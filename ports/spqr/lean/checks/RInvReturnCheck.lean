import Spqr.Proofs.RInvFrame
import Mathlib.Tactic.IntervalCases

namespace Spqr.RInvReturnCheck

open WalkM WalkState

def g : Graph := ⟨3, #[(0, 1), (1, 2), (1, 2), (1, 2), (0, 2)]⟩

def child : DfsTree := .node 2 [
  .back 4 0 (.ret 0 .backEdge),
  .back 2 1 (.ret 1 .backEdge),
  .back 3 1 (.ret 1 .backEdge)]

def childOut : DfsOut := .tree 1 (.ret 0 .type1Child) child

def rootOut : DfsOut := .tree 0 .component (.node 1 [childOut])

theorem dfsForest_eq : g.dfsForest [] [] = [.node 0 [rootOut]] := by cbv

def before : WalkState :=
  let s := { WalkState.init g false with stackVerts := #[0, 0, 0] }
  let s := (walkOutPre 0 0 rootOut false).run s |>.2
  let s := { s with firstOccurrence := s.firstOccurrence.set! 0 g.ne }
  let s := { s with stackVerts := s.stackVerts.set! 1 1 }
  let s := (walkOutPre 1 1 childOut false).run s |>.2
  { s with firstOccurrence := s.firstOccurrence.set! 1 g.ne }

def returned : WalkState := after (walkTree child 2) before

def pEntry : TEntry := ⟨2, 1, 1, ([], [9])⟩

theorem pEntry_mem : pEntry ∈ returned.tstack := by
  cbv
  exact List.mem_cons_self ..

theorem pEntry_piece : 9 ∈ returned.entryPieceItems pEntry := by
  cbv
  exact List.mem_cons_self ..

theorem pEntry_vs : Items.vs returned.items 9 = (some 2, some 1) := by cbv

theorem pending_no_parent {j : Nat} (hj : j = 4 ∨ j = 5) (p : Nat) :
    ¬Items.IsParent returned.items p j := by
  by_cases hp : p < 10
  · interval_cases p <;> rcases hj with rfl | rfl <;> cbv <;>
      change ¬(_ : Nat) ∈ (_ : List Nat) <;> decide
  · have hsize : returned.items.size = 10 := by cbv
    have hle : returned.items.size ≤ p := by omega
    simp [Items.IsParent, Items.ch_of_le _ _ hle]

theorem pending_outside {e : Nat} (he : e = 0 ∨ e = 1) :
    ¬Items.EdgeBelow g returned.items 9 e := by
  intro h
  have hroot : ∀ p, ¬Items.IsParent returned.items p (edgeItem g e) :=
    pending_no_parent (by rcases he with rfl | rfl <;> decide)
  have heq := Items.Below.eq_of_no_parent hroot h
  rcases he with rfl | rfl <;> cases heq

theorem not_sepClass : ¬g.SepClass 2 1 1 0 := by
  rintro (h | ⟨x, y, ⟨z, hx⟩, _, hr⟩)
  · cases h
  · simp [Graph.Joins, g] at hx
    rcases hx with ⟨rfl, _⟩ | ⟨_, rfl⟩
    · exact hr.ok_left.2 rfl
    · exact hr.ok_left.1 rfl

theorem returned_not_rInvAt (dfs : DfsData) : ¬returned.RInvAt dfs 1 := by
  intro h
  have hentry := h.entries pEntry pEntry_mem (by decide)
  have hg : returned.g = g := by cbv
  have hmax := hentry.maximal 9 pEntry_piece 2 1 pEntry_vs
  rw [hg] at hmax
  exact not_sepClass (hmax 1 0 (by decide) (by decide)
    (pending_outside (.inr rfl)) (pending_outside (.inl rfl)))

theorem returned_base :
    returned.tstack.drop (returned.tstack.length - before.tstack.length) = before.tstack := by
  cbv

theorem returned_base_entries (dfs : DfsData) :
    ∀ t ∈ returned.tstack.drop (returned.tstack.length - before.tstack.length),
      t.vStart ≠ 1 → returned.EntryR dfs t := by
  rw [returned_base]
  intro t ht hne
  change t ∈ [⟨1, 1, 0, ([], [2])⟩] at ht
  have ht' : t = ⟨1, 1, 0, ([], [2])⟩ := by simpa using ht
  subst t
  exact False.elim (hne rfl)

#print axioms pEntry_mem
#print axioms dfsForest_eq
#print axioms pEntry_piece
#print axioms pEntry_vs
#print axioms pending_no_parent
#print axioms pending_outside
#print axioms not_sepClass
#print axioms returned_not_rInvAt
#print axioms returned_base
#print axioms returned_base_entries

end Spqr.RInvReturnCheck

open Spqr WalkM

namespace RReturnCheck

def entryKey (t : TEntry) := (t.vStart, t.topDepth, t.firstIdx, t.spans)

def descendants (items : Items) : Nat → Nat → List Nat
  | 0, i => [i]
  | fuel + 1, i => i :: (items.ch i).flatMap (descendants items fuel)

def below (s : WalkState) (i : Nat) : List Nat :=
  ((descendants s.items s.items.size i).filter (1 + s.g.nv ≤ ·)).mergeSort

def entryBelow (s : WalkState) (t : TEntry) : List Nat :=
  (t.spans.1 ++ t.spans.2).flatMap (below s)

def isBlock (g : Graph) : Bool := (List.range g.nv).all fun v => Id.run do
  let mut seen := if g.ne == 0 then [] else [0]
  for _ in List.range g.ne do
    seen := (List.range g.ne).filter fun e => seen.contains e || seen.any fun f =>
      [g.edges[e]!.1, g.edges[e]!.2].any fun u =>
        u != v && (g.edges[f]!.1 == u || g.edges[f]!.2 == u)
  return seen.length == g.ne

/-- Check return shape/disjointness on every graph, and the `EntryR` observation frame on blocks.
The frame check assumes the input's settled-base obligations; empty entries need no terminal. -/
def checkReturn (parent : Nat) (s s' : WalkState) : List String := Id.run do
  let mut bad := []
  let base := s'.tstack.drop (s'.tstack.length - s.tstack.length)
  if s.tstack.length > s'.tstack.length then bad := "size" :: bad
  if base.map entryKey != s.tstack.map entryKey then bad := "base" :: bad
  if s.g.nv != s'.g.nv || s.g.edges != s'.g.edges then bad := "graph" :: bad
  for t in base do
    if isBlock s.g && t.vStart != parent then
      if !(entryBelow s t).isEmpty && s.stackVerts[t.topDepth]! != s'.stackVerts[t.topDepth]! then
        bad := "base terminal" :: bad
      for i in t.spans.1 ++ t.spans.2 do
        if Items.type s.items i != Items.type s'.items i || Items.vs s.items i != Items.vs s'.items i ||
            below s i != below s' i then bad := "base item frame" :: bad
  for i in List.range s'.tstack.length do
    for j in List.range i do
      if (entryBelow s' s'.tstack[i]!).any ((entryBelow s' s'.tstack[j]!).contains) then
        bad := "disjointness" :: bad
  return bad

mutual
partial def tree (t : DfsTree) (d : Nat) (s : WalkState) : WalkState × List String :=
  match t with
  | .node v os =>
    let s := { s with stackVerts := s.stackVerts.set! d v }
    let (hv, s, bad) := outs v d os false s
    let s := if hv then s else (setStackDir d true *> pushVertTstack v d).run s |>.2
    (s, bad)

partial def outs (v d : Nat) (os : List DfsOut) (hv : Bool) (s : WalkState) :
    Bool × WalkState × List String :=
  match os with
  | [] => (hv, s, [])
  | o :: rest =>
    let (hv, s) := (walkOutPre v d o hv).run s
    let orig := s.tstack.length
    let (s, bad) := match o with
      | .back .. => (s, [])
      | .tree _ _ child =>
        let before := { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }
        let (after, bad) := tree child (d + 1) before
        (after, bad ++ (checkReturn v before after).map (s!"v={v} e={o.e} " ++ ·))
    let (hv, s) := (finishEdge v d o orig hv).run s
    let (hv, s, bad') := outs v d rest hv s
    (hv, s, bad ++ bad')
end

def runGraph (g : Graph) (vo eo : List Nat) (tern : Bool) : List String := Id.run do
  let forest := g.dfsForest vo eo
  let mut s := WalkState.init g tern
  let mut bad := []
  for t in forest do
    let (s', bs) := tree t 0 s
    bad := bad ++ bs
    s := (do
      let top ← popTstack
      modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 }).run s' |>.2
  let ref := g.walk tern forest
  if s.items.toList.map (fun i => (i.type, i.vs, i.ch)) !=
      ref.items.toList.map (fun i => (i.type, i.vs, i.ch)) ||
      s.tstack.map entryKey != ref.tstack.map entryKey then
    bad := "instrumented walk differs" :: bad
  return bad

end RReturnCheck

/-- Input is concatenated `gen.py` cases, prefixed by the number of cases. -/
def main : IO UInt32 := do
  let input ← (← IO.getStdin).readToEnd
  let toks := ((input.splitOn " ").flatMap (·.splitOn "\n") |>.filter (· != "") |>.map String.toNat!).toArray
  let mut p := 1
  let mut fails := 0
  let mut blocks := 0
  for seed in [0:toks[0]!] do
    let nv := toks[p]!; let ne := toks[p+1]!
    p := p + 3
    let edges := (List.range ne).map fun i => (toks[p+2*i]!, toks[p+2*i+1]!)
    p := p + 2*ne
    let k := toks[p]!; p := p + 1
    let vo := (List.range k).map fun i => toks[p+i]!
    p := p + k
    let k := toks[p]!; p := p + 1
    let eo := (List.range k).map fun i => toks[p+i]!
    p := p + k
    if RReturnCheck.isBlock ⟨nv, edges.toArray⟩ then blocks := blocks + 1
    for tern in [false, true] do
      let bad := RReturnCheck.runGraph ⟨nv, edges.toArray⟩ vo eo tern
      if !bad.isEmpty then
        fails := fails + 1
        IO.println s!"seed={seed} tern={tern}: {bad}"
  let fixed := RReturnCheck.checkReturn 1 Spqr.RInvReturnCheck.before Spqr.RInvReturnCheck.returned
  if !fixed.isEmpty then
    fails := fails + 1
    IO.println s!"fixed counterexample: {fixed}"
  IO.println s!"cases={toks[0]!} x both ternarize modes; block cases={blocks}; fixed regression; failures={fails}"
  return if fails == 0 then 0 else 1
