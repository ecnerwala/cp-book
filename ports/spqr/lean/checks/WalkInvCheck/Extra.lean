import WalkInvCheck.St
/-!
# The `WalkInv` fields no existing checker evaluates

* `invCheck`: `Inv' D s` (`Spqr/WalkSpec.lean`: `EntryInv'` per open entry with the entries above
  it, `ItemInv` per allocated node);
* `rskelCheck`: `Items.RSkelInv` in its `Items.RThreeConnected` form (`items.rSkeleton g i` of every
  R item is 3-connected, `SpqrTree.ThreeConnected`; `Items.rSkeleton_perm_contract` relates the two);
* the R handoff candidates of PROOF.md §4.5: `e1Check`/`e3Check` at the child-entry site (`entry`,
  `buried_vacuous`), `e2Check` (`bd_free`: no open entry owns an edge below `vertItem v`),
  `e4Check` (span items are roots, every item has at most one parent), `e5Check` (`mid_cur` at the
  reached splits of loop 1);
* `stLiveB`: `StLive g items new b` (`Spqr/StBdPop.lean`) of one stack segment against the block
  `openBlock …` of the truncated reference — the per-segment pairing PROOF.md §7 lists as not yet
  checked.
-/
namespace WalkInvCheck.Extra
open Spqr WalkM

partial def edgesBelow (s : WalkState) (i : ItemId) : List Nat :=
  (if 1 + s.g.nv ≤ i ∧ i < 1 + s.g.nv + s.g.ne then [i - 1 - s.g.nv] else []) ++
    (s.items[i]!.ch).flatMap (edgesBelow s)
def spanItems (t : TEntry) : List ItemId := t.spans.1 ++ t.spans.2
def entryEdges (s : WalkState) (t : TEntry) : List Nat := (spanItems t).flatMap (edgesBelow s)
def inc (s : WalkState) (e v : Nat) : Bool := s.g.edges[e]!.1 == v || s.g.edges[e]!.2 == v
def touches (s : WalkState) (E : List Nat) (v : Nat) : Bool := E.any (inc s · v)
def isParent (s : WalkState) (p c : ItemId) : Bool := (s.items[p]!.ch).contains c
def hasParent (s : WalkState) (c : ItemId) : Bool :=
  (List.range s.items.size).any fun p => isParent s p c
def showT (t : TEntry) : String := s!"({t.vStart},{t.topDepth},{t.firstIdx},{t.spans})"
def showStack (s : WalkState) : String := s!"{s.tstack.map showT}"

/-- `g.ConnEdges E`: the endpoints of the edges of `E` are connected through `E`. -/
def connB (s : WalkState) (E : List Nat) : Bool :=
  match E with
  | [] => true
  | e₀ :: _ => Id.run do
    let mut seen : List Nat := [s.g.edges[e₀]!.1]
    for _ in List.range (s.g.nv + 1) do
      for e in E do
        let (a, b) := s.g.edges[e]!
        if seen.contains a && !seen.contains b then seen := b :: seen
        if seen.contains b && !seen.contains a then seen := a :: seen
    return E.all fun e => seen.contains s.g.edges[e]!.1 && seen.contains s.g.edges[e]!.2

/-- The attachment vertices of `E`: incident to an edge of `E` and to an edge outside `E`. -/
def boundary (s : WalkState) (E : List Nat) : List Nat :=
  (List.range s.g.nv).filter fun v =>
    touches s E v && (List.range s.g.ne).any fun e => !E.contains e && inc s e v
def termB (D : Nat) (s : WalkState) (t : TEntry) (v : Nat) : Bool :=
  v == t.vStart || (List.range (D + 1)).any fun k => t.topDepth ≤ k && s.stackVerts[k]! == v
def term'B (D : Nat) (s : WalkState) (above : List TEntry) (t : TEntry) (v : Nat) : Bool :=
  termB D s t v || above.any fun t' => termB D s t' v
def twoAttB (s : WalkState) (E : List Nat) (a b : Nat) : Bool :=
  (boundary s E).all fun v => v == a || v == b

/-- The failing fields of `s.Inv' D`. -/
def invCheck (seed : Nat) (site : String) (D : Nat) (s : WalkState) : List Viol := Id.run do
  let mut out : List Viol := []
  let bad (k info : String) : Viol :=
    ⟨seed, s.ternarize, s!"inv.{k}", s!"{site} D={D} {info} | sv={s.stackVerts.toList.take (D+2)} stack={showStack s}"⟩
  for i in List.range s.tstack.length do
    let t := s.tstack[i]!
    let above := s.tstack.take i
    let E := entryEdges s t
    if !connB s E then out := bad "conn" (showT t) :: out
    let badB := (boundary s E).filter fun v => !term'B D s above t v
    if badB ≠ [] then out := bad "attached" s!"{showT t} bad={badB}" :: out
  for i in List.range s.items.size do
    if 1 + s.g.nv + s.g.ne ≤ i then
      let E := edgesBelow s i
      if !connB s E then out := bad "node_conn" s!"i={i}" :: out
      match Items.vs s.items i with
      | (some u, some v) => if !twoAttB s E u v then out := bad "node_att" s!"i={i} vs=({u},{v})" :: out
      | _ => pure ()
  return out

/-- `SpqrTree.ThreeConnected n es`, decidably. -/
def threeConnB (n : Nat) (es : List (Nat × Nat)) : Bool :=
  4 ≤ n && (List.range n).all fun a => (List.range n).all fun b =>
    let alive := fun v => v < n && v != a && v != b
    let vs := (List.range n).filter alive
    match vs with
    | [] => true
    | u :: _ => Id.run do
      let mut seen := [u]
      for _ in List.range (n + 1) do
        for (x, y) in es do
          if alive x && alive y then
            if seen.contains x && !seen.contains y then seen := y :: seen
            if seen.contains y && !seen.contains x then seen := x :: seen
      return vs.all seen.contains

/-- `Items.RSkelInv` (as `Items.RThreeConnected`): every R item's skeleton is 3-connected. -/
def rskelCheck (seed : Nat) (site : String) (s : WalkState) : List Viol :=
  (List.range s.items.size).flatMap fun i =>
    if Items.type s.items i != .R then [] else
    let bad (k info : String) : Viol := ⟨seed, s.ternarize, s!"rskel.{k}", s!"{site} i={i} {info}"⟩
    (match Items.vs s.items i with
      | (some _, some _) => []
      | vs => [bad "vs" s!"vs={repr vs}"]) ++
    (if threeConnB (Items.nvList s.g s.items i).length (Items.rSkeleton s.g s.items i) then []
     else [bad "three_conn" s!"nv={Items.nvList s.g s.items i} skel={Items.rSkeleton s.g s.items i} ch={Items.ch s.items i}"])

/-- `g.Reach ok x y` (`Spqr/SepPair.lean`): a walk from `x` to `y` all of whose vertices satisfy `ok`. -/
def reachB (s : WalkState) (ok : Nat → Bool) (x y : Nat) : Bool :=
  if !ok x then false else Id.run do
    let mut seen := [x]
    for _ in List.range (s.g.nv + 1) do
      for e in List.range s.g.ne do
        let (a, b) := s.g.edges[e]!
        if ok a && ok b then
          if seen.contains a && !seen.contains b then seen := b :: seen
          if seen.contains b && !seen.contains a then seen := a :: seen
    return seen.contains y
/-- `g.SepClass a b e e'` = `g.EdgeConn (· ∉ {a, b}) e e'`. -/
def sepClassB (s : WalkState) (a b e e' : Nat) : Bool :=
  e == e' || Id.run do
    let ok := fun v => v != a && v != b
    let (x₁, x₂) := s.g.edges[e]!
    let (y₁, y₂) := s.g.edges[e']!
    return [x₁, x₂].any fun x => [y₁, y₂].any fun y => reachB s ok x y

/-- E1 (`entry`) and the restated E3 (`stab_single`, PROOF.md §4.5) at the child-entry site
`walkTree c D` (`D = d + 1 ≥ 1`, parent `p = stackVerts[d]`): no open entry starts at `p` with
`topDepth > d`; every open entry with `topDepth = D` (the only ones whose `EntryR.single` reads the
updated `stackVerts[D]`) has `EntryR.single` in the updated state:
`TwoAttached (t.edges) t.vStart c ∨ ∀ e e' ∈ t.edges, SepClass t.vStart c e e'`. The original E3
(`buried_vacuous`: single edge or `TwoAttached … vStart vStart`) is refuted in `E3False.lean`. -/
def e3Check (seed D c : Nat) (s : WalkState) : List Viol :=
  if D = 0 then [] else
  let p := s.stackVerts[D - 1]!
  s.tstack.flatMap fun t =>
    let es := entryEdges s t
    (if t.vStart == p && t.topDepth > D - 1 then
      [⟨seed, s.ternarize, "e1.entry", s!"D={D} p={p} {showT t} | {showStack s}"⟩] else []) ++
    (if t.topDepth != D || twoAttB s es t.vStart c ||
        es.all (fun e => es.all fun e' => sepClassB s t.vStart c e e') then [] else
      [⟨seed, s.ternarize, "e3.stab_single", s!"D={D} c={c} {showT t} edges={es} bd={boundary s es} sv={s.stackVerts.toList.take (D+3)} | {showStack s}"⟩])

/-- E2 (`bd_free`): no open entry owns an edge below `vertItem v`. -/
def e2Check (seed : Nat) (site : String) (v : Nat) (s : WalkState) : List Viol :=
  let EV := edgesBelow s (vertItem v)
  s.tstack.flatMap fun t =>
    if (entryEdges s t).any EV.contains then
      [⟨seed, s.ternarize, "e2.bd_free", s!"{site} v={v} {showT t} EV={EV} | {showStack s}"⟩] else []

/-- E4: span items are roots; every item has at most one parent (counted with multiplicity). -/
def e4Check (seed : Nat) (site : String) (s : WalkState) : List Viol := Id.run do
  let mut out : List Viol := []
  for t in s.tstack do
    for i in spanItems t do
      if hasParent s i then out := ⟨seed, s.ternarize, "e4.span_root", s!"{site} {showT t} i={i}"⟩ :: out
  for i in List.range s.items.size do
    let k : Nat := (List.range s.items.size).foldl (fun a p => a + (s.items[p]!.ch).count i) 0
    if k > 1 then out := ⟨seed, s.ternarize, "e4.unique_parent", s!"{site} i={i} parents={k}"⟩ :: out
  return out

def interiorB (s : WalkState) (E : List Nat) (v : Nat) : Bool :=
  (List.range s.g.ne).all fun e => !inc s e v || E.contains e

/-- E5 (`mid_cur`) at the reached splits of loop 1 (the `Loop1Spec` walk of `hi`, the `topDepth ≥ d`
prefix of the child's entries, as in `checks/EarCheck.lean`'s `specWalk`): at a split with the
consumed entries `done`, the child `c = stackVerts[d+1]` is the current bottom
(`l1Bot o done = (done.getLast?.map vStart).getD o.dest`) or interior to the edges consumed so far
(`l1Edges o s done = o.e :: edges of done`). `mid_cur` is the split closing at depth `d`
(`t.topDepth = d`, the P/R branch `RBranch.mid` consumes); `mid_cur_s` the series split. -/
def e5Walk (seed : Nat) (s : WalkState) (d : Nat) (o : DfsOut) : Nat → List TEntry → List TEntry → List Viol
  | 0, _, _ => []
  | _ + 1, _, [] => []
  | fuel + 1, done, t :: rest =>
    let v := (done.getLast?.map TEntry.vStart).getD o.dest
    let E := o.e :: done.flatMap (entryEdges s)
    let c := s.stackVerts[d + 1]!
    let bad (k : String) : List Viol :=
      if v == c || interiorB s E c then [] else
      [⟨seed, s.ternarize, s!"e5.{k}", s!"d={d} e={o.e} c={c} bot={v} E={E} t={showT t} done={done.map showT} | {showStack s}"⟩]
    if t.topDepth == d then bad "mid_cur" ++ e5Walk seed s d o fuel (done ++ [t]) rest
    else bad "mid_cur_s" ++ match rest with
      | t' :: rest' => e5Walk seed s d o fuel (done ++ [t, t']) rest'
      | [] => []
def e5Check (seed v d : Nat) (o : DfsOut) (orig : Nat) (s : WalkState) : List Viol :=
  if !(o.cls.isTree && o.cls.lowval d < d) then [] else
  let sub := s.tstack.take (s.tstack.length - orig)
  let hi := sub.takeWhile fun t => d ≤ t.topDepth
  (e5Walk seed s d o (hi.length + 1) [] hi).map fun x => { x with info := s!"v={v} {x.info}" }

/-- `HangingUnder items x i`: some V/Q item `y ≠ i` below `x` has `i` below it. -/
def hangingUnderB (items : Items) (x i : ItemId) : Bool :=
  (St.descOf items items.size x).any fun y =>
    y != i && (Items.type items y == .V || Items.type items y == .Q) && i ∈ St.descOf items items.size y

/-- The failing clauses of `StLive g items new b`. -/
def stLiveB (g : Graph) (items : Items) (new : List TEntry) (b : StBlock) : List String :=
  (readStack new).flatMap fun x => (St.descOf items items.size x).filterMap fun i =>
    if !St.isSPR items i || hangingUnderB items x i || St.inBlockB g items b i then none
    else some s!"item {i} ({repr (Items.type items i)}) vs {repr (Items.vs items i)} ch {Items.ch items i} leaves {St.leavesB items i} under span item {x} not InBlock root={repr b.root} items={b.items}"

end WalkInvCheck.Extra
