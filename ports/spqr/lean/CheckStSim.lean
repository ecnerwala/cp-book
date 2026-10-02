import Spqr
import Spqr.StSim

open Spqr WalkM

/-! Differential test of the simulation relation `StSim` / `StSimOuts` (`Spqr/StSim.lean`): the walk
is re-run step by step (mirroring `walkTree` / `walkOuts` / `walkOut`, calling the real
`finishEdge`) and the relation is checked at the start of every out-edge, before every `finishEdge`
and at the end of every `walkTree`. -/

instance (p q : Nat × Nat) : Decidable (Items.PairEq p q) := by unfold Items.PairEq; infer_instance

structure Stats where
  checks : Nat := 0
  bad : Nat := 0

instance : BEq TEntry :=
  ⟨fun a b => a.vStart == b.vStart && a.topDepth == b.topDepth && a.firstIdx == b.firstIdx && a.spans == b.spans⟩

def isSPR (items : Items) (i : ItemId) : Bool :=
  match Items.type items i with
  | .S | .P | .R => true
  | _ => false

def isInfixB (l m : List Nat) : Bool :=
  (List.range (m.length + 1)).any fun k => (m.drop k).take l.length == l

def descOf (items : Items) : Nat → ItemId → List ItemId
  | 0, i => [i]
  | fuel + 1, i => i :: (Items.ch items i).flatMap (descOf items fuel)

def leavesB (items : Items) (i : ItemId) : List ItemId := Items.leaves items items.size i

def stReadB (items : Items) (new : List TEntry) (ps : List StPiece) : Bool :=
  (readStack new).flatMap (leavesB items) == stNest ps

def vsOrientedAtB (g : Graph) (items : Items) (b : StBlock) (i : ItemId) : Bool :=
  let lv := leavesB items i
  let ds := descOf items items.size i
  lv.all (· ∈ b.items) &&
  ((List.range g.ne).all fun e =>
    let x := edgeItem g e
    !(x ∈ b.items || b.root.any fun r => decide (Items.PairEq g.edges[e]! r)) || !(x ∈ ds) || x ∈ lv) &&
  decide (Oriented (b.seq g) (Items.vs items i)) &&
  (match Items.vs items i with
    | (some s, some t) => (Items.ch items i).all fun c => Items.type items c ≠ .V ||
        (decide (Precedes (b.seq g) s (c - 1)) && decide (Precedes (b.seq g) (c - 1) t))
    | _ => false) &&
  ((Items.ch items i).all fun c => Items.type items c = .V || decide (Oriented (b.seq g) (Items.vs items c))) &&
  ((Items.ch items i).all fun c => Items.type items c = .V ||
    match Items.vs items c with
    | (some u, some v) =>
      let dc := descOf items items.size c
      (List.range g.ne).all fun e =>
        !(edgeItem g e ∈ b.items || b.root.any fun r => decide (Items.PairEq g.edges[e]! r)) ||
        !(edgeItem g e ∈ dc) ||
        [(g.edges[e]!).1, (g.edges[e]!).2].all fun y =>
          (u == y || decide (Precedes (b.seq g) u y)) && (y == v || decide (Precedes (b.seq g) y v))
    | _ => false)

def inBlockB (g : Graph) (items : Items) (b : StBlock) (i : ItemId) : Bool :=
  isInfixB (leavesB items i) b.items && vsOrientedAtB g items b i

/-- Violations of `StItems`. -/
def stItemsB (g : Graph) (s : WalkState) (blocks : List StBlock) : List String :=
  let rs := readStack s.tstack
  let items := s.items
  let m1 := if rs.all (fun x => (List.range items.size).all fun p => !(x ∈ Items.ch items p)) then []
    else ["roots: a stack item has a parent"]
  let m2 := if decide rs.Nodup then [] else [s!"nodup: {rs}"]
  let m4 := if rs.all (fun x => (descOf items items.size x).all (· < items.size)) then []
    else ["bounded: a stack item's subtree is out of range"]
  let m3 := (List.range items.size).filterMap fun i =>
    if !isSPR items i then none
    else if rs.any (fun x => i ∈ descOf items items.size x) then none
    else if blocks.any (inBlockB g items · i) then none
    else some s!"item {i} ({repr (Items.type items i)}) vs {repr (Items.vs items i)} ch {Items.ch items i} leaves {leavesB items i} not live and in no block {blocks.map (·.items)}"
  let m5 := if (List.range items.size).all (fun p => (Items.ch items p).all (· < items.size)) then []
    else ["chLt: a child is out of range"]
  let m6 := if (List.range items.size).all (fun p => decide (Items.ch items p).Nodup) then []
    else ["chNodup: a child list repeats"]
  m1 ++ m2 ++ m4 ++ m5 ++ m6 ++ m3

def report (st : IO.Ref Stats) (ok : Bool) (msg : String) : IO Unit := do
  st.modify fun x => { x with checks := x.checks + 1, bad := if ok then x.bad else x.bad + 1 }
  unless ok do IO.println s!"BAD {msg}"

def above (s : WalkState) (base : List TEntry) : List TEntry := s.tstack.take (s.tstack.length - base.length)
def baseOk (s : WalkState) (base : List TEntry) : Bool :=
  base.length ≤ s.tstack.length && s.tstack.drop (s.tstack.length - base.length) == base

mutual
partial def chkTree (st : IO.Ref Stats) (g : Graph) (prev : List DfsTree) (fs : List PathFrame)
    (t : DfsTree) (d : Nat) (s : WalkState) : IO WalkState := do
  match t with
  | .node v outs =>
    let s : WalkState := { s with stackVerts := s.stackVerts.set! d v }
    let base := s.tstack
    let (hasVert, s) ← chkOuts st g prev fs v d outs [] false base s
    let s := if hasVert then s else ((do setStackDir d true; pushVertTstack v d : WalkM Unit).run s).2
    let (ps, _) := refTree g t d (DirsOf s d)
    report st (baseOk s base) s!"base changed at end of walkTree {v} d={d}"
    report st (stReadB s.items (above s base) ps)
      s!"read at end of walkTree {v} d={d}: {(readStack (above s base)).flatMap (leavesB s.items)} vs {stNest ps}"
    for m in stItemsB g s (refBlocks g (prev ++ [truncTree fs t])) do
      report st false s!"items at end of walkTree {v} d={d}: {m}"
    return s

partial def chkOuts (st : IO.Ref Stats) (g : Graph) (prev : List DfsTree) (fs : List PathFrame)
    (v d : Nat) (outs done : List DfsOut) (hasVert : Bool) (base : List TEntry) (s : WalkState) :
    IO (Bool × WalkState) := do
  let (ps, _, hv) := refOuts g v d (DirsOf s d) done false
  report st (baseOk s base) s!"base changed at out-edge {done.length} of {v} d={d}"
  report st (stReadB s.items (above s base) ps)
    s!"read at out-edge {done.length} of {v} d={d}: {(readStack (above s base)).flatMap (leavesB s.items)} vs {stNest ps}"
  report st (hv == hasVert) s!"hasVert at out-edge {done.length} of {v} d={d}: ref {hv} walk {hasVert}"
  for m in stItemsB g s (refBlocks g (prev ++ [truncTree fs (.node v done)])) do
    report st false s!"items at out-edge {done.length} of {v} d={d}: {m}"
  match outs with
  | [] => return (hasVert, s)
  | o :: rest =>
    let (hasVert, s) ← chkOut st g prev fs v d o done hasVert s
    chkOuts st g prev fs v d rest (done ++ [o]) hasVert base s

partial def chkOut (st : IO.Ref Stats) (g : Graph) (prev : List DfsTree) (fs : List PathFrame)
    (v d : Nat) (o : DfsOut) (done : List DfsOut) (hasVert : Bool) (s : WalkState) :
    IO (Bool × WalkState) := do
  let lowval := o.cls.lowval d
  let lowDir := s.stackDir[lowval]!
  let s := ((setStackDir d (if lowval ≥ d then false else !lowDir)).run s).2
  let (hasVert, s) := if !hasVert && lowval < d && o.cls.isType1 then (true, ((pushVertTstack v d).run s).2)
    else (hasVert, s)
  let orig := s.tstack
  let s ← match o with
    | .tree _ _ child =>
      let s : WalkState := { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }
      chkTree st g prev (fs ++ [⟨v, done, o⟩]) child (d + 1) s
    | .back .. => pure s
  let psChild := match o with
    | .tree _ _ child => (refTree g child (d + 1) (DirsOf s (d + 1))).1
    | .back .. => []
  report st (baseOk s orig) s!"orig changed before finishEdge {o.e} of {v} d={d}"
  report st (stReadB s.items (above s orig) psChild)
    s!"read before finishEdge {o.e} of {v} d={d}: {(readStack (above s orig)).flatMap (leavesB s.items)} vs {stNest psChild}"
  for m in stItemsB g s (refBlocks g (prev ++ [truncTree fs (.node v (done ++ [o]))])) do
    report st false s!"items before finishEdge {o.e} of {v} d={d}: {m}"
  return (finishEdge v d o orig.length hasVert).run s
end

def main : IO Unit := do
  let input ← (← IO.getStdin).readToEnd
  let toks := (input.splitOn " ").flatMap (·.splitOn "\n") |>.filter (· ≠ "") |>.map String.toNat!
  let toks := toks.toArray
  let mut p := 0
  let next := fun (p : Nat) => (toks[p]!, p + 1)
  let (nv, p1) := next p; let (ne, p2) := next p1; let (tern, p3) := next p2
  p := p3
  let mut edges : Array (Nat × Nat) := #[]
  for _ in [0:ne] do
    edges := edges.push (toks[p]!, toks[p+1]!)
    p := p + 2
  let k := toks[p]!; p := p + 1
  let vertOrder := (List.range k).map fun i => toks[p + i]!
  p := p + k
  let l := toks[p]!; p := p + 1
  let edgeOrder := (List.range l).map fun i => toks[p + i]!
  let g : Graph := ⟨nv, edges⟩
  let forest := g.dfsForest vertOrder edgeOrder
  let st ← IO.mkRef ({} : Stats)
  let mut s := WalkState.init g (tern != 0)
  let mut prev : List DfsTree := []
  for t in forest do
    s ← chkTree st g prev [] t 0 s
    s := ((do let top ← popTstack; modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 } : WalkM Unit).run s).2
    prev := prev ++ [t]
  let real := g.walk (tern != 0) forest
  report st (toString (repr s.items) == toString (repr real.items)) "mirrored walk differs from g.walk"
  let stats ← st.get
  IO.println s!"checks {stats.checks} bad {stats.bad}"
