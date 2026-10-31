import Spqr.Dfs
import Spqr.Items

/-!
# Phase 2: the ear-decomposition walk

Walking the DFS forest, we maintain a stack of partially built "ears" (`TEntry`), each covering a
tree path from `vStart` up to depth `topDepth`, with the items hanging off the two sides of the
path. Returning from a child merges and closes ears into S / P / R nodes.
-/

namespace Spqr

/-- `sides dir a b` has `a` on side `dir` and `b` on the other side. -/
def setSides (dir : Bool) (a b : α) : α × α := if dir then (b, a) else (a, b)
def getSide (p : α × α) (dir : Bool) : α := if dir then p.2 else p.1

structure TEntry where
  vStart : Nat
  topDepth : Nat
  firstIdx : Nat
  spans : List ItemId × List ItemId
deriving Repr, Inhabited

structure WalkState where
  g : Graph
  ternarize : Bool
  items : Array Item
  stackVerts : Array Nat
  stackDir : Array Bool
  /-- Counts back edges seen so far. -/
  nxtEdgeIdx : Nat := 0
  /-- Index of the first back edge into each depth. -/
  firstOccurrence : Array Nat
  /-- Head is the top (`cur`), second is `nxt`. -/
  tstack : List TEntry := []
  totBlocks : Nat := 0
  totSelfLoops : Nat := 0

namespace WalkState

def init (g : Graph) (ternarize : Bool) : WalkState where
  g := g
  ternarize := ternarize
  items := initialItems g
  stackVerts := Array.replicate g.nv 0
  stackDir := Array.replicate g.nv false
  firstOccurrence := Array.replicate g.nv 0

end WalkState

abbrev WalkM := StateM WalkState

namespace WalkM

def modifyItem (item : ItemId) (f : Item → Item) : WalkM Unit :=
  modify fun s => { s with items := s.items.modify item f }

def getItem (item : ItemId) : WalkM Item := do return (← get).items[item]!

def allocItem (type : NodeType) : WalkM ItemId := do
  let s ← get
  set { s with items := s.items.push ⟨type, (none, none), []⟩ }
  return s.items.size

def stackDir (d : Nat) : WalkM Bool := do return (← get).stackDir[d]!
def setStackDir (d : Nat) (b : Bool) : WalkM Unit :=
  modify fun s => { s with stackDir := s.stackDir.set! d b }

def makeVs (vStart topDepth : Nat) : WalkM (Option Nat × Option Nat) := do
  let s ← get
  return setSides s.stackDir[topDepth]! (some s.stackVerts[topDepth]!) (some vStart)

def cur : WalkM TEntry := do return (← get).tstack.head!
def nxt : WalkM TEntry := do return (← get).tstack.tail.head!
def tstackSize : WalkM Nat := do return (← get).tstack.length
def modifyCur (f : TEntry → TEntry) : WalkM Unit :=
  modify fun s => { s with tstack := match s.tstack with | a :: rest => f a :: rest | [] => [] }
def modifyNxt (f : TEntry → TEntry) : WalkM Unit :=
  modify fun s => { s with tstack := match s.tstack with | a :: b :: rest => a :: f b :: rest | l => l }
def popTstack : WalkM TEntry := do
  let s ← get
  set { s with tstack := s.tstack.tail }
  return s.tstack.head!

def pushTstack (vStart topDepth : Nat) (item : ItemId) : WalkM Unit :=
  modify fun s =>
    { s with tstack := ⟨vStart, topDepth, s.nxtEdgeIdx, setSides s.stackDir[topDepth]! [item] []⟩ :: s.tstack }
def pushVertTstack (v topDepth : Nat) : WalkM Unit := pushTstack v topDepth (vertItem v)
def pushEdgeTstack (vStart topDepth e : Nat) : WalkM Unit := do
  pushTstack vStart topDepth (edgeItem (← get).g e)

/-- Merge the top ear into the one below it. -/
def mergeTstackTops : WalkM Unit := do
  let b ← popTstack
  modifyCur fun a =>
    { a with topDepth := min a.topDepth b.topDepth,
             spans := (b.spans.1 ++ a.spans.1, a.spans.2 ++ b.spans.2) }

/-- Allocate a node of the given type for the ear below the top, unless (when not ternarizing)
that ear is a single node of the same type, which is then reopened and reused. -/
def maybeUnwrapNxt (type : NodeType) : WalkM ItemId := do
  if type == .R || (← get).ternarize then return ← allocItem type
  let t ← nxt
  let topDir ← stackDir t.topDepth
  let item := (getSide t.spans topDir).head!
  let it ← getItem item
  if it.type == type then
    modifyNxt fun t => { t with spans := setSides topDir it.ch [] }
    return item
  else
    allocItem type

/-- Close the top ear into `item`: it becomes the ear's single child. -/
def finishTstackTop (item : ItemId) : WalkM Unit := do
  let t ← cur
  let topDir ← stackDir t.topDepth
  let vs ← makeVs t.vStart t.topDepth
  modifyItem item fun it => { it with vs := vs, ch := getSide t.spans topDir }
  modifyCur fun t => { t with spans := setSides topDir [item] [] }

/-- `while cond do body`, with `fuel` bounding the number of iterations. -/
def loop : Nat → WalkM Bool → WalkM Unit → WalkM Unit
  | 0, _, _ => pure ()
  | fuel + 1, cond, body => do
    if ← cond then
      body
      loop fuel cond body

end WalkM

open WalkM in
/-- Finish out-edge `o` of vertex `curV` at depth `d`, after its subtree (if any) has been walked.
`hasVert` says whether `cur`'s own vertex ear has been pushed; returns the updated flag. -/
def finishEdge (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool) : WalkM Bool := do
  let g := (← get).g
  let nxtV := o.dest
  let e := o.e
  let lowval := o.cls.lowval d
  let isTree := o.cls.isTree
  let isType1 := o.cls.isType1
  let edgeDir ← stackDir d
  let qItem := edgeItem g e

  if lowval ≥ d then
    -- A block boundary: bridges, components hanging from a cut vertex, and self-loops.
    modifyItem qItem fun it => { it with vs := (some curV, none) }
    modify fun s => { s with totBlocks := s.totBlocks + 1 }
    if isTree then
      if lowval == d + 1 then
        -- The top ear is just the child vertex; prepend a bridge node.
        let item ← allocItem .I
        let vs ← makeVs nxtV d
        modifyItem item fun it => { it with vs := vs }
        let t ← popTstack
        modifyItem qItem fun it => { it with ch := item :: t.spans.2 }
      else
        -- Below the top: the vertex ear, then the back-edge ear.
        let backedge ← popTstack
        let t ← popTstack
        modifyItem qItem fun it => { it with ch := backedge.spans.1 ++ t.spans.2 }
    else
      modify fun s => { s with totSelfLoops := s.totSelfLoops + 1 }
      let item ← allocItem .O
      modifyItem item fun it => { it with vs := (some curV, none) }
      modifyItem qItem fun it => { it with ch := [item] }
    modifyItem (vertItem curV) fun it => { it with ch := it.ch ++ [qItem] }
    return hasVert

  let vs ← makeVs nxtV d
  modifyItem qItem fun it => { it with vs := vs }

  let mut isSingle := true
  if isTree then
    pushEdgeTstack nxtV d e
    -- Close every ear that stays at or below `cur`.
    loop (← tstackSize) (do return (← tstackSize) ≥ 2 && (← nxt).topDepth ≥ d) do
      let type ← do
        if (← nxt).topDepth > d then
          setStackDir (← nxt).topDepth edgeDir
          mergeTstackTops
          pure NodeType.S
        else if (← nxt).vStart == (← cur).vStart then pure .P
        else pure .R
      let item ← maybeUnwrapNxt type
      mergeTstackTops
      finishTstackTop item
    -- Merge ears whose first back edge is after a back edge of `cur` into the top ear.
    let fo := (← get).firstOccurrence[d]!
    if (← cur).firstIdx > fo then
      loop (← tstackSize) (do return (← cur).firstIdx > fo) mergeTstackTops
      isSingle := false
    if hasVert then
      if !isType1 then
        loop (← tstackSize) (do return (← tstackSize) > origTstack + 3) mergeTstackTops
        isSingle := false
      let item ← if isType1 then some <$> maybeUnwrapNxt (if isSingle then .S else .R) else pure none
      mergeTstackTops  -- with the back edge
      mergeTstackTops  -- with the vertex
      modifyCur fun t => { t with vStart := curV, spans := setSides (!edgeDir) (t.spans.1 ++ t.spans.2) [] }
      if let some item := item then
        finishTstackTop item
        isSingle := true
  else
    pushEdgeTstack curV lowval e
    modify fun s =>
      { s with firstOccurrence := s.firstOccurrence.modify lowval (min · s.nxtEdgeIdx),
               nxtEdgeIdx := s.nxtEdgeIdx + 1 }

  if isType1 && (← tstackSize) ≥ 2 && (← nxt).vStart == curV && (← nxt).topDepth == lowval then
    let item ← maybeUnwrapNxt .P
    mergeTstackTops
    finishTstackTop item

  if !hasVert then
    pushVertTstack curV d
    if !isSingle then mergeTstackTops
    return true
  return hasVert

open WalkM in
mutual
/-- Walk the subtree `t`, which sits at depth `d`. -/
def walkTree (t : DfsTree) (d : Nat) : WalkM Unit := do
  match t with
  | .node v outs =>
    modify fun s => { s with stackVerts := s.stackVerts.set! d v }
    let hasVert ← walkOuts v d outs false
    unless hasVert do
      setStackDir d true
      pushVertTstack v d

def walkOuts (v d : Nat) (outs : List DfsOut) (hasVert : Bool) : WalkM Bool := do
  match outs with
  | [] => return hasVert
  | o :: rest =>
    let hasVert ← walkOut v d o hasVert
    walkOuts v d rest hasVert

def walkOut (v d : Nat) (o : DfsOut) (hasVert : Bool) : WalkM Bool := do
  let lowval := o.cls.lowval d
  let lowDir ← stackDir lowval
  setStackDir d (if lowval ≥ d then false else !lowDir)
  let hasVert ← do
    if !hasVert && lowval < d && o.cls.isType1 then
      pushVertTstack v d
      pure true
    else pure hasVert
  let origTstack ← tstackSize
  match o with
  | .tree _ _ child =>
    modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }
    walkTree child (d + 1)
  | .back .. => pure ()
  finishEdge v d o origTstack hasVert
end

open WalkM in
def walkForest (forest : List DfsTree) : WalkM Unit :=
  forest.forM fun t => do
    walkTree t 0
    let top ← popTstack
    modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 }

/-- Phase 2 entry point. -/
def Graph.walk (g : Graph) (ternarize : Bool) (forest : List DfsTree) : WalkState :=
  (walkForest forest).run (WalkState.init g ternarize) |>.2

end Spqr
