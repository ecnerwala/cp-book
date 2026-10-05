import Spqr.Walk
import Spqr.Planar

/-!
# Phase 2, planar variant: the walk with planarity bookkeeping

The planar walk is the ordinary walk (`Spqr.Walk`) run on the same `WalkState`, with extra state
`PlanarAux` carried alongside: every tstack entry gets a `PlEntry` (the flip bits of its span items
and its `Planarity` data), closed S/P/R items record four quarter-edge matches, and the global
`qem` (quarter-edge matches) accumulates the partial rotation system. Every step below is a base
step lifted by `liftW` interleaved with steps that only touch the auxiliary state, so the base
component of the planar walk is the ordinary walk (`planarWalk_proj`, in `Spqr.PlanarWalkProj`).

Virtual edge `ve = item - (1 + nv)` caps item `item`; its quarter-edges are `4 ve + 2 side + dir`.

## Flip bits

The C++ stores span items as `2 * item + flip` and encodes the flip of each next item relative
to the previous one (`ch_nxt` low bit), so that flipping a whole list only touches its two ends.
Here every span item carries its *absolute* flip as a `Bool`, and flipping a list is `map not`.
The two encodings agree: the C++ reader (`relabel`) recovers the absolute flip of the `k`-th item
as the XOR of the head's bit and the first `k` deltas, concatenation stores the delta
`head(b) ^ tail(a)` of absolute bits, and flipping both end bits of a list flips every recovered
absolute bit while leaving the deltas unchanged — which is exactly `map not`.
-/

namespace Spqr

/-- The exposed back-edge ends of one side: `ends.1`/`depths.1` is the outermost (lowest
return) and `ends.2`/`depths.2` the innermost, with depths increasing inwards. -/
structure PlTop where
  ends : Nat × Nat
  depths : Nat × Nat
deriving Repr, Inhabited

/-- One side of an ear's planarity data: `bot` are the outer / inner exposed quarter-edges along
the tree path (attached to the bottom- / top-most vertex), `top` the exposed back edges. -/
structure PlSide where
  bot : Option (Nat × Nat) := none
  top : Option PlTop := none
deriving Repr, Inhabited

/-- Planarity data of a tstack entry: side 0 holds a minimal return (`sides.1.top.depths.1 =
topDepth` whenever there is one). -/
structure Planarity where
  sides : PlSide × PlSide := (default, default)
deriving Repr, Inhabited

/-- Planarity payload of a tstack entry: per-side flips of the span items (parallel to
`TEntry.spans`) and the planarity data, `none` once nonplanarity has been detected. -/
structure PlEntry where
  flips : List Bool × List Bool
  pl : Option Planarity
deriving Repr, Inhabited

/-- Recorded planarity of a closed node: the four quarter-edge matches of its cap (indexed by
`2 * side + dir`), or nonplanar. `unset` for I/O nodes (never read). -/
inductive NodePlanarity where
  | unset
  | planar (m : Array Nat)
  | nonplanar
deriving Repr, Inhabited

/-- Quarter-edge matches: `qem[q] = some r` when `q` and `r` face each other. -/
abbrev Qem := Array (Option Nat)
abbrev QemM := StateM Qem

namespace QemM
def link (a b : Nat) : QemM Unit := modify fun m => (m.set! a (some b)).set! b (some a)
def take (a : Nat) : QemM (Option Nat) := do
  let m ← get
  set (m.set! a none)
  return m[a]!
def clear (a : Nat) : QemM Unit := modify (·.set! a none)
end QemM

structure PlanarAux where
  qem : Qem
  /-- `topDepth` of the entry each virtual edge was created for. -/
  edgeTopDepths : Array Nat
  /-- Indexed by `item - (1 + nv + ne)`. -/
  nodePlanarity : Array NodePlanarity
  /-- Flips of each item's children, parallel to `Item.ch`. -/
  itemFlips : Array (List Bool)
  /-- Parallel to `WalkState.tstack`. -/
  plStack : List PlEntry := []

structure PlanarWalkState where
  base : WalkState
  aux : PlanarAux

def PlanarWalkState.init (g : Graph) (ternarize : Bool) : PlanarWalkState where
  base := WalkState.init g ternarize
  aux := {
    qem := Array.replicate (8 * g.ne + 4) none
    edgeTopDepths := Array.replicate (2 * g.ne) 0
    nodePlanarity := #[]
    itemFlips := Array.replicate (1 + g.nv + g.ne) [] }

abbrev PlanarWalkM := StateM PlanarWalkState

/-- Index `2 * side + dir` into a node's four cap matches. -/
def qidx (side dir : Bool) : Nat := 2 * side.toNat + dir.toNat

def flipEntry (e : PlEntry) : PlEntry :=
  { flips := (e.flips.1.map not, e.flips.2.map not), pl := e.pl.map fun p => ⟨(p.sides.2, p.sides.1)⟩ }

/-- Planarity data of a fresh single-edge ear on virtual edge `ve`. -/
def makeEdgePlanarity (ve topDepth : Nat) (topDir isTree : Bool) : Planarity :=
  if isTree then
    ⟨({ bot := some (QE.mk ve (!topDir) false, QE.mk ve topDir true) },
      { bot := some (QE.mk ve (!topDir) true, QE.mk ve topDir false) })⟩
  else
    ⟨({ bot := some (QE.mk ve (!topDir) false, QE.mk ve (!topDir) true),
        top := some ⟨(QE.mk ve topDir true, QE.mk ve topDir false), (topDepth, topDepth)⟩ },
      default)⟩

/-- The cap matches recorded when closing an ear. -/
def finishMatches (isTree topDir : Bool) (p : Planarity) : Array Nat :=
  let b0 := p.sides.1.bot.getD (0, 0)
  let b1 := p.sides.2.bot.getD (0, 0)
  let t0 := (p.sides.1.top.map (·.ends)).getD (0, 0)
  let m := Array.replicate 4 0
  if isTree then
    (((m.set! (qidx (!topDir) true) b0.1).set! (qidx topDir false) b0.2).set! (qidx (!topDir) false) b1.1).set! (qidx topDir true) b1.2
  else
    (((m.set! (qidx (!topDir) true) b0.1).set! (qidx (!topDir) false) b0.2).set! (qidx topDir false) t0.1).set! (qidx topDir true) t0.2

/-- Reopen a closed node: its recorded matches replace the ends of its cap edge. -/
def unwrapPlanarity (isTree topDir : Bool) (m : Array Nat) (p : Planarity) : Planarity :=
  let s0 := p.sides.1
  let s1 := p.sides.2
  if isTree then
    ⟨({ s0 with bot := some (m[qidx (!topDir) true]!, m[qidx topDir false]!) },
      { s1 with bot := some (m[qidx (!topDir) false]!, m[qidx topDir true]!) })⟩
  else
    let bot := some (m[qidx (!topDir) true]!, m[qidx (!topDir) false]!)
    let top := s0.top.map fun t => { t with ends := (m[qidx topDir false]!, m[qidx topDir true]!) }
    ⟨({ s0 with bot := bot, top := top }, s1)⟩

/-- Merge the side `b` of an ear on top of side `a` of the ear below: the bottom chains are
joined, and the back-edge chains too unless they nest the wrong way (`none`). -/
def mergeSide (a b : PlSide) : QemM (Option PlSide) := do
  match b.bot, a.bot with
  | none, _ => return some a
  | some _, none => return some b
  | some (bb0, bb1), some (ab0, ab1) =>
    QemM.link ab1 bb0
    let a := { a with bot := some (ab0, bb1) }
    match b.top, a.top with
    | none, _ => return some a
    | some bt, none => return some { a with top := some bt }
    | some bt, some at_ =>
      if at_.depths.2 > bt.depths.1 then return none
      QemM.link at_.ends.2 bt.ends.1
      return some { a with top := some ⟨(at_.ends.1, bt.ends.2), (at_.depths.1, bt.depths.2)⟩ }

def mergePlanarity : Option Planarity → Option Planarity → QemM (Option Planarity)
  | none, _ => return none
  | some _, none => return none
  | some a, some b => do
    match ← mergeSide a.sides.1 b.sides.1 with
    | none => return none
    | some s0 =>
      match ← mergeSide a.sides.2 b.sides.2 with
      | none => return none
      | some s1 => return some ⟨(s0, s1)⟩

/-- Close the back edges of a side into the ear (they all return to the current depth). -/
def closeSide (s : PlSide) : QemM PlSide := do
  match s.top, s.bot with
  | some t, some b =>
    QemM.link b.2 t.ends.2
    return { bot := some (b.1, t.ends.1), top := none }
  | _, _ => return s

/-- Peel off the back edges returning to depth `d` from the inner end of a side. -/
def pruneSide (edgeTopDepths : Array Nat) (d : Nat) : Nat → PlSide → QemM PlSide
  | 0, s => return s
  | fuel + 1, s => do
    match s.top, s.bot with
    | some t, some b =>
      if t.depths.2 != d then return s
      QemM.link b.2 t.ends.2
      let nb := QE.flipDir t.ends.2
      let s := { s with bot := some (b.1, nb) }
      match ← QemM.take nb with
      | some x =>
        QemM.clear x
        pruneSide edgeTopDepths d fuel { s with top := some ⟨(t.ends.1, x), (t.depths.1, edgeTopDepths[QE.edge x]!)⟩ }
      | none => return { s with top := none }
    | _, _ => return s

/-- Leaving a child: fold side 1 (which holds only the `lowval` returns) around onto side 0. -/
def foldPlanarity (lowval : Nat) : Option Planarity → QemM (Option Planarity)
  | none => return none
  | some p => do
    let s0 := p.sides.1
    let s1 := p.sides.2
    match s0.bot, s1.bot with
    | some b0, some b1 =>
      QemM.link b0.1 b1.1
      let s0 := { s0 with bot := some (b1.2, b0.2) }
      match s1.top with
      | none => return some ⟨(s0, default)⟩
      | some t1 =>
        if t1.depths.2 != lowval then return none
        match s0.top with
        | some t0 =>
          QemM.link t0.ends.1 t1.ends.1
          return some ⟨({ s0 with top := some { t0 with ends := (t1.ends.2, t0.ends.2) } }, default)⟩
        | none => return some ⟨(s0, default)⟩
    | _, _ => return some p

namespace PlanarWalkM

def liftW (m : WalkM α) : PlanarWalkM α := fun s =>
  let (a, b) := m s.base
  (a, { s with base := b })

def liftQem (m : QemM α) : PlanarWalkM α := fun s =>
  let (a, q) := m s.aux.qem
  (a, { s with aux := { s.aux with qem := q } })

def getAux : PlanarWalkM PlanarAux := do return (← get).aux
def modifyAux (f : PlanarAux → PlanarAux) : PlanarWalkM Unit := modify fun s => { s with aux := f s.aux }

def curPl : PlanarWalkM PlEntry := do return (← getAux).plStack.head!
def nxtPl : PlanarWalkM PlEntry := do return (← getAux).plStack.tail.head!
def modifyCurPl (f : PlEntry → PlEntry) : PlanarWalkM Unit :=
  modifyAux fun a => { a with plStack := match a.plStack with | x :: rest => f x :: rest | [] => [] }
def modifyNxtPl (f : PlEntry → PlEntry) : PlanarWalkM Unit :=
  modifyAux fun a => { a with plStack := match a.plStack with | x :: y :: rest => x :: f y :: rest | l => l }
def pushPl (e : PlEntry) : PlanarWalkM Unit := modifyAux fun a => { a with plStack := e :: a.plStack }
def popPl : PlanarWalkM PlEntry := do
  let a ← getAux
  modifyAux fun a => { a with plStack := a.plStack.tail }
  return a.plStack.head!
/-- Entry `i` counted from the bottom of the stack. -/
def plAt (i : Nat) : PlanarWalkM PlEntry := do
  let a ← getAux
  return a.plStack[a.plStack.length - 1 - i]!
def modifyPlAt (i : Nat) (f : PlEntry → PlEntry) : PlanarWalkM Unit :=
  modifyAux fun a => { a with plStack := a.plStack.modify (a.plStack.length - 1 - i) f }
def topDepthAt (i : Nat) : PlanarWalkM Nat := do
  let s ← liftW get
  return (s.tstack[s.tstack.length - 1 - i]!).topDepth

def nodeIdx (item : ItemId) : PlanarWalkM Nat := do
  let g := (← liftW get).g
  return item - (1 + g.nv + g.ne)
def setNodePlanarity (item : ItemId) (p : NodePlanarity) : PlanarWalkM Unit := do
  let i ← nodeIdx item
  modifyAux fun a => { a with nodePlanarity := a.nodePlanarity.set! i p }
def getNodePlanarity (item : ItemId) : PlanarWalkM NodePlanarity := do
  let i ← nodeIdx item
  return (← getAux).nodePlanarity[i]!
def getItemFlips (item : ItemId) : PlanarWalkM (List Bool) := do return (← getAux).itemFlips[item]!
def modifyItemFlips (item : ItemId) (f : List Bool → List Bool) : PlanarWalkM Unit :=
  modifyAux fun a => { a with itemFlips := a.itemFlips.modify item f }
def setItemFlips (item : ItemId) (fl : List Bool) : PlanarWalkM Unit := modifyItemFlips item fun _ => fl

def allocItem (type : NodeType) : PlanarWalkM ItemId := do
  let item ← liftW (WalkM.allocItem type)
  modifyAux fun a => { a with nodePlanarity := a.nodePlanarity.push .unset, itemFlips := a.itemFlips.push [] }
  return item

def makeEdgePlanarity (item topDepth : Nat) (isTree : Bool) : PlanarWalkM Planarity := do
  let s ← liftW get
  let ve := item - (1 + s.g.nv)
  modifyAux fun a => { a with edgeTopDepths := a.edgeTopDepths.set! ve topDepth }
  return Spqr.makeEdgePlanarity ve topDepth s.stackDir[topDepth]! isTree

def pushTstack (vStart topDepth : Nat) (item : ItemId) (pl : Option Planarity) : PlanarWalkM Unit := do
  liftW (WalkM.pushTstack vStart topDepth item)
  let topDir ← liftW (WalkM.stackDir topDepth)
  pushPl ⟨setSides topDir [false] [], pl⟩
def pushVertTstack (v topDepth : Nat) : PlanarWalkM Unit := pushTstack v topDepth (vertItem v) (some default)
def pushEdgeTstack (vStart topDepth e : Nat) (isTree : Bool) : PlanarWalkM Unit := do
  let item := edgeItem (← liftW get).g e
  let pl ← makeEdgePlanarity item topDepth isTree
  pushTstack vStart topDepth item (some pl)

def mergeTstackTops : PlanarWalkM Unit := do
  liftW WalkM.mergeTstackTops
  let b ← popPl
  let a ← curPl
  let pl ← liftQem (mergePlanarity a.pl b.pl)
  modifyCurPl fun a => { flips := (b.flips.1 ++ a.flips.1, a.flips.2 ++ b.flips.2), pl := pl }

def maybeUnwrapNxt (type : NodeType) (isTree : Bool) : PlanarWalkM ItemId := do
  if type == .R || (← liftW get).ternarize then return ← allocItem type
  let t ← liftW WalkM.nxt
  let topDir ← liftW (WalkM.stackDir t.topDepth)
  let item := (getSide t.spans topDir).head!
  let it ← liftW (WalkM.getItem item)
  if it.type == type then
    liftW (WalkM.modifyNxt fun t => { t with spans := setSides topDir it.ch [] })
    let flips ← getItemFlips item
    let np ← getNodePlanarity item
    modifyNxtPl fun e =>
      { flips := setSides topDir flips [],
        pl := match np, e.pl with
          | .planar m, some p => some (unwrapPlanarity isTree topDir m p)
          | _, pl => pl }
    return item
  else
    allocItem type

def finishTstackTop (item : ItemId) (isTree : Bool) : PlanarWalkM Unit := do
  let t ← liftW WalkM.cur
  let topDir ← liftW (WalkM.stackDir t.topDepth)
  let e ← curPl
  setNodePlanarity item (match e.pl with
    | some p => .planar (finishMatches isTree topDir p)
    | none => .nonplanar)
  liftW (WalkM.finishTstackTop item)
  setItemFlips item (getSide e.flips topDir)
  let pl ← makeEdgePlanarity item t.topDepth isTree
  modifyCurPl fun _ => ⟨setSides topDir [false] [], some pl⟩

def loop : Nat → PlanarWalkM Bool → PlanarWalkM Unit → PlanarWalkM Unit
  | 0, _, _ => pure ()
  | fuel + 1, cond, body => do
    if ← cond then
      body
      loop fuel cond body

/-- `loop fuel cond body`, where the body also sees whether this is the first iteration. -/
def loopFirst : Nat → PlanarWalkM Bool → (Bool → PlanarWalkM Unit) → Bool → PlanarWalkM Unit
  | 0, _, _, _ => pure ()
  | fuel + 1, cond, body, first => do
    if ← cond then
      body first
      loopFirst fuel cond body false

/-- After closing an S/P/R ear at depth `d`: all its remaining back edges return to `d`. -/
def closeBackedges : PlanarWalkM Unit := do
  let e ← curPl
  if let some p := e.pl then
    let s0 ← liftQem (closeSide p.sides.1)
    let s1 ← liftQem (closeSide p.sides.2)
    modifyCurPl fun e => { e with pl := some ⟨(s0, s1)⟩ }

/-- Before merging `nxt` into `cur` (ears whose first back edge is after a back edge of `cur`):
put the returns to `d` on side 1. -/
def flipBeforeMerge (d fo : Nat) (isSingle : Bool) : PlanarWalkM Unit := do
  let n ← liftW WalkM.nxt
  let c ← liftW WalkM.cur
  if n.firstIdx > fo then
    if n.topDepth == d then modifyNxtPl flipEntry
  else if !isSingle then
    let ne ← nxtPl
    if let some p := ne.pl then
      if (p.sides.1.top.map (·.depths.2)) == some d then
        if c.topDepth < n.topDepth then modifyNxtPl flipEntry else modifyCurPl flipEntry

/-- Peel off the back edges of `cur` returning to `d`. -/
def pruneBackedges (d : Nat) : PlanarWalkM Unit := do
  let e ← curPl
  if let some p := e.pl then
    let a ← getAux
    let s0 ← liftQem (pruneSide a.edgeTopDepths d a.qem.size p.sides.1)
    let s1 ← liftQem (pruneSide a.edgeTopDepths d a.qem.size p.sides.2)
    modifyCurPl fun e => { e with pl := some ⟨(s0, s1)⟩ }

/-- Before merging everything above the vertex ear: the `lowval` returns go on side 1. -/
def flipForLowval (lowval origTstack : Nat) : PlanarWalkM Unit := do
  let e ← plAt (origTstack + 2)
  if let some p := e.pl then
    if (p.sides.1.top.map (·.depths.2)) == some lowval then modifyPlAt (origTstack + 2) flipEntry
  let n ← liftW WalkM.tstackSize
  (List.range' (origTstack + 3) (n - (origTstack + 3))).forM fun i => do
    if (← topDepthAt i) == lowval then modifyPlAt i flipEntry

/-- Leaving a child along an edge of direction `edgeDir`: fold both sides onto `!edgeDir`. -/
def foldSides (edgeDir : Bool) (lowval : Nat) : PlanarWalkM Unit := do
  let e ← curPl
  let pl ← liftQem (foldPlanarity lowval e.pl)
  modifyCurPl fun e => { flips := setSides (!edgeDir) (e.flips.1 ++ e.flips.2) [], pl := pl }

end PlanarWalkM

open PlanarWalkM in
/-- `finishEdge` with planarity bookkeeping: the same steps, each base step lifted by `liftW`,
with the planarity-only steps in between. -/
def planarFinishEdge (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool) : PlanarWalkM Bool := do
  let g := (← liftW get).g
  let nxtV := o.dest
  let e := o.e
  let lowval := o.cls.lowval d
  let isTree := o.cls.isTree
  let isType1 := o.cls.isType1
  let edgeDir ← liftW (WalkM.stackDir d)
  let qItem := edgeItem g e

  if lowval ≥ d then
    liftW (WalkM.modifyItem qItem fun it => { it with vs := (some curV, none) })
    liftW (modify fun s => { s with totBlocks := s.totBlocks + 1 })
    if isTree then
      if lowval == d + 1 then
        let item ← allocItem .I
        let vs ← liftW (WalkM.makeVs nxtV d)
        liftW (WalkM.modifyItem item fun it => { it with vs := vs })
        let t ← liftW WalkM.popTstack
        let tp ← popPl
        liftW (WalkM.modifyItem qItem fun it => { it with ch := item :: t.spans.2 })
        setItemFlips qItem (false :: tp.flips.2)
      else
        let backedge ← liftW WalkM.popTstack
        let bp ← popPl
        let t ← liftW WalkM.popTstack
        let tp ← popPl
        liftW (WalkM.modifyItem qItem fun it => { it with ch := backedge.spans.1 ++ t.spans.2 })
        setItemFlips qItem (bp.flips.1 ++ tp.flips.2)
    else
      liftW (modify fun s => { s with totSelfLoops := s.totSelfLoops + 1 })
      let item ← allocItem .O
      liftW (WalkM.modifyItem item fun it => { it with vs := (some curV, none) })
      liftW (WalkM.modifyItem qItem fun it => { it with ch := [item] })
      setItemFlips qItem [false]
    liftW (WalkM.modifyItem (vertItem curV) fun it => { it with ch := it.ch ++ [qItem] })
    modifyItemFlips (vertItem curV) (· ++ [false])
    return hasVert

  let vs ← liftW (WalkM.makeVs nxtV d)
  liftW (WalkM.modifyItem qItem fun it => { it with vs := vs })

  let mut isSingle := true
  if isTree then
    pushEdgeTstack nxtV d e true
    loop (← liftW WalkM.tstackSize) (do return (← liftW WalkM.tstackSize) ≥ 2 && (← liftW WalkM.nxt).topDepth ≥ d) do
      let type ← do
        if (← liftW WalkM.nxt).topDepth > d then
          liftW (WalkM.setStackDir (← liftW WalkM.nxt).topDepth edgeDir)
          mergeTstackTops
          pure NodeType.S
        else if (← liftW WalkM.nxt).vStart == (← liftW WalkM.cur).vStart then pure .P
        else pure .R
      let item ← maybeUnwrapNxt type (type == .S)
      mergeTstackTops
      closeBackedges
      finishTstackTop item true
    let fo := (← liftW get).firstOccurrence[d]!
    if (← liftW WalkM.cur).firstIdx > fo then
      loopFirst (← liftW WalkM.tstackSize) (do return (← liftW WalkM.cur).firstIdx > fo) (fun first => do
        flipBeforeMerge d fo (isSingle && first)
        mergeTstackTops) true
      isSingle := false
      pruneBackedges d
    if hasVert then
      if !isType1 then
        flipForLowval lowval origTstack
        loop (← liftW WalkM.tstackSize) (do return (← liftW WalkM.tstackSize) > origTstack + 3) mergeTstackTops
        isSingle := false
      let item ← if isType1 then some <$> maybeUnwrapNxt (if isSingle then .S else .R) false else pure none
      mergeTstackTops
      mergeTstackTops
      liftW (WalkM.modifyCur fun t => { t with vStart := curV, spans := setSides (!edgeDir) (t.spans.1 ++ t.spans.2) [] })
      foldSides edgeDir lowval
      if let some item := item then
        finishTstackTop item false
        isSingle := true
  else
    pushEdgeTstack curV lowval e false
    liftW (modify fun s =>
      { s with firstOccurrence := s.firstOccurrence.modify lowval (min · s.nxtEdgeIdx),
               nxtEdgeIdx := s.nxtEdgeIdx + 1 })

  if isType1 && (← liftW WalkM.tstackSize) ≥ 2 && (← liftW WalkM.nxt).vStart == curV && (← liftW WalkM.nxt).topDepth == lowval then
    let item ← maybeUnwrapNxt .P false
    mergeTstackTops
    finishTstackTop item false

  if !hasVert then
    pushVertTstack curV d
    if !isSingle then mergeTstackTops
    return true
  return hasVert

open PlanarWalkM in
mutual
def planarWalkTree (t : DfsTree) (d : Nat) : PlanarWalkM Unit := do
  match t with
  | .node v outs =>
    liftW (modify fun s => { s with stackVerts := s.stackVerts.set! d v })
    let hasVert ← planarWalkOuts v d outs false
    unless hasVert do
      liftW (WalkM.setStackDir d true)
      pushVertTstack v d

def planarWalkOuts (v d : Nat) (outs : List DfsOut) (hasVert : Bool) : PlanarWalkM Bool := do
  match outs with
  | [] => return hasVert
  | o :: rest =>
    let hasVert ← planarWalkOut v d o hasVert
    planarWalkOuts v d rest hasVert

def planarWalkOut (v d : Nat) (o : DfsOut) (hasVert : Bool) : PlanarWalkM Bool := do
  let lowval := o.cls.lowval d
  let lowDir ← liftW (WalkM.stackDir lowval)
  liftW (WalkM.setStackDir d (if lowval ≥ d then false else !lowDir))
  let hasVert ← do
    if !hasVert && lowval < d && o.cls.isType1 then
      pushVertTstack v d
      pure true
    else pure hasVert
  let origTstack ← liftW WalkM.tstackSize
  match o with
  | .tree _ _ child =>
    liftW (modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne })
    planarWalkTree child (d + 1)
  | .back .. => pure ()
  planarFinishEdge v d o origTstack hasVert
end

open PlanarWalkM in
def planarWalkForest (forest : List DfsTree) : PlanarWalkM Unit :=
  forest.forM fun t => do
    planarWalkTree t 0
    let top ← liftW WalkM.popTstack
    let tp ← popPl
    liftW (WalkM.modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 })
    modifyItemFlips rootItem (· ++ tp.flips.2)

/-- Phase 2 entry point, planar variant. -/
def Graph.planarWalk (g : Graph) (ternarize : Bool) (forest : List DfsTree) : PlanarWalkState :=
  (planarWalkForest forest).run (PlanarWalkState.init g ternarize) |>.2

end Spqr
