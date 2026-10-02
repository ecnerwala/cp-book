import Spqr.Relabel
import Spqr.DfsFast
import Spqr.PlanarWalk

/-!
# Phase 3, planar variant: relabel with local rotation systems

The planar relabel is `relabel` (`Spqr.Relabel`) run on the same `RelabelState`, with extra state
`PlanarRelabelAux` alongside. Every node gets a `nodePlanar` flag and a local rotation system on
its node-edge quarter-edges `4 * ne + 2 * side + dir` (`neRotAdj`): S/P/Q/I/O skeletons have a
fixed layout, R skeletons map the walk's quarter-edge matches through `rotEdgeNe` (and are all
`none` when nonplanar). The base component of the result is the ordinary `SpqrTree`
(`planarRelabel_proj`, in `Spqr.PlanarSpec`).
-/

namespace Spqr

deriving instance DecidableEq for SpqrTree

structure PlanarSpqrTree extends SpqrTree where
  nodePlanar : Array Bool
  /-- Facing quarter-edges of node-edges; `none` for nonplanar nodes. -/
  neRotAdj : Array (Option Nat)
deriving Repr

structure PlanarRelabelAux where
  qem : Qem
  nodePlanarity : Array NodePlanarity
  itemFlips : Array (List Bool)
  nodePlanar : Array Bool := #[]
  neRotAdj : Array (Option Nat) := #[]
  /-- Scratch: node-edge of each virtual edge inside the current R node (`2 ne` is the cap). -/
  rotEdgeNe : Array Nat

structure PlanarRelabelState where
  base : RelabelState
  aux : PlanarRelabelAux

abbrev PlanarRelabelM := StateM PlanarRelabelState

/-- Local rotation system of a node laid out by `layoutNode`, `4 * nEdges` entries indexed by
`4 * (ne - neSt) + 2 * side + dir`. `mapRot ve` gives the four entries of the virtual edge `ve`
of an R node (`2 * ne` for the cap). -/
def layoutRot (type : NodeType) (nVerts neSt neEn : Nat) (edgeVes : List Nat)
    (mapRot : Nat → Array (Option Nat)) (capVe : Nat) : Array (Option Nat) :=
  let nEdges := neEn - neSt
  let rot (a b c d : Nat) : Array (Option Nat) := #[some a, some b, some c, some d]
  match type with
  | .F | .V => #[]
  | _ =>
    if nVerts == 1 then
      rot (4 * neSt + 3) (4 * neSt + 2) (4 * neSt + 1) (4 * neSt + 0)
    else if type == .Q || type == .I then
      rot (4 * neSt + 1) (4 * neSt + 0) (4 * neSt + 3) (4 * neSt + 2)
    else if type == .P then
      (List.range nEdges).foldl (init := #[]) fun acc k =>
        let ne := neSt + k
        let prv := (if ne == neSt then neEn else ne) - 1
        let nxt := if ne + 1 == neEn then neSt else ne + 1
        acc ++ rot (4 * prv + 1) (4 * nxt + 0) (4 * nxt + 3) (4 * prv + 2)
    else if type == .S then
      (List.range (nEdges - 1)).foldl
        (init := rot (4 * (neSt + 1) + 1) (4 * (neSt + 1) + 0) (4 * (neEn - 1) + 3) (4 * (neEn - 1) + 2))
        fun acc i =>
          let ne := neSt + 1 + i
          let (a, b) := if ne - 1 == neSt then (4 * neSt + 1, 4 * neSt + 0) else (4 * (ne - 1) + 3, 4 * (ne - 1) + 2)
          let (c, d) := if ne + 1 == neEn then (4 * neSt + 3, 4 * neSt + 2) else (4 * (ne + 1) + 1, 4 * (ne + 1) + 0)
          acc ++ rot a b c d
    else
      edgeVes.foldl (init := mapRot capVe) fun acc ve => acc ++ mapRot ve

namespace PlanarRelabelM

def liftR (m : RelabelM α) : PlanarRelabelM α := fun s =>
  let (a, b) := m s.base
  (a, { s with base := b })

def getAux : PlanarRelabelM PlanarRelabelAux := do return (← get).aux
def modifyAux (f : PlanarRelabelAux → PlanarRelabelAux) : PlanarRelabelM Unit :=
  modify fun s => { s with aux := f s.aux }

/-- Link the cap's four quarter-edges `8 ne + s` to the node's recorded matches; `false` if the
node was found nonplanar. -/
def setupNode (g : Graph) (type : NodeType) (cur : ItemId) : PlanarRelabelM Bool := do
  match type with
  | .S | .P | .R =>
    match (← getAux).nodePlanarity[cur - (1 + g.nv + g.ne)]! with
    | .planar m =>
      for s in [0 : 4] do
        modifyAux fun a => { a with qem := (a.qem.set! (8 * g.ne + s) (some m[s]!)).set! m[s]! (some (8 * g.ne + s)) }
      return true
    | _ => return false
  | _ => return true

/-- Flipped children of an R node have their virtual edge's direction bits swapped. -/
def applyFlips (g : Graph) (it : Item) (flips : List Bool) : PlanarRelabelM Unit := do
  for (c, flip) in it.ch.zip flips do
    if c ≥ 1 + g.nv && flip then
      let ve := c - (1 + g.nv)
      modifyAux fun a => { a with qem := (a.qem.swapIfInBounds (4 * ve + 0) (4 * ve + 1)).swapIfInBounds (4 * ve + 2) (4 * ve + 3) }

end PlanarRelabelM

open PlanarRelabelM in
/-- `relabel` with planarity: the same steps, each base step lifted by `liftR`. -/
def planarRelabel : Nat → ItemId → Option Nat → Option Nat → Option Nat → PlanarRelabelM Unit
  | 0, _, _, _, _ => pure ()
  | fuel + 1, cur, parent, parNv, capTwin => do
    let g := (← liftR get).g
    let it ← liftR (RelabelM.item cur)
    let curIdx := (← liftR get).types.size
    liftR (modify fun s => { s with
      types := s.types.push it.type, par := s.par.push parent, subtreeEnd := s.subtreeEnd.push 0,
      vertParNv := s.vertParNv.push parNv, origId := s.origId.push none })
    match it.type with
    | .V =>
      let v := cur - 1
      liftR (modify fun s => { s with origId := s.origId.set! curIdx (some v), vertIndex := s.vertIndex.set! v (some curIdx) })
    | .Q =>
      let e := cur - 1 - g.nv
      let flipped := it.vs.1 != some (g.edges[e]!).1
      liftR (modify fun s => { s with
        origId := s.origId.set! curIdx (some e)
        edgeIndex := s.edgeIndex.set! e (some curIdx)
        edgeFlipped := s.edgeFlipped.set! e flipped })
    | _ => pure ()
    let planar ← setupNode g it.type cur
    modifyAux fun a => { a with nodePlanar := a.nodePlanar.push planar }

    let nvSt := (← liftR get).nodeVerts.size
    let chSt := (← liftR get).chDat.size
    let vertChildren := it.ch.filter (· < 1 + g.nv)
    let nodeVerts : List NodeVert :=
      (it.vs.1.toList ++ vertChildren.map (· - 1) ++ it.vs.2.toList).map (⟨curIdx, ·⟩)
    liftR (modify fun s => { s with nodeVerts := s.nodeVerts ++ nodeVerts.toArray })
    let nvEn := nvSt + nodeVerts.length
    if it.type == .R then
      liftR (modify fun s => { s with vertPos := (nodeVerts.zipIdx nvSt).foldl (init := s.vertPos) fun a (nv, pos) => a.set! nv.vert pos })
      applyFlips g it (← getAux).itemFlips[cur]!
    let children ← liftR (RelabelM.orderedChildren it nvSt)
    liftR (modify fun s => { s with chDat := s.chDat ++ children.toArray })
    let chEn := chSt + children.length

    let isNode := it.type.isNode
    let hasCap := isNode && !(it.type == .Q && !children.isEmpty)
    let nEdges := (if isNode then children.countP (· ≥ 1 + g.nv) else 0) + (if hasCap then 1 else 0)
    let neSt := (← liftR get).nodeEdges.size
    let neEn := neSt + nEdges

    let vertPos := (← liftR get).vertPos
    let items := (← liftR get).items
    let edgeChildren := (children.filter (· ≥ 1 + g.nv)).map fun c =>
      let cvs := items[c]!.vs
      (vertPos[cvs.1.getD 0]!, vertPos[cvs.2.getD 0]!)
    let l := layoutNode it.type curIdx nvSt nvEn neSt neEn edgeChildren
    let edgeVes := (children.filter (· ≥ 1 + g.nv)).map (· - (1 + g.nv))
    if it.type == .R then
      modifyAux fun a => { a with rotEdgeNe :=
        ((edgeVes.zipIdx (neSt + 1)).foldl (init := a.rotEdgeNe) fun r (ve, ne) => r.set! ve ne).set! (2 * g.ne) neSt }
    let aux ← getAux
    let mapRot (ve : Nat) : Array (Option Nat) :=
      if !planar then Array.replicate 4 none
      else (Array.range 4).map fun z =>
        (aux.qem[4 * ve + z]!).map fun o => 4 * aux.rotEdgeNe[QE.edge o]! + (o &&& 2) + (1 - z % 2)
    let rot := layoutRot it.type (nvEn - nvSt) neSt neEn edgeVes mapRot (2 * g.ne)
    liftR (modify fun s => { s with
      nodeEdges := s.nodeEdges ++ l.edges,
      adjDat := s.adjDat ++ l.adjDat,
      adjBounds := s.adjBounds ++ (l.adjBounds.extract 1 (2 * (nvEn - nvSt) + 1)),
      chBounds := s.chBounds.push chEn, nvBounds := s.nvBounds.push nvEn, neBounds := s.neBounds.push neEn })
    modifyAux fun a => { a with neRotAdj := a.neRotAdj ++ rot }
    if hasCap then
      liftR (modify fun s => { s with nodeEdges := s.nodeEdges.modify neSt fun ne => { ne with twin := capTwin } })

    let mut curNv := nvSt + (if it.vs.1.isSome then 1 else 0)
    let mut curNe := neSt + (if hasCap then 1 else 0)
    let mut k := 0
    for c in children do
      let nxtIdx := (← liftR get).types.size
      let nxtNe := (← liftR get).nodeEdges.size
      liftR (modify fun s => { s with chDat := s.chDat.set! (chSt + k) nxtIdx })
      k := k + 1
      if c < 1 + g.nv then
        planarRelabel fuel c (some curIdx) (some curNv) none
        curNv := curNv + 1
      else if isNode then
        liftR (modify fun s => { s with nodeEdges := s.nodeEdges.modify curNe fun ne => { ne with twin := some nxtNe } })
        planarRelabel fuel c (some curIdx) none (some curNe)
        curNe := curNe + 1
      else
        planarRelabel fuel c (some curIdx) none none
    liftR (modify fun s => { s with subtreeEnd := s.subtreeEnd.set! curIdx s.types.size })

/-- The output tree read off the final relabel state (the tail of `relabelTree`). -/
def SpqrTree.ofState (g : Graph) (s : RelabelState) : SpqrTree :=
  let nodeVerts := s.nodeVerts.map fun nv => { nv with vert := (s.vertIndex[nv.vert]!).getD 0 }
  { nv := g.nv, ne := g.ne, vertIndex := s.vertIndex, edgeIndex := s.edgeIndex, edgeFlipped := s.edgeFlipped,
    par := s.par, subtreeEnd := s.subtreeEnd, types := s.types, origId := s.origId,
    chBounds := s.chBounds, chDat := s.chDat, nodeVerts := nodeVerts, nvBounds := s.nvBounds,
    vertParNv := s.vertParNv, nodeEdges := s.nodeEdges, neBounds := s.neBounds,
    adjBounds := s.adjBounds, adjDat := s.adjDat }

theorem relabelTree_eq_ofState (g : Graph) (items : Array Item) :
    relabelTree g items = SpqrTree.ofState g ((relabel items.size rootItem none none none).run (RelabelState.init g items) |>.2) := rfl

def PlanarRelabelState.init (g : Graph) (w : PlanarWalkState) : PlanarRelabelState where
  base := RelabelState.init g w.base.items
  aux := { qem := w.aux.qem, nodePlanarity := w.aux.nodePlanarity, itemFlips := w.aux.itemFlips,
           rotEdgeNe := Array.replicate (2 * g.ne + 1) 0 }

/-- Phase 3 entry point, planar variant. -/
def planarRelabelTree (g : Graph) (w : PlanarWalkState) : PlanarSpqrTree :=
  let s := (planarRelabel w.base.items.size rootItem none none none).run (PlanarRelabelState.init g w) |>.2
  { SpqrTree.ofState g s.base with nodePlanar := s.aux.nodePlanar, neRotAdj := s.aux.neRotAdj }

/-- Build the planar SPQR tree of `g` (same arguments as `Graph.spqrTree`). -/
def Graph.planarSpqrTree (g : Graph) (ternarize : Bool := false) (vertOrder edgeOrder : List Nat := []) : PlanarSpqrTree :=
  let forest := g.dfsForestFast vertOrder edgeOrder
  let w := g.planarWalk ternarize forest
  planarRelabelTree g w

end Spqr
