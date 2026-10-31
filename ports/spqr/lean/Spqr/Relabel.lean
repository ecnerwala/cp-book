import Spqr.Items

/-!
# Phase 3: relabel the item tree in preorder

The item tree built by the walk is numbered in preorder, and each SPQR node's skeleton is laid
out: its vertices (`NodeVert`), its skeleton edges (`NodeEdge`, each a cap edge shared with the
parent or a virtual edge shared with a child, both ends linked by `twin`), and a per-vertex
adjacency CSR in bracket order.
-/

namespace Spqr

structure NodeVert where
  node : Nat
  vert : Nat
deriving Repr, Inhabited, DecidableEq

structure NodeEdge where
  node : Nat
  twin : Option Nat := none
  /-- Endpoints, as node-vertex ids. -/
  nvs : Nat × Nat
deriving Repr, Inhabited, DecidableEq

structure NodeAdj where
  ne : Nat
  destNv : Nat
deriving Repr, Inhabited, DecidableEq

structure SpqrTree where
  nv : Nat
  ne : Nat
  vertIndex : Array (Option Nat)
  edgeIndex : Array (Option Nat)
  edgeFlipped : Array Bool
  par : Array (Option Nat)
  subtreeEnd : Array Nat
  types : Array NodeType
  origId : Array (Option Nat)
  chBounds : Array Nat
  chDat : Array Nat
  nodeVerts : Array NodeVert
  nvBounds : Array Nat
  vertParNv : Array (Option Nat)
  nodeEdges : Array NodeEdge
  neBounds : Array Nat
  adjBounds : Array Nat
  adjDat : Array NodeAdj
deriving Repr

structure RelabelState where
  g : Graph
  items : Array Item
  vertIndex : Array (Option Nat)
  edgeIndex : Array (Option Nat)
  edgeFlipped : Array Bool
  par : Array (Option Nat) := #[]
  subtreeEnd : Array Nat := #[]
  types : Array NodeType := #[]
  origId : Array (Option Nat) := #[]
  chBounds : Array Nat := #[0]
  chDat : Array Nat := #[]
  nodeVerts : Array NodeVert := #[]
  nvBounds : Array Nat := #[0]
  vertParNv : Array (Option Nat) := #[]
  nodeEdges : Array NodeEdge := #[]
  neBounds : Array Nat := #[0]
  adjBounds : Array Nat := #[0]
  adjDat : Array NodeAdj := #[]
  /-- Scratch: position of each vertex inside the current R node. -/
  vertPos : Array Nat

abbrev RelabelM := StateM RelabelState

namespace NodeType
def isNode : NodeType → Bool
  | F | V => false
  | _ => true
end NodeType

/-- Skeleton layout of one node, relative to its first node-vertex `nvSt` and node-edge `neSt`:
`adjBounds[i]` for `i ∈ [1, 2 * nVerts]` is the global bound `adjBounds[2 * nvSt + i]`, and
`adjDat[j]` is the global slot `2 * neSt + j`. -/
structure Layout where
  edges : Array NodeEdge
  adjBounds : Array Nat
  adjDat : Array NodeAdj

namespace Layout

def empty (nVerts nEdges : Nat) : Layout :=
  ⟨Array.replicate nEdges default, Array.replicate (2 * nVerts + 1) 0, Array.replicate (2 * nEdges) default⟩

/-- Place edge `ne` between node-verts `nvs`, writing its two adjacency entries at slots `nds`. -/
def setNe (l : Layout) (neSt : Nat) (node ne : Nat) (nvs : Nat × Nat) (nds : Nat × Nat) : Layout :=
  { l with
    edges := l.edges.set! (ne - neSt) ⟨node, none, nvs⟩
    adjDat := (l.adjDat.set! (nds.1 - 2 * neSt) ⟨ne, nvs.2⟩).set! (nds.2 - 2 * neSt) ⟨ne, nvs.1⟩ }

end Layout

/-- Lay out the skeleton of a node of type `type`, given its node-verts `[nvSt, nvEn)`, node-edges
`[neSt, neEn)`, and (for R nodes) its ordered children with their node-vert endpoints. -/
def layoutNode (type : NodeType) (node nvSt nvEn neSt neEn : Nat)
    (edgeChildren : List (Nat × Nat)) : Layout := Id.run do
  let nVerts := nvEn - nvSt
  let nEdges := neEn - neSt
  let mut l := Layout.empty nVerts nEdges
  -- Local index of global adjacency bound `i`.
  let b (i : Nat) := i - 2 * nvSt
  match type with
  | .F =>
    for i in [2 * nvSt + 1 : 2 * nvEn + 1] do
      l := { l with adjBounds := l.adjBounds.set! (b i) (2 * neSt) }
  | .V => pure ()
  | _ =>
    if nVerts == 1 then
      -- Q self-loop or O node
      l := { l with adjBounds := (l.adjBounds.set! 1 (2 * neSt + 1)).set! 2 (2 * neSt + 2) }
      l := l.setNe neSt node neSt (nvSt, nvSt) (2 * neSt + 1, 2 * neSt)
    else if type == .Q || type == .I then
      l := { l with adjBounds := (((l.adjBounds.set! 1 (2 * neSt)).set! 2 (2 * neSt + 1)).set! 3 (2 * neSt + 2)).set! 4 (2 * neSt + 2) }
      l := l.setNe neSt node neSt (nvSt, nvSt + 1) (2 * neSt, 2 * neSt + 1)
    else if type == .P then
      l := { l with adjBounds := (((l.adjBounds.set! 1 (2 * neSt)).set! 2 (2 * neSt + nEdges)).set! 3 (2 * neSt + 2 * nEdges)).set! 4 (2 * neSt + 2 * nEdges) }
      for k in [0 : nEdges] do
        l := l.setNe neSt node (neSt + k) (nvSt, nvSt + 1) (2 * neSt + k, 2 * neEn - 1 - k)
    else if type == .S then
      for i in [2 * nvSt + 1 : 2 * nvEn + 1] do
        l := { l with adjBounds := l.adjBounds.set! (b i) (b i + 2 * neSt) }
      l := { l with adjBounds := (l.adjBounds.modify 1 (· - 1)).modify (b (2 * nvEn - 1)) (· + 1) }
      l := l.setNe neSt node neSt (nvSt, nvEn - 1) (2 * neSt, 2 * neEn - 1)
      for i in [1 : nEdges] do
        let ne := neSt + i
        l := l.setNe neSt node ne (nvSt + i - 1, nvSt + i) (2 * ne - 1, 2 * ne)
    else
      -- R: count the entries of each half-row, prefix-sum, then fill in reverse child order.
      let inc (l : Layout) (i : Nat) : Layout := { l with adjBounds := l.adjBounds.modify (b i) (· + 1) }
      l := inc (inc l (2 * nvSt + 2)) (2 * nvEn - 1)
      for (a, c) in edgeChildren do
        l := inc (inc l (2 * a + 2)) (2 * c + 1)
      let mut off := 2 * neSt
      for i in [2 * nvSt + 1 : 2 * nvEn + 1] do
        let cnt := l.adjBounds[b i]!
        l := { l with adjBounds := l.adjBounds.set! (b i) off }
        off := off + cnt
      l := inc l (2 * nvSt + 2)
      let mut nxtNe := neEn
      for (a, c) in edgeChildren.reverse do
        nxtNe := nxtNe - 1
        let s0 := l.adjBounds[b (2 * a + 2)]!
        let s1 := l.adjBounds[b (2 * c + 1)]!
        l := inc (inc l (2 * a + 2)) (2 * c + 1)
        l := l.setNe neSt node nxtNe (a, c) (s0, s1)
      l := l.setNe neSt node neSt (nvSt, nvEn - 1) (2 * neSt, 2 * neEn - 1)
      l := inc l (2 * nvEn - 1)
  return l

namespace RelabelM

def item (i : ItemId) : RelabelM Item := do return (← get).items[i]!

/-- The children of `item` in output order: for R nodes, stably sorted by the sum of their
endpoint positions, so that the adjacency lists come out in bracket order. -/
def orderedChildren (it : Item) (nvSt : Nat) : RelabelM (List ItemId) := do
  if it.type != .R then return it.ch
  let s ← get
  let loc (c : ItemId) : Nat :=
    if c < 1 + s.g.nv then 2 * (s.vertPos[c - 1]! - nvSt)
    else
      let cvs := s.items[c]!.vs
      (s.vertPos[cvs.1.getD 0]! - nvSt) + (s.vertPos[cvs.2.getD 0]! - nvSt)
  return it.ch.mergeSort fun a b => loc a ≤ loc b

end RelabelM

open RelabelM in
/-- Number `cur` as the next preorder index, lay out its skeleton, then recurse into its children.
`parent` / `parNv` / `capTwin` are what the parent recorded for us. `fuel` bounds the recursion
depth (any `fuel ≥ items.size` suffices). -/
def relabel : Nat → ItemId → Option Nat → Option Nat → Option Nat → RelabelM Unit
  | 0, _, _, _, _ => pure ()
  | fuel + 1, cur, parent, parNv, capTwin => do
    let g := (← get).g
    let it ← item cur
    let curIdx := (← get).types.size
    modify fun s => { s with
      types := s.types.push it.type, par := s.par.push parent, subtreeEnd := s.subtreeEnd.push 0,
      vertParNv := s.vertParNv.push parNv, origId := s.origId.push none }
    match it.type with
    | .V =>
      let v := cur - 1
      modify fun s => { s with origId := s.origId.set! curIdx (some v), vertIndex := s.vertIndex.set! v (some curIdx) }
    | .Q =>
      let e := cur - 1 - g.nv
      let flipped := it.vs.1 != some (g.edges[e]!).1
      modify fun s => { s with
        origId := s.origId.set! curIdx (some e)
        edgeIndex := s.edgeIndex.set! e (some curIdx)
        edgeFlipped := s.edgeFlipped.set! e flipped }
    | _ => pure ()

    -- Node-verts: the first cap endpoint, the vertex children, the second cap endpoint.
    let nvSt := (← get).nodeVerts.size
    let chSt := (← get).chDat.size
    let vertChildren := it.ch.filter (· < 1 + g.nv)
    let nodeVerts : List NodeVert :=
      (it.vs.1.toList ++ vertChildren.map (· - 1) ++ it.vs.2.toList).map (⟨curIdx, ·⟩)
    modify fun s => { s with nodeVerts := s.nodeVerts ++ nodeVerts.toArray }
    let nvEn := nvSt + nodeVerts.length
    if it.type == .R then
      modify fun s => { s with vertPos := (nodeVerts.zipIdx nvSt).foldl (init := s.vertPos) fun a (nv, pos) => a.set! nv.vert pos }
    let children ← orderedChildren it nvSt
    modify fun s => { s with chDat := s.chDat ++ children.toArray }
    let chEn := chSt + children.length

    let isNode := it.type.isNode
    let hasCap := isNode && !(it.type == .Q && !children.isEmpty)
    let nEdges := (if isNode then children.countP (· ≥ 1 + g.nv) else 0) + (if hasCap then 1 else 0)
    let neSt := (← get).nodeEdges.size
    let neEn := neSt + nEdges

    let vertPos := (← get).vertPos
    let items := (← get).items
    let edgeChildren := (children.filter (· ≥ 1 + g.nv)).map fun c =>
      let cvs := items[c]!.vs
      (vertPos[cvs.1.getD 0]!, vertPos[cvs.2.getD 0]!)
    let l := layoutNode it.type curIdx nvSt nvEn neSt neEn edgeChildren
    modify fun s => { s with
      nodeEdges := s.nodeEdges ++ l.edges,
      adjDat := s.adjDat ++ l.adjDat,
      adjBounds := s.adjBounds ++ (l.adjBounds.extract 1 (2 * (nvEn - nvSt) + 1)),
      chBounds := s.chBounds.push chEn, nvBounds := s.nvBounds.push nvEn, neBounds := s.neBounds.push neEn }
    if hasCap then
      modify fun s => { s with nodeEdges := s.nodeEdges.modify neSt fun ne => { ne with twin := capTwin } }

    -- Children, in order.
    let mut curNv := nvSt + (if it.vs.1.isSome then 1 else 0)
    let mut curNe := neSt + (if hasCap then 1 else 0)
    let mut k := 0
    for c in children do
      let nxtIdx := (← get).types.size
      let nxtNe := (← get).nodeEdges.size
      modify fun s => { s with chDat := s.chDat.set! (chSt + k) nxtIdx }
      k := k + 1
      if c < 1 + g.nv then
        relabel fuel c (some curIdx) (some curNv) none
        curNv := curNv + 1
      else if isNode then
        modify fun s => { s with nodeEdges := s.nodeEdges.modify curNe fun ne => { ne with twin := some nxtNe } }
        relabel fuel c (some curIdx) none (some curNe)
        curNe := curNe + 1
      else
        relabel fuel c (some curIdx) none none
    modify fun s => { s with subtreeEnd := s.subtreeEnd.set! curIdx s.types.size }

def RelabelState.init (g : Graph) (items : Array Item) : RelabelState where
  g := g
  items := items
  vertIndex := Array.replicate g.nv none
  edgeIndex := Array.replicate g.ne none
  edgeFlipped := Array.replicate g.ne false
  vertPos := Array.replicate g.nv 0

/-- Phase 3 entry point. -/
def relabelTree (g : Graph) (items : Array Item) : SpqrTree :=
  let s := (relabel items.size rootItem none none none).run (RelabelState.init g items) |>.2
  -- Node-verts refer to vertices by the preorder index of their V item.
  let nodeVerts := s.nodeVerts.map fun nv => { nv with vert := (s.vertIndex[nv.vert]!).getD 0 }
  { nv := g.nv, ne := g.ne, vertIndex := s.vertIndex, edgeIndex := s.edgeIndex, edgeFlipped := s.edgeFlipped,
    par := s.par, subtreeEnd := s.subtreeEnd, types := s.types, origId := s.origId,
    chBounds := s.chBounds, chDat := s.chDat, nodeVerts := nodeVerts, nvBounds := s.nvBounds,
    vertParNv := s.vertParNv, nodeEdges := s.nodeEdges, neBounds := s.neBounds,
    adjBounds := s.adjBounds, adjDat := s.adjDat }

end Spqr
