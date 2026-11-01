import Spqr.Relabel
import Spqr.Refine

/-!
# `relabel` with a ghost visiting order

`Ghost.relabel` is `relabel` plus an `order : Array ItemId` recording the items in the order they
are numbered, so that the preorder index of item `i` is `order.idxOf i` — a definition rather
than an existential. Projecting `order` away gives back `relabel` (`sim_relabel`), hence
`relabelTree g items = SpqrTree.ofRelabelState g (Ghost.relabelRun g items).toRelabelState`.
-/

namespace Spqr

/-- `relabelTree`'s packaging of a final `RelabelState`. -/
def SpqrTree.ofRelabelState (g : Graph) (s : RelabelState) : SpqrTree :=
  let nodeVerts := s.nodeVerts.map fun nv => { nv with vert := (s.vertIndex[nv.vert]!).getD 0 }
  { nv := g.nv, ne := g.ne, vertIndex := s.vertIndex, edgeIndex := s.edgeIndex, edgeFlipped := s.edgeFlipped,
    par := s.par, subtreeEnd := s.subtreeEnd, types := s.types, origId := s.origId,
    chBounds := s.chBounds, chDat := s.chDat, nodeVerts := nodeVerts, nvBounds := s.nvBounds,
    vertParNv := s.vertParNv, nodeEdges := s.nodeEdges, neBounds := s.neBounds,
    adjBounds := s.adjBounds, adjDat := s.adjDat }

theorem relabelTree_eq_ofRelabelState (g : Graph) (items : Array Item) :
    relabelTree g items =
      SpqrTree.ofRelabelState g ((relabel items.size rootItem none none none).run (RelabelState.init g items) |>.2) :=
  rfl

namespace Ghost

structure RelabelState extends Spqr.RelabelState where
  /-- The items in the order they were numbered. -/
  order : Array ItemId := #[]

abbrev RelabelM := StateM RelabelState

namespace RelabelM

def item (i : ItemId) : RelabelM Item := do return (← get).items[i]!

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

open Spqr.Ghost.RelabelM in
/-- `Spqr.relabel`, also pushing `cur` onto `order` when it is numbered. -/
def relabel : Nat → ItemId → Option Nat → Option Nat → Option Nat → RelabelM Unit
  | 0, _, _, _, _ => pure ()
  | fuel + 1, cur, parent, parNv, capTwin => do
    let g := (← get).g
    let it ← item cur
    let curIdx := (← get).types.size
    modify fun s => { s with
      order := s.order.push cur,
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

def RelabelState.init (g : Graph) (items : Array Item) : RelabelState :=
  { Spqr.RelabelState.init g items with }

/-- The ghost state after relabeling the whole item tree. -/
def relabelRun (g : Graph) (items : Array Item) : RelabelState :=
  (relabel items.size rootItem none none none).run (RelabelState.init g items) |>.2

/-- The preorder index of item `i` (`items.size` if `i` was never numbered). -/
def relabelIdx (g : Graph) (items : Array Item) (i : ItemId) : Nat :=
  (relabelRun g items).order.toList.idxOf i

/-! ## Refinement -/

namespace RelabelM

theorem sim_item (i : ItemId) :
    Refine.Sim RelabelState.toRelabelState Eq (item i) (Spqr.RelabelM.item i) := fun _ => ⟨rfl, rfl⟩

theorem sim_orderedChildren (it : Item) (n : Nat) :
    Refine.Sim RelabelState.toRelabelState Eq (orderedChildren it n) (Spqr.RelabelM.orderedChildren it n) :=
  fun s => by
  unfold orderedChildren Spqr.RelabelM.orderedChildren
  by_cases h : (it.type != .R) = true
  · simp only [h, ↓reduceIte]
    refine ⟨?_, ?_⟩ <;> first | rfl | trivial
  · simp only [h, Bool.false_eq_true, ↓reduceIte, run_simp]
    refine ⟨?_, ?_⟩ <;> first | rfl | trivial

end RelabelM

macro_rules
  | `(tactic| sim_leaf) => `(tactic| first | apply RelabelM.sim_item | apply RelabelM.sim_orderedChildren)

set_option maxHeartbeats 4000000 in
theorem sim_relabel : ∀ (fuel : Nat) (cur : ItemId) (p pn ct : Option Nat),
    Refine.Sim RelabelState.toRelabelState Eq (relabel fuel cur p pn ct) (Spqr.relabel fuel cur p pn ct)
  | 0, _, _, _, _ => by unfold relabel Spqr.relabel; sim_auto
  | fuel + 1, cur, p, pn, ct => by
    have ih := sim_relabel fuel
    unfold relabel Spqr.relabel
    sim_auto

theorem relabelRun_toSlow (g : Graph) (items : Array Item) :
    (relabelRun g items).toRelabelState =
      ((Spqr.relabel items.size rootItem none none none).run (Spqr.RelabelState.init g items) |>.2) :=
  (sim_relabel items.size rootItem none none none (RelabelState.init g items)).1

theorem relabelTree_eq (g : Graph) (items : Array Item) :
    relabelTree g items = SpqrTree.ofRelabelState g (relabelRun g items).toRelabelState := by
  rw [relabelTree_eq_ofRelabelState, relabelRun_toSlow]

end Ghost

end Spqr
