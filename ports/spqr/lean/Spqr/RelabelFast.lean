import Spqr.Relabel
import Spqr.Refine

/-!
# Phase 3, linear version

`Spqr.relabel` with a step counter and R children sorted by precomputed keys. `relabelTreeFast_eq`
shows it computes the same `SpqrTree` as `relabelTree`, by simulation (`ticks` projected away).
-/

namespace Spqr.Fast

structure RelabelState extends Spqr.RelabelState where
  ticks : Nat := 0

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
  return ((it.ch.map fun c => (loc c, c)).mergeSort fun a b => a.1 ≤ b.1).map (·.2)

end RelabelM

open Spqr.Fast.RelabelM in
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
      ticks := s.ticks + 1 + (if it.type.isNode && !(it.type == .Q && !it.ch.isEmpty) then 1 else 0)
        + (if it.vs.1.isSome then 1 else 0) + (if it.vs.2.isSome then 1 else 0),
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
    modify fun s => { s with ticks := s.ticks + 2 * children.length, chDat := s.chDat ++ children.toArray }
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

def RelabelState.init (g : Graph) (items : Array Item) : RelabelState :=
  { Spqr.RelabelState.init g items with }

/-- The phase 3 state after relabeling the whole item tree. -/
def relabelRun (g : Graph) (items : Array Item) : RelabelState :=
  (relabel items.size rootItem none none none).run (RelabelState.init g items) |>.2

/-- Phase 3 entry point. -/
def relabelTreeFast (g : Graph) (items : Array Item) : SpqrTree :=
  let s := relabelRun g items
  -- Node-verts refer to vertices by the preorder index of their V item.
  let nodeVerts := s.nodeVerts.map fun nv => { nv with vert := (s.vertIndex[nv.vert]!).getD 0 }
  { nv := g.nv, ne := g.ne, vertIndex := s.vertIndex, edgeIndex := s.edgeIndex, edgeFlipped := s.edgeFlipped,
    par := s.par, subtreeEnd := s.subtreeEnd, types := s.types, origId := s.origId,
    chBounds := s.chBounds, chDat := s.chDat, nodeVerts := nodeVerts, nvBounds := s.nvBounds,
    vertParNv := s.vertParNv, nodeEdges := s.nodeEdges, neBounds := s.neBounds,
    adjBounds := s.adjBounds, adjDat := s.adjDat }

/-! ## Refinement -/

theorem mergeSort_key {α : Type} (l : List α) (key : α → Nat) :
    ((l.map fun c => (key c, c)).mergeSort fun a b => a.1 ≤ b.1).map (·.2) =
      l.mergeSort fun a b => key a ≤ key b := by
  rw [← List.map_mergeSort (f := fun c => (key c, c)) (r := fun a b => decide (key a ≤ key b))
    (s := fun a b => decide (a.1 ≤ b.1)) (fun _ _ _ _ => rfl), List.map_map]
  exact List.map_id _

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
    refine ⟨?_, ?_⟩ <;> first | rfl | trivial | exact mergeSort_key _ _

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

theorem relabelTreeFast_eq (g : Graph) (items : Array Item) :
    relabelTreeFast g items = relabelTree g items := by
  simp only [relabelTreeFast, relabelTree, relabelRun_toSlow]

end Spqr.Fast
