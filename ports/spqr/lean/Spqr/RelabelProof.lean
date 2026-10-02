import Spqr.RelabelMono
import Spqr.RelabelWp
import Spqr.LayoutSize

/-!
# `relabel`: per-call specification

`relabel_spec`: from a `CallPre` state, a call `relabel fuel cur p pn ct` (with enough fuel) reaches a
`CallPost` state: every item of `desc cur` is numbered, with its `NodeS` record, and nothing below the
entry bounds is touched.
-/

namespace Spqr

namespace Items

variable {items : Items} {g : Graph}

theorem getElem!_type {j : ItemId} (h : j < items.size) : items[j]!.type = items.type j := by
  simp [type, getElem!_pos, h]

theorem getElem!_vs {j : ItemId} (h : j < items.size) : items[j]!.vs = items.vs j := by
  simp [vs, getElem!_pos, h]

end Items

namespace Ghost

open RelabelM

/-- Discharge one arm of a join point, leaving the relation as a `sorry`. -/
macro "jp_armR_sorry" : tactic =>
  `(tactic| first
      | exact ⟨_, by sorry, arm_modify _ _ _⟩
      | exact ⟨_, by sorry, arm_pure _ _⟩
      | exact ⟨_, by sorry, arm_id _ _⟩)

structure StepNum (cur : ItemId) (ty : NodeType) (p pn : Option Nat) (s σ : RelabelState) : Prop where
  g_eq : σ.g = s.g
  items_eq : σ.items = s.items
  vertIndex : σ.vertIndex = s.vertIndex
  edgeIndex : σ.edgeIndex = s.edgeIndex
  edgeFlipped : σ.edgeFlipped = s.edgeFlipped
  chBounds : σ.chBounds = s.chBounds
  chDat : σ.chDat = s.chDat
  nodeVerts : σ.nodeVerts = s.nodeVerts
  nvBounds : σ.nvBounds = s.nvBounds
  nodeEdges : σ.nodeEdges = s.nodeEdges
  neBounds : σ.neBounds = s.neBounds
  adjBounds : σ.adjBounds = s.adjBounds
  adjDat : σ.adjDat = s.adjDat
  vertPos : σ.vertPos = s.vertPos
  order : σ.order = s.order.push cur
  types : σ.types = s.types.push ty
  par : σ.par = s.par.push p
  subtreeEnd : σ.subtreeEnd = s.subtreeEnd.push 0
  vertParNv : σ.vertParNv = s.vertParNv.push pn
  origId : σ.origId = s.origId.push none
structure StepVQ (g : Graph) (items : Items) (cur curIdx : Nat) (s σ : RelabelState) : Prop where
  g_eq : σ.g = s.g
  items_eq : σ.items = s.items
  par : σ.par = s.par
  subtreeEnd : σ.subtreeEnd = s.subtreeEnd
  types : σ.types = s.types
  chBounds : σ.chBounds = s.chBounds
  chDat : σ.chDat = s.chDat
  nodeVerts : σ.nodeVerts = s.nodeVerts
  nvBounds : σ.nvBounds = s.nvBounds
  vertParNv : σ.vertParNv = s.vertParNv
  nodeEdges : σ.nodeEdges = s.nodeEdges
  neBounds : σ.neBounds = s.neBounds
  adjBounds : σ.adjBounds = s.adjBounds
  adjDat : σ.adjDat = s.adjDat
  vertPos : σ.vertPos = s.vertPos
  order : σ.order = s.order
  origId_size : σ.origId.size = s.origId.size
  origId_ne : ∀ k, k ≠ curIdx → σ.origId[k]? = s.origId[k]?
  origId_cur : σ.origId[curIdx]! = items.origOf g cur
  vertIndex_size : σ.vertIndex.size = s.vertIndex.size
  vertIndex : ∀ v, σ.vertIndex[v]! =
    if items.type cur = NodeType.V ∧ v = cur - 1 then some curIdx else s.vertIndex[v]!
  edgeIndex_size : σ.edgeIndex.size = s.edgeIndex.size
  edgeIndex : ∀ e, σ.edgeIndex[e]! =
    if items.type cur = NodeType.Q ∧ e = cur - 1 - g.nv then some curIdx else s.edgeIndex[e]!
  edgeFlipped_size : σ.edgeFlipped.size = s.edgeFlipped.size
  edgeFlipped : ∀ e, σ.edgeFlipped[e]! =
    if items.type cur = NodeType.Q ∧ e = cur - 1 - g.nv then ((items.vs cur).1 != some (g.edges[e]!).1)
    else s.edgeFlipped[e]!
structure StepNV (nv : List NodeVert) (s σ : RelabelState) : Prop where
  g_eq : σ.g = s.g
  items_eq : σ.items = s.items
  vertIndex : σ.vertIndex = s.vertIndex
  edgeIndex : σ.edgeIndex = s.edgeIndex
  edgeFlipped : σ.edgeFlipped = s.edgeFlipped
  par : σ.par = s.par
  subtreeEnd : σ.subtreeEnd = s.subtreeEnd
  types : σ.types = s.types
  origId : σ.origId = s.origId
  chBounds : σ.chBounds = s.chBounds
  chDat : σ.chDat = s.chDat
  nvBounds : σ.nvBounds = s.nvBounds
  vertParNv : σ.vertParNv = s.vertParNv
  nodeEdges : σ.nodeEdges = s.nodeEdges
  neBounds : σ.neBounds = s.neBounds
  adjBounds : σ.adjBounds = s.adjBounds
  adjDat : σ.adjDat = s.adjDat
  vertPos : σ.vertPos = s.vertPos
  order : σ.order = s.order
  nodeVerts : σ.nodeVerts = s.nodeVerts ++ nv.toArray
structure StepPos (ty : NodeType) (nvs : List NodeVert) (nvSt : Nat) (s σ : RelabelState) : Prop where
  g_eq : σ.g = s.g
  items_eq : σ.items = s.items
  vertIndex : σ.vertIndex = s.vertIndex
  edgeIndex : σ.edgeIndex = s.edgeIndex
  edgeFlipped : σ.edgeFlipped = s.edgeFlipped
  par : σ.par = s.par
  subtreeEnd : σ.subtreeEnd = s.subtreeEnd
  types : σ.types = s.types
  origId : σ.origId = s.origId
  chBounds : σ.chBounds = s.chBounds
  chDat : σ.chDat = s.chDat
  nodeVerts : σ.nodeVerts = s.nodeVerts
  nvBounds : σ.nvBounds = s.nvBounds
  vertParNv : σ.vertParNv = s.vertParNv
  nodeEdges : σ.nodeEdges = s.nodeEdges
  neBounds : σ.neBounds = s.neBounds
  adjBounds : σ.adjBounds = s.adjBounds
  adjDat : σ.adjDat = s.adjDat
  order : σ.order = s.order
  vertPos : σ.vertPos = if (ty == .R) = true then
    (nvs.zipIdx nvSt).foldl (fun a x => a.set! x.1.vert x.2) s.vertPos else s.vertPos
structure StepLay (children : List ItemId) (l : Layout) (chEn nvSt nvEn neEn : Nat) (s σ : RelabelState) : Prop where
  g_eq : σ.g = s.g
  items_eq : σ.items = s.items
  vertIndex : σ.vertIndex = s.vertIndex
  edgeIndex : σ.edgeIndex = s.edgeIndex
  edgeFlipped : σ.edgeFlipped = s.edgeFlipped
  par : σ.par = s.par
  subtreeEnd : σ.subtreeEnd = s.subtreeEnd
  types : σ.types = s.types
  origId : σ.origId = s.origId
  chDat : σ.chDat = s.chDat ++ children.toArray
  nodeVerts : σ.nodeVerts = s.nodeVerts
  vertParNv : σ.vertParNv = s.vertParNv
  vertPos : σ.vertPos = s.vertPos
  order : σ.order = s.order
  nodeEdges : σ.nodeEdges = s.nodeEdges ++ l.edges
  adjDat : σ.adjDat = s.adjDat ++ l.adjDat
  adjBounds : σ.adjBounds = s.adjBounds ++ l.adjBounds.extract 1 (2 * (nvEn - nvSt) + 1)
  chBounds : σ.chBounds = s.chBounds.push chEn
  nvBounds : σ.nvBounds = s.nvBounds.push nvEn
  neBounds : σ.neBounds = s.neBounds.push neEn
structure StepCap (b : Bool) (neSt : Nat) (ct : Option Nat) (s σ : RelabelState) : Prop where
  g_eq : σ.g = s.g
  items_eq : σ.items = s.items
  vertIndex : σ.vertIndex = s.vertIndex
  edgeIndex : σ.edgeIndex = s.edgeIndex
  edgeFlipped : σ.edgeFlipped = s.edgeFlipped
  par : σ.par = s.par
  subtreeEnd : σ.subtreeEnd = s.subtreeEnd
  types : σ.types = s.types
  origId : σ.origId = s.origId
  chBounds : σ.chBounds = s.chBounds
  chDat : σ.chDat = s.chDat
  nodeVerts : σ.nodeVerts = s.nodeVerts
  nvBounds : σ.nvBounds = s.nvBounds
  vertParNv : σ.vertParNv = s.vertParNv
  neBounds : σ.neBounds = s.neBounds
  adjBounds : σ.adjBounds = s.adjBounds
  adjDat : σ.adjDat = s.adjDat
  vertPos : σ.vertPos = s.vertPos
  order : σ.order = s.order
  nodeEdges_size : σ.nodeEdges.size = s.nodeEdges.size
  nodeEdges_ne : ∀ k, k ≠ neSt → σ.nodeEdges[k]? = s.nodeEdges[k]?
  nodeEdges_node : ∀ k : Nat, NodeEdge.node σ.nodeEdges[k]! = NodeEdge.node s.nodeEdges[k]!
  nodeEdges_nvs : ∀ k : Nat, NodeEdge.nvs σ.nodeEdges[k]! = NodeEdge.nvs s.nodeEdges[k]!
  cap : b = true → neSt < s.nodeEdges.size → NodeEdge.twin σ.nodeEdges[neSt]! = ct
structure StepEnd (curIdx n : Nat) (s σ : RelabelState) : Prop where
  g_eq : σ.g = s.g
  items_eq : σ.items = s.items
  vertIndex : σ.vertIndex = s.vertIndex
  edgeIndex : σ.edgeIndex = s.edgeIndex
  edgeFlipped : σ.edgeFlipped = s.edgeFlipped
  par : σ.par = s.par
  types : σ.types = s.types
  origId : σ.origId = s.origId
  chBounds : σ.chBounds = s.chBounds
  chDat : σ.chDat = s.chDat
  nodeVerts : σ.nodeVerts = s.nodeVerts
  nvBounds : σ.nvBounds = s.nvBounds
  vertParNv : σ.vertParNv = s.vertParNv
  nodeEdges : σ.nodeEdges = s.nodeEdges
  neBounds : σ.neBounds = s.neBounds
  adjBounds : σ.adjBounds = s.adjBounds
  adjDat : σ.adjDat = s.adjDat
  vertPos : σ.vertPos = s.vertPos
  order : σ.order = s.order
  subtreeEnd : σ.subtreeEnd = s.subtreeEnd.set! curIdx n

/-- Entry condition of a `relabel` call. -/
structure CallPre (g : Graph) (items : Items) (cur : ItemId) (s : RelabelState) : Prop where
  cons : Consistent g items s
  lt : cur < items.size
  fresh : ∀ j ∈ items.desc cur, j ∉ s.order.toList

/-- Exit condition of `relabel fuel cur p pn ct` started from `s`. -/
structure CallPost (g : Graph) (items : Items) (cur : ItemId) (p pn ct : Option Nat) (s s' : RelabelState) : Prop where
  cons : Consistent g items s'
  agree : Agree Bounds.zero s s'
  idx_cur : s'.idx cur = s.types.size
  mem : ∀ j, j ∈ s'.order.toList ↔ j ∈ s.order.toList ∨ j ∈ items.desc cur
  subtree_end : s'.subtreeEnd[s.types.size]! = s'.types.size
  par : s'.par[s.types.size]! = p
  par_nv : s'.vertParNv[s.types.size]! = pn
  ne_st : s'.neBounds[s.types.size]! = s.nodeEdges.size
  cap : items.hasCap cur → NodeEdge.twin s'.nodeEdges[s.nodeEdges.size]! = ct
  nodes : ∀ j ∈ items.desc cur, NodeS g items s' j ∧ LowB items (Bounds.of s) s' j

set_option pp.deepTerms false
set_option pp.deepTerms.threshold 80
set_option pp.proofs.withType false

theorem relabel_spec {g : Graph} {items : Items} (hwf : items.WF g) :
    ∀ (fuel : Nat) (cur : ItemId) (p pn ct : Option Nat) (s : RelabelState),
      CallPre g items cur s → (items.desc cur).card ≤ fuel →
      wp (relabel fuel cur p pn ct) (fun _ s' => CallPost g items cur p pn ct s s') s
  | 0, cur, p, pn, ct, s, hpre, hf => by
    have := Items.one_le_card_desc hpre.lt
    omega
  | fuel + 1, cur, p, pn, ct, s, hpre, hf => by
    have ih := relabel_spec hwf fuel
    have hc := hpre.cons
    have hcur := hpre.lt
    unfold relabel
    simp only [wp_bind, wp_get, wp_item, wp_modify, hc.items_eq, hc.g_eq, Items.getElem!_ch hcur,
      Items.getElem!_type hcur, Items.getElem!_vs hcur]
    try dsimp only
    -- number `cur`
    refine wp_abs _ (StepNum cur (items.type cur) p pn s) _ (by constructor <;> first | rfl | exact hc.g_eq.symm | exact hc.items_eq.symm) fun s1 h1 => ?_
    -- V/Q bookkeeping
    refine wp_jpR _ _ (StepVQ g items cur s.types.size) _ (by split <;> jp_armR_sorry) fun s2 h2 => ?_
    try simp only [wp_bind, wp_get, wp_modify]
    -- node-verts
    refine wp_abs _ (StepNV ((items.nvList g cur).map (⟨s.types.size, ·⟩)) s2) _ (by constructor <;> first | rfl | exact hc.g_eq.symm | exact hc.items_eq.symm) fun s3 h3 => ?_
    try dsimp only
    -- R: vertex positions
    refine wp_jpR _ _ (StepPos (items.type cur) ((items.nvList g cur).map (⟨s.types.size, ·⟩)) s2.nodeVerts.size) _ (by split <;> jp_armR_sorry) fun s4 h4 => ?_
    try simp only [wp_bind, wp_get, wp_modify]
    refine wp_orderedChildren' _ _ fun children hchildren => ?_
    try simp only [wp_bind, wp_get, wp_modify]
    -- child slots, node edges / bounds
    refine wp_abs _ (StepLay children _ _ _ _ _ s4) _ (by constructor <;> first | rfl | exact hc.g_eq.symm | exact hc.items_eq.symm) fun s6 h6 => ?_
    try simp only [wp_bind, wp_get, wp_modify]
    -- cap twin
    split
    · rename_i hcap
      try simp only [wp_bind, wp_get, wp_modify]
      refine wp_abs _ (StepCap true s4.nodeEdges.size ct s6) _ (by sorry) fun s7 h7 => ?_
      trace_state
      sorry
    · rename_i hcap
      refine wp_abs _ (StepCap false s4.nodeEdges.size ct s6) _ (by sorry) fun s7 h7 => ?_
      sorry

end Ghost

end Spqr
