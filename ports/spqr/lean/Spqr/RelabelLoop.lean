import Spqr.RelabelMono
import Spqr.RelabelWp
import Spqr.LayoutSize

/-!
# `relabel`: abstract steps and the children-loop invariant

Each primitive of the ghost `relabel` body is described by a `Step*` relation between the states
before and after it; `Entry` collects what the prefix of a call wrote before the children loop,
`LoopInv` is the invariant of that loop, and `CallPre` / `CallPost` are the call's contract. The
`wp` glue lives in `RelabelProof.lean`; everything here is about plain states.
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
section helpers
variable {α : Type} [Inhabited α]
theorem push_get!_last (a : Array α) (x : α) : (a.push x)[a.size]! = x := by
  rw [Array.getElem!_eq_getD_getElem?, Array.getElem?_push_size]; rfl
theorem push_get!_lt (a : Array α) (x : α) {k : Nat} (hk : k < a.size) : (a.push x)[k]! = a[k]! := by
  rw [Array.getElem!_eq_getD_getElem?, Array.getElem!_eq_getD_getElem?, Array.getElem?_push_lt hk, Array.getElem?_eq_getElem hk]
theorem append_get!_right (a b : Array α) (k : Nat) : (a ++ b)[a.size + k]! = b[k]! := by
  rw [Array.getElem!_eq_getD_getElem?, Array.getElem!_eq_getD_getElem?,
    Array.getElem?_append_right (Nat.le_add_right _ _), Nat.add_sub_cancel_left]
theorem append_get!_left (a b : Array α) {k : Nat} (hk : k < a.size) : (a ++ b)[k]! = a[k]! := by
  rw [Array.getElem!_eq_getD_getElem?, Array.getElem!_eq_getD_getElem?, Array.getElem?_append_left hk]
theorem set!_get!_self (a : Array α) {i : Nat} (x : α) (hi : i < a.size) : (a.set! i x)[i]! = x := by
  rw [Array.getElem!_eq_getD_getElem?, Array.set!_eq_setIfInBounds, Array.getElem?_setIfInBounds_self_of_lt hi]
  rfl
theorem set!_get!_ne (a : Array α) {i k : Nat} (x : α) (h : i ≠ k) : (a.set! i x)[k]! = a[k]! := by
  rw [Array.getElem!_eq_getD_getElem?, Array.getElem!_eq_getD_getElem?, Array.set!_eq_setIfInBounds,
    Array.getElem?_setIfInBounds_ne h]
theorem modify_get!_ne (a : Array α) {i k : Nat} (f : α → α) (h : i ≠ k) : (a.modify i f)[k]! = a[k]! := by
  rw [Array.getElem!_eq_getD_getElem?, Array.getElem!_eq_getD_getElem?, Array.getElem?_modify, if_neg h]
theorem modify_get!_self (a : Array α) {i : Nat} (f : α → α) (hi : i < a.size) : (a.modify i f)[i]! = f a[i]! := by
  rw [Array.getElem!_eq_getD_getElem?, Array.getElem!_eq_getD_getElem?, Array.getElem?_modify, if_pos rfl,
    Array.getElem?_eq_getElem hi]; rfl
end helpers

theorem modify_twin_node (a : Array NodeEdge) (j : Nat) (t : Option Nat) (k : Nat) :
    (a.modify j fun ne => { ne with twin := t })[k]!.node = a[k]!.node := by
  rw [Array.getElem!_eq_getD_getElem?, Array.getElem!_eq_getD_getElem?, Array.getElem?_modify]
  split
  · cases a[k]? <;> rfl
  · rfl
theorem modify_twin_nvs (a : Array NodeEdge) (j : Nat) (t : Option Nat) (k : Nat) :
    (a.modify j fun ne => { ne with twin := t })[k]!.nvs = a[k]!.nvs := by
  rw [Array.getElem!_eq_getD_getElem?, Array.getElem!_eq_getD_getElem?, Array.getElem?_modify]
  split
  · cases a[k]? <;> rfl
  · rfl

namespace PreFrom
variable {α : Type} {b : Nat} {a a' : Array α}
theorem set!_of (h : PreFrom b a a') (i : Nat) (x : α) (hi : i < b ∨ a.size ≤ i) : PreFrom b a (a'.set! i x) :=
  ⟨by rw [Array.size_set!]; exact h.1, fun k hk hka => by
    rw [Array.set!_eq_setIfInBounds, Array.getElem?_setIfInBounds_ne (by omega)]; exact h.2 k hk hka⟩
theorem modify_of (h : PreFrom b a a') (i : Nat) (f : α → α) (hi : i < b ∨ a.size ≤ i) :
    PreFrom b a (a'.modify i f) :=
  ⟨by rw [Array.size_modify]; exact h.1, fun k hk hka => by
    rw [Array.getElem?_modify, ite_eq_right (by omega)]; exact h.2 k hk hka⟩
end PreFrom

theorem Items.mem_desc_iff_children {items : Items} {a x : ItemId} (ha : a < items.size) :
    x ∈ items.desc a ↔ x = a ∨ ∃ c ∈ items.ch a, x ∈ items.desc c := by
  rw [Items.mem_desc]
  constructor
  · rintro ⟨hx, hb⟩
    rcases Relation.ReflTransGen.cases_head hb with rfl | ⟨c, hc, hcx⟩
    · exact Or.inl rfl
    · exact Or.inr ⟨c, hc, Items.mem_desc.2 ⟨hx, hcx⟩⟩
  · rintro (rfl | ⟨c, hc, hx⟩)
    · exact ⟨ha, Relation.ReflTransGen.refl⟩
    · rw [Items.mem_desc] at hx; exact ⟨hx.1, Relation.ReflTransGen.head hc hx.2⟩

theorem Items.one_le_nEdges_of_hasCap {g : Graph} {items : Items} {i : ItemId} (h : items.hasCap i = true) : 1 ≤ items.nEdges g i := by
  unfold Items.nEdges Items.capCount; rw [if_pos h]; omega

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
  origId : σ.origId =
    if items.type cur = .V then s.origId.set! curIdx (some (cur - 1))
    else if items.type cur = .Q then s.origId.set! curIdx (some (cur - 1 - g.nv)) else s.origId
  vertIndex : σ.vertIndex =
    if items.type cur = .V then s.vertIndex.set! (cur - 1) (some curIdx) else s.vertIndex
  edgeIndex : σ.edgeIndex =
    if items.type cur = .Q then s.edgeIndex.set! (cur - 1 - g.nv) (some curIdx) else s.edgeIndex
  edgeFlipped : σ.edgeFlipped =
    if items.type cur = .Q then
      s.edgeFlipped.set! (cur - 1 - g.nv) ((items.vs cur).1 != some (g.edges[cur - 1 - g.nv]!).1)
    else s.edgeFlipped
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
  nodeEdges : σ.nodeEdges =
    if b = true then s.nodeEdges.modify neSt (fun ne => { ne with twin := ct }) else s.nodeEdges
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

/-- Child slot `i` set to `n`; with `tw = some j`, node-edge `j` is twinned with the next node-edge
to be written. -/
structure StepSlotT (i n : Nat) (tw : Option Nat) (s σ : RelabelState) : Prop where
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
  nodeVerts : σ.nodeVerts = s.nodeVerts
  nvBounds : σ.nvBounds = s.nvBounds
  vertParNv : σ.vertParNv = s.vertParNv
  neBounds : σ.neBounds = s.neBounds
  adjBounds : σ.adjBounds = s.adjBounds
  adjDat : σ.adjDat = s.adjDat
  vertPos : σ.vertPos = s.vertPos
  order : σ.order = s.order
  chDat : σ.chDat = s.chDat.set! i n
  nodeEdges : σ.nodeEdges =
    tw.elim s.nodeEdges fun j => s.nodeEdges.modify j fun ne => { ne with twin := some s.nodeEdges.size }

/-- Entry condition of a `relabel` call. -/
structure CallPre (g : Graph) (items : Items) (cur : ItemId) (s : RelabelState) : Prop where
  cons : Consistent g items s
  lt : cur < items.size
  fresh : ∀ j ∈ items.desc cur, j ∉ s.order.toList

/-- Exit condition of `relabel fuel cur p pn ct` started from `s`. -/
structure CallPost (g : Graph) (items : Items) (cur : ItemId) (p pn ct : Option Nat) (s s' : RelabelState) :
    Prop where
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

/-- The skeleton the call for `cur` appends, in terms of the entry positions. -/
def entryLayout (g : Graph) (items : Items) (cur : ItemId) (curIdx nvSt neSt : Nat) (pos : Nat → Nat)
    (children : List ItemId) : Layout :=
  layoutNode (items.type cur) curIdx nvSt (nvSt + (items.nvList g cur).length) neSt
    (neSt + items.nEdges g cur) (items.edgeChildren g pos children)

/-- What the prefix of the call for `cur` (numbering, V/Q bookkeeping, node-verts, `vertPos`,
skeleton, cap twin) wrote, as seen at the loop-entry state `s7`; `s` is the call's entry state. -/
structure Entry (g : Graph) (items : Items) (cur : ItemId) (p pn ct : Option Nat) (curIdx chSt nvSt neSt : Nat)
    (pos : Nat → Nat) (children : List ItemId) (s s7 : RelabelState) : Prop where
  cons : Consistent g items s7
  agree : Agree Bounds.zero s s7
  curIdx_eq : curIdx = s.types.size
  chSt_eq : chSt = s.chDat.size
  nvSt_eq : nvSt = s.nodeVerts.size
  neSt_eq : neSt = s.nodeEdges.size
  cur_not_mem : cur ∉ s.order.toList
  order : s7.order = s.order.push cur
  types : s7.types = s.types.push (items.type cur)
  par : s7.par = s.par.push p
  vertParNv : s7.vertParNv = s.vertParNv.push pn
  subtreeEnd : s7.subtreeEnd = s.subtreeEnd.push 0
  origId : s7.origId[curIdx]! = items.origOf g cur
  vertIndex : items.type cur = .V → s7.vertIndex[cur - 1]! = some curIdx
  edgeIndex : items.type cur = .Q → s7.edgeIndex[cur - 1 - g.nv]! = some curIdx ∧
    s7.edgeFlipped[cur - 1 - g.nv]! = ((items.vs cur).1 != some (g.edges[cur - 1 - g.nv]!).1)
  children_eq : children = items.ordered g cur nvSt pos
  pos_ok : items.type cur = .R → Items.PosOK nvSt (items.nvList g cur) pos
  chBounds : s7.chBounds = s.chBounds.push (chSt + children.length)
  nvBounds : s7.nvBounds = s.nvBounds.push (nvSt + (items.nvList g cur).length)
  neBounds : s7.neBounds = s.neBounds.push (neSt + items.nEdges g cur)
  chDat : s7.chDat = s.chDat ++ children.toArray
  nodeVerts : s7.nodeVerts = s.nodeVerts ++ ((items.nvList g cur).map (⟨curIdx, ·⟩)).toArray
  nodeEdges_size : s7.nodeEdges.size = neSt + items.nEdges g cur
  edge_node : ∀ k, k < items.nEdges g cur →
    (s7.nodeEdges[neSt + k]!).node = ((entryLayout g items cur curIdx nvSt neSt pos children).edges[k]!).node
  edge_nvs : ∀ k, k < items.nEdges g cur →
    (s7.nodeEdges[neSt + k]!).nvs = ((entryLayout g items cur curIdx nvSt neSt pos children).edges[k]!).nvs
  adj_bounds : ∀ j, 1 ≤ j → j ≤ 2 * (items.nvList g cur).length →
    s7.adjBounds[2 * nvSt + j]! = (entryLayout g items cur curIdx nvSt neSt pos children).adjBounds[j]!
  adj_dat : ∀ j, j < 2 * items.nEdges g cur →
    s7.adjDat[2 * neSt + j]! = (entryLayout g items cur curIdx nvSt neSt pos children).adjDat[j]!
  cap : items.hasCap cur → (s7.nodeEdges[neSt]!).twin = ct

/-- Invariant of the children loop of the call for `cur`: `done ++ rest = children`, the loop
counters `b = (nvCur, neCur, chCur)`, and what the calls for `done` wrote. -/
structure LoopInv (g : Graph) (items : Items) (cur : ItemId) (ct : Option Nat) (curIdx chSt nvSt neSt : Nat)
    (s s7 : RelabelState) (children done rest : List ItemId) (b : Nat × Nat × Nat) (σ : RelabelState) : Prop where
  split : children = done ++ rest
  cons : Consistent g items σ
  agree0 : Agree Bounds.zero s σ
  agree7 : Agree (Bounds.of s7) s7 σ
  idx_cur : σ.idx cur = curIdx
  ch_cur : b.2.2 = done.length
  nv_cur : b.1 = nvSt + (items.vs cur).1.toList.length + done.countP (· < 1 + g.nv)
  ne_cur : b.2.1 = neSt + items.capCount cur +
    (if (items.type cur).isNode then done.countP (· ≥ 1 + g.nv) else 0)
  mem : ∀ j, j ∈ σ.order.toList ↔ j ∈ s7.order.toList ∨ ∃ c ∈ done, j ∈ items.desc c
  size : σ.types.size = match done.getLast? with
    | none => curIdx + 1
    | some c => σ.subtreeEnd[σ.idx c]!
  slots : ∀ k (hk : k < done.length), σ.chDat[chSt + k]! = σ.idx done[k]
  child_idx : ∀ k (hk : k < done.length),
    σ.idx done[k] = if k = 0 then curIdx + 1 else σ.subtreeEnd[σ.idx done[k - 1]]!
  par : ∀ c ∈ done, σ.par[σ.idx c]! = some curIdx
  par_nv : ∀ k (hk : k < (done.filter (· < 1 + g.nv)).length),
    σ.vertParNv[σ.idx (done.filter (· < 1 + g.nv))[k]]! = some (nvSt + (items.vs cur).1.toList.length + k)
  par_nv_none : ∀ c ∈ done, 1 + g.nv ≤ c → σ.vertParNv[σ.idx c]! = none
  twin : (items.type cur).isNode → ∀ k (hk : k < (done.filter (· ≥ 1 + g.nv)).length),
    (σ.tree g).twin (neSt + items.capCount cur + k) =
      some σ.neBounds[σ.idx (done.filter (· ≥ 1 + g.nv))[k]]! ∧
    (items.hasCap (done.filter (· ≥ 1 + g.nv))[k] →
      (σ.tree g).twin σ.neBounds[σ.idx (done.filter (· ≥ 1 + g.nv))[k]]! =
        some (neSt + items.capCount cur + k))
  cap_none : ¬ (items.type cur).isNode → ∀ c ∈ done, items.hasCap c →
    (σ.tree g).twin σ.neBounds[σ.idx c]! = none
  cap : items.hasCap cur → (σ.nodeEdges[neSt]!).twin = ct
  nodes : ∀ c ∈ done, ∀ j ∈ items.desc c, NodeS g items σ j ∧ LowB items (Bounds.of s7) σ j

variable {g : Graph} {items : Items} {cur : ItemId} {p pn ct : Option Nat} {curIdx chSt nvSt neSt : Nat}
  {pos : Nat → Nat} {children done rest : List ItemId} {s s7 σ σ1 σ3 : RelabelState} {b : Nat × Nat × Nat}
  {c : ItemId}

theorem Entry.isEmpty_children (he : Entry g items cur p pn ct curIdx chSt nvSt neSt pos children s s7) :
    children.isEmpty = (items.ch cur).isEmpty := by
  rw [he.children_eq]
  have := (Items.ordered_perm (items := items) (g := g) (i := cur) (nvSt := nvSt) (pos := pos)).length_eq
  exact Bool.eq_iff_iff.2 (by simp only [List.isEmpty_iff_length_eq_zero, this])

theorem Entry.hasCap_eq (he : Entry g items cur p pn ct curIdx chSt nvSt neSt pos children s s7) :
    items.hasCap cur = ((items.type cur).isNode && !(items.type cur == .Q && !children.isEmpty)) := by
  rw [Items.hasCap, he.isEmpty_children]

theorem Entry.nodup_children (hwf : items.WF g)
    (he : Entry g items cur p pn ct curIdx chSt nvSt neSt pos children s s7) : children.Nodup := by
  rw [he.children_eq]; exact Items.ordered_perm.nodup_iff.2 (hwf.tree.ch_nodup cur)

theorem entry_of_steps {s1 s2 s3 s4 s6 : RelabelState} {b : Bool} (hwf : items.WF g) (hpre : CallPre g items cur s)
    (ct : Option Nat)
    (h1 : StepNum cur (items.type cur) p pn s s1) (h2 : StepVQ g items cur s.types.size s1 s2)
    (h3 : StepNV ((items.nvList g cur).map (fun v => (⟨s.types.size, v⟩ : NodeVert))) s2 s3)
    (h4 : StepPos (items.type cur) ((items.nvList g cur).map (fun v => (⟨s.types.size, v⟩ : NodeVert))) s2.nodeVerts.size s3 s4)
    (hch : children = RelabelM.orderedList items[cur]! s2.nodeVerts.size s4)
    (h6 : StepLay children
      (layoutNode (items.type cur) s.types.size s2.nodeVerts.size
        (s2.nodeVerts.size + ((items.nvList g cur).map (fun v => (⟨s.types.size, v⟩ : NodeVert))).length)
        s4.nodeEdges.size
        (s4.nodeEdges.size + ((if (items.type cur).isNode = true then children.countP (· ≥ 1 + g.nv) else 0) +
          if ((items.type cur).isNode && !(items.type cur == NodeType.Q && !children.isEmpty)) = true then 1 else 0))
        ((children.filter (· ≥ 1 + g.nv)).map fun c =>
          (s4.vertPos[s4.items[c]!.vs.1.getD 0]!, s4.vertPos[s4.items[c]!.vs.2.getD 0]!)))
      (s2.chDat.size + children.length) s2.nodeVerts.size
      (s2.nodeVerts.size + ((items.nvList g cur).map (fun v => (⟨s.types.size, v⟩ : NodeVert))).length)
      (s4.nodeEdges.size + ((if (items.type cur).isNode = true then children.countP (· ≥ 1 + g.nv) else 0) +
        if ((items.type cur).isNode && !(items.type cur == NodeType.Q && !children.isEmpty)) = true then 1 else 0))
      s4 s6)
    (h7 : StepCap b s4.nodeEdges.size ct s6 s7)
    (hb : b = true ↔ ((items.type cur).isNode && !(items.type cur == NodeType.Q && !children.isEmpty)) = true) :
    Entry g items cur p pn ct s.types.size s2.chDat.size s2.nodeVerts.size s4.nodeEdges.size
      (fun v => s4.vertPos[v]!) children s s7 := by
  sorry

theorem loop_init (hpre : CallPre g items cur s)
    (he : Entry g items cur p pn ct curIdx chSt nvSt neSt pos children s s7) {nv0 ne0 : Nat}
    (hnv : nv0 = nvSt + if (items.vs cur).1.isSome = true then 1 else 0)
    (hne : ne0 = neSt + if items.hasCap cur then 1 else 0) :
    LoopInv g items cur ct curIdx chSt nvSt neSt s s7 children [] children (nv0, ne0, 0) s7 where
  split := rfl
  cons := he.cons
  agree0 := he.agree
  agree7 := Agree.refl
  idx_cur := by
    rw [RelabelState.idx, he.order, Array.toList_push, List.idxOf_append, if_neg he.cur_not_mem,
      List.idxOf_cons_self, Array.length_toList, hpre.cons.order_size, he.curIdx_eq, Nat.zero_add]
  ch_cur := rfl
  nv_cur := by rw [hnv, List.countP_nil, Nat.add_zero]; cases (items.vs cur).1 <;> rfl
  ne_cur := by simp only [hne, Items.capCount, List.countP_nil, ite_self, Nat.add_zero]
  mem := fun j => by simp
  size := by show s7.types.size = curIdx + 1; rw [he.types, Array.size_push, he.curIdx_eq]
  slots := fun k hk => absurd hk (Nat.not_lt_zero _)
  child_idx := fun k hk => absurd hk (Nat.not_lt_zero _)
  par := fun c hc => nomatch hc
  par_nv := fun k hk => absurd hk (Nat.not_lt_zero _)
  par_nv_none := fun c hc => nomatch hc
  twin := fun _ k hk => absurd hk (Nat.not_lt_zero _)
  cap_none := fun _ c hc => nomatch hc
  cap := he.cap
  nodes := fun c hc => nomatch hc

theorem loop_fuel {fuel : Nat} (ht : items.Tree g) (hcur : cur < items.size) (hc : c ∈ items.ch cur)
    (hf : (items.desc cur).card ≤ fuel + 1) : (items.desc c).card ≤ fuel := by
  have h := Finset.card_le_card (ht.desc_subset hc)
  rw [Finset.card_erase_of_mem (Items.mem_desc_self hcur)] at h
  have := Items.one_le_card_desc hcur
  omega

theorem LoopInv.mem_ch (he : Entry g items cur p pn ct curIdx chSt nvSt neSt pos children s s7)
    (hi : LoopInv g items cur ct curIdx chSt nvSt neSt s s7 children done (c :: rest) b σ) : c ∈ items.ch cur := by
  rw [← Items.mem_ordered_iff (g := g) (nvSt := nvSt) (pos := pos), ← he.children_eq, hi.split]
  simp

/-- Transfer `Consistent` along equations of all fields. -/
macro "cons_transfer" hc:term : tactic =>
  `(tactic| first
    | exact ($hc).g_eq | exact ($hc).items_eq | exact ($hc).order_size | exact ($hc).order_nodup | exact ($hc).par_size
    | exact ($hc).subtreeEnd_size | exact ($hc).origId_size | exact ($hc).vertParNv_size | exact ($hc).chBounds_size
    | exact ($hc).nvBounds_size | exact ($hc).neBounds_size | exact ($hc).ch_last | exact ($hc).nv_last | exact ($hc).ne_last
    | exact ($hc).adjBounds_size | exact ($hc).adjDat_size | exact ($hc).vertIndex_size | exact ($hc).edgeIndex_size
    | exact ($hc).edgeFlipped_size | exact ($hc).vertPos_size | exact ($hc).vert_index_mem | exact ($hc).edge_index_mem)

theorem LoopInv.mem_ch_done (he : Entry g items cur p pn ct curIdx chSt nvSt neSt pos children s s7)
    (hi : LoopInv g items cur ct curIdx chSt nvSt neSt s s7 children done rest b σ) (hc : c ∈ done) :
    c ∈ items.ch cur := by
  rw [← Items.mem_ordered_iff (g := g) (nvSt := nvSt) (pos := pos), ← he.children_eq, hi.split]
  exact List.mem_append_left _ hc

theorem StepSlotT.nodeEdges_size {i n : Nat} {tw : Option Nat} (h : StepSlotT i n tw σ σ1) :
    σ1.nodeEdges.size = σ.nodeEdges.size := by
  rw [h.nodeEdges]; cases tw <;> simp

theorem StepSlotT.nodeEdges_node {i n : Nat} {tw : Option Nat} (h : StepSlotT i n tw σ σ1) (k : Nat) :
    σ1.nodeEdges[k]!.node = σ.nodeEdges[k]!.node := by
  rw [h.nodeEdges]; cases tw
  · rfl
  · exact modify_twin_node _ _ _ _

theorem StepSlotT.nodeEdges_nvs {i n : Nat} {tw : Option Nat} (h : StepSlotT i n tw σ σ1) (k : Nat) :
    σ1.nodeEdges[k]!.nvs = σ.nodeEdges[k]!.nvs := by
  rw [h.nodeEdges]; cases tw
  · rfl
  · exact modify_twin_nvs _ _ _ _

theorem StepSlotT.cons {i n : Nat} {tw : Option Nat} (hc : Consistent g items σ) (h : StepSlotT i n tw σ σ1) :
    Consistent g items σ1 := by
  have hsz := h.nodeEdges_size
  obtain ⟨hg, hit, hvi, hei, hef, hpar, hse, hty, hor, hcb, hnv, hnvb, hvp, hneb, hab, had, hvpos, hord,
    hcd, -⟩ := h
  constructor <;> simp only [hg, hit, hvi, hei, hef, hpar, hse, hty, hor, hcb, hnv, hnvb, hvp, hneb, hab, had,
    hvpos, hord, hcd, hsz, Array.size_set!] <;> cons_transfer hc

theorem StepEnd.cons {curIdx n : Nat} (hc : Consistent g items σ) (h : StepEnd curIdx n σ σ1) :
    Consistent g items σ1 := by
  obtain ⟨hg, hit, hvi, hei, hef, hpar, hty, hor, hcb, hcd, hnv, hnvb, hvp, hne, hneb, hab, had, hvpos, hord,
    hse⟩ := h
  constructor <;> simp only [hg, hit, hvi, hei, hef, hpar, hse, hty, hor, hcb, hnv, hnvb, hvp, hne, hneb, hab, had,
    hvpos, hord, hcd, Array.size_set!] <;> cons_transfer hc

/-- Transfer `Agree B s0 _` along equations of all fields. -/
macro "agree_transfer" ha:term : tactic =>
  `(tactic| first
    | exact ($ha).g_eq | exact ($ha).items_eq | exact ($ha).order | exact ($ha).types | exact ($ha).par
    | exact ($ha).subtreeEnd | exact ($ha).origId | exact ($ha).vertParNv | exact ($ha).chBounds
    | exact ($ha).nvBounds | exact ($ha).neBounds | exact ($ha).nodeVerts | exact ($ha).adjBounds
    | exact ($ha).adjDat | exact ($ha).chDat | exact ($ha).nodeEdges | exact ($ha).nodeEdges_node
    | exact ($ha).nodeEdges_nvs | exact ($ha).vertIndex_size | exact ($ha).vertIndex | exact ($ha).edgeIndex_size
    | exact ($ha).edgeIndex | exact ($ha).edgeFlipped_size | exact ($ha).edgeFlipped | exact ($ha).vertPos_size)

theorem StepSlotT.agree {B : Bounds} {s0 : RelabelState} {i n : Nat} {tw : Option Nat}
    (h : StepSlotT i n tw σ σ1) (ha : Agree B s0 σ) (hi : i < B.chDat ∨ s0.chDat.size ≤ i)
    (hj : ∀ j, tw = some j → j < B.nodeEdges ∨ s0.nodeEdges.size ≤ j) : Agree B s0 σ1 := by
  have hnode := h.nodeEdges_node
  have hnvs := h.nodeEdges_nvs
  obtain ⟨hg, hit, hvi, hei, hef, hpar, hse, hty, hor, hcb, hnv, hnvb, hvp, hneb, hab, had, hvpos, hord,
    hcd, hne⟩ := h
  constructor <;> (try simp only [hg, hit, hvi, hei, hef, hpar, hse, hty, hor, hcb, hnv, hnvb, hvp, hneb, hab, had,
    hvpos, hord, hcd]) <;> first
    | agree_transfer ha
    | exact PreFrom.set!_of ha.chDat _ _ hi
    | (rw [hne]; cases tw with
        | none => exact ha.nodeEdges
        | some j => exact PreFrom.modify_of ha.nodeEdges _ _ (hj j rfl))
    | (intro k hk; rw [hnode]; exact ha.nodeEdges_node k hk)
    | (intro k hk; rw [hnvs]; exact ha.nodeEdges_nvs k hk)

theorem StepEnd.agree {B : Bounds} {s0 : RelabelState} {curIdx n : Nat} (h : StepEnd curIdx n σ σ1)
    (ha : Agree B s0 σ) (hi : curIdx < B.idx ∨ s0.subtreeEnd.size ≤ curIdx) : Agree B s0 σ1 := by
  obtain ⟨hg, hit, hvi, hei, hef, hpar, hty, hor, hcb, hcd, hnv, hnvb, hvp, hne, hneb, hab, had, hvpos, hord,
    hse⟩ := h
  constructor <;> simp only [hg, hit, hvi, hei, hef, hpar, hse, hty, hor, hcb, hnv, hnvb, hvp, hne, hneb, hab, had,
    hvpos, hord, hcd] <;> first
    | agree_transfer ha
    | exact PreFrom.set!_of ha.subtreeEnd _ _ hi

/-- The call for the next child `c` starts from a `CallPre` state. -/
theorem LoopInv.pre (hwf : items.WF g) (hpre : CallPre g items cur s)
    (he : Entry g items cur p pn ct curIdx chSt nvSt neSt pos children s s7)
    (hi : LoopInv g items cur ct curIdx chSt nvSt neSt s s7 children done (c :: rest) b σ)
    (hcons : Consistent g items σ1) (horder : σ1.order = σ.order) : CallPre g items c σ1 := by
  have hcc := hi.mem_ch he
  refine ⟨hcons, hwf.tree.ch_lt cur c hcc, fun j hj hmem => ?_⟩
  have hjc : j ∈ (items.desc cur).erase cur := hwf.tree.desc_subset hcc hj
  rw [Finset.mem_erase] at hjc
  rw [horder, hi.mem, he.order, Array.toList_push, List.mem_append, List.mem_singleton] at hmem
  rcases hmem with (hmem | rfl) | ⟨c', hc', hjc'⟩
  · exact hpre.fresh j hjc.2 hmem
  · exact hjc.1 rfl
  · have hnd := he.nodup_children hwf
    rw [hi.split, List.nodup_append] at hnd
    have hne : c' ≠ c := hnd.2.2 c' hc' c List.mem_cons_self
    exact Finset.disjoint_left.1 (hwf.tree.desc_disjoint (hi.mem_ch_done he hc') hcc hne) hjc' hj

/-- One iteration of the children loop: the slot write (and twin for a non-V child of a node),
then the recursive call for `c`. -/
theorem LoopInv.step (hwf : items.WF g) (hpre : CallPre g items cur s)
    (he : Entry g items cur p pn ct curIdx chSt nvSt neSt pos children s s7)
    (hi : LoopInv g items cur ct curIdx chSt nvSt neSt s s7 children done (c :: rest) b σ)
    {tw pn' ct' : Option Nat} {b' : Nat × Nat × Nat}
    (hcase : (c < 1 + g.nv ∧ tw = none ∧ pn' = some b.1 ∧ ct' = none ∧ b' = (b.1 + 1, b.2.1, b.2.2 + 1)) ∨
      (¬ c < 1 + g.nv ∧ (items.type cur).isNode = true ∧ tw = some b.2.1 ∧ pn' = none ∧ ct' = some b.2.1 ∧
        b' = (b.1, b.2.1 + 1, b.2.2 + 1)) ∨
      (¬ c < 1 + g.nv ∧ ¬ (items.type cur).isNode = true ∧ tw = none ∧ pn' = none ∧ ct' = none ∧
        b' = (b.1, b.2.1, b.2.2 + 1)))
    (h1 : StepSlotT (chSt + b.2.2) σ.types.size tw σ σ1)
    (h3 : CallPost g items c (some curIdx) pn' ct' σ1 σ3) :
    LoopInv g items cur ct curIdx chSt nvSt neSt s s7 children (done ++ [c]) rest b' σ3 := by
  sorry

theorem NodeS.neSt_lt {B : Bounds} {i : ItemId} (hn : NodeS g items σ i) (hl : LowB items B σ i)
    (hcap : items.hasCap i = true) : σ.neBounds[σ.idx i]! < σ.nodeEdges.size := by
  have h1 := hn.node.ne_range
  rw [tree_neRange_fst, tree_neRange_snd] at h1
  have h2 := hl.nodeEdges_le
  have h3 := Items.one_le_nEdges_of_hasCap (g := g) hcap
  omega

/-- After the loop, writing `subtreeEnd[curIdx]` finishes the call. -/
theorem LoopInv.fin (hwf : items.WF g) (hpre : CallPre g items cur s)
    (he : Entry g items cur p pn ct curIdx chSt nvSt neSt pos children s s7)
    (hi : LoopInv g items cur ct curIdx chSt nvSt neSt s s7 children done [] b σ)
    (hend : StepEnd curIdx σ.types.size σ σ1) : CallPost g items cur p pn ct s σ1 := by
  have hdone : children = done := by rw [hi.split, List.append_nil]
  subst hdone
  have hcs := hi.cons
  have hc1 : Consistent g items σ1 := hend.cons hcs
  have hce := he.curIdx_eq
  have hcse := he.chSt_eq
  have hnse := he.nvSt_eq
  have hnese := he.neSt_eq
  have hs7size : s7.types.size = curIdx + 1 := by rw [he.types, Array.size_push, hce]
  have hs7ch : s7.chDat.size = chSt + children.length := by
    rw [he.chDat]; simp [hcse]
  have hs7nv : s7.nodeVerts.size = nvSt + (items.nvList g cur).length := by
    rw [he.nodeVerts]; simp [hnse]
  have hs7ne := he.nodeEdges_size
  have hlen : children.length = (items.ch cur).length := by
    rw [he.children_eq]
    exact (Items.ordered_perm (items := items) (g := g) (i := cur) (nvSt := nvSt) (pos := pos)).length_eq
  have hmemch : ∀ c ∈ items.ch cur, c ∈ children := fun c hc => by
    rw [he.children_eq]; exact Items.mem_ordered_iff.2 hc
  have hchmem : ∀ c ∈ children, c ∈ items.ch cur := fun c hc =>
    Items.mem_ordered_iff.1 (he.children_eq ▸ hc)
  have hcur7 : cur ∈ s7.order.toList := by rw [he.order]; simp
  have hcurσ : cur ∈ σ.order.toList := (hi.mem cur).2 (Or.inl hcur7)
  have hcσ : ∀ c ∈ items.ch cur, c ∈ σ.order.toList := fun c hc =>
    (hi.mem c).2 (Or.inr ⟨c, hmemch c hc, Items.mem_desc_self (hwf.tree.ch_lt cur c hc)⟩)
  have hord1 : σ1.order = σ.order := hend.order
  have hidx1 : ∀ j, σ1.idx j = σ.idx j := fun j => by unfold RelabelState.idx; rw [hord1]
  have hic : σ1.idx cur = curIdx := by rw [hidx1, hi.idx_cur]
  have hcurlt : curIdx < σ.types.size := by
    have := hcs.idx_lt hcurσ; rwa [hi.idx_cur] at this
  have hlow : ∀ c ∈ children, curIdx + 1 ≤ σ.idx c := fun c hc => by
    have h := (hi.nodes c hc c (Items.mem_desc_self (hwf.tree.ch_lt cur c (hchmem c hc)))).2.idx
    have h' : s7.types.size ≤ σ.idx c := h
    omega
  have hendA : Agree (Bounds.of s7) σ σ1 :=
    hend.agree Agree.refl (Or.inl (by show curIdx < s7.types.size; omega))
  have hagree0 : Agree Bounds.zero s σ1 :=
    hend.agree hi.agree0 (Or.inr (by rw [hpre.cons.subtreeEnd_size, hce]))
  have htree_twin : ∀ n, (σ1.tree g).twin n = (σ.tree g).twin n := fun n => by
    rw [tree_twin, tree_twin, hend.nodeEdges]
  have hchB0 : σ.chBounds[curIdx]! = chSt := by
    rw [hi.agree7.chBounds_get! he.cons (by omega), he.chBounds, hce,
      push_get!_lt _ _ (by rw [hpre.cons.chBounds_size]; omega), hpre.cons.ch_last, hcse]
  have hchB1 : σ.chBounds[curIdx + 1]! = chSt + children.length := by
    rw [hi.agree7.chBounds_get! he.cons (by omega), he.chBounds, hce, ← hpre.cons.chBounds_size,
      push_get!_last]
  have hnvB0 : σ.nvBounds[curIdx]! = nvSt := by
    rw [hi.agree7.nvBounds_get! he.cons (by omega), he.nvBounds, hce,
      push_get!_lt _ _ (by rw [hpre.cons.nvBounds_size]; omega), hpre.cons.nv_last, hnse]
  have hnvB1 : σ.nvBounds[curIdx + 1]! = nvSt + (items.nvList g cur).length := by
    rw [hi.agree7.nvBounds_get! he.cons (by omega), he.nvBounds, hce, ← hpre.cons.nvBounds_size,
      push_get!_last]
  have hneB0 : σ.neBounds[curIdx]! = neSt := by
    rw [hi.agree7.neBounds_get! he.cons (by omega), he.neBounds, hce,
      push_get!_lt _ _ (by rw [hpre.cons.neBounds_size]; omega), hpre.cons.ne_last, hnese]
  have hneB1 : σ.neBounds[curIdx + 1]! = neSt + items.nEdges g cur := by
    rw [hi.agree7.neBounds_get! he.cons (by omega), he.neBounds, hce, ← hpre.cons.neBounds_size,
      push_get!_last]
  have hchSt : ((σ1.tree g).chRange curIdx).1 = chSt := by rw [tree_chRange_fst, hend.chBounds, hchB0]
  have hchEn : ((σ1.tree g).chRange curIdx).2 = chSt + children.length := by
    rw [tree_chRange_snd, hend.chBounds, hchB1]
  have hnvSt : ((σ1.tree g).nvRange curIdx).1 = nvSt := by rw [tree_nvRange_fst, hend.nvBounds, hnvB0]
  have hnvEn : ((σ1.tree g).nvRange curIdx).2 = nvSt + (items.nvList g cur).length := by
    rw [tree_nvRange_snd, hend.nvBounds, hnvB1]
  have hneSt : ((σ1.tree g).neRange curIdx).1 = neSt := by rw [tree_neRange_fst, hend.neBounds, hneB0]
  have hneEn : ((σ1.tree g).neRange curIdx).2 = neSt + items.nEdges g cur := by
    rw [tree_neRange_snd, hend.neBounds, hneB1]
  have hty : (σ1.tree g).type curIdx = items.type cur := by
    rw [tree_type, hend.types, hi.agree7.types_getD (by omega), he.types, hce, Array.getElem?_push_size]; rfl
  have hlay : nodeLayout g items (σ1.tree g) σ1.idx cur pos =
      entryLayout g items cur curIdx nvSt neSt pos children := by
    unfold nodeLayout entryLayout
    rw [hic, hnvSt, hnvEn, hneSt, hneEn, he.children_eq]
  have hnv1 : s7.nodeVerts.size ≤ σ1.nodeVerts.size := by
    rw [hend.nodeVerts]; exact hi.agree7.nodeVerts.size_le
  have hraw : ∀ k (hk : k < (items.nvList g cur).length),
      σ1.nodeVerts[nvSt + k]! = ⟨curIdx, (items.nvList g cur)[k]⟩ := by
    intro k hk
    rw [hend.nodeVerts, hi.agree7.nodeVerts_get! (by omega), he.nodeVerts, hnse, append_get!_right,
      getElem!_pos _ _ (by simpa using hk)]
    simp
  have hnodeS : NodeS g items σ1 cur := by
    refine ⟨{
      idx_lt := ?_
      type := ?_
      orig := ?_
      vert_index := ?_
      edge_index := ?_
      ch_range := ?_
      nv_range := ?_
      ne_range := ?_
      node_verts := ?_
      child_lt := ?_
      child_par := ?_
      child_par_nv_none := ?_
      child_cap_twin_none := ?_
      layout := ⟨pos, ?_⟩ }, ?_⟩
    · rw [hic, tree_size, hend.types]; exact hcurlt
    · rw [hic]; exact hty
    · rw [hic, tree_origId, hend.origId, hi.agree7.origId_get! he.cons (by omega), he.origId]
    · intro hV
      rw [hic, tree_vertIndex, hend.vertIndex,
        hi.agree7.vertIndex _ (by rw [he.vertIndex hV]; exact Option.some_ne_none _), he.vertIndex hV]
    · intro hQ
      obtain ⟨h1, h2⟩ := he.edgeIndex hQ
      have hne : s7.edgeIndex[cur - 1 - g.nv]! ≠ none := by rw [h1]; exact Option.some_ne_none _
      refine ⟨?_, ?_⟩
      · rw [hic, tree_edgeIndex, hend.edgeIndex, hi.agree7.edgeIndex _ hne, h1]
      · rw [tree_edgeFlipped, hend.edgeFlipped, hi.agree7.edgeFlipped _ hne, h2]
    · rw [hic, hchEn, hchSt, hlen]
    · rw [hic, hnvEn, hnvSt]
    · rw [hic, hneEn, hneSt]
    · intro k hk
      rw [hic, hnvSt, tree_nodeVerts_getElem! g σ1 _ (by omega), hraw k hk]; rfl
    · intro c hc
      rw [hidx1, tree_size, hend.types]; exact hcs.idx_lt (hcσ c hc)
    · intro c hc
      rw [tree_parent_eq, hidx1, hic, hend.par]; exact hi.par c (hmemch c hc)
    · intro c hc hge
      rw [tree_vertParNv, hidx1, hend.vertParNv]; exact hi.par_nv_none c (hmemch c hc) hge
    · intro hn c hc hcap
      rw [htree_twin, tree_neRange_fst, hidx1, hend.neBounds]; exact hi.cap_none hn c (hmemch c hc) hcap
    · refine {
        pos_ok := ?_
        children := ?_
        edge_node := ?_
        edge_nvs := ?_
        adj_bounds := ?_
        adj_dat := ?_
        vert_par_nv := ?_
        twin := ?_
        child_idx := ?_
        subtree_end := ?_ }
      · intro hR; rw [hic, hnvSt]; exact he.pos_ok hR
      · rw [hic, hnvSt, ← he.children_eq]
        unfold SpqrTree.children
        rw [hchEn, hchSt, Nat.add_sub_cancel_left]
        refine List.ext_getElem (by simp) fun k h1 h2 => ?_
        rw [List.getElem_map, List.getElem_map, List.getElem_range, tree_chDat, hend.chDat]
        have := hi.slots k (by simpa using h2)
        rw [Array.getElem!_eq_getD_getElem?] at this
        exact this.trans (hidx1 _).symm
      · intro k hk
        rw [hic, hneSt, hlay, tree_nodeEdges, hend.nodeEdges, hi.agree7.nodeEdges_node _ (by omega)]
        exact he.edge_node k hk
      · intro k hk
        rw [hic, hneSt, hlay, tree_nodeEdges, hend.nodeEdges, hi.agree7.nodeEdges_nvs _ (by omega)]
        exact he.edge_nvs k hk
      · intro j hj1 hj2
        rw [hic, hnvSt, hlay, tree_adjBounds, hend.adjBounds, hi.agree7.adjBounds_get! he.cons (by omega)]
        exact he.adj_bounds j hj1 hj2
      · intro j hj
        rw [hic, hneSt, hlay, tree_adjDat, hend.adjDat, hi.agree7.adjDat_get! he.cons (by omega)]
        exact he.adj_dat j hj
      · rw [hic, hnvSt, ← he.children_eq]
        intro k hk
        rw [tree_vertParNv, hidx1, hend.vertParNv]; exact hi.par_nv k hk
      · intro hn
        rw [hic, hneSt, hnvSt, ← he.children_eq]
        intro k hk
        simp only [htree_twin, tree_neRange_fst, hidx1, hend.neBounds]
        exact hi.twin hn k hk
      · simp only [hic, hnvSt, ← he.children_eq]
        intro k hk
        simp only [hidx1, tree_subtreeEnd]
        rw [hend.subtreeEnd]
        by_cases h0 : k = 0
        · rw [if_pos h0]; have := hi.child_idx k hk; rwa [if_pos h0] at this
        · rw [if_neg h0, set!_get!_ne _ _ (Nat.ne_of_lt (hlow _ (List.getElem_mem _)))]
          have := hi.child_idx k hk; rwa [if_neg h0] at this
      · rw [hic, hnvSt, ← he.children_eq, tree_subtreeEnd, hend.subtreeEnd,
          set!_get!_self _ _ (by rw [hcs.subtreeEnd_size]; exact hcurlt)]
        have hsz := hi.size
        cases hl : children.getLast? with
        | none => show σ.types.size = curIdx + 1; rw [hsz, hl]
        | some c =>
          show σ.types.size = (σ.subtreeEnd.set! curIdx σ.types.size)[σ1.idx c]!
          rw [hidx1, set!_get!_ne _ _ (Nat.ne_of_lt (hlow c (List.mem_of_getLast? hl))), hsz, hl]
    · intro k hk
      rw [hic, hnvSt, hend.nodeVerts]
      have := hraw k hk; rw [hend.nodeVerts] at this; exact this
  have hBle : (Bounds.of s).le (Bounds.of s7) := by
    refine ⟨?_, ?_, ?_⟩ <;> simp only [Bounds.of] <;> omega
  have hlowB : LowB items (Bounds.of s) σ1 cur := {
    mem := by rw [hord1]; exact hcurσ
    idx := by show s.types.size ≤ σ1.idx cur; rw [hic]; omega
    chDat := by show s.chDat.size ≤ σ1.chBounds[σ1.idx cur]!; rw [hic, hend.chBounds, hchB0]; omega
    chDat_le := by
      rw [hic, hend.chBounds, hchB1, hend.chDat]
      have := hi.agree7.chDat.size_le; omega
    nodeEdges := by show s.nodeEdges.size ≤ σ1.neBounds[σ1.idx cur]!; rw [hic, hend.neBounds, hneB0]; omega
    nodeEdges_le := by
      rw [hic, hend.neBounds, hneB1, hend.nodeEdges]
      have := hi.agree7.nodeEdges.size_le; omega
    nodeVerts_le := by
      rw [hic, hend.nvBounds, hnvB1, hend.nodeVerts]
      have := hi.agree7.nodeVerts.size_le; omega
    ch_mem := fun c hc => by rw [hord1]; exact hcσ c hc
    ch_idx := fun c hc => by
      show s.types.size ≤ σ1.idx c; rw [hidx1]; have := hlow c (hmemch c hc); omega
    ch_ne := fun c hc => by
      rw [hidx1, hend.neBounds, hend.nodeEdges]
      have hn := hi.nodes c (hmemch c hc) c (Items.mem_desc_self (hwf.tree.ch_lt cur c hc))
      refine ⟨?_, fun hcap => hn.1.neSt_lt hn.2 hcap⟩
      have h7 : s7.nodeEdges.size ≤ σ.neBounds[σ.idx c]! := hn.2.nodeEdges
      show s.nodeEdges.size ≤ _; omega }
  exact {
    cons := hc1
    agree := hagree0
    idx_cur := by rw [hic, hce]
    mem := fun j => by
      rw [hord1, hi.mem j, he.order, Array.toList_push, List.mem_append, List.mem_singleton,
        Items.mem_desc_iff_children hpre.lt]
      constructor
      · rintro ((h | h) | ⟨c, hc, hj⟩)
        · exact Or.inl h
        · exact Or.inr (Or.inl h)
        · exact Or.inr (Or.inr ⟨c, hchmem c hc, hj⟩)
      · rintro (h | h | ⟨c, hc, hj⟩)
        · exact Or.inl (Or.inl h)
        · exact Or.inl (Or.inr h)
        · exact Or.inr ⟨c, hmemch c hc, hj⟩
    subtree_end := by
      rw [← hce, hend.subtreeEnd, hend.types, set!_get!_self _ _ (by rw [hcs.subtreeEnd_size]; exact hcurlt)]
    par := by
      rw [← hce, hend.par, hi.agree7.par_get! he.cons (by omega), he.par, hce, ← hpre.cons.par_size,
        push_get!_last]
    par_nv := by
      rw [← hce, hend.vertParNv, hi.agree7.vertParNv_get! he.cons (by omega), he.vertParNv, hce,
        ← hpre.cons.vertParNv_size, push_get!_last]
    ne_st := by rw [← hce, hend.neBounds, hneB0, hnese]
    cap := fun h => by rw [← hnese, hend.nodeEdges]; exact hi.cap h
    nodes := fun j hj => by
      rcases (Items.mem_desc_iff_children hpre.lt).1 hj with rfl | ⟨c, hc, hj'⟩
      · exact ⟨hnodeS, hlowB⟩
      · obtain ⟨hn, hl⟩ := hi.nodes c (hmemch c hc) j hj'
        exact ⟨hn.mono hcs hl hendA, (hl.agree hcs hendA).mono hBle⟩ }

end Ghost

end Spqr
