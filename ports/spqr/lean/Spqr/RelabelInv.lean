import Spqr.RelabelSpec
import Spqr.RelabelGhost
import Spqr.ItemTree

/-!
# State-level invariants of the `relabel` fold

The ghost state `Ghost.RelabelState` is viewed as an `SpqrTree` (`tree`) with the preorder index
`idx := order.idxOf`. `Agree B s s'` says `s'` extends `s` without rewriting anything except the
child slots `chDat[≥ B.chDat]` / twins `nodeEdges[≥ B.nodeEdges]` (and the scratch `vertPos`);
`NodeS` is `RelabelNode` plus the raw node-vert entries, and `LowB` locates a finished node's
slots at or above `B`, so that `NodeS` is preserved along `Agree B` (`NodeS.mono`).
-/

namespace Spqr

namespace Items

variable {g : Graph} {items : Items}

theorem type_eq_getElem! (i : ItemId) : items.type i = items[i]!.type := by
  unfold type; by_cases h : i < items.size
  · simp [h]
  · simp [h]; rfl
theorem ch_eq_getElem! (i : ItemId) : items.ch i = items[i]!.ch := by
  unfold ch; by_cases h : i < items.size
  · simp [h]
  · simp [h]; rfl
theorem vs_eq_getElem! (i : ItemId) : items.vs i = items[i]!.vs := by
  unfold vs; by_cases h : i < items.size
  · simp [h]
  · simp [h]; rfl

end Items

/-- `xs[k]!` from `xs[k]?`. -/
theorem Array.getElem!_eq_getD_getElem? {α : Type} [Inhabited α] (xs : Array α) (k : Nat) :
    xs[k]! = xs[k]?.getD default := by
  rw [Array.getElem!_eq_getD, Array.getD_eq_getD_getElem?]

namespace Ghost

/-- The `SpqrTree` read off a ghost state. -/
def RelabelState.tree (g : Graph) (s : RelabelState) : SpqrTree := SpqrTree.ofRelabelState g s.toRelabelState

/-- Preorder index of item `i` in `s` (`s.order.size` if not yet numbered). -/
def RelabelState.idx (s : RelabelState) (i : ItemId) : Nat := s.order.toList.idxOf i

section tree
variable (g : Graph) (s : RelabelState)

@[simp] theorem tree_size : (s.tree g).size = s.types.size := rfl
@[simp] theorem tree_nv : (s.tree g).nv = g.nv := rfl
@[simp] theorem tree_ne : (s.tree g).ne = g.ne := rfl
@[simp] theorem tree_types : (s.tree g).types = s.types := rfl
@[simp] theorem tree_par : (s.tree g).par = s.par := rfl
@[simp] theorem tree_subtreeEnd : (s.tree g).subtreeEnd = s.subtreeEnd := rfl
@[simp] theorem tree_origId : (s.tree g).origId = s.origId := rfl
@[simp] theorem tree_vertIndex : (s.tree g).vertIndex = s.vertIndex := rfl
@[simp] theorem tree_edgeIndex : (s.tree g).edgeIndex = s.edgeIndex := rfl
@[simp] theorem tree_edgeFlipped : (s.tree g).edgeFlipped = s.edgeFlipped := rfl
@[simp] theorem tree_chBounds : (s.tree g).chBounds = s.chBounds := rfl
@[simp] theorem tree_chDat : (s.tree g).chDat = s.chDat := rfl
@[simp] theorem tree_nvBounds : (s.tree g).nvBounds = s.nvBounds := rfl
@[simp] theorem tree_vertParNv : (s.tree g).vertParNv = s.vertParNv := rfl
@[simp] theorem tree_nodeEdges : (s.tree g).nodeEdges = s.nodeEdges := rfl
@[simp] theorem tree_neBounds : (s.tree g).neBounds = s.neBounds := rfl
@[simp] theorem tree_adjBounds : (s.tree g).adjBounds = s.adjBounds := rfl
@[simp] theorem tree_adjDat : (s.tree g).adjDat = s.adjDat := rfl
@[simp] theorem tree_nodeVerts : (s.tree g).nodeVerts =
    s.nodeVerts.map fun nv => { nv with vert := (s.vertIndex[nv.vert]!).getD 0 } := rfl
theorem tree_type (n : Nat) : (s.tree g).type n = s.types[n]?.getD .F := rfl
theorem tree_parent (n : Nat) : (s.tree g).parent n = s.par[n]?.getD none := rfl
theorem tree_chRange (n : Nat) : (s.tree g).chRange n = (s.chBounds[n]?.getD 0, s.chBounds[n + 1]?.getD 0) := rfl
theorem tree_nvRange (n : Nat) : (s.tree g).nvRange n = (s.nvBounds[n]?.getD 0, s.nvBounds[n + 1]?.getD 0) := rfl
theorem tree_neRange (n : Nat) : (s.tree g).neRange n = (s.neBounds[n]?.getD 0, s.neBounds[n + 1]?.getD 0) := rfl
theorem tree_twin (n : Nat) : (s.tree g).twin n = (s.nodeEdges[n]?).bind (·.twin) := rfl

theorem tree_nodeVerts_getElem! (n : Nat) (h : n < s.nodeVerts.size) :
    (s.tree g).nodeVerts[n]! = ⟨s.nodeVerts[n]!.node, (s.vertIndex[s.nodeVerts[n]!.vert]!).getD 0⟩ := by
  rw [tree_nodeVerts, getElem!_pos (s.nodeVerts.map _) n (by simpa using h), Array.getElem_map,
    getElem!_pos s.nodeVerts n h]

end tree

/-- `a'` agrees with `a` from index `b` on (and is at least as long). -/
def PreFrom {α : Type} (b : Nat) (a a' : Array α) : Prop :=
  a.size ≤ a'.size ∧ ∀ k, b ≤ k → k < a.size → a'[k]? = a[k]?

abbrev Pre {α : Type} (a a' : Array α) : Prop := PreFrom 0 a a'

namespace PreFrom

variable {α : Type} {b : Nat} {a a' a'' : Array α}

theorem refl : PreFrom b a a := ⟨le_rfl, fun _ _ _ => rfl⟩
theorem trans (h1 : PreFrom b a a') (h2 : PreFrom b a' a'') : PreFrom b a a'' :=
  ⟨h1.1.trans h2.1, fun k hb hk => (h2.2 k hb (hk.trans_le h1.1)).trans (h1.2 k hb hk)⟩
theorem mono {b' : Nat} (h : PreFrom b a a') (hb : b ≤ b') : PreFrom b' a a' :=
  ⟨h.1, fun k hb' hk => h.2 k (hb.trans hb') hk⟩
theorem size_le (h : PreFrom b a a') : a.size ≤ a'.size := h.1
theorem get? (h : PreFrom b a a') {k : Nat} (hb : b ≤ k) (hk : k < a.size) : a'[k]? = a[k]? := h.2 k hb hk
theorem get! [Inhabited α] (h : PreFrom b a a') {k : Nat} (hb : b ≤ k) (hk : k < a.size) : a'[k]! = a[k]! := by
  rw [Array.getElem!_eq_getD_getElem?, Array.getElem!_eq_getD_getElem?, h.2 k hb hk]
theorem getD (h : PreFrom b a a') {k : Nat} (hb : b ≤ k) (hk : k < a.size) (d : α) :
    a'[k]?.getD d = a[k]?.getD d := by rw [h.2 k hb hk]

theorem push (x : α) : Pre a (a.push x) :=
  ⟨by simp, fun k _ hk => by rw [Array.getElem?_push, ite_eq_right (by omega)]⟩
theorem append (c : Array α) : Pre a (a ++ c) :=
  ⟨by simp, fun k _ hk => Array.getElem?_append_left hk⟩
theorem set! (i : Nat) (x : α) (hi : i < b) : PreFrom b a (a.set! i x) :=
  ⟨by simp, fun k hk _ => by
    rw [Array.set!_eq_setIfInBounds, Array.getElem?_setIfInBounds_ne]; omega⟩
theorem modify (i : Nat) (f : α → α) (hi : i < b) : PreFrom b a (a.modify i f) :=
  ⟨by simp, fun k hk _ => by rw [Array.getElem?_modify, ite_eq_right (by omega)]⟩
theorem set!_ge (i : Nat) (x : α) (hi : a.size ≤ i) : Pre a (a.set! i x) :=
  ⟨by simp, fun k _ hk => by
    rw [Array.set!_eq_setIfInBounds, Array.getElem?_setIfInBounds_ne]; omega⟩
theorem modify_ge (i : Nat) (f : α → α) (hi : a.size ≤ i) : Pre a (a.modify i f) :=
  ⟨by simp, fun k _ hk => by rw [Array.getElem?_modify, ite_eq_right (by omega)]⟩

end PreFrom

/-- Lower bounds on the child slots / node-edges a step may rewrite. -/
structure Bounds where
  chDat : Nat
  nodeEdges : Nat
  /-- Node index: `subtreeEnd` entries at or above it are kept. -/
  idx : Nat

def Bounds.of (s : RelabelState) : Bounds := ⟨s.chDat.size, s.nodeEdges.size, s.types.size⟩

/-- Nothing may be rewritten. -/
def Bounds.zero : Bounds := ⟨0, 0, 0⟩

def Bounds.le (B B' : Bounds) : Prop := B.chDat ≤ B'.chDat ∧ B.nodeEdges ≤ B'.nodeEdges ∧ B.idx ≤ B'.idx

theorem Bounds.zero_le (B : Bounds) : Bounds.zero.le B := ⟨Nat.zero_le _, Nat.zero_le _, Nat.zero_le _⟩
theorem Bounds.le_refl (B : Bounds) : B.le B := ⟨le_rfl, le_rfl, le_rfl⟩

/-- `s'` extends `s`: every array grows at the end, entries are unchanged except `chDat` below
`B.chDat`, the `twin`s of `nodeEdges` below `B.nodeEdges`, `subtreeEnd` below `B.idx`, and
`vertIndex` / `edgeIndex` / `edgeFlipped` entries that were still unset. -/
structure Agree (B : Bounds) (s s' : RelabelState) : Prop where
  g_eq : s'.g = s.g
  items_eq : s'.items = s.items
  order : s.order.toList <+: s'.order.toList
  types : Pre s.types s'.types
  par : Pre s.par s'.par
  subtreeEnd : PreFrom B.idx s.subtreeEnd s'.subtreeEnd
  origId : Pre s.origId s'.origId
  vertParNv : Pre s.vertParNv s'.vertParNv
  chBounds : Pre s.chBounds s'.chBounds
  nvBounds : Pre s.nvBounds s'.nvBounds
  neBounds : Pre s.neBounds s'.neBounds
  nodeVerts : Pre s.nodeVerts s'.nodeVerts
  adjBounds : Pre s.adjBounds s'.adjBounds
  adjDat : Pre s.adjDat s'.adjDat
  chDat : PreFrom B.chDat s.chDat s'.chDat
  nodeEdges : PreFrom B.nodeEdges s.nodeEdges s'.nodeEdges
  nodeEdges_node : ∀ k, k < s.nodeEdges.size → s'.nodeEdges[k]!.node = s.nodeEdges[k]!.node
  nodeEdges_nvs : ∀ k, k < s.nodeEdges.size → s'.nodeEdges[k]!.nvs = s.nodeEdges[k]!.nvs
  vertIndex_size : s'.vertIndex.size = s.vertIndex.size
  vertIndex : ∀ v : Nat, s.vertIndex[v]! ≠ none → s'.vertIndex[v]! = s.vertIndex[v]!
  edgeIndex_size : s'.edgeIndex.size = s.edgeIndex.size
  edgeIndex : ∀ e : Nat, s.edgeIndex[e]! ≠ none → s'.edgeIndex[e]! = s.edgeIndex[e]!
  edgeFlipped_size : s'.edgeFlipped.size = s.edgeFlipped.size
  edgeFlipped : ∀ e : Nat, s.edgeIndex[e]! ≠ none → s'.edgeFlipped[e]! = s.edgeFlipped[e]!
  vertPos_size : s'.vertPos.size = s.vertPos.size

namespace Agree

variable {B B' : Bounds} {s s' s'' : RelabelState}

theorem refl : Agree B s s where
  g_eq := rfl; items_eq := rfl; order := List.prefix_rfl
  types := PreFrom.refl; par := PreFrom.refl; subtreeEnd := PreFrom.refl; origId := PreFrom.refl
  vertParNv := PreFrom.refl; chBounds := PreFrom.refl; nvBounds := PreFrom.refl
  neBounds := PreFrom.refl; nodeVerts := PreFrom.refl; adjBounds := PreFrom.refl
  adjDat := PreFrom.refl; chDat := PreFrom.refl; nodeEdges := PreFrom.refl
  nodeEdges_node := fun _ _ => rfl; nodeEdges_nvs := fun _ _ => rfl
  vertIndex_size := rfl; vertIndex := fun _ _ => rfl; edgeIndex_size := rfl
  edgeIndex := fun _ _ => rfl; edgeFlipped_size := rfl; edgeFlipped := fun _ _ => rfl
  vertPos_size := rfl

theorem trans (h1 : Agree B s s') (h2 : Agree B s' s'') : Agree B s s'' where
  g_eq := h2.g_eq.trans h1.g_eq
  items_eq := h2.items_eq.trans h1.items_eq
  order := h1.order.trans h2.order
  types := h1.types.trans h2.types
  par := h1.par.trans h2.par
  subtreeEnd := h1.subtreeEnd.trans h2.subtreeEnd
  origId := h1.origId.trans h2.origId
  vertParNv := h1.vertParNv.trans h2.vertParNv
  chBounds := h1.chBounds.trans h2.chBounds
  nvBounds := h1.nvBounds.trans h2.nvBounds
  neBounds := h1.neBounds.trans h2.neBounds
  nodeVerts := h1.nodeVerts.trans h2.nodeVerts
  adjBounds := h1.adjBounds.trans h2.adjBounds
  adjDat := h1.adjDat.trans h2.adjDat
  chDat := h1.chDat.trans h2.chDat
  nodeEdges := h1.nodeEdges.trans h2.nodeEdges
  nodeEdges_node := fun k hk =>
    (h2.nodeEdges_node k (hk.trans_le h1.nodeEdges.size_le)).trans (h1.nodeEdges_node k hk)
  nodeEdges_nvs := fun k hk =>
    (h2.nodeEdges_nvs k (hk.trans_le h1.nodeEdges.size_le)).trans (h1.nodeEdges_nvs k hk)
  vertIndex_size := h2.vertIndex_size.trans h1.vertIndex_size
  vertIndex := fun v hv => by
    have e := h1.vertIndex v hv
    rw [h2.vertIndex v (e ▸ hv), e]
  edgeIndex_size := h2.edgeIndex_size.trans h1.edgeIndex_size
  edgeIndex := fun e he => by
    have e' := h1.edgeIndex e he
    rw [h2.edgeIndex e (e' ▸ he), e']
  edgeFlipped_size := h2.edgeFlipped_size.trans h1.edgeFlipped_size
  edgeFlipped := fun e he => by
    have e' := h1.edgeIndex e he
    rw [h2.edgeFlipped e (e' ▸ he), h1.edgeFlipped e he]
  vertPos_size := h2.vertPos_size.trans h1.vertPos_size

theorem mono (h : Agree B s s') (hB : B.le B') : Agree B' s s' :=
  { h with
    chDat := h.chDat.mono hB.1
    nodeEdges := h.nodeEdges.mono hB.2.1
    subtreeEnd := h.subtreeEnd.mono hB.2.2 }

theorem order_size (h : Agree B s s') : s.order.size ≤ s'.order.size := by
  simpa using h.order.length_le

theorem idx_eq (h : Agree B s s') {i : ItemId} (hi : i ∈ s.order.toList) : s'.idx i = s.idx i := by
  obtain ⟨l, hl⟩ := h.order
  unfold RelabelState.idx
  rw [← hl, List.idxOf_append_of_mem hi]

end Agree

/-- Global shape of a ghost state between steps. -/
structure Consistent (g : Graph) (items : Items) (s : RelabelState) : Prop where
  g_eq : s.g = g
  items_eq : s.items = items
  order_size : s.order.size = s.types.size
  order_nodup : s.order.toList.Nodup
  par_size : s.par.size = s.types.size
  subtreeEnd_size : s.subtreeEnd.size = s.types.size
  origId_size : s.origId.size = s.types.size
  vertParNv_size : s.vertParNv.size = s.types.size
  chBounds_size : s.chBounds.size = s.types.size + 1
  nvBounds_size : s.nvBounds.size = s.types.size + 1
  neBounds_size : s.neBounds.size = s.types.size + 1
  ch_last : s.chBounds[s.types.size]! = s.chDat.size
  nv_last : s.nvBounds[s.types.size]! = s.nodeVerts.size
  ne_last : s.neBounds[s.types.size]! = s.nodeEdges.size
  adjBounds_size : s.adjBounds.size = 2 * s.nodeVerts.size + 1
  adjDat_size : s.adjDat.size = 2 * s.nodeEdges.size
  vertIndex_size : s.vertIndex.size = g.nv
  edgeIndex_size : s.edgeIndex.size = g.ne
  edgeFlipped_size : s.edgeFlipped.size = g.ne
  vertPos_size : s.vertPos.size = g.nv
  vert_index_mem : ∀ v : Nat, s.vertIndex[v]! ≠ none → vertItem v ∈ s.order.toList
  edge_index_mem : ∀ e : Nat, s.edgeIndex[e]! ≠ none → edgeItem g e ∈ s.order.toList

theorem Consistent.idx_lt_iff {g : Graph} {items : Items} {s : RelabelState} (hc : Consistent g items s)
    (i : ItemId) : s.idx i < s.types.size ↔ i ∈ s.order.toList := by
  rw [← hc.order_size, RelabelState.idx, ← Array.length_toList, List.idxOf_lt_length_iff]

/-- `RelabelNode` together with the raw node-vert entries (before re-indexing). -/
structure NodeS (g : Graph) (items : Items) (s : RelabelState) (i : ItemId) : Prop where
  node : RelabelNode g items (s.tree g) s.idx i
  raw : ∀ k (hk : k < (items.nvList g i).length),
    s.nodeVerts[((s.tree g).nvRange (s.idx i)).1 + k]! = ⟨s.idx i, (items.nvList g i)[k]⟩

/-- Item `i` and its children are numbered and all the slots `RelabelNode … i` reads lie at or
above `B` and within the arrays. -/
structure LowB (items : Items) (B : Bounds) (s : RelabelState) (i : ItemId) : Prop where
  mem : i ∈ s.order.toList
  idx : B.idx ≤ s.idx i
  chDat : B.chDat ≤ s.chBounds[s.idx i]!
  chDat_le : s.chBounds[s.idx i + 1]! ≤ s.chDat.size
  nodeEdges : B.nodeEdges ≤ s.neBounds[s.idx i]!
  nodeEdges_le : s.neBounds[s.idx i + 1]! ≤ s.nodeEdges.size
  nodeVerts_le : s.nvBounds[s.idx i + 1]! ≤ s.nodeVerts.size
  ch_mem : ∀ c ∈ items.ch i, c ∈ s.order.toList
  ch_idx : ∀ c ∈ items.ch i, B.idx ≤ s.idx c
  ch_ne : ∀ c ∈ items.ch i, B.nodeEdges ≤ s.neBounds[s.idx c]! ∧
    (items.hasCap c → s.neBounds[s.idx c]! < s.nodeEdges.size)

theorem LowB.mono {items : Items} {B B' : Bounds} {s : RelabelState} {i : ItemId} (h : LowB items B s i)
    (hB : B'.le B) : LowB items B' s i :=
  { h with
    idx := hB.2.2.trans h.idx
    chDat := hB.1.trans h.chDat
    nodeEdges := hB.2.1.trans h.nodeEdges
    ch_idx := fun c hc => hB.2.2.trans (h.ch_idx c hc)
    ch_ne := fun c hc => ⟨hB.2.1.trans (h.ch_ne c hc).1, (h.ch_ne c hc).2⟩ }

end Ghost

end Spqr
