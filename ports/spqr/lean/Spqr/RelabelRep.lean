import Spqr.RelabelSpec
import Spqr.ItemTree
import Spqr.Proofs.Dfs
import Spqr.LayoutShape

/-!
# Relabel phase D: `Represents` transport from the per-node interface

Everything here is proved for an abstract output tree `t` satisfying the phase-3 interface
`RelabelIdx` + `RelabelNode` (`Spqr/RelabelSpec.lean`), packaged as `RelabelOK`; the only
admission reaching `relabelTree` is `relabel_node_spec`.
-/

namespace Spqr

/-- `omega` after unfolding the `ItemId` abbreviation and the item constructors. -/
macro "riomega" : tactic =>
  `(tactic| ((try simp only [ItemId, vertItem, edgeItem, rootItem] at *); omega))

/-- `a[i]! = a[i]?.getD 0` for `Nat` arrays (the two spellings used by `Relabel.lean`). -/
theorem Array.getElem!_nat (a : Array Nat) (i : Nat) : a[i]! = a[i]?.getD 0 := by
  simp only [getElem!_def]; cases a[i]? <;> rfl

theorem Array.getElem!_opt {α : Type} (a : Array (Option α)) (i : Nat) :
    a[i]! = a[i]?.getD none := by
  simp only [getElem!_def]; cases a[i]? <;> rfl

/-- Item-level facts the representation transport needs beyond `Items.WF`. All hold for the
walk's output (`PROOF.md` §5); `Items.WF` alone does not imply them. -/
structure Items.RepOK (g : Graph) (items : Items) : Prop where
  /-- A block-root Q with children `[c]` is a self-loop; with `[c, v]` its edge is `{u, v}`. -/
  q_root : ∀ e, e < g.ne →
    (∀ c, items.ch (edgeItem g e) = [c] → (g.edges[e]!).1 = (g.edges[e]!).2) ∧
    (∀ c w u, items.ch (edgeItem g e) = [c, vertItem w] → (items.vs (edgeItem g e)).1 = some u →
      Items.PairEq (u, w) g.edges[e]!)
  /-- O items hang alone under a (self-loop) Q. -/
  o_parent : ∀ p c, items.IsParent p c → items.type c = .O → items.type p = .Q ∧ items.ch p = [c]
  /-- S children are in path order (positional `s_shape`). -/
  s_order : ∀ i, i < items.size → items.type i = .S → ∃ u v xs,
    items.vs i = (some u, some v) ∧ ((items.ch i).filter fun c => items.type c = .V) = xs.map vertItem ∧
    items.virtualEdges i = List.zip (u :: xs) (xs ++ [v])

/-! ### `layoutNode` facts, in `Layout`-local form (derived from `Spqr.LayoutShape`) -/

namespace LayoutFacts

open LayoutShape LayoutR

theorem one_vert (ty : NodeType) (ht : ty.isNode = true) (node nvSt neSt neEn : Nat)
    (ec : List (Nat × Nat)) (hne : neEn - neSt = 1) :
    ((layoutNode ty node nvSt (nvSt + 1) neSt neEn ec).edges[0]!).nvs = (nvSt, nvSt) := by
  rw [layoutNode_loop_eq ty _ _ _ _ _ ec (by cases ty <;> simp [NodeType.isNode] at ht ⊢)
    (by omega), runLoop_edges_get _ _ _ _ _ hne]

theorem qi_edge (ty : NodeType) (ht : ty = .Q ∨ ty = .I) (node nvSt neSt neEn : Nat)
    (ec : List (Nat × Nat)) (hne : neEn - neSt = 1) :
    ((layoutNode ty node nvSt (nvSt + 2) neSt neEn ec).edges[0]!).nvs = (nvSt, nvSt + 1) := by
  rw [layoutNode_QI_eq ty _ _ _ _ _ ec ht (by omega), runQI_edges_get _ _ _ _ _ hne]

theorem p_edge (node nvSt neSt neEn : Nat) (ec : List (Nat × Nat)) (k : Nat) (hk : k < neEn - neSt) :
    ((layoutNode .P node nvSt (nvSt + 2) neSt neEn ec).edges[k]!).nvs = (nvSt, nvSt + 1) := by
  rw [layoutNode_P_eq _ _ _ _ _ _ (by omega)]
  obtain ⟨-, -, -, ge, -⟩ := runP_spec node nvSt (nvSt + 2) neSt neEn (by omega)
  rw [ge k hk]

theorem s_edge (node nvSt nvEn neSt neEn : Nat) (ec : List (Nat × Nat)) (hv : nvSt + 2 ≤ nvEn)
    (hne : neEn - neSt = nvEn - nvSt) (k : Nat) (hk : k < neEn - neSt) :
    ((layoutNode .S node nvSt nvEn neSt neEn ec).edges[k]!).nvs =
      if k = 0 then (nvSt, nvEn - 1) else (nvSt + k - 1, nvSt + k) := by
  rw [layoutNode_S_eq _ _ _ _ _ _ (by omega)]
  obtain ⟨-, -, -, ge, -⟩ := runS_spec node nvSt nvEn neSt neEn (by omega)
  rw [ge k hk]; split_ifs <;> rfl

theorem foldl_countStep_edges (nvSt : Nat) (L : List (Nat × Nat)) (l : Layout) :
    (L.foldl (countStep nvSt) l).edges = l.edges := by
  induction L generalizing l with
  | nil => rfl
  | cons p L ih => rw [List.foldl_cons, ih, countStep_edges]

/-- `LayoutShape.run_edges` without the orientation hypothesis (the edge slots do not depend on
the adjacency counts). -/
theorem r_edge (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat)) (hv : nvSt + 2 ≤ nvEn)
    (hne : neEn = neSt + E.length + 1) (k : Nat) (hk : k < E.length + 1) :
    ((layoutNode .R node nvSt nvEn neSt neEn E).edges[k]!).nvs =
      if k = 0 then (nvSt, nvEn - 1) else E[k - 1]! := by
  rw [layoutNode_R_eq _ _ _ _ _ _ hv]
  simp only [run]
  set l₀ := inc nvSt (inc nvSt (Layout.empty (nvEn - nvSt) (neEn - neSt)) (2 * nvSt + 2))
    (2 * nvEn - 1) with hl₀
  have hsz0 : l₀.adjBounds.size = 2 * (nvEn - nvSt) + 1 := by simp [hl₀, inc, Layout.empty]
  have he0 : l₀.edges = (Layout.empty (nvEn - nvSt) (neEn - neSt)).edges := rfl
  set l₁ := E.foldl (countStep nvSt) l₀ with hl₁
  have e1 : l₁.edges = l₀.edges := foldl_countStep_edges nvSt E l₀
  have sz1 : l₁.adjBounds.size = l₀.adjBounds.size := (foldl_countStep_sizes nvSt E l₀).1
  obtain ⟨e2, -, -, -, -⟩ := foldl_prefixStep nvSt (2 * nvEn + 1 - (2 * nvSt + 1)) (2 * nvSt + 1)
    l₁ (2 * neSt) (Nat.le_refl _) (by rw [sz1, hsz0]; omega)
  set l₂ := ((List.range' (2 * nvSt + 1) (2 * nvEn + 1 - (2 * nvSt + 1))).foldl (prefixStep nvSt)
    (l₁, 2 * neSt)).1 with hl₂
  have he2 : (inc nvSt l₂ (2 * nvSt + 2)).edges = (Layout.empty (nvEn - nvSt) (neEn - neSt)).edges := by
    show l₂.edges = _
    rw [e2, e1, he0]
  have hesz : (inc nvSt l₂ (2 * nvSt + 2)).edges.size = E.length + 1 := by
    rw [he2, empty_edges_size]; omega
  obtain ⟨sz3, g3⟩ := foldl_fillStep_edges nvSt neSt node E.reverse (inc nvSt l₂ (2 * nvSt + 2)) neEn
    (by simp; omega) (by rw [hesz]; omega)
  set l₃ := (E.reverse.foldl (fillStep nvSt neSt node) (inc nvSt l₂ (2 * nvSt + 2), neEn)).1 with hl₃
  show (l₃.setNe neSt node neSt (nvSt, nvEn - 1) (2 * neSt, 2 * neEn - 1)).edges[k]!.nvs = _
  rw [setNe_edges_get _ _ _ _ _ _ _ (by rw [sz3, hesz]; exact hk), Nat.sub_self]
  by_cases hk0 : k = 0
  · rw [ite_of_pos hk0, ite_of_pos hk0]
  · rw [ite_of_neg hk0, ite_of_neg hk0, g3 k (by rw [hesz]; exact hk)]
    simp only [List.length_reverse]
    rw [ite_of_pos (by omega), List.getD_eq_getElem?_getD, List.getElem?_reverse (by omega),
      show E.length - 1 - (neEn - 1 - neSt - k) = k - 1 by omega]
    simp only [getElem!_def]
    cases E[k - 1]? <;> rfl

end LayoutFacts

/-- Item-level skeleton of `i`: its cap followed by its virtual edges, as pairs of positions in
`items.nvList g i`. -/
def Items.rSkeleton (g : Graph) (items : Items) (i : ItemId) : List (Nat × Nat) :=
  (((items.vs i).1.getD 0, (items.vs i).2.getD 0) :: items.virtualEdges i).map fun q =>
    ((items.nvList g i).idxOf q.1, (items.nvList g i).idxOf q.2)

/-- Item-level R 3-connectivity: the contract the R-content proof must establish; here it is only
transported to the output tree (`RelabelOK.r_three_connected`), not proved. -/
def Items.RThreeConnected (g : Graph) (items : Items) : Prop :=
  ∀ i, i < items.size → items.type i = .R →
    SpqrTree.ThreeConnected (items.nvList g i).length (items.rSkeleton g i)

theorem SpqrTree.ThreeConnected_congr {n : Nat} {es es' : List (Nat × Nat)}
    (h : ∀ p, p ∈ es ↔ p ∈ es') :
    SpqrTree.ThreeConnected n es ↔ SpqrTree.ThreeConnected n es' := by
  simp only [SpqrTree.ThreeConnected, h]

theorem Array.getElem!_eq_getD {α : Type} [Inhabited α] (a : Array α) (i : Nat) :
    a[i]! = a[i]?.getD default := by
  simp only [getElem!_def]; cases a[i]? <;> rfl

/-- The phase-3 interface as a hypothesis on an abstract output `t`. -/
structure RelabelOK (g : Graph) (items : Items) (t : SpqrTree) (idx : ItemId → Nat) : Prop where
  wf : items.WF g
  ridx : RelabelIdx g items t idx
  node : ∀ i, i < items.size → RelabelNode g items t idx i

namespace RelabelOK

variable {g : Graph} {items : Items} {t : SpqrTree} {idx : ItemId → Nat}
  (h : RelabelOK g items t idx)
include h

theorem tree : items.Tree g := h.wf.tree
theorem endpoints : items.Endpoints g := h.wf.endpoints
theorem shapes : items.Shapes g := h.wf.shapes

theorem size : t.size = items.size := h.ridx.size

theorem idx_lt {i : ItemId} (hi : i < items.size) : idx i < t.size := (h.node i hi).idx_lt

theorem idx_inj {i j : ItemId} (hi : i < items.size) (hj : j < items.size)
    (hij : idx i = idx j) : i = j := h.ridx.inj i j hi hj hij

theorem idx_surj {n : Nat} (hn : n < t.size) : ∃ i, i < items.size ∧ idx i = n := by
  have hsz := h.size
  have := Finset.surj_on_of_inj_on_of_card_le (s := Finset.range items.size)
    (t := Finset.range t.size) (fun i _ => idx i)
    (fun i hi => by
      simp only [Finset.mem_range] at hi ⊢; exact h.idx_lt hi)
    (fun i j hi hj hij => by
      simp only [Finset.mem_range] at hi hj; exact h.idx_inj hi hj hij)
    (by simp [hsz])
  obtain ⟨i, hi, rfl⟩ := this n (by simpa using hn)
  exact ⟨i, by simpa using hi, rfl⟩

theorem type_eq {i : ItemId} (hi : i < items.size) : t.type (idx i) = items.type i :=
  (h.node i hi).type

/-! ### Item-level facts from `Items.Tree` -/

theorem ch_lt {i c : ItemId} (_hi : i < items.size) (hc : c ∈ items.ch i) : c < items.size :=
  h.tree.ch_lt i c hc

theorem ch_ne_root {i c : ItemId} (hc : c ∈ items.ch i) : c ≠ rootItem := by
  intro hc0; subst hc0; exact h.tree.root_no_parent i hc

theorem vertItem_lt {v : Nat} (hv : v < g.nv) : vertItem v < items.size := by
  have := h.tree.size; riomega

theorem edgeItem_lt {e : Nat} (he : e < g.ne) : edgeItem g e < items.size := by
  have := h.tree.size; riomega

theorem type_vertItem {v : Nat} (hv : v < g.nv) : items.type (vertItem v) = .V :=
  h.tree.vert v hv

theorem type_edgeItem {e : Nat} (he : e < g.ne) : items.type (edgeItem g e) = .Q :=
  h.tree.edge e he

/-- A child below `1 + g.nv` is a vertex item. -/
theorem ch_lt_vert {i c : ItemId} (hc : c ∈ items.ch i) (hlt : c < 1 + g.nv) :
    c = vertItem (c - 1) ∧ items.type c = .V := by
  have h0 : c ≠ 0 := h.ch_ne_root hc
  have : c = vertItem (c - 1) := by riomega
  refine ⟨this, ?_⟩
  rw [this]; exact h.type_vertItem (by riomega)

theorem type_V_iff {i c : ItemId} (hi : i < items.size) (hc : c ∈ items.ch i) :
    items.type c = .V ↔ c < 1 + g.nv := by
  constructor
  · intro hV
    by_contra hge
    have hcs := h.ch_lt hi hc
    by_cases hq : c < 1 + g.nv + g.ne
    · have : c = edgeItem g (c - (1 + g.nv)) := by riomega
      rw [this] at hV
      rw [h.type_edgeItem (by riomega)] at hV; cases hV
    · have := h.tree.node c (by riomega) hcs
      simp [hV] at this
  · intro hlt; exact (h.ch_lt_vert hc hlt).2

omit h in
theorem ch_of_ge {q : ItemId} (hq : items.size ≤ q) : items.ch q = [] := by
  simp [Items.ch, Array.getElem?_eq_none_iff.2 hq]

omit h in
theorem parent_lt {p c : ItemId} (hc : c ∈ items.ch p) : p < items.size := by
  by_contra hge
  rw [ch_of_ge (by riomega)] at hc
  exact List.not_mem_nil hc

theorem card_desc_child {i c : ItemId} (hi : i < items.size) (hc : c ∈ items.ch i) :
    (items.desc c).card < (items.desc i).card :=
  lt_of_le_of_lt (Finset.card_le_card (h.tree.desc_subset hc))
    (Finset.card_lt_card (Finset.erase_ssubset (Items.mem_desc_self hi)))

/-! ### Ranges and bounds -/

omit h in
theorem nvRange_fst (n : Nat) : (t.nvRange n).1 = t.nvBounds[n]! := by
  simp [SpqrTree.nvRange, Array.getElem!_nat]
omit h in
theorem nvRange_snd (n : Nat) : (t.nvRange n).2 = t.nvBounds[n + 1]! := by
  simp [SpqrTree.nvRange, Array.getElem!_nat]
omit h in
theorem neRange_fst (n : Nat) : (t.neRange n).1 = t.neBounds[n]! := by
  simp [SpqrTree.neRange, Array.getElem!_nat]
omit h in
theorem neRange_snd (n : Nat) : (t.neRange n).2 = t.neBounds[n + 1]! := by
  simp [SpqrTree.neRange, Array.getElem!_nat]
omit h in
theorem inSubtree_iff' (i j : Nat) : t.inSubtree i j ↔ i ≤ j ∧ j < t.subtreeEnd[i]! := by
  simp [SpqrTree.inSubtree, Array.getElem!_nat]

theorem nvBounds_step {n : Nat} (hn : n < t.size) : t.nvBounds[n]! ≤ t.nvBounds[n + 1]! := by
  obtain ⟨i, hi, rfl⟩ := h.idx_surj hn
  have := (h.node i hi).nv_range
  rw [nvRange_fst, nvRange_snd] at this; omega

theorem neBounds_step {n : Nat} (hn : n < t.size) : t.neBounds[n]! ≤ t.neBounds[n + 1]! := by
  obtain ⟨i, hi, rfl⟩ := h.idx_surj hn
  have := (h.node i hi).ne_range
  rw [neRange_fst, neRange_snd] at this; omega

theorem nvBounds_mono {n m : Nat} (hnm : n ≤ m) (hm : m ≤ t.size) :
    t.nvBounds[n]! ≤ t.nvBounds[m]! := by
  induction m with
  | zero => obtain rfl : n = 0 := by omega
            exact le_rfl
  | succ k ih =>
    rcases Nat.eq_or_lt_of_le hnm with rfl | hlt
    · exact le_rfl
    · exact (ih (by omega) (by omega)).trans (h.nvBounds_step (by omega))

theorem neBounds_mono {n m : Nat} (hnm : n ≤ m) (hm : m ≤ t.size) :
    t.neBounds[n]! ≤ t.neBounds[m]! := by
  induction m with
  | zero => obtain rfl : n = 0 := by omega
            exact le_rfl
  | succ k ih =>
    rcases Nat.eq_or_lt_of_le hnm with rfl | hlt
    · exact le_rfl
    · exact (ih (by omega) (by omega)).trans (h.neBounds_step (by omega))

theorem nvEn_le {i : ItemId} (hi : i < items.size) : (t.nvRange (idx i)).2 ≤ t.nodeVerts.size := by
  rw [nvRange_snd, ← h.ridx.nv_last]; exact h.nvBounds_mono (h.idx_lt hi) le_rfl

theorem neEn_le {i : ItemId} (hi : i < items.size) : (t.neRange (idx i)).2 ≤ t.nodeEdges.size := by
  rw [neRange_snd, ← h.ridx.ne_last]; exact h.neBounds_mono (h.idx_lt hi) le_rfl

/-! ### Node-verts and their original vertices -/

theorem nvList_lt {i : ItemId} (_hi : i < items.size) {v : Nat} (hv : v ∈ items.nvList g i) :
    v < g.nv := by
  unfold Items.nvList at hv
  simp only [List.mem_append, List.mem_map, List.mem_filter, Option.mem_toList,
    decide_eq_true_eq] at hv
  rcases hv with (hv | ⟨c, ⟨hc, hlt⟩, rfl⟩) | hv
  · exact h.endpoints.vs_lt i v (Or.inl hv)
  · have h0 : c ≠ 0 := h.ch_ne_root hc
    riomega
  · exact h.endpoints.vs_lt i v (Or.inr hv)

theorem nodeVerts_get {i : ItemId} (hi : i < items.size) {k : Nat}
    (hk : k < (items.nvList g i).length) :
    (t.nvRange (idx i)).1 + k < t.nodeVerts.size ∧
    t.nodeVerts[(t.nvRange (idx i)).1 + k]? =
      some ⟨idx i, (t.vertIndex[(items.nvList g i)[k]]!).getD 0⟩ := by
  have hlt : (t.nvRange (idx i)).1 + k < t.nodeVerts.size := by
    have := (h.node i hi).nv_range; have := h.nvEn_le hi; omega
  refine ⟨hlt, ?_⟩
  rw [Array.getElem?_eq_getElem hlt, ← getElem!_pos]
  exact congrArg some ((h.node i hi).node_verts k hk)

theorem origId_vert {v : Nat} (hv : v < g.nv) : t.origId[(t.vertIndex[v]!).getD 0]! = some v := by
  rw [h.ridx.vert_index v hv, Option.getD_some]
  have horig := (h.node _ (h.vertItem_lt hv)).orig
  rw [Items.origOf, h.type_vertItem hv] at horig
  rw [horig]
  show some (1 + v - 1) = some v
  congr 1; omega

theorem nvOrig_get {i : ItemId} (hi : i < items.size) {k : Nat}
    (hk : k < (items.nvList g i).length) :
    t.nvOrig ((t.nvRange (idx i)).1 + k) = some (items.nvList g i)[k] := by
  obtain ⟨-, hnv⟩ := h.nodeVerts_get hi hk
  have hv := h.nvList_lt hi (List.getElem_mem hk)
  simp only [SpqrTree.nvOrig, hnv, ← Array.getElem!_opt]
  exact h.origId_vert hv

theorem nodeVertsOf_length {i : ItemId} (hi : i < items.size) :
    (t.nodeVertsOf (idx i)).length = (items.nvList g i).length := by
  have := (h.node i hi).nv_range
  simp only [SpqrTree.nodeVertsOf, List.length_map, List.length_range]; omega

theorem nodeVertsOf_get {i : ItemId} (hi : i < items.size) {k : Nat}
    (hk : k < (items.nvList g i).length) :
    (t.nodeVertsOf (idx i))[k]? = some ⟨idx i, (t.vertIndex[(items.nvList g i)[k]]!).getD 0⟩ := by
  have hr := (h.node i hi).nv_range
  obtain ⟨-, hnv⟩ := h.nodeVerts_get hi hk
  simp only [SpqrTree.nodeVertsOf, List.getElem?_map, List.getElem?_range, hr, Nat.add_sub_cancel_left,
    hk, Option.map_some, hnv, Option.getD_some]

theorem nv_orig_map {i : ItemId} (hi : i < items.size) :
    ((t.nodeVertsOf (idx i)).map fun nv => t.origId[nv.vert]!) = (items.nvList g i).map some := by
  apply List.ext_getElem
  · simp [h.nodeVertsOf_length hi]
  · intro k h1 h2
    simp only [List.getElem_map]
    rw [(List.getElem_eq_iff _).2 (h.nodeVertsOf_get hi (by simpa using h2))]
    exact h.origId_vert (h.nvList_lt hi (List.getElem_mem _))

/-! ### Preorder intervals -/

omit h in
/-- Consecutive intervals `[lo c, hi c)` over `L` starting at `a` cover exactly `[a, b)`, `b` the
last `hi`. -/
theorem interval_cover {α : Type} (L : List α) (lo hi : α → Nat) (a : Nat)
    (hlo : ∀ k (hk : k < L.length), lo L[k] = if k = 0 then a else hi L[k - 1])
    (hlt : ∀ c ∈ L, lo c < hi c) (b : Nat)
    (hb : b = match L.getLast? with | none => a | some c => hi c) :
    a ≤ b ∧ ∀ x, (a ≤ x ∧ x < b) ↔ ∃ k, ∃ hk : k < L.length, lo L[k] ≤ x ∧ x < hi L[k] := by
  subst hb
  rcases L with _ | ⟨c0, L'⟩
  · simp
  set L := c0 :: L' with hL
  have hlen : 0 < L.length := by simp [hL]
  have hlast : L.getLast? = some L[L.length - 1] := by
    rw [List.getLast?_eq_getElem?, List.getElem?_eq_getElem (by omega)]
  rw [hlast]
  simp only
  have hA : ∀ k (hk : k < L.length), a ≤ lo L[k] := by
    intro k
    induction k with
    | zero => intro hk; rw [hlo 0 hk]; simp
    | succ k ih =>
      intro hk
      rw [hlo (k + 1) hk]; simp only [Nat.add_one_ne_zero, ite_false, Nat.add_sub_cancel]
      exact (ih (by omega)).trans (hlt _ (List.getElem_mem _)).le
  have hB : ∀ d k (hk : k < L.length), k + d = L.length - 1 → hi L[k] ≤ hi L[L.length - 1] := by
    intro d
    induction d with
    | zero => intro k hk hkd; simp only [Nat.add_zero] at hkd; subst hkd; exact le_rfl
    | succ d ih =>
      intro k hk hkd
      have hk1 : k + 1 < L.length := by omega
      have := hlo (k + 1) hk1
      simp only [Nat.add_one_ne_zero, ite_false, Nat.add_sub_cancel] at this
      calc hi L[k] = lo L[k + 1] := this.symm
        _ ≤ hi L[k + 1] := (hlt _ (List.getElem_mem _)).le
        _ ≤ _ := ih (k + 1) hk1 (by omega)
  have hC : ∀ k (hk : k < L.length) x, a ≤ x → x < hi L[k] →
      ∃ k', ∃ hk' : k' < L.length, lo L[k'] ≤ x ∧ x < hi L[k'] := by
    intro k
    induction k with
    | zero => intro hk x hax hx; exact ⟨0, hk, by rw [hlo 0 hk]; simpa using hax, hx⟩
    | succ k ih =>
      intro hk x hax hx
      by_cases hxk : x < hi L[k]
      · exact ih (by omega) x hax hxk
      · refine ⟨k + 1, hk, ?_, hx⟩
        rw [hlo (k + 1) hk]; simp only [Nat.add_one_ne_zero, ite_false, Nat.add_sub_cancel]; omega
  refine ⟨(hA _ (by omega)).trans (hlt _ (List.getElem_mem _)).le, fun x => ⟨?_, ?_⟩⟩
  · rintro ⟨hax, hx⟩; exact hC _ (by omega) x hax hx
  · rintro ⟨k, hk, h1, h2⟩
    exact ⟨(hA k hk).trans h1, h2.trans_le (hB (L.length - 1 - k) k hk (by omega))⟩

omit h in
theorem ordered_perm (i : ItemId) (nvSt : Nat) (pos : Nat → Nat) :
    (items.ordered g i nvSt pos).Perm (items.ch i) := by
  unfold Items.ordered
  split_ifs
  · exact List.Perm.refl _
  · exact List.mergeSort_perm _ _

theorem below_iff_aux : ∀ n i, i < items.size → (items.desc i).card = n →
    idx i < t.subtreeEnd[idx i]! ∧
    ∀ j, j < items.size → (items.Below i j ↔ idx i ≤ idx j ∧ idx j < t.subtreeEnd[idx i]!) := by
  intro n
  induction n using Nat.strong_induction_on with
  | _ n ih =>
  intro i hi hn
  obtain ⟨pos, hl⟩ := (h.node i hi).layout
  have hmem : ∀ c, c ∈ items.ordered g i (t.nvRange (idx i)).1 pos ↔ c ∈ items.ch i :=
    fun c => (ordered_perm i _ pos).mem_iff
  have ihc : ∀ c ∈ items.ordered g i (t.nvRange (idx i)).1 pos,
      idx c < t.subtreeEnd[idx c]! ∧
      ∀ j, j < items.size → (items.Below c j ↔ idx c ≤ idx j ∧ idx j < t.subtreeEnd[idx c]!) :=
    fun c hc => ih _ (hn ▸ h.card_desc_child hi ((hmem c).1 hc)) c (h.ch_lt hi ((hmem c).1 hc)) rfl
  obtain ⟨hab, hcov⟩ := interval_cover (items.ordered g i (t.nvRange (idx i)).1 pos) idx
    (fun c => t.subtreeEnd[idx c]!) (idx i + 1) hl.child_idx (fun c hc => (ihc c hc).1)
    (t.subtreeEnd[idx i]!)
    (by rw [hl.subtree_end]; cases (items.ordered g i (t.nvRange (idx i)).1 pos).getLast? <;> rfl)
  refine ⟨by omega, fun j hj => ⟨?_, ?_⟩⟩
  · intro hb
    rcases Relation.ReflTransGen.cases_head hb with rfl | ⟨c, hc, hcj⟩
    · exact ⟨le_rfl, by omega⟩
    · have hcL := (hmem c).2 hc
      obtain ⟨k, hk, rfl⟩ := List.getElem_of_mem hcL
      have := (hcov (idx j)).2 ⟨k, hk, ((ihc _ hcL).2 j hj).1 hcj⟩
      omega
  · rintro ⟨h1, h2⟩
    rcases Nat.eq_or_lt_of_le h1 with heq | hlt
    · obtain rfl := h.idx_inj hi hj heq
      exact Relation.ReflTransGen.refl
    · obtain ⟨k, hk, hk1, hk2⟩ := (hcov (idx j)).1 ⟨hlt, h2⟩
      have hcL := List.getElem_mem hk
      exact Relation.ReflTransGen.head ((hmem _).1 hcL) (((ihc _ hcL).2 j hj).2 ⟨hk1, hk2⟩)

theorem below_iff {i j : ItemId} (hi : i < items.size) (hj : j < items.size) :
    items.Below i j ↔ idx i ≤ idx j ∧ idx j < t.subtreeEnd[idx i]! :=
  (h.below_iff_aux _ i hi rfl).2 j hj

theorem inSubtree_iff {i j : ItemId} (hi : i < items.size) (hj : j < items.size) :
    t.inSubtree (idx i) (idx j) ↔ items.Below i j := by
  rw [inSubtree_iff', h.below_iff hi hj]

theorem edgeIn_iff {i : ItemId} (hi : i < items.size) {e : Nat} (he : e < g.ne) :
    t.EdgeIn (idx i) e ↔ items.EdgeBelow g i e := by
  unfold SpqrTree.EdgeIn Items.EdgeBelow
  rw [h.ridx.edge_index e he]
  constructor
  · rintro ⟨j, hj, hs⟩
    cases hj
    exact (h.inSubtree_iff hi (h.edgeItem_lt he)).1 hs
  · intro hb
    exact ⟨_, rfl, (h.inSubtree_iff hi (h.edgeItem_lt he)).2 hb⟩

/-! ### Parent and children transport -/

theorem parent_eq_iff {p c : ItemId} (hp : p < items.size) (hc : c < items.size) :
    t.parent (idx c) = some (idx p) ↔ items.IsParent p c := by
  constructor
  · intro hpar
    by_cases hc0 : c = rootItem
    · subst hc0
      rw [h.ridx.root, h.ridx.root_par] at hpar
      cases hpar
    · obtain ⟨q, hq, -⟩ := h.tree.unique_parent c (by riomega) hc
      have := (h.node q (parent_lt hq)).child_par c hq
      rw [this] at hpar
      obtain rfl := h.idx_inj (parent_lt hq) hp (Option.some_injective _ hpar)
      exact hq
  · exact (h.node p hp).child_par c

theorem lt_size_of_parent {n m : Nat} (hp : t.parent n = some m) : n < t.size := by
  by_contra hge
  have : t.par[n]? = none := Array.getElem?_eq_none_iff.2 (by rw [h.ridx.sizes.par]; omega)
  simp [SpqrTree.parent, this] at hp

theorem parent_cases {n m : Nat} (hp : t.parent n = some m) :
    ∃ c p, c < items.size ∧ p < items.size ∧ idx c = n ∧ idx p = m ∧ items.IsParent p c := by
  obtain ⟨c, hc, rfl⟩ := h.idx_surj (h.lt_size_of_parent hp)
  by_cases hc0 : c = rootItem
  · subst hc0
    rw [h.ridx.root, h.ridx.root_par] at hp
    cases hp
  · obtain ⟨q, hq, -⟩ := h.tree.unique_parent c (by riomega) hc
    have := (h.node q (parent_lt hq)).child_par c hq
    rw [this] at hp
    exact ⟨c, q, hc, parent_lt hq, rfl, Option.some_injective _ hp, hq⟩

theorem mem_children_iff {i : ItemId} (hi : i < items.size) (n : Nat) :
    n ∈ t.children (idx i) ↔ ∃ c ∈ items.ch i, idx c = n := by
  obtain ⟨pos, hl⟩ := (h.node i hi).layout
  rw [hl.children, List.mem_map]
  simp only [(ordered_perm i _ pos).mem_iff]

/-! ### `Represents` fields that need only `Items.WF` -/

theorem canonical : ∀ i p, t.parent i = some p →
    (t.type i = .S → t.type p ≠ .S) ∧ (t.type i = .P → t.type p ≠ .P) := by
  intro n m hp
  obtain ⟨c, p, hc, hp', rfl, rfl, hpar⟩ := h.parent_cases hp
  rw [h.type_eq hc, h.type_eq hp']
  exact h.shapes.canonical p c hpar

theorem interior : ∀ i v, i < t.size → v < g.nv → ∀ j, t.vertIndex[v]! = some j →
    (t.parent j = some i ↔
      (∀ e, SpqrTree.Graph.Incident g v e → t.EdgeIn i e) ∧
      ∀ c ∈ t.children i, ¬ ∀ e, SpqrTree.Graph.Incident g v e → t.EdgeIn c e) := by
  intro n v hn hv j hj
  obtain ⟨a, ha, rfl⟩ := h.idx_surj hn
  rw [h.ridx.vert_index v hv] at hj
  cases hj
  rw [h.parent_eq_iff ha (h.vertItem_lt hv), h.endpoints.interior a v ha hv]
  constructor
  · rintro ⟨h1, h2⟩
    refine ⟨fun e ⟨he, hinc⟩ => (h.edgeIn_iff ha he).2 (h1 e he hinc), ?_⟩
    intro c' hc' hall
    obtain ⟨c, hc, rfl⟩ := (h.mem_children_iff ha c').1 hc'
    exact h2 c hc fun e he hinc => (h.edgeIn_iff (h.ch_lt ha hc) he).1 (hall e ⟨he, hinc⟩)
  · rintro ⟨h1, h2⟩
    refine ⟨fun e he hinc => (h.edgeIn_iff ha he).1 (h1 e ⟨he, hinc⟩), ?_⟩
    intro c hc hall
    exact h2 (idx c) ((h.mem_children_iff ha _).2 ⟨c, hc, rfl⟩)
      fun e ⟨he, hinc⟩ => (h.edgeIn_iff (h.ch_lt ha hc) he).2 (hall e he hinc)

/-- `nv_orig_inj`, from `Items.nvList` having no repeats. -/
theorem nv_orig_inj (hnd : ∀ i, i < items.size → (items.nvList g i).Nodup) :
    ∀ i, i < t.size → t.type i ≠ .O → t.type i ≠ .Q →
      ((t.nodeVertsOf i).map fun nv => t.origId[nv.vert]!).Nodup := by
  intro n hn _ _
  obtain ⟨a, ha, rfl⟩ := h.idx_surj hn
  rw [h.nv_orig_map ha]
  exact (hnd a ha).map (Option.some_injective _)

/-! ### Node-edge transport -/

omit h in
theorem isNode_iff (ty : NodeType) : ty.isNode = true ↔ ty ∉ [NodeType.F, .V] := by
  cases ty <;> simp [NodeType.isNode]

theorem filter_lt_eq {i : ItemId} (hi : i < items.size) :
    (items.ch i).filter (· < 1 + g.nv) = (items.ch i).filter fun c => items.type c = .V := by
  apply List.filter_congr
  intro c hc
  exact decide_eq_decide.2 (h.type_V_iff hi hc).symm

theorem nv_nodup {i : ItemId} (hi : i < items.size) : (items.nvList g i).Nodup := by
  unfold Items.nvList; rw [h.filter_lt_eq hi]; exact h.endpoints.nv_nodup i hi

theorem filter_ge_eq {i : ItemId} (hi : i < items.size) :
    (items.ch i).filter (· ≥ 1 + g.nv) = (items.ch i).filter fun c => items.type c ≠ .V := by
  apply List.filter_congr
  intro c hc
  exact decide_eq_decide.2 (by rw [ne_eq, h.type_V_iff hi hc]; constructor <;> intro <;> riomega)

theorem virt_length {i : ItemId} (hi : i < items.size) :
    (items.virtualEdges i).length = (items.ch i).countP (· ≥ 1 + g.nv) := by
  rw [Items.virtualEdges, List.length_map, ← h.filter_ge_eq hi, List.countP_eq_length_filter]

theorem nEdges_node {i : ItemId} (hi : i < items.size) (hn : (items.type i).isNode = true) :
    items.nEdges g i = (items.virtualEdges i).length + items.capCount i := by
  rw [Items.nEdges, ite_eq_left hn, h.virt_length hi]

omit h in
theorem nEdges_not_node {i : ItemId} (hn : (items.type i).isNode = false) : items.nEdges g i = 0 := by
  simp [Items.nEdges, Items.capCount, Items.hasCap, hn]

omit h in
theorem ordered_filter_perm (i : ItemId) (nvSt : Nat) (pos : Nat → Nat) (p : ItemId → Bool) :
    ((items.ordered g i nvSt pos).filter p).Perm ((items.ch i).filter p) :=
  (ordered_perm i nvSt pos).filter p

theorem type_Q_eq {i : ItemId} (hi : i < items.size) (hQ : items.type i = .Q) :
    ∃ e, e < g.ne ∧ i = edgeItem g e := by
  by_cases h1 : i < 1 + g.nv
  · exfalso
    rcases Nat.eq_zero_or_pos i with h0 | h0
    · subst h0; have := h.tree.root; rw [show (0 : ItemId) = rootItem from rfl, this] at hQ; cases hQ
    · have : i = vertItem (i - 1) := by riomega
      rw [this, h.type_vertItem (by riomega)] at hQ; cases hQ
  · by_cases h2 : i < 1 + g.nv + g.ne
    · exact ⟨i - (1 + g.nv), by riomega, by riomega⟩
    · exfalso
      have := h.tree.node i (by riomega) hi
      simp [hQ] at this

theorem hasCap_eq {i : ItemId} (hi : i < items.size) : t.hasCap (idx i) = items.hasCap i := by
  obtain ⟨pos, hl⟩ := (h.node i hi).layout
  have hperm := ordered_perm (items := items) (g := g) i (t.nvRange (idx i)).1 pos
  by_cases hQ : items.type i = .Q
  · obtain ⟨e, he, rfl⟩ := h.type_Q_eq hi hQ
    have key : (((items.ordered g (edgeItem g e) (t.nvRange (idx (edgeItem g e))).1 pos).map idx).any
        fun c => t.type c != .V) = !(items.ch (edgeItem g e)).isEmpty := by
      rcases h.shapes.q_children e he with h0 | ⟨c, hc, hch⟩
      · rw [h0] at hperm ⊢
        rw [List.perm_nil.1 hperm]; rfl
      · have hne : items.ch (edgeItem g e) ≠ [] := by
          rcases hch with h1 | ⟨v, h1⟩ <;> simp [h1]
        have hcm : c ∈ items.ch (edgeItem g e) := by
          rcases hch with h1 | ⟨v, h1⟩ <;> simp [h1]
        rw [Bool.eq_iff_iff, List.any_eq_true]
        constructor
        · intro _; simpa using hne
        · intro _
          refine ⟨idx c, List.mem_map_of_mem (hperm.mem_iff.2 hcm), ?_⟩
          rw [h.type_eq (h.ch_lt hi hcm)]
          simpa using fun hV => hc (by simp [hV])
    rw [SpqrTree.hasCap, Items.hasCap, h.type_eq hi, hl.children, key]
  · have : (items.type (i) == .Q) = false := by simpa using hQ
    rw [SpqrTree.hasCap, Items.hasCap, h.type_eq hi, this]
    simp

theorem neBounds_cover : ∀ n, n ≤ t.size → ∀ ne, ne < t.neBounds[n]! →
    ∃ m, m < n ∧ t.neBounds[m]! ≤ ne ∧ ne < t.neBounds[m + 1]! := by
  intro n
  induction n with
  | zero => intro _ ne hne; rw [h.ridx.ne_zero] at hne; omega
  | succ n ih =>
    intro hn ne hne
    by_cases hlt : ne < t.neBounds[n]!
    · obtain ⟨m, hm, h1, h2⟩ := ih (by omega) ne hlt
      exact ⟨m, by omega, h1, h2⟩
    · exact ⟨n, by omega, by omega, hne⟩

theorem ne_mem_range {ne : Nat} (hne : ne < t.nodeEdges.size) :
    ∃ a, a < items.size ∧ (t.neRange (idx a)).1 ≤ ne ∧ ne < (t.neRange (idx a)).2 := by
  obtain ⟨m, hm, h1, h2⟩ := h.neBounds_cover t.size le_rfl ne (by rw [h.ridx.ne_last]; exact hne)
  obtain ⟨a, ha, rfl⟩ := h.idx_surj hm
  exact ⟨a, ha, by rw [neRange_fst]; exact h1, by rw [neRange_snd]; exact h2⟩

omit h in
theorem twin_some_lt {ne ne' : Nat} (ht : t.twin ne = some ne') : ne < t.nodeEdges.size := by
  by_contra hge
  simp [SpqrTree.twin, Array.getElem?_eq_none_iff.2 (by omega : t.nodeEdges.size ≤ ne)] at ht

theorem nodeEdges_get {i : ItemId} (hi : i < items.size) {k : Nat} (hk : k < items.nEdges g i) :
    (t.neRange (idx i)).1 + k < t.nodeEdges.size ∧
    t.nodeEdges[(t.neRange (idx i)).1 + k]? = some (t.nodeEdges[(t.neRange (idx i)).1 + k]!) := by
  have hlt : (t.neRange (idx i)).1 + k < t.nodeEdges.size := by
    have := (h.node i hi).ne_range; have := h.neEn_le hi; omega
  exact ⟨hlt, by rw [Array.getElem?_eq_getElem hlt, ← getElem!_pos]⟩

theorem neOrig_eq {i : ItemId} (hi : i < items.size) {k a b : Nat} (hk : k < items.nEdges g i)
    (ha : a < (items.nvList g i).length) (hb : b < (items.nvList g i).length)
    (hnvs : (t.nodeEdges[(t.neRange (idx i)).1 + k]!).nvs =
      ((t.nvRange (idx i)).1 + a, (t.nvRange (idx i)).1 + b)) :
    t.neOrig ((t.neRange (idx i)).1 + k) = some ((items.nvList g i)[a], (items.nvList g i)[b]) := by
  obtain ⟨-, hget⟩ := h.nodeEdges_get hi hk
  have h1 := h.nvOrig_get hi ha
  have h2 := h.nvOrig_get hi hb
  simp [SpqrTree.neOrig, SpqrTree.nvsOf, hget, hnvs, h1, h2]

/-! ### Endpoints and node-vert lists per type -/

theorem vs_fst {i : ItemId} (hi : i < items.size) (hn : (items.type i).isNode = true) :
    ∃ u, (items.vs i).1 = some u := by
  have hs := h.endpoints.vs_shape i hi
  cases hty : items.type i <;> rw [hty] at hs hn <;> simp [NodeType.isNode] at hn
  · exact hs.elim (fun h1 => h1.elim fun v hv => ⟨v, by rw [hv]⟩)
      (fun h2 => h2.elim fun u h2' => h2'.elim fun v hv => ⟨u, by rw [hv]⟩)
  · obtain ⟨u, v, hv⟩ := hs; exact ⟨u, by rw [hv]⟩
  · obtain ⟨v, hv⟩ := hs; exact ⟨v, by rw [hv]⟩
  all_goals obtain ⟨u, v, hv⟩ := hs; exact ⟨u, by rw [hv]⟩

omit h in
theorem nvList_eq (i : ItemId) : items.nvList g i =
    (items.vs i).1.toList ++ ((items.ch i).filter (· < 1 + g.nv)).map (· - 1) ++ (items.vs i).2.toList := rfl

omit h in
theorem vs_head {i : ItemId} {u : Nat} (hu : (items.vs i).1 = some u) :
    (items.nvList g i)[0]? = some u := by
  rw [nvList_eq, hu]; rfl

omit h in
theorem vs_last {i : ItemId} {v : Nat} (hv : (items.vs i).2 = some v) :
    (items.nvList g i)[(items.nvList g i).length - 1]? = some v := by
  rw [nvList_eq, hv, Option.toList_some]
  set A := (items.vs i).1.toList ++ ((items.ch i).filter (· < 1 + g.nv)).map (· - 1)
  have : (A ++ [v]).length - 1 = A.length := by simp
  rw [this, List.getElem?_concat_length]

omit h in
theorem nvList_pos {i : ItemId} {u : Nat} (hu : (items.vs i).1 = some u) : 0 < (items.nvList g i).length := by
  rw [nvList_eq, hu]; simp

theorem filter_V_eq_nil {i : ItemId} (hi : i < items.size)
    (hV : ∀ c, items.IsParent i c → items.type c ≠ .V) :
    ((items.ch i).filter (· < 1 + g.nv)).map (· - 1) = [] := by
  rw [h.filter_lt_eq hi, List.filter_eq_nil_iff.2 (fun c hc => by simpa using hV c hc)]; rfl

omit h in
theorem map_vertItem_sub (xs : List Nat) : (xs.map vertItem).map (· - 1) = xs := by
  rw [List.map_map]
  conv_rhs => rw [← List.map_id xs]
  apply List.map_congr_left
  intro x _; show vertItem x - 1 = x; riomega

theorem nvList_O {i : ItemId} (hi : i < items.size) (hO : items.type i = .O) :
    ∃ v, items.vs i = (some v, none) ∧ items.nvList g i = [v] := by
  have hs := h.endpoints.vs_shape i hi
  rw [hO] at hs
  obtain ⟨v, hv⟩ := hs
  refine ⟨v, hv, ?_⟩
  rw [nvList_eq, hv, h.shapes.i_o_leaf i hi (Or.inr hO)]; rfl

theorem nvList_I {i : ItemId} (hi : i < items.size) (hI : items.type i = .I) :
    ∃ u v, items.vs i = (some u, some v) ∧ items.nvList g i = [u, v] := by
  have hs := h.endpoints.vs_shape i hi
  rw [hI] at hs
  obtain ⟨u, v, hv⟩ := hs
  refine ⟨u, v, hv, ?_⟩
  rw [nvList_eq, hv, h.shapes.i_o_leaf i hi (Or.inl hI)]; rfl

theorem nvList_P {i : ItemId} (hi : i < items.size) (hP : items.type i = .P) :
    ∃ u v, items.vs i = (some u, some v) ∧ items.nvList g i = [u, v] := by
  have hs := h.endpoints.vs_shape i hi
  rw [hP] at hs
  obtain ⟨u, v, hv⟩ := hs
  refine ⟨u, v, hv, ?_⟩
  rw [nvList_eq, hv, h.filter_V_eq_nil hi (h.shapes.p_shape i hi hP).2.1]; rfl

theorem nvList_S (hr : items.RepOK g) {i : ItemId} (hi : i < items.size) (hS : items.type i = .S) :
    ∃ u v xs, items.vs i = (some u, some v) ∧ items.nvList g i = u :: xs ++ [v] ∧ 2 ≤ xs.length ∧
      items.virtualEdges i = List.zip (u :: xs) (xs ++ [v]) := by
  obtain ⟨u, v, xs, hv, hxs, hvirt⟩ := hr.s_order i hi hS
  obtain ⟨u', v', xs', hv', hxs', hlen, -⟩ := h.shapes.s_shape i hi hS
  rw [hv] at hv'; cases hv'
  rw [hxs] at hxs'
  have : xs = xs' := by
    have := congrArg (List.map (· - 1)) hxs'
    rwa [map_vertItem_sub, map_vertItem_sub] at this
  subst this
  refine ⟨u, v, xs, hv, ?_, hlen, hvirt⟩
  rw [nvList_eq, hv, h.filter_lt_eq hi, hxs, map_vertItem_sub]; rfl

theorem nvList_R {i : ItemId} (hi : i < items.size) (hR : items.type i = .R) :
    ∃ u v, items.vs i = (some u, some v) ∧ 4 ≤ (items.nvList g i).length := by
  have hs := h.endpoints.vs_shape i hi
  rw [hR] at hs
  obtain ⟨u, v, hv⟩ := hs
  refine ⟨u, v, hv, ?_⟩
  have := (h.shapes.r_shape i hi hR).1
  rw [nvList_eq, hv, h.filter_lt_eq hi]
  simp only [Option.toList_some, List.length_append, List.length_cons, List.length_nil, List.length_map]
  omega

/-- A Q item: its endpoints `u`, `v` (`v = u` for a loop), its node-vert list, and its edge. -/
theorem q_root_iff {e : Nat} (he : e < g.ne) :
    (items.vs (edgeItem g e)).2 = none ↔ items.ch (edgeItem g e) ≠ [] := by
  have hs := h.endpoints.vs_shape _ (h.edgeItem_lt he)
  rw [h.type_edgeItem he] at hs
  obtain ⟨u, hu⟩ : ∃ u, (items.vs (edgeItem g e)).1 = some u := by
    rcases hs with ⟨v, hv⟩ | ⟨u, v, hv⟩
    · exact ⟨v, by rw [hv]⟩
    · exact ⟨u, by rw [hv]⟩
  exact (h.endpoints.q_vs e he u hu).2.1

theorem nEdges_Q {e : Nat} (he : e < g.ne) : items.nEdges g (edgeItem g e) = 1 := by
  have hq := h.edgeItem_lt he
  have hQ := h.type_edgeItem he
  rw [h.nEdges_node hq (by rw [hQ]; decide), Items.virtualEdges, Items.capCount]
  rcases h.shapes.q_children e he with h0 | ⟨c, hc, hch⟩
  · simp [h0, Items.hasCap, hQ, NodeType.isNode]
  · have hcV : items.type c ≠ .V := fun hV => hc (by simp [hV])
    rcases hch with h1 | ⟨w, hw, h1⟩
    · simp [h1, Items.hasCap, hQ, hcV]
    · simp [h1, Items.hasCap, hQ, hcV, h.type_vertItem hw]

theorem nvList_Q_le {e : Nat} (he : e < g.ne) : (items.nvList g (edgeItem g e)).length ≤ 2 := by
  have hq := h.edgeItem_lt he
  have hs := h.endpoints.vs_shape _ hq
  rw [h.type_edgeItem he] at hs
  unfold Items.nvList
  rw [h.filter_lt_eq hq]
  rcases h.shapes.q_children e he with h0 | ⟨c, hc, hch⟩
  · rw [h0]; rcases hs with ⟨v, hv⟩ | ⟨u, v, hv⟩ <;> simp [hv]
  · have hcV : items.type c ≠ .V := fun hV => hc (by simp [hV])
    have hv2 := (h.q_root_iff he).2 (by rcases hch with h1 | ⟨w, -, h1⟩ <;> simp [h1])
    rcases hs with ⟨v, hv⟩ | ⟨u, v, hv⟩
    · rcases hch with h1 | ⟨w, hw, h1⟩
      · simp [hv, h1, hcV]
      · simp [hv, h1, hcV, h.type_vertItem hw]
    · rw [hv] at hv2; cases hv2

theorem q_pair (hr : items.RepOK g) {e : Nat} (he : e < g.ne) :
    ∃ u v, (items.vs (edgeItem g e)).1 = some u ∧ Items.PairEq (u, v) g.edges[e]! ∧
      ((items.ch (edgeItem g e) = [] ∧ items.vs (edgeItem g e) = (some u, some v) ∧
          items.nvList g (edgeItem g e) = [u, v]) ∨
        (items.ch (edgeItem g e) ≠ [] ∧ items.vs (edgeItem g e) = (some u, none) ∧ ∃ c,
          (items.ch (edgeItem g e) = [c, vertItem v] ∧ items.nvList g (edgeItem g e) = [u, v]) ∨
          (items.ch (edgeItem g e) = [c] ∧ items.nvList g (edgeItem g e) = [u] ∧ v = u))) ∧
      (items.ch (edgeItem g e) = [] → items.hasCap (edgeItem g e) = true) ∧
      items.nEdges g (edgeItem g e) = 1 := by
  have hq := h.edgeItem_lt he
  have hQ := h.type_edgeItem he
  have hs := h.endpoints.vs_shape _ hq
  rw [hQ] at hs
  obtain ⟨hloop, hpair⟩ := hr.q_root e he
  have hroot := h.q_root_iff he
  have hcap0 : items.ch (edgeItem g e) = [] → items.hasCap (edgeItem g e) = true := by
    intro h0; simp [Items.hasCap, hQ, h0, NodeType.isNode]
  rcases h.shapes.q_children e he with h0 | ⟨c, hc, hch⟩
  · -- leaf
    have hcap := hcap0 h0
    have hne : items.nEdges g (edgeItem g e) = 1 := by
      rw [h.nEdges_node hq (by rw [hQ]; decide), Items.virtualEdges, h0, Items.capCount, hcap]; rfl
    rcases hs with ⟨v, hv⟩ | ⟨u, v, hv⟩
    · exfalso
      have : (items.vs (edgeItem g e)).2 = none := by rw [hv]
      exact (hroot.1 this) h0
    · have hu : (items.vs (edgeItem g e)).1 = some u := by rw [hv]
      refine ⟨u, v, hu, (h.endpoints.q_vs e he u hu).2.2.2 v (by rw [hv]), Or.inl ⟨h0, hv, ?_⟩, hcap0, hne⟩
      rw [nvList_eq, hv, h0]; rfl
  · have hne0 : items.ch (edgeItem g e) ≠ [] := by
      rcases hch with h1 | ⟨w, -, h1⟩ <;> simp [h1]
    have hv2 : (items.vs (edgeItem g e)).2 = none := hroot.2 hne0
    have hcapf : items.hasCap (edgeItem g e) = false := by
      simp [Items.hasCap, hQ, hne0]
    have hcV : items.type c ≠ .V := fun hV => hc (by simp [hV])
    have hcm : c ∈ items.ch (edgeItem g e) := by
      rcases hch with h1 | ⟨w, -, h1⟩ <;> simp [h1]
    have hcge : ¬ c < 1 + g.nv := fun hlt => hcV ((h.type_V_iff hq hcm).2 hlt)
    obtain ⟨u, hu⟩ : ∃ u, (items.vs (edgeItem g e)).1 = some u := by
      rcases hs with ⟨v, hv⟩ | ⟨u, v, hv⟩
      · exact ⟨v, by rw [hv]⟩
      · exfalso; rw [hv] at hv2; cases hv2
    have hvs : items.vs (edgeItem g e) = (some u, none) := by
      ext1 <;> simp [hu, hv2]
    have hne0' := h.nEdges_node hq (by rw [hQ]; decide)
    rw [Items.capCount, hcapf, Items.virtualEdges] at hne0'
    rcases hch with h1 | ⟨w, hw, h1⟩
    · have hne : items.nEdges g (edgeItem g e) = 1 := by rw [hne0', h1]; simp [hcV]
      have hl := hloop c h1
      refine ⟨u, u, hu, ?_, Or.inr ⟨hne0, hvs, c, Or.inr ⟨h1, ?_, rfl⟩⟩, hcap0, hne⟩
      · have := (h.endpoints.q_vs e he u hu).1
        left; ext <;> simp <;> omega
      · rw [nvList_eq, hvs, h1]; simp [hcge]
    · have hpw := hpair c w u h1 hu
      have hwV := h.type_vertItem hw
      have hne : items.nEdges g (edgeItem g e) = 1 := by rw [hne0', h1]; simp [hcV, hwV]
      have hwlt : vertItem w < 1 + g.nv := by riomega
      refine ⟨u, w, hu, hpw, Or.inr ⟨hne0, hvs, c, Or.inl ⟨h1, ?_⟩⟩, hcap0, hne⟩
      rw [nvList_eq, hvs, h1]
      simp [hcge, hwlt]
      show 1 + w - 1 = w; omega

/-- The cap (first node-edge) of a node joins its first and last node-vert. -/
theorem edge0_nvs (hr : items.RepOK g) {i : ItemId} (hi : i < items.size)
    (hn : (items.type i).isNode = true) (h1 : 1 ≤ items.nEdges g i) :
    (t.nodeEdges[(t.neRange (idx i)).1]!).nvs =
      ((t.nvRange (idx i)).1, (t.nvRange (idx i)).1 + ((items.nvList g i).length - 1)) := by
  obtain ⟨pos, hl⟩ := (h.node i hi).layout
  have hnv := (h.node i hi).nv_range
  have hne := (h.node i hi).ne_range
  have hE : (items.edgeChildren g pos (items.ordered g i (t.nvRange (idx i)).1 pos)).length =
      (items.virtualEdges i).length := by
    rw [Items.edgeChildren, List.length_map, (ordered_filter_perm i _ pos _).length_eq,
      h.virt_length hi, List.countP_eq_length_filter]
  have hnE := h.nEdges_node hi hn
  have := hl.edge_nvs 0 h1
  rw [Nat.add_zero] at this
  rw [this, nodeLayout, hnv, hne]
  generalize (t.nvRange (idx i)).1 = nvSt at hl hE ⊢
  generalize (t.neRange (idx i)).1 = neSt at hl ⊢
  generalize items.edgeChildren g pos (items.ordered g i nvSt pos) = E at hE ⊢
  cases hty : items.type i
  · rw [hty] at hn; cases hn
  · rw [hty] at hn; cases hn
  · obtain ⟨e, he, rfl⟩ := h.type_Q_eq hi hty
    obtain ⟨u, v, -, -, hl, -, hne1⟩ := h.q_pair hr he
    rcases hl with ⟨-, -, hl2⟩ | ⟨-, -, c, ⟨-, hl2⟩ | ⟨-, hl1, -⟩⟩
    · rw [hl2]; exact LayoutFacts.qi_edge _ (Or.inl rfl) _ _ _ _ _ (by omega)
    · rw [hl2]; exact LayoutFacts.qi_edge _ (Or.inl rfl) _ _ _ _ _ (by omega)
    · rw [hl1]; exact LayoutFacts.one_vert _ (by decide) _ _ _ _ _ (by omega)
  · have hne1 : items.nEdges g i = 1 := by
      rw [h.nEdges_node hi hn, Items.virtualEdges, h.shapes.i_o_leaf i hi (Or.inl hty)]
      simp [Items.capCount, Items.hasCap, hty, NodeType.isNode]
    obtain ⟨u, v, -, hl2⟩ := h.nvList_I hi hty
    rw [hl2]; exact LayoutFacts.qi_edge _ (Or.inr rfl) _ _ _ _ _ (by omega)
  · have hne1 : items.nEdges g i = 1 := by
      rw [h.nEdges_node hi hn, Items.virtualEdges, h.shapes.i_o_leaf i hi (Or.inr hty)]
      simp [Items.capCount, Items.hasCap, hty, NodeType.isNode]
    obtain ⟨v, -, hl1⟩ := h.nvList_O hi hty
    rw [hl1]; exact LayoutFacts.one_vert _ (by decide) _ _ _ _ _ (by omega)
  · obtain ⟨u, v, xs, -, hnl, hxs, hvirt⟩ := h.nvList_S hr hi hty
    have hlen : (items.nvList g i).length = xs.length + 2 := by rw [hnl]; simp
    have hcap : items.capCount i = 1 := by simp [Items.capCount, Items.hasCap, hty, NodeType.isNode]
    have hvl : (items.virtualEdges i).length = xs.length + 1 := by rw [hvirt, List.length_zip]; simp
    rw [hvl, hcap] at hnE
    have := LayoutFacts.s_edge (idx i) nvSt (nvSt + (items.nvList g i).length) neSt
      (neSt + items.nEdges g i) E (by omega) (by omega) 0 (by omega)
    rw [ite_eq_left rfl] at this
    rw [this]; congr 1; omega
  · obtain ⟨u, v, -, hl2⟩ := h.nvList_P hi hty
    rw [hl2]; exact LayoutFacts.p_edge _ _ _ _ _ 0 (by omega)
  · obtain ⟨u, v, -, hlen⟩ := h.nvList_R hi hty
    have hcap : items.capCount i = 1 := by simp [Items.capCount, Items.hasCap, hty, NodeType.isNode]
    rw [hcap] at hnE
    have := LayoutFacts.r_edge (idx i) nvSt (nvSt + (items.nvList g i).length) neSt
      (neSt + items.nEdges g i) E (by omega) (by omega) 0 (by omega)
    rw [ite_eq_left rfl] at this
    rw [this]; congr 1; omega

theorem edges_lt (hr : items.RepOK g) {e : Nat} (he : e < g.ne) :
    (g.edges[e]!).1 < g.nv ∧ (g.edges[e]!).2 < g.nv := by
  obtain ⟨u, v, hu, hpe, hl, -, -⟩ := h.q_pair hr he
  have hu' := h.endpoints.vs_lt _ u (Or.inl hu)
  have hv' : v < g.nv := by
    rcases hl with ⟨-, -, hl⟩ | ⟨-, -, c, ⟨-, hl⟩ | ⟨-, hl, rfl⟩⟩
    · exact h.nvList_lt (h.edgeItem_lt he) (by rw [hl]; simp)
    · exact h.nvList_lt (h.edgeItem_lt he) (by rw [hl]; simp)
    · exact hu'
  rcases hpe with hpe | hpe
  · rw [← hpe]; exact ⟨hu', hv'⟩
  · obtain ⟨h1, h2⟩ := Prod.mk.inj hpe; omega

theorem q_endpoints (hr : items.RepOK g) : ∀ e, e < g.ne → ∀ i, t.edgeIndex[e]! = some i →
    ∀ ne, (t.neRange i).1 = ne → ∀ p, t.neOrig ne = some p →
      (t.edgeFlipped[e]! = false → p = g.edges[e]!) ∧
      (t.edgeFlipped[e]! = true → p = ((g.edges[e]!).2, (g.edges[e]!).1)) := by
  intro e he n hn ne hne p hp
  rw [h.ridx.edge_index e he] at hn
  cases hn; subst hne
  have hq := h.edgeItem_lt he
  have hQ := h.type_edgeItem he
  obtain ⟨u, v, hu, hpe, hl, -, hne1⟩ := h.q_pair hr he
  have h1 : 1 ≤ items.nEdges g (edgeItem g e) := by omega
  have hpos := nvList_pos (g := g) hu
  have hnvs := h.edge0_nvs hr hq (by rw [hQ]; decide) h1
  have hp' := h.neOrig_eq hq (k := 0) (a := 0) (b := (items.nvList g (edgeItem g e)).length - 1) h1
    (by omega) (by omega) (by rw [Nat.add_zero]; exact hnvs)
  rw [Nat.add_zero] at hp'
  rw [hp'] at hp
  cases hp
  have hpv : ((items.nvList g (edgeItem g e))[0]'(by omega),
      (items.nvList g (edgeItem g e))[(items.nvList g (edgeItem g e)).length - 1]'(by omega)) = (u, v) := by
    rcases hl with ⟨-, -, hl⟩ | ⟨-, -, c, ⟨-, hl⟩ | ⟨-, hl, rfl⟩⟩ <;> simp [hl]
  rw [hpv]
  have hfl := ((h.node _ hq).edge_index hQ).2
  have hidx : edgeItem g e - 1 - g.nv = e := by riomega
  rw [hidx] at hfl
  rw [hfl, hu]
  constructor
  · intro hf
    have hue : u = (g.edges[e]!).1 := by simpa using hf
    rcases hpe with hpe | hpe
    · exact hpe
    · obtain ⟨h1, h2⟩ := Prod.mk.inj hpe
      exact Prod.ext (by show u = _; omega) (by show v = _; omega)
  · intro hf
    have hue : u ≠ (g.edges[e]!).1 := by simpa using hf
    rcases hpe with hpe | hpe
    · exact absurd (Prod.mk.inj hpe).1 hue
    · exact hpe

theorem separation (hr : items.RepOK g) : ∀ i ne, i < t.size → t.capNe i = some ne →
    ∀ p, t.neOrig ne = some p →
    ∀ v e e', SpqrTree.Graph.Incident g v e → SpqrTree.Graph.Incident g v e' →
      t.EdgeIn i e → ¬ t.EdgeIn i e' → v = p.1 ∨ v = p.2 := by
  intro n ne hn hcap p hp v e e' ⟨he, hinc⟩ ⟨he', hinc'⟩ hin hnin
  obtain ⟨a, ha, rfl⟩ := h.idx_surj hn
  rw [SpqrTree.capNe, h.hasCap_eq ha] at hcap
  by_cases hc : items.hasCap a = true
  swap
  · rw [ite_eq_right hc] at hcap; cases hcap
  rw [ite_eq_left hc, Option.some.injEq] at hcap
  subst hcap
  have hn' : (items.type a).isNode = true := by
    simp only [Items.hasCap, Bool.and_eq_true] at hc; exact hc.1
  have hv : v < g.nv := by
    have := h.edges_lt hr he; rcases hinc with rfl | rfl <;> omega
  have hcc : items.capCount a = 1 := by simp [Items.capCount, hc]
  have h1 : 1 ≤ items.nEdges g a := by rw [h.nEdges_node ha hn', hcc]; omega
  obtain ⟨u, hu⟩ := h.vs_fst ha hn'
  have hpos := nvList_pos (g := g) hu
  have hnvs := h.edge0_nvs hr ha hn' h1
  have hp' := h.neOrig_eq ha (k := 0) (a := 0) (b := (items.nvList g a).length - 1) h1
    (by omega) (by omega) (by rw [Nat.add_zero]; exact hnvs)
  rw [Nat.add_zero] at hp'
  rw [hp'] at hp
  cases hp
  have hsep := h.endpoints.separation a ha ((isNode_iff _).1 hn') v e e' hv he he' hinc hinc'
    ((h.edgeIn_iff ha he).1 hin) (fun hb => hnin ((h.edgeIn_iff ha he').2 hb))
  rcases hsep with hs | hs
  · left; exact ((List.getElem_eq_iff _).2 (vs_head hs)).symm
  · right; exact ((List.getElem_eq_iff _).2 (vs_last hs)).symm

/-! ### Twin glue -/

omit h in
theorem PairEq.symm' {p q : Nat × Nat} (hpq : Items.PairEq p q) : Items.PairEq q p := by
  rcases hpq with rfl | rfl
  · exact Or.inl rfl
  · right; rfl

omit h in
theorem ne_of_nodup_head_last {x y : Nat} {l : List Nat} (hnd : (x :: l ++ [y]).Nodup) : x ≠ y := by
  rintro rfl
  simp at hnd

theorem child_type_ne_F {i c : ItemId} (hi : i < items.size) (hc : c ∈ items.ch i) : items.type c ≠ .F := by
  intro hF
  have hcs := h.ch_lt hi hc
  have h0 := h.ch_ne_root hc
  by_cases h1 : c < 1 + g.nv
  · exact absurd ((h.ch_lt_vert hc h1).2.symm.trans hF) (by decide)
  · by_cases h2 : c < 1 + g.nv + g.ne
    · have : c = edgeItem g (c - (1 + g.nv)) := by riomega
      rw [this, h.type_edgeItem (by riomega)] at hF; cases hF
    · have := h.tree.node c (by riomega) hcs
      simp [hF] at this

/-- The cap of a capped node joins `vs.1` to `vs.2` (or to itself for a one-vertex node). -/
theorem cap_orig (hr : items.RepOK g) {i : ItemId} (hi : i < items.size)
    (hn : (items.type i).isNode = true) (hc : items.hasCap i = true) :
    ∃ x, (items.vs i).1 = some x ∧ t.neOrig (t.neRange (idx i)).1 = some (x, (items.vs i).2.getD x) := by
  have hcc : items.capCount i = 1 := by simp [Items.capCount, hc]
  have h1 : 1 ≤ items.nEdges g i := by rw [h.nEdges_node hi hn, hcc]; omega
  obtain ⟨x, hx⟩ := h.vs_fst hi hn
  have hpos := nvList_pos (g := g) hx
  have hnvs := h.edge0_nvs hr hi hn h1
  have hp := h.neOrig_eq hi (k := 0) (a := 0) (b := (items.nvList g i).length - 1) h1
    (by omega) (by omega) (by rw [Nat.add_zero]; exact hnvs)
  rw [Nat.add_zero] at hp
  refine ⟨x, hx, ?_⟩
  rw [hp, (List.getElem_eq_iff _).2 (vs_head hx)]
  cases hy : (items.vs i).2 with
  | some y => rw [(List.getElem_eq_iff _).2 (vs_last hy)]; rfl
  | none =>
    have hs := h.endpoints.vs_shape i hi
    cases hty : items.type i <;> rw [hty] at hs hn
    · cases hn
    · cases hn
    · obtain ⟨e, he, rfl⟩ := h.type_Q_eq hi hty
      have hne : items.ch (edgeItem g e) ≠ [] := (h.q_root_iff he).1 hy
      simp [Items.hasCap, hty, hne] at hc
    · obtain ⟨u, v, hv⟩ := hs; rw [hv] at hy; cases hy
    · obtain ⟨v, hv, hl⟩ := h.nvList_O hi hty
      rw [hv] at hx; cases hx
      simp [hl]
    all_goals obtain ⟨u, v, hv⟩ := hs; rw [hv] at hy; cases hy

/-- A non-V child of a node has a cap joining its two endpoints (`O` children: a loop). -/
theorem child_cap_orig (hr : items.RepOK g) {p c : ItemId} (hp : p < items.size)
    (hn : (items.type p).isNode = true) (hc : c ∈ items.ch p) (hcV : items.type c ≠ .V) :
    ∃ x y, (items.vs c).1 = some x ∧ (items.vs c).2.getD x = y ∧
      t.neOrig (t.neRange (idx c)).1 = some (x, y) ∧ (items.type c ≠ .O → (items.vs c).2 = some y) := by
  have hcs := h.ch_lt hp hc
  have hcF := h.child_type_ne_F hp hc
  have hcn : (items.type c).isNode = true := by
    rw [isNode_iff]; simp [hcF, hcV]
  have hcap : items.hasCap c = true := by
    by_cases hQ : items.type c = .Q
    · obtain ⟨e, he, rfl⟩ := h.type_Q_eq hcs hQ
      have h0 : items.ch (edgeItem g e) = [] :=
        h.shapes.q_leaf_of_node p _ hc ((isNode_iff _).1 hn) hQ
      simp [Items.hasCap, hQ, h0, NodeType.isNode]
    · have : (items.type c == .Q) = false := by simpa using hQ
      simp [Items.hasCap, hcn, this]
  obtain ⟨x, hx, horig⟩ := h.cap_orig hr hcs hcn hcap
  refine ⟨x, _, hx, rfl, horig, ?_⟩
  intro hO
  cases hy : (items.vs c).2 with
  | some y => rfl
  | none =>
    exfalso
    have hs := h.endpoints.vs_shape c hcs
    cases hty : items.type c <;> rw [hty] at hs hcn hcF hcV hO
    · exact hcF rfl
    · exact hcV rfl
    · obtain ⟨e, he, rfl⟩ := h.type_Q_eq hcs hty
      exact (h.q_root_iff he).1 hy (h.shapes.q_leaf_of_node p _ hc ((isNode_iff _).1 hn) hty)
    · obtain ⟨u, v, hv⟩ := hs; rw [hv] at hy; cases hy
    · exact hO rfl
    all_goals obtain ⟨u, v, hv⟩ := hs; rw [hv] at hy; cases hy

/-- Endpoints of a non-V child are node-verts of the parent. -/
theorem mem_nvList_of_child {p c : ItemId} (hn : (items.type p).isNode = true)
    (hc : c ∈ items.ch p) (hcV : items.type c ≠ .V) {x : Nat}
    (hx : (items.vs c).1 = some x ∨ (items.vs c).2 = some x) : x ∈ items.nvList g p := by
  have hxlt := h.endpoints.vs_lt c x hx
  rcases h.endpoints.child_vs_in_parent p c hc ((isNode_iff _).1 hn) hcV x hx with h1 | h1 | h1
  · rw [nvList_eq, h1]; simp
  · rw [nvList_eq, h1]; simp
  · rw [nvList_eq]
    refine List.mem_append_left _ (List.mem_append_right _ ?_)
    rw [List.mem_map]
    refine ⟨vertItem x, List.mem_filter.2 ⟨h1, by simp; riomega⟩, by show 1 + x - 1 = x; omega⟩

/-- Virtual edge `capCount + j` of a node and the cap of its `j`-th non-V child carry the same
original vertices. -/
theorem virt_glue (hr : items.RepOK g) {a : ItemId} (ha : a < items.size)
    (hn : (items.type a).isNode = true) {pos : Nat → Nat} (hl : RelabelLayout g items t idx a pos)
    {j : Nat} (hj : j < ((items.ordered g a (t.nvRange (idx a)).1 pos).filter (· ≥ 1 + g.nv)).length) :
    ∃ p q, t.neOrig ((t.neRange (idx a)).1 + items.capCount a + j) = some p ∧
      t.neOrig (t.neRange (idx ((items.ordered g a (t.nvRange (idx a)).1 pos).filter
        (· ≥ 1 + g.nv))[j])).1 = some q ∧ Items.PairEq p q := by
  set F := (items.ordered g a (t.nvRange (idx a)).1 pos).filter (· ≥ 1 + g.nv) with hF
  have hperm : F.Perm ((items.ch a).filter (· ≥ 1 + g.nv)) := ordered_filter_perm a _ pos _
  have hFlen : F.length = (items.virtualEdges a).length := by
    rw [hperm.length_eq, h.virt_length ha, List.countP_eq_length_filter]
  have hcF : F[j] ∈ F := List.getElem_mem hj
  have hcm : F[j] ∈ items.ch a := List.mem_of_mem_filter (hperm.mem_iff.1 hcF)
  have hcge : F[j] ≥ 1 + g.nv := by simpa using List.of_mem_filter (hperm.mem_iff.1 hcF)
  have hcV : items.type F[j] ≠ .V := fun hV =>
    absurd ((h.type_V_iff ha hcm).1 hV) (not_lt.2 hcge)
  obtain ⟨x, y, hx, hy, hcq, hO⟩ := h.child_cap_orig hr ha hn hcm hcV
  have hxm : x ∈ items.nvList g a := h.mem_nvList_of_child hn hcm hcV (Or.inl hx)
  have hcs := h.ch_lt ha hcm
  have hym : y ∈ items.nvList g a := by
    by_cases hO' : items.type F[j] = .O
    · obtain ⟨v, hv, -⟩ := h.nvList_O hcs hO'
      rw [hv] at hy; simp at hy; subst hy; exact hxm
    · exact h.mem_nvList_of_child hn hcm hcV (Or.inr (hO hO'))
  have hnE := h.nEdges_node ha hn
  have hnv := (h.node a ha).nv_range
  have hne := (h.node a ha).ne_range
  have hk : items.capCount a + j < items.nEdges g a := by omega
  have hlay := hl.edge_nvs _ hk
  rw [nodeLayout, hnv, hne, Items.edgeChildren, ← hF] at hlay
  have hnotO : items.type a ≠ .Q → items.type F[j] ≠ .O :=
    fun hQ hO' => hQ (hr.o_parent a _ hcm hO').1
  set nvSt := (t.nvRange (idx a)).1 with hnvSt
  set neSt := (t.neRange (idx a)).1 with hneSt
  set len := (items.nvList g a).length with hlen0
  cases hty : items.type a
  · rw [hty] at hn; cases hn
  · rw [hty] at hn; cases hn
  · -- Q: a block root with one non-V child
    obtain ⟨e, he, rfl⟩ := h.type_Q_eq ha hty
    obtain ⟨u, v, hu, -, hcase, -, hne1⟩ := h.q_pair hr he
    rcases hcase with ⟨h0, -, -⟩ | ⟨hne0, hvs, c, hc⟩
    · rw [h0] at hcm; cases hcm
    have hcap0 : items.capCount (edgeItem g e) = 0 := by
      simp [Items.capCount, Items.hasCap, hty, hne0]
    have hj0 : j = 0 := by omega
    subst hj0
    simp only [hcap0, Nat.add_zero] at hlay hk ⊢
    have hpos := nvList_pos (g := g) hu
    have hnvs := h.edge0_nvs hr ha hn (by omega)
    have hp := h.neOrig_eq ha (k := 0) (a := 0) (b := len - 1) (by omega) (by omega) (by omega)
      (by rw [Nat.add_zero]; exact hnvs)
    rw [Nat.add_zero] at hp
    refine ⟨_, (x, y), hp, hcq, ?_⟩
    rcases hc with ⟨hch, hnl⟩ | ⟨hch, hnl, -⟩
    · have hvlt : v < g.nv := h.nvList_lt ha (by rw [hnl]; simp)
      have hO' : items.type F[0] ≠ .O := fun hO' => by
        have := (hr.o_parent _ _ hcm hO').2
        rw [hch] at this; cases this
      have hy' := hO hO'
      have hnd := h.nv_nodup hcs
      have hxy : x ≠ y := by
        rw [nvList_eq, hx, hy'] at hnd
        exact ne_of_nodup_head_last hnd
      have hpv : ((items.nvList g (edgeItem g e))[0]'(by omega),
          (items.nvList g (edgeItem g e))[len - 1]'(by omega)) = (u, v) := by
        simp [hlen0, hnl]
      rw [hpv]
      rw [hnl] at hxm hym
      simp only [List.mem_cons, List.mem_nil_iff, or_false] at hxm hym
      rcases hxm with rfl | rfl <;> rcases hym with rfl | rfl
      · exact absurd rfl hxy
      · exact Or.inl rfl
      · exact Or.inr rfl
      · exact absurd rfl hxy
    · have hpv : ((items.nvList g (edgeItem g e))[0]'(by omega),
          (items.nvList g (edgeItem g e))[len - 1]'(by omega)) = (u, u) := by
        simp [hlen0, hnl]
      rw [hpv]
      rw [hnl] at hxm hym
      simp only [List.mem_cons, List.mem_nil_iff, or_false] at hxm hym
      subst hxm; subst hym
      exact Or.inl rfl
  · rw [h.shapes.i_o_leaf a ha (Or.inl hty)] at hcm; cases hcm
  · rw [h.shapes.i_o_leaf a ha (Or.inr hty)] at hcm; cases hcm
  · -- S
    rw [hty] at hlay
    obtain ⟨u, v, xs, hvs, hnl, hxs, hvirt⟩ := h.nvList_S hr ha hty
    have hcap1 : items.capCount a = 1 := by
      simp [Items.capCount, Items.hasCap, hty, NodeType.isNode]
    have hFeq : F = (items.ch a).filter fun c => items.type c ≠ .V := by
      rw [hF, Items.ordered, ite_eq_left (by rw [hty]; decide), h.filter_ge_eq ha]
    have hvirt' : items.virtualEdges a =
        F.map fun c => ((items.vs c).1.getD 0, (items.vs c).2.getD 0) := by
      rw [Items.virtualEdges, hFeq]
    have hlen : len = xs.length + 2 := by rw [hlen0, hnl]; simp
    have hvl : (items.virtualEdges a).length = xs.length + 1 := by
      rw [hvirt, List.length_zip]; simp
    rw [hcap1] at hlay hk ⊢
    rw [hvl, hcap1] at hnE
    have hs := LayoutFacts.s_edge (idx a) nvSt (nvSt + len) neSt (neSt + items.nEdges g a)
      (F.map fun c => (pos ((items.vs c).1.getD 0), pos ((items.vs c).2.getD 0)))
      (by omega) (by omega) (1 + j) (by omega)
    rw [ite_eq_right (by omega)] at hs
    rw [hs] at hlay
    have hp := h.neOrig_eq ha (k := 1 + j) (a := j) (b := j + 1) hk (by omega) (by omega)
      (by rw [hlay]; exact Prod.ext (by dsimp only; omega) (by dsimp only; omega))
    rw [← Nat.add_assoc] at hp
    refine ⟨_, (x, y), hp, hcq, Or.inl ?_⟩
    have hy' := hO (hnotO (by rw [hty]; decide))
    have hvj : (items.virtualEdges a)[j]'(by omega) = (x, y) := by
      rw [List.getElem_of_eq hvirt', List.getElem_map, hx, hy']; rfl
    rw [List.getElem_of_eq hvirt, List.getElem_zip] at hvj
    rw [← hvj]
    have hj1 : j < (u :: xs).length := by simp; omega
    have hj2 : j < (xs ++ [v]).length := by simp; omega
    have e1 : (items.nvList g a)[j]? = some ((u :: xs)[j]'hj1) := by
      rw [hnl, List.getElem?_append_left hj1, List.getElem?_eq_getElem]
    have e2 : (items.nvList g a)[j + 1]? = some ((xs ++ [v])[j]'hj2) := by
      rw [hnl, List.cons_append, List.getElem?_cons_succ, List.getElem?_eq_getElem]
    rw [(List.getElem_eq_iff _).2 e1, (List.getElem_eq_iff _).2 e2]
  · -- P
    rw [hty] at hlay
    obtain ⟨u, v, hvs, hnl⟩ := h.nvList_P ha hty
    have hcap1 : items.capCount a = 1 := by
      simp [Items.capCount, Items.hasCap, hty, NodeType.isNode]
    have hlen : len = 2 := by rw [hlen0, hnl]; rfl
    rw [hcap1] at hlay hk ⊢
    rw [hcap1] at hnE
    rw [hlen] at hlay
    have hs := LayoutFacts.p_edge (idx a) nvSt neSt (neSt + items.nEdges g a)
      (F.map fun c => (pos ((items.vs c).1.getD 0), pos ((items.vs c).2.getD 0))) (1 + j) (by omega)
    rw [hs] at hlay
    have hp := h.neOrig_eq ha (k := 1 + j) (a := 0) (b := 1) hk (by omega) (by omega) (by rw [hlay]; rfl)
    rw [← Nat.add_assoc] at hp
    refine ⟨_, (x, y), hp, hcq, ?_⟩
    have hy' := hO (hnotO (by rw [hty]; decide))
    have hvm : (x, y) ∈ items.virtualEdges a := by
      rw [Items.virtualEdges, List.mem_map]
      exact ⟨F[j], List.mem_filter.2 ⟨hcm, by simpa using hcV⟩, by rw [hx, hy']; rfl⟩
    obtain ⟨u', v', hvs', hpq⟩ := (h.shapes.p_shape a ha hty).2.2 _ hvm
    rw [hvs] at hvs'
    cases hvs'
    have hpv : ((items.nvList g a)[0]'(by omega), (items.nvList g a)[1]'(by omega)) = (u, v) := by
      simp [hnl]
    rw [hpv]
    exact PairEq.symm' hpq
  · -- R
    rw [hty] at hlay
    obtain ⟨u, v, hvs, hlen4⟩ := h.nvList_R ha hty
    have hcap1 : items.capCount a = 1 := by
      simp [Items.capCount, Items.hasCap, hty, NodeType.isNode]
    rw [hcap1] at hlay hk ⊢
    rw [hcap1] at hnE
    have hy' := hO (hnotO (by rw [hty]; decide))
    have hpos := hl.pos_ok hty
    obtain ⟨hxle, hxp⟩ := hpos x hxm
    obtain ⟨hyle, hyp⟩ := hpos y hym
    have hxa := (List.getElem?_eq_some_iff.1 hxp).1
    have hyb := (List.getElem?_eq_some_iff.1 hyp).1
    set E := F.map fun c => (pos ((items.vs c).1.getD 0), pos ((items.vs c).2.getD 0)) with hE
    have hElen : E.length = F.length := List.length_map _
    have hs := LayoutFacts.r_edge (idx a) nvSt (nvSt + len) neSt (neSt + items.nEdges g a) E
      (by omega) (by omega) (1 + j) (by omega)
    rw [ite_eq_right (by omega), show 1 + j - 1 = j by omega, getElem!_pos E j (by omega)] at hs
    have hEj : E[j]'(by omega) = (pos x, pos y) := by
      simp only [hE, List.getElem_map, hx, hy', Option.getD_some]
    rw [hEj] at hs
    rw [hs] at hlay
    have hp := h.neOrig_eq ha (k := 1 + j) (a := pos x - nvSt) (b := pos y - nvSt) hk hxa hyb
      (by rw [hlay]; exact Prod.ext (by dsimp only; omega) (by dsimp only; omega))
    rw [← Nat.add_assoc] at hp
    refine ⟨_, (x, y), hp, hcq, Or.inl ?_⟩
    rw [(List.getElem_eq_iff _).2 hxp, (List.getElem_eq_iff _).2 hyp]

theorem twin_glue (hr : items.RepOK g) : ∀ ne ne', t.twin ne = some ne' →
    ∀ p q, t.neOrig ne = some p → t.neOrig ne' = some q → Items.PairEq p q := by
  intro ne ne' ht p q hp hq
  obtain ⟨a, ha, h1, h2⟩ := h.ne_mem_range (twin_some_lt ht)
  obtain ⟨pos, hl⟩ := (h.node a ha).layout
  have hner := (h.node a ha).ne_range
  rw [hner] at h2
  by_cases hn : (items.type a).isNode = true
  swap
  · exfalso
    have := nEdges_not_node (g := g) (i := a) (by simpa using hn)
    omega
  have hnE := h.nEdges_node ha hn
  have hFlen : ((items.ordered g a (t.nvRange (idx a)).1 pos).filter (· ≥ 1 + g.nv)).length =
      (items.virtualEdges a).length := by
    rw [(ordered_filter_perm a _ pos _).length_eq, h.virt_length ha, List.countP_eq_length_filter]
  by_cases hk : (t.neRange (idx a)).1 + items.capCount a ≤ ne
  · -- `ne` is a virtual edge of `a`
    have hj : ne - ((t.neRange (idx a)).1 + items.capCount a) <
        ((items.ordered g a (t.nvRange (idx a)).1 pos).filter (· ≥ 1 + g.nv)).length := by omega
    obtain ⟨p', q', hp', hq', hpq⟩ := h.virt_glue hr ha hn hl hj
    have hne : (t.neRange (idx a)).1 + items.capCount a +
        (ne - ((t.neRange (idx a)).1 + items.capCount a)) = ne := by omega
    rw [hne] at hp'
    have htw := (hl.twin hn _ hj).1
    rw [hne, ht] at htw
    cases htw
    rw [hp] at hp'; rw [hq] at hq'
    cases hp'; cases hq'
    exact hpq
  · -- `ne` is the cap of `a`: twinned with a virtual edge of its parent
    have hc : items.hasCap a = true := by
      by_contra hc
      have : items.capCount a = 0 := by simp [Items.capCount, hc]
      omega
    have hne : ne = (t.neRange (idx a)).1 := by
      have : items.capCount a = 1 := by simp [Items.capCount, hc]
      omega
    subst hne
    have ha0 : a ≠ rootItem := by
      rintro rfl; rw [h.tree.root] at hn; cases hn
    obtain ⟨p0, hp0, -⟩ := h.tree.unique_parent a (Nat.pos_of_ne_zero ha0) ha
    have hp0s := parent_lt hp0
    by_cases hn0 : (items.type p0).isNode = true
    · obtain ⟨pos0, hl0⟩ := (h.node p0 hp0s).layout
      have hperm := ordered_perm (items := items) (g := g) p0 (t.nvRange (idx p0)).1 pos0
      have haV : items.type a ≠ .V := by intro hV; rw [hV] at hn; cases hn
      have hage : ¬ a < 1 + g.nv := fun hlt => haV ((h.type_V_iff hp0s hp0).2 hlt)
      have hmem : a ∈ (items.ordered g p0 (t.nvRange (idx p0)).1 pos0).filter (· ≥ 1 + g.nv) :=
        List.mem_filter.2 ⟨hperm.mem_iff.2 hp0, by simpa using hage⟩
      obtain ⟨j, hj, hja⟩ := List.getElem_of_mem hmem
      obtain ⟨p', q', hp', hq', hpq⟩ := h.virt_glue hr hp0s hn0 hl0 hj
      have htw := (hl0.twin hn0 j hj).2
      rw [hja] at htw hq'
      rw [ht] at htw
      cases htw
      rw [hq] at hp'; rw [hp] at hq'
      cases hp'; cases hq'
      exact PairEq.symm' hpq
    · exfalso
      have := (h.node p0 hp0s).child_cap_twin_none (by simpa using hn0) a hp0 hc
      rw [this] at ht; cases ht

/-! ### R skeleton transport -/

omit h in
theorem skeleton_mem (n : Nat) (p : Nat × Nat) :
    p ∈ t.skeleton n ↔ ∃ k, k < (t.neRange n).2 - (t.neRange n).1 ∧
      (t.nodeEdges[(t.neRange n).1 + k]!).nvs = p := by
  simp [SpqrTree.skeleton, SpqrTree.nodeEdgesOf, Array.getElem!_eq_getD]

theorem idxOf_eq_of_pos {a : ItemId} (ha : a < items.size)
    (hR : items.type a = .R) {pos : Nat → Nat} (hl : RelabelLayout g items t idx a pos)
    {x : Nat} (hx : x ∈ items.nvList g a) :
    (items.nvList g a).idxOf x = pos x - (t.nvRange (idx a)).1 := by
  obtain ⟨-, hp⟩ := hl.pos_ok hR x hx
  obtain ⟨hlt, hget⟩ := List.getElem?_eq_some_iff.1 hp
  have := (h.nv_nodup ha).idxOf_getElem _ hlt
  rwa [hget] at this

theorem r_virt_nvs (hr : items.RepOK g) {a : ItemId} (ha : a < items.size)
    (hR : items.type a = .R) {pos : Nat → Nat} (hl : RelabelLayout g items t idx a pos)
    {j : Nat}
    (hj : j < ((items.ordered g a (t.nvRange (idx a)).1 pos).filter (· ≥ 1 + g.nv)).length) :
    ∃ x y, items.vs ((items.ordered g a (t.nvRange (idx a)).1 pos).filter (· ≥ 1 + g.nv))[j] =
        (some x, some y) ∧
      x ∈ items.nvList g a ∧ y ∈ items.nvList g a ∧
      (t.nodeEdges[(t.neRange (idx a)).1 + (1 + j)]!).nvs = (pos x, pos y) := by
  have hn : (items.type a).isNode = true := by rw [hR]; rfl
  obtain ⟨u, v, -, hlen4⟩ := h.nvList_R ha hR
  set F := (items.ordered g a (t.nvRange (idx a)).1 pos).filter (· ≥ 1 + g.nv) with hF
  have hcF : F[j] ∈ F := List.getElem_mem hj
  have hcm : F[j] ∈ items.ch a := (ordered_perm a _ pos).mem_iff.1 (List.mem_filter.1 hcF).1
  have hcge : F[j] ≥ 1 + g.nv := by simpa using (List.mem_filter.1 hcF).2
  have hcV : items.type F[j] ≠ .V := fun hV =>
    absurd ((h.type_V_iff ha hcm).1 hV) (not_lt.2 hcge)
  have hcO : items.type F[j] ≠ .O := fun hO => by
    have := (hr.o_parent a _ hcm hO).1; rw [hR] at this; cases this
  obtain ⟨x, y, hx, -, -, hO⟩ := h.child_cap_orig hr ha hn hcm hcV
  have hy' := hO hcO
  have hxm := h.mem_nvList_of_child hn hcm hcV (Or.inl hx)
  have hym := h.mem_nvList_of_child hn hcm hcV (Or.inr hy')
  refine ⟨x, y, Prod.ext hx hy', hxm, hym, ?_⟩
  have hFlen : F.length = (items.virtualEdges a).length := by
    rw [(ordered_filter_perm a _ pos _).length_eq, h.virt_length ha, List.countP_eq_length_filter]
  have hnE := h.nEdges_node ha hn
  have hcap1 : items.capCount a = 1 := by
    simp [Items.capCount, Items.hasCap, hR, NodeType.isNode]
  have hk : 1 + j < items.nEdges g a := by omega
  have hlay := hl.edge_nvs _ hk
  rw [nodeLayout, (h.node a ha).nv_range, (h.node a ha).ne_range, Items.edgeChildren, ← hF, hR]
    at hlay
  set nvSt := (t.nvRange (idx a)).1 with hnvSt
  set neSt := (t.neRange (idx a)).1 with hneSt
  set len := (items.nvList g a).length with hlen0
  set E := F.map fun c => (pos ((items.vs c).1.getD 0), pos ((items.vs c).2.getD 0)) with hE
  have hElen : E.length = F.length := List.length_map _
  have hs := LayoutFacts.r_edge (idx a) nvSt (nvSt + len) neSt (neSt + items.nEdges g a) E
    (by omega) (by omega) (1 + j) (by omega)
  rw [ite_eq_right (by omega), show 1 + j - 1 = j by omega, getElem!_pos E j (by omega)] at hs
  have hEj : E[j]'(by omega) = (pos x, pos y) := by
    simp only [hE, List.getElem_map, hx, hy', Option.getD_some]
  rw [hEj] at hs
  rw [hs] at hlay
  exact hlay

theorem r_three_connected (hr : items.RepOK g) (hR : items.RThreeConnected g) :
    ∀ i, i < t.size → t.type i = .R →
      SpqrTree.ThreeConnected (t.nVerts i)
        ((t.skeleton i).map fun p => (p.1 - (t.nvRange i).1, p.2 - (t.nvRange i).1)) := by
  intro n hn hty
  obtain ⟨a, ha, rfl⟩ := h.idx_surj hn
  rw [h.type_eq ha] at hty
  have hnode : (items.type a).isNode = true := by rw [hty]; rfl
  obtain ⟨pos, hl⟩ := (h.node a ha).layout
  obtain ⟨u, v, hvs, hlen4⟩ := h.nvList_R ha hty
  have hnd := h.nv_nodup ha
  have hnv := (h.node a ha).nv_range
  have hne := (h.node a ha).ne_range
  have hnE := h.nEdges_node ha hnode
  have hcap1 : items.capCount a = 1 := by
    simp [Items.capCount, Items.hasCap, hty, NodeType.isNode]
  set F := (items.ordered g a (t.nvRange (idx a)).1 pos).filter (· ≥ 1 + g.nv) with hF
  have hFlen : F.length = (items.virtualEdges a).length := by
    rw [(ordered_filter_perm a _ pos _).length_eq, h.virt_length ha, List.countP_eq_length_filter]
  have hperm : F.Perm ((items.ch a).filter fun c => items.type c ≠ .V) := by
    rw [← h.filter_ge_eq ha]; exact ordered_filter_perm a _ pos _
  have hu0 : (items.nvList g a)[0]'(by omega) = u :=
    (List.getElem_eq_iff _).2 (vs_head (g := g) (show (items.vs a).1 = some u by rw [hvs]))
  have hvl : (items.nvList g a)[(items.nvList g a).length - 1]'(by omega) = v :=
    (List.getElem_eq_iff _).2 (vs_last (g := g) (show (items.vs a).2 = some v by rw [hvs]))
  have hidx_u : (items.nvList g a).idxOf u = 0 := by
    have := hnd.idxOf_getElem 0 (by omega); rwa [hu0] at this
  have hidx_v : (items.nvList g a).idxOf v = (items.nvList g a).length - 1 := by
    have := hnd.idxOf_getElem ((items.nvList g a).length - 1) (by omega); rwa [hvl] at this
  have hcap := h.edge0_nvs hr ha hnode (by omega)
  have hnV : t.nVerts (idx a) = (items.nvList g a).length := by
    rw [SpqrTree.nVerts, hnv]; omega
  rw [hnV, SpqrTree.ThreeConnected_congr (es' := items.rSkeleton g a)]
  · exact hR a ha hty
  intro p
  rw [List.mem_map, Items.rSkeleton, List.mem_map]
  constructor
  · rintro ⟨q, hq, rfl⟩
    rw [skeleton_mem, hne, Nat.add_sub_cancel_left] at hq
    obtain ⟨k, hk, hq⟩ := hq
    rcases Nat.eq_zero_or_pos k with rfl | hkpos
    · rw [Nat.add_zero, hcap] at hq
      refine ⟨((items.vs a).1.getD 0, (items.vs a).2.getD 0), List.mem_cons_self .., ?_⟩
      rw [← hq, hvs]
      simp only [Option.getD_some, hidx_u, hidx_v]
      exact Prod.ext (by dsimp only; omega) (by dsimp only; omega)
    · obtain ⟨j, rfl⟩ : ∃ j, k = 1 + j := ⟨k - 1, by omega⟩
      obtain ⟨x, y, hvsj, hxm, hym, hnvs⟩ := h.r_virt_nvs hr ha hty hl (j := j) (by show j < F.length; omega)
      rw [hnvs] at hq
      refine ⟨((items.vs F[j]).1.getD 0, (items.vs F[j]).2.getD 0), List.mem_cons_of_mem _ ?_, ?_⟩
      · rw [Items.virtualEdges, List.mem_map]
        exact ⟨F[j], hperm.subset (List.getElem_mem _), rfl⟩
      · rw [← hq, hvsj]
        simp only [Option.getD_some]
        rw [h.idxOf_eq_of_pos ha hty hl hxm, h.idxOf_eq_of_pos ha hty hl hym]
  · rintro ⟨q, hq, rfl⟩
    rcases List.mem_cons.1 hq with rfl | hq
    · refine ⟨_, (skeleton_mem _ _).2 ⟨0, by rw [hne, Nat.add_sub_cancel_left]; omega, rfl⟩, ?_⟩
      rw [Nat.add_zero, hcap, hvs]
      simp only [Option.getD_some, hidx_u, hidx_v]
      exact Prod.ext (by dsimp only; omega) (by dsimp only; omega)
    · rw [Items.virtualEdges, List.mem_map] at hq
      obtain ⟨c, hc, rfl⟩ := hq
      have hcF : c ∈ F := hperm.mem_iff.2 hc
      obtain ⟨j, hj, rfl⟩ := List.mem_iff_getElem.1 hcF
      obtain ⟨x, y, hvsj, hxm, hym, hnvs⟩ := h.r_virt_nvs hr ha hty hl hj
      refine ⟨_, (skeleton_mem _ _).2 ⟨1 + j, by rw [hne, Nat.add_sub_cancel_left]; omega, rfl⟩, ?_⟩
      rw [hnvs, hvsj]
      simp only [Option.getD_some]
      rw [h.idxOf_eq_of_pos ha hty hl hxm, h.idxOf_eq_of_pos ha hty hl hym]

/-! ### Assembly -/

/-- `Represents` from the per-node interface, `Items.RepOK`, and item-level R 3-connectivity. -/
theorem represents (hr : items.RepOK g) (hR : items.RThreeConnected g) :
    t.Represents g where
  nv := h.ridx.nv
  ne := h.ridx.ne
  q_endpoints := h.q_endpoints hr
  twin_glue := h.twin_glue hr
  nv_orig_inj := h.nv_orig_inj fun _ hi => h.nv_nodup hi
  separation := h.separation hr
  interior := h.interior
  r_three_connected := h.r_three_connected hr hR
  canonical := h.canonical

/-- `Represents` with the R clause supplied at the output level (the form `Correctness.lean`
assumes). -/
theorem represents_of_r (hr : items.RepOK g)
    (hR : ∀ i, i < t.size → t.type i = .R →
      SpqrTree.ThreeConnected (t.nVerts i)
        ((t.skeleton i).map fun p => (p.1 - (t.nvRange i).1, p.2 - (t.nvRange i).1))) :
    t.Represents g where
  nv := h.ridx.nv
  ne := h.ridx.ne
  q_endpoints := h.q_endpoints hr
  twin_glue := h.twin_glue hr
  nv_orig_inj := h.nv_orig_inj fun _ hi => h.nv_nodup hi
  separation := h.separation hr
  interior := h.interior
  r_three_connected := hR
  canonical := h.canonical

end RelabelOK

theorem relabelOK_of_wf (g : Graph) (items : Items) (h : items.WF g) :
    ∃ idx, RelabelOK g items (relabelTree g items) idx := by
  obtain ⟨idx, -, hidx, hnode⟩ := relabel_node_spec g items h
  exact ⟨idx, ⟨h, hidx, hnode⟩⟩

/-- Phase-D main result: the relabelled tree represents `g`, given the item contracts. -/
theorem relabelTree_represents' (g : Graph) (items : Items) (h : items.WF g)
    (hr : items.RepOK g) (hR : items.RThreeConnected g) :
    (relabelTree g items).Represents g := by
  obtain ⟨idx, hok⟩ := relabelOK_of_wf g items h
  exact hok.represents hr hR

/-- Variant with the output-level R clause, matching `Correctness.relabelTree_represents` up to the
extra `Items.RepOK` hypothesis. -/
theorem relabelTree_represents_of_r (g : Graph) (items : Items) (h : items.WF g)
    (hr : items.RepOK g)
    (hR : ∀ i, i < (relabelTree g items).size → (relabelTree g items).type i = .R →
      SpqrTree.ThreeConnected ((relabelTree g items).nVerts i)
        (((relabelTree g items).skeleton i).map fun p =>
          (p.1 - ((relabelTree g items).nvRange i).1, p.2 - ((relabelTree g items).nvRange i).1))) :
    (relabelTree g items).Represents g := by
  obtain ⟨idx, hok⟩ := relabelOK_of_wf g items h
  exact hok.represents_of_r hr hR

end Spqr
