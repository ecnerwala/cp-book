import Spqr.RelabelSpec
import Spqr.StSpec
import Spqr.ItemTree
import Spqr.Proofs.Dfs

/-!
# Relabel phase D: `Represents` transport from the per-node interface

Everything here is proved for an abstract output tree `t` satisfying the phase-3 interface
`RelabelIdx` + `RelabelNode` (`Spqr/RelabelSpec.lean`), packaged as `RelabelOK`; the only
admission reaching `relabelTree` is `relabel_node_spec`.
-/

namespace Spqr

/-- `omega` after unfolding the `ItemId` abbreviation and the item constructors. -/
macro "iomega" : tactic =>
  `(tactic| ((try simp only [ItemId, vertItem, edgeItem, rootItem] at *); omega))

/-- `a[i]! = a[i]?.getD 0` for `Nat` arrays (the two spellings used by `Relabel.lean`). -/
theorem Array.getElem!_nat (a : Array Nat) (i : Nat) : a[i]! = a[i]?.getD 0 := by
  simp only [getElem!_def]; cases a[i]? <;> rfl

theorem Array.getElem!_opt {α : Type} (a : Array (Option α)) (i : Nat) :
    a[i]! = a[i]?.getD none := by
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
  have := h.tree.size; iomega

theorem edgeItem_lt {e : Nat} (he : e < g.ne) : edgeItem g e < items.size := by
  have := h.tree.size; iomega

theorem type_vertItem {v : Nat} (hv : v < g.nv) : items.type (vertItem v) = .V :=
  h.tree.vert v hv

theorem type_edgeItem {e : Nat} (he : e < g.ne) : items.type (edgeItem g e) = .Q :=
  h.tree.edge e he

/-- A child below `1 + g.nv` is a vertex item. -/
theorem ch_lt_vert {i c : ItemId} (hc : c ∈ items.ch i) (hlt : c < 1 + g.nv) :
    c = vertItem (c - 1) ∧ items.type c = .V := by
  have h0 : c ≠ 0 := h.ch_ne_root hc
  have : c = vertItem (c - 1) := by iomega
  refine ⟨this, ?_⟩
  rw [this]; exact h.type_vertItem (by iomega)

theorem type_V_iff {i c : ItemId} (hi : i < items.size) (hc : c ∈ items.ch i) :
    items.type c = .V ↔ c < 1 + g.nv := by
  constructor
  · intro hV
    by_contra hge
    have hcs := h.ch_lt hi hc
    by_cases hq : c < 1 + g.nv + g.ne
    · have : c = edgeItem g (c - (1 + g.nv)) := by iomega
      rw [this] at hV
      rw [h.type_edgeItem (by iomega)] at hV; cases hV
    · have := h.tree.node c (by iomega) hcs
      simp [hV] at this
  · intro hlt; exact (h.ch_lt_vert hc hlt).2

omit h in
theorem ch_of_ge {q : ItemId} (hq : items.size ≤ q) : items.ch q = [] := by
  simp [Items.ch, Array.getElem?_eq_none_iff.2 hq]

omit h in
theorem parent_lt {p c : ItemId} (hc : c ∈ items.ch p) : p < items.size := by
  by_contra hge
  rw [ch_of_ge (by iomega)] at hc
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
    iomega
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

theorem nvOrig_get {i : ItemId} (hi : i < items.size) {k : Nat}
    (hk : k < (items.nvList g i).length) :
    t.nvOrig ((t.nvRange (idx i)).1 + k) = some (items.nvList g i)[k] := by
  obtain ⟨-, hnv⟩ := h.nodeVerts_get hi hk
  have hv := h.nvList_lt hi (List.getElem_mem hk)
  generalize (items.nvList g i)[k] = v at hnv hv ⊢
  have hvi := h.ridx.vert_index _ hv
  have horig := (h.node _ (h.vertItem_lt hv)).orig
  rw [Items.origOf, h.type_vertItem hv] at horig
  simp only [SpqrTree.nvOrig, hnv, hvi, Option.getD_some, ← Array.getElem!_opt]
  show t.origId[idx (vertItem v)]! = some v
  rw [horig]
  show some (1 + v - 1) = some v
  congr 1; omega

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
    · obtain ⟨q, hq, -⟩ := h.tree.unique_parent c (by iomega) hc
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
  · obtain ⟨q, hq, -⟩ := h.tree.unique_parent c (by iomega) hc
    have := (h.node q (parent_lt hq)).child_par c hq
    rw [this] at hp
    exact ⟨c, q, hc, parent_lt hq, rfl, Option.some_injective _ hp, hq⟩

theorem mem_children_iff {i : ItemId} (hi : i < items.size) (n : Nat) :
    n ∈ t.children (idx i) ↔ ∃ c ∈ items.ch i, idx c = n := by
  obtain ⟨pos, hl⟩ := (h.node i hi).layout
  rw [hl.children, List.mem_map]
  simp only [(ordered_perm i _ pos).mem_iff]

end RelabelOK

end Spqr
