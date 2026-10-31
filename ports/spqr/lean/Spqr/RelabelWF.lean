import Spqr.RelabelAdj

/-!
# `SpqrTree.WF` of the relabel output

Assembles `SpqrTree.WF` (`Spec.lean`) for `relabelTree g items` from the per-node interface
`relabel_node_spec` (taken as a hypothesis, as in `RelabelOwn.lean`/`RelabelAdj.lean`):
`Preorder` from `RelabelLayout.child_idx`/`subtree_end` (the children of a node are numbered
`idx i + 1`, then each after its predecessor's subtree, and the node's subtree ends with its last
child's), `Shape` from `LayoutShape.shape_*`, and the adjacency clauses from `Layout.Local`.
-/

namespace Spqr

namespace Items

variable {g : Graph} {items : Items}

theorem Tree.card_desc_child_lt (ht : items.Tree g) {a c : ItemId} (ha : a < items.size)
    (hc : c ∈ items.ch a) : (items.desc c).card < (items.desc a).card := by
  have h1 := Finset.card_le_card (ht.desc_subset hc)
  rw [Finset.card_erase_of_mem (mem_desc_self ha)] at h1
  have h2 := one_le_card_desc (items := items) ha
  omega

end Items

namespace RelabelAll

variable {g : Graph} {items : Items} {t : SpqrTree} {idx : ItemId → Nat}
variable (H : RelabelAll g items t idx)
include H

/-! ### Preorder: the children chain -/

/-- Where the `k`-th child of `i` is numbered: `idx i + 1` for the first, the end of the previous
child's subtree otherwise. -/
def chainEnd (t : SpqrTree) (idx : ItemId → Nat) (i : ItemId) (L : List ItemId) (k : Nat) : Nat :=
  if k = 0 then idx i + 1 else t.subtreeEnd[idx L[k - 1]!]!

omit H in
theorem chainEnd_zero (i : ItemId) (L : List ItemId) : chainEnd t idx i L 0 = idx i + 1 := rfl

omit H in
theorem chainEnd_succ (i : ItemId) (L : List ItemId) (k : Nat) :
    chainEnd t idx i L (k + 1) = t.subtreeEnd[idx L[k]!]! := rfl

omit H in
theorem child_idx_eq {i : ItemId} {pos : Nat → Nat}
    (hl : RelabelLayout g items t idx i pos) {k : Nat}
    (hk : k < (items.ordered g i (t.nvRange (idx i)).1 pos).length) :
    idx (items.ordered g i (t.nvRange (idx i)).1 pos)[k]! =
      chainEnd t idx i (items.ordered g i (t.nvRange (idx i)).1 pos) k := by
  have := hl.child_idx k hk
  rw [getElem!_pos _ k hk]
  unfold chainEnd
  split_ifs with h0
  · subst h0; simpa using this
  · simp only [h0, ↓reduceIte] at this
    rw [this, getElem!_pos _ (k - 1) (by omega)]

omit H in
theorem subtreeEnd_eq_chain {i : ItemId} {pos : Nat → Nat}
    (hl : RelabelLayout g items t idx i pos) :
    t.subtreeEnd[idx i]! = chainEnd t idx i (items.ordered g i (t.nvRange (idx i)).1 pos)
      (items.ordered g i (t.nvRange (idx i)).1 pos).length := by
  have := hl.subtree_end
  set L := items.ordered g i (t.nvRange (idx i)).1 pos with hL
  rcases hLnil : L with _ | ⟨a, L'⟩
  · rw [hLnil] at this; simpa [chainEnd] using this
  · rw [hLnil] at this
    have hne : ¬ (a :: L').length = 0 := by simp
    rw [List.getLast?_eq_getElem?, getElem?_pos (a :: L') _ (by simp)] at this
    unfold chainEnd; simp only [hne, ↓reduceIte]; rw [getElem!_pos (a :: L') _ (by simp)]
    exact this

omit H in
/-- The chain of a node whose children's subtrees are nonempty: ends increase, every index between
the node and the end of the chain lies in a child's subtree. -/
theorem chain {i : ItemId} {pos : Nat → Nat}
    (hl : RelabelLayout g items t idx i pos)
    (hch : ∀ c ∈ items.ch i, idx c < t.subtreeEnd[idx c]!) :
    ∀ k, k ≤ (items.ordered g i (t.nvRange (idx i)).1 pos).length →
      (∀ k', k' < k → chainEnd t idx i (items.ordered g i (t.nvRange (idx i)).1 pos) k' <
        chainEnd t idx i (items.ordered g i (t.nvRange (idx i)).1 pos) (k' + 1)) ∧
      ∀ j, idx i < j → j < chainEnd t idx i (items.ordered g i (t.nvRange (idx i)).1 pos) k →
        ∃ k', k' < k ∧ idx (items.ordered g i (t.nvRange (idx i)).1 pos)[k']! ≤ j ∧
          j < t.subtreeEnd[idx (items.ordered g i (t.nvRange (idx i)).1 pos)[k']!]! := by
  set L := items.ordered g i (t.nvRange (idx i)).1 pos with hL
  intro k
  induction k with
  | zero => intro _; exact ⟨fun _ h => absurd h (Nat.not_lt_zero _), fun j h1 h2 => by
      rw [chainEnd_zero] at h2; omega⟩
  | succ k ih =>
    intro hk
    obtain ⟨ih1, ih2⟩ := ih (by omega)
    have hmem : L[k]! ∈ items.ch i := by
      rw [getElem!_pos L k hk]; exact (Items.ordered_perm (i := i) (pos := pos)).mem_iff.1 (List.getElem_mem hk)
    have hc := hch _ hmem
    have hidx := child_idx_eq hl hk
    rw [← hL] at hidx
    have hstep : chainEnd t idx i L k < chainEnd t idx i L (k + 1) := by
      rw [chainEnd_succ, ← hidx]; exact hc
    refine ⟨fun k' hk' => ?_, fun j h1 h2 => ?_⟩
    · rcases Nat.lt_succ_iff_lt_or_eq.1 hk' with h | rfl
      · exact ih1 k' h
      · exact hstep
    · by_cases hj : j < chainEnd t idx i L k
      · obtain ⟨k', h1', h2', h3'⟩ := ih2 j h1 hj
        exact ⟨k', by omega, h2', h3'⟩
      · refine ⟨k, Nat.lt_succ_self k, ?_, ?_⟩
        · rw [hidx]; omega
        · rw [chainEnd_succ] at h2; exact h2

omit H in
theorem chainEnd_mono {i : ItemId} {pos : Nat → Nat}
    (hl : RelabelLayout g items t idx i pos)
    (hch : ∀ c ∈ items.ch i, idx c < t.subtreeEnd[idx c]!) {k k' : Nat} (hkk : k ≤ k')
    (hk' : k' ≤ (items.ordered g i (t.nvRange (idx i)).1 pos).length) :
    chainEnd t idx i (items.ordered g i (t.nvRange (idx i)).1 pos) k ≤
      chainEnd t idx i (items.ordered g i (t.nvRange (idx i)).1 pos) k' := by
  induction k', hkk using Nat.le_induction with
  | base => exact le_rfl
  | succ m hm ih =>
    exact le_trans (ih (by omega)) ((chain hl hch (m + 1) hk').1 m (Nat.lt_succ_self m)).le

/-- Strong induction over the item tree: subtrees are nonempty and every index strictly inside a
node's subtree range lies in a child's. -/
theorem subtree_props : ∀ i, i < items.size →
    idx i < t.subtreeEnd[idx i]! ∧
    ∀ j, idx i < j → j < t.subtreeEnd[idx i]! →
      ∃ c ∈ items.ch i, idx c ≤ j ∧ j < t.subtreeEnd[idx c]! := by
  suffices ∀ n i, i < items.size → (items.desc i).card = n →
      idx i < t.subtreeEnd[idx i]! ∧
      ∀ j, idx i < j → j < t.subtreeEnd[idx i]! →
        ∃ c ∈ items.ch i, idx c ≤ j ∧ j < t.subtreeEnd[idx c]! from
    fun i hi => this _ i hi rfl
  intro n
  induction n using Nat.strong_induction_on with
  | _ n ih =>
  intro i hi hcard
  have hch : ∀ c ∈ items.ch i, idx c < t.subtreeEnd[idx c]! := by
    intro c hc
    exact (ih _ (hcard ▸ H.tree.card_desc_child_lt hi hc) c (H.tree.child_lt hc) rfl).1
  obtain ⟨pos, hl⟩ := (H.node i hi).layout
  have hend := subtreeEnd_eq_chain hl
  obtain ⟨-, h2⟩ := chain hl hch _ le_rfl
  refine ⟨?_, fun j h1 hj => ?_⟩
  · rw [hend]
    exact lt_of_lt_of_le (Nat.lt_succ_self _)
      (chainEnd_mono hl hch (Nat.zero_le _) le_rfl)
  · rw [hend] at hj
    obtain ⟨k', hk', h3, h4⟩ := h2 j h1 hj
    refine ⟨_, ?_, h3, h4⟩
    rw [getElem!_pos (items.ordered g i (t.nvRange (idx i)).1 pos) _ hk']
    exact (Items.ordered_perm (i := i) (pos := pos)).mem_iff.1 (List.getElem_mem hk')

theorem child_subtree_pos (i : ItemId) :
    ∀ c ∈ items.ch i, idx c < t.subtreeEnd[idx c]! :=
  fun c hc => (H.subtree_props c (H.tree.child_lt hc)).1

/-- A child is numbered after its parent, inside the parent's subtree range. -/
theorem child_idx_bounds {i : ItemId} (hi : i < items.size) {c : ItemId} (hc : c ∈ items.ch i) :
    idx i < idx c ∧ idx c < t.subtreeEnd[idx i]! := by
  obtain ⟨pos, hl⟩ := (H.node i hi).layout
  have hch := H.child_subtree_pos i
  set L := items.ordered g i (t.nvRange (idx i)).1 pos with hL
  have hcL : c ∈ L := (Items.ordered_perm (i := i) (pos := pos)).mem_iff.2 hc
  obtain ⟨k, hk, hkc⟩ := List.getElem_of_mem hcL
  rw [← hkc, ← getElem!_pos L k hk, child_idx_eq hl hk]
  have h0 := chainEnd_mono hl hch (Nat.zero_le k) hk.le
  have h1 := (chain hl hch (k + 1) hk).1 k (Nat.lt_succ_self k)
  have h2 := chainEnd_mono hl hch hk le_rfl
  have hend := subtreeEnd_eq_chain hl
  simp only [← hL, chainEnd_zero, chainEnd_succ] at h0 h1 h2 hend ⊢
  rw [hend]
  omega

theorem children_pairwise {i : ItemId} {pos : Nat → Nat}
    (hl : RelabelLayout g items t idx i pos) :
    ((items.ordered g i (t.nvRange (idx i)).1 pos).map idx).Pairwise (· < ·) := by
  have hch := H.child_subtree_pos i
  set L := items.ordered g i (t.nvRange (idx i)).1 pos with hL
  rw [List.pairwise_iff_getElem]
  intro a b ha hb hab
  simp only [List.length_map] at ha hb
  simp only [List.getElem_map]
  rw [← getElem!_pos L a ha, ← getElem!_pos L b hb, child_idx_eq hl ha, child_idx_eq hl hb]
  have h1 := (chain hl hch (a + 1) ha).1 a (Nat.lt_succ_self a)
  have h2 := chainEnd_mono hl hch hab hb.le
  simp only [← hL, Nat.succ_eq_add_one] at h1 h2 ⊢
  omega

theorem children_sum {i : ItemId} {pos : Nat → Nat}
    (hl : RelabelLayout g items t idx i pos) :
    ∀ m, m ≤ (items.ordered g i (t.nvRange (idx i)).1 pos).length →
      ((List.range m).map fun k =>
          t.subtreeEnd[idx (items.ordered g i (t.nvRange (idx i)).1 pos)[k]!]! -
            idx (items.ordered g i (t.nvRange (idx i)).1 pos)[k]!).sum + (idx i + 1) =
        chainEnd t idx i (items.ordered g i (t.nvRange (idx i)).1 pos) m := by
  have hch := H.child_subtree_pos i
  set L := items.ordered g i (t.nvRange (idx i)).1 pos with hL
  intro m
  induction m with
  | zero => intro _; simp [chainEnd_zero]
  | succ m ih =>
    intro hm
    rw [List.range_succ, List.map_append, List.sum_append, List.map_singleton, List.sum_singleton]
    have := ih (by omega)
    have h1 := (chain hl hch (m + 1) hm).1 m (Nat.lt_succ_self m)
    have h2 := child_idx_eq hl (show m < L.length by omega)
    simp only [← hL, chainEnd_succ] at h1 h2 this ⊢
    omega

/-- Every node with a parent is the image of a child item, its parent the image of the item's
parent. -/
theorem parent_some {j p : Nat} (h : t.parent j = some p) :
    ∃ q c, q < items.size ∧ c ∈ items.ch q ∧ idx c = j ∧ idx q = p := by
  have hj : j < t.size := by
    by_contra hj
    have : t.par[j]? = none := Array.getElem?_eq_none_iff.2 (by rw [H.gl.sizes.par]; omega)
    simp [SpqrTree.parent, this] at h
  obtain ⟨c, hc, rfl⟩ := H.idx_surj hj
  have hc0 : c ≠ 0 := by
    rintro rfl
    rw [show idx 0 = 0 from H.gl.root, H.gl.root_par] at h; cases h
  obtain ⟨q, hq⟩ := H.tree.parent_exists hc hc0
  have hq' := H.tree.parent_lt hq
  have := (H.node q hq').child_par c hq
  rw [this] at h
  exact ⟨q, c, hq', hq, rfl, Option.some_inj.1 h⟩

theorem parent_eq_iff {i : ItemId} (hi : i < items.size) (j : Nat) :
    t.parent j = some (idx i) ↔ ∃ c ∈ items.ch i, idx c = j := by
  constructor
  · intro h
    obtain ⟨q, c, hq, hc, rfl, hqi⟩ := H.parent_some h
    rw [H.idx_inj hq hi hqi] at hc
    exact ⟨c, hc, rfl⟩
  · rintro ⟨c, hc, rfl⟩
    exact (H.node i hi).child_par c hc

theorem root_lt : rootItem < items.size := by
  have := H.tree.size; simp only [rootItem]; iomega

/-- `Preorder.only_root_F`: non-root nodes are not F (`Items.Tree.type_F_iff` via `idx`). -/
theorem only_root_F : ∀ n, 0 < n → n < t.size → t.type n ≠ .F := by
  intro n hn0 hn
  obtain ⟨c, hc, rfl⟩ := H.idx_surj hn
  have hc0 : c ≠ 0 := by
    rintro rfl; rw [show idx 0 = 0 from H.gl.root] at hn0; exact absurd hn0 (lt_irrefl 0)
  rw [H.type hc]
  exact fun h => hc0 ((H.tree.type_F_iff hc).1 h)

theorem preorder : t.Preorder := by
  have hroot := H.gl.root
  refine ⟨?_, H.gl.root_par, ?_, ?_, ?_, H.ch_mono, ?_, H.only_root_F⟩
  · have := H.type H.root_lt
    rw [hroot] at this; rw [this]; exact H.tree.root
  · intro n hn0 hn
    obtain ⟨c, hc, rfl⟩ := H.idx_surj hn
    have hc0 : c ≠ 0 := by rintro rfl; rw [show idx 0 = 0 from hroot] at hn0; exact absurd hn0 (lt_irrefl 0)
    obtain ⟨q, hq⟩ := H.tree.parent_exists hc hc0
    have hq' := H.tree.parent_lt hq
    exact ⟨idx q, (H.node q hq').child_par c hq, (H.child_idx_bounds hq' hq).1⟩
  · intro j p h
    obtain ⟨q, c, hq, hc, rfl, rfl⟩ := H.parent_some h
    have := H.child_idx_bounds hq hc
    refine ⟨this.1.le, ?_⟩
    have e : t.subtreeEnd[idx q]?.getD 0 = t.subtreeEnd[idx q]! :=
      (Array.getElem!_eq_getD_getElem? _ _).symm
    rw [e]; exact this.2
  · intro n hn
    obtain ⟨i, hi, rfl⟩ := H.idx_surj hn
    obtain ⟨pos, hl⟩ := (H.node i hi).layout
    rw [hl.children]
    have hpw := H.children_pairwise hl
    refine List.Perm.eq_of_pairwise (le := (· < ·)) (fun a b _ _ h1 h2 => by omega) hpw
      (List.pairwise_lt_range.filter _) ?_
    refine LayoutShape.perm_of_mem_iff (hpw.imp Nat.ne_of_lt)
      (List.nodup_range.filter _) fun j => ?_
    rw [List.mem_filter, List.mem_range, List.mem_map, decide_eq_true_eq, H.parent_eq_iff hi]
    constructor
    · rintro ⟨c, hc, rfl⟩
      have hc' := (Items.ordered_perm (i := i) (pos := pos)).mem_iff.1 hc
      exact ⟨H.idx_lt (H.tree.child_lt hc'), c, hc', rfl⟩
    · rintro ⟨-, c, hc, rfl⟩
      exact ⟨c, (Items.ordered_perm (i := i) (pos := pos)).mem_iff.2 hc, rfl⟩
  · intro n hn
    obtain ⟨i, hi, rfl⟩ := H.idx_surj hn
    obtain ⟨pos, hl⟩ := (H.node i hi).layout
    rw [hl.children, subtreeEnd_eq_chain hl, ← H.children_sum hl _ le_rfl]
    set L := items.ordered g i (t.nvRange (idx i)).1 pos with hL
    have : (L.map idx).map (fun c => t.subtreeEnd[c]! - c) =
        (List.range L.length).map fun k => t.subtreeEnd[idx L[k]!]! - idx L[k]! := by
      apply List.ext_getElem (by simp)
      intro k h1 h2
      simp only [List.getElem_map, List.getElem_range]
      rw [getElem!_pos L k (by simpa using h2)]
    rw [this]; omega

/-! ### Shape -/

theorem skeleton_eq {i : ItemId} (hi : i < items.size) {pos : Nat → Nat}
    (hl : RelabelLayout g items t idx i pos)
    (hloc : (nodeLayout g items t idx i pos).Local (idx i) (t.nvRange (idx i)).1 (t.nvRange (idx i)).2
      (t.neRange (idx i)).1 (t.neRange (idx i)).2) :
    (nodeLayout g items t idx i pos).skeleton = t.skeleton (idx i) := by
  have hne := (H.node i hi).ne_range
  have hsz := hloc.edges_size
  unfold Layout.skeleton SpqrTree.skeleton SpqrTree.nodeEdgesOf
  apply LayoutShape.toList_map_eq
  · rw [List.length_map, List.length_map, List.length_range]; omega
  · intro k hk
    have hk' : k < items.nEdges g i := by omega
    rw [List.map_map, getElem!_pos (List.map _ (List.range _)) k
      (by rw [List.length_map, List.length_range]; omega)]
    simp only [List.getElem_map, List.getElem_range, Function.comp]
    rw [← Array.getElem!_eq_getD_getElem?, hl.edge_nvs k hk']

open LayoutShape in
theorem shape_local (hor : items.ROriented g) {i : ItemId} (hi : i < items.size) {pos : Nat → Nat}
    (hl : RelabelLayout g items t idx i pos) :
    (nodeLayout g items t idx i pos).Shape (items.type i) (t.nvRange (idx i)).1
      (t.nvRange (idx i)).2 := by
  have ht := H.tree
  have hnv := (H.node i hi).nv_range
  have hne := (H.node i hi).ne_range
  have hlen := Items.nvList_length (g := g) (items := items) i
  unfold nodeLayout
  by_cases hn : (items.type i).isNode = true
  · have hy := H.wf.layout_hyps hi hn
    have h1 := H.wf.one_le_nvList hi hn
    have hpos := H.wf.nEdges_pos hi hn
    cases hty : items.type i
    all_goals simp only [hty] at hn
    all_goals simp [NodeType.isNode] at hn
    · -- Q
      rcases H.wf.nvList_Q hi hty with hv1 | hv2
      · exact shape_loop _ _ _ _ _ _ _ (Or.inl rfl) (by omega) (by have := hy.1 hv1; omega)
      · exact shape_QI _ _ _ _ _ _ _ (Or.inl rfl) (by omega) (by have := hy.2.1 (Or.inl hty); omega)
    · -- I
      obtain ⟨u, v, huv⟩ := H.wf.vs_two hi (Or.inl hty)
      have h0 := H.wf.shapes.i_o_leaf i hi (Or.inl hty)
      rw [huv, h0] at hlen; simp at hlen
      exact shape_QI _ _ _ _ _ _ _ (Or.inr rfl) (by omega) (by have := hy.2.1 (Or.inr hty); omega)
    · -- O
      have hv1 := hy.2.2.2.2.2 hty
      exact shape_loop _ _ _ _ _ _ _ (Or.inr rfl) (by omega) (by have := hy.1 hv1; omega)
    · -- S
      obtain ⟨u, v, xs, huv, hxs, hk, -⟩ := H.wf.shapes.s_shape i hi hty
      have hV : ((items.ch i).filter (· < 1 + g.nv)).length = xs.length := by
        rw [← ht.filter_V_eq, hxs, List.length_map]
      rw [huv] at hlen; simp at hlen
      exact shape_S _ _ _ _ _ _ (by omega) (by have := hy.2.2.1 hty; omega)
    · -- P
      have hcap := Items.hasCap_of_ne_Q (items := items) (i := i) (by rw [hty]; rfl)
        (by rw [hty]; decide)
      have h2 := hy.2.2.2.2.1 hcap (Or.inr (Or.inr hty))
      obtain ⟨u, v, huv⟩ := H.wf.vs_two hi (Or.inr (Or.inr (Or.inl hty)))
      rw [huv] at hlen; simp at hlen
      have hn0 : (items.type i).isNode = true := by rw [hty]; rfl
      have h3 := (H.wf.shapes.p_shape i hi hty).1
      have hcnt := ht.countP_nonV i
      have hnE := Items.nEdges_eq (g := g) hn0
      simp only [Items.capCount, hcap, ↓reduceIte] at hnE
      exact shape_P _ _ _ _ _ _ (by omega) (by omega)
    · -- R
      have h4 := hy.2.2.2.1 hty
      have hcap := Items.hasCap_of_ne_Q (items := items) (i := i) (by rw [hty]; rfl)
        (by rw [hty]; decide)
      have hec := Items.edgeChildren_length (g := g) (items := items) i (t.nvRange (idx i)).1 pos
      have hE := H.edgeChildren_bounds hor hi hty hl
      have hn0 : (items.type i).isNode = true := by rw [hty]; rfl
      have hcnt := ht.countP_nonV i
      have hnE := Items.nEdges_eq (g := g) hn0
      simp only [Items.capCount, hcap, ↓reduceIte] at hnE
      have hposok := hl.pos_ok hty
      have hnd := H.wf.nv_nodup hi
      obtain ⟨u, v, huv⟩ := H.wf.vs_two hi (Or.inr (Or.inr (Or.inr hty)))
      obtain ⟨-, h6, hnodup, -, hnpar, -⟩ := H.wf.shapes.r_shape i hi hty
      have hvE := H.wf.r_edges_in_nv hi hty
      set E := items.edgeChildren g pos (items.ordered g i (t.nvRange (idx i)).1 pos) with hEdef
      set nvSt := (t.nvRange (idx i)).1 with hnvSt
      have hEperm : E.Perm ((items.virtualEdges i).map fun q => (pos q.1, pos q.2)) := by
        rw [hEdef]
        unfold Items.edgeChildren Items.virtualEdges
        rw [List.map_map]
        refine (((Items.ordered_perm (i := i) (pos := pos)).filter _).map _).trans ?_
        rw [ht.filter_nonV_eq]
        exact List.Perm.refl _
      have hinj : ∀ a ∈ items.nvList g i, ∀ b ∈ items.nvList g i, pos a = pos b → a = b := by
        intro a ha b hb hab
        have h1 := (hposok a ha).2
        have h2 := (hposok b hb).2
        rw [hab, h2] at h1
        exact (Option.some_inj.1 h1).symm
      have hEnd : E.Nodup := by
        refine hEperm.nodup_iff.2 (List.Nodup.map_on ?_ (hnodup.of_map _))
        intro q hq q' hq' h
        obtain ⟨h1, h2⟩ := Prod.mk.inj h
        exact Prod.ext (hinj _ (hvE q hq).1 _ (hvE q' hq').1 h1) (hinj _ (hvE q hq).2 _ (hvE q' hq').2 h2)
      have hu0 : (items.nvList g i)[0]? = some u := by
        simp [Items.nvList, huv]
      have hvl : (items.nvList g i)[(items.nvList g i).length - 1]? = some v := by
        rw [← List.getLast?_eq_getElem?]
        simp only [Items.nvList, huv, Option.toList_some, List.singleton_append]
        rw [List.getLast?_concat]
      have hum : u ∈ items.nvList g i := List.mem_iff_getElem?.2 ⟨_, hu0⟩
      have hvm : v ∈ items.nvList g i := List.mem_iff_getElem?.2 ⟨_, hvl⟩
      have hpu : pos u = nvSt := by
        have := Items.PosOK.pos_sub hnd hposok hum
        rw [Items.idxOf_eq_of_getElem? hnd hu0] at this
        have := (hposok u hum).1; omega
      have hpv : pos v = (t.nvRange (idx i)).2 - 1 := by
        have := Items.PosOK.pos_sub hnd hposok hvm
        rw [Items.idxOf_eq_of_getElem? hnd hvl] at this
        have := (hposok v hvm).1; omega
      have hnotin : (nvSt, (t.nvRange (idx i)).2 - 1) ∉ E := by
        intro hmem
        rw [hEperm.mem_iff] at hmem
        obtain ⟨q, hq, hqe⟩ := List.mem_map.1 hmem
        obtain ⟨h1, h2⟩ := Prod.mk.inj hqe
        rw [← hpu] at h1; rw [← hpv] at h2
        have hq1 := hinj _ (hvE q hq).1 _ hum h1
        have hq2 := hinj _ (hvE q hq).2 _ hvm h2
        exact hnpar u v huv q hq (Or.inl (Prod.ext hq1 hq2))
      exact shape_R _ _ _ _ _ _ (by omega) hE (by omega) (by omega) (by omega)
        (List.nodup_cons.2 ⟨hnotin, hEnd⟩)
  · have hn' : (items.type i).isNode = false := Bool.eq_false_iff.2 hn
    have hne0 := Items.nEdges_eq_zero (g := g) hn'
    cases hty : items.type i
    all_goals simp only [hty] at hn'
    all_goals simp [NodeType.isNode] at hn'
    · exact shape_F _ _ _ _ _ _ (by omega)
    · have h0 := H.wf.nvList_V hi hty
      rw [h0] at hnv; simp at hnv
      exact shape_V _ _ _ _ _ _ (by omega) (by omega)

theorem shape (hor : items.ROriented g) {i : ItemId} (hi : i < items.size) :
    t.Shape (idx i) := by
  obtain ⟨pos, hl, hloc⟩ := H.layout_local hor hi
  have hS := H.shape_local hor hi hl
  have hsk := H.skeleton_eq hi hl hloc
  have hed : t.nEdges (idx i) = (nodeLayout g items t idx i pos).edges.size := by
    have hne := (H.node i hi).ne_range
    unfold SpqrTree.nEdges; rw [hloc.edges_size]
  rcases hr : t.nvRange (idx i) with ⟨s, e⟩
  simp only [SpqrTree.Shape, Layout.Shape, hr] at hS ⊢
  rw [H.type hi, ← hsk, hed]
  exact hS

/-! ### Adjacency -/

omit H in
theorem mono_le (f : Nat → Nat) (lo n : Nat)
    (hm : ∀ r, lo ≤ r → r < lo + n → f r ≤ f (r + 1)) :
    ∀ d a, lo ≤ a → a + d ≤ lo + n → f a ≤ f (a + d) := by
  intro d
  induction d with
  | zero => intros; simp
  | succ d ih =>
    intro a ha hd
    exact (ih a ha (by omega)).trans (by rw [← Nat.add_assoc]; exact hm _ (by omega) (by omega))

omit H in
theorem find_row (f : Nat → Nat) (lo : Nat) :
    ∀ n, (∀ r, lo ≤ r → r < lo + n → f r ≤ f (r + 1)) →
      ∀ x, f lo ≤ x → x < f (lo + n) →
        ∃ r, lo ≤ r ∧ r < lo + n ∧ f r ≤ x ∧ x < f (r + 1) := by
  intro n
  induction n with
  | zero => intro _ x h1 h2; rw [Nat.add_zero] at h2; omega
  | succ n ih =>
    intro hm x h1 h2
    by_cases hx : x < f (lo + n)
    · obtain ⟨r, hr1, hr2, hr3, hr4⟩ := ih (fun r h1 h2 => hm r h1 (by omega)) x h1 hx
      exact ⟨r, hr1, by omega, hr3, hr4⟩
    · exact ⟨lo + n, by omega, by omega, by omega, by rw [Nat.add_assoc]; exact h2⟩

open LayoutR (row rowBound) in
/-- Global adjacency bounds of node `i` are its local row bounds. -/
theorem global_bound (hor : items.ROriented g) {i : ItemId} (hi : i < items.size) {pos : Nat → Nat}
    (hl : RelabelLayout g items t idx i pos) :
    ∀ r, 2 * (t.nvRange (idx i)).1 ≤ r → r ≤ 2 * (t.nvRange (idx i)).2 →
      t.adjBounds[r]! =
        rowBound (t.nvRange (idx i)).1 (t.neRange (idx i)).1 (nodeLayout g items t idx i pos) r := by
  intro r h1 h2
  have hnv := (H.node i hi).nv_range
  unfold rowBound
  split
  · rename_i h; rw [h]; exact (H.adj_spec hor).1 (idx i) (H.idx_lt hi)
  · have hj := hl.adj_bounds (r - 2 * (t.nvRange (idx i)).1) (by omega) (by omega)
    rw [show 2 * (t.nvRange (idx i)).1 + (r - 2 * (t.nvRange (idx i)).1) = r by omega] at hj
    exact hj

theorem adj_bounds_mono (hor : items.ROriented g) :
    ∀ r, r + 1 < t.adjBounds.size → t.adjBounds[r]! ≤ t.adjBounds[r + 1]! := by
  intro r hr
  have hsz := H.gl.sizes.adjBounds
  obtain ⟨i, hi, h1, h2⟩ := H.nv_locate (nv := r / 2) (by omega)
  obtain ⟨pos, hl, hloc⟩ := H.layout_local hor hi
  rw [H.global_bound hor hi hl r (by omega) (by omega),
    H.global_bound hor hi hl (r + 1) (by omega) (by omega)]
  exact hloc.bound_mono r (by omega) (by omega)

omit H in
open LayoutR (row rowBound) in
theorem rowBound_start {i : ItemId} {pos : Nat → Nat} :
    rowBound (t.nvRange (idx i)).1 (t.neRange (idx i)).1 (nodeLayout g items t idx i pos)
      (2 * (t.nvRange (idx i)).1) = 2 * (t.neRange (idx i)).1 := by
  unfold rowBound; simp

open LayoutR (row rowBound) in
/-- Local row bounds of node `i` lie in `[2 neSt, 2 neEn]`. -/
theorem rowBound_range {i : ItemId} (hi : i < items.size) {pos : Nat → Nat}
    (hloc : (nodeLayout g items t idx i pos).Local (idx i) (t.nvRange (idx i)).1 (t.nvRange (idx i)).2
      (t.neRange (idx i)).1 (t.neRange (idx i)).2) :
    ∀ r, 2 * (t.nvRange (idx i)).1 ≤ r → r ≤ 2 * (t.nvRange (idx i)).2 →
      2 * (t.neRange (idx i)).1 ≤
        rowBound (t.nvRange (idx i)).1 (t.neRange (idx i)).1 (nodeLayout g items t idx i pos) r ∧
      rowBound (t.nvRange (idx i)).1 (t.neRange (idx i)).1 (nodeLayout g items t idx i pos) r ≤
        2 * (t.neRange (idx i)).2 := by
  intro r h1 h2
  have hnv := (H.node i hi).nv_range
  have hm := mono_le _ (2 * (t.nvRange (idx i)).1) (2 * (items.nvList g i).length)
    (fun r h1 h2 => hloc.bound_mono r h1 (by omega))
  constructor
  · have := hm (r - 2 * (t.nvRange (idx i)).1) _ le_rfl (by omega)
    rwa [show 2 * (t.nvRange (idx i)).1 + (r - 2 * (t.nvRange (idx i)).1) = r by omega,
      rowBound_start] at this
  · have := hm (2 * (t.nvRange (idx i)).2 - r) r h1 (by omega)
    rwa [show r + (2 * (t.nvRange (idx i)).2 - r) = 2 * (t.nvRange (idx i)).2 by omega,
      hloc.bound_last] at this

open LayoutR (row rowBound) in
/-- Global row `r` of node `i` is its local row. -/
theorem global_row (hor : items.ROriented g) {i : ItemId} (hi : i < items.size) {pos : Nat → Nat}
    (hl : RelabelLayout g items t idx i pos)
    (hloc : (nodeLayout g items t idx i pos).Local (idx i) (t.nvRange (idx i)).1 (t.nvRange (idx i)).2
      (t.neRange (idx i)).1 (t.neRange (idx i)).2) :
    ∀ r, 2 * (t.nvRange (idx i)).1 ≤ r → r < 2 * (t.nvRange (idx i)).2 →
      (List.range (t.adjBounds[r + 1]! - t.adjBounds[r]!)).map (fun k => t.adjDat[t.adjBounds[r]! + k]!) =
        row (t.nvRange (idx i)).1 (t.neRange (idx i)).1 (nodeLayout g items t idx i pos) r := by
  intro r h1 h2
  have hne := (H.node i hi).ne_range
  unfold row
  rw [H.global_bound hor hi hl r h1 h2.le, H.global_bound hor hi hl (r + 1) (by omega) (by omega)]
  have hb := H.rowBound_range hi hloc r h1 h2.le
  have hb' := H.rowBound_range hi hloc (r + 1) (by omega) (by omega)
  apply List.map_congr_left
  intro k hk
  rw [List.mem_range] at hk
  have := hl.adj_dat (rowBound (t.nvRange (idx i)).1 (t.neRange (idx i)).1
    (nodeLayout g items t idx i pos) r + k - 2 * (t.neRange (idx i)).1) (by omega)
  rwa [show 2 * (t.neRange (idx i)).1 + (rowBound (t.nvRange (idx i)).1 (t.neRange (idx i)).1
    (nodeLayout g items t idx i pos) r + k - 2 * (t.neRange (idx i)).1) =
    rowBound (t.nvRange (idx i)).1 (t.neRange (idx i)).1 (nodeLayout g items t idx i pos) r + k
    by omega] at this

open LayoutR (row rowBound) in
theorem adj_dest (hor : items.ROriented g) :
    ∀ k, k < t.adjDat.size → ∀ nvs, t.nvsOf t.adjDat[k]!.ne = some nvs →
      t.adjDat[k]!.destNv = nvs.1 ∨ t.adjDat[k]!.destNv = nvs.2 := by
  intro k hk nvs hnvs
  have hds := H.gl.adj_dat_size
  obtain ⟨i, hi, h1, h2⟩ := H.ne_locate (ne := k / 2) (by omega)
  obtain ⟨pos, hl, hloc⟩ := H.layout_local hor hi
  have hne := (H.node i hi).ne_range
  have hnv := (H.node i hi).nv_range
  have hk' := hl.adj_dat (k - 2 * (t.neRange (idx i)).1) (by omega)
  rw [show 2 * (t.neRange (idx i)).1 + (k - 2 * (t.neRange (idx i)).1) = k by omega] at hk'
  obtain ⟨r, hr1, hr2, hr3, hr4⟩ := find_row
    (rowBound (t.nvRange (idx i)).1 (t.neRange (idx i)).1 (nodeLayout g items t idx i pos))
    (2 * (t.nvRange (idx i)).1) (2 * (items.nvList g i).length)
    (fun r h1 h2 => hloc.bound_mono r h1 (by omega)) k (by rw [rowBound_start]; omega)
    (by rw [show 2 * (t.nvRange (idx i)).1 + 2 * (items.nvList g i).length = 2 * (t.nvRange (idx i)).2
      by omega, hloc.bound_last]; omega)
  have hmem : (nodeLayout g items t idx i pos).adjDat[k - 2 * (t.neRange (idx i)).1]! ∈
      row (t.nvRange (idx i)).1 (t.neRange (idx i)).1 (nodeLayout g items t idx i pos) r := by
    unfold row
    refine List.mem_map.2 ⟨k - rowBound (t.nvRange (idx i)).1 (t.neRange (idx i)).1
      (nodeLayout g items t idx i pos) r, List.mem_range.2 (by omega), ?_⟩
    congr 1; omega
  obtain ⟨ha1, ha2, ha3⟩ := hloc.adj_dest r hr1 (by omega) _ hmem
  rw [hk'] at hnvs ⊢
  have hnvs' : t.nodeEdges[(nodeLayout g items t idx i pos).adjDat[k - 2 * (t.neRange (idx i)).1]!.ne]!.nvs
      = nvs := by
    unfold SpqrTree.nvsOf at hnvs
    obtain ⟨e, he, he2⟩ := Option.map_eq_some_iff.1 hnvs
    rw [Array.getElem!_eq_getD_getElem?, he, Option.getD_some]; exact he2
  have hev := hl.edge_nvs
    ((nodeLayout g items t idx i pos).adjDat[k - 2 * (t.neRange (idx i)).1]!.ne - (t.neRange (idx i)).1)
    (by omega)
  rw [show (t.neRange (idx i)).1 +
    ((nodeLayout g items t idx i pos).adjDat[k - 2 * (t.neRange (idx i)).1]!.ne - (t.neRange (idx i)).1) =
    (nodeLayout g items t idx i pos).adjDat[k - 2 * (t.neRange (idx i)).1]!.ne by omega] at hev
  rw [← hnvs', hev]
  exact ha3

omit H in
theorem filter_range_restrict (p : Nat → Bool) (N a b : Nat) (hab : a ≤ b) (hbN : b ≤ N)
    (hout : ∀ x, x < N → ¬ (a ≤ x ∧ x < b) → p x = false) :
    (List.range N).filter p = ((List.range (b - a)).filter fun k => p (a + k)).map (a + ·) := by
  have h0 : (List.range a).filter p = [] :=
    List.filter_eq_nil_iff.2 fun x hx => by
      rw [List.mem_range] at hx; simp [hout x (by omega) (by omega)]
  have h2 : ((List.range (N - b)).map ((a + (b - a)) + ·)).filter p = [] :=
    List.filter_eq_nil_iff.2 fun x hx => by
      obtain ⟨k, hk, rfl⟩ := List.mem_map.1 hx
      rw [List.mem_range] at hk
      simp only [Bool.not_eq_true]
      exact hout _ (by omega) (by omega)
  conv_lhs => rw [show N = a + (b - a) + (N - b) by omega]
  rw [List.range_add, List.range_add, List.filter_append, List.filter_append, List.filter_map, h0, h2,
    List.nil_append, List.append_nil]
  rfl

/-- The global node-edges of node `i` selected by an endpoint predicate are its local ones. -/
theorem row_filter {i : ItemId} (hi : i < items.size) {pos : Nat → Nat}
    (hl : RelabelLayout g items t idx i pos)
    (hloc : (nodeLayout g items t idx i pos).Local (idx i) (t.nvRange (idx i)).1 (t.nvRange (idx i)).2
      (t.neRange (idx i)).1 (t.neRange (idx i)).2)
    (sel : Nat × Nat → Nat) (nv : Nat)
    (hout : ∀ ne, ne < t.nodeEdges.size →
      ¬ ((t.neRange (idx i)).1 ≤ ne ∧ ne < (t.neRange (idx i)).2) → sel t.nodeEdges[ne]!.nvs ≠ nv) :
    (List.range t.nodeEdges.size).filter (fun ne => decide (sel t.nodeEdges[ne]!.nvs = nv)) =
      ((List.range (nodeLayout g items t idx i pos).edges.size).filter fun k =>
        decide (sel (nodeLayout g items t idx i pos).edges[k]!.nvs = nv)).map
        (· + (t.neRange (idx i)).1) := by
  have hne := (H.node i hi).ne_range
  have hle := H.neEn_le_size (H.idx_lt hi)
  have hc : ∀ k ∈ List.range ((t.neRange (idx i)).2 - (t.neRange (idx i)).1),
      decide (sel t.nodeEdges[(t.neRange (idx i)).1 + k]!.nvs = nv) =
        decide (sel (nodeLayout g items t idx i pos).edges[k]!.nvs = nv) := by
    intro k hk; rw [List.mem_range] at hk; rw [hl.edge_nvs k (by omega)]
  rw [filter_range_restrict (fun ne => decide (sel t.nodeEdges[ne]!.nvs = nv)) t.nodeEdges.size
    (t.neRange (idx i)).1 (t.neRange (idx i)).2 (by omega) hle
    (fun x hx hn => by simpa using hout x hx hn)]
  rw [hloc.edges_size, List.filter_congr hc]
  exact List.map_congr_left fun k _ => Nat.add_comm _ _

open LayoutR (row rowBound) in
/-- Node-edges of a node other than `i` have no endpoint among `i`'s node-vertices. -/
theorem foreign_ne (hor : items.ROriented g) {i : ItemId} (hi : i < items.size) {nv : Nat}
    (h1 : (t.nvRange (idx i)).1 ≤ nv) (h2 : nv < (t.nvRange (idx i)).2) :
    ∀ ne, ne < t.nodeEdges.size →
      ¬ ((t.neRange (idx i)).1 ≤ ne ∧ ne < (t.neRange (idx i)).2) →
      t.nodeEdges[ne]!.nvs.1 ≠ nv ∧ t.nodeEdges[ne]!.nvs.2 ≠ nv := by
  intro ne hne' hnot
  obtain ⟨j, hj, hj1, hj2⟩ := H.ne_locate hne'
  obtain ⟨pos', hl', hloc'⟩ := H.layout_local hor hj
  have hij : idx j ≠ idx i := by
    intro he; rw [he] at hj1 hj2; exact hnot ⟨hj1, hj2⟩
  have hnej := (H.node j hj).ne_range
  have hev := hl'.edge_nvs (ne - (t.neRange (idx j)).1) (by omega)
  rw [show (t.neRange (idx j)).1 + (ne - (t.neRange (idx j)).1) = ne by omega] at hev
  have hb := hloc'.ne_nvs (ne - (t.neRange (idx j)).1) (by rw [hloc'.edges_size]; omega)
  rw [← hev] at hb
  constructor
  · intro h; rw [h] at hb
    exact H.nv_disjoint (H.idx_lt hj) (H.idx_lt hi) hij ⟨hb.1, by omega⟩ ⟨h1, h2⟩
  · intro h; rw [h] at hb
    exact H.nv_disjoint (H.idx_lt hj) (H.idx_lt hi) hij ⟨by omega, hb.2.2⟩ ⟨h1, h2⟩

open LayoutR (row rowBound) in
/-- `WF.adj_incident`: rows `2 nv` and `2 nv + 1` list the node-edges ending, respectively
starting, at `nv`. -/
theorem adj_incident' (hor : items.ROriented g) : ∀ nv, nv < t.nodeVerts.size →
    ((List.range (t.adjBounds[2 * nv + 2]! - t.adjBounds[2 * nv]!)).map fun k =>
        (t.adjDat[t.adjBounds[2 * nv]! + k]!).ne).Perm
      (((List.range t.nodeEdges.size).filter fun ne => t.nodeEdges[ne]!.nvs.2 = nv) ++
       ((List.range t.nodeEdges.size).filter fun ne => t.nodeEdges[ne]!.nvs.1 = nv)) := by
  intro nv hnv
  obtain ⟨i, hi, h1, h2⟩ := H.nv_locate hnv
  obtain ⟨pos, hl, hloc⟩ := H.layout_local hor hi
  have hmono := H.adj_bounds_mono hor
  have hsz := H.gl.sizes.adjBounds
  have hb01 : t.adjBounds[2 * nv]! ≤ t.adjBounds[2 * nv + 1]! := hmono (2 * nv) (by omega)
  have hb12 : t.adjBounds[2 * nv + 1]! ≤ t.adjBounds[2 * nv + 2]! := hmono (2 * nv + 1) (by omega)
  have hsplit : ((List.range (t.adjBounds[2 * nv + 2]! - t.adjBounds[2 * nv]!)).map fun k =>
        (t.adjDat[t.adjBounds[2 * nv]! + k]!).ne) =
      ((List.range (t.adjBounds[2 * nv + 1]! - t.adjBounds[2 * nv]!)).map fun k =>
        t.adjDat[t.adjBounds[2 * nv]! + k]!).map (·.ne) ++
      ((List.range (t.adjBounds[2 * nv + 1 + 1]! - t.adjBounds[2 * nv + 1]!)).map fun k =>
        t.adjDat[t.adjBounds[2 * nv + 1]! + k]!).map (·.ne) := by
    have e : t.adjBounds[2 * nv + 2]! = t.adjBounds[2 * nv + 1 + 1]! :=
      congrArg (fun m => t.adjBounds[m]!) (by omega : 2 * nv + 2 = 2 * nv + 1 + 1)
    rw [show t.adjBounds[2 * nv + 2]! - t.adjBounds[2 * nv]! =
      (t.adjBounds[2 * nv + 1]! - t.adjBounds[2 * nv]!) +
        (t.adjBounds[2 * nv + 1 + 1]! - t.adjBounds[2 * nv + 1]!) by omega,
      List.range_add, List.map_append, List.map_map, List.map_map, List.map_map]
    congr 1
    apply List.map_congr_left
    intro k hk
    simp only [Function.comp]
    congr 2; omega
  rw [hsplit, H.global_row hor hi hl hloc (2 * nv) (by omega) (by omega),
    H.global_row hor hi hl hloc (2 * nv + 1) (by omega) (by omega),
    H.row_filter hi hl hloc Prod.snd nv (fun ne h hn => (H.foreign_ne hor hi h1 h2 ne h hn).2),
    H.row_filter hi hl hloc Prod.fst nv (fun ne h hn => (H.foreign_ne hor hi h1 h2 ne h hn).1)]
  exact (hloc.adj_incident_lo nv h1 h2).append (hloc.adj_incident_hi nv h1 h2)

theorem wf_tree (hor : items.ROriented g) : t.WF :=
  ⟨H.gl.sizes, H.preorder, H.bijections, H.ownership hor, H.twins,
    fun n hn => by obtain ⟨i, hi, rfl⟩ := H.idx_surj hn; exact H.shape hor hi,
    H.adj_bounds_mono hor, (H.adj_spec hor).2, H.adj_incident' hor, H.adj_dest hor⟩

end RelabelAll

end Spqr
