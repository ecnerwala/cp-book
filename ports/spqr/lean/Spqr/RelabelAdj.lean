import Spqr.RelabelOwn

/-!
# `relabel_adj_spec` from the per-node interface

The CSR facts about the output `adjBounds`: every node's rows start at `2 neSt`, and the last
bound is `adjDat.size`. Both follow from `relabel_node_spec` (taken as a hypothesis, as in
`RelabelOwn.lean`): `RelabelLayout.adj_bounds` identifies node `i`'s bounds `2 nvSt + j`,
`j ∈ [1, 2 nVerts]`, with those of `layoutNode`, whose last bound is `2 neEn` by
`Layout.Local.bound_last` (`LayoutShape.lean`; the R instance needs the edge children oriented along
the node-vert order, i.e. `Items.ROriented`). Walking the nodes in index order, bound `2 nvSt` of
node `n + 1` is bound `2 nvEn` of node `n` (`nvBounds[n + 1]`), so an induction on `n` gives
`adjBounds[2 nvBounds[n]] = 2 neBounds[n]` for all `n ≤ size`; `n = size` is the last bound.
-/

namespace Spqr

namespace Items

variable {g : Graph} {items : Items}

/-- A V item has no node-verts: no endpoints, and its children are Q items. -/
theorem WF.nvList_V (hw : items.WF g) {i : ItemId} (hi : i < items.size) (hV : items.type i = .V) :
    items.nvList g i = [] := by
  have hvs := hw.endpoints.vs_shape i hi
  simp only [hV] at hvs
  have hch : (items.ch i).filter (· < 1 + g.nv) = [] := by
    rw [List.filter_eq_nil_iff]
    intro c hc
    have hr := (hw.tree.type_V_iff hi).1 hV
    have hvi : vertItem (i - 1) = i := by simp only [vertItem]; iomega
    have hQ : items.type c = .Q := hw.tree.v_children (i - 1) c (by iomega) (by rw [hvi]; exact hc)
    have := hw.tree.child_ge_of_ne_V hc (by rw [hQ]; decide)
    simp only [decide_eq_true_eq]; omega
  unfold Items.nvList
  rw [hvs, hch]; rfl

/-- Every node has a node-edge: its cap, or (block-root Q) its node child. -/
theorem WF.nEdges_pos (hw : items.WF g) {i : ItemId} (hi : i < items.size)
    (hn : (items.type i).isNode = true) : 1 ≤ items.nEdges g i := by
  rw [nEdges_eq hn]
  by_cases hcap : items.hasCap i = true
  · simp [Items.capCount, hcap]
  · have hQ : items.type i = .Q := by
      by_contra h; exact hcap (hasCap_of_ne_Q hn h)
    rcases hw.q_cases hi hQ with h0 | ⟨c, hc, h1 | ⟨v, hv, h1⟩⟩
    · exact absurd (by rw [hasCap_Q hQ, h0]; rfl) hcap
    · have hcge := hw.tree.child_ge_of_ne_V (p := i) (c := c) (by rw [h1]; simp) (by simp at hc; tauto)
      rw [h1]; simp [hcge]
    · have hcge := hw.tree.child_ge_of_ne_V (p := i) (c := c) (by rw [h1]; simp) (by simp at hc; tauto)
      have hv' : ¬ 1 + g.nv ≤ vertItem v := Nat.not_le.2 (Nat.add_lt_add_left hv 1)
      rw [h1]; simp [hcge, hv']

/-- A block-root Q records one endpoint. -/
theorem WF.q_vs_snd_none (hw : items.WF g) {i : ItemId} (hi : i < items.size)
    (hQ : items.type i = .Q) (hch : items.ch i ≠ []) : (items.vs i).2 = none := by
  have hr := (hw.tree.type_Q_iff hi).1 hQ
  have he : edgeItem g (i - 1 - g.nv) = i := by simp only [edgeItem]; iomega
  obtain ⟨u, hu⟩ := hw.vs_fst hi (by rw [hQ]; rfl)
  have := (hw.endpoints.q_vs (i - 1 - g.nv) (by iomega) u (by rw [he]; exact hu)).2.1
  rw [he] at this
  exact this.2 hch

/-- Node-vert count of a Q item: `1` (block root without a V child, or loop leaf) or `2`. -/
theorem WF.nvList_Q (hw : items.WF g) {i : ItemId} (hi : i < items.size)
    (hQ : items.type i = .Q) :
    (items.nvList g i).length = 1 ∨ (items.nvList g i).length = 2 := by
  have hlen := nvList_length (g := g) (items := items) i
  have hvs := hw.endpoints.vs_shape i hi
  rw [hQ] at hvs
  rcases hw.q_cases hi hQ with h0 | ⟨c, hc, h1 | ⟨v, hv, h1⟩⟩
  · rw [h0] at hlen
    rcases hvs with ⟨v, hv⟩ | ⟨u, v, hv⟩ <;> rw [hv] at hlen <;> simp at hlen <;> omega
  · have h2 := hw.q_vs_snd_none hi hQ (by rw [h1]; simp)
    obtain ⟨u, hu⟩ := hw.vs_fst hi (by rw [hQ]; rfl)
    have hcge := hw.tree.child_ge_of_ne_V (p := i) (c := c) (by rw [h1]; simp) (by simp at hc; tauto)
    have hc' : ¬ c < 1 + g.nv := Nat.not_lt.2 hcge
    rw [h1, hu, h2] at hlen
    simp [hc'] at hlen; omega
  · have h2 := hw.q_vs_snd_none hi hQ (by rw [h1]; simp)
    obtain ⟨u, hu⟩ := hw.vs_fst hi (by rw [hQ]; rfl)
    have hcge := hw.tree.child_ge_of_ne_V (p := i) (c := c) (by rw [h1]; simp) (by simp at hc; tauto)
    have hc' : ¬ c < 1 + g.nv := Nat.not_lt.2 hcge
    have hv' : vertItem v < 1 + g.nv := Nat.add_lt_add_left hv 1
    rw [h1, hu, h2] at hlen
    simp [hc', hv'] at hlen; omega

end Items

namespace RelabelAll

variable {g : Graph} {items : Items} {t : SpqrTree} {idx : ItemId → Nat}
variable (H : RelabelAll g items t idx)
include H

/-- The edge children of an R node lie in its node-vert range and are oriented along it. -/
theorem edgeChildren_bounds (hor : items.ROriented g) {i : ItemId} (hi : i < items.size)
    (hR : items.type i = .R) {pos : Nat → Nat} (hl : RelabelLayout g items t idx i pos) :
    ∀ q ∈ items.edgeChildren g pos (items.ordered g i (t.nvRange (idx i)).1 pos),
      (t.nvRange (idx i)).1 ≤ q.1 ∧ q.1 < q.2 ∧ q.2 < (t.nvRange (idx i)).2 := by
  intro q hq
  obtain ⟨k, hk⟩ := List.mem_iff_getElem?.1 hq
  have hklt := (List.getElem?_eq_some_iff.1 hk).1
  obtain ⟨u, v, hu, hv, huv, heq⟩ := H.edgeChild_at hi hR hklt
  rw [heq] at hk
  obtain rfl := Option.some.inj hk
  have hpos := hl.pos_ok hR
  have hu1 := (hpos u hu).1
  have hv1 := (hpos v hv).1
  have hu2 := Items.PosOK.pos_lt hpos hu
  have hv2 := Items.PosOK.pos_lt hpos hv
  obtain ⟨hnd, hord⟩ := hor i hi hR
  have h3 := hord _ huv
  rw [← Items.PosOK.pos_sub hnd hpos hu, ← Items.PosOK.pos_sub hnd hpos hv] at h3
  have hnv := (H.node i hi).nv_range
  exact ⟨by omega, by omega, by omega⟩

open LayoutShape in
/-- Every item's layout satisfies the local adjacency interface. -/
theorem layout_local (hor : items.ROriented g) {i : ItemId} (hi : i < items.size) :
    ∃ pos, RelabelLayout g items t idx i pos ∧
      (nodeLayout g items t idx i pos).Local (idx i) (t.nvRange (idx i)).1 (t.nvRange (idx i)).2
        (t.neRange (idx i)).1 (t.neRange (idx i)).2 := by
  obtain ⟨pos, hl⟩ := (H.node i hi).layout
  refine ⟨pos, hl, ?_⟩
  have ht := H.tree
  have hnv := (H.node i hi).nv_range
  have hne := (H.node i hi).ne_range
  have hlen := Items.nvList_length (g := g) (items := items) i
  unfold nodeLayout
  by_cases hn : (items.type i).isNode = true
  · have hy := H.wf.layout_hyps hi hn
    have h1 := H.wf.one_le_nvList hi hn
    have hpos := H.wf.nEdges_pos hi hn
    by_cases hv1 : (items.nvList g i).length = 1
    · exact local_loop _ _ _ _ _ _ _
        ⟨by intro h; simp [h, NodeType.isNode] at hn, by intro h; simp [h, NodeType.isNode] at hn⟩
        (by omega) (by have := hy.1 hv1; omega)
    cases hty : items.type i
    all_goals simp only [hty] at hn
    all_goals simp [NodeType.isNode] at hn
    · -- Q
      have h2 := (H.wf.nvList_Q hi hty).resolve_left hv1
      exact local_QI _ _ _ _ _ _ _ (Or.inl rfl) (by omega) (by have := hy.2.1 (Or.inl hty); omega)
    · -- I
      obtain ⟨u, v, huv⟩ := H.wf.vs_two hi (Or.inl hty)
      have h0 := H.wf.shapes.i_o_leaf i hi (Or.inl hty)
      rw [huv, h0] at hlen; simp at hlen
      exact local_QI _ _ _ _ _ _ _ (Or.inr rfl) (by omega) (by have := hy.2.1 (Or.inr hty); omega)
    · -- O
      exact absurd (hy.2.2.2.2.2 hty) hv1
    · -- S
      obtain ⟨u, v, xs, huv, hxs, hk, -⟩ := H.wf.shapes.s_shape i hi hty
      have hV : ((items.ch i).filter (· < 1 + g.nv)).length = xs.length := by
        rw [← ht.filter_V_eq, hxs, List.length_map]
      rw [huv] at hlen; simp at hlen
      exact local_S _ _ _ _ _ _ (by omega) (by have := hy.2.2.1 hty; omega)
    · -- P
      have hcap := Items.hasCap_of_ne_Q (items := items) (i := i) (by rw [hty]; rfl)
        (by rw [hty]; decide)
      have h2 := hy.2.2.2.2.1 hcap (Or.inr (Or.inr hty))
      obtain ⟨u, v, huv⟩ := H.wf.vs_two hi (Or.inr (Or.inr (Or.inl hty)))
      rw [huv] at hlen; simp at hlen
      exact local_P _ _ _ _ _ _ (by omega) (by omega)
    · -- R
      have h4 := hy.2.2.2.1 hty
      have hcap := Items.hasCap_of_ne_Q (items := items) (i := i) (by rw [hty]; rfl)
        (by rw [hty]; decide)
      have hec := Items.edgeChildren_length (g := g) (items := items) i (t.nvRange (idx i)).1 pos
      have hE := H.edgeChildren_bounds hor hi hty hl
      have hn0 : (items.type i).isNode = true := by rw [hty]; rfl
      exact local_R _ _ _ _ _ _ (by omega) hE
        (by rw [hec, hne, Items.nEdges_eq hn0]; simp only [Items.capCount, hcap, ↓reduceIte]; omega)
  · have hn' : (items.type i).isNode = false := Bool.eq_false_iff.2 hn
    have hne0 := Items.nEdges_eq_zero (g := g) hn'
    cases hty : items.type i
    all_goals simp only [hty] at hn'
    all_goals simp [NodeType.isNode] at hn'
    · exact local_F _ _ _ _ _ _ (by omega) (by omega)
    · have h0 := H.wf.nvList_V hi hty
      rw [h0] at hnv; simp at hnv
      exact local_V _ _ _ _ _ _ (by omega) (by omega)

/-- If a node's rows start at `2 neSt`, they end at `2 neEn`. -/
theorem adj_bound_step (hor : items.ROriented g) {i : ItemId} (hi : i < items.size)
    (h0 : t.adjBounds[2 * (t.nvRange (idx i)).1]! = 2 * (t.neRange (idx i)).1) :
    t.adjBounds[2 * (t.nvRange (idx i)).2]! = 2 * (t.neRange (idx i)).2 := by
  obtain ⟨pos, hl, hloc⟩ := H.layout_local hor hi
  have hnv := (H.node i hi).nv_range
  have hne := (H.node i hi).ne_range
  by_cases hz : (items.nvList g i).length = 0
  · have hn : (items.type i).isNode = false := by
      by_contra h
      have := H.wf.one_le_nvList hi (by simpa using h); omega
    rw [hnv, hne, hz, Items.nEdges_eq_zero (g := g) hn, Nat.add_zero, Nat.add_zero]
    exact h0
  · have hb := hloc.bound_last
    have ha := hl.adj_bounds (2 * (items.nvList g i).length) (by omega) le_rfl
    unfold LayoutR.rowBound at hb
    rw [hnv] at hb ⊢
    rw [LayoutShape.ite_of_neg (show ¬ (2 * ((t.nvRange (idx i)).1 + (items.nvList g i).length) =
      2 * (t.nvRange (idx i)).1) by omega), show 2 * ((t.nvRange (idx i)).1 + (items.nvList g i).length) -
      2 * (t.nvRange (idx i)).1 = 2 * (items.nvList g i).length by omega] at hb
    rw [show 2 * ((t.nvRange (idx i)).1 + (items.nvList g i).length) =
      2 * (t.nvRange (idx i)).1 + 2 * (items.nvList g i).length by omega, ha]
    exact hb

/-- `adjBounds[2 nvBounds[n]] = 2 neBounds[n]` for every `n ≤ size`. -/
theorem adj_bounds_start (hor : items.ROriented g) :
    ∀ n, n ≤ t.size → t.adjBounds[2 * t.nvBounds[n]!]! = 2 * t.neBounds[n]! := by
  intro n
  induction n with
  | zero => intro _; rw [H.gl.nv_zero, H.gl.ne_zero]; simpa using H.gl.adj_zero
  | succ n ih =>
    intro hn
    obtain ⟨i, hi, rfl⟩ := H.idx_surj (by omega : n < t.size)
    have := H.adj_bound_step hor hi (by rw [t.nvRange_fst, t.neRange_fst]; exact ih (by omega))
    rwa [t.nvRange_snd, t.neRange_snd] at this

/-- `relabel_adj_spec` for the tree described by `RelabelAll`. -/
theorem adj_spec (hor : items.ROriented g) :
    (∀ n, n < t.size → t.adjBounds[2 * (t.nvRange n).1]! = 2 * (t.neRange n).1) ∧
    t.adjBounds[2 * t.nodeVerts.size]! = t.adjDat.size := by
  refine ⟨fun n hn => ?_, ?_⟩
  · rw [t.nvRange_fst, t.neRange_fst]; exact H.adj_bounds_start hor n hn.le
  · rw [← H.gl.nv_last, H.gl.adj_dat_size, ← H.gl.ne_last]
    exact H.adj_bounds_start hor t.size le_rfl

end RelabelAll

/-- `relabel_adj_spec`: the output `adjBounds` CSR facts, from `relabel_node_spec`. -/
theorem relabelTree_adj (g : Graph) (items : Items) (h : items.WF g) (hor : items.ROriented g) :
    (∀ n, n < (relabelTree g items).size →
      (relabelTree g items).adjBounds[2 * ((relabelTree g items).nvRange n).1]! =
        2 * ((relabelTree g items).neRange n).1) ∧
    (relabelTree g items).adjBounds[2 * (relabelTree g items).nodeVerts.size]! =
      (relabelTree g items).adjDat.size := by
  obtain ⟨idx, -, hidx, hnode⟩ := relabel_node_spec g items h
  exact (RelabelAll.mk h hidx hnode).adj_spec hor

alias relabel_adj_spec := relabelTree_adj

end Spqr
