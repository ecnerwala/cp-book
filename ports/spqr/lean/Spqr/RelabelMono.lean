import Spqr.RelabelInv

/-!
# Transport of a node's record along `Agree`

`RelabelLayout.congr`: the layout record only reads the slots listed in `LowB`, so it transports
to any tree agreeing there. `NodeS.mono`: a finished node's record survives any later step that
rewrites nothing at or above its `LowB` bounds.
-/

namespace Spqr

namespace Items

variable {g : Graph} {items : Items} {i c : ItemId} {nvSt : Nat} {pos : Nat → Nat}

theorem ordered_perm : (items.ordered g i nvSt pos).Perm (items.ch i) := by
  unfold ordered; split
  · exact List.Perm.refl _
  · exact List.mergeSort_perm _ _

theorem mem_ordered_iff : c ∈ items.ordered g i nvSt pos ↔ c ∈ items.ch i := ordered_perm.mem_iff

theorem mem_ch_of_mem_filter_ordered {p : ItemId → Bool} (h : c ∈ (items.ordered g i nvSt pos).filter p) :
    c ∈ items.ch i := mem_ordered_iff.1 (List.mem_of_mem_filter h)

theorem length_filter_ordered (p : ItemId → Bool) :
    ((items.ordered g i nvSt pos).filter p).length = (items.ch i).countP p := by
  rw [← List.countP_eq_length_filter, ordered_perm.countP_eq]

theorem capCount_add_lt_nEdges (hn : (items.type i).isNode) {k : Nat}
    (hk : k < ((items.ordered g i nvSt pos).filter (· ≥ 1 + g.nv)).length) :
    items.capCount i + k < items.nEdges g i := by
  rw [length_filter_ordered] at hk
  unfold nEdges; rw [if_pos hn]; omega

end Items

theorem nodeLayout_congr {g : Graph} {items : Items} {t t' : SpqrTree} {idx idx' : ItemId → Nat} {i : ItemId}
    {pos : Nat → Nat} (hi : idx' i = idx i) (hnv : t'.nvRange (idx i) = t.nvRange (idx i))
    (hne : t'.neRange (idx i) = t.neRange (idx i)) :
    nodeLayout g items t' idx' i pos = nodeLayout g items t idx i pos := by
  unfold nodeLayout; rw [hi, hnv, hne]

/-- `RelabelLayout` only reads the node's ranges, child slots, node-edges, adjacency rows, the
children's `vertParNv` / `neRange` / cap twin / `subtreeEnd`, and the node's own `subtreeEnd`. -/
theorem RelabelLayout.congr {g : Graph} {items : Items} {t t' : SpqrTree} {idx idx' : ItemId → Nat}
    {i : ItemId} {pos : Nat → Nat} (h : RelabelLayout g items t idx i pos)
    (hi : idx' i = idx i) (hch : ∀ c ∈ items.ch i, idx' c = idx c)
    (hnv : t'.nvRange (idx i) = t.nvRange (idx i)) (hne : t'.neRange (idx i) = t.neRange (idx i))
    (hchildren : t'.children (idx i) = t.children (idx i))
    (hedges : ∀ k, k < items.nEdges g i →
      (t'.nodeEdges[(t.neRange (idx i)).1 + k]!).node = (t.nodeEdges[(t.neRange (idx i)).1 + k]!).node ∧
      (t'.nodeEdges[(t.neRange (idx i)).1 + k]!).nvs = (t.nodeEdges[(t.neRange (idx i)).1 + k]!).nvs)
    (hadjb : ∀ j, 1 ≤ j → j ≤ 2 * (items.nvList g i).length →
      t'.adjBounds[2 * (t.nvRange (idx i)).1 + j]! = t.adjBounds[2 * (t.nvRange (idx i)).1 + j]!)
    (hadjd : ∀ j, j < 2 * items.nEdges g i →
      t'.adjDat[2 * (t.neRange (idx i)).1 + j]! = t.adjDat[2 * (t.neRange (idx i)).1 + j]!)
    (hvpn : ∀ c ∈ items.ch i, t'.vertParNv[idx c]! = t.vertParNv[idx c]!)
    (htwin : ∀ k, k < items.nEdges g i → t'.twin ((t.neRange (idx i)).1 + k) = t.twin ((t.neRange (idx i)).1 + k))
    (hcne : ∀ c ∈ items.ch i, t'.neRange (idx c) = t.neRange (idx c))
    (hctwin : ∀ c ∈ items.ch i, items.hasCap c → t'.twin (t.neRange (idx c)).1 = t.twin (t.neRange (idx c)).1)
    (hse : ∀ c ∈ items.ch i, t'.subtreeEnd[idx c]! = t.subtreeEnd[idx c]!)
    (hsei : t'.subtreeEnd[idx i]! = t.subtreeEnd[idx i]!) :
    RelabelLayout g items t' idx' i pos := by
  have hL := nodeLayout_congr (g := g) (items := items) (pos := pos) hi hnv hne
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> simp only [hi, hnv, hne, hL, hchildren, hsei]
  · exact h.pos_ok
  · rw [h.children]
    exact List.map_congr_left fun c hc => (hch c (Items.mem_ordered_iff.1 hc)).symm
  · intro k hk; rw [(hedges k hk).1]; exact h.edge_node k hk
  · intro k hk; rw [(hedges k hk).2]; exact h.edge_nvs k hk
  · intro j h1 h2; rw [hadjb j h1 h2]; exact h.adj_bounds j h1 h2
  · intro j hj; rw [hadjd j hj]; exact h.adj_dat j hj
  · intro k hk
    have hc := Items.mem_ch_of_mem_filter_ordered (List.getElem_mem hk)
    rw [hch _ hc, hvpn _ hc]; exact h.vert_par_nv k hk
  · intro hn k hk
    obtain ⟨h1, h2⟩ := h.twin hn k hk
    have hc := Items.mem_ch_of_mem_filter_ordered (List.getElem_mem hk)
    have ht := htwin _ (Items.capCount_add_lt_nEdges hn hk)
    rw [← Nat.add_assoc] at ht
    rw [hch _ hc, hcne _ hc, ht]
    exact ⟨h1, fun hcap => by rw [hctwin _ hc hcap]; exact h2 hcap⟩
  · intro k hk
    have hc := Items.mem_ordered_iff.1 (List.getElem_mem hk)
    rw [hch _ hc, h.child_idx k hk]
    split
    · rfl
    · have hc' := Items.mem_ordered_iff.1 (List.getElem_mem (l := items.ordered g i (t.nvRange (idx i)).1 pos)
        (n := k - 1) (by omega))
      rw [hch _ hc', hse _ hc']
  · rw [h.subtree_end]
    cases hl : (items.ordered g i (t.nvRange (idx i)).1 pos).getLast?
    · rfl
    · have hc := Items.mem_ordered_iff.1 (List.mem_of_getLast? hl)
      simp only [hch _ hc, hse _ hc]

namespace Ghost

variable {g : Graph} {items : Items} {B : Bounds} {s s' : RelabelState} {i : ItemId}

section tree
variable (g : Graph) (s : RelabelState) (n : Nat)
theorem tree_chRange_fst : ((s.tree g).chRange n).1 = s.chBounds[n]! := by
  rw [tree_chRange, Array.getElem!_eq_getD_getElem?]; rfl
theorem tree_chRange_snd : ((s.tree g).chRange n).2 = s.chBounds[n + 1]! := by
  rw [tree_chRange, Array.getElem!_eq_getD_getElem?]; rfl
theorem tree_nvRange_fst : ((s.tree g).nvRange n).1 = s.nvBounds[n]! := by
  rw [tree_nvRange, Array.getElem!_eq_getD_getElem?]; rfl
theorem tree_nvRange_snd : ((s.tree g).nvRange n).2 = s.nvBounds[n + 1]! := by
  rw [tree_nvRange, Array.getElem!_eq_getD_getElem?]; rfl
theorem tree_neRange_fst : ((s.tree g).neRange n).1 = s.neBounds[n]! := by
  rw [tree_neRange, Array.getElem!_eq_getD_getElem?]; rfl
theorem tree_neRange_snd : ((s.tree g).neRange n).2 = s.neBounds[n + 1]! := by
  rw [tree_neRange, Array.getElem!_eq_getD_getElem?]; rfl
theorem tree_parent_eq : (s.tree g).parent n = s.par[n]! := by
  rw [tree_parent, Array.getElem!_eq_getD_getElem?]; rfl
end tree

theorem Consistent.idx_lt (hc : Consistent g items s) (hi : i ∈ s.order.toList) : s.idx i < s.types.size :=
  (hc.idx_lt_iff i).2 hi

namespace Agree

variable (ha : Agree B s s') (hc : Consistent g items s)
include ha hc

omit hc in
theorem types_getD {n : Nat} (hn : n < s.types.size) : s'.types[n]?.getD .F = s.types[n]?.getD .F :=
  ha.types.getD (Nat.zero_le _) hn _
theorem par_get! {n : Nat} (hn : n < s.types.size) : s'.par[n]! = s.par[n]! :=
  ha.par.get! (Nat.zero_le _) (by rw [hc.par_size]; exact hn)
theorem origId_get! {n : Nat} (hn : n < s.types.size) : s'.origId[n]! = s.origId[n]! :=
  ha.origId.get! (Nat.zero_le _) (by rw [hc.origId_size]; exact hn)
theorem vertParNv_get! {n : Nat} (hn : n < s.types.size) : s'.vertParNv[n]! = s.vertParNv[n]! :=
  ha.vertParNv.get! (Nat.zero_le _) (by rw [hc.vertParNv_size]; exact hn)
theorem subtreeEnd_get! {n : Nat} (hb : B.idx ≤ n) (hn : n < s.types.size) : s'.subtreeEnd[n]! = s.subtreeEnd[n]! :=
  ha.subtreeEnd.get! hb (by rw [hc.subtreeEnd_size]; exact hn)
theorem chBounds_get! {n : Nat} (hn : n ≤ s.types.size) : s'.chBounds[n]! = s.chBounds[n]! :=
  ha.chBounds.get! (Nat.zero_le _) (by rw [hc.chBounds_size]; omega)
theorem nvBounds_get! {n : Nat} (hn : n ≤ s.types.size) : s'.nvBounds[n]! = s.nvBounds[n]! :=
  ha.nvBounds.get! (Nat.zero_le _) (by rw [hc.nvBounds_size]; omega)
theorem neBounds_get! {n : Nat} (hn : n ≤ s.types.size) : s'.neBounds[n]! = s.neBounds[n]! :=
  ha.neBounds.get! (Nat.zero_le _) (by rw [hc.neBounds_size]; omega)
theorem chRange_eq {n : Nat} (hn : n < s.types.size) : (s'.tree g).chRange n = (s.tree g).chRange n := by
  rw [tree_chRange, tree_chRange, ha.chBounds.getD (Nat.zero_le _) (by rw [hc.chBounds_size]; omega),
    ha.chBounds.getD (Nat.zero_le _) (by rw [hc.chBounds_size]; omega)]
theorem nvRange_eq {n : Nat} (hn : n < s.types.size) : (s'.tree g).nvRange n = (s.tree g).nvRange n := by
  rw [tree_nvRange, tree_nvRange, ha.nvBounds.getD (Nat.zero_le _) (by rw [hc.nvBounds_size]; omega),
    ha.nvBounds.getD (Nat.zero_le _) (by rw [hc.nvBounds_size]; omega)]
theorem neRange_eq {n : Nat} (hn : n < s.types.size) : (s'.tree g).neRange n = (s.tree g).neRange n := by
  rw [tree_neRange, tree_neRange, ha.neBounds.getD (Nat.zero_le _) (by rw [hc.neBounds_size]; omega),
    ha.neBounds.getD (Nat.zero_le _) (by rw [hc.neBounds_size]; omega)]
omit hc in
theorem twin_eq {n : Nat} (hb : B.nodeEdges ≤ n) (hn : n < s.nodeEdges.size) :
    (s'.tree g).twin n = (s.tree g).twin n := by
  rw [tree_twin, tree_twin, ha.nodeEdges.get? hb hn]
theorem adjBounds_get! {n : Nat} (hn : n ≤ 2 * s.nodeVerts.size) : s'.adjBounds[n]! = s.adjBounds[n]! :=
  ha.adjBounds.get! (Nat.zero_le _) (by rw [hc.adjBounds_size]; omega)
theorem adjDat_get! {n : Nat} (hn : n < 2 * s.nodeEdges.size) : s'.adjDat[n]! = s.adjDat[n]! :=
  ha.adjDat.get! (Nat.zero_le _) (by rw [hc.adjDat_size]; exact hn)
omit hc in
theorem nodeVerts_get! {n : Nat} (hn : n < s.nodeVerts.size) : s'.nodeVerts[n]! = s.nodeVerts[n]! :=
  ha.nodeVerts.get! (Nat.zero_le _) hn

/-- The children list of a finished node is unchanged. -/
theorem children_eq {n : Nat} (hn : n < s.types.size) (hb : B.chDat ≤ s.chBounds[n]!)
    (hle : s.chBounds[n + 1]! ≤ s.chDat.size) : (s'.tree g).children n = (s.tree g).children n := by
  unfold SpqrTree.children
  rw [ha.chRange_eq hc hn]
  refine List.map_congr_left fun k hk => ?_
  rw [List.mem_range] at hk
  rw [tree_chDat, tree_chDat, tree_chRange_fst, tree_chRange_snd] at *
  rw [ha.chDat.get? (k := s.chBounds[n]! + k) (by omega) (by omega)]

end Agree

theorem LowB.agree {B' : Bounds} (hc : Consistent g items s) (hl : LowB items B s i) (ha : Agree B' s s') :
    LowB items B s' i := by
  have e := ha.idx_eq hl.mem
  have hlt := hc.idx_lt hl.mem
  have hche : ∀ c ∈ items.ch i, s'.idx c = s.idx c := fun c hc' => ha.idx_eq (hl.ch_mem c hc')
  have hchlt : ∀ c ∈ items.ch i, s.idx c < s.types.size := fun c hc' => hc.idx_lt (hl.ch_mem c hc')
  refine ⟨ha.order.subset hl.mem, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · rw [e]; exact hl.idx
  · rw [e, ha.chBounds_get! hc hlt.le]; exact hl.chDat
  · rw [e, ha.chBounds_get! hc hlt]; exact hl.chDat_le.trans ha.chDat.size_le
  · rw [e, ha.neBounds_get! hc hlt.le]; exact hl.nodeEdges
  · rw [e, ha.neBounds_get! hc hlt]; exact hl.nodeEdges_le.trans ha.nodeEdges.size_le
  · rw [e, ha.nvBounds_get! hc hlt]; exact hl.nodeVerts_le.trans ha.nodeVerts.size_le
  · exact fun c hc' => ha.order.subset (hl.ch_mem c hc')
  · intro c hc'; rw [hche c hc']; exact hl.ch_idx c hc'
  · intro c hc'; rw [hche c hc', ha.neBounds_get! hc (hchlt c hc').le]
    exact ⟨(hl.ch_ne c hc').1, fun hcap => ((hl.ch_ne c hc').2 hcap).trans_le ha.nodeEdges.size_le⟩

/-- A finished node's record survives a step that keeps everything at or above its bounds. -/
theorem NodeS.mono (hc : Consistent g items s) (hn : NodeS g items s i) (hl : LowB items B s i)
    (ha : Agree B s s') : NodeS g items s' i := by
  have e := ha.idx_eq hl.mem
  have hlt := hc.idx_lt hl.mem
  have hche : ∀ c ∈ items.ch i, s'.idx c = s.idx c := fun c hc' => ha.idx_eq (hl.ch_mem c hc')
  have hchlt : ∀ c ∈ items.ch i, s.idx c < s.types.size := fun c hc' => hc.idx_lt (hl.ch_mem c hc')
  obtain ⟨node, raw⟩ := hn
  have hnvr := ha.nvRange_eq hc (g := g) hlt
  have hner := ha.neRange_eq hc (g := g) hlt
  have hchr := ha.chRange_eq hc (g := g) hlt
  have nvSt_eq := tree_nvRange_fst g s (s.idx i)
  have nvEn_eq := tree_nvRange_snd g s (s.idx i)
  have neSt_eq := tree_neRange_fst g s (s.idx i)
  have neEn_eq := tree_neRange_snd g s (s.idx i)
  have hnvEn : s.nvBounds[s.idx i + 1]! = s.nvBounds[s.idx i]! + (items.nvList g i).length := by
    have := node.nv_range; rwa [nvSt_eq, nvEn_eq] at this
  have hneEn : s.neBounds[s.idx i + 1]! = s.neBounds[s.idx i]! + items.nEdges g i := by
    have := node.ne_range; rwa [neSt_eq, neEn_eq] at this
  have hnvle := hl.nodeVerts_le
  have hnele := hl.nodeEdges_le
  -- raw node-verts transport
  have raw' : ∀ k (hk : k < (items.nvList g i).length),
      s'.nodeVerts[((s'.tree g).nvRange (s'.idx i)).1 + k]! = ⟨s'.idx i, (items.nvList g i)[k]⟩ := by
    intro k hk
    rw [e, hnvr, ha.nodeVerts_get! (by rw [nvSt_eq]; omega)]
    exact raw k hk
  refine ⟨?_, raw'⟩
  have hnvlt : ∀ k, k < (items.nvList g i).length → ((s.tree g).nvRange (s.idx i)).1 + k < s.nodeVerts.size := by
    intro k hk; rw [nvSt_eq]; omega
  have hnvlt' : ∀ k, k < (items.nvList g i).length → ((s.tree g).nvRange (s.idx i)).1 + k < s'.nodeVerts.size :=
    fun k hk => (hnvlt k hk).trans_le ha.nodeVerts.size_le
  refine
    { idx_lt := ?_, type := ?_, orig := ?_, vert_index := ?_, edge_index := ?_, ch_range := ?_, nv_range := ?_
      ne_range := ?_, node_verts := ?_, child_lt := ?_, child_par := ?_, child_par_nv_none := ?_
      child_cap_twin_none := ?_, layout := ?_ }
  · rw [e, tree_size]; exact hlt.trans_le ha.types.size_le
  · rw [e, tree_type, ha.types_getD hlt]; exact node.type
  · rw [e, tree_origId, ha.origId_get! hc hlt]; exact node.orig
  · intro hV; rw [e, tree_vertIndex]
    have := node.vert_index hV
    rw [tree_vertIndex] at this
    rw [ha.vertIndex _ (by rw [this]; exact Option.some_ne_none _)]; exact this
  · intro hQ; rw [e, tree_edgeIndex, tree_edgeFlipped]
    obtain ⟨h1, h2⟩ := node.edge_index hQ
    rw [tree_edgeIndex] at h1; rw [tree_edgeFlipped] at h2
    have hne : s.edgeIndex[i - 1 - g.nv]! ≠ none := by rw [h1]; exact Option.some_ne_none _
    rw [ha.edgeIndex _ hne, ha.edgeFlipped _ hne]; exact ⟨h1, h2⟩
  · rw [e, hchr]; exact node.ch_range
  · rw [e, hnvr]; exact node.nv_range
  · rw [e, hner]; exact node.ne_range
  · intro k hk
    rw [e, hnvr, tree_nodeVerts_getElem! g s' _ (hnvlt' k hk), ha.nodeVerts_get! (hnvlt k hk)]
    have := raw k hk
    rw [this]; rfl
  · intro c hc'; rw [hche c hc', tree_size]; exact (hchlt c hc').trans_le ha.types.size_le
  · intro c hc'; rw [hche c hc', e, tree_parent_eq, ha.par_get! hc (hchlt c hc'), ← tree_parent_eq]
    exact node.child_par c hc'
  · intro c hc' hnv; rw [hche c hc', tree_vertParNv, ha.vertParNv_get! hc (hchlt c hc')]
    exact node.child_par_nv_none c hc' hnv
  · intro hnn c hc' hcap
    have := node.child_cap_twin_none hnn c hc' hcap
    rw [hche c hc', ha.neRange_eq hc (g := g) (hchlt c hc')]
    rw [tree_neRange_fst] at this ⊢
    rw [ha.twin_eq (hl.ch_ne c hc').1 ((hl.ch_ne c hc').2 hcap)]
    exact this
  · obtain ⟨pos, hL⟩ := node.layout
    refine ⟨pos, hL.congr e hche hnvr hner (ha.children_eq hc hlt ?_ ?_) ?_ ?_ ?_ ?_ ?_ ?_ ?_ ?_ ?_⟩
    · rw [← tree_chRange_fst g s]; rw [tree_chRange_fst]; exact hl.chDat
    · exact hl.chDat_le
    · intro k hk
      rw [tree_nodeEdges, tree_nodeEdges, neSt_eq]
      exact ⟨ha.nodeEdges_node _ (by omega), ha.nodeEdges_nvs _ (by omega)⟩
    · intro j h1 h2
      rw [tree_adjBounds, tree_adjBounds, nvSt_eq, ha.adjBounds_get! hc (by omega)]
    · intro j hj
      rw [tree_adjDat, tree_adjDat, neSt_eq, ha.adjDat_get! hc (by omega)]
    · intro c hc'; rw [tree_vertParNv, tree_vertParNv, ha.vertParNv_get! hc (hchlt c hc')]
    · intro k hk
      rw [neSt_eq, ha.twin_eq (hl.nodeEdges.trans (Nat.le_add_right _ _)) (by omega)]
    · intro c hc'; exact ha.neRange_eq hc (hchlt c hc')
    · intro c hc' hcap
      rw [tree_neRange_fst, ha.twin_eq (hl.ch_ne c hc').1 ((hl.ch_ne c hc').2 hcap)]
    · intro c hc'; rw [tree_subtreeEnd, tree_subtreeEnd, ha.subtreeEnd_get! hc (hl.ch_idx c hc') (hchlt c hc')]
    · rw [tree_subtreeEnd, tree_subtreeEnd, ha.subtreeEnd_get! hc hl.idx hlt]

end Ghost

end Spqr
