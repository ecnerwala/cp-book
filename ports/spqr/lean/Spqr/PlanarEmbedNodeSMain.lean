import Spqr.PlanarEmbedNodeSFold

/-!
# The `S` node: `nodeFold_capped_S`

The whole node fold of an `S` node is `blockFold` over all `n` blocks; blocks `0..n-2` yield the
`InvS` certificate of the chain of all children (`invS_all`), block `n - 1` is a no-op
(`last_block_S`), and the chain's edge list is a permutation of `edgesBelow i` (every child is a
non-`V` child or the `V` item of an interior node-vertex, `PieceSep.nv_child`/`v_child_nv`), so
`Capped.perm` transports the certificate to `pieceBelow g i`.
-/

namespace Spqr.PlanarSpqrTree

variable (t : PlanarSpqrTree)

section S

variable {i : Nat}

theorem sC_mem_sL (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hS : t.toSpqrTree.type i = .S)
    {j m c : Nat} (hj : j ≤ m) (hc : (t.sC i)[j]? = some c) : c ∈ t.sL i m := by
  obtain ⟨hjl, hjv⟩ := List.getElem?_eq_some_iff.1 hc
  have hc' : (t.sC i)[j]! = c := by rw [getElem!_pos (t.sC i) j hjl, hjv]
  subst hc'
  cases j with
  | zero => exact List.mem_cons.2 (Or.inl rfl)
  | succ j => exact List.mem_cons_of_mem _ (List.mem_flatMap.2 ⟨j, List.mem_range.2 (by omega), by simp⟩)

theorem sW_mem_sL {j m : Nat} (hj1 : 1 ≤ j) (hj : j ≤ m) :
    t.sW ((t.toSpqrTree.nvRange i).1 + j) ∈ t.sL i m := by
  obtain ⟨j, rfl⟩ : ∃ j', j = j' + 1 := ⟨j - 1, by omega⟩
  exact List.mem_cons_of_mem _ (List.mem_flatMap.2 ⟨j, List.mem_range.2 (by omega), by simp⟩)

theorem sL_nodup (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hS : t.toSpqrTree.type i = .S) {m : Nat} (hm : m ≤ t.toSpqrTree.nVerts i - 2) :
    (t.sL i m).Nodup := by
  induction m with
  | zero => simp [sL]
  | succ m ih =>
    rw [sL_succ]
    refine (ih (by omega)).append ?_ ?_
    · exact List.nodup_cons.2
        ⟨by simpa using t.W_ne_C (g := g) hwf hsep hi hS (j := m + 1) (by omega) (by omega),
         List.nodup_singleton _⟩
    · intro c' hc' hc2
      simp only [List.mem_cons, List.mem_singleton, List.not_mem_nil, or_false] at hc2
      rcases hc2 with rfl | rfl
      · exact t.ne_sL_W (g := g) hwf hsep hi hS hm hc' rfl
      · exact t.ne_sL_C (g := g) hwf hsep hi hS hm hc' rfl

/-- Both node-vertices of the cap of an `S` node are cap endpoints. -/
theorem capEnd_S (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hS : t.toSpqrTree.type i = .S) :
    t.toSpqrTree.CapEnd i (t.toSpqrTree.nvRange i).1 ∧
    t.toSpqrTree.CapEnd i ((t.toSpqrTree.nvRange i).1 + t.toSpqrTree.nVerts i - 1) := by
  have h3 := (t.shape_S hwf hi hS).1
  have hn := t.nEdges_S hwf hi hS
  have hsz : (t.toSpqrTree.neRange i).1 < t.nodeEdges.size := by
    have := PlanarRot.neRange_le_size t.toSpqrTree hwf i hi; omega
  have h0 := t.nvs_S hwf hi hS (k := 0) (by omega)
  rw [Nat.add_zero, if_pos rfl, getElem!_pos t.nodeEdges _ hsz] at h0
  exact ⟨⟨_, _, t.capNe_S hS, Array.getElem?_eq_getElem hsz, Or.inl (by rw [h0])⟩,
    ⟨_, _, t.capNe_S hS, Array.getElem?_eq_getElem hsz, Or.inr (by rw [h0])⟩⟩

/-- Every child of an `S` node is a non-`V` child or the `V` item of an interior node-vertex. -/
theorem mem_sL_of_children (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hS : t.toSpqrTree.type i = .S) {c : Nat} (hc : c ∈ t.children i) :
    c ∈ t.sL i (t.toSpqrTree.nVerts i - 2) := by
  have h3 := (t.shape_S hwf hi hS).1
  by_cases hV : t.toSpqrTree.type c = .V
  · obtain ⟨nv, d, hnv, hd, hdv⟩ :=
      hsep.v_child_nv i c hi (Or.inl hS) (by rw [t.children_eq]; exact hc) hV
    have hpar := (t.child_data hwf hi hc).2.1
    have hnc : ¬ t.toSpqrTree.CapEnd i nv :=
      (hsep.nv_child i nv d hi (Or.inl hS) hnv hd).1 (by rw [hdv]; exact hpar)
    obtain ⟨hcap0, hcapl⟩ := t.capEnd_S hwf hi hS
    have hnv1 := hnv.1
    have hnv2 := hnv.2
    have hn : t.toSpqrTree.nVerts i = (t.toSpqrTree.nvRange i).2 - (t.toSpqrTree.nvRange i).1 := rfl
    obtain ⟨j, rfl⟩ : ∃ j, nv = (t.toSpqrTree.nvRange i).1 + j :=
      ⟨nv - (t.toSpqrTree.nvRange i).1, by omega⟩
    have hj1 : 1 ≤ j := by
      by_contra h0
      exact hnc (by rw [show j = 0 by omega, Nat.add_zero]; exact hcap0)
    have hj : j ≤ t.toSpqrTree.nVerts i - 2 := by
      by_contra h0
      exact hnc (by
        rw [show (t.toSpqrTree.nvRange i).1 + j =
          (t.toSpqrTree.nvRange i).1 + t.toSpqrTree.nVerts i - 1 by omega]
        exact hcapl)
    have hc' : c = t.sW ((t.toSpqrTree.nvRange i).1 + j) := by
      rw [← hdv]; unfold sW; rw [getElem!_of_getElem? hd]
    rw [hc']
    exact t.sW_mem_sL hj1 hj
  · have hcS : c ∈ t.sC i := List.mem_filter.2 ⟨hc, by simpa using hV⟩
    obtain ⟨j, hj⟩ := List.mem_iff_getElem?.1 hcS
    exact t.sC_mem_sL hwf hi hS (by have := (List.getElem?_eq_some_iff.1 hj).1; rw [t.sC_length hwf hi hS] at this; omega) hj

theorem children_perm_sL (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hS : t.toSpqrTree.type i = .S) : (t.children i).Perm (t.sL i (t.toSpqrTree.nVerts i - 2)) :=
  (List.perm_ext_iff_of_nodup (t.children_nodup hwf hi) (t.sL_nodup (g := g) hwf hsep hi hS le_rfl)).2
    fun _ => ⟨t.mem_sL_of_children (g := g) hwf hsep hi hS, t.mem_sL_children (g := g) hwf hsep hi hS le_rfl⟩

theorem chain_perm_edgesBelow (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size) (hS : t.toSpqrTree.type i = .S) :
    (t.chain g i (t.toSpqrTree.nVerts i - 2)).ves.Perm (t.edgesBelow i) := by
  have he : t.edgesBelow i = (t.children i).flatMap t.edgesBelow := by
    rw [t.edgesBelow_eq hwf hsh i hi, ← t.type_eq_of_lt i hi, hS]; rfl
  rw [he]
  exact (t.children_perm_sL (g := g) hwf hsep hi hS).symm.flatMap_right t.edgesBelow

theorem fold_eq_blockFold_S (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hS : t.toSpqrTree.type i = .S)
    (s : EmbedState) :
    (List.range' (4 * t.neBounds[i]!) (4 * t.neBounds[i + 1]! - 4 * t.neBounds[i]!)).foldl
      (t.nodeStep i t.neBounds[i]!) s =
      t.blockFold i (t.toSpqrTree.neRange i).1 s (t.toSpqrTree.nVerts i - 2 + 1) := by
  have h3 := (t.shape_S hwf hi hS).1
  have hn := t.nEdges_S hwf hi hS
  rw [PlanarRot.neRange_eq] at hn ⊢
  unfold blockFold
  dsimp only at hn ⊢
  rw [show 4 * t.neBounds[i + 1]! - 4 * t.neBounds[i]! = 4 * (t.toSpqrTree.nVerts i - 2 + 1 + 1) by omega]

/-- `nodeFold_capped` for `S` nodes: iterated `Capped.join`/`attachOpen` along the cycle layout. -/
theorem nodeFold_capped_S (g : Graph) (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hrep : t.toSpqrTree.Represents g) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hS : t.toSpqrTree.type i = .S)
    (hlay : ∃ ev mr cv, t.neRotAdj.extract (4 * (t.toSpqrTree.neRange i).1) (4 * (t.toSpqrTree.neRange i).2) =
      layoutRot (t.toSpqrTree.type i) (t.toSpqrTree.nVerts i) (t.toSpqrTree.neRange i).1
        (t.toSpqrTree.neRange i).2 ev mr cv)
    {ne : Nat} (hcne : t.toSpqrTree.capNe i = some ne) {p : Nat × Nat}
    (hp : t.toSpqrTree.neOrig ne = some p)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    let s' := (List.range' (4 * t.neBounds[i]!) (4 * t.neBounds[i + 1]! - 4 * t.neBounds[i]!)).foldl
      (t.nodeStep i t.neBounds[i]!) s
    (∀ q, ¬(t.pieceBelow g i).Mem q → s'.rotAdj[q]? = s.rotAdj[q]?) ∧
    ∃ a0 a1 a2 a3 ρ, s'.outerE[i]? = some #[some a0, some a1, some a2, some a3] ∧
      (t.pieceBelow g i).Capped s'.rotAdj ρ a0 a1 a2 a3 p.1 p.2 := by
  dsimp only
  have hne := hrep.ne
  have h3 := (t.shape_S hwf hi hS).1
  rw [t.fold_eq_blockFold_S hwf hi hS, t.last_block_S hwf hi hS hlay]
  obtain ⟨houter, hrsz, hframe, a0, a1, a2, a3, ρ, ha0, ha1, ha2, ha3, hcap⟩ :=
    t.invS_all hwf hsh hne hsep hi hS hlay s h (t.toSpqrTree.nVerts i - 2) le_rfl
  obtain ⟨-, b0, b1, e2, e3, hb0, hb1, he2, he3, hrow⟩ := t.cap_block_S hwf hne hsep hi hS hlay s h
  rw [ha0] at hb0; rw [ha1] at hb1; rw [ha2] at he2; rw [ha3] at he3
  cases hb0; cases hb1; cases he2; cases he3
  have hperm := t.chain_perm_edgesBelow (g := g) hwf hsh hsep hi hS
  rw [t.capNe_S hS] at hcne
  cases hcne
  have hp' := t.capOrig_S (g := g) hwf hsep hi hS hp
  subst hp'
  refine ⟨fun q hq => hframe q fun hq' => hq (hperm.mem_iff.1 hq'), a0, a1, a2, a3, ?_⟩
  obtain ⟨σ, hσ⟩ := hcap.perm hperm (hperm.nodup_iff.2 (t.edgesBelow_nodup hwf i))
  refine ⟨σ, ?_, ?_⟩
  · rw [t.blockFold_outerE, hrow]
  · rw [show t.toSpqrTree.nVerts i - 2 + 1 = t.toSpqrTree.nVerts i - 1 by omega] at hσ
    exact hσ

end S

end Spqr.PlanarSpqrTree
