import Spqr.PlanarEmbedNodeP

/-!
# The `R` node: skeleton, twins, local rotation and original vertices

Bookkeeping for `nodeFold_capped_R`: the node-edges of an `R` node `i` are the cap followed by
`nEdges i - 1` edges whose twins are the caps of the non-`V` children in order (`child_R`), the
`neRotAdj` entries of the node are `nodeRot i` shifted by `4 * neSt` (`rotR`, using
`NodeRotClosed`), and every node-edge has original endpoints `sV` of its node-vertices
(`neOrig_R`).
-/

namespace Spqr.PlanarSpqrTree

variable (t : PlanarSpqrTree) {g : Graph} {i : Nat}

section R

theorem shape_R (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hR : t.toSpqrTree.type i = .R) :
    4 ≤ t.toSpqrTree.nVerts i ∧ 6 ≤ t.toSpqrTree.nEdges i ∧ (t.toSpqrTree.skeleton i).Nodup ∧
    ∀ p ∈ t.toSpqrTree.skeleton i, p.1 < p.2 := by
  have hsh := hwf.shape i hi
  unfold SpqrTree.Shape at hsh
  rw [hR] at hsh
  obtain ⟨⟨s, e⟩, hnv⟩ : ∃ p, t.toSpqrTree.nvRange i = p := ⟨_, rfl⟩
  rw [hnv] at hsh
  simp only at hsh
  unfold SpqrTree.nVerts
  rw [hnv]
  exact hsh

theorem hasCap_R (hR : t.toSpqrTree.type i = .R) : t.toSpqrTree.hasCap i = true := by
  simp [SpqrTree.hasCap, hR]; decide

theorem capNe_R (hR : t.toSpqrTree.type i = .R) :
    t.toSpqrTree.capNe i = some (t.toSpqrTree.neRange i).1 := by
  simp [SpqrTree.capNe, t.hasCap_R hR]

theorem neEn_R (hwf : t.toSpqrTree.WF) (hi : i < t.size) :
    (t.toSpqrTree.neRange i).2 = (t.toSpqrTree.neRange i).1 + t.toSpqrTree.nEdges i := by
  have := PlanarRot.neRange_mono t.toSpqrTree hwf i hi
  unfold SpqrTree.nEdges; omega

/-- The non-cap node-edges of an `R` node are the twins of the caps of its non-`V` children, in
order. -/
theorem noncap_R (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hR : t.toSpqrTree.type i = .R) :
    (t.sC i).length = t.toSpqrTree.nEdges i - 1 ∧
    ∀ k, k < t.toSpqrTree.nEdges i - 1 →
      ∃ c, (t.sC i)[k]? = some c ∧
        (t.nodeEdges[(t.toSpqrTree.neRange i).1 + (k + 1)]!).twin = t.toSpqrTree.capNe c := by
  have hnc := hwf.twins.noncap_children i hi (by rw [hR]; decide)
  rw [t.hasCap_R hR, t.children_eq] at hnc
  simp only [↓reduceIte] at hnc
  have hlen := congrArg List.length hnc
  simp only [List.length_map, List.length_drop, t.nodeEdgesOf_length] at hlen
  unfold SpqrTree.nEdges
  refine ⟨by unfold sC; omega, ?_⟩
  intro k hk
  have hget := congrArg (fun l => l[k]?) hnc
  simp only [List.getElem?_map, List.getElem?_drop, t.nodeEdgesOf_getElem? i (1 + k) (by omega),
    Option.map_some] at hget
  obtain ⟨c, hc, hce⟩ := Option.map_eq_some_iff.1 hget.symm
  refine ⟨c, hc, ?_⟩
  rw [t.nodeEdges_getElem! hwf hi (by omega), show k + 1 = 1 + k by omega, hce]

theorem sC_length_R (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hR : t.toSpqrTree.type i = .R) :
    (t.sC i).length = t.toSpqrTree.nEdges i - 1 := (t.noncap_R hwf hi hR).1

theorem sC_getElem?_R (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hR : t.toSpqrTree.type i = .R)
    {j : Nat} (hj : j < t.toSpqrTree.nEdges i - 1) : (t.sC i)[j]? = some (t.sC i)[j]! := by
  have := t.sC_length_R hwf hi hR
  rw [getElem!_pos (t.sC i) j (by omega), List.getElem?_eq_getElem (by omega)]

/-- Facts about the `j`-th non-`V` child `c` of an `R` node (cf. `child_S`, `child_P`). -/
theorem child_R (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hR : t.toSpqrTree.type i = .R) {j c : Nat} (hj : j < t.toSpqrTree.nEdges i - 1)
    (hc : (t.sC i)[j]? = some c) :
    c ∈ t.children i ∧ t.toSpqrTree.type c ≠ .V ∧ t.toSpqrTree.hasCap c = true ∧
    t.toSpqrTree.type c ≠ .I ∧ t.toSpqrTree.type c ≠ .O ∧
    t.toSpqrTree.twin ((t.toSpqrTree.neRange i).1 + (j + 1)) = some (t.toSpqrTree.neRange c).1 ∧
    t.toSpqrTree.neOrig ((t.toSpqrTree.neRange i).1 + (j + 1)) =
      t.toSpqrTree.neOrig (t.toSpqrTree.neRange c).1 ∧
    ∀ r, r < 4 → ∀ s : EmbedState,
      t.treeQe s (4 * ((t.toSpqrTree.neRange i).1 + (j + 1)) + r) = s.outerE[c]![r]! := by
  have hmem := List.mem_of_getElem? hc
  rw [sC, List.mem_filter] at hmem
  have hcm : c ∈ t.children i := hmem.1
  have hcV : t.toSpqrTree.type c ≠ .V := by simpa using hmem.2
  obtain ⟨hcap, hcI, hcO⟩ := hsep.node_child_cap i c hi (Or.inr (Or.inr hR)) (by rwa [t.children_eq]) hcV
  obtain ⟨c', hc', htw⟩ := (t.noncap_R hwf hi hR).2 j hj
  rw [hc] at hc'; cases hc'
  have hcne : t.toSpqrTree.capNe c = some (t.toSpqrTree.neRange c).1 := by
    simp [SpqrTree.capNe, hcap]
  rw [hcne] at htw
  have hbnd : (t.toSpqrTree.neRange i).1 + (j + 1) < (t.toSpqrTree.neRange i).2 := by
    unfold SpqrTree.nEdges at hj; omega
  have hsz := PlanarRot.neRange_le_size t.toSpqrTree hwf i hi
  have htw' : t.toSpqrTree.twin ((t.toSpqrTree.neRange i).1 + (j + 1)) =
      some (t.toSpqrTree.neRange c).1 := by
    unfold SpqrTree.twin
    rw [Array.getElem?_eq_getElem (show (t.toSpqrTree.neRange i).1 + (j + 1) < t.nodeEdges.size by omega),
      ← getElem!_pos t.nodeEdges ((t.toSpqrTree.neRange i).1 + (j + 1)) (by omega), Option.bind_some, htw]
  refine ⟨hcm, hcV, hcap, hcI, hcO, htw', hsep.twin_orient _ _ htw', ?_⟩
  intro r hr s
  have hcs := (t.child_data hwf hi hcm).1
  unfold treeQe
  have he : QE.edge (4 * ((t.toSpqrTree.neRange i).1 + (j + 1)) + r) =
      (t.toSpqrTree.neRange i).1 + (j + 1) := by unfold QE.edge; omega
  have hslot : 2 * QE.side (4 * ((t.toSpqrTree.neRange i).1 + (j + 1)) + r) +
      QE.dir (4 * ((t.toSpqrTree.neRange i).1 + (j + 1)) + r) = r := by
    unfold QE.side QE.dir; omega
  rw [he, htw, hslot, Option.getD_some, t.nodeOfNe_cap hwf hcs hcap]

/-- Endpoints of node-edge `k` of an `R` node: inside its node-vertex range, increasing. -/
theorem nvs_R (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hR : t.toSpqrTree.type i = .R)
    {k : Nat} (hk : k < t.toSpqrTree.nEdges i) :
    (t.toSpqrTree.nvRange i).1 ≤ (t.nodeEdges[(t.toSpqrTree.neRange i).1 + k]!).nvs.1 ∧
    (t.nodeEdges[(t.toSpqrTree.neRange i).1 + k]!).nvs.1 <
      (t.nodeEdges[(t.toSpqrTree.neRange i).1 + k]!).nvs.2 ∧
    (t.nodeEdges[(t.toSpqrTree.neRange i).1 + k]!).nvs.2 < (t.toSpqrTree.nvRange i).2 := by
  have hsz := PlanarRot.neRange_le_size t.toSpqrTree hwf i hi
  unfold SpqrTree.nEdges at hk
  have hnode := hwf.own.ne_node i ((t.toSpqrTree.neRange i).1 + k) hi (by omega) (by omega)
  have h1 := hwf.own.ne_nvs ((t.toSpqrTree.neRange i).1 + k) (by omega) i hnode
  have hmem : (t.nodeEdges[(t.toSpqrTree.neRange i).1 + k]!).nvs ∈ t.toSpqrTree.skeleton i := by
    have := t.skeleton_getElem? i k (by omega)
    rw [← t.nodeEdges_getElem! hwf hi (by omega)] at this
    exact List.mem_of_getElem? this
  exact ⟨h1.1, (t.shape_R hwf hi hR).2.2.2 _ hmem, h1.2.2⟩

theorem nvsOf_R (hwf : t.toSpqrTree.WF) (hi : i < t.size) {k : Nat} (hk : k < t.toSpqrTree.nEdges i) :
    t.toSpqrTree.nvsOf ((t.toSpqrTree.neRange i).1 + k) =
      some (t.nodeEdges[(t.toSpqrTree.neRange i).1 + k]!).nvs := by
  have hsz := PlanarRot.neRange_le_size t.toSpqrTree hwf i hi
  unfold SpqrTree.nEdges at hk
  unfold SpqrTree.nvsOf
  rw [Array.getElem?_eq_getElem (show (t.toSpqrTree.neRange i).1 + k < t.nodeEdges.size by omega),
    getElem!_pos t.nodeEdges ((t.toSpqrTree.neRange i).1 + k) (by omega)]
  rfl

theorem localSkeleton_length : (t.localSkeleton i).length = t.toSpqrTree.nEdges i := by
  simp [localSkeleton, PlanarRot.skeleton_length]

theorem localSkeleton_getElem? (hwf : t.toSpqrTree.WF) (hi : i < t.size) {k : Nat}
    (hk : k < t.toSpqrTree.nEdges i) :
    (t.localSkeleton i)[k]? =
      some ((t.nodeEdges[(t.toSpqrTree.neRange i).1 + k]!).nvs.1 - (t.toSpqrTree.nvRange i).1,
        (t.nodeEdges[(t.toSpqrTree.neRange i).1 + k]!).nvs.2 - (t.toSpqrTree.nvRange i).1) := by
  unfold localSkeleton
  unfold SpqrTree.nEdges at hk
  rw [List.getElem?_map, t.skeleton_getElem? i k hk, ← t.nodeEdges_getElem! hwf hi (by omega)]
  rfl

/-- `nodeRot i` reads `neRotAdj` on the node's segment, shifted by `4 * neSt`. -/
theorem nodeRot_get (hwf : t.toSpqrTree.WF) (hi : i < t.size)
    (hsz : 4 * (t.toSpqrTree.neRange i).2 ≤ t.neRotAdj.size) {l : Nat}
    (hl : l < 4 * t.toSpqrTree.nEdges i) :
    (t.nodeRot i).get l =
      (t.neRotAdj[4 * (t.toSpqrTree.neRange i).1 + l]?).bind
        (Option.map (· - 4 * (t.toSpqrTree.neRange i).1)) := by
  have hen := t.neEn_R hwf hi
  unfold nodeRot RotationSystem.get
  obtain ⟨⟨neSt, neEn⟩, hne⟩ : ∃ p, t.toSpqrTree.neRange i = p := ⟨_, rfl⟩
  rw [hne] at hsz hen ⊢
  rw [Array.getElem?_map, Array.getElem?_extract, if_pos (by omega)]
  cases t.neRotAdj[4 * neSt + l]? with
  | none => rfl
  | some o => cases o <;> rfl

theorem nodeRot_size_le (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hR : t.toSpqrTree.type i = .R)
    (hloc : IsPlanarEmbedding (t.localSkeleton i) (t.toSpqrTree.nVerts i) (t.nodeRot i)) :
    4 * (t.toSpqrTree.neRange i).2 ≤ t.neRotAdj.size := by
  have hs := hloc.size
  rw [t.localSkeleton_length] at hs
  have h6 := (t.shape_R hwf hi hR).2.1
  have hen := t.neEn_R hwf hi
  unfold nodeRot RotationSystem.size at hs
  obtain ⟨⟨neSt, neEn⟩, hne⟩ : ∃ p, t.toSpqrTree.neRange i = p := ⟨_, rfl⟩
  rw [hne] at hs hen ⊢
  simp only [Array.size_map, Array.size_extract] at hs ⊢
  omega

/-- The executable rotation entry of the node quarter-edge `4 * neSt + l` of an `R` node is the
`nodeRot` entry shifted back; it stays inside the node's segment. -/
theorem rotR (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hR : t.toSpqrTree.type i = .R)
    (hloc : IsPlanarEmbedding (t.localSkeleton i) (t.toSpqrTree.nVerts i) (t.nodeRot i))
    (hclosed : t.NodeRotClosed) {l : Nat} (hl : l < 4 * t.toSpqrTree.nEdges i) :
    t.neRotAdj[4 * (t.toSpqrTree.neRange i).1 + l]! =
      some (4 * (t.toSpqrTree.neRange i).1 + (t.nodeRot i).rot l) ∧
    (t.nodeRot i).rot l < 4 * t.toSpqrTree.nEdges i := by
  have hsz := t.nodeRot_size_le hwf hi hR hloc
  have hs := hloc.size
  rw [t.localSkeleton_length] at hs
  have hen := t.neEn_R hwf hi
  have htot := hloc.total l (by omega)
  obtain ⟨r, hr⟩ := Option.isSome_iff_exists.1 htot
  have hinv := hloc.involution l (by omega) r hr
  have hget := t.nodeRot_get hwf hi hsz hl
  rw [hr] at hget
  obtain ⟨o, ho, hmap⟩ := Option.bind_eq_some_iff.1 hget.symm
  obtain ⟨tb, rfl, htb⟩ := Option.map_eq_some_iff.1 hmap
  have hcl := hclosed i hi hR _ _ (by omega) (by omega) ho
  have hrot : (t.nodeRot i).rot l = r := by unfold RotationSystem.rot; rw [hr]; rfl
  refine ⟨?_, by omega⟩
  rw [hrot, getElem!_pos t.neRotAdj _ (by omega)]
  rw [Array.getElem?_eq_getElem (by omega)] at ho
  rw [Option.some.inj ho]
  congr 1; omega

/-- Original endpoints of node-edge `k` of an `R` node, as `sV` of its (local) node-vertices. -/
theorem neOrig_R (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hR : t.toSpqrTree.type i = .R) {k : Nat} (hk : k < t.toSpqrTree.nEdges i) :
    t.toSpqrTree.neOrig ((t.toSpqrTree.neRange i).1 + k) =
      some (t.sV i ((t.nodeEdges[(t.toSpqrTree.neRange i).1 + k]!).nvs.1 - (t.toSpqrTree.nvRange i).1),
        t.sV i ((t.nodeEdges[(t.toSpqrTree.neRange i).1 + k]!).nvs.2 - (t.toSpqrTree.nvRange i).1)) := by
  obtain ⟨p, hp⟩ : ∃ p, t.toSpqrTree.neOrig ((t.toSpqrTree.neRange i).1 + k) = some p := by
    cases k with
    | zero => exact hsep.cap_orig i _ hi (t.capNe_R hR)
    | succ j =>
      obtain ⟨hc, -, hcap, -, -, -, horig, -⟩ :=
        t.child_R (g := g) hwf hsep hi hR (j := j) (by omega) (t.sC_getElem?_R hwf hi hR (by omega))
      rw [horig]
      exact hsep.cap_orig _ _ (t.child_data hwf hi hc).1 (by unfold SpqrTree.capNe; rw [hcap]; rfl)
  have hb := t.nvs_R hwf hi hR hk
  rw [hp]
  rw [neOrig_iff] at hp
  obtain ⟨nvs, hnvs, h1, h2⟩ := hp
  rw [t.nvsOf_R hwf hi hk] at hnvs
  cases hnvs
  unfold sV
  rw [Nat.add_sub_cancel' hb.1, Nat.add_sub_cancel' (by omega), h1, h2]
  rfl

end R

end Spqr.PlanarSpqrTree
