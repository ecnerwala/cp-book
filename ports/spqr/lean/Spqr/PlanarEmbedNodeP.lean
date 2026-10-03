import Spqr.PlanarEmbedNodeSMain

/-!
# The `P` node: data layer, bond fold and `nodeFold_capped_P`

A `P` node has node-vertices `s, s + 1`, `k ≥ 3` node-edges all with `nvs = (s, s + 1)`, the cap
first and then the caps' twins of its `k - 1` children in order (`Twins.noncap_children`), and no
`V` child (both node-vertices are cap endpoints, `PieceSep.nv_child`). Its `neRotAdj` segment is
`layoutRot .P`, i.e. `rotP`: block `0` (the cap) exposes slots `0, 3` of the first child and
`1, 2` of the last one; block `j` (`1 ≤ j ≤ k - 2`) links child `j - 1` to child `j` by
`c1 ↔ d0`, `c2 ↔ d3` (`Capped.parJoin`); block `k - 1` is a no-op. The induction invariant `InvP`
carries the `Capped` certificate of the chain of the first `m + 1` children at the cap's original
endpoints; the chain of all children is `pieceBelow g i` itself.
-/

namespace Spqr.PlanarSpqrTree

variable (t : PlanarSpqrTree)

/-- The `neRotAdj` segment of node `i` is the `layoutRot` of its type (`neRotAdj_segment` for
`g.planarTree`). -/
abbrev LayoutAt (i : Nat) : Prop :=
  ∃ ev mr cv, t.neRotAdj.extract (4 * (t.toSpqrTree.neRange i).1) (4 * (t.toSpqrTree.neRange i).2) =
    layoutRot (t.toSpqrTree.type i) (t.toSpqrTree.nVerts i) (t.toSpqrTree.neRange i).1
      (t.toSpqrTree.neRange i).2 ev mr cv

section P

variable {i : Nat}

theorem shape_P (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hP : t.toSpqrTree.type i = .P) :
    t.toSpqrTree.nVerts i = 2 ∧ 3 ≤ t.toSpqrTree.nEdges i ∧
    ∀ p ∈ t.toSpqrTree.skeleton i, p = ((t.toSpqrTree.nvRange i).1, (t.toSpqrTree.nvRange i).1 + 1) := by
  have hsh := hwf.shape i hi
  unfold SpqrTree.Shape at hsh
  rw [hP] at hsh
  obtain ⟨⟨s, e⟩, hnv⟩ : ∃ p, t.toSpqrTree.nvRange i = p := ⟨_, rfl⟩
  rw [hnv] at hsh
  simp only at hsh
  unfold SpqrTree.nVerts
  rw [hnv]
  exact hsh

theorem hasCap_P (hP : t.toSpqrTree.type i = .P) : t.toSpqrTree.hasCap i = true := by
  simp [SpqrTree.hasCap, hP]; decide

theorem capNe_P (hP : t.toSpqrTree.type i = .P) :
    t.toSpqrTree.capNe i = some (t.toSpqrTree.neRange i).1 := by
  simp [SpqrTree.capNe, t.hasCap_P hP]

/-- Every node-edge of a `P` node has endpoints `(s, s + 1)`. -/
theorem nvs_P (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hP : t.toSpqrTree.type i = .P)
    {k : Nat} (hk : k < t.toSpqrTree.nEdges i) :
    (t.nodeEdges[(t.toSpqrTree.neRange i).1 + k]!).nvs =
      ((t.toSpqrTree.nvRange i).1, (t.toSpqrTree.nvRange i).1 + 1) := by
  have h := t.skeleton_getElem? i k (by unfold SpqrTree.nEdges at hk; omega)
  rw [t.nodeEdges_getElem! hwf hi (by unfold SpqrTree.nEdges at hk; omega)]
  exact (t.shape_P hwf hi hP).2.2 _ (List.mem_of_getElem? h)

theorem nvsOf_P (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hP : t.toSpqrTree.type i = .P)
    {k : Nat} (hk : k < t.toSpqrTree.nEdges i) :
    t.toSpqrTree.nvsOf ((t.toSpqrTree.neRange i).1 + k) =
      some ((t.toSpqrTree.nvRange i).1, (t.toSpqrTree.nvRange i).1 + 1) := by
  have hsz := PlanarRot.neRange_le_size t.toSpqrTree hwf i hi
  unfold SpqrTree.nEdges at hk
  rw [← t.nvs_P hwf hi hP (by unfold SpqrTree.nEdges; omega)]
  unfold SpqrTree.nvsOf
  rw [Array.getElem?_eq_getElem (show (t.toSpqrTree.neRange i).1 + k < t.nodeEdges.size by omega),
    getElem!_pos t.nodeEdges ((t.toSpqrTree.neRange i).1 + k) (by omega)]
  rfl

/-- The non-cap node-edges of a `P` node are the twins of the caps of its non-`V` children, in
order. -/
theorem noncap_P (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hP : t.toSpqrTree.type i = .P) :
    (t.sC i).length = t.toSpqrTree.nEdges i - 1 ∧
    ∀ k, k < t.toSpqrTree.nEdges i - 1 →
      ∃ c, (t.sC i)[k]? = some c ∧
        (t.nodeEdges[(t.toSpqrTree.neRange i).1 + (k + 1)]!).twin = t.toSpqrTree.capNe c := by
  have hnc := hwf.twins.noncap_children i hi (by rw [hP]; decide)
  rw [t.hasCap_P hP, t.children_eq] at hnc
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

/-- The `neRotAdj` entries of a `P` node are `rotP`. -/
theorem neRotAdj_P (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hP : t.toSpqrTree.type i = .P)
    (hlay : t.LayoutAt i) {k r : Nat} (hk : k < t.toSpqrTree.nEdges i) (hr : r < 4) :
    t.neRotAdj[4 * ((t.toSpqrTree.neRange i).1 + k) + r]! =
      some (rotP (t.toSpqrTree.nEdges i) (t.toSpqrTree.neRange i).1 k r) := by
  obtain ⟨ev, mr, cv, h⟩ := hlay
  have hmono := PlanarRot.neRange_mono t.toSpqrTree hwf i hi
  have hen : (t.toSpqrTree.neRange i).2 = (t.toSpqrTree.neRange i).1 + t.toSpqrTree.nEdges i := by
    unfold SpqrTree.nEdges; omega
  rw [hP, hen, (t.shape_P hwf hi hP).1] at h
  have hg := layoutRot_P_get (t.toSpqrTree.nEdges i) (t.toSpqrTree.neRange i).1 ev mr cv k r hk hr
  rw [← h, Array.getElem?_extract] at hg
  split at hg
  · have hg' : t.neRotAdj[4 * (t.toSpqrTree.neRange i).1 + (4 * k + r)]? =
        some (some (rotP (t.toSpqrTree.nEdges i) (t.toSpqrTree.neRange i).1 k r)) := hg
    obtain ⟨hb, hv⟩ := Array.getElem?_eq_some_iff.1 hg'
    rw [show 4 * ((t.toSpqrTree.neRange i).1 + k) + r = 4 * (t.toSpqrTree.neRange i).1 + (4 * k + r) by omega,
      getElem!_pos t.neRotAdj (4 * (t.toSpqrTree.neRange i).1 + (4 * k + r)) hb, hv]
  · cases hg

theorem rotP_P (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hP : t.toSpqrTree.type i = .P)
    (hlay : t.LayoutAt i) {k r : Nat} (hk : k < t.toSpqrTree.nEdges i) (hr : r < 4) {ta tb : Nat}
    (hta : ta = 4 * ((t.toSpqrTree.neRange i).1 + k) + r)
    (htb : tb = rotP (t.toSpqrTree.nEdges i) (t.toSpqrTree.neRange i).1 k r) :
    t.neRotAdj[ta]! = some tb := by
  subst hta htb; exact t.neRotAdj_P hwf hi hP hlay hk hr

theorem sC_length_P (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hP : t.toSpqrTree.type i = .P) :
    (t.sC i).length = t.toSpqrTree.nEdges i - 1 := (t.noncap_P hwf hi hP).1

theorem sC_getElem?_P (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hP : t.toSpqrTree.type i = .P)
    {j : Nat} (hj : j < t.toSpqrTree.nEdges i - 1) : (t.sC i)[j]? = some (t.sC i)[j]! := by
  have := t.sC_length_P hwf hi hP
  rw [getElem!_pos (t.sC i) j (by omega), List.getElem?_eq_getElem (by omega)]

/-- Facts about the `j`-th non-`V` child `c` of a `P` node (cf. `child_S`). -/
theorem child_P (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hP : t.toSpqrTree.type i = .P) {j c : Nat} (hj : j < t.toSpqrTree.nEdges i - 1)
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
  obtain ⟨hcap, hcI, hcO⟩ := hsep.node_child_cap i c hi (Or.inr (Or.inl hP)) (by rwa [t.children_eq]) hcV
  obtain ⟨c', hc', htw⟩ := (t.noncap_P hwf hi hP).2 j hj
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

theorem nvOrig_P (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hP : t.toSpqrTree.type i = .P) {j : Nat} (hj : j < 2) :
    t.toSpqrTree.nvOrig ((t.toSpqrTree.nvRange i).1 + j) = some (t.sV i j) := by
  have h3 := (t.shape_P hwf hi hP).2.1
  obtain ⟨hc, -, hcap, -, -, -, horig, -⟩ :=
    t.child_P (g := g) hwf hsep hi hP (j := 0) (by omega) (t.sC_getElem?_P hwf hi hP (by omega))
  have hcne : t.toSpqrTree.capNe (t.sC i)[0]! = some (t.toSpqrTree.neRange (t.sC i)[0]!).1 := by
    simp only [SpqrTree.capNe, hcap, ↓reduceIte]
  obtain ⟨p, hp⟩ := hsep.cap_orig _ _ (t.child_data hwf hi hc).1 hcne
  have hp' : t.toSpqrTree.neOrig ((t.toSpqrTree.neRange i).1 + 1) = some p := horig.trans hp
  rw [neOrig_iff] at hp'
  obtain ⟨nvs, hnvs, h1, h2⟩ := hp'
  rw [t.nvsOf_P hwf hi hP (k := 1) (by omega)] at hnvs
  cases hnvs
  simp only at h1 h2
  unfold sV
  interval_cases j
  · rw [Nat.add_zero, h1]; rfl
  · rw [h2]; rfl

/-- Every non-cap node-edge of a `P` node has original endpoints `sV 0`, `sV 1`. -/
theorem neOrig_P (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hP : t.toSpqrTree.type i = .P) {k : Nat} (hk : k < t.toSpqrTree.nEdges i) :
    t.toSpqrTree.neOrig ((t.toSpqrTree.neRange i).1 + k) = some (t.sV i 0, t.sV i 1) := by
  rw [neOrig_iff]
  refine ⟨_, t.nvsOf_P hwf hi hP hk, ?_, ?_⟩
  · have := t.nvOrig_P (g := g) hwf hsep hi hP (j := 0) (by omega)
    rwa [Nat.add_zero] at this
  · exact t.nvOrig_P (g := g) hwf hsep hi hP (j := 1) (by omega)

theorem capOrig_P (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hP : t.toSpqrTree.type i = .P) {p : Nat × Nat}
    (hp : t.toSpqrTree.neOrig (t.toSpqrTree.neRange i).1 = some p) : p = (t.sV i 0, t.sV i 1) := by
  have h := t.neOrig_P (g := g) hwf hsep hi hP (k := 0) (by have := (t.shape_P hwf hi hP).2.1; omega)
  rw [Nat.add_zero, hp] at h
  exact Option.some.inj h

theorem sV01_ne_P (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hP : t.toSpqrTree.type i = .P) : t.sV i 0 ≠ t.sV i 1 := by
  intro e
  have h1 := t.nvOrig_P (g := g) hwf hsep hi hP (j := 0) (by omega)
  have h2 := t.nvOrig_P (g := g) hwf hsep hi hP (j := 1) (by omega)
  rw [e, ← h2] at h1
  have hn := (t.shape_P hwf hi hP).1
  unfold SpqrTree.nVerts at hn
  have := hsep.nv_orig_inj i _ _ hi (Or.inr (Or.inl hP)) ⟨by omega, by omega⟩ ⟨by omega, by omega⟩ h1
  omega

/-- The `j`-th child of a `P` node carries a `Capped` certificate at `sV 0`, `sV 1`. -/
theorem child_cert_P (hwf : t.toSpqrTree.WF) (hne : t.ne = g.ne) (hsep : t.toSpqrTree.PieceSep g)
    (hi : i < t.size) (hP : t.toSpqrTree.type i = .P) {j : Nat} (hj : j < t.toSpqrTree.nEdges i - 1)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    (t.sC i)[j]! ∈ t.children i ∧
    ∃ a0 a1 a2 a3 ρ, s.outerE[(t.sC i)[j]!]![0]! = some a0 ∧ s.outerE[(t.sC i)[j]!]![1]! = some a1 ∧
      s.outerE[(t.sC i)[j]!]![2]! = some a2 ∧ s.outerE[(t.sC i)[j]!]![3]! = some a3 ∧
      (t.pieceBelow g (t.sC i)[j]!).Capped s.rotAdj ρ a0 a1 a2 a3 (t.sV i 0) (t.sV i 1) := by
  obtain ⟨hc, -, hcap, hcI, hcO, -, horig, -⟩ :=
    t.child_P (g := g) hwf hsep hi hP hj (t.sC_getElem?_P hwf hi hP hj)
  refine ⟨hc, ?_⟩
  have hp : t.toSpqrTree.neOrig (t.toSpqrTree.neRange (t.sC i)[j]!).1 = some (t.sV i 0, t.sV i 1) := by
    rw [← horig]
    exact t.neOrig_P (g := g) hwf hsep hi hP (by unfold SpqrTree.nEdges at hj ⊢; omega)
  exact t.capped_child g hwf hne hsep hi hc hcI hcO hcap hp s h

/-- Both node-vertices of a `P` node are cap endpoints, so it has no `V` child. -/
theorem noV_P (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hP : t.toSpqrTree.type i = .P) {c : Nat} (hc : c ∈ t.children i) : t.toSpqrTree.type c ≠ .V := by
  intro hV
  obtain ⟨nv, d, hnv, hd, hdv⟩ :=
    hsep.v_child_nv i c hi (Or.inr (Or.inl hP)) (by rw [t.children_eq]; exact hc) hV
  have hpar := (t.child_data hwf hi hc).2.1
  have hnc : ¬ t.toSpqrTree.CapEnd i nv :=
    (hsep.nv_child i nv d hi (Or.inr (Or.inl hP)) hnv hd).1 (by rw [hdv]; exact hpar)
  apply hnc
  have hn := (t.shape_P hwf hi hP).1
  unfold SpqrTree.nVerts at hn
  have hsz : (t.toSpqrTree.neRange i).1 < t.nodeEdges.size := by
    have := PlanarRot.neRange_le_size t.toSpqrTree hwf i hi
    have := (t.shape_P hwf hi hP).2.1
    unfold SpqrTree.nEdges at this
    omega
  have h0 := t.nvs_P hwf hi hP (k := 0) (by have := (t.shape_P hwf hi hP).2.1; omega)
  rw [Nat.add_zero, getElem!_pos t.nodeEdges _ hsz] at h0
  refine ⟨_, _, t.capNe_P hP, Array.getElem?_eq_getElem hsz, ?_⟩
  rw [h0]
  obtain ⟨h1, h2⟩ := hnv
  simp only
  omega

theorem sC_eq_children_P (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hP : t.toSpqrTree.type i = .P) : t.sC i = t.children i := by
  unfold sC
  rw [List.filter_eq_self]
  intro c hc
  simpa using t.noV_P (g := g) hwf hsep hi hP hc

/-- Two distinct children of a `P` node meet only at `sV 0`, `sV 1`. -/
theorem attach_P (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hP : t.toSpqrTree.type i = .P) {a b x : Nat} (ha : a ∈ t.children i) (hb : b ∈ t.children i)
    (hab : a ≠ b) (hta : t.toSpqrTree.Touches g a x) (htb : t.toSpqrTree.Touches g b x) :
    x = t.sV i 0 ∨ x = t.sV i 1 := by
  rw [← t.children_eq] at ha hb
  obtain ⟨nv, hnv, ho⟩ := hsep.node_attach i a b x hi (Or.inr (Or.inl hP)) ha hb hab hta htb
  have hn := (t.shape_P hwf hi hP).1
  unfold SpqrTree.nVerts at hn
  obtain ⟨h1, h2⟩ := hnv
  obtain ⟨j, rfl⟩ : ∃ j, nv = (t.toSpqrTree.nvRange i).1 + j := ⟨nv - (t.toSpqrTree.nvRange i).1, by omega⟩
  have hj : j < 2 := by omega
  rw [t.nvOrig_P (g := g) hwf hsep hi hP (j := j) hj] at ho
  cases ho
  interval_cases j
  · exact Or.inl rfl
  · exact Or.inr rfl

end P

end Spqr.PlanarSpqrTree
