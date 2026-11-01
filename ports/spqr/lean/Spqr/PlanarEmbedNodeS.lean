import Spqr.PlanarEmbedNode
import Spqr.PlanarLayout
import Spqr.PlanarRotSpec

/-!
# The `S` node: skeleton, twins and rotation entries

Bookkeeping for `nodeFold_capped_S`: the node-edges of an `S` node `i` with node-vertices
`s, …, s + n - 1` are the cap `(s, s + n - 1)` followed by the path edges `(s + k - 1, s + k)`,
the twin of path edge `k` is the cap of the `k`-th non-`V` child, and the node's `neRotAdj`
segment is `rotS` (`layoutRot .S`).
-/

namespace Spqr.PlanarSpqrTree

variable (t : PlanarSpqrTree)

theorem mem_children_of_parent (hwf : t.toSpqrTree.WF) {i j : Nat} (hi : i < t.size)
    (hj : j < t.size) (hp : t.toSpqrTree.parent j = some i) : j ∈ t.children i := by
  rw [← t.children_eq, hwf.preorder.ch_eq i hi, List.mem_filter]
  exact ⟨List.mem_range.2 hj, by simp [hp]⟩

theorem nodeEdges_getElem! (hwf : t.toSpqrTree.WF) {i k : Nat} (hi : i < t.size)
    (hk : (t.toSpqrTree.neRange i).1 + k < (t.toSpqrTree.neRange i).2) :
    t.nodeEdges[(t.toSpqrTree.neRange i).1 + k]! =
      (t.nodeEdges[(t.toSpqrTree.neRange i).1 + k]?).getD default := by
  have := PlanarRot.neRange_le_size t.toSpqrTree hwf i hi
  rw [getElem!_pos _ _ (by omega), Array.getElem?_eq_getElem (by omega)]; rfl

theorem skeleton_getElem? (i k : Nat)
    (hk : k < (t.toSpqrTree.neRange i).2 - (t.toSpqrTree.neRange i).1) :
    (t.toSpqrTree.skeleton i)[k]? =
      some (t.nodeEdges[(t.toSpqrTree.neRange i).1 + k]?.getD default).nvs := by
  simp [SpqrTree.skeleton, SpqrTree.nodeEdgesOf, List.getElem?_range, hk]

theorem nodeEdgesOf_getElem? (i k : Nat)
    (hk : k < (t.toSpqrTree.neRange i).2 - (t.toSpqrTree.neRange i).1) :
    (t.toSpqrTree.nodeEdgesOf i)[k]? =
      some (t.nodeEdges[(t.toSpqrTree.neRange i).1 + k]?.getD default) := by
  simp [SpqrTree.nodeEdgesOf, List.getElem?_range, hk]

theorem nodeEdgesOf_length (i : Nat) :
    (t.toSpqrTree.nodeEdgesOf i).length = (t.toSpqrTree.neRange i).2 - (t.toSpqrTree.neRange i).1 := by
  simp [SpqrTree.nodeEdgesOf]

section S

variable {i : Nat}

theorem shape_S (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hS : t.toSpqrTree.type i = .S) :
    3 ≤ t.toSpqrTree.nVerts i ∧
    t.toSpqrTree.skeleton i =
      ((t.toSpqrTree.nvRange i).1, (t.toSpqrTree.nvRange i).1 + t.toSpqrTree.nVerts i - 1) ::
        (List.range (t.toSpqrTree.nVerts i - 1)).map fun k =>
          ((t.toSpqrTree.nvRange i).1 + k, (t.toSpqrTree.nvRange i).1 + k + 1) := by
  have hsh := hwf.shape i hi
  unfold SpqrTree.Shape at hsh
  rw [hS] at hsh
  obtain ⟨⟨s, e⟩, hnv⟩ : ∃ p, t.toSpqrTree.nvRange i = p := ⟨_, rfl⟩
  rw [hnv] at hsh
  simp only at hsh
  unfold SpqrTree.nVerts
  rw [hnv]
  simp only
  refine ⟨hsh.1, ?_⟩
  rw [hsh.2]
  congr 2
  omega

theorem nEdges_S (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hS : t.toSpqrTree.type i = .S) :
    (t.toSpqrTree.neRange i).2 - (t.toSpqrTree.neRange i).1 = t.toSpqrTree.nVerts i := by
  have hl := PlanarRot.skeleton_length t.toSpqrTree i
  rw [(t.shape_S hwf hi hS).2] at hl
  simp only [List.length_cons, List.length_map, List.length_range] at hl
  have := (t.shape_S hwf hi hS).1
  unfold SpqrTree.nEdges at hl
  omega

theorem hasCap_S (hS : t.toSpqrTree.type i = .S) : t.toSpqrTree.hasCap i = true := by
  simp [SpqrTree.hasCap, hS]; decide

theorem capNe_S (hS : t.toSpqrTree.type i = .S) :
    t.toSpqrTree.capNe i = some (t.toSpqrTree.neRange i).1 := by
  simp [SpqrTree.capNe, t.hasCap_S hS]

/-- Endpoints of the node-edges of an `S` node: the cap closes the path `s, s + 1, …, s + n - 1`. -/
theorem nvs_S (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hS : t.toSpqrTree.type i = .S)
    {k : Nat} (hk : k < t.toSpqrTree.nVerts i) :
    (t.nodeEdges[(t.toSpqrTree.neRange i).1 + k]!).nvs =
      if k = 0 then ((t.toSpqrTree.nvRange i).1, (t.toSpqrTree.nvRange i).1 + t.toSpqrTree.nVerts i - 1)
      else ((t.toSpqrTree.nvRange i).1 + k - 1, (t.toSpqrTree.nvRange i).1 + k) := by
  have hn := t.nEdges_S hwf hi hS
  have h := t.skeleton_getElem? i k (by omega)
  rw [t.nodeEdges_getElem! hwf hi (by omega), (t.shape_S hwf hi hS).2] at *
  cases k with
  | zero =>
    simp only [List.getElem?_cons_zero, Option.some.injEq, Nat.add_zero] at h
    simp [← h]
  | succ k =>
    simp only [List.getElem?_cons_succ, List.getElem?_map,
      List.getElem?_range (show k < t.toSpqrTree.nVerts i - 1 by omega),
      Option.map_some, Option.some.injEq] at h
    rw [if_neg (Nat.succ_ne_zero k), ← h]
    congr 1 <;> omega

/-- The non-cap node-edges of an `S` node are the twins of the caps of its non-`V` children, in
order. -/
theorem noncap_S (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hS : t.toSpqrTree.type i = .S) :
    ((t.children i).filter fun c => t.toSpqrTree.type c ≠ .V).length = t.toSpqrTree.nVerts i - 1 ∧
    ∀ k, k < t.toSpqrTree.nVerts i - 1 →
      ∃ c, ((t.children i).filter fun c => t.toSpqrTree.type c ≠ .V)[k]? = some c ∧
        (t.nodeEdges[(t.toSpqrTree.neRange i).1 + (k + 1)]!).twin = t.toSpqrTree.capNe c := by
  have hn := t.nEdges_S hwf hi hS
  have hnc := hwf.twins.noncap_children i hi (by rw [hS]; decide)
  rw [t.hasCap_S hS, t.children_eq] at hnc
  simp only [↓reduceIte] at hnc
  have hlen := congrArg List.length hnc
  simp only [List.length_map, List.length_drop, t.nodeEdgesOf_length] at hlen
  refine ⟨by omega, ?_⟩
  intro k hk
  have hget := congrArg (fun l => l[k]?) hnc
  simp only [List.getElem?_map, List.getElem?_drop, t.nodeEdgesOf_getElem? i (1 + k) (by omega),
    Option.map_some] at hget
  obtain ⟨c, hc, hce⟩ := Option.map_eq_some_iff.1 hget.symm
  refine ⟨c, hc, ?_⟩
  rw [t.nodeEdges_getElem! hwf hi (by omega), show k + 1 = 1 + k by omega, hce]

/-- The `neRotAdj` entries of an `S` node are `rotS`. -/
theorem neRotAdj_S (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hS : t.toSpqrTree.type i = .S)
    (hlay : ∃ ev mr cv, t.neRotAdj.extract (4 * (t.toSpqrTree.neRange i).1) (4 * (t.toSpqrTree.neRange i).2) =
      layoutRot (t.toSpqrTree.type i) (t.toSpqrTree.nVerts i) (t.toSpqrTree.neRange i).1
        (t.toSpqrTree.neRange i).2 ev mr cv)
    {k r : Nat} (hk : k < t.toSpqrTree.nVerts i) (hr : r < 4) :
    t.neRotAdj[4 * ((t.toSpqrTree.neRange i).1 + k) + r]! =
      some (rotS (t.toSpqrTree.nVerts i) (t.toSpqrTree.neRange i).1 k r) := by
  obtain ⟨ev, mr, cv, h⟩ := hlay
  have hn := t.nEdges_S hwf hi hS
  have h3 := (t.shape_S hwf hi hS).1
  have hmono := PlanarRot.neRange_mono t.toSpqrTree hwf i hi
  have hen : (t.toSpqrTree.neRange i).2 = (t.toSpqrTree.neRange i).1 + t.toSpqrTree.nVerts i := by omega
  rw [hS, hen] at h
  have hg := layoutRot_S_get (t.toSpqrTree.nVerts i) (t.toSpqrTree.neRange i).1 (by omega) ev mr cv k r hk hr
  rw [← h, Array.getElem?_extract] at hg
  split at hg
  · have hg' : t.neRotAdj[4 * (t.toSpqrTree.neRange i).1 + (4 * k + r)]? =
        some (some (rotS (t.toSpqrTree.nVerts i) (t.toSpqrTree.neRange i).1 k r)) := hg
    obtain ⟨hb, hv⟩ := Array.getElem?_eq_some_iff.1 hg'
    rw [show 4 * ((t.toSpqrTree.neRange i).1 + k) + r = 4 * (t.toSpqrTree.neRange i).1 + (4 * k + r) by omega,
      getElem!_pos t.neRotAdj (4 * (t.toSpqrTree.neRange i).1 + (4 * k + r)) hb, hv]
  · cases hg

end S

end Spqr.PlanarSpqrTree

namespace Spqr.PlanarSpqrTree

variable (t : PlanarSpqrTree)

theorem getElem!_of_getElem? {α : Type} [Inhabited α] {a : Array α} {i : Nat} {x : α}
    (h : a[i]? = some x) : a[i]! = x := by
  obtain ⟨hb, hv⟩ := Array.getElem?_eq_some_iff.1 h
  rw [getElem!_pos a i hb, hv]

/-! ### Unfolding `nodeStep` -/

theorem nodeStep_skip (i neSt : Nat) (s : EmbedState) {ta tb : Nat} (h : t.neRotAdj[ta]! = some tb)
    (hlt : tb < ta) : t.nodeStep i neSt s ta = s := by
  unfold nodeStep; rw [h]; simp [hlt]

theorem nodeStep_cap (i neSt : Nat) (s : EmbedState) {ta tb : Nat} (h : t.neRotAdj[ta]! = some tb)
    (hlt : ¬ tb < ta) (hta : ta < 4 * (neSt + 1)) :
    t.nodeStep i neSt s ta =
      ((EmbedM.setOuter i (2 * QE.side ta + (1 - QE.dir ta)) (t.treeQe s tb)).run s).2 := by
  unfold nodeStep; rw [h]; simp [hlt, hta]

theorem nodeStep_link (i neSt : Nat) (s : EmbedState) {ta tb : Nat} (h : t.neRotAdj[ta]! = some tb)
    (hlt : ¬ tb < ta) (hta : ¬ ta < 4 * (neSt + 1)) (hc : ¬ (ta % 4 = 2 ∧ tb % 4 = 1)) :
    t.nodeStep i neSt s ta = ((EmbedM.link (t.treeQe s ta) (t.treeQe s tb)).run s).2 := by
  unfold nodeStep; rw [h]
  have : (ta % 4 == 2 && tb % 4 == 1) = false := by
    simp only [Bool.and_eq_false_imp, beq_iff_eq, beq_eq_false_iff_ne]
    intro h2; exact fun h1 => hc ⟨h2, h1⟩
  simp [hlt, hta, this]

theorem nodeStep_corner (i neSt : Nat) (s : EmbedState) {ta tb : Nat} (h : t.neRotAdj[ta]! = some tb)
    (hlt : ¬ tb < ta) (hta : ¬ ta < 4 * (neSt + 1)) (h2 : ta % 4 = 2) (h1 : tb % 4 = 1) :
    t.nodeStep i neSt s ta =
      let v := t.nodeVerts[(t.nodeEdges[QE.edge ta]!).nvs.2]!.vert
      if s.outerE[v]![0]!.isSome then
        ((EmbedM.link (t.treeQe s tb) s.outerE[v]![0]!).run
          ((EmbedM.link (t.treeQe s ta) s.outerE[v]![1]!).run s).2).2
      else ((EmbedM.link (t.treeQe s ta) (t.treeQe s tb)).run s).2 := by
  unfold nodeStep; rw [h]
  simp [hlt, hta, h2, h1]

/-! ### Children of an `S` node -/

theorem nEdges_pos_of_hasCap (hwf : t.toSpqrTree.WF) {c : Nat} (hc : c < t.size)
    (hcap : t.toSpqrTree.hasCap c = true) :
    (t.toSpqrTree.neRange c).1 < (t.toSpqrTree.neRange c).2 := by
  have hsh := hwf.shape c hc
  unfold SpqrTree.Shape at hsh
  have hl := PlanarRot.skeleton_length t.toSpqrTree c
  unfold SpqrTree.nEdges at hl hsh
  obtain ⟨⟨a, b⟩, hnv⟩ : ∃ p, t.toSpqrTree.nvRange c = p := ⟨_, rfl⟩
  rw [hnv] at hsh
  simp only at hsh
  have hnode : (t.toSpqrTree.type c).isNode = true := by
    unfold SpqrTree.hasCap at hcap
    exact (Bool.and_eq_true_iff.1 hcap).1
  revert hsh hnode
  cases t.toSpqrTree.type c <;> intro hsh hnode <;> simp [NodeType.isNode] at hnode
  all_goals first
    | (rcases hsh with ⟨_, h⟩ | ⟨_, h⟩ <;> rw [h] at hl <;> simp at hl <;> omega)
    | (rw [hsh.2] at hl; simp at hl; omega)
    | (obtain ⟨_, h2, _⟩ := hsh; omega)

theorem nodeOfNe_cap (hwf : t.toSpqrTree.WF) {c : Nat} (hc : c < t.size)
    (hcap : t.toSpqrTree.hasCap c = true) :
    (t.nodeEdges[(t.toSpqrTree.neRange c).1]!).node = c := by
  have h := hwf.own.ne_node c (t.toSpqrTree.neRange c).1 hc (Nat.le_refl _)
    (t.nEdges_pos_of_hasCap hwf hc hcap)
  unfold SpqrTree.nodeOfNe at h
  obtain ⟨d, hd, hdn⟩ := Option.map_eq_some_iff.1 h
  rw [getElem!_of_getElem? hd, hdn]

section S

variable {i : Nat}

/-- Facts about the `j`-th non-`V` child `c` of an `S` node: it is a capped non-`I`/`O` child, the
twin of path edge `j + 1` is its cap, and `treeQe` on that path edge reads its exposed row. -/
theorem child_S (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hS : t.toSpqrTree.type i = .S) {j c : Nat} (hj : j < t.toSpqrTree.nVerts i - 1)
    (hc : ((t.children i).filter fun c => t.toSpqrTree.type c ≠ .V)[j]? = some c) :
    c ∈ t.children i ∧ t.toSpqrTree.type c ≠ .V ∧ t.toSpqrTree.hasCap c = true ∧
    t.toSpqrTree.type c ≠ .I ∧ t.toSpqrTree.type c ≠ .O ∧
    t.toSpqrTree.twin ((t.toSpqrTree.neRange i).1 + (j + 1)) = some (t.toSpqrTree.neRange c).1 ∧
    t.toSpqrTree.neOrig ((t.toSpqrTree.neRange i).1 + (j + 1)) =
      t.toSpqrTree.neOrig (t.toSpqrTree.neRange c).1 ∧
    ∀ r, r < 4 → ∀ s : EmbedState,
      t.treeQe s (4 * ((t.toSpqrTree.neRange i).1 + (j + 1)) + r) = s.outerE[c]![r]! := by
  have hmem := List.mem_of_getElem? hc
  rw [List.mem_filter] at hmem
  have hcm : c ∈ t.children i := hmem.1
  have hcV : t.toSpqrTree.type c ≠ .V := by simpa using hmem.2
  obtain ⟨hcap, hcI, hcO⟩ := hsep.node_child_cap i c hi (Or.inl hS) (by rwa [t.children_eq]) hcV
  have hn := t.nEdges_S hwf hi hS
  obtain ⟨c', hc', htw⟩ := (t.noncap_S hwf hi hS).2 j hj
  rw [hc] at hc'; cases hc'
  have hcne : t.toSpqrTree.capNe c = some (t.toSpqrTree.neRange c).1 := by
    simp [SpqrTree.capNe, hcap]
  rw [hcne] at htw
  have hbnd : (t.toSpqrTree.neRange i).1 + (j + 1) < (t.toSpqrTree.neRange i).2 := by omega
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

end S

end Spqr.PlanarSpqrTree
