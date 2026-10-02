import Spqr.PlanarEmbedSteps

/-!
# The leaf step of `planarEmbed`, under the hypotheses it needs

`embedItem_step_leaf` (`PlanarEmbedSteps.lean`) is stated for an arbitrary `PlanarSpqrTree` and
is not provable as stated: nothing says the leaf has no items below it (`subtreeEnd`) nor that
`outerE[i]` is still unset when item `i` is processed (`GluedUpTo` has no clause about unprocessed
items). `embedItem_step_leaf_of` proves it from exactly those two facts; the first is
`WF.preorder.subtree_eq` for a childless item, the second is the missing `GluedUpTo` clause.
-/

namespace Spqr

/-- The empty rotation system is a planar embedding of the empty graph. -/
theorem isPlanarEmbedding_nil (n : Nat) : IsPlanarEmbedding [] n ⟨#[]⟩ where
  size := rfl
  verts := by simp
  total := fun q hq => absurd hq (by simp [RotationSystem.size])
  involution := fun q hq => absurd hq (by simp [RotationSystem.size])
  opposite_dir := fun q hq => absurd hq (by simp [RotationSystem.size])
  same_vertex := fun q hq => absurd hq (by simp [RotationSystem.size])
  vertex_orbits := by
    simp [RotationSystem.numVertexOrbits, RotationSystem.size, numOrbits, numNonIsolated, nonIsolated]
  euler := by
    simp [EulerFormula, RotationSystem.numFaceOrbits, RotationSystem.size, numOrbits, numNonIsolated,
      numComponents, nonIsolated]

namespace PlanarSpqrTree

variable (t : PlanarSpqrTree)

/-- An item whose subtree is itself alone and which is not a `Q` item has no edges below it. -/
theorem edgesBelow_leaf (i : Nat) (hsub : t.subtreeEnd[i]! = i + 1) (hq : t.types[i]! ≠ .Q) :
    t.edgesBelow i = [] := by
  unfold edgesBelow
  rw [hsub, Nat.add_sub_cancel_left, List.range'_one, List.filterMap_cons, List.filterMap_nil]
  simp [hq]

theorem maximal_succ (i j : Nat) (hij : i + 1 ≤ j) (h : t.Maximal i j) : t.Maximal (i + 1) j :=
  ⟨hij, h.2.1, fun p hp => Nat.lt_succ_of_lt (h.2.2 p hp)⟩

/-- `embedItem` is the identity on `O`/`I` items. -/
theorem embedItem_leaf (i : Nat) (hty : t.types[i]! = .O ∨ t.types[i]! = .I) (s : EmbedState) :
    ((t.embedItem i).run s).2 = s := by
  rcases hty with hty | hty <;> (unfold embedItem; rw [hty]; rfl)

/-- The `O`/`I` step of `planarEmbed`, from `GluedUpTo (i + 1)` plus the two facts the admitted
`embedItem_step_leaf` lacks: the leaf has nothing below it, and its exposed ends are still unset. -/
theorem embedItem_step_leaf_of (g : Graph) (i : Nat) (hi : i < t.size)
    (hty : t.types[i]! = .O ∨ t.types[i]! = .I)
    (hsub : t.subtreeEnd[i]! = i + 1)
    (s : EmbedState) (houter : ∀ q, ¬ s.exposedAt i q) (h : t.GluedUpTo g (i + 1) s) :
    t.GluedUpTo g i ((t.embedItem i).run s).2 := by
  rw [t.embedItem_leaf i hty s]
  have hq : t.types[i]! ≠ .Q := by rcases hty with hty | hty <;> rw [hty] <;> decide
  have hbelow : t.edgesBelow i = [] := t.edgesBelow_leaf i hsub hq
  refine ⟨h.rot_size, h.outer_size, ?_, ?_⟩
  · intro q hq'
    exact h.unset q fun j hj hjs => hq' j (by omega) hjs
  · intro j hj
    by_cases hji : i + 1 ≤ j
    · exact h.piece j (t.maximal_succ i j hji hj)
    · have hji' : j = i := by have := hj.1; omega
      subst hji'
      have hes : (t.pieceBelow g j).es = [] := by simp [pieceBelow, Piece.es, hbelow]
      have hmem : ∀ q, ¬ (t.pieceBelow g j).Mem q := by
        intro q hm; simp [pieceBelow, Piece.Mem, hbelow] at hm
      refine ⟨⟨#[]⟩, ?_, ?_, ?_, ?_⟩
      · rw [hes]; exact isPlanarEmbedding_nil _
      · intro q _ hm; exact absurd hm (hmem q)
      · intro q hm; exact absurd hm (hmem q)
      · intro q _ hq; exact absurd hq (houter q)

end PlanarSpqrTree

end Spqr
