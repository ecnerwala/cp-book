import Spqr.PlanarEmbedNodeSMain

/-!
# The node fold: `nodeFold_capped` and the `S`/`P`/`R` step

`nodeFold_capped` dispatches on the node type: `S` is `nodeFold_capped_S` (proved,
`PlanarEmbedNodeSMain.lean`); `P` and `R` are the named admissions `nodeFold_capped_P` /
`nodeFold_capped_R` below. `embedItem_step_node_faces` combines it with `embedItem_node` and
`node_capped_glued`.
-/

namespace Spqr.PlanarSpqrTree

variable (t : PlanarSpqrTree)

/-- The `neRotAdj` segment of node `i` is the `layoutRot` of its type (`neRotAdj_segment` for
`g.planarTree`). -/
abbrev LayoutAt (i : Nat) : Prop :=
  ∃ ev mr cv, t.neRotAdj.extract (4 * (t.toSpqrTree.neRange i).1) (4 * (t.toSpqrTree.neRange i).2) =
    layoutRot (t.toSpqrTree.type i) (t.toSpqrTree.nVerts i) (t.toSpqrTree.neRange i).1
      (t.toSpqrTree.neRange i).2 ev mr cv

/-- `nodeFold_capped` for `P` nodes. Admitted: the bond layout `layoutRot .P` links consecutive
children by `c1 ↔ d0`, `c2 ↔ d3` (`Capped.parJoin`, proved) in `nodeRot` order; what is missing is
the executable correspondence (`nodeStep` over the bond segment with `treeQe`/`Piece.loc`
bookkeeping, the analogue of `PlanarEmbedNodeS*.lean` for `S`). The statement is the `P` case of the
previously admitted `nodeFold_capped`; its executable consequences (planarity of every piece and
cofaciality of the exposed slots `0`/`2` of every maximal item after each step) are checked by
`check_piece_sep` (`checkOuter`/`checkCapFace`, seeds 0..300 + tiny cases) and `compare_planar_lean.sh`. -/
theorem nodeFold_capped_P (g : Graph) (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hrep : t.toSpqrTree.Represents g) (hsep : t.toSpqrTree.PieceSep g) {i : Nat} (hi : i < t.size)
    (hP : t.toSpqrTree.type i = .P)
    (hloc : IsPlanarEmbedding (t.localSkeleton i) (t.toSpqrTree.nVerts i) (t.nodeRot i))
    (hlay : t.LayoutAt i)
    {ne : Nat} (hcne : t.toSpqrTree.capNe i = some ne) {p : Nat × Nat}
    (hp : t.toSpqrTree.neOrig ne = some p)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    let s' := (List.range' (4 * t.neBounds[i]!) (4 * t.neBounds[i + 1]! - 4 * t.neBounds[i]!)).foldl
      (t.nodeStep i t.neBounds[i]!) s
    (∀ q, ¬(t.pieceBelow g i).Mem q → s'.rotAdj[q]? = s.rotAdj[q]?) ∧
    ∃ a0 a1 a2 a3 ρ, s'.outerE[i]? = some #[some a0, some a1, some a2, some a3] ∧
      (t.pieceBelow g i).Capped s'.rotAdj ρ a0 a1 a2 a3 p.1 p.2 := by
  sorry

/-- `nodeFold_capped` for `R` nodes. Admitted: the 2-sum of `hloc` (the skeleton embedding
`nodeRot i`, cap included) with each child's `Capped` certificate through its twin
(`TwoSum.splice_isPlanarEmbedding`, `RotationSystem.insert`, `IsPlanarEmbedding.uninsert` and
`cap_sides_distinct` are proved); missing are the executable correspondence of `nodeStep` over the
`mapRot` segment with the splice's links (`treeQe`/`Piece.loc` bookkeeping, the corner `V` items,
`NodeCorners`) and that the cap is not a bridge of the skeleton. The statement is the `R` case of the
previously admitted `nodeFold_capped`; its executable consequences are checked by `check_piece_sep`
(`checkOuter`/`checkCapFace`, seeds 0..300 + tiny cases) and `compare_planar_lean.sh`. -/
theorem nodeFold_capped_R (g : Graph) (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hrep : t.toSpqrTree.Represents g) (hsep : t.toSpqrTree.PieceSep g) {i : Nat} (hi : i < t.size)
    (hR : t.toSpqrTree.type i = .R)
    (hloc : IsPlanarEmbedding (t.localSkeleton i) (t.toSpqrTree.nVerts i) (t.nodeRot i))
    (hlay : t.LayoutAt i)
    {ne : Nat} (hcne : t.toSpqrTree.capNe i = some ne) {p : Nat × Nat}
    (hp : t.toSpqrTree.neOrig ne = some p)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    let s' := (List.range' (4 * t.neBounds[i]!) (4 * t.neBounds[i + 1]! - 4 * t.neBounds[i]!)).foldl
      (t.nodeStep i t.neBounds[i]!) s
    (∀ q, ¬(t.pieceBelow g i).Mem q → s'.rotAdj[q]? = s.rotAdj[q]?) ∧
    ∃ a0 a1 a2 a3 ρ, s'.outerE[i]? = some #[some a0, some a1, some a2, some a3] ∧
      (t.pieceBelow g i).Capped s'.rotAdj ρ a0 a1 a2 a3 p.1 p.2 := by
  sorry

/-- Semantic content of the node fold: starting from `GluedFaces g (i + 1) s`, the fold of
`nodeStep` over the node's quarter-edges leaves `rotAdj` unchanged outside the node's piece, fills
the node's row with four exposed ends, and the piece is `Capped` there — planar witness agreeing
with `rotAdj`, the two exposed pairs facing at the cap's endpoints, slots 0/2 cofacial. `S` proved
(`nodeFold_capped_S`); `P`/`R` admitted (`nodeFold_capped_P`/`_R`); see PROOF.md §8.6. -/
theorem nodeFold_capped (g : Graph) (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hrep : t.toSpqrTree.Represents g) (hsep : t.toSpqrTree.PieceSep g) {i : Nat} (hi : i < t.size)
    (hty : t.toSpqrTree.type i = .S ∨ t.toSpqrTree.type i = .P ∨ t.toSpqrTree.type i = .R)
    (hloc : IsPlanarEmbedding (t.localSkeleton i) (t.toSpqrTree.nVerts i) (t.nodeRot i))
    (hlay : t.LayoutAt i)
    {ne : Nat} (hcne : t.toSpqrTree.capNe i = some ne) {p : Nat × Nat}
    (hp : t.toSpqrTree.neOrig ne = some p)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    let s' := (List.range' (4 * t.neBounds[i]!) (4 * t.neBounds[i + 1]! - 4 * t.neBounds[i]!)).foldl
      (t.nodeStep i t.neBounds[i]!) s
    (∀ q, ¬(t.pieceBelow g i).Mem q → s'.rotAdj[q]? = s.rotAdj[q]?) ∧
    ∃ a0 a1 a2 a3 ρ, s'.outerE[i]? = some #[some a0, some a1, some a2, some a3] ∧
      (t.pieceBelow g i).Capped s'.rotAdj ρ a0 a1 a2 a3 p.1 p.2 := by
  rcases hty with hS | hP | hR
  · exact t.nodeFold_capped_S g hwf hsh hrep hsep hi hS hlay hcne hp s h
  · exact t.nodeFold_capped_P g hwf hsh hrep hsep hi hP hloc hlay hcne hp s h
  · exact t.nodeFold_capped_R g hwf hsh hrep hsep hi hR hloc hlay hcne hp s h

/-- `S`/`P`/`R` step over `GluedFaces`: `embedItem_node` + `nodeFold_capped` + `node_capped_glued`. -/
theorem embedItem_step_node_faces (g : Graph) (hwf : t.toSpqrTree.WF)
    (hsh : t.toSpqrTree.ChildShape) (hrep : t.toSpqrTree.Represents g)
    (hsep : t.toSpqrTree.PieceSep g) (i : Nat) (hi : i < t.size)
    (hty : t.types[i]! = .S ∨ t.types[i]! = .P ∨ t.types[i]! = .R)
    (hloc : IsPlanarEmbedding (t.localSkeleton i) (t.toSpqrTree.nVerts i) (t.nodeRot i))
    (hlay : t.LayoutAt i)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    t.GluedFaces g i ((t.embedItem i).run s).2 := by
  have ht : t.toSpqrTree.type i = .S ∨ t.toSpqrTree.type i = .P ∨ t.toSpqrTree.type i = .R := by
    rw [t.type_eq_of_lt i hi]; exact hty
  have hcne : t.toSpqrTree.capNe i = some (t.toSpqrTree.neRange i).1 := by
    rcases ht with h' | h' | h' <;> simp [SpqrTree.capNe, SpqrTree.hasCap, h'] <;> decide
  obtain ⟨p, hp⟩ := hsep.cap_orig i _ hi hcne
  rw [t.embedItem_node i hty s]
  obtain ⟨hframe, a0, a1, a2, a3, ρ, hrow, hcap⟩ :=
    t.nodeFold_capped g hwf hsh hrep hsep hi ht hloc hlay hcne hp s h
  exact t.node_capped_glued g hwf hsh hrep hsep hi ht hcne hp s _ h
    (t.nodeFold_rotAdj_size _ _ _ _) (t.nodeFold_outerE_size _ _ _ _) hframe
    (fun j hji => t.nodeFold_outerE_ne _ _ _ _ j hji) hrow hcap

end Spqr.PlanarSpqrTree
