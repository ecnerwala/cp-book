import Spqr.Proofs.PlanarGenus

/-!
# Deleting an edge of a planar embedding (`IsPlanarEmbedding.uninsert`)

The inverse of `RotationSystem.insert`: the first edge `p` of `p :: es` is deleted from a planar
embedding `rs` in which `p` is not a bridge (`EdgesConn es p.1 p.2`). Its two sides are then
distinct faces (`cap_sides_distinct`: were they cofacial, `detach_both` would give the deleted
system `F + 2` face orbits, contradicting `genus_le`), the deleted system `σ` is planar, pairs
the former neighbours `rot 0 ↔ rot 1` and `rot 2 ↔ rot 3`, agrees with `rs` elsewhere, and the
two sides merge: `rot 1 - 4` and `rot 3 - 4` are cofacial in `σ` (`sameOrbit_detach_merge`).
-/

namespace Spqr

open RotationSystem

/-- **Edge deletion.** Deleting the first edge `p` (distinct endpoints, not a bridge) of a planar
embedding `rs` of `p :: es`: its sides `2`, `0` are distinct faces, all four neighbours lie above
the edge, and the shortened system `σ` (quarter-edge `4 + q` of `rs` is `q` of `σ`) is a planar
embedding of `es` which pairs `rot 0 ↔ rot 1`, `rot 2 ↔ rot 3`, agrees with `rs` elsewhere, and
has `rot 1 - 4`, `rot 3 - 4` cofacial. -/
theorem IsPlanarEmbedding.uninsert {es : List (Nat × Nat)} {p : Nat × Nat} {n : Nat}
    {rs : RotationSystem} (h : IsPlanarEmbedding (p :: es) n rs) (hne : p.1 ≠ p.2)
    (hconn : EdgesConn es p.1 p.2) :
    ¬SameOrbit (rs.stepC 3) 2 0 ∧
    4 ≤ rs.rot 0 ∧ 4 ≤ rs.rot 1 ∧ 4 ≤ rs.rot 2 ∧ 4 ≤ rs.rot 3 ∧
    ∃ σ : RotationSystem, IsPlanarEmbedding es n σ ∧ σ.size + 4 = rs.size ∧
      (∀ q, q < σ.size → σ.get q = some
        (if 4 + q = rs.rot 0 then rs.rot 1 - 4 else if 4 + q = rs.rot 1 then rs.rot 0 - 4 else
         if 4 + q = rs.rot 2 then rs.rot 3 - 4 else if 4 + q = rs.rot 3 then rs.rot 2 - 4 else
         rs.rot (4 + q) - 4)) ∧
      SameOrbit (σ.stepC 3) (rs.rot 1 - 4) (rs.rot 3 - 4) := by
  have hE := h.toIsEmbedding
  have hm : rs.size = 4 * (es.length + 1) := by rw [hE.size]; rfl
  have ht := hE.total
  have hi := hE.involution
  have ho := hE.opposite_dir
  have hsv := hE.same_vertex
  have hverts : ∀ q ∈ es, q.1 < n ∧ q.2 < n := fun q hq => hE.verts q (List.mem_cons_of_mem _ hq)
  have hp : p.1 < n ∧ p.2 < n := hE.verts p (List.mem_cons_self ..)
  have hu : HasEdge es p.1 := (HasEdge.of_edgesConn hconn hne).1
  have hv : HasEdge es p.2 := (HasEdge.of_edgesConn hconn hne).2
  have hv0 : QE.vert (p :: es) 0 = some p.1 := by rw [vert_cons_lt _ _ (by decide)]; rfl
  have hv1 : QE.vert (p :: es) 1 = some p.1 := by rw [vert_cons_lt _ _ (by decide)]; rfl
  have hv2 : QE.vert (p :: es) 2 = some p.2 := by rw [vert_cons_lt _ _ (by decide)]; rfl
  have hv3 : QE.vert (p :: es) 3 = some p.2 := by rw [vert_cons_lt _ _ (by decide)]; rfl
  have hr0 := rot_lt ht hi (q := 0) (by omega)
  have hr1 := rot_lt ht hi (q := 1) (by omega)
  have hr2 := rot_lt ht hi (q := 2) (by omega)
  have hr3 := rot_lt ht hi (q := 3) (by omega)
  have hp0 := rot_mod2 ht ho (q := 0) (by omega)
  have hp1 := rot_mod2 ht ho (q := 1) (by omega)
  have hp2 := rot_mod2 ht ho (q := 2) (by omega)
  have hp3 := rot_mod2 ht ho (q := 3) (by omega)
  have hrr0 := rot_rot ht hi (q := 0) (by omega)
  have hrr1 := rot_rot ht hi (q := 1) (by omega)
  have hrr2 := rot_rot ht hi (q := 2) (by omega)
  have hrr3 := rot_rot ht hi (q := 3) (by omega)
  have hstep1 : ∀ q, q < rs.size → rs.stepC 1 q = rs.rot (q ^^^ 1) :=
    fun q hq => stepC_eq_rot ht hm (by decide) hq
  have hvr : ∀ a, a < 4 → QE.vert (p :: es) (rs.rot a) = QE.vert (p :: es) a :=
    fun a ha => vert_rot ht hsv (by omega)
  -- no neighbour of the edge lies on the edge itself
  have h0 : 4 ≤ rs.rot 0 := by
    by_contra hlt
    rcases (by omega : rs.rot 0 = 1 ∨ rs.rot 0 = 3) with e | e
    · have e' : rs.rot 1 = 0 := by rw [← e, hrr0]
      exact not_hasEdge_of_pair hE (a := 0) (b := 0) (by omega) (by omega)
        (by rw [hstep1 0 (by omega)]; exact e') (by rw [hstep1 0 (by omega)]; exact e') hv0 hu
    · have := hvr 0 (by omega); rw [e, hv3, hv0] at this; exact hne (Option.some_inj.1 this).symm
  have h1 : 4 ≤ rs.rot 1 := by
    by_contra hlt
    rcases (by omega : rs.rot 1 = 0 ∨ rs.rot 1 = 2) with e | e
    · have e' : rs.rot 0 = 1 := by rw [← e, hrr1]
      omega
    · have := hvr 1 (by omega); rw [e, hv2, hv1] at this; exact hne (Option.some_inj.1 this).symm
  have h2 : 4 ≤ rs.rot 2 := by
    by_contra hlt
    rcases (by omega : rs.rot 2 = 1 ∨ rs.rot 2 = 3) with e | e
    · have := hvr 2 (by omega); rw [e, hv1, hv2] at this; exact hne (Option.some_inj.1 this)
    · have e' : rs.rot 3 = 2 := by rw [← e, hrr2]
      exact not_hasEdge_of_pair hE (a := 2) (b := 2) (by omega) (by omega)
        (by rw [hstep1 2 (by omega)]; exact e') (by rw [hstep1 2 (by omega)]; exact e') hv2 hv
  have h3 : 4 ≤ rs.rot 3 := by
    by_contra hlt
    rcases (by omega : rs.rot 3 = 0 ∨ rs.rot 3 = 2) with e | e
    · have := hvr 3 (by omega); rw [e, hv0, hv3] at this; exact hne (Option.some_inj.1 this)
    · have e' : rs.rot 2 = 3 := by rw [← e, hrr3]
      omega
  obtain ⟨H, hHt, hHi, hHo, hHsv, hHs, hH0, hH1, hH2, hH3, hHV, hHF, hHrot⟩ :=
    detach_both hE h0 h1 h2 h3
  -- `H` is the disjoint union of the isolated edge and the shortened system
  have hhigh : ∀ q, 4 ≤ q → q < H.size → 4 ≤ H.rot q := by
    intro q h4 hq
    rw [hHrot q (by omega) h4]
    have hrq := rot_rot ht hi (q := q) (by omega)
    split_ifs <;> try omega
    by_contra hlt
    interval_cases hv : rs.rot q <;> omega
  have hH4 : 4 ≤ H.size := by omega
  have hlow : ∀ q, q < 4 → H.get q = single.get q := by
    intro q hq
    obtain ⟨s0, s1, s2, s3⟩ := single_rot4
    rw [get_eq_rot hHt (by omega), get_eq_rot single_total (by rw [single_size]; omega)]
    interval_cases q <;> simp only [hH0, hH1, hH2, hH3, s0, s1, s2, s3]
  have hHeq : H = single.union H.shrink := eq_union_shrink single_size single_total hHt hH4 hlow hhigh
  set σ := H.shrink with hσ
  have hσs : σ.size = 4 * es.length := by rw [hσ, shrink_size, hHs, hm]; omega
  have hσt : σ.Total := shrink_total hHt
  have hσi : σ.Involution := shrink_involution hHt hHi hhigh
  have hVu := union_numVertexOrbits single σ single_total single_involution hσt hσi (m₁ := 1)
    single_size hσs
  have hFu := union_numFaceOrbits single σ single_total single_involution hσt hσi (m₁ := 1)
    single_size hσs
  rw [← hHeq, single_numVertexOrbits] at hVu
  rw [← hHeq, single_numFaceOrbits] at hFu
  have hNI : numNonIsolated (p :: es) n = numNonIsolated es n := by
    rw [numNonIsolated_cons_eq hp, if_pos hu, if_pos (Or.inr hv)]; omega
  have hσE : IsEmbedding es n σ :=
    isEmbedding_of_union single_size (by rw [← hHeq, hHs, hE.size]) hverts (by rw [← hHeq]; exact hHt)
      (by rw [← hHeq]; exact hHi) (by rw [← hHeq]; exact hHo) (by rw [← hHeq]; exact hHsv)
      (by have := hE.vertex_orbits; omega)
  have hgen := genus_le es hσE
  have heul := h.euler
  unfold EulerFormula at heul
  rw [List.length_cons] at heul
  have hC := numComponents_cons_le_of_conn hverts hp hconn
  rcases hHF with hF | ⟨hF, hside, hcof⟩
  · exfalso; omega
  refine ⟨hside, h0, h1, h2, h3, σ, ⟨hσE, ?_⟩, by omega, ?_, ?_⟩
  · unfold EulerFormula; omega
  · intro q hq
    rw [hσ, shrink_get, get_eq_rot hHt (by omega), Option.map_some, hHrot (4 + q) (by omega) (by omega)]
    congr 1
    split_ifs <;> rfl
  · have hstep : H.stepC 3 = unionStep (single.stepC 3) (shift single.size (σ.stepC 3))
        (Finset.range single.size) := by
      rw [hHeq]; exact union_stepC single σ (m₁ := 1) single_size (by decide)
    rw [hstep] at hcof
    have hdisj : Disjoint (Finset.range single.size) ((Finset.range σ.size).image (single.size + ·)) := by
      rw [Finset.disjoint_left]
      intro q hq hq'
      rw [Finset.mem_range] at hq
      obtain ⟨a, -, rfl⟩ := Finset.mem_image.1 hq'
      omega
    have hB : rs.rot 1 ∈ (Finset.range σ.size).image (single.size + ·) :=
      Finset.mem_image.2 ⟨rs.rot 1 - 4, Finset.mem_range.2 (by omega), by rw [single_size]; omega⟩
    have := sameOrbit_shift_sub ((sameOrbit_unionStep_right
      (isPermOn_stepC single_total single_involution (m := 1) single_size (by decide))
      (isPermOn_shift (isPermOn_stepC hσt hσi hσs (by decide)) _) hdisj hB).1 hcof)
    rwa [single_size] at this

/-- A non-bridge edge of a planar embedding has its two sides in distinct faces. -/
theorem IsPlanarEmbedding.cap_sides_distinct {es : List (Nat × Nat)} {p : Nat × Nat} {n : Nat}
    {rs : RotationSystem} (h : IsPlanarEmbedding (p :: es) n rs) (hne : p.1 ≠ p.2)
    (hconn : EdgesConn es p.1 p.2) : ¬rs.SameFaceOrbit 0 2 := by
  intro hs
  have hm : rs.size = 4 * (es.length + 1) := by rw [h.size]; rfl
  exact (h.uninsert hne hconn).1
    ((sameFaceOrbit_iff h.total h.involution hm (by omega)).1 hs).symm


end Spqr
