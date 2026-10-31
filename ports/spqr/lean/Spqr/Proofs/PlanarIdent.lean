import Spqr.Proofs.PlanarDegree

/-!
# Identifying two vertices of one planar embedding

`rs.conj a b` with `a` at `v` and `b` at `w` (`v ≠ w`, connected in `es`) whose two face pairs
`(a ^^^ 3, b ^^^ 3)` and `(rot a ^^^ 3, rot b ^^^ 3)` are cofacial embeds `identEdges v w es`:
two vertex orbits merge, two face orbits split off, the component count is unchanged.
This is the second link of a 2-sum (the first one is `IsPlanarEmbedding.oneSum_conj`).
-/

namespace Spqr

open Classical

theorem IsPlanarEmbedding.ident_conj {es : List (Nat × Nat)} {n : Nat} {rs : RotationSystem}
    (h : IsPlanarEmbedding es n rs) {a b v w : Nat}
    (ha : a < rs.size) (hb : b < rs.size) (hab : a % 2 = b % 2)
    (hva : QE.vert es a = some v) (hvb : QE.vert es b = some w) (hvw : v ≠ w)
    (hconn : EdgesConn es v w)
    (h1 : SameOrbit (rs.stepC 3) (a ^^^ 3) (b ^^^ 3))
    (h2 : SameOrbit (rs.stepC 3) (rs.rot a ^^^ 3) (rs.rot b ^^^ 3)) :
    IsPlanarEmbedding (identEdges v w es) n (rs.conj a b) := by
  have hs := h.size
  have hea := hasEdge_of_vert hva
  have heb := hasEdge_of_vert hvb
  have hvn : v < n := hea.lt_of h.verts
  have hwn : w < n := heb.lt_of h.verts
  have hab' : a ≠ b := fun e => hvw (Option.some.inj (hva.symm.trans (e ▸ hvb)))
  have hrota : QE.vert es (rs.rot a) = some v := by
    rw [RotationSystem.vert_rot h.total h.same_vertex ha, hva]
  have hrotb : QE.vert es (rs.rot b) = some w := by
    rw [RotationSystem.vert_rot h.total h.same_vertex hb, hvb]
  have hrab : rs.rot a ≠ b := fun e => hvw (Option.some.inj (hrota.symm.trans (e ▸ hvb)))
  have hvert : ∀ z, QE.vert es (rs.stepC 1 z) = QE.vert es z := by
    intro z
    by_cases hz : z < rs.size
    · have hx : z ^^^ 1 < rs.size := hs ▸ xor_lt_mul4 (hs ▸ hz) (by decide)
      rw [RotationSystem.stepC_eq_rot h.total hs (by decide) hz,
        RotationSystem.vert_rot h.total h.same_vertex hx, vert_xor_one]
    · rw [(RotationSystem.isPermOn_stepC h.total h.involution hs (by decide)).fix z
        (by simpa using hz)]
  have hne : ∀ z, QE.vert es z = some w → ¬QE.vert es z = some v := fun z hz e =>
    hvw (Option.some.inj (e.symm.trans hz))
  have hVO := RotationSystem.conj_numVertexOrbits rs a b h.total h.involution h.opposite_dir
    ha hb hab' hrab hs (fun z => QE.vert es z = some v) (fun z => by rw [hvert])
    (by rw [vert_xor_one, hva]) (by rw [vert_xor_one]; exact hne _ hvb)
    (by rw [vert_xor_one, hrota]) (by rw [vert_xor_one]; exact hne _ hrotb)
  have hFO := RotationSystem.conj_numFaceOrbits_split rs a b h.total h.involution h.opposite_dir
    ha hb hab' hrab hs h1 h2
  have hNI := numNonIsolated_ident (n := n) hvw hvn hwn hea heb
  have hCC := ccCount_ident_of_conn (n := n) hvw hvn hwn hea heb hconn
  have hev : ∀ p ∈ identEdges v w es, p.1 < n ∧ p.2 < n := by
    intro p hp
    obtain ⟨q, hq, rfl⟩ := mem_identEdges.1 hp
    have := h.verts q hq
    simp only [ident]
    split_ifs <;> omega
  exact {
    size := by rw [RotationSystem.conj_size, hs]; simp [identEdges]
    verts := hev
    total := RotationSystem.conj_total rs a b h.total ha hb
    involution := RotationSystem.conj_involution rs a b h.total h.involution ha hb
    opposite_dir := RotationSystem.conj_oppositeDir rs a b h.total h.opposite_dir ha hb hab
    same_vertex := RotationSystem.conj_sameVertex rs a b h.total ha hb h.same_vertex hva hvb
    vertex_orbits := by have := h.vertex_orbits; omega
    euler := by
      unfold EulerFormula
      have hE := h.euler
      unfold EulerFormula at hE
      rw [numComponents_eq_ccCount hev, numComponents_eq_ccCount h.verts] at *
      simp only [identEdges, List.length_map] at *
      omega }

end Spqr
