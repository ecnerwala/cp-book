import Spqr.Proofs.PlanarMap
import Spqr.Proofs.PlanarGlue

namespace Spqr

open Classical

theorem graphCounts_congr {es es' : List (Nat × Nat)} {n : Nat}
    (hm : ∀ p, p ∈ es ↔ p ∈ es') (hv : ∀ p ∈ es, p.1 < n ∧ p.2 < n) :
    numNonIsolated es n = numNonIsolated es' n ∧ numComponents es n = numComponents es' n := by
  have he : HasEdge es = HasEdge es' := by
    funext v
    exact propext (by simp only [HasEdge, hm])
  have hc : EdgesConn es = EdgesConn es' := by
    funext a b
    unfold EdgesConn
    congr 1
    funext u v
    exact propext (by simp only [hm])
  refine ⟨?_, ?_⟩
  · rw [numNonIsolated_eq_card, numNonIsolated_eq_card, he]
  · rw [numComponents_eq_ccCount hv,
      numComponents_eq_ccCount (fun p hp => hv p ((hm p).2 hp))]
    simp only [ccCount, he, hc]

theorem IsPlanarEmbedding.reindex {es es' : List (Nat × Nat)} {n : Nat}
    {ρ σ : RotationSystem} (h : IsPlanarEmbedding es n ρ) (φ : Nat → Nat)
    (hsize : σ.size = 4 * es'.length) (hlen : es.length = es'.length)
    (hmem : ∀ p, p ∈ es ↔ p ∈ es')
    (hmap : ∀ q, q < σ.size → φ q < ρ.size)
    (hinj : ∀ q r, q < σ.size → r < σ.size → φ q = φ r → q = r)
    (hget : ∀ q, q < σ.size → ρ.get (φ q) = (σ.get q).map φ)
    (hbound : ∀ q, q < σ.size → ∀ r ∈ σ.get q, r < σ.size)
    (hdir : ∀ q, q < σ.size → QE.dir (φ q) = QE.dir q)
    (hvert : ∀ q, q < σ.size → QE.vert es (φ q) = QE.vert es' q)
    (hxor : ∀ q, q < σ.size → ∀ c, c < 4 → φ (q ^^^ c) = φ q ^^^ c) :
    IsPlanarEmbedding es' n σ := by
  have ht : σ.Total := by
    intro q hq
    have hh := h.total (φ q) (hmap q hq)
    simpa only [hget q hq, Option.isSome_map] using hh
  have hi : σ.Involution := by
    intro q hq r hr
    have hrb := hbound q hq r hr
    refine ⟨hrb, ?_⟩
    have hρ := (h.involution (φ q) (hmap q hq) (φ r)
      (by rw [hget q hq, hr]; rfl)).2
    rw [hget r hrb] at hρ
    obtain ⟨q', hq', he⟩ := Option.map_eq_some_iff.1 hρ
    have : q' = q := hinj q' q (hbound r hrb q' hq') hq he
    simpa [this] using hq'
  have himg : (Finset.range σ.size).image φ = Finset.range ρ.size := by
    apply Finset.eq_of_subset_of_card_le
    · intro r hr
      obtain ⟨q, hq, rfl⟩ := Finset.mem_image.1 hr
      exact Finset.mem_range.2 (hmap q (Finset.mem_range.1 hq))
    · rw [Finset.card_image_of_injOn (fun q hq r hr =>
        hinj q r (Finset.mem_range.1 hq) (Finset.mem_range.1 hr)), Finset.card_range,
        Finset.card_range, h.size, hsize, hlen]
  have hc : ∀ c, c < 4 → orbitCount (ρ.stepC c) (Finset.range ρ.size) =
      orbitCount (σ.stepC c) (Finset.range σ.size) := by
    intro c hc
    apply orbitCount_congr (σ.isPermOn_stepC ht hi hsize hc)
      (ρ.isPermOn_stepC h.total h.involution h.size hc)
      (fun q hq r hr => hinj q r (Finset.mem_range.1 hq) (Finset.mem_range.1 hr)) himg
    intro q hq
    have hqb := Finset.mem_range.1 hq
    have hqc : q ^^^ c < σ.size := by
      rw [hsize] at hqb ⊢
      exact xor_lt_mul4 hqb hc
    rw [RotationSystem.stepC_apply, RotationSystem.stepC_apply, ← hxor q hqb c hc, hget _ hqc]
    obtain ⟨r, hr⟩ := Option.isSome_iff_exists.1 (ht _ hqc)
    simp [hr]
  have hv : σ.numVertexOrbits = ρ.numVertexOrbits := by
    rw [RotationSystem.numVertexOrbits_eq _ ht hi hsize,
      RotationSystem.numVertexOrbits_eq _ h.total h.involution h.size]
    simpa only [hsize, h.size] using (hc 1 (by decide)).symm
  have hf : σ.numFaceOrbits = ρ.numFaceOrbits := by
    rw [RotationSystem.numFaceOrbits_eq _ ht hi hsize,
      RotationSystem.numFaceOrbits_eq _ h.total h.involution h.size]
    simpa only [hsize, h.size] using (hc 3 (by decide)).symm
  obtain ⟨hn, hcomp⟩ := graphCounts_congr hmem h.verts
  refine ⟨⟨hsize, fun p hp => h.verts p ((hmem p).2 hp), ht, hi, ?_, ?_, ?_⟩, ?_⟩
  · intro q hq r hr
    have hh := h.opposite_dir (φ q) (hmap q hq) (φ r) (by rw [hget q hq, hr]; rfl)
    simpa only [hdir q hq, hdir r (hbound q hq r hr)] using hh
  · intro q hq r hr
    have hh := h.same_vertex (φ q) (hmap q hq) (φ r) (by rw [hget q hq, hr]; rfl)
    simpa only [hvert q hq, hvert r (hbound q hq r hr)] using hh
  · rw [hv, h.vertex_orbits, hn]
  · simpa only [EulerFormula, hf, ← hn, ← hcomp, ← hlen] using h.euler

end Spqr
