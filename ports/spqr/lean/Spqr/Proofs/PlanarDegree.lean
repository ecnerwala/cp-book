import Spqr.Proofs.PlanarInsert
import Spqr.Proofs.PlanarMap
import Spqr.Proofs.TwoSumCount
import Mathlib.Data.Finset.Prod

/-!
# Vertex orbits of an embedding

In an `IsEmbedding` the vertex permutation has exactly two orbits per non-isolated vertex (one per
direction), so two quarter-edges at the same vertex with the same direction lie on the same vertex
orbit. Hence at an end of a non-loop edge whose vertex has another incident edge, the rotation
leaves the edge (`IsEmbedding.rot_div_ne`, the `deg` hypotheses of `TwoSum.WF`).
-/

namespace Spqr

open Classical

theorem sameOrbit_eq_of_fixed {f : Nat → Nat} {S : Finset Nat} (hf : IsPermOn f S) {q r : Nat}
    (hq : q ∈ S) (hfix : f q = q) (h : SameOrbit f q r) : r = q := by
  obtain ⟨k, hk⟩ := IsPermOn.reach_iff.1 (hf.reach_of_sameOrbit hq h)
  rw [← hk]
  clear hk h
  induction k with
  | zero => rfl
  | succ k ih => rw [Function.iterate_succ_apply', ih, hfix]

theorem xor_one_eq (q : Nat) : q ^^^ 1 = if q % 2 = 0 then q + 1 else q - 1 := by
  have hx : ∀ k, k < 4 → k ^^^ 1 = if k % 2 = 0 then k + 1 else k - 1 := by decide
  rw [xor_eq4 q (by decide), hx _ (Nat.mod_lt q (by omega))]
  split_ifs <;> omega

theorem vert_xor_one (es : List (Nat × Nat)) (q : Nat) : QE.vert es (q ^^^ 1) = QE.vert es q := by
  unfold QE.vert QE.edge QE.side
  rw [xor_one_eq]
  have h1 : (if q % 2 = 0 then q + 1 else q - 1) / 4 = q / 4 := by split_ifs <;> omega
  have h2 : (if q % 2 = 0 then q + 1 else q - 1) / 2 % 2 = q / 2 % 2 := by split_ifs <;> omega
  rw [h1, h2]

/-- The vertex-and-direction label of a quarter-edge. -/
def lbl (es : List (Nat × Nat)) (q : Nat) : Nat × Nat := ((QE.vert es q).getD 0, q % 2)

namespace IsEmbedding

variable {es : List (Nat × Nat)} {n : Nat} {rs : RotationSystem}

theorem lbl_stepC (h : IsEmbedding es n rs) (q : Nat) : lbl es (rs.stepC 1 q) = lbl es q := by
  by_cases hq : q < rs.size
  · have hx : q ^^^ 1 < rs.size := h.size ▸ xor_lt_mul4 (h.size ▸ hq) (by decide)
    have h1 : QE.vert es (rs.stepC 1 q) = QE.vert es q := by
      rw [RotationSystem.stepC_eq_rot h.total h.size (by decide) hq,
        RotationSystem.vert_rot h.total h.same_vertex hx, vert_xor_one]
    have h2 := RotationSystem.stepC_mod2 (c := 1) h.total h.opposite_dir h.size (by decide) (by decide) q
    simp only [lbl, h1, h2]
  · have : rs.stepC 1 q = q :=
      (RotationSystem.isPermOn_stepC h.total h.involution h.size (by decide)).fix q
        (by simpa using hq)
    rw [this]

/-- Two quarter-edges at the same vertex in the same direction lie on one vertex orbit. -/
theorem sameOrbit_of_lbl (h : IsEmbedding es n rs) {q r : Nat} (hq : q < rs.size)
    (hr : r < rs.size) (hl : lbl es q = lbl es r) : SameOrbit (rs.stepC 1) q r := by
  set f := rs.stepC 1 with hfdef
  set S := Finset.range rs.size with hSdef
  have hf : IsPermOn f S := RotationSystem.isPermOn_stepC h.total h.involution h.size (by decide)
  have hinv : ∀ a b, SameOrbit f a b → lbl es a = lbl es b := fun a b hab =>
    sameOrbit_invariant (lbl es) h.lbl_stepC hab
  have hcount : orbitCount f S = 2 * numNonIsolated es n := by
    rw [← h.vertex_orbits, RotationSystem.numVertexOrbits_eq rs h.total h.involution h.size,
      ← h.size]
  set V := (Finset.range n).filter (HasEdge es) with hVdef
  set T := V ×ˢ ({0, 1} : Finset Nat) with hTdef
  have hlT : ∀ a ∈ S, lbl es a ∈ T := by
    intro a ha
    rw [hSdef, Finset.mem_range] at ha
    obtain ⟨v, hv⟩ : ∃ v, QE.vert es a = some v := by
      unfold QE.vert QE.edge
      have : a / 4 < es.length := by have := h.size; omega
      rw [List.getElem?_eq_getElem this]; exact ⟨_, rfl⟩
    rw [hTdef, Finset.mem_product]
    refine ⟨?_, ?_⟩
    · rw [hVdef, Finset.mem_filter, Finset.mem_range]
      obtain ⟨p, hp, hpv⟩ := hasEdge_of_vert hv
      have := h.verts p hp
      simp only [lbl, hv, Option.getD_some]
      exact ⟨by omega, hasEdge_of_vert hv⟩
    · simp only [lbl, Finset.mem_insert, Finset.mem_singleton]; omega
  have hTl : ∀ b ∈ T, ∃ a ∈ S, lbl es a = b := by
    intro b hb
    rw [hTdef, Finset.mem_product, hVdef, Finset.mem_filter, Finset.mem_range] at hb
    obtain ⟨⟨hv, p, hp, hpv⟩, hd⟩ := hb
    simp only [Finset.mem_insert, Finset.mem_singleton] at hd
    obtain ⟨k, hk, hpk⟩ := List.getElem_of_mem hp
    by_cases h1 : p.1 = b.1
    · refine ⟨4 * k + b.2, ?_, ?_⟩
      · rw [hSdef, Finset.mem_range, h.size]; omega
      · unfold lbl QE.vert QE.edge QE.side
        have e1 : (4 * k + b.2) / 4 = k := by omega
        have e2 : (4 * k + b.2) / 2 % 2 = 0 := by omega
        have e3 : (4 * k + b.2) % 2 = b.2 := by omega
        rw [e1, e2, e3, List.getElem?_eq_getElem hk, hpk]
        simp only [Option.map_some, Option.getD_some, ↓reduceIte]
        exact Prod.ext h1 rfl
    · have h2 : p.2 = b.1 := by
        rcases hpv with h' | h'
        · exact absurd h' h1
        · exact h'
      refine ⟨4 * k + 2 + b.2, ?_, ?_⟩
      · rw [hSdef, Finset.mem_range, h.size]; omega
      · unfold lbl QE.vert QE.edge QE.side
        have e1 : (4 * k + 2 + b.2) / 4 = k := by omega
        have e2 : (4 * k + 2 + b.2) / 2 % 2 = 1 := by omega
        have e3 : (4 * k + 2 + b.2) % 2 = b.2 := by omega
        rw [e1, e2, e3, List.getElem?_eq_getElem hk, hpk]
        simp only [Option.map_some, Option.getD_some, Nat.one_ne_zero, ↓reduceIte]
        exact Prod.ext h2 rfl
  have hcardT : T.card = 2 * numNonIsolated es n := by
    rw [hTdef, Finset.card_product, Finset.card_pair (by decide), numNonIsolated_eq_card, hVdef]
    omega
  set O := S.image (orbit f S) with hOdef
  have hcardO : O.card = T.card := by
    rw [hcardT, ← hcount]; rfl
  have hrep : ∀ a (ha : a ∈ O), ∀ q ∈ S, orbit f S q = a →
      lbl es (Classical.choose (Finset.mem_image.1 ha)) = lbl es q := by
    intro a ha q hq hqa
    obtain ⟨hc1, hc2⟩ := Classical.choose_spec (Finset.mem_image.1 ha)
    generalize Classical.choose (Finset.mem_image.1 ha) = c at hc1 hc2 ⊢
    have : c ∈ orbit f S q := by
      rw [hqa, ← hc2]; exact mem_orbit.2 ⟨hc1, SameOrbit.refl f _⟩
    exact (hinv _ _ (mem_orbit.1 this).2).symm
  have hF : ∀ a (ha : a ∈ O), lbl es (Classical.choose (Finset.mem_image.1 ha)) ∈ T :=
    fun a ha => hlT _ (Classical.choose_spec (Finset.mem_image.1 ha)).1
  have hsurj : ∀ b ∈ T, ∃ (a : Finset Nat) (ha : a ∈ O), lbl es (Classical.choose (Finset.mem_image.1 ha)) = b := by
    intro b hb
    obtain ⟨a, ha, rfl⟩ := hTl b hb
    have hO : orbit f S a ∈ O := Finset.mem_image_of_mem (orbit f S) ha
    exact ⟨orbit f S a, hO, hrep _ hO a ha rfl⟩
  have hinj := Finset.inj_on_of_surj_on_of_card_le
    (fun a ha => lbl es (Classical.choose (Finset.mem_image.1 ha))) hF hsurj hcardO.le
  have hqS : q ∈ S := Finset.mem_range.2 hq
  have hrS : r ∈ S := Finset.mem_range.2 hr
  have hOq : orbit f S q ∈ O := Finset.mem_image_of_mem (orbit f S) hqS
  have hOr : orbit f S r ∈ O := Finset.mem_image_of_mem (orbit f S) hrS
  have heq : orbit f S q = orbit f S r :=
    hinj hOq hOr (by rw [hrep _ hOq q hqS rfl, hrep _ hOr r hrS rfl, hl])
  have : r ∈ orbit f S q := by rw [heq]; exact mem_orbit.2 ⟨hrS, SameOrbit.refl f _⟩
  exact (mem_orbit.1 this).2

/-- At an end of a non-loop edge whose vertex carries another edge, the rotation leaves the edge. -/
theorem rot_div_ne (h : IsEmbedding es n rs) {q u v k : Nat} (hq : q < rs.size)
    (he : es[q / 4]? = some (u, v)) (huv : u ≠ v) (hk : k ≠ q / 4)
    (hinc : ∃ p, es[k]? = some p ∧ (QE.vert es q = some p.1 ∨ QE.vert es q = some p.2)) :
    rs.rot q / 4 ≠ q / 4 := by
  intro hdiv
  have hrq : rs.rot q < rs.size := RotationSystem.rot_lt h.total h.involution hq
  have hvq : QE.vert es q = some (if q / 2 % 2 = 0 then u else v) := by
    unfold QE.vert QE.edge QE.side; rw [he]; rfl
  have hvr : QE.vert es (rs.rot q) = some (if rs.rot q / 2 % 2 = 0 then u else v) := by
    unfold QE.vert QE.edge QE.side; rw [hdiv, he]; rfl
  have hvv := RotationSystem.vert_rot h.total h.same_vertex hq
  rw [hvr, hvq] at hvv
  have hside : rs.rot q / 2 % 2 = q / 2 % 2 := by
    by_contra hne
    apply huv
    have := Option.some.inj hvv
    split_ifs at this <;> omega
  have hdir := RotationSystem.rot_mod2 h.total h.opposite_dir hq
  have hrot : rs.rot q = q ^^^ 1 := by rw [xor_one_eq]; split_ifs <;> omega
  have hfix : rs.stepC 1 q = q := by
    rw [RotationSystem.stepC_eq_rot h.total h.size (by decide) hq, ← hrot,
      RotationSystem.rot_rot h.total h.involution hq]
  obtain ⟨p, hp, hpv⟩ := hinc
  have hkl : k < es.length := (List.getElem?_eq_some_iff.1 hp).1
  obtain ⟨r, hr, hrk, hlr⟩ : ∃ r, r < rs.size ∧ r / 4 = k ∧ lbl es r = lbl es q := by
    rcases hpv with hpv | hpv
    · refine ⟨4 * k + q % 2, by rw [h.size]; omega, by omega, ?_⟩
      unfold lbl QE.vert QE.edge QE.side
      have e1 : (4 * k + q % 2) / 4 = k := by omega
      have e2 : (4 * k + q % 2) / 2 % 2 = 0 := by omega
      have e3 : (4 * k + q % 2) % 2 = q % 2 := by omega
      rw [e1, e2, e3, hp]
      unfold QE.vert QE.edge QE.side at hpv
      rw [hpv]; rfl
    · refine ⟨4 * k + 2 + q % 2, by rw [h.size]; omega, by omega, ?_⟩
      unfold lbl QE.vert QE.edge QE.side
      have e1 : (4 * k + 2 + q % 2) / 4 = k := by omega
      have e2 : (4 * k + 2 + q % 2) / 2 % 2 = 1 := by omega
      have e3 : (4 * k + 2 + q % 2) % 2 = q % 2 := by omega
      rw [e1, e2, e3, hp]
      unfold QE.vert QE.edge QE.side at hpv
      rw [hpv]; rfl
  have hso := h.sameOrbit_of_lbl hq hr hlr.symm
  have := sameOrbit_eq_of_fixed
    (RotationSystem.isPermOn_stepC h.total h.involution h.size (by decide))
    (Finset.mem_range.2 hq) hfix hso
  omega

end IsEmbedding

end Spqr
