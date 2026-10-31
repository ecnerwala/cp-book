import Spqr.Proofs.OrbitSplit

/-!
# Orbit segments

Small facts about iterates of a permutation used to track where a quarter-edge lands after a
face-splitting conjugation: iterates agree with a modified step as long as the path avoids the
modified points, the first-hit path to a point never returns to its start, and face orbits are
reversed by `rot` (`stepC 3 ∘ rot ∘ stepC 3 = rot`).
-/

namespace Spqr

theorem iterate_eq_of_eqOn {f g : Nat → Nat} {a : Nat} {k : Nat}
    (h : ∀ j, j < k → g (f^[j] a) = f (f^[j] a)) : g^[k] a = f^[k] a := by
  induction k generalizing a with
  | zero => rfl
  | succ k ih =>
    have h0 := h 0 (Nat.succ_pos k)
    simp only [Function.iterate_zero, id_eq] at h0
    rw [Function.iterate_succ_apply, Function.iterate_succ_apply, h0]
    exact ih fun j hj => by
      have := h (j + 1) (by omega)
      rwa [Function.iterate_succ_apply] at this

theorem IsPermOn.exists_first_hit {f : Nat → Nat} {S : Finset Nat} (hf : IsPermOn f S) {a b : Nat}
    (ha : a ∈ S) (h : SameOrbit f a b) :
    ∃ k, f^[k] a = b ∧ ∀ j, j < k → f^[j] a ≠ b ∧ f^[j + 1] a ≠ a := by
  have hex : ∃ k, f^[k] a = b := IsPermOn.reach_iff.1 (hf.reach_of_sameOrbit ha h)
  refine ⟨Nat.find hex, Nat.find_spec hex, fun j hj => ⟨Nat.find_min hex hj, fun hp => ?_⟩⟩
  apply Nat.find_min hex (show Nat.find hex - (j + 1) < Nat.find hex by omega)
  have e : f^[Nat.find hex - (j + 1)] (f^[j + 1] a) = b := by
    rw [← Function.iterate_add_apply, Nat.sub_add_cancel hj]
    exact Nat.find_spec hex
  rwa [hp] at e

namespace RotationSystem

variable {rs : RotationSystem}

theorem stepC_lt_iff (ht : rs.Total) (hi : rs.Involution) {m c : Nat} (hs : rs.size = 4 * m)
    (hc : c < 4) (q : Nat) : rs.stepC c q < rs.size ↔ q < rs.size := by
  have hf := isPermOn_stepC ht hi hs hc
  constructor
  · intro h; exact Finset.mem_range.1 (hf.mem_of_apply_mem (Finset.mem_range.2 h))
  · intro h; exact Finset.mem_range.1 (hf.maps q (Finset.mem_range.2 h))

theorem sameOrbit_stepC_lt_iff (ht : rs.Total) (hi : rs.Involution) {m c : Nat}
    (hs : rs.size = 4 * m) (hc : c < 4) {a b : Nat} (h : SameOrbit (rs.stepC c) a b) :
    a < rs.size ↔ b < rs.size := by
  have := sameOrbit_invariant (fun q => decide (q < rs.size))
    (fun q => by simp only [decide_eq_decide]; exact stepC_lt_iff ht hi hs hc q) h
  simpa using this

theorem stepC3_rot_stepC3 (ht : rs.Total) (hi : rs.Involution) {m : Nat} (hs : rs.size = 4 * m)
    {q : Nat} (hq : q < rs.size) : rs.stepC 3 (rs.rot (rs.stepC 3 q)) = rs.rot q := by
  have h1 : q ^^^ 3 < rs.size := hs ▸ xor_lt_mul4 (hs ▸ hq) (by decide)
  rw [stepC_eq_rot ht hs (by decide) hq, rot_rot ht hi h1, stepC_eq_rot ht hs (by decide) h1,
    xor_xor_self]

theorem sameOrbit_stepC3_rot (ht : rs.Total) (hi : rs.Involution) {m : Nat} (hs : rs.size = 4 * m)
    {a b : Nat} (ha : a < rs.size) (h : SameOrbit (rs.stepC 3) a b) :
    SameOrbit (rs.stepC 3) (rs.rot a) (rs.rot b) := by
  induction h with
  | rel a b e =>
    subst e
    have := SameOrbit.step (f := rs.stepC 3) (rs.rot (rs.stepC 3 a))
    rw [stepC3_rot_stepC3 ht hi hs ha] at this
    exact this.symm
  | refl a => exact SameOrbit.refl _ _
  | symm a b h ih =>
    exact (ih ((sameOrbit_stepC_lt_iff ht hi hs (by decide) h).2 ha)).symm
  | trans a b c h₁ h₂ ih₁ ih₂ =>
    exact (ih₁ ha).trans (ih₂ ((sameOrbit_stepC_lt_iff ht hi hs (by decide) h₁).1 ha))

end RotationSystem

end Spqr
