import Spqr.Proofs.Orbits
import Spqr.Proofs.Planar
import Mathlib.Data.Finset.Max

/-!
# `numOrbits` counts orbits

The marking loop `numOrbits step n` of `Planar.lean` (points that are minimal on their `n`-step
forward orbit) equals `orbitCount (stepFn step) (range n)` when `stepFn step` permutes `range n`:
every point of an orbit is reached within `n` steps, so the marked points are exactly the orbit
minima, one per orbit.
-/

namespace Spqr

open Classical

namespace IsPermOn

variable {f : Nat → Nat} {n : Nat}

/-- A period `0 < p ≤ n` of every point of `range n`. -/
theorem exists_period_le (hf : IsPermOn f (Finset.range n)) {q : Nat} (hq : q < n) :
    ∃ p, 0 < p ∧ p ≤ n ∧ f^[p] q = q := by
  have hq' : q ∈ Finset.range n := Finset.mem_range.2 hq
  have hc : (Finset.range n).card < (Finset.range (n + 1)).card := by simp
  obtain ⟨i, hi, j, hj, hij, h⟩ := Finset.exists_ne_map_eq_of_card_lt_of_maps_to hc
    (f := fun i => f^[i] q) (fun i _ => hf.iterate_mem i hq')
  rw [Finset.mem_range] at hi hj
  rcases Nat.lt_or_gt_of_ne hij with hlt | hlt
  · refine ⟨j - i, by omega, by omega, hf.iterate_inj i (hf.iterate_mem _ hq') hq' ?_⟩
    rw [← Function.iterate_add_apply, show i + (j - i) = j by omega]
    exact h.symm
  · refine ⟨i - j, by omega, by omega, hf.iterate_inj j (hf.iterate_mem _ hq') hq' ?_⟩
    rw [← Function.iterate_add_apply, show j + (i - j) = i by omega]
    exact h

/-- Every iterate is one of the first `n`. -/
theorem iterate_eq_iterate_lt (hf : IsPermOn f (Finset.range n)) {q : Nat} (hq : q < n)
    (k : Nat) : ∃ i, i < n ∧ f^[i] q = f^[k] q := by
  obtain ⟨p, hp0, hpn, hp⟩ := hf.exists_period_le hq
  refine ⟨k % p, by have := Nat.mod_lt k hp0; omega, ?_⟩
  conv_rhs => rw [← Nat.mod_add_div k p, Function.iterate_add_apply, Function.iterate_mul,
    Function.iterate_fixed hp]

end IsPermOn

/-- `isOrbitMin`: `q` is the least point of its orbit in `range n`. -/
theorem isOrbitMin_iff_orbit_min (step : Nat → Option Nat) (n : Nat)
    (hf : IsPermOn (stepFn step) (Finset.range n)) {q : Nat} (hq : q < n) :
    isOrbitMin step n q = true ↔ ∀ r ∈ orbit (stepFn step) (Finset.range n) q, q ≤ r := by
  rw [isOrbitMin_iff]
  constructor
  · intro H r hr
    rw [mem_orbit] at hr
    obtain ⟨k, hk⟩ :=
      IsPermOn.reach_iff.1 (hf.reach_of_sameOrbit (Finset.mem_range.2 hq) hr.2)
    obtain ⟨i, hi, hik⟩ := hf.iterate_eq_iterate_lt hq k
    rw [← hk, ← hik]
    exact H i hi
  · intro H i _
    exact H _ (mem_orbit.2 ⟨hf.iterate_mem i (Finset.mem_range.2 hq),
      IsPermOn.sameOrbit_of_reach (IsPermOn.reach_iff.2 ⟨i, rfl⟩)⟩)

/-- The marking loop `numOrbits` of `Planar.lean` counts the orbits of a step function that
permutes `range n`. -/
theorem numOrbits_eq_orbitCount (step : Nat → Option Nat) (n : Nat)
    (hf : IsPermOn (stepFn step) (Finset.range n)) :
    numOrbits step n = orbitCount (stepFn step) (Finset.range n) := by
  rw [numOrbits_eq_card step n (fun q => ∀ r ∈ orbit (stepFn step) (Finset.range n) q, q ≤ r)
    (fun q hq => isOrbitMin_iff_orbit_min step n hf hq)]
  unfold orbitCount
  have himg : (Finset.range n).image (orbit (stepFn step) (Finset.range n)) =
      ((Finset.range n).filter fun q => ∀ r ∈ orbit (stepFn step) (Finset.range n) q, q ≤ r).image
        (orbit (stepFn step) (Finset.range n)) := by
    ext O
    simp only [Finset.mem_image, Finset.mem_filter]
    constructor
    · rintro ⟨q, hq, rfl⟩
      have hne : (orbit (stepFn step) (Finset.range n) q).Nonempty := ⟨q, mem_orbit_self hq⟩
      have hm := mem_orbit.1 (Finset.min'_mem _ hne)
      refine ⟨Finset.min' _ hne, ⟨hm.1, ?_⟩, (orbit_eq_of_sameOrbit hm.2).symm⟩
      rw [orbit_eq_of_sameOrbit hm.2.symm]
      exact fun r hr => Finset.min'_le _ r hr
    · rintro ⟨q, ⟨hq, -⟩, rfl⟩
      exact ⟨q, hq, rfl⟩
  rw [himg, Finset.card_image_of_injOn]
  intro q hq r hr hqr
  rw [Finset.mem_coe, Finset.mem_filter] at hq hr
  have h1 := sameOrbit_of_orbit_eq hq.1 hqr
  exact le_antisymm (hq.2 r (mem_orbit.2 ⟨hr.1, h1⟩)) (hr.2 q (mem_orbit.2 ⟨hq.1, h1.symm⟩))

end Spqr
