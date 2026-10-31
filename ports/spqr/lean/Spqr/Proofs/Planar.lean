import Mathlib.Data.Finset.Card
import Mathlib.Data.Nat.Bitwise
import Spqr.Planar

/-!
# Lemmas about the rotation-system spec

Orbit counting (`numOrbits`) through orbit minima: a point is counted iff it is the least of
its `n`-step forward orbit, so a label that `step` preserves gives the minima of a class, an
explicit iterate below `q` shows `q` is not a minimum, and a fixed-point-free involution on
`[0, n)` has exactly `n / 2` orbits. Also the arithmetic forms of `QE.flipDir` / `QE.across`.
-/

namespace Spqr

theorem isOrbitMin_foldl (step : Nat → Option Nat) (q : Nat) : ∀ (m i : Nat) (ok : Bool),
    ((List.range' i m).foldl (fun (p : Nat × Bool) _ => (stepFn step p.1, p.2 && decide (q ≤ p.1)))
      ((stepFn step)^[i] q, ok)).2 = true ↔
      ok = true ∧ ∀ j, i ≤ j → j < i + m → q ≤ (stepFn step)^[j] q := by
  intro m
  induction m with
  | zero =>
    intro i ok
    simp only [List.range'_zero, List.foldl_nil, Nat.add_zero]
    constructor
    · intro h; exact ⟨h, fun j h1 h2 => absurd h1 (by omega)⟩
    · intro h; exact h.1
  | succ m ih =>
    intro i ok
    rw [List.range'_succ, List.foldl_cons]
    simp only
    rw [← Function.iterate_succ_apply' (stepFn step) i q, ih (i + 1)]
    constructor
    · rintro ⟨h1, h2⟩
      simp only [Bool.and_eq_true, decide_eq_true_eq] at h1
      refine ⟨h1.1, fun j hj1 hj2 => ?_⟩
      rcases Nat.eq_or_lt_of_le hj1 with rfl | h
      · exact h1.2
      · exact h2 j h (by omega)
    · rintro ⟨h1, h2⟩
      refine ⟨?_, fun j hj1 hj2 => h2 j (by omega) (by omega)⟩
      simp only [Bool.and_eq_true, decide_eq_true_eq]
      exact ⟨h1, h2 i le_rfl (by omega)⟩

theorem isOrbitMin_iff (step : Nat → Option Nat) (n q : Nat) :
    isOrbitMin step n q = true ↔ ∀ i, i < n → q ≤ (stepFn step)^[i] q := by
  unfold isOrbitMin
  rw [List.range_eq_range']
  have := isOrbitMin_foldl step q n 0 true
  simp only [Function.iterate_zero, id_eq, Nat.zero_add] at this
  rw [this]
  simp

theorem isOrbitMin_of_label (step : Nat → Option Nat) (n q : Nat) (c : Nat → Nat)
    (hc : ∀ x, c (stepFn step x) = c x) (hmin : ∀ x, c x = c q → q ≤ x) :
    isOrbitMin step n q = true := by
  rw [isOrbitMin_iff]
  intro i hi
  clear hi
  apply hmin
  induction i with
  | zero => simp
  | succ i ih => rw [Function.iterate_succ_apply', hc]; exact ih

theorem not_isOrbitMin_of_iter (step : Nat → Option Nat) (n q i : Nat) (hi : i < n)
    (h : (stepFn step)^[i] q < q) : isOrbitMin step n q = false := by
  rw [Bool.eq_false_iff, Ne, isOrbitMin_iff]
  intro H
  exact absurd (H i hi) (not_le.mpr h)

theorem length_filter_range_eq_card (n : Nat) (b : Nat → Bool) :
    ((List.range n).filter b).length = ((Finset.range n).filter (fun q => b q = true)).card := by
  rw [← List.toFinset_card_of_nodup ((List.nodup_range).filter b)]
  congr 1
  ext q
  simp

theorem numOrbits_eq_card (step : Nat → Option Nat) (n : Nat) (p : Nat → Prop) [DecidablePred p]
    (h : ∀ q, q < n → (isOrbitMin step n q = true ↔ p q)) :
    numOrbits step n = ((Finset.range n).filter p).card := by
  unfold numOrbits
  rw [length_filter_range_eq_card]
  congr 1
  ext q
  simp only [Finset.mem_filter, Finset.mem_range]
  exact and_congr_right (h q)

/-- A fixed-point-free involution of `[0, n)` has `n / 2` orbits. -/
theorem numOrbits_involution (step : Nat → Option Nat) (n : Nat)
    (hlt : ∀ q, q < n → stepFn step q < n)
    (hinv : ∀ q, q < n → stepFn step (stepFn step q) = q)
    (hne : ∀ q, q < n → stepFn step q ≠ q) :
    numOrbits step n = n / 2 := by
  have hiter : ∀ q, q < n → ∀ i, (stepFn step)^[i] q = if i % 2 = 0 then q else stepFn step q := by
    intro q hq i
    induction i with
    | zero => simp
    | succ i ih =>
      rw [Function.iterate_succ_apply', ih]
      rcases Nat.mod_two_eq_zero_or_one i with h | h
      · have h' : (i + 1) % 2 = 1 := by omega
        simp [h, h']
      · have h' : (i + 1) % 2 = 0 := by omega
        simp [h, h', hinv q hq]
  have hmin : ∀ q, q < n → (isOrbitMin step n q = true ↔ q < stepFn step q) := by
    intro q hq
    have hn : 2 ≤ n := by
      by_contra hn
      have h0 : q = 0 := by omega
      have := hlt q hq
      have := hne q hq
      omega
    rw [isOrbitMin_iff]
    constructor
    · intro H
      have := H 1 (by omega)
      simp only [Function.iterate_one] at this
      exact lt_of_le_of_ne this (hne q hq).symm
    · intro H i _
      rw [hiter q hq i]
      split <;> omega
  rw [numOrbits_eq_card step n (fun q => q < stepFn step q) hmin]
  have hsplit := Finset.card_filter_add_card_filter_not (s := Finset.range n) (fun q => q < stepFn step q)
  have hcomp : (Finset.range n).filter (fun q => ¬ q < stepFn step q) =
      (Finset.range n).filter (fun q => stepFn step q < q) := by
    ext q
    simp only [Finset.mem_filter, Finset.mem_range]
    constructor
    · rintro ⟨hq, h⟩; exact ⟨hq, lt_of_le_of_ne (not_lt.mp h) (hne q hq)⟩
    · rintro ⟨hq, h⟩; exact ⟨hq, not_lt.mpr h.le⟩
  have hbij : ((Finset.range n).filter (fun q => q < stepFn step q)).card =
      ((Finset.range n).filter (fun q => stepFn step q < q)).card := by
    apply Finset.card_bij' (fun q _ => stepFn step q) (fun q _ => stepFn step q)
    · intro q hq
      simp only [Finset.mem_filter, Finset.mem_range] at hq ⊢
      exact ⟨hlt q hq.1, by rw [hinv q hq.1]; exact hq.2⟩
    · intro q hq
      simp only [Finset.mem_filter, Finset.mem_range] at hq ⊢
      exact ⟨hlt q hq.1, by rw [hinv q hq.1]; exact hq.2⟩
    · intro q hq
      simp only [Finset.mem_filter, Finset.mem_range] at hq
      exact hinv q hq.1
    · intro q hq
      simp only [Finset.mem_filter, Finset.mem_range] at hq
      exact hinv q hq.1
  rw [hcomp, ← hbij, Finset.card_range] at hsplit
  omega

/-! ### Quarter-edge arithmetic -/

namespace QE

theorem flipDir_eq (q : Nat) : flipDir q = 2 * (q / 2) + (1 - q % 2) := by
  unfold flipDir
  rcases Nat.even_or_odd q with h | h
  · rw [Nat.xor_one_of_even h]; obtain ⟨k, hk⟩ := h; omega
  · rw [Nat.xor_one_of_odd h]; obtain ⟨k, hk⟩ := h; omega

theorem across_eq (q : Nat) : across q = 4 * (q / 4) + (3 - q % 4) := by
  unfold across
  apply Nat.eq_of_testBit_eq
  intro i
  rw [Nat.testBit_xor]
  match i with
  | 0 =>
    simp only [Nat.testBit_zero]
    by_cases h : q % 2 = 1 <;> simp [h] <;> omega
  | 1 =>
    simp only [Nat.testBit_succ, Nat.testBit_zero]
    by_cases h : q / 2 % 2 = 1 <;> simp [h] <;> omega
  | i + 2 =>
    rw [Nat.testBit_lt_two_pow (x := 3) (by
      have := Nat.one_le_two_pow (n := i)
      rw [Nat.pow_succ, Nat.pow_succ]; omega), Bool.xor_false]
    simp only [Nat.testBit_succ]
    congr 1
    omega

theorem edge_mk (ne r : Nat) (_hr : r < 4) : edge (4 * ne + r) = ne := by unfold edge; omega
theorem side_mk (ne r : Nat) (_hr : r < 4) : side (4 * ne + r) = r / 2 := by unfold side; omega
theorem dir_mk (ne r : Nat) (hr : r < 4) : dir (4 * ne + r) = r % 2 := by unfold dir; omega

end QE

end Spqr
