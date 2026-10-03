import Spqr.Proofs.PlanarOneSum

/-!
# Splitting an orbit

`swapImg f x y` with `x`, `y` on a common orbit of `f` splits that orbit in two
(`orbitCount_swapImg_same`), and `rs.conj a b` with `a ^^^ 3`, `b ^^^ 3` cofacial (and likewise
their rotation partners) gains two face orbits (`conj_numFaceOrbits_split`): inserting an edge
across a face.
-/

namespace Spqr

open Classical

theorem swapImg_swapImg (f : Nat → Nat) (x y : Nat) : swapImg (swapImg f x y) x y = f := by
  funext a
  rw [swapImg_eq, swapImg_eq, Equiv.swap_apply_self]

theorem not_sameOrbit_swapImg_of_sameOrbit {f : Nat → Nat} {S : Finset Nat} (hf : IsPermOn f S)
    {x y : Nat} (hx : x ∈ S) (hne : x ≠ y) (h : SameOrbit f x y) :
    ¬SameOrbit (swapImg f x y) x y := by
  have hy : y ∈ S := (hf.mem_iff_of_sameOrbit h).1 hx
  have hex : ∃ k, f^[k] x = y := IsPermOn.reach_iff.1 (hf.reach_of_sameOrbit hx h)
  set k := Nat.find hex with hkdef
  have hk : f^[k] x = y := Nat.find_spec hex
  have hkmin : ∀ j, j < k → f^[j] x ≠ y := fun j hj => Nat.find_min hex hj
  have hkpos : 0 < k := by
    rcases Nat.eq_zero_or_pos k with h0 | h0
    · rw [h0] at hk; exact absurd hk hne
    · exact h0
  obtain ⟨p, hp0, hp⟩ := hf.exists_iterate_eq hx
  have hpex : ∃ p, 0 < p ∧ f^[p] x = x := ⟨p, hp0, hp⟩
  set p₀ := Nat.find hpex with hp₀def
  have hp₀ : 0 < p₀ ∧ f^[p₀] x = x := Nat.find_spec hpex
  have hp₀min : ∀ j, j < p₀ → ¬(0 < j ∧ f^[j] x = x) := fun j hj => Nat.find_min hpex hj
  have hkp : k < p₀ := by
    rcases Nat.lt_or_ge k p₀ with hlt | hle
    · exact hlt
    exfalso
    apply hkmin (k - p₀) (by omega)
    rw [← hk]
    conv_rhs => rw [show k = (k - p₀) + p₀ by omega, Function.iterate_add_apply, hp₀.2]
  have hiter_ne : ∀ j, 0 < j → j ≤ k → f^[j] x ≠ x :=
    fun j hj hjk heq => hp₀min j (by omega) ⟨hj, heq⟩
  let ℓ : Nat → Prop := fun z => ∃ j, 0 < j ∧ j ≤ k ∧ f^[j] x = z
  have hℓx : ¬ℓ x := fun ⟨j, hj, hjk, heq⟩ => hiter_ne j hj hjk heq
  have hℓy : ℓ y := ⟨k, hkpos, le_refl k, hk⟩
  have hinv : ∀ z, ℓ (swapImg f x y z) ↔ ℓ z := by
    intro z
    by_cases hzx : z = x
    · subst hzx
      rw [swapImg_left]
      refine ⟨fun ⟨j, hj, hjk, heq⟩ => ?_, fun h => absurd h hℓx⟩
      exfalso
      have hsucc : f (f^[k] z) = f^[k + 1] z := (Function.iterate_succ_apply' f k z).symm
      rw [← hk, hsucc, show k + 1 = j + (k + 1 - j) by omega, Function.iterate_add_apply] at heq
      exact hiter_ne (k + 1 - j) (by omega) (by omega)
        (hf.iterate_inj j (hf.iterate_mem _ hx) hx heq.symm)
    by_cases hzy : z = y
    · subst hzy
      rw [swapImg_right f hzx]
      exact ⟨fun _ => hℓy, fun _ => ⟨1, Nat.one_pos, hkpos, rfl⟩⟩
    rw [swapImg_of_ne f hzx hzy]
    constructor
    · rintro ⟨j, hj, hjk, heq⟩
      have hzS : z ∈ S := by
        by_contra hzS
        rw [hf.fix z hzS] at heq
        exact hzS (heq ▸ hf.iterate_mem j hx)
      obtain ⟨j', rfl⟩ : ∃ j', j = j' + 1 := ⟨j - 1, by omega⟩
      have hsucc : f^[j' + 1] x = f (f^[j'] x) := Function.iterate_succ_apply' f j' x
      rw [hsucc] at heq
      have := hf.inj _ (hf.iterate_mem j' hx) z hzS heq
      rcases Nat.eq_zero_or_pos j' with h0 | h0
      · rw [h0] at this; exact absurd this.symm hzx
      · exact ⟨j', h0, by omega, this⟩
    · rintro ⟨j, hj, hjk, rfl⟩
      have hjk' : j ≠ k := fun h => hzy (h ▸ hk)
      exact ⟨j + 1, by omega, by omega, (Function.iterate_succ_apply' f j x)⟩
  intro hxy
  have := sameOrbit_invariant (fun z => ℓ z) (fun z => propext (hinv z)) hxy
  exact hℓx (this ▸ hℓy)

theorem orbitCount_swapImg_same {f : Nat → Nat} {S : Finset Nat} (hf : IsPermOn f S)
    {x y : Nat} (hx : x ∈ S) (hne : x ≠ y) (h : SameOrbit f x y) :
    orbitCount f S + 1 = orbitCount (swapImg f x y) S := by
  have hy : y ∈ S := (hf.mem_iff_of_sameOrbit h).1 hx
  have := orbitCount_swapImg (isPermOn_swapImg hf hx hy) hx hy
    (not_sameOrbit_swapImg_of_sameOrbit hf hx hne h)
  rwa [swapImg_swapImg] at this

namespace RotationSystem

variable (rs : RotationSystem) (a b : Nat)

theorem conj_orbitCount_split (ht : rs.Total) (hi : rs.Involution) (ho : rs.OppositeDir)
    (ha : a < rs.size) (hb : b < rs.size) (hab : a ≠ b) (hrab : rs.rot a ≠ b)
    {m c : Nat} (hs : rs.size = 4 * m) (hc : c < 4) (hc2 : c % 2 = 1)
    (h1 : SameOrbit (rs.stepC c) (a ^^^ c) (b ^^^ c))
    (h2 : SameOrbit (rs.stepC c) (rs.rot a ^^^ c) (rs.rot b ^^^ c)) :
    orbitCount (rs.stepC c) (Finset.range rs.size) + 2 =
      orbitCount ((rs.conj a b).stepC c) (Finset.range rs.size) := by
  rw [conj_stepC rs a b ht hi ho ha hb hab hrab hs hc]
  have hf := isPermOn_stepC ht hi hs hc
  have hpar : ∀ x y, SameOrbit (rs.stepC c) x y → x % 2 = y % 2 := fun x y h =>
    sameOrbit_invariant (· % 2) (stepC_mod2 ht ho hs hc hc2) h
  have hra := rot_lt ht hi ha
  have hrb := rot_lt ht hi hb
  have mem : ∀ z, z < rs.size → z ^^^ c ∈ Finset.range rs.size := fun z hz =>
    Finset.mem_range.2 (hs ▸ xor_lt_mul4 (hs ▸ hz) hc)
  have hne1 : a ^^^ c ≠ b ^^^ c := fun h => hab (xor_right_inj.1 h)
  have hne2 : rs.rot a ^^^ c ≠ rs.rot b ^^^ c :=
    fun h => hab (rot_inj ht hi ha hb (xor_right_inj.1 h))
  have hdis : ¬SameOrbit (rs.stepC c) (a ^^^ c) (rs.rot a ^^^ c) := by
    intro h
    have h1 := hpar _ _ h
    have h2 := xor_mod2 a hc2
    have h3 := xor_mod2 (rs.rot a) hc2
    have h4 := rot_mod2 ht ho ha
    omega
  set f := rs.stepC c with hfdef
  set f' := swapImg f (a ^^^ c) (b ^^^ c) with hf'def
  have hf' := isPermOn_swapImg hf (mem a ha) (mem b hb)
  have hf'a : f' (a ^^^ c) = f (b ^^^ c) := swapImg_left f _ _
  have hf'b : f' (b ^^^ c) = f (a ^^^ c) := swapImg_right f hne1.symm
  have hf'z : ∀ z, z ≠ a ^^^ c → z ≠ b ^^^ c → f' z = f z := fun z h1 h2 => swapImg_of_ne f h1 h2
  have h2' : SameOrbit f' (rs.rot a ^^^ c) (rs.rot b ^^^ c) := by
    refine (sameOrbit_congr (fun z => ¬SameOrbit f (a ^^^ c) z) ?_ ?_ ?_ hdis).1 h2
    · intro z
      exact not_congr ⟨fun h => h.trans (SameOrbit.step z), fun h => h.trans (SameOrbit.step z).symm⟩
    · intro z
      by_cases hz : SameOrbit f (a ^^^ c) z
      · simp only [hz, not_true_eq_false, false_iff, not_not]
        by_cases hza : z = a ^^^ c
        · rw [hza, hf'a]; exact h1.trans (SameOrbit.step _)
        by_cases hzb : z = b ^^^ c
        · rw [hzb, hf'b]; exact SameOrbit.step _
        rw [hf'z z hza hzb]; exact hz.trans (SameOrbit.step z)
      · have hza : z ≠ a ^^^ c := fun h => hz (h ▸ SameOrbit.refl f _)
        have hzb : z ≠ b ^^^ c := fun h => hz (h ▸ h1)
        rw [hf'z z hza hzb]
        exact not_congr ⟨fun h => h.trans (SameOrbit.step z), fun h => h.trans (SameOrbit.step z).symm⟩
    · intro z hz
      have hza : z ≠ a ^^^ c := fun h => hz (h ▸ SameOrbit.refl f _)
      have hzb : z ≠ b ^^^ c := fun h => hz (h ▸ h1)
      rw [hf'z z hza hzb]
  have e1 := orbitCount_swapImg_same hf (mem a ha) hne1 h1
  have e2 := orbitCount_swapImg_same hf' (mem _ hra) hne2 h2'
  simp only [hf'def] at e1 e2 ⊢
  omega

theorem conj_numFaceOrbits_split (ht : rs.Total) (hi : rs.Involution) (ho : rs.OppositeDir)
    (ha : a < rs.size) (hb : b < rs.size) (hab : a ≠ b) (hrab : rs.rot a ≠ b)
    {m : Nat} (hs : rs.size = 4 * m)
    (h1 : SameOrbit (rs.stepC 3) (a ^^^ 3) (b ^^^ 3))
    (h2 : SameOrbit (rs.stepC 3) (rs.rot a ^^^ 3) (rs.rot b ^^^ 3)) :
    rs.numFaceOrbits + 2 = (rs.conj a b).numFaceOrbits := by
  rw [numFaceOrbits_eq _ (conj_total rs a b ht ha hb) (conj_involution rs a b ht hi ha hb)
      (m := m) (by rw [conj_size, hs]), numFaceOrbits_eq _ ht hi hs, ← hs]
  exact conj_orbitCount_split rs a b ht hi ho ha hb hab hrab hs (by decide) (by decide) h1 h2

/-- After conjugating by two cofacial pairs, the two first elements lie in distinct orbits. -/
theorem conj_not_sameOrbit_split (ht : rs.Total) (hi : rs.Involution) (ho : rs.OppositeDir)
    (ha : a < rs.size) (hb : b < rs.size) (hab : a ≠ b) (hrab : rs.rot a ≠ b)
    {m c : Nat} (hs : rs.size = 4 * m) (hc : c < 4) (hc2 : c % 2 = 1)
    (h1 : SameOrbit (rs.stepC c) (a ^^^ c) (b ^^^ c))
    (h2 : SameOrbit (rs.stepC c) (rs.rot a ^^^ c) (rs.rot b ^^^ c)) :
    ¬SameOrbit ((rs.conj a b).stepC c) (a ^^^ c) (b ^^^ c) := by
  rw [conj_stepC rs a b ht hi ho ha hb hab hrab hs hc]
  have hf := isPermOn_stepC ht hi hs hc
  have hpar : ∀ x y, SameOrbit (rs.stepC c) x y → x % 2 = y % 2 := fun x y h =>
    sameOrbit_invariant (· % 2) (stepC_mod2 ht ho hs hc hc2) h
  have hra := rot_lt ht hi ha
  have hrb := rot_lt ht hi hb
  have mem : ∀ z, z < rs.size → z ^^^ c ∈ Finset.range rs.size := fun z hz =>
    Finset.mem_range.2 (hs ▸ xor_lt_mul4 (hs ▸ hz) hc)
  have hne1 : a ^^^ c ≠ b ^^^ c := fun h => hab (xor_right_inj.1 h)
  have hdis : ¬SameOrbit (rs.stepC c) (a ^^^ c) (rs.rot a ^^^ c) := by
    intro h
    have h1 := hpar _ _ h
    have h2 := xor_mod2 a hc2
    have h3 := xor_mod2 (rs.rot a) hc2
    have h4 := rot_mod2 ht ho ha
    omega
  set f := rs.stepC c with hfdef
  set f' := swapImg f (a ^^^ c) (b ^^^ c) with hf'def
  have hf' := isPermOn_swapImg hf (mem a ha) (mem b hb)
  have hf'a : f' (a ^^^ c) = f (b ^^^ c) := swapImg_left f _ _
  have hf'b : f' (b ^^^ c) = f (a ^^^ c) := swapImg_right f hne1.symm
  have hf'z : ∀ z, z ≠ a ^^^ c → z ≠ b ^^^ c → f' z = f z := fun z h1 h2 => swapImg_of_ne f h1 h2
  have h2' : SameOrbit f' (rs.rot a ^^^ c) (rs.rot b ^^^ c) := by
    refine (sameOrbit_congr (fun z => ¬SameOrbit f (a ^^^ c) z) ?_ ?_ ?_ hdis).1 h2
    · intro z
      exact not_congr ⟨fun h => h.trans (SameOrbit.step z), fun h => h.trans (SameOrbit.step z).symm⟩
    · intro z
      by_cases hz : SameOrbit f (a ^^^ c) z
      · simp only [hz, not_true_eq_false, false_iff, not_not]
        by_cases hza : z = a ^^^ c
        · rw [hza, hf'a]; exact h1.trans (SameOrbit.step _)
        by_cases hzb : z = b ^^^ c
        · rw [hzb, hf'b]; exact SameOrbit.step _
        rw [hf'z z hza hzb]; exact hz.trans (SameOrbit.step z)
      · have hza : z ≠ a ^^^ c := fun h => hz (h ▸ SameOrbit.refl f _)
        have hzb : z ≠ b ^^^ c := fun h => hz (h ▸ h1)
        rw [hf'z z hza hzb]
        exact not_congr ⟨fun h => h.trans (SameOrbit.step z), fun h => h.trans (SameOrbit.step z).symm⟩
    · intro z hz
      have hza : z ≠ a ^^^ c := fun h => hz (h ▸ SameOrbit.refl f _)
      have hzb : z ≠ b ^^^ c := fun h => hz (h ▸ h1)
      rw [hf'z z hza hzb]
  have hsplit : ¬SameOrbit f' (a ^^^ c) (b ^^^ c) :=
    not_sameOrbit_swapImg_of_sameOrbit hf (mem a ha) hne1 h1
  have hpar' : ∀ z, f' z % 2 = z % 2 := by
    intro z
    by_cases hza : z = a ^^^ c
    · rw [hza, hf'a, hfdef, stepC_mod2 ht ho hs hc hc2]; exact (hpar _ _ h1).symm
    by_cases hzb : z = b ^^^ c
    · rw [hzb, hf'b, hfdef, stepC_mod2 ht ho hs hc hc2]; exact hpar _ _ h1
    rw [hf'z z hza hzb, hfdef]; exact stepC_mod2 ht ho hs hc hc2 z
  have hdis' : ¬SameOrbit f' (rs.rot a ^^^ c) (a ^^^ c) := by
    intro h
    have h1 := sameOrbit_invariant (· % 2) hpar' h
    have h2 := xor_mod2 a hc2
    have h3 := xor_mod2 (rs.rot a) hc2
    have h4 := rot_mod2 ht ho ha
    omega
  set f'' := swapImg f' (rs.rot a ^^^ c) (rs.rot b ^^^ c) with hf''def
  have hne2 : rs.rot a ^^^ c ≠ rs.rot b ^^^ c :=
    fun h => hab (rot_inj ht hi ha hb (xor_right_inj.1 h))
  have hf''a : f'' (rs.rot a ^^^ c) = f' (rs.rot b ^^^ c) := swapImg_left f' _ _
  have hf''b : f'' (rs.rot b ^^^ c) = f' (rs.rot a ^^^ c) := swapImg_right f' hne2.symm
  have hf''z : ∀ z, z ≠ rs.rot a ^^^ c → z ≠ rs.rot b ^^^ c → f'' z = f' z :=
    fun z h1 h2 => swapImg_of_ne f' h1 h2
  intro h
  apply hsplit
  refine (sameOrbit_congr (fun z => ¬SameOrbit f' (rs.rot a ^^^ c) z) ?_ ?_ ?_ hdis').2 h
  · intro z
    exact not_congr ⟨fun h => h.trans (SameOrbit.step z), fun h => h.trans (SameOrbit.step z).symm⟩
  · intro z
    by_cases hz : SameOrbit f' (rs.rot a ^^^ c) z
    · simp only [hz, not_true_eq_false, false_iff, not_not]
      by_cases hza : z = rs.rot a ^^^ c
      · rw [hza, hf''a]; exact h2'.trans (SameOrbit.step _)
      by_cases hzb : z = rs.rot b ^^^ c
      · rw [hzb, hf''b]; exact SameOrbit.step _
      rw [hf''z z hza hzb]; exact hz.trans (SameOrbit.step z)
    · have hza : z ≠ rs.rot a ^^^ c := fun h => hz (h ▸ SameOrbit.refl f' _)
      have hzb : z ≠ rs.rot b ^^^ c := fun h => hz (h ▸ h2')
      rw [hf''z z hza hzb]
      exact not_congr ⟨fun h => h.trans (SameOrbit.step z), fun h => h.trans (SameOrbit.step z).symm⟩
  · intro z hz
    have hza : z ≠ rs.rot a ^^^ c := fun h => hz (h ▸ SameOrbit.refl f' _)
    have hzb : z ≠ rs.rot b ^^^ c := fun h => hz (h ▸ h2')
    rw [hf''z z hza hzb]

end RotationSystem

end Spqr
