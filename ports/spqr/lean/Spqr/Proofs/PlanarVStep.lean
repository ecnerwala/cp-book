import Spqr.Proofs.PlanarSplice
import Spqr.Proofs.PlanarInsert

/-!
# One-sum with an open piece, pointwise

`vstep`: `IsPlanarEmbedding.splice` (the 1-sum `(ρ₁.union ρ₂).conj a (ρ₁.size + w0)` of two
embeddings sharing exactly the vertex `x`) together with its rotation given pointwise: the
corner `a ↔ ρ₁.rot a` of `ρ₁` is opened and rewired to the boundary pair `w0 ↔ w1` of `ρ₂`
(`a ↔ w1`, `ρ₁.rot a ↔ w0`), everything else is unchanged (the `ρ₂` part shifted by `ρ₁.size`).
-/

namespace Spqr

open RotationSystem

theorem vstep {es₁ es₂ : List (Nat × Nat)} {n : Nat} {ρ₁ ρ₂ : RotationSystem}
    (h₁ : IsPlanarEmbedding es₁ n ρ₁) (h₂ : IsPlanarEmbedding es₂ n ρ₂) {x a w0 w1 : Nat}
    (ha : a < ρ₁.size) (hw0 : ρ₂.get w0 = some w1) (hw0lt : w0 < ρ₂.size) (hab : a % 2 = w0 % 2)
    (hva : QE.vert es₁ a = some x) (hvw : QE.vert es₂ w0 = some x)
    (hsep : ∀ w, HasEdge es₁ w → HasEdge es₂ w → w = x) :
    IsPlanarEmbedding (es₁ ++ es₂) n ((ρ₁.union ρ₂).conj a (ρ₁.size + w0)) ∧
    ((ρ₁.union ρ₂).conj a (ρ₁.size + w0)).size = ρ₁.size + ρ₂.size ∧
    (∀ q, q < ρ₁.size → ((ρ₁.union ρ₂).conj a (ρ₁.size + w0)).rot q =
      if q = a then ρ₁.size + w1 else if q = ρ₁.rot a then ρ₁.size + w0 else ρ₁.rot q) ∧
    (∀ r, r < ρ₂.size → ((ρ₁.union ρ₂).conj a (ρ₁.size + w0)).rot (ρ₁.size + r) =
      if r = w0 then ρ₁.rot a else if r = w1 then a else ρ₁.size + ρ₂.rot r) := by
  have hpl : IsPlanarEmbedding (es₁ ++ es₂) n ((ρ₁.union ρ₂).conj a (ρ₁.size + w0)) :=
    h₁.splice h₂ (by rw [← h₁.size]; exact ha) (by rw [← h₂.size]; exact hw0lt) hab hva hvw hsep
  have hsz : ((ρ₁.union ρ₂).conj a (ρ₁.size + w0)).size = ρ₁.size + ρ₂.size := by
    rw [conj_size, union_size]
  have hw1 : ρ₂.rot w0 = w1 := rot_eq_of_get hw0
  have hw1lt : w1 < ρ₂.size := hw1 ▸ rot_lt h₂.total h₂.involution hw0lt
  have hw10 : w1 ≠ w0 := hw1 ▸ rot_ne h₂.total h₂.opposite_dir hw0lt
  have hw1w : ρ₂.rot w1 = w0 := by rw [← hw1]; exact rot_rot h₂.total h₂.involution hw0lt
  have hra : ρ₁.rot a < ρ₁.size := rot_lt h₁.total h₁.involution ha
  have hraa : ρ₁.rot a ≠ a := rot_ne h₁.total h₁.opposite_dir ha
  have hrra : ρ₁.rot (ρ₁.rot a) = a := rot_rot h₁.total h₁.involution ha
  set b := ρ₁.size + w0 with hb
  have hsw : ∀ y, y ≠ a → y ≠ b → Equiv.swap a b y = y := fun y h1 h2 =>
    Equiv.swap_apply_of_ne_of_ne h1 h2
  refine ⟨hpl, hsz, ?_, ?_⟩
  · intro q hq
    have hq' : q < ((ρ₁.union ρ₂)).size := by rw [union_size]; omega
    apply rot_eq_of_get
    rw [conj_get_lt (hq := hq')]
    by_cases hqa : q = a
    · subst hqa
      rw [Equiv.swap_apply_left, union_get_ge (hq := by omega), hb, Nat.add_sub_cancel_left, hw0]
      simp only [Option.map_some, ite_true]
      rw [hsw _ (by omega) (by omega)]
    rw [ite_eq_right hqa, hsw q hqa (by omega), union_get_lt (hq := hq), get_eq_rot h₁.total hq]
    by_cases hqr : q = ρ₁.rot a
    · subst hqr
      rw [hrra, ite_eq_left rfl]
      simp only [Option.map_some, Equiv.swap_apply_left]
    rw [ite_eq_right hqr]
    simp only [Option.map_some]
    rw [hsw _ (fun h => hqr (by rw [← h, rot_rot h₁.total h₁.involution hq])) (by
      have := rot_lt h₁.total h₁.involution hq; omega)]
  · intro r hr
    have hq' : ρ₁.size + r < ((ρ₁.union ρ₂)).size := by rw [union_size]; omega
    apply rot_eq_of_get
    rw [conj_get_lt (hq := hq')]
    by_cases hr0 : r = w0
    · subst hr0
      rw [← hb, Equiv.swap_apply_right, union_get_lt (hq := ha), get_eq_rot h₁.total ha, ite_eq_left rfl]
      simp only [Option.map_some]
      rw [hsw _ hraa (by omega)]
    rw [ite_eq_right hr0, hsw _ (by omega) (by omega), union_get_ge (hq := by omega),
      Nat.add_sub_cancel_left, get_eq_rot h₂.total hr]
    simp only [Option.map_some]
    by_cases hr1 : r = w1
    · subst hr1
      rw [hw1w, ite_eq_left rfl, ← hb, Equiv.swap_apply_right]
    rw [ite_eq_right hr1, hsw _ (by omega) (fun h => hr1 (by
      have : ρ₂.rot r = w0 := by omega
      rw [← rot_rot h₂.total h₂.involution hr, this, hw1]))]

end Spqr
