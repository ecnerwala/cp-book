import Spqr.Proofs.PieceLoc
import Spqr.Proofs.PlanarAppend

namespace Spqr.Piece

def Agrees (P : Piece) (a : Array (Option Nat)) (ρ : RotationSystem) : Prop :=
  ∀ q r, P.Mem q → a[q]? = some (some r) →
    ∃ lq lr, P.loc q = some lq ∧ P.loc r = some lr ∧ ρ.get lq = some lr

theorem agrees_append {P : Piece} {E : List Nat} {a : Array (Option Nat)}
    {ρ₁ ρ₂ : RotationSystem} (hdis : List.Disjoint P.ves E)
    (hs : ρ₁.size = 4 * P.ves.length)
    (h₁ : P.Agrees a ρ₁) (h₂ : ({P with ves := E} : Piece).Agrees a ρ₂) :
    ({P with ves := P.ves ++ E} : Piece).Agrees a (ρ₁.union ρ₂) := by
  intro q r hq hr
  rcases List.mem_append.1 hq with hq | hq
  · obtain ⟨lq, lr, hq', hr', hρ⟩ := h₁ q r hq hr
    refine ⟨lq, lr, loc_append_left E hq', loc_append_left E hr', ?_⟩
    rw [RotationSystem.union_get_lt _ _ (by rw [hs]; exact loc_lt hq')]
    exact hρ
  · obtain ⟨lq, lr, hq', hr', hρ⟩ := h₂ q r hq hr
    have hnq : QE.edge q ∉ P.ves := fun h => List.disjoint_left.1 hdis h hq
    have hnr : QE.edge r ∉ P.ves := fun h => List.disjoint_left.1 hdis h (mem_of_loc hr')
    refine ⟨4 * P.ves.length + lq, 4 * P.ves.length + lr,
      loc_append_right P.ves hnq hq', loc_append_right P.ves hnr hr', ?_⟩
    rw [RotationSystem.union_get_ge _ _ (by rw [hs]; omega), hs, Nat.add_sub_cancel_left, hρ]
    rfl

end Spqr.Piece
