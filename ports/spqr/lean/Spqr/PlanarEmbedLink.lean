import Spqr.PlanarEmbedSteps

/-!
# Frame lemmas for the primitives of `embedItem`

`link`, `outer` and `setOuter` are the only state operations of `embedItem`; the per-item steps
reason about them through these lemmas on an abstract `EmbedState` instead of unfolding the
`do`-blocks.
-/

namespace Spqr

namespace PlanarSpqrTree

open EmbedM

theorem outer_run (i k : Nat) (s : EmbedState) : (outer i k).run s = (s.outerE[i]![k]!, s) := rfl

theorem link_none_left (b : Option Nat) (s : EmbedState) : (link none b).run s = ((), s) := by
  cases b <;> rfl

theorem link_none_right (a : Option Nat) (s : EmbedState) : (link a none).run s = ((), s) := by
  cases a <;> rfl

theorem link_some_run (a b : Nat) (s : EmbedState) :
    (link (some a) (some b)).run s =
      ((), { s with rotAdj := (s.rotAdj.set! a (some b)).set! b (some a) }) := rfl

theorem link_outerE (a b : Option Nat) (s : EmbedState) : ((link a b).run s).2.outerE = s.outerE := by
  cases a <;> cases b <;> rfl

theorem link_rotAdj_size (a b : Option Nat) (s : EmbedState) :
    ((link a b).run s).2.rotAdj.size = s.rotAdj.size := by
  cases a <;> cases b <;> first | rfl | (rw [link_some_run]; simp)

theorem link_exposedAt (a b : Option Nat) (s : EmbedState) (j q : Nat) :
    ((link a b).run s).2.exposedAt j q ↔ s.exposedAt j q := by
  simp [EmbedState.exposedAt, link_outerE]

/-- The rotation after `link (some a) (some b)` with `a ≠ b` in range. -/
theorem link_rotAdj_get (a b : Nat) (s : EmbedState) (ha : a < s.rotAdj.size) (hb : b < s.rotAdj.size)
    (hab : a ≠ b) (q : Nat) :
    ((link (some a) (some b)).run s).2.rotAdj[q]? =
      if q = a then some (some b) else if q = b then some (some a) else s.rotAdj[q]? := by
  rw [link_some_run]
  simp only [Array.set!, Array.getElem?_setIfInBounds, Array.size_setIfInBounds]
  by_cases hq : q = a
  · subst hq; simp [hab.symm, ha]
  by_cases hq' : q = b
  · subst hq'; simp [hb, Ne.symm hab]
  · simp [hq, hq', Ne.symm hq, Ne.symm hq', ]

theorem setOuter_rotAdj (i k : Nat) (q : Option Nat) (s : EmbedState) :
    ((setOuter i k q).run s).2.rotAdj = s.rotAdj := rfl

theorem setOuter_run (i k : Nat) (q : Option Nat) (s : EmbedState) :
    (setOuter i k q).run s = ((), { s with outerE := s.outerE.modify i (·.set! k q) }) := rfl

theorem setOuter_outerE_size (i k : Nat) (q : Option Nat) (s : EmbedState) :
    ((setOuter i k q).run s).2.outerE.size = s.outerE.size := by
  rw [setOuter_run]; simp

theorem setOuter_outerE_ne (i k : Nat) (q : Option Nat) (s : EmbedState) (j : Nat) (hj : j ≠ i) :
    ((setOuter i k q).run s).2.outerE[j]? = s.outerE[j]? := by
  rw [setOuter_run]; simp [Array.getElem?_modify, Ne.symm hj]

theorem setOuter_outerE_eq (i k : Nat) (q : Option Nat) (s : EmbedState) (hi : i < s.outerE.size) :
    ((setOuter i k q).run s).2.outerE[i]? = some (s.outerE[i].set! k q) := by
  rw [setOuter_run]; simp [Array.getElem_modify_self, hi]

theorem setOuter_exposedAt_ne (i k : Nat) (q : Option Nat) (s : EmbedState) (j r : Nat) (hj : j ≠ i) :
    ((setOuter i k q).run s).2.exposedAt j r ↔ s.exposedAt j r := by
  simp [EmbedState.exposedAt, setOuter_outerE_ne i k q s j hj]

end PlanarSpqrTree

end Spqr
