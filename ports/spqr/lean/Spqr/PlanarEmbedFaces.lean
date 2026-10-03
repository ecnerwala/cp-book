import Spqr.PlanarEmbedTree

/-!
# Cofacial exposed cap ends

`PlanarEmbedFaceCounterexample.lean` shows `GluedUpTo` is too weak for the `Q` step: the two
exposed pairs of a capped piece must lie on a common face of the planar certificate that agrees
with `rotAdj`. `GluedFaces` adds this, for the same witness, on every maximal capped piece.
-/

namespace Spqr

namespace PlanarSpqrTree

variable (t : PlanarSpqrTree)

/-- `ρ` certifies the piece below `j` in state `s`: planar, agreeing with `rotAdj`, with exactly
the exposed ends unset and each exposed pair facing (the body of `GluedPieces.piece`). -/
def PieceWitness (g : Graph) (s : EmbedState) (j : Nat) (ρ : RotationSystem) : Prop :=
  IsPlanarEmbedding (t.pieceBelow g j).es g.nv ρ ∧
  (∀ q r, (t.pieceBelow g j).Mem q → s.rotAdj[q]? = some (some r) →
    ∃ lq lr, (t.pieceBelow g j).loc q = some lq ∧ (t.pieceBelow g j).loc r = some lr ∧
      ρ.get lq = some lr) ∧
  (∀ q, (t.pieceBelow g j).Mem q → (s.rotAdj[q]? = some none ↔ s.exposedAt j q)) ∧
  (∀ k, k < 2 →
    (∀ a, s.outerE[j]?.bind (fun o => o[2 * k]?) = some (some a) →
      ∃ b, s.outerE[j]?.bind (fun o => o[2 * k + 1]?) = some (some b) ∧
        ∃ la lb, (t.pieceBelow g j).loc a = some la ∧ (t.pieceBelow g j).loc b = some lb ∧
          ρ.get la = some lb) ∧
    (∀ b, s.outerE[j]?.bind (fun o => o[2 * k + 1]?) = some (some b) →
      ∃ a, s.outerE[j]?.bind (fun o => o[2 * k]?) = some (some a)))

/-- The exposed ends in slots 0 and 2 of row `j` are on a common face of `ρ`. -/
def CapFace (g : Graph) (s : EmbedState) (j : Nat) (ρ : RotationSystem) : Prop :=
  ∀ a c, s.outerE[j]?.bind (fun o => o[0]?) = some (some a) →
    s.outerE[j]?.bind (fun o => o[2]?) = some (some c) →
    ∃ la lc, (t.pieceBelow g j).loc a = some la ∧ (t.pieceBelow g j).loc c = some lc ∧
      ρ.SameFaceOrbit la lc

/-- `GluedUpTo` plus cofaciality of the two exposed pairs of every maximal capped piece, for the
certificate that agrees with `rotAdj`. -/
structure GluedFaces (g : Graph) (i : Nat) (s : EmbedState) : Prop
    extends t.GluedUpTo g i s where
  cap_face : ∀ j, t.Maximal i j → ∀ ne, t.toSpqrTree.capNe j = some ne →
    ∃ ρ, t.PieceWitness g s j ρ ∧ t.CapFace g s j ρ

theorem pieceWitness_frame {g : Graph} {s s' : EmbedState} {j : Nat} {ρ : RotationSystem}
    (h : t.PieceWitness g s j ρ)
    (hrot : ∀ q, (t.pieceBelow g j).Mem q → s'.rotAdj[q]? = s.rotAdj[q]?)
    (hout : s'.outerE[j]? = s.outerE[j]?) : t.PieceWitness g s' j ρ := by
  obtain ⟨hρ, ha, ho, hp⟩ := h
  refine ⟨hρ, ?_, ?_, ?_⟩
  · intro q r hq hr; rw [hrot q hq] at hr; exact ha q r hq hr
  · intro q hq; rw [hrot q hq, ho q hq]; simp only [EmbedState.exposedAt, hout]
  · simpa only [hout] using hp

theorem capFace_frame {g : Graph} {s s' : EmbedState} {j : Nat} {ρ : RotationSystem}
    (h : t.CapFace g s j ρ) (hout : s'.outerE[j]? = s.outerE[j]?) : t.CapFace g s' j ρ := by
  simpa only [CapFace, hout] using h

theorem maximal_succ_of_ne {i j : Nat} (h : t.Maximal i j) (hij : j ≠ i) :
    t.Maximal (i + 1) j :=
  ⟨by have := h.1; omega, h.2.1, fun p hp => Nat.lt_succ_of_lt (h.2.2 p hp)⟩

/-- A step that frames every other maximal piece preserves `GluedFaces`, given the certificate
for the processed item itself. -/
theorem gluedFaces_of_frame {g : Graph} {i : Nat} {s s' : EmbedState}
    (h : t.GluedFaces g (i + 1) s) (h' : t.GluedUpTo g i s')
    (hself : ∀ ne, t.toSpqrTree.capNe i = some ne →
      ∃ ρ, t.PieceWitness g s' i ρ ∧ t.CapFace g s' i ρ)
    (hrot : ∀ j, t.Maximal i j → j ≠ i → ∀ q, (t.pieceBelow g j).Mem q →
      s'.rotAdj[q]? = s.rotAdj[q]?)
    (hout : ∀ j, j ≠ i → s'.outerE[j]? = s.outerE[j]?) : t.GluedFaces g i s' := by
  refine ⟨h', ?_⟩
  intro j hm ne hne
  by_cases hji : j = i
  · subst hji; exact hself ne hne
  · obtain ⟨ρ, hw, hf⟩ := h.cap_face j (t.maximal_succ_of_ne hm hji) ne hne
    exact ⟨ρ, t.pieceWitness_frame hw (hrot j hm hji) (hout j hji),
      t.capFace_frame hf (hout j hji)⟩

theorem gluedFaces_init (g : Graph) : t.GluedFaces g t.size t.initState := by
  refine ⟨t.gluedUpTo_init g, ?_⟩
  intro j hm
  have := hm.2.1; have := hm.1
  omega

/-- A capped item whose row is unexposed is certified by any piece witness. -/
theorem capFace_of_unexposed {g : Graph} {s : EmbedState} {i : Nat} (ρ : RotationSystem)
    (hun : ∀ q, ¬ s.exposedAt i q) : t.CapFace g s i ρ := by
  intro a c ha _
  exact (hun a ⟨0, ha⟩).elim

theorem maximal_self (hwf : t.toSpqrTree.WF) {i : Nat} (hi : i < t.size) : t.Maximal i i :=
  ⟨le_rfl, hi, fun p hp => t.parent_lt hwf hi ((t.parent_some_iff i p).2 hp)⟩

end PlanarSpqrTree

end Spqr
