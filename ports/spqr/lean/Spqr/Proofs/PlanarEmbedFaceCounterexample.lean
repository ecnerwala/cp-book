import Spqr.Proofs.PlanarEmbedCapCounterexample

open Spqr Spqr.PlanarSpqrTree

namespace Spqr.PlanarEmbedCounterexample

def tripleGraph : Graph := ⟨2, #[(0, 1), (0, 1), (0, 1)]⟩
def tripleTree : PlanarSpqrTree := {
  nv := 2, ne := 3, vertIndex := #[some 1, some 6], edgeIndex := #[some 2, some 5, some 4],
  edgeFlipped := #[false, false, false], par := #[none, some 0, some 1, some 2, some 3, some 3, some 2],
  subtreeEnd := #[7, 7, 7, 6, 5, 6, 7], types := #[.F, .V, .Q, .P, .Q, .Q, .V],
  origId := #[none, some 0, some 0, none, some 2, some 1, some 1],
  chBounds := #[0, 1, 2, 4, 6, 6, 6, 6], chDat := #[1, 2, 3, 6, 4, 5],
  nodeVerts := #[⟨0, 1⟩, ⟨2, 1⟩, ⟨2, 6⟩, ⟨3, 1⟩, ⟨3, 6⟩, ⟨4, 1⟩, ⟨4, 6⟩, ⟨5, 1⟩, ⟨5, 6⟩],
  nvBounds := #[0, 1, 1, 3, 5, 7, 9, 9], vertParNv := #[none, some 0, none, none, none, none, some 2],
  nodeEdges := #[⟨2, some 1, (1, 2)⟩, ⟨3, some 0, (3, 4)⟩, ⟨3, some 4, (3, 4)⟩,
    ⟨3, some 5, (3, 4)⟩, ⟨4, some 2, (5, 6)⟩, ⟨5, some 3, (7, 8)⟩],
  neBounds := #[0, 0, 0, 1, 4, 5, 6, 6],
  adjBounds := #[0, 0, 0, 0, 1, 2, 2, 2, 5, 8, 8, 8, 9, 10, 10, 10, 11, 12, 12],
  adjDat := #[⟨0, 2⟩, ⟨0, 1⟩, ⟨1, 4⟩, ⟨2, 4⟩, ⟨3, 4⟩, ⟨3, 3⟩, ⟨2, 3⟩,
    ⟨1, 3⟩, ⟨4, 6⟩, ⟨4, 5⟩, ⟨5, 8⟩, ⟨5, 7⟩],
  nodePlanar := #[true, true, true, true, true, true, true],
  neRotAdj := #[some 1, some 0, some 3, some 2, some 13, some 8, some 11, some 14,
    some 5, some 12, some 15, some 6, some 9, some 4, some 7, some 10,
    some 17, some 16, some 19, some 18, some 21, some 20, some 23, some 22] }

def doubleRot : RotationSystem := ⟨#[some 5, some 4, some 7, some 6, some 1, some 0, some 3, some 2]⟩

def badFaceState : EmbedState :=
  ⟨#[none, none, none, none, some 9, none, some 11, none, none, some 4, none, some 6],
    #[#[none, none, none, none], #[none, none, none, none], #[none, none, none, none],
      #[some 8, some 5, some 10, some 7], #[some 8, some 9, some 10, some 11],
      #[some 4, some 5, some 6, some 7], #[none, none, none, none]]⟩

theorem doubleRot_planar : IsPlanarEmbedding [(0, 1), (0, 1)] 2 doubleRot := by
  refine ⟨⟨rfl, by simp, ?_, ?_, ?_, ?_, by decide⟩, by decide⟩
  · intro q hq
    change q < 8 at hq
    interval_cases q <;> decide
  all_goals
    intro q hq r hr
    change q < 8 at hq
    interval_cases q <;> simp [RotationSystem.get, doubleRot] at hr <;> subst r <;> decide

set_option maxHeartbeats 1600000 in
theorem badFaceState_before : tripleTree.GluedUpTo tripleGraph 3 badFaceState := by
  have hslots : ∀ j k q, badFaceState.outerE[j]?.bind (fun o => o[k]?) = some (some q) →
      3 ≤ j ∧ j ≤ 5 ∧ k < 4 ∧ q < 12 := by
    intro j k q hk
    have hj : j < 7 := by
      by_contra hn
      rw [Array.getElem?_eq_none (by exact Nat.le_of_not_gt hn)] at hk
      cases hk
    interval_cases j <;> simp [badFaceState] at hk
    all_goals
      have hkl : k < 4 := by
        by_contra hn
        rw [List.getElem?_eq_none (by simpa using Nat.le_of_not_gt hn)] at hk
        cases hk
      interval_cases k <;> simp_all <;> omega
  refine ⟨⟨⟨⟨⟨⟨by decide, by decide, ?_, ?_, ?_⟩, ?_, ?_⟩, ?_⟩, ?_, ?_⟩, ?_⟩, ?_⟩
  · intro j hj q ⟨k, hk⟩
    have := (hslots j k q hk).1
    omega
  · intro q hq
    have h0 := hq 3 (by omega) (by decide)
    change q / 4 ∉ [2, 1] at h0
    by_cases hb : q < 12
    · have hh : q < 4 := by simp at h0; omega
      interval_cases q <;> simp [badFaceState]
    · exact Or.inl (Array.getElem?_eq_none (by exact Nat.le_of_not_gt hb))
  · intro j hj
    have hjs : j < 7 := hj.2.1
    have hj' : j = 3 ∨ j = 6 := by
      interval_cases j <;> simp_all [PlanarSpqrTree.Maximal, tripleTree]
    rcases hj' with rfl | rfl
    · refine ⟨doubleRot, doubleRot_planar, ?_, ?_, ?_⟩
      · intro q r hq hr
        change q / 4 ∈ [2, 1] at hq
        have hqb : q < 12 := by simp at hq; omega
        interval_cases q <;> norm_num at hq <;> simp [badFaceState] at hr <;> subst r
        · exact ⟨4, 1, by decide, by decide, by decide⟩
        · exact ⟨6, 3, by decide, by decide, by decide⟩
        · exact ⟨1, 4, by decide, by decide, by decide⟩
        · exact ⟨3, 6, by decide, by decide, by decide⟩
      · intro q hq
        change q / 4 ∈ [2, 1] at hq
        have hqb : q < 12 := by simp at hq; omega
        interval_cases q <;> norm_num at hq
        all_goals
          simp only [badFaceState, EmbedState.exposedAt]
          change _ ↔ ∃ k, _[k]? = some (some _)
          rw [← Array.mem_iff_getElem?]
          simp
      · intro k hk
        interval_cases k <;> simp [badFaceState]
        · exact ⟨0, by decide, 5, by decide, by decide⟩
        · exact ⟨2, by decide, 7, by decide, by decide⟩
    · refine ⟨⟨#[]⟩, isPlanarEmbedding_nil 2, ?_, ?_, ?_⟩
      · intro q r hq; change QE.edge q ∈ [] at hq; cases hq
      · intro q hq; change QE.edge q ∈ [] at hq; cases hq
      · intro k hk; interval_cases k <;> simp [badFaceState]
  · intro j hj
    change j < 7 at hj
    interval_cases j <;> decide
  · intro j k q hk
    obtain ⟨hjl, hjh, hkb, _⟩ := hslots j k q hk
    refine ⟨hkb, ?_, ?_⟩
    · interval_cases j <;> decide
    · intro p hp hpt
      interval_cases j <;> simp [SpqrTree.parent, tripleTree] at hp <;> subst p <;>
        simp [SpqrTree.type, tripleTree] at hpt
  · intro j p v q hp hpt hv ⟨k, hk⟩
    obtain ⟨hjl, hjh, _, _⟩ := hslots j k q hk
    interval_cases j <;> simp [SpqrTree.parent, tripleTree] at hp <;> subst p <;> cases hpt
  · intro j k q hk
    obtain ⟨hjl, hjh, hkb, _⟩ := hslots j k q hk
    interval_cases j <;> interval_cases k <;> simp [badFaceState] at hk <;> subst q <;> decide
  · intro j hj hjs p hp hpt _
    change j < 7 at hjs
    interval_cases j <;> simp [SpqrTree.parent, tripleTree] at hp <;> subst p <;>
      simp_all [SpqrTree.type, tripleTree]
  · intro j hj hjt v hv
    refine ⟨?_, ?_, ?_⟩
    · intro q ⟨k, hk⟩
      obtain ⟨hjl, hjh, _, _⟩ := hslots j k q hk
      interval_cases j <;> cases hjt
    · intro k q hk
      obtain ⟨hjl, hjh, _, _⟩ := hslots j k q hk
      interval_cases j <;> cases hjt
    · intro hge hne
      change j < 7 at hj
      interval_cases j <;> simp [SpqrTree.type, tripleTree] at hjt
      exact (hne (by rfl)).elim
  · intro j hj ne hn p hp
    change j < 7 at hj
    have hcaps : ∀ j, j < 7 → tripleTree.toSpqrTree.capNe j =
        [none, none, none, some 1, some 4, some 5, none][j]! := by
      intro j hj; interval_cases j <;> decide
    rw [hcaps j hj] at hn
    interval_cases j <;> simp at hn <;> subst ne
    all_goals
      have hp' : p = (0, 1) := by exact Option.some.inj hp.symm
      subst p
    all_goals
      refine ⟨?_, ?_⟩
      · intro k q hk
        have hkb := (hslots _ k q hk).2.2.1
        interval_cases k <;> simp [badFaceState] at hk <;> subst q <;> decide
      · intro _ _ k hk
        interval_cases k <;> exact ⟨_, rfl⟩

def badFaceClosed : RotationSystem :=
  ⟨#[some 9, some 4, some 11, some 6, some 1, some 8, some 3, some 10,
    some 5, some 0, some 7, some 2]⟩

theorem badFaceState_after :
    ¬tripleTree.GluedPieces tripleGraph 2 ((tripleTree.embedItem 2).run badFaceState).2 := by
  intro h
  obtain ⟨ρ, hρ, ha, _, hb⟩ := h.piece 2 ⟨by decide, by decide, by
    intro p hp
    have : p = 1 := by simpa [tripleTree] using hp.symm
    omega⟩
  obtain ⟨b, hb, a', b', hla, hlb, hab⟩ := (hb 0 (by decide)).1 0 (by decide)
  have hb' : b = 5 := Option.some.inj (Option.some.inj hb.symm)
  subst b
  have haa : a' = 0 := Option.some.inj (hla.symm.trans (by decide))
  have hbb : b' = 9 := Option.some.inj (hlb.symm.trans (by decide))
  subst a'; subst b'
  have h9 := (hρ.involution 0 (by rw [hρ.size]; decide) 9 hab).2
  let glob := fun q => if q < 4 then q else if q < 8 then q + 4 else q - 4
  have hgets : ∀ q, q < 12 → ρ.get q = badFaceClosed.get q := by
    intro q hq
    by_cases h0 : q = 0
    · subst q; exact hab
    by_cases hn : q = 9
    · subst q; exact h9
    have hmem : (tripleTree.pieceBelow tripleGraph 2).Mem (glob q) := by
      change glob q / 4 ∈ [0, 2, 1]
      interval_cases q <;> decide
    have hlookup : ((tripleTree.embedItem 2).run badFaceState).2.rotAdj[glob q]? =
        some (some (glob ((badFaceClosed.get q).getD 0))) := by
      interval_cases q <;> first | contradiction | decide
    obtain ⟨a, b, hla, hlb, hab⟩ := ha _ _ hmem hlookup
    have hqa : a = q := by
      apply Option.some.inj
      apply hla.symm.trans
      interval_cases q <;> decide
    have hqb : b = (badFaceClosed.get q).getD 0 := by
      apply Option.some.inj
      apply hlb.symm.trans
      interval_cases q <;> decide
    subst a; subst b
    convert hab using 1
    interval_cases q <;> decide
  have hadj : ρ.rotAdj = badFaceClosed.rotAdj := by
    apply Array.ext
    · exact hρ.size
    · intro q hq hqb
      have hb : q < 12 := by have := hρ.size; change ρ.rotAdj.size = 12 at this; omega
      simpa only [RotationSystem.get, Array.getElem?_eq_getElem hq,
        Array.getElem?_eq_getElem hqb, Option.bind_some, id_eq] using hgets q hb
  have heq : ρ = badFaceClosed := by cases ρ; congr
  have he := hρ.euler
  rw [heq] at he
  exact (by decide : ¬EulerFormula [(0, 1), (0, 1), (0, 1)] 2 badFaceClosed) he

theorem doubleRot_not_cofacial : ¬doubleRot.SameFaceOrbit 0 2 := by
  intro h
  have hpres : ∀ a b, a = 0 ∨ a = 6 → doubleRot.faceStep a = some b → b = 0 ∨ b = 6 := by
    intro a b ha hb
    rcases ha with rfl | rfl
    · have : b = 6 := Option.some.inj hb.symm
      exact Or.inr this
    · have : b = 0 := Option.some.inj hb.symm
      exact Or.inl this
  have hreach : ∀ b, doubleRot.SameFaceOrbit 0 b → b = 0 ∨ b = 6 := by
    intro b hb
    induction hb with
    | refl => exact Or.inl rfl
    | tail h hab ih => exact hpres _ _ ih hab
  have hh := hreach 2 h
  omega

end Spqr.PlanarEmbedCounterexample
