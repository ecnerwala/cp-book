import Spqr.Proofs.PlanarEmbedVCounterexample

open Spqr Spqr.PlanarSpqrTree

namespace Spqr.PlanarEmbedCounterexample

def pathGraph : Graph := ⟨3, #[(0, 1), (1, 2)]⟩

def pathTree : PlanarSpqrTree := { starTree with
  par := #[none, some 0, some 1, some 2, some 2, some 4, some 5, some 5]
  subtreeEnd := #[8, 8, 8, 4, 8, 8, 7, 8]
  chBounds := #[0, 1, 2, 4, 4, 5, 7, 7, 7]
  chDat := #[1, 2, 3, 4, 5, 6, 7]
  nodeVerts := #[⟨0, 1⟩, ⟨2, 1⟩, ⟨2, 4⟩, ⟨3, 1⟩, ⟨3, 4⟩, ⟨5, 4⟩, ⟨5, 7⟩, ⟨6, 4⟩, ⟨6, 7⟩] }

def badQState : EmbedState :=
  ⟨#[none, none, none, none, some 5, some 4, none, none],
    #[#[none, none, none, none], #[none, none, none, none], #[none, none, none, none],
      #[none, none, none, none], #[some 6, some 7, none, none], #[some 4, some 5, none, none],
      #[none, none, none, none], #[none, none, none, none]]⟩

theorem pathEdgeRot_planar : IsPlanarEmbedding [(1, 2)] 3 edgeRot := by
  constructor
  · constructor
    · rfl
    · decide
    · exact edgeRot_planar.total
    · exact edgeRot_planar.involution
    · exact edgeRot_planar.opposite_dir
    · intro q hq r hr
      change q < 4 at hq
      interval_cases q <;> simp [RotationSystem.get, edgeRot] at hr <;> subst hr <;> decide
    · decide
  · decide

theorem badQState_before : pathTree.GluedOriented pathGraph 3 badQState := by
  refine ⟨⟨⟨⟨by decide, by decide, ?_, ?_, ?_⟩, ?_, ?_⟩, ?_⟩, ?_, ?_⟩
  · intro j hj q
    interval_cases j <;> rintro ⟨k, hk⟩
    all_goals
      have hm : some q ∈ [none, none, none, none] := List.mem_of_getElem? (by simpa [badQState] using hk)
      simp at hm
  · intro q hq
    have h1 := hq 4 (by omega) (by decide)
    change QE.edge q ∉ [1] at h1
    by_cases hq' : q < 8
    · simp only [QE.edge, List.mem_singleton] at h1
      interval_cases q <;> norm_num at h1 <;> exact Or.inr rfl
    · exact Or.inl (Array.getElem?_eq_none (by change 8 ≤ q; omega))
  · intro j hj
    have hjs : j < 8 := hj.2.1
    have hjs' : j = 3 ∨ j = 4 := by
      interval_cases j <;> simp_all [PlanarSpqrTree.Maximal, pathTree, starTree]
    rcases hjs' with rfl | rfl
    · refine ⟨⟨#[]⟩, isPlanarEmbedding_nil 3, ?_, ?_, ?_⟩
      · intro q r hq; change QE.edge q ∈ [] at hq; cases hq
      · intro q hq; change QE.edge q ∈ [] at hq; cases hq
      · intro k hk; interval_cases k <;> simp [badQState]
    · refine ⟨edgeRot, pathEdgeRot_planar, ?_, ?_, ?_⟩
      · intro q r hq hr
        change q / 4 ∈ [1] at hq
        simp only [List.mem_singleton] at hq
        have hq' : q < 8 := by omega
        interval_cases q <;> norm_num at hq <;> simp [badQState] at hr <;> subst r
        · exact ⟨0, 1, by decide, by decide, by decide⟩
        · exact ⟨1, 0, by decide, by decide, by decide⟩
      · intro q hq
        change q / 4 ∈ [1] at hq
        simp only [List.mem_singleton] at hq
        have hq' : q < 8 := by omega
        interval_cases q <;> norm_num at hq
        all_goals
          simp only [badQState, EmbedState.exposedAt]
          change _ ↔ ∃ k, _[k]? = some (some _)
          rw [← Array.mem_iff_getElem?]
          simp
      · intro k hk
        interval_cases k <;> simp [badQState]
        exact ⟨2, by decide, 3, by decide, by decide⟩
  · intro j hj
    change j < 8 at hj
    interval_cases j <;> decide
  · intro j k q hk
    have hj : j < 8 := by
      by_contra hn
      rw [Array.getElem?_eq_none (by exact Nat.le_of_not_gt hn)] at hk
      cases hk
    interval_cases j <;> simp [badQState] at hk
    all_goals
      have hkl : k < 4 := by
        by_contra hn
        rw [List.getElem?_eq_none (by simpa using Nat.le_of_not_gt hn)] at hk
        cases hk
      interval_cases k <;> simp [SpqrTree.type, SpqrTree.parent, pathTree, starTree] at hk ⊢
  · intro j p v q hp hpt hv ⟨k, hk⟩
    have hj : j < 8 := by
      by_contra hn
      rw [Array.getElem?_eq_none (by exact Nat.le_of_not_gt hn)] at hk
      cases hk
    interval_cases j <;> simp [badQState] at hk
    all_goals
      have hkl : k < 4 := by
        by_contra hn
        rw [List.getElem?_eq_none (by simpa using Nat.le_of_not_gt hn)] at hk
        cases hk
      interval_cases k <;> simp at hk <;> subst q
      all_goals
        simp [SpqrTree.parent, pathTree, starTree] at hp
        subst p
        simp [SpqrTree.type, pathTree, starTree] at hpt
        all_goals
          simp [pathTree, starTree] at hv
          subst v
          rfl
  · intro j k q hk
    have hj : j < 8 := by
      by_contra hn
      rw [Array.getElem?_eq_none (by exact Nat.le_of_not_gt hn)] at hk
      cases hk
    interval_cases j <;> simp [badQState] at hk
    all_goals
      have hkl : k < 4 := by
        by_contra hn
        rw [List.getElem?_eq_none (by simpa using Nat.le_of_not_gt hn)] at hk
        cases hk
      interval_cases k <;> simp_all <;> omega
  · intro j hj hjs p hp hpt _
    change j < 8 at hjs
    interval_cases j <;> simp [SpqrTree.parent, pathTree, starTree] at hp <;> subst p <;>
      simp_all [SpqrTree.type, pathTree, starTree]
    exact ⟨4, 0, rfl⟩

theorem badQState_after :
    ¬pathTree.GluedPieces pathGraph 2 ((pathTree.embedItem 2).run badQState).2 := by
  intro h
  obtain ⟨ρ, hρ, ha, _, _⟩ := h.piece 2 ⟨by decide, by decide, by
    intro p hp
    have : p = 1 := by simpa [pathTree, starTree] using hp.symm
    omega⟩
  obtain ⟨a, b, hla, hlb, hab⟩ := ha 3 6 (by change 0 ∈ [0,1]; decide) (by decide)
  have ha : a = 3 := Option.some.inj (hla.symm.trans (by decide))
  have hb : b = 6 := Option.some.inj (hlb.symm.trans (by decide))
  subst a; subst b
  have hv := hρ.same_vertex 3 (by rw [hρ.size]; decide) 6 hab
  have : (some 1 : Option Nat) = some 2 := hv
  cases this

theorem badQState_excluded : ¬pathTree.GluedUpTo pathGraph 3 badQState := by
  intro h
  have hv := (h.outer_vertex 4 (by decide) (by decide) 1 (by decide)).1 6 ⟨0, rfl⟩
  have : (some 2 : Option Nat) = some 1 := hv
  cases this

end Spqr.PlanarEmbedCounterexample
