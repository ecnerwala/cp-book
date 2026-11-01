import Spqr.Proofs.PlanarEmbedQCounterexample

open Spqr Spqr.PlanarSpqrTree

namespace Spqr.PlanarEmbedCounterexample

def parallelGraph : Graph := ⟨2, #[(0, 1), (0, 1)]⟩
def parallelTree : PlanarSpqrTree := { edgeTree with
  ne := 2, edgeIndex := #[some 2, some 3], edgeFlipped := #[false, false]
  types := #[.F, .V, .Q, .Q, .V], origId := #[none, some 0, some 0, some 1, some 1] }

def badCapState : EmbedState :=
  ⟨Array.replicate 8 none,
    #[#[none, none, none, none], #[none, none, none, none], #[none, none, none, none],
      #[some 6, some 7, some 4, some 5], #[none, none, none, none]]⟩

set_option maxHeartbeats 1000000 in
theorem badCapState_before : parallelTree.GluedVertex parallelGraph 3 badCapState := by
  have hslots : ∀ j k q, badCapState.outerE[j]?.bind (fun o => o[k]?) = some (some q) →
      j = 3 ∧ ((k = 0 ∧ q = 6) ∨ (k = 1 ∧ q = 7) ∨ (k = 2 ∧ q = 4) ∨ (k = 3 ∧ q = 5)) := by
    intro j k q hk
    have hj : j < 5 := by
      by_contra hn
      rw [Array.getElem?_eq_none (by exact Nat.le_of_not_gt hn)] at hk
      cases hk
    interval_cases j <;> simp [badCapState] at hk
    all_goals
      have hkl : k < 4 := by
        by_contra hn
        rw [List.getElem?_eq_none (by simpa using Nat.le_of_not_gt hn)] at hk
        cases hk
      interval_cases k <;> simp_all
  refine ⟨⟨⟨⟨⟨by decide, by decide, ?_, ?_, ?_⟩, ?_, ?_⟩, ?_⟩, ?_, ?_⟩, ?_⟩
  · intro j hj q ⟨k, hk⟩
    have := (hslots j k q hk).1
    omega
  · intro q _
    by_cases hq : q < 8 <;> simp [badCapState, hq]
  · intro j hj
    have hjs : j < 5 := hj.2.1
    have hjs' : j = 3 ∨ j = 4 := by
      interval_cases j <;> simp_all [PlanarSpqrTree.Maximal, parallelTree, edgeTree]
    rcases hjs' with rfl | rfl
    · refine ⟨edgeRot, edgeRot_planar, ?_, ?_, ?_⟩
      · intro q r hq hr
        by_cases hq' : q < 8 <;> simp [badCapState, hq'] at hr
      · intro q hq
        change q / 4 ∈ [1] at hq
        simp only [List.mem_singleton] at hq
        have hq' : q < 8 := by omega
        interval_cases q <;> norm_num at hq
        all_goals
          simp only [badCapState, EmbedState.exposedAt]
          change _ ↔ ∃ k, _[k]? = some (some _)
          rw [← Array.mem_iff_getElem?]
          simp
      · intro k hk
        interval_cases k <;> simp [badCapState]
        · exact ⟨2, by decide, 3, by decide, by decide⟩
        · exact ⟨0, by decide, 1, by decide, by decide⟩
    · refine ⟨⟨#[]⟩, isPlanarEmbedding_nil 2, ?_, ?_, ?_⟩
      · intro q r hq; change QE.edge q ∈ [] at hq; cases hq
      · intro q hq; change QE.edge q ∈ [] at hq; cases hq
      · intro k hk; interval_cases k <;> simp [badCapState]
  · intro j hj
    change j < 5 at hj
    interval_cases j <;> decide
  · intro j k q hk
    obtain ⟨rfl, hs⟩ := hslots j k q hk
    have hk' : k < 4 := by rcases hs with h | h | h | h <;> omega
    refine ⟨hk', by decide, ?_⟩
    intro p hp hpt
    have hp' : p = 2 := by simpa [SpqrTree.parent, parallelTree, edgeTree] using hp.symm
    subst p
    cases hpt <;> contradiction
  · intro j p v q hp hpt hv ⟨k, hk⟩
    obtain ⟨rfl, _⟩ := hslots j k q hk
    have hp' : p = 2 := by simpa [SpqrTree.parent, parallelTree, edgeTree] using hp.symm
    subst p
    cases hpt
  · intro j k q hk
    obtain ⟨_, hs⟩ := hslots j k q hk
    rcases hs with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ <;> decide
  · intro j hj hjs p hp hpt _
    change j < 5 at hjs
    interval_cases j <;> simp [SpqrTree.parent, parallelTree, edgeTree] at hp <;> subst p <;>
      simp_all [SpqrTree.type, parallelTree, edgeTree]
  · intro j hj hjt v hv
    refine ⟨?_, ?_, ?_⟩
    · intro q ⟨k, hk⟩
      obtain ⟨rfl, _⟩ := hslots j k q hk
      cases hjt
    · intro k q hk
      obtain ⟨rfl, _⟩ := hslots j k q hk
      cases hjt
    · intro hge hne
      change j < 5 at hj
      interval_cases j <;> simp [SpqrTree.type, parallelTree, edgeTree] at hjt
      exact (hne (by rfl)).elim

theorem badCapState_after :
    ¬parallelTree.GluedPieces parallelGraph 2 ((parallelTree.embedItem 2).run badCapState).2 := by
  intro h
  obtain ⟨ρ, hρ, ha, _, _⟩ := h.piece 2 ⟨by decide, by decide, by
    intro p hp
    have : p = 1 := by simpa [parallelTree, edgeTree] using hp.symm
    omega⟩
  obtain ⟨a, b, hla, hlb, hab⟩ := ha 1 6 (by change 0 ∈ [0,1]; decide) (by decide)
  have ha : a = 1 := Option.some.inj (hla.symm.trans (by decide))
  have hb : b = 6 := Option.some.inj (hlb.symm.trans (by decide))
  subst a; subst b
  have hv := hρ.same_vertex 1 (by rw [hρ.size]; decide) 6 hab
  have : (some 0 : Option Nat) = some 1 := hv
  cases this

theorem badCapState_excluded : ¬parallelTree.GluedUpTo parallelGraph 3 badCapState := by
  intro h
  have hv := (h.outer_cap 3 (by decide) 1 (by decide) (0, 1) (by decide)).1 0 6 rfl
  have : (some 1 : Option Nat) = some 0 := hv
  cases this

end Spqr.PlanarEmbedCounterexample
