import Spqr.Proofs.PlanarEmbedCounterexample

open Spqr Spqr.PlanarSpqrTree

namespace Spqr.PlanarEmbedCounterexample

def starGraph : Graph := ⟨3, #[(0, 1), (0, 2)]⟩

def starTree : PlanarSpqrTree := {
  nv := 3, ne := 2, vertIndex := #[some 1, some 4, some 7], edgeIndex := #[some 2, some 5],
  edgeFlipped := #[false, false], par := #[none, some 0, some 1, some 2, some 2, some 1, some 5, some 5],
  subtreeEnd := #[8, 8, 5, 4, 5, 8, 7, 8], types := #[.F, .V, .Q, .I, .V, .Q, .I, .V],
  origId := #[none, some 0, some 0, none, some 1, some 1, none, some 2],
  chBounds := #[0, 1, 3, 5, 5, 5, 7, 7, 7], chDat := #[1, 2, 5, 3, 4, 6, 7],
  nodeVerts := #[⟨0, 1⟩, ⟨2, 1⟩, ⟨2, 4⟩, ⟨3, 1⟩, ⟨3, 4⟩, ⟨5, 1⟩, ⟨5, 7⟩, ⟨6, 1⟩, ⟨6, 7⟩],
  nvBounds := #[0, 1, 1, 3, 5, 5, 7, 9, 9],
  vertParNv := #[none, some 0, none, none, some 2, none, none, some 6],
  nodeEdges := #[⟨2, some 1, (1, 2)⟩, ⟨3, some 0, (3, 4)⟩, ⟨5, some 3, (5, 6)⟩, ⟨6, some 2, (7, 8)⟩],
  neBounds := #[0, 0, 0, 1, 2, 2, 3, 4, 4],
  adjBounds := #[0, 0, 0, 0, 1, 2, 2, 2, 3, 4, 4, 4, 5, 6, 6, 6, 7, 8, 8],
  adjDat := #[⟨0, 2⟩, ⟨0, 1⟩, ⟨1, 4⟩, ⟨1, 3⟩, ⟨2, 6⟩, ⟨2, 5⟩, ⟨3, 8⟩, ⟨3, 7⟩],
  nodePlanar := #[true, true, true, true, true, true, true, true],
  neRotAdj := #[some 1, some 0, some 3, some 2, some 5, some 4, some 7, some 6,
    some 9, some 8, some 11, some 10, some 13, some 12, some 15, some 14] }

def badVState : EmbedState :=
  ⟨#[some 1, some 0, none, none, some 5, some 4, none, none],
    #[#[none, none, none, none], #[none, none, none, none], #[some 2, some 3, none, none],
      #[none, none, none, none], #[none, none, none, none], #[some 6, some 7, none, none],
      #[none, none, none, none], #[none, none, none, none]]⟩

theorem starEdgeRot_planar (e : Nat) (he : e < 2) :
    IsPlanarEmbedding [starGraph.edges[e]!] 3 edgeRot := by
  constructor
  · constructor
    · rfl
    · interval_cases e <;> simp [starGraph]
    · exact edgeRot_planar.total
    · exact edgeRot_planar.involution
    · exact edgeRot_planar.opposite_dir
    · intro q hq r hr
      change q < 4 at hq
      interval_cases q <;> simp [RotationSystem.get, edgeRot] at hr <;> subst hr <;>
        interval_cases e <;> decide
    · interval_cases e <;> decide
  · interval_cases e <;> decide

theorem badVState_before : starTree.GluedSlots starGraph 2 badVState := by
  refine ⟨⟨by decide, by decide, ?_, ?_, ?_⟩, ?_, ?_⟩
  · intro j hj q
    interval_cases j <;> rintro ⟨k, hk⟩
    all_goals
      have hm : some q ∈ [none, none, none, none] := List.mem_of_getElem? (by simpa [badVState] using hk)
      simp at hm
  · intro q hq
    have h0 := hq 2 (by omega) (by decide)
    have h1 := hq 5 (by omega) (by decide)
    change QE.edge q ∉ [0] at h0
    change QE.edge q ∉ [1] at h1
    have hge : 8 ≤ q := by simp [QE.edge] at h0 h1; omega
    exact Or.inl (Array.getElem?_eq_none hge)
  · intro j hj
    have hjs : j < 8 := hj.2.1
    have hj' : j = 2 ∨ j = 5 := by
      interval_cases j <;> simp_all [PlanarSpqrTree.Maximal, starTree]
    rcases hj' with rfl | rfl
    all_goals
      refine ⟨edgeRot, ?_, ?_, ?_, ?_⟩
    case inl.refine_1 => exact starEdgeRot_planar 0 (by decide)
    case inr.refine_1 => exact starEdgeRot_planar 1 (by decide)
    case inl.refine_2 | inr.refine_2 =>
        intro q r hq hr
        first | change q / 4 ∈ [0] at hq | change q / 4 ∈ [1] at hq
        simp only [List.mem_singleton] at hq
        have hql : q < 8 := by omega
        interval_cases q <;> norm_num at hq <;> simp [badVState] at hr <;> subst r
        all_goals first
          | exact ⟨0, 1, by decide, by decide, by decide⟩
          | exact ⟨1, 0, by decide, by decide, by decide⟩
    case inl.refine_3 | inr.refine_3 =>
        intro q hq
        first | change q / 4 ∈ [0] at hq | change q / 4 ∈ [1] at hq
        simp only [List.mem_singleton] at hq
        have hql : q < 8 := by omega
        interval_cases q <;> norm_num at hq
        all_goals
          simp only [badVState, EmbedState.exposedAt]
          change _ ↔ ∃ k, _[k]? = some (some _)
          rw [← Array.mem_iff_getElem?]
          simp
    case inl.refine_4 | inr.refine_4 =>
        intro k hk
        interval_cases k <;> simp [badVState]
        all_goals exact ⟨2, by decide, 3, by decide, by decide⟩
  · intro j hj
    change j < 8 at hj
    interval_cases j <;> decide
  · intro j k q hk
    have hj : j < 8 := by
      by_contra hh
      have hn : badVState.outerE[j]? = none := Array.getElem?_eq_none (by exact Nat.le_of_not_gt hh)
      rw [hn] at hk
      cases hk
    interval_cases j <;> simp [badVState] at hk
    all_goals
      have hkl : k < 4 := by
        by_contra hh
        rw [List.getElem?_eq_none (by simpa using Nat.le_of_not_gt hh)] at hk
        cases hk
      interval_cases k <;> simp [SpqrTree.type, SpqrTree.parent, starTree] at hk ⊢

theorem badVState_after :
    ¬starTree.GluedSlots starGraph 1 ((starTree.embedItem 1).run badVState).2 := by
  intro h
  obtain ⟨ρ, hρ, ha, _, _⟩ := h.piece 1 ⟨by decide, by decide, by
    intro p hp
    have : p = 0 := by simpa [starTree] using hp.symm
    omega⟩
  obtain ⟨a, b, hla, hlb, hab⟩ := ha 3 6 (by change 0 ∈ [0,1]; decide) (by decide)
  have ha : a = 3 := Option.some.inj (hla.symm.trans (by decide))
  have hb : b = 6 := Option.some.inj (hlb.symm.trans (by decide))
  subst a; subst b
  have hv := hρ.same_vertex 3 (by rw [hρ.size]; decide) 6 hab
  have : (some 1 : Option Nat) = some 2 := hv
  cases this

theorem badVState_excluded : ¬starTree.GluedUpTo starGraph 2 badVState := by
  intro h
  have hv := h.outer_at_vertex 2 1 0 2 (by rfl) (by rfl) (by rfl) ⟨0, rfl⟩
  have : (some 1 : Option Nat) = some 0 := hv
  cases this

end Spqr.PlanarEmbedCounterexample
