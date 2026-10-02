import Spqr.PlanarEmbedLeaf

open Spqr Spqr.PlanarSpqrTree

namespace Spqr.PlanarEmbedCounterexample

def edgeGraph : Graph := ⟨2, #[(0, 1)]⟩
def edgeTree : PlanarSpqrTree := {
  nv := 2, ne := 1, vertIndex := #[some 1, some 4], edgeIndex := #[some 2],
  edgeFlipped := #[false], par := #[none, some 0, some 1, some 2, some 2],
  subtreeEnd := #[5, 5, 5, 4, 5], types := #[.F, .V, .Q, .I, .V],
  origId := #[none, some 0, some 0, none, some 1], chBounds := #[0, 1, 2, 4, 4, 4],
  chDat := #[1, 2, 3, 4],
  nodeVerts := #[⟨0, 1⟩, ⟨2, 1⟩, ⟨2, 4⟩, ⟨3, 1⟩, ⟨3, 4⟩],
  nvBounds := #[0, 1, 1, 3, 5, 5], vertParNv := #[none, some 0, none, none, some 2],
  nodeEdges := #[⟨2, some 1, (1, 2)⟩, ⟨3, some 0, (3, 4)⟩],
  neBounds := #[0, 0, 0, 1, 2, 2], adjBounds := #[0, 0, 0, 0, 1, 2, 2, 2, 3, 4, 4],
  adjDat := #[⟨0, 2⟩, ⟨0, 1⟩, ⟨1, 4⟩, ⟨1, 3⟩], nodePlanar := #[true, true, true, true, true],
  neRotAdj := #[some 1, some 0, some 3, some 2, some 5, some 4, some 7, some 6] }

def badState : EmbedState :=
  ⟨#[some 1, some 0, none, none],
    #[#[none, none, none, none], #[none, none, some 2, some 3],
      #[none, none, none, none], #[none, none, none, none], #[none, none, none, none]]⟩

def edgeRot : RotationSystem := ⟨#[some 1, some 0, some 3, some 2]⟩

theorem edgeRot_planar : IsPlanarEmbedding [(0,1)] 2 edgeRot := by
  constructor
  · constructor
    · rfl
    · simp
    · intro q hq
      change q < 4 at hq
      interval_cases q <;> decide
    · intro q hq r hr
      change q < 4 at hq
      interval_cases q <;> simp [RotationSystem.get, edgeRot] at hr <;> subst hr <;> decide
    · intro q hq r hr
      change q < 4 at hq
      interval_cases q <;> simp [RotationSystem.get, edgeRot] at hr <;> subst hr <;> decide
    · intro q hq r hr
      change q < 4 at hq
      interval_cases q <;> simp [RotationSystem.get, edgeRot] at hr <;> subst hr <;> decide
    · decide
  · decide

theorem badState_before : edgeTree.GluedPieces edgeGraph 1 badState := by
  constructor
  · decide
  · decide
  · intro j hj q
    have : j = 0 := by omega
    subst j
    rintro ⟨k, hk⟩
    change (#[none, none, none, none] : Array (Option Nat))[k]? = some (some q) at hk
    have hk' : k < 4 := by
      by_contra hh
      rw [Array.getElem?_eq_none (by simpa using Nat.le_of_not_gt hh)] at hk
      cases hk
    interval_cases k <;> cases hk
  · intro q hq
    have hn := hq 1 (by omega) (by decide)
    change QE.edge q ∉ [0] at hn
    have hge : 4 ≤ q := by simp [QE.edge] at hn; omega
    exact Or.inl (Array.getElem?_eq_none (by exact hge))
  · intro j hj
    have hjs : j < 5 := hj.2.1
    have hj1 : j = 1 := by
      interval_cases j <;> simp_all [PlanarSpqrTree.Maximal, edgeTree]
    subst j
    refine ⟨edgeRot, by exact edgeRot_planar, ?_, ?_, ?_⟩
    · intro q r hq hr
      change QE.edge q ∈ [0] at hq
      have hql : q < 4 := by simp [QE.edge] at hq; omega
      interval_cases q <;> simp [badState] at hr <;> subst r
      · exact ⟨0, 1, by decide, by decide, by decide⟩
      · exact ⟨1, 0, by decide, by decide, by decide⟩
    · intro q hq
      change QE.edge q ∈ [0] at hq
      have hql : q < 4 := by simp [QE.edge] at hq; omega
      interval_cases q <;> simp only [badState, EmbedState.exposedAt]
      · constructor
        · intro h; cases h
        · rintro ⟨k, hk⟩
          have hm : some 0 ∈ [none, none, some 2, some 3] :=
            List.mem_of_getElem? (by simpa using hk)
          simp at hm
      · constructor
        · intro h; cases h
        · rintro ⟨k, hk⟩
          have hm : some 1 ∈ [none, none, some 2, some 3] :=
            List.mem_of_getElem? (by simpa using hk)
          simp at hm
      · exact ⟨fun _ => ⟨2, rfl⟩, fun _ => rfl⟩
      · exact ⟨fun _ => ⟨3, rfl⟩, fun _ => rfl⟩
    · intro k hk
      interval_cases k
      · simp [badState]
      · constructor
        · intro a ha
          have : a = 2 := by simpa [badState] using ha.symm
          subst a
          exact ⟨3, by rfl, 2, 3, by decide, by decide, by decide⟩
        · intro b hb
          have : b = 3 := by simpa [badState] using hb.symm
          subst b
          exact ⟨2, rfl⟩

theorem badState_after : ¬edgeTree.GluedPieces edgeGraph 0 ((edgeTree.embedItem 0).run badState).2 := by
  intro h
  obtain ⟨ρ, _, _, hu, _⟩ := h.piece 0 ⟨by omega, by decide, by intro p hp; cases hp⟩
  have he := (hu 2 (by change 0 ∈ [0]; simp)).1 (by decide)
  obtain ⟨k, hk⟩ := he
  change (#[none, none, none, none] : Array (Option Nat))[k]? = some (some 2) at hk
  have hk' : k < 4 := by
    by_contra hh
    rw [Array.getElem?_eq_none (by simpa using Nat.le_of_not_gt hh)] at hk
    cases hk
  interval_cases k <;> cases hk

theorem badState_excluded : ¬edgeTree.GluedUpTo edgeGraph 1 badState := by
  intro h
  have hh := (h.outer_slots 1 2 2 (by rfl)).2.2 0 (by rfl) (Or.inl rfl)
  omega

end Spqr.PlanarEmbedCounterexample
