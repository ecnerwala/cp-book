import Spqr.Proofs.RunSaturation
import Spqr.ItemAcyc

namespace Spqr.REdgeDomainCheck

def k4 : Graph := ⟨4, #[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]⟩

def state : WalkState := k4.walk false (k4.dfsForest [] [])

theorem self_edge : ¬6 < k4.ne ∧ Items.EdgeBelow k4 state.items 11 6 := by
  exact ⟨by decide, Relation.ReflTransGen.refl⟩

theorem child_leaves : ∀ c ∈ state.rChildren 11, c < 11 ∧ Items.ch state.items c = [] := by
  have hc : state.rChildren 11 = [7, 9, 10, 6, 8] := by cbv
  intro c hc'
  rw [hc] at hc'
  simp only [List.mem_cons, List.not_mem_nil, or_false] at hc'
  rcases hc' with rfl | rfl | rfl | rfl | rfl <;> cbv

theorem not_raw_coverage : ¬(∀ e, Items.EdgeBelow k4 state.items 11 e ↔
    runEdges (Items.EdgeBelow k4 state.items) (state.rChildren 11) e) := by
  intro h
  obtain ⟨c, hc, he⟩ := (h 6).1 self_edge.2
  have heq := Items.below_eq_of_ch_nil (child_leaves c hc).2 he
  have hlt := (child_leaves c hc).1
  change 11 = c at heq
  exact Nat.not_lt_of_ge (Nat.le_of_eq heq) hlt

theorem bounded_coverage : ∀ e, e < k4.ne → (Items.EdgeBelow k4 state.items 11 e ↔
    runEdges (Items.EdgeBelow k4 state.items) (state.rChildren 11) e) := by
  have hch : Items.ch state.items 11 = [7, 4, 9, 10, 6, 3, 8] := by cbv
  have hr : state.rChildren 11 = [7, 9, 10, 6, 8] := by cbv
  have hleaf : ∀ c ∈ Items.ch state.items 11, Items.ch state.items c = [] := by
    intro c hc
    rw [hch] at hc
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hc
    rcases hc with rfl | rfl | rfl | rfl | rfl | rfl | rfl <;> cbv
  intro e he
  constructor
  · intro hb
    rcases Relation.ReflTransGen.cases_head hb with h | ⟨c, hc, hb⟩
    · change 11 = 1 + 4 + e at h
      change e < 6 at he
      omega
    · have heq := Items.below_eq_of_ch_nil (hleaf c hc) hb
      change 1 + 4 + e = c at heq
      refine ⟨c, ?_, hb⟩
      change c ∈ Items.ch state.items 11 at hc
      rw [hch] at hc
      rw [hr]
      simp only [List.mem_cons, List.not_mem_nil, or_false] at hc ⊢
      rcases hc with rfl | rfl | rfl | rfl | rfl | rfl | rfl <;> first | decide | omega
  · rintro ⟨c, hc, hb⟩
    exact (Relation.ReflTransGen.single (List.mem_of_mem_filter hc)).trans hb

#print axioms self_edge
#print axioms child_leaves
#print axioms not_raw_coverage
#print axioms bounded_coverage

end Spqr.REdgeDomainCheck
