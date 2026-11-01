import Spqr.Proofs.RSkelRoot
import Mathlib.Tactic.IntervalCases

/-!
# Kernel-checked counterexample: `(DfsData.ofForest (g.dfsForest vo eo)).Rooted g` fails

`DfsData.ofForest` takes `root` from the *first* tree of the forest, and `dfsForest` roots a tree
at every vertex in `inOrder` order — an isolated vertex ordered first (here `0` in the triple
edge `1 — 2`, a 2-connected graph) becomes `root` with no out-edges, so no edge endpoint is below
it. `Proofs/RSkelRoot.lean`'s `dfsForest_rooted` is therefore stated for `ofForest` re-rooted at
a depth-0 vertex (the root of the tree holding the edges), which is what the R proofs consume.
-/

namespace Spqr.RRootedCounter

def g : Graph := ⟨3, #[(1, 2), (1, 2), (1, 2)]⟩

def t₁ : DfsTree := .node 1 [.tree 0 .component
  (.node 2 [.back 1 1 (.ret 0 .backEdge), .back 2 1 (.ret 0 .backEdge)])]

theorem dfsForest_eq : g.dfsForest [] [] = [.node 0 [], t₁] := by cbv

theorem g_wf : g.WF := by
  intro p hp
  simp only [g, List.mem_cons, List.not_mem_nil, or_false, Array.mem_def] at hp
  rcases hp with rfl | rfl | rfl <;> decide

theorem edges_eq (e : Nat) (he : e < g.ne) : g.edges[e]? = some (1, 2) := by
  have he : e < 3 := he
  interval_cases e <;> rfl

theorem two_connected : g.TwoConnected := by
  intro v e e' he he'
  by_cases hv : v = 1
  · exact .inr ⟨2, 2, ⟨1, .inr (edges_eq e he)⟩, ⟨1, .inr (edges_eq e' he')⟩, .refl (by omega)⟩
  · exact .inr ⟨1, 1, ⟨2, .inl (edges_eq e he)⟩, ⟨2, .inl (edges_eq e' he')⟩, .refl (Ne.symm hv)⟩

theorem root_eq : (DfsData.ofForest (g.dfsForest [] [])).root = 0 := by rw [dfsForest_eq]; rfl

theorem outs_zero : (DfsData.ofForest (g.dfsForest [] [])).outs 0 = [] := by rw [dfsForest_eq]; rfl

theorem not_rooted : ¬ (DfsData.ofForest (g.dfsForest [] [])).Rooted g := by
  intro h
  have h1 := h 0 1 ⟨2, .inl rfl⟩
  rw [root_eq] at h1
  rcases Relation.ReflTransGen.cases_head h1 with h | ⟨c, ⟨o, ho, -, -⟩, -⟩
  · cases h
  · rw [outs_zero] at ho
    simp at ho

/-- The original statement of `dfsForest_rooted` is false. -/
theorem counter : ∃ (g : Graph) (vo eo : List Nat), g.WF ∧ OrderOK g.nv vo ∧ OrderOK g.ne eo ∧
    g.TwoConnected ∧ ¬ (DfsData.ofForest (g.dfsForest vo eo)).Rooted g :=
  ⟨g, [], [], g_wf, ⟨List.nodup_nil, by simp⟩, ⟨List.nodup_nil, by simp⟩, two_connected, not_rooted⟩

#print axioms counter

end Spqr.RRootedCounter
