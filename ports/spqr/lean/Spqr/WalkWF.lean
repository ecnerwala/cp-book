import Spqr.Build
import Spqr.ItemSpec
import Spqr.RangesWF

/-!
# Walk-phase contract

`walk_items_wf` is the walk→relabel interface (`Items.WF`, `PROOF.md` §3–§4) and `spqrTree_eq`
identifies `Graph.spqrTree` with the reference pipeline. They live below `Spqr.Correctness` so that
the st-order layer (`StSpec`/`StWalk`/`StOriented`) can use them and `Correctness` can import that
layer.
-/

namespace Spqr

/-- Admitted (hypothesis-free bridge): `Items.Tree` for the walk. Under `g.WF`/`OrderOK` this is
`WalkCover.walk_tree` (which lives above this module and is itself modulo `walk_sides`). -/
theorem walk_items_tree (g : Graph) (tern : Bool) (vo eo : List Nat) :
    Items.Tree g (g.walk tern (g.dfsForest vo eo)).items := by
  sorry

/-- Admitted (hypothesis-free bridge): the typing facts of `WalkTyping.walk_typing`
(`i_o_leaf`, `vs_shape`, `vs_lt`), which needs `0 < g.nv`, bounded trees and edge coverage. -/
theorem walk_typingFacts (g : Graph) (tern : Bool) (vo eo : List Nat) :
    Items.TypingFacts g (g.walk tern (g.dfsForest vo eo)).items := by
  sorry

/-- Admitted: the walk's items have the postorder range / attachment structure (`PROOF.md` §4.6),
relative to the DFS edge postorder. Checked empirically by `check_ranges`. -/
theorem walk_ranges (g : Graph) (tern : Bool) (vo eo : List Nat) :
    Items.Ranges g (g.walk tern (g.dfsForest vo eo)).items (edgePostorderForest (g.dfsForest vo eo)) := by
  sorry

/-- Admitted, and FALSE for `tern = true`: ternarizing a P or S node leaves P under P / S under S
(`check_ranges`, seeds 0 and 386; `PROOF.md` §4.6). `Spec.Represents.canonical` already says "unless
ternarizing"; `Items.Shapes.canonical` and this statement need the `tern = false` guard. -/
theorem walk_canonical (g : Graph) (tern : Bool) (vo eo : List Nat) :
    ∀ p c, Items.IsParent (g.walk tern (g.dfsForest vo eo)).items p c →
      (Items.type (g.walk tern (g.dfsForest vo eo)).items c = .S →
        Items.type (g.walk tern (g.dfsForest vo eo)).items p ≠ .S) ∧
      (Items.type (g.walk tern (g.dfsForest vo eo)).items c = .P →
        Items.type (g.walk tern (g.dfsForest vo eo)).items p ≠ .P) := by
  sorry

/-- Phase 2: the walk's items satisfy the item-level specification. `Endpoints`/`Shapes` are derived
from `walk_ranges` by `Items.wf_of_ranges` (pure item-level reasoning, no tstack facts). -/
theorem walk_items_wf (g : Graph) (tern : Bool) (vo eo : List Nat) :
    Items.WF g (g.walk tern (g.dfsForest vo eo)).items :=
  Items.wf_of_ranges (walk_items_tree g tern vo eo) (walk_typingFacts g tern vo eo)
    (walk_ranges g tern vo eo) (walk_canonical g tern vo eo)

theorem spqrTree_eq (g : Graph) (tern : Bool) (vo eo : List Nat) :
    g.spqrTree tern vo eo = relabelTree g (g.walk tern (g.dfsForest vo eo)).items := by
  simp [Graph.spqrTree, Graph.dfsForestFast_eq, Graph.walkFast_items, Fast.relabelTreeFast_eq]

end Spqr
