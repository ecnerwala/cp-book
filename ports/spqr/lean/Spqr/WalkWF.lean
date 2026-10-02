import Spqr.Build
import Spqr.ItemSpec

/-!
# Walk-phase contract

`walk_items_wf` is the walk→relabel interface (`Items.WF`, `PROOF.md` §3–§4) and `spqrTree_eq`
identifies `Graph.spqrTree` with the reference pipeline. They live below `Spqr.Correctness` so that
the st-order layer (`StSpec`/`StWalk`/`StOriented`) can use them and `Correctness` can import that
layer.
-/

namespace Spqr

/-- Phase 2: the walk's items satisfy the item-level specification. -/
theorem walk_items_wf (g : Graph) (tern : Bool) (vo eo : List Nat) :
    Items.WF g (g.walk tern (g.dfsForest vo eo)).items := by
  sorry

theorem spqrTree_eq (g : Graph) (tern : Bool) (vo eo : List Nat) :
    g.spqrTree tern vo eo = relabelTree g (g.walk tern (g.dfsForest vo eo)).items := by
  simp [Graph.spqrTree, Graph.dfsForestFast_eq, Graph.walkFast_items, Fast.relabelTreeFast_eq]

end Spqr
