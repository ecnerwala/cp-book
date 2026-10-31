import Spqr.Proofs.Dfs
import Spqr.RelabelRep
import Spqr.Walk

/-!
# R items of the walk are 3-connected (PROOF.md §4.5, item level, output contract)

`Items.RThreeConnected g items` (`Spqr.RelabelRep`) is the item-level contract the relabel phase
transports to the output tree (`RelabelOK.r_three_connected`). This file states it for the walk's
items; `Spqr.Correctness.spqrTree_r_three_connected` is its transport. On a block it is
`walk_items_rThreeConnected_of_twoConnected` (`Spqr.Proofs.RItems`).
-/

namespace Spqr

/-- Admitted (PROOF.md §4.5): every R item of the walk's output has a 3-connected `rSkeleton`
(the cut form, `SpqrTree.ThreeConnected (items.nvList g i).length (items.rSkeleton g i)`). On a
block this is `walk_items_rThreeConnected_of_twoConnected` (`Spqr.Proofs.RItems`, from
`items_r_three_connected` through `Items.RSkel3.rThreeConnected_of_block`). For a general
`g.WF` input the whole-graph complement of `Items.RSkel3` is not the parent piece of an R item of a
non-trivial block, so the bridge needs the block-local complement (edges of the R item's block
outside it); that restatement is the remaining obligation. -/
theorem walk_items_rThreeConnected (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    Items.RThreeConnected g (g.walk tern (g.dfsForest vo eo)).items := by
  sorry

end Spqr
