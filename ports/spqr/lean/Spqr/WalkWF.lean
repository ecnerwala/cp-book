import Spqr.Build
import Spqr.ItemSpec
import Spqr.RangesWF
import Spqr.RangesFinal

/-!
# Walk-phase contract

The walk-side canonicity admission behind the walk→relabel interface `Items.WF`
(`PROOF.md` §3–§4) is `walk_canonical`; `spqrTree_eq`
identifies `Graph.spqrTree` with the reference pipeline. `walk_items_wf` itself is assembled in
`Spqr.WalkItemsWF` (above `WalkCover`/`WalkTyping`, whose `walk_tree`/`walk_typing` it consumes);
this module stays below the st-order layer (`StSpec`/`StWalk`/`StOriented`), which only needs
`spqrTree_eq` and takes `Items.WF` as a hypothesis.
-/

namespace Spqr

/-- Admitted: canonicity of the unternarized walk (`PROOF.md` §4.6; for `tern = true` P under P /
S under S do occur, `check_ranges` seeds 0 and 386). -/
theorem walk_canonical (g : Graph) (vo eo : List Nat) :
    Items.Canonical (g.walk false (g.dfsForest vo eo)).items := by
  sorry

/-- `Graph.spqrTree` (the fast pipeline) is the reference walk followed by `relabelTree`. -/
theorem spqrTree_eq (g : Graph) (tern : Bool) (vo eo : List Nat) :
    g.spqrTree tern vo eo = relabelTree g (g.walk tern (g.dfsForest vo eo)).items := by
  simp [Graph.spqrTree, Graph.dfsForestFast_eq, Graph.walkFast_items, Fast.relabelTreeFast_eq]

end Spqr
