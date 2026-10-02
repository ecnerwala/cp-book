import Spqr.Build
import Spqr.Spec
import Spqr.ItemSpec
import Spqr.Proofs.Dfs
import Spqr.RelabelRep

/-!
# Correctness theorems

Top level: `spqrTree_wf` and `spqrTree_represents`. They are assembled from one theorem per
phase; each phase theorem is proven in its own module (`Spqr.Proofs.*`). Statements below that
are still `sorry` are listed in the README.
-/

namespace Spqr

/-- Phase 1: the DFS forest is a spanning forest of `g` in which every edge appears exactly once
as an out-edge (tree edge from the parent, or back edge from the deeper endpoint / loop). -/
theorem dfsForest_spanning (g : Graph) (hg : g.WF) (vo eo : List Nat) (hvo : OrderOK g.nv vo)
    (heo : OrderOK g.ne eo) :
    let forest := g.dfsForest vo eo
    -- every vertex appears exactly once
    ((forest.flatMap DfsTree.verts).Perm (List.range g.nv)) ∧
    -- every edge appears exactly once
    ((forest.flatMap DfsTree.edges).Perm (List.range g.ne)) :=
  dfsForest_spanning' hg hvo heo

/-- Phase 2: the walk's items satisfy the item-level specification. -/
theorem walk_items_wf (g : Graph) (tern : Bool) (vo eo : List Nat) :
    Items.WF g (g.walk tern (g.dfsForest vo eo)).items := by
  sorry

/-- Phase 3: relabelling a well-formed item tree gives a well-formed output ... -/
theorem relabelTree_wf (g : Graph) (items : Items) (h : items.WF g) :
    (relabelTree g items).WF := by
  sorry

/-- ... that represents `g`. The one hypothesis beyond `Items.WF` is the item-level 3-connectivity
of the R items' skeletons (`Items.RThreeConnected`, `RelabelRep.lean`), the walk-side statement of
`items_r_three_connected` (`Spqr.Proofs.RItems`); transported, not implied by `Items.Shapes`. -/
theorem relabelTree_represents (g : Graph) (items : Items) (h : items.WF g)
    (hR : Items.RThreeConnected g items) :
    (relabelTree g items).Represents g :=
  relabelTree_represents' g items h hR

theorem spqrTree_r_three_connected (g : Graph) (tern : Bool) (vo eo : List Nat) :
    let t := g.spqrTree tern vo eo
    ∀ i, i < t.size → t.type i = .R →
      SpqrTree.ThreeConnected (t.nVerts i)
        ((t.skeleton i).map fun p => (p.1 - (t.nvRange i).1, p.2 - (t.nvRange i).1)) := by
  sorry

theorem spqrTree_eq (g : Graph) (tern : Bool) (vo eo : List Nat) :
    g.spqrTree tern vo eo = relabelTree g (g.walk tern (g.dfsForest vo eo)).items := by
  simp [Graph.spqrTree, Graph.dfsForestFast_eq, Graph.walkFast_items, Fast.relabelTreeFast_eq]

theorem spqrTree_wf (g : Graph) (tern : Bool) (vo eo : List Nat) : (g.spqrTree tern vo eo).WF := by
  rw [spqrTree_eq]; exact relabelTree_wf g _ (walk_items_wf g tern vo eo)

theorem spqrTree_represents (g : Graph) (tern : Bool) (vo eo : List Nat) :
    (g.spqrTree tern vo eo).Represents g := by
  have h := spqrTree_r_three_connected g tern vo eo
  rw [spqrTree_eq] at h ⊢
  exact relabelTree_represents_of_r g _ (walk_items_wf g tern vo eo) h

end Spqr
