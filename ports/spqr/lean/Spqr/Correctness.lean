import Spqr.Build
import Spqr.Spec
import Spqr.ItemSpec
import Spqr.Proofs.Dfs
import Spqr.RelabelRep
import Spqr.WalkWF
import Spqr.RelabelWF
import Spqr.StOriented

/-!
# Correctness theorems

Top level: `spqrTree_wf` and `spqrTree_represents`. They are assembled from one theorem per
phase; each phase theorem is proven in its own module (`Spqr.Proofs.*`; `walk_items_wf`/`spqrTree_eq` in
`Spqr.WalkWF`, `walk_items_rOriented'` in `Spqr.StOriented`). Statements below that are still
`sorry` are listed in the README.
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

/-- Phase 3: relabelling a well-formed, R-oriented item tree gives a well-formed output ...
(`RelabelWF.lean`, modulo `relabel_node_spec`). -/
theorem relabelTree_wf (g : Graph) (items : Items) (h : items.WF g) (hor : items.ROriented g) :
    (relabelTree g items).WF := by
  obtain ⟨idx, -, hidx, hnode⟩ := relabel_node_spec g items h
  exact RelabelAll.wf_tree ⟨h, hidx, hnode⟩ hor

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

theorem spqrTree_wf (g : Graph) (tern : Bool) (vo eo : List Nat) : (g.spqrTree tern vo eo).WF := by
  rw [spqrTree_eq]; exact relabelTree_wf g _ (walk_items_wf g tern vo eo) (walk_items_rOriented' g tern vo eo)

theorem spqrTree_represents (g : Graph) (tern : Bool) (vo eo : List Nat) :
    (g.spqrTree tern vo eo).Represents g := by
  have h := spqrTree_r_three_connected g tern vo eo
  rw [spqrTree_eq] at h ⊢
  exact relabelTree_represents_of_r g _ (walk_items_wf g tern vo eo) h

end Spqr
