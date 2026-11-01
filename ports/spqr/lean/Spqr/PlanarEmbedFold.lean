import Spqr.PlanarEmbedLeaf
import Spqr.PlanarEmbedF
import Spqr.PlanarEmbedV

/-!
# Assembling the per-item steps of `planarEmbed`

The leaf step `embedItem_step_leaf` over `GluedUpTo`; the per-type dispatch and the fold over the
reverse preorder live in `PlanarEmbedFacesSteps.lean` / `PlanarEmbedFacesFold.lean` (over
`GluedFaces`). Takes the structural well-formedness `WF` and the children shapes `ChildShape` of
the tree as hypotheses (obtained for `g.planarTree …` from `spqrTree_wf'` / `spqrTree_childShape`
via `planarRelabel_proj`).
-/

namespace Spqr

namespace PlanarSpqrTree

variable (t : PlanarSpqrTree)

/-- `O`/`I` items have nothing below them. -/
theorem subtreeEnd_leaf (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape) (i : Nat)
    (hi : i < t.size) (hty : t.types[i]! = .O ∨ t.types[i]! = .I) :
    t.subtreeEnd[i]! = i + 1 := by
  have h1 := hwf.preorder.subtree_eq i hi
  rw [hsh.leaf i hi (by rw [t.type_eq_of_lt i hi]; exact hty)] at h1
  simpa using h1

/-- `O`/`I` step: leaves without edges; `embedItem` is the identity. -/
theorem embedItem_step_leaf (g : Graph) (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (i : Nat) (hi : i < t.size) (hty : t.types[i]! = .O ∨ t.types[i]! = .I)
    (s : EmbedState) (h : t.GluedUpTo g (i + 1) s) :
    t.GluedUpTo g i ((t.embedItem i).run s).2 :=
  t.embedItem_step_leaf_of g i hi hty (t.subtreeEnd_leaf hwf hsh i hi hty) s h

end PlanarSpqrTree

end Spqr
