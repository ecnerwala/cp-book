import Spqr.PlanarEmbedFacesSteps
import Spqr.PlanarEmbedQ
import Spqr.PlanarNodeSpec
import Spqr.PlanarEmbedNodeExec
import Spqr.PlanarEmbedNodeFold

/-!
# The fold of `planarEmbed` over the same-witness invariant `GluedFaces`

`embedItem_step_faces` dispatches on the item type (leaf, `F`, `V`, `Q` proved; `S`/`P`/`R` from
`PlanarEmbedNodeFold.lean`, through `nodeFold_capped`) and `gluedFaces_planarEmbed` folds it over the reverse
preorder from `gluedFaces_init`.
-/

namespace Spqr

namespace PlanarSpqrTree

variable (t : PlanarSpqrTree)

theorem embedItem_step_faces (g : Graph) (hg : g.WF) (hwf : t.toSpqrTree.WF)
    (hsh : t.toSpqrTree.ChildShape) (hrep : t.toSpqrTree.Represents g)
    (hsep : t.toSpqrTree.PieceSep g) (i : Nat) (hi : i < t.size)
    (hloc : ∀ i, i < t.size → t.types[i]! = .S ∨ t.types[i]! = .P ∨ t.types[i]! = .R →
      IsPlanarEmbedding (t.localSkeleton i) (t.toSpqrTree.nVerts i) (t.nodeRot i))
    (hlay : ∀ i, i < t.size → t.LayoutAt i) (hclosed : t.NodeRotClosed) (hcor : t.NodeCorners)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    t.GluedFaces g i ((t.embedItem i).run s).2 := by
  match hty : t.types[i]! with
  | .F => exact t.embedItem_step_F_faces g hwf hsh hrep hsep i hi hty s h
  | .V => exact t.embedItem_step_V_faces g hwf hsh hrep hsep i hi hty s h
  | .Q => exact t.embedItem_step_Q_faces g hg hwf hsh hrep hsep i hi hty s h
  | .O => exact t.embedItem_step_leaf_faces g hwf hsh i hi (Or.inl hty) s h
  | .I => exact t.embedItem_step_leaf_faces g hwf hsh i hi (Or.inr hty) s h
  | .S => exact t.embedItem_step_node_faces g hwf hsh hrep hsep i hi (Or.inl hty) (hloc i hi (Or.inl hty)) (hlay i hi) hclosed hcor s h
  | .P =>
    exact t.embedItem_step_node_faces g hwf hsh hrep hsep i hi (Or.inr (Or.inl hty))
      (hloc i hi (Or.inr (Or.inl hty))) (hlay i hi) hclosed hcor s h
  | .R =>
    exact t.embedItem_step_node_faces g hwf hsh hrep hsep i hi (Or.inr (Or.inr hty))
      (hloc i hi (Or.inr (Or.inr hty))) (hlay i hi) hclosed hcor s h

theorem gluedFaces_planarEmbed (g : Graph) (hg : g.WF) (hwf : t.toSpqrTree.WF)
    (hsh : t.toSpqrTree.ChildShape) (hrep : t.toSpqrTree.Represents g)
    (hsep : t.toSpqrTree.PieceSep g)
    (hloc : ∀ i, i < t.size → t.types[i]! = .S ∨ t.types[i]! = .P ∨ t.types[i]! = .R →
      IsPlanarEmbedding (t.localSkeleton i) (t.toSpqrTree.nVerts i) (t.nodeRot i))
    (hlay : ∀ i, i < t.size → t.LayoutAt i) (hclosed : t.NodeRotClosed) (hcor : t.NodeCorners) :
    t.GluedFaces g 0 (((List.range t.size).reverse.forM t.embedItem).run t.initState).2 :=
  t.forM_reverse_range_inv (fun i s => t.GluedFaces g i s) t.size
    (fun i s hi h => t.embedItem_step_faces g hg hwf hsh hrep hsep i hi hloc hlay hclosed hcor s h) _
    (t.gluedFaces_init g)

end PlanarSpqrTree

end Spqr
