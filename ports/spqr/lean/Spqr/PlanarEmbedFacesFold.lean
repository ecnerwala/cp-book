import Spqr.PlanarEmbedFacesSteps
import Spqr.PlanarEmbedQ

/-!
# The fold of `planarEmbed` over the same-witness invariant `GluedFaces`

`embedItem_step_faces` dispatches on the item type (leaf, `F`, `V`, `Q` proved; `S`/`P`/`R`
admitted as `embedItem_step_node_faces`) and `gluedFaces_planarEmbed` folds it over the reverse
preorder from `gluedFaces_init`.
-/

namespace Spqr

namespace PlanarSpqrTree

variable (t : PlanarSpqrTree)

/-- `S`/`P`/`R` step over `GluedFaces`: the node's local rotation (`nodePlanar_sound`) is 2-summed
with each child's piece through the twin virtual edge, the cap's quarter-edges become the exposed
ends, and the cap's two exposed pairs are cofacial in the resulting witness (they are the two
sides of the removed cap edge). Admitted: the executable node branch (`neRotAdj` traversal, twin
lookup, inner `V` items) has no unfolding lemmas yet, and `nodePlanar_sound` is stated for
`g.planarTree …` rather than an arbitrary `t` with `t.nodePlanar.all id = true`. -/
theorem embedItem_step_node_faces (g : Graph) (hwf : t.toSpqrTree.WF)
    (hsh : t.toSpqrTree.ChildShape) (hrep : t.toSpqrTree.Represents g)
    (hsep : t.toSpqrTree.PieceSep g) (i : Nat) (hi : i < t.size)
    (hty : t.types[i]! = .S ∨ t.types[i]! = .P ∨ t.types[i]! = .R)
    (hall : t.nodePlanar.all id = true)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    t.GluedFaces g i ((t.embedItem i).run s).2 := by
  sorry

theorem embedItem_step_faces (g : Graph) (hg : g.WF) (hwf : t.toSpqrTree.WF)
    (hsh : t.toSpqrTree.ChildShape) (hrep : t.toSpqrTree.Represents g)
    (hsep : t.toSpqrTree.PieceSep g) (i : Nat) (hi : i < t.size)
    (hall : t.nodePlanar.all id = true)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    t.GluedFaces g i ((t.embedItem i).run s).2 := by
  match hty : t.types[i]! with
  | .F => exact t.embedItem_step_F_faces g hwf hsh hrep hsep i hi hty s h
  | .V => exact t.embedItem_step_V_faces g hwf hsh hrep hsep i hi hty s h
  | .Q => exact t.embedItem_step_Q_faces g hg hwf hsh hrep hsep i hi hty s h
  | .O => exact t.embedItem_step_leaf_faces g hwf hsh i hi (Or.inl hty) s h
  | .I => exact t.embedItem_step_leaf_faces g hwf hsh i hi (Or.inr hty) s h
  | .S => exact t.embedItem_step_node_faces g hwf hsh hrep hsep i hi (Or.inl hty) hall s h
  | .P =>
    exact t.embedItem_step_node_faces g hwf hsh hrep hsep i hi (Or.inr (Or.inl hty)) hall s h
  | .R =>
    exact t.embedItem_step_node_faces g hwf hsh hrep hsep i hi (Or.inr (Or.inr hty)) hall s h

theorem gluedFaces_planarEmbed (g : Graph) (hg : g.WF) (hwf : t.toSpqrTree.WF)
    (hsh : t.toSpqrTree.ChildShape) (hrep : t.toSpqrTree.Represents g)
    (hsep : t.toSpqrTree.PieceSep g) (hall : t.nodePlanar.all id = true) :
    t.GluedFaces g 0 (((List.range t.size).reverse.forM t.embedItem).run t.initState).2 :=
  t.forM_reverse_range_inv (fun i s => t.GluedFaces g i s) t.size
    (fun i s hi h => t.embedItem_step_faces g hg hwf hsh hrep hsep i hi hall s h) _
    (t.gluedFaces_init g)

end PlanarSpqrTree

end Spqr
