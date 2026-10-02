import Spqr.PlanarEmbedLeaf
import Spqr.PlanarEmbedF

/-!
# Assembling the per-item steps of `planarEmbed`

`embedItem_step` dispatches on the item type; `gluedUpTo_planarEmbed` folds it over the reverse
preorder (`forM_reverse_range_inv`) from `gluedUpTo_init`. Both take the structural
well-formedness `WF` and the children shapes `ChildShape` of the tree as hypotheses (obtained for
`g.planarTree …` from `spqrTree_wf'` / `spqrTree_childShape` via `planarRelabel_proj`).
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

theorem embedItem_step (g : Graph) (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hrep : t.toSpqrTree.Represents g) (hsep : t.toSpqrTree.PieceSep g)
    (i : Nat) (hi : i < t.size) (hall : t.nodePlanar.all id = true)
    (s : EmbedState) (h : t.GluedUpTo g (i + 1) s) :
    t.GluedUpTo g i ((t.embedItem i).run s).2 := by
  match hty : t.types[i]! with
  | .F => exact t.embedItem_step_F g hwf hsh hrep hsep i hi hty s h
  | .V => exact t.embedItem_step_V g hwf hsh hrep hsep i hi hty s h
  | .Q => exact t.embedItem_step_Q g hwf hsh hrep hsep i hi hty s h
  | .O => exact t.embedItem_step_leaf g hwf hsh i hi (Or.inl hty) s h
  | .I => exact t.embedItem_step_leaf g hwf hsh i hi (Or.inr hty) s h
  | .S => exact t.embedItem_step_node g hwf hsh hrep hsep i hi (Or.inl hty) hall s h
  | .P => exact t.embedItem_step_node g hwf hsh hrep hsep i hi (Or.inr (Or.inl hty)) hall s h
  | .R => exact t.embedItem_step_node g hwf hsh hrep hsep i hi (Or.inr (Or.inr hty)) hall s h

theorem gluedUpTo_planarEmbed (g : Graph) (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hrep : t.toSpqrTree.Represents g) (hsep : t.toSpqrTree.PieceSep g)
    (hall : t.nodePlanar.all id = true) :
    t.GluedUpTo g 0 (((List.range t.size).reverse.forM t.embedItem).run t.initState).2 :=
  t.forM_reverse_range_inv (fun i s => t.GluedUpTo g i s) t.size
    (fun i s hi h => t.embedItem_step g hwf hsh hrep hsep i hi hall s h) _ (t.gluedUpTo_init g)

end PlanarSpqrTree

end Spqr
