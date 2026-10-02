import Spqr.PlanarShape
import Spqr.RelabelRep
import Spqr.WalkWF

/-!
# `ChildShape` of the relabelled tree

The children of `idx i` in the output are the images of `items.ch i` (`RelabelOK.mem_children_iff`),
so each `Items.Shapes` clause transfers to the output tree.
-/

namespace Spqr

namespace RelabelOK

variable {g : Graph} {items : Items} {t : SpqrTree} {idx : ItemId → Nat}
  (h : RelabelOK g items t idx)
include h

theorem children_nil {i : ItemId} (hi : i < items.size) (hc : items.ch i = []) :
    t.children (idx i) = [] := by
  rw [List.eq_nil_iff_forall_not_mem]
  intro n hn
  obtain ⟨c, hc', -⟩ := (h.mem_children_iff hi n).1 hn
  simp [hc] at hc'

theorem childShape : t.ChildShape where
  leaf := by
    intro n hn hty
    obtain ⟨i, hi, rfl⟩ := h.idx_surj hn
    rw [h.type_eq hi] at hty
    exact h.children_nil hi (h.shapes.i_o_leaf i hi hty.symm)

end RelabelOK

theorem relabelTree_childShape (g : Graph) (items : Items) (h : items.WF g) :
    (relabelTree g items).ChildShape := by
  obtain ⟨idx, hok⟩ := relabelOK_of_wf g items h
  exact hok.childShape

theorem spqrTree_childShape (g : Graph) (tern : Bool) (vo eo : List Nat) :
    (g.spqrTree tern vo eo).ChildShape := by
  rw [spqrTree_eq]; exact relabelTree_childShape g _ (walk_items_wf g tern vo eo)

end Spqr
