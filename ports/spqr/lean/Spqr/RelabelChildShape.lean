import Spqr.PlanarShape
import Spqr.RelabelRep
import Spqr.WalkItemsWF
import Spqr.RelabelWF

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
  tile := by
    intro n hn k hk
    obtain ⟨i, hi, rfl⟩ := h.idx_surj hn
    have H : RelabelAll g items t idx := ⟨h.wf, h.ridx, h.node⟩
    obtain ⟨pos, hl⟩ := (h.node i hi).layout
    set L := items.ordered g i (t.nvRange (idx i)).1 pos with hL
    rw [hl.children] at hk ⊢
    simp only [List.length_map, ← hL] at hk ⊢
    rw [getElem!_pos _ k (by simpa using hk), List.getElem_map, ← getElem!_pos L k hk,
      RelabelAll.child_idx_eq hl hk, ← hL, ← H.children_sum hl k hk.le]
    have : (((L.map idx).take k).map fun c => t.subtreeEnd[c]! - c) =
        (List.range k).map fun j => t.subtreeEnd[idx L[j]!]! - idx L[j]! := by
      apply List.ext_getElem (by simp; omega)
      intro j h1 h2
      simp only [List.getElem_map, List.getElem_take, List.getElem_range]
      rw [getElem!_pos L j (by simp at h2; omega)]
    rw [this]; simp only [← hL]; omega

end RelabelOK

theorem relabelTree_childShape (g : Graph) (items : Items) (h : items.WF g) :
    (relabelTree g items).ChildShape := by
  obtain ⟨idx, hok⟩ := relabelOK_of_wf g items h
  exact hok.childShape

theorem spqrTree_childShape (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) : (g.spqrTree tern vo eo).ChildShape := by
  rw [spqrTree_eq]; exact relabelTree_childShape g _ (walk_items_wf g hg tern vo eo hvo heo)

end Spqr
