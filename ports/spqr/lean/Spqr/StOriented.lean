import Spqr.StWalk
import Spqr.WalkItemsWF
import Spqr.RelabelOwn

/-!
# `Items.ROriented` from the st-ordering

`Items.StNumbered` orients every R item's edge children along its vertex list, which is
`Items.ROriented` once `vertList` and `nvList` are identified (`Items.Tree`).
`walk_items_rOriented'` (from `walk_st`) supplies it to `Correctness.spqrTree_wf`.
-/

namespace Spqr

theorem Items.vertList_eq_nvList (g : Graph) {items : Items} (ht : items.Tree g) (i : ItemId) :
    items.vertList i = items.nvList g i := by
  unfold Items.vertList Items.nvList
  rw [ht.filter_V_eq]

theorem Items.rOriented_of_stNumbered (g : Graph) {items : Items} (hst : items.StNumbered)
    (h : items.WF g) : items.ROriented g := by
  intro i hi hR
  obtain ⟨s, t, -, -, hor⟩ := hst i hi (Or.inr (Or.inr hR))
  refine ⟨h.nv_nodup hi, ?_⟩
  rw [← Items.vertList_eq_nvList g h.tree]
  exact hor

theorem walk_items_rOriented' (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    Items.ROriented g (g.walk tern (g.dfsForest vo eo)).items :=
  Items.rOriented_of_stNumbered g (walk_st g hg tern vo eo hvo heo (walk_items_wf g hg tern vo eo hvo heo))
    (walk_items_wf g hg tern vo eo hvo heo)

end Spqr
