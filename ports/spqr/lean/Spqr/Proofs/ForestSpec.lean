import Mathlib.Tactic.Tauto
import Spqr.Proofs.Dfs
import Spqr.Proofs.Interval

/-!
# `DfsForestSpec` for the forest `dfsForest` computes

Discharges the phase-1 hypotheses of `Spqr/SepPair.lean` from `Proofs/Dfs.lean`
(`dfsForest_spanning'`, `dfsForest_wf`) and one extra invariant of `dfsVisit`
(`dfsVisit_outsIn`: every out-edge is an adjacency entry of its vertex).
-/

namespace Spqr

/-! ### `depthAt` -/

mutual
theorem DfsTree.depthAt_eq_none : ∀ (t : DfsTree) {x d₀ : Nat}, x ∉ t.verts → t.depthAt x d₀ = none
  | .node v outs, x, d₀, h => by
    have hv : v ≠ x := fun hvx => h (hvx ▸ List.mem_cons_self)
    simp only [DfsTree.depthAt, hv, ↓reduceIte]
    exact DfsOut.depthAtList_eq_none outs fun h' => h (List.mem_cons_of_mem _ h')
theorem DfsOut.depthAtList_eq_none :
    ∀ (outs : List DfsOut) {x d₀ : Nat}, x ∉ DfsOut.vertsList outs →
      DfsOut.depthAtList x d₀ outs = none
  | [], _, _, _ => rfl
  | .back _ _ _ :: rest, x, d₀, h => DfsOut.depthAtList_eq_none rest h
  | .tree e cls child :: rest, x, d₀, h => by
    simp only [DfsOut.vertsList, List.mem_append, not_or] at h
    show (child.depthAt x d₀).orElse _ = none
    rw [DfsTree.depthAt_eq_none child h.1]
    exact DfsOut.depthAtList_eq_none rest h.2
end

theorem DfsOut.depthAtList_eq_of_mem {outs : List DfsOut} (hnd : (DfsOut.vertsList outs).Nodup)
    {e : Nat} {cls : OutClass} {child : DfsTree} (hmem : DfsOut.tree e cls child ∈ outs)
    {x d₀ : Nat} (hx : x ∈ child.verts) : DfsOut.depthAtList x d₀ outs = child.depthAt x d₀ := by
  induction outs with
  | nil => exact absurd hmem List.not_mem_nil
  | cons o rest ih =>
    rcases List.mem_cons.mp hmem with rfl | hmem'
    · simp only [DfsOut.vertsList, List.nodup_append] at hnd
      show (child.depthAt x d₀).orElse _ = _
      rw [DfsOut.depthAtList_eq_none rest fun h => hnd.2.2 _ hx _ h rfl]
      cases child.depthAt x d₀ <;> rfl
    · cases o with
      | back => exact ih hnd hmem'
      | tree e' cls' child' =>
        simp only [DfsOut.vertsList, List.nodup_append] at hnd
        have hx' : x ∈ DfsOut.vertsList rest := DfsOut.mem_vertsList.mpr ⟨_, _, _, hmem', hx⟩
        show (child'.depthAt x d₀).orElse _ = _
        rw [DfsTree.depthAt_eq_none child' fun h => hnd.2.2 _ h _ hx' rfl]
        exact ih hnd.2.1 hmem'

theorem DfsTree.depthAt_v (t : DfsTree) (d₀ : Nat) : t.depthAt t.v d₀ = some d₀ := by
  cases t; simp [DfsTree.depthAt, DfsTree.v]

theorem DfsTree.Sub.depthAt_eq {s t : DfsTree} (hnd : t.verts.Nodup) (h : s.Sub t) {x : Nat}
    (hx : x ∈ s.verts) (d₀ : Nat) :
    ∃ k, t.depthAt s.v d₀ = some k ∧ t.depthAt x d₀ = s.depthAt x k := by
  induction h generalizing d₀ with
  | refl => exact ⟨d₀, DfsTree.depthAt_v _ _, rfl⟩
  | @step v outs e cls child hsub hmem ih =>
    simp only [DfsTree.verts, List.nodup_cons] at hnd
    have hsv : s.v ∈ child.verts := hsub.verts_subset s.v_mem_verts
    have hxc : x ∈ child.verts := hsub.verts_subset hx
    have hne1 : v ≠ s.v := fun h => hnd.1 (h ▸ DfsOut.mem_vertsList.mpr ⟨_, _, _, hmem, hsv⟩)
    have hne2 : v ≠ x := fun h => hnd.1 (h ▸ DfsOut.mem_vertsList.mpr ⟨_, _, _, hmem, hxc⟩)
    simp only [DfsTree.depthAt, hne1, hne2, ↓reduceIte]
    rw [DfsOut.depthAtList_eq_of_mem hnd.2 hmem hsv, DfsOut.depthAtList_eq_of_mem hnd.2 hmem hxc]
    exact ih (hnd.2.sublist (DfsOut.verts_infix_of_mem hmem).sublist) _

theorem DfsTree.Sub.verts_nodup {s t : DfsTree} (hnd : t.verts.Nodup) (h : s.Sub t) :
    s.verts.Nodup := by
  induction h with
  | refl => exact hnd
  | @step v outs e cls child _ hmem ih =>
    simp only [DfsTree.verts, List.nodup_cons] at hnd
    exact ih (hnd.2.sublist (DfsOut.verts_infix_of_mem hmem).sublist)

/-- The depth of a tree child is one more than its parent's. -/
theorem DfsTree.Sub.depthAt_child {s t : DfsTree} (hnd : t.verts.Nodup) (h : s.Sub t)
    {e : Nat} {cls : OutClass} {c : DfsTree} (hc : DfsOut.tree e cls c ∈ s.outs) :
    ∃ k, t.depthAt s.v 0 = some k ∧ t.depthAt c.v 0 = some (k + 1) := by
  have hnd' := h.verts_nodup hnd
  obtain ⟨v, outs⟩ := s
  have hcv : c.v ∈ (DfsTree.node v outs).verts :=
    List.mem_cons_of_mem _ (DfsOut.mem_vertsList.mpr ⟨_, _, _, hc, c.v_mem_verts⟩)
  obtain ⟨k, hk, hx⟩ := h.depthAt_eq hnd hcv 0
  refine ⟨k, hk, hx.trans ?_⟩
  simp only [DfsTree.verts, List.nodup_cons] at hnd'
  have hne : v ≠ c.v := fun h => hnd'.1 (h ▸ DfsOut.mem_vertsList.mpr ⟨_, _, _, hc, c.v_mem_verts⟩)
  simp only [DfsTree.depthAt, hne, ↓reduceIte]
  rw [DfsOut.depthAtList_eq_of_mem hnd'.2 hc c.v_mem_verts, DfsTree.depthAt_v]

/-- `DfsData.ofForest`'s depth is `depthAt` in the tree containing the vertex. -/
theorem DfsData.ofForest_depth_eq {forest : List DfsTree} (hnd : (forest.flatMap DfsTree.verts).Nodup)
    {t : DfsTree} (ht : t ∈ forest) {x k : Nat} (hx : t.depthAt x 0 = some k) :
    (DfsData.ofForest forest).depth x = k := by
  show (forest.findSome? (·.depthAt x 0)).getD 0 = k
  induction forest with
  | nil => exact absurd ht List.not_mem_nil
  | cons t₀ rest ih =>
    rw [List.flatMap_cons, List.nodup_append] at hnd
    rcases List.mem_cons.mp ht with rfl | ht'
    · simp [hx]
    · have hxt : x ∈ t.verts := by
        by_contra hc
        rw [t.depthAt_eq_none hc] at hx
        cases hx
      rw [List.findSome?_cons, t₀.depthAt_eq_none fun h =>
        hnd.2.2 _ h _ ((List.infix_flatMap_of_mem ht').subset hxt) rfl]
      exact ih hnd.2.1 ht'

theorem DfsData.ofForest_depth_parent {forest : List DfsTree}
    (hnd : (forest.flatMap DfsTree.verts).Nodup) {t s : DfsTree} (ht : t ∈ forest) (hs : s.Sub t)
    {e : Nat} {cls : OutClass} {c : DfsTree} (hc : DfsOut.tree e cls c ∈ s.outs) :
    (DfsData.ofForest forest).depth c.v = (DfsData.ofForest forest).depth s.v + 1 := by
  have hnd' : t.verts.Nodup := hnd.sublist (List.infix_flatMap_of_mem ht).sublist
  obtain ⟨k, hk, hk'⟩ := hs.depthAt_child hnd' hc
  rw [DfsData.ofForest_depth_eq hnd ht hk, DfsData.ofForest_depth_eq hnd ht hk']

theorem DfsData.ofForest_depth_root {forest : List DfsTree}
    (hnd : (forest.flatMap DfsTree.verts).Nodup) {t : DfsTree} (ht : t ∈ forest) :
    (DfsData.ofForest forest).depth t.v = 0 :=
  DfsData.ofForest_depth_eq hnd ht (t.depthAt_v 0)

/-! ### All out-edges with their vertices -/

mutual
/-- The pairs `(vertex, out-edge)` of a tree, in the order of `DfsTree.edges`. -/
def DfsTree.allOuts : DfsTree → List (Nat × DfsOut)
  | .node v outs => DfsOut.allOutsList v outs
def DfsOut.allOutsList (v : Nat) : List DfsOut → List (Nat × DfsOut)
  | [] => []
  | .back e dest cls :: rest => (v, .back e dest cls) :: DfsOut.allOutsList v rest
  | .tree e cls child :: rest =>
    (v, .tree e cls child) :: (child.allOuts ++ DfsOut.allOutsList v rest)
end

mutual
theorem DfsTree.allOuts_map_e : ∀ t : DfsTree, t.allOuts.map (·.2.e) = t.edges
  | .node v outs => DfsOut.allOutsList_map_e v outs
theorem DfsOut.allOutsList_map_e :
    ∀ (v : Nat) (outs : List DfsOut), (DfsOut.allOutsList v outs).map (·.2.e) = DfsOut.edgesList outs
  | _, [] => rfl
  | v, .back e dest cls :: rest => by
    simp only [DfsOut.allOutsList, List.map_cons, DfsOut.edgesList]
    rw [DfsOut.allOutsList_map_e v rest]; rfl
  | v, .tree e cls child :: rest => by
    simp only [DfsOut.allOutsList, List.map_cons, List.map_append, DfsOut.edgesList]
    rw [DfsOut.allOutsList_map_e v rest, DfsTree.allOuts_map_e child]; rfl
end

mutual
theorem DfsTree.verts_eq_allOuts : ∀ t : DfsTree,
    t.verts = t.v :: (t.allOuts.filter (·.2.isTree)).map (·.2.dest)
  | .node v outs => by
    simp only [DfsTree.verts, DfsTree.v, DfsTree.allOuts]
    rw [DfsOut.vertsList_eq_allOutsList v outs]
theorem DfsOut.vertsList_eq_allOutsList : ∀ (v : Nat) (outs : List DfsOut),
    DfsOut.vertsList outs = ((DfsOut.allOutsList v outs).filter (·.2.isTree)).map (·.2.dest)
  | _, [] => rfl
  | v, .back e dest cls :: rest => by
    simp only [DfsOut.vertsList, DfsOut.allOutsList, List.filter_cons, DfsOut.isTree,
      Bool.false_eq_true, ↓reduceIte]
    exact DfsOut.vertsList_eq_allOutsList v rest
  | v, .tree e cls child :: rest => by
    simp only [DfsOut.vertsList, DfsOut.allOutsList, List.filter_cons, List.filter_append,
      DfsOut.isTree, ↓reduceIte, List.map_cons, List.map_append]
    rw [DfsOut.vertsList_eq_allOutsList v rest, DfsTree.verts_eq_allOuts child]
    rfl
end

mutual
theorem DfsTree.mem_allOuts : ∀ (t : DfsTree) {x : Nat} {o : DfsOut},
    (x, o) ∈ t.allOuts ↔ o ∈ t.outsAt x
  | .node v outs, x, o => by
    simp only [DfsTree.allOuts, DfsTree.outsAt, List.mem_append]
    rw [DfsOut.mem_allOutsList v outs]
    by_cases hv : v = x
    · subst hv; simp
    · simp [hv, Ne.symm hv]
theorem DfsOut.mem_allOutsList : ∀ (v : Nat) (outs : List DfsOut) {x : Nat} {o : DfsOut},
    (x, o) ∈ DfsOut.allOutsList v outs ↔ (x = v ∧ o ∈ outs) ∨ o ∈ DfsOut.outsAtList x outs
  | _, [], _, _ => by simp [DfsOut.allOutsList, DfsOut.outsAtList]
  | v, .back e dest cls :: rest, x, o => by
    simp only [DfsOut.allOutsList, List.mem_cons, Prod.mk.injEq, DfsOut.outsAtList]
    rw [DfsOut.mem_allOutsList v rest]
    tauto
  | v, .tree e cls child :: rest, x, o => by
    simp only [DfsOut.allOutsList, List.mem_cons, Prod.mk.injEq, List.mem_append,
      DfsOut.outsAtList]
    rw [DfsOut.mem_allOutsList v rest, DfsTree.mem_allOuts child]
    tauto
end

theorem List.flatMap_sublist_flatMap {α β : Type} {f g : α → List β} :
    ∀ {l : List α}, (∀ a ∈ l, (f a).Sublist (g a)) → (l.flatMap f).Sublist (l.flatMap g)
  | [], _ => List.Sublist.refl _
  | a :: l, h => by
    simp only [List.flatMap_cons]
    exact (h a List.mem_cons_self).append
      (List.flatMap_sublist_flatMap fun a' ha' => h a' (List.mem_cons_of_mem _ ha'))

theorem List.flatMap_perm_flatMap {α β : Type} {f g : α → List β} :
    ∀ {l : List α}, (∀ a ∈ l, (f a).Perm (g a)) → (l.flatMap f).Perm (l.flatMap g)
  | [], _ => List.Perm.refl _
  | a :: l, h => by
    simp only [List.flatMap_cons]
    exact (h a List.mem_cons_self).append
      (List.flatMap_perm_flatMap fun a' ha' => h a' (List.mem_cons_of_mem _ ha'))

theorem DfsOut.map_sublist_allOutsList (v : Nat) :
    ∀ outs : List DfsOut, (outs.map (v, ·)).Sublist (DfsOut.allOutsList v outs)
  | [] => List.Sublist.refl _
  | .back _ _ _ :: rest => (DfsOut.map_sublist_allOutsList v rest).cons_cons _
  | .tree _ _ _ :: rest =>
    ((DfsOut.map_sublist_allOutsList v rest).trans (List.sublist_append_right _ _)).cons_cons _

theorem DfsOut.allOuts_infix_of_mem {outs : List DfsOut} {v e : Nat} {cls : OutClass}
    {child : DfsTree} (h : DfsOut.tree e cls child ∈ outs) :
    child.allOuts <:+: DfsOut.allOutsList v outs := by
  induction outs with
  | nil => exact absurd h List.not_mem_nil
  | cons o rest ih =>
    rcases List.mem_cons.mp h with rfl | h'
    · exact ((List.prefix_append _ _).isInfix).trans (List.infix_cons (List.infix_refl _))
    · cases o with
      | back => exact (ih h').trans (List.infix_cons (List.infix_refl _))
      | tree =>
        exact (ih h').trans (((List.suffix_append _ _).isInfix).trans
          (List.infix_cons (List.infix_refl _)))

theorem DfsTree.Sub.allOuts_infix {s t : DfsTree} (h : s.Sub t) : s.allOuts <:+: t.allOuts := by
  induction h with
  | refl => exact List.infix_refl _
  | step _ hmem ih => exact ih.trans (DfsOut.allOuts_infix_of_mem hmem)

mutual
theorem DfsTree.edgePostorder_perm_edges : ∀ t : DfsTree, t.edgePostorder.Perm t.edges
  | .node _ outs => DfsOut.edgePostorderList_perm_edgesList outs
theorem DfsOut.edgePostorderList_perm_edgesList :
    ∀ outs : List DfsOut, (DfsOut.edgePostorderList outs).Perm (DfsOut.edgesList outs)
  | [] => List.Perm.refl _
  | .back e _ _ :: rest => (DfsOut.edgePostorderList_perm_edgesList rest).cons e
  | .tree e _ child :: rest => by
    simp only [DfsOut.edgePostorderList, DfsOut.edgesList]
    exact List.perm_middle.trans (((DfsTree.edgePostorder_perm_edges child).append
      (DfsOut.edgePostorderList_perm_edgesList rest)).cons e)
end

/-! ### The forest -/

namespace DfsData

variable {forest : List DfsTree}

theorem mem_forestAllOuts {x : Nat} {o : DfsOut} :
    (x, o) ∈ forest.flatMap DfsTree.allOuts ↔ o ∈ (ofForest forest).outs x := by
  show _ ↔ o ∈ forest.flatMap (·.outsAt x)
  simp only [List.mem_flatMap, DfsTree.mem_allOuts]

theorem forestAllOuts_map_e :
    (forest.flatMap DfsTree.allOuts).map (·.2.e) = forest.flatMap DfsTree.edges := by
  rw [List.map_flatMap]
  exact List.flatMap_congr fun t _ => t.allOuts_map_e

theorem forestAllOuts_dests_sublist :
    (((forest.flatMap DfsTree.allOuts).filter (·.2.isTree)).map (·.2.dest)).Sublist
      (forest.flatMap DfsTree.verts) := by
  rw [List.filter_flatMap, List.map_flatMap]
  exact List.flatMap_sublist_flatMap fun t _ => by
    rw [t.verts_eq_allOuts]; exact List.sublist_cons_self _ _

theorem edgePostorderForest_perm :
    (edgePostorderForest forest).Perm (forest.flatMap DfsTree.edges) :=
  List.flatMap_perm_flatMap fun t _ => t.edgePostorder_perm_edges

theorem exists_sub_of_mem_outs {v : Nat} {o : DfsOut} (ho : o ∈ (ofForest forest).outs v)
    (hnd : (forest.flatMap DfsTree.verts).Nodup) :
    ∃ t ∈ forest, ∃ s : DfsTree, s.Sub t ∧ s.v = v ∧ (ofForest forest).outs v = s.outs := by
  obtain ⟨t, ht, ho'⟩ := List.mem_flatMap.mp (show o ∈ forest.flatMap (·.outsAt v) from ho)
  have hv : v ∈ t.verts := by
    by_contra h
    rw [t.outsAt_eq_nil h] at ho'
    exact List.not_mem_nil ho'
  obtain ⟨s, hs, rfl⟩ := t.exists_sub_of_mem_verts hv
  exact ⟨t, ht, s, hs, rfl, ofForest_outs hnd ht hs⟩

section Structural

variable (hnd : (forest.flatMap DfsTree.verts).Nodup)
  (hend : (forest.flatMap DfsTree.edges).Nodup)

include hend in
theorem ofForest_out_inj (v v' : Nat) {o o' : DfsOut} (ho : o ∈ (ofForest forest).outs v)
    (ho' : o' ∈ (ofForest forest).outs v') (he : o.e = o'.e) : v = v' ∧ o = o' := by
  have hn : ((forest.flatMap DfsTree.allOuts).map (·.2.e)).Nodup := by
    rw [forestAllOuts_map_e]; exact hend
  have h := List.inj_on_of_nodup_map hn (mem_forestAllOuts.mpr ho) (mem_forestAllOuts.mpr ho') he
  exact Prod.mk.inj h

include hnd hend in
theorem ofForest_outs_nodup (v : Nat) : ((ofForest forest).outs v).Nodup := by
  by_cases h : (ofForest forest).outs v = []
  · rw [h]; exact List.nodup_nil
  · obtain ⟨o, ho⟩ := List.exists_mem_of_ne_nil _ h
    obtain ⟨t, ht, s, hs, rfl, heq⟩ := exists_sub_of_mem_outs ho hnd
    rw [heq]
    obtain ⟨sv, outs⟩ := s
    have hsub : ((outs.map (sv, ·)).map (·.2.e)).Sublist (forest.flatMap DfsTree.edges) := by
      rw [← forestAllOuts_map_e]
      refine List.Sublist.map _ ((DfsOut.map_sublist_allOutsList sv outs).trans ?_)
      exact (hs.allOuts_infix.trans (List.infix_flatMap_of_mem ht)).sublist
    have := (hend.sublist hsub).of_map
    exact this.of_map _

include hnd in
theorem ofForest_tree_dest_inj {v v' : Nat} {o o' : DfsOut} (ho : o ∈ (ofForest forest).outs v)
    (ho' : o' ∈ (ofForest forest).outs v') (ht : o.isTree = true) (ht' : o'.isTree = true)
    (hd : o.dest = o'.dest) : v = v' ∧ o = o' := by
  have hn := hnd.sublist forestAllOuts_dests_sublist
  have h := List.inj_on_of_nodup_map hn (List.mem_filter.mpr ⟨mem_forestAllOuts.mpr ho, ht⟩)
    (List.mem_filter.mpr ⟨mem_forestAllOuts.mpr ho', ht'⟩) hd
  exact Prod.mk.inj h

include hnd in
theorem ofForest_depth_parent' {p c : Nat} (h : (ofForest forest).IsParent p c) :
    (ofForest forest).depth c = (ofForest forest).depth p + 1 := by
  obtain ⟨o, ho, hti, rfl⟩ := h
  obtain ⟨t, ht, s, hs, rfl, heq⟩ := exists_sub_of_mem_outs ho hnd
  rw [heq] at ho
  cases o with
  | back => cases hti
  | tree e cls child => exact ofForest_depth_parent hnd ht hs ho

include hnd in
theorem ofForest_depth_root' : (ofForest forest).depth (ofForest forest).root = 0 := by
  cases forest with
  | nil => rfl
  | cons t rest => exact ofForest_depth_root hnd List.mem_cons_self

end Structural

end DfsData

/-! ### `classify` facts -/

theorem classify_false_ne_component (d : Nat) (n : Lowvals) : classify d false n ≠ .component := by
  unfold classify; split <;> (try split) <;> simp

theorem classify_true_ne_selfLoop (d : Nat) (n : Lowvals) : classify d true n ≠ .selfLoop := by
  unfold classify; split <;> (try split) <;> simp

theorem classify_true_ne_backEdge (d l : Nat) (n : Lowvals) :
    classify d true n ≠ .ret l .backEdge := by
  unfold classify; split <;> (try split) <;> (try split) <;> simp_all

/-! ### `retDepths` -/

mutual
theorem DfsTree.mem_retDepths_iff (D : DfsTree → Nat) : ∀ (s : DfsTree),
    (∀ s' e cls c, s'.Sub s → DfsOut.tree e cls c ∈ s'.outs → D c = D s' + 1) → ∀ x,
    (x ∈ s.retDepths (D s) ↔
      ∃ s' : DfsTree, s'.Sub s ∧ ∃ o : DfsOut, o ∈ s'.outs ∧ o.isTree = false ∧ o.cls.lowval (D s') = x)
  | .node v outs, hD, x => by
    simp only [DfsTree.retDepths]
    rw [DfsOut.mem_retDepthsList_iff D outs (D (.node v outs))
      (fun e cls c h => hD _ e cls c (.refl _) h)
      (fun e cls c h s' e' cls' c' hs' h' => hD s' e' cls' c' (.step hs' h) h') x]
    constructor
    · rintro (⟨o, ho, hb, hx⟩ | ⟨e, cls, c, hmem, s', hs', o, ho, hb, hx⟩)
      · exact ⟨_, .refl _, o, ho, hb, hx⟩
      · exact ⟨s', .step hs' hmem, o, ho, hb, hx⟩
    · rintro ⟨s', hs', o, ho, hb, hx⟩
      cases hs' with
      | refl => exact .inl ⟨o, ho, hb, hx⟩
      | step hs'' hmem => exact .inr ⟨_, _, _, hmem, s', hs'', o, ho, hb, hx⟩
theorem DfsOut.mem_retDepthsList_iff (D : DfsTree → Nat) :
    ∀ (outs : List DfsOut) (dv : Nat),
    (∀ e cls c, DfsOut.tree e cls c ∈ outs → D c = dv + 1) →
    (∀ e cls c, DfsOut.tree e cls c ∈ outs → ∀ s' e' cls' c', s'.Sub c →
      DfsOut.tree e' cls' c' ∈ s'.outs → D c' = D s' + 1) → ∀ x,
    (x ∈ DfsOut.retDepthsList dv outs ↔
      (∃ o : DfsOut, o ∈ outs ∧ o.isTree = false ∧ o.cls.lowval dv = x) ∨
      ∃ e cls c, DfsOut.tree e cls c ∈ outs ∧
        ∃ s' : DfsTree, s'.Sub c ∧ ∃ o : DfsOut, o ∈ s'.outs ∧ o.isTree = false ∧ o.cls.lowval (D s') = x)
  | [], dv, _, _, x => by simp [DfsOut.retDepthsList]
  | .back e dest cls :: rest, dv, h1, h2, x => by
    simp only [DfsOut.retDepthsList, List.mem_cons]
    rw [DfsOut.mem_retDepthsList_iff D rest dv
      (fun e' cls' c' h => h1 e' cls' c' (List.mem_cons_of_mem _ h))
      (fun e' cls' c' h => h2 e' cls' c' (List.mem_cons_of_mem _ h)) x]
    constructor
    · rintro (rfl | ⟨o, ho, hb, hx⟩ | ⟨e', cls', c', hmem, rest'⟩)
      · exact .inl ⟨_, .inl rfl, rfl, rfl⟩
      · exact .inl ⟨o, .inr ho, hb, hx⟩
      · exact .inr ⟨e', cls', c', .inr hmem, rest'⟩
    · rintro (⟨o, rfl | ho, hb, hx⟩ | ⟨e', cls', c', h | hmem, rest'⟩)
      · exact .inl hx.symm
      · exact .inr (.inl ⟨o, ho, hb, hx⟩)
      · cases h
      · exact .inr (.inr ⟨e', cls', c', hmem, rest'⟩)
  | .tree e cls c :: rest, dv, h1, h2, x => by
    simp only [DfsOut.retDepthsList, List.mem_append, List.mem_cons]
    have hDc : D c = dv + 1 := h1 e cls c List.mem_cons_self
    rw [← hDc, DfsTree.mem_retDepths_iff D c (h2 e cls c List.mem_cons_self) x,
      DfsOut.mem_retDepthsList_iff D rest dv
        (fun e' cls' c' h => h1 e' cls' c' (List.mem_cons_of_mem _ h))
        (fun e' cls' c' h => h2 e' cls' c' (List.mem_cons_of_mem _ h)) x]
    constructor
    · rintro (⟨s', hs', o, ho, hb, hx⟩ | ⟨o, ho, hb, hx⟩ | ⟨e', cls', c', hmem, rest'⟩)
      · exact .inr ⟨e, cls, c, .inl rfl, s', hs', o, ho, hb, hx⟩
      · exact .inl ⟨o, .inr ho, hb, hx⟩
      · exact .inr ⟨e', cls', c', .inr hmem, rest'⟩
    · rintro (⟨o, rfl | ho, hb, hx⟩ | ⟨e', cls', c', h | hmem, rest'⟩)
      · cases hb
      · exact .inr (.inl ⟨o, ho, hb, hx⟩)
      · cases h; exact .inl rest'
      · exact .inr (.inr ⟨e', cls', c', hmem, rest'⟩)
end

/-! ### From `DfsTree.WF` to the `Spec` fields -/

namespace DfsData

variable {forest : List DfsTree}

/-- `anc` lists the proper ancestors of `x`, `anc[i]` at depth `i`. -/
def AncList (d : DfsData) (anc : List Nat) (x : Nat) : Prop :=
  anc.length = d.depth x ∧ ∀ i w, anc[i]? = some w → d.depth w = i ∧ d.Anc w x

theorem AncList.index {d : DfsData} {anc : List Nat} {v : Nat} (h : AncList d anc v) {i w : Nat}
    (hi : (anc ++ [v])[i]? = some w) :
    d.depth w = i ∧ d.Anc w v ∧ i ≤ d.depth v ∧ (i = d.depth v → w = v) := by
  have hlen := h.1
  by_cases hlt : i < anc.length
  · rw [List.getElem?_append_left hlt] at hi
    obtain ⟨h1, h2⟩ := h.2 i w hi
    exact ⟨h1, h2, by omega, fun h' => absurd h' (by omega)⟩
  · rw [List.getElem?_append_right (Nat.le_of_not_lt hlt)] at hi
    cases hk : i - anc.length with
    | zero =>
      rw [hk] at hi
      simp only [List.getElem?_cons_zero, Option.some.injEq] at hi
      subst hi
      exact ⟨by omega, Anc.refl _, by omega, fun _ => rfl⟩
    | succ k =>
      rw [hk] at hi
      simp at hi

theorem wf_sub_aux (hnd : (forest.flatMap DfsTree.verts).Nodup) {t₀ : DfsTree} (ht₀ : t₀ ∈ forest)
    {s t : DfsTree} (hs : s.Sub t) :
    ∀ anc, t.Sub t₀ → t.WF anc → AncList (ofForest forest) anc t.v →
      ∃ anc', s.WF anc' ∧ AncList (ofForest forest) anc' s.v := by
  induction hs with
  | refl => exact fun anc _ h1 h2 => ⟨anc, h1, h2⟩
  | @step v outs e cls child hsub hmem ih =>
    intro anc htt₀ hwf hanc
    simp only [DfsTree.WF] at hwf
    have h := hwf.2 _ hmem
    simp only [DfsOut.WF] at h
    have hout : (ofForest forest).outs v = outs := ofForest_outs hnd ht₀ htt₀
    have hpar : (ofForest forest).IsParent v child.v :=
      ⟨.tree e cls child, by rw [hout]; exact hmem, rfl, rfl⟩
    have hdc := ofForest_depth_parent' hnd hpar
    have hlen := hanc.1
    simp only [DfsTree.v] at hlen
    refine ih (anc ++ [v]) ((DfsTree.Sub.child hmem).trans htt₀) h.1 ⟨?_, ?_⟩
    · rw [List.length_append, List.length_singleton, hdc]; omega
    · intro i w hi
      by_cases hlt : i < anc.length
      · rw [List.getElem?_append_left hlt] at hi
        obtain ⟨h1, h2⟩ := hanc.2 i w hi
        exact ⟨h1, h2.tail hpar⟩
      · rw [List.getElem?_append_right (Nat.le_of_not_lt hlt)] at hi
        cases hk : i - anc.length with
        | zero =>
          rw [hk] at hi
          simp only [List.getElem?_cons_zero, Option.some.injEq] at hi
          subst hi
          exact ⟨by omega, .single hpar⟩
        | succ k =>
          rw [hk] at hi
          simp at hi

theorem wf_sub (hnd : (forest.flatMap DfsTree.verts).Nodup) (hwf : ∀ t ∈ forest, t.WF [])
    {t s : DfsTree} (ht : t ∈ forest) (hs : s.Sub t) :
    ∃ anc, s.WF anc ∧ AncList (ofForest forest) anc s.v :=
  wf_sub_aux hnd ht hs [] (.refl _) (hwf t ht)
    ⟨by simp [ofForest_depth_root hnd ht], fun i w h => by simp at h⟩

section WF

variable (hnd : (forest.flatMap DfsTree.verts).Nodup) (hwf : ∀ t ∈ forest, t.WF [])
include hnd hwf

theorem back_out_facts {u : Nat} {e dest : Nat} {cls : OutClass}
    (ho : DfsOut.back e dest cls ∈ (ofForest forest).outs u) :
    (ofForest forest).depth dest ≤ (ofForest forest).depth u ∧
      (ofForest forest).Anc dest u ∧
      cls = classify ((ofForest forest).depth u) false
        ((ofForest forest).depth dest, (ofForest forest).depth u) ∧
      ((ofForest forest).depth dest = (ofForest forest).depth u → dest = u) := by
  obtain ⟨t, ht, s, hs, rfl, heq⟩ := exists_sub_of_mem_outs ho hnd
  obtain ⟨anc, hswf, hanc⟩ := wf_sub hnd hwf ht hs
  obtain ⟨v, outs⟩ := s
  rw [heq] at ho
  simp only [DfsTree.WF] at hswf
  have h := hswf.2 _ ho
  simp only [DfsOut.WF] at h
  obtain ⟨i, hi, hcls⟩ := h
  obtain ⟨h1, h2, h3, h4⟩ := hanc.index hi
  have hlen := hanc.1
  simp only [DfsTree.v] at h1 h2 h3 h4 hlen ⊢
  rw [← hlen] at h3 h4 ⊢
  rw [h1]
  exact ⟨h3, h2, hcls, h4⟩

theorem tree_out_facts {u : Nat} {e : Nat} {cls : OutClass} {c : DfsTree}
    (ho : DfsOut.tree e cls c ∈ (ofForest forest).outs u) :
    cls = classify ((ofForest forest).depth u) true
        (low2 ((ofForest forest).depth u + 1) (c.retDepths ((ofForest forest).depth u + 1))) ∧
      ∀ y, y ∈ c.retDepths ((ofForest forest).depth u + 1) ↔ (ofForest forest).Returns c.v y := by
  obtain ⟨t, ht, s, hs, rfl, heq⟩ := exists_sub_of_mem_outs ho hnd
  obtain ⟨anc, hswf, hanc⟩ := wf_sub hnd hwf ht hs
  have hc : c.Sub t := (DfsTree.Sub.child (by rw [← heq]; exact ho)).trans hs
  have hdc : (ofForest forest).depth c.v = (ofForest forest).depth s.v + 1 :=
    ofForest_depth_parent' hnd ⟨_, ho, rfl, rfl⟩
  obtain ⟨v, outs⟩ := s
  rw [heq] at ho
  simp only [DfsTree.WF] at hswf
  have h := hswf.2 _ ho
  simp only [DfsOut.WF] at h
  have hlen : anc.length = (ofForest forest).depth v := hanc.1
  have hdc' : (ofForest forest).depth c.v = (ofForest forest).depth v + 1 := hdc
  show cls = classify ((ofForest forest).depth v) true _ ∧
    ∀ y, y ∈ c.retDepths ((ofForest forest).depth v + 1) ↔ _
  rw [hlen] at h
  refine ⟨h.2, fun y => ?_⟩
  rw [← hdc']
  rw [DfsTree.mem_retDepths_iff (fun s' => (ofForest forest).depth s'.v) c
    (fun s' e' cls' c' hs' hc' => ofForest_depth_parent hnd ht (hs'.trans hc) hc')]
  constructor
  · rintro ⟨s', hs', o, ho', hb, hx⟩
    refine ⟨s'.v, o, ?_, ofForest_anc_of_sub hnd ht hs' hc, hb, ?_⟩
    · rw [ofForest_outs hnd ht (hs'.trans hc)]; exact ho'
    · cases o with
      | tree => cases hb
      | back e' dest' cls' =>
        have ho'' : DfsOut.back e' dest' cls' ∈ (ofForest forest).outs s'.v := by
          rw [ofForest_outs hnd ht (hs'.trans hc)]; exact ho'
        obtain ⟨hle, -, hcls', -⟩ := back_out_facts hnd hwf ho''
        rw [← hx]
        show (ofForest forest).depth dest' = cls'.lowval ((ofForest forest).depth s'.v)
        rw [hcls', lowval_classify_back hle]
  · rintro ⟨u, o, ho', hanc', hb, hx⟩
    obtain ⟨s', hs', rfl⟩ := c.exists_sub_of_mem_verts ((ofForest_anc_iff hnd ht hc).mp hanc')
    have ho'' := ho'
    rw [ofForest_outs hnd ht (hs'.trans hc)] at ho''
    refine ⟨s', hs', o, ho'', hb, ?_⟩
    cases o with
    | tree => cases hb
    | back e' dest' cls' =>
      obtain ⟨hle, -, hcls', -⟩ := back_out_facts hnd hwf ho'
      rw [← hx]
      show cls'.lowval ((ofForest forest).depth s'.v) = (ofForest forest).depth dest'
      rw [hcls', lowval_classify_back hle]

theorem ofForest_back_anc {v : Nat} {o : DfsOut} (ho : o ∈ (ofForest forest).outs v)
    (hb : o.isTree = false) : (ofForest forest).Anc o.dest v := by
  cases o with
  | tree => cases hb
  | back e dest cls => exact (back_out_facts hnd hwf ho).2.1

theorem ofForest_sorted (v : Nat) :
    ((ofForest forest).outs v).Pairwise fun o o' => o.cls.rank ≤ o'.cls.rank := by
  by_cases h : (ofForest forest).outs v = []
  · rw [h]; exact List.Pairwise.nil
  · obtain ⟨o, ho⟩ := List.exists_mem_of_ne_nil _ h
    obtain ⟨t, ht, s, hs, rfl, heq⟩ := exists_sub_of_mem_outs ho hnd
    obtain ⟨anc, hswf, -⟩ := wf_sub hnd hwf ht hs
    rw [heq]
    obtain ⟨v, outs⟩ := s
    simp only [DfsTree.WF] at hswf
    exact hswf.1

theorem ofForest_cls_bridge {v : Nat} {o : DfsOut} (ho : o ∈ (ofForest forest).outs v) :
    o.cls = .bridge ↔
      o.isTree = true ∧ ∀ l, l ≤ (ofForest forest).depth v → ¬(ofForest forest).Returns o.dest l := by
  cases o with
  | back e dest cls =>
    obtain ⟨hle, -, hcls, -⟩ := back_out_facts hnd hwf ho
    simp only [DfsOut.cls, DfsOut.isTree, Bool.false_eq_true, false_and, iff_false]
    rw [hcls, classify_eq_bridge_iff]; omega
  | tree e cls c =>
    obtain ⟨hcls, hret⟩ := tree_out_facts hnd hwf ho
    simp only [DfsOut.cls, DfsOut.isTree, DfsOut.dest, true_and]
    rw [hcls, classify_child_bridge_iff]
    constructor
    · intro h l hl hr; have := h l ((hret l).mpr hr); omega
    · intro h y hy; by_contra hc; exact h y (by omega) ((hret y).mp hy)

theorem ofForest_cls_component {v : Nat} {o : DfsOut} (ho : o ∈ (ofForest forest).outs v) :
    o.cls = .component ↔
      o.isTree = true ∧ (ofForest forest).Returns o.dest ((ofForest forest).depth v) ∧
        ∀ l, l < (ofForest forest).depth v → ¬(ofForest forest).Returns o.dest l := by
  cases o with
  | back e dest cls =>
    obtain ⟨-, -, hcls, -⟩ := back_out_facts hnd hwf ho
    simp only [DfsOut.cls, DfsOut.isTree, Bool.false_eq_true, false_and, iff_false]
    rw [hcls]; exact classify_false_ne_component _ _
  | tree e cls c =>
    obtain ⟨hcls, hret⟩ := tree_out_facts hnd hwf ho
    simp only [DfsOut.cls, DfsOut.isTree, DfsOut.dest, true_and]
    rw [hcls, classify_child_component_iff]
    constructor
    · rintro ⟨h1, h2⟩
      exact ⟨(hret _).mp h1, fun l hl hr => by have := h2 l ((hret l).mpr hr); omega⟩
    · rintro ⟨h1, h2⟩
      exact ⟨(hret _).mpr h1, fun y hy => by by_contra hc; exact h2 y (by omega) ((hret y).mp hy)⟩

theorem ofForest_cls_selfLoop {v : Nat} {o : DfsOut} (ho : o ∈ (ofForest forest).outs v) :
    o.cls = .selfLoop ↔ o.isTree = false ∧ o.dest = v := by
  cases o with
  | back e dest cls =>
    obtain ⟨-, -, hcls, h4⟩ := back_out_facts hnd hwf ho
    simp only [DfsOut.cls, DfsOut.isTree, DfsOut.dest, true_and]
    rw [hcls, classify_eq_selfLoop_iff]
    exact ⟨h4, fun h => by rw [h]⟩
  | tree e cls c =>
    obtain ⟨hcls, -⟩ := tree_out_facts hnd hwf ho
    simp only [DfsOut.cls, DfsOut.isTree, Bool.true_eq_false, false_and, iff_false]
    rw [hcls]; exact classify_true_ne_selfLoop _ _

theorem ofForest_cls_backEdge {v : Nat} {o : DfsOut} (ho : o ∈ (ofForest forest).outs v) (l : Nat) :
    o.cls = .ret l .backEdge ↔
      o.isTree = false ∧ o.dest ≠ v ∧ (ofForest forest).depth o.dest = l := by
  cases o with
  | back e dest cls =>
    obtain ⟨hle, -, hcls, h4⟩ := back_out_facts hnd hwf ho
    simp only [DfsOut.cls, DfsOut.isTree, DfsOut.dest, true_and]
    rw [hcls, classify_eq_ret_iff_back]
    constructor
    · rintro ⟨h1, h2, -⟩
      exact ⟨fun h => by subst h; omega, h1⟩
    · rintro ⟨h1, h2⟩
      refine ⟨h2, ?_, rfl⟩
      rcases Nat.lt_or_eq_of_le hle with h | h
      · omega
      · exact absurd (h4 h) h1
  | tree e cls c =>
    obtain ⟨hcls, -⟩ := tree_out_facts hnd hwf ho
    simp only [DfsOut.cls, DfsOut.isTree, Bool.true_eq_false, false_and, iff_false]
    rw [hcls]; exact classify_true_ne_backEdge _ _ _

theorem ofForest_cls_type1 {v : Nat} {o : DfsOut} (ho : o ∈ (ofForest forest).outs v) (l : Nat) :
    o.cls = .ret l .type1Child ↔
      o.isTree = true ∧ l < (ofForest forest).depth v ∧ (ofForest forest).Returns o.dest l ∧
        (∀ l', l' < l → ¬(ofForest forest).Returns o.dest l') ∧
        ∀ l', l < l' → l' < (ofForest forest).depth v → ¬(ofForest forest).Returns o.dest l' := by
  cases o with
  | back e dest cls =>
    obtain ⟨-, -, hcls, -⟩ := back_out_facts hnd hwf ho
    simp only [DfsOut.cls, DfsOut.isTree, Bool.false_eq_true, false_and, iff_false]
    rw [hcls, classify_eq_ret_iff_back]; simp
  | tree e cls c =>
    obtain ⟨hcls, hret⟩ := tree_out_facts hnd hwf ho
    simp only [DfsOut.cls, DfsOut.isTree, DfsOut.dest, true_and]
    rw [hcls, classify_child_type1_iff]
    constructor
    · rintro ⟨h1, h2, h3⟩
      refine ⟨h2, (hret _).mp h1, fun l' hl' hr => ?_, fun l' hl hl' hr => ?_⟩
      · have := h3 l' ((hret l').mpr hr); omega
      · have := h3 l' ((hret l').mpr hr); omega
    · rintro ⟨h1, h2, h3, h4⟩
      refine ⟨(hret _).mpr h2, h1, fun y hy => ?_⟩
      have hy' := (hret y).mp hy
      by_contra hc
      rcases Nat.lt_trichotomy y l with h | h | h
      · exact h3 y h hy'
      · exact hc (.inl h)
      · exact h4 y h (by omega) hy'

theorem ofForest_cls_type2 {v : Nat} {o : DfsOut} (ho : o ∈ (ofForest forest).outs v) (l : Nat) :
    o.cls = .ret l .type2Child ↔
      o.isTree = true ∧ l < (ofForest forest).depth v ∧ (ofForest forest).Returns o.dest l ∧
        (∀ l', l' < l → ¬(ofForest forest).Returns o.dest l') ∧
        ∃ l', l < l' ∧ l' < (ofForest forest).depth v ∧ (ofForest forest).Returns o.dest l' := by
  cases o with
  | back e dest cls =>
    obtain ⟨-, -, hcls, -⟩ := back_out_facts hnd hwf ho
    simp only [DfsOut.cls, DfsOut.isTree, Bool.false_eq_true, false_and, iff_false]
    rw [hcls, classify_eq_ret_iff_back]; simp
  | tree e cls c =>
    obtain ⟨hcls, hret⟩ := tree_out_facts hnd hwf ho
    simp only [DfsOut.cls, DfsOut.isTree, DfsOut.dest, true_and]
    rw [hcls, classify_child_type2_iff]
    constructor
    · rintro ⟨h1, h2, h3, y, hy, hyl, hyd⟩
      refine ⟨h2, (hret _).mp h1, fun l' hl' hr => ?_, y, ?_, hyd, (hret _).mp hy⟩
      · have := h3 l' ((hret l').mpr hr); omega
      · have := h3 y hy; omega
    · rintro ⟨h1, h2, h3, y, hyl, hyd, hy⟩
      refine ⟨(hret _).mpr h2, h1, fun y' hy' => ?_, y, (hret _).mpr hy, by omega, hyd⟩
      by_contra hc
      exact h3 y' (by omega) ((hret y').mp hy')

end WF

end DfsData

/-! ### Out-edges are adjacency entries -/

mutual
/-- Every out-edge `(e, dest)` of a vertex `v` of the tree is an entry `(dest, e)` of `adj[v]`. -/
def DfsTree.OutsIn (adj : Array (List (Nat × Nat))) : DfsTree → Prop
  | .node v outs => DfsOut.OutsInList adj v outs
def DfsOut.OutsInList (adj : Array (List (Nat × Nat))) (v : Nat) : List DfsOut → Prop
  | [] => True
  | .back e dest _ :: rest => (dest, e) ∈ adj[v]! ∧ DfsOut.OutsInList adj v rest
  | .tree e _ child :: rest =>
    (child.v, e) ∈ adj[v]! ∧ child.OutsIn adj ∧ DfsOut.OutsInList adj v rest
end

mutual
theorem DfsTree.OutsIn.mem {adj : Array (List (Nat × Nat))} : ∀ {t : DfsTree}, t.OutsIn adj →
    ∀ {x : Nat} {o : DfsOut}, (x, o) ∈ t.allOuts → (o.dest, o.e) ∈ adj[x]!
  | .node _ _, h, _, _, hmem => DfsOut.OutsInList.mem h hmem
theorem DfsOut.OutsInList.mem {adj : Array (List (Nat × Nat))} {v : Nat} :
    ∀ {outs : List DfsOut}, DfsOut.OutsInList adj v outs →
    ∀ {x : Nat} {o : DfsOut}, (x, o) ∈ DfsOut.allOutsList v outs → (o.dest, o.e) ∈ adj[x]!
  | [], _, _, _, hmem => absurd hmem List.not_mem_nil
  | .back e dest cls :: rest, h, x, o, hmem => by
    rcases List.mem_cons.mp hmem with heq | hmem'
    · obtain ⟨rfl, rfl⟩ := Prod.mk.inj heq
      exact h.1
    · exact DfsOut.OutsInList.mem h.2 hmem'
  | .tree e cls child :: rest, h, x, o, hmem => by
    rcases List.mem_cons.mp hmem with heq | hmem'
    · obtain ⟨rfl, rfl⟩ := Prod.mk.inj heq
      exact h.1
    · rcases List.mem_append.mp hmem' with hmem'' | hmem''
      · exact DfsTree.OutsIn.mem h.2.1 hmem''
      · exact DfsOut.OutsInList.mem h.2.2 hmem''
end

/-- The per-out-edge content of `DfsOut.OutsInList`. -/
def DfsOut.OutIn (adj : Array (List (Nat × Nat))) (v : Nat) (o : DfsOut) : Prop :=
  (o.dest, o.e) ∈ adj[v]! ∧ ∀ e cls child, o = .tree e cls child → child.OutsIn adj

theorem DfsOut.outsInList_iff {adj : Array (List (Nat × Nat))} {v : Nat} :
    ∀ {outs : List DfsOut}, DfsOut.OutsInList adj v outs ↔ ∀ o ∈ outs, DfsOut.OutIn adj v o
  | [] => by simp [DfsOut.OutsInList]
  | .back e dest cls :: rest => by
    simp only [DfsOut.OutsInList, DfsOut.outsInList_iff, List.forall_mem_cons]
    constructor
    · rintro ⟨h1, h2⟩; exact ⟨⟨h1, fun _ _ _ h => nomatch h⟩, h2⟩
    · rintro ⟨h1, h2⟩; exact ⟨h1.1, h2⟩
  | .tree e cls child :: rest => by
    simp only [DfsOut.OutsInList, DfsOut.outsInList_iff, List.forall_mem_cons]
    constructor
    · rintro ⟨h1, h2, h3⟩
      exact ⟨⟨h1, fun _ _ _ h => by cases h; exact h2⟩, h3⟩
    · rintro ⟨h1, h2⟩; exact ⟨h1.1, h1.2 _ _ _ rfl, h2⟩

theorem dfsVisit_v (adj : Array (List (Nat × Nat))) :
    ∀ (fuel v d : Nat) (prvE : Option Nat) (depth : Array (Option Nat)),
      (dfsVisit adj fuel v d prvE depth).1.v = v
  | 0, _, _, _, _ => rfl
  | _ + 1, _, _, _, _ => by rw [dfsVisit_succ]; rfl

theorem dfsVisit_outsIn (adj : Array (List (Nat × Nat))) :
    ∀ (fuel v d : Nat) (prvE : Option Nat) (depth : Array (Option Nat)),
      (dfsVisit adj fuel v d prvE depth).1.OutsIn adj
  | 0, v, _, _, _ => by show DfsOut.OutsInList adj v []; trivial
  | fuel + 1, v, d, prvE, depth => by
    rw [dfsVisit_succ]
    show DfsOut.OutsInList adj v _
    rw [DfsOut.outsInList_iff]
    intro o ho
    rw [List.mem_mergeSort, List.mem_reverse] at ho
    have key : ∀ (l : List (Nat × Nat)) (init : List DfsOut × Lowvals × Array (Option Nat)),
        (∀ p ∈ l, p ∈ adj[v]!) → (∀ o ∈ init.1, DfsOut.OutIn adj v o) →
        ∀ o ∈ (l.foldl (dfsStep adj fuel d prvE) init).1, DfsOut.OutIn adj v o := by
      intro l
      induction l with
      | nil => exact fun _ _ h => h
      | cons p l ih =>
        intro init hl hinit
        rw [List.foldl_cons]
        refine ih _ (fun p' hp' => hl p' (List.mem_cons_of_mem _ hp')) ?_
        obtain ⟨outs, lv, cur⟩ := init
        obtain ⟨nxt, e⟩ := p
        have hp : (nxt, e) ∈ adj[v]! := hl _ List.mem_cons_self
        simp only [dfsStep]
        split
        · exact hinit
        · split
          · obtain ⟨⟨child, n, cur'⟩, hr⟩ :
                ∃ r, dfsVisit adj fuel nxt (d + 1) (some e) cur = r := ⟨_, rfl⟩
            have hv : child.v = nxt := by
              have := dfsVisit_v adj fuel nxt (d + 1) (some e) cur
              rw [hr] at this; exact this
            have hin : child.OutsIn adj := by
              have := dfsVisit_outsIn adj fuel nxt (d + 1) (some e) cur
              rw [hr] at this; exact this
            rw [hr]
            intro o ho
            simp only [List.mem_cons] at ho
            rcases ho with rfl | ho
            · refine ⟨?_, fun _ _ _ h => by cases h; exact hin⟩
              simp only [DfsOut.dest, DfsOut.e]
              rw [hv]; exact hp
            · exact hinit o ho
          · intro o ho
            simp only [List.mem_cons] at ho
            rcases ho with rfl | ho
            · exact ⟨hp, fun _ _ _ h => nomatch h⟩
            · exact hinit o ho
    exact key adj[v]! _ (fun _ h => h) (fun _ h => absurd h List.not_mem_nil) o ho

theorem dfsForest_outsIn (g : Graph) (vo eo : List Nat) :
    ∀ t ∈ g.dfsForest vo eo, t.OutsIn (g.adjacency eo) := by
  rw [dfsForest_eq]
  have key : ∀ (l : List Nat) (init : List DfsTree × Array (Option Nat)),
      (∀ t ∈ init.1, t.OutsIn (g.adjacency eo)) →
      ∀ t ∈ (l.foldl (forestStep (g.adjacency eo) g.nv) init).1, t.OutsIn (g.adjacency eo) := by
    intro l
    induction l with
    | nil => exact fun _ h => h
    | cons rt l ih =>
      intro init hinit
      rw [List.foldl_cons]
      refine ih _ ?_
      obtain ⟨roots, depth⟩ := init
      simp only [forestStep]
      split
      · exact hinit
      · obtain ⟨⟨t, n, depth'⟩, hr⟩ :
            ∃ r, dfsVisit (g.adjacency eo) g.nv rt 0 none depth = r := ⟨_, rfl⟩
        have hin : t.OutsIn (g.adjacency eo) := by
          have := dfsVisit_outsIn (g.adjacency eo) g.nv rt 0 none depth
          rw [hr] at this; exact this
        rw [hr]
        intro t' ht'
        simp only [List.mem_cons] at ht'
        rcases ht' with rfl | ht'
        · exact hin
        · exact hinit t' ht'
  intro t ht
  exact key _ _ (fun _ h => absurd h List.not_mem_nil) t (List.mem_reverse.mp ht)

theorem adjacency_mem_iff {g : Graph} (hg : g.WF) {eo : List Nat} (heo : OrderOK g.ne eo)
    (x y e : Nat) :
    (y, e) ∈ (g.adjacency eo)[x]! ↔
      e < g.ne ∧ x < g.nv ∧ (g.edges[e]! = (x, y) ∨ g.edges[e]! = (y, x)) := by
  have hperm := inOrder_perm heo
  have hinv := adjInv_foldl hg (inOrder g.ne eo) [] _ (adjInv_nil g)
    (by simpa using hperm.nodup_iff.2 List.nodup_range) (fun e he => by simpa using hperm.subset he)
  simp only [List.nil_append] at hinv
  obtain ⟨_, hmem, -⟩ := hinv
  have hrev : (g.adjacency eo)[x]! =
      ((inOrder g.ne eo).foldl (adjStep g) (Array.replicate g.nv []))[x]!.reverse := by
    rw [adjacency_eq]; exact getElem!_map' List.reverse rfl _ _
  rw [hrev, List.mem_reverse, hmem]
  constructor
  · rintro ⟨_, h⟩; exact h
  · intro h; exact ⟨hperm.mem_iff.2 (List.mem_range.2 h.1), h⟩

/-! ### Main theorem -/

/-- The forest computed by phase 1 satisfies the hypotheses of `Spqr/SepPair.lean`. -/
theorem dfsForestSpec_of_dfsForest {g : Graph} (hg : g.WF) {vo eo : List Nat}
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) : DfsForestSpec g (g.dfsForest vo eo) := by
  obtain ⟨hvp, hep⟩ := dfsForest_spanning' hg hvo heo
  have hnd : ((g.dfsForest vo eo).flatMap DfsTree.verts).Nodup := hvp.nodup_iff.2 List.nodup_range
  have hend : ((g.dfsForest vo eo).flatMap DfsTree.edges).Nodup := hep.nodup_iff.2 List.nodup_range
  have hwf := dfsForest_wf hg hvo heo
  have hin := dfsForest_outsIn g vo eo
  refine
    { verts_nodup := hnd
      edges_nodup := DfsData.edgePostorderForest_perm.nodup_iff.2 hend
      joins := ?_
      edge_out := ?_
      out_inj := fun v v' o ho o' ho' he => DfsData.ofForest_out_inj hend v v' ho ho' he
      nodup := DfsData.ofForest_outs_nodup hnd hend
      tree_inj := fun v o ho o' ho' ht ht' hd =>
        (DfsData.ofForest_tree_dest_inj hnd ho ho' ht ht' hd).2
      back_anc := fun v o ho hb => DfsData.ofForest_back_anc hnd hwf ho hb
      parent_unique := fun p p' c ⟨o, ho, ht, hd⟩ ⟨o', ho', ht', hd'⟩ =>
        (DfsData.ofForest_tree_dest_inj hnd ho ho' ht ht' (hd.trans hd'.symm)).1
      depth_parent := fun p c h => DfsData.ofForest_depth_parent' hnd h
      depth_root := DfsData.ofForest_depth_root' hnd
      sorted := DfsData.ofForest_sorted hnd hwf
      cls_bridge := fun v o ho => DfsData.ofForest_cls_bridge hnd hwf ho
      cls_component := fun v o ho => DfsData.ofForest_cls_component hnd hwf ho
      cls_selfLoop := fun v o ho => DfsData.ofForest_cls_selfLoop hnd hwf ho
      cls_backEdge := fun v o ho l => DfsData.ofForest_cls_backEdge hnd hwf ho l
      cls_type1 := fun v o ho l => DfsData.ofForest_cls_type1 hnd hwf ho l
      cls_type2 := fun v o ho l => DfsData.ofForest_cls_type2 hnd hwf ho l }
  · intro v o ho
    obtain ⟨t, ht, hto⟩ := List.mem_flatMap.mp (DfsData.mem_forestAllOuts.mpr ho)
    obtain ⟨he, -, h⟩ := (adjacency_mem_iff hg heo _ _ _).mp ((hin t ht).mem hto)
    rw [getElem!_pos g.edges o.e he] at h
    show g.edges[o.e]? = some (v, o.dest) ∨ g.edges[o.e]? = some (o.dest, v)
    rw [getElem?_pos g.edges o.e he]
    rcases h with h | h <;> rw [h] <;> simp
  · intro e he
    have hmem : e ∈ ((g.dfsForest vo eo).flatMap DfsTree.allOuts).map (·.2.e) := by
      rw [DfsData.forestAllOuts_map_e]
      exact hep.mem_iff.2 (List.mem_range.2 he)
    obtain ⟨⟨v, o⟩, hvo', rfl⟩ := List.mem_map.mp hmem
    exact ⟨v, o, DfsData.mem_forestAllOuts.mp hvo', rfl⟩

end Spqr
