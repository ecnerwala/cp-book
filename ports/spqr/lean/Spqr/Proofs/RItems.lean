import Spqr.Proofs.RSkel3
import Spqr.Correctness
import Spqr.Proofs.RSkelRoot

/-!
# R items are 3-connected (PROOF.md §4.5, item level)

`Items.RSkel3 g items i`: the skeleton of the R item `i` (its non-`V` children as pieces, the
complement of its edge set as the parent piece at `i`'s terminals, contracted) is 3-connected.
`RBranch.rSkel3` establishes it for the item Loop 1's R branch closes (`rCloseItems`: `allocItem .R`,
then `finishTstackTop` of the merged entry), from `RBranch.threeConnected`. `items_r_three_connected`
(admitted) is the statement for every R item of the walk's output on a block.
-/

namespace Spqr

theorem Items.WF.q_root_covers_block {g : Graph} {items : Items} (hwf : items.WF g)
    (hg : g.WF) (h2 : g.TwoConnected) {e : Nat} (he : e < g.ne)
    (hch : items.ch (edgeItem g e) ≠ []) :
    ∀ f, f < g.ne → items.EdgeBelow g (edgeItem g e) f := by
  have hi : edgeItem g e < items.size := by
    have := hwf.tree.size
    simp only [edgeItem, ItemId] at *
    omega
  obtain ⟨u, hu⟩ := hwf.vs_fst hi (by simp [hwf.tree.edge e he, NodeType.isNode])
  have hv := (hwf.endpoints.q_vs e he u hu).2.1.2 hch
  have hatt : g.TwoAttached (items.EdgeBelow g (edgeItem g e)) u u := by
    intro v f f' hf hf' hE hE' hfv hf'v
    have hvl : v < g.nv := by
      rcases hfv with h | h
      · exact h ▸ (hg.getElem! hf).1
      · exact h ▸ (hg.getElem! hf).2
    have hvs := hwf.endpoints.separation (edgeItem g e) hi
      (by simp [hwf.tree.edge e he]) v f f' hvl hf hf' hfv hf'v hE hE'
    exact .inl (by simpa [hu, hv, eq_comm] using hvs)
  intro f hf
  exact hatt.edgeConn_mem (fun h => h rfl) (fun h => h rfl)
    Relation.ReflTransGen.refl (h2 u e f he hf)

theorem Items.WF.vertex_child_leaf_of_block {g : Graph} {items : Items} (hwf : items.WF g)
    (hg : g.WF) (h2 : g.TwoConnected) {i v : Nat} (hi : i < items.size) (hv : v < g.nv)
    (hn : items.type i ∈ [NodeType.S, .P, .R]) (hp : items.IsParent i (vertItem v))
    (hroot : ∀ c, items.IsParent (vertItem v) c → items.ch c ≠ []) :
    items.ch (vertItem v) = [] := by
  apply List.eq_nil_iff_forall_not_mem.2
  intro c hc
  have hQ := hwf.tree.v_children v c hv hc
  have hid := (hwf.tree.type_Q_iff (hwf.tree.ch_lt _ _ hc)).1 hQ
  have he : c - 1 - g.nv < g.ne := by unfold ItemId at *; omega
  have heq : edgeItem g (c - 1 - g.nv) = c := by unfold edgeItem ItemId at *; omega
  have hall := hwf.q_root_covers_block hg h2 he (by rw [heq]; exact hroot c hc)
  have hnot := ((hwf.endpoints.interior i v hi hv hn).1 hp).2.2 (vertItem v) hp
  apply hnot
  intro e he _
  exact Relation.ReflTransGen.head hc (heq ▸ hall e he)

theorem Items.WF.nonV_child_cover_of_block {g : Graph} {items : Items} (hwf : items.WF g)
    (hg : g.WF) (h2 : g.TwoConnected)
    (hroot : ∀ v, v < g.nv → ∀ c, items.IsParent (vertItem v) c → items.ch c ≠ [])
    {i : ItemId} (hi : i < items.size) (hn : items.type i ∈ [NodeType.S, .P, .R]) :
    ∀ e, e < g.ne → items.EdgeBelow g i e →
      ∃ c ∈ items.ch i, items.type c ≠ .V ∧ items.EdgeBelow g c e := by
  intro e he hb
  rcases Relation.ReflTransGen.cases_head hb with h | ⟨c, hc, hb⟩
  · rw [h, hwf.tree.edge e he] at hn
    simp at hn
  · refine ⟨c, hc, ?_, hb⟩
    intro hV
    have hid := (hwf.tree.type_V_iff (hwf.tree.ch_lt _ _ hc)).1 hV
    have hv : c - 1 < g.nv := by unfold ItemId at *; omega
    have heq : vertItem (c - 1) = c := by unfold vertItem ItemId at *; omega
    have hleaf := hwf.vertex_child_leaf_of_block hg h2 hi hv hn (by rwa [heq]) (hroot _ hv)
    rw [heq] at hleaf
    have hec : c = edgeItem g e := by
      rcases Relation.ReflTransGen.cases_head hb with h | ⟨j, hj, _⟩
      · exact h
      · simp [Items.IsParent, hleaf] at hj
    rw [hec, hwf.tree.edge e he] at hV
    cases hV

theorem Items.RSkel3.rThreeConnected {g : Graph} {items : Items} {i : ItemId}
    (h : Items.RSkel3 g items i) (hwf : items.WF g) (hi : i < items.size)
    (hR : items.type i = .R)
    (hcover : ∀ e, e < g.ne → items.EdgeBelow g i e →
      ∃ c ∈ items.ch i, items.type c ≠ .V ∧ items.EdgeBelow g c e) :
    SpqrTree.ThreeConnected (items.nvList g i).length (items.rSkeleton g i) := by
  classical
  obtain ⟨s, t, hvs, h3⟩ := h
  let P := (Pieces.ofItems g items ((items.ch i).filter fun c => decide (items.type c ≠ .V))).addParent g
    (items.EdgeBelow g i) s t
  obtain ⟨idx, _, hidx, hnode⟩ := relabel_node_spec g items hwf
  let hr : RelabelOK g items (relabelTree g items) idx := ⟨hwf, hidx, hnode⟩
  have hedges : (P.contract g).edges.toList = items.virtualEdges i ++ [(s, t)] :=
    Pieces.ofItems_addParent_edges s t (fun e he hE => by
      obtain ⟨c, hc, ht, hce⟩ := hcover e he hE
      exact ⟨c, List.mem_filter.2 ⟨hc, by simpa using ht⟩, hce⟩)
  have hcap : ∀ v, (items.vs i).1 = some v ∨ (items.vs i).2 = some v → v = s ∨ v = t := by
    intro v hv
    rcases hvs with hvs | hvs <;> simp only [hvs, Option.some.injEq] at hv <;> tauto
  have hcapmem : s ∈ items.nvList g i ∧ t ∈ items.nvList g i := by
    rcases hvs with hvs | hvs <;> simp [RelabelOK.nvList_eq, hvs]
  have hvirt : ∀ q ∈ items.virtualEdges i, q.1 ∈ items.nvList g i ∧ q.2 ∈ items.nvList g i := by
    intro q hq
    obtain ⟨c, hc, rfl⟩ := List.mem_map.1 hq
    obtain ⟨hc, hcV⟩ := List.mem_filter.1 hc
    obtain ⟨u, v, hcv⟩ := (hwf.shapes.r_shape i hi hR).2.2.2.2.2 c hc (by simpa using hcV)
    simp only [hcv, Option.getD_some]
    exact ⟨hr.mem_nvList_of_child (by simp [hR, NodeType.isNode]) hc (by simpa using hcV)
      (.inl (by rw [hcv])),
      hr.mem_nvList_of_child (by simp [hR, NodeType.isNode]) hc (by simpa using hcV)
      (.inr (by rw [hcv]))⟩
  have hends : ∀ e u v, (P.contract g).Joins e u v →
      u ∈ items.nvList g i ∧ v ∈ items.nvList g i := by
    intro e u v he
    rcases he with he | he
    all_goals
      have hm := Array.mem_of_getElem? he
      rw [← Array.mem_toList_iff, hedges] at hm
      rcases List.mem_append.1 hm with hm | hm
      · first | exact hvirt _ hm | exact (hvirt _ hm).symm
      · simp only [List.mem_singleton, Prod.mk.injEq] at hm
        rcases hm with ⟨rfl, rfl⟩
        first | exact hcapmem | exact hcapmem.symm
  have hactive : ∀ v ∈ items.nvList g i, ∃ e, (P.contract g).IsEnd e v := by
    intro v hv
    have hin : ∃ q ∈ (P.contract g).edges.toList, q.1 = v ∨ q.2 = v := by
      rw [RelabelOK.nvList_eq, hr.filter_lt_eq hi] at hv
      rcases List.mem_append.1 hv with hv | hv
      · rcases List.mem_append.1 hv with hv | hv
        · exact ⟨(s, t), by simp [hedges],
            (hcap v (.inl (by simpa using hv))).imp Eq.symm Eq.symm⟩
        · obtain ⟨c, hc, rfl⟩ := List.mem_map.1 hv
          obtain ⟨hc, hcV⟩ := List.mem_filter.1 hc
          have hclt := (hr.type_V_iff hi hc).1 (by simpa using hcV)
          obtain ⟨heq, _⟩ := hr.ch_lt_vert hc hclt
          have hvlt : c - 1 < g.nv := by
            have hc0 : c ≠ 0 := hr.ch_ne_root hc
            have hcLe : c ≤ g.nv := Nat.lt_succ_iff.1 (by simpa [Nat.add_comm] using hclt)
            omega
          have hvpar : items.IsParent i (vertItem (c - 1)) := by rwa [← heq]
          obtain ⟨⟨e, he, hev⟩, hall, hnot⟩ :=
            (hwf.endpoints.interior i (c - 1) hi hvlt (by simp [hR])).1 hvpar
          obtain ⟨j, hj, hjV, hje⟩ := hcover e he (hall e he hev)
          have hother := hnot j hj
          push Not at hother
          obtain ⟨f, hf, hfv, hjf⟩ := hother
          have hterm := hwf.endpoints.separation j (hwf.tree.ch_lt i j hj)
            (by simp [hr.child_type_ne_F hi hj, hjV])
            (c - 1) e f hvlt he hf hev hfv hje hjf
          obtain ⟨x, y, hjvs⟩ := (hwf.shapes.r_shape i hi hR).2.2.2.2.2 j hj hjV
          refine ⟨(x, y), ?_, ?_⟩
          · rw [hedges]
            apply List.mem_append_left
            exact List.mem_map.2 ⟨j, List.mem_filter.2 ⟨hj, by simpa using hjV⟩, by simp [hjvs]⟩
          · simpa [hjvs] using hterm
      · exact ⟨(s, t), by simp [hedges],
          (hcap v (.inr (by simpa using hv))).imp Eq.symm Eq.symm⟩
    obtain ⟨⟨u, w⟩, hq, hqv⟩ := hin
    obtain ⟨e, he⟩ := Graph.exists_joins_of_mem hq
    rcases hqv with rfl | rfl
    · exact ⟨e, he.isEnd⟩
    · exact ⟨e, he.symm.isEnd⟩
  have hne : 3 ≤ (P.contract g).ne := by
    have hlen := congrArg List.length hedges
    have := (hwf.shapes.r_shape i hi hR).2.1
    simp only [List.length_append, List.length_cons, List.length_nil, Array.length_toList] at hlen
    change 3 ≤ (P.contract g).edges.size
    omega
  have hcut := h3.relabel (h3.twoConnected hne) (hr.nv_nodup hi)
    (hr.nvList_R hi hR).choose_spec.choose_spec.2 hactive hends
  apply (SpqrTree.ThreeConnected_congr_undirected (es' := items.rSkeleton g i) ?_).1 hcut
  intro u v
  rw [hedges]
  rcases hvs with hvs | hvs <;>
    simp only [List.map_append, Items.rSkeleton, hvs, Option.getD_some,
      List.map_cons, List.mem_append, List.mem_cons, Prod.mk.injEq] <;> tauto

theorem Items.RSkel3.rThreeConnected_of_block {g : Graph} {items : Items} {i : ItemId}
    (h : Items.RSkel3 g items i) (hwf : items.WF g) (hg : g.WF) (h2 : g.TwoConnected)
    (hroot : ∀ v, v < g.nv → ∀ c, items.IsParent (vertItem v) c → items.ch c ≠ [])
    (hi : i < items.size) (hR : items.type i = .R) :
    SpqrTree.ThreeConnected (items.nvList g i).length (items.rSkeleton g i) :=
  h.rThreeConnected hwf hi hR (hwf.nonV_child_cover_of_block hg h2 hroot hi (by simp [hR]))


/-- Admitted (PROOF.md §4.5, item level): on a block, every R item of the walk's output has a
3-connected skeleton. Every R item is closed either at Loop 1's R branch (`RBranch.rSkel3` under
`loop1_rBranch`) or as the type-1 R close of `finishEdge` (the `isSingle = false` case, the same
argument with `cur` the type-1 entry); the skeleton is then untouched by later steps (items are
only modified when allocated or closed). -/
theorem items_r_three_connected (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) (h2 : g.TwoConnected) :
    let items := (g.walk tern (g.dfsForest vo eo)).items
    ∀ i, i < items.size → Items.type items i = NodeType.R → Items.RSkel3 g items i := by
  intro items i hi hty
  exact walk_rSkelInv g hg tern vo eo hvo heo h2 i hi hty

/-- On a block, the output contract `Items.RThreeConnected` for the walk's items: `items_r_three_connected`
bridged by `Items.RSkel3.rThreeConnected_of_block`, whose placement input (`Q` children of `V` items
have children) is `Items.Ranges.q_under_v` from `walk_ranges`. -/
theorem walk_items_rThreeConnected_of_twoConnected (g : Graph) (hg : g.WF) (tern : Bool)
    (vo eo : List Nat) (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) (h2 : g.TwoConnected) :
    Items.RThreeConnected g (g.walk tern (g.dfsForest vo eo)).items := by
  intro i hi hR
  have hwf := walk_items_wf g hg tern vo eo hvo heo
  refine (items_r_three_connected g hg tern vo eo hvo heo h2 i hi hR).rThreeConnected_of_block
    hwf hg h2 ?_ hi hR
  intro v hv c hc
  have hnv : 0 < g.nv := by omega
  obtain ⟨hvp, hep⟩ := dfsForest_spanning' hg hvo heo
  have hcov : ∀ e, e < g.ne → ∃ t ∈ g.dfsForest vo eo, e ∈ t.edges := fun e he => by
    simpa [List.mem_flatMap] using hep.mem_iff.2 (List.mem_range.2 he)
  exact (walk_ranges g tern vo eo hg hvo heo hnv (dfsForest_bounded g hg hvo heo) hcov).q_under_v
    v c hv hc

end Spqr
