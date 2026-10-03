import Spqr.EarWalk
import Spqr.WalkWF
import Spqr.RangesWalk
import Spqr.RangesClose
import Spqr.RangesCloseSites
import Spqr.RangesCloseTree

/-! # `Items.WF` for the walk on a DFS forest

`walk_items_wf` assembles `Items.WF` for `g.walk tern (g.dfsForest vo eo)` from the proved tree
(`walk_tree`) and typing (`WalkTyping.walk_typing`) layers and the range layer
(`wf_of_ranges` with `walk_ranges`, from the range invariant and `walk_closeFacts`). It lives above `WalkCover` (which imports
`StWalk`, which imports `WalkWF`), so the DFS prerequisites of `walk_tree`/`walk_typing` are
discharged here from `g.WF` and `OrderOK`. -/

namespace Spqr

theorem WalkTree.toTree {g : Graph} {items : Items} (h : WalkTree g items) : items.Tree g :=
  { size := h.size, root := h.root, vert := h.vert, edge := h.edge, node := h.node, ch_lt := h.ch_lt,
    unique_parent := h.unique_parent, root_no_parent := h.root_no_parent, ch_nodup := h.ch_nodup,
    reach := h.reach, v_children := h.v_children, root_children := h.root_children }

theorem WalkTyping.toTypingFacts {g : Graph} {items : Items} (h : WalkTyping g items) :
    Items.TypingFacts g items :=
  { i_o_leaf := h.i_o_leaf, vs_shape := h.vs_shape, vs_lt := h.vs_lt }

/-! ### `DfsTree.Bounded` from `DfsTree.WF` and vertex/edge ranges -/

mutual
theorem DfsTree.bounded_of_wf {nv ne : Nat} :
    ∀ (anc : List Nat) (t : DfsTree), (∀ x ∈ anc, x < nv) → t.WF anc →
      (∀ v ∈ t.verts, v < nv) → (∀ e ∈ t.edges, e < ne) → t.Bounded nv ne
  | anc, .node v outs, hanc, hwf, hv, he => by
    have hv0 : v < nv := hv v (by simp [DfsTree.verts])
    have hwf' : ∀ o ∈ outs, o.WF anc v := by unfold DfsTree.WF at hwf; exact hwf.2
    refine ⟨hv0, DfsOut.boundedList_of_wf anc v outs hanc hv0 hwf' ?_ ?_⟩
    · intro w hw; exact hv w (by simp [DfsTree.verts, hw])
    · intro e he'; exact he e (by simp [DfsTree.edges, he'])
theorem DfsOut.boundedList_of_wf {nv ne : Nat} :
    ∀ (anc : List Nat) (v : Nat) (outs : List DfsOut), (∀ x ∈ anc, x < nv) → v < nv →
      (∀ o ∈ outs, o.WF anc v) → (∀ w ∈ DfsOut.vertsList outs, w < nv) →
      (∀ e ∈ DfsOut.edgesList outs, e < ne) → DfsOut.BoundedList nv ne outs
  | _, _, [], _, _, _, _, _ => trivial
  | anc, v, .back e dest cls :: rest, hanc, hv, hwf, hvs, hes => by
    refine ⟨⟨hes e (by simp [DfsOut.edgesList]), ?_⟩,
      DfsOut.boundedList_of_wf anc v rest hanc hv (fun o ho => hwf o (by simp [ho]))
        (by simpa [DfsOut.vertsList] using hvs) (fun e' he' => hes e' (by simp [DfsOut.edgesList, he']))⟩
    have hw := hwf (.back e dest cls) (by simp)
    unfold DfsOut.WF at hw
    obtain ⟨i, hi, -⟩ := hw
    have hmem : dest ∈ anc ++ [v] := List.mem_of_getElem? hi
    rcases List.mem_append.1 hmem with h | h
    · exact hanc _ h
    · simp at h; exact h ▸ hv
  | anc, v, .tree e cls child :: rest, hanc, hv, hwf, hvs, hes => by
    refine ⟨⟨hes e (by simp [DfsOut.edgesList]), ?_⟩,
      DfsOut.boundedList_of_wf anc v rest hanc hv (fun o ho => hwf o (by simp [ho]))
        (fun w hw => hvs w (by simp [DfsOut.vertsList, hw]))
        (fun e' he' => hes e' (by simp [DfsOut.edgesList, he']))⟩
    have hw := hwf (.tree e cls child) (by simp)
    unfold DfsOut.WF at hw
    refine DfsTree.bounded_of_wf (anc ++ [v]) child ?_ hw.1
      (fun w hw => hvs w (by simp [DfsOut.vertsList, hw]))
      (fun e' he' => hes e' (by simp [DfsOut.edgesList, he']))
    intro x hx
    rcases List.mem_append.1 hx with h | h
    · exact hanc _ h
    · simp at h; exact h ▸ hv
end

theorem dfsForest_bounded (g : Graph) {vo eo : List Nat} (hg : g.WF) (hvo : OrderOK g.nv vo)
    (heo : OrderOK g.ne eo) : ∀ t ∈ g.dfsForest vo eo, t.Bounded g.nv g.ne := by
  intro t ht
  obtain ⟨hvp, hep⟩ := dfsForest_spanning' hg hvo heo
  refine DfsTree.bounded_of_wf [] t (by simp) (dfsForest_wf hg hvo heo t ht) ?_ ?_
  · intro v hv
    exact List.mem_range.1 (hvp.subset (List.mem_flatMap.2 ⟨t, ht, hv⟩))
  · intro e he
    exact List.mem_range.1 (hep.subset (List.mem_flatMap.2 ⟨t, ht, he⟩))

/-! ### The empty graph: the walk is the initial items -/

theorem dfsForest_nil_of_nv_zero (g : Graph) {vo eo : List Nat} (hnv : g.nv = 0)
    (hvo : OrderOK g.nv vo) : g.dfsForest vo eo = [] := by
  have hvo' : vo = [] := List.eq_nil_iff_forall_not_mem.2 fun x hx => by
    have := hvo.2 x hx; omega
  rw [dfsForest_eq]; simp [inOrder, hvo', hnv]

theorem initialItems_facts (g : Graph) (hnv : g.nv = 0) (hne : g.ne = 0) (i : ItemId) :
    Items.type (initialItems g) i = .F ∧ Items.ch (initialItems g) i = [] ∧
      Items.vs (initialItems g) i = (none, none) := by
  have : initialItems g = #[⟨.F, (none, none), []⟩] := by simp [initialItems, hnv, hne]
  rw [this]
  rcases i with _ | i <;> simp [Items.type, Items.ch, Items.vs]

theorem wf_initialItems (g : Graph) (hnv : g.nv = 0) (hne : g.ne = 0) :
    Items.WF g (initialItems g) := by
  have hit := initialItems_facts g hnv hne
  have hnp : ∀ p c, ¬ Items.IsParent (initialItems g) p c := fun p c hc => by
    simp [Items.IsParent, (hit p).2.1] at hc
  have hsz : (initialItems g).size = 1 := by simp [initialItems, hnv, hne]
  have hF : ∀ i, Items.type (initialItems g) i ∉ [NodeType.F, .V] → False := fun i h => by
    simp [(hit i).1] at h
  have hSPR : ∀ i, Items.type (initialItems g) i ∈ [NodeType.S, .P, .R] → False := fun i h => by
    simp [(hit i).1] at h
  refine Items.wf_of_ranges (σ := []) ?_ ?_ ?_
  · exact {
      size := by rw [hsz]; omega
      root := (hit _).1
      vert := fun v hv => absurd hv (by omega)
      edge := fun e he => absurd he (by omega)
      node := fun i h1 h2 => absurd h2 (by omega)
      ch_lt := fun p c h => (hnp p c h).elim
      unique_parent := fun c h1 h2 => absurd h2 (by omega)
      root_no_parent := fun p h => hnp p _ h
      ch_nodup := fun p => by rw [(hit p).2.1]; exact List.nodup_nil
      reach := fun i hi => by
        rw [hsz] at hi
        obtain rfl : i = 0 := by omega
        exact Relation.ReflTransGen.refl
      v_children := fun v c hv _ => absurd hv (by omega)
      root_children := fun c h => (hnp _ c h).elim }
  · exact {
      i_o_leaf := fun i _ h => by simp [(hit i).1] at h
      vs_shape := fun i _ => by rw [(hit i).1]; exact (hit i).2.2
      vs_lt := fun i u h => by simp [(hit i).2.2] at h }
  · exact {
      convex := fun i _ h => (hF i h).elim
      att_vs := fun i _ h => (hF i h).elim
      vs_att := fun i _ h => by
        rcases h with h | ⟨h, -⟩
        · exact (hSPR i h).elim
        · simp [(hit i).1] at h
      vs_ne := fun i _ h => (hSPR i h).elim
      interior := fun i v _ hv => absurd hv (by omega)
      child_two := fun p c h => (hnp p c h).elim
      io_parent := fun p c h => (hnp p c h).elim
      q_leaf := fun e he => absurd he (by omega)
      q_root := fun e he => absurd he (by omega)
      q_under_v := fun v c hv => absurd hv (by omega)
      p_shape := fun i _ h => by simp [(hit i).1] at h
      s_order := fun i _ h => by simp [(hit i).1] at h
      r_shape := fun i _ h => by simp [(hit i).1] at h }

/-! ### `Items.Ranges` -/

/-- Admitted: at each actual P site, stack ownership covers the interval from the settled
base piece through the child's postorder block (`FinishPCover`); at each unpushed vertex,
its descendants have already been processed (`PushVertR`). The latter follows from `Place`
once its pushed-edge predicate is bounded by the current postorder prefix. This bookkeeping
is proved for each back/tree `CoverOut` by `coverOut_back` / `coverOut_tree`. The remaining
P ownership and enclosing `CoverOuts` / `CoverTree` / `RootsCover` induction are not proved;
no adjacency premise remains. -/
theorem walk_rootsCover (g : Graph) (tern : Bool) (vo eo : List Nat) (hg : g.WF)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    WalkState.RootsCover (edgePostorderForest (g.dfsForest vo eo)) 0
      (g.dfsForest vo eo) (WalkState.init g tern) := by
  sorry

theorem walk_rangesInv (g : Graph) (tern : Bool) (vo eo : List Nat) (hg : g.WF)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    (g.walk tern (g.dfsForest vo eo)).RangesInv (edgePostorderForest (g.dfsForest vo eo))
      (edgePostorderForest (g.dfsForest vo eo)).length 0 := by
  obtain ⟨hvp, hep⟩ := dfsForest_spanning' hg hvo heo
  exact WalkState.walk_rangesInv_of_cover g tern _ (ForestOK.of_perm hvp hep)
    (dfsForest_wf hg hvo heo) (dfsForest_ends g hg hvo heo)
    (fun e he => hep.symm.subset (List.mem_range.2 he)) (walk_rootsCover g tern vo eo hg hvo heo)

theorem walk_g (g : Graph) (tern : Bool) (forest : List DfsTree) (hnv : 0 < g.nv)
    (hb : ∀ t ∈ forest, t.Bounded g.nv g.ne) : (g.walk tern forest).g = g :=
  (WalkM.walkForest_typing forest (WalkState.init_typing g tern hnv) hb).g_eq

/-- `CloseInv` of the forest walk: `walk_closeInv'` runs the six close sites through the walk/forest
induction; the DFS-side site facts are proved (`dsTree`), the remaining per-site admissions are the
named `closeCtx_bd_vert`/`bd_node`/`p_site`/`v_site`/`l1_site` (RangesCloseTree.lean), plus
`walk_rootsCover` for the schedule. -/
theorem walk_closeInv (g : Graph) (tern : Bool) (vo eo : List Nat) (hg : g.WF)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    (g.walk tern (g.dfsForest vo eo)).CloseInv := by
  obtain ⟨hvp, hep⟩ := dfsForest_spanning' hg hvo heo
  exact WalkState.walk_closeInv' g tern _ (ForestOK.of_perm hvp hep)
    (dfsForest_wf hg hvo heo) (dfsForest_ends g hg hvo heo) (walk_rootsCover g tern vo eo hg hvo heo)

theorem walk_closeFacts (g : Graph) (tern : Bool) (vo eo : List Nat) (hg : g.WF)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) (hnv : 0 < g.nv) :
    Items.CloseFacts g (g.walk tern (g.dfsForest vo eo)).items := by
  obtain ⟨hvp, hep⟩ := dfsForest_spanning' hg hvo heo
  have hb := dfsForest_bounded g hg hvo heo
  have hvcov : ∀ v, v < g.nv → v ∈ (g.dfsForest vo eo).flatMap DfsTree.verts :=
    fun v hv => hvp.mem_iff.2 (List.mem_range.2 hv)
  have hecov : ∀ e, e < g.ne → e ∈ (g.dfsForest vo eo).flatMap DfsTree.edges :=
    fun e he => hep.mem_iff.2 (List.mem_range.2 he)
  have hcov : ∀ e, e < g.ne → ∃ t ∈ g.dfsForest vo eo, e ∈ t.edges :=
    fun e he => by simpa [List.mem_flatMap] using hecov e he
  have ht := walk_tree g tern _ hnv hb (ForestOK.of_perm hvp hep) (dfsForest_wf hg hvo heo)
    (dfsForest_ends g hg hvo heo) hvcov hecov
  have heq := walk_g g tern _ hnv hb
  have hc := (walk_closeInv g tern vo eo hg hvo heo).of_tree
    (by rw [heq]; exact ht.toTree) (by rw [heq]; exact walk_typing g tern _ hnv hb hcov)
  rwa [heq] at hc

/-- `Items.Ranges` for the walk: `convex`/node `att_vs` from the final range invariant
(`walk_rangesInv`), the remaining attachment/shape clauses from `walk_closeFacts`. -/
theorem walk_ranges (g : Graph) (tern : Bool) (vo eo : List Nat) (hgf : g.WF)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) (hnv : 0 < g.nv)
    (hb : ∀ t ∈ g.dfsForest vo eo, t.Bounded g.nv g.ne)
    (hcov : ∀ e, e < g.ne → ∃ t ∈ g.dfsForest vo eo, e ∈ t.edges) :
    Items.Ranges g (g.walk tern (g.dfsForest vo eo)).items (edgePostorderForest (g.dfsForest vo eo)) := by
  have hg := walk_g g tern _ hnv hb
  have key := WalkState.ranges_of_rangesInv (walk_rangesInv g tern vo eo hgf hvo heo)
    (by rw [hg]; exact walk_typing g tern _ hnv hb hcov)
    (by rw [hg]; exact walk_closeFacts g tern vo eo hgf hvo heo hnv)
  rwa [hg] at key

/-! ### `Items.WF` -/

/-- Phase 2: the walk's items satisfy the item-level specification. `Tree` is `walk_tree`,
`Endpoints`/`Shapes` are derived from `walk_ranges` by `Items.wf_of_ranges` (pure item-level
reasoning, no tstack facts); the empty graph is handled separately since `walk_tree` needs
`0 < g.nv`. -/
theorem walk_items_wf (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    Items.WF g (g.walk tern (g.dfsForest vo eo)).items := by
  obtain ⟨hvp, hep⟩ := dfsForest_spanning' hg hvo heo
  have hb := dfsForest_bounded g hg hvo heo
  have hwf := dfsForest_wf hg hvo heo
  have hvcov : ∀ v, v < g.nv → v ∈ (g.dfsForest vo eo).flatMap DfsTree.verts :=
    fun v hv => hvp.mem_iff.2 (List.mem_range.2 hv)
  have hecov : ∀ e, e < g.ne → e ∈ (g.dfsForest vo eo).flatMap DfsTree.edges :=
    fun e he => hep.mem_iff.2 (List.mem_range.2 he)
  rcases Nat.eq_zero_or_pos g.nv with hnv | hnv
  · have hne : g.ne = 0 := by
      by_contra hne
      have hpos : 0 < g.edges.size := Nat.pos_of_ne_zero hne
      have := (hg _ (Array.getElem_mem hpos)).1
      omega
    rw [dfsForest_nil_of_nv_zero g hnv hvo]
    exact wf_initialItems g hnv hne
  · have hf : ForestOK g (g.dfsForest vo eo) := ForestOK.of_perm hvp hep
    have ht := walk_tree g tern _ hnv hb hf hwf (dfsForest_ends g hg hvo heo) hvcov hecov
    have hcov : ∀ e, e < g.ne → ∃ t ∈ g.dfsForest vo eo, e ∈ t.edges :=
      fun e he => by simpa [List.mem_flatMap] using hecov e he
    have hty := walk_typing g tern _ hnv hb hcov
    exact Items.wf_of_ranges ht.toTree hty.toTypingFacts (walk_ranges g tern vo eo hg hvo heo hnv hb hcov)

end Spqr
