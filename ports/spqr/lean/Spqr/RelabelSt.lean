import Spqr.RelabelRep
import Spqr.StSpec
import Spqr.RelabelAdj

/-!
# Relabel phase D: st-order transport (`relabel_st'`)

`Items.StNumbered` (per-item st-numbering of the S/P/R item vertex lists, `Spqr/StSpec.lean`)
is transported along the per-node relabel interface (`relabel_node_spec`, `relabelTree_adj`) to
`SpqrTree.StOrder` of `relabelTree g items`:

* `StOrder.st`: the skeleton of each S/P/R node is st-numbered (`skeleton_S/P/R` give the skeleton
  lists from `LayoutFacts`; for R the virtual edges are `(pos x, pos y)` with
  `pos = nvSt + idxOf`, so `Items.StList` transfers verbatim);
* `StOrder.dom`: edge dominance (`ownEdges` of S/P/R from the lists above; small nodes trivially);
* `StOrder.adj`: the per-node-vertex adjacency rows, via the shared row transport
  `adjRow_eq_layout` (global CSR row = `LayoutR.row` of `layoutNode`, from
  `RelabelLayout.adj_bounds/adj_dat`, `Layout.Local`, and `relabelTree_adj`) and the per-type row
  descriptions of `LayoutShape` / `LayoutR.layoutNode_R_bracket`.

See `PROOF.md` §7.3.
-/

namespace Spqr

open LayoutR (row rowBound)

theorem List.pairwise_of_length_le_one {α : Type} {R : α → α → Prop} :
    ∀ (l : List α), l.length ≤ 1 → l.Pairwise R
  | [], _ => List.Pairwise.nil
  | [_], _ => List.pairwise_singleton _ _
  | _ :: _ :: _, h => by simp at h

@[simp] theorem pairwise_true {α : Type} : ∀ (l : List α), l.Pairwise fun _ _ => True
  | [] => List.Pairwise.nil
  | _ :: l => List.Pairwise.cons (fun _ _ => trivial) (pairwise_true l)

/-! ### `rowBound` of a local layout stays inside `[2 neSt, 2 neEn]` -/

theorem Layout.Local.rowBound_ge {node nvSt nvEn neSt neEn : Nat} {L : Layout}
    (hL : L.Local node nvSt nvEn neSt neEn) :
    ∀ r, 2 * nvSt ≤ r → r ≤ 2 * nvEn → 2 * neSt ≤ rowBound nvSt neSt L r := by
  intro r hr1
  induction r, hr1 using Nat.le_induction with
  | base => intro; simp [rowBound]
  | succ r hr ih => intro h2; exact (ih (by omega)).trans (hL.bound_mono r hr (by omega))

theorem Layout.Local.rowBound_le {node nvSt nvEn neSt neEn : Nat} {L : Layout}
    (hL : L.Local node nvSt nvEn neSt neEn) :
    ∀ r, 2 * nvSt ≤ r → r ≤ 2 * nvEn → rowBound nvSt neSt L r ≤ 2 * neEn := by
  intro r hr1 hr2
  have key : ∀ d, d ≤ 2 * nvEn - 2 * nvSt → rowBound nvSt neSt L (2 * nvEn - d) ≤ 2 * neEn := by
    intro d
    induction d with
    | zero => intro; simp only [Nat.sub_zero]; exact hL.bound_last.le
    | succ d ih =>
      intro hd
      have := hL.bound_mono (2 * nvEn - (d + 1)) (by omega) (by omega)
      rw [show 2 * nvEn - (d + 1) + 1 = 2 * nvEn - d by omega] at this
      exact this.trans (ih (by omega))
  have := key (2 * nvEn - r) (by omega)
  rwa [show 2 * nvEn - (2 * nvEn - r) = r by omega] at this

/-! ### Sortedness of `ordered` for R nodes -/

theorem Items.ordered_pairwise_loc (g : Graph) (items : Items) (a : ItemId) (nvSt : Nat)
    (pos : Nat → Nat) (hR : items.type a = .R) :
    (items.ordered g a nvSt pos).Pairwise fun x y =>
      items.loc g nvSt pos x ≤ items.loc g nvSt pos y := by
  rw [Items.ordered, LayoutShape.ite_of_neg (not_not.2 hR)]
  have := List.pairwise_mergeSort
    (le := fun x y => decide (items.loc g nvSt pos x ≤ items.loc g nvSt pos y))
    (fun a b c => by simp only [decide_eq_true_eq]; omega)
    (fun a b => by simp only [Bool.or_eq_true, decide_eq_true_eq]; omega) (items.ch a)
  exact this.imp fun h => by simpa using h

namespace RelabelOK

variable {g : Graph} {items : Items} {t : SpqrTree} {idx : ItemId → Nat}
  (h : RelabelOK g items t idx)
include h

theorem vertList_eq {i : ItemId} (hi : i < items.size) : items.vertList i = items.nvList g i := by
  rw [Items.vertList, nvList_eq, h.filter_lt_eq hi]

/-- `Items.StNumbered` orients the R nodes (`relabelTree_adj`'s hypothesis). -/
theorem rOriented (hst : items.StNumbered) : items.ROriented g := by
  intro i hi hR
  obtain ⟨u, v, -, ⟨hnd, -, -⟩, hor⟩ := hst i hi (Or.inr (Or.inr hR))
  rw [h.vertList_eq hi] at hnd hor
  exact ⟨hnd, hor⟩

/-- A V item has no node-vertices. -/
theorem nvList_V {i : ItemId} (hi : i < items.size) (hV : items.type i = .V) :
    items.nvList g i = [] := by
  have hs := h.endpoints.vs_shape i hi
  rw [hV] at hs
  have hi' : ∃ v, v < g.nv ∧ i = vertItem v := by
    by_cases h0 : i = 0
    · subst h0; exact absurd (hV.symm.trans h.tree.root) (by decide)
    by_cases h1 : i < 1 + g.nv
    · exact ⟨i - 1, by riomega, by riomega⟩
    by_cases h2 : i < 1 + g.nv + g.ne
    · have : i = edgeItem g (i - (1 + g.nv)) := by riomega
      rw [this, h.type_edgeItem (by riomega)] at hV; cases hV
    · have := h.tree.node i (by riomega) hi; simp [hV] at this
  obtain ⟨v, hv, rfl⟩ := hi'
  rw [nvList_eq, hs, h.filter_V_eq_nil hi]
  · rfl
  · intro c hc
    rw [h.tree.v_children v c hv hc]; decide

theorem nEdges_IO {a : ItemId} (ha : a < items.size) (hty : items.type a = .I ∨ items.type a = .O) :
    items.nEdges g a = 1 := by
  have hn : (items.type a).isNode = true := by rcases hty with h | h <;> rw [h] <;> rfl
  rw [h.nEdges_node ha hn, Items.virtualEdges, h.shapes.i_o_leaf a ha hty]
  rcases hty with h | h <;> simp [Items.capCount, Items.hasCap, h, NodeType.isNode]

/-! ### Skeleton lists -/

omit h in
theorem skeleton_eq_of (n : Nat) (L : List (Nat × Nat))
    (hlen : (t.neRange n).2 - (t.neRange n).1 = L.length)
    (hk : ∀ k (hk : k < L.length), (t.nodeEdges[(t.neRange n).1 + k]!).nvs = L[k]) :
    t.skeleton n = L := by
  apply List.ext_getElem
  · simp [SpqrTree.skeleton, SpqrTree.nodeEdgesOf, hlen]
  · intro k h1 h2
    simp only [SpqrTree.skeleton, SpqrTree.nodeEdgesOf, List.getElem_map, List.getElem_range]
    rw [← Array.getElem!_eq_getD]
    exact hk k h2

omit h in
theorem ownEdges_nvs (n : Nat) :
    (t.ownEdges n).map (·.nvs) = (t.skeleton n).drop (if t.hasCap n then 1 else 0) := by
  rw [SpqrTree.ownEdges, SpqrTree.skeleton, List.map_drop]

omit h in
theorem skeleton_length (n : Nat) : (t.skeleton n).length = (t.neRange n).2 - (t.neRange n).1 := by
  simp [SpqrTree.skeleton, SpqrTree.nodeEdgesOf]

theorem nvList_S' {a : ItemId} (ha : a < items.size) (hS : items.type a = .S) :
    ∃ u v xs, items.vs a = (some u, some v) ∧ items.nvList g a = u :: (xs ++ [v]) ∧
      1 ≤ xs.length ∧ (items.virtualEdges a).length = xs.length + 1 := by
  obtain ⟨u, v, xs, hvs, hxs, h2, hperm⟩ := h.shapes.s_shape a ha hS
  refine ⟨u, v, xs, hvs, ?_, h2, ?_⟩
  · rw [nvList_eq, hvs, h.filter_lt_eq ha, hxs, List.map_map]
    have : ((fun c => c - 1) ∘ vertItem) = id := by
      funext x; simp only [Function.comp, vertItem, id]; omega
    rw [this, List.map_id]; rfl
  · rw [hperm.length_eq, List.length_zip]; simp

theorem nEdges_S {a : ItemId} (ha : a < items.size) (hS : items.type a = .S) :
    items.nEdges g a = (items.nvList g a).length := by
  obtain ⟨u, v, xs, -, hnl, -, hvl⟩ := h.nvList_S' ha hS
  rw [h.nEdges_node ha (by rw [hS]; rfl), hvl, hnl]
  simp [Items.capCount, Items.hasCap, hS, NodeType.isNode]

theorem skeleton_S {a : ItemId} (ha : a < items.size) (hS : items.type a = .S) {pos : Nat → Nat}
    (hl : RelabelLayout g items t idx a pos) :
    t.skeleton (idx a) =
      ((t.nvRange (idx a)).1, (t.nvRange (idx a)).1 + ((items.nvList g a).length - 1)) ::
        (List.range ((items.nvList g a).length - 1)).map
          fun k => ((t.nvRange (idx a)).1 + k, (t.nvRange (idx a)).1 + k + 1) := by
  obtain ⟨u, v, xs, -, hnl, h2, -⟩ := h.nvList_S' ha hS
  have hlen : (items.nvList g a).length = xs.length + 2 := by rw [hnl]; simp
  have hnE := h.nEdges_S ha hS
  have hnv := (h.node a ha).nv_range
  have hne := (h.node a ha).ne_range
  apply skeleton_eq_of
  · rw [hne]; simp; omega
  · intro k hk
    simp only [List.length_cons, List.length_map, List.length_range] at hk
    rw [hl.edge_nvs k (by omega), nodeLayout, hS, hnv, hne,
      LayoutFacts.s_edge _ _ _ _ _ _ (by omega) (by omega) k (by omega)]
    cases k with
    | zero =>
      simp only [↓reduceIte, List.getElem_cons_zero, Prod.mk.injEq, true_and]; omega
    | succ k =>
      simp only [Nat.succ_ne_zero, ↓reduceIte, List.getElem_cons_succ, List.getElem_map,
        List.getElem_range, Prod.mk.injEq]
      omega

theorem skeleton_P {a : ItemId} (ha : a < items.size) (hP : items.type a = .P) {pos : Nat → Nat}
    (hl : RelabelLayout g items t idx a pos) :
    t.skeleton (idx a) =
      List.replicate (items.nEdges g a) ((t.nvRange (idx a)).1, (t.nvRange (idx a)).1 + 1) := by
  obtain ⟨u, v, -, hl2⟩ := h.nvList_P ha hP
  have hne := (h.node a ha).ne_range
  apply skeleton_eq_of
  · rw [hne, List.length_replicate]; omega
  · intro k hk
    rw [List.length_replicate] at hk
    rw [List.getElem_replicate, hl.edge_nvs k hk, nodeLayout, hP, (h.node a ha).nv_range, hne, hl2]
    exact LayoutFacts.p_edge _ _ _ _ _ k (by omega)

theorem nEdges_R {a : ItemId} (ha : a < items.size) (hR : items.type a = .R) (pos : Nat → Nat) :
    items.nEdges g a =
      (items.edgeChildren g pos (items.ordered g a (t.nvRange (idx a)).1 pos)).length + 1 := by
  rw [h.nEdges_node ha (by rw [hR]; rfl), Items.edgeChildren, List.length_map,
    (ordered_filter_perm a _ pos _).length_eq, h.virt_length ha, List.countP_eq_length_filter]
  simp [Items.capCount, Items.hasCap, hR, NodeType.isNode]

theorem skeleton_R {a : ItemId} (ha : a < items.size) (hR : items.type a = .R) {pos : Nat → Nat}
    (hl : RelabelLayout g items t idx a pos) :
    t.skeleton (idx a) =
      ((t.nvRange (idx a)).1, (t.nvRange (idx a)).1 + ((items.nvList g a).length - 1)) ::
        items.edgeChildren g pos (items.ordered g a (t.nvRange (idx a)).1 pos) := by
  obtain ⟨-, -, -, hlen4⟩ := h.nvList_R ha hR
  have hnE := h.nEdges_R ha hR pos
  have hnv := (h.node a ha).nv_range
  have hne := (h.node a ha).ne_range
  apply skeleton_eq_of
  · rw [hne]; simp; omega
  · intro k hk
    simp only [List.length_cons] at hk
    rw [hl.edge_nvs k (by omega), nodeLayout, hR, hnv, hne,
      LayoutFacts.r_edge _ _ _ _ _ _ (by omega) (by omega) k (by omega)]
    cases k with
    | zero =>
      simp only [↓reduceIte, List.getElem_cons_zero, Prod.mk.injEq, true_and]; omega
    | succ k =>
      simp only [Nat.succ_ne_zero, ↓reduceIte, Nat.add_sub_cancel, List.getElem_cons_succ]
      exact getElem!_pos _ _ (by omega)

/-! ### R nodes of an st-numbered item tree -/

/-- `pos` on an R node's vertex list is `nvSt + idxOf`. -/
theorem pos_eq {a : ItemId} (_ha : a < items.size) (hR : items.type a = .R)
    (hnd : (items.nvList g a).Nodup) {pos : Nat → Nat} (hl : RelabelLayout g items t idx a pos)
    {x : Nat} (hx : x ∈ items.nvList g a) :
    pos x = (t.nvRange (idx a)).1 + (items.nvList g a).idxOf x := by
  obtain ⟨hle, hp⟩ := hl.pos_ok hR x hx
  obtain ⟨hlt, hget⟩ := List.getElem?_eq_some_iff.1 hp
  have := hnd.idxOf_getElem _ hlt
  rw [hget] at this
  omega

theorem r_facts (hst : items.StNumbered) {a : ItemId} (ha : a < items.size)
    (hR : items.type a = .R) {pos : Nat → Nat} (hl : RelabelLayout g items t idx a pos) :
    (∀ q ∈ items.edgeChildren g pos (items.ordered g a (t.nvRange (idx a)).1 pos),
      (t.nvRange (idx a)).1 ≤ q.1 ∧ q.1 < q.2 ∧ q.2 < (t.nvRange (idx a)).2) ∧
    (items.edgeChildren g pos (items.ordered g a (t.nvRange (idx a)).1 pos)).Nodup ∧
    (items.edgeChildren g pos (items.ordered g a (t.nvRange (idx a)).1 pos)).Pairwise
      fun p q => p.1 + p.2 ≤ q.1 + q.2 := by
  obtain ⟨u, v, -, ⟨hnd, hes, -⟩, hor⟩ := hst a ha (Or.inr (Or.inr hR))
  rw [h.vertList_eq ha] at hnd hes hor
  have hnv := (h.node a ha).nv_range
  have hperm : ((items.ordered g a (t.nvRange (idx a)).1 pos).filter (· ≥ 1 + g.nv)).Perm
      ((items.ch a).filter fun c => items.type c ≠ .V) := by
    rw [← h.filter_ge_eq ha]; exact ordered_filter_perm a _ pos _
  have hvmem : ∀ c ∈ (items.ordered g a (t.nvRange (idx a)).1 pos).filter (· ≥ 1 + g.nv),
      ((items.vs c).1.getD 0, (items.vs c).2.getD 0) ∈ items.virtualEdges a := by
    intro c hc
    rw [Items.virtualEdges]
    exact List.mem_map.2 ⟨c, hperm.subset hc, rfl⟩
  have hpos : ∀ x ∈ items.nvList g a,
      pos x = (t.nvRange (idx a)).1 + (items.nvList g a).idxOf x :=
    fun x hx => h.pos_eq ha hR hnd hl hx
  refine ⟨?_, ?_, ?_⟩
  · intro q hq
    obtain ⟨c, hc, rfl⟩ := List.mem_map.1 hq
    have hq' := hvmem c hc
    obtain ⟨h1, h2, -⟩ := hes _ (List.mem_cons_of_mem _ hq')
    have ho := hor _ hq'
    simp only at h1 h2 ho ⊢
    rw [hpos _ h1, hpos _ h2, hnv]
    have := List.idxOf_lt_length_iff.2 h2
    omega
  · have hvnd : (items.virtualEdges a).Nodup :=
      List.Nodup.of_map _ (h.shapes.r_shape a ha hR).2.2.1
    have hE : items.edgeChildren g pos (items.ordered g a (t.nvRange (idx a)).1 pos) =
        (((items.ordered g a (t.nvRange (idx a)).1 pos).filter (· ≥ 1 + g.nv)).map
          fun c => ((items.vs c).1.getD 0, (items.vs c).2.getD 0)).map
            fun q => (pos q.1, pos q.2) := by
      rw [Items.edgeChildren, List.map_map]; rfl
    have hmap : ∀ q, q ∈ ((items.ordered g a (t.nvRange (idx a)).1 pos).filter (· ≥ 1 + g.nv)).map
        (fun c => ((items.vs c).1.getD 0, (items.vs c).2.getD 0)) → q ∈ items.virtualEdges a := by
      intro q hq
      rw [Items.virtualEdges]
      exact (hperm.map _).subset hq
    rw [hE]
    refine List.Nodup.map_on ?_ ((hperm.map _).nodup_iff.2 (by rwa [Items.virtualEdges] at hvnd))
    intro x hx y hy hxy
    obtain ⟨hx1, hx2, -⟩ := hes _ (List.mem_cons_of_mem _ (hmap x hx))
    obtain ⟨hy1, hy2, -⟩ := hes _ (List.mem_cons_of_mem _ (hmap y hy))
    simp only [Prod.mk.injEq] at hxy
    rw [hpos _ hx1, hpos _ hy1, hpos _ hx2, hpos _ hy2] at hxy
    exact Prod.ext ((List.idxOf_inj hx1).1 (by omega)) ((List.idxOf_inj hx2).1 (by omega))
  · have hord := Items.ordered_pairwise_loc g items a (t.nvRange (idx a)).1 pos hR
    rw [Items.edgeChildren, List.pairwise_map]
    refine (hord.filter (· ≥ 1 + g.nv)).imp_of_mem ?_
    intro c c' hc hc' hle
    have hcge : ¬ c < 1 + g.nv := by have := (List.mem_filter.1 hc).2; simp at this; omega
    have hcge' : ¬ c' < 1 + g.nv := by have := (List.mem_filter.1 hc').2; simp at this; omega
    simp only [Items.loc, LayoutShape.ite_of_neg hcge, LayoutShape.ite_of_neg hcge'] at hle
    obtain ⟨h1, h2, -⟩ := hes _ (List.mem_cons_of_mem _ (hvmem c hc))
    obtain ⟨h1', h2', -⟩ := hes _ (List.mem_cons_of_mem _ (hvmem c' hc'))
    dsimp only at h1 h2 h1' h2'
    have := (hl.pos_ok hR _ h1).1
    have := (hl.pos_ok hR _ h2).1
    have := (hl.pos_ok hR _ h1').1
    have := (hl.pos_ok hR _ h2').1
    dsimp only at hle ⊢
    omega

/-! ### `StOrder.st` -/

theorem st_S {a : ItemId} (ha : a < items.size) (hS : items.type a = .S) {pos : Nat → Nat}
    (hl : RelabelLayout g items t idx a pos) : t.StNumbered (idx a) := by
  obtain ⟨u, v, xs, -, hnl, h2, -⟩ := h.nvList_S' ha hS
  have hlen : (items.nvList g a).length = xs.length + 2 := by rw [hnl]; simp
  have hnv := (h.node a ha).nv_range
  rw [SpqrTree.StNumbered, h.skeleton_S ha hS hl]
  refine ⟨?_, ?_⟩
  · intro p hp
    rcases List.mem_cons.1 hp with rfl | hp
    · simp only; omega
    · obtain ⟨k, hk, rfl⟩ := List.mem_map.1 hp
      simp only; omega
  · intro nv h1 h3
    rw [hnv] at h3
    refine ⟨⟨(nv - 1, nv), List.mem_cons_of_mem _ (List.mem_map.2
        ⟨nv - 1 - (t.nvRange (idx a)).1, List.mem_range.2 (by omega), ?_⟩), rfl⟩,
      ⟨(nv, nv + 1), List.mem_cons_of_mem _ (List.mem_map.2
        ⟨nv - (t.nvRange (idx a)).1, List.mem_range.2 (by omega), ?_⟩), rfl⟩⟩
    · simp only [Prod.mk.injEq]; omega
    · simp only [Prod.mk.injEq]; omega

theorem st_P {a : ItemId} (ha : a < items.size) (hP : items.type a = .P) {pos : Nat → Nat}
    (hl : RelabelLayout g items t idx a pos) : t.StNumbered (idx a) := by
  obtain ⟨u, v, -, hl2⟩ := h.nvList_P ha hP
  have hnv := (h.node a ha).nv_range
  rw [hl2] at hnv
  rw [SpqrTree.StNumbered, h.skeleton_P ha hP hl]
  refine ⟨?_, ?_⟩
  · intro p hp
    rw [List.eq_of_mem_replicate hp]; simp
  · intro nv h1 h3
    rw [hnv] at h3
    simp at h3
    omega

theorem st_R (hst : items.StNumbered) {a : ItemId} (ha : a < items.size)
    (hR : items.type a = .R) {pos : Nat → Nat} (hl : RelabelLayout g items t idx a pos) :
    t.StNumbered (idx a) := by
  obtain ⟨u, v, hvs, ⟨hnd, hes, hint⟩, hor⟩ := hst a ha (Or.inr (Or.inr hR))
  rw [h.vertList_eq ha] at hnd hes hint hor
  obtain ⟨hE, -, -⟩ := h.r_facts hst ha hR hl
  obtain ⟨-, -, -, hlen4⟩ := h.nvList_R ha hR
  have hnv := (h.node a ha).nv_range
  have hperm : ((items.ordered g a (t.nvRange (idx a)).1 pos).filter (· ≥ 1 + g.nv)).Perm
      ((items.ch a).filter fun c => items.type c ≠ .V) := by
    rw [← h.filter_ge_eq ha]; exact ordered_filter_perm a _ pos _
  have hpos : ∀ x ∈ items.nvList g a,
      pos x = (t.nvRange (idx a)).1 + (items.nvList g a).idxOf x :=
    fun x hx => h.pos_eq ha hR hnd hl hx
  have hu0 : (items.nvList g a)[0]'(by omega) = u :=
    (List.getElem_eq_iff _).2 (vs_head (g := g) (show (items.vs a).1 = some u by rw [hvs]))
  have hvl : (items.nvList g a)[(items.nvList g a).length - 1]'(by omega) = v :=
    (List.getElem_eq_iff _).2 (vs_last (g := g) (show (items.vs a).2 = some v by rw [hvs]))
  have hidx_u : (items.nvList g a).idxOf u = 0 := by
    have := hnd.idxOf_getElem 0 (by omega); rwa [hu0] at this
  have hidx_v : (items.nvList g a).idxOf v = (items.nvList g a).length - 1 := by
    have := hnd.idxOf_getElem ((items.nvList g a).length - 1) (by omega); rwa [hvl] at this
  have hor' : ∀ p ∈ (u, v) :: items.virtualEdges a,
      (items.nvList g a).idxOf p.1 < (items.nvList g a).idxOf p.2 := by
    intro p hp
    rcases List.mem_cons.1 hp with rfl | hp
    · simp only; rw [hidx_u, hidx_v]; omega
    · exact hor p hp
  -- the skeleton edge of an item-level edge
  have hφ : ∀ p ∈ (u, v) :: items.virtualEdges a,
      ((t.nvRange (idx a)).1 + (items.nvList g a).idxOf p.1,
        (t.nvRange (idx a)).1 + (items.nvList g a).idxOf p.2) ∈ t.skeleton (idx a) := by
    intro p hp
    rw [h.skeleton_R ha hR hl]
    rcases List.mem_cons.1 hp with rfl | hp
    · simp only [hidx_u, hidx_v, Nat.add_zero]
      exact List.mem_cons.2 (Or.inl rfl)
    · rw [Items.virtualEdges] at hp
      obtain ⟨c, hc, rfl⟩ := List.mem_map.1 hp
      have hcF : c ∈ (items.ordered g a (t.nvRange (idx a)).1 pos).filter (· ≥ 1 + g.nv) :=
        hperm.mem_iff.2 hc
      have hq' : ((items.vs c).1.getD 0, (items.vs c).2.getD 0) ∈ items.virtualEdges a := by
        rw [Items.virtualEdges]; exact List.mem_map.2 ⟨c, hc, rfl⟩
      obtain ⟨h1, h2, -⟩ := hes _ (List.mem_cons_of_mem _ hq')
      refine List.mem_cons_of_mem _ (List.mem_map.2 ⟨c, hcF, ?_⟩)
      show (pos _, pos _) = _
      rw [hpos _ h1, hpos _ h2]
  rw [SpqrTree.StNumbered]
  refine ⟨?_, ?_⟩
  · intro p hp
    rw [h.skeleton_R ha hR hl] at hp
    rcases List.mem_cons.1 hp with rfl | hp
    · simp only; omega
    · exact (hE p hp).2.1
  · intro nv hlo hhi
    rw [hnv] at hhi
    have hklt : nv - (t.nvRange (idx a)).1 < (items.nvList g a).length := by omega
    have hxm := List.getElem_mem hklt
    have hidx := hnd.idxOf_getElem _ hklt
    have hhead : (items.nvList g a).head? ≠ some (items.nvList g a)[nv - (t.nvRange (idx a)).1] := by
      rw [List.head?_eq_getElem?, List.getElem?_eq_getElem (by omega)]
      intro heq
      have := hnd.getElem_inj_iff.1 (Option.some.inj heq)
      omega
    have hlast : (items.nvList g a).getLast? ≠
        some (items.nvList g a)[nv - (t.nvRange (idx a)).1] := by
      rw [List.getLast?_eq_getElem?, List.getElem?_eq_getElem (by omega)]
      intro heq
      have := hnd.getElem_inj_iff.1 (Option.some.inj heq)
      omega
    obtain ⟨⟨p, hp, hp'⟩, ⟨q, hq, hq'⟩⟩ := hint _ hxm hhead hlast
    refine ⟨⟨_, hφ p hp, ?_⟩, ⟨_, hφ q hq, ?_⟩⟩
    · rcases hp' with ⟨h1, h2⟩ | ⟨h1, h2⟩
      · exact absurd (hor' p hp) (by rw [h1]; omega)
      · simp only; rw [h1, hidx]; omega
    · rcases hq' with ⟨h1, h2⟩ | ⟨h1, h2⟩
      · simp only; rw [h1, hidx]; omega
      · exact absurd (hor' q hq) (by rw [h1]; omega)

/-! ### `StOrder.dom` -/

theorem dom_small {a : ItemId} (ha : a < items.size) (hle : items.nEdges g a ≤ 1) :
    t.EdgeDominance (idx a) := by
  rw [SpqrTree.EdgeDominance, ownEdges_nvs]
  apply List.pairwise_of_length_le_one
  rw [List.length_drop, skeleton_length, (h.node a ha).ne_range]
  omega

theorem dom_node (hst : items.StNumbered) {a : ItemId} (ha : a < items.size)
    {pos : Nat → Nat} (hl : RelabelLayout g items t idx a pos) : t.EdgeDominance (idx a) := by
  cases hty : items.type a
  · exact h.dom_small ha (by rw [nEdges_not_node (g := g) (i := a) (by rw [hty]; rfl)]; omega)
  · exact h.dom_small ha (by rw [nEdges_not_node (g := g) (i := a) (by rw [hty]; rfl)]; omega)
  · obtain ⟨e, he, rfl⟩ := h.type_Q_eq ha hty
    have hne1 := h.nEdges_Q he
    exact h.dom_small ha (by omega)
  · exact h.dom_small ha (by rw [h.nEdges_IO ha (Or.inl hty)])
  · exact h.dom_small ha (by rw [h.nEdges_IO ha (Or.inr hty)])
  · rw [SpqrTree.EdgeDominance, ownEdges_nvs, h.skeleton_S ha hty hl, h.hasCap_eq ha,
      show items.hasCap a = true by simp [Items.hasCap, hty, NodeType.isNode],
      LayoutShape.ite_of_pos rfl, List.drop_succ_cons, List.drop_zero, List.pairwise_map]
    refine List.pairwise_lt_range.imp ?_
    intro i j hij _ hle
    simp only at hle
    omega
  · rw [SpqrTree.EdgeDominance, ownEdges_nvs, h.skeleton_P ha hty hl, List.drop_replicate,
      List.pairwise_replicate]
    exact Or.inr fun hne => absurd rfl hne
  · obtain ⟨-, -, hsum⟩ := h.r_facts hst ha hty hl
    rw [SpqrTree.EdgeDominance, ownEdges_nvs, h.skeleton_R ha hty hl, h.hasCap_eq ha,
      show items.hasCap a = true by simp [Items.hasCap, hty, NodeType.isNode],
      LayoutShape.ite_of_pos rfl, List.drop_succ_cons, List.drop_zero]
    exact pairwise_dominance_of_sorted_sum _ hsum

/-! ### `StOrder.adj`: row transport -/

theorem nvBounds_cover : ∀ n, n ≤ t.size → ∀ nv, nv < t.nvBounds[n]! →
    ∃ m, m < n ∧ t.nvBounds[m]! ≤ nv ∧ nv < t.nvBounds[m + 1]! := by
  intro n
  induction n with
  | zero => intro _ nv hnv; rw [h.ridx.nv_zero] at hnv; omega
  | succ n ih =>
    intro hn nv hnv
    by_cases hlt : nv < t.nvBounds[n]!
    · obtain ⟨m, hm, h1, h2⟩ := ih (by omega) nv hlt
      exact ⟨m, by omega, h1, h2⟩
    · exact ⟨n, by omega, by omega, hnv⟩

/-- Which node a node-vertex belongs to: the node whose range contains it. -/
theorem nodeOfNv_range {nv n : Nat} (hnv : nv < t.nodeVerts.size) (hn : t.nodeOfNv nv = some n) :
    ∃ a, a < items.size ∧ idx a = n ∧ (t.nvRange (idx a)).1 ≤ nv ∧ nv < (t.nvRange (idx a)).2 := by
  obtain ⟨m, hm, hPm, hnext⟩ := h.nvBounds_cover t.size le_rfl nv (by rwa [h.ridx.nv_last])
  obtain ⟨a, ha, hidx⟩ := h.idx_surj hm
  have hlo : (t.nvRange (idx a)).1 ≤ nv := by rw [nvRange_fst, hidx]; exact hPm
  have hhi : nv < (t.nvRange (idx a)).2 := by rw [nvRange_snd, hidx]; exact hnext
  refine ⟨a, ha, ?_, hlo, hhi⟩
  have hnvr := (h.node a ha).nv_range
  have hk := (h.node a ha).node_verts (nv - (t.nvRange (idx a)).1) (by omega)
  rw [Nat.add_sub_cancel' hlo] at hk
  have hk' := (getElem!_pos t.nodeVerts nv hnv).symm.trans hk
  rw [SpqrTree.nodeOfNv, Array.getElem?_eq_getElem hnv, Option.map_some, hk'] at hn
  have := Option.some.inj hn
  simp only at this
  exact this

/-- Row transport: a global adjacency row of a node-vertex of `a` is the corresponding row of
`layoutNode` (given `relabelTree_adj`'s row-start fact for `a`). -/
theorem adjRow_eq_layout {a : ItemId} (ha : a < items.size) {pos : Nat → Nat}
    (hl : RelabelLayout g items t idx a pos)
    (hL : (nodeLayout g items t idx a pos).Local (idx a) (t.nvRange (idx a)).1 (t.nvRange (idx a)).2
      (t.neRange (idx a)).1 (t.neRange (idx a)).2)
    (hA : t.adjBounds[2 * (t.nvRange (idx a)).1]! = 2 * (t.neRange (idx a)).1)
    (r : Nat) (hr1 : 2 * (t.nvRange (idx a)).1 ≤ r) (hr2 : r < 2 * (t.nvRange (idx a)).2) :
    t.adjRow r = row (t.nvRange (idx a)).1 (t.neRange (idx a)).1 (nodeLayout g items t idx a pos) r := by
  have hnv := (h.node a ha).nv_range
  have hne := (h.node a ha).ne_range
  have hb : ∀ r', 2 * (t.nvRange (idx a)).1 ≤ r' → r' ≤ 2 * (t.nvRange (idx a)).2 →
      t.adjBounds[r']! = rowBound (t.nvRange (idx a)).1 (t.neRange (idx a)).1
        (nodeLayout g items t idx a pos) r' := by
    intro r' h1 h2
    by_cases h0 : r' = 2 * (t.nvRange (idx a)).1
    · subst h0; rw [hA]; simp [rowBound]
    · rw [rowBound, LayoutShape.ite_of_neg h0]
      have := hl.adj_bounds (r' - 2 * (t.nvRange (idx a)).1) (by omega) (by omega)
      rw [Nat.add_sub_cancel' h1] at this
      exact this
  have hge := hL.rowBound_ge
  have hle := hL.rowBound_le
  rw [SpqrTree.adjRow, row, hb r hr1 (by omega), hb (r + 1) (by omega) (by omega)]
  apply List.map_congr_left
  intro k hk
  rw [List.mem_range] at hk
  have h1 := hge r hr1 (by omega)
  have h2 := hle (r + 1) (by omega) (by omega)
  have := hl.adj_dat (rowBound (t.nvRange (idx a)).1 (t.neRange (idx a)).1
    (nodeLayout g items t idx a pos) r + k - 2 * (t.neRange (idx a)).1) (by omega)
  rw [show 2 * (t.neRange (idx a)).1 + (rowBound (t.nvRange (idx a)).1 (t.neRange (idx a)).1
    (nodeLayout g items t idx a pos) r + k - 2 * (t.neRange (idx a)).1) =
    rowBound (t.nvRange (idx a)).1 (t.neRange (idx a)).1 (nodeLayout g items t idx a pos) r + k
    by omega] at this
  exact this

theorem adj_rows {a : ItemId} (ha : a < items.size) {pos : Nat → Nat}
    (hl : RelabelLayout g items t idx a pos)
    (hL : (nodeLayout g items t idx a pos).Local (idx a) (t.nvRange (idx a)).1 (t.nvRange (idx a)).2
      (t.neRange (idx a)).1 (t.neRange (idx a)).2)
    (hA : t.adjBounds[2 * (t.nvRange (idx a)).1]! = 2 * (t.neRange (idx a)).1)
    {nv : Nat} (hlo : (t.nvRange (idx a)).1 ≤ nv) (hhi : nv < (t.nvRange (idx a)).2) :
    t.AdjBracket nv ↔
      ((∀ x ∈ row (t.nvRange (idx a)).1 (t.neRange (idx a)).1 (nodeLayout g items t idx a pos)
          (2 * nv), x.destNv < nv) ∧
        (∀ x ∈ row (t.nvRange (idx a)).1 (t.neRange (idx a)).1 (nodeLayout g items t idx a pos)
          (2 * nv + 1), nv < x.destNv) ∧
        ((row (t.nvRange (idx a)).1 (t.neRange (idx a)).1 (nodeLayout g items t idx a pos)
          (2 * nv)).map (·.destNv)).Pairwise (· ≥ ·) ∧
        ((row (t.nvRange (idx a)).1 (t.neRange (idx a)).1 (nodeLayout g items t idx a pos)
          (2 * nv + 1)).map (·.destNv)).Pairwise (· ≥ ·)) := by
  rw [SpqrTree.AdjBracket, h.adjRow_eq_layout ha hl hL hA _ (by omega) (by omega),
    h.adjRow_eq_layout ha hl hL hA _ (by omega) (by omega)]

section PerType

variable {a : ItemId} (ha : a < items.size) {pos : Nat → Nat}
  (hl : RelabelLayout g items t idx a pos)
  (hA : t.adjBounds[2 * (t.nvRange (idx a)).1]! = 2 * (t.neRange (idx a)).1)
  {nv : Nat} (hlo : (t.nvRange (idx a)).1 ≤ nv) (hhi : nv < (t.nvRange (idx a)).2)
include ha hl hA hlo hhi

theorem adj_F (hty : items.type a = .F) : t.AdjBracket nv := by
  have hnv := (h.node a ha).nv_range
  have hne := (h.node a ha).ne_range
  have h0 := nEdges_not_node (g := g) (i := a) (by rw [hty]; rfl)
  have hL : (nodeLayout g items t idx a pos).Local (idx a) (t.nvRange (idx a)).1
      (t.nvRange (idx a)).2 (t.neRange (idx a)).1 (t.neRange (idx a)).2 := by
    rw [nodeLayout, hty]; exact LayoutShape.local_F _ _ _ _ _ _ (by omega) (by omega)
  rw [h.adj_rows ha hl hL hA hlo hhi, nodeLayout, hty, LayoutShape.layoutNode_F_eq,
    LayoutShape.runF_row _ _ _ _ _ (by omega) (by omega),
    LayoutShape.runF_row _ _ _ _ _ (by omega) (by omega)]
  simp

theorem adj_QI (hty : items.type a = .Q ∨ items.type a = .I)
    (hlen : (items.nvList g a).length = 2) (hne1 : items.nEdges g a = 1) : t.AdjBracket nv := by
  have hnv := (h.node a ha).nv_range
  have hne := (h.node a ha).ne_range
  have hL : (nodeLayout g items t idx a pos).Local (idx a) (t.nvRange (idx a)).1
      (t.nvRange (idx a)).2 (t.neRange (idx a)).1 (t.neRange (idx a)).2 := by
    rw [nodeLayout]; exact LayoutShape.local_QI _ _ _ _ _ _ _ hty (by omega) (by omega)
  rw [h.adj_rows ha hl hL hA hlo hhi, nodeLayout,
    LayoutShape.layoutNode_QI_eq _ _ _ _ _ _ _ hty (by omega),
    LayoutShape.runQI_row _ _ _ _ _ (by omega) (by omega) _ (by omega) (by omega),
    LayoutShape.runQI_row _ _ _ _ _ (by omega) (by omega) _ (by omega) (by omega)]
  split_ifs <;> simp <;> omega

theorem adj_S (hty : items.type a = .S) : t.AdjBracket nv := by
  have hnv := (h.node a ha).nv_range
  have hne := (h.node a ha).ne_range
  have hnE := h.nEdges_S ha hty
  obtain ⟨u, v, xs, -, hnl, h2, -⟩ := h.nvList_S' ha hty
  have hlen : (items.nvList g a).length = xs.length + 2 := by rw [hnl]; simp
  have hL : (nodeLayout g items t idx a pos).Local (idx a) (t.nvRange (idx a)).1
      (t.nvRange (idx a)).2 (t.neRange (idx a)).1 (t.neRange (idx a)).2 := by
    rw [nodeLayout, hty]; exact LayoutShape.local_S _ _ _ _ _ _ (by omega) (by omega)
  rw [h.adj_rows ha hl hL hA hlo hhi, nodeLayout, hty,
    LayoutShape.layoutNode_S_eq _ _ _ _ _ _ (by omega),
    LayoutShape.runS_row _ _ _ _ _ (by omega) (by omega) _ (by omega) (by omega),
    LayoutShape.runS_row _ _ _ _ _ (by omega) (by omega) _ (by omega) (by omega)]
  split_ifs <;> simp <;> omega

theorem adj_P (hty : items.type a = .P) : t.AdjBracket nv := by
  have hnv := (h.node a ha).nv_range
  have hne := (h.node a ha).ne_range
  obtain ⟨u, v, -, hl2⟩ := h.nvList_P ha hty
  have hlen : (items.nvList g a).length = 2 := by rw [hl2]; rfl
  have hL : (nodeLayout g items t idx a pos).Local (idx a) (t.nvRange (idx a)).1
      (t.nvRange (idx a)).2 (t.neRange (idx a)).1 (t.neRange (idx a)).2 := by
    rw [nodeLayout, hty]; exact LayoutShape.local_P _ _ _ _ _ _ (by omega) (by omega)
  have hv2 : (t.nvRange (idx a)).2 - (t.nvRange (idx a)).1 = 2 := by omega
  have hle : (t.neRange (idx a)).1 ≤ (t.neRange (idx a)).2 := by omega
  have r0 := LayoutShape.runP_row (idx a) (t.nvRange (idx a)).1 (t.nvRange (idx a)).2
    (t.neRange (idx a)).1 (t.neRange (idx a)).2 hv2 hle (2 * nv) (by omega) (by omega)
  have r1 := LayoutShape.runP_row (idx a) (t.nvRange (idx a)).1 (t.nvRange (idx a)).2
    (t.neRange (idx a)).1 (t.neRange (idx a)).2 hv2 hle (2 * nv + 1) (by omega) (by omega)
  have hPeq := LayoutShape.layoutNode_P_eq (idx a) (t.nvRange (idx a)).1 (t.nvRange (idx a)).2
    (t.neRange (idx a)).1 (t.neRange (idx a)).2
    (items.edgeChildren g pos (items.ordered g a (t.nvRange (idx a)).1 pos)) (by omega)
  rw [h.adj_rows ha hl hL hA hlo hhi, nodeLayout, hty, hPeq, r0, r1]
  split_ifs <;> simp [List.pairwise_map] <;> omega

theorem adj_R (hst : items.StNumbered) (hty : items.type a = .R) : t.AdjBracket nv := by
  have hnv := (h.node a ha).nv_range
  have hne := (h.node a ha).ne_range
  obtain ⟨hE, hnd, hsum⟩ := h.r_facts hst ha hty hl
  have hnE := h.nEdges_R ha hty pos
  obtain ⟨-, -, -, hlen4⟩ := h.nvList_R ha hty
  have hL : (nodeLayout g items t idx a pos).Local (idx a) (t.nvRange (idx a)).1
      (t.nvRange (idx a)).2 (t.neRange (idx a)).1 (t.neRange (idx a)).2 := by
    rw [nodeLayout, hty]; exact LayoutShape.local_R _ _ _ _ _ _ (by omega) hE (by omega)
  have hb := LayoutR.layoutNode_R_bracket (idx a) (t.nvRange (idx a)).1 (t.nvRange (idx a)).2
    (t.neRange (idx a)).1 (t.neRange (idx a)).2
    (items.edgeChildren g pos (items.ordered g a (t.nvRange (idx a)).1 pos)) (by omega) hE
    (pairwise_dominance_of_sorted_sum _ hsum) hnd (by omega) nv hlo hhi
  rw [h.adj_rows ha hl hL hA hlo hhi, nodeLayout, hty]
  exact hb

end PerType

theorem adj_node (hst : items.StNumbered) {a : ItemId} (ha : a < items.size)
    {pos : Nat → Nat} (hl : RelabelLayout g items t idx a pos)
    (hA : t.adjBounds[2 * (t.nvRange (idx a)).1]! = 2 * (t.neRange (idx a)).1)
    {nv : Nat} (h1 : 1 < t.nVerts (idx a)) (hlo : (t.nvRange (idx a)).1 ≤ nv)
    (hhi : nv < (t.nvRange (idx a)).2) : t.AdjBracket nv := by
  have hnv := (h.node a ha).nv_range
  have hV : t.nVerts (idx a) = (items.nvList g a).length := by rw [SpqrTree.nVerts, hnv]; omega
  rw [hV] at h1
  cases hty : items.type a
  · exact h.adj_F ha hl hA hlo hhi hty
  · exfalso; rw [h.nvList_V ha hty] at h1; simp at h1
  · obtain ⟨e, he, rfl⟩ := h.type_Q_eq ha hty
    refine h.adj_QI ha hl hA hlo hhi (Or.inl hty) ?_ (h.nEdges_Q he)
    have := h.nvList_Q_le he
    omega
  · obtain ⟨u, v, -, hl2⟩ := h.nvList_I ha hty
    exact h.adj_QI ha hl hA hlo hhi (Or.inr hty) (by rw [hl2]; rfl) (h.nEdges_IO ha (Or.inl hty))
  · obtain ⟨v, -, hl1⟩ := h.nvList_O ha hty
    exfalso; rw [hl1] at h1; simp at h1
  · exact h.adj_S ha hl hA hlo hhi hty
  · exact h.adj_P ha hl hA hlo hhi hty
  · exact h.adj_R ha hl hA hlo hhi hst hty

/-! ### Assembly -/

theorem stOrder (hst : items.StNumbered)
    (hA : ∀ n, n < t.size → t.adjBounds[2 * (t.nvRange n).1]! = 2 * (t.neRange n).1) :
    t.StOrder where
  st := by
    intro n hn hty
    obtain ⟨a, ha, rfl⟩ := h.idx_surj hn
    rw [h.type_eq ha] at hty
    obtain ⟨pos, hl⟩ := (h.node a ha).layout
    rcases hty with hS | hP | hR
    · exact h.st_S ha hS hl
    · exact h.st_P ha hP hl
    · exact h.st_R hst ha hR hl
  dom := by
    intro n hn
    obtain ⟨a, ha, rfl⟩ := h.idx_surj hn
    obtain ⟨pos, hl⟩ := (h.node a ha).layout
    exact h.dom_node hst ha hl
  adj := by
    intro nv n hnv hn h1
    obtain ⟨a, ha, rfl, hlo, hhi⟩ := h.nodeOfNv_range hnv hn
    obtain ⟨pos, hl⟩ := (h.node a ha).layout
    exact h.adj_node hst ha hl (hA _ (h.idx_lt ha)) h1 hlo hhi

end RelabelOK

/-- Phase 3: relabelling well-formed items in s-t order gives st-ordered output. From the per-node
interface, modulo `relabel_node_spec` (via `relabelTree_adj` for the CSR row starts). -/
theorem relabel_st (g : Graph) (items : Items) (hst : items.StNumbered) (h : items.WF g) :
    (relabelTree g items).StOrder := by
  obtain ⟨idx, hok⟩ := relabelOK_of_wf g items h
  exact hok.stOrder hst (relabelTree_adj g items h (hok.rOriented hst)).1

end Spqr
