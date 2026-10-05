import Spqr.EarFrame
import Spqr.EarRoot
import Spqr.EarShape

/-!
# The forest walk: `walk_sides` and its consumers

`walk_sides` (moved here from `WalkCover.lean`, below `EarRoot`/`EarShape` in the import order) is
`walk_sides_of_roots` applied to `RootsBook forest (init g tern)`, which is threaded root by root
(`RootState`): `BookTree` at a root is the admitted `walkTree_book`, `GuardsTree` is
`walkTree_guards'`, `Inv' 0` at a later root follows from `Inv' 0` after the previous root walk
(`invTree`) through the root pop/append (`rootItem` is parentless by `Place.root`, the popped stack is
empty by `RootOK`), and the freshness of a later root's `V`/`Q` items from `Place` (parentless) and the
frame `walkTree_frame` (`EarFrame.lean`, childless; its typing hypothesis is `Place.types`).
-/

namespace Spqr
open WalkM

namespace WalkState

theorem pushed_append {g : Graph} {P : ItemId → Prop} {vs es vs' es' : List Nat} {i : ItemId}
    (h : Pushed g (Pushed g P vs es) vs' es' i) : Pushed g P (vs ++ vs') (es ++ es') i := by
  rcases h with (h | ⟨v, hv, rfl⟩ | ⟨e, he, rfl⟩) | ⟨v, hv, rfl⟩ | ⟨e, he, rfl⟩
  · exact Or.inl h
  · exact Or.inr (Or.inl ⟨v, List.mem_append_left _ hv, rfl⟩)
  · exact Or.inr (Or.inr ⟨e, List.mem_append_left _ he, rfl⟩)
  · exact Or.inr (Or.inl ⟨v, List.mem_append_right _ hv, rfl⟩)
  · exact Or.inr (Or.inr ⟨e, List.mem_append_right _ he, rfl⟩)

/-- The endpoints of an edge of a well-formed tree are vertices of the tree. -/
theorem tree_inc_verts {g : Graph} {t : DfsTree} (hwf : t.WF []) (hends : t.Ends g) {e x : Nat}
    (he : e ∈ t.edges) (hx : g.Inc e x) : x ∈ t.verts := by
  obtain ⟨v, outs⟩ := t
  rw [DfsTree.WF] at hwf
  rw [DfsTree.Ends] at hends
  obtain ⟨o, ho, hs⟩ := mem_subEdges_edgesList.1 he
  rcases endsOut_wf g [] v o (hwf.2 o ho) (hends o ho) e hs x hx with h | rfl | h
  · exact absurd h (List.not_mem_nil)
  · exact List.mem_cons_self ..
  · exact List.mem_cons_of_mem _ (mem_vertsList_of_verts ho h)

/-- Per-tree edge completeness from the forest: every edge of `g` incident to a vertex of a tree
of the forest is an edge of that tree (the edge lies in some tree, whose vertices it joins, and the
trees' vertex sets are disjoint). `walkTree_book`/`walkTree_ear` need it (`EarFalse.lean`). -/
theorem comp_of_forest {g : Graph} {forest : List DfsTree} (hf : ForestOK g forest)
    (hwf : ∀ t ∈ forest, t.WF []) (hends : ∀ t ∈ forest, t.Ends g)
    (hcov : ∀ e, e < g.ne → e ∈ forest.flatMap DfsTree.edges) {t : DfsTree} (ht : t ∈ forest) :
    ∀ e, e < g.ne → ∀ x, g.Inc e x → x ∈ t.verts → e ∈ t.edges := by
  intro e he x hx hxt
  obtain ⟨t', ht', he'⟩ := List.mem_flatMap.1 (hcov e he)
  have hxt' := tree_inc_verts (hwf t' ht') (hends t' ht') he' hx
  obtain ⟨l₁, l₂, rfl⟩ := List.append_of_mem ht
  have hnd := hf.verts_nodup
  simp only [List.flatMap_append, List.flatMap_cons] at hnd
  rcases List.mem_append.1 ht' with h | h
  · exact absurd (List.mem_flatMap.2 ⟨t', h, hxt'⟩) fun hm =>
      List.disjoint_of_nodup_append hnd hm (List.mem_append_left _ hxt)
  · rcases List.mem_cons.1 h with rfl | h
    · exact he'
    · exact absurd hxt fun hm => List.disjoint_of_nodup_append (List.nodup_append.1 hnd).2.1 hm
        (List.mem_flatMap.2 ⟨t', h, hxt'⟩)

/-- The state at the start of a root walk, after the roots `pre`. -/
structure RootState (g : Graph) (pre : List DfsTree) (s : WalkState) : Prop where
  place : s.Place g (Pushed g (fun _ => False) (pre.flatMap DfsTree.verts) (pre.flatMap DfsTree.edges))
    (fun _ => False)
  g_eq : s.g = g
  sv : s.stackVerts.size = g.nv
  sd : s.stackDir.size = g.nv
  fo : s.firstOccurrence.size = g.nv
  shape : Shape s
  tstack : s.tstack = []
  inv : s.Inv' 0
  fresh : ∀ i, 0 < i → i < 1 + g.nv + g.ne → (∀ v ∈ pre.flatMap DfsTree.verts, i ≠ vertItem v) →
    (∀ e ∈ pre.flatMap DfsTree.edges, i ≠ edgeItem g e) → Items.ch s.items i = []

theorem rootState_init (g : Graph) (ternarize : Bool) : RootState g [] (WalkState.init g ternarize) where
  place := (WalkState.init_place g ternarize).mono (fun _ h => Or.inl h) fun _ h => h
  g_eq := rfl
  sv := by simp [WalkState.init]
  sd := by simp [WalkState.init]
  fo := by simp [WalkState.init]
  shape := init_shape g ternarize
  tstack := rfl
  inv := init_inv g ternarize
  fresh i _ _ _ _ := Items.initialItems_ch g i

namespace RootState

variable {g : Graph} {pre rest : List DfsTree} {t : DfsTree} {s : WalkState}

theorem hvlt (hf : ForestOK g (pre ++ t :: rest)) : ∀ v ∈ t.verts, v < g.nv := fun v hv =>
  hf.verts_lt v (by simp [hv])
theorem helt (hf : ForestOK g (pre ++ t :: rest)) : ∀ e ∈ t.edges, e < g.ne := fun e he =>
  hf.edges_lt e (by simp [he])
theorem hvn (hf : ForestOK g (pre ++ t :: rest)) : t.verts.Nodup ∧ ∀ v ∈ t.verts, v ∉ pre.flatMap DfsTree.verts := by
  have h := hf.verts_nodup
  rw [List.flatMap_append, List.flatMap_cons, List.nodup_append, List.nodup_append] at h
  exact ⟨h.2.1.1, fun v hv hv' => h.2.2 v hv' v (List.mem_append_left _ hv) rfl⟩
theorem hen (hf : ForestOK g (pre ++ t :: rest)) : t.edges.Nodup ∧ ∀ e ∈ t.edges, e ∉ pre.flatMap DfsTree.edges := by
  have h := hf.edges_nodup
  rw [List.flatMap_append, List.flatMap_cons, List.nodup_append, List.nodup_append] at h
  exact ⟨h.2.1.1, fun e he he' => h.2.2 e he' e (List.mem_append_left _ he) rfl⟩

theorem hPv (hf : ForestOK g (pre ++ t :: rest)) :
    ∀ v ∈ t.verts, ¬ Pushed g (fun _ => False) (pre.flatMap DfsTree.verts) (pre.flatMap DfsTree.edges) (vertItem v) := by
  rintro v hv (h | ⟨w, hw, hvw⟩ | ⟨e, _, hve⟩)
  · exact h
  · exact (hvn hf).2 v hv (vertItem_inj hvw ▸ hw)
  · exact vertItem_ne_edgeItem (hvlt hf v hv) e hve
theorem hPe (hf : ForestOK g (pre ++ t :: rest)) :
    ∀ e ∈ t.edges, ¬ Pushed g (fun _ => False) (pre.flatMap DfsTree.verts) (pre.flatMap DfsTree.edges) (edgeItem g e) := by
  rintro e he (h | ⟨w, hw, hew⟩ | ⟨e', he', hee⟩)
  · exact h
  · exact vertItem_ne_edgeItem (hf.verts_lt w (by simp [hw])) e hew.symm
  · exact (hen hf).2 e he (edgeItem_inj hee ▸ he')

/-- `BookTree` at a root from the threaded state (the admitted `walkTree_book`). -/
theorem book (h : RootState g pre s) (hf : ForestOK g (pre ++ t :: rest)) (hwf : t.WF [])
    (hends : t.Ends g) (hcomp : ∀ e, e < g.ne → ∀ x, g.Inc e x → x ∈ t.verts → e ∈ t.edges) :
    BookTree t 0 s := by
  obtain ⟨hp, hg, hsv, hsd, hfo, hs, hts, hi, hfresh⟩ := h
  subst hg
  refine walkTree_book t s hwf hends (hvlt hf) (helt hf) (hvn hf).1 (hen hf).1 hcomp hsv hsd hfo hts hi hs
    (fun v hv => ⟨?_, ?_⟩) fun e he => ⟨?_, ?_⟩
  · refine hfresh _ (by show 0 < 1 + v; omega) (by have := hvlt hf v hv; show 1 + v < _; omega) ?_ ?_
    · intro w hw hvw; exact (hvn hf).2 v hv (vertItem_inj hvw ▸ hw)
    · intro e _ hve; exact vertItem_ne_edgeItem (hvlt hf v hv) e hve
  · exact noParent_of_cnt_eq_zero (hp.cnt_eq_zero (by show 0 < 1 + v; omega)
      (by have := hvlt hf v hv; show 1 + v < _; omega) (hPv hf v hv))
  · refine hfresh _ (by show 0 < 1 + s.g.nv + e; omega) (by have := helt hf e he; show 1 + s.g.nv + e < _; omega) ?_ ?_
    · intro w hw hew; exact vertItem_ne_edgeItem (hf.verts_lt w (by simp [hw])) e hew.symm
    · intro e' he' hee; exact (hen hf).2 e he (edgeItem_inj hee ▸ he')
  · exact noParent_of_cnt_eq_zero (hp.cnt_eq_zero (by show 0 < 1 + s.g.nv + e; omega)
      (by have := helt hf e he; show 1 + s.g.nv + e < _; omega) (hPe hf e he))

/-- One root: walk it, then pop its entry onto `rootItem`. -/
theorem step (h : RootState g pre s) (hf : ForestOK g (pre ++ t :: rest)) (hwf : t.WF [])
    (hends : t.Ends g) (hcomp : ∀ e, e < g.ne → ∀ x, g.Inc e x → x ∈ t.verts → e ∈ t.edges) :
    wp (walkTree t 0) (fun _ s₁ =>
      wp (popTstack >>= fun top => modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 })
        (fun _ s₂ => RootState g (pre ++ [t]) s₂) s₁) s := by
  have hb := h.book hf hwf hends hcomp
  have hg := gbTree t 0 s hb
  have hi' : ∀ v outs, t = .node v outs →
      ({ s with stackVerts := s.stackVerts.set! 0 v } : WalkState).Inv' 0 :=
    fun _ _ _ => h.inv.stackVerts_of_nil h.tstack _
  have hnv : 0 < g.nv := by
    obtain ⟨v, outs⟩ := t
    exact Nat.lt_of_le_of_lt (Nat.zero_le _) (hvlt hf v (by simp [DfsTree.verts]))
  have hplace := (walk_place_aux g).1 t 0 _ _ s h.place (hvlt hf) (helt hf) (hvn hf).1 (hen hf).1
    (hPv hf) (hPe hf)
  have hframe := walkTree_frame t 0 s h.place.types
  have hinv := invTree t 0 s hi' h.shape hg hb
  have hrk : wp (walkTree t 0) (fun _ s' => RootOK s') s :=
    walkTree_rootOK t s hi' h.shape hg hb h.tstack (by rw [h.sd]; exact hnv)
  refine wp_mono _ (wp_and hplace (wp_and hframe (wp_and hinv hrk)))
    fun _ s₁ ⟨hp₁, ⟨hg₁, hsv₁, hsd₁, hfo₁, hch₁⟩, ⟨hi₁, hs₁⟩, hrk₁⟩ => ?_
  obtain ⟨tt, hts₁, -, -⟩ := hrk₁
  simp only [wp_bind, wp_popTstack, wp_modifyItem]
  have hroot : ∀ p, ¬ Items.IsParent s₁.items p rootItem := noParent_of_cnt_eq_zero hp₁.root
  have hpop : ({ s₁ with tstack := s₁.tstack.tail } : WalkState).Inv' 0 :=
    ⟨fun _ _ _ h' => by simp [hts₁] at h', fun i h1 h2 => ⟨(hi₁.nodes i h1 h2).conn, (hi₁.nodes i h1 h2).attached⟩⟩
  refine
    { place := ?_
      g_eq := hg₁.trans h.g_eq
      sv := hsv₁.trans h.sv
      sd := hsd₁.trans h.sd
      fo := hfo₁.trans h.fo
      shape := ?_
      tstack := by simp [hts₁]
      inv := ?_
      fresh := ?_ }
  · refine (hp₁.root_append).mono (fun i hi => ?_) fun _ hi => hi
    simpa using pushed_append hi
  · refine (hs₁.tstack (l := s₁.tstack.tail) fun e he => hs₁.span e (List.mem_of_mem_tail he)).modify
      rootItem _ (fun _ => rfl) fun hj c hc => ?_
    simp only [List.mem_append] at hc
    rcases hc with hc | hc
    · exact hs₁.ch_lt rootItem c (by simpa [Items.ch, Items.IsParent, hj] using hc)
    · exact hs₁.span tt (by simp [hts₁]) c (List.mem_append_right _ (by simpa [hts₁] using hc))
  · exact hpop.modifyCh rootItem _ (by show (0 : Nat) < 1 + _ + _; omega) hroot (by simp [hts₁])
  · intro i hi0 hi hv he
    simp only [List.flatMap_append, List.flatMap_cons, List.flatMap_nil, List.append_nil] at hv he
    show Items.ch (s₁.items.modify rootItem _) i = []
    rw [Items.ch_modify_of_ne rootItem _ (Nat.ne_of_gt hi0)]
    rw [hch₁ i hi0 (by rw [h.g_eq]; exact hi) (fun v hv' => hv v (List.mem_append_right _ hv'))
      (fun e he' => by rw [h.g_eq]; exact he e (List.mem_append_right _ he'))]
    exact h.fresh i hi0 hi (fun v hv' => hv v (List.mem_append_left _ hv'))
      fun e he' => he e (List.mem_append_left _ he')

end RootState

theorem rootsBook_of_state {g : Graph} : ∀ (forest pre : List DfsTree) (s : WalkState),
    RootState g pre s → ForestOK g (pre ++ forest) → (∀ t ∈ forest, t.WF []) →
    (∀ t ∈ forest, t.Ends g) →
    (∀ t ∈ forest, ∀ e, e < g.ne → ∀ x, g.Inc e x → x ∈ t.verts → e ∈ t.edges) →
    RootsBook forest s
  | [], _, _, _, _, _, _, _ => trivial
  | t :: rest, pre, s, h, hf, hwf, hends, hcomp => by
    have hb := h.book hf (hwf t (by simp)) (hends t (by simp)) (hcomp t (by simp))
    refine ⟨gbTree t 0 s hb, hb, h.inv, ?_⟩
    refine wp_mono _ (h.step hf (hwf t (by simp)) (hends t (by simp)) (hcomp t (by simp)))
      fun _ s₁ h₁ => ?_
    refine wp_mono _ h₁ fun _ s₂ h₂ => ?_
    exact rootsBook_of_state rest (pre ++ [t]) s₂ h₂ (by simpa using hf)
      (fun t' ht' => hwf t' (by simp [ht'])) (fun t' ht' => hends t' (by simp [ht']))
      fun t' ht' => hcomp t' (by simp [ht'])

end WalkState
open WalkState

/-! ### Admitted: the side discipline of the walk

The one orientation fact this file needs: at every discard site of the walk the discarded side is
empty (and the entries the sites pop exist). PROOF.md §4.4 says which ear facts discharge it:
`EarSpec.walkTree_guards` (the entries exist, so `MergeOK`/`BoundaryOK`/`RootOK`'s shape) and
`chain_stackDir_const` + `TEntry.OnSide` (every ear hangs on the side `stackDir[topDepth]`, so the
side `finishTstackTop`/`maybeUnwrapNxt`/the block branch/`walkForest` discard is `[]`). -/
theorem walk_sides (g : Graph) (ternarize : Bool) (forest : List DfsTree) (hf : ForestOK g forest)
    (hwf : ∀ t ∈ forest, t.WF []) (hends : ∀ t ∈ forest, t.Ends g)
    (hecov : ∀ e, e < g.ne → e ∈ forest.flatMap DfsTree.edges) :
    SidesForest forest (WalkState.init g ternarize) :=
  walk_sides_of_roots g ternarize forest hf
    (rootsBook_of_state forest [] _ (rootState_init g ternarize) (by simpa using hf) hwf hends
      fun t ht => comp_of_forest hf hwf hends hecov ht)

/-- Exact placement at the end of the walk, and the `tstack` is empty. -/
theorem walk_full (g : Graph) (ternarize : Bool) (forest : List DfsTree) (hf : ForestOK g forest)
    (hwf : ∀ t ∈ forest, t.WF []) (hends : ∀ t ∈ forest, t.Ends g)
    (hecov : ∀ e, e < g.ne → e ∈ forest.flatMap DfsTree.edges) :
    (g.walk ternarize forest).Full g
      (WalkM.Pushed g (fun _ => False) (forest.flatMap DfsTree.verts) (forest.flatMap DfsTree.edges))
      (fun _ => False) ∧ (g.walk ternarize forest).tstack = [] :=
  walkForest_full forest (WalkState.init_full g ternarize) rfl hf.verts_lt hf.edges_lt
    hf.verts_nodup hf.edges_nodup (fun _ _ h => h) (fun _ _ h => h)
    (walk_sides g ternarize forest hf hwf hends hecov)

section Consequences

variable (g : Graph) (ternarize : Bool) (forest : List DfsTree) (hf : ForestOK g forest)
  (hwf : ∀ t ∈ forest, t.WF []) (hends : ∀ t ∈ forest, t.Ends g)
include hf hwf hends hends

/-- The walk ends with an empty `tstack`. -/
theorem walk_tstack_nil (hecov : ∀ e, e < g.ne → e ∈ forest.flatMap DfsTree.edges) :
    (g.walk ternarize forest).tstack = [] := (walk_full g ternarize forest hf hwf hends hecov).2

/-- Every allocated non-root item ends up in some `ch` list. -/
theorem walk_covered
    (hvcov : ∀ v, v < g.nv → v ∈ forest.flatMap DfsTree.verts)
    (hecov : ∀ e, e < g.ne → e ∈ forest.flatMap DfsTree.edges) :
    ∀ i, 0 < i → i < (g.walk ternarize forest).items.size →
      ∃ p, Items.IsParent (g.walk ternarize forest).items p i := fun i hi hlt => by
  obtain ⟨h, ht⟩ := walk_full g ternarize forest hf hwf hends hecov
  refine WalkState.exists_parent_of_cnt ht ?_
  by_cases hn : 1 + g.nv + g.ne ≤ i
  · exact h.placed i hn hlt id
  · refine h.pushed i ?_
    by_cases hv : i < 1 + g.nv
    · exact Or.inr (Or.inl ⟨i - 1, hvcov _ (by omega), by show i = 1 + (i - 1); omega⟩)
    · exact Or.inr (Or.inr ⟨i - (1 + g.nv), hecov _ (by omega),
        by show i = 1 + g.nv + (i - (1 + g.nv)); omega⟩)

/-- The root's children are `vertItem`s of real vertices. -/
theorem walk_root_children (hecov : ∀ e, e < g.ne → e ∈ forest.flatMap DfsTree.edges) :
    ∀ c, Items.IsParent (g.walk ternarize forest).items rootItem c →
      Items.type (g.walk ternarize forest).items c = .V := fun c hc => by
  obtain ⟨v, hv, rfl⟩ := (walk_full g ternarize forest hf hwf hends hecov).1.rootch c hc
  exact (walk_full g ternarize forest hf hwf hends hecov).1.place.vert v hv

/-- Every item is below the root: by `Full.acyc` every item is below a parentless one, and by
`walk_covered` only the root is parentless. -/
theorem walk_reach
    (hvcov : ∀ v, v < g.nv → v ∈ forest.flatMap DfsTree.verts)
    (hecov : ∀ e, e < g.ne → e ∈ forest.flatMap DfsTree.edges) :
    ∀ i, i < (g.walk ternarize forest).items.size →
      Items.Below (g.walk ternarize forest).items rootItem i :=
  (walk_full g ternarize forest hf hwf hends hecov).1.acyc.reach fun r hr0 hrlt hnp =>
    let ⟨p, hp⟩ := walk_covered g ternarize forest hf hwf hends hvcov hecov r hr0 hrlt
    hnp p hp

theorem walk_unique_parent
    (hvcov : ∀ v, v < g.nv → v ∈ forest.flatMap DfsTree.verts)
    (hecov : ∀ e, e < g.ne → e ∈ forest.flatMap DfsTree.edges) :
    ∀ c, 0 < c → c < (g.walk ternarize forest).items.size →
      ∃ p, Items.IsParent (g.walk ternarize forest).items p c ∧
        ∀ p', Items.IsParent (g.walk ternarize forest).items p' c → p' = p := fun c hc hlt =>
  let ⟨p, hp⟩ := walk_covered g ternarize forest hf hwf hends hvcov hecov c hc hlt
  ⟨p, hp, fun p' hp' => walk_parent_unique g ternarize forest hf c p p' hp hp'⟩

/-- Completeness: after the walk every edge is below exactly one child chain from the root — the
`Items.Tree` content of `Items.WF` (false without the forest hypotheses, e.g. for `forest = []`). -/
theorem walk_nodes_partition
    (hvcov : ∀ v, v < g.nv → v ∈ forest.flatMap DfsTree.verts)
    (hecov : ∀ e, e < g.ne → e ∈ forest.flatMap DfsTree.edges) :
    ∀ e, e < g.ne →
      ∃ p, Items.IsParent (g.walk ternarize forest).items p (edgeItem g e) ∧
        ∀ p', Items.IsParent (g.walk ternarize forest).items p' (edgeItem g e) → p' = p := fun e he =>
  walk_unique_parent g ternarize forest hf hwf hends hvcov hecov (edgeItem g e) (by show 0 < 1 + g.nv + e; omega)
    (Nat.lt_of_lt_of_le (by show 1 + g.nv + e < _; omega) (walk_place g ternarize forest hf).size)

end Consequences

/-- The fields of `Items.Tree`, with `v_children` restricted to real vertices. -/
structure WalkTree (g : Graph) (items : Items) : Prop where
  size : 1 + g.nv + g.ne ≤ items.size
  root : items.type rootItem = .F
  vert : ∀ v, v < g.nv → items.type (vertItem v) = .V
  edge : ∀ e, e < g.ne → items.type (edgeItem g e) = .Q
  node : ∀ i, 1 + g.nv + g.ne ≤ i → i < items.size → items.type i ∉ [NodeType.F, .V, .Q]
  ch_lt : ∀ p c, items.IsParent p c → c < items.size
  unique_parent : ∀ c, 0 < c → c < items.size →
    ∃ p, items.IsParent p c ∧ ∀ p', items.IsParent p' c → p' = p
  root_no_parent : ∀ p, ¬ items.IsParent p rootItem
  ch_nodup : ∀ p, (items.ch p).Nodup
  reach : ∀ i, i < items.size → items.Below rootItem i
  v_children : ∀ v c, v < g.nv → items.IsParent (vertItem v) c → items.type c = .Q
  root_children : ∀ c, items.IsParent rootItem c → items.type c = .V ∨ items.type c = .Q

/-- `Items.Tree` for the walk (modulo the admitted `walkTree_book`). -/
theorem walk_tree (g : Graph) (ternarize : Bool) (forest : List DfsTree) (hnv : 0 < g.nv)
    (hb : ∀ t ∈ forest, t.Bounded g.nv g.ne) (hf : ForestOK g forest) (hwf : ∀ t ∈ forest, t.WF [])
    (hends : ∀ t ∈ forest, t.Ends g)
    (hvcov : ∀ v, v < g.nv → v ∈ forest.flatMap DfsTree.verts)
    (hecov : ∀ e, e < g.ne → e ∈ forest.flatMap DfsTree.edges) :
    WalkTree g (g.walk ternarize forest).items :=
  have hcov : ∀ e, e < g.ne → ∃ t ∈ forest, e ∈ t.edges := fun e he =>
    List.mem_flatMap.1 (hecov e he)
  let ht := walk_typing g ternarize forest hnv hb hcov
  { size := ht.size
    root := ht.root
    vert := ht.vert
    edge := ht.edge
    node := ht.node
    ch_lt := walk_ch_lt g ternarize forest hf
    unique_parent := walk_unique_parent g ternarize forest hf hwf hends hvcov hecov
    root_no_parent := walk_root_no_parent g ternarize forest hf
    ch_nodup := walk_ch_nodup g ternarize forest hf
    reach := walk_reach g ternarize forest hf hwf hends hvcov hecov
    v_children := ht.v_children
    root_children := fun c hc => Or.inl (walk_root_children g ternarize forest hf hwf hends hecov c hc) }

end Spqr
