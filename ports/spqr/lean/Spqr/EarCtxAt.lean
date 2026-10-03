import Spqr.EarCtx
import Spqr.EarLoop2

/-!
# `earAt_of_ctx`: the site contract `EarAt` from the between-edges invariant `EarCtx`

PROOF.md §4.2b. The next out's `finishEdge` site is reached from `EarCtx v d done (o :: rest) …`
through `walkOutPre` (direction of `d`, maybe the vertex push) and, for a tree edge, the child's
walk. `earAt_of_ctx_back` derives every `EarFinish` field at a back-edge site; the tree-edge site
is `earAt_of_ctx_tree` (mapped fields derived from the parent context and the child's end-of-outs
context, the rest named per-field admissions).
-/

namespace Spqr

theorem setSides_mem_single {α : Type} (b : Bool) (i : α) :
    i ∈ (setSides b [i] []).1 ++ (setSides b [i] []).2 := by
  cases b <;> simp [setSides]

theorem setSides_mem_single_iff {α : Type} (b : Bool) (i j : α) :
    j ∈ (setSides b [i] []).1 ++ (setSides b [i] []).2 ↔ j = i := by
  cases b <;> simp [setSides]

namespace TEntry

theorem edges_of_root_not_mem {g : Graph} {items : Items} {t : TEntry} {e : Nat}
    (hr : ∀ p, ¬ items.IsParent p (edgeItem g e)) (hm : edgeItem g e ∉ t.spans.1 ++ t.spans.2) :
    ¬ t.edges g items e := by
  rintro ⟨i, hi, hb⟩
  exact hm (hb.eq_of_no_parent hr ▸ hi)

theorem edges_single_entry {g : Graph} {items : Items} (vStart topDepth firstIdx : Nat) (b : Bool)
    (i e : Nat) :
    TEntry.edges g items ⟨vStart, topDepth, firstIdx, setSides b [i] []⟩ e ↔
      items.EdgeBelow g i e := by
  cases b <;> simp [TEntry.edges, setSides]

end TEntry

theorem vertItem_ne_edgeItem' {g : Graph} {x e : Nat} (hx : x < g.nv) : vertItem x ≠ edgeItem g e := by
  intro h
  have h' : (1 + x : Nat) = 1 + g.nv + e := h
  omega

theorem getElem!_set!_ne' {α : Type} [Inhabited α] (a : Array α) (i j : Nat) (v : α) (h : j ≠ i) :
    (a.set! i v)[j]! = a[j]! := by
  rw [Array.set!_eq_setIfInBounds]
  simp only [getElem!_def]
  rw [Array.getElem?_setIfInBounds]
  simp [Ne.symm h]

theorem getElem!_set!_self' {α : Type} [Inhabited α] (a : Array α) (i : Nat) (v : α) (h : i < a.size) :
    (a.set! i v)[i]! = v := by
  rw [Array.set!_eq_setIfInBounds]
  simp only [getElem!_def]
  rw [Array.getElem?_setIfInBounds]
  simp [h]

namespace WalkState

theorem ret_of_lowval_lt' {o : DfsOut} {d : Nat} (h : o.cls.lowval d < d) :
    ∃ lv kind, o.cls = .ret lv kind ∧ lv < d := by
  cases hc : o.cls <;> simp only [OutClass.lowval, hc] at h
  all_goals first | omega | exact ⟨_, _, rfl, h⟩

theorem mem_subEdges_edgesList {outs : List DfsOut} {e : Nat} :
    e ∈ DfsOut.edgesList outs ↔ ∃ o ∈ outs, subEdges o e := by
  induction outs with
  | nil => simp [DfsOut.edgesList]
  | cons o rest ih =>
    cases o with
    | back e' dest cls => simp [DfsOut.edgesList, ih, subEdges, DfsOut.e]
    | tree e' cls child => simp [DfsOut.edgesList, ih, subEdges, DfsOut.e, or_assoc]

theorem mem_vertsList {outs : List DfsOut} {x : Nat} :
    x ∈ DfsOut.vertsList outs ↔ ∃ o ∈ outs, ∃ e cls child, o = .tree e cls child ∧ x ∈ child.verts := by
  induction outs with
  | nil => simp [DfsOut.vertsList]
  | cons o rest ih =>
    cases o with
    | back e' dest cls => simp [DfsOut.vertsList, ih]
    | tree e' cls child =>
      simp only [DfsOut.vertsList, List.mem_append, ih]
      constructor
      · rintro (h | ⟨o, ho, e'', cls', c, rfl, hx⟩)
        · exact ⟨_, List.mem_cons_self, e', cls, child, rfl, h⟩
        · exact ⟨_, List.mem_cons_of_mem _ ho, e'', cls', c, rfl, hx⟩
      · rintro ⟨o, ho, e'', cls', c, rfl, hx⟩
        rcases List.mem_cons.1 ho with h | h
        · cases h; exact .inl hx
        · exact .inr ⟨_, h, _, _, _, rfl, hx⟩

theorem mem_afterVert {done : List (DfsOut × Bool)} {o : DfsOut} (h : o ∈ afterVert done) :
    (o, true) ∈ done := by
  obtain ⟨⟨o', b⟩, hx, rfl⟩ := List.mem_map.1 h
  obtain ⟨hmem, hb⟩ := List.mem_filter.1 hx
  simp at hb; subst hb; exact hmem

/-- Rank-sorted outs: the finished outs at the lowval of a type-1 out are type 1. -/
theorem allType1_of_rank {d : Nat} {done : List (DfsOut × Bool)} {cls : OutClass}
    (hlt : cls.lowval d < d) (ht1 : cls.isType1 = true)
    (hrank : ∀ o' ∈ done, o'.1.cls.rank ≤ cls.rank)
    (hret : ∀ o ∈ afterVert done, o.cls.lowval d < d) :
    allType1 d (cls.lowval d) done := by
  intro o ho hl
  have hr := hrank _ (mem_afterVert ho)
  obtain ⟨lv, k, hc, _⟩ := ret_of_lowval_lt' (hret o ho)
  simp only at hr
  cases hc' : cls <;> simp only [hc', OutClass.lowval] at hlt hl ht1 hr <;> try omega
  rw [hc] at hl hr; simp only [OutClass.rank] at hl hr
  subst hl
  rw [hc]
  cases k <;> cases ‹RetKind› <;> simp_all [RetKind.rank, OutClass.isType1] <;> omega

/-- Rank-sorted outs: no returning out is finished before a boundary one (`hcls`: a `ret` class
returns below `d`, as `classify` guarantees). -/
theorem no_ret_before_boundary {d : Nat} {done : List (DfsOut × Bool)} {cls : OutClass}
    (hge : d ≤ cls.lowval d) (hcls : ∀ lv k, cls = .ret lv k → lv < d)
    (hrank : ∀ o' ∈ done, o'.1.cls.rank ≤ cls.rank) :
    ¬ ∃ o ∈ done, o.1.cls.lowval d < d := by
  rintro ⟨o, ho, hlt⟩
  have hr := hrank o ho
  obtain ⟨lv, k, hc, hlv⟩ := ret_of_lowval_lt' hlt
  rw [hc] at hr
  have hk : 0 ≤ k.rank := Nat.zero_le _
  cases cls with
  | bridge => simp [OutClass.rank] at hr
  | component => simp [OutClass.rank] at hr; omega
  | selfLoop => simp [OutClass.rank] at hr; omega
  | ret lv' k' => exact absurd (hcls lv' k' rfl) (Nat.not_lt.2 hge)

theorem Loop1Range.nil {d : Nat} {hi : List TEntry} (h : Loop1Range d [] hi) : hi = [] := by
  obtain ⟨lo, h, -, -⟩ := h
  exact (List.append_eq_nil_iff.1 h.symm).1

/-- The back-edge site: `walkOutPre` set `stackDir[d]` and pushed `V v` (`L = [V v]`) iff the vertex
entry was missing and the edge returns; every `EarFinish` field follows from `EarCtx`. -/
theorem earAt_back_of_ctx {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s : WalkState} {e dest : Nat} {cls : OutClass}
    (hC : EarCtx v d done (.back e dest cls :: rest) hasVert base bE sv sd s)
    (hv : v < s.g.nv) (hsd : d < s.stackDir.size)
    (hrank : ∀ o' ∈ done, o'.1.cls.rank ≤ cls.rank)
    (hinc : ∀ o' ∈ done, o'.1.e < s.g.ne ∧ s.g.Inc o'.1.e v)
    (hnd : ∀ o' ∈ done, ∀ e', subEdges o'.1 e' → e' ≠ e)
    (hb : cls.isTree = false) (hcls : ∀ lv k, cls = .ret lv k → lv < d)
    (L : List TEntry) (push : Bool)
    (hpush : push = true ↔ hasVert = false ∧ cls.lowval d < d ∧ cls.isType1 = true)
    (hL : L = if push then
      [⟨v, d, s.nxtEdgeIdx, setSides (s.stackDir.set! d (if cls.lowval d ≥ d then false
        else !s.stackDir[cls.lowval d]!))[d]! [vertItem v] []⟩] else []) :
    EarAt v d (.back e dest cls) (L ++ s.tstack).length (hasVert || push)
      { s with stackDir := s.stackDir.set! d (if cls.lowval d ≥ d then false
                 else !s.stackDir[cls.lowval d]!),
               tstack := L ++ s.tstack } := by
  set b := (if cls.lowval d ≥ d then false else !s.stackDir[cls.lowval d]!) with hb_def
  set V : TEntry := ⟨v, d, s.nxtEdgeIdx, setSides (s.stackDir.set! d b)[d]! [vertItem v] []⟩ with hV
  have hT : (DfsOut.cls (.back e dest cls)).isTree = true → False := fun h => by
    simp [DfsOut.cls, hb] at h
  have hmemL : ∀ t ∈ L, t = V ∧ push = true := by
    intro t ht; rw [hL] at ht; split at ht <;> simp_all
  have hLnil : push = false → L = [] := fun hp => by rw [hL, hp]; rfl
  have hLpw : ∀ R : TEntry → TEntry → Prop, L.Pairwise R := by
    intro R; rw [hL]; split <;> simp
  have hVedges : ∀ e', V.edges s.g s.items e' ↔ Items.EdgeBelow s.g s.items (vertItem v) e' :=
    fun e' => TEntry.edges_single_entry _ _ _ _ _ _
  have hVspan : ∀ i, i ∈ V.spans.1 ++ V.spans.2 ↔ i = vertItem v :=
    fun i => setSides_mem_single_iff _ _ _
  have hpushV : push = true → hasVert = false ∧ cls.lowval d < d ∧ cls.isType1 = true := hpush.1
  have hsdk : ∀ k, k ≠ d → (s.stackDir.set! d b)[k]! = s.stackDir[k]! :=
    fun k hk => getElem!_set!_ne' _ _ _ _ hk
  have hsdd : (s.stackDir.set! d b)[d]! = b := getElem!_set!_self' _ _ _ hsd
  have hsub : ∀ e', subEdges (.back e dest cls) e' ↔ e' = e := fun e' => by
    simp [subEdges, DfsOut.e]
  have hqo := hC.q_fresh _ List.mem_cons_self e ((hsub e).2 rfl)
  obtain ⟨top, hts, hTop⟩ := hC.top
  refine ⟨[], L ++ s.tstack, rfl, ?_⟩
  refine
    { tstack := rfl
      back_nil := fun _ => rfl
      base_bot := fun h => (hT h).elim
      sub_bot := by simp
      path := hC.path
      disj := ?disj
      span_disj := ?span_disj
      sub_edges := by simp
      base_disj := ?base_disj
      sub_cover := fun e' _ hne hs => (hne ((hsub e').1 hs)).elim
      loop1_side := fun _ hi hr => by rw [hr.nil]; simp
      loop1_touch := fun _ hi hr => by rw [hr.nil]; simp
      touch_bot := ?touch_bot
      vert := ?vert
      vert_free := ?vert_free
      q_free := ?q_free
      q_root := hqo.1
      v_root := hC.v_root
      p_entry := ?p_entry
      loop1 := fun h => (hT h).elim
      bottom := fun h => (hT h).elim
      loops := fun h => (hT h).elim
      late := fun h => (hT h).elim
      late_fo := fun h => (hT h).elim
      close := fun h => (hT h).elim
      sv_d := hC.sv_d
      sv_child := fun h => (hT h).elim
      path_child := fun h => (hT h).elim
      dir_d := ?dir_d
      boundary := by simp
      base_touch := fun h => (hT h).elim
      bd_noVert := ?bd_noVert
      bd_bridge := fun h => (hT h).elim
      bd_comp := fun h => (hT h).elim
      bd_term := fun h => (hT h).elim
      bd_side := fun h => (hT h).elim
      lower := fun h => (hT h).elim }
  case disj =>
    refine List.pairwise_append.2 ⟨hLpw _, hC.disj, ?_⟩
    intro t ht t' ht' e' he' hte
    obtain ⟨rfl, hp⟩ := hmemL t ht
    exact fun h => hC.vert_disj (hpushV hp).1 t' ht' e' he' h ((hVedges e').1 hte)
  case span_disj =>
    refine List.pairwise_append.2 ⟨hLpw _, hC.span_disj, ?_⟩
    intro t ht t' ht' i hi hi'
    obtain ⟨rfl, hp⟩ := hmemL t ht
    rw [(hVspan i).1 hi] at hi'
    exact absurd (hC.vert_free t' ht' hi') (by rw [(hpushV hp).1]; decide)
  case base_disj =>
    intro t ht e' he' hte hs
    obtain rfl := (hsub e').1 hs
    rcases List.mem_append.1 ht with ht | ht
    · obtain ⟨rfl, -⟩ := hmemL t ht
      obtain ⟨o', ho', -, hso'⟩ := (hC.vert_edges e' he').1 ((hVedges e').1 hte)
      exact hnd o' ho' e' hso' rfl
    · exact TEntry.edges_of_root_not_mem hqo.1 (hqo.2.2 t ht) hte
  case touch_bot =>
    intro t ht ⟨e', he', hte⟩
    rcases List.mem_append.1 ht with ht | ht
    · obtain ⟨rfl, -⟩ := hmemL t ht
      obtain ⟨o', ho', hge, -⟩ := (hC.vert_edges e' he').1 ((hVedges e').1 hte)
      refine ⟨o'.1.e, (hinc o' ho').1, (hVedges _).2 ?_, (hinc o' ho').2⟩
      exact (hC.vert_edges _ (hinc o' ho').1).2 ⟨o', ho', hge, .inl rfl⟩
    · exact hC.touch_bot t ht ⟨e', he', hte⟩
  case vert =>
    intro h
    cases hp : push
    · rw [hp, Bool.or_false] at h
      obtain ⟨above, below, hsplit, -⟩ := hTop.split
      rw [h] at hsplit
      obtain ⟨vt, htop, hle, hmem, -⟩ := hsplit
      exact ⟨vt, List.mem_append_right _ (hts ▸ htop ▸ List.mem_append_left _
        (List.mem_append_right _ List.mem_cons_self)), hle, hmem⟩
    · exact ⟨V, List.mem_append_left _ (by rw [hL, hp]; simp), Nat.le_refl _,
        (hVspan _).2 rfl⟩
  case vert_free =>
    intro t ht hi
    refine ⟨?_, ht⟩
    rcases List.mem_append.1 ht with ht | ht
    · rw [(hmemL t ht).2]; simp
    · rw [hC.vert_free t ht hi]; simp
  case q_free =>
    intro t ht hi
    rcases List.mem_append.1 ht with ht | ht
    · obtain ⟨rfl, -⟩ := hmemL t ht
      exact vertItem_ne_edgeItem' hv ((hVspan _).1 hi).symm
    · exact hqo.2.2 t ht hi
  case p_entry =>
    intro hlt ht1 t ht hvs htd
    simp only [DfsOut.cls] at hlt ht1 htd
    have hne : t.topDepth ≠ d := by omega
    rcases List.mem_append.1 ht with ht | ht
    · obtain ⟨rfl, -⟩ := hmemL t ht
      exact absurd rfl hne
    rw [hts] at ht
    rcases List.mem_append.1 ht with ht | ht
    · obtain ⟨above, below, hsplit, hent, hsingle, -, -, -, hbelow⟩ := hTop.split
      have hA : t ∈ above := by
        cases hasVert
        · obtain ⟨rfl, rfl⟩ := hsplit; exact absurd hvs (hbelow t ht)
        · obtain ⟨vt, htop, -, -, hvt, -⟩ := hsplit
          rw [htop] at ht
          rcases List.mem_append.1 ht with ht | ht
          · exact ht
          rcases List.mem_cons.1 ht with rfl | ht
          · exact absurd (hvt hvs).1 hne
          · exact absurd hvs (hbelow t ht)
      have hS := hsingle t hA (htd ▸ allType1_of_rank hlt ht1 hrank hC.afterVert_ret)
      refine ⟨?_, hS.att⟩
      show ∃ i, t.spans = setSides (s.stackDir.set! d b)[t.topDepth]! [i] [] ∧ _
      rw [hsdk _ hne]; exact hS.item
    · exact absurd hvs (hC.base_bot t ht)
  case dir_d =>
    intro hlt
    simp only [DfsOut.cls] at hlt
    show (s.stackDir.set! d b)[d]! = !(s.stackDir.set! d b)[cls.lowval d]!
    rw [hsdd, hsdk _ (by omega), hb_def]
    simp [Nat.not_le.2 hlt]
  case bd_noVert =>
    intro hge
    simp only [DfsOut.cls] at hge
    have hp : push = false := by
      cases hp : push
      · rfl
      · exact absurd (hpushV hp).2.1 (by omega)
    rw [hp, Bool.or_false]
    cases hhv : hasVert
    · rfl
    · exact absurd (hC.hv_ret hhv) (no_ret_before_boundary hge hcls hrank)

theorem vertItem_inj' {x y : Nat} (h : vertItem x = vertItem y) : x = y := by
  have h' : (1 + x : Nat) = 1 + y := h
  omega

theorem _root_.Spqr.Graph.Touches.congr' {g : Graph} {E₁ E₂ : Nat → Prop}
    (h : ∀ e, e < g.ne → (E₁ e ↔ E₂ e)) {x : Nat} : g.Touches E₁ x ↔ g.Touches E₂ x := by
  constructor
  · rintro ⟨e, he, hE, hi⟩; exact ⟨e, he, (h e he).1 hE, hi⟩
  · rintro ⟨e, he, hE, hi⟩; exact ⟨e, he, (h e he).2 hE, hi⟩

/-- The state after the child's end push: `walkTree` pushes `V y` (direction `dir'`) iff the
child's outs left no vertex entry; `D₃` is the direction array afterwards. -/
def pushEnd (sE : WalkState) (D₃ : Array Bool) (L' : List TEntry) : WalkState :=
  { sE with stackDir := D₃, tstack := L' ++ sE.tstack }

/-- Shape facts about the entries above `base` between outs (`ctxCheck`: `ret_hv`, `t1_flag`,
`t1_above_distinct`/`t1_below`/`t1_vt`, `vt_side_strict`): a returning done out means the vertex
entry is pushed; type-1 returning outs are recorded after it; while every returning done out is
type 1 the vertex entry is the bottom of the top and the ears above it have distinct depths
(decreasing upwards: the sorted outs' lowvals, P-merged when equal); the vertex entry sits on the
side opposite to one of the returning outs' lowval directions. -/
structure CtxShape (v d : Nat) (done : List (DfsOut × Bool)) (hasVert : Bool) (base : List TEntry)
    (s : WalkState) : Prop where
  ret_hv : (∃ o ∈ done, o.1.cls.lowval d < d) → hasVert = true
  t1_flag : ∀ o ∈ done, o.1.cls.lowval d < d → o.1.cls.isType1 = true → o.2 = true
  t1 : hasVert = true → (∀ o ∈ done, o.1.cls.lowval d < d → o.1.cls.isType1 = true) →
    ∃ above vt, s.tstack = above ++ vt :: base ∧ vt.vStart = v ∧
      vertItem v ∈ vt.spans.1 ++ vt.spans.2 ∧
      above.Pairwise (fun t t' => t'.topDepth < t.topDepth)
  vt_side : ∀ top, s.tstack = top ++ base → ∀ t ∈ top, t.vStart = v →
    vertItem v ∈ t.spans.1 ++ t.spans.2 → ∃ o ∈ done, o.1.cls.lowval d < d ∧
      t.spans = setSides (!s.stackDir[o.1.cls.lowval d]!) [vertItem v] []

theorem ctxShape_init {v d : Nat} {base : List TEntry} {s : WalkState} (hts : s.tstack = base) :
    CtxShape v d [] false base s where
  ret_hv := fun ⟨_, ho, _⟩ => nomatch ho
  t1_flag := fun _ ho => nomatch ho
  t1 := fun h => nomatch h
  vt_side := fun top htop t ht => by
    rw [hts] at htop
    have : top = [] := by simpa using htop.symm
    rw [this] at ht
    exact absurd ht List.not_mem_nil

theorem ctxShape_same {v d : Nat} {done : List (DfsOut × Bool)} {hasVert : Bool}
    {base : List TEntry} {s s' : WalkState} {o : DfsOut} {b : Bool}
    (hS : CtxShape v d done hasVert base s) (hol : d ≤ o.cls.lowval d)
    (hts : s'.tstack = s.tstack) (hsd : ∀ k, k < d → s'.stackDir[k]! = s.stackDir[k]!) :
    CtxShape v d (done ++ [(o, b)]) hasVert base s' where
  ret_hv := fun ⟨o', ho', hl⟩ => by
    rcases List.mem_append.1 ho' with h | h
    · exact hS.ret_hv ⟨o', h, hl⟩
    · rw [List.mem_singleton] at h; subst h; exact absurd hl (Nat.not_lt.2 hol)
  t1_flag := fun o' ho' hl ht => by
    rcases List.mem_append.1 ho' with h | h
    · exact hS.t1_flag o' h hl ht
    · rw [List.mem_singleton] at h; subst h; exact absurd hl (Nat.not_lt.2 hol)
  t1 := fun hv ht1 => by
    rw [hts]
    exact hS.t1 hv fun o' ho' hl => ht1 o' (List.mem_append_left _ ho') hl
  vt_side := fun top htop t ht htv hm => by
    obtain ⟨o', ho', hl, hsp⟩ := hS.vt_side top (by rw [← hts]; exact htop) t ht htv hm
    exact ⟨o', List.mem_append_left _ ho', hl, by rw [hsd _ hl]; exact hsp⟩

/-- Two decompositions of a span-disjoint stack around an entry holding the item `i` coincide. -/
theorem split_unique {i : ItemId} : ∀ {a a' : List TEntry} {l b b' : List TEntry} {x x' : TEntry},
    l = a ++ x :: b → l = a' ++ x' :: b' →
    i ∈ x.spans.1 ++ x.spans.2 → i ∈ x'.spans.1 ++ x'.spans.2 →
    l.Pairwise (fun t t' => ∀ j ∈ t.spans.1 ++ t.spans.2, j ∉ t'.spans.1 ++ t'.spans.2) →
    a = a' ∧ x = x' ∧ b = b'
  | [], [], _, _, _, _, _, h, h', _, _, _ => by
    subst h
    obtain ⟨rfl, rfl⟩ := List.cons.inj h'
    exact ⟨rfl, rfl, rfl⟩
  | [], _ :: _, _, _, _, _, _, h, h', hx, hx', hd => by
    subst h
    obtain ⟨rfl, rfl⟩ := List.cons.inj h'
    exact absurd hx' ((List.pairwise_cons.1 hd).1 _
      (List.mem_append_right _ (List.mem_cons_self ..)) i hx)
  | _ :: _, [], _, _, _, _, _, h, h', hx, hx', hd => by
    subst h'
    obtain ⟨rfl, rfl⟩ := List.cons.inj h
    exact absurd hx ((List.pairwise_cons.1 hd).1 _
      (List.mem_append_right _ (List.mem_cons_self ..)) i hx')
  | y :: a, y' :: a', _, _, _, _, _, h, h', hx, hx', hd => by
    subst h
    rw [List.cons_append, List.cons_append, List.cons.injEq] at h'
    obtain ⟨rfl, h'⟩ := h'
    rw [List.cons_append] at hd
    obtain ⟨rfl, rfl, rfl⟩ :=
      split_unique (a := a) (a' := a') rfl h' hx hx' (List.pairwise_cons.1 hd).2
    exact ⟨rfl, rfl, rfl⟩

/-- A tree-edge site `o = (e, cls, node y outs)` of `walkOuts v d`: the parent context at `s`
(before `walkOutPre`), the WF facts of the out (rank order, nodup, endpoints), the vertex push
(`L`, `push`) of `walkOutPre`, the child's end-of-outs context at `sE` over the pushed stack, the
frame of the child's walk (graph, the stack arrays below `d + 1`, the two foreign roots), and the
child's end push (`L'`, `push'`, `D₃`). -/
structure TreeSite (v d : Nat) (done : List (DfsOut × Bool)) (rest : List DfsOut) (hasVert : Bool)
    (base : List TEntry) (bE : List (Nat → Prop)) (sv : List Nat) (sd : List Bool) (s : WalkState)
    (e : Nat) (cls : OutClass) (y : Nat) (outs : List DfsOut) (L : List TEntry) (push : Bool)
    (done' : List (DfsOut × Bool)) (hv' : Bool) (bE' : List (Nat → Prop)) (sv' : List Nat)
    (sd' : List Bool) (sE : WalkState) (dir' : Bool) (L' : List TEntry) (push' : Bool)
    (D₃ : Array Bool) : Prop where
  ctx : EarCtx v d done (.tree e cls (.node y outs) :: rest) hasVert base bE sv sd s
  hv : v < s.g.nv
  hy : y < s.g.nv
  e_lt : e < s.g.ne
  hsd : d < s.stackDir.size
  rank : ∀ o' ∈ done, o'.1.cls.rank ≤ cls.rank
  inc : ∀ o' ∈ done, o'.1.e < s.g.ne ∧ s.g.Inc o'.1.e v
  nd : ∀ o' ∈ done, ∀ e', subEdges o'.1 e' → ¬ subEdges (.tree e cls (.node y outs)) e'
  done_nc : ∀ o' ∈ done, ∀ e', subEdges o'.1 e' → ∀ x, s.g.Inc e' x →
    x ∉ (DfsTree.node y outs).verts
  tree : cls.isTree = true
  cls_ret : ∀ lv k, cls = .ret lv k → lv < d
  v_nc : v ∉ (DfsTree.node y outs).verts
  e_ne : e ∉ (DfsTree.node y outs).edges
  comp : ∀ e', e' < s.g.ne → ∀ x, s.g.Inc e' x → x ∈ (DfsTree.node y outs).verts →
      subEdges (.tree e cls (.node y outs)) e'
  ends : d ≤ cls.lowval d → ∀ e', subEdges (.tree e cls (.node y outs)) e' → e' < s.g.ne →
    ∀ x, s.g.Inc e' x → x = v ∨ x ∈ (DfsTree.node y outs).verts
  hpush : push = true ↔ hasVert = false ∧ cls.lowval d < d ∧ cls.isType1 = true
  hL : L = if push then
    [⟨v, d, s.nxtEdgeIdx, setSides (s.stackDir.set! d (if cls.lowval d ≥ d then false
      else !s.stackDir[cls.lowval d]!))[d]! [vertItem v] []⟩] else []
  ctx' : EarCtx y (d + 1) done' [] hv' (L ++ s.tstack) bE' sv' sd' sE
  hdone' : done'.map (·.1) = outs
  inc' : ∀ o' ∈ done', o'.1.e < s.g.ne ∧ s.g.Inc o'.1.e y
  hbE' : ∀ k, k < (L ++ s.tstack).length → ∀ e', e' < s.g.ne →
    (bE'[k]! e' ↔ (L ++ s.tstack)[k]!.edges s.g s.items e')
  gE : sE.g = s.g
  svlo : ∀ k, k ≤ d → sE.stackVerts[k]! = s.stackVerts[k]!
  sdlo : ∀ k, k ≤ d → sE.stackDir[k]! = (s.stackDir.set! d (if cls.lowval d ≥ d then false
    else !s.stackDir[cls.lowval d]!))[k]!
  q_root : ∀ p, ¬ Items.IsParent sE.items p (edgeItem s.g e)
  v_root : ∀ p, ¬ Items.IsParent sE.items p (vertItem v)
  hpush' : push' = true ↔ hv' = false
  hL' : L' = if push' then [⟨y, d + 1, sE.nxtEdgeIdx, setSides dir' [vertItem y] []⟩] else []
  sd₃ : ∀ k, k ≤ d → D₃[k]! = sE.stackDir[k]!
  hdir' : push' = true → dir' = true
  size_le : s.items.size ≤ sE.items.size
  hsz : 1 + s.g.nv + s.g.ne ≤ s.items.size
  verts_lt : ∀ x ∈ (DfsTree.node y outs).verts, x < s.g.nv
  e_ends : Items.PairEq (y, v) s.g.edges[e]!
  rest_nd : ∀ o' ∈ rest, ∀ e', subEdges o' e' →
    ¬ subEdges (.tree e cls (.node y outs)) e' ∧ e' < s.g.ne
  rest_nv : ∀ o' ∈ rest, ∀ e₁ cls₁ c, o' = .tree e₁ cls₁ c →
    ∀ w ∈ c.verts, w ∉ (DfsTree.node y outs).verts ∧ w < s.g.nv
  /-- The done outs' edges avoid the remaining subtrees' vertices. -/
  done_rest : ∀ o' ∈ done, ∀ e', subEdges o'.1 e' → ∀ x, s.g.Inc e' x →
    ∀ o'' ∈ rest, ∀ e₁ cls₁ child, o'' = .tree e₁ cls₁ child → x ∉ child.verts
  /-- The child's walk touches only its own items and the ones it allocates. -/
  items_kept : ∀ j, j < s.items.size → (∀ x ∈ (DfsTree.node y outs).verts, j ≠ vertItem x) →
    (∀ e' ∈ (DfsTree.node y outs).edges, j ≠ edgeItem s.g e') →
    Items.ch sE.items j = Items.ch s.items j ∧
    ∀ p, Items.IsParent sE.items p j ↔ Items.IsParent s.items p j
  /-- The items spanned by the entries the child's walk left above the pushed stack are the
  child's `V`/`Q` items or allocated during it (`freshCheck`). -/
  child_fresh : ∀ top, sE.tstack = top ++ (L ++ s.tstack) → ∀ t ∈ top,
    ∀ i ∈ t.spans.1 ++ t.spans.2,
    (∃ x ∈ (DfsTree.node y outs).verts, i = vertItem x) ∨
    (∃ e' ∈ (DfsTree.node y outs).edges, i = edgeItem s.g e') ∨ s.items.size ≤ i
  bridge_bd : cls = .bridge → ∀ o' ∈ done', d + 1 ≤ o'.1.cls.lowval (d + 1)
  comp_ret : cls = .component → ∀ o' ∈ done', o'.1.cls.lowval (d + 1) < d + 1 →
    o'.1.cls.lowval (d + 1) = d
  comp_ex : cls = .component → ∃ o' ∈ done', o'.1.cls.lowval (d + 1) < d + 1
  comp_t1 : cls = .component → ∀ o' ∈ done', o'.1.cls.lowval (d + 1) < d + 1 →
    o'.1.cls.isType1 = true
  /-- The child's end-of-outs shape (`ctxCheck` at the child's end). -/
  shape' : CtxShape y (d + 1) done' hv' (L ++ s.tstack) sE
  /-- The child's edges are connected to `y` through the child (DFS tree), within any edge set
  containing them. -/
  c_reach : ∀ E : Nat → Prop, (∀ e' ∈ (DfsTree.node y outs).edges, E e') →
    ∀ e' ∈ (DfsTree.node y outs).edges, Relation.ReflTransGen (s.g.AdjIn E) y (s.g.edges[e']!).1
  edges_lt : ∀ e' ∈ (DfsTree.node y outs).edges, e' < s.g.ne

/-- The stack shape after loops 1–3 of `finishEdge` at a returning tree edge (`retCheck`): relative
to the parent's pre-push state `s`, the state `sX` (after `mergeLate`, or after `closeVert'` when
the vertex entry is pushed) has `R ++ L ++ s.tstack`, the items of `L ++ s.tstack` untouched (except
the child's and `Q e`); the entries of `R` are root-spanned by the child's, `Q e` or fresh items, own
exactly the out's edges, start at `v` or in the child, and touch only `v`, the child and the path at
depths `[lowval, d]`. With the vertex entry (`hv`) `R` is one entry `(v, lowval)` (one item if type
1); without it (type 2) its entries start in the child and the top one holds `e`. -/
structure RetTop (v d : Nat) (s : WalkState) (e : Nat) (cls : OutClass) (y : Nat)
    (outs : List DfsOut) (L : List TEntry) (hv : Bool) (sX : WalkState) (R : List TEntry) :
    Prop where
  tstack : sX.tstack = R ++ (L ++ s.tstack)
  g : sX.g = s.g
  sv : ∀ k, k ≤ d → sX.stackVerts[k]! = s.stackVerts[k]!
  sd : ∀ k, k < d → sX.stackDir[k]! = s.stackDir[k]!
  /-- Without the vertex entry, loop 2 set the direction of `d` opposite to the lowval's. -/
  sdd : hv = false → sX.stackDir[d]! = !s.stackDir[cls.lowval d]!
  size : s.items.size ≤ sX.items.size
  nxt : s.nxtEdgeIdx ≤ sX.nxtEdgeIdx
  kept : ∀ j, j < s.items.size → (∀ x ∈ (DfsTree.node y outs).verts, j ≠ vertItem x) →
    (∀ e' ∈ (DfsTree.node y outs).edges, j ≠ edgeItem s.g e') → j ≠ edgeItem s.g e →
    Items.ch sX.items j = Items.ch s.items j ∧
    ∀ p, Items.IsParent sX.items p j ↔ Items.IsParent s.items p j
  ch_lt : ∀ p c, Items.IsParent sX.items p c → c < sX.items.size
  span_lt : ∀ t ∈ R, ∀ i ∈ t.spans.1 ++ t.spans.2, i < sX.items.size
  span_root : ∀ t ∈ R, ∀ i ∈ t.spans.1 ++ t.spans.2, ∀ p, ¬ Items.IsParent sX.items p i
  span_new : ∀ t ∈ R, ∀ i ∈ t.spans.1 ++ t.spans.2,
    (∀ t' ∈ L ++ s.tstack, i ∉ t'.spans.1 ++ t'.spans.2) ∧
    (i < s.items.size → (∃ x ∈ (DfsTree.node y outs).verts, i = vertItem x) ∨
      (∃ e' ∈ (DfsTree.node y outs).edges, i = edgeItem s.g e') ∨ i = edgeItem s.g e)
  touch_bot : ∀ t ∈ R, (∃ e', e' < s.g.ne ∧ t.edges s.g sX.items e') →
    s.g.Touches (t.edges s.g sX.items) t.vStart
  edges : ∀ e', e' < s.g.ne →
    ((∃ t ∈ R, t.edges s.g sX.items e') ↔ subEdges (.tree e cls (.node y outs)) e')
  touch : ∀ t ∈ R, ∀ x, s.g.Touches (t.edges s.g sX.items) x →
    x = v ∨ x ∈ (DfsTree.node y outs).verts ∨ ∃ k, cls.lowval d ≤ k ∧ k ≤ d ∧ x = s.stackVerts[k]!
  disj : R.Pairwise fun t t' => ∀ e', e' < s.g.ne → t.edges s.g sX.items e' →
    ¬ t'.edges s.g sX.items e'
  span_disj : R.Pairwise fun t t' => ∀ i ∈ t.spans.1 ++ t.spans.2, i ∉ t'.spans.1 ++ t'.spans.2
  vstart : ∀ t ∈ R, t.vStart = v ∨ t.vStart ∈ (DfsTree.node y outs).verts
  touch_k : ∀ t ∈ R, ∀ k, k ≤ d → s.g.Touches (t.edges s.g sX.items) s.stackVerts[k]! →
    t.topDepth ≤ k
  vert : hv = true → ∃ x, R = [x] ∧ x.vStart = v ∧ x.topDepth = cls.lowval d ∧
    getSide x.spans (!s.stackDir[cls.lowval d]!) = [] ∧
    s.nxtEdgeIdx ≤ x.firstIdx ∧ x.firstIdx < sX.nxtEdgeIdx ∧
    s.g.Touches (x.edges s.g sX.items) s.stackVerts[cls.lowval d]! ∧
    (∀ w, s.g.Touches (x.edges s.g sX.items) w → w ≠ v →
      ¬ s.g.Interior (x.edges s.g sX.items) w →
      ∃ k, cls.lowval d ≤ k ∧ k ≤ d ∧ w = s.stackVerts[k]!) ∧
    (cls.isType1 = true → (∃ i, x.spans = setSides s.stackDir[cls.lowval d]! [i] []) ∧
      ∀ w, s.g.Touches (x.edges s.g sX.items) w →
        w = v ∨ w = s.stackVerts[cls.lowval d]! ∨ s.g.Interior (x.edges s.g sX.items) w)
  noVert : hv = false → (∃ c R', R = c :: R' ∧ c.edges s.g sX.items e) ∧
    ∀ t ∈ R, t.vStart ≠ v ∧ ∀ k, k ≤ d → s.stackVerts[k]! ≠ t.vStart

theorem tree_comp_shape_of_shape {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s : WalkState} {e : Nat} {cls : OutClass} {y : Nat} {outs : List DfsOut} {L : List TEntry}
    {push : Bool} {done' : List (DfsOut × Bool)} {hv' : Bool} {bE' : List (Nat → Prop)}
    {sv' : List Nat} {sd' : List Bool} {sE : WalkState} {dir' : Bool} {L' : List TEntry}
    {push' : Bool} {D₃ : Array Bool}
    (H : TreeSite v d done rest hasVert base bE sv sd s e cls y outs L push done' hv' bE' sv' sd'
      sE dir' L' push' D₃)
    (hS' : CtxShape y (d + 1) done' hv' (L ++ s.tstack) sE)
    (hex : ∃ o' ∈ done', o'.1.cls.lowval (d + 1) < d + 1)
    (ht1 : ∀ o' ∈ done', o'.1.cls.lowval (d + 1) < d + 1 → o'.1.cls.isType1 = true)
    (hc : cls = .component) :
    hv' = true ∧ ∃ t₁ f₂,
      sE.tstack = t₁ :: ⟨y, d + 1, f₂, ([], [vertItem y])⟩ :: (L ++ s.tstack) ∧
      t₁.vStart = y ∧ t₁.topDepth = d ∧ t₁.spans.2 = [] := by
  have hv't : hv' = true := hS'.ret_hv hex
  have hsdd : sE.stackDir[d]! = false := by
    rw [H.sdlo d (Nat.le_refl _), getElem!_set!_self' _ _ _ H.hsd, hc]
    simp [OutClass.lowval]
  have hC' := H.ctx'
  obtain ⟨top, htop, hCT⟩ := hC'.top
  obtain ⟨above₀, below₀, hsp, hCE, -, -, hLow, -, -⟩ := hCT.split
  simp only [hv't, ↓reduceIte] at hsp
  obtain ⟨vt₀, htv₀, -, hvtm₀, hvtd, -, hcov, -⟩ := hsp
  have hab₀ : sE.tstack = above₀ ++ vt₀ :: (below₀ ++ (L ++ s.tstack)) := by
    rw [htop, htv₀]; simp
  obtain ⟨above, vt, hab, hvtv, hvtm, hpw⟩ := hS'.t1 hv't ht1
  obtain ⟨rfl, rfl, hbelow⟩ := split_unique hab hab₀ hvtm hvtm₀ hC'.span_disj
  have hbelow0 : below₀ = [] := by simpa using hbelow.symm
  subst hbelow0
  have hvtd' : vt.topDepth = d + 1 := (hvtd hvtv).1
  have hdep : ∀ t ∈ above, t.topDepth = d := fun t ht => by
    obtain ⟨o, ho, hto⟩ := hLow t ht
    rw [hto]
    exact H.comp_ret hc (o, true) (mem_afterVert ho) (by rw [← hto]; exact (hCE t ht).depth)
  obtain ⟨o', ho', hl'⟩ := hex
  have ho'a : o'.1 ∈ afterVert done' := by
    have h2 := hS'.t1_flag o' ho' hl' (ht1 o' ho' hl')
    exact List.mem_map.2 ⟨o', List.mem_filter.2 ⟨ho', by simp [h2]⟩, rfl⟩
  have hlo' : o'.1.cls.lowval (d + 1) = d := H.comp_ret hc o' ho' hl'
  rcases above with _ | ⟨t₁, _ | ⟨t₂, above''⟩⟩
  · exfalso
    rcases hcov o'.1 ho'a with ⟨t, ht, -⟩ | h
    · exact absurd ht List.not_mem_nil
    · rw [hvtd', hlo'] at h; omega
  · refine ⟨hv't, t₁, vt.firstIdx, ?_, (hCE t₁ (List.mem_cons_self ..)).vStart,
      hdep t₁ (List.mem_cons_self ..), ?_⟩
    · obtain ⟨o'', ho'', hl'', hsp⟩ := hS'.vt_side top htop vt
        (by rw [htv₀]; simp) hvtv hvtm
      have hlo'' : o''.1.cls.lowval (d + 1) = d := H.comp_ret hc o'' ho'' hl''
      rw [hlo'', hsdd] at hsp
      rw [hab]
      rcases vt with ⟨a, b, c, sp⟩
      simp only at hvtv hvtd' hsp
      subst hvtv hvtd' hsp
      rfl
    · have hside := (hCE t₁ (List.mem_cons_self ..)).side
      rw [hdep t₁ (List.mem_cons_self ..), hsdd] at hside
      exact hside
  · exfalso
    have h12 := (List.pairwise_cons.1 hpw).1 t₂ (List.mem_cons_self ..)
    rw [hdep t₁ (List.mem_cons_self ..),
      hdep t₂ (List.mem_cons_of_mem _ (List.mem_cons_self ..))] at h12
    exact Nat.lt_irrefl _ h12

theorem isParent_lt_size {items : Items} {p c : ItemId} (h : Items.IsParent items p c) :
    p < items.size := by
  by_contra hp
  simp [Items.IsParent, Items.ch, Array.getElem?_eq_none (Nat.le_of_not_lt hp)] at h

theorem mem_l2Cur_spans (c₀ : TEntry) (R : List TEntry) (i : ItemId) :
    i ∈ (l2Cur c₀ R).spans.1 ++ (l2Cur c₀ R).spans.2 ↔
      i ∈ c₀.spans.1 ++ c₀.spans.2 ∨ ∃ t ∈ R, i ∈ t.spans.1 ++ t.spans.2 := by
  induction R generalizing c₀ with
  | nil => simp [l2Cur]
  | cons t R ih =>
    rw [l2Cur_cons, ih, TEntry.mem_mergeInto]
    simp [List.mem_cons, exists_eq_or_imp, or_assoc]

theorem l2Cur_top_le (c₀ : TEntry) (R : List TEntry) :
    (l2Cur c₀ R).topDepth ≤ c₀.topDepth ∧ ∀ t ∈ R, (l2Cur c₀ R).topDepth ≤ t.topDepth := by
  induction R generalizing c₀ with
  | nil => simp [l2Cur]
  | cons t R ih =>
    rw [l2Cur_cons]
    obtain ⟨h1, h2⟩ := ih (TEntry.mergeInto c₀ t)
    refine ⟨Nat.le_trans h1 (Nat.min_le_right _ _), fun u hu => ?_⟩
    rcases List.mem_cons.1 hu with rfl | hu
    · exact Nat.le_trans h1 (Nat.min_le_left _ _)
    · exact h2 u hu

theorem le_l2Cur_top (c₀ : TEntry) (R : List TEntry) (x : Nat) (h0 : x ≤ c₀.topDepth)
    (h : ∀ t ∈ R, x ≤ t.topDepth) : x ≤ (l2Cur c₀ R).topDepth := by
  induction R generalizing c₀ with
  | nil => simpa [l2Cur] using h0
  | cons t R ih =>
    rw [l2Cur_cons]
    exact ih _ (Nat.le_min.2 ⟨h t (List.mem_cons_self ..), h0⟩) fun u hu => h u (List.mem_cons_of_mem _ hu)

section
variable {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut} {hasVert : Bool}
  {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool} {s : WalkState}
  {e : Nat} {cls : OutClass} {y : Nat} {outs : List DfsOut} {L : List TEntry} {push : Bool}
  {done' : List (DfsOut × Bool)} {hv' : Bool} {bE' : List (Nat → Prop)} {sv' : List Nat}
  {sd' : List Bool} {sE : WalkState} {dir' : Bool} {L' : List TEntry} {push' : Bool}
  {D₃ : Array Bool}
  (H : TreeSite v d done rest hasVert base bE sv sd s e cls y outs L push done' hv' bE' sv' sd' sE
    dir' L' push' D₃)
include H

local notation "o₀" => DfsOut.tree e cls (DfsTree.node y outs)
local notation "s₃" => pushEnd sE D₃ L'

/-- At a component edge the child's end-of-outs stack is `[(y, d), V y]` above the parent's, its
vertex entry already pushed (so there is no end push), the `(y, d)` entry on side 1 and `V y` on
side 2 — the two entries `finishBoundary` pops (`compEndCheck`; from the child's end shape: every
returning child out is type 1 at lowval `d`). -/
theorem tree_comp_shape : cls = .component → hv' = true ∧ ∃ t₁ f₂,
    sE.tstack = t₁ :: ⟨y, d + 1, f₂, ([], [vertItem y])⟩ :: (L ++ s.tstack) ∧
    t₁.vStart = y ∧ t₁.topDepth = d ∧ t₁.spans.2 = [] := fun hc =>
  tree_comp_shape_of_shape H H.shape' (H.comp_ex hc) (H.comp_t1 hc) hc

/-! ### Named admissions (PROOF.md §4.2b): the `EarFinish` fields at a tree-edge site not yet
derived from the child's end-of-outs context. Each is the field verbatim, over the site's `sub`. -/
section
variable (sub : List TEntry) (hsub : (pushEnd sE D₃ L').tstack = sub ++ (L ++ s.tstack))
include hsub

theorem earAt_tree_loop1_side : (o₀).cls.lowval d < d → ∀ hi, Loop1Range d sub hi →
    ∀ t ∈ hi, getSide t.spans (!(s₃).stackDir[d]!) = [] := by
  sorry

theorem earAt_tree_loop1_touch : (o₀).cls.lowval d < d → ∀ hi, Loop1Range d sub hi → ∀ t ∈ hi,
    (∃ e, e < (s₃).g.ne ∧ t.edges (s₃).g (s₃).items e) → t.topDepth ≤ d + 1 →
    (s₃).g.Touches (t.edges (s₃).g (s₃).items) (s₃).stackVerts[t.topDepth]! := by
  sorry

theorem earAt_tree_loop1 : (o₀).cls.isTree = true → (o₀).cls.lowval d < d → ∀ hi, Loop1Range d sub hi →
    Loop1Spec d o₀ s₃ hi := by
  sorry

theorem earAt_tree_bottom : (o₀).cls.isTree = true → (o₀).cls.lowval d < d →
    ∃ mid py vy, sub = mid ++ [py, vy] ∧ EarBottom d ((o₀).cls.lowval d) s₃ py vy := by
  sorry

theorem earAt_tree_loops : (o₀).cls.isTree = true → (o₀).cls.lowval d < d →
    ∃ c mid py vy, (feS₂ d o₀ s₃).tstack = c :: mid ++ [py, vy] ++ (L ++ s.tstack) ∧
      (∃ mid₀, sub = mid₀ ++ [py, vy]) ∧
      (o₀).cls.lowval d ≤ c.topDepth ∧ c.topDepth ≤ d ∧
      ((o₀).cls.isType1 = true → mid = [] ∧
        (feS₂ d o₀ s₃).g.Touches (c.edges (feS₂ d o₀ s₃).g (feS₂ d o₀ s₃).items) v ∧
        (feS₂ d o₀ s₃).g.Touches (c.edges (feS₂ d o₀ s₃).g (feS₂ d o₀ s₃).items) (o₀).dest ∧
        ∀ e, e < (s₃).g.ne → subEdges o₀ e →
          c.edges (feS₂ d o₀ s₃).g (feS₂ d o₀ s₃).items e ∨
          py.edges (feS₂ d o₀ s₃).g (feS₂ d o₀ s₃).items e ∨
          vy.edges (feS₂ d o₀ s₃).g (feS₂ d o₀ s₃).items e) := by
  sorry

theorem earAt_tree_late : (o₀).cls.isTree = true → (o₀).cls.lowval d < d →
    ∃ c₀ R, EarLate d (feS₁ d o₀ s₃) c₀ R := by
  sorry

theorem earAt_tree_late_fo : (o₀).cls.isTree = true → (o₀).cls.lowval d < d →
    ∀ t ∈ (feS₁ d o₀ s₃).tstack.getLast?, t.firstIdx ≤ (feS₁ d o₀ s₃).firstOccurrence[d]! := by
  sorry

theorem earAt_tree_close : (o₀).cls.isTree = true → (o₀).cls.lowval d < d →
    ∃ c mid py vy, EarClose v d ((o₀).cls.lowval d) o₀ (hasVert || push) (L ++ s.tstack) s₃
      (feS₂ d o₀ s₃) c mid py vy := by
  sorry

theorem earAt_tree_bd_bridge : (o₀).cls.isTree = true → (o₀).cls.lowval d = d + 1 →
    ∃ t, sub = [t] ∧ t.vStart = (o₀).dest ∧ t.topDepth = d + 1 := by
  intro _ hl
  change OutClass.lowval d cls = d + 1 at hl
  have hb : cls = .bridge := by
    cases hc : cls with
    | bridge => rfl
    | component => rw [hc] at hl; simp [OutClass.lowval] at hl
    | selfLoop => rw [hc] at hl; simp [OutClass.lowval] at hl
    | ret lv k =>
      have := H.cls_ret lv k hc
      rw [hc] at hl; simp [OutClass.lowval] at hl; omega
  have hC' := H.ctx'
  have hbd := H.bridge_bd hb
  have hv'f : hv' = false := by
    cases hh : hv'
    · rfl
    · obtain ⟨o', ho', hl'⟩ := hC'.hv_ret hh
      exact absurd (hbd o' ho') (Nat.not_le.2 hl')
  have hpf : push = false := by
    cases hp : push
    · rfl
    · have := (H.hpush.1 hp).2.1
      omega
  have hp' : push' = true := H.hpush'.2 hv'f
  have hdir : dir' = true := H.hdir' hp'
  have hL : L = [] := by rw [H.hL, hpf]; rfl
  have hL' : L' = [⟨y, d + 1, sE.nxtEdgeIdx, setSides true [vertItem y] []⟩] := by
    rw [H.hL', hp', hdir]; rfl
  have hsEts : sE.tstack = s.tstack := by
    obtain ⟨top', htop', hCT'⟩ := hC'.top
    rw [hL, List.nil_append] at htop'
    cases top' with
    | nil => simpa using htop'
    | cons t tl =>
      obtain ⟨o', ho', hl'⟩ := hCT'.ret (List.cons_ne_nil _ _)
      exact absurd (hbd o' ho') (Nat.not_le.2 hl')
  have h : L' ++ sE.tstack = sub ++ (L ++ s.tstack) := hsub
  rw [hL, hL', hsEts, List.nil_append] at h
  exact ⟨_, (List.append_cancel_right h).symm, rfl, rfl⟩

theorem earAt_tree_bd_comp : (o₀).cls.isTree = true → d ≤ (o₀).cls.lowval d → (o₀).cls.lowval d ≠ d + 1 →
    ∃ t₁ t₂, sub = [t₁, t₂] ∧ t₁.vStart = (o₀).dest ∧ t₁.topDepth = (o₀).cls.lowval d ∧
      t₂.vStart = (o₀).dest ∧ t₂.topDepth = d + 1 := by
  intro ht hge hne
  change OutClass.isTree cls = true at ht
  change d ≤ OutClass.lowval d cls at hge
  change OutClass.lowval d cls ≠ d + 1 at hne
  have hc : cls = .component := by
    cases hc : cls with
    | component => rfl
    | bridge => rw [hc] at hne; simp [OutClass.lowval] at hne
    | selfLoop => rw [hc] at ht; simp [OutClass.isTree] at ht
    | ret lv k =>
      have := H.cls_ret lv k hc
      rw [hc] at hge; simp [OutClass.lowval] at hge; omega
  obtain ⟨hv't, t₁, f₂, hts, hvs, htop, -⟩ := tree_comp_shape H hc
  have hp' : push' = false := by
    cases hp : push'
    · rfl
    · rw [H.hpush'.1 hp] at hv't; cases hv't
  have hL' : L' = [] := by rw [H.hL', hp']; rfl
  have h : L' ++ sE.tstack = sub ++ (L ++ s.tstack) := hsub
  rw [hL', List.nil_append, hts] at h
  have h' : [t₁, ⟨y, d + 1, f₂, ([], [vertItem y])⟩] ++ (L ++ s.tstack) = sub ++ (L ++ s.tstack) := h
  refine ⟨t₁, _, (List.append_cancel_right h').symm, hvs, ?_, rfl, rfl⟩
  rw [htop]; subst hc; rfl

theorem earAt_tree_bd_term : (o₀).cls.isTree = true → d ≤ (o₀).cls.lowval d → ∀ u ∈ (s₃).tstack.tail,
    (s₃).g.Touches (u.edges (s₃).g (s₃).items) v → u.vStart = v ∨ u.topDepth ≤ d := by
  intro _ _ u hu htouch
  have hC := H.ctx
  have hC' := H.ctx'
  obtain ⟨top', htop', hCT'⟩ := hC'.top
  obtain ⟨top, htop, hCT⟩ := hC.top
  have hvd : sE.stackVerts[d]! = v := (H.svlo d (Nat.le_refl _)).trans hC.sv_d
  have hmem : u ∈ L' ++ sE.tstack := List.mem_of_mem_tail hu
  have htE : sE.g.Touches (u.edges sE.g sE.items) sE.stackVerts[d]! := by rw [hvd]; exact htouch
  rcases List.mem_append.1 hmem with h | h
  · exfalso
    rw [H.hL'] at h
    split at h
    · rw [List.mem_singleton] at h
      subst h
      obtain ⟨e', he', ⟨i, hi, hb⟩, hinc⟩ := htE
      refine hC'.vert_touch d (Nat.lt_succ_self d) ⟨e', he', ?_, hinc⟩
      cases dir' <;> simp [setSides] at hi <;> subst hi <;> exact hb
    · nomatch h
  · rw [htop'] at h
    rcases List.mem_append.1 h with h | h
    · exact .inr (hCT'.touch_k u h d (Nat.le_succ d) htE)
    · obtain ⟨k, hk, hku⟩ := List.mem_iff_getElem.1 h
      have hk! : (L ++ s.tstack)[k]! = u := by rw [getElem!_pos _ k hk]; exact hku
      obtain ⟨e', he', hte, hinc⟩ := htE
      have hbE : bE'[k]! e' := (hC'.base_edges.2 k hk e' he').1 (by rw [hk!]; exact hte)
      rw [H.gE] at he' hinc; rw [H.svlo d (Nat.le_refl _)] at hinc
      have hte' : u.edges s.g s.items e' := by rw [← hk!]; exact (H.hbE' k hk e' he').1 hbE
      have hts : s.g.Touches (u.edges s.g s.items) s.stackVerts[d]! := ⟨e', he', hte', hinc⟩
      rcases List.mem_append.1 h with h | h
      · left; rw [H.hL] at h; split at h
        · rw [List.mem_singleton] at h; subst h; rfl
        · nomatch h
      · rw [htop] at h
        rcases List.mem_append.1 h with h | h
        · exact .inr (hCT.touch_k u h d (Nat.le_refl _) hts)
        · obtain ⟨j, hj, hju⟩ := List.mem_iff_getElem.1 h
          have hj! : base[j]! = u := by rw [getElem!_pos base j hj]; exact hju
          exact absurd rfl (hC.base_touch j hj v ⟨e', he',
            (hC.base_edges.2 j hj e' he').1 (by rw [hj!]; exact hte'),
            by rw [← hC.sv_d]; exact hinc⟩).1

theorem earAt_tree_bd_side : (o₀).cls.isTree = true → d ≤ (o₀).cls.lowval d →
    if (o₀).cls.lowval d = d + 1 then ∀ t ∈ (s₃).tstack.head?, t.spans.1 = []
    else (∀ b ∈ (s₃).tstack.head?, b.spans.2 = []) ∧ ∀ t ∈ (s₃).tstack.tail.head?, t.spans.1 = [] := by
  intro ht hge
  change OutClass.isTree cls = true at ht
  change d ≤ OutClass.lowval d cls at hge
  show (if OutClass.lowval d cls = d + 1 then _ else _)
  split
  · rename_i hl
    obtain ⟨t, hsub', -, -⟩ := earAt_tree_bd_bridge H sub hsub ht hl
    intro t' ht'
    have h0 : L' ++ sE.tstack = sub ++ (L ++ s.tstack) := hsub
    change t' ∈ (L' ++ sE.tstack).head? at ht'
    rw [h0, hsub'] at ht'
    simp only [List.cons_append, List.nil_append, List.head?_cons, Option.mem_def,
      Option.some.injEq] at ht'
    subst ht'
    have hb : cls = .bridge := by
      cases hc : cls with
      | bridge => rfl
      | component => rw [hc] at hl; simp [OutClass.lowval] at hl
      | selfLoop => rw [hc] at hl; simp [OutClass.lowval] at hl
      | ret lv k =>
        have := H.cls_ret lv k hc
        rw [hc] at hl; simp [OutClass.lowval] at hl; omega
    have hv'f : hv' = false := by
      cases hh : hv'
      · rfl
      · obtain ⟨o', ho', hl'⟩ := H.ctx'.hv_ret hh
        exact absurd (H.bridge_bd hb o' ho') (Nat.not_le.2 hl')
    have hp' : push' = true := H.hpush'.2 hv'f
    have hdir : dir' = true := H.hdir' hp'
    have hL' : L' = [⟨y, d + 1, sE.nxtEdgeIdx, setSides true [vertItem y] []⟩] := by
      rw [H.hL', hp', hdir]; rfl
    have h : L' ++ sE.tstack = [t] ++ (L ++ s.tstack) := h0.trans (by rw [hsub'])
    rw [hL'] at h
    have := (List.cons.inj h).1
    rw [← this]; rfl
  · rename_i hl
    obtain ⟨t₁, t₂, hsub', -, -, -, -⟩ := earAt_tree_bd_comp H sub hsub ht hge hl
    have hc : cls = .component := by
      cases hc : cls with
      | component => rfl
      | bridge => rw [hc] at hl; simp [OutClass.lowval] at hl
      | selfLoop => rw [hc] at ht; simp [OutClass.isTree] at ht
      | ret lv k =>
        have := H.cls_ret lv k hc
        rw [hc] at hge; simp [OutClass.lowval] at hge; omega
    obtain ⟨hv't, t₁', f₂, hts, -, -, hs2⟩ := tree_comp_shape H hc
    have hp' : push' = false := by
      cases hp : push'
      · rfl
      · rw [H.hpush'.1 hp] at hv't; cases hv't
    have hL' : L' = [] := by rw [H.hL', hp']; rfl
    have h0 : L' ++ sE.tstack = sub ++ (L ++ s.tstack) := hsub
    have h : L' ++ sE.tstack = [t₁, t₂] ++ (L ++ s.tstack) := h0.trans (by rw [hsub'])
    rw [hL', List.nil_append, hts] at h
    have h1 := (List.cons.inj h).1
    have h2 := (List.cons.inj (List.cons.inj h).2).1
    change (∀ b ∈ (L' ++ sE.tstack).head?, b.spans.2 = []) ∧
      ∀ t ∈ (L' ++ sE.tstack).tail.head?, t.spans.1 = []
    rw [hL', List.nil_append, hts]
    simp only [List.head?_cons, List.tail_cons, Option.mem_def, Option.some.injEq]
    exact ⟨fun b hb => by rw [← hb]; exact hs2, fun t ht => by rw [← ht]⟩

theorem earAt_tree_lower : (o₀).cls.isTree = true →
    (after (finishEdge v d o₀ (L ++ s.tstack).length (hasVert || push)) s₃).g = (s₃).g ∧
    (after (finishEdge v d o₀ (L ++ s.tstack).length (hasVert || push)) s₃).stackVerts =
      (s₃).stackVerts ∧
    ∀ above t below,
      (after (finishEdge v d o₀ (L ++ s.tstack).length (hasVert || push)) s₃).tstack =
        above ++ t :: below →
      (s₃).g.Touches (t.edges (s₃).g
        (after (finishEdge v d o₀ (L ++ s.tstack).length (hasVert || push)) s₃).items) (o₀).dest →
      (s₃).g.Interior (t.edges (s₃).g
        (after (finishEdge v d o₀ (L ++ s.tstack).length (hasVert || push)) s₃).items) (o₀).dest ∨
      (o₀).dest = t.vStart ∨ ∃ t' ∈ above, (o₀).dest = t'.vStart := by
  sorry

end

/-- The tree-edge site: every mapped `EarFinish` field from the parent context, the child's
end-of-outs context and the frame; the rest are the named admissions above. -/
theorem earAt_tree_of_ctx :
    EarAt v d (DfsOut.tree e cls (.node y outs)) (L ++ s.tstack).length (hasVert || push)
      (pushEnd sE D₃ L') := by
  set o : DfsOut := .tree e cls (.node y outs) with ho
  set child : DfsTree := .node y outs with hchild
  set b := (if cls.lowval d ≥ d then false else !s.stackDir[cls.lowval d]!) with hb_def
  set V : TEntry := ⟨v, d, s.nxtEdgeIdx, setSides (s.stackDir.set! d b)[d]! [vertItem v] []⟩ with hV
  set Vy : TEntry := ⟨y, d + 1, sE.nxtEdgeIdx, setSides dir' [vertItem y] []⟩ with hVy
  have hgE := H.gE
  have hne : sE.g.ne = s.g.ne := by rw [hgE]
  have hnv : sE.g.nv = s.g.nv := by rw [hgE]
  have hT : o.cls.isTree = false → False := fun h => by simp [ho, DfsOut.cls, H.tree] at h
  have hmemL : ∀ t ∈ L, t = V ∧ push = true := by
    intro t ht; rw [H.hL] at ht; split at ht <;> simp_all
  have hmemL' : ∀ t ∈ L', t = Vy ∧ push' = true := by
    intro t ht; rw [H.hL'] at ht; split at ht <;> simp_all
  have hLpw : ∀ R : TEntry → TEntry → Prop, L.Pairwise R := by
    intro R; rw [H.hL]; split <;> simp
  have hLpw' : ∀ R : TEntry → TEntry → Prop, L'.Pairwise R := by
    intro R; rw [H.hL']; split <;> simp
  have hVedges : ∀ e', V.edges s.g s.items e' ↔ Items.EdgeBelow s.g s.items (vertItem v) e' :=
    fun e' => TEntry.edges_single_entry _ _ _ _ _ _
  have hVyedges : ∀ e', Vy.edges sE.g sE.items e' ↔ Items.EdgeBelow sE.g sE.items (vertItem y) e' :=
    fun e' => TEntry.edges_single_entry _ _ _ _ _ _
  have hVspan : ∀ i, i ∈ V.spans.1 ++ V.spans.2 ↔ i = vertItem v :=
    fun i => setSides_mem_single_iff _ _ _
  have hVyspan : ∀ i, i ∈ Vy.spans.1 ++ Vy.spans.2 ↔ i = vertItem y :=
    fun i => setSides_mem_single_iff _ _ _
  have hsub : ∀ e', subEdges o e' ↔ e' = e ∨ e' ∈ child.edges := fun e' => by
    simp [ho, subEdges, DfsOut.e]
  have hsubChild : ∀ e', (∃ x ∈ done', subEdges x.1 e') → e' ∈ child.edges := by
    rintro e' ⟨x, hx, hs⟩
    exact mem_subEdges_edgesList.2 ⟨x.1, by rw [← H.hdone']; exact List.mem_map_of_mem hx, hs⟩
  have hyv : y ∈ child.verts := List.mem_cons_self
  have hvy : v ≠ y := fun h => H.v_nc (by rw [h]; exact hyv)
  have hfreshY := H.ctx.v_fresh o List.mem_cons_self e cls child rfl
  have hqo := H.ctx.q_fresh o List.mem_cons_self e (.inl rfl)
  have hpushV := H.hpush.1
  have hpushV' := H.hpush'.1
  have hsdk : ∀ k, k ≠ d → (s.stackDir.set! d b)[k]! = s.stackDir[k]! :=
    fun k hk => getElem!_set!_ne' _ _ _ _ hk
  have hsdd : (s.stackDir.set! d b)[d]! = b := getElem!_set!_self' _ _ _ H.hsd
  have hbaseE : ∀ t ∈ L ++ s.tstack, ∀ e', e' < s.g.ne →
      (t.edges sE.g sE.items e' ↔ t.edges s.g s.items e') := by
    intro t ht e' he'
    obtain ⟨k, hk, rfl⟩ := List.getElem_of_mem ht
    rw [← getElem!_pos (L ++ s.tstack) k hk]
    exact (H.ctx'.base_edges.2 k hk e' (by rw [hne]; exact he')).trans (H.hbE' k hk e' he')
  have hbaseE' : ∀ t ∈ L ++ s.tstack, ∀ e', e' < s.g.ne →
      (t.edges s.g sE.items e' ↔ t.edges s.g s.items e') := by
    intro t ht e' he'; have := hbaseE t ht e' he'; rwa [hgE] at this
  have hbaseNT : ∀ u ∈ L ++ s.tstack, ∀ x ∈ child.verts, ¬ s.g.Touches (u.edges s.g s.items) x := by
    intro u hu x hx ⟨e', he', hue, hix⟩
    rcases List.mem_append.1 hu with hu | hu
    · obtain ⟨rfl, -⟩ := hmemL u hu
      obtain ⟨o', ho', -, hs⟩ := (H.ctx.vert_edges e' he').1 ((hVedges e').1 hue)
      exact H.done_nc o' ho' e' hs x hix hx
    · exact (hfreshY x hx).2.2.2.2.2 u hu ⟨e', he', hue, hix⟩
  have hbase_nv : ∀ u ∈ L ++ s.tstack, u.vStart ≠ y := by
    intro u hu
    rcases List.mem_append.1 hu with hu | hu
    · rw [(hmemL u hu).1]; exact hvy
    · exact (hfreshY y hyv).2.2.2.2.1 u hu
  have hsvd : sE.stackVerts[d]! = v := by rw [H.svlo d (Nat.le_refl d)]; exact H.ctx.sv_d
  have hInc : ∀ e' x, s.g.Inc e' x → sE.g.Inc e' x := fun e' x h => by rw [hgE]; exact h
  obtain ⟨top', hts', hTop'⟩ := H.ctx'.top
  have htstack : (pushEnd sE D₃ L').tstack = (L' ++ top') ++ (L ++ s.tstack) := by
    show L' ++ sE.tstack = _; rw [hts', List.append_assoc]
  have hsubE : ∀ t ∈ L' ++ top', ∀ e', e' < sE.g.ne → t.edges sE.g sE.items e' → subEdges o e' := by
    intro t ht e' he' hte
    refine (hsub e').2 (.inr (hsubChild e' ?_))
    rcases List.mem_append.1 ht with ht | ht
    · obtain ⟨rfl, -⟩ := hmemL' t ht
      obtain ⟨x, hx, -, hs⟩ := (H.ctx'.vert_edges e' he').1 ((hVyedges e').1 hte)
      exact ⟨x, hx, hs⟩
    · exact hTop'.edges t ht e' he' hte
  refine ⟨L' ++ top', L ++ s.tstack, rfl, ?_⟩
  refine
    { tstack := htstack
      back_nil := fun h => (hT h).elim
      base_bot := fun _ t ht => hbase_nv t ht
      sub_bot := ?sub_bot
      path := fun k k' hk hk' => H.ctx'.path k k' hk (Nat.le_succ_of_le hk')
      disj := ?disj
      span_disj := ?span_disj
      sub_edges := hsubE
      base_disj := ?base_disj
      sub_cover := ?sub_cover
      loop1_side := earAt_tree_loop1_side H (L' ++ top') htstack
      loop1_touch := earAt_tree_loop1_touch H (L' ++ top') htstack
      touch_bot := ?touch_bot
      vert := ?vert
      vert_free := ?vert_free
      q_free := ?q_free
      q_root := fun p => by
        show ¬ Items.IsParent sE.items p (edgeItem sE.g e); rw [hgE]; exact H.q_root p
      v_root := H.v_root
      p_entry := ?p_entry
      loop1 := earAt_tree_loop1 H (L' ++ top') htstack
      bottom := earAt_tree_bottom H (L' ++ top') htstack
      loops := earAt_tree_loops H (L' ++ top') htstack
      late := earAt_tree_late H (L' ++ top') htstack
      late_fo := earAt_tree_late_fo H (L' ++ top') htstack
      close := earAt_tree_close H (L' ++ top') htstack
      sv_d := hsvd
      sv_child := fun _ => H.ctx'.sv_d
      path_child := fun _ k hk => by
        show sE.stackVerts[k]! ≠ y
        rw [← H.ctx'.sv_d]; exact H.ctx'.path k (d + 1) (by omega) (Nat.le_refl _)
      dir_d := ?dir_d
      boundary := ?boundary
      base_touch := ?base_touch
      bd_noVert := ?bd_noVert
      bd_bridge := earAt_tree_bd_bridge H (L' ++ top') htstack
      bd_comp := earAt_tree_bd_comp H (L' ++ top') htstack
      bd_term := earAt_tree_bd_term H (L' ++ top') htstack
      bd_side := earAt_tree_bd_side H (L' ++ top') htstack
      lower := earAt_tree_lower H (L' ++ top') htstack }
  case sub_bot =>
    intro t ht
    rcases List.mem_append.1 ht with ht | ht
    · rw [(hmemL' t ht).1]; exact fun h => hvy h.symm
    · have := hTop'.bot t ht d (Nat.lt_succ_self d); rwa [hsvd] at this
  case disj =>
    refine List.pairwise_append.2 ⟨hLpw' _, H.ctx'.disj, ?_⟩
    intro t ht t' ht' e' he' hte
    obtain ⟨rfl, hp⟩ := hmemL' t ht
    exact fun h => H.ctx'.vert_disj (hpushV' hp) t' ht' e' he' h ((hVyedges e').1 hte)
  case span_disj =>
    refine List.pairwise_append.2 ⟨hLpw' _, H.ctx'.span_disj, ?_⟩
    intro t ht t' ht' i hi hi'
    obtain ⟨rfl, hp⟩ := hmemL' t ht
    rw [(hVyspan i).1 hi] at hi'
    exact absurd (H.ctx'.vert_free t' ht' hi') (by rw [hpushV' hp]; decide)
  case base_disj =>
    intro t ht e' he' hte hs
    have he's : e' < s.g.ne := by rw [← hne]; exact he'
    have hte' := (hbaseE t ht e' he's).1 hte
    rcases List.mem_append.1 ht with ht | ht
    · obtain ⟨rfl, -⟩ := hmemL t ht
      obtain ⟨o', ho', -, hs'⟩ := (H.ctx.vert_edges e' he's).1 ((hVedges e').1 hte')
      exact H.nd o' ho' e' hs' hs
    · have hq := H.ctx.q_fresh o List.mem_cons_self e' hs
      exact TEntry.edges_of_root_not_mem hq.1 (hq.2.2 t ht) hte'
  case sub_cover =>
    intro e' he' hne' hs
    rcases (hsub e').1 hs with h | h
    · exact (hne' h).elim
    obtain ⟨o', ho', hs'⟩ := mem_subEdges_edgesList.1 h
    rw [← H.hdone'] at ho'
    obtain ⟨x, hx, rfl⟩ := List.mem_map.1 ho'
    by_cases hlt : x.1.cls.lowval (d + 1) < d + 1
    · obtain ⟨t, ht, hte⟩ := hTop'.cover x hx hlt e' he' hs'
      exact ⟨t, List.mem_append_right _ ht, hte⟩
    · have hEB := (H.ctx'.vert_edges e' he').2 ⟨x, hx, Nat.le_of_not_lt hlt, hs'⟩
      cases hhv : hv'
      · have hp : push' = true := H.hpush'.2 hhv
        exact ⟨Vy, List.mem_append_left _ (by rw [H.hL', hp]; simp [hVy]), (hVyedges e').2 hEB⟩
      · obtain ⟨above, below, hsplit, -⟩ := hTop'.split
        rw [hhv] at hsplit
        obtain ⟨vt, htop, -, hmem, -⟩ := hsplit
        exact ⟨vt, List.mem_append_right _ (htop ▸ List.mem_append_right _ List.mem_cons_self),
          vertItem y, hmem, hEB⟩
  case touch_bot =>
    intro t ht ⟨e', he', hte⟩
    rcases List.mem_append.1 ht with ht | ht
    · obtain ⟨rfl, -⟩ := hmemL' t ht
      obtain ⟨x, hx, hge, -⟩ := (H.ctx'.vert_edges e' he').1 ((hVyedges e').1 hte)
      have hi := H.inc' x hx
      have hlt : x.1.e < sE.g.ne := by rw [hne]; exact hi.1
      exact ⟨x.1.e, hlt, (hVyedges _).2 ((H.ctx'.vert_edges _ hlt).2 ⟨x, hx, hge, .inl rfl⟩),
        hInc _ _ hi.2⟩
    · exact H.ctx'.touch_bot t ht ⟨e', he', hte⟩
  case vert =>
    intro h
    cases hp : push
    · rw [hp, Bool.or_false] at h
      obtain ⟨top, hts, hTop⟩ := H.ctx.top
      obtain ⟨above, below, hsplit, -⟩ := hTop.split
      rw [h] at hsplit
      obtain ⟨vt, htop, hle, hmem, -⟩ := hsplit
      exact ⟨vt, List.mem_append_right _ (hts ▸ htop ▸ List.mem_append_left _
        (List.mem_append_right _ List.mem_cons_self)), hle, hmem⟩
    · exact ⟨V, List.mem_append_left _ (by rw [H.hL, hp]; simp [hV, hb_def]), Nat.le_refl _,
        (hVspan _).2 rfl⟩
  case vert_free =>
    intro t ht hi
    rcases List.mem_append.1 ht with ht | ht
    · obtain ⟨rfl, -⟩ := hmemL' t ht
      exact absurd (vertItem_inj' ((hVyspan _).1 hi)) hvy
    rw [hts'] at ht
    rcases List.mem_append.1 ht with ht | ht
    · rcases hTop'.vitems t ht v (by rw [hnv]; exact H.hv) hi with h | h
      · exact absurd h hvy
      · rw [H.hdone'] at h
        exact absurd (List.mem_cons_of_mem _ h) H.v_nc
    · refine ⟨?_, ht⟩
      rcases List.mem_append.1 ht with ht | ht
      · rw [(hmemL t ht).2]; simp
      · rw [H.ctx.vert_free t ht hi]; simp
  case q_free =>
    intro t ht hi
    have hiE : edgeItem sE.g e ∈ t.spans.1 ++ t.spans.2 := hi
    have hi : edgeItem s.g e ∈ t.spans.1 ++ t.spans.2 := by rw [← hgE]; exact hiE
    rcases List.mem_append.1 ht with ht | ht
    · obtain ⟨rfl, -⟩ := hmemL' t ht
      exact vertItem_ne_edgeItem' H.hy ((hVyspan _).1 hi).symm
    rw [hts'] at ht
    rcases List.mem_append.1 ht with ht | ht
    · exact H.e_ne (hsubChild e (hTop'.qitems t ht e (by rw [hne]; exact H.e_lt) hiE))
    · rcases List.mem_append.1 ht with ht | ht
      · obtain ⟨rfl, -⟩ := hmemL t ht
        exact vertItem_ne_edgeItem' H.hv ((hVspan _).1 hi).symm
      · exact hqo.2.2 t ht hi
  case p_entry =>
    intro hlt ht1 t ht hvs htd
    simp only [ho, DfsOut.cls] at hlt ht1 htd
    have hne_d : t.topDepth ≠ d := by omega
    have hS : CtxSingle v s t := by
      rcases List.mem_append.1 ht with ht | ht
      · obtain ⟨rfl, -⟩ := hmemL t ht
        exact absurd rfl hne_d
      obtain ⟨top, hts, hTop⟩ := H.ctx.top
      rw [hts] at ht
      rcases List.mem_append.1 ht with ht | ht
      · obtain ⟨above, below, hsplit, -, hsingle, -, -, -, hbelow⟩ := hTop.split
        have hA : t ∈ above := by
          cases hasVert
          · obtain ⟨rfl, rfl⟩ := hsplit; exact absurd hvs (hbelow t ht)
          · obtain ⟨vt, htop, -, -, hvt, -⟩ := hsplit
            rw [htop] at ht
            rcases List.mem_append.1 ht with ht | ht
            · exact ht
            rcases List.mem_cons.1 ht with rfl | ht
            · exact absurd (hvt hvs).1 hne_d
            · exact absurd hvs (hbelow t ht)
        exact hsingle t hA (htd ▸ allType1_of_rank hlt ht1 H.rank H.ctx.afterVert_ret)
      · exact absurd hvs (H.ctx.base_bot t ht)
    obtain ⟨i, hspans, -⟩ := hS.item
    have htE : t ∈ sE.tstack := hts' ▸ List.mem_append_right _ ht
    refine ⟨⟨i, ?_, fun p => H.ctx'.span_root t htE i (by rw [hspans]; exact setSides_mem_single _ _) p⟩, ?_⟩
    · show t.spans = setSides D₃[t.topDepth]! [i] []
      rw [H.sd₃ _ (by omega), H.sdlo _ (by omega), hsdk _ hne_d]; exact hspans
    · intro x hx
      have hx' : sE.g.Touches (t.edges sE.g sE.items) x := hx
      rw [hgE] at hx'
      rcases hS.att x ((Graph.Touches.congr' (hbaseE' t ht)).1 hx') with h | h | h
      · exact .inl h
      · refine .inr (.inl ?_)
        rw [h]
        show s.stackVerts[t.topDepth]! = sE.stackVerts[t.topDepth]!
        rw [H.svlo _ (by omega)]
      · refine .inr (.inr ?_)
        show sE.g.Interior (t.edges sE.g sE.items) x
        rw [hgE]
        exact (Graph.Interior.congr fun e' he' => (hbaseE' t ht e' he').symm).1 h
  case dir_d =>
    intro hlt
    simp only [ho, DfsOut.cls] at hlt
    show D₃[d]! = !D₃[cls.lowval d]!
    rw [H.sd₃ d (Nat.le_refl d), H.sd₃ _ (Nat.le_of_lt hlt), H.sdlo d (Nat.le_refl d),
      H.sdlo _ (Nat.le_of_lt hlt), hsdd, hsdk _ (by omega), hb_def]
    simp [Nat.not_le.2 hlt]
  case boundary =>
    intro hge t ht u hu x htx hux
    simp only [ho, DfsOut.cls] at hge
    obtain ⟨e', he', hte, hix⟩ := htx
    have he's : e' < s.g.ne := by rw [← hne]; exact he'
    have hix' : s.g.Inc e' x := by rw [← hgE]; exact hix
    rcases H.ends hge e' (hsubE t ht e' he' hte) he's x hix' with rfl | hx
    · rfl
    · exfalso
      have hux' : sE.g.Touches (u.edges sE.g sE.items) x := hux
      rw [hgE] at hux'
      exact hbaseNT u hu x hx ((Graph.Touches.congr' (hbaseE' u hu)).1 hux')
  case base_touch =>
    intro _ u hu hux
    have hux' : sE.g.Touches (u.edges sE.g sE.items) y := hux
    rw [hgE] at hux'
    exact hbaseNT u hu y hyv ((Graph.Touches.congr' (hbaseE' u hu)).1 hux')
  case bd_noVert =>
    intro hge
    simp only [ho, DfsOut.cls] at hge
    have hp : push = false := by
      cases hp : push
      · rfl
      · exact absurd (hpushV hp).2.1 (by omega)
    rw [hp, Bool.or_false]
    cases hhv : hasVert
    · rfl
    · exact absurd (H.ctx.hv_ret hhv)
        (no_ret_before_boundary hge H.cls_ret H.rank)

/-- The frame of loops 1–2 at a returning tree edge without the vertex entry, over the `EarClose`
stack `R`: the `RetTop` fields `EarClose` does not supply (`retCheck`: `sd`, `sdd`, `size`, `nxt`,
`keep_ch`/`keep_par`, `ch_lt`, `span_lt`/`span_root`/`span_old`/`span_owned`, `touch_bot`, `touch`,
`vstart`, `noVert_v`/`noVert_path`). -/
structure RetFrame (v d : Nat) (s : WalkState) (e : Nat) (cls : OutClass) (y : Nat)
    (outs : List DfsOut) (L : List TEntry) (sX : WalkState) (R : List TEntry) : Prop where
  sd : ∀ k, k < d → sX.stackDir[k]! = s.stackDir[k]!
  sdd : sX.stackDir[d]! = !s.stackDir[cls.lowval d]!
  size : s.items.size ≤ sX.items.size
  nxt : s.nxtEdgeIdx ≤ sX.nxtEdgeIdx
  kept : ∀ j, j < s.items.size → (∀ x ∈ (DfsTree.node y outs).verts, j ≠ vertItem x) →
    (∀ e' ∈ (DfsTree.node y outs).edges, j ≠ edgeItem s.g e') → j ≠ edgeItem s.g e →
    Items.ch sX.items j = Items.ch s.items j ∧
    ∀ p, Items.IsParent sX.items p j ↔ Items.IsParent s.items p j
  ch_lt : ∀ p c, Items.IsParent sX.items p c → c < sX.items.size
  span_lt : ∀ t ∈ R, ∀ i ∈ t.spans.1 ++ t.spans.2, i < sX.items.size
  span_root : ∀ t ∈ R, ∀ i ∈ t.spans.1 ++ t.spans.2, ∀ p, ¬ Items.IsParent sX.items p i
  span_new : ∀ t ∈ R, ∀ i ∈ t.spans.1 ++ t.spans.2,
    (∀ t' ∈ L ++ s.tstack, i ∉ t'.spans.1 ++ t'.spans.2) ∧
    (i < s.items.size → (∃ x ∈ (DfsTree.node y outs).verts, i = vertItem x) ∨
      (∃ e' ∈ (DfsTree.node y outs).edges, i = edgeItem s.g e') ∨ i = edgeItem s.g e)
  touch_bot : ∀ t ∈ R, (∃ e', e' < s.g.ne ∧ t.edges s.g sX.items e') →
    s.g.Touches (t.edges s.g sX.items) t.vStart
  touch : ∀ t ∈ R, ∀ x, s.g.Touches (t.edges s.g sX.items) x →
    x = v ∨ x ∈ (DfsTree.node y outs).verts ∨ ∃ k, cls.lowval d ≤ k ∧ k ≤ d ∧ x = s.stackVerts[k]!
  vstart : ∀ t ∈ R, t.vStart = v ∨ t.vStart ∈ (DfsTree.node y outs).verts
  touch_k : ∀ t ∈ R, ∀ k, k ≤ d → s.g.Touches (t.edges s.g sX.items) s.stackVerts[k]! →
    t.topDepth ≤ k

/-- Admitted (dump-checked, `retCheck`): the frame of loops 1–2 (`RetFrame`) over the `EarClose`
stack of the site's `EarFinish`; without the vertex entry its entries start in the child
(`noVert_*`); with it, `vy` was opened during the child's walk (`x_first_pre`), the entries reach
`stackVerts[lowval]` (`x_top_pre`) and lie at depths `≥ lowval` (`x_mid_top`). -/
theorem tree_ret_frame (hr : cls.lowval d < d) (hv : Bool) (c : TEntry) (mid : List TEntry)
    (py vy : TEntry)
    (hC : EarClose v d (cls.lowval d) o₀ hv (L ++ s.tstack) s₃ (feS₂ d o₀ s₃) c mid py vy) :
    RetFrame v d s e cls y outs L (feS₂ d o₀ s₃) (c :: mid ++ [py, vy]) ∧
    (hv = false → ∀ t ∈ c :: mid ++ [py, vy], t.vStart ≠ v ∧
      ∀ k, k ≤ d → s.stackVerts[k]! ≠ t.vStart) ∧
    (hv = true → (s.nxtEdgeIdx ≤ vy.firstIdx ∧ vy.firstIdx < (feS₂ d o₀ s₃).nxtEdgeIdx) ∧
      s.g.Touches (fun e' => ∃ t ∈ c :: mid ++ [py, vy], t.edges s.g (feS₂ d o₀ s₃).items e')
        s.stackVerts[cls.lowval d]! ∧
      ∀ t ∈ c :: mid ++ [py, vy], cls.lowval d ≤ t.topDepth) := by
  sorry

/-- `tree_ret_vert` from the item well-formedness (`Shape`, `Inv'`) of the state loop 3 runs from:
the vertex close is `closeVert'_type2_eq`/`closeVert'_type1_eq` over the `EarClose` stack of the
site's `EarFinish`, whose frame is `tree_ret_frame`. -/
theorem tree_ret_vert_of_wf (hr : cls.lowval d < d) (hb : (hasVert || push) = true)
    (hs₂ : Shape (feS₂ d o₀ s₃)) (hi₂ : (feS₂ d o₀ s₃).Inv' (d + 1)) :
    ∃ R, RetTop v d s e cls y outs L true (feS₃ v d o₀ (L ++ s.tstack).length s₃) R := by
  obtain ⟨sub, base₀, hlen, hE⟩ := earAt_tree_of_ctx H
  have hbase : base₀ = L ++ s.tstack := by
    obtain ⟨top', hts', -⟩ := H.ctx'.top
    have h1 : (s₃).tstack = (L' ++ top') ++ (L ++ s.tstack) := by
      show L' ++ sE.tstack = _; rw [hts', List.append_assoc]
    exact List.append_inj_right' (hE.tstack.symm.trans h1) hlen
  subst hbase
  rw [hb] at hE
  have hg3 : (s₃).g = s.g := H.gE
  obtain ⟨c, mid, py, vy, hC⟩ := hE.close H.tree hr
  obtain ⟨F, -, hX⟩ := tree_ret_frame H hr true c mid py vy hC
  obtain ⟨⟨hf1, hf2⟩, htop, hmidtop⟩ := hX rfl
  have hgs : (feS₂ d o₀ s₃).g = s.g := hC.g.trans hg3
  have hdir : (s₃).stackDir[d]! = !(s₃).stackDir[cls.lowval d]! := hE.dir_d hr
  have hsl : (s₃).stackDir[cls.lowval d]! = s.stackDir[cls.lowval d]! :=
    hC.dir_l.symm.trans (F.sd _ hr)
  have hsvl : (s₃).stackVerts[cls.lowval d]! = s.stackVerts[cls.lowval d]! := H.svlo _ (Nat.le_of_lt hr)
  have hchild : ∀ k, k ≤ d → (s₃).stackVerts[k]! ≠ (s₃).stackVerts[d + 1]! := fun k hk => by
    rw [hE.sv_child H.tree]; exact hE.path_child H.tree k hk
  have hok := closeVertOk_of_close (feSingle d o₀ s₃) hC hs₂ hr hE.path hchild hdir
  have hsh : Shape (feS₃ v d o₀ (L ++ s.tstack).length s₃) :=
    (Step.closeVert' hi₂ hs₂ (v := v) (by rw [hgs]; exact H.hv) hok).shape
  have hsubE : ∀ t ∈ c :: mid ++ [py, vy], ∀ e', e' < s.g.ne →
      t.edges s.g (feS₂ d o₀ s₃).items e' → subEdges o₀ e' := by rw [← hg3]; exact hC.sub_edges
  have hcover : ∀ e', e' < s.g.ne → subEdges o₀ e' →
      ∃ t ∈ c :: mid ++ [py, vy], t.edges s.g (feS₂ d o₀ s₃).items e' := by
    rw [← hg3]; exact hC.sub_cover
  have hce : c.edges s.g (feS₂ d o₀ s₃).items e := by rw [← hg3]; exact hC.c_edge
  have hinc : s.g.Inc e v := by
    rcases H.e_ends with h | h
    · exact Or.inr (congrArg Prod.snd h).symm
    · exact Or.inl (congrArg Prod.snd h).symm
  have hint : ∀ E : Nat → Prop, (∀ e', e' < s.g.ne → subEdges o₀ e' → E e') →
      ∀ x ∈ (DfsTree.node y outs).verts, s.g.Interior E x :=
    fun E hE' x hx e' he' hinc' => hE' e' he' (H.comp e' he' x hinc' hx)
  have hatt : ∀ E : Nat → Prop, (∀ e', e' < s.g.ne → subEdges o₀ e' → E e') →
      (∀ x, s.g.Touches E x → x = v ∨ x ∈ (DfsTree.node y outs).verts ∨
        ∃ k, cls.lowval d ≤ k ∧ k ≤ d ∧ x = s.stackVerts[k]!) →
      ∀ w, s.g.Touches E w → w ≠ v → ¬ s.g.Interior E w →
        ∃ k, cls.lowval d ≤ k ∧ k ≤ d ∧ w = s.stackVerts[k]! := fun E hE' hT w hw hwv hwi => by
    rcases hT w hw with h | h | h
    · exact absurd h hwv
    · exact absurd (hint E hE' w h) hwi
    · exact h
  cases ht1 : cls.isType1
  · have h3 : feS₃ v d o₀ (L ++ s.tstack).length s₃ =
        { feS₂ d o₀ s₃ with
          tstack := { l2Cur c (mid ++ [py, vy]) with
            vStart := v
            spans := setSides (!(s₃).stackDir[d]!)
              ((l2Cur c (mid ++ [py, vy])).spans.1 ++ (l2Cur c (mid ++ [py, vy])).spans.2) [] } ::
            (L ++ s.tstack) } := by
      show after (closeVert' v (s₃).stackDir[d]! cls.isType1 _ _) (feS₂ d o₀ s₃) = _
      rw [ht1]; exact closeVert'_type2_eq _ _ hC
    set m := l2Cur c (mid ++ [py, vy]) with hm
    set r : TEntry := { m with
      vStart := v
      spans := setSides (!(s₃).stackDir[d]!) (m.spans.1 ++ m.spans.2) [] } with hr_def
    rw [h3] at hsh ⊢
    have hrsp : ∀ i, i ∈ r.spans.1 ++ r.spans.2 ↔
        ∃ t ∈ c :: mid ++ [py, vy], i ∈ t.spans.1 ++ t.spans.2 := fun i => by
      show i ∈ (setSides _ _ []).1 ++ (setSides _ _ []).2 ↔ _
      rw [mem_setSides, hm, mem_l2Cur_spans, List.cons_append]
      simp only [List.mem_cons, exists_eq_or_imp]
    have hrE : ∀ e', r.edges s.g (feS₂ d o₀ s₃).items e' ↔
        ∃ t ∈ c :: mid ++ [py, vy], t.edges s.g (feS₂ d o₀ s₃).items e' := fun e' => by
      simp only [TEntry.edges, hrsp]
      constructor
      · rintro ⟨i, ⟨t, ht, hi⟩, hb⟩; exact ⟨t, ht, i, hi, hb⟩
      · rintro ⟨t, ht, i, hi, hb⟩; exact ⟨i, ⟨t, ht, hi⟩, hb⟩
    have hrtop : r.topDepth = cls.lowval d := by
      show m.topDepth = _
      apply Nat.le_antisymm
      · have hpt : py.topDepth = cls.lowval d := hC.py_top
        rw [← hpt]; exact (l2Cur_top_le c _).2 py (by simp)
      · exact le_l2Cur_top c _ _ (hmidtop c (by simp)) fun t ht => hmidtop t (by simp [ht])
    have hrfirst : r.firstIdx = vy.firstIdx := by
      show (l2Cur c (mid ++ [py, vy])).firstIdx = _
      rw [show mid ++ [py, vy] = (mid ++ [py]) ++ [vy] by simp, l2Cur_concat]; rfl
    have htouchr : ∀ t ∈ [r], ∀ x, s.g.Touches (t.edges s.g (feS₂ d o₀ s₃).items) x →
        x = v ∨ x ∈ (DfsTree.node y outs).verts ∨
          ∃ k, cls.lowval d ≤ k ∧ k ≤ d ∧ x = s.stackVerts[k]! := fun t ht x hx => by
      rw [List.mem_singleton] at ht; rw [ht] at hx
      obtain ⟨e', he', hre, hinc'⟩ := hx
      obtain ⟨t', ht', hte⟩ := (hrE e').1 hre
      exact F.touch t' ht' x ⟨e', he', hte, hinc'⟩
    refine ⟨[r], {
      tstack := rfl
      g := hgs
      sv := fun k hk => by show (feS₂ d o₀ s₃).stackVerts[k]! = _; rw [hC.sv]; exact H.svlo k hk
      sd := F.sd
      sdd := fun h => by cases h
      size := F.size
      nxt := F.nxt
      kept := F.kept
      ch_lt := hsh.ch_lt
      span_lt := fun t ht => hsh.span t (by rw [List.mem_singleton] at ht; rw [ht]; exact List.mem_cons_self ..)
      span_root := fun t ht i hi => by
        rw [List.mem_singleton] at ht; rw [ht, hrsp] at hi
        obtain ⟨t', ht', hi'⟩ := hi
        exact F.span_root t' ht' i hi'
      span_new := fun t ht i hi => by
        rw [List.mem_singleton] at ht; rw [ht, hrsp] at hi
        obtain ⟨t', ht', hi'⟩ := hi
        exact F.span_new t' ht' i hi'
      touch_bot := fun t ht _ => by
        rw [List.mem_singleton] at ht; rw [ht]
        exact ⟨e, H.e_lt, (hrE e).2 ⟨c, by simp, hce⟩, hinc⟩
      edges := fun e' he' => by
        simp only [List.mem_singleton, exists_eq_left]
        rw [hrE e']
        exact ⟨fun ⟨t, ht, hte⟩ => hsubE t ht e' he' hte, hcover e' he'⟩
      touch := htouchr
      touch_k := fun t ht k hk hx => by
        rw [List.mem_singleton] at ht; rw [ht] at hx ⊢; rw [hrtop]
        obtain ⟨e', he', hre, hinc'⟩ := hx
        obtain ⟨t', ht', hte⟩ := (hrE e').1 hre
        exact Nat.le_trans ((hX rfl).2.2 t' ht') (F.touch_k t' ht' k hk ⟨e', he', hte, hinc'⟩)
      disj := List.pairwise_singleton ..
      span_disj := List.pairwise_singleton ..
      vstart := fun t ht => by rw [List.mem_singleton] at ht; rw [ht]; exact Or.inl rfl
      vert := fun _ => ⟨r, rfl, rfl, hrtop, ?_, by rw [hrfirst]; exact hf1,
        by show r.firstIdx < (feS₂ d o₀ s₃).nxtEdgeIdx; rw [hrfirst]; exact hf2, ?_,
        hatt _ (fun e' he' hs => (hrE e').2 (hcover e' he' hs)) (htouchr r (List.mem_singleton_self _)),
        fun h => by rw [ht1] at h; cases h⟩
      noVert := fun h => by cases h }⟩
    · show getSide (setSides (!(s₃).stackDir[d]!) _ []) _ = []
      rw [← hsl, hdir, Bool.not_not]
      exact getSide_setSides_not _ _
    · obtain ⟨e', he', ⟨t, ht, hte⟩, hinc'⟩ := htop
      exact ⟨e', he', (hrE e').2 ⟨t, ht, hte⟩, hinc'⟩
  · obtain ⟨item, i, s₂, py', f, h3, hfch, hsdr, hmtop, hts₂, hg₂, hsv₂, hsd₂, hnxt₂, hsize₂, hch₂,
      hfree, hpysp, hiroot, hwhich, hpysp', hallE, hpE', hr'E'⟩ :=
      closeVert'_type1_eq (feSingle d o₀ s₃) hC hi₂ hs₂ ht1 hr hE.path hchild hdir
    have hmid : mid = [] := (hC.type1 ht1).1
    have htouch1 := (hC.type1 ht1).2
    subst hmid
    have h3' : feS₃ v d o₀ (L ++ s.tstack).length s₃ =
        after (closeVert' v (s₃).stackDir[d]! true (L ++ s.tstack).length (feSingle d o₀ s₃))
          (feS₂ d o₀ s₃) := by
      show after (closeVert' v (s₃).stackDir[d]! cls.isType1 _ _) (feS₂ d o₀ s₃) = _
      rw [ht1]
    rw [← h3'] at h3 hr'E'
    set m := TEntry.mergeInto (TEntry.mergeInto c py') vy with hm
    set r : TEntry := { m with
      vStart := v
      spans := setSides s₂.stackDir[m.topDepth]! [item] [] } with hr_def
    set sX := feS₃ v d o₀ (L ++ s.tstack).length s₃ with hsX
    have hXi : sX.items = s₂.items.modify item f := by rw [h3]
    have hXt : sX.tstack = r :: (L ++ s.tstack) := by rw [h3]
    have hXg : sX.g = s.g := by rw [h3]; exact hg₂.trans hgs
    have hXsv : sX.stackVerts = (feS₂ d o₀ s₃).stackVerts := by rw [h3]; exact hsv₂
    have hXsd : sX.stackDir = (feS₂ d o₀ s₃).stackDir := by rw [h3]; exact hsd₂
    have hXn : sX.nxtEdgeIdx = (feS₂ d o₀ s₃).nxtEdgeIdx := by rw [h3]; exact hnxt₂
    have hequiv : ∀ e', e' < s.g.ne → (r.edges s.g sX.items e' ↔
        c.edges s.g (feS₂ d o₀ s₃).items e' ∨ py.edges s.g (feS₂ d o₀ s₃).items e' ∨
          vy.edges s.g (feS₂ d o₀ s₃).items e') :=
      fun e' he' => by rw [← hgs] at he' ⊢; exact hr'E' e' he'
    have hrE : ∀ e', e' < s.g.ne → (r.edges s.g sX.items e' ↔
        ∃ t ∈ c :: [] ++ [py, vy], t.edges s.g (feS₂ d o₀ s₃).items e') := fun e' he' => by
      rw [hequiv e' he']; simp
    have htouch1' : ∀ w, s.g.Touches (fun e' => c.edges s.g (feS₂ d o₀ s₃).items e' ∨
        py.edges s.g (feS₂ d o₀ s₃).items e' ∨ vy.edges s.g (feS₂ d o₀ s₃).items e') w →
        w = v ∨ w = s.stackVerts[cls.lowval d]! ∨
          s.g.Interior (fun e' => c.edges s.g (feS₂ d o₀ s₃).items e' ∨
            py.edges s.g (feS₂ d o₀ s₃).items e' ∨ vy.edges s.g (feS₂ d o₀ s₃).items e') w := by
      rw [← hg3, ← hsvl]; exact htouch1
    have hsdr' : s₂.stackDir[m.topDepth]! = (s₃).stackDir[cls.lowval d]! := hsdr
    have hrsp : ∀ i', i' ∈ r.spans.1 ++ r.spans.2 ↔ i' = item := fun i' => by
      show i' ∈ (setSides _ _ []).1 ++ (setSides _ _ []).2 ↔ _
      rw [mem_setSides, List.mem_singleton]
    have hipy : i ∈ py.spans.1 ++ py.spans.2 := by
      rw [hpysp]; exact (mem_setSides _ _ _).2 (List.mem_singleton_self _)
    have hich : ∀ j, ¬ Items.IsParent s.items i j := fun j hij => by
      have hilt := isParent_lt_size hij
      rcases (F.span_new py (by simp) i hipy).2 hilt with ⟨x, hx, hi⟩ | ⟨e', he', hi⟩ | hi
      · have := (H.ctx.v_fresh o₀ (List.mem_cons_self ..) e cls _ rfl x hx).2.1
        rw [Items.IsParent, hi, this] at hij; simp at hij
      · have := (H.ctx.q_fresh o₀ (List.mem_cons_self ..) e' (Or.inr he')).2.1
        rw [Items.IsParent, hi, this] at hij; simp at hij
      · have := (H.ctx.q_fresh o₀ (List.mem_cons_self ..) e (Or.inl rfl)).2.1
        rw [Items.IsParent, hi, this] at hij; simp at hij
    have hnew : ∀ j, j < s.items.size → (∀ x ∈ (DfsTree.node y outs).verts, j ≠ vertItem x) →
        (∀ e' ∈ (DfsTree.node y outs).edges, j ≠ edgeItem s.g e') → j ≠ edgeItem s.g e →
        ∀ t ∈ c :: [] ++ [py, vy], j ∉ t.spans.1 ++ t.spans.2 := fun j hj hjv hje hjq t ht hjm => by
      rcases (F.span_new t ht j hjm).2 hj with ⟨x, hx, h⟩ | ⟨e', he', h⟩ | h
      · exact hjv x hx h
      · exact hje e' he' h
      · exact hjq h
    have hmsp : ∀ j, j ∈ m.spans.1 ++ m.spans.2 → j ∈ c.spans.1 ++ c.spans.2 ∨
        j ∈ py'.spans.1 ++ py'.spans.2 ∨ j ∈ vy.spans.1 ++ vy.spans.2 := fun j hjm => by
      rcases (TEntry.mem_mergeInto _ _ j).1 hjm with hjm | hjm
      · rcases (TEntry.mem_mergeInto _ _ j).1 hjm with hjm | hjm
        · exact Or.inl hjm
        · exact Or.inr (Or.inl hjm)
      · exact Or.inr (Or.inr hjm)
    have hitem_free : item ∉ m.spans.1 ++ m.spans.2 := fun hm' => by
      rcases hmsp item hm' with h | h | h
      · exact hfree.free c (by rw [hts₂]; simp) h
      · exact hfree.free py' (by rw [hts₂]; simp) h
      · exact hfree.free vy (by rw [hts₂]; simp) h
    have hXroot : ∀ p, ¬ Items.IsParent sX.items p item := fun p hp => by
      rw [hXi] at hp
      rcases Items.IsParent_modify hp with hp | ⟨_, _, hm'⟩
      · exact hfree.root p hp
      · rw [hfch] at hm'; exact hitem_free hm'
    have hitem_new : (∀ t' ∈ L ++ s.tstack, item ∉ t'.spans.1 ++ t'.spans.2) ∧
        (item < s.items.size → (∃ x ∈ (DfsTree.node y outs).verts, item = vertItem x) ∨
          (∃ e' ∈ (DfsTree.node y outs).edges, item = edgeItem s.g e') ∨ item = edgeItem s.g e) := by
      refine ⟨fun t' ht' => hfree.free t' (by rw [hts₂]; simp [ht']), fun hlt => ?_⟩
      rcases hwhich with h | h
      · rw [h] at hlt ⊢; exact (F.span_new py (by simp) i hipy).2 hlt
      · rw [h] at hlt; exact absurd hlt (Nat.not_lt.2 F.size)
    have hkept : ∀ j, j < s.items.size → (∀ x ∈ (DfsTree.node y outs).verts, j ≠ vertItem x) →
        (∀ e' ∈ (DfsTree.node y outs).edges, j ≠ edgeItem s.g e') → j ≠ edgeItem s.g e →
        Items.ch sX.items j = Items.ch s.items j ∧
        ∀ p, Items.IsParent sX.items p j ↔ Items.IsParent s.items p j := fun j hj hjv hje hjq => by
      have hk := F.kept j hj hjv hje hjq
      have hn := hnew j hj hjv hje hjq
      have hj' : j ≠ item := by
        rcases hwhich with h | h
        · rw [h]; intro hji; exact hn py (by simp) (by rw [hji]; exact hipy)
        · rw [h]; exact Nat.ne_of_lt (Nat.lt_of_lt_of_le hj F.size)
      rw [hXi]
      refine ⟨(Items.ch_modify_of_ne item f hj').trans ((hch₂ j).trans hk.1),
        fun p => ⟨fun hp => ?_, fun hp => ?_⟩⟩
      · rcases Items.IsParent_modify hp with hp | ⟨_, _, hjm⟩
        · exact (hk.2 p).1 ((Items.IsParent_congr hch₂).1 hp)
        · exfalso
          rw [hfch] at hjm
          rcases hmsp j hjm with h | h | h
          · exact hn c (by simp) h
          · rcases hpysp' j h with h | h
            · exact hn py (by simp) h
            · exact hich j ((hk.2 i).1 h)
          · exact hn vy (by simp) h
      · have hp' := (Items.IsParent_congr hch₂).2 ((hk.2 p).2 hp)
        have hpne : p ≠ item := fun hpi => by
          rw [hpi] at hp
          rcases hwhich with h | h
          · rw [h] at hp; exact hich j hp
          · rw [h] at hp; exact absurd (isParent_lt_size hp) (Nat.not_lt.2 F.size)
        rw [Items.IsParent, Items.ch_modify_of_ne item f hpne]; exact hp'
    have htouchr : ∀ t ∈ [r], ∀ x, s.g.Touches (t.edges s.g sX.items) x →
        x = v ∨ x ∈ (DfsTree.node y outs).verts ∨
          ∃ k, cls.lowval d ≤ k ∧ k ≤ d ∧ x = s.stackVerts[k]! := fun t ht x hx => by
      rw [List.mem_singleton] at ht; rw [ht] at hx
      obtain ⟨e', he', hre, hinc'⟩ := hx
      obtain ⟨t', ht', hte⟩ := (hrE e' he').1 hre
      exact F.touch t' ht' x ⟨e', he', hte, hinc'⟩
    refine ⟨[r], {
      tstack := hXt
      g := hXg
      sv := fun k hk => by rw [hXsv, hC.sv]; exact H.svlo k hk
      sd := fun k hk => by rw [hXsd]; exact F.sd k hk
      sdd := fun h => by cases h
      size := by rw [hXi, Array.size_modify]; exact Nat.le_trans F.size hsize₂
      nxt := by rw [hXn]; exact F.nxt
      kept := hkept
      ch_lt := hsh.ch_lt
      span_lt := fun t ht => hsh.span t (by rw [hXt, List.mem_singleton] at *; rw [ht]; exact List.mem_cons_self ..)
      span_root := fun t ht i' hi' => by
        rw [List.mem_singleton] at ht; rw [ht, hrsp] at hi'; rw [hi']; exact hXroot
      span_new := fun t ht i' hi' => by
        rw [List.mem_singleton] at ht; rw [ht, hrsp] at hi'; rw [hi']; exact hitem_new
      touch_bot := fun t ht _ => by
        rw [List.mem_singleton] at ht; rw [ht]
        exact ⟨e, H.e_lt, (hrE e H.e_lt).2 ⟨c, by simp, hce⟩, hinc⟩
      edges := fun e' he' => by
        simp only [List.mem_singleton, exists_eq_left]
        rw [hrE e' he']
        exact ⟨fun ⟨t, ht, hte⟩ => hsubE t ht e' he' hte, hcover e' he'⟩
      touch := htouchr
      touch_k := fun t ht k hk hx => by
        rw [List.mem_singleton] at ht; rw [ht] at hx ⊢; rw [hmtop]
        obtain ⟨e', he', hre, hinc'⟩ := hx
        obtain ⟨t', ht', hte⟩ := (hrE e' he').1 hre
        exact Nat.le_trans ((hX rfl).2.2 t' ht') (F.touch_k t' ht' k hk ⟨e', he', hte, hinc'⟩)
      disj := List.pairwise_singleton ..
      span_disj := List.pairwise_singleton ..
      vstart := fun t ht => by rw [List.mem_singleton] at ht; rw [ht]; exact Or.inl rfl
      vert := fun _ => ⟨r, rfl, rfl, hmtop, ?_, hf1, by rw [hXn]; exact hf2, ?_,
        hatt _ (fun e' he' hs => (hrE e' he').2 (hcover e' he' hs)) (htouchr r (List.mem_singleton_self _)),
        fun _ => ⟨⟨item, by show setSides _ [item] [] = _; rw [hsdr', hsl]⟩, fun w hw => ?_⟩⟩
      noVert := fun h => by cases h }⟩
    · show getSide (setSides s₂.stackDir[m.topDepth]! [item] []) _ = []
      rw [hsdr', hsl]
      exact getSide_setSides_not _ _
    · obtain ⟨e', he', ⟨t, ht, hte⟩, hinc'⟩ := htop
      exact ⟨e', he', (hrE e' he').2 ⟨t, ht, hte⟩, hinc'⟩
    · obtain ⟨e', he', hre, hinc'⟩ := hw
      rcases htouch1' w ⟨e', he', (hequiv e' he').1 hre, hinc'⟩ with h | h | h
      · exact Or.inl h
      · exact Or.inr (Or.inl h)
      · exact Or.inr (Or.inr fun e'' he'' hi'' => (hequiv e'' he'').2 (h e'' he'' hi''))

/-- Admitted (plumbing, not a dump-checked invariant): item well-formedness (`Shape`, `Inv'`) of
the state the vertex close runs from. `WalkBackbone.WalkInvOut` carries `Inv'`/`Shape` at every
out site and the `Step` lemmas of `finishEdge` carry them through loops 1–2; `TreeSite` does not
expose them. -/
theorem treeSite_wf (hr : cls.lowval d < d) :
    Shape (feS₂ d o₀ s₃) ∧ (feS₂ d o₀ s₃).Inv' (d + 1) := by
  sorry

/-- The stack after the vertex close (loop 3, the two merges, the retarget and the type-1 close)
at a returning tree edge with the vertex entry (`RetTop`, `hv = true`). -/
theorem tree_ret_vert (hr : cls.lowval d < d) (hb : (hasVert || push) = true) :
    ∃ R, RetTop v d s e cls y outs L true (feS₃ v d o₀ (L ++ s.tstack).length s₃) R :=
  tree_ret_vert_of_wf H hr hb (treeSite_wf H hr).1 (treeSite_wf H hr).2

/-- The stack shape after loops 1–3 at a returning tree edge (`RetTop`), at the state `finishRest`
runs from: without the vertex entry it is the `EarClose` stack of the site's `EarFinish`
(`earAt_tree_of_ctx`) with the frame `tree_ret_frame`; with it, `tree_ret_vert`. -/
theorem tree_ret_shape (hr : cls.lowval d < d) : ∃ R,
    RetTop v d s e cls y outs L (hasVert || push)
      (if hasVert || push then feS₃ v d o₀ (L ++ s.tstack).length s₃ else feS₂ d o₀ s₃) R := by
  obtain ⟨sub, base₀, hlen, hE⟩ := earAt_tree_of_ctx H
  have hbase : base₀ = L ++ s.tstack := by
    obtain ⟨top', hts', -⟩ := H.ctx'.top
    have h1 : (s₃).tstack = (L' ++ top') ++ (L ++ s.tstack) := by
      show L' ++ sE.tstack = _; rw [hts', List.append_assoc]
    exact List.append_inj_right' (hE.tstack.symm.trans h1) hlen
  subst hbase
  have hg3 : (s₃).g = s.g := H.gE
  generalize hb : (hasVert || push) = hv₀ at hE ⊢
  cases hv₀
  · simp only [Bool.false_eq_true, ↓reduceIte]
    obtain ⟨c, mid, py, vy, hC⟩ := hE.close H.tree hr
    obtain ⟨F, hnv, -⟩ := tree_ret_frame H hr false c mid py vy hC
    have hsubE : ∀ t ∈ c :: mid ++ [py, vy], ∀ e', e' < s.g.ne →
        t.edges s.g (feS₂ d o₀ s₃).items e' → subEdges o₀ e' := by
      rw [← hg3]; exact hC.sub_edges
    have hcover : ∀ e', e' < s.g.ne → subEdges o₀ e' →
        ∃ t ∈ c :: mid ++ [py, vy], t.edges s.g (feS₂ d o₀ s₃).items e' := by
      rw [← hg3]; exact hC.sub_cover
    have hdisj : (feS₂ d o₀ s₃).tstack.Pairwise fun t t' => ∀ e', e' < s.g.ne →
        t.edges s.g (feS₂ d o₀ s₃).items e' → ¬ t'.edges s.g (feS₂ d o₀ s₃).items e' := by
      rw [← hg3]; exact hC.disj
    have hce : c.edges s.g (feS₂ d o₀ s₃).items e := by rw [← hg3]; exact hC.c_edge
    have hsdisj := hC.span_disj
    rw [hC.tstack] at hdisj hsdisj
    refine ⟨c :: mid ++ [py, vy],
      { tstack := hC.tstack
        g := hC.g.trans hg3
        sv := fun k hk => by rw [hC.sv]; exact H.svlo k hk
        sd := F.sd
        sdd := fun _ => F.sdd
        size := F.size
        nxt := F.nxt
        kept := F.kept
        ch_lt := F.ch_lt
        span_lt := F.span_lt
        span_root := F.span_root
        span_new := F.span_new
        touch_bot := F.touch_bot
        edges := fun e' he' => ⟨fun ⟨t, ht, hte⟩ => hsubE t ht e' he' hte, hcover e' he'⟩
        touch := F.touch
        disj := (List.pairwise_append.1 hdisj).1
        span_disj := (List.pairwise_append.1 hsdisj).1
        vstart := F.vstart
        touch_k := F.touch_k
        vert := fun h => absurd h (by decide)
        noVert := fun _ => ⟨⟨c, mid ++ [py, vy], rfl, hce⟩, hnv rfl⟩ }⟩
  · simp only [↓reduceIte]
    exact tree_ret_vert H hr hb

end

end WalkState
end Spqr
