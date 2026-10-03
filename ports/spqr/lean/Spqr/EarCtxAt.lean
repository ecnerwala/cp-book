import Spqr.EarCtx
import Spqr.WalkInv

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
  obtain ⟨lv, k, hc, _⟩ := ret_of_lowval_lt (hret o ho)
  simp only at hr
  cases hc' : cls <;> simp only [hc', OutClass.lowval] at hlt hl ht1 hr <;> try omega
  rw [hc] at hl hr; simp only [OutClass.rank] at hl hr
  subst hl
  rw [hc]
  cases k <;> cases ‹RetKind› <;> simp_all [RetKind.rank, OutClass.isType1]

/-- Rank-sorted outs: no returning out is finished before a boundary one (`hcls`: a `ret` class
returns below `d`, as `classify` guarantees). -/
theorem no_ret_before_boundary {d : Nat} {done : List (DfsOut × Bool)} {cls : OutClass}
    (hge : d ≤ cls.lowval d) (hcls : ∀ lv k, cls = .ret lv k → lv < d)
    (hrank : ∀ o' ∈ done, o'.1.cls.rank ≤ cls.rank) :
    ¬ ∃ o ∈ done, o.1.cls.lowval d < d := by
  rintro ⟨o, ho, hlt⟩
  have hr := hrank o ho
  obtain ⟨lv, k, hc, hlv⟩ := ret_of_lowval_lt hlt
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
    (hpush : push = true ↔ hasVert = false ∧ cls.lowval d < d)
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
  have hpushV : push = true → hasVert = false ∧ cls.lowval d < d := hpush.1
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
      · exact absurd (hpushV hp).2 (by omega)
    rw [hp, Bool.or_false]
    cases hhv : hasVert
    · rfl
    · exact absurd (hC.hv_ret hhv (List.cons_ne_nil _ _)) (no_ret_before_boundary hge hcls hrank)

end WalkState
end Spqr
