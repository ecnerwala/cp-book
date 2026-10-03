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
  ends : ∀ e', subEdges (.tree e cls (.node y outs)) e' → e' < s.g.ne → ∀ x, s.g.Inc e' x →
    x = v ∨ x ∈ (DfsTree.node y outs).verts
  hpush : push = true ↔ hasVert = false ∧ cls.lowval d < d
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

/-! ### Named admissions (PROOF.md §4.2b): the `EarFinish` fields at a tree-edge site not yet
derived from the child's end-of-outs context. Each is the field verbatim, over the site's `sub`. -/
section
variable (sub : List TEntry) (hsub : (s₃).tstack = sub ++ (L ++ s.tstack))
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
  sorry

theorem earAt_tree_bd_comp : (o₀).cls.isTree = true → d ≤ (o₀).cls.lowval d → (o₀).cls.lowval d ≠ d + 1 →
    ∃ t₁ t₂, sub = [t₁, t₂] ∧ t₁.vStart = (o₀).dest ∧ t₁.topDepth = (o₀).cls.lowval d ∧
      t₂.vStart = (o₀).dest ∧ t₂.topDepth = d + 1 := by
  sorry

theorem earAt_tree_bd_term : (o₀).cls.isTree = true → d ≤ (o₀).cls.lowval d → ∀ u ∈ (s₃).tstack.tail,
    (s₃).g.Touches (u.edges (s₃).g (s₃).items) v → u.vStart = v ∨ u.topDepth ≤ d := by
  sorry

theorem earAt_tree_bd_side : (o₀).cls.isTree = true → d ≤ (o₀).cls.lowval d →
    if (o₀).cls.lowval d = d + 1 then ∀ t ∈ (s₃).tstack.head?, t.spans.1 = []
    else (∀ b ∈ (s₃).tstack.head?, b.spans.2 = []) ∧ ∀ t ∈ (s₃).tstack.tail.head?, t.spans.1 = [] := by
  sorry

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
    intro _ t ht u hu x htx hux
    obtain ⟨e', he', hte, hix⟩ := htx
    have he's : e' < s.g.ne := by rw [← hne]; exact he'
    have hix' : s.g.Inc e' x := by rw [← hgE]; exact hix
    rcases H.ends e' (hsubE t ht e' he' hte) he's x hix' with rfl | hx
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
      · exact absurd (hpushV hp).2 (by omega)
    rw [hp, Bool.or_false]
    cases hhv : hasVert
    · rfl
    · exact absurd (H.ctx.hv_ret hhv (List.cons_ne_nil _ _))
        (no_ret_before_boundary hge H.cls_ret H.rank)

end

end WalkState
end Spqr
