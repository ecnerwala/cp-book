import Spqr.RangesCloseContent

/-!
# Close invariant through the walk

`CloseInv` is threaded through the mutual walk induction (the shape of `rgTree`/`rgOuts`/`rgOut`):
each `finishEdge` is discharged by `finishEdge_closeInv` from a `CloseCtx`, whose range/ear parts
come from the enclosing induction, whose DFS-side parts (`DfsSite`) are supplied in the `wp` shape of
`RgTree` (`CsTree`/`CsOuts`/`CsOut`) and discharged on the real walk by `dsTree`, and whose remaining
site fields are the named admissions `closeCtx_bd_vert`/`bd_node`/`p_site`/`v_site`/`l1_site`.
-/

namespace Spqr
open WalkM
namespace WalkState
variable {s : WalkState} {σ : List Nat} {n D : Nat}

/-- The DFS-side facts of `CloseCtx` at a `finishEdge` call: the site's position in the schedule,
the edge's endpoints, the open path, and the child's second edge; all follow from `DfsTree.WF`/`Ends`
and the frame (`dsTree`). -/
structure DfsSite (σ : List Nat) (n curV d : Nat) (o : DfsOut) (s : WalkState) : Prop where
  pos : σ[n]? = some o.e
  block : o.block <:+: σ
  ends : o.Ends s.g curV
  path : ∀ k, k < d → ∃ e, e < s.g.ne ∧ n < σ.idxOf e ∧
    Items.PairEq (s.stackVerts[k]!, s.stackVerts[k + 1]!) s.g.edges[e]!
  dest_edge : o.cls.isTree = true → o.cls.lowval d < d →
    ∃ e, e < s.g.ne ∧ e ≠ o.e ∧ s.g.Inc e o.dest
  dest_edges : o.cls.isTree = true → ∀ e, e < s.g.ne → s.g.Inc e o.dest → subEdges o e
  dest_lt : o.dest < s.g.nv
  bd_loop : o.cls.isTree = false → d ≤ o.cls.lowval d → o.dest = curV

/-- Everything the walk induction exports at a `finishEdge` call: the range/ear/frontier/R-side
postconditions (`RgTree`, `gbTree`, `BookTree`, `FrontiersTree`, `scheduleTree`), the close invariant,
and the DFS-side site facts. -/
structure CloseBase (σ : List Nat) (n D curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool)
    (s : WalkState) : Prop where
  hD : D = if o.cls.isTree then d + 1 else d
  nodup : σ.Nodup
  rgs : RgS σ n D s
  guards : FinishGuards d o origTstack hasVert s
  book : FinishBook curV d o origTstack hasVert s
  frontier : Frontier (o := o) d origTstack s
  finishR : FinishR σ n curV d o origTstack hasVert s
  close : s.CloseInv
  site : DfsSite σ n curV d o s


section Sites
variable {curV d : Nat} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}

/-! ### The range invariant inside `finishEdge`, from `CloseBase` -/

/-- `RgStep` to the `k`-th iterate of a loop whose first `k` conditions held. -/
theorem RgStep.iter {v : Nat} (cond : WalkM Bool) (body : WalkM Unit) (Ok Adj : WalkState → Prop)
    (hbody : ∀ s, v < s.g.nv → s.RangesInv σ n D → Shape s → (∀ e ∈ σ, e < s.g.ne) →
      (cond.run s).1 = true → Ok s → Adj s → RgStep σ n D v s (body.run s).2)
    (hv : v < s.g.nv) (h : s.RangesInv σ n D) (hs : Shape s) (hσ : ∀ e ∈ σ, e < s.g.ne)
    (hok : ∀ k, (∀ j, j ≤ k → (cond.run (iter body j s)).1 = true) → Ok (iter body k s))
    (hadj : ∀ k, (∀ j, j ≤ k → (cond.run (iter body j s)).1 = true) → Adj (iter body k s)) :
    ∀ k, (∀ j, j < k → (cond.run (iter body j s)).1 = true) → RgStep σ n D v s (iter body k s)
  | 0, _ => RgStep.refl h hs
  | k + 1, hk => by
    have st := RgStep.iter cond body Ok Adj hbody hv h hs hσ hok hadj k fun j hj => hk j (by omega)
    rw [iter_succ']
    exact st.trans (hbody _ (by rw [st.step.g]; exact hv) st.ranges st.step.shape (st.hσ hσ)
      (hk k (Nat.lt_succ_self k)) (hok k fun j hj => hk j (by omega)) (hadj k fun j hj => hk j (by omega)))

theorem CloseBase.finishOk (h : CloseBase σ n D curV d o origTstack hasVert s) {lv : Nat} {kind : RetKind}
    (ho : o.cls = .ret lv kind) (hl : lv < d) : FinishOk D curV d lv o origTstack hasVert s := by
  obtain ⟨sub, base, hlen, hE⟩ := h.book.ear
  exact finishOk_of_guards ho hl h.guards hE hlen h.rgs.1.inv h.rgs.2.1 h.hD h.book.v_lt h.book.e_lt
    h.book.q (h.book.ends lv kind ho) h.book.vert

theorem CloseBase.step₀ (h : CloseBase σ n D curV d o origTstack hasVert s) :
    RgStep σ n D curV s (feS₀ d o s) :=
  have hj : edgeItem s.g o.e < 1 + s.g.nv + s.g.ne := by
    show 1 + s.g.nv + o.e < _; have := h.book.e_lt; omega
  ⟨Step.modifyVs h.rgs.1.inv h.rgs.2.1 (edgeItem s.g o.e) _ hj, h.rgs.1.modifyVs (edgeItem s.g o.e) _ hj⟩

theorem CloseBase.etype₀ (h : CloseBase σ n D curV d o origTstack hasVert s) :
    Items.type (feS₀ d o s).items (edgeItem (feS₀ d o s).g o.e) ≠ .V := by
  rw [h.step₀.step.shape.edge o.e (by rw [h.step₀.step.g]; exact h.book.e_lt)]; decide

/-- The range invariant (with `o.e` processed) at the `k`-th loop-1 iterate of `closeEars`
(`l1Iter d o s k`), for a returning tree edge whose first `k` loop conditions held. -/
theorem rangesInv_l1Iter (h : CloseBase σ n D curV d o origTstack hasVert s)
    (ht : o.cls.isTree = true) (hlow : o.cls.lowval d < d) (k : Nat)
    (hk : ∀ j, j < k → result (loop1Cond d) (l1Iter d o s j) = true) :
    RgStep σ (n + 1) D curV s (l1Iter d o s k) := by
  obtain ⟨lv, kind, ho, hl⟩ := ret_of_lowval_lt hlow
  have hok := h.finishOk ho hl
  have st₀ := h.step₀
  have hσ₀ := st₀.hσ h.rgs.2.2
  have hv₀ : curV < (feS₀ d o s).g.nv := by rw [st₀.step.g]; exact h.book.v_lt
  have ha := closeEarsAdj_of_frontier h.frontier ht hlow st₀.ranges st₀.step.shape hv₀ h.nodup hσ₀
    h.site.block h.site.pos h.etype₀ (hok.ears ht)
  have st₁ : RgStep σ (n + 1) D curV (feS₀ d o s) (ceS₁ o.dest d o.e (feS₀ d o s)) :=
    RgStep.pushEdge st₀.ranges st₀.step.shape h.nodup o.dest d o.e (hok.ears ht).e_lt (hok.ears ht).q
      ha.etype (hok.ears ht).ends (hok.ears ht).d_le ha.pos
  exact st₀.trans (st₁.trans (RgStep.iter (loop1Cond d) (Spqr.loop1Body d s.stackDir[d]!)
    (Loop1BodyOk D d s.stackDir[d]!) (Loop1BodyAdj σ d s.stackDir[d]!)
    (fun _ hv h' hs hσ _ hok hadj => RgStep.loop1Body h' hs h.nodup hσ hv hok hadj)
    (by rw [st₁.step.g]; exact hv₀) st₁.ranges st₁.step.shape (st₁.hσ hσ₀) (hok.ears ht).body ha.body k hk))

/-- The range invariant (with `o.e` processed) at `feS₂ d o s`, the state after `closeEars` and
`mergeLate` of a returning tree edge, where a type-1 vertex close starts. -/
theorem rangesInv_feS₂ (h : CloseBase σ n D curV d o origTstack hasVert s)
    (ht : o.cls.isTree = true) (hlow : o.cls.lowval d < d) :
    RgStep σ (n + 1) D curV s (feS₂ d o s) := by
  obtain ⟨lv, kind, ho, hl⟩ := ret_of_lowval_lt hlow
  have hok := h.finishOk ho hl
  have st₀ := h.step₀
  have hσ₀ := st₀.hσ h.rgs.2.2
  have hv₀ : curV < (feS₀ d o s).g.nv := by rw [st₀.step.g]; exact h.book.v_lt
  have ha₁ := closeEarsAdj_of_frontier h.frontier ht hlow st₀.ranges st₀.step.shape hv₀ h.nodup hσ₀
    h.site.block h.site.pos h.etype₀ (hok.ears ht)
  have st₁ : RgStep σ (n + 1) D curV s (feS₁ d o s) :=
    st₀.trans (RgStep.closeEars st₀.ranges st₀.step.shape h.nodup hσ₀ hv₀ (hok.ears ht) ha₁)
  have hσ₁ := st₁.hσ h.rgs.2.2
  have ha₂ := mergeLateAdj_of_frontier h.frontier ht hlow st₁.ranges h.nodup hσ₁ h.site.block (hok.late ht)
  exact st₁.trans (RgStep.mergeLate st₁.ranges st₁.step.shape h.nodup hσ₁
    (by rw [st₁.step.g]; exact h.book.v_lt) (hok.late ht) ha₂)


/-- The attachments of a top entry about to be closed (`FinishTopOk`): its bottom, the path vertex
at its `topDepth`, or an interior vertex (`AttachedIn` + `FinishTopOk.mid`). -/
theorem att_of_finishTop {r : WalkState} {t : TEntry} {rest : List TEntry} (hi : r.Inv' D)
    (hts : r.tstack = t :: rest) (hok : FinishTopOk D r) :
    ∀ w, r.g.Touches (t.edges r.g r.items) w →
      w = t.vStart ∨ w = r.stackVerts[t.topDepth]! ∨ r.g.Interior (t.edges r.g r.items) w := by
  intro w hw
  by_cases hint : r.g.Interior (t.edges r.g r.items) w
  · exact .inr (.inr hint)
  have hcur : curE r = t := by simp [curE, hts]
  obtain ⟨e, he, hte, hinc⟩ := hw
  have hint' := hint
  simp only [Graph.Interior, not_forall] at hint'
  obtain ⟨e', he', hinc', hne⟩ := hint'
  have hterm := (hi.entries [] t rest hts).attached w e e' he he' hte hne hinc hinc'
  rw [TEntry.Term'_nil] at hterm
  rcases hterm with h | ⟨k, hk1, hk2, hwk⟩
  · exact .inl h
  rcases Nat.eq_or_lt_of_le hk1 with h | h
  · exact .inr (.inl (by rw [hwk, h]))
  rcases hok.mid k (by rw [hcur]; exact h) hk2 with h1 | h1 | h1
  · rw [hcur] at h1; exact .inl (hwk.trans h1)
  · rw [hcur] at h1; exact absurd (by rw [hwk]; exact h1) hint
  · rw [hcur] at h1; exact absurd (by rw [← hwk]; exact ⟨e, he, hte, hinc⟩) h1

/-- The path edge `(stackVerts[k], stackVerts[k+1])`, `k < d`, is held by no entry of a state
reached from a `finishEdge` pre-state with `o.e` processed. -/
theorem CloseBase.path_pend (h : CloseBase σ n D curV d o origTstack hasVert s) {r : WalkState}
    (st : RgStep σ (n + 1) D curV s r) {k : Nat} (hk : k < d) :
    ∃ e, e < r.g.ne ∧ r.g.Inc e s.stackVerts[k]! ∧ r.g.Inc e s.stackVerts[k + 1]! ∧
      ∀ t ∈ r.tstack, ¬ t.edges r.g r.items e := by
  obtain ⟨e, he, hn, hpe⟩ := h.site.path k hk
  have hg := st.step.g
  refine ⟨e, by rw [hg]; exact he, ?_, ?_, fun t ht hte => ?_⟩
  · rw [hg]; exact (Graph.inc_of_pairEq hpe).1
  · rw [hg]; exact (Graph.inc_of_pairEq hpe).2
  · have := st.ranges.processed t ht e (by rw [hg]; exact he) hte
    omega

/-- The range invariant (with `o.e` processed) at the state `finishP` runs from. -/
theorem CloseBase.rgFeRest (h : CloseBase σ n D curV d o origTstack hasVert s) {lv : Nat} {kind : RetKind}
    (ho : o.cls = .ret lv kind) (hl : lv < d) :
    RgStep σ (n + 1) D curV s (feRest curV d o origTstack hasVert s) := by
  have hlv : o.cls.lowval d = lv := by rw [ho]; rfl
  have hlow : o.cls.lowval d < d := by rwa [hlv]
  have hok := h.finishOk ho hl
  have hadj := h.finishR.1 lv kind ho hl
  have hnd := h.nodup
  by_cases ht : o.cls.isTree = true
  · have st₂ := rangesInv_feS₂ h ht hlow
    have hσ₂ := st₂.hσ h.rgs.2.2
    have hv₂ : curV < (feS₂ d o s).g.nv := by rw [st₂.step.g]; exact h.book.v_lt
    cases hasVert
    · simpa [feRest, ht] using st₂
    · have st₃ : RgStep σ (n + 1) D curV _ (feS₃ curV d o origTstack s) :=
        RgStep.closeVert' st₂.ranges st₂.step.shape hnd hσ₂ hv₂ (hok.vert ht rfl) (hadj.vert ht rfl)
      simpa [feRest, ht] using st₂.trans st₃
  · have ht' : o.cls.isTree = false := Bool.eq_false_iff.2 ht
    have st₀ := h.step₀
    have hq : Items.ch (feS₀ d o s).items (edgeItem (feS₀ d o s).g o.e) = [] := by
      show Items.ch (s.items.modify _ _) (edgeItem s.g o.e) = []
      rw [Items.ch_modify_ch_eq (edgeItem s.g o.e)
        (fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) (fun _ => rfl)]
      exact hok.q ht'
    have st₁ : RgStep σ (n + 1) D curV _ (after (pushEdgeTstack curV lv o.e) (feS₀ d o s)) :=
      RgStep.pushEdge st₀.ranges st₀.step.shape hnd curV lv o.e hok.e_lt hq (hadj.etype ht') (hok.ends ht')
        (hok.lv_le ht') (hadj.pos ht')
    have st₂ : RgStep σ (n + 1) D curV _ (feBack curV lv d o s) :=
      ⟨Step.frame st₁.step.inv st₁.step.shape rfl rfl rfl rfl, st₁.ranges.frame rfl rfl rfl rfl⟩
    have hr : feRest curV d o origTstack hasVert s = feBack curV lv d o s := by simp [feRest, ht', hlv]
    rw [hr]; exact st₀.trans (st₁.trans st₂)

/-- The range invariant through the merges and the retarget of `closeVert'` (before the type-1
`finishTstackTop`). -/
theorem RgStep.cvS₅ {v : Nat} (h : s.RangesInv σ n D) (hs : Shape s) (hnd : σ.Nodup)
    (hσ : ∀ e ∈ σ, e < s.g.ne) (hv : v < s.g.nv) {curV : Nat} {edgeDir isType1 : Bool}
    {origTstack : Nat} {isSingle : Bool}
    (hok : CloseVertOk D curV edgeDir isType1 origTstack isSingle s)
    (hadj : CloseVertAdj σ isType1 origTstack isSingle s) :
    RgStep σ n D v s (cvS₅ curV edgeDir isType1 origTstack isSingle s) := by
  have st₁ : RgStep σ n D v s (cvS₁ isType1 origTstack isSingle s) :=
    RgStep.vertPre h hs hnd hσ hv hok.loop3 hadj.loop3
  have hσ₁ := st₁.hσ hσ
  obtain ⟨st₂, -⟩ := RgStep.vertUnwrap (v := v) st₁.ranges st₁.step.shape hσ₁ (isType1 := isType1)
    (isSingle := cvB₁ isType1 origTstack isSingle s) (fun h => by subst h; exact hok.unwrap rfl)
  have st₂ : RgStep σ n D v (cvS₁ isType1 origTstack isSingle s) (cvS₂ isType1 origTstack isSingle s) := st₂
  have hσ₂ := st₂.hσ hσ₁
  have st₃ : RgStep σ n D v _ (cvS₃ isType1 origTstack isSingle s) :=
    RgStep.mergeTop st₂.ranges st₂.step.shape hnd hσ₂ hok.merge₁ hadj.merge₁
  have hσ₃ := st₃.hσ hσ₂
  have st₄ : RgStep σ n D v _ (cvS₄ isType1 origTstack isSingle s) :=
    RgStep.mergeTop st₃.ranges st₃.step.shape hnd hσ₃ hok.merge₂ hadj.merge₂
  exact st₁.trans (st₂.trans (st₃.trans (st₄.trans
    (RgStep.retarget st₄.ranges st₄.step.shape (st₄.hσ hσ₃) curV edgeDir hok.retarget))))

/-- `maybeUnwrapNxt` keeps the stack's entries and their `topDepth` (it only rewrites `nxt`'s spans). -/
theorem unwrap_topDepth {ty : NodeType} {a b : TEntry} {rest : List TEntry} (hts : s.tstack = a :: b :: rest) :
    ∃ b', (after (maybeUnwrapNxt ty) s).tstack = a :: b' :: rest ∧ b'.topDepth = b.topDepth := by
  rw [show after (maybeUnwrapNxt ty) s = ((maybeUnwrapNxt ty).run s).2 from rfl,
    maybeUnwrapNxt_run_eq ty s a b rest hts _ rfl _ rfl]
  split_ifs
  · exact ⟨b, by rw [run_allocItem]; exact hts, rfl⟩
  · exact ⟨_, rfl, rfl⟩
  · exact ⟨b, by rw [run_allocItem]; exact hts, rfl⟩

/-- At a type-1 vertex close, the entry retargeted at `cvS₅` has `topDepth = lowval`: the unwrap
keeps `py`'s depth, the two merges take the minimum with `c` (`≥ lowval`) and `vy` (`> d`). -/
theorem cvS₅_topDepth (ht : o.cls.isTree = true) (hlow : o.cls.lowval d < d) (h1 : o.cls.isType1 = true)
    {sub base : List TEntry} (hE : s.EarFinish curV d o hasVert sub base) {t : TEntry} {rest : List TEntry}
    (hts : (cvS₅ curV s.stackDir[d]! true origTstack (feSingle d o s) (feS₂ d o s)).tstack = t :: rest) :
    t.topDepth = o.cls.lowval d := by
  obtain ⟨c, mid, py, vy, hcl⟩ := hE.close ht hlow
  have hmid : mid = [] := (hcl.type1 h1).1
  have hts₂ : (feS₂ d o s).tstack = c :: py :: vy :: base := by simp [hcl.tstack, hmid]
  have hS₂ : cvS₂ true origTstack (feSingle d o s) (feS₂ d o s) =
      after (maybeUnwrapNxt (if feSingle d o s then .S else .R)) (feS₂ d o s) := by
    simp only [cvS₂, cvS₁, cvB₁, after, result, vertPre, vertUnwrap, Bool.not_true, Bool.false_eq_true,
      ↓reduceIte, WalkM.pure_run, WalkM.map_run]
  obtain ⟨py', hts₂', hpy'⟩ := unwrap_topDepth (ty := if feSingle d o s then .S else .R) hts₂
  have hts₃ : (cvS₃ true origTstack (feSingle d o s) (feS₂ d o s)).tstack =
      TEntry.mergeInto c py' :: vy :: base := by
    simp only [cvS₃, hS₂]; rw [after_mergeTstackTops_eq hts₂']
  have hts₄ : (cvS₄ true origTstack (feSingle d o s) (feS₂ d o s)).tstack =
      TEntry.mergeInto (TEntry.mergeInto c py') vy :: base := by
    simp only [cvS₄]; rw [after_mergeTstackTops_eq hts₃]
  have hts₅ : (cvS₅ curV s.stackDir[d]! true origTstack (feSingle d o s) (feS₂ d o s)).tstack =
      { TEntry.mergeInto (TEntry.mergeInto c py') vy with
        vStart := curV,
        spans := setSides (!s.stackDir[d]!)
          ((TEntry.mergeInto (TEntry.mergeInto c py') vy).spans.1 ++
            (TEntry.mergeInto (TEntry.mergeInto c py') vy).spans.2) [] } :: base := by
    show ((WalkState.retarget _ _).run _).2.tstack = _
    rw [retarget_run_eq _ _ _ _ _ hts₄]
  rw [hts₅] at hts
  obtain ⟨rfl, -⟩ := List.cons.inj hts
  have h₁ := hcl.c_top.1
  have h₂ := hcl.py_top
  have h₃ := hcl.vy_top
  simp only [TEntry.mergeInto, hpy', h₂]
  omega

end Sites

section Admissions
variable {curV d : Nat} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}

/-- Block entry of a boundary edge carries exactly `vertItem o.dest` on its active side
(checker: `closeCtx_bd_vert`, kind `ctx_bd_vert`). -/
theorem closeCtx_bd_vert (_h : CloseBase σ n D curV d o origTstack hasVert s)
    (hc : CloseContent curV d o origTstack hasVert s) :
    o.cls.isTree = true → d ≤ o.cls.lowval d →
    ∀ t ∈ (if o.cls.lowval d = d + 1 then s.tstack.head? else s.tstack.tail.head?),
      t.spans.2 = [vertItem o.dest] :=
  hc.bd_vert

/-- Completed block: the top entry's passive side is one non-`F`/`V` node (a leaf if `Q`) with
terminals `{curV, o.dest}` (checker: `closeCtx_bd_node`, kind `ctx_bd_node`). -/
theorem closeCtx_bd_node (_h : CloseBase σ n D curV d o origTstack hasVert s)
    (hc : CloseContent curV d o origTstack hasVert s) :
    o.cls.isTree = true → d ≤ o.cls.lowval d → o.cls.lowval d ≠ d + 1 →
    ∀ b ∈ s.tstack.head?, ∃ c, b.spans.1 = [c] ∧ Items.type s.items c ∉ [NodeType.F, .V] ∧
      (Items.type s.items c = .Q → Items.ch s.items c = []) ∧
      ∃ a b', Items.vs s.items c = (some a, some b') ∧ Items.PairEq (a, b') (curV, o.dest) :=
  hc.bd_node

/-- The P-merge site (checker: `closeCtx_p_site`, kinds `psite_*`). -/
theorem closeCtx_p_site (h : CloseBase σ n D curV d o origTstack hasVert s)
    (hc : CloseContent curV d o origTstack hasVert s) :
    o.cls.lowval d < d → o.cls.isType1 = true →
    result (condP curV (o.cls.lowval d) true) (feRest curV d o origTstack hasVert s) = true →
    PSite curV (o.cls.lowval d) (feRest curV d o origTstack hasVert s) := by
  intro hlow h1 hp
  obtain ⟨lv, kind, ho, hl⟩ := ret_of_lowval_lt hlow
  have st := h.rgFeRest ho hl
  obtain ⟨sub, base, hlen, hE⟩ := h.book.ear
  have pc := hc.p_site hlow h1 hp
  have hsv := st.step.sv
  exact {
    shape := st.step.shape, stack := pc.stack, single := pc.single, att := pc.att,
    touch := pc.touch, kinds := pc.kinds, once := pc.once,
    ne := by
      rw [hsv]; intro heq
      exact hE.path _ _ hlow le_rfl (hE.sv_d.trans heq).symm,
    pend := by
      intro w hw
      rcases hw with rfl | rfl
      · obtain ⟨e, he, -, hinc, hpend⟩ := h.path_pend st (k := d - 1) (by omega)
        refine ⟨e, he, ?_, fun t ht => hpend t (List.mem_of_mem_take ht)⟩
        rwa [show d - 1 + 1 = d by omega, hE.sv_d] at hinc
      · obtain ⟨e, he, hinc, -, hpend⟩ := h.path_pend st hlow
        exact ⟨e, he, by rw [hsv]; exact hinc, fun t ht => hpend t (List.mem_of_mem_take ht)⟩ }

/-- The type-1 vertex-close site (checker: `closeCtx_v_site`, kinds `vsite_*`). -/
theorem closeCtx_v_site (h : CloseBase σ n D curV d o origTstack hasVert s)
    (hc : CloseContent curV d o origTstack hasVert s) :
    o.cls.isTree = true → o.cls.lowval d < d → hasVert = true → o.cls.isType1 = true →
    ∃ t, VSite curV ((maybeUnwrapNxt (if feSingle d o s then NodeType.S else .R)).run (feS₂ d o s)).1 t
      (cvS₅ curV s.stackDir[d]! true origTstack (feSingle d o s) (feS₂ d o s)) := by
  intro ht hlow hv h1
  obtain ⟨lv, kind, ho, hl⟩ := ret_of_lowval_lt hlow
  obtain ⟨t, vc⟩ := hc.v_site ht hlow hv h1
  obtain ⟨sub, base, hlen, hE⟩ := h.book.ear
  have hok := (h.finishOk ho hl).vert ht hv
  have hadj := (h.finishR.1 lv kind ho hl).vert ht hv
  rw [h1] at hok hadj
  have st₂ := rangesInv_feS₂ h ht hlow
  have st₅ := RgStep.cvS₅ st₂.ranges st₂.step.shape h.nodup (st₂.hσ h.rgs.2.2)
    (by rw [st₂.step.g]; exact h.book.v_lt) hok hadj
  have st := st₂.trans st₅
  obtain ⟨rest, hts⟩ := vc.stack
  have htop := cvS₅_topDepth ht hlow h1 hE hts
  have hsv := st.step.sv
  refine ⟨t, {
    shape := st₅.step.shape, stack := vc.stack, vstart := vc.vstart, side := vc.side,
    free := vc.free, kinds := vc.kinds, two := vc.two, touch := vc.touch, inner := vc.inner,
    s_order := vc.s_order, r_shape := vc.r_shape, p_shape := vc.p_shape,
    ne := by
      rw [htop, hsv]; intro heq
      exact hE.path _ _ hlow le_rfl (hE.sv_d.trans heq).symm,
    att := fun w hw => ?_,
    pend := by
      intro w hw
      rcases hw with rfl | rfl
      · obtain ⟨e, he, -, hinc, hpend⟩ := h.path_pend st (k := d - 1) (by omega)
        refine ⟨e, he, ?_, hpend t (hts ▸ List.mem_cons_self ..)⟩
        rwa [show d - 1 + 1 = d by omega, hE.sv_d] at hinc
      · obtain ⟨e, he, hinc, -, hpend⟩ := h.path_pend st hlow
        exact ⟨e, he, by rw [htop, hsv]; exact hinc, hpend t (hts ▸ List.mem_cons_self ..)⟩ }⟩
  rcases att_of_finishTop st₅.ranges.inv hts (hok.finish rfl) w hw with h' | h' | h'
  exacts [.inl (h'.trans vc.vstart), .inr (.inl h'), .inr (.inr h')]

/-- The loop-1 iteration sites (checker: `closeCtx_l1_site`, kinds `l1site_*`). -/
theorem closeCtx_l1_site (h : CloseBase σ n D curV d o origTstack hasVert s)
    (hc : CloseContent curV d o origTstack hasVert s) :
    o.cls.isTree = true → o.cls.lowval d < d → ∀ k,
    (∀ j, j ≤ k → result (loop1Cond d) (l1Iter d o s j) = true) →
    Shape (l1S₁ d s.stackDir[d]! (l1Iter d o s k)) ∧
    (∃ a b rest, (l1S₁ d s.stackDir[d]! (l1Iter d o s k)).tstack = a :: b :: rest) ∧
    ∃ t, VSite t.vStart
      (result (maybeUnwrapNxt (l1Ty d s.stackDir[d]! (l1Iter d o s k))) (l1S₁ d s.stackDir[d]! (l1Iter d o s k)))
      t (after mergeTstackTops (l1S₂ d s.stackDir[d]! (l1Iter d o s k))) := by
  intro ht hlow k hk
  obtain ⟨lv, kind, ho, hl⟩ := ret_of_lowval_lt hlow
  have hok : Loop1BodyOk D d s.stackDir[d]! (l1Iter d o s k) := ((h.finishOk ho hl).ears ht).body k hk
  have hadj : Loop1BodyAdj σ d s.stackDir[d]! (l1Iter d o s k) :=
    ((h.finishR.1 lv kind ho hl).ears ht).body k hk
  have st := rangesInv_l1Iter h ht hlow k (fun j hj => hk j hj.le)
  have hσ := st.hσ h.rgs.2.2
  have st₁ : RgStep σ (n + 1) D curV _ (l1S₁ d s.stackDir[d]! (l1Iter d o s k)) :=
    RgStep.loop1Type st.ranges st.step.shape h.nodup hσ hok.mergeS hadj.mergeS
  have hσ₁ := st₁.hσ hσ
  have r := RgStep.unwrapRes (v := curV) st₁.ranges st₁.step.shape hσ₁ (loop1Type_result d _ _) hok.unwrap
  have st₂ : RgStep σ (n + 1) D curV _ (l1S₂ d s.stackDir[d]! (l1Iter d o s k)) := r.step
  have st₃ : RgStep σ (n + 1) D curV _ (l1Pre d o s k) :=
    RgStep.mergeTop st₂.ranges st₂.step.shape h.nodup (st₂.hσ hσ₁) hok.close.merge hadj.close
  obtain ⟨t, hne, hpend, vc⟩ := hc.l1_site ht hlow k hk
  obtain ⟨rest, hts⟩ := vc.stack
  refine ⟨st₁.step.shape, ?_, t, {
    shape := st₃.step.shape, stack := vc.stack, vstart := vc.vstart,
    side := vc.side, free := vc.free, ne := hne, kinds := vc.kinds, two := vc.two, touch := vc.touch,
    pend := hpend, inner := vc.inner, s_order := vc.s_order, r_shape := vc.r_shape, p_shape := vc.p_shape,
    att := att_of_finishTop st₃.ranges.inv hts hok.close.finish }⟩
  have h2 := hok.unwrap.two
  rcases hts₁ : (l1S₁ d s.stackDir[d]! (l1Iter d o s k)).tstack with _ | ⟨a, _ | ⟨b, rest⟩⟩
  · rw [hts₁] at h2; simp at h2
  · rw [hts₁] at h2; simp at h2
  · exact ⟨a, b, rest, rfl⟩

/-- `CloseCtx` from the exports and the block contents. -/
theorem CloseCtx.of_exports (h : CloseBase σ n D curV d o origTstack hasVert s)
    (hc : CloseContent curV d o origTstack hasVert s) :
    CloseCtx σ n D curV d o origTstack hasVert s where
  hD := h.hD
  nodup := h.nodup
  lt := h.rgs.2.2
  pos := h.site.pos
  block := h.site.block
  ranges := h.rgs.1
  shape := h.rgs.2.1
  guards := h.guards
  book := h.book
  frontier := h.frontier
  finishR := h.finishR
  close := h.close
  ends := h.site.ends
  path := h.site.path
  dest_edge := h.site.dest_edge
  dest_edges := h.site.dest_edges
  dest_lt := h.site.dest_lt
  bd_loop := h.site.bd_loop
  bd_vert := closeCtx_bd_vert h hc
  bd_node := closeCtx_bd_node h hc
  p_site := closeCtx_p_site h hc
  v_site := closeCtx_v_site h hc
  l1_site := closeCtx_l1_site h hc

/-- The block contents at every `finishEdge` pre-state; supplied by the walk induction
(`WalkBackbone.lean`, `WalkInvEnd`), a named hypothesis here. -/
theorem closeBase_content (h : CloseBase σ n D curV d o origTstack hasVert s) :
    CloseContent curV d o origTstack hasVert s := by
  sorry

theorem CloseCtx.of_base (h : CloseBase σ n D curV d o origTstack hasVert s) :
    CloseCtx σ n D curV d o origTstack hasVert s :=
  CloseCtx.of_exports h (closeBase_content h)

end Admissions

mutual
def CsTree (σ : List Nat) (n : Nat) (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  match t with
  | .node v outs => CsOuts σ n v d outs false { s with stackVerts := s.stackVerts.set! d v }

def CsOuts (σ : List Nat) (n v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  match outs with
  | [] => True
  | o :: rest => CsOut σ n v d o hasVert s ∧
      wp (walkOut v d o hasVert) (fun hasVert' s' => CsOuts σ (n + o.block.length) v d rest hasVert' s') s

def CsOut (σ : List Nat) (n v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  wp (walkOutPre v d o hasVert) (fun hasVert' s₁ =>
    match o with
    | .tree _ _ child =>
      wp (modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }) (fun _ s₂ =>
        CsTree σ n child (d + 1) s₂ ∧
        wp (walkTree child (d + 1)) (fun _ s₃ =>
          DfsSite σ (n + child.edgePostorder.length) v d o s₃) s₂) s₁
    | .back .. => DfsSite σ n v d o s₁) s
end

theorem walkOutPre_closeInv' (h : s.CloseInv) (hs : Shape s) {v : Nat} (d : Nat)
    (o : DfsOut) {hasVert : Bool} (hb : VertBook v hasVert s) :
    wp (walkOutPre v d o hasVert) (fun _ s' => s'.CloseInv) s := by
  have h' : ({ s with stackDir := s.stackDir.set! d (if o.cls.lowval d ≥ d then false else !s.stackDir[o.cls.lowval d]!) } : WalkState).CloseInv :=
    h.frame rfl rfl (fun _ hi => hi)
  unfold walkOutPre
  simp only [wp_bind, wp_stackDir, wp_setStackDir, wp_ite]
  split
  · rename_i hc
    have hf : hasVert = false := by cases hasVert <;> simp_all
    obtain ⟨hv, -, -⟩ := hb hf
    simp only [wp_bind, wp_pure]
    exact h'.pushVert v d hv (by have := hs.size; show s.g.nv < s.items.size; omega)
  · simp only [wp_pure]; exact h'

abbrev CcInvTree (σ : List Nat) (n : Nat) (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  (∀ v outs, t = .node v outs → ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).RangesInv σ n d) →
  Shape s → σ.Nodup → (∀ e ∈ σ, e < s.g.ne) → GuardsTree t d s → BookTree t d s → FrontiersTree t d s →
  RgTree σ n t d s → CsTree σ n t d s → s.CloseInv →
  wp (walkTree t d) (fun _ s' => s'.CloseInv) s

abbrev CcInvOuts (σ : List Nat) (n v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  RgS σ n d s → σ.Nodup → GuardsOuts v d outs hasVert s → BookOuts v d outs hasVert s →
  FrontiersOuts v d outs hasVert s → RgOuts σ n v d outs hasVert s → CsOuts σ n v d outs hasVert s →
  s.CloseInv → wp (walkOuts v d outs hasVert) (fun _ s' => s'.CloseInv) s

abbrev CcInvOut (σ : List Nat) (n v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  RgS σ n d s → σ.Nodup → GuardsOut v d o hasVert s → BookOut v d o hasVert s →
  FrontiersOut v d o hasVert s → RgOut σ n v d o hasVert s → CsOut σ n v d o hasVert s →
  s.CloseInv → wp (walkOut v d o hasVert) (fun _ s' => s'.CloseInv) s

mutual
theorem ccTree : ∀ (σ : List Nat) (n : Nat) (t : DfsTree) (d : Nat) (s : WalkState), CcInvTree σ n t d s
  | σ, n, .node v outs, d, s => fun hi hs hnd hσ hg hb hf hr hc hcl => by
    unfold walkTree
    simp only [wp_bind, wp_modify]
    unfold GuardsTree at hg; unfold BookTree at hb; unfold FrontiersTree at hf; unfold RgTree at hr
    unfold CsTree at hc
    refine wp_imp (wp_imp (wp_of_forall fun hv s' ⟨⟨hi', hs', hσ'⟩, hvb, hpr⟩ hcl' => ?_)
      (rgOuts σ n v d outs false _ ⟨hi v outs rfl, hs.frame', hσ⟩ hnd hg hb hr))
      (ccOuts σ n v d outs false _ ⟨hi v outs rfl, hs.frame', hσ⟩ hnd hg hb hf hr hc
        (hcl.frame rfl rfl (fun _ h => h)))
    cases hv
    · simp only [Bool.false_eq_true, ↓reduceIte, wp_bind, wp_setStackDir]
      obtain ⟨hv, -, -⟩ := hvb rfl
      exact (hcl'.frame (s' := { s' with stackDir := s'.stackDir.set! d true }) rfl rfl (fun _ h => h)).pushVert
        v d hv (by have := hs'.size; show s'.g.nv < s'.items.size; omega)
    · exact hcl'

theorem ccOuts : ∀ (σ : List Nat) (n v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState),
    CcInvOuts σ n v d outs hasVert s
  | σ, n, v, d, [], hasVert, s => fun _ _ _ _ _ _ _ hcl => by
    unfold walkOuts
    simp only [wp_pure]
    exact hcl
  | σ, n, v, d, o :: rest, hasVert, s => fun hrs hnd hg hb hf hr hc hcl => by
    unfold GuardsOuts at hg; unfold BookOuts at hb; unfold FrontiersOuts at hf; unfold RgOuts at hr
    unfold CsOuts at hc
    unfold walkOuts
    simp only [wp_bind]
    exact wp_imp (wp_imp (wp_imp (wp_imp (wp_imp (wp_imp (wp_imp (wp_of_forall
      fun hv' s' hrs' hg' hb' hf' hr' hc' hcl' =>
        ccOuts σ (n + o.block.length) v d rest hv' s' hrs' hnd hg' hb' hf' hr' hc' hcl')
      (rgOut σ n v d o hasVert s hrs hnd hg.1 hb.1 hr.1)) hg.2) hb.2) hf.2) hr.2) hc.2)
      (ccOut σ n v d o hasVert s hrs hnd hg.1 hb.1 hf.1 hr.1 hc.1 hcl)

theorem ccOut : ∀ (σ : List Nat) (n v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState), CcInvOut σ n v d o hasVert s
  | σ, n, v, d, o, hasVert, s => fun ⟨hi, hs, hσ⟩ hnd hg hb hf hr hc hcl => by
    rw [walkOut_eq, wp_bind]
    unfold GuardsOut at hg; unfold BookOut at hb; unfold FrontiersOut at hf; unfold RgOut at hr
    unfold CsOut at hc
    refine wp_imp (wp_of_forall fun hv' s₁ ⟨⟨hi₁, hs₁, hσ₁⟩, hg₁, hb₁, hf₁, hr₁, hc₁, hcl₁⟩ => ?_)
      (wp_and (walkOutPre_ranges hi hs hσ hb.1 hr.1) (wp_and hg (wp_and hb.2 (wp_and hf (wp_and hr.2
        (wp_and hc (walkOutPre_closeInv' hcl hs d o hb.1)))))))
    unfold walkOutRest
    rw [wp_bind, wp_tstackSize]
    cases o with
    | back e cls dest =>
      try simp only [wp_pure] at hg₁ hb₁ hf₁ hr₁ hc₁
      simp only [wp_bind, wp_pure]
      exact finishEdge_closeInv (CloseCtx.of_base (D := d)
        ⟨by rw [if_neg fun h => by obtain ⟨_, _, _, h⟩ := hb₁.tree.1 h; cases h],
          hnd, ⟨hi₁, hs₁, hσ₁⟩, hg₁, hb₁, hf₁, hr₁, hcl₁, hc₁⟩)
    | tree e cls child =>
      try simp only [wp_bind, wp_modify] at hg₁ hb₁ hf₁ hr₁ hc₁
      simp only [wp_bind, wp_modify]
      have pre : ∀ w outs, child = .node w outs →
          ({ { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne } with
            stackVerts := s₁.stackVerts.set! (d + 1) w } : WalkState).RangesInv σ n (d + 1) :=
        fun w outs _ =>
          (hi₁.frame' (s' := { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne })).setSv w
      refine wp_imp (wp_imp (wp_imp (wp_of_forall fun _ s₃ ⟨hg₃, hb₃, hf₃, hr₃, hc₃⟩ ⟨hi₃, hs₃, hσ₃⟩ hcl₃ => ?_)
        (wp_and hg₁.2 (wp_and hb₁.2 (wp_and hf₁.2 (wp_and hr₁.2 hc₁.2)))))
        (rgTree σ n child (d + 1) _ pre hs₁.frame' hnd hσ₁ hg₁.1 hb₁.1 hr₁.1))
        (ccTree σ n child (d + 1) _ pre hs₁.frame' hnd hσ₁ hg₁.1 hb₁.1 hf₁.1 hr₁.1 hc₁.1
          (hcl₁.frame rfl rfl (fun _ h => h)))
      exact finishEdge_closeInv (CloseCtx.of_base (D := d + 1)
        ⟨by rw [if_pos (hb₃.tree.2 ⟨_, _, _, rfl⟩)], hnd, ⟨hi₃, hs₃, hσ₃⟩, hg₃, hb₃, hf₃, hr₃, hcl₃, hc₃⟩)
end

/-- `walkTree` preserves `CloseInv`, given the range/ear hypotheses of `walkTree_rangesInv`, the
frontiers, and the per-site facts `CsTree`. -/
theorem walkTree_closeInv (t : DfsTree) (d : Nat) (s : WalkState)
    (hi : ∀ v outs, t = .node v outs → ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).RangesInv σ n d)
    (hs : Shape s) (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne) (hg : GuardsTree t d s) (hb : BookTree t d s)
    (hf : FrontiersTree t d s) (hr : RgTree σ n t d s) (hc : CsTree σ n t d s) (hcl : s.CloseInv) :
    ((walkTree t d).run s).2.CloseInv :=
  ccTree σ n t d s hi hs hnd hσ hg hb hf hr hc hcl

/-! ### The site facts on the real walk -/

/-- The open path `sv` (`sv[k]` at depth `k`) is joined by edges scheduled at or after `m`. -/
def AncPath (g : Graph) (σ : List Nat) (m : Nat) (sv : List Nat) : Prop :=
  ∀ k, k + 1 < sv.length → ∃ e, e < g.ne ∧ m ≤ σ.idxOf e ∧ Items.PairEq (sv[k]!, sv[k + 1]!) g.edges[e]!

theorem AncPath.mono {g : Graph} {m m' : Nat} {sv : List Nat} (h : m' ≤ m) (hp : AncPath g σ m sv) :
    AncPath g σ m' sv := fun k hk => by
  obtain ⟨e, he, hm, hpe⟩ := hp k hk
  exact ⟨e, he, le_trans h hm, hpe⟩

theorem getElem!_of_getElem? {l : List Nat} {k x : Nat} (h : l[k]? = some x) : l[k]! = x := by
  rw [List.getElem!_eq_getElem?_getD, h]; rfl

theorem getElem!_append_left' {l₁ l₂ : List Nat} {k : Nat} (h : k < l₁.length) : (l₁ ++ l₂)[k]! = l₁[k]! := by
  simp [List.getElem!_eq_getElem?_getD, List.getElem?_append_left h]

theorem getElem!_concat_length' (l : List Nat) (a : Nat) : (l ++ [a])[l.length]! = a :=
  getElem!_of_getElem? List.getElem?_concat_length

theorem PairEq.flip {a b : Nat} {q : Nat × Nat} (h : Items.PairEq (a, b) q) : Items.PairEq (b, a) q := by
  rcases h with h | h
  · subst h; exact .inr rfl
  · obtain ⟨rfl, rfl⟩ := Prod.mk.inj h; exact .inl rfl

theorem classify_back {d i : Nat} (hi : i < d + 1) :
    classify d false (i, d) = if i = d then .selfLoop else .ret i .backEdge := by
  unfold classify
  by_cases h : i = d
  · subst h; simp
  · have : ¬ (i ≥ d) := by omega
    simp [h, this]

abbrev DsTree (σ : List Nat) (n : Nat) (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  ∀ (g : Graph) (anc : List Nat) (pe : Nat → Prop), Types g s → d = anc.length → t.WF anc → t.Ends g →
    (∀ e', e' < g.ne → ∀ x, g.Inc e' x → x ∈ t.verts → e' ∈ t.edges ∨ pe e') →
    (∀ e', pe e' → ∀ x, g.Inc e' x → x ∈ anc ++ [t.v]) →
    (anc ++ t.verts).Nodup → (∀ a ∈ anc, a < g.nv) → (∀ v ∈ t.verts, v < g.nv) →
    (∀ e ∈ t.edges, e < g.ne) → t.edges.Nodup → s.stackVerts.size = g.nv →
    (∀ k, k < d → anc[k]? = some s.stackVerts[k]!) →
    σ.Nodup → PostAt σ n t.edgePostorder →
    (∀ v outs, t = .node v outs → AncPath g σ (n + t.edgePostorder.length) (anc ++ [v])) →
    CsTree σ n t d s
abbrev DsOuts (σ : List Nat) (n v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  ∀ (g : Graph) (anc : List Nat) (pe : Nat → Prop), Types g s → d = anc.length → (∀ o ∈ outs, o.WF anc v) →
    (∀ o ∈ outs, DfsOut.Ends g v o) →
    (∀ e', e' < g.ne → ∀ x, g.Inc e' x → x ∈ DfsOut.vertsList outs → e' ∈ DfsOut.edgesList outs ∨ pe e') →
    (∀ e', pe e' → ∀ x, g.Inc e' x → x ∈ anc ++ [v]) →
    (anc ++ v :: DfsOut.vertsList outs).Nodup →
    (∀ a ∈ anc, a < g.nv) → v < g.nv → (∀ w ∈ DfsOut.vertsList outs, w < g.nv) →
    (∀ e ∈ DfsOut.edgesList outs, e < g.ne) → (DfsOut.edgesList outs).Nodup →
    s.stackVerts.size = g.nv → (∀ k, k < d → anc[k]? = some s.stackVerts[k]!) → s.stackVerts[d]! = v →
    σ.Nodup → PostAt σ n (DfsOut.edgePostorderList outs) →
    AncPath g σ (n + (DfsOut.edgePostorderList outs).length) (anc ++ [v]) →
    CsOuts σ n v d outs hasVert s
abbrev DsOut (σ : List Nat) (n v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  ∀ (g : Graph) (anc : List Nat) (pe : Nat → Prop), Types g s → d = anc.length → o.WF anc v →
    DfsOut.Ends g v o →
    (∀ e', e' < g.ne → ∀ x, g.Inc e' x → x ∈ o.verts → subEdges o e' ∨ pe e') →
    (∀ e', pe e' → ∀ x, g.Inc e' x → x ∈ anc ++ [v]) →
    (anc ++ v :: o.verts).Nodup →
    (∀ a ∈ anc, a < g.nv) → v < g.nv → (∀ w ∈ o.verts, w < g.nv) →
    (∀ e ∈ o.edges, e < g.ne) → o.edges.Nodup →
    s.stackVerts.size = g.nv → (∀ k, k < d → anc[k]? = some s.stackVerts[k]!) → s.stackVerts[d]! = v →
    σ.Nodup → PostAt σ n o.block →
    AncPath g σ (n + o.block.length) (anc ++ [v]) →
    CsOut σ n v d o hasVert s

mutual
theorem dsTree : ∀ (σ : List Nat) (n : Nat) (t : DfsTree) (d : Nat) (s : WalkState), DsTree σ n t d s
  | σ, n, .node v outs, d, s => by
    intro g anc pe hT hd hwf hends hcomp hpeA hnd hanc hvlt helt hen hsz hsvk hσn hat hpath
    simp only [DfsTree.verts, DfsTree.edges] at hnd hvlt helt hen
    simp only [DfsTree.WF] at hwf
    simp only [DfsTree.Ends] at hends
    simp only [DfsTree.edgePostorder] at hat
    have hp := hpath v outs rfl
    simp only [DfsTree.edgePostorder] at hp
    unfold CsTree
    have hv : v < g.nv := hvlt v (List.mem_cons_self ..)
    have hdlt : d < g.nv := by
      have hnd' : (anc ++ [v]).Nodup :=
        List.Nodup.sublist (List.Sublist.append_left (List.cons_sublist_cons.2 (List.nil_sublist _)) _) hnd
      have hsub : anc ++ [v] ⊆ List.range g.nv := fun w hw => by
        rw [List.mem_range]
        rcases List.mem_append.1 hw with hw | hw
        · exact hanc w hw
        · exact (List.mem_singleton.1 hw) ▸ hv
      have := (List.Nodup.subperm hnd' hsub).length_le
      simp at this; omega
    exact dsOuts σ n v d outs false _ g anc pe ⟨hT.g_eq, hT.size, hT.root, hT.vert, hT.edge⟩ hd hwf.2 hends
      (fun e' he' x hx hxm => hcomp e' he' x hx (List.mem_cons_of_mem _ hxm)) hpeA
      hnd hanc hv (fun w hw => hvlt w (List.mem_cons_of_mem _ hw)) helt hen (by simp [hsz])
      (fun k hk => by rw [getElem!_set!_ne _ _ _ _ (by omega)]; exact hsvk k hk)
      (getElem!_set!_self _ _ _ (by rw [hsz]; exact hdlt)) hσn hat hp

theorem dsOuts : ∀ (σ : List Nat) (n v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState),
    DsOuts σ n v d outs hasVert s
  | σ, n, v, d, [], hasVert, s => by
    intro g anc pe hT hd hwf hends hcomp hpeA hnd hanc hv hw he hen hsz hsvk hsvd hσn hat hpath
    unfold CsOuts; trivial
  | σ, n, v, d, o :: rest, hasVert, s => by
    intro g anc pe hT hd hwf hends hcomp hpeA hnd hanc hv hw he hen hsz hsvk hsvd hσn hat hpath
    have hnd' := hnd
    rw [DfsOut.vertsList_eq, List.flatMap_cons, ← DfsOut.vertsList_eq] at hnd'
    have hndA : ∀ x, x ∈ anc → x ∈ v :: (o.verts ++ DfsOut.vertsList rest) → False :=
      fun x hx hx' => List.disjoint_of_nodup_append hnd' hx hx'
    have hndV : v ∉ o.verts ++ DfsOut.vertsList rest :=
      (List.nodup_cons.1 (List.nodup_append.1 hnd').2.1).1
    have hndOR : ∀ x, x ∈ o.verts → x ∈ DfsOut.vertsList rest → False :=
      fun x hx hx' =>
        List.disjoint_of_nodup_append (List.nodup_cons.1 (List.nodup_append.1 hnd').2.1).2 hx hx'
    have hmemR : ∀ x, x ∈ DfsOut.vertsList rest → x ∈ DfsOut.vertsList (o :: rest) := fun x hx => by
      rw [DfsOut.vertsList_eq, List.flatMap_cons, ← DfsOut.vertsList_eq]
      exact List.mem_append_right _ hx
    have hcompO : ∀ e', e' < g.ne → ∀ x, g.Inc e' x → x ∈ o.verts → subEdges o e' ∨ pe e' := by
      intro e' he' x hx hxo
      rcases hcomp e' he' x hx (mem_vertsList_of_verts (List.mem_cons_self ..) hxo) with h | h
      · obtain ⟨o', ho', hs⟩ := mem_subEdges_edgesList.1 h
        rcases List.mem_cons.1 ho' with heq | ho'
        · subst heq; exact .inl hs
        · exfalso
          rcases endsOut_wf g anc v o' (hwf o' (List.mem_cons_of_mem _ ho'))
            (hends o' (List.mem_cons_of_mem _ ho')) e' hs x hx with h' | rfl | h'
          · exact hndA x h' (List.mem_cons_of_mem _ (List.mem_append_left _ hxo))
          · exact hndV (List.mem_append_left _ hxo)
          · exact hndOR x hxo (mem_vertsList_of_verts ho' h')
      · exact .inr h
    have hcompR : ∀ e', e' < g.ne → ∀ x, g.Inc e' x → x ∈ DfsOut.vertsList rest →
        e' ∈ DfsOut.edgesList rest ∨ pe e' := by
      intro e' he' x hx hxr
      rcases hcomp e' he' x hx (hmemR x hxr) with h | h
      · obtain ⟨o', ho', hs⟩ := mem_subEdges_edgesList.1 h
        rcases List.mem_cons.1 ho' with heq | ho'
        · exfalso
          subst heq
          rcases endsOut_wf g anc v o' (hwf o' (List.mem_cons_self ..))
            (hends o' (List.mem_cons_self ..)) e' hs x hx with h' | rfl | h'
          · exact hndA x h' (List.mem_cons_of_mem _ (List.mem_append_right _ hxr))
          · exact hndV (List.mem_append_right _ hxr)
          · exact hndOR x h' hxr
        · exact .inl (mem_subEdges_edgesList.2 ⟨o', ho', hs⟩)
      · exact .inr h
    rw [DfsOut.vertsList_eq, List.flatMap_cons, ← DfsOut.vertsList_eq] at hnd hw
    rw [DfsOut.edgesList_eq, List.flatMap_cons, ← DfsOut.edgesList_eq] at he hen
    rw [DfsOut.edgePostorderList_cons] at hat hpath
    unfold CsOuts
    refine ⟨dsOut σ n v d o hasVert s g anc pe hT hd (hwf o (List.mem_cons_self ..)) (hends o (List.mem_cons_self ..))
      hcompO hpeA
      (List.Nodup.sublist (List.Sublist.append_left (List.Sublist.cons_cons v (List.sublist_append_left _ _)) _) hnd)
      hanc hv (fun w hw' => hw w (List.mem_append_left _ hw')) (fun e he' => he e (List.mem_append_left _ he'))
      (List.Nodup.sublist (List.sublist_append_left _ _) hen) hsz hsvk hsvd hσn hat.left
      (hpath.mono (by rw [List.length_append]; omega)), ?_⟩
    have hK := kOut v d o hasVert s g (d + 1) 0 s hT (Nat.le_refl _) (by omega) (vertItem_ne_zero v)
      (fun w _ => vertItem_ne_zero w) (fun e _ => edgeItem_ne_zero g e) Keep.refl
    refine wp_mono _ hK fun hv' s' hK' => ?_
    exact dsOuts σ (n + o.block.length) v d rest hv' s' g anc pe (hT.of_keep hK') hd
      (fun o' ho' => hwf o' (List.mem_cons_of_mem _ ho')) (fun o' ho' => hends o' (List.mem_cons_of_mem _ ho'))
      hcompR hpeA
      (List.Nodup.sublist (List.Sublist.append_left (List.Sublist.cons_cons v (List.sublist_append_right _ _)) _) hnd)
      hanc hv (fun w hw' => hw w (List.mem_append_right _ hw')) (fun e he' => he e (List.mem_append_right _ he'))
      (List.Nodup.sublist (List.sublist_append_right _ _) hen) (hK'.sv.trans hsz)
      (fun k hk => by rw [hK'.svlo k (by omega)]; exact hsvk k hk)
      (by rw [hK'.svlo d (Nat.lt_succ_self _)]; exact hsvd) hσn hat.right
      (hpath.mono (by rw [List.length_append, Nat.add_assoc]))

theorem dsOut : ∀ (σ : List Nat) (n v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState), DsOut σ n v d o hasVert s
  | σ, n, v, d, o, hasVert, s => by
    intro g anc pe hT hd hwf hends hcomp hpeA hnd hanc hv hw he hen hsz hsvk hsvd hσn hat hpath
    subst hd
    unfold CsOut
    have hK₀ := keep_walkOutPre (D := anc.length + 1) (j := 0) v anc.length o hasVert (Keep.refl (s := s))
    refine wp_mono _ hK₀ fun hv' s₁ hK₁ => ?_
    clear hK₀
    have hg₁ : s₁.g = g := hK₁.g.trans hT.g_eq
    have hsv₁ : ∀ k, k ≤ anc.length → s₁.stackVerts[k]! = (anc ++ [v])[k]! := fun k hk => by
      rw [hK₁.svlo k (by omega)]
      rcases Nat.lt_or_ge k anc.length with hk' | hk'
      · exact (getElem!_of_getElem? ((List.getElem?_append_left hk').trans (hsvk k hk'))).symm
      · have : k = anc.length := by omega
        subst this
        rw [hsvd, getElem!_concat_length']
    cases o with
    | back e dest cls =>
      simp only [DfsOut.WF] at hwf
      obtain ⟨i, hi, hcls⟩ := hwf
      simp only [DfsOut.block] at hat hpath
      have hil : i < anc.length + 1 := by
        have := (List.getElem?_eq_some_iff.1 hi).1; simpa using this
      subst hcls
      rw [classify_back hil] at hends ⊢
      refine ⟨hat.singleton, hat.block_infix, by rw [hg₁]; exact hends, ?_, ?_, ?_, ?_, ?_⟩
      · intro k hk
        obtain ⟨e', he', hm, hpe⟩ := hpath k (by simp; omega)
        refine ⟨e', hg₁ ▸ he', by simp at hm; omega, ?_⟩
        rw [hsv₁ k (by omega), hsv₁ (k + 1) (by omega), hg₁]; exact hpe
      · intro ht
        simp only [DfsOut.cls] at ht
        split at ht <;> simp [OutClass.isTree] at ht
      · intro ht
        simp only [DfsOut.cls] at ht
        split at ht <;> simp [OutClass.isTree] at ht
      · show dest < s₁.g.nv
        rw [hg₁]
        rcases List.mem_append.1 (List.mem_of_getElem? hi) with h | h
        · exact hanc _ h
        · exact (List.mem_singleton.1 h) ▸ hv
      · intro _ hge
        simp only [DfsOut.cls, DfsOut.dest] at hge ⊢
        split at hge
        · rename_i hid
          subst hid
          rw [List.getElem?_concat_length] at hi
          exact (Option.some.inj hi).symm
        · simp [OutClass.lowval] at hge; omega
    | tree e cls child =>
      simp only [DfsOut.WF] at hwf
      obtain ⟨hwfc, hcls⟩ := hwf
      simp only [DfsOut.Ends] at hends
      obtain ⟨hpe, hendsc⟩ := hends
      simp only [DfsOut.verts] at hnd hw
      simp only [DfsOut.edges, List.mem_cons, forall_eq_or_imp] at he
      simp only [DfsOut.edges] at hen
      simp only [DfsOut.block] at hat hpath
      have htr : cls.isTree = true := by
        cases hc : cls with
        | bridge => rfl
        | component => rfl
        | selfLoop => rw [hc] at hcls; unfold classify at hcls; split at hcls <;> split at hcls <;> cases hcls
        | ret lv kind =>
          rw [hc] at hcls
          obtain ⟨_, _, hk⟩ := classify_eq_ret_iff_tree.1 hcls.symm
          cases kind with
          | backEdge => split at hk <;> cases hk
          | _ => rfl
      simp only [wp_modify]
      have hT₁ := hT.of_keep hK₁
      have hT₂ : Types g { s₁ with firstOccurrence := s₁.firstOccurrence.set! anc.length s₁.g.ne } :=
        ⟨hT₁.g_eq, hT₁.size, hT₁.root, hT₁.vert, hT₁.edge⟩
      have hnd' : (anc ++ [v] ++ child.verts).Nodup := by simpa using hnd
      have hpos : σ.idxOf e = n + child.edgePostorder.length := by
        have h1 := hat.right.singleton
        have h2 := (List.getElem?_eq_some_iff.1 h1).1
        rw [← getElem!_of_getElem? h1]; exact idxOf_getElem! hσn h2
      have hpc : ∀ w outs, child = .node w outs →
          AncPath g σ (n + child.edgePostorder.length) (anc ++ [v] ++ [w]) := by
        rintro w outs rfl k hk
        simp only [List.length_append, List.length_singleton] at hk
        rcases Nat.lt_or_ge (k + 1) (anc.length + 1) with hk' | hk'
        · obtain ⟨e', he', hm, hpe'⟩ := hpath k (by simp; omega)
          refine ⟨e', he', by simp at hm; omega, ?_⟩
          rw [getElem!_append_left' (l₁ := anc ++ [v]) (l₂ := [w]) (k := k) (by simp; omega),
            getElem!_append_left' (l₁ := anc ++ [v]) (l₂ := [w]) (k := k + 1) (by simp; omega)]
          exact hpe'
        · have : k = anc.length := by omega
          subst this
          refine ⟨e, he.1, hpos.symm.le, ?_⟩
          have hl : anc.length + 1 = (anc ++ [v]).length := by simp
          rw [getElem!_append_left' (l₁ := anc ++ [v]) (l₂ := [w]) (by simp), getElem!_concat_length', hl,
            getElem!_concat_length']
          simp only [DfsTree.v] at hpe
          exact PairEq.flip hpe
      have hcompC : ∀ e', e' < g.ne → ∀ x, g.Inc e' x → x ∈ child.verts →
          e' ∈ child.edges ∨ (e' = e ∨ pe e') := by
        intro e' he' x hx hxc
        rcases hcomp e' he' x hx hxc with h | h
        · rcases h with h | h
          · exact .inr (.inl h)
          · exact .inl h
        · exact .inr (.inr h)
      have hpeC : ∀ e', (e' = e ∨ pe e') → ∀ x, g.Inc e' x → x ∈ anc ++ [v] ++ [child.v] := by
        intro e' h x hx
        rcases h with rfl | h
        · rcases eq_of_inc_pairEq hpe hx with h | h <;> rw [h]
          · exact List.mem_append_right _ (List.mem_singleton_self _)
          · exact List.mem_append_left _ (List.mem_append_right _ (List.mem_singleton_self _))
        · exact List.mem_append_left _ (hpeA e' h x hx)
      have hcv : child.v ∈ child.verts := by
        obtain ⟨w, outs⟩ := child; exact List.mem_cons_self ..
      refine ⟨dsTree σ n child (anc.length + 1) _ g (anc ++ [v]) (fun e' => e' = e ∨ pe e') hT₂ (by simp)
        hwfc hendsc hcompC hpeC hnd'
        (fun a ha => by
          rcases List.mem_append.1 ha with ha | ha
          · exact hanc a ha
          · exact (List.mem_singleton.1 ha) ▸ hv)
        hw he.2 (List.nodup_cons.1 hen).2 (hK₁.sv.trans hsz)
        (fun k hk => by
          rcases Nat.lt_or_ge k anc.length with hk' | hk'
          · rw [List.getElem?_append_left hk', hsvk k hk']
            show _ = some s₁.stackVerts[k]!
            rw [hK₁.svlo k (by omega)]
          · have : k = anc.length := by omega
            subst this
            rw [List.getElem?_concat_length]
            show _ = some s₁.stackVerts[anc.length]!
            rw [hK₁.svlo _ (Nat.lt_succ_self _), hsvd])
        hσn hat.left hpc, ?_⟩
      have hKc := kTree child (anc.length + 1) _ g (anc.length + 1) 0 _ hT₂ (Nat.le_refl _) (by omega)
        (fun w _ => vertItem_ne_zero w) (fun e' _ => edgeItem_ne_zero g e') Keep.refl
      refine wp_mono _ hKc fun _ s₃ hK₃ => ?_
      have hg₃ : s₃.g = g := hK₃.g.trans hg₁
      have hsv₃ : ∀ k, k ≤ anc.length → s₃.stackVerts[k]! = (anc ++ [v])[k]! := fun k hk => by
        rw [hK₃.svlo k (by omega)]; exact hsv₁ k hk
      refine ⟨hat.right.singleton, hat.block_infix, by rw [hg₃]; simp only [DfsOut.Ends]; exact ⟨hpe, hendsc⟩, ?_, ?_, ?_, ?_, ?_⟩
      · intro k hk
        obtain ⟨e', he', hm, hpe'⟩ := hpath k (by simp; omega)
        refine ⟨e', hg₃ ▸ he', by simp at hm; omega, ?_⟩
        rw [hsv₃ k (by omega), hsv₃ (k + 1) (by omega), hg₃]; exact hpe'
      · intro _ hlt
        simp only [DfsOut.cls] at hlt
        rw [hcls] at hlt
        obtain ⟨w, outs⟩ := child
        simp only [DfsTree.Ends] at hendsc
        cases outs with
        | nil =>
          exfalso
          simp [DfsTree.retDepths, DfsOut.retDepthsList, low2, classify, OutClass.lowval] at hlt
        | cons o' rest' =>
          have hmem : o'.e ∈ (DfsTree.node w (o' :: rest')).edges := by
            cases o' <;> simp [DfsTree.edges, DfsOut.edgesList, DfsOut.e]
          refine ⟨o'.e, by rw [hg₃]; exact he.2 _ hmem, fun h => (List.nodup_cons.1 hen).1 (by simp only [DfsOut.e] at h; exact h ▸ hmem), ?_⟩
          show s₃.g.Inc o'.e w
          rw [hg₃]
          have hE := hendsc o' (List.mem_cons_self ..)
          cases o' with
          | back e₀ d₀ c₀ => simp only [DfsOut.Ends] at hE; exact (Graph.inc_of_pairEq hE).1
          | tree e₀ c₀ ch => simp only [DfsOut.Ends] at hE; exact (Graph.inc_of_pairEq hE.1).2
      · intro _ e' he' hinc
        rw [hg₃] at he' hinc
        rcases hcomp e' he' _ hinc hcv with h | h
        · exact h
        · exact (List.disjoint_of_nodup_append hnd' (hpeA e' h _ hinc) hcv).elim
      · show child.v < s₃.g.nv
        rw [hg₃]
        obtain ⟨w, outs⟩ := child
        exact hw w (by simp [DfsTree.verts])
      · intro h
        simp only [DfsOut.cls] at h
        rw [h] at htr; cases htr
end

/-- `CsTree` on the real walk: a DFS tree below the ancestor path `anc`, scheduled at `n`. -/
theorem csTree_of_dfs {g : Graph} {anc : List Nat} (t : DfsTree) (d : Nat) (s : WalkState)
    (hT : Types g s) (hd : d = anc.length) (hwf : t.WF anc) (hends : t.Ends g)
    (hcomp : ∀ e, e < g.ne → ∀ x, g.Inc e x → x ∈ t.verts → e ∈ t.edges)
    (hnd : (anc ++ t.verts).Nodup) (hanc : ∀ a ∈ anc, a < g.nv) (hvlt : ∀ v ∈ t.verts, v < g.nv)
    (helt : ∀ e ∈ t.edges, e < g.ne) (hen : t.edges.Nodup) (hsz : s.stackVerts.size = g.nv)
    (hsvk : ∀ k, k < d → anc[k]? = some s.stackVerts[k]!) (hσn : σ.Nodup) (hat : PostAt σ n t.edgePostorder)
    (hpath : ∀ v outs, t = .node v outs → AncPath g σ (n + t.edgePostorder.length) (anc ++ [v])) :
    CsTree σ n t d s :=
  dsTree σ n t d s g anc (fun _ => False) hT hd hwf hends
    (fun e he x hx hxt => .inl (hcomp e he x hx hxt)) (fun _ h => h.elim)
    hnd hanc hvlt helt hen hsz hsvk hσn hat hpath

theorem forest_closeInv {g : Graph} {σ : List Nat} (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < g.ne) :
    ∀ forest pre n s, RootState g pre s → s.RangesInv σ n 0 → ForestOK g (pre ++ forest) →
      (∀ t ∈ forest, t.WF []) → (∀ t ∈ forest, t.Ends g) →
      (∀ t ∈ forest, ∀ e, e < g.ne → ∀ x, g.Inc e x → x ∈ t.verts → e ∈ t.edges) →
      RootsCover σ n forest s →
      PostAt σ n (edgePostorderForest forest) → s.CloseInv →
      wp (walkForest forest) (fun _ s' => s'.CloseInv) s
  | [], _, _, _, _, _, _, _, _, _, _, _, hcl => by simpa [walkForest, wp_pure] using hcl
  | t :: rest, pre, n, s, h, hr, hf, hwf, hends, hcomp, hc, hat, hcl => by
    have hb := h.book hf (hwf t (by simp)) (hends t (by simp)) (hcomp t (by simp))
    have hg := gbTree t 0 s hb
    have hi : ∀ v outs, t = .node v outs →
        ({ s with stackVerts := s.stackVerts.set! 0 v } : WalkState).RangesInv σ n 0 :=
      fun v _ _ => hr.stackVerts_of_nil h.tstack _
    have hσ' : ∀ e ∈ σ, e < s.g.ne := by rwa [h.g_eq]
    have hfront := walkTree_frontiers t 0 s (fun v outs ht => (hi v outs ht).inv) h.shape hg hb
    have hat' : PostAt σ n (t.edgePostorder ++ edgePostorderForest rest) := hat
    have hsched := scheduleTree σ n t 0 s hi h.shape hnd hσ' hg hb hfront hc.1 hat'.left
    have hrg := rgTree σ n t 0 s hi h.shape hnd hσ' hg hb hsched
    have hcs : CsTree σ n t 0 s :=
      csTree_of_dfs (anc := []) t 0 s h.place.types rfl (hwf t (by simp)) (hends t (by simp))
        (hcomp t (by simp))
        (by simpa using (RootState.hvn hf).1) (by simp) (RootState.hvlt hf) (RootState.helt hf)
        (RootState.hen hf).1 h.sv (fun k hk => absurd hk (Nat.not_lt_zero _)) hnd hat'.left
        (fun v outs _ k hk => absurd hk (by simp))
    have hcc := ccTree σ n t 0 s hi h.shape hnd hσ' hg hb hfront hsched hcs hcl
    have hnv : 0 < s.g.nv := by
      obtain ⟨v, outs⟩ := t
      have hv := RootState.hvlt hf v (by simp [DfsTree.verts])
      rw [h.g_eq]; omega
    have hk : wp (walkTree t 0) (fun _ s' => RootOK s') s :=
      walkTree_rootOK t s (fun v outs ht => (hi v outs ht).inv) h.shape hg hb h.tstack
      (by rw [h.sd, ← h.g_eq]; exact hnv)
    have hp := (walk_place_aux g).1 t 0 _ _ s h.place (RootState.hvlt hf) (RootState.helt hf)
      (RootState.hvn hf).1 (RootState.hen hf).1 (RootState.hPv hf) (RootState.hPe hf)
    have hst := h.step hf (hwf t (by simp)) (hends t (by simp)) (hcomp t (by simp))
    show wp ((walkTree t 0 >>= fun _ => popTstack >>= fun top =>
      modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 }) >>= fun _ => walkForest rest) _ s
    simp only [wp_bind]
    refine wp_imp (wp_of_forall fun _ s₁ ⟨hrs, hk, hp, hst, hc, hcl⟩ => ?_)
      (wp_and hrg (wp_and hk (wp_and hp (wp_and hst (wp_and hc.2 hcc)))))
    have hrpop := hrs.1.root_append hrs.2.1 hk (noParent_of_cnt_eq_zero hp.root)
    have hp' : s₁.Place s₁.g (Pushed g (Pushed g (fun _ => False) (pre.flatMap DfsTree.verts)
        (pre.flatMap DfsTree.edges)) t.verts t.edges) (fun _ => False) := by
      rw [hp.g_eq]; exact hp
    have hclpop := rootAppend_closeInv hcl hp' hk
    exact wp_mono (popTstack >>= fun top => modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 })
      (wp_and hrpop (wp_and hst (wp_and hc hclpop))) fun _ s₂ ⟨hr₂, hs₂, hc₂, hcl₂⟩ =>
      forest_closeInv hnd hσ rest (pre ++ [t]) (n + t.edgePostorder.length) s₂ hs₂ hr₂
        (by simpa using hf) (fun t' ht' => hwf t' (by simp [ht']))
        (fun t' ht' => hends t' (by simp [ht'])) (fun t' ht' => hcomp t' (by simp [ht'])) hc₂ hat'.right hcl₂

/-- `CloseInv` for the walk of a DFS forest, from `RootsCover` (the range-side schedule). -/
theorem walk_closeInv' (g : Graph) (tern : Bool) (forest : List DfsTree)
    (hf : ForestOK g forest) (hwf : ∀ t ∈ forest, t.WF []) (hends : ∀ t ∈ forest, t.Ends g)
    (hecov : ∀ e, e < g.ne → e ∈ forest.flatMap DfsTree.edges)
    (hc : RootsCover (edgePostorderForest forest) 0 forest (WalkState.init g tern)) :
    (g.walk tern forest).CloseInv := by
  have hnd := DfsData.edgePostorderForest_perm.nodup_iff.2 hf.edges_nodup
  have hσ : ∀ e ∈ edgePostorderForest forest, e < g.ne :=
    fun e he => hf.edges_lt e (DfsData.edgePostorderForest_perm.subset he)
  have h := forest_closeInv hnd hσ forest [] 0 (WalkState.init g tern)
    (rootState_init g tern) (init_rangesInv g tern hnd) (by simpa using hf) hwf hends
    (fun t ht => comp_of_forest hf hwf hends hecov ht) hc
    ⟨[], [], rfl, by simp⟩ (init_closeInv g tern)
  simpa only [wp, Graph.walk] using h

end WalkState
end Spqr
