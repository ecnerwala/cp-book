import Spqr.WalkCover
import Spqr.EarFrontier

/-!
# Pop-site conditions at `finishEdge` from the ear contract

`finishEdge_sides` establishes `FinishSides` — the `MergeOK`/`CloseOK`/`UnwrapOK`/`BoundaryOK`
site conditions that `walk_sides` needs at a `finishEdge` — from `FinishGuards`, `FinishBook`
(whose ear contract gives `FinishOk` through `finishOk_of_guards`) and the schedule frontier
(`finishEdge_frontier`, for the loop-2/3 stack-length bounds). `walkTree_sides` propagates it
through the walk (`SidesTree`), like `walkTree_frontiers`.
-/

namespace Spqr
open WalkM
namespace WalkState

variable {s : WalkState}

/-- `LoopSides` from a per-iterate statement: `B` at every iterate reached with the condition true. -/
theorem loopSides_of_iter {B : WalkState → Prop} {cond : WalkM Bool} {body : WalkM Unit}
    (hc : ∀ s, (cond.run s).2 = s) :
    ∀ (n : Nat) (s : WalkState),
      (∀ k, (∀ j, j < k → (cond.run (iter body j s)).1 = true) →
        (cond.run (iter body k s)).1 = true → B (iter body k s)) →
      LoopSides B cond body n s
  | 0, _, _ => trivial
  | n + 1, s, h => by
    show (cond.run s).1 = true → B (cond.run s).2 ∧ LoopSides B cond body n (body.run (cond.run s).2).2
    rw [hc s]
    intro hb
    refine ⟨h 0 (fun j hj => absurd hj (Nat.not_lt_zero _)) hb, loopSides_of_iter hc n _ fun k hj hk => ?_⟩
    exact h (k + 1) (fun j hj' => by
      cases j with
      | zero => exact hb
      | succ j => exact hj j (Nat.lt_of_succ_lt_succ hj')) hk

theorem exists_cons_cons_of_two (h : 2 ≤ s.tstack.length) : ∃ a b rest, s.tstack = a :: b :: rest := by
  match hts : s.tstack with
  | a :: b :: rest => exact ⟨a, b, rest, rfl⟩
  | [] | [_] => rw [hts] at h; simp at h

theorem length_after_merge (h : 2 ≤ s.tstack.length) :
    (mergeTstackTops.run s).2.tstack.length = s.tstack.length - 1 := by
  obtain ⟨a, b, rest, hs⟩ := exists_cons_cons_of_two h
  rw [mergeTstackTops_run_eq s a b rest hs]; simp [hs]

theorem length_after_unwrap (ty : NodeType) (h : 2 ≤ s.tstack.length) :
    ((maybeUnwrapNxt ty).run s).2.tstack.length = s.tstack.length := by
  obtain ⟨a, b, rest, hs⟩ := exists_cons_cons_of_two h
  rw [maybeUnwrapNxt_run_eq ty s a b rest hs _ rfl _ rfl]
  split_ifs <;> simp [run_allocItem, hs]

/-- Merging while more than `m ≥ 1` entries remain never drops below `m`. -/
theorem iter_merge_length_ge (m : Nat) (hm : 1 ≤ m) : ∀ (k : Nat) (s : WalkState), m ≤ s.tstack.length →
    (∀ j, j < k → m < (iter mergeTstackTops j s).tstack.length) →
    m ≤ (iter mergeTstackTops k s).tstack.length
  | 0, _, h0, _ => h0
  | k + 1, s, _, hj => by
    have h1 : m < s.tstack.length := hj 0 (Nat.succ_pos _)
    have hlen := length_after_merge (s := s) (by omega)
    exact iter_merge_length_ge m hm k _ (by rw [hlen]; omega)
      fun j hj' => hj (j + 1) (Nat.succ_lt_succ hj')

theorem closeOK_of_finishTopOk {D : Nat} (h : FinishTopOk D s) : CloseOK s := by
  obtain ⟨t, rest, hts⟩ : ∃ t rest, s.tstack = t :: rest := by
    cases hts : s.tstack with
    | nil => exact absurd hts h.nonempty
    | cons t rest => exact ⟨t, rest, rfl⟩
  refine ⟨t, rest, hts, ?_⟩
  have := h.side
  rw [curE, hts, List.head!_cons] at this
  exact this

theorem unwrapOK_of_ok {ty : NodeType} (hty : ty ≠ .F) (h : UnwrapOk ty s) : UnwrapOK s ty := by
  intro hne a t rest hts x xs hx htx
  have hne' : ¬ (ty = .R ∨ s.ternarize = true) := by
    simp only [Bool.or_eq_false_iff, beq_eq_false_iff_ne, ne_eq] at hne
    rintro (h | h)
    · exact hne.1 h
    · rw [h] at hne; exact Bool.noConfusion hne.2
  have hxlt : x < s.items.size := by
    by_contra hge
    have : Items.type s.items x = .F := by
      simp [Items.type, Array.getElem?_eq_none_iff.2 (Nat.le_of_not_lt hge)]
    exact hty (htx ▸ this)
  have hhead : nxtHead s = x := by
    simp only [nxtHead, nxtE, nxtDir, hts, List.tail_cons, List.head!_cons, hx]
  have hat := h.unwrap hne' (by
    rw [hhead, getElem!_pos s.items x hxlt, ← Items.type_eq_getElem hxlt]; exact htx)
  have hs1 := hat.single
  have hs2 := hat.side
  rw [hhead] at hs1
  simp only [nxtE, nxtDir, hts, List.tail_cons, List.head!_cons] at hs1 hs2
  rw [hx] at hs1
  exact ⟨(List.cons.inj hs1).2, hs2⟩

theorem umc_of_ok {D : Nat} {ty : NodeType} (hty : ty ≠ .F) (hu : UnwrapOk ty s)
    (hc : CloseTwoOk D (after (maybeUnwrapNxt ty) s)) : UnwrapMergeClose ty s := by
  refine ⟨unwrapOK_of_ok hty hu, ?_⟩
  show MergeOK (after (maybeUnwrapNxt ty) s) ∧ CloseOK (after mergeTstackTops (after (maybeUnwrapNxt ty) s))
  exact ⟨by show 2 ≤ _; rw [after, length_after_unwrap ty hu.two]; exact hu.two,
    closeOK_of_finishTopOk hc.finish⟩

theorem finishPSides_of_ok {D curV lowval : Nat} {isType1 : Bool} (h : FinishPOk D curV lowval isType1 s) :
    FinishPSides curV lowval isType1 s := by
  show result (condP curV lowval isType1) s = true → UnwrapMergeClose .P s
  intro hc
  exact umc_of_ok (by decide) (h.ok hc).1 (h.ok hc).2

theorem after_finishP_of_not {curV lowval : Nat} {isType1 : Bool}
    (h : result (condP curV lowval isType1) s = false) : after (finishP curV lowval isType1) s = s := by
  have h' : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == curV) &&
      (s.tstack.tail.head!.topDepth == lowval)) = false := h
  simp only [after, Spqr.finishP, WalkM.run_bind, run_condP, h', Bool.false_eq_true, ↓reduceIte, WalkM.pure_run]

/-- `FinishRestSides` from `FinishRestOk`; the vertex push (`!hasVert`, `!isSingle`) needs a nonempty
stack, which `finishP` leaves alone when its condition fails. -/
theorem finishRestSides_of_ok {D curV d lowval : Nat} {isType1 hasVert isSingle : Bool}
    (h : FinishRestOk D curV d lowval isType1 hasVert isSingle s)
    (htail : hasVert = false → isSingle = false →
      s.tstack ≠ [] ∧ result (condP curV lowval isType1) s = false) :
    FinishRestSides curV d lowval isType1 hasVert isSingle s := by
  refine ⟨finishPSides_of_ok h.p, ?_⟩
  show FinishTailSides curV d hasVert isSingle (after (finishP curV lowval isType1) s)
  intro hv hsg
  obtain ⟨hne, hc⟩ := htail hv hsg
  rw [after_finishP_of_not hc]
  obtain ⟨t, rest, hts⟩ := List.exists_cons_of_ne_nil hne
  show 2 ≤ (_ :: s.tstack).length
  rw [hts]; simp

/-- Site conditions of `finishEdge` from the guards, the bookkeeping (ear contract) and the
frontier. -/
theorem finishEdge_sides {D curV d : Nat} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}
    (hg : FinishGuards d o origTstack hasVert s) (hb : FinishBook curV d o origTstack hasVert s)
    (hi : s.Inv' D) (hs : Shape s) (hD : D = if o.cls.isTree then d + 1 else d) :
    FinishSides curV d o origTstack hasVert s := by
  obtain ⟨sub, base, hlen, hE⟩ := hb.ear
  refine ⟨fun hge ht => ?_, fun hlt => ?_⟩
  · -- boundary: `bd_side` + freshness of `vertItem curV`
    have hside := hE.bd_side ht hge
    have hnb : ∀ t ∈ s.tstack, ∀ c ∈ t.spans.1 ++ t.spans.2, ¬ Items.Below s.items c (vertItem curV) := by
      intro t ht' c hc hbel
      have hc' := Items.Below.eq_of_no_parent hE.v_root hbel
      subst hc'
      have := (hE.vert_free t ht' hc).1
      rw [hE.bd_noVert hge] at this
      exact Bool.noConfusion this
    by_cases h1 : o.cls.lowval d = d + 1
    · rw [if_pos h1] at hside
      simp only [h1, beq_self_eq_true, ↓reduceIte]
      intro t ht'
      have ht'' := List.mem_of_mem_head? ht'
      exact ⟨hside t ht', fun c hc => hnb t ht'' c (List.mem_append_right _ hc)⟩
    · rw [if_neg h1] at hside
      simp only [beq_eq_false_iff_ne.2 h1, Bool.false_eq_true, ↓reduceIte]
      refine ⟨fun b hb' => ⟨hside.1 b hb', fun c hc =>
        hnb b (List.mem_of_mem_head? hb') c (List.mem_append_left _ hc)⟩, fun t ht' => ?_⟩
      have ht'' := List.mem_of_mem_tail (List.mem_of_mem_head? ht')
      exact ⟨hside.2 t ht', fun c hc => hnb t ht'' c (List.mem_append_right _ hc)⟩
  · obtain ⟨lv, kind, ho, hl⟩ := ret_of_lowval_lt hlt
    have hlv : o.cls.lowval d = lv := by rw [ho]; rfl
    have hok := finishOk_of_guards ho hl hg hE hlen hi hs hD hb.v_lt hb.e_lt hb.q (hb.ends lv kind ho) hb.vert
    have hfr := finishEdge_frontier (D := D) hb hi hs hD
    show (o.cls.isTree = true → FinishTreeSides curV d o origTstack hasVert s.stackDir[d]! (feS₀ d o s)) ∧
      (o.cls.isTree = false → FinishBackSides curV d o hasVert (feS₀ d o s))
    unfold FinishTreeSides FinishBackSides
    rw [hlv]
    refine ⟨fun ht => ⟨?_, ?_⟩, fun ht' => ?_⟩
    · -- loop 1
      show LoopSides (Loop1Sides d s.stackDir[d]!) (loop1Cond d) (loop1Body d s.stackDir[d]!)
        (ceS₁ o.dest d o.e (feS₀ d o s)).tstack.length (ceS₁ o.dest d o.e (feS₀ d o s))
      refine loopSides_of_iter (fun _ => rfl) _ _ fun k hj hk => ?_
      have hbody := (hok.ears ht).body k fun j hj' => by
        rcases Nat.lt_or_eq_of_le hj' with h | h
        · exact hj j h
        · subst h; exact hk
      have hres := loop1Type_result d s.stackDir[d]! (iter (loop1Body d s.stackDir[d]!) k (ceS₁ o.dest d o.e (feS₀ d o s)))
      refine ⟨fun _ => ((result_loop1Cond_iff d _).1 hk).1, ?_⟩
      show UnwrapMergeClose (l1Ty d s.stackDir[d]! _) (l1S₁ d s.stackDir[d]! _)
      exact umc_of_ok (fun h => hres (by simp only [l1Ty] at h; simp only [h, List.mem_cons, true_or]))
        hbody.unwrap hbody.close
    · -- loops 2, 3 and the rest
      show MergeLateSides d (feS₁ d o s) ∧
        (hasVert = true → CloseVertSides curV s.stackDir[d]! o.cls.isType1 origTstack (feSingle d o s)
          (fun b => FinishRestSides curV d lv o.cls.isType1 hasVert b) (feS₂ d o s)) ∧
        (hasVert = false → FinishRestSides curV d lv o.cls.isType1 hasVert (feSingle d o s) (feS₂ d o s))
      refine ⟨fun _ => loopSides_of_iter (fun _ => rfl) _ _ fun k hj hk => ?_, fun hv => ?_, fun hv => ?_⟩
      · have h2 := hfr.loop2 ht hlt k hj
        dsimp only at h2
        show 2 ≤ _
        have := h2.2 hk
        omega
      · -- closeVert
        have hg' : Inv1 d (feS₀ d o s).tstack ∧ (Inv2 (feS₁ d o s).firstOccurrence[d]! (feS₁ d o s).tstack ∧
            (hasVert = true → 3 ≤ (feS₂ d o s).tstack.length ∧
              (o.cls.isType1 = false → origTstack + 3 ≤ (feS₂ d o s).tstack.length))) := hg.2 hlt ht
        have hg3 := hg'.2.2 hv
        have hcv := hok.vert ht hv
        have hrest := hok.rest_vert ht hv
        subst hv
        cases h1 : o.cls.isType1
        · have hS₃ : feS₃ curV d o origTstack s = cvS₅ curV s.stackDir[d]! false origTstack (feSingle d o s) (feS₂ d o s) := by
            simp only [feS₃, h1, closeVert', cvS₅, cvS₄, cvS₃, cvS₂, cvS₁, after, vertPre,
              vertUnwrap, vertFinish, Bool.not_false, Bool.false_eq_true, ↓reduceIte, WalkM.run_bind,
              WalkM.pure_run, run_tstackSize]
          have hB₃ : feB₃ curV d o origTstack s = false := by
            simp only [feB₃, h1, closeVert', result, vertPre, vertUnwrap, vertFinish, Bool.not_false,
              Bool.false_eq_true, ↓reduceIte, WalkM.run_bind, WalkM.pure_run, run_tstackSize]
          have hpre : cvS₁ false origTstack (feSingle d o s) (feS₂ d o s) =
              after (loop (feS₂ d o s).tstack.length (loop3Cond origTstack) mergeTstackTops) (feS₂ d o s) := by
            simp only [cvS₁, after, vertPre, Bool.not_false, ↓reduceIte, WalkM.run_bind, WalkM.pure_run,
              run_tstackSize]
          rw [h1] at hcv hrest
          rw [hS₃, hB₃] at hrest
          refine ⟨fun _ => ⟨loopSides_of_iter (fun _ => rfl) _ _ fun k hj hk => ?_, ?_⟩, fun h => nomatch h⟩
          · have h3 := hfr.loop3 ht hlt h1 k hj
            dsimp only at h3
            show 2 ≤ _
            have := h3.2 hk
            omega
          · show CloseVertTailSides curV s.stackDir[d]! none
              (FinishRestSides curV d lv false true false)
              (after (loop (feS₂ d o s).tstack.length (loop3Cond origTstack) mergeTstackTops) (feS₂ d o s))
            rw [← hpre]
            have hL : origTstack + 3 ≤ (cvS₁ false origTstack (feSingle d o s) (feS₂ d o s)).tstack.length := by
              rw [hpre]
              obtain ⟨k, -, hk, hj, -⟩ := loop_run_iter (feS₂ d o s).tstack.length (loop3Cond origTstack)
                mergeTstackTops (feS₂ d o s) (fun _ => rfl)
              show origTstack + 3 ≤ ((loop _ _ _).run _).2.tstack.length
              rw [hk]
              exact iter_merge_length_ge _ (by omega) k _ (hg3.2 h1) fun j hj' => by
                have := hj j hj'
                rw [run_loop3Cond] at this
                simpa using this
            refine ⟨by show 2 ≤ _; omega, ?_⟩
            show MergeOK (after mergeTstackTops _) ∧
              FinishRestSides curV d lv false true false (cvS₅ curV s.stackDir[d]! false origTstack (feSingle d o s) (feS₂ d o s))
            refine ⟨by show 2 ≤ _; rw [after, length_after_merge (by omega)]; omega, ?_⟩
            exact finishRestSides_of_ok hrest fun h => nomatch h
        · have hS₃ : feS₃ curV d o origTstack s =
              after (finishTstackTop ((maybeUnwrapNxt (if feSingle d o s then .S else .R)).run (feS₂ d o s)).1)
                (cvS₅ curV s.stackDir[d]! true origTstack (feSingle d o s) (feS₂ d o s)) := by
            simp only [feS₃, h1, closeVert', cvS₅, cvS₄, cvS₃, cvS₂, cvS₁, cvB₁, after, result, vertPre,
              vertUnwrap, vertFinish, Bool.not_true, Bool.false_eq_true, ↓reduceIte, WalkM.run_bind,
              WalkM.pure_run, WalkM.map_run]
            rfl
          have hB₃ : feB₃ curV d o origTstack s = true := by
            simp only [feB₃, h1, closeVert', result, vertPre, vertUnwrap, vertFinish, Bool.not_true,
              Bool.false_eq_true, ↓reduceIte, WalkM.run_bind, WalkM.pure_run, WalkM.map_run]
          have hU : cvS₂ true origTstack (feSingle d o s) (feS₂ d o s) =
              after (maybeUnwrapNxt (if feSingle d o s then .S else .R)) (feS₂ d o s) := by
            simp only [cvS₂, cvS₁, cvB₁, after, result, vertPre, vertUnwrap, Bool.not_true, Bool.false_eq_true,
              ↓reduceIte, WalkM.pure_run, WalkM.map_run]
          rw [h1] at hcv hrest
          rw [hS₃, hB₃] at hrest
          refine ⟨(fun h => nomatch h), fun _ => ⟨unwrapOK_of_ok (by split <;> decide) (hcv.unwrap rfl), ?_⟩⟩
          show CloseVertTailSides curV s.stackDir[d]! (some _) (FinishRestSides curV d lv true true true)
            (after (maybeUnwrapNxt (if feSingle d o s then .S else .R)) (feS₂ d o s))
          rw [← hU]
          have hlenU : (cvS₂ true origTstack (feSingle d o s) (feS₂ d o s)).tstack.length = (feS₂ d o s).tstack.length := by
            rw [hU, after, length_after_unwrap _ (by omega)]
          refine ⟨by show 2 ≤ _; omega, ?_⟩
          show MergeOK (after mergeTstackTops _) ∧
            (CloseOK (cvS₅ curV s.stackDir[d]! true origTstack (feSingle d o s) (feS₂ d o s)) ∧
              FinishRestSides curV d lv true true true
                (after (finishTstackTop _) (cvS₅ curV s.stackDir[d]! true origTstack (feSingle d o s) (feS₂ d o s))))
          refine ⟨by show 2 ≤ _; rw [after, length_after_merge (by omega)]; omega, ?_, ?_⟩
          · exact closeOK_of_finishTopOk (hcv.finish rfl)
          · exact finishRestSides_of_ok hrest fun h => nomatch h
      · refine finishRestSides_of_ok (hok.rest_tree ht hv) fun _ _ => ?_
        obtain ⟨c, mid, py, vy, hts, -⟩ := hE.loops ht hlt
        exact ⟨(by rw [hts]; exact List.cons_ne_nil _ _), ear_condP_tree ho hl ht hg hE hi hs hv⟩
    · show FinishRestSides curV d lv o.cls.isType1 hasVert true (feBack curV lv d o s)
      exact finishRestSides_of_ok (hok.rest_back ht') fun _ h => nomatch h

/-! `SidesTree` by the walk induction (the shape of `walkTree_frontiers`): `Inv'`/`Shape` are threaded
to each `finishEdge` site, where `finishEdge_sides` applies. -/

abbrev SdTree (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  (∀ v outs, t = .node v outs → ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).Inv' d) →
  Shape s → GuardsTree t d s → BookTree t d s → SidesTree t d s

abbrev SdOuts (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  s.Inv' d → Shape s → GuardsOuts v d outs hasVert s → BookOuts v d outs hasVert s →
  SidesOuts v d outs hasVert s

abbrev SdOut (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  s.Inv' d → Shape s → GuardsOut v d o hasVert s → BookOut v d o hasVert s →
  SidesOut v d o hasVert s

mutual
theorem sdTree : ∀ (t : DfsTree) (d : Nat) (s : WalkState), SdTree t d s
  | .node v outs, d, s => fun hi hs hg hb => by
    unfold SidesTree; unfold GuardsTree at hg; unfold BookTree at hb
    exact sdOuts v d outs false _ (hi v outs rfl) hs.frame' hg hb

theorem sdOuts : ∀ (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState),
    SdOuts v d outs hasVert s
  | v, d, [], hasVert, s => fun _ _ _ _ => by
    unfold SidesOuts; trivial
  | v, d, o :: rest, hasVert, s => fun hi hs hg hb => by
    unfold SidesOuts; unfold GuardsOuts at hg; unfold BookOuts at hb
    refine ⟨sdOut v d o hasVert s hi hs hg.1 hb.1, ?_⟩
    exact wp_imp (wp_imp (wp_imp (wp_of_forall fun hv' s' ⟨hi', hs'⟩ hg' hb' =>
      sdOuts v d rest hv' s' hi' hs' hg' hb') (invOut v d o hasVert s hi hs hg.1 hb.1)) hg.2) hb.2

theorem sdOut : ∀ (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState), SdOut v d o hasVert s
  | v, d, o, hasVert, s => fun hi hs hg hb => by
    unfold SidesOut; unfold GuardsOut at hg; unfold BookOut at hb
    refine wp_imp (wp_imp (wp_imp (wp_of_forall fun hv' s₁ ⟨hi₁, hs₁⟩ hg₁ hb₁ => ?_)
      (walkOutPre_inv hi hs hb.1)) hg) hb.2
    cases o with
    | back e cls dest =>
      exact finishEdge_sides (D := d) hg₁ hb₁ hi₁ hs₁
        (by rw [if_neg fun h => by obtain ⟨_, _, _, h⟩ := hb₁.tree.1 h; cases h])
    | tree e cls child =>
      try simp only [wp_modify] at hg₁ hb₁
      simp only [wp_modify]
      have pre : ∀ w outs, child = .node w outs →
          ({ { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne } with
            stackVerts := s₁.stackVerts.set! (d + 1) w } : WalkState).Inv' (d + 1) :=
        fun w outs _ => (hi₁.frame' (s' := { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne })).setSv w
      refine ⟨sdTree child (d + 1) _ pre hs₁.frame' hg₁.1 hb₁.1, ?_⟩
      exact wp_imp (wp_imp (wp_of_forall fun _ s₃ ⟨hg₃, hb₃⟩ ⟨hi₃, hs₃⟩ =>
        finishEdge_sides (D := d + 1) hg₃ hb₃ hi₃ hs₃ (by rw [if_pos (hb₃.tree.2 ⟨_, _, _, rfl⟩)]))
        (wp_and hg₁.2 hb₁.2)) (invTree child (d + 1) _ pre hs₁.frame' hg₁.1 hb₁.1)
end

/-- `FinishSides` holds at every `finishEdge` of the walk of `t`, given the ear guards and the
bookkeeping facts (the hypotheses of `walkTree_inv'`). -/
theorem walkTree_sides (t : DfsTree) (d : Nat) (s : WalkState)
    (hi : ∀ v outs, t = .node v outs → ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).Inv' d)
    (hs : Shape s) (hg : GuardsTree t d s) (hb : BookTree t d s) : SidesTree t d s :=
  sdTree t d s hi hs hg hb

end WalkState
end Spqr
