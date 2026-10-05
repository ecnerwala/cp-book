import Spqr.Proofs.RInvBase
import Spqr.Proofs.Blocks
import Spqr.EarFrontier
import Spqr.EarFrame
import Spqr.WalkPlace

/-!
# The child-return contract `RReturn` by the walk induction

`finishEdge_bot_keep`: a `finishEdge` whose boundary `origTstack` lies above the bottom `n₀`
(`hasVert = true → n₀ + 1 ≤ origTstack`) leaves the bottom `n₀` entries unchanged (`BotKeep`).
Loops 1–2 are read off the ear shape at `feS₂` (`EarFinish.close`: `c :: mid ++ [py, vy] ++ base`);
`closeVert'`, the P-check and the first-edge vertex push touch only the top one or two entries.
-/

namespace Spqr
open WalkM
namespace WalkState
variable {s : WalkState} {dfs : DfsData}

/-- The bottom `n` entries of the stack are unchanged. -/
def BotKeep (n : Nat) (s s' : WalkState) : Prop :=
  n ≤ s'.tstack.length ∧
    s'.tstack.drop (s'.tstack.length - n) = s.tstack.drop (s.tstack.length - n)

theorem BotKeep.refl {n : Nat} (hn : n ≤ s.tstack.length) : BotKeep n s s := ⟨hn, rfl⟩

theorem BotKeep.trans {n : Nat} {s₁ s₂ : WalkState} (h₁ : BotKeep n s s₁) (h₂ : BotKeep n s₁ s₂) :
    BotKeep n s s₂ := ⟨h₂.1, h₂.2.trans h₁.2⟩

theorem BotKeep.of_tstack {n : Nat} {s' : WalkState} (h : s'.tstack = s.tstack)
    (hn : n ≤ s.tstack.length) : BotKeep n s s' := ⟨by rw [h]; exact hn, by rw [h]⟩

theorem botKeep_cons {n : Nat} {s' : WalkState} {t : TEntry} (h : s'.tstack = t :: s.tstack)
    (hn : n ≤ s.tstack.length) : BotKeep n s s' := by
  refine ⟨by rw [h]; simp only [List.length_cons]; omega, ?_⟩
  rw [h, List.length_cons, Nat.succ_sub hn, List.drop_succ_cons]

theorem botKeep_head {n : Nat} {s' : WalkState} {a a' : TEntry} {rest : List TEntry}
    (hs : s.tstack = a :: rest) (hs' : s'.tstack = a' :: rest) (hn : n ≤ rest.length) :
    BotKeep n s s' := by
  refine ⟨by rw [hs']; simp only [List.length_cons]; omega, ?_⟩
  rw [hs, hs']
  simp only [List.length_cons]
  rw [show rest.length + 1 - n = (rest.length - n) + 1 by omega, List.drop_succ_cons, List.drop_succ_cons]

theorem botKeep_head2 {n : Nat} {s' : WalkState} {a b a' b' : TEntry} {rest : List TEntry}
    (hs : s.tstack = a :: b :: rest) (hs' : s'.tstack = a' :: b' :: rest) (hn : n ≤ rest.length) :
    BotKeep n s s' := by
  refine ⟨by rw [hs']; simp only [List.length_cons]; omega, ?_⟩
  rw [hs, hs']
  simp only [List.length_cons]
  rw [show rest.length + 1 + 1 - n = (rest.length - n) + 1 + 1 by omega]
  simp only [List.drop_succ_cons]

theorem botKeep_merge {n : Nat} {s' : WalkState} {a b m : TEntry} {rest : List TEntry}
    (hs : s.tstack = a :: b :: rest) (hs' : s'.tstack = m :: rest) (hn : n ≤ rest.length) :
    BotKeep n s s' := by
  refine ⟨by rw [hs']; simp only [List.length_cons]; omega, ?_⟩
  rw [hs, hs']
  simp only [List.length_cons]
  rw [show rest.length + 1 + 1 - n = (rest.length - n) + 1 + 1 by omega,
    show rest.length + 1 - n = (rest.length - n) + 1 by omega, List.drop_succ_cons,
    List.drop_succ_cons, List.drop_succ_cons]

theorem botKeep_mergeTop {n : Nat} (hlen : n + 2 ≤ s.tstack.length) :
    BotKeep n s (after mergeTstackTops s) := by
  match hts : s.tstack with
  | [] | [_] => rw [hts] at hlen; simp at hlen
  | a :: b :: rest =>
    have hts' : (after mergeTstackTops s).tstack = TEntry.mergeInto a b :: rest := by
      show (mergeTstackTops.run s).2.tstack = _
      rw [mergeTstackTops_run_eq s a b rest hts]
    exact botKeep_merge hts hts' (by rw [hts] at hlen; simp only [List.length_cons] at hlen; omega)

theorem botKeep_loop3 (origTstack : Nat) {n : Nat} (hn : n + 1 ≤ origTstack) : ∀ (k : Nat) (s : WalkState),
    (∀ j, j < k → result (loop3Cond origTstack) (iter mergeTstackTops j s) = true) →
    origTstack + 3 ≤ s.tstack.length → BotKeep n s (iter mergeTstackTops k s)
  | 0, s, _, h0 => BotKeep.refl (by omega)
  | k + 1, s, hj, hl => by
    have hc : decide (s.tstack.length > origTstack + 3) = true := hj 0 (Nat.succ_pos k)
    have hc' := of_decide_eq_true hc
    show BotKeep n s (iter mergeTstackTops k (mergeTstackTops.run s).2)
    refine (botKeep_mergeTop (by omega)).trans
      (botKeep_loop3 origTstack hn k (mergeTstackTops.run s).2 (fun j hjk => hj (j + 1) (Nat.succ_lt_succ hjk)) ?_)
    match hts : s.tstack with
    | [] | [_] => rw [hts] at hc'; simp at hc'
    | a :: b :: rest =>
      show origTstack + 3 ≤ (mergeTstackTops.run s).2.tstack.length
      rw [mergeTstackTops_run_eq s a b rest hts]; rw [hts] at hc'
      simp only [List.length_cons] at hc' ⊢; omega

/-- A fuelled loop with a pure condition is some iterate of its body, with the condition true on
all earlier iterates. -/
theorem loop_iter_exists (cond : WalkM Bool) (body : WalkM Unit) : ∀ (fuel : Nat) (s : WalkState),
    (∀ s, (cond.run s).2 = s) →
    ∃ k, after (loop fuel cond body) s = iter body k s ∧
      (∀ j, j < k → result cond (iter body j s) = true) ∧ True
  | 0, s, _ => ⟨0, rfl, fun j hj => absurd hj (Nat.not_lt_zero j), trivial⟩
  | fuel + 1, s, hcond => by
    unfold after
    rw [loop_succ_run fuel cond body s (hcond s)]
    by_cases hc : (cond.run s).1 = true
    · simp only [hc, ↓reduceIte]
      obtain ⟨k, heq, hj, -⟩ := loop_iter_exists cond body fuel (body.run s).2 hcond
      refine ⟨k + 1, heq, fun j hjk => ?_, trivial⟩
      cases j with
      | zero => exact hc
      | succ j => exact hj j (Nat.lt_of_succ_lt_succ hjk)
    · simp only [hc, Bool.false_eq_true, ↓reduceIte]
      exact ⟨0, rfl, fun j hj => absurd hj (Nat.not_lt_zero j), trivial⟩

theorem botKeep_unwrap {n : Nat} (ty : NodeType) (hlen : n + 2 ≤ s.tstack.length) :
    BotKeep n s (after (maybeUnwrapNxt ty) s) := by
  match hts : s.tstack with
  | [] | [_] => rw [hts] at hlen; simp at hlen
  | a :: b :: rest =>
    obtain ⟨b', hts', -, -⟩ := maybeUnwrapNxt_tstack ty a b rest hts
    exact botKeep_head2 hts hts' (by rw [hts] at hlen; simp only [List.length_cons] at hlen; omega)

/-- `closeVert'` keeps the bottom `n` when `n + 1 ≤ origTstack` and the stack has at least
`origTstack + 3` entries. -/
theorem botKeep_closeVert {n u : Nat} {edgeDir isType1 : Bool} {origTstack : Nat} {isSingle : Bool}
    (hn : n + 1 ≤ origTstack) (hlen : origTstack + 3 ≤ s.tstack.length) :
    BotKeep n s (after (WalkState.closeVert' u edgeDir isType1 origTstack isSingle) s) ∧
      origTstack + 1 ≤ (after (WalkState.closeVert' u edgeDir isType1 origTstack isSingle) s).tstack.length := by
  have K₁ : BotKeep n s (cvS₁ isType1 origTstack isSingle s) ∧
      origTstack + 3 ≤ (cvS₁ isType1 origTstack isSingle s).tstack.length := by
    cases isType1
    · show BotKeep n s ((vertPre false origTstack isSingle).run s).2 ∧
        origTstack + 3 ≤ ((vertPre false origTstack isSingle).run s).2.tstack.length
      simp only [WalkState.vertPre, Bool.not_false, ↓reduceIte, WalkM.run_bind, run_tstackSize, WalkM.pure_run]
      obtain ⟨k, heq, hj, -⟩ := loop_iter_exists (loop3Cond origTstack) mergeTstackTops s.tstack.length s
        (fun _ => rfl)
      have key : BotKeep n s (after (loop s.tstack.length (loop3Cond origTstack) mergeTstackTops) s) ∧
          origTstack + 3 ≤ (after (loop s.tstack.length (loop3Cond origTstack) mergeTstackTops) s).tstack.length := by
        rw [heq]; exact ⟨botKeep_loop3 origTstack hn k s hj hlen, loop3_iter_len origTstack k s hj hlen⟩
      exact key
    · exact ⟨BotKeep.refl (by omega), hlen⟩
  have K₂ : BotKeep n s (cvS₂ isType1 origTstack isSingle s) ∧
      origTstack + 3 ≤ (cvS₂ isType1 origTstack isSingle s).tstack.length := by
    cases isType1
    · exact K₁
    · have hl₁ : origTstack + 3 ≤ s.tstack.length := K₁.2
      show BotKeep n s ((some <$> maybeUnwrapNxt (if isSingle then NodeType.S else .R)).run s).2 ∧
        origTstack + 3 ≤ ((some <$> maybeUnwrapNxt (if isSingle then NodeType.S else .R)).run s).2.tstack.length
      rw [WalkM.map_run]
      refine ⟨botKeep_unwrap _ (by omega), ?_⟩
      match hts : s.tstack with
      | [] | [_] => rw [hts] at hl₁; simp at hl₁
      | a :: b :: rest =>
        obtain ⟨b', hts', -, -⟩ := maybeUnwrapNxt_tstack (if isSingle then NodeType.S else .R) a b rest hts
        have hts'' : ((maybeUnwrapNxt (if isSingle then NodeType.S else .R)).run s).2.tstack = a :: b' :: rest := hts'
        rw [hts'']; rw [hts] at hl₁; simpa using hl₁
  have hl₂ := K₂.2
  have K₃ : BotKeep n (cvS₂ isType1 origTstack isSingle s) (cvS₃ isType1 origTstack isSingle s) :=
    botKeep_mergeTop (by omega)
  have hl₃ : (cvS₃ isType1 origTstack isSingle s).tstack.length + 1 =
      (cvS₂ isType1 origTstack isSingle s).tstack.length := by
    match hts : (cvS₂ isType1 origTstack isSingle s).tstack with
    | [] | [_] => rw [hts] at hl₂; simp at hl₂
    | a :: b :: rest =>
      show (mergeTstackTops.run _).2.tstack.length + 1 = _
      rw [mergeTstackTops_run_eq _ a b rest hts]; rfl
  have K₄ : BotKeep n (cvS₃ isType1 origTstack isSingle s) (cvS₄ isType1 origTstack isSingle s) :=
    botKeep_mergeTop (by omega)
  have hl₄ : (cvS₄ isType1 origTstack isSingle s).tstack.length + 1 =
      (cvS₃ isType1 origTstack isSingle s).tstack.length := by
    match hts : (cvS₃ isType1 origTstack isSingle s).tstack with
    | [] | [_] => rw [hts] at hl₃; simp at hl₃; omega
    | a :: b :: rest =>
      show (mergeTstackTops.run _).2.tstack.length + 1 = _
      rw [mergeTstackTops_run_eq _ a b rest hts]; rfl
  match hts₄ : (cvS₄ isType1 origTstack isSingle s).tstack with
  | [] => rw [hts₄] at hl₄; simp at hl₄; omega
  | t :: rest =>
    have hl₄' : rest.length + 1 = (cvS₄ isType1 origTstack isSingle s).tstack.length := by rw [hts₄]; rfl
    have hts₅ : (cvS₅ u edgeDir isType1 origTstack isSingle s).tstack =
        { t with vStart := u, spans := setSides (!edgeDir) (t.spans.1 ++ t.spans.2) [] } :: rest := by
      show ((WalkState.retarget u edgeDir).run _).2.tstack = _
      rw [retarget_run_eq u edgeDir _ t rest hts₄]
    have K₅ : BotKeep n (cvS₄ isType1 origTstack isSingle s) (cvS₅ u edgeDir isType1 origTstack isSingle s) :=
      botKeep_head hts₄ hts₅ (by omega)
    have hl₅ : origTstack + 1 ≤ (cvS₅ u edgeDir isType1 origTstack isSingle s).tstack.length := by
      rw [hts₅]; simp only [List.length_cons]; omega
    have K := K₂.1.trans (K₃.trans (K₄.trans K₅))
    cases isType1
    · exact ⟨K, hl₅⟩
    · set item := ((maybeUnwrapNxt (if isSingle then NodeType.S else .R)).run s).1
      obtain ⟨x', it', hF, -, -, -⟩ := finishTstackTop_run item (cvS₅ u edgeDir true origTstack isSingle s) hts₅
      have hts₆ : ((finishTstackTop item).run (cvS₅ u edgeDir true origTstack isSingle s)).2.tstack = x' :: rest := by
        rw [hF]
      show BotKeep n s ((finishTstackTop item).run (cvS₅ u edgeDir true origTstack isSingle s)).2 ∧
        origTstack + 1 ≤ ((finishTstackTop item).run (cvS₅ u edgeDir true origTstack isSingle s)).2.tstack.length
      refine ⟨K.trans (botKeep_head hts₅ hts₆ (by omega)), ?_⟩
      rw [hts₆]; simp only [List.length_cons]; omega

theorem botKeep_finishP {n u lowval : Nat} {isType1 : Bool} (hlen : n + 2 ≤ s.tstack.length) :
    BotKeep n s (after (Spqr.finishP u lowval isType1) s) ∧
      s.tstack.length ≤ (after (Spqr.finishP u lowval isType1) s).tstack.length + 1 := by
  show BotKeep n s ((Spqr.finishP u lowval isType1).run s).2 ∧
    s.tstack.length ≤ ((Spqr.finishP u lowval isType1).run s).2.tstack.length + 1
  simp only [Spqr.finishP, WalkM.run_bind, run_condP]
  by_cases hc : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == u) &&
      (s.tstack.tail.head!.topDepth == lowval)) = true
  · simp only [hc, ↓reduceIte, WalkM.run_bind]
    match hts : s.tstack with
    | [] | [_] => rw [hts] at hlen; simp at hlen
    | a :: b :: rest =>
      obtain ⟨b', hts₁, -, -⟩ := maybeUnwrapNxt_tstack .P a b rest hts
      set s₁ := ((maybeUnwrapNxt .P).run s).2 with hs₁
      have hts₁' : s₁.tstack = a :: b' :: rest := hts₁
      have hn : n ≤ rest.length := by rw [hts] at hlen; simp only [List.length_cons] at hlen; omega
      have hts₂ : (mergeTstackTops.run s₁).2.tstack = TEntry.mergeInto a b' :: rest := by
        rw [mergeTstackTops_run_eq s₁ a b' rest hts₁']
      obtain ⟨x', it', hF, -, -, -⟩ := finishTstackTop_run ((maybeUnwrapNxt .P).run s).1 (mergeTstackTops.run s₁).2 hts₂
      have hts₃ : ((finishTstackTop ((maybeUnwrapNxt .P).run s).1).run (mergeTstackTops.run s₁).2).2.tstack =
          x' :: rest := by rw [hF]
      refine ⟨(botKeep_head2 hts hts₁' hn).trans ((botKeep_merge hts₁' hts₂ hn).trans
        (botKeep_head hts₂ hts₃ hn)), ?_⟩
      rw [hts₃]; simp
  · simp only [hc, Bool.false_eq_true, ↓reduceIte]
    exact ⟨BotKeep.refl (by omega), Nat.le_succ _⟩

theorem mergeTop_len (hlen : 2 ≤ s.tstack.length) :
    (after mergeTstackTops s).tstack.length + 1 = s.tstack.length := by
  match hts : s.tstack with
  | [] | [_] => rw [hts] at hlen; simp at hlen
  | a :: b :: rest =>
    have hts' : (after mergeTstackTops s).tstack = TEntry.mergeInto a b :: rest := by
      show (mergeTstackTops.run s).2.tstack = _
      rw [mergeTstackTops_run_eq s a b rest hts]
    rw [hts']; simp

theorem botKeep_finishTail {n u k : Nat} {hasVert isSingle : Bool}
    (hlen : hasVert = false → n + 1 ≤ s.tstack.length) (hlen' : n ≤ s.tstack.length) :
    BotKeep n s (after (Spqr.finishTail u k hasVert isSingle) s) ∧
      s.tstack.length ≤ (after (Spqr.finishTail u k hasVert isSingle) s).tstack.length := by
  cases hasVert
  · have K₁ : BotKeep n s (after (pushVertTstack u k) s) := botKeep_cons rfl hlen'
    have hl₁ : (after (pushVertTstack u k) s).tstack.length = s.tstack.length + 1 := rfl
    cases isSingle
    · have K₂ := botKeep_mergeTop (s := after (pushVertTstack u k) s) (n := n)
        (by have := hlen rfl; omega)
      have hl₂ := mergeTop_len (s := after (pushVertTstack u k) s) (by have := hlen rfl; omega)
      refine ⟨K₁.trans K₂, ?_⟩
      show s.tstack.length ≤ (after mergeTstackTops (after (pushVertTstack u k) s)).tstack.length
      omega
    · exact ⟨K₁, by show s.tstack.length ≤ s.tstack.length + 1; omega⟩
  · exact ⟨BotKeep.refl hlen', le_refl _⟩

theorem botKeep_finishRest {n u k lowval : Nat} {isType1 hasVert isSingle : Bool}
    (hlen : n + 2 ≤ s.tstack.length) :
    BotKeep n s (after (Spqr.finishRest u k lowval isType1 hasVert isSingle) s) ∧
      s.tstack.length ≤ (after (Spqr.finishRest u k lowval isType1 hasVert isSingle) s).tstack.length + 1 := by
  obtain ⟨K₁, hl₁⟩ := botKeep_finishP (u := u) (lowval := lowval) (isType1 := isType1) hlen
  obtain ⟨K₂, hl₂⟩ := botKeep_finishTail (s := after (Spqr.finishP u lowval isType1) s) (n := n)
    (u := u) (k := k) (hasVert := hasVert) (isSingle := isSingle) (fun _ => by omega) (by omega)
  refine ⟨K₁.trans K₂, ?_⟩
  show s.tstack.length ≤
    (after (Spqr.finishTail u k hasVert isSingle) (after (Spqr.finishP u lowval isType1) s)).tstack.length + 1
  omega

theorem botKeep_of_base {n : Nat} {s' : WalkState} {A A' base : List TEntry} (hs : s.tstack = A ++ base)
    (hs' : s'.tstack = A' ++ base) (hn : n ≤ base.length) : BotKeep n s s' := by
  have key : ∀ (B : List TEntry), (B ++ base).drop ((B ++ base).length - n) = base.drop (base.length - n) := by
    intro B
    rw [show (B ++ base).length - n = B.length + (base.length - n) by simp; omega, ← List.drop_drop,
      List.drop_left]
  refine ⟨by rw [hs']; simp; omega, ?_⟩
  rw [hs, hs', key, key]

/-- A `finishEdge` at `(curV, d)` whose boundary lies above the bottom `n₀` leaves those entries
unchanged. -/
theorem finishEdge_bot_keep (curV d lv : Nat) (kind : RetKind) (o : DfsOut) (origTstack n₀ : Nat)
    (hasVert : Bool) {sub base : List TEntry} (ho : o.cls = .ret lv kind) (hlow : lv < d)
    (hE : s.EarFinish curV d o hasVert sub base) (hl : base.length = origTstack)
    (hback : kind = .backEdge → hasVert = true)
    (hn : n₀ ≤ origTstack) (hnv : hasVert = true → n₀ + 1 ≤ origTstack) :
    BotKeep n₀ s (after (finishEdge curV d o origTstack hasVert) s) ∧
      origTstack + (if hasVert then 0 else 1) ≤
        (after (finishEdge curV d o origTstack hasVert) s).tstack.length := by
  have hlv : o.cls.lowval d = lv := by rw [ho]; rfl
  have hge : ¬ (lv ≥ d) := by omega
  have hlow' : o.cls.lowval d < d := by rw [hlv]; exact hlow
  have hsz : origTstack ≤ s.tstack.length := by rw [hE.tstack, List.length_append, hl]; omega
  have K₀ : BotKeep n₀ s (feS₀ d o s) := BotKeep.of_tstack rfl (by omega)
  by_cases hk : kind = .backEdge
  · subst hk
    have hhv := hback rfl
    subst hhv
    have ht' : o.cls.isTree = false := by rw [ho]; rfl
    have h1 : o.cls.isType1 = true := by rw [ho]; rfl
    have hl₂ : (feBack curV lv d o s).tstack.length = s.tstack.length + 1 := rfl
    have K₂ : BotKeep n₀ s (feBack curV lv d o s) := botKeep_cons (s := s) rfl (by omega)
    have hsz' : s.tstack.length = origTstack := by rw [hE.tstack, hE.back_nil ht']; simpa using hl
    show BotKeep n₀ s ((finishEdge curV d o origTstack true).run s).2 ∧
      origTstack + (if true then 0 else 1) ≤ ((finishEdge curV d o origTstack true).run s).2.tstack.length
    rw [finishEdge_eq]
    simp only [finishEdge', finishBack, hlv, hge, ↓reduceIte, WalkM.run_bind, WalkM.get_run,
      run_stackDir, run_makeVs, run_modifyItem, ht', Bool.false_eq_true, h1]
    obtain ⟨K₃, hl₃⟩ := botKeep_finishRest (s := feBack curV lv d o s) (n := n₀) (u := curV) (k := d)
      (lowval := lv) (isType1 := true) (hasVert := true) (isSingle := true) (by have := hnv rfl; omega)
    refine ⟨K₂.trans K₃, ?_⟩
    show origTstack + 0 ≤ (after (finishRest curV d lv true true true) (feBack curV lv d o s)).tstack.length
    omega
  · have ht : o.cls.isTree = true := by
      rw [ho]; cases kind with
      | backEdge => exact absurd rfl hk
      | type1Child => rfl
      | type2Child => rfl
    obtain ⟨c, mid, py, vy, hC⟩ := hE.close ht hlow'
    have K₂ : BotKeep n₀ s (feS₂ d o s) := botKeep_of_base hE.tstack hC.tstack (by omega)
    have hlen₂ : origTstack + 3 ≤ (feS₂ d o s).tstack.length := by
      rw [hC.tstack]; simp only [List.length_cons, List.length_append, hl]; omega
    show BotKeep n₀ s ((finishEdge curV d o origTstack hasVert).run s).2 ∧
      origTstack + (if hasVert then 0 else 1) ≤ ((finishEdge curV d o origTstack hasVert).run s).2.tstack.length
    rw [finishEdge_eq]
    simp only [finishEdge', finishTree, hlv, hge, ↓reduceIte, ht, closeVert_eq]
    cases hasVert
    · simp only [Bool.false_eq_true, ↓reduceIte]
      obtain ⟨K₃, hl₃⟩ := botKeep_finishRest (s := feS₂ d o s) (n := n₀) (u := curV) (k := d) (lowval := lv)
        (isType1 := o.cls.isType1) (hasVert := false) (isSingle := feSingle d o s) (by omega)
      refine ⟨K₂.trans K₃, ?_⟩
      show origTstack + 1 ≤
        (after (finishRest curV d lv o.cls.isType1 false (feSingle d o s)) (feS₂ d o s)).tstack.length
      omega
    · simp only [↓reduceIte]
      obtain ⟨K₃, hl₃⟩ := botKeep_closeVert (s := feS₂ d o s) (n := n₀) (u := curV)
        (edgeDir := s.stackDir[d]!) (isType1 := o.cls.isType1) (isSingle := feSingle d o s)
        (hnv rfl) hlen₂
      have hl₃' : origTstack + 1 ≤ (feS₃ curV d o origTstack s).tstack.length := hl₃
      have K₃' : BotKeep n₀ (feS₂ d o s) (feS₃ curV d o origTstack s) := K₃
      obtain ⟨K₄, hl₄⟩ := botKeep_finishRest (s := feS₃ curV d o origTstack s) (n := n₀) (u := curV) (k := d)
        (lowval := lv) (isType1 := o.cls.isType1) (hasVert := true) (isSingle := feB₃ curV d o origTstack s)
        (by have := hnv rfl; omega)
      refine ⟨K₂.trans (K₃'.trans K₄), ?_⟩
      show origTstack + 0 ≤ (after (finishRest curV d lv o.cls.isType1 true (feB₃ curV d o origTstack s))
        (feS₃ curV d o origTstack s)).tstack.length
      omega


theorem BotKeep.mono {n m : Nat} {s s' : WalkState} (h : BotKeep n s s') (hm : m ≤ n) :
    BotKeep m s s' := by
  obtain ⟨hl, hd⟩ := h
  have hlen := congrArg List.length hd
  simp only [List.length_drop] at hlen
  refine ⟨by omega, ?_⟩
  rw [show s'.tstack.length - m = (s'.tstack.length - n) + (n - m) by omega, ← List.drop_drop, hd,
    List.drop_drop]
  congr 1; omega

theorem earFinish_close_len {curV d lv : Nat} {kind : RetKind} {o : DfsOut} {origTstack : Nat}
    {hasVert : Bool} {sub base : List TEntry} (ho : o.cls = .ret lv kind) (hlow : lv < d)
    (hk : kind ≠ .backEdge) (hE : s.EarFinish curV d o hasVert sub base) (hl : base.length = origTstack) :
    origTstack + 3 ≤ (feS₂ d o s).tstack.length := by
  have hlow' : o.cls.lowval d < d := by rw [ho]; exact hlow
  have ht : o.cls.isTree = true := by
    rw [ho]; cases kind with
    | backEdge => exact absurd rfl hk
    | type1Child => rfl
    | type2Child => rfl
  obtain ⟨c, mid, py, vy, hC⟩ := hE.close ht hlow'
  rw [hC.tstack]; simp only [List.length_cons, List.length_append, hl]; omega

/-- The ancestor chain at `(v, d)`: `stackVerts[d] = v` and `stackVerts[k]` is the depth-`k`
ancestor of `v` for `k ≤ d`. -/
def AncChain (dfs : DfsData) (v d : Nat) (s : WalkState) : Prop :=
  s.stackVerts[d]! = v ∧ ∀ k, k ≤ d → dfs.Anc s.stackVerts[k]! v ∧ dfs.depth s.stackVerts[k]! = k

/-- An ancestor frame `(p, dp, n₀)`: the bottom `n₀` entries are settled at `(p, dp)`. -/
abbrev RFrame := Nat × Nat × Nat

/-- The state carried through the walk at `(v, d)`: every ancestor frame of `F` holds
positionally and the whole stack is settled at `(v, d)`. -/
structure RWalk (dfs : DfsData) (F : List RFrame) (v d : Nat) (s : WalkState) : Prop where
  frames : ∀ f ∈ F, f.2.2 ≤ s.tstack.length ∧ s.RInvG dfs f.1 f.2.1 f.2.2
  top : s.RInvTop dfs v d

mutual
/-- The side hypotheses of the child-return induction at the `walkTree` entries and `finishEdge`
sites of the walk of `t` from `s` (checked on seeds 0..300 in `checks/RFinishEdgeCheck.lean`):
at each entry the ancestor chain, the stability of `EntryR` under `stackVerts.set! d v`, and that
no entry starting at the parent tops out above `d - 1`; at each out-edge that it returns
(`lowval < d`), that no entry owns an edge of `v`'s vertex item before a vertex push, and the
R-maximality content `FinishRShape` at tree-edge sites. -/
def RSideTree (dfs : DfsData) (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  match t with
  | .node v outs =>
    AncChain dfs v d { s with stackVerts := s.stackVerts.set! d v } ∧
    (∀ t ∈ s.tstack, s.EntryR dfs t →
      ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).EntryR dfs t) ∧
    (∀ t ∈ s.tstack, t.vStart = s.stackVerts[d - 1]! → t.topDepth ≤ d - 1) ∧
    RSideOuts dfs v d outs false { s with stackVerts := s.stackVerts.set! d v }

def RSideOuts (dfs : DfsData) (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) :
    Prop :=
  match outs with
  | [] => hasVert = false → VertFree v s
  | o :: rest => RSideOut dfs v d o hasVert s ∧
      wp (walkOut v d o hasVert) (fun hasVert' s' => RSideOuts dfs v d rest hasVert' s') s

def RSideOut (dfs : DfsData) (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  o.cls.lowval d < d ∧ (hasVert = false → VertFree v s) ∧
  wp (walkOutPre v d o hasVert) (fun hasVert' s₁ =>
    match o with
    | .tree _ _ child =>
      wp (modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }) (fun _ s₂ =>
        RSideTree dfs child (d + 1) s₂ ∧
        wp (walkTree child (d + 1))
          (fun _ s₃ => FinishRShape dfs v d o s₁.tstack.length hasVert' s₃) s₂) s₁
    | .back .. => True) s
end

theorem RWalk.of_eq {F : List RFrame} {v d : Nat} {s' : WalkState} (hg : s'.g = s.g)
    (hi : s'.items = s.items) (hsv : s'.stackVerts = s.stackVerts) (hts : s'.tstack = s.tstack)
    (h : RWalk dfs F v d s) : RWalk dfs F v d s' :=
  ⟨fun f hf => ⟨by rw [hts]; exact (h.frames f hf).1, RInvG.of_eq' hg hi hsv hts (h.frames f hf).2⟩,
    RInvTop.of_eq hg hi hsv hts h.top⟩

theorem RWalk.pushVert {F : List RFrame} {v d k : Nat} (hown : VertFree v s)
    (h : RWalk dfs F v d s) : RWalk dfs F v d (after (pushVertTstack v k) s) := by
  refine ⟨fun f hf => ⟨?_, RInvG.pushVert_pos v k (h.frames f hf).1 hown (h.frames f hf).2⟩, ?_⟩
  · show f.2.2 ≤ (_ :: s.tstack).length
    have := (h.frames f hf).1; simp only [List.length_cons]; omega
  · exact (h.top.pushVert_own hown).toTop fun _ _ hne => (hne rfl).elim

theorem walkOutPre_r {v d : Nat} {o : DfsOut} {hasVert : Bool} {F : List RFrame} {B : Nat}
    (hfree : hasVert = false → VertFree v s) (hW : RWalk dfs F v d s) (hBl : B ≤ s.tstack.length)
    (hnv : hasVert = true → B + 1 ≤ s.tstack.length) :
    wp (walkOutPre v d o hasVert) (fun hv' s₁ => RWalk dfs F v d s₁ ∧ BotKeep B s s₁ ∧
      (hv' = true → B + 1 ≤ s₁.tstack.length) ∧ (hasVert = true → hv' = true) ∧
      (o.cls.lowval d < d → o.cls.isType1 = true → hv' = true) ∧
      s₁.g = s.g ∧ s₁.stackVerts = s.stackVerts ∧ s₁.items = s.items) s := by
  unfold walkOutPre
  simp only [wp_bind, wp_stackDir, wp_setStackDir, wp_ite]
  have hW' : RWalk dfs F v d ({ s with stackDir := s.stackDir.set! d (if o.cls.lowval d ≥ d then false else !s.stackDir[o.cls.lowval d]!) } : WalkState) :=
    RWalk.of_eq (s := s) rfl rfl rfl rfl hW
  split
  · rename_i hc
    have hf : hasVert = false := by cases hasVert <;> simp_all
    subst hf
    show RWalk dfs F v d _ ∧ BotKeep B s _ ∧ (true = true → B + 1 ≤ _) ∧ (false = true → true = true) ∧
      (o.cls.lowval d < d → o.cls.isType1 = true → true = true) ∧ _ ∧ _ ∧ _
    refine ⟨hW'.pushVert (hfree rfl), botKeep_cons (s := s) rfl hBl, fun _ => ?_, fun h => absurd h Bool.false_ne_true,
      fun _ _ => rfl, rfl, rfl, rfl⟩
    show B + 1 ≤ (_ :: s.tstack).length
    simp only [List.length_cons]; omega
  · rename_i hc
    show RWalk dfs F v d _ ∧ BotKeep B s _ ∧ (hasVert = true → B + 1 ≤ _) ∧ (hasVert = true → hasVert = true) ∧
      (o.cls.lowval d < d → o.cls.isType1 = true → hasVert = true) ∧ _ ∧ _ ∧ _
    refine ⟨hW', BotKeep.refl hBl, hnv, fun h => h, fun h1 h2 => ?_, rfl, rfl, rfl⟩
    cases hasVert
    · exact absurd (by simp [h1, h2]) hc
    · rfl

/-- The child-return induction, `walkTree` level (`d = dp + 1`, parent `stackVerts[dp]` settled on
the whole stack at entry): the ancestor frames `F` hold positionally at exit, the walked vertex's
own `RInvTop` holds, the entry stack is kept as the bottom, and the graph and `stackVerts[< d]`
are unchanged. -/
abbrev RRTree (dfs : DfsData) (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  ∀ (F : List RFrame) (dp : Nat), d = dp + 1 →
    (∀ v outs, t = .node v outs → ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).Inv' d) →
    Shape s → GuardsTree t d s → BookTree t d s → RSideTree dfs t d s →
    s.g.TwoConnected → dfs.Spec s.g → dfs.Rooted s.g →
    (∀ f ∈ F, f.2.2 ≤ s.tstack.length ∧ s.RInvG dfs f.1 f.2.1 f.2.2) →
    s.RInvTop dfs s.stackVerts[dp]! dp →
    wp (walkTree t d) (fun _ s' => (∀ v outs, t = .node v outs → RWalk dfs F v d s') ∧
      BotKeep s.tstack.length s s' ∧ s'.g = s.g ∧
      ∀ k, k < d → s'.stackVerts[k]! = s.stackVerts[k]!) s

abbrev RROuts (dfs : DfsData) (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) :
    Prop :=
  ∀ (F : List RFrame) (B : Nat),
    s.Inv' d → Shape s → GuardsOuts v d outs hasVert s → BookOuts v d outs hasVert s →
    RSideOuts dfs v d outs hasVert s → s.g.TwoConnected → dfs.Spec s.g → dfs.Rooted s.g →
    AncChain dfs v d s → RWalk dfs F v d s → (∀ f ∈ F, f.2.2 ≤ B) → B ≤ s.tstack.length →
    (hasVert = true → B + 1 ≤ s.tstack.length) →
    wp (walkOuts v d outs hasVert) (fun hasVert' s' => RWalk dfs F v d s' ∧ BotKeep B s s' ∧
      (hasVert' = true → B + 1 ≤ s'.tstack.length) ∧ (hasVert' = false → VertFree v s') ∧
      s'.g = s.g ∧ ∀ k, k ≤ d → s'.stackVerts[k]! = s.stackVerts[k]!) s

abbrev RROut (dfs : DfsData) (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  ∀ (F : List RFrame) (B : Nat),
    s.Inv' d → Shape s → GuardsOut v d o hasVert s → BookOut v d o hasVert s →
    RSideOut dfs v d o hasVert s → s.g.TwoConnected → dfs.Spec s.g → dfs.Rooted s.g →
    AncChain dfs v d s → RWalk dfs F v d s → (∀ f ∈ F, f.2.2 ≤ B) → B ≤ s.tstack.length →
    (hasVert = true → B + 1 ≤ s.tstack.length) →
    wp (walkOut v d o hasVert) (fun hasVert' s' => RWalk dfs F v d s' ∧ BotKeep B s s' ∧
      (hasVert' = true → B + 1 ≤ s'.tstack.length) ∧
      s'.g = s.g ∧ ∀ k, k ≤ d → s'.stackVerts[k]! = s.stackVerts[k]!) s

mutual
theorem rrTree : ∀ (t : DfsTree) (d : Nat) (s : WalkState), RRTree dfs t d s
  | .node v outs, d, s => fun F dp hdp hi hs hg hb hr h2 hsp hrt hF hpar => by
    subst hdp
    unfold walkTree
    simp only [wp_bind, wp_modify]
    unfold GuardsTree at hg; unfold BookTree at hb; unfold RSideTree at hr
    obtain ⟨hanc, hstab, hvs, hr⟩ := hr
    simp only [Nat.add_sub_cancel] at hvs
    have hW₀ : RWalk dfs F v (dp + 1) { s with stackVerts := s.stackVerts.set! (dp + 1) v } :=
      ⟨fun f hf => ⟨(hF f hf).1, ⟨fun t ht hd hne => hstab t (List.mem_of_mem_drop ht) ((hF f hf).2.entries t ht hd hne),
          (hF f hf).2.disj⟩⟩,
        ⟨fun t ht hd hne => hstab t ht (hpar.entries t ht (by omega) fun h => by
          have := hvs t ht h; omega), hpar.disj⟩⟩
    refine wp_imp (wp_imp (wp_of_forall fun hv s₁ ⟨hi₁, hs₁, hvb⟩
        ⟨hW₁, hK₁, hn₁, hfree₁, hg₁, hsv₁⟩ => ?_)
      (invOuts v (dp + 1) outs false _ (hi v outs rfl) hs.frame' hg hb))
      (rrOuts v (dp + 1) outs false _ F s.tstack.length (hi v outs rfl) hs.frame' hg hb hr h2 hsp hrt
        hanc hW₀ (fun f hf => (hF f hf).1) (le_refl _) (fun h => nomatch h))
    have hsv : ∀ k, k < dp + 1 → s₁.stackVerts[k]! = s.stackVerts[k]! := fun k hk => by
      rw [hsv₁ k (by omega)]
      exact getElem!_set!_ne _ _ _ _ (by omega)
    cases hv
    · simp only [Bool.false_eq_true, ↓reduceIte, wp_bind, wp_setStackDir]
      have hW₂ : RWalk dfs F v (dp + 1) ({ s₁ with stackDir := s₁.stackDir.set! (dp + 1) true } : WalkState) :=
        RWalk.of_eq (s := s₁) rfl rfl rfl rfl hW₁
      refine ⟨fun v' outs' h => ?_, ?_, hg₁, hsv⟩
      · cases h; exact hW₂.pushVert (hfree₁ rfl)
      · exact ((BotKeep.of_tstack (s := s) rfl (le_refl _)).trans hK₁).trans
          (botKeep_cons (s := { s₁ with stackDir := s₁.stackDir.set! (dp + 1) true }) rfl hK₁.1)
    · exact ⟨fun v' outs' h => by cases h; exact hW₁,
        (BotKeep.of_tstack (s := s) rfl (le_refl _)).trans hK₁, hg₁, hsv⟩

theorem rrOuts : ∀ (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState),
    RROuts dfs v d outs hasVert s
  | v, d, [], hasVert, s => fun F B hi hs hg hb hr h2 hsp hrt hanc hW hB hBl hnv => by
    unfold RSideOuts at hr
    unfold walkOuts
    rw [wp_pure]
    exact ⟨hW, BotKeep.refl hBl, hnv, hr, rfl, fun _ _ => rfl⟩
  | v, d, o :: rest, hasVert, s => fun F B hi hs hg hb hr h2 hsp hrt hanc hW hB hBl hnv => by
    unfold GuardsOuts at hg; unfold BookOuts at hb; unfold RSideOuts at hr
    unfold walkOuts
    simp only [wp_bind]
    refine wp_imp (wp_imp (wp_imp (wp_imp (wp_imp (wp_of_forall fun hv' s' ⟨hi', hs'⟩ hg' hb' hr'
        ⟨hW', hK', hn', hg'', hsv'⟩ => ?_)
      (invOut v d o hasVert s hi hs hg.1 hb.1)) hg.2) hb.2) hr.2)
      (rrOut v d o hasVert s F B hi hs hg.1 hb.1 hr.1 h2 hsp hrt hanc hW hB hBl hnv)
    have hanc' : AncChain dfs v d s' :=
      ⟨(hsv' d (le_refl _)).trans hanc.1, fun k hk => by rw [hsv' k hk]; exact hanc.2 k hk⟩
    refine wp_imp (wp_of_forall fun hv'' s'' ⟨R1, R2, R3, R4, R5, R6⟩ =>
      ⟨R1, hK'.trans R2, R3, R4, R5.trans hg'', fun k hk => (R6 k hk).trans (hsv' k hk)⟩)
      (rrOuts v d rest hv' s' F B hi' hs' hg' hb' hr' (by rw [hg'']; exact h2) (by rw [hg'']; exact hsp)
        (by rw [hg'']; exact hrt) hanc' hW' hB hK'.1 hn')

theorem rrOut : ∀ (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState), RROut dfs v d o hasVert s
  | v, d, o, hasVert, s => fun F B hi hs hg hb hr h2 hsp hrt hanc hW hB hBl hnv => by
    rw [walkOut_eq, wp_bind]
    unfold GuardsOut at hg; unfold BookOut at hb; unfold RSideOut at hr
    obtain ⟨hlow, hfree, hr⟩ := hr
    refine wp_imp (wp_imp (wp_imp (wp_imp (wp_imp (wp_of_forall fun hv₁ s₁ ⟨hi₁, hs₁⟩ hg₁ hb₁ hr₁
        ⟨hW₁, hK₁, hn₁, hhv₁, hpush₁, hg₁', hsv₁, hit₁⟩ => ?_)
      (walkOutPre_inv hi hs hb.1)) hg) hb.2) hr) (walkOutPre_r hfree hW hBl hnv)
    have hanc₁ : AncChain dfs v d s₁ := by rw [AncChain, hsv₁]; exact hanc
    have h2₁ : s₁.g.TwoConnected := by rw [hg₁']; exact h2
    have hsp₁ : dfs.Spec s₁.g := by rw [hg₁']; exact hsp
    have hrt₁ : dfs.Rooted s₁.g := by rw [hg₁']; exact hrt
    obtain ⟨lv, kind, ho, hl⟩ := ret_of_lowval_lt hlow
    unfold walkOutRest
    rw [wp_bind, wp_tstackSize]
    cases o with
    | back e cls dest =>
      try simp only at hg₁ hb₁
      simp only
      have hk : kind = .backEdge := by
        cases kind with
        | backEdge => rfl
        | type1Child => exact absurd (hb₁.tree.1 (by rw [ho]; rfl)) (by simp)
        | type2Child => exact absurd (hb₁.tree.1 (by rw [ho]; rfl)) (by simp)
      subst hk
      have hhv : hv₁ = true := hpush₁ hlow (by rw [ho]; rfl)
      subst hhv
      have hD : d = if (DfsOut.back e cls dest).cls.isTree then d + 1 else d := by
        rw [if_neg fun h => by obtain ⟨_, _, _, h⟩ := hb₁.tree.1 h; cases h]
      obtain ⟨sub, base, hlen, hE⟩ := hb₁.ear
      have hok := finishOk_of_guards ho hl hg₁ hE hlen hi₁ hs₁ hD hb₁.v_lt hb₁.e_lt hb₁.q
        (hb₁.ends lv .backEdge ho) hb₁.vert
      have hfront := finishEdge_frontier hb₁ hi₁ hs₁ hD
      have hhd : dfs.depth v = d := by
        have := (hanc₁.2 d (le_refl _)).2; rwa [hanc₁.1] at this
      have hR : s₁.RInvFront dfs v d s₁.tstack.length := hW₁.top.toFront _
      have hbot := finishEdge_bot_keep v d lv .backEdge _ s₁.tstack.length B true ho hl hE hlen
        (fun _ => rfl) hK₁.1 hn₁
      have hstep := finishEdge_inv v d lv .backEdge _ s₁.tstack.length true ho hl hb₁.v_lt hi₁ hs₁ hok
      refine ⟨⟨fun f hf => ⟨?_, ?_⟩, ?_⟩, hK₁.trans hbot.1, fun _ => ?_, hstep.g.trans hg₁',
        fun k hk => by rw [hstep.sv]; exact congrArg (·[k]!) hsv₁⟩
      · show f.2.2 ≤ (after (finishEdge v d _ s₁.tstack.length true) s₁).tstack.length
        have := hbot.1.1; have := hB f hf; omega
      · exact finishEdge_rInvG_base f.1 f.2.1 v d lv .backEdge _ s₁.tstack.length f.2.2 true ho hl
          hb₁.v_lt hi₁ hs₁ hok hfront (fun h => absurd rfl h) (fun _ => ⟨rfl, le_refl _⟩)
          (fun h => absurd rfl h) (by have := hB f hf; have := hK₁.1; omega)
          (fun _ => by have := hB f hf; have := hn₁ rfl; omega)
          (hW₁.frames f hf).2
      · exact finishEdge_rInvTop v d lv .backEdge _ s₁.tstack.length true ho hl hb₁.v_lt hi₁ hs₁ hok
          hg₁ hfront (fun h => absurd rfl h) (fun _ => ⟨rfl, le_refl _⟩) h2₁ hsp₁ hrt₁ hhd hanc₁.1
          hanc₁.2 hR
      · show B + 1 ≤ (after (finishEdge v d _ s₁.tstack.length true) s₁).tstack.length
        have h1 := hbot.2; have h2 := hn₁ rfl; simp only [↓reduceIte] at h1; omega
    | tree e cls child =>
      obtain ⟨c, couts⟩ := child
      try simp only [wp_modify] at hg₁ hb₁ hr₁
      simp only [wp_bind, wp_modify]
      have hW₂ : RWalk dfs F v d ({ s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne } : WalkState) :=
        RWalk.of_eq (s := s₁) rfl rfl rfl rfl hW₁
      have hpar : ({ s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne } : WalkState).RInvTop dfs
          ({ s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne } : WalkState).stackVerts[d]! d :=
        RInvTop.of_eq (s := s₁) rfl rfl rfl rfl (by
          show s₁.RInvTop dfs s₁.stackVerts[d]! d
          rw [hanc₁.1]; exact hW₁.top)
      refine wp_imp (wp_imp (wp_imp (wp_imp (wp_of_forall fun _ s₃ ⟨hg₃, hb₃⟩ hr₃ ⟨hi₃, hs₃⟩
          ⟨hW₃, hK₃, hg₃', hsv₃⟩ => ?body) (wp_and hg₁.2 hb₁.2)) hr₁.2)
        (invTree (.node c couts) (d + 1) _ ?pre hs₁.frame' hg₁.1 hb₁.1))
        (rrTree (.node c couts) (d + 1) _ ((v, d, s₁.tstack.length) :: F) d rfl ?pre hs₁.frame' hg₁.1
          hb₁.1 hr₁.1 h2₁ hsp₁ hrt₁ ?frames hpar)
      case pre => exact fun w outs _ => (hi₁.frame' (s' := { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne })).setSv w
      case frames =>
        intro f hf
        simp only [List.mem_cons] at hf
        rcases hf with rfl | hf
        · exact ⟨le_refl _, hW₂.top.toG _⟩
        · exact hW₂.frames f hf
      case body =>
        have hW₃ := hW₃ c couts rfl
        have hanc₃ : AncChain dfs v d s₃ :=
          ⟨(hsv₃ d (by omega)).trans hanc₁.1, fun k hk => by rw [hsv₃ k (by omega)]; exact hanc₁.2 k hk⟩
        have h2₃ : s₃.g.TwoConnected := by rw [hg₃']; exact h2₁
        have hsp₃ : dfs.Spec s₃.g := by rw [hg₃']; exact hsp₁
        have hrt₃ : dfs.Rooted s₃.g := by rw [hg₃']; exact hrt₁
        have hk : kind ≠ .backEdge := fun hk => by
          subst hk
          have := hb₃.tree.2 ⟨_, _, _, rfl⟩
          rw [ho] at this; cases this
        have hD : d + 1 = if (DfsOut.tree e cls (.node c couts)).cls.isTree then d + 1 else d := by
          rw [if_pos (hb₃.tree.2 ⟨_, _, _, rfl⟩)]
        obtain ⟨sub, base, hlen, hE⟩ := hb₃.ear
        have hok := finishOk_of_guards ho hl hg₃ hE hlen hi₃ hs₃ hD hb₃.v_lt hb₃.e_lt hb₃.q
          (hb₃.ends lv kind ho) hb₃.vert
        have hfront := finishEdge_frontier hb₃ hi₃ hs₃ hD
        have hhd : dfs.depth v = d := by
          have := (hanc₃.2 d (le_refl _)).2; rwa [hanc₃.1] at this
        have hclose := earFinish_close_len ho hl hk hE hlen
        have hfr := (hW₃.frames (v, d, s₁.tstack.length) (List.mem_cons_self ..)).2
        have hR : s₃.RInvFront dfs v d s₁.tstack.length := ⟨hfr.entries, hfr.disj⟩
        have hbot := finishEdge_bot_keep v d lv kind _ s₁.tstack.length B hv₁ ho hl hE hlen
          (fun h => absurd h hk) hK₁.1 hn₁
        have hstep := finishEdge_inv v d lv kind _ s₁.tstack.length hv₁ ho hl hb₃.v_lt hi₃ hs₃ hok
        refine ⟨⟨fun f hf => ⟨?_, ?_⟩, ?_⟩, (hK₁.trans (hK₃.mono hK₁.1)).trans hbot.1, fun hv => ?_,
          hstep.g.trans (hg₃'.trans hg₁'),
          fun k hk => by rw [hstep.sv, hsv₃ k (by omega)]; exact congrArg (·[k]!) hsv₁⟩
        · show f.2.2 ≤ (after (finishEdge v d _ s₁.tstack.length hv₁) s₃).tstack.length
          have := hbot.1.1; have := hB f hf; omega
        · exact finishEdge_rInvG_base f.1 f.2.1 v d lv kind _ s₁.tstack.length f.2.2 hv₁ ho hl
            hb₃.v_lt hi₃ hs₃ hok hfront (fun _ => hr₃) (fun h => absurd h hk) (fun _ => hclose)
            (by have := hB f hf; have := hK₁.1; omega) (fun h => by have := hB f hf; have := hn₁ h; omega)
            (hW₃.frames f (List.mem_cons_of_mem _ hf)).2
        · exact finishEdge_rInvTop v d lv kind _ s₁.tstack.length hv₁ ho hl hb₃.v_lt hi₃ hs₃ hok
            hg₃ hfront (fun _ => hr₃) (fun h => absurd h hk) h2₃ hsp₃ hrt₃ hhd hanc₃.1 hanc₃.2 hR
        · show B + 1 ≤ (after (finishEdge v d _ s₁.tstack.length hv₁) s₃).tstack.length
          have := hbot.2
          cases hv₁
          · simp only [Bool.false_eq_true, ↓reduceIte] at this; have := hK₁.1; omega
          · simp only [↓reduceIte] at this; have := hn₁ rfl; omega
end

/-- The provisional child-return contract from the induction, with the side facts (`RSideTree`)
and the ear bookkeeping (`BookTree`) as explicit hypotheses. -/
theorem walkTree_rReturn (s : WalkState) (d c : Nat) (outs : List DfsOut)
    (hi : s.Inv' d) (hs : Shape s) (hg : GuardsTree (.node c outs) (d + 1) s)
    (hb : BookTree (.node c outs) (d + 1) s) (hr : RSideTree dfs (.node c outs) (d + 1) s)
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hR : s.RInvTop dfs s.stackVerts[d]! d) :
    RReturn s (after (walkTree (.node c outs) (d + 1)) s) dfs s.stackVerts[d]! d := by
  have h := rrTree (.node c outs) (d + 1) s [(s.stackVerts[d]!, d, s.tstack.length)] d rfl
    (fun _ _ h => by cases h; exact hi.setSv _) hs hg hb hr h2 hsp hrt
    (fun f hf => by
      simp only [List.mem_singleton] at hf
      subst hf; exact ⟨le_refl _, hR.toG _⟩) hR
  obtain ⟨hW, hK, -, -⟩ := h
  have hW := hW c outs rfl
  have hfr := (hW.frames _ (List.mem_singleton_self _)).2
  refine ⟨hK.1, ?_, hfr.entries, hW.top.disj⟩
  have := hK.2; rwa [Nat.sub_self, List.drop_zero] at this

/-- The DFS-layer part of `RSideOut` (`ret`): every out-edge of a non-root vertex `c` returns,
`o.cls.lowval (depth c) < depth c`. Tree edges: `bridge`/`component` children contradict
`child_returns_above`; `ret` children by `ret_lowpt`. Back edges: a self-loop at `c` is not
`EdgeConn (· ≠ c)` to the parent edge; a `ret _ .backEdge` lands at a proper ancestor. -/
theorem _root_.Spqr.DfsData.Spec.outs_lowval_lt {g : Graph} {dfs : DfsData} (hs : dfs.Spec g)
    (h2 : g.TwoConnected) {p c : Nat} (hp : dfs.IsParent p c) {o : DfsOut} (ho : o ∈ dfs.outs c) :
    o.cls.lowval (dfs.depth c) < dfs.depth c := by
  have hdc := hs.depth_parent _ _ hp
  cases ht : o.isTree with
  | true =>
    rcases DfsData.tree_cls hs ho ht with h | h | ⟨l, k, h⟩
    · obtain ⟨l, hl, hr⟩ := DfsData.child_returns_above hs h2 hp (DfsData.isParent_of_tree ho ht)
      exact absurd hr (((hs.cls_bridge c o ho).mp h).2 l (Nat.le_of_lt hl))
    · obtain ⟨l, hl, hr⟩ := DfsData.child_returns_above hs h2 hp (DfsData.isParent_of_tree ho ht)
      exact absurd hr (((hs.cls_component c o ho).mp h).2.2 l hl)
    · rw [h]; exact (DfsData.ret_lowpt hs ho ht h).1
  | false =>
    cases hc : o.cls with
    | bridge => have := ((hs.cls_bridge c o ho).mp hc).1; rw [ht] at this; cases this
    | component => have := ((hs.cls_component c o ho).mp hc).1; rw [ht] at this; cases this
    | selfLoop =>
      exfalso
      have hdest := ((hs.cls_selfLoop c o ho).mp hc).2
      have hj := hs.joins c o ho
      rw [hdest] at hj
      obtain ⟨e', he'⟩ := hp.joins hs
      rcases h2 c o.e e' hj.lt he'.lt with h | ⟨x, y, hx, hy, hr⟩
      · subst h
        rcases hj.eq_or he' with ⟨h1, -⟩ | ⟨-, h1⟩ <;> (subst h1; omega)
      · obtain ⟨z, hz⟩ := hx
        rcases hj.eq_or hz with ⟨h1, -⟩ | ⟨-, h1⟩ <;> exact hr.ok_left h1.symm
    | ret l k =>
      cases k with
      | type1Child => have := ((hs.cls_type1 c o ho l).mp hc).1; rw [ht] at this; cases this
      | type2Child => have := ((hs.cls_type2 c o ho l).mp hc).1; rw [ht] at this; cases this
      | backEdge =>
        obtain ⟨-, hne, hl⟩ := (hs.cls_backEdge c o ho l).mp hc
        have hanc := hs.back_anc c o ho ht
        show l < dfs.depth c
        rw [← hl]; exact hanc.depth_lt hs hne

/-- The ancestor chain at a child entry `(c, d + 1)` from the parent's chain (`RSideTree`'s first
conjunct; `d + 1 < stackVerts.size` is a layout fact, `EarWalk`'s `sv : stackVerts.size = g.nv`). -/
theorem ancChain_child {d c : Nat} (hsp : dfs.Spec s.g) (hsize : d + 1 < s.stackVerts.size)
    (hp : dfs.IsParent s.stackVerts[d]! c)
    (hchain : ∀ k, k ≤ d → dfs.Anc s.stackVerts[k]! s.stackVerts[d]! ∧ dfs.depth s.stackVerts[k]! = k) :
    AncChain dfs c (d + 1) { s with stackVerts := s.stackVerts.set! (d + 1) c } := by
  have hself : (s.stackVerts.set! (d + 1) c)[d + 1]! = c := getElem!_set!_self _ _ _ hsize
  refine ⟨hself, fun k hk => ?_⟩
  rcases Nat.lt_or_eq_of_le hk with hk | rfl
  · show dfs.Anc (s.stackVerts.set! (d + 1) c)[k]! c ∧ dfs.depth (s.stackVerts.set! (d + 1) c)[k]! = k
    rw [getElem!_set!_ne _ _ _ _ (by omega)]
    exact ⟨(hchain k (by omega)).1.trans hp.anc, (hchain k (by omega)).2⟩
  · show dfs.Anc (s.stackVerts.set! (d + 1) c)[d + 1]! c ∧ dfs.depth (s.stackVerts.set! (d + 1) c)[d + 1]! = d + 1
    rw [hself]
    exact ⟨.refl _, by rw [hsp.depth_parent _ _ hp, (hchain d (le_refl _)).2]⟩

/-- **Named admission** — the standalone form `WalkTreeRReturnSpec` needs, kept only because that
statement is fixed; nothing else consumes it (the root glue uses `walkTree_rSide`). It is
under-hypothesized rather than false-by-counterexample: `BookTree` at a non-root entry is ear
content that only the root-threaded `bTree` produces, and `AncChain` for `k < d` constrains
`stackVerts` below `d`, which none of the hypotheses mention. -/
theorem walkTree_rSide_spec (s : WalkState) (d c : Nat) (outs : List DfsOut)
    (hi : s.Inv' d) (hs : Shape s) (hg : GuardsTree (.node c outs) (d + 1) s)
    (hf : FrontiersTree (.node c outs) (d + 1) s) (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g)
    (hrt : dfs.Rooted s.g) (hp : dfs.IsParent s.stackVerts[d]! c) (ho : outs = dfs.outs c)
    (hR : s.RInvTop dfs s.stackVerts[d]! d) :
    BookTree (.node c outs) (d + 1) s ∧ RSideTree dfs (.node c outs) (d + 1) s := by
  sorry

theorem walkTreeRReturnSpec : WalkTreeRReturnSpec dfs :=
  fun s d c outs hi hs hg hf h2 hsp hrt hp ho hR =>
    have ⟨hb, hr⟩ := walkTree_rSide_spec s d c outs hi hs hg hf h2 hsp hrt hp ho hR
    walkTree_rReturn s d c outs hi hs hg hb hr h2 hsp hrt hR


end WalkState
end Spqr
