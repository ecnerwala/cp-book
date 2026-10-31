import Spqr.Proofs.RInvTree

/-!
# The parent's base across a `finishEdge` of the child (`finishEdge_rInvG_base`)

The frame statement of the `WalkTreeRReturnSpec` induction design (PROOF.md §4.5): a `finishEdge`
at `(curV, d)` whose schedule boundary `origTstack` lies above the bottom `n₀` entries preserves
`RInvG dfs p dp n₀` — the bottom `n₀` entries are untouched, so their `EntryR` at `(p, dp)` is
kept, and the whole stack stays edge-disjoint. The content is purely positional: every
`RInvG.*` primitive takes its exemption hypothesis only when the touched entry lies in the
bottom `n₀`, so what is needed per primitive is a stack-length bound. These come from
`Frontier.loop1/2/3` (the loop iterates), `hclose` (`origTstack + 3 ≤ length` after loop 2 — the
child's entries `c :: mid ++ [py, vy]` of `EarFinish.loops`; `FinishGuards` gives this only for
type-2 edges with `hasVert`), and `hasVert = true → n₀ + 1 ≤ origTstack` (the vertex entry of
`curV` was pushed after the descent into `p`'s child), with `n₀ ≤ origTstack`.
-/

namespace Spqr
open WalkM
namespace WalkState
variable {s : WalkState} {dfs : DfsData}

theorem RInvG.of_eq' {s' : WalkState} {v d n : Nat} (hg : s'.g = s.g) (hi : s'.items = s.items)
    (hsv : s'.stackVerts = s.stackVerts) (hts : s'.tstack = s.tstack) (h : s.RInvG dfs v d n) :
    s'.RInvG dfs v d n :=
  RInvG.congr (s := s) (s' := s') hg hsv hts (fun _ _ _ _ => by rw [hi]) (fun _ _ _ _ => by rw [hi])
    (fun _ _ _ _ _ => by rw [hi]) h

/-! ### Positional primitives: operations strictly above the bottom `n` -/

/-- Pushing the vertex entry of any `u` on top of the bottom `n`. -/
theorem RInvG.pushVert_pos {v d n : Nat} (u k : Nat) (hn : n ≤ s.tstack.length)
    (hown : ∀ t ∈ s.tstack, ∀ e, t.edges s.g s.items e → ¬ Items.EdgeBelow s.g s.items (vertItem u) e)
    (h : s.RInvG dfs v d n) : (after (pushVertTstack u k) s).RInvG dfs v d n := by
  show RInvG { s with tstack := ⟨u, k, s.nxtEdgeIdx, setSides s.stackDir[k]! [vertItem u] []⟩ :: s.tstack } dfs v d n
  have hed : ∀ e, (⟨u, k, s.nxtEdgeIdx, setSides s.stackDir[k]! [vertItem u] []⟩ : TEntry).edges
      s.g s.items e → Items.EdgeBelow s.g s.items (vertItem u) e := by
    rintro e ⟨i, hi, hb⟩
    rw [List.mem_singleton.1 ((mem_setSides _ _ i).1 hi)] at hb
    exact hb
  refine ⟨fun t ht hd hne => ?_, List.pairwise_cons.2 ⟨fun t' ht' e he hte => hown t' ht' e hte (hed e he), h.disj⟩⟩
  simp only [List.length_cons] at ht
  rw [Nat.succ_sub hn, List.drop_succ_cons] at ht
  exact EntryR.congr (s := s) rfl rfl (fun _ _ => rfl) (fun _ _ => rfl) (fun _ _ _ => Iff.rfl)
    (h.entries t ht hd hne)

/-- A merge of two entries above the bottom `n`. -/
theorem RInvG.mergeTop_pos {v d n : Nat} (hlen : n + 2 ≤ s.tstack.length) (h : s.RInvG dfs v d n) :
    (after mergeTstackTops s).RInvG dfs v d n ∧
      (after mergeTstackTops s).tstack.length + 1 = s.tstack.length := by
  match hts : s.tstack with
  | [] | [_] => rw [hts] at hlen; simp at hlen
  | a :: b :: rest =>
    refine ⟨h.mergeTop a b rest hts (fun hl => absurd hl (by omega)), ?_⟩
    show (mergeTstackTops.run s).2.tstack.length + 1 = _
    rw [mergeTstackTops_run_eq s a b rest hts]; rfl

/-- Re-targeting the top to any `u` above the bottom `n`. -/
theorem RInvG.retarget_pos {v d n : Nat} (u : Nat) (edgeDir : Bool) (t : TEntry) (rest : List TEntry)
    (hts : s.tstack = t :: rest) (hlen : n + 1 ≤ s.tstack.length) (h : s.RInvG dfs v d n) :
    (after (WalkState.retarget u edgeDir) s).RInvG dfs v d n := by
  show ((WalkState.retarget u edgeDir).run s).2.RInvG dfs v d n
  rw [retarget_run_eq u edgeDir s t rest hts]
  have hE : ∀ e, TEntry.edges s.g s.items
      { t with vStart := u, spans := setSides (!edgeDir) (t.spans.1 ++ t.spans.2) [] } e ↔
      t.edges s.g s.items e := by
    intro e; constructor
    · rintro ⟨i, hi, hb⟩; exact ⟨i, (mem_setSides _ _ i).1 hi, hb⟩
    · rintro ⟨i, hi, hb⟩; exact ⟨i, (mem_setSides _ _ i).2 hi, hb⟩
  have hd := h.disj; rw [hts] at hd
  obtain ⟨htop, hrest⟩ := List.pairwise_cons.1 hd
  have hlen' : rest.length + 1 = s.tstack.length := by rw [hts]; rfl
  refine ⟨fun x hx hdep hne => ?_, List.pairwise_cons.2 ⟨fun x hx e he => htop x hx e ((hE e).1 he), hrest⟩⟩
  simp only [List.length_cons] at hx
  rcases mem_drop_cons hx with ⟨hk, rfl⟩ | hx
  · exact absurd hk (by omega)
  · exact EntryR.congr (s := s) rfl rfl (fun _ _ => rfl) (fun _ _ => rfl) (fun _ _ _ => Iff.rfl)
      (h.entries x (by rw [hts]; simp only [List.length_cons]; exact mem_drop_cons_of hx) hdep hne)

/-- A merge loop whose every firing iterate has at least two entries above the bottom `n`. -/
theorem RInvG.mergeLoop {D v w d n : Nat} (cond : WalkM Bool) (fuel : Nat)
    (hcond : ∀ s, (cond.run s).2 = s) (hv : w < s.g.nv) (hi : s.Inv' D) (hs : Shape s)
    (hok : ∀ k, (∀ j, j ≤ k → (cond.run (iter mergeTstackTops j s)).1 = true) →
      MergeTopOk D (iter mergeTstackTops k s))
    (hlen : ∀ k, (∀ j, j < k → (cond.run (iter mergeTstackTops j s)).1 = true) →
      (cond.run (iter mergeTstackTops k s)).1 = true → n + 2 ≤ (iter mergeTstackTops k s).tstack.length)
    (h : s.RInvG dfs v d n) :
    ∃ k, after (loop fuel cond mergeTstackTops) s = iter mergeTstackTops k s ∧
      (∀ j, j < k → (cond.run (iter mergeTstackTops j s)).1 = true) ∧
      (iter mergeTstackTops k s).RInvG dfs v d n ∧ Step D w s (iter mergeTstackTops k s) := by
  induction fuel generalizing s with
  | zero => exact ⟨0, rfl, fun j hj => absurd hj (Nat.not_lt_zero j), h, Step.refl hi hs⟩
  | succ fuel ih =>
    unfold after
    rw [loop_succ_run fuel cond mergeTstackTops s (hcond s)]
    by_cases hc : (cond.run s).1 = true
    · simp only [hc, ↓reduceIte]
      have h0 : ∀ j, j ≤ 0 → (cond.run (iter mergeTstackTops j s)).1 = true :=
        fun j hj => by rw [Nat.le_zero.1 hj]; exact hc
      have st := Step.mergeTop (v := w) hi hs (hok 0 h0)
      have R := (h.mergeTop_pos (hlen 0 (fun j hj => absurd hj (Nat.not_lt_zero j)) hc)).1
      obtain ⟨k, heq, hj, R', st'⟩ := ih (by rw [st.g]; exact hv) st.inv st.shape
        (fun k hk => hok (k + 1) fun j hj => by
          cases j with
          | zero => exact hc
          | succ j => exact hk j (Nat.le_of_succ_le_succ hj))
        (fun k hk hck => hlen (k + 1) (fun j hj => by
          cases j with
          | zero => exact hc
          | succ j => exact hk j (Nat.lt_of_succ_lt_succ hj)) hck) R
      refine ⟨k + 1, heq, fun j hjk => ?_, R', st.trans st'⟩
      cases j with
      | zero => exact hc
      | succ j => exact hj j (Nat.lt_of_succ_lt_succ hjk)
    · simp only [hc, Bool.false_eq_true, ↓reduceIte]
      exact ⟨0, rfl, fun j hj => absurd hj (Nat.not_lt_zero j), h, Step.refl hi hs⟩

/-- Loop 2 above the bottom `n` (`hlen` from `Frontier.loop2`). -/
theorem RInvG.mergeLate_pos {D v w d n dl : Nat} (hv : w < s.g.nv) (hi : s.Inv' D) (hs : Shape s)
    (hok : MergeLateOk D dl s)
    (hlen : ∀ k, (∀ j, j < k →
      result (loop2Cond s.firstOccurrence[dl]!) (iter mergeTstackTops j s) = true) →
      result (loop2Cond s.firstOccurrence[dl]!) (iter mergeTstackTops k s) = true →
      n + 2 ≤ (iter mergeTstackTops k s).tstack.length)
    (h : s.RInvG dfs v d n) : (after (Spqr.mergeLate dl) s).RInvG dfs v d n := by
  show ((Spqr.mergeLate dl).run s).2.RInvG dfs v d n
  rw [mergeLate_run]
  by_cases hf : (curE s).firstIdx > s.firstOccurrence[dl]!
  · simp only [hf, ↓reduceIte]
    obtain ⟨k, heq, -, R, -⟩ := h.mergeLoop _ _ (fun _ => rfl) hv hi hs (hok.body hf) hlen
    have key : (after (loop s.tstack.length (loop2Cond s.firstOccurrence[dl]!) mergeTstackTops) s).RInvG
        dfs v d n := by rw [heq]; exact R
    exact key
  · simp only [hf, ↓reduceIte]; exact h

/-- Loop 3 keeps at least `origTstack + 3` entries: it only fires above that. -/
theorem loop3_iter_len (origTstack : Nat) : ∀ (k : Nat) (s : WalkState),
    (∀ j, j < k → result (loop3Cond origTstack) (iter mergeTstackTops j s) = true) →
    origTstack + 3 ≤ s.tstack.length → origTstack + 3 ≤ (iter mergeTstackTops k s).tstack.length
  | 0, _, _, h0 => h0
  | k + 1, s, hj, _ => by
    have hc : decide (s.tstack.length > origTstack + 3) = true := hj 0 (Nat.succ_pos k)
    have hc' := of_decide_eq_true hc
    show origTstack + 3 ≤ (iter mergeTstackTops k (mergeTstackTops.run s).2).tstack.length
    apply loop3_iter_len origTstack k (mergeTstackTops.run s).2
      (fun j hjk => hj (j + 1) (Nat.succ_lt_succ hjk))
    match hts : s.tstack with
    | [] | [_] => rw [hts] at hc'; simp at hc'
    | a :: b :: rest =>
      rw [mergeTstackTops_run_eq s a b rest hts]; rw [hts] at hc'
      simp only [List.length_cons] at hc' ⊢; omega

/-- `closeVert'` above the bottom `n`: loop 3 fires only above `origTstack + 3`, the unwrap, the
two merges, the re-targeting and the type-1 close all stay above `origTstack + 1 ≥ n + 2`. -/
theorem RInvG.closeVert_pos {D v w d n : Nat} {u : Nat} {edgeDir isType1 : Bool} {origTstack : Nat}
    {isSingle : Bool} (hv : w < s.g.nv) (hi : s.Inv' D) (hs : Shape s)
    (hok : CloseVertOk D u edgeDir isType1 origTstack isSingle s)
    (hn : n + 1 ≤ origTstack) (hlen : origTstack + 3 ≤ s.tstack.length)
    (h : s.RInvG dfs v d n) :
    (after (WalkState.closeVert' u edgeDir isType1 origTstack isSingle) s).RInvG dfs v d n ∧
      origTstack + 1 ≤ (after (WalkState.closeVert' u edgeDir isType1 origTstack isSingle) s).tstack.length := by
  have st₁ : Step D w s (cvS₁ isType1 origTstack isSingle s) := Step.vertPre hi hs hv hok.loop3
  have R₁ : (cvS₁ isType1 origTstack isSingle s).RInvG dfs v d n ∧
      origTstack + 3 ≤ (cvS₁ isType1 origTstack isSingle s).tstack.length := by
    cases isType1
    · show ((vertPre false origTstack isSingle).run s).2.RInvG dfs v d n ∧
        origTstack + 3 ≤ ((vertPre false origTstack isSingle).run s).2.tstack.length
      simp only [WalkState.vertPre, Bool.not_false, ↓reduceIte, WalkM.run_bind, run_tstackSize, WalkM.pure_run]
      obtain ⟨k, heq, hj, R, -⟩ := h.mergeLoop _ _ (fun _ => rfl) hv hi hs (hok.loop3 rfl)
        (fun k _ hck => by
          have hck' : decide ((iter mergeTstackTops k s).tstack.length > origTstack + 3) = true := hck
          have := of_decide_eq_true hck'; omega)
      have key : (after (loop s.tstack.length (loop3Cond origTstack) mergeTstackTops) s).RInvG dfs v d n ∧
          origTstack + 3 ≤ (after (loop s.tstack.length (loop3Cond origTstack) mergeTstackTops) s).tstack.length := by
        rw [heq]; exact ⟨R, loop3_iter_len origTstack k s hj hlen⟩
      exact key
    · exact ⟨h, hlen⟩
  obtain ⟨st₂, hfree⟩ := Step.vertUnwrap (v := w) st₁.inv st₁.shape (isType1 := isType1)
    (isSingle := cvB₁ isType1 origTstack isSingle s) (fun h => by subst h; exact hok.unwrap rfl)
  have st₂ : Step D w (cvS₁ isType1 origTstack isSingle s) (cvS₂ isType1 origTstack isSingle s) := st₂
  have R₂ : (cvS₂ isType1 origTstack isSingle s).RInvG dfs v d n ∧
      origTstack + 3 ≤ (cvS₂ isType1 origTstack isSingle s).tstack.length := by
    cases isType1
    · exact R₁
    · have hl₁ : origTstack + 3 ≤ s.tstack.length := R₁.2
      show ((some <$> maybeUnwrapNxt (if isSingle then NodeType.S else .R)).run s).2.RInvG dfs v d n ∧
        origTstack + 3 ≤ ((some <$> maybeUnwrapNxt (if isSingle then NodeType.S else .R)).run s).2.tstack.length
      rw [WalkM.map_run]
      refine ⟨R₁.1.unwrapNxt st₁.shape (hok.unwrap rfl) (fun hl => absurd hl (by omega)), ?_⟩
      match hts : s.tstack with
      | [] | [_] => rw [hts] at hl₁; simp at hl₁
      | a :: b :: rest =>
        obtain ⟨b', hts', -, -⟩ := maybeUnwrapNxt_tstack (if isSingle then NodeType.S else .R) a b rest hts
        have hts'' : ((maybeUnwrapNxt (if isSingle then NodeType.S else .R)).run s).2.tstack = a :: b' :: rest := hts'
        rw [hts'']; rw [hts] at hl₁; simpa using hl₁
  obtain ⟨R₃, hl₃⟩ := R₂.1.mergeTop_pos (by have := R₂.2; omega)
  have hl₃' : (cvS₃ isType1 origTstack isSingle s).tstack.length + 1 =
      (cvS₂ isType1 origTstack isSingle s).tstack.length := hl₃
  obtain ⟨R₄, hl₄⟩ := R₃.mergeTop_pos (by have := R₂.2; omega)
  have hl₄' : (cvS₄ isType1 origTstack isSingle s).tstack.length + 1 =
      (cvS₃ isType1 origTstack isSingle s).tstack.length := hl₄
  match hts₄ : (cvS₄ isType1 origTstack isSingle s).tstack with
  | [] => exact absurd hts₄ hok.retarget.nonempty
  | t :: rest =>
    have hl₄'' : (cvS₄ isType1 origTstack isSingle s).tstack.length = rest.length + 1 := by rw [hts₄]; rfl
    have R₅ : (cvS₅ u edgeDir isType1 origTstack isSingle s).RInvG dfs v d n :=
      R₄.retarget_pos u edgeDir t rest hts₄ (by have := R₂.2; omega)
    have hts₅ : (cvS₅ u edgeDir isType1 origTstack isSingle s).tstack =
        { t with vStart := u, spans := setSides (!edgeDir) (t.spans.1 ++ t.spans.2) [] } :: rest := by
      show ((WalkState.retarget u edgeDir).run _).2.tstack = _
      rw [retarget_run_eq u edgeDir _ t rest hts₄]
    have hl₅ : origTstack + 1 ≤ (cvS₅ u edgeDir isType1 origTstack isSingle s).tstack.length := by
      rw [hts₅]; simp only [List.length_cons]; have := R₂.2; omega
    cases isType1
    · exact ⟨R₅, hl₅⟩
    · have hf : ItemFree (cvS₂ true origTstack isSingle s)
          ((maybeUnwrapNxt (if isSingle then NodeType.S else .R)).run s).1 := hfree _ rfl
      have hf₅ := ((hf.merge).merge).retarget u edgeDir
      have hside := (hok.finish rfl).side
      have hcur : curE (cvS₅ u edgeDir true origTstack isSingle s) =
          { t with vStart := u, spans := setSides (!edgeDir) (t.spans.1 ++ t.spans.2) [] } := by
        rw [curE, hts₅]; rfl
      rw [hcur] at hside
      have hl₆ : ((finishTstackTop ((maybeUnwrapNxt (if isSingle then NodeType.S else .R)).run s).1).run
          (cvS₅ u edgeDir true origTstack isSingle s)).2.tstack.length =
          (cvS₅ u edgeDir true origTstack isSingle s).tstack.length := by
        rw [finishTstackTop_run_eq _ _ _ rest hts₅, hts₅]; rfl
      show ((finishTstackTop ((maybeUnwrapNxt (if isSingle then NodeType.S else .R)).run s).1).run
          (cvS₅ u edgeDir true origTstack isSingle s)).2.RInvG dfs v d n ∧
        origTstack + 1 ≤ ((finishTstackTop ((maybeUnwrapNxt (if isSingle then NodeType.S else .R)).run s).1).run
          (cvS₅ u edgeDir true origTstack isSingle s)).2.tstack.length
      refine ⟨R₅.finishTop _ _ rest hts₅ (fun hl => absurd hl (by omega)) hf₅.lt hf₅.root hf₅.free hside, ?_⟩
      rw [hl₆]; exact hl₅

/-- The P-check above the bottom `n`: when it fires, it unwraps and merges the two top entries and
closes the result, all above the bottom `n`; the stack loses at most one entry. -/
theorem RInvG.finishP_pos {D v d n : Nat} {u lowval : Nat} {isType1 : Bool} (hi : s.Inv' D) (hs : Shape s)
    (hok : FinishPOk D u lowval isType1 s) (hlen : n + 2 ≤ s.tstack.length) (h : s.RInvG dfs v d n) :
    (after (Spqr.finishP u lowval isType1) s).RInvG dfs v d n ∧
      s.tstack.length ≤ (after (Spqr.finishP u lowval isType1) s).tstack.length + 1 := by
  show ((Spqr.finishP u lowval isType1).run s).2.RInvG dfs v d n ∧
    s.tstack.length ≤ ((Spqr.finishP u lowval isType1).run s).2.tstack.length + 1
  simp only [Spqr.finishP, WalkM.run_bind, run_condP]
  by_cases hc : result (condP u lowval isType1) s = true
  · have h' : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == u) &&
        (s.tstack.tail.head!.topDepth == lowval)) = true := hc
    simp only [h', ↓reduceIte, WalkM.run_bind]
    obtain ⟨hu, hc2⟩ := hok.ok hc
    have r := maybeUnwrapNxt_spec (v := u) hi hs (by decide) hu
    have R₁ := h.unwrapNxt hs hu (fun hl => absurd hl (by omega))
    match hts : s.tstack with
    | [] | [_] => have := hu.two; rw [hts] at this; simp at this
    | a :: b :: rest =>
      obtain ⟨b', hts₁, -, -⟩ := maybeUnwrapNxt_tstack .P a b rest hts
      set s₁ := ((maybeUnwrapNxt .P).run s).2 with hs₁
      have hts₁' : s₁.tstack = a :: b' :: rest := hts₁
      have hl₁ : n + 2 ≤ s₁.tstack.length := by rw [hts₁']; rw [hts] at hlen; exact hlen
      obtain ⟨R₂, hl₂⟩ := R₁.mergeTop_pos hl₁
      have hl₂' : (mergeTstackTops.run s₁).2.tstack.length + 1 = s₁.tstack.length := hl₂
      have hts₂ : (mergeTstackTops.run s₁).2.tstack = TEntry.mergeInto a b' :: rest := by
        rw [mergeTstackTops_run_eq s₁ a b' rest hts₁']
      have hf := r.free.merge
      have hside := hc2.finish.side
      have hcur : curE (after mergeTstackTops (after (maybeUnwrapNxt .P) s)) = TEntry.mergeInto a b' := by
        show (mergeTstackTops.run s₁).2.tstack.head! = _
        rw [hts₂]; rfl
      rw [hcur] at hside
      refine ⟨R₂.finishTop _ (TEntry.mergeInto a b') rest hts₂ (fun hl => absurd hl (by omega))
        hf.lt hf.root hf.free hside, ?_⟩
      rw [finishTstackTop_run_eq _ _ _ rest hts₂]
      simp
  · have h' : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == u) &&
        (s.tstack.tail.head!.topDepth == lowval)) = false := Bool.eq_false_iff.2 hc
    simp only [h', Bool.false_eq_true, ↓reduceIte]
    exact ⟨h, Nat.le_succ _⟩

/-- The first-edge vertex push of any `u` above the bottom `n`. -/
theorem RInvG.finishTail_pos {v d n : Nat} {u k : Nat} {hasVert isSingle : Bool}
    (hown : hasVert = false → ∀ t ∈ s.tstack, ∀ e, t.edges s.g s.items e →
      ¬ Items.EdgeBelow s.g s.items (vertItem u) e)
    (hlen : hasVert = false → n + 1 ≤ s.tstack.length)
    (h : s.RInvG dfs v d n) : (after (Spqr.finishTail u k hasVert isSingle) s).RInvG dfs v d n := by
  cases hasVert
  · have R₁ := h.pushVert_pos u k (by have := hlen rfl; omega) (hown rfl)
    cases isSingle
    · exact (R₁.mergeTop_pos (by show n + 2 ≤ s.tstack.length + 1; have := hlen rfl; omega)).1
    · exact R₁
  · exact h

theorem RInvG.finishRest_pos {D v d n : Nat} {u k lowval : Nat} {isType1 hasVert isSingle : Bool}
    (hi : s.Inv' D) (hs : Shape s)
    (hok : FinishRestOk D u k lowval isType1 hasVert isSingle s)
    (hown : hasVert = false → ∀ t ∈ (after (Spqr.finishP u lowval isType1) s).tstack, ∀ e,
      t.edges (after (Spqr.finishP u lowval isType1) s).g (after (Spqr.finishP u lowval isType1) s).items e →
      ¬ Items.EdgeBelow (after (Spqr.finishP u lowval isType1) s).g
        (after (Spqr.finishP u lowval isType1) s).items (vertItem u) e)
    (hlen : n + 2 ≤ s.tstack.length) (h : s.RInvG dfs v d n) :
    (after (Spqr.finishRest u k lowval isType1 hasVert isSingle) s).RInvG dfs v d n := by
  obtain ⟨R₁, hl₁⟩ := h.finishP_pos hi hs hok.p hlen
  have hl : n + 1 ≤ (after (Spqr.finishP u lowval isType1) s).tstack.length := by omega
  exact R₁.finishTail_pos hown (fun _ => hl)

/-! ### The frame -/

/-- A `finishEdge` at `(curV, d)` whose boundary `origTstack` lies above the bottom `n₀` preserves
`RInvG dfs p dp n₀`. Tree edge: `feS₀` is free of the stack (`FinishRShape.pend`), loop 1 by
`RInvG.closeEars` with `Frontier.loop1`, loop 2 with `Frontier.loop2`, `closeVert'` and the tail
above `origTstack + 3` (`hclose`) resp. `origTstack + 1`. Back edge: the vertex entry of `curV`
is on top of the whole base (`hback`), so the push and the P-check stay above `n₀ + 1`. -/
theorem finishEdge_rInvG_base {D : Nat} (p dp curV d lv : Nat) (kind : RetKind) (o : DfsOut)
    (origTstack n₀ : Nat) (hasVert : Bool) (ho : o.cls = .ret lv kind) (hlow : lv < d)
    (hv : curV < s.g.nv) (hi : s.Inv' D) (hs : Shape s)
    (hok : FinishOk D curV d lv o origTstack hasVert s)
    (hfront : Frontier (o := o) d origTstack s)
    (hshape : kind ≠ .backEdge → FinishRShape dfs curV d o origTstack hasVert s)
    (hback : kind = .backEdge → hasVert = true ∧ s.tstack.length ≤ origTstack)
    (hclose : kind ≠ .backEdge → origTstack + 3 ≤ (feS₂ d o s).tstack.length)
    (hn : n₀ ≤ origTstack) (hnv : hasVert = true → n₀ + 1 ≤ origTstack)
    (hR : s.RInvG dfs p dp n₀) :
    (after (finishEdge curV d o origTstack hasVert) s).RInvG dfs p dp n₀ := by
  have hlv : o.cls.lowval d = lv := by rw [ho]; rfl
  have hge : ¬ (lv ≥ d) := by omega
  have hlow' : o.cls.lowval d < d := by rw [hlv]; exact hlow
  have hsub : subEdges o o.e := (hfront.owns o.e hok.e_lt).1 (Or.inl rfl)
  have st₀ : Step D curV s (feS₀ d o s) :=
    Step.modifyVs hi hs (edgeItem s.g o.e) _ (by show 1 + s.g.nv + o.e < _; have := hok.e_lt; omega)
  have hv₀ : curV < (feS₀ d o s).g.nv := by rw [st₀.g]; exact hv
  set f : Item → Item := fun it =>
    { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) } with hf
  have hn₀ : origTstack ≤ (feS₀ d o s).tstack.length := hfront.size
  by_cases hk : kind = .backEdge
  · subst hk
    obtain ⟨hhv, hlen⟩ := hback rfl
    subst hhv
    have ht' : o.cls.isTree = false := by rw [ho]; rfl
    have h1 : o.cls.isType1 = true := by rw [ho]; rfl
    have hown : ∀ t ∈ s.tstack, ¬ t.edges s.g s.items o.e := fun t ht hte =>
      hfront.base_disj t (by rw [Nat.sub_eq_zero_of_le hlen, List.drop_zero]; exact ht) o.e hok.e_lt hte hsub
    have hfree₀ : ∀ t ∈ s.tstack, edgeItem s.g o.e ∉ t.spans.1 ++ t.spans.2 := fun t ht hmem =>
      hown t ht ⟨_, hmem, .refl⟩
    have R₀ : (feS₀ d o s).RInvG dfs p dp n₀ :=
      hR.modifyVs_free (edgeItem s.g o.e) f (fun _ => rfl) (fun _ => rfl) hfree₀
    have hq₀ : Items.ch (feS₀ d o s).items (edgeItem (feS₀ d o s).g o.e) = [] := by
      show Items.ch (s.items.modify _ f) (edgeItem s.g o.e) = []
      rw [Items.ch_modify_ch_eq (edgeItem s.g o.e) f (fun _ => rfl)]; exact hok.q ht'
    have hown₀ : ∀ t ∈ (feS₀ d o s).tstack, ¬ t.edges (feS₀ d o s).g (feS₀ d o s).items o.e :=
      fun t ht hte => hown t ht ((TEntry.edges_congr (fun i _ e =>
        Items.Below_modify_ch_eq (edgeItem s.g o.e) f (fun _ => rfl)) o.e).1 hte)
    have st₁ : Step D curV _ (after (pushEdgeTstack curV lv o.e) (feS₀ d o s)) :=
      Step.pushEdge st₀.inv st₀.shape curV lv o.e hok.e_lt hq₀ (hok.ends ht') (hok.lv_le ht')
    have R₁ : (after (pushEdgeTstack curV lv o.e) (feS₀ d o s)).RInvG dfs p dp n₀ :=
      R₀.pushEdge curV lv o.e (by omega) hq₀ hown₀
    have st₂ : Step D curV _ (feBack curV lv d o s) := Step.frame st₁.inv st₁.shape rfl rfl rfl rfl
    have R₂ : (feBack curV lv d o s).RInvG dfs p dp n₀ :=
      RInvG.of_eq' (s := after (pushEdgeTstack curV lv o.e) (feS₀ d o s)) rfl rfl rfl rfl R₁
    have hl₂ : (feBack curV lv d o s).tstack.length = s.tstack.length + 1 := rfl
    have hsz := hfront.size
    show ((finishEdge curV d o origTstack true).run s).2.RInvG dfs p dp n₀
    rw [finishEdge_eq]
    simp only [finishEdge', finishBack, hlv, hge, ↓reduceIte, WalkM.run_bind, WalkM.get_run,
      run_stackDir, run_makeVs, run_modifyItem, ht', Bool.false_eq_true, h1]
    have hrest := hok.rest_back ht'
    rw [h1] at hrest
    exact RInvG.finishRest_pos (s := feBack curV lv d o s) st₂.inv st₂.shape hrest (fun h => nomatch h)
      (by have := hnv rfl; omega) R₂
  · have hshape := hshape hk
    have hclose := hclose hk
    have ht : o.cls.isTree = true := by
      rw [ho]; cases kind with
      | backEdge => exact absurd rfl hk
      | type1Child => rfl
      | type2Child => rfl
    have hfree₀ : ∀ t ∈ s.tstack, edgeItem s.g o.e ∉ t.spans.1 ++ t.spans.2 := fun t ht hmem =>
      hshape.pend t ht ⟨_, hmem, .refl⟩
    have R₀ : (feS₀ d o s).RInvG dfs p dp n₀ :=
      hR.modifyVs_free (edgeItem s.g o.e) f (fun _ => rfl) (fun _ => rfl) hfree₀
    have hown₀ : ∀ t ∈ (feS₀ d o s).tstack, ¬ t.edges (feS₀ d o s).g (feS₀ d o s).items o.e :=
      fun t ht hte => hshape.pend t ht ((TEntry.edges_congr (fun i _ e =>
        Items.Below_modify_ch_eq (edgeItem s.g o.e) f (fun _ => rfl)) o.e).1 hte)
    obtain ⟨k, heq, -, -, R₁, -⟩ := R₀.closeEars st₀.inv st₀.shape hv₀ (by omega) (hok.ears ht) hown₀
      (fun k hk hck => Nat.le_trans (Nat.add_le_add_right hn _) ((hfront.loop1 ht hlow' k hk).2 hck))
    have st₁ : Step D curV _ (feS₁ d o s) := Step.closeEars st₀.inv st₀.shape hv₀ (hok.ears ht)
    have hv₁ : curV < (feS₁ d o s).g.nv := by rw [st₁.g]; exact hv₀
    have R₁' : (feS₁ d o s).RInvG dfs p dp n₀ := by
      have heq' : feS₁ d o s =
          iter (Spqr.loop1Body d s.stackDir[d]!) k (ceS₁ o.dest d o.e (feS₀ d o s)) := heq
      rw [heq']; exact R₁
    have R₂ : (feS₂ d o s).RInvG dfs p dp n₀ := R₁'.mergeLate_pos hv₁ st₁.inv st₁.shape (hok.late ht)
      (fun k hk hck => Nat.le_trans (Nat.add_le_add_right hn 2) ((hfront.loop2 ht hlow' k hk).2 hck))
    have st₂ : Step D curV _ (feS₂ d o s) := Step.mergeLate st₁.inv st₁.shape hv₁ (hok.late ht)
    have hv₂ : curV < (feS₂ d o s).g.nv := by rw [st₂.g]; exact hv₁
    show ((finishEdge curV d o origTstack hasVert).run s).2.RInvG dfs p dp n₀
    rw [finishEdge_eq]
    simp only [finishEdge', finishTree, hlv, hge, ↓reduceIte, ht, closeVert_eq]
    cases hasVert
    · simp only [Bool.false_eq_true, ↓reduceIte]
      have hown := hshape.vert_own rfl
      simp only [feP, hlv] at hown
      exact RInvG.finishRest_pos (s := feS₂ d o s) st₂.inv st₂.shape (hok.rest_tree ht rfl) (fun _ => hown)
        (by omega) R₂
    · simp only [↓reduceIte]
      have st₃ : Step D curV _ (feS₃ curV d o origTstack s) :=
        Step.closeVert' st₂.inv st₂.shape hv₂ (hok.vert ht rfl)
      obtain ⟨R₃, hl₃⟩ := R₂.closeVert_pos hv₂ st₂.inv st₂.shape (hok.vert ht rfl) (hnv rfl) hclose
      have hl₃' : origTstack + 1 ≤ (feS₃ curV d o origTstack s).tstack.length := hl₃
      exact RInvG.finishRest_pos (s := feS₃ curV d o origTstack s) st₃.inv st₃.shape (hok.rest_vert ht rfl)
        (fun h => nomatch h) (by have := hnv rfl; omega) R₃

end WalkState
end Spqr
