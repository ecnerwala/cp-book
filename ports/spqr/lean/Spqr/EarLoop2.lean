import Spqr.EarLoop1

/-!
# Loop 2 from `EarFinish.late`

`mergeLate` folds the top entries into `cur` while `firstIdx > firstOccurrence[d]` (`loop2Cond`).
`iter_merge_eq` computes the `k`-th iterate of `mergeTstackTops` on `c₀ :: R` as
`l2Cur c₀ (R.take k) :: R.drop k`; `mergeLateOk_of_late` derives `MergeLateOk` from `EarLate`: the
condition at iterate `j` is `firstIdx` of the `j`-th merged entry (`mergeInto` keeps `nxt.firstIdx`),
so every merged entry passed it and `EarLate.merge` applies; the rest is pairwise disjoint.
-/
namespace Spqr
open WalkM WalkState

theorem l2Cur_cons (c₀ t : TEntry) (R : List TEntry) :
    l2Cur c₀ (t :: R) = l2Cur (TEntry.mergeInto c₀ t) R := rfl

theorem l2Cur_concat (c₀ : TEntry) (done : List TEntry) (t : TEntry) :
    l2Cur c₀ (done ++ [t]) = TEntry.mergeInto (l2Cur c₀ done) t := by
  simp [l2Cur, List.foldl_append]

theorem mergeInto_firstIdx (cur nxt : TEntry) : (TEntry.mergeInto cur nxt).firstIdx = nxt.firstIdx := rfl
theorem mergeInto_vStart (cur nxt : TEntry) : (TEntry.mergeInto cur nxt).vStart = nxt.vStart := rfl

theorem l2Cur_edges {g : Graph} {items : Items} (c₀ : TEntry) (done : List TEntry) (e : Nat) :
    (l2Cur c₀ done).edges g items e ↔ c₀.edges g items e ∨ ∃ t ∈ done, t.edges g items e := by
  induction done generalizing c₀ with
  | nil => simp [l2Cur]
  | cons t R ih => rw [l2Cur_cons, ih, TEntry.edges_mergeInto]; simp [or_assoc]

theorem iter_add (body : WalkM Unit) (m n : Nat) (s : WalkState) :
    iter body (n + m) s = iter body n (iter body m s) := by
  induction m generalizing s with
  | zero => rfl
  | succ m ih => exact ih _

theorem mergeTstackTops_short (s : WalkState) (h : s.tstack.length ≤ 1) :
    (mergeTstackTops.run s).2.tstack.length ≤ 1 := by
  obtain ⟨g, tern, items, sv, sd, nei, fo, tstack, tb, tsl⟩ := s
  match tstack, h with
  | [], _ => exact Nat.zero_le _
  | [_], _ => exact Nat.zero_le _

theorem iter_merge_short (s : WalkState) (h : s.tstack.length ≤ 1) (k : Nat) :
    (iter mergeTstackTops k s).tstack.length ≤ 1 := by
  induction k generalizing s with
  | zero => exact h
  | succ k ih => exact ih _ (mergeTstackTops_short s h)

theorem iter_merge_eq (s : WalkState) (c₀ : TEntry) (R : List TEntry) (hs : s.tstack = c₀ :: R) (k : Nat)
    (hk : k ≤ R.length) :
    iter mergeTstackTops k s = { s with tstack := l2Cur c₀ (R.take k) :: R.drop k } := by
  induction k generalizing s c₀ R with
  | zero => show s = { s with tstack := c₀ :: R }; rw [← hs]
  | succ k ih =>
    cases R with
    | nil => simp at hk
    | cons t R =>
      show iter mergeTstackTops k (mergeTstackTops.run s).2 = _
      rw [mergeTstackTops_run_eq s c₀ t R hs, ih _ (TEntry.mergeInto c₀ t) R rfl (by simpa using hk)]
      rfl

theorem result_loop2Cond_iff (fo : Nat) (st : WalkState) :
    result (loop2Cond fo) st = true ↔ fo < st.tstack.head!.firstIdx := by
  unfold result; rw [run_loop2Cond]; simp

/-- Every entry merged so far passed `loop2Cond`. -/
theorem late_done_gt (s : WalkState) (c₀ : TEntry) (R : List TEntry) (hs : s.tstack = c₀ :: R) (fo : Nat)
    (k : Nat) (hk : k ≤ R.length)
    (hc : ∀ j, j ≤ k → result (loop2Cond fo) (iter mergeTstackTops j s) = true) :
    ∀ u ∈ R.take k, fo < u.firstIdx := by
  induction k with
  | zero => simp
  | succ k ih =>
    rw [List.take_succ_eq_append_getElem (by omega)]
    intro u hu
    rw [List.mem_append, List.mem_singleton] at hu
    rcases hu with hu | rfl
    · exact ih (by omega) (fun j hj => hc j (by omega)) u hu
    · have := (result_loop2Cond_iff fo _).1 (hc (k + 1) (Nat.le_refl _))
      rw [iter_merge_eq s c₀ R hs (k + 1) hk, List.take_succ_eq_append_getElem (by omega),
        l2Cur_concat] at this
      exact this

theorem mergeLateOk_of_late {d : Nat} {s : WalkState} {c₀ : TEntry} {R : List TEntry}
    (hL : EarLate d s c₀ R) : MergeLateOk (d + 1) d s := by
  refine ⟨fun hfo k hc => ?_⟩
  have hc₀ : s.firstOccurrence[d]! < c₀.firstIdx := by simpa [curE, hL.tstack] using hfo
  by_cases hk : k ≤ R.length
  · rw [iter_merge_eq s c₀ R hL.tstack k hk]
    intro cur nxt rest h
    simp only [List.cons.injEq] at h
    obtain ⟨rfl, hd⟩ := h
    have hR : R = R.take k ++ nxt :: rest := by rw [← hd, List.take_append_drop]
    have hdone := late_done_gt s c₀ R hL.tstack _ k hk hc
    refine ⟨MergeOk.congr (s := s) rfl rfl (fun _ _ => Iff.rfl) (fun _ _ => Iff.rfl) (hL.merge (R.take k) nxt rest hR hc₀ hdone), ?_⟩
    intro t ht e he hte hor
    have hp := hL.disj
    rw [hL.tstack, hR, List.pairwise_cons, List.pairwise_append, List.pairwise_cons] at hp
    obtain ⟨h0, -, ⟨hn, -⟩, hcross⟩ := hp
    rw [l2Cur_edges] at hor
    rcases hor with (h | ⟨u, hu, h⟩) | h
    · exact h0 t (by simp [ht]) e he h hte
    · exact hcross u hu t (by simp [ht]) e he h hte
    · exact hn t ht e he h hte
  · intro cur nxt rest h
    have h1 : (iter mergeTstackTops R.length s).tstack.length ≤ 1 := by
      rw [iter_merge_eq s c₀ R hL.tstack R.length (Nat.le_refl _)]; simp
    have h2 := iter_merge_short _ h1 (k - R.length)
    rw [← iter_add, Nat.sub_add_cancel (by omega), h] at h2
    simp at h2


theorem getSide_setSides_self {α : Type} (dir : Bool) (a b : α) : getSide (setSides dir a b) dir = a := by
  cases dir <;> rfl

/-- `maybeUnwrapNxt` on a one-sided single-item `nxt` keeps the stack shape and every entry's edge
set (the unwrapped entry holds the children of its item, hence the same edges). -/
theorem maybeUnwrapNxt_edges {ty : NodeType} (hs : Shape s) (hty : ty ∉ [NodeType.F, .V, .Q])
    {a b : TEntry} {rest : List TEntry} (hts : s.tstack = a :: b :: rest) {i : ItemId}
    (hside : getSide b.spans (!s.stackDir[b.topDepth]!) = [])
    (hsingle : getSide b.spans s.stackDir[b.topDepth]! = [i]) :
    ∃ b', (after (maybeUnwrapNxt ty) s).tstack = a :: b' :: rest ∧
      (after (maybeUnwrapNxt ty) s).g = s.g ∧ (after (maybeUnwrapNxt ty) s).stackVerts = s.stackVerts ∧
      (after (maybeUnwrapNxt ty) s).stackDir = s.stackDir ∧
      b'.vStart = b.vStart ∧ b'.topDepth = b.topDepth ∧
      getSide b'.spans (!s.stackDir[b.topDepth]!) = [] ∧
      (∀ e, e < s.g.ne → (b'.edges s.g (after (maybeUnwrapNxt ty) s).items e ↔ b.edges s.g s.items e)) ∧
      (∀ t : TEntry, ∀ e, t.edges s.g (after (maybeUnwrapNxt ty) s).items e ↔ t.edges s.g s.items e) := by
  have hh : i = (getSide b.spans s.stackDir[b.topDepth]!).head! := by rw [hsingle]; rfl
  unfold after
  rw [maybeUnwrapNxt_run_eq ty s a b rest hts _ rfl i hh]
  have hpush : ∀ t : TEntry, ∀ e, t.edges s.g (s.items.push ⟨ty, (none, none), []⟩) e ↔ t.edges s.g s.items e :=
    fun _ e => TEntry.edges_congr (fun _ _ _ => Items.Below_push_nil _ rfl) e
  by_cases h1 : ty = .R ∨ s.ternarize = true
  · rw [if_pos h1, run_allocItem]
    exact ⟨b, hts, rfl, rfl, rfl, rfl, rfl, hside, fun e _ => hpush b e, hpush⟩
  rw [if_neg h1]
  by_cases h2 : s.items[i]!.type = ty
  · rw [if_pos h2]
    refine ⟨_, rfl, rfl, rfl, rfl, rfl, rfl, getSide_setSides_not _ _, fun e he => ?_, fun _ _ => Iff.rfl⟩
    have hib : i ∈ b.spans.1 ++ b.spans.2 := by
      rw [mem_of_getSide_nil _ b.spans hside, hsingle]; exact List.mem_singleton_self i
    have hilt : i < s.items.size := hs.span b (by simp [hts]) i hib
    have hget : s.items[i]! = s.items[i] := getElem!_pos s.items i hilt
    have hity : Items.type s.items i = ty := by rw [Items.type_eq_getElem hilt, ← hget]; exact h2
    have hie : ∀ e, e < s.g.ne → edgeItem s.g e ≠ i := fun e he =>
      hs.edgeItem_ne he (hs.node_of_type hilt hty hity)
    have hch : s.items[i]!.ch = Items.ch s.items i := by rw [Items.ch_eq_getElem hilt, hget]
    rw [TEntry.edges_single _ i hside hsingle, hch]
    exact TEntry.edges_unwrap _ b.vStart b.topDepth b.firstIdx i hie he
  · rw [if_neg h2, run_allocItem]
    exact ⟨b, hts, rfl, rfl, rfl, rfl, rfl, hside, fun e _ => hpush b e, hpush⟩


/-! ### The vertex close (`closeVert'`) from `EarClose` -/

/-- A `loop` with a state-independent condition stops at the first iterate failing it, or at the fuel. -/
theorem loop_run_iter (fuel : Nat) (cond : WalkM Bool) (body : WalkM Unit) (s : WalkState)
    (hcond : ∀ s, (cond.run s).2 = s) :
    ∃ k, k ≤ fuel ∧ ((loop fuel cond body).run s).2 = iter body k s ∧
      (∀ j, j < k → (cond.run (iter body j s)).1 = true) ∧
      (k = fuel ∨ (cond.run (iter body k s)).1 = false) := by
  induction fuel generalizing s with
  | zero => exact ⟨0, Nat.le_refl _, rfl, fun j hj => absurd hj (Nat.not_lt_zero _), .inl rfl⟩
  | succ fuel ih =>
    rw [loop_succ_run fuel cond body s (hcond s)]
    by_cases hc : (cond.run s).1 = true
    · rw [if_pos hc]
      obtain ⟨k, hk, heq, hall, hend⟩ := ih (body.run s).2
      refine ⟨k + 1, by omega, heq, fun j hj => ?_, ?_⟩
      · cases j with
        | zero => exact hc
        | succ j => exact hall j (by omega)
      · rcases hend with h | h
        · exact .inl (by omega)
        · exact .inr h
    · rw [if_neg hc]
      exact ⟨0, Nat.zero_le _, rfl, fun j hj => absurd hj (Nat.not_lt_zero _), .inr (Bool.eq_false_iff.2 hc)⟩

/-- `MergeTopOk` at the `k`-th merge of a fold over `R'` (`k < R'.length`), from `FoldSpec` and the
pairwise edge-disjointness of the stack `c₀ :: R' ++ tail`. -/
theorem mergeTopOk_iter_fold {D : Nat} {s : WalkState} {c₀ : TEntry} {R' tail : List TEntry}
    (hs : s.tstack = c₀ :: (R' ++ tail))
    (hdisj : s.tstack.Pairwise fun t t' => ∀ e, e < s.g.ne → t.edges s.g s.items e → ¬ t'.edges s.g s.items e)
    (hfold : FoldSpec D s c₀ R') (k : Nat) (hk : k < R'.length) :
    MergeTopOk D (iter mergeTstackTops k s) := by
  rw [iter_merge_eq s c₀ _ hs k (by simp; omega)]
  intro cur nxt rest h
  simp only [List.cons.injEq] at h
  obtain ⟨rfl, hd⟩ := h
  have hdrop : R'.drop k = R'[k] :: R'.drop (k + 1) := List.drop_eq_getElem_cons hk
  rw [List.drop_append_of_le_length (Nat.le_of_lt hk), hdrop] at hd
  simp only [List.cons_append, List.cons.injEq] at hd
  obtain ⟨rfl, rfl⟩ := hd
  rw [List.take_append_of_le_length (Nat.le_of_lt hk)]
  have hR : R' = R'.take k ++ R'[k] :: R'.drop (k + 1) := by
    conv_lhs => rw [← List.take_append_drop k R']
    rw [hdrop]
  refine ⟨MergeOk.congr (s := s) rfl rfl (fun _ _ => Iff.rfl) (fun _ _ => Iff.rfl) (hfold _ _ _ hR), ?_⟩
  have hsplit := List.pairwise_split hdisj (l₁ := c₀ :: (R'.take k ++ [R'[k]])) (l₂ := R'.drop (k + 1) ++ tail)
    (by rw [hs]; (conv_lhs => rw [hR]); simp only [List.append_assoc, List.cons_append, List.nil_append])
  intro t ht e he hte hor
  rw [l2Cur_edges] at hor
  rcases hor with (h | ⟨u, hu, h⟩) | h
  · exact hsplit c₀ List.mem_cons_self t ht e he h hte
  · exact hsplit u (List.mem_cons_of_mem _ (List.mem_append_left _ hu)) t ht e he h hte
  · exact hsplit R'[k] (List.mem_cons_of_mem _ (List.mem_append_right _ (List.mem_singleton_self _))) t ht e he h hte

/-- `MergeOk` transports along entries with the same `vStart`s, merged `topDepth` and edge sets. -/
theorem MergeOk.congr_entry {D : Nat} {s s' : WalkState} {cur nxt cur' nxt' : TEntry}
    (hg : s'.g = s.g) (hsv : s'.stackVerts = s.stackVerts)
    (hcv : cur'.vStart = cur.vStart) (hnv : nxt'.vStart = nxt.vStart)
    (htd : (TEntry.mergeInto cur' nxt').topDepth = (TEntry.mergeInto cur nxt).topDepth)
    (hc : ∀ e, e < s.g.ne → (cur'.edges s.g s'.items e ↔ cur.edges s.g s.items e))
    (hn : ∀ e, e < s.g.ne → (nxt'.edges s.g s'.items e ↔ nxt.edges s.g s.items e))
    (h : MergeOk D s cur nxt) : MergeOk D s' cur' nxt' := by
  refine ⟨?_, ?_⟩
  · rw [hg]
    rintro ⟨e, he, hce⟩ ⟨e', he', hne⟩
    obtain ⟨x, ⟨e₁, he₁, h₁, hi₁⟩, ⟨e₂, he₂, h₂, hi₂⟩⟩ := h.share ⟨e, he, (hc e he).1 hce⟩ ⟨e', he', (hn e' he').1 hne⟩
    exact ⟨x, ⟨e₁, he₁, (hc e₁ he₁).2 h₁, hi₁⟩, ⟨e₂, he₂, (hn e₂ he₂).2 h₂, hi₂⟩⟩
  · rw [hcv, TEntry.Term_congr hsv (show (TEntry.mergeInto cur' nxt').vStart = (TEntry.mergeInto cur nxt).vStart from hnv) htd]
    rcases h.bottom with hb | hb
    · exact .inl hb
    · refine .inr ?_
      rw [hg]
      intro e he hi
      rcases hb e he hi with h₁ | h₁
      · exact .inl ((hc e he).2 h₁)
      · exact .inr ((hn e he).2 h₁)

/-- The type-2 vertex close: loop 3 and the two merges fold the sub-stack `mid ++ [py, vy]` into `c`. -/
theorem closeVertOk_type2 {curV d l : Nat} {o : DfsOut} {base : List TEntry} {s st : WalkState}
    {c : TEntry} {mid : List TEntry} {py vy : TEntry} (isSingle edgeDir : Bool)
    (hC : EarClose curV d l o true base s st c mid py vy) :
    CloseVertOk (d + 1) curV edgeDir false base.length isSingle st := by
  set R := mid ++ [py, vy] ++ base with hR
  have hts : st.tstack = c :: R := hC.tstack
  have hRl : R.length = mid.length + 2 + base.length := by simp [hR]; omega
  have hdisj : st.tstack.Pairwise fun t t' => ∀ e, e < st.g.ne → t.edges st.g st.items e → ¬ t'.edges st.g st.items e := by
    rw [hC.g]; exact hC.disj
  have hsplit : ∀ l₁ l₂, st.tstack = l₁ ++ l₂ → ∀ u ∈ l₁, ∀ t ∈ l₂, ∀ e, e < st.g.ne →
      u.edges st.g st.items e → ¬ t.edges st.g st.items e :=
    fun l₁ l₂ hl u hu t ht => List.pairwise_split hdisj hl u hu t ht
  have hiter : ∀ k, k ≤ R.length → (iter mergeTstackTops k st).tstack.length = R.length + 1 - k := fun k hk => by
    rw [iter_merge_eq st c R hts k hk]; simp [List.length_drop]; omega
  have hcond : ∀ k, result (loop3Cond base.length) (iter mergeTstackTops k st) = true ↔
      base.length + 3 < (iter mergeTstackTops k st).tstack.length := fun k => by
    unfold result; rw [run_loop3Cond]; simp
  have hlt : ∀ k, result (loop3Cond base.length) (iter mergeTstackTops k st) = true → k < mid.length := fun k hk => by
    rw [hcond] at hk
    by_cases hkR : k ≤ R.length
    · rw [hiter k hkR] at hk; omega
    · have h1 : (iter mergeTstackTops R.length st).tstack.length ≤ 1 := by rw [hiter _ (Nat.le_refl _)]; omega
      have h2 := iter_merge_short _ h1 (k - R.length)
      rw [← iter_add, Nat.sub_add_cancel (by omega)] at h2
      omega
  have hmto : ∀ k, k < mid.length → MergeTopOk (d + 1) (iter mergeTstackTops k st) := fun k hk =>
    mergeTopOk_iter_fold (R' := mid ++ [py, vy]) (tail := base) hts hdisj hC.fold k (by simp; omega)
  have hs₁ : cvS₁ false base.length isSingle st = iter mergeTstackTops mid.length st := by
    have : cvS₁ false base.length isSingle st =
        ((loop st.tstack.length (loop3Cond base.length) mergeTstackTops).run st).2 := by
      simp only [cvS₁, after, vertPre, Bool.not_false, ↓reduceIte, WalkM.run_bind, run_tstackSize, WalkM.pure_run]
    rw [this]
    obtain ⟨k, hk, heq, hall, hend⟩ :=
      loop_run_iter st.tstack.length (loop3Cond base.length) mergeTstackTops st (fun _ => rfl)
    rw [heq]
    have hk1 : k ≤ mid.length := by
      by_contra h
      have := hlt mid.length (hall mid.length (by omega))
      omega
    have hk2 : mid.length ≤ k := by
      rcases hend with h | h
      · rw [hts] at h; simp at h; omega
      · by_contra h'
        have : result (loop3Cond base.length) (iter mergeTstackTops k st) = true := by
          rw [hcond, hiter k (by omega)]; omega
        have h' : result (loop3Cond base.length) (iter mergeTstackTops k st) = false := h
        rw [this] at h'; cases h'
    rw [Nat.le_antisymm hk1 hk2]
  have hs₂ : cvS₂ false base.length isSingle st = { st with tstack := l2Cur c mid :: py :: vy :: base } := by
    show cvS₁ false base.length isSingle st = _
    rw [hs₁, iter_merge_eq st c R hts mid.length (by omega)]
    simp only [hR, List.append_assoc]
    rw [List.take_left, List.drop_left]; rfl
  have hs₃ : cvS₃ false base.length isSingle st = { st with tstack := l2Cur c (mid ++ [py]) :: vy :: base } := by
    show (mergeTstackTops.run (cvS₂ false base.length isSingle st)).2 = _
    rw [hs₂, mergeTstackTops_run_eq _ _ _ _ rfl, l2Cur_concat]
  have hs₄ : cvS₄ false base.length isSingle st = { st with tstack := l2Cur c (mid ++ [py, vy]) :: base } := by
    show (mergeTstackTops.run (cvS₃ false base.length isSingle st)).2 = _
    rw [hs₃, mergeTstackTops_run_eq _ _ _ _ rfl]
    simp only [show mid ++ [py, vy] = (mid ++ [py]) ++ [vy] by simp, l2Cur_concat]
  have hR₁ : st.tstack = (c :: (mid ++ [py])) ++ (vy :: base) := by rw [hts, hR]; simp
  have hR₂ : st.tstack = (c :: (mid ++ [py, vy])) ++ base := by rw [hts, hR]; simp
  refine ⟨fun _ k hk => hmto k (hlt k (hk k (Nat.le_refl _))), fun h => Bool.noConfusion h, ?_, ?_, ?_, fun h => Bool.noConfusion h⟩
  · rw [hs₂]
    intro cur nxt rest h
    simp only [List.cons.injEq] at h
    obtain ⟨rfl, rfl, rfl⟩ := h
    refine ⟨MergeOk.congr (s := st) rfl rfl (fun _ _ => Iff.rfl) (fun _ _ => Iff.rfl) (hC.fold mid py [vy] (by simp)), ?_⟩
    intro t ht e he hte hor
    rw [l2Cur_edges] at hor
    rcases hor with (h | ⟨u, hu, h⟩) | h
    · exact hsplit _ _ hR₁ c (by simp) t ht e he h hte
    · exact hsplit _ _ hR₁ u (by simp [hu]) t ht e he h hte
    · exact hsplit _ _ hR₁ py (by simp) t ht e he h hte
  · rw [hs₃]
    intro cur nxt rest h
    simp only [List.cons.injEq] at h
    obtain ⟨rfl, rfl, rfl⟩ := h
    refine ⟨MergeOk.congr (s := st) rfl rfl (fun _ _ => Iff.rfl) (fun _ _ => Iff.rfl) (hC.fold (mid ++ [py]) vy [] (by simp)), ?_⟩
    intro t ht e he hte hor
    rw [l2Cur_edges] at hor
    rcases hor with (h | ⟨u, hu, h⟩) | h
    · exact hsplit _ _ hR₂ c (by simp) t ht e he h hte
    · exact hsplit _ _ hR₂ u (List.mem_cons_of_mem _ (by
        simp only [List.mem_append, List.mem_cons, List.not_mem_nil, or_false] at hu ⊢
        rcases hu with h | h <;> simp [h])) t ht e he h hte
    · exact hsplit _ _ hR₂ vy (by simp) t ht e he h hte
  · rw [hs₄]
    refine ⟨by simp, ?_, ?_⟩
    · have hmv : (l2Cur c (mid ++ [py, vy])).vStart = vy.vStart := by
        rw [show mid ++ [py, vy] = (mid ++ [py]) ++ [vy] by simp, l2Cur_concat]; rfl
      refine .inr (.inl ?_)
      show st.g.Interior ((l2Cur c (mid ++ [py, vy])).edges st.g st.items) (l2Cur c (mid ++ [py, vy])).vStart
      rw [hmv, hC.g]
      intro e he hinc
      obtain ⟨t, ht, hte⟩ := hC.sub_cover e he (hC.y_edges e he hinc)
      rw [l2Cur_edges]
      simp only [List.cons_append, List.mem_cons] at ht
      rcases ht with rfl | ht
      · exact .inl hte
      · exact .inr ⟨t, ht, hte⟩
    · intro u hu e he hue hme
      have hme' : (l2Cur c (mid ++ [py, vy])).edges st.g st.items e := hme
      rw [l2Cur_edges] at hme'
      rcases hme' with h | ⟨t, ht, h⟩
      · exact hsplit _ _ hR₂ c (by simp) u hu e he h hue
      · exact hsplit _ _ hR₂ t (List.mem_cons_of_mem _ ht) u hu e he h hue

/-- The type-1 vertex close: `py` is unwrapped, then `py` and `vy` are merged into `c`, which is
retargeted to `curV`; the closed entry sits at depth `l` and touches the path only at `curV`. -/
theorem closeVertOk_type1 {curV d l : Nat} {o : DfsOut} {base : List TEntry} {s st : WalkState}
    {c : TEntry} {mid : List TEntry} {py vy : TEntry} (isSingle : Bool)
    (hC : EarClose curV d l o true base s st c mid py vy) (hst : Shape st) (ht1 : o.cls.isType1 = true)
    (hl : l < d) (hpath : ∀ k k', k < k' → k' ≤ d → s.stackVerts[k]! ≠ s.stackVerts[k']!)
    (hchild : ∀ k, k ≤ d → s.stackVerts[k]! ≠ s.stackVerts[d + 1]!)
    (hdir : s.stackDir[d]! = !s.stackDir[l]!) :
    CloseVertOk (d + 1) curV s.stackDir[d]! true base.length isSingle st := by
  obtain ⟨hmid, htouch⟩ := hC.type1 ht1
  subst hmid
  have hts : st.tstack = c :: py :: vy :: base := hC.tstack
  obtain ⟨i, hpysp, hroot⟩ := hC.py_item
  have hdl : st.stackDir[py.topDepth]! = s.stackDir[l]! := by rw [hC.py_top]; exact hC.dir_l
  have hside : getSide py.spans (!st.stackDir[py.topDepth]!) = [] := by
    rw [hdl, hpysp]; exact getSide_setSides_not _ _
  have hsingle : getSide py.spans st.stackDir[py.topDepth]! = [i] := by
    rw [hdl, hpysp]; cases s.stackDir[l]! <;> rfl
  have hib : i ∈ py.spans.1 ++ py.spans.2 := by
    rw [mem_of_getSide_nil _ py.spans hside, hsingle]; exact List.mem_singleton_self i
  have htyn : (if isSingle then NodeType.S else NodeType.R) ∉ [NodeType.F, .V, .Q] := by
    cases isSingle <;> decide
  have hdisj : st.tstack.Pairwise fun t t' => ∀ e, e < st.g.ne → t.edges st.g st.items e → ¬ t'.edges st.g st.items e := by
    rw [hC.g]; exact hC.disj
  have hsplit : ∀ l₁ l₂, st.tstack = l₁ ++ l₂ → ∀ u ∈ l₁, ∀ t ∈ l₂, ∀ e, e < st.g.ne →
      u.edges st.g st.items e → ¬ t.edges st.g st.items e :=
    fun l₁ l₂ hl u hu t ht => List.pairwise_split hdisj hl u hu t ht
  have hR₁ : st.tstack = [c, py] ++ (vy :: base) := hts
  have hR₂ : st.tstack = [c, py, vy] ++ base := hts
  have hs₂ : cvS₂ true base.length isSingle st = after (maybeUnwrapNxt (if isSingle then .S else .R)) st := by
    simp only [cvS₂, cvS₁, cvB₁, vertUnwrap, vertPre, Bool.not_true, Bool.false_eq_true, ↓reduceIte, after, result,
      WalkM.map_run, WalkM.pure_run]
  obtain ⟨py', hts₂, hg₂, hsv₂, hsd₂, hpv', hpt', hside', hpE', hallE⟩ :=
    maybeUnwrapNxt_edges (s := st) (ty := if isSingle then .S else .R) hst htyn hts (i := i) hside hsingle
  set s₂ := after (maybeUnwrapNxt (if isSingle then NodeType.S else NodeType.R)) st with hs₂def
  have hs₃ : cvS₃ true base.length isSingle st = { s₂ with tstack := TEntry.mergeInto c py' :: vy :: base } := by
    show (mergeTstackTops.run (cvS₂ true base.length isSingle st)).2 = _
    rw [hs₂, mergeTstackTops_run_eq _ _ _ _ hts₂]
  set m := TEntry.mergeInto (TEntry.mergeInto c py') vy with hm
  have hs₄ : cvS₄ true base.length isSingle st = { s₂ with tstack := m :: base } := by
    show (mergeTstackTops.run (cvS₃ true base.length isSingle st)).2 = _
    rw [hs₃, mergeTstackTops_run_eq _ _ _ _ rfl]
  have hmE : ∀ e, e < st.g.ne → (m.edges st.g s₂.items e ↔
      c.edges st.g st.items e ∨ py.edges st.g st.items e ∨ vy.edges st.g st.items e) := fun e he => by
    rw [hm, TEntry.edges_mergeInto, TEntry.edges_mergeInto, hallE c, hpE' e he, hallE vy, or_assoc]
  have hmtop : m.topDepth = l := by
    show min vy.topDepth (min py'.topDepth c.topDepth) = l
    rw [hpt', hC.py_top]; have := hC.c_top; have := hC.vy_top; omega
  have hmv : m.vStart = vy.vStart := rfl
  have hu : UnwrapOk (if isSingle then .S else .R) st := by
    refine ⟨by rw [hts]; simp, fun _ _ => ?_⟩
    have hn : nxtE st = py := by rw [nxtE, hts]; rfl
    have hd : nxtDir st = st.stackDir[py.topDepth]! := by rw [nxtDir, hn]
    have hh : nxtHead st = i := by rw [nxtHead, hn, hd, hsingle]; rfl
    refine ⟨?_, ?_, ?_, ?_⟩
    · rw [hn, hd]; exact hside
    · rw [hn, hd, hh]; exact hsingle
    · rw [hh]; exact hroot
    · rw [hh]; intro t ht hmem
      have hcur : curE st = c := by rw [curE, hts]; rfl
      rw [hcur, hts] at ht
      simp only [List.tail_cons, List.mem_cons] at ht
      have hpws := hC.span_disj; rw [hts] at hpws
      rcases ht with rfl | ht
      · exact (List.pairwise_cons.1 hpws).1 py (by simp) i hmem hib
      · exact (List.pairwise_cons.1 (List.pairwise_cons.1 hpws).2).1 t (List.mem_cons.2 ht) i hib hmem
  refine ⟨fun h => Bool.noConfusion h, fun _ => hu, ?_, ?_, ?_, fun _ => ?_⟩
  · rw [hs₂]
    intro cur nxt rest h
    rw [hts₂] at h
    simp only [List.cons.injEq] at h
    obtain ⟨rfl, rfl, rfl⟩ := h
    refine ⟨MergeOk.congr_entry (s := st) hg₂ hsv₂ rfl hpv'
      (show min py'.topDepth c.topDepth = min py.topDepth c.topDepth by rw [hpt'])
      (fun e _ => hallE c e) hpE' (hC.fold [] py [vy] rfl), ?_⟩
    intro t ht e he hte hor
    rw [hg₂] at he hte hor
    rw [hallE] at hte
    rcases hor with h | h
    · exact hsplit _ _ hR₁ c (by simp) t ht e he ((hallE c e).1 h) hte
    · exact hsplit _ _ hR₁ py (by simp) t ht e he ((hpE' e he).1 h) hte
  · rw [hs₃]
    intro cur nxt rest h
    simp only [List.cons.injEq] at h
    obtain ⟨rfl, rfl, rfl⟩ := h
    refine ⟨MergeOk.congr_entry (s := st) (cur' := TEntry.mergeInto c py') (cur := TEntry.mergeInto c py) hg₂ hsv₂ hpv' rfl
      (show min vy.topDepth (min py'.topDepth c.topDepth) = min vy.topDepth (min py.topDepth c.topDepth) by rw [hpt'])
      (fun e he => by
        show (TEntry.mergeInto c py').edges st.g s₂.items e ↔ (TEntry.mergeInto c py).edges st.g st.items e
        rw [TEntry.edges_mergeInto, TEntry.edges_mergeInto, hallE c, hpE' e he])
      (fun e _ => hallE vy e) (hC.fold [py] vy [] rfl), ?_⟩
    intro t ht e he hte hor
    simp only [hg₂] at he hte hor
    simp only [hallE] at hte hor
    rw [TEntry.edges_mergeInto] at hor
    have hpy : py'.edges st.g st.items e ↔ py.edges st.g st.items e := (hallE py' e).symm.trans (hpE' e he)
    rcases hor with (h | h) | h
    · exact hsplit _ _ hR₂ c (by simp) t ht e he h hte
    · exact hsplit _ _ hR₂ py (by simp) t ht e he (hpy.1 h) hte
    · exact hsplit _ _ hR₂ vy (by simp) t ht e he h hte
  · rw [hs₄]
    refine ⟨by simp, ?_, ?_⟩
    · refine .inr (.inl ?_)
      show s₂.g.Interior (m.edges s₂.g s₂.items) m.vStart
      rw [hg₂, hmv]
      intro e he hinc
      rw [hmE e he]
      rw [hC.g] at he hinc
      obtain ⟨t, ht, hte⟩ := hC.sub_cover e he (hC.y_edges e he hinc)
      rw [← hC.g] at hte
      simp only [List.cons_append, List.nil_append, List.mem_cons, List.not_mem_nil, or_false] at ht
      rcases ht with rfl | rfl | rfl
      · exact .inl hte
      · exact .inr (.inl hte)
      · exact .inr (.inr hte)
    · intro u hu e he hue hme
      simp only [hg₂] at he hue hme
      simp only [hallE] at hue
      rcases (hmE e he).1 hme with h | h | h
      · exact hsplit _ _ hR₂ c (by simp) u hu e he h hue
      · exact hsplit _ _ hR₂ py (by simp) u hu e he h hue
      · exact hsplit _ _ hR₂ vy (by simp) u hu e he h hue
  · have hs₅ : cvS₅ curV s.stackDir[d]! true base.length isSingle st =
        { s₂ with tstack :=
          { m with vStart := curV, spans := setSides (!s.stackDir[d]!) (m.spans.1 ++ m.spans.2) [] } :: base } := by
      show ((retarget curV s.stackDir[d]!).run (cvS₄ true base.length isSingle st)).2 = _
      rw [hs₄, retarget_run_eq _ _ _ _ _ rfl]
    rw [hs₅]
    set r : TEntry := { m with vStart := curV, spans := setSides (!s.stackDir[d]!) (m.spans.1 ++ m.spans.2) [] }
      with hr
    have hcur : curE { s₂ with tstack := r :: base } = r := rfl
    have hrE : ∀ e, r.edges st.g s₂.items e ↔ m.edges st.g s₂.items e := fun e => by
      simp only [TEntry.edges, hr, mem_setSides]
    have hrtop : r.topDepth = l := hmtop
    refine ⟨by simp, ?_, ?_⟩
    · rw [hcur, hrtop]
      show getSide (setSides (!s.stackDir[d]!) (m.spans.1 ++ m.spans.2) []) (!s₂.stackDir[l]!) = []
      rw [hsd₂, hC.dir_l, hdir, Bool.not_not]
      exact getSide_setSides_not _ _
    · intro k hk hkD
      rw [hcur, hrtop] at hk
      rw [hcur]
      show s₂.stackVerts[k]! = curV ∨ s₂.g.Interior (r.edges s₂.g s₂.items) s₂.stackVerts[k]! ∨
        ¬ s₂.g.Touches (r.edges s₂.g s₂.items) s₂.stackVerts[k]!
      rw [hsv₂, hg₂, hC.sv]
      have hne : s.stackVerts[k]! ≠ s.stackVerts[l]! := by
        rcases Nat.lt_or_ge k (d + 1) with h | h
        · exact fun h' => hpath l k hk (by omega) h'.symm
        · have : k = d + 1 := by omega
          subst this; exact fun h' => hchild l (by omega) h'.symm
      have hEq : ∀ e, e < st.g.ne → (r.edges st.g s₂.items e ↔
          c.edges st.g st.items e ∨ py.edges st.g st.items e ∨ vy.edges st.g st.items e) :=
        fun e he => (hrE e).trans (hmE e he)
      by_cases htc : st.g.Touches (r.edges st.g s₂.items) s.stackVerts[k]!
      · have htc' : s.g.Touches (fun e => c.edges s.g st.items e ∨ py.edges s.g st.items e ∨
            vy.edges s.g st.items e) s.stackVerts[k]! := by
          rw [← hC.g]; exact (Graph.Touches.congr hEq).1 htc
        rcases htouch _ htc' with h | h | h
        · exact .inl h
        · exact absurd h hne
        · refine .inr (.inl ?_)
          rw [← hC.g] at h
          exact (Graph.Interior.congr hEq).2 h
      · exact .inr (.inr htc)

/-- `CloseVertOk` at `feS₂` from the close contract. -/
theorem closeVertOk_of_close {curV d l : Nat} {o : DfsOut} {base : List TEntry} {s st : WalkState}
    {c : TEntry} {mid : List TEntry} {py vy : TEntry} (isSingle : Bool)
    (hC : EarClose curV d l o true base s st c mid py vy) (hst : Shape st)
    (hl : l < d) (hpath : ∀ k k', k < k' → k' ≤ d → s.stackVerts[k]! ≠ s.stackVerts[k']!)
    (hchild : ∀ k, k ≤ d → s.stackVerts[k]! ≠ s.stackVerts[d + 1]!)
    (hdir : s.stackDir[d]! = !s.stackDir[l]!) :
    CloseVertOk (d + 1) curV s.stackDir[d]! o.cls.isType1 base.length isSingle st := by
  cases h : o.cls.isType1
  · exact closeVertOk_type2 isSingle _ hC
  · exact closeVertOk_type1 isSingle hC hst h hl hpath hchild hdir

end Spqr
