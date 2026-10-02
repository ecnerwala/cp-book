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

end Spqr
