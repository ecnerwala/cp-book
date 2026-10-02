import Spqr.StEars

/-!
# Merge loops keep the st-reading

Loops 2 and 3 of `finishEdge` only run `mergeTstackTops`. A merge of the top two entries keeps the
left and right readings; as long as the loop ends with at least one entry above the suffix `B`, every
merge happened above `B`, so `L1StInv` is kept.
-/

namespace Spqr
open WalkM WalkState

def mergeTopsN : Nat → List TEntry → List TEntry
  | 0, l => l
  | k + 1, l => mergeTopsN k (mergeTops l)

theorem iter_mergeTstackTops (k : Nat) (s : WalkState) :
    iter mergeTstackTops k s = { s with tstack := mergeTopsN k s.tstack } := by
  induction k generalizing s with
  | zero => rfl
  | succ k ih =>
    show iter mergeTstackTops k (mergeTstackTops.run s).2 = _
    rw [ih, run_mergeTstackTops]
    rfl

theorem mergeTopsN_nil (k : Nat) : mergeTopsN k [] = [] := by
  induction k with
  | zero => rfl
  | succ k ih => exact ih

theorem mergeTopsN_length (k : Nat) (l : List TEntry) (h : 1 ≤ (mergeTopsN k l).length) :
    (mergeTopsN k l).length + k = l.length := by
  induction k generalizing l with
  | zero => rfl
  | succ k ih =>
    match l with
    | [] => rw [mergeTopsN_nil] at h; simp at h
    | [a] => rw [show mergeTopsN (k + 1) [a] = mergeTopsN k [] from rfl, mergeTopsN_nil] at h; simp at h
    | b :: a :: rest =>
      have h1 := ih (mergeTops (b :: a :: rest)) h
      have h2 : (mergeTops (b :: a :: rest)).length = rest.length + 1 := by simp [mergeTops]
      show (mergeTopsN k (mergeTops (b :: a :: rest))).length + (k + 1) = rest.length + 2
      omega

theorem mergeTopsN_above (k : Nat) (new B : List TEntry) (hk : k < new.length) :
    ∃ new', mergeTopsN k (new ++ B) = new' ++ B ∧ readL new' = readL new ∧ readR new' = readR new ∧
      new'.length + k = new.length := by
  induction k generalizing new with
  | zero => exact ⟨new, rfl, rfl, rfl, rfl⟩
  | succ k ih =>
    match new, hk with
    | [], hk => simp at hk
    | [a], hk => simp at hk
    | b :: a :: rest, hk =>
      have hm : mergeTops (b :: a :: rest ++ B) = mergeTops (b :: a :: rest) ++ B := rfl
      have hl : (mergeTops (b :: a :: rest)).length = rest.length + 1 := by simp [mergeTops]
      obtain ⟨new', h1, h2, h3, h4⟩ := ih (mergeTops (b :: a :: rest)) (by simp at hk; omega)
      refine ⟨new', ?_, ?_, ?_, ?_⟩
      · show mergeTopsN k (mergeTops (b :: a :: rest ++ B)) = new' ++ B
        rw [hm, h1]
      · rw [h2]; simp [mergeTops, readL]
      · rw [h3]; simp [mergeTops, readR]
      · simp; omega

theorem readStack_append_congr {new new' B : List TEntry} (hL : readL new' = readL new)
    (hR : readR new' = readR new) : readStack (new' ++ B) = readStack (new ++ B) := by
  simp only [readStack, readL_append, readR_append, hL, hR]

theorem L1StInv.mergeTopsN {g : Graph} {s : WalkState} {d : Nat} {ps : List StPiece}
    {blocks : List StBlock} {B : List TEntry} {st : WalkState} (hJ : L1StInv g s d ps blocks B st)
    (k : Nat) (hk : k + B.length < st.tstack.length) :
    L1StInv g s d ps blocks B (iter mergeTstackTops k st) := by
  obtain ⟨new, hts, hR⟩ := hJ.read
  have hlen : new.length + B.length = st.tstack.length := by rw [hts, List.length_append]
  obtain ⟨new', h1, h2, h3, -⟩ := mergeTopsN_above k new B (by omega)
  rw [iter_mergeTstackTops]
  refine ⟨⟨new', by dsimp only; rw [hts, h1], ?_⟩, ?_, hJ.dirs⟩
  · unfold StRead at hR ⊢
    dsimp only
    rw [readStack, h2, h3]
    exact hR
  · exact StItems.congr hJ.items (by dsimp only; rw [hts, h1]; exact readStack_append_congr h2 h3) rfl

/-- A merge loop that ends with an entry above `B` keeps the reading. -/
theorem L1StInv.mergeLoop {g : Graph} {s : WalkState} {d : Nat} {ps : List StPiece}
    {blocks : List StBlock} {B : List TEntry} (cond : WalkM Bool) (hcond : ∀ s, (cond.run s).2 = s)
    (fuel : Nat) {st : WalkState} (hJ : L1StInv g s d ps blocks B st)
    (hlen : B.length + 1 ≤ ((loop fuel cond mergeTstackTops).run st).2.tstack.length) :
    L1StInv g s d ps blocks B ((loop fuel cond mergeTstackTops).run st).2 := by
  obtain ⟨k, hk, -⟩ := loop_run_iter cond mergeTstackTops hcond fuel st
  rw [hk] at hlen ⊢
  rw [iter_mergeTstackTops] at hlen
  dsimp only at hlen
  have := mergeTopsN_length k st.tstack (by omega)
  exact hJ.mergeTopsN k (by omega)

theorem L1StInv.mergeLate {g : Graph} {s : WalkState} {d : Nat} {ps : List StPiece}
    {blocks : List StBlock} {B : List TEntry} (d' : Nat) {st : WalkState}
    (hJ : L1StInv g s d ps blocks B st)
    (hlen : B.length + 1 ≤ ((mergeLate d').run st).2.tstack.length) :
    L1StInv g s d ps blocks B ((mergeLate d').run st).2 := by
  rw [mergeLate_run] at hlen ⊢
  by_cases h : (curE st).firstIdx > st.firstOccurrence[d']!
  · simp only [h, ↓reduceIte] at hlen ⊢
    exact hJ.mergeLoop _ (fun _ => rfl) _ hlen
  · simp only [h, ↓reduceIte] at hlen ⊢
    exact hJ

end Spqr
