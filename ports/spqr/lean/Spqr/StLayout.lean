import Mathlib.Algebra.BigOperators.Group.List.Basic
import Mathlib.Data.List.Nodup
import Spqr.Relabel

/-!
# The R-node layout as a fold

`layoutNode .R` written with explicit folds (`LayoutR.run`), and the bracket property of its
adjacency rows (`layoutNode_r_bracket`).
-/

namespace Spqr

namespace LayoutR

/-- `adjBounds[i - 2 nvSt] += 1`. -/
def inc (nvSt : Nat) (l : Layout) (i : Nat) : Layout :=
  { l with adjBounds := l.adjBounds.modify (i - 2 * nvSt) (· + 1) }

/-- Count the two half-row entries of the edge `p`. -/
def countStep (nvSt : Nat) (l : Layout) (p : Nat × Nat) : Layout :=
  inc nvSt (inc nvSt l (2 * p.1 + 2)) (2 * p.2 + 1)

/-- Prefix-sum step: replace the count at bound `i` by the running offset. -/
def prefixStep (nvSt : Nat) (s : Layout × Nat) (i : Nat) : Layout × Nat :=
  ({ s.1 with adjBounds := s.1.adjBounds.set! (i - 2 * nvSt) s.2 },
    s.2 + s.1.adjBounds[i - 2 * nvSt]!)

/-- Reverse-fill step: place the edge `p` (node-edge `s.2 - 1`) at the current cursors. -/
def fillStep (nvSt neSt node : Nat) (s : Layout × Nat) (p : Nat × Nat) : Layout × Nat :=
  ((countStep nvSt s.1 p).setNe neSt node (s.2 - 1) (p.1, p.2)
      (s.1.adjBounds[2 * p.1 + 2 - 2 * nvSt]!, s.1.adjBounds[2 * p.2 + 1 - 2 * nvSt]!),
    s.2 - 1)

def run (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat)) : Layout :=
  let l₁ := E.foldl (countStep nvSt)
    (inc nvSt (inc nvSt (Layout.empty (nvEn - nvSt) (neEn - neSt)) (2 * nvSt + 2)) (2 * nvEn - 1))
  let l₂ := ((List.range' (2 * nvSt + 1) (2 * nvEn + 1 - (2 * nvSt + 1))).foldl (prefixStep nvSt)
    (l₁, 2 * neSt)).1
  let l₃ := (E.reverse.foldl (fillStep nvSt neSt node) (inc nvSt l₂ (2 * nvSt + 2), neEn)).1
  inc nvSt (l₃.setNe neSt node neSt (nvSt, nvEn - 1) (2 * neSt, 2 * neEn - 1)) (2 * nvEn - 1)

theorem forIn_id_yield {α β : Type} (l : List α) (f : α → β → β) (init : β) :
    (forIn (m := Id) l init fun a b => ForInStep.yield (f a b)) = l.foldl (fun b a => f a b) init := by
  induction l generalizing init with
  | nil => rfl
  | cons a l ih => simp [List.forIn_cons, ih]; rfl

theorem layoutNode_R_eq (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat))
    (hv : nvSt + 2 ≤ nvEn) :
    layoutNode .R node nvSt nvEn neSt neEn E = run node nvSt nvEn neSt neEn E := by
  simp only [layoutNode, Id.run, Std.Legacy.Range.forIn_eq_forIn_range', Std.Legacy.Range.size,
    bind, pure, Id.instMonad, beq_self_eq_true, Bool.or_self, BEq.beq, reduceCtorEq,
    forIn_id_yield]
  have : ¬ (nvEn - nvSt = 1) := by omega
  simp only [this, decide_false, Bool.false_eq_true, ↓reduceIte, Bool.or_self, Nat.add_sub_cancel,
    Nat.div_one]
  rfl


/-! ### Counting -/

/-- The edge `p` has an entry at adjacency bound `i`: `i = 2 p.1 + 2` (its entry in the higher row
of `p.1`) or `i = 2 p.2 + 1` (its entry in the lower row of `p.2`). -/
def hits (i : Nat) (p : Nat × Nat) : Bool := 2 * p.1 + 2 == i || 2 * p.2 + 1 == i

/-- Number of entries of the edges `P` at bound `i`. -/
def cnt (i : Nat) (P : List (Nat × Nat)) : Nat := (P.filter (hits i)).length

/-- The destination of the entry of `p` at bound `i`. -/
def other (i : Nat) (p : Nat × Nat) : Nat := if i % 2 = 0 then p.2 else p.1

/-- Global slot of the first entry at bound `i`: `2 neSt` plus the entries at bounds below `i`. -/
def start (nvSt neSt : Nat) (A : List (Nat × Nat)) (i : Nat) : Nat :=
  2 * neSt + ((List.range' (2 * nvSt + 1) (i - (2 * nvSt + 1))).map fun j => cnt j A).sum

theorem cnt_nil (i : Nat) : cnt i [] = 0 := rfl

theorem cnt_cons (i : Nat) (p : Nat × Nat) (P : List (Nat × Nat)) :
    cnt i (p :: P) = (if hits i p then 1 else 0) + cnt i P := by
  unfold cnt; rw [List.filter_cons]; split <;> simp <;> omega

theorem cnt_append (i : Nat) (P Q : List (Nat × Nat)) : cnt i (P ++ Q) = cnt i P + cnt i Q := by
  unfold cnt; simp [List.filter_append]

theorem cnt_reverse (i : Nat) (P : List (Nat × Nat)) : cnt i P.reverse = cnt i P := by
  unfold cnt; simp [List.filter_reverse]

theorem hits_iff (i : Nat) (p : Nat × Nat) : hits i p = true ↔ 2 * p.1 + 2 = i ∨ 2 * p.2 + 1 = i := by
  simp [hits]

theorem start_succ (nvSt neSt : Nat) (A : List (Nat × Nat)) (i : Nat) (hi : 2 * nvSt + 1 ≤ i) :
    start nvSt neSt A (i + 1) = start nvSt neSt A i + cnt i A := by
  unfold start
  rw [show i + 1 - (2 * nvSt + 1) = (i - (2 * nvSt + 1)) + 1 by omega, List.range'_1_concat,
    show 2 * nvSt + 1 + (i - (2 * nvSt + 1)) = i by omega]
  simp [List.sum_append]; omega

theorem start_le (nvSt neSt : Nat) (A : List (Nat × Nat)) {i j : Nat} (hij : i ≤ j) :
    start nvSt neSt A i ≤ start nvSt neSt A j := by
  induction j with
  | zero => cases Nat.le_zero.1 hij; exact Nat.le_refl _
  | succ j ih =>
    rcases Nat.lt_or_eq_of_le hij with h | rfl
    · refine Nat.le_trans (ih (by omega)) ?_
      by_cases hj : 2 * nvSt + 1 ≤ j
      · rw [start_succ _ _ _ _ hj]; omega
      · unfold start
        rw [show j - (2 * nvSt + 1) = 0 by omega, show j + 1 - (2 * nvSt + 1) = 0 by omega]
        exact Nat.le_refl _
    · exact Nat.le_refl _

theorem start_low (nvSt neSt : Nat) (A : List (Nat × Nat)) (i : Nat) (hi : i ≤ 2 * nvSt + 1) :
    start nvSt neSt A i = 2 * neSt := by
  unfold start; rw [show i - (2 * nvSt + 1) = 0 by omega]; rfl

theorem sum_map_ite_hits (p : Nat × Nat) (l : List Nat) :
    (l.map fun j => if hits j p then 1 else 0).sum = l.countP fun j => hits j p := by
  induction l with
  | nil => rfl
  | cons j l ih =>
    rw [List.map_cons, List.sum_cons, ih, List.countP_cons]
    by_cases h : hits j p = true <;> simp [h] <;> omega

theorem countP_hits_range' (p : Nat × Nat) (s n : Nat) (h1 : s ≤ 2 * p.1 + 2) (h2 : 2 * p.1 + 2 < s + n)
    (h3 : s ≤ 2 * p.2 + 1) (h4 : 2 * p.2 + 1 < s + n) (hp : p.1 < p.2) :
    ((List.range' s n).countP fun j => hits j p) = 2 := by
  have key : ∀ (l : List Nat), (l.countP fun j => hits j p) = l.count (2 * p.1 + 2) + l.count (2 * p.2 + 1) := by
    intro l
    induction l with
    | nil => rfl
    | cons j l ih =>
      rw [List.countP_cons, ih, List.count_cons, List.count_cons]
      by_cases ha : j = 2 * p.1 + 2 <;> by_cases hb : j = 2 * p.2 + 1 <;> simp [hits, ha, hb, show ¬ (2 * p.1 + 1 = 2 * p.2) by omega,
          show ¬ (2 * p.2 = 2 * p.1 + 1) by omega] <;> omega
  rw [key, List.count_eq_one_of_mem (List.nodup_range' _ _) (by simp; omega),
    List.count_eq_one_of_mem (List.nodup_range' _ _) (by simp; omega)]
  all_goals exact Nat.one_pos

/-- Every edge with `nvSt ≤ p.1 < p.2 < nvEn` has exactly two entries among the bounds
`[2 nvSt + 1, 2 nvEn]`. -/
theorem sum_cnt_range' (nvSt nvEn : Nat) (A : List (Nat × Nat))
    (hA : ∀ p ∈ A, nvSt ≤ p.1 ∧ p.1 < p.2 ∧ p.2 < nvEn) :
    ((List.range' (2 * nvSt + 1) (2 * nvEn - 2 * nvSt)).map fun j => cnt j A).sum = 2 * A.length := by
  induction A with
  | nil => simp [cnt_nil]
  | cons p A ih =>
    simp only [cnt_cons]
    rw [List.sum_map_add, ih (fun q hq => hA q (List.mem_cons_of_mem _ hq)), sum_map_ite_hits]
    obtain ⟨h1, h2, h3⟩ := hA p List.mem_cons_self
    rw [countP_hits_range' p _ _ (by omega) (by omega) (by omega) (by omega) h2]
    simp; omega

theorem start_top (nvSt nvEn neSt : Nat) (A : List (Nat × Nat))
    (hA : ∀ p ∈ A, nvSt ≤ p.1 ∧ p.1 < p.2 ∧ p.2 < nvEn) (hv : nvSt ≤ nvEn) :
    start nvSt neSt A (2 * nvEn + 1) = 2 * neSt + 2 * A.length := by
  unfold start
  rw [show 2 * nvEn + 1 - (2 * nvSt + 1) = 2 * nvEn - 2 * nvSt by omega, sum_cnt_range' _ _ _ hA]


/-! ### Array helpers -/

theorem getElem!_modify_self' {α : Type} [Inhabited α] (xs : Array α) (i : Nat) (f : α → α)
    (hi : i < xs.size) : (xs.modify i f)[i]! = f xs[i]! := by
  simp [Array.getElem!_eq_getD, Array.getD_eq_getD_getElem?, Array.getElem?_modify,
    Array.getElem?_eq_getElem hi]

theorem getElem!_modify_ne' {α : Type} [Inhabited α] (xs : Array α) (i j : Nat) (f : α → α)
    (hij : i ≠ j) : (xs.modify i f)[j]! = xs[j]! := by
  simp [Array.getElem!_eq_getD, Array.getD_eq_getD_getElem?, Array.getElem?_modify, hij]

/-! ### The counting pass -/

theorem countStep_get (nvSt : Nat) (l : Layout) (p : Nat × Nat) (hp : nvSt ≤ p.1 ∧ p.1 < p.2) (k : Nat)
    (hk1 : 2 * p.1 + 2 - 2 * nvSt < l.adjBounds.size) (hk2 : 2 * p.2 + 1 - 2 * nvSt < l.adjBounds.size) :
    (countStep nvSt l p).adjBounds[k]! =
      l.adjBounds[k]! + if hits (k + 2 * nvSt) p then 1 else 0 := by
  unfold countStep inc
  simp only
  by_cases h2 : k = 2 * p.2 + 1 - 2 * nvSt
  · subst h2
    rw [getElem!_modify_self' _ _ _ (by simpa using hk2),
      getElem!_modify_ne' _ _ _ _ (by omega)]
    simp [hits]; omega
  · rw [getElem!_modify_ne' _ _ _ _ (Ne.symm h2)]
    by_cases h1 : k = 2 * p.1 + 2 - 2 * nvSt
    · subst h1
      rw [getElem!_modify_self' _ _ _ hk1]
      simp [hits]; omega
    · rw [getElem!_modify_ne' _ _ _ _ (Ne.symm h1)]
      have : hits (k + 2 * nvSt) p = false := by
        simp [hits]; omega
      simp [this]

theorem countStep_edges (nvSt : Nat) (l : Layout) (p : Nat × Nat) :
    (countStep nvSt l p).edges = l.edges := rfl
theorem countStep_adjDat (nvSt : Nat) (l : Layout) (p : Nat × Nat) :
    (countStep nvSt l p).adjDat = l.adjDat := rfl
theorem countStep_size (nvSt : Nat) (l : Layout) (p : Nat × Nat) :
    (countStep nvSt l p).adjBounds.size = l.adjBounds.size := by
  simp [countStep, inc]

theorem foldl_countStep (nvSt nvEn : Nat) (P : List (Nat × Nat))
    (hP : ∀ p ∈ P, nvSt ≤ p.1 ∧ p.1 < p.2 ∧ p.2 < nvEn) :
    ∀ (l : Layout), l.adjBounds.size = 2 * (nvEn - nvSt) + 1 →
      (P.foldl (countStep nvSt) l).edges = l.edges ∧
      (P.foldl (countStep nvSt) l).adjDat = l.adjDat ∧
      (P.foldl (countStep nvSt) l).adjBounds.size = l.adjBounds.size ∧
      ∀ k, (P.foldl (countStep nvSt) l).adjBounds[k]! = l.adjBounds[k]! + cnt (k + 2 * nvSt) P := by
  induction P with
  | nil => intro l _; simp [cnt_nil]
  | cons p P ih =>
    intro l hl
    obtain ⟨h1, h2, h3⟩ := hP p List.mem_cons_self
    have hl' : (countStep nvSt l p).adjBounds.size = 2 * (nvEn - nvSt) + 1 := by
      rw [countStep_size]; exact hl
    obtain ⟨e, d, sz, g⟩ := ih (fun q hq => hP q (List.mem_cons_of_mem _ hq)) _ hl'
    rw [List.foldl_cons]
    refine ⟨e.trans (countStep_edges _ _ _), d.trans (countStep_adjDat _ _ _),
      sz.trans (countStep_size _ _ _), fun k => ?_⟩
    rw [g k, countStep_get nvSt l p ⟨h1, h2⟩ k (by omega) (by omega), cnt_cons]
    omega

/-! ### The prefix-sum pass -/

theorem foldl_prefixStep (nvSt : Nat) :
    ∀ (m s : Nat) (l : Layout) (off : Nat), 2 * nvSt + 1 ≤ s →
      s + m ≤ l.adjBounds.size + 2 * nvSt →
      ((List.range' s m).foldl (prefixStep nvSt) (l, off)).1.edges = l.edges ∧
      ((List.range' s m).foldl (prefixStep nvSt) (l, off)).1.adjDat = l.adjDat ∧
      ((List.range' s m).foldl (prefixStep nvSt) (l, off)).1.adjBounds.size = l.adjBounds.size ∧
      ((List.range' s m).foldl (prefixStep nvSt) (l, off)).2 =
        off + ((List.range' s m).map fun j => l.adjBounds[j - 2 * nvSt]!).sum ∧
      ∀ k, (k + 2 * nvSt < s ∨ s + m ≤ k + 2 * nvSt →
          ((List.range' s m).foldl (prefixStep nvSt) (l, off)).1.adjBounds[k]! = l.adjBounds[k]!) ∧
        (s ≤ k + 2 * nvSt → k + 2 * nvSt < s + m →
          ((List.range' s m).foldl (prefixStep nvSt) (l, off)).1.adjBounds[k]! =
            off + ((List.range' s (k + 2 * nvSt - s)).map fun j => l.adjBounds[j - 2 * nvSt]!).sum) := by
  intro m
  induction m with
  | zero =>
    intro s l off _ _
    refine ⟨rfl, rfl, rfl, by simp, fun k => ⟨fun _ => rfl, fun h1 h2 => by omega⟩⟩
  | succ m ih =>
    intro s l off hs hsz
    rw [List.range'_succ, List.foldl_cons]
    set l' : Layout := { l with adjBounds := l.adjBounds.set! (s - 2 * nvSt) off } with hl'
    have hstep : prefixStep nvSt (l, off) s = (l', off + l.adjBounds[s - 2 * nvSt]!) := rfl
    rw [hstep]
    have hsz' : l'.adjBounds.size = l.adjBounds.size := by simp [hl']
    obtain ⟨e, d, sz, o, g⟩ := ih (s + 1) l' (off + l.adjBounds[s - 2 * nvSt]!) (by omega) (by omega)
    have hsame : ∀ j, s + 1 ≤ j → l'.adjBounds[j - 2 * nvSt]! = l.adjBounds[j - 2 * nvSt]! := by
      intro j hj
      simp only [hl']
      exact Array.getElem!_set!_ne _ _ _ _ (by omega)
    have hmap : ∀ n, ((List.range' (s + 1) n).map fun j => l'.adjBounds[j - 2 * nvSt]!) =
        (List.range' (s + 1) n).map fun j => l.adjBounds[j - 2 * nvSt]! := by
      intro n
      apply List.map_congr_left
      intro j hj
      exact hsame j (List.mem_range'_1.1 hj).1
    refine ⟨e, d, sz.trans hsz', ?_, fun k => ⟨fun hk => ?_, fun hk1 hk2 => ?_⟩⟩
    · rw [o, hmap, List.map_cons, List.sum_cons]; omega
    · rw [(g k).1 (by omega)]
      simp only [hl']
      exact Array.getElem!_set!_ne _ _ _ _ (by omega)
    · rcases Nat.lt_or_ge (k + 2 * nvSt) (s + 1) with h | h
      · have hk : k = s - 2 * nvSt := by omega
        rw [(g k).1 (by omega), hk, show s - 2 * nvSt + 2 * nvSt - s = 0 by omega]
        simp only [hl', List.range'_zero, List.map_nil, List.sum_nil, Nat.add_zero]
        exact Array.getElem!_set!_self _ _ _ (by omega)
      · rw [(g k).2 h (by omega), hmap,
          show k + 2 * nvSt - s = (k + 2 * nvSt - (s + 1)) + 1 by omega,
          List.range'_succ, List.map_cons, List.sum_cons]
        omega


/-! ### The fill pass -/

/-- The cap edge followed by the children: all edges of the node. -/
def allE (nvSt nvEn : Nat) (E : List (Nat × Nat)) : List (Nat × Nat) := (nvSt, nvEn - 1) :: E

/-- The cap occupies the first slot at bound `2 nvSt + 2`, before the children's entries. -/
def bump (nvSt i : Nat) : Nat := if i = 2 * nvSt + 2 then 1 else 0

theorem cnt_allE (nvSt nvEn : Nat) (E : List (Nat × Nat)) (hv : nvSt + 2 ≤ nvEn) (i : Nat) :
    cnt i (allE nvSt nvEn E) = bump nvSt i + (if i = 2 * nvEn - 1 then 1 else 0) + cnt i E := by
  unfold allE bump
  rw [cnt_cons]
  simp only [hits, beq_iff_eq, Bool.or_eq_true]
  split_ifs <;> omega

theorem start_ge (nvSt neSt : Nat) (A : List (Nat × Nat)) (i : Nat) :
    2 * neSt ≤ start nvSt neSt A i := Nat.le_add_right _ _

theorem allE_bounds (nvSt nvEn : Nat) (E : List (Nat × Nat)) (hv : nvSt + 2 ≤ nvEn)
    (hE : ∀ q ∈ E, nvSt ≤ q.1 ∧ q.1 < q.2 ∧ q.2 < nvEn) :
    ∀ q ∈ allE nvSt nvEn E, nvSt ≤ q.1 ∧ q.1 < q.2 ∧ q.2 < nvEn := by
  intro q hq
  rcases List.mem_cons.1 hq with rfl | hq
  · simp; omega
  · exact hE q hq

theorem allE_length (nvSt nvEn : Nat) (E : List (Nat × Nat)) :
    (allE nvSt nvEn E).length = E.length + 1 := rfl

/-- Slots of different bounds are different. -/
theorem slot_ne (nvSt neSt : Nat) (A : List (Nat × Nat)) {i j a b : Nat}
    (hi : 2 * nvSt + 1 ≤ i) (hj : 2 * nvSt + 1 ≤ j) (hij : i ≠ j)
    (ha : a < cnt i A) (hb : b < cnt j A) :
    start nvSt neSt A i + a ≠ start nvSt neSt A j + b := by
  rcases Nat.lt_or_gt_of_ne hij with h | h
  · have h1 := start_succ nvSt neSt A i hi
    have h2 := start_le nvSt neSt A (show i + 1 ≤ j by omega)
    omega
  · have h1 := start_succ nvSt neSt A j hj
    have h2 := start_le nvSt neSt A (show j + 1 ≤ i by omega)
    omega

theorem slot_lt_top (nvSt nvEn neSt : Nat) (E : List (Nat × Nat)) (hv : nvSt + 2 ≤ nvEn)
    (hE : ∀ q ∈ E, nvSt ≤ q.1 ∧ q.1 < q.2 ∧ q.2 < nvEn) {i a : Nat}
    (hi : 2 * nvSt + 1 ≤ i) (hi2 : i ≤ 2 * nvEn) (ha : a < cnt i (allE nvSt nvEn E)) :
    start nvSt neSt (allE nvSt nvEn E) i + a < 2 * neSt + 2 * (E.length + 1) := by
  have h1 := start_succ nvSt neSt (allE nvSt nvEn E) i hi
  have h2 := start_le nvSt neSt (allE nvSt nvEn E) (show i + 1 ≤ 2 * nvEn + 1 by omega)
  have h3 := start_top nvSt nvEn neSt (allE nvSt nvEn E) (allE_bounds nvSt nvEn E hv hE) (by omega)
  rw [allE_length] at h3
  omega

theorem cnt_split (i : Nat) (P R E : List (Nat × Nat)) (hPR : P ++ R = E.reverse) :
    cnt i P + cnt i R = cnt i E := by
  rw [← cnt_append, hPR, cnt_reverse]

/-- Invariant of the reverse fill: the children `P` have been placed, in order, right after the
cap's slot (at bound `2 nvSt + 2`) resp. at the beginning of their bounds' slot ranges. -/
structure FillInv (nvSt nvEn neSt : Nat) (E P : List (Nat × Nat)) (l : Layout) : Prop where
  size_b : l.adjBounds.size = 2 * (nvEn - nvSt) + 1
  size_d : l.adjDat.size = 2 * (E.length + 1)
  bounds : ∀ i, 2 * nvSt + 1 ≤ i → i ≤ 2 * nvEn →
    l.adjBounds[i - 2 * nvSt]! = start nvSt neSt (allE nvSt nvEn E) i + bump nvSt i + cnt i P
  dat : ∀ i, 2 * nvSt + 1 ≤ i → i ≤ 2 * nvEn → ∀ k q, (P.filter (hits i))[k]? = some q →
    (l.adjDat[start nvSt neSt (allE nvSt nvEn E) i + bump nvSt i + k - 2 * neSt]!).destNv =
      other i q

theorem other_even (i : Nat) (p : Nat × Nat) (h : i % 2 = 0) : other i p = p.2 := by
  simp [other, h]
theorem other_odd (i : Nat) (p : Nat × Nat) (h : i % 2 = 1) : other i p = p.1 := by
  simp [other, h]

theorem fillStep_inv (nvSt nvEn neSt node : Nat) (E P R : List (Nat × Nat)) (p : Nat × Nat)
    (l : Layout) (nxt : Nat) (hv : nvSt + 2 ≤ nvEn)
    (hE : ∀ q ∈ E, nvSt ≤ q.1 ∧ q.1 < q.2 ∧ q.2 < nvEn)
    (hPR : P ++ p :: R = E.reverse) (h : FillInv nvSt nvEn neSt E P l) :
    FillInv nvSt nvEn neSt E (P ++ [p]) (fillStep nvSt neSt node (l, nxt) p).1 := by
  have hp : nvSt ≤ p.1 ∧ p.1 < p.2 ∧ p.2 < nvEn := hE p (by
    have : p ∈ E.reverse := by rw [← hPR]; simp
    simpa using this)
  have hcntP : ∀ i, cnt i (P ++ [p]) ≤ cnt i E := by
    intro i
    have := cnt_split i P (p :: R) E hPR
    rw [cnt_cons] at this
    rw [cnt_append, cnt_cons, cnt_nil]
    omega
  have hcntA : ∀ i, bump nvSt i + cnt i E ≤ cnt i (allE nvSt nvEn E) := by
    intro i; rw [cnt_allE _ _ _ hv]; omega
  have hits0 : hits (2 * p.1 + 2) p = true := by simp [hits]
  have hits1 : hits (2 * p.2 + 1) p = true := by simp [hits]
  have hs0lt : bump nvSt (2 * p.1 + 2) + cnt (2 * p.1 + 2) P < cnt (2 * p.1 + 2) (allE nvSt nvEn E) := by
    have h1 := hcntP (2 * p.1 + 2)
    rw [cnt_append, cnt_cons, cnt_nil, hits0] at h1
    have h2 := hcntA (2 * p.1 + 2)
    simp at h1; omega
  have hs1lt : bump nvSt (2 * p.2 + 1) + cnt (2 * p.2 + 1) P < cnt (2 * p.2 + 1) (allE nvSt nvEn E) := by
    have h1 := hcntP (2 * p.2 + 1)
    rw [cnt_append, cnt_cons, cnt_nil, hits1] at h1
    have h2 := hcntA (2 * p.2 + 1)
    simp at h1; omega
  have hb0 := h.bounds (2 * p.1 + 2) (by omega) (by omega)
  have hb1 := h.bounds (2 * p.2 + 1) (by omega) (by omega)
  have hge0 := start_ge nvSt neSt (allE nvSt nvEn E) (2 * p.1 + 2)
  have hge1 := start_ge nvSt neSt (allE nvSt nvEn E) (2 * p.2 + 1)
  have hne01 := slot_ne nvSt neSt (allE nvSt nvEn E) (i := 2 * p.1 + 2) (j := 2 * p.2 + 1)
    (by omega) (by omega) (by omega) hs0lt hs1lt
  have hlt0 := slot_lt_top nvSt nvEn neSt E hv hE (i := 2 * p.1 + 2) (by omega) (by omega) hs0lt
  have hlt1 := slot_lt_top nvSt nvEn neSt E hv hE (i := 2 * p.2 + 1) (by omega) (by omega) hs1lt
  simp only [fillStep, Layout.setNe, countStep_adjDat, hb0, hb1]
  refine ⟨?_, ?_, ?_, ?_⟩
  · rw [countStep_size]; exact h.size_b
  · simp [Array.size_set!, h.size_d]
  · intro i hi1 hi2
    rw [countStep_get nvSt l p ⟨hp.1, hp.2.1⟩ _ (by rw [h.size_b]; omega) (by rw [h.size_b]; omega),
      show i - 2 * nvSt + 2 * nvSt = i by omega, h.bounds i hi1 hi2, cnt_append, cnt_cons, cnt_nil]
    omega
  · intro i hi1 hi2 k q hq
    have hk := (List.getElem?_eq_some_iff.1 hq).1
    rw [List.filter_append] at hq
    have hgei := start_ge nvSt neSt (allE nvSt nvEn E) i
    by_cases hhit : hits i p = true
    · have hf : [p].filter (hits i) = [p] := by simp [hhit]
      rw [hf] at hq
      have hki : k < cnt i P + 1 := by
        rw [List.filter_append, hf, List.length_append, List.length_singleton] at hk
        exact hk
      by_cases hk' : k < (P.filter (hits i)).length
      · rw [List.getElem?_append_left hk'] at hq
        have hkA : bump nvSt i + k < cnt i (allE nvSt nvEn E) := by
          have h1 := hcntP i; have h2 := hcntA i
          rw [cnt_append] at h1
          change k < cnt i P at hk'
          omega
        rcases (hits_iff i p).1 hhit with hi0 | hi1'
        · subst hi0
          rw [Array.getElem!_set!_ne _ _ _ _ (by
              have := slot_ne nvSt neSt (allE nvSt nvEn E) (i := 2 * p.2 + 1) (j := 2 * p.1 + 2)
                (by omega) (by omega) (by omega) hs1lt hkA
              omega),
            Array.getElem!_set!_ne _ _ _ _ (by change k < cnt (2 * p.1 + 2) P at hk'; omega)]
          exact h.dat _ hi1 hi2 k q hq
        · subst hi1'
          rw [Array.getElem!_set!_ne _ _ _ _ (by change k < cnt (2 * p.2 + 1) P at hk'; omega),
            Array.getElem!_set!_ne _ _ _ _ (by
              have := slot_ne nvSt neSt (allE nvSt nvEn E) (i := 2 * p.1 + 2) (j := 2 * p.2 + 1)
                (by omega) (by omega) (by omega) hs0lt hkA
              omega)]
          exact h.dat _ hi1 hi2 k q hq
      · rw [List.getElem?_append_right (Nat.le_of_not_lt hk'), List.getElem?_singleton] at hq
        have hkeq : k = cnt i P := by
          change ¬ k < cnt i P at hk'; omega
        rw [show k - (P.filter (hits i)).length = 0 from by change k - cnt i P = 0; omega] at hq
        simp only [↓reduceIte, Option.some.injEq] at hq
        subst hq
        rw [hkeq]
        rcases (hits_iff i p).1 hhit with hi0 | hi1'
        · subst hi0
          rw [Array.getElem!_set!_ne _ _ _ _ (by omega),
            Array.getElem!_set!_self _ _ _ (by rw [h.size_d]; omega), other_even _ _ (by omega)]
        · subst hi1'
          rw [Array.getElem!_set!_self _ _ _ (by rw [Array.size_set!, h.size_d]; omega),
            other_odd _ _ (by omega)]
    · have hf : [p].filter (hits i) = [] := by simp [hhit]
      rw [hf, List.append_nil] at hq
      have hk' : k < cnt i P := (List.getElem?_eq_some_iff.1 hq).1
      have hkA : bump nvSt i + k < cnt i (allE nvSt nvEn E) := by
        have h1 := hcntP i; have h2 := hcntA i
        rw [cnt_append] at h1
        omega
      have hni : 2 * p.1 + 2 ≠ i ∧ 2 * p.2 + 1 ≠ i := by
        have := (hits_iff i p).not.1 hhit
        omega
      have hne0 := slot_ne nvSt neSt (allE nvSt nvEn E) (i := 2 * p.1 + 2) (j := i)
        (by omega) hi1 hni.1 hs0lt hkA
      have hne1 := slot_ne nvSt neSt (allE nvSt nvEn E) (i := 2 * p.2 + 1) (j := i)
        (by omega) hi1 hni.2 hs1lt hkA
      rw [Array.getElem!_set!_ne _ _ _ _ (by omega), Array.getElem!_set!_ne _ _ _ _ (by omega)]
      exact h.dat i hi1 hi2 k q hq

theorem foldl_fillStep_inv (nvSt nvEn neSt node : Nat) (E : List (Nat × Nat))
    (hv : nvSt + 2 ≤ nvEn) (hE : ∀ q ∈ E, nvSt ≤ q.1 ∧ q.1 < q.2 ∧ q.2 < nvEn) :
    ∀ (P R : List (Nat × Nat)) (l : Layout) (nxt : Nat), P ++ R = E.reverse →
      FillInv nvSt nvEn neSt E P l →
      FillInv nvSt nvEn neSt E (P ++ R) (R.foldl (fillStep nvSt neSt node) (l, nxt)).1 := by
  intro P R
  induction R generalizing P with
  | nil => intro l nxt _ h; simpa using h
  | cons p R ih =>
    intro l nxt hPR h
    rw [List.foldl_cons]
    have := ih (P ++ [p]) (fillStep nvSt neSt node (l, nxt) p).1 (fillStep nvSt neSt node (l, nxt) p).2
      (by simpa using hPR) (fillStep_inv nvSt nvEn neSt node E P R p l nxt hv hE hPR h)
    simpa using this


/-! ### The finished layout -/

theorem getElem!_replicate' {α : Type} [Inhabited α] (n : Nat) (a : α) (k : Nat) :
    (Array.replicate n a)[k]! = if k < n then a else default := by
  simp [Array.getElem!_eq_getD, Array.getD_eq_getD_getElem?, Array.getElem?_replicate]
  split <;> rfl

theorem cnt_low_zero (nvSt nvEn : Nat) (E : List (Nat × Nat)) (hv : nvSt + 2 ≤ nvEn)
    (hE : ∀ q ∈ E, nvSt ≤ q.1 ∧ q.1 < q.2 ∧ q.2 < nvEn) :
    cnt (2 * nvSt + 1) (allE nvSt nvEn E) = 0 := by
  unfold cnt
  rw [List.length_eq_zero_iff, List.filter_eq_nil_iff]
  intro q hq
  have := allE_bounds nvSt nvEn E hv hE q hq
  simp [hits]; omega

theorem cnt_high_zero (nvSt nvEn : Nat) (E : List (Nat × Nat)) (hv : nvSt + 2 ≤ nvEn)
    (hE : ∀ q ∈ E, nvSt ≤ q.1 ∧ q.1 < q.2 ∧ q.2 < nvEn) :
    cnt (2 * nvEn) (allE nvSt nvEn E) = 0 := by
  unfold cnt
  rw [List.length_eq_zero_iff, List.filter_eq_nil_iff]
  intro q hq
  have := allE_bounds nvSt nvEn E hv hE q hq
  simp [hits]; omega

theorem start_end (nvSt nvEn neSt : Nat) (E : List (Nat × Nat)) (hv : nvSt + 2 ≤ nvEn)
    (hE : ∀ q ∈ E, nvSt ≤ q.1 ∧ q.1 < q.2 ∧ q.2 < nvEn) :
    start nvSt neSt (allE nvSt nvEn E) (2 * nvEn) = 2 * neSt + 2 * (E.length + 1) := by
  have h1 := start_top nvSt nvEn neSt (allE nvSt nvEn E) (allE_bounds nvSt nvEn E hv hE) (by omega)
  have h2 := start_succ nvSt neSt (allE nvSt nvEn E) (2 * nvEn) (by omega)
  rw [allE_length] at h1
  rw [cnt_high_zero nvSt nvEn E hv hE] at h2
  omega

theorem start_two (nvSt nvEn neSt : Nat) (E : List (Nat × Nat)) (hv : nvSt + 2 ≤ nvEn)
    (hE : ∀ q ∈ E, nvSt ≤ q.1 ∧ q.1 < q.2 ∧ q.2 < nvEn) :
    start nvSt neSt (allE nvSt nvEn E) (2 * nvSt + 2) = 2 * neSt := by
  rw [start_succ _ _ _ _ (Nat.le_refl _), start_low _ _ _ _ (Nat.le_refl _), cnt_low_zero nvSt nvEn E hv hE]
  all_goals omega

theorem slot_pos (nvSt nvEn neSt : Nat) (E : List (Nat × Nat)) (hv : nvSt + 2 ≤ nvEn) {i : Nat}
    (hi : 2 * nvSt + 3 ≤ i) : 2 * neSt + 1 ≤ start nvSt neSt (allE nvSt nvEn E) i := by
  have h1 : start nvSt neSt (allE nvSt nvEn E) (2 * nvSt + 3) =
      start nvSt neSt (allE nvSt nvEn E) (2 * nvSt + 2) + cnt (2 * nvSt + 2) (allE nvSt nvEn E) :=
    start_succ nvSt neSt (allE nvSt nvEn E) (2 * nvSt + 2) (by omega)
  have h2 := start_le nvSt neSt (allE nvSt nvEn E) hi
  have h3 := start_ge nvSt neSt (allE nvSt nvEn E) (2 * nvSt + 2)
  have h4 : 1 ≤ cnt (2 * nvSt + 2) (allE nvSt nvEn E) := by
    rw [cnt_allE _ _ _ hv]; unfold bump; simp only [↓reduceIte]; omega
  omega

theorem slot_lt_top' (nvSt nvEn neSt : Nat) (E : List (Nat × Nat)) (hv : nvSt + 2 ≤ nvEn)
    (hE : ∀ q ∈ E, nvSt ≤ q.1 ∧ q.1 < q.2 ∧ q.2 < nvEn) {i a : Nat}
    (hi : 2 * nvSt + 1 ≤ i) (hi2 : i ≤ 2 * nvEn)
    (ha : a + (if i = 2 * nvEn - 1 then 1 else 0) < cnt i (allE nvSt nvEn E)) :
    start nvSt neSt (allE nvSt nvEn E) i + a + 1 < 2 * neSt + 2 * (E.length + 1) := by
  have hend := start_end nvSt nvEn neSt E hv hE
  have h1 := start_succ nvSt neSt (allE nvSt nvEn E) i hi
  by_cases hil : i = 2 * nvEn - 1
  · subst hil
    rw [if_pos rfl] at ha
    rw [show 2 * nvEn - 1 + 1 = 2 * nvEn by omega] at h1
    omega
  · rw [if_neg hil] at ha
    rcases Nat.lt_or_ge i (2 * nvEn - 1) with h | h
    · have h2 := start_le nvSt neSt (allE nvSt nvEn E) (show i + 1 ≤ 2 * nvEn - 1 by omega)
      have h3 := start_succ nvSt neSt (allE nvSt nvEn E) (2 * nvEn - 1) (by omega)
      have h4 : 1 ≤ cnt (2 * nvEn - 1) (allE nvSt nvEn E) := by
        rw [cnt_allE _ _ _ hv]; simp only [↓reduceIte]; omega
      rw [show 2 * nvEn - 1 + 1 = 2 * nvEn by omega] at h3
      omega
    · have : i = 2 * nvEn := by omega
      subst this
      rw [cnt_high_zero nvSt nvEn E hv hE] at ha
      omega

/-- The destinations of the adjacency entries at bound `i`, in slot order: the cap first (at
bound `2 nvSt + 2`), then the children in reverse order, then the cap (at bound `2 nvEn - 1`). -/
def rowDest (nvSt nvEn : Nat) (E : List (Nat × Nat)) (i : Nat) : List Nat :=
  (if i = 2 * nvSt + 2 then [nvEn - 1] else []) ++ (E.reverse.filter (hits i)).map (other i) ++
    (if i = 2 * nvEn - 1 then [nvSt] else [])

theorem rowDest_length (nvSt nvEn : Nat) (E : List (Nat × Nat)) (hv : nvSt + 2 ≤ nvEn) (i : Nat) :
    (rowDest nvSt nvEn E i).length = cnt i (allE nvSt nvEn E) := by
  rw [cnt_allE _ _ _ hv, ← cnt_reverse]
  unfold rowDest bump cnt
  simp only [List.length_append, List.length_map]
  split_ifs <;> simp <;> omega

theorem run_spec (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat)) (hv : nvSt + 2 ≤ nvEn)
    (hE : ∀ q ∈ E, nvSt ≤ q.1 ∧ q.1 < q.2 ∧ q.2 < nvEn) (hne : neEn = neSt + E.length + 1) :
    (∀ i, 2 * nvSt + 1 ≤ i → i ≤ 2 * nvEn →
      (run node nvSt nvEn neSt neEn E).adjBounds[i - 2 * nvSt]! =
        start nvSt neSt (allE nvSt nvEn E) (i + 1)) ∧
    (∀ i, 2 * nvSt + 1 ≤ i → i ≤ 2 * nvEn → ∀ k q, (rowDest nvSt nvEn E i)[k]? = some q →
      ((run node nvSt nvEn neSt neEn E).adjDat[start nvSt neSt (allE nvSt nvEn E) i + k -
        2 * neSt]!).destNv = q) := by
  have hA' := allE_bounds nvSt nvEn E hv hE
  have hend := start_end nvSt nvEn neSt E hv hE
  have htwo := start_two nvSt nvEn neSt E hv hE
  have hlow := cnt_low_zero nvSt nvEn E hv hE
  simp only [run]
  set l₀ := inc nvSt (inc nvSt (Layout.empty (nvEn - nvSt) (neEn - neSt)) (2 * nvSt + 2))
    (2 * nvEn - 1) with hl₀
  have hl₀' : l₀ = countStep nvSt (Layout.empty (nvEn - nvSt) (neEn - neSt)) (nvSt, nvEn - 1) := by
    simp only [hl₀, countStep]
    rw [show 2 * (nvEn - 1) + 1 = 2 * nvEn - 1 by omega]
  have hsz0 : l₀.adjBounds.size = 2 * (nvEn - nvSt) + 1 := by simp [hl₀, inc, Layout.empty]
  have hget0 : ∀ k, l₀.adjBounds[k]! = cnt (k + 2 * nvSt) [(nvSt, nvEn - 1)] := by
    intro k
    rw [hl₀', countStep_get nvSt _ _ ⟨Nat.le_refl _, by simp; omega⟩ k
      (by simp [Layout.empty]; omega) (by simp [Layout.empty]; omega), cnt_cons, cnt_nil]
    simp only [Layout.empty, getElem!_replicate']
    split <;> simp
  set l₁ := E.foldl (countStep nvSt) l₀ with hl₁
  obtain ⟨e1, d1, sz1, g1⟩ := foldl_countStep nvSt nvEn E hE l₀ hsz0
  have hget1 : ∀ k, l₁.adjBounds[k]! = cnt (k + 2 * nvSt) (allE nvSt nvEn E) := by
    intro k; rw [g1, hget0, allE, cnt_cons, cnt_cons, cnt_nil]; omega
  have hsz1 : l₁.adjBounds.size = 2 * (nvEn - nvSt) + 1 := sz1.trans hsz0
  obtain ⟨e2, d2, sz2, o2, g2⟩ := foldl_prefixStep nvSt (2 * nvEn + 1 - (2 * nvSt + 1)) (2 * nvSt + 1)
    l₁ (2 * neSt) (Nat.le_refl _) (by omega)
  set l₂ := ((List.range' (2 * nvSt + 1) (2 * nvEn + 1 - (2 * nvSt + 1))).foldl (prefixStep nvSt)
    (l₁, 2 * neSt)).1 with hl₂
  have hget2 : ∀ i, 2 * nvSt + 1 ≤ i → i ≤ 2 * nvEn →
      l₂.adjBounds[i - 2 * nvSt]! = start nvSt neSt (allE nvSt nvEn E) i := by
    intro i hi1 hi2
    rw [(g2 (i - 2 * nvSt)).2 (by omega) (by omega), show i - 2 * nvSt + 2 * nvSt = i by omega]
    unfold start
    congr 2
    apply List.map_congr_left
    intro j hj
    have := (List.mem_range'_1.1 hj).1
    rw [hget1, show j - 2 * nvSt + 2 * nvSt = j by omega]
  have hsz2 : l₂.adjBounds.size = 2 * (nvEn - nvSt) + 1 := sz2.trans hsz1
  have hdat2 : l₂.adjDat = Array.replicate (2 * (neEn - neSt)) default := by
    rw [d2, d1]; simp [hl₀, inc, Layout.empty]
  have hinv0 : FillInv nvSt nvEn neSt E [] (inc nvSt l₂ (2 * nvSt + 2)) := by
    refine ⟨by simp [inc, hsz2], by simp [inc, hdat2]; omega, ?_, ?_⟩
    · intro i hi1 hi2
      unfold inc
      simp only
      rw [cnt_nil, Nat.add_zero]
      by_cases hi : i = 2 * nvSt + 2
      · subst hi
        rw [getElem!_modify_self' _ _ _ (by rw [hsz2]; omega), hget2 _ (by omega) (by omega)]
        simp [bump]
      · rw [getElem!_modify_ne' _ _ _ _ (by omega), hget2 _ hi1 hi2]
        simp [bump, hi]
    · intro i _ _ k q hq; simp at hq
  have hinv := foldl_fillStep_inv nvSt nvEn neSt node E hv hE [] E.reverse _ neEn (by simp) hinv0
  rw [List.nil_append] at hinv
  set l₄ := (E.reverse.foldl (fillStep nvSt neSt node) (inc nvSt l₂ (2 * nvSt + 2), neEn)).1 with hl₄
  have hsz4 := hinv.size_d
  have hcntA : ∀ i, cnt i (allE nvSt nvEn E) =
      bump nvSt i + (if i = 2 * nvEn - 1 then 1 else 0) + cnt i E.reverse := by
    intro i; rw [cnt_allE _ _ _ hv, cnt_reverse]
  constructor
  · intro i hi1 hi2
    simp only [inc, Layout.setNe]
    have hb := hinv.bounds i hi1 hi2
    have hs := start_succ nvSt neSt (allE nvSt nvEn E) i hi1
    have hc := hcntA i
    by_cases hi : i = 2 * nvEn - 1
    · subst hi
      rw [getElem!_modify_self' _ _ _ (by rw [hinv.size_b]; omega), hb]
      rw [if_pos rfl] at hc
      omega
    · rw [getElem!_modify_ne' _ _ _ _ (by omega), hb]
      rw [if_neg hi] at hc
      omega
  · intro i hi1 hi2 k q hq
    simp only [inc, Layout.setNe, Nat.sub_self]
    have hge := start_ge nvSt neSt (allE nvSt nvEn E) i
    unfold rowDest at hq
    have hmid : ∀ k' q', ((E.reverse.filter (hits i)).map (other i))[k']? = some q' →
        k' < cnt i E.reverse ∧
        (l₄.adjDat[start nvSt neSt (allE nvSt nvEn E) i + bump nvSt i + k' - 2 * neSt]!).destNv = q' := by
      intro k' q' hq'
      rw [List.getElem?_map, Option.map_eq_some_iff] at hq'
      obtain ⟨p, hp, rfl⟩ := hq'
      exact ⟨(List.getElem?_eq_some_iff.1 hp).1, hinv.dat i hi1 hi2 k' p hp⟩
    by_cases hi' : i = 2 * nvSt + 2
    · subst hi'
      rw [if_pos rfl, if_neg (by omega), List.append_nil, List.singleton_append] at hq
      have hb : bump nvSt (2 * nvSt + 2) = 1 := by simp [bump]
      cases k with
      | zero =>
        rw [List.getElem?_cons_zero, Option.some.injEq] at hq
        subst hq
        rw [htwo, Nat.add_zero, Nat.sub_self, Array.getElem!_set!_ne _ _ _ _ (by omega),
          Array.getElem!_set!_self _ _ _ (by rw [hsz4]; omega)]
      | succ k' =>
        rw [List.getElem?_cons_succ] at hq
        obtain ⟨hk', hd⟩ := hmid k' q hq
        have hlt := slot_lt_top' nvSt nvEn neSt E hv hE (i := 2 * nvSt + 2) (a := bump nvSt (2 * nvSt + 2) + k')
          (by omega) (by omega) (by rw [if_neg (by omega), hcntA, if_neg (by omega)]; omega)
        rw [Array.getElem!_set!_ne _ _ _ _ (by omega), Array.getElem!_set!_ne _ _ _ _ (by omega),
          show start nvSt neSt (allE nvSt nvEn E) (2 * nvSt + 2) + (k' + 1) - 2 * neSt =
            start nvSt neSt (allE nvSt nvEn E) (2 * nvSt + 2) + bump nvSt (2 * nvSt + 2) + k' - 2 * neSt by omega]
        exact hd
    · rw [if_neg hi', List.nil_append] at hq
      have hb : bump nvSt i = 0 := by simp [bump, hi']
      by_cases hi'' : i = 2 * nvEn - 1
      · subst hi''
        rw [if_pos rfl] at hq
        by_cases hk : k < ((E.reverse.filter (hits (2 * nvEn - 1))).map (other (2 * nvEn - 1))).length
        · rw [List.getElem?_append_left hk] at hq
          obtain ⟨hk', hd⟩ := hmid k q hq
          have hlt := slot_lt_top' nvSt nvEn neSt E hv hE (i := 2 * nvEn - 1) (a := k)
            (by omega) (by omega) (by rw [if_pos rfl, hcntA, if_pos rfl]; omega)
          have hpos := slot_pos nvSt nvEn neSt E hv (i := 2 * nvEn - 1) (by omega)
          rw [Array.getElem!_set!_ne _ _ _ _ (by omega), Array.getElem!_set!_ne _ _ _ _ (by omega),
            show start nvSt neSt (allE nvSt nvEn E) (2 * nvEn - 1) + k - 2 * neSt =
              start nvSt neSt (allE nvSt nvEn E) (2 * nvEn - 1) + bump nvSt (2 * nvEn - 1) + k - 2 * neSt by omega]
          exact hd
        · have hk0 := (List.getElem?_eq_some_iff.1 hq).1
          rw [List.getElem?_append_right (Nat.le_of_not_lt hk), List.getElem?_singleton] at hq
          rw [List.length_append, List.length_map, List.length_singleton] at hk0
          rw [List.length_map] at hk
          have hlen : k = cnt (2 * nvEn - 1) E.reverse := by
            change k = (E.reverse.filter (hits (2 * nvEn - 1))).length
            omega
          rw [show k - ((E.reverse.filter (hits (2 * nvEn - 1))).map (other (2 * nvEn - 1))).length = 0 by
            rw [List.length_map]; omega] at hq
          simp only [↓reduceIte, Option.some.injEq] at hq
          subst hq
          have hs := start_succ nvSt neSt (allE nvSt nvEn E) (2 * nvEn - 1) (by omega)
          rw [show 2 * nvEn - 1 + 1 = 2 * nvEn by omega, hend] at hs
          have hc := hcntA (2 * nvEn - 1)
          rw [if_pos rfl] at hc
          rw [show start nvSt neSt (allE nvSt nvEn E) (2 * nvEn - 1) + k - 2 * neSt = 2 * neEn - 1 - 2 * neSt by omega,
            Array.getElem!_set!_self _ _ _ (by rw [Array.size_set!, hsz4]; omega)]
      · rw [if_neg hi'', List.append_nil] at hq
        obtain ⟨hk', hd⟩ := hmid k q hq
        have hi3 : 2 * nvSt + 3 ≤ i := by
          rcases Nat.lt_or_ge i (2 * nvSt + 3) with h | h
          · have : i = 2 * nvSt + 1 := by omega
            subst this
            rw [cnt_reverse] at hk'
            have := hcntA (2 * nvSt + 1)
            rw [hlow, cnt_reverse] at this
            omega
          · exact h
        have hpos := slot_pos nvSt nvEn neSt E hv hi3
        have hlt := slot_lt_top' nvSt nvEn neSt E hv hE (i := i) (a := k) hi1 hi2
          (by rw [if_neg hi'', hcntA, if_neg hi'']; omega)
        rw [Array.getElem!_set!_ne _ _ _ _ (by omega), Array.getElem!_set!_ne _ _ _ _ (by omega),
          show start nvSt neSt (allE nvSt nvEn E) i + k - 2 * neSt =
            start nvSt neSt (allE nvSt nvEn E) i + bump nvSt i + k - 2 * neSt by omega]
        exact hd


/-! ### The bracket property -/

theorem rowDest_mem (nvSt nvEn : Nat) (E : List (Nat × Nat)) (i x : Nat)
    (hx : x ∈ rowDest nvSt nvEn E i) :
    (i = 2 * nvSt + 2 ∧ x = nvEn - 1) ∨ (∃ p ∈ E, hits i p = true ∧ x = other i p) ∨
      (i = 2 * nvEn - 1 ∧ x = nvSt) := by
  unfold rowDest at hx
  simp only [List.mem_append, List.mem_map, List.mem_filter, List.mem_reverse] at hx
  rcases hx with (hx | ⟨p, ⟨hp, hh⟩, rfl⟩) | hx
  · left
    split_ifs at hx with h
    · simp at hx; exact ⟨h, hx⟩
    · simp at hx
  · right; left; exact ⟨p, hp, hh, rfl⟩
  · right; right
    split_ifs at hx with h
    · simp at hx; exact ⟨h, hx⟩
    · simp at hx

theorem rowDest_lt (nvSt nvEn : Nat) (E : List (Nat × Nat)) (hv : nvSt + 2 ≤ nvEn)
    (hE : ∀ q ∈ E, nvSt ≤ q.1 ∧ q.1 < q.2 ∧ q.2 < nvEn) (nv x : Nat)
    (hx : x ∈ rowDest nvSt nvEn E (2 * nv + 1)) : x < nv := by
  rcases rowDest_mem nvSt nvEn E _ x hx with ⟨h, _⟩ | ⟨p, hp, hh, rfl⟩ | ⟨h, rfl⟩
  · omega
  · have := hE p hp
    rw [hits_iff] at hh
    rw [other_odd _ _ (by omega)]
    omega
  · omega

theorem rowDest_gt (nvSt nvEn : Nat) (E : List (Nat × Nat)) (hv : nvSt + 2 ≤ nvEn)
    (hE : ∀ q ∈ E, nvSt ≤ q.1 ∧ q.1 < q.2 ∧ q.2 < nvEn) (nv x : Nat) (hnv : nv < nvEn)
    (hx : x ∈ rowDest nvSt nvEn E (2 * nv + 2)) : nv < x := by
  rcases rowDest_mem nvSt nvEn E _ x hx with ⟨h, rfl⟩ | ⟨p, hp, hh, rfl⟩ | ⟨h, _⟩
  · omega
  · have := hE p hp
    rw [hits_iff] at hh
    rw [other_even _ _ (by omega)]
    omega
  · omega

theorem rowDest_pairwise (nvSt nvEn : Nat) (E : List (Nat × Nat))
    (hE : ∀ q ∈ E, nvSt ≤ q.1 ∧ q.1 < q.2 ∧ q.2 < nvEn)
    (hdom : E.Pairwise fun p q => p ≠ q → ¬ (q.1 ≤ p.1 ∧ q.2 ≤ p.2))
    (hnd : E.Nodup) (i : Nat) : (rowDest nvSt nvEn E i).Pairwise (· ≥ ·) := by
  have hdom' : E.Pairwise fun p q => ¬ (q.1 ≤ p.1 ∧ q.2 ≤ p.2) :=
    (hdom.and hnd).imp fun ⟨h, hne⟩ => h hne
  have hrev : E.reverse.Pairwise fun a b => ¬ (a.1 ≤ b.1 ∧ a.2 ≤ b.2) :=
    (List.pairwise_reverse (R := fun (a b : Nat × Nat) => ¬ (a.1 ≤ b.1 ∧ a.2 ≤ b.2))).2 hdom'
  have hM : ((E.reverse.filter (hits i)).map (other i)).Pairwise (· ≥ ·) := by
    rw [List.pairwise_map]
    refine List.Pairwise.imp_of_mem ?_ (hrev.filter (hits i))
    intro p q hp hq hpq
    rw [List.mem_filter] at hp hq
    have h1 := (hits_iff i p).1 hp.2
    have h2 := (hits_iff i q).1 hq.2
    unfold other
    split_ifs <;> omega
  have hmem : ∀ x ∈ (E.reverse.filter (hits i)).map (other i), nvSt ≤ x ∧ x ≤ nvEn - 1 := by
    intro x hx
    simp only [List.mem_map, List.mem_filter, List.mem_reverse] at hx
    obtain ⟨p, ⟨hp, _⟩, rfl⟩ := hx
    have := hE p hp
    unfold other
    split_ifs <;> omega
  unfold rowDest
  rw [List.pairwise_append, List.pairwise_append]
  refine ⟨⟨?_, hM, ?_⟩, ?_, ?_⟩
  · split_ifs <;> simp
  · intro a ha b hb
    split_ifs at ha with h
    · simp at ha; subst ha
      exact (hmem b hb).2
    · simp at ha
  · split_ifs <;> simp
  · intro a ha b hb
    split_ifs at hb with h
    · simp at hb; subst hb
      rw [List.mem_append] at ha
      rcases ha with ha | ha
      · split_ifs at ha with h'
        · simp at ha; omega
        · simp at ha
      · exact (hmem a ha).1
    · simp at hb

/-- Start of row `r` in a node's local layout: the local `adjBounds` entry, except that the
first row starts at the node's first node-edge slot. -/
def rowBound (nvSt neSt : Nat) (l : Layout) (r : Nat) : Nat :=
  if r = 2 * nvSt then 2 * neSt else l.adjBounds[r - 2 * nvSt]!

/-- Row `r` of a node's local layout. -/
def row (nvSt neSt : Nat) (l : Layout) (r : Nat) : List NodeAdj :=
  (List.range (rowBound nvSt neSt l (r + 1) - rowBound nvSt neSt l r)).map fun k =>
    l.adjDat[rowBound nvSt neSt l r + k - 2 * neSt]!

theorem row_destNv (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat)) (hv : nvSt + 2 ≤ nvEn)
    (hE : ∀ q ∈ E, nvSt ≤ q.1 ∧ q.1 < q.2 ∧ q.2 < nvEn) (hne : neEn = neSt + E.length + 1)
    (r : Nat) (hr1 : 2 * nvSt ≤ r) (hr2 : r ≤ 2 * nvEn - 1) :
    (row nvSt neSt (layoutNode .R node nvSt nvEn neSt neEn E) r).map (·.destNv) =
      rowDest nvSt nvEn E (r + 1) := by
  rw [layoutNode_R_eq _ _ _ _ _ _ hv]
  obtain ⟨hb, hd⟩ := run_spec node nvSt nvEn neSt neEn E hv hE hne
  have hbnd : rowBound nvSt neSt (run node nvSt nvEn neSt neEn E) r =
      start nvSt neSt (allE nvSt nvEn E) (r + 1) := by
    unfold rowBound
    split_ifs with h
    · subst h; rw [start_low _ _ _ _ (Nat.le_refl _)]
    · rw [hb r (by omega) (by omega)]
  have hbnd' : rowBound nvSt neSt (run node nvSt nvEn neSt neEn E) (r + 1) =
      start nvSt neSt (allE nvSt nvEn E) (r + 2) := by
    unfold rowBound
    rw [if_neg (by omega), hb (r + 1) (by omega) (by omega)]
  have hlen : rowBound nvSt neSt (run node nvSt nvEn neSt neEn E) (r + 1) -
      rowBound nvSt neSt (run node nvSt nvEn neSt neEn E) r = cnt (r + 1) (allE nvSt nvEn E) := by
    rw [hbnd, hbnd', start_succ _ _ _ _ (by omega)]; omega
  apply List.ext_getElem
  · simp [row, hlen, rowDest_length _ _ _ hv]
  · intro k h1 h2
    simp only [row, List.getElem_map, List.getElem_range]
    rw [hbnd]
    exact hd (r + 1) (by omega) (by omega) k _ (List.getElem?_eq_getElem h2)

/-- The reverse fill of `layoutNode` for an R node: with `E` in dominance order and each `(a, c)`
oriented `a < c`, row `2 nv` of the layout holds the lower incidences of `nv`, row `2 nv + 1` the
higher ones, each with non-increasing destination. -/
theorem layoutNode_R_bracket (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat))
    (hv : nvSt + 2 ≤ nvEn)
    (hE : ∀ q ∈ E, nvSt ≤ q.1 ∧ q.1 < q.2 ∧ q.2 < nvEn)
    (hdom : E.Pairwise fun p q => p ≠ q → ¬ (q.1 ≤ p.1 ∧ q.2 ≤ p.2))
    (hnd : E.Nodup) (hne : neEn = neSt + E.length + 1)
    (nv : Nat) (hnv1 : nvSt ≤ nv) (hnv2 : nv < nvEn) :
    let l := layoutNode .R node nvSt nvEn neSt neEn E
    (∀ a ∈ row nvSt neSt l (2 * nv), a.destNv < nv) ∧
    (∀ a ∈ row nvSt neSt l (2 * nv + 1), nv < a.destNv) ∧
    ((row nvSt neSt l (2 * nv)).map (·.destNv)).Pairwise (· ≥ ·) ∧
    ((row nvSt neSt l (2 * nv + 1)).map (·.destNv)).Pairwise (· ≥ ·) := by
  intro l
  have h0 := row_destNv node nvSt nvEn neSt neEn E hv hE hne (2 * nv) (by omega) (by omega)
  have h1 := row_destNv node nvSt nvEn neSt neEn E hv hE hne (2 * nv + 1) (by omega) (by omega)
  refine ⟨?_, ?_, ?_, ?_⟩
  · intro a ha
    exact rowDest_lt nvSt nvEn E hv hE nv a.destNv (h0 ▸ List.mem_map_of_mem ha)
  · intro a ha
    exact rowDest_gt nvSt nvEn E hv hE nv a.destNv hnv2 (h1 ▸ List.mem_map_of_mem ha)
  · rw [h0]; exact rowDest_pairwise nvSt nvEn E hE hdom hnd _
  · rw [h1]; exact rowDest_pairwise nvSt nvEn E hE hdom hnd _

end LayoutR

end Spqr
