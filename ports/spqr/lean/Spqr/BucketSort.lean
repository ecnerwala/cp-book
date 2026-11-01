/-!
# Stable counting sort

`bucketSort key bound l` distributes `l` into `bound + 1` buckets in one pass and concatenates them:
bucket `j < bound` holds the elements with `key x = j` in their original order, and bucket `bound`
holds every element with `key x ≥ bound` (also in original order). For keys below `bound` this is
exactly `List.mergeSort` by key, in time `O(bound + l.length)`.
-/

namespace Spqr

variable {α : Type u}

def bucketLists (key : α → Nat) (bound : Nat) (l : List α) : Array (List α) :=
  (l.foldl (fun bs x => bs.modify (min (key x) bound) (x :: ·)) (Array.replicate (bound + 1) [])).map
    List.reverse

def bucketSort (key : α → Nat) (bound : Nat) (l : List α) : List α :=
  (bucketLists key bound l).toList.flatten

section spec
variable {key : α → Nat} {bound : Nat}

private theorem foldl_modify_size (l : List α) (bs : Array (List α)) :
    (l.foldl (fun bs x => bs.modify (min (key x) bound) (x :: ·)) bs).size = bs.size := by
  induction l generalizing bs with
  | nil => rfl
  | cons x l ih => rw [List.foldl_cons, ih, Array.size_modify]

private theorem foldl_modify_getElem (l : List α) (bs : Array (List α)) (j : Nat) (hj : j < bs.size) :
    (l.foldl (fun bs x => bs.modify (min (key x) bound) (x :: ·)) bs)[j]'(by rwa [foldl_modify_size]) =
      (l.filter (fun x => min (key x) bound == j)).reverse ++ bs[j] := by
  induction l generalizing bs with
  | nil => simp
  | cons x l ih =>
    simp only [List.foldl_cons]
    rw [ih _ (by simpa using hj), Array.getElem_modify]
    by_cases hx : min (key x) bound = j <;> simp [hx]

theorem size_bucketLists (l : List α) : (bucketLists key bound l).size = bound + 1 := by
  simp [bucketLists, foldl_modify_size]

theorem getElem_bucketLists (l : List α) (j : Nat) (hj : j < (bucketLists key bound l).size) :
    (bucketLists key bound l)[j] = l.filter (fun x => min (key x) bound == j) := by
  rw [size_bucketLists] at hj
  simp only [bucketLists, Array.getElem_map]
  rw [foldl_modify_getElem (key := key) l _ j (by simpa using hj)]
  simp

theorem toList_bucketLists (l : List α) :
    (bucketLists key bound l).toList =
      (List.range (bound + 1)).map (fun j => l.filter (fun x => min (key x) bound == j)) := by
  apply List.ext_getElem
  · simp [size_bucketLists]
  · intro i h₁ h₂
    simp only [Array.getElem_toList, List.getElem_map, List.getElem_range]
    exact getElem_bucketLists l i _

theorem bucketSort_eq_flatMap (l : List α) :
    bucketSort key bound l =
      (List.range (bound + 1)).flatMap (fun j => l.filter (fun x => min (key x) bound == j)) := by
  rw [bucketSort, toList_bucketLists, List.flatMap_def]

private theorem flatMap_eq_single {β : Type v} {L : List Nat} {g : Nat → List β} {a : Nat}
    (ha : a ∈ L) (hn : L.Nodup) (h : ∀ b ∈ L, b ≠ a → g b = []) : L.flatMap g = g a := by
  induction L with
  | nil => simp at ha
  | cons b L ih =>
    rw [List.nodup_cons] at hn
    rw [List.flatMap_cons]
    by_cases hb : b = a
    · subst hb
      have : L.flatMap g = [] := List.flatMap_eq_nil_iff.2 fun c hc =>
        h c (List.mem_cons_of_mem _ hc) (fun e => hn.1 (e ▸ hc))
      simp [this]
    · rw [h b (List.mem_cons_self ..) hb, List.nil_append]
      exact ih (List.mem_of_ne_of_mem (Ne.symm hb) ha) hn.2 fun c hc => h c (List.mem_cons_of_mem _ hc)

/-- `bucketSort` is stable: it preserves the relative order of the elements of every class of a
predicate that is constant on (clamped) key classes. -/
theorem bucketSort_filter_of_const (l : List α) (p : α → Bool)
    (hp : ∀ x y, p x → p y → min (key x) bound = min (key y) bound) :
    (bucketSort key bound l).filter p = l.filter p := by
  rw [bucketSort_eq_flatMap, List.filter_flatMap]
  simp only [List.filter_filter]
  by_cases hl : l.filter p = []
  · rw [hl, List.flatMap_eq_nil_iff]
    intro j _
    rw [List.filter_eq_nil_iff] at hl ⊢
    intro x hx
    simp [hl x hx]
  · obtain ⟨x₀, hx₀, hpx₀⟩ := List.exists_mem_of_ne_nil _ hl |>.imp fun x h => List.mem_filter.1 h
    rw [flatMap_eq_single (L := List.range (bound + 1)) (a := min (key x₀) bound)
      (by simp only [List.mem_range]; omega) List.nodup_range]
    · apply List.filter_congr
      intro x _
      by_cases hx : p x
      · simp [hx, hp x x₀ hx hpx₀]
      · simp [hx]
    · intro b _ hb
      rw [List.filter_eq_nil_iff]
      intro x _
      by_cases hx : p x
      · simp [hx, hp x x₀ hx hpx₀, Ne.symm hb]
      · simp [hx]

theorem bucketSort_filter (l : List α) (k : Nat) :
    (bucketSort key bound l).filter (fun x => key x == k) = l.filter (fun x => key x == k) :=
  bucketSort_filter_of_const l _ fun x y hx hy => by
    simp only [beq_iff_eq] at hx hy
    rw [hx, hy]

theorem bucketSort_pairwise_min (l : List α) :
    (bucketSort key bound l).Pairwise (fun x y => min (key x) bound ≤ min (key y) bound) := by
  rw [bucketSort_eq_flatMap, List.pairwise_flatMap]
  refine ⟨fun a _ => ?_, ?_⟩
  · rw [List.pairwise_filter]
    apply List.pairwise_of_forall
    intro x y hx hy
    simp only [beq_iff_eq] at hx hy
    omega
  · refine List.pairwise_lt_range.imp fun {a b} hab x hx y hy => ?_
    simp only [List.mem_filter, beq_iff_eq] at hx hy
    omega

end spec

/-- A stable sort is unique: two lists sorted by `key` with the same elements of every key, in the
same order, are equal. -/
theorem eq_of_pairwise_of_filter_eq {key : α → Nat} :
    ∀ {l₁ l₂ : List α}, l₁.Pairwise (fun x y => key x ≤ key y) → l₂.Pairwise (fun x y => key x ≤ key y) →
      (∀ k, l₁.filter (fun x => key x == k) = l₂.filter (fun x => key x == k)) → l₁ = l₂
  | [], [], _, _, _ => rfl
  | [], y :: _, _, _, h => by simpa using h (key y)
  | x :: _, [], _, _, h => by simpa using h (key x)
  | x :: l₁, y :: l₂, h₁, h₂, h => by
    have le₁ : key x ≤ key y := by
      have hy : y ∈ (x :: l₁).filter (fun z => key z == key y) := by
        rw [h]; exact List.mem_filter.2 ⟨List.mem_cons_self .., by simp⟩
      rcases List.mem_cons.1 (List.mem_filter.1 hy).1 with hyx | hy
      · exact hyx ▸ Nat.le_refl _
      · exact (List.pairwise_cons.1 h₁).1 y hy
    have le₂ : key y ≤ key x := by
      have hx : x ∈ (y :: l₂).filter (fun z => key z == key x) := by
        rw [← h]; exact List.mem_filter.2 ⟨List.mem_cons_self .., by simp⟩
      rcases List.mem_cons.1 (List.mem_filter.1 hx).1 with hxy | hx
      · exact hxy ▸ Nat.le_refl _
      · exact (List.pairwise_cons.1 h₂).1 x hx
    have hk : key x = key y := Nat.le_antisymm le₁ le₂
    have hx := h (key x)
    rw [List.filter_cons_of_pos (by simp), List.filter_cons_of_pos (by simp [hk])] at hx
    obtain ⟨hxy, htl⟩ := List.cons.inj hx
    subst hxy
    congr
    refine eq_of_pairwise_of_filter_eq h₁.of_cons h₂.of_cons fun k => ?_
    by_cases hkx : k = key x
    · exact hkx ▸ htl
    · have := h k
      rwa [List.filter_cons_of_neg (by simp [Ne.symm hkx]), List.filter_cons_of_neg (by simp [Ne.symm hkx])] at this

/-- `mergeSort` is stable: it preserves the relative order of the elements of every class of a
predicate whose elements are pairwise `le`. -/
theorem _root_.List.mergeSort_filter_of_pairwise {le : α → α → Bool} (l : List α) (p : α → Bool)
    (hp : ∀ x y, p x → p y → le x y)
    (trans : ∀ a b c, le a b → le b c → le a c) (total : ∀ a b, le a b || le b a) :
    (l.mergeSort le).filter p = l.filter p := by
  have hsub : (l.filter p).Sublist ((l.mergeSort le).filter p) := by
    have := (List.sublist_mergeSort trans total
      (List.pairwise_filter.2 (List.pairwise_of_forall fun x y hx hy => hp x y hx hy))
      (List.filter_sublist (p := p) (l := l))).filter p
    rwa [List.filter_filter, List.filter_congr (fun x _ => Bool.and_self (p x))] at this
  exact (hsub.eq_of_length ((List.mergeSort_perm l le).filter p).length_eq.symm).symm

section mergeSort
variable {key : α → Nat} {bound : Nat}

theorem decide_le_trans (f : α → Nat) : ∀ a b c : α, decide (f a ≤ f b) → decide (f b ≤ f c) → decide (f a ≤ f c) := by
  intro a b c; simp only [decide_eq_true_eq]; omega
theorem decide_le_total (f : α → Nat) : ∀ a b : α, (decide (f a ≤ f b) || decide (f b ≤ f a)) = true := by
  intro a b; simp only [Bool.or_eq_true, decide_eq_true_eq]; omega

theorem bucketSort_eq_mergeSort_min (l : List α) :
    bucketSort key bound l = l.mergeSort fun x y => min (key x) bound ≤ min (key y) bound := by
  refine eq_of_pairwise_of_filter_eq (key := fun x => min (key x) bound) (bucketSort_pairwise_min l)
    ((List.pairwise_mergeSort (decide_le_trans _) (decide_le_total _) l).imp fun h => by simpa using h) fun k => ?_
  rw [List.mergeSort_filter_of_pairwise _ _ (fun x y hx hy => ?_) (decide_le_trans _) (decide_le_total _)]
  · exact bucketSort_filter_of_const l _ fun x y hx hy => by
      simp only [beq_iff_eq] at hx hy; omega
  · simp only [beq_iff_eq] at hx hy; simp; omega

theorem bucketSort_perm (l : List α) : (bucketSort key bound l).Perm l := by
  rw [bucketSort_eq_mergeSort_min]; exact List.mergeSort_perm _ _

theorem bucketSort_sorted (l : List α) (h : ∀ x ∈ l, key x < bound) :
    (bucketSort key bound l).Pairwise (fun x y => key x ≤ key y) :=
  (bucketSort_pairwise_min l).imp_of_mem fun {a b} ha hb hab => by
    have ha := h a ((bucketSort_perm l).mem_iff.1 ha)
    have hb := h b ((bucketSort_perm l).mem_iff.1 hb)
    omega

theorem bucketSort_stable (l : List α) (k : Nat) :
    (bucketSort key bound l).filter (fun x => key x == k) = l.filter (fun x => key x == k) :=
  bucketSort_filter l k

theorem bucketSort_eq_mergeSort (l : List α) (h : ∀ x ∈ l, key x < bound) :
    bucketSort key bound l = l.mergeSort fun x y => key x ≤ key y :=
  eq_of_pairwise_of_filter_eq (key := key) (bucketSort_sorted l h)
    ((List.pairwise_mergeSort (decide_le_trans _) (decide_le_total _) l).imp fun h => by simpa using h) fun k => by
    rw [bucketSort_stable, List.mergeSort_filter_of_pairwise _ _ (fun x y hx hy => ?_) (decide_le_trans _) (decide_le_total _)]
    simp only [beq_iff_eq] at hx hy; simp; omega

end mergeSort

end Spqr
