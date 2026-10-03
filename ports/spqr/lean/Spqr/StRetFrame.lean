import Spqr.StTree
import Spqr.StSimLemmas
import Spqr.StIParent
import Spqr.WalkCover

/-! # Frame facts of a returning `finishEdge` (`PROOF.md` §7)

A returning `finishEdge` (`o.cls = .ret lv kind`, `lv < d`) only touches the stack above the
untouched base `B` and the items it pulls out of those entries. `Fr U s₀ s` says the protected
items `U` (a set closed under children at `s₀`) kept their type and children, and no unprotected
item acquired a protected child; `Out U B s` says the stack is `above ++ B` with no protected id in
`above`. Both hold for `U` = the items below the spans of `B` (plus `vertItem curV` before the
vertex push), because every item on the stack has at most one owner (`cnt ≤ 1`). -/

namespace Spqr
open WalkM WalkState

theorem spansCount_pos_of_mem_readStack : ∀ {ts : List TEntry} {x : ItemId}, x ∈ readStack ts →
    0 < spansCount ts x
  | [], _, h => by simp [readStack, readL, readR] at h
  | t :: ts, x, h => by
    rw [spansCount_cons]
    rcases mem_readStack_cons.1 h with h | h
    · exact Nat.lt_of_lt_of_le (List.count_pos_iff.2 h) (Nat.le_add_right _ _)
    · exact Nat.lt_of_lt_of_le (spansCount_pos_of_mem_readStack h) (Nat.le_add_left _ _)

theorem lt_of_mem_ch {items : Items} {j c : ItemId} (hc : c ∈ Items.ch items j) : j < items.size := by
  by_contra h
  rw [Items.ch_of_le _ _ (Nat.le_of_not_lt h)] at hc
  cases hc

theorem chCount_pos_of_mem_ch {items : Items} {j c : ItemId} (hc : c ∈ Items.ch items j) :
    0 < chCount items c :=
  Nat.lt_of_lt_of_le (List.count_pos_iff.2 hc) (Items.count_le_chCount _ (lt_of_mem_ch hc) c)

theorem mem_getSide' {p : List ItemId × List ItemId} {dir : Bool} {c : ItemId} (hc : c ∈ getSide p dir) :
    c ∈ p.1 ++ p.2 := by
  cases dir <;> simp [getSide] at hc <;> simp [hc]

/-- The protected set: allocated at `s₀`, and not the root. -/
structure UOk (U : ItemId → Prop) (s₀ : WalkState) : Prop where
  lt : ∀ y : Nat, U y → y < s₀.items.size
  root : ¬ U 0

theorem UOk.head {U : ItemId → Prop} {s₀ : WalkState} (hU : UOk U s₀) {l : List ItemId}
    (hl : ∀ i ∈ l, ¬ U i) : ¬ U l.head! := by
  cases l with
  | nil => exact hU.root
  | cons x _ => exact hl x List.mem_cons_self

theorem UOk.head_getSide {U : ItemId → Prop} {s₀ : WalkState} (hU : UOk U s₀)
    {p : List ItemId × List ItemId} (dir : Bool) (hp : ∀ i ∈ p.1 ++ p.2, ¬ U i) :
    ¬ U (getSide p dir).head! := by
  cases dir
  · exact hU.head fun i hi => hp i (List.mem_append_left _ hi)
  · exact hU.head fun i hi => hp i (List.mem_append_right _ hi)

/-- The protected items keep type and children; no unprotected item gets a protected child. -/
structure Fr (U : ItemId → Prop) (s₀ s : WalkState) : Prop where
  size : s₀.items.size ≤ s.items.size
  type : ∀ y : Nat, U y → Items.type s.items y = Items.type s₀.items y
  ch : ∀ y : Nat, U y → Items.ch s.items y = Items.ch s₀.items y
  out : ∀ i : Nat, ¬ U i → ∀ c ∈ Items.ch s.items i, ¬ U c

/-- The stack is `above ++ B` with no protected id in `above`. -/
def Out (U : ItemId → Prop) (B : List TEntry) (s : WalkState) : Prop :=
  ∃ above, s.tstack = above ++ B ∧ ∀ t ∈ above, ∀ i ∈ t.spans.1 ++ t.spans.2, ¬ U i

structure FrO (U : ItemId → Prop) (B : List TEntry) (s₀ s : WalkState) : Prop where
  fr : Fr U s₀ s
  out : Out U B s

variable {U : ItemId → Prop} {B : List TEntry} {s₀ s s' : WalkState}

theorem Fr.refl (h : ∀ i : Nat, ¬ U i → ∀ c ∈ Items.ch s₀.items i, ¬ U c) : Fr U s₀ s₀ :=
  ⟨Nat.le_refl _, fun _ _ => rfl, fun _ _ => rfl, h⟩

theorem Fr.of_items (h : Fr U s₀ s) (he : s'.items = s.items) : Fr U s₀ s' := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> rw [he]
  exacts [h.size, h.type, h.ch, h.out]

theorem Fr.push (hU : UOk U s₀) (h : Fr U s₀ s) (x : Item) (hx : x.ch = []) :
    Fr U s₀ { s with items := s.items.push x } := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · show s₀.items.size ≤ (s.items.push x).size
    simp only [Array.size_push]; exact Nat.le_succ_of_le h.size
  · intro y hy
    have hne : y ≠ s.items.size := by have h1 := hU.lt y hy; have h2 := h.size; omega
    show Items.type (s.items.push x) y = _
    rw [Items.type_push, if_neg hne]
    exact h.type y hy
  · intro y hy
    have hne : y ≠ s.items.size := by have h1 := hU.lt y hy; have h2 := h.size; omega
    show Items.ch (s.items.push x) y = _
    rw [Items.ch_push, if_neg hne]
    exact h.ch y hy
  · intro i hi c hc
    change c ∈ Items.ch (s.items.push x) i at hc
    rw [Items.ch_push] at hc
    split at hc
    · rw [hx] at hc; cases hc
    · exact h.out i hi c hc

theorem Fr.modify_keep (h : Fr U s₀ s) (j : ItemId) (f : Item → Item)
    (ht : ∀ it, (f it).type = it.type) (hc : ∀ it, (f it).ch = it.ch) :
    Fr U s₀ { s with items := s.items.modify j f } := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · show s₀.items.size ≤ (s.items.modify j f).size
    simp only [Array.size_modify]; exact h.size
  · intro y hy; show Items.type (s.items.modify j f) y = _
    rw [Items.type_modify _ _ _ _ ht]; exact h.type y hy
  · intro y hy; show Items.ch (s.items.modify j f) y = _
    rw [Items.ch_modify_of_ch _ _ _ _ hc]; exact h.ch y hy
  · intro i hi c hc'; change c ∈ Items.ch (s.items.modify j f) i at hc'
    rw [Items.ch_modify_of_ch _ _ _ _ hc] at hc'; exact h.out i hi c hc'

theorem Fr.modify_ch (h : Fr U s₀ s) (j : ItemId) (f : Item → Item) (ht : ∀ it, (f it).type = it.type)
    (hj : ¬ U j) (hnew : ∀ c ∈ (f s.items[j]!).ch, ¬ U c) :
    Fr U s₀ { s with items := s.items.modify j f } := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · show s₀.items.size ≤ (s.items.modify j f).size
    simp only [Array.size_modify]; exact h.size
  · intro y hy; show Items.type (s.items.modify j f) y = _
    rw [Items.type_modify _ _ _ _ ht]; exact h.type y hy
  · intro y hy; show Items.ch (s.items.modify j f) y = _
    rw [Items.ch_modify_ne _ _ _ _ (by rintro rfl; exact hj hy)]
    exact h.ch y hy
  · intro i hi c hc
    change c ∈ Items.ch (s.items.modify j f) i at hc
    rcases Items.IsParent_modify hc with hc | ⟨rfl, hlt, hc⟩
    · exact h.out i hi c hc
    · exact hnew c (by rwa [getElem!_pos s.items i hlt])

theorem Out.of_tstack (h : Out U B s) (he : s'.tstack = s.tstack) : Out U B s' := by
  obtain ⟨a, ha, hu⟩ := h
  exact ⟨a, he.trans ha, hu⟩

theorem Out.length (h : Out U B s) : B.length ≤ s.tstack.length := by
  obtain ⟨a, ha, -⟩ := h
  rw [ha, List.length_append]; omega

theorem Out.cons (h : Out U B s) {t : TEntry} (he : s'.tstack = t :: s.tstack)
    (ht : ∀ i ∈ t.spans.1 ++ t.spans.2, ¬ U i) : Out U B s' := by
  obtain ⟨a, ha, hu⟩ := h
  refine ⟨t :: a, by rw [he, ha]; rfl, fun u hu' => ?_⟩
  rcases List.mem_cons.1 hu' with rfl | hu'
  · exact ht
  · exact hu u hu'

theorem Out.merge (h : Out U B s) (hl : B.length + 2 ≤ s.tstack.length)
    (he : s'.tstack = mergeTop s.tstack) : Out U B s' := by
  obtain ⟨abv, ha, hu⟩ := h
  rw [ha] at hl
  rcases abv with _ | ⟨b, _ | ⟨a, rest⟩⟩
  · rw [List.nil_append] at hl; omega
  · rw [List.singleton_append, List.length_cons] at hl; omega
  · refine ⟨{ a with topDepth := min a.topDepth b.topDepth, spans := (b.spans.1 ++ a.spans.1, a.spans.2 ++ b.spans.2) } :: rest, by rw [he, ha]; rfl, fun u hu' => ?_⟩
    rcases List.mem_cons.1 hu' with rfl | hu'
    · intro i hi
      simp only [List.mem_append] at hi
      rcases hi with (hi | hi) | (hi | hi)
      · exact hu b (by simp) i (List.mem_append_left _ hi)
      · exact hu a (by simp) i (List.mem_append_left _ hi)
      · exact hu a (by simp) i (List.mem_append_right _ hi)
      · exact hu b (by simp) i (List.mem_append_right _ hi)
    · exact hu u (by simp [hu'])

theorem Out.modifyCur (h : Out U B s) (hl : B.length + 1 ≤ s.tstack.length) (f : TEntry → TEntry)
    (he : s'.tstack = match s.tstack with | a :: rest => f a :: rest | [] => [])
    (hf : ∀ t, (∀ i ∈ t.spans.1 ++ t.spans.2, ¬ U i) → ∀ i ∈ (f t).spans.1 ++ (f t).spans.2, ¬ U i) :
    Out U B s' := by
  obtain ⟨abv, ha, hu⟩ := h
  rcases abv with _ | ⟨a, rest⟩
  · rw [ha, List.nil_append] at hl; omega
  · rw [ha] at he
    refine ⟨f a :: rest, he, fun u hu' => ?_⟩
    rcases List.mem_cons.1 hu' with rfl | hu'
    · exact hf a (hu a (by simp))
    · exact hu u (by simp [hu'])

theorem Out.modifyNxt (h : Out U B s) (hl : B.length + 2 ≤ s.tstack.length) (f : TEntry → TEntry)
    (he : s'.tstack = match s.tstack with | a :: b :: rest => a :: f b :: rest | l => l)
    (hf : ∀ t, (∀ i ∈ t.spans.1 ++ t.spans.2, ¬ U i) → ∀ i ∈ (f t).spans.1 ++ (f t).spans.2, ¬ U i) :
    Out U B s' := by
  obtain ⟨abv, ha, hu⟩ := h
  rw [ha] at hl
  rcases abv with _ | ⟨a, _ | ⟨b, rest⟩⟩
  · rw [List.nil_append] at hl; omega
  · rw [List.singleton_append, List.length_cons] at hl; omega
  · rw [ha] at he
    refine ⟨a :: f b :: rest, he, fun u hu' => ?_⟩
    simp only [List.mem_cons] at hu'
    rcases hu' with rfl | rfl | hu'
    · exact hu u (by simp)
    · exact hf b (hu b (by simp))
    · exact hu u (by simp [hu'])

theorem Out.head (h : Out U B s) (hl : B.length + 1 ≤ s.tstack.length) :
    ∀ i ∈ s.tstack.head!.spans.1 ++ s.tstack.head!.spans.2, ¬ U i := by
  obtain ⟨abv, ha, hu⟩ := h
  rcases abv with _ | ⟨a, rest⟩
  · rw [ha, List.nil_append] at hl; omega
  · rw [ha]; exact hu a (by simp)

theorem Out.nxt (h : Out U B s) (hl : B.length + 2 ≤ s.tstack.length) :
    ∀ i ∈ s.tstack.tail.head!.spans.1 ++ s.tstack.tail.head!.spans.2, ¬ U i := by
  obtain ⟨abv, ha, hu⟩ := h
  rw [ha] at hl
  rcases abv with _ | ⟨a, _ | ⟨b, rest⟩⟩
  · rw [List.nil_append] at hl; omega
  · rw [List.singleton_append, List.length_cons] at hl; omega
  · rw [ha]; exact hu b (by simp)

theorem length_mergeTop (l : List TEntry) : (mergeTop l).length = l.length - 1 := by
  rcases l with _ | ⟨b, _ | ⟨a, rest⟩⟩ <;> simp [mergeTop]

theorem length_modifyCur (f : TEntry → TEntry) (l : List TEntry) :
    (match l with | a :: rest => f a :: rest | [] => []).length = l.length := by
  cases l <;> rfl

theorem length_modifyNxt (f : TEntry → TEntry) (l : List TEntry) :
    (match l with | a :: b :: rest => a :: f b :: rest | l => l).length = l.length := by
  rcases l with _ | ⟨a, _ | ⟨b, rest⟩⟩ <;> rfl

/-! ## Primitives -/

theorem FrO.of_eq (h : FrO U B s₀ s) (hi : s'.items = s.items) (ht : s'.tstack = s.tstack) :
    FrO U B s₀ s' := ⟨h.fr.of_items hi, h.out.of_tstack ht⟩

theorem FrO.merge (h : FrO U B s₀ s) (hl : B.length + 2 ≤ s.tstack.length) :
    FrO U B s₀ { s with tstack := mergeTop s.tstack } := ⟨h.fr.of_items rfl, h.out.merge hl rfl⟩

theorem FrO.merge' (h : FrO U B s₀ s) (hl : B.length + 2 ≤ s.tstack.length) (hi : s'.items = s.items)
    (ht : s'.tstack = mergeTop s.tstack) : FrO U B s₀ s' := ⟨h.fr.of_items hi, h.out.merge hl ht⟩

theorem FrO.cons (h : FrO U B s₀ s) {t : TEntry} (ht : ∀ i ∈ t.spans.1 ++ t.spans.2, ¬ U i) :
    FrO U B s₀ { s with tstack := t :: s.tstack } := ⟨h.fr.of_items rfl, h.out.cons rfl ht⟩

theorem length_maybeUnwrapNxt (ty : NodeType) :
    wp (maybeUnwrapNxt ty) (fun _ s' => s'.tstack.length = s.tstack.length) s := by
  unfold maybeUnwrapNxt
  simp only [wp_bind, wp_ite, wp_get, wp_allocItem, wp_pure, wp_nxt, wp_stackDir, wp_getItem,
    wp_modifyNxt]
  split
  · trivial
  · split
    · exact length_modifyNxt _ _
    · trivial

theorem length_finishTstackTop (item : ItemId) :
    wp (finishTstackTop item) (fun _ s' => s'.tstack.length = s.tstack.length) s := by
  unfold finishTstackTop
  simp only [wp_bind, wp_cur, wp_stackDir, wp_makeVs, wp_modifyItem, wp_modifyCur]
  exact length_modifyCur _ _

theorem length_loop1Type (d : Nat) (dir : Bool) :
    wp (loop1Type d dir) (fun _ s' => s'.tstack.length ≤ s.tstack.length) s := by
  unfold loop1Type
  simp only [wp_bind, wp_ite, wp_nxt, wp_cur, wp_setStackDir, wp_mergeTstackTops, wp_pure]
  split
  · simp only [length_mergeTop]; omega
  · split <;> exact Nat.le_refl _

theorem length_loop1Body (d : Nat) (dir : Bool) :
    wp (loop1Body d dir) (fun _ s' => s'.tstack.length ≤ s.tstack.length) s := by
  unfold loop1Body
  simp only [wp_bind]
  refine wp_mono _ (length_loop1Type d dir) fun ty s₁ h₁ => ?_
  refine wp_mono _ (length_maybeUnwrapNxt ty) fun item s₂ h₂ => ?_
  rw [wp_mergeTstackTops]
  refine wp_mono _ (length_finishTstackTop item) fun _ s' h' => ?_
  rw [h']; simp only [length_mergeTop]; omega

theorem FrO.unwrapNxt (hU : UOk U s₀) (h : FrO U B s₀ s) (hl : B.length + 2 ≤ s.tstack.length)
    (ty : NodeType) :
    wp (maybeUnwrapNxt ty)
      (fun item s' => FrO U B s₀ s' ∧ ¬ U item ∧ s'.tstack.length = s.tstack.length) s := by
  have alloc : FrO U B s₀ { s with items := s.items.push ⟨ty, (none, none), []⟩ } ∧ ¬ U s.items.size :=
    ⟨⟨h.fr.push hU _ rfl, h.out.of_tstack rfl⟩,
      fun hu => by have h1 := hU.lt _ hu; have h2 := h.fr.size; omega⟩
  unfold maybeUnwrapNxt
  simp only [wp_bind, wp_ite, wp_get, wp_allocItem, wp_pure, wp_nxt, wp_stackDir, wp_getItem,
    wp_modifyNxt]
  split
  · exact ⟨alloc.1, alloc.2, trivial⟩
  · split
    · have hb := h.out.nxt hl
      have hitem := hU.head_getSide s.stackDir[s.tstack.tail.head!.topDepth]! hb
      refine ⟨⟨h.fr.of_items rfl, h.out.modifyNxt hl _ rfl (fun t _ c hc => ?_)⟩, hitem,
        length_modifyNxt _ _⟩
      simp only [mem_setSides_append, List.append_nil] at hc
      refine h.fr.out _ hitem c ?_
      by_cases hlt : (getSide s.tstack.tail.head!.spans s.stackDir[s.tstack.tail.head!.topDepth]!).head!
          < s.items.size
      · rwa [Items.ch_getElem! _ _ hlt] at hc
      · rw [getElem!_neg s.items _ hlt, show (default : Item).ch = [] from rfl] at hc
        cases hc
    · exact ⟨alloc.1, alloc.2, trivial⟩

theorem FrO.finishTop (h : FrO U B s₀ s) (hl : B.length + 1 ≤ s.tstack.length) {item : ItemId}
    (hitem : ¬ U item) :
    wp (finishTstackTop item) (fun _ s' => FrO U B s₀ s' ∧ s'.tstack.length = s.tstack.length) s := by
  unfold finishTstackTop
  simp only [wp_bind, wp_cur, wp_stackDir, wp_makeVs, wp_modifyItem, wp_modifyCur]
  have hhead := h.out.head hl
  have hfr : Fr U s₀ { s with items := s.items.modify item fun it =>
      { it with vs := setSides s.stackDir[s.tstack.head!.topDepth]! (some s.stackVerts[s.tstack.head!.topDepth]!)
                        (some s.tstack.head!.vStart),
                ch := getSide s.tstack.head!.spans s.stackDir[s.tstack.head!.topDepth]! } } :=
    h.fr.modify_ch item _ (fun _ => rfl) hitem (fun c hc => hhead c (mem_getSide' hc))
  refine ⟨⟨hfr.of_items rfl, h.out.modifyCur hl _ rfl (fun t _ i hi => ?_)⟩, length_modifyCur _ _⟩
  simp only [mem_setSides_append, List.append_nil, List.mem_singleton] at hi
  exact hi ▸ hitem

theorem FrO.l1Type (h : FrO U B s₀ s) (d : Nat) (dir : Bool) :
    wp (loop1Type d dir) (fun _ s' => (FrO U B s₀ s' ∧ B.length + 2 ≤ s'.tstack.length) ∨
      s'.tstack.length ≤ B.length + 1) s := by
  unfold loop1Type
  simp only [wp_bind, wp_ite, wp_nxt, wp_cur, wp_setStackDir, wp_mergeTstackTops, wp_pure]
  split
  · by_cases hl : B.length + 3 ≤ s.tstack.length
    · refine Or.inl ⟨h.merge' (by omega) rfl rfl, ?_⟩
      simp only [length_mergeTop]; omega
    · exact Or.inr (by simp only [length_mergeTop]; omega)
  · split
    · exact if hl : B.length + 2 ≤ s.tstack.length then Or.inl ⟨h, hl⟩ else Or.inr (by omega)
    · exact if hl : B.length + 2 ≤ s.tstack.length then Or.inl ⟨h, hl⟩ else Or.inr (by omega)

theorem FrO.l1Body (hU : UOk U s₀) (h : FrO U B s₀ s) (d : Nat) (dir : Bool) :
    wp (loop1Body d dir) (fun _ s' => FrO U B s₀ s' ∨ s'.tstack.length ≤ B.length) s := by
  unfold loop1Body
  simp only [wp_bind]
  refine wp_mono _ (h.l1Type d dir) fun ty s₁ h₁ => ?_
  rcases h₁ with ⟨h₁, hl₁⟩ | hl₁
  · refine wp_mono _ (h₁.unwrapNxt hU hl₁ ty) fun item s₂ ⟨h₂, hitem, hl₂⟩ => ?_
    rw [wp_mergeTstackTops]
    refine wp_mono _ ((h₂.merge (by omega)).finishTop (by simp only [length_mergeTop]; omega) hitem)
      fun _ s' h' => Or.inl h'.1
  · refine wp_mono _ (length_maybeUnwrapNxt ty) fun item s₂ hl₂ => ?_
    rw [wp_mergeTstackTops]
    refine wp_mono _ (length_finishTstackTop item) fun _ s' h' => Or.inr ?_
    rw [h']; simp only [length_mergeTop]; omega

/-! ## Loops -/

theorem wp_cond_result (cond : WalkM Bool) (h : wp cond (fun _ s' => s' = s) s) :
    wp cond (fun b s' => s' = s ∧ b = result cond s) s := ⟨h, rfl⟩

theorem wp_loop_cond (I : WalkState → Prop) (n : Nat) (cond : WalkM Bool) (body : WalkM Unit)
    (hc : ∀ s, wp cond (fun _ s' => s' = s) s)
    (hb : ∀ s, I s → result cond s = true → wp body (fun _ s' => I s') s) (h : I s) :
    wp (loop n cond body) (fun _ s' => I s') s := by
  induction n generalizing s with
  | zero => exact h
  | succ n ih =>
    simp only [loop, wp_bind, wp_ite, wp_pure]
    refine wp_mono _ (wp_cond_result cond (hc s)) fun b s₁ ⟨e, hb'⟩ => ?_
    subst e; subst hb'
    split
    · exact wp_mono _ (hb s₁ h ‹_›) fun _ s₂ h₂ => ih h₂
    · exact h

theorem wp_loop_or (I : WalkState → Prop) (m n : Nat) (cond : WalkM Bool) (body : WalkM Unit)
    (hc : ∀ s, wp cond (fun _ s' => s' = s) s)
    (hlen : ∀ s, wp body (fun _ s' => s'.tstack.length ≤ s.tstack.length) s)
    (hb : ∀ s, I s → wp body (fun _ s' => I s' ∨ s'.tstack.length ≤ m) s)
    (h : I s ∨ s.tstack.length ≤ m) :
    wp (loop n cond body) (fun _ s' => I s' ∨ s'.tstack.length ≤ m) s :=
  wp_loop (fun s => I s ∨ s.tstack.length ≤ m) n cond body
    (fun s hs => wp_mono _ (hc s) fun _ _ e => e ▸ hs)
    (fun s hs => hs.elim (hb s) fun hm => wp_mono _ (hlen s) fun _ _ h' => Or.inr (h'.trans hm))
    h (fun _ => id)

theorem loop3Cond_result (o : Nat) (s : WalkState) (h : result (loop3Cond o) s = true) :
    o + 3 < s.tstack.length := by
  have : result (loop3Cond o) s = decide (s.tstack.length > o + 3) := rfl
  rw [this, decide_eq_true_eq] at h
  exact h

theorem FrO.ears (hU : UOk U s₀) (h : FrO U B s₀ s) {nxtV d e : Nat} {dir : Bool}
    (hq : ¬ U (edgeItem s.g e)) :
    wp (closeEars nxtV d e dir) (fun _ s' => FrO U B s₀ s' ∨ s'.tstack.length ≤ B.length) s := by
  unfold closeEars
  simp only [wp_bind, wp_pushEdgeTstack, wp_tstackSize]
  refine wp_loop_or (FrO U B s₀) B.length _ _ _ (fun _ => loop1Cond_pure d) (fun s => length_loop1Body d dir)
    (fun s hs => hs.l1Body hU d dir) (Or.inl (h.cons fun i hi => ?_))
  rw [mem_setSides, List.mem_singleton] at hi
  exact hi ▸ hq

theorem length_mergeLate (d : Nat) :
    wp (mergeLate d) (fun _ s' => s'.tstack.length ≤ s.tstack.length) s := by
  unfold mergeLate
  simp only [wp_bind, wp_get, wp_cur, wp_ite, wp_tstackSize, wp_pure]
  split
  · refine wp_loop (fun s' => s'.tstack.length ≤ s.tstack.length) _ _ _
      (fun s hs => wp_mono _ (loop2Cond_pure _) fun _ _ e => e ▸ hs)
      (fun s hs => by rw [wp_mergeTstackTops]; simp only [length_mergeTop]; omega) (Nat.le_refl _)
      (fun _ => id)
  · exact Nat.le_refl _

theorem FrO.late (h : FrO U B s₀ s) (d : Nat) :
    wp (mergeLate d) (fun _ s' => FrO U B s₀ s' ∨ s'.tstack.length ≤ B.length) s := by
  unfold mergeLate
  simp only [wp_bind, wp_get, wp_cur, wp_ite, wp_tstackSize, wp_pure]
  split
  · refine wp_loop_or (FrO U B s₀) B.length _ _ _ (fun _ => loop2Cond_pure _)
      (fun s => by rw [wp_mergeTstackTops]; simp only [length_mergeTop]; omega)
      (fun s hs => ?_) (Or.inl h)
    rw [wp_mergeTstackTops]
    by_cases hl : B.length + 2 ≤ s.tstack.length
    · exact Or.inl (hs.merge hl)
    · exact Or.inr (by simp only [length_mergeTop]; omega)
  · exact Or.inl h

/-! ## `finishP`, `finishTail`, `closeVert'` -/

theorem FrO.finP (hU : UOk U s₀) (h : FrO U B s₀ s) (hl : B.length + 1 ≤ s.tstack.length)
    {curV lv : Nat} {ty1 : Bool} (hB : ∀ t ∈ B, t.vStart ≠ curV) :
    wp (finishP curV lv ty1) (fun _ s' => FrO U B s₀ s' ∧ B.length + 1 ≤ s'.tstack.length) s := by
  unfold finishP condP
  simp only [wp_bind, wp_ite, wp_tstackSize, wp_nxt, wp_pure]
  split
  · rename_i hc
    simp only [Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq] at hc
    have hl2 : B.length + 2 ≤ s.tstack.length := by
      by_contra hlt
      obtain ⟨abv, ha, -⟩ := h.out
      rcases abv with _ | ⟨a, _ | ⟨b, rest⟩⟩
      · rw [ha, List.nil_append] at hl; omega
      · rcases B with _ | ⟨t, B⟩
        · have := hc.1.1.2; rw [ha] at this; simp at this
        · refine hB t List.mem_cons_self ?_
          have : s.tstack.tail.head! = t := by rw [ha]; rfl
          rw [this] at hc; exact hc.1.2
      · rw [ha] at hlt; simp at hlt
    refine wp_mono _ (h.unwrapNxt hU hl2 .P) fun item s₁ ⟨h₁, hitem, hl₁⟩ => ?_
    rw [wp_mergeTstackTops]
    refine wp_mono _ ((h₁.merge (by omega)).finishTop (by simp only [length_mergeTop]; omega) hitem)
      fun _ s' ⟨h', hl'⟩ => ⟨h', by rw [hl']; simp only [length_mergeTop]; omega⟩
  · exact ⟨h, hl⟩

theorem FrO.tail (h : FrO U B s₀ s) (hl : B.length + 1 ≤ s.tstack.length) {curV d : Nat}
    (isSingle : Bool) (hv : ¬ U (vertItem curV)) :
    wp (finishTail curV d false isSingle) (fun _ s' => FrO U B s₀ s') s := by
  unfold finishTail
  simp only [Bool.not_false, ↓reduceIte, wp_bind, wp_pushVertTstack, wp_ite, wp_mergeTstackTops, wp_pure]
  have h₁ := h.cons (t := ⟨curV, d, s.nxtEdgeIdx, setSides s.stackDir[d]! [vertItem curV] []⟩)
    fun i hi => by rw [mem_setSides, List.mem_singleton] at hi; exact hi ▸ hv
  split
  · exact h₁.merge (by simp only [List.length_cons]; omega)
  · exact h₁

theorem after_finishTail_true {curV d : Nat} (isSingle : Bool) (s : WalkState) :
    after (finishTail curV d true isSingle) s = s := rfl

theorem FrO.vPre (h : FrO U B s₀ s) {orig : Nat} (hB : B.length ≤ orig)
    (hl : orig + 3 ≤ s.tstack.length) (ty1 isSingle : Bool) :
    wp (vertPre ty1 orig isSingle) (fun _ s' => FrO U B s₀ s' ∧ orig + 3 ≤ s'.tstack.length) s := by
  unfold vertPre
  simp only [wp_ite, wp_bind, wp_tstackSize, wp_pure]
  split
  · refine wp_loop_cond (fun s => FrO U B s₀ s ∧ orig + 3 ≤ s.tstack.length) _ _ _ (fun _ => loop3Cond_pure orig)
      (fun s ⟨hs, hls⟩ hc => ?_) ⟨h, hl⟩
    have hc' := loop3Cond_result orig s hc
    rw [wp_mergeTstackTops]
    exact ⟨hs.merge (by omega), by simp only [length_mergeTop]; omega⟩
  · exact ⟨h, hl⟩

theorem FrO.vUnwrap (hU : UOk U s₀) (h : FrO U B s₀ s) (hl : B.length + 2 ≤ s.tstack.length)
    (ty1 isSingle : Bool) :
    wp (vertUnwrap ty1 isSingle)
      (fun item s' => FrO U B s₀ s' ∧ (∀ i ∈ item, ¬ U i) ∧ s'.tstack.length = s.tstack.length) s := by
  unfold vertUnwrap
  simp only [wp_ite, wp_map, wp_pure]
  split
  · exact wp_mono _ (h.unwrapNxt hU hl _) fun item s' ⟨h', hi, hl'⟩ =>
      ⟨h', fun j hj => by rw [Option.mem_some_iff] at hj; exact hj ▸ hi, hl'⟩
  · exact ⟨h, fun _ h => (Option.not_mem_none _ h).elim, trivial⟩

theorem FrO.closeV (hU : UOk U s₀) (h : FrO U B s₀ s) {orig : Nat} (hB : B.length ≤ orig)
    (hl : orig + 3 ≤ s.tstack.length) (curV : Nat) (dir ty1 isSingle : Bool) :
    wp (closeVert' curV dir ty1 orig isSingle)
      (fun _ s' => FrO U B s₀ s' ∧ B.length + 1 ≤ s'.tstack.length) s := by
  unfold closeVert'
  simp only [wp_bind]
  refine wp_mono _ (h.vPre hB hl ty1 isSingle) fun b s₁ ⟨h₁, hl₁⟩ => ?_
  refine wp_mono _ (h₁.vUnwrap hU (by omega) ty1 b) fun item s₂ ⟨h₂, hitem, hl₂⟩ => ?_
  rw [wp_mergeTstackTops, wp_mergeTstackTops]
  have h₃ : FrO U B s₀ { s₂ with tstack := mergeTop (mergeTop s₂.tstack) } :=
    (h₂.merge (by omega)).merge' (by simp only [length_mergeTop]; omega) rfl rfl
  have hl₃ : B.length + 1 ≤ (mergeTop (mergeTop s₂.tstack)).length := by
    simp only [length_mergeTop]; omega
  unfold retarget
  rw [wp_modifyCur]
  have h₄ : FrO U B s₀ { s₂ with tstack := match mergeTop (mergeTop s₂.tstack) with
      | a :: rest => { a with vStart := curV, spans := setSides (!dir) (a.spans.1 ++ a.spans.2) [] } :: rest
      | [] => [] } :=
    ⟨h₃.fr.of_items rfl, h₃.out.modifyCur hl₃ _ rfl (fun t ht i hi => by
      rw [mem_setSides] at hi; exact ht i hi)⟩
  have hl₄ : B.length + 1 ≤ (match mergeTop (mergeTop s₂.tstack) with
      | a :: rest => { a with vStart := curV, spans := setSides (!dir) (a.spans.1 ++ a.spans.2) [] } :: rest
      | [] => []).length := by
    rw [length_modifyCur]; exact hl₃
  cases item with
  | some item =>
    simp only [vertFinish, wp_bind, wp_pure]
    exact wp_mono _ (h₄.finishTop hl₄ (hitem item rfl)) fun _ s' ⟨h', hl'⟩ => ⟨h', by rw [hl']; exact hl₄⟩
  | none => exact ⟨h₄, hl₄⟩

/-! ## The returning `finishEdge` -/

/-- A returning `finishEdge` keeps `FrO` up to the vertex push (`fePState`) and, if `vertItem curV`
is unprotected, to the end. -/
theorem FrO.finishEdge_ret {d lv : Nat} {kind : RetKind} {o : DfsOut} {curV : Nat} {hasVert : Bool}
    {sub pre : List TEntry} (hU : UOk U s) (h : FrO U B s s)
    (hE : s.EarFinish curV d o hasVert sub (pre ++ B)) (ho : o.cls = .ret lv kind) (hlow : lv < d)
    (hq : ¬ U (edgeItem s.g o.e)) (hB : ∀ t ∈ B, t.vStart ≠ curV) :
    FrO U B s (fePState curV lv d o s) ∧
    ((hasVert = false → ¬ U (vertItem curV)) →
      FrO U B s ((finishEdge curV d o (pre ++ B).length hasVert).run s).2) := by
  have hlv : o.cls.lowval d = lv := by rw [ho]; rfl
  have hge : ¬ (lv ≥ d) := by omega
  have hr : (finishEdge curV d o (pre ++ B).length hasVert).run s =
      if o.cls.isTree = true then
        (finishTree curV d o (pre ++ B).length hasVert s.stackDir[d]!).run (feS₀ d o s)
      else (finishBack curV d o hasVert).run (feS₀ d o s) := by
    rw [finishEdge_eq]
    simp only [finishEdge', hlv, hge, ↓reduceIte, WalkM.run_bind, WalkM.get_run, run_stackDir, run_makeVs,
      run_modifyItem]
    by_cases ht : o.cls.isTree = true <;> simp only [ht, ↓reduceIte, Bool.false_eq_true] <;> rfl
  rw [hr]
  have hlenB : B.length ≤ s.tstack.length := h.out.length
  have hfr₀ : Fr U s { s with items := s.items.modify (edgeItem s.g o.e) fun it =>
      { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) } } :=
    h.fr.modify_keep _ _ (fun _ => rfl) (fun _ => rfl)
  have h₀ : FrO U B s (feS₀ d o s) := ⟨hfr₀.of_items rfl, h.out.of_tstack rfl⟩
  by_cases ht : o.cls.isTree = true
  · rw [if_pos ht]
    simp only [fePState, ht, ↓reduceIte]
    have h₁ : FrO U B s (feS₁ d o s) ∨ (feS₁ d o s).tstack.length ≤ B.length := h₀.ears hU hq
    obtain ⟨c, mid, py, vy, hts₂, -⟩ := hE.loops ht (by rw [hlv]; exact hlow)
    have hlen₂' : (pre ++ B).length + 3 ≤ (feS₂ d o s).tstack.length := by
      rw [hts₂]; simp only [List.length_cons, List.length_append]; omega
    have hlen₂ : B.length + 3 ≤ (feS₂ d o s).tstack.length := by
      have := hlen₂'; rw [List.length_append] at this; omega
    have hlen₁ : (feS₂ d o s).tstack.length ≤ (feS₁ d o s).tstack.length := length_mergeLate d
    have h₁' : FrO U B s (feS₁ d o s) := h₁.resolve_right (by omega)
    have h₂ : FrO U B s (feS₂ d o s) := (h₁'.late d).resolve_right fun hk => by
      have hk' : (feS₂ d o s).tstack.length ≤ B.length := hk
      omega
    have h₃ := h₂.finP hU (curV := curV) (lv := lv) (ty1 := o.cls.isType1) (by omega) hB
    have htree : ((finishTree curV d o (pre ++ B).length hasVert s.stackDir[d]!).run (feS₀ d o s)).2 =
        if hasVert = true then
          ((closeVert' curV s.stackDir[d]! o.cls.isType1 (pre ++ B).length (feSingle d o s) >>=
            finishRest curV d lv o.cls.isType1 hasVert).run (feS₂ d o s)).2
        else ((finishRest curV d lv o.cls.isType1 hasVert (feSingle d o s)).run (feS₂ d o s)).2 := by
      simp only [finishTree, WalkM.run_bind, closeVert_eq, hlv]
      split <;> rfl
    refine ⟨h₃.1, fun hvi => ?_⟩
    rw [htree]
    cases hasVert
    · simp only [Bool.false_eq_true, ↓reduceIte, finishRest, WalkM.run_bind]
      exact h₃.1.tail h₃.2 _ (hvi rfl)
    · simp only [↓reduceIte, finishRest, WalkM.run_bind]
      have h₄ := h₂.closeV hU (orig := (pre ++ B).length)
        (by simp only [List.length_append]; omega) hlen₂' curV s.stackDir[d]! o.cls.isType1 (feSingle d o s)
      have h₅ := h₄.1.finP hU (curV := curV) (lv := lv) (ty1 := o.cls.isType1) h₄.2 hB
      exact h₅.1
  · rw [if_neg ht]
    have ht' : o.cls.isTree = false := Bool.eq_false_iff.2 ht
    simp only [fePState, ht', Bool.false_eq_true, ↓reduceIte]
    have hb : FrO U B s (feBack curV lv d o s) :=
      (h₀.cons (t := ⟨curV, lv, (feS₀ d o s).nxtEdgeIdx, setSides (feS₀ d o s).stackDir[lv]! [edgeItem (feS₀ d o s).g o.e] []⟩)
        fun i hi => by rw [mem_setSides, List.mem_singleton] at hi; exact hi ▸ hq).of_eq rfl rfl
    have hlb : (feBack curV lv d o s).tstack.length = s.tstack.length + 1 := rfl
    have h₃ := hb.finP hU (curV := curV) (lv := lv) (ty1 := o.cls.isType1) (by omega) hB
    refine ⟨h₃.1, fun hvi => ?_⟩
    simp only [finishBack, WalkM.run_bind, hlv, finishRest]
    cases hasVert
    · exact h₃.1.tail h₃.2 _ (hvi rfl)
    · exact h₃.1



theorem spansCount_append_eq (A B : List TEntry) (i : ItemId) :
    spansCount (A ++ B) i = spansCount A i + spansCount B i := by
  induction A with
  | nil => simp [spansCount_nil]
  | cons t A ih => rw [List.cons_append, spansCount_cons, spansCount_cons, ih]; omega

theorem below_cases {items : Items} {x y : ItemId} (h : Items.Below items x y) :
    y = x ∨ ∃ p, Items.Below items x p ∧ Items.IsParent items p y :=
  Relation.ReflTransGen.cases_tail h

/-- The frame of a returning `finishEdge`: `V curV` (when `hasVert = false`) stays a root off the
stack, and every item below the untouched base `B` keeps its type and children. -/
theorem finishRet_frame {g : Graph} {P X : ItemId → Prop} {blocks : List StBlock} {d lv : Nat}
    {kind : RetKind} {o : DfsOut} {s : WalkState} {curV : Nat} {hasVert : Bool} {sub pre B : List TEntry}
    (hE : s.EarFinish curV d o hasVert sub (pre ++ B)) (hfull : s.Full g P X) (hI : StItems g s blocks)
    (hv : curV < s.g.nv) (ho : o.cls = .ret lv kind) (hlow : lv < d) (hB : ∀ t ∈ B, t.vStart ≠ curV) :
    (hasVert = false →
      (∀ p, ¬ Items.IsParent (fePState curV lv d o s).items p (vertItem curV)) ∧
      vertItem curV ∉ readStack (fePState curV lv d o s).tstack) ∧
    ∀ x ∈ readStack B, ∀ y, Items.Below s.items x y →
      Items.type ((finishEdge curV d o (pre ++ B).length hasVert).run s).2.items y = Items.type s.items y ∧
      Items.ch ((finishEdge curV d o (pre ++ B).length hasVert).run s).2.items y = Items.ch s.items y := by
  have hts : s.tstack = (sub ++ pre) ++ B := by rw [hE.tstack, List.append_assoc]
  have hBsub : ∀ t ∈ B, t ∈ s.tstack := fun t ht => by rw [hts]; exact List.mem_append_right _ ht
  have hAsub : ∀ t ∈ sub ++ pre, t ∈ s.tstack := fun t ht => by rw [hts]; exact List.mem_append_left _ ht
  have hmemB : ∀ x ∈ readStack B, x ∈ readStack s.tstack := fun x hx => by
    rw [hts]; exact mem_readStack_append.2 (Or.inr hx)
  have hmemA : ∀ x ∈ readStack (sub ++ pre), x ∈ readStack s.tstack := fun x hx => by
    rw [hts]; exact mem_readStack_append.2 (Or.inl hx)
  have hle := hfull.place.le
  have hdisj : ∀ x ∈ readStack (sub ++ pre), x ∉ readStack B := fun x hA hB' => by
    have h1 := spansCount_pos_of_mem_readStack hA
    have h2 := spansCount_pos_of_mem_readStack hB'
    have h3 := hle x
    simp only [WalkState.cnt, hts, spansCount_append_eq] at h1 h3
    omega
  have huniq : ∀ c p p', Items.IsParent s.items p c → Items.IsParent s.items p' c → p = p' := by
    intro c p p' hp hp'
    by_contra hne
    have h2 := Items.two_le_chCount _ (lt_of_mem_ch hp) (lt_of_mem_ch hp') hne hp hp'
    have h1 := hle c
    simp only [WalkState.cnt] at h1
    omega
  have hroot0 : s.cnt 0 = 0 := hfull.place.root
  have h0np : ∀ p, ¬ Items.IsParent s.items p 0 := noParent_of_cnt_eq_zero hroot0
  have h0st : (0 : ItemId) ∉ readStack s.tstack := fun h => by
    have := spansCount_pos_of_mem_readStack h
    simp only [WalkState.cnt] at hroot0
    omega
  have hvU : hasVert = false → ¬ (∃ x ∈ readStack B, Items.Below s.items x (vertItem curV)) :=
    fun hhv ⟨x, hx, hb⟩ => (below_cases hb).elim
      (fun e => by
        obtain ⟨t, ht, hi⟩ := mem_readStack_exists (e ▸ hx)
        have h1 := (hE.vert_free t (hBsub t ht) hi).1
        rw [hhv] at h1; exact Bool.false_ne_true h1)
      fun ⟨p, _, hp⟩ => hE.v_root p hp
  have hmainU := FrO.finishEdge_ret (U := fun i => ∃ x ∈ readStack B, Items.Below s.items x i) (B := B)
    ⟨fun y ⟨x, hx, hb⟩ => hI.bounded x (hmemB x hx) y hb,
     fun ⟨x, hx, hb⟩ => (below_cases hb).elim (fun e => h0st (e ▸ hmemB x hx)) fun ⟨p, _, hp⟩ => h0np p hp⟩
    ⟨Fr.refl fun i hi c hc ⟨x, hx, hb⟩ => (below_cases hb).elim
        (fun e => hI.roots c (e ▸ hmemB x hx) i hc)
        (fun ⟨p, hxp, hp⟩ => hi ⟨x, hx, huniq c p i hp hc ▸ hxp⟩),
     ⟨sub ++ pre, hts, fun t ht i hi ⟨x, hx, hb⟩ => (below_cases hb).elim
        (fun e => hdisj i (mem_readStack_of_mem ht hi) (e ▸ hx))
        (fun ⟨p, _, hp⟩ => hI.roots i (hmemA i (mem_readStack_of_mem ht hi)) p hp)⟩⟩
    hE ho hlow
    (fun ⟨x, hx, hb⟩ => (below_cases hb).elim
      (fun e => by
        obtain ⟨t, ht, hi⟩ := mem_readStack_exists (e ▸ hx)
        exact hE.q_free t (hBsub t ht) hi)
      fun ⟨p, _, hp⟩ => hE.q_root p hp)
    hB
  refine ⟨fun hhv => ?_, fun x hx y hb => ⟨(hmainU.2 hvU).fr.type y ⟨x, hx, hb⟩, (hmainU.2 hvU).fr.ch y ⟨x, hx, hb⟩⟩⟩
  have hmainV := FrO.finishEdge_ret (U := fun i => i = vertItem curV) (B := B)
    ⟨fun y hy => by
        have h1 := hfull.place.size; have h2 := hfull.place.g_eq
        rw [hy]; show 1 + curV < _; rw [← h2] at h1; omega,
     fun h => by have : (0 : Nat) = 1 + curV := h; omega⟩
    ⟨Fr.refl fun i _ c hc hcv => hE.v_root i (hcv ▸ hc),
     ⟨sub ++ pre, hts, fun t ht i hi hiv => by
        have h1 := (hE.vert_free t (hAsub t ht) (hiv ▸ hi)).1
        rw [hhv] at h1; exact Bool.false_ne_true h1⟩⟩
    hE ho hlow (fun h => vertItem_ne_edgeItem (g := s.g) hv o.e h.symm) hB
  refine ⟨fun p hp => ?_, fun hmem => ?_⟩
  · by_cases hpv : p = vertItem curV
    · have h1 := hmainV.1.fr.ch (vertItem curV) rfl
      rw [hpv, Items.IsParent, h1] at hp
      exact hE.v_root _ hp
    · exact hmainV.1.fr.out p hpv _ hp rfl
  · obtain ⟨above, hab, hfree⟩ := hmainV.1.out
    rw [hab, mem_readStack_append] at hmem
    rcases hmem with hmem | hmem
    · obtain ⟨t, ht, hi⟩ := mem_readStack_exists hmem
      exact hfree t ht _ hi rfl
    · obtain ⟨t, ht, hi⟩ := mem_readStack_exists hmem
      have h1 := (hE.vert_free t (hBsub t ht) hi).1
      rw [hhv] at h1; exact Bool.false_ne_true h1

end Spqr
