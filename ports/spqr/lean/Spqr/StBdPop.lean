import Spqr.StClose
import Spqr.StSimLemmas
import Spqr.StVert
/-!
# Item-array bookkeeping of a boundary pop (`finishBoundary_st`)

At a boundary edge the block root `q = Q e` (a root not on the stack) receives as children the span
items of the popped entries `sub` plus fresh `I`/`O` items, and `q` is appended to the children of
`v = V curV` (also a root not on the stack); every other item is untouched.  `BdPop` records these
facts about the final item array; `StItems.bdPop` turns them, together with the block completion of
the popped items (`hbd`), into `StItems` of the new state.
-/
namespace Spqr

theorem ExpandsList.nil_inv {items : Items} {L : List ItemId} (h : ExpandsList items [] L) : L = [] := by
  cases h; rfl

theorem Expands.leaf_inv {items : Items} {x : ItemId} {L : List ItemId} (h : Expands items x L)
    (hx : Items.type items x = .V ∨ Items.type items x = .Q) : L = [x] := by
  cases h with
  | leaf _ h' => rw [ExpandsList.nil_inv h']
  | node h' _ => exact absurd hx h'

theorem Expands.node_inv {items : Items} {x : ItemId} {L : List ItemId} (h : Expands items x L)
    (hx : ¬ (Items.type items x = .V ∨ Items.type items x = .Q)) :
    ExpandsList items (Items.ch items x) L := by
  cases h with
  | leaf h' _ => exact absurd h' hx
  | node _ h' => simpa using h'

structure BdPop (s : WalkState) (sub : List TEntry) (q v : ItemId) (items' : Items) : Prop where
  size : s.items.size ≤ items'.size
  frame : ∀ y, y ≠ q → y ≠ v → y < s.items.size →
    Items.type items' y = Items.type s.items y ∧ Items.ch items' y = Items.ch s.items y ∧
    Items.vs items' y = Items.vs s.items y
  type_q : Items.type items' q = Items.type s.items q
  type_v : Items.type items' v = Items.type s.items v
  fresh_type : ∀ y, s.items.size ≤ y →
    ¬ (Items.type items' y = .S ∨ Items.type items' y = .P ∨ Items.type items' y = .R)
  fresh_ch : ∀ y, s.items.size ≤ y → Items.ch items' y = []
  ch_v : Items.ch items' v = Items.ch s.items v ++ [q]
  ch_q : ∀ c ∈ Items.ch items' q, c ∈ readStack sub ∨ (s.items.size ≤ c ∧ c < items'.size)
  ch_q_nodup : (Items.ch items' q).Nodup

theorem readStack_nodup_right {sub base : List TEntry} (h : (readStack (sub ++ base)).Nodup) :
    (readStack base).Nodup := by
  refine h.sublist ?_
  simp only [readStack, readL_append, readR_append]
  exact (List.sublist_append_right _ _).append (List.sublist_append_left _ _)

section
variable {g : Graph} {s s' : WalkState} {sub base : List TEntry} {q v : ItemId}
  (hI : StItems g s (blocks := blocks)) (hP : BdPop s sub q v s'.items)
  (hq_root : ∀ p, ¬ Items.IsParent s.items p q) (hv_root : ∀ p, ¬ Items.IsParent s.items p v)
include hI hP hq_root hv_root

theorem BdPop.below_old {x y : ItemId} (hxq : x ≠ q) (hxv : x ≠ v) (hx : x < s.items.size)
    (h : Items.Below s'.items x y) : Items.Below s.items x y := by
  revert hxq hxv hx
  induction h using Relation.ReflTransGen.head_induction_on with
  | refl => intros; exact .refl
  | @head a c hac _ ih =>
    intro haq hav ha
    have hac' : Items.IsParent s.items a c := by
      unfold Items.IsParent at hac ⊢; rwa [(hP.frame a haq hav ha).2.1] at hac
    have hcq : c ≠ q := fun e => hq_root a (e ▸ hac')
    have hcv : c ≠ v := fun e => hv_root a (e ▸ hac')
    exact .head hac' (ih hcq hcv (hI.chLt a c hac'))

theorem BdPop.below_new {x y : ItemId} (hxq : x ≠ q) (hxv : x ≠ v) (hx : x < s.items.size)
    (h : Items.Below s.items x y) : Items.Below s'.items x y := by
  revert hxq hxv hx
  induction h using Relation.ReflTransGen.head_induction_on with
  | refl => intros; exact .refl
  | @head a c hac _ ih =>
    intro haq hav ha
    have hac' : Items.IsParent s'.items a c := by
      unfold Items.IsParent at hac ⊢; rwa [(hP.frame a haq hav ha).2.1]
    have hcq : c ≠ q := fun e => hq_root a (e ▸ hac)
    have hcv : c ≠ v := fun e => hv_root a (e ▸ hac)
    exact .head hac' (ih hcq hcv (hI.chLt a c hac))

/-- Items below a sub/base span item are below it in `s`; `q`, `v` and fresh items are never below. -/
theorem BdPop.inBlock {b : StBlock} {i : ItemId} (hiq : i ≠ q) (hiv : i ≠ v) (hi : i < s.items.size)
    (h : InBlock g s.items b i) : InBlock g s'.items b i := by
  have hnb : ∀ y, Items.Below s.items i y → y ≠ q ∧ y ≠ v ∧ y < s.items.size := fun y hy =>
    ⟨fun e => Items.not_below_root_of_ne hq_root hiq (e ▸ hy),
     fun e => Items.not_below_root_of_ne hv_root hiv (e ▸ hy),
     Items.Below.lt_of_chLt hI.chLt hi hy⟩
  refine h.congr (fun y hy => ?_) (hP.frame i hiq hiv hi).2.2 (fun c hc => ?_)
    (fun e he => hP.below_old hI hq_root hv_root hiq hiv hi he) (fun c hc e he _ => ?_)
  · obtain ⟨h1, h2, h3⟩ := hnb y hy
    exact ⟨(hP.frame y h1 h2 h3).1, (hP.frame y h1 h2 h3).2.1⟩
  · obtain ⟨h1, h2, h3⟩ := hnb c (.single hc)
    exact (hP.frame c h1 h2 h3).2.2
  · obtain ⟨h1, h2, h3⟩ := hnb c (.single hc)
    exact hP.below_old hI hq_root hv_root h1 h2 h3 he

end

theorem StItems.bdPop {g : Graph} {s s' : WalkState} {sub base : List TEntry} {q v : ItemId}
    {blocks newBlocks : List StBlock}
    (hI : StItems g s blocks) (hts : s.tstack = sub ++ base) (hts' : s'.tstack = base)
    (hP : BdPop s sub q v s'.items)
    (hq_root : ∀ p, ¬ Items.IsParent s.items p q) (hq_free : q ∉ readStack s.tstack)
    (hv_root : ∀ p, ¬ Items.IsParent s.items p v) (hv_free : v ∉ readStack s.tstack)
    (hq_ty : Items.type s.items q = .Q) (hv_ty : Items.type s.items v = .V)
    (hq_lt : q < s.items.size) (_hv_lt : v < s.items.size)
    (hbd : ∀ i, (∃ x ∈ readStack sub, Items.Below s.items x i) →
      Items.type s.items i = .S ∨ Items.type s.items i = .P ∨ Items.type s.items i = .R →
      ∃ b ∈ blocks ++ newBlocks, InBlock g s.items b i) :
    StItems g s' (blocks ++ newBlocks) := by
  have hqv : q ≠ v := fun e => by rw [e, hv_ty] at hq_ty; cases hq_ty
  have hmem : ∀ x, x ∈ readStack s.tstack ↔ x ∈ readStack sub ∨ x ∈ readStack base := fun x => by
    rw [hts]; exact mem_readStack_append
  have hsub_ne : ∀ x ∈ readStack sub, x ≠ q ∧ x ≠ v ∧ x < s.items.size := fun x hx =>
    ⟨fun e => hq_free (e ▸ (hmem x).2 (.inl hx)), fun e => hv_free (e ▸ (hmem x).2 (.inl hx)),
     hI.bounded x ((hmem x).2 (.inl hx)) x .refl⟩
  have hbase_ne : ∀ x ∈ readStack base, x ≠ q ∧ x ≠ v ∧ x < s.items.size := fun x hx =>
    ⟨fun e => hq_free (e ▸ (hmem x).2 (.inr hx)), fun e => hv_free (e ▸ (hmem x).2 (.inr hx)),
     hI.bounded x ((hmem x).2 (.inr hx)) x .refl⟩
  -- the S/P/R items of `s'` are old, not `q`/`v`, and S/P/R in `s`
  have hspr : ∀ i, Items.type s'.items i = .S ∨ Items.type s'.items i = .P ∨ Items.type s'.items i = .R →
      i ≠ q ∧ i ≠ v ∧ i < s.items.size ∧
      (Items.type s.items i = .S ∨ Items.type s.items i = .P ∨ Items.type s.items i = .R) := by
    intro i hi
    have hlt : i < s.items.size := by
      by_contra h; exact hP.fresh_type i (Nat.le_of_not_lt h) hi
    have hiq : i ≠ q := fun e => by rw [e, hP.type_q, hq_ty] at hi; simp at hi
    have hiv : i ≠ v := fun e => by rw [e, hP.type_v, hv_ty] at hi; simp at hi
    exact ⟨hiq, hiv, hlt, by rwa [(hP.frame i hiq hiv hlt).1] at hi⟩
  -- descendants of `q` in `s'`
  have hq_below : ∀ y, Items.Below s'.items q y → y = q ∨ s.items.size ≤ y ∨
      ∃ c ∈ readStack sub, Items.Below s.items c y := by
    intro y hy
    rcases Relation.ReflTransGen.cases_head hy with e | ⟨c, hc, hcy⟩
    · exact .inl e.symm
    · rcases hP.ch_q c hc with hc' | ⟨hc', _⟩
      · obtain ⟨h1, h2, h3⟩ := hsub_ne c hc'
        exact .inr (.inr ⟨c, hc', hP.below_old hI hq_root hv_root h1 h2 h3 hcy⟩)
      · rcases Relation.ReflTransGen.cases_head hcy with e | ⟨c', hc'', _⟩
        · exact .inr (.inl (e ▸ hc'))
        · unfold Items.IsParent at hc''; rw [hP.fresh_ch c hc'] at hc''; simp at hc''
  -- descendants of `v` in `s'`
  have hv_below : ∀ y, Items.Below s'.items v y → y = v ∨ Items.Below s.items v y ∨
      Items.Below s'.items q y := by
    intro y hy
    rcases Relation.ReflTransGen.cases_head hy with e | ⟨c, hc, hcy⟩
    · exact .inl e.symm
    · unfold Items.IsParent at hc
      rw [hP.ch_v, List.mem_append, List.mem_singleton] at hc
      rcases hc with hc | rfl
      · have hcq : c ≠ q := fun e => hq_root v (e ▸ hc)
        have hcv : c ≠ v := fun e => hv_root v (e ▸ hc)
        exact .inr (.inl (.head hc (hP.below_old hI hq_root hv_root hcq hcv (hI.chLt v c hc) hcy)))
      · exact .inr (.inr hcy)
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · -- roots
    intro x hx p hp
    rw [hts'] at hx
    obtain ⟨hxq, hxv, hxlt⟩ := hbase_ne x hx
    unfold Items.IsParent at hp
    by_cases hpv : p = v
    · rw [hpv, hP.ch_v, List.mem_append, List.mem_singleton] at hp
      rcases hp with hp | hp
      · exact hI.roots x ((hmem x).2 (.inr hx)) v hp
      · exact hxq hp
    by_cases hpq : p = q
    · rw [hpq] at hp
      rcases hP.ch_q x hp with hx' | ⟨hx', _⟩
      · have hnd := hI.nodup
        rw [hts] at hnd
        simp only [readStack, readL_append, readR_append] at hnd hx hx'
        rw [List.nodup_append] at hnd
        obtain ⟨hL, hR, hLR⟩ := hnd
        rw [List.nodup_append] at hL hR
        simp only [List.mem_append] at hx hx'
        rcases hx with hx | hx <;> rcases hx' with hx' | hx'
        · exact hL.2.2 x hx' x hx rfl
        · exact hLR x (List.mem_append_right _ hx) x (List.mem_append_right _ hx') rfl
        · exact hLR x (List.mem_append_left _ hx') x (List.mem_append_left _ hx) rfl
        · exact hR.2.2 x hx x hx' rfl
      · exact absurd hxlt (Nat.not_lt.2 hx')
    by_cases hplt : p < s.items.size
    · rw [(hP.frame p hpq hpv hplt).2.1] at hp
      exact hI.roots x ((hmem x).2 (.inr hx)) p hp
    · rw [hP.fresh_ch p (Nat.le_of_not_lt hplt)] at hp; simp at hp
  · -- nodup
    rw [hts']
    have := hI.nodup
    rw [hts] at this
    exact readStack_nodup_right this
  · -- bounded
    intro x hx y hy
    rw [hts'] at hx
    obtain ⟨hxq, hxv, hxlt⟩ := hbase_ne x hx
    exact Nat.lt_of_lt_of_le
      (hI.bounded x ((hmem x).2 (.inr hx)) y (hP.below_old hI hq_root hv_root hxq hxv hxlt hy)) hP.size
  · -- chLt
    intro p c hp
    unfold Items.IsParent at hp
    by_cases hpv : p = v
    · rw [hpv, hP.ch_v, List.mem_append, List.mem_singleton] at hp
      rcases hp with hp | rfl
      · exact Nat.lt_of_lt_of_le (hI.chLt _ _ hp) hP.size
      · exact Nat.lt_of_lt_of_le hq_lt hP.size
    by_cases hpq : p = q
    · rw [hpq] at hp
      rcases hP.ch_q c hp with hc | ⟨_, hc⟩
      · exact Nat.lt_of_lt_of_le (hsub_ne c hc).2.2 hP.size
      · exact hc
    by_cases hplt : p < s.items.size
    · rw [(hP.frame p hpq hpv hplt).2.1] at hp
      exact Nat.lt_of_lt_of_le (hI.chLt _ _ hp) hP.size
    · rw [hP.fresh_ch p (Nat.le_of_not_lt hplt)] at hp; simp at hp
  · -- chNodup
    intro p
    by_cases hpv : p = v
    · rw [hpv, hP.ch_v, List.nodup_append]
      refine ⟨hI.chNodup v, List.nodup_singleton q, fun a ha b hb e => ?_⟩
      rw [List.mem_singleton] at hb
      subst hb; subst e
      exact hq_root v ha
    by_cases hpq : p = q
    · rw [hpq]; exact hP.ch_q_nodup
    by_cases hplt : p < s.items.size
    · rw [(hP.frame p hpq hpv hplt).2.1]; exact hI.chNodup p
    · rw [hP.fresh_ch p (Nat.le_of_not_lt hplt)]; exact List.nodup_nil
  · -- closed
    intro i _ hty
    obtain ⟨hiq, hiv, hilt, hty'⟩ := hspr i hty
    rcases hI.closed i hilt hty' with ⟨x, hx, hxi⟩ | ⟨b, hb, hbi⟩
    · rcases (hmem x).1 hx with hx' | hx'
      · obtain ⟨b, hb, hbi⟩ := hbd i ⟨x, hx', hxi⟩ hty'
        exact .inr ⟨b, hb, hP.inBlock hI hq_root hv_root hiq hiv hilt hbi⟩
      · obtain ⟨hxq, hxv, hxlt⟩ := hbase_ne x hx'
        exact .inl ⟨x, hts' ▸ hx', hP.below_new hI hq_root hv_root hxq hxv hxlt hxi⟩
    · exact .inr ⟨b, List.mem_append_left _ hb, hP.inBlock hI hq_root hv_root hiq hiv hilt hbi⟩
  · -- finished
    intro x i hx hxi hty
    obtain ⟨hiq, hiv, hilt, hty'⟩ := hspr i hty
    have fromQ : Items.Below s'.items q i → ∃ b ∈ blocks ++ newBlocks, InBlock g s'.items b i := by
      intro h
      rcases hq_below i h with e | h | ⟨c, hc, hci⟩
      · exact absurd e hiq
      · exact absurd hilt (Nat.not_lt.2 h)
      · obtain ⟨b, hb, hbi⟩ := hbd i ⟨c, hc, hci⟩ hty'
        exact ⟨b, hb, hP.inBlock hI hq_root hv_root hiq hiv hilt hbi⟩
    by_cases hxv : x = v
    · rw [hxv] at hxi
      rcases hv_below i hxi with e | h | h
      · exact absurd e hiv
      · obtain ⟨b, hb, hbi⟩ := hI.finished v i (.inl hv_ty) h hty'
        exact ⟨b, List.mem_append_left _ hb, hP.inBlock hI hq_root hv_root hiq hiv hilt hbi⟩
      · exact fromQ h
    by_cases hxq : x = q
    · rw [hxq] at hxi; exact fromQ hxi
    by_cases hxlt : x < s.items.size
    · rw [(hP.frame x hxq hxv hxlt).1] at hx
      obtain ⟨b, hb, hbi⟩ := hI.finished x i hx (hP.below_old hI hq_root hv_root hxq hxv hxlt hxi) hty'
      exact ⟨b, List.mem_append_left _ hb, hP.inBlock hI hq_root hv_root hiq hiv hilt hbi⟩
    · rcases Relation.ReflTransGen.cases_head hxi with e | ⟨c, hc, _⟩
      · exact absurd hilt (e ▸ hxlt)
      · unfold Items.IsParent at hc; rw [hP.fresh_ch x (Nat.le_of_not_lt hxlt)] at hc; simp at hc

/-- A member of an expanded list owns a contiguous segment of the expansion. -/
theorem ExpandsList.segment_of_mem {items : Items} {xs M : List ItemId} (h : ExpandsList items xs M)
    {x : ItemId} (hx : x ∈ xs) : ∃ L A B, Expands items x L ∧ M = A ++ L ++ B := by
  obtain ⟨a, b, rfl⟩ := List.append_of_mem hx
  obtain ⟨A, M', rfl, _, h₂⟩ := ExpandsList.split h
  obtain ⟨L, B, rfl, hL, _⟩ := ExpandsList.cons_iff.1 h₂
  exact ⟨L, A, B, hL, (List.append_assoc _ _ _).symm⟩

/-- An item below a member `x` of an expanded list either lies below a leaf-type item (which
`Expands` does not enter) or owns a segment of the expansion. -/
theorem Items.Below.expands_cases {items : Items} {x i : ItemId} (hb : Items.Below items x i) :
    ∀ {xs M : List ItemId}, ExpandsList items xs M → x ∈ xs →
    (∃ y, y ≠ i ∧ (Items.type items y = .V ∨ Items.type items y = .Q) ∧ Items.Below items y i) ∨
    (∃ L, Expands items i L ∧ ∃ A B, M = A ++ L ++ B) := by
  induction hb using Relation.ReflTransGen.head_induction_on with
  | refl =>
    intro xs M h hx
    obtain ⟨L, A, B, hL, hM⟩ := h.segment_of_mem hx
    exact .inr ⟨L, hL, A, B, hM⟩
  | @head a c hac hcb ih =>
    intro xs M h ha
    obtain ⟨L, A, B, hL, hM⟩ := h.segment_of_mem ha
    by_cases hai : a = i
    · subst hai; exact .inr ⟨L, hL, A, B, hM⟩
    by_cases hleaf : Items.type items a = .V ∨ Items.type items a = .Q
    · exact .inl ⟨a, hai, hleaf, .head hac hcb⟩
    rcases ih (hL.node_inv hleaf) (show c ∈ _ from hac) with h' | ⟨L', hL', A', B', hLeq⟩
    · exact .inl h'
    · refine .inr ⟨L', hL', A ++ A', B' ++ B, ?_⟩
      rw [hM, hLeq]; simp only [List.append_assoc]

/-- Block completion of the popped entries: every S / P / R item below a span item of `sub` is
`InBlock` of the new block `B` (its leaves are a segment of `stNest ps`, the orientation clauses are
`hor`) or, when the path to it passes a leaf-type item, of an old block (`StItems.finished`). -/
theorem StRead.complete {g : Graph} {items : Items} {sub : List TEntry} {ps : List StPiece}
    {blocks : List StBlock} {B : StBlock} (hR : StRead items sub ps) (hB : B.items = stNest ps)
    (hfin : ∀ x i, (Items.type items x = .V ∨ Items.type items x = .Q) → Items.Below items x i →
      Items.type items i = .S ∨ Items.type items i = .P ∨ Items.type items i = .R →
      ∃ b ∈ blocks, InBlock g items b i)
    (hor : ∀ i L, (∃ x ∈ readStack sub, Items.Below items x i) →
      Items.type items i = .S ∨ Items.type items i = .P ∨ Items.type items i = .R →
      Expands items i L → VsOrientedAt g items B i L)
    {i : ItemId} (hx : ∃ x ∈ readStack sub, Items.Below items x i)
    (hty : Items.type items i = .S ∨ Items.type items i = .P ∨ Items.type items i = .R) :
    ∃ b ∈ blocks ++ [B], InBlock g items b i := by
  obtain ⟨x, hxs, hxi⟩ := hx
  have key : ∀ {xs M : List ItemId}, ExpandsList items xs M → x ∈ xs →
      (∃ A' B', stNest ps = A' ++ M ++ B') → ∃ b ∈ blocks ++ [B], InBlock g items b i := by
    intro xs M h hxm hM
    obtain ⟨A', B', hM⟩ := hM
    rcases hxi.expands_cases h hxm with ⟨y, _, hy, hyi⟩ | ⟨L, hL, A, Bb, hLeq⟩
    · obtain ⟨b, hb, hbi⟩ := hfin y i hy hyi hty
      exact ⟨b, List.mem_append_left _ hb, hbi⟩
    · refine ⟨B, List.mem_append_right _ (List.mem_singleton_self _), L, hL, ⟨A' ++ A, Bb ++ B', ?_⟩,
        hor i L ⟨x, hxs, hxi⟩ hty hL⟩
      rw [hB, hM, hLeq]; simp only [List.append_assoc]
  rcases List.mem_append.1 (show x ∈ readL sub ++ readR sub from hxs) with hxs' | hxs'
  · exact key hR.1 hxs' ⟨[], stNestR ps, by simp [stNest]⟩
  · exact key hR.2 hxs' ⟨stNestL ps, [], by simp [stNest]⟩

theorem readL_sublist_of_mem {ts : List TEntry} {t : TEntry} (h : t ∈ ts) : t.spans.1.Sublist (readL ts) := by
  induction ts with
  | nil => exact absurd h List.not_mem_nil
  | cons u ts ih =>
    rcases List.mem_cons.1 h with rfl | h
    · exact List.sublist_append_left _ _
    · exact (ih h).trans (List.sublist_append_right _ _)

theorem readR_sublist_of_mem {ts : List TEntry} {t : TEntry} (h : t ∈ ts) : t.spans.2.Sublist (readR ts) := by
  induction ts with
  | nil => exact absurd h List.not_mem_nil
  | cons u ts ih =>
    rcases List.mem_cons.1 h with rfl | h
    · exact List.sublist_append_right _ _
    · exact (ih h).trans (List.sublist_append_left _ _)

theorem readStack_nodup_mix {ts : List TEntry} (h : (readStack ts).Nodup) {t₁ t₂ : TEntry}
    (h₁ : t₁ ∈ ts) (h₂ : t₂ ∈ ts) : (t₁.spans.1 ++ t₂.spans.2).Nodup := by
  obtain ⟨hL, hR, hLR⟩ := List.nodup_append.1 (show (readL ts ++ readR ts).Nodup from h)
  exact List.nodup_append.2 ⟨hL.sublist (readL_sublist_of_mem h₁), hR.sublist (readR_sublist_of_mem h₂),
    fun a ha b hb => hLR a ((readL_sublist_of_mem h₁).subset ha) b ((readR_sublist_of_mem h₂).subset hb)⟩

/-- `BdPop` for the item array `mid` (the pre-state with `q`'s `vs` set and fresh items pushed)
followed by the two child-list updates of `q` and `v`. -/
theorem BdPop.mk' {s : WalkState} {sub : List TEntry} {q v : ItemId} {mid : Items} {C : List ItemId}
    (hqv : q ≠ v) (hq_lt : q < s.items.size) (hv_lt : v < s.items.size)
    (hsize : s.items.size ≤ mid.size)
    (hmid : ∀ y, y < s.items.size → Items.type mid y = Items.type s.items y ∧
      Items.ch mid y = Items.ch s.items y ∧ (y ≠ q → Items.vs mid y = Items.vs s.items y))
    (hfresh : ∀ y, s.items.size ≤ y →
      ¬ (Items.type mid y = .S ∨ Items.type mid y = .P ∨ Items.type mid y = .R) ∧ Items.ch mid y = [])
    (hC : ∀ c ∈ C, c ∈ readStack sub ∨ (s.items.size ≤ c ∧ c < mid.size)) (hCnd : C.Nodup) :
    BdPop s sub q v ((mid.modify q fun it => { it with ch := C }).modify v
      fun it => { it with ch := it.ch ++ [q] }) := by
  have hty : ∀ y, Items.type ((mid.modify q fun it => { it with ch := C }).modify v
      fun it => { it with ch := it.ch ++ [q] }) y = Items.type mid y := fun y => by
    rw [Items.type_modify _ _ _ (fun it => { it with ch := it.ch ++ [q] }) (fun _ => rfl),
      Items.type_modify _ _ _ (fun it => { it with ch := C }) (fun _ => rfl)]
  have hvs : ∀ y, Items.vs ((mid.modify q fun it => { it with ch := C }).modify v
      fun it => { it with ch := it.ch ++ [q] }) y = Items.vs mid y := fun y => by
    rw [Items.vs_modify_of_vs _ _ _ (fun it => { it with ch := it.ch ++ [q] }) (fun _ => rfl),
      Items.vs_modify_of_vs _ _ _ (fun it => { it with ch := C }) (fun _ => rfl)]
  have hch : ∀ y, y ≠ q → y ≠ v → Items.ch ((mid.modify q fun it => { it with ch := C }).modify v
      fun it => { it with ch := it.ch ++ [q] }) y = Items.ch mid y := fun y hq hv => by
    rw [Items.ch_modify_ne _ _ _ _ hv.symm, Items.ch_modify_ne _ _ _ _ hq.symm]
  have hq_lt' : q < mid.size := Nat.lt_of_lt_of_le hq_lt hsize
  have hv_lt' : v < (mid.modify q fun it => { it with ch := C }).size := by
    rw [Array.size_modify]; exact Nat.lt_of_lt_of_le hv_lt hsize
  have hchq : Items.ch ((mid.modify q fun it => { it with ch := C }).modify v
      fun it => { it with ch := it.ch ++ [q] }) q = C := by
    rw [Items.ch_modify_ne _ _ _ _ hqv.symm, Items.ch_modify_self _ _ _ hq_lt']
  refine ⟨by simp only [Array.size_modify]; exact hsize, fun y hyq hyv hy => ?_, ?_, ?_, fun y hy => ?_,
    fun y hy => ?_, ?_, fun c hc => ?_, ?_⟩
  · rw [hty, hch y hyq hyv, hvs]
    exact ⟨(hmid y hy).1, (hmid y hy).2.1, (hmid y hy).2.2 hyq⟩
  · rw [hty]; exact (hmid q hq_lt).1
  · rw [hty]; exact (hmid v hv_lt).1
  · rw [hty]; exact (hfresh y hy).1
  · have h1 : y ≠ q := Nat.ne_of_gt (Nat.lt_of_lt_of_le hq_lt hy)
    have h2 : y ≠ v := Nat.ne_of_gt (Nat.lt_of_lt_of_le hv_lt hy)
    rw [hch y h1 h2]; exact (hfresh y hy).2
  · rw [Items.ch_modify_self _ _ _ hv_lt']
    show (mid.modify q fun it => { it with ch := C })[v].ch ++ [q] = _
    rw [← Items.ch_eq_getElem hv_lt', Items.ch_modify_ne _ _ _ _ hqv, (hmid v hv_lt).2.1]
  · rw [hchq] at hc
    rcases hC c hc with h | ⟨h1, h2⟩
    · exact .inl h
    · exact .inr ⟨h1, by simpa only [Array.size_modify] using h2⟩
  · rw [hchq]; exact hCnd

end Spqr
