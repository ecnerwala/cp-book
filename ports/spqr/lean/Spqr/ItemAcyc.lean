import Spqr.WalkTyping

/-!
# Acyclicity of `Items.IsParent`, as the walk maintains it

Every allocated item is below some parentless item (`Acyc`). The walk's writes to child lists are
of two shapes: appending items under a node that is below none of them (`Acyc.append`),
and clearing a node's child list when its children have no other parent (`Acyc.dropCh`,
`maybeUnwrapNxt`). At the end every non-root item has a parent, so the only parentless item is the
root (`Acyc.reach`).
-/

namespace Spqr.Items

variable {items items' : Items}

/-- `r` is in no child list. -/
def NoParent (items : Items) (r : ItemId) : Prop := ∀ p, ¬ items.IsParent p r

/-- Every allocated item is below a parentless item. -/
def Acyc (items : Items) : Prop := ∀ i, i < items.size → ∃ r, items.NoParent r ∧ items.Below r i

theorem Below.mono (h : ∀ p c, items.IsParent p c → items'.IsParent p c) {r i : ItemId}
    (hb : items.Below r i) : items'.Below r i := by
  induction hb with
  | refl => exact Relation.ReflTransGen.refl
  | tail _ hs ih => exact ih.tail (h _ _ hs)

theorem NoParent.below_eq {x c : ItemId} (hx : items.NoParent x) (h : items.Below c x) : c = x := by
  rcases Relation.ReflTransGen.cases_tail h with rfl | ⟨b, -, hb⟩
  · rfl
  · exact (hx b hb).elim

theorem below_eq_of_ch_nil {x c : ItemId} (hx : items.ch x = []) (h : items.Below x c) : c = x := by
  rcases Relation.ReflTransGen.cases_head h with rfl | ⟨b, hb, -⟩
  · rfl
  · have : b ∈ items.ch x := hb
    rw [hx] at this
    exact (List.not_mem_nil this).elim

theorem Acyc.of_ch_eq (ha : items.Acyc) (hsz : items'.size = items.size)
    (hch : ∀ j, items'.ch j = items.ch j) : items'.Acyc := by
  have hp : ∀ p c, items'.IsParent p c ↔ items.IsParent p c := fun p c => by simp [IsParent, hch]
  intro i hi
  obtain ⟨r, hr, hb⟩ := ha i (hsz ▸ hi)
  exact ⟨r, fun p hp' => hr p ((hp p r).1 hp'), Below.mono (fun p c => (hp p c).2) hb⟩

/-- A fresh item: no children, in no child list. -/
theorem Acyc.push (ha : items.Acyc) (hsz : items'.size = items.size + 1)
    (hch : ∀ j, items'.ch j = if j = items.size then [] else items.ch j)
    (hfree : items.NoParent items.size) : items'.Acyc := by
  have hp : ∀ p c, items'.IsParent p c ↔ items.IsParent p c := fun p c => by
    simp only [IsParent, hch]
    split
    · subst p; rw [ch_of_le _ _ (Nat.le_refl _)]
    · exact Iff.rfl
  intro i hi
  rw [hsz] at hi
  rcases Nat.lt_succ_iff_lt_or_eq.1 hi with hi | rfl
  · obtain ⟨r, hr, hb⟩ := ha i hi
    exact ⟨r, fun p hp' => hr p ((hp p r).1 hp'), Below.mono (fun p c => (hp p c).2) hb⟩
  · exact ⟨_, fun p hp' => hfree p ((hp p _).1 hp'), Relation.ReflTransGen.refl⟩

theorem isParent_append_iff {a : ItemId} {L : List ItemId}
    (hch : ∀ j, items'.ch j = if j = a then items.ch a ++ L else items.ch j) (p c : ItemId) :
    items'.IsParent p c ↔ items.IsParent p c ∨ (p = a ∧ c ∈ L) := by
  simp only [IsParent, hch]
  split
  · subst p; simp
  · simp [*]

/-- Along `items'` (which only adds edges out of `a`), a path from `c ≠ a` cannot pass through the
parentless `a`. -/
theorem below_of_append {a : ItemId} {L : List ItemId}
    (hch : ∀ j, items'.ch j = if j = a then items.ch a ++ L else items.ch j) (ha : items.NoParent a)
    {c j : ItemId} (hca : c ≠ a) (h : items'.Below c j) : items.Below c j := by
  induction h with
  | refl => exact Relation.ReflTransGen.refl
  | @tail j j' _ hjj' ih =>
    rcases (isParent_append_iff hch j j').1 hjj' with h | ⟨rfl, _⟩
    · exact ih.tail h
    · exact absurd (ha.below_eq ih) hca

/-- Appending items under `a`, which is below none of them. -/
theorem Acyc.append (ha : items.Acyc) {a : ItemId} {L : List ItemId} (hsz : items'.size = items.size)
    (hch : ∀ j, items'.ch j = if j = a then items.ch a ++ L else items.ch j) (halt : a < items.size)
    (hLa : ∀ c ∈ L, ¬ items.Below c a) : items'.Acyc := by
  have hp := isParent_append_iff hch
  have hmono : ∀ {r i}, items.Below r i → items'.Below r i := fun hb =>
    Below.mono (fun p c h => (hp p c).2 (Or.inl h)) hb
  intro i hi
  obtain ⟨r, hr, hb⟩ := ha i (hsz ▸ hi)
  by_cases hrL : r ∈ L
  · obtain ⟨r', hr', hb'⟩ := ha a halt
    have hr'L : r' ∉ L := fun h => hLa r' h hb'
    refine ⟨r', fun p hp' => ?_, ?_⟩
    · rcases (hp p r').1 hp' with h | ⟨_, h⟩
      · exact hr' p h
      · exact hr'L h
    · exact (hmono hb').trans (Relation.ReflTransGen.head ((hp a r).2 (Or.inr ⟨rfl, hrL⟩)) (hmono hb))
  · refine ⟨r, fun p hp' => ?_, hmono hb⟩
    rcases (hp p r).1 hp' with h | ⟨_, h⟩
    · exact hr p h
    · exact hrL h

/-- Clearing the children of `x`, none of which has another parent. -/
theorem Acyc.dropCh (ha : items.Acyc) {x : ItemId} (hsz : items'.size = items.size)
    (hch : ∀ j, items'.ch j = if j = x then [] else items.ch j)
    (huniq : ∀ p c, c ∈ items.ch x → items.IsParent p c → p = x) : items'.Acyc := by
  have hp : ∀ p c, items'.IsParent p c ↔ items.IsParent p c ∧ p ≠ x := fun p c => by
    simp only [IsParent, hch]
    split
    · simp [*]
    · simp [*]
  have split : ∀ {r i}, items.Below r i → items'.Below r i ∨ ∃ c ∈ items.ch x, items'.Below c i := by
    intro r i hb
    induction hb with
    | refl => exact Or.inl Relation.ReflTransGen.refl
    | @tail j i _ hji ih =>
      by_cases hjx : j = x
      · subst hjx; exact Or.inr ⟨i, hji, Relation.ReflTransGen.refl⟩
      · have hji' : items'.IsParent j i := (hp j i).2 ⟨hji, hjx⟩
        rcases ih with h | ⟨c, hc, h⟩
        · exact Or.inl (h.tail hji')
        · exact Or.inr ⟨c, hc, h.tail hji'⟩
  intro i hi
  obtain ⟨r, hr, hb⟩ := ha i (hsz ▸ hi)
  rcases split hb with h | ⟨c, hc, h⟩
  · exact ⟨r, fun p hp' => hr p ((hp p r).1 hp').1, h⟩
  · refine ⟨c, fun p hp' => ?_, h⟩
    obtain ⟨h1, h2⟩ := (hp p c).1 hp'
    exact h2 (huniq p c hc h1)

/-- With every non-root item parented, the root reaches everything. -/
theorem Acyc.reach (ha : items.Acyc) (hroot : ∀ r, 0 < r → r < items.size → ¬ items.NoParent r) :
    ∀ i, i < items.size → items.Below rootItem i := by
  intro i hi
  obtain ⟨r, hr, hb⟩ := ha i hi
  have hrlt : r < items.size := by
    rcases Relation.ReflTransGen.cases_head hb with rfl | ⟨c, hrc, -⟩
    · exact hi
    · refine Nat.lt_of_not_le fun h => ?_
      have : c ∈ items.ch r := hrc
      rw [ch_of_le _ _ h] at this
      exact List.not_mem_nil this
  rcases Nat.eq_zero_or_pos r with rfl | hpos
  · exact (hb : items.Below 0 i)
  · exact absurd hr (hroot r hpos hrlt)

end Spqr.Items
