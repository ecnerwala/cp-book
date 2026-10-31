import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Data.Finset.Card
import Spqr.ItemSpec

/-!
# Counting consequences of `Items.Tree`

A rooted item tree is acyclic, its ancestor chains are linear, and the descendant sets of distinct
children are disjoint; hence the descendant sets of the children of `a` together have fewer
elements than the descendant set of `a`. This is what bounds the preorder traversal of phase 3.
-/

namespace Spqr.Items

variable {items : Items} {g : Graph}

theorem getElem!_ch {j : ItemId} (h : j < items.size) : items[j]!.ch = items.ch j := by
  simp [ch, getElem!_pos, h]

theorem Tree.parent_eq (ht : items.Tree g) {p p' c : ItemId} (h : items.IsParent p c)
    (h' : items.IsParent p' c) : p = p' := by
  have hc : c < items.size := ht.ch_lt p c h
  have hc0 : 0 < c := by
    rcases Nat.eq_zero_or_pos c with rfl | hpos
    · exact absurd h (ht.root_no_parent p)
    · exact hpos
  obtain ⟨q, _, hq⟩ := ht.unique_parent c hc0 hc
  rw [hq p h, hq p' h']

theorem Tree.acyclic (ht : items.Tree g) (i : ItemId) : ¬ Relation.TransGen items.IsParent i i := by
  have key : ∀ j, items.Below rootItem j → ¬ Relation.TransGen items.IsParent j j := by
    intro j hj
    induction hj with
    | refl =>
      intro hc
      obtain ⟨b, -, hb⟩ := Relation.TransGen.tail'_iff.1 hc
      exact ht.root_no_parent b hb
    | @tail p i' _ hpi ih =>
      intro hc
      obtain ⟨b, hib, hb⟩ := Relation.TransGen.tail'_iff.1 hc
      obtain rfl := ht.parent_eq hb hpi
      exact ih (Relation.TransGen.head' hpi hib)
  intro hc
  obtain ⟨b, -, hb⟩ := Relation.TransGen.tail'_iff.1 hc
  exact key i (ht.reach i (ht.ch_lt b i hb)) hc

theorem Tree.below_linear (ht : items.Tree g) {a b x : ItemId} (hax : items.Below a x)
    (hbx : items.Below b x) : items.Below a b ∨ items.Below b a := by
  induction hax with
  | refl => exact Or.inr hbx
  | @tail x' x hax' hx'x ih =>
    rcases Relation.ReflTransGen.cases_tail hbx with rfl | ⟨c, hbc, hcx⟩
    · exact Or.inl (Relation.ReflTransGen.tail hax' hx'x)
    · obtain rfl := ht.parent_eq hcx hx'x
      exact ih hbc

open scoped Classical in
/-- The descendants-or-self of `a` (restricted to valid ids). -/
noncomputable def desc (items : Items) (a : ItemId) : Finset ItemId :=
  (Finset.range items.size).filter (items.Below a)

theorem mem_desc {a x : ItemId} : x ∈ items.desc a ↔ x < items.size ∧ items.Below a x := by
  simp [desc]

theorem card_desc_le {a : ItemId} : (items.desc a).card ≤ items.size :=
  (Finset.card_le_card fun _ hx => Finset.mem_range.2 (mem_desc.1 hx).1).trans_eq
    (Finset.card_range _)

theorem mem_desc_self {a : ItemId} (h : a < items.size) : a ∈ items.desc a :=
  mem_desc.2 ⟨h, Relation.ReflTransGen.refl⟩

theorem one_le_card_desc {a : ItemId} (h : a < items.size) : 1 ≤ (items.desc a).card :=
  Finset.card_pos.2 ⟨a, mem_desc_self h⟩

theorem Tree.desc_subset (ht : items.Tree g) {a c : ItemId} (h : items.IsParent a c) :
    items.desc c ⊆ (items.desc a).erase a := by
  intro x hx
  rw [mem_desc] at hx
  rw [Finset.mem_erase, mem_desc]
  refine ⟨?_, hx.1, Relation.ReflTransGen.head h hx.2⟩
  rintro rfl
  exact ht.acyclic x (Relation.TransGen.head' h hx.2)

theorem Tree.desc_disjoint (ht : items.Tree g) {a c₁ c₂ : ItemId} (h1 : items.IsParent a c₁)
    (h2 : items.IsParent a c₂) (hne : c₁ ≠ c₂) : Disjoint (items.desc c₁) (items.desc c₂) := by
  have aux : ∀ {c₁ c₂}, items.IsParent a c₁ → items.IsParent a c₂ → c₁ ≠ c₂ →
      items.Below c₁ c₂ → False := by
    intro c₁ c₂ h1 h2 hne h12
    rcases Relation.ReflTransGen.cases_tail h12 with rfl | ⟨p, h1p, hp2⟩
    · exact hne rfl
    · obtain rfl := ht.parent_eq hp2 h2
      exact ht.acyclic p (Relation.TransGen.head' h1 h1p)
  rw [Finset.disjoint_left]
  intro x hx1 hx2
  rw [mem_desc] at hx1 hx2
  rcases ht.below_linear hx1.2 hx2.2 with h | h
  · exact aux h1 h2 hne h
  · exact aux h2 h1 (Ne.symm hne) h

/-- The descendant sets of the children of `a` are disjoint and omit `a`. -/
theorem Tree.sum_card_desc_children (ht : items.Tree g) {a : ItemId} (ha : a < items.size) :
    ((items.ch a).map fun c => (items.desc c).card).sum + 1 ≤ (items.desc a).card := by
  classical
  have hch : ∀ c ∈ items.ch a, items.IsParent a c := fun _ hc => hc
  have hsub : ((items.ch a).toFinset.biUnion items.desc) ⊆ (items.desc a).erase a :=
    Finset.biUnion_subset.2 fun c hc => ht.desc_subset (hch c (List.mem_toFinset.1 hc))
  have hdisj : (↑(items.ch a).toFinset : Set ItemId).PairwiseDisjoint items.desc :=
    fun c₁ hc₁ c₂ hc₂ hne =>
      ht.desc_disjoint (hch c₁ (List.mem_toFinset.1 hc₁)) (hch c₂ (List.mem_toFinset.1 hc₂)) hne
  have h1 := Finset.card_le_card hsub
  rw [Finset.card_biUnion hdisj, Finset.card_erase_of_mem (mem_desc_self ha),
    List.sum_toFinset _ (ht.ch_nodup a)] at h1
  have := one_le_card_desc ha
  omega

end Spqr.Items
