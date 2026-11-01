import Spqr.Proofs.CompCount

/-!
# Component counts under disjoint union and vertex identification

`ccCount` (number of `EdgesConn` classes containing an edge) is invariant under edge lists with
the same connectivity/incidence, adds under disjoint union (`shiftEdges`), and drops by one when
two non-isolated vertices in different components are identified (`identEdges`).
-/

namespace Spqr

open Classical

/-! ### Congruence -/

theorem ccCount_congr {es es' : List (Nat × Nat)} (n : Nat)
    (h1 : ∀ v, HasEdge es v ↔ HasEdge es' v)
    (h2 : ∀ x y, EdgesConn es x y ↔ EdgesConn es' x y) :
    ccCount es n = ccCount es' n := by
  have e1 : HasEdge es = HasEdge es' := funext fun v => propext (h1 v)
  have e2 : EdgesConn es = EdgesConn es' := funext fun x => funext fun y => propext (h2 x y)
  unfold ccCount
  rw [e1, e2]

theorem edgesConn_lift {es es' : List (Nat × Nat)} (f : Nat → Nat)
    (h : ∀ p ∈ es, (f p.1, f p.2) ∈ es') {x y : Nat} (hc : EdgesConn es x y) :
    EdgesConn es' (f x) (f y) := by
  induction hc with
  | refl => exact Relation.ReflTransGen.refl
  | tail _ hbc ih => exact Relation.ReflTransGen.tail ih (hbc.imp (h _) (h _))

theorem edgesConn_mono {es es' : List (Nat × Nat)} (h : ∀ p ∈ es, p ∈ es') {x y : Nat}
    (hc : EdgesConn es x y) : EdgesConn es' x y := by
  induction hc with
  | refl => exact Relation.ReflTransGen.refl
  | tail _ hbc ih => exact Relation.ReflTransGen.tail ih (hbc.imp (h _) (h _))

/-- Adding edges between already-connected vertices does not change connectivity. -/
theorem edgesConn_congr_of_conn {es es' : List (Nat × Nat)} (hsub : ∀ p ∈ es, p ∈ es')
    (hext : ∀ p ∈ es', p ∈ es ∨ EdgesConn es p.1 p.2) (x y : Nat) :
    EdgesConn es' x y ↔ EdgesConn es x y := by
  refine ⟨fun h => ?_, edgesConn_mono hsub⟩
  induction h with
  | refl => exact Relation.ReflTransGen.refl
  | @tail b c _ hbc ih =>
    refine Relation.ReflTransGen.trans ih ?_
    rcases hbc with h | h
    · rcases hext _ h with h | h
      · exact Relation.ReflTransGen.single (Or.inl h)
      · exact h
    · rcases hext _ h with h | h
      · exact Relation.ReflTransGen.single (Or.inr h)
      · exact edgesConn_symm h

theorem ccCount_congr_of_conn {es es' : List (Nat × Nat)} (n : Nat) (hsub : ∀ p ∈ es, p ∈ es')
    (hext : ∀ p ∈ es', p ∈ es ∨ (p.1 ≠ p.2 ∧ EdgesConn es p.1 p.2)) :
    ccCount es n = ccCount es' n := by
  refine ccCount_congr n (fun v => ⟨fun ⟨p, hp, hpv⟩ => ⟨p, hsub _ hp, hpv⟩, ?_⟩)
    (fun x y => (edgesConn_congr_of_conn hsub (fun p hp => (hext p hp).imp id And.right) x y).symm)
  rintro ⟨p, hp, hpv⟩
  rcases hext p hp with h | ⟨hne, hc⟩
  · exact ⟨p, h, hpv⟩
  · obtain ⟨h1, h2⟩ := HasEdge.of_edgesConn hc hne
    rcases hpv with rfl | rfl
    · exact h1
    · exact h2

/-! ### Disjoint union -/

def shiftEdges (n₁ : Nat) (es : List (Nat × Nat)) : List (Nat × Nat) :=
  es.map fun p => (n₁ + p.1, n₁ + p.2)

theorem mem_shiftEdges {n₁ : Nat} {es : List (Nat × Nat)} {p : Nat × Nat} :
    p ∈ shiftEdges n₁ es ↔ ∃ q ∈ es, p = (n₁ + q.1, n₁ + q.2) := by
  simp only [shiftEdges, List.mem_map]
  exact ⟨fun ⟨q, hq, h⟩ => ⟨q, hq, h.symm⟩, fun ⟨q, hq, h⟩ => ⟨q, hq, h.symm⟩⟩

section Union

variable {es₁ es₂ : List (Nat × Nat)} {n₁ : Nat} (h₁ : ∀ p ∈ es₁, p.1 < n₁ ∧ p.2 < n₁)
include h₁

theorem edgesConn_union_left {x y : Nat} (hx : x < n₁)
    (h : EdgesConn (es₁ ++ shiftEdges n₁ es₂) x y) : y < n₁ ∧ EdgesConn es₁ x y := by
  induction h with
  | refl => exact ⟨hx, Relation.ReflTransGen.refl⟩
  | @tail b c _ hbc ih =>
    obtain ⟨hb, hxb⟩ := ih
    rcases hbc with h | h <;> rw [List.mem_append, mem_shiftEdges] at h
    · rcases h with h | ⟨q, -, hq⟩
      · exact ⟨(h₁ _ h).2, Relation.ReflTransGen.tail hxb (Or.inl h)⟩
      · simp only [Prod.mk.injEq] at hq; omega
    · rcases h with h | ⟨q, -, hq⟩
      · exact ⟨(h₁ _ h).1, Relation.ReflTransGen.tail hxb (Or.inr h)⟩
      · simp only [Prod.mk.injEq] at hq; omega

theorem edgesConn_union_right {x y : Nat}
    (h : EdgesConn (es₁ ++ shiftEdges n₁ es₂) (n₁ + x) y) :
    ∃ y', y = n₁ + y' ∧ EdgesConn es₂ x y' := by
  induction h with
  | refl => exact ⟨x, rfl, Relation.ReflTransGen.refl⟩
  | @tail b c _ hbc ih =>
    obtain ⟨b', rfl, hxb⟩ := ih
    rcases hbc with h | h <;> rw [List.mem_append, mem_shiftEdges] at h
    · rcases h with h | ⟨q, hq, hq'⟩
      · have : n₁ + b' < n₁ := (h₁ _ h).1; omega
      · simp only [Prod.mk.injEq] at hq'
        obtain ⟨h1, h2⟩ := hq'
        refine ⟨q.2, h2, Relation.ReflTransGen.tail hxb (Or.inl ?_)⟩
        have : b' = q.1 := by omega
        rw [this]; exact hq
    · rcases h with h | ⟨q, hq, hq'⟩
      · have : n₁ + b' < n₁ := (h₁ _ h).2; omega
      · simp only [Prod.mk.injEq] at hq'
        obtain ⟨h1, h2⟩ := hq'
        refine ⟨q.1, h1, Relation.ReflTransGen.tail hxb (Or.inr ?_)⟩
        have : b' = q.2 := by omega
        rw [this]; exact hq

omit h₁ in
theorem edgesConn_union_of_left {x y : Nat} (h : EdgesConn es₁ x y) :
    EdgesConn (es₁ ++ shiftEdges n₁ es₂) x y :=
  edgesConn_mono (fun _ hp => List.mem_append_left _ hp) h

omit h₁ in
theorem edgesConn_union_of_right {x y : Nat} (h : EdgesConn es₂ x y) :
    EdgesConn (es₁ ++ shiftEdges n₁ es₂) (n₁ + x) (n₁ + y) :=
  edgesConn_lift (fun v => n₁ + v)
    (fun p hp => List.mem_append_right _ (mem_shiftEdges.2 ⟨p, hp, rfl⟩)) h

omit h₁ in
theorem hasEdge_union_iff {v : Nat} :
    HasEdge (es₁ ++ shiftEdges n₁ es₂) v ↔
      HasEdge es₁ v ∨ ∃ v', v = n₁ + v' ∧ HasEdge es₂ v' := by
  constructor
  · rintro ⟨p, hp, hpv⟩
    rw [List.mem_append, mem_shiftEdges] at hp
    rcases hp with hp | ⟨q, hq, rfl⟩
    · exact Or.inl ⟨p, hp, hpv⟩
    · right
      simp only at hpv
      rcases hpv with h | h
      · exact ⟨q.1, h.symm, q, hq, Or.inl rfl⟩
      · exact ⟨q.2, h.symm, q, hq, Or.inr rfl⟩
  · rintro (⟨p, hp, hpv⟩ | ⟨v', rfl, p, hp, hpv⟩)
    · exact ⟨p, List.mem_append_left _ hp, hpv⟩
    · refine ⟨(n₁ + p.1, n₁ + p.2), List.mem_append_right _ (mem_shiftEdges.2 ⟨p, hp, rfl⟩), ?_⟩
      rcases hpv with h | h
      · exact Or.inl (by simp [h])
      · exact Or.inr (by simp [h])

theorem ccCount_union (n₂ : Nat) :
    ccCount (es₁ ++ shiftEdges n₁ es₂) (n₁ + n₂) = ccCount es₁ n₁ + ccCount es₂ n₂ := by
  rw [ccCount_eq_card_min, ccCount_eq_card_min, ccCount_eq_card_min]
  have hsplit : (Finset.range (n₁ + n₂)).filter
      (fun v => HasEdge (es₁ ++ shiftEdges n₁ es₂) v ∧
        IsMin (es₁ ++ shiftEdges n₁ es₂) (n₁ + n₂) v) =
      (Finset.range n₁).filter (fun v => HasEdge es₁ v ∧ IsMin es₁ n₁ v) ∪
        ((Finset.range n₂).filter (fun v => HasEdge es₂ v ∧ IsMin es₂ n₂ v)).image
          (fun v => n₁ + v) := by
    ext v
    simp only [Finset.mem_filter, Finset.mem_range, Finset.mem_union, Finset.mem_image]
    by_cases hv : v < n₁
    · have hL : HasEdge (es₁ ++ shiftEdges n₁ es₂) v ↔ HasEdge es₁ v := by
        rw [hasEdge_union_iff]
        exact ⟨fun h => h.elim id (fun ⟨v', hv', _⟩ => by omega), Or.inl⟩
      have hM : IsMin (es₁ ++ shiftEdges n₁ es₂) (n₁ + n₂) v ↔ IsMin es₁ n₁ v := by
        constructor
        · intro h w hw hvw
          exact h w (by omega) (edgesConn_union_of_left hvw)
        · intro h w hw hvw
          obtain ⟨hw', hvw'⟩ := edgesConn_union_left h₁ hv hvw
          exact h w hw' hvw'
      rw [hL, hM]
      constructor
      · rintro ⟨-, h⟩; exact Or.inl ⟨hv, h⟩
      · rintro (⟨-, h⟩ | ⟨a, -, ha⟩)
        · exact ⟨by omega, h⟩
        · omega
    · obtain ⟨v', rfl⟩ : ∃ v', v = n₁ + v' := ⟨v - n₁, by omega⟩
      have hL : HasEdge (es₁ ++ shiftEdges n₁ es₂) (n₁ + v') ↔ HasEdge es₂ v' := by
        rw [hasEdge_union_iff]
        constructor
        · rintro (⟨p, hp, hpv⟩ | ⟨w, hw, h⟩)
          · have := h₁ p hp; omega
          · have : w = v' := by omega
            exact this ▸ h
        · intro h; exact Or.inr ⟨v', rfl, h⟩
      have hM : IsMin (es₁ ++ shiftEdges n₁ es₂) (n₁ + n₂) (n₁ + v') ↔ IsMin es₂ n₂ v' := by
        constructor
        · intro h w hw hvw
          have := h (n₁ + w) (by omega) (edgesConn_union_of_right hvw); omega
        · intro h w hw hvw
          obtain ⟨w', rfl, hvw'⟩ := edgesConn_union_right h₁ hvw
          have := h w' (by omega) hvw'; omega
      rw [hL, hM]
      constructor
      · rintro ⟨hlt, h⟩; exact Or.inr ⟨v', ⟨by omega, h⟩, rfl⟩
      · rintro (⟨h, -⟩ | ⟨a, ⟨ha, h⟩, hav⟩)
        · omega
        · have : a = v' := by omega
          subst this; exact ⟨by omega, h⟩
  rw [hsplit, Finset.card_union_of_disjoint,
    Finset.card_image_of_injective _ (fun a b h => Nat.add_left_cancel h)]
  rw [Finset.disjoint_left]
  intro v hv hv'
  rw [Finset.mem_filter, Finset.mem_range] at hv
  rw [Finset.mem_image] at hv'
  obtain ⟨w, -, rfl⟩ := hv'
  omega

end Union

/-! ### Identifying vertex `b` with `a` -/

def ident (a b x : Nat) : Nat := if x = b then a else x

def identEdges (a b : Nat) (es : List (Nat × Nat)) : List (Nat × Nat) :=
  es.map fun p => (ident a b p.1, ident a b p.2)

theorem mem_identEdges {a b : Nat} {es : List (Nat × Nat)} {p : Nat × Nat} :
    p ∈ identEdges a b es ↔ ∃ q ∈ es, p = (ident a b q.1, ident a b q.2) := by
  simp only [identEdges, List.mem_map]
  exact ⟨fun ⟨q, hq, h⟩ => ⟨q, hq, h.symm⟩, fun ⟨q, hq, h⟩ => ⟨q, hq, h.symm⟩⟩

section Ident

variable {a b : Nat} (hab : a ≠ b)
include hab

theorem ident_ne_b (x : Nat) : ident a b x ≠ b := by
  unfold ident; split_ifs <;> omega

omit hab in
theorem ident_of_ne {x : Nat} (h : x ≠ b) : ident a b x = x := by simp [ident, h]

omit hab in
theorem ident_b : ident a b b = a := by simp [ident]

theorem ident_eq_iff {x y : Nat} :
    ident a b x = ident a b y ↔ x = y ∨ (x = a ∧ y = b) ∨ (x = b ∧ y = a) := by
  unfold ident; split_ifs <;> omega

omit hab in
theorem edgesConn_ident_of {es : List (Nat × Nat)} {x y : Nat} (h : EdgesConn es x y) :
    EdgesConn (identEdges a b es) (ident a b x) (ident a b y) :=
  edgesConn_lift (ident a b) (fun p hp => mem_identEdges.2 ⟨p, hp, rfl⟩) h

theorem edgesConn_ident_iff {es : List (Nat × Nat)} (x y : Nat) :
    EdgesConn (identEdges a b es) (ident a b x) (ident a b y) ↔
      EdgesConn es x y ∨ (EdgesConn es x a ∧ EdgesConn es b y) ∨
        (EdgesConn es x b ∧ EdgesConn es a y) := by
  constructor
  · intro h
    suffices ∀ p q, EdgesConn (identEdges a b es) p q →
        ∀ x y, ident a b x = p → ident a b y = q →
        EdgesConn es x y ∨ (EdgesConn es x a ∧ EdgesConn es b y) ∨
          (EdgesConn es x b ∧ EdgesConn es a y) from
      this _ _ h x y rfl rfl
    intro p q h
    induction h with
    | refl =>
      intro x y hx hy
      rw [← hy] at hx
      rcases (ident_eq_iff hab).1 hx with rfl | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩
      · exact Or.inl Relation.ReflTransGen.refl
      · exact Or.inr (Or.inl ⟨Relation.ReflTransGen.refl, Relation.ReflTransGen.refl⟩)
      · exact Or.inr (Or.inr ⟨Relation.ReflTransGen.refl, Relation.ReflTransGen.refl⟩)
    | @tail q r _ hqr ih =>
      intro x y hx hy
      obtain ⟨s, t, hst, hs, ht⟩ : ∃ s t, ((s, t) ∈ es ∨ (t, s) ∈ es) ∧
          ident a b s = q ∧ ident a b t = r := by
        rcases hqr with h | h <;> obtain ⟨u, hu, hu'⟩ := mem_identEdges.1 h <;>
          simp only [Prod.mk.injEq] at hu'
        · exact ⟨u.1, u.2, Or.inl hu, hu'.1.symm, hu'.2.symm⟩
        · exact ⟨u.2, u.1, Or.inr hu, hu'.2.symm, hu'.1.symm⟩
      have h1 := ih x s hx hs
      have hst' : EdgesConn es s t := Relation.ReflTransGen.single hst
      have h2 : EdgesConn es x t ∨ (EdgesConn es x a ∧ EdgesConn es b t) ∨
          (EdgesConn es x b ∧ EdgesConn es a t) := by
        rcases h1 with h1 | ⟨h1, h1'⟩ | ⟨h1, h1'⟩
        · exact Or.inl (Relation.ReflTransGen.trans h1 hst')
        · exact Or.inr (Or.inl ⟨h1, Relation.ReflTransGen.trans h1' hst'⟩)
        · exact Or.inr (Or.inr ⟨h1, Relation.ReflTransGen.trans h1' hst'⟩)
      rw [← hy] at ht
      rcases (ident_eq_iff hab).1 ht with rfl | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩
      · exact h2
      · rcases h2 with h2 | ⟨h2, -⟩ | ⟨h2, -⟩
        · exact Or.inr (Or.inl ⟨h2, Relation.ReflTransGen.refl⟩)
        · exact Or.inr (Or.inl ⟨h2, Relation.ReflTransGen.refl⟩)
        · exact Or.inl h2
      · rcases h2 with h2 | ⟨h2, -⟩ | ⟨h2, -⟩
        · exact Or.inr (Or.inr ⟨h2, Relation.ReflTransGen.refl⟩)
        · exact Or.inl h2
        · exact Or.inr (Or.inr ⟨h2, Relation.ReflTransGen.refl⟩)
  · rintro (h | ⟨h1, h2⟩ | ⟨h1, h2⟩)
    · exact edgesConn_ident_of h
    · have e1 := edgesConn_ident_of (a := a) (b := b) h1
      have e2 := edgesConn_ident_of (a := a) (b := b) h2
      rw [ident_b] at e2
      rw [ident_of_ne hab] at e1
      exact Relation.ReflTransGen.trans e1 e2
    · have e1 := edgesConn_ident_of (a := a) (b := b) h1
      have e2 := edgesConn_ident_of (a := a) (b := b) h2
      rw [ident_b] at e1
      rw [ident_of_ne hab] at e2
      exact Relation.ReflTransGen.trans e1 e2

omit hab in
theorem hasEdge_ident_iff {es : List (Nat × Nat)} (v : Nat) :
    HasEdge (identEdges a b es) v ↔ ∃ s, ident a b s = v ∧ HasEdge es s := by
  constructor
  · rintro ⟨p, hp, hpv⟩
    obtain ⟨q, hq, rfl⟩ := mem_identEdges.1 hp
    simp only at hpv
    rcases hpv with h | h
    · exact ⟨q.1, h, q, hq, Or.inl rfl⟩
    · exact ⟨q.2, h, q, hq, Or.inr rfl⟩
  · rintro ⟨s, rfl, p, hp, hpv⟩
    refine ⟨(ident a b p.1, ident a b p.2), mem_identEdges.2 ⟨p, hp, rfl⟩, ?_⟩
    rcases hpv with h | h
    · exact Or.inl (by rw [h])
    · exact Or.inr (by rw [h])

theorem not_hasEdge_ident_b {es : List (Nat × Nat)} : ¬HasEdge (identEdges a b es) b := by
  intro h
  obtain ⟨s, hs, -⟩ := (hasEdge_ident_iff b).1 h
  exact ident_ne_b hab s hs

omit hab in
theorem hasEdge_ident_of_ne {es : List (Nat × Nat)} {v : Nat} (hva : v ≠ a) (hvb : v ≠ b) :
    HasEdge (identEdges a b es) v ↔ HasEdge es v := by
  rw [hasEdge_ident_iff]
  constructor
  · rintro ⟨s, hs, h⟩
    have : s = v := by unfold ident at hs; split_ifs at hs <;> omega
    exact this ▸ h
  · intro h; exact ⟨v, ident_of_ne hvb, h⟩

theorem hasEdge_ident_a {es : List (Nat × Nat)} (h : HasEdge es a) :
    HasEdge (identEdges a b es) a :=
  (hasEdge_ident_iff a).2 ⟨a, ident_of_ne hab, h⟩

theorem not_edgesConn_ident_b {es : List (Nat × Nat)} {v : Nat} (hv : v ≠ b) :
    ¬EdgesConn (identEdges a b es) v b := fun h =>
  not_hasEdge_ident_b hab (HasEdge.of_edgesConn h hv).2

theorem edgesConn_ident_a_iff {es : List (Nat × Nat)} {v : Nat} (hv : v ≠ b) :
    EdgesConn (identEdges a b es) v a ↔ EdgesConn es v a ∨ EdgesConn es v b := by
  have := edgesConn_ident_iff hab (es := es) v a
  rw [ident_of_ne hv, ident_of_ne hab] at this
  rw [this]
  constructor
  · rintro (h | ⟨h, -⟩ | ⟨h, -⟩)
    · exact Or.inl h
    · exact Or.inl h
    · exact Or.inr h
  · rintro (h | h)
    · exact Or.inl h
    · exact Or.inr (Or.inr ⟨h, Relation.ReflTransGen.refl⟩)

theorem edgesConn_ident_of_ne_iff {es : List (Nat × Nat)} {v w : Nat} (hv : v ≠ b) (hw : w ≠ b)
    (hva : ¬EdgesConn es v a) (hvb : ¬EdgesConn es v b) :
    EdgesConn (identEdges a b es) v w ↔ EdgesConn es v w := by
  have := edgesConn_ident_iff hab (es := es) v w
  rw [ident_of_ne hv, ident_of_ne hw] at this
  rw [this]
  constructor
  · rintro (h | ⟨h, -⟩ | ⟨h, -⟩)
    · exact h
    · exact absurd h hva
    · exact absurd h hvb
  · exact Or.inl

/-- Least vertices of components touching neither `a` nor `b`. -/
noncomputable def restSet (es : List (Nat × Nat)) (n a b : Nat) : Finset Nat :=
  (Finset.range n).filter fun v =>
    HasEdge es v ∧ IsMin es n v ∧ ¬EdgesConn es v a ∧ ¬EdgesConn es v b

theorem ccCount_ident_eq {es : List (Nat × Nat)} {n : Nat} (ha : a < n) (hea : HasEdge es a) :
    ccCount (identEdges a b es) n = (restSet es n a b).card + 1 := by
  rw [ccCount_eq_card_min]
  rw [← Finset.card_filter_add_card_filter_not
    (s := (Finset.range n).filter
      (fun v => HasEdge (identEdges a b es) v ∧ IsMin (identEdges a b es) n v))
    (fun v => ¬EdgesConn (identEdges a b es) v a)]
  rw [Finset.filter_filter, Finset.filter_filter]
  congr 1
  · unfold restSet
    congr 1
    apply Finset.filter_congr
    intro v hv
    rw [Finset.mem_range] at hv
    by_cases hvb : v = b
    · subst hvb
      constructor
      · rintro ⟨⟨h, -⟩, -⟩; exact absurd h (not_hasEdge_ident_b hab)
      · rintro ⟨-, -, -, h⟩; exact absurd Relation.ReflTransGen.refl h
    · rw [edgesConn_ident_a_iff hab hvb, not_or]
      constructor
      · rintro ⟨⟨he, hm⟩, hva, hvb'⟩
        have hva' : v ≠ a := fun h => hva (h ▸ Relation.ReflTransGen.refl)
        refine ⟨(hasEdge_ident_of_ne hva' hvb).1 he, fun w hw hvw => ?_, hva, hvb'⟩
        have hwb : w ≠ b := fun h => hvb' (h ▸ hvw)
        exact hm w hw ((edgesConn_ident_of_ne_iff hab hvb hwb hva hvb').2 hvw)
      · rintro ⟨he, hm, hva, hvb'⟩
        have hva' : v ≠ a := fun h => hva (h ▸ Relation.ReflTransGen.refl)
        refine ⟨⟨(hasEdge_ident_of_ne hva' hvb).2 he, fun w hw hvw => ?_⟩, hva, hvb'⟩
        by_cases hwb : w = b
        · subst hwb; exact absurd hvw (not_edgesConn_ident_b hab hvb)
        · exact hm w hw ((edgesConn_ident_of_ne_iff hab hvb hwb hva hvb').1 hvw)
  · rw [Finset.card_eq_one]
    have hne : ((Finset.range n).filter (EdgesConn (identEdges a b es) a)).Nonempty :=
      ⟨a, Finset.mem_filter.2 ⟨Finset.mem_range.2 ha, Relation.ReflTransGen.refl⟩⟩
    refine ⟨Finset.min' _ hne, ?_⟩
    obtain ⟨hmn, ham⟩ := Finset.mem_filter.1 (Finset.min'_mem _ hne)
    rw [Finset.mem_range] at hmn
    ext v
    simp only [Finset.mem_filter, Finset.mem_range, Finset.mem_singleton, not_not]
    constructor
    · rintro ⟨hv, ⟨-, hm⟩, hva⟩
      apply le_antisymm
      · exact hm _ hmn (Relation.ReflTransGen.trans hva ham)
      · exact Finset.min'_le _ _
          (Finset.mem_filter.2 ⟨Finset.mem_range.2 hv, edgesConn_symm hva⟩)
    · rintro rfl
      refine ⟨hmn, ⟨hasEdge_of_edgesConn ham (hasEdge_ident_a hab hea), fun w hw hmw => ?_⟩,
        edgesConn_symm ham⟩
      exact Finset.min'_le _ _
        (Finset.mem_filter.2 ⟨Finset.mem_range.2 hw, Relation.ReflTransGen.trans ham hmw⟩)

omit hab in
theorem ccCount_eq_rest_add {es : List (Nat × Nat)} {n : Nat} (ha : a < n) (hb : b < n)
    (hea : HasEdge es a) (heb : HasEdge es b) :
    ∃ ma mb, (ma = mb ↔ EdgesConn es a b) ∧
      ccCount es n = (restSet es n a b).card + ({ma, mb} : Finset Nat).card := by
  have hnea : ((Finset.range n).filter (EdgesConn es a)).Nonempty :=
    ⟨a, Finset.mem_filter.2 ⟨Finset.mem_range.2 ha, Relation.ReflTransGen.refl⟩⟩
  have hneb : ((Finset.range n).filter (EdgesConn es b)).Nonempty :=
    ⟨b, Finset.mem_filter.2 ⟨Finset.mem_range.2 hb, Relation.ReflTransGen.refl⟩⟩
  obtain ⟨hman, hama⟩ := Finset.mem_filter.1 (Finset.min'_mem _ hnea)
  obtain ⟨hmbn, hbmb⟩ := Finset.mem_filter.1 (Finset.min'_mem _ hneb)
  rw [Finset.mem_range] at hman hmbn
  refine ⟨Finset.min' _ hnea, Finset.min' _ hneb, ⟨fun h => ?_, fun h => ?_⟩, ?_⟩
  · exact Relation.ReflTransGen.trans (h ▸ hama) (edgesConn_symm hbmb)
  · apply le_antisymm
    · exact Finset.min'_le _ _ (Finset.mem_filter.2
        ⟨Finset.mem_range.2 hmbn, Relation.ReflTransGen.trans h hbmb⟩)
    · exact Finset.min'_le _ _ (Finset.mem_filter.2
        ⟨Finset.mem_range.2 hman, Relation.ReflTransGen.trans (edgesConn_symm h) hama⟩)
  rw [ccCount_eq_card_min]
  rw [← Finset.card_filter_add_card_filter_not
    (s := (Finset.range n).filter (fun v => HasEdge es v ∧ IsMin es n v))
    (fun v => ¬(EdgesConn es v a ∨ EdgesConn es v b))]
  rw [Finset.filter_filter, Finset.filter_filter]
  congr 1
  · unfold restSet
    congr 1
    apply Finset.filter_congr
    intro v _
    rw [not_or, and_assoc]
  · congr 1
    ext v
    simp only [Finset.mem_filter, Finset.mem_range, Finset.mem_insert, Finset.mem_singleton,
      not_not]
    constructor
    · rintro ⟨hv, ⟨-, hm⟩, hva | hvb⟩
      · left
        apply le_antisymm
        · exact hm _ hman (Relation.ReflTransGen.trans hva hama)
        · exact Finset.min'_le _ _
            (Finset.mem_filter.2 ⟨Finset.mem_range.2 hv, edgesConn_symm hva⟩)
      · right
        apply le_antisymm
        · exact hm _ hmbn (Relation.ReflTransGen.trans hvb hbmb)
        · exact Finset.min'_le _ _
            (Finset.mem_filter.2 ⟨Finset.mem_range.2 hv, edgesConn_symm hvb⟩)
    · rintro (rfl | rfl)
      · refine ⟨hman, ⟨hasEdge_of_edgesConn hama hea, fun w hw hmw => ?_⟩,
          Or.inl (edgesConn_symm hama)⟩
        exact Finset.min'_le _ _
          (Finset.mem_filter.2 ⟨Finset.mem_range.2 hw, Relation.ReflTransGen.trans hama hmw⟩)
      · refine ⟨hmbn, ⟨hasEdge_of_edgesConn hbmb heb, fun w hw hmw => ?_⟩,
          Or.inr (edgesConn_symm hbmb)⟩
        exact Finset.min'_le _ _
          (Finset.mem_filter.2 ⟨Finset.mem_range.2 hw, Relation.ReflTransGen.trans hbmb hmw⟩)

theorem ccCount_ident_of_conn {es : List (Nat × Nat)} {n : Nat} (ha : a < n) (hb : b < n)
    (hea : HasEdge es a) (heb : HasEdge es b) (hc : EdgesConn es a b) :
    ccCount (identEdges a b es) n = ccCount es n := by
  obtain ⟨ma, mb, hiff, h⟩ := ccCount_eq_rest_add ha hb hea heb
  rw [ccCount_ident_eq hab ha hea, h, hiff.2 hc,
    Finset.insert_eq_of_mem (Finset.mem_singleton_self mb), Finset.card_singleton]

theorem ccCount_ident_of_not_conn {es : List (Nat × Nat)} {n : Nat} (ha : a < n) (hb : b < n)
    (hea : HasEdge es a) (heb : HasEdge es b) (hc : ¬EdgesConn es a b) :
    ccCount (identEdges a b es) n + 1 = ccCount es n := by
  obtain ⟨ma, mb, hiff, h⟩ := ccCount_eq_rest_add ha hb hea heb
  rw [ccCount_ident_eq hab ha hea, h, Finset.card_pair (fun e => hc (hiff.1 e))]

end Ident

end Spqr
