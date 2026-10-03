import Spqr.Proofs.PlanarInsert
/-!
# Component / non-isolated bookkeeping of one edge (PROOF.md §8.6)
How `numComponents` and `numNonIsolated` change when the edge `p` is prepended to `es`:
components drop by at most one, do not drop when an endpoint is isolated in `es`, and grow by
one when both are; non-isolated vertices grow by at most the number of endpoints isolated in `es`.
-/
namespace Spqr
open Classical

theorem edgesConn_cons_cases {es : List (Nat × Nat)} {p : Nat × Nat} {a b : Nat}
    (h : EdgesConn (p :: es) a b) :
    EdgesConn es a b ∨ (EdgesConn es a p.1 ∧ EdgesConn es p.2 b) ∨
      (EdgesConn es a p.2 ∧ EdgesConn es p.1 b) := by
  induction h with
  | refl => exact Or.inl Relation.ReflTransGen.refl
  | @tail c b _ hcb ih =>
    rcases hcb with hm | hm <;> rcases List.mem_cons.1 hm with heq | hm
    · subst heq
      rcases ih with h | ⟨h1, _⟩ | ⟨h1, _⟩
      · exact Or.inr (Or.inl ⟨h, Relation.ReflTransGen.refl⟩)
      · exact Or.inr (Or.inl ⟨h1, Relation.ReflTransGen.refl⟩)
      · exact Or.inl h1
    · rcases ih with h | ⟨h1, h2⟩ | ⟨h1, h2⟩
      · exact Or.inl (h.tail (Or.inl hm))
      · exact Or.inr (Or.inl ⟨h1, h2.tail (Or.inl hm)⟩)
      · exact Or.inr (Or.inr ⟨h1, h2.tail (Or.inl hm)⟩)
    · subst heq
      rcases ih with h | ⟨h1, _⟩ | ⟨h1, _⟩
      · exact Or.inr (Or.inr ⟨h, Relation.ReflTransGen.refl⟩)
      · exact Or.inl h1
      · exact Or.inr (Or.inr ⟨h1, Relation.ReflTransGen.refl⟩)
    · rcases ih with h | ⟨h1, h2⟩ | ⟨h1, h2⟩
      · exact Or.inl (h.tail (Or.inr hm))
      · exact Or.inr (Or.inl ⟨h1, h2.tail (Or.inr hm)⟩)
      · exact Or.inr (Or.inr ⟨h1, h2.tail (Or.inr hm)⟩)

theorem hasEdge_of_edgesConn_ne {es : List (Nat × Nat)} {a b : Nat} (h : EdgesConn es a b)
    (ha : HasEdge es a) : HasEdge es b := hasEdge_of_edgesConn h ha

theorem isMin_unique {es : List (Nat × Nat)} {n a b : Nat} (ha : IsMin es n a) (hb : IsMin es n b)
    (han : a < n) (hbn : b < n) (h : EdgesConn es a b) : a = b :=
  le_antisymm (ha b hbn h) (hb a han (edgesConn_symm h))

/-- The least vertex of `w`'s component (`0` for `w` out of range). -/
noncomputable def compMin (es : List (Nat × Nat)) (n w : Nat) : Nat :=
  if h : ((Finset.range n).filter (EdgesConn es w)).Nonempty then Finset.min' _ h else 0

theorem compMin_spec {es : List (Nat × Nat)} {n w : Nat} (hw : w < n) :
    compMin es n w < n ∧ EdgesConn es w (compMin es n w) ∧ IsMin es n (compMin es n w) := by
  have hne : ((Finset.range n).filter (EdgesConn es w)).Nonempty :=
    ⟨w, Finset.mem_filter.2 ⟨Finset.mem_range.2 hw, Relation.ReflTransGen.refl⟩⟩
  unfold compMin
  rw [dif_pos hne]
  obtain ⟨h1, h2⟩ := Finset.mem_filter.1 (Finset.min'_mem _ hne)
  refine ⟨Finset.mem_range.1 h1, h2, fun x hx hmx => ?_⟩
  exact Finset.min'_le _ x (Finset.mem_filter.2 ⟨Finset.mem_range.2 hx, h2.trans hmx⟩)

/-- Least vertices of the components with an edge. -/
noncomputable def minSet (es : List (Nat × Nat)) (n : Nat) : Finset Nat :=
  (Finset.range n).filter (fun v => HasEdge es v ∧ IsMin es n v)

theorem ccCount_eq_card_minSet (es : List (Nat × Nat)) (n : Nat) :
    ccCount es n = (minSet es n).card := ccCount_eq_card_min es n

theorem mem_minSet {es : List (Nat × Nat)} {n v : Nat} :
    v ∈ minSet es n ↔ v < n ∧ HasEdge es v ∧ IsMin es n v := by
  simp [minSet]

section Cons
variable {es : List (Nat × Nat)} {p : Nat × Nat} {n : Nat}

theorem compMin_cons_mem {w : Nat} (hw : w ∈ minSet es n) :
    compMin (p :: es) n w ∈ minSet (p :: es) n := by
  obtain ⟨hwn, hwe, _⟩ := mem_minSet.1 hw
  obtain ⟨h1, h2, h3⟩ := compMin_spec (es := p :: es) hwn
  exact mem_minSet.2 ⟨h1, hasEdge_of_edgesConn h2
    (let ⟨q, hq, hz⟩ := hwe; ⟨q, List.mem_cons_of_mem _ hq, hz⟩), h3⟩

theorem compMin_cons_eq {w₁ w₂ : Nat} (h1 : w₁ ∈ minSet es n) (h2 : w₂ ∈ minSet es n)
    (h : compMin (p :: es) n w₁ = compMin (p :: es) n w₂) :
    w₁ = w₂ ∨ (EdgesConn es w₁ p.1 ∧ EdgesConn es p.2 w₂) ∨
      (EdgesConn es w₁ p.2 ∧ EdgesConn es p.1 w₂) := by
  obtain ⟨hn1, _, hm1⟩ := mem_minSet.1 h1
  obtain ⟨hn2, _, hm2⟩ := mem_minSet.1 h2
  have c : EdgesConn (p :: es) w₁ w₂ :=
    (compMin_spec (es := p :: es) hn1).2.1.trans (h ▸ edgesConn_symm (compMin_spec hn2).2.1)
  rcases edgesConn_cons_cases c with h | h | h
  · exact Or.inl (isMin_unique hm1 hm2 hn1 hn2 h)
  · exact Or.inr (Or.inl h)
  · exact Or.inr (Or.inr h)

theorem not_edgesConn_of_isolated {w x : Nat} (hw : w ∈ minSet es n) (hx : ¬HasEdge es x)
    (h : EdgesConn es w x) : False := by
  obtain ⟨_, hwe, _⟩ := mem_minSet.1 hw
  exact hx (hasEdge_of_edgesConn h hwe)

theorem card_minSet_le_succ (hp : p.1 < n ∧ p.2 < n) :
    (minSet es n).card ≤ (minSet (p :: es) n).card + 1 := by
  have hinj : Set.InjOn (compMin (p :: es) n) ((minSet es n).erase (compMin es n p.2)) := by
    intro w₁ hw₁ w₂ hw₂ h
    rw [Finset.coe_erase, Set.mem_diff, Finset.mem_coe, Set.mem_singleton_iff] at hw₁ hw₂
    rcases compMin_cons_eq hw₁.1 hw₂.1 h with h | ⟨_, h2⟩ | ⟨h1, _⟩
    · exact h
    · exfalso
      obtain ⟨hn2, _, hm2⟩ := mem_minSet.1 hw₂.1
      obtain ⟨hvn, hv1, hv2⟩ := compMin_spec (es := es) hp.2
      exact hw₂.2 (isMin_unique hm2 hv2 hn2 hvn ((edgesConn_symm h2).trans hv1))
    · exfalso
      obtain ⟨hn1, _, hm1⟩ := mem_minSet.1 hw₁.1
      obtain ⟨hvn, hv1, hv2⟩ := compMin_spec (es := es) hp.2
      exact hw₁.2 (isMin_unique hm1 hv2 hn1 hvn (h1.trans hv1))
  have := Finset.card_le_card_of_injOn _ (fun w hw => compMin_cons_mem (Finset.mem_of_mem_erase hw))
    hinj
  have := Finset.pred_card_le_card_erase (s := minSet es n) (a := compMin es n p.2)
  omega

theorem card_minSet_le_of_isolated (hu : ¬HasEdge es p.1 ∨ ¬HasEdge es p.2) :
    (minSet es n).card ≤ (minSet (p :: es) n).card := by
  refine Finset.card_le_card_of_injOn _ (fun w hw => compMin_cons_mem hw) ?_
  intro w₁ hw₁ w₂ hw₂ h
  rw [Finset.mem_coe] at hw₁ hw₂
  rcases compMin_cons_eq hw₁ hw₂ h with h | ⟨h1, h2⟩ | ⟨h1, h2⟩
  · exact h
  · exfalso
    rcases hu with hu | hv
    · exact not_edgesConn_of_isolated hw₁ hu h1
    · exact not_edgesConn_of_isolated hw₂ hv (edgesConn_symm h2)
  · exfalso
    rcases hu with hu | hv
    · exact not_edgesConn_of_isolated hw₂ hu (edgesConn_symm h2)
    · exact not_edgesConn_of_isolated hw₁ hv h1

theorem card_minSet_le_of_loop (hl : p.1 = p.2) :
    (minSet es n).card ≤ (minSet (p :: es) n).card := by
  refine Finset.card_le_card_of_injOn _ (fun w hw => compMin_cons_mem hw) ?_
  intro w₁ hw₁ w₂ hw₂ h
  rw [Finset.mem_coe] at hw₁ hw₂
  obtain ⟨hn1, _, hm1⟩ := mem_minSet.1 hw₁
  obtain ⟨hn2, _, hm2⟩ := mem_minSet.1 hw₂
  rcases compMin_cons_eq hw₁ hw₂ h with h | ⟨h1, h2⟩ | ⟨h1, h2⟩
  · exact h
  · exact isMin_unique hm1 hm2 hn1 hn2 (h1.trans (by rw [hl]; exact h2))
  · exact isMin_unique hm1 hm2 hn1 hn2 (h1.trans (by rw [← hl]; exact h2))

theorem card_minSet_succ_le (hp : p.1 < n ∧ p.2 < n) (hu : ¬HasEdge es p.1)
    (hv : ¬HasEdge es p.2) :
    (minSet es n).card + 1 ≤ (minSet (p :: es) n).card := by
  set m₀ := compMin (p :: es) n p.1
  obtain ⟨h0n, h0c, h0m⟩ := compMin_spec (es := p :: es) hp.1
  have hm₀ : m₀ ∈ minSet (p :: es) n :=
    mem_minSet.2 ⟨h0n, hasEdge_of_edgesConn h0c ⟨p, List.mem_cons_self .., Or.inl rfl⟩, h0m⟩
  have hsub : ∀ w ∈ minSet es n, compMin (p :: es) n w ∈ (minSet (p :: es) n).erase m₀ := by
    intro w hw
    refine Finset.mem_erase.2 ⟨fun h => ?_, compMin_cons_mem hw⟩
    obtain ⟨hwn, _, _⟩ := mem_minSet.1 hw
    have c : EdgesConn (p :: es) w p.1 :=
      (compMin_spec (es := p :: es) hwn).2.1.trans (h ▸ edgesConn_symm h0c)
    rcases edgesConn_cons_cases c with h | ⟨h1, _⟩ | ⟨h1, _⟩
    · exact not_edgesConn_of_isolated hw hu h
    · exact not_edgesConn_of_isolated hw hu h1
    · exact not_edgesConn_of_isolated hw hv h1
  have hinj : Set.InjOn (compMin (p :: es) n) (minSet es n) := by
    intro w₁ hw₁ w₂ hw₂ h
    rw [Finset.mem_coe] at hw₁ hw₂
    rcases compMin_cons_eq hw₁ hw₂ h with h | ⟨h1, _⟩ | ⟨_, h2⟩
    · exact h
    · exact (not_edgesConn_of_isolated hw₁ hu h1).elim
    · exact (not_edgesConn_of_isolated hw₂ hu (edgesConn_symm h2)).elim
  have := Finset.card_le_card_of_injOn _ hsub hinj
  rw [Finset.card_erase_of_mem hm₀] at this
  have := Finset.card_pos.2 ⟨m₀, hm₀⟩
  omega

theorem hes_cons (hes : ∀ q ∈ es, q.1 < n ∧ q.2 < n) (hp : p.1 < n ∧ p.2 < n) :
    ∀ q ∈ p :: es, q.1 < n ∧ q.2 < n := fun q hq => by
  rcases List.mem_cons.1 hq with rfl | hq
  · exact hp
  · exact hes q hq

/-- Deleting an edge creates at most one component. -/
theorem numComponents_cons_le_succ (hes : ∀ q ∈ es, q.1 < n ∧ q.2 < n)
    (hp : p.1 < n ∧ p.2 < n) : numComponents es n ≤ numComponents (p :: es) n + 1 := by
  rw [numComponents_eq_ccCount hes, numComponents_eq_ccCount (hes_cons hes hp),
    ccCount_eq_card_minSet, ccCount_eq_card_minSet]
  exact card_minSet_le_succ hp

/-- Deleting a pendant edge (an endpoint isolated afterwards) creates no component. -/
theorem numComponents_cons_le_of_isolated (hes : ∀ q ∈ es, q.1 < n ∧ q.2 < n)
    (hp : p.1 < n ∧ p.2 < n) (hu : ¬HasEdge es p.1 ∨ ¬HasEdge es p.2) :
    numComponents es n ≤ numComponents (p :: es) n := by
  rw [numComponents_eq_ccCount hes, numComponents_eq_ccCount (hes_cons hes hp),
    ccCount_eq_card_minSet, ccCount_eq_card_minSet]
  exact card_minSet_le_of_isolated hu

/-- Deleting a loop creates no component. -/
theorem numComponents_cons_le_of_loop (hes : ∀ q ∈ es, q.1 < n ∧ q.2 < n)
    (hp : p.1 < n ∧ p.2 < n) (hl : p.1 = p.2) :
    numComponents es n ≤ numComponents (p :: es) n := by
  rw [numComponents_eq_ccCount hes, numComponents_eq_ccCount (hes_cons hes hp),
    ccCount_eq_card_minSet, ccCount_eq_card_minSet]
  exact card_minSet_le_of_loop hl

/-- Deleting an isolated edge removes a component. -/
theorem numComponents_cons_succ_le (hes : ∀ q ∈ es, q.1 < n ∧ q.2 < n)
    (hp : p.1 < n ∧ p.2 < n) (hu : ¬HasEdge es p.1) (hv : ¬HasEdge es p.2) :
    numComponents es n + 1 ≤ numComponents (p :: es) n := by
  rw [numComponents_eq_ccCount hes, numComponents_eq_ccCount (hes_cons hes hp),
    ccCount_eq_card_minSet, ccCount_eq_card_minSet]
  exact card_minSet_succ_le hp hu hv

theorem numNonIsolated_cons_eq (hp : p.1 < n ∧ p.2 < n) :
    numNonIsolated (p :: es) n = numNonIsolated es n + (if HasEdge es p.1 then 0 else 1) +
      (if p.2 = p.1 ∨ HasEdge es p.2 then 0 else 1) := by
  rw [numNonIsolated_eq_card, numNonIsolated_eq_card]
  let T₁ : Finset Nat := if HasEdge es p.1 then ∅ else {p.1}
  let T₂ : Finset Nat := if p.2 = p.1 ∨ HasEdge es p.2 then ∅ else {p.2}
  have hT₁ : T₁.card = if HasEdge es p.1 then 0 else 1 := by
    simp only [T₁]; split_ifs <;> simp
  have hT₂ : T₂.card = if p.2 = p.1 ∨ HasEdge es p.2 then 0 else 1 := by
    simp only [T₂]; split_ifs <;> simp
  have heq : (Finset.range n).filter (HasEdge (p :: es)) =
      (Finset.range n).filter (HasEdge es) ∪ T₁ ∪ T₂ := by
    ext z
    constructor
    · intro hz
      obtain ⟨hzn, hz'⟩ := Finset.mem_filter.1 hz
      clear hz
      obtain ⟨q, hq, hz⟩ := hz'
      simp only [Finset.mem_union, Finset.mem_filter, T₁, T₂]
      rcases List.mem_cons.1 hq with rfl | hq
      · by_cases h1 : HasEdge es q.1
        · by_cases h2 : HasEdge es q.2
          · rcases hz with rfl | rfl
            · exact Or.inl (Or.inl ⟨hzn, h1⟩)
            · exact Or.inl (Or.inl ⟨hzn, h2⟩)
          · rcases hz with rfl | rfl
            · exact Or.inl (Or.inl ⟨hzn, h1⟩)
            · by_cases h21 : q.2 = q.1
              · exact Or.inl (Or.inl ⟨hzn, h21 ▸ h1⟩)
              · exact Or.inr (by
                  rw [ite_eq_right (not_or.2 ⟨h21, h2⟩)]; exact Finset.mem_singleton_self _)
        · rcases hz with rfl | rfl
          · exact Or.inl (Or.inr (by rw [ite_eq_right h1]; exact Finset.mem_singleton_self _))
          · by_cases h2 : HasEdge es q.2
            · exact Or.inl (Or.inl ⟨hzn, h2⟩)
            · by_cases h21 : q.2 = q.1
              · exact Or.inl (Or.inr (by
                  rw [h21, ite_eq_right h1]; exact Finset.mem_singleton_self _))
              · exact Or.inr (by
                  rw [ite_eq_right (not_or.2 ⟨h21, h2⟩)]; exact Finset.mem_singleton_self _)
      · exact Or.inl (Or.inl ⟨hzn, q, hq, hz⟩)
    · intro hz
      simp only [Finset.mem_union, Finset.mem_filter, T₁, T₂] at hz
      refine Finset.mem_filter.2 ?_
      rcases hz with (⟨hzn, q, hq, hz⟩ | hz) | hz
      · exact ⟨hzn, q, List.mem_cons_of_mem _ hq, hz⟩
      · split_ifs at hz with h1
        · simp at hz
        · rw [Finset.mem_singleton] at hz
          subst hz
          exact ⟨Finset.mem_range.2 hp.1, p, List.mem_cons_self .., Or.inl rfl⟩
      · split_ifs at hz with h2
        · simp at hz
        · rw [Finset.mem_singleton] at hz
          subst hz
          exact ⟨Finset.mem_range.2 hp.2, p, List.mem_cons_self .., Or.inr rfl⟩
  have hd1 : Disjoint ((Finset.range n).filter (HasEdge es)) T₁ := by
    simp only [T₁]
    split_ifs with h1
    · exact Finset.disjoint_empty_right _
    · rw [Finset.disjoint_singleton_right, Finset.mem_filter]
      exact fun h => h1 h.2
  have hd2 : Disjoint ((Finset.range n).filter (HasEdge es) ∪ T₁) T₂ := by
    simp only [T₂]
    split_ifs with h2
    · exact Finset.disjoint_empty_right _
    · rw [Finset.disjoint_singleton_right, Finset.mem_union, Finset.mem_filter]
      push_neg at h2
      simp only [T₁]
      split_ifs with h1 <;> simp [h2.1, h2.2]
  rw [heq, Finset.card_union_of_disjoint hd2, Finset.card_union_of_disjoint hd1, hT₁, hT₂]

end Cons
end Spqr
