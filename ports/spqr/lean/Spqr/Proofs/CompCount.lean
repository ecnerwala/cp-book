import Spqr.PlanarLayout
import Spqr.Proofs.TwoSumCount
import Mathlib.Data.Finset.Max

/-!
# Component counting (PROOF.md §8.5)

`numComponents es n` (`Planar.lean`: `n` rounds of min-label relaxation, then count the non-isolated
vertices that are their own label) equals `ccCount es n`, the number of connected components of
`es` (as vertex sets, under `EdgesConn`) containing an edge, whenever all edge endpoints are `< n`.

Proof: labels always stay in the component of the vertex and below it (`LabInv`); a vertex is
*converged* when its label is a lower bound of its component (`Conv`), i.e. equals the component
minimum. The minimum is converged from the start, and each round converges at least one more vertex
of every not-yet-converged component (an edge from a converged to an unconverged vertex), so after
`n` rounds every vertex is converged, and the counted vertices are exactly the component minima.
-/

namespace Spqr

open Classical

theorem edgesConn_symm {es : List (Nat × Nat)} {a b : Nat} (h : EdgesConn es a b) :
    EdgesConn es b a := by
  induction h with
  | refl => exact Relation.ReflTransGen.refl
  | tail _ hbc ih => exact Relation.ReflTransGen.trans (Relation.ReflTransGen.single hbc.symm) ih

/-- Connected components (as vertex sets) that contain an edge. -/
noncomputable def ccCount (es : List (Nat × Nat)) (n : Nat) : Nat :=
  (((Finset.range n).filter (HasEdge es)).image
    fun v => (Finset.range n).filter (EdgesConn es v)).card

theorem hasEdge_of_edgesConn {es : List (Nat × Nat)} {a b : Nat} (h : EdgesConn es a b)
    (ha : HasEdge es a) : HasEdge es b := by
  by_cases hab : a = b
  · exact hab ▸ ha
  · exact (HasEdge.of_edgesConn h hab).2

theorem cls_eq_of_edgesConn {es : List (Nat × Nat)} {n a b : Nat} (h : EdgesConn es a b) :
    (Finset.range n).filter (EdgesConn es a) = (Finset.range n).filter (EdgesConn es b) := by
  ext w
  simp only [Finset.mem_filter]
  exact and_congr_right fun _ =>
    ⟨fun h' => Relation.ReflTransGen.trans (edgesConn_symm h) h',
     fun h' => Relation.ReflTransGen.trans h h'⟩

/-- `v` is the least vertex of its component. -/
def IsMin (es : List (Nat × Nat)) (n v : Nat) : Prop := ∀ w, w < n → EdgesConn es v w → v ≤ w

/-- Components with an edge are counted by their least vertices. -/
theorem ccCount_eq_card_min (es : List (Nat × Nat)) (n : Nat) :
    ccCount es n = ((Finset.range n).filter (fun v => HasEdge es v ∧ IsMin es n v)).card := by
  symm
  rw [Finset.filter_congr (s := Finset.range n) (p := fun v => HasEdge es v ∧ IsMin es n v)
    (q := fun v => HasEdge es v ∧ ∀ w, w < n → EdgesConn es v w → v ≤ w) (fun _ _ => Iff.rfl)]
  unfold ccCount
  have himg : ((Finset.range n).filter (HasEdge es)).image
        (fun v => (Finset.range n).filter (EdgesConn es v)) =
      ((Finset.range n).filter
        (fun v => HasEdge es v ∧ ∀ w, w < n → EdgesConn es v w → v ≤ w)).image
        (fun v => (Finset.range n).filter (EdgesConn es v)) := by
    ext C
    simp only [Finset.mem_image, Finset.mem_filter, Finset.mem_range]
    constructor
    · rintro ⟨v, ⟨hv, he⟩, rfl⟩
      have hne : ((Finset.range n).filter (EdgesConn es v)).Nonempty :=
        ⟨v, Finset.mem_filter.2 ⟨Finset.mem_range.2 hv, Relation.ReflTransGen.refl⟩⟩
      obtain ⟨hmn, hvm⟩ := Finset.mem_filter.1 (Finset.min'_mem _ hne)
      rw [Finset.mem_range] at hmn
      refine ⟨_, ⟨hmn, hasEdge_of_edgesConn hvm he, fun w hw hmw => ?_⟩,
        (cls_eq_of_edgesConn hvm).symm⟩
      exact Finset.min'_le _ w
        (Finset.mem_filter.2 ⟨Finset.mem_range.2 hw, Relation.ReflTransGen.trans hvm hmw⟩)
    · rintro ⟨v, ⟨hv, he, -⟩, rfl⟩
      exact ⟨v, ⟨hv, he⟩, rfl⟩
  rw [himg]
  symm
  apply Finset.card_image_of_injOn
  intro a ha b hb hab
  have hab : (Finset.range n).filter (EdgesConn es a) =
      (Finset.range n).filter (EdgesConn es b) := hab
  rw [Finset.mem_coe, Finset.mem_filter, Finset.mem_range] at ha hb
  have hab' : EdgesConn es a b := by
    have : b ∈ (Finset.range n).filter (EdgesConn es a) := by
      rw [hab]; exact Finset.mem_filter.2 ⟨Finset.mem_range.2 hb.1, Relation.ReflTransGen.refl⟩
    exact (Finset.mem_filter.1 this).2
  exact le_antisymm (ha.2.2 b hb.1 hab') (hb.2.2 a ha.1 (edgesConn_symm hab'))

/-! ### Label lists -/

/-- Label of `v` in a label list (`v` itself when out of range). -/
def lab (l : List Nat) (v : Nat) : Nat := l[v]?.getD v

theorem lab_range (n v : Nat) : lab (List.range n) v = v := by
  unfold lab
  by_cases h : v < n
  · rw [List.getElem?_range h]; rfl
  · rw [List.getElem?_eq_none (by simpa using h)]; rfl

theorem lab_relaxEdge {l : List Nat} {p : Nat × Nat} (hp : p.1 < l.length ∧ p.2 < l.length)
    (v : Nat) :
    lab (relaxEdge l p) v =
      if v = p.1 ∨ v = p.2 then min (lab l p.1) (lab l p.2) else lab l v := by
  simp only [relaxEdge, lab]
  by_cases hvb : v = p.2
  · rw [hvb, List.getElem?_set_self (by rw [List.length_set]; exact hp.2)]; simp
  · rw [List.getElem?_set_ne (Ne.symm hvb)]
    by_cases hva : v = p.1
    · rw [hva, List.getElem?_set_self hp.1]; simp
    · rw [List.getElem?_set_ne (Ne.symm hva)]; simp [hva, hvb]

theorem lab_relaxEdge_le {l : List Nat} {p : Nat × Nat} (hp : p.1 < l.length ∧ p.2 < l.length)
    (v : Nat) : lab (relaxEdge l p) v ≤ lab l v := by
  rw [lab_relaxEdge hp]
  split_ifs with h
  · rcases h with rfl | rfl
    · exact min_le_left _ _
    · exact min_le_right _ _
  · exact le_rfl

theorem foldl_relaxEdge_length :
    ∀ (es : List (Nat × Nat)) (l : List Nat), (es.foldl relaxEdge l).length = l.length := by
  intro es
  induction es with
  | nil => intro l; rfl
  | cons p es ih => intro l; rw [List.foldl_cons, ih, relaxEdge_length]

/-- Labels only decrease. -/
theorem lab_foldl_le (n : Nat) : ∀ (es : List (Nat × Nat)), (∀ p ∈ es, p.1 < n ∧ p.2 < n) →
    ∀ (l : List Nat), l.length = n → ∀ v, lab (es.foldl relaxEdge l) v ≤ lab l v := by
  intro es
  induction es with
  | nil => intro _ l _ v; exact le_rfl
  | cons p es ih =>
    intro hes l hl v
    rw [List.foldl_cons]
    have hp := hes p (List.mem_cons_self ..)
    refine le_trans (ih (fun q hq => hes q (List.mem_cons.2 (Or.inr hq))) _
      (by rw [relaxEdge_length, hl]) v) ?_
    exact lab_relaxEdge_le ⟨by omega, by omega⟩ v

/-- After a round, each endpoint's label is at most the other endpoint's label before it. -/
theorem lab_foldl_edge (n : Nat) : ∀ (es : List (Nat × Nat)), (∀ p ∈ es, p.1 < n ∧ p.2 < n) →
    ∀ (l : List Nat), l.length = n → ∀ p ∈ es,
      lab (es.foldl relaxEdge l) p.2 ≤ lab l p.1 ∧ lab (es.foldl relaxEdge l) p.1 ≤ lab l p.2 := by
  intro es
  induction es with
  | nil => intro _ _ _ p hp; simp at hp
  | cons q es ih =>
    intro hes l hl p hp
    rw [List.foldl_cons]
    have hq := hes q (List.mem_cons_self ..)
    have hes' : ∀ p ∈ es, p.1 < n ∧ p.2 < n := fun r hr => hes r (List.mem_cons.2 (Or.inr hr))
    have hl' : (relaxEdge l q).length = n := by rw [relaxEdge_length, hl]
    have hmono := lab_foldl_le n es hes' _ hl'
    have hql : q.1 < l.length ∧ q.2 < l.length := ⟨by omega, by omega⟩
    rcases List.mem_cons.1 hp with hpq | hp
    · rw [hpq]
      constructor
      · refine le_trans (hmono q.2) ?_
        rw [lab_relaxEdge hql, if_pos (Or.inr rfl)]; exact min_le_left _ _
      · refine le_trans (hmono q.1) ?_
        rw [lab_relaxEdge hql, if_pos (Or.inl rfl)]; exact min_le_right _ _
    · obtain ⟨h1, h2⟩ := ih hes' _ hl' p hp
      exact ⟨le_trans h1 (lab_relaxEdge_le hql _), le_trans h2 (lab_relaxEdge_le hql _)⟩

/-- Labels lie in the component of the vertex and below it. -/
def LabInv (es : List (Nat × Nat)) (l : List Nat) : Prop :=
  ∀ v, EdgesConn es v (lab l v) ∧ lab l v ≤ v

theorem labInv_range (es : List (Nat × Nat)) (n : Nat) : LabInv es (List.range n) := fun v => by
  rw [lab_range]; exact ⟨Relation.ReflTransGen.refl, le_rfl⟩

theorem labInv_relaxEdge {es : List (Nat × Nat)} {l : List Nat} {p : Nat × Nat} (hp : p ∈ es)
    (hl : p.1 < l.length ∧ p.2 < l.length) (h : LabInv es l) : LabInv es (relaxEdge l p) := by
  intro v
  rw [lab_relaxEdge hl]
  split_ifs with hv
  · have hab : EdgesConn es p.1 p.2 := Relation.ReflTransGen.single (Or.inl hp)
    have h1 := h p.1; have h2 := h p.2
    rcases le_total (lab l p.1) (lab l p.2) with hle | hle
    · rw [min_eq_left hle]
      rcases hv with rfl | rfl
      · exact h1
      · exact ⟨Relation.ReflTransGen.trans (edgesConn_symm hab) h1.1, by omega⟩
    · rw [min_eq_right hle]
      rcases hv with rfl | rfl
      · exact ⟨Relation.ReflTransGen.trans hab h2.1, by omega⟩
      · exact h2
  · exact h v

theorem labInv_foldl (n : Nat) (es₀ : List (Nat × Nat)) :
    ∀ (es : List (Nat × Nat)), (∀ p ∈ es, p ∈ es₀) → (∀ p ∈ es, p.1 < n ∧ p.2 < n) →
      ∀ (l : List Nat), l.length = n → LabInv es₀ l → LabInv es₀ (es.foldl relaxEdge l) := by
  intro es
  induction es with
  | nil => intro _ _ l _ h; exact h
  | cons p es ih =>
    intro hsub hes l hl h
    rw [List.foldl_cons]
    have hp := hes p (List.mem_cons_self ..)
    exact ih (fun q hq => hsub q (List.mem_cons.2 (Or.inr hq)))
      (fun q hq => hes q (List.mem_cons.2 (Or.inr hq))) _ (by rw [relaxEdge_length, hl])
      (labInv_relaxEdge (hsub p (List.mem_cons_self ..)) ⟨by omega, by omega⟩ h)

/-- `v` is converged: its label is a lower bound of its component. -/
def Conv (es : List (Nat × Nat)) (n : Nat) (l : List Nat) (v : Nat) : Prop :=
  ∀ w, w < n → EdgesConn es v w → lab l v ≤ w

theorem conv_mono {es : List (Nat × Nat)} {n : Nat} {l l' : List Nat}
    (h : ∀ v, lab l' v ≤ lab l v) {v : Nat} (hc : Conv es n l v) : Conv es n l' v :=
  fun w hw hvw => le_trans (h v) (hc w hw hvw)

/-- A path from a `P`-vertex to a non-`P`-vertex crosses an edge. -/
theorem exists_crossing {es : List (Nat × Nat)} (P : Nat → Prop) {m v : Nat}
    (h : EdgesConn es m v) (hm : P m) (hv : ¬P v) :
    ∃ a b, P a ∧ ¬P b ∧ ((a, b) ∈ es ∨ (b, a) ∈ es) := by
  induction h with
  | refl => exact absurd hm hv
  | @tail b c _ hbc ih =>
    by_cases hb : P b
    · exact ⟨b, c, hb, hv, hbc⟩
    · exact ih hb

/-- A round converges a new vertex of every component that has an unconverged vertex. -/
theorem conv_round (n : Nat) (es : List (Nat × Nat)) (hes : ∀ p ∈ es, p.1 < n ∧ p.2 < n)
    (l : List Nat) (hl : l.length = n) (hinv : LabInv es l) {v : Nat} (hv : v < n)
    (hnc : ¬Conv es n l v) :
    ∃ b, b < n ∧ ¬Conv es n l b ∧ Conv es n (es.foldl relaxEdge l) b := by
  have hmem : ∀ w, w ∈ (Finset.range n).filter (EdgesConn es v) ↔ w < n ∧ EdgesConn es v w :=
    fun w => by rw [Finset.mem_filter, Finset.mem_range]
  have hne : ((Finset.range n).filter (EdgesConn es v)).Nonempty :=
    ⟨v, (hmem v).2 ⟨hv, Relation.ReflTransGen.refl⟩⟩
  obtain ⟨hmn, hvm⟩ := (hmem _).1 (Finset.min'_mem _ hne)
  have hconv_m : Conv es n l (Finset.min' _ hne) := by
    intro w hw hmw
    exact le_trans (hinv _).2
      (Finset.min'_le _ w ((hmem w).2 ⟨hw, Relation.ReflTransGen.trans hvm hmw⟩))
  obtain ⟨a, b, ha, hb, hab⟩ := exists_crossing (Conv es n l) (edgesConn_symm hvm) hconv_m hnc
  have hbn : b < n := by
    rcases hab with h | h
    · exact (hes _ h).2
    · exact (hes _ h).1
  refine ⟨b, hbn, hb, fun w hw hbw => ?_⟩
  have hconn_ab : EdgesConn es a b := Relation.ReflTransGen.single hab
  have h1 : lab (es.foldl relaxEdge l) b ≤ lab l a := by
    rcases hab with h | h
    · exact (lab_foldl_edge n es hes l hl (a, b) h).1
    · exact (lab_foldl_edge n es hes l hl (b, a) h).2
  exact le_trans h1 (ha w hw (Relation.ReflTransGen.trans hconn_ab hbw))

/-! ### Rounds -/

def rounds (es : List (Nat × Nat)) (n k : Nat) : List Nat := (relaxLabels es)^[k] (List.range n)

theorem rounds_succ (es : List (Nat × Nat)) (n k : Nat) :
    rounds es n (k + 1) = es.foldl relaxEdge (rounds es n k) := by
  unfold rounds; rw [Function.iterate_succ_apply', relaxLabels_eq]

theorem rounds_length (es : List (Nat × Nat)) (n : Nat) : ∀ k, (rounds es n k).length = n := by
  intro k
  induction k with
  | zero => simp [rounds]
  | succ k ih => rw [rounds_succ, foldl_relaxEdge_length, ih]

theorem compLabels_eq_rounds (es : List (Nat × Nat)) (n : Nat) :
    compLabels es n = rounds es n n := rfl

section
variable {es : List (Nat × Nat)} {n : Nat} (hes : ∀ p ∈ es, p.1 < n ∧ p.2 < n)
include hes

theorem rounds_labInv : ∀ k, LabInv es (rounds es n k) := by
  intro k
  induction k with
  | zero => exact labInv_range es n
  | succ k ih =>
    rw [rounds_succ]; exact labInv_foldl n es es (fun _ h => h) hes _ (rounds_length es n k) ih

theorem rounds_mono (k v : Nat) : lab (rounds es n (k + 1)) v ≤ lab (rounds es n k) v := by
  rw [rounds_succ]; exact lab_foldl_le n es hes _ (rounds_length es n k) v

theorem rounds_conv : ∀ k, (∀ v, v < n → Conv es n (rounds es n k) v) ∨
    k ≤ ((Finset.range n).filter (Conv es n (rounds es n k))).card := by
  intro k
  induction k with
  | zero => exact Or.inr (Nat.zero_le _)
  | succ k ih =>
    by_cases hall : ∀ v, v < n → Conv es n (rounds es n k) v
    · exact Or.inl fun v hv => conv_mono (rounds_mono hes k) (hall v hv)
    · right
      have hk := ih.resolve_left hall
      push_neg at hall
      obtain ⟨v, hv, hnc⟩ := hall
      obtain ⟨b, hb, hnb, hcb⟩ :=
        conv_round n es hes _ (rounds_length es n k) (rounds_labInv hes k) hv hnc
      rw [← rounds_succ] at hcb
      have hsub : (Finset.range n).filter (Conv es n (rounds es n k)) ⊆
          (Finset.range n).filter (Conv es n (rounds es n (k + 1))) := by
        intro x hx
        rw [Finset.mem_filter] at hx ⊢
        exact ⟨hx.1, conv_mono (rounds_mono hes k) hx.2⟩
      have hlt := Finset.card_lt_card ((Finset.ssubset_iff_of_subset hsub).2
        ⟨b, Finset.mem_filter.2 ⟨Finset.mem_range.2 hb, hcb⟩,
          fun h => hnb (Finset.mem_filter.1 h).2⟩)
      omega

theorem compLabels_conv (v : Nat) (hv : v < n) : Conv es n (compLabels es n) v := by
  rw [compLabels_eq_rounds]
  rcases rounds_conv hes n with h | h
  · exact h v hv
  · have hcard : ((Finset.range n).filter (Conv es n (rounds es n n))).card = n :=
      le_antisymm
        (by simpa using Finset.card_filter_le (Finset.range n) (Conv es n (rounds es n n))) h
    have heq := Finset.eq_of_subset_of_card_le (Finset.filter_subset _ (Finset.range n))
      (by rw [hcard, Finset.card_range])
    have hv' : v ∈ (Finset.range n).filter (Conv es n (rounds es n n)) := by
      rw [heq]; exact Finset.mem_range.2 hv
    exact (Finset.mem_filter.1 hv').2

/-- After `n` rounds a vertex is its own label iff it is the minimum of its component. -/
theorem lab_compLabels_eq_iff (v : Nat) (hv : v < n) :
    lab (compLabels es n) v = v ↔ ∀ w, w < n → EdgesConn es v w → v ≤ w := by
  have hc := compLabels_conv hes v hv
  have hi := rounds_labInv hes n v
  rw [← compLabels_eq_rounds] at hi
  constructor
  · intro h w hw hvw; rw [← h]; exact hc w hw hvw
  · intro h
    exact le_antisymm hi.2 (h _ (by omega) hi.1)

/-- The executable component counter counts the components (vertex classes) with an edge. -/
theorem numComponents_eq_ccCount : numComponents es n = ccCount es n := by
  have hfilt : ∀ v, v < n →
      ((nonIsolated es v && (compLabels es n)[v]?.getD v == v) = true ↔
        HasEdge es v ∧ ∀ w, w < n → EdgesConn es v w → v ≤ w) := by
    intro v hv
    rw [Bool.and_eq_true, beq_iff_eq, nonIsolated_iff]
    exact and_congr_right fun _ => lab_compLabels_eq_iff hes v hv
  show ((List.range n).filter
    fun v => nonIsolated es v && (compLabels es n)[v]?.getD v == v).length = _
  rw [length_filter_range_eq_card, ccCount_eq_card_min]
  congr 1
  exact Finset.filter_congr fun v hv => hfilt v (Finset.mem_range.1 hv)

end

end Spqr
