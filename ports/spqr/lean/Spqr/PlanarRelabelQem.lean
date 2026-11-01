import Spqr.PlanarRelabelSem

/-!
# `qem` bookkeeping for `capLinked`/`flipped`

Slot-wise agreement of two `qem` arrays on a set of slots (`QemAgree`) is preserved by `set!` (the
written slot joins the set) and by `swapIfInBounds` of two slots in the set; `capLinked` and
`flipped` leave every slot they do not address unchanged.
-/

namespace Spqr

/-- `a` and `b` have the same size and agree on every slot satisfying `S`. -/
def QemAgree (S : Nat → Prop) (a b : Qem) : Prop := a.size = b.size ∧ ∀ q, S q → a[q]! = b[q]!

theorem QemAgree.mono {S T : Nat → Prop} {a b : Qem} (h : ∀ q, T q → S q) (ha : QemAgree S a b) :
    QemAgree T a b := ⟨ha.1, fun q hq => ha.2 q (h q hq)⟩

theorem QemAgree.set! {S : Nat → Prop} {a b : Qem} (h : QemAgree S a b) (i : Nat) (v : Option Nat) :
    QemAgree (fun q => S q ∨ q = i) (a.set! i v) (b.set! i v) := by
  refine ⟨by simp only [Array.set!_eq_setIfInBounds, Array.size_setIfInBounds, h.1], fun q hq => ?_⟩
  by_cases hqi : q = i
  · subst hqi
    by_cases hlt : q < a.size
    · rw [Ghost.set!_get!_self a v hlt, Ghost.set!_get!_self b v (h.1 ▸ hlt)]
    · rw [getElem!_neg _ _ (by simpa [Array.set!_eq_setIfInBounds] using hlt),
        getElem!_neg _ _ (by simpa [Array.set!_eq_setIfInBounds, ← h.1] using hlt)]
  · rw [Ghost.set!_get!_ne a v (Ne.symm hqi), Ghost.set!_get!_ne b v (Ne.symm hqi)]
    exact h.2 q (hq.resolve_right hqi)

theorem swapIfInBounds_get! (a : Qem) (i j q : Nat) :
    (a.swapIfInBounds i j)[q]! =
      if i < a.size ∧ j < a.size then (if q = i then a[j]! else if q = j then a[i]! else a[q]!) else a[q]! := by
  rw [Array.swapIfInBounds_def]
  by_cases hi : i < a.size
  · by_cases hj : j < a.size
    · rw [dif_pos hi, dif_pos hj, if_pos ⟨hi, hj⟩]
      by_cases hq : q < a.size
      · rw [getElem!_pos _ _ (by simpa using hq), Array.getElem_swap, getElem!_pos a q hq]
        split_ifs with h1 h2
        · subst h1; rw [getElem!_pos a j hj]
        · subst h2; rw [getElem!_pos a i hi]
        · rfl
      · rw [getElem!_neg _ _ (by simpa using hq), getElem!_neg a q hq]
        split_ifs with h1 h2
        · subst h1; exact absurd hi hq
        · subst h2; exact absurd hj hq
        · rfl
    · rw [dif_pos hi, dif_neg hj, if_neg (fun h => hj h.2)]
  · rw [dif_neg hi, if_neg (fun h => hi h.1)]

theorem QemAgree.swap {S : Nat → Prop} {a b : Qem} (h : QemAgree S a b) {i j : Nat} (hi : S i) (hj : S j) :
    QemAgree S (a.swapIfInBounds i j) (b.swapIfInBounds i j) := by
  refine ⟨by simp [h.1], fun q hq => ?_⟩
  rw [swapIfInBounds_get!, swapIfInBounds_get!, h.1]
  split_ifs <;> first | exact h.2 _ hi | exact h.2 _ hj | exact h.2 _ hq

theorem swapIfInBounds_get!_of_ne (a : Qem) (i j q : Nat) (hi : q ≠ i) (hj : q ≠ j) :
    (a.swapIfInBounds i j)[q]! = a[q]! := by
  rw [swapIfInBounds_get!]
  split_ifs <;> rfl

theorem size_capLinked (g : Graph) (a : Qem) (m : Array Nat) : (capLinked g a m).size = a.size := by
  rw [capLinked_eq]; simp only [Array.set!_eq_setIfInBounds, Array.size_setIfInBounds]

theorem size_flipped (g : Graph) (ch : List ItemId) (fl : List Bool) (a : Qem) :
    (flipped g ch fl a).size = a.size := by
  unfold flipped
  generalize ch.zip fl = l
  induction l generalizing a with
  | nil => rfl
  | cons x l ih =>
    obtain ⟨c, f⟩ := x
    rw [List.foldl_cons, ih]
    dsimp only
    split <;> simp

theorem capLinked_get!_of (g : Graph) (a : Qem) (m : Array Nat) (q : Nat) (hq : q < 8 * g.ne)
    (hm : ∀ s, s < 4 → q ≠ m[s]!) : (capLinked g a m)[q]! = a[q]! := by
  rw [capLinked_eq]
  simp only
  rw [Ghost.set!_get!_ne _ _ (Ne.symm (hm 3 (by omega))), Ghost.set!_get!_ne _ _ (show 8 * g.ne + 3 ≠ q by omega),
    Ghost.set!_get!_ne _ _ (Ne.symm (hm 2 (by omega))), Ghost.set!_get!_ne _ _ (show 8 * g.ne + 2 ≠ q by omega),
    Ghost.set!_get!_ne _ _ (Ne.symm (hm 1 (by omega))), Ghost.set!_get!_ne _ _ (show 8 * g.ne + 1 ≠ q by omega),
    Ghost.set!_get!_ne _ _ (Ne.symm (hm 0 (by omega))), Ghost.set!_get!_ne _ _ (show 8 * g.ne + 0 ≠ q by omega)]

theorem flipped_get!_of (g : Graph) (ch : List ItemId) (fl : List Bool) (a : Qem) (q : Nat)
    (hq : ∀ c ∈ ch, 1 + g.nv ≤ c → ¬ (4 * (c - (1 + g.nv)) ≤ q ∧ q < 4 * (c - (1 + g.nv)) + 4)) :
    (flipped g ch fl a)[q]! = a[q]! := by
  unfold flipped
  have key : ∀ (l : List (ItemId × Bool)), (∀ x ∈ l, x.1 ∈ ch) → ∀ a : Qem,
      (l.foldl (fun q (c, f) =>
        if c ≥ 1 + g.nv && f then
          (q.swapIfInBounds (4 * (c - (1 + g.nv)) + 0) (4 * (c - (1 + g.nv)) + 1)).swapIfInBounds
            (4 * (c - (1 + g.nv)) + 2) (4 * (c - (1 + g.nv)) + 3)
        else q) a)[q]! = a[q]! := by
    intro l
    induction l with
    | nil => intro _ a; rfl
    | cons x l ih =>
      obtain ⟨c, f⟩ := x
      intro hl a
      rw [List.foldl_cons, ih (fun y hy => hl y (List.mem_cons_of_mem _ hy))]
      dsimp only
      split
      · rename_i hc
        have hc' : 1 + g.nv ≤ c := by simpa using (Bool.and_eq_true_iff.1 hc).1
        have := hq c (hl (c, f) (List.mem_cons_self ..)) hc'
        have hq' : q < 4 * (c - (1 + g.nv)) ∨ 4 * (c - (1 + g.nv)) + 4 ≤ q := by
          rcases Nat.lt_or_ge q (4 * (c - (1 + g.nv))) with h1 | h1
          · exact Or.inl h1
          · rcases Nat.lt_or_ge q (4 * (c - (1 + g.nv)) + 4) with h2 | h2
            · exact absurd ⟨h1, h2⟩ this
            · exact Or.inr h2
        rw [swapIfInBounds_get!_of_ne _ (4 * (c - (1 + g.nv)) + 2) (4 * (c - (1 + g.nv)) + 3) q (by simp only [ItemId] at hq' ⊢; omega) (by simp only [ItemId] at hq' ⊢; omega),
          swapIfInBounds_get!_of_ne _ (4 * (c - (1 + g.nv)) + 0) (4 * (c - (1 + g.nv)) + 1) q (by simp only [ItemId] at hq' ⊢; omega) (by simp only [ItemId] at hq' ⊢; omega)]
      · rfl
  exact key _ (fun x hx => by obtain ⟨a, b⟩ := x; exact (List.of_mem_zip hx).1) a

theorem capLinked_agree (g : Graph) {S : Nat → Prop} {a b : Qem} (h : QemAgree S a b) (m : Array Nat) :
    QemAgree (fun q => S q ∨ (8 * g.ne ≤ q ∧ q < 8 * g.ne + 4) ∨ ∃ s, s < 4 ∧ q = m[s]!)
      (capLinked g a m) (capLinked g b m) := by
  rw [capLinked_eq, capLinked_eq]
  simp only
  refine QemAgree.mono ?_ ((((((((h.set! _ _).set! _ _).set! _ _).set! _ _).set! _ _).set! _ _).set! _ _).set! _ _)
  intro q hq
  rcases hq with hq | ⟨h1, h2⟩ | ⟨s, hs, rfl⟩
  · simp only [hq, true_or]
  · have : q = 8 * g.ne + 0 ∨ q = 8 * g.ne + 1 ∨ q = 8 * g.ne + 2 ∨ q = 8 * g.ne + 3 := by omega
    rcases this with rfl | rfl | rfl | rfl <;> simp
  · match s, hs with
    | 0, _ => simp
    | 1, _ => simp
    | 2, _ => simp
    | 3, _ => simp

theorem flipped_agree (g : Graph) (ch : List ItemId) (fl : List Bool) {S : Nat → Prop} {a b : Qem}
    (hS : ∀ c ∈ ch, 1 + g.nv ≤ c → ∀ z, z < 4 → S (4 * (c - (1 + g.nv)) + z)) (h : QemAgree S a b) :
    QemAgree S (flipped g ch fl a) (flipped g ch fl b) := by
  unfold flipped
  have key : ∀ (l : List (ItemId × Bool)), (∀ x ∈ l, x.1 ∈ ch) → ∀ a b : Qem, QemAgree S a b →
      QemAgree S (l.foldl (fun q (c, f) =>
        if c ≥ 1 + g.nv && f then
          (q.swapIfInBounds (4 * (c - (1 + g.nv)) + 0) (4 * (c - (1 + g.nv)) + 1)).swapIfInBounds
            (4 * (c - (1 + g.nv)) + 2) (4 * (c - (1 + g.nv)) + 3)
        else q) a)
      (l.foldl (fun q (c, f) =>
        if c ≥ 1 + g.nv && f then
          (q.swapIfInBounds (4 * (c - (1 + g.nv)) + 0) (4 * (c - (1 + g.nv)) + 1)).swapIfInBounds
            (4 * (c - (1 + g.nv)) + 2) (4 * (c - (1 + g.nv)) + 3)
        else q) b) := by
    intro l
    induction l with
    | nil => intro _ a b h; exact h
    | cons x l ih =>
      obtain ⟨c, f⟩ := x
      intro hl a b h
      rw [List.foldl_cons, List.foldl_cons]
      refine ih (fun y hy => hl y (List.mem_cons_of_mem _ hy)) _ _ ?_
      dsimp only
      split
      · rename_i hc
        have hc' : 1 + g.nv ≤ c := by simpa using (Bool.and_eq_true_iff.1 hc).1
        have hx : c ∈ ch := hl (c, f) (List.mem_cons_self ..)
        exact (h.swap (hS _ hx hc' 0 (by omega)) (hS _ hx hc' 1 (by omega))).swap
          (hS _ hx hc' 2 (by omega)) (hS _ hx hc' 3 (by omega))
      · exact h
  exact key _ (fun x hx => by obtain ⟨a, b⟩ := x; exact (List.of_mem_zip hx).1) a b h

/-- `capLinked` writes the same values at the same slots, so agreement spreads to the cap slots. -/
theorem qemAgree_capLinked {S : Nat → Prop} {a b : Qem} (g : Graph) (m : Array Nat) (h : QemAgree S a b) :
    QemAgree (fun q => S q ∨ ∃ s, s < 4 ∧ q = 4 * capVe g + s) (capLinked g a m) (capLinked g b m) := by
  unfold capLinked
  have e : List.range 4 = [0, 1, 2, 3] := rfl
  rw [e]
  simp only [List.foldl_cons, List.foldl_nil]
  refine (((((((h.set! _ _).set! _ _).set! _ _).set! _ _).set! _ _).set! _ _).set! _ _).set! _ _ |>.mono ?_
  rintro q (hq | ⟨s, hs, rfl⟩)
  · exact Or.inl (Or.inl (Or.inl (Or.inl (Or.inl (Or.inl (Or.inl (Or.inl hq)))))))
  · rcases (show s = 0 ∨ s = 1 ∨ s = 2 ∨ s = 3 by omega) with rfl | rfl | rfl | rfl
    · exact Or.inl (Or.inl (Or.inl (Or.inl (Or.inl (Or.inl (Or.inl (Or.inr rfl)))))))
    · exact Or.inl (Or.inl (Or.inl (Or.inl (Or.inl (Or.inr rfl)))))
    · exact Or.inl (Or.inl (Or.inl (Or.inr rfl)))
    · exact Or.inl (Or.inr rfl)

/-- `flipped` only swaps slots of the edge children; if those are in `S`, agreement is kept. -/
theorem qemAgree_flipped {S : Nat → Prop} {a b : Qem} (g : Graph) (ch : List ItemId) (fl : List Bool)
    (hS : ∀ c ∈ ch, 1 + g.nv ≤ c → ∀ z, z < 4 → S (4 * (c - (1 + g.nv)) + z))
    (h : QemAgree S a b) : QemAgree S (flipped g ch fl a) (flipped g ch fl b) := by
  unfold flipped
  induction ch generalizing fl a b with
  | nil => simpa using h
  | cons c ch ih =>
    cases fl with
    | nil => simpa using h
    | cons f fl =>
      simp only [List.zip_cons_cons, List.foldl_cons]
      refine ih fl (fun c' hc' => hS c' (List.mem_cons_of_mem _ hc')) ?_
      try dsimp only
      split
      · rename_i hcf
        have hge : 1 + g.nv ≤ c := by
          have h' := hcf
          simp only [Bool.and_eq_true, decide_eq_true_eq, ge_iff_le] at h'
          exact h'.1
        have hc : c ∈ c :: ch := List.mem_cons_self ..
        exact (h.swap (hS c hc hge 0 (by omega)) (hS c hc hge 1 (by omega))).swap
          (hS c hc hge 2 (by omega)) (hS c hc hge 3 (by omega))
      · exact h

theorem range4_get! (z : Nat) (hz : z < 4) : (Array.range 4)[z]! = z := by
  rw [getElem!_pos _ _ (by simpa using hz)]; simp

end Spqr
