import Spqr.StBdPop
/-!
# Growth of the open block

The block an open stack segment belongs to grows by inserting fresh items (the pieces of a new
out-edge, the vertex item of the current vertex) around the items already present: a list
`A ++ M ++ B` becomes `A ++ P ++ M ++ Q ++ B`.  `InBlock` / `StLive` are monotone under such an
insertion as long as the inserted items are fresh (`InBlock.insert`, `StLive.insert`).
-/

namespace Spqr

/-! ### Segments and insertion -/

theorem prefix_of_prefix_append_middle {α : Type} {L P Y : List α} :
    ∀ (X : List α), L <+: X ++ P ++ Y → (∀ x ∈ L, x ∉ P) → L <+: X ++ Y
  | [], h, hd => by
    simp only [List.nil_append] at h ⊢
    cases P with
    | nil => simpa using h
    | cons p P =>
      rcases List.prefix_cons_iff.1 h with rfl | ⟨t, rfl, -⟩
      · exact List.nil_prefix
      · exact absurd List.mem_cons_self (hd p List.mem_cons_self)
  | x :: X, h, hd => by
    simp only [List.cons_append] at h ⊢
    rcases List.prefix_cons_iff.1 h with rfl | ⟨t, rfl, ht⟩
    · exact List.nil_prefix
    · exact List.prefix_cons_iff.2 (Or.inr ⟨t, rfl,
        prefix_of_prefix_append_middle X ht fun y hy => hd y (List.mem_cons_of_mem _ hy)⟩)

theorem infix_of_infix_append_middle {α : Type} {L P Y : List α} :
    ∀ (X : List α), L <:+: X ++ P ++ Y → (∀ x ∈ L, x ∉ P) → L <:+: X ++ Y
  | [], h, hd => by
    simp only [List.nil_append] at h ⊢
    induction P with
    | nil => simpa using h
    | cons p P ih =>
      rcases List.infix_cons_iff.1 h with h | h
      · rcases List.prefix_cons_iff.1 h with rfl | ⟨t, rfl, -⟩
        · exact List.nil_infix
        · exact absurd List.mem_cons_self (hd p List.mem_cons_self)
      · exact ih h fun y hy hP => hd y hy (List.mem_cons_of_mem _ hP)
  | x :: X, h, hd => by
    simp only [List.cons_append] at h ⊢
    rcases List.infix_cons_iff.1 h with h | h
    · exact (prefix_of_prefix_append_middle (x :: X) (by simpa using h) hd).isInfix
    · exact List.infix_cons_iff.2 (Or.inr (infix_of_infix_append_middle X h hd))

/-- A segment of the grown list consisting of old items is a segment of the old list. -/
theorem segment_of_insert {A P M Q B L : List ItemId}
    (hd : ∀ x ∈ L, x ∉ P ++ Q) (h : ∃ A' B', A ++ P ++ M ++ Q ++ B = A' ++ L ++ B') :
    ∃ A' B', A ++ M ++ B = A' ++ L ++ B' := by
  obtain ⟨A', B', h⟩ := h
  have h1 : L <:+: (A ++ P ++ M) ++ Q ++ B := ⟨A', B', h.symm⟩
  have h2 := infix_of_infix_append_middle _ h1 fun x hx => fun hP => hd x hx (List.mem_append_right _ hP)
  have h3 : L <:+: A ++ P ++ (M ++ B) := by rw [← List.append_assoc]; exact h2
  obtain ⟨s, t, hst⟩ := infix_of_infix_append_middle _ h3
    fun x hx => fun hP => hd x hx (List.mem_append_left _ hP)
  exact ⟨s, t, by rw [List.append_assoc]; exact hst.symm⟩

theorem segment_insert {A P M Q B L : List ItemId}
    (hL : L <:+: A ∨ L <:+: M ∨ L <:+: B) :
    ∃ A' B', A ++ P ++ M ++ Q ++ B = A' ++ L ++ B' := by
  rcases hL with ⟨s, t, h⟩ | ⟨s, t, h⟩ | ⟨s, t, h⟩
  · exact ⟨s, t ++ P ++ M ++ Q ++ B, by rw [← h]; simp only [List.append_assoc]⟩
  · exact ⟨A ++ P ++ s, t ++ Q ++ B, by rw [← h]; simp only [List.append_assoc]⟩
  · exact ⟨A ++ P ++ M ++ Q ++ s, t, by rw [← h]; simp only [List.append_assoc]⟩

/-! ### `Precedes` under insertion -/

theorem idxOf_lt_length_of_mem_left {a : Nat} {X Y : List Nat} (h : a ∈ X) :
    (X ++ Y).idxOf a < X.length := by
  rw [List.idxOf_append_of_mem h]; exact List.idxOf_lt_length_iff.2 h

theorem idxOf_append_of_not_mem_left {a : Nat} {X Y : List Nat} (h : a ∉ X) :
    (X ++ Y).idxOf a = X.length + Y.idxOf a := by
  rw [List.idxOf_append, if_neg h, Nat.add_comm]

theorem Precedes.insert {A P M Q B : List Nat} {a b : Nat}
    (hP : ∀ x ∈ P, x ∉ A ++ M ++ B) (hQ : ∀ x ∈ Q, x ∉ A ++ M ++ B)
    (h : Precedes (A ++ M ++ B) a b) : Precedes (A ++ P ++ M ++ Q ++ B) a b := by
  obtain ⟨ha, hb, hlt⟩ := h
  have mem : ∀ x ∈ A ++ M ++ B, x ∈ A ++ P ++ M ++ Q ++ B := fun x hx => by
    simp only [List.mem_append] at hx ⊢; tauto
  have idx : ∀ x ∈ A ++ M ++ B, (A ++ P ++ M ++ Q ++ B).idxOf x =
      (A ++ M ++ B).idxOf x + (if x ∈ A then 0 else if x ∈ M then P.length else P.length + Q.length) := by
    intro x hx
    have hxP : x ∉ P := fun h => hP x h hx
    have hxQ : x ∉ Q := fun h => hQ x h hx
    by_cases hA : x ∈ A
    · simp only [List.append_assoc] at *
      rw [List.idxOf_append_of_mem hA, List.idxOf_append_of_mem hA, if_pos hA]; rfl
    · simp only [List.append_assoc] at *
      rw [idxOf_append_of_not_mem_left hA, idxOf_append_of_not_mem_left hA, idxOf_append_of_not_mem_left hxP,
        if_neg hA]
      by_cases hM : x ∈ M
      · rw [List.idxOf_append_of_mem hM, List.idxOf_append_of_mem hM, if_pos hM]; omega
      · rw [idxOf_append_of_not_mem_left hM, idxOf_append_of_not_mem_left hM, idxOf_append_of_not_mem_left hxQ,
          if_neg hM]; omega
  refine ⟨mem a ha, mem b hb, ?_⟩
  rw [idx a ha, idx b hb]
  have hpos : ∀ x ∈ A ++ M ++ B, (x ∉ A → A.length ≤ (A ++ M ++ B).idxOf x) ∧
      (x ∈ A → (A ++ M ++ B).idxOf x < A.length) ∧
      (x ∉ A → x ∉ M → A.length + M.length ≤ (A ++ M ++ B).idxOf x) ∧
      (x ∉ A → x ∈ M → (A ++ M ++ B).idxOf x < A.length + M.length) := by
    intro x hx
    simp only [List.append_assoc]
    refine ⟨fun hA => ?_, fun hA => ?_, fun hA hM => ?_, fun hA hM => ?_⟩
    · rw [idxOf_append_of_not_mem_left hA]; omega
    · exact idxOf_lt_length_of_mem_left hA
    · rw [idxOf_append_of_not_mem_left hA, idxOf_append_of_not_mem_left hM]; omega
    · rw [idxOf_append_of_not_mem_left hA]
      have := idxOf_lt_length_of_mem_left (Y := B) hM; omega
  obtain ⟨a1, a2, a3, a4⟩ := hpos a ha
  obtain ⟨b1, b2, b3, b4⟩ := hpos b hb
  by_cases haA : a ∈ A
  · rw [if_pos haA]; have := a2 haA
    by_cases hbA : b ∈ A
    · rw [if_pos hbA]; omega
    · rw [if_neg hbA]; have := b1 hbA; split_ifs <;> omega
  · rw [if_neg haA]; have := a1 haA
    by_cases hbA : b ∈ A
    · rw [if_pos hbA]; have := b2 hbA; omega
    · rw [if_neg hbA]; have := b1 hbA
      by_cases haM : a ∈ M
      · rw [if_pos haM]; have := a4 haA haM
        by_cases hbM : b ∈ M
        · rw [if_pos hbM]; omega
        · rw [if_neg hbM]; have := b3 hbA hbM; omega
      · rw [if_neg haM]; have := a3 haA haM
        by_cases hbM : b ∈ M
        · rw [if_pos hbM]; have := b4 hbA hbM; omega
        · rw [if_neg hbM]; omega

theorem Oriented.insert {A P M Q B : List Nat} {vs : Option Nat × Option Nat}
    (hP : ∀ x ∈ P, x ∉ A ++ M ++ B) (hQ : ∀ x ∈ Q, x ∉ A ++ M ++ B)
    (h : Oriented (A ++ M ++ B) vs) : Oriented (A ++ P ++ M ++ Q ++ B) vs := by
  obtain ⟨a | _, b | _⟩ := vs <;> simp only [Oriented] at h ⊢
  exact h.insert hP hQ

/-! ### The vertex sequence of a grown block -/

/-- The vertices of an item list. -/
def seqOf (g : Graph) (L : List ItemId) : List Nat :=
  L.filterMap fun x => if 1 ≤ x ∧ x < 1 + g.nv then some (x - 1) else none

theorem seqOf_append (g : Graph) (L L' : List ItemId) : seqOf g (L ++ L') = seqOf g L ++ seqOf g L' :=
  List.filterMap_append

theorem StBlock.seq_seqOf (g : Graph) (r : Option (Nat × Nat)) (L : List ItemId) :
    StBlock.seq g ⟨r, L⟩ = (r.map (·.1)).toList ++ seqOf g L := rfl

theorem mem_seqOf {g : Graph} {L : List ItemId} {y : Nat} :
    y ∈ seqOf g L ↔ vertItem y ∈ L ∧ y < g.nv := by
  simp only [seqOf, List.mem_filterMap, Option.ite_none_right_eq_some, Option.some.injEq, vertItem]
  constructor
  · rintro ⟨x, hx, ⟨h1, h2⟩, rfl⟩
    refine ⟨?_, by omega⟩
    have : 1 + (x - 1) = x := by omega
    rwa [this]
  · rintro ⟨h, hy⟩
    exact ⟨1 + y, h, ⟨by omega, by omega⟩, by omega⟩

theorem seqOf_fresh {g : Graph} {P A : List ItemId} (hP : ∀ x ∈ P, x ∉ A) :
    ∀ y ∈ seqOf g P, y ∉ seqOf g A := fun y hy hA =>
  hP _ (mem_seqOf.1 hy).1 (mem_seqOf.1 hA).1

theorem seq_insert (g : Graph) (r : Option (Nat × Nat)) (A P M Q B : List ItemId) :
    StBlock.seq g ⟨r, A ++ P ++ M ++ Q ++ B⟩ =
      ((r.map (·.1)).toList ++ seqOf g A) ++ seqOf g P ++ seqOf g M ++ seqOf g Q ++ seqOf g B ∧
    StBlock.seq g ⟨r, A ++ M ++ B⟩ = ((r.map (·.1)).toList ++ seqOf g A) ++ seqOf g M ++ seqOf g B := by
  simp [StBlock.seq_seqOf, seqOf_append, List.append_assoc]

/-- Freshness of inserted items, at the level of the vertex sequence. -/
theorem seqOf_fresh_root {g : Graph} {r : Option (Nat × Nat)} {P A M B : List ItemId}
    (hP : ∀ x ∈ P, x ∉ A ++ M ++ B) (hr : ∀ x ∈ P, ∀ a, r = some a → x ≠ vertItem a.1) :
    ∀ y ∈ seqOf g P, y ∉ ((r.map (·.1)).toList ++ seqOf g A) ++ seqOf g M ++ seqOf g B := by
  intro y hy h
  obtain ⟨hyP, -⟩ := mem_seqOf.1 hy
  simp only [List.mem_append, Option.mem_toList, Option.map_eq_some_iff] at h
  rcases h with ((⟨a, ha, hay⟩ | hA) | hM) | hB
  · exact hr _ hyP a ha (by rw [hay])
  · exact hP _ hyP (List.mem_append_left _ (List.mem_append_left _ (mem_seqOf.1 hA).1))
  · exact hP _ hyP (List.mem_append_left _ (List.mem_append_right _ (mem_seqOf.1 hM).1))
  · exact hP _ hyP (List.mem_append_right _ (mem_seqOf.1 hB).1)

/-! ### `InBlock` / `StLive` under insertion -/

theorem ExpandsList.below_of_mem {items : Items} {xs L : List ItemId} (h : ExpandsList items xs L)
    {y : ItemId} (hy : y ∈ L) : ∃ x ∈ xs, Items.Below items x y := by
  induction h with
  | nil => exact absurd hy List.not_mem_nil
  | leaf _ _ ih =>
    rcases List.mem_cons.1 hy with rfl | hy
    · exact ⟨_, List.mem_cons_self, .refl⟩
    · obtain ⟨x, hx, hb⟩ := ih hy
      exact ⟨x, List.mem_cons_of_mem _ hx, hb⟩
  | node _ _ ih =>
    obtain ⟨x, hx, hb⟩ := ih hy
    rcases List.mem_append.1 hx with hx | hx
    · exact ⟨_, List.mem_cons_self, .head hx hb⟩
    · exact ⟨x, List.mem_cons_of_mem _ hx, hb⟩

theorem Expands.below_of_mem {items : Items} {i : ItemId} {L : List ItemId} (h : Expands items i L)
    {y : ItemId} (hy : y ∈ L) : Items.Below items i y := by
  obtain ⟨x, hx, hb⟩ := ExpandsList.below_of_mem h hy
  rw [List.mem_singleton.1 hx] at hb; exact hb

/-- Growing the block by fresh items (not below `i`) keeps `i` in it, provided its leaves stay a
segment. -/
theorem InBlock.insert {g : Graph} {items : Items} {r : Option (Nat × Nat)} {A P M Q B : List ItemId}
    {i : ItemId}
    (hP : ∀ x ∈ P ++ Q, x ∉ A ++ M ++ B) (hr : ∀ x ∈ P ++ Q, ∀ a, r = some a → x ≠ vertItem a.1)
    (hbelow : ∀ y, Items.Below items i y → y ∉ P ++ Q)
    (h : InBlock g items ⟨r, A ++ M ++ B⟩ i) {L : List ItemId} (hL : Expands items i L)
    (hseg : ∃ A' B', A ++ P ++ M ++ Q ++ B = A' ++ L ++ B') :
    InBlock g items ⟨r, A ++ P ++ M ++ Q ++ B⟩ i := by
  obtain ⟨L', hL', -, h1, h2, h3, h4, h5, h6⟩ := h
  have hLe : L = L' := ExpandsList.unique hL hL'
  subst hLe
  obtain ⟨hs1, hs2⟩ := seq_insert g r A P M Q B
  have hP' : ∀ x ∈ P, x ∉ A ++ M ++ B := fun x hx => hP x (List.mem_append_left _ hx)
  have hQ' : ∀ x ∈ Q, x ∉ A ++ M ++ B := fun x hx => hP x (List.mem_append_right _ hx)
  have fP := seqOf_fresh_root (g := g) hP' fun x hx => hr x (List.mem_append_left _ hx)
  have fQ := seqOf_fresh_root (g := g) hQ' fun x hx => hr x (List.mem_append_right _ hx)
  have mem : ∀ x ∈ A ++ M ++ B, x ∈ A ++ P ++ M ++ Q ++ B := fun x hx => by
    simp only [List.mem_append] at hx ⊢; tauto
  have old : ∀ e, edgeItem g e ∈ A ++ P ++ M ++ Q ++ B → Items.Below items i (edgeItem g e) →
      edgeItem g e ∈ A ++ M ++ B := by
    intro e he hb
    have := hbelow _ hb
    simp only [List.mem_append] at he this ⊢; tauto
  have pre : ∀ a b, Precedes (StBlock.seq g ⟨r, A ++ M ++ B⟩) a b →
      Precedes (StBlock.seq g ⟨r, A ++ P ++ M ++ Q ++ B⟩) a b := by
    intro a b hab; rw [hs1]; rw [hs2] at hab; exact hab.insert fP fQ
  refine ⟨L, hL, hseg, fun x hx => mem x (h1 x hx), ?_, ?_, ?_, ?_, ?_⟩
  · intro e he hb hbe
    exact h2 e he (hb.imp_left (old e · hbe)) hbe
  · show Oriented _ _
    rw [hs1]; rw [hs2] at h3; exact h3.insert fP fQ
  · intro s t hst c hc hV
    exact ⟨pre _ _ (h4 s t hst c hc hV).1, pre _ _ (h4 s t hst c hc hV).2⟩
  · intro c hc hV
    have := h5 c hc hV
    show Oriented _ _
    rw [hs1]; rw [hs2] at this; exact this.insert fP fQ
  · intro c hc hV u v huv e he hb hbe y hy
    have hbe' : Items.Below items i (edgeItem g e) := Relation.ReflTransGen.head hc hbe
    obtain ⟨l, r'⟩ := h6 c hc hV u v huv e he (hb.imp_left (old e · hbe')) hbe y hy
    exact ⟨l.imp_right (pre _ _), r'.imp_right (pre _ _)⟩

/-! ### Frames -/

/-- The items below the span items of `new` keep type, children and `vs`. -/
def ItemsFrame (items items' : Items) (new : List TEntry) : Prop :=
  ∀ x ∈ readStack new, ∀ y, Items.Below items x y →
    Items.type items' y = Items.type items y ∧ Items.ch items' y = Items.ch items y ∧
      Items.vs items' y = Items.vs items y

theorem Items.Below.of_ch_frame {items items' : Items} {x : ItemId}
    (hf : ∀ y, Items.Below items x y → Items.ch items' y = Items.ch items y) :
    ∀ {z : ItemId}, Items.Below items x z → Items.Below items' x z := by
  intro z h
  induction h with
  | refl => exact .refl
  | tail hxy hyz ih => exact .tail ih (show _ ∈ Items.ch items' _ by rw [hf _ hxy]; exact hyz)

theorem Items.Below.to_of_ch_frame {items items' : Items} {x : ItemId}
    (hf : ∀ y, Items.Below items x y → Items.ch items' y = Items.ch items y) :
    ∀ {z : ItemId}, Items.Below items' x z → Items.Below items x z := by
  intro z h
  induction h with
  | refl => exact .refl
  | tail _ hyz ih => exact .tail ih (show _ ∈ Items.ch items _ by rw [← hf _ ih]; exact hyz)

theorem HangingUnder.iff_of_frame {items items' : Items} {x i : ItemId}
    (hf : ∀ y, Items.Below items x y →
      Items.type items' y = Items.type items y ∧ Items.ch items' y = Items.ch items y) :
    HangingUnder items' x i ↔ HangingUnder items x i := by
  have hch : ∀ y, Items.Below items x y → Items.ch items' y = Items.ch items y := fun y hy => (hf y hy).2
  constructor
  · rintro ⟨y, h1, h2, h3, h4⟩
    have h1' := Items.Below.to_of_ch_frame hch h1
    exact ⟨y, h1', h2, by rwa [(hf y h1').1] at h3,
      Items.Below.to_of_ch_frame (fun z hz => hch z (h1'.trans hz)) h4⟩
  · rintro ⟨y, h1, h2, h3, h4⟩
    exact ⟨y, Items.Below.of_ch_frame hch h1, h2, by rwa [(hf y h1).1],
      Items.Below.of_ch_frame (fun z hz => hch z (h1.trans hz)) h4⟩

theorem StLive.frame {g : Graph} {items items' : Items} {new : List TEntry} {b : StBlock}
    (hf : ItemsFrame items items' new) (h : StLive g items new b) : StLive g items' new b := by
  intro x hx i hxi' hty' hnh'
  have hf' := hf x hx
  have hch : ∀ y, Items.Below items x y → Items.ch items' y = Items.ch items y :=
    fun y hy => (hf' y hy).2.1
  have hxi := Items.Below.to_of_ch_frame hch hxi'
  have hty : Items.type items i = .S ∨ Items.type items i = .P ∨ Items.type items i = .R := by
    rwa [(hf' i hxi).1] at hty'
  have hnh : ¬ HangingUnder items x i := fun hh =>
    hnh' ((HangingUnder.iff_of_frame fun y hy => ⟨(hf' y hy).1, (hf' y hy).2.1⟩).2 hh)
  refine InBlock.congr (fun y hy => ⟨(hf' y (hxi.trans hy)).1, (hf' y (hxi.trans hy)).2.1⟩)
    (hf' i hxi).2.2 (fun c hc => (hf' c (hxi.trans (.single hc))).2.2) ?_ ?_ (h x hx i hxi hty hnh)
  · intro e he
    exact Items.Below.to_of_ch_frame (fun z hz => hch z (hxi.trans hz)) he
  · intro c hc e he _
    exact Items.Below.to_of_ch_frame
      (fun z hz => hch z (hxi.trans ((Relation.ReflTransGen.single hc).trans hz))) he

theorem StLive.nil {g : Graph} {items : Items} {b : StBlock} : StLive g items [] b := by
  intro x hx
  simp [readStack, readL, readR] at hx

theorem StLive.of_subset {g : Graph} {items : Items} {new new' : List TEntry} {b : StBlock}
    (hs : ∀ x ∈ readStack new', x ∈ readStack new) (h : StLive g items new b) :
    StLive g items new' b := fun x hx => h x (hs x hx)

theorem StLive.append {g : Graph} {items : Items} {new new' : List TEntry} {b : StBlock}
    (h : StLive g items new b) (h' : StLive g items new' b) : StLive g items (new ++ new') b := by
  intro x hx
  rcases mem_readStack_append.1 hx with hx | hx
  · exact h x hx
  · exact h' x hx

/-- Inserting fresh items around the pieces `ps` the segment `new` reads as keeps it live: the
leaves of a live item are a segment of `stNest ps` (`expands_cases`, not hanging), hence of the
grown list. -/
theorem StLive.insert {g : Graph} {items : Items} {new : List TEntry} {ps : List StPiece}
    {r : Option (Nat × Nat)} {A P M Q B : List ItemId}
    (hR : StRead items new ps) (hM : M = stNest ps)
    (hP : ∀ x ∈ P ++ Q, x ∉ A ++ M ++ B) (hr : ∀ x ∈ P ++ Q, ∀ a, r = some a → x ≠ vertItem a.1)
    (hbelow : ∀ x ∈ readStack new, ∀ y, Items.Below items x y → y ∉ P ++ Q)
    (h : StLive g items new ⟨r, A ++ M ++ B⟩) : StLive g items new ⟨r, A ++ P ++ M ++ Q ++ B⟩ := by
  intro x hx i hxi hty hnh
  have hI := h x hx i hxi hty hnh
  have hI' := hI
  obtain ⟨L, hL, -, -⟩ := hI'
  have key : ∀ {xs N : List ItemId}, ExpandsList items xs N → x ∈ xs → N <:+: M → L <:+: M := by
    intro xs N hN hxN hNM
    rcases hxi.expands_cases hN hxN with hh | ⟨L', hL', A', B', hEq⟩
    · exact absurd hh hnh
    · rw [ExpandsList.unique hL hL']
      exact List.IsInfix.trans (show L' <:+: N from ⟨A', B', hEq.symm⟩) hNM
  have hLM : L <:+: M := by
    rcases List.mem_append.1 (show x ∈ readL new ++ readR new from hx) with hx' | hx'
    · exact key hR.1 hx' (by rw [hM]; exact (List.prefix_append (stNestL ps) (stNestR ps)).isInfix)
    · exact key hR.2 hx' (by rw [hM]; exact (List.suffix_append (stNestL ps) (stNestR ps)).isInfix)
  exact hI.insert hP hr (fun y hy => hbelow x hx y (hxi.trans hy)) hL
    (segment_insert (Or.inr (Or.inl hLM)))

end Spqr

namespace Spqr

/-- A stack segment whose span items are all `V`/`Q` is live in any block: every S/P/R item below
such an item hangs under it. -/
theorem StLive.of_leafTypes {g : Graph} {items : Items} {new : List TEntry} {b : StBlock}
    (h : ∀ x ∈ readStack new, Items.type items x = .V ∨ Items.type items x = .Q) :
    StLive g items new b := by
  intro x hx i hxi hty hnh
  exfalso
  refine hnh ⟨x, .refl, fun e => ?_, h x hx, hxi⟩
  subst e
  rcases h x hx with h' | h' <;> rcases hty with hty | hty | hty <;> rw [h'] at hty <;> cases hty

theorem StLive.cons {g : Graph} {items : Items} {new : List TEntry} {t : TEntry} {b : StBlock}
    (ht : StLive g items [t] b) (h : StLive g items new b) : StLive g items (t :: new) b := by
  rw [← List.singleton_append]; exact ht.append h

end Spqr
