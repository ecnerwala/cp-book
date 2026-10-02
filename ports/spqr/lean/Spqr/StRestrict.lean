import Spqr.StRefEt
import Spqr.ItemTree

/-!
# Restricting a block's st-order to an item

`stItem_of_refOrder`: an S/P/R item of the walk whose children are in reference order
(`walk_st'`) and which is placed in a reference block as `VsOriented` says is `Items.StItem`,
given that the block is st-numbered (`refBlocks_st`). The item's vertex list is the block order
restricted to the item's vertices; a lower (higher) block neighbour of an interior vertex is found
under one child of the item, and that child's virtual edge is the item-level neighbour. See
`PROOF.md` §7.6.
-/

namespace Spqr.StRefEt

/-! ### `Before`, `collapseRuns` -/

namespace Before

variable {α : Type} {L : List α} {a b : α}

theorem or_of_mem_mem (ha : a ∈ L) (hb : b ∈ L) (hne : a ≠ b) : Before L a b ∨ Before L b a := by
  obtain ⟨l₁, l₂, rfl⟩ := List.append_of_mem ha
  rcases List.mem_append.mp hb with hb | hb
  · exact Or.inr (of_mem_mem hb (List.mem_cons_self ..))
  · rcases List.mem_cons.mp hb with rfl | hb
    · exact absurd rfl hne
    · obtain ⟨s, t, rfl⟩ := List.append_of_mem hb
      exact Or.inl ⟨l₁, s, t, by simp⟩

theorem filter {p : α → Bool} (h : Before (L.filter p) a b) : Before L a b :=
  h.sub List.filter_sublist

theorem cons (h : Before L a b) (x : α) : Before (x :: L) a b :=
  h.sub (List.sublist_cons_self x L)

theorem append_left (h : Before L a b) (L' : List α) : Before (L' ++ L) a b :=
  h.sub (List.sublist_append_right L' L)

theorem append_right (h : Before L a b) (L' : List α) : Before (L ++ L') a b :=
  h.sub (List.sublist_append_left L L')

theorem cons_of_mem (hb : b ∈ L) : Before (a :: L) a b :=
  of_mem_mem (List.mem_singleton_self a) hb

theorem filter_of {p : α → Bool} (h : Before L a b) (ha : p a = true) (hb : p b = true) :
    Before (L.filter p) a b := by
  obtain ⟨l₁, l₂, l₃, rfl⟩ := h
  exact ⟨l₁.filter p, l₂.filter p, l₃.filter p, by
    simp [List.filter_append, ha, hb]⟩

theorem map {β : Type} (f : α → β) (h : Before L a b) : Before (L.map f) (f a) (f b) := by
  obtain ⟨l₁, l₂, l₃, rfl⟩ := h
  exact ⟨l₁.map f, l₂.map f, l₃.map f, by simp⟩

theorem of_cons [DecidableEq α] (h : Before (a :: L) a b) (hnd : (a :: L).Nodup) : b ∈ L := by
  have hb := h.mem_right
  rcases List.mem_cons.mp hb with rfl | hb
  · exact absurd rfl (h.ne hnd)
  · exact hb

end Before

theorem mem_collapseRuns : ∀ (L : List ItemId) (a : ItemId), a ∈ collapseRuns L ↔ a ∈ L
  | [], a => by simp [collapseRuns]
  | [x], a => by simp [collapseRuns]
  | x :: y :: xs, a => by
    rw [collapseRuns]
    split_ifs with h
    · rw [mem_collapseRuns (y :: xs)]; simp [h]
    · simp [mem_collapseRuns]

theorem collapseRuns_sublist : ∀ (L : List ItemId), List.Sublist (collapseRuns L) L
  | [] => by simp [collapseRuns]
  | [x] => by simp [collapseRuns]
  | x :: y :: xs => by
    rw [collapseRuns]
    split_ifs with h
    · exact (collapseRuns_sublist (y :: xs)).cons x
    · exact (collapseRuns_sublist (y :: xs)).cons_cons x

/-- Collapsing runs keeps the order of distinct elements. -/
theorem before_collapseRuns : ∀ (L : List ItemId) (a b : ItemId), Before L a b → a ≠ b →
    Before (collapseRuns L) a b
  | [], a, b, h, _ => absurd h.mem_left List.not_mem_nil
  | [x], a, b, h, _ => by
    obtain ⟨l₁, l₂, l₃, hL⟩ := h
    have := congrArg List.length hL; simp at this; omega
  | x :: y :: xs, a, b, h, hne => by
    rw [collapseRuns]
    obtain ⟨l₁, l₂, l₃, hL⟩ := h
    split_ifs with hxy
    · apply before_collapseRuns _ _ _ _ hne
      cases l₁ with
      | nil =>
        simp only [List.nil_append, List.cons_append, List.cons.injEq] at hL
        obtain ⟨rfl, hL⟩ := hL
        have hb : b ∈ y :: xs := by rw [hL]; simp
        rw [hxy]
        rcases List.mem_cons.mp hb with rfl | hb
        · exact absurd hxy hne
        · exact Before.cons_of_mem hb
      | cons z l₁ =>
        simp only [List.cons_append, List.cons.injEq, List.append_assoc] at hL
        exact ⟨l₁, l₂, l₃, by simpa using hL.2⟩
    · cases l₁ with
      | nil =>
        simp only [List.nil_append, List.cons_append, List.cons.injEq] at hL
        obtain ⟨rfl, hL⟩ := hL
        exact Before.cons_of_mem ((mem_collapseRuns _ _).mpr (by rw [hL]; simp))
      | cons z l₁ =>
        simp only [List.cons_append, List.cons.injEq, List.append_assoc] at hL
        exact (before_collapseRuns _ _ _ ⟨l₁, l₂, l₃, by simpa using hL.2⟩ hne).cons x

end Spqr.StRefEt

namespace Spqr

open StRefEt

/-! ### Leaves and the item tree -/

namespace Items

variable {g : Graph} {items : Items}

theorem below_of_mem_leaves : ∀ (fuel : Nat) (i x : ItemId), x ∈ Items.leaves items fuel i →
    items.Below i x
  | 0, i, x, h => by
    simp [Items.leaves] at h; subst h; exact Relation.ReflTransGen.refl
  | fuel + 1, i, x, h => by
    unfold Items.leaves at h
    split at h
    · simp at h; subst h; exact Relation.ReflTransGen.refl
    · simp at h; subst h; exact Relation.ReflTransGen.refl
    · obtain ⟨c, hc, hx⟩ := List.mem_flatMap.mp h
      exact Relation.ReflTransGen.head hc (below_of_mem_leaves fuel c x hx)

theorem leaves_of_type (fuel : Nat) (i : ItemId) (h : items.type i = .V ∨ items.type i = .Q) :
    Items.leaves items fuel i = [i] := by
  cases fuel with
  | zero => rfl
  | succ fuel => unfold Items.leaves; rcases h with h | h <;> rw [h]

theorem leaves_succ (fuel : Nat) (i : ItemId)
    (h : items.type i = .S ∨ items.type i = .P ∨ items.type i = .R) :
    Items.leaves items (fuel + 1) i = (items.ch i).flatMap (Items.leaves items fuel) := by
  unfold Items.leaves; rcases h with h | h | h <;> rw [h]

theorem Tree.child_below_child (ht : items.Tree g) {a c₁ c₂ : ItemId} (h1 : items.IsParent a c₁)
    (h2 : items.IsParent a c₂) (h12 : items.Below c₁ c₂) : c₁ = c₂ := by
  rcases Relation.ReflTransGen.cases_tail h12 with rfl | ⟨p, h1p, hp2⟩
  · rfl
  · obtain rfl := ht.parent_eq hp2 h2
    exact absurd (Relation.TransGen.head' h1 h1p) (ht.acyclic _)

theorem find?_leaves_vertItem (ht : items.Tree g) {i : ItemId} {y : Nat} (fuel : Nat)
    (hy : items.IsParent i (vertItem y)) (hV : items.type (vertItem y) = .V) :
    (items.ch i).find? (fun c => vertItem y ∈ Items.leaves items fuel c) = some (vertItem y) := by
  have hsome : ((items.ch i).find? (fun c => vertItem y ∈ Items.leaves items fuel c)).isSome := by
    rw [List.find?_isSome]
    exact ⟨vertItem y, hy, by simp [leaves_of_type fuel _ (Or.inl hV)]⟩
  obtain ⟨c, hc⟩ := Option.isSome_iff_exists.mp hsome
  have hmem := List.mem_of_find?_eq_some hc
  have hp := List.find?_some hc
  simp only [decide_eq_true_eq] at hp
  rw [hc, ht.child_below_child hmem hy (below_of_mem_leaves fuel c _ hp)]

/-- `walk_st'` transports the order of two V children from the reference order to `ch`. -/
theorem before_ch_of_before_order (ht : items.Tree g) {i : ItemId} {order : List ItemId}
    {fuel : Nat} (hch : items.ch i = restrictCh items fuel order i) {y y' : Nat} (hne : y ≠ y')
    (hy : items.IsParent i (vertItem y)) (hy' : items.IsParent i (vertItem y'))
    (hV : items.type (vertItem y) = .V) (hV' : items.type (vertItem y') = .V)
    (h : Before order (vertItem y) (vertItem y')) :
    Before (items.ch i) (vertItem y) (vertItem y') := by
  rw [hch]
  unfold restrictCh
  exact before_collapseRuns _ _ _
    (h.filterMap (find?_leaves_vertItem ht fuel hy hV) (find?_leaves_vertItem ht fuel hy' hV'))
    (fun h => hne (vertItem_inj h))

end Items

theorem Precedes.of_before {xs : List Nat} {a b : Nat} (hnd : xs.Nodup) (h : Before xs a b) :
    Precedes xs a b :=
  ⟨h.mem_left, h.mem_right, h.idxOf_lt hnd⟩

theorem Precedes.before {xs : List Nat} {a b : Nat} (hnd : xs.Nodup) (h : Precedes xs a b) :
    Before xs a b := by
  obtain ⟨ha, hb, hlt⟩ := h
  rcases Before.or_of_mem_mem ha hb (fun h => by subst h; exact lt_irrefl _ hlt) with h | h
  · exact h
  · have := h.idxOf_lt hnd; omega

theorem Precedes.ne {xs : List Nat} {a b : Nat} (h : Precedes xs a b) : a ≠ b := by
  rintro rfl; exact lt_irrefl _ h.2.2

theorem Precedes.asymm {xs : List Nat} {a b : Nat} (h : Precedes xs a b) (h' : Precedes xs b a) :
    False := by
  have := h.2.2; have := h'.2.2; omega

theorem Precedes.trans {xs : List Nat} {a b c : Nat} (h : Precedes xs a b) (h' : Precedes xs b c) :
    Precedes xs a c :=
  ⟨h.1, h'.2.1, lt_trans h.2.2 h'.2.2⟩

theorem idxOf_eq_zero_of_head? {L : List Nat} {x : Nat} (h : L.head? = some x) : L.idxOf x = 0 := by
  cases L with
  | nil => cases h
  | cons y L => simp at h; subst h; simp

theorem idxOf_le_of_getLast? {L : List Nat} {x : Nat} (hnd : L.Nodup) (h : L.getLast? = some x)
    {y : Nat} (hy : y ∈ L) : L.idxOf y ≤ L.idxOf x := by
  obtain ⟨L', rfl⟩ := List.getLast?_eq_some_iff.mp h
  have hx : x ∉ L' := (List.nodup_cons.mp (List.nodup_append_comm.mp hnd)).1
  rw [List.idxOf_append_of_notMem hx, List.idxOf_cons_self]
  have := List.idxOf_lt_length_of_mem hy
  simp at this; omega

theorem vertItem_sub_one (y : Nat) : vertItem y - 1 = y := by
  show 1 + y - 1 = y; omega

theorem edgeOf_eq_some {g : Graph} {y : ItemId} {p : Nat × Nat} (h : edgeOf g y = some p) :
    ∃ e, e < g.ne ∧ y = edgeItem g e ∧ g.edges[e]! = p := by
  obtain ⟨n, rfl⟩ : ∃ n : Nat, n = y := ⟨y, rfl⟩
  unfold edgeOf at h
  split_ifs at h with hc
  · have h1 : 1 + g.nv ≤ n := hc.1
    have h2 : n < 1 + g.nv + g.ne := hc.2
    have h3 : n - (1 + g.nv) < g.ne := by omega
    have h4 : n = edgeItem g (n - (1 + g.nv)) := by simp only [edgeItem]; omega
    exact ⟨n - (1 + g.nv), h3, h4, by simpa using h⟩

namespace Items

variable {g : Graph} {items : Items}

theorem Tree.child_class (ht : items.Tree g) {p c : ItemId} (hp : items.IsParent p c) :
    (∃ y, y < g.nv ∧ c = vertItem y) ∨ (∃ e, e < g.ne ∧ c = edgeItem g e) ∨
      (1 + g.nv + g.ne ≤ c ∧ c < items.size) := by
  have hc := ht.ch_lt p c hp
  obtain ⟨n, rfl⟩ : ∃ n : Nat, n = c := ⟨c, rfl⟩
  rcases Nat.eq_zero_or_pos n with h0 | h0
  · subst h0; exact absurd hp (ht.root_no_parent p)
  rcases Nat.lt_or_ge n (1 + g.nv) with h1 | h1
  · have h3 : n - 1 < g.nv := by omega
    have h4 : n = vertItem (n - 1) := by simp only [vertItem]; omega
    exact Or.inl ⟨n - 1, h3, h4⟩
  rcases Nat.lt_or_ge n (1 + g.nv + g.ne) with h2 | h2
  · have h3 : n - (1 + g.nv) < g.ne := by omega
    have h4 : n = edgeItem g (n - (1 + g.nv)) := by simp only [edgeItem]; omega
    exact Or.inr (Or.inl ⟨n - (1 + g.nv), h3, h4⟩)
  · exact Or.inr (Or.inr ⟨h2, hc⟩)

theorem Tree.v_child_eq (ht : items.Tree g) {p c : ItemId} (hp : items.IsParent p c)
    (hV : items.type c = .V) : ∃ y, y < g.nv ∧ c = vertItem y := by
  rcases ht.child_class hp with h | ⟨e, he, rfl⟩ | ⟨h1, h2⟩
  · exact h
  · rw [ht.edge e he] at hV; cases hV
  · have := ht.node c h1 h2; rw [hV] at this; simp at this

theorem Tree.type_ne_F (ht : items.Tree g) {p c : ItemId} (hp : items.IsParent p c) :
    items.type c ≠ .F := by
  rcases ht.child_class hp with ⟨y, hy, rfl⟩ | ⟨e, he, rfl⟩ | ⟨h1, h2⟩
  · rw [ht.vert y hy]; simp
  · rw [ht.edge e he]; simp
  · have := ht.node c h1 h2; intro h; rw [h] at this; simp at this

end Items

/-- The restriction argument for one block `b` of the reference order: an S/P/R item whose
children are in `b`'s order, whose leaves lie in `b`, and which is oriented along `b.seq`
(the `VsOriented` clauses) inherits `Items.StItem` from `b.St`. -/
theorem stItem_of_block (g : Graph) (items : Items)
    (hwf : Items.WF g items) (i : ItemId) (hi : i < items.size)
    (ht : items.type i = .S ∨ items.type i = .P ∨ items.type i = .R)
    (b : StBlock) (order : List ItemId) (L₁ L₂ : List ItemId)
    (hord : order = L₁ ++ b.items ++ L₂)
    (hch : items.ch i = restrictCh items items.size order i)
    (hb : b.St g)
    (hre : ∀ r ∈ b.root, ∃ e, e < g.ne ∧ Items.PairEq g.edges[e]! r)
    (hleaf : ∀ x ∈ Items.leaves items items.size i, x ∈ b.items)
    (hedge : ∀ e, e < g.ne →
      (edgeItem g e ∈ b.items ∨ ∃ r ∈ b.root, Items.PairEq g.edges[e]! r) →
      items.Below i (edgeItem g e) → edgeItem g e ∈ Items.leaves items items.size i)
    (hor : Oriented (b.seq g) (items.vs i))
    (hbet : ∀ s t, items.vs i = (some s, some t) → ∀ c ∈ items.ch i, items.type c = .V →
      Precedes (b.seq g) s (c - 1) ∧ Precedes (b.seq g) (c - 1) t)
    (hcor : ∀ c ∈ items.ch i, items.type c ≠ .V → Oriented (b.seq g) (items.vs c))
    (hspan : ∀ c ∈ items.ch i, items.type c ≠ .V → ∀ u v, items.vs c = (some u, some v) →
      ∀ e, e < g.ne → (edgeItem g e ∈ b.items ∨ ∃ r ∈ b.root, Items.PairEq g.edges[e]! r) →
      items.Below c (edgeItem g e) →
      ∀ y, (g.edges[e]!).1 = y ∨ (g.edges[e]!).2 = y →
        (u = y ∨ Precedes (b.seq g) u y) ∧ (y = v ∨ Precedes (b.seq g) y v)) :
    Items.StItem items i := by
  set sq := b.seq g with hsq
  have ht' : items.type i ∉ [NodeType.F, .V] := by rcases ht with h | h | h <;> simp [h]
  obtain ⟨s, t, hst⟩ : ∃ s t, items.vs i = (some s, some t) := by
    have := hwf.endpoints.vs_shape i hi
    rcases ht with h | h | h <;> rw [h] at this <;> exact this
  have hsqnd : sq.Nodup := hb.1
  have hsqed := hb.2.1
  have hsqint := hb.2.2
  have hst' : Precedes sq s t := by rw [hst] at hor; exact hor
  obtain ⟨n, hn⟩ : ∃ n, items.size = n + 1 := ⟨items.size - 1, by have := hwf.tree.size; omega⟩
  have hleaves : Items.leaves items items.size i = (items.ch i).flatMap (Items.leaves items n) := by
    rw [hn]; exact Items.leaves_succ n i ht
  have hchild_leaf : ∀ c ∈ items.ch i, (items.type c = .V ∨ items.type c = .Q) →
      c ∈ Items.leaves items items.size i := by
    intro c hc hcT
    rw [hleaves]
    exact List.mem_flatMap.mpr ⟨c, hc, by rw [Items.leaves_of_type n c hcT]; simp⟩
  -- the vertex list
  set xs := ((items.ch i).filter fun c => items.type c = .V).map (· - 1) with hxs
  have hV : items.vertList i = s :: (xs ++ [t]) := by
    simp [Items.vertList, hst, hxs]
  have hVnd : (items.vertList i).Nodup := hwf.endpoints.nv_nodup i hi
  rw [hV] at hVnd
  set V := s :: (xs ++ [t]) with hVdef
  have hUmem : ∀ y, y ∈ V ↔ y = s ∨ y ∈ xs ∨ y = t := by intro y; simp [hVdef]
  have hxs_iff : ∀ x, x ∈ xs ↔ x < g.nv ∧ items.IsParent i (vertItem x) := by
    intro x
    simp only [hxs, List.mem_map, List.mem_filter, decide_eq_true_eq]
    constructor
    · rintro ⟨c, ⟨hc, hcV⟩, rfl⟩
      obtain ⟨y, hy, rfl⟩ := hwf.tree.v_child_eq hc hcV
      rw [vertItem_sub_one]; exact ⟨hy, hc⟩
    · rintro ⟨hx, hp⟩
      exact ⟨vertItem x, ⟨hp, hwf.tree.vert x hx⟩, vertItem_sub_one x⟩
  have hbet' : ∀ x ∈ xs, Precedes sq s x ∧ Precedes sq x t := by
    intro x hx
    obtain ⟨hxn, hp⟩ := (hxs_iff x).mp hx
    have := hbet s t hst (vertItem x) hp (hwf.tree.vert x hxn)
    rwa [vertItem_sub_one] at this
  -- order transfer on the V children
  have hxs_order : ∀ x ∈ xs, ∀ x' ∈ xs, Precedes sq x x' → Before xs x x' := by
    intro x hx x' hx' hlt
    obtain ⟨hxn, hpx⟩ := (hxs_iff x).mp hx
    obtain ⟨hxn', hpx'⟩ := (hxs_iff x').mp hx'
    have hne : x ≠ x' := hlt.ne
    have hbx : vertItem x ∈ b.items := hleaf _ (hchild_leaf _ hpx (Or.inl (hwf.tree.vert x hxn)))
    have hbx' : vertItem x' ∈ b.items :=
      hleaf _ (hchild_leaf _ hpx' (Or.inl (hwf.tree.vert x' hxn')))
    rcases Before.or_of_mem_mem hbx hbx' (fun h => hne (vertItem_inj h)) with hB | hB
    · have hB' : Before order (vertItem x) (vertItem x') := by
        rw [hord]; exact (hB.append_left L₁).append_right L₂
      have hch' := Items.before_ch_of_before_order hwf.tree hch hne hpx hpx'
        (hwf.tree.vert x hxn) (hwf.tree.vert x' hxn') hB'
      have := (hch'.filter_of (p := fun c => decide (items.type c = .V))
        (by simp [hwf.tree.vert x hxn]) (by simp [hwf.tree.vert x' hxn'])).map (· - 1)
      simpa [vertItem_sub_one, hxs] using this
    · exfalso
      have h1 : Before (b.items.filterMap (vertOf g)) x' x :=
        hB.filterMap (vertOf_vertItem hxn') (vertOf_vertItem hxn)
      have h2 : Before sq x' x := by rw [hsq, StBlock.seq_eq]; exact h1.append_left _
      exact hlt.asymm (Precedes.of_before hsqnd h2)
  -- order transfer from sq to the item's vertex list
  have hΦ : ∀ y ∈ V, ∀ y' ∈ V, Precedes sq y y' → Precedes V y y' := by
    intro y hy y' hy' hlt
    refine Precedes.of_before hVnd ?_
    have hne := hlt.ne
    rcases (hUmem y).mp hy with h | h | h <;> rcases (hUmem y').mp hy' with h' | h' | h'
    · exact absurd (h.trans h'.symm) hne
    · rw [h]; exact Before.cons_of_mem (List.mem_append_left _ h')
    · rw [h, h']; exact Before.cons_of_mem (List.mem_append_right _ (List.mem_singleton_self _))
    · rw [h'] at hlt; exact absurd (hbet' y h).1 hlt.asymm
    · exact ((hxs_order y h y' h' hlt).append_right [t]).cons s
    · rw [h']
      exact (Before.of_mem_mem (List.mem_cons_of_mem s h) (List.mem_singleton_self t) :
        Before ((s :: xs) ++ [t]) y t)
    · rw [h, h'] at hlt; exact absurd hst' hlt.asymm
    · rw [h] at hlt; exact absurd (hbet' y' h').2 hlt.asymm
    · exact absurd (h.trans h'.symm) hne
  -- the non-V children
  have hchild : ∀ c ∈ items.ch i, items.type c ≠ .V → ∃ u v, items.vs c = (some u, some v) ∧
      Precedes sq u v ∧ u ∈ V ∧ v ∈ V ∧ (u, v) ∈ items.virtualEdges i := by
    intro c hc hcV
    have hor_c := hcor c hc hcV
    obtain ⟨u, v, huv⟩ : ∃ u v, items.vs c = (some u, some v) := by
      rcases h : items.vs c with ⟨_ | u, _ | v⟩ <;> rw [h] at hor_c <;>
        first | exact ⟨u, v, rfl⟩ | (simp only [Oriented] at hor_c)
    have huv' : Precedes sq u v := by rw [huv] at hor_c; exact hor_c
    have hcvs := hwf.endpoints.child_vs_in_parent i c hc ht' hcV
    have hinV : ∀ w, (items.vs c).1 = some w ∨ (items.vs c).2 = some w → w ∈ V := by
      intro w hw
      have hwn : w < g.nv := hwf.endpoints.vs_lt c w hw
      rcases hcvs w hw with h | h | h
      · rw [hst] at h; simp at h; rw [← h]; simp [hVdef]
      · rw [hst] at h; simp at h; rw [← h]; simp [hVdef]
      · have : w ∈ xs := (hxs_iff w).mpr ⟨hwn, h⟩
        simp [hVdef, this]
    refine ⟨u, v, huv, huv', hinV u (Or.inl (by rw [huv])), hinV v (Or.inr (by rw [huv])), ?_⟩
    unfold Items.virtualEdges
    exact List.mem_map.mpr ⟨c, List.mem_filter.mpr ⟨hc, by simpa using hcV⟩, by simp [huv]⟩
  -- a block edge at a V child of `i` is below a non-V child whose endpoints bracket it
  have hnb : ∀ x ∈ xs, ∀ p ∈ b.edges g, (p.1 = x ∨ p.2 = x) →
      ∃ u v, (u, v) ∈ items.virtualEdges i ∧ u ∈ V ∧ v ∈ V ∧ Precedes sq u v ∧ (x = u ∨ x = v) ∧
        ∀ y, (p.1 = y ∨ p.2 = y) → (u = y ∨ Precedes sq u y) ∧ (y = v ∨ Precedes sq y v) := by
    intro x hx p hp hpx
    obtain ⟨hxn, hpx'⟩ := (hxs_iff x).mp hx
    obtain ⟨e, he, hpe, hbe⟩ : ∃ e, e < g.ne ∧ Items.PairEq g.edges[e]! p ∧
        (edgeItem g e ∈ b.items ∨ ∃ r ∈ b.root, Items.PairEq g.edges[e]! r) := by
      rw [StBlock.edges_eq] at hp
      rcases List.mem_append.mp hp with hp | hp
      · have hr : p ∈ b.root := by simpa using hp
        obtain ⟨e, he, hpe⟩ := hre p hr
        exact ⟨e, he, hpe, Or.inr ⟨p, hr, hpe⟩⟩
      · obtain ⟨y, hy, hye⟩ := List.mem_filterMap.mp hp
        obtain ⟨e, he, rfl, hep⟩ := edgeOf_eq_some hye
        exact ⟨e, he, Or.inl hep, Or.inl hy⟩
    have hinc : ∀ y, (p.1 = y ∨ p.2 = y) → (g.edges[e]!).1 = y ∨ (g.edges[e]!).2 = y := by
      intro y hy
      rcases hpe with hpe | hpe <;> rw [hpe] <;> simpa [or_comm] using hy
    have hbelow : items.Below i (edgeItem g e) :=
      ((hwf.endpoints.interior i x hi hxn).mp hpx').1 e he (hinc x hpx)
    have hlf := hedge e he hbe hbelow
    rw [hleaves] at hlf
    obtain ⟨c, hc, hcl⟩ := List.mem_flatMap.mp hlf
    have hcB : items.Below c (edgeItem g e) := Items.below_of_mem_leaves n c _ hcl
    have hcV : items.type c ≠ .V := by
      intro hV
      rw [Items.leaves_of_type n c (Or.inl hV)] at hcl
      simp at hcl
      rw [← hcl, hwf.tree.edge e he] at hV; cases hV
    obtain ⟨u, v, huv, huv', huV, hvV, hve⟩ := hchild c hc hcV
    have hsp := hspan c hc hcV u v huv e he hbe hcB
    have hxuv : x = u ∨ x = v := by
      have hnot := ((hwf.endpoints.interior i x hi hxn).mp hpx').2 c hc
      push Not at hnot
      obtain ⟨e', he', hinc', hnb'⟩ := hnot
      have hcT : items.type c ∉ [NodeType.F, .V] := by
        simp [hwf.tree.type_ne_F hc, hcV]
      have hsep := hwf.endpoints.separation c (hwf.tree.ch_lt i c hc) hcT x e e' hxn he he'
        (hinc x hpx) hinc' hcB hnb'
      rw [huv] at hsep
      simpa [eq_comm] using hsep
    exact ⟨u, v, hve, huV, hvV, huv', hxuv, fun y hy => hsp y (hinc y hy)⟩
  -- assemble
  refine ⟨s, t, hst, ?_, ?_⟩
  · rw [hV]
    refine ⟨hVnd, ?_, ?_⟩
    · intro p hp
      rcases List.mem_cons.mp hp with rfl | hp
      · exact ⟨by simp [hVdef], by simp [hVdef], hst'.ne⟩
      · unfold Items.virtualEdges at hp
        obtain ⟨c, hcf, rfl⟩ := List.mem_map.mp hp
        obtain ⟨hc, hcV⟩ := List.mem_filter.mp hcf
        simp at hcV
        obtain ⟨u, v, huv, hlt, huV, hvV, -⟩ := hchild c hc hcV
        rw [huv]
        exact ⟨huV, hvV, hlt.ne⟩
    · intro x hxV hhd hlast
      have hx : x ∈ xs := by
        rcases (hUmem x).mp hxV with h | h | h
        · exact absurd (by rw [h]; rfl) hhd
        · exact h
        · exact absurd (by rw [h]; exact List.getLast?_eq_some_iff.mpr ⟨s :: xs, rfl⟩) hlast
      obtain ⟨hsx, hxt⟩ := hbet' x hx
      have hxsq : x ∈ sq := hsx.2.1
      have hxhd : sq.head? ≠ some x := by
        intro h; have := idxOf_eq_zero_of_head? h; have := hsx.2.2; omega
      have hxls : sq.getLast? ≠ some x := by
        intro h; have := idxOf_le_of_getLast? hsqnd h hxt.2.1; have := hxt.2.2; omega
      obtain ⟨⟨p, hp, hplo⟩, ⟨q, hq, hqhi⟩⟩ := hsqint x hxsq hxhd hxls
      constructor
      · have hpx : p.1 = x ∨ p.2 = x := by
          rcases hplo with ⟨h, _⟩ | ⟨h, _⟩
          · exact Or.inl h
          · exact Or.inr h
        obtain ⟨u, v, hve, huV, hvV, huv, hxuv, hsp⟩ := hnb x hx p hp hpx
        obtain ⟨z, hzx, hz⟩ : ∃ z, (p.1 = z ∨ p.2 = z) ∧ Precedes sq z x := by
          rcases hplo with ⟨h1, h2⟩ | ⟨h1, h2⟩
          · exact ⟨p.2, Or.inr rfl, (hsqed p hp).2.1, hxsq, h2⟩
          · exact ⟨p.1, Or.inl rfl, (hsqed p hp).1, hxsq, h2⟩
        have hzsp := hsp z hzx
        have hxv : x = v := by
          rcases hxuv with h | h
          · exfalso
            rcases hzsp.1 with h' | h'
            · rw [h'.symm.trans h.symm] at hz; exact hz.ne rfl
            · rw [← h] at h'; exact hz.asymm h'
          · exact h
        rw [hxv] at hxV ⊢
        refine ⟨(u, v), List.mem_cons_of_mem _ hve, Or.inr ⟨rfl, ?_⟩⟩
        exact (hΦ u huV v hxV huv).2.2
      · have hqx : q.1 = x ∨ q.2 = x := by
          rcases hqhi with ⟨h, _⟩ | ⟨h, _⟩
          · exact Or.inl h
          · exact Or.inr h
        obtain ⟨u, v, hve, huV, hvV, huv, hxuv, hsp⟩ := hnb x hx q hq hqx
        obtain ⟨z, hzx, hz⟩ : ∃ z, (q.1 = z ∨ q.2 = z) ∧ Precedes sq x z := by
          rcases hqhi with ⟨h1, h2⟩ | ⟨h1, h2⟩
          · exact ⟨q.2, Or.inr rfl, hxsq, (hsqed q hq).2.1, h2⟩
          · exact ⟨q.1, Or.inl rfl, hxsq, (hsqed q hq).1, h2⟩
        have hzsp := hsp z hzx
        have hxu : x = u := by
          rcases hxuv with h | h
          · exact h
          · exfalso
            rcases hzsp.2 with h' | h'
            · rw [h'.trans h.symm] at hz; exact hz.ne rfl
            · rw [← h] at h'; exact hz.asymm h'
        rw [hxu] at hxV ⊢
        refine ⟨(u, v), List.mem_cons_of_mem _ hve, Or.inl ⟨rfl, ?_⟩⟩
        exact (hΦ u hxV v hvV huv).2.2
  · intro p hp
    rw [hV]
    unfold Items.virtualEdges at hp
    obtain ⟨c, hcf, rfl⟩ := List.mem_map.mp hp
    obtain ⟨hc, hcV⟩ := List.mem_filter.mp hcf
    simp at hcV
    obtain ⟨u, v, huv, hlt, huV, hvV, -⟩ := hchild c hc hcV
    rw [huv]
    exact (hΦ u huV v hvV hlt).2.2

/-- An S/P/R item of the walk whose children are in reference order (`walk_st'`) and which is
oriented along its block (`walk_vsOriented`) is `Items.StItem`: `stItem_of_block` on the block
supplied by `VsOriented`, using `refBlocks_st` and `refBlocks_root_edge`. -/
theorem stItem_of_refOrder (g : Graph) (tern : Bool) (vo eo : List Nat) (i : ItemId)
    (hi : i < (g.walk tern (g.dfsForest vo eo)).items.size)
    (ht : Items.type (g.walk tern (g.dfsForest vo eo)).items i = .S ∨
      Items.type (g.walk tern (g.dfsForest vo eo)).items i = .P ∨
      Items.type (g.walk tern (g.dfsForest vo eo)).items i = .R)
    (hwf : Items.WF g (g.walk tern (g.dfsForest vo eo)).items)
    (hch : Items.ch (g.walk tern (g.dfsForest vo eo)).items i =
      restrictCh (g.walk tern (g.dfsForest vo eo)).items
        (g.walk tern (g.dfsForest vo eo)).items.size (refOrder g (g.dfsForest vo eo)) i)
    (hvs : VsOriented g (g.walk tern (g.dfsForest vo eo)).items (refBlocks g (g.dfsForest vo eo)))
    (hbl : ∀ b ∈ refBlocks g (g.dfsForest vo eo), b.St g)
    (hre : ∀ b ∈ refBlocks g (g.dfsForest vo eo), ∀ r ∈ b.root,
      ∃ e, e < g.ne ∧ Items.PairEq g.edges[e]! r) :
    Items.StItem (g.walk tern (g.dfsForest vo eo)).items i := by
  obtain ⟨b, hb, hleaf, hedge, hor, hbet, hcor, hspan⟩ := hvs i hi ht
  obtain ⟨L₁, L₂, hord⟩ : ∃ L₁ L₂, refOrder g (g.dfsForest vo eo) = L₁ ++ b.items ++ L₂ := by
    obtain ⟨B₁, B₂, hB⟩ := List.append_of_mem hb
    refine ⟨B₁.flatMap StBlock.items, B₂.flatMap StBlock.items, ?_⟩
    unfold refOrder; rw [hB]; simp [List.flatMap_append, List.flatMap_cons, List.append_assoc]
  exact stItem_of_block g _ hwf i hi ht b _ L₁ L₂ hord hch (hbl b hb) (hre b hb) hleaf hedge hor
    hbet hcor hspan

end Spqr
