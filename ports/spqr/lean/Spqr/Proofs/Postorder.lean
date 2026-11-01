import Spqr.SepPair

/-!
# The postorder edge order: blocks are intervals and form a laminar family (Fact D, structure)
-/

namespace Spqr

theorem List.infix_flatMap_of_mem {α β : Type} {f : α → List β} {a : α} {l : List α} (h : a ∈ l) :
    f a <:+: l.flatMap f := by
  induction l with
  | nil => exact absurd h (List.not_mem_nil)
  | cons b rest ih =>
    rw [List.flatMap_cons]
    rcases List.mem_cons.mp h with rfl | h
    · exact (List.prefix_append _ _).isInfix
    · exact (ih h).trans (List.suffix_append _ _).isInfix

theorem List.eq_or_disjoint_of_nodup_flatMap {α β : Type} {f : α → List β} {l : List α}
    (hnd : (l.flatMap f).Nodup) {a b : α} (ha : a ∈ l) (hb : b ∈ l) :
    a = b ∨ ∀ x, x ∈ f a → x ∈ f b → False := by
  induction l with
  | nil => exact absurd ha (List.not_mem_nil)
  | cons z rest ih =>
    rw [List.flatMap_cons, List.nodup_append] at hnd
    obtain ⟨-, hnd, hdisj⟩ := hnd
    rcases List.mem_cons.mp ha with rfl | ha' <;> rcases List.mem_cons.mp hb with rfl | hb'
    · exact .inl rfl
    · exact .inr fun x hx hx' => hdisj _ hx _ ((List.infix_flatMap_of_mem hb').subset hx') rfl
    · exact .inr fun x hx hx' => hdisj _ hx' _ ((List.infix_flatMap_of_mem ha').subset hx) rfl
    · exact ih hnd ha' hb'

theorem List.eq_cons_of_mem_of_not_sublist {α : Type} {l : List α} {x : α} (hx : x ∈ l)
    (h : ∀ y, [y, x].Sublist l → False) : ∃ rest, l = x :: rest := by
  cases l with
  | nil => exact absurd hx (List.not_mem_nil)
  | cons y rest =>
    rcases List.mem_cons.mp hx with rfl | hx'
    · exact ⟨rest, rfl⟩
    · exact absurd (List.Sublist.cons_cons y (List.singleton_sublist.mpr hx')) (h y)

theorem List.mem_dropWhile_sorted {α : Type} {f : α → Nat} {k : Nat} {l : List α}
    (hs : l.Pairwise fun x y => f x ≤ f y) {x : α} :
    x ∈ l.dropWhile (fun y => decide (f y ≤ k)) ↔ x ∈ l ∧ k < f x := by
  induction l with
  | nil => simp
  | cons y rest ih =>
    rw [List.pairwise_cons] at hs
    rw [List.dropWhile_cons]
    by_cases hy : f y ≤ k
    · rw [ite_eq_left_of_eq_true _ _ (eq_true (decide_eq_true hy)), ih hs.2, List.mem_cons]
      constructor
      · rintro ⟨h1, h2⟩; exact ⟨.inr h1, h2⟩
      · rintro ⟨rfl | h1, h2⟩
        · omega
        · exact ⟨h1, h2⟩
    · rw [ite_eq_right_of_eq_false _ _ (eq_false (by simpa using hy))]
      constructor
      · intro h
        refine ⟨h, ?_⟩
        rcases List.mem_cons.mp h with rfl | h
        · omega
        · have := hs.1 x h; omega
      · exact fun h => h.1

theorem DfsOut.edgePostorderList_eq_flatMap (outs : List DfsOut) :
    DfsOut.edgePostorderList outs = outs.flatMap DfsOut.block := by
  induction outs with
  | nil => rfl
  | cons o rest ih => cases o <;> simp [DfsOut.edgePostorderList, DfsOut.block, ih]

theorem DfsOut.block_infix_of_mem {o : DfsOut} {outs : List DfsOut} (h : o ∈ outs) :
    o.block <:+: DfsOut.edgePostorderList outs := by
  rw [DfsOut.edgePostorderList_eq_flatMap]
  exact List.infix_flatMap_of_mem h

theorem DfsOut.edgePostorderList_tree_cons (e : Nat) (cls : OutClass) (child : DfsTree)
    (rest : List DfsOut) :
    DfsOut.edgePostorderList (.tree e cls child :: rest) =
      (child.edgePostorder ++ [e]) ++ DfsOut.edgePostorderList rest := by
  simp [DfsOut.edgePostorderList]

mutual
/-- Every block is an interval of the postorder. -/
theorem DfsTree.blocks_infix : ∀ (t : DfsTree), ∀ B ∈ t.blocks, B <:+: t.edgePostorder
  | .node _ outs => DfsOut.blocksList_infix outs
theorem DfsOut.blocksList_infix :
    ∀ (outs : List DfsOut), ∀ B ∈ DfsOut.blocksList outs, B <:+: DfsOut.edgePostorderList outs
  | [], _, h => absurd h (List.not_mem_nil)
  | .back _ _ _ :: rest, B, h =>
    (DfsOut.blocksList_infix rest B h).trans (List.suffix_append [_] _).isInfix
  | .tree e cls child :: rest, B, h => by
    rw [DfsOut.edgePostorderList_tree_cons]
    simp only [DfsOut.blocksList, List.mem_cons, List.mem_append] at h
    rcases h with rfl | h | h
    · exact (List.prefix_append _ _).isInfix
    · exact ((DfsTree.blocks_infix child B h).trans (List.prefix_append _ _).isInfix).trans
        (List.prefix_append _ _).isInfix
    · exact (DfsOut.blocksList_infix rest B h).trans (List.suffix_append _ _).isInfix
end

mutual
/-- The blocks below a tree with distinct edges are pairwise nested or disjoint. -/
theorem DfsTree.blocks_laminar :
    ∀ (t : DfsTree), t.edgePostorder.Nodup → ∀ B₁ ∈ t.blocks, ∀ B₂ ∈ t.blocks, Laminar B₁ B₂
  | .node _ outs => DfsOut.blocksList_laminar outs
theorem DfsOut.blocksList_laminar :
    ∀ (outs : List DfsOut), (DfsOut.edgePostorderList outs).Nodup →
      ∀ B₁ ∈ DfsOut.blocksList outs, ∀ B₂ ∈ DfsOut.blocksList outs, Laminar B₁ B₂
  | [], _, _, h, _, _ => absurd h (List.not_mem_nil)
  | .back _ _ _ :: rest, hnd, B₁, h₁, B₂, h₂ =>
    DfsOut.blocksList_laminar rest (List.nodup_cons.mp hnd).2 B₁ h₁ B₂ h₂
  | .tree e cls child :: rest, hnd, B₁, h₁, B₂, h₂ => by
    rw [DfsOut.edgePostorderList_tree_cons, List.nodup_append] at hnd
    obtain ⟨hndb, hndr, hdisj⟩ := hnd
    have hndc : child.edgePostorder.Nodup := (List.nodup_append.mp hndb).1
    have hsubB : ∀ B, B = child.edgePostorder ++ [e] ∨ B ∈ child.blocks →
        B ⊆ child.edgePostorder ++ [e] := fun B h =>
      h.elim (fun h => h ▸ List.Subset.refl _) fun h =>
        (DfsTree.blocks_infix child B h).subset.trans (List.prefix_append _ _).subset
    have hsubR : ∀ B, B ∈ DfsOut.blocksList rest → B ⊆ DfsOut.edgePostorderList rest :=
      fun B h => (DfsOut.blocksList_infix rest B h).subset
    simp only [DfsOut.blocksList, List.mem_cons, List.mem_append] at h₁ h₂
    rcases h₁ with h₁ | h₁ | h₁ <;> rcases h₂ with h₂ | h₂ | h₂
    · exact .inl (h₁ ▸ h₂ ▸ List.Subset.refl _)
    · exact .inr (.inl (h₁ ▸ hsubB B₂ (.inr h₂)))
    · exact .inr (.inr fun x hx hx' => hdisj _ (hsubB B₁ (.inl h₁) hx) _ (hsubR B₂ h₂ hx') rfl)
    · exact .inl (h₂ ▸ hsubB B₁ (.inr h₁))
    · exact DfsTree.blocks_laminar child hndc B₁ h₁ B₂ h₂
    · exact .inr (.inr fun x hx hx' => hdisj _ (hsubB B₁ (.inr h₁) hx) _ (hsubR B₂ h₂ hx') rfl)
    · exact .inr (.inr fun x hx hx' => hdisj _ (hsubB B₂ (.inl h₂) hx') _ (hsubR B₁ h₁ hx) rfl)
    · exact .inr (.inr fun x hx hx' => hdisj _ (hsubB B₂ (.inr h₂) hx') _ (hsubR B₁ h₁ hx) rfl)
    · exact DfsOut.blocksList_laminar rest hndr B₁ h₁ B₂ h₂
end

theorem DfsOut.blocks_subset_blocksList {outs : List DfsOut} {e : Nat} {cls : OutClass}
    {child : DfsTree} (hmem : DfsOut.tree e cls child ∈ outs) :
    child.blocks ⊆ DfsOut.blocksList outs := by
  induction outs with
  | nil => exact absurd hmem (List.not_mem_nil)
  | cons o rest ih =>
    rcases List.mem_cons.mp hmem with rfl | hmem'
    · exact fun B hB => by simp [DfsOut.blocksList, hB]
    · exact (ih hmem').trans fun B hB => by cases o <;> simp [DfsOut.blocksList, hB]

theorem DfsOut.block_mem_blocksList {outs : List DfsOut} {e : Nat} {cls : OutClass}
    {child : DfsTree} (hmem : DfsOut.tree e cls child ∈ outs) :
    (DfsOut.tree e cls child).block ∈ DfsOut.blocksList outs := by
  induction outs with
  | nil => exact absurd hmem (List.not_mem_nil)
  | cons o rest ih =>
    rcases List.mem_cons.mp hmem with rfl | hmem'
    · simp [DfsOut.blocksList, DfsOut.block]
    · cases o <;> simp [DfsOut.blocksList, ih hmem']

namespace DfsTree

theorem Sub.trans {s t u : DfsTree} (h : s.Sub t) (h' : t.Sub u) : s.Sub u := by
  induction h' with
  | refl => exact h
  | step _ hmem ih => exact .step ih hmem

theorem Sub.child {s : DfsTree} {e : Nat} {cls : OutClass} {child : DfsTree}
    (ho : DfsOut.tree e cls child ∈ s.outs) : child.Sub s := by
  cases s; exact .step (.refl _) ho

theorem Sub.edgePostorder_infix {s t : DfsTree} (h : s.Sub t) :
    s.edgePostorder <:+: t.edgePostorder := by
  induction h with
  | refl => exact List.infix_rfl
  | step _ hmem ih =>
    exact (ih.trans (List.prefix_append _ _).isInfix).trans (DfsOut.block_infix_of_mem hmem)

theorem Sub.blocks_subset {s t : DfsTree} (h : s.Sub t) : s.blocks ⊆ t.blocks := by
  induction h with
  | refl => exact List.Subset.refl _
  | step _ hmem ih => exact ih.trans (DfsOut.blocks_subset_blocksList hmem)

/-- The block of a tree out-edge of a subtree is one of the tree's blocks. -/
theorem block_mem_blocks {s t : DfsTree} (h : s.Sub t) {e : Nat} {cls : OutClass} {child : DfsTree}
    (hmem : DfsOut.tree e cls child ∈ s.outs) : (DfsOut.tree e cls child).block ∈ t.blocks :=
  h.blocks_subset (by cases s; exact DfsOut.block_mem_blocksList hmem)

end DfsTree

/-- Fact D, structure: the block of every tree edge is an interval of the walk's edge order. -/
theorem block_infix_forest {forest : List DfsTree} {t s : DfsTree} (ht : t ∈ forest) (h : s.Sub t)
    {e : Nat} {cls : OutClass} {child : DfsTree} (hmem : DfsOut.tree e cls child ∈ s.outs) :
    (DfsOut.tree e cls child).block <:+: edgePostorderForest forest :=
  (t.blocks_infix _ (DfsTree.block_mem_blocks h hmem)).trans (List.infix_flatMap_of_mem ht)

/-- Fact D, structure: blocks of tree edges are pairwise nested or disjoint. -/
theorem blocks_laminar_forest {forest : List DfsTree} (hnd : (edgePostorderForest forest).Nodup)
    {t₁ t₂ : DfsTree} (h₁ : t₁ ∈ forest) (h₂ : t₂ ∈ forest) {B₁ B₂ : List Nat}
    (hB₁ : B₁ ∈ t₁.blocks) (hB₂ : B₂ ∈ t₂.blocks) : Laminar B₁ B₂ := by
  rcases List.eq_or_disjoint_of_nodup_flatMap hnd h₁ h₂ with rfl | hdisj
  · exact t₁.blocks_laminar ((hnd.sublist (List.infix_flatMap_of_mem h₁).sublist)) B₁ hB₁ B₂ hB₂
  · exact .inr (.inr fun x hx hx' =>
      hdisj _ ((t₁.blocks_infix _ hB₁).subset hx) ((t₂.blocks_infix _ hB₂).subset hx'))

end Spqr

/-! ### Vertices, out-lists and edges of subtrees -/

namespace Spqr

theorem DfsTree.v_mem_verts (t : DfsTree) : t.v ∈ t.verts := by
  cases t; simp [DfsTree.verts, DfsTree.v]

theorem DfsOut.mem_vertsList {outs : List DfsOut} {x : Nat} :
    x ∈ DfsOut.vertsList outs ↔
      ∃ e cls child, DfsOut.tree e cls child ∈ outs ∧ x ∈ child.verts := by
  induction outs with
  | nil => simp [DfsOut.vertsList]
  | cons o rest ih =>
    cases o with
    | back e dest cls =>
      simp only [DfsOut.vertsList, ih, List.mem_cons]
      constructor
      · rintro ⟨e', c', ch, h, hx⟩
        exact ⟨e', c', ch, .inr h, hx⟩
      · rintro ⟨e', c', ch, h | h, hx⟩
        · cases h
        · exact ⟨e', c', ch, h, hx⟩
    | tree e cls child =>
      simp only [DfsOut.vertsList, ih, List.mem_cons, List.mem_append]
      constructor
      · rintro (hx | ⟨e', c', ch, h, hx⟩)
        · exact ⟨e, cls, child, .inl rfl, hx⟩
        · exact ⟨e', c', ch, .inr h, hx⟩
      · rintro ⟨e', c', ch, h | h, hx⟩
        · cases h
          exact .inl hx
        · exact .inr ⟨e', c', ch, h, hx⟩

theorem DfsOut.verts_infix_of_mem {outs : List DfsOut} {e : Nat} {cls : OutClass} {child : DfsTree}
    (hmem : DfsOut.tree e cls child ∈ outs) : child.verts <:+: DfsOut.vertsList outs := by
  induction outs with
  | nil => exact absurd hmem (List.not_mem_nil)
  | cons o rest ih =>
    rcases List.mem_cons.mp hmem with rfl | hmem'
    · exact (List.prefix_append _ _).isInfix
    · cases o with
      | back => exact ih hmem'
      | tree => exact (ih hmem').trans (List.suffix_append _ _).isInfix

theorem DfsTree.Sub.verts_subset {s t : DfsTree} (h : s.Sub t) : s.verts ⊆ t.verts := by
  induction h with
  | refl => exact List.Subset.refl _
  | step _ hmem ih =>
    exact ih.trans fun x hx =>
      List.mem_cons_of_mem _ (DfsOut.mem_vertsList.mpr ⟨_, _, _, hmem, hx⟩)

mutual
theorem DfsTree.exists_sub_of_mem_verts :
    ∀ (t : DfsTree) {x : Nat}, x ∈ t.verts → ∃ s : DfsTree, s.Sub t ∧ s.v = x
  | .node v outs, x, h => by
    rcases List.mem_cons.mp h with rfl | h
    · exact ⟨_, .refl _, rfl⟩
    · obtain ⟨e, cls, child, hmem, s, hs, rfl⟩ := DfsOut.exists_sub_of_mem_vertsList outs h
      exact ⟨s, .step hs hmem, rfl⟩
theorem DfsOut.exists_sub_of_mem_vertsList :
    ∀ (outs : List DfsOut) {x : Nat}, x ∈ DfsOut.vertsList outs →
      ∃ e cls child, DfsOut.tree e cls child ∈ outs ∧ ∃ s : DfsTree, s.Sub child ∧ s.v = x
  | [], _, h => absurd h (List.not_mem_nil)
  | .back _ _ _ :: rest, x, h => by
    obtain ⟨e, cls, child, hmem, hs⟩ := DfsOut.exists_sub_of_mem_vertsList rest h
    exact ⟨e, cls, child, List.mem_cons_of_mem _ hmem, hs⟩
  | .tree e cls child :: rest, x, h => by
    rcases List.mem_append.mp h with h | h
    · exact ⟨e, cls, child, List.mem_cons_self, DfsTree.exists_sub_of_mem_verts child h⟩
    · obtain ⟨e', cls', child', hmem, hs⟩ := DfsOut.exists_sub_of_mem_vertsList rest h
      exact ⟨e', cls', child', List.mem_cons_of_mem _ hmem, hs⟩
end

mutual
theorem DfsTree.outsAt_eq_nil : ∀ (t : DfsTree) {x : Nat}, x ∉ t.verts → t.outsAt x = []
  | .node v outs, x, h => by
    have hv : v ≠ x := fun hvx => h (hvx ▸ List.mem_cons_self)
    simp only [DfsTree.outsAt, hv, ↓reduceIte, List.nil_append]
    exact DfsOut.outsAtList_eq_nil outs fun h' => h (List.mem_cons_of_mem _ h')
theorem DfsOut.outsAtList_eq_nil :
    ∀ (outs : List DfsOut) {x : Nat}, x ∉ DfsOut.vertsList outs → DfsOut.outsAtList x outs = []
  | [], _, _ => rfl
  | .back _ _ _ :: rest, x, h => DfsOut.outsAtList_eq_nil rest h
  | .tree e cls child :: rest, x, h => by
    simp only [DfsOut.vertsList, List.mem_append, not_or] at h
    simp only [DfsOut.outsAtList, DfsTree.outsAt_eq_nil child h.1,
      DfsOut.outsAtList_eq_nil rest h.2, List.nil_append]
end

theorem DfsOut.outsAtList_eq_of_mem {outs : List DfsOut} (hnd : (DfsOut.vertsList outs).Nodup)
    {e : Nat} {cls : OutClass} {child : DfsTree} (hmem : DfsOut.tree e cls child ∈ outs) {x : Nat}
    (hx : x ∈ child.verts) : DfsOut.outsAtList x outs = child.outsAt x := by
  induction outs with
  | nil => exact absurd hmem (List.not_mem_nil)
  | cons o rest ih =>
    rcases List.mem_cons.mp hmem with rfl | hmem'
    · simp only [DfsOut.vertsList, List.nodup_append] at hnd
      simp only [DfsOut.outsAtList,
        DfsOut.outsAtList_eq_nil rest fun h => hnd.2.2 _ hx _ h rfl, List.append_nil]
    · cases o with
      | back => exact ih hnd hmem'
      | tree e' cls' child' =>
        simp only [DfsOut.vertsList, List.nodup_append] at hnd
        have hx' : x ∈ DfsOut.vertsList rest := DfsOut.mem_vertsList.mpr ⟨_, _, _, hmem', hx⟩
        simp only [DfsOut.outsAtList,
          DfsTree.outsAt_eq_nil child' fun h => hnd.2.2 _ h _ hx' rfl, List.nil_append]
        exact ih hnd.2.1 hmem'

theorem DfsTree.Sub.outsAt_eq {s t : DfsTree} (hnd : t.verts.Nodup) (h : s.Sub t) :
    t.outsAt s.v = s.outs := by
  induction h with
  | refl =>
    obtain ⟨v, outs⟩ := s
    simp only [DfsTree.verts, List.nodup_cons] at hnd
    simp [DfsTree.outsAt, DfsTree.v, DfsTree.outs, DfsOut.outsAtList_eq_nil outs hnd.1]
  | @step v outs e cls child hsub hmem ih =>
    simp only [DfsTree.verts, List.nodup_cons] at hnd
    have hsv : s.v ∈ child.verts := hsub.verts_subset s.v_mem_verts
    have hne : v ≠ s.v := fun h => hnd.1 (h ▸ DfsOut.mem_vertsList.mpr ⟨_, _, _, hmem, hsv⟩)
    simp only [DfsTree.outsAt, hne, ↓reduceIte, List.nil_append]
    rw [DfsOut.outsAtList_eq_of_mem hnd.2 hmem hsv]
    exact ih (hnd.2.sublist (DfsOut.verts_infix_of_mem hmem).sublist)

theorem outsAt_forest_eq {forest : List DfsTree} (hnd : (forest.flatMap DfsTree.verts).Nodup)
    {t : DfsTree} (ht : t ∈ forest) {x : Nat} (hx : x ∈ t.verts) :
    forest.flatMap (·.outsAt x) = t.outsAt x := by
  induction forest with
  | nil => exact absurd ht (List.not_mem_nil)
  | cons t₀ rest ih =>
    rw [List.flatMap_cons, List.nodup_append] at hnd
    rw [List.flatMap_cons]
    rcases List.mem_cons.mp ht with rfl | ht'
    · rw [List.flatMap_eq_nil_iff.mpr fun t' ht' => t'.outsAt_eq_nil fun h =>
        hnd.2.2 _ hx _ ((List.infix_flatMap_of_mem ht').subset h) rfl, List.append_nil]
    · rw [t₀.outsAt_eq_nil fun h => hnd.2.2 _ h _ ((List.infix_flatMap_of_mem ht').subset hx) rfl,
        List.nil_append]
      exact ih hnd.2.1 ht'

/-- The abstract out-list of a vertex of the forest is its node's out-list. -/
theorem DfsData.ofForest_outs {forest : List DfsTree} (hnd : (forest.flatMap DfsTree.verts).Nodup)
    {t s : DfsTree} (ht : t ∈ forest) (hsub : s.Sub t) :
    (DfsData.ofForest forest).outs s.v = s.outs := by
  show forest.flatMap (·.outsAt s.v) = s.outs
  rw [outsAt_forest_eq hnd ht (hsub.verts_subset s.v_mem_verts)]
  exact hsub.outsAt_eq (hnd.sublist (List.infix_flatMap_of_mem ht).sublist)

theorem DfsTree.edgePostorder_eq (t : DfsTree) :
    t.edgePostorder = DfsOut.edgePostorderList t.outs := by
  cases t; rfl

theorem DfsOut.edgePostorderList_takeWhile_dropWhile (p : DfsOut → Bool) (l : List DfsOut) :
    DfsOut.edgePostorderList l =
      DfsOut.edgePostorderList (l.takeWhile p) ++ DfsOut.edgePostorderList (l.dropWhile p) := by
  rw [DfsOut.edgePostorderList_eq_flatMap, DfsOut.edgePostorderList_eq_flatMap,
    DfsOut.edgePostorderList_eq_flatMap, ← List.flatMap_append, List.takeWhile_append_dropWhile]

theorem DfsOut.edgePostorderList_subset_of_sublist {l l' : List DfsOut} (h : l.Sublist l') :
    DfsOut.edgePostorderList l ⊆ DfsOut.edgePostorderList l' := by
  rw [DfsOut.edgePostorderList_eq_flatMap, DfsOut.edgePostorderList_eq_flatMap]
  intro x hx
  obtain ⟨o, ho, hx⟩ := List.mem_flatMap.mp hx
  exact List.mem_flatMap.mpr ⟨o, h.subset ho, hx⟩

theorem DfsTree.mem_type2Block {c sb : DfsTree} (hnd : c.edgePostorder.Nodup)
    (hpre : sb.edgePostorder <+: c.edgePostorder) {l e₀ x : Nat} :
    x ∈ DfsTree.type2Block l e₀ c sb ↔
      x = e₀ ∨ (x ∈ c.edgePostorder ∧ x ∉ sb.edgePostorder) ∨
        x ∈ DfsOut.edgePostorderList
          (sb.outs.dropWhile fun o => decide (o.cls.rank ≤ (OutClass.ret l .backEdge).rank)) := by
  obtain ⟨R, hR⟩ := hpre
  unfold DfsTree.type2Block
  rw [← hR, List.drop_left]
  rw [← hR, List.nodup_append] at hnd
  simp only [List.mem_append, List.mem_singleton]
  constructor
  · rintro ((h | h) | h)
    · exact .inr (.inr h)
    · exact .inr (.inl ⟨.inr h, fun h' => hnd.2.2 _ h' _ h rfl⟩)
    · exact .inl h
  · rintro (h | ⟨h | h, h'⟩ | h)
    · exact .inr h
    · exact absurd h h'
    · exact .inl (.inr h)
    · exact .inl (.inl h)

theorem DfsTree.type2Block_disjoint_takeWhile {c sb : DfsTree} {e₀ : Nat}
    (hnd : (c.edgePostorder ++ [e₀]).Nodup) (hpre : sb.edgePostorder <+: c.edgePostorder)
    {l x : Nat}
    (hx : x ∈ DfsOut.edgePostorderList
      (sb.outs.takeWhile fun o => decide (o.cls.rank ≤ (OutClass.ret l .backEdge).rank))) :
    x ∉ DfsTree.type2Block l e₀ c sb := by
  obtain ⟨R, hR⟩ := hpre
  unfold DfsTree.type2Block
  rw [← hR, List.drop_left]
  rw [← hR, sb.edgePostorder_eq, DfsOut.edgePostorderList_takeWhile_dropWhile
    (fun o => decide (o.cls.rank ≤ (OutClass.ret l .backEdge).rank)), List.append_assoc,
    List.append_assoc, List.nodup_append] at hnd
  rw [List.append_assoc]
  exact fun h => hnd.2.2 _ hx _ h rfl

theorem DfsOut.e_mem_block (o : DfsOut) : o.e ∈ o.block := by
  cases o <;> simp [DfsOut.block, DfsOut.e]

theorem DfsTree.Sub.e_mem_edgePostorder {s t : DfsTree} (h : s.Sub t) {o : DfsOut}
    (ho : o ∈ s.outs) : o.e ∈ t.edgePostorder :=
  h.edgePostorder_infix.subset
    ((by cases s; exact DfsOut.block_infix_of_mem ho : o.block <:+: s.edgePostorder).subset
      o.e_mem_block)

mutual
theorem DfsTree.mem_edgePostorder :
    ∀ (t : DfsTree) {e : Nat}, e ∈ t.edgePostorder → ∃ s : DfsTree, s.Sub t ∧ ∃ o ∈ s.outs, o.e = e
  | .node _ outs, e, h => by
    rcases DfsOut.mem_edgePostorderList outs h with ⟨o, ho, he⟩ | ⟨_, _, _, hmem, s, hs, ho⟩
    · exact ⟨_, .refl _, o, ho, he⟩
    · exact ⟨s, .step hs hmem, ho⟩
theorem DfsOut.mem_edgePostorderList :
    ∀ (outs : List DfsOut) {e : Nat}, e ∈ DfsOut.edgePostorderList outs →
      (∃ o ∈ outs, o.e = e) ∨
        ∃ e₁ cls child, DfsOut.tree e₁ cls child ∈ outs ∧ ∃ s : DfsTree, s.Sub child ∧ ∃ o ∈ s.outs, o.e = e
  | [], _, h => absurd h (List.not_mem_nil)
  | .back _ _ _ :: rest, e, h => by
    rcases List.mem_cons.mp h with rfl | h
    · exact .inl ⟨_, List.mem_cons_self, rfl⟩
    · rcases DfsOut.mem_edgePostorderList rest h with ⟨o, ho, he⟩ | ⟨e₁, cls, child, hmem, hs⟩
      · exact .inl ⟨o, List.mem_cons_of_mem _ ho, he⟩
      · exact .inr ⟨e₁, cls, child, List.mem_cons_of_mem _ hmem, hs⟩
  | .tree e₁ cls child :: rest, e, h => by
    simp only [DfsOut.edgePostorderList, List.mem_append, List.mem_cons] at h
    rcases h with h | rfl | h
    · obtain ⟨s, hs, ho⟩ := DfsTree.mem_edgePostorder child h
      exact .inr ⟨e₁, cls, child, List.mem_cons_self, s, hs, ho⟩
    · exact .inl ⟨_, List.mem_cons_self, rfl⟩
    · rcases DfsOut.mem_edgePostorderList rest h with ⟨o, ho, he⟩ | ⟨e₂, cls', child', hmem, hs⟩
      · exact .inl ⟨o, List.mem_cons_of_mem _ ho, he⟩
      · exact .inr ⟨e₂, cls', child', List.mem_cons_of_mem _ hmem, hs⟩
end

end Spqr
