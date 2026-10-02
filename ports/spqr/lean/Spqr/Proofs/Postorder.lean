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
