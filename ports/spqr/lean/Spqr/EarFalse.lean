import Spqr.WalkInv
import Spqr.EarRoot

/-!
# `walkTree_ear`/`walkTree_book` are false without edge completeness

The hypotheses of `WalkState.walkTree_ear` (and of `walkTree_book`, `walkTree_guards`,
`walk_sides`) — `t.WF []`, `t.Ends g`, bounds, nodup, the start state — allow graph edges that
are incident to a vertex of `t` but are not edges of `t` (`WF` only classifies the tree's own
edges). `VertBook`'s `TwoAttached` at the end of a vertex's outs is then false: a path
`0 - 1 - 2` (two bridges) with an extra edge `2 - 3` outside the tree. After the bridge `1 - 2`
is finished, the edge set below `V 1` is `{1 - 2}`, attached at vertex `2` by the edge `2 - 3`,
which is neither of the two terminals (`1`, `1`). Kernel-checked below (`decide +kernel` runs the
walk on the concrete state). The statements below are the pre-correction ones; the correction is
the per-tree edge completeness hypothesis `hcomp` (every edge of `g` incident to a vertex of `t` is
an edge of `t`) on `walkTree_ear`/`walkTree_book`/`walkTree_guards`/`RootState.book`, from the
forest's edge coverage via `comp_of_forest` (`EarWalk.lean`); see PROOF.md §4.2b.
-/

namespace Spqr
namespace WalkState
open WalkM


instance (items : Items) (p c : ItemId) : Decidable (items.IsParent p c) := by
  unfold Items.IsParent; infer_instance

theorem noParent_of_ball {items : Items} {c : ItemId}
    (h : ∀ p, p < items.size → ¬ items.IsParent p c) : ∀ p, ¬ items.IsParent p c := fun p hp =>
  if hlt : p < items.size then h p hlt hp else by
    rw [Items.IsParent, ch_of_ge (Nat.le_of_not_lt hlt)] at hp; exact List.not_mem_nil hp

theorem walkTree_ear_false :
    ∃ (t : DfsTree) (s : WalkState), t.WF [] ∧ t.Ends s.g ∧
      (∀ v ∈ t.verts, v < s.g.nv) ∧ (∀ e ∈ t.edges, e < s.g.ne) ∧
      t.verts.Nodup ∧ t.edges.Nodup ∧
      s.stackVerts.size = s.g.nv ∧ s.stackDir.size = s.g.nv ∧ s.firstOccurrence.size = s.g.nv ∧
      s.tstack = [] ∧ s.Inv' 0 ∧ Shape s ∧
      (∀ v ∈ t.verts, Items.ch s.items (vertItem v) = [] ∧
        ∀ p, ¬ Items.IsParent s.items p (vertItem v)) ∧
      (∀ e ∈ t.edges, Items.ch s.items (edgeItem s.g e) = [] ∧
        ∀ p, ¬ Items.IsParent s.items p (edgeItem s.g e)) ∧
      ¬ EarTree t 0 s := by
  let g : Graph := ⟨4, #[(1, 0), (2, 1), (2, 3)]⟩
  let t : DfsTree := .node 0 [.tree 0 .bridge (.node 1 [.tree 1 .bridge (.node 2 [])])]
  refine ⟨t, WalkState.init g false, ?_, ?_, ?_, ?_, ?_, ?_, rfl, rfl, rfl, rfl, init_inv g false,
    init_shape g false, ?_, ?_, ?_⟩
  · unfold DfsTree.WF
    refine ⟨List.pairwise_singleton .., fun o ho => ?_⟩
    simp only [List.mem_singleton] at ho; subst ho
    unfold DfsOut.WF DfsTree.WF
    refine ⟨⟨List.pairwise_singleton .., fun o ho => ?_⟩, rfl⟩
    simp only [List.mem_singleton] at ho; subst ho
    unfold DfsOut.WF DfsTree.WF
    exact ⟨⟨List.Pairwise.nil, fun o ho => nomatch ho⟩, rfl⟩
  · unfold DfsTree.Ends
    intro o ho
    simp only [List.mem_singleton] at ho; subst ho
    unfold DfsOut.Ends DfsTree.Ends
    refine ⟨Or.inl rfl, fun o ho => ?_⟩
    simp only [List.mem_singleton] at ho; subst ho
    unfold DfsOut.Ends DfsTree.Ends
    exact ⟨Or.inl rfl, fun o ho => nomatch ho⟩
  · decide
  · decide
  · decide
  · decide
  · intro v _
    refine ⟨Items.initialItems_ch g _, fun p h => ?_⟩
    rw [Items.IsParent, show (init g false).items = initialItems g from rfl,
      Items.initialItems_ch] at h
    exact absurd h (List.not_mem_nil)
  · intro e _
    refine ⟨Items.initialItems_ch g _, fun p h => ?_⟩
    rw [Items.IsParent, show (init g false).items = initialItems g from rfl,
      Items.initialItems_ch] at h
    exact absurd h (List.not_mem_nil)
  · intro h
    unfold EarTree EarOuts EarOut at h
    have h2 := h.1.2
    simp only [WalkM.wp] at h2
    unfold EarTree EarOuts at h2
    have h3 := h2.1.2
    simp only [WalkM.wp, EarOuts] at h3
    obtain ⟨-, -, hTA⟩ := h3 (by decide +kernel)
    refine absurd (hTA 2 1 2 (by decide +kernel) (by decide +kernel)
      (Relation.ReflTransGen.single (by decide +kernel)) (fun hb => ?_) (Or.inl rfl) (Or.inl rfl))
      (by decide)
    have := Items.Below.eq_of_no_parent (noParent_of_ball (by decide +kernel)) hb
    exact absurd this (by decide +kernel)


theorem walkTree_book_false :
    ∃ (t : DfsTree) (s : WalkState), t.WF [] ∧ t.Ends s.g ∧
      (∀ v ∈ t.verts, v < s.g.nv) ∧ (∀ e ∈ t.edges, e < s.g.ne) ∧
      t.verts.Nodup ∧ t.edges.Nodup ∧
      s.stackVerts.size = s.g.nv ∧ s.stackDir.size = s.g.nv ∧ s.firstOccurrence.size = s.g.nv ∧
      s.tstack = [] ∧ s.Inv' 0 ∧ Shape s ∧
      (∀ v ∈ t.verts, Items.ch s.items (vertItem v) = [] ∧
        ∀ p, ¬ Items.IsParent s.items p (vertItem v)) ∧
      (∀ e ∈ t.edges, Items.ch s.items (edgeItem s.g e) = [] ∧
        ∀ p, ¬ Items.IsParent s.items p (edgeItem s.g e)) ∧
      ¬ BookTree t 0 s := by
  let g : Graph := ⟨4, #[(1, 0), (2, 1), (2, 3)]⟩
  let t : DfsTree := .node 0 [.tree 0 .bridge (.node 1 [.tree 1 .bridge (.node 2 [])])]
  refine ⟨t, WalkState.init g false, ?_, ?_, ?_, ?_, ?_, ?_, rfl, rfl, rfl, rfl, init_inv g false,
    init_shape g false, ?_, ?_, ?_⟩
  · unfold DfsTree.WF
    refine ⟨List.pairwise_singleton .., fun o ho => ?_⟩
    simp only [List.mem_singleton] at ho; subst ho
    unfold DfsOut.WF DfsTree.WF
    refine ⟨⟨List.pairwise_singleton .., fun o ho => ?_⟩, rfl⟩
    simp only [List.mem_singleton] at ho; subst ho
    unfold DfsOut.WF DfsTree.WF
    exact ⟨⟨List.Pairwise.nil, fun o ho => nomatch ho⟩, rfl⟩
  · unfold DfsTree.Ends
    intro o ho
    simp only [List.mem_singleton] at ho; subst ho
    unfold DfsOut.Ends DfsTree.Ends
    refine ⟨Or.inl rfl, fun o ho => ?_⟩
    simp only [List.mem_singleton] at ho; subst ho
    unfold DfsOut.Ends DfsTree.Ends
    exact ⟨Or.inl rfl, fun o ho => nomatch ho⟩
  · decide
  · decide
  · decide
  · decide
  · intro v _
    refine ⟨Items.initialItems_ch g _, fun p h => ?_⟩
    rw [Items.IsParent, show (init g false).items = initialItems g from rfl,
      Items.initialItems_ch] at h
    exact absurd h (List.not_mem_nil)
  · intro e _
    refine ⟨Items.initialItems_ch g _, fun p h => ?_⟩
    rw [Items.IsParent, show (init g false).items = initialItems g from rfl,
      Items.initialItems_ch] at h
    exact absurd h (List.not_mem_nil)
  · intro h
    unfold BookTree BookOuts BookOut at h
    have h2 := h.1.2
    simp only [WalkM.wp] at h2
    unfold BookTree BookOuts at h2
    have h3 := h2.1.2
    simp only [WalkM.wp, BookOuts] at h3
    obtain ⟨-, -, hTA⟩ := h3 (by decide +kernel)
    refine absurd (hTA 2 1 2 (by decide +kernel) (by decide +kernel)
      (Relation.ReflTransGen.single (by decide +kernel)) (fun hb => ?_) (Or.inl rfl) (Or.inl rfl))
      (by decide)
    have := Items.Below.eq_of_no_parent (noParent_of_ball (by decide +kernel)) hb
    exact absurd this (by decide +kernel)


#print axioms walkTree_ear_false
#print axioms walkTree_book_false

end WalkState
end Spqr
