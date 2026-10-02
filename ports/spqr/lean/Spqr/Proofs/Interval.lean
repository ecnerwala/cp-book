import Spqr.Proofs.Postorder
import Spqr.Proofs.SepPair

/-!
# Fact D: separation classes as intervals of the walk's edge order

The abstract facts of `Spqr.Proofs.SepPair` are transported to the forest actually computed:
`DfsData.ofForest forest` has the forest's out-lists, its ancestor relation is subtree
membership, and the edges with an endpoint in `T_c` are exactly the block of the tree edge into
`c`.
-/

namespace Spqr

namespace DfsData

variable {forest : List DfsTree} (hnd : (forest.flatMap DfsTree.verts).Nodup)
include hnd

theorem ofForest_anc_of_sub {t s s' : DfsTree} (ht : t ∈ forest) (h : s'.Sub s) (hs : s.Sub t) :
    (ofForest forest).Anc s.v s'.v := by
  induction h with
  | refl => exact Anc.refl _
  | @step v outs e cls child hsub hmem ih =>
    have hchild : child.Sub t := (DfsTree.Sub.step (.refl child) hmem).trans hs
    refine Anc.head ⟨DfsOut.tree e cls child, ?_, rfl, rfl⟩ (ih hchild)
    rw [ofForest_outs hnd ht hs]
    exact hmem

theorem ofForest_mem_verts_of_anc {t s : DfsTree} (ht : t ∈ forest) (hs : s.Sub t) {x : Nat}
    (h : (ofForest forest).Anc s.v x) : x ∈ s.verts := by
  induction h with
  | refl => exact s.v_mem_verts
  | tail _ hp ih =>
    obtain ⟨s', hs', rfl⟩ := s.exists_sub_of_mem_verts ih
    obtain ⟨o, ho, hti, rfl⟩ := hp
    rw [ofForest_outs hnd ht (hs'.trans hs)] at ho
    cases o with
    | back => cases hti
    | tree e cls child =>
      refine hs'.verts_subset ?_
      obtain ⟨v, outs⟩ := s'
      exact List.mem_cons_of_mem _ (DfsOut.mem_vertsList.mpr ⟨_, _, _, ho, child.v_mem_verts⟩)

/-- `T_c` is the vertex set of the subtree at `c`. -/
theorem ofForest_anc_iff {t s : DfsTree} (ht : t ∈ forest) (hs : s.Sub t) {x : Nat} :
    (ofForest forest).Anc s.v x ↔ x ∈ s.verts :=
  ⟨ofForest_mem_verts_of_anc hnd ht hs, fun hx =>
    let ⟨_, h, hv⟩ := s.exists_sub_of_mem_verts hx
    hv ▸ ofForest_anc_of_sub hnd ht h hs⟩

theorem ofForest_e_mem_edgePostorder {t s : DfsTree} (ht : t ∈ forest) (hs : s.Sub t) {x : Nat}
    (hx : (ofForest forest).Anc s.v x) {o : DfsOut} (ho : o ∈ (ofForest forest).outs x) :
    o.e ∈ s.edgePostorder := by
  obtain ⟨s', hs', rfl⟩ := s.exists_sub_of_mem_verts (ofForest_mem_verts_of_anc hnd ht hs hx)
  rw [ofForest_outs hnd ht (hs'.trans hs)] at ho
  exact hs'.e_mem_edgePostorder ho

end DfsData

variable {g : Graph} {forest : List DfsTree} (hs : DfsForestSpec g forest)
include hs

open DfsData in
/-- The edges with an endpoint in `T_c` are the block of the tree edge into `c`. -/
theorem endIn_iff_mem_block {t s : DfsTree} (ht : t ∈ forest) (hsub : s.Sub t) {e : Nat}
    {cls : OutClass} {c : DfsTree} (ho : DfsOut.tree e cls c ∈ s.outs) {e' : Nat} :
    (ofForest forest).EndIn c.v e' g ↔ e' ∈ (DfsOut.tree e cls c).block := by
  have hcsub : c.Sub t := (DfsTree.Sub.child ho).trans hsub
  have ho' : DfsOut.tree e cls c ∈ (ofForest forest).outs s.v := by
    rw [ofForest_outs hs.verts_nodup ht hsub]; exact ho
  have hpar : (ofForest forest).IsParent s.v c.v := ⟨_, ho', rfl, rfl⟩
  simp only [DfsOut.block, List.mem_append, List.mem_singleton]
  constructor
  · rintro ⟨x, ⟨y, hj⟩, hcx⟩
    obtain ⟨v, o', ho'', rfl, ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩⟩ := edge_out_cases hs.toSpec hj
    · exact .inl (ofForest_e_mem_edgePostorder hs.verts_nodup ht hcsub hcx ho'')
    · by_cases hcv : (ofForest forest).Anc c.v v
      · exact .inl (ofForest_e_mem_edgePostorder hs.verts_nodup ht hcsub hcv ho'')
      · cases hti : o'.isTree
        · exact absurd (hcx.trans (hs.back_anc _ _ ho'' hti)) hcv
        · have hp := isParent_of_tree ho'' hti
          have hdest : o'.dest = c.v := IsParent.eq_of_anc hs.toSpec hcx hp hcv
          obtain rfl : v = s.v := hs.parent_unique _ _ _ hp (hdest ▸ hpar)
          exact .inr (congrArg DfsOut.e (hs.tree_inj _ _ ho'' _ ho' hti rfl hdest))
  · rintro (h | rfl)
    · obtain ⟨s', hs', o', ho', rfl⟩ := c.mem_edgePostorder h
      have hmem : o' ∈ (ofForest forest).outs s'.v := by
        rw [ofForest_outs hs.verts_nodup ht (hs'.trans hcsub)]; exact ho'
      exact ⟨s'.v, (hs.joins _ _ hmem).isEnd, ofForest_anc_of_sub hs.verts_nodup ht hs' hcsub⟩
    · exact ⟨c.v, (hs.joins _ _ ho').symm.isEnd, Anc.refl _⟩

open DfsData in
/-- Fact D for type-1 classes: the class of a type-1 child `c` of `b` returning to `a` is the
block of `b → c`, an interval of the walk's edge order. -/
theorem type1_class_interval {t s : DfsTree} (ht : t ∈ forest) (hsub : s.Sub t) {a : Nat}
    (hab : (ofForest forest).Anc a s.v) {e : Nat} {cls : OutClass} {c : DfsTree}
    (ho : DfsOut.tree e cls c ∈ s.outs)
    (hcls : cls = .ret ((ofForest forest).depth a) .type1Child) :
    (DfsOut.tree e cls c).block <:+: edgePostorderForest forest ∧
      ∀ e', g.SepClass a s.v e e' ↔ e' ∈ (DfsOut.tree e cls c).block := by
  have ho' : DfsOut.tree e cls c ∈ (ofForest forest).outs s.v := by
    rw [ofForest_outs hs.verts_nodup ht hsub]; exact ho
  refine ⟨block_infix_forest ht hsub ho, fun e' => ?_⟩
  have he : (ofForest forest).EndIn c.v e g := ⟨c.v, (hs.joins _ _ ho').symm.isEnd, Anc.refl _⟩
  exact (type1_class hs.toSpec hab ho' hcls he).trans (endIn_iff_mem_block hs ht hsub ho)

/-- Fact D, laminarity: the type-1 classes (blocks of tree edges) of all pairs are pairwise nested
or disjoint. -/
theorem type1_classes_laminar {t₁ t₂ : DfsTree} (h₁ : t₁ ∈ forest) (h₂ : t₂ ∈ forest)
    {B₁ B₂ : List Nat} (hB₁ : B₁ ∈ t₁.blocks) (hB₂ : B₂ ∈ t₂.blocks) : Laminar B₁ B₂ :=
  blocks_laminar_forest hs.edges_nodup h₁ h₂ hB₁ hB₂

end Spqr
