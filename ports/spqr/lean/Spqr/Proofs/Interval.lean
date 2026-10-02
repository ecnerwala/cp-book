import Spqr.Proofs.Postorder
import Spqr.Proofs.SepPair
import Spqr.Proofs.Type2

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

open DfsData in
/-- Out-edges of a vertex of `T_{a'}` lie in the postorder of the subtree at `a'`. -/
theorem ofForest_mem_edgePostorder_iff {t c : DfsTree} (ht : t ∈ forest) (hc : c.Sub t)
    {e' : Nat} :
    e' ∈ c.edgePostorder ↔
      ∃ v, ∃ o ∈ (ofForest forest).outs v, o.e = e' ∧ (ofForest forest).Anc c.v v := by
  constructor
  · intro h
    obtain ⟨s', hs', o, ho, rfl⟩ := c.mem_edgePostorder h
    refine ⟨s'.v, o, ?_, rfl, ofForest_anc_of_sub hs.verts_nodup ht hs' hc⟩
    rw [ofForest_outs hs.verts_nodup ht (hs'.trans hc)]; exact ho
  · rintro ⟨v, o, ho, rfl, hv⟩
    exact ofForest_e_mem_edgePostorder hs.verts_nodup ht hc hv ho

open DfsData in
/-- Fact D for type-2 pairs: with `a = sa.v` not the root, `a' = c.v` its child toward `b = sb.v`,
and the *above*/*between* parts separated by `{a, b}`, the class of `a → a'` is `type2Block`, an
interval of the walk's edge order. -/
theorem type2_class_interval (h2 : g.TwoConnected) {t sa : DfsTree} (ht : t ∈ forest)
    (hsa : sa.Sub t) {q : Nat} (hq : (ofForest forest).IsParent q sa.v)
    {e₀ : Nat} {cls₀ : OutClass} {c : DfsTree} (ho₀ : DfsOut.tree e₀ cls₀ c ∈ sa.outs)
    {sb : DfsTree} (hsb : sb.Sub c) (hne : c.v ≠ sb.v)
    (hsep : ∀ e e', (ofForest forest).Above c.v sa.v sb.v e g →
      (ofForest forest).Between c.v sb.v e' g → ¬g.SepClass sa.v sb.v e e') :
    DfsTree.type2Block ((ofForest forest).depth sa.v) e₀ c sb <:+: edgePostorderForest forest ∧
      ∀ e', g.SepClass sa.v sb.v e₀ e' ↔
        e' ∈ DfsTree.type2Block ((ofForest forest).depth sa.v) e₀ c sb := by
  have hnd := hs.verts_nodup
  have hct : c.Sub t := (DfsTree.Sub.child ho₀).trans hsa
  have hbt : sb.Sub t := hsb.trans hct
  have ho₀' : DfsOut.tree e₀ cls₀ c ∈ (ofForest forest).outs sa.v := by
    rw [ofForest_outs hnd ht hsa]; exact ho₀
  have hpa : (ofForest forest).IsParent sa.v c.v := ⟨_, ho₀', rfl, rfl⟩
  have ha'b : (ofForest forest).Anc c.v sb.v := ofForest_anc_of_sub hnd ht hsb hct
  have hblock := block_infix_forest ht hsa ho₀
  have hndc : c.edgePostorder.Nodup :=
    hs.edges_nodup.sublist ((List.prefix_append _ _).isInfix.trans hblock).sublist
  have hpre : ∀ s : DfsTree, sb.Sub s → s.Sub c → sb.edgePostorder <+: s.edgePostorder := by
    intro s h
    induction h with
    | refl => exact fun _ => List.prefix_refl _
    | @step v outs e cls child hsub hmem ih =>
      intro hsc
      have hchild : child.Sub c := (DfsTree.Sub.child hmem).trans hsc
      refine (ih hchild).trans ?_
      have houts : (ofForest forest).outs v = outs := ofForest_outs hnd ht (hsc.trans hct)
      obtain ⟨rest, hrest⟩ := List.eq_cons_of_mem_of_not_sublist hmem fun o hsl =>
        type2_first_out hs.toSpec h2 hq hpa ha'b hsep (o := o) (o' := DfsOut.tree e cls child)
          (ofForest_anc_of_sub hnd ht hsc hct) (houts ▸ hsl) rfl
          (ofForest_anc_of_sub hnd ht hsub (hchild.trans hct))
      rw [hrest]
      exact List.prefix_append _ _
  obtain ⟨R, hR⟩ := hpre c hsb (.refl c)
  have hdrop : c.edgePostorder.drop sb.edgePostorder.length = R := by
    rw [← hR, List.drop_left]
  have hsorted : sb.outs.Pairwise fun o o' => o.cls.rank ≤ o'.cls.rank := by
    rw [← ofForest_outs hnd ht hbt]; exact hs.sorted _
  have hiff : ∀ e', g.SepClass sa.v sb.v e₀ e' ↔ _ :=
    fun e' => type2_class_iff hs.toSpec h2 hpa ha'b hne hsep ho₀' rfl rfl (e' := e')
  unfold DfsTree.type2Block
  rw [hdrop]
  constructor
  · have hsb' : ∀ (p : DfsOut → Bool) (l : List DfsOut), DfsOut.edgePostorderList l =
        DfsOut.edgePostorderList (l.takeWhile p) ++ DfsOut.edgePostorderList (l.dropWhile p) := by
      intro p l
      rw [DfsOut.edgePostorderList_eq_flatMap, DfsOut.edgePostorderList_eq_flatMap,
        DfsOut.edgePostorderList_eq_flatMap, ← List.flatMap_append,
        List.takeWhile_append_dropWhile]
    have hblock' : c.edgePostorder ++ [e₀] <:+: edgePostorderForest forest := hblock
    rw [← hR, sb.edgePostorder_eq, hsb' (fun o =>
      decide (o.cls.rank ≤ (OutClass.ret ((ofForest forest).depth sa.v) .backEdge).rank)),
      List.append_assoc, List.append_assoc] at hblock'
    rw [List.append_assoc]
    exact (List.suffix_append _ _).isInfix.trans hblock'
  · intro e'
    rw [hiff e']
    simp only [List.mem_append, List.mem_singleton]
    have hRmem : e' ∈ R ↔ e' ∈ c.edgePostorder ∧ e' ∉ sb.edgePostorder := by
      rw [← hR] at hndc ⊢
      rw [List.nodup_append] at hndc
      simp only [List.mem_append]
      constructor
      · exact fun h => ⟨.inr h, fun h' => hndc.2.2 _ h' _ h rfl⟩
      · rintro ⟨h | h, h'⟩
        · exact absurd h h'
        · exact h
    have hQ : e' ∈ DfsOut.edgePostorderList (sb.outs.dropWhile fun o =>
        decide (o.cls.rank ≤ (OutClass.ret ((ofForest forest).depth sa.v) .backEdge).rank)) ↔
        ∃ o ∈ (ofForest forest).outs sb.v,
          (OutClass.ret ((ofForest forest).depth sa.v) .backEdge).rank < o.cls.rank ∧
          ((o.isTree = true ∧ (ofForest forest).EndIn o.dest e' g) ∨
            (o.isTree = false ∧ o.e = e')) := by
      rw [DfsOut.edgePostorderList_eq_flatMap, List.mem_flatMap, ofForest_outs hnd ht hbt]
      constructor
      · rintro ⟨o, ho, he⟩
        obtain ⟨ho, hrank⟩ := (List.mem_dropWhile_sorted hsorted).mp ho
        refine ⟨o, ho, hrank, ?_⟩
        cases o with
        | back e dest cls => exact .inr ⟨rfl, (List.mem_singleton.mp he).symm⟩
        | tree e cls child => exact .inl ⟨rfl, (endIn_iff_mem_block hs ht hbt ho).mpr he⟩
      · rintro ⟨o, ho, hrank, h⟩
        refine ⟨o, (List.mem_dropWhile_sorted hsorted).mpr ⟨ho, hrank⟩, ?_⟩
        cases o with
        | back e dest cls =>
          rcases h with ⟨h, -⟩ | ⟨-, rfl⟩
          · cases h
          · exact List.mem_singleton.mpr rfl
        | tree e cls child =>
          rcases h with ⟨-, h⟩ | ⟨h, -⟩
          · exact (endIn_iff_mem_block hs ht hbt ho).mp h
          · cases h
    have hRmem' : e' ∈ R ↔ ∃ v, ∃ o ∈ (ofForest forest).outs v, o.e = e' ∧
        (ofForest forest).Anc c.v v ∧ ¬(ofForest forest).Anc sb.v v := by
      rw [hRmem, ofForest_mem_edgePostorder_iff hs ht hct, ofForest_mem_edgePostorder_iff hs ht hbt]
      constructor
      · rintro ⟨⟨v, o, ho, rfl, hv⟩, hnot⟩
        exact ⟨v, o, ho, rfl, hv, fun hbv => hnot ⟨v, o, ho, rfl, hbv⟩⟩
      · rintro ⟨v, o, ho, rfl, hv, hbv⟩
        refine ⟨⟨v, o, ho, rfl, hv⟩, ?_⟩
        rintro ⟨v', o', ho', he, hbv'⟩
        obtain ⟨rfl, -⟩ := hs.out_inj _ _ _ ho' _ ho he
        exact hbv hbv'
    rw [hQ, hRmem']
    constructor
    · rintro (h | h | h)
      · exact .inr h
      · exact .inl (.inr h)
      · exact .inl (.inl h)
    · rintro ((h | h) | h)
      · exact .inr (.inr h)
      · exact .inr (.inl h)
      · exact .inl h

/-- Fact D, laminarity: the type-1 classes (blocks of tree edges) of all pairs are pairwise nested
or disjoint. -/
theorem type1_classes_laminar {t₁ t₂ : DfsTree} (h₁ : t₁ ∈ forest) (h₂ : t₂ ∈ forest)
    {B₁ B₂ : List Nat} (hB₁ : B₁ ∈ t₁.blocks) (hB₂ : B₂ ∈ t₂.blocks) : Laminar B₁ B₂ :=
  blocks_laminar_forest hs.edges_nodup h₁ h₂ hB₁ hB₂

end Spqr
