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

theorem Laminar.symm {B₁ B₂ : List Nat} (h : Laminar B₁ B₂) : Laminar B₂ B₁ := by
  rcases h with h | h | h
  · exact .inr (.inl h)
  · exact .inl h
  · exact .inr (.inr fun x hx hx' => h x hx' hx)

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
/-- Along the tree path from `a'` to `b` of a type-2 pair every tree edge is first in its
out-list (`type2_first_out`), so `σ(T_b)` is a prefix of `σ(T_{a'})`. -/
theorem type2_chain_prefix (h2 : g.TwoConnected) {t sa : DfsTree} (ht : t ∈ forest)
    (hsa : sa.Sub t) {q : Nat} (hq : (ofForest forest).IsParent q sa.v)
    {e₀ : Nat} {cls₀ : OutClass} {c : DfsTree} (ho₀ : DfsOut.tree e₀ cls₀ c ∈ sa.outs)
    {sb : DfsTree} (hsb : sb.Sub c)
    (hsep : ∀ e e', (ofForest forest).Above c.v sa.v sb.v e g →
      (ofForest forest).Between c.v sb.v e' g → ¬g.SepClass sa.v sb.v e e') :
    sb.edgePostorder <+: c.edgePostorder := by
  have hnd := hs.verts_nodup
  have hct : c.Sub t := (DfsTree.Sub.child ho₀).trans hsa
  have ho₀' : DfsOut.tree e₀ cls₀ c ∈ (ofForest forest).outs sa.v := by
    rw [ofForest_outs hnd ht hsa]; exact ho₀
  have hpa : (ofForest forest).IsParent sa.v c.v := ⟨_, ho₀', rfl, rfl⟩
  have ha'b : (ofForest forest).Anc c.v sb.v := ofForest_anc_of_sub hnd ht hsb hct
  suffices h : ∀ s : DfsTree, sb.Sub s → s.Sub c → sb.edgePostorder <+: s.edgePostorder from
    h c hsb (.refl c)
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
  obtain ⟨R, hR⟩ := type2_chain_prefix hs h2 ht hsa hq ho₀ hsb hsep
  have hdrop : c.edgePostorder.drop sb.edgePostorder.length = R := by
    rw [← hR, List.drop_left]
  have hsorted : sb.outs.Pairwise fun o o' => o.cls.rank ≤ o'.cls.rank := by
    rw [← ofForest_outs hnd ht hbt]; exact hs.sorted _
  have hiff : ∀ e', g.SepClass sa.v sb.v e₀ e' ↔ _ :=
    fun e' => type2_class_iff hs.toSpec h2 hpa ha'b hne hsep ho₀' rfl rfl (e' := e')
  unfold DfsTree.type2Block
  rw [hdrop]
  constructor
  · have hblock' : c.edgePostorder ++ [e₀] <:+: edgePostorderForest forest := hblock
    rw [← hR, sb.edgePostorder_eq, DfsOut.edgePostorderList_takeWhile_dropWhile (fun o =>
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

open DfsData in
theorem ofForest_anc_of_e_mem_edgePostorder {t s : DfsTree} (ht : t ∈ forest) (hst : s.Sub t)
    {v : Nat} {o : DfsOut} (ho : o ∈ (ofForest forest).outs v) (h : o.e ∈ s.edgePostorder) :
    (ofForest forest).Anc s.v v := by
  obtain ⟨v', o', ho', he, hv'⟩ := (ofForest_mem_edgePostorder_iff hs ht hst).mp h
  obtain ⟨rfl, -⟩ := hs.out_inj _ _ _ ho' _ ho he
  exact hv'

open DfsData in
theorem edgePostorder_subset_of_anc {t s t' s' : DfsTree} (ht : t ∈ forest) (hst : s.Sub t)
    (ht' : t' ∈ forest) (hs't' : s'.Sub t') (h : (ofForest forest).Anc s.v s'.v) :
    s'.edgePostorder ⊆ s.edgePostorder := by
  intro x hx
  obtain ⟨v, o, ho, rfl, hv⟩ := (ofForest_mem_edgePostorder_iff hs ht' hs't').mp hx
  exact (ofForest_mem_edgePostorder_iff hs ht hst).mpr ⟨v, o, ho, rfl, h.trans hv⟩

open DfsData in
theorem edgePostorder_disjoint_of_not_anc {t s t' s' : DfsTree} (ht : t ∈ forest) (hst : s.Sub t)
    (ht' : t' ∈ forest) (hs't' : s'.Sub t') (h : ¬(ofForest forest).Anc s.v s'.v)
    (h' : ¬(ofForest forest).Anc s'.v s.v) {x : Nat} (hx : x ∈ s.edgePostorder)
    (hx' : x ∈ s'.edgePostorder) : False := by
  obtain ⟨v, o, ho, rfl, hv⟩ := (ofForest_mem_edgePostorder_iff hs ht hst).mp hx
  have hv' := ofForest_anc_of_e_mem_edgePostorder hs ht' hs't' ho hx'
  rcases hv.comparable hs.toSpec hv' with h₁ | h₁
  · exact h h₁
  · exact h' h₁

open DfsData in
/-- A vertex is entered by one tree edge: equal destinations give equal tree out-edges. -/
theorem tree_out_eq_of_v_eq {v v' : Nat} {e e' : Nat} {cls cls' : OutClass} {w w' : DfsTree}
    (ho : DfsOut.tree e cls w ∈ (ofForest forest).outs v)
    (ho' : DfsOut.tree e' cls' w' ∈ (ofForest forest).outs v') (h : w.v = w'.v) :
    v = v' ∧ e = e' ∧ w = w' := by
  obtain rfl : v = v' := hs.parent_unique _ _ _ ⟨_, ho, rfl, h⟩ ⟨_, ho', rfl, rfl⟩
  have := hs.tree_inj _ _ ho _ ho' rfl rfl h
  cases this
  exact ⟨rfl, rfl, rfl⟩

open DfsData in
/-- Fact D, laminarity (type-2 vs type-1): the type-2 class `type2Block` of `{a, b}` (with `a'` the
child of `a` toward `b`) and the block of any tree edge `p → w` are nested or disjoint, unless `w`
lies strictly below `a'` on the tree path to `b` — the S caveat: such a block cuts the same cycle
at another place. `hpre` is `type2_chain_prefix`. -/
theorem type2Block_laminar_block {t sa : DfsTree} (ht : t ∈ forest) (hsa : sa.Sub t)
    {e₀ : Nat} {cls₀ : OutClass} {c : DfsTree} (ho₀ : DfsOut.tree e₀ cls₀ c ∈ sa.outs)
    {sb : DfsTree} (hsb : sb.Sub c) (hpre : sb.edgePostorder <+: c.edgePostorder) {l : Nat}
    {t' s' : DfsTree} (ht' : t' ∈ forest) (hs' : s'.Sub t')
    {e : Nat} {cls : OutClass} {w : DfsTree} (ho : DfsOut.tree e cls w ∈ s'.outs)
    (hchain : ¬((ofForest forest).Anc c.v w.v ∧ w.v ≠ c.v ∧ (ofForest forest).Anc w.v sb.v)) :
    Laminar (DfsTree.type2Block l e₀ c sb) (DfsOut.tree e cls w).block := by
  have hnd := hs.verts_nodup
  have hct : c.Sub t := (DfsTree.Sub.child ho₀).trans hsa
  have hbt : sb.Sub t := hsb.trans hct
  have hwt : w.Sub t' := (DfsTree.Sub.child ho).trans hs'
  have ho₀' : DfsOut.tree e₀ cls₀ c ∈ (ofForest forest).outs sa.v := by
    rw [ofForest_outs hnd ht hsa]; exact ho₀
  have ho' : DfsOut.tree e cls w ∈ (ofForest forest).outs s'.v := by
    rw [ofForest_outs hnd ht' hs']; exact ho
  have hpa : (ofForest forest).IsParent sa.v c.v := ⟨_, ho₀', rfl, rfl⟩
  have hpw : (ofForest forest).IsParent s'.v w.v := ⟨_, ho', rfl, rfl⟩
  have hndc : (c.edgePostorder ++ [e₀]).Nodup :=
    hs.edges_nodup.sublist (block_infix_forest ht hsa ho₀).sublist
  have hmemT : ∀ x, x ∈ DfsTree.type2Block l e₀ c sb ↔ _ := fun x =>
    DfsTree.mem_type2Block (List.nodup_append.mp hndc).1 hpre (l := l) (e₀ := e₀) (x := x)
  have hmemB : ∀ x, x ∈ (DfsOut.tree e cls w).block ↔ x ∈ w.edgePostorder ∨ x = e := fun x => by
    show x ∈ w.edgePostorder ++ [e] ↔ _
    rw [List.mem_append, List.mem_singleton]
  have hdrop_sub : DfsOut.edgePostorderList
      (sb.outs.dropWhile fun o => decide (o.cls.rank ≤ (OutClass.ret l .backEdge).rank)) ⊆
        sb.edgePostorder := by
    rw [sb.edgePostorder_eq]; exact DfsOut.edgePostorderList_subset_of_sublist (List.dropWhile_sublist _)
  have hsbc : sb.edgePostorder ⊆ c.edgePostorder := hpre.subset
  have hTsub : ∀ x, x ∈ DfsTree.type2Block l e₀ c sb → x = e₀ ∨ x ∈ c.edgePostorder := by
    intro x hx
    rcases (hmemT x).mp hx with h | ⟨h, -⟩ | h
    · exact .inl h
    · exact .inr h
    · exact .inr (hsbc (hdrop_sub h))
  have he₀_mem : ∀ {s : DfsTree} {t'' : DfsTree}, t'' ∈ forest → s.Sub t'' →
      (ofForest forest).Anc s.v sa.v → e₀ ∈ s.edgePostorder := fun ht'' hst'' h =>
    (ofForest_mem_edgePostorder_iff hs ht'' hst'').mpr ⟨_, _, ho₀', rfl, h⟩
  have he_mem : ∀ {s : DfsTree} {t'' : DfsTree}, t'' ∈ forest → s.Sub t'' →
      (ofForest forest).Anc s.v s'.v → e ∈ s.edgePostorder := fun ht'' hst'' h =>
    (ofForest_mem_edgePostorder_iff hs ht'' hst'').mpr ⟨_, _, ho', rfl, h⟩
  by_cases hwc : (ofForest forest).Anc w.v c.v
  · -- `T_{a'} ⊆ T_w`: the type-2 class is nested in the block of `p → w`
    left
    intro x hx
    rw [hmemB]
    by_cases hvv : w.v = c.v
    · obtain ⟨-, rfl, rfl⟩ := tree_out_eq_of_v_eq hs ho' ho₀' hvv
      rcases hTsub x hx with rfl | h
      · exact .inr rfl
      · exact .inl h
    · rcases hTsub x hx with rfl | h
      · exact .inl (he₀_mem ht' hwt (hwc.to_parent hs.toSpec hvv hpa))
      · exact .inl (edgePostorder_subset_of_anc hs ht' hwt ht hct hwc h)
  by_cases hcw : (ofForest forest).Anc c.v w.v
  · have hvv : w.v ≠ c.v := fun h => hwc (h ▸ Anc.refl _)
    have hcs' : (ofForest forest).Anc c.v s'.v := hcw.to_parent hs.toSpec (Ne.symm hvv) hpw
    have hwb : ¬(ofForest forest).Anc w.v sb.v := fun h => hchain ⟨hcw, hvv, h⟩
    have hBc : ∀ x, x ∈ (DfsOut.tree e cls w).block → x ∈ c.edgePostorder := by
      intro x hx
      rcases (hmemB x).mp hx with h | rfl
      · exact edgePostorder_subset_of_anc hs ht hct ht' hwt hcw h
      · exact he_mem ht hct hcs'
    by_cases hbw : (ofForest forest).Anc sb.v w.v
    · -- `w` is strictly below `b`: `block(p → w)` lies in the block of one out-edge of `b`
      have hwb' : w.v ≠ sb.v := fun h => hwb (h ▸ Anc.refl _)
      obtain ⟨c', hpc', hc'w⟩ := hbw.child (Ne.symm hwb')
      obtain ⟨o₁, ho₁, ht₁, rfl⟩ := hpc'
      have ho₁' := ho₁
      rw [ofForest_outs hnd ht hbt] at ho₁'
      have hBo : ∀ x, x ∈ (DfsOut.tree e cls w).block → x ∈ o₁.block := by
        cases o₁ with
        | back => cases ht₁
        | tree e₁ cls₁ c₁ =>
          intro x hx
          have hc₁t : c₁.Sub t := (DfsTree.Sub.child ho₁').trans hbt
          show x ∈ c₁.edgePostorder ++ [e₁]
          rw [List.mem_append, List.mem_singleton]
          by_cases hc₁w : c₁.v = w.v
          · obtain ⟨-, rfl, rfl⟩ := tree_out_eq_of_v_eq hs ho₁ ho' hc₁w
            rcases (hmemB x).mp hx with h | rfl
            · exact .inl h
            · exact .inr rfl
          · rcases (hmemB x).mp hx with h | rfl
            · exact .inl (edgePostorder_subset_of_anc hs ht hc₁t ht' hwt hc'w h)
            · exact .inl (he_mem ht hc₁t (hc'w.to_parent hs.toSpec hc₁w hpw))
      rcases List.mem_append.mp ((List.takeWhile_append_dropWhile
          (p := fun o : DfsOut => decide (o.cls.rank ≤ (OutClass.ret l .backEdge).rank))
          (l := sb.outs)).symm ▸
            ho₁') with htake | hdrop
      · right; right
        intro x hxT hxB
        exact DfsTree.type2Block_disjoint_takeWhile hndc hpre
          ((DfsOut.block_infix_of_mem htake).subset (hBo x hxB)) hxT
      · right; left
        intro x hx
        exact (hmemT x).mpr (.inr (.inr ((DfsOut.block_infix_of_mem hdrop).subset (hBo x hx))))
    · -- `w` is off the path to `b`: `block(p → w) ⊆ σ(T_{a'}) − σ(T_b)`
      right; left
      intro x hx
      refine (hmemT x).mpr (.inr (.inl ⟨hBc x hx, fun hxb => ?_⟩))
      rcases (hmemB x).mp hx with h | rfl
      · exact edgePostorder_disjoint_of_not_anc hs ht' hwt ht hbt hwb hbw h hxb
      · exact hbw ((ofForest_anc_of_e_mem_edgePostorder hs ht hbt ho' hxb).tail hpw)
  · -- `T_w` and `T_{a'}` are disjoint
    right; right
    intro x hxT hxB
    rcases hTsub x hxT with rfl | hxc
    · rcases (hmemB x).mp hxB with h | rfl
      · exact hwc ((ofForest_anc_of_e_mem_edgePostorder hs ht' hwt ho₀' h).tail hpa)
      · obtain ⟨-, heq⟩ := hs.out_inj _ _ _ ho' _ ho₀' rfl
        cases heq
        exact hwc (Anc.refl _)
    · rcases (hmemB x).mp hxB with h | rfl
      · exact edgePostorder_disjoint_of_not_anc hs ht hct ht' hwt hcw hwc hxc h
      · exact hcw ((ofForest_anc_of_e_mem_edgePostorder hs ht hct ho' hxc).tail hpw)

/-- Fact D, laminarity: the type-1 classes (blocks of tree edges) of all pairs are pairwise nested
or disjoint. -/
theorem type1_classes_laminar {t₁ t₂ : DfsTree} (h₁ : t₁ ∈ forest) (h₂ : t₂ ∈ forest)
    {B₁ B₂ : List Nat} (hB₁ : B₁ ∈ t₁.blocks) (hB₂ : B₂ ∈ t₂.blocks) : Laminar B₁ B₂ :=
  blocks_laminar_forest hs.edges_nodup h₁ h₂ hB₁ hB₂

open DfsData in
/-- Type-2 vs type-2, `T_{a₂'} ∌ a₁'`: the class of `{a₁, b₁}` is laminar with the block of
`a₂ → a₂'` (unless `a₂'` is strictly below `a₁'` on the path to `b₁`), and the class of `{a₂, b₂}`
is contained in that block but the class of `{a₁, b₁}` is not. -/
theorem type2Block_laminar_type2Block_aux (h2 : g.TwoConnected)
    {t₁ sa₁ : DfsTree} (ht₁ : t₁ ∈ forest) (hsa₁ : sa₁.Sub t₁)
    {q₁ : Nat} (hq₁ : (ofForest forest).IsParent q₁ sa₁.v)
    {e₁ : Nat} {cls₁ : OutClass} {c₁ : DfsTree} (ho₁ : DfsOut.tree e₁ cls₁ c₁ ∈ sa₁.outs)
    {sb₁ : DfsTree} (hsb₁ : sb₁.Sub c₁)
    (hsep₁ : ∀ e e', (ofForest forest).Above c₁.v sa₁.v sb₁.v e g →
      (ofForest forest).Between c₁.v sb₁.v e' g → ¬g.SepClass sa₁.v sb₁.v e e')
    {t₂ sa₂ : DfsTree} (ht₂ : t₂ ∈ forest) (hsa₂ : sa₂.Sub t₂)
    {q₂ : Nat} (hq₂ : (ofForest forest).IsParent q₂ sa₂.v)
    {e₂ : Nat} {cls₂ : OutClass} {c₂ : DfsTree} (ho₂ : DfsOut.tree e₂ cls₂ c₂ ∈ sa₂.outs)
    {sb₂ : DfsTree} (hsb₂ : sb₂.Sub c₂)
    (hsep₂ : ∀ e e', (ofForest forest).Above c₂.v sa₂.v sb₂.v e g →
      (ofForest forest).Between c₂.v sb₂.v e' g → ¬g.SepClass sa₂.v sb₂.v e e')
    (hchain : ¬((ofForest forest).Anc c₁.v c₂.v ∧ c₂.v ≠ c₁.v ∧ (ofForest forest).Anc c₂.v sb₁.v))
    (hnot : ¬(ofForest forest).Anc c₂.v c₁.v) :
    Laminar (DfsTree.type2Block ((ofForest forest).depth sa₁.v) e₁ c₁ sb₁)
      (DfsTree.type2Block ((ofForest forest).depth sa₂.v) e₂ c₂ sb₂) := by
  have hnd := hs.verts_nodup
  have hpre₁ := type2_chain_prefix hs h2 ht₁ hsa₁ hq₁ ho₁ hsb₁ hsep₁
  have hpre₂ := type2_chain_prefix hs h2 ht₂ hsa₂ hq₂ ho₂ hsb₂ hsep₂
  have hlam := type2Block_laminar_block hs ht₁ hsa₁ ho₁ hsb₁ hpre₁
    (l := (ofForest forest).depth sa₁.v) ht₂ hsa₂ ho₂ hchain
  have hc₂t : c₂.Sub t₂ := (DfsTree.Sub.child ho₂).trans hsa₂
  have ho₁' : DfsOut.tree e₁ cls₁ c₁ ∈ (ofForest forest).outs sa₁.v := by
    rw [ofForest_outs hnd ht₁ hsa₁]; exact ho₁
  have ho₂' : DfsOut.tree e₂ cls₂ c₂ ∈ (ofForest forest).outs sa₂.v := by
    rw [ofForest_outs hnd ht₂ hsa₂]; exact ho₂
  have hndc₂ : (c₂.edgePostorder ++ [e₂]).Nodup :=
    hs.edges_nodup.sublist (block_infix_forest ht₂ hsa₂ ho₂).sublist
  have hsub : DfsTree.type2Block ((ofForest forest).depth sa₂.v) e₂ c₂ sb₂ ⊆
      (DfsOut.tree e₂ cls₂ c₂).block := by
    intro x hx
    show x ∈ c₂.edgePostorder ++ [e₂]
    rw [List.mem_append, List.mem_singleton]
    rcases (DfsTree.mem_type2Block (List.nodup_append.mp hndc₂).1 hpre₂).mp hx with h | ⟨h, -⟩ | h
    · exact .inr h
    · exact .inl h
    · refine .inl (hpre₂.subset ?_)
      rw [sb₂.edgePostorder_eq]
      exact DfsOut.edgePostorderList_subset_of_sublist (List.dropWhile_sublist _) h
  have he₁ : e₁ ∈ DfsTree.type2Block ((ofForest forest).depth sa₁.v) e₁ c₁ sb₁ := by
    unfold DfsTree.type2Block; simp
  have hne₁ : e₁ ∉ (DfsOut.tree e₂ cls₂ c₂).block := by
    intro h
    rcases List.mem_append.mp h with h | h
    · exact hnot ((ofForest_anc_of_e_mem_edgePostorder hs ht₂ hc₂t ho₁' h).tail ⟨_, ho₁', rfl, rfl⟩)
    · obtain ⟨-, heq⟩ := hs.out_inj _ _ _ ho₁' _ ho₂' (List.mem_singleton.mp h)
      cases heq
      exact hnot (Anc.refl _)
  rcases hlam with h | h | h
  · exact absurd (h he₁) hne₁
  · exact .inr (.inl fun x hx => h (hsub hx))
  · exact .inr (.inr fun x hx hx' => h x hx (hsub hx'))

open DfsData in
/-- Type-2 vs type-2 with the same `a, a'` and `b₂` weakly below `b₁`: the class of `{a, b₁}` is
contained in the class of `{a, b₂}`. The out-edges of `b₁` sorted after `ret l backEdge` do not
return above `a`, while `T_{b₂}` does (`subtree_returns_above`), so none of them contains `b₂`. -/
theorem type2Block_subset_type2Block (h2 : g.TwoConnected)
    {t₁ sa₁ : DfsTree} (ht₁ : t₁ ∈ forest) (hsa₁ : sa₁.Sub t₁)
    {q₁ : Nat} (hq₁ : (ofForest forest).IsParent q₁ sa₁.v)
    {e₁ : Nat} {cls₁ : OutClass} {c₁ : DfsTree} (ho₁ : DfsOut.tree e₁ cls₁ c₁ ∈ sa₁.outs)
    {sb₁ : DfsTree} (hsb₁ : sb₁.Sub c₁) (hne₁ : c₁.v ≠ sb₁.v)
    (hsep₁ : ∀ e e', (ofForest forest).Above c₁.v sa₁.v sb₁.v e g →
      (ofForest forest).Between c₁.v sb₁.v e' g → ¬g.SepClass sa₁.v sb₁.v e e')
    {t₂ sa₂ : DfsTree} (ht₂ : t₂ ∈ forest) (hsa₂ : sa₂.Sub t₂)
    {q₂ : Nat} (hq₂ : (ofForest forest).IsParent q₂ sa₂.v)
    {e₂ : Nat} {cls₂ : OutClass} {c₂ : DfsTree} (ho₂ : DfsOut.tree e₂ cls₂ c₂ ∈ sa₂.outs)
    {sb₂ : DfsTree} (hsb₂ : sb₂.Sub c₂)
    (hsep₂ : ∀ e e', (ofForest forest).Above c₂.v sa₂.v sb₂.v e g →
      (ofForest forest).Between c₂.v sb₂.v e' g → ¬g.SepClass sa₂.v sb₂.v e e')
    (hvv : c₁.v = c₂.v) (hb : (ofForest forest).Anc sb₁.v sb₂.v) :
    DfsTree.type2Block ((ofForest forest).depth sa₁.v) e₁ c₁ sb₁ ⊆
      DfsTree.type2Block ((ofForest forest).depth sa₂.v) e₂ c₂ sb₂ := by
  have hnd := hs.verts_nodup
  have hpre₁ := type2_chain_prefix hs h2 ht₁ hsa₁ hq₁ ho₁ hsb₁ hsep₁
  have hpre₂ := type2_chain_prefix hs h2 ht₂ hsa₂ hq₂ ho₂ hsb₂ hsep₂
  have ho₁' : DfsOut.tree e₁ cls₁ c₁ ∈ (ofForest forest).outs sa₁.v := by
    rw [ofForest_outs hnd ht₁ hsa₁]; exact ho₁
  have ho₂' : DfsOut.tree e₂ cls₂ c₂ ∈ (ofForest forest).outs sa₂.v := by
    rw [ofForest_outs hnd ht₂ hsa₂]; exact ho₂
  obtain ⟨hv, rfl, rfl⟩ := tree_out_eq_of_v_eq hs ho₁' ho₂' hvv
  have hl : (ofForest forest).depth sa₁.v = (ofForest forest).depth sa₂.v := by rw [hv]
  have hc₁t : c₁.Sub t₁ := (DfsTree.Sub.child ho₁).trans hsa₁
  have hc₂t : c₁.Sub t₂ := (DfsTree.Sub.child ho₂).trans hsa₂
  have hb₁t : sb₁.Sub t₁ := hsb₁.trans hc₁t
  have hb₂t : sb₂.Sub t₂ := hsb₂.trans hc₂t
  have hpa₂ : (ofForest forest).IsParent sa₂.v c₁.v := ⟨_, ho₂', rfl, rfl⟩
  have ha'b₁ : (ofForest forest).Anc c₁.v sb₁.v := ofForest_anc_of_sub hnd ht₁ hsb₁ hc₁t
  have ha'b₂ : (ofForest forest).Anc c₁.v sb₂.v := ofForest_anc_of_sub hnd ht₂ hsb₂ hc₂t
  have hndc₁ : c₁.edgePostorder.Nodup :=
    hs.edges_nodup.sublist ((List.prefix_append _ _).isInfix.trans
      (block_infix_forest ht₁ hsa₁ ho₁)).sublist
  have hsorted₁ : sb₁.outs.Pairwise fun o o' => o.cls.rank ≤ o'.cls.rank := by
    rw [← ofForest_outs hnd ht₁ hb₁t]; exact hs.sorted _
  have hdrop_sub₁ : DfsOut.edgePostorderList (sb₁.outs.dropWhile fun o =>
      decide (o.cls.rank ≤ (OutClass.ret ((ofForest forest).depth sa₁.v) .backEdge).rank)) ⊆
        sb₁.edgePostorder := by
    rw [sb₁.edgePostorder_eq]
    exact DfsOut.edgePostorderList_subset_of_sublist (List.dropWhile_sublist _)
  intro x hx
  rw [DfsTree.mem_type2Block hndc₁ hpre₁] at hx
  rw [DfsTree.mem_type2Block hndc₁ hpre₂]
  rcases hx with rfl | ⟨hxc, hxb₁⟩ | hx
  · exact .inl rfl
  · exact .inr (.inl ⟨hxc, fun h => hxb₁ (edgePostorder_subset_of_anc hs ht₁ hb₁t ht₂ hb₂t hb h)⟩)
  by_cases hxb₂ : x ∈ sb₂.edgePostorder
  · by_cases hbb : sb₁.v = sb₂.v
    · right; right
      have houts : sb₁.outs = sb₂.outs := by
        rw [← ofForest_outs hnd ht₁ hb₁t, ← ofForest_outs hnd ht₂ hb₂t, hbb]
      rw [← houts, ← hl]; exact hx
    exfalso
    rw [DfsOut.edgePostorderList_eq_flatMap] at hx
    obtain ⟨o, ho, hxo⟩ := List.mem_flatMap.mp hx
    obtain ⟨ho, hrank⟩ := (List.mem_dropWhile_sorted hsorted₁).mp ho
    have ho' : o ∈ (ofForest forest).outs sb₁.v := by rw [ofForest_outs hnd ht₁ hb₁t]; exact ho
    obtain ⟨v, o', ho'', rfl, hb₂v⟩ := (ofForest_mem_edgePostorder_iff hs ht₂ hb₂t).mp hxb₂
    cases o with
    | back e dest cls =>
      obtain ⟨rfl, -⟩ := hs.out_inj _ _ _ ho'' _ ho' (List.mem_singleton.mp hxo)
      exact hbb (hb.antisymm hs.toSpec hb₂v)
    | tree e' cls' c' =>
      rcases List.mem_append.mp hxo with hxc' | hxe
      · have hc't : c'.Sub t₁ := (DfsTree.Sub.child ho).trans hb₁t
        have hc'v : (ofForest forest).Anc c'.v v :=
          ofForest_anc_of_e_mem_edgePostorder hs ht₁ hc't ho'' hxc'
        have hpc' : (ofForest forest).IsParent sb₁.v c'.v := ⟨_, ho', rfl, rfl⟩
        have hc'b₂ : (ofForest forest).Anc c'.v sb₂.v := by
          rcases hc'v.comparable hs.toSpec hb₂v with h | h
          · exact h
          · have hd := hs.depth_parent _ _ hpc'
            have h1 := h.depth_le hs.toSpec
            have h2 := hb.depth_lt hs.toSpec hbb
            have heq := h.eq_of_depth_eq hs.toSpec (Anc.refl _) (by omega)
            rw [heq]
            exact Anc.refl _
        obtain ⟨l'', hl'', hret⟩ :=
          subtree_returns_above hs.toSpec h2 hq₂ hpa₂ ha'b₂ hsep₂ ha'b₂ (Anc.refl _)
        have hret' : (ofForest forest).Returns c'.v l'' := by
          obtain ⟨u, o, ho, hu, hback, hl⟩ := hret
          exact ⟨u, o, ho, hc'b₂.trans hu, hback, hl⟩
        obtain ⟨p, -, hpb₁⟩ := ha'b₁.cases_tail.resolve_left fun h => hne₁ h.symm
        obtain ⟨lx, kx, hcx, -⟩ := cls_ret_of_tree hs.toSpec h2 hpb₁ ho' rfl
        obtain ⟨-, -, hmin⟩ := ret_lowpt hs.toSpec ho' rfl hcx
        have hk := kx.rank_le
        rw [hcx] at hrank
        simp only [OutClass.rank, RetKind.rank] at hrank hk
        exact hmin l'' (by omega) hret'
      · obtain ⟨rfl, -⟩ := hs.out_inj _ _ _ ho'' _ ho' (List.mem_singleton.mp hxe)
        exact hbb (hb.antisymm hs.toSpec hb₂v)
  · exact .inr (.inl ⟨hpre₁.subset (hdrop_sub₁ hx), hxb₂⟩)

open DfsData in
/-- Fact D, laminarity (type-2 vs type-2): the type-2 classes of two pairs `{a₁, b₁}`, `{a₂, b₂}`
(each under the hypotheses of `type2_class_interval`) are nested or disjoint, unless one pair's
`a'` lies strictly below the other's `a'` on the other's tree path to `b` — the S caveat: both
pairs then cut the same cycle. -/
theorem type2Block_laminar_type2Block (h2 : g.TwoConnected)
    {t₁ sa₁ : DfsTree} (ht₁ : t₁ ∈ forest) (hsa₁ : sa₁.Sub t₁)
    {q₁ : Nat} (hq₁ : (ofForest forest).IsParent q₁ sa₁.v)
    {e₁ : Nat} {cls₁ : OutClass} {c₁ : DfsTree} (ho₁ : DfsOut.tree e₁ cls₁ c₁ ∈ sa₁.outs)
    {sb₁ : DfsTree} (hsb₁ : sb₁.Sub c₁) (hne₁ : c₁.v ≠ sb₁.v)
    (hsep₁ : ∀ e e', (ofForest forest).Above c₁.v sa₁.v sb₁.v e g →
      (ofForest forest).Between c₁.v sb₁.v e' g → ¬g.SepClass sa₁.v sb₁.v e e')
    {t₂ sa₂ : DfsTree} (ht₂ : t₂ ∈ forest) (hsa₂ : sa₂.Sub t₂)
    {q₂ : Nat} (hq₂ : (ofForest forest).IsParent q₂ sa₂.v)
    {e₂ : Nat} {cls₂ : OutClass} {c₂ : DfsTree} (ho₂ : DfsOut.tree e₂ cls₂ c₂ ∈ sa₂.outs)
    {sb₂ : DfsTree} (hsb₂ : sb₂.Sub c₂) (hne₂ : c₂.v ≠ sb₂.v)
    (hsep₂ : ∀ e e', (ofForest forest).Above c₂.v sa₂.v sb₂.v e g →
      (ofForest forest).Between c₂.v sb₂.v e' g → ¬g.SepClass sa₂.v sb₂.v e e')
    (hchain₁ : ¬((ofForest forest).Anc c₁.v c₂.v ∧ c₂.v ≠ c₁.v ∧ (ofForest forest).Anc c₂.v sb₁.v))
    (hchain₂ : ¬((ofForest forest).Anc c₂.v c₁.v ∧ c₁.v ≠ c₂.v ∧ (ofForest forest).Anc c₁.v sb₂.v)) :
    Laminar (DfsTree.type2Block ((ofForest forest).depth sa₁.v) e₁ c₁ sb₁)
      (DfsTree.type2Block ((ofForest forest).depth sa₂.v) e₂ c₂ sb₂) := by
  by_cases hvv : c₁.v = c₂.v
  · by_cases hb : (ofForest forest).Anc sb₁.v sb₂.v
    · exact .inl (type2Block_subset_type2Block hs h2 ht₁ hsa₁ hq₁ ho₁ hsb₁ hne₁ hsep₁
        ht₂ hsa₂ hq₂ ho₂ hsb₂ hsep₂ hvv hb)
    by_cases hb' : (ofForest forest).Anc sb₂.v sb₁.v
    · exact .inr (.inl (type2Block_subset_type2Block hs h2 ht₂ hsa₂ hq₂ ho₂ hsb₂ hne₂ hsep₂
        ht₁ hsa₁ hq₁ ho₁ hsb₁ hsep₁ hvv.symm hb'))
    exfalso
    have hnd := hs.verts_nodup
    have ho₁' : DfsOut.tree e₁ cls₁ c₁ ∈ (ofForest forest).outs sa₁.v := by
      rw [ofForest_outs hnd ht₁ hsa₁]; exact ho₁
    have ho₂' : DfsOut.tree e₂ cls₂ c₂ ∈ (ofForest forest).outs sa₂.v := by
      rw [ofForest_outs hnd ht₂ hsa₂]; exact ho₂
    obtain ⟨hv, rfl, rfl⟩ := tree_out_eq_of_v_eq hs ho₁' ho₂' hvv
    have hc₁t : c₁.Sub t₁ := (DfsTree.Sub.child ho₁).trans hsa₁
    have hc₂t : c₁.Sub t₂ := (DfsTree.Sub.child ho₂).trans hsa₂
    have hpa₁ : (ofForest forest).IsParent sa₁.v c₁.v := ⟨_, ho₁', rfl, rfl⟩
    have hpa₂ : (ofForest forest).IsParent sa₂.v c₁.v := ⟨_, ho₂', rfl, rfl⟩
    have ha'b₁ : (ofForest forest).Anc c₁.v sb₁.v := ofForest_anc_of_sub hnd ht₁ hsb₁ hc₁t
    have ha'b₂ : (ofForest forest).Anc c₁.v sb₂.v := ofForest_anc_of_sub hnd ht₂ hsb₂ hc₂t
    obtain ⟨l'', hl'', u, o, ho, hb₁u, -, hdepth⟩ :=
      subtree_returns_above hs.toSpec h2 hq₁ hpa₁ ha'b₁ hsep₁ ha'b₁ (Anc.refl _)
    have hj := hs.joins _ _ ho
    have hda := hs.depth_parent _ _ hpa₂
    have hdb₂ := ha'b₂.depth_le hs.toSpec
    rw [hv] at hl''
    refine hsep₂ o.e o.e ⟨o.dest, hj.symm.isEnd, ?_, ?_, ?_⟩
      ⟨u, hj.isEnd, ha'b₁.trans hb₁u, ?_⟩ (.inl rfl)
    · intro h; rw [h] at hdepth; omega
    · intro h; rw [h] at hdepth; omega
    · intro h; have := h.depth_le hs.toSpec; omega
    · intro h
      rcases hb₁u.comparable hs.toSpec h with h' | h'
      · exact hb h'
      · exact hb' h'
  · by_cases h21 : (ofForest forest).Anc c₂.v c₁.v
    · have h12 : ¬(ofForest forest).Anc c₁.v c₂.v := fun h => hvv (h.antisymm hs.toSpec h21)
      exact (type2Block_laminar_type2Block_aux hs h2 ht₂ hsa₂ hq₂ ho₂ hsb₂ hsep₂
        ht₁ hsa₁ hq₁ ho₁ hsb₁ hsep₁ hchain₂ h12).symm
    · exact type2Block_laminar_type2Block_aux hs h2 ht₁ hsa₁ hq₁ ho₁ hsb₁ hsep₁
        ht₂ hsa₂ hq₂ ho₂ hsb₂ hsep₂ hchain₁ h21

end Spqr
