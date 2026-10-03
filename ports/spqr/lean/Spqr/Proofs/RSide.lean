import Spqr.Proofs.RInvWalk
import Spqr.RangesSites

/-!
# `RSideTree` by the walk induction

`walkTree_rSide` (the side hypotheses of the child-return induction `rrTree`) by a mutual
induction mirroring `rrTree`/`rrOuts`/`rrOut`: the positional parts are proved here (the
ancestor chain `ancChain_child` with `d + 1 < stackVerts.size` from `ancChain_lt_size`, the
DFS-layer `ret` from `DfsData.Spec.outs_lowval_lt`, the children's `dfs.outs` from `DfsTree.Sub`),
and the R content is isolated at its three sites as the named admissions `rSide_entry_site`
(`EntryR` stability and the parent-start bound at a child entry), `rSide_vertFree_site`
(`VertFree` where a vertex entry may be pushed) and `rSide_finish_content_site`
(`FinishRShape`'s `settled`/`unwrap`/`vert_own` at a tree-edge `finishEdge` site, with the child's
`RWalk`, `FinishGuards`/`FinishBook`/`Frontier` and the parent's chain in hand; `pend`/`ear` are
`rSide_finish_site`'s bookkeeping).
-/

namespace Spqr
open WalkM
namespace WalkState

variable {dfs : DfsData} {s : WalkState}

/-- The depth-`k` ancestors `stackVerts[k]`, `k ≤ d`, together with the child `c` are `d + 2`
distinct vertices of `g`, so `d + 1 < g.nv = stackVerts.size`. -/
theorem ancChain_lt_size {d c : Nat} (hsp : dfs.Spec s.g) (hwf : s.g.WF)
    (hsz : s.stackVerts.size = s.g.nv) (hp : dfs.IsParent s.stackVerts[d]! c)
    (hchain : ∀ k, k ≤ d → dfs.Anc s.stackVerts[k]! s.stackVerts[d]! ∧ dfs.depth s.stackVerts[k]! = k) :
    d + 1 < s.stackVerts.size := by
  have hend : ∀ {e x y : Nat}, s.g.Joins e x y → x < s.g.nv := by
    intro e x y h
    rcases h with h | h
    · obtain ⟨hi, he⟩ := Array.getElem?_eq_some_iff.1 h
      exact (hwf _ (he ▸ Array.getElem_mem hi)).1
    · obtain ⟨hi, he⟩ := Array.getElem?_eq_some_iff.1 h
      exact (hwf _ (he ▸ Array.getElem_mem hi)).2
  let f : Nat → Nat := fun k => if k = d + 1 then c else s.stackVerts[k]!
  have hdep : ∀ k, k ≤ d + 1 → dfs.depth (f k) = k := by
    intro k hk
    rcases Nat.lt_or_eq_of_le hk with hk | rfl
    · simp only [f, Nat.ne_of_lt hk, ↓reduceIte]; exact (hchain k (by omega)).2
    · simp only [f, ↓reduceIte]; rw [hsp.depth_parent _ _ hp, (hchain d (le_refl _)).2]
  have hbound : ∀ k, k ≤ d + 1 → f k < s.g.nv := by
    intro k hk
    rcases Nat.lt_or_eq_of_le hk with hk | rfl
    · simp only [f, Nat.ne_of_lt hk, ↓reduceIte]
      rcases Nat.lt_or_eq_of_le (Nat.le_of_lt_succ hk) with hk' | rfl
      · rcases (hchain k (by omega)).1.cases_head with heq | ⟨w, hw, -⟩
        · exfalso
          have h1 := (hchain k (by omega)).2
          rw [heq, (hchain d (le_refl _)).2] at h1; omega
        · obtain ⟨e, he⟩ := hw.joins hsp; exact hend he
      · obtain ⟨e, he⟩ := hp.joins hsp; exact hend he
    · simp only [f, ↓reduceIte]; obtain ⟨e, he⟩ := hp.joins hsp; exact hend he.symm
  have hnd : ((List.range (d + 2)).map f).Nodup := by
    refine List.Nodup.map_on (fun x hx y hy hxy => ?_) List.nodup_range
    rw [List.mem_range] at hx hy
    have h1 := hdep x (by omega)
    rw [hxy, hdep y (by omega)] at h1
    exact h1.symm
  have hsub : (List.range (d + 2)).map f ⊆ List.range s.g.nv := by
    intro x hx
    obtain ⟨k, hk, rfl⟩ := List.mem_map.1 hx
    rw [List.mem_range] at hk ⊢
    exact hbound k (by omega)
  have hlen := (hnd.subperm hsub).length_le
  simp only [List.length_map, List.length_range] at hlen
  omega

/-- A tree out-edge's class is a tree class (`Spec.cls_selfLoop`/`cls_backEdge`). -/
theorem _root_.Spqr.DfsData.Spec.cls_isTree {g : Graph} {v : Nat} {o : DfsOut} (hsp : dfs.Spec g)
    (hmem : o ∈ dfs.outs v) (ht : o.isTree = true) : o.cls.isTree = true := by
  cases hc : o.cls with
  | selfLoop => have := ((hsp.cls_selfLoop v o hmem).1 hc).1; rw [ht] at this; cases this
  | ret l k =>
    cases k with
    | backEdge => have := ((hsp.cls_backEdge v o hmem l).1 hc).1; rw [ht] at this; cases this
    | type1Child => rfl
    | type2Child => rfl
  | bridge => rfl
  | component => rfl

/-- **Named admission** (R content, PROOF.md §4.5): `REntryContent` at a child entry `(c, d + 1)`
(parent `stackVerts[d]` at depth `d`, its settled top `RInvTop`). Checked field by field
(`rc.entry_stab`, `rc.entry_parent_top`, `checks/WalkInvCheck/RContent.lean`), 0 violations over
seeds 0..3000 × both modes. -/
theorem entry_rContent {d c : Nat} (hi : s.Inv' d) (hs : Shape s)
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hp : dfs.IsParent s.stackVerts[d]! c)
    (hchain : ∀ k, k ≤ d → dfs.Anc s.stackVerts[k]! s.stackVerts[d]! ∧ dfs.depth s.stackVerts[k]! = k)
    (hR : s.RInvTop dfs s.stackVerts[d]! d) : s.REntryContent dfs d c := by
  sorry

/-- **Named admission** (R content, PROOF.md §4.5): `RCloseContent` at the pre-state of a returning
tree-edge `finishEdge` site, from the Ranges site record `CloseBase` (the backbone proves it at
the sites from the primitives). Checked field by field (`rc.settled`/`s2_top`/`vert_own`/`l1_mid`/
`l1_fields`/`l1_top`/`v_*`, `checks/WalkInvCheck/RContent.lean`), 0 violations over seeds
0..3000 × both modes. -/
theorem closeBase_rContent {σ : List Nat} {n D curV d : Nat} {o : DfsOut} {origTstack : Nat}
    {hasVert : Bool} (hcb : CloseBase σ n D curV d o origTstack hasVert s)
    (ht : o.isTree = true) (hlow : o.cls.lowval d < d)
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g) :
    s.RCloseContent dfs curV d o origTstack hasVert := by
  sorry

/-- **Named admission** (R content, PROOF.md §4.5). from `entry_rContent`: at a child entry
`(c, d + 1)` of the walk (parent `stackVerts[d]` at depth `d`, its settled top `RInvTop`), every
`EntryR` entry stays `EntryR` after `stackVerts.set! (d + 1) c` (its `EntryR` does not read
`stackVerts[d + 1]`: `EntryR.set_stackVerts` given `topDepth ≠ d + 1`; note `∀ t ∈ tstack,
topDepth ≤ d` is FALSE — a buried `V y` entry of a finished deeper vertex may remain, cf.
`EarInv.lean`'s dropped `base_top`), and every entry starting at the parent tops out at depth
`≤ d`. Checked at every entry of seeds 0..400 × both modes + 6000 random multigraphs
(`checks/RFinishEdgeCheck.lean` `stab`/`entry` lines, 0 failures). -/
theorem rSide_entry_site {d c : Nat} (hi : s.Inv' d) (hs : Shape s)
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hp : dfs.IsParent s.stackVerts[d]! c)
    (hchain : ∀ k, k ≤ d → dfs.Anc s.stackVerts[k]! s.stackVerts[d]! ∧ dfs.depth s.stackVerts[k]! = k)
    (hR : s.RInvTop dfs s.stackVerts[d]! d) :
    (∀ t ∈ s.tstack, s.EntryR dfs t →
      ({ s with stackVerts := s.stackVerts.set! (d + 1) c } : WalkState).EntryR dfs t) ∧
    (∀ t ∈ s.tstack, t.vStart = s.stackVerts[d]! → t.topDepth ≤ d) := by
  have hrc := entry_rContent hi hs h2 hsp hrt hp hchain hR
  refine ⟨fun t ht hE => ?_, hrc.parent_top⟩
  by_cases hk : t.topDepth = d + 1
  · exact hrc.stab t ht hk hE
  · exact hE.set_stackVerts hk

/-- **Named admission** (R content, PROOF.md §4.5). Exact obligation: while `v`'s vertex entry has
not been pushed (`VertBook v false`: `vertItem v` holds the blocks already closed at `v`), no open
entry owns an edge below `vertItem v` (`VertFree`; the `vert_disj` shape of `EarClose` at the
`finishEdge` sites, here at the `walkOut` entries and the end of `walkOuts`). Checked at every
site of seeds 0..400 × both modes + 6000 random multigraphs (`vertown` line, 0 failures). -/
theorem rSide_vertFree_site {F : List RFrame} {v d : Nat} (hi : s.Inv' d) (hs : Shape s)
    (hb : VertBook v false s) (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hanc : AncChain dfs v d s) (hW : RWalk dfs F v d s) : VertFree v s := by
  sorry

/-- **Named admission** (R content, PROOF.md §4.5). Exact obligation: the content fields of
`FinishRShape` at a tree-edge `finishEdge` site of a non-root vertex `v = stackVerts[d]` (`o` the
tree edge to the child `c`, `B` the stack size before the child, the child walked: `RWalk` with
the parent frame `(v, d, B)` and the child's settled top `RInvTop c (d + 1)`, `FinishGuards`/
`FinishBook`/`Frontier` of the site): `settled` (the entries loop 1 leaves above the base are
`EntryR` at `feS₁`), `unwrap` (type 1 with a vertex entry: `Exempt` of `feS₂`'s next entry),
`vert_own` (`VertFree` at `feP`). Checked at every site of seeds 0..400 × both modes + 6000
random multigraphs (`FinishRShape` lines, 0 failures). -/
theorem rSide_finish_content_site {F : List RFrame} {v d B e c : Nat} {cls : OutClass}
    {couts : List DfsOut} {o : DfsOut} {hasVert : Bool} (ho : o = .tree e cls (.node c couts))
    (hlow : o.cls.lowval d < d) (hmem : o ∈ dfs.outs v)
    (hi : s.Inv' (d + 1)) (hs : Shape s) (hg : FinishGuards d o B hasVert s)
    (hb : FinishBook v d o B hasVert s) (hfr : Frontier (o := o) d B s)
    (hW : RWalk dfs ((v, d, B) :: F) c (d + 1) s) (hanc : AncChain dfs v d s)
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hrc : s.RCloseContent dfs v d o B hasVert) :
    (∀ t ∈ (feS₁ d o s).tstack.tail.take ((feS₁ d o s).tstack.length - 1 - B),
      d ≤ t.topDepth → t.vStart ≠ v → (feS₁ d o s).EntryR dfs t) ∧
    (hasVert = true → o.cls.isType1 = true → Exempt v d (nxtE (feS₂ d o s))) ∧
    (hasVert = false → ∀ t ∈ (feP v d o s).tstack, ∀ e,
      t.edges (feP v d o s).g (feP v d o s).items e →
      ¬ Items.EdgeBelow (feP v d o s).g (feP v d o s).items (vertItem v) e) := by
  refine ⟨hrc.settled, fun hhv h1 => ?_, hrc.vert_own⟩
  obtain ⟨sub, base, -, hE⟩ := hb.ear
  have htree : o.cls.isTree = true := hsp.cls_isTree hmem (by rw [ho]; rfl)
  obtain ⟨c', mid, py, vy, hC⟩ := hE.close htree hlow
  obtain ⟨rfl, -⟩ := hC.type1 h1
  show Exempt v d (feS₂ d o s).tstack.tail.head!
  rw [hC.tstack]
  exact Exempt.of_lt (by show py.topDepth < d; rw [hC.py_top]; exact hlow)

/-- `FinishRShape` at a tree-edge site: `pend` (no entry holds `o.e`) from `EarFinish`'s
`q_root`/`q_free`, `ear` from `FinishBook.ear`, the content fields from
`rSide_finish_content_site`. -/
theorem rSide_finish_site {F : List RFrame} {v d B e c : Nat} {cls : OutClass}
    {couts : List DfsOut} {o : DfsOut} {hasVert : Bool} (ho : o = .tree e cls (.node c couts))
    (hlow : o.cls.lowval d < d) (hmem : o ∈ dfs.outs v)
    (hi : s.Inv' (d + 1)) (hs : Shape s) (hg : FinishGuards d o B hasVert s)
    (hb : FinishBook v d o B hasVert s) (hfr : Frontier (o := o) d B s)
    (hW : RWalk dfs ((v, d, B) :: F) c (d + 1) s) (hanc : AncChain dfs v d s)
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hrc : s.RCloseContent dfs v d o B hasVert) :
    FinishRShape dfs v d o B hasVert s := by
  obtain ⟨hsettled, hunwrap, hvert⟩ :=
    rSide_finish_content_site ho hlow hmem hi hs hg hb hfr hW hanc h2 hsp hrt hrc
  refine ⟨?_, hsettled, hunwrap, hvert, hb.ear, hrc⟩
  intro t ht ⟨i, hi', hbel⟩
  obtain ⟨sub, base, -, hE⟩ := hb.ear
  exact hE.q_free t ht (Items.Below.eq_of_no_parent hE.q_root hbel ▸ hi')

/-- The parent of a non-root walked vertex `v = stackVerts[d]`, `d = dp + 1`, from its chain. -/
theorem AncChain.parent {v dp : Nat} (hanc : AncChain dfs v (dp + 1) s) :
    ∃ p, dfs.IsParent p v := by
  rcases (hanc.2 dp (Nat.le_succ _)).1.cases_tail with heq | ⟨p, -, hp⟩
  · exfalso
    have h1 := (hanc.2 dp (Nat.le_succ _)).2
    have h2 := (hanc.2 (dp + 1) (le_refl _)).2
    rw [hanc.1, heq] at h2; omega
  · exact ⟨p, hanc.1 ▸ hp⟩

abbrev RSTree (dfs : DfsData) (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  ∀ (σ : List Nat) (n dp : Nat), d = dp + 1 →
    s.Inv' dp → Shape s → GuardsTree t d s → BookTree t d s → CbTree σ n t d s →
    s.g.TwoConnected → s.g.WF → dfs.Spec s.g → dfs.Rooted s.g → s.stackVerts.size = s.g.nv →
    (∀ v outs, t = .node v outs → dfs.IsParent s.stackVerts[dp]! v) →
    (∀ t' : DfsTree, t'.Sub t → dfs.outs t'.v = t'.outs) →
    (∀ k, k ≤ dp → dfs.Anc s.stackVerts[k]! s.stackVerts[dp]! ∧ dfs.depth s.stackVerts[k]! = k) →
    s.RInvTop dfs s.stackVerts[dp]! dp →
    RSideTree dfs t d s

abbrev RSOuts (dfs : DfsData) (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) :
    Prop :=
  ∀ (σ : List Nat) (n : Nat) (F : List RFrame) (B dp : Nat), d = dp + 1 →
    s.Inv' d → Shape s → GuardsOuts v d outs hasVert s → BookOuts v d outs hasVert s →
    CbOuts σ n v d outs hasVert s →
    s.g.TwoConnected → s.g.WF → dfs.Spec s.g → dfs.Rooted s.g → s.stackVerts.size = s.g.nv →
    (∀ o ∈ outs, o ∈ dfs.outs v) →
    (∀ o ∈ outs, ∀ e cls child, o = .tree e cls child → ∀ t' : DfsTree, t'.Sub child → dfs.outs t'.v = t'.outs) →
    AncChain dfs v d s → RWalk dfs F v d s → (∀ f ∈ F, f.2.2 ≤ B) → B ≤ s.tstack.length →
    (hasVert = true → B + 1 ≤ s.tstack.length) →
    RSideOuts dfs v d outs hasVert s

abbrev RSOut (dfs : DfsData) (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  ∀ (σ : List Nat) (n : Nat) (F : List RFrame) (B dp : Nat), d = dp + 1 →
    s.Inv' d → Shape s → GuardsOut v d o hasVert s → BookOut v d o hasVert s →
    CbOut σ n v d o hasVert s →
    s.g.TwoConnected → s.g.WF → dfs.Spec s.g → dfs.Rooted s.g → s.stackVerts.size = s.g.nv →
    o ∈ dfs.outs v →
    (∀ e cls child, o = .tree e cls child → ∀ t' : DfsTree, t'.Sub child → dfs.outs t'.v = t'.outs) →
    AncChain dfs v d s → RWalk dfs F v d s → (∀ f ∈ F, f.2.2 ≤ B) → B ≤ s.tstack.length →
    (hasVert = true → B + 1 ≤ s.tstack.length) →
    RSideOut dfs v d o hasVert s

mutual
theorem rsTree : ∀ (t : DfsTree) (d : Nat) (s : WalkState), RSTree dfs t d s
  | .node v outs, d, s => fun σ n dp hdp hi hs hg hb hcb h2 hwf hsp hrt hsz hp hout hchain hpar => by
    subst hdp
    unfold GuardsTree at hg; unfold BookTree at hb; unfold CbTree at hcb
    unfold RSideTree
    have hp := hp v outs rfl
    have hsize : dp + 1 < s.stackVerts.size := ancChain_lt_size hsp hwf hsz hp hchain
    have hanc := ancChain_child hsp hsize hp hchain
    obtain ⟨hstab, hvs⟩ := rSide_entry_site hi hs h2 hsp hrt hp hchain hpar
    refine ⟨hanc, hstab, by simpa only [Nat.add_sub_cancel] using hvs, ?_⟩
    have hW : RWalk dfs [] v (dp + 1) { s with stackVerts := s.stackVerts.set! (dp + 1) v } :=
      ⟨fun f hf => by simp at hf,
        ⟨fun t ht hd hne => hstab t ht (hpar.entries t ht (by omega) fun hv => by
          have := hvs t ht hv; omega), hpar.disj⟩⟩
    have hvo : dfs.outs v = outs := hout _ (DfsTree.Sub.refl _)
    exact rsOuts v (dp + 1) outs false _ σ n [] s.tstack.length dp rfl (hi.setSv v) hs.frame' hg hb hcb h2 hwf
      hsp hrt (by show (s.stackVerts.set! (dp + 1) v).size = s.g.nv; rw [Array.size_set!]; exact hsz)
      (fun o ho => by rw [hvo]; exact ho)
      (fun o ho e cls child hoe t' ht' => hout t' (DfsTree.Sub.step ht' (hoe ▸ ho)))
      hanc hW (fun f hf => by simp at hf) (le_refl _) (fun h => by cases h)

theorem rsOuts : ∀ (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState),
    RSOuts dfs v d outs hasVert s
  | v, d, [], hasVert, s => fun _ _ F B dp hdp hi hs hg hb _ h2 hwf hsp hrt hsz hmem hout hanc hW hB hBl hnv => by
    unfold BookOuts at hb
    unfold RSideOuts
    intro hv
    subst hv
    exact rSide_vertFree_site hi hs hb h2 hsp hrt hanc hW
  | v, d, o :: rest, hasVert, s => fun σ n F B dp hdp hi hs hg hb hcb h2 hwf hsp hrt hsz hmem hout hanc hW hB hBl hnv => by
    unfold GuardsOuts at hg; unfold BookOuts at hb; unfold CbOuts at hcb
    have hro : RSideOut dfs v d o hasVert s :=
      rsOut v d o hasVert s σ n F B dp hdp hi hs hg.1 hb.1 hcb.1 h2 hwf hsp hrt hsz (hmem o (List.mem_cons_self ..))
        (hout o (List.mem_cons_self ..)) hanc hW hB hBl hnv
    unfold RSideOuts
    refine ⟨hro, ?_⟩
    have hT : Types s.g s := ⟨rfl, hs.size, hs.root, hs.vert, hs.edge⟩
    have hk := kOut v d o hasVert s s.g (d + 1) rootItem s hT (le_refl _) (by show 0 < 1 + _ + _; omega)
      (by show 1 + v ≠ 0; omega) (fun w _ => by show 1 + w ≠ 0; omega)
      (fun e _ => by show 1 + s.g.nv + e ≠ 0; omega) Keep.refl
    refine wp_imp (wp_imp (wp_imp (wp_imp (wp_imp (wp_imp (wp_of_forall fun hv' s' ⟨hi', hs'⟩ hg' hb' hcb'
        ⟨hW', hK', hn', hg'', hsv'⟩ hk' => ?_)
      (invOut v d o hasVert s hi hs hg.1 hb.1)) hg.2) hb.2) hcb.2)
      (rrOut v d o hasVert s F B hi hs hg.1 hb.1 hro h2 hsp hrt hanc hW hB hBl hnv)) hk
    have hanc' : AncChain dfs v d s' :=
      ⟨(hsv' d (le_refl _)).trans hanc.1, fun k hk => by rw [hsv' k hk]; exact hanc.2 k hk⟩
    exact rsOuts v d rest hv' s' σ _ F B dp hdp hi' hs' hg' hb' hcb' (by rw [hg'']; exact h2)
      (by rw [hg'']; exact hwf) (by rw [hg'']; exact hsp) (by rw [hg'']; exact hrt)
      (by rw [hk'.sv, hg'']; exact hsz) (fun o ho => hmem o (List.mem_cons_of_mem _ ho))
      (fun o ho => hout o (List.mem_cons_of_mem _ ho)) hanc' hW' hB hK'.1 hn'

theorem rsOut : ∀ (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState), RSOut dfs v d o hasVert s
  | v, d, o, hasVert, s => fun σ n F B dp hdp hi hs hg hb hcb h2 hwf hsp hrt hsz hmem hout hanc hW hB hBl hnv => by
    have hfr := frOut v d o hasVert s hi hs hg hb
    unfold GuardsOut at hg; unfold BookOut at hb; unfold FrontiersOut at hfr; unfold CbOut at hcb
    have hlow : o.cls.lowval d < d := by
      subst hdp
      obtain ⟨p, hpv⟩ := hanc.parent
      have h := DfsData.Spec.outs_lowval_lt hsp h2 hpv hmem
      have hd : dfs.depth v = dp + 1 := by
        have := (hanc.2 (dp + 1) (le_refl _)).2; rwa [hanc.1] at this
      rwa [hd] at h
    have hfree : hasVert = false → VertFree v s := fun hv =>
      rSide_vertFree_site hi hs (hv ▸ hb.1) h2 hsp hrt hanc hW
    unfold RSideOut
    refine ⟨hlow, hfree, ?_⟩
    refine wp_imp (wp_imp (wp_imp (wp_imp (wp_imp (wp_imp (wp_of_forall fun hv₁ s₁ ⟨hi₁, hs₁⟩ hg₁ hb₁ hfr₁ hcb₁
        ⟨hW₁, hK₁, hn₁, hhv₁, hpush₁, hg₁', hsv₁, hit₁⟩ => ?_)
      (walkOutPre_inv hi hs hb.1)) hg) hb.2) hfr) hcb) (walkOutPre_r hfree hW hBl hnv)
    have hanc₁ : AncChain dfs v d s₁ := by rw [AncChain, hsv₁]; exact hanc
    have h2₁ : s₁.g.TwoConnected := by rw [hg₁']; exact h2
    have hwf₁ : s₁.g.WF := by rw [hg₁']; exact hwf
    have hsp₁ : dfs.Spec s₁.g := by rw [hg₁']; exact hsp
    have hrt₁ : dfs.Rooted s₁.g := by rw [hg₁']; exact hrt
    cases o with
    | back e dest cls => trivial
    | tree e cls child =>
      obtain ⟨c, couts⟩ := child
      try simp only [wp_modify] at hg₁ hb₁ hfr₁ hcb₁
      simp only [wp_modify]
      have hW₂ : RWalk dfs F v d ({ s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne } : WalkState) :=
        RWalk.of_eq (s := s₁) rfl rfl rfl rfl hW₁
      have hpar : ({ s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne } : WalkState).RInvTop dfs
          ({ s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne } : WalkState).stackVerts[d]! d :=
        RInvTop.of_eq (s := s₁) rfl rfl rfl rfl (by
          show s₁.RInvTop dfs s₁.stackVerts[d]! d
          rw [hanc₁.1]; exact hW₁.top)
      have hpv : dfs.IsParent s₁.stackVerts[d]! c := by
        rw [hanc₁.1]; exact ⟨_, hmem, rfl, rfl⟩
      have hside : RSideTree dfs (.node c couts) (d + 1)
          ({ s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne } : WalkState) :=
        rsTree (.node c couts) (d + 1) _ σ n d rfl hi₁.frame' hs₁.frame' hg₁.1 hb₁.1 hcb₁.1 h2₁ hwf₁ hsp₁ hrt₁
          (by show s₁.stackVerts.size = s₁.g.nv; rw [hsv₁, hg₁']; exact hsz)
          (fun w wouts h => by cases h; exact hpv) (hout e cls _ rfl)
          (fun k hk => by
            show dfs.Anc s₁.stackVerts[k]! s₁.stackVerts[d]! ∧ dfs.depth s₁.stackVerts[k]! = k
            rw [hanc₁.1]; exact hanc₁.2 k hk) hpar
      refine ⟨hside, ?_⟩
      refine wp_imp (wp_imp (wp_imp (wp_imp (wp_imp (wp_of_forall fun _ s₃ ⟨hg₃, hb₃⟩ hfr₃ hcb₃ ⟨hi₃, hs₃⟩
          ⟨hW₃, hK₃, hg₃', hsv₃⟩ => ?body) (wp_and hg₁.2 hb₁.2)) hfr₁.2) hcb₁.2)
        (invTree (.node c couts) (d + 1) _ ?pre hs₁.frame' hg₁.1 hb₁.1))
        (rrTree (.node c couts) (d + 1) _ ((v, d, s₁.tstack.length) :: F) d rfl ?pre hs₁.frame' hg₁.1
          hb₁.1 hside h2₁ hsp₁ hrt₁ ?frames hpar)
      case pre => exact fun w outs _ => (hi₁.frame' (s' := { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne })).setSv w
      case frames =>
        intro f hf
        simp only [List.mem_cons] at hf
        rcases hf with rfl | hf
        · exact ⟨le_refl _, hW₂.top.toG _⟩
        · exact hW₂.frames f hf
      case body =>
        have hW₃ := hW₃ c couts rfl
        have hanc₃ : AncChain dfs v d s₃ :=
          ⟨(hsv₃ d (by omega)).trans hanc₁.1, fun k hk => by rw [hsv₃ k (by omega)]; exact hanc₁.2 k hk⟩
        have h2₃ : s₃.g.TwoConnected := by rw [hg₃']; exact h2₁
        have hsp₃ : dfs.Spec s₃.g := by rw [hg₃']; exact hsp₁
        have hrt₃ : dfs.Rooted s₃.g := by rw [hg₃']; exact hrt₁
        exact rSide_finish_site rfl hlow hmem hi₃ hs₃ hg₃ hb₃ hfr₃ hW₃ hanc₃ h2₃ hsp₃ hrt₃
          (closeBase_rContent hcb₃ rfl hlow h2₃ hsp₃ hrt₃)
end

/-- The side facts of the child-return induction at a block's non-root `walkTree` entry
`(c, d + 1)` (`RSideTree`), given the ear bookkeeping `BookTree` of the entry, the parent's
ancestor chain, the layout `stackVerts.size = g.nv`, and the DFS out-lists of the subtree
(`DfsData.ofForest_outs`) and the Ranges site record `CbTree` (for `closeBase_rContent`). -/
theorem walkTree_rSide (s : WalkState) (d c : Nat) (outs : List DfsOut) (σ : List Nat) (n : Nat)
    (hi : s.Inv' d) (hs : Shape s) (hg : GuardsTree (.node c outs) (d + 1) s)
    (hb : BookTree (.node c outs) (d + 1) s) (hcb : CbTree σ n (.node c outs) (d + 1) s)
    (h2 : s.g.TwoConnected) (hwf : s.g.WF) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hsz : s.stackVerts.size = s.g.nv) (hp : dfs.IsParent s.stackVerts[d]! c)
    (hout : ∀ t : DfsTree, t.Sub (.node c outs) → dfs.outs t.v = t.outs)
    (hchain : ∀ k, k ≤ d → dfs.Anc s.stackVerts[k]! s.stackVerts[d]! ∧ dfs.depth s.stackVerts[k]! = k)
    (hR : s.RInvTop dfs s.stackVerts[d]! d) :
    RSideTree dfs (.node c outs) (d + 1) s :=
  rsTree (.node c outs) (d + 1) s σ n d rfl hi hs hg hb hcb h2 hwf hsp hrt hsz
    (fun _ _ h => by cases h; exact hp) hout hchain hR

end WalkState
end Spqr
