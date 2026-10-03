import Spqr.Proofs.RSkelFinish
import Spqr.Proofs.RInvWalk

/-!
# `Items.RSkelInv` through the walk (PROOF.md §4.5)

The R-skeleton invariant along the child-return induction of `RInvWalk`: at every `finishEdge`
site with `lv < d`, `keepsR_finishEdge` with its two R-creating sites discharged — the `.R`
iterate of Loop 1 by `loop1_r_keepsR` (`loop1_rBranch` + `rCloseItems_rSkelInv`, side condition
from `FinishTopOk.side`), the type-1 vertex close by the named admission `closeVert_type1_rSkel3`.
-/

namespace Spqr
open WalkState WalkM

namespace WalkState
variable {s : WalkState} {dfs : DfsData}

theorem RInvFront.congr {s' : WalkState} {v d n : Nat} (hg : s'.g = s.g)
    (hsv : s'.stackVerts = s.stackVerts) (hts : s'.tstack = s.tstack)
    (hty : ∀ t ∈ s.tstack, ∀ i ∈ t.spans.1 ++ t.spans.2, Items.type s'.items i = Items.type s.items i)
    (hvs : ∀ t ∈ s.tstack, ∀ i ∈ t.spans.1 ++ t.spans.2, Items.vs s'.items i = Items.vs s.items i)
    (hE : ∀ t ∈ s.tstack, ∀ i ∈ t.spans.1 ++ t.spans.2, ∀ e,
      Items.EdgeBelow s.g s'.items i e ↔ Items.EdgeBelow s.g s.items i e)
    (h : s.RInvFront dfs v d n) : s'.RInvFront dfs v d n := by
  refine ⟨fun t ht hd hv => ?_, ?_⟩
  · rw [hts] at ht
    exact EntryR.congr (s := s) hg hsv (hty t (List.mem_of_mem_drop ht))
      (hvs t (List.mem_of_mem_drop ht)) (hE t (List.mem_of_mem_drop ht)) (h.base t ht hd hv)
  · rw [hts, hg]
    exact List.Pairwise.imp_of_mem (fun {a b} ha hb hab e he hbe =>
      hab e ((TEntry.edges_congr (hE a ha) e).1 he) ((TEntry.edges_congr (hE b hb) e).1 hbe)) h.disj

theorem RInvFront.modifyVs_free {v d n : Nat} (j : ItemId) (f : Item → Item)
    (hch : ∀ it, (f it).ch = it.ch) (hty : ∀ it, (f it).type = it.type)
    (hfree : ∀ t ∈ s.tstack, j ∉ t.spans.1 ++ t.spans.2)
    (h : s.RInvFront dfs v d n) : (after (modifyItem j f) s).RInvFront dfs v d n := by
  refine RInvFront.congr (s := s) (s' := after (modifyItem j f) s) rfl rfl rfl
    (fun t ht i hi => ?_) (fun t ht i hi => ?_) (fun t ht i hi e => ?_) h
  · exact Items.type_modify_type_eq j f hty i
  · exact Items.vs_modify_of_ne j f fun hij => hfree t ht (hij ▸ hi)
  · exact Items.Below_modify_ch_eq j f hch

/-- The `.R` iterate of Loop 1 keeps `Items.RSkelInv`: its fresh R item is
`rCloseItems_rSkelInv`'s (shape and `RTop` from `loop1_rBranch`), with the empty-side condition
read off `FinishTopOk.side` of the close. -/
theorem loop1_r_keepsR {v nxtV d e : Nat} {edgeDir : Bool} (hi : s.Inv' (d + 1)) (hs : Shape s)
    (hok : CloseEarsOk (d + 1) nxtV d e edgeDir s) (hv : v < s.g.nv)
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hc : nxtV = s.stackVerts[d + 1]!) {origTstack : Nat}
    (hR : s.RInvFront dfs s.stackVerts[d]! d origTstack) (k : Nat)
    (hk : ∀ j, j ≤ k → result (loop1Cond d) (iter (loop1Body d edgeDir) j (ceS₁ nxtV d e s)) = true)
    (hty : l1Ty d edgeDir (iter (loop1Body d edgeDir) k (ceS₁ nxtV d e s)) = .R) :
    KeepsR (loop1Body d edgeDir) (iter (loop1Body d edgeDir) k (ceS₁ nxtV d e s)) := by
  obtain ⟨cur, nxt, rest, hb, hR'⟩ := loop1_rBranch hi hs hok h2 hsp hrt hc hR k hk hty
  have st := closeEars_iter_step (v := v) hi hs hv hok k hk
  set sk := iter (loop1Body d edgeDir) k (ceS₁ nxtV d e s) with hsk
  have hg : sk.g = s.g := st.g
  have hside : getSide (TEntry.mergeInto cur nxt).spans (!sk.stackDir[d]!) = [] := by
    have h := (hok.body k hk).close.finish
    have hl2 : l1S₂ d edgeDir sk = { sk with items := sk.items.push ⟨.R, (none, none), []⟩ } := by
      unfold l1S₂; rw [hty, l1S₁_of_R hty]
      unfold after
      rw [maybeUnwrapNxt_run_eq .R sk cur nxt rest hb.tstack _ rfl _ rfl]
      simp only [true_or, ↓reduceIte, run_allocItem]
    rw [hl2] at h
    unfold after at h
    rw [mergeTstackTops_run_eq { sk with items := sk.items.push ⟨.R, (none, none), []⟩ } cur nxt rest
      hb.tstack] at h
    have h := h.side
    simp only [curE, List.head!_cons] at h
    have htop : (TEntry.mergeInto cur nxt).topDepth = d := by
      simp [TEntry.mergeInto, hb.cur_top, hb.nxt_top]
    rw [htop] at h
    exact h
  unfold KeepsR
  intro hinv
  rw [loop1Body_run_of_R hb.tstack hb.cur_top hb.nxt_top hty]
  exact rCloseItems_rSkelInv hinv st.shape hb hR' st.inv (hg ▸ h2) (hg ▸ hsp) (hg ▸ hrt) _ hside

/-- Admitted (PROOF.md §4.5): the type-1 vertex close of a returning tree edge with
`isSingle = false` (`maybeUnwrapNxt .R` at `cvS₁`, then two merges, the retarget to `curV` and
`finishTstackTop`) builds an R item whose skeleton is 3-connected: the pieces are the items of
the merged ear (the whole subtree ear at `curV` plus the back edge and the vertex entry), the
parent piece its complement at `(curV, stackVerts[lv])`. The vertex-level analogue of
`RBranch.rSkel3`: it needs the merged entry's items to be pairwise edge-disjoint maximal pieces
(`RInvH` at `feS₂`, `FinishRShape.unwrap`) and the HT argument that no separation pair survives
(`RCloseShape.threeConnected`). -/
theorem closeVert_type1_rSkel3 {D : Nat} (curV d lv : Nat) (kind : RetKind) (o : DfsOut)
    (origTstack : Nat) (ho : o.cls = .ret lv kind) (hk : kind ≠ .backEdge) (hlow : lv < d)
    (hv : curV < s.g.nv) (hi : s.Inv' D) (hs : Shape s)
    (hok : FinishOk D curV d lv o origTstack true s)
    (hfront : Frontier (o := o) d origTstack s)
    (hshape : FinishRShape dfs curV d o origTstack true s)
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hd : dfs.depth curV = d) (hcur : s.stackVerts[d]! = curV)
    (hanc : ∀ k, k ≤ d → dfs.Anc s.stackVerts[k]! curV ∧ dfs.depth s.stackVerts[k]! = k)
    (hR : s.RInvFront dfs curV d origTstack)
    (h1 : o.cls.isType1 = true) (hsingle : feSingle d o s = false) :
    Items.RSkel3 s.g (feS₃ curV d o origTstack s).items
      (cvS₁ o.cls.isType1 origTstack (feSingle d o s) (feS₂ d o s)).items.size := by
  sorry

/-- `KeepsR` at a returning `finishEdge` site of the walk induction (`lv < d`), from the site
facts of `rrOut`: the Loop-1 R close by `loop1_r_keepsR`, the type-1 vertex close by
`closeVert_type1_rSkel3`. -/
theorem keepsR_finishEdge_site {D : Nat} (curV d lv : Nat) (kind : RetKind) (o : DfsOut)
    (origTstack : Nat) (hasVert : Bool) (ho : o.cls = .ret lv kind) (hlow : lv < d)
    (hv : curV < s.g.nv) (hi : s.Inv' D) (hs : Shape s)
    (hok : FinishOk D curV d lv o origTstack hasVert s)
    (hD : D = if o.cls.isTree then d + 1 else d)
    (hfront : Frontier (o := o) d origTstack s)
    (hshape : kind ≠ .backEdge → FinishRShape dfs curV d o origTstack hasVert s)
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hd : dfs.depth curV = d) (hcur : s.stackVerts[d]! = curV)
    (hanc : ∀ k, k ≤ d → dfs.Anc s.stackVerts[k]! curV ∧ dfs.depth s.stackVerts[k]! = k)
    (hR : s.RInvFront dfs curV d origTstack)
    (hq : ∀ p, ¬ Items.IsParent s.items p (edgeItem s.g o.e))
    (hchild : o.cls.isTree = true → s.stackVerts[d + 1]! = o.dest) :
    KeepsR (finishEdge curV d o origTstack hasVert) s := by
  have hk : o.cls.isTree = true → kind ≠ .backEdge := fun ht h => by
    subst h; rw [ho] at ht; cases ht
  refine keepsR_finishEdge curV d lv kind o origTstack hasVert ho hlow hv hi hs hok hq ?_ ?_
  · intro ht k hk' hty
    have hshape := hshape (hk ht)
    obtain rfl : D = d + 1 := by simp [hD, ht]
    have st₀ : Step (d + 1) curV s (feS₀ d o s) :=
      Step.modifyVs hi hs (edgeItem s.g o.e) _ (by show 1 + s.g.nv + o.e < _; have := hok.e_lt; omega)
    have hfree₀ : ∀ t ∈ s.tstack, edgeItem s.g o.e ∉ t.spans.1 ++ t.spans.2 := fun t ht hmem =>
      hshape.pend t ht ⟨_, hmem, .refl⟩
    have R₀ : (feS₀ d o s).RInvFront dfs (feS₀ d o s).stackVerts[d]! d origTstack := by
      rw [st₀.sv, hcur]
      exact hR.modifyVs_free (edgeItem s.g o.e) _ (fun _ => rfl) (fun _ => rfl) hfree₀
    exact loop1_r_keepsR st₀.inv st₀.shape (hok.ears ht) (by rw [st₀.g]; exact hv)
      (by rw [st₀.g]; exact h2) (by rw [st₀.g]; exact hsp) (by rw [st₀.g]; exact hrt)
      (by rw [st₀.sv]; exact (hchild ht).symm) R₀ k hk' hty
  · intro ht hhv h1 hsingle
    subst hhv
    exact closeVert_type1_rSkel3 curV d lv kind o origTstack ho (hk ht) hlow hv hi hs hok hfront
      (hshape (hk ht)) h2 hsp hrt hd hcur hanc hR h1 hsingle

/-! ## The walk induction

The shape of `rrTree`/`rrOuts`/`rrOut` (`RInvWalk`), whose outputs (`RWalk`, `BotKeep`, graph and
`stackVerts` frames) supply the site facts of `keepsR_finishEdge_site`. -/

abbrev RKTree (dfs : DfsData) (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  ∀ (F : List RFrame) (dp : Nat), d = dp + 1 →
    (∀ v outs, t = .node v outs → ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).Inv' d) →
    Shape s → GuardsTree t d s → BookTree t d s → RSideTree dfs t d s →
    s.g.TwoConnected → dfs.Spec s.g → dfs.Rooted s.g →
    (∀ f ∈ F, f.2.2 ≤ s.tstack.length ∧ s.RInvG dfs f.1 f.2.1 f.2.2) →
    s.RInvTop dfs s.stackVerts[dp]! dp → Items.RSkelInv s.g s.items →
    wp (walkTree t d) (fun _ s' => Items.RSkelInv s.g s'.items) s

abbrev RKOuts (dfs : DfsData) (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) :
    Prop :=
  ∀ (F : List RFrame) (B : Nat),
    s.Inv' d → Shape s → GuardsOuts v d outs hasVert s → BookOuts v d outs hasVert s →
    RSideOuts dfs v d outs hasVert s → s.g.TwoConnected → dfs.Spec s.g → dfs.Rooted s.g →
    AncChain dfs v d s → RWalk dfs F v d s → (∀ f ∈ F, f.2.2 ≤ B) → B ≤ s.tstack.length →
    (hasVert = true → B + 1 ≤ s.tstack.length) → Items.RSkelInv s.g s.items →
    wp (walkOuts v d outs hasVert) (fun _ s' => Items.RSkelInv s.g s'.items) s

abbrev RKOut (dfs : DfsData) (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  ∀ (F : List RFrame) (B : Nat),
    s.Inv' d → Shape s → GuardsOut v d o hasVert s → BookOut v d o hasVert s →
    RSideOut dfs v d o hasVert s → s.g.TwoConnected → dfs.Spec s.g → dfs.Rooted s.g →
    AncChain dfs v d s → RWalk dfs F v d s → (∀ f ∈ F, f.2.2 ≤ B) → B ≤ s.tstack.length →
    (hasVert = true → B + 1 ≤ s.tstack.length) → Items.RSkelInv s.g s.items →
    wp (walkOut v d o hasVert) (fun _ s' => Items.RSkelInv s.g s'.items) s

mutual
theorem rkTree : ∀ (t : DfsTree) (d : Nat) (s : WalkState), RKTree dfs t d s
  | .node v outs, d, s => fun F dp hdp hi hs hg hb hr h2 hsp hrt hF hpar hinv => by
    subst hdp
    unfold walkTree
    simp only [wp_bind, wp_modify]
    unfold GuardsTree at hg; unfold BookTree at hb; unfold RSideTree at hr
    obtain ⟨hanc, hstab, hvs, hr⟩ := hr
    simp only [Nat.add_sub_cancel] at hvs
    have hW₀ : RWalk dfs F v (dp + 1) { s with stackVerts := s.stackVerts.set! (dp + 1) v } :=
      ⟨fun f hf => ⟨(hF f hf).1, ⟨fun t ht hd hne => hstab t (List.mem_of_mem_drop ht) ((hF f hf).2.entries t ht hd hne),
          (hF f hf).2.disj⟩⟩,
        ⟨fun t ht hd hne => hstab t ht (hpar.entries t ht (by omega) fun h => by
          have := hvs t ht h; omega), hpar.disj⟩⟩
    refine wp_imp (wp_of_forall fun hv s₁ hK₁ => ?_)
      (rkOuts v (dp + 1) outs false _ F s.tstack.length (hi v outs rfl) hs.frame' hg hb hr h2 hsp hrt
        hanc hW₀ (fun f hf => (hF f hf).1) (le_refl _) (fun h => nomatch h) hinv)
    cases hv
    · simp only [Bool.false_eq_true, ↓reduceIte, wp_bind, wp_setStackDir, wp_pushVertTstack]
      exact hK₁
    · exact hK₁

theorem rkOuts : ∀ (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState),
    RKOuts dfs v d outs hasVert s
  | v, d, [], hasVert, s => fun F B hi hs hg hb hr h2 hsp hrt hanc hW hB hBl hnv hinv => by
    unfold walkOuts
    rw [wp_pure]
    exact hinv
  | v, d, o :: rest, hasVert, s => fun F B hi hs hg hb hr h2 hsp hrt hanc hW hB hBl hnv hinv => by
    unfold GuardsOuts at hg; unfold BookOuts at hb; unfold RSideOuts at hr
    unfold walkOuts
    simp only [wp_bind]
    refine wp_imp (wp_imp (wp_imp (wp_imp (wp_imp (wp_imp (wp_of_forall fun hv' s' ⟨hi', hs'⟩ hg' hb' hr'
        ⟨hW', hK', hn', hg'', hsv'⟩ hinv' => ?_)
      (invOut v d o hasVert s hi hs hg.1 hb.1)) hg.2) hb.2) hr.2)
      (rrOut v d o hasVert s F B hi hs hg.1 hb.1 hr.1 h2 hsp hrt hanc hW hB hBl hnv))
      (rkOut v d o hasVert s F B hi hs hg.1 hb.1 hr.1 h2 hsp hrt hanc hW hB hBl hnv hinv)
    have hanc' : AncChain dfs v d s' :=
      ⟨(hsv' d (le_refl _)).trans hanc.1, fun k hk => by rw [hsv' k hk]; exact hanc.2 k hk⟩
    refine wp_imp (wp_of_forall fun hv'' s'' h => ?_)
      (rkOuts v d rest hv' s' F B hi' hs' hg' hb' hr' (by rw [hg'']; exact h2) (by rw [hg'']; exact hsp)
        (by rw [hg'']; exact hrt) hanc' hW' hB hK'.1 hn' (by rw [hg'']; exact hinv'))
    rw [hg''] at h; exact h

theorem rkOut : ∀ (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState), RKOut dfs v d o hasVert s
  | v, d, o, hasVert, s => fun F B hi hs hg hb hr h2 hsp hrt hanc hW hB hBl hnv hinv => by
    rw [walkOut_eq, wp_bind]
    unfold GuardsOut at hg; unfold BookOut at hb; unfold RSideOut at hr
    obtain ⟨hlow, hfree, hr⟩ := hr
    refine wp_imp (wp_imp (wp_imp (wp_imp (wp_imp (wp_of_forall fun hv₁ s₁ ⟨hi₁, hs₁⟩ hg₁ hb₁ hr₁
        ⟨hW₁, hK₁, hn₁, hhv₁, hpush₁, hg₁', hsv₁, hit₁⟩ => ?_)
      (walkOutPre_inv hi hs hb.1)) hg) hb.2) hr) (walkOutPre_r hfree hW hBl hnv)
    have hinv₁ : Items.RSkelInv s₁.g s₁.items := by rw [hg₁', hit₁]; exact hinv
    have hanc₁ : AncChain dfs v d s₁ := by rw [AncChain, hsv₁]; exact hanc
    have h2₁ : s₁.g.TwoConnected := by rw [hg₁']; exact h2
    have hsp₁ : dfs.Spec s₁.g := by rw [hg₁']; exact hsp
    have hrt₁ : dfs.Rooted s₁.g := by rw [hg₁']; exact hrt
    obtain ⟨lv, kind, ho, hl⟩ := ret_of_lowval_lt hlow
    unfold walkOutRest
    rw [wp_bind, wp_tstackSize]
    cases o with
    | back e cls dest =>
      try simp only at hg₁ hb₁
      simp only
      have hk : kind = .backEdge := by
        cases kind with
        | backEdge => rfl
        | type1Child => exact absurd (hb₁.tree.1 (by rw [ho]; rfl)) (by simp)
        | type2Child => exact absurd (hb₁.tree.1 (by rw [ho]; rfl)) (by simp)
      subst hk
      have hhv : hv₁ = true := hpush₁ hlow (by rw [ho]; rfl)
      subst hhv
      have hnt : (DfsOut.back e cls dest).cls.isTree = false := by
        cases h : (DfsOut.back e cls dest).cls.isTree
        · rfl
        · obtain ⟨_, _, _, h'⟩ := hb₁.tree.1 h; cases h'
      have hD : d = if (DfsOut.back e cls dest).cls.isTree then d + 1 else d := by simp [hnt]
      obtain ⟨sub, base, hlen, hE⟩ := hb₁.ear
      have hok := finishOk_of_guards ho hl hg₁ hE hlen hi₁ hs₁ hD hb₁.v_lt hb₁.e_lt hb₁.q
        (hb₁.ends lv .backEdge ho) hb₁.vert
      have hfront := finishEdge_frontier hb₁ hi₁ hs₁ hD
      have hhd : dfs.depth v = d := by
        have := (hanc₁.2 d (le_refl _)).2; rwa [hanc₁.1] at this
      have hR : s₁.RInvFront dfs v d s₁.tstack.length := hW₁.top.toFront _
      have hkeep := keepsR_finishEdge_site v d lv .backEdge _ s₁.tstack.length true ho hl hb₁.v_lt
        hi₁ hs₁ hok hD hfront (fun h => absurd rfl h) h2₁ hsp₁ hrt₁ hhd hanc₁.1 hanc₁.2 hR hE.q_root
        hE.sv_child
      have := hkeep hinv₁
      rw [hg₁'] at this
      exact this
    | tree e cls child =>
      obtain ⟨c, couts⟩ := child
      try simp only [wp_modify] at hg₁ hb₁ hr₁
      simp only [wp_bind, wp_modify]
      set s₂ : WalkState := { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne } with hs₂
      have hW₂ : RWalk dfs F v d s₂ := RWalk.of_eq (s := s₁) rfl rfl rfl rfl hW₁
      have hpar : s₂.RInvTop dfs s₂.stackVerts[d]! d :=
        RInvTop.of_eq (s := s₁) rfl rfl rfl rfl (by
          show s₁.RInvTop dfs s₁.stackVerts[d]! d
          rw [hanc₁.1]; exact hW₁.top)
      have hpre := fun w outs (_ : DfsTree.node c couts = .node w outs) =>
        (hi₁.frame' (s' := s₂)).setSv w
      have hframes : ∀ f ∈ (v, d, s₁.tstack.length) :: F,
          f.2.2 ≤ s₂.tstack.length ∧ s₂.RInvG dfs f.1 f.2.1 f.2.2 := by
        intro f hf
        simp only [List.mem_cons] at hf
        rcases hf with rfl | hf
        · exact ⟨le_refl _, hW₂.top.toG _⟩
        · exact hW₂.frames f hf
      refine wp_imp (wp_imp (wp_imp (wp_imp (wp_imp (wp_of_forall fun _ s₃ ⟨hg₃, hb₃⟩ hr₃ ⟨hi₃, hs₃⟩
          ⟨hW₃, hK₃, hg₃', hsv₃⟩ hinv₃ => ?body) (wp_and hg₁.2 hb₁.2)) hr₁.2)
        (invTree (.node c couts) (d + 1) s₂ hpre hs₁.frame' hg₁.1 hb₁.1))
        (rrTree (.node c couts) (d + 1) s₂ ((v, d, s₁.tstack.length) :: F) d rfl hpre hs₁.frame' hg₁.1
          hb₁.1 hr₁.1 h2₁ hsp₁ hrt₁ hframes hpar))
        (rkTree (.node c couts) (d + 1) s₂ ((v, d, s₁.tstack.length) :: F) d rfl hpre hs₁.frame' hg₁.1
          hb₁.1 hr₁.1 h2₁ hsp₁ hrt₁ hframes hpar hinv₁)
      case body =>
        have hW₃ := hW₃ c couts rfl
        have hanc₃ : AncChain dfs v d s₃ :=
          ⟨(hsv₃ d (by omega)).trans hanc₁.1, fun k hk => by rw [hsv₃ k (by omega)]; exact hanc₁.2 k hk⟩
        have h2₃ : s₃.g.TwoConnected := by rw [hg₃']; exact h2₁
        have hsp₃ : dfs.Spec s₃.g := by rw [hg₃']; exact hsp₁
        have hrt₃ : dfs.Rooted s₃.g := by rw [hg₃']; exact hrt₁
        have hD : d + 1 = if (DfsOut.tree e cls (.node c couts)).cls.isTree then d + 1 else d := by
          simp [hb₃.tree.2 ⟨_, _, _, rfl⟩]
        obtain ⟨sub, base, hlen, hE⟩ := hb₃.ear
        have hok := finishOk_of_guards ho hl hg₃ hE hlen hi₃ hs₃ hD hb₃.v_lt hb₃.e_lt hb₃.q
          (hb₃.ends lv kind ho) hb₃.vert
        have hfront := finishEdge_frontier hb₃ hi₃ hs₃ hD
        have hhd : dfs.depth v = d := by
          have := (hanc₃.2 d (le_refl _)).2; rwa [hanc₃.1] at this
        have hfr := (hW₃.frames (v, d, s₁.tstack.length) (List.mem_cons_self ..)).2
        have hR : s₃.RInvFront dfs v d s₁.tstack.length := ⟨hfr.entries, hfr.disj⟩
        have hinv₃' : Items.RSkelInv s₃.g s₃.items := by rw [hg₃']; exact hinv₃
        have hkeep := keepsR_finishEdge_site v d lv kind _ s₁.tstack.length hv₁ ho hl hb₃.v_lt hi₃ hs₃
          hok hD hfront (fun _ => hr₃) h2₃ hsp₃ hrt₃ hhd hanc₃.1 hanc₃.2 hR hE.q_root hE.sv_child
        have := hkeep hinv₃'
        rw [hg₃'.trans hg₁'] at this
        exact this
end

end WalkState
end Spqr
