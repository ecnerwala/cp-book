import Spqr.Proofs.RSkelFinish
import Spqr.Proofs.RInvWalk
import Spqr.Proofs.RLoop1

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

/-- `loop1_r_keepsR` under the ear context of the site (`L1Ctx`/`L1Inv`), via `loop1_rBranch_ctx`. -/
theorem loop1_r_keepsR' {D v d : Nat} {o : DfsOut} {hi lo base : List TEntry}
    (hc : L1Ctx D d o s hi lo base)
    (h0 : L1Inv D d o s hi lo base (ceS₁ o.dest d o.e (feS₀ d o s)))
    (hD : D = d + 1) (he : o.e < s.g.ne)
    (hi₀ : (feS₀ d o s).Inv' (d + 1)) (hs₀ : Shape (feS₀ d o s))
    (hok : CloseEarsOk (d + 1) o.dest d o.e s.stackDir[d]! (feS₀ d o s)) (hv : v < s.g.nv)
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hchild : o.dest = s.stackVerts[d + 1]!) {origTstack : Nat}
    (hR : (feS₀ d o s).RInvFront dfs s.stackVerts[d]! d origTstack) (k : Nat)
    (hk : ∀ j, j ≤ k → result (loop1Cond d) (rl1Iter d o s j) = true)
    (hty : l1Ty d s.stackDir[d]! (rl1Iter d o s k) = .R) :
    KeepsR (loop1Body d s.stackDir[d]!) (rl1Iter d o s k) := by
  subst hD
  obtain ⟨cur, nxt, rest, hb, hR'⟩ := loop1_rBranch_ctx hc h0 hv he hi₀ hs₀ hok h2 hsp hrt hchild hR k hk hty
  have st := closeEars_iter_step (v := v) hi₀ hs₀ hv hok k hk
  set sk := rl1Iter d o s k with hsk
  have hg : sk.g = s.g := st.g
  have hside : getSide (TEntry.mergeInto cur nxt).spans (!sk.stackDir[d]!) = [] := by
    have h := (hok.body k hk).close.finish
    have hl2 : l1S₂ d s.stackDir[d]! sk = { sk with items := sk.items.push ⟨.R, (none, none), []⟩ } := by
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
    obtain ⟨sub, base, -, hE⟩ := hshape.ear
    have hlow' : o.cls.lowval d < d := by
      have : o.cls.lowval d = lv := by rw [ho]; rfl
      omega
    obtain ⟨hi', lo, hrange, hctx⟩ := L1Ctx.ofEar hE hs rfl ht hlow'
    have hq' : Items.ch s.items (edgeItem s.g o.e) = [] := by
      have h := (hok.ears ht).q
      exact (Items.ch_modify_ch_eq (edgeItem s.g o.e) (fun it =>
        { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) })
        (fun _ => rfl) _).symm.trans h
    have h0 := l1_init hE hi hs rfl hok.e_lt hq' (hok.ears ht).ends hrange
    exact loop1_r_keepsR' hctx h0 rfl hok.e_lt st₀.inv st₀.shape (hok.ears ht) hv h2 hsp hrt
      (hchild ht).symm (by rw [st₀.sv] at R₀; exact R₀) k hk' hty
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

/-! ## The boundary branch (`finishBoundary`): no R item is created or re-parented -/

theorem Items.noParent_modify {items : Items} {j k : ItemId} (f : Item → Item)
    (hk : ∀ p, ¬ Items.IsParent items p k)
    (hnew : ∀ hj : j < items.size, k ∉ (f items[j]).ch) :
    ∀ p, ¬ Items.IsParent (items.modify j f) p k :=
  fun p hp => (Items.IsParent_modify hp).elim (hk p) fun ⟨_, hj, hmem⟩ => hnew hj hmem

theorem Items.noParent_push_size {items : Items} (x : Item) (hx : x.ch = [])
    (hc : ∀ p c, Items.IsParent items p c → c < items.size) :
    ∀ p, ¬ Items.IsParent (items.push x) p items.size :=
  fun p hp => Nat.lt_irrefl _ (hc p _ (Items.IsParent_push hx hp))

theorem Items.chLt_modify {items : Items} (j : ItemId) (f : Item → Item) (hf : ∀ it, (f it).ch = it.ch)
    (hc : ∀ p c, Items.IsParent items p c → c < items.size) :
    ∀ p c, Items.IsParent (items.modify j f) p c → c < (items.modify j f).size := by
  intro p c hp
  rw [Array.size_modify]
  unfold Items.IsParent at hp
  rw [Items.ch_modify_ch_eq j f hf] at hp
  exact hc p c hp

theorem Items.type_modify_ch {items : Items} (j : ItemId) (L : Item → List ItemId) (p : ItemId) :
    Items.type (items.modify j fun it => { it with ch := L it }) p = items.type p :=
  Items.type_modify_type_eq j (fun it => { it with ch := L it }) (fun _ => rfl) p

theorem Items.type_modify_vs {items : Items} (j : ItemId) (v : Option Nat × Option Nat) (p : ItemId) :
    Items.type (items.modify j fun it => { it with vs := v }) p = items.type p :=
  Items.type_modify_type_eq j (fun it => { it with vs := v }) (fun _ => rfl) p

theorem keepsR_finishBoundary {D : Nat} (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool)
    (hge : d ≤ o.cls.lowval d) (hs : Shape s) (hv : curV < s.g.nv) (he : o.e < s.g.ne)
    (hok : BoundaryOk D curV d o s) :
    KeepsR (finishEdge curV d o origTstack hasVert) s := by
  intro h0
  have hge' : o.cls.lowval d ≥ d := hge
  have hqt : Items.type s.items (edgeItem s.g o.e) = .Q := hs.edge _ he
  have hvt : Items.type s.items (vertItem curV) = .V := hs.vert _ hv
  have hqlt : edgeItem s.g o.e < s.items.size := by
    have := hs.size; show 1 + s.g.nv + o.e < _; omega
  have hvlt : vertItem curV < s.items.size := by
    have := hs.size; show 1 + curV < _; omega
  have hc0 := hs.ch_lt
  show wp (finishEdge curV d o origTstack hasVert) (fun _ s' => Items.RSkelInv s.g s'.items) s
  rw [finishEdge_eq]
  simp only [finishEdge', wp_bind, wp_get, wp_stackDir, hge', ↓reduceIte]
  unfold finishBoundary
  simp only [wp_bind, wp_modifyItem, wp_modify, wp_ite, wp_allocItem, wp_makeVs, wp_popTstack,
    wp_pure]
  obtain ⟨I₁, hI₁⟩ : ∃ I₁ : Items,
      I₁ = s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) } :=
    ⟨_, rfl⟩
  rw [← hI₁]
  have h1 : Items.RSkelInv s.g I₁ := by
    rw [hI₁]; exact h0.modify_root_of_ne hok.q_root _ (fun _ => rfl) (by rw [hqt]; decide)
  have hq1 : ∀ p, ¬ Items.IsParent I₁ p (edgeItem s.g o.e) := by
    rw [hI₁]; exact fun p hp => hok.q_root p ((isParent_modifyVs_iff ..).1 hp)
  have hv1 : ∀ p, ¬ Items.IsParent I₁ p (vertItem curV) := by
    rw [hI₁]; exact fun p hp => hok.v_root p ((isParent_modifyVs_iff ..).1 hp)
  have hc1 : ∀ p c, Items.IsParent I₁ p c → c < I₁.size := by
    rw [hI₁]; exact Items.chLt_modify (edgeItem s.g o.e) _ (fun _ => rfl) hc0
  have hsz1 : I₁.size = s.items.size := by rw [hI₁]; exact Array.size_modify ..
  have hqt1 : Items.type I₁ (edgeItem s.g o.e) = .Q := by rw [hI₁, Items.type_modify_vs, hqt]
  have hvt1 : Items.type I₁ (vertItem curV) = .V := by rw [hI₁, Items.type_modify_vs, hvt]
  clear hI₁
  have hqne : edgeItem s.g o.e ≠ I₁.size := Nat.ne_of_lt (by rw [hsz1]; exact hqlt)
  have hvne : vertItem curV ≠ I₁.size := Nat.ne_of_lt (by rw [hsz1]; exact hvlt)
  split
  · rename_i hT
    have hpops := hok.pops hT
    split
    · rename_i hL
      simp only [hL, ↓reduceIte] at hpops
      obtain ⟨t, rest, hts⟩ : ∃ t rest, s.tstack = t :: rest := by
        match h : s.tstack, hpops with
        | t :: rest, _ => exact ⟨t, rest, rfl⟩
      have ht : t ∈ s.tstack := by rw [hts]; exact List.mem_cons_self ..
      simp only [hts, List.head!_cons]
      refine Items.RSkelInv.modify_root_of_ne ?h4 ?hv4 _ (fun _ => rfl) ?hvt4
      case h4 =>
        refine Items.RSkelInv.modify_root_of_ne ?h3 ?hq3 _ (fun _ => rfl) ?hqt3
        case h3 =>
          refine Items.RSkelInv.modify_root_of_ne ?h2 ?hn2 _ (fun _ => rfl) ?hnt2
          case h2 => exact h1.push_nil (fun i _ c hc => hc1 i c hc) _ rfl (by decide)
          case hn2 => exact Items.noParent_push_size _ rfl hc1
          case hnt2 => rw [Items.type_push_size]; decide
        case hq3 =>
          intro p hp
          exact hq1 p (Items.IsParent_push rfl ((isParent_modifyVs_iff ..).1 hp))
        case hqt3 => rw [Items.type_modify_vs, Items.type_push_of_ne _ hqne, hqt1]; decide
      case hv4 =>
        refine Items.noParent_modify _ ?_ fun _ hmem => ?_
        · intro p hp
          exact hv1 p (Items.IsParent_push rfl ((isParent_modifyVs_iff ..).1 hp))
        · simp only [List.mem_cons] at hmem
          rcases hmem with h | h
          · exact hvne h
          · exact hok.v_free t ht (List.mem_append_right _ h)
      case hvt4 =>
        rw [Items.type_modify_ch, Items.type_modify_vs, Items.type_push_of_ne _ hvne, hvt1]; decide
    · rename_i hL
      simp only [hL] at hpops
      obtain ⟨t₁, t₂, rest, hts⟩ : ∃ t₁ t₂ rest, s.tstack = t₁ :: t₂ :: rest := by
        rcases hl : s.tstack with _ | ⟨t₁, _ | ⟨t₂, rest⟩⟩
        · rw [hl] at hpops; simp at hpops
        · rw [hl] at hpops; simp at hpops
        · exact ⟨t₁, t₂, rest, rfl⟩
      simp only [hts, List.head!_cons, List.tail_cons]
      refine Items.RSkelInv.modify_root_of_ne ?h3 ?hv3 _ (fun _ => rfl) ?hvt3
      case h3 => exact h1.modify_root_of_ne hq1 _ (fun _ => rfl) (by rw [hqt1]; decide)
      case hv3 =>
        refine Items.noParent_modify _ hv1 fun _ hmem => ?_
        simp only [List.mem_append] at hmem
        rcases hmem with h | h
        · exact hok.v_free t₁ (by rw [hts]; exact List.mem_cons_self ..) (List.mem_append_left _ h)
        · exact hok.v_free t₂ (by rw [hts]; exact List.mem_cons_of_mem _ (List.mem_cons_self ..))
            (List.mem_append_right _ h)
      case hvt3 => rw [Items.type_modify_ch, hvt1]; decide
  · rename_i hT
    refine Items.RSkelInv.modify_root_of_ne ?h4 ?hv4 _ (fun _ => rfl) ?hvt4
    case h4 =>
      refine Items.RSkelInv.modify_root_of_ne ?h3 ?hq3 _ (fun _ => rfl) ?hqt3
      case h3 =>
        refine Items.RSkelInv.modify_root_of_ne ?h2 ?hn2 _ (fun _ => rfl) ?hnt2
        case h2 => exact h1.push_nil (fun i _ c hc => hc1 i c hc) _ rfl (by decide)
        case hn2 => exact Items.noParent_push_size _ rfl hc1
        case hnt2 => rw [Items.type_push_size]; decide
      case hq3 =>
        intro p hp
        exact hq1 p (Items.IsParent_push rfl ((isParent_modifyVs_iff ..).1 hp))
      case hqt3 => rw [Items.type_modify_vs, Items.type_push_of_ne _ hqne, hqt1]; decide
    case hv4 =>
      refine Items.noParent_modify _ ?_ fun _ hmem => ?_
      · intro p hp
        exact hv1 p (Items.IsParent_push rfl ((isParent_modifyVs_iff ..).1 hp))
      · simp only [List.mem_singleton] at hmem
        exact hvne hmem
    case hvt4 =>
      rw [Items.type_modify_ch, Items.type_modify_vs, Items.type_push_of_ne _ hvne, hvt1]; decide

end WalkState
end Spqr
