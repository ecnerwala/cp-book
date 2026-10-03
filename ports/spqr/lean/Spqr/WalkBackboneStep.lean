import Spqr.WalkBackbone

/-!
# Backbone: the `walkOut` step from the `finishEdge` primitives

`bbOut_step` opens `walkOut` (`walkOutPre` site → child entry → child's `WalkInvEnd` → `finishEdge`)
and assembles the post-`finishEdge` `WalkInvOut` from the per-primitive preservation lemmas
(`finishEdge_step`/`finishEdge_ranges`/`finishEdge_closeInv`/`finishEdge_ownedD`/`finishEdge_vertCover`/
`finishEdge_full`/`keep_finishEdge`, the R `finishEdge_rInvG_base`/`finishEdge_rInvTop`/
`keepsR_finishEdge_site`/`keepsR_finishBoundary`, the ST `finishBoundary_st`/`stRet_finish`).
-/

namespace Spqr
open WalkM StRefEt
namespace WalkState

variable {G : TreeGhost} {B v d n : Nat} {outs₀ rest : List DfsOut} {done : List (DfsOut × Bool)}
  {hasVert : Bool} {o : DfsOut} {P : ItemId → Prop} {s : WalkState}

/-- The R content export at a tree-edge finish site (`FinishRShape`), in the shape of the recursion;
vacuous outside the 2-connected non-root case. -/
theorem WalkInvOut.rshape (h : WalkInvOut G B v d outs₀ done (o :: rest) hasVert n P s) :
    wp (walkOutPre v d o hasVert) (fun hv₁ s₁ => ∀ e cls child, o = .tree e cls child →
      wp (modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }) (fun _ s₂ =>
        wp (walkTree child (d + 1)) (fun _ s₃ => s.g.TwoConnected → ∀ dp, d = dp + 1 →
          FinishRShape G.dfs v d o s₁.tstack.length hv₁ s₃) s₂) s₁) s := by
  by_cases h2 : s.g.TwoConnected
  · rcases d with _ | dp
    · exact wp_of_forall fun hv₁ s₁ e cls child _ => by
        simp only [wp_modify]; exact wp_of_forall fun _ _ _ dp h => absurd h (by omega)
    · have hr := h.sites.rside h2 dp rfl
      unfold RSideOut at hr
      refine wp_mono _ hr.2.2 fun hv₁ s₁ hr e cls child ho => ?_
      subst ho
      simp only [wp_modify] at hr ⊢
      exact wp_mono _ hr.2 fun _ _ h _ _ _ => h
  · exact wp_of_forall fun hv₁ s₁ e cls child _ => by
      simp only [wp_modify]; exact wp_of_forall fun _ _ h2' => absurd h2' h2

/-- Admitted (Ranges content at the `finishEdge` site): the `FinishCanon` facts at the three
`finishTstackTop` sites of one `finishEdge` call follow from `CloseBase` there (checker field
`canon_*`, 0 violations). Obligation listed in `PROOF.md` §4.7. -/
theorem closeBase_canon {σ : List Nat} {n D curV d : Nat} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}
    {s : WalkState} (h : CloseBase σ n D curV d o origTstack hasVert s) :
    CloseCanon curV d o origTstack hasVert s := by
  sorry

/-- The layer-independent part of one `finishEdge` step, from `CloseBase` at the site. -/
theorem finish_core {σ : List Nat} {g : Graph} {n D m : Nat} {hv₁ : Bool} {s₀ s₃ : WalkState}
    {sts origs : List Nat} {P P' X : ItemId → Prop}
    (hcb : CloseBase σ n D v d o m hv₁ s₃) (hi : s₃.Inv' D) (hc : s₃.CanonInv)
    (hfull : s₃.Full g P' X) (how : OwnedD σ sts origs P d n s₃)
    (hv : v < g.nv) (he : o.e < g.ne) (hPe : ¬ P' (edgeItem g o.e))
    (hPv : hv₁ = false → ¬ P' (vertItem v)) (hPv' : hv₁ = true → P' (vertItem v))
    (hdest : o.cls.isTree = true → P (vertItem o.dest))
    (hsts : ∀ k, k ≤ d → sts[k]! ≤ sts[d]!) (hstn : sts[d]! ≤ n)
    (horigs : ∀ k, k ≤ d → origs[k]! ≤ origs[d]!) (horig : origs[d]! ≤ m)
    (hpath : ∀ k, k < d → s₃.stackVerts[k]! ≠ v) (hsvlt : ∀ k, k ≤ d → s₃.stackVerts[k]! < s₃.g.nv)
    (hT : Types g s₀) (hK : Keep (d + 1) 0 s₀ s₃) :
    wp (finishEdge v d o m hv₁) (fun hv' s' =>
      (hv₁ = true → hv' = true) ∧ s'.Inv' d ∧ Shape s' ∧ s'.RangesInv σ (n + 1) d ∧ s'.g = s₃.g ∧
      s'.CloseInv ∧ OwnedD σ sts origs P d (n + 1) s' ∧ (hv' = true → VertCover v s') ∧
      s'.Full g (fun i => P' i ∨ i = edgeItem g o.e ∨ (hv' = true ∧ i = vertItem v)) X ∧
      Keep (d + 1) 0 s₀ s' ∧ s'.CanonInv ∧ s'.ternarize = s₃.ternarize) s₃ := by
  have hD := hcb.hD
  have hstep := finishEdge_step hD hi hcb.rgs.2.1 hcb.guards hcb.book
  have hrg := finishEdge_ranges hD hcb.rgs.1 hcb.rgs.2.1 hcb.nodup hcb.rgs.2.2 hcb.guards hcb.book
    hcb.finishR
  have hcl := finishEdge_closeInv (CloseCtx.of_exports hcb (closeBase_content hcb))
  have hos : OwnSite σ n D v d o m hv₁ s₃ :=
    ⟨hcb.hD, hcb.nodup, hcb.rgs, hcb.guards, hcb.book, hcb.site.pos⟩
  have how' := finishEdge_ownedD hos how hsts hstn horigs horig hpath hsvlt hdest
  have hvc := finishEdge_vertCover hos
  have hF := finishEdge_full hfull hv he hPe hPv hPv'
    (finishEdge_sides hcb.guards hcb.book hi hcb.rgs.2.1 hD)
  have hK' := keep_finishEdge hT (j := 0) (by omega) v d o m hv₁ (vertItem_ne_zero v)
    (edgeItem_ne_zero g o.e) hK
  have hcn := finishEdge_canon hcb (closeBase_canon hcb) hc
  have htn := tern_finishEdge (b := s₃.ternarize) v d o m hv₁ rfl
  exact ⟨hF.2, hstep.1, hstep.2, hrg.1, hrg.2.2, hcl, how', hvc, hF.1, hK', hcn, htn⟩

theorem pushed_back_iff {g : Graph} {P : ItemId → Prop} {v e : Nat} {hv₁ hv' : Bool}
    (hh : hv₁ = true → hv' = true) (i : ItemId) :
    ((P i ∨ (hv₁ = true ∧ i = vertItem v)) ∨ i = edgeItem g e ∨ (hv' = true ∧ i = vertItem v)) ↔
      Pushed g (fun i => P i ∨ (hv' = true ∧ i = vertItem v)) [] [e] i := by
  unfold Pushed
  constructor
  · rintro ((h | ⟨h1, h2⟩) | h | h)
    · exact Or.inl (Or.inl h)
    · exact Or.inl (Or.inr ⟨hh h1, h2⟩)
    · exact Or.inr (Or.inr ⟨e, List.mem_singleton_self e, h⟩)
    · exact Or.inl (Or.inr h)
  · rintro ((h | h) | ⟨w, hw, _⟩ | ⟨e', he', h⟩)
    · exact Or.inl (Or.inl h)
    · exact Or.inr (Or.inr h)
    · exact nomatch hw
    · rw [List.mem_singleton.1 he'] at h; exact Or.inr (Or.inl h)

/-- Site `walkOut v d o hasVert`, from the primitives. -/
theorem pushed_tree_iff {g : Graph} {P : ItemId → Prop} {b₁ b₂ : Bool} {v e : Nat} {vs es : List Nat}
    (hb : b₁ = true → b₂ = true) (i : ItemId) :
    (Pushed g (fun i => P i ∨ (b₁ = true ∧ i = vertItem v)) vs es i ∨ i = edgeItem g e ∨
        (b₂ = true ∧ i = vertItem v)) ↔
      Pushed g (fun i => P i ∨ (b₂ = true ∧ i = vertItem v)) vs (e :: es) i := by
  constructor
  · rintro (((hp | ⟨h1, hi⟩) | hv | ⟨e', he', hi⟩) | he | ⟨h2, hi⟩)
    · exact Or.inl (Or.inl hp)
    · exact Or.inl (Or.inr ⟨hb h1, hi⟩)
    · exact Or.inr (Or.inl hv)
    · exact Or.inr (Or.inr ⟨e', List.mem_cons_of_mem _ he', hi⟩)
    · exact Or.inr (Or.inr ⟨e, List.mem_cons_self, he⟩)
    · exact Or.inl (Or.inr ⟨h2, hi⟩)
  · rintro ((hp | ⟨h2, hi⟩) | hv | ⟨e', he', hi⟩)
    · exact Or.inl (Or.inl (Or.inl hp))
    · exact Or.inr (Or.inr ⟨h2, hi⟩)
    · exact Or.inl (Or.inr (Or.inl hv))
    · rcases List.mem_cons.1 he' with rfl | he'
      · exact Or.inr (Or.inl hi)
      · exact Or.inl (Or.inr (Or.inr ⟨e', he', hi⟩))

theorem bbOut_step (o : DfsOut) (v d : Nat) (hasVert : Bool)
    (ih : match o with
      | .tree _ _ child => BbTree child (d + 1)
      | .back .. => True) :
    BbOut v d o hasVert := by
  intro G B outs₀ done rest n P s h
  obtain ⟨hndO, hdisjV, hvo, hvr, hendO, hdisjE⟩ := h.split_facts
  have S := h.sites
  have hE := h.earOut
  have hg := S.guards; unfold GuardsOut at hg
  have hb := S.book; unfold BookOut at hb
  have hf := S.frontiers; unfold FrontiersOut at hf
  have hcb := S.cb; unfold CbOut at hcb
  have hRw := h.rshape
  have hE₂ := hE.2
  have hL := walkOut_stLive h
  have hgs : s.g = G.g := h.g_eq
  have hoe : o.e ∈ o.edges := by cases o <;> simp [DfsOut.edges, DfsOut.e]
  have helt : o.e < G.g.ne := h.e_lt_o o.e hoe
  have hPeo : ¬ P (edgeItem G.g o.e) := h.Pe_o o.e hoe
  have hd₀ : d < s.stackDir.size := by have := h.height; omega
  have hB : ∀ t ∈ segsStack G.segs, t.vStart ≠ v := fun t ht h' =>
    h.base_out t ht (h' ▸ List.mem_cons_self)
  have hsvne : ∀ k, k < d → s.stackVerts[k]! ≠ v := by
    intro k hk heq
    have hm : s.stackVerts[k]! ∈ G.anc := List.mem_of_getElem? (h.anc_sv k hk)
    rw [heq] at hm
    exact (List.nodup_append.1 h.nodup).2.2 _ hm _ List.mem_cons_self rfl
  have hsvlt : ∀ k, k ≤ d → s.stackVerts[k]! < G.g.nv := by
    intro k hk
    rcases Nat.lt_or_eq_of_le hk with hk | rfl
    · exact h.anc_lt _ (List.mem_of_getElem? (h.anc_sv k hk))
    · rw [h.sv_d]; exact h.v_lt
  obtain ⟨new, hts, hR, hnew⟩ := h.stPre.read
  rw [walkOut_eq, wp_bind] at hE₂ hL ⊢
  refine wp_mono _ (wp_and (WalkInvOut.walkOutPre_out h) (wp_and hg (wp_and hb.2 (wp_and hf
    (wp_and hcb (wp_and hRw (wp_and hE₂ hL))))))) ?_
  intro hv₁ s₁ ⟨⟨push, hp⟩, hg₁, hb₁, hf₁, hcb₁, hr₁, hE₁, hL₁⟩
  unfold walkOutRest at hE₁ hL₁ ⊢
  rw [wp_bind, wp_tstackSize] at hE₁ hL₁ ⊢
  obtain ⟨new₁, hts₁, hR₁, hI₁⟩ := hp.pre.read
  have hhv₁ := hp.hv
  have hsd₁ := hp.pre.sd
  have hnopush := hp.pre.nopush
  have hdirs₁ : DirsOf s₁ d = DirsOf s d := hp.pre.dirs
  have hvs₁ : ∀ t ∈ s₁.tstack,
      t.vStart ∈ v :: DfsOut.vertsList (done.map (·.1)) ∨ ∃ t₀ ∈ segsStack G.segs, t.vStart = t₀.vStart :=
    fun t ht => (hp.pre.tstack t ht).elim (fun h => Or.inl (h ▸ List.mem_cons_self)) (h.stPre.vStart t)
  have hseg₁ : SegRead s₁.items G.segs := by rw [hp.pre.items]; exact h.stPre.segRead
  have hgeq₁ : s₁.g = s.g := hp.pre.geq
  have hsgeq : s₁.g = G.g := hgeq₁.trans hgs
  have hit₁ : s₁.items = s.items := hp.pre.items
  have hpushb : (!hasVert && decide (o.cls.lowval d < d) && o.cls.isType1) = push := by
    cases hpu : push
    · apply Bool.eq_false_iff.2
      intro hb'
      simp only [Bool.and_eq_true, Bool.not_eq_true', decide_eq_true_eq] at hb'
      exact absurd (hp.hpush.2 ⟨hb'.1.1, hb'.1.2, hb'.2⟩) (by rw [hpu]; decide)
    · obtain ⟨h1, h2, h3⟩ := hp.hpush.1 hpu
      simp [h1, h2, h3]
  rw [hpushb] at hR₁ hnopush
  generalize hx : (if d ≤ o.cls.lowval d then false else !s.stackDir[o.cls.lowval d]!) = x at hsd₁ hR₁
  have hxd : s₁.stackDir[d]! = x := by rw [hsd₁]; exact Array.getElem!_set!_self _ _ _ hd₀
  have hnew₁ : push = false → hasVert = false → new₁ = [] := fun hpf h0 =>
    (List.append_cancel_right (hts₁.symm.trans ((hnopush hpf).trans hts))).trans (hnew h0)
  have hPv : hv₁ = false → ¬ (P (vertItem v) ∨ (hv₁ = true ∧ vertItem v = vertItem v)) := fun h0 hp' =>
    hp'.elim (h.Pcur (by rw [hhv₁] at h0; exact (Bool.or_eq_false_iff.1 h0).1))
      (fun h1 => by rw [h0] at h1; exact Bool.false_ne_true h1.1)
  have hPv' : hv₁ = true → (P (vertItem v) ∨ (hv₁ = true ∧ vertItem v = vertItem v)) :=
    fun h1 => Or.inr ⟨h1, rfl⟩
  have hPe' : ¬ (P (edgeItem G.g o.e) ∨ (hv₁ = true ∧ edgeItem G.g o.e = vertItem v)) := fun hp' =>
    hp'.elim hPeo fun h1 => vertItem_ne_edgeItem h.v_lt o.e h1.2.symm
  have hpath₁ : ∀ k, k < d → s₁.stackVerts[k]! ≠ v := fun k hk => by
    rw [hp.keep.svlo k (by omega)]; exact hsvne k hk
  have hsvlt₁ : ∀ k, k ≤ d → s₁.stackVerts[k]! < s₁.g.nv := fun k hk => by
    rw [hp.keep.svlo k (by omega), hsgeq]; exact hsvlt k hk
  cases o with
  | back e dest cls =>
    try simp only at hg₁ hb₁ hf₁ hcb₁ hr₁ hE₁ hL₁
    simp only
    obtain ⟨sub, base, hlen, hEf⟩ := hb₁.ear
    have hnt : (DfsOut.back e dest cls).cls.isTree = false := by
      cases hh : (DfsOut.back e dest cls).cls.isTree
      · rfl
      · obtain ⟨_, _, _, h'⟩ := hb₁.tree.1 hh; cases h'
    have hsub : sub = [] := hEf.back_nil hnt
    subst hsub
    have hbase : base = new₁ ++ segsStack G.segs := by simpa [hts₁] using hEf.tstack.symm
    subst hbase
    have hD : d = if (DfsOut.back e dest cls).cls.isTree then d + 1 else d := by simp [hnt]
    have hcore := finish_core hcb₁ hp.inv hp.canon hp.full hp.owned h.v_lt helt hPe' hPv hPv'
      (fun ht => absurd ht (by rw [hnt]; decide)) h.sts_le h.sts_n h.origs_le
      (hp.owned.len d (le_refl _)) hpath₁ hsvlt₁ h.types hp.keep
    rw [hts₁] at hg₁ hb₁ hf₁ hcb₁ hcore hE₁ hL₁ ⊢
    obtain ⟨hhv', hi', hs', hrg', hgs', hc', ho', hvc', hF', hK', hcn', htn'⟩ := hcore
    have hgR : ((finishEdge v d (.back e dest cls) (new₁ ++ segsStack G.segs).length hv₁).run s₁).2.g =
      s₁.g := hgs'
    obtain ⟨hvF, hC'⟩ := hE₁
    obtain ⟨hlc, hll⟩ := hL₁
    have hvl : DfsOut.vertsList (done.map (·.1) ++ [DfsOut.back e dest cls]) =
        DfsOut.vertsList (done.map (·.1)) := by
      simp [DfsOut.vertsList_append, DfsOut.vertsList]
    have key : StPre G.g G.prev G.fs G.segs v d (done.map (·.1) ++ [.back e dest cls])
          ((finishEdge v d (.back e dest cls) (new₁ ++ segsStack G.segs).length hv₁).run s₁).1
          ((finishEdge v d (.back e dest cls) (new₁ ++ segsStack G.segs).length hv₁).run s₁).2 ∧
        (s.g.TwoConnected → ROutCtx G.dfs G.F B v d outs₀
          ((finishEdge v d (.back e dest cls) (new₁ ++ segsStack G.segs).length hv₁).run s₁).1
          ((finishEdge v d (.back e dest cls) (new₁ ++ segsStack G.segs).length hv₁).run s₁).2) := by
      by_cases hge : d ≤ (DfsOut.back e dest cls).cls.lowval d
      · have hge' : d ≤ cls.lowval d := hge
        have hok := ear_boundary hge hg₁ hp.inv hp.shape hb₁ hD
        obtain ⟨hr1, hrg, hrsd, hrts, hrI, hrfr⟩ :=
          finishBoundary_st hEf hp.inv hp.shape hD hok hb₁ hge hsgeq StRead.nil hI₁
            (fun ht => absurd ht (by rw [hnt]; decide))
        have hpf : push = false := by rw [← hpushb]; simp [Nat.not_lt.mpr hge]
        have hhv : hv₁ = hasVert := by rw [hhv₁, hpf, Bool.or_false]
        have hdr : DirsOf ((finishEdge v d (.back e dest cls) (new₁ ++ segsStack G.segs).length hv₁).run
            s₁).2 d = DirsOf s d := by
          rw [DirsOf_congr (s := s₁) (fun k _ => by rw [hrsd]), hdirs₁]
        rw [hpf] at hR₁
        simp only [Bool.false_eq_true, ↓reduceIte, List.append_nil] at hR₁
        rw [hnt] at hrI
        simp only [Bool.false_eq_true, ↓reduceIte, List.append_nil] at hrI
        refine ⟨⟨⟨new₁, hrts, ?_, fun h0 => hnew₁ hpf ?_⟩, ?_, ?_, ?_, ?_⟩, fun h2 => ?_⟩
        · rw [hdr, refOuts_snoc, h.stPre.hv, refOut_boundary_back hge']
          simp only [List.append_nil]
          exact StRead.congr (fun x hx => fun y hb =>
            ((hrfr x (mem_readStack_append.2 (Or.inl hx)) y hb).imp_right And.left)) hR₁
        · rw [hr1, hhv] at h0; exact h0
        · rw [hdr, refOuts_snoc, h.stPre.hv, refOut_boundary_back hge', hr1, hhv]
        · exact hseg₁.congr fun x hx => fun y hb =>
            ((hrfr x (mem_readStack_append.2 (Or.inr hx)) y hb).imp_right And.left)
        · rw [hrts, ← hts₁, hvl]; exact hvs₁
        · rw [hdr, refOuts_snoc, h.stPre.hv, refOut_boundary_back hge']
          simpa only [List.append_nil] using hrI
        · -- R: a boundary back edge only occurs at the root in the 2-connected case
          have R := h.r h2
          have hd0 : d = 0 := by
            by_contra hne
            obtain ⟨dp, rfl⟩ : ∃ dp, d = dp + 1 := ⟨d - 1, by omega⟩
            have hr := S.rside h2 dp rfl
            unfold RSideOut at hr
            exact absurd hr.1 (Nat.not_lt.2 hge)
          obtain ⟨hts0, hhv0, hF0⟩ := R.root hd0
          have hts₁' : s₁.tstack = [] := by rw [hnopush hpf, hts0]
          have hnil := List.append_eq_nil_iff.1 (hts₁'.symm.trans hts₁).symm
          have hts' : ((finishEdge v d (.back e dest cls) (new₁ ++ segsStack G.segs).length hv₁).run
              s₁).2.tstack = [] := by rw [hrts, hnil.1, hnil.2]; rfl
          have hhvF : ((finishEdge v d (.back e dest cls) (new₁ ++ segsStack G.segs).length hv₁).run
              s₁).1 = false := by rw [hr1, hhv, hhv0]
          have hg'' := hgR.trans hgeq₁
          exact {
            wf := by rw [hg'']; exact R.wf
            spec := by rw [hg'']; exact R.spec
            rooted := by rw [hg'']; exact R.rooted
            outs_v := R.outs_v
            sub := R.sub
            chain := ⟨by rw [hK'.svlo d (by omega)]; exact R.chain.1,
              fun k hk => by rw [hK'.svlo k (by omega)]; exact R.chain.2 k hk⟩
            rwalk := ⟨fun f hf => absurd hf (by rw [hF0]; exact List.not_mem_nil),
              ⟨fun t ht => absurd ht (by rw [hts']; exact List.not_mem_nil),
                by rw [hts']; exact List.Pairwise.nil⟩⟩
            fB := R.fB
            B_le := by rw [hts']; have := R.B_le; rw [hts0] at this; exact this
            hvB := fun hh => absurd hh (by rw [hhvF]; decide)
            root := fun _ => ⟨hts', hhvF, hF0⟩
            skel := by
              have := keepsR_finishBoundary (D := d) v d (.back e dest cls)
                (new₁ ++ segsStack G.segs).length hv₁ hge hp.shape
                (by rw [hgeq₁]; exact hgs ▸ h.v_lt) (by rw [hsgeq]; exact helt) hok
                (by rw [hgeq₁, hit₁]; exact R.skel)
              rw [hgR]; exact this }
      · have hlt : (DfsOut.back e dest cls).cls.lowval d < d := Nat.lt_of_not_le hge
        have hlt' : cls.lowval d < d := hlt
        have hge' : ¬ d ≤ cls.lowval d := hge
        have hpr : push = (!hasVert && cls.isType1) := by rw [← hpushb]; simp [DfsOut.cls, hlt']
        have hhv : hv₁ = (hasVert || cls.isType1) := by
          rw [hhv₁, hpr]; cases hasVert <;> cases cls.isType1 <;> rfl
        have hxv : x = !s.stackDir[cls.lowval d]! := by rw [← hx]; simp [DfsOut.cls, hge']
        have hpre' : hv₁ = false → new₁ = [] ∧
            ((refOuts G.g v d (DirsOf s d) (done.map (·.1)) false).1 ++
              if push = true then [⟨x, [vertItem v]⟩] else []) = [] := by
          intro h0
          rw [hhv, Bool.or_eq_false_iff] at h0
          have hpf : push = false := by simp [hpr, h0.1, h0.2]
          refine ⟨hnew₁ hpf h0.1, ?_⟩
          rw [hpf, (refOuts_hv_false _ _ (h.stPre.hv.trans h0.1)).2]; rfl
        obtain ⟨hr1, hrg, hdr', ⟨new', hrts, hR'⟩, hrI, hseg', hvs'⟩ :=
          stRet_finish (D := d) hEf hp.inv hp.shape hD hg₁ hb₁ hlt hB hpre' StRead.nil hR₁ hI₁ hseg₁
            hp.full
        rw [hdirs₁] at hdr'
        rw [hxd, hxv, Bool.not_not, hsgeq, hpr, hnt] at hR'
        simp only [Bool.false_eq_true, ↓reduceIte] at hR'
        refine ⟨?_, fun h2 => ?_⟩
        · subst hhv
          refine ⟨⟨new', hrts, ?_, fun h0 => by rw [hr1] at h0; cases h0⟩, ?_, hseg', ?_, ?_⟩
          · rw [hdr', refOuts_snoc, h.stPre.hv, refOut_ret_back' hlt', DirsOf_getD s hlt']
            simpa only [List.append_assoc, DfsOut.e] using hR'
          · rw [hdr', refOuts_snoc, h.stPre.hv, refOut_ret_back' hlt', hr1]
          · intro t ht
            rw [hvl]
            rcases hvs' t ht with h' | ⟨hT, -, -⟩ | ⟨t₀, ht₀, h'⟩
            · exact Or.inl (h' ▸ List.mem_cons_self)
            · exact Bool.noConfusion ((show cls.isTree = true from hT).symm.trans
                (show cls.isTree = false from hnt))
            · exact (hvs₁ t₀ ht₀).imp (h' ▸ id) fun ⟨t₁, ht₁, h₁⟩ => ⟨t₁, ht₁, h'.trans h₁⟩
          · rw [hdr', refOuts_snoc, h.stPre.hv, refOut_ret_back' hlt']
            simpa only [List.append_nil] using hrI
        · have R := h.r h2
          obtain ⟨lv, kind, ho, hl⟩ := ret_of_lowval_lt hlt
          have hk : kind = .backEdge := by
            cases kind with
            | backEdge => rfl
            | type1Child => exact absurd (hb₁.tree.1 (by rw [ho]; rfl)) (by simp)
            | type2Child => exact absurd (hb₁.tree.1 (by rw [ho]; rfl)) (by simp)
          subst hk
          have hhvt : hv₁ = true := by
            rw [hhv₁]; cases hhV : hasVert
            · rw [hp.hpush.2 ⟨hhV, hlt, by rw [ho]; rfl⟩]; rfl
            · rfl
          have hanc₁ : AncChain G.dfs v d s₁ :=
            ⟨by rw [hp.keep.svlo d (by omega)]; exact R.chain.1,
              fun k hk => by rw [hp.keep.svlo k (by omega)]; exact R.chain.2 k hk⟩
          have h2₁ : s₁.g.TwoConnected := by rw [hgeq₁]; exact h2
          have hsp₁ : G.dfs.Spec s₁.g := by rw [hgeq₁]; exact R.spec
          have hrt₁ : G.dfs.Rooted s₁.g := by rw [hgeq₁]; exact R.rooted
          obtain ⟨hW₁, hK₁, hn₁⟩ := hp.r h2
          have hK₁' := hK₁.1
          rw [hts₁] at hK₁' hn₁
          have hok := finishOk_of_guards ho hl hg₁ hEf rfl hp.inv hp.shape hD hb₁.v_lt hb₁.e_lt hb₁.q
            (hb₁.ends lv .backEdge ho) hb₁.vert
          have hfront := finishEdge_frontier hb₁ hp.inv hp.shape hD
          have hhd : G.dfs.depth v = d := by
            have := (hanc₁.2 d (le_refl _)).2; rwa [hanc₁.1] at this
          have hRf : s₁.RInvFront G.dfs v d (new₁ ++ segsStack G.segs).length := hW₁.top.toFront _
          have hbot := finishEdge_bot_keep v d lv .backEdge _ (new₁ ++ segsStack G.segs).length B hv₁ ho
            hl hEf rfl (fun _ => hhvt) hK₁' hn₁
          have hBr : B ≤ ((finishEdge v d (.back e dest cls) (new₁ ++ segsStack G.segs).length hv₁).run
            s₁).2.tstack.length := hbot.1.1
          have hBr' : (new₁ ++ segsStack G.segs).length + (if hv₁ then 0 else 1) ≤
            ((finishEdge v d (.back e dest cls) (new₁ ++ segsStack G.segs).length hv₁).run
              s₁).2.tstack.length := hbot.2
          have hle₁ : s₁.tstack.length ≤ (new₁ ++ segsStack G.segs).length := Nat.le_of_eq (congrArg List.length hts₁)
          have hg'' := hgR.trans hgeq₁
          exact {
            wf := by rw [hg'']; exact R.wf
            spec := by rw [hg'']; exact R.spec
            rooted := by rw [hg'']; exact R.rooted
            outs_v := R.outs_v
            sub := R.sub
            chain := ⟨by rw [hK'.svlo d (by omega)]; exact R.chain.1,
              fun k hk => by rw [hK'.svlo k (by omega)]; exact R.chain.2 k hk⟩
            rwalk := ⟨fun f hf => ⟨by have := R.fB f hf; omega,
                finishEdge_rInvG_base f.1 f.2.1 v d lv .backEdge _ _ f.2.2 hv₁ ho hl hb₁.v_lt hp.inv
                  hp.shape hok hfront (fun h' => absurd rfl h') (fun _ => ⟨hhvt, hle₁⟩)
                  (fun h' => absurd rfl h') (by have := R.fB f hf; omega)
                  (fun _ => by have := R.fB f hf; have := hn₁ hhvt; omega) (hW₁.frames f hf).2⟩,
              finishEdge_rInvTop v d lv .backEdge _ _ hv₁ ho hl hb₁.v_lt hp.inv hp.shape hok hg₁ hfront
                (fun h' => absurd rfl h') (fun _ => ⟨hhvt, hle₁⟩) h2₁ hsp₁ hrt₁ hhd hanc₁.1
                hanc₁.2 hRf⟩
            fB := R.fB
            B_le := hBr
            hvB := fun _ => by
              have h2 := hn₁ hhvt
              have h0 : (if hv₁ then 0 else 1) = 0 := by rw [hhvt]; rfl
              rw [h0] at hBr'
              omega
            root := fun hd0 => absurd hlt (by rw [hd0]; exact Nat.not_lt_zero _)
            skel := by
              have := keepsR_finishEdge_site (dfs := G.dfs) (D := d) v d lv .backEdge _ _ hv₁ ho hl
                hb₁.v_lt hp.inv hp.shape hok hD hfront (fun h' => absurd rfl h') h2₁ hsp₁ hrt₁ hhd
                hanc₁.1 hanc₁.2 hRf hEf.q_root (fun ht => absurd ht (by rw [hnt]; decide)) hcb₁
                (by rw [hgeq₁, hit₁]; exact R.skel)
              rw [hgR]; exact this }
    obtain ⟨hSt, hRR⟩ := key
    refine ⟨fun h0 => hhv' (by rw [hhv₁, h0]; rfl), hK', hvF, ?_⟩
    have hhf : ((finishEdge v d (.back e dest cls) (new₁ ++ segsStack G.segs).length hv₁).run s₁).1 =
        false → hasVert = false := fun h0 => by
      cases hhV : hasVert
      · rfl
      · exact absurd (hhv' (by rw [hhv₁, hhV]; rfl)) (by rw [h0]; decide)
    exact {
      full := hF'.monoP (pushed_back_iff hhv')
      d_anc := h.d_anc
      split := by rw [List.map_append, List.append_assoc]; exact h.split
      sorted := h.sorted
      wf := h.wf
      ends := h.ends
      nodup := h.nodup
      anc_lt := h.anc_lt
      v_lt := h.v_lt
      w_lt := h.w_lt
      e_lt := h.e_lt
      enodup := h.enodup
      comp := h.comp
      pe_anc := h.pe_anc
      sv_size := hK'.sv.trans h.sv_size
      sd_size := hK'.sd.trans h.sd_size
      anc_sv := fun k hk => by rw [hK'.svlo k (by omega)]; exact h.anc_sv k hk
      sv_d := by rw [hK'.svlo d (by omega)]; exact h.sv_d
      ear := hC'
      inv := hi'
      shape := hs'
      ranges := hrg'
      σ_nodup := h.σ_nodup
      σ_lt := by rw [hgs', hgeq₁]; exact h.σ_lt
      post := by have := h.post; rw [DfsOut.edgePostorderList_cons] at this; exact this.right
      path := by
        have := h.path; rw [DfsOut.edgePostorderList_cons, List.length_append, ← Nat.add_assoc] at this
        exact this
      close := hc'
      canon := hcn'
      tern := htn'.trans hp.tern
      P_past := pushed_past (fun e he hp => by
          rcases hp with hp | ⟨_, hp⟩
          · exact h.P_past e he hp
          · exact absurd hp.symm (vertItem_ne_edgeItem h.v_lt e))
        h.w_lt_o (DfsOut.edges_subset_block _) h.post_o h.σ_nodup
      owned := ho'.mono fun i hi => Or.inl (Or.inl hi)
      sts_len := h.sts_len
      origs_len := h.origs_len
      sts_le := h.sts_le
      sts_n := Nat.le_trans h.sts_n (Nat.le_add_right _ _)
      origs_le := h.origs_le
      Pv := by
        rintro w hw ((hp | ⟨_, hp⟩) | ⟨w', hw', hp⟩ | ⟨e, _, hp⟩)
        · exact h.Pv_rest w hw hp
        · exact hvr (WalkM.vertItem_inj hp ▸ hw)
        · exact hdisjV w' hw' (WalkM.vertItem_inj hp ▸ hw)
        · exact vertItem_ne_edgeItem (h.w_lt_rest w hw) e hp
      Pe := by
        rintro e he ((hp | ⟨_, hp⟩) | ⟨w', hw', hp⟩ | ⟨e', he', hp⟩)
        · exact h.Pe_rest e he hp
        · exact vertItem_ne_edgeItem h.v_lt e hp.symm
        · exact vertItem_ne_edgeItem (h.w_lt_o w' hw') e hp.symm
        · exact hdisjE e' he' (edgeItem_inj hp ▸ he)
      Pcur := by
        rintro h0 ((hp | ⟨h1, _⟩) | ⟨w', hw', hp⟩ | ⟨e, _, hp⟩)
        · exact h.Pcur (hhf h0) hp
        · rw [h0] at h1; cases h1
        · exact hvo (WalkM.vertItem_inj hp ▸ hw')
        · exact vertItem_ne_edgeItem h.v_lt e hp
      Pcur' := fun h1 => Or.inl (Or.inr ⟨h1, rfl⟩)
      vcover := hvc'
      r := fun h2 => hRR (by rw [hK'.g] at h2; exact h2)
      d_fs := h.d_fs
      height := by rw [hK'.sd]; exact h.height
      base_out := h.base_out
      stPre := by rw [List.map_append]; exact hSt
      segs_len := h.segs_len
      live_cur := by rw [List.map_append]; exact hlc
      live_lower := hll }
  | tree e cls child =>
    cases child with
    | node c couts =>
    try simp only [wp_bind, wp_modify] at hg₁ hb₁ hcb₁ hE₁ hL₁
    simp only [wp_bind, wp_modify]
    have hr₁' := hr₁ e cls (.node c couts) rfl
    simp only [wp_modify] at hr₁'
    obtain ⟨sv', sd', new₁', hts₁', hW⟩ := PreOut.child h hp
    have hnew_eq : new₁ = new₁' := List.append_cancel_right (hts₁.symm.trans hts₁')
    subst hnew_eq
    have hih := ih _ _ hW
    have hkb : wp (walkTree (.node c couts) (d + 1))
        (fun _ s₃ => ∀ k, k < d + 1 → s₃.stackDir[k]! = s₁.stackDir[k]!)
        { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne } :=
      fun k hk => walk_keepsBelow.1 (.node c couts) (d + 1) _ k hk
    refine wp_mono _ (wp_and hih (wp_and hg₁.2 (wp_and hb₁.2 (wp_and hcb₁.2
      (wp_and hr₁' (wp_and hE₁ (wp_and hL₁ hkb))))))) ?_
    intro _ s₃ ⟨hend, hg₃, hb₃, hcb₃, hr₃, hE₃, hL₃, hkb₃⟩
    have hit : (DfsOut.tree e cls (.node c couts)).cls.isTree = true := hb₃.tree.2 ⟨_, _, _, rfl⟩
    have hD : d + 1 = if (DfsOut.tree e cls (.node c couts)).cls.isTree then d + 1 else d := by simp [hit]
    have hK₂ : Keep (d + 1) 0 s₁ { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne } :=
      Keep.frame rfl rfl rfl (by simp) (fun _ _ => rfl) rfl
    have hK₃ : Keep (d + 1) 0 s s₃ := (hp.keep.trans hK₂).trans hend.keep
    have hpath₃ : ∀ k, k < d → s₃.stackVerts[k]! ≠ v := fun k hk => by
      rw [hK₃.svlo k (by omega)]; exact hsvne k hk
    have hsvlt₃ : ∀ k, k ≤ d → s₃.stackVerts[k]! < s₃.g.nv := fun k hk => by
      rw [hK₃.svlo k (by omega), hK₃.g, hgs]; exact hsvlt k hk
    have hgeq₃ : s₃.g = s.g := hK₃.g
    have hsgeq₃ : s₃.g = G.g := hgeq₃.trans hgs
    have hcv : c ∈ (DfsTree.node c couts).verts := List.mem_cons_self
    have hnσ : n + (DfsTree.node c couts).edgePostorder.length ≤ G.σ.length := by
      obtain ⟨pre, post, hpre, hσ⟩ := h.post
      rw [hσ]
      simp only [List.length_append, DfsOut.edgePostorderList, List.length_cons, hpre]
      omega
    have how₃ := hend.owned.exit h.sts_len h.origs_len (Nat.le_add_right _ _) hend.ranges.2.2 hnσ
      (by rw [hend.sv_v]; exact hend.vertCover) (by rw [hend.sv_v]; exact Or.inr (Or.inl ⟨c, hcv, rfl⟩))
    have hendO' : (e :: (DfsTree.node c couts).edges).Nodup := hendO
    have hPe₃ : ¬ Pushed G.g (fun i => P i ∨ (hv₁ = true ∧ i = vertItem v)) (DfsTree.node c couts).verts
        (DfsTree.node c couts).edges (edgeItem G.g e) := by
      rintro (hp' | ⟨w, hw, hw'⟩ | ⟨e', he', he''⟩)
      · exact hPe' hp'
      · exact vertItem_ne_edgeItem (h.w_lt_o w hw) e hw'.symm
      · exact (List.nodup_cons.1 hendO').1 ((edgeItem_inj he'').symm ▸ he')
    have hPv₃ : hv₁ = false → ¬ Pushed G.g (fun i => P i ∨ (hv₁ = true ∧ i = vertItem v))
        (DfsTree.node c couts).verts (DfsTree.node c couts).edges (vertItem v) := by
      rintro h0 (hp' | ⟨w, hw, hw'⟩ | ⟨e', _, he''⟩)
      · exact hPv h0 hp'
      · exact hvo (WalkM.vertItem_inj hw' ▸ hw)
      · exact vertItem_ne_edgeItem h.v_lt e' he''
    have hPv'₃ : hv₁ = true → Pushed G.g (fun i => P i ∨ (hv₁ = true ∧ i = vertItem v))
        (DfsTree.node c couts).verts (DfsTree.node c couts).edges (vertItem v) :=
      fun h1 => Or.inl (hPv' h1)
    have hcore := finish_core hcb₃ hend.inv hend.canon hend.full how₃ h.v_lt helt hPe₃ hPv₃ hPv'₃
      (fun _ => Or.inr (Or.inl ⟨c, hcv, rfl⟩)) h.sts_le (Nat.le_trans h.sts_n (Nat.le_add_right _ _))
      h.origs_le (hp.owned.len d (le_refl _)) hpath₃ hsvlt₃ h.types hK₃
    have hR₂ : s.g.TwoConnected →
        RWalk G.dfs ((v, d, s₁.tstack.length) :: G.F) c (d + 1) s₃ ∧
          BotKeep s₁.tstack.length s₁ s₃ ∧ Items.RSkelInv s₃.g s₃.items := fun h2 => by
      have hr := hend.r (by show s₁.g.TwoConnected; rw [hgeq₁]; exact h2) d rfl
      exact ⟨hr.1, hr.2.1, hr.2.2.2⟩
    rw [hts₁] at hg₃ hb₃ hcb₃ hcore hE₃ hL₃ hr₃ hR₂ ⊢
    obtain ⟨hhv', hi', hs', hrg', hgs', hc', ho', hvc', hF', hK', hcn', htn'⟩ := hcore
    have hgR : ((finishEdge v d (.tree e cls (.node c couts)) (new₁ ++ segsStack G.segs).length hv₁).run
      s₃).2.g = s₃.g := hgs'
    obtain ⟨hvF, hC'⟩ := hE₃
    obtain ⟨hlc, hll⟩ := hL₃
    have hvl : DfsOut.vertsList (done.map (·.1) ++ [DfsOut.tree e cls (.node c couts)]) =
        DfsOut.vertsList (done.map (·.1)) ++ (DfsTree.node c couts).verts := by
      simp [DfsOut.vertsList_append, DfsOut.vertsList]
    have hdirs₃ : DirsOf s₃ (d + 1) = DirsOf s d ++ [x] := by
      rw [DirsOf_congr (s := s₁) hkb₃, DirsOf_succ, hxd, hdirs₁]
    have hdirs₃' : DirsOf s₃ d = DirsOf s d := by
      rw [DirsOf_congr (s := s₁) (fun k hk => hkb₃ k (Nat.lt_succ_of_lt hk))]; exact hdirs₁
    have hxd₃ : s₃.stackDir[d]! = x := by rw [hkb₃ d (Nat.lt_succ_self d)]; exact hxd
    obtain ⟨sub, base, hlen, hEf⟩ := hb₃.ear
    obtain ⟨new₃, hts₃, hR₃⟩ := hend.st.read
    have hI₃ := hend.st.items
    have hlive₃ := hend.live
    simp only [List.length_append, List.length_singleton, ← h.d_fs, hdirs₃, segsStack_cons]
      at hts₃ hR₃ hI₃ hlive₃
    rw [simBlocks_frame h.d_fs] at hI₃
    have hRq₃ : StRead s₃.items new₁ ((refOuts G.g v d (DirsOf s d) (done.map (·.1)) false).1 ++
        if push = true then [⟨x, [vertItem v]⟩] else []) := by
      have := hend.segRead _ List.mem_cons_self; rwa [hxd] at this
    have hseg₃' : SegRead s₃.items G.segs := fun sg hsg => hend.segRead sg (List.mem_cons_of_mem _ hsg)
    obtain ⟨hsub, hbase⟩ := List.append_inj' (hEf.tstack.symm.trans hts₃) hlen
    subst hsub hbase
    have hvs₃ := hend.vStart
    simp only [segsStack_cons] at hvs₃
    have hvs₃' : ∀ t ∈ s₃.tstack,
        t.vStart ∈ v :: DfsOut.vertsList (done.map (·.1) ++ [DfsOut.tree e cls (.node c couts)]) ∨
          ∃ t₀ ∈ segsStack G.segs, t.vStart = t₀.vStart := by
      intro t ht
      rw [hvl]
      rcases hvs₃ t ht with h' | ⟨t₀, ht₀, h'⟩
      · exact Or.inl (List.mem_cons_of_mem _ (List.mem_append_right _ h'))
      · rw [← hts₁] at ht₀
        exact (hvs₁ t₀ ht₀).imp (fun h'' => h' ▸ mem_cons_append_left h'')
          fun ⟨t₁, ht₁, h₁⟩ => ⟨t₁, ht₁, h'.trans h₁⟩
    have key : StPre G.g G.prev G.fs G.segs v d (done.map (·.1) ++ [.tree e cls (.node c couts)])
          ((finishEdge v d (.tree e cls (.node c couts)) (new₁ ++ segsStack G.segs).length hv₁).run s₃).1
          ((finishEdge v d (.tree e cls (.node c couts)) (new₁ ++ segsStack G.segs).length hv₁).run s₃).2 ∧
        (s.g.TwoConnected → ROutCtx G.dfs G.F B v d outs₀
          ((finishEdge v d (.tree e cls (.node c couts)) (new₁ ++ segsStack G.segs).length hv₁).run s₃).1
          ((finishEdge v d (.tree e cls (.node c couts)) (new₁ ++ segsStack G.segs).length hv₁).run s₃).2) := by
      by_cases hge : d ≤ (DfsOut.tree e cls (.node c couts)).cls.lowval d
      · have hge' : d ≤ cls.lowval d := hge
        have hxv : x = false := by rw [← hx]; simp [DfsOut.cls, hge']
        have hok := ear_boundary hge hg₃ hend.inv hend.shape hb₃ hD
        subst hxv
        have hlive := hlive₃ sub hts₃
        simp only [openBlock, openBlockP_snoc_bd G.g _ ⟨v, done.map (·.1), .tree e cls (.node c couts)⟩ _
          G.fs 0 (by rw [Nat.zero_add, ← h.d_fs]; exact hge')] at hlive
        obtain ⟨hr1, hrg, hrsd, hrts, hrI, hrfr⟩ :=
          finishBoundary_st hEf hend.inv hend.shape hD hok hb₃ hge hsgeq₃ hR₃ hI₃ (fun _ => hlive)
        have hpf : push = false := by rw [← hpushb]; simp [Nat.not_lt.mpr hge]
        have hhv : hv₁ = hasVert := by rw [hhv₁, hpf, Bool.or_false]
        have hdr : DirsOf ((finishEdge v d (.tree e cls (.node c couts)) (new₁ ++ segsStack G.segs).length
            hv₁).run s₃).2 d = DirsOf s d := by
          rw [DirsOf_congr (s := s₃) (fun k _ => by rw [hrsd]), hdirs₃']
        rw [hpf] at hRq₃
        simp only [Bool.false_eq_true, ↓reduceIte, List.append_nil] at hRq₃
        rw [hit] at hrI
        simp only [↓reduceIte] at hrI
        refine ⟨⟨⟨new₁, hrts, ?_, fun h0 => hnew₁ hpf ?_⟩, ?_, ?_, ?_, ?_⟩, fun h2 => ?_⟩
        · rw [hdr, refOuts_snoc, h.stPre.hv, refOut_boundary_tree hge']
          simp only [List.append_nil]
          exact StRead.congr (fun x hx => fun y hb =>
            ((hrfr x (mem_readStack_append.2 (Or.inl hx)) y hb).imp_right And.left)) hRq₃
        · rw [hr1, hhv] at h0; exact h0
        · rw [hdr, refOuts_snoc, h.stPre.hv, refOut_boundary_tree hge', hr1, hhv]
        · exact hseg₃'.congr fun x hx => fun y hb =>
            ((hrfr x (mem_readStack_append.2 (Or.inr hx)) y hb).imp_right And.left)
        · rw [hrts, ← hts₁]
          intro t ht; rw [hvl]; exact (hvs₁ t ht).imp_left mem_cons_append_left
        · rw [hdr, refOuts_snoc, h.stPre.hv, refOut_boundary_tree hge']
          simpa only [List.append_assoc, DfsOut.dest] using hrI
        · have R := h.r h2
          have hd0 : d = 0 := by
            by_contra hne
            obtain ⟨dp, rfl⟩ : ∃ dp, d = dp + 1 := ⟨d - 1, by omega⟩
            have hr := S.rside h2 dp rfl
            unfold RSideOut at hr
            exact absurd hr.1 (Nat.not_lt.2 hge)
          obtain ⟨hts0, hhv0, hF0⟩ := R.root hd0
          have hts₁' : s₁.tstack = [] := by rw [hnopush hpf, hts0]
          have hnil := List.append_eq_nil_iff.1 (hts₁'.symm.trans hts₁).symm
          have hts' : ((finishEdge v d (.tree e cls (.node c couts)) (new₁ ++ segsStack G.segs).length
              hv₁).run s₃).2.tstack = [] := by rw [hrts, hnil.1, hnil.2]; rfl
          have hhvF : ((finishEdge v d (.tree e cls (.node c couts)) (new₁ ++ segsStack G.segs).length
              hv₁).run s₃).1 = false := by rw [hr1, hhv, hhv0]
          have hg'' := hgR.trans hgeq₃
          exact {
            wf := by rw [hg'']; exact R.wf
            spec := by rw [hg'']; exact R.spec
            rooted := by rw [hg'']; exact R.rooted
            outs_v := R.outs_v
            sub := R.sub
            chain := ⟨by rw [hK'.svlo d (by omega)]; exact R.chain.1,
              fun k hk => by rw [hK'.svlo k (by omega)]; exact R.chain.2 k hk⟩
            rwalk := ⟨fun f hf => absurd hf (by rw [hF0]; exact List.not_mem_nil),
              ⟨fun t ht => absurd ht (by rw [hts']; exact List.not_mem_nil),
                by rw [hts']; exact List.Pairwise.nil⟩⟩
            fB := R.fB
            B_le := by rw [hts']; have := R.B_le; rw [hts0] at this; exact this
            hvB := fun hh => absurd hh (by rw [hhvF]; decide)
            root := fun _ => ⟨hts', hhvF, hF0⟩
            skel := by
              have := keepsR_finishBoundary (D := d + 1) v d (.tree e cls (.node c couts))
                (new₁ ++ segsStack G.segs).length hv₁ hge hend.shape
                (by rw [hgeq₃]; exact hgs ▸ h.v_lt) (by rw [hsgeq₃]; exact helt) hok (hR₂ h2).2.2
              rw [hgR]; exact this }
      · have hlt : (DfsOut.tree e cls (.node c couts)).cls.lowval d < d := Nat.lt_of_not_le hge
        have hlt' : cls.lowval d < d := hlt
        have hge' : ¬ d ≤ cls.lowval d := hge
        have hpr : push = (!hasVert && cls.isType1) := by rw [← hpushb]; simp [DfsOut.cls, hlt']
        have hhv : hv₁ = (hasVert || cls.isType1) := by
          rw [hhv₁, hpr]; cases hasVert <;> cases cls.isType1 <;> rfl
        have hxv : x = !s.stackDir[cls.lowval d]! := by rw [← hx]; simp [DfsOut.cls, hge']
        subst hxv
        have hpre' : hv₁ = false → new₁ = [] ∧
            ((refOuts G.g v d (DirsOf s d) (done.map (·.1)) false).1 ++
              if push = true then [⟨!s.stackDir[cls.lowval d]!, [vertItem v]⟩] else []) = [] := by
          intro h0
          rw [hhv, Bool.or_eq_false_iff] at h0
          have hpf : push = false := by simp [hpr, h0.1, h0.2]
          refine ⟨hnew₁ hpf h0.1, ?_⟩
          rw [hpf, (refOuts_hv_false _ _ (h.stPre.hv.trans h0.1)).2]; rfl
        obtain ⟨hr1, hrg, hdr', ⟨new', hrts, hR'⟩, hrI, hseg', hvs'⟩ :=
          stRet_finish (D := d + 1) hEf hend.inv hend.shape hD hg₃ hb₃ hlt hB hpre' hR₃ hRq₃ hI₃ hseg₃'
            hend.full
        rw [hdirs₃'] at hdr'
        rw [hxd₃, hsgeq₃, hit, hpr] at hR'
        simp only [↓reduceIte, Bool.not_not] at hR'
        refine ⟨?_, fun h2 => ?_⟩
        · subst hhv
          refine ⟨⟨new', hrts, ?_, fun h0 => by rw [hr1] at h0; cases h0⟩, ?_, hseg', ?_, ?_⟩
          · rw [hdr', refOuts_snoc, h.stPre.hv, refOut_ret_tree' hlt', DirsOf_getD s hlt']
            simpa only [List.append_assoc, DfsOut.e] using hR'
          · rw [hdr', refOuts_snoc, h.stPre.hv, refOut_ret_tree' hlt', hr1]
          · intro t ht
            rcases hvs' t ht with h' | ⟨-, -, h'⟩ | ⟨t₀, ht₀, h'⟩
            · exact Or.inl (h' ▸ List.mem_cons_self)
            · refine Or.inl (List.mem_cons_of_mem _ ?_)
              rw [h', DfsOut.vertsList_append]
              exact List.mem_append_right _
                (by simp [DfsOut.dest, DfsOut.vertsList, DfsTree.verts, DfsTree.v])
            · exact (hvs₃' t₀ ht₀).imp (h' ▸ id) fun ⟨t₁, ht₁, h₁⟩ => ⟨t₁, ht₁, h'.trans h₁⟩
          · rw [hdr', refOuts_snoc, h.stPre.hv, refOut_ret_tree' hlt', DirsOf_getD s hlt']
            simpa only [List.append_assoc, List.append_nil] using hrI
        · have R := h.r h2
          obtain ⟨lv, kind, ho, hl⟩ := ret_of_lowval_lt hlt
          have hk : kind ≠ .backEdge := fun hk => by
            subst hk
            have := hit; rw [ho] at this; cases this
          obtain ⟨hW₃, hK₃', hskel₃⟩ := hR₂ h2
          have hanc₃ : AncChain G.dfs v d s₃ :=
            ⟨by rw [hK₃.svlo d (by omega)]; exact R.chain.1,
              fun k hk => by rw [hK₃.svlo k (by omega)]; exact R.chain.2 k hk⟩
          have h2₃ : s₃.g.TwoConnected := by rw [hgeq₃]; exact h2
          have hsp₃ : G.dfs.Spec s₃.g := by rw [hgeq₃]; exact R.spec
          have hrt₃ : G.dfs.Rooted s₃.g := by rw [hgeq₃]; exact R.rooted
          obtain ⟨hW₁, hK₁, hn₁⟩ := hp.r h2
          have hK₁' := hK₁.1
          rw [hts₁] at hK₁' hn₁
          have hok := finishOk_of_guards ho hl hg₃ hEf rfl hend.inv hend.shape hD hb₃.v_lt hb₃.e_lt hb₃.q
            (hb₃.ends lv kind ho) hb₃.vert
          have hfront := finishEdge_frontier hb₃ hend.inv hend.shape hD
          have hhd : G.dfs.depth v = d := by
            have := (hanc₃.2 d (le_refl _)).2; rwa [hanc₃.1] at this
          have hclose := earFinish_close_len ho hl hk hEf rfl
          have hfr := (hW₃.frames (v, d, (new₁ ++ segsStack G.segs).length) List.mem_cons_self).2
          have hRf : s₃.RInvFront G.dfs v d (new₁ ++ segsStack G.segs).length := ⟨hfr.entries, hfr.disj⟩
          have hsh : FinishRShape G.dfs v d (.tree e cls (.node c couts)) (new₁ ++ segsStack G.segs).length
              hv₁ s₃ := hr₃ h2 (d - 1) (by omega)
          have hbot := finishEdge_bot_keep v d lv kind _ (new₁ ++ segsStack G.segs).length B hv₁ ho hl hEf
            rfl (fun h' => absurd h' hk) hK₁' hn₁
          have hBr : B ≤ ((finishEdge v d (.tree e cls (.node c couts)) (new₁ ++ segsStack G.segs).length
            hv₁).run s₃).2.tstack.length := hbot.1.1
          have hBr' : (new₁ ++ segsStack G.segs).length + (if hv₁ then 0 else 1) ≤
            ((finishEdge v d (.tree e cls (.node c couts)) (new₁ ++ segsStack G.segs).length
              hv₁).run s₃).2.tstack.length := hbot.2
          have hg'' := hgR.trans hgeq₃
          exact {
            wf := by rw [hg'']; exact R.wf
            spec := by rw [hg'']; exact R.spec
            rooted := by rw [hg'']; exact R.rooted
            outs_v := R.outs_v
            sub := R.sub
            chain := ⟨by rw [hK'.svlo d (by omega)]; exact R.chain.1,
              fun k hk => by rw [hK'.svlo k (by omega)]; exact R.chain.2 k hk⟩
            rwalk := ⟨fun f hf => ⟨by have := R.fB f hf; omega,
                finishEdge_rInvG_base f.1 f.2.1 v d lv kind _ _ f.2.2 hv₁ ho hl hb₃.v_lt hend.inv
                  hend.shape hok hfront (fun _ => hsh) (fun h' => absurd h' hk) (fun _ => hclose)
                  (by have := R.fB f hf; omega)
                  (fun h' => by have := R.fB f hf; have := hn₁ h'; omega)
                  (hW₃.frames f (List.mem_cons_of_mem _ hf)).2⟩,
              finishEdge_rInvTop v d lv kind _ _ hv₁ ho hl hb₃.v_lt hend.inv hend.shape hok hg₃ hfront
                (fun _ => hsh) (fun h' => absurd h' hk) h2₃ hsp₃ hrt₃ hhd hanc₃.1 hanc₃.2 hRf⟩
            fB := R.fB
            B_le := hBr
            hvB := fun _ => by
              cases hv₁
              · simp only [Bool.false_eq_true, ↓reduceIte] at hBr'; omega
              · simp only [↓reduceIte] at hBr'; have := hn₁ rfl; omega
            root := fun hd0 => absurd hlt (by rw [hd0]; exact Nat.not_lt_zero _)
            skel := by
              have := keepsR_finishEdge_site (dfs := G.dfs) (D := d + 1) v d lv kind _ _ hv₁ ho hl
                hb₃.v_lt hend.inv hend.shape hok hD hfront (fun _ => hsh) h2₃ hsp₃ hrt₃ hhd
                hanc₃.1 hanc₃.2 hRf hEf.q_root hEf.sv_child hcb₃ hskel₃
              rw [hgR]; exact this }
    obtain ⟨hSt, hRR⟩ := key
    refine ⟨fun h0 => hhv' (by rw [hhv₁, h0]; rfl), hK', hvF, ?_⟩
    have hhf : ((finishEdge v d (.tree e cls (.node c couts)) (new₁ ++ segsStack G.segs).length hv₁).run
        s₃).1 = false → hasVert = false := fun h0 => by
      cases hhV : hasVert
      · rfl
      · exact absurd (hhv' (by rw [hhv₁, hhV]; rfl)) (by rw [h0]; decide)
    have hnb : n + (DfsTree.node c couts).edgePostorder.length + 1 =
        n + (DfsOut.tree e cls (.node c couts)).block.length := by
      simp [DfsOut.block, Nat.add_assoc]
    rw [hnb] at hrg' ho'
    exact {
      full := hF'.monoP (pushed_tree_iff hhv')
      d_anc := h.d_anc
      split := by rw [List.map_append, List.append_assoc]; exact h.split
      sorted := h.sorted
      wf := h.wf
      ends := h.ends
      nodup := h.nodup
      anc_lt := h.anc_lt
      v_lt := h.v_lt
      w_lt := h.w_lt
      e_lt := h.e_lt
      enodup := h.enodup
      comp := h.comp
      pe_anc := h.pe_anc
      sv_size := hK'.sv.trans h.sv_size
      sd_size := hK'.sd.trans h.sd_size
      anc_sv := fun k hk => by rw [hK'.svlo k (by omega)]; exact h.anc_sv k hk
      sv_d := by rw [hK'.svlo d (by omega)]; exact h.sv_d
      ear := hC'
      inv := hi'
      shape := hs'
      ranges := hrg'
      σ_nodup := h.σ_nodup
      σ_lt := by rw [hgs', hgeq₃]; exact h.σ_lt
      post := by have := h.post; rw [DfsOut.edgePostorderList_cons] at this; exact this.right
      path := by
        have := h.path; rw [DfsOut.edgePostorderList_cons, List.length_append, ← Nat.add_assoc] at this
        exact this
      close := hc'
      canon := hcn'
      tern := htn'.trans hend.tern
      P_past := pushed_past (fun e he hp => by
          rcases hp with hp | ⟨_, hp⟩
          · exact h.P_past e he hp
          · exact absurd hp.symm (vertItem_ne_edgeItem h.v_lt e))
        h.w_lt_o (DfsOut.edges_subset_block _) h.post_o h.σ_nodup
      owned := ho'.mono fun i hi => (pushed_tree_iff hhv' i).1 (Or.inl hi)
      sts_len := h.sts_len
      origs_len := h.origs_len
      sts_le := h.sts_le
      sts_n := Nat.le_trans h.sts_n (Nat.le_add_right _ _)
      origs_le := h.origs_le
      Pv := by
        rintro w hw ((hp | ⟨_, hp⟩) | ⟨w', hw', hp⟩ | ⟨e, _, hp⟩)
        · exact h.Pv_rest w hw hp
        · exact hvr (WalkM.vertItem_inj hp ▸ hw)
        · exact hdisjV w' hw' (WalkM.vertItem_inj hp ▸ hw)
        · exact vertItem_ne_edgeItem (h.w_lt_rest w hw) e hp
      Pe := by
        rintro e he ((hp | ⟨_, hp⟩) | ⟨w', hw', hp⟩ | ⟨e', he', hp⟩)
        · exact h.Pe_rest e he hp
        · exact vertItem_ne_edgeItem h.v_lt e hp.symm
        · exact vertItem_ne_edgeItem (h.w_lt_o w' hw') e hp.symm
        · exact hdisjE e' he' (edgeItem_inj hp ▸ he)
      Pcur := by
        rintro h0 ((hp | ⟨h1, _⟩) | ⟨w', hw', hp⟩ | ⟨e, _, hp⟩)
        · exact h.Pcur (hhf h0) hp
        · rw [h0] at h1; cases h1
        · exact hvo (WalkM.vertItem_inj hp ▸ hw')
        · exact vertItem_ne_edgeItem h.v_lt e hp
      Pcur' := fun h1 => Or.inl (Or.inr ⟨h1, rfl⟩)
      vcover := hvc'
      r := fun h2 => hRR (by rw [hK'.g] at h2; exact h2)
      d_fs := h.d_fs
      height := by rw [hK'.sd]; exact h.height
      base_out := h.base_out
      stPre := by rw [List.map_append]; exact hSt
      segs_len := h.segs_len
      live_cur := by rw [List.map_append]; exact hlc
      live_lower := hll }

variable {t : DfsTree}

/-- **The backbone induction.** -/
theorem backbone :
    (∀ (t : DfsTree) (d : Nat), BbTree t d) ∧
    (∀ (v d : Nat) (outs : List DfsOut) (hasVert : Bool), BbOuts v d outs hasVert) ∧
    (∀ (v d : Nat) (o : DfsOut) (hasVert : Bool), BbOut v d o hasVert) :=
  walkTree.mutual_induct _ _ _
    (fun d v outs ih => bbTree_node v outs d ih)
    bbOut_step
    bbOuts_nil
    bbOuts_cons

/-- **The backbone theorem**: the conjunction at the entry gives the conjunction at the exit. -/
theorem walkTree_inv (h : WalkInv G t d s) :
    wp (walkTree t d) (fun _ s' => WalkInvEnd G t d s s') s :=
  backbone.1 t d G s h

theorem walkOuts_inv (h : WalkInvOut G B v d outs₀ done rest hasVert n P s) :
    wp (walkOuts v d rest hasVert) (fun hv' s' => (hasVert = true → hv' = true) ∧ Keep (d + 1) 0 s s' ∧
      ∃ done', WalkInvOut G B v d outs₀ done' [] hv' (n + (DfsOut.edgePostorderList rest).length)
        (Pushed G.g (fun i => P i ∨ (hv' = true ∧ i = vertItem v))
          (DfsOut.vertsList rest) (DfsOut.edgesList rest)) s') s :=
  backbone.2.1 v d rest hasVert G B outs₀ done n P s h

end WalkState
end Spqr
