import Spqr.RangesCloseSites

/-!
# Close invariant through the walk

`CloseInv` is threaded through the mutual walk induction (the shape of `rgTree`/`rgOuts`/`rgOut`):
each `finishEdge` is discharged by `finishEdge_closeInv` from a `CloseCtx`, whose range/ear parts
come from the enclosing induction and whose per-site parts (`CloseSite`) are supplied as
hypotheses in the `wp` shape of `RgTree` (`CsTree`/`CsOuts`/`CsOut`).
-/

namespace Spqr
open WalkM
namespace WalkState
variable {s : WalkState} {σ : List Nat} {n D : Nat}

/-- The per-site facts of `CloseCtx` not provided by the range/ear layers. -/
structure CloseSite (σ : List Nat) (n curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool)
    (s : WalkState) : Prop where
  pos : σ[n]? = some o.e
  block : o.block <:+: σ
  ends : o.Ends s.g curV
  path : ∀ k, k < d → ∃ e, e < s.g.ne ∧ n < σ.idxOf e ∧
    Items.PairEq (s.stackVerts[k]!, s.stackVerts[k + 1]!) s.g.edges[e]!
  dest_edge : o.cls.isTree = true → o.cls.lowval d < d →
    ∃ e, e < s.g.ne ∧ e ≠ o.e ∧ s.g.Inc e o.dest
  dest_lt : o.dest < s.g.nv
  bd_loop : o.cls.isTree = false → d ≤ o.cls.lowval d → o.dest = curV
  bd_vert : o.cls.isTree = true → d ≤ o.cls.lowval d →
    ∀ t ∈ (if o.cls.lowval d = d + 1 then s.tstack.head? else s.tstack.tail.head?),
      t.spans.2 = [vertItem o.dest]
  bd_node : o.cls.isTree = true → d ≤ o.cls.lowval d → o.cls.lowval d ≠ d + 1 →
    ∀ b ∈ s.tstack.head?, ∃ c, b.spans.1 = [c] ∧ Items.type s.items c ∉ [NodeType.F, .V] ∧
      (Items.type s.items c = .Q → Items.ch s.items c = []) ∧
      ∃ a b', Items.vs s.items c = (some a, some b') ∧ Items.PairEq (a, b') (curV, o.dest)
  p_site : o.cls.lowval d < d → o.cls.isType1 = true →
    result (condP curV (o.cls.lowval d) true) (feRest curV d o origTstack hasVert s) = true →
    PSite curV (o.cls.lowval d) (feRest curV d o origTstack hasVert s)
  v_site : o.cls.isTree = true → o.cls.lowval d < d → hasVert = true → o.cls.isType1 = true →
    ∃ t, VSite curV ((maybeUnwrapNxt (if feSingle d o s then NodeType.S else .R)).run (feS₂ d o s)).1 t
      (cvS₅ curV s.stackDir[d]! true origTstack (feSingle d o s) (feS₂ d o s))
  l1_site : o.cls.isTree = true → o.cls.lowval d < d → ∀ k,
    (∀ j, j ≤ k → result (loop1Cond d) (l1Iter d o s j) = true) →
    Shape (l1S₁ d s.stackDir[d]! (l1Iter d o s k)) ∧
    (∃ a b rest, (l1S₁ d s.stackDir[d]! (l1Iter d o s k)).tstack = a :: b :: rest) ∧
    ∃ t, VSite t.vStart
      (result (maybeUnwrapNxt (l1Ty d s.stackDir[d]! (l1Iter d o s k))) (l1S₁ d s.stackDir[d]! (l1Iter d o s k)))
      t (after mergeTstackTops (l1S₂ d s.stackDir[d]! (l1Iter d o s k)))

theorem CloseCtx.of_site {curV d : Nat} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}
    (hD : D = if o.cls.isTree then d + 1 else d) (hrs : RgS σ n D s) (hnd : σ.Nodup)
    (hg : FinishGuards d o origTstack hasVert s) (hb : FinishBook curV d o origTstack hasVert s)
    (hf : Frontier (o := o) d origTstack s) (hr : FinishR σ n curV d o origTstack hasVert s)
    (hcl : s.CloseInv) (hc : CloseSite σ n curV d o origTstack hasVert s) :
    CloseCtx σ n D curV d o origTstack hasVert s where
  hD := hD
  nodup := hnd
  lt := hrs.2.2
  pos := hc.pos
  block := hc.block
  ranges := hrs.1
  shape := hrs.2.1
  guards := hg
  book := hb
  frontier := hf
  finishR := hr
  close := hcl
  ends := hc.ends
  path := hc.path
  dest_edge := hc.dest_edge
  dest_lt := hc.dest_lt
  bd_loop := hc.bd_loop
  bd_vert := hc.bd_vert
  bd_node := hc.bd_node
  p_site := hc.p_site
  v_site := hc.v_site
  l1_site := hc.l1_site

mutual
def CsTree (σ : List Nat) (n : Nat) (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  match t with
  | .node v outs => CsOuts σ n v d outs false { s with stackVerts := s.stackVerts.set! d v }

def CsOuts (σ : List Nat) (n v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  match outs with
  | [] => True
  | o :: rest => CsOut σ n v d o hasVert s ∧
      wp (walkOut v d o hasVert) (fun hasVert' s' => CsOuts σ (n + o.block.length) v d rest hasVert' s') s

def CsOut (σ : List Nat) (n v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  wp (walkOutPre v d o hasVert) (fun hasVert' s₁ =>
    match o with
    | .tree _ _ child =>
      wp (modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }) (fun _ s₂ =>
        CsTree σ n child (d + 1) s₂ ∧
        wp (walkTree child (d + 1)) (fun _ s₃ =>
          CloseSite σ (n + child.edgePostorder.length) v d o s₁.tstack.length hasVert' s₃) s₂) s₁
    | .back .. => CloseSite σ n v d o s₁.tstack.length hasVert' s₁) s
end

theorem walkOutPre_closeInv' (h : s.CloseInv) (hs : Shape s) {v : Nat} (d : Nat)
    (o : DfsOut) {hasVert : Bool} (hb : VertBook v hasVert s) :
    wp (walkOutPre v d o hasVert) (fun _ s' => s'.CloseInv) s := by
  have h' : ({ s with stackDir := s.stackDir.set! d (if o.cls.lowval d ≥ d then false else !s.stackDir[o.cls.lowval d]!) } : WalkState).CloseInv :=
    h.frame rfl rfl (fun _ hi => hi)
  unfold walkOutPre
  simp only [wp_bind, wp_stackDir, wp_setStackDir, wp_ite]
  split
  · rename_i hc
    have hf : hasVert = false := by cases hasVert <;> simp_all
    obtain ⟨hv, -, -⟩ := hb hf
    simp only [wp_bind, wp_pure]
    exact h'.pushVert v d hv (by have := hs.size; show s.g.nv < s.items.size; omega)
  · simp only [wp_pure]; exact h'

abbrev CcInvTree (σ : List Nat) (n : Nat) (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  (∀ v outs, t = .node v outs → ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).RangesInv σ n d) →
  Shape s → σ.Nodup → (∀ e ∈ σ, e < s.g.ne) → GuardsTree t d s → BookTree t d s → FrontiersTree t d s →
  RgTree σ n t d s → CsTree σ n t d s → s.CloseInv →
  wp (walkTree t d) (fun _ s' => s'.CloseInv) s

abbrev CcInvOuts (σ : List Nat) (n v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  RgS σ n d s → σ.Nodup → GuardsOuts v d outs hasVert s → BookOuts v d outs hasVert s →
  FrontiersOuts v d outs hasVert s → RgOuts σ n v d outs hasVert s → CsOuts σ n v d outs hasVert s →
  s.CloseInv → wp (walkOuts v d outs hasVert) (fun _ s' => s'.CloseInv) s

abbrev CcInvOut (σ : List Nat) (n v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  RgS σ n d s → σ.Nodup → GuardsOut v d o hasVert s → BookOut v d o hasVert s →
  FrontiersOut v d o hasVert s → RgOut σ n v d o hasVert s → CsOut σ n v d o hasVert s →
  s.CloseInv → wp (walkOut v d o hasVert) (fun _ s' => s'.CloseInv) s

mutual
theorem ccTree : ∀ (σ : List Nat) (n : Nat) (t : DfsTree) (d : Nat) (s : WalkState), CcInvTree σ n t d s
  | σ, n, .node v outs, d, s => fun hi hs hnd hσ hg hb hf hr hc hcl => by
    unfold walkTree
    simp only [wp_bind, wp_modify]
    unfold GuardsTree at hg; unfold BookTree at hb; unfold FrontiersTree at hf; unfold RgTree at hr
    unfold CsTree at hc
    refine wp_imp (wp_imp (wp_of_forall fun hv s' ⟨⟨hi', hs', hσ'⟩, hvb, hpr⟩ hcl' => ?_)
      (rgOuts σ n v d outs false _ ⟨hi v outs rfl, hs.frame', hσ⟩ hnd hg hb hr))
      (ccOuts σ n v d outs false _ ⟨hi v outs rfl, hs.frame', hσ⟩ hnd hg hb hf hr hc
        (hcl.frame rfl rfl (fun _ h => h)))
    cases hv
    · simp only [Bool.false_eq_true, ↓reduceIte, wp_bind, wp_setStackDir]
      obtain ⟨hv, -, -⟩ := hvb rfl
      exact (hcl'.frame (s' := { s' with stackDir := s'.stackDir.set! d true }) rfl rfl (fun _ h => h)).pushVert
        v d hv (by have := hs'.size; show s'.g.nv < s'.items.size; omega)
    · exact hcl'

theorem ccOuts : ∀ (σ : List Nat) (n v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState),
    CcInvOuts σ n v d outs hasVert s
  | σ, n, v, d, [], hasVert, s => fun _ _ _ _ _ _ _ hcl => by
    unfold walkOuts
    simp only [wp_pure]
    exact hcl
  | σ, n, v, d, o :: rest, hasVert, s => fun hrs hnd hg hb hf hr hc hcl => by
    unfold GuardsOuts at hg; unfold BookOuts at hb; unfold FrontiersOuts at hf; unfold RgOuts at hr
    unfold CsOuts at hc
    unfold walkOuts
    simp only [wp_bind]
    exact wp_imp (wp_imp (wp_imp (wp_imp (wp_imp (wp_imp (wp_imp (wp_of_forall
      fun hv' s' hrs' hg' hb' hf' hr' hc' hcl' =>
        ccOuts σ (n + o.block.length) v d rest hv' s' hrs' hnd hg' hb' hf' hr' hc' hcl')
      (rgOut σ n v d o hasVert s hrs hnd hg.1 hb.1 hr.1)) hg.2) hb.2) hf.2) hr.2) hc.2)
      (ccOut σ n v d o hasVert s hrs hnd hg.1 hb.1 hf.1 hr.1 hc.1 hcl)

theorem ccOut : ∀ (σ : List Nat) (n v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState), CcInvOut σ n v d o hasVert s
  | σ, n, v, d, o, hasVert, s => fun ⟨hi, hs, hσ⟩ hnd hg hb hf hr hc hcl => by
    rw [walkOut_eq, wp_bind]
    unfold GuardsOut at hg; unfold BookOut at hb; unfold FrontiersOut at hf; unfold RgOut at hr
    unfold CsOut at hc
    refine wp_imp (wp_of_forall fun hv' s₁ ⟨⟨hi₁, hs₁, hσ₁⟩, hg₁, hb₁, hf₁, hr₁, hc₁, hcl₁⟩ => ?_)
      (wp_and (walkOutPre_ranges hi hs hσ hb.1 hr.1) (wp_and hg (wp_and hb.2 (wp_and hf (wp_and hr.2
        (wp_and hc (walkOutPre_closeInv' hcl hs d o hb.1)))))))
    unfold walkOutRest
    rw [wp_bind, wp_tstackSize]
    cases o with
    | back e cls dest =>
      try simp only [wp_pure] at hg₁ hb₁ hf₁ hr₁ hc₁
      simp only [wp_bind, wp_pure]
      exact finishEdge_closeInv (CloseCtx.of_site (D := d)
        (by rw [if_neg fun h => by obtain ⟨_, _, _, h⟩ := hb₁.tree.1 h; cases h])
        ⟨hi₁, hs₁, hσ₁⟩ hnd hg₁ hb₁ hf₁ hr₁ hcl₁ hc₁)
    | tree e cls child =>
      try simp only [wp_bind, wp_modify] at hg₁ hb₁ hf₁ hr₁ hc₁
      simp only [wp_bind, wp_modify]
      have pre : ∀ w outs, child = .node w outs →
          ({ { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne } with
            stackVerts := s₁.stackVerts.set! (d + 1) w } : WalkState).RangesInv σ n (d + 1) :=
        fun w outs _ =>
          (hi₁.frame' (s' := { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne })).setSv w
      refine wp_imp (wp_imp (wp_imp (wp_of_forall fun _ s₃ ⟨hg₃, hb₃, hf₃, hr₃, hc₃⟩ ⟨hi₃, hs₃, hσ₃⟩ hcl₃ => ?_)
        (wp_and hg₁.2 (wp_and hb₁.2 (wp_and hf₁.2 (wp_and hr₁.2 hc₁.2)))))
        (rgTree σ n child (d + 1) _ pre hs₁.frame' hnd hσ₁ hg₁.1 hb₁.1 hr₁.1))
        (ccTree σ n child (d + 1) _ pre hs₁.frame' hnd hσ₁ hg₁.1 hb₁.1 hf₁.1 hr₁.1 hc₁.1
          (hcl₁.frame rfl rfl (fun _ h => h)))
      exact finishEdge_closeInv (CloseCtx.of_site (D := d + 1)
        (by rw [if_pos (hb₃.tree.2 ⟨_, _, _, rfl⟩)]) ⟨hi₃, hs₃, hσ₃⟩ hnd hg₃ hb₃ hf₃ hr₃ hcl₃ hc₃)
end

/-- `walkTree` preserves `CloseInv`, given the range/ear hypotheses of `walkTree_rangesInv`, the
frontiers, and the per-site facts `CsTree`. -/
theorem walkTree_closeInv (t : DfsTree) (d : Nat) (s : WalkState)
    (hi : ∀ v outs, t = .node v outs → ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).RangesInv σ n d)
    (hs : Shape s) (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne) (hg : GuardsTree t d s) (hb : BookTree t d s)
    (hf : FrontiersTree t d s) (hr : RgTree σ n t d s) (hc : CsTree σ n t d s) (hcl : s.CloseInv) :
    ((walkTree t d).run s).2.CloseInv :=
  ccTree σ n t d s hi hs hnd hσ hg hb hf hr hc hcl

/-- The per-site facts for a forest, in the shape of `RootsCover`. -/
def CsForest (σ : List Nat) (n : Nat) : List DfsTree → WalkState → Prop
  | [], _ => True
  | t :: rest, s => CsTree σ n t 0 s ∧
      wp (walkTree t 0) (fun _ s₁ =>
        wp (popTstack >>= fun top => modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 })
          (fun _ s₂ => CsForest σ (n + t.edgePostorder.length) rest s₂) s₁) s

theorem forest_closeInv_of_sites {g : Graph} {σ : List Nat} (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < g.ne) :
    ∀ forest pre n s, RootState g pre s → s.RangesInv σ n 0 → ForestOK g (pre ++ forest) →
      (∀ t ∈ forest, t.WF []) → (∀ t ∈ forest, t.Ends g) → RootsCover σ n forest s →
      PostAt σ n (edgePostorderForest forest) → CsForest σ n forest s → s.CloseInv →
      wp (walkForest forest) (fun _ s' => s'.CloseInv) s
  | [], _, _, _, _, _, _, _, _, _, _, _, hcl => by simpa [walkForest, wp_pure] using hcl
  | t :: rest, pre, n, s, h, hr, hf, hwf, hends, hc, hat, hcs, hcl => by
    have hb := h.book hf (hwf t (by simp)) (hends t (by simp))
    have hg := gbTree t 0 s hb
    have hi : ∀ v outs, t = .node v outs →
        ({ s with stackVerts := s.stackVerts.set! 0 v } : WalkState).RangesInv σ n 0 :=
      fun v _ _ => hr.stackVerts_of_nil h.tstack _
    have hσ' : ∀ e ∈ σ, e < s.g.ne := by rwa [h.g_eq]
    have hfront := walkTree_frontiers t 0 s (fun v outs ht => (hi v outs ht).inv) h.shape hg hb
    have hat' : PostAt σ n (t.edgePostorder ++ edgePostorderForest rest) := hat
    have hsched := scheduleTree σ n t 0 s hi h.shape hnd hσ' hg hb hfront hc.1 hat'.left
    have hrg := rgTree σ n t 0 s hi h.shape hnd hσ' hg hb hsched
    have hcc := ccTree σ n t 0 s hi h.shape hnd hσ' hg hb hfront hsched hcs.1 hcl
    have hnv : 0 < s.g.nv := by
      obtain ⟨v, outs⟩ := t
      have hv := RootState.hvlt hf v (by simp [DfsTree.verts])
      rw [h.g_eq]; omega
    have hk : wp (walkTree t 0) (fun _ s' => RootOK s') s :=
      walkTree_rootOK t s (fun v outs ht => (hi v outs ht).inv) h.shape hg hb h.tstack
      (by rw [h.sd, ← h.g_eq]; exact hnv)
    have hp := (walk_place_aux g).1 t 0 _ _ s h.place (RootState.hvlt hf) (RootState.helt hf)
      (RootState.hvn hf).1 (RootState.hen hf).1 (RootState.hPv hf) (RootState.hPe hf)
    have hst := h.step hf (hwf t (by simp)) (hends t (by simp))
    show wp ((walkTree t 0 >>= fun _ => popTstack >>= fun top =>
      modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 }) >>= fun _ => walkForest rest) _ s
    simp only [wp_bind]
    refine wp_imp (wp_of_forall fun _ s₁ ⟨hrs, hk, hp, hst, hc, hcs, hcl⟩ => ?_)
      (wp_and hrg (wp_and hk (wp_and hp (wp_and hst (wp_and hc.2 (wp_and hcs.2 hcc))))))
    have hrpop := hrs.1.root_append hrs.2.1 hk (noParent_of_cnt_eq_zero hp.root)
    have hp' : s₁.Place s₁.g (Pushed g (Pushed g (fun _ => False) (pre.flatMap DfsTree.verts)
        (pre.flatMap DfsTree.edges)) t.verts t.edges) (fun _ => False) := by
      rw [hp.g_eq]; exact hp
    have hclpop := rootAppend_closeInv hcl hp' hk
    exact wp_mono (popTstack >>= fun top => modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 })
      (wp_and hrpop (wp_and hst (wp_and hc (wp_and hcs hclpop)))) fun _ s₂ ⟨hr₂, hs₂, hc₂, hcs₂, hcl₂⟩ =>
      forest_closeInv_of_sites hnd hσ rest (pre ++ [t]) (n + t.edgePostorder.length) s₂ hs₂ hr₂
        (by simpa using hf) (fun t' ht' => hwf t' (by simp [ht']))
        (fun t' ht' => hends t' (by simp [ht'])) hc₂ hat'.right hcs₂ hcl₂

/-- `CloseInv` for the walk of a forest, from `RootsCover` (the range-side schedule) and the
per-site facts `CsForest`. -/
theorem walk_closeInv_of_sites (g : Graph) (tern : Bool) (forest : List DfsTree)
    (hf : ForestOK g forest) (hwf : ∀ t ∈ forest, t.WF []) (hends : ∀ t ∈ forest, t.Ends g)
    (hc : RootsCover (edgePostorderForest forest) 0 forest (WalkState.init g tern))
    (hcs : CsForest (edgePostorderForest forest) 0 forest (WalkState.init g tern)) :
    (g.walk tern forest).CloseInv := by
  have hnd := DfsData.edgePostorderForest_perm.nodup_iff.2 hf.edges_nodup
  have hσ : ∀ e ∈ edgePostorderForest forest, e < g.ne :=
    fun e he => hf.edges_lt e (DfsData.edgePostorderForest_perm.subset he)
  have h := forest_closeInv_of_sites hnd hσ forest [] 0 (WalkState.init g tern)
    (rootState_init g tern) (init_rangesInv g tern hnd) (by simpa using hf) hwf hends hc
    ⟨[], [], rfl, by simp⟩ hcs (init_closeInv g tern)
  simpa only [wp, Graph.walk] using h

end WalkState
end Spqr
