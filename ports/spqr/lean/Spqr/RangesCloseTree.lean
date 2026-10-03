import Spqr.RangesCloseContent

/-!
# Close invariant through the walk

`CloseInv` is threaded through the mutual walk induction (the shape of `rgTree`/`rgOuts`/`rgOut`):
each `finishEdge` is discharged by `finishEdge_closeInv` from a `CloseCtx`, whose range/ear parts
come from the enclosing induction, whose DFS-side parts (`DfsSite`) are supplied in the `wp` shape of
`RgTree` (`CsTree`/`CsOuts`/`CsOut`) and discharged on the real walk by `dsTree`, and whose remaining
site fields are the named admissions `closeCtx_bd_vert`/`bd_node`/`p_site`/`v_site`/`l1_site`.
-/

namespace Spqr
open WalkM
namespace WalkState
variable {s : WalkState} {σ : List Nat} {n D : Nat}

/-- The DFS-side facts of `CloseCtx` at a `finishEdge` call: the site's position in the schedule,
the edge's endpoints, the open path, and the child's second edge; all follow from `DfsTree.WF`/`Ends`
and the frame (`dsTree`). -/
structure DfsSite (σ : List Nat) (n curV d : Nat) (o : DfsOut) (s : WalkState) : Prop where
  pos : σ[n]? = some o.e
  block : o.block <:+: σ
  ends : o.Ends s.g curV
  path : ∀ k, k < d → ∃ e, e < s.g.ne ∧ n < σ.idxOf e ∧
    Items.PairEq (s.stackVerts[k]!, s.stackVerts[k + 1]!) s.g.edges[e]!
  dest_edge : o.cls.isTree = true → o.cls.lowval d < d →
    ∃ e, e < s.g.ne ∧ e ≠ o.e ∧ s.g.Inc e o.dest
  dest_edges : o.cls.isTree = true → ∀ e, e < s.g.ne → s.g.Inc e o.dest → subEdges o e
  dest_lt : o.dest < s.g.nv
  bd_loop : o.cls.isTree = false → d ≤ o.cls.lowval d → o.dest = curV

/-- Everything the walk induction exports at a `finishEdge` call: the range/ear/frontier/R-side
postconditions (`RgTree`, `gbTree`, `BookTree`, `FrontiersTree`, `scheduleTree`), the close invariant,
and the DFS-side site facts. -/
structure CloseBase (σ : List Nat) (n D curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool)
    (s : WalkState) : Prop where
  hD : D = if o.cls.isTree then d + 1 else d
  nodup : σ.Nodup
  rgs : RgS σ n D s
  guards : FinishGuards d o origTstack hasVert s
  book : FinishBook curV d o origTstack hasVert s
  frontier : Frontier (o := o) d origTstack s
  finishR : FinishR σ n curV d o origTstack hasVert s
  close : s.CloseInv
  site : DfsSite σ n curV d o s

section Admissions
variable {curV d : Nat} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}

/-- Block entry of a boundary edge carries exactly `vertItem o.dest` on its active side
(checker: `closeCtx_bd_vert`, kind `ctx_bd_vert`). -/
theorem closeCtx_bd_vert (h : CloseBase σ n D curV d o origTstack hasVert s)
    (hc : CloseContent curV d o origTstack hasVert s) :
    o.cls.isTree = true → d ≤ o.cls.lowval d →
    ∀ t ∈ (if o.cls.lowval d = d + 1 then s.tstack.head? else s.tstack.tail.head?),
      t.spans.2 = [vertItem o.dest] := by
  sorry

/-- Completed block: the top entry's passive side is one non-`F`/`V` node (a leaf if `Q`) with
terminals `{curV, o.dest}` (checker: `closeCtx_bd_node`, kind `ctx_bd_node`). -/
theorem closeCtx_bd_node (h : CloseBase σ n D curV d o origTstack hasVert s)
    (hc : CloseContent curV d o origTstack hasVert s) :
    o.cls.isTree = true → d ≤ o.cls.lowval d → o.cls.lowval d ≠ d + 1 →
    ∀ b ∈ s.tstack.head?, ∃ c, b.spans.1 = [c] ∧ Items.type s.items c ∉ [NodeType.F, .V] ∧
      (Items.type s.items c = .Q → Items.ch s.items c = []) ∧
      ∃ a b', Items.vs s.items c = (some a, some b') ∧ Items.PairEq (a, b') (curV, o.dest) := by
  sorry

/-- The P-merge site (checker: `closeCtx_p_site`, kinds `psite_*`). -/
theorem closeCtx_p_site (h : CloseBase σ n D curV d o origTstack hasVert s)
    (hc : CloseContent curV d o origTstack hasVert s) :
    o.cls.lowval d < d → o.cls.isType1 = true →
    result (condP curV (o.cls.lowval d) true) (feRest curV d o origTstack hasVert s) = true →
    PSite curV (o.cls.lowval d) (feRest curV d o origTstack hasVert s) := by
  sorry

/-- The type-1 vertex-close site (checker: `closeCtx_v_site`, kinds `vsite_*`). -/
theorem closeCtx_v_site (h : CloseBase σ n D curV d o origTstack hasVert s)
    (hc : CloseContent curV d o origTstack hasVert s) :
    o.cls.isTree = true → o.cls.lowval d < d → hasVert = true → o.cls.isType1 = true →
    ∃ t, VSite curV ((maybeUnwrapNxt (if feSingle d o s then NodeType.S else .R)).run (feS₂ d o s)).1 t
      (cvS₅ curV s.stackDir[d]! true origTstack (feSingle d o s) (feS₂ d o s)) := by
  sorry

/-- The loop-1 iteration sites (checker: `closeCtx_l1_site`, kinds `l1site_*`). -/
theorem closeCtx_l1_site (h : CloseBase σ n D curV d o origTstack hasVert s)
    (hc : CloseContent curV d o origTstack hasVert s) :
    o.cls.isTree = true → o.cls.lowval d < d → ∀ k,
    (∀ j, j ≤ k → result (loop1Cond d) (l1Iter d o s j) = true) →
    Shape (l1S₁ d s.stackDir[d]! (l1Iter d o s k)) ∧
    (∃ a b rest, (l1S₁ d s.stackDir[d]! (l1Iter d o s k)).tstack = a :: b :: rest) ∧
    ∃ t, VSite t.vStart
      (result (maybeUnwrapNxt (l1Ty d s.stackDir[d]! (l1Iter d o s k))) (l1S₁ d s.stackDir[d]! (l1Iter d o s k)))
      t (after mergeTstackTops (l1S₂ d s.stackDir[d]! (l1Iter d o s k))) := by
  sorry

/-- `CloseCtx` from the exports and the block contents. -/
theorem CloseCtx.of_exports (h : CloseBase σ n D curV d o origTstack hasVert s)
    (hc : CloseContent curV d o origTstack hasVert s) :
    CloseCtx σ n D curV d o origTstack hasVert s where
  hD := h.hD
  nodup := h.nodup
  lt := h.rgs.2.2
  pos := h.site.pos
  block := h.site.block
  ranges := h.rgs.1
  shape := h.rgs.2.1
  guards := h.guards
  book := h.book
  frontier := h.frontier
  finishR := h.finishR
  close := h.close
  ends := h.site.ends
  path := h.site.path
  dest_edge := h.site.dest_edge
  dest_edges := h.site.dest_edges
  dest_lt := h.site.dest_lt
  bd_loop := h.site.bd_loop
  bd_vert := closeCtx_bd_vert h hc
  bd_node := closeCtx_bd_node h hc
  p_site := closeCtx_p_site h hc
  v_site := closeCtx_v_site h hc
  l1_site := closeCtx_l1_site h hc

/-- The block contents at every `finishEdge` pre-state; supplied by the walk induction
(`WalkBackbone.lean`, `WalkInvEnd`), a named hypothesis here. -/
theorem closeBase_content (h : CloseBase σ n D curV d o origTstack hasVert s) :
    CloseContent curV d o origTstack hasVert s := by
  sorry

theorem CloseCtx.of_base (h : CloseBase σ n D curV d o origTstack hasVert s) :
    CloseCtx σ n D curV d o origTstack hasVert s :=
  CloseCtx.of_exports h (closeBase_content h)

end Admissions

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
          DfsSite σ (n + child.edgePostorder.length) v d o s₃) s₂) s₁
    | .back .. => DfsSite σ n v d o s₁) s
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
      exact finishEdge_closeInv (CloseCtx.of_base (D := d)
        ⟨by rw [if_neg fun h => by obtain ⟨_, _, _, h⟩ := hb₁.tree.1 h; cases h],
          hnd, ⟨hi₁, hs₁, hσ₁⟩, hg₁, hb₁, hf₁, hr₁, hcl₁, hc₁⟩)
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
      exact finishEdge_closeInv (CloseCtx.of_base (D := d + 1)
        ⟨by rw [if_pos (hb₃.tree.2 ⟨_, _, _, rfl⟩)], hnd, ⟨hi₃, hs₃, hσ₃⟩, hg₃, hb₃, hf₃, hr₃, hcl₃, hc₃⟩)
end

/-- `walkTree` preserves `CloseInv`, given the range/ear hypotheses of `walkTree_rangesInv`, the
frontiers, and the per-site facts `CsTree`. -/
theorem walkTree_closeInv (t : DfsTree) (d : Nat) (s : WalkState)
    (hi : ∀ v outs, t = .node v outs → ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).RangesInv σ n d)
    (hs : Shape s) (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne) (hg : GuardsTree t d s) (hb : BookTree t d s)
    (hf : FrontiersTree t d s) (hr : RgTree σ n t d s) (hc : CsTree σ n t d s) (hcl : s.CloseInv) :
    ((walkTree t d).run s).2.CloseInv :=
  ccTree σ n t d s hi hs hnd hσ hg hb hf hr hc hcl

/-! ### The site facts on the real walk -/

/-- The open path `sv` (`sv[k]` at depth `k`) is joined by edges scheduled at or after `m`. -/
def AncPath (g : Graph) (σ : List Nat) (m : Nat) (sv : List Nat) : Prop :=
  ∀ k, k + 1 < sv.length → ∃ e, e < g.ne ∧ m ≤ σ.idxOf e ∧ Items.PairEq (sv[k]!, sv[k + 1]!) g.edges[e]!

theorem AncPath.mono {g : Graph} {m m' : Nat} {sv : List Nat} (h : m' ≤ m) (hp : AncPath g σ m sv) :
    AncPath g σ m' sv := fun k hk => by
  obtain ⟨e, he, hm, hpe⟩ := hp k hk
  exact ⟨e, he, le_trans h hm, hpe⟩

theorem getElem!_of_getElem? {l : List Nat} {k x : Nat} (h : l[k]? = some x) : l[k]! = x := by
  rw [List.getElem!_eq_getElem?_getD, h]; rfl

theorem getElem!_append_left' {l₁ l₂ : List Nat} {k : Nat} (h : k < l₁.length) : (l₁ ++ l₂)[k]! = l₁[k]! := by
  simp [List.getElem!_eq_getElem?_getD, List.getElem?_append_left h]

theorem getElem!_concat_length' (l : List Nat) (a : Nat) : (l ++ [a])[l.length]! = a :=
  getElem!_of_getElem? List.getElem?_concat_length

theorem PairEq.flip {a b : Nat} {q : Nat × Nat} (h : Items.PairEq (a, b) q) : Items.PairEq (b, a) q := by
  rcases h with h | h
  · subst h; exact .inr rfl
  · obtain ⟨rfl, rfl⟩ := Prod.mk.inj h; exact .inl rfl

theorem classify_back {d i : Nat} (hi : i < d + 1) :
    classify d false (i, d) = if i = d then .selfLoop else .ret i .backEdge := by
  unfold classify
  by_cases h : i = d
  · subst h; simp
  · have : ¬ (i ≥ d) := by omega
    simp [h, this]

abbrev DsTree (σ : List Nat) (n : Nat) (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  ∀ (g : Graph) (anc : List Nat) (pe : Nat → Prop), Types g s → d = anc.length → t.WF anc → t.Ends g →
    (∀ e', e' < g.ne → ∀ x, g.Inc e' x → x ∈ t.verts → e' ∈ t.edges ∨ pe e') →
    (∀ e', pe e' → ∀ x, g.Inc e' x → x ∈ anc ++ [t.v]) →
    (anc ++ t.verts).Nodup → (∀ a ∈ anc, a < g.nv) → (∀ v ∈ t.verts, v < g.nv) →
    (∀ e ∈ t.edges, e < g.ne) → t.edges.Nodup → s.stackVerts.size = g.nv →
    (∀ k, k < d → anc[k]? = some s.stackVerts[k]!) →
    σ.Nodup → PostAt σ n t.edgePostorder →
    (∀ v outs, t = .node v outs → AncPath g σ (n + t.edgePostorder.length) (anc ++ [v])) →
    CsTree σ n t d s
abbrev DsOuts (σ : List Nat) (n v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  ∀ (g : Graph) (anc : List Nat) (pe : Nat → Prop), Types g s → d = anc.length → (∀ o ∈ outs, o.WF anc v) →
    (∀ o ∈ outs, DfsOut.Ends g v o) →
    (∀ e', e' < g.ne → ∀ x, g.Inc e' x → x ∈ DfsOut.vertsList outs → e' ∈ DfsOut.edgesList outs ∨ pe e') →
    (∀ e', pe e' → ∀ x, g.Inc e' x → x ∈ anc ++ [v]) →
    (anc ++ v :: DfsOut.vertsList outs).Nodup →
    (∀ a ∈ anc, a < g.nv) → v < g.nv → (∀ w ∈ DfsOut.vertsList outs, w < g.nv) →
    (∀ e ∈ DfsOut.edgesList outs, e < g.ne) → (DfsOut.edgesList outs).Nodup →
    s.stackVerts.size = g.nv → (∀ k, k < d → anc[k]? = some s.stackVerts[k]!) → s.stackVerts[d]! = v →
    σ.Nodup → PostAt σ n (DfsOut.edgePostorderList outs) →
    AncPath g σ (n + (DfsOut.edgePostorderList outs).length) (anc ++ [v]) →
    CsOuts σ n v d outs hasVert s
abbrev DsOut (σ : List Nat) (n v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  ∀ (g : Graph) (anc : List Nat) (pe : Nat → Prop), Types g s → d = anc.length → o.WF anc v →
    DfsOut.Ends g v o →
    (∀ e', e' < g.ne → ∀ x, g.Inc e' x → x ∈ o.verts → subEdges o e' ∨ pe e') →
    (∀ e', pe e' → ∀ x, g.Inc e' x → x ∈ anc ++ [v]) →
    (anc ++ v :: o.verts).Nodup →
    (∀ a ∈ anc, a < g.nv) → v < g.nv → (∀ w ∈ o.verts, w < g.nv) →
    (∀ e ∈ o.edges, e < g.ne) → o.edges.Nodup →
    s.stackVerts.size = g.nv → (∀ k, k < d → anc[k]? = some s.stackVerts[k]!) → s.stackVerts[d]! = v →
    σ.Nodup → PostAt σ n o.block →
    AncPath g σ (n + o.block.length) (anc ++ [v]) →
    CsOut σ n v d o hasVert s

mutual
theorem dsTree : ∀ (σ : List Nat) (n : Nat) (t : DfsTree) (d : Nat) (s : WalkState), DsTree σ n t d s
  | σ, n, .node v outs, d, s => by
    intro g anc pe hT hd hwf hends hcomp hpeA hnd hanc hvlt helt hen hsz hsvk hσn hat hpath
    simp only [DfsTree.verts, DfsTree.edges] at hnd hvlt helt hen
    simp only [DfsTree.WF] at hwf
    simp only [DfsTree.Ends] at hends
    simp only [DfsTree.edgePostorder] at hat
    have hp := hpath v outs rfl
    simp only [DfsTree.edgePostorder] at hp
    unfold CsTree
    have hv : v < g.nv := hvlt v (List.mem_cons_self ..)
    have hdlt : d < g.nv := by
      have hnd' : (anc ++ [v]).Nodup :=
        List.Nodup.sublist (List.Sublist.append_left (List.cons_sublist_cons.2 (List.nil_sublist _)) _) hnd
      have hsub : anc ++ [v] ⊆ List.range g.nv := fun w hw => by
        rw [List.mem_range]
        rcases List.mem_append.1 hw with hw | hw
        · exact hanc w hw
        · exact (List.mem_singleton.1 hw) ▸ hv
      have := (List.Nodup.subperm hnd' hsub).length_le
      simp at this; omega
    exact dsOuts σ n v d outs false _ g anc pe ⟨hT.g_eq, hT.size, hT.root, hT.vert, hT.edge⟩ hd hwf.2 hends
      (fun e' he' x hx hxm => hcomp e' he' x hx (List.mem_cons_of_mem _ hxm)) hpeA
      hnd hanc hv (fun w hw => hvlt w (List.mem_cons_of_mem _ hw)) helt hen (by simp [hsz])
      (fun k hk => by rw [getElem!_set!_ne _ _ _ _ (by omega)]; exact hsvk k hk)
      (getElem!_set!_self _ _ _ (by rw [hsz]; exact hdlt)) hσn hat hp

theorem dsOuts : ∀ (σ : List Nat) (n v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState),
    DsOuts σ n v d outs hasVert s
  | σ, n, v, d, [], hasVert, s => by
    intro g anc pe hT hd hwf hends hcomp hpeA hnd hanc hv hw he hen hsz hsvk hsvd hσn hat hpath
    unfold CsOuts; trivial
  | σ, n, v, d, o :: rest, hasVert, s => by
    intro g anc pe hT hd hwf hends hcomp hpeA hnd hanc hv hw he hen hsz hsvk hsvd hσn hat hpath
    have hnd' := hnd
    rw [DfsOut.vertsList_eq, List.flatMap_cons, ← DfsOut.vertsList_eq] at hnd'
    have hndA : ∀ x, x ∈ anc → x ∈ v :: (o.verts ++ DfsOut.vertsList rest) → False :=
      fun x hx hx' => List.disjoint_of_nodup_append hnd' hx hx'
    have hndV : v ∉ o.verts ++ DfsOut.vertsList rest :=
      (List.nodup_cons.1 (List.nodup_append.1 hnd').2.1).1
    have hndOR : ∀ x, x ∈ o.verts → x ∈ DfsOut.vertsList rest → False :=
      fun x hx hx' =>
        List.disjoint_of_nodup_append (List.nodup_cons.1 (List.nodup_append.1 hnd').2.1).2 hx hx'
    have hmemR : ∀ x, x ∈ DfsOut.vertsList rest → x ∈ DfsOut.vertsList (o :: rest) := fun x hx => by
      rw [DfsOut.vertsList_eq, List.flatMap_cons, ← DfsOut.vertsList_eq]
      exact List.mem_append_right _ hx
    have hcompO : ∀ e', e' < g.ne → ∀ x, g.Inc e' x → x ∈ o.verts → subEdges o e' ∨ pe e' := by
      intro e' he' x hx hxo
      rcases hcomp e' he' x hx (mem_vertsList_of_verts (List.mem_cons_self ..) hxo) with h | h
      · obtain ⟨o', ho', hs⟩ := mem_subEdges_edgesList.1 h
        rcases List.mem_cons.1 ho' with heq | ho'
        · subst heq; exact .inl hs
        · exfalso
          rcases endsOut_wf g anc v o' (hwf o' (List.mem_cons_of_mem _ ho'))
            (hends o' (List.mem_cons_of_mem _ ho')) e' hs x hx with h' | rfl | h'
          · exact hndA x h' (List.mem_cons_of_mem _ (List.mem_append_left _ hxo))
          · exact hndV (List.mem_append_left _ hxo)
          · exact hndOR x hxo (mem_vertsList_of_verts ho' h')
      · exact .inr h
    have hcompR : ∀ e', e' < g.ne → ∀ x, g.Inc e' x → x ∈ DfsOut.vertsList rest →
        e' ∈ DfsOut.edgesList rest ∨ pe e' := by
      intro e' he' x hx hxr
      rcases hcomp e' he' x hx (hmemR x hxr) with h | h
      · obtain ⟨o', ho', hs⟩ := mem_subEdges_edgesList.1 h
        rcases List.mem_cons.1 ho' with heq | ho'
        · exfalso
          subst heq
          rcases endsOut_wf g anc v o' (hwf o' (List.mem_cons_self ..))
            (hends o' (List.mem_cons_self ..)) e' hs x hx with h' | rfl | h'
          · exact hndA x h' (List.mem_cons_of_mem _ (List.mem_append_right _ hxr))
          · exact hndV (List.mem_append_right _ hxr)
          · exact hndOR x h' hxr
        · exact .inl (mem_subEdges_edgesList.2 ⟨o', ho', hs⟩)
      · exact .inr h
    rw [DfsOut.vertsList_eq, List.flatMap_cons, ← DfsOut.vertsList_eq] at hnd hw
    rw [DfsOut.edgesList_eq, List.flatMap_cons, ← DfsOut.edgesList_eq] at he hen
    rw [DfsOut.edgePostorderList_cons] at hat hpath
    unfold CsOuts
    refine ⟨dsOut σ n v d o hasVert s g anc pe hT hd (hwf o (List.mem_cons_self ..)) (hends o (List.mem_cons_self ..))
      hcompO hpeA
      (List.Nodup.sublist (List.Sublist.append_left (List.Sublist.cons_cons v (List.sublist_append_left _ _)) _) hnd)
      hanc hv (fun w hw' => hw w (List.mem_append_left _ hw')) (fun e he' => he e (List.mem_append_left _ he'))
      (List.Nodup.sublist (List.sublist_append_left _ _) hen) hsz hsvk hsvd hσn hat.left
      (hpath.mono (by rw [List.length_append]; omega)), ?_⟩
    have hK := kOut v d o hasVert s g (d + 1) 0 s hT (Nat.le_refl _) (by omega) (vertItem_ne_zero v)
      (fun w _ => vertItem_ne_zero w) (fun e _ => edgeItem_ne_zero g e) Keep.refl
    refine wp_mono _ hK fun hv' s' hK' => ?_
    exact dsOuts σ (n + o.block.length) v d rest hv' s' g anc pe (hT.of_keep hK') hd
      (fun o' ho' => hwf o' (List.mem_cons_of_mem _ ho')) (fun o' ho' => hends o' (List.mem_cons_of_mem _ ho'))
      hcompR hpeA
      (List.Nodup.sublist (List.Sublist.append_left (List.Sublist.cons_cons v (List.sublist_append_right _ _)) _) hnd)
      hanc hv (fun w hw' => hw w (List.mem_append_right _ hw')) (fun e he' => he e (List.mem_append_right _ he'))
      (List.Nodup.sublist (List.sublist_append_right _ _) hen) (hK'.sv.trans hsz)
      (fun k hk => by rw [hK'.svlo k (by omega)]; exact hsvk k hk)
      (by rw [hK'.svlo d (Nat.lt_succ_self _)]; exact hsvd) hσn hat.right
      (hpath.mono (by rw [List.length_append, Nat.add_assoc]))

theorem dsOut : ∀ (σ : List Nat) (n v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState), DsOut σ n v d o hasVert s
  | σ, n, v, d, o, hasVert, s => by
    intro g anc pe hT hd hwf hends hcomp hpeA hnd hanc hv hw he hen hsz hsvk hsvd hσn hat hpath
    subst hd
    unfold CsOut
    have hK₀ := keep_walkOutPre (D := anc.length + 1) (j := 0) v anc.length o hasVert (Keep.refl (s := s))
    refine wp_mono _ hK₀ fun hv' s₁ hK₁ => ?_
    clear hK₀
    have hg₁ : s₁.g = g := hK₁.g.trans hT.g_eq
    have hsv₁ : ∀ k, k ≤ anc.length → s₁.stackVerts[k]! = (anc ++ [v])[k]! := fun k hk => by
      rw [hK₁.svlo k (by omega)]
      rcases Nat.lt_or_ge k anc.length with hk' | hk'
      · exact (getElem!_of_getElem? ((List.getElem?_append_left hk').trans (hsvk k hk'))).symm
      · have : k = anc.length := by omega
        subst this
        rw [hsvd, getElem!_concat_length']
    cases o with
    | back e dest cls =>
      simp only [DfsOut.WF] at hwf
      obtain ⟨i, hi, hcls⟩ := hwf
      simp only [DfsOut.block] at hat hpath
      have hil : i < anc.length + 1 := by
        have := (List.getElem?_eq_some_iff.1 hi).1; simpa using this
      subst hcls
      rw [classify_back hil] at hends ⊢
      refine ⟨hat.singleton, hat.block_infix, by rw [hg₁]; exact hends, ?_, ?_, ?_, ?_, ?_⟩
      · intro k hk
        obtain ⟨e', he', hm, hpe⟩ := hpath k (by simp; omega)
        refine ⟨e', hg₁ ▸ he', by simp at hm; omega, ?_⟩
        rw [hsv₁ k (by omega), hsv₁ (k + 1) (by omega), hg₁]; exact hpe
      · intro ht
        simp only [DfsOut.cls] at ht
        split at ht <;> simp [OutClass.isTree] at ht
      · intro ht
        simp only [DfsOut.cls] at ht
        split at ht <;> simp [OutClass.isTree] at ht
      · show dest < s₁.g.nv
        rw [hg₁]
        rcases List.mem_append.1 (List.mem_of_getElem? hi) with h | h
        · exact hanc _ h
        · exact (List.mem_singleton.1 h) ▸ hv
      · intro _ hge
        simp only [DfsOut.cls, DfsOut.dest] at hge ⊢
        split at hge
        · rename_i hid
          subst hid
          rw [List.getElem?_concat_length] at hi
          exact (Option.some.inj hi).symm
        · simp [OutClass.lowval] at hge; omega
    | tree e cls child =>
      simp only [DfsOut.WF] at hwf
      obtain ⟨hwfc, hcls⟩ := hwf
      simp only [DfsOut.Ends] at hends
      obtain ⟨hpe, hendsc⟩ := hends
      simp only [DfsOut.verts] at hnd hw
      simp only [DfsOut.edges, List.mem_cons, forall_eq_or_imp] at he
      simp only [DfsOut.edges] at hen
      simp only [DfsOut.block] at hat hpath
      have htr : cls.isTree = true := by
        cases hc : cls with
        | bridge => rfl
        | component => rfl
        | selfLoop => rw [hc] at hcls; unfold classify at hcls; split at hcls <;> split at hcls <;> cases hcls
        | ret lv kind =>
          rw [hc] at hcls
          obtain ⟨_, _, hk⟩ := classify_eq_ret_iff_tree.1 hcls.symm
          cases kind with
          | backEdge => split at hk <;> cases hk
          | _ => rfl
      simp only [wp_modify]
      have hT₁ := hT.of_keep hK₁
      have hT₂ : Types g { s₁ with firstOccurrence := s₁.firstOccurrence.set! anc.length s₁.g.ne } :=
        ⟨hT₁.g_eq, hT₁.size, hT₁.root, hT₁.vert, hT₁.edge⟩
      have hnd' : (anc ++ [v] ++ child.verts).Nodup := by simpa using hnd
      have hpos : σ.idxOf e = n + child.edgePostorder.length := by
        have h1 := hat.right.singleton
        have h2 := (List.getElem?_eq_some_iff.1 h1).1
        rw [← getElem!_of_getElem? h1]; exact idxOf_getElem! hσn h2
      have hpc : ∀ w outs, child = .node w outs →
          AncPath g σ (n + child.edgePostorder.length) (anc ++ [v] ++ [w]) := by
        rintro w outs rfl k hk
        simp only [List.length_append, List.length_singleton] at hk
        rcases Nat.lt_or_ge (k + 1) (anc.length + 1) with hk' | hk'
        · obtain ⟨e', he', hm, hpe'⟩ := hpath k (by simp; omega)
          refine ⟨e', he', by simp at hm; omega, ?_⟩
          rw [getElem!_append_left' (l₁ := anc ++ [v]) (l₂ := [w]) (k := k) (by simp; omega),
            getElem!_append_left' (l₁ := anc ++ [v]) (l₂ := [w]) (k := k + 1) (by simp; omega)]
          exact hpe'
        · have : k = anc.length := by omega
          subst this
          refine ⟨e, he.1, hpos.symm.le, ?_⟩
          have hl : anc.length + 1 = (anc ++ [v]).length := by simp
          rw [getElem!_append_left' (l₁ := anc ++ [v]) (l₂ := [w]) (by simp), getElem!_concat_length', hl,
            getElem!_concat_length']
          simp only [DfsTree.v] at hpe
          exact PairEq.flip hpe
      have hcompC : ∀ e', e' < g.ne → ∀ x, g.Inc e' x → x ∈ child.verts →
          e' ∈ child.edges ∨ (e' = e ∨ pe e') := by
        intro e' he' x hx hxc
        rcases hcomp e' he' x hx hxc with h | h
        · rcases h with h | h
          · exact .inr (.inl h)
          · exact .inl h
        · exact .inr (.inr h)
      have hpeC : ∀ e', (e' = e ∨ pe e') → ∀ x, g.Inc e' x → x ∈ anc ++ [v] ++ [child.v] := by
        intro e' h x hx
        rcases h with rfl | h
        · rcases eq_of_inc_pairEq hpe hx with h | h <;> rw [h]
          · exact List.mem_append_right _ (List.mem_singleton_self _)
          · exact List.mem_append_left _ (List.mem_append_right _ (List.mem_singleton_self _))
        · exact List.mem_append_left _ (hpeA e' h x hx)
      have hcv : child.v ∈ child.verts := by
        obtain ⟨w, outs⟩ := child; exact List.mem_cons_self ..
      refine ⟨dsTree σ n child (anc.length + 1) _ g (anc ++ [v]) (fun e' => e' = e ∨ pe e') hT₂ (by simp)
        hwfc hendsc hcompC hpeC hnd'
        (fun a ha => by
          rcases List.mem_append.1 ha with ha | ha
          · exact hanc a ha
          · exact (List.mem_singleton.1 ha) ▸ hv)
        hw he.2 (List.nodup_cons.1 hen).2 (hK₁.sv.trans hsz)
        (fun k hk => by
          rcases Nat.lt_or_ge k anc.length with hk' | hk'
          · rw [List.getElem?_append_left hk', hsvk k hk']
            show _ = some s₁.stackVerts[k]!
            rw [hK₁.svlo k (by omega)]
          · have : k = anc.length := by omega
            subst this
            rw [List.getElem?_concat_length]
            show _ = some s₁.stackVerts[anc.length]!
            rw [hK₁.svlo _ (Nat.lt_succ_self _), hsvd])
        hσn hat.left hpc, ?_⟩
      have hKc := kTree child (anc.length + 1) _ g (anc.length + 1) 0 _ hT₂ (Nat.le_refl _) (by omega)
        (fun w _ => vertItem_ne_zero w) (fun e' _ => edgeItem_ne_zero g e') Keep.refl
      refine wp_mono _ hKc fun _ s₃ hK₃ => ?_
      have hg₃ : s₃.g = g := hK₃.g.trans hg₁
      have hsv₃ : ∀ k, k ≤ anc.length → s₃.stackVerts[k]! = (anc ++ [v])[k]! := fun k hk => by
        rw [hK₃.svlo k (by omega)]; exact hsv₁ k hk
      refine ⟨hat.right.singleton, hat.block_infix, by rw [hg₃]; simp only [DfsOut.Ends]; exact ⟨hpe, hendsc⟩, ?_, ?_, ?_, ?_, ?_⟩
      · intro k hk
        obtain ⟨e', he', hm, hpe'⟩ := hpath k (by simp; omega)
        refine ⟨e', hg₃ ▸ he', by simp at hm; omega, ?_⟩
        rw [hsv₃ k (by omega), hsv₃ (k + 1) (by omega), hg₃]; exact hpe'
      · intro _ hlt
        simp only [DfsOut.cls] at hlt
        rw [hcls] at hlt
        obtain ⟨w, outs⟩ := child
        simp only [DfsTree.Ends] at hendsc
        cases outs with
        | nil =>
          exfalso
          simp [DfsTree.retDepths, DfsOut.retDepthsList, low2, classify, OutClass.lowval] at hlt
        | cons o' rest' =>
          have hmem : o'.e ∈ (DfsTree.node w (o' :: rest')).edges := by
            cases o' <;> simp [DfsTree.edges, DfsOut.edgesList, DfsOut.e]
          refine ⟨o'.e, by rw [hg₃]; exact he.2 _ hmem, fun h => (List.nodup_cons.1 hen).1 (by simp only [DfsOut.e] at h; exact h ▸ hmem), ?_⟩
          show s₃.g.Inc o'.e w
          rw [hg₃]
          have hE := hendsc o' (List.mem_cons_self ..)
          cases o' with
          | back e₀ d₀ c₀ => simp only [DfsOut.Ends] at hE; exact (Graph.inc_of_pairEq hE).1
          | tree e₀ c₀ ch => simp only [DfsOut.Ends] at hE; exact (Graph.inc_of_pairEq hE.1).2
      · intro _ e' he' hinc
        rw [hg₃] at he' hinc
        rcases hcomp e' he' _ hinc hcv with h | h
        · exact h
        · exact (List.disjoint_of_nodup_append hnd' (hpeA e' h _ hinc) hcv).elim
      · show child.v < s₃.g.nv
        rw [hg₃]
        obtain ⟨w, outs⟩ := child
        exact hw w (by simp [DfsTree.verts])
      · intro h
        simp only [DfsOut.cls] at h
        rw [h] at htr; cases htr
end

/-- `CsTree` on the real walk: a DFS tree below the ancestor path `anc`, scheduled at `n`. -/
theorem csTree_of_dfs {g : Graph} {anc : List Nat} (t : DfsTree) (d : Nat) (s : WalkState)
    (hT : Types g s) (hd : d = anc.length) (hwf : t.WF anc) (hends : t.Ends g)
    (hcomp : ∀ e, e < g.ne → ∀ x, g.Inc e x → x ∈ t.verts → e ∈ t.edges)
    (hnd : (anc ++ t.verts).Nodup) (hanc : ∀ a ∈ anc, a < g.nv) (hvlt : ∀ v ∈ t.verts, v < g.nv)
    (helt : ∀ e ∈ t.edges, e < g.ne) (hen : t.edges.Nodup) (hsz : s.stackVerts.size = g.nv)
    (hsvk : ∀ k, k < d → anc[k]? = some s.stackVerts[k]!) (hσn : σ.Nodup) (hat : PostAt σ n t.edgePostorder)
    (hpath : ∀ v outs, t = .node v outs → AncPath g σ (n + t.edgePostorder.length) (anc ++ [v])) :
    CsTree σ n t d s :=
  dsTree σ n t d s g anc (fun _ => False) hT hd hwf hends
    (fun e he x hx hxt => .inl (hcomp e he x hx hxt)) (fun _ h => h.elim)
    hnd hanc hvlt helt hen hsz hsvk hσn hat hpath

theorem forest_closeInv {g : Graph} {σ : List Nat} (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < g.ne) :
    ∀ forest pre n s, RootState g pre s → s.RangesInv σ n 0 → ForestOK g (pre ++ forest) →
      (∀ t ∈ forest, t.WF []) → (∀ t ∈ forest, t.Ends g) →
      (∀ t ∈ forest, ∀ e, e < g.ne → ∀ x, g.Inc e x → x ∈ t.verts → e ∈ t.edges) →
      RootsCover σ n forest s →
      PostAt σ n (edgePostorderForest forest) → s.CloseInv →
      wp (walkForest forest) (fun _ s' => s'.CloseInv) s
  | [], _, _, _, _, _, _, _, _, _, _, _, hcl => by simpa [walkForest, wp_pure] using hcl
  | t :: rest, pre, n, s, h, hr, hf, hwf, hends, hcomp, hc, hat, hcl => by
    have hb := h.book hf (hwf t (by simp)) (hends t (by simp)) (hcomp t (by simp))
    have hg := gbTree t 0 s hb
    have hi : ∀ v outs, t = .node v outs →
        ({ s with stackVerts := s.stackVerts.set! 0 v } : WalkState).RangesInv σ n 0 :=
      fun v _ _ => hr.stackVerts_of_nil h.tstack _
    have hσ' : ∀ e ∈ σ, e < s.g.ne := by rwa [h.g_eq]
    have hfront := walkTree_frontiers t 0 s (fun v outs ht => (hi v outs ht).inv) h.shape hg hb
    have hat' : PostAt σ n (t.edgePostorder ++ edgePostorderForest rest) := hat
    have hsched := scheduleTree σ n t 0 s hi h.shape hnd hσ' hg hb hfront hc.1 hat'.left
    have hrg := rgTree σ n t 0 s hi h.shape hnd hσ' hg hb hsched
    have hcs : CsTree σ n t 0 s :=
      csTree_of_dfs (anc := []) t 0 s h.place.types rfl (hwf t (by simp)) (hends t (by simp))
        (hcomp t (by simp))
        (by simpa using (RootState.hvn hf).1) (by simp) (RootState.hvlt hf) (RootState.helt hf)
        (RootState.hen hf).1 h.sv (fun k hk => absurd hk (Nat.not_lt_zero _)) hnd hat'.left
        (fun v outs _ k hk => absurd hk (by simp))
    have hcc := ccTree σ n t 0 s hi h.shape hnd hσ' hg hb hfront hsched hcs hcl
    have hnv : 0 < s.g.nv := by
      obtain ⟨v, outs⟩ := t
      have hv := RootState.hvlt hf v (by simp [DfsTree.verts])
      rw [h.g_eq]; omega
    have hk : wp (walkTree t 0) (fun _ s' => RootOK s') s :=
      walkTree_rootOK t s (fun v outs ht => (hi v outs ht).inv) h.shape hg hb h.tstack
      (by rw [h.sd, ← h.g_eq]; exact hnv)
    have hp := (walk_place_aux g).1 t 0 _ _ s h.place (RootState.hvlt hf) (RootState.helt hf)
      (RootState.hvn hf).1 (RootState.hen hf).1 (RootState.hPv hf) (RootState.hPe hf)
    have hst := h.step hf (hwf t (by simp)) (hends t (by simp)) (hcomp t (by simp))
    show wp ((walkTree t 0 >>= fun _ => popTstack >>= fun top =>
      modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 }) >>= fun _ => walkForest rest) _ s
    simp only [wp_bind]
    refine wp_imp (wp_of_forall fun _ s₁ ⟨hrs, hk, hp, hst, hc, hcl⟩ => ?_)
      (wp_and hrg (wp_and hk (wp_and hp (wp_and hst (wp_and hc.2 hcc)))))
    have hrpop := hrs.1.root_append hrs.2.1 hk (noParent_of_cnt_eq_zero hp.root)
    have hp' : s₁.Place s₁.g (Pushed g (Pushed g (fun _ => False) (pre.flatMap DfsTree.verts)
        (pre.flatMap DfsTree.edges)) t.verts t.edges) (fun _ => False) := by
      rw [hp.g_eq]; exact hp
    have hclpop := rootAppend_closeInv hcl hp' hk
    exact wp_mono (popTstack >>= fun top => modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 })
      (wp_and hrpop (wp_and hst (wp_and hc hclpop))) fun _ s₂ ⟨hr₂, hs₂, hc₂, hcl₂⟩ =>
      forest_closeInv hnd hσ rest (pre ++ [t]) (n + t.edgePostorder.length) s₂ hs₂ hr₂
        (by simpa using hf) (fun t' ht' => hwf t' (by simp [ht']))
        (fun t' ht' => hends t' (by simp [ht'])) (fun t' ht' => hcomp t' (by simp [ht'])) hc₂ hat'.right hcl₂

/-- `CloseInv` for the walk of a DFS forest, from `RootsCover` (the range-side schedule). -/
theorem walk_closeInv' (g : Graph) (tern : Bool) (forest : List DfsTree)
    (hf : ForestOK g forest) (hwf : ∀ t ∈ forest, t.WF []) (hends : ∀ t ∈ forest, t.Ends g)
    (hecov : ∀ e, e < g.ne → e ∈ forest.flatMap DfsTree.edges)
    (hc : RootsCover (edgePostorderForest forest) 0 forest (WalkState.init g tern)) :
    (g.walk tern forest).CloseInv := by
  have hnd := DfsData.edgePostorderForest_perm.nodup_iff.2 hf.edges_nodup
  have hσ : ∀ e ∈ edgePostorderForest forest, e < g.ne :=
    fun e he => hf.edges_lt e (DfsData.edgePostorderForest_perm.subset he)
  have h := forest_closeInv hnd hσ forest [] 0 (WalkState.init g tern)
    (rootState_init g tern) (init_rangesInv g tern hnd) (by simpa using hf) hwf hends
    (fun t ht => comp_of_forest hf hwf hends hecov ht) hc
    ⟨[], [], rfl, by simp⟩ (init_closeInv g tern)
  simpa only [wp, Graph.walk] using h

end WalkState
end Spqr
