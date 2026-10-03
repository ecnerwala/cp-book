import Spqr.RangesStep
import Spqr.WalkInv

/-!
# `walkTree_rangesInv`: the range invariant along the walk

Mirrors `walkTree_inv'` (`WalkInv.lean`): the same mutual induction over `walkTree`/`walkOuts`/
`walkOut`, with the ear guards (`GuardsTree`), the bookkeeping (`BookTree`) and, new here, the
range-side hypotheses `RgTree σ n t d s`: at every `finishEdge` the `FinishAdj`/`BoundaryAdj`
bundle (edge position in `σ`, adjacency at each merge), and at every vertex push `PushVertR`. The
edge counter `n` advances by one per `finishEdge`, i.e. by `o.block.length` per out-edge and by
`t.edgePostorder.length` per subtree.

`finishBoundary_rangesInv` (the block-closing edge, `d ≤ o.cls.lowval d`) is the one admission of
this file; its ear-side counterpart `ear_boundary` is admitted in the ear layer.
-/

namespace Spqr
open WalkM

theorem DfsOut.edgePostorderList_cons (o : DfsOut) (rest : List DfsOut) :
    DfsOut.edgePostorderList (o :: rest) = o.block ++ DfsOut.edgePostorderList rest := by
  cases o <;> simp [DfsOut.edgePostorderList, DfsOut.block]

namespace WalkState
variable {s : WalkState} {σ : List Nat} {n D : Nat}

theorem RangesInv.withInv {D' : Nat} (h : s.RangesInv σ n D) (hi : s.Inv' D') : s.RangesInv σ n D' :=
  ⟨hi, h.processed, h.ordered, h.convex, h.closed⟩

theorem RangesInv.setSv {d : Nat} (x : Nat) (h : s.RangesInv σ n d) :
    ({ s with stackVerts := s.stackVerts.set! (d + 1) x } : WalkState).RangesInv σ n (d + 1) :=
  ⟨h.inv.setSv x, h.processed, h.ordered, h.convex, h.closed⟩

theorem RangesInv.frame' {s' : WalkState} (h : s.RangesInv σ n D) (hg : s'.g = s.g := by rfl)
    (hsv : s'.stackVerts = s.stackVerts := by rfl) (hi : s'.items = s.items := by rfl)
    (hts : s'.tstack = s.tstack := by rfl) : s'.RangesInv σ n D :=
  h.frame hg hsv hi hts

/-- Range-side hypotheses of `finishBoundary` (the block-closing edge `o.e`, `d ≤ o.cls.lowval d`):
the edge sits at position `n` of `σ`, the entries popped into its `Q` item carry the block on the
side that is moved (`side`), and their edges fill `σ` up to `n` from the first piece edge (`fill`),
so the closed `Q` item is convex. -/
structure BoundaryAdj (σ : List Nat) (n d : Nat) (o : DfsOut) (s : WalkState) : Prop where
  pos : σ[n]? = some o.e
  bridge : o.cls.isTree = true → o.cls.lowval d = d + 1 → ∀ t rest, s.tstack = t :: rest →
    t.spans.1 = [] ∧ ∀ a b, a ≤ b → b < n → t.piece s.g s.items σ[a]! → t.edges s.g s.items σ[b]!
  component : o.cls.isTree = true → o.cls.lowval d ≠ d + 1 → ∀ t₁ t₂ rest, s.tstack = t₁ :: t₂ :: rest →
    t₁.spans.2 = [] ∧ t₂.spans.1 = [] ∧ ∀ a b, a ≤ b → b < n →
      (t₁.piece s.g s.items σ[a]! ∨ t₂.piece s.g s.items σ[a]!) →
      (t₁.edges s.g s.items σ[b]! ∨ t₂.edges s.g s.items σ[b]!)

/-- ADMITTED: `finishBoundary` preserves the range invariant, advancing `n`, under `BoundaryOk`
(ear layer) and `BoundaryAdj`. -/
theorem finishBoundary_rangesInv {curV d : Nat} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}
    (hge : d ≤ o.cls.lowval d) (h : s.RangesInv σ n D) (hs : Shape s) (hnd : σ.Nodup)
    (hσ : ∀ e ∈ σ, e < s.g.ne) (hb : FinishBook curV d o origTstack hasVert s)
    (hok : BoundaryOk D curV d o s) (hadj : BoundaryAdj σ n d o s) :
    (after (finishEdge curV d o origTstack hasVert) s).RangesInv σ (n + 1) D ∧
      (after (finishEdge curV d o origTstack hasVert) s).g = s.g := by
  sorry

/-- The range hypotheses of one `finishEdge` call, at the state where it runs. -/
def FinishR (σ : List Nat) (n curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool) (s : WalkState) :
    Prop :=
  (∀ lv kind, o.cls = .ret lv kind → lv < d → FinishAdj σ n curV d lv o origTstack hasVert s) ∧
  (d ≤ o.cls.lowval d → BoundaryAdj σ n d o s)

/-- `finishEdge` preserves the range invariant (at depth `d`, like `finishEdge_step`) and advances
`n`, under the guards, the bookkeeping and `FinishR`. -/
theorem finishEdge_ranges {v d : Nat} {o : DfsOut} {m : Nat} {hasVert : Bool} {D : Nat}
    (hD : D = if o.cls.isTree then d + 1 else d) (h : s.RangesInv σ n D) (hs : Shape s) (hnd : σ.Nodup)
    (hσ : ∀ e ∈ σ, e < s.g.ne) (hg : FinishGuards d o m hasVert s) (hb : FinishBook v d o m hasVert s)
    (hr : FinishR σ n v d o m hasVert s) :
    (after (finishEdge v d o m hasVert) s).RangesInv σ (n + 1) d ∧ Shape (after (finishEdge v d o m hasVert) s) ∧
      (after (finishEdge v d o m hasVert) s).g = s.g := by
  obtain ⟨hi', hs'⟩ := finishEdge_step hD h.inv hs hg hb
  by_cases hge : d ≤ o.cls.lowval d
  · obtain ⟨hr', hgg⟩ := finishBoundary_rangesInv hge h hs hnd hσ hb (ear_boundary hge hg h.inv hs hb hD) (hr.2 hge)
    exact ⟨hr'.withInv hi', hs', hgg⟩
  · obtain ⟨lv, kind, ho, hl⟩ := ret_of_lowval_lt (Nat.lt_of_not_le hge)
    obtain ⟨sub, base, hlen, hE⟩ := hb.ear
    have hok := finishOk_of_guards ho hl hg hE hlen h.inv hs hD hb.v_lt hb.e_lt hb.q (hb.ends lv kind ho) hb.vert
    have hr' := finishEdge_rangesInv v d lv kind o m hasVert ho hl hb.v_lt h hs hnd hσ hok (hr.1 lv kind ho hl)
    have st := finishEdge_inv v d lv kind o m hasVert ho hl hb.v_lt h.inv hs hok
    exact ⟨hr'.withInv hi', hs', st.g⟩

mutual
/-- The range-side hypotheses of the walk, in the style of `GuardsTree`. -/
def RgTree (σ : List Nat) (n : Nat) (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  match t with
  | .node v outs => RgOuts σ n v d outs false { s with stackVerts := s.stackVerts.set! d v }

def RgOuts (σ : List Nat) (n v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  match outs with
  | [] => hasVert = false → PushVertR σ n v s
  | o :: rest => RgOut σ n v d o hasVert s ∧
      wp (walkOut v d o hasVert) (fun hasVert' s' => RgOuts σ (n + o.block.length) v d rest hasVert' s') s

def RgOut (σ : List Nat) (n v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  (hasVert = false → PushVertR σ n v s) ∧
  wp (walkOutPre v d o hasVert) (fun hasVert' s₁ =>
    match o with
    | .tree _ _ child =>
      wp (modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }) (fun _ s₂ =>
        RgTree σ n child (d + 1) s₂ ∧
        wp (walkTree child (d + 1)) (fun _ s₃ =>
          FinishR σ (n + child.edgePostorder.length) v d o s₁.tstack.length hasVert' s₃) s₂) s₁
    | .back .. => FinishR σ n v d o s₁.tstack.length hasVert' s₁) s
end

theorem lt_ne_of_g {s s' : WalkState} (hg : s'.g = s.g) (hσ : ∀ e ∈ σ, e < s.g.ne) :
    ∀ e ∈ σ, e < s'.g.ne := by
  rw [hg]; exact hσ

theorem walkOutPre_ranges {v d : Nat} {o : DfsOut} {hasVert : Bool} (hi : s.RangesInv σ n d) (hs : Shape s)
    (hσ : ∀ e ∈ σ, e < s.g.ne) (hb : VertBook v hasVert s) (hr : hasVert = false → PushVertR σ n v s) :
    wp (walkOutPre v d o hasVert) (fun _ s₁ => s₁.RangesInv σ n d ∧ Shape s₁ ∧ ∀ e ∈ σ, e < s₁.g.ne) s := by
  have hi' : ({ s with stackDir := s.stackDir.set! d (if o.cls.lowval d ≥ d then false else !s.stackDir[o.cls.lowval d]!) } : WalkState).RangesInv σ n d :=
    hi.frame'
  have hs' : Shape ({ s with stackDir := s.stackDir.set! d (if o.cls.lowval d ≥ d then false else !s.stackDir[o.cls.lowval d]!) } : WalkState) :=
    hs.frame'
  unfold walkOutPre
  simp only [wp_bind, wp_stackDir, wp_setStackDir, wp_ite]
  split
  · rename_i hc
    have hf : hasVert = false := by cases hasVert <;> simp_all
    obtain ⟨hv, hc, ha⟩ := hb hf
    have st := RgStep.pushVert (v := v) hi' hs' d hv hc ha ⟨(hr hf).vtype, (hr hf).below⟩
    exact ⟨st.ranges, st.step.shape, lt_ne_of_g st.step.g hσ⟩
  · exact ⟨hi', hs', hσ⟩

abbrev RgS (σ : List Nat) (n d : Nat) (s : WalkState) : Prop :=
  s.RangesInv σ n d ∧ Shape s ∧ ∀ e ∈ σ, e < s.g.ne

abbrev RgInvTree (σ : List Nat) (n : Nat) (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  (∀ v outs, t = .node v outs → ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).RangesInv σ n d) →
  Shape s → σ.Nodup → (∀ e ∈ σ, e < s.g.ne) → GuardsTree t d s → BookTree t d s → RgTree σ n t d s →
  wp (walkTree t d) (fun _ s' => RgS σ (n + t.edgePostorder.length) d s') s

abbrev RgInvOuts (σ : List Nat) (n v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  RgS σ n d s → σ.Nodup → GuardsOuts v d outs hasVert s → BookOuts v d outs hasVert s →
  RgOuts σ n v d outs hasVert s →
  wp (walkOuts v d outs hasVert) (fun hasVert' s' =>
    RgS σ (n + (DfsOut.edgePostorderList outs).length) d s' ∧ VertBook v hasVert' s' ∧
      (hasVert' = false → PushVertR σ (n + (DfsOut.edgePostorderList outs).length) v s')) s

abbrev RgInvOut (σ : List Nat) (n v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  RgS σ n d s → σ.Nodup → GuardsOut v d o hasVert s → BookOut v d o hasVert s → RgOut σ n v d o hasVert s →
  wp (walkOut v d o hasVert) (fun _ s' => RgS σ (n + o.block.length) d s') s

mutual
theorem rgTree : ∀ (σ : List Nat) (n : Nat) (t : DfsTree) (d : Nat) (s : WalkState), RgInvTree σ n t d s
  | σ, n, .node v outs, d, s => fun hi hs hnd hσ hg hb hr => by
    unfold walkTree
    simp only [wp_bind, wp_modify]
    unfold GuardsTree at hg; unfold BookTree at hb; unfold RgTree at hr
    refine wp_imp (wp_of_forall fun hv s' ⟨⟨hi', hs', hσ'⟩, hvb, hpr⟩ => ?_)
      (rgOuts σ n v d outs false _ ⟨hi v outs rfl, hs.frame', hσ⟩ hnd hg hb hr)
    cases hv
    · simp only [Bool.false_eq_true, ↓reduceIte, wp_bind, wp_setStackDir]
      obtain ⟨hv, hc, ha⟩ := hvb rfl
      have hi₂ : ({ s' with stackDir := s'.stackDir.set! d true } : WalkState).RangesInv σ _ d := hi'.frame'
      have hs₂ : Shape ({ s' with stackDir := s'.stackDir.set! d true } : WalkState) := hs'.frame'
      have st := RgStep.pushVert (v := v) hi₂ hs₂ d hv hc ha ⟨(hpr rfl).vtype, (hpr rfl).below⟩
      exact ⟨st.ranges, st.step.shape, lt_ne_of_g st.step.g hσ'⟩
    · exact ⟨hi', hs', hσ'⟩

theorem rgOuts : ∀ (σ : List Nat) (n v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState),
    RgInvOuts σ n v d outs hasVert s
  | σ, n, v, d, [], hasVert, s => fun ⟨hi, hs, hσ⟩ _ _ hb hr => by
    unfold BookOuts at hb; unfold RgOuts at hr
    unfold walkOuts
    simp only [wp_pure, DfsOut.edgePostorderList, List.length_nil, Nat.add_zero]
    exact ⟨⟨hi, hs, hσ⟩, hb, hr⟩
  | σ, n, v, d, o :: rest, hasVert, s => fun hrs hnd hg hb hr => by
    unfold GuardsOuts at hg; unfold BookOuts at hb; unfold RgOuts at hr
    unfold walkOuts
    simp only [wp_bind]
    rw [DfsOut.edgePostorderList_cons, List.length_append, ← Nat.add_assoc]
    exact wp_imp (wp_imp (wp_imp (wp_imp (wp_of_forall fun hv' s' hrs' hg' hb' hr' =>
      rgOuts σ (n + o.block.length) v d rest hv' s' hrs' hnd hg' hb' hr') (rgOut σ n v d o hasVert s hrs hnd hg.1 hb.1 hr.1)) hg.2) hb.2) hr.2

theorem rgOut : ∀ (σ : List Nat) (n v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState), RgInvOut σ n v d o hasVert s
  | σ, n, v, d, o, hasVert, s => fun ⟨hi, hs, hσ⟩ hnd hg hb hr => by
    rw [walkOut_eq, wp_bind]
    unfold GuardsOut at hg; unfold BookOut at hb; unfold RgOut at hr
    refine wp_imp (wp_imp (wp_imp (wp_imp (wp_of_forall fun hv' s₁ ⟨hi₁, hs₁, hσ₁⟩ hg₁ hb₁ hr₁ => ?_)
      (walkOutPre_ranges hi hs hσ hb.1 hr.1)) hg) hb.2) hr.2
    unfold walkOutRest
    rw [wp_bind, wp_tstackSize]
    cases o with
    | back e cls dest =>
      try simp only [wp_pure] at hg₁ hb₁ hr₁
      simp only [wp_bind, wp_pure, DfsOut.block, List.length_singleton]
      obtain ⟨h₁, h₂, h₃⟩ := finishEdge_ranges (D := d)
        (by rw [if_neg fun h => by obtain ⟨_, _, _, h⟩ := hb₁.tree.1 h; cases h]) hi₁ hs₁ hnd hσ₁ hg₁ hb₁ hr₁
      exact ⟨h₁, h₂, lt_ne_of_g h₃ hσ₁⟩
    | tree e cls child =>
      try simp only [wp_bind, wp_modify] at hg₁ hb₁ hr₁
      simp only [wp_bind, wp_modify, DfsOut.block, List.length_append, List.length_singleton, ← Nat.add_assoc]
      refine wp_imp (wp_imp (wp_of_forall fun _ s₃ ⟨hg₃, hb₃, hr₃⟩ ⟨hi₃, hs₃, hσ₃⟩ => ?body)
        (wp_and hg₁.2 (wp_and hb₁.2 hr₁.2)))
        (rgTree σ n child (d + 1) _ ?pre hs₁.frame' hnd hσ₁ hg₁.1 hb₁.1 hr₁.1)
      case body =>
        obtain ⟨h₁, h₂, h₃⟩ := finishEdge_ranges (D := d + 1)
          (by rw [if_pos (hb₃.tree.2 ⟨_, _, _, rfl⟩)]) hi₃ hs₃ hnd hσ₃ hg₃ hb₃ hr₃
        exact ⟨h₁, h₂, lt_ne_of_g h₃ hσ₃⟩
      case pre =>
        exact fun w outs _ =>
          (hi₁.frame' (s' := { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne })).setSv w
end

/-- `walkTree` preserves the range invariant and advances `n` by the number of edges below `t`,
given the ear guards, the bookkeeping facts and the range-side hypotheses `RgTree`. -/
theorem walkTree_rangesInv (t : DfsTree) (d : Nat) (s : WalkState)
    (hi : ∀ v outs, t = .node v outs → ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).RangesInv σ n d)
    (hs : Shape s) (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne) (hg : GuardsTree t d s) (hb : BookTree t d s)
    (hr : RgTree σ n t d s) :
    ((walkTree t d).run s).2.RangesInv σ (n + t.edgePostorder.length) d ∧ Shape ((walkTree t d).run s).2 :=
  let r := rgTree σ n t d s hi hs hnd hσ hg hb hr
  ⟨r.1, r.2.1⟩

end WalkState
end Spqr
