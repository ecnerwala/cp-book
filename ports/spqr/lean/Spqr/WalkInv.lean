import Spqr.WalkSpec
import Spqr.Frame

/-!
# From the tstack guards to `FinishOk`

`finishOk_of_guards` assembles `WalkState.FinishOk` (the per-block hypotheses of `finishEdge_inv`)
from `FinishGuards`, the invariant, bookkeeping facts about the edge being finished, and the ear
facts `ear_*` below, which `FinishGuards` does not provide and are left to the ear invariant
(`EarShape`).
-/

namespace Spqr
open WalkM

namespace WalkState

variable {D : Nat} {s : WalkState}

theorem edgeBelow_vert_nil {v : Nat} (hv : v < s.g.nv) (hch : Items.ch s.items (vertItem v) = []) (e : Nat) :
    ¬ Items.EdgeBelow s.g s.items (vertItem v) e := by
  intro h
  rcases Relation.ReflTransGen.cases_head h with h | ⟨c, hc, -⟩
  · exact absurd h (by show (1 + v : Nat) ≠ 1 + s.g.nv + e; omega)
  · rw [Items.IsParent, hch] at hc; exact List.not_mem_nil hc

theorem finishTailOk_of_nil {curV d : Nat} {hasVert isSingle : Bool} (hv : curV < s.g.nv)
    (hch : hasVert = false → Items.ch s.items (vertItem curV) = [])
    (hm : hasVert = false → isSingle = false → MergeTopOk D (after (pushVertTstack curV d) s)) :
    FinishTailOk D curV d hasVert isSingle s :=
  ⟨fun h => Graph.ConnEdges.empty fun e _ => edgeBelow_vert_nil hv (hch h) e,
   fun h => Graph.TwoAttached.empty fun e _ => edgeBelow_vert_nil hv (hch h) e, hm⟩

section Ear

variable {curV d lv : Nat} {kind : RetKind} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}

/-! Ear facts (Invariant W / `EarShape`), one per `FinishOk` field that `FinishGuards` does not
cover. Each is stated at the state where the block runs. -/

/-- Loop 1: every iteration merges/unwraps/closes a finished sub-ear (`Loop1BodyOk`). -/
theorem ear_loop1 (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hi : s.Inv D) (hs : Shape s) (k : Nat)
    (hk : ∀ j, j ≤ k → result (loop1Cond d)
      (iter (Spqr.loop1Body d s.stackDir[d]!) j (ceS₁ o.dest d o.e (feS₀ d o s))) = true) :
    Loop1BodyOk D d s.stackDir[d]! (iter (Spqr.loop1Body d s.stackDir[d]!) k (ceS₁ o.dest d o.e (feS₀ d o s))) := by
  sorry

/-- Loop 2: every late merge joins entries sharing a terminal (`MergeTopOk`). -/
theorem ear_mergeLate (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hi : s.Inv D) (hs : Shape s) :
    MergeLateOk D d (feS₁ d o s) := by
  sorry

/-- The vertex close: loop 3 merges, the unwrap, the two merges, the retarget and the type-1 close. -/
theorem ear_closeVert (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hi : s.Inv D) (hs : Shape s) (hv : hasVert = true) :
    CloseVertOk D curV s.stackDir[d]! o.cls.isType1 origTstack (feSingle d o s) (feS₂ d o s) := by
  sorry

/-- The P-check after the vertex close. -/
theorem ear_finishP_vert (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hi : s.Inv D) (hs : Shape s) (hv : hasVert = true) :
    FinishPOk D curV lv o.cls.isType1 (feS₃ curV d o origTstack s) := by
  sorry

/-- The P-check of a first tree edge (no vertex entry yet). -/
theorem ear_finishP_tree (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hi : s.Inv D) (hs : Shape s) (hv : hasVert = false) :
    FinishPOk D curV lv o.cls.isType1 (feS₂ d o s) := by
  sorry

/-- The P-check of a back edge. -/
theorem ear_finishP_back (ho : o.cls = .ret lv kind) (hlow : lv < d) (hb : o.cls.isTree = false)
    (hg : FinishGuards d o origTstack hasVert s) (hi : s.Inv D) (hs : Shape s) :
    FinishPOk D curV lv o.cls.isType1 (feBack curV lv d o s) := by
  sorry

/-- The merge of the vertex entry into the type-2 first-edge entry. -/
theorem ear_tail_tree (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hi : s.Inv D) (hs : Shape s) (hv : hasVert = false)
    (hsingle : feSingle d o s = false) :
    MergeTopOk D (after (pushVertTstack curV d) (after (finishP curV lv o.cls.isType1) (feS₂ d o s))) := by
  sorry

/-! Bookkeeping frame facts through the blocks: `g` is constant and the vertex item of `curV` gains
no children before its entry is pushed (only node items' `ch` change). -/

theorem tail_frame_tree (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hi : s.Inv D) (hs : Shape s) (hv : hasVert = false) :
    (after (finishP curV lv o.cls.isType1) (feS₂ d o s)).g = s.g ∧
    (Items.ch s.items (vertItem curV) = [] →
      Items.ch (after (finishP curV lv o.cls.isType1) (feS₂ d o s)).items (vertItem curV) = []) := by
  sorry

theorem tail_frame_back (ho : o.cls = .ret lv kind) (hlow : lv < d) (hb : o.cls.isTree = false)
    (hg : FinishGuards d o origTstack hasVert s) (hi : s.Inv D) (hs : Shape s) :
    (after (finishP curV lv o.cls.isType1) (feBack curV lv d o s)).g = s.g ∧
    (Items.ch s.items (vertItem curV) = [] →
      Items.ch (after (finishP curV lv o.cls.isType1) (feBack curV lv d o s)).items (vertItem curV) = []) := by
  sorry

/-- `FinishOk` from the guards, the invariant, the bookkeeping facts of the finished edge
(`he`, `hq`, `hends`, `hch`) and the ear facts. -/
theorem finishOk_of_guards (ho : o.cls = .ret lv kind) (hlow : lv < d)
    (hg : FinishGuards d o origTstack hasVert s) (hi : s.Inv D) (hs : Shape s)
    (hdD : d ≤ D) (hv : curV < s.g.nv) (he : o.e < s.g.ne)
    (hq : Items.ch s.items (edgeItem s.g o.e) = [])
    (hends : Items.PairEq (if o.cls.isTree then (o.dest, s.stackVerts[d]!) else (curV, s.stackVerts[lv]!))
      s.g.edges[o.e]!)
    (hch : hasVert = false → Items.ch s.items (vertItem curV) = []) :
    FinishOk D curV d lv o origTstack hasVert s where
  e_lt := he
  ears ht := ⟨he,
    by show Items.ch (s.items.modify _ _) _ = []; rw [Items.ch_modify_ch_eq (edgeItem s.g o.e)
        (fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) (fun _ => rfl)]; exact hq,
    by show Items.PairEq (o.dest, s.stackVerts[d]!) s.g.edges[o.e]!; simpa [ht] using hends,
    hdD, ear_loop1 ho hlow ht hg hi hs⟩
  late ht := ear_mergeLate ho hlow ht hg hi hs
  vert ht hv' := ear_closeVert ho hlow ht hg hi hs hv'
  rest_vert ht hv' := ⟨ear_finishP_vert ho hlow ht hg hi hs hv',
    fun h => by simp [hv'] at h, fun h => by simp [hv'] at h, fun h => by simp [hv'] at h⟩
  rest_tree ht hv' :=
    have hf := tail_frame_tree (curV := curV) ho hlow ht hg hi hs hv'
    ⟨ear_finishP_tree ho hlow ht hg hi hs hv',
     finishTailOk_of_nil (by rw [hf.1]; exact hv) (fun _ => hf.2 (hch hv')) fun _ => ear_tail_tree ho hlow ht hg hi hs hv'⟩
  q _ := hq
  ends hb := by simpa [hb] using hends
  lv_le _ := by omega
  rest_back hb :=
    have hf := tail_frame_back (curV := curV) ho hlow hb hg hi hs
    ⟨ear_finishP_back ho hlow hb hg hi hs,
     finishTailOk_of_nil (by rw [hf.1]; exact hv) (fun h => hf.2 (hch h)) fun _ h => by cases h⟩

end Ear

/-! ## `walkTree_inv` by mutual induction

Hypotheses are supplied in the style of `GuardsTree`: the ear facts through `FinishGuards` and the
bookkeeping facts (vertex/edge bounds, fresh `Q`/`V` items, edge endpoints) through `BookTree`. -/

section Walk

variable {α : Type}

theorem Shape.frame' {s' : WalkState} (h : Shape s) (hg : s'.g = s.g := by rfl)
    (hi : s'.items = s.items := by rfl) (hts : s'.tstack = s.tstack := by rfl) : Shape s' :=
  h.frame hg hi hts

theorem Inv.frame' {D : Nat} {s' : WalkState} (h : s.Inv D) (hg : s'.g = s.g := by rfl)
    (hsv : s'.stackVerts = s.stackVerts := by rfl) (hi : s'.items = s.items := by rfl)
    (hts : s'.tstack = s.tstack := by rfl) : s'.Inv D :=
  h.frame hg hsv hi hts

theorem wp_of_forall {m : WalkM α} {Q : α → WalkState → Prop} (h : ∀ a s', Q a s') : wp m Q s := h _ _

theorem wp_imp {m : WalkM α} {P Q : α → WalkState → Prop} (h : wp m (fun a s' => P a s' → Q a s') s)
    (hp : wp m P s) : wp m Q s := h hp

theorem wp_and {m : WalkM α} {P Q : α → WalkState → Prop} (h₁ : wp m P s) (h₂ : wp m Q s) :
    wp m (fun a s' => P a s' ∧ Q a s') s := ⟨h₁, h₂⟩

theorem Inv.setSv {d : Nat} (x : Nat) (h : s.Inv d) :
    ({ s with stackVerts := s.stackVerts.set! (d + 1) x } : WalkState).Inv (d + 1) := by
  refine ⟨fun t ht => ?_, fun i hi hsz => ItemInv.congr (s := s) rfl rfl (fun _ _ => Iff.rfl) (h.nodes i hi hsz)⟩
  obtain ⟨conn, att⟩ := h.entries t ht
  refine ⟨conn, att.mono fun w hw => ?_⟩
  rcases hw with hw | ⟨k, h1, h2, h3⟩
  · exact .inl hw
  · refine .inr ⟨k, h1, by omega, ?_⟩
    have hk : k ≠ d + 1 := by omega
    have : (s.stackVerts.set! (d + 1) x)[k]! = s.stackVerts[k]! := by
      simp [Array.set!, getElem!_def, Array.getElem?_setIfInBounds, hk, Ne.symm hk]
    rw [h3]; exact this.symm

theorem ret_of_lowval_lt {o : DfsOut} {d : Nat} (h : o.cls.lowval d < d) :
    ∃ lv kind, o.cls = .ret lv kind ∧ lv < d := by
  cases hc : o.cls <;> simp only [OutClass.lowval, hc] at h
  all_goals first | omega | exact ⟨_, _, rfl, h⟩

/-- Bookkeeping facts needed by `finishEdge_inv` at the state where `finishEdge` runs. -/
structure FinishBook (curV d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop where
  v_lt : curV < s.g.nv
  e_lt : o.e < s.g.ne
  q : Items.ch s.items (edgeItem s.g o.e) = []
  ends : ∀ lv kind, o.cls = .ret lv kind →
    Items.PairEq (if o.cls.isTree then (o.dest, s.stackVerts[d]!) else (curV, s.stackVerts[lv]!))
      s.g.edges[o.e]!
  vch : hasVert = false → Items.ch s.items (vertItem curV) = []
  tree : o.cls.isTree = true ↔ ∃ e cls c, o = .tree e cls c

/-- Before the vertex entry of `v` is pushed, `v` is in range and its item has no children. -/
def VertBook (v : Nat) (hasVert : Bool) (s : WalkState) : Prop :=
  hasVert = false → v < s.g.nv ∧ Items.ch s.items (vertItem v) = []

mutual
def BookTree (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  match t with
  | .node v outs => BookOuts v d outs false { s with stackVerts := s.stackVerts.set! d v }

def BookOuts (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  match outs with
  | [] => VertBook v hasVert s
  | o :: rest => BookOut v d o hasVert s ∧
      wp (walkOut v d o hasVert) (fun hasVert' s' => BookOuts v d rest hasVert' s') s

def BookOut (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  VertBook v hasVert s ∧
  wp (walkOutPre v d o hasVert) (fun hasVert' s₁ =>
    match o with
    | .tree _ _ child =>
      wp (modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }) (fun _ s₂ =>
        BookTree child (d + 1) s₂ ∧
        wp (walkTree child (d + 1)) (fun _ s₃ => FinishBook v d o hasVert' s₃) s₂) s₁
    | .back .. => FinishBook v d o hasVert' s₁) s
end

/-- Bookkeeping of the walk, to be discharged from the DFS well-formedness and endpoint facts
(`DfsTree.WF`, the `dfsForest` endpoint lemma) and the fresh-items start state; cf. `walkTree_guards`. -/
theorem walkTree_book (t : DfsTree) (d : Nat) (s : WalkState)
    (hfresh : ∀ i, i < 1 + s.g.nv + s.g.ne → Items.ch s.items i = []) :
    BookTree t d s := by
  sorry

/-- Ear fact: after `finishEdge` of a returning tree edge at depth `d`, no open entry is attached at
`stackVerts[d+1]`, so `Inv (d+1)` lowers to `Inv d`. -/
theorem ear_lower {curV d lv : Nat} {kind : RetKind} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}
    (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hi : s.Inv (d + 1)) (hs : Shape s)
    (hi' : (after (finishEdge curV d o origTstack hasVert) s).Inv (d + 1)) :
    (after (finishEdge curV d o origTstack hasVert) s).Inv d := by
  sorry

/-- Boundary edges (`lowval ≥ d`: bridges, components, self-loops) close a block via
`finishBoundary`; the invariant is kept. -/
theorem finishBoundary_inv {curV d : Nat} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}
    (hge : d ≤ o.cls.lowval d) (hi : s.Inv (d + 1)) (hs : Shape s) (hb : FinishBook curV d o hasVert s) :
    Step d curV s (after (finishEdge curV d o origTstack hasVert) s) := by
  sorry

theorem walkOutPre_inv {v d : Nat} {o : DfsOut} {hasVert : Bool} (hi : s.Inv d) (hs : Shape s)
    (hb : VertBook v hasVert s) :
    wp (walkOutPre v d o hasVert) (fun _ s₁ => s₁.Inv d ∧ Shape s₁) s := by
  have hi' : ({ s with stackDir := s.stackDir.set! d (if o.cls.lowval d ≥ d then false else !s.stackDir[o.cls.lowval d]!) } : WalkState).Inv d :=
    hi.frame'
  have hs' : Shape ({ s with stackDir := s.stackDir.set! d (if o.cls.lowval d ≥ d then false else !s.stackDir[o.cls.lowval d]!) } : WalkState) :=
    hs.frame'
  unfold walkOutPre
  simp only [wp_bind, wp_stackDir, wp_setStackDir, wp_ite]
  split
  · rename_i hc
    have hf : hasVert = false := by cases hasVert <;> simp_all
    obtain ⟨hv, hch⟩ := hb hf
    exact ⟨(Step.pushVert hi' hs' d hv (Graph.ConnEdges.empty fun e _ => edgeBelow_vert_nil hv hch e)
        (Graph.TwoAttached.empty fun e _ => edgeBelow_vert_nil hv hch e)).inv,
      (Step.pushVert hi' hs' d hv (Graph.ConnEdges.empty fun e _ => edgeBelow_vert_nil hv hch e)
        (Graph.TwoAttached.empty fun e _ => edgeBelow_vert_nil hv hch e)).shape⟩
  · exact ⟨hi', hs'⟩

theorem finishEdge_step {v d : Nat} {o : DfsOut} {n : Nat} {hasVert : Bool} {D : Nat}
    (hD : D = if o.cls.isTree then d + 1 else d) (hi : s.Inv D) (hs : Shape s)
    (hg : FinishGuards d o n hasVert s) (hb : FinishBook v d o hasVert s) :
    (after (finishEdge v d o n hasVert) s).Inv d ∧ Shape (after (finishEdge v d o n hasVert) s) := by
  have hdD : d ≤ D := by split at hD <;> omega
  have hD' : D ≤ d + 1 := by split at hD <;> omega
  by_cases hge : d ≤ o.cls.lowval d
  · have hst := finishBoundary_inv (origTstack := n) hge (hi.mono hD') hs hb
    exact ⟨hst.inv, hst.shape⟩
  · obtain ⟨lv, kind, ho, hl⟩ := ret_of_lowval_lt (Nat.lt_of_not_le hge)
    have hst := finishEdge_inv v d lv kind o n hasVert ho hl hb.v_lt hi hs
      (finishOk_of_guards ho hl hg hi hs hdD hb.v_lt hb.e_lt hb.q (hb.ends lv kind ho) hb.vch)
    refine ⟨?_, hst.shape⟩
    by_cases ht : o.cls.isTree = true
    · rw [if_pos ht] at hD; subst hD
      exact ear_lower ho hl ht hg hi hs hst.inv
    · rw [if_neg ht] at hD; subst hD; exact hst.inv

abbrev InvTree (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  (∀ v outs, t = .node v outs → ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).Inv d) →
  Shape s → GuardsTree t d s → BookTree t d s → wp (walkTree t d) (fun _ s' => s'.Inv d ∧ Shape s') s

abbrev InvOuts (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  s.Inv d → Shape s → GuardsOuts v d outs hasVert s → BookOuts v d outs hasVert s →
  wp (walkOuts v d outs hasVert) (fun hasVert' s' => s'.Inv d ∧ Shape s' ∧ VertBook v hasVert' s') s

abbrev InvOut (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  s.Inv d → Shape s → GuardsOut v d o hasVert s → BookOut v d o hasVert s →
  wp (walkOut v d o hasVert) (fun _ s' => s'.Inv d ∧ Shape s') s

mutual
theorem invTree : ∀ (t : DfsTree) (d : Nat) (s : WalkState), InvTree t d s
  | .node v outs, d, s => fun hi hs hg hb => by
    unfold walkTree
    simp only [wp_bind, wp_modify]
    unfold GuardsTree at hg; unfold BookTree at hb
    refine wp_imp (wp_of_forall fun hv s' ⟨hi', hs', hvb⟩ => ?_)
      (invOuts v d outs false _ (hi v outs rfl) hs.frame' hg hb)
    cases hv
    · simp only [Bool.false_eq_true, ↓reduceIte, wp_bind, wp_setStackDir]
      obtain ⟨hv, hch⟩ := hvb rfl
      have hi₂ : ({ s' with stackDir := s'.stackDir.set! d true } : WalkState).Inv d := hi'.frame'
      have hs₂ : Shape ({ s' with stackDir := s'.stackDir.set! d true } : WalkState) := hs'.frame'
      exact ⟨(Step.pushVert hi₂ hs₂ d hv (Graph.ConnEdges.empty fun e _ => edgeBelow_vert_nil hv hch e)
          (Graph.TwoAttached.empty fun e _ => edgeBelow_vert_nil hv hch e)).inv,
        (Step.pushVert hi₂ hs₂ d hv (Graph.ConnEdges.empty fun e _ => edgeBelow_vert_nil hv hch e)
          (Graph.TwoAttached.empty fun e _ => edgeBelow_vert_nil hv hch e)).shape⟩
    · exact ⟨hi', hs'⟩

theorem invOuts : ∀ (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState),
    InvOuts v d outs hasVert s
  | v, d, [], hasVert, s => fun hi hs _ hb => by
    unfold BookOuts at hb
    unfold walkOuts
    simp only [wp_pure]
    exact ⟨hi, hs, hb⟩
  | v, d, o :: rest, hasVert, s => fun hi hs hg hb => by
    unfold GuardsOuts at hg; unfold BookOuts at hb
    unfold walkOuts
    simp only [wp_bind]
    exact wp_imp (wp_imp (wp_imp (wp_of_forall fun hv' s' ⟨hi', hs'⟩ hg' hb' =>
      invOuts v d rest hv' s' hi' hs' hg' hb') (invOut v d o hasVert s hi hs hg.1 hb.1)) hg.2) hb.2

theorem invOut : ∀ (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState), InvOut v d o hasVert s
  | v, d, o, hasVert, s => fun hi hs hg hb => by
    rw [walkOut_eq, wp_bind]
    unfold GuardsOut at hg; unfold BookOut at hb
    refine wp_imp (wp_imp (wp_imp (wp_of_forall fun hv' s₁ ⟨hi₁, hs₁⟩ hg₁ hb₁ => ?_)
      (walkOutPre_inv hi hs hb.1)) hg) hb.2
    unfold walkOutRest
    rw [wp_bind, wp_tstackSize]
    cases o with
    | back e cls dest =>
      try simp only [wp_pure] at hg₁ hb₁
      simp only [wp_bind, wp_pure]
      exact finishEdge_step (D := d)
        (by rw [if_neg fun h => by obtain ⟨_, _, _, h⟩ := hb₁.tree.1 h; cases h]) hi₁ hs₁ hg₁ hb₁
    | tree e cls child =>
      try simp only [wp_bind, wp_modify] at hg₁ hb₁
      simp only [wp_bind, wp_modify]
      refine wp_imp (wp_imp (wp_of_forall fun _ s₃ ⟨hg₃, hb₃⟩ ⟨hi₃, hs₃⟩ => ?body) (wp_and hg₁.2 hb₁.2))
        (invTree child (d + 1) _ ?pre hs₁.frame' hg₁.1 hb₁.1)
      case body =>
        exact finishEdge_step (D := d + 1) (by rw [if_pos (hb₃.tree.2 ⟨_, _, _, rfl⟩)]) hi₃ hs₃ hg₃ hb₃
      case pre => exact fun w outs _ => (hi₁.frame' (s' := { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne })).setSv w
end

/-- `walkTree` preserves the invariant, given the ear guards and the bookkeeping facts. -/
theorem walkTree_inv' (t : DfsTree) (d : Nat) (s : WalkState)
    (hi : ∀ v outs, t = .node v outs → ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).Inv d)
    (hs : Shape s) (hg : GuardsTree t d s) (hb : BookTree t d s) :
    ((walkTree t d).run s).2.Inv d ∧ Shape ((walkTree t d).run s).2 :=
  invTree t d s hi hs hg hb

end Walk

end WalkState
end Spqr
