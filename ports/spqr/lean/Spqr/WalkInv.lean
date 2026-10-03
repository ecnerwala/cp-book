import Spqr.WalkSpec
import Spqr.Frame
import Spqr.EarInv
import Spqr.EarLoop1
import Spqr.EarLoop2
import Spqr.Proofs.Dfs
import Spqr.EarFrame
import Spqr.EarCtxAt
import Spqr.StFrame

/-!
# From the tstack guards to `FinishOk`

`finishOk_of_guards` assembles `WalkState.FinishOk` (the per-block hypotheses of `finishEdge_inv`)
from `FinishGuards`, the invariant, bookkeeping facts about the edge being finished, and the ear
facts `ear_*` below, which `FinishGuards` does not provide and are left to the ear invariant
(`EarShape`).
-/

namespace Spqr

/-! Lowval / return-depth facts about a well-formed out-edge (the class of a tree out is determined
by the return depths of its subtree). -/
namespace EarDfs

theorem lowval_classify_tree {d : Nat} {n : Lowvals} (h : n.1 ≤ d + 1) :
    (classify d true n).lowval d = n.1 := by
  unfold classify
  split
  · next h1 =>
    split
    · next h2 => simp at h2; simp [h2, OutClass.lowval]
    · next h2 => simp at h1 h2; simp [OutClass.lowval]; omega
  · rfl

/-- The lowval of a well-formed out-edge at depth `d` is the least depth its piece returns to
(`d + 1` for a bridge). -/
theorem lowval_eq_lmin {anc : List Nat} {v : Nat} {o : DfsOut} (hwf : o.WF anc v) :
    o.cls.lowval anc.length = lmin (anc.length + 1) (o.retDepths anc.length) := by
  cases o with
  | back e dest cls =>
    rw [DfsOut.WF] at hwf
    obtain ⟨i, hi, rfl⟩ := hwf
    have hi' : i ≤ anc.length := by
      have := (List.getElem?_eq_some_iff.1 hi).1; simp at this; omega
    simp only [DfsOut.cls, DfsOut.retDepths, lmin_cons, lmin_nil, lowval_classify_back hi']
    omega
  | tree e cls child =>
    rw [DfsOut.WF] at hwf
    obtain ⟨-, rfl⟩ := hwf
    simp only [DfsOut.cls, DfsOut.retDepths]
    exact lowval_classify_tree (lmin_le _ _)

theorem lowval_mem_retDepths {anc : List Nat} {v : Nat} {o : DfsOut} (hwf : o.WF anc v)
    (hret : o.cls.lowval anc.length < anc.length) :
    o.cls.lowval anc.length ∈ o.retDepths anc.length := by
  rw [lowval_eq_lmin hwf] at hret ⊢
  rcases lmin_eq_or_mem (anc.length + 1) (o.retDepths anc.length) with h | h
  · omega
  · exact h

theorem retDepths_sub {d y : Nat} {outs : List DfsOut} {o : DfsOut} (ho : o ∈ outs) :
    ∀ x ∈ o.retDepths d, x ∈ (DfsTree.node y outs).retDepths d := by
  intro x hx
  show x ∈ DfsOut.retDepthsList d outs
  rw [DfsOut.retDepthsList_eq]
  exact List.mem_flatMap.2 ⟨o, ho, hx⟩

/-- A bridge's child has only boundary outs. -/
theorem bridge_outs_bd {anc : List Nat} {v y : Nat} {cls : OutClass} {outs : List DfsOut}
    (hwf : (DfsTree.node y outs).WF (anc ++ [v]))
    (hcls : cls = classify anc.length true
      (low2 (anc.length + 1) ((DfsTree.node y outs).retDepths (anc.length + 1))))
    (hb : cls = .bridge) : ∀ o ∈ outs, anc.length + 1 ≤ o.cls.lowval (anc.length + 1) := by
  intro o ho
  rw [DfsTree.WF] at hwf
  have hwo := hwf.2 o ho
  have hlen : (anc ++ [v]).length = anc.length + 1 := by simp
  rw [hcls, classify_eq_bridge_iff, low2_fst] at hb
  by_contra hlt
  rw [Nat.not_le] at hlt
  have hm := lowval_mem_retDepths hwo
  rw [hlen] at hm
  have := lmin_le_of_mem (d := anc.length + 1) (retDepths_sub (y := y) ho _ (hm hlt))
  omega

/-- A component's child returns exactly to the parent depth. -/
theorem comp_outs_ret {anc : List Nat} {v y : Nat} {cls : OutClass} {outs : List DfsOut}
    (hwf : (DfsTree.node y outs).WF (anc ++ [v]))
    (hcls : cls = classify anc.length true
      (low2 (anc.length + 1) ((DfsTree.node y outs).retDepths (anc.length + 1))))
    (hc : cls = .component) : ∀ o ∈ outs, o.cls.lowval (anc.length + 1) < anc.length + 1 →
      o.cls.lowval (anc.length + 1) = anc.length := by
  intro o ho hlt
  rw [DfsTree.WF] at hwf
  have hwo := hwf.2 o ho
  have hlen : (anc ++ [v]).length = anc.length + 1 := by simp
  rw [hcls, classify_eq_component_iff, low2_fst] at hc
  have hm := lowval_mem_retDepths hwo
  rw [hlen] at hm
  have := lmin_le_of_mem (d := anc.length + 1) (retDepths_sub (y := y) ho _ (hm hlt))
  omega

theorem comp_outs_ex {anc : List Nat} {v y : Nat} {cls : OutClass} {outs : List DfsOut}
    (hwf : (DfsTree.node y outs).WF (anc ++ [v]))
    (hcls : cls = classify anc.length true
      (low2 (anc.length + 1) ((DfsTree.node y outs).retDepths (anc.length + 1))))
    (hc : cls = .component) : ∃ o ∈ outs, o.cls.lowval (anc.length + 1) < anc.length + 1 := by
  rw [DfsTree.WF] at hwf
  have hlen : (anc ++ [v]).length = anc.length + 1 := by simp
  rw [hcls, classify_eq_component_iff, low2_fst] at hc
  rcases lmin_eq_or_mem (anc.length + 1) ((DfsTree.node y outs).retDepths (anc.length + 1))
    with h | h
  · omega
  · rw [hc] at h
    change anc.length ∈ DfsOut.retDepthsList _ outs at h
    rw [DfsOut.retDepthsList_eq, List.mem_flatMap] at h
    obtain ⟨o, ho, hm⟩ := h
    refine ⟨o, ho, ?_⟩
    have := lowval_eq_lmin (hwf.2 o ho)
    rw [hlen] at this
    rw [this]
    have := lmin_le_of_mem (d := anc.length + 1 + 1) hm
    omega

theorem comp_outs_t1 {anc : List Nat} {v y : Nat} {cls : OutClass} {outs : List DfsOut}
    (hwf : (DfsTree.node y outs).WF (anc ++ [v]))
    (hcls : cls = classify anc.length true
      (low2 (anc.length + 1) ((DfsTree.node y outs).retDepths (anc.length + 1))))
    (hc : cls = .component) : ∀ o ∈ outs, o.cls.isType1 = true := by
  intro o ho
  rw [DfsTree.WF] at hwf
  have hwo := hwf.2 o ho
  have hlen : (anc ++ [v]).length = anc.length + 1 := by simp
  rw [hcls, classify_eq_component_iff, low2_fst] at hc
  cases o with
  | back e dest cls' =>
    rw [DfsOut.WF] at hwo
    obtain ⟨i, -, rfl⟩ := hwo
    simp only [DfsOut.cls, classify]
    split
    · split <;> rfl
    · rfl
  | tree e cls' child =>
    rw [DfsOut.WF] at hwo
    obtain ⟨-, rfl⟩ := hwo
    simp only [DfsOut.cls, hlen]
    by_contra hne
    have hk : ∃ l, classify (anc.length + 1) true
        (low2 (anc.length + 1 + 1) (child.retDepths (anc.length + 1 + 1))) = .ret l .type2Child := by
      generalize classify (anc.length + 1) true
        (low2 (anc.length + 1 + 1) (child.retDepths (anc.length + 1 + 1))) = c at hne
      cases c with
      | ret l k => cases k <;> simp_all [OutClass.isType1]
      | _ => simp [OutClass.isType1] at hne
    obtain ⟨l, hl⟩ := hk
    rw [classify_eq_ret_iff_tree] at hl
    obtain ⟨h1, hld, h2⟩ := hl
    have h2' : (low2 (anc.length + 1 + 1) (child.retDepths (anc.length + 1 + 1))).2 < anc.length + 1 := by
      by_contra hge
      rw [ite_eq_right hge] at h2
      cases h2
    rw [low2_snd] at h2'
    rw [low2_fst] at h1
    rcases lmin_eq_or_mem (anc.length + 1 + 1)
      ((child.retDepths (anc.length + 1 + 1)).filter (· ≠ lmin (anc.length + 1 + 1) (child.retDepths (anc.length + 1 + 1))))
      with h | h
    · omega
    · rw [List.mem_filter] at h
      obtain ⟨hm, hne'⟩ := h
      have hsub := retDepths_sub (d := anc.length + 1) (y := y) ho _ hm
      have hge := lmin_le_of_mem (d := anc.length + 1) hsub
      rw [hc] at hge
      have hl' : lmin (anc.length + 1 + 1) (child.retDepths (anc.length + 1 + 1)) = l := h1
      have hlm : l ∈ child.retDepths (anc.length + 1 + 1) := by
        rcases lmin_eq_or_mem (anc.length + 1 + 1) (child.retDepths (anc.length + 1 + 1)) with h3 | h3
        · omega
        · rwa [hl'] at h3
      have := lmin_le_of_mem (d := anc.length + 1) (retDepths_sub (d := anc.length + 1) (y := y) ho _ hlm)
      rw [hc] at this
      have hne2 := of_decide_eq_true hne'
      omega

end EarDfs

open WalkM

namespace WalkState

variable {D : Nat} {s : WalkState}

theorem edgeBelow_vert_nil {v : Nat} (hv : v < s.g.nv) (hch : Items.ch s.items (vertItem v) = []) (e : Nat) :
    ¬ Items.EdgeBelow s.g s.items (vertItem v) e := by
  intro h
  rcases Relation.ReflTransGen.cases_head h with h | ⟨c, hc, -⟩
  · exact absurd h (by show (1 + v : Nat) ≠ 1 + s.g.nv + e; omega)
  · rw [Items.IsParent, hch] at hc; exact List.not_mem_nil hc

theorem finishTailOk_of_vert {curV d : Nat} {hasVert isSingle : Bool}
    (hc : hasVert = false → s.g.ConnEdges (Items.EdgeBelow s.g s.items (vertItem curV)))
    (ha : hasVert = false → s.g.TwoAttached (Items.EdgeBelow s.g s.items (vertItem curV)) curV curV)
    (hm : hasVert = false → isSingle = false → MergeTopOk D (after (pushVertTstack curV d) s)) :
    FinishTailOk D curV d hasVert isSingle s :=
  ⟨hc, ha, hm⟩

/-- A `Step` keeps the graph and the subtree of the vertex item, hence its connectivity and
2-attachment. -/
theorem Step.vertTransport {v : Nat} {s' : WalkState} (st : Step D v s s')
    (hc : s.g.ConnEdges (Items.EdgeBelow s.g s.items (vertItem v)))
    (ha : s.g.TwoAttached (Items.EdgeBelow s.g s.items (vertItem v)) v v) :
    s'.g.ConnEdges (Items.EdgeBelow s'.g s'.items (vertItem v)) ∧
    s'.g.TwoAttached (Items.EdgeBelow s'.g s'.items (vertItem v)) v v := by
  rw [st.g]
  have hE : ∀ e, e < s.g.ne →
      (Items.EdgeBelow s.g s'.items (vertItem v) e ↔ Items.EdgeBelow s.g s.items (vertItem v) e) :=
    fun e _ => st.below _
  exact ⟨(Graph.ConnEdges.congr hE).2 hc, (Graph.TwoAttached.congr hE).2 ha⟩

section Ear

variable {curV d lv : Nat} {kind : RetKind} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}
  {sub base : List TEntry}

/-! Ear facts (Invariant W / `EarShape`), one per `FinishOk` field that `FinishGuards` does not
cover. Each is stated at the state where the block runs. -/

/-- Loop 1: every iteration merges/unwraps/closes a finished sub-ear (`Loop1BodyOk`). -/
theorem ear_loop1 (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hE : s.EarFinish curV d o hasVert sub base) (hi : s.Inv' D) (hs : Shape s)
    (hD : D = d + 1) (hv : curV < s.g.nv) (he : o.e < s.g.ne) (hq : Items.ch s.items (edgeItem s.g o.e) = [])
    (hends : Items.PairEq (o.dest, s.stackVerts[d]!) s.g.edges[o.e]!) (k : Nat)
    (hk : ∀ j, j ≤ k → result (loop1Cond d)
      (iter (Spqr.loop1Body d s.stackDir[d]!) j (ceS₁ o.dest d o.e (feS₀ d o s))) = true) :
    Loop1BodyOk D d s.stackDir[d]! (iter (Spqr.loop1Body d s.stackDir[d]!) k (ceS₁ o.dest d o.e (feS₀ d o s))) := by
  have hlow' : o.cls.lowval d < d := by rw [ho]; exact hlow
  obtain ⟨hi', lo, hsplit, hc⟩ := L1Ctx.ofEar hE hs hD ht hlow'
  exact loop1_ok hc (l1_init hE hi hs hD he hq hends hsplit) hv k hk

/-- Loop 2: every late merge joins entries sharing a terminal (`MergeTopOk`). -/
theorem ear_mergeLate (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hE : s.EarFinish curV d o hasVert sub base) (hi : s.Inv' D) (hs : Shape s) (hD : D = d + 1) :
    MergeLateOk D d (feS₁ d o s) := by
  obtain ⟨c₀, R, hL⟩ := hE.late ht (by rw [ho]; exact hlow)
  exact hD ▸ mergeLateOk_of_late hL

/-- The vertex close: loop 3 merges, the unwrap, the two merges, the retarget and the type-1 close. -/
theorem ear_closeVert (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hE : s.EarFinish curV d o hasVert sub base) (hi : s.Inv' D) (hs : Shape s) (hv : hasVert = true)
    (hD : D = d + 1) (hlen : base.length = origTstack) (hs₂ : Shape (feS₂ d o s)) :
    CloseVertOk D curV s.stackDir[d]! o.cls.isType1 origTstack (feSingle d o s) (feS₂ d o s) := by
  have hlow' : o.cls.lowval d < d := by rw [ho]; exact hlow
  obtain ⟨c, mid, py, vy, hC⟩ := hE.close ht hlow'
  subst hv hD hlen
  exact closeVertOk_of_close _ hC hs₂ hlow' hE.path
    (fun k hk => by rw [hE.sv_child ht]; exact hE.path_child ht k hk) (hE.dir_d hlow')

/-- The P-check after the vertex close. -/
theorem ear_finishP_vert (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hE : s.EarFinish curV d o hasVert sub base) (hi : s.Inv' D) (hs : Shape s) (hv : hasVert = true)
    (hD : D = d + 1) (hlen : base.length = origTstack) (hi₂ : (feS₂ d o s).Inv' D) (hs₂ : Shape (feS₂ d o s))
    (hs₃ : Shape (feS₃ curV d o origTstack s)) (he : o.e < s.g.ne)
    (hends : Items.PairEq (o.dest, s.stackVerts[d]!) s.g.edges[o.e]!) :
    FinishPOk D curV lv o.cls.isType1 (feS₃ curV d o origTstack s) := by
  have hlv : o.cls.lowval d = lv := by rw [ho]; rfl
  cases ht1 : o.cls.isType1
  · exact ⟨fun h => by simp [result, run_condP] at h⟩
  obtain ⟨c, mid, py, vy, hC⟩ := hE.close ht (by rw [hlv]; exact hlow)
  subst hv hD hlen
  unfold feS₃ at hs₃ ⊢
  rw [ht1] at hs₃ ⊢
  exact hlv ▸ finishPOk_type1_of_close _ hE hC rfl hi₂ hs₂ hs₃ ht ht1 (by rw [hlv]; exact hlow) he hends

/-- The P-check of a first tree edge (no vertex entry yet). -/
theorem ear_condP_tree (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hE : s.EarFinish curV d o hasVert sub base) (hi : s.Inv' D) (hs : Shape s) (hv : hasVert = false) :
    result (condP curV lv o.cls.isType1) (feS₂ d o s) = false := by
  refine Bool.eq_false_iff.2 fun h => ?_
  have hlow' : o.cls.lowval d < d := by rw [ho]; exact hlow
  obtain ⟨c, mid, py, vy, hts, ⟨mid₀, hsub⟩, -, -, h1⟩ := hE.loops ht hlow'
  have h' : (o.cls.isType1 && decide ((feS₂ d o s).tstack.length ≥ 2) &&
      ((feS₂ d o s).tstack.tail.head!.vStart == curV) && ((feS₂ d o s).tstack.tail.head!.topDepth == lv)) = true := h
  simp only [Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq] at h'
  obtain ⟨⟨⟨ht1, -⟩, hvs⟩, -⟩ := h'
  obtain ⟨hmid, -⟩ := h1 ht1
  subst hmid
  rw [hts] at hvs
  simp only [List.nil_append, List.cons_append, List.tail_cons, List.head!_cons] at hvs
  exact hE.sub_bot py (by rw [hsub]; simp) hvs

theorem ear_finishP_tree (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hE : s.EarFinish curV d o hasVert sub base) (hi : s.Inv' D) (hs : Shape s) (hv : hasVert = false) :
    FinishPOk D curV lv o.cls.isType1 (feS₂ d o s) :=
  ⟨fun h => by rw [ear_condP_tree ho hlow ht hg hE hi hs hv] at h; cases h⟩

/-- A back edge: the P-check merges the fresh `(curV, lv)` edge entry into the `(curV, lv)` entry of
`base` (`EarFinish.p_entry`); the result is one-sided on `stackDir[lv]`, attached at `curV`,
`stackVerts[lv]` and interior vertices only. -/
theorem ear_finishP_back (ho : o.cls = .ret lv kind) (hlow : lv < d) (hb : o.cls.isTree = false)
    (hg : FinishGuards d o origTstack hasVert s) (hE : s.EarFinish curV d o hasVert sub base) (hi : s.Inv' D) (hs : Shape s)
    (hD : D = d) (he : o.e < s.g.ne) (hq : Items.ch s.items (edgeItem s.g o.e) = [])
    (hends : Items.PairEq (curV, s.stackVerts[lv]!) s.g.edges[o.e]!)
    (hsB : Shape (feBack curV lv d o s)) :
    FinishPOk D curV lv o.cls.isType1 (feBack curV lv d o s) := by
  have hlv : o.cls.lowval d = lv := by rw [ho]; rfl
  set q := edgeItem s.g o.e with hqdef
  set f : Item → Item := fun it =>
    { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) } with hf
  set dir := s.stackDir[lv]! with hdir
  set c : TEntry := ⟨curV, lv, s.nxtEdgeIdx, setSides dir [q] []⟩ with hc
  set s₂ := feBack curV lv d o s with hs₂
  have hts₂ : s₂.tstack = c :: s.tstack := rfl
  have hit₂ : s₂.items = s.items.modify q f := rfl
  refine ⟨fun hcond => ?_⟩
  have hcond' : (o.cls.isType1 && decide (s₂.tstack.length ≥ 2) && (s₂.tstack.tail.head!.vStart == curV) &&
      (s₂.tstack.tail.head!.topDepth == lv)) = true := hcond
  simp only [Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq] at hcond'
  obtain ⟨⟨⟨ht1, hlen⟩, hbv⟩, hbt⟩ := hcond'
  rw [hts₂] at hlen hbv hbt
  obtain ⟨b, rest, hts⟩ : ∃ b rest, s.tstack = b :: rest := by
    match h : s.tstack with
    | [] => rw [h] at hlen; simp at hlen
    | b :: rest => exact ⟨b, rest, rfl⟩
  rw [hts] at hbv hbt
  simp only [List.tail_cons, List.head!_cons] at hbv hbt
  have hsub : sub = [] := hE.back_nil hb
  have hbase : base = s.tstack := by rw [hE.tstack, hsub]; rfl
  have hbmem : b ∈ s.tstack := by rw [hts]; exact List.mem_cons_self
  have hbbase : b ∈ base := by rw [hbase]; exact hbmem
  obtain ⟨⟨i, hbsp, hroot⟩, hatt⟩ :=
    hE.p_entry (by rw [hlv]; exact hlow) ht1 b hbbase hbv (by rw [hbt, hlv])
  rw [hbt] at hbsp hatt
  have hbside : getSide b.spans (!dir) = [] := by rw [hbsp]; exact getSide_setSides_not _ _
  have hbsingle : getSide b.spans dir = [i] := by rw [hbsp]; exact getSide_setSides_self _ _ _
  have hib : i ∈ b.spans.1 ++ b.spans.2 := by
    rw [mem_of_getSide_nil dir b.spans hbside, hbsingle]; exact List.mem_singleton_self i
  have hiq : i ≠ q := fun h => hE.q_free b hbmem (by rw [show edgeItem s.g o.e = i from h.symm]; exact hib)
  have hchq : ∀ p, Items.ch (s.items.modify q f) p = Items.ch s.items p :=
    Items.ch_modify_ch_eq q f fun _ => rfl
  have hE₁ : ∀ t ∈ s.tstack, ∀ e, t.edges s.g (s.items.modify q f) e ↔ t.edges s.g s.items e :=
    fun t ht e => TEntry.edges_modify_of_not_mem q f hE.q_root (hE.q_free t ht) e
  have hq₁ : Items.ch (s.items.modify q f) q = [] := by rw [hchq]; exact hq
  have hcE : ∀ e', c.edges s.g (s.items.modify q f) e' ↔ e' = o.e :=
    TEntry.edges_edgeEntry dir curV lv s.nxtEdgeIdx o.e hq₁
  have hts₂' : s₂.tstack = c :: b :: rest := by rw [hts₂, hts]
  have hpw := hE.disj
  have hpws := hE.span_disj
  rw [hts] at hpw hpws
  -- the unwrap
  have hu : UnwrapOk .P s₂ := by
    refine ⟨by rw [hts₂']; simp, fun _ _ => ?_⟩
    have hn : nxtE s₂ = b := by rw [nxtE, hts₂']; rfl
    have hd : nxtDir s₂ = dir := by rw [nxtDir, hn, hbt]; rfl
    have hh : nxtHead s₂ = i := by rw [nxtHead, hn, hd, hbsingle]; rfl
    refine ⟨?_, ?_, ?_, ?_⟩
    · rw [hn, hd]; exact hbside
    · rw [hn, hd, hh]; exact hbsingle
    · rw [hh]; intro p hp; exact hroot p ((Items.IsParent_congr hchq).1 hp)
    · rw [hh]; intro t ht hmem
      have hcur : curE s₂ = c := by rw [curE, hts₂']; rfl
      rw [hcur, hts₂'] at ht
      simp only [List.tail_cons, List.mem_cons] at ht
      rcases ht with rfl | ht
      · exact hiq (List.mem_singleton.1 ((mem_setSides dir [q] i).1 hmem))
      · exact (List.pairwise_cons.1 hpws).1 t ht i hib hmem
  -- the close
  obtain ⟨b', hts₃, hg₃, hsv₃, hsd₃, hbv', hbt', hside', hbE', hallE, -⟩ :=
    maybeUnwrapNxt_edges (s := s₂) (ty := .P) hsB (by decide) hts₂' (i := i)
      (by rw [hbt]; exact hbside) (by rw [hbt]; exact hbsingle)
  set s₃ := after (maybeUnwrapNxt .P) s₂ with hs₃
  have hg₃' : s₃.g = s.g := hg₃
  have hsv₃' : s₃.stackVerts = s.stackVerts := hsv₃
  have hsd₃' : s₃.stackDir = s.stackDir := hsd₃
  rw [hbt] at hside'
  have hbE₃ : ∀ e, e < s.g.ne → (b'.edges s.g s₃.items e ↔ b.edges s.g s.items e) :=
    fun e he => (hbE' e he).trans (hE₁ b hbmem e)
  have hcE₃ : ∀ e, c.edges s.g s₃.items e ↔ e = o.e := fun e => (hallE c e).trans (hcE e)
  have hrE₃ : ∀ t ∈ rest, ∀ e, t.edges s.g s₃.items e ↔ t.edges s.g s.items e :=
    fun t ht e => (hallE t e).trans (hE₁ t (by simp [hts, ht]) e)
  have hbtouch : ∀ e, e < s.g.ne → b.edges s.g s.items e →
      s.g.Touches (b.edges s.g s.items) curV := fun e he hbe => hbv ▸ hE.touch_bot b hbmem ⟨e, he, hbe⟩
  have hmerge : MergeTopOk D s₃ := by
    intro cur nxt rest' h
    rw [hts₃] at h
    simp only [List.cons.injEq] at h
    obtain ⟨rfl, rfl, rfl⟩ := h
    refine ⟨⟨?_, ?_⟩, ?_⟩
    · rintro - ⟨e₂, he₂, hbe⟩
      rw [hg₃'] at he₂ hbe ⊢
      refine ⟨curV, ⟨o.e, he, (hcE₃ o.e).2 rfl, (Graph.inc_of_pairEq hends).1⟩, ?_⟩
      obtain ⟨e, he', hbe', hinc⟩ := hbtouch e₂ he₂ ((hbE₃ e₂ he₂).1 hbe)
      exact ⟨e, he', (hbE₃ e he').2 hbe', hinc⟩
    · exact Or.inl (Or.inl (show c.vStart = b'.vStart by rw [hc, hbv', hbv]))
    · intro t ht e he' hte hcb
      rw [hg₃'] at he' hte hcb
      rcases hcb with hce | hbe
      · obtain rfl := (hcE₃ e).1 hce
        obtain ⟨j, hj, hjb⟩ := (hrE₃ t ht _).1 hte
        exact hE.q_free t (by simp [hts, ht]) (Items.Below.eq_of_no_parent hE.q_root hjb ▸ hj)
      · exact (List.pairwise_cons.1 hpw).1 t ht e he' ((hbE₃ e he').1 hbe) ((hrE₃ t ht e).1 hte)
  have hfin : FinishTopOk D (after mergeTstackTops s₃) := by
    have hrun : after mergeTstackTops s₃ = { s₃ with tstack := TEntry.mergeInto c b' :: rest } := by
      show (mergeTstackTops.run s₃).2 = _; rw [mergeTstackTops_run_eq s₃ c b' rest hts₃]
    rw [hrun]
    set m := TEntry.mergeInto c b' with hm
    have hcur : curE { s₃ with tstack := m :: rest } = m := rfl
    have hmtop : m.topDepth = lv := by
      show min b'.topDepth c.topDepth = lv; rw [hbt', hbt]; exact Nat.min_self lv
    have hmv : m.vStart = curV := by show b'.vStart = curV; rw [hbv', hbv]
    have hmE : ∀ e, m.edges s.g s₃.items e ↔ (e = o.e ∨ b'.edges s.g s₃.items e) := fun e => by
      rw [hm, TEntry.edges_mergeInto, hcE₃]
    refine ⟨by simp, ?_, ?_⟩
    · rw [hcur, hmtop]
      show getSide (c.spans.1 ++ b'.spans.1, b'.spans.2 ++ c.spans.2) (!s₃.stackDir[lv]!) = []
      rw [hsd₃']
      exact getSide_merge_nil (!dir) c.spans b'.spans (getSide_setSides_not dir [q]) hside'
    · intro k hk hkD
      rw [hcur, hmtop] at hk
      rw [hcur, hmv]
      show s₃.stackVerts[k]! = curV ∨ s₃.g.Interior (m.edges s₃.g s₃.items) s₃.stackVerts[k]! ∨
        ¬ s₃.g.Touches (m.edges s₃.g s₃.items) s₃.stackVerts[k]!
      rw [hg₃', hsv₃']
      have hk' : s.stackVerts[k]! ≠ s.stackVerts[lv]! := fun h => hE.path lv k hk (by omega) h.symm
      by_cases ht : s.g.Touches (m.edges s.g s₃.items) s.stackVerts[k]!
      · obtain ⟨e, he', hme, hinc⟩ := ht
        rcases (hmE e).1 hme with rfl | hbe
        · left
          rcases hends with h | h <;> rcases hinc with h' | h'
          · exact h'.symm.trans (congrArg Prod.fst h).symm
          · exact absurd (h'.symm.trans (congrArg Prod.snd h).symm) hk'
          · exact absurd (h'.symm.trans (congrArg Prod.snd h).symm) hk'
          · exact h'.symm.trans (congrArg Prod.fst h).symm
        · rcases hatt _ ⟨e, he', (hbE₃ e he').1 hbe, hinc⟩ with h | h | h
          · exact .inl h
          · exact absurd h hk'
          · exact .inr (.inl fun e' he'' hinc' => (hmE e').2 (.inr ((hbE₃ e' he'').2 (h e' he'' hinc'))))
      · exact .inr (.inr ht)
  exact ⟨hu, hmerge, hfin⟩

/-- The merge of the vertex entry into the type-2 first-edge entry. -/
theorem ear_tail_tree (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hE : s.EarFinish curV d o hasVert sub base) (hi : s.Inv' D) (hs : Shape s) (hv : hasVert = false)
    (hsingle : feSingle d o s = false) (hD : D = d + 1) (he : o.e < s.g.ne)
    (hends : Items.PairEq (o.dest, s.stackVerts[d]!) s.g.edges[o.e]!) :
    MergeTopOk D (after (pushVertTstack curV d) (after (finishP curV lv o.cls.isType1) (feS₂ d o s))) := by
  have hlow' : o.cls.lowval d < d := by rw [ho]; exact hlow
  have hcond := ear_condP_tree ho hlow ht hg hE hi hs hv
  have hfin : after (finishP curV lv o.cls.isType1) (feS₂ d o s) = feS₂ d o s := by
    have h' : (o.cls.isType1 && decide ((feS₂ d o s).tstack.length ≥ 2) &&
        ((feS₂ d o s).tstack.tail.head!.vStart == curV) &&
        ((feS₂ d o s).tstack.tail.head!.topDepth == lv)) = false := hcond
    show ((finishP curV lv o.cls.isType1).run (feS₂ d o s)).2 = _
    simp only [Spqr.finishP, WalkM.run_bind, run_condP, h', Bool.false_eq_true, ↓reduceIte]
    rfl
  rw [hfin]
  obtain ⟨c, mid, py, vy, hC⟩ := hE.close ht hlow'
  subst hv hD
  set s₂ := feS₂ d o s with hs₂
  set ve : TEntry := ⟨curV, d, s₂.nxtEdgeIdx, setSides s₂.stackDir[d]! [vertItem curV] []⟩ with hve
  have hrun : after (pushVertTstack curV d) s₂ = { s₂ with tstack := ve :: s₂.tstack } := by
    show ((pushTstack curV d (vertItem curV)).run s₂).2 = _
    rw [run_pushTstack]
  rw [hrun]
  set s₃ : WalkState := { s₂ with tstack := ve :: s₂.tstack } with hs₃
  have hg₃ : s₃.g = s.g := hC.g
  have hsv₃ : s₃.stackVerts = s.stackVerts := hC.sv
  have hit₃ : s₃.items = s₂.items := rfl
  have hveE : ∀ e, ve.edges s.g s₂.items e ↔ Items.EdgeBelow s.g s₂.items (vertItem curV) e := fun e =>
    TEntry.edges_single s₂.stackDir[d]! (vertItem curV) (getSide_setSides_not _ _) (getSide_setSides_self _ _ _) e
  have hinc : s.g.Inc o.e curV := by rw [← hE.sv_d]; exact (Graph.inc_of_pairEq hends).2
  intro cur nxt rest h
  have h' : ve :: (c :: mid ++ [py, vy] ++ base) = cur :: nxt :: rest := by rw [← hC.tstack]; exact h
  simp only [List.cons.injEq] at h'
  obtain ⟨rfl, rfl, rfl⟩ := h'
  have hpw := hC.disj
  rw [hC.tstack] at hpw
  have hcrest := (List.pairwise_cons.1 hpw).1
  refine ⟨⟨?_, ?_⟩, ?_⟩
  · rintro ⟨e₁, he₁, hve₁⟩ -
    rw [hg₃] at he₁
    rw [hg₃, hit₃] at hve₁
    rw [hg₃, hit₃]
    obtain ⟨e, he', hb, hinc'⟩ := hC.vert_touch rfl ⟨e₁, he₁, (hveE e₁).1 hve₁⟩
    exact ⟨curV, ⟨e, he', (hveE e).2 hb, hinc'⟩, ⟨o.e, he, hC.c_edge, hinc⟩⟩
  · refine .inl (.inr ⟨d, ?_, by omega, ?_⟩)
    · show min c.topDepth d ≤ d
      exact Nat.min_le_right _ _
    · show curV = s₃.stackVerts[d]!
      rw [hsv₃, hE.sv_d]
  · intro t ht e he' hte hor
    rw [hg₃] at he'
    rw [hg₃, hit₃] at hte hor
    rcases hor with h | h
    · exact hC.vert_disj rfl t (by rw [hC.tstack]; exact List.mem_cons_of_mem _ ht) e he' hte ((hveE e).1 h)
    · exact hcrest t ht e he' h hte

/-- `FinishOk` from the guards, the invariant, the bookkeeping facts of the finished edge
(`he`, `hq`, `hends`, `hvert`) and the ear facts. The vertex item's connectivity/2-attachment is
transported to the tail through the `Step`s of the preceding blocks. -/
theorem finishOk_of_guards (ho : o.cls = .ret lv kind) (hlow : lv < d)
    (hg : FinishGuards d o origTstack hasVert s) (hE : s.EarFinish curV d o hasVert sub base)
    (hlen : base.length = origTstack) (hi : s.Inv' D) (hs : Shape s)
    (hD : D = if o.cls.isTree then d + 1 else d) (hv : curV < s.g.nv) (he : o.e < s.g.ne)
    (hq : Items.ch s.items (edgeItem s.g o.e) = [])
    (hends : Items.PairEq (if o.cls.isTree then (o.dest, s.stackVerts[d]!) else (curV, s.stackVerts[lv]!))
      s.g.edges[o.e]!)
    (hvert : hasVert = false → s.g.ConnEdges (Items.EdgeBelow s.g s.items (vertItem curV)) ∧
      s.g.TwoAttached (Items.EdgeBelow s.g s.items (vertItem curV)) curV curV) :
    FinishOk D curV d lv o origTstack hasVert s := by
  have hdD : d ≤ D := by split at hD <;> omega
  have hq₀ : Items.ch (feS₀ d o s).items (edgeItem (feS₀ d o s).g o.e) = [] := by
    show Items.ch (s.items.modify _ _) (edgeItem s.g o.e) = []
    rw [Items.ch_modify_ch_eq (edgeItem s.g o.e)
      (fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) (fun _ => rfl)]
    exact hq
  have hears : o.cls.isTree = true → CloseEarsOk D o.dest d o.e s.stackDir[d]! (feS₀ d o s) := fun ht =>
    ⟨he, hq₀, by show Items.PairEq (o.dest, s.stackVerts[d]!) s.g.edges[o.e]!; simpa [ht] using hends,
     hdD, ear_loop1 ho hlow ht hg hE hi hs (by rw [hD, if_pos ht]) hv he hq (by simpa [ht] using hends)⟩
  have st₀ : Step D curV s (feS₀ d o s) :=
    Step.modifyVs hi hs (edgeItem s.g o.e) _ (by show 1 + s.g.nv + o.e < _; omega)
  have hv₀ : curV < (feS₀ d o s).g.nv := by rw [st₀.g]; exact hv
  have st₁ : o.cls.isTree = true → Step D curV (feS₀ d o s) (feS₁ d o s) := fun ht =>
    Step.closeEars st₀.inv st₀.shape hv₀ (hears ht)
  have st₂ : o.cls.isTree = true → Step D curV (feS₁ d o s) (feS₂ d o s) := fun ht =>
    Step.mergeLate (st₁ ht).inv (st₁ ht).shape (by rw [(st₁ ht).g]; exact hv₀)
      (ear_mergeLate ho hlow ht hg hE hi hs (by rw [hD, if_pos ht]))
  have hv₂ : o.cls.isTree = true → curV < (feS₂ d o s).g.nv := fun ht => by
    rw [(st₂ ht).g, (st₁ ht).g]; exact hv₀
  have st₃ : o.cls.isTree = true → hasVert = true → Step D curV (feS₂ d o s) (feS₃ curV d o origTstack s) :=
    fun ht hv' => Step.closeVert' (st₂ ht).inv (st₂ ht).shape (hv₂ ht)
      (ear_closeVert ho hlow ht hg hE hi hs hv' (by rw [hD, if_pos ht]) hlen (st₂ ht).shape)
  refine
    { e_lt := he
      ears := hears
      late := fun ht => ear_mergeLate ho hlow ht hg hE hi hs (by rw [hD, if_pos ht])
      vert := fun ht hv' => ear_closeVert ho hlow ht hg hE hi hs hv' (by rw [hD, if_pos ht]) hlen (st₂ ht).shape
      rest_vert := fun ht hv' => ⟨ear_finishP_vert ho hlow ht hg hE hi hs hv' (by rw [hD, if_pos ht]) hlen
          (st₂ ht).inv (st₂ ht).shape (st₃ ht hv').shape he (by simpa [ht] using hends),
        fun h => by simp [hv'] at h, fun h => by simp [hv'] at h, fun h => by simp [hv'] at h⟩
      rest_tree := fun ht hv' => ?_
      q := fun _ => hq
      ends := fun hb => by simpa [hb] using hends
      lv_le := fun _ => by omega
      rest_back := fun hb => ?_ }
  · have st₁ := st₁ ht
    have st₂ := st₂ ht
    have hv₂ : curV < (feS₂ d o s).g.nv := by rw [st₂.g, st₁.g]; exact hv₀
    have hp := ear_finishP_tree (curV := curV) ho hlow ht hg hE hi hs hv'
    have st₃ := Step.finishP (v := curV) st₂.inv st₂.shape hv₂ hp
    have st := st₀.trans (st₁.trans (st₂.trans st₃))
    obtain ⟨hc, ha⟩ := hvert hv'
    exact ⟨hp, finishTailOk_of_vert (fun _ => (st.vertTransport hc ha).1) (fun _ => (st.vertTransport hc ha).2)
      fun _ hsg => ear_tail_tree ho hlow ht hg hE hi hs hv' hsg (by rw [hD, if_pos ht]) he (by simpa [ht] using hends)⟩
  · have st₁ : Step D curV _ (after (pushEdgeTstack curV lv o.e) (feS₀ d o s)) :=
      Step.pushEdge st₀.inv st₀.shape curV lv o.e he hq₀
        (by show Items.PairEq (curV, s.stackVerts[lv]!) s.g.edges[o.e]!; simpa [hb] using hends) (by omega)
    have st₂ : Step D curV _ (feBack curV lv d o s) := Step.frame st₁.inv st₁.shape rfl rfl rfl rfl
    have hv₂ : curV < (feBack curV lv d o s).g.nv := by rw [st₂.g, st₁.g]; exact hv₀
    have hp := ear_finishP_back (curV := curV) ho hlow hb hg hE hi hs (by simpa [hb] using hD) he hq
      (by simpa [hb] using hends) st₂.shape
    have st₃ := Step.finishP (v := curV) st₂.inv st₂.shape hv₂ hp
    have st := st₀.trans (st₁.trans (st₂.trans st₃))
    exact ⟨hp, finishTailOk_of_vert (fun h => (st.vertTransport (hvert h).1 (hvert h).2).1)
      (fun h => (st.vertTransport (hvert h).1 (hvert h).2).2) fun _ h => by cases h⟩

end Ear

/-! ## `walkTree_inv` by mutual induction

Hypotheses are supplied in the style of `GuardsTree`: the ear facts through `FinishGuards` and the
bookkeeping facts (vertex/edge bounds, fresh `Q`/`V` items, edge endpoints) through `BookTree`. -/

section Walk

variable {α : Type}

theorem Shape.frame' {s' : WalkState} (h : Shape s) (hg : s'.g = s.g := by rfl)
    (hi : s'.items = s.items := by rfl) (hts : s'.tstack = s.tstack := by rfl) : Shape s' :=
  h.frame hg hi hts

theorem Inv'.frame' {D : Nat} {s' : WalkState} (h : s.Inv' D) (hg : s'.g = s.g := by rfl)
    (hsv : s'.stackVerts = s.stackVerts := by rfl) (hi : s'.items = s.items := by rfl)
    (hts : s'.tstack = s.tstack := by rfl) : s'.Inv' D :=
  h.frame hg hsv hi hts

theorem wp_of_forall {m : WalkM α} {Q : α → WalkState → Prop} (h : ∀ a s', Q a s') : wp m Q s := h _ _

theorem wp_imp {m : WalkM α} {P Q : α → WalkState → Prop} (h : wp m (fun a s' => P a s' → Q a s') s)
    (hp : wp m P s) : wp m Q s := h hp

theorem wp_and {m : WalkM α} {P Q : α → WalkState → Prop} (h₁ : wp m P s) (h₂ : wp m Q s) :
    wp m (fun a s' => P a s' ∧ Q a s') s := ⟨h₁, h₂⟩

theorem Inv'.setSv {d : Nat} (x : Nat) (h : s.Inv' d) :
    ({ s with stackVerts := s.stackVerts.set! (d + 1) x } : WalkState).Inv' (d + 1) := by
  have hsv : ∀ k, k ≤ d → (s.stackVerts.set! (d + 1) x)[k]! = s.stackVerts[k]! := by
    intro k hk
    have hk' : k ≠ d + 1 := by omega
    simp [Array.set!, getElem!_def, Ne.symm hk']
  have hT : ∀ (t : TEntry) v, t.Term d s v →
      t.Term (d + 1) { s with stackVerts := s.stackVerts.set! (d + 1) x } v := by
    rintro t v (hv | ⟨k, h1, h2, h3⟩)
    · exact .inl hv
    · exact .inr ⟨k, h1, by omega, by rw [h3]; exact (hsv k h2).symm⟩
  refine Inv'.ofStack ?_ fun i hi hsz => ItemInv.congr (s := s) rfl rfl (fun _ _ => Iff.rfl) (h.nodes i hi hsz)
  suffices ∀ A l, s.Stack d A l → Stack (d + 1) { s with stackVerts := s.stackVerts.set! (d + 1) x } A l from
    this [] _ h.stack
  intro A l hl
  induction l generalizing A with
  | nil => exact Stack_nil
  | cons t rest ih =>
    obtain ⟨ht, hrest⟩ := Stack_cons.1 hl
    refine Stack_cons.2 ⟨⟨ht.conn, ht.attached.mono fun v hv => ?_⟩, ih _ hrest⟩
    rcases hv with hv | ⟨t', ht', hT'⟩
    · exact .inl (hT t v hv)
    · exact .inr ⟨t', ht', hT t' v hT'⟩

theorem ret_of_lowval_lt {o : DfsOut} {d : Nat} (h : o.cls.lowval d < d) :
    ∃ lv kind, o.cls = .ret lv kind ∧ lv < d := by
  cases hc : o.cls <;> simp only [OutClass.lowval, hc] at h
  all_goals first | omega | exact ⟨_, _, rfl, h⟩

/-- Bookkeeping facts needed by `finishEdge_inv` at the state where `finishEdge` runs. -/
structure FinishBook (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool) (s : WalkState) : Prop where
  v_lt : curV < s.g.nv
  e_lt : o.e < s.g.ne
  q : Items.ch s.items (edgeItem s.g o.e) = []
  ends : ∀ lv kind, o.cls = .ret lv kind →
    Items.PairEq (if o.cls.isTree then (o.dest, s.stackVerts[d]!) else (curV, s.stackVerts[lv]!))
      s.g.edges[o.e]!
  vert : hasVert = false → s.g.ConnEdges (Items.EdgeBelow s.g s.items (vertItem curV)) ∧
    s.g.TwoAttached (Items.EdgeBelow s.g s.items (vertItem curV)) curV curV
  tree : o.cls.isTree = true ↔ ∃ e cls c, o = .tree e cls c
  /-- The ear content of the stack (`EarInv.lean`). -/
  ear : s.EarAt curV d o origTstack hasVert

/-- Before the vertex entry of `v` is pushed, `v` is in range and its item is a connected piece
attached only at `v` (empty before the first out-edge; the blocks closed by the boundary edges of
`v` afterwards — `hasVert = false → ch (vertItem v) = []` is false after a bridge/component edge,
e.g. edges `0-1, 1-2, 1-0`: at vertex 1 the bridge `1-2` is finished first, leaving
`ch (vertItem 1) = [Q(1-2)]` with `hasVert = false`). -/
def VertBook (v : Nat) (hasVert : Bool) (s : WalkState) : Prop :=
  hasVert = false → v < s.g.nv ∧ s.g.ConnEdges (Items.EdgeBelow s.g s.items (vertItem v)) ∧
    s.g.TwoAttached (Items.EdgeBelow s.g s.items (vertItem v)) v v

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
        wp (walkTree child (d + 1)) (fun _ s₃ => FinishBook v d o s₁.tstack.length hasVert' s₃) s₂) s₁
    | .back .. => FinishBook v d o s₁.tstack.length hasVert' s₁) s
end

mutual
/-- The out-edges of `v` are edges of `g` with the DFS endpoints: a tree edge joins the child's
vertex and `v`, a back edge joins `v` and its destination. -/
def _root_.Spqr.DfsOut.Ends (g : Graph) (v : Nat) : DfsOut → Prop
  | .back e dest _ => Items.PairEq (v, dest) g.edges[e]!
  | .tree e _ child => Items.PairEq (child.v, v) g.edges[e]! ∧ DfsTree.Ends g child
def _root_.Spqr.DfsTree.Ends (g : Graph) : DfsTree → Prop
  | .node v outs => ∀ o ∈ outs, DfsOut.Ends g v o
end

theorem _root_.Spqr.dfsVisit_root_v (adj : Array (List (Nat × Nat))) :
    ∀ (fuel v d : Nat) (prvE : Option Nat) (depth : Array (Option Nat)),
      (dfsVisit adj fuel v d prvE depth).1.v = v
  | 0, _, _, _, _ => rfl
  | _ + 1, _, _, _, _ => by rw [dfsVisit_succ]; exact rfl

theorem _root_.Spqr.dfsVisit_ends {g : Graph} {adj : Array (List (Nat × Nat))}
    (hadj : ∀ x y e : Nat, (y, e) ∈ adj[x]! → Items.PairEq (x, y) g.edges[e]!) :
    ∀ (fuel v d : Nat) (prvE : Option Nat) (depth : Array (Option Nat)),
      (dfsVisit adj fuel v d prvE depth).1.Ends g
  | 0, v, _, _, _ => by simp [dfsVisit, DfsTree.Ends]
  | fuel + 1, v, d, prvE, depth => by
    rw [dfsVisit_succ]
    have key : ∀ (l : List (Nat × Nat)) (acc : List DfsOut × Lowvals × Array (Option Nat)),
        (∀ p ∈ l, p ∈ adj[v]!) → (∀ o ∈ acc.1, DfsOut.Ends g v o) →
        ∀ o ∈ (l.foldl (dfsStep adj fuel d prvE) acc).1, DfsOut.Ends g v o := by
      intro l
      induction l with
      | nil => intro acc _ h; simpa using h
      | cons p l ih =>
        rintro ⟨outs, lv, cur⟩ hl h
        refine ih _ (fun q hq => hl q (List.mem_cons_of_mem _ hq)) ?_
        obtain ⟨nxt, e⟩ := p
        have hpe : Items.PairEq (v, nxt) g.edges[e]! := hadj _ _ _ (hl _ (List.mem_cons_self ..))
        simp only [dfsStep]
        split
        · exact h
        · split
          · obtain ⟨⟨c, n, dp⟩, hr⟩ : ∃ r, dfsVisit adj fuel nxt (d + 1) (some e) cur = r := ⟨_, rfl⟩
            have hv := dfsVisit_root_v adj fuel nxt (d + 1) (some e) cur
            have he := dfsVisit_ends hadj fuel nxt (d + 1) (some e) cur
            rw [hr] at hv he ⊢
            intro o ho
            rcases List.mem_cons.1 ho with rfl | ho
            · simp only [DfsOut.Ends] at hv ⊢
              refine ⟨?_, he⟩
              rw [hv]
              rcases hpe with hpe | hpe
              · exact .inr (by rw [Prod.ext_iff] at hpe ⊢; exact ⟨hpe.2, hpe.1⟩)
              · exact .inl (by rw [Prod.ext_iff] at hpe ⊢; exact ⟨hpe.2, hpe.1⟩)
            · exact h o ho
          · intro o ho
            rcases List.mem_cons.1 ho with rfl | ho
            · simp only [DfsOut.Ends]; exact hpe
            · exact h o ho
    simp only [DfsTree.Ends]
    intro o ho
    exact key adj[v]! _ (fun _ h => h) (by simp) o
      (List.mem_reverse.1 ((List.mergeSort_perm _ _).subset ho))

set_option linter.unusedVariables false in
/-- Every tree of `dfsForest` has the DFS endpoints (`DfsTree.Ends`): its tree edges join child
and parent, its back edges join the vertex and the destination. -/
theorem _root_.Spqr.dfsForest_ends (g : Graph) (hg : g.WF) {vo eo : List Nat}
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) : ∀ t ∈ g.dfsForest vo eo, t.Ends g := by
  have hadj : ∀ x y e : Nat, (y, e) ∈ (g.adjacency eo)[x]! → Items.PairEq (x, y) g.edges[e]! :=
    fun x y e h => by
      rcases adjacency_ends hg heo x y e h with h | h
      · exact .inl h.symm
      · exact .inr (by rw [h])
  have key : ∀ (l : List Nat) (acc : List DfsTree × Array (Option Nat)), (∀ t ∈ acc.1, t.Ends g) →
      ∀ t ∈ (l.foldl (forestStep (g.adjacency eo) g.nv) acc).1, t.Ends g := by
    intro l
    induction l with
    | nil => intro acc h; simpa using h
    | cons rt l ih =>
      rintro ⟨roots, depth⟩ h
      refine ih _ ?_
      simp only [forestStep]
      split
      · exact h
      · obtain ⟨⟨t', n, depth'⟩, hr⟩ : ∃ r, dfsVisit (g.adjacency eo) g.nv rt 0 none depth = r :=
          ⟨_, rfl⟩
        have he := dfsVisit_ends hadj g.nv rt 0 none depth
        rw [hr] at he ⊢
        intro t ht
        rcases List.mem_cons.1 ht with rfl | ht
        · exact he
        · exact h t ht
  rw [dfsForest_eq]
  intro t ht
  exact key _ _ (by simp) t (List.mem_reverse.1 ht)

/-- The ear fields of `FinishBook`: the content of a `finishEdge` that only the ear preservation
induction supplies (the vertex item's shape before the vertex entry, and the `EarAt` contract). -/
structure FinishEar (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool) (s : WalkState) : Prop where
  vert : hasVert = false → s.g.ConnEdges (Items.EdgeBelow s.g s.items (vertItem curV)) ∧
    s.g.TwoAttached (Items.EdgeBelow s.g s.items (vertItem curV)) curV curV
  ear : s.EarAt curV d o origTstack hasVert

mutual
/-- `BookTree` restricted to the ear fields (`FinishEar` instead of `FinishBook`). -/
def EarTree (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  match t with
  | .node v outs => EarOuts v d outs false { s with stackVerts := s.stackVerts.set! d v }

def EarOuts (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  match outs with
  | [] => VertBook v hasVert s
  | o :: rest => EarOut v d o hasVert s ∧
      wp (walkOut v d o hasVert) (fun hasVert' s' => EarOuts v d rest hasVert' s') s

def EarOut (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  VertBook v hasVert s ∧
  wp (walkOutPre v d o hasVert) (fun hasVert' s₁ =>
    match o with
    | .tree _ _ child =>
      wp (modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }) (fun _ s₂ =>
        EarTree child (d + 1) s₂ ∧
        wp (walkTree child (d + 1)) (fun _ s₃ => FinishEar v d o s₁.tstack.length hasVert' s₃) s₂) s₁
    | .back .. => FinishEar v d o s₁.tstack.length hasVert' s₁) s
end

theorem edgeItem_lt {g : Graph} {e : Nat} (he : e < g.ne) : edgeItem g e < 1 + g.nv + g.ne := by
  show (1 + g.nv + e : Nat) < 1 + g.nv + g.ne; omega

theorem edgeItem_inj {g : Graph} {e e' : Nat} (h : edgeItem g e = edgeItem g e') : e = e' := by
  have h' : (1 + g.nv + e : Nat) = 1 + g.nv + e' := h; omega

theorem vertItem_ne_zero (v : Nat) : vertItem v ≠ 0 := by simp [vertItem]
theorem edgeItem_ne_zero (g : Graph) (e : Nat) : edgeItem g e ≠ 0 := by simp [edgeItem]

/-! ### Assembly of `walkTree_ear` from `EarCtx` (PROOF.md §4.3)

The mutual induction over `walkTree`/`walkOuts`/`walkOut` carries the between-edges context
`EarCtx` (EarCtx.lean); the site contracts come from `earAt_back_of_ctx`/`earAt_tree_of_ctx`
(EarCtxAt.lean). The step facts that are not yet proved are the named admissions `ctx_init_root`,
`ctx_init_child`, `ctx_step_back`, `ctx_step_tree`, `walkTree_qroot_kept`, `walkTree_vroot_kept`,
`walkTree_below_kept` and the DFS fact `ends_of_wf_boundary`; each is checked literally in
`checks/EarCheck.lean` (`ctxCheck` at every `iOuts` entry / after every out, `kept_*`). -/

theorem subEdges_e (o : DfsOut) : subEdges o o.e := by cases o <;> simp [subEdges]

theorem mem_edges_of_subEdges {o : DfsOut} {e : Nat} (h : subEdges o e) : e ∈ o.edges := by
  cases o <;> simp_all [subEdges, DfsOut.edges, DfsOut.e]

theorem mem_vertsList_of_verts {outs : List DfsOut} {o : DfsOut} {x : Nat} (ho : o ∈ outs)
    (hx : x ∈ o.verts) : x ∈ DfsOut.vertsList outs := by
  rw [DfsOut.vertsList_eq]; exact List.mem_flatMap.2 ⟨o, ho, hx⟩

theorem eq_of_inc_pairEq {g : Graph} {e a b x : Nat} (h : Items.PairEq (a, b) g.edges[e]!)
    (hx : g.Inc e x) : x = a ∨ x = b := by
  unfold Items.PairEq at h; unfold Graph.Inc at hx
  rcases h with h | h
  · rw [← h] at hx
    rcases hx with hx | hx
    · exact .inl hx.symm
    · exact .inr hx.symm
  · have ha : a = (g.edges[e]!).2 := congrArg Prod.fst h
    have hb : b = (g.edges[e]!).1 := congrArg Prod.snd h
    rcases hx with hx | hx
    · exact .inr (by rw [hb]; exact hx.symm)
    · exact .inl (by rw [ha]; exact hx.symm)

theorem inc_e_of_ends {g : Graph} {v : Nat} {o : DfsOut} (h : DfsOut.Ends g v o) : g.Inc o.e v := by
  cases o with
  | back e dest cls => rw [DfsOut.Ends] at h; exact (Graph.inc_of_pairEq h).1
  | tree e cls child => rw [DfsOut.Ends] at h; exact (Graph.inc_of_pairEq h.1).2

-- Endpoints of the edges below a well-formed out: the vertex, its subtree, or an ancestor.
mutual
theorem endsOut_wf : ∀ (g : Graph) (anc : List Nat) (v : Nat) (o : DfsOut), o.WF anc v →
    DfsOut.Ends g v o → ∀ e, subEdges o e → ∀ x, g.Inc e x → x ∈ anc ∨ x = v ∨ x ∈ o.verts
  | g, anc, v, .back e' dest cls, hwf, hE, e, he, x, hx => by
    rw [DfsOut.WF] at hwf
    obtain ⟨i, hi, -⟩ := hwf
    rw [DfsOut.Ends] at hE
    have he' : e = e' := by simpa [subEdges, DfsOut.e] using he
    subst he'
    rcases eq_of_inc_pairEq hE hx with rfl | rfl
    · exact .inr (.inl rfl)
    · rcases List.mem_append.1 (List.mem_of_getElem? hi) with h | h
      · exact .inl h
      · exact .inr (.inl (List.mem_singleton.1 h))
  | g, anc, v, .tree e' cls (.node y outs), hwf, hE, e, he, x, hx => by
    rw [DfsOut.WF] at hwf
    obtain ⟨hwf', -⟩ := hwf
    rw [DfsTree.WF] at hwf'
    rw [DfsOut.Ends] at hE
    obtain ⟨hpe, hE'⟩ := hE
    rw [DfsTree.Ends] at hE'
    rcases he with rfl | he
    · rcases eq_of_inc_pairEq hpe hx with rfl | rfl
      · exact .inr (.inr (List.mem_cons_self ..))
      · exact .inr (.inl rfl)
    · rcases endsOuts_wf g (anc ++ [v]) y outs hwf'.2 hE' e he x hx with h | rfl | h
      · rcases List.mem_append.1 h with h | h
        · exact .inl h
        · exact .inr (.inl (List.mem_singleton.1 h))
      · exact .inr (.inr (List.mem_cons_self ..))
      · exact .inr (.inr (List.mem_cons_of_mem _ h))
theorem endsOuts_wf : ∀ (g : Graph) (anc : List Nat) (v : Nat) (outs : List DfsOut),
    (∀ o ∈ outs, o.WF anc v) → (∀ o ∈ outs, DfsOut.Ends g v o) →
    ∀ e, e ∈ DfsOut.edgesList outs → ∀ x, g.Inc e x → x ∈ anc ∨ x = v ∨ x ∈ DfsOut.vertsList outs
  | g, anc, v, [], _, _, e, he, x, hx => by simp [DfsOut.edgesList] at he
  | g, anc, v, o :: rest, hwf, hE, e, he, x, hx => by
    obtain ⟨o', ho', hs⟩ := mem_subEdges_edgesList.1 he
    rcases List.mem_cons.1 ho' with h' | ho'
    · rw [h'] at hs
      rcases endsOut_wf g anc v o (hwf _ (List.mem_cons_self ..)) (hE _ (List.mem_cons_self ..))
        e hs x hx with h | h | h
      · exact .inl h
      · exact .inr (.inl h)
      · exact .inr (.inr (mem_vertsList_of_verts (List.mem_cons_self ..) h))
    · rcases endsOuts_wf g anc v rest (fun o h => hwf o (List.mem_cons_of_mem _ h))
        (fun o h => hE o (List.mem_cons_of_mem _ h)) e (mem_subEdges_edgesList.2 ⟨o', ho', hs⟩) x hx
        with h | h | h
      · exact .inl h
      · exact .inr (.inl h)
      · rw [DfsOut.vertsList_eq] at h
        obtain ⟨o'', ho'', hx'⟩ := List.mem_flatMap.1 h
        exact .inr (.inr (mem_vertsList_of_verts (List.mem_cons_of_mem _ ho'') hx'))
end

theorem back_wf_facts {anc : List Nat} {v e dest : Nat} {cls : OutClass}
    (h : DfsOut.WF anc v (.back e dest cls)) :
    cls.isTree = false ∧ ∀ lv k, cls = .ret lv k → lv < anc.length := by
  rw [DfsOut.WF] at h
  obtain ⟨i, hi, rfl⟩ := h
  have hi' : i < anc.length + 1 := by
    have := (List.getElem?_eq_some_iff.1 hi).1; simpa using this
  refine ⟨?_, fun lv k h => (classify_eq_ret_iff.1 h).2.1⟩
  unfold classify
  by_cases h1 : anc.length ≤ i
  · have : i = anc.length := by omega
    subst this; simp [OutClass.isTree]
  · simp [h1, OutClass.isTree]

theorem pairEq_swap {p : Nat × Nat} {a b : Nat} (h : Items.PairEq (a, b) p) :
    Items.PairEq (b, a) p := by
  rcases h with h | h
  · exact Or.inr (by rw [← h])
  · exact Or.inl (Prod.ext (Prod.mk.inj h).2 (Prod.mk.inj h).1)

theorem reach_of_pairEq {g : Graph} {E : Nat → Prop} {e v w : Nat} (he : e < g.ne) (hE : E e)
    (h : Items.PairEq (v, w) g.edges[e]!) :
    Relation.ReflTransGen (g.AdjIn E) v (g.edges[e]!).1 := by
  rcases h with h | h
  · rw [← h]
  · exact .single ⟨e, he, hE, Or.inr (by rw [(Prod.mk.inj h).1])⟩

-- The edges below an out are connected to its vertex through the out's subtree.
mutual
theorem out_reach : ∀ (g : Graph) (E : Nat → Prop) (v : Nat) (o : DfsOut), DfsOut.Ends g v o →
    (∀ e', subEdges o e' → e' < g.ne ∧ E e') → ∀ e', subEdges o e' →
    Relation.ReflTransGen (g.AdjIn E) v (g.edges[e']!).1
  | g, E, v, .back e dest cls, hE, hb, e', he' => by
    rw [DfsOut.Ends] at hE
    have h : e' = e := by simpa [subEdges, DfsOut.e] using he'
    subst h
    exact reach_of_pairEq (hb _ he').1 (hb _ he').2 hE
  | g, E, v, .tree e cls (.node y outs), hE, hb, e', he' => by
    rw [DfsOut.Ends] at hE
    obtain ⟨hpe, hE'⟩ := hE
    rw [DfsTree.Ends] at hE'
    rcases he' with he' | he'
    · subst he'
      exact reach_of_pairEq (hb _ (Or.inl rfl)).1 (hb _ (Or.inl rfl)).2 (pairEq_swap hpe)
    · exact (Relation.ReflTransGen.single ⟨e, (hb e (Or.inl rfl)).1, (hb e (Or.inl rfl)).2,
        pairEq_swap hpe⟩).trans (outs_reach g E y outs hE' (fun e'' h => hb _ (Or.inr h)) _ he')
theorem outs_reach : ∀ (g : Graph) (E : Nat → Prop) (v : Nat) (outs : List DfsOut),
    (∀ o ∈ outs, DfsOut.Ends g v o) → (∀ e' ∈ DfsOut.edgesList outs, e' < g.ne ∧ E e') →
    ∀ e' ∈ DfsOut.edgesList outs, Relation.ReflTransGen (g.AdjIn E) v (g.edges[e']!).1
  | g, E, v, [], _, _, e', he' => by simp [DfsOut.edgesList] at he'
  | g, E, v, o :: rest, hE, hb, e', he' => by
    obtain ⟨o', ho', hs⟩ := mem_subEdges_edgesList.1 he'
    rcases List.mem_cons.1 ho' with h' | ho'
    · rw [h'] at hs
      exact out_reach g E v o (hE _ (List.mem_cons_self ..))
        (fun e'' h => hb _ (mem_subEdges_edgesList.2 ⟨_, List.mem_cons_self .., h⟩)) e' hs
    · exact outs_reach g E v rest (fun o ho => hE _ (List.mem_cons_of_mem _ ho))
        (fun e'' h => hb _ (mem_subEdges_edgesList.2
          (let ⟨o, ho, hs⟩ := mem_subEdges_edgesList.1 h; ⟨o, List.mem_cons_of_mem _ ho, hs⟩)))
        e' (mem_subEdges_edgesList.2 ⟨o', ho', hs⟩)
end

theorem tree_wf_facts {anc : List Nat} {v e : Nat} {cls : OutClass} {child : DfsTree}
    (h : DfsOut.WF anc v (.tree e cls child)) :
    cls.isTree = true ∧ ∀ lv k, cls = .ret lv k → lv < anc.length := by
  rw [DfsOut.WF] at h
  obtain ⟨-, rfl⟩ := h
  refine ⟨?_, fun lv k h => (classify_eq_ret_iff.1 h).2.1⟩
  unfold classify
  split
  · split <;> simp [OutClass.isTree]
  · simp only [↓reduceIte]
    split <;> simp [OutClass.isTree]

theorem length_le_nv {l : List Nat} {n : Nat} (hnd : l.Nodup) (h : ∀ a ∈ l, a < n) : l.length ≤ n := by
  have hsub : l ⊆ List.range n := fun w hw => List.mem_range.2 (h w hw)
  simpa using (List.Nodup.subperm hnd hsub).length_le

theorem getElem!_map_fn {α : Type} [Inhabited α] (l : List α) (f : α → Nat → Prop) (k : Nat)
    (hk : k < l.length) : (l.map f)[k]! = f l[k]! := by
  rw [getElem!_pos (l.map f) k (by simpa using hk), getElem!_pos l k hk, List.getElem_map]

mutual
/-- Below the ancestor path `anc ++ pre`, if every return depth of an out of `w` is `≥ anc.length`,
every edge of the out has its endpoints in `pre`, `w` or the out's vertices (a back edge's
destination `(anc ++ pre ++ [w])[i]` has `i = lowval ≥ anc.length`). -/
theorem bdOut_wf : ∀ (g : Graph) (anc pre : List Nat) (w : Nat) (o : DfsOut), o.WF (anc ++ pre) w →
    DfsOut.Ends g w o → (∀ r ∈ o.retDepths (anc ++ pre).length, anc.length ≤ r) →
    ∀ e, subEdges o e → ∀ x, g.Inc e x → x ∈ pre ∨ x = w ∨ x ∈ o.verts
  | g, anc, pre, w, .back e' dest cls, hwf, hE, hret, e, he, x, hx => by
    rw [DfsOut.WF] at hwf
    obtain ⟨i, hi, hcls⟩ := hwf
    rw [DfsOut.Ends] at hE
    have he' : e = e' := by simpa [subEdges, DfsOut.e] using he
    subst he'
    rcases eq_of_inc_pairEq hE hx with rfl | rfl
    · exact .inr (.inl rfl)
    · have hi' : i < (anc ++ pre ++ [w]).length := (List.getElem?_eq_some_iff.1 hi).1
      have hile : i ≤ (anc ++ pre).length := by simp at hi' ⊢; omega
      have hr := hret (cls.lowval (anc ++ pre).length) (by simp [DfsOut.retDepths])
      rw [hcls, lowval_classify_back hile] at hr
      rw [List.append_assoc, List.getElem?_append_right hr] at hi
      rcases List.mem_append.1 (List.mem_of_getElem? hi) with h | h
      · exact .inl h
      · exact .inr (.inl (List.mem_singleton.1 h))
  | g, anc, pre, w, .tree e' cls (.node y outs), hwf, hE, hret, e, he, x, hx => by
    rw [DfsOut.WF] at hwf
    obtain ⟨hwf', -⟩ := hwf
    rw [DfsTree.WF] at hwf'
    rw [DfsOut.Ends] at hE
    obtain ⟨hpe, hE'⟩ := hE
    rw [DfsTree.Ends] at hE'
    rcases he with rfl | he
    · rcases eq_of_inc_pairEq hpe hx with rfl | rfl
      · exact .inr (.inr (List.mem_cons_self ..))
      · exact .inr (.inl rfl)
    · have hret' : ∀ r ∈ DfsOut.retDepthsList (anc ++ (pre ++ [w])).length outs, anc.length ≤ r := by
        intro r hr
        apply hret
        simp only [DfsOut.retDepths, DfsTree.retDepths]
        have hl : (anc ++ (pre ++ [w])).length = (anc ++ pre).length + 1 := by simp [Nat.add_assoc]
        rw [hl] at hr
        exact hr
      rcases bdOuts_wf g anc (pre ++ [w]) y outs (by rw [← List.append_assoc]; exact hwf'.2) hE'
        hret' e he x hx with h | rfl | h
      · rcases List.mem_append.1 h with h | h
        · exact .inl h
        · exact .inr (.inl (List.mem_singleton.1 h))
      · exact .inr (.inr (List.mem_cons_self ..))
      · exact .inr (.inr (List.mem_cons_of_mem _ h))

theorem bdOuts_wf : ∀ (g : Graph) (anc pre : List Nat) (w : Nat) (outs : List DfsOut),
    (∀ o ∈ outs, o.WF (anc ++ pre) w) → (∀ o ∈ outs, DfsOut.Ends g w o) →
    (∀ r ∈ DfsOut.retDepthsList (anc ++ pre).length outs, anc.length ≤ r) →
    ∀ e, e ∈ DfsOut.edgesList outs → ∀ x, g.Inc e x → x ∈ pre ∨ x = w ∨ x ∈ DfsOut.vertsList outs
  | g, anc, pre, w, [], _, _, _, e, he, x, hx => by simp [DfsOut.edgesList] at he
  | g, anc, pre, w, o :: rest, hwf, hE, hret, e, he, x, hx => by
    obtain ⟨o', ho', hs⟩ := mem_subEdges_edgesList.1 he
    rw [DfsOut.retDepthsList_eq] at hret
    rcases List.mem_cons.1 ho' with h' | ho'
    · rw [h'] at hs
      rcases bdOut_wf g anc pre w o (hwf _ (List.mem_cons_self ..)) (hE _ (List.mem_cons_self ..))
        (fun r hr => hret r (List.mem_flatMap.2 ⟨o, List.mem_cons_self .., hr⟩)) e hs x hx with h | h | h
      · exact .inl h
      · exact .inr (.inl h)
      · exact .inr (.inr (mem_vertsList_of_verts (List.mem_cons_self ..) h))
    · rcases bdOuts_wf g anc pre w rest (fun o h => hwf o (List.mem_cons_of_mem _ h))
        (fun o h => hE o (List.mem_cons_of_mem _ h))
        (fun r hr => by
          rw [DfsOut.retDepthsList_eq] at hr
          obtain ⟨o'', ho'', hr'⟩ := List.mem_flatMap.1 hr
          exact hret r (List.mem_flatMap.2 ⟨o'', List.mem_cons_of_mem _ ho'', hr'⟩))
        e (mem_subEdges_edgesList.2 ⟨o', ho', hs⟩) x hx with h | h | h
      · exact .inl h
      · exact .inr (.inl h)
      · rw [DfsOut.vertsList_eq] at h
        obtain ⟨o'', ho'', hx'⟩ := List.mem_flatMap.1 h
        exact .inr (.inr (mem_vertsList_of_verts (List.mem_cons_of_mem _ ho'') hx'))
end

/-- Below a boundary tree edge (`d ≤ lowval`, i.e. no return above `v`), every subtree edge has both
endpoints in `v :: child.verts` (dump-checked `kept_ends`): `bdOut_wf` with `pre = []`, the lowval
being the least return depth (`lowval_eq_lmin`). -/
theorem ends_of_wf_boundary {g : Graph} {anc : List Nat} {v e : Nat} {cls : OutClass} {child : DfsTree}
    (hwf : DfsOut.WF anc v (.tree e cls child)) (hE : DfsOut.Ends g v (.tree e cls child))
    (hge : anc.length ≤ cls.lowval anc.length) :
    ∀ e', subEdges (.tree e cls child) e' → ∀ x, g.Inc e' x → x = v ∨ x ∈ child.verts := by
  intro e' he' x hx
  have hret : ∀ r ∈ (DfsOut.tree e cls child).retDepths (anc ++ []).length, anc.length ≤ r := by
    intro r hr
    rw [List.append_nil] at hr
    have h := EarDfs.lowval_eq_lmin hwf
    rw [DfsOut.cls] at h
    rw [h] at hge
    exact Nat.le_trans hge (lmin_le_of_mem hr)
  rcases bdOut_wf g anc [] v _ (by rw [List.append_nil]; exact hwf) hE hret e' he' x hx with h | h | h
  · simp at h
  · exact .inl h
  · exact .inr h

/-! ### The item frame of one out (dump-checked `siteKeptCheck`/`freshCheck`) -/

/-- The items allocated at `s₀` other than the `V`/`Q` items of `vs`/`es` keep their children and
their parents in `s`. -/
def KeptRel (g : Graph) (s₀ : WalkState) (vs es : List Nat) (s : WalkState) : Prop :=
  ∀ j, j < s₀.items.size → (∀ x ∈ vs, j ≠ vertItem x) → (∀ e ∈ es, j ≠ edgeItem g e) →
    Items.ch s.items j = Items.ch s₀.items j ∧
    ∀ p, Items.IsParent s.items p j ↔ Items.IsParent s₀.items p j

/-- Every item in the spans of the entries above `base` is a `V`/`Q` item of the tree
(`v`, the vertices and edges of `outs₀`) or was allocated at or after index `sz₀`. -/
def Fresh (g : Graph) (v : Nat) (outs₀ : List DfsOut) (sz₀ : Nat) (base : List TEntry)
    (s : WalkState) : Prop :=
  ∀ top, s.tstack = top ++ base → ∀ t ∈ top, ∀ i ∈ t.spans.1 ++ t.spans.2,
    (∃ x ∈ v :: DfsOut.vertsList outs₀, i = vertItem x) ∨
    (∃ e ∈ DfsOut.edgesList outs₀, i = edgeItem g e) ∨ sz₀ ≤ i

/-- The items `j` an out `o` of `v` may not touch: allocated at `s`, not `V v`, not a `V`/`Q` item
of `o`, not in the spans of the entries above `base`. -/
def Untouched (v : Nat) (o : DfsOut) (base : List TEntry) (s : WalkState) (j : ItemId) : Prop :=
  j < s.items.size ∧ j ≠ vertItem v ∧ (∀ x ∈ o.verts, j ≠ vertItem x) ∧
  (∀ e' ∈ o.edges, j ≠ edgeItem s.g e') ∧
  ∀ top, s.tstack = top ++ base → ∀ t ∈ top, j ∉ t.spans.1 ++ t.spans.2

/-- Item frame of one out `o` of `v` from `s` to `s'` (over the base `base`): the `Untouched` items
keep their children and parents, and every item in the spans of the entries above `base` at `s'` is
`V v`, a `V`/`Q` item of `o`, allocated after `s`, or in the spans of an entry above `base` at
`s`. -/
def OutFrame (v : Nat) (o : DfsOut) (base : List TEntry) (s s' : WalkState) : Prop :=
  s.items.size ≤ s'.items.size ∧ s'.g = s.g ∧
  (∀ j, Untouched v o base s j →
    Items.ch s'.items j = Items.ch s.items j ∧
    ∀ p, Items.IsParent s'.items p j ↔ Items.IsParent s.items p j) ∧
  (∀ top', s'.tstack = top' ++ base → ∀ t ∈ top', ∀ i ∈ t.spans.1 ++ t.spans.2,
    i = vertItem v ∨ (∃ x ∈ o.verts, i = vertItem x) ∨ (∃ e' ∈ o.edges, i = edgeItem s.g e') ∨
    s.items.size ≤ i ∨ ∃ top, s.tstack = top ++ base ∧ ∃ t ∈ top, i ∈ t.spans.1 ++ t.spans.2)

theorem subEdges_of_mem_edges {o : DfsOut} {e' : Nat} (h : e' ∈ o.edges) : subEdges o e' := by
  cases o with
  | back e dest cls => simp [DfsOut.edges] at h; simp [subEdges, DfsOut.e, h]
  | tree e cls child =>
    simp only [DfsOut.edges, List.mem_cons] at h
    rcases h with rfl | h
    · exact Or.inl rfl
    · exact Or.inr h

/-- `OutFrame` from the items whose children changed (`Mod`): no `Untouched` item is in `Mod`, the
`Mod` items did not gain or lose an `Untouched` child, and no entry above `base` is new. -/
theorem outFrame_of_mod {v : Nat} {o : DfsOut} {base : List TEntry} {s s' : WalkState}
    (Mod : ItemId → Prop) (hsz : s.items.size ≤ s'.items.size) (hg : s'.g = s.g)
    (hmod : ∀ p, ¬ Mod p → Items.ch s'.items p = Items.ch s.items p)
    (hnew : ∀ p j, Mod p → Untouched v o base s j →
      (Items.IsParent s'.items p j ↔ Items.IsParent s.items p j))
    (hPj : ∀ j, Untouched v o base s j → ¬ Mod j)
    (hsub : ∀ top', s'.tstack = top' ++ base → ∀ t ∈ top',
      ∃ top, s.tstack = top ++ base ∧ t ∈ top) :
    OutFrame v o base s s' := by
  refine ⟨hsz, hg, fun j hj => ⟨hmod j (hPj j hj), fun p => ?_⟩, fun top' htop' t ht i hi => ?_⟩
  · by_cases hp : Mod p
    · exact hnew p j hp hj
    · rw [Items.IsParent, Items.IsParent, hmod p hp]
  · obtain ⟨top, htop, ht'⟩ := hsub top' htop' t ht
    exact Or.inr (Or.inr (Or.inr (Or.inr ⟨top, htop, t, ht', hi⟩)))

theorem OutFrame.trans {v : Nat} {o : DfsOut} {base : List TEntry} {s s₁ s' : WalkState}
    (h₁ : OutFrame v o base s s₁) (h₂ : OutFrame v o base s₁ s') : OutFrame v o base s s' := by
  obtain ⟨hsz₁, hg₁, hk₁, hsp₁⟩ := h₁
  obtain ⟨hsz₂, hg₂, hk₂, hsp₂⟩ := h₂
  refine ⟨Nat.le_trans hsz₁ hsz₂, hg₂.trans hg₁, fun j hj => ?_, fun top' htop' t ht i hi => ?_⟩
  · obtain ⟨hjlt, hjv, hxv, hje, htop⟩ := hj
    have h1 := hk₁ j ⟨hjlt, hjv, hxv, hje, htop⟩
    have h2 := hk₂ j ⟨Nat.lt_of_lt_of_le hjlt hsz₁, hjv, hxv, by rw [hg₁]; exact hje,
      fun top₁ htop₁ t ht hi => by
      rcases hsp₁ top₁ htop₁ t ht j hi with h | ⟨x, hx, h⟩ | ⟨e', he', h⟩ | h | ⟨top, htop₀, t', ht', hi'⟩
      · exact hjv h
      · exact hxv x hx h
      · exact hje e' he' h
      · exact absurd h (Nat.not_le.2 hjlt)
      · exact htop top htop₀ t' ht' hi'⟩
    rw [h2.1, h1.1]
    exact ⟨rfl, fun p => (h2.2 p).trans (h1.2 p)⟩
  · rcases hsp₂ top' htop' t ht i hi with h | h | ⟨e', he', h⟩ | h | ⟨top₁, htop₁, t', ht', hi'⟩
    · exact Or.inl h
    · exact Or.inr (Or.inl h)
    · exact Or.inr (Or.inr (Or.inl ⟨e', he', by rw [← hg₁]; exact h⟩))
    · exact Or.inr (Or.inr (Or.inr (Or.inl (Nat.le_trans hsz₁ h))))
    · exact hsp₁ top₁ htop₁ t' ht' i hi'

/-- Frame of the returning back-edge push (`Q e` entry on top of `L ++ s.tstack`, children
unchanged). -/
theorem outFrame_pushBack {s s' : WalkState} {v d lv e dest n n' : Nat} {cls : OutClass}
    {base top L : List TEntry} {push dirD dir : Bool}
    (htop : s.tstack = top ++ base)
    (hL : L = if push then [⟨v, d, n, setSides dirD [vertItem v] []⟩] else [])
    (hg : s'.g = s.g) (hsz' : s'.items.size = s.items.size)
    (hch : ∀ j, Items.ch s'.items j = Items.ch s.items j)
    (hts : s'.tstack = ⟨v, lv, n', setSides dir [edgeItem s.g e] []⟩ :: (L ++ s.tstack)) :
    OutFrame v (.back e dest cls) base s s' := by
  refine ⟨hsz'.ge, hg, fun j _ => ⟨hch j, fun p => by rw [Items.IsParent, Items.IsParent, hch p]⟩,
    fun top' htop' t ht i hi => ?_⟩
  rw [hts, htop, ← List.append_assoc, ← List.cons_append] at htop'
  obtain rfl := List.append_cancel_right htop'
  simp only [List.mem_cons, List.mem_append] at ht
  rcases ht with rfl | ht | ht
  · simp only [mem_setSides, List.mem_singleton] at hi
    exact Or.inr (Or.inr (Or.inl ⟨e, by simp [DfsOut.edges], hi⟩))
  · subst hL
    cases push
    · simp at ht
    · simp only [↓reduceIte, List.mem_singleton] at ht
      subst ht
      simp only [mem_setSides, List.mem_singleton] at hi
      exact Or.inl hi
  · exact Or.inr (Or.inr (Or.inr (Or.inr ⟨top, htop, t, ht, hi⟩)))

/-- Frame of a P-merge over the two top entries `c` (single item `q`) and `a` (holding `i`): only
`p`'s children change, by items that are `q`, `i` or already below `p`. -/
theorem outFrame_mergeP {s s' : WalkState} {v lv n : Nat} {o : DfsOut} {base tl' : List TEntry}
    {c a : TEntry} {q i p : ItemId} {S : List ItemId} {dir : Bool}
    (hts : s.tstack = c :: a :: (tl' ++ base))
    (hcs : c.spans.1 ++ c.spans.2 = [q]) (hai : i ∈ a.spans.1 ++ a.spans.2)
    (hpa : p = i ∨ s.items.size ≤ p)
    (hS : ∀ r ∈ S, r = q ∨ (r = i ∧ i ≠ p) ∨ Items.IsParent s.items p r)
    (hP : ∀ p' c', Items.IsParent s'.items p' c' ↔
      Items.IsParent s.items p' c' ∨ (p' = p ∧ c' ∈ S))
    (hch : ∀ j, j ≠ p → Items.ch s'.items j = Items.ch s.items j)
    (hg : s'.g = s.g) (hsz : s.items.size ≤ s'.items.size)
    (hts' : s'.tstack = ⟨v, lv, n, setSides dir [p] []⟩ :: (tl' ++ base)) :
    OutFrame v o base s s' := by
  have hstop : s.tstack = (c :: a :: tl') ++ base := hts
  refine ⟨hsz, hg, fun j hj => ?_, fun top' htop' t ht i' hi' => ?_⟩
  · obtain ⟨hjlt, -, -, -, hT⟩ := hj
    have hjc : j ≠ q := fun h => hT _ hstop c (List.mem_cons_self ..) (by rw [hcs, h]; simp)
    have hji : j ≠ i := fun h => hT _ hstop a (by simp) (h ▸ hai)
    have hjp : j ≠ p := by
      rcases hpa with rfl | hp
      · exact hji
      · exact fun h => absurd hp (Nat.not_le.2 (h ▸ hjlt))
    refine ⟨hch j hjp, fun p' => ?_⟩
    rw [hP]
    constructor
    · rintro (h | ⟨rfl, hjS⟩)
      · exact h
      · rcases hS j hjS with h | ⟨h, -⟩ | h
        · exact absurd h hjc
        · exact absurd h hji
        · exact h
    · exact Or.inl
  · rw [hts', ← List.cons_append] at htop'
    obtain rfl := List.append_cancel_right htop'
    rcases List.mem_cons.1 ht with rfl | ht
    · simp only [mem_setSides, List.mem_singleton] at hi'
      subst hi'
      rcases hpa with rfl | hp
      · exact Or.inr (Or.inr (Or.inr (Or.inr ⟨_, hstop, a, by simp, hai⟩)))
      · exact Or.inr (Or.inr (Or.inr (Or.inl hp)))
    · exact Or.inr (Or.inr (Or.inr (Or.inr ⟨_, hstop, t, by simp [ht], hi'⟩)))

/-- One out's frame composes with the kept/fresh invariants of the enclosing tree walk. -/
theorem outFrame_step {g : Graph} {v : Nat} {o : DfsOut} {outs₀ : List DfsOut} {sz₀ : Nat}
    {base : List TEntry} {s₀ s s' : WalkState} (ho : o ∈ outs₀) (hg : s.g = g)
    (hsz₀ : sz₀ ≤ s.items.size) (hs₀ : s₀.items.size ≤ sz₀)
    (hK : KeptRel g s₀ (v :: DfsOut.vertsList outs₀) (DfsOut.edgesList outs₀) s)
    (hF : Fresh g v outs₀ sz₀ base s) (hO : OutFrame v o base s s') :
    sz₀ ≤ s'.items.size ∧ s'.g = g ∧
    KeptRel g s₀ (v :: DfsOut.vertsList outs₀) (DfsOut.edgesList outs₀) s' ∧
    Fresh g v outs₀ sz₀ base s' := by
  obtain ⟨hsz, hg', hk, hsp⟩ := hO
  have hvs : ∀ x ∈ o.verts, x ∈ v :: DfsOut.vertsList outs₀ :=
    fun x hx => List.mem_cons_of_mem _ (mem_vertsList_of_verts ho hx)
  have hes : ∀ e' ∈ o.edges, e' ∈ DfsOut.edgesList outs₀ :=
    fun e' he' => mem_subEdges_edgesList.2 ⟨o, ho, subEdges_of_mem_edges he'⟩
  refine ⟨Nat.le_trans hsz₀ hsz, hg'.trans hg, fun j hj hjv hje => ?_, ?_⟩
  · have hj' := hK j hj hjv hje
    have hjlt : j < s.items.size := by omega
    have h := hk j ⟨hjlt, hjv v (List.mem_cons_self ..), fun x hx => hjv x (hvs x hx),
      fun e' he' => by rw [hg]; exact hje e' (hes e' he'),
      fun top htop t ht hi => by
        rcases hF top htop t ht j hi with ⟨x, hx, hjx⟩ | ⟨e', he', hje'⟩ | hle
        · exact hjv x hx hjx
        · exact hje e' he' hje'
        · omega⟩
    rw [h.1, hj'.1]
    refine ⟨rfl, fun p => (h.2 p).trans (hj'.2 p)⟩
  · intro top' htop' t ht i hi
    rcases hsp top' htop' t ht i hi with h | ⟨x, hx, h⟩ | ⟨e', he', h⟩ | h | ⟨top, htop, t', ht', hi'⟩
    · exact Or.inl ⟨v, List.mem_cons_self .., h⟩
    · exact Or.inl ⟨x, hvs x hx, h⟩
    · exact Or.inr (Or.inl ⟨e', hes e' he', by rw [← hg]; exact h⟩)
    · exact Or.inr (Or.inr (Nat.le_trans hsz₀ h))
    · exact hF top htop t' ht' i hi'


/-- Frame of a returning tree out from `RetTop`: the items are `sX`'s, the stack is `T` (spans in
`R`'s or `V v`) over `L ++ s.tstack`. -/
theorem outFrame_retTop {v d e y : Nat} {cls : OutClass} {outs : List DfsOut}
    {L T base top R : List TEntry} {hv₀ : Bool} {s sX s' : WalkState}
    (HR : RetTop v d s e cls y outs L hv₀ sX R) (htop : s.tstack = top ++ base)
    (hLv : ∀ t ∈ L, ∀ i ∈ t.spans.1 ++ t.spans.2, i = vertItem v)
    (hTs : ∀ t ∈ T, ∀ i ∈ t.spans.1 ++ t.spans.2,
      i = vertItem v ∨ ∃ t' ∈ R, i ∈ t'.spans.1 ++ t'.spans.2)
    (hts : s'.tstack = T ++ (L ++ s.tstack)) (hg : s'.g = sX.g) (hit : s'.items = sX.items) :
    OutFrame v (.tree e cls (.node y outs)) base s s' := by
  refine ⟨by rw [hit]; exact HR.size, hg.trans HR.g, fun j hj => ?_, fun top' htop' t ht i hi => ?_⟩
  · obtain ⟨hjlt, -, hxv, hje, -⟩ := hj
    rw [hit]
    exact HR.kept j hjlt hxv (fun e' he' => hje e' (List.mem_cons_of_mem _ he'))
      (hje e (List.mem_cons_self ..))
  · rw [hts, htop, ← List.append_assoc, ← List.append_assoc] at htop'
    obtain rfl := List.append_cancel_right htop'
    simp only [List.mem_append] at ht
    rcases ht with (ht | ht) | ht
    · rcases hTs t ht i hi with h | ⟨t', ht', hi'⟩
      · exact Or.inl h
      · by_cases hlt : i < s.items.size
        · rcases (HR.span_new t' ht' i hi').2 hlt with ⟨x, hx, h⟩ | ⟨e', he', h⟩ | h
          · exact Or.inr (Or.inl ⟨x, hx, h⟩)
          · exact Or.inr (Or.inr (Or.inl ⟨e', List.mem_cons_of_mem _ he', h⟩))
          · exact Or.inr (Or.inr (Or.inl ⟨e, List.mem_cons_self .., h⟩))
        · exact Or.inr (Or.inr (Or.inr (Or.inl (Nat.le_of_not_lt hlt))))
    · exact Or.inl (hLv t ht i hi)
    · exact Or.inr (Or.inr (Or.inr (Or.inr ⟨top, htop, t, ht, hi⟩)))

/-- Frame of the child's walk at a tree out, up to the end push: the child's kept items and the
freshness of its top entries. -/
theorem outFrame_child {v e y : Nat} {cls : OutClass} {outs : List DfsOut}
    {L L' topc base top : List TEntry} {s sE : WalkState} {D₃ : Array Bool}
    (htop : s.tstack = top ++ base) (htopE : sE.tstack = topc ++ (L ++ s.tstack))
    (hsz : s.items.size ≤ sE.items.size) (hg : sE.g = s.g)
    (hk : KeptRel s.g s (DfsTree.node y outs).verts (DfsTree.node y outs).edges sE)
    (hfr : ∀ t ∈ topc, ∀ i ∈ t.spans.1 ++ t.spans.2,
      (∃ x ∈ (DfsTree.node y outs).verts, i = vertItem x) ∨
      (∃ e' ∈ (DfsTree.node y outs).edges, i = edgeItem s.g e') ∨ s.items.size ≤ i)
    (hLv : ∀ t ∈ L, ∀ i ∈ t.spans.1 ++ t.spans.2, i = vertItem v)
    (hL'v : ∀ t ∈ L', ∀ i ∈ t.spans.1 ++ t.spans.2, i = vertItem y) :
    OutFrame v (.tree e cls (.node y outs)) base s (pushEnd sE D₃ L') := by
  refine ⟨hsz, hg, fun j hj => ?_, fun top' htop' t ht i hi => ?_⟩
  · obtain ⟨hjlt, -, hxv, hje, -⟩ := hj
    exact hk j hjlt hxv (fun e' he' => hje e' (List.mem_cons_of_mem _ he'))
  · simp only [pushEnd] at htop'
    rw [htopE, htop, ← List.append_assoc, ← List.append_assoc, ← List.append_assoc] at htop'
    obtain rfl := List.append_cancel_right htop'
    simp only [List.mem_append] at ht
    rcases ht with ((ht | ht) | ht) | ht
    · exact Or.inr (Or.inl ⟨y, by simp [DfsOut.verts, DfsTree.verts], hL'v t ht i hi⟩)
    · rcases hfr t ht i hi with ⟨x, hx, h⟩ | ⟨e', he', h⟩ | h
      · exact Or.inr (Or.inl ⟨x, hx, h⟩)
      · exact Or.inr (Or.inr (Or.inl ⟨e', List.mem_cons_of_mem _ he', h⟩))
      · exact Or.inr (Or.inr (Or.inr (Or.inl h)))
    · exact Or.inl (hLv t ht i hi)
    · exact Or.inr (Or.inr (Or.inr (Or.inr ⟨top, htop, t, ht, hi⟩)))

/-- Frame of `finishBoundary` at a self-loop: the items only gain `Q e` under `V v` and the fresh
`O` item under `Q e`; the tstack is unchanged. -/
theorem outFrame_selfLoop {s s' : WalkState} {v e dest : Nat} {cls : OutClass} {base : List TEntry}
    (hv : v < s.g.nv) (he : e < s.g.ne) (hsz : 1 + s.g.nv + s.g.ne ≤ s.items.size)
    (hq : Items.ch s.items (edgeItem s.g e) = [])
    (hg : s'.g = s.g) (hts : s'.tstack = s.tstack) (hsz' : s'.items.size = s.items.size + 1)
    (hch : ∀ j, Items.ch s'.items j =
      if j = vertItem v then Items.ch s.items (vertItem v) ++ [edgeItem s.g e]
      else if j = edgeItem s.g e then [s.items.size]
      else if j = s.items.size then [] else Items.ch s.items j) :
    OutFrame v (.back e dest cls) base s s' := by
  refine outFrame_of_mod (fun p => p = vertItem v ∨ p = edgeItem s.g e ∨ p = s.items.size)
    (by omega) hg ?_ ?_ ?_ (fun top' htop' t ht => ⟨top', hts ▸ htop', ht⟩)
  · intro p hp
    rw [hch p]
    simp only [not_or] at hp
    simp [hp.1, hp.2.1, hp.2.2]
  · intro p j hp hj
    obtain ⟨hjlt, hjv, -, hje, -⟩ := hj
    have hjq : j ≠ edgeItem s.g e := hje e (by simp [DfsOut.edges])
    have hvq := vertItem_ne_edgeItem' (g := s.g) (e := e) hv
    have h1 : s.items.size ≠ vertItem v := by show _ ≠ 1 + v; omega
    have h2 : s.items.size ≠ edgeItem s.g e := by show _ ≠ 1 + s.g.nv + e; omega
    simp only [Items.IsParent]
    rw [hch p]
    rcases hp with rfl | rfl | rfl
    · simp [hjq]
    · simp [hq, Nat.ne_of_lt hjlt, hvq.symm]
    · simp [Items.ch_of_le s.items _ (Nat.le_refl _), h1, h2]
  · intro j hj
    obtain ⟨hjlt, hjv, -, hje, -⟩ := hj
    simp only [not_or]
    exact ⟨hjv, hje e (by simp [DfsOut.edges]), Nat.ne_of_lt hjlt⟩

/-- A child walk never gives a parent to the `Q` item of an edge outside the subtree that is a root
before it. -/
theorem walkTree_qroot_kept (t : DfsTree) (d : Nat) (s : WalkState)
    (hk : wp (walkTree t d) (fun _ s' => KeptRel s.g s t.verts t.edges s') s) (e : Nat) (he : e ∉ t.edges)
    (he' : e < s.g.ne) (hsz : 1 + s.g.nv + s.g.ne ≤ s.items.size) (hvne : ∀ x ∈ t.verts, x < s.g.nv)
    (hr : ∀ p, ¬ Items.IsParent s.items p (edgeItem s.g e)) :
    wp (walkTree t d) (fun _ s' => ∀ p, ¬ Items.IsParent s'.items p (edgeItem s.g e)) s := by
  refine wp_mono _ hk fun _ s' h p hp => ?_
  refine hr p (((h (edgeItem s.g e) (by show 1 + s.g.nv + e < _; omega) ?_ ?_).2 p).1 hp)
  · exact fun x hx h => (vertItem_ne_edgeItem' (hvne x hx)) h.symm
  · exact fun e₁ he₁ h => he ((edgeItem_inj h) ▸ he₁)

/-- A child walk never gives a parent to the `V` item of a vertex outside the subtree that is a root
before it. -/
theorem walkTree_vroot_kept (t : DfsTree) (d : Nat) (s : WalkState)
    (hk : wp (walkTree t d) (fun _ s' => KeptRel s.g s t.verts t.edges s') s) (v : Nat) (hv : v ∉ t.verts)
    (hv' : v < s.g.nv) (hsz : 1 + s.g.nv + s.g.ne ≤ s.items.size)
    (hr : ∀ p, ¬ Items.IsParent s.items p (vertItem v)) :
    wp (walkTree t d) (fun _ s' => ∀ p, ¬ Items.IsParent s'.items p (vertItem v)) s := by
  refine wp_mono _ hk fun _ s' h p hp => ?_
  refine hr p (((h (vertItem v) (by show 1 + v < _; omega) ?_ ?_).2 p).1 hp)
  · exact fun x hx h => hv ((vertItem_inj' h) ▸ hx)
  · exact fun e₁ _ h => vertItem_ne_edgeItem' hv' h

theorem getElem!_range_map {β : Type} [Inhabited β] (n : Nat) (f : Nat → β) (k : Nat) (hk : k < n) :
    ((List.range n).map f)[k]! = f k := by
  rw [getElem!_pos ((List.range n).map f) k (by simpa using hk), List.getElem_map, List.getElem_range]

/-- The context at the entry of the root's `walkOuts` (`ctxCheck` at the root's `iOuts` entry,
`done = []`). -/
theorem ctx_init_root (t : DfsTree) (s : WalkState) (hwf : t.WF []) (hends : t.Ends s.g)
    (hvlt : ∀ v ∈ t.verts, v < s.g.nv) (helt : ∀ e ∈ t.edges, e < s.g.ne)
    (hvn : t.verts.Nodup) (hen : t.edges.Nodup)
    (hsv : s.stackVerts.size = s.g.nv) (hsd : s.stackDir.size = s.g.nv)
    (hfo : s.firstOccurrence.size = s.g.nv)
    (hts : s.tstack = []) (hi : s.Inv' 0) (hs : Shape s)
    (hvfresh : ∀ v ∈ t.verts, Items.ch s.items (vertItem v) = [] ∧
      ∀ p, ¬ Items.IsParent s.items p (vertItem v))
    (hefresh : ∀ e ∈ t.edges, Items.ch s.items (edgeItem s.g e) = [] ∧
      ∀ p, ¬ Items.IsParent s.items p (edgeItem s.g e)) :
    ∃ sv sd, EarCtx t.v 0 [] t.outs false [] [] sv sd { s with stackVerts := s.stackVerts.set! 0 t.v } := by
  obtain ⟨v, outs⟩ := t
  have hv : v < s.g.nv := hvlt v (List.mem_cons_self ..)
  have hsv0 : (s.stackVerts.set! 0 v)[0]! = v := getElem!_set!_self' s.stackVerts 0 v (by omega)
  have hfr := hvfresh v (List.mem_cons_self ..)
  have hnil : ∀ t, t ∉ ({ s with stackVerts := s.stackVerts.set! 0 v } : WalkState).tstack :=
    fun t ht => by simp [hts] at ht
  refine ⟨[v], [], ?_⟩
  exact
    { top := ⟨[], by simp [hts],
        { split := ⟨[], [], by simp, fun t ht => (nomatch ht), fun t ht => (nomatch ht), List.Pairwise.nil,
            fun t ht => (nomatch ht), fun t ht => (nomatch ht), fun t ht => (nomatch ht)⟩
          bot := fun t ht => (nomatch ht), vitems := fun t ht => (nomatch ht)
          qitems := fun t ht => (nomatch ht), edges := fun t ht => (nomatch ht)
          cover := fun o ho => (nomatch ho), ret := fun h => absurd rfl h }⟩
      base_edges := ⟨rfl, fun k hk => absurd hk (by simp)⟩
      base_bot := fun t ht => (nomatch ht)
      noVert_after := fun _ => rfl
      vfirst := fun t ht _ _ => by simp [hts] at ht
      sv := fun k hk => by
        obtain rfl : k = 0 := (by omega)
        show (s.stackVerts.set! 0 v)[0]! = [v][0]!
        rw [hsv0]; simp
      sd := fun k hk => absurd hk (Nat.not_lt_zero k)
      sv_d := hsv0
      path := fun k k' h h' => absurd (Nat.lt_of_lt_of_le h h') (Nat.not_lt_zero k)
      v_root := hfr.2
      vert_free := fun t ht => absurd ht (hnil t)
      afterVert_ret := fun o ho => (nomatch ho)
      hv_ret := fun h => (nomatch h)
      vert_book := fun _ => ⟨fun e _ _ _ hE _ => (edgeBelow_vert_nil hv hfr.1 e hE).elim,
        fun _ e _ _ _ hE _ _ _ => (edgeBelow_vert_nil hv hfr.1 e hE).elim⟩
      vert_disj := fun _ t ht => absurd ht (hnil t)
      vert_edges := fun e _ => ⟨fun h => (edgeBelow_vert_nil hv hfr.1 e h).elim,
        fun ⟨_, ho, _⟩ => (nomatch ho)⟩
      touch_bot := fun t ht => absurd ht (hnil t)
      span_root := fun t ht => absurd ht (hnil t)
      span_lt := fun t ht => absurd ht (hnil t)
      ch_lt := fun p c h => hs.ch_lt p c h
      disj := by simp [hts]
      span_disj := by simp [hts]
      q_fresh := fun o ho e he => by
        have := hefresh e (mem_subEdges_edgesList.2 ⟨o, ho, he⟩)
        exact ⟨this.2, this.1, fun t ht => absurd ht (hnil t)⟩
      v_fresh := fun o ho e cls child ho' y hy => by
        have hyv : y ∈ DfsOut.vertsList outs := mem_vertsList.2 ⟨o, ho, e, cls, child, ho', hy⟩
        have hf := hvfresh y (List.mem_cons_of_mem _ hyv)
        refine ⟨hf.2, hf.1, fun t ht => absurd ht (hnil t), fun k hk => ?_,
          fun t ht => absurd ht (hnil t), fun t ht => absurd ht (hnil t)⟩
        obtain rfl : k = 0 := by omega
        show (s.stackVerts.set! 0 v)[0]! ≠ y
        rw [hsv0]
        exact fun h => (List.nodup_cons.1 hvn).1 (h ▸ hyv) }

/-- The context at the entry of the child's `walkOuts` (`ctxCheck` at a child's `iOuts` entry,
`done = []`), from the parent's context at the tree-edge site (after `walkOutPre`, the
`firstOccurrence` reset and the child's `stackVerts` write). The base is the whole parent tstack
(with the parent's own `V` entry when `walkOutPre` pushed it); `hinc`/`hnc` (the done outs' edges
touch `v` and avoid the child's vertices) make that entry attached at `v` and foreign to the child. -/
theorem ctx_init_child {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s : WalkState} {e : Nat} {cls : OutClass} {y : Nat} {outs : List DfsOut}
    (hC : EarCtx v d done (.tree e cls (.node y outs) :: rest) hasVert base bE sv sd s)
    (hv : v < s.g.nv) (hy : y < s.g.nv) (hsd : d < s.stackDir.size) (hsv : d + 1 < s.stackVerts.size)
    (hinc : ∀ o' ∈ done, o'.1.e < s.g.ne ∧ s.g.Inc o'.1.e v)
    (hnc : ∀ e', e' < s.g.ne → Items.EdgeBelow s.g s.items (vertItem v) e' →
      ∀ x ∈ (DfsTree.node y outs).verts, ¬ s.g.Inc e' x)
    (hndc : (DfsTree.node y outs).verts.Nodup) (hsz : 1 + s.g.nv ≤ s.items.size)
    (push : Bool) (L : List TEntry)
    (hpush : push = true ↔ hasVert = false ∧ cls.lowval d < d ∧ cls.isType1 = true)
    (hL : L = if push then
      [⟨v, d, s.nxtEdgeIdx, setSides (s.stackDir.set! d (if cls.lowval d ≥ d then false
        else !s.stackDir[cls.lowval d]!))[d]! [vertItem v] []⟩] else []) :
    ∃ sv' sd', EarCtx y (d + 1) [] outs false (L ++ s.tstack)
      ((L ++ s.tstack).map fun t e' => t.edges s.g s.items e') sv' sd'
      { s with stackDir := s.stackDir.set! d (if cls.lowval d ≥ d then false
                 else !s.stackDir[cls.lowval d]!),
               tstack := L ++ s.tstack,
               firstOccurrence := s.firstOccurrence.set! d s.g.ne,
               stackVerts := s.stackVerts.set! (d + 1) y } := by
  have hfr := hC.v_fresh _ (List.mem_cons_self ..) e cls (.node y outs) rfl
  have hfy := hfr y (List.mem_cons_self ..)
  have hvy : v ≠ y := fun h => hfy.2.2.2.1 d (Nat.le_refl _) (hC.sv_d.trans h)
  have hsvk : ∀ k, k ≤ d → (s.stackVerts.set! (d + 1) y)[k]! = s.stackVerts[k]! :=
    fun k hk => getElem!_set!_ne' _ _ _ _ (by omega)
  have hsvd : (s.stackVerts.set! (d + 1) y)[d + 1]! = y := getElem!_set!_self' _ _ _ hsv
  have hLmem : ∀ t ∈ L, t = ⟨v, d, s.nxtEdgeIdx, setSides (s.stackDir.set! d (if cls.lowval d ≥ d
      then false else !s.stackDir[cls.lowval d]!))[d]! [vertItem v] []⟩ ∧ hasVert = false := by
    intro t ht
    subst hL
    split at ht
    · next hp => simp only [List.mem_singleton] at ht; exact ⟨ht, (hpush.1 hp).1⟩
    · simp at ht
  have hLspans : ∀ t ∈ L, t.spans.1 ++ t.spans.2 = [vertItem v] := fun t ht => by
    rw [(hLmem t ht).1]; cases (s.stackDir.set! d (if cls.lowval d ≥ d
      then false else !s.stackDir[cls.lowval d]!))[d]! <;> simp [setSides]
  have hLedges : ∀ t ∈ L, ∀ e', TEntry.edges s.g s.items t e' ↔
      Items.EdgeBelow s.g s.items (vertItem v) e' := fun t ht e' => by
    rw [(hLmem t ht).1]; exact TEntry.edges_single_entry ..
  have hLv : ∀ t ∈ L, t.vStart = v := fun t ht => by rw [(hLmem t ht).1]
  have hL1 : ∀ R : TEntry → TEntry → Prop, L.Pairwise R := fun R => by
    subst hL; split <;> simp
  refine ⟨(List.range (d + 2)).map fun k => (s.stackVerts.set! (d + 1) y)[k]!,
    (List.range (d + 1)).map fun k => (s.stackDir.set! d (if cls.lowval d ≥ d then false
      else !s.stackDir[cls.lowval d]!))[k]!, ?_⟩
  exact
    { top := ⟨[], rfl,
        { split := ⟨[], [], by simp, fun t ht => (nomatch ht), fun t ht => (nomatch ht), List.Pairwise.nil,
            fun t ht => (nomatch ht), fun t ht => (nomatch ht), fun t ht => (nomatch ht)⟩
          bot := fun t ht => (nomatch ht), vitems := fun t ht => (nomatch ht)
          qitems := fun t ht => (nomatch ht), edges := fun t ht => (nomatch ht)
          cover := fun o ho => (nomatch ho), ret := fun h => absurd rfl h }⟩
      base_edges := ⟨by simp, fun k hk e' _ => by rw [getElem!_map_fn _ _ k hk]⟩
      base_bot := fun t ht => by
        rcases List.mem_append.1 ht with h | h
        · rw [hLv t h]; exact hvy
        · exact hfy.2.2.2.2.1 t h
      noVert_after := fun _ => rfl
      vfirst := fun t ht hvt _ => absurd hvt (by
        rcases List.mem_append.1 ht with h | h
        · rw [hLv t h]; exact hvy
        · exact hfy.2.2.2.2.1 t h)
      sv := fun k hk => (getElem!_range_map _ _ k (by omega)).symm
      sd := fun k hk => (getElem!_range_map _ _ k hk).symm
      sv_d := hsvd
      path := fun k k' hkk' hk' => by
        show (s.stackVerts.set! (d + 1) y)[k]! ≠ (s.stackVerts.set! (d + 1) y)[k']!
        rcases Nat.lt_or_ge k' (d + 1) with h | h
        · rw [hsvk k (by omega), hsvk k' (by omega)]; exact hC.path k k' hkk' (by omega)
        · obtain rfl : k' = d + 1 := by omega
          rw [hsvk k (by omega), hsvd]; exact hfy.2.2.2.1 k (by omega)
      v_root := hfy.1
      vert_free := fun t ht hmem => by
        rcases List.mem_append.1 ht with h | h
        · rw [hLspans t h] at hmem
          have h' : (1 + y : Nat) = 1 + v := List.mem_singleton.1 hmem
          exact absurd (by omega) hvy.symm
        · exact absurd hmem (hfy.2.2.1 t h)
      afterVert_ret := fun o ho => (nomatch ho)
      hv_ret := fun h => (nomatch h)
      vert_book := fun _ => ⟨fun e _ _ _ hE _ => (edgeBelow_vert_nil hy hfy.2.1 e hE).elim,
        fun _ e _ _ _ hE _ _ _ => (edgeBelow_vert_nil hy hfy.2.1 e hE).elim⟩
      vert_disj := fun _ _ _ e' _ _ hE => edgeBelow_vert_nil hy hfy.2.1 e' hE
      vert_edges := fun e' _ => ⟨fun h => (edgeBelow_vert_nil hy hfy.2.1 e' h).elim,
        fun ⟨_, ho, _⟩ => (nomatch ho)⟩
      touch_bot := fun t ht ⟨e', he', hte⟩ => by
        rcases List.mem_append.1 ht with h | h
        · obtain ⟨o, ho, hlo, -⟩ := (hC.vert_edges e' he').1 ((hLedges t h e').1 hte)
          rw [hLv t h]
          exact ⟨o.1.e, (hinc o ho).1, (hLedges t h _).2
            ((hC.vert_edges _ (hinc o ho).1).2 ⟨o, ho, hlo, subEdges_e _⟩), (hinc o ho).2⟩
        · exact hC.touch_bot t h ⟨e', he', hte⟩
      span_root := fun t ht i hi => by
        rcases List.mem_append.1 ht with h | h
        · rw [hLspans t h] at hi; rw [List.mem_singleton.1 hi]; exact hC.v_root
        · exact hC.span_root t h i hi
      span_lt := fun t ht i hi => by
        rcases List.mem_append.1 ht with h | h
        · rw [hLspans t h] at hi; rw [List.mem_singleton.1 hi]; show 1 + v < s.items.size; omega
        · exact hC.span_lt t h i hi
      ch_lt := hC.ch_lt
      disj := List.pairwise_append.2 ⟨hL1 _, hC.disj, fun t ht t' ht' e' he' hte hte' =>
        hC.vert_disj (hLmem t ht).2 t' ht' e' he' hte' ((hLedges t ht e').1 hte)⟩
      span_disj := List.pairwise_append.2 ⟨hL1 _, hC.span_disj, fun t ht t' ht' i hi hi' => by
        rw [hLspans t ht] at hi
        rw [List.mem_singleton.1 hi] at hi'
        have := hC.vert_free t' ht' hi'
        rw [(hLmem t ht).2] at this
        exact Bool.false_ne_true this⟩
      q_fresh := fun o ho e' he' => by
        have hq := hC.q_fresh _ (List.mem_cons_self ..) e'
          (Or.inr (mem_subEdges_edgesList.2 ⟨o, ho, he'⟩))
        refine ⟨hq.1, hq.2.1, fun t ht => ?_⟩
        rcases List.mem_append.1 ht with h | h
        · rw [hLspans t h]
          intro hm
          have h' : (1 + s.g.nv + e' : Nat) = 1 + v := List.mem_singleton.1 hm
          omega
        · exact hq.2.2 t h
      v_fresh := fun o ho e₁ cls₁ child₁ ho₁ y' hy' => by
        have hy'o : y' ∈ DfsOut.vertsList outs := mem_vertsList.2 ⟨o, ho, e₁, cls₁, child₁, ho₁, hy'⟩
        have hy'v : y' ∈ (DfsTree.node y outs).verts := List.mem_cons_of_mem _ hy'o
        have hf := hfr y' hy'v
        refine ⟨hf.1, hf.2.1, fun t ht => ?_, fun k hk => ?_, fun t ht => ?_, fun t ht => ?_⟩
        · rcases List.mem_append.1 ht with h | h
          · rw [hLspans t h]
            intro hm
            have h' : (1 + y' : Nat) = 1 + v := List.mem_singleton.1 hm
            exact hf.2.2.2.1 d (Nat.le_refl _) (hC.sv_d.trans (by omega))
          · exact hf.2.2.1 t h
        · show (s.stackVerts.set! (d + 1) y)[k]! ≠ y'
          rcases Nat.lt_or_ge k (d + 1) with h | h
          · rw [hsvk k (by omega)]; exact hf.2.2.2.1 k (by omega)
          · obtain rfl : k = d + 1 := by omega
            rw [hsvd]
            exact fun h => (List.nodup_cons.1 hndc).1 (h ▸ hy'o)
        · rcases List.mem_append.1 ht with h | h
          · rw [hLv t h]; exact fun h' => hf.2.2.2.1 d (Nat.le_refl _) (hC.sv_d.trans h')
          · exact hf.2.2.2.2.1 t h
        · rcases List.mem_append.1 ht with h | h
          · rintro ⟨e', he', hte, hinc'⟩
            exact hnc e' he' ((hLedges t h e').1 hte) y' hy'v hinc'
          · exact hf.2.2.2.2.2 t h }

/-- `finishBoundary` at a back-edge (self-loop) site, over an abstract post-state: the items gain
`edgeItem e` under `vertItem v` and the fresh `O` item `s.items.size` under `edgeItem e`; the
stack, path and everything else are unchanged. -/
theorem earCtx_selfLoop {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s s' : WalkState} {e dest : Nat} {cls : OutClass}
    (hC : EarCtx v d done (.back e dest cls :: rest) false base bE sv sd s)
    (hv : v < s.g.nv) (he : e < s.g.ne) (hge : d ≤ cls.lowval d)
    (hsz : 1 + s.g.nv + s.g.ne ≤ s.items.size) (hend : s.g.edges[e]! = (v, v))
    (hinc : ∀ o' ∈ done, o'.1.e < s.g.ne ∧ s.g.Inc o'.1.e v)
    (hrest_e : ∀ o' ∈ rest, ∀ e', subEdges o' e' → e' ≠ e ∧ e' < s.g.ne)
    (hrest_v : ∀ o' ∈ rest, ∀ e₁ cls₁ child, o' = .tree e₁ cls₁ child →
      ∀ y ∈ child.verts, y < s.g.nv)
    (hg : s'.g = s.g) (hts : s'.tstack = s.tstack) (hsv : s'.stackVerts = s.stackVerts)
    (hsd : ∀ k, k < d → s'.stackDir[k]! = s.stackDir[k]!)
    (hsz' : s'.items.size = s.items.size + 1)
    (hch : ∀ j, Items.ch s'.items j =
      if j = vertItem v then Items.ch s.items (vertItem v) ++ [edgeItem s.g e]
      else if j = edgeItem s.g e then [s.items.size]
      else if j = s.items.size then [] else Items.ch s.items j) :
    EarCtx v d (done ++ [(.back e dest cls, false)]) rest false base bE sv sd s' := by
  have hqf := hC.q_fresh _ (List.mem_cons_self ..) e (Or.inl rfl)
  have hvq : vertItem v ≠ edgeItem s.g e := vertItem_ne_edgeItem' hv
  have hvn : vertItem v ≠ s.items.size := by show 1 + v ≠ _; omega
  have hqn : edgeItem s.g e ≠ s.items.size := by show 1 + s.g.nv + e ≠ _; omega
  have hchn : Items.ch s.items s.items.size = [] := by simp [Items.ch]
  have hP : ∀ p c, Items.IsParent s'.items p c ↔
      Items.IsParent s.items p c ∨ (p = vertItem v ∧ c = edgeItem s.g e) ∨
        (p = edgeItem s.g e ∧ c = s.items.size) := by
    intro p c
    simp only [Items.IsParent, hch]
    by_cases h1 : p = vertItem v
    · subst h1; simp [hvq]
    · by_cases h2 : p = edgeItem s.g e
      · subst h2; simp [h1, hqf.2.1]
      · by_cases h3 : p = s.items.size
        · subst h3; simp [h1, h2, hchn]
        · simp [h1, h2, h3]
  have hmono : ∀ {a i}, Items.Below s.items a i → Items.Below s'.items a i := by
    intro a i h
    induction h with
    | refl => exact .refl
    | tail _ hs ih => exact ih.tail ((hP _ _).2 (Or.inl hs))
  have hbelow : ∀ a i, ¬ Items.Below s.items a (vertItem v) → ¬ Items.Below s.items a (edgeItem s.g e) →
      (Items.Below s'.items a i ↔ Items.Below s.items a i) := by
    intro a i hav haq
    refine ⟨fun h => ?_, hmono⟩
    induction h with
    | refl => exact .refl
    | tail _ hji ih =>
      rcases (hP _ _).1 hji with h | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩
      · exact ih.tail h
      · exact absurd ih hav
      · exact absurd ih haq
  have hbelow_v : ∀ i, Items.Below s'.items (vertItem v) i ↔
      Items.Below s.items (vertItem v) i ∨ i = edgeItem s.g e ∨ i = s.items.size := by
    intro i
    constructor
    · intro h
      induction h with
      | refl => exact .inl .refl
      | tail _ hji ih =>
        rcases (hP _ _).1 hji with h | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩
        · rcases ih with ih | rfl | rfl
          · exact .inl (ih.tail h)
          · exact absurd h (by rw [Items.IsParent, hqf.2.1]; simp)
          · exact absurd h (by rw [Items.IsParent, hchn]; simp)
        · exact .inr (.inl rfl)
        · exact .inr (.inr rfl)
    · rintro (h | rfl | rfl)
      · exact hmono h
      · exact Relation.ReflTransGen.single ((hP _ _).2 (.inr (.inl ⟨rfl, rfl⟩)))
      · exact (Relation.ReflTransGen.single ((hP _ _).2 (.inr (.inl ⟨rfl, rfl⟩)))).tail
          ((hP _ _).2 (.inr (.inr ⟨rfl, rfl⟩)))
  have hEB_v : ∀ e', e' < s.g.ne → (Items.EdgeBelow s.g s'.items (vertItem v) e' ↔
      Items.EdgeBelow s.g s.items (vertItem v) e' ∨ e' = e) := by
    intro e' he'
    rw [Items.EdgeBelow, Items.EdgeBelow, hbelow_v]
    have hne : edgeItem s.g e' ≠ s.items.size := by show 1 + s.g.nv + e' ≠ _; omega
    constructor
    · rintro (h | h | h)
      · exact .inl h
      · exact .inr (edgeItem_inj h)
      · exact absurd h hne
    · rintro (h | rfl)
      · exact .inl h
      · exact .inr (.inl rfl)
  have hspanB : ∀ t ∈ s.tstack, ∀ i ∈ t.spans.1 ++ t.spans.2, ∀ j,
      (Items.Below s'.items i j ↔ Items.Below s.items i j) := by
    intro t ht i hi j
    apply hbelow
    · intro h
      have := Items.Below.eq_of_no_parent hC.v_root h
      subst this
      exact Bool.false_ne_true (hC.vert_free t ht hi)
    · intro h
      have := Items.Below.eq_of_no_parent hqf.1 h
      subst this
      exact hqf.2.2 t ht hi
  have hedges : ∀ t ∈ s.tstack, ∀ e',
      (TEntry.edges s.g s'.items t e' ↔ TEntry.edges s.g s.items t e') := fun t ht e' =>
    TEntry.edges_congr (fun i hi e' => hspanB t ht i hi _) e'
  have hAV : afterVert (done ++ [(.back e dest cls, false)]) = afterVert done := by simp [afterVert]
  have hVL : ∀ x, x ∈ DfsOut.vertsList (done.map (·.1)) →
      x ∈ DfsOut.vertsList ((done ++ [((DfsOut.back e dest cls, false) : DfsOut × Bool)]).map (·.1)) := by
    intro x hx
    rw [DfsOut.vertsList_eq] at hx ⊢
    rw [List.map_append, List.flatMap_append]
    exact List.mem_append_left _ hx
  obtain ⟨top, htop, hCT⟩ := hC.top
  have hmemT : ∀ t ∈ top, t ∈ s.tstack := fun t ht => by
    rw [htop]; exact List.mem_append_left _ ht
  have hbase_mem : ∀ k, k < base.length → base[k]! ∈ s.tstack := fun k hk => by
    rw [htop, getElem!_pos base k hk]; exact List.mem_append_right _ (List.getElem_mem hk)
  obtain ⟨above, below, hsp, -, -, -, -, -, hB⟩ := hCT.split
  simp only [Bool.false_eq_true, ↓reduceIte] at hsp
  obtain ⟨rfl, rfl⟩ := hsp
  have hbook := hC.vert_book rfl
  exact
    { top := ⟨top, by rw [hts, htop],
        { split := ⟨[], top, by simp, fun t ht => (nomatch ht), fun t ht => (nomatch ht),
            List.Pairwise.nil, fun t ht => (nomatch ht), fun t ht => (nomatch ht), hB⟩
          bot := fun t ht k hk => by rw [hsv]; exact hCT.bot t ht k hk
          vitems := fun t ht x hx hm =>
            (hCT.vitems t ht x (by rwa [hg] at hx) hm).imp id (hVL x)
          qitems := fun t ht e' he' hm =>
            let ⟨o, ho, hs⟩ := hCT.qitems t ht e' (by rwa [hg] at he') (by rwa [hg] at hm)
            ⟨o, List.mem_append_left _ ho, hs⟩
          edges := fun t ht e' he' hte => by
            rw [hg] at he' hte
            obtain ⟨o, ho, hs⟩ := hCT.edges t ht e' he' ((hedges t (hmemT t ht) e').1 hte)
            exact ⟨o, List.mem_append_left _ ho, hs⟩
          cover := fun o ho hlt e' he' hs => by
            rw [hg] at he' ⊢
            rcases List.mem_append.1 ho with h | h
            · obtain ⟨t, ht, hte⟩ := hCT.cover o h hlt e' he' hs
              exact ⟨t, ht, (hedges t (hmemT t ht) e').2 hte⟩
            · rw [List.mem_singleton] at h
              subst h
              exact absurd hlt (Nat.not_lt.2 hge)
          ret := fun h => (hCT.ret h).imp fun o ⟨ho, hl⟩ => ⟨List.mem_append_left _ ho, hl⟩ }⟩
      base_edges := ⟨hC.base_edges.1, fun k hk e' he' => by
        rw [hg] at he' ⊢
        exact (hedges _ (hbase_mem k hk) e').trans (hC.base_edges.2 k hk e' he')⟩
      base_bot := hC.base_bot
      noVert_after := fun _ => by
        have h := hC.noVert_after rfl
        simp only [afterVert, List.map_eq_nil_iff] at h
        simp [afterVert, List.filter_append, h]
      vfirst := fun t ht hvt _ => absurd hvt (by
        obtain ⟨top, htop, hCT⟩ := hC.top
        obtain ⟨above, below, hsplit, -, -, -, -, -, hbelow⟩ := hCT.split
        rw [if_neg Bool.false_ne_true] at hsplit
        rw [hts, htop] at ht
        rcases List.mem_append.1 ht with h | h
        · exact hbelow t (hsplit.2 ▸ h)
        · exact hC.base_bot t h)
      sv := fun k hk => by rw [hsv]; exact hC.sv k hk
      sd := fun k hk => by rw [hsd k hk]; exact hC.sd k hk
      sv_d := by rw [hsv]; exact hC.sv_d
      path := fun k k' h h' => by rw [hsv]; exact hC.path k k' h h'
      v_root := fun p h => by
        rcases (hP p _).1 h with h | ⟨_, h⟩ | ⟨_, h⟩
        · exact hC.v_root p h
        · exact hvq h
        · exact hvn h
      vert_free := fun t ht hm => by rw [hts] at ht; exact hC.vert_free t ht hm
      afterVert_ret := fun o ho => by rw [hAV] at ho; exact hC.afterVert_ret o ho
      hv_ret := fun h => (nomatch h)
      vert_book := fun _ => by
        rw [hg]
        have hadj : ∀ u w, s.g.AdjIn (Items.EdgeBelow s.g s.items (vertItem v)) u w →
            s.g.AdjIn (Items.EdgeBelow s.g s'.items (vertItem v)) u w :=
          fun u w ⟨e₁, h1, h2, h3⟩ => ⟨e₁, h1, (hEB_v e₁ h1).2 (.inl h2), h3⟩
        have hreach : ∀ e₁, e₁ < s.g.ne → Items.EdgeBelow s.g s.items (vertItem v) e₁ →
            Relation.ReflTransGen (s.g.AdjIn (Items.EdgeBelow s.g s'.items (vertItem v))) v
              (s.g.edges[e₁]!).1 := by
          intro e₁ h1 hE1
          obtain ⟨o, ho, hlo, -⟩ := (hC.vert_edges e₁ h1).1 hE1
          have he0 : Items.EdgeBelow s.g s.items (vertItem v) o.1.e :=
            (hC.vert_edges _ (hinc o ho).1).2 ⟨o, ho, hlo, subEdges_e _⟩
          have step : Relation.ReflTransGen (s.g.AdjIn (Items.EdgeBelow s.g s'.items (vertItem v))) v
              (s.g.edges[o.1.e]!).1 := by
            rcases (hinc o ho).2 with h | h
            · rw [h]
            · exact .single ⟨o.1.e, (hinc o ho).1, (hEB_v _ (hinc o ho).1).2 (.inl he0),
                .inr (by rw [h])⟩
          exact step.trans (Relation.ReflTransGen.mono hadj _ _ (hbook.1 _ _ (hinc o ho).1 h1 he0 hE1))
        refine ⟨fun e₁ e₂ h1 h2 hE1 hE2 => ?_, fun x e₁ e₂ h1 h2 hE1 hnE2 hi1 hi2 => ?_⟩
        · rcases (hEB_v e₁ h1).1 hE1 with hE1 | hE1 <;> rcases (hEB_v e₂ h2).1 hE2 with hE2 | hE2
          · exact Relation.ReflTransGen.mono hadj _ _ (hbook.1 e₁ e₂ h1 h2 hE1 hE2)
          · subst hE2
            simp only [hend]
            have hsymm : ∀ a b, Relation.ReflTransGen
                (s.g.AdjIn (Items.EdgeBelow s.g s'.items (vertItem v))) a b →
                Relation.ReflTransGen (s.g.AdjIn (Items.EdgeBelow s.g s'.items (vertItem v))) b a := by
              intro a b h
              induction h with
              | refl => exact .refl
              | tail _ hs ih => exact (Relation.ReflTransGen.single (Graph.AdjIn.symm hs)).trans ih
            exact hsymm _ _ (hreach e₁ h1 hE1)
          · subst hE1
            simp only [hend]
            exact hreach e₂ h2 hE2
          · subst hE1 hE2; exact .refl
        · rcases (hEB_v e₁ h1).1 hE1 with hE1 | hE1
          · exact hbook.2 x e₁ e₂ h1 h2 hE1 (fun h => hnE2 ((hEB_v e₂ h2).2 (.inl h))) hi1 hi2
          · subst hE1
            left
            rcases hi1 with h | h <;> (rw [hend] at h; exact h.symm)
      vert_disj := fun _ t ht e' he' hte hE => by
        rw [hts] at ht
        rw [hg] at he' hte hE
        rcases (hEB_v e' he').1 hE with h | rfl
        · exact hC.vert_disj rfl t ht e' he' ((hedges t ht e').1 hte) h
        · obtain ⟨i, hi, hb⟩ := (hedges t ht e').1 hte
          exact hqf.2.2 t ht (Items.Below.eq_of_no_parent hqf.1 hb ▸ hi)
      vert_edges := fun e' he' => by
        rw [hg] at he' ⊢
        rw [hEB_v e' he', hC.vert_edges e' he']
        constructor
        · rintro (⟨o, ho, hlo, hs⟩ | h)
          · exact ⟨o, List.mem_append_left _ ho, hlo, hs⟩
          · rw [h]
            exact ⟨(DfsOut.back e dest cls, false),
              List.mem_append_right _ (List.mem_singleton_self _), hge, Or.inl rfl⟩
        · rintro ⟨o, ho, hlo, hs⟩
          rcases List.mem_append.1 ho with h | h
          · exact .inl ⟨o, h, hlo, hs⟩
          · rw [List.mem_singleton] at h
            subst h
            rcases hs with hs | hs
            · exact .inr hs
            · exact hs.elim
      touch_bot := fun t ht ⟨e', he', hte⟩ => by
        rw [hts] at ht
        rw [hg] at he' hte ⊢
        obtain ⟨e₀, h0, hE0, hi0⟩ := hC.touch_bot t ht ⟨e', he', (hedges t ht e').1 hte⟩
        exact ⟨e₀, h0, (hedges t ht e₀).2 hE0, hi0⟩
      span_root := fun t ht i hi p h => by
        rw [hts] at ht
        rcases (hP p i).1 h with h | ⟨_, rfl⟩ | ⟨_, rfl⟩
        · exact hC.span_root t ht i hi p h
        · exact hqf.2.2 t ht hi
        · exact Nat.lt_irrefl _ (hC.span_lt t ht _ hi)
      span_lt := fun t ht i hi => by
        rw [hts] at ht
        rw [hsz']
        exact Nat.lt_succ_of_lt (hC.span_lt t ht i hi)
      ch_lt := fun p c h => by
        rw [hsz']
        rcases (hP p c).1 h with h | ⟨_, rfl⟩ | ⟨_, rfl⟩
        · exact Nat.lt_succ_of_lt (hC.ch_lt p c h)
        · show 1 + s.g.nv + e < _; omega
        · exact Nat.lt_succ_self _
      disj := by
        rw [hts]
        refine List.Pairwise.imp_of_mem (fun {t t'} ht ht' h e' he' hte hte' => ?_) hC.disj
        rw [hg] at he' hte hte'
        exact h e' he' ((hedges t ht e').1 hte) ((hedges t' ht' e').1 hte')
      span_disj := by rw [hts]; exact hC.span_disj
      q_fresh := fun o ho e' hs => by
        obtain ⟨hne, hlt⟩ := hrest_e o ho e' hs
        obtain ⟨hr, hc, hm⟩ := hC.q_fresh o (List.mem_cons_of_mem _ ho) e' hs
        rw [hg]
        have hq' : edgeItem s.g e' ≠ edgeItem s.g e := fun h => hne (edgeItem_inj h)
        have hn' : edgeItem s.g e' ≠ s.items.size := by show 1 + s.g.nv + e' ≠ _; omega
        have hv' : edgeItem s.g e' ≠ vertItem v := (vertItem_ne_edgeItem' hv).symm
        refine ⟨fun p h => ?_, ?_, fun t ht => ?_⟩
        · rcases (hP p _).1 h with h | ⟨_, h⟩ | ⟨_, h⟩
          · exact hr p h
          · exact hq' h
          · exact hn' h
        · rw [hch, if_neg hv', if_neg hq', if_neg hn']; exact hc
        · rw [hts] at ht; exact hm t ht
      v_fresh := fun o ho e₁ cls₁ child ho₁ y hy => by
        have hyv := hrest_v o ho e₁ cls₁ child ho₁ y hy
        obtain ⟨hr, hc, hm, hsvk, hvs, htch⟩ :=
          hC.v_fresh o (List.mem_cons_of_mem _ ho) e₁ cls₁ child ho₁ y hy
        have hyv' : vertItem y ≠ vertItem v := fun h =>
          hsvk d (Nat.le_refl _) (hC.sv_d.trans (vertItem_inj' h).symm)
        have hyq : vertItem y ≠ edgeItem s.g e := vertItem_ne_edgeItem' hyv
        have hyn : vertItem y ≠ s.items.size := by show 1 + y ≠ _; omega
        refine ⟨fun p h => ?_, ?_, fun t ht => ?_, fun k hk => ?_, fun t ht => ?_, fun t ht => ?_⟩
        · rcases (hP p _).1 h with h | ⟨_, h⟩ | ⟨_, h⟩
          · exact hr p h
          · exact hyq h
          · exact hyn h
        · rw [hch, if_neg hyv', if_neg hyq, if_neg hyn]; exact hc
        · rw [hts] at ht; exact hm t ht
        · rw [hsv]; exact hsvk k hk
        · rw [hts] at ht; exact hvs t ht
        · rw [hts] at ht
          rw [hg]
          rintro ⟨e', he', hte, hi⟩
          exact htch t ht ⟨e', he', (hedges t ht e').1 hte, hi⟩ }

/-- `ctx_step_back`, the boundary (self-loop, `d ≤ lowval`) case: `hasVert = false` here
(`hv_ret` + rank order), nothing is pushed, and `finishBoundary` only relinks items. -/
theorem ctx_step_back_boundary {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s : WalkState} {e dest : Nat} {cls : OutClass}
    (hC : EarCtx v d done (.back e dest cls :: rest) hasVert base bE sv sd s)
    (hv : v < s.g.nv) (hsd : d < s.stackDir.size)
    (hrank : ∀ o' ∈ done, o'.1.cls.rank ≤ cls.rank)
    (hinc : ∀ o' ∈ done, o'.1.e < s.g.ne ∧ s.g.Inc o'.1.e v)
    (hnd : ∀ o' ∈ done, ∀ e', subEdges o'.1 e' → e' ≠ e)
    (hb : cls.isTree = false) (hcls : ∀ lv k, cls = .ret lv k → lv < d)
    (he : e < s.g.ne) (hsz : 1 + s.g.nv + s.g.ne ≤ s.items.size)
    (hend : Items.PairEq (v, dest) s.g.edges[e]!) (hself : d ≤ cls.lowval d → dest = v)
    (hrest_e : ∀ o' ∈ rest, ∀ e', subEdges o' e' → e' ≠ e ∧ e' < s.g.ne)
    (hrest_v : ∀ o' ∈ rest, ∀ e₁ cls₁ child, o' = .tree e₁ cls₁ child →
      ∀ y ∈ child.verts, y < s.g.nv)
    (L : List TEntry) (push : Bool)
    (hpush : push = true ↔ hasVert = false ∧ cls.lowval d < d ∧ cls.isType1 = true)
    (hL : L = if push then
      [⟨v, d, s.nxtEdgeIdx, setSides (s.stackDir.set! d (if cls.lowval d ≥ d then false
        else !s.stackDir[cls.lowval d]!))[d]! [vertItem v] []⟩] else [])
    (hge : d ≤ cls.lowval d) :
    wp (finishEdge v d (.back e dest cls) (L ++ s.tstack).length (hasVert || push))
      (fun hv' s' => EarCtx v d (done ++ [(.back e dest cls, hasVert || push)]) rest hv' base bE sv sd s' ∧
        OutFrame v (.back e dest cls) base s s' ∧
        (CtxShape v d done hasVert base s →
          CtxShape v d (done ++ [(.back e dest cls, hasVert || push)]) hv' base s'))
      { s with stackDir := s.stackDir.set! d (if cls.lowval d ≥ d then false
                 else !s.stackDir[cls.lowval d]!),
               tstack := L ++ s.tstack } := by
  have hhv : hasVert = false := by
    cases hhv : hasVert
    · rfl
    · exact absurd (hC.hv_ret hhv) (no_ret_before_boundary hge hcls hrank)
  have hpf : push = false := by
    cases hp : push
    · rfl
    · exact absurd (hpush.1 hp).2.1 (Nat.not_lt.2 hge)
  subst hhv hpf
  have hL' : L = [] := by rw [hL]; simp
  subst hL'
  have hend' : s.g.edges[e]! = (v, v) := by
    have := hself hge
    subst this
    rcases hend with h | h
    · exact h.symm
    · exact Prod.ext (congrArg Prod.snd h).symm (congrArg Prod.fst h).symm
  have hge' : cls.lowval d ≥ d := hge
  rw [finishEdge_eq]
  simp only [finishEdge', wp_bind, wp_get, wp_stackDir, DfsOut.cls, hge', ↓reduceIte]
  unfold finishBoundary
  simp only [wp_bind, wp_modifyItem, wp_modify, wp_allocItem, wp_pure, DfsOut.cls, hb,
    Bool.false_eq_true, ↓reduceIte, DfsOut.e, Bool.or_self, List.nil_append]
  refine (fun hch => ⟨earCtx_selfLoop hC hv he hge hsz hend' hinc hrest_e hrest_v rfl rfl rfl
    (fun k hk => getElem!_set!_ne' _ _ _ _ (by omega)) (by simp) hch,
    outFrame_selfLoop hv he hsz (hC.q_fresh _ (List.mem_cons_self ..) e (subEdges_e _)).2.1
      rfl rfl (by simp) hch,
    fun hS => ctxShape_same hS hge (by simp) fun k hk => getElem!_set!_ne' _ _ _ _ (Nat.ne_of_lt hk)⟩) ?_
  intro j
  have hqlt : edgeItem s.g e < s.items.size := by show 1 + s.g.nv + e < _; omega
  have hvlt : vertItem v < s.items.size := by show 1 + v < _; omega
  have hvq : vertItem v ≠ edgeItem s.g e := vertItem_ne_edgeItem' hv
  simp only [Items.ch, Array.getElem?_modify, Array.getElem?_push, Array.size_modify]
  by_cases h1 : j = vertItem v
  · subst h1; simp [hvq, hvq.symm, Array.getElem?_eq_getElem hvlt, Nat.ne_of_lt hvlt, Nat.ne_of_gt hvlt]
  · by_cases h2 : j = edgeItem s.g e
    · subst h2; simp [h1, Ne.symm h1, Array.getElem?_eq_getElem hqlt, Nat.ne_of_lt hqlt, Nat.ne_of_gt hqlt]
    · by_cases h3 : j = s.items.size
      · subst h3; simp [h1, h2, Ne.symm h1, Ne.symm h2]
      · simp [h1, h2, h3, Ne.symm h1, Ne.symm h2, Ne.symm h3]

/-- `Items.Below` from an item without children: only the item itself. -/
theorem below_of_ch_nil {items : Items} {q j : ItemId} (hq : Items.ch items q = [])
    (h : Items.Below items q j) : j = q := by
  induction h with
  | refl => rfl
  | tail _ hs ih => subst ih; exact absurd hs (by rw [Items.IsParent, hq]; simp)

theorem RetKind.rank_le_two (k : RetKind) : k.rank ≤ 2 := by cases k <;> decide

/-- Rank order of two returning classes orders their lowvals. -/
theorem lowval_le_of_rank {d : Nat} {c c' : OutClass} (h : c'.lowval d < d) (h' : c.lowval d < d)
    (hr : c'.rank ≤ c.rank) : c'.lowval d ≤ c.lowval d := by
  cases c <;> cases c' <;> simp only [OutClass.lowval, OutClass.rank] at h h' hr ⊢ <;> try omega
  rename_i lv k lv' k'
  have := RetKind.rank_le_two k; have := RetKind.rank_le_two k'; omega

theorem spans_setSides_single (dir : Bool) (i : ItemId) :
    (setSides dir [i] []).1 ++ (setSides dir [i] []).2 = [i] := by
  unfold setSides; cases dir <;> rfl

theorem mem_mergeInto_spans (u c : TEntry) (i : ItemId) :
    i ∈ (TEntry.mergeInto u c).spans.1 ++ (TEntry.mergeInto u c).spans.2 ↔
      i ∈ u.spans.1 ++ u.spans.2 ∨ i ∈ c.spans.1 ++ c.spans.2 := by
  simp only [TEntry.mergeInto, List.mem_append]
  tauto

/-- The shape after the `Q e` push (and the first-edge `V v` push) at a returning back edge, when
the P-check fails (`hnoP`). -/
theorem ctxShape_pushBack {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s s' : WalkState} {e dest : Nat} {cls : OutClass}
    (hC : EarCtx v d done (.back e dest cls :: rest) hasVert base bE sv sd s)
    (hS : CtxShape v d done hasVert base s)
    (hv : v < s.g.nv) (hlt : cls.lowval d < d)
    (hrank : ∀ o' ∈ done, o'.1.cls.rank ≤ cls.rank)
    (L : List TEntry) (push dirD : Bool)
    (hpush : push = true ↔ hasVert = false)
    (hL : L = if push then [⟨v, d, s.nxtEdgeIdx, setSides dirD [vertItem v] []⟩] else [])
    (hdirD : push = true → dirD = !s.stackDir[cls.lowval d]!)
    (hnoP : push = false → ∀ a tl, s.tstack = a :: tl → a.vStart = v → a.topDepth ≠ cls.lowval d)
    (hsd : ∀ k, k < d → s'.stackDir[k]! = s.stackDir[k]!)
    (hts : s'.tstack = ⟨v, cls.lowval d, s.nxtEdgeIdx,
      setSides s.stackDir[cls.lowval d]! [edgeItem s.g e] []⟩ :: (L ++ s.tstack)) :
    CtxShape v d (done ++ [(.back e dest cls, true)]) true base s' := by
  set lv := cls.lowval d with hlv
  set c : TEntry := ⟨v, lv, s.nxtEdgeIdx, setSides s.stackDir[lv]! [edgeItem s.g e] []⟩ with hc
  have hcV : vertItem v ∉ c.spans.1 ++ c.spans.2 := by
    show vertItem v ∉ (setSides s.stackDir[lv]! [edgeItem s.g e] []).1 ++
      (setSides s.stackDir[lv]! [edgeItem s.g e] []).2
    rw [spans_setSides_single, List.mem_singleton]; exact vertItem_ne_edgeItem' hv
  obtain ⟨top, htop, hCT⟩ := hC.top
  refine
    { ret_hv := fun _ => rfl
      t1_flag := fun o ho hl ht => by
        rcases List.mem_append.1 ho with h | h
        · exact hS.t1_flag o h hl ht
        · rw [List.mem_singleton] at h; subst h; rfl
      t1 := fun _ hT => ?_
      vt_side := fun top' htop' t ht htv hm => ?_ }
  · cases hp : push
    · have hhv : hasVert = true := by
        cases h : hasVert
        · exact absurd (hpush.2 h) (by rw [hp]; decide)
        · rfl
      simp only [hp, Bool.false_eq_true, ↓reduceIte] at hL
      rw [hL, List.nil_append] at hts
      obtain ⟨above, vt, hab, hvtv, hvtm, hpw⟩ :=
        hS.t1 hhv fun o ho hl => hT o (List.mem_append_left _ ho) hl
      refine ⟨c :: above, vt, by rw [hts, hab]; rfl, hvtv, hvtm, List.pairwise_cons.2 ⟨?_, hpw⟩⟩
      obtain ⟨above₀, below₀, hsp, hCE, -, -, hLow, -, -⟩ := hCT.split
      simp only [hhv, ↓reduceIte] at hsp
      obtain ⟨vt₀, htv₀, -, hvtm₀, -, -, -, -⟩ := hsp
      have hab₀ : s.tstack = above₀ ++ vt₀ :: (below₀ ++ base) := by rw [htop, htv₀]; simp
      obtain ⟨rfl, rfl, -⟩ := split_unique hab hab₀ hvtm hvtm₀ hC.span_disj
      have hle : ∀ t ∈ above, t.topDepth ≤ lv := fun t ht => by
        obtain ⟨o, ho, hto⟩ := hLow t ht
        have hr := hrank (o, true) (mem_afterVert ho)
        rw [hto]
        exact lowval_le_of_rank (by rw [← hto]; exact (hCE t ht).depth) hlt hr
      rcases above with _ | ⟨a, above'⟩
      · intro t ht; exact absurd ht List.not_mem_nil
      · have ha : a.topDepth < lv := Nat.lt_of_le_of_ne (hle a (List.mem_cons_self ..))
          (hnoP hp a _ (by rw [hab]; rfl) (hCE a (List.mem_cons_self ..)).vStart)
        intro t ht
        rcases List.mem_cons.1 ht with rfl | ht'
        · exact ha
        · exact Nat.lt_trans ((List.pairwise_cons.1 hpw).1 t ht') ha
    · have hhv : hasVert = false := hpush.1 hp
      have htop0 : top = [] := by
        by_contra h
        exact absurd (hS.ret_hv (hCT.ret h)) (by rw [hhv]; decide)
      simp only [hp, ↓reduceIte] at hL
      rw [hL] at hts
      refine ⟨[c], ⟨v, d, s.nxtEdgeIdx, setSides dirD [vertItem v] []⟩, ?_, rfl, ?_,
        List.pairwise_singleton _ _⟩
      · rw [hts, htop, htop0]; rfl
      · show vertItem v ∈ (setSides dirD [vertItem v] []).1 ++ (setSides dirD [vertItem v] []).2
        rw [spans_setSides_single]; exact List.mem_singleton_self _
  · have htop'' : top' = c :: (L ++ top) := by
      apply List.append_cancel_right (bs := base)
      rw [← htop', hts, htop]; simp
    rw [htop''] at ht
    rcases List.mem_cons.1 ht with rfl | ht
    · exact absurd hm hcV
    · rcases List.mem_append.1 ht with ht | ht
      · rw [hL] at ht
        split at ht
        · rw [List.mem_singleton] at ht; subst ht
          refine ⟨(.back e dest cls, true), List.mem_append_right _ (List.mem_singleton_self _),
            hlt, ?_⟩
          show setSides dirD [vertItem v] [] = setSides (!s'.stackDir[lv]!) [vertItem v] []
          rw [hdirD ‹_›, hsd _ hlt]
        · exact absurd ht List.not_mem_nil
      · obtain ⟨o, ho, hl, hsp⟩ := hS.vt_side top htop t ht htv hm
        exact ⟨o, List.mem_append_left _ ho, hl, by rw [hsd _ hl]; exact hsp⟩

/-- The shape after a P-merge: the returning type-1 out's entry `(v, lowval)` merged into the
previous ear at the same lowval (the top of `above`), leaving one entry `a'` in its place. -/
theorem ctxShape_mergeP {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s s' : WalkState} {o : DfsOut} {a a' : TEntry} {tl : List TEntry} {p : ItemId} {dir : Bool}
    (hC : EarCtx v d done (o :: rest) true base bE sv sd s)
    (hS : CtxShape v d done true base s)
    (hv : v < s.g.nv) (hlt : o.cls.lowval d < d)
    (hsz : 1 + s.g.nv + s.g.ne ≤ s.items.size)
    (hts : s.tstack = a :: tl) (hav : a.vStart = v) (had : a.topDepth = o.cls.lowval d)
    (hp : p ∈ a.spans.1 ++ a.spans.2 ∨ s.items.size ≤ p)
    (ha'v : a'.vStart = v) (ha'd : a'.topDepth = o.cls.lowval d)
    (ha's : a'.spans = setSides dir [p] [])
    (hsd : ∀ k, k < d → s'.stackDir[k]! = s.stackDir[k]!)
    (hts' : s'.tstack = a' :: tl) :
    CtxShape v d (done ++ [(o, true)]) true base s' := by
  obtain ⟨top, htop, hCT⟩ := hC.top
  obtain ⟨above₀, below₀, hsp, -, -, -, -, -, -⟩ := hCT.split
  simp only [↓reduceIte] at hsp
  obtain ⟨vt₀, htv₀, -, hvtm₀, hvtd, hvt4, -, -⟩ := hsp
  have hab₀ : s.tstack = above₀ ++ vt₀ :: (below₀ ++ base) := by rw [htop, htv₀]; simp
  obtain ⟨above', rfl⟩ : ∃ above', above₀ = a :: above' := by
    rcases above₀ with _ | ⟨a₀, above'⟩
    · exfalso
      have h := hab₀; rw [hts] at h
      rw [← (List.cons.inj h).1] at hvtd
      exact absurd (hvtd hav).1 (by rw [had]; exact Nat.ne_of_lt hlt)
    · have h := hab₀; rw [hts] at h
      exact ⟨above', by rw [(List.cons.inj h).1]⟩
  have htl : tl = above' ++ vt₀ :: (below₀ ++ base) := by
    have h := hab₀; rw [hts] at h; exact (List.cons.inj h).2
  have haV : vertItem v ∉ a.spans.1 ++ a.spans.2 := hvt4 a (List.mem_cons_self ..)
  have hpV : p ≠ vertItem v := by
    rintro rfl
    rcases hp with h | h
    · exact haV h
    · exact absurd h (Nat.not_le.2 (by show 1 + v < _; omega))
  have ha'V : vertItem v ∉ a'.spans.1 ++ a'.spans.2 := by
    rw [ha's, spans_setSides_single, List.mem_singleton]; exact fun h => hpV h.symm
  refine
    { ret_hv := fun _ => rfl
      t1_flag := fun o' ho' hl ht => by
        rcases List.mem_append.1 ho' with h | h
        · exact hS.t1_flag o' h hl ht
        · rw [List.mem_singleton] at h; subst h; rfl
      t1 := fun _ hT => ?_
      vt_side := fun top' htop' t ht htv hm => ?_ }
  · obtain ⟨above, vt, hab, hvtv, hvtm, hpw⟩ :=
      hS.t1 rfl fun o' ho' hl => hT o' (List.mem_append_left _ ho') hl
    obtain ⟨rfl, rfl, -⟩ := split_unique hab hab₀ hvtm hvtm₀ hC.span_disj
    refine ⟨a' :: above', vt, ?_, hvtv, hvtm, ?_⟩
    · have h := hab; rw [hts] at h
      rw [hts', (List.cons.inj h).2]; rfl
    · rw [List.pairwise_cons] at hpw ⊢
      refine ⟨fun t ht => ?_, hpw.2⟩
      rw [ha'd, ← had]; exact hpw.1 t ht
  · have htop'' : top' = a' :: (above' ++ vt₀ :: below₀) := by
      apply List.append_cancel_right (bs := base)
      rw [← htop', hts', htl]; simp
    rw [htop''] at ht
    rcases List.mem_cons.1 ht with rfl | ht
    · exact absurd hm ha'V
    · have ht₀ : t ∈ top := by rw [htv₀]; exact List.mem_cons_of_mem _ ht
      obtain ⟨o', ho', hl, hsp⟩ := hS.vt_side top htop t ht₀ htv hm
      exact ⟨o', List.mem_append_left _ ho', hl, by rw [hsd _ hl]; exact hsp⟩


/-- The returning back-edge push (no P-merge): the post-state has the new one-edge entry
`(v, lowval)` holding `Q e` on top of `L ++ s.tstack`, the items' children unchanged (`Q e` only got
its `vs`), and `nxtEdgeIdx` bumped. -/
theorem earCtx_pushBack {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s s' : WalkState} {e dest : Nat} {cls : OutClass}
    (hC : EarCtx v d done (.back e dest cls :: rest) hasVert base bE sv sd s)
    (hv : v < s.g.nv) (he : e < s.g.ne) (hlt : cls.lowval d < d)
    (hsz : 1 + s.g.nv + s.g.ne ≤ s.items.size)
    (hrank : ∀ o' ∈ done, o'.1.cls.rank ≤ cls.rank)
    (hinc : ∀ o' ∈ done, o'.1.e < s.g.ne ∧ s.g.Inc o'.1.e v)
    (hend : Items.PairEq (v, dest) s.g.edges[e]!) (hdest : s.stackVerts[cls.lowval d]! = dest)
    (hrest_e : ∀ o' ∈ rest, ∀ e', subEdges o' e' → e' ≠ e ∧ e' < s.g.ne)
    (hrest_v : ∀ o' ∈ rest, ∀ e₁ cls₁ child, o' = .tree e₁ cls₁ child →
      ∀ y ∈ child.verts, y < s.g.nv)
    (hnc : ∀ o' ∈ done, ∀ e', subEdges o'.1 e' → ∀ x, s.g.Inc e' x →
      ∀ o'' ∈ rest, ∀ e₁ cls₁ child, o'' = .tree e₁ cls₁ child → x ∉ child.verts)
    (L : List TEntry) (push dirD : Bool)
    (hpush : push = true ↔ hasVert = false)
    (hL : L = if push then [⟨v, d, s.nxtEdgeIdx, setSides dirD [vertItem v] []⟩] else [])
    (hg : s'.g = s.g) (hsv : s'.stackVerts = s.stackVerts)
    (hsd : ∀ k, k < d → s'.stackDir[k]! = s.stackDir[k]!)
    (hsz' : s'.items.size = s.items.size)
    (hch : ∀ j, Items.ch s'.items j = Items.ch s.items j)
    (hnx : s'.nxtEdgeIdx = s.nxtEdgeIdx + 1)
    (hts : s'.tstack = ⟨v, cls.lowval d, s.nxtEdgeIdx,
      setSides s.stackDir[cls.lowval d]! [edgeItem s.g e] []⟩ :: (L ++ s.tstack)) :
    EarCtx v d (done ++ [(.back e dest cls, true)]) rest true base bE sv sd s' := by
  set lv := cls.lowval d with hlv
  set q := edgeItem s.g e with hq
  set c : TEntry := ⟨v, lv, s.nxtEdgeIdx, setSides s.stackDir[lv]! [q] []⟩ with hc
  set o : DfsOut := .back e dest cls with ho
  have hqf := hC.q_fresh _ (List.mem_cons_self ..) e (Or.inl rfl)
  have hvq : vertItem v ≠ q := vertItem_ne_edgeItem' hv
  have hP : ∀ p c, Items.IsParent s'.items p c ↔ Items.IsParent s.items p c := by
    intro p c; simp only [Items.IsParent, hch]
  have hB : ∀ a i, Items.Below s'.items a i ↔ Items.Below s.items a i := by
    intro a i
    constructor
    · intro h
      induction h with
      | refl => exact .refl
      | tail _ hs ih => exact ih.tail ((hP _ _).1 hs)
    · intro h
      induction h with
      | refl => exact .refl
      | tail _ hs ih => exact ih.tail ((hP _ _).2 hs)
  have hEB : ∀ i e', Items.EdgeBelow s'.g s'.items i e' ↔ Items.EdgeBelow s.g s.items i e' := by
    intro i e'; rw [Items.EdgeBelow, Items.EdgeBelow, hg, hB]
  have hedges : ∀ t e', TEntry.edges s'.g s'.items t e' ↔ TEntry.edges s.g s.items t e' := by
    intro t e'; simp only [TEntry.edges, hEB]
  have hedgesF : ∀ t, TEntry.edges s'.g s'.items t = TEntry.edges s.g s.items t :=
    fun t => funext fun e' => propext (hedges t e')
  have hcs : c.spans.1 ++ c.spans.2 = [q] := by
    show (setSides s.stackDir[lv]! [q] []).1 ++ (setSides s.stackDir[lv]! [q] []).2 = [q]
    unfold setSides; split <;> rfl
  have hce : ∀ e', TEntry.edges s'.g s'.items c e' ↔ e' = e := by
    intro e'
    rw [hedges, TEntry.edges, hcs]
    constructor
    · rintro ⟨i, hi, hb⟩
      rw [List.mem_singleton] at hi; subst hi
      exact edgeItem_inj (below_of_ch_nil hqf.2.1 hb)
    · rintro rfl; exact ⟨q, List.mem_singleton_self _, .refl⟩
  have hinc_e : ∀ x, s.g.Inc e x → x = v ∨ x = dest := by
    intro x hx
    rcases hend with h | h
    · rw [Graph.Inc, ← h] at hx; exact hx.imp Eq.symm Eq.symm
    · obtain ⟨h1, h2⟩ := Prod.mk.inj h
      rcases hx with h' | h'
      · exact .inr (h'.symm.trans h2.symm)
      · exact .inl (h'.symm.trans h1.symm)
  have hIv : s.g.Inc e v := by
    rcases hend with h | h
    · exact .inl (by rw [← h])
    · exact .inr (Prod.mk.inj h).1.symm
  have hId : s.g.Inc e dest := by
    rcases hend with h | h
    · exact .inr (by rw [← h])
    · exact .inl (Prod.mk.inj h).2.symm
  have hTc : ∀ x, s'.g.Touches (TEntry.edges s'.g s'.items c) x → x = v ∨ x = dest := by
    rintro x ⟨e', -, hce', hx⟩
    rw [hce] at hce'; subst hce'
    rw [hg] at hx; exact hinc_e x hx
  have hTcv : s'.g.Touches (TEntry.edges s'.g s'.items c) v :=
    ⟨e, by rw [hg]; exact he, (hce e).2 rfl, by rw [hg]; exact hIv⟩
  have hTcd : s'.g.Touches (TEntry.edges s'.g s'.items c) dest :=
    ⟨e, by rw [hg]; exact he, (hce e).2 rfl, by rw [hg]; exact hId⟩
  have hvd : s.stackVerts[d]! = v := hC.sv_d
  have hLprop : ∀ t ∈ L, t.vStart = v ∧ t.topDepth = d ∧ t.firstIdx = s.nxtEdgeIdx ∧
      t.spans.1 ++ t.spans.2 = [vertItem v] ∧ hasVert = false := by
    intro t ht
    rw [hL] at ht
    split at ht
    · rw [List.mem_singleton] at ht; subst ht
      refine ⟨rfl, rfl, rfl, ?_, hpush.1 ‹_›⟩
      show (setSides dirD [vertItem v] []).1 ++ (setSides dirD [vertItem v] []).2 = _
      unfold setSides; split <;> rfl
    · exact absurd ht (List.not_mem_nil)
  have hLe : ∀ t ∈ L, ∀ e', TEntry.edges s'.g s'.items t e' ↔
      Items.EdgeBelow s.g s.items (vertItem v) e' := by
    intro t ht e'
    rw [hedges, TEntry.edges, (hLprop t ht).2.2.2.1]
    constructor
    · rintro ⟨i, hi, h⟩; rw [List.mem_singleton] at hi; subst hi; exact h
    · intro h; exact ⟨_, List.mem_singleton_self _, h⟩
  have hmem : ∀ t ∈ s'.tstack, t = c ∨ t ∈ L ∨ t ∈ s.tstack := by
    intro t ht
    rw [hts] at ht
    rcases List.mem_cons.1 ht with h | h
    · exact .inl h
    · exact .inr (List.mem_append.1 h)
  have hAV' : afterVert (done ++ [(o, true)]) = afterVert done ++ [o] := by simp [afterVert]
  have hoD : (o, true) ∈ done ++ [(o, true)] := List.mem_append_right _ (List.mem_singleton_self _)
  have hoA : o ∈ afterVert (done ++ [(o, true)]) := by
    rw [hAV']; exact List.mem_append_right _ (List.mem_singleton_self _)
  have hsub : ∀ e', subEdges o e' → e' = e := by
    intro e' h
    rcases h with h | h
    · exact h
    · exact h.elim
  have hVL : ∀ x, x ∈ DfsOut.vertsList (done.map (·.1)) →
      x ∈ DfsOut.vertsList ((done ++ [(o, true)]).map (·.1)) := by
    intro x hx
    rw [DfsOut.vertsList_eq] at hx ⊢
    rw [List.map_append, List.flatMap_append]
    exact List.mem_append_left _ hx
  have hvne : ∀ k, k < d → v ≠ s.stackVerts[k]! := fun k hk h =>
    hC.path k d hk (Nat.le_refl _) (by rw [hvd, h])
  have hqlt : q < s'.items.size := by rw [hsz']; show 1 + s.g.nv + e < _; omega
  have hvlt : vertItem v < s'.items.size := by rw [hsz']; show 1 + v < _; omega
  obtain ⟨top, htop, hCT⟩ := hC.top
  have hmemT : ∀ t ∈ top, t ∈ s.tstack := fun t ht => by
    rw [htop]; exact List.mem_append_left _ ht
  have hbase_mem : ∀ k, k < base.length → base[k]! ∈ s.tstack := fun k hk => by
    rw [htop, getElem!_pos base k hk]; exact List.mem_append_right _ (List.getElem_mem hk)
  obtain ⟨above, below, hsp, hCE, hCS, hPW, hLow, hEdg, hBel⟩ := hCT.split
  have habm : ∀ t ∈ above, t ∈ s.tstack := by
    intro t ht
    apply hmemT
    cases hasVert with
    | true =>
      obtain ⟨vt, htv, -⟩ := hsp
      rw [htv]; exact List.mem_append_left _ ht
    | false => exact absurd ht (hsp.1 ▸ List.not_mem_nil)
  have habv : ∀ t ∈ above, vertItem v ∉ t.spans.1 ++ t.spans.2 := by
    intro t ht
    cases hasVert with
    | true =>
      obtain ⟨vt, -, -, -, -, h4, -⟩ := hsp
      exact h4 t ht
    | false => exact absurd ht (hsp.1 ▸ List.not_mem_nil)
  have hCE' : ∀ t, CtxEntry v d s t → CtxEntry v d s' t := by
    intro t h
    obtain ⟨e', h1, h2⟩ := h.nonempty
    exact
      { vStart := h.vStart
        depth := h.depth
        nonempty := ⟨e', by rw [hg]; exact h1, (hedges t e').2 h2⟩
        touch_bot := by rw [hedgesF, hg]; exact h.touch_bot
        touch_top := by rw [hedgesF, hg, hsv]; exact h.touch_top
        side := by rw [hsd _ h.depth]; exact h.side
        att := by rw [hedgesF, hg, hsv]; exact h.att }
  have hCS' : ∀ t, t.topDepth < d → CtxSingle v s t → CtxSingle v s' t := by
    intro t hd h
    obtain ⟨i, hi, hroot⟩ := h.item
    exact
      { item := ⟨i, by rw [hsd _ hd]; exact hi, fun p hp => hroot p ((hP _ _).1 hp)⟩
        att := by rw [hedgesF, hg, hsv]; exact h.att }
  have hcCE : CtxEntry v d s' c :=
    { vStart := rfl
      depth := hlt
      nonempty := ⟨e, by rw [hg]; exact he, (hce e).2 rfl⟩
      touch_bot := hTcv
      touch_top := by show s'.g.Touches _ s'.stackVerts[lv]!; rw [hsv, hdest]; exact hTcd
      side := by
        show getSide (setSides s.stackDir[lv]! [q] []) (!s'.stackDir[lv]!) = []
        rw [hsd lv hlt]; unfold setSides getSide; cases s.stackDir[lv]! <;> rfl
      att := fun x hx hxv _ => by
        rcases hTc x hx with h | h
        · exact absurd h hxv
        · exact ⟨lv, Nat.le_refl _, hlt.le, by rw [hsv, hdest]; exact h⟩ }
  have hcCS : CtxSingle v s' c :=
    { item := ⟨q, by rw [hsd lv hlt], fun p hp => hqf.1 p ((hP _ _).1 hp)⟩
      att := fun x hx => by
        rcases hTc x hx with h | h
        · exact .inl h
        · exact .inr (.inl (by show x = s'.stackVerts[lv]!; rw [hsv, hdest]; exact h)) }
  have hfalse : ∀ t ∈ L, hasVert = false := fun t ht => (hLprop t ht).2.2.2.2
  have hAVnil : hasVert = false → afterVert done = [] := hC.noVert_after
  have hLspan : ∀ t ∈ L, ∀ i, i ∈ t.spans.1 ++ t.spans.2 → i = vertItem v := by
    intro t ht i hi
    rw [(hLprop t ht).2.2.2.1] at hi
    exact List.mem_singleton.1 hi
  have hcspan : ∀ i, i ∈ c.spans.1 ++ c.spans.2 → i = q := by
    intro i hi; rw [hcs] at hi; exact List.mem_singleton.1 hi
  have hTL : ∀ t ∈ L, ∀ x, s'.g.Touches (TEntry.edges s'.g s'.items t) x →
      ∃ o' ∈ done, ∃ e', e' < s.g.ne ∧ subEdges o'.1 e' ∧ s.g.Inc e' x := by
    rintro t ht x ⟨e', he', hte, hx⟩
    rw [hg] at he' hx
    obtain ⟨o', ho', -, hs⟩ := (hC.vert_edges e' he').1 ((hLe t ht e').1 hte)
    exact ⟨o', ho', e', he', hs, hx⟩
  refine
    { top := ⟨c :: L ++ top, by rw [hts, htop]; simp,
        { split := ⟨c :: above, below, ?_, ?_, ?_, ?_, ?_, ?_, hBel⟩
          bot := fun t ht k hk => by
            rw [hsv]
            rcases List.mem_cons.1 ht with rfl | ht
            · exact hvne k hk
            rcases List.mem_append.1 ht with ht | ht
            · rw [(hLprop t ht).1]; exact hvne k hk
            · exact hCT.bot t ht k hk
          vitems := fun t ht x hx hm => by
            rw [hg] at hx
            rcases List.mem_cons.1 ht with rfl | ht
            · exact absurd (hcspan _ hm) (vertItem_ne_edgeItem' hx)
            rcases List.mem_append.1 ht with ht | ht
            · exact .inl (vertItem_inj' (hLspan t ht _ hm))
            · exact (hCT.vitems t ht x hx hm).imp id (hVL x)
          qitems := fun t ht e' he' hm => by
            rw [hg] at he' hm
            rcases List.mem_cons.1 ht with rfl | ht
            · exact ⟨(o, true), hoD, Or.inl (edgeItem_inj (hcspan _ hm))⟩
            rcases List.mem_append.1 ht with ht | ht
            · exact absurd (hLspan t ht _ hm) (vertItem_ne_edgeItem' hv).symm
            · obtain ⟨o', ho', hs⟩ := hCT.qitems t ht e' he' hm
              exact ⟨o', List.mem_append_left _ ho', hs⟩
          edges := fun t ht e' he' hte => by
            rw [hg] at he'
            rcases List.mem_cons.1 ht with rfl | ht
            · exact ⟨(o, true), hoD, Or.inl ((hce e').1 hte)⟩
            rcases List.mem_append.1 ht with ht | ht
            · obtain ⟨o', ho', -, hs⟩ := (hC.vert_edges e' he').1 ((hLe t ht e').1 hte)
              exact ⟨o', List.mem_append_left _ ho', hs⟩
            · obtain ⟨o', ho', hs⟩ := hCT.edges t ht e' he' ((hedges t e').1 hte)
              exact ⟨o', List.mem_append_left _ ho', hs⟩
          cover := fun o' ho' hl e' he' hs => by
            rw [hg] at he'
            rcases List.mem_append.1 ho' with h | h
            · obtain ⟨t, ht, hte⟩ := hCT.cover o' h hl e' he' hs
              exact ⟨t, List.mem_cons_of_mem _ (List.mem_append_right _ ht), (hedges t e').2 hte⟩
            · rw [List.mem_singleton] at h
              subst h
              exact ⟨c, List.mem_cons_self .., (hce e').2 (hsub e' hs)⟩
          ret := fun _ => ⟨(.back e dest cls, true), List.mem_append_right _ (List.mem_singleton_self _), hlt⟩ }⟩
      base_edges := ⟨hC.base_edges.1, fun k hk e' he' => by
        rw [hg] at he'
        exact (hedges _ e').trans (hC.base_edges.2 k hk e' he')⟩
      base_bot := hC.base_bot
      noVert_after := fun h => nomatch h
      vfirst := fun t ht hvt hnv => by
        rw [hnx]
        rcases hmem t ht with rfl | ht | ht
        · exact Nat.lt_succ_self _
        · exact absurd (by rw [(hLprop t ht).2.2.2.1]; exact List.mem_singleton_self _) hnv
        · exact Nat.lt_succ_of_lt (hC.vfirst t ht hvt hnv)
      sv := fun k hk => by rw [hsv]; exact hC.sv k hk
      sd := fun k hk => by rw [hsd k hk]; exact hC.sd k hk
      sv_d := by rw [hsv]; exact hvd
      path := fun k k' h h' => by rw [hsv]; exact hC.path k k' h h'
      v_root := fun p h => hC.v_root p ((hP _ _).1 h)
      vert_free := fun _ _ _ => rfl
      afterVert_ret := fun o' ho' => by
        rw [hAV'] at ho'
        rcases List.mem_append.1 ho' with h | h
        · exact hC.afterVert_ret o' h
        · rw [List.mem_singleton] at h; subst h; exact hlt
      hv_ret := fun _ => ⟨(o, true), hoD, hlt⟩
      vert_book := fun h => nomatch h
      vert_disj := fun h => nomatch h
      vert_edges := fun e' he' => by
        rw [hg] at he'
        rw [hEB, hC.vert_edges e' he']
        constructor
        · rintro ⟨o', ho', hl, hs⟩; exact ⟨o', List.mem_append_left _ ho', hl, hs⟩
        · rintro ⟨o', ho', hl, hs⟩
          rcases List.mem_append.1 ho' with h | h
          · exact ⟨o', h, hl, hs⟩
          · rw [List.mem_singleton] at h; subst h; exact absurd hl (Nat.not_le.2 hlt)
      touch_bot := fun t ht hne => by
        rcases hmem t ht with rfl | ht | ht
        · exact hTcv
        · obtain ⟨e', he', hte⟩ := hne
          rw [hg] at he'
          obtain ⟨o', ho', hlo, -⟩ := (hC.vert_edges e' he').1 ((hLe t ht e').1 hte)
          rw [(hLprop t ht).1]
          exact ⟨o'.1.e, by rw [hg]; exact (hinc o' ho').1,
            (hLe t ht _).2 ((hC.vert_edges _ (hinc o' ho').1).2 ⟨o', ho', hlo, subEdges_e _⟩),
            by rw [hg]; exact (hinc o' ho').2⟩
        · rw [hedgesF, hg]
          rw [hedgesF, hg] at hne
          exact hC.touch_bot t ht hne
      span_root := fun t ht i hi p hp => by
        rw [hP] at hp
        rcases hmem t ht with rfl | ht | ht
        · rw [hcspan i hi] at hp; exact hqf.1 p hp
        · rw [hLspan t ht i hi] at hp; exact hC.v_root p hp
        · exact hC.span_root t ht i hi p hp
      span_lt := fun t ht i hi => by
        rcases hmem t ht with rfl | ht | ht
        · rw [hcspan i hi]; exact hqlt
        · rw [hLspan t ht i hi]; exact hvlt
        · rw [hsz']; exact hC.span_lt t ht i hi
      ch_lt := fun p c h => by rw [hsz']; exact hC.ch_lt p c ((hP _ _).1 h)
      disj := by
        rw [hts]
        refine List.pairwise_cons.2 ⟨fun t' ht' e' he' h1 h2 => ?_, ?_⟩
        · rw [hce] at h1; subst h1
          rw [hedges] at h2
          obtain ⟨i, hi, hb⟩ := h2
          have := Items.Below.eq_of_no_parent hqf.1 hb
          subst this
          rcases List.mem_append.1 ht' with ht' | ht'
          · exact hvq (hLspan t' ht' _ hi).symm
          · exact hqf.2.2 t' ht' hi
        · refine List.pairwise_append.2 ⟨by rw [hL]; split <;> simp, ?_,
            fun t ht t' ht' e' he' h1 h2 => ?_⟩
          · refine hC.disj.imp fun {t t'} h e' he' h1 h2 => ?_
            rw [hg] at he'
            exact h e' he' ((hedges t e').1 h1) ((hedges t' e').1 h2)
          · rw [hg] at he'
            exact hC.vert_disj (hfalse t ht) t' ht' e' he' ((hedges t' e').1 h2)
              ((hLe t ht e').1 h1)
      span_disj := by
        rw [hts]
        refine List.pairwise_cons.2 ⟨fun t' ht' i hi hi' => ?_, ?_⟩
        · rw [hcspan i hi] at hi'
          rcases List.mem_append.1 ht' with ht' | ht'
          · exact hvq (hLspan t' ht' _ hi').symm
          · exact hqf.2.2 t' ht' hi'
        · refine List.pairwise_append.2 ⟨by rw [hL]; split <;> simp, hC.span_disj,
            fun t ht t' ht' i hi hi' => ?_⟩
          · rw [hLspan t ht i hi] at hi'
            exact Bool.false_ne_true ((hC.vert_free t' ht' hi').symm.trans (hfalse t ht)).symm
      q_fresh := fun o' ho' e' hs => by
        obtain ⟨h1, h2, h3⟩ := hC.q_fresh o' (List.mem_cons_of_mem _ ho') e' hs
        rw [hg]
        refine ⟨fun p hp => h1 p ((hP _ _).1 hp), by rw [hch]; exact h2, fun t ht hm => ?_⟩
        rcases hmem t ht with rfl | ht | ht
        · exact (hrest_e o' ho' e' hs).1 (edgeItem_inj (hcspan _ hm))
        · exact (vertItem_ne_edgeItem' hv) (hLspan t ht _ hm).symm
        · exact h3 t ht hm
      v_fresh := fun o' ho' e₁ cls₁ child h y hy => by
        obtain ⟨h1, h2, h3, h4, h5, h6⟩ := hC.v_fresh o' (List.mem_cons_of_mem _ ho') e₁ cls₁ child h y hy
        have hyv : y < s.g.nv := hrest_v o' ho' e₁ cls₁ child h y hy
        have hyne : v ≠ y := fun hvy => h4 d (Nat.le_refl _) (hvd.trans hvy)
        refine ⟨fun p hp => h1 p ((hP _ _).1 hp), by rw [hch]; exact h2, fun t ht hm => ?_,
          fun k hk => by rw [hsv]; exact h4 k hk, fun t ht => ?_, fun t ht hT => ?_⟩
        · rcases hmem t ht with rfl | ht | ht
          · exact (vertItem_ne_edgeItem' hyv) (hcspan _ hm)
          · exact hyne (vertItem_inj' (hLspan t ht _ hm)).symm
          · exact h3 t ht hm
        · rcases hmem t ht with rfl | ht | ht
          · exact hyne
          · rw [(hLprop t ht).1]; exact hyne
          · exact h5 t ht
        · rcases hmem t ht with rfl | ht | ht
          · rcases hTc y hT with hh | hh
            · exact hyne hh.symm
            · exact h4 lv hlt.le (hdest.trans hh.symm)
          · obtain ⟨o'', ho'', e', he', hs, hx⟩ := hTL t ht y hT
            exact hnc o'' ho'' e' hs y hx o' ho' e₁ cls₁ child h hy
          · rw [hedgesF, hg] at hT
            exact h6 t ht hT }
  · rw [if_pos rfl]
    by_cases hhv : hasVert = true
    · rw [if_pos hhv] at hsp
      obtain ⟨vt, htv, h1, h2, h3, h4, h5, h6⟩ := hsp
      have hpf : push = false := by
        cases hp : push
        · rfl
        · exact absurd (hpush.1 hp) (by rw [hhv]; decide)
      have hL0 : L = [] := by rw [hL, hpf]; rfl
      refine ⟨vt, by rw [hL0, htv]; rfl, h1, h2, h3, fun t ht => ?_, fun o' ho' => ?_,
        fun o' ho' e' he' hs => ?_⟩
      · rcases List.mem_cons.1 ht with rfl | ht
        · intro hm; exact hvq (hcspan _ hm)
        · exact h4 t ht
      · rw [hAV'] at ho'
        rcases List.mem_append.1 ho' with ho' | ho'
        · exact (h5 o' ho').imp (fun ⟨t, ht, h⟩ => ⟨t, List.mem_cons_of_mem _ ht, h⟩) id
        · rw [List.mem_singleton] at ho'; subst ho'; exact .inl ⟨c, List.mem_cons_self .., rfl⟩
      · rw [hg] at he'
        rw [hAV'] at ho'
        rcases List.mem_append.1 ho' with ho' | ho'
        · exact (h6 o' ho' e' he' hs).imp
            (fun ⟨t, ht, h⟩ => ⟨t, List.mem_cons_of_mem _ ht, (hedges t e').2 h⟩) (hedges vt e').2
        · rw [List.mem_singleton] at ho'; subst ho'
          exact .inl ⟨c, List.mem_cons_self .., (hce e').2 (hsub e' hs)⟩
    · rw [if_neg hhv] at hsp
      have hhf : hasVert = false := by
        cases h : hasVert with
        | false => rfl
        | true => exact absurd h hhv
      obtain ⟨hab, htb⟩ := hsp
      have hpt : push = true := hpush.2 hhf
      have hL1 : L = [⟨v, d, s.nxtEdgeIdx, setSides dirD [vertItem v] []⟩] := by rw [hL, hpt]; rfl
      have hAV0 := hAVnil hhf
      refine ⟨⟨v, d, s.nxtEdgeIdx, setSides dirD [vertItem v] []⟩, by rw [hL1, hab, htb]; rfl,
        Nat.le_refl _, ?_, fun _ => ⟨rfl, ?_⟩, fun t ht => ?_, fun o' ho' => ?_,
        fun o' ho' e' he' hs => ?_⟩
      · show vertItem v ∈ (setSides dirD [vertItem v] []).1 ++ (setSides dirD [vertItem v] []).2
        unfold setSides; cases dirD <;> simp
      · show (setSides dirD [vertItem v] []).1 ++ (setSides dirD [vertItem v] []).2 = _
        unfold setSides; cases dirD <;> rfl
      · rcases List.mem_cons.1 ht with rfl | ht
        · intro hm; exact hvq (hcspan _ hm)
        · exact absurd ht (hab ▸ List.not_mem_nil)
      · rw [hAV', hAV0, List.nil_append, List.mem_singleton] at ho'; subst ho'
        exact .inl ⟨c, List.mem_cons_self .., rfl⟩
      · rw [hAV', hAV0, List.nil_append, List.mem_singleton] at ho'; subst ho'
        exact .inl ⟨c, List.mem_cons_self .., (hce e').2 (hsub e' hs)⟩
  · intro t ht
    rcases List.mem_cons.1 ht with rfl | ht
    · exact hcCE
    · exact hCE' t (hCE t ht)
  · intro t ht hA
    rcases List.mem_cons.1 ht with rfl | ht
    · exact hcCS
    · exact hCS' t (hCE t ht).depth
        (hCS t ht fun o' ho' hl => hA o' (by rw [hAV']; exact List.mem_append_left _ ho') hl)
  · refine List.pairwise_cons.2 ⟨fun t' ht' => ⟨?_, ?_⟩, hPW⟩
    · obtain ⟨o', ho', hto⟩ := hLow t' ht'
      show t'.topDepth ≤ lv
      rw [hto]
      exact lowval_le_of_rank (hC.afterVert_ret o' ho') hlt (hrank _ (mem_afterVert ho'))
    · show t'.firstIdx < s.nxtEdgeIdx
      exact hC.vfirst t' (habm t' ht') (hCE t' ht').vStart (habv t' ht')
  · intro t ht
    rcases List.mem_cons.1 ht with rfl | ht
    · exact ⟨o, hoA, rfl⟩
    · obtain ⟨o', ho', h⟩ := hLow t ht
      exact ⟨o', by rw [hAV']; exact List.mem_append_left _ ho', h⟩
  · intro t ht e' he' hte
    rw [hg] at he'
    rcases List.mem_cons.1 ht with rfl | ht
    · exact ⟨o, hoA, Or.inl ((hce e').1 hte)⟩
    · obtain ⟨o', ho', h⟩ := hEdg t ht e' he' ((hedges t e').1 hte)
      exact ⟨o', by rw [hAV']; exact List.mem_append_left _ ho', h⟩

theorem ch_of_ge {items : Items} {p : ItemId} (h : items.size ≤ p) : Items.ch items p = [] := by
  simp [Items.ch, Array.getElem?_eq_none_iff.2 h]

theorem below_addChildren {I I' : Items} {p : ItemId} {S : List ItemId}
    (hP : ∀ p' c', Items.IsParent I' p' c' ↔ Items.IsParent I p' c' ∨ (p' = p ∧ c' ∈ S))
    (x i : ItemId) :
    Items.Below I' x i ↔ Items.Below I x i ∨ (Items.Below I x p ∧ ∃ r ∈ S, Items.Below I r i) := by
  have hm : ∀ {a b}, Items.Below I a b → Items.Below I' a b :=
    fun h => Relation.ReflTransGen.mono (fun a b hab => (hP a b).2 (.inl hab)) _ _ h
  constructor
  · intro h
    induction h with
    | refl => exact .inl .refl
    | @tail b c _ hs ih =>
      rcases (hP _ _).1 hs with hs | ⟨hb, hc⟩
      · rcases ih with ih | ⟨h1, r, hr, h2⟩
        · exact .inl (ih.tail hs)
        · exact .inr ⟨h1, r, hr, h2.tail hs⟩
      · subst hb
        rcases ih with ih | ⟨h1, -⟩
        · exact .inr ⟨ih, c, hc, .refl⟩
        · exact .inr ⟨h1, c, hc, .refl⟩
  · rintro (h | ⟨h1, r, hr, h2⟩)
    · exact hm h
    · exact (hm h1).trans ((Relation.ReflTransGen.single ((hP _ _).2 (.inr ⟨rfl, hr⟩))).trans (hm h2))


/-- When the P-check fires after the `Q e` push, the entry below the new top is a type-1 ear
`(v, lowval)` whose active side is a single item with no parent (`CtxTop`/`vfirst`/`allType1`). -/
theorem mergeP_single {v d lv : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s : WalkState} {c a : TEntry} {tl : List TEntry}
    (hC : EarCtx v d done rest true base bE sv sd s)
    (hts : s.tstack = c :: a :: tl) (hcV : vertItem v ∉ c.spans.1 ++ c.spans.2)
    (hav : a.vStart = v) (had : a.topDepth = lv) (hlt : lv < d)
    (hA1 : allType1 d lv done) :
    ∃ i, a.spans = setSides s.stackDir[lv]! [i] [] ∧ ∀ p, ¬ Items.IsParent s.items p i := by
  obtain ⟨top, htop, hCT⟩ := hC.top
  obtain ⟨above, below, hsp, hCE, hCS, -, -, -, -⟩ := hCT.split
  rw [if_pos rfl] at hsp
  obtain ⟨vt, htv, -, hvt2, hvt3, -, -, -⟩ := hsp
  have ha : a ∈ above := by
    have h := htop
    rw [hts, htv] at h
    rcases above with _ | ⟨x, _ | ⟨y, r⟩⟩
    · simp only [List.nil_append, List.cons_append, List.cons.injEq] at h
      rw [← h.1] at hvt2; exact absurd hvt2 hcV
    · simp only [List.cons_append, List.nil_append, List.cons.injEq] at h
      obtain ⟨-, h2, -⟩ := h
      rw [← h2] at hvt3
      exact absurd (hvt3 hav).1 (by rw [had]; exact Nat.ne_of_lt hlt)
    · simp only [List.cons_append, List.cons.injEq] at h
      obtain ⟨-, h2, -⟩ := h
      rw [← h2]; exact List.mem_cons_of_mem _ (List.mem_cons_self ..)
  have hA1a : allType1 d a.topDepth done := by rw [had]; exact hA1
  obtain ⟨i, hi, hroot⟩ := (hCS a ha hA1a).item
  rw [had] at hi
  exact ⟨i, hi, hroot⟩

/-- Abstract post-state of the P-merge: the top `Q e` and the `(v, lowval)` entry below it are
replaced by a single entry whose side is `[p]`, where `p` (fresh or the reused `P` item) gains the
children `S` (`Q e` plus the old side / its children). All `EarCtx` clauses are transported. -/
theorem earCtx_mergeP {v d lv : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s s' : WalkState} {c a : TEntry} {tl : List TEntry} {p q i : ItemId}
    {S : List ItemId}
    (hC : EarCtx v d done rest true base bE sv sd s)
    (hv : v < s.g.nv) (hlt : lv < d)
    (hsz : 1 + s.g.nv + s.g.ne ≤ s.items.size)
    (hts : s.tstack = c :: a :: tl)
    (hcd : c.topDepth = lv) (hcs : c.spans.1 ++ c.spans.2 = [q]) (hcv : c.vStart = v)
    (hav : a.vStart = v) (had : a.topDepth = lv)
    (hai : a.spans = setSides s.stackDir[lv]! [i] [])
    (hA1 : allType1 d lv done)
    (hrest_e : ∀ o' ∈ rest, ∀ e', subEdges o' e' → e' < s.g.ne)
    (hrest_v : ∀ o' ∈ rest, ∀ e₁ cls₁ child, o' = .tree e₁ cls₁ child →
      ∀ y ∈ child.verts, y < s.g.nv)
    (hpa : p = i ∨ s.items.size ≤ p)
    (hS : ∀ r ∈ S, r = q ∨ (r = i ∧ i ≠ p) ∨ Items.IsParent s.items p r)
    (hSq : q ∈ S) (hSi : i ≠ p → i ∈ S)
    (hP : ∀ p' c', Items.IsParent s'.items p' c' ↔
      Items.IsParent s.items p' c' ∨ (p' = p ∧ c' ∈ S))
    (hch : ∀ j, j ≠ p → Items.ch s'.items j = Items.ch s.items j)
    (hsz' : s.items.size ≤ s'.items.size) (hplt : p < s'.items.size)
    (hg : s'.g = s.g) (hsv : s'.stackVerts = s.stackVerts) (hsd : s'.stackDir = s.stackDir)
    (hnx : s'.nxtEdgeIdx = s.nxtEdgeIdx)
    (hts' : s'.tstack = ⟨v, lv, a.firstIdx, setSides s.stackDir[lv]! [p] []⟩ :: tl) :
    EarCtx v d done rest true base bE sv sd s' := by
  set a₃ : TEntry := ⟨v, lv, a.firstIdx, setSides s.stackDir[lv]! [p] []⟩ with ha₃
  have hcm : c ∈ s.tstack := by rw [hts]; exact List.mem_cons_self ..
  have ham : a ∈ s.tstack := by rw [hts]; exact List.mem_cons_of_mem _ (List.mem_cons_self ..)
  have htlm : ∀ t ∈ tl, t ∈ s.tstack := fun t ht => by
    rw [hts]; exact List.mem_cons_of_mem _ (List.mem_cons_of_mem _ ht)
  have has : a.spans.1 ++ a.spans.2 = [i] := by
    rw [hai]; unfold setSides; split <;> rfl
  have ha₃s : a₃.spans.1 ++ a₃.spans.2 = [p] := by
    show (setSides s.stackDir[lv]! [p] []).1 ++ (setSides s.stackDir[lv]! [p] []).2 = [p]
    unfold setSides; split <;> rfl
  have hcspan : ∀ j, j ∈ c.spans.1 ++ c.spans.2 → j = q := fun j hj => by
    rw [hcs] at hj; exact List.mem_singleton.1 hj
  have hqm : q ∈ c.spans.1 ++ c.spans.2 := by rw [hcs]; exact List.mem_singleton_self _
  have him : i ∈ a.spans.1 ++ a.spans.2 := by rw [has]; exact List.mem_singleton_self _
  have hqlt : q < s.items.size := hC.span_lt c hcm q hqm
  have hilt : i < s.items.size := hC.span_lt a ham i him
  have hvlt : vertItem v < s.items.size := by show 1 + v < _; omega
  have hdisj := hC.disj
  rw [hts] at hdisj
  have hsdisj := hC.span_disj
  rw [hts] at hsdisj
  have hca_sd := (List.pairwise_cons.1 hsdisj).1 a (List.mem_cons_self ..)
  have hct_sd : ∀ t ∈ tl, ∀ j ∈ c.spans.1 ++ c.spans.2, j ∉ t.spans.1 ++ t.spans.2 :=
    fun t ht => (List.pairwise_cons.1 hsdisj).1 t (List.mem_cons_of_mem _ ht)
  have hat_sd : ∀ t ∈ tl, ∀ j ∈ a.spans.1 ++ a.spans.2, j ∉ t.spans.1 ++ t.spans.2 :=
    fun t ht => (List.pairwise_cons.1 (List.pairwise_cons.1 hsdisj).2).1 t ht
  have htl_sd := (List.pairwise_cons.1 (List.pairwise_cons.1 hsdisj).2).2
  have hct_d : ∀ t ∈ tl, ∀ e', e' < s.g.ne → TEntry.edges s.g s.items c e' →
      ¬ TEntry.edges s.g s.items t e' :=
    fun t ht => (List.pairwise_cons.1 hdisj).1 t (List.mem_cons_of_mem _ ht)
  have hat_d : ∀ t ∈ tl, ∀ e', e' < s.g.ne → TEntry.edges s.g s.items a e' →
      ¬ TEntry.edges s.g s.items t e' :=
    fun t ht => (List.pairwise_cons.1 (List.pairwise_cons.1 hdisj).2).1 t ht
  have htl_d := (List.pairwise_cons.1 (List.pairwise_cons.1 hdisj).2).2
  have hqi : q ≠ i := fun h => hca_sd q hqm (h ▸ him)
  have hiroot : ∀ p', ¬ Items.IsParent s.items p' i := hC.span_root a ham i him
  have hproot : ∀ p', ¬ Items.IsParent s.items p' p := by
    rcases hpa with hpi | hp
    · rw [hpi]; exact hiroot
    · intro p' h; exact absurd (hC.ch_lt p' p h) (Nat.not_lt.2 hp)
  have hpS : p ∉ S := fun hp => by
    rcases hS p hp with h | ⟨h, hne⟩ | h
    · rcases hpa with hpi | hp'
      · exact hqi (h.symm.trans hpi)
      · exact absurd hqlt (by rw [← h]; exact Nat.not_lt.2 hp')
    · exact hne h.symm
    · exact hproot p h
  have hpT : ∀ t ∈ tl, p ∉ t.spans.1 ++ t.spans.2 := fun t ht hp => by
    rcases hpa with hpi | hp'
    · rw [hpi] at hp; exact hat_sd t ht i him hp
    · exact absurd (hC.span_lt t (htlm t ht) p hp) (Nat.not_lt.2 hp')
  obtain ⟨top, htop, hCT⟩ := hC.top
  obtain ⟨above, below, hsp, hCE, hCS, hPW, hLow, hEdg, hBel⟩ := hCT.split
  rw [if_pos rfl] at hsp
  obtain ⟨vt, htv, hvt1, hvt2, hvt3, hvt4, hvt5, hvt6⟩ := hsp
  obtain ⟨above'', hab⟩ : ∃ above'', above = c :: a :: above'' := by
    have h := htop
    rw [hts, htv] at h
    rcases above with _ | ⟨x, _ | ⟨y, r⟩⟩
    · simp only [List.nil_append, List.cons_append, List.cons.injEq] at h
      rw [← h.1] at hvt3; exact absurd (hvt3 hcv).1 (by rw [hcd]; exact Nat.ne_of_lt hlt)
    · simp only [List.cons_append, List.nil_append, List.cons.injEq] at h
      obtain ⟨-, h2, -⟩ := h
      rw [← h2] at hvt3
      exact absurd (hvt3 hav).1 (by rw [had]; exact Nat.ne_of_lt hlt)
    · simp only [List.cons_append, List.cons.injEq] at h
      obtain ⟨h1, h2, -⟩ := h
      exact ⟨r, by rw [h1, h2]⟩
  have htl : tl = (above'' ++ vt :: below) ++ base := by
    have h := htop
    rw [hts, htv, hab] at h
    simp only [List.cons_append, List.cons.injEq] at h
    exact h.2.2
  have hmemT' : ∀ t ∈ above'' ++ vt :: below, t ∈ tl := fun t ht => by
    rw [htl]; exact List.mem_append_left _ ht
  have hab_tl : ∀ t ∈ above'', t ∈ tl := fun t ht => hmemT' t (List.mem_append_left _ ht)
  have hvt_tl : vt ∈ tl := hmemT' vt (List.mem_append_right _ (List.mem_cons_self ..))
  have hbase_tl : ∀ k, k < base.length → base[k]! ∈ tl := fun k hk => by
    rw [htl, getElem!_pos base k hk]; exact List.mem_append_right _ (List.getElem_mem hk)
  have hcmA : c ∈ above := by rw [hab]; exact List.mem_cons_self ..
  have hamA : a ∈ above := by rw [hab]; exact List.mem_cons_of_mem _ (List.mem_cons_self ..)
  have habA : ∀ t ∈ above'', t ∈ above := fun t ht => by
    rw [hab]; exact List.mem_cons_of_mem _ (List.mem_cons_of_mem _ ht)
  have hmemTop : ∀ t ∈ above'' ++ vt :: below, t ∈ top := fun t ht => by
    rw [htv]
    rcases List.mem_append.1 ht with ht | ht
    · exact List.mem_append_left _ (habA t ht)
    · exact List.mem_append_right _ ht
  have hcmT : c ∈ top := by rw [htv]; exact List.mem_append_left _ hcmA
  have hamT : a ∈ top := by rw [htv]; exact List.mem_append_left _ hamA
  have hvq : vertItem v ≠ q := fun h => hvt4 c hcmA (h ▸ hqm)
  have hva : vertItem v ∉ a.spans.1 ++ a.spans.2 := hvt4 a hamA
  have hpv : p ≠ vertItem v := by
    rcases hpa with hpi | hp
    · exact fun h => hva (by rw [← h, hpi]; exact him)
    · exact fun h => absurd hvlt (by rw [← h]; exact Nat.not_lt.2 hp)
  have hB := below_addChildren hP
  have hBp : ∀ x, Items.Below s.items x p → x = p := fun x h =>
    Items.Below.eq_of_no_parent hproot h
  have hm : ∀ {x j}, Items.Below s.items x j → Items.Below s'.items x j :=
    fun h => (hB _ _).2 (.inl h)
  have hpch : s.items.size ≤ p → ∀ j, Items.Below s.items p j → j = p :=
    fun hp j h => below_of_ch_nil (ch_of_ge hp) h
  have hBpj : ∀ j, Items.Below s'.items p j ↔
      Items.Below s.items p j ∨ ∃ r ∈ S, Items.Below s.items r j := by
    intro j; rw [hB]
    constructor
    · rintro (h | ⟨-, h⟩)
      · exact .inl h
      · exact .inr h
    · rintro (h | h)
      · exact .inl h
      · exact .inr ⟨.refl, h⟩
  have hBx : ∀ x j, x ≠ p → (Items.Below s'.items x j ↔ Items.Below s.items x j) := by
    intro x j hx; rw [hB]
    constructor
    · rintro (h | ⟨h, -⟩)
      · exact h
      · exact absurd (hBp x h) hx
    · exact .inl
  have hEBv : ∀ e', Items.EdgeBelow s'.g s'.items (vertItem v) e' ↔
      Items.EdgeBelow s.g s.items (vertItem v) e' := by
    intro e'; rw [Items.EdgeBelow, Items.EdgeBelow, hg]; exact hBx _ _ hpv.symm
  have hedgesT : ∀ t ∈ tl, ∀ e', TEntry.edges s'.g s'.items t e' ↔
      TEntry.edges s.g s.items t e' := by
    intro t ht e'
    simp only [TEntry.edges, Items.EdgeBelow, hg]
    constructor
    · rintro ⟨j, hj, h⟩; exact ⟨j, hj, (hBx j _ (fun h' => hpT t ht (h' ▸ hj))).1 h⟩
    · rintro ⟨j, hj, h⟩; exact ⟨j, hj, hm h⟩
  have hedgesTF : ∀ t ∈ tl, TEntry.edges s'.g s'.items t = TEntry.edges s.g s.items t :=
    fun t ht => funext fun e' => propext (hedgesT t ht e')
  have hce : ∀ e', TEntry.edges s.g s.items c e' ↔ Items.Below s.items q (edgeItem s.g e') := by
    intro e'; rw [TEntry.edges, hcs]
    constructor
    · rintro ⟨j, hj, h⟩; rw [List.mem_singleton] at hj; subst hj; exact h
    · intro h; exact ⟨q, List.mem_singleton_self _, h⟩
  have hae : ∀ e', TEntry.edges s.g s.items a e' ↔ Items.Below s.items i (edgeItem s.g e') := by
    intro e'; rw [TEntry.edges, has]
    constructor
    · rintro ⟨j, hj, h⟩; rw [List.mem_singleton] at hj; subst hj; exact h
    · intro h; exact ⟨i, List.mem_singleton_self _, h⟩
  have hce' : ∀ e', e' < s.g.ne → (TEntry.edges s'.g s'.items a₃ e' ↔
      TEntry.edges s.g s.items a e' ∨ TEntry.edges s.g s.items c e') := by
    intro e' he'
    have hq'lt : edgeItem s.g e' < s.items.size := by show 1 + s.g.nv + e' < _; omega
    rw [TEntry.edges, ha₃s, hae]
    constructor
    · rintro ⟨j, hj, h⟩
      rw [List.mem_singleton] at hj; subst hj
      rw [Items.EdgeBelow, hg, hBpj] at h
      rcases h with h | ⟨r, hr, h⟩
      · rcases hpa with hpi | hp
        · rw [hpi] at h; exact .inl h
        · exact absurd hq'lt (by rw [hpch hp _ h]; exact Nat.not_lt.2 hp)
      · rcases hS r hr with hrq | ⟨hri, -⟩ | hpr
        · rw [hrq] at h; exact .inr ((hce e').2 h)
        · rw [hri] at h; exact .inl h
        · rcases hpa with hpi | hp
          · rw [hpi] at hpr; exact .inl (Relation.ReflTransGen.head hpr h)
          · exact absurd hpr (by rw [Items.IsParent, ch_of_ge hp]; simp)
    · rintro (h | h)
      · refine ⟨p, List.mem_singleton_self _, ?_⟩
        rw [Items.EdgeBelow, hg, hBpj]
        rcases hpa with hpi | hp
        · rw [hpi]; exact .inl h
        · exact .inr ⟨i, hSi (fun h' => absurd hilt (by rw [h']; exact Nat.not_lt.2 hp)), h⟩
      · exact ⟨p, List.mem_singleton_self _, by
          rw [Items.EdgeBelow, hg, hBpj]; exact .inr ⟨q, hSq, (hce e').1 h⟩⟩
  have hTa₃ : ∀ x, s'.g.Touches (TEntry.edges s'.g s'.items a₃) x →
      s.g.Touches (TEntry.edges s.g s.items a) x ∨ s.g.Touches (TEntry.edges s.g s.items c) x := by
    rintro x ⟨e', he', hte, hx⟩
    rw [hg] at he' hx
    rcases (hce' e' he').1 hte with h | h
    · exact .inl ⟨e', he', h, hx⟩
    · exact .inr ⟨e', he', h, hx⟩
  have hTc : ∀ x, s.g.Touches (TEntry.edges s.g s.items c) x →
      s'.g.Touches (TEntry.edges s'.g s'.items a₃) x := by
    rintro x ⟨e', he', hte, hx⟩
    exact ⟨e', by rw [hg]; exact he', (hce' e' he').2 (.inr hte), by rw [hg]; exact hx⟩
  have hIntc : ∀ x, s.g.Interior (TEntry.edges s.g s.items c) x →
      s'.g.Interior (TEntry.edges s'.g s'.items a₃) x := by
    intro x h e' he' hx
    rw [hg] at he' hx
    exact (hce' e' he').2 (.inr (h e' he' hx))
  have hInt : ∀ x, s.g.Interior (TEntry.edges s.g s.items a) x →
      s'.g.Interior (TEntry.edges s'.g s'.items a₃) x := by
    intro x h e' he' hx
    rw [hg] at he' hx
    exact (hce' e' he').2 (.inl (h e' he' hx))
  have hCEa := hCE a hamA
  have hA1a : allType1 d a.topDepth done := by rw [had]; exact hA1
  have hCSa := hCS a hamA hA1a
  have hCEc := hCE c hcmA
  have hCSc := hCS c hcmA (by rw [hcd]; exact hA1)
  have hCE' : ∀ t ∈ tl, CtxEntry v d s t → CtxEntry v d s' t := by
    intro t ht h
    obtain ⟨e', h1, h2⟩ := h.nonempty
    exact
      { vStart := h.vStart
        depth := h.depth
        nonempty := ⟨e', by rw [hg]; exact h1, (hedgesT t ht e').2 h2⟩
        touch_bot := by rw [hedgesTF t ht, hg]; exact h.touch_bot
        touch_top := by rw [hedgesTF t ht, hg, hsv]; exact h.touch_top
        side := by rw [hsd]; exact h.side
        att := by rw [hedgesTF t ht, hg, hsv]; exact h.att }
  have hSnot : ∀ t ∈ tl, ∀ j ∈ t.spans.1 ++ t.spans.2, j ∉ S := fun t ht j hj hjS => by
    rcases hS j hjS with h | ⟨h, -⟩ | h
    · rw [h] at hj; exact hct_sd t ht _ hqm hj
    · rw [h] at hj; exact hat_sd t ht _ him hj
    · exact hC.span_root t (htlm t ht) j hj p h
  have hCS' : ∀ t ∈ tl, CtxSingle v s t → CtxSingle v s' t := by
    intro t ht h
    obtain ⟨j, hj, hroot⟩ := h.item
    have hjm : j ∈ t.spans.1 ++ t.spans.2 := by
      rw [hj]; unfold setSides; split <;> simp
    exact
      { item := ⟨j, by rw [hsd]; exact hj, fun p' hp' => by
          rcases (hP _ _).1 hp' with h' | ⟨-, hjS⟩
          · exact hroot p' h'
          · exact hSnot t ht j hjm hjS⟩
        att := by rw [hedgesTF t ht, hg, hsv]; exact h.att }
  have hproot' : ∀ p', ¬ Items.IsParent s'.items p' p := fun p' h => by
    rcases (hP _ _).1 h with h | ⟨-, h⟩
    · exact hproot p' h
    · exact hpS h
  have hCEa₃ : CtxEntry v d s' a₃ :=
    { vStart := rfl
      depth := hlt
      nonempty := by
        obtain ⟨e', h1, h2⟩ := hCEc.nonempty
        exact ⟨e', by rw [hg]; exact h1, (hce' e' h1).2 (.inr h2)⟩
      touch_bot := hTc v hCEc.touch_bot
      touch_top := by
        show s'.g.Touches _ s'.stackVerts[lv]!; rw [hsv, ← hcd]; exact hTc _ hCEc.touch_top
      side := by
        show getSide (setSides s.stackDir[lv]! [p] []) (!s'.stackDir[lv]!) = []
        rw [hsd]; unfold setSides getSide; cases s.stackDir[lv]! <;> rfl
      att := fun x hx hxv hxi => by
        rcases hTa₃ x hx with h | h
        · obtain ⟨k, hk1, hk2, hk3⟩ := hCEa.att x h hxv (fun hI => hxi (hInt x hI))
          exact ⟨k, by rw [had] at hk1; exact hk1, hk2, by rw [hsv]; exact hk3⟩
        · obtain ⟨k, hk1, hk2, hk3⟩ := hCEc.att x h hxv (fun hI => hxi (hIntc x hI))
          exact ⟨k, by rw [hcd] at hk1; exact hk1, hk2, by rw [hsv]; exact hk3⟩ }
  have hCSa₃ : CtxSingle v s' a₃ :=
    { item := ⟨p, by
        show setSides s.stackDir[lv]! [p] [] = setSides s'.stackDir[lv]! [p] []
        rw [hsd], hproot'⟩
      att := fun x hx => by
        rcases hTa₃ x hx with h | h
        · rcases hCSa.att x h with h | h | h
          · exact .inl h
          · exact .inr (.inl (by show x = s'.stackVerts[lv]!; rw [hsv, ← had]; exact h))
          · exact .inr (.inr (hInt x h))
        · rcases hCSc.att x h with h | h | h
          · exact .inl h
          · exact .inr (.inl (by show x = s'.stackVerts[lv]!; rw [hsv, ← hcd]; exact h))
          · exact .inr (.inr (hIntc x h)) }
  have hvd : s.stackVerts[d]! = v := hC.sv_d
  have hvne : ∀ k, k < d → v ≠ s.stackVerts[k]! := fun k hk h =>
    hC.path k d hk (Nat.le_refl _) (by rw [hvd, h])
  have hmem' : ∀ t ∈ s'.tstack, t = a₃ ∨ t ∈ tl := fun t ht => by
    rw [hts'] at ht; exact List.mem_cons.1 ht
  refine
    { top := ⟨a₃ :: (above'' ++ vt :: below), by rw [hts', htl]; rfl,
        { split := ⟨a₃ :: above'', below, ?_, ?_, ?_, ?_, ?_, ?_, hBel⟩
          bot := fun t ht k hk => by
            rw [hsv]
            rcases List.mem_cons.1 ht with rfl | ht
            · exact hvne k hk
            · exact hCT.bot t (hmemTop t ht) k hk
          vitems := fun t ht x hx hm' => by
            rw [hg] at hx
            rcases List.mem_cons.1 ht with rfl | ht
            · rw [ha₃s, List.mem_singleton] at hm'
              rcases hpa with hpi | hp
              · exact hCT.vitems a hamT x hx (by rw [has, List.mem_singleton]; exact hm'.trans hpi)
              · exact absurd (show vertItem x < s.items.size by show 1 + x < _; omega)
                  (by rw [hm']; exact Nat.not_lt.2 hp)
            · exact hCT.vitems t (hmemTop t ht) x hx hm'
          qitems := fun t ht e' he' hm' => by
            rw [hg] at he' hm'
            rcases List.mem_cons.1 ht with rfl | ht
            · rw [ha₃s, List.mem_singleton] at hm'
              rcases hpa with hpi | hp
              · exact hCT.qitems a hamT e' he' (by rw [has, List.mem_singleton]; exact hm'.trans hpi)
              · exact absurd (show edgeItem s.g e' < s.items.size by show 1 + s.g.nv + e' < _; omega)
                  (by rw [hm']; exact Nat.not_lt.2 hp)
            · exact hCT.qitems t (hmemTop t ht) e' he' hm'
          edges := fun t ht e' he' hte => by
            rw [hg] at he'
            rcases List.mem_cons.1 ht with rfl | ht
            · rcases (hce' e' he').1 hte with h | h
              · exact hCT.edges a hamT e' he' h
              · exact hCT.edges c hcmT e' he' h
            · exact hCT.edges t (hmemTop t ht) e' he' ((hedgesT t (hmemT' t ht) e').1 hte)
          cover := fun o' ho' hl e' he' hs => by
            rw [hg] at he'
            obtain ⟨t, ht, hte⟩ := hCT.cover o' ho' hl e' he' hs
            rw [htv, hab, List.cons_append, List.cons_append] at ht
            rcases List.mem_cons.1 ht with h | ht
            · rw [h] at hte
              exact ⟨a₃, List.mem_cons_self .., (hce' e' he').2 (.inr hte)⟩
            rcases List.mem_cons.1 ht with h | ht
            · rw [h] at hte
              exact ⟨a₃, List.mem_cons_self .., (hce' e' he').2 (.inl hte)⟩
            · exact ⟨t, List.mem_cons_of_mem _ ht, (hedgesT t (hmemT' t ht) e').2 hte⟩
          ret := fun _ => hC.hv_ret rfl }⟩
      base_edges := ⟨hC.base_edges.1, fun k hk e' he' => by
        rw [hg] at he'
        exact (hedgesT _ (hbase_tl k hk) e').trans (hC.base_edges.2 k hk e' he')⟩
      base_bot := hC.base_bot
      noVert_after := fun h => nomatch h
      vfirst := fun t ht hvt hnv => by
        rw [hnx]
        rcases hmem' t ht with rfl | ht
        · exact hC.vfirst a ham hav hva
        · exact hC.vfirst t (htlm t ht) hvt hnv
      sv := fun k hk => by rw [hsv]; exact hC.sv k hk
      sd := fun k hk => by rw [hsd]; exact hC.sd k hk
      sv_d := by rw [hsv]; exact hvd
      path := fun k k' h h' => by rw [hsv]; exact hC.path k k' h h'
      v_root := fun p' h => by
        rcases (hP _ _).1 h with h | ⟨-, h⟩
        · exact hC.v_root p' h
        · rcases hS _ h with h | ⟨h, -⟩ | h
          · exact hvq h
          · exact hva (by rw [h]; exact him)
          · exact hC.v_root p h
      vert_free := fun _ _ _ => rfl
      afterVert_ret := hC.afterVert_ret
      hv_ret := hC.hv_ret
      vert_book := fun h => nomatch h
      vert_disj := fun h => nomatch h
      vert_edges := fun e' he' => by rw [hg] at he'; rw [hEBv]; exact hC.vert_edges e' he'
      touch_bot := fun t ht hne => by
        rcases hmem' t ht with rfl | ht
        · exact hTc v hCEc.touch_bot
        · rw [hedgesTF t ht, hg]
          rw [hedgesTF t ht, hg] at hne
          exact hC.touch_bot t (htlm t ht) hne
      span_root := fun t ht j hj p' hp' => by
        rcases hmem' t ht with rfl | ht
        · rw [ha₃s, List.mem_singleton] at hj; subst hj; exact hproot' p' hp'
        · rcases (hP _ _).1 hp' with h | ⟨-, h⟩
          · exact hC.span_root t (htlm t ht) j hj p' h
          · exact hSnot t ht j hj h
      span_lt := fun t ht j hj => by
        rcases hmem' t ht with rfl | ht
        · rw [ha₃s, List.mem_singleton] at hj; subst hj; exact hplt
        · exact Nat.lt_of_lt_of_le (hC.span_lt t (htlm t ht) j hj) hsz'
      ch_lt := fun p' c' h => by
        rcases (hP _ _).1 h with h | ⟨-, h⟩
        · exact Nat.lt_of_lt_of_le (hC.ch_lt p' c' h) hsz'
        · rcases hS _ h with h | ⟨h, -⟩ | h
          · rw [h]; exact Nat.lt_of_lt_of_le hqlt hsz'
          · rw [h]; exact Nat.lt_of_lt_of_le hilt hsz'
          · exact Nat.lt_of_lt_of_le (hC.ch_lt p c' h) hsz'
      disj := by
        rw [hts']
        refine List.pairwise_cons.2 ⟨fun t' ht' e' he' h1 h2 => ?_, ?_⟩
        · rw [hg] at he'
          rw [hedgesT t' ht'] at h2
          rcases (hce' e' he').1 h1 with h | h
          · exact hat_d t' ht' e' he' h h2
          · exact hct_d t' ht' e' he' h h2
        · refine List.Pairwise.imp_of_mem (fun {t t'} ht ht' h e' he' h1 h2 => ?_) htl_d
          rw [hg] at he'
          exact h e' he' ((hedgesT t ht e').1 h1) ((hedgesT t' ht' e').1 h2)
      span_disj := by
        rw [hts']
        refine List.pairwise_cons.2 ⟨fun t' ht' j hj hj' => ?_, htl_sd⟩
        rw [ha₃s, List.mem_singleton] at hj; subst hj
        exact hpT t' ht' hj'
      q_fresh := fun o' ho' e' hs => by
        obtain ⟨h1, h2, h3⟩ := hC.q_fresh o' ho' e' hs
        have he'lt := hrest_e o' ho' e' hs
        have hq'lt : edgeItem s.g e' < s.items.size := by show 1 + s.g.nv + e' < _; omega
        have hq'p : edgeItem s.g e' ≠ p := by
          rcases hpa with hpi | hp
          · intro h; exact h3 a ham (by rw [h, hpi]; exact him)
          · intro h; exact absurd hq'lt (by rw [h]; exact Nat.not_lt.2 hp)
        rw [hg]
        refine ⟨fun p' hp' => ?_, by rw [hch _ hq'p]; exact h2, fun t ht hm' => ?_⟩
        · rcases (hP _ _).1 hp' with h | ⟨-, h⟩
          · exact h1 p' h
          · rcases hS _ h with h | ⟨h, -⟩ | h
            · exact h3 c hcm (by rw [h]; exact hqm)
            · exact h3 a ham (by rw [h]; exact him)
            · exact h1 p h
        · rcases hmem' t ht with rfl | ht
          · rw [ha₃s, List.mem_singleton] at hm'; exact hq'p hm'
          · exact h3 t (htlm t ht) hm'
      v_fresh := fun o' ho' e₁ cls₁ child h y hy => by
        obtain ⟨h1, h2, h3, h4, h5, h6⟩ := hC.v_fresh o' ho' e₁ cls₁ child h y hy
        have hyv : y < s.g.nv := hrest_v o' ho' e₁ cls₁ child h y hy
        have hyne : v ≠ y := fun hvy => h4 d (Nat.le_refl _) (hvd.trans hvy)
        have hylt : vertItem y < s.items.size := by show 1 + y < _; omega
        have hyp : vertItem y ≠ p := by
          rcases hpa with hpi | hp
          · intro h'; exact h3 a ham (by rw [h', hpi]; exact him)
          · intro h'; exact absurd hylt (by rw [h']; exact Nat.not_lt.2 hp)
        refine ⟨fun p' hp' => ?_, by rw [hch _ hyp]; exact h2, fun t ht hm' => ?_,
          fun k hk => by rw [hsv]; exact h4 k hk, fun t ht => ?_, fun t ht hT => ?_⟩
        · rcases (hP _ _).1 hp' with h' | ⟨-, h'⟩
          · exact h1 p' h'
          · rcases hS _ h' with h' | ⟨h', -⟩ | h'
            · exact h3 c hcm (by rw [h']; exact hqm)
            · exact h3 a ham (by rw [h']; exact him)
            · exact h1 p h'
        · rcases hmem' t ht with rfl | ht
          · rw [ha₃s, List.mem_singleton] at hm'; exact hyp hm'
          · exact h3 t (htlm t ht) hm'
        · rcases hmem' t ht with rfl | ht
          · exact hyne
          · exact h5 t (htlm t ht)
        · rcases hmem' t ht with rfl | ht
          · rcases hTa₃ y hT with hh | hh
            · exact h6 a ham hh
            · exact h6 c hcm hh
          · rw [hedgesTF t ht, hg] at hT
            exact h6 t (htlm t ht) hT }
  · rw [if_pos rfl]
    refine ⟨vt, rfl, hvt1, hvt2, hvt3, fun t ht => ?_, fun o' ho' => ?_,
      fun o' ho' e' he' hs => ?_⟩
    · rcases List.mem_cons.1 ht with rfl | ht
      · rw [ha₃s, List.mem_singleton]; exact fun h => hpv h.symm
      · exact hvt4 t (habA t ht)
    · rcases hvt5 o' ho' with ⟨t, ht, h⟩ | h
      · rw [hab] at ht
        rcases List.mem_cons.1 ht with h' | ht
        · rw [h'] at h; exact .inl ⟨a₃, List.mem_cons_self .., by show lv = _; rw [← hcd]; exact h⟩
        rcases List.mem_cons.1 ht with h' | ht
        · rw [h'] at h; exact .inl ⟨a₃, List.mem_cons_self .., by show lv = _; rw [← had]; exact h⟩
        · exact .inl ⟨t, List.mem_cons_of_mem _ ht, h⟩
      · exact .inr h
    · rw [hg] at he'
      rcases hvt6 o' ho' e' he' hs with ⟨t, ht, h⟩ | h
      · rw [hab] at ht
        rcases List.mem_cons.1 ht with h' | ht
        · rw [h'] at h
          exact .inl ⟨a₃, List.mem_cons_self .., (hce' e' he').2 (.inr h)⟩
        rcases List.mem_cons.1 ht with h' | ht
        · rw [h'] at h
          exact .inl ⟨a₃, List.mem_cons_self .., (hce' e' he').2 (.inl h)⟩
        · exact .inl ⟨t, List.mem_cons_of_mem _ ht, (hedgesT t (hab_tl t ht) e').2 h⟩
      · exact .inr ((hedgesT vt hvt_tl e').2 h)
  · intro t ht
    rcases List.mem_cons.1 ht with rfl | ht
    · exact hCEa₃
    · exact hCE' t (hab_tl t ht) (hCE t (habA t ht))
  · intro t ht hA
    rcases List.mem_cons.1 ht with rfl | ht
    · exact hCSa₃
    · exact hCS' t (hab_tl t ht) (hCS t (habA t ht) hA)
  · have hPW' := hPW
    rw [hab] at hPW'
    refine List.pairwise_cons.2 ⟨fun t' ht' => ?_,
      (List.pairwise_cons.1 (List.pairwise_cons.1 hPW').2).2⟩
    obtain ⟨h1, h2⟩ := (List.pairwise_cons.1 (List.pairwise_cons.1 hPW').2).1 t' ht'
    exact ⟨by show t'.topDepth ≤ lv; rw [← had]; exact h1, h2⟩
  · intro t ht
    rcases List.mem_cons.1 ht with rfl | ht
    · obtain ⟨o', ho', h⟩ := hLow a hamA
      exact ⟨o', ho', by show lv = _; rw [← had]; exact h⟩
    · exact hLow t (habA t ht)
  · intro t ht e' he' hte
    rw [hg] at he'
    rcases List.mem_cons.1 ht with rfl | ht
    · rcases (hce' e' he').1 hte with h | h
      · exact hEdg a hamA e' he' h
      · exact hEdg c hcmA e' he' h
    · exact hEdg t (habA t ht) e' he' ((hedgesT t (hab_tl t ht) e').1 hte)

/-- After a returning back-edge out whose P-check fires: from the state after the `Q e` push
(`s₂`, described by its equations), when the entry below the new top is the parallel ear
`(v, lowval)` (`condP`), `finishRest` (P-merge; no vertex push since `hasVert || push = true`)
re-establishes the context. `condP` forces `push = false` (the pushed `V v` entry would sit at
depth `d ≠ lowval`), `mergeP_single` gives the merged side, `earCtx_mergeP` the post-state. -/
theorem ctx_step_back_ret_P {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s s₂ : WalkState} {e dest : Nat} {cls : OutClass}
    (hC : EarCtx v d done (.back e dest cls :: rest) hasVert base bE sv sd s)
    (hv : v < s.g.nv) (he : e < s.g.ne) (hlt : cls.lowval d < d)
    (hsz : 1 + s.g.nv + s.g.ne ≤ s.items.size)
    (hrank : ∀ o' ∈ done, o'.1.cls.rank ≤ cls.rank)
    (hinc : ∀ o' ∈ done, o'.1.e < s.g.ne ∧ s.g.Inc o'.1.e v)
    (hend : Items.PairEq (v, dest) s.g.edges[e]!) (hdest : s.stackVerts[cls.lowval d]! = dest)
    (hrest_e : ∀ o' ∈ rest, ∀ e', subEdges o' e' → e' ≠ e ∧ e' < s.g.ne)
    (hrest_v : ∀ o' ∈ rest, ∀ e₁ cls₁ child, o' = .tree e₁ cls₁ child →
      ∀ y ∈ child.verts, y < s.g.nv)
    (hnc : ∀ o' ∈ done, ∀ e', subEdges o'.1 e' → ∀ x, s.g.Inc e' x →
      ∀ o'' ∈ rest, ∀ e₁ cls₁ child, o'' = .tree e₁ cls₁ child → x ∉ child.verts)
    (L : List TEntry) (push dirD : Bool)
    (hpush : push = true ↔ hasVert = false)
    (hL : L = if push then [⟨v, d, s.nxtEdgeIdx, setSides dirD [vertItem v] []⟩] else [])
    (hg : s₂.g = s.g) (hsv : s₂.stackVerts = s.stackVerts)
    (hsd : ∀ k, k < d → s₂.stackDir[k]! = s.stackDir[k]!)
    (hsz' : s₂.items.size = s.items.size)
    (hch : ∀ j, Items.ch s₂.items j = Items.ch s.items j)
    (hnx : s₂.nxtEdgeIdx = s.nxtEdgeIdx + 1)
    (hts : s₂.tstack = ⟨v, cls.lowval d, s.nxtEdgeIdx,
      setSides s.stackDir[cls.lowval d]! [edgeItem s.g e] []⟩ :: (L ++ s.tstack))
    (hc : result (condP v (cls.lowval d) cls.isType1) s₂ = true) :
    wp (finishRest v d (cls.lowval d) cls.isType1 true true)
      (fun hv' s' => EarCtx v d (done ++ [(.back e dest cls, true)]) rest true base bE sv sd s₂ →
        EarCtx v d (done ++ [(.back e dest cls, true)]) rest hv' base bE sv sd s' ∧
        OutFrame v (.back e dest cls) base s₂ s' ∧
        (CtxShape v d done hasVert base s →
          CtxShape v d (done ++ [(.back e dest cls, true)]) hv' base s')) s₂ := by
  have hc' : (cls.isType1 && decide (s₂.tstack.length ≥ 2) && (s₂.tstack.tail.head!.vStart == v) &&
      (s₂.tstack.tail.head!.topDepth == cls.lowval d)) = true := by
    rw [result, run_condP] at hc; exact hc
  have hc2 := hc'
  simp only [Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq] at hc2
  obtain ⟨⟨⟨ht1, hlen⟩, hnv⟩, hnd⟩ := hc2
  have hL0 : L = [] := by
    rw [hL]
    cases hp : push
    · rfl
    · exfalso
      rw [hL, if_pos hp] at hts
      rw [hts] at hnd
      simp only [List.tail_cons, List.head!_cons] at hnd
      exact absurd hnd (Nat.ne_of_gt hlt)
  rw [hL0, List.nil_append] at hts
  have hp0 : push = false := by
    cases hp : push
    · rfl
    · exfalso; simp [hp, hL0] at hL
  have hhvT : hasVert = true := by
    cases hhv : hasVert
    · exact absurd (hpush.2 hhv) (by rw [hp0]; decide)
    · rfl
  subst hhvT
  cases hst : s.tstack with
  | nil =>
    exfalso
    rw [hts, hst] at hlen
    simp at hlen
  | cons a tl =>
  rw [hst] at hts
  rw [hts] at hnv hnd
  simp only [List.tail_cons, List.head!_cons] at hnv hnd
  obtain ⟨tl', rfl⟩ : ∃ tl', tl = tl' ++ base := by
    obtain ⟨top, htop, -⟩ := hC.top
    cases top with
    | nil =>
      exfalso
      rw [hst] at htop
      have hb : base = a :: tl := by simpa using htop.symm
      exact hC.base_bot a (by rw [hb]; exact List.mem_cons_self ..) hnv
    | cons a' tl' =>
      rw [hst, List.cons_append, List.cons.injEq] at htop
      exact ⟨tl', htop.2⟩
  set lv := cls.lowval d with hlv
  set q := edgeItem s.g e with hq
  set c : TEntry := ⟨v, lv, s.nxtEdgeIdx, setSides s.stackDir[lv]! [q] []⟩ with hcdef
  have hsdl : s₂.stackDir[lv]! = s.stackDir[lv]! := hsd lv hlt
  have hcs : c.spans.1 ++ c.spans.2 = [q] := by
    show (setSides s.stackDir[lv]! [q] []).1 ++ (setSides s.stackDir[lv]! [q] []).2 = [q]
    unfold setSides; split <;> rfl
  have hA1 : allType1 d lv (done ++ [(.back e dest cls, true)]) := by
    intro o ho hol
    unfold afterVert at ho
    rw [List.filter_append, List.map_append] at ho
    rcases List.mem_append.1 ho with ho | ho
    · have hr := hrank (o, true) (mem_afterVert ho)
      obtain ⟨lv', k', hk', -⟩ := ret_of_lowval_lt (o := o) (d := d) (by rw [hol]; exact hlt)
      obtain ⟨lv'', k'', hk'', -⟩ := ret_of_lowval_lt (o := DfsOut.back e dest cls) (d := d) hlt
      simp only [DfsOut.cls] at hk''
      rw [hk'] at hol hr ⊢
      rw [hk''] at hr hlv ht1
      simp only [OutClass.lowval] at hol hlv
      simp only [OutClass.rank] at hr
      cases k' <;> cases k'' <;> simp only [RetKind.rank, OutClass.isType1] at hr ht1 ⊢ <;> omega
    · simp [afterVert] at ho
      rw [ho]; exact ht1
  unfold finishRest finishP finishTail condP maybeUnwrapNxt finishTstackTop
  simp only [wp_bind, wp_tstackSize, wp_nxt, wp_pure, wp_ite, hc', ↓reduceIte, Bool.not_true,
    Bool.false_eq_true, wp_get, wp_allocItem, wp_stackDir, wp_getItem, wp_modifyNxt,
    wp_mergeTstackTops, wp_cur, wp_makeVs, wp_modifyItem, wp_modifyCur]
  have hts₂ : s₂.tstack = c :: a :: (tl' ++ base) := hts
  have hcd : c.topDepth = lv := rfl
  rw [hts₂]
  simp only [mergeTop, List.tail_cons, List.head!_cons, hnd, hnv, hcd, Nat.min_self]
  have hcsp : c.spans = setSides s₂.stackDir[lv]! [q] [] := by
    show setSides s.stackDir[lv]! [q] [] = _
    rw [hsdl]
  have hv₂ : v < s₂.g.nv := by rw [hg]; exact hv
  have hsz₂ : 1 + s₂.g.nv + s₂.g.ne ≤ s₂.items.size := by rw [hg, hsz']; exact hsz
  have hrest_e₂ : ∀ o' ∈ rest, ∀ e', subEdges o' e' → e' < s₂.g.ne :=
    fun o' ho' e' h => by rw [hg]; exact (hrest_e o' ho' e' h).2
  have hrest_v₂ : ∀ o' ∈ rest, ∀ e₁ cls₁ child, o' = .tree e₁ cls₁ child →
      ∀ y ∈ child.verts, y < s₂.g.nv :=
    fun o' ho' e₁ cls₁ child h y hy => by rw [hg]; exact hrest_v o' ho' e₁ cls₁ child h y hy
  have hcV : vertItem v ∉ c.spans.1 ++ c.spans.2 := by
    rw [hcs, List.mem_singleton, hq]; exact vertItem_ne_edgeItem' hv
  have ham : a ∈ s₂.tstack := by rw [hts₂]; exact List.mem_cons_of_mem _ (List.mem_cons_self ..)
  split_ifs with h1 h2 <;> intro hC₂ <;>
    (obtain ⟨i, hai, hiroot⟩ := mergeP_single hC₂ hts₂ hcV hnv hnd hlt hA1
     have hilt : i < s₂.items.size :=
       hC₂.span_lt a ham i (by rw [hai]; unfold setSides; split <;> simp)
     have hhead : (getSide a.spans s₂.stackDir[lv]!).head! = i := by
       rw [hai]; unfold setSides getSide; cases s₂.stackDir[lv]! <;> rfl
     have hchi : s₂.items[i]!.ch = Items.ch s₂.items i := by
       rw [Items.ch_eq_getElem hilt, getElem!_pos s₂.items i hilt])
  on_goal 2 =>
    rw [hhead, hchi]
    set S := getSide (c.spans.1 ++ (setSides s₂.stackDir[lv]! (Items.ch s₂.items i) []).1,
      (setSides s₂.stackDir[lv]! (Items.ch s₂.items i) []).2 ++ c.spans.2) s₂.stackDir[lv]! with hSdef
    have hSmem : ∀ r, r ∈ S ↔ r = q ∨ r ∈ Items.ch s₂.items i := by
      intro r; rw [hSdef, hcsp]; unfold setSides getSide
      cases s₂.stackDir[lv]! <;> simp <;> tauto
    refine (fun hS hP => ⟨earCtx_mergeP hC₂ hv₂ hlt hsz₂ hts₂ hcd hcs rfl hnv hnd hai hA1
      hrest_e₂ hrest_v₂ (p := i) (S := S) (.inl rfl) hS ((hSmem q).2 (.inl rfl))
      (fun h => absurd rfl h) hP (fun j hj => Items.ch_modify_of_ne _ _ hj)
      (by simp) (by rw [Array.size_modify]; exact hilt) rfl rfl rfl rfl rfl,
      outFrame_mergeP hts₂ hcs (by rw [hai, mem_setSides]; exact List.mem_singleton_self _)
        (.inl rfl) hS hP (fun j hj => Items.ch_modify_of_ne _ _ hj) rfl (by simp) rfl,
      fun hS0 => ctxShape_mergeP hC hS0 hv hlt hsz hst hnv hnd (p := i) (hts' := rfl)
        (hsd := fun k hk => hsd k hk)
        (Or.inl (by rw [hai, mem_setSides]; exact List.mem_singleton_self _)) rfl rfl rfl⟩) ?_ ?_
    · intro r hr
      rcases (hSmem r).1 hr with h | h
      · exact .inl h
      · exact .inr (.inr h)
    · intro p' c'
      unfold Items.IsParent
      by_cases hp' : p' = i
      · subst hp'
        rw [Items.ch_modify_at _ _ hilt]
        simp only [hSmem, true_and]
        tauto
      · rw [Items.ch_modify_of_ne _ _ hp']
        simp [hp']
  all_goals
    set S := getSide (c.spans.1 ++ a.spans.1, a.spans.2 ++ c.spans.2) s₂.stackDir[lv]! with hSdef
    have hSmem : ∀ r, r ∈ S ↔ r = q ∨ r = i := by
      intro r; rw [hSdef, hcsp, hai]; unfold setSides getSide
      cases s₂.stackDir[lv]! <;> simp <;> tauto
    refine (fun hS hP => ⟨earCtx_mergeP hC₂ hv₂ hlt hsz₂ hts₂ hcd hcs rfl hnv hnd hai hA1
      hrest_e₂ hrest_v₂ (p := s₂.items.size) (S := S) (.inr (Nat.le_refl _)) hS
      ((hSmem q).2 (.inl rfl)) (fun _ => (hSmem i).2 (.inr rfl)) hP
      (fun j hj => by rw [Items.ch_modify_of_ne _ _ hj, Items.ch_push, if_neg hj])
      (by simp) (by simp) rfl rfl rfl rfl rfl,
      outFrame_mergeP hts₂ hcs (by rw [hai, mem_setSides]; exact List.mem_singleton_self _)
        (.inr (Nat.le_refl _)) hS hP
        (fun j hj => by rw [Items.ch_modify_of_ne _ _ hj, Items.ch_push, if_neg hj])
        rfl (by simp) rfl,
      fun hS0 => ctxShape_mergeP hC hS0 hv hlt hsz hst hnv hnd (p := s₂.items.size) (hts' := rfl)
        (hsd := fun k hk => hsd k hk) (Or.inr (Nat.le_of_eq hsz'.symm)) rfl rfl rfl⟩) ?_ ?_
    · intro r hr
      rcases (hSmem r).1 hr with h | h
      · exact .inl h
      · exact .inr (.inl ⟨h, Nat.ne_of_lt hilt⟩)
    · intro p' c'
      unfold Items.IsParent
      by_cases hp' : p' = s₂.items.size
      · subst hp'
        rw [Items.ch_modify_at _ _ (by simp), ch_of_ge (Nat.le_refl _), Array.getElem_push_eq]
        simp
      · rw [Items.ch_modify_of_ne _ _ hp', Items.ch_push, if_neg hp']
        simp [hp']

/-- `finishRest` after the `Q e` push of a returning back edge: `by_cases` on the P-check; the
no-merge case is `earCtx_pushBack`, the merge case `ctx_step_back_ret_P`. -/
theorem ctx_step_back_ret_rest {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s s₂ : WalkState} {e dest : Nat} {cls : OutClass}
    (hC : EarCtx v d done (.back e dest cls :: rest) hasVert base bE sv sd s)
    (hv : v < s.g.nv) (he : e < s.g.ne) (hlt : cls.lowval d < d)
    (hsz : 1 + s.g.nv + s.g.ne ≤ s.items.size)
    (hrank : ∀ o' ∈ done, o'.1.cls.rank ≤ cls.rank)
    (hinc : ∀ o' ∈ done, o'.1.e < s.g.ne ∧ s.g.Inc o'.1.e v)
    (hend : Items.PairEq (v, dest) s.g.edges[e]!) (hdest : s.stackVerts[cls.lowval d]! = dest)
    (hrest_e : ∀ o' ∈ rest, ∀ e', subEdges o' e' → e' ≠ e ∧ e' < s.g.ne)
    (hrest_v : ∀ o' ∈ rest, ∀ e₁ cls₁ child, o' = .tree e₁ cls₁ child →
      ∀ y ∈ child.verts, y < s.g.nv)
    (hnc : ∀ o' ∈ done, ∀ e', subEdges o'.1 e' → ∀ x, s.g.Inc e' x →
      ∀ o'' ∈ rest, ∀ e₁ cls₁ child, o'' = .tree e₁ cls₁ child → x ∉ child.verts)
    (L : List TEntry) (push dirD : Bool) (ht1 : cls.isType1 = true)
    (hdirD : push = true → dirD = !s.stackDir[cls.lowval d]!)
    (hpush : push = true ↔ hasVert = false)
    (hL : L = if push then [⟨v, d, s.nxtEdgeIdx, setSides dirD [vertItem v] []⟩] else [])
    (hg : s₂.g = s.g) (hsv : s₂.stackVerts = s.stackVerts)
    (hsd : ∀ k, k < d → s₂.stackDir[k]! = s.stackDir[k]!)
    (hsz' : s₂.items.size = s.items.size)
    (hch : ∀ j, Items.ch s₂.items j = Items.ch s.items j)
    (hnx : s₂.nxtEdgeIdx = s.nxtEdgeIdx + 1)
    (hts : s₂.tstack = ⟨v, cls.lowval d, s.nxtEdgeIdx,
      setSides s.stackDir[cls.lowval d]! [edgeItem s.g e] []⟩ :: (L ++ s.tstack)) :
    wp (finishRest v d (cls.lowval d) cls.isType1 (hasVert || push) true)
      (fun hv' s' => EarCtx v d (done ++ [(.back e dest cls, hasVert || push)]) rest hv'
        base bE sv sd s' ∧ OutFrame v (.back e dest cls) base s s' ∧
        (CtxShape v d done hasVert base s →
          CtxShape v d (done ++ [(.back e dest cls, hasVert || push)]) hv' base s')) s₂ := by
  have hhv' : (hasVert || push) = true := by
    cases h : hasVert
    · rw [hpush.2 h]; rfl
    · rfl
  rw [hhv']
  have hpb := earCtx_pushBack hC hv he hlt hsz hrank hinc hend hdest hrest_e hrest_v hnc L push dirD
    hpush hL hg hsv hsd hsz' hch hnx hts
  have hPB : OutFrame v (.back e dest cls) base s s₂ := by
    obtain ⟨top, htop, -⟩ := hC.top
    exact outFrame_pushBack htop hL hg hsz' hch hts
  by_cases hc : result (condP v (cls.lowval d) cls.isType1) s₂ = true
  · exact wp_mono _ (ctx_step_back_ret_P hC hv he hlt hsz hrank hinc hend hdest hrest_e hrest_v hnc
      L push dirD hpush hL hg hsv hsd hsz' hch hnx hts hc) fun _ _ h =>
        ⟨(h hpb).1, hPB.trans (h hpb).2.1, (h hpb).2.2⟩
  · have hc' : (cls.isType1 && decide (s₂.tstack.length ≥ 2) && (s₂.tstack.tail.head!.vStart == v) &&
        (s₂.tstack.tail.head!.topDepth == cls.lowval d)) = false := Bool.eq_false_iff.2 hc
    unfold finishRest finishP finishTail condP
    simp only [wp_bind, wp_tstackSize, wp_nxt, wp_pure, wp_ite, hc', Bool.false_eq_true,
      ↓reduceIte, Bool.not_true]
    have hnoP : push = false → ∀ a tl, s.tstack = a :: tl → a.vStart = v →
        a.topDepth ≠ cls.lowval d := by
      intro hp a tl hst hav hd
      rw [hts, hL, hp, ite_eq_right Bool.false_ne_true, List.nil_append, hst] at hc'
      simp [hav, hd, ht1] at hc'
    exact ⟨hpb, hPB, fun hS => ctxShape_pushBack hC hS hv hlt hrank L push dirD hpush hL hdirD hnoP hsd hts⟩

/-- After a returning back-edge out (`lowval < d`): `finishBack` (the `Q` push, the P-check, the
first-edge vertex push) re-establishes the context. -/
theorem ctx_step_back_ret {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s : WalkState} {e dest : Nat} {cls : OutClass}
    (hC : EarCtx v d done (.back e dest cls :: rest) hasVert base bE sv sd s)
    (hv : v < s.g.nv) (hsd : d < s.stackDir.size)
    (hrank : ∀ o' ∈ done, o'.1.cls.rank ≤ cls.rank)
    (hinc : ∀ o' ∈ done, o'.1.e < s.g.ne ∧ s.g.Inc o'.1.e v)
    (hnd : ∀ o' ∈ done, ∀ e', subEdges o'.1 e' → e' ≠ e)
    (hb : cls.isTree = false) (hcls : ∀ lv k, cls = .ret lv k → lv < d)
    (he : e < s.g.ne) (hsz : 1 + s.g.nv + s.g.ne ≤ s.items.size)
    (hend : Items.PairEq (v, dest) s.g.edges[e]!) (hself : d ≤ cls.lowval d → dest = v)
    (hrest_e : ∀ o' ∈ rest, ∀ e', subEdges o' e' → e' ≠ e ∧ e' < s.g.ne)
    (hrest_v : ∀ o' ∈ rest, ∀ e₁ cls₁ child, o' = .tree e₁ cls₁ child →
      ∀ y ∈ child.verts, y < s.g.nv)
    (hdest : s.stackVerts[cls.lowval d]! = dest)
    (hnc : ∀ o' ∈ done, ∀ e', subEdges o'.1 e' → ∀ x, s.g.Inc e' x →
      ∀ o'' ∈ rest, ∀ e₁ cls₁ child, o'' = .tree e₁ cls₁ child → x ∉ child.verts)
    (L : List TEntry) (push : Bool)
    (hpush : push = true ↔ hasVert = false ∧ cls.lowval d < d ∧ cls.isType1 = true)
    (hL : L = if push then
      [⟨v, d, s.nxtEdgeIdx, setSides (s.stackDir.set! d (if cls.lowval d ≥ d then false
        else !s.stackDir[cls.lowval d]!))[d]! [vertItem v] []⟩] else [])
    (hlt : cls.lowval d < d) :
    wp (finishEdge v d (.back e dest cls) (L ++ s.tstack).length (hasVert || push))
      (fun hv' s' => EarCtx v d (done ++ [(.back e dest cls, hasVert || push)]) rest hv' base bE sv sd s' ∧
        OutFrame v (.back e dest cls) base s s' ∧
        (CtxShape v d done hasVert base s →
          CtxShape v d (done ++ [(.back e dest cls, hasVert || push)]) hv' base s'))
      { s with stackDir := s.stackDir.set! d (if cls.lowval d ≥ d then false
                 else !s.stackDir[cls.lowval d]!),
               tstack := L ++ s.tstack } := by
  have ht1 : cls.isType1 = true := by
    obtain ⟨lv', k, hc', -⟩ := ret_of_lowval_lt (o := DfsOut.back e dest cls) hlt
    simp only [DfsOut.cls] at hc'
    subst hc'
    cases k
    · rfl
    · rfl
    · exact nomatch hb
  have hpush' : push = true ↔ hasVert = false :=
    ⟨fun h => (hpush.1 h).1, fun h => hpush.2 ⟨h, hlt, ht1⟩⟩
  have hge' : ¬ (cls.lowval d ≥ d) := Nat.not_le.2 hlt
  rw [finishEdge_eq]
  simp only [finishEdge', wp_bind, wp_get, wp_stackDir, DfsOut.cls, DfsOut.e, DfsOut.dest, hge',
    ↓reduceIte, wp_makeVs, wp_modifyItem, hb, Bool.false_eq_true]
  unfold finishBack
  simp only [wp_bind, wp_pushEdgeTstack, wp_modify, DfsOut.cls, DfsOut.e]
  refine ctx_step_back_ret_rest hC hv he hlt hsz hrank hinc hend hdest hrest_e hrest_v hnc L push _ ht1
    (fun _ => by rw [getElem!_set!_self' _ _ _ hsd, ite_eq_right hge']) hpush' hL rfl rfl (fun k hk => getElem!_set!_ne' _ _ _ _ (Nat.ne_of_lt hk)) (by simp) ?_ rfl ?_
  · intro j
    simp only [Items.ch, Array.getElem?_modify]
    by_cases hj : edgeItem s.g e = j
    · subst hj; cases s.items[edgeItem s.g e]? <;> simp
    · simp [hj]
  · rw [getElem!_set!_ne' _ _ _ _ (Nat.ne_of_lt hlt)]

/-- `finishEdge` at a back-edge site re-establishes the context with the out appended to `done`
(`ctxCheck` after a back-edge out): the boundary case is proved, the returning case is
`ctx_step_back_ret`. -/
theorem ctx_step_back_frame {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s : WalkState} {e dest : Nat} {cls : OutClass}
    (hC : EarCtx v d done (.back e dest cls :: rest) hasVert base bE sv sd s)
    (hv : v < s.g.nv) (hsd : d < s.stackDir.size)
    (hrank : ∀ o' ∈ done, o'.1.cls.rank ≤ cls.rank)
    (hinc : ∀ o' ∈ done, o'.1.e < s.g.ne ∧ s.g.Inc o'.1.e v)
    (hnd : ∀ o' ∈ done, ∀ e', subEdges o'.1 e' → e' ≠ e)
    (hb : cls.isTree = false) (hcls : ∀ lv k, cls = .ret lv k → lv < d)
    (he : e < s.g.ne) (hsz : 1 + s.g.nv + s.g.ne ≤ s.items.size)
    (hend : Items.PairEq (v, dest) s.g.edges[e]!) (hself : d ≤ cls.lowval d → dest = v)
    (hrest_e : ∀ o' ∈ rest, ∀ e', subEdges o' e' → e' ≠ e ∧ e' < s.g.ne)
    (hrest_v : ∀ o' ∈ rest, ∀ e₁ cls₁ child, o' = .tree e₁ cls₁ child →
      ∀ y ∈ child.verts, y < s.g.nv)
    (hdest : cls.lowval d < d → s.stackVerts[cls.lowval d]! = dest)
    (hnc : ∀ o' ∈ done, ∀ e', subEdges o'.1 e' → ∀ x, s.g.Inc e' x →
      ∀ o'' ∈ rest, ∀ e₁ cls₁ child, o'' = .tree e₁ cls₁ child → x ∉ child.verts)
    (L : List TEntry) (push : Bool)
    (hpush : push = true ↔ hasVert = false ∧ cls.lowval d < d ∧ cls.isType1 = true)
    (hL : L = if push then
      [⟨v, d, s.nxtEdgeIdx, setSides (s.stackDir.set! d (if cls.lowval d ≥ d then false
        else !s.stackDir[cls.lowval d]!))[d]! [vertItem v] []⟩] else []) :
    wp (finishEdge v d (.back e dest cls) (L ++ s.tstack).length (hasVert || push))
      (fun hv' s' => EarCtx v d (done ++ [(.back e dest cls, hasVert || push)]) rest hv' base bE sv sd s' ∧
        OutFrame v (.back e dest cls) base s s' ∧
        (CtxShape v d done hasVert base s →
          CtxShape v d (done ++ [(.back e dest cls, hasVert || push)]) hv' base s'))
      { s with stackDir := s.stackDir.set! d (if cls.lowval d ≥ d then false
                 else !s.stackDir[cls.lowval d]!),
               tstack := L ++ s.tstack } := by
  by_cases hge : d ≤ cls.lowval d
  · exact ctx_step_back_boundary hC hv hsd hrank hinc hnd hb hcls he hsz hend hself hrest_e hrest_v
      L push hpush hL hge
  · exact ctx_step_back_ret hC hv hsd hrank hinc hnd hb hcls he hsz hend hself hrest_e hrest_v
      (hdest (Nat.lt_of_not_le hge)) hnc L push hpush hL (Nat.lt_of_not_le hge)

/-- `ctx_step_back_frame` without the item frame. -/
theorem ctx_step_back {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s : WalkState} {e dest : Nat} {cls : OutClass}
    (hC : EarCtx v d done (.back e dest cls :: rest) hasVert base bE sv sd s)
    (hv : v < s.g.nv) (hsd : d < s.stackDir.size)
    (hrank : ∀ o' ∈ done, o'.1.cls.rank ≤ cls.rank)
    (hinc : ∀ o' ∈ done, o'.1.e < s.g.ne ∧ s.g.Inc o'.1.e v)
    (hnd : ∀ o' ∈ done, ∀ e', subEdges o'.1 e' → e' ≠ e)
    (hb : cls.isTree = false) (hcls : ∀ lv k, cls = .ret lv k → lv < d)
    (he : e < s.g.ne) (hsz : 1 + s.g.nv + s.g.ne ≤ s.items.size)
    (hend : Items.PairEq (v, dest) s.g.edges[e]!) (hself : d ≤ cls.lowval d → dest = v)
    (hrest_e : ∀ o' ∈ rest, ∀ e', subEdges o' e' → e' ≠ e ∧ e' < s.g.ne)
    (hrest_v : ∀ o' ∈ rest, ∀ e₁ cls₁ child, o' = .tree e₁ cls₁ child →
      ∀ y ∈ child.verts, y < s.g.nv)
    (hdest : cls.lowval d < d → s.stackVerts[cls.lowval d]! = dest)
    (hnc : ∀ o' ∈ done, ∀ e', subEdges o'.1 e' → ∀ x, s.g.Inc e' x →
      ∀ o'' ∈ rest, ∀ e₁ cls₁ child, o'' = .tree e₁ cls₁ child → x ∉ child.verts)
    (L : List TEntry) (push : Bool)
    (hpush : push = true ↔ hasVert = false ∧ cls.lowval d < d ∧ cls.isType1 = true)
    (hL : L = if push then
      [⟨v, d, s.nxtEdgeIdx, setSides (s.stackDir.set! d (if cls.lowval d ≥ d then false
        else !s.stackDir[cls.lowval d]!))[d]! [vertItem v] []⟩] else []) :
    wp (finishEdge v d (.back e dest cls) (L ++ s.tstack).length (hasVert || push))
      (fun hv' s' => EarCtx v d (done ++ [(.back e dest cls, hasVert || push)]) rest hv' base bE sv sd s')
      { s with stackDir := s.stackDir.set! d (if cls.lowval d ≥ d then false
                 else !s.stackDir[cls.lowval d]!),
               tstack := L ++ s.tstack } :=
  wp_mono _ (ctx_step_back_frame hC hv hsd hrank hinc hnd hb hcls he hsz hend hself hrest_e hrest_v
    hdest hnc L push hpush hL) fun _ _ h => h.1

/-- `Items.Below` from an item of a set `Old` closed under parenthood in `I` is unchanged when the
children of every `Old` item are unchanged. -/
theorem below_kept {I I' : Items} (Old : ItemId → Prop)
    (hcl : ∀ a c, Old a → Items.IsParent I a c → Old c)
    (hk : ∀ j, Old j → Items.ch I' j = Items.ch I j) {a : ItemId} (ha : Old a) (i : ItemId) :
    Items.Below I' a i ↔ Items.Below I a i := by
  have h1 : ∀ i, Items.Below I' a i → Items.Below I a i ∧ Old i := by
    intro i h
    induction h with
    | refl => exact ⟨.refl, ha⟩
    | tail _ hs ih =>
      have hs' := (Iff.of_eq (congrArg (_ ∈ ·) (hk _ ih.2))).1 hs
      exact ⟨ih.1.tail hs', hcl _ _ ih.2 hs'⟩
  have h2 : ∀ i, Items.Below I a i → Items.Below I' a i ∧ Old i := by
    intro i h
    induction h with
    | refl => exact ⟨.refl, ha⟩
    | tail _ hs ih =>
      refine ⟨ih.1.tail ?_, hcl _ _ ih.2 hs⟩
      rw [Items.IsParent, hk _ ih.2]; exact hs
  exact ⟨fun h => (h1 i h).1, fun h => (h2 i h).1⟩

/-- A child walk does not change the edge set below the `V` item of a vertex outside the subtree
(dump-checked `kept_v_below`): from the kept frame `hk`, since the items below `V v` before the
walk are neither the tree's items (those are roots) nor allocated by it (`below_kept`). -/
theorem walkTree_below_kept (t : DfsTree) (d : Nat) (s : WalkState)
    (hk : wp (walkTree t d) (fun _ s' => KeptRel s.g s t.verts t.edges s') s) {g : Graph} (hT : Types g s)
    (v : Nat) (hv : v ∉ t.verts) (hv' : v < s.g.nv)
    (hch : ∀ p c, Items.IsParent s.items p c → c < s.items.size)
    (hvr : ∀ x ∈ t.verts, ∀ p, ¬ Items.IsParent s.items p (vertItem x))
    (her : ∀ e ∈ t.edges, ∀ p, ¬ Items.IsParent s.items p (edgeItem s.g e)) :
    wp (walkTree t d) (fun _ s' => ∀ e, e < s.g.ne →
      (Items.EdgeBelow s'.g s'.items (vertItem v) e ↔ Items.EdgeBelow s.g s.items (vertItem v) e)) s := by
  refine wp_mono _ (wp_and (walkTree_frame t d s hT) hk)
    fun _ s' h e he => ?_
  obtain ⟨⟨hg, -⟩, hk⟩ := h
  unfold Items.EdgeBelow
  rw [hg]
  refine below_kept (fun j => j < s.items.size ∧ (∀ x ∈ t.verts, j ≠ vertItem x) ∧
      ∀ e ∈ t.edges, j ≠ edgeItem s.g e)
    (fun a c _ hp => ⟨hch a c hp, fun x hx hc => hvr x hx a (hc ▸ hp),
      fun e' he' hc => her e' he' a (hc ▸ hp)⟩)
    (fun j hj => (hk j hj.1 hj.2.1 hj.2.2).1) ⟨?_, fun x hx h => hv (by rw [vertItem_inj' h]; exact hx),
      fun e' _ h => vertItem_ne_edgeItem' hv' h⟩ _
  have := hT.size
  rw [hT.g_eq] at hv'
  show 1 + v < _
  omega

/-- The vertex entry pushed by `walkOutPre` spans only `V v`. -/
theorem TreeSite.L_spans {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s : WalkState} {e : Nat} {cls : OutClass} {y : Nat} {outs : List DfsOut} {L : List TEntry}
    {push : Bool} {done' : List (DfsOut × Bool)} {hv' : Bool} {bE' : List (Nat → Prop)}
    {sv' : List Nat} {sd' : List Bool} {sE : WalkState} {dir' : Bool} {L' : List TEntry}
    {push' : Bool} {D₃ : Array Bool}
    (H : TreeSite v d done rest hasVert base bE sv sd s e cls y outs L push done' hv' bE' sv' sd'
      sE dir' L' push' D₃) :
    ∀ t ∈ L, ∀ i ∈ t.spans.1 ++ t.spans.2, i = vertItem v := by
  intro t ht i hi
  rw [H.hL] at ht
  by_cases hp : push = true
  · rw [if_pos hp, List.mem_singleton] at ht
    subst ht
    rw [mem_setSides] at hi
    simpa using hi
  · rw [if_neg hp] at ht
    simp at ht

/-- The child's end push spans only `V y`. -/
theorem TreeSite.L'_spans {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s : WalkState} {e : Nat} {cls : OutClass} {y : Nat} {outs : List DfsOut} {L : List TEntry}
    {push : Bool} {done' : List (DfsOut × Bool)} {hv' : Bool} {bE' : List (Nat → Prop)}
    {sv' : List Nat} {sd' : List Bool} {sE : WalkState} {dir' : Bool} {L' : List TEntry}
    {push' : Bool} {D₃ : Array Bool}
    (H : TreeSite v d done rest hasVert base bE sv sd s e cls y outs L push done' hv' bE' sv' sd'
      sE dir' L' push' D₃) :
    ∀ t ∈ L', ∀ i ∈ t.spans.1 ++ t.spans.2, i = vertItem y := by
  intro t ht i hi
  rw [H.hL'] at ht
  by_cases hp : push' = true
  · rw [if_pos hp, List.mem_singleton] at ht
    subst ht
    rw [mem_setSides] at hi
    simpa using hi
  · rw [if_neg hp] at ht
    simp at ht

/-- The item frame of the child's walk at a tree site, up to its end push. -/
theorem TreeSite.frame_child {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s : WalkState} {e : Nat} {cls : OutClass} {y : Nat} {outs : List DfsOut} {L : List TEntry}
    {push : Bool} {done' : List (DfsOut × Bool)} {hv' : Bool} {bE' : List (Nat → Prop)}
    {sv' : List Nat} {sd' : List Bool} {sE : WalkState} {dir' : Bool} {L' : List TEntry}
    {push' : Bool} {D₃ : Array Bool}
    (H : TreeSite v d done rest hasVert base bE sv sd s e cls y outs L push done' hv' bE' sv' sd'
      sE dir' L' push' D₃) :
    OutFrame v (.tree e cls (.node y outs)) base s (pushEnd sE D₃ L') := by
  obtain ⟨top, htop, -⟩ := H.ctx.top
  obtain ⟨topc, htopc, -⟩ := H.ctx'.top
  exact outFrame_child htop htopc H.size_le H.gE H.items_kept (H.child_fresh topc htopc) H.L_spans
    H.L'_spans

/-- `finishBoundary` at a bridge tree out, over an abstract post-state: the child's outs are all
boundary (`bridge_bd`), so its end-of-outs `top` is empty and `hasVert` false; the end push `V y`
is the only entry popped; the items gain `Q e` under `V v`, and the fresh `I` item
`sE.items.size` and `V y` under `Q e`. Everything the parent's context sees is transported through
`items_kept` (`below_kept`) and the three new links. -/
theorem earCtx_bridge {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s : WalkState} {e : Nat} {cls : OutClass} {y : Nat} {outs : List DfsOut} {L : List TEntry}
    {push : Bool} {done' : List (DfsOut × Bool)} {hv' : Bool} {bE' : List (Nat → Prop)}
    {sv' : List Nat} {sd' : List Bool} {sE : WalkState} {dir' : Bool} {L' : List TEntry}
    {push' : Bool} {D₃ : Array Bool} {s' : WalkState}
    (H : TreeSite v d done rest hasVert base bE sv sd s e cls y outs L push done' hv' bE' sv' sd'
      sE dir' L' push' D₃)
    (hb : cls = .bridge)
    (hg : s'.g = s.g) (hts : s'.tstack = s.tstack)
    (hsv : ∀ k, k ≤ d → s'.stackVerts[k]! = s.stackVerts[k]!)
    (hsd : ∀ k, k < d → s'.stackDir[k]! = s.stackDir[k]!)
    (hsz' : s'.items.size = sE.items.size + 1)
    (hch : ∀ j, Items.ch s'.items j =
      if j = vertItem v then Items.ch sE.items (vertItem v) ++ [edgeItem s.g e]
      else if j = edgeItem s.g e then [sE.items.size, vertItem y]
      else if j = sE.items.size then [] else Items.ch sE.items j) :
    EarCtx v d (done ++ [(.tree e cls (.node y outs), false)]) rest false base bE sv sd s' := by
  have hC := H.ctx
  have hC' := H.ctx'
  have hge : d ≤ cls.lowval d := by rw [hb]; exact Nat.le_succ d
  have hhv : hasVert = false := by
    cases hhv : hasVert
    · rfl
    · exact absurd (hC.hv_ret hhv) (no_ret_before_boundary hge H.cls_ret H.rank)
  have hbd := H.bridge_bd hb
  have hv'f : hv' = false := by
    cases hh : hv'
    · rfl
    · obtain ⟨o', ho', hl⟩ := hC'.hv_ret hh
      exact absurd (hbd o' ho') (Nat.not_le.2 hl)
  have hpf : push = false := by
    cases hp : push
    · rfl
    · exact absurd (H.hpush.1 hp).2.1 (Nat.not_lt.2 hge)
  subst hhv hv'f hpf
  have hgE := H.gE
  have hv := H.hv
  have hy := H.hy
  have hsz := H.hsz
  have he := H.e_lt
  have hqf := hC.q_fresh _ (List.mem_cons_self ..) e (Or.inl rfl)
  have hyf := hC.v_fresh _ (List.mem_cons_self ..) e cls _ rfl y (List.mem_cons_self ..)
  have hedgesC : ∀ e', e' ∈ (DfsTree.node y outs).edges ↔ ∃ o' ∈ outs, subEdges o' e' :=
    fun e' => mem_subEdges_edgesList
  -- items untouched by the child's walk
  have hOldcl : ∀ a c, (a < s.items.size ∧ (∀ x ∈ (DfsTree.node y outs).verts, a ≠ vertItem x) ∧
      ∀ e' ∈ (DfsTree.node y outs).edges, a ≠ edgeItem s.g e') → Items.IsParent s.items a c →
      c < s.items.size ∧ (∀ x ∈ (DfsTree.node y outs).verts, c ≠ vertItem x) ∧
      ∀ e' ∈ (DfsTree.node y outs).edges, c ≠ edgeItem s.g e' := by
    intro a c _ h
    refine ⟨hC.ch_lt a c h, fun x hx hc => ?_, fun e' he' hc => ?_⟩
    · exact (hC.v_fresh _ (List.mem_cons_self ..) e cls _ rfl x hx).1 a (by rw [← hc]; exact h)
    · exact (hC.q_fresh _ (List.mem_cons_self ..) e' (Or.inr he')).1 a (by rw [← hc]; exact h)
  have hbelowE : ∀ a, (a < s.items.size ∧ (∀ x ∈ (DfsTree.node y outs).verts, a ≠ vertItem x) ∧
      ∀ e' ∈ (DfsTree.node y outs).edges, a ≠ edgeItem s.g e') → ∀ i,
      (Items.Below sE.items a i ↔ Items.Below s.items a i) :=
    fun a ha i => below_kept _ hOldcl (fun j hj => (H.items_kept j hj.1 hj.2.1 hj.2.2).1) ha i
  have hvlt : vertItem v < s.items.size := by show 1 + v < _; omega
  have hqlt : edgeItem s.g e < s.items.size := by show 1 + s.g.nv + e < _; omega
  have hvq : vertItem v ≠ edgeItem s.g e := vertItem_ne_edgeItem' hv
  have hOld_v : vertItem v < s.items.size ∧
      (∀ x ∈ (DfsTree.node y outs).verts, vertItem v ≠ vertItem x) ∧
      ∀ e' ∈ (DfsTree.node y outs).edges, vertItem v ≠ edgeItem s.g e' :=
    ⟨hvlt, fun x hx h => H.v_nc (by rw [vertItem_inj' h]; exact hx),
      fun _ _ => vertItem_ne_edgeItem' hv⟩
  have hOld_q : edgeItem s.g e < s.items.size ∧
      (∀ x ∈ (DfsTree.node y outs).verts, edgeItem s.g e ≠ vertItem x) ∧
      ∀ e' ∈ (DfsTree.node y outs).edges, edgeItem s.g e ≠ edgeItem s.g e' :=
    ⟨hqlt, fun x hx h => vertItem_ne_edgeItem' (H.verts_lt x hx) h.symm,
      fun e' he' h => H.e_ne (by rw [edgeItem_inj h]; exact he')⟩
  have hOld_span : ∀ t ∈ s.tstack, ∀ a ∈ t.spans.1 ++ t.spans.2, a < s.items.size ∧
      (∀ x ∈ (DfsTree.node y outs).verts, a ≠ vertItem x) ∧
      ∀ e' ∈ (DfsTree.node y outs).edges, a ≠ edgeItem s.g e' := fun t ht a ha =>
    ⟨hC.span_lt t ht a ha,
     fun x hx h => (hC.v_fresh _ (List.mem_cons_self ..) e cls _ rfl x hx).2.2.1 t ht
       (by rw [← h]; exact ha),
     fun e' he' h => (hC.q_fresh _ (List.mem_cons_self ..) e' (Or.inr he')).2.2 t ht
       (by rw [← h]; exact ha)⟩
  have hchq : Items.ch sE.items (edgeItem s.g e) = [] :=
    (H.items_kept _ hqlt hOld_q.2.1 hOld_q.2.2).1.trans hqf.2.1
  have hchn : Items.ch sE.items sE.items.size = [] := by simp [Items.ch]
  have hsle := H.size_le
  have hvn : vertItem v ≠ sE.items.size := by show 1 + v ≠ _; omega
  have hqn : edgeItem s.g e ≠ sE.items.size := by show 1 + s.g.nv + e ≠ _; omega
  have hyn : vertItem y ≠ sE.items.size := by show 1 + y ≠ _; omega
  have hyq : vertItem y ≠ edgeItem s.g e := vertItem_ne_edgeItem' hy
  have hyv : vertItem y ≠ vertItem v := fun h =>
    H.v_nc (by rw [← vertItem_inj' h]; exact List.mem_cons_self ..)
  have hP : ∀ p c, Items.IsParent s'.items p c ↔
      Items.IsParent sE.items p c ∨ (p = vertItem v ∧ c = edgeItem s.g e) ∨
        (p = edgeItem s.g e ∧ (c = sE.items.size ∨ c = vertItem y)) := by
    intro p c
    simp only [Items.IsParent, hch]
    by_cases h1 : p = vertItem v
    · subst h1; simp [hvq]
    · by_cases h2 : p = edgeItem s.g e
      · subst h2; simp [h1, hchq]
      · by_cases h3 : p = sE.items.size
        · subst h3; simp [h1, h2, hchn]
        · simp [h1, h2, h3]
  have hmono : ∀ {a i}, Items.Below sE.items a i → Items.Below s'.items a i := by
    intro a i h
    induction h with
    | refl => exact .refl
    | tail _ hs ih => exact ih.tail ((hP _ _).2 (Or.inl hs))
  have hbelow : ∀ a i, ¬ Items.Below sE.items a (vertItem v) →
      ¬ Items.Below sE.items a (edgeItem s.g e) →
      (Items.Below s'.items a i ↔ Items.Below sE.items a i) := by
    intro a i hav haq
    refine ⟨fun h => ?_, hmono⟩
    induction h with
    | refl => exact .refl
    | tail _ hji ih =>
      rcases (hP _ _).1 hji with h | ⟨rfl, rfl⟩ | ⟨rfl, -⟩
      · exact ih.tail h
      · exact absurd ih hav
      · exact absurd ih haq
  have hbelow_v : ∀ i, Items.Below s'.items (vertItem v) i ↔
      Items.Below sE.items (vertItem v) i ∨ i = edgeItem s.g e ∨ i = sE.items.size ∨
        Items.Below sE.items (vertItem y) i := by
    intro i
    constructor
    · intro h
      induction h with
      | refl => exact .inl .refl
      | tail _ hji ih =>
        rcases (hP _ _).1 hji with h | ⟨rfl, rfl⟩ | ⟨rfl, rfl | rfl⟩
        · rcases ih with ih | rfl | rfl | ih
          · exact .inl (ih.tail h)
          · exact absurd h (by rw [Items.IsParent, hchq]; simp)
          · exact absurd h (by rw [Items.IsParent, hchn]; simp)
          · exact .inr (.inr (.inr (ih.tail h)))
        · exact .inr (.inl rfl)
        · exact .inr (.inr (.inl rfl))
        · exact .inr (.inr (.inr .refl))
    · rintro (h | rfl | rfl | h)
      · exact hmono h
      · exact .single ((hP _ _).2 (.inr (.inl ⟨rfl, rfl⟩)))
      · exact (Relation.ReflTransGen.single ((hP _ _).2 (.inr (.inl ⟨rfl, rfl⟩)))).tail
          ((hP _ _).2 (.inr (.inr ⟨rfl, .inl rfl⟩)))
      · exact ((Relation.ReflTransGen.single ((hP _ _).2 (.inr (.inl ⟨rfl, rfl⟩)))).tail
          ((hP _ _).2 (.inr (.inr ⟨rfl, .inr rfl⟩)))).trans (hmono h)
  have hEBy : ∀ e', e' < s.g.ne →
      (Items.EdgeBelow s.g sE.items (vertItem y) e' ↔ e' ∈ (DfsTree.node y outs).edges) := by
    intro e' he'
    have h := hC'.vert_edges e' (by rw [hgE]; exact he')
    rw [hgE] at h
    rw [h, hedgesC]
    constructor
    · rintro ⟨o', ho', -, hs⟩
      exact ⟨o'.1, by rw [← H.hdone']; exact List.mem_map_of_mem ho', hs⟩
    · rintro ⟨o', ho', hs⟩
      rw [← H.hdone'] at ho'
      obtain ⟨p, hp, rfl⟩ := List.mem_map.1 ho'
      exact ⟨p, hp, hbd p hp, hs⟩
  have hEB_v : ∀ e', e' < s.g.ne → (Items.EdgeBelow s.g s'.items (vertItem v) e' ↔
      Items.EdgeBelow s.g s.items (vertItem v) e' ∨ subEdges (.tree e cls (.node y outs)) e') := by
    intro e' he'
    have hne : edgeItem s.g e' ≠ sE.items.size := by show 1 + s.g.nv + e' ≠ _; omega
    rw [Items.EdgeBelow, hbelow_v, hbelowE _ hOld_v]
    constructor
    · rintro (h | h | h | h)
      · exact .inl h
      · exact .inr (Or.inl (edgeItem_inj h))
      · exact absurd h hne
      · exact .inr (Or.inr ((hEBy e' he').1 h))
    · rintro (h | h)
      · exact .inl h
      rcases (show e' = e ∨ e' ∈ (DfsTree.node y outs).edges from h) with rfl | h
      · exact .inr (.inl rfl)
      · exact .inr (.inr (.inr ((hEBy e' he').2 h)))
  have hspanB : ∀ t ∈ s.tstack, ∀ i ∈ t.spans.1 ++ t.spans.2, ∀ j,
      (Items.Below s'.items i j ↔ Items.Below s.items i j) := by
    intro t ht i hi j
    rw [← hbelowE _ (hOld_span t ht i hi)]
    apply hbelow
    · intro h
      have := Items.Below.eq_of_no_parent H.v_root h
      subst this
      exact Bool.false_ne_true (hC.vert_free t ht hi)
    · intro h
      have := Items.Below.eq_of_no_parent H.q_root h
      subst this
      exact hqf.2.2 t ht hi
  have hedges : ∀ t ∈ s.tstack, ∀ e',
      (TEntry.edges s.g s'.items t e' ↔ TEntry.edges s.g s.items t e') := fun t ht e' =>
    TEntry.edges_congr (fun i hi e' => hspanB t ht i hi _) e'
  have hAV : afterVert (done ++ [(.tree e cls (.node y outs), false)]) = afterVert done := by
    simp [afterVert]
  have hVL : ∀ x, x ∈ DfsOut.vertsList (done.map (·.1)) →
      x ∈ DfsOut.vertsList ((done ++ [((DfsOut.tree e cls (.node y outs), false) :
        DfsOut × Bool)]).map (·.1)) := by
    intro x hx
    rw [DfsOut.vertsList_eq] at hx ⊢
    rw [List.map_append, List.flatMap_append]
    exact List.mem_append_left _ hx
  obtain ⟨top, htop, hCT⟩ := hC.top
  have hmemT : ∀ t ∈ top, t ∈ s.tstack := fun t ht => by
    rw [htop]; exact List.mem_append_left _ ht
  have hbase_mem : ∀ k, k < base.length → base[k]! ∈ s.tstack := fun k hk => by
    rw [htop, getElem!_pos base k hk]; exact List.mem_append_right _ (List.getElem_mem hk)
  obtain ⟨above, below, hsp, -, -, -, -, -, hB⟩ := hCT.split
  simp only [Bool.false_eq_true, ↓reduceIte] at hsp
  obtain ⟨rfl, rfl⟩ := hsp
  have hbook := hC.vert_book rfl
  exact
    { top := ⟨top, by rw [hts, htop],
        { split := ⟨[], top, by simp, fun t ht => (nomatch ht), fun t ht => (nomatch ht),
            List.Pairwise.nil, fun t ht => (nomatch ht), fun t ht => (nomatch ht), hB⟩
          bot := fun t ht k hk => by rw [hsv k (by omega)]; exact hCT.bot t ht k hk
          vitems := fun t ht x hx hm =>
            (hCT.vitems t ht x (by rwa [hg] at hx) hm).imp id (hVL x)
          qitems := fun t ht e' he' hm =>
            let ⟨o, ho, hs⟩ := hCT.qitems t ht e' (by rwa [hg] at he') (by rwa [hg] at hm)
            ⟨o, List.mem_append_left _ ho, hs⟩
          edges := fun t ht e' he' hte => by
            rw [hg] at he' hte
            obtain ⟨o, ho, hs⟩ := hCT.edges t ht e' he' ((hedges t (hmemT t ht) e').1 hte)
            exact ⟨o, List.mem_append_left _ ho, hs⟩
          cover := fun o ho hlt e' he' hs => by
            rw [hg] at he' ⊢
            rcases List.mem_append.1 ho with h | h
            · obtain ⟨t, ht, hte⟩ := hCT.cover o h hlt e' he' hs
              exact ⟨t, ht, (hedges t (hmemT t ht) e').2 hte⟩
            · rw [List.mem_singleton] at h
              subst h
              exact absurd hlt (Nat.not_lt.2 hge)
          ret := fun h => (hCT.ret h).imp fun o ⟨ho, hl⟩ => ⟨List.mem_append_left _ ho, hl⟩ }⟩
      base_edges := ⟨hC.base_edges.1, fun k hk e' he' => by
        rw [hg] at he' ⊢
        exact (hedges _ (hbase_mem k hk) e').trans (hC.base_edges.2 k hk e' he')⟩
      base_bot := hC.base_bot
      noVert_after := fun _ => by
        have h := hC.noVert_after rfl
        simp only [afterVert, List.map_eq_nil_iff] at h
        simp [afterVert, List.filter_append, h]
      vfirst := fun t ht hvt _ => absurd hvt (by
        rw [hts, htop] at ht
        rcases List.mem_append.1 ht with h | h
        · exact hB t h
        · exact hC.base_bot t h)
      sv := fun k hk => by rw [hsv k hk]; exact hC.sv k hk
      sd := fun k hk => by rw [hsd k hk]; exact hC.sd k hk
      sv_d := by rw [hsv d (Nat.le_refl _)]; exact hC.sv_d
      path := fun k k' h h' => by rw [hsv k (by omega), hsv k' h']; exact hC.path k k' h h'
      v_root := fun p h => by
        rcases (hP p _).1 h with h | ⟨_, h⟩ | ⟨_, h | h⟩
        · exact H.v_root p h
        · exact hvq h
        · exact hvn h
        · exact hyv h.symm
      vert_free := fun t ht hm => by rw [hts] at ht; exact hC.vert_free t ht hm
      afterVert_ret := fun o ho => by rw [hAV] at ho; exact hC.afterVert_ret o ho
      hv_ret := fun h => (nomatch h)
      vert_book := fun _ => by
        rw [hg]
        have hadjO : ∀ u w, s.g.AdjIn (Items.EdgeBelow s.g s.items (vertItem v)) u w →
            s.g.AdjIn (Items.EdgeBelow s.g s'.items (vertItem v)) u w :=
          fun u w ⟨e₁, h1, h2, h3⟩ => ⟨e₁, h1, (hEB_v e₁ h1).2 (.inl h2), h3⟩
        have hadjC : ∀ u w, s.g.AdjIn (Items.EdgeBelow s.g sE.items (vertItem y)) u w →
            s.g.AdjIn (Items.EdgeBelow s.g s'.items (vertItem v)) u w :=
          fun u w ⟨e₁, h1, h2, h3⟩ =>
            ⟨e₁, h1, (hEB_v e₁ h1).2 (.inr (Or.inr ((hEBy e₁ h1).1 h2))), h3⟩
        have hsymm : ∀ a b, Relation.ReflTransGen
            (s.g.AdjIn (Items.EdgeBelow s.g s'.items (vertItem v))) a b →
            Relation.ReflTransGen (s.g.AdjIn (Items.EdgeBelow s.g s'.items (vertItem v))) b a := by
          intro a b h
          induction h with
          | refl => exact .refl
          | tail _ hs ih => exact (Relation.ReflTransGen.single (Graph.AdjIn.symm hs)).trans ih
        have hEe : Items.EdgeBelow s.g s'.items (vertItem v) e := (hEB_v e he).2 (.inr (Or.inl rfl))
        have hyv_step : Relation.ReflTransGen
            (s.g.AdjIn (Items.EdgeBelow s.g s'.items (vertItem v))) v y :=
          .single ⟨e, he, hEe, by
            rcases H.e_ends with h | h
            · exact Or.inr (by rw [← h])
            · exact Or.inl (Prod.ext (Prod.mk.inj h).2 (Prod.mk.inj h).1)⟩
        have hreachO : ∀ e₁, e₁ < s.g.ne → Items.EdgeBelow s.g s.items (vertItem v) e₁ →
            Relation.ReflTransGen (s.g.AdjIn (Items.EdgeBelow s.g s'.items (vertItem v))) v
              (s.g.edges[e₁]!).1 := by
          intro e₁ h1 hE1
          obtain ⟨o, ho, hlo, -⟩ := (hC.vert_edges e₁ h1).1 hE1
          have he0 : Items.EdgeBelow s.g s.items (vertItem v) o.1.e :=
            (hC.vert_edges _ (H.inc o ho).1).2 ⟨o, ho, hlo, subEdges_e _⟩
          have step : Relation.ReflTransGen (s.g.AdjIn (Items.EdgeBelow s.g s'.items (vertItem v)))
              v (s.g.edges[o.1.e]!).1 := by
            rcases (H.inc o ho).2 with h | h
            · rw [h]
            · exact .single ⟨o.1.e, (H.inc o ho).1, (hEB_v _ (H.inc o ho).1).2 (.inl he0),
                .inr (by rw [h])⟩
          exact step.trans (Relation.ReflTransGen.mono hadjO _ _
            (hbook.1 _ _ (H.inc o ho).1 h1 he0 hE1))
        have hreachC : ∀ e₁, e₁ < s.g.ne → e₁ ∈ (DfsTree.node y outs).edges →
            Relation.ReflTransGen (s.g.AdjIn (Items.EdgeBelow s.g s'.items (vertItem v))) v
              (s.g.edges[e₁]!).1 := by
          intro e₁ h1 hm
          obtain ⟨o', ho', hs⟩ := (hedgesC e₁).1 hm
          rw [← H.hdone'] at ho'
          obtain ⟨p, hp, rfl⟩ := List.mem_map.1 ho'
          have hpe := H.inc' p hp
          have hE0 : Items.EdgeBelow s.g sE.items (vertItem y) p.1.e :=
            (hEBy _ hpe.1).2 ((hedgesC _).2
              ⟨p.1, by rw [← H.hdone']; exact List.mem_map_of_mem hp, subEdges_e _⟩)
          have hE1 : Items.EdgeBelow s.g sE.items (vertItem y) e₁ := (hEBy _ h1).2 hm
          have step : Relation.ReflTransGen (s.g.AdjIn (Items.EdgeBelow s.g s'.items (vertItem v)))
              y (s.g.edges[p.1.e]!).1 := by
            rcases hpe.2 with h | h
            · rw [h]
            · exact .single ⟨p.1.e, hpe.1,
                (hEB_v _ hpe.1).2 (.inr (Or.inr ((hEBy _ hpe.1).1 hE0))), .inr (by rw [h])⟩
          have hbookC := hC'.vert_book rfl
          rw [hgE] at hbookC
          exact hyv_step.trans (step.trans (Relation.ReflTransGen.mono hadjC _ _
            (hbookC.1 _ _ hpe.1 h1 hE0 hE1)))
        have hreach : ∀ e₁, e₁ < s.g.ne → Items.EdgeBelow s.g s'.items (vertItem v) e₁ →
            Relation.ReflTransGen (s.g.AdjIn (Items.EdgeBelow s.g s'.items (vertItem v))) v
              (s.g.edges[e₁]!).1 := by
          intro e₁ h1 hE1
          rcases (hEB_v e₁ h1).1 hE1 with h | h
          · exact hreachO e₁ h1 h
          rcases (show e₁ = e ∨ e₁ ∈ (DfsTree.node y outs).edges from h) with rfl | h
          · rcases H.e_ends with h | h
            · rw [← h]; exact hyv_step
            · rw [(Prod.mk.inj h).2.symm]
          · exact hreachC e₁ h1 h
        refine ⟨fun e₁ e₂ h1 h2 hE1 hE2 => (hsymm _ _ (hreach e₁ h1 hE1)).trans (hreach e₂ h2 hE2),
          fun x e₁ e₂ h1 h2 hE1 hnE2 hi1 hi2 => ?_⟩
        rcases (hEB_v e₁ h1).1 hE1 with h | h
        · exact hbook.2 x e₁ e₂ h1 h2 h (fun h' => hnE2 ((hEB_v e₂ h2).2 (.inl h'))) hi1 hi2
        · rcases H.ends hge e₁ h h1 x hi1 with hx | hx
          · exact .inl hx
          · exact absurd ((hEB_v e₂ h2).2 (.inr (H.comp e₂ h2 x hi2 hx))) hnE2
      vert_disj := fun _ t ht e' he' hte hE => by
        rw [hts] at ht
        rw [hg] at he' hte hE
        rcases (hEB_v e' he').1 hE with h | h
        · exact hC.vert_disj rfl t ht e' he' ((hedges t ht e').1 hte) h
        · obtain ⟨i, hi, hb'⟩ := (hedges t ht e').1 hte
          have hqf' := hC.q_fresh _ (List.mem_cons_self ..) e' h
          exact hqf'.2.2 t ht (Items.Below.eq_of_no_parent hqf'.1 hb' ▸ hi)
      vert_edges := fun e' he' => by
        rw [hg] at he' ⊢
        rw [hEB_v e' he', hC.vert_edges e' he']
        constructor
        · rintro (⟨o, ho, hlo, hs⟩ | h)
          · exact ⟨o, List.mem_append_left _ ho, hlo, hs⟩
          · exact ⟨(.tree e cls (.node y outs), false),
              List.mem_append_right _ (List.mem_singleton_self _), hge, h⟩
        · rintro ⟨o, ho, hlo, hs⟩
          rcases List.mem_append.1 ho with h | h
          · exact .inl ⟨o, h, hlo, hs⟩
          · rw [List.mem_singleton] at h
            subst h
            exact .inr hs
      touch_bot := fun t ht ⟨e', he', hte⟩ => by
        rw [hts] at ht
        rw [hg] at he' hte ⊢
        obtain ⟨e₀, h0, hE0, hi0⟩ := hC.touch_bot t ht ⟨e', he', (hedges t ht e').1 hte⟩
        exact ⟨e₀, h0, (hedges t ht e₀).2 hE0, hi0⟩
      span_root := fun t ht i hi p h => by
        rw [hts] at ht
        rcases (hP p i).1 h with h | ⟨_, rfl⟩ | ⟨_, rfl | rfl⟩
        · have hO := hOld_span t ht i hi
          exact hC.span_root t ht i hi p (((H.items_kept i hO.1 hO.2.1 hO.2.2).2 p).1 h)
        · exact hqf.2.2 t ht hi
        · exact Nat.lt_irrefl _ (Nat.lt_of_lt_of_le (hC.span_lt t ht _ hi) hsle)
        · exact hyf.2.2.1 t ht hi
      span_lt := fun t ht i hi => by
        rw [hts] at ht
        rw [hsz']
        exact Nat.lt_succ_of_lt (Nat.lt_of_lt_of_le (hC.span_lt t ht i hi) hsle)
      ch_lt := fun p c h => by
        rw [hsz']
        rcases (hP p c).1 h with h | ⟨_, rfl⟩ | ⟨_, rfl | rfl⟩
        · exact Nat.lt_succ_of_lt (hC'.ch_lt p c h)
        · exact Nat.lt_succ_of_lt (Nat.lt_of_lt_of_le hqlt hsle)
        · exact Nat.lt_succ_self _
        · exact Nat.lt_succ_of_lt (Nat.lt_of_lt_of_le (by show 1 + y < _; omega) hsle)
      disj := by
        rw [hts]
        refine List.Pairwise.imp_of_mem (fun {t t'} ht ht' h e' he' hte hte' => ?_) hC.disj
        rw [hg] at he' hte hte'
        exact h e' he' ((hedges t ht e').1 hte) ((hedges t' ht' e').1 hte')
      span_disj := by rw [hts]; exact hC.span_disj
      q_fresh := fun o ho e' hs => by
        obtain ⟨hns, hlt⟩ := H.rest_nd o ho e' hs
        obtain ⟨hr, hc, hm⟩ := hC.q_fresh o (List.mem_cons_of_mem _ ho) e' hs
        rw [hg]
        have hOld : edgeItem s.g e' < s.items.size ∧
            (∀ x ∈ (DfsTree.node y outs).verts, edgeItem s.g e' ≠ vertItem x) ∧
            ∀ e'' ∈ (DfsTree.node y outs).edges, edgeItem s.g e' ≠ edgeItem s.g e'' :=
          ⟨by show 1 + s.g.nv + e' < _; omega,
           fun x hx h => vertItem_ne_edgeItem' (H.verts_lt x hx) h.symm,
           fun e'' he'' h => hns (Or.inr (by rw [edgeItem_inj h]; exact he''))⟩
        have hK := H.items_kept _ hOld.1 hOld.2.1 hOld.2.2
        have hq' : edgeItem s.g e' ≠ edgeItem s.g e := fun h => hns (Or.inl (edgeItem_inj h))
        have hn' : edgeItem s.g e' ≠ sE.items.size := by show 1 + s.g.nv + e' ≠ _; omega
        have hv' : edgeItem s.g e' ≠ vertItem v := (vertItem_ne_edgeItem' hv).symm
        have hy' : edgeItem s.g e' ≠ vertItem y := (vertItem_ne_edgeItem' hy).symm
        refine ⟨fun p h => ?_, ?_, fun t ht => ?_⟩
        · rcases (hP p _).1 h with h | ⟨_, h⟩ | ⟨_, h | h⟩
          · exact hr p ((hK.2 p).1 h)
          · exact hq' h
          · exact hn' h
          · exact hy' h
        · rw [hch, if_neg hv', if_neg hq', if_neg hn', hK.1]; exact hc
        · rw [hts] at ht; exact hm t ht
      v_fresh := fun o ho e₁ cls₁ child ho₁ w hw => by
        obtain ⟨hwC, hwv⟩ := H.rest_nv o ho e₁ cls₁ child ho₁ w hw
        obtain ⟨hr, hc, hm, hsvk, hvs, htch⟩ :=
          hC.v_fresh o (List.mem_cons_of_mem _ ho) e₁ cls₁ child ho₁ w hw
        have hOld : vertItem w < s.items.size ∧
            (∀ x ∈ (DfsTree.node y outs).verts, vertItem w ≠ vertItem x) ∧
            ∀ e'' ∈ (DfsTree.node y outs).edges, vertItem w ≠ edgeItem s.g e'' :=
          ⟨by show 1 + w < _; omega, fun x hx h => hwC (by rw [vertItem_inj' h]; exact hx),
           fun _ _ => vertItem_ne_edgeItem' hwv⟩
        have hK := H.items_kept _ hOld.1 hOld.2.1 hOld.2.2
        have hwv' : vertItem w ≠ vertItem v := fun h =>
          hsvk d (Nat.le_refl _) (hC.sv_d.trans (vertItem_inj' h).symm)
        have hwq : vertItem w ≠ edgeItem s.g e := vertItem_ne_edgeItem' hwv
        have hwn : vertItem w ≠ sE.items.size := by show 1 + w ≠ _; omega
        have hwy : vertItem w ≠ vertItem y := fun h =>
          hwC (by rw [vertItem_inj' h]; exact List.mem_cons_self ..)
        refine ⟨fun p h => ?_, ?_, fun t ht => ?_, fun k hk => ?_, fun t ht => ?_, fun t ht => ?_⟩
        · rcases (hP p _).1 h with h | ⟨_, h⟩ | ⟨_, h | h⟩
          · exact hr p ((hK.2 p).1 h)
          · exact hwq h
          · exact hwn h
          · exact hwy h
        · rw [hch, if_neg hwv', if_neg hwq, if_neg hwn, hK.1]; exact hc
        · rw [hts] at ht; exact hm t ht
        · rw [hsv k hk]; exact hsvk k hk
        · rw [hts] at ht; exact hvs t ht
        · rw [hts] at ht
          rw [hg]
          rintro ⟨e'', he'', hte, hi⟩
          exact htch t ht ⟨e'', he'', (hedges t ht e'').1 hte, hi⟩ }

/-- The item frame of `finishBoundary` at a bridge: the child's frame, then `Q e` under `V v`, the
fresh `O` item and `V y` under `Q e`. -/
theorem outFrame_bridge {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s : WalkState} {e : Nat} {cls : OutClass} {y : Nat} {outs : List DfsOut} {L : List TEntry}
    {push : Bool} {done' : List (DfsOut × Bool)} {hv' : Bool} {bE' : List (Nat → Prop)}
    {sv' : List Nat} {sd' : List Bool} {sE : WalkState} {dir' : Bool} {L' : List TEntry}
    {push' : Bool} {D₃ : Array Bool} {s' : WalkState}
    (H : TreeSite v d done rest hasVert base bE sv sd s e cls y outs L push done' hv' bE' sv' sd'
      sE dir' L' push' D₃)
    (hb : cls = .bridge)
    (hg : s'.g = s.g) (hts : s'.tstack = s.tstack)
    (hsz' : s'.items.size = sE.items.size + 1)
    (hch : ∀ j, Items.ch s'.items j =
      if j = vertItem v then Items.ch sE.items (vertItem v) ++ [edgeItem s.g e]
      else if j = edgeItem s.g e then [sE.items.size, vertItem y]
      else if j = sE.items.size then [] else Items.ch sE.items j) :
    OutFrame v (.tree e cls (.node y outs)) base s s' := by
  have hC := H.ctx
  have hC' := H.ctx'
  have hge : d ≤ cls.lowval d := by rw [hb]; exact Nat.le_succ d
  have hbd := H.bridge_bd hb
  have hpf : push = false := by
    cases hp : push
    · rfl
    · exact absurd (H.hpush.1 hp).2.1 (Nat.not_lt.2 hge)
  have hL : L = [] := by rw [H.hL, hpf]; rfl
  have hsEts : sE.tstack = s.tstack := by
    obtain ⟨top', htop', hCT'⟩ := hC'.top
    rw [hL, List.nil_append] at htop'
    cases top' with
    | nil => simpa using htop'
    | cons t tl =>
      obtain ⟨o', ho', hl⟩ := hCT'.ret (List.cons_ne_nil _ _)
      exact absurd (hbd o' ho') (Nat.not_le.2 hl)
  have hsz := H.hsz
  have hsle := H.size_le
  have hv := H.hv
  have he := H.e_lt
  have hgE := H.gE
  have hqlt : edgeItem s.g e < s.items.size := by show 1 + s.g.nv + e < _; omega
  have hqE : Items.ch sE.items (edgeItem s.g e) = [] := by
    rw [(H.items_kept _ hqlt (fun x hx => (vertItem_ne_edgeItem' (H.verts_lt x hx)).symm)
      (fun e' he' h => H.e_ne ((edgeItem_inj h) ▸ he'))).1]
    exact (hC.q_fresh _ (List.mem_cons_self ..) e (subEdges_e _)).2.1
  have hvq : vertItem v ≠ edgeItem s.g e := vertItem_ne_edgeItem' hv
  have hjq' : ∀ j, (∀ e' ∈ (DfsOut.tree e cls (.node y outs)).edges,
      j ≠ edgeItem (pushEnd sE D₃ L').g e') → j ≠ edgeItem s.g e := fun j hje => by
    have := hje e (by simp [DfsOut.edges])
    rwa [show (pushEnd sE D₃ L').g = s.g from hgE] at this
  have h1 : sE.items.size ≠ vertItem v := by show _ ≠ 1 + v; omega
  have h2 : sE.items.size ≠ edgeItem s.g e := by show _ ≠ 1 + s.g.nv + e; omega
  refine H.frame_child.trans (outFrame_of_mod
    (fun p => p = vertItem v ∨ p = edgeItem s.g e ∨ p = sE.items.size)
    (by simp only [pushEnd]; omega) (hg.trans hgE.symm) ?_ ?_ ?_ ?_)
  · intro p hp
    rw [hch p]
    simp only [not_or] at hp
    simp [pushEnd, hp.1, hp.2.1, hp.2.2]
  · intro p j hp hj
    obtain ⟨hjlt, hjv, hxv, hje, -⟩ := hj
    have hjq : j ≠ edgeItem s.g e := hjq' j hje
    have hjy : j ≠ vertItem y := hxv y (by simp [DfsOut.verts, DfsTree.verts])
    have hjlt' : j < sE.items.size := hjlt
    simp only [Items.IsParent, pushEnd]
    rw [hch p]
    rcases hp with rfl | rfl | rfl
    · simp [hjq]
    · simp [hvq.symm, hqE, Nat.ne_of_lt hjlt', hjy]
    · simp [Items.ch_of_le sE.items _ (Nat.le_refl _), h1, h2]
  · intro j hj
    obtain ⟨hjlt, hjv, -, hje, -⟩ := hj
    simp only [not_or]
    exact ⟨hjv, hjq' j hje, Nat.ne_of_lt hjlt⟩
  · intro top' htop' t ht
    rw [hts] at htop'
    exact ⟨L' ++ top', by simp [pushEnd, hsEts, htop'], List.mem_append_right _ ht⟩

/-- The tree-edge step for a bridge (`lowval = d + 1`): the child's outs are all boundary, so its
end push `V y` is the only entry popped; `earCtx_bridge` over the concrete post-state. -/
theorem ctx_step_tree_bridge {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s : WalkState} {e : Nat} {cls : OutClass} {y : Nat} {outs : List DfsOut} {L : List TEntry}
    {push : Bool} {done' : List (DfsOut × Bool)} {hv' : Bool} {bE' : List (Nat → Prop)}
    {sv' : List Nat} {sd' : List Bool} {sE : WalkState} {dir' : Bool} {L' : List TEntry}
    {push' : Bool} {D₃ : Array Bool}
    (H : TreeSite v d done rest hasVert base bE sv sd s e cls y outs L push done' hv' bE' sv' sd'
      sE dir' L' push' D₃)
    (hb : cls = .bridge) :
    wp (finishEdge v d (.tree e cls (.node y outs)) (L ++ s.tstack).length (hasVert || push))
      (fun hv'' s' => EarCtx v d (done ++ [(.tree e cls (.node y outs), hasVert || push)]) rest hv''
        base bE sv sd s' ∧ OutFrame v (.tree e cls (.node y outs)) base s s' ∧
        (CtxShape v d done hasVert base s →
          CtxShape v d (done ++ [(.tree e cls (.node y outs), hasVert || push)]) hv'' base s'))
      (pushEnd sE D₃ L') := by
  have hC := H.ctx
  have hC' := H.ctx'
  have hge : d ≤ cls.lowval d := by rw [hb]; exact Nat.le_succ d
  have hhv : hasVert = false := by
    cases hhv : hasVert
    · rfl
    · exact absurd (hC.hv_ret hhv) (no_ret_before_boundary hge H.cls_ret H.rank)
  have hbd := H.bridge_bd hb
  have hv'f : hv' = false := by
    cases hh : hv'
    · rfl
    · obtain ⟨o', ho', hl⟩ := hC'.hv_ret hh
      exact absurd (hbd o' ho') (Nat.not_le.2 hl)
  have hpf : push = false := by
    cases hp : push
    · rfl
    · exact absurd (H.hpush.1 hp).2.1 (Nat.not_lt.2 hge)
  have hp' : push' = true := H.hpush'.2 hv'f
  have hdir : dir' = true := H.hdir' hp'
  have hL : L = [] := by rw [H.hL, hpf]; rfl
  have hL' : L' = [⟨y, d + 1, sE.nxtEdgeIdx, setSides true [vertItem y] []⟩] := by
    rw [H.hL', hp', hdir]; rfl
  have hsEts : sE.tstack = s.tstack := by
    obtain ⟨top', htop', hCT'⟩ := hC'.top
    rw [hL, List.nil_append] at htop'
    cases top' with
    | nil => simpa using htop'
    | cons t tl =>
      obtain ⟨o', ho', hl⟩ := hCT'.ret (List.cons_ne_nil _ _)
      exact absurd (hbd o' ho') (Nat.not_le.2 hl)
  have hsz := H.hsz
  have hsle := H.size_le
  have hv := H.hv
  have he := H.e_lt
  have hvlt : vertItem v < sE.items.size := by show 1 + v < _; omega
  have hqlt : edgeItem s.g e < sE.items.size := by show 1 + s.g.nv + e < _; omega
  have hvq : vertItem v ≠ edgeItem s.g e := vertItem_ne_edgeItem' H.hv
  have hgE := H.gE
  have hsdk : ∀ k, k < d → D₃[k]! = s.stackDir[k]! := fun k hk => by
    rw [H.sd₃ k (Nat.le_of_lt hk), H.sdlo k (Nat.le_of_lt hk)]
    exact getElem!_set!_ne' _ _ _ _ (by omega)
  subst hhv hpf
  simp only [Bool.or_self]
  have hge' : cls.lowval d ≥ d := hge
  rw [finishEdge_eq]
  simp only [finishEdge', wp_bind, wp_get, wp_stackDir, DfsOut.cls, hge', ↓reduceIte]
  unfold finishBoundary
  subst hb
  simp only [wp_bind, wp_modifyItem, wp_modify, wp_allocItem, wp_makeVs, wp_popTstack, wp_pure,
    DfsOut.cls, DfsOut.e, DfsOut.dest, OutClass.isTree, OutClass.lowval, beq_self_eq_true,
    ↓reduceIte, pushEnd, hL', hgE, List.singleton_append, List.head!_cons, List.tail_cons, setSides]
  refine (fun hch => ⟨earCtx_bridge H rfl rfl hsEts (fun k hk => H.svlo k hk) hsdk (by simp) hch,
    outFrame_bridge H rfl rfl hsEts (by simp) hch,
    fun hS => ctxShape_same hS hge' (by simp [hsEts]) hsdk⟩) ?_
  intro j
  simp only [Items.ch, Array.getElem?_modify, Array.getElem?_push, Array.size_modify]
  by_cases h1 : j = vertItem v
  · subst h1
    simp [hvq.symm, Array.getElem?_eq_getElem hvlt, Nat.ne_of_lt hvlt, Nat.ne_of_gt hvlt]
  · by_cases h2 : j = edgeItem s.g e
    · subst h2
      simp [h1, Ne.symm h1, Array.getElem?_eq_getElem hqlt, Nat.ne_of_lt hqlt, Nat.ne_of_gt hqlt]
    · by_cases h3 : j = sE.items.size
      · subst h3; simp [h1, h2, Ne.symm h1, Ne.symm h2]
      · simp [h1, h2, h3, Ne.symm h1, Ne.symm h2, Ne.symm h3]

/-- The parent context after a component tree out (`lowval = d`), from the post-state equations:
`sE` (the child's end-of-outs state, its stack `t₁ :: V y :: s.tstack`) with `t₁`'s side-1 items
and `V y` moved under `Q e` and `Q e` appended under `V v`; the two entries popped. -/
theorem earCtx_comp {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s : WalkState} {e : Nat} {cls : OutClass} {y : Nat} {outs : List DfsOut} {L : List TEntry}
    {push : Bool} {done' : List (DfsOut × Bool)} {hv' : Bool} {bE' : List (Nat → Prop)}
    {sv' : List Nat} {sd' : List Bool} {sE : WalkState} {dir' : Bool} {L' : List TEntry}
    {push' : Bool} {D₃ : Array Bool} {s' : WalkState}
    (H : TreeSite v d done rest hasVert base bE sv sd s e cls y outs L push done' hv' bE' sv' sd'
      sE dir' L' push' D₃)
    (hc : cls = .component) {t₁ : TEntry} {f₂ : Nat}
    (hsh : sE.tstack = t₁ :: ⟨y, d + 1, f₂, ([], [vertItem y])⟩ :: s.tstack)
    (ht₁ : t₁.spans.2 = [])
    (hg : s'.g = s.g) (hts : s'.tstack = s.tstack)
    (hsv : ∀ k, k ≤ d → s'.stackVerts[k]! = s.stackVerts[k]!)
    (hsd : ∀ k, k < d → s'.stackDir[k]! = s.stackDir[k]!)
    (hsz' : s'.items.size = sE.items.size)
    (hch : ∀ j, Items.ch s'.items j =
      if j = vertItem v then Items.ch sE.items (vertItem v) ++ [edgeItem s.g e]
      else if j = edgeItem s.g e then t₁.spans.1 ++ [vertItem y]
      else Items.ch sE.items j) :
    EarCtx v d (done ++ [(.tree e cls (.node y outs), false)]) rest false base bE sv sd s' := by
  have hC := H.ctx
  have hC' := H.ctx'
  have hge : d ≤ cls.lowval d := by rw [hc]; exact Nat.le_refl d
  have hhv : hasVert = false := by
    cases hhv : hasVert
    · rfl
    · exact absurd (hC.hv_ret hhv) (no_ret_before_boundary hge H.cls_ret H.rank)
  have hpf : push = false := by
    cases hp : push
    · rfl
    · exact absurd (H.hpush.1 hp).2.1 (Nat.not_lt.2 hge)
  subst hhv hpf
  have hL : L = [] := by rw [H.hL]; rfl
  have hgE := H.gE
  have hv := H.hv
  have hy := H.hy
  have hsz := H.hsz
  have he := H.e_lt
  have hqf := hC.q_fresh _ (List.mem_cons_self ..) e (Or.inl rfl)
  have hyf := hC.v_fresh _ (List.mem_cons_self ..) e cls _ rfl y (List.mem_cons_self ..)
  have hedgesC : ∀ e', e' ∈ (DfsTree.node y outs).edges ↔ ∃ o' ∈ outs, subEdges o' e' :=
    fun e' => mem_subEdges_edgesList
  -- the child's top: `[t₁, V y]`
  have ht1 : t₁ ∈ sE.tstack := by rw [hsh]; exact List.mem_cons_self ..
  have ht2 : (⟨y, d + 1, f₂, ([], [vertItem y])⟩ : TEntry) ∈ sE.tstack := by
    rw [hsh]; exact List.mem_cons_of_mem _ (List.mem_cons_self ..)
  obtain ⟨top', htop', hCT'⟩ := hC'.top
  rw [hL, List.nil_append] at htop'
  have htop'eq : top' = [t₁, ⟨y, d + 1, f₂, ([], [vertItem y])⟩] :=
    (List.append_inj_left' (show [t₁, _] ++ s.tstack = top' ++ s.tstack from
      hsh.symm.trans htop') rfl).symm
  have hmem1 : t₁ ∈ top' := by rw [htop'eq]; exact List.mem_cons_self ..
  have hsp1 : ∀ a ∈ t₁.spans.1, (∀ p, ¬ Items.IsParent sE.items p a) ∧ a < sE.items.size ∧
      ∀ t ∈ s.tstack, a ∉ t.spans.1 ++ t.spans.2 := by
    intro a ha
    have ha' : a ∈ t₁.spans.1 ++ t₁.spans.2 := List.mem_append_left _ ha
    refine ⟨hC'.span_root t₁ ht1 a ha', hC'.span_lt t₁ ht1 a ha', fun t ht => ?_⟩
    have hpw := hC'.span_disj
    rw [hsh] at hpw
    exact (List.pairwise_cons.1 hpw).1 t (List.mem_cons_of_mem _ ht) a ha'
  have hvy_lt : vertItem y < sE.items.size := hC'.span_lt _ ht2 (vertItem y) (by simp)
  have hv1 : ∀ x, x < s.g.nv → vertItem x ∈ t₁.spans.1 → x = y ∨ x ∈ DfsOut.vertsList outs := by
    intro x hx h
    have := hCT'.vitems t₁ hmem1 x (by rw [hgE]; exact hx) (List.mem_append_left _ h)
    rwa [H.hdone'] at this
  have hq1 : ∀ e', e' < s.g.ne → edgeItem s.g e' ∈ t₁.spans.1 → e' ∈ (DfsTree.node y outs).edges := by
    intro e' he' h
    obtain ⟨o', ho', hs⟩ := hCT'.qitems t₁ hmem1 e' (by rw [hgE]; exact he')
      (by rw [hgE]; exact List.mem_append_left _ h)
    exact (hedgesC e').2 ⟨o'.1, by rw [← H.hdone']; exact List.mem_map_of_mem ho', hs⟩
  -- items untouched by the child's walk
  have hOldcl : ∀ a c, (a < s.items.size ∧ (∀ x ∈ (DfsTree.node y outs).verts, a ≠ vertItem x) ∧
      ∀ e' ∈ (DfsTree.node y outs).edges, a ≠ edgeItem s.g e') → Items.IsParent s.items a c →
      c < s.items.size ∧ (∀ x ∈ (DfsTree.node y outs).verts, c ≠ vertItem x) ∧
      ∀ e' ∈ (DfsTree.node y outs).edges, c ≠ edgeItem s.g e' := by
    intro a c _ h
    refine ⟨hC.ch_lt a c h, fun x hx hc => ?_, fun e' he' hc => ?_⟩
    · exact (hC.v_fresh _ (List.mem_cons_self ..) e cls _ rfl x hx).1 a (by rw [← hc]; exact h)
    · exact (hC.q_fresh _ (List.mem_cons_self ..) e' (Or.inr he')).1 a (by rw [← hc]; exact h)
  have hbelowE : ∀ a, (a < s.items.size ∧ (∀ x ∈ (DfsTree.node y outs).verts, a ≠ vertItem x) ∧
      ∀ e' ∈ (DfsTree.node y outs).edges, a ≠ edgeItem s.g e') → ∀ i,
      (Items.Below sE.items a i ↔ Items.Below s.items a i) :=
    fun a ha i => below_kept _ hOldcl (fun j hj => (H.items_kept j hj.1 hj.2.1 hj.2.2).1) ha i
  have hvlt : vertItem v < s.items.size := by show 1 + v < _; omega
  have hqlt : edgeItem s.g e < s.items.size := by show 1 + s.g.nv + e < _; omega
  have hvq : vertItem v ≠ edgeItem s.g e := vertItem_ne_edgeItem' hv
  have hOld_v : vertItem v < s.items.size ∧
      (∀ x ∈ (DfsTree.node y outs).verts, vertItem v ≠ vertItem x) ∧
      ∀ e' ∈ (DfsTree.node y outs).edges, vertItem v ≠ edgeItem s.g e' :=
    ⟨hvlt, fun x hx h => H.v_nc (by rw [vertItem_inj' h]; exact hx),
      fun _ _ => vertItem_ne_edgeItem' hv⟩
  have hOld_q : edgeItem s.g e < s.items.size ∧
      (∀ x ∈ (DfsTree.node y outs).verts, edgeItem s.g e ≠ vertItem x) ∧
      ∀ e' ∈ (DfsTree.node y outs).edges, edgeItem s.g e ≠ edgeItem s.g e' :=
    ⟨hqlt, fun x hx h => vertItem_ne_edgeItem' (H.verts_lt x hx) h.symm,
      fun e' he' h => H.e_ne (by rw [edgeItem_inj h]; exact he')⟩
  have hOld_span : ∀ t ∈ s.tstack, ∀ a ∈ t.spans.1 ++ t.spans.2, a < s.items.size ∧
      (∀ x ∈ (DfsTree.node y outs).verts, a ≠ vertItem x) ∧
      ∀ e' ∈ (DfsTree.node y outs).edges, a ≠ edgeItem s.g e' := fun t ht a ha =>
    ⟨hC.span_lt t ht a ha,
     fun x hx h => (hC.v_fresh _ (List.mem_cons_self ..) e cls _ rfl x hx).2.2.1 t ht
       (by rw [← h]; exact ha),
     fun e' he' h => (hC.q_fresh _ (List.mem_cons_self ..) e' (Or.inr he')).2.2 t ht
       (by rw [← h]; exact ha)⟩
  have hchq : Items.ch sE.items (edgeItem s.g e) = [] :=
    (H.items_kept _ hqlt hOld_q.2.1 hOld_q.2.2).1.trans hqf.2.1
  have hsle := H.size_le
  have hyq : vertItem y ≠ edgeItem s.g e := vertItem_ne_edgeItem' hy
  have hyv : vertItem y ≠ vertItem v := fun h =>
    H.v_nc (by rw [← vertItem_inj' h]; exact List.mem_cons_self ..)
  have hP : ∀ p c, Items.IsParent s'.items p c ↔
      Items.IsParent sE.items p c ∨ (p = vertItem v ∧ c = edgeItem s.g e) ∨
        (p = edgeItem s.g e ∧ (c ∈ t₁.spans.1 ∨ c = vertItem y)) := by
    intro p c
    simp only [Items.IsParent, hch]
    by_cases h1 : p = vertItem v
    · subst h1; simp [hvq]
    · by_cases h2 : p = edgeItem s.g e
      · subst h2; simp [h1, hchq]
      · simp [h1, h2]
  have hmono : ∀ {a i}, Items.Below sE.items a i → Items.Below s'.items a i := by
    intro a i h
    induction h with
    | refl => exact .refl
    | tail _ hs ih => exact ih.tail ((hP _ _).2 (Or.inl hs))
  have hbelow : ∀ a i, ¬ Items.Below sE.items a (vertItem v) →
      ¬ Items.Below sE.items a (edgeItem s.g e) →
      (Items.Below s'.items a i ↔ Items.Below sE.items a i) := by
    intro a i hav haq
    refine ⟨fun h => ?_, hmono⟩
    induction h with
    | refl => exact .refl
    | tail _ hji ih =>
      rcases (hP _ _).1 hji with h | ⟨rfl, rfl⟩ | ⟨rfl, -⟩
      · exact ih.tail h
      · exact absurd ih hav
      · exact absurd ih haq
  have hbelow_v : ∀ i, Items.Below s'.items (vertItem v) i ↔
      Items.Below sE.items (vertItem v) i ∨ i = edgeItem s.g e ∨
        (∃ a ∈ t₁.spans.1, Items.Below sE.items a i) ∨ Items.Below sE.items (vertItem y) i := by
    intro i
    constructor
    · intro h
      induction h with
      | refl => exact .inl .refl
      | tail _ hji ih =>
        rcases (hP _ _).1 hji with h | ⟨rfl, rfl⟩ | ⟨rfl, ha | rfl⟩
        · rcases ih with ih | rfl | ⟨a, ha, hb⟩ | ih
          · exact .inl (ih.tail h)
          · exact absurd h (by rw [Items.IsParent, hchq]; simp)
          · exact .inr (.inr (.inl ⟨a, ha, hb.tail h⟩))
          · exact .inr (.inr (.inr (ih.tail h)))
        · exact .inr (.inl rfl)
        · exact .inr (.inr (.inl ⟨_, ha, .refl⟩))
        · exact .inr (.inr (.inr .refl))
    · rintro (h | rfl | ⟨a, ha, h⟩ | h)
      · exact hmono h
      · exact .single ((hP _ _).2 (.inr (.inl ⟨rfl, rfl⟩)))
      · exact ((Relation.ReflTransGen.single ((hP _ _).2 (.inr (.inl ⟨rfl, rfl⟩)))).tail
          ((hP _ _).2 (.inr (.inr ⟨rfl, .inl ha⟩)))).trans (hmono h)
      · exact ((Relation.ReflTransGen.single ((hP _ _).2 (.inr (.inl ⟨rfl, rfl⟩)))).tail
          ((hP _ _).2 (.inr (.inr ⟨rfl, .inr rfl⟩)))).trans (hmono h)
  -- the child's edges: below `t₁`'s items or below `V y`
  have hEBy : ∀ e', e' < s.g.ne →
      (((∃ a ∈ t₁.spans.1, Items.Below sE.items a (edgeItem s.g e')) ∨
        Items.EdgeBelow s.g sE.items (vertItem y) e') ↔ e' ∈ (DfsTree.node y outs).edges) := by
    intro e' he'
    have hve := hC'.vert_edges e' (by rw [hgE]; exact he')
    rw [hgE] at hve
    rw [hedgesC]
    constructor
    · rintro (⟨a, ha, hb⟩ | h)
      · obtain ⟨o', ho', hs⟩ := hCT'.edges t₁ hmem1 e' (by rw [hgE]; exact he')
          (by rw [hgE]; exact ⟨a, List.mem_append_left _ ha, hb⟩)
        exact ⟨o'.1, by rw [← H.hdone']; exact List.mem_map_of_mem ho', hs⟩
      · obtain ⟨o', ho', -, hs⟩ := hve.1 h
        exact ⟨o'.1, by rw [← H.hdone']; exact List.mem_map_of_mem ho', hs⟩
    · rintro ⟨o', ho', hs⟩
      rw [← H.hdone'] at ho'
      obtain ⟨p, hp, rfl⟩ := List.mem_map.1 ho'
      by_cases hlt : p.1.cls.lowval (d + 1) < d + 1
      · obtain ⟨t, ht, hte⟩ := hCT'.cover p hp hlt e' (by rw [hgE]; exact he') hs
        rw [hgE] at hte
        rw [htop'eq] at ht
        rcases List.mem_cons.1 ht with rfl | ht
        · obtain ⟨a, ha, hb⟩ := hte
          rw [ht₁, List.append_nil] at ha
          exact .inl ⟨a, ha, hb⟩
        · rw [List.mem_singleton] at ht
          subst ht
          obtain ⟨a, ha, hb⟩ := hte
          simp at ha
          subst ha
          exact .inr hb
      · exact .inr (hve.2 ⟨p, hp, Nat.le_of_not_lt hlt, hs⟩)
  have hEB_v : ∀ e', e' < s.g.ne → (Items.EdgeBelow s.g s'.items (vertItem v) e' ↔
      Items.EdgeBelow s.g s.items (vertItem v) e' ∨ subEdges (.tree e cls (.node y outs)) e') := by
    intro e' he'
    rw [Items.EdgeBelow, hbelow_v, hbelowE _ hOld_v]
    constructor
    · rintro (h | h | h | h)
      · exact .inl h
      · exact .inr (Or.inl (edgeItem_inj h))
      · exact .inr (Or.inr ((hEBy e' he').1 (.inl h)))
      · exact .inr (Or.inr ((hEBy e' he').1 (.inr h)))
    · rintro (h | h)
      · exact .inl h
      rcases (show e' = e ∨ e' ∈ (DfsTree.node y outs).edges from h) with rfl | h
      · exact .inr (.inl rfl)
      · rcases (hEBy e' he').2 h with h | h
        · exact .inr (.inr (.inl h))
        · exact .inr (.inr (.inr h))
  have hspanB : ∀ t ∈ s.tstack, ∀ i ∈ t.spans.1 ++ t.spans.2, ∀ j,
      (Items.Below s'.items i j ↔ Items.Below s.items i j) := by
    intro t ht i hi j
    rw [← hbelowE _ (hOld_span t ht i hi)]
    apply hbelow
    · intro h
      have := Items.Below.eq_of_no_parent H.v_root h
      subst this
      exact Bool.false_ne_true (hC.vert_free t ht hi)
    · intro h
      have := Items.Below.eq_of_no_parent H.q_root h
      subst this
      exact hqf.2.2 t ht hi
  have hedges : ∀ t ∈ s.tstack, ∀ e',
      (TEntry.edges s.g s'.items t e' ↔ TEntry.edges s.g s.items t e') := fun t ht e' =>
    TEntry.edges_congr (fun i hi e' => hspanB t ht i hi _) e'
  have hAV : afterVert (done ++ [(.tree e cls (.node y outs), false)]) = afterVert done := by
    simp [afterVert]
  have hVL : ∀ x, x ∈ DfsOut.vertsList (done.map (·.1)) →
      x ∈ DfsOut.vertsList ((done ++ [((DfsOut.tree e cls (.node y outs), false) :
        DfsOut × Bool)]).map (·.1)) := by
    intro x hx
    rw [DfsOut.vertsList_eq] at hx ⊢
    rw [List.map_append, List.flatMap_append]
    exact List.mem_append_left _ hx
  obtain ⟨top, htop, hCT⟩ := hC.top
  have hmemT : ∀ t ∈ top, t ∈ s.tstack := fun t ht => by
    rw [htop]; exact List.mem_append_left _ ht
  have hbase_mem : ∀ k, k < base.length → base[k]! ∈ s.tstack := fun k hk => by
    rw [htop, getElem!_pos base k hk]; exact List.mem_append_right _ (List.getElem_mem hk)
  obtain ⟨above, below, hsp, -, -, -, -, -, hB⟩ := hCT.split
  simp only [Bool.false_eq_true, ↓reduceIte] at hsp
  obtain ⟨rfl, rfl⟩ := hsp
  have hbook := hC.vert_book rfl
  exact
    { top := ⟨top, by rw [hts, htop],
        { split := ⟨[], top, by simp, fun t ht => (nomatch ht), fun t ht => (nomatch ht),
            List.Pairwise.nil, fun t ht => (nomatch ht), fun t ht => (nomatch ht), hB⟩
          bot := fun t ht k hk => by rw [hsv k (by omega)]; exact hCT.bot t ht k hk
          vitems := fun t ht x hx hm =>
            (hCT.vitems t ht x (by rwa [hg] at hx) hm).imp id (hVL x)
          qitems := fun t ht e' he' hm =>
            let ⟨o, ho, hs⟩ := hCT.qitems t ht e' (by rwa [hg] at he') (by rwa [hg] at hm)
            ⟨o, List.mem_append_left _ ho, hs⟩
          edges := fun t ht e' he' hte => by
            rw [hg] at he' hte
            obtain ⟨o, ho, hs⟩ := hCT.edges t ht e' he' ((hedges t (hmemT t ht) e').1 hte)
            exact ⟨o, List.mem_append_left _ ho, hs⟩
          cover := fun o ho hlt e' he' hs => by
            rw [hg] at he' ⊢
            rcases List.mem_append.1 ho with h | h
            · obtain ⟨t, ht, hte⟩ := hCT.cover o h hlt e' he' hs
              exact ⟨t, ht, (hedges t (hmemT t ht) e').2 hte⟩
            · rw [List.mem_singleton] at h
              subst h
              exact absurd hlt (Nat.not_lt.2 hge)
          ret := fun h => (hCT.ret h).imp fun o ⟨ho, hl⟩ => ⟨List.mem_append_left _ ho, hl⟩ }⟩
      base_edges := ⟨hC.base_edges.1, fun k hk e' he' => by
        rw [hg] at he' ⊢
        exact (hedges _ (hbase_mem k hk) e').trans (hC.base_edges.2 k hk e' he')⟩
      base_bot := hC.base_bot
      noVert_after := fun _ => by
        have h := hC.noVert_after rfl
        simp only [afterVert, List.map_eq_nil_iff] at h
        simp [afterVert, List.filter_append, h]
      vfirst := fun t ht hvt _ => absurd hvt (by
        rw [hts, htop] at ht
        rcases List.mem_append.1 ht with h | h
        · exact hB t h
        · exact hC.base_bot t h)
      sv := fun k hk => by rw [hsv k hk]; exact hC.sv k hk
      sd := fun k hk => by rw [hsd k hk]; exact hC.sd k hk
      sv_d := by rw [hsv d (Nat.le_refl _)]; exact hC.sv_d
      path := fun k k' h h' => by rw [hsv k (by omega), hsv k' h']; exact hC.path k k' h h'
      v_root := fun p h => by
        rcases (hP p _).1 h with h | ⟨_, h⟩ | ⟨_, h | h⟩
        · exact H.v_root p h
        · exact hvq h
        · rcases hv1 v hv h with h | h
          · exact H.v_nc (by rw [h]; exact List.mem_cons_self ..)
          · exact H.v_nc (List.mem_cons_of_mem _ h)
        · exact hyv h.symm
      vert_free := fun t ht hm => by rw [hts] at ht; exact hC.vert_free t ht hm
      afterVert_ret := fun o ho => by rw [hAV] at ho; exact hC.afterVert_ret o ho
      hv_ret := fun h => (nomatch h)
      vert_book := fun _ => by
        rw [hg]
        have hadjO : ∀ u w, s.g.AdjIn (Items.EdgeBelow s.g s.items (vertItem v)) u w →
            s.g.AdjIn (Items.EdgeBelow s.g s'.items (vertItem v)) u w :=
          fun u w ⟨e₁, h1, h2, h3⟩ => ⟨e₁, h1, (hEB_v e₁ h1).2 (.inl h2), h3⟩
        have hsymm : ∀ a b, Relation.ReflTransGen
            (s.g.AdjIn (Items.EdgeBelow s.g s'.items (vertItem v))) a b →
            Relation.ReflTransGen (s.g.AdjIn (Items.EdgeBelow s.g s'.items (vertItem v))) b a := by
          intro a b h
          induction h with
          | refl => exact .refl
          | tail _ hs ih => exact (Relation.ReflTransGen.single (Graph.AdjIn.symm hs)).trans ih
        have hEe : Items.EdgeBelow s.g s'.items (vertItem v) e := (hEB_v e he).2 (.inr (Or.inl rfl))
        have hyv_step : Relation.ReflTransGen
            (s.g.AdjIn (Items.EdgeBelow s.g s'.items (vertItem v))) v y :=
          .single ⟨e, he, hEe, by
            rcases H.e_ends with h | h
            · exact Or.inr (by rw [← h])
            · exact Or.inl (Prod.ext (Prod.mk.inj h).2 (Prod.mk.inj h).1)⟩
        have hreachO : ∀ e₁, e₁ < s.g.ne → Items.EdgeBelow s.g s.items (vertItem v) e₁ →
            Relation.ReflTransGen (s.g.AdjIn (Items.EdgeBelow s.g s'.items (vertItem v))) v
              (s.g.edges[e₁]!).1 := by
          intro e₁ h1 hE1
          obtain ⟨o, ho, hlo, -⟩ := (hC.vert_edges e₁ h1).1 hE1
          have he0 : Items.EdgeBelow s.g s.items (vertItem v) o.1.e :=
            (hC.vert_edges _ (H.inc o ho).1).2 ⟨o, ho, hlo, subEdges_e _⟩
          have step : Relation.ReflTransGen (s.g.AdjIn (Items.EdgeBelow s.g s'.items (vertItem v)))
              v (s.g.edges[o.1.e]!).1 := by
            rcases (H.inc o ho).2 with h | h
            · rw [h]
            · exact .single ⟨o.1.e, (H.inc o ho).1, (hEB_v _ (H.inc o ho).1).2 (.inl he0),
                .inr (by rw [h])⟩
          exact step.trans (Relation.ReflTransGen.mono hadjO _ _
            (hbook.1 _ _ (H.inc o ho).1 h1 he0 hE1))
        have hreachC : ∀ e₁, e₁ < s.g.ne → e₁ ∈ (DfsTree.node y outs).edges →
            Relation.ReflTransGen (s.g.AdjIn (Items.EdgeBelow s.g s'.items (vertItem v))) v
              (s.g.edges[e₁]!).1 := fun e₁ _ hm =>
          hyv_step.trans (H.c_reach _ (fun e₂ h2 =>
            (hEB_v e₂ (H.edges_lt e₂ h2)).2 (.inr (.inr h2))) e₁ hm)
        have hreach : ∀ e₁, e₁ < s.g.ne → Items.EdgeBelow s.g s'.items (vertItem v) e₁ →
            Relation.ReflTransGen (s.g.AdjIn (Items.EdgeBelow s.g s'.items (vertItem v))) v
              (s.g.edges[e₁]!).1 := by
          intro e₁ h1 hE1
          rcases (hEB_v e₁ h1).1 hE1 with h | h
          · exact hreachO e₁ h1 h
          rcases (show e₁ = e ∨ e₁ ∈ (DfsTree.node y outs).edges from h) with rfl | h
          · rcases H.e_ends with h | h
            · rw [← h]; exact hyv_step
            · rw [(Prod.mk.inj h).2.symm]
          · exact hreachC e₁ h1 h
        refine ⟨fun e₁ e₂ h1 h2 hE1 hE2 => (hsymm _ _ (hreach e₁ h1 hE1)).trans (hreach e₂ h2 hE2),
          fun x e₁ e₂ h1 h2 hE1 hnE2 hi1 hi2 => ?_⟩
        rcases (hEB_v e₁ h1).1 hE1 with h | h
        · exact hbook.2 x e₁ e₂ h1 h2 h (fun h' => hnE2 ((hEB_v e₂ h2).2 (.inl h'))) hi1 hi2
        · rcases H.ends hge e₁ h h1 x hi1 with hx | hx
          · exact .inl hx
          · exact absurd ((hEB_v e₂ h2).2 (.inr (H.comp e₂ h2 x hi2 hx))) hnE2
      vert_disj := fun _ t ht e' he' hte hE => by
        rw [hts] at ht
        rw [hg] at he' hte hE
        rcases (hEB_v e' he').1 hE with h | h
        · exact hC.vert_disj rfl t ht e' he' ((hedges t ht e').1 hte) h
        · obtain ⟨i, hi, hb'⟩ := (hedges t ht e').1 hte
          have hqf' := hC.q_fresh _ (List.mem_cons_self ..) e' h
          exact hqf'.2.2 t ht (Items.Below.eq_of_no_parent hqf'.1 hb' ▸ hi)
      vert_edges := fun e' he' => by
        rw [hg] at he' ⊢
        rw [hEB_v e' he', hC.vert_edges e' he']
        constructor
        · rintro (⟨o, ho, hlo, hs⟩ | h)
          · exact ⟨o, List.mem_append_left _ ho, hlo, hs⟩
          · exact ⟨(.tree e cls (.node y outs), false),
              List.mem_append_right _ (List.mem_singleton_self _), hge, h⟩
        · rintro ⟨o, ho, hlo, hs⟩
          rcases List.mem_append.1 ho with h | h
          · exact .inl ⟨o, h, hlo, hs⟩
          · rw [List.mem_singleton] at h
            subst h
            exact .inr hs
      touch_bot := fun t ht ⟨e', he', hte⟩ => by
        rw [hts] at ht
        rw [hg] at he' hte ⊢
        obtain ⟨e₀, h0, hE0, hi0⟩ := hC.touch_bot t ht ⟨e', he', (hedges t ht e').1 hte⟩
        exact ⟨e₀, h0, (hedges t ht e₀).2 hE0, hi0⟩
      span_root := fun t ht i hi p h => by
        rw [hts] at ht
        rcases (hP p i).1 h with h | ⟨_, rfl⟩ | ⟨_, h | rfl⟩
        · have hO := hOld_span t ht i hi
          exact hC.span_root t ht i hi p (((H.items_kept i hO.1 hO.2.1 hO.2.2).2 p).1 h)
        · exact hqf.2.2 t ht hi
        · exact (hsp1 i h).2.2 t ht hi
        · exact hyf.2.2.1 t ht hi
      span_lt := fun t ht i hi => by
        rw [hts] at ht
        rw [hsz']
        exact Nat.lt_of_lt_of_le (hC.span_lt t ht i hi) hsle
      ch_lt := fun p c h => by
        rw [hsz']
        rcases (hP p c).1 h with h | ⟨_, rfl⟩ | ⟨_, h | rfl⟩
        · exact hC'.ch_lt p c h
        · exact Nat.lt_of_lt_of_le hqlt hsle
        · exact (hsp1 c h).2.1
        · exact hvy_lt
      disj := by
        rw [hts]
        refine List.Pairwise.imp_of_mem (fun {t t'} ht ht' h e' he' hte hte' => ?_) hC.disj
        rw [hg] at he' hte hte'
        exact h e' he' ((hedges t ht e').1 hte) ((hedges t' ht' e').1 hte')
      span_disj := by rw [hts]; exact hC.span_disj
      q_fresh := fun o ho e' hs => by
        obtain ⟨hns, hlt⟩ := H.rest_nd o ho e' hs
        obtain ⟨hr, hc', hm⟩ := hC.q_fresh o (List.mem_cons_of_mem _ ho) e' hs
        rw [hg]
        have hOld : edgeItem s.g e' < s.items.size ∧
            (∀ x ∈ (DfsTree.node y outs).verts, edgeItem s.g e' ≠ vertItem x) ∧
            ∀ e'' ∈ (DfsTree.node y outs).edges, edgeItem s.g e' ≠ edgeItem s.g e'' :=
          ⟨by show 1 + s.g.nv + e' < _; omega,
           fun x hx h => vertItem_ne_edgeItem' (H.verts_lt x hx) h.symm,
           fun e'' he'' h => hns (Or.inr (by rw [edgeItem_inj h]; exact he''))⟩
        have hK := H.items_kept _ hOld.1 hOld.2.1 hOld.2.2
        have hq' : edgeItem s.g e' ≠ edgeItem s.g e := fun h => hns (Or.inl (edgeItem_inj h))
        have hv' : edgeItem s.g e' ≠ vertItem v := (vertItem_ne_edgeItem' hv).symm
        have hy' : edgeItem s.g e' ≠ vertItem y := (vertItem_ne_edgeItem' hy).symm
        refine ⟨fun p h => ?_, ?_, fun t ht => ?_⟩
        · rcases (hP p _).1 h with h | ⟨_, h⟩ | ⟨_, h | h⟩
          · exact hr p ((hK.2 p).1 h)
          · exact hq' h
          · exact hns (Or.inr (hq1 e' hlt h))
          · exact hy' h
        · rw [hch, if_neg hv', if_neg hq', hK.1]; exact hc'
        · rw [hts] at ht; exact hm t ht
      v_fresh := fun o ho e₁ cls₁ child ho₁ w hw => by
        obtain ⟨hwC, hwv⟩ := H.rest_nv o ho e₁ cls₁ child ho₁ w hw
        obtain ⟨hr, hc', hm, hsvk, hvs, htch⟩ :=
          hC.v_fresh o (List.mem_cons_of_mem _ ho) e₁ cls₁ child ho₁ w hw
        have hOld : vertItem w < s.items.size ∧
            (∀ x ∈ (DfsTree.node y outs).verts, vertItem w ≠ vertItem x) ∧
            ∀ e'' ∈ (DfsTree.node y outs).edges, vertItem w ≠ edgeItem s.g e'' :=
          ⟨by show 1 + w < _; omega, fun x hx h => hwC (by rw [vertItem_inj' h]; exact hx),
           fun _ _ => vertItem_ne_edgeItem' hwv⟩
        have hK := H.items_kept _ hOld.1 hOld.2.1 hOld.2.2
        have hwv' : vertItem w ≠ vertItem v := fun h =>
          hsvk d (Nat.le_refl _) (hC.sv_d.trans (vertItem_inj' h).symm)
        have hwq : vertItem w ≠ edgeItem s.g e := vertItem_ne_edgeItem' hwv
        have hwy : vertItem w ≠ vertItem y := fun h =>
          hwC (by rw [vertItem_inj' h]; exact List.mem_cons_self ..)
        refine ⟨fun p h => ?_, ?_, fun t ht => ?_, fun k hk => ?_, fun t ht => ?_, fun t ht => ?_⟩
        · rcases (hP p _).1 h with h | ⟨_, h⟩ | ⟨_, h | h⟩
          · exact hr p ((hK.2 p).1 h)
          · exact hwq h
          · rcases hv1 w hwv h with h | h
            · exact hwC (by rw [h]; exact List.mem_cons_self ..)
            · exact hwC (List.mem_cons_of_mem _ h)
          · exact hwy h
        · rw [hch, if_neg hwv', if_neg hwq, hK.1]; exact hc'
        · rw [hts] at ht; exact hm t ht
        · rw [hsv k hk]; exact hsvk k hk
        · rw [hts] at ht; exact hvs t ht
        · rw [hts] at ht
          rw [hg]
          rintro ⟨e'', he'', hte, hi⟩
          exact htch t ht ⟨e'', he'', (hedges t ht e'').1 hte, hi⟩ }

/-- The item frame of `finishBoundary` at a component edge: the child's frame, then `t₁`'s side-1
items and `V y` moved under `Q e`, `Q e` under `V v`. -/
theorem outFrame_comp {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s : WalkState} {e : Nat} {cls : OutClass} {y : Nat} {outs : List DfsOut} {L : List TEntry}
    {push : Bool} {done' : List (DfsOut × Bool)} {hv' : Bool} {bE' : List (Nat → Prop)}
    {sv' : List Nat} {sd' : List Bool} {sE : WalkState} {dir' : Bool} {L' : List TEntry}
    {push' : Bool} {D₃ : Array Bool} {s' : WalkState}
    (H : TreeSite v d done rest hasVert base bE sv sd s e cls y outs L push done' hv' bE' sv' sd'
      sE dir' L' push' D₃)
    {t₁ : TEntry} {f₂ : Nat}
    (hsh : sE.tstack = t₁ :: ⟨y, d + 1, f₂, ([], [vertItem y])⟩ :: s.tstack)
    (hg : s'.g = s.g) (hts : s'.tstack = s.tstack)
    (hsz' : s'.items.size = sE.items.size)
    (hch : ∀ j, Items.ch s'.items j =
      if j = vertItem v then Items.ch sE.items (vertItem v) ++ [edgeItem s.g e]
      else if j = edgeItem s.g e then t₁.spans.1 ++ [vertItem y]
      else Items.ch sE.items j) :
    OutFrame v (.tree e cls (.node y outs)) base s s' := by
  have hC := H.ctx
  have hsz := H.hsz
  have hsle := H.size_le
  have hv := H.hv
  have he := H.e_lt
  have hgE := H.gE
  have hqlt : edgeItem s.g e < s.items.size := by show 1 + s.g.nv + e < _; omega
  have hqE : Items.ch sE.items (edgeItem s.g e) = [] := by
    rw [(H.items_kept _ hqlt (fun x hx => (vertItem_ne_edgeItem' (H.verts_lt x hx)).symm)
      (fun e' he' h => H.e_ne ((edgeItem_inj h) ▸ he'))).1]
    exact (hC.q_fresh _ (List.mem_cons_self ..) e (subEdges_e _)).2.1
  have hvq : vertItem v ≠ edgeItem s.g e := vertItem_ne_edgeItem' hv
  have hjq' : ∀ j, (∀ e' ∈ (DfsOut.tree e cls (.node y outs)).edges,
      j ≠ edgeItem (pushEnd sE D₃ L').g e') → j ≠ edgeItem s.g e := fun j hje => by
    have := hje e (by simp [DfsOut.edges])
    rwa [show (pushEnd sE D₃ L').g = s.g from hgE] at this
  obtain ⟨top₀, htop₀, -⟩ := hC.top
  refine H.frame_child.trans (outFrame_of_mod
    (fun p => p = vertItem v ∨ p = edgeItem s.g e)
    (by simp only [pushEnd]; omega) (hg.trans hgE.symm) ?_ ?_ ?_ ?_)
  · intro p hp
    rw [hch p]
    simp only [not_or] at hp
    simp [pushEnd, hp.1, hp.2]
  · intro p j hp hj
    obtain ⟨hjlt, hjv, hxv, hje, hjt⟩ := hj
    have hjq : j ≠ edgeItem s.g e := hjq' j hje
    have hjy : j ≠ vertItem y := hxv y (by simp [DfsOut.verts, DfsTree.verts])
    have hjt1 : j ∉ t₁.spans.1 := fun h =>
      hjt (L' ++ t₁ :: ⟨y, d + 1, f₂, ([], [vertItem y])⟩ :: top₀) (by simp [pushEnd, hsh, htop₀])
        t₁ (List.mem_append_right _ (List.mem_cons_self ..)) (List.mem_append_left _ h)
    simp only [Items.IsParent, pushEnd]
    rw [hch p]
    rcases hp with rfl | rfl
    · simp [hjq]
    · simp [hvq.symm, hqE, hjt1, hjy]
  · intro j hj
    obtain ⟨hjlt, hjv, -, hje, -⟩ := hj
    simp only [not_or]
    exact ⟨hjv, hjq' j hje⟩
  · intro top' htop' t ht
    rw [hts] at htop'
    exact ⟨L' ++ t₁ :: ⟨y, d + 1, f₂, ([], [vertItem y])⟩ :: top', by simp [pushEnd, hsh, htop'],
      List.mem_append_right _ (List.mem_cons_of_mem _ (List.mem_cons_of_mem _ ht))⟩

/-- The tree-edge step for a component (`lowval = d`): `finishBoundary` pops the child's `(y, d)`
entry and its vertex entry (`tree_comp_shape`); `earCtx_comp` over the concrete post-state. -/
theorem ctx_step_tree_comp {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s : WalkState} {e : Nat} {cls : OutClass} {y : Nat} {outs : List DfsOut} {L : List TEntry}
    {push : Bool} {done' : List (DfsOut × Bool)} {hv' : Bool} {bE' : List (Nat → Prop)}
    {sv' : List Nat} {sd' : List Bool} {sE : WalkState} {dir' : Bool} {L' : List TEntry}
    {push' : Bool} {D₃ : Array Bool}
    (H : TreeSite v d done rest hasVert base bE sv sd s e cls y outs L push done' hv' bE' sv' sd'
      sE dir' L' push' D₃)
    (hc : cls = .component) :
    wp (finishEdge v d (.tree e cls (.node y outs)) (L ++ s.tstack).length (hasVert || push))
      (fun hv'' s' => EarCtx v d (done ++ [(.tree e cls (.node y outs), hasVert || push)]) rest hv''
        base bE sv sd s' ∧ OutFrame v (.tree e cls (.node y outs)) base s s' ∧
        (CtxShape v d done hasVert base s →
          CtxShape v d (done ++ [(.tree e cls (.node y outs), hasVert || push)]) hv'' base s'))
      (pushEnd sE D₃ L') := by
  have hC := H.ctx
  have hge : d ≤ cls.lowval d := by rw [hc]; exact Nat.le_refl d
  have hhv : hasVert = false := by
    cases hhv : hasVert
    · rfl
    · exact absurd (hC.hv_ret hhv) (no_ret_before_boundary hge H.cls_ret H.rank)
  have hpf : push = false := by
    cases hp : push
    · rfl
    · exact absurd (H.hpush.1 hp).2.1 (Nat.not_lt.2 hge)
  obtain ⟨hv't, t₁, f₂, hsh, -, -, ht₁⟩ := tree_comp_shape H hc
  have hp' : push' = false := by
    cases hp : push'
    · rfl
    · exact absurd (H.hpush'.1 hp) (by rw [hv't]; simp)
  have hL : L = [] := by rw [H.hL, hpf]; rfl
  have hL' : L' = [] := by rw [H.hL', hp']; rfl
  rw [hL, List.nil_append] at hsh
  have hsz := H.hsz
  have hsle := H.size_le
  have hv := H.hv
  have he := H.e_lt
  have hvlt : vertItem v < sE.items.size := by show 1 + v < _; omega
  have hqlt : edgeItem s.g e < sE.items.size := by show 1 + s.g.nv + e < _; omega
  have hvq : vertItem v ≠ edgeItem s.g e := vertItem_ne_edgeItem' H.hv
  have hgE := H.gE
  have hsdk : ∀ k, k < d → D₃[k]! = s.stackDir[k]! := fun k hk => by
    rw [H.sd₃ k (Nat.le_of_lt hk), H.sdlo k (Nat.le_of_lt hk)]
    exact getElem!_set!_ne' _ _ _ _ (by omega)
  subst hhv hpf
  simp only [Bool.or_self]
  have hge' : cls.lowval d ≥ d := hge
  rw [finishEdge_eq]
  simp only [finishEdge', wp_bind, wp_get, wp_stackDir, DfsOut.cls, hge', ↓reduceIte]
  unfold finishBoundary
  subst hc
  have hne : (d == d + 1) = false := by simp
  simp only [wp_bind, wp_modifyItem, wp_modify, wp_popTstack, wp_pure, DfsOut.cls, DfsOut.e,
    OutClass.isTree, OutClass.lowval, hne, Bool.false_eq_true, ↓reduceIte, pushEnd,
    hL', List.nil_append, hsh, List.head!_cons, List.tail_cons]
  refine (fun hch => ⟨earCtx_comp H rfl hsh ht₁ hgE (by simp) (fun k hk => H.svlo k hk) hsdk (by simp) hch,
    outFrame_comp H hsh hgE (by simp) (by simp) hch,
    fun hS => ctxShape_same hS hge' (by simp [hsh]) hsdk⟩) ?_
  intro j
  simp only [Items.ch, Array.getElem?_modify]
  by_cases h1 : j = vertItem v
  · subst h1
    simp [hgE, hvq.symm, Array.getElem?_eq_getElem hvlt]
  · by_cases h2 : j = edgeItem s.g e
    · subst h2
      simp [hgE, h1, hvq, Array.getElem?_eq_getElem hqlt]
    · simp [hgE, h1, h2, Ne.symm h1, Ne.symm h2]

/-- `finishEdge` at a returning tree edge, block by block: `finishRest` runs from `feS₃` (vertex
entry pushed) or from `feS₂`. -/
theorem wp_finishEdge_tree (curV d : Nat) (o : DfsOut) (orig : Nat) (hasVert : Bool)
    (s : WalkState) (Q : Bool → WalkState → Prop) (ht : o.cls.isTree = true)
    (hr : o.cls.lowval d < d) :
    wp (finishEdge curV d o orig hasVert) Q s ↔
      if hasVert then
        wp (finishRest curV d (o.cls.lowval d) o.cls.isType1 true (feB₃ curV d o orig s)) Q
          (feS₃ curV d o orig s)
      else wp (finishRest curV d (o.cls.lowval d) o.cls.isType1 false (feSingle d o s)) Q
        (feS₂ d o s) := by
  have hge : ¬ (o.cls.lowval d ≥ d) := Nat.not_le.2 hr
  rw [finishEdge_eq]
  cases hasVert <;>
    simp only [finishEdge', finishTree, wp_bind, wp_get, wp_stackDir, hge, ↓reduceIte, ht,
      closeVert_eq, Bool.false_eq_true, wp_makeVs, wp_modifyItem] <;>
    exact Iff.rfl

/-- The parent context after a returning tree out from `RetTop` and the final top `T`: `R` itself
(the vertex entry was already pushed), or `V v` pushed above `R` (`isSingle`), or merged into `R`'s
top entry. -/
theorem earCtx_ret {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s : WalkState} {e : Nat} {cls : OutClass} {y : Nat} {outs : List DfsOut} {L : List TEntry}
    {push : Bool} {done' : List (DfsOut × Bool)} {hv' : Bool} {bE' : List (Nat → Prop)}
    {sv' : List Nat} {sd' : List Bool} {sE : WalkState} {dir' : Bool} {L' : List TEntry}
    {push' : Bool} {D₃ : Array Bool}
    (H : TreeSite v d done rest hasVert base bE sv sd s e cls y outs L push done' hv' bE' sv' sd'
      sE dir' L' push' D₃)
    (hr : cls.lowval d < d) {hv₀ : Bool}
    (hb : (hasVert || push) = hv₀) {sX : WalkState} {R : List TEntry}
    (HR : RetTop v d s e cls y outs L hv₀ sX R) (T : List TEntry) {s' : WalkState}
    (hT : (hv₀ = true ∧ T = R) ∨ (hv₀ = false ∧ ∃ f dir,
      T = ⟨v, d, f, setSides dir [vertItem v] []⟩ :: R ∨
      ∃ c R', R = c :: R' ∧ T = TEntry.mergeInto ⟨v, d, f, setSides dir [vertItem v] []⟩ c :: R'))
    (hts : s'.tstack = T ++ (L ++ s.tstack)) (hg : s'.g = sX.g) (hit : s'.items = sX.items)
    (hsv : s'.stackVerts = sX.stackVerts) (hsd : ∀ k, k < d → s'.stackDir[k]! = sX.stackDir[k]!)
    (hnx : s'.nxtEdgeIdx = sX.nxtEdgeIdx) :
    EarCtx v d (done ++ [(.tree e cls (.node y outs), hv₀)]) rest true base bE sv sd s' := by
  have hC := H.ctx
  have hnc := H.done_rest
  have hv : v < s.g.nv := H.hv
  have he : e < s.g.ne := H.e_lt
  have hsz : 1 + s.g.nv + s.g.ne ≤ s.items.size := H.hsz
  set C : DfsTree := .node y outs with hCdef
  set o : DfsOut := .tree e cls C with ho
  obtain ⟨-, hRg, hRsv, hRsd, -, hRsz, hRnx, hRkept, hRch, hRlt, hRroot, hRnew, hRtb, hRE, hRtouch,
    hRdisj, hRsdisj, hRvs, hRvert, hRnoV⟩ := HR
  rw [← hit] at hRsz hRkept hRch hRlt hRroot hRtb hRE hRtouch hRdisj hRvert hRnoV
  rw [← hnx] at hRnx hRvert
  have hg' : s'.g = s.g := hg.trans hRg
  have hsvX : ∀ k, k ≤ d → s'.stackVerts[k]! = s.stackVerts[k]! := fun k hk => by
    rw [hsv]; exact hRsv k hk
  have hsdX : ∀ k, k < d → s'.stackDir[k]! = s.stackDir[k]! := fun k hk =>
    (hsd k hk).trans (hRsd k hk)
  have hvd : s.stackVerts[d]! = v := hC.sv_d
  have hvne : ∀ k, k < d → v ≠ s.stackVerts[k]! := fun k hk h =>
    hC.path k d hk (Nat.le_refl _) (by rw [hvd, h])
  have hoR : o ∈ o :: rest := List.mem_cons_self ..
  have hqf := fun e' (hs : subEdges o e') => hC.q_fresh o hoR e' hs
  have hvf := fun x (hx : x ∈ C.verts) => hC.v_fresh o hoR e cls C rfl x hx
  have hCv : ∀ x ∈ C.verts, x < s.g.nv := H.verts_lt
  have hoe : ∀ e', subEdges o e' → e' < s.g.ne := by
    rintro e' (h | h)
    · exact h ▸ he
    · exact H.edges_lt e' h
  have hsube : subEdges o e := Or.inl rfl
  -- L
  have hLprop : ∀ t ∈ L, t.vStart = v ∧ t.topDepth = d ∧ t.firstIdx = s.nxtEdgeIdx ∧
      t.spans.1 ++ t.spans.2 = [vertItem v] ∧ hasVert = false ∧ push = true := by
    intro t ht
    rw [H.hL] at ht
    split at ht
    · rw [List.mem_singleton] at ht; subst ht
      exact ⟨rfl, rfl, rfl, spans_setSides_single _ _, (H.hpush.1 ‹_›).1, ‹_›⟩
    · exact absurd ht List.not_mem_nil
  have hhv0 : hv₀ = false → hasVert = false := fun h => by
    subst h; exact (Bool.or_eq_false_iff.1 hb).1
  have hpush0 : hv₀ = false → push = false := fun h => by
    subst h; exact (Bool.or_eq_false_iff.1 hb).2
  have hL0 : hv₀ = false → L = [] := fun h => by rw [H.hL, hpush0 h]; rfl
  -- items untouched by the child
  let Old : ItemId → Prop := fun j => j < s.items.size ∧ (∀ x ∈ C.verts, j ≠ vertItem x) ∧
    (∀ e' ∈ C.edges, j ≠ edgeItem s.g e') ∧ j ≠ edgeItem s.g e
  have hOldcl : ∀ a c, Old a → Items.IsParent s.items a c → Old c := by
    intro a c _ hp
    refine ⟨hC.ch_lt a c hp, fun x hx h => (hvf x hx).1 a (h ▸ hp),
      fun e' he' h => (hqf e' (Or.inr he')).1 a (h ▸ hp), fun h => (hqf e hsube).1 a (h ▸ hp)⟩
  have hkept : ∀ j, Old j → Items.ch s'.items j = Items.ch s.items j :=
    fun j hj => (hRkept j hj.1 hj.2.1 hj.2.2.1 hj.2.2.2).1
  have hPO : ∀ j, Old j → ∀ p, Items.IsParent s'.items p j ↔ Items.IsParent s.items p j :=
    fun j hj => (hRkept j hj.1 hj.2.1 hj.2.2.1 hj.2.2.2).2
  have hbelowO : ∀ a, Old a → ∀ i, Items.Below s'.items a i ↔ Items.Below s.items a i :=
    fun a ha i => below_kept Old hOldcl hkept ha i
  have hOld_v : Old (vertItem v) :=
    ⟨by show 1 + v < _; omega, fun x hx h => H.v_nc ((vertItem_inj' h) ▸ hx),
      fun _ _ => vertItem_ne_edgeItem' hv, vertItem_ne_edgeItem' hv⟩
  have hOld_span : ∀ t ∈ s.tstack, ∀ i ∈ t.spans.1 ++ t.spans.2, Old i := by
    intro t ht i hi
    refine ⟨hC.span_lt t ht i hi, fun x hx h => (hvf x hx).2.2.1 t ht (h ▸ hi),
      fun e' he' h => (hqf e' (Or.inr he')).2.2 t ht (h ▸ hi), fun h => (hqf e hsube).2.2 t ht (h ▸ hi)⟩
  have hedges₀ : ∀ t ∈ s.tstack, TEntry.edges s.g s'.items t = TEntry.edges s.g s.items t := by
    intro t ht
    funext e'
    exact propext (TEntry.edges_congr (fun i hi e'' => hbelowO i (hOld_span t ht i hi) _) e')
  have hedges : ∀ t ∈ s.tstack, TEntry.edges s'.g s'.items t = TEntry.edges s.g s.items t := by
    intro t ht; rw [hg']; exact hedges₀ t ht
  have hEBv : ∀ e', Items.EdgeBelow s.g s'.items (vertItem v) e' ↔
      Items.EdgeBelow s.g s.items (vertItem v) e' := fun e' => hbelowO _ hOld_v _
  have hVE : ∀ e', e' < s.g.ne → (Items.EdgeBelow s.g s'.items (vertItem v) e' ↔
      ∃ o' ∈ done, d ≤ o'.1.cls.lowval d ∧ subEdges o'.1 e') := fun e' he' => by
    rw [hEBv]; exact hC.vert_edges e' he'
  have hUe : ∀ t : TEntry, t.spans.1 ++ t.spans.2 = [vertItem v] → ∀ e',
      TEntry.edges s.g s'.items t e' ↔ Items.EdgeBelow s.g s'.items (vertItem v) e' := by
    intro t hts' e'
    unfold TEntry.edges
    rw [hts']
    constructor
    · rintro ⟨i, hi, h⟩; rw [List.mem_singleton] at hi; subst hi; exact h
    · intro h; exact ⟨_, List.mem_singleton_self _, h⟩
  have hLe : ∀ t ∈ L, ∀ e', TEntry.edges s.g s'.items t e' ↔
      Items.EdgeBelow s.g s'.items (vertItem v) e' := fun t ht => hUe t (hLprop t ht).2.2.2.1
  have hTVv : ∀ e', e' < s.g.ne → Items.EdgeBelow s.g s'.items (vertItem v) e' →
      s.g.Touches (Items.EdgeBelow s.g s'.items (vertItem v)) v := by
    intro e' he' h
    obtain ⟨o', ho', hlo, -⟩ := (hVE e' he').1 h
    exact ⟨o'.1.e, (H.inc o' ho').1, (hVE _ (H.inc o' ho').1).2 ⟨o', ho', hlo, subEdges_e _⟩,
      (H.inc o' ho').2⟩
  have hTV : ∀ x, s.g.Touches (Items.EdgeBelow s.g s'.items (vertItem v)) x →
      ∃ o' ∈ done, ∃ e', e' < s.g.ne ∧ subEdges o'.1 e' ∧ s.g.Inc e' x := by
    rintro x ⟨e', he', hb', hx⟩
    obtain ⟨o', ho', -, hs⟩ := (hVE e' he').1 hb'
    exact ⟨o', ho', e', he', hs, hx⟩
  have hRspan_cls : ∀ c ∈ R, ∀ i ∈ c.spans.1 ++ c.spans.2, i < s.items.size →
      (∃ x ∈ C.verts, i = vertItem x) ∨ (∃ e' ∈ C.edges, i = edgeItem s.g e') ∨ i = edgeItem s.g e :=
    fun c hc i hi => (hRnew c hc i hi).2
  have hRnoVv : ∀ c ∈ R, vertItem v ∉ c.spans.1 ++ c.spans.2 := by
    intro c hc hi
    rcases hRspan_cls c hc _ hi hOld_v.1 with ⟨x, hx, h⟩ | ⟨e', -, h⟩ | h
    · exact H.v_nc ((vertItem_inj' h) ▸ hx)
    · exact vertItem_ne_edgeItem' hv h
    · exact vertItem_ne_edgeItem' hv h
  have hRVdisj : ∀ e', e' < s.g.ne → (∃ c ∈ R, TEntry.edges s.g s'.items c e') →
      ¬ Items.EdgeBelow s.g s'.items (vertItem v) e' := by
    intro e' he' hc hb'
    obtain ⟨o', ho', -, hs⟩ := (hVE e' he').1 hb'
    exact H.nd o' ho' e' hs ((hRE e' he').1 hc)
  have hRolddisj : ∀ e', e' < s.g.ne → (∃ c ∈ R, TEntry.edges s.g s'.items c e') →
      ∀ t' ∈ s.tstack, ¬ TEntry.edges s.g s.items t' e' := by
    intro e' he' hc t' ht' ⟨i, hi, hb'⟩
    have hs := (hRE e' he').1 hc
    have := Items.Below.eq_of_no_parent (hqf e' hs).1 hb'
    subst this
    exact (hqf e' hs).2.2 t' ht' hi
  have hoD : (o, hv₀) ∈ done ++ [(o, hv₀)] := List.mem_append_right _ (List.mem_singleton_self _)
  have hAVm : ∀ o', o' ∈ afterVert (done ++ [(o, hv₀)]) ↔
      o' ∈ afterVert done ∨ (hv₀ = true ∧ o' = o) := by
    intro o'
    cases hv₀ <;> simp [afterVert, List.filter_append]
  -- T
  have hus : ∀ f dir, (⟨v, d, f, setSides dir [vertItem v] []⟩ : TEntry).spans.1 ++
      (⟨v, d, f, setSides dir [vertItem v] []⟩ : TEntry).spans.2 = [vertItem v] :=
    fun f dir => spans_setSides_single dir (vertItem v)
  have hTmem : ∀ t ∈ T, t ∈ R ∨
      (hv₀ = false ∧ t.vStart = v ∧ t.topDepth = d ∧ t.spans.1 ++ t.spans.2 = [vertItem v]) ∨
      (hv₀ = false ∧ ∃ c ∈ R, TEntry.edges s.g s'.items c e ∧ t.vStart = c.vStart ∧
        t.topDepth ≤ d ∧
        (∀ i, i ∈ t.spans.1 ++ t.spans.2 ↔ i = vertItem v ∨ i ∈ c.spans.1 ++ c.spans.2) ∧
        (∀ e', TEntry.edges s.g s'.items t e' ↔
          Items.EdgeBelow s.g s'.items (vertItem v) e' ∨ TEntry.edges s.g s'.items c e')) := by
    intro t ht
    rcases hT with ⟨-, rfl⟩ | ⟨h0, f, dir, rfl | ⟨c, R', hR, rfl⟩⟩
    · exact .inl ht
    · rcases List.mem_cons.1 ht with rfl | ht
      · exact .inr (.inl ⟨h0, rfl, rfl, hus f dir⟩)
      · exact .inl ht
    · obtain ⟨⟨c₀, R₀, hR₀, hce⟩, -⟩ := hRnoV h0
      rw [hR] at hR₀
      obtain ⟨rfl, rfl⟩ := List.cons.inj hR₀
      rcases List.mem_cons.1 ht with rfl | ht
      · refine .inr (.inr ⟨h0, c, by rw [hR]; exact List.mem_cons_self .., hce, rfl,
          Nat.min_le_right _ _, fun i => ?_, fun e' => ?_⟩)
        · rw [mem_mergeInto_spans, hus f dir, List.mem_singleton]
        · rw [TEntry.edges_mergeInto, hUe _ (hus f dir)]
      · exact .inl (by rw [hR]; exact List.mem_cons_of_mem _ ht)
  have hTcov : ∀ c ∈ R, ∃ t ∈ T, ∀ e', TEntry.edges s.g s'.items c e' → TEntry.edges s.g s'.items t e' := by
    intro c hc
    rcases hT with ⟨-, rfl⟩ | ⟨-, f, dir, rfl | ⟨c₀, R', hR, rfl⟩⟩
    · exact ⟨c, hc, fun _ h => h⟩
    · exact ⟨c, List.mem_cons_of_mem _ hc, fun _ h => h⟩
    · rw [hR] at hc
      rcases List.mem_cons.1 hc with rfl | hc
      · exact ⟨_, List.mem_cons_self .., fun e' h => (TEntry.edges_mergeInto _ _ _).2 (.inr h)⟩
      · exact ⟨c, List.mem_cons_of_mem _ hc, fun _ h => h⟩
  have hTvs : ∀ t ∈ T, t.vStart = v ∨ t.vStart ∈ C.verts := by
    intro t ht
    rcases hTmem t ht with h | ⟨-, h, -⟩ | ⟨-, c, hc, -, h, -⟩
    · exact hRvs t h
    · exact .inl h
    · rw [h]; exact hRvs c hc
  have hTspan : ∀ t ∈ T, ∀ i ∈ t.spans.1 ++ t.spans.2, (∃ c ∈ R, i ∈ c.spans.1 ++ c.spans.2) ∨ i = vertItem v := by
    intro t ht i hi
    rcases hTmem t ht with h | ⟨-, -, -, h⟩ | ⟨-, c, hc, -, -, -, h, -⟩
    · exact .inl ⟨t, h, hi⟩
    · rw [h, List.mem_singleton] at hi; exact .inr hi
    · rcases (h i).1 hi with h' | h'
      · exact .inr h'
      · exact .inl ⟨c, hc, h'⟩
  have hTedge : ∀ t ∈ T, ∀ e', TEntry.edges s.g s'.items t e' →
      (∃ c ∈ R, TEntry.edges s.g s'.items c e') ∨ Items.EdgeBelow s.g s'.items (vertItem v) e' := by
    intro t ht e' h
    rcases hTmem t ht with h₁ | ⟨-, -, -, h₁⟩ | ⟨-, c, hc, -, -, -, -, h₁⟩
    · exact .inl ⟨t, h₁, h⟩
    · exact .inr ((hUe t h₁ e').1 h)
    · rcases (h₁ e').1 h with h' | h'
      · exact .inr h'
      · exact .inl ⟨c, hc, h'⟩
  have hTpwE : T.Pairwise fun t t' => ∀ e', e' < s.g.ne → TEntry.edges s.g s'.items t e' →
      ¬ TEntry.edges s.g s'.items t' e' := by
    rcases hT with ⟨-, rfl⟩ | ⟨h0, f, dir, rfl | ⟨c, R', hR, rfl⟩⟩
    · exact hRdisj
    · refine List.pairwise_cons.2 ⟨fun t' ht' e' he' h1 h2 => ?_, hRdisj⟩
      exact hRVdisj e' he' ⟨t', ht', h2⟩ ((hUe _ (hus f dir) e').1 h1)
    · subst hR
      have hcR' := List.pairwise_cons.1 hRdisj
      refine List.pairwise_cons.2 ⟨fun t' ht' e' he' h1 h2 => ?_, hcR'.2⟩
      rcases (TEntry.edges_mergeInto _ _ _).1 h1 with h1 | h1
      · exact hRVdisj e' he' ⟨t', List.mem_cons_of_mem _ ht', h2⟩ ((hUe _ (hus f dir) e').1 h1)
      · exact hcR'.1 t' ht' e' he' h1 h2
  have hTpwS : T.Pairwise fun t t' => ∀ i ∈ t.spans.1 ++ t.spans.2, i ∉ t'.spans.1 ++ t'.spans.2 := by
    rcases hT with ⟨-, rfl⟩ | ⟨h0, f, dir, rfl | ⟨c, R', hR, rfl⟩⟩
    · exact hRsdisj
    · refine List.pairwise_cons.2 ⟨fun t' ht' i hi hi' => ?_, hRsdisj⟩
      rw [hus f dir, List.mem_singleton] at hi; subst hi
      exact hRnoVv t' ht' hi'
    · subst hR
      have hcR' := List.pairwise_cons.1 hRsdisj
      refine List.pairwise_cons.2 ⟨fun t' ht' i hi hi' => ?_, hcR'.2⟩
      rcases (mem_mergeInto_spans _ _ _).1 hi with h1 | h1
      · rw [hus f dir, List.mem_singleton] at h1; subst h1
        exact hRnoVv t' (List.mem_cons_of_mem _ ht') hi'
      · exact hcR'.1 t' ht' i h1 hi'
  -- the old top
  obtain ⟨top, htop, hCT⟩ := hC.top
  have hmemT : ∀ t ∈ top, t ∈ s.tstack := fun t ht => by
    rw [htop]; exact List.mem_append_left _ ht
  have hmemS : ∀ t ∈ s'.tstack, t ∈ T ∨ t ∈ L ∨ t ∈ s.tstack := by
    intro t ht
    rw [hts] at ht
    rcases List.mem_append.1 ht with h | h
    · exact .inl h
    · exact .inr (List.mem_append.1 h)
  have hmemTop : ∀ t ∈ T ++ L ++ top, t ∈ T ∨ t ∈ L ∨ t ∈ top := by
    intro t ht
    rcases List.mem_append.1 ht with h | h
    · exact (List.mem_append.1 h).imp id Or.inl
    · exact .inr (.inr h)
  have hmemTop' : ∀ t ∈ T ++ L ++ top, t ∈ s'.tstack := by
    intro t ht
    rw [hts, htop]
    rcases hmemTop t ht with h | h | h
    · exact List.mem_append_left _ h
    · exact List.mem_append_right _ (List.mem_append_left _ h)
    · exact List.mem_append_right _ (List.mem_append_right _ (List.mem_append_left _ h))
  have hVL : ∀ x, x ∈ DfsOut.vertsList (done.map (·.1)) →
      x ∈ DfsOut.vertsList ((done ++ [(o, hv₀)]).map (·.1)) := by
    intro x hx
    rw [DfsOut.vertsList_eq] at hx ⊢
    rw [List.map_append, List.flatMap_append]
    exact List.mem_append_left _ hx
  have hCVL : ∀ x ∈ C.verts, x ∈ DfsOut.vertsList ((done ++ [(o, hv₀)]).map (·.1)) := by
    intro x hx
    exact mem_vertsList.2 ⟨o, by rw [List.map_append]; exact List.mem_append_right _ (List.mem_singleton_self _),
      e, cls, C, rfl, hx⟩
  have hCE' : ∀ t ∈ s.tstack, CtxEntry v d s t → CtxEntry v d s' t := by
    intro t ht h
    obtain ⟨e', h1, h2⟩ := h.nonempty
    exact
      { vStart := h.vStart
        depth := h.depth
        nonempty := ⟨e', by rw [hg']; exact h1, by rw [hedges t ht]; exact h2⟩
        touch_bot := by rw [hedges t ht, hg']; exact h.touch_bot
        touch_top := by rw [hedges t ht, hg', hsvX _ h.depth.le]; exact h.touch_top
        side := by rw [hsdX _ h.depth]; exact h.side
        att := fun x hx hxv hni => by
          rw [hedges t ht, hg'] at hx hni
          obtain ⟨k, hk1, hk2, hk3⟩ := h.att x hx hxv hni
          exact ⟨k, hk1, hk2, by rw [hsvX k hk2]; exact hk3⟩ }
  have hCS' : ∀ t ∈ s.tstack, t.topDepth < d → CtxSingle v s t → CtxSingle v s' t := by
    intro t ht hd h
    obtain ⟨i, hi, hroot⟩ := h.item
    exact
      { item := ⟨i, by rw [hsdX _ hd]; exact hi, fun p hp =>
          hroot p ((hPO i (hOld_span t ht i (by rw [hi]; unfold setSides; split <;> simp)) p).1 hp)⟩
        att := by rw [hedges t ht, hg', hsvX _ hd.le]; exact h.att }
  obtain ⟨above, below, hsp, hCE, hCS, hPW, hLow, hEdg, hBel⟩ := hCT.split
  refine
    { top := ⟨T ++ L ++ top, by rw [hts, htop]; simp only [List.append_assoc],
        { split := ?_
          bot := fun t ht k hk => by
            rw [hsvX k hk.le]
            rcases hmemTop t ht with h | h | h
            · rcases hTvs t h with h' | h'
              · rw [h']; exact hvne k hk
              · exact fun h'' => (hvf _ h').2.2.2.1 k hk.le h''.symm
            · rw [(hLprop t h).1]; exact hvne k hk
            · exact hCT.bot t h k hk
          vitems := fun t ht x hx hm => by
            rw [hg'] at hx
            rcases hmemTop t ht with h | h | h
            · rcases hTspan t h _ hm with ⟨c, hc, hm'⟩ | h'
              · rcases hRspan_cls c hc _ hm' (by show 1 + x < _; omega) with ⟨x', hx', h'⟩ | ⟨e', -, h'⟩ | h'
                · exact .inr (hCVL x ((vertItem_inj' h') ▸ hx'))
                · exact absurd h' (vertItem_ne_edgeItem' hx)
                · exact absurd h' (vertItem_ne_edgeItem' hx)
              · exact .inl (vertItem_inj' h')
            · rw [(hLprop t h).2.2.2.1, List.mem_singleton] at hm
              exact .inl (vertItem_inj' hm)
            · exact (hCT.vitems t h x hx hm).imp id (hVL x)
          qitems := fun t ht e' he' hm => by
            rw [hg'] at he' hm
            rcases hmemTop t ht with h | h | h
            · rcases hTspan t h _ hm with ⟨c, hc, hm'⟩ | h'
              · rcases hRspan_cls c hc _ hm' (by show 1 + s.g.nv + e' < _; omega) with
                  ⟨x', hx', h'⟩ | ⟨e'', he'', h'⟩ | h'
                · exact absurd h'.symm (vertItem_ne_edgeItem' (hCv x' hx'))
                · exact ⟨(o, hv₀), hoD, Or.inr ((edgeItem_inj h') ▸ he'')⟩
                · exact ⟨(o, hv₀), hoD, Or.inl (edgeItem_inj h')⟩
              · exact absurd h'.symm (vertItem_ne_edgeItem' hv)
            · rw [(hLprop t h).2.2.2.1, List.mem_singleton] at hm
              exact absurd hm.symm (vertItem_ne_edgeItem' hv)
            · obtain ⟨o', ho', hs⟩ := hCT.qitems t h e' he' hm
              exact ⟨o', List.mem_append_left _ ho', hs⟩
          edges := fun t ht e' he' hte => by
            rw [hg'] at he' hte
            rcases hmemTop t ht with h | h | h
            · rcases hTedge t h e' hte with h' | h'
              · exact ⟨(o, hv₀), hoD, (hRE e' he').1 h'⟩
              · obtain ⟨o', ho', -, hs⟩ := (hVE e' he').1 h'
                exact ⟨o', List.mem_append_left _ ho', hs⟩
            · obtain ⟨o', ho', -, hs⟩ := (hVE e' he').1 ((hLe t h e').1 hte)
              exact ⟨o', List.mem_append_left _ ho', hs⟩
            · rw [hedges₀ t (hmemT t h)] at hte
              obtain ⟨o', ho', hs⟩ := hCT.edges t h e' he' hte
              exact ⟨o', List.mem_append_left _ ho', hs⟩
          cover := fun o' ho' hl e' he' hs => by
            rw [hg'] at he' ⊢
            rcases List.mem_append.1 ho' with h | h
            · obtain ⟨t, ht, hte⟩ := hCT.cover o' h hl e' he' hs
              exact ⟨t, List.mem_append_right _ ht,
                by rw [hedges₀ t (hmemT t ht)]; exact hte⟩
            · rw [List.mem_singleton] at h; subst h
              obtain ⟨c, hc, hce⟩ := (hRE e' he').2 hs
              obtain ⟨t, ht, hte⟩ := hTcov c hc
              exact ⟨t, List.mem_append_left _ (List.mem_append_left _ ht), hte e' hce⟩
          ret := fun _ => ⟨(o, hv₀), hoD, hr⟩ }⟩
      base_edges := ⟨hC.base_edges.1, fun k hk e' he' => by
        rw [hg'] at he'
        have hbm : base[k]! ∈ s.tstack := by
          rw [htop, getElem!_pos base k hk]; exact List.mem_append_right _ (List.getElem_mem hk)
        rw [hedges _ hbm]
        exact hC.base_edges.2 k hk e' he'⟩
      base_bot := hC.base_bot
      sv := fun k hk => by rw [hsvX k hk]; exact hC.sv k hk
      sd := fun k hk => by rw [hsdX k hk]; exact hC.sd k hk
      sv_d := by rw [hsvX d (Nat.le_refl _)]; exact hvd
      path := fun k k' h h' => by rw [hsvX k (Nat.le_of_lt (Nat.lt_of_lt_of_le h h')), hsvX k' h']; exact hC.path k k' h h'
      v_root := fun p h => hC.v_root p ((hPO _ hOld_v p).1 h)
      vert_free := fun _ _ _ => rfl
      afterVert_ret := fun o' ho' => by
        rcases (hAVm o').1 ho' with h | ⟨-, rfl⟩
        · exact hC.afterVert_ret o' h
        · exact hr
      hv_ret := fun _ => ⟨(o, hv₀), hoD, hr⟩
      noVert_after := fun h => nomatch h
      vfirst := fun t ht hvt hnv => by
        rcases hmemS t ht with h | h | h
        · rcases hTmem t h with h₁ | ⟨-, -, -, h₁⟩ | ⟨-, c, hc, -, -, -, h₁, -⟩
          · cases hh : hv₀ with
            | true =>
              obtain ⟨x, hRx, -, -, -, -, hxf2, -⟩ := hRvert hh
              subst hRx
              rw [List.mem_singleton] at h₁; subst h₁
              exact hxf2
            | false => exact absurd hvt ((hRnoV hh).2 t h₁).1
          · exact absurd (by rw [h₁]; exact List.mem_singleton_self _) hnv
          · exact absurd ((h₁ _).2 (.inl rfl)) hnv
        · exact absurd (by rw [(hLprop t h).2.2.2.1]; exact List.mem_singleton_self _) hnv
        · exact Nat.lt_of_lt_of_le (hC.vfirst t h hvt hnv) hRnx
      vert_book := fun h => nomatch h
      vert_disj := fun h => nomatch h
      vert_edges := fun e' he' => by
        rw [hg'] at he' ⊢
        rw [hVE e' he']
        constructor
        · rintro ⟨o', ho', hl, hs⟩; exact ⟨o', List.mem_append_left _ ho', hl, hs⟩
        · rintro ⟨o', ho', hl, hs⟩
          rcases List.mem_append.1 ho' with h | h
          · exact ⟨o', h, hl, hs⟩
          · rw [List.mem_singleton] at h; subst h; exact absurd hl (Nat.not_le.2 hr)
      touch_bot := fun t ht hne => by
        rw [hg'] at hne ⊢
        rcases hmemS t ht with h | h | h
        · rcases hTmem t h with h₁ | ⟨-, hvs, -, h₁⟩ | ⟨-, c, hc, hce, hvs, -, -, h₁⟩
          · exact hRtb t h₁ hne
          · rw [hvs]
            obtain ⟨e', he', hte⟩ := hne
            obtain ⟨e₁, he₁, hb₁, hI₁⟩ := hTVv e' he' ((hUe t h₁ e').1 hte)
            exact ⟨e₁, he₁, (hUe t h₁ e₁).2 hb₁, hI₁⟩
          · rw [hvs]
            obtain ⟨e₁, he₁, hb₁, hI₁⟩ := hRtb c hc ⟨e, he, hce⟩
            exact ⟨e₁, he₁, (h₁ e₁).2 (.inr hb₁), hI₁⟩
        · rw [(hLprop t h).1]
          obtain ⟨e', he', hte⟩ := hne
          obtain ⟨e₁, he₁, hb₁, hI₁⟩ := hTVv e' he' ((hLe t h e').1 hte)
          exact ⟨e₁, he₁, (hLe t h e₁).2 hb₁, hI₁⟩
        · rw [hedges₀ t h] at hne ⊢
          exact hC.touch_bot t h hne
      span_root := fun t ht i hi p hp => by
        rcases hmemS t ht with h | h | h
        · rcases hTspan t h i hi with ⟨c, hc, hi'⟩ | rfl
          · exact hRroot c hc i hi' p hp
          · exact hC.v_root p ((hPO _ hOld_v p).1 hp)
        · rw [(hLprop t h).2.2.2.1, List.mem_singleton] at hi; subst hi
          exact hC.v_root p ((hPO _ hOld_v p).1 hp)
        · exact hC.span_root t h i hi p ((hPO i (hOld_span t h i hi) p).1 hp)
      span_lt := fun t ht i hi => by
        rcases hmemS t ht with h | h | h
        · rcases hTspan t h i hi with ⟨c, hc, hi'⟩ | rfl
          · exact hRlt c hc i hi'
          · exact Nat.lt_of_lt_of_le hOld_v.1 hRsz
        · rw [(hLprop t h).2.2.2.1, List.mem_singleton] at hi; subst hi
          exact Nat.lt_of_lt_of_le hOld_v.1 hRsz
        · exact Nat.lt_of_lt_of_le (hC.span_lt t h i hi) hRsz
      ch_lt := hRch
      disj := by
        rw [hts, hg']
        refine List.pairwise_append.2 ⟨hTpwE, ?_, fun t ht t' ht' e' he' h1 h2 => ?_⟩
        · refine List.pairwise_append.2 ⟨?_, ?_, fun t ht t' ht' e' he' h1 h2 => ?_⟩
          · rw [H.hL]; split <;> simp
          · refine hC.disj.imp_of_mem fun {t t'} ht ht' h e' he' h1 h2 => ?_
            rw [hedges₀ t ht] at h1
            rw [hedges₀ t' ht'] at h2
            exact h e' he' h1 h2
          · rw [hedges₀ t' ht'] at h2
            exact hC.vert_disj (hLprop t ht).2.2.2.2.1 t' ht' e' he' h2 ((hEBv e').1 ((hLe t ht e').1 h1))
        · rcases hTedge t ht e' h1 with hR | hV
          · rcases List.mem_append.1 ht' with h' | h'
            · exact hRVdisj e' he' hR ((hLe t' h' e').1 h2)
            · rw [hedges₀ t' h'] at h2
              exact hRolddisj e' he' hR t' h' h2
          · have h0 : hv₀ = false := by
              rcases hTmem t ht with h₁ | ⟨h₁, -⟩ | ⟨h₁, -⟩
              · exact absurd hV (hRVdisj e' he' ⟨t, h₁, h1⟩)
              · exact h₁
              · exact h₁
            rcases List.mem_append.1 ht' with h' | h'
            · exact absurd ((hLprop t' h').2.2.2.2.2.symm.trans (hpush0 h0)) (by decide)
            · rw [hedges₀ t' h'] at h2
              exact hC.vert_disj (hhv0 h0) t' h' e' he' h2 ((hEBv e').1 hV)
      span_disj := by
        rw [hts]
        refine List.pairwise_append.2 ⟨hTpwS, ?_, fun t ht t' ht' i hi hi' => ?_⟩
        · refine List.pairwise_append.2 ⟨by rw [H.hL]; split <;> simp, hC.span_disj,
            fun t ht t' ht' i hi hi' => ?_⟩
          rw [(hLprop t ht).2.2.2.1, List.mem_singleton] at hi; subst hi
          exact Bool.false_ne_true ((hLprop t ht).2.2.2.2.1.symm.trans (hC.vert_free t' ht' hi'))
        · rcases hTspan t ht i hi with ⟨c, hc, hic⟩ | rfl
          · exact (hRnew c hc i hic).1 t' ht' hi'
          · have h0 : hv₀ = false := by
              rcases hTmem t ht with h₁ | ⟨h₁, -⟩ | ⟨h₁, -⟩
              · exact absurd hi (hRnoVv t h₁)
              · exact h₁
              · exact h₁
            rcases List.mem_append.1 ht' with h' | h'
            · exact absurd ((hLprop t' h').2.2.2.2.2.symm.trans (hpush0 h0)) (by decide)
            · exact Bool.false_ne_true ((hhv0 h0).symm.trans (hC.vert_free t' h' hi'))
      q_fresh := fun o' ho' e' hs => by
        obtain ⟨hns, he'⟩ := H.rest_nd o' ho' e' hs
        obtain ⟨h1, h2, h3⟩ := hC.q_fresh o' (List.mem_cons_of_mem _ ho') e' hs
        have hOq : Old (edgeItem s.g e') :=
          ⟨by show 1 + s.g.nv + e' < _; omega, fun x hx h => vertItem_ne_edgeItem' (hCv x hx) h.symm,
            fun e'' he'' h => hns (Or.inr ((edgeItem_inj h) ▸ he'')),
            fun h => hns (Or.inl (edgeItem_inj h))⟩
        rw [hg']
        refine ⟨fun p hp => h1 p ((hPO _ hOq p).1 hp), by rw [hkept _ hOq]; exact h2, fun t ht hm => ?_⟩
        rcases hmemS t ht with h | h | h
        · rcases hTspan t h _ hm with ⟨c, hc, hm'⟩ | h'
          · rcases hRspan_cls c hc _ hm' hOq.1 with ⟨x', hx', h'⟩ | ⟨e'', he'', h'⟩ | h'
            · exact vertItem_ne_edgeItem' (hCv x' hx') h'.symm
            · exact hns (Or.inr ((edgeItem_inj h') ▸ he''))
            · exact hns (Or.inl (edgeItem_inj h'))
          · exact vertItem_ne_edgeItem' hv h'.symm
        · rw [(hLprop t h).2.2.2.1, List.mem_singleton] at hm
          exact vertItem_ne_edgeItem' hv hm.symm
        · exact h3 t h hm
      v_fresh := fun o' ho' e₁ cls₁ child hoc w hw => by
        obtain ⟨h1, h2, h3, h4, h5, h6⟩ :=
          hC.v_fresh o' (List.mem_cons_of_mem _ ho') e₁ cls₁ child hoc w hw
        obtain ⟨hwC, hwv⟩ := H.rest_nv o' ho' e₁ cls₁ child hoc w hw
        have hwne : v ≠ w := fun h => h4 d (Nat.le_refl _) (hvd.trans h)
        have hOw : Old (vertItem w) :=
          ⟨by show 1 + w < _; omega, fun x hx h => hwC ((vertItem_inj' h) ▸ hx),
            fun _ _ => vertItem_ne_edgeItem' hwv, vertItem_ne_edgeItem' hwv⟩
        refine ⟨fun p hp => h1 p ((hPO _ hOw p).1 hp), by rw [hkept _ hOw]; exact h2,
          fun t ht hm => ?_, fun k hk => by rw [hsvX k hk]; exact h4 k hk, fun t ht => ?_,
          fun t ht hT => ?_⟩
        · rcases hmemS t ht with h | h | h
          · rcases hTspan t h _ hm with ⟨c, hc, hm'⟩ | h'
            · rcases hRspan_cls c hc _ hm' hOw.1 with ⟨x', hx', h'⟩ | ⟨e'', -, h'⟩ | h'
              · exact hwC ((vertItem_inj' h') ▸ hx')
              · exact vertItem_ne_edgeItem' hwv h'
              · exact vertItem_ne_edgeItem' hwv h'
            · exact hwne (vertItem_inj' h').symm
          · rw [(hLprop t h).2.2.2.1, List.mem_singleton] at hm
            exact hwne (vertItem_inj' hm).symm
          · exact h3 t h hm
        · rcases hmemS t ht with h | h | h
          · rcases hTvs t h with h' | h'
            · rw [h']; exact hwne
            · exact fun hh => hwC (hh ▸ h')
          · rw [(hLprop t h).1]; exact hwne
          · exact h5 t h
        · rw [hg'] at hT
          rcases hmemS t ht with h | h | h
          · obtain ⟨e', he', hte, hI⟩ := hT
            rcases hTedge t h e' hte with ⟨c, hc, hce⟩ | hV
            · rcases hRtouch c hc w ⟨e', he', hce, hI⟩ with hh | hh | ⟨k, -, hk, hh⟩
              · exact hwne hh.symm
              · exact hwC hh
              · exact h4 k hk hh.symm
            · obtain ⟨o'', ho'', e'', -, hs, hx⟩ := hTV w ⟨e', he', hV, hI⟩
              exact hnc o'' ho'' e'' hs w hx o' ho' e₁ cls₁ child hoc hw
          · obtain ⟨e', he', hte, hI⟩ := hT
            obtain ⟨o'', ho'', e'', -, hs, hx⟩ := hTV w ⟨e', he', (hLe t h e').1 hte, hI⟩
            exact hnc o'' ho'' e'' hs w hx o' ho' e₁ cls₁ child hoc hw
          · rw [hedges₀ t h] at hT
            exact h6 t h hT }
  -- split
  simp only [↓reduceIte]
  cases hv₀ with
  | true =>
    obtain ⟨x, hRx, hxv, hxd, hxside, hxf1, hxf2, hxtt, hxatt, hxt1⟩ := hRvert rfl
    have hTR : T = R := by
      rcases hT with ⟨-, h⟩ | ⟨h, -⟩
      · exact h
      · exact absurd h (by decide)
    subst hTR
    subst hRx
    have hxm : x ∈ [x] := List.mem_singleton_self _
    have hxe : TEntry.edges s.g s'.items x e := by
      obtain ⟨t, ht, h⟩ := (hRE e he).2 hsube
      rw [List.mem_singleton] at ht; subst ht; exact h
    have hxCE : CtxEntry v d s' x :=
      { vStart := hxv
        depth := hxd ▸ hr
        nonempty := ⟨e, by rw [hg']; exact he, by rw [hg']; exact hxe⟩
        touch_bot := by rw [hg', ← hxv]; exact hRtb x hxm ⟨e, he, hxe⟩
        touch_top := by rw [hg', hsvX _ (hxd ▸ hr.le), hxd]; exact hxtt
        side := by rw [hxd, hsdX _ hr]; exact hxside
        att := fun w hw hwv hni => by
          rw [hg'] at hw hni
          obtain ⟨k, hk1, hk2, hk3⟩ := hxatt w hw hwv hni
          exact ⟨k, hxd ▸ hk1, hk2, by rw [hsvX k hk2]; exact hk3⟩ }
    have hxCS : allType1 d x.topDepth (done ++ [(o, true)]) → CtxSingle v s' x := by
      intro hA
      have ht1 : cls.isType1 = true := hA o ((hAVm o).2 (.inr ⟨rfl, rfl⟩)) hxd.symm
      obtain ⟨⟨i, hi⟩, hatt⟩ := hxt1 ht1
      exact
        { item := ⟨i, by rw [hxd, hsdX _ hr]; exact hi, fun p hp =>
            hRroot x hxm i (by rw [hi]; unfold setSides; split <;> simp) p hp⟩
          att := fun w hw => by
            rw [hg'] at hw ⊢
            rcases hatt w hw with h | h | h
            · exact .inl h
            · exact .inr (.inl (by rw [hxd, hsvX _ hr.le]; exact h))
            · exact .inr (.inr h) }
    have habm : ∀ t ∈ above, t ∈ top := by
      intro t ht
      by_cases hhv : hasVert = true
      · rw [if_pos hhv] at hsp
        obtain ⟨vt, htv, -⟩ := hsp
        rw [htv]; exact List.mem_append_left _ ht
      · rw [if_neg hhv] at hsp
        exact absurd ht (hsp.1 ▸ List.not_mem_nil)
    have habove : ∀ t' ∈ above, t'.topDepth ≤ x.topDepth ∧ t'.firstIdx < x.firstIdx := by
      intro t' ht'
      have hhv : hasVert = true := by
        cases h : hasVert
        · exact absurd ht' (by rw [if_neg (by rw [h]; decide)] at hsp; exact hsp.1 ▸ List.not_mem_nil)
        · rfl
      rw [if_pos hhv] at hsp
      obtain ⟨vt, htv, -, -, -, h4, -⟩ := hsp
      obtain ⟨o', ho', hto⟩ := hLow t' ht'
      refine ⟨?_, ?_⟩
      · rw [hxd, hto]
        exact lowval_le_of_rank (hC.afterVert_ret o' ho') hr (H.rank _ (mem_afterVert ho'))
      · exact Nat.lt_of_lt_of_le (hC.vfirst t' (hmemT t' (by rw [htv]; exact List.mem_append_left _ ht'))
          (hCE t' ht').vStart (h4 t' ht')) hxf1
    refine ⟨x :: above, below, ?_, ?_, ?_, ?_, ?_, ?_, hBel⟩
    · by_cases hhv : hasVert = true
      · rw [if_pos hhv] at hsp
        obtain ⟨vt, htv, h1, h2, h3, h4, h5, h6⟩ := hsp
        have hpf : push = false := by
          cases hp : push
          · rfl
          · exact absurd (H.hpush.1 hp).1 (by rw [hhv]; decide)
        have hL0 : L = [] := by rw [H.hL, hpf]; rfl
        refine ⟨vt, by rw [hL0, htv]; rfl, h1, h2, h3, fun t ht => ?_, fun o' ho' => ?_,
          fun o' ho' e' he' hs => ?_⟩
        · rcases List.mem_cons.1 ht with rfl | ht
          · exact hRnoVv t hxm
          · exact h4 t ht
        · rcases (hAVm o').1 ho' with ho' | ⟨-, rfl⟩
          · exact (h5 o' ho').imp (fun ⟨t, ht, h⟩ => ⟨t, List.mem_cons_of_mem _ ht, h⟩) id
          · exact .inl ⟨x, List.mem_cons_self .., hxd⟩
        · rw [hg'] at he'
          rcases (hAVm o').1 ho' with ho' | ⟨-, rfl⟩
          · exact (h6 o' ho' e' he' hs).imp
              (fun ⟨t, ht, h⟩ => ⟨t, List.mem_cons_of_mem _ ht,
                by rw [hedges t (hmemT t (by rw [htv]; exact List.mem_append_left _ ht))]; exact h⟩)
              (by rw [hedges vt (hmemT vt (by rw [htv]; exact List.mem_append_right _ (List.mem_cons_self ..)))]; exact id)
          · obtain ⟨t, ht, h⟩ := (hRE e' he').2 hs
            rw [List.mem_singleton] at ht; subst ht
            exact .inl ⟨t, List.mem_cons_self .., by rw [hg']; exact h⟩
      · have hhf : hasVert = false := by
          cases h : hasVert with
          | false => rfl
          | true => exact absurd h hhv
        rw [if_neg (by rw [hhf]; decide)] at hsp
        obtain ⟨hab, htb⟩ := hsp
        have hpt : push = true := by
          cases hp : push
          · exact absurd hb (by rw [hhf, hp]; decide)
          · rfl
        obtain ⟨u₀, hL1⟩ : ∃ u₀, L = [u₀] := by rw [H.hL, hpt]; exact ⟨_, rfl⟩
        have hu := hLprop u₀ (by rw [hL1]; exact List.mem_singleton_self _)
        have hAV0 := hC.noVert_after hhf
        refine ⟨u₀, by rw [hL1, hab, htb]; rfl, hu.2.1 ▸ Nat.le_refl _,
          by rw [hu.2.2.2.1]; exact List.mem_singleton_self _, fun _ => ⟨hu.2.1, hu.2.2.2.1⟩,
          fun t ht => ?_, fun o' ho' => ?_, fun o' ho' e' he' hs => ?_⟩
        · rcases List.mem_cons.1 ht with rfl | ht
          · exact hRnoVv t hxm
          · exact absurd ht (hab ▸ List.not_mem_nil)
        · rcases (hAVm o').1 ho' with ho' | ⟨-, rfl⟩
          · exact absurd ho' (hAV0 ▸ List.not_mem_nil)
          · exact .inl ⟨x, List.mem_cons_self .., hxd⟩
        · rw [hg'] at he'
          rcases (hAVm o').1 ho' with ho' | ⟨-, rfl⟩
          · exact absurd ho' (hAV0 ▸ List.not_mem_nil)
          · obtain ⟨t, ht, h⟩ := (hRE e' he').2 hs
            rw [List.mem_singleton] at ht; subst ht
            exact .inl ⟨t, List.mem_cons_self .., by rw [hg']; exact h⟩
    · intro t ht
      rcases List.mem_cons.1 ht with rfl | ht
      · exact hxCE
      · exact hCE' t (hmemT t (habm t ht)) (hCE t ht)
    · intro t ht hA
      rcases List.mem_cons.1 ht with rfl | ht
      · exact hxCS hA
      · exact hCS' t (hmemT t (habm t ht)) (hCE t ht).depth
          (hCS t ht fun o' ho' hl => hA o' ((hAVm o').2 (.inl ho')) hl)
    · exact List.pairwise_cons.2 ⟨habove, hPW⟩
    · intro t ht
      rcases List.mem_cons.1 ht with rfl | ht
      · exact ⟨o, (hAVm o).2 (.inr ⟨rfl, rfl⟩), hxd⟩
      · obtain ⟨o', ho', h⟩ := hLow t ht
        exact ⟨o', (hAVm o').2 (.inl ho'), h⟩
    · intro t ht e' he' hte
      rw [hg'] at he' hte
      rcases List.mem_cons.1 ht with rfl | ht
      · exact ⟨o, (hAVm o).2 (.inr ⟨rfl, rfl⟩), (hRE e' he').1 ⟨t, hxm, hte⟩⟩
      · rw [hedges₀ t (hmemT t (habm t ht))] at hte
        obtain ⟨o', ho', h⟩ := hEdg t ht e' he' hte
        exact ⟨o', (hAVm o').2 (.inl ho'), h⟩
  | false =>
    have hhf := hhv0 rfl
    have hL0' := hL0 rfl
    have hAV0 := hC.noVert_after hhf
    rw [if_neg (by rw [hhf]; decide)] at hsp
    obtain ⟨hab, htb⟩ := hsp
    have hAV' : afterVert (done ++ [(o, false)]) = [] := by
      rw [List.eq_nil_iff_forall_not_mem]
      intro o' ho'
      rcases (hAVm o').1 ho' with h | ⟨h, -⟩
      · exact (hAV0 ▸ List.not_mem_nil) h
      · exact Bool.false_ne_true h
    have hRnv : ∀ t ∈ R, t.vStart ≠ v := fun t ht => ((hRnoV rfl).2 t ht).1
    rcases hT with ⟨h, -⟩ | ⟨-, f, dir, rfl | ⟨c, R', hR, rfl⟩⟩
    · exact absurd h (by decide)
    · refine ⟨[], R ++ top, ⟨⟨v, d, f, setSides dir [vertItem v] []⟩, by rw [hL0']; simp,
        Nat.le_refl _, by rw [hus]; exact List.mem_singleton_self _, fun _ => ⟨rfl, hus f dir⟩,
        fun t ht => absurd ht List.not_mem_nil, fun o' ho' => absurd ho' (hAV' ▸ List.not_mem_nil),
        fun o' ho' => absurd ho' (hAV' ▸ List.not_mem_nil)⟩,
        fun t ht => absurd ht List.not_mem_nil, fun t ht => absurd ht List.not_mem_nil,
        List.Pairwise.nil, fun t ht => absurd ht List.not_mem_nil,
        fun t ht => absurd ht List.not_mem_nil, fun t ht => ?_⟩
      rcases List.mem_append.1 ht with h | h
      · exact hRnv t h
      · exact hBel t (htb ▸ h)
    · subst hR
      have hRnv' := hRnv
      refine ⟨[], R' ++ top, ⟨TEntry.mergeInto ⟨v, d, f, setSides dir [vertItem v] []⟩ c,
        by rw [hL0']; simp, Nat.min_le_right _ _,
        (mem_mergeInto_spans _ _ _).2 (.inl (by rw [hus]; exact List.mem_singleton_self _)),
        fun h => absurd h (hRnv' c (List.mem_cons_self ..)),
        fun t ht => absurd ht List.not_mem_nil, fun o' ho' => absurd ho' (hAV' ▸ List.not_mem_nil),
        fun o' ho' => absurd ho' (hAV' ▸ List.not_mem_nil)⟩,
        fun t ht => absurd ht List.not_mem_nil, fun t ht => absurd ht List.not_mem_nil,
        List.Pairwise.nil, fun t ht => absurd ht List.not_mem_nil,
        fun t ht => absurd ht List.not_mem_nil, fun t ht => ?_⟩
      rcases List.mem_append.1 ht with h | h
      · exact hRnv' t (List.mem_cons_of_mem _ h)
      · exact hBel t (htb ▸ h)

theorem ctxShape_ret {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s : WalkState} {e : Nat} {cls : OutClass} {y : Nat} {outs : List DfsOut} {L : List TEntry}
    {push : Bool} {done' : List (DfsOut × Bool)} {hv' : Bool} {bE' : List (Nat → Prop)}
    {sv' : List Nat} {sd' : List Bool} {sE : WalkState} {dir' : Bool} {L' : List TEntry}
    {push' : Bool} {D₃ : Array Bool}
    (H : TreeSite v d done rest hasVert base bE sv sd s e cls y outs L push done' hv' bE' sv' sd'
      sE dir' L' push' D₃)
    (hS : CtxShape v d done hasVert base s)
    (hr : cls.lowval d < d) {hv₀ : Bool}
    (hb : (hasVert || push) = hv₀) {sX : WalkState} {R : List TEntry}
    (HR : RetTop v d s e cls y outs L hv₀ sX R) (T : List TEntry) {s' : WalkState}
    (hT : (hv₀ = true ∧ T = R) ∨ (hv₀ = false ∧ ∃ f,
      T = ⟨v, d, f, setSides (!s.stackDir[cls.lowval d]!) [vertItem v] []⟩ :: R ∨
      ∃ c R', R = c :: R' ∧
        T = TEntry.mergeInto ⟨v, d, f, setSides (!s.stackDir[cls.lowval d]!) [vertItem v] []⟩ c
          :: R'))
    (hnoP : hasVert = true → cls.isType1 = true → ∀ a tl, s.tstack = a :: tl → a.vStart = v →
      a.topDepth ≠ cls.lowval d)
    (hts : s'.tstack = T ++ (L ++ s.tstack))
    (hsd : ∀ k, k < d → s'.stackDir[k]! = sX.stackDir[k]!) :
    CtxShape v d (done ++ [(.tree e cls (.node y outs), hv₀)]) true base s' := by
  have hC := H.ctx
  have hv : v < s.g.nv := H.hv
  set lv := cls.lowval d with hlv
  have hsdX : ∀ k, k < d → s'.stackDir[k]! = s.stackDir[k]! := fun k hk =>
    (hsd k hk).trans (HR.sd k hk)
  have ht1hv : cls.isType1 = true → hv₀ = true := fun ht1 => by
    cases hhv : hasVert
    · have hp : push = true := H.hpush.2 ⟨hhv, hr, ht1⟩
      rw [← hb, hhv, hp]; rfl
    · rw [← hb, hhv]; rfl
  have hVlt : vertItem v < s.items.size := by
    have := H.hsz; show 1 + v < _; omega
  have hVR : ∀ t ∈ R, vertItem v ∉ t.spans.1 ++ t.spans.2 := fun t ht hm => by
    rcases (HR.span_new t ht _ hm).2 hVlt with ⟨x, hx, hxe⟩ | ⟨e', -, he'⟩ | he'
    · exact H.v_nc (by rw [vertItem_inj' hxe]; exact hx)
    · exact vertItem_ne_edgeItem' hv he'
    · exact vertItem_ne_edgeItem' hv he'
  have hLmem : ∀ t ∈ L, t.vStart = v ∧ t.spans = setSides (!s.stackDir[lv]!) [vertItem v] [] := by
    intro t ht
    rw [H.hL] at ht
    split at ht
    · rw [List.mem_singleton] at ht; subst ht
      refine ⟨rfl, ?_⟩
      show setSides _ [vertItem v] [] = _
      rw [getElem!_set!_self' _ _ _ H.hsd, ite_eq_right (Nat.not_le.2 hr)]
    · exact absurd ht List.not_mem_nil
  obtain ⟨top, htop, hCT⟩ := hC.top
  refine
    { ret_hv := fun _ => rfl
      t1_flag := fun o ho hl ht => by
        rcases List.mem_append.1 ho with h | h
        · exact hS.t1_flag o h hl ht
        · rw [List.mem_singleton] at h; subst h; exact ht1hv ht
      t1 := fun _ hT1 => ?_
      vt_side := fun top' htop' t ht htv hm => ?_ }
  · have ht1 : cls.isType1 = true :=
      hT1 _ (List.mem_append_right _ (List.mem_singleton_self _)) hr
    have hv₀t : hv₀ = true := ht1hv ht1
    have hTR : T = R := by
      rcases hT with ⟨-, h⟩ | ⟨h, -⟩
      · exact h
      · rw [hv₀t] at h; exact absurd h (by decide)
    obtain ⟨x, hRx, hxv, hxd, -⟩ := HR.vert hv₀t
    rw [hTR, hRx] at hts
    cases hp : push
    · have hhv : hasVert = true := by
        cases h : hasVert
        · exact absurd (H.hpush.2 ⟨h, hr, ht1⟩) (by rw [hp]; decide)
        · rfl
      have hL0 : L = [] := by rw [H.hL, hp, ite_eq_right Bool.false_ne_true]
      rw [hL0, List.nil_append] at hts
      obtain ⟨above, vt, hab, hvtv, hvtm, hpw⟩ :=
        hS.t1 hhv fun o ho hl => hT1 o (List.mem_append_left _ ho) hl
      refine ⟨x :: above, vt, by rw [hts, hab]; rfl, hvtv, hvtm, List.pairwise_cons.2 ⟨?_, hpw⟩⟩
      obtain ⟨above₀, below₀, hsp, hCE, -, -, hLow, -, -⟩ := hCT.split
      simp only [hhv, ↓reduceIte] at hsp
      obtain ⟨vt₀, htv₀, -, hvtm₀, -, -, -, -⟩ := hsp
      have hab₀ : s.tstack = above₀ ++ vt₀ :: (below₀ ++ base) := by rw [htop, htv₀]; simp
      obtain ⟨rfl, rfl, -⟩ := split_unique hab hab₀ hvtm hvtm₀ hC.span_disj
      have hle : ∀ t ∈ above, t.topDepth ≤ lv := fun t ht => by
        obtain ⟨o, ho, hto⟩ := hLow t ht
        have hrk := H.rank (o, true) (mem_afterVert ho)
        rw [hto]
        exact lowval_le_of_rank (by rw [← hto]; exact (hCE t ht).depth) hr hrk
      rw [hxd]
      rcases above with _ | ⟨a, above'⟩
      · intro t ht; exact absurd ht List.not_mem_nil
      · have ha : a.topDepth < lv := Nat.lt_of_le_of_ne (hle a (List.mem_cons_self ..))
          (hnoP hhv ht1 a _ (by rw [hab]; rfl) (hCE a (List.mem_cons_self ..)).vStart)
        intro t ht
        rcases List.mem_cons.1 ht with rfl | ht'
        · exact ha
        · exact Nat.lt_trans ((List.pairwise_cons.1 hpw).1 t ht') ha
    · have hhv : hasVert = false := (H.hpush.1 hp).1
      have htop0 : top = [] := by
        by_contra h
        exact absurd (hS.ret_hv (hCT.ret h)) (by rw [hhv]; decide)
      rw [H.hL, hp, ite_eq_left rfl] at hts
      refine ⟨[x], ⟨v, d, s.nxtEdgeIdx, setSides (s.stackDir.set! d (if lv ≥ d then false
        else !s.stackDir[lv]!))[d]! [vertItem v] []⟩, ?_, rfl, ?_, List.pairwise_singleton _ _⟩
      · rw [hts, htop, htop0]; rfl
      · show vertItem v ∈ (setSides _ [vertItem v] []).1 ++ (setSides _ [vertItem v] []).2
        rw [mem_setSides]; exact List.mem_singleton_self _
  · have htop'' : top' = T ++ (L ++ top) := by
      apply List.append_cancel_right (bs := base)
      rw [← htop', hts, htop]; simp
    rw [htop''] at ht
    rcases List.mem_append.1 ht with ht | ht
    · rcases hT with ⟨hv₀t, rfl⟩ | ⟨hv₀f, f, hT | ⟨c, R', hRc, hT⟩⟩
      · exact absurd hm (hVR t ht)
      · subst hT
        rcases List.mem_cons.1 ht with rfl | ht
        · exact ⟨(.tree e cls (.node y outs), hv₀),
            List.mem_append_right _ (List.mem_singleton_self _), hr, by
              show setSides (!s.stackDir[lv]!) [vertItem v] [] =
                setSides (!s'.stackDir[lv]!) [vertItem v] []
              rw [hsdX _ hr]⟩
        · exact absurd htv ((HR.noVert hv₀f).2 t ht).1
      · subst hT
        rcases List.mem_cons.1 ht with rfl | ht
        · exact absurd htv ((HR.noVert hv₀f).2 c (by rw [hRc]; exact List.mem_cons_self ..)).1
        · exact absurd htv ((HR.noVert hv₀f).2 t (by rw [hRc]; exact List.mem_cons_of_mem _ ht)).1
    · rcases List.mem_append.1 ht with ht | ht
      · obtain ⟨-, hsp⟩ := hLmem t ht
        exact ⟨(.tree e cls (.node y outs), hv₀),
          List.mem_append_right _ (List.mem_singleton_self _), hr, by
            show _ = setSides (!s'.stackDir[lv]!) [vertItem v] []
            rw [hsp, hsdX _ hr]⟩
      · obtain ⟨o, ho, hl, hsp⟩ := hS.vt_side top htop t ht htv hm
        exact ⟨o, List.mem_append_left _ ho, hl, by rw [hsdX _ hl]; exact hsp⟩


/-- The P-merge at a returning tree edge with the vertex entry (`condP` fires): `condP` forces
`push = false` (a pushed `V v` entry would sit at depth `d ≠ lowval`), hence `hasVert = true` and
`L = []`, and the entry below the retained `(v, lowval)` entry is the enclosing parallel ear.
`earCtx_ret` (with `T = R`) gives the context at `sX`, `mergeP_single` the merged side,
`earCtx_mergeP` the post-state. -/
theorem ctx_step_tree_ret_P {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s : WalkState} {e : Nat} {cls : OutClass} {y : Nat} {outs : List DfsOut} {L : List TEntry}
    {push : Bool} {done' : List (DfsOut × Bool)} {hv' : Bool} {bE' : List (Nat → Prop)}
    {sv' : List Nat} {sd' : List Bool} {sE : WalkState} {dir' : Bool} {L' : List TEntry}
    {push' : Bool} {D₃ : Array Bool}
    (H : TreeSite v d done rest hasVert base bE sv sd s e cls y outs L push done' hv' bE' sv' sd'
      sE dir' L' push' D₃)
    (hr : cls.lowval d < d) {sX : WalkState} {R : List TEntry}
    (HR : RetTop v d s e cls y outs L true sX R) (b : Bool)
    (hc : result (condP v (cls.lowval d) cls.isType1) sX = true) :
    wp (finishRest v d (cls.lowval d) cls.isType1 true b)
      (fun hv'' s' => EarCtx v d (done ++ [(.tree e cls (.node y outs), true)]) rest hv''
        base bE sv sd s' ∧ OutFrame v (.tree e cls (.node y outs)) base s s' ∧
        (CtxShape v d done hasVert base s →
          CtxShape v d (done ++ [(.tree e cls (.node y outs), true)]) hv'' base s')) sX := by
  have hc' : (cls.isType1 && decide (sX.tstack.length ≥ 2) && (sX.tstack.tail.head!.vStart == v) &&
      (sX.tstack.tail.head!.topDepth == cls.lowval d)) = true := by
    rw [result, run_condP] at hc; exact hc
  have hc2 := hc'
  simp only [Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq] at hc2
  obtain ⟨⟨⟨ht1, hlen⟩, hnv⟩, hnd⟩ := hc2
  set lv := cls.lowval d with hlv
  obtain ⟨x, hR, hxv, hxd, -, -, -, -, -, hxt1⟩ := HR.vert rfl
  obtain ⟨⟨q, hxq⟩, -⟩ := hxt1 ht1
  have hts₀ := HR.tstack
  rw [hR] at hts₀
  have hL0 : L = [] := by
    have hL := H.hL
    cases hp : push
    · rw [hp] at hL; exact hL
    · exfalso
      simp only [hp, ↓reduceIte] at hL
      rw [hL] at hts₀
      rw [hts₀] at hnd
      simp only [List.cons_append, List.tail_cons] at hnd
      exact absurd hnd (Nat.ne_of_gt hr)
  have hp0 : push = false := by
    cases hp : push
    · rfl
    · exfalso
      have hL := H.hL
      simp [hp, hL0] at hL
  have hb : (hasVert || push) = true := by
    cases hhv : hasVert
    · exact absurd (H.hpush.2 ⟨hhv, hr, ht1⟩) (by rw [hp0]; decide)
    · rfl
  have hhvT : hasVert = true := by
    cases hhv : hasVert
    · exact absurd (H.hpush.2 ⟨hhv, hr, ht1⟩) (by rw [hp0]; decide)
    · rfl
  subst hhvT
  rw [hL0, List.nil_append] at hts₀
  cases hst : s.tstack with
  | nil =>
    exfalso
    rw [hts₀, hst] at hlen
    simp at hlen
  | cons a tl =>
  rw [hst] at hts₀
  have hts₂ : sX.tstack = x :: a :: tl := hts₀
  rw [hts₂] at hnv hnd
  simp only [List.tail_cons, List.head!_cons] at hnv hnd
  obtain ⟨tl', rfl⟩ : ∃ tl', tl = tl' ++ base := by
    obtain ⟨top, htop, -⟩ := H.ctx.top
    cases top with
    | nil =>
      exfalso
      rw [hst] at htop
      have hb0 : base = a :: tl := by simpa using htop.symm
      exact H.ctx.base_bot a (by rw [hb0]; exact List.mem_cons_self ..) hnv
    | cons a' tl' =>
      rw [hst, List.cons_append, List.cons.injEq] at htop
      exact ⟨tl', htop.2⟩
  have hF0 : OutFrame v (.tree e cls (.node y outs)) base s sX :=
    outFrame_retTop HR (top := a :: tl') (by rw [hst]; rfl) H.L_spans
      (fun t ht i hi => Or.inr ⟨t, ht, hi⟩) HR.tstack rfl rfl
  have hC : EarCtx v d (done ++ [(.tree e cls (.node y outs), true)]) rest true base bE sv sd sX :=
    earCtx_ret H hr hb HR R (Or.inl ⟨rfl, rfl⟩) HR.tstack rfl rfl rfl (fun _ _ => rfl) rfl
  have hsdl : sX.stackDir[lv]! = s.stackDir[lv]! := HR.sd lv hr
  have hcsp : x.spans = setSides sX.stackDir[lv]! [q] [] := by rw [hsdl]; exact hxq
  have hcs : x.spans.1 ++ x.spans.2 = [q] := by rw [hxq]; exact spans_setSides_single _ _
  have hA1 : allType1 d lv (done ++ [(.tree e cls (.node y outs), true)]) := by
    intro o ho hol
    unfold afterVert at ho
    rw [List.filter_append, List.map_append] at ho
    rcases List.mem_append.1 ho with ho | ho
    · have hrk := H.rank (o, true) (mem_afterVert ho)
      obtain ⟨lv', k', hk', -⟩ := ret_of_lowval_lt (o := o) (d := d) (by rw [hol]; exact hr)
      obtain ⟨lv'', k'', hk'', -⟩ :=
        ret_of_lowval_lt (o := DfsOut.tree e cls (.node y outs)) (d := d) hr
      simp only [DfsOut.cls] at hk''
      rw [hk'] at hol hrk ⊢
      rw [hk''] at hrk hlv ht1
      simp only [OutClass.lowval] at hol hlv
      simp only [OutClass.rank] at hrk
      cases k' <;> cases k'' <;> simp only [RetKind.rank, OutClass.isType1] at hrk ht1 ⊢ <;> omega
    · simp at ho
      rw [ho]; exact ht1
  have hcV : vertItem v ∉ x.spans.1 ++ x.spans.2 := by
    obtain ⟨top, htop, hCT⟩ := hC.top
    obtain ⟨above, below, hsp, -, -, -, -, -, -⟩ := hCT.split
    rw [if_pos rfl] at hsp
    obtain ⟨vt, htv, -, -, hvt3, hvt4, -, -⟩ := hsp
    have h := htop
    rw [hts₂, htv] at h
    rcases above with _ | ⟨x', r⟩
    · simp only [List.nil_append, List.cons_append, List.cons.injEq] at h
      rw [← h.1] at hvt3
      exact absurd (hvt3 hxv).1 (by rw [hxd]; exact Nat.ne_of_lt hr)
    · simp only [List.cons_append, List.cons.injEq] at h
      rw [← h.1] at hvt4
      exact hvt4 x (List.mem_cons_self ..)
  have hv₂ : v < sX.g.nv := by rw [HR.g]; exact H.hv
  have hsz₂ : 1 + sX.g.nv + sX.g.ne ≤ sX.items.size := by
    rw [HR.g]; exact Nat.le_trans H.hsz HR.size
  have hrest_e₂ : ∀ o' ∈ rest, ∀ e', subEdges o' e' → e' < sX.g.ne :=
    fun o' ho' e' h => by rw [HR.g]; exact (H.rest_nd o' ho' e' h).2
  have hrest_v₂ : ∀ o' ∈ rest, ∀ e₁ cls₁ child, o' = .tree e₁ cls₁ child →
      ∀ w ∈ child.verts, w < sX.g.nv :=
    fun o' ho' e₁ cls₁ child h w hw => by
      rw [HR.g]; exact (H.rest_nv o' ho' e₁ cls₁ child h w hw).2
  have ham : a ∈ sX.tstack := by rw [hts₂]; exact List.mem_cons_of_mem _ (List.mem_cons_self ..)
  unfold finishRest finishP finishTail condP maybeUnwrapNxt finishTstackTop
  simp only [wp_bind, wp_tstackSize, wp_nxt, wp_pure, wp_ite, hc', ↓reduceIte, Bool.not_true,
    Bool.false_eq_true, wp_get, wp_allocItem, wp_stackDir, wp_getItem, wp_modifyNxt,
    wp_mergeTstackTops, wp_cur, wp_makeVs, wp_modifyItem, wp_modifyCur]
  rw [hts₂]
  simp only [mergeTop, List.tail_cons, List.head!_cons, hnd, hnv, hxd, ← hlv, Nat.min_self]
  split_ifs with h1 h2 <;>
    (obtain ⟨i, hai, hiroot⟩ := mergeP_single hC hts₂ hcV hnv hnd hr hA1
     have hilt : i < sX.items.size :=
       hC.span_lt a ham i (by rw [hai]; unfold setSides; split <;> simp)
     have hhead : (getSide a.spans sX.stackDir[lv]!).head! = i := by
       rw [hai]; unfold setSides getSide; cases sX.stackDir[lv]! <;> rfl
     have hchi : sX.items[i]!.ch = Items.ch sX.items i := by
       rw [Items.ch_eq_getElem hilt, getElem!_pos sX.items i hilt])
  on_goal 2 =>
    rw [hhead, hchi]
    set S := getSide (x.spans.1 ++ (setSides sX.stackDir[lv]! (Items.ch sX.items i) []).1,
      (setSides sX.stackDir[lv]! (Items.ch sX.items i) []).2 ++ x.spans.2) sX.stackDir[lv]! with hSdef
    have hSmem : ∀ r, r ∈ S ↔ r = q ∨ r ∈ Items.ch sX.items i := by
      intro r; rw [hSdef, hcsp]; unfold setSides getSide
      cases sX.stackDir[lv]! <;> simp <;> tauto
    refine (fun hS hP => ⟨earCtx_mergeP hC hv₂ hr hsz₂ hts₂ hxd hcs hxv hnv hnd hai hA1
      hrest_e₂ hrest_v₂ (p := i) (S := S) (.inl rfl) hS ((hSmem q).2 (.inl rfl))
      (fun h => absurd rfl h) hP (fun j hj => Items.ch_modify_of_ne _ _ hj)
      (by simp) (by rw [Array.size_modify]; exact hilt) rfl rfl rfl rfl rfl,
      hF0.trans (outFrame_mergeP hts₂ hcs (by rw [hai, mem_setSides]; exact List.mem_singleton_self _)
        (.inl rfl) hS hP (fun j hj => Items.ch_modify_of_ne _ _ hj) rfl (by simp) rfl),
      fun hS0 => ctxShape_mergeP H.ctx hS0 H.hv hr H.hsz hst hnv hnd (p := i) (hts' := rfl)
        (hsd := fun k hk => HR.sd k hk)
        (Or.inl (by rw [hai, mem_setSides]; exact List.mem_singleton_self _)) rfl rfl rfl⟩) ?_ ?_
    · intro r hrS
      rcases (hSmem r).1 hrS with h | h
      · exact .inl h
      · exact .inr (.inr h)
    · intro p' c'
      unfold Items.IsParent
      by_cases hp' : p' = i
      · subst hp'
        rw [Items.ch_modify_at _ _ hilt]
        simp only [hSmem, true_and]
        tauto
      · rw [Items.ch_modify_of_ne _ _ hp']
        simp [hp']
  all_goals
    set S := getSide (x.spans.1 ++ a.spans.1, a.spans.2 ++ x.spans.2) sX.stackDir[lv]! with hSdef
    have hSmem : ∀ r, r ∈ S ↔ r = q ∨ r = i := by
      intro r; rw [hSdef, hcsp, hai]; unfold setSides getSide
      cases sX.stackDir[lv]! <;> simp <;> tauto
    refine (fun hS hP => ⟨earCtx_mergeP hC hv₂ hr hsz₂ hts₂ hxd hcs hxv hnv hnd hai hA1
      hrest_e₂ hrest_v₂ (p := sX.items.size) (S := S) (.inr (Nat.le_refl _)) hS
      ((hSmem q).2 (.inl rfl)) (fun _ => (hSmem i).2 (.inr rfl)) hP
      (fun j hj => by rw [Items.ch_modify_of_ne _ _ hj, Items.ch_push, if_neg hj])
      (by simp) (by simp) rfl rfl rfl rfl rfl,
      hF0.trans (outFrame_mergeP hts₂ hcs (by rw [hai, mem_setSides]; exact List.mem_singleton_self _)
        (.inr (Nat.le_refl _)) hS hP
        (fun j hj => by rw [Items.ch_modify_of_ne _ _ hj, Items.ch_push, if_neg hj])
        rfl (by simp) rfl),
      fun hS0 => ctxShape_mergeP H.ctx hS0 H.hv hr H.hsz hst hnv hnd (p := sX.items.size) (hts' := rfl)
        (hsd := fun k hk => HR.sd k hk) (Or.inr HR.size) rfl rfl rfl⟩) ?_ ?_
    · intro r hrS
      rcases (hSmem r).1 hrS with h | h
      · exact .inl h
      · exact .inr (.inl ⟨h, Nat.ne_of_lt hilt⟩)
    · intro p' c'
      unfold Items.IsParent
      by_cases hp' : p' = sX.items.size
      · subst hp'
        rw [Items.ch_modify_at _ _ (by simp), ch_of_ge (Nat.le_refl _), Array.getElem_push_eq]
        simp
      · rw [Items.ch_modify_of_ne _ _ hp', Items.ch_push, if_neg hp']
        simp [hp']

/-- The tree-edge step for a returning child (`lowval < d`): `tree_ret_shape` after loops 1–3, then
the P-check (`ctx_step_tree_ret_P`) and the first-edge vertex push, `earCtx_ret` over the final
top. -/
theorem ctx_step_tree_ret {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s : WalkState} {e : Nat} {cls : OutClass} {y : Nat} {outs : List DfsOut} {L : List TEntry}
    {push : Bool} {done' : List (DfsOut × Bool)} {hv' : Bool} {bE' : List (Nat → Prop)}
    {sv' : List Nat} {sd' : List Bool} {sE : WalkState} {dir' : Bool} {L' : List TEntry}
    {push' : Bool} {D₃ : Array Bool}
    (H : TreeSite v d done rest hasVert base bE sv sd s e cls y outs L push done' hv' bE' sv' sd'
      sE dir' L' push' D₃)
    (hr : cls.lowval d < d) :
    wp (finishEdge v d (.tree e cls (.node y outs))
      (L ++ s.tstack).length (hasVert || push))
      (fun hv'' s' =>
        EarCtx v d
          (done ++ [(.tree e cls (.node y outs), hasVert || push)])
          rest hv'' base bE sv sd s' ∧ OutFrame v (.tree e cls (.node y outs)) base s s' ∧
          (CtxShape v d done hasVert base s →
            CtxShape v d (done ++ [(.tree e cls (.node y outs), hasVert || push)]) hv'' base s'))
      (pushEnd sE D₃ L') := by
  obtain ⟨R, HR⟩ := tree_ret_shape H hr
  obtain ⟨top, htop, -⟩ := H.ctx.top
  rw [wp_finishEdge_tree v d (.tree e cls (.node y outs)) _ _ _ _ H.tree hr]
  generalize hb : (hasVert || push) = hv₀ at HR ⊢
  cases hv₀
  · simp only [Bool.false_eq_true, ↓reduceIte] at HR ⊢
    have hhp := Bool.or_eq_false_iff.1 hb
    have ht1 : cls.isType1 = false := by
      cases h : cls.isType1
      · rfl
      · exact absurd (H.hpush.2 ⟨hhp.1, hr, h⟩) (by rw [hhp.2]; decide)
    set sX := feS₂ d (.tree e cls (.node y outs)) (pushEnd sE D₃ L') with hsX
    unfold finishRest finishP finishTail condP
    simp only [wp_bind, wp_tstackSize, wp_nxt, wp_pure, wp_ite, DfsOut.cls, ht1, Bool.false_and,
      Bool.false_eq_true, ↓reduceIte, Bool.not_false, wp_pushVertTstack, wp_mergeTstackTops]
    split
    · obtain ⟨⟨c, R', hR, -⟩, -⟩ := HR.noVert rfl
      refine (fun hts => ⟨earCtx_ret H hr hb HR
        (TEntry.mergeInto ⟨v, d, sX.nxtEdgeIdx, setSides sX.stackDir[d]! [vertItem v] []⟩ c :: R')
        (Or.inr ⟨rfl, sX.nxtEdgeIdx, sX.stackDir[d]!, Or.inr ⟨c, R', hR, rfl⟩⟩) hts rfl rfl rfl
        (fun _ _ => rfl) rfl, outFrame_retTop HR htop H.L_spans ?_ hts rfl rfl,
        fun hS => ctxShape_ret H hS hr hb HR _
          (Or.inr ⟨rfl, sX.nxtEdgeIdx, by rw [← HR.sdd rfl]; exact Or.inr ⟨c, R', hR, rfl⟩⟩)
          (fun hhv => absurd hhv (by rw [hhp.1]; decide)) hts (fun _ _ => rfl)⟩) ?_
      · simp only [HR.tstack, hR, List.cons_append]
        rfl
      · intro t ht i hi
        rcases List.mem_cons.1 ht with rfl | ht
        · rcases (mem_mergeInto_spans _ _ _).1 hi with h | h
          · exact Or.inl (by (try dsimp only at h); rw [mem_setSides] at h; simpa using h)
          · exact Or.inr ⟨c, by rw [hR]; exact List.mem_cons_self .., h⟩
        · exact Or.inr ⟨t, by rw [hR]; exact List.mem_cons_of_mem _ ht, hi⟩
    · refine (fun hts => ⟨earCtx_ret H hr hb HR
        (⟨v, d, sX.nxtEdgeIdx, setSides sX.stackDir[d]! [vertItem v] []⟩ :: R)
        (Or.inr ⟨rfl, sX.nxtEdgeIdx, sX.stackDir[d]!, Or.inl rfl⟩) hts rfl rfl rfl
        (fun _ _ => rfl) rfl, outFrame_retTop HR htop H.L_spans ?_ hts rfl rfl,
        fun hS => ctxShape_ret H hS hr hb HR _
          (Or.inr ⟨rfl, sX.nxtEdgeIdx, by rw [← HR.sdd rfl]; exact Or.inl rfl⟩)
          (fun hhv => absurd hhv (by rw [hhp.1]; decide)) hts (fun _ _ => rfl)⟩) ?_
      · simp only [HR.tstack, List.cons_append]
      · intro t ht i hi
        rcases List.mem_cons.1 ht with rfl | ht
        · exact Or.inl (by rw [mem_setSides] at hi; simpa using hi)
        · exact Or.inr ⟨t, ht, hi⟩
  · simp only [↓reduceIte] at HR ⊢
    set sX := feS₃ v d (.tree e cls (.node y outs)) (L ++ s.tstack).length (pushEnd sE D₃ L')
      with hsX
    by_cases hc : result (condP v (cls.lowval d) cls.isType1) sX = true
    · exact ctx_step_tree_ret_P H hr HR _ hc
    · have hc' : (cls.isType1 && decide (sX.tstack.length ≥ 2) && (sX.tstack.tail.head!.vStart == v) &&
          (sX.tstack.tail.head!.topDepth == cls.lowval d)) = false := Bool.eq_false_iff.2 hc
      unfold finishRest finishP finishTail condP
      simp only [wp_bind, wp_tstackSize, wp_nxt, wp_pure, wp_ite, DfsOut.cls, hc',
        Bool.false_eq_true, ↓reduceIte, Bool.not_true]
      have hnoP : hasVert = true → cls.isType1 = true → ∀ a tl, s.tstack = a :: tl → a.vStart = v →
          a.topDepth ≠ cls.lowval d := by
        intro hhv ht1 a tl hst hav hd
        have hp : push = false := by
          cases hp : push
          · rfl
          · exact absurd (H.hpush.1 hp).1 (by rw [hhv]; decide)
        have hL0 : L = [] := by rw [H.hL, hp]; rfl
        obtain ⟨x, hRx, -⟩ := HR.vert rfl
        rw [HR.tstack, hRx, hL0, hst] at hc'
        simp [hav, hd, ht1] at hc'
      exact ⟨earCtx_ret H hr hb HR R (Or.inl ⟨rfl, rfl⟩) HR.tstack rfl rfl rfl (fun _ _ => rfl) rfl,
        outFrame_retTop HR htop H.L_spans (fun t ht i hi => Or.inr ⟨t, ht, hi⟩) HR.tstack rfl rfl,
        fun hS => ctxShape_ret H hS hr hb HR R (Or.inl ⟨rfl, rfl⟩) hnoP HR.tstack (fun _ _ => rfl)⟩

/-- `finishEdge` at a tree-edge site re-establishes the context with the out appended to `done`,
from the `TreeSite` data of `earAt_tree_of_ctx`: a case split on the class (bridge / component /
returning) over the three named admissions above. -/
theorem ctx_step_tree_frame {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s : WalkState} {e : Nat} {cls : OutClass} {y : Nat} {outs : List DfsOut} {L : List TEntry}
    {push : Bool} {done' : List (DfsOut × Bool)} {hv' : Bool} {bE' : List (Nat → Prop)}
    {sv' : List Nat} {sd' : List Bool} {sE : WalkState} {dir' : Bool} {L' : List TEntry}
    {push' : Bool} {D₃ : Array Bool}
    (H : TreeSite v d done rest hasVert base bE sv sd s e cls y outs L push done' hv' bE' sv' sd'
      sE dir' L' push' D₃) :
    wp (finishEdge v d (.tree e cls (.node y outs)) (L ++ s.tstack).length (hasVert || push))
      (fun hv'' s' => EarCtx v d (done ++ [(.tree e cls (.node y outs), hasVert || push)]) rest hv''
        base bE sv sd s' ∧ OutFrame v (.tree e cls (.node y outs)) base s s' ∧
        (CtxShape v d done hasVert base s →
          CtxShape v d (done ++ [(.tree e cls (.node y outs), hasVert || push)]) hv'' base s'))
      (pushEnd sE D₃ L') := by
  have key : cls = .bridge ∨ cls = .component ∨ cls.lowval d < d := by
    have ht := H.tree
    have hr := H.cls_ret
    cases cls with
    | bridge => exact .inl rfl
    | component => exact .inr (.inl rfl)
    | selfLoop => exact absurd ht (by decide)
    | ret lv k => exact .inr (.inr (hr lv k rfl))
  rcases key with hb | hc | hr
  · exact ctx_step_tree_bridge H hb
  · exact ctx_step_tree_comp H hc
  · exact ctx_step_tree_ret H hr

/-- `ctx_step_tree_frame` without the item frame. -/
theorem ctx_step_tree {v d : Nat} {done : List (DfsOut × Bool)} {rest : List DfsOut}
    {hasVert : Bool} {base : List TEntry} {bE : List (Nat → Prop)} {sv : List Nat} {sd : List Bool}
    {s : WalkState} {e : Nat} {cls : OutClass} {y : Nat} {outs : List DfsOut} {L : List TEntry}
    {push : Bool} {done' : List (DfsOut × Bool)} {hv' : Bool} {bE' : List (Nat → Prop)}
    {sv' : List Nat} {sd' : List Bool} {sE : WalkState} {dir' : Bool} {L' : List TEntry}
    {push' : Bool} {D₃ : Array Bool}
    (H : TreeSite v d done rest hasVert base bE sv sd s e cls y outs L push done' hv' bE' sv' sd'
      sE dir' L' push' D₃) :
    wp (finishEdge v d (.tree e cls (.node y outs)) (L ++ s.tstack).length (hasVert || push))
      (fun hv'' s' => EarCtx v d (done ++ [(.tree e cls (.node y outs), hasVert || push)]) rest hv''
        base bE sv sd s')
      (pushEnd sE D₃ L') :=
  wp_mono _ (ctx_step_tree_frame H) fun _ _ h => h.1

theorem wp_walkOutPre {v d : Nat} {o : DfsOut} {hasVert : Bool} {s : WalkState}
    {Q : Bool → WalkState → Prop}
    (h : ∀ (push : Bool) (L : List TEntry),
      (push = true ↔ hasVert = false ∧ o.cls.lowval d < d ∧ o.cls.isType1 = true) →
      L = (if push then [⟨v, d, s.nxtEdgeIdx, setSides (s.stackDir.set! d (if o.cls.lowval d ≥ d then
        false else !s.stackDir[o.cls.lowval d]!))[d]! [vertItem v] []⟩] else []) →
      Q (hasVert || push) { s with stackDir := s.stackDir.set! d (if o.cls.lowval d ≥ d then false
        else !s.stackDir[o.cls.lowval d]!), tstack := L ++ s.tstack }) :
    wp (walkOutPre v d o hasVert) Q s := by
  unfold walkOutPre
  simp only [wp_bind, wp_stackDir, wp_setStackDir, wp_ite, wp_pushVertTstack, wp_pure]
  split
  · next hc =>
    simp only [Bool.and_eq_true, Bool.not_eq_true', decide_eq_true_eq] at hc
    obtain ⟨⟨h0, h1⟩, h2⟩ := hc
    subst h0
    exact h true [⟨v, d, s.nxtEdgeIdx, setSides (s.stackDir.set! d (if o.cls.lowval d ≥ d then false
      else !s.stackDir[o.cls.lowval d]!))[d]! [vertItem v] []⟩] (by simp [h1, h2]) (by simp)
  · next hc =>
    have := h false [] (by
      simp only [Bool.false_eq_true, false_iff]
      rintro ⟨h0, h1, h2⟩
      exact hc (by simp [h0, h1, h2])) (by simp)
    rwa [Bool.or_false, List.nil_append] at this

/-- End-of-walk data of `walkTree (.node y outs) d` used by the parent's tree site: the child's
end-of-outs context, the push of the child's `V` entry (`pushEnd`), the stack directions below. -/
def TreeEnd (y d : Nat) (outs : List DfsOut) (base : List TEntry) (bE : List (Nat → Prop))
    (sv : List Nat) (sd : List Bool) (s' : WalkState) : Prop :=
  ∃ (done' : List (DfsOut × Bool)) (hv' : Bool) (sE : WalkState) (dir' : Bool) (L' : List TEntry)
    (push' : Bool) (D₃ : Array Bool),
    s' = pushEnd sE D₃ L' ∧ EarCtx y d done' [] hv' base bE sv sd sE ∧ done'.map (·.1) = outs ∧
    (push' = true ↔ hv' = false) ∧
    L' = (if push' then [⟨y, d, sE.nxtEdgeIdx, setSides dir' [vertItem y] []⟩] else []) ∧
    (∀ k, k < d → D₃[k]! = sE.stackDir[k]!) ∧
    (push' = true → d < D₃.size → dir' = true)

/-- `TreeEnd` with the end-of-outs shape `CtxShape` of the finished vertex. -/
def TreeEndS (y d : Nat) (outs : List DfsOut) (base : List TEntry) (bE : List (Nat → Prop))
    (sv : List Nat) (sd : List Bool) (s' : WalkState) : Prop :=
  ∃ (done' : List (DfsOut × Bool)) (hv' : Bool) (sE : WalkState) (dir' : Bool) (L' : List TEntry)
    (push' : Bool) (D₃ : Array Bool),
    s' = pushEnd sE D₃ L' ∧ EarCtx y d done' [] hv' base bE sv sd sE ∧
    CtxShape y d done' hv' base sE ∧ done'.map (·.1) = outs ∧
    (push' = true ↔ hv' = false) ∧
    L' = (if push' then [⟨y, d, sE.nxtEdgeIdx, setSides dir' [vertItem y] []⟩] else []) ∧
    (∀ k, k < d → D₃[k]! = sE.stackDir[k]!) ∧
    (push' = true → d < D₃.size → dir' = true)

/-- The per-out part of `EarOut` after `walkOutPre` (literally its inner match). -/
def OutMid (v d : Nat) (o : DfsOut) (hasVert' : Bool) (s₁ : WalkState) : Prop :=
  match o with
  | .tree _ _ child =>
    wp (modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }) (fun _ s₂ =>
      EarTree child (d + 1) s₂ ∧
      wp (walkTree child (d + 1)) (fun _ s₃ => FinishEar v d o s₁.tstack.length hasVert' s₃) s₂) s₁
  | .back .. => FinishEar v d o s₁.tstack.length hasVert' s₁

abbrev CTree (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  ∀ (g : Graph) (anc : List Nat) (base : List TEntry) (bE : List (Nat → Prop)) (sv : List Nat)
    (sd : List Bool) (pe : Nat → Prop), Types g s → d = anc.length → t.WF anc → t.Ends g →
    (anc ++ t.verts).Nodup →
    (∀ a ∈ anc, a < g.nv) → (∀ v ∈ t.verts, v < g.nv) → (∀ e ∈ t.edges, e < g.ne) → t.edges.Nodup →
    (∀ e', e' < g.ne → ∀ x, g.Inc e' x → x ∈ t.verts → e' ∈ t.edges ∨ pe e') →
    (∀ e', pe e' → ∀ x, g.Inc e' x → x ∈ anc ++ [t.v]) →
    s.stackVerts.size = g.nv → s.stackDir.size = g.nv →
    (∀ k, k < d → anc[k]? = some s.stackVerts[k]!) →
    EarCtx t.v d [] t.outs false base bE sv sd { s with stackVerts := s.stackVerts.set! d t.v } →
    EarTree t d s ∧ wp (walkTree t d) (fun _ s' => TreeEnd t.v d t.outs base bE sv sd s') s ∧
    wp (walkTree t d) (fun _ s' => KeptRel g s t.verts t.edges s' ∧
      ∀ top, s'.tstack = top ++ base → ∀ t' ∈ top, ∀ i ∈ t'.spans.1 ++ t'.spans.2,
        (∃ x ∈ t.verts, i = vertItem x) ∨ (∃ e ∈ t.edges, i = edgeItem g e) ∨
        s.items.size ≤ i) s ∧
    wp (walkTree t d) (fun _ s' => TreeEndS t.v d t.outs base bE sv sd s') s

abbrev COuts (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  ∀ (g : Graph) (anc : List Nat) (outs₀ : List DfsOut) (done : List (DfsOut × Bool))
    (base : List TEntry) (bE : List (Nat → Prop)) (sv : List Nat) (sd : List Bool) (pe : Nat → Prop),
    Types g s → d = anc.length → done.map (·.1) ++ outs = outs₀ →
    outs₀.Pairwise (fun a b => a.cls.rank ≤ b.cls.rank) → (∀ o ∈ outs₀, o.WF anc v) →
    (∀ o ∈ outs₀, DfsOut.Ends g v o) → (anc ++ v :: DfsOut.vertsList outs₀).Nodup →
    (∀ a ∈ anc, a < g.nv) → v < g.nv → (∀ w ∈ DfsOut.vertsList outs₀, w < g.nv) →
    (∀ e ∈ DfsOut.edgesList outs₀, e < g.ne) → (DfsOut.edgesList outs₀).Nodup →
    (∀ e', e' < g.ne → ∀ x, g.Inc e' x → x ∈ v :: DfsOut.vertsList outs₀ → e' ∈ DfsOut.edgesList outs₀ ∨ pe e') →
    (∀ e', pe e' → ∀ x, g.Inc e' x → x ∈ anc ++ [v]) →
    s.stackVerts.size = g.nv → s.stackDir.size = g.nv →
    (∀ k, k < d → anc[k]? = some s.stackVerts[k]!) →
    EarCtx v d done outs hasVert base bE sv sd s →
    EarOuts v d outs hasVert s ∧ wp (walkOuts v d outs hasVert)
      (fun hv' s' => ∃ done', done'.map (·.1) = outs₀ ∧ EarCtx v d done' [] hv' base bE sv sd s') s ∧
    (∀ (s₀ : WalkState) (sz₀ : Nat), s₀.items.size ≤ sz₀ → sz₀ ≤ s.items.size →
      KeptRel g s₀ (v :: DfsOut.vertsList outs₀) (DfsOut.edgesList outs₀) s →
      Fresh g v outs₀ sz₀ base s →
      wp (walkOuts v d outs hasVert) (fun _ s' => sz₀ ≤ s'.items.size ∧ s'.g = g ∧
        KeptRel g s₀ (v :: DfsOut.vertsList outs₀) (DfsOut.edgesList outs₀) s' ∧
        Fresh g v outs₀ sz₀ base s') s) ∧
    (CtxShape v d done hasVert base s → wp (walkOuts v d outs hasVert)
      (fun hv' s' => ∃ done', done'.map (·.1) = outs₀ ∧ EarCtx v d done' [] hv' base bE sv sd s' ∧
        CtxShape v d done' hv' base s') s)

abbrev COut (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  ∀ (g : Graph) (anc : List Nat) (outs₀ rest : List DfsOut) (done : List (DfsOut × Bool))
    (base : List TEntry) (bE : List (Nat → Prop)) (sv : List Nat) (sd : List Bool) (pe : Nat → Prop),
    Types g s → d = anc.length → done.map (·.1) ++ o :: rest = outs₀ →
    outs₀.Pairwise (fun a b => a.cls.rank ≤ b.cls.rank) → (∀ o ∈ outs₀, o.WF anc v) →
    (∀ o ∈ outs₀, DfsOut.Ends g v o) → (anc ++ v :: DfsOut.vertsList outs₀).Nodup →
    (∀ a ∈ anc, a < g.nv) → v < g.nv → (∀ w ∈ DfsOut.vertsList outs₀, w < g.nv) →
    (∀ e ∈ DfsOut.edgesList outs₀, e < g.ne) → (DfsOut.edgesList outs₀).Nodup →
    (∀ e', e' < g.ne → ∀ x, g.Inc e' x → x ∈ v :: DfsOut.vertsList outs₀ → e' ∈ DfsOut.edgesList outs₀ ∨ pe e') →
    (∀ e', pe e' → ∀ x, g.Inc e' x → x ∈ anc ++ [v]) →
    s.stackVerts.size = g.nv → s.stackDir.size = g.nv →
    (∀ k, k < d → anc[k]? = some s.stackVerts[k]!) →
    EarCtx v d done (o :: rest) hasVert base bE sv sd s →
    EarOut v d o hasVert s ∧ wp (walkOut v d o hasVert)
      (fun hv' s' => ∃ hvF, EarCtx v d (done ++ [(o, hvF)]) rest hv' base bE sv sd s') s ∧
    wp (walkOut v d o hasVert) (fun _ s' => OutFrame v o base s s') s ∧
    (CtxShape v d done hasVert base s → wp (walkOut v d o hasVert)
      (fun hv' s' => ∃ hvF, EarCtx v d (done ++ [(o, hvF)]) rest hv' base bE sv sd s' ∧
        CtxShape v d (done ++ [(o, hvF)]) hv' base s') s)

mutual
theorem cTree : ∀ (t : DfsTree) (d : Nat) (s : WalkState), CTree t d s
  | .node v outs, d, s => by
    intro g anc base bE sv sd pe hT hd hwf hends hnd hanc hvlt helt hen hcomp hpe hsz hsdz hancsv hC
    simp only [DfsTree.verts, DfsTree.edges, DfsTree.v, DfsTree.outs] at hnd hvlt helt hen hcomp hpe hC ⊢
    simp only [DfsTree.WF] at hwf
    simp only [DfsTree.Ends] at hends
    have hv : v < g.nv := hvlt v (List.mem_cons_self ..)
    obtain ⟨hEO, hpost, hfr, hpostS⟩ := cOuts v d outs false { s with stackVerts := s.stackVerts.set! d v } g anc
      outs [] base bE sv sd pe ⟨hT.g_eq, hT.size, hT.root, hT.vert, hT.edge⟩ hd rfl hwf.1 hwf.2 hends hnd
      hanc hv (fun w hw => hvlt w (List.mem_cons_of_mem _ hw)) helt hen hcomp hpe (by simp [hsz]) hsdz
      (fun k hk => by
        show anc[k]? = some (s.stackVerts.set! d v)[k]!
        rw [getElem!_set!_ne' _ _ _ _ (Nat.ne_of_lt hk)]; exact hancsv k hk) hC
    refine ⟨by unfold EarTree; exact hEO, ?_, ?_, ?_⟩
    · unfold walkTree
      simp only [wp_bind, wp_modify]
      refine wp_mono _ hpost fun hv' s' ⟨done', hmap, hC'⟩ => ?_
      cases hv'
      · simp only [Bool.false_eq_true, ↓reduceIte, wp_bind, wp_setStackDir, wp_pushVertTstack]
        exact ⟨done', false, s', (s'.stackDir.set! d true)[d]!, _, true, s'.stackDir.set! d true, rfl,
          hC', hmap, by simp, rfl, fun k hk => getElem!_set!_ne' s'.stackDir d k true (Nat.ne_of_lt hk),
          fun _ hd => getElem!_set!_self' _ _ _ (by simpa using hd)⟩
      · exact ⟨done', true, s', false, [], false, s'.stackDir, rfl, hC', hmap, by simp, rfl,
          fun _ _ => rfl, fun h => absurd h Bool.false_ne_true⟩
    · have hF0 : Fresh g v outs s.items.size base { s with stackVerts := s.stackVerts.set! d v } := by
        intro top htop t ht
        obtain ⟨top₀, htop₀, hCT⟩ := hC.top
        have h0 : top₀ = [] := by
          by_contra h
          obtain ⟨o, ho, -⟩ := hCT.ret h
          simp at ho
        rw [h0] at htop₀
        have : top = [] := List.append_cancel_right (htop.symm.trans htop₀)
        simp [this] at ht
      have hfr' := hfr s s.items.size (Nat.le_refl _) (Nat.le_refl _)
        (fun j _ _ _ => ⟨rfl, fun p => Iff.rfl⟩) hF0
      unfold walkTree
      simp only [wp_bind, wp_modify]
      refine wp_mono _ (wp_and hpost hfr') fun hv' s' ⟨⟨_, _, hC'⟩, hsz', hg', hK', hF'⟩ => ?_
      cases hv'
      · simp only [Bool.false_eq_true, ↓reduceIte, wp_bind, wp_setStackDir, wp_pushVertTstack]
        refine ⟨hK', ?_⟩
        intro top htop t' ht i hi
        cases top with
        | nil =>
          exfalso
          obtain ⟨top₁, htop₁, -⟩ := hC'.top
          have h := congrArg List.length htop
          simp only [List.length_cons, htop₁, List.length_append, List.nil_append] at h
          omega
        | cons t₀ top₁ =>
          simp only [List.cons_append, List.cons.injEq] at htop
          obtain ⟨rfl, htop₁⟩ := htop
          rcases List.mem_cons.1 ht with rfl | ht
          · exact Or.inl ⟨v, List.mem_cons_self .., by rw [mem_setSides] at hi; simpa using hi⟩
          · exact hF' top₁ htop₁ t' ht i hi
      · exact ⟨hK', hF'⟩
    · have hts0 : ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).tstack = base := by
        obtain ⟨top₀, htop₀, hCT⟩ := hC.top
        have h0 : top₀ = [] := by
          by_contra h
          obtain ⟨o, ho, -⟩ := hCT.ret h
          simp at ho
        rw [h0] at htop₀; exact htop₀
      unfold walkTree
      simp only [wp_bind, wp_modify]
      refine wp_mono _ (hpostS (ctxShape_init hts0)) fun hv' s' ⟨done', hmap, hC', hS'⟩ => ?_
      cases hv'
      · simp only [Bool.false_eq_true, ↓reduceIte, wp_bind, wp_setStackDir, wp_pushVertTstack]
        exact ⟨done', false, s', (s'.stackDir.set! d true)[d]!, _, true, s'.stackDir.set! d true, rfl,
          hC', hS', hmap, by simp, rfl,
          fun k hk => getElem!_set!_ne' s'.stackDir d k true (Nat.ne_of_lt hk),
          fun _ hd => getElem!_set!_self' _ _ _ (by simpa using hd)⟩
      · exact ⟨done', true, s', false, [], false, s'.stackDir, rfl, hC', hS', hmap, by simp, rfl,
          fun _ _ => rfl, fun h => absurd h Bool.false_ne_true⟩

theorem cOuts : ∀ (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState),
    COuts v d outs hasVert s
  | v, d, [], hasVert, s => by
    intro g anc outs₀ done base bE sv sd pe hT hd hmap hsort hwf hends hnd hanc hv hw he hen hcomp hpe hsz hsdz
      hancsv hC
    refine ⟨?_, ?_, ?_, ?_⟩
    · unfold EarOuts VertBook
      intro h
      exact ⟨by rw [hT.g_eq]; exact hv, hC.vert_book h⟩
    · unfold walkOuts; simp only [wp_pure]
      exact ⟨done, by simpa using hmap, hC⟩
    · intro s₀ sz₀ hs₀ hsz₀ hK hF
      unfold walkOuts; simp only [wp_pure]
      exact ⟨hsz₀, hT.g_eq, hK, hF⟩
    · intro hS
      unfold walkOuts; simp only [wp_pure]
      exact ⟨done, by simpa using hmap, hC, hS⟩
  | v, d, o :: rest, hasVert, s => by
    intro g anc outs₀ done base bE sv sd pe hT hd hmap hsort hwf hends hnd hanc hv hw he hen hcomp hpe hsz hsdz
      hancsv hC
    obtain ⟨hEo, hstep, hOF, hstepS⟩ := cOut v d o hasVert s g anc outs₀ rest done base bE sv sd pe hT hd hmap hsort
      hwf hends hnd hanc hv hw he hen hcomp hpe hsz hsdz hancsv hC
    have hK0 : wp (walkOut v d o hasVert) (fun _ s' => Keep (d + 1) 0 s s') s :=
      kOut v d o hasVert s g (d + 1) 0 s hT (Nat.le_refl _) (by have := hT.size; omega)
        (vertItem_ne_zero v) (fun w _ => vertItem_ne_zero w) (fun e _ => edgeItem_ne_zero g e) Keep.refl
    have hrest : wp (walkOut v d o hasVert) (fun hv' s' => EarOuts v d rest hv' s' ∧
        wp (walkOuts v d rest hv') (fun hv'' s'' => ∃ done', done'.map (·.1) = outs₀ ∧
          EarCtx v d done' [] hv'' base bE sv sd s'') s') s := by
      refine wp_mono _ (wp_and hstep hK0) fun hv' s' ⟨⟨hvF, hC'⟩, hK⟩ => ?_
      have h := cOuts v d rest hv' s' g anc outs₀ (done ++ [(o, hvF)]) base bE sv sd pe (hT.of_keep hK) hd
        (by simpa using hmap) hsort hwf hends hnd hanc hv hw he hen hcomp hpe (hK.sv.trans hsz)
        (hK.sd.trans hsdz)
        (fun k hk => by rw [hancsv k hk, hC.sv k hk.le, hC'.sv k hk.le]) hC'
      exact ⟨h.1, h.2.1⟩
    refine ⟨?_, ?_, ?_, ?_⟩
    · unfold EarOuts
      exact ⟨hEo, wp_mono _ hrest fun _ _ h => h.1⟩
    · unfold walkOuts; simp only [wp_bind]
      exact wp_mono _ hrest fun _ _ h => h.2
    · intro s₀ sz₀ hs₀ hsz₀ hK hF
      have ho : o ∈ outs₀ := by rw [← hmap]; exact List.mem_append_right _ (List.mem_cons_self ..)
      have hrest' : wp (walkOut v d o hasVert) (fun hv' s' => wp (walkOuts v d rest hv')
          (fun _ s'' => sz₀ ≤ s''.items.size ∧ s''.g = g ∧
            KeptRel g s₀ (v :: DfsOut.vertsList outs₀) (DfsOut.edgesList outs₀) s'' ∧
            Fresh g v outs₀ sz₀ base s'') s') s := by
        refine wp_mono _ (wp_and hstep (wp_and hK0 hOF)) fun hv' s' ⟨⟨hvF, hC'⟩, hK', hOF'⟩ => ?_
        obtain ⟨hsz₁, -, hK₁, hF₁⟩ := outFrame_step ho hT.g_eq hsz₀ hs₀ hK hF hOF'
        exact (cOuts v d rest hv' s' g anc outs₀ (done ++ [(o, hvF)]) base bE sv sd pe
          (hT.of_keep hK') hd (by simpa using hmap) hsort hwf hends hnd hanc hv hw he hen hcomp hpe
          (hK'.sv.trans hsz) (hK'.sd.trans hsdz)
          (fun k hk => by rw [hancsv k hk, hC.sv k hk.le, hC'.sv k hk.le]) hC').2.2.1
          s₀ sz₀ hs₀ hsz₁ hK₁ hF₁
      unfold walkOuts; simp only [wp_bind]
      exact hrest'
    · intro hS
      unfold walkOuts; simp only [wp_bind]
      refine wp_mono _ (wp_and (hstepS hS) hK0) fun hv' s' ⟨⟨hvF, hC', hS'⟩, hK⟩ => ?_
      exact (cOuts v d rest hv' s' g anc outs₀ (done ++ [(o, hvF)]) base bE sv sd pe (hT.of_keep hK) hd
        (by simpa using hmap) hsort hwf hends hnd hanc hv hw he hen hcomp hpe (hK.sv.trans hsz)
        (hK.sd.trans hsdz)
        (fun k hk => by rw [hancsv k hk, hC.sv k hk.le, hC'.sv k hk.le]) hC').2.2.2 hS'

theorem cOut : ∀ (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState), COut v d o hasVert s
  | v, d, o, hasVert, s => by
    intro g anc outs₀ rest done base bE sv sd pe hT hd hmap hsort hwf hends hnd hanc hv hw he hen hcomp hpe hsz hsdz
      hancsv hC
    subst hd hmap
    have hgs : s.g = g := hT.g_eq
    have hv' : v < s.g.nv := by rw [hgs]; exact hv
    have hndv : (anc ++ [v]).Nodup := List.Nodup.sublist
      (List.Sublist.append_left (List.cons_sublist_cons.2 (List.nil_sublist _)) _) hnd
    have hdlt : anc.length < g.nv := by
      have := length_le_nv hndv fun w hw => by
        rcases List.mem_append.1 hw with hw | hw
        · exact hanc w hw
        · exact (List.mem_singleton.1 hw) ▸ hv
      simp at this; omega
    have hsd : anc.length < s.stackDir.size := by rw [hsdz]; exact hdlt
    have ho : o ∈ done.map (·.1) ++ o :: rest := by simp
    have hrank : ∀ o' ∈ done, o'.1.cls.rank ≤ o.cls.rank := fun o' ho' =>
      (List.pairwise_append.1 hsort).2.2 _ (List.mem_map_of_mem ho') o (List.mem_cons_self ..)
    have hinc : ∀ o' ∈ done, o'.1.e < s.g.ne ∧ s.g.Inc o'.1.e v := fun o' ho' => by
      have hm : o'.1 ∈ done.map (·.1) ++ o :: rest :=
        List.mem_append_left _ (List.mem_map_of_mem ho')
      rw [hgs]
      exact ⟨he _ (mem_subEdges_edgesList.2 ⟨_, hm, subEdges_e _⟩), inc_e_of_ends (hends _ hm)⟩
    have hen' := hen
    rw [DfsOut.edgesList_eq, List.flatMap_append, List.flatMap_cons] at hen'
    have hdisjE := List.disjoint_of_nodup_append hen'
    have hnd_o : ∀ o' ∈ done, ∀ e', subEdges o'.1 e' → ¬ subEdges o e' := fun o' ho' e' hs hs' =>
      hdisjE (by rw [← DfsOut.edgesList_eq]
                 exact mem_subEdges_edgesList.2 ⟨_, List.mem_map_of_mem ho', hs⟩)
        (List.mem_append_left _ (mem_edges_of_subEdges hs'))
    have hnd_vo : (o.edges ++ rest.flatMap DfsOut.edges).Nodup := (List.nodup_append.1 hen').2.1
    have hndV := hnd
    rw [DfsOut.vertsList_eq, List.flatMap_append, List.flatMap_cons] at hndV
    have hwf_o := hwf o ho
    have hends_o := hends o ho
    unfold EarOut
    rw [walkOut_eq]
    simp only [wp_bind]
    have hvb0 : VertBook v hasVert s := fun h => ⟨hv', hC.vert_book h⟩
    have hcore : wp (walkOutPre v anc.length o hasVert) (fun hvF s₁ =>
        OutMid v anc.length o hvF s₁ ∧ wp (walkOutRest v anc.length o hvF)
          (fun hv' s' => (∃ hvF', EarCtx v anc.length (done ++ [(o, hvF')]) rest hv' base bE sv sd s') ∧
            OutFrame v o base s s' ∧
            (CtxShape v anc.length done hasVert base s → ∃ hvF',
              EarCtx v anc.length (done ++ [(o, hvF')]) rest hv' base bE sv sd s' ∧
              CtxShape v anc.length (done ++ [(o, hvF')]) hv' base s'))
          s₁) s := by
      refine wp_walkOutPre fun push L hpush hL => ?_
      cases o with
      | back e dest cls =>
        obtain ⟨hb, hcls⟩ := back_wf_facts hwf_o
        have hnd_e : ∀ o' ∈ done, ∀ e', subEdges o'.1 e' → e' ≠ e := fun o' ho' e' hs h =>
          hnd_o o' ho' e' hs (by rw [h]; exact subEdges_e _)
        have he_e : e < s.g.ne := by
          rw [hgs]; exact he _ (mem_subEdges_edgesList.2 ⟨_, ho, subEdges_e _⟩)
        have hszS : 1 + s.g.nv + s.g.ne ≤ s.items.size := by rw [hgs]; exact hT.size
        have hend' : Items.PairEq (v, dest) s.g.edges[e]! := by
          have h := hends_o
          rw [DfsOut.Ends] at h
          rw [hgs]
          exact h
        have hself : anc.length ≤ cls.lowval anc.length → dest = v := by
          intro hge
          have h := hwf_o
          rw [DfsOut.WF] at h
          obtain ⟨i, hi, hcls'⟩ := h
          have hi' : i < anc.length + 1 := by
            have := (List.getElem?_eq_some_iff.1 hi).1; simpa using this
          by_cases hid : i = anc.length
          · subst hid
            rw [List.getElem?_append_right (Nat.le_refl _), Nat.sub_self] at hi
            simp at hi
            omega
          · exfalso
            have hlt : i < anc.length := by omega
            rw [hcls'] at hge
            simp [classify, Nat.not_le.2 hlt, OutClass.lowval] at hge
        have hrest_e : ∀ o' ∈ rest, ∀ e', subEdges o' e' → e' ≠ e ∧ e' < s.g.ne := by
          intro o' ho' e' hs
          have hm : o' ∈ done.map (·.1) ++ .back e dest cls :: rest :=
            List.mem_append_right _ (List.mem_cons_of_mem _ ho')
          refine ⟨fun h => ?_, by rw [hgs]; exact he _ (mem_subEdges_edgesList.2 ⟨_, hm, hs⟩)⟩
          exact List.disjoint_of_nodup_append hnd_vo (by simp [DfsOut.edges])
            (List.mem_flatMap.2 ⟨o', ho', h ▸ mem_edges_of_subEdges hs⟩)
        have hrest_v : ∀ o' ∈ rest, ∀ e₁ cls₁ child, o' = .tree e₁ cls₁ child →
            ∀ y ∈ child.verts, y < s.g.nv := by
          intro o' ho' e₁ cls₁ child ho₁ y hy
          rw [hgs]
          exact hw _ (mem_vertsList.2 ⟨o', List.mem_append_right _ (List.mem_cons_of_mem _ ho'),
            e₁, cls₁, child, ho₁, hy⟩)
        have hdest : cls.lowval anc.length < anc.length →
            s.stackVerts[cls.lowval anc.length]! = dest := by
          intro hlt
          have h := hwf_o
          rw [DfsOut.WF] at h
          obtain ⟨i, hi, hcls'⟩ := h
          obtain ⟨lv, k, hret, -⟩ := ret_of_lowval_lt (o := DfsOut.back e dest cls) hlt
          simp only [DfsOut.cls] at hret
          rw [hret] at hcls' hlt ⊢
          have hil : i = lv := (classify_eq_ret_iff.1 hcls'.symm).1
          simp only [OutClass.lowval] at hlt ⊢
          subst hil
          rw [List.getElem?_append_left hlt, hancsv i hlt] at hi
          exact Option.some.inj hi
        have hnc : ∀ o' ∈ done, ∀ e', subEdges o'.1 e' → ∀ x, s.g.Inc e' x →
            ∀ o'' ∈ rest, ∀ e₁ cls₁ child, o'' = .tree e₁ cls₁ child → x ∉ child.verts := by
          intro o' ho' e' hs x hx o'' ho'' e₁ cls₁ child ho₁ hxC
          have hm : o'.1 ∈ done.map (·.1) ++ DfsOut.back e dest cls :: rest :=
            List.mem_append_left _ (List.mem_map_of_mem ho')
          have hxR : x ∈ rest.flatMap DfsOut.verts :=
            List.mem_flatMap.2 ⟨o'', ho'', by rw [ho₁]; exact hxC⟩
          rw [hgs] at hx
          rcases endsOut_wf g anc v o'.1 (hwf _ hm) (hends _ hm) e' hs x hx with h | rfl | h
          · exact List.disjoint_of_nodup_append hndV h (List.mem_cons_of_mem _
              (List.mem_append_right _ (List.mem_append_right _ hxR)))
          · exact (List.nodup_cons.1 (List.nodup_append.1 hndV).2.1).1
              (List.mem_append_right _ (List.mem_append_right _ hxR))
          · exact List.disjoint_of_nodup_append (List.nodup_cons.1 (List.nodup_append.1 hndV).2.1).2
              (List.mem_flatMap.2 ⟨o'.1, List.mem_map_of_mem ho', h⟩) (List.mem_append_right _ hxR)
        refine ⟨⟨fun h => hC.vert_book (Bool.or_eq_false_iff.1 h).1,
          earAt_back_of_ctx hC hv' hsd hrank hinc hnd_e hb hcls L push hpush hL⟩, ?_⟩
        unfold walkOutRest
        rw [wp_bind, wp_tstackSize]
        exact wp_mono _ (ctx_step_back_frame hC hv' hsd hrank hinc hnd_e hb hcls he_e hszS hend' hself
          hrest_e hrest_v hdest hnc L push hpush hL) fun _ _ h =>
            ⟨⟨_, h.1⟩, h.2.1, fun hS => ⟨_, h.1, h.2.2 hS⟩⟩
      | tree e cls child =>
        obtain ⟨y, outs'⟩ := child
        obtain ⟨ht, hcls⟩ := tree_wf_facts hwf_o
        have hwf_c := hwf_o
        rw [DfsOut.WF] at hwf_c
        have hends_c := hends_o
        rw [DfsOut.Ends] at hends_c
        have hEt := hends_c.2
        rw [DfsTree.Ends] at hEt
        simp only [DfsOut.verts] at hndV
        have hyv : y ∈ (DfsTree.node y outs').verts := List.mem_cons_self ..
        have hy : y < g.nv := hw y (mem_vertsList_of_verts ho hyv)
        have he_lt : e < g.ne := he e (mem_subEdges_edgesList.2 ⟨_, ho, subEdges_e _⟩)
        have hnd_vC : (v :: ((done.map (·.1)).flatMap DfsOut.verts ++
            ((DfsTree.node y outs').verts ++ rest.flatMap DfsOut.verts))).Nodup :=
          (List.nodup_append.1 hndV).2.1
        have hv_nc : v ∉ (DfsTree.node y outs').verts := fun h =>
          (List.nodup_cons.1 hnd_vC).1 (List.mem_append_right _ (List.mem_append_left _ h))
        have he_ne : e ∉ (DfsTree.node y outs').edges :=
          (List.nodup_cons.1 (List.nodup_append.1 hnd_vo).1).1
        have hen_c : (DfsTree.node y outs').edges.Nodup :=
          (List.nodup_cons.1 (List.nodup_append.1 hnd_vo).1).2
        have hnd_c : (anc ++ [v] ++ (DfsTree.node y outs').verts).Nodup := by
          rw [List.append_assoc, List.singleton_append]
          exact List.Nodup.sublist (List.Sublist.append_left (List.cons_sublist_cons.2
            ((List.sublist_append_left _ _).trans (List.sublist_append_right _ _))) _) hndV
        have hdlt2 : anc.length + 1 < g.nv := by
          have hnd2 : (anc ++ [v] ++ [y]).Nodup := List.Nodup.sublist
            (List.Sublist.append_left (List.cons_sublist_cons.2 (List.nil_sublist _)) _) hnd_c
          have := length_le_nv hnd2 fun w hw => by
            rcases List.mem_append.1 hw with hw | hw
            · rcases List.mem_append.1 hw with hw | hw
              · exact hanc w hw
              · exact (List.mem_singleton.1 hw) ▸ hv
            · exact (List.mem_singleton.1 hw) ▸ hy
          simp at this; omega
        have hdone_nc : ∀ o' ∈ done, ∀ e', subEdges o'.1 e' → ∀ x, s.g.Inc e' x →
            x ∉ (DfsTree.node y outs').verts := fun o' ho' e' hs x hx hxC => by
          have hm : o'.1 ∈ done.map (·.1) ++ DfsOut.tree e cls (.node y outs') :: rest :=
            List.mem_append_left _ (List.mem_map_of_mem ho')
          rw [hgs] at hx
          rcases endsOut_wf g anc v o'.1 (hwf _ hm) (hends _ hm) e' hs x hx with h | rfl | h
          · exact List.disjoint_of_nodup_append hndV h
              (List.mem_cons_of_mem _ (List.mem_append_right _ (List.mem_append_left _ hxC)))
          · exact hv_nc hxC
          · exact List.disjoint_of_nodup_append (List.nodup_cons.1 hnd_vC).2
              (List.mem_flatMap.2 ⟨o'.1, List.mem_map_of_mem ho', h⟩) (List.mem_append_left _ hxC)
        have hrest_nc : ∀ o' ∈ rest, ∀ e', subEdges o' e' → ∀ x, s.g.Inc e' x →
            x ∉ (DfsTree.node y outs').verts := fun o' ho' e' hs x hx hxC => by
          have hm : o' ∈ done.map (·.1) ++ DfsOut.tree e cls (.node y outs') :: rest :=
            List.mem_append_right _ (List.mem_cons_of_mem _ ho')
          rw [hgs] at hx
          rcases endsOut_wf g anc v o' (hwf _ hm) (hends _ hm) e' hs x hx with h | rfl | h
          · exact List.disjoint_of_nodup_append hndV h
              (List.mem_cons_of_mem _ (List.mem_append_right _ (List.mem_append_left _ hxC)))
          · exact hv_nc hxC
          · exact List.disjoint_of_nodup_append
              (List.nodup_append.1 (List.nodup_cons.1 hnd_vC).2).2.1 hxC
              (List.mem_flatMap.2 ⟨o', ho', h⟩)
        have hcomp_c : ∀ e', e' < g.ne → ∀ x, g.Inc e' x → x ∈ (DfsTree.node y outs').verts →
            e' ∈ (DfsTree.node y outs').edges ∨ e' = e := by
          intro e' he' x hx hxC
          have hx' : x ∈ v :: DfsOut.vertsList
              (done.map (·.1) ++ DfsOut.tree e cls (.node y outs') :: rest) :=
            List.mem_cons_of_mem _ (mem_vertsList_of_verts ho hxC)
          rcases hcomp e' he' x hx hx' with h | h
          · obtain ⟨o', ho', hs⟩ := mem_subEdges_edgesList.1 h
            rcases List.mem_append.1 ho' with ho' | ho'
            · obtain ⟨p, hp, rfl⟩ := List.mem_map.1 ho'
              exact absurd hxC (hdone_nc p hp e' hs x (by rw [hgs]; exact hx))
            · rcases List.mem_cons.1 ho' with rfl | ho'
              · rcases (show e' = e ∨ e' ∈ (DfsTree.node y outs').edges from hs) with h | h
                · exact Or.inr h
                · exact Or.inl h
              · exact absurd hxC (hrest_nc o' ho' e' hs x (by rw [hgs]; exact hx))
          · rcases List.mem_append.1 (hpe e' h x hx) with h | h
            · exact absurd (List.mem_cons_of_mem _ (List.mem_append_right _
                (List.mem_append_left _ hxC))) (List.disjoint_of_nodup_append hndV h)
            · have hxv := List.mem_singleton.1 h
              subst hxv
              exact absurd hxC hv_nc
        have hpe_c : ∀ e', e' = e → ∀ x, g.Inc e' x → x ∈ anc ++ [v] ++ [y] := by
          rintro e' rfl x hx
          rcases eq_of_inc_pairEq hends_c.1 hx with rfl | rfl
          · simp [DfsTree.v]
          · simp
        obtain ⟨sv', sd', hC₂⟩ := ctx_init_child hC hv' (by rw [hgs]; exact hy) hsd
          (by rw [hsz]; exact hdlt2) hinc
          (fun e' he' hb x hx hi => by
            obtain ⟨o', ho', -, hs⟩ := (hC.vert_edges e' he').1 hb
            exact hdone_nc o' ho' e' hs x hi hx)
          (List.nodup_append.1 (List.nodup_append.1 (List.nodup_cons.1 hnd_vC).2).2.1).1
          (by rw [hgs]; have := hT.size; omega) push L hpush hL
        obtain ⟨S₂, hS₂⟩ : ∃ S₂ : WalkState, S₂ = { s with
            stackDir := (s.stackDir.set! anc.length (if cls.lowval anc.length ≥ anc.length then false
              else !s.stackDir[cls.lowval anc.length]!)),
            tstack := L ++ s.tstack,
            firstOccurrence := s.firstOccurrence.set! anc.length s.g.ne } := ⟨_, rfl⟩
        have hT₂ : Types g S₂ := by subst hS₂; exact ⟨hT.g_eq, hT.size, hT.root, hT.vert, hT.edge⟩
        have hC₂' : EarCtx (DfsTree.node y outs').v (anc.length + 1) [] (DfsTree.node y outs').outs false
            (L ++ s.tstack) ((L ++ s.tstack).map fun (t : TEntry) (e' : Nat) => t.edges s.g s.items e')
            sv' sd' { S₂ with stackVerts := S₂.stackVerts.set! (anc.length + 1) (DfsTree.node y outs').v } := by
          subst hS₂; exact hC₂
        obtain ⟨hET, -, hfr, hendS⟩ := cTree (.node y outs') (anc.length + 1) S₂ g (anc ++ [v]) (L ++ s.tstack)
          ((L ++ s.tstack).map fun (t : TEntry) (e' : Nat) => t.edges s.g s.items e') sv' sd' (· = e) hT₂
          (by simp) hwf_c.1 hends_c.2 hnd_c (fun a ha => by
            rcases List.mem_append.1 ha with ha | ha
            · exact hanc a ha
            · exact (List.mem_singleton.1 ha) ▸ hv)
          (fun w hw' => hw w (mem_vertsList_of_verts ho hw'))
          (fun e' h => he e' (mem_subEdges_edgesList.2 ⟨_, ho, Or.inr h⟩)) hen_c hcomp_c hpe_c
          (by subst hS₂; simp [hsz]) (by subst hS₂; simp [hsdz])
          (fun k hk => by
            subst hS₂
            show (anc ++ [v])[k]? = some s.stackVerts[k]!
            rcases Nat.lt_succ_iff_lt_or_eq.1 hk with hk | rfl
            · rw [List.getElem?_append_left hk]; exact hancsv k hk
            · rw [List.getElem?_append_right (Nat.le_refl _), hC.sv_d]; simp) hC₂'
        have hK := kTree (.node y outs') (anc.length + 1) S₂ g (anc.length + 1) 0 S₂ hT₂ (Nat.le_refl _)
          (by have := hT.size; omega) (fun w _ => vertItem_ne_zero w) (fun e _ => edgeItem_ne_zero g e)
          Keep.refl
        have hsdb : wp (walkTree (.node y outs') (anc.length + 1)) (fun _ s₃ => ∀ k, k < anc.length + 1 →
            s₃.stackDir[k]! = S₂.stackDir[k]!) S₂ :=
          fun k hk => walkTree_stackDir_below _ _ S₂ k hk
        have hik : wp (walkTree (.node y outs') (anc.length + 1)) (fun _ s₃ =>
            KeptRel S₂.g S₂ (DfsTree.node y outs').verts (DfsTree.node y outs').edges s₃) S₂ :=
          wp_mono _ hfr fun _ _ h => by rw [hT₂.g_eq]; exact h.1
        have hq := walkTree_qroot_kept (.node y outs') (anc.length + 1) S₂ hik e he_ne
          (by subst hS₂; rw [hgs]; exact he_lt) (by subst hS₂; rw [hgs]; exact hT.size)
          (by subst hS₂; rw [hgs]; exact fun x hx => hw x (mem_vertsList_of_verts ho hx))
          (by subst hS₂; exact (hC.q_fresh _ (List.mem_cons_self ..) e (subEdges_e _)).1)
        have hvr := walkTree_vroot_kept (.node y outs') (anc.length + 1) S₂ hik v hv_nc
          (by subst hS₂; exact hv') (by subst hS₂; rw [hgs]; exact hT.size)
          (by subst hS₂; exact hC.v_root)
        have hbel := walkTree_below_kept (.node y outs') (anc.length + 1) S₂ hik hT₂ v hv_nc
          (by subst hS₂; exact hv') (by subst hS₂; exact hC.ch_lt)
          (fun x hx => by
            subst hS₂; exact (hC.v_fresh _ (List.mem_cons_self ..) e cls (.node y outs') rfl x hx).1)
          (fun e' he' => by
            subst hS₂; exact (hC.q_fresh _ (List.mem_cons_self ..) e' (Or.inr he')).1)
        have hboth : wp (walkTree (.node y outs') (anc.length + 1)) (fun _ s₃ =>
            FinishEar v anc.length (.tree e cls (.node y outs')) (L ++ s.tstack).length
              (hasVert || push) s₃ ∧
            wp (finishEdge v anc.length (.tree e cls (.node y outs')) (L ++ s.tstack).length
                (hasVert || push))
              (fun hv' s' => (∃ hvF', EarCtx v anc.length
                (done ++ [(.tree e cls (.node y outs'), hvF')]) rest hv' base bE sv sd s') ∧
                OutFrame v (.tree e cls (.node y outs')) base s s' ∧
                (CtxShape v anc.length done hasVert base s → ∃ hvF',
                  EarCtx v anc.length (done ++ [(.tree e cls (.node y outs'), hvF')]) rest hv'
                    base bE sv sd s' ∧
                  CtxShape v anc.length (done ++ [(.tree e cls (.node y outs'), hvF')]) hv' base s'))
                s₃) S₂ := by
          refine wp_mono _ (wp_and hendS (wp_and hK (wp_and hsdb (wp_and hq (wp_and hvr
            (wp_and hbel (wp_and hik hfr))))))) ?_
          rintro _ s₃ ⟨⟨done', hv'', sE, dir', L', push', D₃, rfl, hC', hS', hdone', hpush', hL', hsd₃,
            hdir₃⟩,
            hK3, hsdb3, hq3, hvr3, hbel3, hik3, hfr3⟩
          have H : TreeSite v anc.length done rest hasVert base bE sv sd s e cls y outs' L push done'
              hv'' _ sv' sd' sE dir' L' push' D₃ :=
            { ctx := hC, hv := hv', hy := by rw [hgs]; exact hy, e_lt := by rw [hgs]; exact he_lt
              hsd := hsd, rank := hrank, inc := hinc, nd := hnd_o, done_nc := hdone_nc, tree := ht
              cls_ret := hcls, v_nc := hv_nc, e_ne := he_ne
              comp := fun e' he' x hx hxC => by
                rw [hgs] at he' hx
                rcases hcomp_c e' he' x hx hxC with h | h
                · exact Or.inr h
                · exact Or.inl h
              ends := fun hge e' hs _ x hx =>
                ends_of_wf_boundary hwf_o hends_o hge e' hs x (by rw [hgs] at hx; exact hx)
              hpush := hpush, hL := hL, ctx' := hC', hdone' := hdone'
              inc' := fun o' ho' => by
                have hm : o'.1 ∈ (DfsTree.node y outs').outs := by
                  rw [← hdone']; exact List.mem_map_of_mem ho'
                rw [hgs]
                exact ⟨he _ (mem_subEdges_edgesList.2 ⟨_, ho, Or.inr
                  (mem_subEdges_edgesList.2 ⟨_, hm, subEdges_e _⟩)⟩), inc_e_of_ends (hEt _ hm)⟩
              hbE' := fun k hk e' _ => by rw [getElem!_map_fn _ _ k hk]
              gE := hK3.g.trans (by rw [hS₂])
              svlo := fun k hk => (hK3.svlo k (Nat.lt_succ_of_le hk)).trans (by rw [hS₂])
              sdlo := fun k hk => (hsd₃ k (Nat.lt_succ_of_le hk)).symm.trans
                ((hsdb3 k (Nat.lt_succ_of_le hk)).trans (by rw [hS₂]))
              q_root := by rw [hS₂] at hq3; exact hq3
              v_root := hvr3, hpush' := hpush', hL' := hL'
              sd₃ := fun k hk => hsd₃ k (Nat.lt_succ_of_le hk)
              size_le := by have h := hK3.size; subst hS₂; simpa [pushEnd] using h
              hsz := by rw [hgs]; exact hT.size
              verts_lt := fun x hx => by rw [hgs]; exact hw x (mem_vertsList_of_verts ho hx)
              e_ends := by rw [hgs]; exact hends_c.1
              rest_nd := fun o' ho' e' hs => by
                have hm : o' ∈ done.map (·.1) ++ DfsOut.tree e cls (.node y outs') :: rest :=
                  List.mem_append_right _ (List.mem_cons_of_mem _ ho')
                refine ⟨fun hs' => List.disjoint_of_nodup_append hnd_vo (mem_edges_of_subEdges hs')
                  (List.mem_flatMap.2 ⟨o', ho', mem_edges_of_subEdges hs⟩), ?_⟩
                rw [hgs]; exact he _ (mem_subEdges_edgesList.2 ⟨_, hm, hs⟩)
              rest_nv := fun o' ho' e₁ cls₁ c ho₁ w hw' => by
                subst ho₁
                refine ⟨fun hwC => List.disjoint_of_nodup_append
                  (List.nodup_append.1 (List.nodup_cons.1 hnd_vC).2).2.1 hwC
                  (List.mem_flatMap.2 ⟨_, ho', hw'⟩), ?_⟩
                rw [hgs]
                exact hw _ (mem_vertsList.2 ⟨_, List.mem_append_right _ (List.mem_cons_of_mem _ ho'),
                  e₁, cls₁, c, rfl, hw'⟩)
              done_rest := fun o' ho' e' hs x hx o'' ho'' e₁ cls₁ child ho₁ hxC => by
                have hm : o'.1 ∈ done.map (·.1) ++ DfsOut.tree e cls (.node y outs') :: rest :=
                  List.mem_append_left _ (List.mem_map_of_mem ho')
                have hxR : x ∈ rest.flatMap DfsOut.verts :=
                  List.mem_flatMap.2 ⟨o'', ho'', by rw [ho₁]; exact hxC⟩
                rw [hgs] at hx
                rcases endsOut_wf g anc v o'.1 (hwf _ hm) (hends _ hm) e' hs x hx with h | rfl | h
                · exact List.disjoint_of_nodup_append hndV h (List.mem_cons_of_mem _
                    (List.mem_append_right _ (List.mem_append_right _ hxR)))
                · exact (List.nodup_cons.1 hnd_vC).1
                    (List.mem_append_right _ (List.mem_append_right _ hxR))
                · exact List.disjoint_of_nodup_append (List.nodup_cons.1 hnd_vC).2
                    (List.mem_flatMap.2 ⟨o'.1, List.mem_map_of_mem ho', h⟩)
                    (List.mem_append_right _ hxR)
              items_kept := fun j hj hjv hje => by
                have h := hik3 j (by subst hS₂; exact hj) hjv (by subst hS₂; exact hje)
                subst hS₂; simpa [pushEnd] using h
              child_fresh := fun top htop t ht i hi => by
                have h := hfr3.2 (L' ++ top) (by simp [pushEnd, htop]) t (List.mem_append_right _ ht) i hi
                rcases h with h | ⟨e', he', h⟩ | h
                · exact Or.inl h
                · exact Or.inr (Or.inl ⟨e', he', by rw [hgs]; exact h⟩)
                · exact Or.inr (Or.inr (by subst hS₂; exact h))
              hdir' := fun hp => hdir₃ hp (by
                have h1 := hK3.sd
                simp only [pushEnd, hS₂, Array.size_set!] at h1
                rw [h1, hsdz]; omega)
              bridge_bd := fun hb o' ho' => EarDfs.bridge_outs_bd hwf_c.1 hwf_c.2 hb o'.1
                (by have h := List.mem_map_of_mem (f := (·.1)) ho'; rw [hdone'] at h; exact h)
              comp_ret := fun hc o' ho' => EarDfs.comp_outs_ret hwf_c.1 hwf_c.2 hc o'.1
                (by have h := List.mem_map_of_mem (f := (·.1)) ho'; rw [hdone'] at h; exact h)
              comp_ex := fun hc => by
                obtain ⟨o', ho', hl⟩ := EarDfs.comp_outs_ex hwf_c.1 hwf_c.2 hc
                have ho'' : o' ∈ done'.map (·.1) := by rw [hdone']; exact ho'
                obtain ⟨o'', hmem, hEq⟩ := List.mem_map.1 ho''
                exact ⟨o'', hmem, by rw [hEq]; exact hl⟩
              comp_t1 := fun hc o' ho' _ => EarDfs.comp_outs_t1 hwf_c.1 hwf_c.2 hc o'.1
                (by have h := List.mem_map_of_mem (f := (·.1)) ho'; rw [hdone'] at h; exact h)
              shape' := hS'
              c_reach := fun E hE e' he' => by
                rw [hgs]
                exact outs_reach g E y outs' hEt (fun e'' h =>
                  ⟨he _ (mem_subEdges_edgesList.2 ⟨_, ho, Or.inr h⟩), hE e'' h⟩) e' he'
              edges_lt := fun e' he' => by
                rw [hgs]; exact he _ (mem_subEdges_edgesList.2 ⟨_, ho, Or.inr he'⟩) }
          refine ⟨⟨fun h => ?_, earAt_tree_of_ctx H⟩, wp_mono _ (ctx_step_tree_frame H) fun _ _ h =>
            ⟨⟨_, h.1⟩, h.2.1, fun hS => ⟨_, h.1, h.2.2 hS⟩⟩⟩
          have h0 : hasVert = false := (Bool.or_eq_false_iff.1 h).1
          obtain ⟨h1, h2⟩ := hC.vert_book h0
          have hg3 : (pushEnd sE D₃ L').g = s.g := hK3.g.trans (by rw [hS₂])
          rw [hS₂] at hbel3
          rw [hg3] at hbel3 ⊢
          exact ⟨(Graph.ConnEdges.congr fun e' he' => hbel3 e' he').2 h1,
            (Graph.TwoAttached.congr fun e' he' => hbel3 e' he').2 h2⟩
        subst hS₂
        unfold OutMid
        simp only [wp_modify]
        refine ⟨⟨hET, wp_mono _ hboth fun _ _ h => h.1⟩, ?_⟩
        unfold walkOutRest
        rw [wp_bind, wp_tstackSize]
        simp only [wp_bind, wp_modify]
        exact wp_mono _ hboth fun _ _ h => h.2
    exact ⟨⟨hvb0, wp_mono _ hcore fun _ _ h => h.1⟩,
      wp_mono _ hcore fun _ _ h => wp_mono _ h.2 fun _ _ h' => h'.1,
      wp_mono _ hcore fun _ _ h => wp_mono _ h.2 fun _ _ h' => h'.2.1,
      fun hS => wp_mono _ hcore fun _ _ h => wp_mono _ h.2 fun _ _ h' => h'.2.2 hS⟩
end

/-- Admitted: the ear content of a root walk — `FinishBook.ear`/`vert` at every `finishEdge`
(PROOF.md §4.1–4.3, every field 0 violations on 3000 seeds), from the DFS well-formedness and
endpoint facts, the start state (empty tstack, `Inv' 0`, `Shape`, the arrays sized `g.nv`) and the
freshness of the tree's own `V`/`Q` items (childless and parentless). The remaining fields of
`FinishBook` are derived in `walkTree_book`. -/
theorem walkTree_ear (t : DfsTree) (s : WalkState) (hwf : t.WF []) (hends : t.Ends s.g)
    (hvlt : ∀ v ∈ t.verts, v < s.g.nv) (helt : ∀ e ∈ t.edges, e < s.g.ne)
    (hvn : t.verts.Nodup) (hen : t.edges.Nodup)
    (hcomp : ∀ e, e < s.g.ne → ∀ x, s.g.Inc e x → x ∈ t.verts → e ∈ t.edges)
    (hsv : s.stackVerts.size = s.g.nv) (hsd : s.stackDir.size = s.g.nv)
    (hfo : s.firstOccurrence.size = s.g.nv)
    (hts : s.tstack = []) (hi : s.Inv' 0) (hs : Shape s)
    (hvfresh : ∀ v ∈ t.verts, Items.ch s.items (vertItem v) = [] ∧
      ∀ p, ¬ Items.IsParent s.items p (vertItem v))
    (hefresh : ∀ e ∈ t.edges, Items.ch s.items (edgeItem s.g e) = [] ∧
      ∀ p, ¬ Items.IsParent s.items p (edgeItem s.g e)) :
    EarTree t 0 s := by
  obtain ⟨sv, sd, hC⟩ := ctx_init_root t s hwf hends hvlt helt hvn hen hsv hsd hfo hts hi hs hvfresh hefresh
  exact (cTree t 0 s s.g [] [] [] sv sd (fun _ => False) ⟨rfl, hs.size, hs.root, hs.vert, hs.edge⟩ rfl
    hwf hends (by simpa using hvn) (fun a ha => nomatch ha) hvlt helt hen
    (fun e' he' x hx hv => Or.inl (hcomp e' he' x hx hv)) (fun _ h => h.elim) hsv hsd
    (fun k hk => absurd hk (Nat.not_lt_zero _)) hC).1

/-! `BookTree` from `EarTree`: the endpoint field `ends` comes from `DfsTree.Ends`/`WF` with the
stack vertices below the current depth equal to the ancestor path, and `q` from the freshness of
the `Q` items, both carried by the frame `Keep` of `EarFrame.lean` (a `Q` item is written only at
its own `finishEdge`). -/

abbrev BTree (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  ∀ (g : Graph) (anc : List Nat), Types g s → d = anc.length → t.WF anc → t.Ends g →
    (anc ++ t.verts).Nodup → (∀ a ∈ anc, a < g.nv) → (∀ v ∈ t.verts, v < g.nv) →
    (∀ e ∈ t.edges, e < g.ne) → t.edges.Nodup → s.stackVerts.size = g.nv →
    (∀ k, k < d → anc[k]? = some s.stackVerts[k]!) →
    (∀ e ∈ t.edges, Items.ch s.items (edgeItem g e) = []) → EarTree t d s → BookTree t d s

abbrev BOuts (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  ∀ (g : Graph) (anc : List Nat), Types g s → d = anc.length → (∀ o ∈ outs, o.WF anc v) →
    (∀ o ∈ outs, DfsOut.Ends g v o) → (anc ++ v :: DfsOut.vertsList outs).Nodup →
    (∀ a ∈ anc, a < g.nv) → v < g.nv → (∀ w ∈ DfsOut.vertsList outs, w < g.nv) →
    (∀ e ∈ DfsOut.edgesList outs, e < g.ne) → (DfsOut.edgesList outs).Nodup →
    s.stackVerts.size = g.nv → (∀ k, k < d → anc[k]? = some s.stackVerts[k]!) →
    s.stackVerts[d]! = v → (∀ e ∈ DfsOut.edgesList outs, Items.ch s.items (edgeItem g e) = []) →
    EarOuts v d outs hasVert s → BookOuts v d outs hasVert s

abbrev BOut (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  ∀ (g : Graph) (anc : List Nat), Types g s → d = anc.length → o.WF anc v →
    DfsOut.Ends g v o → (anc ++ v :: o.verts).Nodup →
    (∀ a ∈ anc, a < g.nv) → v < g.nv → (∀ w ∈ o.verts, w < g.nv) →
    (∀ e ∈ o.edges, e < g.ne) → o.edges.Nodup →
    s.stackVerts.size = g.nv → (∀ k, k < d → anc[k]? = some s.stackVerts[k]!) →
    s.stackVerts[d]! = v → (∀ e ∈ o.edges, Items.ch s.items (edgeItem g e) = []) →
    EarOut v d o hasVert s → BookOut v d o hasVert s

mutual
theorem bTree : ∀ (t : DfsTree) (d : Nat) (s : WalkState), BTree t d s
  | .node v outs, d, s => by
    intro g anc hT hd hwf hends hnd hanc hvlt helt hen hsz hsvk hfresh hE
    simp only [DfsTree.verts, DfsTree.edges] at hnd hvlt helt hen hfresh
    simp only [DfsTree.WF] at hwf
    simp only [DfsTree.Ends] at hends
    unfold BookTree
    unfold EarTree at hE
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
    exact bOuts v d outs false _ g anc ⟨hT.g_eq, hT.size, hT.root, hT.vert, hT.edge⟩ hd hwf.2 hends
      hnd hanc hv (fun w hw => hvlt w (List.mem_cons_of_mem _ hw)) helt hen (by simp [hsz])
      (fun k hk => by rw [getElem!_set!_ne _ _ _ _ (by omega)]; exact hsvk k hk)
      (getElem!_set!_self _ _ _ (by rw [hsz]; exact hdlt)) hfresh hE

theorem bOuts : ∀ (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState),
    BOuts v d outs hasVert s
  | v, d, [], hasVert, s => by
    intro g anc _ _ _ _ _ _ _ _ _ _ _ _ _ _ hE
    unfold EarOuts at hE
    unfold BookOuts
    exact hE
  | v, d, o :: rest, hasVert, s => by
    intro g anc hT hd hwf hends hnd hanc hv hw he hen hsz hsvk hsvd hfresh hE
    unfold EarOuts at hE
    unfold BookOuts
    obtain ⟨hEo, hErest⟩ := hE
    rw [DfsOut.vertsList_eq, List.flatMap_cons] at hnd hw
    rw [DfsOut.edgesList_eq, List.flatMap_cons] at he hen hfresh
    simp only [List.mem_append, ← DfsOut.vertsList_eq, ← DfsOut.edgesList_eq] at hnd hw he hen hfresh
    have hnd_o : (anc ++ v :: o.verts).Nodup :=
      List.Nodup.sublist
        (List.Sublist.append_left (List.cons_sublist_cons.2 (List.sublist_append_left _ _)) _) hnd
    have hnd_r : (anc ++ v :: DfsOut.vertsList rest).Nodup :=
      List.Nodup.sublist
        (List.Sublist.append_left (List.cons_sublist_cons.2 (List.sublist_append_right _ _)) _) hnd
    have hdisj := List.disjoint_of_nodup_append hen
    refine ⟨bOut v d o hasVert s g anc hT hd (hwf o (by simp)) (hends o (by simp)) hnd_o hanc hv
      (fun w h => hw w (.inl h)) (fun e h => he e (.inl h))
      (List.Nodup.sublist (List.sublist_append_left _ _) hen) hsz hsvk hsvd
      (fun e h => hfresh e (.inl h)) hEo, ?_⟩
    have hK0 : wp (walkOut v d o hasVert) (fun _ s' => Keep (d + 1) 0 s s') s :=
      kOut v d o hasVert s g (d + 1) 0 s hT (Nat.le_refl _) (by have := hT.size; omega)
        (vertItem_ne_zero v) (fun w _ => vertItem_ne_zero w) (fun e _ => edgeItem_ne_zero g e) Keep.refl
    have hK : wp (walkOut v d o hasVert)
        (fun _ s' => ∀ e ∈ DfsOut.edgesList rest, Keep (d + 1) (edgeItem g e) s s') s :=
      wp_forall fun e he' =>
        kOut v d o hasVert s g (d + 1) (edgeItem g e) s hT (Nat.le_refl _) (edgeItem_lt (he e (.inr he')))
          (vertItem_ne_edgeItem hv e) (fun w hw' => vertItem_ne_edgeItem (hw w (.inl hw')) e)
          (fun e' he'' h => hdisj he'' (edgeItem_inj h ▸ he')) Keep.refl
    refine wp_mono _ (wp_and hErest (wp_and hK0 hK)) fun hv' s' ⟨hE', hK0, hK'⟩ => ?_
    exact bOuts v d rest hv' s' g anc (hT.of_keep hK0) hd (fun o h => hwf o (by simp [h]))
      (fun o h => hends o (by simp [h])) hnd_r hanc hv (fun w h => hw w (.inr h))
      (fun e h => he e (.inr h)) (List.Nodup.sublist (List.sublist_append_right _ _) hen)
      (hK0.sv.trans hsz) (fun k hk => by rw [hK0.svlo k (by omega)]; exact hsvk k hk)
      (by rw [hK0.svlo d (Nat.lt_succ_self d)]; exact hsvd)
      (fun e h => by rw [(hK' e h).ch]; exact hfresh e (.inr h)) hE'

theorem bOut : ∀ (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState), BOut v d o hasVert s
  | v, d, o, hasVert, s => by
    intro g anc hT hd hwf hends hnd hanc hv hw he hen hsz hsvk hsvd hfresh hE
    subst hd
    unfold EarOut at hE
    unfold BookOut
    obtain ⟨hvb, hpre⟩ := hE
    refine ⟨hvb, ?_⟩
    have hK0 : wp (walkOutPre v anc.length o hasVert) (fun _ s₁ => Keep (anc.length + 1) 0 s s₁) s :=
      keep_walkOutPre _ _ _ _ Keep.refl
    have hKE : wp (walkOutPre v anc.length o hasVert)
        (fun _ s₁ => ∀ e ∈ o.edges, Keep (anc.length + 1) (edgeItem g e) s s₁) s :=
      wp_forall fun e _ => keep_walkOutPre _ _ _ _ Keep.refl
    refine wp_mono _ (wp_and hpre (wp_and hK0 hKE)) fun hv' s₁ ⟨hE₁, hK1, hKE1⟩ => ?_
    have hT₁ := hT.of_keep hK1
    have hg₁ : s₁.g = g := hT₁.g_eq
    cases o with
    | back e dest cls =>
      simp only at hE₁ ⊢
      simp only [DfsOut.WF] at hwf
      simp only [DfsOut.Ends] at hends
      simp only [DfsOut.edges, List.mem_singleton, forall_eq] at he hfresh hKE1
      obtain ⟨i, hi, hcls⟩ := hwf
      have hi_lt : i < anc.length + 1 := by
        have := (List.getElem?_eq_some_iff.1 hi).1; simpa using this
      have hbk : cls.isTree = false := by
        cases hc : cls with
        | bridge => rw [hc] at hcls; have := classify_eq_bridge_iff.1 hcls.symm; simp at this; omega
        | component => rw [hc] at hcls; unfold classify at hcls; split at hcls <;> split at hcls <;> cases hcls
        | selfLoop => rfl
        | ret lv kind =>
          rw [hc] at hcls
          obtain ⟨_, _, rfl⟩ := classify_eq_ret_iff_back.1 hcls.symm
          rfl
      refine ⟨hg₁ ▸ hv, hg₁ ▸ he, ?_, ?_, hE₁.vert, ?_, hE₁.ear⟩
      · show Items.ch s₁.items (edgeItem s₁.g e) = []
        rw [hg₁, hKE1.ch]; exact hfresh
      · intro lv kind hret
        simp only [DfsOut.cls] at hret
        subst hret
        obtain ⟨h1, hlv, rfl⟩ := classify_eq_ret_iff_back.1 hcls.symm
        have h1' : i = lv := h1
        subst h1'
        simp only [DfsOut.cls, DfsOut.e, hbk, Bool.false_eq_true, ↓reduceIte]
        rw [hg₁, hK1.svlo i (by omega)]
        have hk := hsvk i hlv
        rw [List.getElem?_append_left hlv, hk] at hi
        obtain rfl := Option.some.inj hi
        exact hends
      · simp only [DfsOut.cls, hbk, Bool.false_eq_true, false_iff]
        rintro ⟨_, _, _, h⟩; cases h
    | tree e cls child =>
      simp only at hE₁ ⊢
      simp only [wp_modify] at hE₁ ⊢
      obtain ⟨hEc, hfin⟩ := hE₁
      simp only [DfsOut.WF] at hwf
      simp only [DfsOut.Ends] at hends
      simp only [DfsOut.verts] at hnd hw
      simp only [DfsOut.edges, List.mem_cons, forall_eq_or_imp] at he hfresh hKE1
      obtain ⟨hwfc, hcls⟩ := hwf
      obtain ⟨hpe, hendsc⟩ := hends
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
      have hT₂ : Types g { s₁ with firstOccurrence := s₁.firstOccurrence.set! anc.length s₁.g.ne } :=
        ⟨hT₁.g_eq, hT₁.size, hT₁.root, hT₁.vert, hT₁.edge⟩
      have hnd' : (anc ++ [v] ++ child.verts).Nodup := by simpa using hnd
      have hKc : wp (walkTree child (anc.length + 1))
          (fun _ s₃ => Keep (anc.length + 1) (edgeItem g e)
            { s₁ with firstOccurrence := s₁.firstOccurrence.set! anc.length s₁.g.ne } s₃) _ :=
        kTree child _ _ g _ _ _ hT₂ (Nat.le_refl _) (edgeItem_lt he.1)
          (fun w hw' => vertItem_ne_edgeItem (hw w hw') e)
          (fun e' he' h => (List.nodup_cons.1 hen).1 (edgeItem_inj h ▸ he')) Keep.refl
      refine ⟨bTree child (anc.length + 1) _ g (anc ++ [v]) hT₂ (by simp) hwfc hendsc hnd'
        (fun a ha => by
          rcases List.mem_append.1 ha with ha | ha
          · exact hanc a ha
          · exact (List.mem_singleton.1 ha) ▸ hv)
        hw he.2 (List.nodup_cons.1 hen).2 (show s₁.stackVerts.size = g.nv from hK1.sv.trans hsz)
        (fun k hk => by
          rcases Nat.lt_or_ge k anc.length with hk' | hk'
          · rw [List.getElem?_append_left hk', hsvk k hk']
            show _ = some s₁.stackVerts[k]!
            rw [hK1.svlo k (by omega)]
          · have : k = anc.length := by omega
            subst this
            rw [List.getElem?_concat_length]
            show _ = some s₁.stackVerts[anc.length]!
            rw [hK1.svlo _ (Nat.lt_succ_self _), hsvd])
        (fun e' he' => by
          show Items.ch s₁.items _ = []
          rw [(hKE1.2 e' he').ch]; exact hfresh.2 e' he') hEc, ?_⟩
      refine wp_mono _ (wp_and hfin hKc) fun _ s₃ ⟨hF, hK₃⟩ => ?_
      have hg₃ : s₃.g = g := hK₃.g.trans hg₁
      refine ⟨hg₃ ▸ hv, hg₃ ▸ he.1, ?_, ?_, hF.vert, ?_, hF.ear⟩
      · show Items.ch s₃.items (edgeItem s₃.g e) = []
        rw [hg₃, hK₃.ch]
        show Items.ch s₁.items _ = []
        rw [hKE1.1.ch]; exact hfresh.1
      · intro lv kind _
        simp only [DfsOut.cls, htr, ↓reduceIte, DfsOut.dest, DfsOut.e]
        rw [hg₃, hK₃.svlo _ (Nat.lt_succ_self _)]
        show Items.PairEq (child.v, s₁.stackVerts[anc.length]!) _
        rw [hK1.svlo _ (Nat.lt_succ_self _), hsvd]
        exact hpe
      · simp only [DfsOut.cls, htr, true_iff]
        exact ⟨e, cls, child, rfl⟩
end

/-- The bookkeeping of a root walk (`BookTree`): the ear content from `walkTree_ear`, the endpoint
and `Q`-item fields from the DFS facts and the frame. -/
theorem walkTree_book (t : DfsTree) (s : WalkState) (hwf : t.WF []) (hends : t.Ends s.g)
    (hvlt : ∀ v ∈ t.verts, v < s.g.nv) (helt : ∀ e ∈ t.edges, e < s.g.ne)
    (hvn : t.verts.Nodup) (hen : t.edges.Nodup)
    (hcomp : ∀ e, e < s.g.ne → ∀ x, s.g.Inc e x → x ∈ t.verts → e ∈ t.edges)
    (hsv : s.stackVerts.size = s.g.nv) (hsd : s.stackDir.size = s.g.nv)
    (hfo : s.firstOccurrence.size = s.g.nv)
    (hts : s.tstack = []) (hi : s.Inv' 0) (hs : Shape s)
    (hvfresh : ∀ v ∈ t.verts, Items.ch s.items (vertItem v) = [] ∧
      ∀ p, ¬ Items.IsParent s.items p (vertItem v))
    (hefresh : ∀ e ∈ t.edges, Items.ch s.items (edgeItem s.g e) = [] ∧
      ∀ p, ¬ Items.IsParent s.items p (edgeItem s.g e)) :
    BookTree t 0 s :=
  bTree t 0 s s.g [] ⟨rfl, hs.size, hs.root, hs.vert, hs.edge⟩ rfl hwf hends (by simpa using hvn)
    (fun a ha => nomatch ha) hvlt helt hen hsv (fun k hk => nomatch hk) (fun e he => (hefresh e he).1)
    (walkTree_ear t s hwf hends hvlt helt hvn hen hcomp hsv hsd hfo hts hi hs hvfresh hefresh)

/-- Admitted (ear content; PROOF.md §4.2b correction). After `finishEdge` of a tree edge at depth
`d`, run under `Inv' (d+1)`, the result satisfies `Inv' d`: every entry still attached at the child
`stackVerts[d+1]` has it as its own terminal or lies below an entry whose `vStart` it is (Loop 1
leaves the entries with `topDepth > d` with `vStart = stackVerts[d+1]`). Validated on 9000 random
multigraphs (`checks/InvCheck.lean`); replaces the false `ear_lower` (`Term d` without the context,
counterexample: the type-2 chain `2→3→4→5→6` of the cycle `0..6` with chords `6-1`, `5-2`). -/
theorem ear_lower' {curV d : Nat} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}
    (ht : o.cls.isTree = true) (hg : FinishGuards d o origTstack hasVert s) (hi : s.Inv' (d + 1))
    (hs : Shape s) (hb : FinishBook curV d o origTstack hasVert s)
    (hpost : (after (finishEdge curV d o origTstack hasVert) s).Inv' (d + 1)) :
    (after (finishEdge curV d o origTstack hasVert) s).Inv' d := by
  obtain ⟨sub, base, hlen, hE⟩ := hb.ear
  subst hlen
  obtain ⟨hg', hsv, hlow⟩ := hE.lower ht
  have hchild := hE.sv_child ht
  set fe := after (finishEdge curV d o base.length hasVert) s with hfe
  refine ⟨fun above t below hts => ?_, hpost.nodes⟩
  have h := hpost.entries above t below hts
  refine ⟨h.conn, fun v e e' he he' hEe hEe' hv hv' => ?_⟩
  have hT := h.attached v e e' he he' hEe hEe' hv hv'
  rw [hg'] at he he' hEe hEe' hv hv'
  have hdest : v = o.dest → t.Term' d fe above v := by
    rintro rfl
    rcases hlow above t below hts ⟨e, he, hEe, hv⟩ with hint | hb | ⟨t', ht', hb⟩
    · exact absurd (hint e' he' hv') hEe'
    · exact .inl (.inl hb)
    · exact .inr ⟨t', ht', .inl hb⟩
  have hk : ∀ k, k ≤ d + 1 → v = fe.stackVerts[k]! → k ≤ d ∨ v = o.dest := by
    intro k hk hvk
    rcases Nat.lt_or_ge k (d + 1) with h | h
    · exact .inl (Nat.le_of_lt_succ h)
    · right
      have : k = d + 1 := by omega
      subst this
      rw [hvk, hsv, hchild]
  rcases hT with (hvS | ⟨k, hk1, hk2, hk3⟩) | ⟨t', ht', (hvS | ⟨k, hk1, hk2, hk3⟩)⟩
  · exact .inl (.inl hvS)
  · rcases hk k hk2 hk3 with hkd | hvd
    · exact .inl (.inr ⟨k, hk1, hkd, hk3⟩)
    · exact hdest hvd
  · exact .inr ⟨t', ht', .inl hvS⟩
  · rcases hk k hk2 hk3 with hkd | hvd
    · exact .inr ⟨t', ht', .inr ⟨k, hk1, hkd, hk3⟩⟩
    · exact hdest hvd

/-! ### Boundary edges -/

/-- Recording `vs` on a childless node item (a fresh `I`/`O` leaf). -/
theorem Inv'.modifyVs_leaf (j : ItemId) (vsv : Option Nat × Option Nat) (h : s.Inv' D)
    (hj : Items.ch s.items j = []) :
    Inv' D { s with items := s.items.modify j fun it => { it with vs := vsv } } := by
  have hB : ∀ a i, Items.Below (s.items.modify j fun it => { it with vs := vsv }) a i ↔ Items.Below s.items a i :=
    fun _ _ => Items.Below_modify_ch_eq j (fun it => { it with vs := vsv }) fun _ => rfl
  refine Inv'.ofStack (h.stack.congr (s := s) rfl rfl fun _ _ e _ => TEntry.edges_congr (fun i _ e => hB i _) e)
    fun i hi hsz => ?_
  have hi : 1 + s.g.nv + s.g.ne ≤ i := hi
  have hsz' : i < s.items.size := by simpa using hsz
  by_cases hne : i = j
  · subst hne
    have hE : ∀ e, e < s.g.ne → ¬ Items.EdgeBelow s.g (s.items.modify i fun it => { it with vs := vsv }) i e := by
      intro e he hb
      rcases ((hB i _).1 hb).head_cases with heq | ⟨c, hc, _⟩
      · have : i = 1 + s.g.nv + e := heq
        omega
      · rw [Items.IsParent, hj] at hc; exact List.not_mem_nil hc
    exact ⟨Graph.ConnEdges.empty hE, fun _ _ _ => Graph.TwoAttached.empty hE⟩
  · exact ItemInv.congr (s := s) rfl (Items.vs_modify_of_ne _ _ hne) (fun e _ => hB i _) (h.nodes i hi hsz')

/-- Writing `ch` of a non-node item (`Q`, `V`) that has no parent and lies in no span changes no
entry's edge set and no node's subtree. -/
theorem Inv'.modifyCh (j : ItemId) (f : Item → Item) (h : s.Inv' D) (hj : j < 1 + s.g.nv + s.g.ne)
    (hroot : ∀ p, ¬ Items.IsParent s.items p j) (hfree : ∀ t ∈ s.tstack, j ∉ t.spans.1 ++ t.spans.2) :
    Inv' D { s with items := s.items.modify j f } := by
  refine Inv'.ofStack (h.stack.congr (s := s) rfl rfl fun t ht e _ => TEntry.edges_modify_of_not_mem j f hroot (hfree t ht) e)
    fun i hi hsz => ?_
  have hi : 1 + s.g.nv + s.g.ne ≤ i := hi
  have hsz' : i < s.items.size := by simpa using hsz
  have hne : i ≠ j := by intro h; subst h; exact absurd hi (Nat.not_le.2 hj)
  refine ItemInv.congr (s := s) rfl (Items.vs_modify_of_ne _ _ hne)
    (fun e _ => Items.Below_modify_of_not_below j f fun hb => hne (Items.Below.eq_of_no_parent hroot hb))
    (h.nodes i hi hsz')

/-- The item facts `finishBoundary` relies on: the popped entries exist, the `Q` item of the edge and
the vertex item of `curV` are roots outside every span (so writing their `ch` touches no entry and
no node). From the ear invariant: boundary edges are never pushed, and the vertex entry is pushed
only after all boundary edges of `curV`. -/
structure BoundaryOk (D curV d : Nat) (o : DfsOut) (s : WalkState) : Prop where
  pops : o.cls.isTree = true → if o.cls.lowval d == d + 1 then s.tstack ≠ [] else 2 ≤ s.tstack.length
  /-- Popping the top entry (bridge: the child's block; component: the child's block, then the
  vertex entry of `curV`) loses no attachment of the entries below: a remaining entry touching a
  terminal of the popped entry has it as its own terminal. -/
  gone : o.cls.isTree = true → ∀ t rest, s.tstack = t :: rest → ∀ u ∈ rest, ∀ v, t.Term D s v →
    s.g.Touches (u.edges s.g s.items) v → u.Term D s v
  gone₂ : o.cls.isTree = true → o.cls.lowval d ≠ d + 1 →
    ∀ t₁ t₂ rest, s.tstack = t₁ :: t₂ :: rest → ∀ u ∈ rest, ∀ v, t₂.Term D s v →
    s.g.Touches (u.edges s.g s.items) v → u.Term D s v
  q_root : ∀ p, ¬ Items.IsParent s.items p (edgeItem s.g o.e)
  q_free : ∀ t ∈ s.tstack, edgeItem s.g o.e ∉ t.spans.1 ++ t.spans.2
  v_root : ∀ p, ¬ Items.IsParent s.items p (vertItem curV)
  v_free : ∀ t ∈ s.tstack, vertItem curV ∉ t.spans.1 ++ t.spans.2

/-- The ear fact behind `BoundaryOk` (span ownership: `Q` items are pushed only by returning edges,
the vertex item only after the boundary edges; the popped block is separated from the entries below
by the articulation vertex `curV`, so its terminals are touched by none of them). -/
theorem ear_boundary {curV d : Nat} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}
    (hge : d ≤ o.cls.lowval d) (hg : FinishGuards d o origTstack hasVert s) (hi : s.Inv' D)
    (hs : Shape s) (hb : FinishBook curV d o origTstack hasVert s)
    (hD : D = if o.cls.isTree then d + 1 else d) : BoundaryOk D curV d o s := by
  obtain ⟨sub, base, hlen, hE⟩ := hb.ear
  have hts := hE.tstack
  have hnd := hE.base_touch
  have hvf : ∀ t ∈ s.tstack, vertItem curV ∉ t.spans.1 ++ t.spans.2 := by
    intro t ht hmem
    have := (hE.vert_free t ht hmem).1
    rw [hE.bd_noVert hge] at this
    exact Bool.false_ne_true this
  refine ⟨?pops, ?gone, ?gone₂, hE.q_root, hE.q_free, hE.v_root, hvf⟩
  case pops =>
    intro ht
    rw [hts]
    split
    · rename_i hL
      obtain ⟨t, hsub, -, -⟩ := hE.bd_bridge ht (by simpa using hL)
      rw [hsub]; simp
    · rename_i hL
      obtain ⟨t₁, t₂, hsub, -, -, -, -⟩ := hE.bd_comp ht hge (by simpa using hL)
      rw [hsub]; simp
  case gone =>
    intro ht t rest hts' u hu v hT htouch
    rw [if_pos ht] at hD; subst hD
    have hu' : u ∈ s.tstack.tail := by rw [hts']; exact hu
    by_cases hL : o.cls.lowval d = d + 1
    · obtain ⟨t', hsub, hbot, htop⟩ := hE.bd_bridge ht hL
      rw [hts, hsub] at hts'
      simp only [List.singleton_append, List.cons.injEq] at hts'
      obtain ⟨rfl, hrest⟩ := hts'
      subst hrest
      rcases hT with hv | ⟨k, hk1, hk2, hk3⟩
      · subst hv; rw [hbot] at htouch; exact absurd htouch (hnd ht u hu)
      · subst hk3; rw [htop] at hk1
        have : k = d + 1 := by omega
        subst this
        rw [hE.sv_child ht] at htouch; exact absurd htouch (hnd ht u hu)
    · obtain ⟨t₁, t₂, hsub, hbot₁, htop₁, hbot₂, htop₂⟩ := hE.bd_comp ht hge hL
      rw [hts, hsub] at hts'
      simp only [List.cons_append, List.nil_append, List.cons.injEq] at hts'
      obtain ⟨rfl, hrest⟩ := hts'
      subst hrest
      rcases hT with hv | ⟨k, hk1, hk2, hk3⟩
      · subst hv
        rcases List.mem_cons.1 hu with rfl | hu
        · exact .inl (by rw [hbot₁, hbot₂])
        · rw [hbot₁] at htouch; exact absurd htouch (hnd ht u hu)
      · subst hk3; rw [htop₁] at hk1
        have : k = d ∨ k = d + 1 := by omega
        rcases this with hkd | hkd
        · rw [hkd, hE.sv_d] at htouch
          rcases hE.bd_term ht hge u hu' htouch with h | h
          · exact .inl (by rw [hkd, hE.sv_d, h])
          · exact .inr ⟨d, h, Nat.le_succ d, by rw [hkd]⟩
        · rw [hkd, hE.sv_child ht] at htouch
          rcases List.mem_cons.1 hu with rfl | hu
          · exact .inl (by rw [hkd, hE.sv_child ht, hbot₂])
          · exact absurd htouch (hnd ht u hu)
  case gone₂ =>
    intro ht hL t₁ t₂ rest hts' u hu v hT htouch
    rw [if_pos ht] at hD; subst hD
    obtain ⟨t₁', t₂', hsub, hbot₁, htop₁, hbot₂, htop₂⟩ := hE.bd_comp ht hge hL
    rw [hts, hsub] at hts'
    simp only [List.cons_append, List.nil_append, List.cons.injEq] at hts'
    obtain ⟨rfl, rfl, hrest⟩ := hts'
    subst hrest
    rcases hT with hv | ⟨k, hk1, hk2, hk3⟩
    · subst hv; rw [hbot₂] at htouch; exact absurd htouch (hnd ht u hu)
    · subst hk3; rw [htop₂] at hk1
      have : k = d + 1 := by omega
      subst this
      rw [hE.sv_child ht] at htouch; exact absurd htouch (hnd ht u hu)

section Boundary

variable {curV d : Nat} {o : DfsOut}

/-- `Inv ∧ Shape`, with the `g`/`stackVerts` frame and the `ch`-equality needed to carry the
`BoundaryOk` root facts. -/
structure BStep (D : Nat) (s s' : WalkState) : Prop where
  inv : s'.Inv' D
  shape : Shape s'
  g : s'.g = s.g
  sv : s'.stackVerts = s.stackVerts
  size : s.items.size ≤ s'.items.size

theorem BStep.refl (hi : s.Inv' D) (hs : Shape s) : BStep D s s := ⟨hi, hs, rfl, rfl, Nat.le_refl _⟩

theorem BStep.trans {s₁ s₂ s₃ : WalkState} (h₁ : BStep D s₁ s₂) (h₂ : BStep D s₂ s₃) : BStep D s₁ s₃ :=
  ⟨h₂.inv, h₂.shape, h₂.g.trans h₁.g, h₂.sv.trans h₁.sv, Nat.le_trans h₁.size h₂.size⟩

theorem BStep.ofStep {v : Nat} {s' : WalkState} (st : Step D v s s') (hsz : s.items.size ≤ s'.items.size) :
    BStep D s s' := ⟨st.inv, st.shape, st.g, st.sv, hsz⟩

theorem BStep.modifyVs (hi : s.Inv' D) (hs : Shape s) (j : ItemId) (vsv : Option Nat × Option Nat)
    (hj : j < 1 + s.g.nv + s.g.ne) :
    BStep D s { s with items := s.items.modify j fun it => { it with vs := vsv } } :=
  .ofStep (Step.modifyVs (v := 0) hi hs j vsv hj) (by simp)

theorem BStep.alloc (hi : s.Inv' D) (hs : Shape s) (ty : NodeType) :
    BStep D s { s with items := s.items.push ⟨ty, (none, none), []⟩ } :=
  .ofStep (Step.alloc (v := 0) hi hs ty) (by simp)

theorem BStep.modifyVs_leaf (hi : s.Inv' D) (hs : Shape s) (j : ItemId) (vsv : Option Nat × Option Nat)
    (hj : Items.ch s.items j = []) :
    BStep D s { s with items := s.items.modify j fun it => { it with vs := vsv } } := by
  refine ⟨hi.modifyVs_leaf j vsv hj, ?_, rfl, rfl, by simp⟩
  refine hs.modify j (fun it => { it with vs := vsv }) (fun _ => rfl) fun hj' c hc => hs.ch_lt j c ?_
  rw [Items.IsParent, Items.ch_eq_getElem hj']
  exact hc

/-- Popping the top entry `t`; the `gone` hypothesis is stated over a state `s₀` with the same
graph, `stackVerts` and item subtrees (the `vs` writes of `finishBoundary` keep `Items.Below`). -/
theorem BStep.pop' (hi : s.Inv' D) (hs : Shape s) {s₀ : WalkState} (t : TEntry) (rest : List TEntry)
    (hts : s.tstack = t :: rest) (hg : s.g = s₀.g) (hsv : s.stackVerts = s₀.stackVerts)
    (hB : ∀ a i, Items.Below s.items a i ↔ Items.Below s₀.items a i)
    (hgone : ∀ u ∈ rest, ∀ v, t.Term D s₀ v → s₀.g.Touches (u.edges s₀.g s₀.items) v → u.Term D s₀ v) :
    BStep D s { s with tstack := s.tstack.tail } := by
  refine ⟨?_, hs.tstack fun t ht => hs.span t (List.mem_of_mem_tail ht), rfl, rfl, Nat.le_refl _⟩
  have hT : ∀ (t : TEntry) v, t.Term D s v ↔ t.Term D s₀ v := by
    intro t v; unfold TEntry.Term; rw [hsv]
  have hE : ∀ (u : TEntry) e, u.edges s.g s.items e ↔ u.edges s₀.g s₀.items e := by
    intro u e; rw [hg]; exact TEntry.edges_congr (fun i _ e => hB i _) e
  have := hi.pop t rest hts fun u hu v hT' htouch => by
    obtain ⟨e, he, hue, hinc⟩ := htouch
    refine (hT u v).2 (hgone u hu v ((hT t v).1 hT') ⟨e, ?_, (hE u e).1 hue, ?_⟩)
    · rwa [← hg]
    · rwa [← hg]
  simpa only [hts, List.tail_cons] using this

theorem BStep.frame (hi : s.Inv' D) (hs : Shape s) {s' : WalkState} (hg : s'.g = s.g := by rfl)
    (hsv : s'.stackVerts = s.stackVerts := by rfl) (hitems : s'.items = s.items := by rfl)
    (hts : s'.tstack = s.tstack := by rfl) : BStep D s s' :=
  ⟨hi.frame hg hsv hitems hts, hs.frame hg hitems hts, hg, hsv,
    Nat.le_of_eq (congrArg Array.size hitems).symm⟩

theorem BStep.modifyCh (hi : s.Inv' D) (hs : Shape s) (j : ItemId) (f : Item → Item)
    (hj : j < 1 + s.g.nv + s.g.ne) (hty : ∀ it, (f it).type = it.type)
    (hroot : ∀ p, ¬ Items.IsParent s.items p j) (hfree : ∀ t ∈ s.tstack, j ∉ t.spans.1 ++ t.spans.2)
    (hch : ∀ hj : j < s.items.size, ∀ c ∈ (f s.items[j]).ch, c < s.items.size) :
    BStep D s { s with items := s.items.modify j f } :=
  ⟨hi.modifyCh j f hj hroot hfree, hs.modify j f hty hch, rfl, rfl, by simp⟩

/-- `IsParent` is unchanged by `vs` writes and by pushing a childless item. -/
theorem isParent_modifyVs_iff (items : Items) (j : ItemId) (vsv : Option Nat × Option Nat) (p c : ItemId) :
    Items.IsParent (items.modify j fun it => { it with vs := vsv }) p c ↔ Items.IsParent items p c :=
  Items.IsParent_congr (Items.ch_modify_ch_eq j (fun it => { it with vs := vsv }) fun _ => rfl)

theorem isParent_push_iff (items : Items) (ty : NodeType) (p c : ItemId) :
    Items.IsParent (items.push ⟨ty, (none, none), []⟩) p c ↔ Items.IsParent items p c :=
  Items.IsParent_congr (Items.ch_push_nil _ rfl)

end Boundary

/-- Boundary edges (`lowval ≥ d`: bridges, components, self-loops) close a block via
`finishBoundary`: `Inv' D` and `Shape` are kept. (Not a `Step d curV`: the `Q` item is appended to
`vertItem curV`, so `Items.Below (vertItem curV)` grows.) -/
theorem finishBoundary_inv {curV d : Nat} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}
    (hge : d ≤ o.cls.lowval d) (hi : s.Inv' D) (hs : Shape s) (hb : FinishBook curV d o origTstack hasVert s)
    (hok : BoundaryOk D curV d o s) :
    (after (finishEdge curV d o origTstack hasVert) s).Inv' D ∧
      Shape (after (finishEdge curV d o origTstack hasVert) s) := by
  have hge' : o.cls.lowval d ≥ d := hge
  have hqlt : edgeItem s.g o.e < 1 + s.g.nv + s.g.ne := by
    show 1 + s.g.nv + o.e < _; have := hb.e_lt; omega
  have hvlt : vertItem curV < 1 + s.g.nv + s.g.ne := by
    show 1 + curV < _; have := hb.v_lt; omega
  have hvq : vertItem curV ≠ edgeItem s.g o.e := by
    show 1 + curV ≠ 1 + s.g.nv + o.e; have := hb.v_lt; omega
  have hvsz : vertItem curV < s.items.size := by
    show 1 + curV < _; have := hb.v_lt; have := hs.size; omega
  -- the Q item's `vs`, the block counter
  have b₀ : BStep D s { s with items := s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) } } :=
    BStep.modifyVs hi hs _ _ hqlt
  set s₀ := { s with items := s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) } }
    with hs₀
  have b₁ : BStep D s { s₀ with totBlocks := s₀.totBlocks + 1 } := b₀.trans (BStep.frame b₀.inv b₀.shape)
  set s₁ := { s₀ with totBlocks := s₀.totBlocks + 1 } with hs₁
  have hq_root₁ : ∀ p, ¬ Items.IsParent s₁.items p (edgeItem s.g o.e) := fun p h =>
    hok.q_root p ((isParent_modifyVs_iff _ _ _ _ _).1 h)
  have hv_root₁ : ∀ p, ¬ Items.IsParent s₁.items p (vertItem curV) := fun p h =>
    hok.v_root p ((isParent_modifyVs_iff _ _ _ _ _).1 h)
  have hts₁ : s₁.tstack = s.tstack := rfl
  have hsz₁ : s₁.items.size = s.items.size := by simp [hs₁, hs₀]
  -- the final vertex write, generic in the state after the branch
  have fin : ∀ s₂ : WalkState, BStep D s s₂ → s₂.g = s.g →
      (∀ p, ¬ Items.IsParent s₂.items p (vertItem curV)) →
      (∀ t ∈ s₂.tstack, vertItem curV ∉ t.spans.1 ++ t.spans.2) →
      let s₃ := { s₂ with items := s₂.items.modify (vertItem curV) fun it => { it with ch := it.ch ++ [edgeItem s.g o.e] } }
      s₃.Inv' D ∧ Shape s₃ := by
    intro s₂ b hg hroot hfree
    have hq : edgeItem s.g o.e < s₂.items.size := Nat.lt_of_lt_of_le hqlt (Nat.le_trans hs.size b.size)
    have b' := BStep.modifyCh b.inv b.shape (vertItem curV)
      (fun it => { it with ch := it.ch ++ [edgeItem s.g o.e] }) (by rw [hg]; exact hvlt) (fun _ => rfl) hroot hfree
      (fun hj c hc => by
        rcases List.mem_append.1 hc with hc | hc
        · exact b.shape.ch_lt _ c (by rw [Items.IsParent, Items.ch_eq_getElem hj]; exact hc)
        · rw [List.mem_singleton] at hc; subst hc; exact hq)
    exact ⟨b'.inv, b'.shape⟩
  show wp (finishEdge curV d o origTstack hasVert) (fun _ s' => s'.Inv' D ∧ Shape s') s
  rw [finishEdge_eq]
  simp only [finishEdge', wp_bind, wp_get, wp_stackDir, hge', ↓reduceIte]
  unfold finishBoundary
  simp only [wp_bind, wp_modifyItem, wp_modify, wp_ite, wp_allocItem, wp_makeVs, wp_popTstack, wp_pure]
  split
  · rename_i hT
    have hpops := hok.pops hT
    split
    · -- bridge
      rename_i hL
      simp only [hL, ↓reduceIte] at hpops
      obtain ⟨t, rest, hts⟩ : ∃ t rest, s.tstack = t :: rest := by
        match h : s.tstack, hpops with
        | t :: rest, _ => exact ⟨t, rest, rfl⟩
      have ht : t ∈ s.tstack := by rw [hts]; exact List.mem_cons_self ..
      have b₂ := b₁.trans (BStep.alloc b₁.inv b₁.shape .I)
      have b₃ := b₂.trans (BStep.modifyVs_leaf b₂.inv b₂.shape s₁.items.size
        (setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest)) (Items.ch_push_size _ rfl))
      have b₄ := b₃.trans (BStep.pop' b₃.inv b₃.shape (s₀ := s) t rest hts rfl rfl
        (fun a i => (Items.Below_modify_ch_eq s₁.items.size
            (fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) fun _ => rfl).trans
          ((Items.Below_push_nil ⟨.I, (none, none), []⟩ rfl).trans
            (Items.Below_modify_ch_eq (edgeItem s.g o.e) (fun it => { it with vs := (some curV, none) }) fun _ => rfl)))
        (hok.gone hT t rest hts))
      have hsz₄ : ((s₁.items.push ⟨.I, (none, none), []⟩).modify s₁.items.size fun it =>
          { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }).size = s.items.size + 1 := by
        simp [hs₁, hs₀]
      have hroot₄ : ∀ p, ¬ Items.IsParent ((s₁.items.push ⟨.I, (none, none), []⟩).modify s₁.items.size fun it =>
          { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) p (edgeItem s.g o.e) :=
        fun p h => hq_root₁ p ((isParent_push_iff _ _ _ _).1 ((isParent_modifyVs_iff _ _ _ _ _).1 h))
      have hvroot₄ : ∀ p, ¬ Items.IsParent ((s₁.items.push ⟨.I, (none, none), []⟩).modify s₁.items.size fun it =>
          { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) p (vertItem curV) :=
        fun p h => hv_root₁ p ((isParent_push_iff _ _ _ _).1 ((isParent_modifyVs_iff _ _ _ _ _).1 h))
      have b₅ := b₄.trans (BStep.modifyCh b₄.inv b₄.shape (edgeItem s.g o.e)
        (fun it => { it with ch := s₁.items.size :: s.tstack.head!.spans.2 }) hqlt (fun _ => rfl) hroot₄
        (fun t' ht' => hok.q_free t' (List.mem_of_mem_tail ht'))
        (fun hj c hc => by
          show c < ((s₁.items.push _).modify _ _).size
          rw [hsz₄]
          rcases List.mem_cons.1 hc with hc | hc
          · rw [hc, hsz₁]; exact Nat.lt_succ_self _
          · exact Nat.lt_succ_of_lt (hs.span t ht c (List.mem_append_right _ (by simpa [hts] using hc)))))
      refine fin _ b₅ rfl (fun p h => ?_) (fun t' ht' => hok.v_free t' (List.mem_of_mem_tail ht'))
      · rcases Items.IsParent_modify h with h | ⟨rfl, hj, hc⟩
        · exact hvroot₄ p h
        · rcases List.mem_cons.1 hc with hc | hc
          · exact absurd (hc.trans hsz₁) (Nat.ne_of_lt hvsz)
          · exact hok.v_free t ht (List.mem_append_right _ (by simpa [hts] using hc))
    · -- component
      rename_i hL
      simp only [hL, Bool.false_eq_true, ↓reduceIte] at hpops
      obtain ⟨t₁, t₂, rest, hts⟩ : ∃ t₁ t₂ rest, s.tstack = t₁ :: t₂ :: rest := by
        match h : s.tstack, hpops with
        | t₁ :: t₂ :: rest, _ => exact ⟨t₁, t₂, rest, rfl⟩
      have ht₁ : t₁ ∈ s.tstack := by rw [hts]; exact List.mem_cons_self ..
      have ht₂ : t₂ ∈ s.tstack := by rw [hts]; exact List.mem_cons_of_mem _ (List.mem_cons_self ..)
      have hB₁ : ∀ a i, Items.Below s₁.items a i ↔ Items.Below s.items a i :=
        fun a i => Items.Below_modify_ch_eq (edgeItem s.g o.e) (fun it => { it with vs := (some curV, none) }) fun _ => rfl
      have b₂ := b₁.trans (BStep.pop' b₁.inv b₁.shape (s₀ := s) t₁ (t₂ :: rest) hts rfl rfl hB₁
        (hok.gone hT t₁ (t₂ :: rest) hts))
      have b₃ := b₂.trans (BStep.pop' b₂.inv b₂.shape (s₀ := s) t₂ rest
        (by show s.tstack.tail = _; rw [hts, List.tail_cons]) rfl rfl hB₁ (hok.gone₂ hT (fun h => hL (by simp [h])) t₁ t₂ rest hts))
      have b₄ := b₃.trans (BStep.modifyCh b₃.inv b₃.shape (edgeItem s.g o.e)
        (fun it => { it with ch := s.tstack.head!.spans.1 ++ s.tstack.tail.head!.spans.2 }) hqlt (fun _ => rfl) hq_root₁
        (fun t' ht' => hok.q_free t' (List.mem_of_mem_tail (List.mem_of_mem_tail ht')))
        (fun hj c hc => by
          show c < s₁.items.size
          rw [hsz₁]
          rcases List.mem_append.1 hc with hc | hc
          · exact hs.span t₁ ht₁ c (List.mem_append_left _ (by simpa [hts] using hc))
          · exact hs.span t₂ ht₂ c (List.mem_append_right _ (by simpa [hts] using hc))))
      refine fin _ b₄ rfl (fun p h => ?_)
        (fun t' ht' => hok.v_free t' (List.mem_of_mem_tail (List.mem_of_mem_tail ht')))
      · rcases Items.IsParent_modify h with h | ⟨rfl, hj, hc⟩
        · exact hv_root₁ p h
        · rcases List.mem_append.1 hc with hc | hc
          · exact hok.v_free t₁ ht₁ (List.mem_append_left _ (by simpa [hts] using hc))
          · exact hok.v_free t₂ ht₂ (List.mem_append_right _ (by simpa [hts] using hc))
  · -- self-loop
    have b₂ := b₁.trans (BStep.frame (s' := { s₁ with totSelfLoops := s₁.totSelfLoops + 1 }) b₁.inv b₁.shape)
    set s₂ := { s₁ with totSelfLoops := s₁.totSelfLoops + 1 } with hs₂
    have hsz₂ : s₂.items.size = s.items.size := by simp [hs₂, hs₁, hs₀]
    have b₃ := b₂.trans (BStep.alloc b₂.inv b₂.shape .O)
    have b₄ := b₃.trans (BStep.modifyVs_leaf b₃.inv b₃.shape s₂.items.size (some curV, none) (Items.ch_push_size _ rfl))
    have hsz₄ : ((s₂.items.push ⟨.O, (none, none), []⟩).modify s₂.items.size fun it =>
        { it with vs := (some curV, none) }).size = s.items.size + 1 := by simp [hs₂, hs₁, hs₀]
    have hroot₄ : ∀ p, ¬ Items.IsParent ((s₂.items.push ⟨.O, (none, none), []⟩).modify s₂.items.size fun it =>
        { it with vs := (some curV, none) }) p (edgeItem s.g o.e) :=
      fun p h => hq_root₁ p ((isParent_push_iff _ _ _ _).1 ((isParent_modifyVs_iff _ _ _ _ _).1 h))
    have hvroot₄ : ∀ p, ¬ Items.IsParent ((s₂.items.push ⟨.O, (none, none), []⟩).modify s₂.items.size fun it =>
        { it with vs := (some curV, none) }) p (vertItem curV) :=
      fun p h => hv_root₁ p ((isParent_push_iff _ _ _ _).1 ((isParent_modifyVs_iff _ _ _ _ _).1 h))
    have b₅ := b₄.trans (BStep.modifyCh b₄.inv b₄.shape (edgeItem s.g o.e)
      (fun it => { it with ch := [s₂.items.size] }) hqlt (fun _ => rfl) hroot₄
      (fun t' ht' => hok.q_free t' ht')
      (fun hj c hc => by
        show c < ((s₂.items.push _).modify _ _).size
        rw [hsz₄, List.mem_singleton.1 hc, hsz₂]; exact Nat.lt_succ_self _))
    refine fin _ b₅ rfl (fun p h => ?_) (fun t' ht' => hok.v_free t' ht')
    · rcases Items.IsParent_modify h with h | ⟨rfl, hj, hc⟩
      · exact hvroot₄ p h
      · rw [List.mem_singleton] at hc
        exact absurd (hc.trans hsz₂) (Nat.ne_of_lt hvsz)

theorem walkOutPre_inv {v d : Nat} {o : DfsOut} {hasVert : Bool} (hi : s.Inv' d) (hs : Shape s)
    (hb : VertBook v hasVert s) :
    wp (walkOutPre v d o hasVert) (fun _ s₁ => s₁.Inv' d ∧ Shape s₁) s := by
  have hi' : ({ s with stackDir := s.stackDir.set! d (if o.cls.lowval d ≥ d then false else !s.stackDir[o.cls.lowval d]!) } : WalkState).Inv' d :=
    hi.frame'
  have hs' : Shape ({ s with stackDir := s.stackDir.set! d (if o.cls.lowval d ≥ d then false else !s.stackDir[o.cls.lowval d]!) } : WalkState) :=
    hs.frame'
  unfold walkOutPre
  simp only [wp_bind, wp_stackDir, wp_setStackDir, wp_ite]
  split
  · rename_i hc
    have hf : hasVert = false := by cases hasVert <;> simp_all
    obtain ⟨hv, hc, ha⟩ := hb hf
    exact ⟨(Step.pushVert hi' hs' d hv hc ha).inv, (Step.pushVert hi' hs' d hv hc ha).shape⟩
  · exact ⟨hi', hs'⟩

theorem finishEdge_step {v d : Nat} {o : DfsOut} {n : Nat} {hasVert : Bool} {D : Nat}
    (hD : D = if o.cls.isTree then d + 1 else d) (hi : s.Inv' D) (hs : Shape s)
    (hg : FinishGuards d o n hasVert s) (hb : FinishBook v d o n hasVert s) :
    (after (finishEdge v d o n hasVert) s).Inv' d ∧ Shape (after (finishEdge v d o n hasVert) s) := by
  have hdD : d ≤ D := by split at hD <;> omega
  have hst : (after (finishEdge v d o n hasVert) s).Inv' D ∧ Shape (after (finishEdge v d o n hasVert) s) := by
    by_cases hge : d ≤ o.cls.lowval d
    · exact finishBoundary_inv (origTstack := n) hge hi hs hb (ear_boundary hge hg hi hs hb hD)
    · obtain ⟨lv, kind, ho, hl⟩ := ret_of_lowval_lt (Nat.lt_of_not_le hge)
      obtain ⟨sub, base, hlen, hE⟩ := hb.ear
      have st := finishEdge_inv v d lv kind o n hasVert ho hl hb.v_lt hi hs
        (finishOk_of_guards ho hl hg hE hlen hi hs hD hb.v_lt hb.e_lt hb.q (hb.ends lv kind ho) hb.vert)
      exact ⟨st.inv, st.shape⟩
  refine ⟨?_, hst.2⟩
  by_cases ht : o.cls.isTree = true
  · rw [if_pos ht] at hD; subst hD; exact ear_lower' ht hg hi hs hb hst.1
  · rw [if_neg ht] at hD; subst hD; exact hst.1

abbrev InvTree (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  (∀ v outs, t = .node v outs → ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).Inv' d) →
  Shape s → GuardsTree t d s → BookTree t d s → wp (walkTree t d) (fun _ s' => s'.Inv' d ∧ Shape s') s

abbrev InvOuts (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  s.Inv' d → Shape s → GuardsOuts v d outs hasVert s → BookOuts v d outs hasVert s →
  wp (walkOuts v d outs hasVert) (fun hasVert' s' => s'.Inv' d ∧ Shape s' ∧ VertBook v hasVert' s') s

abbrev InvOut (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  s.Inv' d → Shape s → GuardsOut v d o hasVert s → BookOut v d o hasVert s →
  wp (walkOut v d o hasVert) (fun _ s' => s'.Inv' d ∧ Shape s') s

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
      obtain ⟨hv, hc, ha⟩ := hvb rfl
      have hi₂ : ({ s' with stackDir := s'.stackDir.set! d true } : WalkState).Inv' d := hi'.frame'
      have hs₂ : Shape ({ s' with stackDir := s'.stackDir.set! d true } : WalkState) := hs'.frame'
      exact ⟨(Step.pushVert hi₂ hs₂ d hv hc ha).inv, (Step.pushVert hi₂ hs₂ d hv hc ha).shape⟩
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

/-- `walkTree` preserves the (contextual) invariant `Inv' d`, given the ear guards and the
bookkeeping facts. (The same statement for `Inv' d` is false: at a vertex whose first out-edge is
type 2 an entry stays attached at a vertex that is only the `vStart` of an entry above it; see
`ear_lower'`.) -/
theorem walkTree_inv' (t : DfsTree) (d : Nat) (s : WalkState)
    (hi : ∀ v outs, t = .node v outs → ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).Inv' d)
    (hs : Shape s) (hg : GuardsTree t d s) (hb : BookTree t d s) :
    ((walkTree t d).run s).2.Inv' d ∧ Shape ((walkTree t d).run s).2 :=
  invTree t d s hi hs hg hb

end Walk

end WalkState
end Spqr
