import Spqr.WalkSpec
import Spqr.Frame
import Spqr.EarInv

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
    (hg : FinishGuards d o origTstack hasVert s) (hE : s.EarFinish curV d o hasVert sub base) (hi : s.Inv' D) (hs : Shape s) (k : Nat)
    (hk : ∀ j, j ≤ k → result (loop1Cond d)
      (iter (Spqr.loop1Body d s.stackDir[d]!) j (ceS₁ o.dest d o.e (feS₀ d o s))) = true) :
    Loop1BodyOk D d s.stackDir[d]! (iter (Spqr.loop1Body d s.stackDir[d]!) k (ceS₁ o.dest d o.e (feS₀ d o s))) := by
  sorry

/-- Loop 2: every late merge joins entries sharing a terminal (`MergeTopOk`). -/
theorem ear_mergeLate (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hE : s.EarFinish curV d o hasVert sub base) (hi : s.Inv' D) (hs : Shape s) :
    MergeLateOk D d (feS₁ d o s) := by
  sorry

/-- The vertex close: loop 3 merges, the unwrap, the two merges, the retarget and the type-1 close. -/
theorem ear_closeVert (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hE : s.EarFinish curV d o hasVert sub base) (hi : s.Inv' D) (hs : Shape s) (hv : hasVert = true) :
    CloseVertOk D curV s.stackDir[d]! o.cls.isType1 origTstack (feSingle d o s) (feS₂ d o s) := by
  sorry

/-- The P-check after the vertex close. -/
theorem ear_finishP_vert (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hE : s.EarFinish curV d o hasVert sub base) (hi : s.Inv' D) (hs : Shape s) (hv : hasVert = true) :
    FinishPOk D curV lv o.cls.isType1 (feS₃ curV d o origTstack s) := by
  sorry

/-- The P-check of a first tree edge (no vertex entry yet). -/
theorem ear_finishP_tree (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hE : s.EarFinish curV d o hasVert sub base) (hi : s.Inv' D) (hs : Shape s) (hv : hasVert = false) :
    FinishPOk D curV lv o.cls.isType1 (feS₂ d o s) := by
  sorry

/-- The P-check of a back edge. -/
theorem getSide_setSides_self {α : Type} (dir : Bool) (a b : α) : getSide (setSides dir a b) dir = a := by
  cases dir <;> rfl

/-- `maybeUnwrapNxt` on a one-sided single-item `nxt` keeps the stack shape and every entry's edge
set (the unwrapped entry holds the children of its item, hence the same edges). -/
theorem maybeUnwrapNxt_edges {ty : NodeType} (hs : Shape s) (hty : ty ∉ [NodeType.F, .V, .Q])
    {a b : TEntry} {rest : List TEntry} (hts : s.tstack = a :: b :: rest) {i : ItemId}
    (hside : getSide b.spans (!s.stackDir[b.topDepth]!) = [])
    (hsingle : getSide b.spans s.stackDir[b.topDepth]! = [i]) :
    ∃ b', (after (maybeUnwrapNxt ty) s).tstack = a :: b' :: rest ∧
      (after (maybeUnwrapNxt ty) s).g = s.g ∧ (after (maybeUnwrapNxt ty) s).stackVerts = s.stackVerts ∧
      (after (maybeUnwrapNxt ty) s).stackDir = s.stackDir ∧
      b'.vStart = b.vStart ∧ b'.topDepth = b.topDepth ∧
      getSide b'.spans (!s.stackDir[b.topDepth]!) = [] ∧
      (∀ e, e < s.g.ne → (b'.edges s.g (after (maybeUnwrapNxt ty) s).items e ↔ b.edges s.g s.items e)) ∧
      (∀ t : TEntry, ∀ e, t.edges s.g (after (maybeUnwrapNxt ty) s).items e ↔ t.edges s.g s.items e) := by
  have hh : i = (getSide b.spans s.stackDir[b.topDepth]!).head! := by rw [hsingle]; rfl
  unfold after
  rw [maybeUnwrapNxt_run_eq ty s a b rest hts _ rfl i hh]
  have hpush : ∀ t : TEntry, ∀ e, t.edges s.g (s.items.push ⟨ty, (none, none), []⟩) e ↔ t.edges s.g s.items e :=
    fun _ e => TEntry.edges_congr (fun _ _ _ => Items.Below_push_nil _ rfl) e
  by_cases h1 : ty = .R ∨ s.ternarize = true
  · rw [if_pos h1, run_allocItem]
    exact ⟨b, hts, rfl, rfl, rfl, rfl, rfl, hside, fun e _ => hpush b e, hpush⟩
  rw [if_neg h1]
  by_cases h2 : s.items[i]!.type = ty
  · rw [if_pos h2]
    refine ⟨_, rfl, rfl, rfl, rfl, rfl, rfl, getSide_setSides_not _ _, fun e he => ?_, fun _ _ => Iff.rfl⟩
    have hib : i ∈ b.spans.1 ++ b.spans.2 := by
      rw [mem_of_getSide_nil _ b.spans hside, hsingle]; exact List.mem_singleton_self i
    have hilt : i < s.items.size := hs.span b (by simp [hts]) i hib
    have hget : s.items[i]! = s.items[i] := getElem!_pos s.items i hilt
    have hity : Items.type s.items i = ty := by rw [Items.type_eq_getElem hilt, ← hget]; exact h2
    have hie : ∀ e, e < s.g.ne → edgeItem s.g e ≠ i := fun e he =>
      hs.edgeItem_ne he (hs.node_of_type hilt hty hity)
    have hch : s.items[i]!.ch = Items.ch s.items i := by rw [Items.ch_eq_getElem hilt, hget]
    rw [TEntry.edges_single _ i hside hsingle, hch]
    exact TEntry.edges_unwrap _ b.vStart b.topDepth b.firstIdx i hie he
  · rw [if_neg h2, run_allocItem]
    exact ⟨b, hts, rfl, rfl, rfl, rfl, rfl, hside, fun e _ => hpush b e, hpush⟩

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
  obtain ⟨b', hts₃, hg₃, hsv₃, hsd₃, hbv', hbt', hside', hbE', hallE⟩ :=
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
    (hsingle : feSingle d o s = false) :
    MergeTopOk D (after (pushVertTstack curV d) (after (finishP curV lv o.cls.isType1) (feS₂ d o s))) := by
  sorry

/-- `FinishOk` from the guards, the invariant, the bookkeeping facts of the finished edge
(`he`, `hq`, `hends`, `hvert`) and the ear facts. The vertex item's connectivity/2-attachment is
transported to the tail through the `Step`s of the preceding blocks. -/
theorem finishOk_of_guards (ho : o.cls = .ret lv kind) (hlow : lv < d)
    (hg : FinishGuards d o origTstack hasVert s) (hE : s.EarFinish curV d o hasVert sub base)
    (hi : s.Inv' D) (hs : Shape s)
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
     hdD, ear_loop1 ho hlow ht hg hE hi hs⟩
  have st₀ : Step D curV s (feS₀ d o s) :=
    Step.modifyVs hi hs (edgeItem s.g o.e) _ (by show 1 + s.g.nv + o.e < _; omega)
  have hv₀ : curV < (feS₀ d o s).g.nv := by rw [st₀.g]; exact hv
  refine
    { e_lt := he
      ears := hears
      late := fun ht => ear_mergeLate ho hlow ht hg hE hi hs
      vert := fun ht hv' => ear_closeVert ho hlow ht hg hE hi hs hv'
      rest_vert := fun ht hv' => ⟨ear_finishP_vert ho hlow ht hg hE hi hs hv',
        fun h => by simp [hv'] at h, fun h => by simp [hv'] at h, fun h => by simp [hv'] at h⟩
      rest_tree := fun ht hv' => ?_
      q := fun _ => hq
      ends := fun hb => by simpa [hb] using hends
      lv_le := fun _ => by omega
      rest_back := fun hb => ?_ }
  · have st₁ : Step D curV _ (feS₁ d o s) := Step.closeEars st₀.inv st₀.shape hv₀ (hears ht)
    have hv₁ : curV < (feS₁ d o s).g.nv := by rw [st₁.g]; exact hv₀
    have st₂ : Step D curV _ (feS₂ d o s) :=
      Step.mergeLate st₁.inv st₁.shape hv₁ (ear_mergeLate ho hlow ht hg hE hi hs)
    have hv₂ : curV < (feS₂ d o s).g.nv := by rw [st₂.g]; exact hv₁
    have hp := ear_finishP_tree (curV := curV) ho hlow ht hg hE hi hs hv'
    have st₃ := Step.finishP (v := curV) st₂.inv st₂.shape hv₂ hp
    have st := st₀.trans (st₁.trans (st₂.trans st₃))
    obtain ⟨hc, ha⟩ := hvert hv'
    exact ⟨hp, finishTailOk_of_vert (fun _ => (st.vertTransport hc ha).1) (fun _ => (st.vertTransport hc ha).2)
      fun _ => ear_tail_tree ho hlow ht hg hE hi hs hv'⟩
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

/-- Bookkeeping of the walk, to be discharged from the DFS well-formedness and endpoint facts
(`DfsTree.WF`, the `dfsForest` endpoint lemma), the fresh-items start state and the ear invariant
(the `VertBook` connectivity/2-attachment of the vertex item after its boundary edges needs the
popped block to be a connected piece attached only at `v`); cf. `walkTree_guards`. -/
theorem walkTree_book (t : DfsTree) (d : Nat) (s : WalkState)
    (hfresh : ∀ i, i < 1 + s.g.nv + s.g.ne → Items.ch s.items i = []) :
    BookTree t d s := by
  sorry

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
  sorry

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
  gone : ∀ t rest, s.tstack = t :: rest → ∀ u ∈ rest, ∀ v, t.Term D s v →
    s.g.Touches (u.edges s.g s.items) v → u.Term D s v
  gone₂ : ∀ t₁ t₂ rest, s.tstack = t₁ :: t₂ :: rest → ∀ u ∈ rest, ∀ v, t₂.Term D s v →
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
    (hs : Shape s) (hb : FinishBook curV d o origTstack hasVert s) : BoundaryOk D curV d o s := by
  sorry

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
  ⟨hi.frame hg hsv hitems hts, hs.frame hg hitems hts, hg, hsv, by rw [hitems]; exact Nat.le_refl _⟩

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
        (hok.gone t rest hts))
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
        (hok.gone t₁ (t₂ :: rest) hts))
      have b₃ := b₂.trans (BStep.pop' b₂.inv b₂.shape (s₀ := s) t₂ rest
        (by show s.tstack.tail = _; rw [hts, List.tail_cons]) rfl rfl hB₁ (hok.gone₂ t₁ t₂ rest hts))
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
    · exact finishBoundary_inv (origTstack := n) hge hi hs hb (ear_boundary hge hg hi hs hb)
    · obtain ⟨lv, kind, ho, hl⟩ := ret_of_lowval_lt (Nat.lt_of_not_le hge)
      obtain ⟨sub, base, -, hE⟩ := hb.ear
      have st := finishEdge_inv v d lv kind o n hasVert ho hl hb.v_lt hi hs
        (finishOk_of_guards ho hl hg hE hi hs hD hb.v_lt hb.e_lt hb.q (hb.ends lv kind ho) hb.vert)
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
