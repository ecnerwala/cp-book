import Spqr.EarSpec
import Spqr.WalkInv

/-!
# The tstack guards from the ear contract

`FinishGuards` (the shape `finishEdge` pops: boundary entries, `Inv1` for loop 1, `Inv2` for
loop 2, the vertex ear's three entries) follows from the ear contract `EarAt` read at the same
site: `finishGuards_of_ear`. `GuardsTree` then follows from `BookTree` by the walk induction
(`guards_of_book`, the shape of `walkTree_frontiers`): `walkTree_guards'`.

The predecessor's standalone shape `EarShape` (depth/first-index bounds on a `hi ++ lo` split) is
too weak for `FinishGuards`: `EarShape.finishGuards` as stated is false (`EarShape.finishGuards_false`:
a single entry returning above `d` satisfies every field, but `Inv1 d` needs an entry below `d`).
The shape the walk carries is `EarFinish` (`Spqr/EarInv.lean`), via `BookTree`.
-/

namespace Spqr
open WalkM

/-- The depth-range/first-index shape of an active ear (`lam < top ≤ bot`, the stack split, head
ownership, first-occurrence bounds). Kept for the counterexample below. -/
structure Ear where
  lam : Nat
  top : Nat
  bot : Nat
  v : Nat
  fo : Nat
  orig : Nat

structure EarShape (E : Ear) (s : WalkState) : Prop where
  lam_lt : E.lam < E.top
  top_le : E.top ≤ E.bot
  split : ∃ hi lo, s.tstack = hi ++ lo ∧ lo.length = E.orig ∧
    (∀ e ∈ hi, E.lam ≤ e.topDepth ∧ E.fo ≤ e.firstIdx) ∧
    (∀ e ∈ lo, e.topDepth < E.top ∧ e.vStart ≠ E.v)
  head : ∀ e ∈ s.tstack.head?, e.vStart = E.v ∧ E.lam ≤ e.topDepth
  firstOcc : ∀ k, E.top ≤ k → k ≤ E.bot → E.fo ≤ s.firstOccurrence[k]!
  nxt : s.nxtEdgeIdx ≤ s.g.ne

theorem EarShape.below (h : EarShape E s) :
    ∃ hi lo, s.tstack = hi ++ lo ∧ Below E.top E.v lo := by
  obtain ⟨hi, lo, hl, -, -, hlo⟩ := h.split
  exact ⟨hi, lo, hl, hlo⟩

/-- `EarShape.finishGuards` (`FinishGuards` from `EarShape` + `E.bot = d + 1`, `E.top ≤ d`,
`lowval = E.lam`, a type-2 tree edge whose child is the head's `vStart`) is false: the stack
`[(0, d + 1, 0, _)]` at `d = 1` satisfies `EarShape ⟨0, 1, 2, 0, 0, 0⟩`, but `Inv1 1` needs an
entry with `topDepth < 1`. -/
theorem EarShape.finishGuards_false :
    ∃ (E : Ear) (s : WalkState) (d : Nat) (o : DfsOut) (hasVert : Bool),
      EarShape E s ∧ E.bot = d + 1 ∧ E.top ≤ d ∧ o.cls.lowval d = E.lam ∧ o.cls.isTree = true ∧
      o.cls.isType1 = false ∧ (∀ e ∈ s.tstack.head?, o.dest = e.vStart) ∧
      ¬ FinishGuards d o E.orig hasVert s := by
  let t : TEntry := ⟨0, 2, 0, ([], [])⟩
  let s : WalkState :=
    { g := ⟨1, #[]⟩
      ternarize := false
      items := #[]
      stackVerts := #[0, 0]
      stackDir := #[false, false]
      nxtEdgeIdx := 0
      firstOccurrence := #[0, 0, 0]
      tstack := [t] }
  refine ⟨⟨0, 1, 2, 0, 0, 0⟩, s, 1, .tree 0 (.ret 0 .type2Child) (.node 0 []), false,
    ⟨by decide, by decide, ⟨[t], [], rfl, rfl, by simp [t], by simp⟩, by simp [s, t], fun _ _ _ => Nat.zero_le _,
      Nat.zero_le _⟩,
    rfl, Nat.le_refl _, rfl, rfl, rfl, by simp [s, t, DfsOut.dest, DfsTree.v], fun h => ?_⟩
  have h1 : Inv1 1 (WalkState.feS₀ 1 (.tree 0 (.ret 0 .type2Child) (.node 0 [])) s).tstack := (h.2 (by decide) rfl).1
  obtain ⟨hi, b, lo, hl, hb, -⟩ := h1
  change [t] = hi ++ b :: lo at hl
  match hi, hl with
  | [], hl => cases hl; exact absurd hb (by decide)
  | [_], hl => simp at hl
  | _ :: _ :: _, hl => simp at hl

namespace WalkState

/-- Loop 1's range is `Units d`: every reached split continues with a depth-`d` entry (closed alone)
or a deeper one followed by another entry (`Loop1Spec`). -/
theorem units_of_reach {d : Nat} {o : DfsOut} {s : WalkState} {hi : List TEntry}
    (hhi : ∀ t ∈ hi, d ≤ t.topDepth) (hspec : Loop1Spec d o s hi) :
    ∀ (n : Nat) (rest done : List TEntry), rest.length ≤ n → L1Reach d hi done rest → Units d rest
  | _, [], _, _, _ => .nil
  | 0, _ :: _, _, hn, _ => absurd hn (by simp)
  | n + 1, t :: rest, done, hn, hr => by
    have hmem : ∀ u ∈ t :: rest, u ∈ hi := fun u hu => by
      rw [hr.append]; exact List.mem_append_right _ hu
    by_cases htd : t.topDepth = d
    · exact .single htd (units_of_reach hhi hspec n rest (done ++ [t]) (by simp at hn; omega) (.close hr htd))
    · have hlt : d < t.topDepth := Nat.lt_of_le_of_ne (hhi t (hmem t (by simp))) (Ne.symm htd)
      obtain ⟨-, t', rest', rfl, -, -⟩ := (hspec done t rest hr).2 hlt
      exact .pair hlt (hhi t' (hmem t' (by simp)))
        (units_of_reach hhi hspec n rest' (done ++ [t, t']) (by simp at hn; omega) (.series hr hlt))

/-- `Inv1 d` of the stack at a returning tree edge: loop 1's range (`Units d`) over the first entry
returning above `d` (which exists: the `(y, lowval)` piece). -/
theorem inv1_of_ear {curV d : Nat} {o : DfsOut} {hasVert : Bool} {sub base : List TEntry} {s : WalkState}
    (hE : s.EarFinish curV d o hasVert sub base) (ht : o.cls.isTree = true) (hlow : o.cls.lowval d < d) :
    Inv1 d s.tstack := by
  have hhi : ∀ t ∈ sub.takeWhile (fun t => decide (d ≤ t.topDepth)), d ≤ t.topDepth := fun t ht => by
    simpa using List.mem_takeWhile' ht
  have hlo : ∀ t ∈ (sub.dropWhile fun t => decide (d ≤ t.topDepth)).head?, t.topDepth < d := by
    intro t ht
    have := List.head?_dropWhile_not (fun t : TEntry => decide (d ≤ t.topDepth)) sub
    rw [Option.mem_def] at ht
    rw [ht] at this
    simpa using this
  have hrange : Loop1Range d sub (sub.takeWhile fun t => decide (d ≤ t.topDepth)) :=
    ⟨_, List.takeWhile_append_dropWhile.symm, hhi, hlo⟩
  have hspec := hE.loop1 ht hlow _ hrange
  obtain ⟨l, lo', hlo'⟩ : ∃ l lo', sub.dropWhile (fun t => decide (d ≤ t.topDepth)) = l :: lo' := by
    cases h : sub.dropWhile (fun t => decide (d ≤ t.topDepth)) with
    | cons l lo' => exact ⟨l, lo', rfl⟩
    | nil =>
      exfalso
      obtain ⟨mid, py, vy, hsub, hB⟩ := hE.bottom ht hlow
      have hsub' : sub.takeWhile (fun t => decide (d ≤ t.topDepth)) = sub := by
        have := @List.takeWhile_append_dropWhile _ (fun t => decide (d ≤ t.topDepth)) sub
        rwa [h, List.append_nil] at this
      have hpy : py ∈ sub.takeWhile fun t => decide (d ≤ t.topDepth) := by rw [hsub', hsub]; simp
      have := hhi py hpy
      rw [hB.py_top] at this
      omega
  refine ⟨_, l, lo' ++ base, ?_, hlo l (by rw [hlo']; rfl),
    units_of_reach hhi hspec _ _ _ (Nat.le_refl _) (L1Reach.nil _)⟩
  have h := @List.takeWhile_append_dropWhile _ (fun t => decide (d ≤ t.topDepth)) sub
  rw [hlo'] at h
  rw [hE.tstack]
  conv_lhs => rw [← h]
  simp

theorem feS₀_tstack (d : Nat) (o : DfsOut) (s : WalkState) : (feS₀ d o s).tstack = s.tstack := rfl

/-- The tstack guards of `finishEdge` from the ear contract at the same site. -/
theorem finishGuards_of_ear {curV d origTstack : Nat} {o : DfsOut} {hasVert : Bool} {s : WalkState}
    (hb : s.EarAt curV d o origTstack hasVert) : FinishGuards d o origTstack hasVert s := by
  obtain ⟨sub, base, hlen, hE⟩ := hb
  refine ⟨fun hge ht => ?_, fun hlt ht => ?_⟩
  · split
    · rename_i h
      obtain ⟨t, hsub, -, -⟩ := hE.bd_bridge ht (by simpa using h)
      rw [hE.tstack, hsub]; simp
    · rename_i h
      obtain ⟨t₁, t₂, hsub, -⟩ := hE.bd_comp ht hge (by simpa using h)
      rw [hE.tstack, hsub]; simp
  · show Inv1 d (feS₀ d o s).tstack ∧ (Inv2 (feS₁ d o s).firstOccurrence[d]! (feS₁ d o s).tstack ∧
      (hasVert = true → 3 ≤ (feS₂ d o s).tstack.length ∧
        (o.cls.isType1 = false → origTstack + 3 ≤ (feS₂ d o s).tstack.length)))
    refine ⟨by rw [feS₀_tstack]; exact inv1_of_ear hE ht hlt, hE.late_fo ht hlt, fun _ => ?_⟩
    obtain ⟨c, mid, py, vy, hts, -⟩ := hE.loops ht hlt
    rw [hts]
    simp only [List.length_append, List.length_cons, List.length_nil, hlen]
    omega

/-! ### `GuardsTree` from `BookTree` (the shape of `walkTree_frontiers`) -/

mutual
theorem gbTree : ∀ (t : DfsTree) (d : Nat) (s : WalkState), BookTree t d s → GuardsTree t d s
  | .node v outs, d, s => fun hb => by
    unfold GuardsTree; unfold BookTree at hb
    exact gbOuts v d outs false _ hb

theorem gbOuts : ∀ (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState),
    BookOuts v d outs hasVert s → GuardsOuts v d outs hasVert s
  | v, d, [], hasVert, s => fun _ => by unfold GuardsOuts; trivial
  | v, d, o :: rest, hasVert, s => fun hb => by
    unfold GuardsOuts; unfold BookOuts at hb
    exact ⟨gbOut v d o hasVert s hb.1,
      wp_imp (wp_of_forall fun hv' s' hb' => gbOuts v d rest hv' s' hb') hb.2⟩

theorem gbOut : ∀ (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState),
    BookOut v d o hasVert s → GuardsOut v d o hasVert s
  | v, d, o, hasVert, s => fun hb => by
    unfold GuardsOut; unfold BookOut at hb
    refine wp_imp (wp_of_forall fun hv' s₁ hb₁ => ?_) hb.2
    cases o with
    | back e cls dest => exact finishGuards_of_ear hb₁.ear
    | tree e cls child =>
      try simp only [wp_modify] at hb₁
      simp only [wp_modify]
      exact ⟨gbTree child (d + 1) _ hb₁.1,
        wp_imp (wp_of_forall fun _ s₃ hb₃ => finishGuards_of_ear hb₃.ear) hb₁.2⟩
end

/-- The tstack guards hold at every `finishEdge` of the walk of `t`, given the bookkeeping facts
(`FinishBook.ear` at every site). This is `EarSpec.walkTree_guards`/`walkEarTree_guards` with the
walk invariant made explicit: from an empty tstack and the DFS hypotheses alone the guards are
exactly the (admitted) ear contract `walkTree_book`. -/
theorem walkTree_guards' (t : DfsTree) (d : Nat) (s : WalkState) (hb : BookTree t d s) :
    GuardsTree t d s :=
  gbTree t d s hb

end WalkState

/-- The walk invariant (PROOF.md §4) at a root: along the walk of a DFS tree from the start state of
a root walk, the tstack-shape guards of every `finishEdge` (`FinishGuards`: boundary pops find their
entries, `Inv1`/`Inv2` for the merge loops, the vertex ear has its three entries) hold. The
hypotheses are those of `WalkState.walkTree_book` (the DFS facts, the start state and the freshness of
the tree's own `V`/`Q` items); the guards come from the ear contract (`finishGuards_of_ear`). -/
theorem walkTree_guards (t : DfsTree) (s : WalkState) (hwf : t.WF []) (hends : t.Ends s.g)
    (hvlt : ∀ v ∈ t.verts, v < s.g.nv) (helt : ∀ e ∈ t.edges, e < s.g.ne)
    (hvn : t.verts.Nodup) (hen : t.edges.Nodup)
    (hsv : s.stackVerts.size = s.g.nv) (hsd : s.stackDir.size = s.g.nv)
    (hfo : s.firstOccurrence.size = s.g.nv)
    (hts : s.tstack = []) (hi : s.Inv' 0) (hs : WalkState.Shape s)
    (hvfresh : ∀ v ∈ t.verts, Items.ch s.items (vertItem v) = [] ∧
      ∀ p, ¬ Items.IsParent s.items p (vertItem v))
    (hefresh : ∀ e ∈ t.edges, Items.ch s.items (edgeItem s.g e) = [] ∧
      ∀ p, ¬ Items.IsParent s.items p (edgeItem s.g e)) :
    GuardsTree t 0 s :=
  WalkState.walkTree_guards' t 0 s
    (WalkState.walkTree_book t s hwf hends hvlt helt hvn hen hsv hsd hfo hts hi hs hvfresh hefresh)

end Spqr
