import Spqr.EarSpec

/-!
# Per-ear stack shape (PROOF.md §4.3/§4.4) — WIP

`EarShape E s`: the tstack of `s` is this ear's entries over the `E.orig` entries that were there when
the ear started. It is the structural (depth/index/owner) part of the stack-shape hypotheses of
`finishEdge`; it is meant to discharge `FinishGuards` (Frame.lean) and the length/`vStart` fields of
`WalkSpec.lean`'s `FinishOk`, not its graph/span fields.
-/

namespace Spqr
open WalkM

/-- The ear being walked: return depth `lam`, depth `top` of the chain top, current chain bottom
`v` at depth `bot`, index `fo` of the ear's first edge, tstack size `orig` under the ear. -/
structure Ear where
  lam : Nat
  top : Nat
  bot : Nat
  v : Nat
  fo : Nat
  orig : Nat

/-- Stack shape relative to the ear `E` (chain `top … bot`, current bottom vertex `E.v`). -/
structure EarShape (E : Ear) (s : WalkState) : Prop where
  lam_lt : E.lam < E.top
  top_le : E.top ≤ E.bot
  split : ∃ hi lo, s.tstack = hi ++ lo ∧ lo.length = E.orig ∧
    (∀ e ∈ hi, E.lam ≤ e.topDepth ∧ E.fo ≤ e.firstIdx) ∧
    (∀ e ∈ lo, e.topDepth < E.top ∧ e.vStart ≠ E.v)
  head : ∀ e ∈ s.tstack.head?, e.vStart = E.v ∧ E.lam ≤ e.topDepth
  firstOcc : ∀ k, E.top ≤ k → k ≤ E.bot → E.fo ≤ s.firstOccurrence[k]!
  nxt : s.nxtEdgeIdx ≤ s.g.ne

theorem EarShape.below (h : EarShape E s) : ∃ hi lo, s.tstack = hi ++ lo ∧ Below E.top E.v lo := by
  obtain ⟨hi, lo, hl, -, -, hlo⟩ := h.split
  exact ⟨hi, lo, hl, hlo⟩

/-- (b) At a chain frame `(v, d)` with chain edge `o` (type 2, `lowval = E.lam`) whose subtree has been
walked and whose deeper frames have been finished, `EarShape` (with `E.v` the child at `d + 1` and
`E.orig` the frame's `origTstack`) gives the guards of this `finishEdge`. -/
theorem EarShape.finishGuards (E : Ear) (s : WalkState) (d : Nat) (o : DfsOut) (hasVert : Bool)
    (h : EarShape E s) (hbot : E.bot = d + 1) (htop : E.top ≤ d)
    (ho : o.cls.lowval d = E.lam) (ht : o.cls.isTree = true) (h2 : o.cls.isType1 = false)
    (hv : ∀ e ∈ s.tstack.head?, o.dest = e.vStart) :
    FinishGuards d o E.orig hasVert s := by
  sorry

/-- (b') `finishEdge` at a chain frame re-establishes the shape one frame up: the chain bottom becomes
`(v, d)`, and the entries above the frame collapse to the frame's pieces. -/
theorem EarShape.finishEdge (E : Ear) (s : WalkState) (v d : Nat) (o : DfsOut) (hasVert : Bool)
    (h : EarShape E s) (hbot : E.bot = d + 1) (htop : E.top ≤ d)
    (ho : o.cls.lowval d = E.lam) (ht : o.cls.isTree = true) (h2 : o.cls.isType1 = false) :
    wp (finishEdge v d o E.orig hasVert) (fun _ s' =>
      EarShape { E with bot := d, v := v } s') s := by
  sorry

/-- (c) Ear view of `walkTree_guards`: from an empty tstack the walk of a well-formed subtree keeps
`EarShape` along every chain and hence `FinishGuards` at every `finishEdge`; transported along
`walkEarTree_eq_walkTree`. -/
theorem walkEarTree_guards (t : DfsTree) (anc : List Nat) (s : WalkState) (hwf : t.WF anc)
    (hs : s.tstack = []) (hne : s.nxtEdgeIdx ≤ s.g.ne)
    (hsz : anc.length + t.height ≤ s.firstOccurrence.size) : GuardsTree t anc.length s := by
  sorry

end Spqr
