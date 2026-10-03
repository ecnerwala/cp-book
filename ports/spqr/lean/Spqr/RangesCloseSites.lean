import Spqr.RangesWalk

/-! # The close records at the six close sites of `finishEdge`

`CloseInv` preservation through one `finishEdge`, split into one named statement per block that
creates or completes an item record (PROOF.md §4.6). Each is stated at the state where its block
runs, under `CloseCtx`: what the walk induction knows at a `finishEdge` call (the range
invariant, the ear book-keeping, the frontier, and the scheduled adjacency `FinishR`). The
blocks between them (`feS₀`, `mergeLate`, the `closeVert'` merges/unwrap/retarget, `finishTail`)
only move spans and are covered by the proved `CloseInv` frame lemmas. The empirical checker
(`checks/RangesInvCheck.lean`, `closeSites`) evaluates every `CloseAt` clause at every one of
these block boundaries, under the statement's name. -/

namespace Spqr.WalkState
open WalkM

/-- The context of a `finishEdge curV d o origTstack hasVert` call of the walk at state `s`. -/
structure CloseCtx (σ : List Nat) (n D curV d : Nat) (o : DfsOut) (origTstack : Nat)
    (hasVert : Bool) (s : WalkState) : Prop where
  hD : D = if o.cls.isTree then d + 1 else d
  nodup : σ.Nodup
  lt : ∀ e ∈ σ, e < s.g.ne
  pos : σ[n]? = some o.e
  block : o.block <:+: σ
  ranges : s.RangesInv σ n D
  shape : Shape s
  guards : FinishGuards d o origTstack hasVert s
  book : FinishBook curV d o origTstack hasVert s
  frontier : Frontier (o := o) d origTstack s
  finishR : FinishR σ n curV d o origTstack hasVert s
  close : s.CloseInv

/-- The state `finishRest` (hence `finishP`) runs from: after `closeVert'` for a tree edge with
a vertex ear, after `mergeLate` for a tree edge without one, after the push for a back edge. -/
def feRest (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool) (s : WalkState) :
    WalkState :=
  if o.cls.isTree then
    if hasVert then feS₃ curV d o origTstack s else feS₂ d o s
  else feBack curV (o.cls.lowval d) d o s

variable {σ : List Nat} {n D curV d : Nat} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}
  {s : WalkState}

/-- Leaf Q at a returning tree edge: `pushEdgeTstack o.dest d o.e` after `feS₀` keeps every
record (`CloseInv.pushEdge`: distinct endpoints and both terminal `Att` witnesses of
`edgeItem o.e`, whose other incident edges are not below it). -/
theorem closeEars_closeAt (hc : CloseCtx σ n D curV d o origTstack hasVert s)
    (ht : o.cls.isTree = true) (hlow : o.cls.lowval d < d) :
    (ceS₁ o.dest d o.e (feS₀ d o s)).CloseInv := by
  sorry

/-- One iteration of loop 1: the S/P/R record created (or reopened) by `maybeUnwrapNxt` and
closed by `finishTstackTop` over the merged entry — attachment/terminal facts, the interior
V-child equivalence, two terminals of every non-V child, and the S/P/R shape clause. -/
theorem loop1Body_closeAt (hc : CloseCtx σ n D curV d o origTstack hasVert s)
    (ht : o.cls.isTree = true) (hlow : o.cls.lowval d < d) (k : Nat)
    (hk : ∀ j, j ≤ k → result (loop1Cond d)
      (iter (loop1Body d s.stackDir[d]!) j (ceS₁ o.dest d o.e (feS₀ d o s))) = true)
    (h : (iter (loop1Body d s.stackDir[d]!) k (ceS₁ o.dest d o.e (feS₀ d o s))).CloseInv) :
    (after (loop1Body d s.stackDir[d]!)
      (iter (loop1Body d s.stackDir[d]!) k (ceS₁ o.dest d o.e (feS₀ d o s)))).CloseInv := by
  sorry

/-- The `some item` arm of `closeVertTail` (`vertFinish`): the S or R record of the vertex ear,
from the state `cvS₅` after the merges and the retarget (the `none` arm is the identity). -/
theorem closeVertTail_closeAt (hc : CloseCtx σ n D curV d o origTstack hasVert s)
    (ht : o.cls.isTree = true) (hlow : o.cls.lowval d < d) (hv : hasVert = true)
    (h : (cvS₅ curV s.stackDir[d]! o.cls.isType1 origTstack (feSingle d o s) (feS₂ d o s)).CloseInv) :
    (feS₃ curV d o origTstack s).CloseInv := by
  sorry

/-- `finishP`: when `condP` holds, the new or reused P record over the merged `(curV,
stackVerts[lowval])` entry — at least two virtual edges, all equal to its terminals, no V child. -/
theorem finishP_closeAt (hc : CloseCtx σ n D curV d o origTstack hasVert s)
    (hlow : o.cls.lowval d < d) (h : (feRest curV d o origTstack hasVert s).CloseInv) :
    (after (finishP curV (o.cls.lowval d) o.cls.isType1)
      (feRest curV d o origTstack hasVert s)).CloseInv := by
  sorry

/-- Leaf Q at a returning back edge: `pushEdgeTstack curV lowval o.e` after `feS₀`
(`CloseInv.pushEdge` as for `closeEars_closeAt`; the counter writes are frames). -/
theorem finishBack_closeAt (hc : CloseCtx σ n D curV d o origTstack hasVert s)
    (ht : o.cls.isTree = false) (hlow : o.cls.lowval d < d) :
    (feBack curV (o.cls.lowval d) d o s).CloseInv := by
  sorry

/-- `finishBoundary`: the nonempty Q record of `o.e` with children `I :: t.spans.2` (bridge),
`backedge.spans.1 ++ t.spans.2` (completed block), or `[O]` (self-loop), then appended under
`vertItem curV` (`CloseInv.vertex_append`). -/
theorem finishBoundary_closeAt (hc : CloseCtx σ n D curV d o origTstack hasVert s)
    (hge : d ≤ o.cls.lowval d) :
    (after (finishEdge curV d o origTstack hasVert) s).CloseInv := by
  sorry

end Spqr.WalkState
