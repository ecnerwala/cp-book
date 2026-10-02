import Spqr.WalkSpec

/-!
# Ear invariant (PROOF.md §4.1–4.3): content of the stack at a `finishEdge`

The corrected attachment invariant `WalkState.Inv'` and its primitive steps live in
`Spqr/WalkSpec.lean`. This file states the *ear content* of the tstack at the moment `finishEdge
curV d o origTstack hasVert` runs (`WalkState.EarFinish`), i.e. the facts beyond `FinishGuards`
(length/`topDepth` shape, `Frame.lean`) and `Inv'` that the per-block hypotheses of
`WalkSpec.FinishOk` and `WalkInv.BoundaryOk` need. Nothing here is proved yet: `EarFinish` is the
statement the preservation induction (`EarShape.finishEdge`/`walkEarTree_guards`) has to carry, and
the `ear_*` admissions of `WalkInv.lean` are to be derived from it.

Consumers, field by field (`sub` = entries created by the child's walk, `base` = the `origTstack`
entries that were there when the out-edge was entered):
* `tstack`, `back_nil`, `base_top`, `base_bot`: the `Below d curV`/`Inv1` shapes of `FinishGuards`.
* `disj`: `MergeTopOk`/`RetargetOk.disj` (entries below a merge are edge-disjoint from it) at every
  merge of loops 1–3, the vertex close and the P-check; `RInvAt.disj`.
* `sub_edges`, `base_disj`, `sub_cover`: span ownership — the entries with `topDepth ≥ d` hold
  exactly the tree edge plus the child's subtree edges (`loop1_rBranch`'s `interior`/`cur_c`,
  `FinishTopOk.mid` of loop 1's closes: the child is finished, so it is interior to the union);
  the entries below hold none of them (`ear_lower'`: after loop 1 nothing below `cur` is attached
  at `stackVerts[d+1]`).
* `loop1_bot`: `topDepth > d ⇒ vStart = o.dest` in loop 1's range (`loop1_rBranch`'s `cur_c`, with
  `mergeInto` keeping the S entry's `vStart`; `MergeOk.bottom` of the S merges).
* `loop1_side`: the `origTstack`-indexed side fact (`walk_sides`' `CloseOK`/`UnwrapOK` at loop 1's
  closes; `FinishTopOk.side`); together with `StFrame.walkTree_stackDir_below` and
  `StEar.chain_stackDir_step` it is the "whole ear on one side" fact.
* `touch_bot`, `touch_top`: `MergeOk.share` of every merge (adjacent entries share a terminal:
  `cur.vStart = nxt.vStart` in the P case, `stackVerts[d]` in the R/late cases, the chain vertex in
  the S case).
* `vert`, `vert_free`, `q_free`, `q_root`, `v_root`: `UnwrapAt.fresh`/`root`, `ItemFree`,
  `BoundaryOk.q_root/q_free/v_root/v_free`, `FinishTailOk`'s vertex item.
* `boundary`: `BoundaryOk.gone/gone₂` (the popped block is separated from the rest by `curV`).

Not covered (needs the chain structure of §4.3 and is left for the preservation session): the exact
three-entry shape `[cur, (y, l) back-edge piece, V y]` above `origTstack` at a type-1 vertex close
(`CloseVertOk.merge₁/merge₂/retarget.old`) and the `firstIdx` order behind loop 2's `MergeOk`.
-/

namespace Spqr
namespace WalkState

/-- The edges of the sub-ear finished by the tree edge `o`: the tree edge and the child's subtree
edges (tree and back). For a back edge: the edge itself. -/
def subEdges (o : DfsOut) (e : Nat) : Prop :=
  e = o.e ∨ match o with
    | .tree _ _ child => e ∈ child.edges
    | .back .. => False

/-- Loop 1's range: the top segment of `sub` whose entries have `topDepth ≥ d`, down to (excluding)
the first entry returning above `d`. -/
def Loop1Range (d : Nat) (sub hi : List TEntry) : Prop :=
  ∃ lo, sub = hi ++ lo ∧ (∀ t ∈ hi, d ≤ t.topDepth) ∧ (∀ t ∈ lo.head?, t.topDepth < d)

/-- Ear content of the tstack `sub ++ base` when `finishEdge curV d o _ hasVert` runs (`o` an
out-edge of `curV` at depth `d`, `stackVerts[d] = curV`, and `stackVerts[d+1] = o.dest` for a tree
edge). -/
structure EarFinish (curV d : Nat) (o : DfsOut) (hasVert : Bool) (sub base : List TEntry)
    (s : WalkState) : Prop where
  tstack : s.tstack = sub ++ base
  /-- A back edge creates no entries before `finishEdge`. -/
  back_nil : o.cls.isTree = false → sub = []
  /-- The enclosing entries return to `curV` or above and were not started at the child. -/
  base_top : ∀ t ∈ base, t.topDepth ≤ d
  base_bot : o.cls.isTree = true → ∀ t ∈ base, t.vStart ≠ o.dest
  /-- Open entries are pairwise edge-disjoint. -/
  disj : s.tstack.Pairwise fun t t' => ∀ e, e < s.g.ne → t.edges s.g s.items e → ¬ t'.edges s.g s.items e
  /-- Span ownership: the child's entries hold only sub-ear edges, the enclosing ones none, and
  every subtree edge is held by some child entry. -/
  sub_edges : ∀ t ∈ sub, ∀ e, e < s.g.ne → t.edges s.g s.items e → subEdges o e
  base_disj : ∀ t ∈ base, ∀ e, e < s.g.ne → t.edges s.g s.items e → ¬ subEdges o e
  sub_cover : ∀ e, e < s.g.ne → e ≠ o.e → subEdges o e → ∃ t ∈ sub, t.edges s.g s.items e
  /-- In loop 1's range an entry returning strictly below `curV` hangs at the child. -/
  loop1_bot : ∀ hi, Loop1Range d sub hi → ∀ t ∈ hi, d < t.topDepth → t.vStart = o.dest
  /-- Loop 1's range lies on the side `stackDir[d]` (the ear's side). -/
  loop1_side : ∀ hi, Loop1Range d sub hi → ∀ t ∈ hi, getSide t.spans (!s.stackDir[d]!) = []
  /-- A non-empty entry touches both its terminals (its top only while it is on the open path). -/
  touch_bot : ∀ t ∈ s.tstack, (∃ e, e < s.g.ne ∧ t.edges s.g s.items e) →
    s.g.Touches (t.edges s.g s.items) t.vStart
  touch_top : ∀ t ∈ s.tstack, (∃ e, e < s.g.ne ∧ t.edges s.g s.items e) → t.topDepth ≤ d + 1 →
    s.g.Touches (t.edges s.g s.items) s.stackVerts[t.topDepth]!
  /-- The vertex item of `curV` is on the stack iff `hasVert` (in an enclosing entry started at
  `curV`), and never elsewhere; the edge's `Q` item is on no entry; both are roots. -/
  vert : hasVert = true → ∃ t ∈ base, t.vStart = curV ∧ vertItem curV ∈ t.spans.1 ++ t.spans.2
  vert_free : ∀ t ∈ s.tstack, vertItem curV ∈ t.spans.1 ++ t.spans.2 → hasVert = true ∧ t ∈ base ∧ t.vStart = curV
  q_free : ∀ t ∈ s.tstack, edgeItem s.g o.e ∉ t.spans.1 ++ t.spans.2
  q_root : ∀ p, ¬ Items.IsParent s.items p (edgeItem s.g o.e)
  v_root : ∀ p, ¬ Items.IsParent s.items p (vertItem curV)
  /-- A block boundary: the child's block meets the rest only at the articulation vertex `curV`. -/
  boundary : d ≤ o.cls.lowval d → ∀ t ∈ sub, ∀ u ∈ base, ∀ v,
    s.g.Touches (t.edges s.g s.items) v → s.g.Touches (u.edges s.g s.items) v → v = curV

/-- `EarFinish` with the split existentially quantified; `origTstack` is the size of `base`. -/
def EarAt (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool) (s : WalkState) : Prop :=
  ∃ sub base, base.length = origTstack ∧ s.EarFinish curV d o hasVert sub base

theorem EarAt.length {curV d origTstack : Nat} {o : DfsOut} {hasVert : Bool} {s : WalkState}
    (h : s.EarAt curV d o origTstack hasVert) : origTstack ≤ s.tstack.length := by
  obtain ⟨sub, base, hl, he⟩ := h
  rw [he.tstack, List.length_append, hl]; omega

/-- A back edge: the stack is exactly `base`. -/
theorem EarAt.back {curV d origTstack : Nat} {o : DfsOut} {hasVert : Bool} {s : WalkState}
    (h : s.EarAt curV d o origTstack hasVert) (hb : o.cls.isTree = false) :
    s.tstack.length = origTstack ∧ ∀ t ∈ s.tstack, t.topDepth ≤ d := by
  obtain ⟨sub, base, hl, he⟩ := h
  have hs := he.back_nil hb
  subst hs
  rw [he.tstack]
  exact ⟨by simpa using hl, fun t ht => he.base_top t (by simpa using ht)⟩

end WalkState
end Spqr
