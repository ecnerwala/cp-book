import Spqr.WalkSpec

/-!
# Ear invariant (PROOF.md §4.1–4.3): content of the stack at a `finishEdge`

The corrected attachment invariant `WalkState.Inv'` and its primitive steps live in
`Spqr/WalkSpec.lean`. This file states the *ear content* of the tstack at the moment `finishEdge
curV d o origTstack hasVert` runs (`WalkState.EarFinish`), i.e. the facts beyond `FinishGuards`
(length/`topDepth` shape, `Frame.lean`) and `Inv'` that the per-block hypotheses of
`WalkSpec.FinishOk` and `WalkInv.BoundaryOk` need. `EarFinish` is the statement the preservation
induction (`EarShape.finishEdge`/`walkEarTree_guards`) has to carry, and the `ear_*` admissions of
`WalkInv.lean` are to be derived from it; every field has 0 violations on 3000 random multigraphs
(`checks/EarCheck.lean`).

Consumers, field by field (`sub` = entries created by the child's walk, `base` = the `origTstack`
entries that were there when the out-edge was entered):
* `tstack`, `back_nil`, `base_bot`, `sub_bot`: the `Below d curV`/`Inv1` shapes of `FinishGuards`,
  the P-check condition (`condP` looks at `nxt.vStart = curV`).
* `path`: the open path has distinct vertices (`FinishTopOk.mid`: an edge to `stackVerts[lv]` is
  not at `stackVerts[k]` for `lv < k ≤ d`).
* `disj`, `span_disj`: `MergeTopOk`/`RetargetOk.disj` (entries below a merge are edge-disjoint from
  it) at every merge of loops 1–3, the vertex close and the P-check; `UnwrapAt.fresh`; `RInvAt.disj`.
* `sub_edges`, `base_disj`, `sub_cover`: span ownership — the entries with `topDepth ≥ d` hold
  exactly the tree edge plus the child's subtree edges (`loop1_rBranch`'s `interior`/`cur_c`,
  `FinishTopOk.mid` of loop 1's closes: the child is finished, so it is interior to the union);
  the entries below hold none of them (`ear_lower'`: after loop 1 nothing below `cur` is attached
  at `stackVerts[d+1]`).
* `loop1_side`: the `origTstack`-indexed side fact (`walk_sides`' `CloseOK`/`UnwrapOK` at loop 1's
  closes; `FinishTopOk.side`); together with `StFrame.walkTree_stackDir_below` and
  `StEar.chain_stackDir_step` it is the "whole ear on one side" fact.
* `touch_bot`, `loop1_touch`: `MergeOk.share` of every merge (adjacent entries share a terminal:
  `cur.vStart = nxt.vStart` in the P case, `stackVerts[d]` in the R/late cases, the chain vertex in
  the S case).
* `vert`, `vert_free`, `q_free`, `q_root`, `v_root`: `UnwrapAt.fresh`/`root`, `ItemFree`,
  `BoundaryOk.q_root/q_free/v_root/v_free`, `FinishTailOk`'s vertex item.
* `p_entry`: the P-check of a type-1 edge/back edge merges into an entry `(curV, lowval)` of `base`
  that is a single root item on the side `stackDir[lowval]`, attached only at `curV`,
  `stackVerts[lowval]` (`UnwrapAt`, `MergeOk`, `FinishTopOk.side/mid` of the P close).
* `bottom` (the author's anchor): for a returning tree edge `sub` ends, bottom-up, with the vertex
  entry of the chain bottom `y` and the `(y, lowval)` piece (`y`'s lowval back edge, or the type-1
  close of `y`'s first child, P-merged with `y`'s other lowval edges), and no loop touches them
  until the fold at the ear's top (`CloseVertOk.merge₁/merge₂/retarget.old`, `MergeOk.share` of the
  fold: both hang at `y`).
* `loops`: the shape after loops 1–2 (`feS₂`): `[c, py, vy]` for a type-1 edge (at most one entry
  of the child's type-1 edges by the P merges, at most one entry reaching `curV` by the
  first-occurrence merges, both present: the lowval edge and the tree edge), `c :: mid ++ [py, vy]`
  for type 2, `lowval ≤ c.topDepth ≤ d` (`CloseVertOk`, `ear_finishP_vert`, `ear_tail_tree`).
* `boundary`: `BoundaryOk.gone/gone₂` (the popped block is separated from the rest by `curV`).

Statements that were FALSE for the real walk (`checks/EarCheck.lean`, 3000 random multigraphs,
every `finishEdge`) and were dropped or weakened: `base_top` (`∀ t ∈ base, topDepth ≤ d`: a vertex
entry `V y` of a finished deeper vertex stays buried under a lower piece), `loop1_bot`
(`topDepth > d ⇒ vStart = o.dest` in loop 1's range: an S-merged entry keeps a deeper chain
vertex as `vStart`), `touch_top` for every entry with `topDepth ≤ d + 1` (buried `V` entries hold
block `Q` items not at `stackVerts[topDepth]`), `vert`/`vert_free` with `vStart = curV` (after a
late merge the entry holding `V curV` has a chain vertex as `vStart`), `loop1_side` without
`lowval < d`, and the candidates "loop 1's range has one entry / entries of depth `≤ d+1` only /
`firstIdx`-ordered". Still missing: the `firstIdx` order behind loop 2's `MergeOk`.
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

/-- The bottom two entries of an ear with lowval `l` and chain bottom `y`: the vertex entry `vy` of
`y` (holding only `vertItem y`, hence only `y`'s boundary blocks) under the `(y, l)` piece `py`, a
single root item on the side `stackDir[l]` attached at `y` and `stackVerts[l]`. -/
structure EarBottom (d l : Nat) (s : WalkState) (py vy : TEntry) : Prop where
  vy_spans : ∃ dir, vy.spans = setSides dir [vertItem vy.vStart] []
  vy_top : d < vy.topDepth
  vy_bd : ∀ v, s.g.Touches (vy.edges s.g s.items) v → v = vy.vStart ∨ s.g.Interior (vy.edges s.g s.items) v
  py_bot : py.vStart = vy.vStart
  py_top : py.topDepth = l
  py_item : ∃ i, py.spans = setSides s.stackDir[l]! [i] [] ∧ ∀ p, ¬ Items.IsParent s.items p i
  py_touch : s.g.Touches (py.edges s.g s.items) vy.vStart ∧ s.g.Touches (py.edges s.g s.items) s.stackVerts[l]!

/-- Ear content of the tstack `sub ++ base` when `finishEdge curV d o _ hasVert` runs (`o` an
out-edge of `curV` at depth `d`, `stackVerts[d] = curV`, and `stackVerts[d+1] = o.dest` for a tree
edge). -/
structure EarFinish (curV d : Nat) (o : DfsOut) (hasVert : Bool) (sub base : List TEntry)
    (s : WalkState) : Prop where
  tstack : s.tstack = sub ++ base
  /-- A back edge creates no entries before `finishEdge`. -/
  back_nil : o.cls.isTree = false → sub = []
  /-- The enclosing entries were not started at the child, the child's not at `curV`. -/
  base_bot : o.cls.isTree = true → ∀ t ∈ base, t.vStart ≠ o.dest
  sub_bot : ∀ t ∈ sub, t.vStart ≠ curV
  /-- The open path has distinct vertices. -/
  path : ∀ k k', k < k' → k' ≤ d → s.stackVerts[k]! ≠ s.stackVerts[k']!
  /-- Open entries are pairwise edge-disjoint and span-disjoint. -/
  disj : s.tstack.Pairwise fun t t' => ∀ e, e < s.g.ne → t.edges s.g s.items e → ¬ t'.edges s.g s.items e
  span_disj : s.tstack.Pairwise fun t t' => ∀ i ∈ t.spans.1 ++ t.spans.2, i ∉ t'.spans.1 ++ t'.spans.2
  /-- Span ownership: the child's entries hold only sub-ear edges, the enclosing ones none, and
  every subtree edge is held by some child entry. -/
  sub_edges : ∀ t ∈ sub, ∀ e, e < s.g.ne → t.edges s.g s.items e → subEdges o e
  base_disj : ∀ t ∈ base, ∀ e, e < s.g.ne → t.edges s.g s.items e → ¬ subEdges o e
  sub_cover : ∀ e, e < s.g.ne → e ≠ o.e → subEdges o e → ∃ t ∈ sub, t.edges s.g s.items e
  /-- Loop 1's range lies on the side `stackDir[d]` (the ear's side), and its entries returning to
  `curV` or the child touch that vertex. -/
  loop1_side : o.cls.lowval d < d → ∀ hi, Loop1Range d sub hi → ∀ t ∈ hi, getSide t.spans (!s.stackDir[d]!) = []
  loop1_touch : o.cls.lowval d < d → ∀ hi, Loop1Range d sub hi → ∀ t ∈ hi,
    (∃ e, e < s.g.ne ∧ t.edges s.g s.items e) → t.topDepth ≤ d + 1 →
    s.g.Touches (t.edges s.g s.items) s.stackVerts[t.topDepth]!
  /-- A non-empty entry touches its bottom. -/
  touch_bot : ∀ t ∈ s.tstack, (∃ e, e < s.g.ne ∧ t.edges s.g s.items e) →
    s.g.Touches (t.edges s.g s.items) t.vStart
  /-- The vertex item of `curV` is on the stack iff `hasVert` (in an enclosing entry returning to
  `curV` or above), and never elsewhere; the edge's `Q` item is on no entry; both are roots. -/
  vert : hasVert = true → ∃ t ∈ base, t.topDepth ≤ d ∧ vertItem curV ∈ t.spans.1 ++ t.spans.2
  vert_free : ∀ t ∈ s.tstack, vertItem curV ∈ t.spans.1 ++ t.spans.2 → hasVert = true ∧ t ∈ base
  q_free : ∀ t ∈ s.tstack, edgeItem s.g o.e ∉ t.spans.1 ++ t.spans.2
  q_root : ∀ p, ¬ Items.IsParent s.items p (edgeItem s.g o.e)
  v_root : ∀ p, ¬ Items.IsParent s.items p (vertItem curV)
  /-- The target of a type-1 P merge: an enclosing `(curV, lowval)` entry is a single root item on
  the side `stackDir[lowval]`, attached only at `curV` and `stackVerts[lowval]`. -/
  p_entry : o.cls.lowval d < d → o.cls.isType1 = true → ∀ t ∈ base, t.vStart = curV →
    t.topDepth = o.cls.lowval d →
    (∃ i, t.spans = setSides s.stackDir[t.topDepth]! [i] [] ∧ ∀ p, ¬ Items.IsParent s.items p i) ∧
    ∀ v, s.g.Touches (t.edges s.g s.items) v →
      v = curV ∨ v = s.stackVerts[t.topDepth]! ∨ s.g.Interior (t.edges s.g s.items) v
  /-- The ear's bottom two entries (`EarBottom`) end `sub`. -/
  bottom : o.cls.isTree = true → o.cls.lowval d < d →
    ∃ mid py vy, sub = mid ++ [py, vy] ∧ EarBottom d (o.cls.lowval d) s py vy
  /-- After loops 1–2 the child's entries are `c :: mid ++ [py, vy]` with the bottom two untouched
  and `lowval ≤ c.topDepth ≤ d`; for a type-1 edge `mid = []`, `c` reaches `curV` and `o.dest`, and
  `c, py, vy` hold all sub-ear edges. -/
  loops : o.cls.isTree = true → o.cls.lowval d < d →
    ∃ c mid py vy, (feS₂ d o s).tstack = c :: mid ++ [py, vy] ++ base ∧
      (∃ mid₀, sub = mid₀ ++ [py, vy]) ∧
      o.cls.lowval d ≤ c.topDepth ∧ c.topDepth ≤ d ∧
      (o.cls.isType1 = true → mid = [] ∧
        (feS₂ d o s).g.Touches (c.edges (feS₂ d o s).g (feS₂ d o s).items) curV ∧
        (feS₂ d o s).g.Touches (c.edges (feS₂ d o s).g (feS₂ d o s).items) o.dest ∧
        ∀ e, e < s.g.ne → subEdges o e →
          c.edges (feS₂ d o s).g (feS₂ d o s).items e ∨ py.edges (feS₂ d o s).g (feS₂ d o s).items e ∨
          vy.edges (feS₂ d o s).g (feS₂ d o s).items e)
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
    s.tstack.length = origTstack := by
  obtain ⟨sub, base, hl, he⟩ := h
  have hs := he.back_nil hb
  subst hs
  rw [he.tstack]; simpa using hl

end WalkState
end Spqr
