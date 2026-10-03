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
`firstIdx`-ordered".

Ear session 5 (`late`, `close`, `sv_d`, `sv_child`, `path_child`, `dir_d`): loop 2 is read at `feS₁`
(`EarLate`: `MergeOk` at every split the loop reaches, i.e. while every merged entry has
`firstIdx > firstOccurrence[d]`, plus pairwise edge-disjointness there); the state after loops 1–2
is described by `EarClose` (`c :: mid ++ [py, vy]` over the untouched `base`, span ownership of the
sub-ear edges, every edge at the chain bottom `y` is a sub-ear edge, `FoldSpec`: `MergeOk` at every
step of folding the child's entries top-down — loop 3 and the two merges of the vertex close, or the
`!hasVert` merge — and for type 1 the touch set `{curV, stackVerts[lowval]} ∪ interior`); all fields
0 violations on 3000 random multigraphs (`checks/EarCheck.lean`, `closeCheck`). `vy.vStart` is
*not* `o.dest` in general (type-2 chain), so `RetargetOk.old`/`MergeOk.bottom` go through
`Interior` (`y_edges`), never through the path.
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

/-- Loop 1 read on the original stack: after consuming the entries `done` (top first) of its range,
the current piece has bottom `l1Bot o done` (the bottom of the last consumed entry, initially the
child) and edge set `l1Edges o s done` (the tree edge plus the consumed entries' edges). -/
def l1Bot (o : DfsOut) (done : List TEntry) : Nat := ((done.getLast?).map TEntry.vStart).getD o.dest
def l1Edges (o : DfsOut) (s : WalkState) (done : List TEntry) (e : Nat) : Prop :=
  e = o.e ∨ ∃ t ∈ done, t.edges s.g s.items e

/-- Merging the loop-1 piece `(v, E)` (bottom `v`, edges `E`, top `d`) with the entry `t`:
`MergeOk.share`/`bottom` at depth `d + 1`. -/
structure L1Merge (d : Nat) (s : WalkState) (v : Nat) (E : Nat → Prop) (t : TEntry) : Prop where
  share : (∃ e, e < s.g.ne ∧ t.edges s.g s.items e) →
    ∃ x, s.g.Touches E x ∧ s.g.Touches (t.edges s.g s.items) x
  bottom : v = t.vStart ∨ (∃ k, d ≤ k ∧ k ≤ d + 1 ∧ v = s.stackVerts[k]!) ∨
    s.g.Interior (fun e => E e ∨ t.edges s.g s.items e) v

/-- Closing the loop-1 piece `(v, E)` with the entry `t`: the merge, and `FinishTopOk.mid` at the
only intermediate depth `d + 1`. -/
structure L1Close (d : Nat) (s : WalkState) (v : Nat) (E : Nat → Prop) (t : TEntry) : Prop
    extends L1Merge d s v E t where
  mid : s.stackVerts[d + 1]! = t.vStart ∨
    s.g.Interior (fun e => E e ∨ t.edges s.g s.items e) s.stackVerts[d + 1]! ∨
    ¬ s.g.Touches (fun e => E e ∨ t.edges s.g s.items e) s.stackVerts[d + 1]!

/-- `t` unwraps as a `ty` node (`maybeUnwrapNxt ty`): if the head of a side of `t` is a `ty` item,
`t` is that single root item, whose children are on no entry. -/
def L1Unwrap (s : WalkState) (ty : NodeType) (t : TEntry) : Prop :=
  ∀ dir h, (getSide t.spans dir).head! = h → Items.type s.items h = ty →
    getSide t.spans dir = [h] ∧ (∀ p, ¬ Items.IsParent s.items p h) ∧
    ∀ c ∈ Items.ch s.items h, ∀ u ∈ s.tstack, c ∉ u.spans.1 ++ u.spans.2

/-- The splits `done ++ rest` of loop 1's range reached at iteration boundaries: an entry at depth
`d` is closed alone, a deeper one is S-merged and closed together with the entry under it. -/
inductive L1Reach (d : Nat) : List TEntry → List TEntry → List TEntry → Prop
  | nil (hi : List TEntry) : L1Reach d hi [] hi
  | close {hi done : List TEntry} {t : TEntry} {rest : List TEntry} :
    L1Reach d hi done (t :: rest) → t.topDepth = d → L1Reach d hi (done ++ [t]) rest
  | series {hi done : List TEntry} {t t' : TEntry} {rest : List TEntry} :
    L1Reach d hi done (t :: t' :: rest) → d < t.topDepth → L1Reach d hi (done ++ [t, t']) rest

/-- What loop 1 needs of its range `hi` (top first), read in the state before `finishEdge`: at every
reached split, the next entry `t` closes (depth `d`; as a `P` node when it shares the piece's
bottom), or is S-merged and the entry under it closes (as an `S` node). -/
def Loop1Spec (d : Nat) (o : DfsOut) (s : WalkState) (hi : List TEntry) : Prop :=
  ∀ done t rest, L1Reach d hi done (t :: rest) →
    (t.topDepth = d → L1Close d s (l1Bot o done) (l1Edges o s done) t ∧
      (t.vStart = l1Bot o done → L1Unwrap s .P t)) ∧
    (d < t.topDepth → L1Merge d s (l1Bot o done) (l1Edges o s done) t ∧
      ∃ t' rest', rest = t' :: rest' ∧ L1Close d s t.vStart (l1Edges o s (done ++ [t])) t' ∧
        L1Unwrap s .S t')

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

/-- The piece obtained by merging the entries `done` (top first) one by one into `c₀`
(`mergeTstackTops` folded: `cur := mergeInto cur t`). -/
def l2Cur (c₀ : TEntry) (done : List TEntry) : TEntry := done.foldl TEntry.mergeInto c₀

/-- Merging the entries `R` (top first) one by one into the piece `c₀`: `MergeOk` at every split. -/
def FoldSpec (D : Nat) (st : WalkState) (c₀ : TEntry) (R : List TEntry) : Prop :=
  ∀ done t rest, R = done ++ t :: rest → MergeOk D st (l2Cur c₀ done) t

/-- Loop 2 (`mergeLate`), read at `feS₁` (after loop 1): the stack is `c₀ :: R`, pairwise
edge-disjoint, and merging `R` into `c₀` is `MergeOk` (at depth bound `d + 1`) at every split the
loop reaches, i.e. while `c₀` and every merged entry were opened after `firstOccurrence[d]`
(`loop2Cond`: `cur.firstIdx > fo`, and `mergeInto` keeps `nxt.firstIdx`). -/
structure EarLate (d : Nat) (st : WalkState) (c₀ : TEntry) (R : List TEntry) : Prop where
  tstack : st.tstack = c₀ :: R
  disj : st.tstack.Pairwise fun t t' => ∀ e, e < st.g.ne → t.edges st.g st.items e → ¬ t'.edges st.g st.items e
  merge : ∀ done t rest, R = done ++ t :: rest → st.firstOccurrence[d]! < c₀.firstIdx →
    (∀ u ∈ done, st.firstOccurrence[d]! < u.firstIdx) → MergeOk (d + 1) st (l2Cur c₀ done) t

/-- The stack after loops 1–2 (`st = feS₂ d o s`, `l = lowval`): the child's entries are
`c :: mid ++ [py, vy]` over the untouched `base`; they hold exactly the sub-ear edges, every edge at
the chain bottom `y = vy.vStart` is one of them, and folding them top-down (loop 3, then the two
merges of the vertex close, or the single `!hasVert` merge) is `MergeOk` at every step. For a type-1
edge the three entries touch only `curV`, `stackVerts[l]` and interior vertices
(`FinishTopOk.mid`/`RetargetOk.old` of the vertex close and the P-check). When the vertex item of
`curV` is not on the stack yet, its edges are on no entry and it touches `curV` if non-empty
(`ear_tail_tree`). -/
structure EarClose (curV d l : Nat) (o : DfsOut) (hasVert : Bool) (base : List TEntry) (s st : WalkState)
    (c : TEntry) (mid : List TEntry) (py vy : TEntry) : Prop where
  tstack : st.tstack = c :: mid ++ [py, vy] ++ base
  g : st.g = s.g
  sv : st.stackVerts = s.stackVerts
  dir_l : st.stackDir[l]! = s.stackDir[l]!
  c_top : l ≤ c.topDepth ∧ c.topDepth ≤ d
  base_edges : ∀ t ∈ base, ∀ e, t.edges s.g st.items e ↔ t.edges s.g s.items e
  base_root : ∀ t ∈ base, ∀ i ∈ t.spans.1 ++ t.spans.2, (∀ p, ¬ Items.IsParent s.items p i) →
    ∀ p, ¬ Items.IsParent st.items p i
  disj : st.tstack.Pairwise fun t t' => ∀ e, e < s.g.ne → t.edges s.g st.items e → ¬ t'.edges s.g st.items e
  span_disj : st.tstack.Pairwise fun t t' => ∀ i ∈ t.spans.1 ++ t.spans.2, i ∉ t'.spans.1 ++ t'.spans.2
  sub_edges : ∀ t ∈ c :: mid ++ [py, vy], ∀ e, e < s.g.ne → t.edges s.g st.items e → subEdges o e
  sub_cover : ∀ e, e < s.g.ne → subEdges o e → ∃ t ∈ c :: mid ++ [py, vy], t.edges s.g st.items e
  y_edges : ∀ e, e < s.g.ne → s.g.Inc e vy.vStart → subEdges o e
  vy_top : d < vy.topDepth
  c_edge : c.edges s.g st.items o.e
  py_item : ∃ i, py.spans = setSides s.stackDir[l]! [i] [] ∧ ∀ p, ¬ Items.IsParent st.items p i
  py_bot : py.vStart = vy.vStart
  py_top : py.topDepth = l
  fold : FoldSpec (d + 1) st c (mid ++ [py, vy])
  type1 : o.cls.isType1 = true → mid = [] ∧
    ∀ v, s.g.Touches (fun e => c.edges s.g st.items e ∨ py.edges s.g st.items e ∨ vy.edges s.g st.items e) v →
      v = curV ∨ v = s.stackVerts[l]! ∨
      s.g.Interior (fun e => c.edges s.g st.items e ∨ py.edges s.g st.items e ∨ vy.edges s.g st.items e) v
  vert_disj : hasVert = false → ∀ t ∈ st.tstack, ∀ e, e < s.g.ne → t.edges s.g st.items e →
    ¬ Items.EdgeBelow s.g st.items (vertItem curV) e
  vert_touch : hasVert = false → (∃ e, e < s.g.ne ∧ Items.EdgeBelow s.g st.items (vertItem curV) e) →
    s.g.Touches (Items.EdgeBelow s.g st.items (vertItem curV)) curV

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
  /-- Loop 1's per-iteration merge/unwrap/close facts (`Loop1Spec`). -/
  loop1 : o.cls.isTree = true → o.cls.lowval d < d → ∀ hi, Loop1Range d sub hi → Loop1Spec d o s hi
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
  /-- Loop 2 read after loop 1 (`EarLate`), and the stack after loops 1–2 (`EarClose`). -/
  late : o.cls.isTree = true → o.cls.lowval d < d → ∃ c₀ R, EarLate d (feS₁ d o s) c₀ R
  close : o.cls.isTree = true → o.cls.lowval d < d →
    ∃ c mid py vy, EarClose curV d (o.cls.lowval d) o hasVert base s (feS₂ d o s) c mid py vy
  /-- The open path: `stackVerts[d] = curV`, `stackVerts[d+1]` is the child (not on the path up to
  `d`), and the ear's side `stackDir[d]` is opposite to `stackDir[lowval]` (`finishSetup`). -/
  sv_d : s.stackVerts[d]! = curV
  sv_child : o.cls.isTree = true → s.stackVerts[d + 1]! = o.dest
  path_child : o.cls.isTree = true → ∀ k, k ≤ d → s.stackVerts[k]! ≠ o.dest
  dir_d : o.cls.lowval d < d → s.stackDir[d]! = !s.stackDir[o.cls.lowval d]!
  /-- A block boundary: the child's block meets the rest only at the articulation vertex `curV`. -/
  boundary : d ≤ o.cls.lowval d → ∀ t ∈ sub, ∀ u ∈ base, ∀ v,
    s.g.Touches (t.edges s.g s.items) v → s.g.Touches (u.edges s.g s.items) v → v = curV
  /-- Every edge at the child is a sub-ear edge. -/
  dest_edges : o.cls.isTree = true → ∀ e, e < s.g.ne → s.g.Inc e o.dest → subEdges o e
  /-- A boundary edge (`d ≤ lowval`) comes before the vertex entry of `curV`; the child's entries are
  the block `(o.dest, d + 1)` (bridge) or `[(o.dest, lowval), (o.dest, d + 1)]` (component: the block
  on the child's vertex entry); an entry under the top one touching `curV` has it as a terminal. -/
  bd_noVert : d ≤ o.cls.lowval d → hasVert = false
  bd_bridge : o.cls.isTree = true → o.cls.lowval d = d + 1 →
    ∃ t, sub = [t] ∧ t.vStart = o.dest ∧ t.topDepth = d + 1
  bd_comp : o.cls.isTree = true → d ≤ o.cls.lowval d → o.cls.lowval d ≠ d + 1 →
    ∃ t₁ t₂, sub = [t₁, t₂] ∧ t₁.vStart = o.dest ∧ t₁.topDepth = o.cls.lowval d ∧
      t₂.vStart = o.dest ∧ t₂.topDepth = d + 1
  bd_term : o.cls.isTree = true → d ≤ o.cls.lowval d → ∀ u ∈ s.tstack.tail,
    s.g.Touches (u.edges s.g s.items) curV → u.vStart = curV ∨ u.topDepth ≤ d
  /-- The popped boundary entries lie on the side the walk fixed for them (`BoundaryOK`): the
  block entry has `spans.1 = []`; at a component edge the `(o.dest, lowval)` entry above it has
  `spans.2 = []`. -/
  bd_side : o.cls.isTree = true → d ≤ o.cls.lowval d →
    if o.cls.lowval d = d + 1 then ∀ t ∈ s.tstack.head?, t.spans.1 = []
    else (∀ b ∈ s.tstack.head?, b.spans.2 = []) ∧ ∀ t ∈ s.tstack.tail.head?, t.spans.1 = []
  /-- After a tree edge the child is attached only as a bottom: an entry of the final stack touching
  `o.dest` has it interior, or as its `vStart` or that of an entry above it (`ear_lower'`); the
  frame facts are included until the preservation induction supplies them. -/
  lower : o.cls.isTree = true →
    (after (finishEdge curV d o base.length hasVert) s).g = s.g ∧
    (after (finishEdge curV d o base.length hasVert) s).stackVerts = s.stackVerts ∧
    ∀ above t below, (after (finishEdge curV d o base.length hasVert) s).tstack = above ++ t :: below →
      s.g.Touches (t.edges s.g (after (finishEdge curV d o base.length hasVert) s).items) o.dest →
      s.g.Interior (t.edges s.g (after (finishEdge curV d o base.length hasVert) s).items) o.dest ∨
      o.dest = t.vStart ∨ ∃ t' ∈ above, o.dest = t'.vStart

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
