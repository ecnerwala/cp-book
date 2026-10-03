import Spqr.EarInv

/-!
# `EarCtx`: the between-edges invariant of `walkOuts v d` (PROOF.md §4.2b)

Statement only — no proof attempt. Checked literally in `checks/EarCheck.lean` (`ctxCheck`, every
clause, at the start of `walkOuts` and after every `walkOut` return): 0 violations on seeds 0..3000
with `ternarize = false` and `true`. It is meant to be the induction hypothesis of the `EarTree`
induction (`WalkInv.walkTree_ear`): the next site's `EarAt v d o origTstack hasVert` is to be
derived from it (plus the child's own `EarCtx` at the end of its outs for the loop-1/2/close
fields), and `walkOut` re-establishes it.

Parameters: `done` = the outs of `v` already finished, each tagged with whether the vertex entry
had been pushed when its `finishEdge` ran (`hasVert ∨ (lowval < d ∧ isType1)` before the walk);
`rest` = the outs still to come; `hasVert` = the vertex entry of `v` has been pushed (by
`walkOutPre` or by `finishEdge` itself); `base` = the tstack at entry to `v` (compared as values),
`bE` = its entries' edge sets then; `sv` = `stackVerts[0..d]` at entry.

Shape of the stack: `above ++ vt :: below ++ base` once `hasVert`, `top ++ base` before, where
* `vt` is the unique entry holding `vertItem v`, `vt.topDepth ≤ d` — it is `V v` itself until a
  type-2 out closes with the vertex entry, after which it is the merged entry with a chain vertex
  as `vStart` (seed 1, `v = 4`, `d = 3`: `(1,1,[29,7,15,23,30,5])` holds `V 4 = 5`); it may also
  sit *above* the chain child's pieces (seed 1, `v = 5`, `d = 5`:
  `[V 5, (6,5), (6,3), (6,1), V 6] ++ base`), which is why the per-out entries are `above`, not
  `top`;
* `above` are the entries of the outs finished after the vertex push: each `(v, l)` with `l` the
  lowval of such an out, non-empty, touching `v` and `stackVerts[l]`, on side `stackDir[l]`,
  attached only at `v` and at `stackVerts[k]`, `l ≤ k ≤ d`; when every such out at `l` is type 1
  it is a single root item attached only at `v` and `stackVerts[l]` (type-2 outs leave a
  multi-item entry and do not P-merge: seed 1, `v = 1`, `d = 4` has `(1,1,[29,7,15,23])` above
  `(1,1,[16])`); `topDepth` is non-increasing and `firstIdx` strictly decreasing top-down; every
  lowval of those outs has an entry in `above` (or is `vt`'s depth), and `above` holds exactly
  their `subEdges` (up to `vt`);
* `below` (and all of `top` before the push) has no entry started at `v`: it is the open ear of
  the chain child.

Which `EarFinish` field each clause is meant to re-establish at the next `finishEdge` site
(`sub` = what the next out's child leaves, `base'` = the whole stack now, with `V v` pushed on top
first when the next out is the first type-1 returning one):
* `base_eq`/`base_edges`/`sv`/`sv_d`/`path` → `path`, `sv_d`, `sv_child`/`path_child` (with the
  child's freshness), `base_bot`; they are the `Keep` frame of `walkTree_frame` restricted to `v`;
* `v_root`/`vert_free`/`vert_le` → `vert`, `vert_free`, `v_root` (and the `UnwrapAt`-style
  freshness of the vertex close in `ear_closeVert`);
* `above_*`/`open_bot`/`base_bot` → `p_entry` (the `(v, lowval)` entries of `base'` are exactly
  those of `above` at that lowval; the single-root-item/attachment part is the `allT1` case),
  `touch_bot` for the enclosing entries, and `MergeOk.share/bottom` of the P merge in
  `ear_finishP_vert`;
* `above_sub_edges`/`above_cover`/`vert_edges` → `base_disj` (no enclosing entry holds an edge of
  a remaining subtree), `dest_edges`, and `sub_edges`/`sub_cover` of the *parent's* site when `v`
  is finished (the entries `v` leaves hold exactly `v`'s subtree edges);
* `disj`/`span_disj` → `disj`/`span_disj`;
* `q_fresh`/`v_fresh` → `q_free`/`q_root`, the child sites' `v_root`/`vert_free`, `path_child`,
  `base_bot`.
Not covered — they need the child's end-of-outs `EarCtx` plus the chain anchor (`bottom`: the
chain bottom's `V y`/`(y, lowval)` pair, never merged until the ear's top): `loop1`,
`loop1_side`/`loop1_touch`, `late`/`late_fo`, `loops`, `close`, `bottom`, `boundary`/`bd_*`,
`lower`, `dir_d`.
-/

namespace Spqr
namespace WalkState

/-- An entry of `above`: `(v, l)` with `l < d`, non-empty, touching `v` and `stackVerts[l]`, on
side `stackDir[l]`, attached only at `v` and at path vertices `stackVerts[k]`, `l ≤ k ≤ d`. -/
structure CtxEntry (v d : Nat) (s : WalkState) (t : TEntry) : Prop where
  vStart : t.vStart = v
  depth : t.topDepth < d
  nonempty : ∃ e, e < s.g.ne ∧ t.edges s.g s.items e
  touch_bot : s.g.Touches (t.edges s.g s.items) v
  touch_top : s.g.Touches (t.edges s.g s.items) s.stackVerts[t.topDepth]!
  side : getSide t.spans (!s.stackDir[t.topDepth]!) = []
  att : ∀ x, s.g.Touches (t.edges s.g s.items) x → x ≠ v → ¬ s.g.Interior (t.edges s.g s.items) x →
    ∃ k, t.topDepth ≤ k ∧ k ≤ d ∧ x = s.stackVerts[k]!

/-- The P-merge target shape (`EarFinish.p_entry`): a single root item attached only at `v` and
`stackVerts[topDepth]`. -/
structure CtxSingle (v : Nat) (s : WalkState) (t : TEntry) : Prop where
  item : ∃ i, t.spans = setSides s.stackDir[t.topDepth]! [i] [] ∧ ∀ p, ¬ Items.IsParent s.items p i
  att : ∀ x, s.g.Touches (t.edges s.g s.items) x →
    x = v ∨ x = s.stackVerts[t.topDepth]! ∨ s.g.Interior (t.edges s.g s.items) x

/-- The outs finished after the vertex push. -/
def afterVert (done : List (DfsOut × Bool)) : List DfsOut := (done.filter (·.2)).map (·.1)

/-- Every out finished after the vertex push at lowval `l` is type 1. -/
def allType1 (d l : Nat) (done : List (DfsOut × Bool)) : Prop :=
  ∀ o ∈ afterVert done, o.cls.lowval d = l → o.cls.isType1 = true

structure EarCtx (v d : Nat) (done : List (DfsOut × Bool)) (rest : List DfsOut) (hasVert : Bool)
    (base : List TEntry) (bE : List (Nat → Prop)) (sv : List Nat) (s : WalkState) : Prop where
  /-- The stack: `above ++ vt :: below ++ base` once `hasVert`, `below ++ base` before. -/
  split : ∃ above below,
    (if hasVert then
      ∃ vt, s.tstack = above ++ vt :: below ++ base ∧ vt.topDepth ≤ d ∧
        vertItem v ∈ vt.spans.1 ++ vt.spans.2 ∧
        (∀ t ∈ above, vertItem v ∉ t.spans.1 ++ t.spans.2) ∧
        (∀ o ∈ afterVert done, (∃ t ∈ above, t.topDepth = o.cls.lowval d) ∨
          vt.topDepth = o.cls.lowval d) ∧
        (∀ o ∈ afterVert done, ∀ e, e < s.g.ne → subEdges o e →
          (∃ t ∈ above, t.edges s.g s.items e) ∨ vt.edges s.g s.items e)
    else above = [] ∧ s.tstack = below ++ base) ∧
    (∀ t ∈ above, CtxEntry v d s t) ∧
    (∀ t ∈ above, allType1 d t.topDepth done → CtxSingle v s t) ∧
    above.Pairwise (fun t t' => t'.topDepth ≤ t.topDepth ∧ t'.firstIdx < t.firstIdx) ∧
    (∀ t ∈ above, ∃ o ∈ afterVert done, t.topDepth = o.cls.lowval d) ∧
    (∀ t ∈ above, ∀ e, e < s.g.ne → t.edges s.g s.items e → ∃ o ∈ afterVert done, subEdges o e) ∧
    (∀ t ∈ below, t.vStart ≠ v)
  /-- The enclosing entries are untouched (values and edge sets), none was started at `v`. -/
  base_edges : bE.length = base.length ∧ ∀ k, k < base.length → ∀ e, e < s.g.ne →
    (base[k]!.edges s.g s.items e ↔ bE[k]! e)
  base_bot : ∀ t ∈ base, t.vStart ≠ v
  /-- The path. -/
  sv : ∀ k, k ≤ d → s.stackVerts[k]! = sv[k]!
  sv_d : s.stackVerts[d]! = v
  path : ∀ k k', k < k' → k' ≤ d → s.stackVerts[k]! ≠ s.stackVerts[k']!
  /-- The vertex item of `v`: a root, on the stack only in `vt`. -/
  v_root : ∀ p, ¬ Items.IsParent s.items p (vertItem v)
  vert_free : ∀ t ∈ s.tstack, vertItem v ∈ t.spans.1 ++ t.spans.2 → hasVert = true
  afterVert_ret : ∀ o ∈ afterVert done, o.cls.lowval d < d
  /-- The boundary blocks (`d ≤ lowval`) are exactly what hangs under `vertItem v`. -/
  vert_edges : ∀ e, e < s.g.ne →
    (Items.EdgeBelow s.g s.items (vertItem v) e ↔ ∃ o ∈ done, d ≤ o.1.cls.lowval d ∧ subEdges o.1 e)
  /-- Open entries are pairwise edge- and span-disjoint. -/
  disj : s.tstack.Pairwise fun t t' => ∀ e, e < s.g.ne → t.edges s.g s.items e → ¬ t'.edges s.g s.items e
  span_disj : s.tstack.Pairwise fun t t' => ∀ i ∈ t.spans.1 ++ t.spans.2, i ∉ t'.spans.1 ++ t'.spans.2
  /-- Freshness of the remaining outs: their `Q` items and subtree `V` items are childless roots on
  no span, and their subtree vertices are off the path. -/
  q_fresh : ∀ o ∈ rest, ∀ e, subEdges o e →
    (∀ p, ¬ Items.IsParent s.items p (edgeItem s.g e)) ∧ Items.ch s.items (edgeItem s.g e) = [] ∧
    ∀ t ∈ s.tstack, edgeItem s.g e ∉ t.spans.1 ++ t.spans.2
  v_fresh : ∀ o ∈ rest, ∀ e cls child, o = .tree e cls child → ∀ y ∈ child.verts,
    (∀ p, ¬ Items.IsParent s.items p (vertItem y)) ∧ Items.ch s.items (vertItem y) = [] ∧
    (∀ t ∈ s.tstack, vertItem y ∉ t.spans.1 ++ t.spans.2) ∧
    ∀ k, k ≤ d → s.stackVerts[k]! ≠ y

end WalkState
end Spqr
