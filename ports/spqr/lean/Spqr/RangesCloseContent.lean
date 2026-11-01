import Spqr.RangesCloseSites

/-! # Block contents at the close sites (`CloseContent`)

The item-level contents of the stack entries at the sites of `closeCtx_bd_vert`/`bd_node`/
`p_site`/`v_site`/`l1_site` (PROOF.md §4.6), stated at the deterministic site states of one
`finishEdge curV d o origTstack hasVert` call from `s`. `CloseContent` is what the walk induction
(`WalkBackbone.lean`) supplies at every `finishEdge` call; the five `closeCtx_*` theorems
(`RangesCloseTree.lean`) derive the full `PSite`/`VSite` records from it and the exports of
`CloseBase` (the graph-level clauses `shape`/`ne`/`att`/`pend` come from the range invariant, the
ear contracts and the open path). The executable checker (`checks/WalkInvCheck/Ranges.lean`,
`checkContent`) evaluates every field separately (`content_*`). -/

namespace Spqr
open WalkM
namespace WalkState

/-- The top two entries at a P site (`PSite` without `shape`/`ne`/`pend`): `cur` is `(curV, lv)`;
each entry is one item on the side `stackDir[lv]` with terminals `{stackVerts[lv], curV}`, touches
both and is attached nowhere else; the items, and the children of a P item among them, are S/P/R or
leaf Q; each is a root listed once on the stack. -/
structure PContent (curV lv : Nat) (r : WalkState) : Prop where
  stack : ∃ cur rest, r.tstack = cur :: rest ∧ cur.vStart = curV ∧ cur.topDepth = lv
  single : ∀ t ∈ r.tstack.take 2, ∃ j, t.spans = setSides r.stackDir[lv]! [j] [] ∧
    ∃ a b, Items.vs r.items j = (some a, some b) ∧ Items.PairEq (a, b) (r.stackVerts[lv]!, curV)
  att : ∀ t ∈ r.tstack.take 2, ∀ w, r.g.Touches (t.edges r.g r.items) w →
    w = curV ∨ w = r.stackVerts[lv]! ∨ r.g.Interior (t.edges r.g r.items) w
  touch : ∀ t ∈ r.tstack.take 2, r.g.Touches (t.edges r.g r.items) curV ∧
    r.g.Touches (t.edges r.g r.items) r.stackVerts[lv]!
  kinds : ∀ t ∈ r.tstack.take 2, ∀ j ∈ t.spans.1 ++ t.spans.2, ∀ j',
    j' = j ∨ (Items.type r.items j = .P ∧ Items.IsParent r.items j j') →
    Items.type r.items j' ∈ [NodeType.S, .P, .R, .Q] ∧
      (Items.type r.items j' = .Q → Items.ch r.items j' = [])
  once : ∀ t ∈ r.tstack.take 2, ∀ j ∈ t.spans.1 ++ t.spans.2,
    (∀ p, ¬ Items.IsParent r.items p j) ∧ spansCount r.tstack j = 1

/-- The top entry `t` about to be closed into the free node `x` (`VSite` without
`shape`/`ne`/`att`/`pend`): `t` bottoms at `curV` with everything on its `stackDir[topDepth]` side;
the side's children are S/P/R/Q/V, non-V ones with two distinct terminals; both terminals are
touched; a vertex item is a child iff all its edges lie in the side and in no single child; and the
S path / R / P shape clauses hold for the side's virtual edges. -/
structure VContent (curV : Nat) (x : ItemId) (t : TEntry) (r : WalkState) : Prop where
  stack : ∃ rest, r.tstack = t :: rest
  vstart : t.vStart = curV
  side : getSide t.spans (!r.stackDir[t.topDepth]!) = []
  free : ItemFree r x
  kinds : ∀ c ∈ vKids t r, Items.type r.items c ∈ [NodeType.S, .P, .R, .Q, .V]
  two : ∀ c ∈ vKids t r, Items.type r.items c ≠ .V →
    ∃ a b, Items.vs r.items c = (some a, some b) ∧ a ≠ b
  touch : r.g.Touches (t.edges r.g r.items) curV ∧
    r.g.Touches (t.edges r.g r.items) r.stackVerts[t.topDepth]!
  inner : ∀ w, w < r.g.nv → (vertItem w ∈ vKids t r ↔
    r.g.Touches (t.edges r.g r.items) w ∧ r.g.Interior (t.edges r.g r.items) w ∧
    ∀ c ∈ vKids t r, ¬ r.g.Interior (Items.EdgeBelow r.g r.items c) w)
  s_order : Items.type r.items x = .S → ∃ xs,
    ((vKids t r).filter fun c => Items.type r.items c = .V) = xs.map vertItem ∧ 1 ≤ xs.length ∧
    vVirt t r = List.zip ((vTerms curV t r).1 :: xs) (xs ++ [(vTerms curV t r).2])
  r_shape : Items.type r.items x = .R →
    2 ≤ ((vKids t r).filter fun c => Items.type r.items c = .V).length ∧ 5 ≤ (vVirt t r).length ∧
    ((vVirt t r).map fun q => min q.1 q.2 + r.g.nv * max q.1 q.2).Nodup ∧
    ∀ q ∈ vVirt t r, ¬ Items.PairEq q (vTerms curV t r)
  p_shape : Items.type r.items x = .P → 2 ≤ (vVirt t r).length ∧
    (∀ c ∈ vKids t r, Items.type r.items c ≠ .V) ∧ ∀ q ∈ vVirt t r, Items.PairEq q (vTerms curV t r)

/-- The merged pre-close state of the `k`-th loop-1 iteration. -/
def l1Pre (d : Nat) (o : DfsOut) (s : WalkState) (k : Nat) : WalkState :=
  after mergeTstackTops (l1S₂ d s.stackDir[d]! (l1Iter d o s k))

/-- The node `maybeUnwrapNxt` returns in the `k`-th loop-1 iteration. -/
def l1Node (d : Nat) (o : DfsOut) (s : WalkState) (k : Nat) : ItemId :=
  result (maybeUnwrapNxt (l1Ty d s.stackDir[d]! (l1Iter d o s k))) (l1S₁ d s.stackDir[d]! (l1Iter d o s k))

/-- The block contents of one `finishEdge curV d o origTstack hasVert` call from `s`:
* `bd_vert`: the block entry of a boundary tree edge (`(o.dest, d + 1)`, the top entry when
  `lowval = d + 1`, else the one under the completed block) holds exactly `vertItem o.dest`;
* `bd_node`: at a completed block (`d ≤ lowval ≠ d + 1`) the top entry's passive side is one closed
  node, a leaf if `Q`, with terminals `{curV, o.dest}`;
* `p_site`: `PContent` at the state `finishP` runs from whenever `condP` holds there;
* `v_site`: `VContent` at `cvS₅` of a type-1 vertex close, for the unwrapped node;
* `l1_site`: at every loop-1 iteration whose first `k` conditions held, `VContent` (with
  `curV := t.vStart`) at the merged pre-close state for the unwrapped node, the terminals distinct,
  each with an edge outside the entry. -/
structure CloseContent (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool)
    (s : WalkState) : Prop where
  bd_vert : o.cls.isTree = true → d ≤ o.cls.lowval d →
    ∀ t ∈ (if o.cls.lowval d = d + 1 then s.tstack.head? else s.tstack.tail.head?),
      t.spans.2 = [vertItem o.dest]
  bd_node : o.cls.isTree = true → d ≤ o.cls.lowval d → o.cls.lowval d ≠ d + 1 →
    ∀ b ∈ s.tstack.head?, ∃ c, b.spans.1 = [c] ∧ Items.type s.items c ∉ [NodeType.F, .V] ∧
      (Items.type s.items c = .Q → Items.ch s.items c = []) ∧
      ∃ a b', Items.vs s.items c = (some a, some b') ∧ Items.PairEq (a, b') (curV, o.dest)
  p_site : o.cls.lowval d < d → o.cls.isType1 = true →
    result (condP curV (o.cls.lowval d) true) (feRest curV d o origTstack hasVert s) = true →
    PContent curV (o.cls.lowval d) (feRest curV d o origTstack hasVert s)
  v_site : o.cls.isTree = true → o.cls.lowval d < d → hasVert = true → o.cls.isType1 = true →
    ∃ t, VContent curV ((maybeUnwrapNxt (if feSingle d o s then NodeType.S else .R)).run (feS₂ d o s)).1 t
      (cvS₅ curV s.stackDir[d]! true origTstack (feSingle d o s) (feS₂ d o s))
  l1_site : o.cls.isTree = true → o.cls.lowval d < d → ∀ k,
    (∀ j, j ≤ k → result (loop1Cond d) (l1Iter d o s j) = true) →
    ∃ t, t.vStart ≠ (l1Pre d o s k).stackVerts[t.topDepth]! ∧
      (∀ w, w = t.vStart ∨ w = (l1Pre d o s k).stackVerts[t.topDepth]! →
        ∃ e, e < (l1Pre d o s k).g.ne ∧ (l1Pre d o s k).g.Inc e w ∧
          ¬ t.edges (l1Pre d o s k).g (l1Pre d o s k).items e) ∧
      VContent t.vStart (l1Node d o s k) t (l1Pre d o s k)

variable {σ : List Nat} {n D : Nat} {s : WalkState}

theorem after_mergeTstackTops_eq {a b : TEntry} {rest : List TEntry} (hts : s.tstack = a :: b :: rest) :
    after mergeTstackTops s = { s with tstack := TEntry.mergeInto a b :: rest } := by
  show (mergeTstackTops.run s).2 = _; rw [mergeTstackTops_run_eq s a b rest hts]


end WalkState
end Spqr
