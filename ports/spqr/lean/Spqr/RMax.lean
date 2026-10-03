import Spqr.Contract
import Spqr.SepPairExhaust

/-!
# R-maximality: the shape of an R close (PROOF.md §4.3 "R case", §4.5)

`finishEdge` closes an R item from tstack entries `E_1..E_k` (loop 1: equal `topDepth`, distinct
`vStart`s) whose union `U` is connected and 2-attached at `{s, t}` (`s = vStart` of the top entry,
`t = stackVerts[topDepth]`); the merged entries' items are the maximal sub-pieces `P` of `U`, and
the R skeleton is `(P.addParent g U s t).contract g` — `U`'s complement is one more piece (the
parent virtual edge).

`RCloseShape` names every walk-side fact the graph theory consumes; it is a hypothesis here, to
be discharged by the walk invariant (`WalkSpec.lean`/`Frame.lean`), and makes no reference to
`WalkState`.
-/

namespace Spqr

namespace Pieces

variable (P : Pieces)

/-- `K` is laminar with `U` and its sub-pieces: inside one sub-piece, disjoint from `U`, or
containing `U`. -/
def LaminarWith (U K : Nat → Prop) : Prop :=
  (∃ i, i < P.k ∧ ∀ e, K e → P.Mem i e) ∨ (∀ e, K e → ¬U e) ∨ ∀ e, U e → K e

open Classical in
/-- The pieces of the R skeleton: the sub-pieces `P` of `U` plus the complement of `U` as piece
`P.k` with terminals `s`, `t`. -/
noncomputable def addParent (g : Graph) (U : Nat → Prop) (s t : Nat) : Pieces where
  k := P.k + 1
  piece e := if U e then P.piece e else if e < g.ne then some P.k else none
  x i := if i = P.k then s else P.x i
  y i := if i = P.k then t else P.y i

/-- A pair of vertices of `U` on the R skeleton (interior to no sub-piece), other than a
sub-piece's terminal pair or `{s, t}`: the pairs that must not separate the block. -/
def SkelPair (g : Graph) (U : Nat → Prop) (s t a b : Nat) : Prop :=
  g.Touches U a ∧ g.Touches U b ∧ P.Skel g a ∧ P.Skel g b ∧ ¬P.TermPair a b ∧
    ¬((a = s ∧ b = t) ∨ (a = t ∧ b = s))

end Pieces

/-- The facts about an R close of `finishEdge`, for the block `g` with sorted DFS tree `d`:
* `wf`, `sub`: the merged entries' items are well-formed pieces inside `U`;
* `conn`, `attached`, `touch_s`, `touch_t`, `ne`, `proper`: `U` (the union of the closed
  entries) is a connected proper edge set 2-attached at its two terminals;
* `single`: `U` is one separation class of `{s, t}` (the close is an R, not a P);
* `maximal`: each sub-piece is maximal — its complement is one class of its terminals;
* `type1`: a type-1 class `T_c ∪ {b–c}` at a skeleton pair `{a, b}` of `U` was closed by the
  type-1 close / P-check at `b` before `U` formed: it is inside a sub-piece, outside `U`, or
  contains `U`;
* `bond`: two parallel edges of `U` lie in a common sub-piece (parallel pieces were merged by the
  P-check, so the skeleton has no multi-edge);
* `type2`: the type-2 class at a skeleton pair `{a, b}` of `U` (the class of the tree edge
  `a → a'`) was a `finishTstackTop` (`firstIdx > firstOccurrence[d]`) before `U` formed: it is
  laminar with `U` in the same sense. -/
structure RCloseShape (g : Graph) (d : DfsData) (P : Pieces) (U : Nat → Prop) (s t : Nat) :
    Prop where
  wf : P.WF g
  sub : ∀ i e, P.Mem i e → U e
  conn : g.ConnEdges U
  attached : g.TwoAttached U s t
  touch_s : g.Touches U s
  touch_t : g.Touches U t
  ne : s ≠ t
  proper : ∃ e, e < g.ne ∧ ¬U e
  single : ∀ e e', e < g.ne → e' < g.ne → U e → U e' → g.SepClass s t e e'
  maximal : ∀ i, i < P.k → ∀ e e', e < g.ne → e' < g.ne → ¬P.Mem i e → ¬P.Mem i e' →
    g.SepClass (P.x i) (P.y i) e e'
  type1 : ∀ a b, P.SkelPair g U s t a b → d.Anc a b → ∀ o ∈ d.outs b,
    o.cls = .ret (d.depth a) .type1Child →
      P.LaminarWith (fun e => e < g.ne ∧ U e) (d.EndIn o.dest · g)
  bond : ∀ a b e₁ e₂, e₁ ≠ e₂ → g.Joins e₁ a b → g.Joins e₂ a b → U e₁ → U e₂ →
    ∃ i, P.Mem i e₁ ∧ P.Mem i e₂
  type2 : ∀ a b, P.SkelPair g U s t a b → d.Type2Pair a b → ∀ o ∈ d.outs a, o.isTree = true →
    d.Anc o.dest b → P.LaminarWith (fun e => e < g.ne ∧ U e) (g.SepClass a b o.e)

end Spqr
