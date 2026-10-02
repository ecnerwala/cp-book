import Spqr.RClose

/-!
# The R-maximality walk invariant (PROOF.md §4.5, walk side)

`WalkState.RContent` and the non-shape fields of `WalkState.RStep` are facts about the *history*
of the walk: what earlier closes (type-1 closes, P-checks, Loop 2's `firstIdx > firstOccurrence`
merges) did to the entries that Loop 1 now merges into an R item. This file states them as an
invariant of the tstack:

* `EntryR s dfs t` — the content facts of one open entry `t`, with terminals `a = t.vStart`,
  `b = stackVerts[t.topDepth]`: its closed items are pairwise edge-disjoint 2-terminal pieces,
  each a *maximal* piece (its complement is one class of its terminal pair: every P-check / Loop-1
  P merge at that pair has fired), parallel edges of `t` are bonded into a common piece (the
  P-check ran), an entry that is still attached at a third vertex (an open type-2 sub-ear) is a
  single `{a, b}`-class (Loop 2 glued every piece that returned to `b` from inside the sub-ear),
  and the type-1 / type-2 classes at skeleton pairs of `t` are laminar with `t` (the type-1 close
  at that vertex, resp. the `firstIdx > firstOccurrence` split, ran before `t` formed).
* `RInv s dfs` — every entry of the tstack satisfies `EntryR`, and distinct entries are
  edge-disjoint.
* `RBranch s d cur nxt rest` — the stack shape at Loop 1's R branch beyond what `Inv (d+1)` says:
  `cur`'s bottom is the child `stackVerts[d+1]` being finished and is interior to the union
  (statement 3 of Invariant W), `cur`'s items are between that child and `stackVerts[d]`, `nxt`
  touches both its terminals, no edge of `nxt` joins `cur`'s terminals (those are `cur`'s own
  parallel edges, P-merged into `cur`), both entries are nonempty and the union is proper.

`Proofs/RInv.lean` derives `RStep` and `RContent` from `RTop ∧ RBranch ∧ Inv (d+1)` (`RTop`: the
two top entries' `EntryR` and disjointness, all the R branch reads), so `RStep.threeConnected'`
gives the 3-connected R skeleton. `RInv` is not an invariant of every intermediate state (the
fresh tree-edge entry of `pushEdgeTstack`, or `cur` just before a P merge, is not `maximal`);
`RInvAt` is the form settled at the current vertex.
-/

namespace Spqr

namespace WalkState

/-- The sub-pieces of one entry: its non-`V` items. -/
def entryPieceItems (s : WalkState) (t : TEntry) : List ItemId :=
  (t.spans.1 ++ t.spans.2).filter fun i => decide (Items.type s.items i ≠ .V)

/-- `K` is laminar with the entry `t`: inside one piece of `t`, disjoint from `t`, or containing
`t` (`Pieces.LaminarWith` at item level). -/
def EntryLaminar (s : WalkState) (t : TEntry) (K : Nat → Prop) : Prop :=
  (∃ i ∈ s.entryPieceItems t, ∀ e, K e → Items.EdgeBelow s.g s.items i e) ∨
    (∀ e, K e → ¬t.edges s.g s.items e) ∨ ∀ e, t.edges s.g s.items e → K e

/-- A skeleton pair relative to the entry `t`: `b` is interior to no piece of `t`, and `{a, b}` is
neither a piece's terminal pair nor `t`'s own terminals (`Pieces.SkelPair` at item level). Neither
`a` nor `b` need be a vertex of `t`: a class at a pair outside `t` is trivially laminar with a
sub-ear (it contains the sub-ear or misses it), so the laminarity fields below are stated for all
such pairs. -/
def EntrySkelPair (s : WalkState) (t : TEntry) (a b : Nat) : Prop :=
  (∀ i ∈ s.entryPieceItems t, ∀ x y, Items.vs s.items i = (some x, some y) →
      s.g.Touches (Items.EdgeBelow s.g s.items i) b → b = x ∨ b = y) ∧
    (∀ i ∈ s.entryPieceItems t, ∀ x y, Items.vs s.items i = (some x, some y) →
      ¬((a = x ∧ b = y) ∨ (a = y ∧ b = x))) ∧
    ¬((a = t.vStart ∧ b = s.stackVerts[t.topDepth]!) ∨
      (a = s.stackVerts[t.topDepth]! ∧ b = t.vStart))

/-- The content facts of one open entry `t` (terminals `t.vStart`, `stackVerts[t.topDepth]`),
with `dfs` the sorted DFS tree of the block `s.g`. -/
structure EntryR (s : WalkState) (dfs : DfsData) (t : TEntry) : Prop where
  /-- The closed items of `t` are pairwise edge-disjoint 2-terminal pieces. -/
  pieces : s.PieceItems (s.entryPieceItems t)
  /-- Each closed item is a maximal piece: its complement is one class of its terminal pair. -/
  maximal : ∀ i ∈ s.entryPieceItems t, ∀ x y, Items.vs s.items i = (some x, some y) →
    ∀ e e', e < s.g.ne → e' < s.g.ne → ¬Items.EdgeBelow s.g s.items i e →
      ¬Items.EdgeBelow s.g s.items i e' → s.g.SepClass x y e e'
  /-- Parallel edges of `t` lie in a common closed item (the P-check bonded them). -/
  bond : ∀ a b e₁ e₂, e₁ ≠ e₂ → s.g.Joins e₁ a b → s.g.Joins e₂ a b →
    t.edges s.g s.items e₁ → t.edges s.g s.items e₂ →
    ∃ i ∈ s.entryPieceItems t, Items.EdgeBelow s.g s.items i e₁ ∧ Items.EdgeBelow s.g s.items i e₂
  /-- An entry still attached at a third vertex (an open type-2 sub-ear) is one class of its
  terminal pair: Loop 2 glued every piece that returned to the top from inside the sub-ear. (A
  2-attached entry may be a union of classes — the P case — so nothing is claimed for it.) -/
  single : ¬s.g.TwoAttached (t.edges s.g s.items) t.vStart s.stackVerts[t.topDepth]! →
    ∀ e e', e < s.g.ne → e' < s.g.ne → t.edges s.g s.items e → t.edges s.g s.items e' →
      s.g.SepClass t.vStart s.stackVerts[t.topDepth]! e e'
  /-- The class of a type-1 child edge at a skeleton pair of `t` is laminar with `t`: the type-1
  close / P-check at `b` ran before `t` formed. -/
  type1 : ∀ a b, s.EntrySkelPair t a b → dfs.Anc a b → ∀ o ∈ dfs.outs b,
    o.cls = .ret (dfs.depth a) .type1Child → s.EntryLaminar t (dfs.EndIn o.dest · s.g)
  /-- The class of the tree edge towards `b` at a type-2 skeleton pair of `t` is laminar with
  `t`: the `firstIdx > firstOccurrence` split fired before `t` formed. -/
  type2 : ∀ a b, s.EntrySkelPair t a b → dfs.Type2Pair a b → ∀ o ∈ dfs.outs a,
    o.isTree = true → dfs.Anc o.dest b → s.EntryLaminar t (s.g.SepClass a b o.e)

/-- The R-maximality invariant of the tstack. -/
structure RInv (s : WalkState) (dfs : DfsData) : Prop where
  entries : ∀ t ∈ s.tstack, s.EntryR dfs t
  disj : s.tstack.Pairwise fun t t' => ∀ e, t.edges s.g s.items e → ¬t'.edges s.g s.items e

/-- The content facts Loop 1's R branch reads: `EntryR` of its two top entries and their
edge-disjointness. -/
structure RTop (s : WalkState) (dfs : DfsData) (cur nxt : TEntry) : Prop where
  entry_cur : s.EntryR dfs cur
  entry_nxt : s.EntryR dfs nxt
  disj : ∀ e, cur.edges s.g s.items e → ¬nxt.edges s.g s.items e

/-- `RInv` settled at the vertex `v` being walked. An entry with bottom `v` is still collecting
`v`'s classes at its top: a `(v, l)` entry's complement contains the `(v, l)` classes of `v`'s
later out-edges until the P-check merges them (two parallel back edges `v → stackVerts[l]`: after
the first is pushed its single-edge entry is not `maximal`), so such entries are exempt until `v`
is finished (`walkTree_rInvAt`). -/
structure RInvAt (s : WalkState) (dfs : DfsData) (v : Nat) : Prop where
  entries : ∀ t ∈ s.tstack, t.vStart ≠ v → s.EntryR dfs t
  disj : s.tstack.Pairwise fun t t' => ∀ e, t.edges s.g s.items e → ¬t'.edges s.g s.items e

/-- The stack shape at Loop 1's R branch (`loop1Type` returns `.R` while `finishEdge` runs at
depth `d` for the tree edge `stackVerts[d] → stackVerts[d+1]`), beyond `Inv (d+1)`. -/
structure RBranch (s : WalkState) (d : Nat) (cur nxt : TEntry) (rest : List TEntry) : Prop where
  tstack : s.tstack = cur :: nxt :: rest
  cur_top : cur.topDepth = d
  nxt_top : nxt.topDepth = d
  ne : nxt.vStart ≠ cur.vStart
  /-- `cur`'s bottom is the child being finished. -/
  cur_c : cur.vStart = s.stackVerts[d + 1]!
  /-- `cur` has a closed item, and all its items are between the child and `stackVerts[d]`. -/
  cur_piece : ∃ i, i ∈ s.entryPieceItems cur
  cur_vs : ∀ i ∈ s.entryPieceItems cur,
    Items.vs s.items i = (some cur.vStart, some s.stackVerts[d]!) ∨
      Items.vs s.items i = (some s.stackVerts[d]!, some cur.vStart)
  /-- Statement 3 of Invariant W: the child is finished, all its edges are in the union. -/
  interior : s.g.Interior (s.rU cur nxt) cur.vStart
  cur_ne : ∃ e, e < s.g.ne ∧ cur.edges s.g s.items e
  nxt_ne : ∃ e, e < s.g.ne ∧ nxt.edges s.g s.items e
  proper : ∃ e, e < s.g.ne ∧ ¬s.rU cur nxt e
  /-- `nxt` touches both its terminals (it was pushed at its bottom and returned to its top). -/
  nxt_touch_top : s.g.Touches (nxt.edges s.g s.items) s.stackVerts[d]!
  nxt_touch_bot : s.g.Touches (nxt.edges s.g s.items) nxt.vStart
  /-- Edges joining the child to `stackVerts[d]` are the child's own back edges, P-merged into
  `cur`, never in `nxt`. -/
  nxt_no_cu : ∀ e, nxt.edges s.g s.items e → ¬s.g.Joins e cur.vStart s.stackVerts[d]!

end WalkState

end Spqr
