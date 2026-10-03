import Spqr.RMax
import Spqr.Frame
import Spqr.WalkSpec

/-!
# The R case of loop 1, walk side (PROOF.md §4.3 "R case", §4.5)

`loop1Body` closes an R item when `nxt` returns to the same depth `d` as `cur` from a different
`vStart`: a fresh `R` item is allocated, `cur` is merged into `nxt`, and the merged entry is closed
into the item (`finishTstackTop`). The closed edge set is `U = edges cur ∪ edges nxt`
(`WalkState.rU`), its terminals are `nxt.vStart` and `stackVerts[d]`, and its sub-pieces are the
non-`V` items of the merged spans (`WalkState.rPieceItems`, read as `Pieces.ofItems`).

`WalkState.RStep` is the stack shape at such a step, as far as `Inv` and the loop-1 guards describe
it; `WalkState.RContent` names the five content facts of `RCloseShape` on that state. The content
facts are history facts of the walk (what earlier closes did), so they are hypotheses here; the
walk invariant that discharges them is the ear-structure invariant of `EarSpec.lean`.
-/

namespace Spqr

open Classical in
/-- The pieces given by a list of items: edge `e` is in the piece of the first item of `L` with
`e` below it, with terminals the item's `vs`. (For pairwise edge-disjoint items the order is
irrelevant: `Pieces.ofItems_mem_iff`.) -/
noncomputable def Pieces.ofItems (g : Graph) (items : Items) (L : List ItemId) : Pieces where
  k := L.length
  piece e := if e < g.ne then L.findIdx? fun i => decide (items.EdgeBelow g i e) else none
  x i := (Items.vs items L[i]!).1.getD 0
  y i := (Items.vs items L[i]!).2.getD 0

namespace WalkState

/-- The edges closed by the R case: the union of the two top entries. -/
def rU (s : WalkState) (cur nxt : TEntry) (e : Nat) : Prop :=
  cur.edges s.g s.items e ∨ nxt.edges s.g s.items e

/-- All items of `cur` merged into `nxt`, both sides. -/
def rItems (cur nxt : TEntry) : List ItemId :=
  (cur.spans.1 ++ nxt.spans.1) ++ (nxt.spans.2 ++ cur.spans.2)

/-- The sub-pieces of the R close: the merged items other than vertex items. -/
def rPieceItems (s : WalkState) (cur nxt : TEntry) : List ItemId :=
  (rItems cur nxt).filter fun i => decide (Items.type s.items i ≠ .V)

/-- A list of items that are pairwise edge-disjoint 2-terminal pieces. For node items `conn` and
`attached` are `Inv.nodes`; for the `Q` items of back edges they are the single-edge facts. -/
structure PieceItems (s : WalkState) (L : List ItemId) : Prop where
  nodup : L.Nodup
  vs : ∀ i ∈ L, ∃ x y, Items.vs s.items i = (some x, some y)
  ne : ∀ i ∈ L, ∃ e, e < s.g.ne ∧ Items.EdgeBelow s.g s.items i e
  conn : ∀ i ∈ L, s.g.ConnEdges (Items.EdgeBelow s.g s.items i)
  attached : ∀ i ∈ L, ∀ x y, Items.vs s.items i = (some x, some y) →
    s.g.TwoAttached (Items.EdgeBelow s.g s.items i) x y
  disj : ∀ i ∈ L, ∀ j ∈ L, i ≠ j → ∀ e, e < s.g.ne → Items.EdgeBelow s.g s.items i e →
    ¬Items.EdgeBelow s.g s.items j e

/-- The stack at an R step of loop 1 (`loop1Type` returns `.R`): `cur` and `nxt` both return to
depth `d` from different vertices, `cur`'s bottom has become interior to the union (PROOF.md §4.3,
statement 3: `cur.vStart` is no longer on the DFS path), both entries are nonempty, the union is
not the whole edge set, and the merged items are pieces. -/
structure RStep (s : WalkState) (d : Nat) (cur nxt : TEntry) (rest : List TEntry) : Prop where
  tstack : s.tstack = cur :: nxt :: rest
  cur_top : cur.topDepth = d
  nxt_top : nxt.topDepth = d
  ne : nxt.vStart ≠ cur.vStart
  interior : s.g.Interior (s.rU cur nxt) cur.vStart
  cur_ne : ∃ e, e < s.g.ne ∧ cur.edges s.g s.items e
  nxt_ne : ∃ e, e < s.g.ne ∧ nxt.edges s.g s.items e
  proper : ∃ e, e < s.g.ne ∧ ¬s.rU cur nxt e
  pieces : s.PieceItems (s.rPieceItems cur nxt)

/-- The content fields of `RCloseShape` at an R step (`U = rU`, `P = ofItems (rPieceItems)`,
`s = nxt.vStart`, `t = stackVerts[d]`), with `dfs` the sorted DFS tree of the block `s.g`. -/
structure RContent (s : WalkState) (dfs : DfsData) (d : Nat) (cur nxt : TEntry) : Prop where
  single : ∀ e e', e < s.g.ne → e' < s.g.ne → s.rU cur nxt e → s.rU cur nxt e' →
    s.g.SepClass nxt.vStart s.stackVerts[d]! e e'
  maximal : ∀ i, i < (Pieces.ofItems s.g s.items (s.rPieceItems cur nxt)).k →
    ∀ e e', e < s.g.ne → e' < s.g.ne →
      ¬(Pieces.ofItems s.g s.items (s.rPieceItems cur nxt)).Mem i e →
      ¬(Pieces.ofItems s.g s.items (s.rPieceItems cur nxt)).Mem i e' →
      s.g.SepClass ((Pieces.ofItems s.g s.items (s.rPieceItems cur nxt)).x i)
        ((Pieces.ofItems s.g s.items (s.rPieceItems cur nxt)).y i) e e'
  type1 : ∀ a b, (Pieces.ofItems s.g s.items (s.rPieceItems cur nxt)).SkelPair s.g (s.rU cur nxt)
      nxt.vStart s.stackVerts[d]! a b →
    dfs.Anc a b → ∀ o ∈ dfs.outs b, o.cls = .ret (dfs.depth a) .type1Child →
      (Pieces.ofItems s.g s.items (s.rPieceItems cur nxt)).LaminarWith
        (fun e => e < s.g.ne ∧ s.rU cur nxt e)
        (dfs.EndIn o.dest · s.g)
  bond : ∀ a b e₁ e₂, e₁ ≠ e₂ → s.g.Joins e₁ a b → s.g.Joins e₂ a b →
    s.rU cur nxt e₁ → s.rU cur nxt e₂ →
    ∃ i, (Pieces.ofItems s.g s.items (s.rPieceItems cur nxt)).Mem i e₁ ∧
      (Pieces.ofItems s.g s.items (s.rPieceItems cur nxt)).Mem i e₂
  type2 : ∀ a b, (Pieces.ofItems s.g s.items (s.rPieceItems cur nxt)).SkelPair s.g (s.rU cur nxt)
      nxt.vStart s.stackVerts[d]! a b →
    dfs.Type2Pair a b → ∀ o ∈ dfs.outs a, o.isTree = true → dfs.Anc o.dest b →
      (Pieces.ofItems s.g s.items (s.rPieceItems cur nxt)).LaminarWith
        (fun e => e < s.g.ne ∧ s.rU cur nxt e)
        (s.g.SepClass a b o.e)

end WalkState

end Spqr
