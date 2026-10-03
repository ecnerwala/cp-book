import Spqr.SepPair

/-!
# The two kinds of separation pair (PROOF.md §3, §4.5)

`Type1Pair` and `Type2Pair` are the lowpoint-level descriptions of the `isType1` / type-2 splits
of `finishEdge`; `Spqr.Proofs.SepPairExhaust` proves they are exactly the separation pairs of a
block (`sepPair_iff`), so a piece with no split left is 3-connected (`three_connected_of_no_split`).
-/

namespace Spqr

namespace Graph

/-- No separation pair (in the SPQR convention of `SeparationPair`). -/
def ThreeConnected (g : Graph) : Prop := ∀ a b, ¬g.SeparationPair a b

end Graph

namespace DfsData

variable (d : DfsData)

/-- A type-1 split at `{a, b}`, `a` a proper ancestor of `b`: a type-1 child `b → c` to `depth a`
(its class is `T_c` with `b–c` and `T_c`'s back edges) leaving at least two other edges, or two
parallel `a–b` edges (back edges `b → a`, or the tree edge when `a` is `b`'s parent) next to a
third edge — a bond. -/
def Type1Pair (a b : Nat) (g : Graph) : Prop :=
  d.Anc a b ∧ a ≠ b ∧
    ((∃ o ∈ d.outs b, o.cls = .ret (d.depth a) .type1Child ∧
        ∃ e e', e < g.ne ∧ e' < g.ne ∧ e ≠ e' ∧ ¬d.EndIn o.dest e g ∧ ¬d.EndIn o.dest e' g) ∨
      (∃ e₁ e₂ e₃, e₃ < g.ne ∧ e₁ ≠ e₂ ∧ e₁ ≠ e₃ ∧ e₂ ≠ e₃ ∧ g.Joins e₁ a b ∧ g.Joins e₂ a b))

/-- No child subtree of `b` attaches both strictly above `a` and strictly between `a` and `b`. -/
def NoBothSides (a b : Nat) : Prop :=
  ∀ c, d.IsParent b c → (∃ l, l < d.depth a ∧ d.Returns c l) →
    ∀ l, d.depth a < l → l < d.depth b → ¬d.Returns c l

/-- No back edge out of `T_{a'} − T_b` lands strictly above `a`. -/
def BetweenStays (a a' b : Nat) : Prop :=
  ∀ u, d.Anc a' u → ¬d.Anc b u → ∀ o ∈ d.outs u, o.isTree = false → d.depth a ≤ d.depth o.dest

/-- A type-2 split at `{a, b}`: `a` is not the root, `b` is a proper descendant of `a`'s child
`a'`, and the *above* part (outside `T_{a'} ∪ {a}`, with the children of `b` returning above `a`)
is not joined to the *between* part (`T_{a'} − T_b`, with the other children of `b`). -/
def Type2Pair (a b : Nat) : Prop :=
  ∃ q a', d.IsParent q a ∧ d.IsParent a a' ∧ d.Anc a' b ∧ a' ≠ b ∧
    d.NoBothSides a b ∧ d.BetweenStays a a' b

end DfsData

end Spqr
