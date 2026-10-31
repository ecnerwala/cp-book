import Spqr.GraphLemmas
import Spqr.SepPairExhaust

/-!
# Contracting pieces to a skeleton (PROOF.md §4.5, §5)

A family of pairwise edge-disjoint *pieces* of `g` — each a connected edge set 2-attached at its
two terminals — is contracted to the skeleton `Pieces.contract`: the edges of `g` in no piece, plus
one virtual edge per piece between its terminals. A laminar family is handled by passing its
maximal members (nested pieces are absorbed, and their terminals are then interior or terminal
vertices of the enclosing piece), so disjointness is built into the representation
(`piece : Nat → Option Nat`). `Spqr.Proofs.Contract` relates separation pairs of the skeleton to
those of `g`.
-/

namespace Spqr

/-- `k` pairwise edge-disjoint pieces of a graph: `piece e = some i` puts edge `e` into piece `i`,
whose terminals are `x i` and `y i`. -/
structure Pieces where
  k : Nat
  piece : Nat → Option Nat
  x : Nat → Nat
  y : Nat → Nat

namespace Pieces

variable (P : Pieces) (g : Graph)

/-- Edge `e` lies in piece `i`. -/
def Mem (i e : Nat) : Prop := P.piece e = some i

/-- The skeleton's edges by origin: the edges of `g` in no piece (in order), then one virtual edge
per piece. -/
def origins : List (Nat ⊕ Nat) :=
  ((List.range g.ne).filter fun e => (P.piece e).isNone).map Sum.inl ++
    (List.range P.k).map Sum.inr

def originEnds : Nat ⊕ Nat → Nat × Nat
  | .inl e => g.edges[e]!
  | .inr i => (P.x i, P.y i)

/-- The skeleton: `g` with every piece replaced by a virtual edge between its terminals. -/
def contract : Graph := ⟨g.nv, ((P.origins g).map (P.originEnds g)).toArray⟩

/-- Skeleton edge `f` stands for `s`: the edge `e` of `g` (`.inl e`) or the piece `i` (`.inr i`). -/
def Orig (f : Nat) (s : Nat ⊕ Nat) : Prop := (P.origins g)[f]? = some s

/-- Skeleton edge `f` is the image of edge `e` of `g`: its copy if `e` is in no piece, else the
virtual edge of `e`'s piece. -/
def Img (e f : Nat) : Prop :=
  (P.piece e = none ∧ P.Orig g f (.inl e)) ∨ ∃ i, P.Mem i e ∧ P.Orig g f (.inr i)

/-- `w` is a vertex of piece `i` other than its terminals. -/
def Int (i w : Nat) : Prop := i < P.k ∧ g.Touches (P.Mem i) w ∧ w ≠ P.x i ∧ w ≠ P.y i

/-- A skeleton vertex: interior to no piece. -/
def Skel (w : Nat) : Prop := ∀ i, ¬P.Int g i w

/-- `{a, b}` is the terminal pair of some piece. -/
def TermPair (a b : Nat) : Prop :=
  ∃ i, i < P.k ∧ ((a = P.x i ∧ b = P.y i) ∨ (a = P.y i ∧ b = P.x i))

/-- Well-formed pieces: in range, connected, 2-attached at their terminals, touching both
terminals, with distinct terminals and not the whole edge set. -/
structure WF : Prop where
  lt : ∀ i e, P.Mem i e → i < P.k ∧ e < g.ne
  conn : ∀ i, i < P.k → g.ConnEdges (P.Mem i)
  attached : ∀ i, i < P.k → g.TwoAttached (P.Mem i) (P.x i) (P.y i)
  touch : ∀ i, i < P.k → g.Touches (P.Mem i) (P.x i) ∧ g.Touches (P.Mem i) (P.y i)
  ne : ∀ i, i < P.k → P.x i ≠ P.y i
  proper : ∀ i, i < P.k → ∃ e, e < g.ne ∧ ¬P.Mem i e

end Pieces

end Spqr
