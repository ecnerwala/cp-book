import Spqr.Graph

namespace Spqr

/-- `F` = forest root, `V` = original vertex, `Q` = original edge, `I` = bridge, `O` = self-loop,
`S` = series (cycle), `P` = parallel (bond), `R` = rigid (3-connected). -/
inductive NodeType where
  | F | V | Q | I | O | S | P | R
deriving DecidableEq, Repr, Inhabited

/-- Items are the SPQR nodes together with the original vertices, all living in one tree.
Item ids: `0` is the root, `1 + v` is vertex `v`, `1 + nv + e` is edge `e`, and the remaining
nodes are allocated after that. -/
abbrev ItemId := Nat

def rootItem : ItemId := 0
def vertItem (v : Nat) : ItemId := 1 + v
def edgeItem (g : Graph) (e : Nat) : ItemId := 1 + g.nv + e

structure Item where
  type : NodeType
  /-- The endpoints of the node's cap edge (`(some v, none)` for the one-vertex cases). -/
  vs : Option Nat × Option Nat := (none, none)
  /-- Children, in insertion order. -/
  ch : List ItemId := []
deriving Repr, Inhabited

def initialItems (g : Graph) : Array Item :=
  #[⟨.F, (none, none), []⟩]
    ++ Array.replicate g.nv ⟨.V, (none, none), []⟩
    ++ Array.replicate g.ne ⟨.Q, (none, none), []⟩

end Spqr
