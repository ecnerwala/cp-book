import Spqr.Spec

namespace Spqr.SpqrTree

variable (t : SpqrTree) (g : Graph)

/-- A vertex incident to an original edge below an item. -/
def Touches (i v : Nat) : Prop :=
  ∃ e, Graph.Incident g v e ∧ t.EdgeIn i e

/-- Separation of the uncapped pieces used in the F/V/Q gluing steps. -/
structure PieceSep : Prop where
  root_disjoint : ∀ a ∈ t.children 0, ∀ b ∈ t.children 0, a ≠ b →
    ∀ v, t.Touches g a v → t.Touches g b v → False
  v_attach : ∀ i v, i < t.size → t.type i = .V → t.origId[i]! = some v →
    ∀ a ∈ t.children i, ∀ b ∈ t.children i, a ≠ b →
      ∀ w, t.Touches g a w → t.Touches g b w → w = v
  v_nonempty : ∀ i, i < t.size → t.type i = .V →
    ∀ c ∈ t.children i, ∃ e, t.EdgeIn c e
  /-- Every cap edge has original endpoints. -/
  cap_orig : ∀ c ne, c < t.size → t.capNe c = some ne → ∃ p, t.neOrig ne = some p
  /-- A capped item other than an `I` or `O` leaf has an edge below it. -/
  cap_nonempty : ∀ c, c < t.size → t.hasCap c = true → t.type c ≠ .I → t.type c ≠ .O →
    ∃ e, t.EdgeIn c e
  q_root_attach : ∀ i c w e, i < t.size → t.type i = .Q →
    t.children i = [c, w] → t.type w = .V → t.origId[i]! = some e →
      ∀ v e', t.Touches g i v → Graph.Incident g v e' → ¬ t.EdgeIn i e' →
        v = (g.edges[e]!).1 ∨ v = (g.edges[e]!).2
  /-- A `Q` item is a capped leaf, the root of a loop block (one `O` child), or the root of a
  block with a capped child `c` (an `I` leaf, a capped `Q` leaf or a node) followed by the lower
  `V` item `w`. -/
  q_shape : ∀ i, i < t.size → t.type i = .Q →
    t.children i = [] ∨ (∃ c, t.type c = .O ∧ t.children i = [c]) ∨
      ∃ c w, t.type c ≠ .V ∧ t.type c ≠ .O ∧ t.hasCap c = true ∧
        t.type w = .V ∧ t.children i = [c, w]
  /-- A capped `Q` leaf hangs under a node (never directly under the root or a `V` item). -/
  q_leaf_parent : ∀ i p, i < t.size → t.type i = .Q → t.children i = [] → t.parent i = some p →
    t.type p ≠ .F ∧ t.type p ≠ .V
  /-- A `Q` edge is a loop exactly when it has an `O` child. -/
  q_loop : ∀ i e, i < t.size → t.type i = .Q → t.origId[i]! = some e →
    ((g.edges[e]!).1 = (g.edges[e]!).2 ↔ ∃ c ∈ t.children i, t.type c = .O)
  /-- A `Q` hanging under a `V` item hangs at the edge's first endpoint in walk orientation. -/
  q_upper : ∀ i p e, i < t.size → t.type i = .Q → t.parent i = some p → t.type p = .V →
    t.origId[i]! = some e →
    t.origId[p]! = some (if t.edgeFlipped[e]! then (g.edges[e]!).2 else (g.edges[e]!).1)
  /-- The lower `V` item of a block-root `Q` is the edge's second endpoint in walk orientation. -/
  q_lower : ∀ i c w e, i < t.size → t.type i = .Q → t.children i = [c, w] → t.origId[i]! = some e →
    t.origId[w]! = some (if t.edgeFlipped[e]! then (g.edges[e]!).1 else (g.edges[e]!).2)
  /-- The cap of a block-root `Q`'s node child has the edge's endpoints in walk orientation. -/
  q_cap_orient : ∀ i c e ne p, i < t.size → t.type i = .Q → c ∈ t.children i →
    t.origId[i]! = some e → t.capNe c = some ne → t.neOrig ne = some p →
    p = if t.edgeFlipped[e]! then ((g.edges[e]!).2, (g.edges[e]!).1) else g.edges[e]!
  /-- The lower `V` piece of a block-root `Q` meets the edge and the node child's piece only at
  the lower endpoint. -/
  q_lower_attach : ∀ i c w e, i < t.size → t.type i = .Q → t.children i = [c, w] →
    t.origId[i]! = some e → ∀ v, (t.Touches g c v ∨ Graph.Incident g v e) → t.Touches g w v →
      v = if t.edgeFlipped[e]! then (g.edges[e]!).1 else (g.edges[e]!).2

end Spqr.SpqrTree
