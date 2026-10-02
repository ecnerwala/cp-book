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
  q_root_attach : ∀ i c w e, i < t.size → t.type i = .Q →
    t.children i = [c, w] → t.type w = .V → t.origId[i]! = some e →
      ∀ v e', t.Touches g i v → Graph.Incident g v e' → ¬ t.EdgeIn i e' →
        v = (g.edges[e]!).1 ∨ v = (g.edges[e]!).2

end Spqr.SpqrTree
