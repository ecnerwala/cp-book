import Spqr.Spec

namespace Spqr.SpqrTree

variable (t : SpqrTree) (g : Graph)

/-- A vertex incident to an original edge below an item. -/
def Touches (i v : Nat) : Prop :=
  ∃ e, Graph.Incident g v e ∧ t.EdgeIn i e

/-- Node-vertex `nv` is an endpoint of the cap of node `i`. -/
def CapEnd (i nv : Nat) : Prop :=
  ∃ ne d, t.capNe i = some ne ∧ t.nodeEdges[ne]? = some d ∧ (d.nvs.1 = nv ∨ d.nvs.2 = nv)

/-- `nv` is a node-vertex of item `i`. -/
def NvOf (i nv : Nat) : Prop := (t.nvRange i).1 ≤ nv ∧ nv < (t.nvRange i).2

/-- Child `c` of node `i` is attached at node-vertex `nv`: `c` is the `V` item of `nv`, or the twin
of `c`'s cap is a node-edge of `i` with endpoint `nv`. -/
def NvInc (i c nv : Nat) : Prop :=
  (∃ d, t.nodeVerts[nv]? = some d ∧ d.vert = c) ∨
    ∃ ne tw d, t.nodeOfNe ne = some i ∧ t.twin ne = some tw ∧ t.capNe c = some tw ∧
      t.nodeEdges[ne]? = some d ∧ (d.nvs.1 = nv ∨ d.nvs.2 = nv)

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
  /-- Twin node-edges list their original endpoints in the same order (the node step reads the
  child's exposed slot `2 * side + dir` for its own quarter-edge `(side, dir)`). -/
  twin_orient : ∀ ne tw, t.twin ne = some tw → t.neOrig ne = t.neOrig tw
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
  /-- The V item of a node-vertex of an S/P/R node is a child of the node exactly when the
  node-vertex is not a cap endpoint (the node step splices it at that vertex's corner). -/
  nv_child : ∀ i nv d, i < t.size → (t.type i = .S ∨ t.type i = .P ∨ t.type i = .R) →
    t.NvOf i nv → t.nodeVerts[nv]? = some d → (t.parent d.vert = some i ↔ ¬ t.CapEnd i nv)
  /-- Every V child of an S/P/R node is the V item of one of its node-vertices. -/
  v_child_nv : ∀ i c, i < t.size → (t.type i = .S ∨ t.type i = .P ∨ t.type i = .R) →
    c ∈ t.children i → t.type c = .V → ∃ nv d, t.NvOf i nv ∧ t.nodeVerts[nv]? = some d ∧ d.vert = c
  /-- An S/P/R node hangs under a `Q` (block root) or another node, never directly under the root
  or a `V` item (so its exposed row has four slots). -/
  node_parent : ∀ i p, i < t.size → (t.type i = .S ∨ t.type i = .P ∨ t.type i = .R) →
    t.parent i = some p → t.type p ≠ .F ∧ t.type p ≠ .V
  /-- Distinct node-vertices of an S/P/R node have distinct V items. -/
  nv_vert_inj : ∀ i nv nv' d d', i < t.size → (t.type i = .S ∨ t.type i = .P ∨ t.type i = .R) →
    t.NvOf i nv → t.NvOf i nv' → t.nodeVerts[nv]? = some d → t.nodeVerts[nv']? = some d' →
    d.vert = d'.vert → nv = nv'
  /-- Every non-`V` child of an S/P/R node is a capped item that is not an `I`/`O` leaf. -/
  node_child_cap : ∀ i c, i < t.size → (t.type i = .S ∨ t.type i = .P ∨ t.type i = .R) →
    c ∈ t.children i → t.type c ≠ .V → t.hasCap c = true ∧ t.type c ≠ .I ∧ t.type c ≠ .O
  /-- Distinct node-vertices of an S/P/R node have distinct original vertices. -/
  nv_orig_inj : ∀ i nv nv', i < t.size → (t.type i = .S ∨ t.type i = .P ∨ t.type i = .R) →
    t.NvOf i nv → t.NvOf i nv' → t.nvOrig nv = t.nvOrig nv' → nv = nv'
  /-- A child of an S/P/R node touches the original vertex of a node-vertex only where it is
  attached (`NvInc`). -/
  node_touch : ∀ i c w nv, i < t.size → (t.type i = .S ∨ t.type i = .P ∨ t.type i = .R) →
    c ∈ t.children i → t.Touches g c w → t.NvOf i nv → t.nvOrig nv = some w → t.NvInc i c nv
  /-- Two distinct children of an S/P/R node meet only at original vertices of its node-vertices. -/
  node_attach : ∀ i a b w, i < t.size → (t.type i = .S ∨ t.type i = .P ∨ t.type i = .R) →
    a ∈ t.children i → b ∈ t.children i → a ≠ b → t.Touches g a w → t.Touches g b w →
    ∃ nv, t.NvOf i nv ∧ t.nvOrig nv = some w

end Spqr.SpqrTree
