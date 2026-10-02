import Mathlib.Logic.Relation
import Spqr.Relabel

/-!
# Specification: what a correct SPQR tree is

`SpqrTree.WF` is the structural well-formedness of the output (preorder tree, CSR bounds, index
bijections, skeleton shapes); `SpqrTree.Represents g` says the skeletons glue back to `g` along
the twin edges and every S/P/R node is a genuine triconnected component. The correctness theorems
are stated in `Spqr.Correctness`.
-/

namespace Spqr

namespace SpqrTree

variable (t : SpqrTree)

/-- Number of items. -/
def size : Nat := t.types.size

def type (i : Nat) : NodeType := t.types[i]?.getD .F
def parent (i : Nat) : Option Nat := t.par[i]?.getD none
def chRange (i : Nat) : Nat × Nat := (t.chBounds[i]?.getD 0, t.chBounds[i + 1]?.getD 0)
def nvRange (i : Nat) : Nat × Nat := (t.nvBounds[i]?.getD 0, t.nvBounds[i + 1]?.getD 0)
def neRange (i : Nat) : Nat × Nat := (t.neBounds[i]?.getD 0, t.neBounds[i + 1]?.getD 0)
/-- The children of item `i`, as listed in `chDat`. -/
def children (i : Nat) : List Nat :=
  ((List.range ((t.chRange i).2 - (t.chRange i).1)).map fun k => t.chDat[(t.chRange i).1 + k]?.getD 0)
def nodeVertsOf (i : Nat) : List NodeVert :=
  (List.range ((t.nvRange i).2 - (t.nvRange i).1)).map fun k => t.nodeVerts[(t.nvRange i).1 + k]?.getD default
def nodeEdgesOf (i : Nat) : List NodeEdge :=
  (List.range ((t.neRange i).2 - (t.neRange i).1)).map fun k => t.nodeEdges[(t.neRange i).1 + k]?.getD default
def nVerts (i : Nat) : Nat := (t.nvRange i).2 - (t.nvRange i).1
def nEdges (i : Nat) : Nat := (t.neRange i).2 - (t.neRange i).1

/-- Original vertex of node-vertex `nv` (via its V item). -/
def nvOrig (nv : Nat) : Option Nat := do
  let v ← t.nodeVerts[nv]?
  t.origId[v.vert]?.getD none
def nodeOfNv (nv : Nat) : Option Nat := (t.nodeVerts[nv]?).map (·.node)
def nodeOfNe (ne : Nat) : Option Nat := (t.nodeEdges[ne]?).map (·.node)
def twin (ne : Nat) : Option Nat := (t.nodeEdges[ne]?).bind (·.twin)
def nvsOf (ne : Nat) : Option (Nat × Nat) := (t.nodeEdges[ne]?).map (·.nvs)

/-- `i` is a proper ancestor-or-self of `j`: preorder with subtree ranges. -/
def inSubtree (i j : Nat) : Prop := i ≤ j ∧ j < t.subtreeEnd[i]?.getD 0

/-- A node (anything but F / V) has a cap edge — its first node-edge — unless it is a Q node that
is the root of its block (then it has children and no parent node). -/
def hasCap (i : Nat) : Bool :=
  (t.type i).isNode && !(t.type i == .Q && (t.children i).any fun c => t.type c != .V)

/-- Is `i` the root Q of a block (an edge whose block has no other edges)? -/
def isBlockRootQ (i : Nat) : Bool := t.type i == .Q && !t.hasCap i

/-- The cap node-edge of `i`. -/
def capNe (i : Nat) : Option Nat := if t.hasCap i then some (t.neRange i).1 else none

/-! ### Structural well-formedness -/

/-- Array sizes agree. -/
structure Sizes : Prop where
  par : t.par.size = t.size
  subtreeEnd : t.subtreeEnd.size = t.size
  origId : t.origId.size = t.size
  chBounds : t.chBounds.size = t.size + 1
  nvBounds : t.nvBounds.size = t.size + 1
  neBounds : t.neBounds.size = t.size + 1
  adjBounds : t.adjBounds.size = 2 * t.nodeVerts.size + 1
  vertParNv : t.vertParNv.size = t.size
  vertIndex : t.vertIndex.size = t.nv
  edgeIndex : t.edgeIndex.size = t.ne
  edgeFlipped : t.edgeFlipped.size = t.ne

/-- The items form a tree stored in preorder: parents precede children, children lists are the
exact fibers of `par`, and subtree ranges are contiguous. -/
structure Preorder : Prop where
  root_type : t.type 0 = .F
  root_par : t.parent 0 = none
  par_lt : ∀ i, 0 < i → i < t.size → ∃ p, t.parent i = some p ∧ p < i
  par_subtree : ∀ i p, t.parent i = some p → t.inSubtree p i
  ch_eq : ∀ i, i < t.size → t.children i = (List.range t.size).filter fun j => t.parent j = some i
  ch_monotone : ∀ i, i + 1 < t.chBounds.size → t.chBounds[i]! ≤ t.chBounds[i + 1]!
  subtree_eq : ∀ i, i < t.size → t.subtreeEnd[i]! = i + 1 + ((t.children i).map fun c => t.subtreeEnd[c]! - c).sum
  only_root_F : ∀ i, 0 < i → t.type i ≠ .F

/-- V items ↔ vertices and Q items ↔ edges, via `vertIndex` / `edgeIndex` / `origId`. -/
structure Bijections : Prop where
  vert_index : ∀ v, v < t.nv → ∃ i, t.vertIndex[v]! = some i ∧ t.type i = .V ∧ t.origId[i]! = some v
  vert_orig : ∀ i, i < t.size → t.type i = .V → ∃ v, t.origId[i]! = some v ∧ t.vertIndex[v]! = some i
  edge_index : ∀ e, e < t.ne → ∃ i, t.edgeIndex[e]! = some i ∧ t.type i = .Q ∧ t.origId[i]! = some e
  edge_orig : ∀ i, i < t.size → t.type i = .Q → ∃ e, t.origId[i]! = some e ∧ t.edgeIndex[e]! = some i
  other_orig : ∀ i, i < t.size → t.type i ≠ .V → t.type i ≠ .Q → t.origId[i]! = none

/-- Every node-vertex / node-edge is owned by the node whose range contains it and refers to a V
item / node-vertices of that node; `vertParNv` points each V item at its slot in its parent. -/
structure Ownership : Prop where
  nv_bounds_mono : ∀ i, i + 1 < t.nvBounds.size → t.nvBounds[i]! ≤ t.nvBounds[i + 1]!
  ne_bounds_mono : ∀ i, i + 1 < t.neBounds.size → t.neBounds[i]! ≤ t.neBounds[i + 1]!
  nv_last : t.nvBounds[t.size]! = t.nodeVerts.size
  ne_last : t.neBounds[t.size]! = t.nodeEdges.size
  nv_node : ∀ i nv, i < t.size → (t.nvRange i).1 ≤ nv → nv < (t.nvRange i).2 → t.nodeOfNv nv = some i
  nv_vert : ∀ nv, nv < t.nodeVerts.size → t.type (t.nodeVerts[nv]!.vert) = .V
  nv_distinct : ∀ i, i < t.size → ((t.nodeVertsOf i).map (·.vert)).Nodup
  ne_node : ∀ i ne, i < t.size → (t.neRange i).1 ≤ ne → ne < (t.neRange i).2 → t.nodeOfNe ne = some i
  ne_nvs : ∀ ne, ne < t.nodeEdges.size → ∀ i, t.nodeOfNe ne = some i →
    (t.nvRange i).1 ≤ t.nodeEdges[ne]!.nvs.1 ∧ t.nodeEdges[ne]!.nvs.1 ≤ t.nodeEdges[ne]!.nvs.2 ∧
    t.nodeEdges[ne]!.nvs.2 < (t.nvRange i).2
  vert_par_nv : ∀ i p, i < t.size → t.type i = .V → t.parent i = some p →
    ∃ nv, t.vertParNv[i]! = some nv ∧ (t.nvRange p).1 ≤ nv ∧ nv < (t.nvRange p).2 ∧ t.nodeVerts[nv]!.vert = i
  vert_par_nv_none : ∀ i, i < t.size → t.type i ≠ .V → t.vertParNv[i]! = none
  /-- Node-vertices of a node are: its first cap endpoint, its V children in order, its second cap
  endpoint (so V children are exactly the vertices "interior" to the node). -/
  nv_layout : ∀ i, i < t.size → ∃ a b : List NodeVert,
    (t.nodeVertsOf i) = a ++ ((t.children i).filter fun c => t.type c = .V).map (fun c => ⟨i, c⟩) ++ b ∧
    a.length ≤ 1 ∧ b.length ≤ 1 ∧
    (t.hasCap i → ∀ nvs, t.nvsOf (t.neRange i).1 = some nvs →
      (nvs.1 = (t.nvRange i).1 ∧ nvs.2 = (t.nvRange i).2 - 1))

/-- Twin edges pair the cap of every non-block-root node with one edge of its parent, bijectively
over the parent's non-cap edges; twins are an involution on node-edges. -/
structure Twins : Prop where
  twin_invol : ∀ ne ne', t.twin ne = some ne' → t.twin ne' = some ne
  twin_ne : ∀ ne ne', t.twin ne = some ne' → ne ≠ ne'
  twin_parent : ∀ i ne, i < t.size → t.capNe i = some ne →
    (∃ p, t.parent i = some p ∧ (t.type p).isNode ∧ ∃ ne', t.twin ne = some ne' ∧ t.nodeOfNe ne' = some p) ∨
    (t.twin ne = none ∧ ∃ p, t.parent i = some p ∧ ¬ (t.type p).isNode)
  /-- Every non-cap edge of a node is the twin of the cap of exactly one child, in children order
  (children that are V items are skipped). -/
  noncap_children : ∀ i, i < t.size → (t.type i).isNode →
    ((t.nodeEdgesOf i).drop (if t.hasCap i then 1 else 0)).map (·.twin) =
      ((t.children i).filter fun c => t.type c ≠ .V).map fun c => t.capNe c
  cap_none : ∀ i, i < t.size → ¬ (t.type i).isNode → t.nEdges i = 0

/-- Edges of node `i` as node-vertex pairs. -/
def skeleton (i : Nat) : List (Nat × Nat) := (t.nodeEdgesOf i).map (·.nvs)

/-- The expected skeleton of each node type. `S`: the cycle `nvSt, nvSt+1, …, nvEn-1` closed by the
cap; `P`: a bond of ≥ 3 edges on 2 vertices; `Q`/`I`: one edge; `O`: one loop; `R`: a simple graph
on ≥ 4 vertices with ≥ 6 edges (3-connectivity is in `Represents`). -/
def Shape (i : Nat) : Prop :=
  let (s, e) := t.nvRange i
  let n := e - s
  match t.type i with
  | .F => t.nEdges i = 0
  | .V => n = 0 ∧ t.nEdges i = 0
  | .Q => (n = 1 ∧ t.skeleton i = [(s, s)]) ∨ (n = 2 ∧ t.skeleton i = [(s, s + 1)])
  | .I => n = 2 ∧ t.skeleton i = [(s, s + 1)]
  | .O => n = 1 ∧ t.skeleton i = [(s, s)]
  | .P => n = 2 ∧ 3 ≤ t.nEdges i ∧ ∀ p ∈ t.skeleton i, p = (s, s + 1)
  | .S => 3 ≤ n ∧ t.skeleton i = (s, e - 1) :: (List.range (n - 1)).map fun k => (s + k, s + k + 1)
  | .R => 4 ≤ n ∧ 6 ≤ t.nEdges i ∧ (t.skeleton i).Nodup ∧ ∀ p ∈ t.skeleton i, p.1 < p.2

structure WF : Prop where
  sizes : t.Sizes
  preorder : t.Preorder
  bij : t.Bijections
  own : t.Ownership
  twins : t.Twins
  shape : ∀ i, i < t.size → t.Shape i
  adj_bounds_mono : ∀ i, i + 1 < t.adjBounds.size → t.adjBounds[i]! ≤ t.adjBounds[i + 1]!
  adj_last : t.adjBounds[2 * t.nodeVerts.size]! = t.adjDat.size
  /-- Rows `2 nv` and `2 nv + 1` of `adjDat` hold exactly the incidences of node-vertex `nv`. -/
  adj_incident : ∀ nv, nv < t.nodeVerts.size →
    ((List.range (t.adjBounds[2 * nv + 2]! - t.adjBounds[2 * nv]!)).map fun k =>
        (t.adjDat[t.adjBounds[2 * nv]! + k]!).ne).Perm
      ((List.range t.nodeEdges.size).filter fun ne =>
        t.nodeEdges[ne]!.nvs.1 = nv ∨ t.nodeEdges[ne]!.nvs.2 = nv)
  adj_dest : ∀ k, k < t.adjDat.size → ∀ nvs, t.nvsOf t.adjDat[k]!.ne = some nvs →
    t.adjDat[k]!.destNv = nvs.1 ∨ t.adjDat[k]!.destNv = nvs.2

/-! ### Semantic correctness: the tree is an SPQR decomposition of `g` -/

/-- Original endpoints of node-edge `ne`, in node-vertex order. -/
def neOrig (ne : Nat) : Option (Nat × Nat) := do
  let nvs ← t.nvsOf ne
  return (← t.nvOrig nvs.1, ← t.nvOrig nvs.2)

/-- Unordered equality of vertex pairs. -/
def PairEq (p q : Nat × Nat) : Prop := p = q ∨ p = (q.2, q.1)

/-- Original edge `e` lies in the subtree of item `i`. -/
def EdgeIn (i e : Nat) : Prop := ∃ j, t.edgeIndex[e]! = some j ∧ t.inSubtree i j

/-- Vertex `v` is incident to edge `e` of `g`. -/
def Graph.Incident (g : Graph) (v e : Nat) : Prop := e < g.ne ∧ ((g.edges[e]!).1 = v ∨ (g.edges[e]!).2 = v)

/-- A simple graph on `Fin n` given as a list of edges is 3-connected: `n ≥ 4` and removing any
two vertices leaves the rest connected. -/
def ThreeConnected (n : Nat) (es : List (Nat × Nat)) : Prop :=
  4 ≤ n ∧ ∀ a b, a < n → b < n →
    let alive := fun v => v < n ∧ v ≠ a ∧ v ≠ b
    let adj := fun u w => alive u ∧ alive w ∧ ((u, w) ∈ es ∨ (w, u) ∈ es)
    ∀ u w, alive u → alive w → Relation.ReflTransGen adj u w

/-- `t` is an SPQR decomposition of `g`. -/
structure Represents (g : Graph) : Prop where
  nv : t.nv = g.nv
  ne : t.ne = g.ne
  /-- Q items carry their edge's endpoints (possibly flipped). -/
  q_endpoints : ∀ e, e < g.ne → ∀ i, t.edgeIndex[e]! = some i →
    ∀ ne, (t.neRange i).1 = ne → ∀ p, t.neOrig ne = some p →
      (t.edgeFlipped[e]! = false → p = g.edges[e]!) ∧ (t.edgeFlipped[e]! = true → p = ((g.edges[e]!).2, (g.edges[e]!).1))
  /-- Twin edges connect the same two original vertices. -/
  twin_glue : ∀ ne ne', t.twin ne = some ne' → ∀ p q, t.neOrig ne = some p → t.neOrig ne' = some q → PairEq p q
  /-- Every node-vertex pair of a node refers to distinct original vertices, except for loops. -/
  nv_orig_inj : ∀ i, i < t.size → t.type i ≠ .O → t.type i ≠ .Q →
    ((t.nodeVertsOf i).map fun nv => t.origId[nv.vert]!).Nodup
  /-- Separation: the subtree below a capped node touches the rest of `g` only through its cap's
  endpoints; the subtree below a block-root Q (a whole block) only through cut vertices, i.e. every
  original vertex of `g` touching the subtree's edges and other edges lies in the cap. -/
  separation : ∀ i ne, i < t.size → t.capNe i = some ne → ∀ p, t.neOrig ne = some p →
    ∀ v e e', Graph.Incident g v e → Graph.Incident g v e' → t.EdgeIn i e → ¬ t.EdgeIn i e' →
      v = p.1 ∨ v = p.2
  /-- The vertex children of a node are exactly the original vertices interior to it: the
  vertices all of whose edges lie in its subtree but not in a single child's subtree. -/
  interior : ∀ i v, i < t.size → v < g.nv → ∀ j, t.vertIndex[v]! = some j →
    (t.parent j = some i ↔
      (∀ e, Graph.Incident g v e → t.EdgeIn i e) ∧
      ∀ c ∈ t.children i, ¬ ∀ e, Graph.Incident g v e → t.EdgeIn c e)
  /-- R skeletons are 3-connected (as graphs on their node-vertices, relabelled from `nvSt`). -/
  r_three_connected : ∀ i, i < t.size → t.type i = .R →
    ThreeConnected (t.nVerts i) ((t.skeleton i).map fun p => (p.1 - (t.nvRange i).1, p.2 - (t.nvRange i).1))
  /-- Canonicity: no two adjacent S nodes and no two adjacent P nodes unless ternarizing. -/
  canonical : ∀ i p, t.parent i = some p → (t.type i = .S → t.type p ≠ .S) ∧ (t.type i = .P → t.type p ≠ .P)

end SpqrTree

end Spqr
