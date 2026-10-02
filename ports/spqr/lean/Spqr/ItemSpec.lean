import Mathlib.Logic.Relation
import Spqr.Walk

/-!
# Item-level specification (the interface between phases 2 and 3)

Phase 2 produces an `Array Item`: the root `F`, one `V` per vertex, one `Q` per edge, and the
allocated S/P/R/I/O nodes, linked by `Item.ch`. These predicates say what the walk must
guarantee so that `relabelTree` produces a well-formed `SpqrTree` representing `g`.
-/

namespace Spqr

abbrev Items := Array Item

namespace Items

variable (g : Graph) (items : Items)

def type (i : ItemId) : NodeType := (items[i]?.map (·.type)).getD .F
def ch (i : ItemId) : List ItemId := (items[i]?.map (·.ch)).getD []
def vs (i : ItemId) : Option Nat × Option Nat := (items[i]?.map (·.vs)).getD (none, none)

/-- `p` is the parent of `c`. -/
def IsParent (p c : ItemId) : Prop := c ∈ items.ch p

/-- `i` is a descendant-or-self of `a`. -/
def Below (a i : ItemId) : Prop := Relation.ReflTransGen (items.IsParent) a i

/-- The items form a rooted tree with the expected leaves. -/
structure Tree : Prop where
  size : 1 + g.nv + g.ne ≤ items.size
  root : items.type rootItem = .F
  vert : ∀ v, v < g.nv → items.type (vertItem v) = .V
  edge : ∀ e, e < g.ne → items.type (edgeItem g e) = .Q
  node : ∀ i, 1 + g.nv + g.ne ≤ i → i < items.size → items.type i ∉ [NodeType.F, .V, .Q]
  ch_lt : ∀ p c, items.IsParent p c → c < items.size
  unique_parent : ∀ c, 0 < c → c < items.size → ∃ p, items.IsParent p c ∧ ∀ p', items.IsParent p' c → p' = p
  root_no_parent : ∀ p, ¬ items.IsParent p rootItem
  ch_nodup : ∀ p, (items.ch p).Nodup
  reach : ∀ i, i < items.size → items.Below rootItem i
  /-- V items have only Q children; the root has V and Q children only. -/
  v_children : ∀ v c, v < g.nv → items.IsParent (vertItem v) c → items.type c = .Q
  root_children : ∀ c, items.IsParent rootItem c → items.type c = .V ∨ items.type c = .Q

/-- Original edge `e` is below item `i`. -/
def EdgeBelow (i e : ItemId) : Prop := items.Below i (edgeItem g e)

def PairEq (p q : Nat × Nat) : Prop := p = q ∨ p = (q.2, q.1)

/-- Endpoint information: each node's `vs`, the children's `vs`, and the orientation of children
relative to the parent (`vs` of a child are both vertices of the parent's node-vertex list). -/
structure Endpoints : Prop where
  /-- Capped nodes have two endpoints; one-vertex nodes (`O`, loop `Q`) have `(some v, none)`;
  `F` and `V` have none. -/
  vs_shape : ∀ i, i < items.size →
    match items.type i with
    | .F | .V => items.vs i = (none, none)
    | .O => ∃ v, items.vs i = (some v, none)
    | .Q => (∃ v, items.vs i = (some v, none)) ∨ (∃ u v, items.vs i = (some u, some v))
    | _ => ∃ u v, items.vs i = (some u, some v)
  vs_lt : ∀ i u, (items.vs i).1 = some u ∨ (items.vs i).2 = some u → u < g.nv
  /-- A Q item's endpoints are its edge's endpoints. -/
  q_vs : ∀ e, e < g.ne → ∀ u, (items.vs (edgeItem g e)).1 = some u →
    (u = (g.edges[e]!).1 ∨ u = (g.edges[e]!).2) ∧
    ((items.vs (edgeItem g e)).2 = none ↔ (g.edges[e]!).1 = (g.edges[e]!).2) ∧
    ∀ v, (items.vs (edgeItem g e)).2 = some v → PairEq (u, v) g.edges[e]!
  /-- Each child node's endpoints are vertices of the parent node: the parent's own endpoints or
  its V children. -/
  child_vs_in_parent : ∀ p c, items.IsParent p c → items.type p ∉ [NodeType.F, .V] → items.type c ≠ .V →
    ∀ u, (items.vs c).1 = some u ∨ (items.vs c).2 = some u →
      (items.vs p).1 = some u ∨ (items.vs p).2 = some u ∨ items.IsParent p (vertItem u)
  /-- Separation: the subtree of a node touches the rest of the graph only at its endpoints. -/
  separation : ∀ i, i < items.size → items.type i ∉ [NodeType.F, .V] →
    ∀ v e e', v < g.nv → e < g.ne → e' < g.ne →
      ((g.edges[e]!).1 = v ∨ (g.edges[e]!).2 = v) → ((g.edges[e']!).1 = v ∨ (g.edges[e']!).2 = v) →
      items.EdgeBelow g i e → ¬ items.EdgeBelow g i e' →
      (items.vs i).1 = some v ∨ (items.vs i).2 = some v
  /-- Interior: the V children of a node are the vertices whose every edge lies below the node but
  below none of its children. -/
  interior : ∀ i v, i < items.size → v < g.nv →
    (items.IsParent i (vertItem v) ↔
      (∀ e, e < g.ne → ((g.edges[e]!).1 = v ∨ (g.edges[e]!).2 = v) → items.EdgeBelow g i e) ∧
      ∀ c, items.IsParent i c → ¬ ∀ e, e < g.ne → ((g.edges[e]!).1 = v ∨ (g.edges[e]!).2 = v) → items.EdgeBelow g c e)

/-- The non-V children of a node, i.e. its virtual edges, as endpoint pairs. -/
def virtualEdges (i : ItemId) : List (Nat × Nat) :=
  ((items.ch i).filter fun c => items.type c ≠ .V).map fun c => ((items.vs c).1.getD 0, (items.vs c).2.getD 0)

/-- Per-type shape of the children lists. -/
structure Shapes : Prop where
  /-- `I`, `O`: a single leaf under a Q. -/
  i_o_leaf : ∀ i, i < items.size → items.type i = .I ∨ items.type i = .O → items.ch i = []
  /-- `Q`: either a leaf (cap to parent), or a block root with children `[node]` or `[node, v]`. -/
  q_children : ∀ e, e < g.ne → items.ch (edgeItem g e) = [] ∨
    ∃ c, items.type c ∉ [NodeType.F, .V, .Q] ∧
      (items.ch (edgeItem g e) = [c] ∨ ∃ v, items.ch (edgeItem g e) = [c, vertItem v])
  /-- `P`: 2 endpoints, ≥ 2 non-V children all on the same endpoints, no V children. -/
  p_shape : ∀ i, i < items.size → items.type i = .P →
    2 ≤ (items.virtualEdges i).length ∧ (∀ c, items.IsParent i c → items.type c ≠ .V) ∧
    ∀ q ∈ items.virtualEdges i, ∃ u v, items.vs i = (some u, some v) ∧ PairEq q (u, v)
  /-- `S`: with endpoints `u = x₀, x₁, …, x_k = v` (the V children in order) the non-V children
  are exactly the path edges `(x_j, x_{j+1})`, `k ≥ 2`. -/
  s_shape : ∀ i, i < items.size → items.type i = .S → ∃ u v xs,
    items.vs i = (some u, some v) ∧ ((items.ch i).filter fun c => items.type c = .V) = xs.map vertItem ∧
    2 ≤ xs.length ∧
    (items.virtualEdges i).Perm (List.zip (u :: xs) (xs ++ [v]))
  /-- `R`: ≥ 2 V children, simple skeleton with ≥ 6 edges, 3-connected (checked on the output). -/
  r_shape : ∀ i, i < items.size → items.type i = .R →
    2 ≤ ((items.ch i).filter fun c => items.type c = .V).length ∧ 6 ≤ (items.virtualEdges i).length ∧
    ((items.virtualEdges i).map fun q => min q.1 q.2 + g.nv * max q.1 q.2).Nodup ∧
    ∀ q ∈ items.virtualEdges i, q.1 ≠ q.2
  canonical : ∀ p c, items.IsParent p c → (items.type c = .S → items.type p ≠ .S) ∧ (items.type c = .P → items.type p ≠ .P)

structure WF : Prop where
  tree : items.Tree g
  endpoints : items.Endpoints g
  shapes : items.Shapes g

end Items

end Spqr
