import Spqr.ItemSpec
import Spqr.GraphLemmas
import Spqr.SepPair

/-!
# Postorder ranges and attachments of the walk's items

`Items.Ranges g σ items` is a final-state property of the walk's items (PROOF.md §4.6): every node
item owns a contiguous interval of the edge order `σ` (R1), and the attachment vertices of the
edge set below a node are exactly its recorded endpoints `vs`, together with the few per-type
attachment facts about children (R2). `RangesWF.lean` derives `Items.Endpoints` and `Items.Shapes`
from it by pure item-level reasoning, with no reference to the tstack.
-/

namespace Spqr
namespace Items

variable (g : Graph) (items : Items)

/-- `v` is one of the endpoints recorded in `vs i`. -/
def IsVs (i v : Nat) : Prop := (items.vs i).1 = some v ∨ (items.vs i).2 = some v

/-- `v` is an attachment vertex of the edge set below `i`: it has an edge below `i` and one not
below `i`. -/
def Att (i v : Nat) : Prop :=
  ∃ e e', e < g.ne ∧ e' < g.ne ∧ g.Inc e v ∧ g.Inc e' v ∧ items.EdgeBelow g i e ∧ ¬ items.EdgeBelow g i e'

/-- `v` is a non-isolated vertex all of whose edges lie below `i`. -/
def Inner (i v : Nat) : Prop :=
  (∃ e, e < g.ne ∧ g.Inc e v) ∧ ∀ e, e < g.ne → g.Inc e v → items.EdgeBelow g i e

/-- `i` reaches `j` by parent steps that never enter a V item: `j` belongs to `i`'s own block. -/
def BelowNoV (a i : ItemId) : Prop :=
  Relation.ReflTransGen (fun p c => items.IsParent p c ∧ items.type c ≠ .V) a i

/-- Original edge `e` is a piece edge of `i`: below `i` within `i`'s block (not in a block hanging
from one of its interior vertices). -/
def PieceEdge (i e : Nat) : Prop := items.BelowNoV i (edgeItem g e)

/-- Range/attachment structure of the final items relative to the edge order `σ` (for the walk,
`edgePostorderForest forest`: each edge at its deeper DFS endpoint, in `walkOut` order). -/
structure Ranges (σ : List Nat) : Prop where
  /-- R1: the piece edges of a node item form a contiguous interval of `σ` once the blocks hanging
  from its interior vertices (processed in the middle of it) are skipped; the blocks hanging from
  interior vertices need not be contiguous with it (PROOF.md §4.6). -/
  convex : ∀ i, i < items.size → items.type i ∉ [NodeType.F, .V] →
    ∀ a b c, a ≤ b → b ≤ c → c < σ.length →
      items.PieceEdge g i σ[a]! → items.PieceEdge g i σ[c]! → items.EdgeBelow g i σ[b]!
  /-- R2: the attachments of a node item are among its endpoints … -/
  att_vs : ∀ i, i < items.size → items.type i ∉ [NodeType.F, .V] →
    ∀ v, items.Att g i v → items.IsVs i v
  /-- … and at an S/P/R node or a leaf Q both endpoints are attachments. -/
  vs_att : ∀ i, i < items.size → items.type i ∈ [NodeType.S, .P, .R] ∨ (items.type i = .Q ∧ items.ch i = []) →
    ∀ v, items.IsVs i v → items.Att g i v
  /-- An S/P/R node's two endpoints are distinct. -/
  vs_ne : ∀ i, i < items.size → items.type i ∈ [NodeType.S, .P, .R] →
    ∀ u v, items.vs i = (some u, some v) → u ≠ v
  /-- The V children of an S/P/R node are its interior vertices that are interior to no child. -/
  interior : ∀ i v, i < items.size → v < g.nv → items.type i ∈ [NodeType.S, .P, .R] →
    (items.IsParent i (vertItem v) ↔
      items.Inner g i v ∧ ∀ c, items.IsParent i c → ¬ ∀ e, e < g.ne → g.Inc e v → items.EdgeBelow g c e)
  /-- A non-V child of an S/P/R node has two distinct endpoints. -/
  child_two : ∀ p c, items.IsParent p c → items.type p ∈ [NodeType.S, .P, .R] → items.type c ≠ .V →
    ∃ a b, items.vs c = (some a, some b) ∧ a ≠ b
  /-- `I`/`O` items hang under a Q. -/
  io_parent : ∀ p c, items.IsParent p c → items.type c = .I ∨ items.type c = .O → items.type p = .Q
  /-- A leaf Q is the cap of a non-loop edge. -/
  q_leaf : ∀ e, e < g.ne → items.ch (edgeItem g e) = [] →
    ∃ a b, items.vs (edgeItem g e) = (some a, some b) ∧ a ≠ b ∧ PairEq (a, b) g.edges[e]!
  /-- A block-root Q records its upper endpoint `u`; a self-loop has one child `c` with `vs c =
  (some u, none)`, otherwise the children are `[c, vertItem w]` with `{u, w}` the edge and `vs c`
  the edge's endpoints (`c` is a leaf Q for a digon). -/
  q_root : ∀ e, e < g.ne → items.ch (edgeItem g e) ≠ [] →
    ∃ u c, items.vs (edgeItem g e) = (some u, none) ∧ g.Inc e u ∧ items.type c ∉ [NodeType.F, .V] ∧
      (items.type c = .Q → items.ch c = []) ∧
      ((g.edges[e]!).1 = (g.edges[e]!).2 →
        items.ch (edgeItem g e) = [c] ∧ items.vs c = (some u, none)) ∧
      ((g.edges[e]!).1 ≠ (g.edges[e]!).2 → ∃ w, w < g.nv ∧ PairEq (u, w) g.edges[e]! ∧
        items.ch (edgeItem g e) = [c, vertItem w] ∧
        ∃ a b, items.vs c = (some a, some b) ∧ PairEq (a, b) (u, w))
  /-- A child of a V item is a block root (leaf Qs hang only under nodes and block-root Qs). -/
  q_under_v : ∀ v c, v < g.nv → items.IsParent (vertItem v) c → items.ch c ≠ []
  /-- `P`: ≥ 2 virtual edges, all on the node's endpoints, no V children. -/
  p_shape : ∀ i, i < items.size → items.type i = .P →
    2 ≤ (items.virtualEdges i).length ∧ (∀ c, items.IsParent i c → items.type c ≠ .V) ∧
    ∀ q ∈ items.virtualEdges i, ∃ u v, items.vs i = (some u, some v) ∧ PairEq q (u, v)
  /-- `S`: the non-V children in `ch` order are the path edges `(x_j, x_{j+1})` from `u` through
  the V children (in `ch` order, at least one) to `v`. -/
  s_order : ∀ i, i < items.size → items.type i = .S → ∃ u v xs,
    items.vs i = (some u, some v) ∧ ((items.ch i).filter fun c => items.type c = .V) = xs.map vertItem ∧
    1 ≤ xs.length ∧ items.virtualEdges i = List.zip (u :: xs) (xs ++ [v])
  /-- `R`: ≥ 2 V children and ≥ 5 pairwise non-parallel virtual edges, none parallel to `vs`. -/
  r_shape : ∀ i, i < items.size → items.type i = .R →
    2 ≤ ((items.ch i).filter fun c => items.type c = .V).length ∧ 5 ≤ (items.virtualEdges i).length ∧
    ((items.virtualEdges i).map fun q => min q.1 q.2 + g.nv * max q.1 q.2).Nodup ∧
    ∀ u v, items.vs i = (some u, some v) → ∀ q ∈ items.virtualEdges i, ¬ PairEq q (u, v)

end Items
end Spqr
