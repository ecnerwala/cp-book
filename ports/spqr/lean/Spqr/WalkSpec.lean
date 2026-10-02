import Spqr.ItemSpec

/-!
# The walk invariant: soundness and completeness per step

Stated on `walkTree` directly (no ear structure). Every tstack entry `t` describes a connected set
of original edges `t.edges` whose attachment to the rest of the block is exactly the two terminals
`t.vStart` and `stackVerts[t.topDepth]`. *Soundness*: `mergeTstackTops` joins two entries sharing a
terminal, so the union is again connected and 2-attached. *Completeness*: every item `finishEdge`
closes is 2-attached and connected, i.e. cut off by a genuine 2-separation, and the closed items
partition the edges. The stack-shape facts these steps rely on (`size ≥ origTstack + 3`, `nxt`
being the chain's entry) are proved on the ear layer (`EarSpec.lean`).
-/

namespace Spqr

/-- `u` and `w` are joined by an edge of `E`. -/
def Graph.AdjIn (g : Graph) (E : Nat → Prop) (u w : Nat) : Prop :=
  ∃ e, e < g.ne ∧ E e ∧ Items.PairEq (u, w) g.edges[e]!

/-- The edge set `E` is connected (as a subgraph; isolated vertices are not part of it). -/
def Graph.ConnEdges (g : Graph) (E : Nat → Prop) : Prop :=
  ∀ e e', e < g.ne → e' < g.ne → E e → E e' →
    Relation.ReflTransGen (g.AdjIn E) (g.edges[e]!).1 (g.edges[e']!).1

/-- Every vertex incident to both an edge in `E` and an edge outside `E` is `a` or `b`. -/
def Graph.TwoAttached (g : Graph) (E : Nat → Prop) (a b : Nat) : Prop :=
  ∀ v e e', e < g.ne → e' < g.ne → E e → ¬ E e' →
    ((g.edges[e]!).1 = v ∨ (g.edges[e]!).2 = v) → ((g.edges[e']!).1 = v ∨ (g.edges[e']!).2 = v) →
    v = a ∨ v = b

namespace TEntry

/-- The original edges covered by an entry: everything below either side's items. -/
def edges (g : Graph) (items : Items) (t : TEntry) (e : Nat) : Prop :=
  ∃ i ∈ t.spans.1 ++ t.spans.2, items.EdgeBelow g i e

/-- Second terminal, read off the vertex stack. -/
def top (s : WalkState) (t : TEntry) : Nat := s.stackVerts[t.topDepth]!

end TEntry

namespace WalkState

/-- Invariant on one tstack entry. -/
structure EntryInv (s : WalkState) (t : TEntry) : Prop where
  conn : s.g.ConnEdges (t.edges s.g s.items)
  attached : s.g.TwoAttached (t.edges s.g s.items) t.vStart (t.top s)

/-- Invariant on a closed item (its subtree is a 2-attached connected piece). -/
structure ItemInv (s : WalkState) (i : ItemId) : Prop where
  conn : s.g.ConnEdges (Items.EdgeBelow s.g s.items i)
  attached : ∀ u v, Items.vs s.items i = (some u, some v) →
    s.g.TwoAttached (Items.EdgeBelow s.g s.items i) u v

/-- The walk invariant: every open entry and every allocated node is a 2-attached connected piece. -/
structure Inv (s : WalkState) : Prop where
  entries : ∀ t ∈ s.tstack, s.EntryInv t
  nodes : ∀ i, 1 + s.g.nv + s.g.ne ≤ i → i < s.items.size → s.ItemInv i

end WalkState

open WalkM

/-- Soundness of a merge: two adjacent entries sharing the terminal `cur.top = nxt.vStart`
(or `cur.vStart = nxt.vStart`, the P case) merge into a 2-attached connected entry. -/
theorem mergeTstackTops_sound (s : WalkState) (a b : TEntry) (rest : List TEntry)
    (hs : s.tstack = a :: b :: rest) (ha : s.EntryInv a) (hb : s.EntryInv b)
    (hshare : a.top s = b.vStart ∨ a.vStart = b.vStart ∨ a.top s = b.top s) :
    let s' := (mergeTstackTops.run s).2
    ∀ t ∈ s'.tstack, s'.EntryInv t := by
  sorry

/-- Closing the top entry into `item` records a 2-attached connected piece. -/
theorem finishTstackTop_complete (s : WalkState) (item : ItemId) (h : s.Inv)
    (hcur : ∀ t, s.tstack.head? = some t → s.EntryInv t) :
    ((finishTstackTop item).run s).2.Inv := by
  sorry

/-- `finishEdge` preserves the invariant. -/
theorem finishEdge_inv (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool)
    (s : WalkState) (h : s.Inv) :
    ((finishEdge curV d o origTstack hasVert).run s).2.Inv := by
  sorry

/-- The whole walk preserves the invariant. -/
theorem walkTree_inv (t : DfsTree) (d : Nat) (s : WalkState) (h : s.Inv) :
    ((walkTree t d).run s).2.Inv := by
  sorry

/-- Completeness: after the walk, the allocated nodes partition the edges (every edge is below
exactly one child chain from the root) — the `Items.Tree` content of `Items.WF`. -/
theorem walk_nodes_partition (g : Graph) (tern : Bool) (forest : List DfsTree) :
    let s := g.walk tern forest
    ∀ e, e < g.ne → ∃ p, Items.IsParent s.items p (edgeItem g e) ∧
      ∀ p', Items.IsParent s.items p' (edgeItem g e) → p' = p := by
  sorry

end Spqr
