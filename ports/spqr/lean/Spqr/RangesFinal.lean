import Spqr.RangesInv

/-!
# From the final range invariant to `Items.Ranges`

`ranges_of_rangesInv`: `Items.Ranges` for a walk state satisfying `RangesInv` (any `n`, `D`),
`WalkTyping`, and `Items.CloseFacts` — the clauses of `Items.Ranges` that `RangesInv` does not
carry. From `RangesInv` come `convex` (= `closed`) and `att_vs` for the allocated node items
(`Inv'.nodes`: `TwoAttached` at their `vs`; `I`/`O` leaves own no edge). Everything else —
`att_vs` for the `Q` items and all the `vs_att`/shape/placement clauses — is about the state of an
item at its close (its `vs` were just written from the entries' terminals, its children's shape),
which `RangesInv` deliberately does not record; `CloseFacts` names exactly that remainder.
-/

namespace Spqr
namespace Items
variable {g : Graph} {items : Items}

theorem Below_of_ch_nil {i j : ItemId} (h : items.Below i j) (hnil : items.ch i = []) : j = i := by
  rcases Relation.ReflTransGen.cases_head h with rfl | ⟨c, hc, -⟩
  · rfl
  · simp [IsParent, hnil] at hc

/-- The clauses of `Items.Ranges` not derivable from the final `RangesInv` (+ `WalkTyping`):
`att_vs` restricted to the `Q` items, and the remaining attachment/shape clauses verbatim. -/
structure CloseFacts (g : Graph) (items : Items) : Prop where
  /-- `att_vs` for the `Q` items (`RangesInv.inv` only covers the allocated node items). -/
  q_att_vs : ∀ e, e < g.ne → ∀ v, items.Att g (edgeItem g e) v → items.IsVs (edgeItem g e) v
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

namespace WalkState
variable {s : WalkState} {σ : List Nat} {n D : Nat}

theorem ranges_of_rangesInv (h : s.RangesInv σ n D) (ht : WalkTyping s.g s.items)
    (hc : Items.CloseFacts s.g s.items) : Items.Ranges s.g s.items σ where
  convex := h.closed
  att_vs := fun i hi hty v hatt => by
    by_cases hlo : 1 + s.g.nv + s.g.ne ≤ i
    · have hnode := ht.node i hlo hi
      have hvs := ht.vs_shape i hi
      obtain ⟨e, e', he, he', hev, he'v, hb, hnb⟩ := hatt
      have two : ∀ u w, Items.vs s.items i = (some u, some w) → Items.IsVs s.items i v := by
        intro u w huw
        rcases (h.inv.nodes i hlo hi).attached u w huw v e e' he he' hb hnb hev he'v with rfl | rfl
        · exact Or.inl (by rw [huw])
        · exact Or.inr (by rw [huw])
      cases hT : Items.type s.items i <;> rw [hT] at hvs hnode hty
      all_goals first
        | (simp at hty; done)
        | (simp at hnode; done)
        | (obtain ⟨u, w, huw⟩ := hvs; exact two u w huw)
        | (have hnil := ht.i_o_leaf i hi (Or.inr hT)
           have h1 : 1 + s.g.nv + e = i := Items.Below_of_ch_nil hb hnil
           omega)
    · by_cases h0 : i = 0
      · subst h0
        have hr := ht.root
        unfold rootItem at hr
        rw [hr] at hty; simp at hty
      · have hi1 : 1 ≤ i := Nat.pos_of_ne_zero h0
        by_cases hv : i < 1 + s.g.nv
        · have hvi : vertItem (i - 1) = i := Nat.add_sub_cancel' hi1
          have hlt : i - 1 < s.g.nv := by omega
          rw [← hvi, ht.vert (i - 1) hlt] at hty; simp at hty
        · have hle : 1 + s.g.nv ≤ i := Nat.le_of_not_lt hv
          have hei : edgeItem s.g (i - (1 + s.g.nv)) = i := Nat.add_sub_cancel' hle
          have hlt : i - (1 + s.g.nv) < s.g.ne := by omega
          rw [← hei] at hatt ⊢
          exact hc.q_att_vs _ hlt v hatt
  vs_att := hc.vs_att
  vs_ne := hc.vs_ne
  interior := hc.interior
  child_two := hc.child_two
  io_parent := hc.io_parent
  q_leaf := hc.q_leaf
  q_root := hc.q_root
  q_under_v := hc.q_under_v
  p_shape := hc.p_shape
  s_order := hc.s_order
  r_shape := hc.r_shape

end WalkState
end Spqr
