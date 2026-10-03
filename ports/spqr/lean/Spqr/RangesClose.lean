import Spqr.RangesFinal
import Spqr.WalkPlace

namespace Spqr
namespace Items

/-- The attachment and shape record of one closed item. -/
structure CloseAt (g : Graph) (items : Items) (i : ItemId) : Prop where
  att_vs : items.type i ∈ [NodeType.S, .P, .R, .Q] →
    ∀ v, items.Att g i v → items.IsVs i v
  vs_att : items.type i ∈ [NodeType.S, .P, .R] ∨ (items.type i = .Q ∧ items.ch i = []) →
    ∀ v, items.IsVs i v → items.Att g i v
  vs_ne : items.type i ∈ [NodeType.S, .P, .R] →
    ∀ u v, items.vs i = (some u, some v) → u ≠ v
  interior : items.type i ∈ [NodeType.S, .P, .R] → ∀ v, v < g.nv →
    (items.IsParent i (vertItem v) ↔
      items.Inner g i v ∧ ∀ c, items.IsParent i c →
        ¬ ∀ e, e < g.ne → g.Inc e v → items.EdgeBelow g c e)
  child_two : ∀ c, items.IsParent i c → items.type i ∈ [NodeType.S, .P, .R] →
    items.type c ≠ .V → ∃ a b, items.vs c = (some a, some b) ∧ a ≠ b
  io_parent : ∀ c, items.IsParent i c → items.type c = .I ∨ items.type c = .O → items.type i = .Q
  q_leaf : ∀ e, i = edgeItem g e → e < g.ne → items.ch i = [] →
    ∃ a b, items.vs i = (some a, some b) ∧ a ≠ b ∧ PairEq (a, b) g.edges[e]!
  q_root : ∀ e, i = edgeItem g e → e < g.ne → items.ch i ≠ [] →
    ∃ u c, items.vs i = (some u, none) ∧ g.Inc e u ∧ items.type c ∉ [NodeType.F, .V] ∧
      (items.type c = .Q → items.ch c = []) ∧
      ((g.edges[e]!).1 = (g.edges[e]!).2 → items.ch i = [c] ∧ items.vs c = (some u, none)) ∧
      ((g.edges[e]!).1 ≠ (g.edges[e]!).2 → ∃ w, w < g.nv ∧ PairEq (u, w) g.edges[e]! ∧
        items.ch i = [c, vertItem w] ∧
        ∃ a b, items.vs c = (some a, some b) ∧ PairEq (a, b) (u, w))
  q_under_v : ∀ v, i = vertItem v → v < g.nv → ∀ c, items.IsParent i c → items.ch c ≠ []
  p_shape : items.type i = .P →
    2 ≤ (items.virtualEdges i).length ∧ (∀ c, items.IsParent i c → items.type c ≠ .V) ∧
    ∀ q ∈ items.virtualEdges i, ∃ u v, items.vs i = (some u, some v) ∧ PairEq q (u, v)
  s_order : items.type i = .S → ∃ u v xs,
    items.vs i = (some u, some v) ∧ ((items.ch i).filter fun c => items.type c = .V) = xs.map vertItem ∧
    1 ≤ xs.length ∧ items.virtualEdges i = List.zip (u :: xs) (xs ++ [v])
  r_shape : items.type i = .R →
    2 ≤ ((items.ch i).filter fun c => items.type c = .V).length ∧ 5 ≤ (items.virtualEdges i).length ∧
    ((items.virtualEdges i).map fun q => min q.1 q.2 + g.nv * max q.1 q.2).Nodup ∧
    ∀ u v, items.vs i = (some u, some v) → ∀ q ∈ items.virtualEdges i, ¬ PairEq q (u, v)

theorem CloseAt.root {g : Graph} {items : Items} (ht : items.type rootItem = .F)
    (hn : items.ch rootItem = []) : CloseAt g items rootItem := by
  constructor
  · simp [ht]
  · simp [ht]
  · simp [ht]
  · simp [ht]
  · simp [ht]
  · simp [hn, IsParent]
  · intro e he; have he' : 0 = 1 + g.nv + e := he; omega
  · intro e he; have he' : 0 = 1 + g.nv + e := he; omega
  · intro v hv; have hv' : 0 = 1 + v := hv; omega
  · simp [ht]
  · simp [ht]
  · simp [ht]

theorem parent_lt {items : Items} {p c : ItemId} (h : items.IsParent p c) : p < items.size := by
  by_contra hn
  simp [IsParent, ch_of_le _ _ (Nat.le_of_not_lt hn)] at h

end Items

namespace WalkState
open WalkM

/-- Closed records are required for the root and for items occurring in a span or child list.
Freshly allocated and reopened items are loose until their next close. -/
structure CloseInv (s : WalkState) : Prop where
  closed : ∀ i, i < s.items.size → i = rootItem ∨ 0 < s.cnt i → Items.CloseAt s.g s.items i

theorem init_closeInv (g : Graph) (tern : Bool) : (init g tern).CloseInv := by
  constructor
  intro i hi hc
  rcases hc with rfl | hc
  · exact Items.CloseAt.root (by simp [init, Items.initialItems_type, rootItem])
      (by simp [init, Items.initialItems_ch])
  · rw [init_cnt] at hc; omega

theorem CloseInv.frame {s s' : WalkState} (h : s.CloseInv) (hg : s'.g = s.g)
    (hi : s'.items = s.items) (hc : ∀ i, 0 < s'.cnt i → 0 < s.cnt i) : s'.CloseInv := by
  constructor
  intro i hil hl
  rw [hg, hi]
  exact h.closed i (by rwa [hi] at hil) (hl.imp_right (hc i))

theorem CloseInv.of_tree {s : WalkState} (h : s.CloseInv) (ht : Items.Tree s.g s.items)
    (hty : WalkTyping s.g s.items) : Items.CloseFacts s.g s.items := by
  have hall : ∀ i, i < s.items.size → Items.CloseAt s.g s.items i := by
    intro i hi
    apply h.closed i hi
    by_cases hz : i = rootItem
    · exact Or.inl hz
    · obtain ⟨p, hp, -⟩ := ht.unique_parent i (Nat.pos_of_ne_zero hz) hi
      have hpos := List.count_pos_iff.mpr hp
      have hle := Items.count_le_chCount s.items (Items.parent_lt hp) i
      exact Or.inr (by dsimp [cnt]; omega)
  have he : ∀ e, e < s.g.ne → edgeItem s.g e < s.items.size := by
    intro e he; have := hty.size; change 1 + s.g.nv + e < s.items.size; omega
  refine ⟨?_, fun i hi => (hall i hi).vs_att, fun i hi => (hall i hi).vs_ne,
    fun i v hi hv ht => (hall i hi).interior ht v hv,
    fun p c hp => (hall p (Items.parent_lt hp)).child_two c hp,
    fun p c hp => (hall p (Items.parent_lt hp)).io_parent c hp,
    fun e he' => (hall _ (he e he')).q_leaf e rfl he',
    fun e he' => (hall _ (he e he')).q_root e rfl he', ?_,
    fun i hi => (hall i hi).p_shape, fun i hi => (hall i hi).s_order,
    fun i hi => (hall i hi).r_shape⟩
  · intro e he' v
    exact (hall _ (he e he')).att_vs (by rw [hty.edge e he']; simp) v
  · intro v c hv hp
    exact (hall _ (Items.parent_lt hp)).q_under_v v rfl hv c hp

end WalkState
end Spqr
