import Spqr.PlanarGlue
import Spqr.Proofs.Planar
import Mathlib.Data.Finset.Card

/-!
# Vertex and component counts of a 2-sum (PROOF.md §8.5)

The computable counters `numNonIsolated` / `numComponents` of `Planar.lean` on the edge list of a
2-sum `TwoSum.edges`, in terms of the two sides: the non-isolated vertices are those of `G₁` plus
those of `G₂` other than its two terminals (identified with `u₁, v₁`).
-/

namespace Spqr

open Classical

/-- `v` has an incident edge in `es` (the proposition behind `nonIsolated`). -/
def HasEdge (es : List (Nat × Nat)) (v : Nat) : Prop := ∃ p ∈ es, p.1 = v ∨ p.2 = v

theorem nonIsolated_iff (es : List (Nat × Nat)) (v : Nat) :
    nonIsolated es v = true ↔ HasEdge es v := by
  simp [nonIsolated, HasEdge, List.any_eq_true]

theorem numNonIsolated_eq_card (es : List (Nat × Nat)) (n : Nat) :
    numNonIsolated es n = ((Finset.range n).filter (HasEdge es)).card := by
  unfold numNonIsolated
  rw [length_filter_range_eq_card]
  congr 1
  exact Finset.filter_congr fun v _ => nonIsolated_iff es v

/-- Two distinct connected vertices both have an edge. -/
theorem HasEdge.of_edgesConn {es : List (Nat × Nat)} {a b : Nat} (h : EdgesConn es a b)
    (hab : a ≠ b) : HasEdge es a ∧ HasEdge es b := by
  constructor
  · rcases h.cases_head with rfl | ⟨c, hac, -⟩
    · exact (hab rfl).elim
    · rcases hac with h | h
      · exact ⟨(a, c), h, Or.inl rfl⟩
      · exact ⟨(c, a), h, Or.inr rfl⟩
  · rcases h.cases_tail with rfl | ⟨c, -, hcb⟩
    · exact (hab rfl).elim
    · rcases hcb with h | h
      · exact ⟨(c, b), h, Or.inr rfl⟩
      · exact ⟨(b, c), h, Or.inl rfl⟩

namespace TwoSum

variable (T : TwoSum) {rs₁ rs₂ : RotationSystem} (W : T.WF rs₁ rs₂)
include W

theorem mem_e₁ : (T.u₁, T.v₁) ∈ T.es₁ := List.mem_iff_getElem?.2 ⟨_, W.e₁⟩
theorem mem_e₂ : (T.u₂, T.v₂) ∈ T.es₂ := List.mem_iff_getElem?.2 ⟨_, W.e₂⟩
theorem u₁_lt' : T.u₁ < T.n₁ := (W.emb₁.verts _ (T.mem_e₁ W)).1
theorem v₁_lt' : T.v₁ < T.n₁ := (W.emb₁.verts _ (T.mem_e₁ W)).2
theorem u₂_lt' : T.u₂ < T.n₂ := (W.emb₂.verts _ (T.mem_e₂ W)).1
theorem v₂_lt' : T.v₂ < T.n₂ := (W.emb₂.verts _ (T.mem_e₂ W)).2

omit W in
theorem vert₂_u' : T.vert₂ T.u₂ = T.u₁ := by simp [vert₂]
theorem vert₂_v' : T.vert₂ T.v₂ = T.v₁ := by simp [vert₂, W.uv₂.symm]

omit W in
theorem vert₂_cases (x : Nat) : T.vert₂ x = T.u₁ ∨ T.vert₂ x = T.v₁ ∨ T.n₁ ≤ T.vert₂ x := by
  unfold vert₂; split_ifs <;> omega

omit W in
theorem vert₂_eq_add (hu : T.u₁ < T.n₁) (hv : T.v₁ < T.n₁) {w x : Nat}
    (h : T.vert₂ x = T.n₁ + w) : x = w ∧ x ≠ T.u₂ ∧ x ≠ T.v₂ := by
  unfold vert₂ at h; split_ifs at h <;> omega

omit W in
theorem hasEdge_edges_of_left {v : Nat} (h : HasEdge (T.es₁.eraseIdx T.e₁) v) :
    HasEdge T.edges v := by
  obtain ⟨p, hp, hpv⟩ := h
  unfold edges
  exact ⟨p, List.mem_append_left _ hp, hpv⟩

omit W in
theorem hasEdge_edges_of_right {w : Nat} (h : HasEdge (T.es₂.eraseIdx T.e₂) w) :
    HasEdge T.edges (T.vert₂ w) := by
  obtain ⟨p, hp, hpv⟩ := h
  unfold edges
  refine ⟨(T.vert₂ p.1, T.vert₂ p.2), List.mem_append_right _ (List.mem_map.2 ⟨p, hp, rfl⟩), ?_⟩
  rcases hpv with h | h
  · exact Or.inl (by rw [h])
  · exact Or.inr (by rw [h])

/-- A vertex of `G₁` has an edge in the 2-sum iff it has one in `G₁` (the terminals keep an edge
by `conn`). -/
theorem hasEdge_edges_lt {v : Nat} (hv : v < T.n₁) : HasEdge T.edges v ↔ HasEdge T.es₁ v := by
  have hu := T.u₁_lt' W; have hv₁ := T.v₁_lt' W
  constructor
  · rintro ⟨p, hp, hpv⟩
    unfold edges at hp
    rw [List.mem_append] at hp
    rcases hp with hp | hp
    · exact ⟨p, List.mem_of_mem_eraseIdx hp, hpv⟩
    · rw [List.mem_map] at hp
      obtain ⟨q, -, rfl⟩ := hp
      simp only at hpv
      have h3 : v = T.u₁ ∨ v = T.v₁ := by
        rcases hpv with h | h
        · have := T.vert₂_cases q.1; omega
        · have := T.vert₂_cases q.2; omega
      rcases h3 with rfl | rfl
      · exact ⟨_, T.mem_e₁ W, Or.inl rfl⟩
      · exact ⟨_, T.mem_e₁ W, Or.inr rfl⟩
  · rintro ⟨p, hp, hpv⟩
    obtain ⟨i, hi⟩ := List.mem_iff_getElem?.1 hp
    by_cases hie : i = T.e₁
    · subst hie
      rw [W.e₁] at hi
      obtain rfl := Option.some.inj hi
      simp only at hpv
      rcases W.conn with hc | hc
      · have h2 := HasEdge.of_edgesConn hc W.uv₁
        rcases hpv with rfl | rfl
        · exact T.hasEdge_edges_of_left h2.1
        · exact T.hasEdge_edges_of_left h2.2
      · have h2 := HasEdge.of_edgesConn hc W.uv₂
        rcases hpv with rfl | rfl
        · have := T.hasEdge_edges_of_right h2.1; rwa [T.vert₂_u'] at this
        · have := T.hasEdge_edges_of_right h2.2; rwa [T.vert₂_v' W] at this
    · exact T.hasEdge_edges_of_left ⟨p, List.mem_eraseIdx_iff_getElem?.2 ⟨i, hie, hi⟩, hpv⟩

/-- A shifted vertex `n₁ + w` of `G₂` has an edge in the 2-sum iff `w` has one in `G₂` and is not
a terminal. -/
theorem hasEdge_edges_add {w : Nat} (hw : w < T.n₂) :
    HasEdge T.edges (T.n₁ + w) ↔ HasEdge T.es₂ w ∧ w ≠ T.u₂ ∧ w ≠ T.v₂ := by
  have hu := T.u₁_lt' W; have hv₁ := T.v₁_lt' W
  constructor
  · rintro ⟨p, hp, hpv⟩
    unfold edges at hp
    rw [List.mem_append] at hp
    rcases hp with hp | hp
    · have := W.emb₁.verts p (List.mem_of_mem_eraseIdx hp); omega
    · rw [List.mem_map] at hp
      obtain ⟨q, hq, rfl⟩ := hp
      simp only at hpv
      rcases hpv with h | h
      · obtain ⟨rfl, h1, h2⟩ := T.vert₂_eq_add hu hv₁ h
        exact ⟨⟨q, List.mem_of_mem_eraseIdx hq, Or.inl rfl⟩, h1, h2⟩
      · obtain ⟨rfl, h1, h2⟩ := T.vert₂_eq_add hu hv₁ h
        exact ⟨⟨q, List.mem_of_mem_eraseIdx hq, Or.inr rfl⟩, h1, h2⟩
  · rintro ⟨⟨p, hp, hpv⟩, hu2, hv2⟩
    obtain ⟨i, hi⟩ := List.mem_iff_getElem?.1 hp
    have hie : i ≠ T.e₂ := by
      rintro rfl
      rw [W.e₂] at hi
      obtain rfl := Option.some.inj hi
      simp only at hpv
      omega
    have h := T.hasEdge_edges_of_right ⟨p, List.mem_eraseIdx_iff_getElem?.2 ⟨i, hie, hi⟩, hpv⟩
    have hw' : T.vert₂ w = T.n₁ + w := by simp [vert₂, hu2, hv2]
    rwa [hw'] at h

/-- Non-isolated vertices of the 2-sum: those of `G₁` plus those of `G₂` other than its two
terminals. -/
theorem numNonIsolated_edges :
    numNonIsolated T.edges T.nVerts + 2 = numNonIsolated T.es₁ T.n₁ + numNonIsolated T.es₂ T.n₂ := by
  have hu := T.u₂_lt' W; have hv := T.v₂_lt' W
  rw [numNonIsolated_eq_card, numNonIsolated_eq_card, numNonIsolated_eq_card]
  set A₂ := (Finset.range T.n₂).filter (HasEdge T.es₂) with hA₂
  have hsub : ({T.u₂, T.v₂} : Finset Nat) ⊆ A₂ := by
    intro x hx
    rw [Finset.mem_insert, Finset.mem_singleton] at hx
    rw [hA₂, Finset.mem_filter, Finset.mem_range]
    rcases hx with rfl | rfl
    · exact ⟨hu, _, T.mem_e₂ W, Or.inl rfl⟩
    · exact ⟨hv, _, T.mem_e₂ W, Or.inr rfl⟩
  have hsplit : (Finset.range T.nVerts).filter (HasEdge T.edges) =
      (Finset.range T.n₁).filter (HasEdge T.es₁) ∪ (A₂ \ {T.u₂, T.v₂}).image (T.n₁ + ·) := by
    ext v
    simp only [Finset.mem_filter, Finset.mem_range, Finset.mem_union, Finset.mem_image,
      Finset.mem_sdiff, Finset.mem_insert, Finset.mem_singleton, hA₂]
    constructor
    · rintro ⟨hv, he⟩
      by_cases h1 : v < T.n₁
      · exact Or.inl ⟨h1, (T.hasEdge_edges_lt W h1).1 he⟩
      · right
        have hw : v - T.n₁ < T.n₂ := by unfold nVerts at hv; omega
        have hv' : v = T.n₁ + (v - T.n₁) := by omega
        rw [hv'] at he
        obtain ⟨h2, h3, h4⟩ := (T.hasEdge_edges_add W hw).1 he
        exact ⟨v - T.n₁, ⟨⟨hw, h2⟩, by omega⟩, by omega⟩
    · rintro (⟨h1, he⟩ | ⟨w, ⟨⟨hw, he⟩, hne⟩, rfl⟩)
      · exact ⟨by unfold nVerts; omega, (T.hasEdge_edges_lt W h1).2 he⟩
      · exact ⟨by unfold nVerts; omega,
          (T.hasEdge_edges_add W hw).2 ⟨he, fun h => hne (Or.inl h), fun h => hne (Or.inr h)⟩⟩
  rw [hsplit, Finset.card_union_of_disjoint, Finset.card_image_of_injective _ (fun a b h => Nat.add_left_cancel h),
    ← Finset.card_sdiff_add_card_eq_card hsub, Finset.card_pair W.uv₂]
  · omega
  · rw [Finset.disjoint_left]
    intro v hv hv'
    rw [Finset.mem_filter, Finset.mem_range] at hv
    rw [Finset.mem_image] at hv'
    obtain ⟨w, -, rfl⟩ := hv'
    omega

end TwoSum

end Spqr
