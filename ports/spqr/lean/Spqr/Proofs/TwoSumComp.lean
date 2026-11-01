import Spqr.Proofs.CompGlue

/-!
# Components of a 2-sum

`numComponents_edges`: `C + 1 = C₁ + C₂`. The 2-sum with both virtual edges kept is the disjoint
union `G₁ ⊔ (G₂ + n₁)` with `n₁ + u₂ ↦ u₁` and `n₁ + v₂ ↦ v₁` identified; the first
identification merges two components, the second does not (`u₁ ~ v₁` through the kept virtual
edge), and dropping the virtual edges changes nothing because `conn` keeps `u₁ ~ v₁`.
-/

namespace Spqr

namespace TwoSum

variable (T : TwoSum)

/-- Vertex map from the disjoint union `G₁ ⊔ (G₂ + n₁)` onto the 2-sum. -/
def vmap (x : Nat) : Nat := ident T.v₁ (T.n₁ + T.v₂) (ident T.u₁ (T.n₁ + T.u₂) x)

/-- The 2-sum with both virtual edges kept. -/
def glued : List (Nat × Nat) :=
  identEdges T.v₁ (T.n₁ + T.v₂)
    (identEdges T.u₁ (T.n₁ + T.u₂) (T.es₁ ++ shiftEdges T.n₁ T.es₂))

theorem mem_glued {p : Nat × Nat} :
    p ∈ T.glued ↔ ∃ q ∈ T.es₁ ++ shiftEdges T.n₁ T.es₂, p = (T.vmap q.1, T.vmap q.2) := by
  unfold glued
  constructor
  · intro h
    obtain ⟨q, hq, rfl⟩ := mem_identEdges.1 h
    obtain ⟨r, hr, rfl⟩ := mem_identEdges.1 hq
    exact ⟨r, hr, rfl⟩
  · rintro ⟨q, hq, rfl⟩
    exact mem_identEdges.2 ⟨_, mem_identEdges.2 ⟨q, hq, rfl⟩, rfl⟩

variable {rs₁ rs₂ : RotationSystem} (W : T.WF rs₁ rs₂)
include W

theorem vmap_lt {x : Nat} (hx : x < T.n₁) : T.vmap x = x := by
  have := T.u₁_lt' W; have := T.v₁_lt' W
  unfold vmap ident; split_ifs <;> omega

theorem vmap_add (w : Nat) : T.vmap (T.n₁ + w) = T.vert₂ w := by
  have := T.u₁_lt' W; have := T.v₁_lt' W; have := W.uv₂
  unfold vmap ident vert₂; split_ifs <;> omega

theorem vert₂_lt {w : Nat} (hw : w < T.n₂) : T.vert₂ w < T.n₁ + T.n₂ := by
  have := T.u₁_lt' W; have := T.v₁_lt' W
  unfold vert₂; split_ifs <;> omega

theorem edges_verts' : ∀ p ∈ T.edges, p.1 < T.n₁ + T.n₂ ∧ p.2 < T.n₁ + T.n₂ := by
  intro p hp
  unfold edges at hp
  rw [List.mem_append, List.mem_map] at hp
  rcases hp with hp | ⟨q, hq, rfl⟩
  · have := W.emb₁.verts p (List.mem_of_mem_eraseIdx hp)
    omega
  · have := W.emb₂.verts q (List.mem_of_mem_eraseIdx hq)
    exact ⟨T.vert₂_lt W this.1, T.vert₂_lt W this.2⟩

theorem edgesConn_edges_uv : EdgesConn T.edges T.u₁ T.v₁ := by
  have hu : T.vert₂ T.u₂ = T.u₁ := by simp [vert₂]
  have hv : T.vert₂ T.v₂ = T.v₁ := by simp [vert₂, W.uv₂.symm]
  rcases W.conn with h | h
  · exact edgesConn_mono (es' := T.edges)
      (fun p hp => by unfold edges; exact List.mem_append_left _ hp) h
  · have := edgesConn_lift (es' := T.edges) T.vert₂ (fun p hp => by
      unfold edges; exact List.mem_append_right _ (List.mem_map.2 ⟨p, hp, rfl⟩)) h
    rwa [hu, hv] at this

theorem mem_es₁_cases {q : Nat × Nat} (hq : q ∈ T.es₁) :
    q ∈ T.es₁.eraseIdx T.e₁ ∨ q = (T.u₁, T.v₁) := by
  obtain ⟨i, hi⟩ := List.mem_iff_getElem?.1 hq
  by_cases hie : i = T.e₁
  · subst hie; right; rw [W.e₁] at hi; exact (Option.some.inj hi).symm
  · left; exact List.mem_eraseIdx_iff_getElem?.2 ⟨i, hie, hi⟩

theorem mem_es₂_cases {q : Nat × Nat} (hq : q ∈ T.es₂) :
    q ∈ T.es₂.eraseIdx T.e₂ ∨ q = (T.u₂, T.v₂) := by
  obtain ⟨i, hi⟩ := List.mem_iff_getElem?.1 hq
  by_cases hie : i = T.e₂
  · subst hie; right; rw [W.e₂] at hi; exact (Option.some.inj hi).symm
  · left; exact List.mem_eraseIdx_iff_getElem?.2 ⟨i, hie, hi⟩

theorem mem_edges_of_mem_glued {p : Nat × Nat} (hp : p ∈ T.glued) :
    p ∈ T.edges ∨ (p.1 ≠ p.2 ∧ EdgesConn T.edges p.1 p.2) := by
  have hu : T.vert₂ T.u₂ = T.u₁ := by simp [vert₂]
  have hv : T.vert₂ T.v₂ = T.v₁ := by simp [vert₂, W.uv₂.symm]
  obtain ⟨q, hq, rfl⟩ := (T.mem_glued).1 hp
  rw [List.mem_append, mem_shiftEdges] at hq
  rcases hq with hq | ⟨r, hr, rfl⟩
  · have hlt := W.emb₁.verts q hq
    rw [T.vmap_lt W hlt.1, T.vmap_lt W hlt.2]
    rcases T.mem_es₁_cases W hq with h | rfl
    · left; unfold edges; exact List.mem_append_left _ h
    · exact Or.inr ⟨W.uv₁, T.edgesConn_edges_uv W⟩
  · simp only [T.vmap_add W]
    rcases T.mem_es₂_cases W hr with h | rfl
    · left; unfold edges; exact List.mem_append_right _ (List.mem_map.2 ⟨r, h, rfl⟩)
    · right
      simp only [hu, hv]
      exact ⟨W.uv₁, T.edgesConn_edges_uv W⟩

theorem mem_glued_of_mem_edges {p : Nat × Nat} (hp : p ∈ T.edges) : p ∈ T.glued := by
  rw [T.mem_glued]
  unfold edges at hp
  rw [List.mem_append, List.mem_map] at hp
  rcases hp with hp | ⟨q, hq, rfl⟩
  · have hlt := W.emb₁.verts p (List.mem_of_mem_eraseIdx hp)
    refine ⟨p, List.mem_append_left _ (List.mem_of_mem_eraseIdx hp), ?_⟩
    exact Prod.ext (T.vmap_lt W hlt.1).symm (T.vmap_lt W hlt.2).symm
  · refine ⟨(T.n₁ + q.1, T.n₁ + q.2),
      List.mem_append_right _ (mem_shiftEdges.2 ⟨q, List.mem_of_mem_eraseIdx hq, rfl⟩), ?_⟩
    exact Prod.ext (T.vmap_add W q.1).symm (T.vmap_add W q.2).symm

/-- Components of the 2-sum, using `conn` (the virtual edge is not a bridge on at least one
side): `C = C₁ + C₂ − 1`. -/
theorem numComponents_edges :
    numComponents T.edges T.nVerts + 1 =
      numComponents T.es₁ T.n₁ + numComponents T.es₂ T.n₂ := by
  have hu₁ := T.u₁_lt' W; have hv₁ := T.v₁_lt' W
  have hu₂ := T.u₂_lt' W; have hv₂ := T.v₂_lt' W
  have huv₂ := W.uv₂
  have h₁ : ∀ p ∈ T.es₁, p.1 < T.n₁ ∧ p.2 < T.n₁ := W.emb₁.verts
  have hmem₁ : (T.u₁, T.v₁) ∈ T.es₁ ++ shiftEdges T.n₁ T.es₂ :=
    List.mem_append_left _ (T.mem_e₁ W)
  have hmem₂ : (T.n₁ + T.u₂, T.n₁ + T.v₂) ∈ T.es₁ ++ shiftEdges T.n₁ T.es₂ :=
    List.mem_append_right _ (mem_shiftEdges.2 ⟨_, T.mem_e₂ W, rfl⟩)
  have hmem₁' : (T.u₁, T.v₁) ∈ identEdges T.u₁ (T.n₁ + T.u₂) (T.es₁ ++ shiftEdges T.n₁ T.es₂) := by
    refine mem_identEdges.2 ⟨_, hmem₁, ?_⟩
    simp only
    rw [ident_of_ne (by omega), ident_of_ne (by omega)]
  have hmem₂' : (T.u₁, T.n₁ + T.v₂) ∈
      identEdges T.u₁ (T.n₁ + T.u₂) (T.es₁ ++ shiftEdges T.n₁ T.es₂) := by
    refine mem_identEdges.2 ⟨_, hmem₂, ?_⟩
    simp only
    rw [ident_b, ident_of_ne (by omega)]
  have hea : HasEdge (T.es₁ ++ shiftEdges T.n₁ T.es₂) T.u₁ := ⟨_, hmem₁, Or.inl rfl⟩
  have heb : HasEdge (T.es₁ ++ shiftEdges T.n₁ T.es₂) (T.n₁ + T.u₂) := ⟨_, hmem₂, Or.inl rfl⟩
  have hea' : HasEdge (identEdges T.u₁ (T.n₁ + T.u₂) (T.es₁ ++ shiftEdges T.n₁ T.es₂)) T.v₁ :=
    ⟨_, hmem₁', Or.inr rfl⟩
  have heb' : HasEdge (identEdges T.u₁ (T.n₁ + T.u₂) (T.es₁ ++ shiftEdges T.n₁ T.es₂))
      (T.n₁ + T.v₂) := ⟨_, hmem₂', Or.inr rfl⟩
  have hcU : ccCount (T.es₁ ++ shiftEdges T.n₁ T.es₂) (T.n₁ + T.n₂) =
      numComponents T.es₁ T.n₁ + numComponents T.es₂ T.n₂ := by
    rw [ccCount_union h₁, numComponents_eq_ccCount h₁, numComponents_eq_ccCount W.emb₂.verts]
  have hcU' : ccCount (identEdges T.u₁ (T.n₁ + T.u₂) (T.es₁ ++ shiftEdges T.n₁ T.es₂))
      (T.n₁ + T.n₂) + 1 = ccCount (T.es₁ ++ shiftEdges T.n₁ T.es₂) (T.n₁ + T.n₂) := by
    apply ccCount_ident_of_not_conn (by omega) (by omega) (by omega) hea heb
    intro h
    have := (edgesConn_union_left h₁ hu₁ h).1
    omega
  have hcG : ccCount T.glued (T.n₁ + T.n₂) =
      ccCount (identEdges T.u₁ (T.n₁ + T.u₂) (T.es₁ ++ shiftEdges T.n₁ T.es₂)) (T.n₁ + T.n₂) := by
    apply ccCount_ident_of_conn (by omega) (by omega) (by omega) hea' heb'
    exact Relation.ReflTransGen.tail (Relation.ReflTransGen.single (Or.inr hmem₁'))
      (Or.inl hmem₂')
  have hcE : numComponents T.edges (T.n₁ + T.n₂) = ccCount T.glued (T.n₁ + T.n₂) := by
    rw [numComponents_eq_ccCount (T.edges_verts' W)]
    exact ccCount_congr_of_conn _ (fun p hp => T.mem_glued_of_mem_edges W hp)
      (fun p hp => T.mem_edges_of_mem_glued W hp)
  show numComponents T.edges (T.n₁ + T.n₂) + 1 = _
  omega

end TwoSum

end Spqr
