import Spqr.Proofs.PlanarUnion

namespace Spqr

/-- Union of vertex-disjoint edge lists in the same vertex numbering. -/
theorem IsPlanarEmbedding.append {es₁ es₂ : List (Nat × Nat)} {n : Nat}
    {rs₁ rs₂ : RotationSystem} (h₁ : IsPlanarEmbedding es₁ n rs₁)
    (h₂ : IsPlanarEmbedding es₂ n rs₂)
    (hdis : ∀ v, HasEdge es₁ v → HasEdge es₂ v → False) :
    IsPlanarEmbedding (es₁ ++ es₂) n (rs₁.union rs₂) := by
  let f := fun x => if x < n then x else x - n
  have hf₁ : ∀ x, HasEdge es₁ x → f x = x := by
    intro x hx
    exact ite_eq_left (hx.lt_of h₁.verts)
  have hf₂ : ∀ x, f (n + x) = x := by
    intro x
    simp [f]
  have hmap : mapEdges f (es₁ ++ shiftEdges n es₂) = es₁ ++ es₂ := by
    simp only [mapEdges, List.map_append]
    congr 1
    · trans es₁.map id
      · apply List.map_congr_left
        intro p hp
        change (f p.1, f p.2) = p
        rw [hf₁ p.1 ⟨p, hp, Or.inl rfl⟩, hf₁ p.2 ⟨p, hp, Or.inr rfl⟩]
      · exact List.map_id es₁
    · simp only [shiftEdges, List.map_map, Function.comp_def, hf₂, Prod.mk.eta,
        List.map_id']
  have hU := h₁.union h₂
  have hinj : ∀ x y, HasEdge (es₁ ++ shiftEdges n es₂) x →
      HasEdge (es₁ ++ shiftEdges n es₂) y → f x = f y → x = y := by
    intro x y hx hy hxy
    rw [hasEdge_union_iff] at hx hy
    rcases hx with hx | ⟨x', rfl, hx⟩ <;> rcases hy with hy | ⟨y', rfl, hy⟩
    · rwa [hf₁ x hx, hf₁ y hy] at hxy
    · rw [hf₁ x hx, hf₂] at hxy
      exact False.elim (hdis x hx (hxy ▸ hy))
    · rw [hf₂, hf₁ y hy] at hxy
      exact False.elim (hdis y hy (hxy ▸ hx))
    · rw [hf₂, hf₂] at hxy
      omega
  have hf : ∀ p ∈ es₁ ++ shiftEdges n es₂, f p.1 < n ∧ f p.2 < n := by
    intro p hp
    have hm : (f p.1, f p.2) ∈ es₁ ++ es₂ := by
      rw [← hmap]
      exact mem_mapEdges.2 ⟨p, hp, rfl⟩
    rcases List.mem_append.1 hm with hm | hm
    · exact h₁.verts _ hm
    · exact h₂.verts _ hm
  rw [← hmap]
  exact IsPlanarEmbedding.map hU.verts hf hinj hU

end Spqr
