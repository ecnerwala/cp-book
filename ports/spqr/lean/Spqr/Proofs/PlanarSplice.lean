import Spqr.Proofs.PlanarOneSum

namespace Spqr

theorem IsPlanarEmbedding.splice {es₁ es₂ : List (Nat × Nat)} {n : Nat}
    {rs₁ rs₂ : RotationSystem} (h₁ : IsPlanarEmbedding es₁ n rs₁)
    (h₂ : IsPlanarEmbedding es₂ n rs₂) {v a b : Nat}
    (ha : a < 4 * es₁.length) (hb : b < 4 * es₂.length) (hab : a % 2 = b % 2)
    (hva : QE.vert es₁ a = some v) (hvb : QE.vert es₂ b = some v)
    (hsep : ∀ w, HasEdge es₁ w → HasEdge es₂ w → w = v) :
    IsPlanarEmbedding (es₁ ++ es₂) n ((rs₁.union rs₂).conj a (rs₁.size + b)) := by
  have hv : v < n := (hasEdge_of_vert hva).lt_of h₁.verts
  let f := fun x => if x < n then x else x - n
  have hf₁ : ∀ x, HasEdge es₁ x → f x = x := by
    intro x hx
    exact ite_eq_left (hx.lt_of h₁.verts)
  have hf₂ : ∀ x, f (n + x) = x := by
    intro x
    simp [f]
  have hf_ident : ∀ x, f (ident v (n + v) x) = f x := by
    intro x
    by_cases hx : x = n + v
    · subst x
      simp [ident, f, hv]
    · simp [ident, hx]
  have hmap : mapEdges f (identEdges v (n + v) (es₁ ++ shiftEdges n es₂)) = es₁ ++ es₂ := by
    simp only [mapEdges, identEdges, List.map_map, Function.comp_def, hf_ident, List.map_append]
    congr 1
    · trans es₁.map id
      · apply List.map_congr_left
        intro p hp
        change (f p.1, f p.2) = p
        rw [hf₁ p.1 ⟨p, hp, Or.inl rfl⟩, hf₁ p.2 ⟨p, hp, Or.inr rfl⟩]
      · exact List.map_id es₁
    · simp only [shiftEdges, List.map_map, Function.comp_def, hf₂, Prod.mk.eta,
        List.map_id']
  have hinj : ∀ x y, HasEdge (identEdges v (n + v) (es₁ ++ shiftEdges n es₂)) x →
      HasEdge (identEdges v (n + v) (es₁ ++ shiftEdges n es₂)) y → f x = f y → x = y := by
    intro x y hx hy hxy
    obtain ⟨x', rfl, hx'⟩ := (hasEdge_ident_iff x).1 hx
    obtain ⟨y', rfl, hy'⟩ := (hasEdge_ident_iff y).1 hy
    rw [hf_ident, hf_ident] at hxy
    rw [hasEdge_union_iff] at hx' hy'
    rcases hx' with hx' | ⟨x', rfl, hx'⟩ <;> rcases hy' with hy' | ⟨y', rfl, hy'⟩
    · rw [hf₁ x' hx', hf₁ y' hy'] at hxy
      rw [hxy]
    · rw [hf₁ x' hx', hf₂] at hxy
      subst y'
      have hxv := hsep x' hx' hy'
      subst x'
      simp [ident]
    · rw [hf₂, hf₁ y' hy'] at hxy
      subst x'
      have hyv := hsep y' hy' hx'
      subst y'
      simp [ident]
    · rw [hf₂, hf₂] at hxy
      rw [hxy]
  have hf : ∀ p ∈ identEdges v (n + v) (es₁ ++ shiftEdges n es₂),
      f p.1 < n ∧ f p.2 < n := by
    intro p hp
    have hm : (f p.1, f p.2) ∈ es₁ ++ es₂ := by
      rw [← hmap]
      exact mem_mapEdges.2 ⟨p, hp, rfl⟩
    rcases List.mem_append.1 hm with hm | hm
    · exact h₁.verts _ hm
    · exact h₂.verts _ hm
  have hC := h₁.oneSum_conj h₂ ha hb hab hva hvb
  rw [← hmap]
  exact IsPlanarEmbedding.map hC.verts hf hinj hC

end Spqr
