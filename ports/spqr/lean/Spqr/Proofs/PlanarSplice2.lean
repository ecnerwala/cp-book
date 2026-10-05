import Spqr.Proofs.PlanarIdent

/-!
# Concrete 2-sum of planar embeddings sharing two vertices

Two planar embeddings on the same vertex set, whose edge sets meet exactly in the vertices `u`
and `v` (each connected to the other inside both), glued by a conjugation at `u` (as in
`IsPlanarEmbedding.splice`) and then one at `v` whose two face pairs are cofacial after the first
gluing.  This is `oneSum_conj` followed by `ident_conj`, relabelled back to `[0, n)`.
-/

namespace Spqr

theorem IsPlanarEmbedding.splice2 {es₁ es₂ : List (Nat × Nat)} {n : Nat}
    {rs₁ rs₂ : RotationSystem} (h₁ : IsPlanarEmbedding es₁ n rs₁)
    (h₂ : IsPlanarEmbedding es₂ n rs₂) {u v a b c d : Nat}
    (ha : a < 4 * es₁.length) (hb : b < 4 * es₂.length) (hab : a % 2 = b % 2)
    (hc : c < 4 * es₁.length) (hd : d < 4 * es₂.length) (hcd : c % 2 = d % 2)
    (hua : QE.vert es₁ a = some u) (hub : QE.vert es₂ b = some u)
    (hvc : QE.vert es₁ c = some v) (hvd : QE.vert es₂ d = some v) (huv : u ≠ v)
    (hconn₁ : EdgesConn es₁ u v) (hconn₂ : EdgesConn es₂ u v)
    (hsep : ∀ w, HasEdge es₁ w → HasEdge es₂ w → w = u ∨ w = v)
    (h1 : SameOrbit (((rs₁.union rs₂).conj a (rs₁.size + b)).stepC 3)
      (c ^^^ 3) ((rs₁.size + d) ^^^ 3))
    (h2 : SameOrbit (((rs₁.union rs₂).conj a (rs₁.size + b)).stepC 3)
      (((rs₁.union rs₂).conj a (rs₁.size + b)).rot c ^^^ 3)
      (((rs₁.union rs₂).conj a (rs₁.size + b)).rot (rs₁.size + d) ^^^ 3)) :
    IsPlanarEmbedding (es₁ ++ es₂) n
      (((rs₁.union rs₂).conj a (rs₁.size + b)).conj c (rs₁.size + d)) := by
  have hu : u < n := (hasEdge_of_vert hua).lt_of h₁.verts
  have hv : v < n := (hasEdge_of_vert hvc).lt_of h₁.verts
  have hs₁ := h₁.size
  have hs₂ := h₂.size
  set R := (rs₁.union rs₂).conj a (rs₁.size + b) with hR
  set E := identEdges u (n + u) (es₁ ++ shiftEdges n es₂) with hE
  have hC : IsPlanarEmbedding E (n + n) R := h₁.oneSum_conj h₂ ha hb hab hua hub
  have hRs : R.size = rs₁.size + rs₂.size := by
    rw [hR, RotationSystem.conj_size, RotationSystem.union_size]
  have hcR : c < R.size := by omega
  have hdR : rs₁.size + d < R.size := by omega
  have hvcE : QE.vert E c = some v := by
    rw [hE, RotationSystem.vert_identEdges, vert_append_left (by omega), hvc]
    simp [ident, show v ≠ n + u by omega]
  have hvdE : QE.vert E (rs₁.size + d) = some (n + v) := by
    rw [hE, RotationSystem.vert_identEdges, hs₁, vert_shiftEdges_add, hvd]
    simp [ident, huv.symm]
  have hconnE : EdgesConn E v (n + v) := by
    have e1 : EdgesConn (es₁ ++ shiftEdges n es₂) v u :=
      edgesConn_union_of_left (edgesConn_symm hconn₁)
    have e2 : EdgesConn (es₁ ++ shiftEdges n es₂) (n + u) (n + v) :=
      edgesConn_union_of_right hconn₂
    have f1 := edgesConn_ident_of (a := u) (b := n + u) e1
    have f2 := edgesConn_ident_of (a := u) (b := n + u) e2
    rw [ident_of_ne (by omega), ident_of_ne (by omega)] at f1
    rw [ident_b, ident_of_ne (by omega)] at f2
    exact f1.trans f2
  have hD : IsPlanarEmbedding (identEdges v (n + v) E) (n + n) (R.conj c (rs₁.size + d)) :=
    hC.ident_conj hcR hdR (by omega) hvcE hvdE (by omega) hconnE h1 h2
  let f := fun x => if x < n then x else x - n
  have hf₁ : ∀ x, HasEdge es₁ x → f x = x := by
    intro x hx
    exact ite_eq_left (hx.lt_of h₁.verts)
  have hf₂ : ∀ x, f (n + x) = x := by
    intro x
    simp [f]
  have hf_ident : ∀ x, f (ident v (n + v) (ident u (n + u) x)) = f x := by
    intro x
    by_cases hx : x = n + u
    · subst x
      simp [ident, f, hu, show u ≠ n + v by omega]
    · by_cases hx' : x = n + v
      · subst x
        simp [ident, f, hv, huv.symm]
      · simp [ident, hx, hx']
  have hmap : mapEdges f (identEdges v (n + v) E) = es₁ ++ es₂ := by
    simp only [hE, mapEdges, identEdges, List.map_map, Function.comp_def, hf_ident,
      List.map_append]
    congr 1
    · trans es₁.map id
      · apply List.map_congr_left
        intro p hp
        change (f p.1, f p.2) = p
        rw [hf₁ p.1 ⟨p, hp, Or.inl rfl⟩, hf₁ p.2 ⟨p, hp, Or.inr rfl⟩]
      · exact List.map_id es₁
    · simp only [shiftEdges, List.map_map, Function.comp_def, hf₂, Prod.mk.eta,
        List.map_id']
  have hinj : ∀ x y, HasEdge (identEdges v (n + v) E) x →
      HasEdge (identEdges v (n + v) E) y → f x = f y → x = y := by
    intro x y hx hy hxy
    obtain ⟨x', rfl, hx'⟩ := (hasEdge_ident_iff x).1 hx
    obtain ⟨y', rfl, hy'⟩ := (hasEdge_ident_iff y).1 hy
    obtain ⟨x'', rfl, hx''⟩ := (hasEdge_ident_iff x').1 hx'
    obtain ⟨y'', rfl, hy''⟩ := (hasEdge_ident_iff y').1 hy'
    rw [hf_ident, hf_ident] at hxy
    rw [hasEdge_union_iff] at hx'' hy''
    rcases hx'' with hx'' | ⟨x'', rfl, hx''⟩ <;> rcases hy'' with hy'' | ⟨y'', rfl, hy''⟩
    · rw [hf₁ x'' hx'', hf₁ y'' hy''] at hxy
      rw [hxy]
    · rw [hf₁ x'' hx'', hf₂] at hxy
      subst y''
      rcases hsep x'' hx'' hy'' with h | h <;> subst x''
      · simp [ident]
      · simp [ident, huv.symm, show v ≠ n + u by omega]
    · rw [hf₂, hf₁ y'' hy''] at hxy
      subst x''
      rcases hsep y'' hy'' hx'' with h | h <;> subst y''
      · simp [ident]
      · simp [ident, huv.symm, show v ≠ n + u by omega]
    · rw [hf₂, hf₂] at hxy
      rw [hxy]
  have hf : ∀ p ∈ identEdges v (n + v) E, f p.1 < n ∧ f p.2 < n := by
    intro p hp
    have hm : (f p.1, f p.2) ∈ es₁ ++ es₂ := by
      rw [← hmap]
      exact mem_mapEdges.2 ⟨p, hp, rfl⟩
    rcases List.mem_append.1 hm with hm | hm
    · exact h₁.verts _ hm
    · exact h₂.verts _ hm
  rw [← hmap]
  exact IsPlanarEmbedding.map hD.verts hf hinj hD

end Spqr
