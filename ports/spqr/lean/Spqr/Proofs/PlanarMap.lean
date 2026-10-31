import Spqr.Proofs.CompGlue

/-!
# Planar embeddings under vertex relabelling

A rotation system only sees quarter-edges, so an `IsPlanarEmbedding` of `es` is one of
`mapEdges f es` for any `f` that is injective on the non-isolated vertices: the non-isolated
vertex count and the component count are transported by `f`.
-/

namespace Spqr

open Classical

def mapEdges (f : Nat → Nat) (es : List (Nat × Nat)) : List (Nat × Nat) :=
  es.map fun p => (f p.1, f p.2)

theorem mem_mapEdges {f : Nat → Nat} {es : List (Nat × Nat)} {p : Nat × Nat} :
    p ∈ mapEdges f es ↔ ∃ q ∈ es, p = (f q.1, f q.2) := by
  simp only [mapEdges, List.mem_map]
  exact ⟨fun ⟨q, hq, h⟩ => ⟨q, hq, h.symm⟩, fun ⟨q, hq, h⟩ => ⟨q, hq, h.symm⟩⟩

theorem length_mapEdges (f : Nat → Nat) (es : List (Nat × Nat)) :
    (mapEdges f es).length = es.length := List.length_map ..

theorem getElem?_mapEdges (f : Nat → Nat) (es : List (Nat × Nat)) (k : Nat) :
    (mapEdges f es)[k]? = es[k]?.map fun p => (f p.1, f p.2) := List.getElem?_map ..

theorem vert_mapEdges (f : Nat → Nat) (es : List (Nat × Nat)) (q : Nat) :
    QE.vert (mapEdges f es) q = (QE.vert es q).map f := by
  unfold QE.vert
  rw [getElem?_mapEdges]
  cases es[QE.edge q]? with
  | none => rfl
  | some p => simp only [Option.map_some]; split <;> rfl

theorem hasEdge_mapEdges {f : Nat → Nat} {es : List (Nat × Nat)} {v : Nat} :
    HasEdge (mapEdges f es) v ↔ ∃ u, HasEdge es u ∧ f u = v := by
  constructor
  · rintro ⟨p, hp, hpv⟩
    obtain ⟨q, hq, rfl⟩ := mem_mapEdges.1 hp
    rcases hpv with h | h
    · exact ⟨q.1, ⟨q, hq, Or.inl rfl⟩, h⟩
    · exact ⟨q.2, ⟨q, hq, Or.inr rfl⟩, h⟩
  · rintro ⟨u, ⟨q, hq, hqu⟩, rfl⟩
    refine ⟨(f q.1, f q.2), mem_mapEdges.2 ⟨q, hq, rfl⟩, ?_⟩
    rcases hqu with rfl | rfl
    · exact Or.inl rfl
    · exact Or.inr rfl

theorem hasEdge_of_vert {es : List (Nat × Nat)} {q v : Nat} (h : QE.vert es q = some v) :
    HasEdge es v := by
  unfold QE.vert at h
  obtain ⟨p, hp, hpv⟩ := Option.map_eq_some_iff.1 h
  refine ⟨p, List.mem_iff_getElem?.2 ⟨_, hp⟩, ?_⟩
  split at hpv
  · exact Or.inl hpv
  · exact Or.inr hpv

theorem sameVertex_mapEdges {rs : RotationSystem} {es : List (Nat × Nat)} (f : Nat → Nat)
    (h : rs.SameVertex es) : rs.SameVertex (mapEdges f es) := by
  intro q hq r hr
  rw [vert_mapEdges, vert_mapEdges, h q hq r hr]

section Inj

variable {f : Nat → Nat} {es : List (Nat × Nat)} {n n' : Nat}
  (hv : ∀ p ∈ es, p.1 < n ∧ p.2 < n) (hf : ∀ p ∈ es, f p.1 < n' ∧ f p.2 < n')
  (hinj : ∀ u v, HasEdge es u → HasEdge es v → f u = f v → u = v)

include hv in
theorem HasEdge.lt_of {u : Nat} (h : HasEdge es u) : u < n := by
  obtain ⟨p, hp, hpu⟩ := h
  rcases hpu with rfl | rfl
  · exact (hv p hp).1
  · exact (hv p hp).2

include hf in
theorem HasEdge.map_lt {u : Nat} (h : HasEdge es u) : f u < n' := by
  obtain ⟨p, hp, hpu⟩ := h
  rcases hpu with rfl | rfl
  · exact (hf p hp).1
  · exact (hf p hp).2

include hv hf in
theorem filter_hasEdge_mapEdges :
    (Finset.range n').filter (HasEdge (mapEdges f es)) =
      ((Finset.range n).filter (HasEdge es)).image f := by
  ext v
  simp only [Finset.mem_filter, Finset.mem_range, Finset.mem_image, hasEdge_mapEdges]
  constructor
  · rintro ⟨-, u, hu, rfl⟩
    exact ⟨u, ⟨hu.lt_of hv, hu⟩, rfl⟩
  · rintro ⟨u, ⟨-, hu⟩, rfl⟩
    exact ⟨hu.map_lt hf, u, hu, rfl⟩

include hv hf hinj in
theorem numNonIsolated_mapEdges :
    numNonIsolated (mapEdges f es) n' = numNonIsolated es n := by
  rw [numNonIsolated_eq_card, numNonIsolated_eq_card, filter_hasEdge_mapEdges hv hf,
    Finset.card_image_of_injOn]
  intro u hu v hv' h
  exact hinj u v (Finset.mem_filter.1 hu).2 (Finset.mem_filter.1 hv').2 h

include hinj in
theorem edgesConn_mapEdges_iff {x y : Nat} (hx : HasEdge es x) :
    EdgesConn (mapEdges f es) (f x) y ↔ ∃ y', EdgesConn es x y' ∧ f y' = y := by
  constructor
  · intro h
    induction h with
    | refl => exact ⟨x, Relation.ReflTransGen.refl, rfl⟩
    | @tail b c _ hbc ih =>
      obtain ⟨b', hb', rfl⟩ := ih
      have hb'e : HasEdge es b' := hasEdge_of_edgesConn hb' hx
      rcases hbc with h | h
      · obtain ⟨q, hq, hqe⟩ := mem_mapEdges.1 h
        simp only [Prod.mk.injEq] at hqe
        obtain ⟨h1, rfl⟩ := hqe
        have := hinj b' q.1 hb'e ⟨q, hq, Or.inl rfl⟩ h1
        subst this
        exact ⟨q.2, Relation.ReflTransGen.tail hb' (Or.inl hq), rfl⟩
      · obtain ⟨q, hq, hqe⟩ := mem_mapEdges.1 h
        simp only [Prod.mk.injEq] at hqe
        obtain ⟨rfl, h2⟩ := hqe
        have := hinj b' q.2 hb'e ⟨q, hq, Or.inr rfl⟩ h2
        subst this
        exact ⟨q.1, Relation.ReflTransGen.tail hb' (Or.inr hq), rfl⟩
  · rintro ⟨y', hy', rfl⟩
    exact edgesConn_lift f (fun p hp => mem_mapEdges.2 ⟨p, hp, rfl⟩) hy'

include hv hf hinj in
theorem ccCount_mapEdges : ccCount (mapEdges f es) n' = ccCount es n := by
  unfold ccCount
  rw [filter_hasEdge_mapEdges hv hf, Finset.image_image]
  have hcls : ∀ u ∈ (Finset.range n).filter (HasEdge es),
      (Finset.range n').filter (EdgesConn (mapEdges f es) (f u)) =
        ((Finset.range n).filter (EdgesConn es u)).image f := by
    intro u hu
    have hue := (Finset.mem_filter.1 hu).2
    ext y
    simp only [Finset.mem_filter, Finset.mem_range, Finset.mem_image,
      edgesConn_mapEdges_iff hinj hue]
    constructor
    · rintro ⟨-, y', hy', rfl⟩
      exact ⟨y', ⟨(hasEdge_of_edgesConn hy' hue).lt_of hv, hy'⟩, rfl⟩
    · rintro ⟨y', ⟨-, hy'⟩, rfl⟩
      exact ⟨(hasEdge_of_edgesConn hy' hue).map_lt hf, y', hy', rfl⟩
  have h1 : ((Finset.range n).filter (HasEdge es)).image
        ((fun v => (Finset.range n').filter (EdgesConn (mapEdges f es) v)) ∘ f) =
      ((Finset.range n).filter (HasEdge es)).image
        (fun u => ((Finset.range n).filter (EdgesConn es u)).image f) :=
    Finset.image_congr fun u hu => hcls u hu
  have h2 : ((Finset.range n).filter (HasEdge es)).image
        (fun u => ((Finset.range n).filter (EdgesConn es u)).image f) =
      (((Finset.range n).filter (HasEdge es)).image
        (fun u => (Finset.range n).filter (EdgesConn es u))).image (Finset.image f) := by
    rw [Finset.image_image]; rfl
  rw [h1, h2, Finset.card_image_of_injOn]
  intro C hC D hD hCD
  obtain ⟨u, hu, rfl⟩ := Finset.mem_image.1 hC
  obtain ⟨v, hv', rfl⟩ := Finset.mem_image.1 hD
  have hue := (Finset.mem_filter.1 hu).2
  have hve := (Finset.mem_filter.1 hv').2
  have key : ∀ (a b : Nat), HasEdge es a → HasEdge es b →
      ∀ w ∈ (Finset.range n).filter (EdgesConn es a),
      ((Finset.range n).filter (EdgesConn es a)).image f =
        ((Finset.range n).filter (EdgesConn es b)).image f →
      w ∈ (Finset.range n).filter (EdgesConn es b) := by
    intro a b hae _ w hw heq
    have : f w ∈ ((Finset.range n).filter (EdgesConn es b)).image f := by
      rw [← heq]; exact Finset.mem_image_of_mem f hw
    obtain ⟨w', hw', hww'⟩ := Finset.mem_image.1 this
    have hwe : HasEdge es w := hasEdge_of_edgesConn (Finset.mem_filter.1 hw).2 hae
    have hw'e : HasEdge es w' := hasEdge_of_edgesConn (Finset.mem_filter.1 hw').2 ‹HasEdge es b›
    rwa [hinj w' w hw'e hwe hww'] at hw'
  ext w
  exact ⟨fun hw => key u v hue hve w hw hCD, fun hw => key v u hve hue w hw hCD.symm⟩

include hv hf hinj in
theorem IsPlanarEmbedding.map {rs : RotationSystem} (h : IsPlanarEmbedding es n rs) :
    IsPlanarEmbedding (mapEdges f es) n' rs where
  size := by rw [length_mapEdges]; exact h.size
  verts := by
    intro p hp
    obtain ⟨q, hq, rfl⟩ := mem_mapEdges.1 hp
    exact hf q hq
  total := h.total
  involution := h.involution
  opposite_dir := h.opposite_dir
  same_vertex := sameVertex_mapEdges f h.same_vertex
  vertex_orbits := by rw [numNonIsolated_mapEdges hv hf hinj]; exact h.vertex_orbits
  euler := by
    have hes' : ∀ p ∈ mapEdges f es, p.1 < n' ∧ p.2 < n' := by
      intro p hp
      obtain ⟨q, hq, rfl⟩ := mem_mapEdges.1 hp
      exact hf q hq
    unfold EulerFormula
    rw [numNonIsolated_mapEdges hv hf hinj, numComponents_eq_ccCount hes',
      ccCount_mapEdges hv hf hinj, ← numComponents_eq_ccCount hv, length_mapEdges]
    exact h.euler

include hv hf hinj in
theorem Planar.map (h : Planar es n) : Planar (mapEdges f es) n' :=
  let ⟨rs, hrs⟩ := h
  ⟨rs, hrs.map hv hf hinj⟩

end Inj

end Spqr
