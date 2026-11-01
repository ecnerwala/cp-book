import Spqr.Proofs.PlanarGlue
import Spqr.Proofs.PlanarMap

/-!
# Disjoint union of planar embeddings

`rs₁.union rs₂` places `rs₂`'s quarter-edges after `rs₁`'s; it is a planar embedding of
`es₁ ++ shiftEdges n₁ es₂` on `n₁ + n₂` vertices (faces, vertices and components add).
-/

namespace Spqr

open Classical

theorem vert_append_left {es es' : List (Nat × Nat)} {q : Nat} (h : q / 4 < es.length) :
    QE.vert (es ++ es') q = QE.vert es q := by
  unfold QE.vert QE.edge
  rw [List.getElem?_append_left h]

theorem vert_shiftEdges_add {es₁ es₂ : List (Nat × Nat)} {n₁ : Nat} (r : Nat) :
    QE.vert (es₁ ++ shiftEdges n₁ es₂) (4 * es₁.length + r) = (QE.vert es₂ r).map (n₁ + ·) := by
  unfold QE.vert QE.edge QE.side
  have h1 : (4 * es₁.length + r) / 4 = es₁.length + r / 4 := by omega
  have h2 : (4 * es₁.length + r) / 2 % 2 = r / 2 % 2 := by omega
  rw [h1, h2, List.getElem?_append_right (by omega), Nat.add_sub_cancel_left]
  unfold shiftEdges
  rw [List.getElem?_map]
  cases es₂[r / 4]? with
  | none => rfl
  | some p => simp only [Option.map_some]; split <;> rfl

theorem xor_sub_mul4 {q m c : Nat} (hq : 4 * m ≤ q) (hc : c < 4) :
    (q ^^^ c) - 4 * m = (q - 4 * m) ^^^ c := by
  obtain ⟨r, rfl⟩ : ∃ r, q = 4 * m + r := ⟨q - 4 * m, by omega⟩
  rw [Nat.add_sub_cancel_left, xor_eq4 (4 * m + r) hc, xor_eq4 r hc]
  have h1 : (4 * m + r) / 4 = m + r / 4 := by omega
  have h2 : (4 * m + r) % 4 = r % 4 := by omega
  rw [h1, h2]
  omega

namespace RotationSystem

/-- Disjoint union: `rs₂`'s quarter-edges are shifted by `rs₁.size`. -/
def union (rs₁ rs₂ : RotationSystem) : RotationSystem :=
  ⟨rs₁.rotAdj ++ rs₂.rotAdj.map (Option.map (rs₁.size + ·))⟩

variable (rs₁ rs₂ : RotationSystem)

theorem union_size : (rs₁.union rs₂).size = rs₁.size + rs₂.size := by
  simp [union, size]

theorem union_get_lt {q : Nat} (hq : q < rs₁.size) : (rs₁.union rs₂).get q = rs₁.get q := by
  unfold union get size at *
  rw [Array.getElem?_append_left hq]

theorem union_get_ge {q : Nat} (hq : rs₁.size ≤ q) :
    (rs₁.union rs₂).get q = (rs₂.get (q - rs₁.size)).map (rs₁.size + ·) := by
  unfold union get size at *
  rw [Array.getElem?_append_right hq, Array.getElem?_map]
  cases rs₂.rotAdj[q - rs₁.rotAdj.size]? with
  | none => rfl
  | some o => cases o <;> rfl

theorem union_stepC {m₁ : Nat} (hs : rs₁.size = 4 * m₁) {c : Nat} (hc : c < 4) :
    (rs₁.union rs₂).stepC c =
      unionStep (rs₁.stepC c) (shift rs₁.size (rs₂.stepC c)) (Finset.range rs₁.size) := by
  funext q
  simp only [stepC_apply, unionStep, Finset.mem_range]
  split_ifs with hq
  · rw [union_get_lt _ _ (hs ▸ xor_lt_mul4 (hs ▸ hq) hc)]
  · have hq' := not_lt.1 hq
    rw [union_get_ge _ _ (hs ▸ ge_of_xor_ge (hs ▸ hq') hc)]
    have hx : (q ^^^ c) - rs₁.size = (q - rs₁.size) ^^^ c := by
      rw [hs] at hq' ⊢; exact xor_sub_mul4 hq' hc
    rw [hx]
    unfold shift
    simp only [hq, ↓reduceIte, stepC_apply]
    cases rs₂.get ((q - rs₁.size) ^^^ c) with
    | none => simp only [Option.map_none, Option.getD_none]; omega
    | some s => rfl

theorem union_total (h₁ : rs₁.Total) (h₂ : rs₂.Total) : (rs₁.union rs₂).Total := by
  intro q hq
  rw [union_size] at hq
  by_cases h : q < rs₁.size
  · rw [union_get_lt _ _ h]; exact h₁ q h
  · rw [union_get_ge _ _ (not_lt.1 h), Option.isSome_map]
    exact h₂ (q - rs₁.size) (by omega)

theorem union_involution (h₁ : rs₁.Involution) (h₂ : rs₂.Involution) :
    (rs₁.union rs₂).Involution := by
  intro q hq r hr
  rw [union_size] at hq
  by_cases h : q < rs₁.size
  · rw [union_get_lt _ _ h] at hr
    obtain ⟨hr1, hr2⟩ := h₁ q h r hr
    rw [union_size, union_get_lt _ _ hr1]
    exact ⟨by omega, hr2⟩
  · rw [union_get_ge _ _ (not_lt.1 h)] at hr
    obtain ⟨r', hr', rfl⟩ := Option.mem_map.1 hr
    obtain ⟨hr1, hr2⟩ := h₂ (q - rs₁.size) (by omega) r' hr'
    rw [union_size, union_get_ge _ _ (by omega), Nat.add_sub_cancel_left, hr2]
    refine ⟨by omega, ?_⟩
    simp only [Option.map_some, Option.some.injEq]
    omega

theorem union_oppositeDir {m₁ : Nat} (hs : rs₁.size = 4 * m₁) (h₁ : rs₁.OppositeDir)
    (h₂ : rs₂.OppositeDir) : (rs₁.union rs₂).OppositeDir := by
  intro q hq r hr
  rw [union_size] at hq
  by_cases h : q < rs₁.size
  · rw [union_get_lt _ _ h] at hr
    exact h₁ q h r hr
  · rw [union_get_ge _ _ (not_lt.1 h)] at hr
    obtain ⟨r', hr', rfl⟩ := Option.mem_map.1 hr
    have := h₂ (q - rs₁.size) (by omega) r' hr'
    unfold QE.dir at *
    omega

theorem union_sameVertex {es₁ es₂ : List (Nat × Nat)} {n₁ : Nat} (hs : rs₁.size = 4 * es₁.length)
    (hi₁ : rs₁.Involution) (h₁ : rs₁.SameVertex es₁) (h₂ : rs₂.SameVertex es₂) :
    (rs₁.union rs₂).SameVertex (es₁ ++ shiftEdges n₁ es₂) := by
  intro q hq r hr
  rw [union_size] at hq
  by_cases h : q < rs₁.size
  · rw [union_get_lt _ _ h] at hr
    have hr1 := (h₁ q h r hr)
    have hrlt : r < rs₁.size := (hi₁ q h r hr).1
    rw [vert_append_left (by omega), vert_append_left (by omega)]
    exact hr1
  · rw [union_get_ge _ _ (not_lt.1 h)] at hr
    obtain ⟨r', hr', rfl⟩ := Option.mem_map.1 hr
    obtain ⟨q', rfl⟩ : ∃ q', q = rs₁.size + q' := ⟨q - rs₁.size, by omega⟩
    rw [Nat.add_sub_cancel_left] at hr'
    have := h₂ q' (by omega) r' hr'
    rw [hs, vert_shiftEdges_add, vert_shiftEdges_add, this]

theorem union_orbitCount (ht₁ : rs₁.Total) (hi₁ : rs₁.Involution) (ht₂ : rs₂.Total)
    (hi₂ : rs₂.Involution) {m₁ m₂ : Nat} (hs₁ : rs₁.size = 4 * m₁) (hs₂ : rs₂.size = 4 * m₂)
    {c : Nat} (hc : c < 4) :
    orbitCount ((rs₁.union rs₂).stepC c) (Finset.range (4 * (m₁ + m₂))) =
      orbitCount (rs₁.stepC c) (Finset.range (4 * m₁)) +
        orbitCount (rs₂.stepC c) (Finset.range (4 * m₂)) := by
  rw [union_stepC _ _ hs₁ hc, show 4 * (m₁ + m₂) = rs₁.size + rs₂.size by omega, ← hs₁, ← hs₂]
  have hT : Finset.range (rs₁.size + rs₂.size) =
      Finset.range rs₁.size ∪ (Finset.range rs₂.size).image (rs₁.size + ·) := by
    ext q
    simp only [Finset.mem_range, Finset.mem_union, Finset.mem_image]
    constructor
    · intro h
      by_cases hq : q < rs₁.size
      · exact Or.inl hq
      · exact Or.inr ⟨q - rs₁.size, by omega, by omega⟩
    · rintro (h | ⟨a, ha, rfl⟩) <;> omega
  have hdisj : Disjoint (Finset.range rs₁.size)
      ((Finset.range rs₂.size).image (rs₁.size + ·)) := by
    rw [Finset.disjoint_left]
    intro q hq hq'
    rw [Finset.mem_range] at hq
    obtain ⟨a, -, rfl⟩ := Finset.mem_image.1 hq'
    omega
  rw [hT, orbitCount_unionStep (isPermOn_stepC ht₁ hi₁ hs₁ hc)
    (isPermOn_shift (isPermOn_stepC ht₂ hi₂ hs₂ hc) _) hdisj,
    orbitCount_shift (isPermOn_stepC ht₂ hi₂ hs₂ hc)]

theorem union_numFaceOrbits (ht₁ : rs₁.Total) (hi₁ : rs₁.Involution) (ht₂ : rs₂.Total)
    (hi₂ : rs₂.Involution) {m₁ m₂ : Nat} (hs₁ : rs₁.size = 4 * m₁) (hs₂ : rs₂.size = 4 * m₂) :
    (rs₁.union rs₂).numFaceOrbits = rs₁.numFaceOrbits + rs₂.numFaceOrbits := by
  rw [numFaceOrbits_eq _ (union_total _ _ ht₁ ht₂) (union_involution _ _ hi₁ hi₂)
      (m := m₁ + m₂) (by rw [union_size, hs₁, hs₂]; omega),
    numFaceOrbits_eq _ ht₁ hi₁ hs₁, numFaceOrbits_eq _ ht₂ hi₂ hs₂]
  exact union_orbitCount _ _ ht₁ hi₁ ht₂ hi₂ hs₁ hs₂ (by decide)

theorem union_numVertexOrbits (ht₁ : rs₁.Total) (hi₁ : rs₁.Involution) (ht₂ : rs₂.Total)
    (hi₂ : rs₂.Involution) {m₁ m₂ : Nat} (hs₁ : rs₁.size = 4 * m₁) (hs₂ : rs₂.size = 4 * m₂) :
    (rs₁.union rs₂).numVertexOrbits = rs₁.numVertexOrbits + rs₂.numVertexOrbits := by
  rw [numVertexOrbits_eq _ (union_total _ _ ht₁ ht₂) (union_involution _ _ hi₁ hi₂)
      (m := m₁ + m₂) (by rw [union_size, hs₁, hs₂]; omega),
    numVertexOrbits_eq _ ht₁ hi₁ hs₁, numVertexOrbits_eq _ ht₂ hi₂ hs₂]
  exact union_orbitCount _ _ ht₁ hi₁ ht₂ hi₂ hs₁ hs₂ (by decide)

end RotationSystem

theorem union_verts {es₁ es₂ : List (Nat × Nat)} {n₁ n₂ : Nat}
    (h₁ : ∀ p ∈ es₁, p.1 < n₁ ∧ p.2 < n₁) (h₂ : ∀ p ∈ es₂, p.1 < n₂ ∧ p.2 < n₂) :
    ∀ p ∈ es₁ ++ shiftEdges n₁ es₂, p.1 < n₁ + n₂ ∧ p.2 < n₁ + n₂ := by
  intro p hp
  rw [List.mem_append, mem_shiftEdges] at hp
  rcases hp with hp | ⟨q, hq, rfl⟩
  · have := h₁ p hp; omega
  · have := h₂ q hq; simp only; omega

theorem numNonIsolated_union {es₁ es₂ : List (Nat × Nat)} {n₁ n₂ : Nat}
    (h₁ : ∀ p ∈ es₁, p.1 < n₁ ∧ p.2 < n₁) (h₂ : ∀ p ∈ es₂, p.1 < n₂ ∧ p.2 < n₂) :
    numNonIsolated (es₁ ++ shiftEdges n₁ es₂) (n₁ + n₂) =
      numNonIsolated es₁ n₁ + numNonIsolated es₂ n₂ := by
  rw [numNonIsolated_eq_card, numNonIsolated_eq_card, numNonIsolated_eq_card]
  have hset : (Finset.range (n₁ + n₂)).filter (HasEdge (es₁ ++ shiftEdges n₁ es₂)) =
      (Finset.range n₁).filter (HasEdge es₁) ∪
        ((Finset.range n₂).filter (HasEdge es₂)).image (n₁ + ·) := by
    ext v
    simp only [Finset.mem_filter, Finset.mem_range, Finset.mem_union, Finset.mem_image,
      hasEdge_union_iff]
    constructor
    · rintro ⟨hv, h | ⟨v', rfl, h⟩⟩
      · exact Or.inl ⟨h.lt_of h₁, h⟩
      · exact Or.inr ⟨v', ⟨h.lt_of h₂, h⟩, rfl⟩
    · rintro (⟨hv, h⟩ | ⟨v', ⟨hv', h⟩, rfl⟩)
      · exact ⟨by omega, Or.inl h⟩
      · exact ⟨by omega, Or.inr ⟨v', rfl, h⟩⟩
  rw [hset, Finset.card_union_of_disjoint, Finset.card_image_of_injective _ (add_right_injective n₁)]
  rw [Finset.disjoint_left]
  intro v hv hv'
  obtain ⟨v', -, rfl⟩ := Finset.mem_image.1 hv'
  have := (Finset.mem_filter.1 hv).1
  rw [Finset.mem_range] at this
  omega

theorem IsPlanarEmbedding.union {es₁ es₂ : List (Nat × Nat)} {n₁ n₂ : Nat}
    {rs₁ rs₂ : RotationSystem} (h₁ : IsPlanarEmbedding es₁ n₁ rs₁)
    (h₂ : IsPlanarEmbedding es₂ n₂ rs₂) :
    IsPlanarEmbedding (es₁ ++ shiftEdges n₁ es₂) (n₁ + n₂) (rs₁.union rs₂) where
  size := by
    rw [RotationSystem.union_size, h₁.size, h₂.size]
    simp only [List.length_append, shiftEdges, List.length_map]
    omega
  verts := union_verts h₁.verts h₂.verts
  total := RotationSystem.union_total _ _ h₁.total h₂.total
  involution := RotationSystem.union_involution _ _ h₁.involution h₂.involution
  opposite_dir := RotationSystem.union_oppositeDir _ _ h₁.size h₁.opposite_dir h₂.opposite_dir
  same_vertex :=
    RotationSystem.union_sameVertex _ _ h₁.size h₁.involution h₁.same_vertex h₂.same_vertex
  vertex_orbits := by
    rw [RotationSystem.union_numVertexOrbits _ _ h₁.total h₁.involution h₂.total h₂.involution
      h₁.size h₂.size, numNonIsolated_union h₁.verts h₂.verts, h₁.vertex_orbits, h₂.vertex_orbits]
    omega
  euler := by
    unfold EulerFormula
    rw [RotationSystem.union_numFaceOrbits _ _ h₁.total h₁.involution h₂.total h₂.involution
      h₁.size h₂.size, numNonIsolated_union h₁.verts h₂.verts,
      numComponents_eq_ccCount (union_verts h₁.verts h₂.verts), ccCount_union h₁.verts n₂,
      ← numComponents_eq_ccCount h₁.verts, ← numComponents_eq_ccCount h₂.verts]
    simp only [List.length_append, shiftEdges, List.length_map]
    have e₁ := h₁.euler; have e₂ := h₂.euler
    unfold EulerFormula at e₁ e₂
    omega

theorem Planar.union {es₁ es₂ : List (Nat × Nat)} {n₁ n₂ : Nat} (h₁ : Planar es₁ n₁)
    (h₂ : Planar es₂ n₂) : Planar (es₁ ++ shiftEdges n₁ es₂) (n₁ + n₂) :=
  let ⟨_, h₁⟩ := h₁
  let ⟨_, h₂⟩ := h₂
  ⟨_, h₁.union h₂⟩

end Spqr
