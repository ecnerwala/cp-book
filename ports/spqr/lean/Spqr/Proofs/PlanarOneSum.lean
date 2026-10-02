import Spqr.Proofs.PlanarUnion

/-!
# 1-sum of planar embeddings

`rs.conj a b` conjugates the rotation by the transposition `(a b)` of two quarter-edges of the
same direction: the two rotations at the vertices of `a` and `b` are spliced into one.  On the
disjoint union `rs₁.union rs₂` with `a` at `v₁` and `b` at `n₁ + v₂` this embeds
`identEdges v₁ (n₁ + v₂) (es₁ ++ shiftEdges n₁ es₂)`: two face orbits and two vertex orbits
merge, one component disappears.
-/

namespace Spqr

open Classical

theorem xor_right_inj {p x c : Nat} : p ^^^ c = x ^^^ c ↔ p = x :=
  ⟨fun h => by rw [← xor_xor_self p c, h, xor_xor_self], fun h => by rw [h]⟩

namespace RotationSystem

/-- Conjugate the rotation by the transposition `(a b)`. -/
def conj (rs : RotationSystem) (a b : Nat) : RotationSystem :=
  ⟨(Array.range rs.size).map fun q => (rs.get (Equiv.swap a b q)).map (Equiv.swap a b)⟩

variable (rs : RotationSystem) (a b : Nat)

theorem conj_size : (rs.conj a b).size = rs.size := by
  simp [conj, size]

theorem conj_get_lt {q : Nat} (hq : q < rs.size) :
    (rs.conj a b).get q = (rs.get (Equiv.swap a b q)).map (Equiv.swap a b) := by
  unfold conj get
  rw [Array.getElem?_map, Array.getElem?_range]
  simp [hq]

theorem conj_get_ge {q : Nat} (hq : rs.size ≤ q) : (rs.conj a b).get q = none := by
  unfold conj get
  rw [Array.getElem?_map, Array.getElem?_range]
  simp [not_lt.2 hq]

section

variable (ht : rs.Total) (hi : rs.Involution) (ho : rs.OppositeDir)
  (ha : a < rs.size) (hb : b < rs.size) (hab : a ≠ b) (hab2 : a % 2 = b % 2)
  (hrab : rs.rot a ≠ b)

include ha hb

theorem swap_lt {q : Nat} (hq : q < rs.size) : Equiv.swap a b q < rs.size := by
  rw [Equiv.swap_apply_def]
  split_ifs <;> assumption

theorem swap_lt_iff {q : Nat} : Equiv.swap a b q < rs.size ↔ q < rs.size := by
  constructor
  · intro h
    have := swap_lt rs a b ha hb h
    rwa [Equiv.swap_apply_self] at this
  · exact swap_lt rs a b ha hb

include hab2 in
theorem swap_mod2 (q : Nat) : Equiv.swap a b q % 2 = q % 2 := by
  rw [Equiv.swap_apply_def]
  split_ifs with h1 h2
  · subst h1; omega
  · subst h2; omega
  · rfl

include ht in
theorem conj_total : (rs.conj a b).Total := by
  intro q hq
  rw [conj_size] at hq
  rw [conj_get_lt _ _ _ hq, Option.isSome_map]
  exact ht _ (swap_lt rs a b ha hb hq)

include ht hi in
theorem conj_involution : (rs.conj a b).Involution := by
  intro q hq r hr
  rw [conj_size] at hq ⊢
  rw [conj_get_lt _ _ _ hq, get_eq_rot ht (swap_lt rs a b ha hb hq)] at hr
  simp only [Option.map_some, Option.mem_def, Option.some.injEq] at hr
  subst hr
  have hlt : rs.rot (Equiv.swap a b q) < rs.size := rot_lt ht hi (swap_lt rs a b ha hb hq)
  refine ⟨swap_lt rs a b ha hb hlt, ?_⟩
  rw [conj_get_lt _ _ _ (swap_lt rs a b ha hb hlt), Equiv.swap_apply_self,
    get_eq_rot ht hlt, rot_rot ht hi (swap_lt rs a b ha hb hq)]
  simp

include ht ho hab2 in
theorem conj_oppositeDir : (rs.conj a b).OppositeDir := by
  intro q hq r hr
  rw [conj_size] at hq
  rw [conj_get_lt _ _ _ hq, get_eq_rot ht (swap_lt rs a b ha hb hq)] at hr
  simp only [Option.map_some, Option.mem_def, Option.some.injEq] at hr
  subst hr
  unfold QE.dir
  rw [swap_mod2 rs a b ha hb hab2]
  have := rot_mod2 ht ho (swap_lt rs a b ha hb hq)
  rw [swap_mod2 rs a b ha hb hab2] at this
  exact this

omit ha hb in
theorem vert_identEdges (v w : Nat) (es : List (Nat × Nat)) (q : Nat) :
    QE.vert (identEdges v w es) q = (QE.vert es q).map (ident v w) :=
  vert_mapEdges (ident v w) es q

include ht in
theorem conj_sameVertex {es : List (Nat × Nat)} {v w : Nat}
    (hsv : rs.SameVertex es) (hva : QE.vert es a = some v) (hwb : QE.vert es b = some w) :
    (rs.conj a b).SameVertex (identEdges v w es) := by
  have key : ∀ x, QE.vert (identEdges v w es) (Equiv.swap a b x) = QE.vert (identEdges v w es) x := by
    intro x
    rw [Equiv.swap_apply_def]
    split_ifs with h1 h2
    · subst h1; rw [vert_identEdges, vert_identEdges, hva, hwb]; simp [ident]
    · subst h2; rw [vert_identEdges, vert_identEdges, hva, hwb]; simp [ident]
    · rfl
  intro q hq r hr
  rw [conj_size] at hq
  rw [conj_get_lt _ _ _ hq, get_eq_rot ht (swap_lt rs a b ha hb hq)] at hr
  simp only [Option.map_some, Option.mem_def, Option.some.injEq] at hr
  subst hr
  rw [key, ← key q, vert_identEdges, vert_identEdges,
    hsv _ (swap_lt rs a b ha hb hq) _ (by rw [get_eq_rot ht (swap_lt rs a b ha hb hq)]; rfl)]

include ht hi ho hab hrab in
theorem conj_stepC {m c : Nat} (hs : rs.size = 4 * m) (hc : c < 4) :
    (rs.conj a b).stepC c =
      swapImg (swapImg (rs.stepC c) (a ^^^ c) (b ^^^ c)) (rs.rot a ^^^ c) (rs.rot b ^^^ c) := by
  have hra := rot_lt ht hi ha
  have hrb := rot_lt ht hi hb
  have hraa := rot_ne ht ho ha
  have hrbb := rot_ne ht ho hb
  have hrba : rs.rot b ≠ a := fun h => hrab (by rw [← h, rot_rot ht hi hb])
  have hrr : rs.rot a ≠ rs.rot b := fun h => hab (rot_inj ht hi ha hb h)
  have hσ : ∀ z, z < rs.size → rs.stepC c (z ^^^ c) = rs.rot z := by
    intro z hz
    rw [stepC_eq_rot ht hs hc (hs ▸ xor_lt_mul4 (hs ▸ hz) hc), xor_xor_self]
  funext q
  by_cases hq : q < rs.size
  · obtain ⟨p, rfl⟩ : ∃ p, q = p ^^^ c := ⟨q ^^^ c, (xor_xor_self q c).symm⟩
    have hp : p < rs.size := by
      have := hs ▸ xor_lt_mul4 (hs ▸ hq) hc
      rwa [xor_xor_self] at this
    rw [stepC_apply, xor_xor_self, conj_get_lt _ _ _ hp, get_eq_rot ht (swap_lt rs a b ha hb hp)]
    simp only [Option.map_some, Option.getD_some, swapImg, xor_right_inj, hσ _ ha, hσ _ hb,
      hσ _ hra, hσ _ hrb, hσ _ hp, rot_rot ht hi ha, rot_rot ht hi hb, hrba, hrbb, hraa, hrab,
      ↓reduceIte]
    by_cases h1 : p = a
    · subst h1
      simp [Equiv.swap_apply_left, Equiv.swap_apply_of_ne_of_ne hrba hrbb, hraa.symm, hrba.symm]
    by_cases h2 : p = b
    · subst h2
      simp [Equiv.swap_apply_right, Equiv.swap_apply_of_ne_of_ne hraa hrab, hrab.symm, hrbb.symm,
        hab.symm]
    by_cases h3 : p = rs.rot a
    · subst h3
      simp [Equiv.swap_apply_of_ne_of_ne hraa hrab, rot_rot ht hi ha, Equiv.swap_apply_left]
    by_cases h4 : p = rs.rot b
    · subst h4
      simp [Equiv.swap_apply_of_ne_of_ne hrba hrbb, rot_rot ht hi hb, Equiv.swap_apply_right,
        hrr.symm]
    have h5 : rs.rot p ≠ a := fun h => h3 (by rw [← h, rot_rot ht hi hp])
    have h6 : rs.rot p ≠ b := fun h => h4 (by rw [← h, rot_rot ht hi hp])
    simp [Equiv.swap_apply_of_ne_of_ne h1 h2, Equiv.swap_apply_of_ne_of_ne h5 h6, h1, h2, h3, h4]
  · have hq' := not_lt.1 hq
    have h1 : q ≠ a ^^^ c := fun h => hq (h ▸ hs ▸ xor_lt_mul4 (hs ▸ ha) hc)
    have h2 : q ≠ b ^^^ c := fun h => hq (h ▸ hs ▸ xor_lt_mul4 (hs ▸ hb) hc)
    have h3 : q ≠ rs.rot a ^^^ c := fun h => hq (h ▸ hs ▸ xor_lt_mul4 (hs ▸ hra) hc)
    have h4 : q ≠ rs.rot b ^^^ c := fun h => hq (h ▸ hs ▸ xor_lt_mul4 (hs ▸ hrb) hc)
    rw [stepC_apply, conj_get_ge _ _ _ (hs ▸ ge_of_xor_ge (hs ▸ hq') hc)]
    simp [swapImg, h1, h2, h3, h4, stepC_of_ge hs hc hq']

include ht hi ho hab hrab in
theorem conj_orbitCount {m c : Nat} (hs : rs.size = 4 * m) (hc : c < 4) (hc2 : c % 2 = 1)
    (C : Nat → Prop) (hC : ∀ z, C z ↔ C (rs.stepC c z)) (hCa : C (a ^^^ c)) (hCb : ¬C (b ^^^ c))
    (hCra : C (rs.rot a ^^^ c)) (hCrb : ¬C (rs.rot b ^^^ c)) :
    orbitCount ((rs.conj a b).stepC c) (Finset.range rs.size) + 2 =
      orbitCount (rs.stepC c) (Finset.range rs.size) := by
  rw [conj_stepC rs a b ht hi ho ha hb hab hrab hs hc]
  have hf := isPermOn_stepC ht hi hs hc
  have hinv : ∀ x y, SameOrbit (rs.stepC c) x y → (C x ↔ C y) := fun x y h =>
    Iff.of_eq (sameOrbit_invariant (fun z => C z) (fun z => propext (hC z).symm) h)
  have hpar : ∀ x y, SameOrbit (rs.stepC c) x y → x % 2 = y % 2 := fun x y h =>
    sameOrbit_invariant (· % 2) (stepC_mod2 ht ho hs hc hc2) h
  have hra := rot_lt ht hi ha
  have hrb := rot_lt ht hi hb
  have mem : ∀ z, z < rs.size → z ^^^ c ∈ Finset.range rs.size := fun z hz =>
    Finset.mem_range.2 (hs ▸ xor_lt_mul4 (hs ▸ hz) hc)
  have n1 : ¬SameOrbit (rs.stepC c) (a ^^^ c) (b ^^^ c) := fun h => hCb ((hinv _ _ h).1 hCa)
  have n2 : ¬SameOrbit (rs.stepC c) (rs.rot a ^^^ c) (rs.rot b ^^^ c) :=
    fun h => hCrb ((hinv _ _ h).1 hCra)
  have nxy : ¬SameOrbit (rs.stepC c) (a ^^^ c) (rs.rot b ^^^ c) :=
    fun h => hCrb ((hinv _ _ h).1 hCa)
  have nyx : ¬SameOrbit (rs.stepC c) (b ^^^ c) (rs.rot a ^^^ c) :=
    fun h => hCb ((hinv _ _ h).2 hCra)
  have hpa : ¬(SameOrbit (rs.stepC c) (a ^^^ c) (rs.rot a ^^^ c) ∧
      SameOrbit (rs.stepC c) (b ^^^ c) (rs.rot b ^^^ c)) := by
    rintro ⟨h, -⟩
    have h1 := hpar _ _ h
    have h2 := xor_mod2 a hc2
    have h3 := xor_mod2 (rs.rot a) hc2
    have h4 := rot_mod2 ht ho ha
    omega
  have e1 := orbitCount_swapImg hf (mem a ha) (mem b hb) n1
  have e2 := orbitCount_swapImg (isPermOn_swapImg hf (mem a ha) (mem b hb)) (mem _ hra) (mem _ hrb)
    (not_sameOrbit_swapImg_of hf (mem a ha) (mem b hb) n1 n2 nxy nyx hpa)
  omega

include ht hi ho hab hrab in
theorem conj_numFaceOrbits {m : Nat} (hs : rs.size = 4 * m)
    (C : Nat → Prop) (hC : ∀ z, C z ↔ C (rs.stepC 3 z)) (hCa : C (a ^^^ 3)) (hCb : ¬C (b ^^^ 3))
    (hCra : C (rs.rot a ^^^ 3)) (hCrb : ¬C (rs.rot b ^^^ 3)) :
    (rs.conj a b).numFaceOrbits + 2 = rs.numFaceOrbits := by
  rw [numFaceOrbits_eq _ (conj_total rs a b ht ha hb) (conj_involution rs a b ht hi ha hb)
      (m := m) (by rw [conj_size, hs]), numFaceOrbits_eq _ ht hi hs, ← hs]
  exact conj_orbitCount rs a b ht hi ho ha hb hab hrab hs (by decide) (by decide) C hC hCa hCb
    hCra hCrb

include ht hi ho hab hrab in
theorem conj_numVertexOrbits {m : Nat} (hs : rs.size = 4 * m)
    (C : Nat → Prop) (hC : ∀ z, C z ↔ C (rs.stepC 1 z)) (hCa : C (a ^^^ 1)) (hCb : ¬C (b ^^^ 1))
    (hCra : C (rs.rot a ^^^ 1)) (hCrb : ¬C (rs.rot b ^^^ 1)) :
    (rs.conj a b).numVertexOrbits + 2 = rs.numVertexOrbits := by
  rw [numVertexOrbits_eq _ (conj_total rs a b ht ha hb) (conj_involution rs a b ht hi ha hb)
      (m := m) (by rw [conj_size, hs]), numVertexOrbits_eq _ ht hi hs, ← hs]
  exact conj_orbitCount rs a b ht hi ho ha hb hab hrab hs (by decide) (by decide) C hC hCa hCb
    hCra hCrb

end

end RotationSystem

theorem HasEdge.exists_quarter {es : List (Nat × Nat)} {v : Nat} (h : HasEdge es v) :
    ∃ a, a < 4 * es.length ∧ a % 2 = 0 ∧ QE.vert es a = some v := by
  obtain ⟨p, hp, hpv⟩ := h
  obtain ⟨i, hi, rfl⟩ := List.mem_iff_getElem.1 hp
  rcases hpv with h | h
  · refine ⟨4 * i, by omega, by omega, ?_⟩
    unfold QE.vert
    rw [show QE.edge (4 * i) = i by unfold QE.edge; omega,
      show QE.side (4 * i) = 0 by unfold QE.side; omega, List.getElem?_eq_getElem hi]
    simp [h]
  · refine ⟨4 * i + 2, by omega, by omega, ?_⟩
    unfold QE.vert
    rw [show QE.edge (4 * i + 2) = i by unfold QE.edge; omega,
      show QE.side (4 * i + 2) = 1 by unfold QE.side; omega, List.getElem?_eq_getElem hi]
    simp [h]

theorem numNonIsolated_ident {es : List (Nat × Nat)} {n a b : Nat} (hab : a ≠ b) (ha : a < n)
    (hb : b < n) (hea : HasEdge es a) (heb : HasEdge es b) :
    numNonIsolated (identEdges a b es) n + 1 = numNonIsolated es n := by
  rw [numNonIsolated_eq_card, numNonIsolated_eq_card]
  have : (Finset.range n).filter (HasEdge es) =
      insert b ((Finset.range n).filter (HasEdge (identEdges a b es))) := by
    ext v
    simp only [Finset.mem_insert, Finset.mem_filter, Finset.mem_range]
    by_cases hv : v = b
    · subst hv; simp [hb, heb]
    · simp only [hv, false_or]
      by_cases hva : v = a
      · subst hva; simp [ha, hea, hasEdge_ident_a hab hea]
      · rw [hasEdge_ident_of_ne hva hv]
  rw [this, Finset.card_insert_of_notMem]
  rw [Finset.mem_filter]
  exact fun h => not_hasEdge_ident_b hab h.2

/-- The 1-sum on the disjoint union: identify `v₁` (with an edge in `es₁`) and `n₁ + v₂` (with
an edge in `es₂`) by conjugating `rs₁.union rs₂` with a transposition of two of their
quarter-edges. -/
theorem IsPlanarEmbedding.oneSum {es₁ es₂ : List (Nat × Nat)} {n₁ n₂ : Nat}
    {rs₁ rs₂ : RotationSystem} (h₁ : IsPlanarEmbedding es₁ n₁ rs₁)
    (h₂ : IsPlanarEmbedding es₂ n₂ rs₂) {v₁ v₂ : Nat} (hv₁ : HasEdge es₁ v₁)
    (hv₂ : HasEdge es₂ v₂) :
    Planar (identEdges v₁ (n₁ + v₂) (es₁ ++ shiftEdges n₁ es₂)) (n₁ + n₂) := by
  obtain ⟨a, ha, ha2, hva⟩ := hv₁.exists_quarter
  obtain ⟨b₂, hb₂, hb2, hvb⟩ := hv₂.exists_quarter
  have hU := h₁.union h₂
  have hs₁ := h₁.size
  have hs₂ := h₂.size
  have hsz := hU.size
  have hv₁n : v₁ < n₁ := hv₁.lt_of h₁.verts
  have hv₂n : v₂ < n₂ := hv₂.lt_of h₂.verts
  have hw : v₁ ≠ n₁ + v₂ := by omega
  set rs := rs₁.union rs₂ with hrs
  set es := es₁ ++ shiftEdges n₁ es₂ with hes
  set b := rs₁.size + b₂ with hbdef
  have hsU : rs.size = rs₁.size + rs₂.size := RotationSystem.union_size rs₁ rs₂
  have ha' : a < rs.size := by omega
  have hb : b < rs.size := by omega
  have hab : a ≠ b := by omega
  have hab2 : a % 2 = b % 2 := by omega
  have hrota : rs.rot a < rs₁.size := by
    have := rs.get_eq_rot hU.total ha'
    rw [RotationSystem.union_get_lt _ _ (by omega)] at this
    exact (h₁.involution a (by omega) _ (by rw [this]; rfl)).1
  have hrotb : rs₁.size ≤ rs.rot b := by
    have := rs.get_eq_rot hU.total hb
    rw [RotationSystem.union_get_ge _ _ (by omega), Option.map_eq_some_iff] at this
    obtain ⟨r', -, hr'⟩ := this
    omega
  have hrab : rs.rot a ≠ b := by omega
  have hC : ∀ c, c < 4 → ∀ z, z < rs₁.size ↔ rs.stepC c z < rs₁.size := by
    intro c hc z
    rw [RotationSystem.union_stepC _ _ hs₁ hc]
    unfold unionStep shift
    simp only [Finset.mem_range]
    by_cases hz : z < rs₁.size
    · simp only [hz, ↓reduceIte, true_iff]
      exact RotationSystem.stepC_lt h₁.total h₁.involution hs₁ hc hz
    · simp only [hz, ↓reduceIte, false_iff, not_lt]
      omega
  have hCa : ∀ c, c < 4 → a ^^^ c < rs₁.size := fun c hc => hs₁ ▸ xor_lt_mul4 (hs₁ ▸ ha) hc
  have hCb : ∀ c, c < 4 → ¬b ^^^ c < rs₁.size := fun c hc =>
    not_lt.2 (hs₁ ▸ ge_of_xor_ge (hs₁ ▸ (by omega : rs₁.size ≤ b)) hc)
  have hCra : ∀ c, c < 4 → rs.rot a ^^^ c < rs₁.size := fun c hc =>
    hs₁ ▸ xor_lt_mul4 (hs₁ ▸ hrota) hc
  have hCrb : ∀ c, c < 4 → ¬rs.rot b ^^^ c < rs₁.size := fun c hc =>
    not_lt.2 (hs₁ ▸ ge_of_xor_ge (hs₁ ▸ hrotb) hc)
  have hva' : QE.vert es a = some v₁ := by
    rw [hes, vert_append_left (by omega)]; exact hva
  have hwb : QE.vert es b = some (n₁ + v₂) := by
    rw [hes, hbdef, hs₁, vert_shiftEdges_add, hvb]; rfl
  have hev : ∀ p ∈ identEdges v₁ (n₁ + v₂) es, p.1 < n₁ + n₂ ∧ p.2 < n₁ + n₂ := by
    intro p hp
    obtain ⟨q, hq, rfl⟩ := mem_identEdges.1 hp
    have := hU.verts q hq
    simp only [ident]
    split_ifs <;> omega
  have hea : HasEdge es v₁ := hasEdge_union_iff.2 (Or.inl hv₁)
  have heb : HasEdge es (n₁ + v₂) := hasEdge_union_iff.2 (Or.inr ⟨v₂, rfl, hv₂⟩)
  have hnc : ¬EdgesConn es v₁ (n₁ + v₂) := fun h =>
    by have := (edgesConn_union_left h₁.verts hv₁n h).1; omega
  have hVO := RotationSystem.conj_numVertexOrbits rs a b hU.total hU.involution hU.opposite_dir
    ha' hb hab hrab hsz (fun z => z < rs₁.size) (hC 1 (by decide)) (hCa 1 (by decide))
    (hCb 1 (by decide)) (hCra 1 (by decide)) (hCrb 1 (by decide))
  have hFO := RotationSystem.conj_numFaceOrbits rs a b hU.total hU.involution hU.opposite_dir
    ha' hb hab hrab hsz (fun z => z < rs₁.size) (hC 3 (by decide)) (hCa 3 (by decide))
    (hCb 3 (by decide)) (hCra 3 (by decide)) (hCrb 3 (by decide))
  have hNI := numNonIsolated_ident (n := n₁ + n₂) hw (by omega) (by omega) hea heb
  have hCC := ccCount_ident_of_not_conn (n := n₁ + n₂) hw (by omega) (by omega) hea heb hnc
  refine ⟨rs.conj a b, ?_⟩
  exact {
    size := by
      rw [RotationSystem.conj_size, hsz]
      simp [identEdges]
    verts := hev
    total := RotationSystem.conj_total rs a b hU.total ha' hb
    involution := RotationSystem.conj_involution rs a b hU.total hU.involution ha' hb
    opposite_dir := RotationSystem.conj_oppositeDir rs a b hU.total hU.opposite_dir ha' hb hab2
    same_vertex := RotationSystem.conj_sameVertex rs a b hU.total ha' hb hU.same_vertex hva' hwb
    vertex_orbits := by
      have := hU.vertex_orbits
      omega
    euler := by
      unfold EulerFormula
      have hE := hU.euler
      unfold EulerFormula at hE
      rw [numComponents_eq_ccCount hev, numComponents_eq_ccCount hU.verts] at *
      simp only [identEdges, List.length_map] at *
      omega }

end Spqr
