import Spqr.Proofs.PlanarInsert
import Spqr.Proofs.OrbitSplit
import Spqr.Proofs.PlanarDegree
import Spqr.Proofs.PlanarCompCons
/-!
# Genus inequality (PROOF.md §8.6)
Every rotation system of `es` (an `IsEmbedding`) satisfies `F + 2V ≤ 2(2C + E)`, with equality
exactly for the planar ones.  Proof by deleting the first edge `p` of `p :: es`: its four
quarter-edges are detached from their vertex rotations by `conj`s (each pairs `p`'s two quarter-edges
of one end with each other and the two former neighbours with each other), after which the system is
`single.union ρ` (or `loop.union ρ`) for the rotation system `ρ = shrink _` of `es`.  Each `conj`
changes the vertex-orbit count by `+2` and the face-orbit count by `±2`; the component/non-isolated
bookkeeping is done on `ccCount` through its component minima.
-/
namespace Spqr
open Classical

theorem sameOrbit_step2 {f : Nat → Nat} {x y z : Nat} (h1 : f x = y) (h2 : f y = z) :
    SameOrbit f x z :=
  (RotationSystem.sameOrbit_of_eq h1).trans (RotationSystem.sameOrbit_of_eq h2)

/-- `swapImg f x y` with `f y = x` is `f` with `x` skipped. -/
theorem swapImg_eq_skip {f : Nat → Nat} {S : Finset Nat} (hf : IsPermOn f S) {x y : Nat}
    (hy : y ∈ S) (hfy : f y = x) : swapImg f x y = skip f [x] := by
  funext a
  unfold swapImg skip
  by_cases hax : a = x
  · simp [hax, hfy]
  · simp only [List.mem_singleton, hax, ↓reduceIte]
    by_cases hay : a = y
    · simp [hay, hfy]
    · have : f a ≠ x := fun h => hay (hf.inj a (by
        by_contra ha
        exact hax (by rw [← hf.fix a ha, h])) y hy (h.trans hfy.symm))
      simp [hay, this]

theorem sameOrbit_swapImg_of_fix {f : Nat → Nat} {S : Finset Nat} (hf : IsPermOn f S) {x y : Nat}
    (hx : x ∈ S) (hy : y ∈ S) (hxy : x ≠ y) (hfy : f y = x) {a b : Nat} (ha : a ≠ x) (hb : b ≠ x)
    (h : SameOrbit f a b) : SameOrbit (swapImg f x y) a b := by
  rw [swapImg_eq_skip hf hy hfy]
  refine (sameOrbit_skip [x] ?_ (List.nodup_singleton x) (by simpa using ha)).2
    ⟨h, by simpa using hb⟩
  intro p hp
  simp only [List.mem_singleton] at hp ⊢
  subst hp
  intro h'
  exact hxy (hf.inj p hx y hy (h'.trans hfy.symm))

theorem swapImg_comm (f : Nat → Nat) (x y : Nat) : swapImg f x y = swapImg f y x := by
  funext a
  unfold swapImg
  by_cases hx : a = x <;> by_cases hy : a = y <;> simp [hx, hy] <;> intro h <;> subst h <;> rfl

theorem sameOrbit_mem_pair {f : Nat → Nat} {x y : Nat} (h1 : f x = y) (h2 : f y = x) {w : Nat}
    (h : Reach f x w) : w = x ∨ w = y := by
  obtain ⟨k, rfl⟩ := IsPermOn.reach_iff.1 h
  clear h
  induction k with
  | zero => exact Or.inl rfl
  | succ k ih =>
    rw [Function.iterate_succ_apply']
    rcases ih with h | h <;> rw [h] <;> simp [h1, h2]

namespace RotationSystem

theorem conj_rot (rs : RotationSystem) (a b : Nat) (ht : rs.Total) (ha : a < rs.size)
    (hb : b < rs.size) {q : Nat} (hq : q < rs.size) :
    (rs.conj a b).rot q = Equiv.swap a b (rs.rot (Equiv.swap a b q)) := by
  unfold rot
  rw [conj_get_lt rs a b hq, get_eq_rot ht (swap_lt rs a b ha hb hq)]
  rfl

/-- Drop the first four quarter-edges (the first edge). -/
def shrink (rs : RotationSystem) : RotationSystem :=
  ⟨(rs.rotAdj.extract 4 rs.size).map (Option.map (· - 4))⟩

theorem shrink_size (rs : RotationSystem) : rs.shrink.size = rs.size - 4 := by
  simp [shrink, size]

theorem shrink_get (rs : RotationSystem) (q : Nat) :
    rs.shrink.get q = (rs.get (4 + q)).map (· - 4) := by
  unfold shrink get
  rw [Array.getElem?_map, Array.getElem?_extract]
  split_ifs with h
  · rw [Nat.add_comm]
    generalize rs.rotAdj[q + 4]? = o
    cases o with
    | none => rfl
    | some o => cases o <;> rfl
  · rw [Array.getElem?_eq_none (by simp [size] at h; omega)]
    rfl

theorem ext_get {rs rs' : RotationSystem} (hs : rs.size = rs'.size) (ht : rs.Total)
    (ht' : rs'.Total) (h : ∀ q, rs.get q = rs'.get q) : rs = rs' := by
  obtain ⟨xs⟩ := rs
  obtain ⟨ys⟩ := rs'
  congr
  apply Array.ext hs
  intro i hi hi'
  have h1 := ht i hi
  have h2 := ht' i hi'
  have := h i
  simp only [get, Array.getElem?_eq_getElem hi, Array.getElem?_eq_getElem hi',
    Option.bind_some, id] at h1 h2 this
  exact this

theorem eq_union_shrink {K rs : RotationSystem} (hK : K.size = 4) (hKt : K.Total) (ht : rs.Total)
    (hs : 4 ≤ rs.size) (hlow : ∀ q, q < 4 → rs.get q = K.get q)
    (hhigh : ∀ q, 4 ≤ q → q < rs.size → 4 ≤ rs.rot q) : rs = K.union rs.shrink := by
  have hst : rs.shrink.Total := by
    intro q hq
    rw [shrink_size] at hq
    rw [shrink_get, Option.isSome_map]
    exact ht (4 + q) (by omega)
  refine ext_get ?_ ht (union_total _ _ hKt hst) fun q => ?_
  · rw [union_size, hK, shrink_size]; omega
  by_cases hq : q < 4
  · rw [union_get_lt _ _ (hK ▸ hq), hlow q hq]
  · rw [union_get_ge _ _ (hK ▸ not_lt.1 hq), hK, shrink_get, Nat.add_sub_cancel' (not_lt.1 hq)]
    by_cases hq' : q < rs.size
    · rw [get_eq_rot ht hq']
      have := hhigh q (not_lt.1 hq) hq'
      simp only [Option.map_some]
      congr 1; omega
    · rw [get_eq_none (not_lt.1 hq')]; rfl

theorem map_add_eq_some {o : Option Nat} {k q : Nat} (h : o.map (k + ·) = some (k + q)) :
    o = some q := by
  cases o with
  | none => simp at h
  | some r => simp only [Option.map_some, Option.some.injEq] at h ⊢; omega

/-- The second summand of a disjoint union with a one-edge system is a rotation system of the tail. -/
theorem isEmbedding_of_union {K ρ : RotationSystem} {p : Nat × Nat} {es : List (Nat × Nat)}
    {n : Nat} (hK : K.size = 4) (hs : (K.union ρ).size = 4 * (p :: es).length)
    (hverts : ∀ q ∈ es, q.1 < n ∧ q.2 < n) (ht : (K.union ρ).Total)
    (hi : (K.union ρ).Involution) (ho : (K.union ρ).OppositeDir)
    (hsv : (K.union ρ).SameVertex (p :: es))
    (hV : ρ.numVertexOrbits = 2 * numNonIsolated es n) : IsEmbedding es n ρ where
  size := by rw [union_size, hK] at hs; simp at hs; omega
  verts := hverts
  total q hq := by
    have := ht (4 + q) (by rw [union_size, hK]; omega)
    rwa [union_get_ge _ _ (by omega), hK, Nat.add_sub_cancel_left, Option.isSome_map] at this
  involution q hq r hr := by
    have h1 := hi (4 + q) (by rw [union_size, hK]; omega) (4 + r)
      (by rw [union_get_ge _ _ (by omega), hK, Nat.add_sub_cancel_left]; exact Option.mem_map_of_mem _ hr)
    rw [union_size, hK] at h1
    refine ⟨by omega, ?_⟩
    have h2 := h1.2
    rw [union_get_ge _ _ (by omega), hK, Nat.add_sub_cancel_left] at h2
    exact map_add_eq_some h2
  opposite_dir q hq r hr := by
    have h1 := ho (4 + q) (by rw [union_size, hK]; omega) (4 + r)
      (by rw [union_get_ge _ _ (by omega), hK, Nat.add_sub_cancel_left]; exact Option.mem_map_of_mem _ hr)
    unfold QE.dir at h1 ⊢
    omega
  same_vertex q hq r hr := by
    have h1 := hsv (4 + q) (by rw [union_size, hK]; omega) (4 + r)
      (by rw [union_get_ge _ _ (by omega), hK, Nat.add_sub_cancel_left]; exact Option.mem_map_of_mem _ hr)
    rwa [vert_cons_add, vert_cons_add] at h1
  vertex_orbits := hV

/-- Detaching the end `{s, t}` (two quarter-edges of the same edge and vertex, opposite directions)
of the first edge from its vertex rotation, when both their rotation partners are outside the first
edge: `rs.conj s (rs.rot t)` pairs `s ↔ t` and the two former partners with each other. -/
structure DetachAt (rs : RotationSystem) (s t : Nat) : Prop where
  total : rs.Total
  involution : rs.Involution
  opposite_dir : rs.OppositeDir
  size4 : 4 ≤ rs.size
  s_lt : s < 4
  t_lt : t < 4
  par : s % 2 ≠ t % 2
  y_ge : 4 ≤ rs.rot t
  z_ge : 4 ≤ rs.rot s

namespace DetachAt

variable {rs : RotationSystem} {s t : Nat} (d : DetachAt rs s t)
include d

theorem s_size : s < rs.size := by have := d.size4; have := d.s_lt; omega
theorem t_size : t < rs.size := by have := d.size4; have := d.t_lt; omega
theorem y_size : rs.rot t < rs.size := rot_lt d.total d.involution d.t_size
theorem z_size : rs.rot s < rs.size := rot_lt d.total d.involution d.s_size
theorem rot_y : rs.rot (rs.rot t) = t := rot_rot d.total d.involution d.t_size
theorem rot_z : rs.rot (rs.rot s) = s := rot_rot d.total d.involution d.s_size
theorem s_ne_t : s ≠ t := fun h => d.par (by rw [h])
theorem y_ne_z : rs.rot t ≠ rs.rot s := fun h => d.s_ne_t (rot_inj d.total d.involution d.s_size d.t_size h.symm)
theorem s_ne_y : s ≠ rs.rot t := by have := d.y_ge; have := d.s_lt; omega
theorem par_y : s % 2 = rs.rot t % 2 := by
  have h1 := rot_mod2 d.total d.opposite_dir d.t_size
  have := d.par; omega
theorem par_z : rs.rot s % 2 = t % 2 := by
  have h1 := rot_mod2 d.total d.opposite_dir d.s_size
  have := d.par; omega

/-- The rotation of the detached system. -/
theorem rot_eq (q : Nat) (hq : q < rs.size) :
    (rs.conj s (rs.rot t)).rot q =
      if q = t then s else if q = s then t else
      if q = rs.rot t then rs.rot s else if q = rs.rot s then rs.rot t else rs.rot q := by
  have hy := d.y_ge; have hz := d.z_ge; have hs := d.s_lt; have ht := d.t_lt
  have hyz := d.y_ne_z
  rw [conj_rot rs s (rs.rot t) d.total d.s_size d.y_size hq]
  split_ifs with h1 h2 h3 h4
  · subst h1
    rw [Equiv.swap_apply_of_ne_of_ne d.s_ne_t.symm (by omega), Equiv.swap_apply_right]
  · subst h2
    rw [Equiv.swap_apply_left, d.rot_y, Equiv.swap_apply_of_ne_of_ne d.s_ne_t.symm (by omega)]
  · subst h3
    rw [Equiv.swap_apply_right, Equiv.swap_apply_of_ne_of_ne (by omega) hyz.symm]
  · subst h4
    rw [Equiv.swap_apply_of_ne_of_ne (by omega) hyz.symm, d.rot_z, Equiv.swap_apply_left]
  · rw [Equiv.swap_apply_of_ne_of_ne h2 h3]
    refine Equiv.swap_apply_of_ne_of_ne ?_ ?_
    · intro h; exact h4 (by rw [← rot_rot d.total d.involution hq, h])
    · intro h; exact h1 (by rw [← rot_rot d.total d.involution hq, h, d.rot_y])

theorem conj_total' : (rs.conj s (rs.rot t)).Total := conj_total rs _ _ d.total d.s_size d.y_size
theorem conj_involution' : (rs.conj s (rs.rot t)).Involution :=
  conj_involution rs _ _ d.total d.involution d.s_size d.y_size
theorem conj_oppositeDir' : (rs.conj s (rs.rot t)).OppositeDir :=
  conj_oppositeDir rs _ _ d.total d.opposite_dir d.s_size d.y_size d.par_y
theorem conj_size' : (rs.conj s (rs.rot t)).size = rs.size := conj_size _ _ _

theorem conj_sameVertex' {es : List (Nat × Nat)} (hsv : rs.SameVertex es)
    (hv : QE.vert es s = QE.vert es t) {u : Nat} (hu : QE.vert es t = some u) :
    (rs.conj s (rs.rot t)).SameVertex es := by
  have := conj_sameVertex rs s (rs.rot t) d.total d.s_size d.y_size hsv (hv.trans hu)
    ((vert_rot d.total hsv d.t_size).trans hu)
  rwa [identEdges_self] at this

theorem numVertexOrbits_eq {m : Nat} (hm : rs.size = 4 * m)
    (h1 : SameOrbit (rs.stepC 1) (s ^^^ 1) (rs.rot t ^^^ 1))
    (h2 : SameOrbit (rs.stepC 1) (rs.rot s ^^^ 1) (t ^^^ 1)) :
    rs.numVertexOrbits + 2 = (rs.conj s (rs.rot t)).numVertexOrbits := by
  rw [RotationSystem.numVertexOrbits_eq _ d.conj_total' d.conj_involution' (m := m)
      (by rw [conj_size, hm]), RotationSystem.numVertexOrbits_eq _ d.total d.involution hm, ← hm]
  exact conj_orbitCount_split rs s (rs.rot t) d.total d.involution d.opposite_dir d.s_size
    d.y_size d.s_ne_y d.y_ne_z.symm hm (by decide) (by decide) h1 (by rw [d.rot_y]; exact h2)

/-- In an embedding, both pairs of quarter-edges exchanged by the detach lie on common vertex
orbits (same vertex, same direction). -/
theorem vertexOrbits_of_embedding {es : List (Nat × Nat)} {n : Nat} (h : IsEmbedding es n rs)
    (hv : QE.vert es s = QE.vert es t) :
    SameOrbit (rs.stepC 1) (s ^^^ 1) (rs.rot t ^^^ 1) ∧
      SameOrbit (rs.stepC 1) (rs.rot s ^^^ 1) (t ^^^ 1) := by
  have hlt : ∀ q, q < rs.size → q ^^^ 1 < rs.size := fun q hq =>
    h.size ▸ xor_lt_mul4 (h.size ▸ hq) (by decide)
  have hvy := vert_rot h.total h.same_vertex d.t_size
  have hvz := vert_rot h.total h.same_vertex d.s_size
  have p1 := xor_mod2 s (c := 1) (by decide)
  have p2 := xor_mod2 t (c := 1) (by decide)
  have p3 := xor_mod2 (rs.rot t) (c := 1) (by decide)
  have p4 := xor_mod2 (rs.rot s) (c := 1) (by decide)
  have py := d.par_y
  have pz := d.par_z
  constructor
  · refine h.sameOrbit_of_lbl (hlt _ d.s_size) (hlt _ d.y_size) ?_
    unfold lbl
    rw [vert_xor_one, vert_xor_one, hvy, hv]
    congr 1; omega
  · refine h.sameOrbit_of_lbl (hlt _ d.z_size) (hlt _ d.t_size) ?_
    unfold lbl
    rw [vert_xor_one, vert_xor_one, hvz, hv]
    congr 1; omega

/-- Vertex-orbit facts about points other than `s`, `t` survive a detach of a direction pair
`s = t ^^^ 1`: the detach only skips `s` and `t` out of their vertex orbits. -/
theorem sameOrbit_stepC1_conj {m : Nat} (hm : rs.size = 4 * m) (hst' : s = t ^^^ 1) {a b : Nat}
    (ha : a ≠ s) (ha' : a ≠ t) (hb : b ≠ s) (hb' : b ≠ t) (h : SameOrbit (rs.stepC 1) a b) :
    SameOrbit ((rs.conj s (rs.rot t)).stepC 1) a b := by
  rw [conj_stepC rs s (rs.rot t) d.total d.involution d.opposite_dir d.s_size d.y_size d.s_ne_y
    d.y_ne_z.symm hm (by decide), d.rot_y]
  have hf := isPermOn_stepC d.total d.involution hm (c := 1) (by decide)
  have hts : t ^^^ 1 = s := by rw [hst']
  have hx : s ^^^ 1 = t := by rw [hst', xor_xor_self]
  have mem : ∀ q, q < rs.size → q ^^^ 1 ∈ Finset.range rs.size := fun q hq =>
    Finset.mem_range.2 (hm ▸ xor_lt_mul4 (hm ▸ hq) (by decide))
  have hfy : rs.stepC 1 (rs.rot t ^^^ 1) = s ^^^ 1 := by
    rw [stepC_eq_rot d.total hm (by decide) (hm ▸ xor_lt_mul4 (hm ▸ d.y_size) (by decide)),
      xor_xor_self, d.rot_y, hx]
  have h1 := sameOrbit_swapImg_of_fix hf (mem _ d.s_size) (mem _ d.y_size)
    (fun h => d.s_ne_y (xor_right_inj.1 h)) hfy (by rw [hx]; exact ha') (by rw [hx]; exact hb') h
  have hg : IsPermOn (swapImg (rs.stepC 1) (s ^^^ 1) (rs.rot t ^^^ 1)) (Finset.range rs.size) :=
    isPermOn_swapImg hf (mem _ d.s_size) (mem _ d.y_size)
  have hzs : rs.rot s ^^^ 1 ≠ s ^^^ 1 := fun h =>
    absurd (xor_right_inj.1 h) (by have := d.z_ge; have := d.s_lt; omega)
  have hzy : rs.rot s ^^^ 1 ≠ rs.rot t ^^^ 1 := fun h => d.y_ne_z (xor_right_inj.1 h).symm
  have hgz : swapImg (rs.stepC 1) (s ^^^ 1) (rs.rot t ^^^ 1) (rs.rot s ^^^ 1) = t ^^^ 1 := by
    unfold swapImg
    rw [ite_eq_right hzs, ite_eq_right hzy,
      stepC_eq_rot d.total hm (by decide) (hm ▸ xor_lt_mul4 (hm ▸ d.z_size) (by decide)),
      xor_xor_self, d.rot_z, hts]
  rw [swapImg_comm]
  exact sameOrbit_swapImg_of_fix hg (mem _ d.t_size) (mem _ d.z_size)
    (fun h => by have := xor_right_inj.1 h; have := d.z_ge; have := d.t_lt; omega) hgz (by rw [hts]; exact ha) (by rw [hts]; exact hb) h1

theorem stepC3_y {m : Nat} (hm : rs.size = 4 * m) : rs.stepC 3 (rs.rot t ^^^ 3) = t := by
  rw [stepC_eq_rot d.total hm (by decide) (hm ▸ xor_lt_mul4 (hm ▸ d.y_size) (by decide)),
    xor_xor_self, d.rot_y]

theorem stepC3_z {m : Nat} (hm : rs.size = 4 * m) : rs.stepC 3 (rs.rot s ^^^ 3) = s := by
  rw [stepC_eq_rot d.total hm (by decide) (hm ▸ xor_lt_mul4 (hm ▸ d.z_size) (by decide)),
    xor_xor_self, d.rot_z]

theorem par_sameOrbit {m : Nat} (hm : rs.size = 4 * m) {x y : Nat}
    (h : SameOrbit (rs.stepC 3) x y) : x % 2 = y % 2 :=
  sameOrbit_invariant (· % 2) (stepC_mod2 d.total d.opposite_dir hm (by decide) (by decide)) h

/-- Face orbits when the two sides of the detached end are cofacial: two orbits split off. -/
theorem numFaceOrbits_split {m : Nat} (hm : rs.size = 4 * m)
    (hside : SameOrbit (rs.stepC 3) (s ^^^ 3) t) :
    rs.numFaceOrbits + 2 = (rs.conj s (rs.rot t)).numFaceOrbits := by
  refine conj_numFaceOrbits_split rs s (rs.rot t) d.total d.involution d.opposite_dir d.s_size
    d.y_size d.s_ne_y d.y_ne_z.symm hm ?_ ?_
  · exact hside.trans (sameOrbit_of_eq (d.stepC3_y hm)).symm
  · rw [d.rot_y]
    have := sameOrbit_stepC_xor d.total d.involution hm (by decide) hside
    rw [xor_xor_self] at this
    exact (sameOrbit_of_eq (d.stepC3_z hm)).trans this

/-- Face orbits when the two sides of the detached end are not cofacial: two orbits merge. -/
theorem numFaceOrbits_merge {m : Nat} (hm : rs.size = 4 * m)
    (hside : ¬SameOrbit (rs.stepC 3) (s ^^^ 3) t) :
    (rs.conj s (rs.rot t)).numFaceOrbits + 2 = rs.numFaceOrbits := by
  have hmir : ∀ {x y}, SameOrbit (rs.stepC 3) x y →
      SameOrbit (rs.stepC 3) (x ^^^ 3) (y ^^^ 3) :=
    fun h => sameOrbit_stepC_xor d.total d.involution hm (by decide) h
  have hy := d.stepC3_y hm
  have hz := d.stepC3_z hm
  refine conj_numFaceOrbits rs s (rs.rot t) d.total d.involution d.opposite_dir d.s_size d.y_size
    d.s_ne_y d.y_ne_z.symm hm
    (fun w => SameOrbit (rs.stepC 3) w (s ^^^ 3) ∨ SameOrbit (rs.stepC 3) w (rs.rot s ^^^ 3))
    ?_ (Or.inl (SameOrbit.refl _ _)) ?_ (Or.inr (SameOrbit.refl _ _)) ?_
  · intro w
    have h := sameOrbit_of_eq (f := rs.stepC 3) (q := w) rfl
    exact or_congr ⟨fun h' => h.symm.trans h', fun h' => h.trans h'⟩
      ⟨fun h' => h.symm.trans h', fun h' => h.trans h'⟩
  · rintro (h | h)
    · exact hside ((sameOrbit_of_eq hy).symm.trans h).symm
    · have h1 := (sameOrbit_of_eq hy).symm.trans (h.trans (sameOrbit_of_eq hz))
      exact d.par (d.par_sameOrbit hm h1).symm
  · rw [d.rot_y]
    rintro (h | h)
    · have := d.par_sameOrbit hm h
      have := xor_mod2 t (c := 3) (by decide)
      have := xor_mod2 s (c := 3) (by decide)
      have := d.par
      omega
    · have h1 := hmir (h.trans (sameOrbit_of_eq hz))
      rw [xor_xor_self] at h1
      exact hside h1.symm

end DetachAt

/-- The edge of a quarter-edge `w ≥ 4` of `p :: es` lies in `es` and joins the vertices of `w` and
of its across quarter-edge. -/
theorem edgesConn_tail_across {es : List (Nat × Nat)} {p : Nat × Nat} {w a b : Nat} (hw : 4 ≤ w)
    (ha : QE.vert (p :: es) w = some a) (hb : QE.vert (p :: es) (w ^^^ 3) = some b) :
    EdgesConn es a b := by
  rw [vert_across] at hb
  unfold QE.vert QE.edge QE.side at ha
  have hidx : (p :: es)[w / 4]? = es[w / 4 - 1]? := by
    obtain ⟨k, hk⟩ : ∃ k, w / 4 = k + 1 := ⟨w / 4 - 1, by omega⟩
    rw [hk, List.getElem?_cons_succ]; simp
  rw [hidx] at ha hb
  obtain ⟨q, hq⟩ : ∃ q, es[w / 4 - 1]? = some q := by
    cases h : es[w / 4 - 1]? with
    | none => simp [h] at ha
    | some q => exact ⟨q, rfl⟩
  have hmem : q ∈ es := List.mem_of_getElem? hq
  rw [hq] at ha hb
  simp only [Option.map_some, Option.some.injEq] at ha hb
  subst ha hb
  split_ifs
  · exact Relation.ReflTransGen.single (Or.inl hmem)
  · exact Relation.ReflTransGen.single (Or.inr hmem)

/-- If the two sides `0`, `2` of the first edge are not cofacial, the face walk from `0` returns
to `0` through the tail only, so the endpoints of the first edge are connected in the tail. -/
theorem edgesConn_of_not_cofacial {es : List (Nat × Nat)} {p : Nat × Nat} {n : Nat}
    {rs : RotationSystem} (h : IsEmbedding (p :: es) n rs) (h0 : 4 ≤ rs.rot 0) (h3 : 4 ≤ rs.rot 3)
    (hside : ¬SameOrbit (rs.stepC 3) 2 0) : EdgesConn es p.1 p.2 := by
  have hm : rs.size = 4 * (es.length + 1) := by rw [h.size]; rfl
  have hsz : 4 ≤ rs.size := by omega
  set f := rs.stepC 3 with hf
  have hstep : ∀ q, q < rs.size → f q = rs.rot (q ^^^ 3) := fun q hq =>
    stepC_eq_rot h.total hm (by decide) hq
  have hpar : ∀ w, SameOrbit f 0 w → w % 2 = 0 := fun w hw =>
    (sameOrbit_invariant (· % 2)
      (stepC_mod2 h.total h.opposite_dir hm (by decide) (by decide)) hw).symm
  have horb : ∀ k, SameOrbit f 0 (f^[k] 0) := by
    intro k
    induction k with
    | zero => exact SameOrbit.refl _ _
    | succ k ih => rw [Function.iterate_succ_apply']; exact ih.trans (sameOrbit_of_eq rfl)
  have hlt : ∀ k, f^[k] 0 < rs.size := by
    intro k
    induction k with
    | zero => simp only [Function.iterate_zero, id]; omega
    | succ k ih =>
      rw [Function.iterate_succ_apply']
      exact stepC_lt h.total h.involution hm (by decide) ih
  have hv3 : QE.vert (p :: es) 3 = some p.2 := by rw [vert_cons_lt _ _ (by decide)]; rfl
  have hv0 : QE.vert (p :: es) 0 = some p.1 := by rw [vert_cons_lt _ _ (by decide)]; rfl
  have claim : ∀ k, f^[k] 0 = 0 ∨ (4 ≤ f^[k] 0 ∧
      ∀ vw, QE.vert (p :: es) (f^[k] 0) = some vw → EdgesConn es p.2 vw) := by
    intro k
    induction k with
    | zero => exact Or.inl rfl
    | succ k ih =>
      have hpar' := hpar _ (horb (k + 1))
      have hne2 : f^[k + 1] 0 ≠ 2 := fun h2 => hside (by
        have := horb (k + 1); rw [h2] at this; exact this.symm)
      rw [Function.iterate_succ_apply'] at hpar' hne2 ⊢
      rcases ih with hz | ⟨hw, hconn⟩
      · rw [hz, hstep 0 (by omega)]
        have h03 : (0 : Nat) ^^^ 3 = 3 := by decide
        rw [h03]
        refine Or.inr ⟨h3, fun vw hvw => ?_⟩
        rw [vert_rot h.total h.same_vertex (by omega), hv3] at hvw
        cases hvw
        exact Relation.ReflTransGen.refl
      · by_cases hz : f (f^[k] 0) = 0
        · exact Or.inl hz
        · refine Or.inr ⟨by omega, fun vw hvw => ?_⟩
          rw [hstep _ (hlt k),
            vert_rot h.total h.same_vertex (hm ▸ xor_lt_mul4 (hm ▸ hlt k) (by decide))] at hvw
          obtain ⟨va, hva⟩ := vert_some_of_lt (es := p :: es) (q := f^[k] 0)
            (by rw [← h.size]; exact hlt k)
          exact (hconn va hva).trans (edgesConn_tail_across hw hva hvw)
  have hr0 : rs.rot 0 < rs.size := rot_lt h.total h.involution (by omega)
  have hr03 : rs.rot 0 ^^^ 3 < rs.size := hm ▸ xor_lt_mul4 (hm ▸ hr0) (by decide)
  have hpred : f (rs.rot 0 ^^^ 3) = 0 := by
    rw [hstep _ hr03, xor_xor_self, rot_rot h.total h.involution (by omega)]
  have hperm := isPermOn_stepC h.total h.involution hm (c := 3) (by decide)
  obtain ⟨k, hk⟩ := IsPermOn.reach_iff.1
    (hperm.reach_of_sameOrbit (Finset.mem_range.2 (by omega)) (sameOrbit_of_eq hpred).symm)
  have hge : 4 * 1 ≤ rs.rot 0 ^^^ 3 := ge_of_xor_ge (by omega) (by decide)
  rcases claim k with hc | ⟨_, hconn⟩
  · rw [hk] at hc; omega
  · rw [hk] at hconn
    obtain ⟨vb, hvb⟩ := vert_some_of_lt (es := p :: es) (q := rs.rot 0 ^^^ 3)
      (by rw [← h.size]; exact hr03)
    have hva : QE.vert (p :: es) (rs.rot 0) = some p.1 := by
      rw [vert_rot h.total h.same_vertex (by omega), hv0]
    exact (edgesConn_tail_across h0 hva hvb).trans (edgesConn_symm (hconn vb hvb))

/-! ### Deleting the first edge -/

theorem loop_total : loop.Total := by decide
theorem loop_involution : loop.Involution := by decide
theorem loop_numFaceOrbits : loop.numFaceOrbits = 4 := by decide
theorem loop_numVertexOrbits : loop.numVertexOrbits = 2 := by decide
theorem single_rot4 :
    single.rot 0 = 1 ∧ single.rot 1 = 0 ∧ single.rot 2 = 3 ∧ single.rot 3 = 2 := by decide
theorem loop_rot4 : loop.rot 0 = 3 ∧ loop.rot 1 = 2 ∧ loop.rot 2 = 1 ∧ loop.rot 3 = 0 := by decide

theorem shrink_total {H : RotationSystem} (ht : H.Total) : H.shrink.Total := by
  intro q hq
  rw [shrink_size] at hq
  rw [shrink_get, Option.isSome_map]
  exact ht (4 + q) (by omega)

theorem shrink_involution {H : RotationSystem} (ht : H.Total) (hi : H.Involution)
    (hhigh : ∀ q, 4 ≤ q → q < H.size → 4 ≤ H.rot q) : H.shrink.Involution := by
  intro q hq r hr
  rw [shrink_size] at hq
  rw [Option.mem_def, shrink_get, get_eq_rot ht (by omega)] at hr
  simp only [Option.map_some, Option.some.injEq] at hr
  have hge := hhigh (4 + q) (by omega) (by omega)
  have hlt := rot_lt ht hi (q := 4 + q) (by omega)
  have h4r : 4 + r = H.rot (4 + q) := by omega
  refine ⟨by rw [shrink_size]; omega, ?_⟩
  rw [shrink_get, h4r, get_eq_rot ht hlt, rot_rot ht hi (by omega)]
  simp

theorem exists_qe_of_hasEdge {es : List (Nat × Nat)} {u : Nat} (h : HasEdge es u) {d : Nat}
    (hd : d < 2) : ∃ q, q < 4 * es.length ∧ QE.vert es q = some u ∧ q % 2 = d := by
  obtain ⟨e, he, hu⟩ := h
  obtain ⟨i, hi, hei⟩ := List.mem_iff_getElem.1 he
  rcases hu with hu | hu
  · refine ⟨4 * i + d, by omega, ?_, by omega⟩
    unfold QE.vert QE.edge QE.side
    rw [show (4 * i + d) / 4 = i by omega, show (4 * i + d) / 2 % 2 = 0 by omega,
      List.getElem?_eq_getElem hi, hei]
    simp [hu]
  · refine ⟨4 * i + 2 + d, by omega, ?_, by omega⟩
    unfold QE.vert QE.edge QE.side
    rw [show (4 * i + 2 + d) / 4 = i by omega, show (4 * i + 2 + d) / 2 % 2 = 1 by omega,
      List.getElem?_eq_getElem hi, hei]
    simp [hu]

/-- A vertex orbit `{a, b}` of size at most two contains every quarter-edge with the label of
`a`. -/
theorem mem_pair_of_lbl {es : List (Nat × Nat)} {n : Nat} {rs : RotationSystem}
    (h : IsEmbedding es n rs) {a b w : Nat} (ha : a < rs.size) (hw : w < rs.size)
    (hab : rs.stepC 1 a = b) (hba : rs.stepC 1 b = a) (hl : lbl es w = lbl es a) :
    w = a ∨ w = b :=
  sameOrbit_mem_pair hab hba
    ((isPermOn_stepC h.total h.involution h.size (c := 1) (by decide)).reach_of_sameOrbit
      (Finset.mem_range.2 ha) (h.sameOrbit_of_lbl ha hw hl.symm))

/-- If the vertex orbit of a quarter-edge of the first edge stays inside the first edge, its
vertex has no other edge. -/
theorem not_hasEdge_of_pair {es : List (Nat × Nat)} {p : Nat × Nat} {n : Nat}
    {rs : RotationSystem} (h : IsEmbedding (p :: es) n rs) {a b u : Nat} (ha : a < 4) (hb : b < 4)
    (hab : rs.stepC 1 a = b) (hba : rs.stepC 1 b = a) (hu : QE.vert (p :: es) a = some u) :
    ¬HasEdge es u := by
  intro he
  obtain ⟨q, hq, hv, hd⟩ := exists_qe_of_hasEdge he (Nat.mod_lt a (by omega) : a % 2 < 2)
  have hm : rs.size = 4 * (es.length + 1) := by rw [h.size]; rfl
  have hl : lbl (p :: es) (4 + q) = lbl (p :: es) a := by
    unfold lbl
    rw [vert_cons_add, hv, hu]
    simp only [Option.getD_some, Prod.mk.injEq, true_and]
    omega
  rcases mem_pair_of_lbl h (by omega) (by omega) hab hba hl with h1 | h1 <;> omega

theorem hasEdge_of_vert_ge {es : List (Nat × Nat)} {p : Nat × Nat} {w u : Nat} (hw : 4 ≤ w)
    (hv : QE.vert (p :: es) w = some u) : HasEdge es u := by
  obtain ⟨q, rfl⟩ : ∃ q, w = 4 + q := ⟨w - 4, by omega⟩
  rw [vert_cons_add] at hv
  exact hasEdge_of_vert hv

/-- After detaching, `H` is `K.union H.shrink` with `H.shrink` an embedding of the tail; the
inductive bound for the tail transfers. -/
theorem genus_cons_bound {es : List (Nat × Nat)} {p : Nat × Nat} {n : Nat}
    {H K : RotationSystem}
    (IH : ∀ ρ, IsEmbedding es n ρ →
      ρ.numFaceOrbits + 2 * numNonIsolated es n ≤ 2 * (2 * numComponents es n + es.length))
    (hverts : ∀ q ∈ es, q.1 < n ∧ q.2 < n)
    (hK4 : K.size = 4) (hKt : K.Total) (hKi : K.Involution)
    (hHt : H.Total) (hHi : H.Involution) (hHo : H.OppositeDir) (hHsv : H.SameVertex (p :: es))
    (hHs : H.size = 4 * (p :: es).length) (hrot : ∀ q, q < 4 → H.rot q = K.rot q)
    (hV : H.numVertexOrbits = K.numVertexOrbits + 2 * numNonIsolated es n) :
    H.numFaceOrbits + 2 * numNonIsolated es n ≤
      K.numFaceOrbits + 2 * (2 * numComponents es n + es.length) := by
  have hs4 : 4 ≤ H.size := by rw [hHs, List.length_cons]; omega
  have hlow : ∀ q, q < 4 → H.get q = K.get q := fun q hq => by
    rw [get_eq_rot hHt (by omega), get_eq_rot hKt (by omega), hrot q hq]
  have hhigh : ∀ q, 4 ≤ q → q < H.size → 4 ≤ H.rot q := by
    intro q hq4 hq
    by_contra hlt
    push_neg at hlt
    have h1 := rot_rot hHt hHi hq
    have h2 : H.rot (H.rot q) < 4 := by
      rw [hrot _ hlt, ← hK4]; exact rot_lt hKt hKi (by omega)
    omega
  have heq : H = K.union H.shrink := eq_union_shrink hK4 hKt hHt hs4 hlow hhigh
  have hss : H.shrink.size = 4 * es.length := by rw [shrink_size, hHs, List.length_cons]; omega
  have hst := shrink_total hHt
  have hsi := shrink_involution hHt hHi hhigh
  have hVu := union_numVertexOrbits K H.shrink hKt hKi hst hsi (m₁ := 1) hK4 hss
  rw [← heq] at hVu
  have hFu := union_numFaceOrbits K H.shrink hKt hKi hst hsi (m₁ := 1) hK4 hss
  rw [← heq] at hFu
  have hemb : IsEmbedding es n H.shrink :=
    isEmbedding_of_union hK4 (by rw [← heq]; exact hHs) hverts (by rw [← heq]; exact hHt)
      (by rw [← heq]; exact hHi) (by rw [← heq]; exact hHo) (by rw [← heq]; exact hHsv)
      (by omega)
  have := IH _ hemb
  omega

/-- One detach (`s`, `t` quarter-edges of the first edge at a common vertex, cofacial sides):
the result pairs `s ↔ t`, keeps the other small rotations, and gains two vertex and two face
orbits. -/
theorem detach_one {es : List (Nat × Nat)} {p : Nat × Nat} {n : Nat} {rs : RotationSystem}
    {s t : Nat} (h : IsEmbedding (p :: es) n rs) (d : DetachAt rs s t)
    (hv : QE.vert (p :: es) s = QE.vert (p :: es) t) {u : Nat} (hu : QE.vert (p :: es) t = some u)
    (hside : SameOrbit (rs.stepC 3) (s ^^^ 3) t) :
    ∃ H : RotationSystem, H.Total ∧ H.Involution ∧ H.OppositeDir ∧ H.SameVertex (p :: es) ∧
      H.size = rs.size ∧ H.rot t = s ∧ H.rot s = t ∧
      (∀ q, q < 4 → q ≠ s → q ≠ t → H.rot q = rs.rot q) ∧
      H.numVertexOrbits = rs.numVertexOrbits + 2 ∧ H.numFaceOrbits = rs.numFaceOrbits + 2 := by
  have hm : rs.size = 4 * (es.length + 1) := by rw [h.size]; rfl
  obtain ⟨h1, h2⟩ := d.vertexOrbits_of_embedding h hv
  refine ⟨rs.conj s (rs.rot t), d.conj_total', d.conj_involution', d.conj_oppositeDir',
    d.conj_sameVertex' h.same_vertex hv hu, d.conj_size', ?_, ?_, ?_,
    (d.numVertexOrbits_eq hm h1 h2).symm, (d.numFaceOrbits_split hm hside).symm⟩
  · rw [d.rot_eq t d.t_size]; simp
  · rw [d.rot_eq s d.s_size]; simp [d.s_ne_t]
  · intro q hq hqs hqt
    rw [d.rot_eq q (by have := d.size4; omega)]
    have := d.y_ge
    have := d.z_ge
    split_ifs <;> omega

/-- Both ends detached (all four neighbours outside the first edge): the first edge becomes an
isolated `single`, vertex orbits gain four, face orbits gain four or stay (the latter exactly
when its two sides are not cofacial). -/
theorem detach_both {es : List (Nat × Nat)} {p : Nat × Nat} {n : Nat} {rs : RotationSystem}
    (h : IsEmbedding (p :: es) n rs)
    (h0 : 4 ≤ rs.rot 0) (h1 : 4 ≤ rs.rot 1) (h2 : 4 ≤ rs.rot 2) (h3 : 4 ≤ rs.rot 3) :
    ∃ H : RotationSystem, H.Total ∧ H.Involution ∧ H.OppositeDir ∧ H.SameVertex (p :: es) ∧
      H.size = rs.size ∧ H.rot 0 = 1 ∧ H.rot 1 = 0 ∧ H.rot 2 = 3 ∧ H.rot 3 = 2 ∧
      H.numVertexOrbits = rs.numVertexOrbits + 4 ∧
      (H.numFaceOrbits = rs.numFaceOrbits + 4 ∨
        (H.numFaceOrbits = rs.numFaceOrbits ∧ ¬SameOrbit (rs.stepC 3) 2 0)) := by
  have hm : rs.size = 4 * (es.length + 1) := by rw [h.size]; rfl
  have ht := h.total
  have hi := h.involution
  have ho := h.opposite_dir
  have hsv := h.same_vertex
  have hv0 : QE.vert (p :: es) 0 = some p.1 := by rw [vert_cons_lt _ _ (by decide)]; rfl
  have hv1 : QE.vert (p :: es) 1 = some p.1 := by rw [vert_cons_lt _ _ (by decide)]; rfl
  have hv2 : QE.vert (p :: es) 2 = some p.2 := by rw [vert_cons_lt _ _ (by decide)]; rfl
  have hv3 : QE.vert (p :: es) 3 = some p.2 := by rw [vert_cons_lt _ _ (by decide)]; rfl
  have d1 : DetachAt rs 1 0 := ⟨ht, hi, ho, by omega, by omega, by omega, by decide, h0, h1⟩
  set H1 := rs.conj 1 (rs.rot 0) with hH1def
  have hs1 : H1.size = rs.size := d1.conj_size'
  have hm1 : H1.size = 4 * (es.length + 1) := hs1.trans hm
  have hH1t : H1.Total := d1.conj_total'
  have hH1i : H1.Involution := d1.conj_involution'
  have hH1o : H1.OppositeDir := d1.conj_oppositeDir'
  have hH1sv : H1.SameVertex (p :: es) := d1.conj_sameVertex' hsv (hv1.trans hv0.symm) hv0
  have hH1 : ∀ q, q < 4 → H1.rot q = if q = 0 then 1 else if q = 1 then 0 else rs.rot q := by
    intro q hq
    rw [d1.rot_eq q (by omega)]
    split_ifs <;> omega
  have hH1_0 : H1.rot 0 = 1 := by rw [hH1 0 (by omega)]; simp
  have hH1_1 : H1.rot 1 = 0 := by rw [hH1 1 (by omega)]; simp
  have hH1_2 : H1.rot 2 = rs.rot 2 := by rw [hH1 2 (by omega)]; simp
  have hH1_3 : H1.rot 3 = rs.rot 3 := by rw [hH1 3 (by omega)]; simp
  have d2 : DetachAt H1 3 2 := ⟨hH1t, hH1i, hH1o, by omega, by omega, by omega, by decide,
    by rw [hH1_2]; exact h2, by rw [hH1_3]; exact h3⟩
  set H := H1.conj 3 (H1.rot 2) with hHdef
  have hs2 : H.size = H1.size := d2.conj_size'
  have hH : ∀ q, q < 4 → H.rot q = if q = 2 then 3 else if q = 3 then 2 else H1.rot q := by
    intro q hq
    rw [d2.rot_eq q (by omega)]
    have := d2.y_ge
    have := d2.z_ge
    split_ifs <;> omega
  obtain ⟨a1, a2⟩ := d1.vertexOrbits_of_embedding h (hv1.trans hv0.symm)
  have hV1 : rs.numVertexOrbits + 2 = H1.numVertexOrbits := d1.numVertexOrbits_eq hm a1 a2
  have d2rs : DetachAt rs 3 2 := ⟨ht, hi, ho, by omega, by omega, by omega, by decide, h2, h3⟩
  obtain ⟨b1, b2⟩ := d2rs.vertexOrbits_of_embedding h (hv3.trans hv2.symm)
  have e31 : (3:Nat) ^^^ 1 = 2 := by decide
  have e21 : (2:Nat) ^^^ 1 = 3 := by decide
  have hx2 : 4 * 1 ≤ rs.rot 2 ^^^ 1 := ge_of_xor_ge (by omega) (by decide)
  have hx3 : 4 * 1 ≤ rs.rot 3 ^^^ 1 := ge_of_xor_ge (by omega) (by decide)
  have b1' : SameOrbit (H1.stepC 1) (3 ^^^ 1) (rs.rot 2 ^^^ 1) :=
    d1.sameOrbit_stepC1_conj hm (by decide) (by rw [e31]; decide) (by rw [e31]; decide)
      (by omega) (by omega) b1
  have b2' : SameOrbit (H1.stepC 1) (rs.rot 3 ^^^ 1) (2 ^^^ 1) :=
    d1.sameOrbit_stepC1_conj hm (by decide) (by omega) (by omega) (by rw [e21]; decide)
      (by rw [e21]; decide) b2
  rw [← hH1_2] at b1'
  rw [← hH1_3] at b2'
  have hV2 : H1.numVertexOrbits + 2 = H.numVertexOrbits := d2.numVertexOrbits_eq hm1 b1' b2'
  have e33 : (3:Nat) ^^^ 3 = 0 := by decide
  have e13 : (1:Nat) ^^^ 3 = 2 := by decide
  have e23 : (2:Nat) ^^^ 3 = 1 := by decide
  have hside2 : SameOrbit (H1.stepC 3) (3 ^^^ 3) 2 := by
    rw [e33]
    have : H1.stepC 3 2 = 0 := by
      rw [stepC_eq_rot hH1t hm1 (by decide) (by omega), e23, hH1_1]
    exact (sameOrbit_of_eq this).symm
  have hF2 : H1.numFaceOrbits + 2 = H.numFaceOrbits := d2.numFaceOrbits_split hm1 hside2
  refine ⟨H, d2.conj_total', d2.conj_involution', d2.conj_oppositeDir',
    d2.conj_sameVertex' hH1sv (hv3.trans hv2.symm) hv2, hs2.trans hs1, ?_, ?_, ?_, ?_,
    by omega, ?_⟩
  · rw [hH 0 (by omega)]; simp [hH1_0]
  · rw [hH 1 (by omega)]; simp [hH1_1]
  · rw [hH 2 (by omega)]; simp
  · rw [hH 3 (by omega)]; simp
  · by_cases hside1 : SameOrbit (rs.stepC 3) (1 ^^^ 3) 0
    · have hF1 : rs.numFaceOrbits + 2 = H1.numFaceOrbits := d1.numFaceOrbits_split hm hside1
      left; omega
    · have hF1 : H1.numFaceOrbits + 2 = rs.numFaceOrbits := d1.numFaceOrbits_merge hm hside1
      right
      rw [e13] at hside1
      exact ⟨by omega, hside1⟩

/-- **Genus inequality**: every rotation system of `es` satisfies `F + 2V ≤ 2(2C + E)`. -/
theorem genus_le : ∀ (es : List (Nat × Nat)) {n : Nat} {rs : RotationSystem},
    IsEmbedding es n rs →
    rs.numFaceOrbits + 2 * numNonIsolated es n ≤ 2 * (2 * numComponents es n + es.length)
  | [], n, rs, h => by
    have hs : rs.size = 0 := by rw [h.size]; rfl
    have hF : rs.numFaceOrbits = 0 := by
      unfold numFaceOrbits numOrbits; rw [hs]; rfl
    have hNI : numNonIsolated [] n = 0 := by
      rw [numNonIsolated_eq_card]; simp [HasEdge]
    omega
  | p :: es, n, rs, h => by
    have IH : ∀ ρ, IsEmbedding es n ρ →
        ρ.numFaceOrbits + 2 * numNonIsolated es n ≤ 2 * (2 * numComponents es n + es.length) :=
      fun ρ hρ => genus_le es hρ
    have hm : rs.size = 4 * (es.length + 1) := by rw [h.size]; rfl
    have ht := h.total
    have hi := h.involution
    have ho := h.opposite_dir
    have hsv := h.same_vertex
    have hverts : ∀ q ∈ es, q.1 < n ∧ q.2 < n := fun q hq =>
      h.verts q (List.mem_cons_of_mem _ hq)
    have hp : p.1 < n ∧ p.2 < n := h.verts p (List.mem_cons_self ..)
    have hV0 := h.vertex_orbits
    have hNI := numNonIsolated_cons_eq (es := es) hp
    have hv0 : QE.vert (p :: es) 0 = some p.1 := by rw [vert_cons_lt _ _ (by decide)]; rfl
    have hv1 : QE.vert (p :: es) 1 = some p.1 := by rw [vert_cons_lt _ _ (by decide)]; rfl
    have hv2 : QE.vert (p :: es) 2 = some p.2 := by rw [vert_cons_lt _ _ (by decide)]; rfl
    have hv3 : QE.vert (p :: es) 3 = some p.2 := by rw [vert_cons_lt _ _ (by decide)]; rfl
    have hrot : ∀ q, q < rs.size → QE.vert (p :: es) (rs.rot q) = QE.vert (p :: es) q :=
      fun q hq => vert_rot ht hsv hq
    have hr0 := rot_lt ht hi (q := 0) (by omega)
    have hr1 := rot_lt ht hi (q := 1) (by omega)
    have hr2 := rot_lt ht hi (q := 2) (by omega)
    have hr3 := rot_lt ht hi (q := 3) (by omega)
    have hp0 := rot_mod2 ht ho (q := 0) (by omega)
    have hp1 := rot_mod2 ht ho (q := 1) (by omega)
    have hp2 := rot_mod2 ht ho (q := 2) (by omega)
    have hp3 := rot_mod2 ht ho (q := 3) (by omega)
    have hrr0 := rot_rot ht hi (q := 0) (by omega)
    have hrr1 := rot_rot ht hi (q := 1) (by omega)
    have hrr2 := rot_rot ht hi (q := 2) (by omega)
    have hrr3 := rot_rot ht hi (q := 3) (by omega)
    have hstep1 : ∀ q, q < rs.size → rs.stepC 1 q = rs.rot (q ^^^ 1) :=
      fun q hq => stepC_eq_rot ht hm (by decide) hq
    have hstep3 : ∀ q, q < rs.size → rs.stepC 3 q = rs.rot (q ^^^ 3) :=
      fun q hq => stepC_eq_rot ht hm (by decide) hq
    have e01 : (0:Nat) ^^^ 1 = 1 := by decide
    have e21 : (2:Nat) ^^^ 1 = 3 := by decide
    have e03 : (0:Nat) ^^^ 3 = 3 := by decide
    have e13 : (1:Nat) ^^^ 3 = 2 := by decide
    have e23 : (2:Nat) ^^^ 3 = 1 := by decide
    have e33 : (3:Nat) ^^^ 3 = 0 := by decide
    have hvv : ∀ q r, q < 4 → r < 4 → rs.rot q = r →
        QE.vert (p :: es) r = QE.vert (p :: es) q := by
      intro q r hq hr hqr
      have := hrot q (by omega)
      rwa [hqr] at this
    have hbig : ∀ q, q < 4 → 4 ≤ rs.rot q → ∀ u, QE.vert (p :: es) q = some u →
        HasEdge es u := by
      intro q hq hge u hu
      exact hasEdge_of_vert_ge hge (by rw [hrot q (by omega)]; exact hu)
    have hpair : ∀ a b u, a < 4 → b < 4 → rs.rot (a ^^^ 1) = b → rs.rot (b ^^^ 1) = a →
        QE.vert (p :: es) a = some u → ¬HasEdge es u := by
      intro a b u ha hb hab hba hu
      exact not_hasEdge_of_pair h ha hb (by rw [hstep1 a (by omega)]; exact hab)
        (by rw [hstep1 b (by omega)]; exact hba) hu
    obtain ⟨s0, s1, s2, s3⟩ := single_rot4
    obtain ⟨l0, l1, l2, l3⟩ := loop_rot4
    have hE : (p :: es).length = es.length + 1 := rfl
    rw [hE]
    rcases Decidable.em (p.1 = p.2) with huv | huv
    · -- loop
      have hl20 : lbl (p :: es) 2 = lbl (p :: es) 0 := by simp [lbl, hv2, hv0, huv]
      have h01 : rs.rot 0 ≠ 1 := fun e => by
        have h10 : rs.rot 1 = 0 := by have := hrr0; rw [e] at this; exact this
        have := mem_pair_of_lbl h (a := 0) (b := 0) (w := 2) (by omega) (by omega)
          (by rw [hstep1 0 (by omega), e01]; exact h10)
          (by rw [hstep1 0 (by omega), e01]; exact h10) hl20
        omega
      have h32 : rs.rot 3 ≠ 2 := fun e => by
        have := mem_pair_of_lbl h (a := 2) (b := 2) (w := 0) (by omega) (by omega)
          (by rw [hstep1 2 (by omega), e21]; exact e)
          (by rw [hstep1 2 (by omega), e21]; exact e) hl20.symm
        omega
      have h10 : rs.rot 1 ≠ 0 := fun e => h01 (by have := hrr1; rw [e] at this; exact this)
      have h23 : rs.rot 2 ≠ 3 := fun e => h32 (by have := hrr2; rw [e] at this; exact this)
      rcases (show rs.rot 0 = 3 ∨ 4 ≤ rs.rot 0 by omega) with h03 | h0ge <;>
      rcases (show rs.rot 1 = 2 ∨ 4 ≤ rs.rot 1 by omega) with h12 | h1ge
      · have h30 : rs.rot 3 = 0 := by have := hrr0; rw [h03] at this; exact this
        have h21 : rs.rot 2 = 1 := by have := hrr1; rw [h12] at this; exact this
        have hu : ¬HasEdge es p.1 := hpair 0 2 p.1 (by omega) (by omega)
          (by rw [e01]; exact h12) (by rw [e21]; exact h30) hv0
        rw [ite_eq_right hu, ite_eq_left (Or.inl huv.symm)] at hNI
        have hC := numComponents_cons_succ_le hverts hp hu (by rw [← huv]; exact hu)
        have hb := genus_cons_bound IH hverts loop_size loop_total loop_involution ht hi ho hsv
          h.size (by intro q hq; interval_cases q <;> omega)
          (by rw [loop_numVertexOrbits]; omega)
        rw [loop_numFaceOrbits] at hb
        omega
      · have h30 : rs.rot 3 = 0 := by have := hrr0; rw [h03] at this; exact this
        have h2ge : 4 ≤ rs.rot 2 := by
          by_contra hlt
          have h21 : rs.rot 2 = 1 := by omega
          have := hrr2; rw [h21] at this; omega
        have d : DetachAt rs 1 2 :=
          ⟨ht, hi, ho, by omega, by omega, by omega, by decide, h2ge, h1ge⟩
        obtain ⟨H, hHt, hHi, hHo, hHsv, hHs, hH2, hH1, hHq, hHV, hHF⟩ :=
          detach_one h d (by rw [hv1, hv2, huv]) hv2 (by rw [e13]; exact SameOrbit.refl _ _)
        have hH0 : H.rot 0 = 3 := by rw [hHq 0 (by omega) (by omega) (by omega)]; exact h03
        have hH3 : H.rot 3 = 0 := by rw [hHq 3 (by omega) (by omega) (by omega)]; exact h30
        have hu : HasEdge es p.1 := hbig 1 (by omega) h1ge p.1 hv1
        rw [ite_eq_left hu, ite_eq_left (Or.inl huv.symm)] at hNI
        have hC := numComponents_cons_le_of_loop hverts hp huv
        have hb := genus_cons_bound IH hverts loop_size loop_total loop_involution hHt hHi hHo
          hHsv (hHs.trans h.size) (by intro q hq; interval_cases q <;> omega)
          (by rw [loop_numVertexOrbits]; omega)
        rw [loop_numFaceOrbits] at hb
        omega
      · have h21 : rs.rot 2 = 1 := by have := hrr1; rw [h12] at this; exact this
        have h3ge : 4 ≤ rs.rot 3 := by
          by_contra hlt
          have h30 : rs.rot 3 = 0 := by omega
          have := hrr3; rw [h30] at this; omega
        have d : DetachAt rs 3 0 :=
          ⟨ht, hi, ho, by omega, by omega, by omega, by decide, h0ge, h3ge⟩
        obtain ⟨H, hHt, hHi, hHo, hHsv, hHs, hH0, hH3, hHq, hHV, hHF⟩ :=
          detach_one h d (by rw [hv3, hv0, huv]) hv0 (by rw [e33]; exact SameOrbit.refl _ _)
        have hH1 : H.rot 1 = 2 := by rw [hHq 1 (by omega) (by omega) (by omega)]; exact h12
        have hH2 : H.rot 2 = 1 := by rw [hHq 2 (by omega) (by omega) (by omega)]; exact h21
        have hu : HasEdge es p.1 := hbig 0 (by omega) h0ge p.1 hv0
        rw [ite_eq_left hu, ite_eq_left (Or.inl huv.symm)] at hNI
        have hC := numComponents_cons_le_of_loop hverts hp huv
        have hb := genus_cons_bound IH hverts loop_size loop_total loop_involution hHt hHi hHo
          hHsv (hHs.trans h.size) (by intro q hq; interval_cases q <;> omega)
          (by rw [loop_numVertexOrbits]; omega)
        rw [loop_numFaceOrbits] at hb
        omega
      · have h2ge : 4 ≤ rs.rot 2 := by
          by_contra hlt
          have h21 : rs.rot 2 = 1 := by omega
          have := hrr2; rw [h21] at this; omega
        have h3ge : 4 ≤ rs.rot 3 := by
          by_contra hlt
          have h30 : rs.rot 3 = 0 := by omega
          have := hrr3; rw [h30] at this; omega
        obtain ⟨H, hHt, hHi, hHo, hHsv, hHs, hH0, hH1, hH2, hH3, hHV, hHF⟩ :=
          detach_both h h0ge h1ge h2ge h3ge
        have hu : HasEdge es p.1 := hbig 0 (by omega) h0ge p.1 hv0
        rw [ite_eq_left hu, ite_eq_left (Or.inl huv.symm)] at hNI
        have hC := numComponents_cons_le_of_loop hverts hp huv
        have hb := genus_cons_bound IH hverts single_size single_total single_involution hHt hHi
          hHo hHsv (hHs.trans h.size) (by intro q hq; interval_cases q <;> omega)
          (by rw [single_numVertexOrbits]; omega)
        rw [single_numFaceOrbits] at hb
        rcases hHF with hHF | ⟨hHF, -⟩ <;> omega
    · -- two distinct endpoints
      have huv' : ¬p.2 = p.1 := fun e => huv e.symm
      have h03 : rs.rot 0 ≠ 3 := fun e => huv (by
        have := hvv 0 3 (by omega) (by omega) e
        rw [hv3, hv0] at this
        exact (Option.some.inj this).symm)
      have h12 : rs.rot 1 ≠ 2 := fun e => huv (by
        have := hvv 1 2 (by omega) (by omega) e
        rw [hv2, hv1] at this
        exact (Option.some.inj this).symm)
      have h21 : rs.rot 2 ≠ 1 := fun e => huv (by
        have := hvv 2 1 (by omega) (by omega) e
        rw [hv1, hv2] at this
        exact Option.some.inj this)
      have h30 : rs.rot 3 ≠ 0 := fun e => huv (by
        have := hvv 3 0 (by omega) (by omega) e
        rw [hv0, hv3] at this
        exact Option.some.inj this)
      rcases (show rs.rot 0 = 1 ∨ 4 ≤ rs.rot 0 by omega) with h01 | h0ge <;>
      rcases (show rs.rot 2 = 3 ∨ 4 ≤ rs.rot 2 by omega) with h23 | h2ge
      · have h10 : rs.rot 1 = 0 := by have := hrr0; rw [h01] at this; exact this
        have h32 : rs.rot 3 = 2 := by have := hrr2; rw [h23] at this; exact this
        have hu : ¬HasEdge es p.1 := hpair 0 0 p.1 (by omega) (by omega)
          (by rw [e01]; exact h10) (by rw [e01]; exact h10) hv0
        have hv : ¬HasEdge es p.2 := hpair 2 2 p.2 (by omega) (by omega)
          (by rw [e21]; exact h32) (by rw [e21]; exact h32) hv2
        rw [ite_eq_right hu, ite_eq_right (not_or.2 ⟨huv', hv⟩)] at hNI
        have hC := numComponents_cons_succ_le hverts hp hu hv
        have hb := genus_cons_bound IH hverts single_size single_total single_involution ht hi ho
          hsv h.size (by intro q hq; interval_cases q <;> omega)
          (by rw [single_numVertexOrbits]; omega)
        rw [single_numFaceOrbits] at hb
        omega
      · have h10 : rs.rot 1 = 0 := by have := hrr0; rw [h01] at this; exact this
        have h3ge : 4 ≤ rs.rot 3 := by
          by_contra hlt
          have h32 : rs.rot 3 = 2 := by omega
          have := hrr3; rw [h32] at this; omega
        have d : DetachAt rs 3 2 :=
          ⟨ht, hi, ho, by omega, by omega, by omega, by decide, h2ge, h3ge⟩
        obtain ⟨H, hHt, hHi, hHo, hHsv, hHs, hH2, hH3, hHq, hHV, hHF⟩ :=
          detach_one h d (hv3.trans hv2.symm) hv2 (by
            rw [e33]
            exact (sameOrbit_of_eq
              (by rw [hstep3 2 (by omega), e23, h10] : rs.stepC 3 2 = 0)).symm)
        have hH0 : H.rot 0 = 1 := by rw [hHq 0 (by omega) (by omega) (by omega)]; exact h01
        have hH1 : H.rot 1 = 0 := by rw [hHq 1 (by omega) (by omega) (by omega)]; exact h10
        have hu : ¬HasEdge es p.1 := hpair 0 0 p.1 (by omega) (by omega)
          (by rw [e01]; exact h10) (by rw [e01]; exact h10) hv0
        have hv : HasEdge es p.2 := hbig 2 (by omega) h2ge p.2 hv2
        rw [ite_eq_right hu, ite_eq_left (Or.inr hv)] at hNI
        have hC := numComponents_cons_le_of_isolated hverts hp (Or.inl hu)
        have hb := genus_cons_bound IH hverts single_size single_total single_involution hHt hHi
          hHo hHsv (hHs.trans h.size) (by intro q hq; interval_cases q <;> omega)
          (by rw [single_numVertexOrbits]; omega)
        rw [single_numFaceOrbits] at hb
        omega
      · have h32 : rs.rot 3 = 2 := by have := hrr2; rw [h23] at this; exact this
        have h1ge : 4 ≤ rs.rot 1 := by
          by_contra hlt
          have h10 : rs.rot 1 = 0 := by omega
          have := hrr1; rw [h10] at this; omega
        have d : DetachAt rs 1 0 :=
          ⟨ht, hi, ho, by omega, by omega, by omega, by decide, h0ge, h1ge⟩
        obtain ⟨H, hHt, hHi, hHo, hHsv, hHs, hH0, hH1, hHq, hHV, hHF⟩ :=
          detach_one h d (hv1.trans hv0.symm) hv0 (by
            rw [e13]
            exact (sameOrbit_of_eq
              (by rw [hstep3 0 (by omega), e03, h32] : rs.stepC 3 0 = 2)).symm)
        have hH2 : H.rot 2 = 3 := by rw [hHq 2 (by omega) (by omega) (by omega)]; exact h23
        have hH3 : H.rot 3 = 2 := by rw [hHq 3 (by omega) (by omega) (by omega)]; exact h32
        have hu : HasEdge es p.1 := hbig 0 (by omega) h0ge p.1 hv0
        have hv : ¬HasEdge es p.2 := hpair 2 2 p.2 (by omega) (by omega)
          (by rw [e21]; exact h32) (by rw [e21]; exact h32) hv2
        rw [ite_eq_left hu, ite_eq_right (not_or.2 ⟨huv', hv⟩)] at hNI
        have hC := numComponents_cons_le_of_isolated hverts hp (Or.inr hv)
        have hb := genus_cons_bound IH hverts single_size single_total single_involution hHt hHi
          hHo hHsv (hHs.trans h.size) (by intro q hq; interval_cases q <;> omega)
          (by rw [single_numVertexOrbits]; omega)
        rw [single_numFaceOrbits] at hb
        omega
      · have h1ge : 4 ≤ rs.rot 1 := by
          by_contra hlt
          have h10 : rs.rot 1 = 0 := by omega
          have := hrr1; rw [h10] at this; omega
        have h3ge : 4 ≤ rs.rot 3 := by
          by_contra hlt
          have h32 : rs.rot 3 = 2 := by omega
          have := hrr3; rw [h32] at this; omega
        obtain ⟨H, hHt, hHi, hHo, hHsv, hHs, hH0, hH1, hH2, hH3, hHV, hHF⟩ :=
          detach_both h h0ge h1ge h2ge h3ge
        have hu : HasEdge es p.1 := hbig 0 (by omega) h0ge p.1 hv0
        have hv : HasEdge es p.2 := hbig 2 (by omega) h2ge p.2 hv2
        rw [ite_eq_left hu, ite_eq_left (Or.inr hv)] at hNI
        have hb := genus_cons_bound IH hverts single_size single_total single_involution hHt hHi
          hHo hHsv (hHs.trans h.size) (by intro q hq; interval_cases q <;> omega)
          (by rw [single_numVertexOrbits]; omega)
        rw [single_numFaceOrbits] at hb
        rcases hHF with hHF | ⟨hHF, hside⟩
        · have hC := numComponents_cons_le_succ hverts hp
          omega
        · have hC := numComponents_cons hverts hp huv (edgesConn_of_not_cofacial h h0ge h3ge hside)
          omega

end RotationSystem
end Spqr
