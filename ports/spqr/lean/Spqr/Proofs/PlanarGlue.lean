import Spqr.Proofs.Orbits

/-!
# 2-sum gluing (PROOF.md §8.5)

The splice `TwoSum.splice` of two planar embeddings along a virtual edge is a planar embedding of
the 2-sum: the face orbits are `faces₁ + faces₂ − 4` (two faces per direction merge) and the
vertex orbits `verts₁ + verts₂ − 4`.
-/

namespace Spqr

open Classical

/-! ### Quarter-edge arithmetic -/

theorem xor_div4 (q : Nat) {c : Nat} (hc : c < 4) : (q ^^^ c) / 4 = q / 4 := by
  have := @Nat.xor_div_two_pow q c 2
  have h0 : c / 4 = 0 := Nat.div_eq_of_lt hc
  simpa [h0] using this

theorem xor_mod4 (q : Nat) {c : Nat} (hc : c < 4) : (q ^^^ c) % 4 = q % 4 ^^^ c := by
  have := @Nat.xor_mod_two_pow q c 2
  have h0 : c % 4 = c := Nat.mod_eq_of_lt hc
  simpa [h0] using this

theorem xor_lt4 {j c : Nat} (hj : j < 4) (hc : c < 4) : j ^^^ c < 4 :=
  @Nat.xor_lt_two_pow j c 2 hj hc

theorem xor_eq4 (q : Nat) {c : Nat} (hc : c < 4) : q ^^^ c = 4 * (q / 4) + (q % 4 ^^^ c) := by
  conv_lhs => rw [← Nat.div_add_mod (q ^^^ c) 4]
  rw [xor_div4 q hc, xor_mod4 q hc]

theorem xor_xor_self (q c : Nat) : (q ^^^ c) ^^^ c = q := by
  rw [Nat.xor_assoc, Nat.xor_self, Nat.xor_zero]

theorem xor_mod2 (q : Nat) {c : Nat} (hc : c % 2 = 1) : (q ^^^ c) % 2 ≠ q % 2 := by
  have := @Nat.xor_mod_two_pow q c 1
  simp only [Nat.pow_one, hc] at this
  rw [this]
  rcases Nat.mod_two_eq_zero_or_one q with h | h <;> simp [h]

theorem xor_lt_mul4 {q m c : Nat} (hq : q < 4 * m) (hc : c < 4) : q ^^^ c < 4 * m := by
  rw [xor_eq4 q hc]
  have := xor_lt4 (Nat.mod_lt q (by omega)) hc
  have : q / 4 < m := by omega
  omega

theorem mul4_add_xor (m j : Nat) {c : Nat} (hj : j < 4) (hc : c < 4) :
    (4 * m + j) ^^^ c = 4 * m + (j ^^^ c) := by
  rw [xor_eq4 _ hc]
  have h1 : (4 * m + j) / 4 = m := by omega
  have h2 : (4 * m + j) % 4 = j := by omega
  rw [h1, h2]

theorem ge_of_xor_ge {q m c : Nat} (hq : 4 * m ≤ q) (hc : c < 4) : 4 * m ≤ q ^^^ c := by
  rw [xor_eq4 q hc]; omega

/-! ### Index bookkeeping -/

theorem insIdx_delIdx {e f : Nat} (h : f ≠ e) : insIdx e (delIdx e f) = f := by
  simp only [insIdx, delIdx]; split_ifs <;> omega

theorem delIdx_insIdx (e g : Nat) : delIdx e (insIdx e g) = g := by
  simp only [insIdx, delIdx]; split_ifs <;> omega

theorem insIdx_ne (e g : Nat) : insIdx e g ≠ e := by
  simp only [insIdx]; split_ifs <;> omega

theorem delIdx_lt {e f m : Nat} (hf : f < m) (h : f ≠ e) (he : e < m) : delIdx e f < m - 1 := by
  simp only [delIdx]; split_ifs <;> omega

theorem insIdx_lt {e g m : Nat} (hg : g < m - 1) (he : e < m) : insIdx e g < m := by
  simp only [insIdx]; split_ifs <;> omega

theorem delIdx_inj {e f f' : Nat} (hf : f ≠ e) (hf' : f' ≠ e) (h : delIdx e f = delIdx e f') :
    f = f' := by
  simp only [delIdx] at h; split_ifs at h <;> omega

/-! ### Rotation systems -/

namespace RotationSystem

variable {rs : RotationSystem}

/-- Totalised rotation. -/
def rot (rs : RotationSystem) (q : Nat) : Nat := (rs.get q).getD q

theorem get_eq_none {q : Nat} (hq : rs.size ≤ q) : rs.get q = none := by
  unfold get
  rw [Array.getElem?_eq_none hq]
  rfl

theorem rot_of_ge {q : Nat} (hq : rs.size ≤ q) : rs.rot q = q := by
  unfold rot; rw [get_eq_none hq]; rfl

theorem get_eq_rot (ht : rs.Total) {q : Nat} (hq : q < rs.size) : rs.get q = some (rs.rot q) := by
  obtain ⟨r, hr⟩ := Option.isSome_iff_exists.1 (ht q hq)
  unfold rot; rw [hr]; rfl

theorem rot_lt (ht : rs.Total) (hi : rs.Involution) {q : Nat} (hq : q < rs.size) :
    rs.rot q < rs.size :=
  (hi q hq _ (by rw [get_eq_rot ht hq]; rfl)).1

theorem rot_rot (ht : rs.Total) (hi : rs.Involution) {q : Nat} (hq : q < rs.size) :
    rs.rot (rs.rot q) = q := by
  have h := (hi q hq _ (by rw [get_eq_rot ht hq]; rfl)).2
  show (rs.get (rs.rot q)).getD _ = q
  rw [h]; rfl

theorem rot_mod2 (ht : rs.Total) (ho : rs.OppositeDir) {q : Nat} (hq : q < rs.size) :
    rs.rot q % 2 ≠ q % 2 :=
  ho q hq _ (by rw [get_eq_rot ht hq]; rfl)

theorem rot_ne (ht : rs.Total) (ho : rs.OppositeDir) {q : Nat} (hq : q < rs.size) :
    rs.rot q ≠ q := fun h => rot_mod2 ht ho hq (by rw [h])

theorem vert_rot {es : List (Nat × Nat)} (ht : rs.Total) (hv : rs.SameVertex es) {q : Nat}
    (hq : q < rs.size) : QE.vert es (rs.rot q) = QE.vert es q :=
  (hv q hq _ (by rw [get_eq_rot ht hq]; rfl)).symm

theorem rot_inj (ht : rs.Total) (hi : rs.Involution) {q r : Nat} (hq : q < rs.size)
    (hr : r < rs.size) (h : rs.rot q = rs.rot r) : q = r := by
  rw [← rot_rot ht hi hq, h, rot_rot ht hi hr]

theorem stepC_apply (c q : Nat) : rs.stepC c q = (rs.get (q ^^^ c)).getD q := rfl

theorem stepC_eq_rot (ht : rs.Total) {m c q : Nat} (hs : rs.size = 4 * m) (hc : c < 4)
    (hq : q < rs.size) : rs.stepC c q = rs.rot (q ^^^ c) := by
  rw [stepC_apply, get_eq_rot ht (hs ▸ xor_lt_mul4 (hs ▸ hq) hc)]; rfl

theorem stepC_of_ge {m c q : Nat} (hs : rs.size = 4 * m) (hc : c < 4) (hq : rs.size ≤ q) :
    rs.stepC c q = q := by
  rw [stepC_apply, get_eq_none (hs ▸ ge_of_xor_ge (hs ▸ hq) hc)]; rfl

theorem isPermOn_stepC (ht : rs.Total) (hi : rs.Involution) {m c : Nat} (hs : rs.size = 4 * m)
    (hc : c < 4) : IsPermOn (rs.stepC c) (Finset.range rs.size) where
  maps q hq := by
    rw [Finset.mem_range] at hq ⊢
    rw [stepC_eq_rot ht hs hc hq]
    exact rot_lt ht hi (hs ▸ xor_lt_mul4 (hs ▸ hq) hc)
  inj q hq r hr h := by
    rw [Finset.mem_range] at hq hr
    rw [stepC_eq_rot ht hs hc hq, stepC_eq_rot ht hs hc hr] at h
    have := rot_inj ht hi (hs ▸ xor_lt_mul4 (hs ▸ hq) hc) (hs ▸ xor_lt_mul4 (hs ▸ hr) hc) h
    rw [← xor_xor_self q c, this, xor_xor_self]
  fix q hq := stepC_of_ge hs hc (by simpa using hq)

theorem stepC_mod2 (ht : rs.Total) (ho : rs.OppositeDir) {m c : Nat} (hs : rs.size = 4 * m)
    (hc : c < 4) (hc2 : c % 2 = 1) (q : Nat) : rs.stepC c q % 2 = q % 2 := by
  by_cases hq : q < rs.size
  · rw [stepC_eq_rot ht hs hc hq]
    have h1 := rot_mod2 ht ho (hs ▸ xor_lt_mul4 (hs ▸ hq) hc)
    have h2 := xor_mod2 q hc2
    omega
  · rw [stepC_of_ge hs hc (not_lt.1 hq)]

theorem stepC_lt (ht : rs.Total) (hi : rs.Involution) {m c q : Nat} (hs : rs.size = 4 * m)
    (hc : c < 4) (hq : q < rs.size) : rs.stepC c q < rs.size :=
  by simpa using (isPermOn_stepC ht hi hs hc).maps q (by simpa using hq)

/-- The face/vertex orbit relation of `Planar.lean` agrees with `SameOrbit` of the total step. -/
theorem sameOrbit_stepC_iff (ht : rs.Total) (hi : rs.Involution) {m c q r : Nat}
    (hs : rs.size = 4 * m) (hc : c < 4) (hq : q < rs.size) :
    SameOrbit (rs.stepC c) q r ↔
      Relation.ReflTransGen (fun a b => rs.get (a ^^^ c) = some b) q r := by
  have hperm := isPermOn_stepC ht hi hs hc
  constructor
  · intro h
    have h' := hperm.reach_of_sameOrbit (by simpa using hq) h
    unfold Reach at h'
    clear h
    induction h' with
    | refl => exact Relation.ReflTransGen.refl
    | @tail b d hb hd ih =>
      refine ih.tail ?_
      have hbm : b < rs.size := by
        have := hperm.mem_iff_of_sameOrbit (IsPermOn.sameOrbit_of_reach hb)
        simpa using this.1 (by simpa using hq)
      rw [stepC_eq_rot ht hs hc hbm] at hd
      rw [← hd, get_eq_rot ht (hs ▸ xor_lt_mul4 (hs ▸ hbm) hc)]
  · intro h
    induction h with
    | refl => exact SameOrbit.refl _ _
    | @tail b d _ hd ih =>
      refine ih.trans ?_
      have : rs.stepC c b = d := by rw [stepC_apply, hd]; rfl
      exact this ▸ SameOrbit.step b

theorem sameFaceOrbit_iff (ht : rs.Total) (hi : rs.Involution) {m q r : Nat}
    (hs : rs.size = 4 * m) (hq : q < rs.size) :
    rs.SameFaceOrbit q r ↔ SameOrbit (rs.stepC 3) q r :=
  (sameOrbit_stepC_iff ht hi hs (by decide) hq).symm

/-- `q ↦ q ^^^ c` conjugates the step to its inverse, so it preserves orbits. -/
theorem sameOrbit_stepC_xor (ht : rs.Total) (hi : rs.Involution) {m c : Nat}
    (hs : rs.size = 4 * m) (hc : c < 4) {q r : Nat} (h : SameOrbit (rs.stepC c) q r) :
    SameOrbit (rs.stepC c) (q ^^^ c) (r ^^^ c) := by
  induction h with
  | rel a b h =>
    subst h
    by_cases ha : a < rs.size
    · have : rs.stepC c (rs.stepC c a ^^^ c) = a ^^^ c := by
        have hlt : rs.rot (a ^^^ c) ^^^ c < rs.size :=
          hs ▸ xor_lt_mul4 (hs ▸ rot_lt ht hi (hs ▸ xor_lt_mul4 (hs ▸ ha) hc)) hc
        rw [stepC_eq_rot ht hs hc ha, stepC_eq_rot ht hs hc hlt, xor_xor_self,
          rot_rot ht hi (hs ▸ xor_lt_mul4 (hs ▸ ha) hc)]
      exact (this ▸ SameOrbit.step (f := rs.stepC c) (rs.stepC c a ^^^ c)).symm
    · rw [stepC_of_ge hs hc (not_lt.1 ha)]; exact SameOrbit.refl _ _
  | refl => exact SameOrbit.refl _ _
  | symm _ _ _ ih => exact ih.symm
  | trans _ _ _ _ _ ih₁ ih₂ => exact ih₁.trans ih₂

theorem vert_stepC_one {es : List (Nat × Nat)} (ht : rs.Total) (hv : rs.SameVertex es) {m : Nat}
    (hs : rs.size = 4 * m) (q : Nat) : QE.vert es (rs.stepC 1 q) = QE.vert es q := by
  by_cases hq : q < rs.size
  · rw [stepC_eq_rot ht hs (by decide) hq, vert_rot ht hv (hs ▸ xor_lt_mul4 (hs ▸ hq) (by decide))]
    unfold QE.vert QE.edge QE.side
    rw [xor_div4 q (by decide)]
    congr 1
    funext p
    have : (q ^^^ 1) / 2 = q / 2 := by
      have := @Nat.xor_div_two_pow q 1 1
      simpa using this
    rw [this]
  · rw [stepC_of_ge hs (by decide) (not_lt.1 hq)]

end RotationSystem

namespace TwoSum

variable (T : TwoSum) {rs₁ rs₂ : RotationSystem}

/-- WIP (PROOF.md §8.5): the splice of two planar embeddings along the virtual edge is a planar
embedding of the 2-sum. Admitted; the embedding fields and the orbit counts are in progress. -/
theorem splice_isPlanarEmbedding (W : T.WF rs₁ rs₂) :
    IsPlanarEmbedding T.edges T.nVerts (T.splice rs₁ rs₂) := by
  sorry

theorem planar (W : T.WF rs₁ rs₂) : Planar T.edges T.nVerts :=
  ⟨_, T.splice_isPlanarEmbedding W⟩

/-- Converse (not attempted): a planar embedding of the 2-sum yields one of `G₁` by contracting
the `G₂` side onto `e₁`. -/
theorem planar_left (he₁ : T.es₁[T.e₁]? = some (T.u₁, T.v₁))
    (he₂ : T.es₂[T.e₂]? = some (T.u₂, T.v₂)) (h : Planar T.edges T.nVerts) :
    Planar T.es₁ T.n₁ := by
  sorry

end TwoSum

end Spqr
