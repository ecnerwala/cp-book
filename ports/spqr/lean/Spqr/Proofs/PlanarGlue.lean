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

theorem getElem?_eraseIdx_delIdx {α : Type} (l : List α) {e f : Nat} (h : f ≠ e) :
    (l.eraseIdx e)[delIdx e f]? = l[f]? := by
  rw [List.getElem?_eraseIdx]; unfold delIdx
  split_ifs <;> first | rfl | (congr 1; omega)

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

/-! ### Index maps -/

theorem qe₁_div (q : Nat) : T.qe₁ q / 4 = delIdx T.e₁ (q / 4) := by unfold qe₁; omega
theorem qe₁_mod (q : Nat) : T.qe₁ q % 4 = q % 4 := by unfold qe₁; omega
theorem qe₂_div (q : Nat) : (T.qe₂ q - T.off) / 4 = delIdx T.e₂ (q / 4) := by unfold qe₂; omega
theorem qe₂_mod (q : Nat) : (T.qe₂ q - T.off) % 4 = q % 4 := by unfold qe₂; omega
theorem off_le_qe₂ (q : Nat) : T.off ≤ T.qe₂ q := by unfold qe₂; omega

theorem pre₁_qe₁ {q : Nat} (h : q / 4 ≠ T.e₁) : T.pre₁ (T.qe₁ q) = q := by
  unfold pre₁; rw [qe₁_div, qe₁_mod, insIdx_delIdx h]; omega
theorem qe₁_pre₁ (r : Nat) : T.qe₁ (T.pre₁ r) = r := by
  unfold qe₁ pre₁
  have h1 : (4 * insIdx T.e₁ (r / 4) + r % 4) / 4 = insIdx T.e₁ (r / 4) := by omega
  have h2 : (4 * insIdx T.e₁ (r / 4) + r % 4) % 4 = r % 4 := by omega
  rw [h1, h2, delIdx_insIdx]; omega
theorem pre₂_qe₂ {q : Nat} (h : q / 4 ≠ T.e₂) : T.pre₂ (T.qe₂ q) = q := by
  unfold pre₂; rw [qe₂_div, qe₂_mod, insIdx_delIdx h]; omega
theorem qe₂_pre₂ {r : Nat} (hr : T.off ≤ r) : T.qe₂ (T.pre₂ r) = r := by
  unfold qe₂ pre₂
  have h1 : (4 * insIdx T.e₂ ((r - T.off) / 4) + (r - T.off) % 4) / 4 =
    insIdx T.e₂ ((r - T.off) / 4) := by omega
  have h2 : (4 * insIdx T.e₂ ((r - T.off) / 4) + (r - T.off) % 4) % 4 = (r - T.off) % 4 := by omega
  rw [h1, h2, delIdx_insIdx]; omega

theorem pre₁_div_ne (r : Nat) : T.pre₁ r / 4 ≠ T.e₁ := by
  unfold pre₁
  have : (4 * insIdx T.e₁ (r / 4) + r % 4) / 4 = insIdx T.e₁ (r / 4) := by omega
  rw [this]; exact insIdx_ne _ _
theorem pre₂_div_ne (r : Nat) : T.pre₂ r / 4 ≠ T.e₂ := by
  unfold pre₂
  have : (4 * insIdx T.e₂ ((r - T.off) / 4) + (r - T.off) % 4) / 4 =
    insIdx T.e₂ ((r - T.off) / 4) := by omega
  rw [this]; exact insIdx_ne _ _

theorem qe₁_xor {q c : Nat} (hc : c < 4) : T.qe₁ (q ^^^ c) = T.qe₁ q ^^^ c := by
  unfold qe₁
  rw [xor_div4 q hc, xor_mod4 q hc, mul4_add_xor _ _ (Nat.mod_lt _ (by omega)) hc]
theorem qe₂_xor {q c : Nat} (hc : c < 4) : T.qe₂ (q ^^^ c) = T.qe₂ q ^^^ c := by
  unfold qe₂ off
  rw [xor_div4 q hc, xor_mod4 q hc]
  have : 4 * (T.m₁ - 1) + 4 * delIdx T.e₂ (q / 4) + q % 4 =
      4 * (T.m₁ - 1 + delIdx T.e₂ (q / 4)) + q % 4 := by omega
  rw [this, mul4_add_xor _ _ (Nat.mod_lt _ (by omega)) hc]
  omega

theorem qe₁_mod2 (q : Nat) : T.qe₁ q % 2 = q % 2 := by unfold qe₁; omega
theorem qe₂_mod2 (q : Nat) : T.qe₂ q % 2 = q % 2 := by unfold qe₂ off; omega

/-- Neighbour of the `k`-th quarter of the virtual edge in `G₁`. -/
def ρ₁ (rs₁ : RotationSystem) (k : Nat) : Nat := rs₁.rot (4 * T.e₁ + k)
/-- Neighbour of the `k`-th quarter of the virtual edge in `G₂`. -/
def ρ₂ (rs₂ : RotationSystem) (k : Nat) : Nat := rs₂.rot (4 * T.e₂ + k)

theorem ρ₁_def (rs₁ : RotationSystem) (k : Nat) : T.ρ₁ rs₁ k = rs₁.rot (4 * T.e₁ + k) := rfl
theorem ρ₂_def (rs₂ : RotationSystem) (k : Nat) : T.ρ₂ rs₂ k = rs₂.rot (4 * T.e₂ + k) := rfl

/-- Total version of `spliceFn`. -/
def glue (rs₁ rs₂ : RotationSystem) (r : Nat) : Nat :=
  if r < T.off then
    if rs₁.rot (T.pre₁ r) / 4 = T.e₁ then T.qe₂ (T.ρ₂ rs₂ (rs₁.rot (T.pre₁ r) % 4 ^^^ 1))
    else T.qe₁ (rs₁.rot (T.pre₁ r))
  else
    if rs₂.rot (T.pre₂ r) / 4 = T.e₂ then T.qe₁ (T.ρ₁ rs₁ (rs₂.rot (T.pre₂ r) % 4 ^^^ 1))
    else T.qe₂ (rs₂.rot (T.pre₂ r))

/-! ### Consequences of well-formedness -/

section WF
variable (W : T.WF rs₁ rs₂)
include W

theorem e₁_lt : T.e₁ < T.m₁ := List.getElem?_eq_some_iff.1 W.e₁ |>.1
theorem e₂_lt : T.e₂ < T.m₂ := List.getElem?_eq_some_iff.1 W.e₂ |>.1
theorem size₁ : rs₁.size = 4 * T.m₁ := W.emb₁.size
theorem size₂ : rs₂.size = 4 * T.m₂ := W.emb₂.size

theorem qe₁_lt_off {q : Nat} (hq : q < 4 * T.m₁) (h : q / 4 ≠ T.e₁) : T.qe₁ q < T.off := by
  have := delIdx_lt (e := T.e₁) (by omega : q / 4 < T.m₁) h (T.e₁_lt W)
  unfold qe₁ off; omega
theorem qe₂_lt {q : Nat} (hq : q < 4 * T.m₂) (h : q / 4 ≠ T.e₂) :
    T.qe₂ q < 4 * (T.m₁ + T.m₂ - 2) := by
  have := delIdx_lt (e := T.e₂) (by omega : q / 4 < T.m₂) h (T.e₂_lt W)
  have := T.e₁_lt W
  unfold qe₂ off; omega
theorem pre₁_lt {r : Nat} (hr : r < T.off) : T.pre₁ r < 4 * T.m₁ := by
  have := insIdx_lt (e := T.e₁) (g := r / 4) (m := T.m₁) (by unfold off at hr; omega) (T.e₁_lt W)
  unfold pre₁; omega
theorem pre₂_lt {r : Nat} (hr : r < 4 * (T.m₁ + T.m₂ - 2)) (hr' : T.off ≤ r) :
    T.pre₂ r < 4 * T.m₂ := by
  have := T.e₁_lt W
  have := insIdx_lt (e := T.e₂) (g := (r - T.off) / 4) (m := T.m₂)
    (by unfold off at hr' ⊢; omega)
    (T.e₂_lt W)
  unfold pre₂; omega

theorem virt₁_lt {k : Nat} (hk : k < 4) : 4 * T.e₁ + k < rs₁.size := by
  rw [T.size₁ W]; have := T.e₁_lt W; omega
theorem virt₂_lt {k : Nat} (hk : k < 4) : 4 * T.e₂ + k < rs₂.size := by
  rw [T.size₂ W]; have := T.e₂_lt W; omega

theorem ρ₁_lt {k : Nat} (hk : k < 4) : T.ρ₁ rs₁ k < 4 * T.m₁ :=
  T.size₁ W ▸ RotationSystem.rot_lt W.emb₁.total W.emb₁.involution (T.virt₁_lt W hk)
theorem ρ₂_lt {k : Nat} (hk : k < 4) : T.ρ₂ rs₂ k < 4 * T.m₂ :=
  T.size₂ W ▸ RotationSystem.rot_lt W.emb₂.total W.emb₂.involution (T.virt₂_lt W hk)
theorem ρ₁_div_ne {k : Nat} (hk : k < 4) : T.ρ₁ rs₁ k / 4 ≠ T.e₁ :=
  W.deg₁ k hk _ (by rw [RotationSystem.get_eq_rot W.emb₁.total (T.virt₁_lt W hk)]; rfl)
theorem ρ₂_div_ne {k : Nat} (hk : k < 4) : T.ρ₂ rs₂ k / 4 ≠ T.e₂ :=
  W.deg₂ k hk _ (by rw [RotationSystem.get_eq_rot W.emb₂.total (T.virt₂_lt W hk)]; rfl)
theorem rot_ρ₁ {k : Nat} (hk : k < 4) : rs₁.rot (T.ρ₁ rs₁ k) = 4 * T.e₁ + k :=
  RotationSystem.rot_rot W.emb₁.total W.emb₁.involution (T.virt₁_lt W hk)
theorem rot_ρ₂ {k : Nat} (hk : k < 4) : rs₂.rot (T.ρ₂ rs₂ k) = 4 * T.e₂ + k :=
  RotationSystem.rot_rot W.emb₂.total W.emb₂.involution (T.virt₂_lt W hk)
theorem ρ₁_mod2 {k : Nat} (hk : k < 4) : T.ρ₁ rs₁ k % 2 ≠ k % 2 := by
  have := RotationSystem.rot_mod2 W.emb₁.total W.emb₁.opposite_dir (T.virt₁_lt W hk)
  rw [ρ₁_def]; omega
theorem ρ₂_mod2 {k : Nat} (hk : k < 4) : T.ρ₂ rs₂ k % 2 ≠ k % 2 := by
  have := RotationSystem.rot_mod2 W.emb₂.total W.emb₂.opposite_dir (T.virt₂_lt W hk)
  rw [ρ₂_def]; omega
theorem ρ₁_inj {k k' : Nat} (hk : k < 4) (hk' : k' < 4) (h : T.ρ₁ rs₁ k = T.ρ₁ rs₁ k') :
    k = k' := by
  have := RotationSystem.rot_inj W.emb₁.total W.emb₁.involution (T.virt₁_lt W hk)
    (T.virt₁_lt W hk') h
  omega
theorem ρ₂_inj {k k' : Nat} (hk : k < 4) (hk' : k' < 4) (h : T.ρ₂ rs₂ k = T.ρ₂ rs₂ k') :
    k = k' := by
  have := RotationSystem.rot_inj W.emb₂.total W.emb₂.involution (T.virt₂_lt W hk)
    (T.virt₂_lt W hk') h
  omega

theorem vert_virt₁ {k : Nat} (hk : k < 4) :
    QE.vert T.es₁ (4 * T.e₁ + k) = some (if k / 2 % 2 = 0 then T.u₁ else T.v₁) := by
  unfold QE.vert QE.edge QE.side
  have h1 : (4 * T.e₁ + k) / 4 = T.e₁ := by omega
  have h2 : (4 * T.e₁ + k) / 2 % 2 = k / 2 % 2 := by omega
  rw [h1, h2, W.e₁]; rfl
theorem vert_virt₂ {k : Nat} (hk : k < 4) :
    QE.vert T.es₂ (4 * T.e₂ + k) = some (if k / 2 % 2 = 0 then T.u₂ else T.v₂) := by
  unfold QE.vert QE.edge QE.side
  have h1 : (4 * T.e₂ + k) / 4 = T.e₂ := by omega
  have h2 : (4 * T.e₂ + k) / 2 % 2 = k / 2 % 2 := by omega
  rw [h1, h2, W.e₂]; rfl
theorem vert_ρ₁ {k : Nat} (hk : k < 4) :
    QE.vert T.es₁ (T.ρ₁ rs₁ k) = some (if k / 2 % 2 = 0 then T.u₁ else T.v₁) := by
  rw [ρ₁_def, RotationSystem.vert_rot W.emb₁.total W.emb₁.same_vertex (T.virt₁_lt W hk),
    T.vert_virt₁ W hk]
theorem vert_ρ₂ {k : Nat} (hk : k < 4) :
    QE.vert T.es₂ (T.ρ₂ rs₂ k) = some (if k / 2 % 2 = 0 then T.u₂ else T.v₂) := by
  rw [ρ₂_def, RotationSystem.vert_rot W.emb₂.total W.emb₂.same_vertex (T.virt₂_lt W hk),
    T.vert_virt₂ W hk]

/-! ### The glued rotation as a total function -/

theorem spliceFn_eq {r : Nat} (hr : r < 4 * (T.m₁ + T.m₂ - 2)) :
    T.spliceFn rs₁ rs₂ r = some (T.glue rs₁ rs₂ r) := by
  unfold spliceFn glue
  split
  · rename_i h
    rw [RotationSystem.get_eq_rot W.emb₁.total (T.size₁ W ▸ T.pre₁_lt W h), Option.bind_some]
    split
    · rw [RotationSystem.get_eq_rot W.emb₂.total
        (T.virt₂_lt W (xor_lt4 (Nat.mod_lt _ (by omega)) (by decide)))]
      rfl
    · rfl
  · rename_i h
    rw [RotationSystem.get_eq_rot W.emb₂.total (T.size₂ W ▸ T.pre₂_lt W hr (not_lt.1 h)), Option.bind_some]
    split
    · rw [RotationSystem.get_eq_rot W.emb₁.total
        (T.virt₁_lt W (xor_lt4 (Nat.mod_lt _ (by omega)) (by decide)))]
      rfl
    · rfl

omit W in
theorem splice_size : (T.splice rs₁ rs₂).size = 4 * (T.m₁ + T.m₂ - 2) := by
  simp [splice, RotationSystem.size]

omit W in
theorem splice_get {r : Nat} (hr : r < 4 * (T.m₁ + T.m₂ - 2)) :
    (T.splice rs₁ rs₂).get r = T.spliceFn rs₁ rs₂ r := by
  simp [splice, RotationSystem.get, List.getElem?_range hr]

theorem splice_get_eq {r : Nat} (hr : r < 4 * (T.m₁ + T.m₂ - 2)) :
    (T.splice rs₁ rs₂).get r = some (T.glue rs₁ rs₂ r) := by
  rw [T.splice_get hr, T.spliceFn_eq W hr]

theorem splice_rot {r : Nat} (hr : r < 4 * (T.m₁ + T.m₂ - 2)) :
    (T.splice rs₁ rs₂).rot r = T.glue rs₁ rs₂ r := by
  unfold RotationSystem.rot; rw [T.splice_get_eq W hr]; rfl

theorem glue_lt {r : Nat} (hr : r < 4 * (T.m₁ + T.m₂ - 2)) :
    T.glue rs₁ rs₂ r < 4 * (T.m₁ + T.m₂ - 2) := by
  have hm := T.e₁_lt W; have hm₂ := T.e₂_lt W
  have hk : ∀ s : Nat, s % 4 ^^^ 1 < 4 := fun s => xor_lt4 (Nat.mod_lt _ (by omega)) (by decide)
  unfold glue
  split
  · rename_i h
    have hq := T.pre₁_lt W h
    have hs := T.size₁ W ▸ RotationSystem.rot_lt W.emb₁.total W.emb₁.involution (T.size₁ W ▸ hq)
    split
    · exact T.qe₂_lt W (T.ρ₂_lt W (hk (rs₁.rot (T.pre₁ r)))) (T.ρ₂_div_ne W (hk (rs₁.rot (T.pre₁ r))))
    · have := T.qe₁_lt_off W hs ‹_›
      unfold off at this; omega
  · rename_i h
    have hq := T.pre₂_lt W hr (not_lt.1 h)
    have hs := T.size₂ W ▸ RotationSystem.rot_lt W.emb₂.total W.emb₂.involution (T.size₂ W ▸ hq)
    split
    · have := T.qe₁_lt_off W (T.ρ₁_lt W (hk (rs₂.rot (T.pre₂ r))))
        (T.ρ₁_div_ne W (hk (rs₂.rot (T.pre₂ r))))
      unfold off at this; omega
    · exact T.qe₂_lt W hs ‹_›

theorem glue_glue {r : Nat} (hr : r < 4 * (T.m₁ + T.m₂ - 2)) :
    T.glue rs₁ rs₂ (T.glue rs₁ rs₂ r) = r := by
  have ht₁ := W.emb₁.total; have hi₁ := W.emb₁.involution
  have ht₂ := W.emb₂.total; have hi₂ := W.emb₂.involution
  have hk : ∀ s : Nat, s % 4 ^^^ 1 < 4 := fun s => xor_lt4 (Nat.mod_lt _ (by omega)) (by decide)
  unfold glue
  split
  · rename_i h
    have hq := T.pre₁_lt W h
    have hs := T.size₁ W ▸ RotationSystem.rot_lt ht₁ hi₁ (T.size₁ W ▸ hq)
    split
    · rename_i he
      have hk1 := hk (rs₁.rot (T.pre₁ r))
      have h1 : ¬T.qe₂ (T.ρ₂ rs₂ (rs₁.rot (T.pre₁ r) % 4 ^^^ 1)) < T.off :=
        not_lt.2 (T.off_le_qe₂ _)
      rw [if_neg h1, T.pre₂_qe₂ (T.ρ₂_div_ne W hk1), T.rot_ρ₂ W hk1]
      have h2 : (4 * T.e₂ + (rs₁.rot (T.pre₁ r) % 4 ^^^ 1)) / 4 = T.e₂ := by omega
      have h3 : (4 * T.e₂ + (rs₁.rot (T.pre₁ r) % 4 ^^^ 1)) % 4 ^^^ 1 = rs₁.rot (T.pre₁ r) % 4 := by
        have : (4 * T.e₂ + (rs₁.rot (T.pre₁ r) % 4 ^^^ 1)) % 4 = rs₁.rot (T.pre₁ r) % 4 ^^^ 1 := by
          omega
        rw [this, xor_xor_self]
      rw [if_pos h2, h3, ρ₁_def]
      have h4 : 4 * T.e₁ + rs₁.rot (T.pre₁ r) % 4 = rs₁.rot (T.pre₁ r) := by omega
      rw [h4, RotationSystem.rot_rot ht₁ hi₁ (T.size₁ W ▸ hq), T.qe₁_pre₁]
    · rename_i he
      rw [if_pos (T.qe₁_lt_off W hs he), T.pre₁_qe₁ he, RotationSystem.rot_rot ht₁ hi₁ (T.size₁ W ▸ hq),
        if_neg (T.pre₁_div_ne r), T.qe₁_pre₁]
  · rename_i h
    have hq := T.pre₂_lt W hr (not_lt.1 h)
    have hs := T.size₂ W ▸ RotationSystem.rot_lt ht₂ hi₂ (T.size₂ W ▸ hq)
    split
    · rename_i he
      have hk1 := hk (rs₂.rot (T.pre₂ r))
      have h1 : T.qe₁ (T.ρ₁ rs₁ (rs₂.rot (T.pre₂ r) % 4 ^^^ 1)) < T.off :=
        T.qe₁_lt_off W (T.ρ₁_lt W hk1) (T.ρ₁_div_ne W hk1)
      rw [if_pos h1, T.pre₁_qe₁ (T.ρ₁_div_ne W hk1), T.rot_ρ₁ W hk1]
      have h2 : (4 * T.e₁ + (rs₂.rot (T.pre₂ r) % 4 ^^^ 1)) / 4 = T.e₁ := by omega
      have h3 : (4 * T.e₁ + (rs₂.rot (T.pre₂ r) % 4 ^^^ 1)) % 4 ^^^ 1 = rs₂.rot (T.pre₂ r) % 4 := by
        have : (4 * T.e₁ + (rs₂.rot (T.pre₂ r) % 4 ^^^ 1)) % 4 = rs₂.rot (T.pre₂ r) % 4 ^^^ 1 := by
          omega
        rw [this, xor_xor_self]
      rw [if_pos h2, h3, ρ₂_def]
      have h4 : 4 * T.e₂ + rs₂.rot (T.pre₂ r) % 4 = rs₂.rot (T.pre₂ r) := by omega
      rw [h4, RotationSystem.rot_rot ht₂ hi₂ (T.size₂ W ▸ hq), T.qe₂_pre₂ (not_lt.1 h)]
    · rename_i he
      rw [if_neg (not_lt.2 (T.off_le_qe₂ _)), T.pre₂_qe₂ he,
        RotationSystem.rot_rot ht₂ hi₂ (T.size₂ W ▸ hq), if_neg (T.pre₂_div_ne r),
        T.qe₂_pre₂ (not_lt.1 h)]

theorem glue_mod2 {r : Nat} (hr : r < 4 * (T.m₁ + T.m₂ - 2)) :
    T.glue rs₁ rs₂ r % 2 ≠ r % 2 := by
  have hk : ∀ s : Nat, s % 4 ^^^ 1 < 4 := fun s => xor_lt4 (Nat.mod_lt _ (by omega)) (by decide)
  have hx : ∀ s : Nat, (s % 4 ^^^ 1) % 2 ≠ s % 2 := fun s => by
    have := xor_mod2 (s % 4) (c := 1) rfl; omega
  unfold glue
  split
  · rename_i h
    have hq := T.pre₁_lt W h
    have h1 := RotationSystem.rot_mod2 W.emb₁.total W.emb₁.opposite_dir (T.size₁ W ▸ hq)
    have h2 : T.pre₁ r % 2 = r % 2 := by unfold pre₁; omega
    split
    · rw [T.qe₂_mod2]
      have h3 := T.ρ₂_mod2 W (hk (rs₁.rot (T.pre₁ r)))
      have h4 := hx (rs₁.rot (T.pre₁ r))
      omega
    · rw [T.qe₁_mod2]; omega
  · rename_i h
    have hq := T.pre₂_lt W hr (not_lt.1 h)
    have h1 := RotationSystem.rot_mod2 W.emb₂.total W.emb₂.opposite_dir (T.size₂ W ▸ hq)
    have h2 : T.pre₂ r % 2 = r % 2 := by unfold pre₂; unfold off at h ⊢; omega
    split
    · rw [T.qe₁_mod2]
      have h3 := T.ρ₁_mod2 W (hk (rs₂.rot (T.pre₂ r)))
      have h4 := hx (rs₂.rot (T.pre₂ r))
      omega
    · rw [T.qe₂_mod2]; omega

/-! ### Vertices of the 2-sum -/

theorem edges_length : T.edges.length = T.m₁ + T.m₂ - 2 := by
  have := T.e₁_lt W; have := T.e₂_lt W
  simp only [edges, List.length_append, List.length_map, List.length_eraseIdx]
  unfold m₁ m₂ at *
  rw [if_pos ‹T.e₁ < _›, if_pos ‹T.e₂ < _›]; omega

theorem length_eraseIdx₁ : (T.es₁.eraseIdx T.e₁).length = T.m₁ - 1 := by
  have := T.e₁_lt W; unfold m₁ at *; rw [List.length_eraseIdx, if_pos this]

theorem vert_qe₁ {q : Nat} (hq : q < 4 * T.m₁) (h : q / 4 ≠ T.e₁) :
    QE.vert T.edges (T.qe₁ q) = QE.vert T.es₁ q := by
  unfold QE.vert QE.edge QE.side
  rw [T.qe₁_div]
  have hlt : delIdx T.e₁ (q / 4) < (T.es₁.eraseIdx T.e₁).length := by
    rw [T.length_eraseIdx₁ W]; exact delIdx_lt (by unfold m₁ at hq ⊢; omega) h (T.e₁_lt W)
  have h2 : T.qe₁ q / 2 % 2 = q / 2 % 2 := by unfold qe₁; omega
  rw [edges, List.getElem?_append_left hlt, getElem?_eraseIdx_delIdx _ h, h2]

theorem vert_qe₂ {q : Nat} (_hq : q < 4 * T.m₂) (h : q / 4 ≠ T.e₂) :
    QE.vert T.edges (T.qe₂ q) = (QE.vert T.es₂ q).map T.vert₂ := by
  unfold QE.vert QE.edge QE.side
  have hd : T.qe₂ q / 4 = (T.m₁ - 1) + delIdx T.e₂ (q / 4) := by unfold qe₂ off; omega
  rw [hd, edges, List.getElem?_append_right (by rw [T.length_eraseIdx₁ W]; omega),
    T.length_eraseIdx₁ W, Nat.add_sub_cancel_left, List.getElem?_map,
    getElem?_eraseIdx_delIdx _ h]
  have h2 : T.qe₂ q / 2 % 2 = q / 2 % 2 := by unfold qe₂ off; omega
  rw [h2]
  cases T.es₂[q / 4]? with
  | none => rfl
  | some p => simp only [Option.map_some]; split <;> rfl

omit W in
theorem vert₂_u : T.vert₂ T.u₂ = T.u₁ := by simp [vert₂]
theorem vert₂_v : T.vert₂ T.v₂ = T.v₁ := by simp [vert₂, W.uv₂.symm]

theorem glue_vert {r : Nat} (hr : r < 4 * (T.m₁ + T.m₂ - 2)) :
    QE.vert T.edges (T.glue rs₁ rs₂ r) = QE.vert T.edges r := by
  have hk : ∀ s : Nat, s % 4 ^^^ 1 < 4 := fun s => xor_lt4 (Nat.mod_lt _ (by omega)) (by decide)
  have hside : ∀ s : Nat, (s % 4 ^^^ 1) / 2 % 2 = s / 2 % 2 := fun s => by
    have := @Nat.xor_div_two_pow (s % 4) 1 1
    simp only [Nat.pow_one] at this
    rw [this, show (1 : Nat) / 2 = 0 from rfl, Nat.xor_zero]
    omega
  unfold glue
  split
  · rename_i h
    have hq := T.pre₁_lt W h
    have hs := T.size₁ W ▸ RotationSystem.rot_lt W.emb₁.total W.emb₁.involution (T.size₁ W ▸ hq)
    have hv := RotationSystem.vert_rot W.emb₁.total W.emb₁.same_vertex (T.size₁ W ▸ hq)
    conv_rhs => rw [← T.qe₁_pre₁ r, T.vert_qe₁ W hq (T.pre₁_div_ne r), ← hv]
    split
    · rename_i he
      rw [T.vert_qe₂ W (T.ρ₂_lt W (hk _)) (T.ρ₂_div_ne W (hk _)), T.vert_ρ₂ W (hk _)]
      have h4 : rs₁.rot (T.pre₁ r) = 4 * T.e₁ + rs₁.rot (T.pre₁ r) % 4 := by omega
      conv_rhs => rw [h4, T.vert_virt₁ W (Nat.mod_lt _ (by omega))]
      have h5 : rs₁.rot (T.pre₁ r) % 4 / 2 % 2 = rs₁.rot (T.pre₁ r) / 2 % 2 := by omega
      rw [hside, h5]
      simp only [Option.map_some]
      split <;> simp [vert₂_u, T.vert₂_v W]
    · rename_i he
      rw [T.vert_qe₁ W hs he]
  · rename_i h
    have hq := T.pre₂_lt W hr (not_lt.1 h)
    have hs := T.size₂ W ▸ RotationSystem.rot_lt W.emb₂.total W.emb₂.involution (T.size₂ W ▸ hq)
    have hv := RotationSystem.vert_rot W.emb₂.total W.emb₂.same_vertex (T.size₂ W ▸ hq)
    conv_rhs => rw [← T.qe₂_pre₂ (not_lt.1 h), T.vert_qe₂ W hq (T.pre₂_div_ne r), ← hv]
    split
    · rename_i he
      rw [T.vert_qe₁ W (T.ρ₁_lt W (hk _)) (T.ρ₁_div_ne W (hk _)), T.vert_ρ₁ W (hk _)]
      have h4 : rs₂.rot (T.pre₂ r) = 4 * T.e₂ + rs₂.rot (T.pre₂ r) % 4 := by omega
      conv_rhs => rw [h4, T.vert_virt₂ W (Nat.mod_lt _ (by omega))]
      have h5 : rs₂.rot (T.pre₂ r) % 4 / 2 % 2 = rs₂.rot (T.pre₂ r) / 2 % 2 := by omega
      rw [hside, h5]
      simp only [Option.map_some]
      split <;> simp [vert₂_u, T.vert₂_v W]
    · rename_i he
      rw [T.vert_qe₂ W hs he]

theorem u₁_lt : T.u₁ < T.n₁ :=
  (W.emb₁.verts _ (List.getElem?_eq_some_iff.1 W.e₁ |>.2 ▸ List.getElem_mem _)).1
theorem v₁_lt : T.v₁ < T.n₁ :=
  (W.emb₁.verts _ (List.getElem?_eq_some_iff.1 W.e₁ |>.2 ▸ List.getElem_mem _)).2

theorem edges_verts : ∀ p ∈ T.edges, p.1 < T.nVerts ∧ p.2 < T.nVerts := by
  intro p hp
  unfold edges at hp
  rw [List.mem_append] at hp
  unfold nVerts
  rcases hp with hp | hp
  · have := W.emb₁.verts p (List.mem_of_mem_eraseIdx hp)
    omega
  · rw [List.mem_map] at hp
    obtain ⟨q, hq, rfl⟩ := hp
    have := W.emb₂.verts q (List.mem_of_mem_eraseIdx hq)
    have hu := T.u₁_lt W; have hv := T.v₁_lt W
    constructor <;> (unfold vert₂; split_ifs <;> omega)

/-- All embedding fields of the splice except the vertex-orbit count. -/
theorem splice_total : (T.splice rs₁ rs₂).Total := fun r hr => by
  rw [T.splice_get_eq W (T.splice_size ▸ hr)]; rfl
theorem splice_involution : (T.splice rs₁ rs₂).Involution := fun r hr s hs => by
  rw [T.splice_size] at hr
  rw [T.splice_get_eq W hr, Option.mem_def, Option.some.injEq] at hs
  subst hs
  rw [T.splice_size, T.splice_get_eq W (T.glue_lt W hr), T.glue_glue W hr]
  exact ⟨T.glue_lt W hr, rfl⟩
theorem splice_oppositeDir : (T.splice rs₁ rs₂).OppositeDir := fun r hr s hs => by
  rw [T.splice_size] at hr
  rw [T.splice_get_eq W hr, Option.mem_def, Option.some.injEq] at hs
  subst hs
  exact T.glue_mod2 W hr
theorem splice_sameVertex : (T.splice rs₁ rs₂).SameVertex T.edges := fun r hr s hs => by
  rw [T.splice_size] at hr
  rw [T.splice_get_eq W hr, Option.mem_def, Option.some.injEq] at hs
  subst hs
  exact (T.glue_vert W hr).symm

end WF

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
