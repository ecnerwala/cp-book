import Spqr.Proofs.PieceInsert

/-!
# One-sum of two capped pieces

`Capped.join`: two capped pieces meeting exactly in the vertex `v` (cap endpoints `u, v` and
`v, w`), glued by the node step's two links `c2 ↔ d1`, `d0 ↔ c3`, form a capped piece with ends
`c0 c1` at `u` and `d2 d3` at `w`. The witness is `(ρ.union σ).conj l2 (ρ.size + m0)`; `c0` and
`d2` stay cofacial because the two outer faces merge through the glued corner (both image swaps of
`conj_stepC` join orbits of distinct components).
-/

namespace Spqr

open Classical RotationSystem

namespace Piece

theorem Capped.join {P : Piece} {E : List Nat} {A A' : Array (Option Nat)} {ρ σ : RotationSystem}
    {c0 c1 c2 c3 d0 d1 d2 d3 u v w : Nat}
    (h : P.Capped A ρ c0 c1 c2 c3 u v)
    (hE : ({P with ves := E} : Piece).Capped A σ d0 d1 d2 d3 v w)
    (huv : u ≠ v) (hvw : v ≠ w)
    (hdis : List.Disjoint P.ves E)
    (hsep : ∀ x, HasEdge P.es x → HasEdge ({P with ves := E} : Piece).es x → x = v)
    (hA' : ∀ q, A'[q]? = if q = c2 then some (some d1) else if q = d1 then some (some c2)
      else if q = d0 then some (some c3) else if q = c3 then some (some d0) else A[q]?) :
    ∃ ρ', ({P with ves := P.ves ++ E} : Piece).Capped A' ρ' c0 c1 d2 d3 u w := by
  let U : Piece := {P with ves := P.ves ++ E}
  obtain ⟨l0, l1, hl0, hl1, h01, hvu⟩ := h.pair0
  obtain ⟨l2, l3, hl2, hl3, h23, hvv⟩ := h.pair2
  obtain ⟨l0', l2', hl0', hl2', hface⟩ := h.face
  rw [hl0] at hl0'; rw [hl2] at hl2'; cases hl0'; cases hl2'
  obtain ⟨m0, m1, hm0, hm1, hm01, hvv'⟩ := hE.pair0
  obtain ⟨m2, m3, hm2, hm3, hm23, hvw'⟩ := hE.pair2
  obtain ⟨m0', m2', hm0', hm2', hfaceE⟩ := hE.face
  rw [hm0] at hm0'; rw [hm2] at hm2'; cases hm0'; cases hm2'
  have hpl := h.planar; have hplE := hE.planar
  have hs₁ : ρ.size = 4 * P.ves.length := by rw [hpl.size, es, List.length_map]
  have hs₂ : σ.size = 4 * E.length := by rw [hplE.size, es, List.length_map]
  have hl0lt : l0 < ρ.size := by rw [hs₁]; exact loc_lt hl0
  have hl1lt : l1 < ρ.size := by rw [hs₁]; exact loc_lt hl1
  have hl2lt : l2 < ρ.size := by rw [hs₁]; exact loc_lt hl2
  have hl3lt : l3 < ρ.size := by rw [hs₁]; exact loc_lt hl3
  have hm0lt : m0 < σ.size := by rw [hs₂]; exact loc_lt hm0
  have hm1lt : m1 < σ.size := by rw [hs₂]; exact loc_lt hm1
  have hm2lt : m2 < σ.size := by rw [hs₂]; exact loc_lt hm2
  have hm3lt : m3 < σ.size := by rw [hs₂]; exact loc_lt hm3
  have h10 : ρ.get l1 = some l0 := (hpl.involution l0 hl0lt l1 h01).2
  have h32 : ρ.get l3 = some l2 := (hpl.involution l2 hl2lt l3 h23).2
  have hm10 : σ.get m1 = some m0 := (hplE.involution m0 hm0lt m1 hm01).2
  have hm32 : σ.get m3 = some m2 := (hplE.involution m2 hm2lt m3 hm23).2
  have hpl0 : l0 % 2 = 0 := by rw [loc_mod_two hl0]; exact h.dir0
  have hpl1 : l1 % 2 = 1 := by rw [loc_mod_two hl1]; exact h.dir1
  have hpl2 : l2 % 2 = 0 := by rw [loc_mod_two hl2]; exact h.dir2
  have hpl3 : l3 % 2 = 1 := by rw [loc_mod_two hl3]; exact h.dir3
  have hpm0 : m0 % 2 = 0 := by rw [loc_mod_two hm0]; exact hE.dir0
  have hpm1 : m1 % 2 = 1 := by rw [loc_mod_two hm1]; exact hE.dir1
  have hpm2 : m2 % 2 = 0 := by rw [loc_mod_two hm2]; exact hE.dir2
  have hpm3 : m3 % 2 = 1 := by rw [loc_mod_two hm3]; exact hE.dir3
  have hvu1 : QE.vert P.es l1 = some u := (hpl.same_vertex l0 hl0lt l1 h01).symm.trans hvu
  have hvv3 : QE.vert P.es l3 = some v := (hpl.same_vertex l2 hl2lt l3 h23).symm.trans hvv
  have hvv1 : QE.vert ({P with ves := E} : Piece).es m1 = some v :=
    (hplE.same_vertex m0 hm0lt m1 hm01).symm.trans hvv'
  have hvw3 : QE.vert ({P with ves := E} : Piece).es m3 = some w :=
    (hplE.same_vertex m2 hm2lt m3 hm23).symm.trans hvw'
  have hl02 : l0 ≠ l2 := fun e => by subst e; exact huv (Option.some.inj (hvu.symm.trans hvv))
  have hl13 : l1 ≠ l3 := fun e => by subst e; exact huv (Option.some.inj (hvu1.symm.trans hvv3))
  have hm02 : m0 ≠ m2 := fun e => by subst e; exact hvw (Option.some.inj (hvv'.symm.trans hvw'))
  have hm13 : m1 ≠ m3 := fun e => by subst e; exact hvw (Option.some.inj (hvv1.symm.trans hvw3))
  have hplan := hpl.splice hplE (by simpa [es] using loc_lt hl2) (by simpa [es] using loc_lt hm0)
    (by rw [hpl2, hpm0]) hvv hvv' hsep
  let R := ρ.union σ
  let D := ρ.size + m0
  let M1 := ρ.size + m1
  let M2 := ρ.size + m2
  let M3 := ρ.size + m3
  let ρ' := R.conj l2 D
  have hρ' : IsPlanarEmbedding U.es U.nVerts ρ' := by
    simpa only [U, es, List.map_append] using hplan
  have hnot : ∀ q, ({P with ves := E} : Piece).Mem q → QE.edge q ∉ P.ves := by
    intro q hq hh
    exact List.disjoint_left.1 hdis hh hq
  have hPE : ∀ q, P.Mem q → ∀ q', ({P with ves := E} : Piece).Mem q' → q ≠ q' :=
    fun q hq q' hq' heq => hdis hq (heq ▸ hq')
  have hmc0 := mem_of_loc hl0; have hmc1 := mem_of_loc hl1
  have hmc2 := mem_of_loc hl2; have hmc3 := mem_of_loc hl3
  have hmd0 := mem_of_loc hm0; have hmd1 := mem_of_loc hm1
  have hmd2 := mem_of_loc hm2; have hmd3 := mem_of_loc hm3
  have hl0U : U.loc c0 = some l0 := loc_append_left E hl0
  have hl1U : U.loc c1 = some l1 := loc_append_left E hl1
  have hl2U : U.loc c2 = some l2 := loc_append_left E hl2
  have hl3U : U.loc c3 = some l3 := loc_append_left E hl3
  have hDU : U.loc d0 = some D := by
    simpa only [D, hs₁] using loc_append_right P.ves (hnot d0 hmd0) hm0
  have hM1U : U.loc d1 = some M1 := by
    simpa only [M1, hs₁] using loc_append_right P.ves (hnot d1 hmd1) hm1
  have hM2U : U.loc d2 = some M2 := by
    simpa only [M2, hs₁] using loc_append_right P.ves (hnot d2 hmd2) hm2
  have hM3U : U.loc d3 = some M3 := by
    simpa only [M3, hs₁] using loc_append_right P.ves (hnot d3 hmd3) hm3
  have hUsize : R.size = ρ.size + σ.size := RotationSystem.union_size ..
  have hRs : R.size = 4 * (P.ves.length + E.length) := by rw [hUsize, hs₁, hs₂]; omega
  have hRt : R.Total := RotationSystem.union_total _ _ hpl.total hplE.total
  have hRi : R.Involution := RotationSystem.union_involution _ _ hpl.involution hplE.involution
  have hRo : R.OppositeDir :=
    RotationSystem.union_oppositeDir _ _ hs₁ hpl.opposite_dir hplE.opposite_dir
  have hR0 : R.get l0 = some l1 := by rw [RotationSystem.union_get_lt _ _ hl0lt]; exact h01
  have hR2 : R.get l2 = some l3 := by rw [RotationSystem.union_get_lt _ _ hl2lt]; exact h23
  have hR3 : R.get l3 = some l2 := by rw [RotationSystem.union_get_lt _ _ hl3lt]; exact h32
  have hRD : R.get D = some M1 := by
    rw [RotationSystem.union_get_ge _ _ (by dsimp [D]; omega)]
    simp only [D, M1, Nat.add_sub_cancel_left, hm01, Option.map_some]
  have hRM1 : R.get M1 = some D := by
    rw [RotationSystem.union_get_ge _ _ (by dsimp [M1]; omega)]
    simp only [D, M1, Nat.add_sub_cancel_left, hm10, Option.map_some]
  have hRM2 : R.get M2 = some M3 := by
    rw [RotationSystem.union_get_ge _ _ (by dsimp [M2]; omega)]
    simp only [M2, M3, Nat.add_sub_cancel_left, hm23, Option.map_some]
  have hUloc : ∀ q l, U.loc q = some l → l < R.size := by
    intro q l hl
    rw [hUsize, hs₁, hs₂]
    have := loc_lt hl
    simpa only [U, List.length_append, Nat.mul_add] using this
  have hl2R := hUloc _ _ hl2U; have hDR := hUloc _ _ hDU
  have hl2D : l2 ≠ D := by dsimp [D]; omega
  have hl3D : l3 ≠ D := by dsimp [D]; omega
  have hl0D : l0 ≠ D := by dsimp [D]; omega
  have hl1D : l1 ≠ D := by dsimp [D]; omega
  have hM1l2 : M1 ≠ l2 := by dsimp [M1]; omega
  have hM1D : M1 ≠ D := by dsimp [M1, D]; omega
  have hM2l2 : M2 ≠ l2 := by dsimp [M2]; omega
  have hM2D : M2 ≠ D := by dsimp [M2, D]; omega
  have hM3l2 : M3 ≠ l2 := by dsimp [M3]; omega
  have hM3D : M3 ≠ D := by dsimp [M3, D]; omega
  have hl32 : l3 ≠ l2 := by omega
  have hl12 : l1 ≠ l2 := by omega
  have hg2 : ρ'.get l2 = some M1 := by
    rw [RotationSystem.conj_get_lt _ _ _ hl2R]
    simp only [Equiv.swap_apply_left, hRD, Option.map_some,
      Equiv.swap_apply_of_ne_of_ne hM1l2 hM1D]
  have hgM1 : ρ'.get M1 = some l2 := by
    rw [RotationSystem.conj_get_lt _ _ _ (hUloc _ _ hM1U)]
    simp only [Equiv.swap_apply_of_ne_of_ne hM1l2 hM1D, hRM1, Option.map_some,
      Equiv.swap_apply_right]
  have hgD : ρ'.get D = some l3 := by
    rw [RotationSystem.conj_get_lt _ _ _ hDR]
    simp only [Equiv.swap_apply_right, hR2, Option.map_some,
      Equiv.swap_apply_of_ne_of_ne hl32 hl3D]
  have hg3 : ρ'.get l3 = some D := by
    rw [RotationSystem.conj_get_lt _ _ _ (hUloc _ _ hl3U)]
    simp only [Equiv.swap_apply_of_ne_of_ne hl32 hl3D, hR3, Option.map_some,
      Equiv.swap_apply_left]
  have hg0 : ρ'.get l0 = some l1 := by
    rw [RotationSystem.conj_get_lt _ _ _ (hUloc _ _ hl0U)]
    simp only [Equiv.swap_apply_of_ne_of_ne hl02 hl0D, hR0, Option.map_some,
      Equiv.swap_apply_of_ne_of_ne hl12 hl1D]
  have hgM2 : ρ'.get M2 = some M3 := by
    rw [RotationSystem.conj_get_lt _ _ _ (hUloc _ _ hM2U)]
    simp only [Equiv.swap_apply_of_ne_of_ne hM2l2 hM2D, hRM2, Option.map_some,
      Equiv.swap_apply_of_ne_of_ne hM3l2 hM3D]
  have hAg : U.Agrees A R := agrees_append hdis hs₁ h.agrees hE.agrees
  have hc20 : c2 ≠ c0 := fun e => hl02 (Option.some.inj (hl0.symm.trans (e ▸ hl2)))
  have hc31 : c3 ≠ c1 := fun e => hl13 (Option.some.inj (hl1.symm.trans (e ▸ hl3)))
  have hd02 : d0 ≠ d2 := fun e => hm02 (Option.some.inj ((e ▸ hm0).symm.trans hm2))
  have hd13 : d1 ≠ d3 := fun e => hm13 (Option.some.inj ((e ▸ hm1).symm.trans hm3))
  have hd0 := h.dir0; have hd1 := h.dir1; have hd2 := h.dir2; have hd3 := h.dir3
  have he0 := hE.dir0; have he1 := hE.dir1; have he2 := hE.dir2; have he3 := hE.dir3
  refine ⟨ρ', hρ', ?_, ?_, ⟨l0, l1, hl0U, hl1U, hg0, ?_⟩, ⟨M2, M3, hM2U, hM3U, hgM2, ?_⟩,
    ⟨l0, M2, hl0U, hM2U, ?_⟩, h.dir0, h.dir1, hE.dir2, hE.dir3⟩
  · intro q r hq hr
    rw [hA'] at hr
    split_ifs at hr with h1 h2 h3 h4
    · subst h1; cases hr; exact ⟨l2, M1, hl2U, hM1U, hg2⟩
    · subst h2; cases hr; exact ⟨M1, l2, hM1U, hl2U, hgM1⟩
    · subst h3; cases hr; exact ⟨D, l3, hDU, hl3U, hgD⟩
    · subst h4; cases hr; exact ⟨l3, D, hl3U, hDU, hg3⟩
    · obtain ⟨lq, lr, hlq, hlr, hqr⟩ := hAg q r hq hr
      have hqB : lq ≠ l2 := fun hh => h1 (loc_injective hlq (hh ▸ hl2U))
      have hqD : lq ≠ D := fun hh => h3 (loc_injective hlq (hh ▸ hDU))
      have hrB : lr ≠ l2 := by
        intro hh
        subst lr
        have := (hRi lq (hUloc q lq hlq) l2 hqr).2
        exact h4 (loc_injective hlq (Option.some.inj (this.symm.trans hR2) ▸ hl3U))
      have hrD : lr ≠ D := by
        intro hh
        subst lr
        have := (hRi lq (hUloc q lq hlq) D hqr).2
        exact h2 (loc_injective hlq (Option.some.inj (this.symm.trans hRD) ▸ hM1U))
      refine ⟨lq, lr, hlq, hlr, ?_⟩
      rw [RotationSystem.conj_get_lt _ _ _ (hUloc q lq hlq)]
      simp only [Equiv.swap_apply_of_ne_of_ne hqB hqD, hqr, Option.map_some,
        Equiv.swap_apply_of_ne_of_ne hrB hrD]
  · intro q hq
    rw [hA']
    split_ifs with h1 h2 h3 h4
    · subst h1
      exact iff_of_false (by simp) (by
        rintro (e | e | e | e)
        · exact hc20 e
        · omega
        · exact hPE _ hmc2 _ hmd2 e
        · exact hPE _ hmc2 _ hmd3 e)
    · subst h2
      exact iff_of_false (by simp) (by
        rintro (e | e | e | e)
        · exact hPE _ hmc0 _ hmd1 e.symm
        · exact hPE _ hmc1 _ hmd1 e.symm
        · omega
        · exact hd13 e)
    · subst h3
      exact iff_of_false (by simp) (by
        rintro (e | e | e | e)
        · exact hPE _ hmc0 _ hmd0 e.symm
        · exact hPE _ hmc1 _ hmd0 e.symm
        · exact hd02 e
        · omega)
    · subst h4
      exact iff_of_false (by simp) (by
        rintro (e | e | e | e)
        · omega
        · exact hc31 e
        · exact hPE _ hmc3 _ hmd2 e
        · exact hPE _ hmc3 _ hmd3 e)
    · rcases List.mem_append.1 hq with hqP | hqE
      · rw [h.unset q hqP]
        constructor
        · rintro (e | e | e | e)
          · exact Or.inl e
          · exact Or.inr (Or.inl e)
          · exact absurd e h1
          · exact absurd e h4
        · rintro (e | e | e | e)
          · exact Or.inl e
          · exact Or.inr (Or.inl e)
          · exact absurd e (hPE _ hqP _ hmd2)
          · exact absurd e (hPE _ hqP _ hmd3)
      · rw [hE.unset q hqE]
        constructor
        · rintro (e | e | e | e)
          · exact absurd e h3
          · exact absurd e h2
          · exact Or.inr (Or.inr (Or.inl e))
          · exact Or.inr (Or.inr (Or.inr e))
        · rintro (e | e | e | e)
          · exact absurd e.symm (hPE _ hmc0 _ hqE)
          · exact absurd e.symm (hPE _ hmc1 _ hqE)
          · exact Or.inr (Or.inr (Or.inl e))
          · exact Or.inr (Or.inr (Or.inr e))
  · rw [vert_loc hl0U, ← vert_loc hl0]; exact hvu
  · rw [vert_loc hM2U, ← vert_loc hm2]; exact hvw'
  · have hl0R := hUloc _ _ hl0U
    have hM2R := hUloc _ _ hM2U
    have hM1R := hUloc _ _ hM1U
    have hl3R := hUloc _ _ hl3U
    rw [RotationSystem.sameFaceOrbit_iff hρ'.total hρ'.involution
      (by rw [RotationSystem.conj_size, hRs]) (by rw [RotationSystem.conj_size]; exact hl0R)]
    have hrD : R.rot l2 ≠ D := by rw [rot_eq_of_get hR2]; exact hl3D
    rw [conj_stepC R l2 D hRt hRi hRo hl2R hDR hl2D hrD hRs (by decide), rot_eq_of_get hR2,
      rot_eq_of_get hRD]
    set f := R.stepC 3 with hfdef
    have hf : IsPermOn f (Finset.range R.size) := isPermOn_stepC hRt hRi hRs (by decide)
    have hfρ : IsPermOn (ρ.stepC 3) (Finset.range ρ.size) :=
      isPermOn_stepC hpl.total hpl.involution hs₁ (by decide)
    have hfU : f = unionStep (ρ.stepC 3) (shift ρ.size (σ.stepC 3)) (Finset.range ρ.size) :=
      union_stepC ρ σ hs₁ (by decide)
    have hpar : ∀ a b, SameOrbit f a b → a % 2 = b % 2 := fun a b hab =>
      sameOrbit_invariant (· % 2) (stepC_mod2 hRt hRo hRs (by decide) (by decide)) hab
    have hcomp : ∀ a, (f a < ρ.size ↔ a < ρ.size) := by
      intro a
      rw [hfU]
      unfold unionStep
      split_ifs with ha
      · rw [Finset.mem_range] at ha
        exact ⟨fun _ => ha, fun _ => Finset.mem_range.1 (hfρ.maps a (Finset.mem_range.2 ha))⟩
      · rw [Finset.mem_range] at ha
        unfold shift
        split_ifs
        omega
    have hcompS : ∀ a b, SameOrbit f a b → (a < ρ.size ↔ b < ρ.size) := by
      intro a b hab
      have := sameOrbit_invariant (fun q => decide (q < ρ.size))
        (fun a => by simp only [decide_eq_decide]; exact hcomp a) hab
      simpa using this
    have hmonoρ : ∀ a b, a < ρ.size → SameOrbit (ρ.stepC 3) a b → SameOrbit f a b := by
      intro a b ha hab
      refine sameOrbit_mono_of_eqOn (fun q => q < ρ.size) ?_ ?_ ha hab
      · intro q
        exact ⟨fun hq => Finset.mem_range.1 (hfρ.maps q (Finset.mem_range.2 hq)),
          fun hq => Finset.mem_range.1 (hfρ.mem_of_apply_mem (Finset.mem_range.2 hq))⟩
      · intro q hq
        rw [hfU]
        unfold unionStep
        simp only [Finset.mem_range.2 hq, ite_true]
    have hmonoσ : ∀ a b, SameOrbit (σ.stepC 3) a b →
        SameOrbit f (ρ.size + a) (ρ.size + b) := by
      intro a b hab
      have hsh := (sameOrbit_shift ρ.size (σ.stepC 3)).2 hab
      refine sameOrbit_mono_of_eqOn (fun q => ρ.size ≤ q) ?_ ?_ (by omega) hsh
      · intro q
        unfold shift
        split_ifs with hq <;> omega
      · intro q hq
        rw [hfU]
        unfold unionStep
        simp only [Finset.mem_range, not_lt.2 hq, ite_false]
    have hface' : SameOrbit (ρ.stepC 3) l0 l2 :=
      (RotationSystem.sameFaceOrbit_iff hpl.total hpl.involution hs₁ hl0lt).1 hface
    have hfaceE' : SameOrbit (σ.stepC 3) m0 m2 :=
      (RotationSystem.sameFaceOrbit_iff hplE.total hplE.involution hs₂ hm0lt).1 hfaceE
    have hxρ : l2 ^^^ 3 < ρ.size := by
      rw [hs₁]; exact xor_lt_mul4 (by rw [← hs₁]; exact hl2lt) (by decide)
    have hx'ρ : l3 ^^^ 3 < ρ.size := by
      rw [hs₁]; exact xor_lt_mul4 (by rw [← hs₁]; exact hl3lt) (by decide)
    have hxR : l2 ^^^ 3 < R.size := by rw [hUsize]; omega
    have hx'R : l3 ^^^ 3 < R.size := by rw [hUsize]; omega
    have hyge : ρ.size ≤ D ^^^ 3 := by
      rw [hs₁]; exact ge_of_xor_ge (by rw [← hs₁]; dsimp [D]; omega) (by decide)
    have hy'ge : ρ.size ≤ M1 ^^^ 3 := by
      rw [hs₁]; exact ge_of_xor_ge (by rw [← hs₁]; dsimp [M1]; omega) (by decide)
    have hyR : D ^^^ 3 < R.size := by
      rw [hRs]; exact xor_lt_mul4 (by rw [← hRs]; exact hDR) (by decide)
    have hy'R : M1 ^^^ 3 < R.size := by
      rw [hRs]; exact xor_lt_mul4 (by rw [← hRs]; exact hM1R) (by decide)
    have hxpar := xor_mod2 l2 (show 3 % 2 = 1 by decide)
    have hx'par := xor_mod2 l3 (show 3 % 2 = 1 by decide)
    have hfx' : f (l3 ^^^ 3) = l2 := by
      rw [hfdef, stepC_eq_rot hRt hRs (by decide) hx'R, xor_xor_self, rot_eq_of_get hR3]
    have hfy' : f (M1 ^^^ 3) = D := by
      rw [hfdef, stepC_eq_rot hRt hRs (by decide) hy'R, xor_xor_self, rot_eq_of_get hRM1]
    have hfx0 : SameOrbit f (l3 ^^^ 3) l0 :=
      (sameOrbit_of_eq hfx').trans (hmonoρ l2 l0 hl2lt hface'.symm)
    have hfy2 : SameOrbit f (M1 ^^^ 3) M2 :=
      (sameOrbit_of_eq hfy').trans (hmonoσ m0 m2 hfaceE')
    have hxy : ¬SameOrbit f (l2 ^^^ 3) (D ^^^ 3) := fun hh => by
      have := (hcompS _ _ hh).1 hxρ; omega
    have hx := Finset.mem_range.2 hxR
    have hy := Finset.mem_range.2 hyR
    have hg : IsPermOn (swapImg f (l2 ^^^ 3) (D ^^^ 3)) (Finset.range R.size) :=
      isPermOn_swapImg hf hx hy
    have hx'y' : ¬SameOrbit (swapImg f (l2 ^^^ 3) (D ^^^ 3)) (l3 ^^^ 3) (M1 ^^^ 3) := by
      rw [sameOrbit_swapImg hf hx hy hxy]
      rintro (hh | ⟨hh | hh, _⟩)
      · have := (hcompS _ _ hh).1 hx'ρ; omega
      · have := hpar _ _ hh; omega
      · have := (hcompS _ _ hh).2 hx'ρ; omega
    rw [sameOrbit_swapImg hg (Finset.mem_range.2 hx'R) (Finset.mem_range.2 hy'R) hx'y']
    exact Or.inr ⟨Or.inl ((sameOrbit_swapImg hf hx hy hxy).2 (Or.inl hfx0)),
      Or.inr ((sameOrbit_swapImg hf hx hy hxy).2 (Or.inl hfy2))⟩

end Piece

end Spqr
