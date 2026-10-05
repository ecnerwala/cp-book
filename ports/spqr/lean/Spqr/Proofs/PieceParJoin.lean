import Spqr.Proofs.PieceJoin
import Spqr.Proofs.PlanarSplice2
import Spqr.Proofs.OrbitSegment

/-!
# Parallel join of two capped pieces

`Capped.parJoin`: two capped pieces with the same cap endpoints `u, v`, meeting exactly in `u` and
`v`, glued by the P-node step's two links `c1 ↔ d0` (at `u`) and `c2 ↔ d3` (at `v`), form a capped
piece with pairs `c0 d1` at `u` and `d2 c3` at `v`.  The witness is
`((ρ.union σ).conj l1 (ρ.size + m1)).conj l2 (ρ.size + m2)`: the first conjugation merges the two
outer faces of each parity, the second splits them again (`IsPlanarEmbedding.splice2`).  `c0` and
`d2` end up on the same face because the face walk from `l0` reaches `l3 ^^^ 3` before any point
moved by the two conjugations, and that point is then sent to `ρ.size + m2`.
-/

namespace Spqr

open Classical RotationSystem

namespace Piece

theorem Capped.parJoin {P : Piece} {E : List Nat} {A A' : Array (Option Nat)}
    {ρ σ : RotationSystem} {c0 c1 c2 c3 d0 d1 d2 d3 u v : Nat}
    (h : P.Capped A ρ c0 c1 c2 c3 u v)
    (hE : ({P with ves := E} : Piece).Capped A σ d0 d1 d2 d3 u v)
    (huv : u ≠ v)
    (hdis : List.Disjoint P.ves E)
    (hsep : ∀ x, HasEdge P.es x → HasEdge ({P with ves := E} : Piece).es x → x = u ∨ x = v)
    (hA' : ∀ q, A'[q]? = if q = c1 then some (some d0) else if q = d0 then some (some c1)
      else if q = c2 then some (some d3) else if q = d3 then some (some c2) else A[q]?) :
    ∃ ρ', ({P with ves := P.ves ++ E} : Piece).Capped A' ρ' c0 d1 d2 c3 u v := by
  let U : Piece := {P with ves := P.ves ++ E}
  obtain ⟨l0, l1, hl0, hl1, h01, hvu⟩ := h.pair0
  obtain ⟨l2, l3, hl2, hl3, h23, hvv⟩ := h.pair2
  obtain ⟨l0', l2', hl0', hl2', hface⟩ := h.face
  rw [hl0] at hl0'; rw [hl2] at hl2'; cases hl0'; cases hl2'
  obtain ⟨m0, m1, hm0, hm1, hm01, hvu'⟩ := hE.pair0
  obtain ⟨m2, m3, hm2, hm3, hm23, hvv'⟩ := hE.pair2
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
  have hvu1' : QE.vert ({P with ves := E} : Piece).es m1 = some u :=
    (hplE.same_vertex m0 hm0lt m1 hm01).symm.trans hvu'
  have hvv3' : QE.vert ({P with ves := E} : Piece).es m3 = some v :=
    (hplE.same_vertex m2 hm2lt m3 hm23).symm.trans hvv'
  have hl02 : l0 ≠ l2 := fun e => by subst e; exact huv (Option.some.inj (hvu.symm.trans hvv))
  have hl13 : l1 ≠ l3 := fun e => by subst e; exact huv (Option.some.inj (hvu1.symm.trans hvv3))
  have hm02 : m0 ≠ m2 := fun e => by subst e; exact huv (Option.some.inj (hvu'.symm.trans hvv'))
  have hm13 : m1 ≠ m3 := fun e => by
    subst e; exact huv (Option.some.inj (hvu1'.symm.trans hvv3'))
  have hface' : SameOrbit (ρ.stepC 3) l0 l2 :=
    (RotationSystem.sameFaceOrbit_iff hpl.total hpl.involution hs₁ hl0lt).1 hface
  have hfaceE' : SameOrbit (σ.stepC 3) m0 m2 :=
    (RotationSystem.sameFaceOrbit_iff hplE.total hplE.involution hs₂ hm0lt).1 hfaceE
  have hface13 : SameOrbit (ρ.stepC 3) l1 l3 := by
    have := sameOrbit_stepC3_rot hpl.total hpl.involution hs₁ hl0lt hface'
    rwa [rot_eq_of_get h01, rot_eq_of_get h23] at this
  have hfaceE13 : SameOrbit (σ.stepC 3) m1 m3 := by
    have := sameOrbit_stepC3_rot hplE.total hplE.involution hs₂ hm0lt hfaceE'
    rwa [rot_eq_of_get hm01, rot_eq_of_get hm23] at this
  have hconn₁ : EdgesConn P.es u v :=
    edgesConn_of_sameOrbit hpl.total hpl.involution hpl.same_vertex hpl.size hface' hvu hvv
  have hconn₂ : EdgesConn ({P with ves := E} : Piece).es u v :=
    edgesConn_of_sameOrbit hplE.total hplE.involution hplE.same_vertex hplE.size hfaceE' hvu' hvv'
  let R := ρ.union σ
  let M0 := ρ.size + m0
  let M1 := ρ.size + m1
  let M2 := ρ.size + m2
  let M3 := ρ.size + m3
  let R1 := R.conj l1 M1
  let ρ' := R1.conj l2 M2
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
  have hM0U : U.loc d0 = some M0 := by
    simpa only [M0, hs₁] using loc_append_right P.ves (hnot d0 hmd0) hm0
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
  have hR1 : R.get l1 = some l0 := by rw [RotationSystem.union_get_lt _ _ hl1lt]; exact h10
  have hR2 : R.get l2 = some l3 := by rw [RotationSystem.union_get_lt _ _ hl2lt]; exact h23
  have hR3 : R.get l3 = some l2 := by rw [RotationSystem.union_get_lt _ _ hl3lt]; exact h32
  have hRM0 : R.get M0 = some M1 := by
    rw [RotationSystem.union_get_ge _ _ (by dsimp [M0]; omega)]
    simp only [M0, M1, Nat.add_sub_cancel_left, hm01, Option.map_some]
  have hRM1 : R.get M1 = some M0 := by
    rw [RotationSystem.union_get_ge _ _ (by dsimp [M1]; omega)]
    simp only [M0, M1, Nat.add_sub_cancel_left, hm10, Option.map_some]
  have hRM2 : R.get M2 = some M3 := by
    rw [RotationSystem.union_get_ge _ _ (by dsimp [M2]; omega)]
    simp only [M2, M3, Nat.add_sub_cancel_left, hm23, Option.map_some]
  have hRM3 : R.get M3 = some M2 := by
    rw [RotationSystem.union_get_ge _ _ (by dsimp [M3]; omega)]
    simp only [M2, M3, Nat.add_sub_cancel_left, hm32, Option.map_some]
  have hUloc : ∀ q l, U.loc q = some l → l < R.size := by
    intro q l hl
    rw [hUsize, hs₁, hs₂]
    have := loc_lt hl
    simpa only [U, List.length_append, Nat.mul_add] using this
  have hl0R := hUloc _ _ hl0U; have hl1R := hUloc _ _ hl1U
  have hl2R := hUloc _ _ hl2U; have hl3R := hUloc _ _ hl3U
  have hM0R := hUloc _ _ hM0U; have hM1R := hUloc _ _ hM1U
  have hM2R := hUloc _ _ hM2U; have hM3R := hUloc _ _ hM3U
  have hM0p : M0 % 2 = 0 := by dsimp [M0]; omega
  have hM1p : M1 % 2 = 1 := by dsimp [M1]; omega
  have hM2p : M2 % 2 = 0 := by dsimp [M2]; omega
  have hM3p : M3 % 2 = 1 := by dsimp [M3]; omega
  have hM0ge : ρ.size ≤ M0 := by dsimp [M0]; omega
  have hM1ge : ρ.size ≤ M1 := by dsimp [M1]; omega
  have hM2ge : ρ.size ≤ M2 := by dsimp [M2]; omega
  have hM3ge : ρ.size ≤ M3 := by dsimp [M3]; omega
  have hM0M2 : M0 ≠ M2 := by dsimp [M0, M2]; omega
  have hM1M3 : M1 ≠ M3 := by dsimp [M1, M3]; omega
  have hl01 : l0 ≠ l1 := by omega
  have hl03 : l0 ≠ l3 := by omega
  have hl12 : l1 ≠ l2 := by omega
  have hl21 : l2 ≠ l1 := by omega
  have hl31 : l3 ≠ l1 := by omega
  have hl32 : l3 ≠ l2 := by omega
  have hl0M1 : l0 ≠ M1 := by omega
  have hl0M2 : l0 ≠ M2 := by omega
  have hl1M1 : l1 ≠ M1 := by omega
  have hl1M2 : l1 ≠ M2 := by omega
  have hl2M1 : l2 ≠ M1 := by omega
  have hl2M2 : l2 ≠ M2 := by omega
  have hl3M1 : l3 ≠ M1 := by omega
  have hl3M2 : l3 ≠ M2 := by omega
  have hM0l1 : M0 ≠ l1 := by omega
  have hM0l2 : M0 ≠ l2 := by omega
  have hM0M1 : M0 ≠ M1 := by omega
  have hM1l2 : M1 ≠ l2 := by omega
  have hM1M2 : M1 ≠ M2 := by omega
  have hM2l1 : M2 ≠ l1 := by omega
  have hM2M1 : M2 ≠ M1 := by omega
  have hM3l1 : M3 ≠ l1 := by omega
  have hM3l2 : M3 ≠ l2 := by omega
  have hM3M1 : M3 ≠ M1 := by omega
  have hM3M2 : M3 ≠ M2 := by omega
  -- the intermediate system `R1 = R.conj l1 M1`
  have hR1s : R1.size = R.size := RotationSystem.conj_size ..
  have hR1s' : R1.size = 4 * (P.ves.length + E.length) := by rw [hR1s, hRs]
  have hR1g : ∀ q r, q < R.size → R.get q = some r → q ≠ l1 → q ≠ M1 → r ≠ l1 → r ≠ M1 →
      R1.get q = some r := by
    intro q r hq hqr hq1 hqM hr1 hrM
    rw [RotationSystem.conj_get_lt _ _ _ hq]
    simp only [Equiv.swap_apply_of_ne_of_ne hq1 hqM, hqr, Option.map_some,
      Equiv.swap_apply_of_ne_of_ne hr1 hrM]
  have hR1l1 : R1.get l1 = some M0 := by
    rw [RotationSystem.conj_get_lt _ _ _ hl1R]
    simp only [Equiv.swap_apply_left, hRM1, Option.map_some,
      Equiv.swap_apply_of_ne_of_ne hM0l1 hM0M1]
  have hR1M0 : R1.get M0 = some l1 := by
    rw [RotationSystem.conj_get_lt _ _ _ hM0R]
    simp only [Equiv.swap_apply_of_ne_of_ne hM0l1 hM0M1, hRM0, Option.map_some,
      Equiv.swap_apply_right]
  have hR1l0 : R1.get l0 = some M1 := by
    rw [RotationSystem.conj_get_lt _ _ _ hl0R]
    simp only [Equiv.swap_apply_of_ne_of_ne hl01 hl0M1, hR0, Option.map_some,
      Equiv.swap_apply_left]
  have hR1M1 : R1.get M1 = some l0 := by
    rw [RotationSystem.conj_get_lt _ _ _ hM1R]
    simp only [Equiv.swap_apply_right, hR1, Option.map_some,
      Equiv.swap_apply_of_ne_of_ne hl01 hl0M1]
  have hR1l2 : R1.get l2 = some l3 := hR1g _ _ hl2R hR2 hl21 hl2M1 hl31 hl3M1
  have hR1l3 : R1.get l3 = some l2 := hR1g _ _ hl3R hR3 hl31 hl3M1 hl21 hl2M1
  have hR1M2 : R1.get M2 = some M3 := hR1g _ _ hM2R hRM2 hM2l1 hM2M1 hM3l1 hM3M1
  have hR1M3 : R1.get M3 = some M2 := hR1g _ _ hM3R hRM3 hM3l1 hM3M1 hM2l1 hM2M1
  have hR1t : R1.Total := RotationSystem.conj_total R l1 M1 hRt hl1R hM1R
  have hR1i : R1.Involution := RotationSystem.conj_involution R l1 M1 hRt hRi hl1R hM1R
  have hR1o : R1.OppositeDir :=
    RotationSystem.conj_oppositeDir R l1 M1 hRt hRo hl1R hM1R (by omega)
  -- the final system `ρ' = R1.conj l2 M2`
  have hρ's : ρ'.size = R.size := by rw [RotationSystem.conj_size, hR1s]
  have hg' : ∀ q r, q < R.size → R1.get q = some r → q ≠ l2 → q ≠ M2 → r ≠ l2 → r ≠ M2 →
      ρ'.get q = some r := by
    intro q r hq hqr hq1 hqM hr1 hrM
    rw [RotationSystem.conj_get_lt _ _ _ (by rw [hR1s]; exact hq)]
    simp only [Equiv.swap_apply_of_ne_of_ne hq1 hqM, hqr, Option.map_some,
      Equiv.swap_apply_of_ne_of_ne hr1 hrM]
  have hg2 : ρ'.get l2 = some M3 := by
    rw [RotationSystem.conj_get_lt _ _ _ (by rw [hR1s]; exact hl2R)]
    simp only [Equiv.swap_apply_left, hR1M2, Option.map_some,
      Equiv.swap_apply_of_ne_of_ne hM3l2 hM3M2]
  have hgM3 : ρ'.get M3 = some l2 := by
    rw [RotationSystem.conj_get_lt _ _ _ (by rw [hR1s]; exact hM3R)]
    simp only [Equiv.swap_apply_of_ne_of_ne hM3l2 hM3M2, hR1M3, Option.map_some,
      Equiv.swap_apply_right]
  have hgM2 : ρ'.get M2 = some l3 := by
    rw [RotationSystem.conj_get_lt _ _ _ (by rw [hR1s]; exact hM2R)]
    simp only [Equiv.swap_apply_right, hR1l2, Option.map_some,
      Equiv.swap_apply_of_ne_of_ne hl32 hl3M2]
  have hg3 : ρ'.get l3 = some M2 := by
    rw [RotationSystem.conj_get_lt _ _ _ (by rw [hR1s]; exact hl3R)]
    simp only [Equiv.swap_apply_of_ne_of_ne hl32 hl3M2, hR1l3, Option.map_some,
      Equiv.swap_apply_left]
  have hg1 : ρ'.get l1 = some M0 := hg' _ _ hl1R hR1l1 hl12 hl1M2 hM0l2 hM0M2
  have hgM0 : ρ'.get M0 = some l1 := hg' _ _ hM0R hR1M0 hM0l2 hM0M2 hl12 hl1M2
  have hg0 : ρ'.get l0 = some M1 := hg' _ _ hl0R hR1l0 hl02 hl0M2 hM1l2 hM1M2
  have hgM1 : ρ'.get M1 = some l0 := hg' _ _ hM1R hR1M1 hM1l2 hM1M2 hl02 hl0M2
  -- face orbits of `R`
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
  have hfρq : ∀ q, q < ρ.size → f q = ρ.stepC 3 q := by
    intro q hq
    rw [hfU]
    unfold unionStep
    simp only [Finset.mem_range.2 hq, ite_true]
  have hxlt : ∀ q, q < ρ.size → q ^^^ 3 < ρ.size := fun q hq => by
    rw [hs₁]; exact xor_lt_mul4 (by rw [← hs₁]; exact hq) (by decide)
  have hxR : ∀ q, q < R.size → q ^^^ 3 < R.size := fun q hq => by
    rw [hRs]; exact xor_lt_mul4 (by rw [← hRs]; exact hq) (by decide)
  have hxge : ∀ q, ρ.size ≤ q → ρ.size ≤ q ^^^ 3 := fun q hq => by
    rw [hs₁]; exact ge_of_xor_ge (by rw [← hs₁]; exact hq) (by decide)
  have hfv : ∀ q r, q < R.size → R.get q = some r → f (q ^^^ 3) = r := fun q r hq hqr => by
    rw [hfdef, stepC_eq_rot hRt hRs (by decide) (hxR q hq), xor_xor_self, rot_eq_of_get hqr]
  have hx0 := xor_mod2 l0 (show 3 % 2 = 1 by decide)
  have hx1 := xor_mod2 l1 (show 3 % 2 = 1 by decide)
  have hx2 := xor_mod2 l2 (show 3 % 2 = 1 by decide)
  have hx3 := xor_mod2 l3 (show 3 % 2 = 1 by decide)
  have hxM0 := xor_mod2 M0 (show 3 % 2 = 1 by decide)
  have hxM1 := xor_mod2 M1 (show 3 % 2 = 1 by decide)
  have hxM2 := xor_mod2 M2 (show 3 % 2 = 1 by decide)
  have hxM3 := xor_mod2 M3 (show 3 % 2 = 1 by decide)
  have hgeM0 := hxge M0 hM0ge; have hgeM1 := hxge M1 hM1ge
  have hgeM2 := hxge M2 hM2ge; have hgeM3 := hxge M3 hM3ge
  have hltl0 := hxlt l0 hl0lt; have hltl1 := hxlt l1 hl1lt
  have hltl2 := hxlt l2 hl2lt; have hltl3 := hxlt l3 hl3lt
  -- the four outer face orbits of `R`
  have hOe1 : SameOrbit f (l1 ^^^ 3) (l3 ^^^ 3) :=
    ((sameOrbit_of_eq (hfv l1 l0 hl1R hR1)).trans (hmonoρ l0 l2 hl0lt hface')).trans
      (sameOrbit_of_eq (hfv l3 l2 hl3R hR3)).symm
  have hOo1 : SameOrbit f (l0 ^^^ 3) (l2 ^^^ 3) :=
    ((sameOrbit_of_eq (hfv l0 l1 hl0R hR0)).trans (hmonoρ l1 l3 hl1lt hface13)).trans
      (sameOrbit_of_eq (hfv l2 l3 hl2R hR2)).symm
  have hOe2 : SameOrbit f (M1 ^^^ 3) (M3 ^^^ 3) :=
    ((sameOrbit_of_eq (hfv M1 M0 hM1R hRM1)).trans (hmonoσ m0 m2 hfaceE')).trans
      (sameOrbit_of_eq (hfv M3 M2 hM3R hRM3)).symm
  have hOo2 : SameOrbit f (M0 ^^^ 3) (M2 ^^^ 3) :=
    ((sameOrbit_of_eq (hfv M0 M1 hM0R hRM0)).trans (hmonoσ m1 m3 hfaceE13)).trans
      (sameOrbit_of_eq (hfv M2 M3 hM2R hRM2)).symm
  -- the first conjugation merges `Oe1` with `Oe2` and `Oo1` with `Oo2`
  have hx := Finset.mem_range.2 (hxR l1 hl1R)
  have hy := Finset.mem_range.2 (hxR M1 hM1R)
  have hx' := Finset.mem_range.2 (hxR l0 hl0R)
  have hy' := Finset.mem_range.2 (hxR M0 hM0R)
  have hxy : ¬SameOrbit f (l1 ^^^ 3) (M1 ^^^ 3) := fun hh => by
    have := (hcompS _ _ hh).1 hltl1; omega
  have hg1p : IsPermOn (swapImg f (l1 ^^^ 3) (M1 ^^^ 3)) (Finset.range R.size) :=
    isPermOn_swapImg hf hx hy
  have hx'y' : ¬SameOrbit (swapImg f (l1 ^^^ 3) (M1 ^^^ 3)) (l0 ^^^ 3) (M0 ^^^ 3) := by
    rw [sameOrbit_swapImg hf hx hy hxy]
    rintro (hh | ⟨hh | hh, _⟩)
    · have := (hcompS _ _ hh).1 hltl0; omega
    · have := hpar _ _ hh; omega
    · have := (hcompS _ _ hh).2 hltl0; omega
  have hgdef : R1.stepC 3 =
      swapImg (swapImg f (l1 ^^^ 3) (M1 ^^^ 3)) (l0 ^^^ 3) (M0 ^^^ 3) := by
    rw [conj_stepC R l1 M1 hRt hRi hRo hl1R hM1R hl1M1
      (by rw [rot_eq_of_get hR1]; exact hl0M1) hRs (by decide), rot_eq_of_get hR1,
      rot_eq_of_get hRM1]
  have h1 : SameOrbit (R1.stepC 3) (l2 ^^^ 3) (M2 ^^^ 3) := by
    rw [hgdef, sameOrbit_swapImg hg1p hx' hy' hx'y']
    exact Or.inr ⟨Or.inl ((sameOrbit_swapImg hf hx hy hxy).2 (Or.inl hOo1)),
      Or.inr ((sameOrbit_swapImg hf hx hy hxy).2 (Or.inl hOo2))⟩
  have h2 : SameOrbit (R1.stepC 3) (R1.rot l2 ^^^ 3) (R1.rot M2 ^^^ 3) := by
    rw [rot_eq_of_get hR1l2, rot_eq_of_get hR1M2, hgdef, sameOrbit_swapImg hg1p hx' hy' hx'y']
    exact Or.inl ((sameOrbit_swapImg hf hx hy hxy).2 (Or.inr ⟨Or.inl hOe1, Or.inr hOe2⟩))
  have hplan := hpl.splice2 hplE (a := l1) (b := m1) (c := l2) (d := m2)
    (by simpa [es] using loc_lt hl1) (by simpa [es] using loc_lt hm1) (by rw [hpl1, hpm1])
    (by simpa [es] using loc_lt hl2) (by simpa [es] using loc_lt hm2) (by rw [hpl2, hpm2])
    hvu1 hvu1' hvv hvv' huv hconn₁ hconn₂ hsep h1 h2
  have hρ' : IsPlanarEmbedding U.es U.nVerts ρ' := by
    simpa only [U, es, List.map_append] using hplan
  -- the face of `l0` after the second conjugation
  have hh' : ρ'.stepC 3 =
      swapImg (swapImg (R1.stepC 3) (l2 ^^^ 3) (M2 ^^^ 3)) (l3 ^^^ 3) (M3 ^^^ 3) := by
    rw [conj_stepC R1 l2 M2 hR1t hR1i hR1o (by rw [hR1s]; exact hl2R) (by rw [hR1s]; exact hM2R)
      hl2M2 (by rw [rot_eq_of_get hR1l2]; exact hl3M2) hR1s' (by decide),
      rot_eq_of_get hR1l2, rot_eq_of_get hR1M2]
  have hF32 : SameOrbit (ρ.stepC 3) l0 (l3 ^^^ 3) := by
    refine hface'.trans (sameOrbit_of_eq ?_).symm
    rw [stepC_eq_rot hpl.total hs₁ (by decide) hltl3, xor_xor_self, rot_eq_of_get h32]
  have hFx : ρ.stepC 3 (l1 ^^^ 3) = l0 := by
    rw [stepC_eq_rot hpl.total hs₁ (by decide) hltl1, xor_xor_self, rot_eq_of_get h10]
  obtain ⟨k, hk, hkmin⟩ := hfρ.exists_first_hit (Finset.mem_range.2 hl0lt) hF32
  have hiter : (ρ'.stepC 3)^[k] l0 = (ρ.stepC 3)^[k] l0 := by
    apply iterate_eq_of_eqOn
    intro j hj
    set q := (ρ.stepC 3)^[j] l0 with hqdef
    have hq : q < ρ.size := Finset.mem_range.1 (hfρ.iterate_mem j (Finset.mem_range.2 hl0lt))
    have hqpar : q % 2 = 0 := by
      have := sameOrbit_invariant (· % 2)
        (stepC_mod2 hpl.total hpl.opposite_dir hs₁ (by decide) (by decide))
        (IsPermOn.sameOrbit_of_reach (IsPermOn.reach_iff.2 ⟨j, hqdef.symm⟩) :
          SameOrbit (ρ.stepC 3) l0 q)
      omega
    have hq3 : q ≠ l3 ^^^ 3 := (hkmin j hj).1
    have hqx : q ≠ l1 ^^^ 3 := fun e =>
      (hkmin j hj).2 (by rw [Function.iterate_succ_apply', ← hqdef, e, hFx])
    rw [hh', swapImg_of_ne _ hq3 (by omega), swapImg_of_ne _ (by omega) (by omega), hgdef,
      swapImg_of_ne _ (by omega) (by omega), swapImg_of_ne _ hqx (by omega), hfρq q hq]
  have hreach : SameOrbit (ρ'.stepC 3) l0 (l3 ^^^ 3) :=
    IsPermOn.sameOrbit_of_reach (IsPermOn.reach_iff.2 ⟨k, hiter.trans hk⟩)
  have hlast : ρ'.stepC 3 (l3 ^^^ 3) = M2 := by
    have hM13 : M3 ^^^ 3 ≠ M1 ^^^ 3 := fun e => by
      have := congrArg (· ^^^ 3) e
      simp only [xor_xor_self] at this
      exact hM3M1 this
    rw [hh', swapImg_left, swapImg_of_ne _ (by omega) (by omega), hgdef,
      swapImg_of_ne _ (by omega) (by omega), swapImg_of_ne _ (by omega) hM13,
      hfv M3 M2 hM3R hRM3]
  have hfaceU : SameOrbit (ρ'.stepC 3) l0 M2 := by
    have := SameOrbit.step (f := ρ'.stepC 3) (l3 ^^^ 3)
    rw [hlast] at this
    exact hreach.trans this
  -- assembly
  have hAg : U.Agrees A R := agrees_append hdis hs₁ h.agrees hE.agrees
  have hc20 : c2 ≠ c0 := fun e => hl02 (Option.some.inj (hl0.symm.trans (e ▸ hl2)))
  have hc31 : c3 ≠ c1 := fun e => hl13 (Option.some.inj (hl1.symm.trans (e ▸ hl3)))
  have hd02 : d0 ≠ d2 := fun e => hm02 (Option.some.inj ((e ▸ hm0).symm.trans hm2))
  have hd13 : d1 ≠ d3 := fun e => hm13 (Option.some.inj ((e ▸ hm1).symm.trans hm3))
  have hd0 := h.dir0; have hd1 := h.dir1; have hd2 := h.dir2; have hd3 := h.dir3
  have he0 := hE.dir0; have he1 := hE.dir1; have he2 := hE.dir2; have he3 := hE.dir3
  refine ⟨ρ', hρ', ?_, ?_, ⟨l0, M1, hl0U, hM1U, hg0, ?_⟩, ⟨M2, l3, hM2U, hl3U, hgM2, ?_⟩,
    ⟨l0, M2, hl0U, hM2U, ?_⟩, h.dir0, hE.dir1, hE.dir2, h.dir3⟩
  · intro q r hq hr
    rw [hA'] at hr
    split_ifs at hr with h1 h2 h3 h4
    · subst h1; cases hr; exact ⟨l1, M0, hl1U, hM0U, hg1⟩
    · subst h2; cases hr; exact ⟨M0, l1, hM0U, hl1U, hgM0⟩
    · subst h3; cases hr; exact ⟨l2, M3, hl2U, hM3U, hg2⟩
    · subst h4; cases hr; exact ⟨M3, l2, hM3U, hl2U, hgM3⟩
    · obtain ⟨lq, lr, hlq, hlr, hqr⟩ := hAg q r hq hr
      have hqlt := hUloc q lq hlq
      have back : ∀ s t, R.get s = some t → lr = s → lq = t := fun s t hst e => by
        subst e
        have := (hRi lq hqlt lr hqr).2
        exact Option.some.inj (this.symm.trans hst)
      have hq1 : lq ≠ l1 := fun hh => h1 (loc_injective hlq (hh ▸ hl1U))
      have hqM1 : lq ≠ M1 := fun hh => by
        rw [loc_injective hlq (hh ▸ hM1U), (hE.unset d1 hmd1).2 (Or.inr (Or.inl rfl))] at hr
        cases hr
      have hq2 : lq ≠ l2 := fun hh => h3 (loc_injective hlq (hh ▸ hl2U))
      have hqM2 : lq ≠ M2 := fun hh => by
        rw [loc_injective hlq (hh ▸ hM2U),
          (hE.unset d2 hmd2).2 (Or.inr (Or.inr (Or.inl rfl)))] at hr
        cases hr
      have hr1 : lr ≠ l1 := fun e => by
        rw [loc_injective hlq ((back _ _ hR1 e) ▸ hl0U), (h.unset c0 hmc0).2 (Or.inl rfl)] at hr
        cases hr
      have hrM1 : lr ≠ M1 := fun e => h2 (loc_injective hlq ((back _ _ hRM1 e) ▸ hM0U))
      have hr2 : lr ≠ l2 := fun e => by
        rw [loc_injective hlq ((back _ _ hR2 e) ▸ hl3U),
          (h.unset c3 hmc3).2 (Or.inr (Or.inr (Or.inr rfl)))] at hr
        cases hr
      have hrM2 : lr ≠ M2 := fun e => h4 (loc_injective hlq ((back _ _ hRM2 e) ▸ hM3U))
      exact ⟨lq, lr, hlq, hlr,
        hg' lq lr hqlt (hR1g lq lr hqlt hqr hq1 hqM1 hr1 hrM1) hq2 hqM2 hr2 hrM2⟩
  · intro q hq
    rw [hA']
    split_ifs with h1 h2 h3 h4
    · subst h1
      exact iff_of_false (by simp) (by
        rintro (e | e | e | e)
        · omega
        · exact hPE _ hmc1 _ hmd1 e
        · exact hPE _ hmc1 _ hmd2 e
        · exact hc31 e.symm)
    · subst h2
      exact iff_of_false (by simp) (by
        rintro (e | e | e | e)
        · exact hPE _ hmc0 _ hmd0 e.symm
        · omega
        · exact hd02 e
        · exact hPE _ hmc3 _ hmd0 e.symm)
    · subst h3
      exact iff_of_false (by simp) (by
        rintro (e | e | e | e)
        · exact hc20 e
        · exact hPE _ hmc2 _ hmd1 e
        · exact hPE _ hmc2 _ hmd2 e
        · omega)
    · subst h4
      exact iff_of_false (by simp) (by
        rintro (e | e | e | e)
        · exact hPE _ hmc0 _ hmd3 e.symm
        · exact hd13 e.symm
        · omega
        · exact hPE _ hmc3 _ hmd3 e.symm)
    · rcases List.mem_append.1 hq with hqP | hqE
      · rw [h.unset q hqP]
        constructor
        · rintro (e | e | e | e)
          · exact Or.inl e
          · exact absurd e h1
          · exact absurd e h3
          · exact Or.inr (Or.inr (Or.inr e))
        · rintro (e | e | e | e)
          · exact Or.inl e
          · exact absurd e (hPE _ hqP _ hmd1)
          · exact absurd e (hPE _ hqP _ hmd2)
          · exact Or.inr (Or.inr (Or.inr e))
      · rw [hE.unset q hqE]
        constructor
        · rintro (e | e | e | e)
          · exact absurd e h2
          · exact Or.inr (Or.inl e)
          · exact Or.inr (Or.inr (Or.inl e))
          · exact absurd e h4
        · rintro (e | e | e | e)
          · exact absurd e.symm (hPE _ hmc0 _ hqE)
          · exact Or.inr (Or.inl e)
          · exact Or.inr (Or.inr (Or.inl e))
          · exact absurd e.symm (hPE _ hmc3 _ hqE)
  · rw [vert_loc hl0U, ← vert_loc hl0]; exact hvu
  · rw [vert_loc hM2U, ← vert_loc hm2]; exact hvv'
  · rw [RotationSystem.sameFaceOrbit_iff hρ'.total hρ'.involution
      (by rw [hρ's, hRs]) (by rw [hρ's]; exact hl0R)]
    exact hfaceU

end Piece

end Spqr
