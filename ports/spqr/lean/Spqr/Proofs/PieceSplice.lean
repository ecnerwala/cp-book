import Spqr.PlanarEmbedClose
import Spqr.Proofs.PlanarSplice

namespace Spqr.Piece

structure OpenEmbedding (P : Piece) (A : Array (Option Nat)) (ρ : RotationSystem)
    (a b v : Nat) : Prop where
  planar : IsPlanarEmbedding P.es P.nVerts ρ
  agrees : P.Agrees A ρ
  boundary : ∃ la lb, P.loc a = some la ∧ P.loc b = some lb ∧
    ρ.get la = some lb ∧ QE.vert P.es la = some v
  unset : ∀ q, P.Mem q → (A[q]? = some none ↔ q = a ∨ q = b)
  left_dir : a % 2 = 0
  right_dir : b % 2 = 1

theorem OpenEmbedding.frame {P : Piece} {A B : Array (Option Nat)} {ρ : RotationSystem}
    {a b v : Nat} (h : P.OpenEmbedding A ρ a b v)
    (hframe : ∀ q, P.Mem q → B[q]? = A[q]?) : P.OpenEmbedding B ρ a b v where
  planar := h.planar
  agrees := fun q r hq hr => h.agrees q r hq (by rwa [hframe q hq] at hr)
  boundary := h.boundary
  unset := by intro q hq; rw [hframe q hq]; exact h.unset q hq
  left_dir := h.left_dir
  right_dir := h.right_dir

theorem OpenEmbedding.splice {P : Piece} {E : List Nat} {ρ₁ ρ₂ : RotationSystem}
    {s : PlanarSpqrTree.EmbedState} {a b c d v : Nat}
    (h₁ : P.OpenEmbedding s.rotAdj ρ₁ a b v)
    (h₂ : ({P with ves := E} : Piece).OpenEmbedding s.rotAdj ρ₂ c d v)
    (hdis : List.Disjoint P.ves E)
    (hsep : ∀ w, HasEdge P.es w → HasEdge ({P with ves := E} : Piece).es w → w = v)
    (hbound : ∀ q, ({P with ves := P.ves ++ E} : Piece).Mem q → q < s.rotAdj.size) :
    ∃ ρ, ({P with ves := P.ves ++ E} : Piece).OpenEmbedding
      ((PlanarSpqrTree.EmbedM.link (some b) (some c)).run s).2.rotAdj ρ a d v := by
  let U : Piece := {P with ves := P.ves ++ E}
  obtain ⟨la, lb, hla, hlb, hab, hva⟩ := h₁.boundary
  obtain ⟨lc, ld, hlc, hld, hcd, hvc⟩ := h₂.boundary
  have hs₁ : ρ₁.size = 4 * P.ves.length := by simpa [es] using h₁.planar.size
  have hs₂ : ρ₂.size = 4 * E.length := by simpa [es] using h₂.planar.size
  have ha : la < ρ₁.size := by rw [hs₁]; exact loc_lt hla
  have hb : lb < ρ₁.size := by rw [hs₁]; exact loc_lt hlb
  have hc : lc < ρ₂.size := by rw [hs₂]; exact loc_lt hlc
  have hd : ld < ρ₂.size := by rw [hs₂]; exact loc_lt hld
  have hba := (h₁.planar.involution la ha lb hab).2
  have hdc := (h₂.planar.involution lc hc ld hcd).2
  have hvb : QE.vert P.es lb = some v := (h₁.planar.same_vertex la ha lb hab).symm.trans hva
  have hvd : QE.vert ({P with ves := E} : Piece).es ld = some v :=
    (h₂.planar.same_vertex lc hc ld hcd).symm.trans hvc
  have hpar : lb % 2 = ld % 2 := by
    rw [loc_mod_two hlb, loc_mod_two hld, h₁.right_dir, h₂.right_dir]
  have hplan := h₁.planar.splice h₂.planar (by simpa [es] using loc_lt hlb)
    (by simpa [es] using loc_lt hld) hpar hvb hvd hsep
  let R := ρ₁.union ρ₂
  let C := ρ₁.size + lc
  let D := ρ₁.size + ld
  let ρ := R.conj lb D
  have hρ : IsPlanarEmbedding U.es U.nVerts ρ := by
    simpa only [U, es, List.map_append] using hplan
  have hlaU : U.loc a = some la := loc_append_left E hla
  have hlbU : U.loc b = some lb := loc_append_left E hlb
  have hnot : ∀ q, ({P with ves := E} : Piece).Mem q → QE.edge q ∉ P.ves := by
    intro q hq hh
    exact List.disjoint_left.1 hdis hh hq
  have hlcU : U.loc c = some C := by
    simpa only [C, hs₁] using loc_append_right P.ves (hnot c (mem_of_loc hlc)) hlc
  have hldU : U.loc d = some D := by
    simpa only [D, hs₁] using loc_append_right P.ves (hnot d (mem_of_loc hld)) hld
  have hRa : R.get la = some lb := by rw [RotationSystem.union_get_lt _ _ ha]; exact hab
  have hRb : R.get lb = some la := by rw [RotationSystem.union_get_lt _ _ hb]; exact hba
  have hRc : R.get C = some D := by
    rw [RotationSystem.union_get_ge _ _ (by dsimp [C]; omega)]
    simp only [C, D, Nat.add_sub_cancel_left, hcd, Option.map_some]
  have hRd : R.get D = some C := by
    rw [RotationSystem.union_get_ge _ _ (by dsimp [D]; omega)]
    simp only [C, D, Nat.add_sub_cancel_left, hdc, Option.map_some]
  have hUsize : R.size = ρ₁.size + ρ₂.size := RotationSystem.union_size ..
  have hUloc : ∀ q l, U.loc q = some l → l < R.size := by
    intro q l hl
    rw [hUsize, hs₁, hs₂]
    have := loc_lt hl
    simpa only [U, List.length_append, Nat.mul_add] using this
  have hdirs : la ≠ lb ∧ C ≠ D := by
    have habdir := loc_mod_two hla
    have hbbdir := loc_mod_two hlb
    have hccdir := loc_mod_two hlc
    have hdddir := loc_mod_two hld
    have := h₁.left_dir; have := h₁.right_dir; have := h₂.left_dir; have := h₂.right_dir
    dsimp [C, D]
    omega
  have had : la ≠ D := by dsimp [D]; omega
  have hbc : lb ≠ C := by dsimp [C]; omega
  have hbd : lb ≠ D := by dsimp [D]; omega
  have hpair : ρ.get la = some D := by
    rw [RotationSystem.conj_get_lt _ _ _ (hUloc a la hlaU)]
    simp only [Equiv.swap_apply_of_ne_of_ne hdirs.1 had, hRa, Option.map_some,
      Equiv.swap_apply_left]
  have hlinkb : ρ.get lb = some C := by
    rw [RotationSystem.conj_get_lt _ _ _ (hUloc b lb hlbU)]
    simp only [Equiv.swap_apply_left, hRd, Option.map_some,
      Equiv.swap_apply_of_ne_of_ne hbc.symm hdirs.2]
  have hlinkc : ρ.get C = some lb := by
    rw [RotationSystem.conj_get_lt _ _ _ (hUloc c C hlcU)]
    simp only [Equiv.swap_apply_of_ne_of_ne hbc.symm hdirs.2, hRc, Option.map_some,
      Equiv.swap_apply_right]
  have hneq : b ≠ c := by
    intro heq
    have : lb = C := Option.some.inj (hlbU.symm.trans (heq ▸ hlcU))
    exact hbc this
  have hget := PlanarSpqrTree.link_rotAdj_get b c s (hbound b (mem_of_loc hlbU))
    (hbound c (mem_of_loc hlcU)) hneq
  have hAg : U.Agrees s.rotAdj R := agrees_append hdis hs₁ h₁.agrees h₂.agrees
  have hnone : ∀ q, U.Mem q → (s.rotAdj[q]? = some none ↔ q = a ∨ q = b ∨ q = c ∨ q = d) := by
    intro q hq
    rcases List.mem_append.1 hq with hq | hq
    · rw [h₁.unset q hq]
      have hqc : q ≠ c := fun hh => hnot c (mem_of_loc hlc) (hh ▸ hq)
      have hqd : q ≠ d := fun hh => hnot d (mem_of_loc hld) (hh ▸ hq)
      simp [hqc, hqd]
    · rw [h₂.unset q hq]
      have hqa : q ≠ a := fun hh => hnot a (hh ▸ hq) (mem_of_loc hla)
      have hqb : q ≠ b := fun hh => hnot b (hh ▸ hq) (mem_of_loc hlb)
      simp [hqa, hqb]
  refine ⟨ρ, hρ, ?_, ⟨la, D, hlaU, hldU, hpair, ?_⟩, ?_, h₁.left_dir, h₂.right_dir⟩
  · intro q r hq hr
    rw [hget] at hr
    split_ifs at hr with hqb hqc
    · subst q
      have : r = c := (Option.some.inj (Option.some.inj hr)).symm
      subst r
      exact ⟨lb, C, hlbU, hlcU, hlinkb⟩
    · subst q
      have : r = b := (Option.some.inj (Option.some.inj hr)).symm
      subst r
      exact ⟨C, lb, hlcU, hlbU, hlinkc⟩
    · obtain ⟨lq, lr, hlq, hlr, hqr⟩ := hAg q r hq hr
      have hqx : q ≠ a ∧ q ≠ d := by
        constructor <;> intro hh
        all_goals have hn := (hnone q hq).2 (by simp [hh])
        all_goals rw [hr] at hn; cases hn
      have hqB : lq ≠ lb := fun hh => hqb (loc_injective hlq (hh ▸ hlbU))
      have hqD : lq ≠ D := fun hh => hqx.2 (loc_injective hlq (hh ▸ hldU))
      have hRI : R.Involution := RotationSystem.union_involution ρ₁ ρ₂
        h₁.planar.involution h₂.planar.involution
      have hrB : lr ≠ lb := by
        intro hh
        subst lr
        have hbq := (hRI lq (hUloc q lq hlq) lb hqr).2
        have : lq = la := Option.some.inj (hbq.symm.trans hRb)
        exact hqx.1 (loc_injective hlq (this ▸ hlaU))
      have hrD : lr ≠ D := by
        intro hh
        subst lr
        have hdq := (hRI lq (hUloc q lq hlq) D hqr).2
        have : lq = C := Option.some.inj (hdq.symm.trans hRd)
        exact hqc (loc_injective hlq (this ▸ hlcU))
      refine ⟨lq, lr, hlq, hlr, ?_⟩
      rw [RotationSystem.conj_get_lt _ _ _ (hUloc q lq hlq)]
      simp only [Equiv.swap_apply_of_ne_of_ne hqB hqD, hqr, Option.map_some,
        Equiv.swap_apply_of_ne_of_ne hrB hrD]
  · rw [vert_loc hlaU, ← vert_loc hla]
    exact hva
  · intro q hq
    rw [hget q]
    have hba' : b ≠ a := by have := h₁.left_dir; have := h₁.right_dir; omega
    have hcd' : c ≠ d := by have := h₂.left_dir; have := h₂.right_dir; omega
    have hbd' : b ≠ d := fun hh => hnot d (mem_of_loc hld) (hh ▸ mem_of_loc hlb)
    have hca' : c ≠ a := fun hh => hnot a (hh ▸ mem_of_loc hlc) (mem_of_loc hla)
    split_ifs with hqb hqc
    · simp [hqb, hba', hbd']
    · simp [hqc, hca', hcd']
    · rw [hnone q hq]
      simp [hqb, hqc]

end Spqr.Piece
