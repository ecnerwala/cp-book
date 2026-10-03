import Spqr.PlanarEmbedNodeS
import Spqr.PlanarEmbedNodeP

/-!
# `NodeCorners` at `S` and `P` nodes

Read off the fixed layouts `layoutRot .S` / `layoutRot .P` (`rotS`, `rotP`) and the skeleton
shapes: an `S` node's corner at inner node-vertex `v0 + k` is slot `2` of edge `k`, facing slot `1`
of edge `k + 1`; a `P` node has no corner (slot `2` always faces a slot `3`).
-/

namespace Spqr
namespace PlanarSpqrTree

variable (t : PlanarSpqrTree)

theorem get?_of_get! {a : Array (Option Nat)} {i x : Nat} (h : a[i]! = some x) :
    a[i]? = some (some x) := by
  by_cases hi : i < a.size
  · rw [getElem!_pos a i hi] at h; rw [Array.getElem?_eq_getElem hi, h]
  · rw [getElem!_neg a i hi] at h; cases h

theorem cornerAt_decomp {i nv ta : Nat} (h : t.CornerAt i nv ta) :
    ∃ k, 1 ≤ k ∧ (t.toSpqrTree.neRange i).1 + k < (t.toSpqrTree.neRange i).2 ∧
      ta = 4 * ((t.toSpqrTree.neRange i).1 + k) + 2 := by
  obtain ⟨h1, h2, h3, -, -⟩ := h
  exact ⟨ta / 4 - (t.toSpqrTree.neRange i).1, by omega, by omega, by omega⟩

theorem rotS_corner (n s k : Nat) (hk1 : 1 ≤ k) :
    rotS n s k 2 = if k + 1 = n then 4 * s + 3 else 4 * (s + k + 1) + 1 := by
  unfold rotS; split_ifs <;> first | contradiction | omega

theorem rotP_two_mod (k' s k : Nat) : rotP k' s k 2 % 4 = 3 := by
  unfold rotP; dsimp only; split_ifs <;> first | contradiction | omega

/-- The conclusion of `NodeCorners` at one node. -/
def CornersAt (i : Nat) : Prop :=
  (∀ nv ta, t.toSpqrTree.CapEnd i nv → ¬ t.CornerAt i nv ta) ∧
  ∀ nv, t.toSpqrTree.NvOf i nv → ¬ t.toSpqrTree.CapEnd i nv →
    ∃ ta, t.CornerAt i nv ta ∧ (∀ tb, t.neRotAdj[ta]? = some (some tb) → ta < tb) ∧
      ∀ ta', t.CornerAt i nv ta' → ta' = ta

section S

variable {i : Nat} (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hS : t.toSpqrTree.type i = .S)
include hwf hi hS

theorem neEn_S : (t.toSpqrTree.neRange i).2 = (t.toSpqrTree.neRange i).1 + t.toSpqrTree.nVerts i := by
  have := t.nEdges_S hwf hi hS
  have := PlanarRot.neRange_mono t.toSpqrTree hwf i hi
  omega

theorem nvEn_S : (t.toSpqrTree.nvRange i).2 = (t.toSpqrTree.nvRange i).1 + t.toSpqrTree.nVerts i := by
  have h3 := (t.shape_S hwf hi hS).1
  unfold SpqrTree.nVerts at h3 ⊢
  omega

theorem capEnd_S_iff {nv : Nat} : t.toSpqrTree.CapEnd i nv ↔
    nv = (t.toSpqrTree.nvRange i).1 ∨
      nv = (t.toSpqrTree.nvRange i).1 + t.toSpqrTree.nVerts i - 1 := by
  have h3 := (t.shape_S hwf hi hS).1
  have hnvs := t.nvs_S hwf hi hS (k := 0) (by omega)
  rw [ite_eq_left rfl, Nat.add_zero] at hnvs
  have hle := PlanarRot.neRange_le_size t.toSpqrTree hwf i hi
  have hen := t.neEn_S hwf hi hS
  have hsz : (t.toSpqrTree.neRange i).1 < t.nodeEdges.size := by omega
  constructor
  · rintro ⟨ne, d, hne, hd, hor⟩
    rw [t.capNe_S hS] at hne
    cases hne
    rw [getElem!_of_getElem? hd] at hnvs
    rw [hnvs] at hor
    simp only at hor
    omega
  · intro h
    refine ⟨_, _, t.capNe_S hS, Array.getElem?_eq_getElem hsz, ?_⟩
    rw [getElem!_pos t.nodeEdges _ hsz] at hnvs
    rw [hnvs]
    simp only
    omega

theorem cornerAt_S {nv ta : Nat} (h : t.CornerAt i nv ta) :
    ∃ k, 1 ≤ k ∧ k < t.toSpqrTree.nVerts i ∧ ta = 4 * ((t.toSpqrTree.neRange i).1 + k) + 2 ∧
      nv = (t.toSpqrTree.nvRange i).1 + k := by
  obtain ⟨k, hk1, hk2, hta⟩ := t.cornerAt_decomp h
  have hen := t.neEn_S hwf hi hS
  obtain ⟨-, -, -, ⟨d, hd, hnv⟩, -⟩ := h
  subst hta
  have hnvs := t.nvs_S hwf hi hS (k := k) (by omega)
  rw [ite_eq_right (by omega)] at hnvs
  rw [show (4 * ((t.toSpqrTree.neRange i).1 + k) + 2) / 4 = (t.toSpqrTree.neRange i).1 + k by
    omega] at hd
  rw [getElem!_of_getElem? hd] at hnvs
  rw [hnvs] at hnv
  simp only at hnv
  exact ⟨k, hk1, by omega, rfl, by omega⟩

theorem cornersAt_S (hlay : t.LayoutAt i) : t.CornersAt i := by
  have h3 := (t.shape_S hwf hi hS).1
  have hen := t.neEn_S hwf hi hS
  have hnv := t.nvEn_S hwf hi hS
  have hle := PlanarRot.neRange_le_size t.toSpqrTree hwf i hi
  constructor
  · intro nv ta hcap hcor
    obtain ⟨k, hk1, hk, hta, hnvk⟩ := t.cornerAt_S hwf hi hS hcor
    rw [t.capEnd_S_iff hwf hi hS] at hcap
    have hk' : k + 1 = t.toSpqrTree.nVerts i := by omega
    obtain ⟨-, -, -, -, tb, htb, h1⟩ := hcor
    subst hta
    have hr := t.neRotAdj_S hwf hi hS hlay (k := k) (r := 2) hk (by omega)
    rw [get?_of_get! hr] at htb
    simp only [Option.some.injEq] at htb
    rw [rotS_corner _ _ _ hk1, ite_eq_left hk'] at htb
    omega
  · intro nv hnvof hncap
    rw [t.capEnd_S_iff hwf hi hS] at hncap
    obtain ⟨hlo, hhi⟩ := hnvof
    have hk1 : 1 ≤ nv - (t.toSpqrTree.nvRange i).1 := by omega
    have hk2 : nv - (t.toSpqrTree.nvRange i).1 + 1 < t.toSpqrTree.nVerts i := by omega
    have hr := t.neRotAdj_S hwf hi hS hlay (k := nv - (t.toSpqrTree.nvRange i).1) (r := 2)
      (by omega) (by omega)
    rw [rotS_corner _ _ _ hk1, ite_eq_right (by omega)] at hr
    have hnvs := t.nvs_S hwf hi hS (k := nv - (t.toSpqrTree.nvRange i).1) (by omega)
    rw [ite_eq_right (by omega)] at hnvs
    have hlt : (t.toSpqrTree.neRange i).1 + (nv - (t.toSpqrTree.nvRange i).1) < t.nodeEdges.size := by
      omega
    refine ⟨4 * ((t.toSpqrTree.neRange i).1 + (nv - (t.toSpqrTree.nvRange i).1)) + 2, ?_, ?_, ?_⟩
    · refine ⟨by omega, by omega, by omega,
        ⟨t.nodeEdges[(t.toSpqrTree.neRange i).1 + (nv - (t.toSpqrTree.nvRange i).1)]!, ?_, ?_⟩,
        _, get?_of_get! hr, by omega⟩
      · rw [show (4 * ((t.toSpqrTree.neRange i).1 + (nv - (t.toSpqrTree.nvRange i).1)) + 2) / 4 =
          (t.toSpqrTree.neRange i).1 + (nv - (t.toSpqrTree.nvRange i).1) by omega,
          getElem!_pos t.nodeEdges _ hlt]
        exact Array.getElem?_eq_getElem hlt
      · rw [hnvs]
        simp only
        omega
    · intro tb htb
      rw [get?_of_get! hr] at htb
      simp only [Option.some.injEq] at htb
      omega
    · intro ta' h'
      obtain ⟨k', -, -, hta', hnv'⟩ := t.cornerAt_S hwf hi hS h'
      omega

end S

section P

variable {i : Nat} (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hP : t.toSpqrTree.type i = .P)
include hwf hi hP

theorem cornersAt_P (hlay : t.LayoutAt i) : t.CornersAt i := by
  obtain ⟨hn2, hne3, -⟩ := t.shape_P hwf hi hP
  have hle := PlanarRot.neRange_le_size t.toSpqrTree hwf i hi
  have hne3' := hne3
  unfold SpqrTree.nEdges at hne3'
  have hno : ∀ nv ta, ¬ t.CornerAt i nv ta := by
    intro nv ta h
    obtain ⟨k, hk1, hk2, hta⟩ := t.cornerAt_decomp h
    obtain ⟨-, -, -, -, tb, htb, h1⟩ := h
    subst hta
    have hr := t.neRotAdj_P hwf hi hP hlay (k := k) (r := 2) (by unfold SpqrTree.nEdges; omega)
      (by omega)
    rw [get?_of_get! hr] at htb
    simp only [Option.some.injEq] at htb
    have := rotP_two_mod (t.toSpqrTree.nEdges i) (t.toSpqrTree.neRange i).1 k
    omega
  refine ⟨fun nv ta _ => hno nv ta, fun nv hnvof hncap => ?_⟩
  exfalso
  apply hncap
  obtain ⟨hlo, hhi⟩ := hnvof
  have hnv : (t.toSpqrTree.nvRange i).2 = (t.toSpqrTree.nvRange i).1 + 2 := by
    unfold SpqrTree.nVerts at hn2; omega
  have hnvs := t.nvs_P hwf hi hP (k := 0) (by omega)
  rw [Nat.add_zero] at hnvs
  have hsz : (t.toSpqrTree.neRange i).1 < t.nodeEdges.size := by omega
  refine ⟨_, _, t.capNe_P hP, Array.getElem?_eq_getElem hsz, ?_⟩
  rw [getElem!_pos t.nodeEdges _ hsz] at hnvs
  rw [hnvs]
  simp only
  omega

end P

theorem nodeCorners_of_cornersAt
    (h : ∀ i, i < t.size →
      (t.toSpqrTree.type i = .S ∨ t.toSpqrTree.type i = .P ∨ t.toSpqrTree.type i = .R) →
      t.CornersAt i) : t.NodeCorners :=
  ⟨fun i nv ta hi hty hcap => (h i hi hty).1 nv ta hcap,
   fun i nv hi hty hnv hncap => (h i hi hty).2 nv hnv hncap⟩

end PlanarSpqrTree
end Spqr
