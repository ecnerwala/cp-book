import Spqr.Proofs.PlanarSubst
import Spqr.Proofs.PlanarInsert

/-!
# Substituting a capped piece for edge `1`

`subst_step`: in the planar embedding `ρ₁` of `es₁`, edge `1 = (u, v)` is replaced by the capped
piece `(esc, ρc)` (exposed pairs `l0 ↔ l1` at `u`, `l2 ↔ l3` at `v`, `l0`/`l2` cofacial): the cap
`(u, v)` is hung into `ρc` by `RotationSystem.insert` and the two systems are 2-summed by
`IsPlanarEmbedding.subst`. The rotation of the result is given pointwise: the piece's slots sit at
offset `off = 4 * (es₁.length - 1)`, the other quarter-edges of `es₁` shift down by `4` past edge
`1` (`sh1`), and a corner `q ↔ 4 + s` of `es₁` becomes `q ↔ off + l_s`.
-/

namespace Spqr

open RotationSystem

/-- Index of a quarter-edge of `es₁` after deleting edge `1`. -/
def sh1 (s : Nat) : Nat := if s < 4 then s else s - 4

/-- The exposed slot of a capped piece facing the quarter `k % 4` of its cap. -/
def pick4 (l0 l1 l2 l3 k : Nat) : Nat :=
  if k % 4 = 0 then l0 else if k % 4 = 1 then l1 else if k % 4 = 2 then l2 else l3

theorem pick4_mod (l0 l1 l2 l3 k : Nat) : pick4 l0 l1 l2 l3 (k % 4) = pick4 l0 l1 l2 l3 k := by
  simp [pick4]

theorem subst_step {es₁ esc : List (Nat × Nat)} {n u v : Nat} {ρ₁ ρc : RotationSystem}
    (h₁ : IsPlanarEmbedding es₁ n ρ₁) (he : es₁[1]? = some (u, v)) (huv : u ≠ v)
    (hdeg : ∀ k, k < 4 → ρ₁.rot (4 + k) / 4 ≠ 1)
    (hc : IsPlanarEmbedding esc n ρc) {l0 l1 l2 l3 : Nat}
    (hl0 : ρc.get l0 = some l1) (hl2 : ρc.get l2 = some l3)
    (h0 : l0 % 2 = 0) (h2 : l2 % 2 = 0) (hl0lt : l0 < ρc.size) (hl2lt : l2 < ρc.size)
    (hvu : QE.vert esc l0 = some u) (hvv : QE.vert esc l2 = some v)
    (hface : SameOrbit (ρc.stepC 3) l0 l2)
    (hsep : ∀ w, HasEdge (es₁.eraseIdx 1) w → HasEdge esc w → w = u ∨ w = v) :
    IsPlanarEmbedding (es₁.eraseIdx 1 ++ esc) n
      ((substSum es₁ ((u, v) :: esc) n 1 u v).splice ρ₁ (ρc.insert 0 l0 l2)) ∧
    ((substSum es₁ ((u, v) :: esc) n 1 u v).splice ρ₁ (ρc.insert 0 l0 l2)).size =
      4 * (es₁.length - 1) + ρc.size ∧
    (∀ r, r < 4 * (es₁.length - 1) →
      ((substSum es₁ ((u, v) :: esc) n 1 u v).splice ρ₁ (ρc.insert 0 l0 l2)).rot r =
        if ρ₁.rot (if r < 4 then r else r + 4) / 4 = 1 then
          4 * (es₁.length - 1) + pick4 l0 l1 l2 l3 (ρ₁.rot (if r < 4 then r else r + 4))
        else sh1 (ρ₁.rot (if r < 4 then r else r + 4))) ∧
    (∀ l, l < ρc.size →
      ((substSum es₁ ((u, v) :: esc) n 1 u v).splice ρ₁ (ρc.insert 0 l0 l2)).rot
          (4 * (es₁.length - 1) + l) =
        if l = l0 then sh1 (ρ₁.rot 4) else if l = l1 then sh1 (ρ₁.rot 5)
        else if l = l2 then sh1 (ρ₁.rot 6) else if l = l3 then sh1 (ρ₁.rot 7)
        else 4 * (es₁.length - 1) + ρc.rot l) := by
  have hm₁ : 1 < es₁.length := (List.getElem?_eq_some_iff.1 he).1
  have hl02 : l0 ≠ l2 := fun h => huv (by rw [h, hvv] at hvu; exact (Option.some.inj hvu).symm)
  have H : InsertSetting ρc 0 l0 l2 :=
    ⟨hc.total, hc.involution, hc.opposite_dir, by rw [hc.size]; omega, Or.inl rfl, hl0lt, hl2lt,
      h0, h2, hl02⟩
  have hrot0 : ρc.rot l0 = l1 := rot_eq_of_get hl0
  have hrot2 : ρc.rot l2 = l3 := rot_eq_of_get hl2
  have hl1lt : l1 < ρc.size := hrot0 ▸ H.a1_lt
  have hl3lt : l3 < ρc.size := hrot2 ▸ H.a3_lt
  have hins : IsPlanarEmbedding ((u, v) :: esc) n (ρc.insert 0 l0 l2) :=
    hc.insert (Or.inl rfl) hl0lt hl2lt h0 h2 hl02 hvu hvv huv hface (p := (u, v))
      (by rw [vert_cons_lt _ _ (by decide)]; rfl) (by rw [vert_cons_lt _ _ (by decide)]; rfl)
  have hconn : EdgesConn esc u v :=
    edgesConn_of_sameOrbit hc.total hc.involution hc.same_vertex hc.size hface hvu hvv
  have hface' : ¬(ρc.insert 0 l0 l2).SameFaceOrbit 0 2 := H.insert_not_sameFaceOrbit hface
  have hdeg' : ∀ k, k < 4 → ∀ s ∈ ρ₁.get (4 * 1 + k), s / 4 ≠ 1 := by
    intro k hk s hs
    rw [← rot_eq_of_get hs, show 4 * 1 + k = 4 + k by omega]
    exact hdeg k hk
  have he₂ : ((u, v) :: esc)[0]? = some (u, v) := rfl
  have hconn' : EdgesConn (((u, v) :: esc).eraseIdx 0) u v := by simpa using hconn
  have hsep' : ∀ w, HasEdge (es₁.eraseIdx 1) w → HasEdge (((u, v) :: esc).eraseIdx 0) w →
      w = u ∨ w = v := by simpa using hsep
  have hpl := IsPlanarEmbedding.subst h₁ hins he he₂ huv hdeg' hface' hconn' hsep'
  have W := substSum_wf h₁ hins he he₂ huv hdeg' hface' hconn'
  simp only [List.eraseIdx_cons_zero] at hpl
  set T := substSum es₁ ((u, v) :: esc) n 1 u v with hT
  have hsz : (T.splice ρ₁ (ρc.insert 0 l0 l2)).size = 4 * (es₁.length - 1) + ρc.size := by
    rw [substSum_size, hc.size]
    simp only [List.length_cons]
    omega
  have hoff : T.off = 4 * (es₁.length - 1) := rfl
  have hm₁' : T.m₁ = es₁.length := rfl
  have hm₂' : T.m₂ = esc.length + 1 := rfl
  have he₁' : T.e₁ = 1 := rfl
  have he₂' : T.e₂ = 0 := rfl
  have hu₁ : T.u₁ = u := rfl
  have hv₁ : T.v₁ = v := rfl
  have hcs : ρc.size = 4 * esc.length := hc.size
  -- the inserted cap's rotation
  have hI : ∀ k, k < 4 → (ρc.insert 0 l0 l2).rot (k ^^^ 1) = 4 + pick4 l0 l1 l2 l3 k := by
    intro k hk
    have e0 := H.insert_get_x
    have e1 := H.insert_get_x1
    have e2 := H.insert_get_x2
    have e3 := H.insert_get_x3
    simp only [Nat.zero_add, Nat.sub_zero, hrot0, hrot2] at e0 e1 e2 e3
    unfold pick4
    rcases (by omega : k = 0 ∨ k = 1 ∨ k = 2 ∨ k = 3) with rfl | rfl | rfl | rfl
    · simpa using rot_eq_of_get e1
    · simpa using rot_eq_of_get e0
    · simpa using rot_eq_of_get e3
    · simpa using rot_eq_of_get e2
  have hA0 := H.insert_get_a0
  have hA1 := H.insert_get_a1
  have hA2 := H.insert_get_a2
  have hA3 := H.insert_get_a3
  simp only [Nat.zero_add, Nat.sub_zero, hrot0, hrot2] at hA0 hA1 hA2 hA3
  refine ⟨hpl, hsz, ?_, ?_⟩
  · intro r hr
    have hr' : r < 4 * (T.m₁ + T.m₂ - 2) := by rw [hm₁', hm₂']; omega
    rw [substSum_rot W hr']
    have hpre : T.pre₁ r = if r < 4 then r else r + 4 := by
      unfold TwoSum.pre₁ insIdx; rw [he₁']; split_ifs <;> omega
    unfold TwoSum.glue
    rw [ite_eq_left (by rw [hoff]; exact hr), hpre, he₁']
    set s := ρ₁.rot (if r < 4 then r else r + 4) with hs_def
    have hpk : pick4 l0 l1 l2 l3 s < ρc.size := by unfold pick4; split_ifs <;> assumption
    by_cases hs : s / 4 = 1
    · simp only [hs, ↓reduceIte]
      unfold TwoSum.ρ₂
      rw [he₂', Nat.mul_zero, Nat.zero_add, hI _ (Nat.mod_lt _ (by omega)), pick4_mod]
      unfold TwoSum.qe₂ delIdx
      rw [he₂', hoff]
      split_ifs <;> omega
    · simp only [hs, ↓reduceIte]
      unfold TwoSum.qe₁ delIdx sh1
      rw [he₁']
      split_ifs <;> omega
  · intro l hl
    have hr' : 4 * (es₁.length - 1) + l < 4 * (T.m₁ + T.m₂ - 2) := by rw [hm₁', hm₂']; omega
    rw [substSum_rot W hr']
    have hpre : T.pre₂ (4 * (es₁.length - 1) + l) = 4 + l := by
      unfold TwoSum.pre₂ insIdx; rw [he₂', hoff]; split_ifs <;> omega
    unfold TwoSum.glue
    rw [ite_eq_right (by rw [hoff]; omega), hpre, he₂']
    have hqe₁ : ∀ s, s / 4 ≠ 1 → T.qe₁ s = sh1 s := by
      intro s hs; unfold TwoSum.qe₁ delIdx sh1; rw [he₁']; split_ifs <;> omega
    have hρ₁ : ∀ k, T.ρ₁ ρ₁ k = ρ₁.rot (4 + k) := by
      intro k; unfold TwoSum.ρ₁; rw [he₁']
    by_cases hl0' : l = l0
    · subst hl0'
      rw [rot_eq_of_get hA0, ite_eq_left (by decide), hρ₁, show (1 % 4 ^^^ 1 : Nat) = 0 by decide,
        hqe₁ _ (hdeg 0 (by decide)), ite_eq_left rfl]
    rw [ite_eq_right hl0']
    by_cases hl1' : l = l1
    · subst hl1'
      rw [rot_eq_of_get hA1, ite_eq_left (by decide), hρ₁, show (0 % 4 ^^^ 1 : Nat) = 1 by decide,
        hqe₁ _ (hdeg 1 (by decide)), ite_eq_left rfl]
    rw [ite_eq_right hl1']
    by_cases hl2' : l = l2
    · subst hl2'
      rw [rot_eq_of_get hA2, ite_eq_left (by decide), hρ₁, show (3 % 4 ^^^ 1 : Nat) = 2 by decide,
        hqe₁ _ (hdeg 2 (by decide)), ite_eq_left rfl]
    rw [ite_eq_right hl2']
    by_cases hl3' : l = l3
    · subst hl3'
      rw [rot_eq_of_get hA3, ite_eq_left (by decide), hρ₁, show (2 % 4 ^^^ 1 : Nat) = 3 by decide,
        hqe₁ _ (hdeg 3 (by decide)), ite_eq_left rfl]
    rw [ite_eq_right hl3']
    rw [rot_eq_of_get (H.insert_get_add hl hl0' (hrot0 ▸ hl1') hl2' (hrot2 ▸ hl3'))]
    have hlt := rot_lt hc.total hc.involution hl
    rw [ite_eq_right (by omega)]
    unfold TwoSum.qe₂ delIdx
    rw [he₂', hoff]
    split_ifs <;> omega

end Spqr
