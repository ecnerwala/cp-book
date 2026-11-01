import Spqr.PlanarGlue
import Mathlib.Data.Finset.Card
import Mathlib.Data.Finset.Image
import Mathlib.Logic.Function.Iterate
import Mathlib.Logic.Equiv.Basic

/-!
# Orbit counting for permutations of a finite set

`orbitCount f S` is the number of orbits of `f` on `S`.  We prove it is invariant under
relabelling, additive over disjoint unions, preserved by skipping a non-fixed point, and
drops by one when the images of two points on different orbits are exchanged.
-/

namespace Spqr

open Classical

theorem SameOrbit.refl (f : Nat → Nat) (q : Nat) : SameOrbit f q q := Relation.EqvGen.refl q
theorem SameOrbit.symm {f : Nat → Nat} {q r : Nat} (h : SameOrbit f q r) : SameOrbit f r q :=
  Relation.EqvGen.symm _ _ h
theorem SameOrbit.trans {f : Nat → Nat} {q r s : Nat} (h : SameOrbit f q r)
    (h' : SameOrbit f r s) : SameOrbit f q s :=
  Relation.EqvGen.trans _ _ _ h h'
theorem SameOrbit.step {f : Nat → Nat} (q : Nat) : SameOrbit f q (f q) :=
  Relation.EqvGen.rel _ _ rfl

theorem mem_orbit {f : Nat → Nat} {S : Finset Nat} {q r : Nat} :
    r ∈ orbit f S q ↔ r ∈ S ∧ SameOrbit f q r :=
  Finset.mem_filter

theorem orbit_subset {f : Nat → Nat} {S : Finset Nat} {q : Nat} : orbit f S q ⊆ S :=
  Finset.filter_subset _ _

theorem mem_orbit_self {f : Nat → Nat} {S : Finset Nat} {q : Nat} (hq : q ∈ S) :
    q ∈ orbit f S q :=
  mem_orbit.2 ⟨hq, SameOrbit.refl f q⟩

theorem orbit_eq_of_sameOrbit {f : Nat → Nat} {S : Finset Nat} {q r : Nat}
    (h : SameOrbit f q r) : orbit f S q = orbit f S r := by
  ext s
  simp only [mem_orbit]
  exact ⟨fun ⟨hs, h'⟩ => ⟨hs, h.symm.trans h'⟩, fun ⟨hs, h'⟩ => ⟨hs, h.trans h'⟩⟩

theorem sameOrbit_of_orbit_eq {f : Nat → Nat} {S : Finset Nat} {q r : Nat} (hq : q ∈ S)
    (h : orbit f S q = orbit f S r) : SameOrbit f q r := by
  have := mem_orbit_self (f := f) (S := S) hq
  rw [h, mem_orbit] at this
  exact this.2.symm

theorem orbit_eq_or_disjoint {f : Nat → Nat} {S : Finset Nat} (q r : Nat) :
    orbit f S q = orbit f S r ∨ Disjoint (orbit f S q) (orbit f S r) := by
  by_cases h : ∃ s, s ∈ orbit f S q ∧ s ∈ orbit f S r
  · obtain ⟨s, hs, hs'⟩ := h
    left
    exact orbit_eq_of_sameOrbit ((mem_orbit.1 hs).2.trans (mem_orbit.1 hs').2.symm)
  · right
    rw [Finset.disjoint_left]
    intro s hs hs'
    exact h ⟨s, hs, hs'⟩

namespace IsPermOn

variable {f : Nat → Nat} {S : Finset Nat}

theorem mem_of_apply_mem (hf : IsPermOn f S) {q : Nat} (h : f q ∈ S) : q ∈ S := by
  by_contra hq
  rw [hf.fix q hq] at h
  exact hq h

theorem mem_iff_of_sameOrbit (hf : IsPermOn f S) {q r : Nat} (h : SameOrbit f q r) :
    q ∈ S ↔ r ∈ S := by
  induction h with
  | rel a b h => subst h; exact ⟨hf.maps a, hf.mem_of_apply_mem⟩
  | refl a => exact Iff.rfl
  | symm a b _ ih => exact ih.symm
  | trans a b c _ _ ih₁ ih₂ => exact ih₁.trans ih₂

theorem iterate_mem (hf : IsPermOn f S) (k : Nat) {q : Nat} (hq : q ∈ S) : f^[k] q ∈ S := by
  induction k generalizing q with
  | zero => exact hq
  | succ k ih => rw [Function.iterate_succ_apply]; exact ih (hf.maps q hq)

theorem iterate_inj (hf : IsPermOn f S) (k : Nat) {q r : Nat} (hq : q ∈ S) (hr : r ∈ S)
    (h : f^[k] q = f^[k] r) : q = r := by
  induction k generalizing q r with
  | zero => exact h
  | succ k ih =>
    rw [Function.iterate_succ_apply, Function.iterate_succ_apply] at h
    exact hf.inj q hq r hr (ih (hf.maps q hq) (hf.maps r hr) h)

theorem exists_iterate_eq (hf : IsPermOn f S) {q : Nat} (hq : q ∈ S) :
    ∃ k, 0 < k ∧ f^[k] q = q := by
  have hc : S.card < (Finset.range (S.card + 1)).card := by simp
  obtain ⟨i, -, j, -, hij, h⟩ := Finset.exists_ne_map_eq_of_card_lt_of_maps_to hc
    (f := fun i => f^[i] q) (fun i _ => hf.iterate_mem i hq)
  rcases Nat.lt_or_gt_of_ne hij with hlt | hlt
  · refine ⟨j - i, by omega, hf.iterate_inj i (hf.iterate_mem _ hq) hq ?_⟩
    rw [← Function.iterate_add_apply, show i + (j - i) = j by omega]
    exact h.symm
  · refine ⟨i - j, by omega, hf.iterate_inj j (hf.iterate_mem _ hq) hq ?_⟩
    rw [← Function.iterate_add_apply, show j + (i - j) = i by omega]
    exact h

theorem reach_iff {q r : Nat} : Reach f q r ↔ ∃ k, f^[k] q = r := by
  constructor
  · intro h
    induction h with
    | refl => exact ⟨0, rfl⟩
    | tail _ h ih =>
      obtain ⟨k, rfl⟩ := ih
      exact ⟨k + 1, by rw [Function.iterate_succ_apply']; exact h⟩
  · rintro ⟨k, rfl⟩
    induction k with
    | zero => exact Relation.ReflTransGen.refl
    | succ k ih =>
      rw [Function.iterate_succ_apply']
      exact ih.tail rfl

theorem reach_symm (hf : IsPermOn f S) {q r : Nat} (hq : q ∈ S) (h : Reach f q r) :
    Reach f r q := by
  obtain ⟨j, rfl⟩ := reach_iff.1 h
  obtain ⟨k, hk, hfix⟩ := hf.exists_iterate_eq hq
  refine reach_iff.2 ⟨k * (j + 1) - j, ?_⟩
  have hle : j + 1 ≤ k * (j + 1) := Nat.le_mul_of_pos_left _ hk
  rw [← Function.iterate_add_apply, show k * (j + 1) - j + j = k * (j + 1) by omega,
    Function.iterate_mul]
  exact Function.iterate_fixed hfix _

theorem sameOrbit_of_reach {q r : Nat} (h : Reach f q r) : SameOrbit f q r := by
  induction h with
  | refl => exact SameOrbit.refl f q
  | tail _ h ih => exact ih.trans (h ▸ SameOrbit.step _)

theorem reach_of_sameOrbit (hf : IsPermOn f S) {q r : Nat} (hq : q ∈ S) (h : SameOrbit f q r) :
    Reach f q r := by
  induction h with
  | rel a b h => exact Relation.ReflTransGen.single h
  | refl a => exact Relation.ReflTransGen.refl
  | symm a b h ih => exact hf.reach_symm ((hf.mem_iff_of_sameOrbit h).2 hq) (ih ((hf.mem_iff_of_sameOrbit h).2 hq))
  | trans a b c h₁ _ ih₁ ih₂ => exact (ih₁ hq).trans (ih₂ ((hf.mem_iff_of_sameOrbit h₁).1 hq))

end IsPermOn

/-! ### Relabelling -/

theorem orbitCount_congr {f g : Nat → Nat} {S T : Finset Nat} {φ : Nat → Nat}
    (hf : IsPermOn f S) (hg : IsPermOn g T) (hinj : Set.InjOn φ S) (himg : S.image φ = T)
    (hcomm : ∀ a ∈ S, g (φ a) = φ (f a)) : orbitCount g T = orbitCount f S := by
  have memT : ∀ c ∈ T, ∃ a ∈ S, φ a = c := fun c hc => by
    rw [← himg] at hc; exact Finset.mem_image.1 hc
  have fwd : ∀ a b, SameOrbit f a b → a ∈ S → SameOrbit g (φ a) (φ b) := by
    intro a b h
    induction h with
    | rel a b h => intro ha; subst h; exact (hcomm a ha) ▸ SameOrbit.step _
    | refl a => intro _; exact SameOrbit.refl _ _
    | symm a b h ih => intro hb; exact (ih ((hf.mem_iff_of_sameOrbit h).2 hb)).symm
    | trans a b c h₁ _ ih₁ ih₂ =>
      intro ha; exact (ih₁ ha).trans (ih₂ ((hf.mem_iff_of_sameOrbit h₁).1 ha))
  have bwd : ∀ c d, SameOrbit g c d → ∀ a ∈ S, ∀ b ∈ S, φ a = c → φ b = d → SameOrbit f a b := by
    intro c d h
    induction h with
    | rel c d h =>
      intro a ha b hb hac hbd
      subst hac; subst hbd
      rw [hcomm a ha] at h
      exact hinj (hf.maps a ha) hb h ▸ SameOrbit.step _
    | refl c => intro a ha b hb hac hbc; exact (hinj ha hb (hac.trans hbc.symm)) ▸ SameOrbit.refl _ _
    | symm c d _ ih => intro a ha b hb hac hbd; exact (ih b hb a ha hbd hac).symm
    | trans c d e h₁ _ ih₁ ih₂ =>
      intro a ha b hb hac hbe
      have hdT : d ∈ T := (hg.mem_iff_of_sameOrbit h₁).1 (hac ▸ himg ▸ Finset.mem_image_of_mem φ ha)
      obtain ⟨a', ha', hd⟩ := memT d hdT
      exact (ih₁ a ha a' ha' hac hd).trans (ih₂ a' ha' b hb hd hbe)
  have horb : ∀ a ∈ S, orbit g T (φ a) = (orbit f S a).image φ := by
    intro a ha
    ext r
    simp only [mem_orbit, Finset.mem_image]
    constructor
    · rintro ⟨hr, h⟩
      obtain ⟨b, hb, rfl⟩ := memT r hr
      exact ⟨b, ⟨hb, bwd _ _ h a ha b hb rfl rfl⟩, rfl⟩
    · rintro ⟨b, ⟨hb, h⟩, rfl⟩
      exact ⟨himg ▸ Finset.mem_image_of_mem φ hb, fwd a b h ha⟩
  have hinjC : Set.InjOn (Finset.image φ) ↑(S.image (orbit f S)) := by
    intro A hA B hB hAB
    simp only [Finset.coe_image, Set.mem_image, Finset.mem_coe] at hA hB
    obtain ⟨a, ha, rfl⟩ := hA
    obtain ⟨b, hb, rfl⟩ := hB
    ext s
    constructor
    · intro hs
      have := Finset.mem_image_of_mem φ hs
      rw [hAB, Finset.mem_image] at this
      obtain ⟨s', hs', h⟩ := this
      rwa [hinj (orbit_subset hs') (orbit_subset hs) h] at hs'
    · intro hs
      have := Finset.mem_image_of_mem φ hs
      rw [← hAB, Finset.mem_image] at this
      obtain ⟨s', hs', h⟩ := this
      rwa [hinj (orbit_subset hs') (orbit_subset hs) h] at hs'
  calc orbitCount g T = ((S.image φ).image (orbit g T)).card := by rw [himg]; rfl
    _ = (S.image (fun a => (orbit f S a).image φ)).card := by
        rw [Finset.image_image]
        exact congrArg _ (Finset.image_congr fun a ha => horb a ha)
    _ = ((S.image (orbit f S)).image (Finset.image φ)).card := by rw [Finset.image_image]; rfl
    _ = orbitCount f S := Finset.card_image_of_injOn hinjC

/-! ### Disjoint unions -/

theorem isPermOn_unionStep {f g : Nat → Nat} {S T : Finset Nat} (hf : IsPermOn f S)
    (hg : IsPermOn g T) (hST : Disjoint S T) : IsPermOn (unionStep f g S) (S ∪ T) where
  maps q hq := by
    unfold unionStep
    rcases Finset.mem_union.1 hq with h | h
    · rw [if_pos h]; exact Finset.mem_union_left _ (hf.maps q h)
    · rw [if_neg (Finset.disjoint_right.1 hST h)]; exact Finset.mem_union_right _ (hg.maps q h)
  inj q hq r hr h := by
    unfold unionStep at h
    rcases Finset.mem_union.1 hq with hq | hq <;> rcases Finset.mem_union.1 hr with hr | hr
    · rw [if_pos hq, if_pos hr] at h; exact hf.inj q hq r hr h
    · rw [if_pos hq, if_neg (Finset.disjoint_right.1 hST hr)] at h
      exact absurd (h ▸ hf.maps q hq) (Finset.disjoint_right.1 hST (hg.maps r hr))
    · rw [if_neg (Finset.disjoint_right.1 hST hq), if_pos hr] at h
      exact absurd (h ▸ hg.maps q hq) (Finset.disjoint_left.1 hST (hf.maps r hr))
    · rw [if_neg (Finset.disjoint_right.1 hST hq), if_neg (Finset.disjoint_right.1 hST hr)] at h
      exact hg.inj q hq r hr h
  fix q hq := by
    unfold unionStep
    rw [if_neg (fun h => hq (Finset.mem_union_left _ h))]
    exact hg.fix q (fun h => hq (Finset.mem_union_right _ h))

theorem sameOrbit_unionStep_left {f g : Nat → Nat} {S T : Finset Nat} (hf : IsPermOn f S)
    (hg : IsPermOn g T) (hST : Disjoint S T) {a b : Nat} (ha : a ∈ S) :
    SameOrbit (unionStep f g S) a b ↔ SameOrbit f a b := by
  have memS : ∀ c d, SameOrbit (unionStep f g S) c d → (c ∈ S ↔ d ∈ S) := by
    intro c d h
    induction h with
    | rel c d h =>
      subst h
      unfold unionStep
      constructor
      · intro hc; rw [if_pos hc]; exact hf.maps c hc
      · intro hd
        by_contra hc
        rw [if_neg hc] at hd
        by_cases hcT : c ∈ T
        · exact Finset.disjoint_right.1 hST (hg.maps c hcT) hd
        · rw [hg.fix c hcT] at hd; exact hc hd
    | refl c => exact Iff.rfl
    | symm _ _ _ ih => exact ih.symm
    | trans _ _ _ _ _ ih₁ ih₂ => exact ih₁.trans ih₂
  have fwd : ∀ c d, SameOrbit (unionStep f g S) c d → c ∈ S → SameOrbit f c d := by
    intro c d h
    induction h with
    | rel c d h => intro hc; subst h; unfold unionStep; rw [if_pos hc]; exact SameOrbit.step _
    | refl c => intro _; exact SameOrbit.refl _ _
    | symm c d h ih => intro hd; exact (ih ((memS _ _ h).2 hd)).symm
    | trans c d e h₁ _ ih₁ ih₂ => intro hc; exact (ih₁ hc).trans (ih₂ ((memS _ _ h₁).1 hc))
  have bwd : ∀ c d, SameOrbit f c d → c ∈ S → SameOrbit (unionStep f g S) c d := by
    intro c d h
    induction h with
    | rel c d h =>
      intro hc; subst h
      have : unionStep f g S c = f c := by unfold unionStep; rw [if_pos hc]
      exact this ▸ SameOrbit.step _
    | refl c => intro _; exact SameOrbit.refl _ _
    | symm c d h ih => intro hd; exact (ih ((hf.mem_iff_of_sameOrbit h).2 hd)).symm
    | trans c d e h₁ _ ih₁ ih₂ =>
      intro hc; exact (ih₁ hc).trans (ih₂ ((hf.mem_iff_of_sameOrbit h₁).1 hc))
  exact ⟨fun h => fwd _ _ h ha, fun h => bwd _ _ h ha⟩

theorem unionStep_comm {f g : Nat → Nat} {S T : Finset Nat} (hST : Disjoint S T) :
    ∀ q ∈ S ∪ T, unionStep f g S q = unionStep g f T q := by
  intro q hq
  unfold unionStep
  rcases Finset.mem_union.1 hq with h | h
  · rw [if_pos h, if_neg (Finset.disjoint_left.1 hST h)]
  · rw [if_neg (Finset.disjoint_right.1 hST h), if_pos h]

theorem orbit_unionStep_left {f g : Nat → Nat} {S T : Finset Nat} (hf : IsPermOn f S)
    (hg : IsPermOn g T) (hST : Disjoint S T) {a : Nat} (ha : a ∈ S) :
    orbit (unionStep f g S) (S ∪ T) a = orbit f S a := by
  ext r
  simp only [mem_orbit]
  constructor
  · rintro ⟨-, h⟩
    have h' := (sameOrbit_unionStep_left hf hg hST ha).1 h
    exact ⟨(hf.mem_iff_of_sameOrbit h').1 ha, h'⟩
  · rintro ⟨hr, h⟩
    exact ⟨Finset.mem_union_left _ hr, (sameOrbit_unionStep_left hf hg hST ha).2 h⟩

theorem orbitCount_unionStep {f g : Nat → Nat} {S T : Finset Nat} (hf : IsPermOn f S)
    (hg : IsPermOn g T) (hST : Disjoint S T) :
    orbitCount (unionStep f g S) (S ∪ T) = orbitCount f S + orbitCount g T := by
  have hsym : ∀ a ∈ T, orbit (unionStep f g S) (S ∪ T) a = orbit g T a := by
    intro a ha
    have h1 : orbit (unionStep f g S) (S ∪ T) a = orbit (unionStep g f T) (T ∪ S) a := by
      ext r
      simp only [mem_orbit, Finset.union_comm S T]
      have hperm := isPermOn_unionStep hf hg hST
      have hperm' := isPermOn_unionStep hg hf hST.symm
      rw [Finset.union_comm] at hperm'
      have key : ∀ c d, SameOrbit (unionStep f g S) c d → c ∈ S ∪ T →
          SameOrbit (unionStep g f T) c d := by
        intro c d h
        induction h with
        | rel c d h => intro hc; subst h; rw [unionStep_comm hST c hc]; exact SameOrbit.step _
        | refl c => intro _; exact SameOrbit.refl _ _
        | symm c d h ih => intro hd; exact (ih ((hperm.mem_iff_of_sameOrbit h).2 hd)).symm
        | trans c d e h₁ _ ih₁ ih₂ =>
          intro hc; exact (ih₁ hc).trans (ih₂ ((hperm.mem_iff_of_sameOrbit h₁).1 hc))
      have key' : ∀ c d, SameOrbit (unionStep g f T) c d → c ∈ S ∪ T →
          SameOrbit (unionStep f g S) c d := by
        intro c d h
        induction h with
        | rel c d h => intro hc; subst h; rw [← unionStep_comm hST c hc]; exact SameOrbit.step _
        | refl c => intro _; exact SameOrbit.refl _ _
        | symm c d h ih => intro hd; exact (ih ((hperm'.mem_iff_of_sameOrbit h).2 hd)).symm
        | trans c d e h₁ _ ih₁ ih₂ =>
          intro hc; exact (ih₁ hc).trans (ih₂ ((hperm'.mem_iff_of_sameOrbit h₁).1 hc))
      exact and_congr_right fun _ => ⟨fun h => key _ _ h (Finset.mem_union_right _ ha),
        fun h => key' _ _ h (Finset.mem_union_right _ ha)⟩
    rw [h1]
    exact orbit_unionStep_left hg hf hST.symm ha
  unfold orbitCount
  rw [Finset.image_union,
    Finset.image_congr (g := orbit f S) (fun a ha => orbit_unionStep_left hf hg hST ha),
    Finset.image_congr (g := orbit g T) (fun a ha => hsym a ha)]
  apply Finset.card_union_of_disjoint
  rw [Finset.disjoint_left]
  intro A hA hA'
  obtain ⟨a, ha, rfl⟩ := Finset.mem_image.1 hA
  obtain ⟨b, hb, hab⟩ := Finset.mem_image.1 hA'
  have := mem_orbit_self (f := g) (S := T) hb
  rw [hab] at this
  exact Finset.disjoint_left.1 hST (orbit_subset this) hb

/-! ### Skipping a point -/

theorem isPermOn_delete {f : Nat → Nat} {S : Finset Nat} (hf : IsPermOn f S) {p : Nat}
    (hp : p ∈ S) (hfp : f p ≠ p) : IsPermOn (delete f p) (S.erase p) where
  maps q hq := by
    obtain ⟨hqp, hqS⟩ := Finset.mem_erase.1 hq
    unfold delete
    rw [if_neg hqp]
    split
    · exact Finset.mem_erase.2 ⟨hfp, hf.maps p hp⟩
    · exact Finset.mem_erase.2 ⟨‹_›, hf.maps q hqS⟩
  inj q hq r hr h := by
    obtain ⟨hqp, hqS⟩ := Finset.mem_erase.1 hq
    obtain ⟨hrp, hrS⟩ := Finset.mem_erase.1 hr
    unfold delete at h
    rw [if_neg hqp, if_neg hrp] at h
    split at h <;> split at h
    · exact hf.inj q hqS r hrS (‹f q = p›.trans ‹f r = p›.symm)
    · exact absurd (hf.inj p hp r hrS h) (Ne.symm hrp)
    · exact absurd (hf.inj q hqS p hp h) hqp
    · exact hf.inj q hqS r hrS h
  fix q hq := by
    unfold delete
    by_cases hqp : q = p
    · rw [if_pos hqp, hqp]
    · rw [if_neg hqp]
      have hqS : q ∉ S := fun h => hq (Finset.mem_erase.2 ⟨hqp, h⟩)
      rw [hf.fix q hqS, if_neg hqp]

theorem sameOrbit_delete {f : Nat → Nat} {p : Nat} (hfp : f p ≠ p) {a b : Nat} (ha : a ≠ p) :
    SameOrbit (delete f p) a b ↔ SameOrbit f a b ∧ b ≠ p := by
  have hne : ∀ c, c ≠ p → delete f p c ≠ p := by
    intro c hc
    unfold delete
    rw [if_neg hc]
    split
    · exact hfp
    · exact ‹_›
  have memp : ∀ c d, SameOrbit (delete f p) c d → (c = p ↔ d = p) := by
    intro c d h
    induction h with
    | rel c d h =>
      subst h
      constructor
      · intro hc; unfold delete; rw [if_pos hc]
      · intro hd; by_contra hc; exact hne c hc hd
    | refl c => exact Iff.rfl
    | symm _ _ _ ih => exact ih.symm
    | trans _ _ _ _ _ ih₁ ih₂ => exact ih₁.trans ih₂
  have fwd : ∀ c d, SameOrbit (delete f p) c d → c ≠ p → SameOrbit f c d := by
    intro c d h
    induction h with
    | rel c d h =>
      intro hc
      subst h
      unfold delete
      rw [if_neg hc]
      split
      · exact (SameOrbit.step (f := f) c).trans (by rw [‹f c = p›]; exact SameOrbit.step p)
      · exact SameOrbit.step c
    | refl c => intro _; exact SameOrbit.refl _ _
    | symm c d h ih => intro hd; exact (ih (fun hc => hd ((memp _ _ h).1 hc))).symm
    | trans c d e h₁ _ ih₁ ih₂ => intro hc; exact (ih₁ hc).trans (ih₂ (fun hd => hc ((memp _ _ h₁).2 hd)))
  constructor
  · intro h
    exact ⟨fwd _ _ h ha, fun hb => ha ((memp _ _ h).2 hb)⟩
  · rintro ⟨h, hb⟩
    let g : Nat → Nat := fun x => if x = p then f p else x
    have key : ∀ c d, SameOrbit f c d → SameOrbit (delete f p) (g c) (g d) := by
      intro c d h
      induction h with
      | rel c d h =>
        subst h
        by_cases hc : c = p
        · subst hc
          simp only [g, if_pos rfl, if_neg hfp]
          exact SameOrbit.refl _ _
        · by_cases hd : f c = p
          · simp only [g, if_neg hc, if_pos hd]
            have : delete f p c = f p := by unfold delete; rw [if_neg hc, if_pos hd]
            exact this ▸ SameOrbit.step c
          · simp only [g, if_neg hc, if_neg hd]
            have : delete f p c = f c := by unfold delete; rw [if_neg hc, if_neg hd]
            exact this ▸ SameOrbit.step c
      | refl c => exact SameOrbit.refl _ _
      | symm _ _ _ ih => exact ih.symm
      | trans _ _ _ _ _ ih₁ ih₂ => exact ih₁.trans ih₂
    have := key a b h
    simpa only [g, if_neg ha, if_neg hb] using this

theorem orbit_delete {f : Nat → Nat} {S : Finset Nat} {p : Nat} (hfp : f p ≠ p) {a : Nat} (ha : a ≠ p) :
    orbit (delete f p) (S.erase p) a = (orbit f S a).erase p := by
  ext r
  simp only [mem_orbit, Finset.mem_erase, sameOrbit_delete hfp ha]
  tauto

theorem orbitCount_delete {f : Nat → Nat} {S : Finset Nat} (hf : IsPermOn f S) {p : Nat}
    (hp : p ∈ S) (hfp : f p ≠ p) : orbitCount (delete f p) (S.erase p) = orbitCount f S := by
  unfold orbitCount
  have h1 : (S.erase p).image (orbit (delete f p) (S.erase p)) =
      ((S.erase p).image (orbit f S)).image (Finset.erase · p) := by
    rw [Finset.image_image]
    exact Finset.image_congr fun a ha => orbit_delete hfp (Finset.mem_erase.1 ha).1
  have h2 : (S.erase p).image (orbit f S) = S.image (orbit f S) := by
    apply Finset.Subset.antisymm (Finset.image_subset_image (Finset.erase_subset _ _))
    intro A hA
    obtain ⟨a, ha, rfl⟩ := Finset.mem_image.1 hA
    by_cases hap : a = p
    · subst hap
      rw [orbit_eq_of_sameOrbit (SameOrbit.step a)]
      exact Finset.mem_image_of_mem _ (Finset.mem_erase.2 ⟨hfp, hf.maps a hp⟩)
    · exact Finset.mem_image_of_mem _ (Finset.mem_erase.2 ⟨hap, ha⟩)
  rw [h1, h2]
  apply Finset.card_image_of_injOn
  intro A hA B hB hAB
  simp only [Finset.coe_image, Set.mem_image, Finset.mem_coe] at hA hB
  obtain ⟨a, ha, rfl⟩ := hA
  obtain ⟨b, hb, rfl⟩ := hB
  simp only at hAB
  rcases orbit_eq_or_disjoint (f := f) (S := S) a b with h | h
  · exact h
  · exfalso
    have hsub : (orbit f S a).erase p ⊆ orbit f S a ∩ orbit f S b := by
      intro s hs
      exact Finset.mem_inter.2 ⟨Finset.erase_subset _ _ hs, hAB ▸ hs |> Finset.erase_subset _ _⟩
    rw [Finset.disjoint_iff_inter_eq_empty.1 h, Finset.subset_empty] at hsub
    have ha' : a ∈ orbit f S a := mem_orbit_self ha
    have hap : a = p := by
      by_contra hap
      have : a ∈ (orbit f S a).erase p := Finset.mem_erase.2 ⟨hap, ha'⟩
      rw [hsub] at this
      exact absurd this (Finset.notMem_empty _)
    subst hap
    have : f a ∈ (orbit f S a).erase a :=
      Finset.mem_erase.2 ⟨hfp, mem_orbit.2 ⟨hf.maps a hp, SameOrbit.step a⟩⟩
    rw [hsub] at this
    exact absurd this (Finset.notMem_empty _)

/-! ### Exchanging two images -/

theorem swapImg_eq (f : Nat → Nat) (x y a : Nat) : swapImg f x y a = f (Equiv.swap x y a) := by
  unfold swapImg
  rw [Equiv.swap_apply_def]
  split
  · rfl
  split <;> rfl

theorem isPermOn_swapImg {f : Nat → Nat} {S : Finset Nat} (hf : IsPermOn f S) {x y : Nat}
    (hx : x ∈ S) (hy : y ∈ S) : IsPermOn (swapImg f x y) S := by
  have hsw : ∀ q ∈ S, Equiv.swap x y q ∈ S := by
    intro q hq
    rw [Equiv.swap_apply_def]
    split
    · exact hy
    split
    · exact hx
    · exact hq
  exact
    { maps := fun q hq => by rw [swapImg_eq]; exact hf.maps _ (hsw q hq)
      inj := fun q hq r hr h => by
        rw [swapImg_eq, swapImg_eq] at h
        exact (Equiv.swap x y).injective (hf.inj _ (hsw q hq) _ (hsw r hr) h)
      fix := fun q hq => by
        rw [swapImg_eq, Equiv.swap_apply_of_ne_of_ne (fun h : q = x => hq (by rw [h]; exact hx))
          (fun h : q = y => hq (by rw [h]; exact hy))]
        exact hf.fix q hq }

theorem sameOrbit_swapImg_self {f : Nat → Nat} {S : Finset Nat} (hf : IsPermOn f S) {x y : Nat}
    (hy : y ∈ S) (hxy : ¬SameOrbit f x y) : SameOrbit (swapImg f x y) x y := by
  set f' := swapImg f x y with hf'
  have hfy : f' x = f y := by simp [hf', swapImg]
  have key : ∀ c, Reach f (f y) c → Reach f' (f y) c ∨ Reach f' (f y) y := by
    intro c h
    induction h with
    | refl => exact Or.inl Relation.ReflTransGen.refl
    | @tail b c hc h ih =>
      rcases ih with ih | ih
      · by_cases hby : b = y
        · subst hby; exact Or.inr ih
        · left
          have hbx : b ≠ x := by
            intro hbx
            apply hxy
            have : SameOrbit f y b := (SameOrbit.step y).trans (IsPermOn.sameOrbit_of_reach hc)
            rw [hbx] at this
            exact this.symm
          have : f' b = f b := by simp [hf', swapImg, hbx, hby]
          exact ih.tail (this.trans h)
      · exact Or.inr ih
  have hret : Reach f (f y) y := hf.reach_symm hy (Relation.ReflTransGen.single rfl)
  have : Reach f' (f y) y := (key y hret).elim id id
  exact (hfy ▸ SameOrbit.step (f := f') x).trans (IsPermOn.sameOrbit_of_reach this)

theorem sameOrbit_swapImg {f : Nat → Nat} {S : Finset Nat} (hf : IsPermOn f S) {x y : Nat}
    (hx : x ∈ S) (hy : y ∈ S) (hxy : ¬SameOrbit f x y) {a b : Nat} :
    SameOrbit (swapImg f x y) a b ↔
      SameOrbit f a b ∨ ((SameOrbit f x a ∨ SameOrbit f y a) ∧ (SameOrbit f x b ∨ SameOrbit f y b)) := by
  set f' := swapImg f x y with hf'
  have hxy' := sameOrbit_swapImg_self hf hy hxy
  have hfx : f' x = f y := by simp [hf', swapImg]
  have hfy : f' y = f x := by simp [hf', swapImg, Ne.symm (show x ≠ y from fun h => hxy (h ▸ SameOrbit.refl f x))]
  have hfx' : f x = f' y := hfy.symm
  have hfy' : f y = f' x := hfx.symm
  constructor
  · intro h
    induction h with
    | rel c d h =>
      subst h
      by_cases hcx : c = x
      · subst hcx; right; exact ⟨Or.inl (SameOrbit.refl _ _), Or.inr (hfx ▸ SameOrbit.step (f := f) y)⟩
      by_cases hcy : c = y
      · subst hcy; right; exact ⟨Or.inr (SameOrbit.refl _ _), Or.inl (hfy ▸ SameOrbit.step (f := f) x)⟩
      · left
        have : f' c = f c := by simp [hf', swapImg, hcx, hcy]
        exact this ▸ SameOrbit.step c
    | refl c => left; exact SameOrbit.refl _ _
    | symm c d _ ih =>
      rcases ih with ih | ⟨h₁, h₂⟩
      · exact Or.inl ih.symm
      · exact Or.inr ⟨h₂, h₁⟩
    | trans c d e _ _ ih₁ ih₂ =>
      rcases ih₁ with ih₁ | ⟨h₁, h₂⟩ <;> rcases ih₂ with ih₂ | ⟨h₃, h₄⟩
      · exact Or.inl (ih₁.trans ih₂)
      · exact Or.inr ⟨h₃.imp (fun h => h.trans ih₁.symm) (fun h => h.trans ih₁.symm), h₄⟩
      · exact Or.inr ⟨h₁, h₂.imp (fun h => h.trans ih₂) (fun h => h.trans ih₂)⟩
      · exact Or.inr ⟨h₁, h₄⟩
  · have stepx : SameOrbit f' x (f x) := hfy ▸ hxy'.trans (SameOrbit.step (f := f') y)
    have stepy : SameOrbit f' y (f y) := hfx ▸ hxy'.symm.trans (SameOrbit.step (f := f') x)
    have lift : ∀ c d, SameOrbit f c d → SameOrbit f' c d := by
      intro c d h
      induction h with
      | rel c d h =>
        subst h
        by_cases hcx : c = x
        · subst hcx; exact stepx
        by_cases hcy : c = y
        · subst hcy; exact stepy
        · have : f' c = f c := by simp [hf', swapImg, hcx, hcy]
          exact this ▸ SameOrbit.step c
      | refl c => exact SameOrbit.refl _ _
      | symm _ _ _ ih => exact ih.symm
      | trans _ _ _ _ _ ih₁ ih₂ => exact ih₁.trans ih₂
    rintro (h | ⟨h₁, h₂⟩)
    · exact lift _ _ h
    · have ha : SameOrbit f' x a := h₁.elim (fun h => lift _ _ h) (fun h => hxy'.trans (lift _ _ h))
      have hb : SameOrbit f' x b := h₂.elim (fun h => lift _ _ h) (fun h => hxy'.trans (lift _ _ h))
      exact ha.symm.trans hb

theorem orbitCount_swapImg {f : Nat → Nat} {S : Finset Nat} (hf : IsPermOn f S) {x y : Nat}
    (hx : x ∈ S) (hy : y ∈ S) (hxy : ¬SameOrbit f x y) :
    orbitCount (swapImg f x y) S + 1 = orbitCount f S := by
  set f' := swapImg f x y with hf'
  set M := orbit f S x ∪ orbit f S y with hM
  have hMS : M ⊆ S := Finset.union_subset orbit_subset orbit_subset
  have memM : ∀ r, r ∈ M ↔ r ∈ S ∧ (SameOrbit f x r ∨ SameOrbit f y r) := by
    intro r; simp only [hM, Finset.mem_union, mem_orbit]; tauto
  have horb : ∀ q ∈ S, orbit f' S q = if q ∈ M then M else orbit f S q := by
    intro q hq
    ext r
    rw [mem_orbit, sameOrbit_swapImg hf hx hy hxy]
    split
    · rename_i hqM
      rw [memM] at hqM ⊢
      constructor
      · rintro ⟨hr, h | ⟨-, h⟩⟩
        · exact ⟨hr, hqM.2.imp (fun h' => h'.trans h) (fun h' => h'.trans h)⟩
        · exact ⟨hr, h⟩
      · rintro ⟨hr, h⟩
        exact ⟨hr, Or.inr ⟨hqM.2, h⟩⟩
    · rename_i hqM
      rw [memM] at hqM
      rw [mem_orbit]
      constructor
      · rintro ⟨hr, h | ⟨h, -⟩⟩
        · exact ⟨hr, h⟩
        · exact absurd ⟨hq, h⟩ hqM
      · rintro ⟨hr, h⟩
        exact ⟨hr, Or.inl h⟩
  have hsplit : S = (S \ M) ∪ M := (Finset.sdiff_union_of_subset hMS).symm
  have hA : ∀ q ∈ S \ M, orbit f' S q = orbit f S q := by
    intro q hq
    rw [horb q (Finset.mem_sdiff.1 hq).1, if_neg (Finset.mem_sdiff.1 hq).2]
  have hB : ∀ q ∈ M, orbit f' S q = M := by
    intro q hq
    rw [horb q (hMS hq), if_pos hq]
  have hxM : x ∈ M := Finset.mem_union_left _ (mem_orbit_self hx)
  have hyM : y ∈ M := Finset.mem_union_right _ (mem_orbit_self hy)
  have hMne : M.Nonempty := ⟨x, hxM⟩
  have img' : S.image (orbit f' S) = (S \ M).image (orbit f S) ∪ {M} := by
    have : S.image (orbit f' S) = ((S \ M) ∪ M).image (orbit f' S) := by rw [← hsplit]
    rw [this, Finset.image_union, Finset.image_congr (g := orbit f S) (fun q hq => hA q hq),
      Finset.image_congr (g := fun _ => M) (fun q hq => hB q hq), Finset.image_const hMne]
  have img : S.image (orbit f S) = (S \ M).image (orbit f S) ∪ {orbit f S x, orbit f S y} := by
    have : S.image (orbit f S) = ((S \ M) ∪ M).image (orbit f S) := by rw [← hsplit]
    rw [this, Finset.image_union]
    congr 1
    ext A
    simp only [Finset.mem_image, Finset.mem_insert, Finset.mem_singleton]
    constructor
    · rintro ⟨q, hq, rfl⟩
      rcases (memM q).1 hq with ⟨-, h | h⟩
      · exact Or.inl (orbit_eq_of_sameOrbit h).symm
      · exact Or.inr (orbit_eq_of_sameOrbit h).symm
    · rintro (rfl | rfl)
      · exact ⟨x, hxM, rfl⟩
      · exact ⟨y, hyM, rfl⟩
  have notA : ∀ q, orbit f S q ∈ (S \ M).image (orbit f S) → q ∈ S → q ∉ M := by
    intro q hq hqS hqM
    obtain ⟨q', hq', h⟩ := Finset.mem_image.1 hq
    have : SameOrbit f q' q := sameOrbit_of_orbit_eq (Finset.mem_sdiff.1 hq').1 h
    apply (Finset.mem_sdiff.1 hq').2
    rw [memM] at hqM ⊢
    exact ⟨(Finset.mem_sdiff.1 hq').1, hqM.2.imp (fun h => h.trans this.symm) (fun h => h.trans this.symm)⟩
  unfold orbitCount
  rw [img', img, Finset.card_union_of_disjoint, Finset.card_union_of_disjoint, Finset.card_singleton,
    Finset.card_pair]
  · intro h
    exact hxy (sameOrbit_of_orbit_eq hx h)
  · rw [Finset.disjoint_left]
    intro A hA hA'
    simp only [Finset.mem_insert, Finset.mem_singleton] at hA'
    rcases hA' with rfl | rfl
    · exact notA x hA hx hxM
    · exact notA y hA hy hyM
  · rw [Finset.disjoint_left]
    intro A hA hA'
    rw [Finset.mem_singleton] at hA'
    subst hA'
    rw [hM, ← orbit_eq_of_sameOrbit (f := f) (S := S) (SameOrbit.refl f x)] at hA
    obtain ⟨q', hq', h⟩ := Finset.mem_image.1 hA
    have hq'M : q' ∈ orbit f S x ∪ orbit f S y := h ▸ mem_orbit_self (Finset.mem_sdiff.1 hq').1
    exact (Finset.mem_sdiff.1 hq').2 hq'M

/-! ### Invariants, congruence, two swaps -/

theorem sameOrbit_invariant {f : Nat → Nat} {α : Type} (ℓ : Nat → α) (hf : ∀ a, ℓ (f a) = ℓ a)
    {a b : Nat} (h : SameOrbit f a b) : ℓ a = ℓ b := by
  induction h with
  | rel a b h => subst h; exact (hf a).symm
  | refl => rfl
  | symm _ _ _ ih => exact ih.symm
  | trans _ _ _ _ _ ih₁ ih₂ => exact ih₁.trans ih₂

theorem sameOrbit_mono_of_eqOn {f g : Nat → Nat} (C : Nat → Prop) (hf : ∀ a, C a ↔ C (f a))
    (heq : ∀ a, C a → f a = g a) {a b : Nat} (ha : C a) (h : SameOrbit f a b) :
    SameOrbit g a b := by
  have memC : ∀ c d, SameOrbit f c d → (C c ↔ C d) := by
    intro c d h
    induction h with
    | rel c d h => subst h; exact hf c
    | refl => exact Iff.rfl
    | symm _ _ _ ih => exact ih.symm
    | trans _ _ _ _ _ ih₁ ih₂ => exact ih₁.trans ih₂
  have key : ∀ c d, SameOrbit f c d → C c → SameOrbit g c d := by
    intro c d h
    induction h with
    | rel c d h => intro hc; subst h; rw [heq c hc]; exact SameOrbit.step _
    | refl c => intro _; exact SameOrbit.refl _ _
    | symm c d h ih => intro hd; exact (ih ((memC _ _ h).2 hd)).symm
    | trans c d e h₁ _ ih₁ ih₂ => intro hc; exact (ih₁ hc).trans (ih₂ ((memC _ _ h₁).1 hc))
  exact key a b h ha

theorem sameOrbit_congr {f g : Nat → Nat} (C : Nat → Prop) (hf : ∀ a, C a ↔ C (f a))
    (hg : ∀ a, C a ↔ C (g a)) (heq : ∀ a, C a → f a = g a) {a b : Nat} (ha : C a) :
    SameOrbit f a b ↔ SameOrbit g a b :=
  ⟨sameOrbit_mono_of_eqOn C hf heq ha,
    sameOrbit_mono_of_eqOn C hg (fun a h => (heq a h).symm) ha⟩

theorem unionStep_swap {f g : Nat → Nat} {S T : Finset Nat} (hf : IsPermOn f S)
    (hg : IsPermOn g T) (hST : Disjoint S T) : unionStep f g S = unionStep g f T := by
  funext q
  unfold unionStep
  by_cases hS : q ∈ S
  · rw [if_pos hS, if_neg (Finset.disjoint_left.1 hST hS)]
  by_cases hT : q ∈ T
  · rw [if_neg hS, if_pos hT]
  · rw [if_neg hS, if_neg hT, hf.fix q hS, hg.fix q hT]

theorem sameOrbit_unionStep_right {f g : Nat → Nat} {S T : Finset Nat} (hf : IsPermOn f S)
    (hg : IsPermOn g T) (hST : Disjoint S T) {a b : Nat} (ha : a ∈ T) :
    SameOrbit (unionStep f g S) a b ↔ SameOrbit g a b := by
  rw [unionStep_swap hf hg hST]
  exact sameOrbit_unionStep_left hg hf hST.symm ha

theorem not_sameOrbit_swapImg_of {f : Nat → Nat} {S : Finset Nat} (hf : IsPermOn f S)
    {x₁ y₁ x₂ y₂ : Nat} (hx₁ : x₁ ∈ S) (hy₁ : y₁ ∈ S) (h₁ : ¬SameOrbit f x₁ y₁)
    (h₂ : ¬SameOrbit f x₂ y₂) (hxy : ¬SameOrbit f x₁ y₂) (hyx : ¬SameOrbit f y₁ x₂)
    (h : ¬(SameOrbit f x₁ x₂ ∧ SameOrbit f y₁ y₂)) : ¬SameOrbit (swapImg f x₁ y₁) x₂ y₂ := by
  rw [sameOrbit_swapImg hf hx₁ hy₁ h₁]
  rintro (h' | ⟨hx | hx, hy | hy⟩)
  · exact h₂ h'
  · exact hxy hy
  · exact h ⟨hx, hy⟩
  · exact hyx hx
  · exact hyx hx

theorem swapImg_left (f : Nat → Nat) (x y : Nat) : swapImg f x y x = f y := by simp [swapImg]
theorem swapImg_right (f : Nat → Nat) {x y : Nat} (h : y ≠ x) : swapImg f x y y = f x := by
  simp [swapImg, h]
theorem swapImg_of_ne (f : Nat → Nat) {x y a : Nat} (hx : a ≠ x) (hy : a ≠ y) :
    swapImg f x y a = f a := by simp [swapImg, hx, hy]

/-! ### Skipping a list of points -/

theorem skip_nil (f : Nat → Nat) : skip f [] = f := by funext a; simp [skip]

theorem skip_cons (f : Nat → Nat) (p : Nat) (P : List Nat) (hp : f p ∉ p :: P)
    (hP : ∀ q ∈ P, f q ∉ p :: P) (hnd : p ∉ P) : skip f (p :: P) = skip (delete f p) P := by
  funext a
  unfold skip delete
  by_cases hap : a = p
  · subst hap
    simp [hnd]
  by_cases haP : a ∈ P
  · simp [haP]
  have ha : a ∉ p :: P := by simp [hap, haP]
  simp only [ha, if_false, haP, hap]
  by_cases hfa : f a = p
  · have : f p ∉ P := fun h => hp (List.mem_cons_of_mem _ h)
    simp [hfa, this]
  by_cases hfaP : f a ∈ P
  · have h1 : f a ∈ p :: P := List.mem_cons_of_mem _ hfaP
    have h2 : f (f a) ≠ p := fun h => hP _ hfaP (h ▸ List.mem_cons_self)
    simp [h1, hfaP, hfa, h2]
  · have h1 : f a ∉ p :: P := by simp [hfa, hfaP]
    simp [h1, hfaP, hfa]

theorem orbitCount_skip {f : Nat → Nat} {S : Finset Nat} (hf : IsPermOn f S) (P : List Nat)
    (hP : ∀ p ∈ P, p ∈ S) (hfP : ∀ p ∈ P, f p ∉ P) (hnd : P.Nodup) :
    IsPermOn (skip f P) (S \ P.toFinset) ∧
      orbitCount (skip f P) (S \ P.toFinset) = orbitCount f S := by
  induction P generalizing f S with
  | nil => simpa [skip_nil] using hf
  | cons p P ih =>
    have hpP : p ∉ P := (List.nodup_cons.1 hnd).1
    have hfp : f p ∉ p :: P := hfP p List.mem_cons_self
    have hfp' : f p ≠ p := fun h => hfp (by rw [h]; exact List.mem_cons_self)
    rw [skip_cons f p P hfp (fun q hq => hfP q (List.mem_cons_of_mem _ hq)) hpP]
    have hS : S.erase p \ P.toFinset = S \ (p :: P).toFinset := by
      ext a; simp only [Finset.mem_sdiff, Finset.mem_erase, List.toFinset_cons,
        Finset.mem_insert, List.mem_toFinset]; tauto
    rw [← hS]
    have hperm := isPermOn_delete hf (hP p List.mem_cons_self) hfp'
    obtain ⟨h1, h2⟩ := ih hperm
      (fun q hq => Finset.mem_erase.2 ⟨fun h => hpP (h ▸ hq), hP q (List.mem_cons_of_mem _ hq)⟩)
      (fun q hq h => by
        unfold delete at h
        have hqp : q ≠ p := fun h => hpP (h ▸ hq)
        rw [if_neg hqp] at h
        split at h
        · exact hfp (List.mem_cons_of_mem _ h)
        · exact hfP q (List.mem_cons_of_mem _ hq) (List.mem_cons_of_mem _ h))
      (List.nodup_cons.1 hnd).2
    exact ⟨h1, h2.trans (orbitCount_delete hf (hP p List.mem_cons_self) hfp')⟩

theorem sameOrbit_skip {f : Nat → Nat} (P : List Nat) (hfP : ∀ p ∈ P, f p ∉ P) (hnd : P.Nodup)
    {a b : Nat} (ha : a ∉ P) : SameOrbit (skip f P) a b ↔ SameOrbit f a b ∧ b ∉ P := by
  induction P generalizing f with
  | nil => simp [skip_nil]
  | cons p P ih =>
    have hpP : p ∉ P := (List.nodup_cons.1 hnd).1
    have hfp : f p ∉ p :: P := hfP p List.mem_cons_self
    have hfp' : f p ≠ p := fun h => hfp (by rw [h]; exact List.mem_cons_self)
    rw [skip_cons f p P hfp (fun q hq => hfP q (List.mem_cons_of_mem _ hq)) hpP]
    have haP : a ∉ P := fun h => ha (List.mem_cons_of_mem _ h)
    have hap : a ≠ p := fun h => ha (h ▸ List.mem_cons_self)
    rw [ih (fun q hq h => by
        unfold delete at h
        have hqp : q ≠ p := fun h => hpP (h ▸ hq)
        rw [if_neg hqp] at h
        split at h
        · exact hfp (List.mem_cons_of_mem _ h)
        · exact hfP q (List.mem_cons_of_mem _ hq) (List.mem_cons_of_mem _ h))
      (List.nodup_cons.1 hnd).2 haP, sameOrbit_delete hfp' hap]
    simp only [List.mem_cons, not_or]
    tauto

end Spqr
