import Spqr.PlanarEmbedLink
import Spqr.Proofs.PieceAppend

namespace Spqr.PlanarSpqrTree

open EmbedM

def closeOuter (j : Nat) : EmbedM Unit := do
  let a ← outer j 0
  if a.isSome then link a (← outer j 1)

theorem outer_some_iff (s : EmbedState) (j k q : Nat) :
    s.outerE[j]![k]! = some q ↔
      s.outerE[j]?.bind (fun o => o[k]?) = some (some q) := by
  simp only [getElem!_def]
  cases hj : s.outerE[j]? with
  | none => simp [show (default : Array (Option Nat)) = #[] from rfl]
  | some o =>
    simp only [Option.bind_some]
    cases hk : o[k]? <;> simp

theorem closeOuter_run (j : Nat) (s : EmbedState) :
    (closeOuter j).run s = (link s.outerE[j]![0]! s.outerE[j]![1]!).run s := by
  unfold closeOuter
  simp only [StateT.run_bind, outer_run]
  cases s.outerE[j]![0]! <;> rfl

theorem closeOuter_outerE (j : Nat) (s : EmbedState) :
    ((closeOuter j).run s).2.outerE = s.outerE := by
  rw [closeOuter_run, link_outerE]

theorem closeOuter_rotAdj_size (j : Nat) (s : EmbedState) :
    ((closeOuter j).run s).2.rotAdj.size = s.rotAdj.size := by
  rw [closeOuter_run, link_rotAdj_size]

theorem link_rotAdj_get_ne (a b : Option Nat) (s : EmbedState) (q : Nat)
    (ha : a ≠ some q) (hb : b ≠ some q) :
    ((link a b).run s).2.rotAdj[q]? = s.rotAdj[q]? := by
  cases a with
  | none => rfl
  | some a =>
    cases b with
    | none => rfl
    | some b =>
      rw [link_some_run]
      simp only [Array.set!, Array.getElem?_setIfInBounds, Array.size_setIfInBounds]
      simp_all

theorem closeOuter_frame (j : Nat) (s : EmbedState) (q : Nat)
    (h : ¬s.exposedAt j q) : ((closeOuter j).run s).2.rotAdj[q]? = s.rotAdj[q]? := by
  rw [closeOuter_run]
  exact link_rotAdj_get_ne _ _ _ _
    (fun hh => h ⟨0, (outer_some_iff s j 0 q).1 hh⟩)
    (fun hh => h ⟨1, (outer_some_iff s j 1 q).1 hh⟩)

theorem closeOuter_spec (P : Piece) (ρ : RotationSystem) (j : Nat) (s : EmbedState)
    (hρ : IsPlanarEmbedding P.es P.nVerts ρ)
    (hbound : ∀ q, P.Mem q → q < s.rotAdj.size)
    (hagree : P.Agrees s.rotAdj ρ)
    (hopen : ∀ q, P.Mem q → (s.rotAdj[q]? = some none ↔ s.exposedAt j q))
    (hslots : ∀ k q, s.outerE[j]?.bind (fun o => o[k]?) = some (some q) → k < 2)
    (hpair : (∀ a, s.outerE[j]?.bind (fun o => o[0]?) = some (some a) →
      ∃ b, s.outerE[j]?.bind (fun o => o[1]?) = some (some b) ∧
        ∃ la lb, P.loc a = some la ∧ P.loc b = some lb ∧ ρ.get la = some lb) ∧
      (∀ b, s.outerE[j]?.bind (fun o => o[1]?) = some (some b) →
        ∃ a, s.outerE[j]?.bind (fun o => o[0]?) = some (some a))) :
    P.Agrees ((closeOuter j).run s).2.rotAdj ρ ∧
    (∀ q, P.Mem q → ((closeOuter j).run s).2.rotAdj[q]? ≠ some none) ∧
    (∀ q, ¬P.Mem q → ((closeOuter j).run s).2.rotAdj[q]? = s.rotAdj[q]?) := by
  have hem : ∀ q, s.exposedAt j q → P.Mem q := by
    rintro q ⟨k, hk⟩
    have hkl := hslots k q hk
    interval_cases k
    · obtain ⟨b, _, la, lb, ha, _, _⟩ := hpair.1 q hk
      exact Piece.mem_of_loc ha
    · obtain ⟨a, ha⟩ := hpair.2 q hk
      obtain ⟨b, hb, la, lb, _, hlb, _⟩ := hpair.1 a ha
      have : b = q := Option.some.inj (Option.some.inj (hb.symm.trans hk))
      subst b
      exact Piece.mem_of_loc hlb
  refine ⟨?_, ?_, fun q hq => closeOuter_frame j s q (fun hh => hq (hem q hh))⟩
  all_goals
    cases ha : s.outerE[j]![0]! with
    | none =>
      have hn : ∀ q, ¬s.exposedAt j q := by
        rintro q ⟨k, hk⟩
        have hkl := hslots k q hk
        interval_cases k
        · have hh := (outer_some_iff s j 0 q).2 hk
          rw [ha] at hh; cases hh
        · obtain ⟨a, ha'⟩ := hpair.2 q hk
          have hh := (outer_some_iff s j 0 a).2 ha'
          rw [ha] at hh; cases hh
      first
      | intro q r hq hr
        rw [closeOuter_frame j s q (hn q)] at hr
        exact hagree q r hq hr
      | intro q hq hh
        rw [closeOuter_frame j s q (hn q)] at hh
        exact hn q ((hopen q hq).1 hh)
    | some a =>
      obtain ⟨b, hb, la, lb, hla, hlb, hab⟩ := hpair.1 a ((outer_some_iff s j 0 a).1 ha)
      have hb' := (outer_some_iff s j 1 b).2 hb
      have hal : la < ρ.size := by rw [hρ.size]; simpa [Piece.es] using Piece.loc_lt hla
      have hba := (hρ.involution la hal lb hab).2
      have hne : a ≠ b := by
        intro heq
        subst b
        have hll : la = lb := Option.some.inj (hla.symm.trans hlb)
        exact hρ.opposite_dir la hal lb hab (congrArg QE.dir hll.symm)
      have hget : ∀ q, ((closeOuter j).run s).2.rotAdj[q]? =
          if q = a then some (some b) else if q = b then some (some a) else s.rotAdj[q]? := by
        intro q
        rw [closeOuter_run, ha, hb']
        exact link_rotAdj_get a b s (hbound a (Piece.mem_of_loc hla))
          (hbound b (Piece.mem_of_loc hlb)) hne q
      first
      | intro q r hq hr
        rw [hget] at hr
        split_ifs at hr with hqa hqb
        · subst q
          have : r = b := (Option.some.inj (Option.some.inj hr)).symm
          subst r
          exact ⟨la, lb, hla, hlb, hab⟩
        · subst q
          have : r = a := (Option.some.inj (Option.some.inj hr)).symm
          subst r
          exact ⟨lb, la, hlb, hla, hba⟩
        · exact hagree q r hq hr
      | intro q hq hh
        rw [hget] at hh
        split_ifs at hh with hqa hqb
        · cases hh
        · cases hh
        · obtain ⟨k, hk⟩ := (hopen q hq).1 hh
          have hkl := hslots k q hk
          interval_cases k
          · have hqa' := (outer_some_iff s j 0 q).2 hk
            exact hqa (Option.some.inj (hqa'.symm.trans ha))
          · have hqb' := (outer_some_iff s j 1 q).2 hk
            exact hqb (Option.some.inj (hqb'.symm.trans hb'))

end Spqr.PlanarSpqrTree
