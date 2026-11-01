import Spqr.PlanarEmbedVBoundary

namespace Spqr.PlanarSpqrTree

open EmbedM

def vLoop : List Nat → EmbedState → Option Nat → Option Nat →
    (Option Nat × Option Nat) × EmbedState
  | [], s, a, b => ((a, b), s)
  | j :: js, s, a, b =>
    if a.isNone then vLoop js s s.outerE[j]![0]! s.outerE[j]![1]!
    else
      let s' := ((link b s.outerE[j]![0]!).run s).2
      vLoop js s' a s'.outerE[j]![1]!

theorem vLoop_outerE (js : List Nat) (s : EmbedState) (a b : Option Nat) :
    (vLoop js s a b).2.outerE = s.outerE := by
  induction js generalizing s a b with
  | nil => rfl
  | cons j js ih =>
    simp only [vLoop]
    split <;> rw [ih]
    exact link_outerE ..

theorem vLoop_rotAdj_size (js : List Nat) (s : EmbedState) (a b : Option Nat) :
    (vLoop js s a b).2.rotAdj.size = s.rotAdj.size := by
  induction js generalizing s a b with
  | nil => rfl
  | cons j js ih =>
    simp only [vLoop]
    split <;> rw [ih]
    exact link_rotAdj_size ..

theorem vLoop_left (js : List Nat) (s : EmbedState) (a : Nat) (b : Option Nat) :
    (vLoop js s (some a) b).1.1 = some a := by
  induction js generalizing s b with
  | nil => rfl
  | cons j js ih => simp only [vLoop, Option.isNone_some, Bool.false_eq_true, ↓reduceIte]; exact ih ..

def setOuterPair (i : Nat) (a b : Option Nat) (s : EmbedState) : EmbedState :=
  ((setOuter i 1 b).run ((setOuter i 0 a).run s).2).2

theorem setOuterPair_rotAdj (i : Nat) (a b : Option Nat) (s : EmbedState) :
    (setOuterPair i a b s).rotAdj = s.rotAdj := rfl

theorem setOuterPair_outer_size (i : Nat) (a b : Option Nat) (s : EmbedState) :
    (setOuterPair i a b s).outerE.size = s.outerE.size := by
  simp only [setOuterPair, setOuter_outerE_size]

theorem setOuterPair_outer_ne (i : Nat) (a b : Option Nat) (s : EmbedState)
    (j : Nat) (hji : j ≠ i) : (setOuterPair i a b s).outerE[j]? = s.outerE[j]? := by
  simp only [setOuterPair, setOuter_outerE_ne i 1 b _ j hji, setOuter_outerE_ne i 0 a s j hji]

theorem setOuterPair_outer_eq (i : Nat) (a b : Option Nat) (s : EmbedState)
    (hi : i < s.outerE.size) :
    (setOuterPair i a b s).outerE[i]? = some ((s.outerE[i]!.set! 0 a).set! 1 b) := by
  simp [setOuterPair, setOuter_run, Array.getElem_modify_self, hi]

theorem setOuterPair_lookup (i : Nat) (a b : Option Nat) (s : EmbedState)
    (hi : i < s.outerE.size) (hs : s.outerE[i]!.size = 4) (j k : Nat) :
    (setOuterPair i a b s).outerE[j]?.bind (fun o => o[k]?) =
      if j = i then
        if k = 0 then some a else if k = 1 then some b
        else s.outerE[j]?.bind (fun o => o[k]?)
      else s.outerE[j]?.bind (fun o => o[k]?) := by
  by_cases hji : j = i
  · subst j
    rw [setOuterPair_outer_eq _ _ _ _ hi]
    simp only [Option.bind_some, ↓reduceIte, Array.set!, Array.getElem?_setIfInBounds,
      Array.size_setIfInBounds, hs]
    rw [Array.getElem?_eq_getElem hi]
    simp only [Option.bind_some, getElem!_pos s.outerE i hi]
    by_cases hk0 : k = 0 <;> by_cases hk1 : k = 1 <;> simp_all [eq_comm] <;> omega
  · rw [setOuterPair_outer_ne _ _ _ _ j hji, ite_eq_right hji]

theorem setOuterPair_row_size (i : Nat) (a b : Option Nat) (s : EmbedState) (j : Nat) :
    (setOuterPair i a b s).outerE[j]!.size = s.outerE[j]!.size := by
  by_cases hji : j = i
  · subst j
    by_cases hi : i < s.outerE.size
    · simp only [getElem!_def, setOuterPair_outer_eq _ _ _ _ hi]
      simp
    · simp [getElem!_def, Array.getElem?_eq_none (by omega : s.outerE.size ≤ i),
        Array.getElem?_eq_none (by rw [setOuterPair_outer_size]; omega :
          (setOuterPair i a b s).outerE.size ≤ i)]
  · simp only [getElem!_def, setOuterPair_outer_ne _ _ _ _ j hji]

theorem embedItem_V (t : PlanarSpqrTree) (i : Nat) (ht : t.types[i]! = .V) (s : EmbedState) :
    ((t.embedItem i).run s).2 =
      let r := vLoop (t.children i) s none none
      setOuterPair i r.1.1 r.1.2 r.2 := by
  have idbind {α β : Type} (a : Id α) (f : α → Id β) : (a >>= f) = f a := rfl
  have idpure {α : Type} (a : α) : (pure a : Id α) = a := rfl
  have idmap {α β : Type} (a : Id α) (f : α → Id β) : (f <$> a) = f a := rfl
  have loop : ∀ (js : List Nat) (s : EmbedState) (qes : Option Nat × Option Nat),
      ((do
        let mut qes := qes
        for j in js do
          if qes.1.isNone then qes := (← outer j 0, ← outer j 1)
          else
            link qes.2 (← outer j 0)
            qes := (qes.1, ← outer j 1)
        pure qes : EmbedM (Option Nat × Option Nat)).run s) = vLoop js s qes.1 qes.2 := by
    intro js
    induction js with
    | nil => intros; rfl
    | cons j js ih =>
      intro s qes
      cases qes with
      | mk a b =>
        have hh := ih (if a.isNone then s else ((link b s.outerE[j]![0]!).run s).2)
          (if a.isNone then (s.outerE[j]![0]!, s.outerE[j]![1]!)
            else (a, s.outerE[j]![1]!))
        cases a <;> cases b <;> cases h0 : s.outerE[j]![0]! <;>
          simpa [List.forIn_cons, vLoop, StateT.run_bind, outer_run, link_none_left,
            link_none_right, link_some_run, h0, idbind, idpure, idmap] using hh
  unfold embedItem
  rw [ht]
  change ((do
    let qes ← (do
      let mut qes : Option Nat × Option Nat := (none, none)
      for j in t.children i do
        if qes.1.isNone then qes := (← outer j 0, ← outer j 1)
        else
          link qes.2 (← outer j 0)
          qes := (qes.1, ← outer j 1)
      pure qes)
    setOuter i 0 qes.1
    setOuter i 1 qes.2).run s).2 = _
  simp only [StateT.run_bind, loop, idbind, setOuterPair]

theorem vLoop_open_spec (P : Piece) (f : Nat → List Nat) (js : List Nat)
    (s : EmbedState) (a b v : Nat) (ρ₀ : RotationSystem)
    (hcur : P.OpenEmbedding s.rotAdj ρ₀ a b v)
    (hchild : ∀ j ∈ js, ∃ c d ρ, s.outerE[j]![0]! = some c ∧ s.outerE[j]![1]! = some d ∧
      ({P with ves := f j} : Piece).OpenEmbedding s.rotAdj ρ c d v)
    (hcursep : ∀ j ∈ js, List.Disjoint P.ves (f j) ∧
      ∀ w, HasEdge P.es w → HasEdge ({P with ves := f j} : Piece).es w → w = v)
    (hpairs : js.Pairwise fun j k => List.Disjoint (f j) (f k) ∧
      ∀ w, HasEdge ({P with ves := f j} : Piece).es w →
        HasEdge ({P with ves := f k} : Piece).es w → w = v)
    (hbound : ∀ q, ({P with ves := P.ves ++ js.flatMap f} : Piece).Mem q → q < s.rotAdj.size) :
    let r := vLoop js s (some a) (some b)
    (∃ d ρ, r.1 = (some a, some d) ∧
      ({P with ves := P.ves ++ js.flatMap f} : Piece).OpenEmbedding r.2.rotAdj ρ a d v) ∧
    (∀ q, ¬({P with ves := P.ves ++ js.flatMap f} : Piece).Mem q →
      r.2.rotAdj[q]? = s.rotAdj[q]?) := by
  induction js generalizing P s b ρ₀ with
  | nil =>
    simp only [vLoop, List.flatMap_nil, List.append_nil]
    exact ⟨⟨b, ρ₀, rfl, hcur⟩, fun _ _ => trivial⟩
  | cons j js ih =>
    obtain ⟨c, d, ρj, hc, hd, hj⟩ := hchild j (by simp)
    have hsep := hcursep j (by simp)
    obtain ⟨hpj, hpjs⟩ := List.pairwise_cons.1 hpairs
    let P' : Piece := {P with ves := P.ves ++ f j}
    let s' := ((link (some b) (some c)).run s).2
    have hsmall : ∀ q, P'.Mem q → q < s.rotAdj.size := by
      intro q hq
      apply hbound q
      simp only [Piece.Mem, List.flatMap_cons, List.mem_append, P'] at hq ⊢
      tauto
    obtain ⟨ρ', hρ'⟩ := hcur.splice hj hsep.1 hsep.2 hsmall
    have hbcMem : P.Mem b ∧ ({P with ves := f j} : Piece).Mem c := by
      obtain ⟨_, _, _, hb, _⟩ := hcur.boundary
      obtain ⟨_, _, hc', _⟩ := hj.boundary
      exact ⟨Piece.mem_of_loc hb, Piece.mem_of_loc hc'⟩
    have hnext : ∀ k ∈ js, ∃ c d ρ, s'.outerE[k]![0]! = some c ∧ s'.outerE[k]![1]! = some d ∧
        ({P' with ves := f k} : Piece).OpenEmbedding s'.rotAdj ρ c d v := by
      intro k hk
      obtain ⟨ck, dk, ρk, hck, hdk, hek⟩ := hchild k (by simp [hk])
      refine ⟨ck, dk, ρk, ?_, ?_, hek.frame ?_⟩
      · simpa only [s', link_outerE] using hck
      · simpa only [s', link_outerE] using hdk
      · intro q hq
        apply link_rotAdj_get_ne
        · intro hh
          have : b = q := Option.some.inj hh
          exact List.disjoint_left.1 (hcursep k (by simp [hk])).1 hbcMem.1 (this ▸ hq)
        · intro hh
          have : c = q := Option.some.inj hh
          exact List.disjoint_left.1 (hpj k hk).1 hbcMem.2 (this ▸ hq)
    have hnextsep : ∀ k ∈ js, List.Disjoint P'.ves (f k) ∧
        ∀ w, HasEdge P'.es w → HasEdge ({P' with ves := f k} : Piece).es w → w = v := by
      intro k hk
      refine ⟨?_, ?_⟩
      · rw [List.disjoint_left]
        intro e he hek
        rcases List.mem_append.1 he with he | he
        · exact List.disjoint_left.1 (hcursep k (by simp [hk])).1 he hek
        · exact List.disjoint_left.1 (hpj k hk).1 he hek
      · intro w hw hwk
        obtain ⟨p, hp, hpw⟩ := hw
        simp only [P', Piece.es, List.map_append, List.mem_append] at hp
        rcases hp with hp | hp
        · exact (hcursep k (by simp [hk])).2 w ⟨p, hp, hpw⟩ hwk
        · exact (hpj k hk).2 w ⟨p, hp, hpw⟩ hwk
    have hnextbound : ∀ q, ({P' with ves := P'.ves ++ js.flatMap f} : Piece).Mem q →
        q < s'.rotAdj.size := by
      intro q hq
      change q < ((link (some b) (some c)).run s).2.rotAdj.size
      rw [link_rotAdj_size]
      apply hbound q
      simpa only [Piece.Mem, P', List.flatMap_cons, List.append_assoc] using hq
    obtain ⟨hout, hframe⟩ := ih P' s' d ρ' hρ' hnext hnextsep hpjs hnextbound
    simp only [vLoop, Option.isNone_some, Bool.false_eq_true, ↓reduceIte, hc,
      link_outerE, hd]
    constructor
    · simpa only [P', List.flatMap_cons, List.append_assoc] using hout
    · intro q hq
      rw [hframe q (by simpa only [P', List.flatMap_cons, List.append_assoc] using hq)]
      apply link_rotAdj_get_ne
      · intro hh
        have : b = q := Option.some.inj hh
        apply hq
        exact List.mem_append_left _ (this ▸ hbcMem.1)
      · intro hh
        have : c = q := Option.some.inj hh
        apply hq
        exact List.mem_append_right _ (List.mem_flatMap.2 ⟨j, by simp, this ▸ hbcMem.2⟩)

end Spqr.PlanarSpqrTree
