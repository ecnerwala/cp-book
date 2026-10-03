import Spqr.PlanarEmbedFacesSteps
import Spqr.PlanarEmbedQExec
import Spqr.PlanarEmbedVBoundary
import Spqr.Proofs.PieceInsert

/-!
# The `Q` step: certificates of the child pieces

The capped child of a block-root `Q` item carries a `Capped` certificate (its two exposed pairs,
cofacial for the witness agreeing with `rotAdj`); the lower `V` child carries an `OpenEmbedding`
or is empty.
-/

namespace Spqr.PlanarSpqrTree

variable (t : PlanarSpqrTree)

/-- The upper endpoint of edge `e` in walk orientation. -/
def qUpperV (g : Graph) (e : Nat) : Nat :=
  if t.edgeFlipped[e]! then (g.edges[e]!).2 else (g.edges[e]!).1

/-- The lower endpoint of edge `e` in walk orientation. -/
def qLowerV (g : Graph) (e : Nat) : Nat :=
  if t.edgeFlipped[e]! then (g.edges[e]!).1 else (g.edges[e]!).2

/-- The side offset of `qes`: slots `0`, `1` are the quarter-edges `4 e + qx`, `4 e + qx + 1`. -/
def qx (e : Nat) : Nat := if t.edgeFlipped[e]! then 2 else 0

theorem qx_cases (e : Nat) : t.qx e = 0 ∨ t.qx e = 2 := by
  unfold qx; cases t.edgeFlipped[e]! <;> simp

theorem qes_zero {i e : Nat} (he : t.origId[i]! = some e) : t.qes i 0 = some (4 * e + t.qx e) := by
  unfold qes qx; rw [he]; cases hf : t.edgeFlipped[e]! <;> simp [QE.mk, hf]

theorem qes_one {i e : Nat} (he : t.origId[i]! = some e) :
    t.qes i 1 = some (4 * e + (t.qx e + 1)) := by
  unfold qes qx; rw [he]; cases hf : t.edgeFlipped[e]! <;> simp [QE.mk, hf]

theorem qes_two {i e : Nat} (he : t.origId[i]! = some e) :
    t.qes i 2 = some (4 * e + (2 - t.qx e)) := by
  unfold qes qx; rw [he]; cases hf : t.edgeFlipped[e]! <;> simp [QE.mk, hf]

theorem qes_three {i e : Nat} (he : t.origId[i]! = some e) :
    t.qes i 3 = some (4 * e + (3 - t.qx e)) := by
  unfold qes qx; rw [he]; cases hf : t.edgeFlipped[e]! <;> simp [QE.mk, hf]

theorem q_lower_boundary (g : Graph) (hwf : t.toSpqrTree.WF) (hne : t.ne = g.ne)
    {i w v : Nat} (hi : i < t.size) (hw : w ∈ t.children i)
    (htw : t.toSpqrTree.type w = .V) (hv : t.origId[w]! = some v)
    (s : EmbedState) (h : t.GluedUpTo g (i + 1) s) :
    (s.outerE[w]![0]! = none ∧ t.edgesBelow w = []) ∨
    ∃ w0 w1 ρ, s.outerE[w]![0]! = some w0 ∧ s.outerE[w]![1]! = some w1 ∧
      (t.pieceBelow g w).OpenEmbedding s.rotAdj ρ w0 w1 v := by
  have hm := t.child_maximal hwf hi hw
  obtain ⟨hws, _, hiw⟩ := t.child_data hwf hi hw
  obtain ⟨ρ, hρ, ha, ho, hpair⟩ := h.piece w hm
  obtain ⟨hvert, hslots, hpres⟩ := h.outer_vertex w hws htw v hv
  have hpr := hpair 0 (by omega)
  cases h0 : s.outerE[w]![0]! with
  | none =>
    left
    refine ⟨rfl, ?_⟩
    by_contra hne'
    obtain ⟨q, k, hk⟩ := hpres (by omega) hne'
    have hkl : k < 2 := hslots k q hk
    interval_cases k
    · rw [← outer_some_iff, h0] at hk; cases hk
    · obtain ⟨a, ha0⟩ := hpr.2 q hk
      rw [← outer_some_iff, h0] at ha0; cases ha0
  | some a =>
    right
    have ha0 := (outer_some_iff s w 0 a).1 h0
    obtain ⟨b, hb1, la, lb, hla, hlb, hab⟩ := hpr.1 a ha0
    refine ⟨a, b, ρ, rfl, (outer_some_iff s w 1 b).2 hb1,
      hρ, ha, ⟨la, lb, hla, hlb, hab, ?_⟩, ?_, h.outer_dir w 0 a ha0, h.outer_dir w 1 b hb1⟩
    · rw [t.pieceBelow_loc_vert hwf hne hla]
      exact hvert a ⟨0, ha0⟩
    · intro q hq
      rw [ho q hq]
      constructor
      · rintro ⟨k, hk⟩
        have hk' := hslots k q hk
        interval_cases k
        · exact Or.inl (Option.some.inj (Option.some.inj (hk.symm.trans ha0)))
        · exact Or.inr (Option.some.inj (Option.some.inj (hk.symm.trans hb1)))
      · rintro (rfl | rfl)
        · exact ⟨0, ha0⟩
        · exact ⟨1, hb1⟩

theorem q_capped_child (g : Graph) (hwf : t.toSpqrTree.WF) (hne : t.ne = g.ne)
    (hsep : t.toSpqrTree.PieceSep g) {i c e : Nat} (hi : i < t.size)
    (ht : t.toSpqrTree.type i = .Q) (hc : c ∈ t.children i)
    (hcI : t.toSpqrTree.type c ≠ .I) (hcO : t.toSpqrTree.type c ≠ .O)
    (hcap : t.toSpqrTree.hasCap c = true) (he : t.origId[i]! = some e)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    ∃ c0 c1 c2 c3 ρ, s.outerE[c]![0]! = some c0 ∧ s.outerE[c]![1]! = some c1 ∧
      s.outerE[c]![2]! = some c2 ∧ s.outerE[c]![3]! = some c3 ∧
      (t.pieceBelow g c).Capped s.rotAdj ρ c0 c1 c2 c3 (t.qUpperV g e) (t.qLowerV g e) := by
  have hm := t.child_maximal hwf hi hc
  obtain ⟨hcs, _, hic⟩ := t.child_data hwf hi hc
  have hcne : t.toSpqrTree.capNe c = some (t.toSpqrTree.neRange c).1 := by
    simp [SpqrTree.capNe, hcap]
  obtain ⟨p, hp⟩ := hsep.cap_orig c _ hcs hcne
  have hpe : p = if t.edgeFlipped[e]! then ((g.edges[e]!).2, (g.edges[e]!).1) else g.edges[e]! :=
    hsep.q_cap_orient i c e _ p hi ht (by rwa [t.children_eq]) he hcne hp
  have hp1 : p.1 = t.qUpperV g e := by
    rw [hpe]; unfold qUpperV; cases t.edgeFlipped[e]! <;> rfl
  have hp2 : p.2 = t.qLowerV g e := by
    rw [hpe]; unfold qLowerV; cases t.edgeFlipped[e]! <;> rfl
  obtain ⟨ρ, ⟨hρ, ha, ho, hpair⟩, hf⟩ := h.cap_face c hm _ hcne
  have hnil : t.edgesBelow c ≠ [] := by
    obtain ⟨e', he'⟩ := hsep.cap_nonempty c hcs hcap hcI hcO
    have := t.mem_edgesBelow_of_edgeIn hwf he'
    intro h0; rw [h0] at this; exact List.not_mem_nil this
  obtain ⟨hvert, hpres⟩ := h.outer_cap c hcs _ hcne p hp
  have hall := hpres (by omega) hnil
  obtain ⟨c0, h0⟩ := hall 0 (by omega)
  obtain ⟨c1, h1⟩ := hall 1 (by omega)
  obtain ⟨c2, h2⟩ := hall 2 (by omega)
  obtain ⟨c3, h3⟩ := hall 3 (by omega)
  have hpr0 := hpair 0 (by omega)
  have hpr1 := hpair 1 (by omega)
  obtain ⟨b, hb, l0, l1, hl0, hl1, hg01⟩ := hpr0.1 c0 h0
  have hbc : b = c1 := Option.some.inj (Option.some.inj (hb.symm.trans h1))
  subst hbc
  obtain ⟨d, hd, l2, l3, hl2, hl3, hg23⟩ := hpr1.1 c2 h2
  have hdc : d = c3 := Option.some.inj (Option.some.inj (hd.symm.trans h3))
  subst hdc
  have hslot : ∀ k q, s.outerE[c]?.bind (fun o => o[k]?) = some (some q) → k < 4 :=
    fun k q hk => (h.outer_slots c k q hk).1
  refine ⟨c0, b, c2, d, ρ, (outer_some_iff s c 0 c0).2 h0, (outer_some_iff s c 1 b).2 h1,
    (outer_some_iff s c 2 c2).2 h2, (outer_some_iff s c 3 d).2 h3, hρ, ha, ?_, ?_, ?_,
    hf c0 c2 h0 h2, h.outer_dir c 0 c0 h0, h.outer_dir c 1 b h1, h.outer_dir c 2 c2 h2,
    h.outer_dir c 3 d h3⟩
  · intro q hq
    rw [ho q hq]
    constructor
    · rintro ⟨k, hk⟩
      have hk' := hslot k q hk
      interval_cases k
      · exact Or.inl (Option.some.inj (Option.some.inj (hk.symm.trans h0)))
      · exact Or.inr (Or.inl (Option.some.inj (Option.some.inj (hk.symm.trans h1))))
      · exact Or.inr (Or.inr (Or.inl (Option.some.inj (Option.some.inj (hk.symm.trans h2)))))
      · exact Or.inr (Or.inr (Or.inr (Option.some.inj (Option.some.inj (hk.symm.trans h3)))))
    · rintro (rfl | rfl | rfl | rfl)
      · exact ⟨0, h0⟩
      · exact ⟨1, h1⟩
      · exact ⟨2, h2⟩
      · exact ⟨3, h3⟩
  · refine ⟨l0, l1, hl0, hl1, hg01, ?_⟩
    rw [t.pieceBelow_loc_vert hwf hne hl0, hvert 0 c0 h0]
    simp [hp1]
  · refine ⟨l2, l3, hl2, hl3, hg23, ?_⟩
    rw [t.pieceBelow_loc_vert hwf hne hl2, hvert 2 c2 h2]
    simp [hp2]

/-- Assemble `GluedFaces g i` for a block-root `Q` item from an open embedding of its whole piece
at the upper endpoint, exposed as the pair `(a, b)` in slots `0`, `1`. -/
theorem q_open_glued (g : Graph) (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hrep : t.toSpqrTree.Represents g) (hsep : t.toSpqrTree.PieceSep g)
    {i e : Nat} (hi : i < t.size) (ht : t.toSpqrTree.type i = .Q) (he : t.origId[i]! = some e)
    (hcap : t.toSpqrTree.capNe i = none)
    (s s₁ : EmbedState) (h : t.GluedFaces g (i + 1) s)
    (houter : s₁.outerE = s.outerE) (hsize : s₁.rotAdj.size = s.rotAdj.size)
    (hframe : ∀ q, ¬(t.pieceBelow g i).Mem q → s₁.rotAdj[q]? = s.rotAdj[q]?)
    {a b : Nat} {ρ : RotationSystem}
    (hemb : (t.pieceBelow g i).OpenEmbedding s₁.rotAdj ρ a b (t.qUpperV g e)) :
    t.GluedFaces g i (setOuterPair i (some a) (some b) s₁) := by
  let s' := setOuterPair i (some a) (some b) s₁
  have hrot : s'.rotAdj = s₁.rotAdj := rfl
  have hlookup : ∀ j k,
      s'.outerE[j]?.bind (fun o => o[k]?) =
        if j = i then if k = 0 then some (some a) else if k = 1 then some (some b)
        else s.outerE[j]?.bind (fun o => o[k]?) else s.outerE[j]?.bind (fun o => o[k]?) := by
    intro j k
    rw [setOuterPair_lookup i _ _ s₁ (by rw [houter, h.outer_size]; exact hi)
      (by rw [houter]; exact h.outer_row_size i hi), houter]
  have hnewslot : ∀ k q, s'.outerE[i]?.bind (fun o => o[k]?) = some (some q) →
      (k = 0 ∧ a = q) ∨ (k = 1 ∧ b = q) := by
    intro k q hk
    rw [hlookup] at hk
    simp only [↓reduceIte] at hk
    by_cases hk0 : k = 0
    · simp only [hk0, ↓reduceIte, Option.some.injEq] at hk
      exact Or.inl ⟨hk0, hk⟩
    · by_cases hk1 : k = 1
      · simp only [hk1, ite_eq_right (by decide : ¬(1 : Nat) = 0), ↓reduceIte,
          Option.some.injEq] at hk
        exact Or.inr ⟨hk1, hk⟩
      · simp only [ite_eq_right hk0, ite_eq_right hk1] at hk
        exact (h.outer_unprocessed i (by omega) q ⟨k, hk⟩).elim
  have hexpose : ∀ q, s'.exposedAt i q ↔ q = a ∨ q = b := by
    intro q
    constructor
    · rintro ⟨k, hk⟩
      rcases hnewslot k q hk with ⟨_, rfl⟩ | ⟨_, rfl⟩
      · exact Or.inl rfl
      · exact Or.inr rfl
    · rintro (rfl | rfl)
      · exact ⟨0, by simp [hlookup]⟩
      · exact ⟨1, by simp [hlookup]⟩
  have hoth : ∀ j, j ≠ i → ∀ q, s'.exposedAt j q ↔ s.exposedAt j q := by
    intro j hji q
    simp only [EmbedState.exposedAt, hlookup, ite_eq_right hji]
  have himax := t.maximal_self hwf hi
  have hvb : QE.vert g.edges.toList a = some (t.qUpperV g e) ∧
      QE.vert g.edges.toList b = some (t.qUpperV g e) := by
    obtain ⟨la, lb, hla, hlb, hab, hva⟩ := hemb.boundary
    refine ⟨(t.pieceBelow_loc_vert hwf hrep.ne hla).symm.trans hva, ?_⟩
    rw [← t.pieceBelow_loc_vert hwf hrep.ne hlb]
    exact (hemb.planar.same_vertex la (by
      rw [hemb.planar.size]
      simpa only [Piece.es, List.length_map] using Piece.loc_lt hla) lb hab).symm.trans hva
  have h' : t.GluedUpTo g i s' := by
    refine ⟨⟨⟨⟨⟨⟨?_, ?_, ?_, ?_, ?_⟩, ?_, ?_⟩, ?_⟩, ?_, ?_⟩, ?_⟩, ?_⟩
    · rw [hrot, hsize, h.rot_size]
    · rw [setOuterPair_outer_size, houter, h.outer_size]
    · intro j hj q he'
      exact h.outer_unprocessed j (by omega) q ((hoth j (by omega) q).1 he')
    · intro q hq
      rw [hrot, hframe q (hq i le_rfl hi)]
      exact h.unset q (fun j hj hjs => hq j (by omega) hjs)
    · intro j hj
      by_cases hji : j = i
      · subst j
        refine ⟨ρ, hemb.planar, hemb.agrees, ?_, ?_⟩
        · intro q hq
          rw [hexpose]
          exact hemb.unset q hq
        · intro k hk
          interval_cases k
          · constructor
            · intro a' ha'
              have haa : a = a' := by simpa [hlookup] using ha'
              subst haa
              obtain ⟨la, lb, hla, hlb, hab, _⟩ := hemb.boundary
              exact ⟨b, by simp [hlookup], la, lb, hla, hlb, hab⟩
            · intro b' _
              exact ⟨a, by simp [hlookup]⟩
          · constructor
            · intro a' ha'
              rcases hnewslot 2 a' ha' with ⟨hf, _⟩ | ⟨hf, _⟩ <;> omega
            · intro b' hb'
              rcases hnewslot 3 b' hb' with ⟨hf, _⟩ | ⟨hf, _⟩ <;> omega
      · have hjsucc : t.Maximal (i + 1) j := t.maximal_succ i j (by have := hj.1; omega) hj
        obtain ⟨ρj, hρj, haj, hoj, hpj⟩ := h.piece j hjsucc
        have hf : ∀ q, (t.pieceBelow g j).Mem q → s'.rotAdj[q]? = s.rotAdj[q]? := by
          intro q hq
          rw [hrot]
          apply hframe q
          intro hqi
          exact List.disjoint_left.1 (t.maximal_pieces_disjoint hwf hsh hj himax hji) hq hqi
        refine ⟨ρj, hρj, ?_, ?_, ?_⟩
        · intro q r hq hr; exact haj q r hq (by rwa [hf q hq] at hr)
        · intro q hq; rw [hf q hq, hoth j hji]; exact hoj q hq
        · simpa only [hlookup, ite_eq_right hji] using hpj
    · intro j hj
      rw [setOuterPair_row_size, houter]
      exact h.outer_row_size j hj
    · intro j k q hk
      by_cases hji : j = i
      · subst j
        have hkl : k < 2 := by rcases hnewslot k q hk with ⟨rfl, _⟩ | ⟨rfl, _⟩ <;> omega
        exact ⟨by omega, by rw [ht]; decide, fun _ _ _ => hkl⟩
      · rw [hlookup, ite_eq_right hji] at hk
        exact h.outer_slots j k q hk
    · intro j p v q hp hpt hv he'
      by_cases hji : j = i
      · subst j
        have hu := hsep.q_upper i p e hi ht hp hpt he
        have hvu : v = t.qUpperV g e := Option.some.inj (hv.symm.trans hu)
        subst hvu
        rcases (hexpose q).1 he' with rfl | rfl
        · exact hvb.1
        · exact hvb.2
      · exact h.outer_at_vertex j p v q hp hpt hv ((hoth j hji q).1 he')
    · intro j k q hk
      by_cases hji : j = i
      · subst j
        rcases hnewslot k q hk with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩
        · exact hemb.left_dir
        · exact hemb.right_dir
      · rw [hlookup, ite_eq_right hji] at hk
        exact h.outer_dir j k q hk
    · intro j hj hjs p hp hpt hne
      by_cases hji : j = i
      · subst j
        exact ⟨a, (hexpose a).2 (Or.inl rfl)⟩
      · obtain ⟨q, hq⟩ := h.outer_present j (by omega) hjs p hp hpt hne
        exact ⟨q, (hoth j hji q).2 hq⟩
    · intro j hjs hjt v hv
      by_cases hji : j = i
      · subst j; rw [ht] at hjt; cases hjt
      · obtain ⟨hvert, hslots, hpres⟩ := h.outer_vertex j hjs hjt v hv
        refine ⟨?_, ?_, ?_⟩
        · intro q hq; exact hvert q ((hoth j hji q).1 hq)
        · intro k q hk
          rw [hlookup, ite_eq_right hji] at hk
          exact hslots k q hk
        · intro hj hne
          obtain ⟨q, hq⟩ := hpres (by omega) hne
          exact ⟨q, (hoth j hji q).2 hq⟩
    · intro j hjs ne hc p hp
      have hji : j ≠ i := by
        intro he'; subst j; rw [hcap] at hc; cases hc
      obtain ⟨hvert, hpres⟩ := h.outer_cap j hjs ne hc p hp
      simp only [hlookup, ite_eq_right hji]
      refine ⟨hvert, ?_⟩
      intro hj hne
      exact hpres (by omega) hne
  refine t.gluedFaces_of_frame h h' ?_ ?_ ?_
  · intro ne hne; rw [hcap] at hne; cases hne
  · intro j hm hji q hq
    rw [setOuterPair_rotAdj]
    apply hframe q
    intro hqi
    exact t.maximal_pieces_disjoint hwf hsh hm himax hji hq hqi
  · intro j hji
    rw [setOuterPair_outer_ne _ _ _ _ j hji, houter]

end Spqr.PlanarSpqrTree
