import Spqr.PlanarEmbedFacesSteps
import Spqr.PlanarEmbedQExec
import Spqr.PlanarEmbedVBoundary
import Spqr.Proofs.PieceInsert
import Spqr.Proofs.Dfs

/-!
# The `Q` step: certificates of the child pieces

The capped child of a block-root `Q` item carries a `Capped` certificate (its two exposed pairs,
cofacial for the witness agreeing with `rotAdj`); the lower `V` child carries an `OpenEmbedding`
or is empty.
-/

namespace Spqr.PlanarSpqrTree

variable (t : PlanarSpqrTree)

open EmbedM

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

/-- The edge of `Q` item `i` lies below no item strictly after `i`. -/
theorem q_edge_not_below (hwf : t.toSpqrTree.WF) {i e j : Nat} (hidx : t.edgeIndex[e]! = some i)
    (hij : i < j) : e ∉ t.edgesBelow j := by
  intro hm
  obtain ⟨j', hjj', _, _, _, _, hidx'⟩ := t.mem_edgesBelow_data hwf hm
  have : j' = i := Option.some.inj (hidx'.symm.trans hidx)
  omega

/-- Before the `Q` step, the four quarter-edges of its edge are unset. -/
theorem q_fresh (g : Graph) (hwf : t.toSpqrTree.WF) {i e : Nat}
    (hidx : t.edgeIndex[e]! = some i) (s : EmbedState) (h : t.GluedUpTo g (i + 1) s)
    {k : Nat} (hk : k < 4) : s.rotAdj[4 * e + k]? = some none := by
  have he : e < t.ne := t.edgeIndex_some_lt hwf hidx
  have hlt : 4 * e + k < s.rotAdj.size := by rw [h.rot_size]; omega
  rcases h.unset (4 * e + k) (fun j hj _ => by
      have : QE.edge (4 * e + k) = e := by unfold QE.edge; omega
      rw [this]; exact t.q_edge_not_below hwf hidx (by omega)) with h0 | h0
  · rw [Array.getElem?_eq_none_iff] at h0; omega
  · exact h0

/-- The four quarter-edges of the edge of a one-edge piece. -/
theorem mem_single {P : Piece} {e k : Nat} (hk : k < 4) :
    ({P with ves := [e]} : Piece).Mem (4 * e + k) := by
  simp only [Piece.Mem, QE.edge, List.mem_singleton]; omega

/-- A `Q` item with an `O` child: the one-edge loop piece with witness `loop`, after the
executable `link` of slots `1` and `2`. -/
theorem loop_open (P : Piece) {s : EmbedState} {u e x : Nat} (hx : x = 0 ∨ x = 2)
    (hends : P.ends e = (u, u)) (hu : u < P.nVerts)
    (hA : ∀ k, k < 4 → s.rotAdj[4 * e + k]? = some none) (hbound : 4 * e + 3 < s.rotAdj.size) :
    ({P with ves := [e]} : Piece).OpenEmbedding
      ((link (some (4 * e + (x + 1))) (some (4 * e + (2 - x)))).run s).2.rotAdj
      RotationSystem.loop (4 * e + x) (4 * e + (3 - x)) u := by
  have hloc : ∀ k, k < 4 → ({P with ves := [e]} : Piece).loc (4 * e + k) = some k :=
    fun k hk => Piece.loc_cons_self {P with ves := []} e k hk
  have hmem : ∀ q, ({P with ves := [e]} : Piece).Mem q → q = 4 * e + q % 4 ∧ q % 4 < 4 := by
    intro q hq
    simp only [Piece.Mem, QE.edge, List.mem_singleton] at hq
    exact ⟨by omega, Nat.mod_lt _ (by decide)⟩
  have hes : ({P with ves := [e]} : Piece).es = [P.ends e] := by simp [Piece.es]
  have hget := link_rotAdj_get (4 * e + (x + 1)) (4 * e + (2 - x)) s (by omega) (by omega)
    (by omega)
  refine ⟨?_, ?_, ?_, ?_, by omega, by omega⟩
  · rw [hes, hends]; exact isPlanarEmbedding_loop hu
  · intro q r hq hr
    obtain ⟨hq', hk⟩ := hmem q hq
    rw [hget] at hr
    split_ifs at hr with h1 h2
    · subst h1
      refine ⟨x + 1, 2 - x, hloc _ (by omega), ?_, ?_⟩
      · rw [← Option.some.inj (Option.some.inj hr)]; exact hloc _ (by omega)
      · rcases hx with rfl | rfl <;> rfl
    · subst h2
      refine ⟨2 - x, x + 1, hloc _ (by omega), ?_, ?_⟩
      · rw [← Option.some.inj (Option.some.inj hr)]; exact hloc _ (by omega)
      · rcases hx with rfl | rfl <;> rfl
    · rw [hq', hA _ hk] at hr; cases hr
  · refine ⟨x, 3 - x, hloc x (by omega), hloc (3 - x) (by omega), ?_, ?_⟩
    · rcases hx with rfl | rfl <;> rfl
    · rw [hes, vert_cons_lt _ _ (by omega), hends]; rcases hx with rfl | rfl <;> simp
  · intro q hq
    obtain ⟨hq', hk⟩ := hmem q hq
    rw [hget]
    split_ifs with h1 h2
    · simp only [Option.some.injEq, reduceCtorEq, false_iff]; omega
    · simp only [Option.some.injEq, reduceCtorEq, false_iff]; omega
    · rw [hq', hA _ hk]; simp only [true_iff]; omega

/-- The `Q` step for the root of a loop block (one `O` child). -/
theorem q_loop_glued (g : Graph) (hg : g.WF) (hwf : t.toSpqrTree.WF)
    (hsh : t.toSpqrTree.ChildShape) (hrep : t.toSpqrTree.Represents g)
    (hsep : t.toSpqrTree.PieceSep g) {i c e : Nat} (hi : i < t.size)
    (ht : t.toSpqrTree.type i = .Q) (he : t.origId[i]! = some e)
    (hidx : t.edgeIndex[e]! = some i) (hc : t.children i = [c])
    (hcO : t.toSpqrTree.type c = .O) (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    t.GluedFaces g i
      (setOuterPair i (t.qes i 0) (t.qes i 3) ((link (t.qes i 1) (t.qes i 2)).run s).2) := by
  have hty : t.types[i]! = .Q := by rw [← t.type_eq_of_lt i hi]; exact ht
  obtain ⟨hcs, _, _⟩ := t.child_data hwf hi (by rw [hc]; exact List.mem_singleton_self c)
  have hcO' : t.types[c]! = .O := by rw [← t.type_eq_of_lt c hcs]; exact hcO
  have hes : t.edgesBelow i = [e] := by
    rw [t.edgesBelow_eq hwf hsh i hi, hty, he, hc, List.flatMap_cons, List.flatMap_nil,
      t.edgesBelow_leaf c (t.subtreeEnd_leaf hwf hsh c hcs (Or.inl hcO')) (by rw [hcO']; decide)]
    rfl
  have hcap : t.toSpqrTree.capNe i = none := by
    simp [SpqrTree.capNe, SpqrTree.hasCap, ht, t.children_eq, hc, hcO]
  have hloop : (g.edges[e]!).1 = (g.edges[e]!).2 :=
    (hsep.q_loop i e hi ht he).2
      ⟨c, by rw [t.children_eq, hc]; exact List.mem_singleton_self c, hcO⟩
  have hu : t.qUpperV g e = (g.edges[e]!).1 := by
    unfold qUpperV; split
    · exact hloop.symm
    · rfl
  have hene : e < g.ne := by rw [← hrep.ne]; exact t.edgeIndex_some_lt hwf hidx
  have hult : (g.edges[e]!).1 < g.nv := (hg.getElem! hene).1
  have hx := t.qx_cases e
  rw [t.qes_zero he, t.qes_one he, t.qes_two he, t.qes_three he]
  have hA : ∀ k, k < 4 → s.rotAdj[4 * e + k]? = some none :=
    fun k hk => t.q_fresh g hwf hidx s h.toGluedUpTo hk
  have hbound : 4 * e + 3 < s.rotAdj.size := by
    rw [h.rot_size]; have := t.edgeIndex_some_lt hwf hidx; omega
  have hopen := loop_open (t.pieceBelow g i) hx (u := (g.edges[e]!).1)
    (Prod.ext rfl hloop.symm) hult hA hbound
  have hPi : ({t.pieceBelow g i with ves := [e]} : Piece) = t.pieceBelow g i := by
    simp [pieceBelow, hes]
  rw [hPi, ← hu] at hopen
  have hmem : ∀ k, k < 4 → (t.pieceBelow g i).Mem (4 * e + k) := fun k hk => by
    simp only [Piece.Mem, pieceBelow, hes, QE.edge, List.mem_singleton]; omega
  refine t.q_open_glued g hwf hsh hrep hsep hi ht he hcap s _ h (link_outerE _ _ _)
    (link_rotAdj_size _ _ _) ?_ hopen
  intro q hq
  rw [link_rotAdj_get _ _ _ (by omega) (by omega) (by omega),
    ite_eq_right (fun h' => hq (by rw [h']; exact hmem _ (by omega))),
    ite_eq_right (fun h' => hq (by rw [h']; exact hmem _ (by omega)))]

theorem setOuter4_rotAdj (i : Nat) (s : EmbedState) : (t.setOuter4 i s).rotAdj = s.rotAdj := rfl

theorem setOuter4_outer_size (i : Nat) (s : EmbedState) :
    (t.setOuter4 i s).outerE.size = s.outerE.size := by
  simp only [setOuter4, setOuter_outerE_size]

theorem setOuter4_outer_ne (i : Nat) (s : EmbedState) (j : Nat) (hji : j ≠ i) :
    (t.setOuter4 i s).outerE[j]? = s.outerE[j]? := by
  simp only [setOuter4, setOuter_outerE_ne i _ _ _ j hji]

theorem setOuter4_outer_eq (i : Nat) (s : EmbedState) (hi : i < s.outerE.size) :
    (t.setOuter4 i s).outerE[i]? =
      some ((((s.outerE[i]!.set! 0 (t.qes i 0)).set! 1 (t.qes i 1)).set! 2 (t.qes i 2)).set! 3
        (t.qes i 3)) := by
  simp [setOuter4, setOuter_run, Array.getElem_modify_self, hi]

theorem setOuter4_lookup (i : Nat) (s : EmbedState) (hi : i < s.outerE.size)
    (hs : s.outerE[i]!.size = 4) (j k : Nat) :
    (t.setOuter4 i s).outerE[j]?.bind (fun o => o[k]?) =
      if j = i then (if k < 4 then some (t.qes i k) else none)
      else s.outerE[j]?.bind (fun o => o[k]?) := by
  by_cases hji : j = i
  · subst j
    rw [t.setOuter4_outer_eq _ _ hi]
    simp only [Option.bind_some, ↓reduceIte, Array.set!, Array.getElem?_setIfInBounds,
      Array.size_setIfInBounds, hs]
    by_cases hk : k < 4
    · by_cases hk0 : k = 0 <;> by_cases hk1 : k = 1 <;> by_cases hk2 : k = 2 <;>
        by_cases hk3 : k = 3 <;> simp_all [eq_comm] <;> omega
    · rw [Array.getElem?_eq_none_iff.2 (by omega : s.outerE[i]!.size ≤ k)]
      simp [show 3 ≠ k by omega, show 2 ≠ k by omega, show 1 ≠ k by omega,
        show 0 ≠ k by omega, hk]
  · rw [t.setOuter4_outer_ne _ _ j hji, ite_eq_right hji]

theorem setOuter4_row_size (i : Nat) (s : EmbedState) (j : Nat) :
    (t.setOuter4 i s).outerE[j]!.size = s.outerE[j]!.size := by
  by_cases hji : j = i
  · subst j
    by_cases hi : i < s.outerE.size
    · simp only [getElem!_def, t.setOuter4_outer_eq _ _ hi]
      simp
    · simp [getElem!_def, Array.getElem?_eq_none (by omega : s.outerE.size ≤ i),
        Array.getElem?_eq_none (by rw [setOuter4_outer_size]; omega :
          (t.setOuter4 i s).outerE.size ≤ i)]
  · simp only [getElem!_def, t.setOuter4_outer_ne _ _ j hji]

/-- `link` leaves every other entry unchanged. -/
theorem link_rotAdj_get_other (a b : Nat) (s : EmbedState) {q : Nat} (hqa : q ≠ a) (hqb : q ≠ b) :
    ((link (some a) (some b)).run s).2.rotAdj[q]? = s.rotAdj[q]? := by
  rw [link_some_run]
  simp only [Array.set!, Array.getElem?_setIfInBounds, Array.size_setIfInBounds]
  simp [Ne.symm hqa, Ne.symm hqb]

/-- The `Q` step for a block root: the capped child `c` (an `I` leaf or a capped node/`Q`) and
the lower `V` item `w`. -/
theorem q_block_glued (g : Graph) (hg : g.WF) (hwf : t.toSpqrTree.WF)
    (hsh : t.toSpqrTree.ChildShape) (hrep : t.toSpqrTree.Represents g)
    (hsep : t.toSpqrTree.PieceSep g) {i c w e : Nat} (hi : i < t.size)
    (ht : t.toSpqrTree.type i = .Q) (he : t.origId[i]! = some e)
    (hidx : t.edgeIndex[e]! = some i) (hc : t.children i = [c, w])
    (hcV : t.toSpqrTree.type c ≠ .V) (hcO : t.toSpqrTree.type c ≠ .O)
    (hcap : t.toSpqrTree.hasCap c = true) (hwV : t.toSpqrTree.type w = .V)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    t.GluedFaces g i
      (setOuterPair i (t.qes i 0) (t.qUpper c (t.qes i 1) (t.qes i 2) s).1.1
        (qLower w (t.qes i 3) (t.qUpper c (t.qes i 1) (t.qes i 2) s).1.2
          (t.qUpper c (t.qes i 1) (t.qes i 2) s).2)) := by
  have hty : t.types[i]! = .Q := by rw [← t.type_eq_of_lt i hi]; exact ht
  have hcm : c ∈ t.children i := by rw [hc]; simp
  have hwm : w ∈ t.children i := by rw [hc]; simp
  obtain ⟨hcs, _, hic⟩ := t.child_data hwf hi hcm
  obtain ⟨hws, _, hiw⟩ := t.child_data hwf hi hwm
  have hcw : c ≠ w := fun h' => hcV (h' ▸ hwV)
  have hx := t.qx_cases e
  have hene : e < g.ne := by rw [← hrep.ne]; exact t.edgeIndex_some_lt hwf hidx
  have hult : t.qUpperV g e < g.nv := by
    unfold qUpperV; split
    · exact (hg.getElem! hene).2
    · exact (hg.getElem! hene).1
  have hvlt : t.qLowerV g e < g.nv := by
    unfold qLowerV; split
    · exact (hg.getElem! hene).1
    · exact (hg.getElem! hene).2
  have hno : ¬ (g.edges[e]!).1 = (g.edges[e]!).2 := by
    intro hl
    obtain ⟨o, ho, hoO⟩ := (hsep.q_loop i e hi ht he).1 hl
    rw [t.children_eq, hc] at ho
    simp only [List.mem_cons, List.not_mem_nil, or_false] at ho
    rcases ho with rfl | rfl
    · exact hcO hoO
    · rw [hwV] at hoO; cases hoO
  have huv : t.qUpperV g e ≠ t.qLowerV g e := by
    unfold qUpperV qLowerV; split
    · exact fun h' => hno h'.symm
    · exact hno
  have hends : (t.pieceBelow g c).ends e =
      if t.qx e = 0 then (t.qUpperV g e, t.qLowerV g e) else (t.qLowerV g e, t.qUpperV g e) := by
    unfold qx qUpperV qLowerV
    show g.edges[e]! = _
    cases t.edgeFlipped[e]! <;> simp
  have hes : t.edgesBelow i = e :: (t.edgesBelow c ++ t.edgesBelow w) := by
    rw [t.edgesBelow_eq hwf hsh i hi, hty, he, hc]
    simp
  have hcapi : t.toSpqrTree.capNe i = none := by
    simp [SpqrTree.capNe, SpqrTree.hasCap, ht, t.children_eq, hc, hcV]
  have henc : e ∉ t.edgesBelow c := t.q_edge_not_below hwf hidx hic
  have henw : e ∉ t.edgesBelow w := t.q_edge_not_below hwf hidx hiw
  have hdis : List.Disjoint (t.edgesBelow c) (t.edgesBelow w) :=
    t.maximal_pieces_disjoint hwf hsh (t.child_maximal hwf hi hcm) (t.child_maximal hwf hi hwm)
      hcw
  have hA : ∀ k, k < 4 → s.rotAdj[4 * e + k]? = some none :=
    fun k hk => t.q_fresh g hwf hidx s h.toGluedUpTo hk
  set Pc : Piece := {t.pieceBelow g c with ves := e :: t.edgesBelow c} with hPc
  have hPi : t.pieceBelow g i = {Pc with ves := Pc.ves ++ t.edgesBelow w} := by
    simp [pieceBelow, hes, hPc]
  have hPw : ({Pc with ves := t.edgesBelow w} : Piece) = t.pieceBelow g w := rfl
  have hsubc : ∀ q, Pc.Mem q → (t.pieceBelow g i).Mem q := by
    intro q hq
    simp only [Piece.Mem, pieceBelow, hes, hPc, List.mem_cons, List.mem_append] at hq ⊢
    tauto
  have hsubw : ∀ q, (t.pieceBelow g w).Mem q → (t.pieceBelow g i).Mem q := by
    intro q hq
    simp only [Piece.Mem, pieceBelow, hes, List.mem_cons, List.mem_append] at hq ⊢
    tauto
  have hPcmem : ∀ q, (t.pieceBelow g c).Mem q → Pc.Mem q :=
    fun q hq => List.mem_cons_of_mem e hq
  have hPce : ∀ k, k < 4 → Pc.Mem (4 * e + k) := by
    intro k hk
    have : QE.edge (4 * e + k) = e := by unfold QE.edge; omega
    show QE.edge (4 * e + k) ∈ e :: t.edgesBelow c
    rw [this]; simp
  have hbound : ∀ q, (t.pieceBelow g i).Mem q → q < s.rotAdj.size := by
    intro q hq; rw [h.rot_size]; exact t.mem_pieceBelow_bound hwf g hq
  have hsepv : ∀ y, HasEdge Pc.es y → HasEdge (t.pieceBelow g w).es y → y = t.qLowerV g e := by
    intro y hy hyw
    have hyw' := t.touches_of_hasEdge hwf g hrep.ne hyw
    have hy' : t.toSpqrTree.Touches g c y ∨ SpqrTree.Graph.Incident g y e := by
      obtain ⟨p, hp, hpy⟩ := hy
      simp only [hPc, Piece.es, List.map_cons, List.mem_cons] at hp
      rcases hp with rfl | hp
      · exact Or.inr ⟨hene, hpy⟩
      · exact Or.inl (t.touches_of_hasEdge hwf g hrep.ne ⟨p, hp, hpy⟩)
    exact hsep.q_lower_attach i c w e hi ht (by rwa [t.children_eq]) he y hy' hyw'
  have hwnmem : ∀ q, (t.pieceBelow g w).Mem q → ¬ Pc.Mem q := by
    intro q hq hq'
    simp only [Piece.Mem, pieceBelow, hPc, List.mem_cons] at hq hq'
    rcases hq' with hq' | hq'
    · exact henw (hq' ▸ hq)
    · exact List.disjoint_left.1 hdis hq' hq
  have hupper : ∃ b a' ρ₂ s₂, t.qUpper c (t.qes i 1) (t.qes i 2) s = ((some b, some a'), s₂) ∧
      Pc.Open2 s₂.rotAdj ρ₂ (4 * e + t.qx e) b a' (4 * e + (3 - t.qx e)) (t.qUpperV g e)
        (t.qLowerV g e) ∧
      s₂.outerE = s.outerE ∧ s₂.rotAdj.size = s.rotAdj.size ∧
      ∀ q, ¬ Pc.Mem q → s₂.rotAdj[q]? = s.rotAdj[q]? := by
    by_cases hcI : t.toSpqrTree.type c = .I
    · have hcI' : t.types[c]! = .I := by rw [← t.type_eq_of_lt c hcs]; exact hcI
      have hcnil : t.edgesBelow c = [] :=
        t.edgesBelow_leaf c (t.subtreeEnd_leaf hwf hsh c hcs (Or.inr hcI')) (by rw [hcI']; decide)
      refine ⟨4 * e + (t.qx e + 1), 4 * e + (2 - t.qx e), RotationSystem.single, s, ?_, ?_, rfl,
        rfl, fun _ _ => rfl⟩
      · simp [qUpper, hcI', t.qes_one he, t.qes_two he]
      · have := Piece.single_open2 (t.pieceBelow g c) hx hends huv hult hvlt hA
        rw [hPc, hcnil]; exact this
    · obtain ⟨c0, c1, c2, c3, ρ, h0, h1, h2, h3, hcap'⟩ :=
        t.q_capped_child g hwf hrep.ne hsep hi ht hcm hcI hcO hcap he s h
      have hcI' : (t.types[c]! != .I) = true := by
        rw [← t.type_eq_of_lt c hcs]; simpa using hcI
      obtain ⟨ρ', hρ'⟩ := hcap'.insert hx henc hends huv hA (fun q hq => hbound q (hsubc q hq))
      obtain ⟨l0, _, hl0, _⟩ := hcap'.pair0
      obtain ⟨_, l3, _, hl3, _⟩ := hcap'.pair2
      have hm0 : Pc.Mem c0 := hPcmem c0 (Piece.mem_of_loc hl0)
      have hm3 : Pc.Mem c3 := hPcmem c3 (Piece.mem_of_loc hl3)
      refine ⟨c1, c2, ρ', _, ?_, hρ', ?_, ?_, ?_⟩
      · simp only [qUpper, hcI', ↓reduceIte, h0, h1, h2, h3, t.qes_one he, t.qes_two he]
      · simp only [link_outerE]
      · simp only [link_rotAdj_size]
      · intro q hq
        rw [link_rotAdj_get_other _ _ _ (fun h' => hq (by subst h'; exact hPce _ (by omega))) (fun h' => hq (by subst h'; exact hm3)),
          link_rotAdj_get_other _ _ _ (fun h' => hq (by subst h'; exact hPce _ (by omega))) (fun h' => hq (by subst h'; exact hm0))]
  obtain ⟨b, a', ρ₂, s₂, hqU, hO2, hout2, hsz2, hfr2⟩ := hupper
  rw [hqU, t.qes_zero he, t.qes_three he]
  simp only
  have hbound2 : ∀ q, (t.pieceBelow g i).Mem q → q < s₂.rotAdj.size := by
    rw [hsz2]; exact hbound
  obtain ⟨_, lb, _, hlb, _⟩ := hO2.pair
  obtain ⟨la', _, hla', _⟩ := hO2.pair'
  have hmemb : Pc.Mem b := Piece.mem_of_loc hlb
  have hmema' : Pc.Mem a' := Piece.mem_of_loc hla'
  have hmem3 : Pc.Mem (4 * e + (3 - t.qx e)) := hPce _ (by omega)
  have hv' : t.origId[w]! = some (t.qLowerV g e) :=
    hsep.q_lower i c w e hi ht (by rwa [t.children_eq]) he
  rcases t.q_lower_boundary g hwf hrep.ne hi hwm hwV hv' s h.toGluedUpTo with
    ⟨hw0, hwnil⟩ | ⟨w0, w1, ρw, hw0, hw1, hwopen⟩
  · have hqL : qLower w (some (4 * e + (3 - t.qx e))) (some a') s₂ =
        ((link (some (4 * e + (3 - t.qx e))) (some a')).run s₂).2 := by
      simp [qLower, hout2, hw0]
    rw [hqL]
    have hPi' : t.pieceBelow g i = Pc := by rw [hPi, hwnil, List.append_nil]
    have hopen := hO2.close huv (fun q hq => hbound2 q (hsubc q hq))
    refine t.q_open_glued g hwf hsh hrep hsep hi ht he hcapi s _ h ?_ ?_ ?_ (by rw [hPi']; exact hopen)
    · rw [link_outerE, hout2]
    · rw [link_rotAdj_size, hsz2]
    · intro q hq
      rw [link_rotAdj_get_other _ _ _ (fun h' => hq (by subst h'; exact hsubc _ hmem3))
        (fun h' => hq (by subst h'; exact hsubc _ hmema')), hfr2 q (fun h' => hq (hsubc q h'))]
  · have hqL : qLower w (some (4 * e + (3 - t.qx e))) (some a') s₂ =
        ((link (some a') (some w1)).run
          ((link (some (4 * e + (3 - t.qx e))) (some w0)).run s₂).2).2 := by
      simp [qLower, hout2, hw0, hw1]
    rw [hqL]
    obtain ⟨lw0, lw1, hlw0, hlw1, _, _⟩ := hwopen.boundary
    have hmw0 : (t.pieceBelow g w).Mem w0 := Piece.mem_of_loc hlw0
    have hmw1 : (t.pieceBelow g w).Mem w1 := Piece.mem_of_loc hlw1
    have hwopen' : ({Pc with ves := t.edgesBelow w} : Piece).OpenEmbedding s₂.rotAdj ρw w0 w1
        (t.qLowerV g e) := by
      rw [hPw]; exact hwopen.frame (fun q hq => hfr2 q (hwnmem q hq))
    obtain ⟨ρ', hρ'⟩ := hO2.splice huv hwopen'
      (by rw [hPc]; exact List.disjoint_cons_left.2 ⟨henw, hdis⟩) hsepv
      (fun q hq => hbound2 q (by rw [hPi]; exact hq))
    refine t.q_open_glued g hwf hsh hrep hsep hi ht he hcapi s _ h ?_ ?_ ?_ (by rw [hPi]; exact hρ')
    · rw [link_outerE, link_outerE, hout2]
    · rw [link_rotAdj_size, link_rotAdj_size, hsz2]
    · intro q hq
      rw [link_rotAdj_get_other _ _ _ (fun h' => hq (by subst h'; exact hsubc _ hmema'))
        (fun h' => hq (by subst h'; exact hsubw _ hmw1)),
        link_rotAdj_get_other _ _ _ (fun h' => hq (by subst h'; exact hsubc _ hmem3))
        (fun h' => hq (by subst h'; exact hsubw _ hmw0)), hfr2 q (fun h' => hq (hsubc q h'))]

theorem vert_edge (g : Graph) {e : Nat} (he : e < g.ne) {k : Nat} (hk : k < 4) :
    QE.vert g.edges.toList (4 * e + k) =
      some (if k < 2 then (g.edges[e]!).1 else (g.edges[e]!).2) := by
  have he' : e < g.edges.size := he
  have hside : (4 * e + k) / 2 % 2 = k / 2 := by omega
  simp only [QE.vert, QE.edge, QE.side, Array.getElem?_toList, show (4 * e + k) / 4 = e by omega,
    hside, Array.getElem?_eq_getElem he', Option.map_some, getElem!_pos g.edges e he']
  congr 1
  interval_cases k <;> simp

/-- The `Q` step for a capped `Q` leaf: all four quarter-edges become exposed. -/
theorem q_leaf_glued (g : Graph) (hg : g.WF) (hwf : t.toSpqrTree.WF)
    (hsh : t.toSpqrTree.ChildShape) (hrep : t.toSpqrTree.Represents g)
    (hsep : t.toSpqrTree.PieceSep g) {i e : Nat} (hi : i < t.size)
    (ht : t.toSpqrTree.type i = .Q) (he : t.origId[i]! = some e)
    (hidx : t.edgeIndex[e]! = some i) (hc : t.children i = [])
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    t.GluedFaces g i (t.setOuter4 i s) := by
  have hty : t.types[i]! = .Q := by rw [← t.type_eq_of_lt i hi]; exact ht
  have hx := t.qx_cases e
  have hene : e < g.ne := by rw [← hrep.ne]; exact t.edgeIndex_some_lt hwf hidx
  have hult : t.qUpperV g e < g.nv := by
    unfold qUpperV; split
    · exact (hg.getElem! hene).2
    · exact (hg.getElem! hene).1
  have hvlt : t.qLowerV g e < g.nv := by
    unfold qLowerV; split
    · exact (hg.getElem! hene).1
    · exact (hg.getElem! hene).2
  have hno : ¬ (g.edges[e]!).1 = (g.edges[e]!).2 := by
    intro hl
    obtain ⟨o, ho, _⟩ := (hsep.q_loop i e hi ht he).1 hl
    rw [t.children_eq, hc] at ho
    cases ho
  have huv : t.qUpperV g e ≠ t.qLowerV g e := by
    unfold qUpperV qLowerV; split
    · exact fun h' => hno h'.symm
    · exact hno
  have hends : (t.pieceBelow g i).ends e =
      if t.qx e = 0 then (t.qUpperV g e, t.qLowerV g e) else (t.qLowerV g e, t.qUpperV g e) := by
    unfold qx qUpperV qLowerV
    show g.edges[e]! = _
    cases t.edgeFlipped[e]! <;> simp
  have hes : t.edgesBelow i = [e] := by
    rw [t.edgesBelow_eq hwf hsh i hi, hty, he, hc]
    simp
  have hA : ∀ k, k < 4 → s.rotAdj[4 * e + k]? = some none :=
    fun k hk => t.q_fresh g hwf hidx s h.toGluedUpTo hk
  have hnoparent : ∀ p, t.toSpqrTree.parent i = some p →
      t.toSpqrTree.type p ≠ .F ∧ t.toSpqrTree.type p ≠ .V :=
    fun p hp => hsep.q_leaf_parent i p hi ht (by rw [t.children_eq]; exact hc) hp
  set P₁ : Piece := {t.pieceBelow g i with ves := [e]} with hP₁
  have hPi : t.pieceBelow g i = P₁ := by simp [pieceBelow, hes, hP₁]
  have ho := Piece.single_open2 (t.pieceBelow g i) hx hends huv hult hvlt hA
  have hloc : ∀ k, k < 4 → P₁.loc (4 * e + k) = some k :=
    fun k hk => Piece.loc_cons_self {t.pieceBelow g i with ves := []} e k hk
  have hmem : ∀ q, P₁.Mem q → q = 4 * e + q % 4 ∧ q % 4 < 4 := by
    intro q hq
    simp only [hP₁, Piece.Mem, QE.edge, List.mem_singleton] at hq
    exact ⟨by omega, Nat.mod_lt _ (by decide)⟩
  have hvq : ∀ k, k < 4 → QE.vert g.edges.toList (4 * e + k) =
      some (if k < 2 then (g.edges[e]!).1 else (g.edges[e]!).2) :=
    fun k hk => vert_edge g hene hk
  set s' := t.setOuter4 i s with hs'
  have hrot : s'.rotAdj = s.rotAdj := rfl
  have hlookup : ∀ j k, s'.outerE[j]?.bind (fun o => o[k]?) =
      if j = i then (if k < 4 then some (t.qes i k) else none)
      else s.outerE[j]?.bind (fun o => o[k]?) :=
    fun j k => t.setOuter4_lookup i s (by rw [h.outer_size]; exact hi) (h.outer_row_size i hi) j k
  have hslot : ∀ k q, s'.outerE[i]?.bind (fun o => o[k]?) = some (some q) →
      (k = 0 ∧ q = 4 * e + t.qx e) ∨ (k = 1 ∧ q = 4 * e + (t.qx e + 1)) ∨
      (k = 2 ∧ q = 4 * e + (2 - t.qx e)) ∨ (k = 3 ∧ q = 4 * e + (3 - t.qx e)) := by
    intro k q hk
    rw [hlookup, ite_eq_left rfl] at hk
    by_cases hk4 : k < 4
    · rw [ite_eq_left hk4] at hk
      interval_cases k
      · rw [t.qes_zero he] at hk; exact Or.inl ⟨rfl, (Option.some.inj (Option.some.inj hk)).symm⟩
      · rw [t.qes_one he] at hk
        exact Or.inr (Or.inl ⟨rfl, (Option.some.inj (Option.some.inj hk)).symm⟩)
      · rw [t.qes_two he] at hk
        exact Or.inr (Or.inr (Or.inl ⟨rfl, (Option.some.inj (Option.some.inj hk)).symm⟩))
      · rw [t.qes_three he] at hk
        exact Or.inr (Or.inr (Or.inr ⟨rfl, (Option.some.inj (Option.some.inj hk)).symm⟩))
    · rw [ite_eq_right hk4] at hk; cases hk
  have hs0 : s'.outerE[i]?.bind (fun o => o[0]?) = some (some (4 * e + t.qx e)) := by
    rw [hlookup, ite_eq_left rfl, ite_eq_left (by decide), t.qes_zero he]
  have hs1 : s'.outerE[i]?.bind (fun o => o[1]?) = some (some (4 * e + (t.qx e + 1))) := by
    rw [hlookup, ite_eq_left rfl, ite_eq_left (by decide), t.qes_one he]
  have hs2 : s'.outerE[i]?.bind (fun o => o[2]?) = some (some (4 * e + (2 - t.qx e))) := by
    rw [hlookup, ite_eq_left rfl, ite_eq_left (by decide), t.qes_two he]
  have hs3 : s'.outerE[i]?.bind (fun o => o[3]?) = some (some (4 * e + (3 - t.qx e))) := by
    rw [hlookup, ite_eq_left rfl, ite_eq_left (by decide), t.qes_three he]
  have hexpose : ∀ q, s'.exposedAt i q ↔
      q = 4 * e + t.qx e ∨ q = 4 * e + (t.qx e + 1) ∨ q = 4 * e + (2 - t.qx e) ∨
        q = 4 * e + (3 - t.qx e) := by
    intro q
    constructor
    · rintro ⟨k, hk⟩
      rcases hslot k q hk with ⟨_, rfl⟩ | ⟨_, rfl⟩ | ⟨_, rfl⟩ | ⟨_, rfl⟩ <;> simp
    · rintro (rfl | rfl | rfl | rfl)
      · exact ⟨0, hs0⟩
      · exact ⟨1, hs1⟩
      · exact ⟨2, hs2⟩
      · exact ⟨3, hs3⟩
  have hoth : ∀ j, j ≠ i → ∀ q, s'.exposedAt j q ↔ s.exposedAt j q := by
    intro j hji q
    simp only [EmbedState.exposedAt, hlookup, ite_eq_right hji]
  have himax := t.maximal_self hwf hi
  have hw : t.PieceWitness g s' i RotationSystem.single := by
    rw [PieceWitness, hPi]
    refine ⟨ho.planar, ho.agrees, ?_, ?_⟩
    · intro q hq
      rw [hrot, ho.unset q hq, hexpose]
    · intro k hk
      interval_cases k
      · refine ⟨fun a ha => ?_, fun b _ => ⟨_, hs0⟩⟩
        have haa : a = 4 * e + t.qx e := by
          rcases hslot 0 a ha with ⟨_, rfl⟩ | ⟨hf, _⟩ | ⟨hf, _⟩ | ⟨hf, _⟩ <;> first | rfl | omega
        subst haa
        obtain ⟨la, lb, hla, hlb, hab, _⟩ := ho.pair
        exact ⟨_, hs1, la, lb, hla, hlb, hab⟩
      · refine ⟨fun a ha => ?_, fun b _ => ⟨_, hs2⟩⟩
        have haa : a = 4 * e + (2 - t.qx e) := by
          rcases hslot 2 a ha with ⟨hf, _⟩ | ⟨hf, _⟩ | ⟨_, rfl⟩ | ⟨hf, _⟩ <;> first | rfl | omega
        subst haa
        obtain ⟨la, lb, hla, hlb, hab, _⟩ := ho.pair'
        exact ⟨_, hs3, la, lb, hla, hlb, hab⟩
  have hface : t.CapFace g s' i RotationSystem.single := by
    intro a c ha hc'
    have haa : a = 4 * e + t.qx e := by
      rcases hslot 0 a ha with ⟨_, rfl⟩ | ⟨hf, _⟩ | ⟨hf, _⟩ | ⟨hf, _⟩ <;> first | rfl | omega
    have hcc : c = 4 * e + (2 - t.qx e) := by
      rcases hslot 2 c hc' with ⟨hf, _⟩ | ⟨hf, _⟩ | ⟨_, rfl⟩ | ⟨hf, _⟩ <;> first | rfl | omega
    subst haa hcc
    rw [hPi]
    refine ⟨t.qx e, 2 - t.qx e, hloc _ (by omega), hloc _ (by omega), ?_⟩
    rcases hx with hx | hx <;> rw [hx]
    · exact Relation.ReflTransGen.single (by rfl)
    · exact Relation.ReflTransGen.single (by rfl)
  have hvert : ∀ ne, t.toSpqrTree.capNe i = some ne → ∀ p, t.toSpqrTree.neOrig ne = some p →
      ∀ k q, s'.outerE[i]?.bind (fun o => o[k]?) = some (some q) →
        QE.vert g.edges.toList q = some (if k < 2 then p.1 else p.2) := by
    intro ne hne p hp k q hk
    have hne' : (t.toSpqrTree.neRange i).1 = ne := by
      unfold SpqrTree.capNe at hne
      split at hne
      · exact Option.some.inj hne
      · cases hne
    have hpe := hrep.q_endpoints e hene i hidx ne hne' p hp
    have hslot' := hslot k q hk
    cases hf : t.edgeFlipped[e]!
    · have hp' := hpe.1 hf
      have hx0 : t.qx e = 0 := by simp [qx, hf]
      rw [hx0] at hslot'
      rcases hslot' with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ <;>
        (rw [hvq _ (by omega)]; simp [hp'])
    · have hp' := hpe.2 hf
      have hx2 : t.qx e = 2 := by simp [qx, hf]
      rw [hx2] at hslot'
      rcases hslot' with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ <;>
        (rw [hvq _ (by omega)]; simp [hp'])
  have h' : t.GluedUpTo g i s' := by
    refine ⟨⟨⟨⟨⟨⟨?_, ?_, ?_, ?_, ?_⟩, ?_, ?_⟩, ?_⟩, ?_, ?_⟩, ?_⟩, ?_⟩
    · rw [hrot, h.rot_size]
    · rw [hs', setOuter4_outer_size, h.outer_size]
    · intro j hj q he'
      exact h.outer_unprocessed j (by omega) q ((hoth j (by omega) q).1 he')
    · intro q hq
      rw [hrot]
      exact h.unset q (fun j hj hjs => hq j (by omega) hjs)
    · intro j hj
      by_cases hji : j = i
      · subst j
        exact ⟨RotationSystem.single, hw⟩
      · have hjsucc : t.Maximal (i + 1) j := t.maximal_succ i j (by have := hj.1; omega) hj
        obtain ⟨ρj, hρj, haj, hoj, hpj⟩ := h.piece j hjsucc
        refine ⟨ρj, hρj, ?_, ?_, ?_⟩
        · intro q r hq hr; exact haj q r hq (by rwa [hrot] at hr)
        · intro q hq; rw [hrot, hoth j hji]; exact hoj q hq
        · simpa only [hlookup, ite_eq_right hji] using hpj
    · intro j hj
      rw [hs', setOuter4_row_size]
      exact h.outer_row_size j hj
    · intro j k q hk
      by_cases hji : j = i
      · subst j
        have hk4 : k < 4 := by
          rcases hslot k q hk with ⟨rfl, _⟩ | ⟨rfl, _⟩ | ⟨rfl, _⟩ | ⟨rfl, _⟩ <;> omega
        refine ⟨hk4, by rw [ht]; decide, fun p hp hpt => ?_⟩
        rcases hpt with hpt | hpt
        · exact ((hnoparent p hp).1 hpt).elim
        · exact ((hnoparent p hp).2 hpt).elim
      · rw [hlookup, ite_eq_right hji] at hk
        exact h.outer_slots j k q hk
    · intro j p v q hp hpt hv he'
      by_cases hji : j = i
      · subst j
        exact ((hnoparent p hp).2 hpt).elim
      · exact h.outer_at_vertex j p v q hp hpt hv ((hoth j hji q).1 he')
    · intro j k q hk
      by_cases hji : j = i
      · subst j
        rcases hslot k q hk with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ <;> omega
      · rw [hlookup, ite_eq_right hji] at hk
        exact h.outer_dir j k q hk
    · intro j hj hjs p hp hpt hne
      by_cases hji : j = i
      · subst j
        exact ((hnoparent p hp).2 hpt).elim
      · obtain ⟨q, hq⟩ := h.outer_present j (by omega) hjs p hp hpt hne
        exact ⟨q, (hoth j hji q).2 hq⟩
    · intro j hjs hjt v hv
      by_cases hji : j = i
      · subst j; rw [ht] at hjt; cases hjt
      · obtain ⟨hvert', hslots, hpres⟩ := h.outer_vertex j hjs hjt v hv
        refine ⟨?_, ?_, ?_⟩
        · intro q hq; exact hvert' q ((hoth j hji q).1 hq)
        · intro k q hk
          rw [hlookup, ite_eq_right hji] at hk
          exact hslots k q hk
        · intro hj hne
          obtain ⟨q, hq⟩ := hpres (by omega) hne
          exact ⟨q, (hoth j hji q).2 hq⟩
    · intro j hjs ne hc' p hp
      by_cases hji : j = i
      · subst j
        refine ⟨hvert ne hc' p hp, fun _ _ k hk => ?_⟩
        interval_cases k
        · exact ⟨_, hs0⟩
        · exact ⟨_, hs1⟩
        · exact ⟨_, hs2⟩
        · exact ⟨_, hs3⟩
      · obtain ⟨hvert', hpres⟩ := h.outer_cap j hjs ne hc' p hp
        simp only [hlookup, ite_eq_right hji]
        refine ⟨hvert', ?_⟩
        intro hj hne
        exact hpres (by omega) hne
  refine t.gluedFaces_of_frame h h' ?_ ?_ ?_
  · intro ne _; exact ⟨RotationSystem.single, hw, hface⟩
  · intro j _ _ q _; rfl
  · intro j hji
    rw [hs', t.setOuter4_outer_ne i s j hji]

/-- The `Q` step of `planarEmbed`, preserving the same-witness cofacial certificate. -/
theorem embedItem_step_Q_faces (g : Graph) (hg : g.WF) (hwf : t.toSpqrTree.WF)
    (hsh : t.toSpqrTree.ChildShape) (hrep : t.toSpqrTree.Represents g)
    (hsep : t.toSpqrTree.PieceSep g) (i : Nat) (hi : i < t.size) (hty : t.types[i]! = .Q)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    t.GluedFaces g i ((t.embedItem i).run s).2 := by
  have ht : t.toSpqrTree.type i = .Q := by rw [t.type_eq_of_lt i hi]; exact hty
  obtain ⟨e, he, hidx⟩ := hwf.bij.edge_orig i hi ht
  rcases hsep.q_shape i hi ht with hc | ⟨c, hcO, hc⟩ | ⟨c, w, hcV, hcO, hcap, hwV, hc⟩
  · have hc' : t.children i = [] := by rw [← t.children_eq]; exact hc
    rw [t.embedItem_Q_nil i hty hc']
    exact t.q_leaf_glued g hg hwf hsh hrep hsep hi ht he hidx hc' s h
  · have hc' : t.children i = [c] := by rw [← t.children_eq]; exact hc
    obtain ⟨hcs, _, _⟩ := t.child_data (j := c) hwf hi (by rw [hc']; simp)
    have hcO' : t.types[c]! = .O := by rw [← t.type_eq_of_lt c hcs]; exact hcO
    rw [t.embedItem_Q_O i c [] hty hc' hcO']
    exact t.q_loop_glued g hg hwf hsh hrep hsep hi ht he hidx hc' hcO s h
  · have hc' : t.children i = [c, w] := by rw [← t.children_eq]; exact hc
    obtain ⟨hcs, _, _⟩ := t.child_data (j := c) hwf hi (by rw [hc']; simp)
    have hcO' : t.types[c]! ≠ .O := by rw [← t.type_eq_of_lt c hcs]; exact hcO
    rw [t.embedItem_Q_cons i c [w] hty hc' hcO']
    exact t.q_block_glued g hg hwf hsh hrep hsep hi ht he hidx hc' hcV hcO hcap hwV s h

end Spqr.PlanarSpqrTree
