import Spqr.PlanarEmbedVLoop
import Spqr.PlanarEmbedLeaf

namespace Spqr.PlanarSpqrTree

variable (t : PlanarSpqrTree)

theorem vLoop_children (g : Graph) (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hne : t.ne = g.ne) (hsep : t.toSpqrTree.PieceSep g) {i v : Nat} (hi : i < t.size)
    (ht : t.toSpqrTree.type i = .V) (hv : t.origId[i]! = some v)
    (s : EmbedState) (h : t.GluedUpTo g (i + 1) s) :
    let r := vLoop (t.children i) s none none
    ((r.1 = (none, none) ∧ t.edgesBelow i = []) ∨
      ∃ a b ρ, r.1 = (some a, some b) ∧ (t.pieceBelow g i).OpenEmbedding r.2.rotAdj ρ a b v) ∧
    (∀ q, ¬(t.pieceBelow g i).Mem q → r.2.rotAdj[q]? = s.rotAdj[q]?) := by
  have hty : t.types[i]! = .V := by rwa [t.type_eq_of_lt i hi] at ht
  have hes : t.edgesBelow i = (t.children i).flatMap t.edgesBelow := by
    simpa [hty] using t.edgesBelow_eq hwf hsh i hi
  have hnd := t.children_nodup hwf hi
  have hpairs : (t.children i).Pairwise fun j k =>
      List.Disjoint (t.edgesBelow j) (t.edgesBelow k) ∧
        ∀ w, HasEdge (t.pieceBelow g j).es w → HasEdge (t.pieceBelow g k).es w → w = v := by
    apply hnd.imp_of_mem
    intro j k hj hk hjk
    exact ⟨t.maximal_pieces_disjoint hwf hsh (t.child_maximal hwf hi hj)
      (t.child_maximal hwf hi hk) hjk, t.v_pieces_meet hwf hne hsep hi ht hv hj hk hjk⟩
  cases hcs : t.children i with
  | nil =>
    rw [hcs] at hes
    exact ⟨Or.inl ⟨rfl, hes⟩, fun _ _ => rfl⟩
  | cons j js =>
    have hj : j ∈ t.children i := by rw [hcs]; simp
    obtain ⟨a, b, ρ, ha, hb, he⟩ := t.v_child_boundary g hwf hne hsep hi ht hv s h hj
    rw [hcs, List.pairwise_cons] at hpairs
    have hbound : ∀ q, ({t.pieceBelow g j with ves := t.edgesBelow j ++ js.flatMap t.edgesBelow} : Piece).Mem q →
        q < s.rotAdj.size := by
      intro q hq
      rw [h.rot_size]
      apply t.mem_pieceBelow_bound hwf g (i := i)
      simpa only [pieceBelow, Piece.Mem, hes, hcs, List.flatMap_cons] using hq
    have hchildren : ∀ k ∈ js, ∃ c d ρ, s.outerE[k]![0]! = some c ∧ s.outerE[k]![1]! = some d ∧
        ({t.pieceBelow g j with ves := t.edgesBelow k} : Piece).OpenEmbedding s.rotAdj ρ c d v := by
      intro k hk
      exact t.v_child_boundary g hwf hne hsep hi ht hv s h (by rw [hcs]; simp [hk])
    obtain ⟨hout, hf⟩ := vLoop_open_spec (t.pieceBelow g j) t.edgesBelow js s a b v ρ he hchildren
      hpairs.1 hpairs.2 hbound
    simp only [vLoop, Option.isNone_none, ↓reduceIte, ha, hb]
    obtain ⟨d, ρ', hr, he'⟩ := hout
    refine ⟨Or.inr ⟨a, d, ρ', hr, ?_⟩, ?_⟩
    · simpa only [pieceBelow, hes, hcs, List.flatMap_cons] using he'
    · simpa only [pieceBelow, hes, hcs, List.flatMap_cons] using hf

theorem vLoop_piece (g : Graph) (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hne : t.ne = g.ne) (hsep : t.toSpqrTree.PieceSep g) {i : Nat} (hi : i < t.size)
    (ht : t.toSpqrTree.type i = .V) (s : EmbedState) (h : t.GluedUpTo g (i + 1) s) :
    let r := vLoop (t.children i) s none none
    (∃ ρ, IsPlanarEmbedding (t.pieceBelow g i).es g.nv ρ ∧
      (t.pieceBelow g i).Agrees r.2.rotAdj ρ ∧
      (∀ q, (t.pieceBelow g i).Mem q → (r.2.rotAdj[q]? = some none ↔
        r.1.1 = some q ∨ r.1.2 = some q)) ∧
      (∀ a, r.1.1 = some a → ∃ b, r.1.2 = some b ∧ ∃ la lb,
        (t.pieceBelow g i).loc a = some la ∧ (t.pieceBelow g i).loc b = some lb ∧ ρ.get la = some lb) ∧
      (∀ b, r.1.2 = some b → ∃ a, r.1.1 = some a) ∧
      (∀ a, r.1.1 = some a → a % 2 = 0) ∧
      (∀ b, r.1.2 = some b → b % 2 = 1)) ∧
    (∀ q, ¬(t.pieceBelow g i).Mem q → r.2.rotAdj[q]? = s.rotAdj[q]?) := by
  obtain ⟨v, hv, _⟩ := hwf.bij.vert_orig i hi ht
  obtain ⟨hp, hf⟩ := t.vLoop_children g hwf hsh hne hsep hi ht hv s h
  refine ⟨?_, hf⟩
  rcases hp with ⟨hn, he⟩ | ⟨a, b, ρ, hr, hp⟩
  · have hm : ∀ q, ¬(t.pieceBelow g i).Mem q := by simp [pieceBelow, Piece.Mem, he]
    refine ⟨⟨#[]⟩, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
    · simpa [pieceBelow, Piece.es, he] using isPlanarEmbedding_nil g.nv
    · intro q r hq; exact (hm q hq).elim
    · intro q hq; exact (hm q hq).elim
    all_goals simp [hn]
  · refine ⟨ρ, hp.planar, hp.agrees, ?_, ?_, ?_, ?_, ?_⟩
    · intro q hq
      simpa only [hr, Option.some.injEq, eq_comm] using hp.unset q hq
    · intro a' ha'
      have : a = a' := Option.some.inj (by simpa only [hr] using ha')
      subst a'
      obtain ⟨la, lb, hla, hlb, hab, _⟩ := hp.boundary
      exact ⟨b, by rw [hr], la, lb, hla, hlb, hab⟩
    · intro _ _; exact ⟨a, by rw [hr]⟩
    · intro a' ha'
      have : a = a' := Option.some.inj (by simpa only [hr] using ha')
      exact this ▸ hp.left_dir
    · intro b' hb'
      have : b = b' := Option.some.inj (by simpa only [hr] using hb')
      exact this ▸ hp.right_dir

theorem embedItem_step_V (g : Graph) (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hrep : t.toSpqrTree.Represents g) (hsep : t.toSpqrTree.PieceSep g)
    (i : Nat) (hi : i < t.size) (hty : t.types[i]! = .V)
    (s : EmbedState) (h : t.GluedUpTo g (i + 1) s) :
    t.GluedUpTo g i ((t.embedItem i).run s).2 := by
  have ht : t.toSpqrTree.type i = .V := by rw [t.type_eq_of_lt i hi]; exact hty
  let r := vLoop (t.children i) s none none
  let s' := setOuterPair i r.1.1 r.1.2 r.2
  rw [t.embedItem_V i hty s]
  change t.GluedUpTo g i s'
  have hrot : s'.rotAdj = r.2.rotAdj := rfl
  have houter : r.2.outerE = s.outerE := vLoop_outerE ..
  obtain ⟨⟨ρ, hρ, haρ, hoρ, hpaira, hpairb, hdira, hdirb⟩, hframe⟩ :=
    t.vLoop_piece g hwf hsh hrep.ne hsep hi ht s h
  have hlookup : ∀ j k,
      s'.outerE[j]?.bind (fun o => o[k]?) =
        if j = i then if k = 0 then some r.1.1 else if k = 1 then some r.1.2
        else s.outerE[j]?.bind (fun o => o[k]?) else s.outerE[j]?.bind (fun o => o[k]?) := by
    intro j k
    rw [setOuterPair_lookup i _ _ r.2 (by rw [houter, h.outer_size]; exact hi)
      (by rw [houter]; exact h.outer_row_size i hi), houter]
  have hnewslot : ∀ k q, s'.outerE[i]?.bind (fun o => o[k]?) = some (some q) →
      (k = 0 ∧ r.1.1 = some q) ∨ (k = 1 ∧ r.1.2 = some q) := by
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
  have hexpose : ∀ q, s'.exposedAt i q ↔ r.1.1 = some q ∨ r.1.2 = some q := by
    intro q
    constructor
    · rintro ⟨k, hk⟩
      exact (hnewslot k q hk).imp And.right And.right
    · rintro (ha | hb)
      · exact ⟨0, by simp [hlookup, ha]⟩
      · exact ⟨1, by simp [hlookup, hb]⟩
  have hoth : ∀ j, j ≠ i → ∀ q, s'.exposedAt j q ↔ s.exposedAt j q := by
    intro j hji q
    simp only [EmbedState.exposedAt, hlookup, ite_eq_right hji]
  have himax : t.Maximal i i := ⟨le_rfl, hi, fun p hp =>
    t.parent_lt hwf hi ((t.parent_some_iff i p).2 hp)⟩
  refine ⟨⟨⟨⟨⟨?_, ?_, ?_, ?_, ?_⟩, ?_, ?_⟩, ?_⟩, ?_, ?_⟩, ?_⟩
  · rw [hrot, vLoop_rotAdj_size, h.rot_size]
  · rw [setOuterPair_outer_size, houter, h.outer_size]
  · intro j hj q he
    exact h.outer_unprocessed j (by omega) q ((hoth j (by omega) q).1 he)
  · intro q hq
    rw [hrot, hframe q (hq i le_rfl hi)]
    exact h.unset q (fun j hj hjs => hq j (by omega) hjs)
  · intro j hj
    by_cases hji : j = i
    · subst j
      refine ⟨ρ, hρ, haρ, ?_, ?_⟩
      · intro q hq
        rw [hexpose]
        exact hoρ q hq
      · intro k hk
        interval_cases k
        · constructor
          · intro a ha
            have ha' : r.1.1 = some a := by simpa [hlookup] using ha
            obtain ⟨b, hb, hp⟩ := hpaira a ha'
            exact ⟨b, by simpa [hlookup] using hb, hp⟩
          · intro b hb
            have hb' : r.1.2 = some b := by simpa [hlookup] using hb
            obtain ⟨a, ha⟩ := hpairb b hb'
            exact ⟨a, by simpa [hlookup] using ha⟩
        · constructor
          · intro a ha
            rcases hnewslot 2 a ha with ⟨hf, _⟩ | ⟨hf, _⟩ <;> omega
          · intro b hb
            rcases hnewslot 3 b hb with ⟨hf, _⟩ | ⟨hf, _⟩ <;> omega
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
  · intro j p v q hp hpt hv he
    by_cases hji : j = i
    · subst j
      exact (t.v_parent_not_v hwf hi ht hp hpt).elim
    · exact h.outer_at_vertex j p v q hp hpt hv ((hoth j hji q).1 he)
  · intro j k q hk
    by_cases hji : j = i
    · subst j
      rcases hnewslot k q hk with ⟨rfl, hq⟩ | ⟨rfl, hq⟩
      · exact hdira q hq
      · exact hdirb q hq
    · rw [hlookup, ite_eq_right hji] at hk
      exact h.outer_dir j k q hk
  · intro j hj hjs p hp hpt hne
    by_cases hji : j = i
    · subst j
      exact (t.v_parent_not_v hwf hi ht hp hpt).elim
    · obtain ⟨q, hq⟩ := h.outer_present j (by omega) hjs p hp hpt hne
      exact ⟨q, (hoth j hji q).2 hq⟩
  · intro j hjs hjt v hv
    by_cases hji : j = i
    · subst j
      obtain ⟨hc, _⟩ := t.vLoop_children g hwf hsh hrep.ne hsep hi ht hv s h
      refine ⟨?_, ?_, ?_⟩
      · intro q hq
        rw [hexpose] at hq
        rcases hc with ⟨hr, _⟩ | ⟨a, b, ρ, hr, he⟩
        · change r.1 = (none, none) at hr
          simp [hr] at hq
        · change r.1 = (some a, some b) at hr
          obtain ⟨la, lb, hla, hlb, hab, hva⟩ := he.boundary
          have hvl : QE.vert g.edges.toList a = some v :=
            (t.pieceBelow_loc_vert hwf hrep.ne hla).symm.trans hva
          have hvr : QE.vert g.edges.toList b = some v := by
            rw [← t.pieceBelow_loc_vert hwf hrep.ne hlb]
            exact (he.planar.same_vertex la (by
              rw [he.planar.size]
              simpa only [Piece.es, List.length_map] using Piece.loc_lt hla) lb hab).symm.trans hva
          rcases hq with hq | hq
          · have heq : a = q := Option.some.inj (by simpa only [hr] using hq)
            exact heq ▸ hvl
          · have heq : b = q := Option.some.inj (by simpa only [hr] using hq)
            exact heq ▸ hvr
      · intro k q hk
        rcases hnewslot k q hk with ⟨rfl, _⟩ | ⟨rfl, _⟩ <;> omega
      · intro _ hne
        rcases hc with ⟨_, he⟩ | ⟨a, b, ρ, hr, he⟩
        · exact (hne he).elim
        · exact ⟨a, (hexpose a).2 (Or.inl (by rw [hr]))⟩
    · obtain ⟨hvert, hslots, hpres⟩ := h.outer_vertex j hjs hjt v hv
      refine ⟨?_, ?_, ?_⟩
      · intro q hq; exact hvert q ((hoth j hji q).1 hq)
      · intro k q hk
        rw [hlookup, ite_eq_right hji] at hk
        exact hslots k q hk
      · intro hj hne
        obtain ⟨q, hq⟩ := hpres (by omega) hne
        exact ⟨q, (hoth j hji q).2 hq⟩

end Spqr.PlanarSpqrTree
