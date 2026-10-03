import Spqr.PlanarEmbedQ
import Spqr.PlanarNodeSpec
import Spqr.PlanarEmbedNodeExec

/-!
# The `S`/`P`/`R` step of `planarEmbed`

`node_capped_glued` assembles `GluedFaces g i` from a `Capped` certificate of the node's piece
whose four exposed ends fill the node's row, framed outside the piece; `nodeFold_capped` (the
semantic content of the node fold, admitted) supplies it; `embedItem_step_node_faces` combines them
through the executable unfolding `embedItem_node`.
-/

namespace Spqr.PlanarSpqrTree

variable (t : PlanarSpqrTree)

theorem node_capped_glued (g : Graph) (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hrep : t.toSpqrTree.Represents g) (hsep : t.toSpqrTree.PieceSep g) {i : Nat} (hi : i < t.size)
    (hty : t.toSpqrTree.type i = .S ∨ t.toSpqrTree.type i = .P ∨ t.toSpqrTree.type i = .R)
    {ne : Nat} (hcne : t.toSpqrTree.capNe i = some ne) {p : Nat × Nat}
    (hp : t.toSpqrTree.neOrig ne = some p)
    (s s' : EmbedState) (h : t.GluedFaces g (i + 1) s)
    (hsize : s'.rotAdj.size = s.rotAdj.size) (hosize : s'.outerE.size = s.outerE.size)
    (hframe : ∀ q, ¬(t.pieceBelow g i).Mem q → s'.rotAdj[q]? = s.rotAdj[q]?)
    (hout : ∀ j, j ≠ i → s'.outerE[j]? = s.outerE[j]?)
    {a0 a1 a2 a3 : Nat} (hrow : s'.outerE[i]? = some #[some a0, some a1, some a2, some a3])
    {ρ : RotationSystem} (hcap : (t.pieceBelow g i).Capped s'.rotAdj ρ a0 a1 a2 a3 p.1 p.2) :
    t.GluedFaces g i s' := by
  have htV : t.toSpqrTree.type i ≠ .V := by rcases hty with h' | h' | h' <;> rw [h'] <;> decide
  have htF : t.toSpqrTree.type i ≠ .F := by rcases hty with h' | h' | h' <;> rw [h'] <;> decide
  have hnoparent : ∀ q, t.toSpqrTree.parent i = some q →
      t.toSpqrTree.type q ≠ .F ∧ t.toSpqrTree.type q ≠ .V :=
    fun q hq => hsep.node_parent i q hi hty hq
  have hlookup : ∀ j k, s'.outerE[j]?.bind (fun o => o[k]?) =
      if j = i then (if k = 0 then some (some a0) else if k = 1 then some (some a1)
        else if k = 2 then some (some a2) else if k = 3 then some (some a3) else none)
      else s.outerE[j]?.bind (fun o => o[k]?) := by
    intro j k
    by_cases hji : j = i
    · subst hji
      rw [ite_eq_left rfl, hrow, Option.bind_some]
      rcases k with _ | _ | _ | _ | k <;> simp
    · rw [ite_eq_right hji, hout j hji]
  have hslot : ∀ k q, s'.outerE[i]?.bind (fun o => o[k]?) = some (some q) →
      (k = 0 ∧ q = a0) ∨ (k = 1 ∧ q = a1) ∨ (k = 2 ∧ q = a2) ∨ (k = 3 ∧ q = a3) := by
    intro k q hk
    rw [hlookup, ite_eq_left rfl] at hk
    split_ifs at hk with h0 h1 h2 h3
    · exact Or.inl ⟨h0, (Option.some.inj (Option.some.inj hk)).symm⟩
    · exact Or.inr (Or.inl ⟨h1, (Option.some.inj (Option.some.inj hk)).symm⟩)
    · exact Or.inr (Or.inr (Or.inl ⟨h2, (Option.some.inj (Option.some.inj hk)).symm⟩))
    · exact Or.inr (Or.inr (Or.inr ⟨h3, (Option.some.inj (Option.some.inj hk)).symm⟩))
  have hs0 : s'.outerE[i]?.bind (fun o => o[0]?) = some (some a0) := by rw [hlookup]; simp
  have hs1 : s'.outerE[i]?.bind (fun o => o[1]?) = some (some a1) := by rw [hlookup]; simp
  have hs2 : s'.outerE[i]?.bind (fun o => o[2]?) = some (some a2) := by rw [hlookup]; simp
  have hs3 : s'.outerE[i]?.bind (fun o => o[3]?) = some (some a3) := by rw [hlookup]; simp
  have hexpose : ∀ q, s'.exposedAt i q ↔ q = a0 ∨ q = a1 ∨ q = a2 ∨ q = a3 := by
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
  have hvl : ∀ {q l}, (t.pieceBelow g i).loc q = some l →
      QE.vert (t.pieceBelow g i).es l = QE.vert g.edges.toList q :=
    fun hl => t.pieceBelow_loc_vert hwf hrep.ne hl
  have hlt : ∀ {q l}, (t.pieceBelow g i).loc q = some l → l < ρ.size := by
    intro q l hl
    rw [hcap.planar.size]
    simpa only [Piece.es, List.length_map] using Piece.loc_lt hl
  have hvert4 : QE.vert g.edges.toList a0 = some p.1 ∧ QE.vert g.edges.toList a1 = some p.1 ∧
      QE.vert g.edges.toList a2 = some p.2 ∧ QE.vert g.edges.toList a3 = some p.2 := by
    obtain ⟨l0, l1, hl0, hl1, h01, hv0⟩ := hcap.pair0
    obtain ⟨l2, l3, hl2, hl3, h23, hv2⟩ := hcap.pair2
    refine ⟨(hvl hl0).symm.trans hv0, ?_, (hvl hl2).symm.trans hv2, ?_⟩
    · rw [← hvl hl1, ← hcap.planar.same_vertex l0 (hlt hl0) l1 h01]; exact hv0
    · rw [← hvl hl3, ← hcap.planar.same_vertex l2 (hlt hl2) l3 h23]; exact hv2
  have hw : t.PieceWitness g s' i ρ := by
    refine ⟨hcap.planar, hcap.agrees, ?_, ?_⟩
    · intro q hq; rw [hcap.unset q hq, hexpose]
    · intro k hk
      interval_cases k
      · refine ⟨fun a ha => ?_, fun b _ => ⟨_, hs0⟩⟩
        have haa : a = a0 := by
          rcases hslot 0 a ha with ⟨_, rfl⟩ | ⟨hf, _⟩ | ⟨hf, _⟩ | ⟨hf, _⟩ <;> first | rfl | omega
        subst haa
        obtain ⟨la, lb, hla, hlb, hab, _⟩ := hcap.pair0
        exact ⟨_, hs1, la, lb, hla, hlb, hab⟩
      · refine ⟨fun a ha => ?_, fun b _ => ⟨_, hs2⟩⟩
        have haa : a = a2 := by
          rcases hslot 2 a ha with ⟨hf, _⟩ | ⟨hf, _⟩ | ⟨_, rfl⟩ | ⟨hf, _⟩ <;> first | rfl | omega
        subst haa
        obtain ⟨la, lb, hla, hlb, hab, _⟩ := hcap.pair2
        exact ⟨_, hs3, la, lb, hla, hlb, hab⟩
  have hface : t.CapFace g s' i ρ := by
    intro a c ha hc'
    have haa : a = a0 := by
      rcases hslot 0 a ha with ⟨_, rfl⟩ | ⟨hf, _⟩ | ⟨hf, _⟩ | ⟨hf, _⟩ <;> first | rfl | omega
    have hcc : c = a2 := by
      rcases hslot 2 c hc' with ⟨hf, _⟩ | ⟨hf, _⟩ | ⟨_, rfl⟩ | ⟨hf, _⟩ <;> first | rfl | omega
    subst haa hcc
    exact hcap.face
  have hd0 := hcap.dir0
  have hd1 := hcap.dir1
  have hd2 := hcap.dir2
  have hd3 := hcap.dir3
  have hrow! : s'.outerE[i]! = #[some a0, some a1, some a2, some a3] := by
    rw [getElem!_def, hrow]
  have hdisj : ∀ j, t.Maximal i j → j ≠ i → ∀ q, (t.pieceBelow g j).Mem q →
      s'.rotAdj[q]? = s.rotAdj[q]? := by
    intro j hj hji q hq
    apply hframe q
    intro hqi
    exact List.disjoint_left.1 (t.maximal_pieces_disjoint hwf hsh hj himax hji) hq hqi
  have h' : t.GluedUpTo g i s' := by
    refine ⟨⟨⟨⟨⟨⟨?_, ?_, ?_, ?_, ?_⟩, ?_, ?_⟩, ?_⟩, ?_, ?_⟩, ?_⟩, ?_⟩
    · rw [hsize, h.rot_size]
    · rw [hosize, h.outer_size]
    · intro j hj q he'
      exact h.outer_unprocessed j (by omega) q ((hoth j (by omega) q).1 he')
    · intro q hq
      rw [hframe q (hq i le_rfl hi)]
      exact h.unset q (fun j hj hjs => hq j (by omega) hjs)
    · intro j hj
      by_cases hji : j = i
      · subst j
        exact ⟨ρ, hw⟩
      · have hjsucc : t.Maximal (i + 1) j := t.maximal_succ i j (by have := hj.1; omega) hj
        obtain ⟨ρj, hρj, haj, hoj, hpj⟩ := h.piece j hjsucc
        refine ⟨ρj, hρj, ?_, ?_, ?_⟩
        · intro q r hq hr; exact haj q r hq (by rwa [hdisj j hj hji q hq] at hr)
        · intro q hq; rw [hdisj j hj hji q hq, hoth j hji]; exact hoj q hq
        · simpa only [hlookup, ite_eq_right hji] using hpj
    · intro j hj
      by_cases hji : j = i
      · subst j; rw [hrow!]; rfl
      · rw [getElem!_def, hout j hji, ← getElem!_def]
        exact h.outer_row_size j hj
    · intro j k q hk
      by_cases hji : j = i
      · subst j
        have hk4 : k < 4 := by
          rcases hslot k q hk with ⟨rfl, _⟩ | ⟨rfl, _⟩ | ⟨rfl, _⟩ | ⟨rfl, _⟩ <;> omega
        refine ⟨hk4, htF, fun p hp hpt => ?_⟩
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
      · subst j; exact absurd hjt htV
      · obtain ⟨hvert', hslots, hpres⟩ := h.outer_vertex j hjs hjt v hv
        refine ⟨?_, ?_, ?_⟩
        · intro q hq; exact hvert' q ((hoth j hji q).1 hq)
        · intro k q hk
          rw [hlookup, ite_eq_right hji] at hk
          exact hslots k q hk
        · intro hj hne
          obtain ⟨q, hq⟩ := hpres (by omega) hne
          exact ⟨q, (hoth j hji q).2 hq⟩
    · intro j hjs ne' hc' p' hp'
      by_cases hji : j = i
      · subst j
        rw [hcne] at hc'
        cases hc'
        rw [hp] at hp'
        cases hp'
        refine ⟨fun k q hk => ?_, fun _ _ k hk => ?_⟩
        · rcases hslot k q hk with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ <;>
            simp [hvert4.1, hvert4.2.1, hvert4.2.2.1, hvert4.2.2.2]
        · interval_cases k
          · exact ⟨_, hs0⟩
          · exact ⟨_, hs1⟩
          · exact ⟨_, hs2⟩
          · exact ⟨_, hs3⟩
      · obtain ⟨hvert', hpres⟩ := h.outer_cap j hjs ne' hc' p' hp'
        simp only [hlookup, ite_eq_right hji]
        refine ⟨hvert', ?_⟩
        intro hj hne
        exact hpres (by omega) hne
  exact t.gluedFaces_of_frame h h' (fun _ _ => ⟨ρ, hw, hface⟩) hdisj hout

/-- Semantic content of the node fold: starting from `GluedFaces g (i + 1) s`, the fold of
`nodeStep` over the node's quarter-edges leaves `rotAdj` unchanged outside the node's piece, fills
the node's row with four exposed ends, and the piece is `Capped` there — planar witness agreeing
with `rotAdj`, the two exposed pairs facing at the cap's endpoints, slots 0/2 cofacial. Admitted:
(N1) the 2-sum of `hloc` with each child's `Capped` certificate through its twin (`Twins.noncap_children`,
`PieceSep.twin_orient`), (N2) iteration in `nodeRot` order with `Piece.loc` bookkeeping, (N3)
the inner `V` items at corners (1-sum) and cofaciality of the cap's sides; see PROOF.md §8.6. -/
theorem nodeFold_capped (g : Graph) (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hrep : t.toSpqrTree.Represents g) (hsep : t.toSpqrTree.PieceSep g) {i : Nat} (hi : i < t.size)
    (hty : t.toSpqrTree.type i = .S ∨ t.toSpqrTree.type i = .P ∨ t.toSpqrTree.type i = .R)
    (hloc : IsPlanarEmbedding (t.localSkeleton i) (t.toSpqrTree.nVerts i) (t.nodeRot i))
    {ne : Nat} (hcne : t.toSpqrTree.capNe i = some ne) {p : Nat × Nat}
    (hp : t.toSpqrTree.neOrig ne = some p)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    let s' := (List.range' (4 * t.neBounds[i]!) (4 * t.neBounds[i + 1]! - 4 * t.neBounds[i]!)).foldl
      (t.nodeStep i t.neBounds[i]!) s
    (∀ q, ¬(t.pieceBelow g i).Mem q → s'.rotAdj[q]? = s.rotAdj[q]?) ∧
    ∃ a0 a1 a2 a3 ρ, s'.outerE[i]? = some #[some a0, some a1, some a2, some a3] ∧
      (t.pieceBelow g i).Capped s'.rotAdj ρ a0 a1 a2 a3 p.1 p.2 := by
  sorry

/-- `S`/`P`/`R` step over `GluedFaces`: `embedItem_node` + `nodeFold_capped` + `node_capped_glued`. -/
theorem embedItem_step_node_faces (g : Graph) (hwf : t.toSpqrTree.WF)
    (hsh : t.toSpqrTree.ChildShape) (hrep : t.toSpqrTree.Represents g)
    (hsep : t.toSpqrTree.PieceSep g) (i : Nat) (hi : i < t.size)
    (hty : t.types[i]! = .S ∨ t.types[i]! = .P ∨ t.types[i]! = .R)
    (hloc : IsPlanarEmbedding (t.localSkeleton i) (t.toSpqrTree.nVerts i) (t.nodeRot i))
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    t.GluedFaces g i ((t.embedItem i).run s).2 := by
  have ht : t.toSpqrTree.type i = .S ∨ t.toSpqrTree.type i = .P ∨ t.toSpqrTree.type i = .R := by
    rw [t.type_eq_of_lt i hi]; exact hty
  have hcne : t.toSpqrTree.capNe i = some (t.toSpqrTree.neRange i).1 := by
    rcases ht with h' | h' | h' <;> simp [SpqrTree.capNe, SpqrTree.hasCap, h'] <;> decide
  obtain ⟨p, hp⟩ := hsep.cap_orig i _ hi hcne
  rw [t.embedItem_node i hty s]
  obtain ⟨hframe, a0, a1, a2, a3, ρ, hrow, hcap⟩ :=
    t.nodeFold_capped g hwf hsh hrep hsep hi ht hloc hcne hp s h
  exact t.node_capped_glued g hwf hsh hrep hsep hi ht hcne hp s _ h
    (t.nodeFold_rotAdj_size _ _ _ _) (t.nodeFold_outerE_size _ _ _ _) hframe
    (fun j hji => t.nodeFold_outerE_ne _ _ _ _ j hji) hrow hcap

end Spqr.PlanarSpqrTree
