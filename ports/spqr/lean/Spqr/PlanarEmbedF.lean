import Spqr.PlanarEmbedCloseList
import Spqr.PlanarEmbedEdges
import Spqr.Proofs.PieceUnion

namespace Spqr.PlanarSpqrTree

variable (t : PlanarSpqrTree)

theorem root_child_exposed_mem (g : Graph) (hwf : t.toSpqrTree.WF) (hi : 0 < t.size)
    (s : EmbedState) (h : t.GluedUpTo g 1 s) {j q : Nat}
    (hj : j ∈ t.children 0) (he : s.exposedAt j q) : (t.pieceBelow g j).Mem q := by
  obtain ⟨ρ, _, _, _, hp⟩ := h.piece j (t.child_maximal_one hwf hi hj)
  obtain ⟨k, hk⟩ := he
  have hkl := (h.outer_slots j k q hk).2.2 0 (t.child_data hwf hi hj).2.1
    (Or.inl hwf.preorder.root_type)
  interval_cases k
  · obtain ⟨b, _, la, lb, hla, _, _⟩ := (hp 0 (by omega)).1 q hk
    exact Piece.mem_of_loc hla
  · obtain ⟨a, ha⟩ := (hp 0 (by omega)).2 q hk
    obtain ⟨b, hb, la, lb, _, hlb, _⟩ := (hp 0 (by omega)).1 a ha
    have : b = q := Option.some.inj (Option.some.inj (hb.symm.trans hk))
    subst b
    exact Piece.mem_of_loc hlb

/-- Close the component boundary pairs and union their vertex-disjoint embeddings. -/
theorem embedItem_step_F (g : Graph) (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hrep : t.toSpqrTree.Represents g) (hsep : t.toSpqrTree.PieceSep g)
    (i : Nat) (hi : i < t.size) (hty : t.types[i]! = .F)
    (s : EmbedState) (h : t.GluedUpTo g (i + 1) s) :
    t.GluedUpTo g i ((t.embedItem i).run s).2 := by
  have hi0 : i = 0 := by
    by_contra hn
    exact hwf.preorder.only_root_F i (by omega) hi (by rw [t.type_eq_of_lt i hi]; exact hty)
  subst i
  rw [t.embedItem_F 0 hty s]
  have hnd := t.children_nodup hwf hi
  have hothers : ∀ j ∈ t.children 0, ∀ q, (t.pieceBelow g j).Mem q →
      ∀ k ∈ t.children 0, k ≠ j → ¬s.exposedAt k q := by
    intro j hj q hq k hk hkj he
    have hmem := t.root_child_exposed_mem g hwf hi s h hk he
    exact List.disjoint_left.1 (t.root_pieces_disjoint hwf g hrep.ne hsep hk hj hkj).1 hmem hq
  have hclosed : ∀ j ∈ t.children 0, ∃ ρ,
      IsPlanarEmbedding (t.pieceBelow g j).es g.nv ρ ∧
      (t.pieceBelow g j).Agrees (closeList (t.children 0) s).rotAdj ρ ∧
      ∀ q, (t.pieceBelow g j).Mem q → (closeList (t.children 0) s).rotAdj[q]? ≠ some none := by
    intro j hj
    obtain ⟨ρ, hρ, ha, ho, hp⟩ := h.piece j (t.child_maximal_one hwf hi hj)
    obtain ⟨hc, ht, _⟩ := closeOuter_spec (t.pieceBelow g j) ρ j s hρ
      (fun q hq => by rw [h.rot_size]; exact t.mem_pieceBelow_bound hwf g hq) ha ho
      (fun k q hk => (h.outer_slots j k q hk).2.2 0 (t.child_data hwf hi hj).2.1
        (Or.inl hwf.preorder.root_type)) (hp 0 (by omega))
    refine ⟨ρ, hρ, ?_, ?_⟩
    · intro q r hq hr
      rw [closeList_get _ hnd j hj s q (hothers j hj q hq)] at hr
      exact hc q r hq hr
    · intro q hq
      rw [closeList_get _ hnd j hj s q (hothers j hj q hq)]
      exact ht q hq
  have hes : t.edgesBelow 0 = (t.children 0).flatMap t.edgesBelow := by
    simpa [hty] using t.edgesBelow_eq hwf hsh 0 hi
  have hpdis : (t.children 0).Pairwise fun j k => List.Disjoint (t.edgesBelow j) (t.edgesBelow k) ∧
      ∀ v, HasEdge (t.pieceBelow g j).es v → HasEdge (t.pieceBelow g k).es v → False := by
    apply hnd.imp_of_mem
    intro j k hj hk hjk
    exact t.root_pieces_disjoint hwf g hrep.ne hsep hj hk hjk
  obtain ⟨ρ, hρ, haρ⟩ := Piece.union_embeddings (t.pieceBelow g 0) (t.children 0) t.edgesBelow
    (closeList (t.children 0) s).rotAdj
    (fun j hj => by obtain ⟨ρ, hρ, haρ, _⟩ := hclosed j hj; exact ⟨ρ, hρ, haρ⟩) hpdis
  have hρ' : IsPlanarEmbedding (t.pieceBelow g 0).es g.nv ρ := by
    simpa only [pieceBelow, ← hes] using hρ
  have haρ' : (t.pieceBelow g 0).Agrees (closeList (t.children 0) s).rotAdj ρ := by
    simpa only [pieceBelow, ← hes] using haρ
  have hroot : ∀ q, ¬(closeList (t.children 0) s).exposedAt 0 q := by
    intro q
    simpa only [EmbedState.exposedAt, closeList_outerE] using h.outer_unprocessed 0 (by omega) q
  refine ⟨⟨?_, ?_, ?_, ?_, ?_⟩, ?_, ?_⟩
  · rw [closeList_rotAdj_size, h.rot_size]
  · rw [closeList_outerE, h.outer_size]
  · intro j hj; omega
  · intro q hq
    rw [closeList_frame]
    · exact h.unset q (fun j hj hjs => hq j (by omega) hjs)
    · intro j hj he
      exact hq j (by omega) (t.child_data hwf hi hj).1
        (t.root_child_exposed_mem g hwf hi s h hj he)
  · intro j hj
    have hj0 := t.maximal_zero_eq hwf hj
    subst j
    refine ⟨ρ, hρ', haρ', ?_, ?_⟩
    · intro q hq
      have hq' : QE.edge q ∈ (t.children 0).flatMap t.edgesBelow := by
        rw [← hes]; exact hq
      obtain ⟨j, hj, hqj⟩ := List.mem_flatMap.1 hq'
      obtain ⟨_, _, _, ht⟩ := hclosed j hj
      exact ⟨fun hh => False.elim (ht q hqj hh), fun hh => False.elim (hroot q hh)⟩
    · intro k hk
      exact ⟨fun a ha => False.elim (hroot a ⟨_, ha⟩),
        fun b hb => False.elim (hroot b ⟨_, hb⟩)⟩
  · intro j hj; rw [closeList_outerE]; exact h.outer_row_size j hj
  · intro j k q hk
    rw [closeList_outerE] at hk
    exact h.outer_slots j k q hk

end Spqr.PlanarSpqrTree
