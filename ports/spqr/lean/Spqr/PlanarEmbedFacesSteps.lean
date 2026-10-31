import Spqr.PlanarEmbedFaces
import Spqr.PlanarEmbedFold

/-!
# `GluedFaces` through the proved steps

Leaves, `F` and `V` carry no cap, and each frames every other maximal piece, so
`gluedFaces_of_frame` lifts their `GluedUpTo` steps.
-/

namespace Spqr

namespace PlanarSpqrTree

variable (t : PlanarSpqrTree)

theorem capNe_eq_none_of_not_node {i : Nat}
    (h : t.toSpqrTree.type i = .F ∨ t.toSpqrTree.type i = .V) : t.toSpqrTree.capNe i = none := by
  unfold SpqrTree.capNe SpqrTree.hasCap
  rcases h with h | h <;> rw [h] <;> rfl

theorem embedItem_step_leaf_faces (g : Graph) (hwf : t.toSpqrTree.WF)
    (hsh : t.toSpqrTree.ChildShape) (i : Nat) (hi : i < t.size)
    (hty : t.types[i]! = .O ∨ t.types[i]! = .I)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    t.GluedFaces g i ((t.embedItem i).run s).2 := by
  have h' := t.embedItem_step_leaf g hwf hsh i hi hty s h.toGluedUpTo
  rw [t.embedItem_leaf i hty s] at h' ⊢
  refine t.gluedFaces_of_frame h h' ?_ (fun _ _ _ _ _ => rfl) (fun _ _ => rfl)
  intro ne _
  obtain ⟨ρ, hw⟩ := h'.piece i (t.maximal_self hwf hi)
  exact ⟨ρ, hw, t.capFace_of_unexposed ρ (h.outer_unprocessed i (by omega))⟩

theorem embedItem_step_F_faces (g : Graph) (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hrep : t.toSpqrTree.Represents g) (hsep : t.toSpqrTree.PieceSep g)
    (i : Nat) (hi : i < t.size) (hty : t.types[i]! = .F)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    t.GluedFaces g i ((t.embedItem i).run s).2 := by
  have h' := t.embedItem_step_F g hwf hsh hrep hsep i hi hty s h.toGluedUpTo
  have hi0 : i = 0 := by
    by_contra hn
    exact hwf.preorder.only_root_F i (by omega) hi (by rw [t.type_eq_of_lt i hi]; exact hty)
  subst hi0
  refine t.gluedFaces_of_frame h h' ?_ ?_ ?_
  · intro ne hne
    rw [t.capNe_eq_none_of_not_node (Or.inl (by rw [t.type_eq_of_lt 0 hi]; exact hty))] at hne
    cases hne
  · intro j hm hj0
    exact absurd (t.maximal_zero_eq hwf hm) hj0
  · intro j _
    rw [t.embedItem_F 0 hty s, closeList_outerE]

theorem embedItem_step_V_faces (g : Graph) (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hrep : t.toSpqrTree.Represents g) (hsep : t.toSpqrTree.PieceSep g)
    (i : Nat) (hi : i < t.size) (hty : t.types[i]! = .V)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    t.GluedFaces g i ((t.embedItem i).run s).2 := by
  have h' := t.embedItem_step_V g hwf hsh hrep hsep i hi hty s h.toGluedUpTo
  have ht : t.toSpqrTree.type i = .V := by rw [t.type_eq_of_lt i hi]; exact hty
  obtain ⟨_, hframe⟩ := t.vLoop_piece g hwf hsh hrep.ne hsep hi ht s h.toGluedUpTo
  refine t.gluedFaces_of_frame h h' ?_ ?_ ?_
  · intro ne hne
    rw [t.capNe_eq_none_of_not_node (Or.inr ht)] at hne
    cases hne
  · intro j hm hji q hq
    rw [t.embedItem_V i hty s, setOuterPair_rotAdj]
    apply hframe q
    intro hqi
    exact t.maximal_pieces_disjoint hwf hsh hm (t.maximal_self hwf hi) hji hq hqi
  · intro j hji
    rw [t.embedItem_V i hty s, setOuterPair_outer_ne _ _ _ _ j hji, vLoop_outerE]

end PlanarSpqrTree

end Spqr
