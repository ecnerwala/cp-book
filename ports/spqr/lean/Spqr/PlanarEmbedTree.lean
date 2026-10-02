import Spqr.PlanarEmbedEdges
import Spqr.Proofs.PieceLoc

namespace Spqr.PlanarSpqrTree

variable (t : PlanarSpqrTree)

theorem mem_edgesBelow_of_edgeIn (hwf : t.toSpqrTree.WF) {i e : Nat}
    (h : t.toSpqrTree.EdgeIn i e) : e ∈ t.edgesBelow i := by
  obtain ⟨j, he, hji, hje⟩ := h
  obtain ⟨j', he', ht, ho⟩ := hwf.bij.edge_index e (t.edgeIndex_some_lt hwf he)
  have hjj : j' = j := Option.some.inj (he'.symm.trans he)
  subst j'
  have hjs := t.orig_some_lt hwf ho
  refine List.mem_filterMap.2 ⟨j, ?_, ?_⟩
  · rw [← getElem!_nat] at hje
    exact List.mem_range'.2 ⟨j - i, by omega, by omega⟩
  · rw [t.type_eq_of_lt j hjs] at ht
    simp [ht, ho]

theorem parent_lt (hwf : t.toSpqrTree.WF) {i p : Nat} (hi : i < t.size)
    (hp : t.toSpqrTree.parent i = some p) : p < i := by
  have hi0 : 0 < i := by
    by_contra hn
    have : i = 0 := by omega
    subst i
    rw [hwf.preorder.root_par] at hp
    cases hp
  obtain ⟨p', hp', hpi⟩ := hwf.preorder.par_lt i hi0 hi
  have : p' = p := Option.some.inj (hp'.symm.trans hp)
  omega

theorem child_maximal (hwf : t.toSpqrTree.WF) {i j : Nat} (hi : i < t.size)
    (hj : j ∈ t.children i) : t.Maximal (i + 1) j := by
  obtain ⟨hjs, hp, hij⟩ := t.child_data hwf hi hj
  refine ⟨by omega, hjs, ?_⟩
  intro p hp'
  have := (t.parent_some_iff j p).2 hp'
  have : p = i := Option.some.inj (this.symm.trans hp)
  omega

theorem v_children_nonempty (hwf : t.toSpqrTree.WF) {g : Graph}
    (hsep : t.toSpqrTree.PieceSep g) {i j : Nat} (hi : i < t.size)
    (ht : t.toSpqrTree.type i = .V) (hj : j ∈ t.children i) : t.edgesBelow j ≠ [] := by
  obtain ⟨e, he⟩ := hsep.v_nonempty i hi ht j (by rwa [t.children_eq])
  have hm := t.mem_edgesBelow_of_edgeIn hwf he
  intro hn
  rw [hn] at hm
  cases hm

theorem v_parent_not_v (hwf : t.toSpqrTree.WF) {i p : Nat} (hi : i < t.size)
    (ht : t.toSpqrTree.type i = .V) (hp : t.toSpqrTree.parent i = some p) :
    t.toSpqrTree.type p ≠ .V := by
  intro hpt
  obtain ⟨nv, _, hstart, hend, _⟩ := hwf.own.vert_par_nv i p hi ht hp
  have hpi := t.parent_lt hwf hi hp
  have hs := hwf.shape p (lt_trans hpi hi)
  simp only [SpqrTree.Shape, hpt] at hs
  omega

theorem maximal_not_inSubtree (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    {i j k : Nat} (hij : i ≤ j) (hjs : j < t.size) (hjk : j < k) (hk : t.Maximal i k) :
    ¬k < t.subtreeEnd[j]! := by
  generalize hd : k - j = d
  induction d using Nat.strong_induction_on generalizing j with
  | h d ih =>
    intro hke
    have hm : k ∈ List.range' (j + 1) (t.subtreeEnd[j]! - (j + 1)) :=
      List.mem_range'.2 ⟨k - (j + 1), by omega, by omega⟩
    rw [t.range'_children hwf hsh j hjs] at hm
    obtain ⟨c, hc, hck⟩ := List.mem_flatMap.1 hm
    obtain ⟨hcs, hcp, hjc⟩ := t.child_data hwf hjs hc
    obtain ⟨r, hr, hkr⟩ := List.mem_range'.1 hck
    simp only [Nat.one_mul] at hkr
    by_cases hkc : k = c
    · subst c
      have := hk.2.2 j ((t.parent_some_iff k j).1 hcp)
      omega
    · exact ih (k - c) (by omega) (by omega) hcs (by omega) rfl (by omega)

theorem maximal_pieces_disjoint (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    {i a b : Nat} (ha : t.Maximal i a) (hb : t.Maximal i b) (hab : a ≠ b) :
    List.Disjoint (t.edgesBelow a) (t.edgesBelow b) := by
  rw [List.disjoint_left]
  intro e hea heb
  obtain ⟨j, haj, hje, _, _, _, hej⟩ := t.mem_edgesBelow_data hwf hea
  obtain ⟨k, hbk, hke, _, _, _, hek⟩ := t.mem_edgesBelow_data hwf heb
  have : k = j := Option.some.inj (hek.symm.trans hej)
  subst k
  rcases lt_or_gt_of_ne hab with hab | hba
  · exact t.maximal_not_inSubtree hwf hsh ha.1 ha.2.1 hab hb (by omega)
  · exact t.maximal_not_inSubtree hwf hsh hb.1 hb.2.1 hba ha (by omega)

theorem pieceBelow_loc_vert (hwf : t.toSpqrTree.WF) {g : Graph} (hne : t.ne = g.ne)
    {i q l : Nat} (h : (t.pieceBelow g i).loc q = some l) :
    QE.vert (t.pieceBelow g i).es l = QE.vert g.edges.toList q := by
  rw [Piece.vert_loc h]
  have he : QE.edge q < g.edges.size := by
    change QE.edge q < g.ne
    rw [← hne]
    exact t.mem_edgesBelow_lt hwf (Piece.mem_of_loc h)
  simp only [QE.vert, Array.getElem?_toList, Array.getElem?_eq_getElem he, Option.map_some,
    pieceBelow, getElem!_pos g.edges (QE.edge q) he]

theorem v_pieces_meet (hwf : t.toSpqrTree.WF) {g : Graph} (hne : t.ne = g.ne)
    (hsep : t.toSpqrTree.PieceSep g) {i v a b : Nat} (hi : i < t.size)
    (ht : t.toSpqrTree.type i = .V) (hv : t.origId[i]! = some v)
    (ha : a ∈ t.children i) (hb : b ∈ t.children i) (hab : a ≠ b) :
    ∀ w, HasEdge (t.pieceBelow g a).es w → HasEdge (t.pieceBelow g b).es w → w = v := by
  intro w hwa hwb
  exact hsep.v_attach i v hi ht hv a (by rwa [t.children_eq]) b (by rwa [t.children_eq])
    hab w (t.touches_of_hasEdge hwf g hne hwa) (t.touches_of_hasEdge hwf g hne hwb)

end Spqr.PlanarSpqrTree
