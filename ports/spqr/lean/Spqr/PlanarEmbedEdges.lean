import Spqr.PlanarEmbedBelow

namespace Spqr.PlanarSpqrTree

variable (t : PlanarSpqrTree)

theorem type_eq_of_lt (i : Nat) (hi : i < t.size) : t.toSpqrTree.type i = t.types[i]! := by
  unfold SpqrTree.type
  rw [Array.getElem?_eq_getElem hi, getElem!_pos t.types i hi]
  rfl

theorem orig_some_lt (hwf : t.toSpqrTree.WF) {j e : Nat} (h : t.origId[j]! = some e) :
    j < t.size := by
  by_contra hn
  have hj : ¬j < t.origId.size := by rwa [hwf.sizes.origId]
  simp [getElem!_neg, hj] at h

theorem edgeIndex_some_lt (hwf : t.toSpqrTree.WF) {e j : Nat} (h : t.edgeIndex[e]! = some j) :
    e < t.ne := by
  by_contra hn
  have he : ¬e < t.edgeIndex.size := by rwa [hwf.sizes.edgeIndex]
  simp [getElem!_neg, he] at h

theorem mem_edgesBelow_data (hwf : t.toSpqrTree.WF) {i e : Nat} (h : e ∈ t.edgesBelow i) :
    ∃ j, i ≤ j ∧ j < t.subtreeEnd[i]! ∧ j < t.size ∧
      t.types[j]! = .Q ∧ t.origId[j]! = some e ∧ t.edgeIndex[e]! = some j := by
  obtain ⟨j, hj, hje⟩ := List.mem_filterMap.1 h
  split at hje <;> rename_i hty
  · have hjs := t.orig_some_lt hwf hje
    have htype : t.types[j]! = .Q := by simpa using hty
    obtain ⟨e', he', hidx⟩ := hwf.bij.edge_orig j hjs (by rw [t.type_eq_of_lt j hjs]; exact htype)
    have heq : e' = e := Option.some.inj (he'.symm.trans hje)
    subst e'
    obtain ⟨k, hk, hjk⟩ := List.mem_range'.1 hj
    simp only [Nat.one_mul] at hjk
    exact ⟨j, by omega, by omega, hjs, htype, hje, hidx⟩
  · cases hje

theorem mem_edgesBelow_lt (hwf : t.toSpqrTree.WF) {i e : Nat} (h : e ∈ t.edgesBelow i) :
    e < t.ne := by
  obtain ⟨_, _, _, _, _, _, hidx⟩ := t.mem_edgesBelow_data hwf h
  exact t.edgeIndex_some_lt hwf hidx

theorem edgeIn_of_mem_edgesBelow (hwf : t.toSpqrTree.WF) {i e : Nat}
    (h : e ∈ t.edgesBelow i) : t.toSpqrTree.EdgeIn i e := by
  obtain ⟨j, hij, hje, _, _, _, hidx⟩ := t.mem_edgesBelow_data hwf h
  exact ⟨j, hidx, hij, by rwa [← getElem!_nat]⟩

theorem edgesBelow_nodup (hwf : t.toSpqrTree.WF) (i : Nat) : (t.edgesBelow i).Nodup := by
  apply List.Nodup.filterMap ?_ (List.nodup_range')
  intro j k e hj hk
  split at hj <;> rename_i htyj
  · split at hk <;> rename_i htyk
    · have hj' := t.orig_some_lt hwf hj
      have hk' := t.orig_some_lt hwf hk
      obtain ⟨ej, hej, hij⟩ := hwf.bij.edge_orig j hj'
        (by rw [t.type_eq_of_lt j hj']; simpa using htyj)
      obtain ⟨ek, hek, hik⟩ := hwf.bij.edge_orig k hk'
        (by rw [t.type_eq_of_lt k hk']; simpa using htyk)
      have hje : ej = e := Option.some.inj (hej.symm.trans hj)
      have hke : ek = e := Option.some.inj (hek.symm.trans hk)
      subst ej; subst ek
      exact Option.some.inj (hij.symm.trans hik)
    · cases hk
  · cases hj

theorem child_data (hwf : t.toSpqrTree.WF) {i j : Nat} (hi : i < t.size)
    (hj : j ∈ t.children i) : j < t.size ∧ t.toSpqrTree.parent j = some i ∧ i < j := by
  rw [← t.children_eq, hwf.preorder.ch_eq i hi, List.mem_filter] at hj
  have hjs : j < t.size := List.mem_range.1 hj.1
  have hjp : t.toSpqrTree.parent j = some i := by simpa using hj.2
  have hj0 : 0 < j := by
    by_contra hn
    have : j = 0 := by omega
    subst j
    rw [hwf.preorder.root_par] at hjp
    cases hjp
  obtain ⟨p, hp, hpi⟩ := hwf.preorder.par_lt j hj0 hjs
  have : p = i := Option.some.inj (hp.symm.trans hjp)
  subst p
  exact ⟨hjs, hjp, hpi⟩

theorem children_nodup (hwf : t.toSpqrTree.WF) {i : Nat} (hi : i < t.size) :
    (t.children i).Nodup := by
  rw [← t.children_eq, hwf.preorder.ch_eq i hi]
  exact List.nodup_range.filter _

theorem parent_some_iff (j p : Nat) :
    t.toSpqrTree.parent j = some p ↔ t.par[j]? = some (some p) := by
  unfold SpqrTree.parent
  cases t.par[j]? <;> simp

theorem child_maximal_one (hwf : t.toSpqrTree.WF) (hi : 0 < t.size) {j : Nat}
    (hj : j ∈ t.children 0) : t.Maximal 1 j := by
  obtain ⟨hjs, hjp, hj0⟩ := t.child_data hwf hi hj
  refine ⟨hj0, hjs, ?_⟩
  intro p hp
  have h0 := (t.parent_some_iff j 0).1 hjp
  have hp0 : p = 0 := Option.some.inj (Option.some.inj (hp.symm.trans h0))
  omega

theorem maximal_zero_eq (hwf : t.toSpqrTree.WF) {j : Nat} (h : t.Maximal 0 j) : j = 0 := by
  by_contra hn
  obtain ⟨p, hp, _⟩ := hwf.preorder.par_lt j (by omega) h.2.1
  have := h.2.2 p ((t.parent_some_iff j p).1 hp)
  omega

theorem mem_pieceBelow_bound (hwf : t.toSpqrTree.WF) (g : Graph) {i q : Nat}
    (h : (t.pieceBelow g i).Mem q) : q < 4 * t.ne := by
  have he := t.mem_edgesBelow_lt hwf h
  unfold QE.edge at he
  omega

theorem touches_of_hasEdge (hwf : t.toSpqrTree.WF) (g : Graph)
    (hne : t.ne = g.ne) {i v : Nat} (h : HasEdge (t.pieceBelow g i).es v) :
    t.toSpqrTree.Touches g i v := by
  obtain ⟨p, hp, hv⟩ := h
  obtain ⟨e, he, rfl⟩ := List.mem_map.1 hp
  refine ⟨e, ?_, t.edgeIn_of_mem_edgesBelow hwf he⟩
  exact ⟨by rw [← hne]; exact t.mem_edgesBelow_lt hwf he, hv⟩

theorem root_pieces_disjoint (hwf : t.toSpqrTree.WF) (g : Graph)
    (hne : t.ne = g.ne) (hsep : t.toSpqrTree.PieceSep g)
    {a b : Nat} (ha : a ∈ t.children 0) (hb : b ∈ t.children 0) (hab : a ≠ b) :
    List.Disjoint (t.edgesBelow a) (t.edgesBelow b) ∧
      ∀ v, HasEdge (t.pieceBelow g a).es v → HasEdge (t.pieceBelow g b).es v → False := by
  have hv : ∀ v, HasEdge (t.pieceBelow g a).es v → HasEdge (t.pieceBelow g b).es v → False := by
    intro v hva hvb
    apply hsep.root_disjoint a (by rwa [t.children_eq]) b (by rwa [t.children_eq]) hab v
    · exact t.touches_of_hasEdge hwf g hne hva
    · exact t.touches_of_hasEdge hwf g hne hvb
  refine ⟨?_, hv⟩
  rw [List.disjoint_left]
  intro e hea heb
  exact hv (g.edges[e]!).1
    ⟨g.edges[e]!, List.mem_map.2 ⟨e, hea, rfl⟩, Or.inl rfl⟩
    ⟨g.edges[e]!, List.mem_map.2 ⟨e, heb, rfl⟩, Or.inl rfl⟩

end Spqr.PlanarSpqrTree
