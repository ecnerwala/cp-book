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

end Spqr.PlanarSpqrTree
