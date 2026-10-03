import Spqr.PlanarRelabelBridge
import Spqr.PlanarSoundR

/-!
# `NodeRotClosed` of the planar relabel output

The rows of a planar R node (`NodeRow`) are `some (4 * (neSt + idxOf …) + …)` wherever the `qem`
entry read is a quarter-edge of the node's skeleton, which `PlanarFinish.closed` guarantees; so
every entry of the node's `neRotAdj` segment points back into the segment.
-/

namespace Spqr

/-- `planarRelabel` from `init`: the tree is read off a final state satisfying `RowInv`. -/
theorem planarRelabelTree_rowInv (g : Graph) (w : PlanarWalkState)
    (hwf : Items.WF g w.base.items) (hpf : PlanarFinish g w) :
    ∃ s, RowInv g w s ∧ planarRelabelTree g w =
      { SpqrTree.ofState g s.base with nodePlanar := s.aux.nodePlanar, neRotAdj := s.aux.neRotAdj } := by
  have hroot : rootItem < w.base.items.size :=
    Nat.lt_of_lt_of_le (by omega : 0 < 1 + g.nv + g.ne) hwf.tree.size
  have hfr : Fresh g w rootItem (PlanarRelabelState.init g w) := fun _ _ _ _ _ _ _ _ => rfl
  obtain ⟨hinv, -⟩ := planarRelabel_rowInv g w hwf hpf w.base.items.size rootItem none none none _ hroot
    (RowInv.init g w) hfr
  exact ⟨_, hinv, rfl⟩

/-- `NodeRow` holds at every `R` node of the output tree (under the all-planar flag). -/
theorem planarRelabelTree_nodeRow (g : Graph) (w : PlanarWalkState)
    (hwf : Items.WF g w.base.items) (hpf : PlanarFinish g w)
    (hall : (planarRelabelTree g w).nodePlanar.all id = true) :
    ∃ s, RowInv g w s ∧ planarRelabelTree g w =
      { SpqrTree.ofState g s.base with nodePlanar := s.aux.nodePlanar, neRotAdj := s.aux.neRotAdj } ∧
      ∀ n, n < s.base.types.size → (SpqrTree.ofState g s.base).type n = .R → NodeRow g w s n := by
  obtain ⟨s, hinv, hT⟩ := planarRelabelTree_rowInv g w hwf hpf
  refine ⟨s, hinv, hT, fun n hn' hR => ?_⟩
  rw [hT] at hall
  have hR' : s.base.types[n]! = .R := by
    rw [getElem!_pos s.base.types n hn']
    have : s.base.types[n]?.getD .F = .R := hR
    rwa [Array.getElem?_eq_getElem hn'] at this
  have hn2 : n < s.aux.nodePlanar.size := by rw [hinv.np_size]; exact hn'
  have hpl' : s.aux.nodePlanar[n]! = true := by
    rw [getElem!_pos s.aux.nodePlanar n hn2]
    have : s.aux.nodePlanar[n]?.getD false = true := PlanarSpqrTree.isPlanar_of_all hall hn2
    rwa [Array.getElem?_eq_getElem hn2] at this
  exact hinv.rows n hn' hR' hpl'

theorem planarRelabelTree_nodeRotClosed (g : Graph) (w : PlanarWalkState)
    (hwf : Items.WF g w.base.items) (hpf : PlanarFinish g w)
    (hall : (planarRelabelTree g w).nodePlanar.all id = true) :
    (planarRelabelTree g w).NodeRotClosed := by
  obtain ⟨s, hinv, hT, hrows⟩ := planarRelabelTree_nodeRow g w hwf hpf hall
  rw [hT]
  intro n hn hR ta tb hlo hhi hget
  have hrow := hrows n hn hR
  unfold NodeRow at hrow
  obtain ⟨it, m, children, pos, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12⟩ := hrow
  have ene0 : ((SpqrTree.ofState g s.base).neRange n).1 = s.base.neBounds[n]! := getD_zero_eq_get! _ _
  have ene1 : ((SpqrTree.ofState g s.base).neRange n).2 = s.base.neBounds[n + 1]! := getD_zero_eq_get! _ _
  simp only at hlo hhi hget
  rw [ene0] at hlo
  rw [ene1] at hhi
  rw [ene0, ene1]
  generalize hneSt : s.base.neBounds[n]! = neSt at *
  generalize hneEn : s.base.neBounds[n + 1]! = neEn at *
  set items := w.base.items with hitems
  set edgeVes := (children.filter (· ≥ 1 + g.nv)).map (· - (1 + g.nv)) with hev
  set ves := capVe g :: edgeVes with hves
  set Q := flipped g items[it]!.ch w.aux.itemFlips[it]! (capLinked g w.aux.qem m) with hQ
  set P := nodeSkel g items it with hP
  have hget' : s.aux.neRotAdj[ta]! = some tb := by
    rw [getElem!_def, hget]
  have hlen : ves.length = neEn - neSt := by
    simp only [hves, List.length_cons]; omega
  obtain ⟨l, hl⟩ : ∃ l, ta = 4 * neSt + l := ⟨ta - 4 * neSt, by omega⟩
  subst hl
  have hl4 : l < 4 * ves.length := by rw [hlen]; omega
  obtain ⟨hnone, hsome⟩ := h12 l hl4
  -- children is a permutation of ch
  have hlt : it < items.size := h2
  have hRt : Items.type items it = .R := by rw [Items.type_eq' hlt]; exact h3
  have hcht : Items.ch items it = items[it]!.ch := Items.ch_eq' hlt
  have hperm_ch : items[it]!.ch.Perm children := by
    rw [h5, Items.ordered, ite_eq_right (by rw [hRt]; simp), hcht]
    exact (List.mergeSort_perm _ _).symm
  have hperm : P.ves.Perm ves := by
    show (capVe g :: (items[it]!.ch.filter (· ≥ 1 + g.nv)).map (· - (1 + g.nv))).Perm
      (capVe g :: edgeVes)
    exact List.Perm.cons _ ((hperm_ch.filter _).map _)
  have hit : 1 + g.nv + g.ne + (it - (1 + g.nv + g.ne)) = it := by omega
  have hty : items[1 + g.nv + g.ne + (it - (1 + g.nv + g.ne))]!.type ∈ [NodeType.S, .P, .R] := by
    rw [hit, h3]; simp
  have hcl := hpf.closed _ m h4 hty
  simp only [hit] at hcl
  change ∀ ve ∈ P.ves, ∀ z, z < 4 → ∃ o, Q[4 * ve + z]! = some o ∧ QE.edge o ∈ P.ves at hcl
  have hvl : l / 4 < ves.length := by omega
  have hmem : ves[l / 4]! ∈ P.ves := by
    rw [hperm.mem_iff, getElem!_pos ves (l / 4) hvl]
    exact List.getElem_mem hvl
  obtain ⟨o, ho, hoe⟩ := hcl _ hmem (l % 4) (Nat.mod_lt _ (by omega))
  have hoe' : QE.edge o ∈ ves := hperm.mem_iff.mp hoe
  have hrow := hsome o ho hoe'
  rw [hget'] at hrow
  have hidx : ves.idxOf (QE.edge o) < ves.length := List.idxOf_lt_length_of_mem hoe'
  have hand : o &&& 2 ≤ 2 := Nat.and_le_right
  have hmod : 1 - l % 2 ≤ 1 := Nat.sub_le _ _
  have htb : tb = 4 * (neSt + ves.idxOf (QE.edge o)) + (o &&& 2) + (1 - l % 2) := by
    injection hrow
  omega

end Spqr
