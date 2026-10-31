import Spqr.PlanarRelabelClosed
import Spqr.PlanarRelabelCornersSP
import Spqr.PlanarEmbedNodeR

/-!
# `NodeCorners` of the planar relabel output

`R` nodes: the node's `nodeEdges` rows and `neRotAdj` rows are `NodeRow`
(`planarRelabelTree_nodeRow`), so a corner at node-vertex `nvSt + j` is exactly an edge child in
`cornerVes` of the `j`-th node vertex, which `PlanarFinish.corners` pins down: none at the cap
endpoints, exactly one elsewhere, facing an edge that sorts later in `Items.ordered` (so the
executable `ta < tb`). `S`/`P` nodes: `PlanarRelabelCornersSP.lean`.
-/

namespace Spqr

open PlanarSpqrTree

theorem get?_ne_of_get!_none {a : Array (Option Nat)} {i x : Nat} (h : a[i]! = none) :
    a[i]? ≠ some (some x) := by
  by_cases hi : i < a.size
  · rw [getElem!_pos a i hi] at h; rw [Array.getElem?_eq_getElem hi, h]; simp
  · rw [Array.getElem?_eq_none_iff.2 (by omega)]; simp

theorem and_two_cases (o : Nat) : o &&& 2 = 0 ∨ o &&& 2 = 2 := by
  have h := Nat.and_two_pow o 1
  rw [Nat.pow_one] at h
  rw [h]
  cases o.testBit 1 <;> simp

theorem nvsArr_get! (s : PlanarRelabelState) {k : Nat} (hk : k < s.base.nodeEdges.size) :
    (nvsArr s)[k]! = s.base.nodeEdges[k]!.nvs := by
  unfold nvsArr
  rw [getElem!_pos (s.base.nodeEdges.map (·.nvs)) k (by simpa using hk),
    getElem!_pos s.base.nodeEdges k hk, Array.getElem_map]

theorem edgeChildren_get! (g : Graph) (items : Items) (pos : Nat → Nat) (children : List ItemId)
    {r : Nat} (hr : r < (children.filter (· ≥ 1 + g.nv)).length) :
    (Items.edgeChildren g items pos children)[r]! =
      (pos ((Items.vs items (children.filter (· ≥ 1 + g.nv))[r]).1.getD 0),
       pos ((Items.vs items (children.filter (· ≥ 1 + g.nv))[r]).2.getD 0)) := by
  unfold Items.edgeChildren
  rw [List.getElem!_eq_getElem?_getD, List.getElem?_map, List.getElem?_eq_getElem hr]
  rfl

/-- An edge child of an `R` item has two endpoints, both node vertices of the item. -/
theorem edgeChild_vs {g : Graph} {items : Items} (hwf : Items.WF g items) {it c : ItemId}
    (hlt : it < items.size) (hRt : Items.type items it = .R) (hc : c ∈ Items.ch items it)
    (hge : 1 + g.nv ≤ c) :
    ∃ a b, items[c]!.vs = (some a, some b) ∧ a ∈ Items.nvList g items it ∧
      b ∈ Items.nvList g items it ∧ c < items.size := by
  have hpar : Items.IsParent items it c := hc
  have hcs := hwf.tree.ch_lt it c hpar
  have hnV := child_not_V hwf hpar hge
  obtain ⟨a, b, hab⟩ := (hwf.shapes.r_shape it hlt hRt).2.2.2.2.2 c hpar hnV
  have hin : ∀ x, (Items.vs items c).1 = some x ∨ (Items.vs items c).2 = some x →
      x ∈ Items.nvList g items it := by
    intro x hx
    have hx' := hwf.endpoints.child_vs_in_parent it c hpar (by rw [hRt]; simp) hnV x hx
    have hxlt := hwf.endpoints.vs_lt c x hx
    unfold Items.nvList
    simp only [List.mem_append, List.mem_map, List.mem_filter, decide_eq_true_eq]
    rcases hx' with h1 | h1 | h1
    · exact Or.inl (Or.inl (by rw [h1]; simp))
    · exact Or.inr (by rw [h1]; simp)
    · refine Or.inl (Or.inr ⟨vertItem x, ⟨h1, ?_⟩, ?_⟩) <;> unfold vertItem <;> omega
  refine ⟨a, b, ?_, hin a (Or.inl (by rw [hab])), hin b (Or.inr (by rw [hab])), hcs⟩
  rw [← Items.vs_eq' hcs]; exact hab

theorem planarRelabelTree_cornersAt_R (g : Graph) (w : PlanarWalkState)
    (hwf : Items.WF g w.base.items) (hpf : PlanarFinish g w)
    (hall : (planarRelabelTree g w).nodePlanar.all id = true) {n : Nat}
    (hn : n < (planarRelabelTree g w).size)
    (hR : (planarRelabelTree g w).toSpqrTree.type n = .R) :
    (planarRelabelTree g w).CornersAt n := by
  obtain ⟨s, hinv, hT, hrows⟩ := planarRelabelTree_nodeRow g w hwf hpf hall
  rw [hT] at hn hR ⊢
  have hrow := hrows n hn hR
  unfold NodeRow at hrow
  obtain ⟨it, m, children, pos, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12⟩ := hrow
  set T : PlanarSpqrTree :=
    { SpqrTree.ofState g s.base with nodePlanar := s.aux.nodePlanar, neRotAdj := s.aux.neRotAdj }
    with hTdef
  have hcapne : T.toSpqrTree.capNe n = some (T.toSpqrTree.neRange n).1 := T.capNe_R hR
  have ene0 : (T.toSpqrTree.neRange n).1 = s.base.neBounds[n]! := getD_zero_eq_get! _ _
  have ene1 : (T.toSpqrTree.neRange n).2 = s.base.neBounds[n + 1]! := getD_zero_eq_get! _ _
  have env0 : (T.toSpqrTree.nvRange n).1 = s.base.nvBounds[n]! := getD_zero_eq_get! _ _
  have env1 : (T.toSpqrTree.nvRange n).2 = s.base.nvBounds[n + 1]! := getD_zero_eq_get! _ _
  have hedges : T.toSpqrTree.nodeEdges = s.base.nodeEdges := rfl
  have hrot : T.neRotAdj = s.aux.neRotAdj := rfl
  generalize hneSt : s.base.neBounds[n]! = neSt at *
  generalize hneEn : s.base.neBounds[n + 1]! = neEn at *
  generalize hnvSt : s.base.nvBounds[n]! = nvSt at *
  generalize hnvEn : s.base.nvBounds[n + 1]! = nvEn at *
  set items := w.base.items with hitems
  set Q := flipped g items[it]!.ch w.aux.itemFlips[it]! (capLinked g w.aux.qem m) with hQ
  set P := nodeSkel g items it with hP
  have hlt : it < items.size := h2
  have hRt : Items.type items it = .R := by rw [Items.type_eq' hlt]; exact h3
  have hcht : Items.ch items it = items[it]!.ch := Items.ch_eq' hlt
  have hperm_ch : items[it]!.ch.Perm children := by
    rw [h5, Items.ordered, ite_eq_right (by rw [hRt]; simp), hcht]
    exact (List.mergeSort_perm _ _).symm
  have hnd_vl : (Items.nvList g items it).Nodup := nvList_nodup hwf hlt
  have hit : 1 + g.nv + g.ne + (it - (1 + g.nv + g.ne)) = it := by omega
  have hty : items[1 + g.nv + g.ne + (it - (1 + g.nv + g.ne))]!.type ∈ [NodeType.S, .P, .R] := by
    rw [hit, h3]; simp
  have hcl := hpf.closed _ m h4 hty
  simp only [hit] at hcl
  change ∀ ve ∈ P.ves, ∀ z, z < 4 → ∃ o, Q[4 * ve + z]! = some o ∧ QE.edge o ∈ P.ves at hcl
  have hcorn := hpf.corners _ m h4 (by rw [hit]; exact h3)
  simp only [hit] at hcorn
  change ∀ j, j < (Items.nvList g items it).length →
      ((j = 0 ∨ j + 1 = (Items.nvList g items it).length) →
        cornerVes g items it Q (Items.nvList g items it)[j]! = []) ∧
      (0 < j → j + 1 < (Items.nvList g items it).length → ∃ ve o,
        cornerVes g items it Q (Items.nvList g items it)[j]! = [ve] ∧
        Q[4 * ve + 2]! = some o ∧ QE.edge o ≠ capVe g ∧
        loc₀ g items it (1 + g.nv + ve) < loc₀ g items it (1 + g.nv + QE.edge o)) at hcorn
  -- facts about the sorted edge children, before they are generalized away
  have hcs_mem : ∀ r (hr : r < (children.filter (· ≥ 1 + g.nv)).length),
      (children.filter (· ≥ 1 + g.nv))[r] ∈ items[it]!.ch ∧
        1 + g.nv ≤ (children.filter (· ≥ 1 + g.nv))[r] := by
    intro r hr
    have hm' := List.mem_filter.1 (List.getElem_mem hr)
    simp only [ge_iff_le, decide_eq_true_eq] at hm'
    exact ⟨hperm_ch.mem_iff.2 hm'.1, hm'.2⟩
  have hev_r : ∀ r (hr : r < (children.filter (· ≥ 1 + g.nv)).length),
      ((children.filter (· ≥ 1 + g.nv)).map (· - (1 + g.nv)))[r]? =
        some ((children.filter (· ≥ 1 + g.nv))[r] - (1 + g.nv)) := by
    intro r hr
    rw [List.getElem?_map, List.getElem?_eq_getElem hr]; rfl
  have hec_r : ∀ r (hr : r < (children.filter (· ≥ 1 + g.nv)).length),
      (Items.edgeChildren g items pos children)[r]! =
        (pos ((Items.vs items (children.filter (· ≥ 1 + g.nv))[r]).1.getD 0),
         pos ((Items.vs items (children.filter (· ≥ 1 + g.nv))[r]).2.getD 0)) :=
    fun r hr => edgeChildren_get! g items pos children hr
  have hnd_cs : (children.filter (· ≥ 1 + g.nv)).Nodup :=
    (hperm_ch.nodup_iff.1 (hcht ▸ hwf.tree.ch_nodup it)).filter _
  have hpw : children.Pairwise
      (fun a b => Items.loc g items nvSt pos a ≤ Items.loc g items nvSt pos b) := by
    rw [h5, Items.ordered, ite_eq_right (by rw [hRt]; simp)]
    refine List.Pairwise.imp (fun h => by simpa using h) (List.pairwise_mergeSort ?_ ?_ _)
    · intro a b c h1 h2; simp only [decide_eq_true_eq] at h1 h2 ⊢; omega
    · intro a b; simp only [Bool.or_eq_true, decide_eq_true_eq]; omega
  have hpw_cs := hpw.filter (· ≥ 1 + g.nv)
  generalize hcs : children.filter (· ≥ 1 + g.nv) = cs at *
  generalize hev : cs.map (· - (1 + g.nv)) = edgeVes at *
  set ves := capVe g :: edgeVes with hves
  have hperm : P.ves.Perm ves := by
    show (capVe g :: (items[it]!.ch.filter (· ≥ 1 + g.nv)).map (· - (1 + g.nv))).Perm
      (capVe g :: edgeVes)
    rw [← hev, ← hcs]
    exact List.Perm.cons _ ((hperm_ch.filter _).map _)
  have hlen_ev : edgeVes.length = cs.length := by rw [← hev]; simp
  have hlen_ves : ves.length = cs.length + 1 := by simp [hves, hlen_ev]
  have hne_size : neEn ≤ s.base.nodeEdges.size := by simpa [nvsArr] using h9
  have hrow_k : ∀ k, k < neEn - neSt → s.base.nodeEdges[neSt + k]!.nvs =
      ((nvSt, nvEn - 1) :: Items.edgeChildren g items pos children)[k]! := by
    intro k hk
    rw [← nvsArr_get! s (by omega)]; exact h10 k hk
  have hcap_nvs : s.base.nodeEdges[neSt]!.nvs = (nvSt, nvEn - 1) := by
    have := hrow_k 0 (by omega)
    rwa [Nat.add_zero, List.getElem!_cons_zero] at this
  have hrow_r : ∀ r (hr : r < cs.length), s.base.nodeEdges[neSt + 1 + r]!.nvs =
      (pos ((Items.vs items cs[r]).1.getD 0), pos ((Items.vs items cs[r]).2.getD 0)) := by
    intro r hr
    have := hrow_k (1 + r) (by omega)
    rw [show neSt + (1 + r) = neSt + 1 + r by omega, Nat.add_comm 1 r, List.getElem!_cons_succ,
      hec_r _ hr] at this
    exact this
  have hves_r : ∀ r (hr : r < cs.length), ves[1 + r]! = cs[r] - (1 + g.nv) := by
    intro r hr
    rw [hves, Nat.add_comm 1 r, List.getElem!_cons_succ, List.getElem!_eq_getElem?_getD, hev_r r hr]
    rfl
  have hev_mem : ∀ r (hr : r < cs.length), cs[r] - (1 + g.nv) ∈ edgeVes := fun r hr =>
    List.mem_iff_getElem?.2 ⟨r, hev_r r hr⟩
  have hQrow : ∀ r (hr : r < cs.length) o, Q[4 * (cs[r] - (1 + g.nv)) + 2]! = some o →
      QE.edge o ∈ ves ∧ s.aux.neRotAdj[4 * (neSt + 1 + r) + 2]! =
        some (4 * (neSt + ves.idxOf (QE.edge o)) + (o &&& 2) + 1) := by
    intro r hr o ho
    have hl : 4 * (1 + r) + 2 < 4 * ves.length := by rw [hlen_ves]; omega
    obtain ⟨-, hsome⟩ := h12 _ hl
    rw [show (4 * (1 + r) + 2) / 4 = 1 + r by omega, show (4 * (1 + r) + 2) % 4 = 2 by omega,
      hves_r r hr] at hsome
    have hmem : cs[r] - (1 + g.nv) ∈ P.ves := by
      rw [hperm.mem_iff, hves]; exact List.mem_cons_of_mem _ (hev_mem r hr)
    obtain ⟨o', ho', hoe⟩ := hcl _ hmem 2 (by omega)
    rw [ho] at ho'
    injection ho' with ho'
    subst ho'
    have hoe' : QE.edge o ∈ ves := hperm.mem_iff.mp hoe
    have := hsome o ho hoe'
    rw [show 4 * neSt + (4 * (1 + r) + 2) = 4 * (neSt + 1 + r) + 2 by omega,
      show 1 - (4 * (1 + r) + 2) % 2 = 1 by omega] at this
    exact ⟨hoe', this⟩
  have hQnone : ∀ r (hr : r < cs.length), Q[4 * (cs[r] - (1 + g.nv)) + 2]! = none →
      s.aux.neRotAdj[4 * (neSt + 1 + r) + 2]! = none := by
    intro r hr ho
    have hl : 4 * (1 + r) + 2 < 4 * ves.length := by rw [hlen_ves]; omega
    obtain ⟨hnone, -⟩ := h12 _ hl
    rw [show (4 * (1 + r) + 2) / 4 = 1 + r by omega, show (4 * (1 + r) + 2) % 4 = 2 by omega,
      hves_r r hr] at hnone
    rw [show 4 * (neSt + 1 + r) + 2 = 4 * neSt + (4 * (1 + r) + 2) by omega]
    exact hnone ho
  -- corners of `T` at node `n` are exactly the edge children with a side-0 partner at slot 2
  have hcorner : ∀ nv ta, T.CornerAt n nv ta ↔ ∃ r, ∃ hr : r < cs.length,
      ta = 4 * (neSt + 1 + r) + 2 ∧ pos ((Items.vs items cs[r]).2.getD 0) = nv ∧
      ∃ o, Q[4 * (cs[r] - (1 + g.nv)) + 2]! = some o ∧ o &&& 2 = 0 := by
    intro nv ta
    constructor
    · intro h
      obtain ⟨k, hk1, hk2, hta⟩ := T.cornerAt_decomp h
      rw [ene0, ene1] at hk2
      rw [ene0] at hta
      obtain ⟨-, -, -, ⟨d, hd, hdnv⟩, tb, htb, htb1⟩ := h
      refine ⟨k - 1, by omega, by omega, ?_, ?_⟩
      · rw [hedges, show ta / 4 = neSt + 1 + (k - 1) by omega] at hd
        rw [← getElem!_of_getElem? hd, hrow_r (k - 1) (by omega)] at hdnv
        exact hdnv
      · rw [hrot, hta, show 4 * (neSt + k) + 2 = 4 * (neSt + 1 + (k - 1)) + 2 by omega] at htb
        cases hQ' : Q[4 * (cs[k - 1] - (1 + g.nv)) + 2]! with
        | none => exact absurd htb (get?_ne_of_get!_none (hQnone _ (by omega) hQ'))
        | some o =>
          obtain ⟨-, hrw⟩ := hQrow _ (by omega) o hQ'
          rw [get?_of_get! hrw] at htb
          injection htb with htb
          injection htb with htb
          refine ⟨o, rfl, ?_⟩
          rcases and_two_cases o with h0 | h2
          · exact h0
          · exfalso; omega
    · rintro ⟨r, hr, hta, hnv, o, ho, h0⟩
      obtain ⟨hoe, hrw⟩ := hQrow r hr o ho
      have hlt' : neSt + 1 + r < s.base.nodeEdges.size := by omega
      refine ⟨by rw [ene0]; omega, by rw [ene1]; omega, by omega,
        ⟨s.base.nodeEdges[neSt + 1 + r]!, ?_, ?_⟩,
        4 * (neSt + ves.idxOf (QE.edge o)) + (o &&& 2) + 1, ?_, ?_⟩
      · rw [hedges, hta, show (4 * (neSt + 1 + r) + 2) / 4 = neSt + 1 + r by omega,
          getElem!_pos s.base.nodeEdges _ hlt']
        exact Array.getElem?_eq_getElem hlt'
      · rw [hrow_r r hr]; exact hnv
      · rw [hrot, hta]; exact get?_of_get! hrw
      · rw [h0]; omega
  have hcv : ∀ v ve, ve ∈ cornerVes g items it Q v ↔
      (∃ c ∈ items[it]!.ch, 1 + g.nv ≤ c ∧ ve = c - (1 + g.nv)) ∧
        (items[1 + g.nv + ve]!.vs).2 = some v ∧ ∃ o, Q[4 * ve + 2]! = some o ∧ o &&& 2 = 0 := by
    intro v ve
    unfold cornerVes
    rw [List.mem_filter, List.mem_map]
    simp only [List.mem_filter, ge_iff_le, decide_eq_true_eq, Bool.and_eq_true, beq_iff_eq]
    constructor
    · rintro ⟨⟨c, ⟨hc, hge⟩, rfl⟩, h2, h3⟩
      refine ⟨⟨c, hc, hge, rfl⟩, h2, ?_⟩
      generalize Q[4 * (c - (1 + g.nv)) + 2]! = q at h3 ⊢
      cases q with
      | none => simp at h3
      | some o => exact ⟨o, rfl, by simpa using h3⟩
    · rintro ⟨⟨c, hc, hge, rfl⟩, h2, o, ho, h0⟩
      refine ⟨⟨c, ⟨hc, hge⟩, rfl⟩, h2, ?_⟩
      rw [ho]; simpa using h0
  have hidx_iff : ∀ b ∈ Items.nvList g items it, ∀ j, j < (Items.nvList g items it).length →
      ((Items.nvList g items it).idxOf b = j ↔ b = (Items.nvList g items it)[j]!) := by
    intro b hb j hj
    constructor
    · intro h
      subst h
      rw [getElem!_pos (Items.nvList g items it) _ hj]
      exact (List.getElem_idxOf hj).symm
    · intro h
      exact Items.idxOf_eq_of_getElem? hnd_vl
        (by rw [h, getElem!_pos (Items.nvList g items it) _ hj]; exact List.getElem?_eq_getElem hj)
  have hloc_eq : ∀ c ∈ items[it]!.ch, 1 + g.nv ≤ c →
      Items.loc g items nvSt pos c = loc₀ g items it c := by
    intro c hc hge
    obtain ⟨a, b, hab, ha, hb, hclt⟩ := edgeChild_vs hwf hlt hRt (by rw [hcht]; exact hc) hge
    have hnc : ¬ c < 1 + g.nv := Nat.not_lt.2 hge
    unfold loc₀ Items.loc
    rw [ite_eq_right hnc, ite_eq_right hnc, Items.vs_eq' hclt, hab]
    simp only [Option.getD_some, Nat.sub_zero]
    rw [Items.PosOK.pos_sub hnd_vl h6 ha, Items.PosOK.pos_sub hnd_vl h6 hb]
  have hcapEnd : ∀ nv, T.toSpqrTree.CapEnd n nv ↔
      nv = nvSt ∨ nv = nvSt + (Items.nvList g items it).length - 1 := by
    intro nv
    constructor
    · rintro ⟨ne, d, hne, hd, hor⟩
      rw [hcapne, ene0] at hne
      injection hne with hne
      subst hne
      rw [hedges] at hd
      rw [← getElem!_of_getElem? hd, hcap_nvs] at hor
      simp only at hor
      omega
    · intro h
      refine ⟨neSt, s.base.nodeEdges[neSt]!, by rw [hcapne, ene0], ?_, ?_⟩
      · rw [hedges, getElem!_pos s.base.nodeEdges neSt (by omega)]
        exact Array.getElem?_eq_getElem (by omega)
      · rw [hcap_nvs]; simp only; omega
  refine ⟨fun nv ta hcap hcor => ?_, fun nv hnvof hncap => ?_⟩
  · rw [hcapEnd] at hcap
    rw [hcorner] at hcor
    obtain ⟨r, hr, -, hnv, o, ho, h0⟩ := hcor
    obtain ⟨hcch, hge⟩ := hcs_mem r hr
    obtain ⟨a, b, hab, ha, hb, hclt⟩ := edgeChild_vs hwf hlt hRt (by rw [hcht]; exact hcch) hge
    rw [Items.vs_eq' hclt, hab] at hnv
    simp only [Option.getD_some] at hnv
    have hpb := Items.PosOK.pos_sub hnd_vl h6 hb
    have hpb' := (h6 b hb).1
    have hjlt : (Items.nvList g items it).idxOf b < (Items.nvList g items it).length :=
      List.idxOf_lt_length_iff.2 hb
    have hj : (Items.nvList g items it).idxOf b = 0 ∨
        (Items.nvList g items it).idxOf b + 1 = (Items.nvList g items it).length := by omega
    have hempty := (hcorn _ hjlt).1 hj
    have hmem : cs[r] - (1 + g.nv) ∈
        cornerVes g items it Q (Items.nvList g items it)[(Items.nvList g items it).idxOf b]! := by
      rw [hcv]
      refine ⟨⟨cs[r], hcch, hge, rfl⟩, ?_, o, ho, h0⟩
      rw [show 1 + g.nv + (cs[r] - (1 + g.nv)) = cs[r] by omega, hab, getElem!_pos (Items.nvList g items it) _ hjlt,
        List.getElem_idxOf hjlt]
    rw [hempty] at hmem
    simp at hmem
  · unfold SpqrTree.NvOf at hnvof
    rw [env0, env1] at hnvof
    rw [hcapEnd] at hncap
    obtain ⟨hlo, hhi⟩ := hnvof
    have hj : nv - nvSt < (Items.nvList g items it).length := by omega
    have hj0 : 0 < nv - nvSt := by omega
    have hj1 : nv - nvSt + 1 < (Items.nvList g items it).length := by omega
    obtain ⟨ve, o, hcvs, ho, hne, hloc⟩ := (hcorn _ hj).2 hj0 hj1
    have hve : ve ∈ cornerVes g items it Q (Items.nvList g items it)[nv - nvSt]! := by
      rw [hcvs]; simp
    rw [hcv] at hve
    obtain ⟨⟨c, hcch, hge, rfl⟩, hvs2, o', ho', h0'⟩ := hve
    rw [ho] at ho'
    injection ho' with ho'
    subst ho'
    have hccs : c ∈ cs := by
      rw [← hcs]; exact List.mem_filter.2 ⟨hperm_ch.mem_iff.1 hcch, by simpa using hge⟩
    have hr : cs.idxOf c < cs.length := List.idxOf_lt_length_iff.2 hccs
    have hcr : cs[cs.idxOf c] = c := List.getElem_idxOf hr
    obtain ⟨a, b, hab, ha, hb, hclt⟩ := edgeChild_vs hwf hlt hRt (by rw [hcht]; exact hcch) hge
    rw [show 1 + g.nv + (c - (1 + g.nv)) = c by omega, hab] at hvs2
    simp only [Option.some.injEq] at hvs2
    have hidxb := (hidx_iff b hb _ hj).2 hvs2
    have hpb := Items.PosOK.pos_sub hnd_vl h6 hb
    have hpb' := (h6 b hb).1
    refine ⟨4 * (neSt + 1 + cs.idxOf c) + 2, ?_, ?_, ?_⟩
    · rw [hcorner]
      refine ⟨cs.idxOf c, hr, rfl, ?_, o, ?_, h0'⟩
      · rw [hcr, Items.vs_eq' hclt, hab]
        simp only [Option.getD_some]
        omega
      · rw [hcr]; exact ho
    · intro tb htb
      obtain ⟨hoe, hrw⟩ := hQrow _ hr o (by rw [hcr]; exact ho)
      rw [hrot, get?_of_get! hrw] at htb
      injection htb with htb
      injection htb with htb
      have hoe' : QE.edge o ∈ edgeVes := by
        rw [hves] at hoe
        rcases List.mem_cons.1 hoe with h | h
        · exact absurd h hne
        · exact h
      have hves_idx : ves.idxOf (QE.edge o) = edgeVes.idxOf (QE.edge o) + 1 := by
        rw [hves, List.idxOf_cons, ite_eq_right (by simpa using Ne.symm hne)]
      have hq : edgeVes.idxOf (QE.edge o) < cs.length := by
        rw [← hlen_ev]; exact List.idxOf_lt_length_iff.2 hoe'
      have hq' : edgeVes[edgeVes.idxOf (QE.edge o)]? = some (QE.edge o) := by
        rw [List.getElem?_eq_getElem (by rw [hlen_ev]; exact hq),
          List.getElem_idxOf (by rw [hlen_ev]; exact hq)]
      rw [hev_r _ hq] at hq'
      injection hq' with hq'
      have hgeq := (hcs_mem _ hq).2
      have hlt_idx : cs.idxOf c < edgeVes.idxOf (QE.edge o) := by
        by_contra hcon
        rcases Nat.lt_or_eq_of_le (Nat.le_of_not_lt hcon) with hlt' | heq
        · have hp := List.pairwise_iff_getElem.1 hpw_cs _ _ hq hr hlt'
          rw [hcr, hloc_eq _ (hcs_mem _ hq).1 hgeq, hloc_eq _ hcch hge] at hp
          have e1 : 1 + g.nv + QE.edge o = cs[edgeVes.idxOf (QE.edge o)] := by omega
          have e2 : 1 + g.nv + (c - (1 + g.nv)) = c := by omega
          rw [e1, e2] at hloc
          omega
        · have e : cs[edgeVes.idxOf (QE.edge o)]? = cs[cs.idxOf c]? := by rw [heq]
          rw [List.getElem?_eq_getElem hq, List.getElem?_eq_getElem hr, hcr] at e
          injection e with e
          rw [e] at hq'
          rw [← hq'] at hloc
          exact Nat.lt_irrefl _ hloc
      omega
    · intro ta' h'
      rw [hcorner] at h'
      obtain ⟨r', hr', hta', hnv', o₂, ho₂, h0₂⟩ := h'
      obtain ⟨hcch', hge'⟩ := hcs_mem r' hr'
      obtain ⟨a', b', hab', ha', hb', hclt'⟩ :=
        edgeChild_vs hwf hlt hRt (by rw [hcht]; exact hcch') hge'
      rw [Items.vs_eq' hclt', hab'] at hnv'
      simp only [Option.getD_some] at hnv'
      have hpb2 := Items.PosOK.pos_sub hnd_vl h6 hb'
      have hpb2' := (h6 b' hb').1
      have hidx : (Items.nvList g items it).idxOf b' = nv - nvSt := by omega
      have hmem' : cs[r'] - (1 + g.nv) ∈
          cornerVes g items it Q (Items.nvList g items it)[nv - nvSt]! := by
        rw [hcv]
        refine ⟨⟨cs[r'], hcch', hge', rfl⟩, ?_, o₂, ho₂, h0₂⟩
        rw [show 1 + g.nv + (cs[r'] - (1 + g.nv)) = cs[r'] by omega, hab',
          ← (hidx_iff b' hb' _ hj).1 hidx]
      rw [hcvs] at hmem'
      simp only [List.mem_singleton] at hmem'
      have hc' : cs[r'] = c := by omega
      have : cs.idxOf c = r' := by
        rw [← hc']; exact Items.idxOf_eq_of_getElem? hnd_cs (List.getElem?_eq_getElem hr')
      omega

/-- `NodeCorners` of the planar relabel output, under the all-planar flag. -/
theorem planarRelabelTree_nodeCorners (g : Graph) (w : PlanarWalkState)
    (hwf : Items.WF g w.base.items) (hpf : PlanarFinish g w)
    (hwfT : (planarRelabelTree g w).toSpqrTree.WF)
    (hlay : ∀ i, i < (planarRelabelTree g w).size → (planarRelabelTree g w).LayoutAt i)
    (hall : (planarRelabelTree g w).nodePlanar.all id = true) :
    (planarRelabelTree g w).NodeCorners := by
  refine nodeCorners_of_cornersAt _ fun i hi hty => ?_
  rcases hty with hS | hP | hR
  · exact (planarRelabelTree g w).cornersAt_S hwfT hi hS (hlay i hi)
  · exact (planarRelabelTree g w).cornersAt_P hwfT hi hP (hlay i hi)
  · exact planarRelabelTree_cornersAt_R g w hwf hpf hall hi hR

end Spqr
