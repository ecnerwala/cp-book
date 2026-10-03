import Spqr.PlanarRelabelNodeR
import Spqr.Proofs.PlanarReindex
import Spqr.Proofs.PlanarMap
import Spqr.PlanarWalkProj
import Spqr.WalkItemsWF
import Spqr.DfsFast

/-!
# `nodePlanar_sound_R` from `PlanarFinish` and `RelabelNodeR`

`PlanarFinish.embedded` is a planar embedding of the node's skeleton read in `ch` order;
`RelabelNodeR` lays the edge children out in `Items.ordered` order and maps the node's vertices
to positions `pos - nvSt`. The two are related by reordering `ves` (`readRot_perm`, a quarter-edge
reindexing) and relabelling the vertices (`IsPlanarEmbedding.map`).
-/

namespace Spqr

/-- `omega` after unfolding `ItemId` (arithmetic on `ItemId`-typed children is opaque to it). -/
macro "iomega" : tactic => `(tactic| ((try dsimp only [ItemId] at *); omega))

theorem readRot_size (P : Piece) (Q : Qem) : (readRot P Q).size = 4 * P.ves.length := by
  simp [readRot, RotationSystem.size]

theorem readRot_get (P : Piece) (Q : Qem) (l : Nat) (hl : l < 4 * P.ves.length) :
    (readRot P Q).get l =
      (Q[4 * P.ves[l / 4]! + l % 4]!).map fun o =>
        4 * P.ves.idxOf (QE.edge o) + (o &&& 2) + (1 - l % 2) := by
  simp [readRot, RotationSystem.get, hl]

theorem readRot_get_of_ge (P : Piece) (Q : Qem) (l : Nat) (hl : 4 * P.ves.length ≤ l) :
    (readRot P Q).get l = none := by
  simp [readRot, RotationSystem.get, hl]

/-- `readRot` reconstructs a planar embedding independently of the order of `ves`. -/
theorem readRot_perm {P : Piece} {E : List Nat} {Q : Qem} {n : Nat}
    (h : IsPlanarEmbedding P.es n (readRot P Q)) (hperm : P.ves.Perm E) (hnd : P.ves.Nodup)
    (hclosed : ∀ ve ∈ P.ves, ∀ z, z < 4 → ∃ o, Q[4 * ve + z]! = some o ∧ QE.edge o ∈ P.ves) :
    IsPlanarEmbedding ({P with ves := E} : Piece).es n (readRot {P with ves := E} Q) := by
  set P' : Piece := {P with ves := E} with hP'
  have hndE : E.Nodup := hperm.nodup_iff.1 hnd
  have hlen : P.ves.length = E.length := hperm.length_eq
  have hves' : P'.ves = E := rfl
  have hes : P.es.length = P.ves.length := by simp [Piece.es]
  have hes' : P'.es.length = E.length := by simp [Piece.es, P']
  have hρs : (readRot P Q).size = 4 * P.ves.length := readRot_size P Q
  have hσs : (readRot P' Q).size = 4 * E.length := readRot_size P' Q
  let φ : Nat → Nat := fun q => 4 * P.ves.idxOf E[q / 4]! + q % 4
  have hEk : ∀ q, q < 4 * E.length → E[q / 4]! ∈ E := by
    intro q hq
    rw [List.getElem!_eq_getElem?_getD, List.getElem?_eq_getElem (by omega)]
    exact List.getElem_mem _
  have hidx : ∀ e ∈ E, P.ves.idxOf e < P.ves.length := fun e he =>
    List.idxOf_lt_length_iff.2 (hperm.mem_iff.2 he)
  have hidxE : ∀ e ∈ E, E.idxOf e < E.length := fun e he => List.idxOf_lt_length_iff.2 he
  have hgetE : ∀ q, q < 4 * E.length → E[E.idxOf E[q / 4]!]! = E[q / 4]! := by
    intro q hq
    have := hidxE _ (hEk q hq)
    rw [List.getElem!_eq_getElem?_getD (i := E.idxOf _), List.getElem?_eq_getElem this]
    simp
  have hgetP : ∀ e ∈ E, P.ves[P.ves.idxOf e]! = e := by
    intro e he
    have := hidx e he
    rw [List.getElem!_eq_getElem?_getD, List.getElem?_eq_getElem this]
    simp [List.getElem_idxOf this]
  have hclosedE : ∀ q, q < 4 * E.length → ∃ o, Q[4 * E[q / 4]! + q % 4]! = some o ∧ QE.edge o ∈ E := by
    intro q hq
    obtain ⟨o, ho, hm⟩ := hclosed _ (hperm.mem_iff.2 (hEk q hq)) (q % 4) (Nat.mod_lt _ (by decide))
    exact ⟨o, ho, hperm.mem_iff.1 hm⟩
  have hy3 : ∀ o q : Nat, (o &&& 2) + (1 - q % 2) < 4 := by
    intro o q
    have := @Nat.and_le_right o 2
    omega
  refine h.reindex φ (by rw [hσs, hes']) (by rw [hes, hes', hlen]) ?_ ?_ ?_ ?_ ?_ ?_ ?_ ?_
  · intro p
    simp only [Piece.es, hves']
    exact (hperm.map P.ends).mem_iff
  · intro q hq
    rw [hσs] at hq
    rw [hρs]
    have := hidx _ (hEk q hq)
    show 4 * P.ves.idxOf E[q / 4]! + q % 4 < 4 * P.ves.length
    omega
  · intro q r hq hr he
    rw [hσs] at hq hr
    change 4 * P.ves.idxOf E[q / 4]! + q % 4 = 4 * P.ves.idxOf E[r / 4]! + r % 4 at he
    have h1 : P.ves.idxOf E[q / 4]! = P.ves.idxOf E[r / 4]! := by omega
    have h2 : E[q / 4]! = E[r / 4]! := by
      rw [← hgetP _ (hEk q hq), h1]
      exact hgetP _ (hEk r hr)
    have h3 : q / 4 = r / 4 := by
      rw [List.getElem!_eq_getElem?_getD, List.getElem!_eq_getElem?_getD,
        List.getElem?_eq_getElem (by omega), List.getElem?_eq_getElem (by omega)] at h2
      exact (List.Nodup.getElem_inj_iff hndE).1 h2
    omega
  · intro q hq
    rw [hσs] at hq
    have hφ : φ q < 4 * P.ves.length := by
      have := hidx _ (hEk q hq)
      show 4 * P.ves.idxOf E[q / 4]! + q % 4 < 4 * P.ves.length
      omega
    rw [readRot_get P Q _ hφ, readRot_get P' Q q (by rw [hves', ← hlen]; rw [hlen]; exact hq)]
    have hd : φ q / 4 = P.ves.idxOf E[q / 4]! := by
      show (4 * P.ves.idxOf E[q / 4]! + q % 4) / 4 = _; omega
    have hm : φ q % 4 = q % 4 := by
      show (4 * P.ves.idxOf E[q / 4]! + q % 4) % 4 = _; omega
    have hm2 : φ q % 2 = q % 2 := by
      show (4 * P.ves.idxOf E[q / 4]! + q % 4) % 2 = _; omega
    rw [hd, hm, hm2, hgetP _ (hEk q hq), hves']
    obtain ⟨o, ho, hoE⟩ := hclosedE q hq
    rw [ho]
    simp only [Option.map_some, Option.some.injEq]
    show _ = 4 * P.ves.idxOf E[(4 * E.idxOf (QE.edge o) + (o &&& 2) + (1 - q % 2)) / 4]! +
      (4 * E.idxOf (QE.edge o) + (o &&& 2) + (1 - q % 2)) % 4
    have h4 : (4 * E.idxOf (QE.edge o) + (o &&& 2) + (1 - q % 2)) / 4 = E.idxOf (QE.edge o) := by
      have := hy3 o q; omega
    have h5 : (4 * E.idxOf (QE.edge o) + (o &&& 2) + (1 - q % 2)) % 4 = (o &&& 2) + (1 - q % 2) := by
      have := hy3 o q; omega
    rw [h4, h5]
    have h6 : E[E.idxOf (QE.edge o)]! = QE.edge o := by
      have := hidxE _ hoE
      rw [List.getElem!_eq_getElem?_getD, List.getElem?_eq_getElem this]
      simp [List.getElem_idxOf this]
    rw [h6]; omega
  · intro q hq r hr
    rw [hσs] at hq
    rw [readRot_get P' Q q (by rw [hves']; exact hq)] at hr
    obtain ⟨o, ho, rfl⟩ := Option.map_eq_some_iff.1 hr
    rw [hves'] at ho
    obtain ⟨o', ho', hoE⟩ := hclosedE q hq
    rw [ho] at ho'
    cases ho'
    have := hidxE _ hoE
    have := hy3 o q
    rw [hσs]
    change 4 * E.idxOf (QE.edge o) + (o &&& 2) + (1 - q % 2) < 4 * E.length
    omega
  · intro q hq
    show (4 * P.ves.idxOf E[q / 4]! + q % 4) % 2 = q % 2
    omega
  · intro q hq
    rw [hσs] at hq
    unfold QE.vert QE.edge QE.side
    have hd : φ q / 4 = P.ves.idxOf E[q / 4]! := by
      show (4 * P.ves.idxOf E[q / 4]! + q % 4) / 4 = _; omega
    have hs : φ q / 2 % 2 = q / 2 % 2 := by
      show (4 * P.ves.idxOf E[q / 4]! + q % 4) / 2 % 2 = _; omega
    rw [hd, hs]
    simp only [Piece.es, hves', List.getElem?_map]
    rw [List.getElem?_eq_getElem (hidx _ (hEk q hq)), List.getElem?_eq_getElem (by omega)]
    have := hgetP _ (hEk q hq)
    rw [List.getElem!_eq_getElem?_getD, List.getElem?_eq_getElem (hidx _ (hEk q hq))] at this
    simp only [Option.getD_some] at this
    rw [this]
    rw [List.getElem!_eq_getElem?_getD, List.getElem?_eq_getElem (by omega)]
    rfl
  · intro q hq c hc
    show 4 * P.ves.idxOf E[(q ^^^ c) / 4]! + (q ^^^ c) % 4 = (4 * P.ves.idxOf E[q / 4]! + q % 4) ^^^ c
    rw [xor_div4 q hc, xor_mod4 q hc, mul4_add_xor _ _ (Nat.mod_lt _ (by decide)) hc]

theorem Items.type_eq' {items : Items} {i : ItemId} (h : i < items.size) :
    Items.type items i = items[i]!.type := by
  simp [Items.type, getElem!_pos, h]

theorem Items.ch_eq' {items : Items} {i : ItemId} (h : i < items.size) :
    Items.ch items i = items[i]!.ch := by
  simp [Items.ch, getElem!_pos, h]

theorem Items.vs_eq' {items : Items} {i : ItemId} (h : i < items.size) :
    Items.vs items i = items[i]!.vs := by
  simp [Items.vs, getElem!_pos, h]

/-- A child `c ≥ 1 + nv` of an item is not a V item (it is a Q or a node). -/
theorem child_not_V {g : Graph} {items : Items} (hwf : Items.WF g items) {p c : ItemId}
    (hpc : Items.IsParent items p c) (hc : 1 + g.nv ≤ c) : Items.type items c ≠ .V := by
  have hlt := hwf.tree.ch_lt p c hpc
  by_cases h : c < 1 + g.nv + g.ne
  · have := hwf.tree.edge (c - (1 + g.nv)) (by iomega)
    have he : edgeItem g (c - (1 + g.nv)) = c := by unfold edgeItem; iomega
    rw [he] at this; rw [this]; decide
  · have := hwf.tree.node c (by iomega) hlt
    intro hv; rw [hv] at this; simp at this

/-- Among the children of an item, `c < 1 + nv` and `type c = V` agree. -/
theorem child_lt_iff_V {g : Graph} {items : Items} (hwf : Items.WF g items) {p c : ItemId}
    (hpc : Items.IsParent items p c) : c < 1 + g.nv ↔ Items.type items c = .V := by
  constructor
  · intro h
    have hc0 : c ≠ 0 := by
      rintro rfl; exact hwf.tree.root_no_parent p hpc
    have := hwf.tree.vert (c - 1) (by iomega)
    have hv : vertItem (c - 1) = c := by unfold vertItem; iomega
    rwa [hv] at this
  · intro h
    by_contra hc
    exact child_not_V hwf hpc (by iomega) h

/-- `nvList` lists the first endpoint, the V children (`ch` order), the second endpoint; it is
`Nodup` (`Endpoints.nv_nodup`). -/
theorem nvList_nodup {g : Graph} {items : Items} (hwf : Items.WF g items) {i : ItemId}
    (hi : i < items.size) : (Items.nvList g items i).Nodup := by
  have h := hwf.endpoints.nv_nodup i hi
  unfold Items.nvList
  have hf : ((Items.ch items i).filter fun c => Items.type items c = NodeType.V) =
      (Items.ch items i).filter (· < 1 + g.nv) := by
    apply List.filter_congr
    intro c hc
    simp only [decide_eq_decide]
    exact (child_lt_iff_V hwf hc).symm
  rwa [hf] at h

/-- `nodePlanar_sound_R` from the two records. -/
theorem nodePlanar_sound_R_of (g : Graph) (w : PlanarWalkState) (T : PlanarSpqrTree)
    (hwf : Items.WF g w.base.items) (hpf : PlanarFinish g w) (n : Nat)
    (hrec : RelabelNodeR g w T n) :
    IsPlanarEmbedding (T.localSkeleton n) (T.toSpqrTree.nVerts n) (T.nodeRot n) := by
  obtain ⟨it, m, children, pos, hge, hlt, hR, hnp, hch, hpos, hnv, hne, hsk, hsz, hrows⟩ := hrec
  set items := w.base.items with hitems
  set nvSt := (T.toSpqrTree.nvRange n).1 with hnvSt
  set nvEn := (T.toSpqrTree.nvRange n).2 with hnvEn
  set neSt := (T.toSpqrTree.neRange n).1 with hneSt
  set neEn := (T.toSpqrTree.neRange n).2 with hneEn
  set edgeVes := (children.filter (· ≥ 1 + g.nv)).map (· - (1 + g.nv)) with hedgeVes
  set ves := capVe g :: edgeVes with hves
  set Q := flipped g items[it]!.ch w.aux.itemFlips[it]! (capLinked g w.aux.qem m) with hQ
  set P := nodeSkel g items it with hP
  set nvl := Items.nvList g items it with hnvl
  have hRt : Items.type items it = .R := by rw [Items.type_eq' hlt]; exact hR
  have hcht : Items.ch items it = items[it]!.ch := Items.ch_eq' hlt
  have hvst : Items.vs items it = items[it]!.vs := Items.vs_eq' hlt
  -- PlanarFinish at `it`
  have hit : 1 + g.nv + g.ne + (it - (1 + g.nv + g.ne)) = it := by iomega
  have hty : items[1 + g.nv + g.ne + (it - (1 + g.nv + g.ne))]!.type ∈ [NodeType.S, .P, .R] := by
    rw [hit, hR]; simp
  have hcl := hpf.closed _ m hnp hty
  have hemb := hpf.embedded _ m hnp hty
  have hchlt := hpf.ch_lt _ m hnp hty
  simp only [hit] at hcl hemb hchlt
  change ∀ ve ∈ P.ves, ∀ z, z < 4 → ∃ o, Q[4 * ve + z]! = some o ∧ QE.edge o ∈ P.ves at hcl
  change IsPlanarEmbedding P.es g.nv (readRot P Q) at hemb
  -- children is a permutation of ch
  have hperm_ch : items[it]!.ch.Perm children := by
    rw [hch, Items.ordered, if_neg (by rw [hRt]; simp), hcht]
    exact (List.mergeSort_perm _ _).symm
  have hperm : P.ves.Perm ves := by
    show (capVe g :: (items[it]!.ch.filter (· ≥ 1 + g.nv)).map (· - (1 + g.nv))).Perm
      (capVe g :: edgeVes)
    exact List.Perm.cons _ ((hperm_ch.filter _).map _)
  have hnd : P.ves.Nodup := by
    show (capVe g :: (items[it]!.ch.filter (· ≥ 1 + g.nv)).map (· - (1 + g.nv))).Nodup
    refine List.nodup_cons.2 ⟨?_, ?_⟩
    · intro hmem
      obtain ⟨c, hc, hce⟩ := List.mem_map.1 hmem
      have hc' := List.mem_filter.1 hc
      have := hchlt c hc'.1
      simp only [ge_iff_le, decide_eq_true_eq] at hc'
      unfold capVe at hce
      iomega
    · apply List.Nodup.map_on
      · intro x hx y hy hxy
        have hx' := List.mem_filter.1 hx
        have hy' := List.mem_filter.1 hy
        simp only [ge_iff_le, decide_eq_true_eq] at hx' hy'
        iomega
      · have := hwf.tree.ch_nodup it
        rw [hcht] at this
        exact this.filter _
  have h' := readRot_perm hemb hperm hnd hcl
  set P' : Piece := {P with ves := ves} with hP'
  have hves' : P'.ves = ves := rfl
  -- the rotation system of the node
  have hrot : T.nodeRot n = readRot P' Q := by
    have hlen : neEn - neSt = ves.length := by rw [hne]; simp [ves]; iomega
    cases hT : T.nodeRot n with
    | mk arr =>
    cases hr : readRot P' Q with
    | mk arr' =>
    congr
    apply Array.ext'
    apply List.ext_getElem?
    intro l
    rw [Array.getElem?_toList, Array.getElem?_toList]
    have e1 : (⟨arr⟩ : RotationSystem).rotAdj = arr := rfl
    have e2 : (⟨arr'⟩ : RotationSystem).rotAdj = arr' := rfl
    rw [← e1, ← e2, ← hT, ← hr]
    simp only [PlanarSpqrTree.nodeRot, readRot, Array.getElem?_map, Array.getElem?_extract,
      Array.getElem?_range]
    rw [← hneSt, ← hneEn, Nat.min_eq_left hsz, show 4 * neEn - 4 * neSt = 4 * ves.length by iomega]
    by_cases hl : l < 4 * ves.length
    · rw [if_pos hl, if_pos hl]
      have hb : 4 * neSt + l < T.neRotAdj.size := by iomega
      have hsome : T.neRotAdj[4 * neSt + l]? = some T.neRotAdj[4 * neSt + l]! := by
        rw [getElem!_pos T.neRotAdj (4 * neSt + l) hb, Array.getElem?_eq_getElem hb]
      rw [hsome]
      have hvk : ves[l / 4]! ∈ ves := by
        rw [List.getElem!_eq_getElem?_getD, List.getElem?_eq_getElem (by iomega)]
        exact List.getElem_mem _
      obtain ⟨o, ho, hoe⟩ := hcl _ (hperm.mem_iff.2 hvk) (l % 4) (Nat.mod_lt _ (by decide))
      rw [(hrows l hl).2 o ho (hperm.mem_iff.1 hoe), hves']
      simp only [Option.map_some]
      rw [ho]
      simp only [Option.map_some, Option.some.injEq]
      iomega
    · rw [if_neg hl, if_neg hl]; rfl
  -- the vertex relabelling
  have hvs : ∃ u v, items[it]!.vs = (some u, some v) := by
    have := hwf.endpoints.vs_shape it hlt
    rw [hRt] at this
    simpa [hvst] using this
  obtain ⟨u, v, huv⟩ := hvs
  have hnvl_eq : nvl = (some u).toList ++ ((items[it]!.ch.filter (· < 1 + g.nv)).map (· - 1)) ++
      (some v).toList := by
    rw [hnvl, Items.nvList, hvst, hcht, huv]
  have hu_mem : u ∈ nvl := by rw [hnvl_eq]; simp
  have hv_mem : v ∈ nvl := by rw [hnvl_eq]; simp
  have hnvl_nd : nvl.Nodup := nvList_nodup hwf hlt
  have hchild_ends : ∀ c ∈ children.filter (· ≥ 1 + g.nv),
      ∃ a b, items[c]!.vs = (some a, some b) ∧ a ∈ nvl ∧ b ∈ nvl ∧ c < items.size ∧
        1 + g.nv ≤ c ∧ c < 1 + g.nv + 2 * g.ne := by
    intro c hc
    have hc' := List.mem_filter.1 hc
    simp only [ge_iff_le, decide_eq_true_eq] at hc'
    have hcch : c ∈ items[it]!.ch := hperm_ch.mem_iff.2 hc'.1
    have hpar : Items.IsParent items it c := by
      show c ∈ Items.ch items it; rw [hcht]; exact hcch
    have hcs := hwf.tree.ch_lt it c hpar
    have hnV := child_not_V hwf hpar hc'.2
    obtain ⟨a, b, hab⟩ := (hwf.shapes.r_shape it hlt hRt).2.2.2.2.2 c hpar hnV
    have hvsc : Items.vs items c = items[c]!.vs := Items.vs_eq' hcs
    have hin : ∀ x, (Items.vs items c).1 = some x ∨ (Items.vs items c).2 = some x → x ∈ nvl := by
      intro x hx
      have hx' := hwf.endpoints.child_vs_in_parent it c hpar (by rw [hRt]; simp)
        hnV x hx
      have hxlt := hwf.endpoints.vs_lt c x hx
      rw [hvst, huv] at hx'
      rcases hx' with h1 | h1 | h1
      · simp at h1; subst h1; exact hu_mem
      · simp at h1; subst h1; exact hv_mem
      · rw [hnvl_eq]
        simp only [List.mem_append, List.mem_map, List.mem_filter, decide_eq_true_eq]
        refine Or.inl (Or.inr ⟨vertItem x, ⟨?_, ?_⟩, ?_⟩)
        · have h1' : vertItem x ∈ Items.ch items it := h1
          rw [hcht] at h1'; exact h1'
        · unfold vertItem; iomega
        · unfold vertItem; iomega
    refine ⟨a, b, hvsc ▸ hab, hin a (Or.inl (by rw [hab])), hin b (Or.inr (by rw [hab])),
      hcs, hc'.2, hchlt c hcch⟩
  have hidxlt : ∀ x ∈ nvl, pos x - nvSt < nvl.length := by
    intro x hx
    obtain ⟨_, h2⟩ := hpos x hx
    exact (List.getElem?_eq_some_iff.1 h2).1
  have hends_mem : ∀ p ∈ P'.es, p.1 ∈ nvl ∧ p.2 ∈ nvl := by
    intro p hp
    obtain ⟨ve, hve, rfl⟩ := List.mem_map.1 hp
    rw [hves'] at hve
    rcases List.mem_cons.1 hve with rfl | hve
    · show ((if capVe g = capVe g then items[it]!.vs else _).1.getD 0 ∈ nvl) ∧
        ((if capVe g = capVe g then items[it]!.vs else _).2.getD 0 ∈ nvl)
      rw [if_pos rfl, huv]
      exact ⟨hu_mem, hv_mem⟩
    · obtain ⟨c, hc, rfl⟩ := List.mem_map.1 hve
      obtain ⟨a, b, hab, ha, hb, _, hc1, hc2⟩ := hchild_ends c hc
      have hne' : c - (1 + g.nv) ≠ capVe g := by unfold capVe; iomega
      have hc3 : 1 + g.nv + (c - (1 + g.nv)) = c := by iomega
      show ((if c - (1 + g.nv) = capVe g then _ else items[1 + g.nv + (c - (1 + g.nv))]!.vs).1.getD 0 ∈ nvl) ∧
        ((if c - (1 + g.nv) = capVe g then _ else items[1 + g.nv + (c - (1 + g.nv))]!.vs).2.getD 0 ∈ nvl)
      rw [if_neg hne', hc3, hab]
      exact ⟨ha, hb⟩
  have hnV' : T.toSpqrTree.nVerts n = nvl.length := by
    show nvEn - nvSt = nvl.length; iomega
  have hf : ∀ p ∈ P'.es, pos p.1 - nvSt < T.toSpqrTree.nVerts n ∧
      pos p.2 - nvSt < T.toSpqrTree.nVerts n := by
    intro p hp
    rw [hnV']
    exact ⟨hidxlt _ (hends_mem p hp).1, hidxlt _ (hends_mem p hp).2⟩
  have hinj : ∀ x y, HasEdge P'.es x → HasEdge P'.es y → pos x - nvSt = pos y - nvSt → x = y := by
    intro x y hx hy hxy
    have hxm : x ∈ nvl := by
      obtain ⟨p, hp, h⟩ := hx
      rcases h with rfl | rfl
      · exact (hends_mem p hp).1
      · exact (hends_mem p hp).2
    have hym : y ∈ nvl := by
      obtain ⟨p, hp, h⟩ := hy
      rcases h with rfl | rfl
      · exact (hends_mem p hp).1
      · exact (hends_mem p hp).2
    have h1 := (hpos x hxm).2
    have h2 := (hpos y hym).2
    rw [hxy, h2] at h1
    exact (Option.some.inj h1).symm
  have hmap := IsPlanarEmbedding.map (f := fun x => pos x - nvSt) (n' := T.toSpqrTree.nVerts n)
    h'.verts hf hinj h'
  -- the mapped edge list is the local skeleton
  have hpu : pos u = nvSt := by
    obtain ⟨h1, h2⟩ := hpos u hu_mem
    have h0 : nvl[0]? = some u := by rw [hnvl_eq]; simp
    rw [← h0] at h2
    have hlt0 : 0 < nvl.length := List.length_pos_of_mem hu_mem
    have hlt1 := hidxlt u hu_mem
    rw [List.getElem?_eq_getElem hlt1, List.getElem?_eq_getElem hlt0] at h2
    have := (List.Nodup.getElem_inj_iff hnvl_nd).1 (Option.some.inj h2)
    iomega
  have hpv : pos v = nvSt + nvl.length - 1 := by
    obtain ⟨h1, h2⟩ := hpos v hv_mem
    have hL : nvl.length = ((some u).toList ++ ((items[it]!.ch.filter (· < 1 + g.nv)).map (· - 1))).length + 1 := by
      rw [hnvl_eq]; simp
    have h0 : nvl[nvl.length - 1]? = some v := by
      rw [hL, Nat.add_sub_cancel, hnvl_eq]
      exact List.getElem?_concat_length
    rw [← h0] at h2
    have hlt0 : nvl.length - 1 < nvl.length := by iomega
    have hlt1 := hidxlt v hv_mem
    rw [List.getElem?_eq_getElem hlt1, List.getElem?_eq_getElem hlt0] at h2
    have := (List.Nodup.getElem_inj_iff hnvl_nd).1 (Option.some.inj h2)
    iomega
  have hskel : mapEdges (fun x => pos x - nvSt) P'.es = T.localSkeleton n := by
    unfold PlanarSpqrTree.localSkeleton
    rw [hsk, ← hnvSt]
    simp only [mapEdges, Piece.es, hves', List.map_cons, List.map_map, Items.edgeChildren]
    rw [hves, List.map_cons]
    congr 1
    · show (pos ((if capVe g = capVe g then items[it]!.vs else _).1.getD 0) - nvSt,
        pos ((if capVe g = capVe g then items[it]!.vs else _).2.getD 0) - nvSt) = _
      rw [if_pos rfl, huv]
      simp only [Option.getD_some, hpu, hpv, Nat.sub_self]
      congr 1
      iomega
    · rw [hedgeVes, List.map_map]
      apply List.map_congr_left
      intro c hc
      obtain ⟨a, b, hab, _, _, hcs, hc1, hc2⟩ := hchild_ends c hc
      have hne' : c - (1 + g.nv) ≠ capVe g := by unfold capVe; iomega
      have hc3 : 1 + g.nv + (c - (1 + g.nv)) = c := by iomega
      show (pos ((if c - (1 + g.nv) = capVe g then _ else items[1 + g.nv + (c - (1 + g.nv))]!.vs).1.getD 0) - nvSt,
        pos ((if c - (1 + g.nv) = capVe g then _ else items[1 + g.nv + (c - (1 + g.nv))]!.vs).2.getD 0) - nvSt) =
        (pos ((Items.vs items c).1.getD 0) - nvSt, pos ((Items.vs items c).2.getD 0) - nvSt)
      rw [if_neg hne', hc3, Items.vs_eq' hcs]
  rw [hskel, ← hrot] at hmap
  exact hmap

/-- `Items.WF` of the planar walk's items under `g.WF`/`OrderOK` (`walk_items_wf`). -/
theorem planarWalk_items_wf' (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    Items.WF g (g.planarWalk tern (g.dfsForestFast vo eo)).base.items := by
  rw [planarWalk_base, Graph.dfsForestFast_eq]
  exact walk_items_wf g hg tern vo eo hvo heo

/-- Hypothesis-free form of `planarWalk_items_wf'`, the planar layer's convention (cf.
`spqrTree_wf`, whose `Items.WF` ingredient this is). Admitted: open for the same reason as
`spqrTree_wf` (`walk_items_wf` needs `g.WF` and `OrderOK`; PROOF.md §7.6). -/
theorem planarWalk_items_wf (g : Graph) (tern : Bool) (vo eo : List Nat) :
    Items.WF g (g.planarWalk tern (g.dfsForestFast vo eo)).base.items := by
  sorry

end Spqr
