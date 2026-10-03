import Spqr.PlanarEmbedNodeRGlue

/-!
# `nodeFold_capped_R`: assembling the `R` node

The node's piece in gluing order (`pieceR`: corner `V` items, then capped children, as flattened by
`RG.vflat`/`RG.cflat`), the transport of `Piece.loc` from each child into it (`locV_sub`, `locC_sub`,
`locQ_R`, `locW_R`), the correspondence between the generic glued rotation `RG.Hyp.glue` (with its
`tgtV`/`res`/`skt`/`fin` indices) and the executable fold `nodeFoldR` (`agrees_V_R`, `agrees_C_R`,
`unset_R`), and the final certificate `main_R`: the `nodeStep` fold yields a `Piece.Capped` witness
of `pieceBelow g i` at the cap's endpoints, framed outside the piece.
-/

namespace Spqr.PlanarSpqrTree

variable (t : PlanarSpqrTree) {g : Graph} {i : Nat}

/-- Edges below the first `k` corner `V` items. -/
noncomputable def vesV (i : Nat) (s : EmbedState) (k : Nat) : List Nat :=
  (List.range k).flatMap fun v => t.edgesBelow (t.vItem i s v)

/-- Edges below the first `k` capped children. -/
def vesC (i k : Nat) : List Nat :=
  (List.range k).flatMap fun j => t.edgesBelow (t.sC i)[j]!

/-- The node's piece in gluing order. -/
noncomputable def pieceR (g : Graph) (i : Nat) (s : EmbedState) : Piece :=
  ⟨t.vesV i s (t.cornersR i s).length ++ t.vesC i (t.toSpqrTree.nEdges i - 1),
    fun e => g.edges[e]!, g.nv, 0, 0⟩

theorem vesV_succ (s : EmbedState) (k : Nat) :
    t.vesV i s (k + 1) = t.vesV i s k ++ t.edgesBelow (t.vItem i s k) := by
  simp [vesV, List.range_succ]

theorem vesC_succ (k : Nat) : t.vesC i (k + 1) = t.vesC i k ++ t.edgesBelow (t.sC i)[k]! := by
  simp [vesC, List.range_succ]

theorem mem_vesV {s : EmbedState} {k e : Nat} :
    e ∈ t.vesV i s k ↔ ∃ v, v < k ∧ e ∈ t.edgesBelow (t.vItem i s v) := by
  simp [vesV]

theorem mem_vesC {k e : Nat} : e ∈ t.vesC i k ↔ ∃ j, j < k ∧ e ∈ t.edgesBelow (t.sC i)[j]! := by
  simp [vesC]

theorem vesV_split {s : EmbedState} {k v : Nat} (hv : v < k) :
    ∃ R, t.vesV i s k = t.vesV i s v ++ (t.edgesBelow (t.vItem i s v) ++ R) := by
  induction k with
  | zero => omega
  | succ k ih =>
    rcases Nat.lt_succ_iff_lt_or_eq.1 hv with h | rfl
    · obtain ⟨R, hR⟩ := ih h
      exact ⟨R ++ t.edgesBelow (t.vItem i s k), by rw [vesV_succ, hR]; simp [List.append_assoc]⟩
    · exact ⟨[], by rw [vesV_succ]; simp⟩

theorem vesC_split {k j : Nat} (hj : j < k) :
    ∃ R, t.vesC i k = t.vesC i j ++ (t.edgesBelow (t.sC i)[j]! ++ R) := by
  induction k with
  | zero => omega
  | succ k ih =>
    rcases Nat.lt_succ_iff_lt_or_eq.1 hj with h | rfl
    · obtain ⟨R, hR⟩ := ih h
      exact ⟨R ++ t.edgesBelow (t.sC i)[k]!, by rw [vesC_succ, hR]; simp [List.append_assoc]⟩
    · exact ⟨[], by rw [vesC_succ]; simp⟩

theorem vflat_eq (g : Graph) (s : EmbedState) (k : Nat) :
    (t.rgR g i s).vflat k = (t.vesV i s k).map fun e => g.edges[e]! := by
  show (List.range k).flatMap (fun v => (t.pieceBelow g (t.vItem i s v)).es) = _
  simp only [vesV, List.map_flatMap]; rfl

theorem cflat_eq (g : Graph) (s : EmbedState) (k : Nat) :
    (t.rgR g i s).cflat k = (t.vesC i k).map fun e => g.edges[e]! := by
  show (List.range k).flatMap (fun j => (t.pieceBelow g (t.sC i)[j]!).es) = _
  simp only [vesC, List.map_flatMap]; rfl

theorem vIdx_eq (g : Graph) (s : EmbedState) (v : Nat) :
    (t.rgR g i s).vIdx v = 4 * ((t.rgR g i s).m + 1) + 4 * (t.vesV i s v).length := by
  unfold RG.vIdx; rw [vflat_eq, List.length_map]

theorem cIdx_eq (g : Graph) (s : EmbedState) (j : Nat) :
    (t.rgR g i s).cIdx j = 4 * ((t.rgR g i s).m + 1) +
      4 * (t.vesV i s (t.cornersR i s).length).length + 4 * (t.vesC i j).length := by
  unfold RG.cIdx; rw [vIdx_eq, cflat_eq, List.length_map]; rfl

theorem pieceR_es (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hR : t.toSpqrTree.type i = .R)
    (g : Graph) (s : EmbedState) :
    (t.pieceR g i s).es = (t.rgR g i s).vflat (t.rgR g i s).nV ++ (t.rgR g i s).cflat (t.rgR g i s).m := by
  have hm := t.rgR_m hwf hi hR g s
  rw [vflat_eq, cflat_eq, show (t.rgR g i s).m = t.toSpqrTree.nEdges i - 1 by omega]
  simp only [pieceR, Piece.es, List.map_append]; rfl

theorem pieceR_mem {g : Graph} {s : EmbedState} {q : Nat} (h : (t.pieceR g i s).Mem q) :
    (∃ v, v < (t.cornersR i s).length ∧ (t.pieceBelow g (t.vItem i s v)).Mem q) ∨
    ∃ j, j < t.toSpqrTree.nEdges i - 1 ∧ (t.pieceBelow g (t.sC i)[j]!).Mem q := by
  rcases List.mem_append.1 h with h | h
  · exact Or.inl (t.mem_vesV.1 h)
  · exact Or.inr (t.mem_vesC.1 h)

section Loc

variable (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape) (hrep : t.toSpqrTree.Represents g)
  (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size) (hR : t.toSpqrTree.type i = .R)
  (hloc : IsPlanarEmbedding (t.localSkeleton i) (t.toSpqrTree.nVerts i) (t.nodeRot i))
  (hclosed : t.NodeRotClosed) (hcor : t.NodeCorners) (s : EmbedState)
  (h : t.GluedFaces g (i + 1) s)
include hwf hi hR

omit hR in
theorem pieces_disjoint_R (hsh : t.toSpqrTree.ChildShape) {a b : Nat} (ha : a ∈ t.children i)
    (hb : b ∈ t.children i) (hab : a ≠ b) {q : Nat} (hqa : (t.pieceBelow g a).Mem q)
    (hqb : (t.pieceBelow g b).Mem q) : False :=
  List.disjoint_left.1 (t.maximal_pieces_disjoint hwf hsh (t.child_maximal hwf hi ha)
    (t.child_maximal hwf hi hb) hab) hqa hqb

include hrep hsep hloc hclosed hcor h in
theorem vItem_data {v : Nat} (hv : v < (t.cornersR i s).length) :
    t.vItem i s v ∈ t.children i ∧ t.toSpqrTree.type (t.vItem i s v) = .V := by
  have hc := t.cornersR_corner (i := i) (s := s) hv
  rw [← getElem!_pos (t.cornersR i s) v hv] at hc
  obtain ⟨h1, h2, -⟩ := t.corner_V_R hwf hrep.ne hsep hi hR hloc hclosed hcor s h hc
  exact ⟨h1, h2⟩

theorem sC_data {j : Nat} (hj : j < t.toSpqrTree.nEdges i - 1) :
    (t.sC i)[j]! ∈ t.children i ∧ t.toSpqrTree.type (t.sC i)[j]! ≠ .V := by
  have hl := t.sC_length_R hwf hi hR
  have hm : (t.sC i)[j]! ∈ t.sC i := by
    rw [getElem!_pos (t.sC i) j (by omega)]; exact List.getElem_mem _
  have := List.mem_filter.1 hm
  exact ⟨this.1, by simpa using this.2⟩

include hsh hrep hsep hloc hclosed hcor h in
theorem locV_R {v q r : Nat} (hv : v < (t.cornersR i s).length)
    (hq : (t.pieceBelow g (t.vItem i s v)).loc q = some r) :
    (t.pieceR g i s).loc q = some (4 * (t.vesV i s v).length + r) := by
  obtain ⟨R, hR'⟩ := t.vesV_split (i := i) (s := s) hv
  have heq : t.pieceR g i s = ⟨t.vesV i s v ++ (t.edgesBelow (t.vItem i s v) ++
      (R ++ t.vesC i (t.toSpqrTree.nEdges i - 1))), fun e => g.edges[e]!, g.nv, 0, 0⟩ := by
    simp only [pieceR, hR', List.append_assoc]
  rw [heq]
  have h1 := Piece.loc_append_left (P := t.pieceBelow g (t.vItem i s v))
    (R ++ t.vesC i (t.toSpqrTree.nEdges i - 1)) hq
  refine Piece.loc_append_right (P := { t.pieceBelow g (t.vItem i s v) with
    ves := t.edgesBelow (t.vItem i s v) ++ (R ++ t.vesC i (t.toSpqrTree.nEdges i - 1)) })
    (t.vesV i s v) ?_ h1
  intro hmem
  obtain ⟨v', hv', hmem⟩ := t.mem_vesV.1 hmem
  exact t.pieces_disjoint_R hwf hi hsh (t.vItem_data hwf hrep hsep hi hR hloc hclosed hcor s h (by omega)).1
    (t.vItem_data hwf hrep hsep hi hR hloc hclosed hcor s h hv).1
    (fun heq => by have := t.vItem_inj hwf hsep hi hR hloc hclosed hcor s (by omega) hv heq; omega)
    hmem (Piece.mem_of_loc hq)

include hsh hrep hsep hloc hclosed hcor h in
theorem locC_R {j q r : Nat} (hj : j < t.toSpqrTree.nEdges i - 1)
    (hq : (t.pieceBelow g (t.sC i)[j]!).loc q = some r) :
    (t.pieceR g i s).loc q =
      some (4 * (t.vesV i s (t.cornersR i s).length).length + (4 * (t.vesC i j).length + r)) := by
  obtain ⟨R, hR'⟩ := t.vesC_split (i := i) hj
  have heq : t.pieceR g i s = ⟨t.vesV i s (t.cornersR i s).length ++ (t.vesC i j ++
      (t.edgesBelow (t.sC i)[j]! ++ R)), fun e => g.edges[e]!, g.nv, 0, 0⟩ := by
    simp only [pieceR, hR']
  rw [heq]
  have h1 := Piece.loc_append_left (P := t.pieceBelow g (t.sC i)[j]!) R hq
  have h2 := Piece.loc_append_right (P := { t.pieceBelow g (t.sC i)[j]! with
    ves := t.edgesBelow (t.sC i)[j]! ++ R }) (t.vesC i j) ?_ h1
  · refine Piece.loc_append_right (P := { t.pieceBelow g (t.sC i)[j]! with
      ves := t.vesC i j ++ (t.edgesBelow (t.sC i)[j]! ++ R) })
      (t.vesV i s (t.cornersR i s).length) ?_ h2
    intro hmem
    obtain ⟨v, hv, hmem⟩ := t.mem_vesV.1 hmem
    have hd := t.vItem_data hwf hrep hsep hi hR hloc hclosed hcor s h hv
    have hd' := t.sC_data hwf hi hR hj
    exact t.pieces_disjoint_R hwf hi hsh hd.1 hd'.1 (fun heq => hd'.2 (heq ▸ hd.2)) hmem
      (Piece.mem_of_loc hq)
  · intro hmem
    obtain ⟨j', hj', hmem⟩ := t.mem_vesC.1 hmem
    exact t.pieces_disjoint_R hwf hi hsh (t.sC_data hwf hi hR (by omega)).1 (t.sC_data hwf hi hR hj).1
      (fun heq => by have := t.sC_inj_R hwf hi hR (by omega) hj heq; omega) hmem (Piece.mem_of_loc hq)

include hsh hrep hsep hloc hclosed hcor h in
/-- `pieceBelow`'s edges are listed by `pieceR` up to order (empty `V` children aside). -/
theorem pieceR_perm : (t.pieceR g i s).ves.Perm (t.edgesBelow i) := by
  have hne := hrep.ne
  have hl := t.sC_length_R hwf hi hR
  set L : List Nat := (List.range (t.cornersR i s).length).map (t.vItem i s) ++ t.sC i with hL
  have hves : (t.pieceR g i s).ves = L.flatMap t.edgesBelow := by
    show t.vesV i s _ ++ t.vesC i _ = _
    rw [hL, List.flatMap_append, List.flatMap_map]
    congr 1
    unfold vesC
    rw [← hl]
    clear hL hl L
    generalize t.sC i = C
    induction C with
    | nil => rfl
    | cons c C ih =>
      rw [List.length_cons, List.range_succ_eq_map, List.flatMap_cons, List.flatMap_map,
        List.flatMap_cons, ← ih]
      simp
  have hLsub : ∀ c ∈ L, c ∈ t.children i := by
    intro c hc
    rcases List.mem_append.1 hc with hc | hc
    · obtain ⟨v, hv, rfl⟩ := List.mem_map.1 hc
      exact (t.vItem_data hwf hrep hsep hi hR hloc hclosed hcor s h (List.mem_range.1 hv)).1
    · exact (List.mem_filter.1 hc).1
  have hLnd : L.Nodup := by
    rw [hL, List.nodup_append]
    refine ⟨?_, t.sC_nodup hwf hi, ?_⟩
    · refine List.Nodup.map_on ?_ List.nodup_range
      intro v hv v' hv' heq
      exact t.vItem_inj hwf hsep hi hR hloc hclosed hcor s (List.mem_range.1 hv) (List.mem_range.1 hv') heq
    · intro a ha b hb hab
      obtain ⟨v, hv, rfl⟩ := List.mem_map.1 ha
      have h1 := (t.vItem_data hwf hrep hsep hi hR hloc hclosed hcor s h (List.mem_range.1 hv)).2
      have h2 := (List.mem_filter.1 hb).2
      rw [← hab] at h2
      simp [h1] at h2
  set Z := (t.children i).filter (fun c => decide (c ∉ L)) with hZ
  have hperm : (t.children i).Perm (L ++ Z) := by
    rw [List.perm_ext_iff_of_nodup (t.children_nodup hwf hi)]
    · intro c
      rw [List.mem_append, hZ, List.mem_filter]
      constructor
      · intro hc
        by_cases hcL : c ∈ L
        · exact Or.inl hcL
        · exact Or.inr ⟨hc, by simpa using hcL⟩
      · rintro (hc | ⟨hc, -⟩)
        · exact hLsub c hc
        · exact hc
    · rw [List.nodup_append]
      refine ⟨hLnd, (t.children_nodup hwf hi).filter _, ?_⟩
      intro a ha b hb hab
      have := (List.mem_filter.1 hb).2
      rw [← hab] at this
      simp [ha] at this
  have hZnil : Z.flatMap t.edgesBelow = [] := by
    rw [List.flatMap_eq_nil_iff]
    intro c hc
    obtain ⟨hc, hcL⟩ := List.mem_filter.1 hc
    have hcL : c ∉ L := by simpa using hcL
    by_contra hcne
    have htV : t.toSpqrTree.type c = .V := by
      by_contra htV
      apply hcL
      rw [hL, List.mem_append]
      right
      unfold sC
      rw [List.mem_filter]
      exact ⟨hc, by simpa using htV⟩
    obtain ⟨hcs, hpar, -⟩ := t.child_data hwf hi hc
    obtain ⟨nv, d, hnv, hd, hdv⟩ := hsep.v_child_nv i c hi (Or.inr (Or.inr hR))
      (by rw [t.children_eq]; exact hc) htV
    have hnce : ¬ t.toSpqrTree.CapEnd i nv :=
      (hsep.nv_child i nv d hi (Or.inr (Or.inr hR)) hnv hd).1 (by rw [hdv]; exact hpar)
    obtain ⟨-, hw, -, horig⟩ := t.sW_R (g := g) hwf hsep hi hR hnv hnce
    have hsw : t.sW nv = c := by
      unfold sW; rw [getElem!_of_getElem? hd, hdv]
    rw [hsw] at horig
    rcases t.q_lower_boundary g hwf hne hi hc htV horig s h.toGluedUpTo with ⟨-, h0⟩ | ⟨w0, w1, ρ, h0, -, -⟩
    · exact hcne h0
    obtain ⟨ta, ⟨hta1, hta2, hta3, ⟨d', hd', hd'nv⟩, tb, htb, htb1⟩, -, -⟩ :=
      hcor.inner i nv hi (Or.inr (Or.inr hR)) hnv hnce
    have hen := t.neEn_R hwf hi
    set l := ta - 4 * (t.toSpqrTree.neRange i).1 with hl'
    have hta : ta = 4 * (t.toSpqrTree.neRange i).1 + l := by omega
    have hrot := t.rotR hwf hi hR hloc hclosed (l := l) (by omega)
    rw [← hta] at hrot
    have hvof : t.vOf ta = c := by
      unfold vOf
      rw [show QE.edge ta = ta / 4 from rfl, getElem!_of_getElem? hd', hd'nv, getElem!_of_getElem? hd, hdv]
    have hcorner : t.Corner (t.toSpqrTree.neRange i).1 (t.toSpqrTree.nEdges i) s (t.nodeRot i).rot l := by
      refine ⟨by omega, by omega, by omega, ?_, ?_⟩
      · have := getElem!_of_getElem? htb
        rw [hrot.1] at this
        have : tb = 4 * (t.toSpqrTree.neRange i).1 + (t.nodeRot i).rot l := by
          injection this with this; omega
        omega
      · rw [← hta, hvof, h0]; rfl
    have hmem : l ∈ t.cornersR i s := t.mem_cornersR.2 hcorner
    obtain ⟨v, hv, hvl⟩ := List.mem_iff_getElem.1 hmem
    apply hcL
    rw [hL, List.mem_append, List.mem_map]
    refine Or.inl ⟨v, List.mem_range.2 hv, ?_⟩
    unfold vItem
    rw [getElem!_pos (t.cornersR i s) v hv, hvl, ← hta, hvof]
  have hty : t.types[i]! = .R := by rw [← t.type_eq_of_lt i hi]; exact hR
  rw [hves, t.edgesBelow_eq hwf hsh i hi, hty, ite_eq_right (by decide), Option.toList_none, List.nil_append]
  refine ((hperm.flatMap_right t.edgesBelow).trans ?_).symm
  rw [List.flatMap_append, hZnil, List.append_nil]

end Loc

end Spqr.PlanarSpqrTree

namespace Spqr.RG

variable {G : RG}

theorem res_lo {k x : Nat} (hx : 4 ≤ x) (hxk : x < 4 * (k + 1)) :
    G.res k x = G.cIdx (x / 4 - 1) + G.cl (x / 4 - 1) (x % 4) := by
  unfold res; rw [ite_eq_left ⟨hx, hxk⟩]

theorem res_hi {k x : Nat} (hx : 4 * (k + 1) ≤ x) : G.res k x = x := by
  unfold res; rw [ite_eq_right (by omega)]

theorem fin_ge {x : Nat} (hx : 4 ≤ x) : G.fin x = x := by
  unfold fin
  rw [ite_eq_right (by omega), ite_eq_right (by omega), ite_eq_right (by omega), ite_eq_right (by omega)]

namespace Hyp

variable (H : G.Hyp)
include H

theorem tgtV_ta_full {v k : Nat} (hv : v < G.nV) (hvk : v < k) (hk : k ≤ G.nV) :
    G.tgtV k (G.vta v) = G.vIdx v + G.vw1 v := by
  induction k with
  | zero => omega
  | succ k ih =>
    rcases Nat.lt_succ_iff_lt_or_eq.1 hvk with hlt | rfl
    · have hne := H.corner_ne hv (by omega) (by omega : v ≠ k)
      show (if G.vta v = G.vta k then _ else if G.vta v = G.ρ₀.rot (G.vta k) then _
        else G.tgtV k (G.vta v)) = _
      rw [ite_eq_right hne.1, ite_eq_right hne.2]
      exact ih hlt (by omega)
    · show (if G.vta v = G.vta v then _ else _) = _
      rw [ite_eq_left rfl]

theorem tgtV_tb_full {v k : Nat} (hv : v < G.nV) (hvk : v < k) (hk : k ≤ G.nV) :
    G.tgtV k (G.ρ₀.rot (G.vta v)) = G.vIdx v + G.vw0 v := by
  induction k with
  | zero => omega
  | succ k ih =>
    rcases Nat.lt_succ_iff_lt_or_eq.1 hvk with hlt | rfl
    · have hne := H.corner_ne (by omega : k < G.nV) hv (by omega : k ≠ v)
      have hinj : G.ρ₀.rot (G.vta v) ≠ G.ρ₀.rot (G.vta k) := fun heq =>
        hne.1 (RotationSystem.rot_inj H.planar₀.total H.planar₀.involution
          (by rw [H.ρ₀_size]; exact H.ta_lt v hv)
          (by rw [H.ρ₀_size]; exact H.ta_lt k (by omega)) heq).symm
      show (if G.ρ₀.rot (G.vta v) = G.vta k then _
        else if G.ρ₀.rot (G.vta v) = G.ρ₀.rot (G.vta k) then _ else G.tgtV k _) = _
      rw [ite_eq_right (Ne.symm hne.2), ite_eq_right hinj]
      exact ih hlt (by omega)
    · have hne : G.ρ₀.rot (G.vta v) ≠ G.vta v := by
        have := H.deg₀ (G.vta v) (H.ta_lt v hv)
        intro heq; apply this; rw [heq]
      show (if G.ρ₀.rot (G.vta v) = G.vta v then _
        else if G.ρ₀.rot (G.vta v) = G.ρ₀.rot (G.vta v) then _ else _) = _
      rw [ite_eq_right hne, ite_eq_left rfl]

omit H in
theorem tgtV_other {l : Nat} (hl : ∀ v, v < G.nV → l ≠ G.vta v ∧ l ≠ G.ρ₀.rot (G.vta v))
    {k : Nat} (hk : k ≤ G.nV) : G.tgtV k l = G.ρ₀.rot l := by
  induction k with
  | zero => rfl
  | succ k ih =>
    have := hl k (by omega)
    show (if l = G.vta k then _ else if l = G.ρ₀.rot (G.vta k) then _ else G.tgtV k l) = _
    rw [ite_eq_right this.1, ite_eq_right this.2]
    exact ih (by omega)

end Hyp

/-- The conclusion of `Hyp.glue` for a fixed final rotation system `σ`. -/
def GlueOut (G : RG) (σ : RotationSystem) : Prop :=
  IsPlanarEmbedding (G.vflat G.nV ++ G.cflat G.m) G.n σ ∧
  σ.size + 4 * (G.m + 1) = G.cIdx G.m ∧
  (∀ v, v < G.nV → ∀ r, r < (G.vρ v).size → σ.rot (G.vIdx v + r - 4 * (G.m + 1)) =
    G.fin (if r = G.vw0 v then G.res G.m (G.ρ₀.rot (G.vta v)) else if r = G.vw1 v then G.res G.m (G.vta v)
      else G.vIdx v + (G.vρ v).rot r) - 4 * (G.m + 1)) ∧
  (∀ j, j < G.m → ∀ r, r < (G.cρ j).size → σ.rot (G.cIdx j + r - 4 * (G.m + 1)) =
    G.fin (if r = G.cl j 0 then G.skt G.m (4 * (j + 1)) else if r = G.cl j 1 then G.skt G.m (4 * (j + 1) + 1)
      else if r = G.cl j 2 then G.skt G.m (4 * (j + 1) + 2) else if r = G.cl j 3 then G.skt G.m (4 * (j + 1) + 3)
      else G.cIdx j + (G.cρ j).rot r) - 4 * (G.m + 1)) ∧
  σ.rot (G.skt G.m 0 - 4 * (G.m + 1)) = G.skt G.m 1 - 4 * (G.m + 1) ∧
  σ.rot (G.skt G.m 2 - 4 * (G.m + 1)) = G.skt G.m 3 - 4 * (G.m + 1) ∧
  SameOrbit (σ.stepC 3) (G.skt G.m 1 - 4 * (G.m + 1)) (G.skt G.m 3 - 4 * (G.m + 1))

theorem Hyp.glue' (H : G.Hyp) : ∃ σ, G.GlueOut σ := H.glue

end Spqr.RG

namespace Spqr.PlanarSpqrTree

variable (t : PlanarSpqrTree) {g : Graph} {i : Nat}

section Assemble

variable (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape) (hrep : t.toSpqrTree.Represents g)
  (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size) (hR : t.toSpqrTree.type i = .R)
  (hloc : IsPlanarEmbedding (t.localSkeleton i) (t.toSpqrTree.nVerts i) (t.nodeRot i))
  (hclosed : t.NodeRotClosed) (hcor : t.NodeCorners) (s : EmbedState)
  (h : t.GluedFaces g (i + 1) s)
include hwf hi hR

omit hwf hi hR in
theorem corner_iff {l : Nat} :
    t.Corner (t.toSpqrTree.neRange i).1 (t.toSpqrTree.nEdges i) s (t.nodeRot i).rot l ↔
      ∃ v, v < (t.cornersR i s).length ∧ (t.cornersR i s)[v]! = l := by
  rw [← t.mem_cornersR, List.mem_iff_getElem]
  constructor
  · rintro ⟨v, hv, hvl⟩; exact ⟨v, hv, by rw [getElem!_pos _ v hv]; exact hvl⟩
  · rintro ⟨v, hv, hvl⟩; exact ⟨v, hv, by rw [← getElem!_pos _ v hv]; exact hvl⟩

include hrep hloc in
theorem T_facts {l : Nat} (hl : l < 4 * t.toSpqrTree.nEdges i) :
    (t.nodeRot i).rot l < 4 * t.toSpqrTree.nEdges i ∧
    (t.nodeRot i).rot ((t.nodeRot i).rot l) = l ∧
    (t.nodeRot i).rot l / 4 ≠ l / 4 ∧ (t.nodeRot i).rot l % 2 ≠ l % 2 := by
  have hs : (t.nodeRot i).size = 4 * t.toSpqrTree.nEdges i := by
    rw [hloc.size, t.localSkeleton_length]
  have hl' : l < (t.nodeRot i).size := by omega
  refine ⟨by rw [← hs]; exact RotationSystem.rot_lt hloc.total hloc.involution hl',
    RotationSystem.rot_rot hloc.total hloc.involution hl', t.deg_R hwf hrep hi hR hloc hl, ?_⟩
  exact hloc.opposite_dir l hl' _ (RotationSystem.get_eq_rot hloc.total hl')


omit hR in
include hsh in
theorem mem_pieceBelow_of_child {c q : Nat} (hc : c ∈ t.children i) (hq : (t.pieceBelow g c).Mem q) :
    (t.pieceBelow g i).Mem q := by
  show QE.edge q ∈ t.edgesBelow i
  rw [t.edgesBelow_eq hwf hsh i hi]
  exact List.mem_append_right _ (List.mem_flatMap.2 ⟨c, hc, hq⟩)

include hsh hrep hsep hloc hclosed hcor h in
theorem locV_sub {v q r : Nat} (hv : v < (t.cornersR i s).length)
    (hq : (t.pieceBelow g (t.vItem i s v)).loc q = some r) :
    (t.pieceR g i s).loc q = some ((t.rgR g i s).vIdx v + r - 4 * ((t.rgR g i s).m + 1)) := by
  rw [t.locV_R hwf hsh hrep hsep hi hR hloc hclosed hcor s h hv hq, t.vIdx_eq]
  congr 1; omega

include hsh hrep hsep hloc hclosed hcor h in
theorem locC_sub {j q r : Nat} (hj : j < t.toSpqrTree.nEdges i - 1)
    (hq : (t.pieceBelow g (t.sC i)[j]!).loc q = some r) :
    (t.pieceR g i s).loc q = some ((t.rgR g i s).cIdx j + r - 4 * ((t.rgR g i s).m + 1)) := by
  rw [t.locC_R hwf hsh hrep hsep hi hR hloc hclosed hcor s h hj hq, t.cIdx_eq]
  congr 1; omega

include hsh hrep hsep hloc hclosed hcor h in
theorem vert_V_sub {v q r : Nat} (hv : v < (t.cornersR i s).length)
    (hq : (t.pieceBelow g (t.vItem i s v)).loc q = some r) :
    QE.vert (t.pieceR g i s).es ((t.rgR g i s).vIdx v + r - 4 * ((t.rgR g i s).m + 1)) =
      QE.vert (t.pieceBelow g (t.vItem i s v)).es r := by
  rw [Piece.vert_loc (t.locV_sub hwf hsh hrep hsep hi hR hloc hclosed hcor s h hv hq), Piece.vert_loc hq]
  rfl

include hsh hrep hsep hloc hclosed hcor h in
theorem vert_C_sub {j q r : Nat} (hj : j < t.toSpqrTree.nEdges i - 1)
    (hq : (t.pieceBelow g (t.sC i)[j]!).loc q = some r) :
    QE.vert (t.pieceR g i s).es ((t.rgR g i s).cIdx j + r - 4 * ((t.rgR g i s).m + 1)) =
      QE.vert (t.pieceBelow g (t.sC i)[j]!).es r := by
  rw [Piece.vert_loc (t.locC_sub hwf hsh hrep hsep hi hR hloc hclosed hcor s h hj hq), Piece.vert_loc hq]
  rfl

include hrep hsep h in
theorem memQ_R {l : Nat} (hl4 : 4 ≤ l) (hl : l < 4 * t.toSpqrTree.nEdges i) :
    (t.pieceBelow g (t.sC i)[l / 4 - 1]!).Mem (t.rQ i s l) ∧ s.rotAdj[t.rQ i s l]? = some none := by
  have := t.slot_R hwf hrep.ne hsep hi hR s h (j := l / 4 - 1) (r := l % 4) (by omega) (by omega)
  rwa [show 4 * (l / 4 - 1 + 1) + l % 4 = l by omega] at this

include hsh hrep hsep hloc hclosed hcor h in
theorem locQ_R {l : Nat} (hl4 : 4 ≤ l) (hl : l < 4 * t.toSpqrTree.nEdges i) :
    (t.pieceR g i s).loc (t.rQ i s l) =
      some ((t.rgR g i s).res (t.rgR g i s).m l - 4 * ((t.rgR g i s).m + 1)) := by
  have hm := t.rgR_m hwf hi hR g s
  have hc := t.rgR_cl hwf hi hR g s hrep.ne hsep h (j := l / 4 - 1) (k := l % 4) (by omega) (by omega)
  rw [show 4 * (l / 4 - 1 + 1) + l % 4 = l by omega] at hc
  rw [t.locC_sub hwf hsh hrep hsep hi hR hloc hclosed hcor s h (by omega) hc, RG.res_lo hl4 (by omega)]

include hrep hsep hloc hclosed hcor h in
theorem memW_R {v : Nat} (hv : v < (t.cornersR i s).length) :
    (t.pieceBelow g (t.vItem i s v)).Mem (t.w0 (t.toSpqrTree.neRange i).1 s (t.cornersR i s)[v]!) ∧
    (t.pieceBelow g (t.vItem i s v)).Mem (t.w1 (t.toSpqrTree.neRange i).1 s (t.cornersR i s)[v]!) :=
  ⟨Piece.mem_of_loc (t.rgR_vw hwf hi hR g s hrep.ne hsep h hloc hclosed hcor hv).1,
   Piece.mem_of_loc (t.rgR_vw hwf hi hR g s hrep.ne hsep h hloc hclosed hcor hv).2⟩

include hsh hrep hsep hloc hclosed hcor h in
theorem locW_R {v : Nat} (hv : v < (t.cornersR i s).length) :
    (t.pieceR g i s).loc (t.w0 (t.toSpqrTree.neRange i).1 s (t.cornersR i s)[v]!) =
      some ((t.rgR g i s).vIdx v + (t.rgR g i s).vw0 v - 4 * ((t.rgR g i s).m + 1)) ∧
    (t.pieceR g i s).loc (t.w1 (t.toSpqrTree.neRange i).1 s (t.cornersR i s)[v]!) =
      some ((t.rgR g i s).vIdx v + (t.rgR g i s).vw1 v - 4 * ((t.rgR g i s).m + 1)) := by
  obtain ⟨h0, h1⟩ := t.rgR_vw hwf hi hR g s hrep.ne hsep h hloc hclosed hcor hv
  exact ⟨t.locV_sub hwf hsh hrep hsep hi hR hloc hclosed hcor s h hv h0,
    t.locV_sub hwf hsh hrep hsep hi hR hloc hclosed hcor s h hv h1⟩

include hsh hrep hsep hloc hclosed hcor h in
theorem disj_VC {v j q : Nat} (hv : v < (t.cornersR i s).length) (hj : j < t.toSpqrTree.nEdges i - 1)
    (hqv : (t.pieceBelow g (t.vItem i s v)).Mem q) (hqc : (t.pieceBelow g (t.sC i)[j]!).Mem q) : False := by
  have hd := t.vItem_data hwf hrep hsep hi hR hloc hclosed hcor s h hv
  have hd' := t.sC_data hwf hi hR hj
  exact t.pieces_disjoint_R hwf hi hsh hd.1 hd'.1 (fun heq => hd'.2 (heq ▸ hd.2)) hqv hqc

include hsh hrep hsep hloc hclosed hcor h in
theorem disj_VV {v v' q : Nat} (hv : v < (t.cornersR i s).length) (hv' : v' < (t.cornersR i s).length)
    (hne : v ≠ v') (hqv : (t.pieceBelow g (t.vItem i s v)).Mem q)
    (hqv' : (t.pieceBelow g (t.vItem i s v')).Mem q) : False :=
  t.pieces_disjoint_R hwf hi hsh (t.vItem_data hwf hrep hsep hi hR hloc hclosed hcor s h hv).1
    (t.vItem_data hwf hrep hsep hi hR hloc hclosed hcor s h hv').1
    (fun heq => hne (t.vItem_inj hwf hsep hi hR hloc hclosed hcor s hv hv' heq)) hqv hqv'

include hsh in
theorem disj_CC {j j' q : Nat} (hj : j < t.toSpqrTree.nEdges i - 1) (hj' : j' < t.toSpqrTree.nEdges i - 1)
    (hne : j ≠ j') (hqj : (t.pieceBelow g (t.sC i)[j]!).Mem q)
    (hqj' : (t.pieceBelow g (t.sC i)[j']!).Mem q) : False :=
  t.pieces_disjoint_R hwf hi hsh (t.sC_data hwf hi hR hj).1 (t.sC_data hwf hi hR hj').1
    (fun heq => hne (t.sC_inj_R hwf hi hR hj hj' heq)) hqj hqj'

include hsh hrep hsep hloc hclosed hcor h in
theorem notQ_of_V {v q : Nat} (hv : v < (t.cornersR i s).length)
    (hq : (t.pieceBelow g (t.vItem i s v)).Mem q) :
    ∀ l, 4 ≤ l → l < 4 * t.toSpqrTree.nEdges i → q ≠ t.rQ i s l := by
  intro l hl4 hl heq
  exact t.disj_VC hwf hsh hrep hsep hi hR hloc hclosed hcor s h hv (j := l / 4 - 1) (by omega) hq
    (by rw [heq]; exact (t.memQ_R hwf hrep hsep hi hR s h hl4 hl).1)

include hsh hrep hsep hloc hclosed hcor h in
theorem notW_of_C {j q : Nat} (hj : j < t.toSpqrTree.nEdges i - 1)
    (hq : (t.pieceBelow g (t.sC i)[j]!).Mem q) :
    ∀ l, t.Corner (t.toSpqrTree.neRange i).1 (t.toSpqrTree.nEdges i) s (t.nodeRot i).rot l →
      q ≠ t.w0 (t.toSpqrTree.neRange i).1 s l ∧ q ≠ t.w1 (t.toSpqrTree.neRange i).1 s l := by
  intro l hc
  obtain ⟨v, hv, rfl⟩ := (t.corner_iff (s := s)).1 hc
  obtain ⟨m0, m1⟩ := t.memW_R hwf hrep hsep hi hR hloc hclosed hcor s h hv
  exact ⟨fun heq => t.disj_VC hwf hsh hrep hsep hi hR hloc hclosed hcor s h hv hj (by rw [heq]; exact m0) hq,
    fun heq => t.disj_VC hwf hsh hrep hsep hi hR hloc hclosed hcor s h hv hj (by rw [heq]; exact m1) hq⟩

include hsh hrep hsep hloc hclosed hcor h in
theorem notW_of_V {v q : Nat} (hv : v < (t.cornersR i s).length)
    (hq : (t.pieceBelow g (t.vItem i s v)).Mem q)
    (h0 : q ≠ t.w0 (t.toSpqrTree.neRange i).1 s (t.cornersR i s)[v]!)
    (h1 : q ≠ t.w1 (t.toSpqrTree.neRange i).1 s (t.cornersR i s)[v]!) :
    ∀ l, t.Corner (t.toSpqrTree.neRange i).1 (t.toSpqrTree.nEdges i) s (t.nodeRot i).rot l →
      q ≠ t.w0 (t.toSpqrTree.neRange i).1 s l ∧ q ≠ t.w1 (t.toSpqrTree.neRange i).1 s l := by
  intro l hc
  obtain ⟨v', hv', rfl⟩ := (t.corner_iff (s := s)).1 hc
  by_cases hvv : v' = v
  · subst hvv; exact ⟨h0, h1⟩
  · obtain ⟨m0, m1⟩ := t.memW_R hwf hrep hsep hi hR hloc hclosed hcor s h hv'
    exact ⟨fun heq => t.disj_VV hwf hsh hrep hsep hi hR hloc hclosed hcor s h hv hv' (Ne.symm hvv) hq
        (by rw [heq]; exact m0),
      fun heq => t.disj_VV hwf hsh hrep hsep hi hR hloc hclosed hcor s h hv hv' (Ne.symm hvv) hq
        (by rw [heq]; exact m1)⟩

include hsh hrep hsep h in
theorem notQ_of_C {j q : Nat} (hj : j < t.toSpqrTree.nEdges i - 1)
    (hq : (t.pieceBelow g (t.sC i)[j]!).Mem q) (hk : ∀ k, k < 4 → q ≠ t.rQ i s (4 * (j + 1) + k)) :
    ∀ l, 4 ≤ l → l < 4 * t.toSpqrTree.nEdges i → q ≠ t.rQ i s l := by
  intro l hl4 hl heq
  by_cases hjl : l / 4 - 1 = j
  · exact hk (l % 4) (by omega) (by rw [show 4 * (j + 1) + l % 4 = l by omega]; exact heq)
  · exact t.disj_CC hwf hsh hi hR hj (by omega) (Ne.symm hjl) hq
      (by rw [heq]; exact (t.memQ_R hwf hrep hsep hi hR s h hl4 hl).1)

include hsh hrep hsep hloc hclosed hcor h in
theorem vert_slot_R {l : Nat} (hl4 : 4 ≤ l) (hl : l < 4 * t.toSpqrTree.nEdges i) :
    QE.vert (t.pieceR g i s).es ((t.rgR g i s).res (t.rgR g i s).m l - 4 * ((t.rgR g i s).m + 1)) =
      QE.vert (t.skR i) l := by
  have hm := t.rgR_m hwf hi hR g s
  have hC := t.rgR_cρ hwf hi hR g s hrep.ne hsep h (j := l / 4 - 1) (by omega)
  have hcl := t.rgR_cl hwf hi hR g s hrep.ne hsep h (j := l / 4 - 1) (k := l % 4) (by omega) (by omega)
  have hskE := t.rgR_skE hwf hi g s (j := l / 4 - 1) (by omega)
  rw [RG.res_lo hl4 (by omega), t.vert_C_sub hwf hsh hrep hsep hi hR hloc hclosed hcor s h (by omega) hcl,
    t.vert_skR hwf hi hl, show (t.toSpqrTree.neRange i).1 + l / 4 = (t.toSpqrTree.neRange i).1 + (l / 4 - 1 + 1) by omega]
  obtain ⟨⟨l0, e0, v0⟩, ⟨l1, e1, v1⟩, ⟨l2, e2, v2⟩, ⟨l3, e3, v3⟩⟩ := hC.slots_vert
  rw [hskE] at v0 v1 v2 v3
  rcases (by omega : l % 4 = 0 ∨ l % 4 = 1 ∨ l % 4 = 2 ∨ l % 4 = 3) with h4 | h4 | h4 | h4
  · rw [h4] at hcl ⊢; rw [Nat.add_zero] at hcl
    obtain rfl := Option.some.inj (e0.symm.trans hcl)
    rw [ite_eq_left (by omega)]; exact v0
  · rw [h4] at hcl ⊢
    obtain rfl := Option.some.inj (e1.symm.trans hcl)
    rw [ite_eq_left (by omega)]; exact v1
  · rw [h4] at hcl ⊢
    obtain rfl := Option.some.inj (e2.symm.trans hcl)
    rw [ite_eq_right (by omega)]; exact v2
  · rw [h4] at hcl ⊢
    obtain rfl := Option.some.inj (e3.symm.trans hcl)
    rw [ite_eq_right (by omega)]; exact v3

include hrep hsep h in
theorem dir_slot_R {l : Nat} (hl4 : 4 ≤ l) (hl : l < 4 * t.toSpqrTree.nEdges i) :
    t.rQ i s l % 2 = l % 2 := by
  obtain ⟨-, -, -, -, ρ, hC⟩ := t.child_cert_R hwf hrep.ne hsep hi hR (j := l / 4 - 1) (by omega) s h
  have d0 := hC.dir0; have d1 := hC.dir1; have d2 := hC.dir2; have d3 := hC.dir3
  have hl' : l = 4 * (l / 4 - 1 + 1) + l % 4 := by omega
  rcases (by omega : l % 4 = 0 ∨ l % 4 = 1 ∨ l % 4 = 2 ∨ l % 4 = 3) with h4 | h4 | h4 | h4
  · rw [h4, Nat.add_zero] at hl'; rw [hl']; omega
  · rw [h4] at hl'; rw [hl']; omega
  · rw [h4] at hl'; rw [hl']; omega
  · rw [h4] at hl'; rw [hl']; omega

include hrep hloc in
theorem capT_R {c : Nat} (hc : c < 4) :
    4 ≤ (t.nodeRot i).rot c ∧ (t.nodeRot i).rot c < 4 * t.toSpqrTree.nEdges i := by
  have h6 := (t.shape_R hwf hi hR).2.1
  obtain ⟨h1, -, h2, -⟩ := t.T_facts hwf hrep hi hR hloc (l := c) (by omega)
  exact ⟨by omega, h1⟩


omit hwf hi hR in
theorem mem_pieceR_iff {q : Nat} : (t.pieceR g i s).Mem q ↔
    (∃ v, v < (t.cornersR i s).length ∧ (t.pieceBelow g (t.vItem i s v)).Mem q) ∨
    (∃ j, j < t.toSpqrTree.nEdges i - 1 ∧ (t.pieceBelow g (t.sC i)[j]!).Mem q) := by
  show QE.edge q ∈ _ ++ _ ↔ _
  rw [List.mem_append, t.mem_vesV, t.mem_vesC]
  exact Iff.rfl

omit hwf hi hR in
theorem vesV_length_le {k k' : Nat} (hk : k ≤ k') : (t.vesV i s k).length ≤ (t.vesV i s k').length := by
  induction k' with
  | zero => rw [Nat.le_zero.1 hk]
  | succ k' ih =>
    rcases Nat.of_le_succ hk with hk | hk
    · exact (ih hk).trans (by rw [t.vesV_succ, List.length_append]; omega)
    · rw [hk]

omit hwf hi hR in
theorem vesC_length_le {k k' : Nat} (hk : k ≤ k') : (t.vesC i k).length ≤ (t.vesC i k').length := by
  induction k' with
  | zero => rw [Nat.le_zero.1 hk]
  | succ k' ih =>
    rcases Nat.of_le_succ hk with hk | hk
    · exact (ih hk).trans (by rw [t.vesC_succ, List.length_append]; omega)
    · rw [hk]

include hrep hsep hloc hclosed hcor h in
theorem vρ_size {v : Nat} (hv : v < (t.cornersR i s).length) :
    ((t.rgR g i s).vρ v).size = 4 * (t.edgesBelow (t.vItem i s v)).length := by
  rw [(t.rgR_vρ hwf hi hR g s hrep.ne hsep h hloc hclosed hcor hv).planar.size]
  simp only [Piece.es, List.length_map]
  rfl

include hrep hsep h in
theorem cρ_size {j : Nat} (hj : j < t.toSpqrTree.nEdges i - 1) :
    ((t.rgR g i s).cρ j).size = 4 * (t.edgesBelow (t.sC i)[j]!).length := by
  rw [(t.rgR_cρ hwf hi hR g s hrep.ne hsep h hj).planar.size]
  simp only [Piece.es, List.length_map]
  rfl

include hrep hsep hloc hclosed hcor h in
theorem vIdx_r_lt {σ : RotationSystem} {v r : Nat} (hv : v < (t.cornersR i s).length)
    (hr : r < ((t.rgR g i s).vρ v).size)
    (hσsz : σ.size + 4 * ((t.rgR g i s).m + 1) = (t.rgR g i s).cIdx (t.rgR g i s).m) :
    (t.rgR g i s).vIdx v + r - 4 * ((t.rgR g i s).m + 1) < σ.size := by
  rw [t.vρ_size hwf hrep hsep hi hR hloc hclosed hcor s h hv] at hr
  have h1 := t.vIdx_eq (i := i) g s v
  have h2 := t.cIdx_eq (i := i) g s (t.rgR g i s).m
  have h3 := t.vesV_length_le (i := i) (s := s) (k := v + 1) hv
  rw [t.vesV_succ, List.length_append] at h3
  omega

include hrep hsep h in
theorem cIdx_r_lt {σ : RotationSystem} {j r : Nat} (hj : j < (t.rgR g i s).m)
    (hr : r < ((t.rgR g i s).cρ j).size)
    (hσsz : σ.size + 4 * ((t.rgR g i s).m + 1) = (t.rgR g i s).cIdx (t.rgR g i s).m) :
    (t.rgR g i s).cIdx j + r - 4 * ((t.rgR g i s).m + 1) < σ.size := by
  have hm := t.rgR_m hwf hi hR g s
  rw [t.cρ_size hwf hrep hsep hi hR s h (by omega)] at hr
  have h1 := t.cIdx_eq (i := i) g s j
  have h2 := t.cIdx_eq (i := i) g s (t.rgR g i s).m
  have h3 := t.vesC_length_le (i := i) (k := j + 1) hj
  rw [t.vesC_succ, List.length_append] at h3
  omega

include hrep hsep h in
theorem cl_lt {j k : Nat} (hj : j < t.toSpqrTree.nEdges i - 1) (hk : k < 4) :
    (t.rgR g i s).cl j k < ((t.rgR g i s).cρ j).size := by
  have hC := t.rgR_cρ hwf hi hR g s hrep.ne hsep h hj
  have hcl := t.rgR_cl hwf hi hR g s hrep.ne hsep h hj hk
  obtain ⟨l0, l1, e0, e1, g0, -⟩ := hC.pair0
  obtain ⟨l2, l3, e2, e3, g2, -⟩ := hC.pair2
  have b0 : l0 < ((t.rgR g i s).cρ j).size := by
    obtain ⟨a, ha, -⟩ := Option.bind_eq_some_iff.1 g0
    exact (Array.getElem?_eq_some_iff.1 ha).1
  have b2 : l2 < ((t.rgR g i s).cρ j).size := by
    obtain ⟨a, ha, -⟩ := Option.bind_eq_some_iff.1 g2
    exact (Array.getElem?_eq_some_iff.1 ha).1
  have b1 := (hC.planar.involution l0 b0 l1 g0).1
  have b3 := (hC.planar.involution l2 b2 l3 g2).1
  rcases (by omega : k = 0 ∨ k = 1 ∨ k = 2 ∨ k = 3) with rfl | rfl | rfl | rfl
  · rw [Nat.add_zero] at hcl; rw [← Option.some.inj (e0.symm.trans hcl)]; exact b0
  · rw [← Option.some.inj (e1.symm.trans hcl)]; exact b1
  · rw [← Option.some.inj (e2.symm.trans hcl)]; exact b2
  · rw [← Option.some.inj (e3.symm.trans hcl)]; exact b3

include hsh hrep hsep hloc hclosed hcor h in
theorem cl_ne {j k k' : Nat} (hj : j < t.toSpqrTree.nEdges i - 1) (hk : k < 4) (hk' : k' < 4) (hkk : k ≠ k') :
    (t.rgR g i s).cl j k ≠ (t.rgR g i s).cl j k' := by
  intro heq
  have hRF := t.rfold_R hwf hsh hrep hsep hi hR hloc hclosed hcor s h
  have h1 := t.rgR_cl hwf hi hR g s hrep.ne hsep h hj hk
  have h2 := t.rgR_cl hwf hi hR g s hrep.ne hsep h hj hk'
  rw [heq] at h1
  have := hRF.Q_inj _ _ (by omega) (by omega) (by omega) (by omega) (Piece.loc_injective h1 h2)
  omega

include hsh hrep hsep hloc hclosed hcor h in
theorem agrees_V_R (σ : RotationSystem) (A : Array (Option Nat)) (hσ : (t.rgR g i s).GlueOut σ)
    (hAw : ∀ l, t.Corner (t.toSpqrTree.neRange i).1 (t.toSpqrTree.nEdges i) s (t.nodeRot i).rot l →
      A[t.w1 (t.toSpqrTree.neRange i).1 s l]? = some (some (t.rQ i s l)) ∧
      A[t.w0 (t.toSpqrTree.neRange i).1 s l]? = some (some (t.rQ i s ((t.nodeRot i).rot l))))
    (hAo : ∀ q, (∀ l, 4 ≤ l → l < 4 * t.toSpqrTree.nEdges i → q ≠ t.rQ i s l) →
      (∀ l, t.Corner (t.toSpqrTree.neRange i).1 (t.toSpqrTree.nEdges i) s (t.nodeRot i).rot l →
        q ≠ t.w0 (t.toSpqrTree.neRange i).1 s l ∧ q ≠ t.w1 (t.toSpqrTree.neRange i).1 s l) →
      A[q]? = s.rotAdj[q]?)
    {v q r : Nat} (hv : v < (t.cornersR i s).length) (hq : (t.pieceBelow g (t.vItem i s v)).Mem q)
    (hr : A[q]? = some (some r)) :
    ∃ lq lr, (t.pieceR g i s).loc q = some lq ∧ (t.pieceR g i s).loc r = some lr ∧ σ.get lq = some lr := by
  obtain ⟨hpl, hsz, hσv, -, -, -, -⟩ := hσ
  have H := t.hyp_R hwf hrep hsep hi hR hloc hclosed hcor s hrep.ne h
  have hm := t.rgR_m hwf hi hR g s
  have hV := t.rgR_vρ hwf hi hR g s hrep.ne hsep h hloc hclosed hcor hv
  have hRF := t.rfold_R hwf hsh hrep hsep hi hR hloc hclosed hcor s h
  have hc := t.cornersR_corner (i := i) (s := s) hv
  rw [← getElem!_pos (t.cornersR i s) v hv] at hc
  have hl4 : 4 ≤ (t.cornersR i s)[v]! := hc.1
  have hlT : 4 ≤ (t.nodeRot i).rot (t.cornersR i s)[v]! := by have := hRF.v_lt _ hc; omega
  have hTlt := (t.T_facts hwf hrep hi hR hloc hc.2.1).1
  obtain ⟨lw0, lw1⟩ := t.locW_R hwf hsh hrep hsep hi hR hloc hclosed hcor s h hv
  obtain ⟨c0, c1⟩ := t.rgR_vw hwf hi hR g s hrep.ne hsep h hloc hclosed hcor hv
  have hw0lt := H.v_lt v hv
  have hw1lt := (hV.planar.involution _ hw0lt _ (H.v_get v hv)).1
  have hw01 : (t.rgR g i s).vw0 v ≠ (t.rgR g i s).vw1 v := by
    have := H.v_even v hv; have := Piece.loc_mod_two c1; have := hV.right_dir; omega
  have hget : ∀ a b, a < σ.size → σ.rot a = b → σ.get a = some b := fun a b ha hab => by
    rw [RotationSystem.get_eq_rot hpl.total ha, hab]
  have hvta : (t.rgR g i s).vta v = (t.cornersR i s)[v]! := rfl
  have hρ₀ : (t.rgR g i s).ρ₀ = t.nodeRot i := rfl
  have hcI := t.cIdx_eq (i := i) g s ((t.nodeRot i).rot (t.cornersR i s)[v]! / 4 - 1)
  have hvI := t.vIdx_eq (i := i) g s v
  have hta_lt := H.ta_lt v hv
  have hta_ge := H.ta_ge v hv
  by_cases h0 : q = t.w0 (t.toSpqrTree.neRange i).1 s (t.cornersR i s)[v]!
  · subst h0
    rw [(hAw _ hc).2] at hr
    obtain rfl := Option.some.inj (Option.some.inj hr)
    refine ⟨_, _, lw0, t.locQ_R hwf hsh hrep hsep hi hR hloc hclosed hcor s h hlT hTlt,
      hget _ _ (t.vIdx_r_lt hwf hrep hsep hi hR hloc hclosed hcor s h hv hw0lt hsz) ?_⟩
    rw [hσv v hv _ hw0lt, ite_eq_left rfl, hρ₀, hvta,
      RG.fin_ge (by rw [RG.res_lo hlT (by omega)]; omega)]
  by_cases h1 : q = t.w1 (t.toSpqrTree.neRange i).1 s (t.cornersR i s)[v]!
  · subst h1
    rw [(hAw _ hc).1] at hr
    obtain rfl := Option.some.inj (Option.some.inj hr)
    refine ⟨_, _, lw1, t.locQ_R hwf hsh hrep hsep hi hR hloc hclosed hcor s h hl4 hc.2.1,
      hget _ _ (t.vIdx_r_lt hwf hrep hsep hi hR hloc hclosed hcor s h hv hw1lt hsz) ?_⟩
    have hcI' := t.cIdx_eq (i := i) g s ((t.cornersR i s)[v]! / 4 - 1)
    rw [hσv v hv _ hw1lt, ite_eq_right (Ne.symm hw01), ite_eq_left rfl, hvta,
      RG.fin_ge (by rw [RG.res_lo hl4 (by omega)]; omega)]
  rw [hAo q (t.notQ_of_V hwf hsh hrep hsep hi hR hloc hclosed hcor s h hv hq)
    (t.notW_of_V hwf hsh hrep hsep hi hR hloc hclosed hcor s h hv hq h0 h1)] at hr
  obtain ⟨lq, lr, hlq, hlr, hg⟩ := hV.agrees q r hq hr
  have hlq_lt : lq < ((t.rgR g i s).vρ v).size := by
    obtain ⟨a, ha, -⟩ := Option.bind_eq_some_iff.1 hg
    exact (Array.getElem?_eq_some_iff.1 ha).1
  have hrot : ((t.rgR g i s).vρ v).rot lq = lr := by unfold RotationSystem.rot; rw [hg]; rfl
  have hne0 : lq ≠ (t.rgR g i s).vw0 v := fun heq => h0 (Piece.loc_injective hlq (by rw [heq]; exact c0))
  have hne1 : lq ≠ (t.rgR g i s).vw1 v := fun heq => h1 (Piece.loc_injective hlq (by rw [heq]; exact c1))
  refine ⟨_, _, t.locV_sub hwf hsh hrep hsep hi hR hloc hclosed hcor s h hv hlq,
    t.locV_sub hwf hsh hrep hsep hi hR hloc hclosed hcor s h hv hlr,
    hget _ _ (t.vIdx_r_lt hwf hrep hsep hi hR hloc hclosed hcor s h hv hlq_lt hsz) ?_⟩
  rw [hσv v hv lq hlq_lt, ite_eq_right hne0, ite_eq_right hne1, hrot, RG.fin_ge (by omega)]


include hsh hrep hsep hloc hclosed hcor h in
theorem agrees_C_R (σ : RotationSystem) (A : Array (Option Nat)) (hσ : (t.rgR g i s).GlueOut σ)
    (hA : ∀ l, 4 ≤ l → l < 4 * t.toSpqrTree.nEdges i → 4 ≤ (t.nodeRot i).rot l →
      A[t.rQ i s l]? = some (some (t.partner (t.toSpqrTree.neRange i).1 (t.toSpqrTree.nEdges i) s
        (t.nodeRot i).rot (t.rQ i s) l)))
    (hAcap : ∀ l, 4 ≤ l → l < 4 * t.toSpqrTree.nEdges i → (t.nodeRot i).rot l < 4 →
      A[t.rQ i s l]? = s.rotAdj[t.rQ i s l]?)
    (hAo : ∀ q, (∀ l, 4 ≤ l → l < 4 * t.toSpqrTree.nEdges i → q ≠ t.rQ i s l) →
      (∀ l, t.Corner (t.toSpqrTree.neRange i).1 (t.toSpqrTree.nEdges i) s (t.nodeRot i).rot l →
        q ≠ t.w0 (t.toSpqrTree.neRange i).1 s l ∧ q ≠ t.w1 (t.toSpqrTree.neRange i).1 s l) →
      A[q]? = s.rotAdj[q]?)
    {j q r : Nat} (hj : j < t.toSpqrTree.nEdges i - 1) (hq : (t.pieceBelow g (t.sC i)[j]!).Mem q)
    (hr : A[q]? = some (some r)) :
    ∃ lq lr, (t.pieceR g i s).loc q = some lq ∧ (t.pieceR g i s).loc r = some lr ∧ σ.get lq = some lr := by
  obtain ⟨hpl, hsz, -, hσc, -, -, -⟩ := hσ
  have H := t.hyp_R hwf hrep hsep hi hR hloc hclosed hcor s hrep.ne h
  have hm := t.rgR_m hwf hi hR g s
  have hC := t.rgR_cρ hwf hi hR g s hrep.ne hsep h hj
  have hget : ∀ a b, a < σ.size → σ.rot a = b → σ.get a = some b := fun a b ha hab => by
    rw [RotationSystem.get_eq_rot hpl.total ha, hab]
  have hjm : j < (t.rgR g i s).m := by omega
  have hcl : ∀ k, k < 4 → (t.pieceBelow g (t.sC i)[j]!).loc (t.rQ i s (4 * (j + 1) + k)) =
      some ((t.rgR g i s).cl j k) := fun k hk => t.rgR_cl hwf hi hR g s hrep.ne hsep h hj hk
  have hcl_lt : ∀ k, k < 4 → (t.rgR g i s).cl j k < ((t.rgR g i s).cρ j).size :=
    fun k hk => t.cl_lt hwf hrep hsep hi hR s h hj hk
  have hcne : ∀ k k', k < 4 → k' < 4 → k ≠ k' → (t.rgR g i s).cl j k ≠ (t.rgR g i s).cl j k' :=
    fun k k' hk hk' hkk => t.cl_ne hwf hsh hrep hsep hi hR hloc hclosed hcor s h hj hk hk' hkk
  have hρ₀ : (t.rgR g i s).ρ₀ = t.nodeRot i := rfl
  have hcI := t.cIdx_eq (i := i) g s j
  by_cases hs : ∃ k, k < 4 ∧ q = t.rQ i s (4 * (j + 1) + k)
  · obtain ⟨k, hk, rfl⟩ := hs
    have hl4 : 4 ≤ 4 * (j + 1) + k := by omega
    have hl' : 4 * (j + 1) + k < 4 * t.toSpqrTree.nEdges i := by omega
    obtain ⟨hTlt, hTT, hTe, -⟩ := t.T_facts hwf hrep hi hR hloc hl'
    by_cases hT : (t.nodeRot i).rot (4 * (j + 1) + k) < 4
    · rw [hAcap _ hl4 hl' hT, (t.memQ_R hwf hrep hsep hi hR s h hl4 hl').2] at hr
      cases hr
    push Not at hT
    rw [hA _ hl4 hl' hT] at hr
    obtain rfl := Option.some.inj (Option.some.inj hr)
    have hlq := t.locQ_R hwf hsh hrep hsep hi hR hloc hclosed hcor s h hl4 hl'
    rw [RG.res_lo hl4 (by omega), show (4 * (j + 1) + k) / 4 - 1 = j by omega,
      show (4 * (j + 1) + k) % 4 = k by omega] at hlq
    have hrot := hσc j hjm _ (hcl_lt k hk)
    have hchain : (if (t.rgR g i s).cl j k = (t.rgR g i s).cl j 0 then
          (t.rgR g i s).skt (t.rgR g i s).m (4 * (j + 1))
        else if (t.rgR g i s).cl j k = (t.rgR g i s).cl j 1 then
          (t.rgR g i s).skt (t.rgR g i s).m (4 * (j + 1) + 1)
        else if (t.rgR g i s).cl j k = (t.rgR g i s).cl j 2 then
          (t.rgR g i s).skt (t.rgR g i s).m (4 * (j + 1) + 2)
        else if (t.rgR g i s).cl j k = (t.rgR g i s).cl j 3 then
          (t.rgR g i s).skt (t.rgR g i s).m (4 * (j + 1) + 3)
        else (t.rgR g i s).cIdx j + ((t.rgR g i s).cρ j).rot ((t.rgR g i s).cl j k)) =
        (t.rgR g i s).skt (t.rgR g i s).m (4 * (j + 1) + k) := by
      rcases (by omega : k = 0 ∨ k = 1 ∨ k = 2 ∨ k = 3) with rfl | rfl | rfl | rfl
      · simp
      · rw [ite_eq_right (hcne 1 0 (by omega) (by omega) (by omega)), ite_eq_left rfl]
      · rw [ite_eq_right (hcne 2 0 (by omega) (by omega) (by omega)),
          ite_eq_right (hcne 2 1 (by omega) (by omega) (by omega)), ite_eq_left rfl]
      · rw [ite_eq_right (hcne 3 0 (by omega) (by omega) (by omega)),
          ite_eq_right (hcne 3 1 (by omega) (by omega) (by omega)),
          ite_eq_right (hcne 3 2 (by omega) (by omega) (by omega)), ite_eq_left rfl]
    rw [hchain] at hrot
    have hskt : ∃ lr, (t.pieceR g i s).loc (t.partner (t.toSpqrTree.neRange i).1 (t.toSpqrTree.nEdges i) s
          (t.nodeRot i).rot (t.rQ i s) (4 * (j + 1) + k)) = some lr ∧
        (t.rgR g i s).fin ((t.rgR g i s).skt (t.rgR g i s).m (4 * (j + 1) + k)) -
          4 * ((t.rgR g i s).m + 1) = lr := by
      unfold partner
      by_cases hc : t.Corner (t.toSpqrTree.neRange i).1 (t.toSpqrTree.nEdges i) s (t.nodeRot i).rot
          (4 * (j + 1) + k)
      · rw [ite_eq_left hc]
        obtain ⟨v, hv, hvl⟩ := (t.corner_iff (s := s)).1 hc
        have hvI := t.vIdx_eq (i := i) g s v
        rw [← hvl]
        refine ⟨_, (t.locW_R hwf hsh hrep hsep hi hR hloc hclosed hcor s h hv).2, ?_⟩
        unfold RG.skt
        rw [show (t.cornersR i s)[v]! = (t.rgR g i s).vta v from rfl, H.tgtV_ta_full (k := (t.rgR g i s).nV) hv hv le_rfl,
          RG.res_hi (by omega), RG.fin_ge (by omega)]
      rw [ite_eq_right hc]
      by_cases hc' : t.Corner (t.toSpqrTree.neRange i).1 (t.toSpqrTree.nEdges i) s (t.nodeRot i).rot
          ((t.nodeRot i).rot (4 * (j + 1) + k))
      · rw [ite_eq_left hc']
        obtain ⟨v, hv, hvl⟩ := (t.corner_iff (s := s)).1 hc'
        have hvI := t.vIdx_eq (i := i) g s v
        rw [← hvl]
        refine ⟨_, (t.locW_R hwf hsh hrep hsep hi hR hloc hclosed hcor s h hv).1, ?_⟩
        have hl_eq : 4 * (j + 1) + k = (t.rgR g i s).ρ₀.rot ((t.rgR g i s).vta v) := by
          show _ = (t.nodeRot i).rot (t.cornersR i s)[v]!
          rw [hvl, hTT]
        unfold RG.skt
        rw [hl_eq, H.tgtV_tb_full (k := (t.rgR g i s).nV) hv hv le_rfl, RG.res_hi (by omega), RG.fin_ge (by omega)]
      rw [ite_eq_right hc']
      refine ⟨_, t.locQ_R hwf hsh hrep hsep hi hR hloc hclosed hcor s h hT hTlt, ?_⟩
      have hoth : ∀ v, v < (t.rgR g i s).nV → 4 * (j + 1) + k ≠ (t.rgR g i s).vta v ∧
          4 * (j + 1) + k ≠ (t.rgR g i s).ρ₀.rot ((t.rgR g i s).vta v) := by
        intro v hv
        have hcv := t.cornersR_corner (i := i) (s := s) hv
        rw [← getElem!_pos (t.cornersR i s) v hv] at hcv
        constructor
        · intro heq
          exact hc (by rw [show 4 * (j + 1) + k = (t.cornersR i s)[v]! from heq]; exact hcv)
        · intro heq
          apply hc'
          have : (t.nodeRot i).rot (4 * (j + 1) + k) = (t.cornersR i s)[v]! := by
            rw [show 4 * (j + 1) + k = (t.nodeRot i).rot (t.cornersR i s)[v]! from heq]
            exact (t.T_facts hwf hrep hi hR hloc hcv.2.1).2.1
          rw [this]; exact hcv
      have hcIT := t.cIdx_eq (i := i) g s ((t.nodeRot i).rot (4 * (j + 1) + k) / 4 - 1)
      unfold RG.skt
      rw [RG.Hyp.tgtV_other hoth le_rfl, hρ₀, RG.fin_ge (by rw [RG.res_lo hT (by omega)]; omega)]
    obtain ⟨lr, hlr, hlr'⟩ := hskt
    refine ⟨_, lr, hlq, hlr, hget _ _ (t.cIdx_r_lt hwf hrep hsep hi hR s h hjm (hcl_lt k hk) hsz) ?_⟩
    rw [hrot]; exact hlr'
  push Not at hs
  rw [hAo q (t.notQ_of_C hwf hsh hrep hsep hi hR s h hj hq hs)
    (t.notW_of_C hwf hsh hrep hsep hi hR hloc hclosed hcor s h hj hq)] at hr
  obtain ⟨lq, lr, hlq, hlr, hg⟩ := hC.agrees q r hq hr
  have hlq_lt : lq < ((t.rgR g i s).cρ j).size := by
    obtain ⟨a, ha, -⟩ := Option.bind_eq_some_iff.1 hg
    exact (Array.getElem?_eq_some_iff.1 ha).1
  have hrot : ((t.rgR g i s).cρ j).rot lq = lr := by unfold RotationSystem.rot; rw [hg]; rfl
  have hne : ∀ k, k < 4 → lq ≠ (t.rgR g i s).cl j k := fun k hk heq =>
    hs k hk (Piece.loc_injective hlq (by rw [heq]; exact hcl k hk))
  refine ⟨_, _, t.locC_sub hwf hsh hrep hsep hi hR hloc hclosed hcor s h hj hlq,
    t.locC_sub hwf hsh hrep hsep hi hR hloc hclosed hcor s h hj hlr,
    hget _ _ (t.cIdx_r_lt hwf hrep hsep hi hR s h hjm hlq_lt hsz) ?_⟩
  rw [hσc j hjm lq hlq_lt, ite_eq_right (hne 0 (by omega)), ite_eq_right (hne 1 (by omega)),
    ite_eq_right (hne 2 (by omega)), ite_eq_right (hne 3 (by omega)), hrot, RG.fin_ge (by omega)]

include hsh hrep hsep hloc hclosed hcor h in
theorem unset_R (A : Array (Option Nat))
    (hA : ∀ l, 4 ≤ l → l < 4 * t.toSpqrTree.nEdges i → 4 ≤ (t.nodeRot i).rot l →
      A[t.rQ i s l]? = some (some (t.partner (t.toSpqrTree.neRange i).1 (t.toSpqrTree.nEdges i) s
        (t.nodeRot i).rot (t.rQ i s) l)))
    (hAcap : ∀ l, 4 ≤ l → l < 4 * t.toSpqrTree.nEdges i → (t.nodeRot i).rot l < 4 →
      A[t.rQ i s l]? = s.rotAdj[t.rQ i s l]?)
    (hAw : ∀ l, t.Corner (t.toSpqrTree.neRange i).1 (t.toSpqrTree.nEdges i) s (t.nodeRot i).rot l →
      A[t.w1 (t.toSpqrTree.neRange i).1 s l]? = some (some (t.rQ i s l)) ∧
      A[t.w0 (t.toSpqrTree.neRange i).1 s l]? = some (some (t.rQ i s ((t.nodeRot i).rot l))))
    (hAo : ∀ q, (∀ l, 4 ≤ l → l < 4 * t.toSpqrTree.nEdges i → q ≠ t.rQ i s l) →
      (∀ l, t.Corner (t.toSpqrTree.neRange i).1 (t.toSpqrTree.nEdges i) s (t.nodeRot i).rot l →
        q ≠ t.w0 (t.toSpqrTree.neRange i).1 s l ∧ q ≠ t.w1 (t.toSpqrTree.neRange i).1 s l) →
      A[q]? = s.rotAdj[q]?) :
    ∀ q, (t.pieceR g i s).Mem q → (A[q]? = some none ↔
      q = t.rQ i s ((t.nodeRot i).rot 1) ∨ q = t.rQ i s ((t.nodeRot i).rot 0) ∨
      q = t.rQ i s ((t.nodeRot i).rot 3) ∨ q = t.rQ i s ((t.nodeRot i).rot 2)) := by
  intro q hq
  have hRF := t.rfold_R hwf hsh hrep hsep hi hR hloc hclosed hcor s h
  have hm := t.rgR_m hwf hi hR g s
  have h6 := (t.shape_R hwf hi hR).2.1
  have hcT : ∀ c, c < 4 → 4 ≤ (t.nodeRot i).rot c ∧ (t.nodeRot i).rot c < 4 * t.toSpqrTree.nEdges i :=
    fun c hc => t.capT_R hwf hrep hi hR hloc hc
  have hmT : ∀ c, c < 4 → (t.pieceBelow g (t.sC i)[(t.nodeRot i).rot c / 4 - 1]!).Mem
      (t.rQ i s ((t.nodeRot i).rot c)) :=
    fun c hc => (t.memQ_R hwf hrep hsep hi hR s h (hcT c hc).1 (hcT c hc).2).1
  have hTT : ∀ c, c < 4 → (t.nodeRot i).rot ((t.nodeRot i).rot c) = c :=
    fun c hc => (t.T_facts hwf hrep hi hR hloc (by omega)).2.1
  have c0 := hcT 0 (by omega); have c1 := hcT 1 (by omega)
  have c2 := hcT 2 (by omega); have c3 := hcT 3 (by omega)
  have t0 := hTT 0 (by omega); have t1 := hTT 1 (by omega)
  have t2 := hTT 2 (by omega); have t3 := hTT 3 (by omega)
  rcases (t.mem_pieceR_iff (s := s)).1 hq with ⟨v, hv, hqv⟩ | ⟨j, hj, hqc⟩
  · have hV := t.rgR_vρ hwf hi hR g s hrep.ne hsep h hloc hclosed hcor hv
    have hc := t.cornersR_corner (i := i) (s := s) hv
    rw [← getElem!_pos (t.cornersR i s) v hv] at hc
    have hr : ¬ A[q]? = some none := by
      by_cases h0 : q = t.w0 (t.toSpqrTree.neRange i).1 s (t.cornersR i s)[v]!
      · rw [h0, (hAw _ hc).2]; intro h; cases h
      by_cases h1 : q = t.w1 (t.toSpqrTree.neRange i).1 s (t.cornersR i s)[v]!
      · rw [h1, (hAw _ hc).1]; intro h; cases h
      rw [hAo q (t.notQ_of_V hwf hsh hrep hsep hi hR hloc hclosed hcor s h hv hqv)
        (t.notW_of_V hwf hsh hrep hsep hi hR hloc hclosed hcor s h hv hqv h0 h1), hV.unset q hqv]
      rintro (h | h)
      · exact h0 h
      · exact h1 h
    refine ⟨fun h => (hr h).elim, ?_⟩
    rintro (rfl | rfl | rfl | rfl)
    · exact (t.disj_VC hwf hsh hrep hsep hi hR hloc hclosed hcor s h hv (by omega) hqv (hmT 1 (by omega))).elim
    · exact (t.disj_VC hwf hsh hrep hsep hi hR hloc hclosed hcor s h hv (by omega) hqv (hmT 0 (by omega))).elim
    · exact (t.disj_VC hwf hsh hrep hsep hi hR hloc hclosed hcor s h hv (by omega) hqv (hmT 3 (by omega))).elim
    · exact (t.disj_VC hwf hsh hrep hsep hi hR hloc hclosed hcor s h hv (by omega) hqv (hmT 2 (by omega))).elim
  · have hC := t.rgR_cρ hwf hi hR g s hrep.ne hsep h hj
    by_cases hs : ∃ k, k < 4 ∧ q = t.rQ i s (4 * (j + 1) + k)
    · obtain ⟨k, hk, rfl⟩ := hs
      have hl4 : 4 ≤ 4 * (j + 1) + k := by omega
      have hl' : 4 * (j + 1) + k < 4 * t.toSpqrTree.nEdges i := by omega
      obtain ⟨hTlt, hTT', -, -⟩ := t.T_facts hwf hrep hi hR hloc hl'
      by_cases hT : (t.nodeRot i).rot (4 * (j + 1) + k) < 4
      · rw [hAcap _ hl4 hl' hT, (t.memQ_R hwf hrep hsep hi hR s h hl4 hl').2]
        refine ⟨fun _ => ?_, fun _ => rfl⟩
        rw [← hTT']
        rcases (by omega : (t.nodeRot i).rot (4 * (j + 1) + k) = 0 ∨ (t.nodeRot i).rot (4 * (j + 1) + k) = 1 ∨
          (t.nodeRot i).rot (4 * (j + 1) + k) = 2 ∨ (t.nodeRot i).rot (4 * (j + 1) + k) = 3) with
          h4 | h4 | h4 | h4 <;> rw [h4] <;> simp
      · push Not at hT
        rw [hA _ hl4 hl' hT]
        refine ⟨(fun h => by cases h), ?_⟩
        rintro (h | h | h | h) <;>
          (have := hRF.Q_inj _ _ hl4 hl' (hcT _ (by omega)).1 (hcT _ (by omega)).2 h; rw [this] at hT; omega)
    · push Not at hs
      rw [hAo q (t.notQ_of_C hwf hsh hrep hsep hi hR s h hj hqc hs)
        (t.notW_of_C hwf hsh hrep hsep hi hR hloc hclosed hcor s h hj hqc), hC.unset q hqc]
      constructor
      · rintro (h | h | h | h)
        · exact (hs 0 (by omega) (by rw [Nat.add_zero]; exact h)).elim
        · exact (hs 1 (by omega) h).elim
        · exact (hs 2 (by omega) h).elim
        · exact (hs 3 (by omega) h).elim
      · rintro (rfl | rfl | rfl | rfl)
        · exact (t.notQ_of_C hwf hsh hrep hsep hi hR s h hj hqc hs _ c1.1 c1.2 rfl).elim
        · exact (t.notQ_of_C hwf hsh hrep hsep hi hR s h hj hqc hs _ c0.1 c0.2 rfl).elim
        · exact (t.notQ_of_C hwf hsh hrep hsep hi hR s h hj hqc hs _ c3.1 c3.2 rfl).elim
        · exact (t.notQ_of_C hwf hsh hrep hsep hi hR s h hj hqc hs _ c2.1 c2.2 rfl).elim


omit hwf hi hR in
theorem fold_eq_R :
    (List.range' (4 * t.neBounds[i]!) (4 * t.neBounds[i + 1]! - 4 * t.neBounds[i]!)).foldl
      (t.nodeStep i t.neBounds[i]!) s =
    (List.range' (4 * (t.toSpqrTree.neRange i).1) (4 * t.toSpqrTree.nEdges i)).foldl
      (t.nodeStep i (t.toSpqrTree.neRange i).1) s := by
  unfold SpqrTree.nEdges
  rw [PlanarRot.neRange_eq]
  dsimp only
  rw [show 4 * t.neBounds[i + 1]! - 4 * t.neBounds[i]! = 4 * (t.neBounds[i + 1]! - t.neBounds[i]!) by omega]

omit hwf hi hR in
theorem array4_ext {r : Array (Option Nat)} {x0 x1 x2 x3 : Option Nat} (hr : r.size = 4)
    (h0 : r[0]? = some x0) (h1 : r[1]? = some x1) (h2 : r[2]? = some x2) (h3 : r[3]? = some x3) :
    r = #[x0, x1, x2, x3] := by
  apply Array.ext (by simp [hr])
  intro k hk hk'
  obtain ⟨_, e0⟩ := Array.getElem?_eq_some_iff.1 h0
  obtain ⟨_, e1⟩ := Array.getElem?_eq_some_iff.1 h1
  obtain ⟨_, e2⟩ := Array.getElem?_eq_some_iff.1 h2
  obtain ⟨_, e3⟩ := Array.getElem?_eq_some_iff.1 h3
  rcases (by omega : k = 0 ∨ k = 1 ∨ k = 2 ∨ k = 3) with rfl | rfl | rfl | rfl
  · simp [e0]
  · simp [e1]
  · simp [e2]
  · simp [e3]

include hsh hrep hsep hloc hclosed hcor h in
theorem main_R {ne : Nat} (hcne : t.toSpqrTree.capNe i = some ne) {p : Nat × Nat}
    (hp : t.toSpqrTree.neOrig ne = some p) :
    let s' := (List.range' (4 * t.neBounds[i]!) (4 * t.neBounds[i + 1]! - 4 * t.neBounds[i]!)).foldl
      (t.nodeStep i t.neBounds[i]!) s
    (∀ q, ¬(t.pieceBelow g i).Mem q → s'.rotAdj[q]? = s.rotAdj[q]?) ∧
    ∃ a0 a1 a2 a3 ρ, s'.outerE[i]? = some #[some a0, some a1, some a2, some a3] ∧
      (t.pieceBelow g i).Capped s'.rotAdj ρ a0 a1 a2 a3 p.1 p.2 := by
  dsimp only
  rw [t.fold_eq_R]
  set s' := (List.range' (4 * (t.toSpqrTree.neRange i).1) (4 * t.toSpqrTree.nEdges i)).foldl
    (t.nodeStep i (t.toSpqrTree.neRange i).1) s with hs'
  have hRF := t.rfold_R hwf hsh hrep hsep hi hR hloc hclosed hcor s h
  obtain ⟨-, -, hA, hAcap, hAw, hAo, r, hrow, hr4, hrk⟩ := t.nodeFoldR hRF
  have H := t.hyp_R hwf hrep hsep hi hR hloc hclosed hcor s hrep.ne h
  obtain ⟨σ, hσ⟩ := H.glue'
  have hσ' := hσ
  obtain ⟨hpl, hsz, -, -, hcap0, hcap2, hface⟩ := hσ'
  have hm := t.rgR_m hwf hi hR g s
  have h6 := (t.shape_R hwf hi hR).2.1
  rw [t.capNe_R hR] at hcne
  cases hcne
  have hsk0 : (t.skR i)[0]? = some p := by
    rw [t.skR_getElem? hwf hi (by omega), Nat.add_zero]
    have hne0 := t.neOrig_R hwf hsep hi hR (k := 0) (by omega)
    rw [Nat.add_zero] at hne0
    rw [hne0] at hp
    exact hp
  have hcT : ∀ c, c < 4 → 4 ≤ (t.nodeRot i).rot c ∧ (t.nodeRot i).rot c < 4 * t.toSpqrTree.nEdges i :=
    fun c hc => t.capT_R hwf hrep hi hR hloc hc
  have hρ₀ : (t.rgR g i s).ρ₀ = t.nodeRot i := rfl
  have hskt : ∀ c, c < 4 → (t.rgR g i s).skt (t.rgR g i s).m c =
      (t.rgR g i s).res (t.rgR g i s).m ((t.nodeRot i).rot c) := by
    intro c hc
    have hoth : ∀ v, v < (t.rgR g i s).nV → c ≠ (t.rgR g i s).vta v ∧
        c ≠ (t.rgR g i s).ρ₀.rot ((t.rgR g i s).vta v) := by
      intro v hv
      have hcv := t.cornersR_corner (i := i) (s := s) hv
      rw [← getElem!_pos (t.cornersR i s) v hv] at hcv
      have := hRF.v_lt _ hcv
      have h4 := hcv.1
      constructor
      · show c ≠ (t.cornersR i s)[v]!; omega
      · show c ≠ (t.nodeRot i).rot (t.cornersR i s)[v]!; omega
    unfold RG.skt
    rw [RG.Hyp.tgtV_other hoth le_rfl, hρ₀]
  have hloc_c : ∀ c, c < 4 → (t.pieceR g i s).loc (t.rQ i s ((t.nodeRot i).rot c)) =
      some ((t.rgR g i s).skt (t.rgR g i s).m c - 4 * ((t.rgR g i s).m + 1)) := by
    intro c hc
    rw [hskt c hc]
    exact t.locQ_R hwf hsh hrep hsep hi hR hloc hclosed hcor s h (hcT c hc).1 (hcT c hc).2
  have hSk := t.rgR_Sk hwf hi hR g s
  have hvert_c : ∀ c, c < 4 → QE.vert (t.pieceR g i s).es
      ((t.rgR g i s).skt (t.rgR g i s).m c - 4 * ((t.rgR g i s).m + 1)) = QE.vert (t.skR i) c := by
    intro c hc
    rw [hskt c hc, t.vert_slot_R hwf hsh hrep hsep hi hR hloc hclosed hcor s h (hcT c hc).1 (hcT c hc).2]
    have hcs : c < (t.rgR g i s).ρ₀.size := by rw [H.ρ₀_size]; omega
    have := H.planar₀.same_vertex c hcs ((t.nodeRot i).rot c)
      (by rw [RotationSystem.get_eq_rot H.planar₀.total hcs]; rfl)
    rw [hSk] at this
    exact this.symm
  have hσs : σ.size = 4 * (t.pieceR g i s).ves.length := by
    rw [hpl.size, ← t.pieceR_es hwf hi hR g s]; simp [Piece.es]
  have hb : ∀ c, c < 4 → (t.rgR g i s).skt (t.rgR g i s).m c - 4 * ((t.rgR g i s).m + 1) < σ.size := by
    intro c hc; rw [hσs]; exact Piece.loc_lt (hloc_c c hc)
  have hget : ∀ a b, a < σ.size → σ.rot a = b → σ.get a = some b := fun a b ha hab => by
    rw [RotationSystem.get_eq_rot hpl.total ha, hab]
  have hrot1 : σ.rot ((t.rgR g i s).skt (t.rgR g i s).m 1 - 4 * ((t.rgR g i s).m + 1)) =
      (t.rgR g i s).skt (t.rgR g i s).m 0 - 4 * ((t.rgR g i s).m + 1) := by
    rw [← hcap0, RotationSystem.rot_rot hpl.total hpl.involution (hb 0 (by omega))]
  have hrot3 : σ.rot ((t.rgR g i s).skt (t.rgR g i s).m 3 - 4 * ((t.rgR g i s).m + 1)) =
      (t.rgR g i s).skt (t.rgR g i s).m 2 - 4 * ((t.rgR g i s).m + 1) := by
    rw [← hcap2, RotationSystem.rot_rot hpl.total hpl.involution (hb 2 (by omega))]
  have hdir : ∀ c, c < 4 → t.rQ i s ((t.nodeRot i).rot c) % 2 = (c + 1) % 2 := by
    intro c hc
    rw [t.dir_slot_R hwf hrep hsep hi hR s h (hcT c hc).1 (hcT c hc).2]
    have := (t.T_facts hwf hrep hi hR hloc (l := c) (by omega)).2.2.2
    omega
  have hcapped : (t.pieceR g i s).Capped s'.rotAdj σ (t.rQ i s ((t.nodeRot i).rot 1))
      (t.rQ i s ((t.nodeRot i).rot 0)) (t.rQ i s ((t.nodeRot i).rot 3)) (t.rQ i s ((t.nodeRot i).rot 2))
      p.1 p.2 := by
    refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
    · rw [t.pieceR_es hwf hi hR g s]; exact hpl
    · intro q r hq hr
      rcases (t.mem_pieceR_iff (s := s)).1 hq with ⟨v, hv, hqv⟩ | ⟨j, hj, hqc⟩
      · exact t.agrees_V_R hwf hsh hrep hsep hi hR hloc hclosed hcor s h σ s'.rotAdj hσ hAw hAo hv hqv hr
      · exact t.agrees_C_R hwf hsh hrep hsep hi hR hloc hclosed hcor s h σ s'.rotAdj hσ hA hAcap hAo hj hqc hr
    · exact t.unset_R hwf hsh hrep hsep hi hR hloc hclosed hcor s h s'.rotAdj hA hAcap hAw hAo
    · refine ⟨_, _, hloc_c 1 (by omega), hloc_c 0 (by omega), hget _ _ (hb 1 (by omega)) hrot1, ?_⟩
      rw [hvert_c 1 (by omega)]
      simp [QE.vert, QE.edge, QE.side, hsk0]
    · refine ⟨_, _, hloc_c 3 (by omega), hloc_c 2 (by omega), hget _ _ (hb 3 (by omega)) hrot3, ?_⟩
      rw [hvert_c 3 (by omega)]
      simp [QE.vert, QE.edge, QE.side, hsk0]
    · exact ⟨_, _, hloc_c 1 (by omega), hloc_c 3 (by omega),
        (RotationSystem.sameFaceOrbit_iff hpl.total hpl.involution hσs (hb 1 (by omega))).2 hface⟩
    · have := hdir 1 (by omega); omega
    · have := hdir 0 (by omega); omega
    · have := hdir 3 (by omega); omega
    · have := hdir 2 (by omega); omega
  have hperm := t.pieceR_perm hwf hsh hrep hsep hi hR hloc hclosed hcor s h
  refine ⟨fun q hq => hAo q ?_ ?_, t.rQ i s ((t.nodeRot i).rot 1), t.rQ i s ((t.nodeRot i).rot 0),
    t.rQ i s ((t.nodeRot i).rot 3), t.rQ i s ((t.nodeRot i).rot 2), ?_⟩
  · intro l hl4 hl heq
    apply hq
    rw [heq]
    exact t.mem_pieceBelow_of_child hwf hsh hi (t.sC_data hwf hi hR (j := l / 4 - 1) (by omega)).1
      (t.memQ_R hwf hrep hsep hi hR s h hl4 hl).1
  · intro l hc
    obtain ⟨v, hv, hvl⟩ := (t.corner_iff (s := s)).1 hc
    subst hvl
    obtain ⟨m0, m1⟩ := t.memW_R hwf hrep hsep hi hR hloc hclosed hcor s h hv
    have hvc := (t.vItem_data hwf hrep hsep hi hR hloc hclosed hcor s h hv).1
    exact ⟨fun heq => hq (heq ▸ t.mem_pieceBelow_of_child hwf hsh hi hvc m0),
      fun heq => hq (heq ▸ t.mem_pieceBelow_of_child hwf hsh hi hvc m1)⟩
  obtain ⟨ρ, hρ⟩ := hcapped.perm hperm (hperm.nodup_iff.2 (t.edgesBelow_nodup hwf i))
  refine ⟨ρ, ?_, hρ⟩
  have e0 := hrk 0 (by omega); have e1 := hrk 1 (by omega)
  have e2 := hrk 2 (by omega); have e3 := hrk 3 (by omega)
  simp only [QE.side, QE.dir] at e0 e1 e2 e3
  norm_num at e0 e1 e2 e3
  rw [hrow, array4_ext hr4 e1 e0 e3 e2]

end Assemble

end Spqr.PlanarSpqrTree
