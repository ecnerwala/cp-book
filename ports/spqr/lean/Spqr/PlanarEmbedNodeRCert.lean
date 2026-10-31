import Spqr.PlanarEmbedNodeRExec
import Spqr.PlanarEmbedNode
import Spqr.PlanarEmbedQ
import Spqr.PlanarEmbedNodeSChain
import Spqr.Proofs.ThreeConn
import Spqr.Proofs.PlanarRGlue

/-!
# The `R` node: node-vertex and child data

Node-vertex facts of an `R` node (cap endpoints, original vertices, `V` items of the inner
node-vertices) and the `RFold` instance of the executable fold.
-/

namespace Spqr.PlanarSpqrTree

open EmbedM

variable (t : PlanarSpqrTree) {g : Graph} {i : Nat}

/-- The child slot read through `treeQe` at the node quarter-edge `4 * neSt + l` of an `R` node. -/
def rQ (i : Nat) (s : EmbedState) (l : Nat) : Nat :=
  (s.outerE[(t.sC i)[l / 4 - 1]!]![l % 4]!).getD 0

theorem capNvs_R (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hR : t.toSpqrTree.type i = .R) :
    (t.nodeEdges[(t.toSpqrTree.neRange i).1]!).nvs =
      ((t.toSpqrTree.nvRange i).1, (t.toSpqrTree.nvRange i).2 - 1) := by
  obtain ⟨_, _, -, -, -, hcap⟩ := hwf.own.nv_layout i hi
  have h6 := (t.shape_R hwf hi hR).2.1
  have := hcap (t.hasCap_R hR) _ (t.nvsOf_R hwf hi (k := 0) (by omega))
  rw [Nat.add_zero] at this
  exact Prod.ext this.1 this.2

theorem capEnd_R (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hR : t.toSpqrTree.type i = .R)
    (nv : Nat) :
    t.toSpqrTree.CapEnd i nv ↔ nv = (t.toSpqrTree.nvRange i).1 ∨ nv = (t.toSpqrTree.nvRange i).2 - 1 := by
  have h6 := (t.shape_R hwf hi hR).2.1
  have hen := t.neEn_R hwf hi
  have hsz : (t.toSpqrTree.neRange i).1 < t.nodeEdges.size := by
    have := PlanarRot.neRange_le_size t.toSpqrTree hwf i hi; omega
  have hc := t.capNvs_R hwf hi hR
  constructor
  · rintro ⟨ne, d, hne, hd, hend⟩
    rw [t.capNe_R hR] at hne
    cases hne
    rw [← getElem!_of_getElem? hd, hc] at hend
    rcases hend with h | h
    · exact Or.inl h.symm
    · exact Or.inr h.symm
  · intro h
    refine ⟨_, _, t.capNe_R hR, Array.getElem?_eq_getElem hsz, ?_⟩
    rw [← getElem!_pos t.nodeEdges _ hsz, hc]
    rcases h with h | h
    · exact Or.inl h.symm
    · exact Or.inr h.symm

/-- Every node-vertex of a node has an original vertex, `sV`. -/
theorem nvOrig_nv (hwf : t.toSpqrTree.WF) (hi : i < t.size) {j : Nat}
    (hj : j < t.toSpqrTree.nVerts i) :
    t.toSpqrTree.nvOrig ((t.toSpqrTree.nvRange i).1 + j) = some (t.sV i j) := by
  have hnv : t.toSpqrTree.NvOf i ((t.toSpqrTree.nvRange i).1 + j) := by
    constructor <;> unfold SpqrTree.nVerts at hj <;> omega
  obtain ⟨d, hd, -⟩ := t.nodeVerts_nv hwf hi hnv
  have hbd : (t.toSpqrTree.nvRange i).1 + j < t.nodeVerts.size := (Array.getElem?_eq_some_iff.1 hd).1
  have hV := hwf.own.nv_vert _ hbd
  rw [getElem!_of_getElem? hd] at hV
  have hds : d.vert < t.size := by
    unfold SpqrTree.type at hV
    cases h : t.types[d.vert]? with
    | none => rw [h] at hV; simp at hV
    | some _ => exact (Array.getElem?_eq_some_iff.1 h).1
  obtain ⟨v, hv, -⟩ := hwf.bij.vert_orig d.vert hds hV
  have ho : t.toSpqrTree.nvOrig ((t.toSpqrTree.nvRange i).1 + j) = some v := by
    unfold SpqrTree.nvOrig
    rw [hd]
    show t.origId[d.vert]?.getD none = some v
    cases h : t.origId[d.vert]? with
    | none => rw [getElem!_neg t.origId d.vert (by simpa using h)] at hv; simp at hv
    | some o => rw [getElem!_of_getElem? h] at hv; simpa using hv
  rw [ho]; unfold sV; rw [ho]; rfl

/-- The `V` item of a node-vertex of an `R` node that is not a cap endpoint. -/
theorem sW_R (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hR : t.toSpqrTree.type i = .R) {nv : Nat} (hnv : t.toSpqrTree.NvOf i nv)
    (hce : ¬ t.toSpqrTree.CapEnd i nv) :
    (∃ d, t.nodeVerts[nv]? = some d ∧ d.vert = t.sW nv) ∧
    t.sW nv ∈ t.children i ∧ t.toSpqrTree.type (t.sW nv) = .V ∧
    t.origId[t.sW nv]! = some (t.sV i (nv - (t.toSpqrTree.nvRange i).1)) := by
  obtain ⟨d, hd, -⟩ := t.nodeVerts_nv hwf hi hnv
  have hdv : d.vert = t.sW nv := by
    unfold sW; rw [getElem!_of_getElem? hd]
  have hpar := (hsep.nv_child i _ d hi (Or.inr (Or.inr hR)) hnv hd).2 hce
  rw [hdv] at hpar
  have hws : t.sW nv < t.size := by
    unfold SpqrTree.parent at hpar
    cases hp : t.par[t.sW nv]? with
    | none => rw [hp] at hpar; simp at hpar
    | some _ => exact lt_of_lt_of_eq (Array.getElem?_eq_some_iff.1 hp).1 hwf.sizes.par
  have hbd : nv < t.nodeVerts.size := (Array.getElem?_eq_some_iff.1 hd).1
  refine ⟨⟨d, hd, hdv⟩, t.mem_children_of_parent hwf hi hws hpar, hwf.own.nv_vert _ hbd, ?_⟩
  have ho := t.nvOrig_nv hwf hi (j := nv - (t.toSpqrTree.nvRange i).1)
    (by unfold SpqrTree.nVerts; have := hnv.1; have := hnv.2; omega)
  rw [show (t.toSpqrTree.nvRange i).1 + (nv - (t.toSpqrTree.nvRange i).1) = nv by have := hnv.1; omega] at ho
  unfold SpqrTree.nvOrig at ho
  rw [hd] at ho
  have ho2 : t.origId[d.vert]?.getD none = some (t.sV i (nv - (t.toSpqrTree.nvRange i).1)) := ho
  rw [hdv] at ho2
  cases ho' : t.origId[t.sW nv]? with
  | none => rw [ho'] at ho2; simp at ho2
  | some o =>
    rw [ho'] at ho2
    simp only [Option.getD_some] at ho2
    rw [getElem!_of_getElem? ho', ho2]


theorem sC_inj_R (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hR : t.toSpqrTree.type i = .R)
    {j j' : Nat} (hj : j < t.toSpqrTree.nEdges i - 1) (hj' : j' < t.toSpqrTree.nEdges i - 1)
    (h : (t.sC i)[j]! = (t.sC i)[j']!) : j = j' := by
  have hl := t.sC_length_R hwf hi hR
  have h1 := t.sC_getElem?_R hwf hi hR hj
  have h2 := t.sC_getElem?_R hwf hi hR hj'
  rw [List.getElem?_eq_getElem (by omega)] at h1 h2
  exact (List.Nodup.getElem_inj_iff (t.sC_nodup hwf hi)).1
    ((Option.some.inj h1).trans (h.trans (Option.some.inj h2).symm))

theorem rQ_eq (i : Nat) (s : EmbedState) {j r : Nat} (hr : r < 4) :
    t.rQ i s (4 * (j + 1) + r) = (s.outerE[(t.sC i)[j]!]![r]!).getD 0 := by
  have e1 : (4 * (j + 1) + r) / 4 - 1 = j := by omega
  have e2 : (4 * (j + 1) + r) % 4 = r := by omega
  unfold rQ; rw [e1, e2]

/-- The `j`-th non-`V` child of an `R` node: a capped child whose four exposed slots are the `rQ`
values of node-edge `j + 1`, with a same-witness `Capped` certificate at `sV` of its node-vertices. -/
theorem child_cert_R (hwf : t.toSpqrTree.WF) (hne : t.ne = g.ne) (hsep : t.toSpqrTree.PieceSep g)
    (hi : i < t.size) (hR : t.toSpqrTree.type i = .R) {j : Nat} (hj : j < t.toSpqrTree.nEdges i - 1)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    (t.sC i)[j]! ∈ t.children i ∧ i < (t.sC i)[j]! ∧ (t.sC i)[j]! < t.size ∧
    (∀ r, r < 4 → s.outerE[(t.sC i)[j]!]![r]! = some (t.rQ i s (4 * (j + 1) + r))) ∧
    ∃ ρ, (t.pieceBelow g (t.sC i)[j]!).Capped s.rotAdj ρ (t.rQ i s (4 * (j + 1)))
      (t.rQ i s (4 * (j + 1) + 1)) (t.rQ i s (4 * (j + 1) + 2)) (t.rQ i s (4 * (j + 1) + 3))
      (t.sV i ((t.nodeEdges[(t.toSpqrTree.neRange i).1 + (j + 1)]!).nvs.1 - (t.toSpqrTree.nvRange i).1))
      (t.sV i ((t.nodeEdges[(t.toSpqrTree.neRange i).1 + (j + 1)]!).nvs.2 - (t.toSpqrTree.nvRange i).1)) := by
  have hc := t.sC_getElem?_R hwf hi hR hj
  obtain ⟨hcm, -, hcap, hcI, hcO, -, horig, -⟩ := t.child_R (g := g) hwf hsep hi hR hj hc
  obtain ⟨hcs, -, hic⟩ := t.child_data hwf hi hcm
  have hp := t.neOrig_R (g := g) hwf hsep hi hR (k := j + 1) (by omega)
  rw [horig] at hp
  obtain ⟨c0, c1, c2, c3, ρ, h0, h1, h2, h3, hcap⟩ := t.capped_child g hwf hne hsep hi hcm hcI hcO hcap hp s h
  have q0 : t.rQ i s (4 * (j + 1)) = c0 := by
    rw [show 4 * (j + 1) = 4 * (j + 1) + 0 by rfl, t.rQ_eq i s (by omega), h0]; rfl
  have q1 : t.rQ i s (4 * (j + 1) + 1) = c1 := by rw [t.rQ_eq i s (by omega), h1]; rfl
  have q2 : t.rQ i s (4 * (j + 1) + 2) = c2 := by rw [t.rQ_eq i s (by omega), h2]; rfl
  have q3 : t.rQ i s (4 * (j + 1) + 3) = c3 := by rw [t.rQ_eq i s (by omega), h3]; rfl
  refine ⟨hcm, hic, hcs, ?_, ρ, ?_⟩
  · intro r hr
    rw [t.rQ_eq i s hr]
    interval_cases r
    · rw [h0]; rfl
    · rw [h1]; rfl
    · rw [h2]; rfl
    · rw [h3]; rfl
  · rw [q0, q1, q2, q3]; exact hcap

/-- The `V` item read at a corner quarter-edge `4 * neSt + l` is the `V` item of the second
node-vertex of node-edge `l / 4`. -/
theorem vOf_R (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hR : t.toSpqrTree.type i = .R) {l : Nat}
    (hl : l < 4 * t.toSpqrTree.nEdges i) :
    t.vOf (4 * (t.toSpqrTree.neRange i).1 + l) =
      t.sW (t.nodeEdges[(t.toSpqrTree.neRange i).1 + l / 4]!).nvs.2 ∧
    t.toSpqrTree.NvOf i (t.nodeEdges[(t.toSpqrTree.neRange i).1 + l / 4]!).nvs.2 := by
  have hnvs := t.nvs_R hwf hi hR (k := l / 4) (by omega)
  refine ⟨?_, ⟨by omega, hnvs.2.2⟩⟩
  unfold vOf sW QE.edge
  rw [show (4 * (t.toSpqrTree.neRange i).1 + l) / 4 = (t.toSpqrTree.neRange i).1 + l / 4 by omega]


end Spqr.PlanarSpqrTree

/-! ### Slots of a `Capped` certificate -/

namespace Spqr

section CappedSlots

variable {P : Piece} {A : Array (Option Nat)} {ρ : RotationSystem} {c0 c1 c2 c3 u v : Nat}

theorem Piece.Capped.mem_slots (hc : P.Capped A ρ c0 c1 c2 c3 u v) :
    P.Mem c0 ∧ P.Mem c1 ∧ P.Mem c2 ∧ P.Mem c3 := by
  obtain ⟨l0, l1, h0, h1, -, -⟩ := hc.pair0
  obtain ⟨l2, l3, h2, h3, -, -⟩ := hc.pair2
  exact ⟨Piece.mem_of_loc h0, Piece.mem_of_loc h1, Piece.mem_of_loc h2, Piece.mem_of_loc h3⟩

theorem Piece.Capped.slots_lt (hc : P.Capped A ρ c0 c1 c2 c3 u v) :
    c0 < A.size ∧ c1 < A.size ∧ c2 < A.size ∧ c3 < A.size := by
  obtain ⟨m0, m1, m2, m3⟩ := hc.mem_slots
  refine ⟨?_, ?_, ?_, ?_⟩
  · exact (Array.getElem?_eq_some_iff.1 ((hc.unset c0 m0).2 (Or.inl rfl))).1
  · exact (Array.getElem?_eq_some_iff.1 ((hc.unset c1 m1).2 (Or.inr (Or.inl rfl)))).1
  · exact (Array.getElem?_eq_some_iff.1 ((hc.unset c2 m2).2 (Or.inr (Or.inr (Or.inl rfl))))).1
  · exact (Array.getElem?_eq_some_iff.1 ((hc.unset c3 m3).2 (Or.inr (Or.inr (Or.inr rfl))))).1

theorem Piece.Capped.slots_unset (hc : P.Capped A ρ c0 c1 c2 c3 u v) :
    A[c0]? = some none ∧ A[c1]? = some none ∧ A[c2]? = some none ∧ A[c3]? = some none := by
  obtain ⟨m0, m1, m2, m3⟩ := hc.mem_slots
  exact ⟨(hc.unset c0 m0).2 (Or.inl rfl), (hc.unset c1 m1).2 (Or.inr (Or.inl rfl)),
    (hc.unset c2 m2).2 (Or.inr (Or.inr (Or.inl rfl))), (hc.unset c3 m3).2 (Or.inr (Or.inr (Or.inr rfl)))⟩

theorem Piece.Capped.slots_distinct (hc : P.Capped A ρ c0 c1 c2 c3 u v) (huv : u ≠ v) :
    c0 ≠ c1 ∧ c0 ≠ c2 ∧ c0 ≠ c3 ∧ c1 ≠ c2 ∧ c1 ≠ c3 ∧ c2 ≠ c3 := by
  have d0 := hc.dir0; have d1 := hc.dir1; have d2 := hc.dir2; have d3 := hc.dir3
  obtain ⟨l0, l1, h0, h1, g0, v0⟩ := hc.pair0
  obtain ⟨l2, l3, h2, h3, g2, v2⟩ := hc.pair2
  have h02 : c0 ≠ c2 := by
    rintro rfl
    rw [h0] at h2; cases h2
    rw [v0] at v2; exact huv (Option.some.inj v2)
  refine ⟨by omega, h02, by omega, by omega, ?_, by omega⟩
  rintro rfl
  rw [h1] at h3; cases h3
  obtain ⟨a0, ha0, -⟩ := Option.bind_eq_some_iff.1 g0
  obtain ⟨a2, ha2, -⟩ := Option.bind_eq_some_iff.1 g2
  have hl0 : l0 < ρ.size := (Array.getElem?_eq_some_iff.1 ha0).1
  have hl2 : l2 < ρ.size := (Array.getElem?_eq_some_iff.1 ha2).1
  have i0 := (hc.planar.involution l0 hl0 l1 g0).2
  have i2 := (hc.planar.involution l2 hl2 l1 g2).2
  rw [i0] at i2
  exact h02 (Piece.loc_injective h0 (i2 ▸ h2))

/-- The vertex of every slot of a `Capped` certificate. -/
theorem Piece.Capped.slots_vert (hc : P.Capped A ρ c0 c1 c2 c3 u v) :
    (∃ l, P.loc c0 = some l ∧ QE.vert P.es l = some u) ∧
    (∃ l, P.loc c1 = some l ∧ QE.vert P.es l = some u) ∧
    (∃ l, P.loc c2 = some l ∧ QE.vert P.es l = some v) ∧
    (∃ l, P.loc c3 = some l ∧ QE.vert P.es l = some v) := by
  obtain ⟨l0, l1, h0, h1, g0, v0⟩ := hc.pair0
  obtain ⟨l2, l3, h2, h3, g2, v2⟩ := hc.pair2
  obtain ⟨a0, ha0, -⟩ := Option.bind_eq_some_iff.1 g0
  obtain ⟨a2, ha2, -⟩ := Option.bind_eq_some_iff.1 g2
  have hl0 : l0 < ρ.size := (Array.getElem?_eq_some_iff.1 ha0).1
  have hl2 : l2 < ρ.size := (Array.getElem?_eq_some_iff.1 ha2).1
  exact ⟨⟨l0, h0, v0⟩, ⟨l1, h1, (hc.planar.same_vertex l0 hl0 l1 g0).symm.trans v0⟩,
    ⟨l2, h2, v2⟩, ⟨l3, h3, (hc.planar.same_vertex l2 hl2 l3 g2).symm.trans v2⟩⟩

end CappedSlots

end Spqr

namespace Spqr.PlanarSpqrTree

open EmbedM

variable (t : PlanarSpqrTree) {g : Graph} {i : Nat}

/-! ### Corners of an `R` node -/

/-- A corner of the executable fold is a `CornerAt` of the second node-vertex of its edge, which
is therefore not a cap endpoint. -/
theorem corner_R (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hR : t.toSpqrTree.type i = .R)
    (hloc : IsPlanarEmbedding (t.localSkeleton i) (t.toSpqrTree.nVerts i) (t.nodeRot i))
    (hclosed : t.NodeRotClosed) (s : EmbedState) {l : Nat}
    (hc : t.Corner (t.toSpqrTree.neRange i).1 (t.toSpqrTree.nEdges i) s (t.nodeRot i).rot l) :
    t.CornerAt i (t.nodeEdges[(t.toSpqrTree.neRange i).1 + l / 4]!).nvs.2
      (4 * (t.toSpqrTree.neRange i).1 + l) := by
  obtain ⟨hl4, hl, hl2, hT, -⟩ := hc
  have hsz := t.nodeRot_size_le hwf hi hR hloc
  have hen := t.neEn_R hwf hi
  have hne := PlanarRot.neRange_le_size t.toSpqrTree hwf i hi
  have hrot := (t.rotR hwf hi hR hloc hclosed hl).1
  refine ⟨by omega, by omega, by omega, ⟨_, ?_, rfl⟩, 4 * (t.toSpqrTree.neRange i).1 + (t.nodeRot i).rot l, ?_, by omega⟩
  · rw [show (4 * (t.toSpqrTree.neRange i).1 + l) / 4 = (t.toSpqrTree.neRange i).1 + l / 4 by omega]
    rw [getElem!_pos t.nodeEdges ((t.toSpqrTree.neRange i).1 + l / 4) (by omega)]
    exact Array.getElem?_eq_getElem (by omega)
  · rw [getElem!_pos t.neRotAdj _ (by omega)] at hrot
    rw [Array.getElem?_eq_getElem (by omega), hrot]

theorem corner_nv_R (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hR : t.toSpqrTree.type i = .R)
    (hloc : IsPlanarEmbedding (t.localSkeleton i) (t.toSpqrTree.nVerts i) (t.nodeRot i))
    (hclosed : t.NodeRotClosed) (hcor : t.NodeCorners) (s : EmbedState) {l : Nat}
    (hc : t.Corner (t.toSpqrTree.neRange i).1 (t.toSpqrTree.nEdges i) s (t.nodeRot i).rot l) :
    t.toSpqrTree.NvOf i (t.nodeEdges[(t.toSpqrTree.neRange i).1 + l / 4]!).nvs.2 ∧
    ¬ t.toSpqrTree.CapEnd i (t.nodeEdges[(t.toSpqrTree.neRange i).1 + l / 4]!).nvs.2 ∧
    l < (t.nodeRot i).rot l := by
  have hca := t.corner_R hwf hi hR hloc hclosed s hc
  have hnv := (t.vOf_R hwf hi hR hc.2.1).2
  have hnce : ¬ t.toSpqrTree.CapEnd i _ := fun h => hcor.capEnd i _ _ hi (Or.inr (Or.inr hR)) h hca
  obtain ⟨ta, -, hlt, huniq⟩ := hcor.inner i _ hi (Or.inr (Or.inr hR)) hnv hnce
  have := huniq _ hca
  subst this
  obtain ⟨-, -, -, -, tb, htb, -⟩ := hca
  have h1 := hlt tb htb
  have hsz := t.nodeRot_size_le hwf hi hR hloc
  have hl := hc.2.1
  have hrot := (t.rotR hwf hi hR hloc hclosed hl).1
  rw [getElem!_pos t.neRotAdj _ (by have := t.neEn_R hwf hi; omega)] at hrot
  rw [Array.getElem?_eq_getElem (by have := t.neEn_R hwf hi; omega), hrot] at htb
  cases Option.some.inj (Option.some.inj htb)
  exact ⟨hnv, hnce, by omega⟩

/-- The `V` item at a corner and its open embedding. -/
theorem corner_V_R (hwf : t.toSpqrTree.WF) (hne : t.ne = g.ne) (hsep : t.toSpqrTree.PieceSep g)
    (hi : i < t.size) (hR : t.toSpqrTree.type i = .R)
    (hloc : IsPlanarEmbedding (t.localSkeleton i) (t.toSpqrTree.nVerts i) (t.nodeRot i))
    (hclosed : t.NodeRotClosed) (hcor : t.NodeCorners) (s : EmbedState)
    (h : t.GluedFaces g (i + 1) s) {l : Nat}
    (hc : t.Corner (t.toSpqrTree.neRange i).1 (t.toSpqrTree.nEdges i) s (t.nodeRot i).rot l) :
    t.vOf (4 * (t.toSpqrTree.neRange i).1 + l) ∈ t.children i ∧
    t.toSpqrTree.type (t.vOf (4 * (t.toSpqrTree.neRange i).1 + l)) = .V ∧
    ∃ w0 w1 ρ, s.outerE[t.vOf (4 * (t.toSpqrTree.neRange i).1 + l)]![0]! = some w0 ∧
      s.outerE[t.vOf (4 * (t.toSpqrTree.neRange i).1 + l)]![1]! = some w1 ∧
      (t.pieceBelow g (t.vOf (4 * (t.toSpqrTree.neRange i).1 + l))).OpenEmbedding s.rotAdj ρ w0 w1
        (t.sV i ((t.nodeEdges[(t.toSpqrTree.neRange i).1 + l / 4]!).nvs.2 - (t.toSpqrTree.nvRange i).1)) := by
  obtain ⟨hnv, hnce, -⟩ := t.corner_nv_R hwf hi hR hloc hclosed hcor s hc
  have hv := (t.vOf_R hwf hi hR hc.2.1).1
  obtain ⟨-, hw, htw, ho⟩ := t.sW_R (g := g) hwf hsep hi hR hnv hnce
  rw [hv]
  refine ⟨hw, htw, ?_⟩
  rcases t.q_lower_boundary g hwf hne hi hw htw ho s h.toGluedUpTo with ⟨h0, -⟩ | ⟨w0, w1, ρ, h0, h1, hop⟩
  · have := hc.2.2.2.2
    rw [hv, h0] at this
    simp at this
  · exact ⟨w0, w1, ρ, h0, h1, hop⟩


/-! ### The `RFold` instance -/

theorem sV_ne_R (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hR : t.toSpqrTree.type i = .R) {k : Nat} (hk : k < t.toSpqrTree.nEdges i) :
    t.sV i ((t.nodeEdges[(t.toSpqrTree.neRange i).1 + k]!).nvs.1 - (t.toSpqrTree.nvRange i).1) ≠
      t.sV i ((t.nodeEdges[(t.toSpqrTree.neRange i).1 + k]!).nvs.2 - (t.toSpqrTree.nvRange i).1) := by
  intro heq
  have hnvs := t.nvs_R hwf hi hR hk
  have h1 := t.nvOrig_nv hwf hi (j := (t.nodeEdges[(t.toSpqrTree.neRange i).1 + k]!).nvs.1 - (t.toSpqrTree.nvRange i).1)
    (by unfold SpqrTree.nVerts; omega)
  have h2 := t.nvOrig_nv hwf hi (j := (t.nodeEdges[(t.toSpqrTree.neRange i).1 + k]!).nvs.2 - (t.toSpqrTree.nvRange i).1)
    (by unfold SpqrTree.nVerts; omega)
  rw [Nat.add_sub_of_le hnvs.1] at h1
  rw [Nat.add_sub_of_le (by omega)] at h2
  have := hsep.nv_orig_inj i _ _ hi (Or.inr (Or.inr hR)) ⟨hnvs.1, by omega⟩ ⟨by omega, hnvs.2.2⟩
    (by rw [h1, h2, heq])
  omega

theorem slot_R (hwf : t.toSpqrTree.WF) (hne : t.ne = g.ne) (hsep : t.toSpqrTree.PieceSep g)
    (hi : i < t.size) (hR : t.toSpqrTree.type i = .R) (s : EmbedState) (h : t.GluedFaces g (i + 1) s)
    {j r : Nat} (hj : j < t.toSpqrTree.nEdges i - 1) (hr : r < 4) :
    (t.pieceBelow g (t.sC i)[j]!).Mem (t.rQ i s (4 * (j + 1) + r)) ∧
    s.rotAdj[t.rQ i s (4 * (j + 1) + r)]? = some none := by
  obtain ⟨-, -, -, -, ρ, hc⟩ := t.child_cert_R hwf hne hsep hi hR hj s h
  obtain ⟨m0, m1, m2, m3⟩ := hc.mem_slots
  obtain ⟨u0, u1, u2, u3⟩ := hc.slots_unset
  interval_cases r
  · rw [Nat.add_zero]; exact ⟨m0, u0⟩
  · exact ⟨m1, u1⟩
  · exact ⟨m2, u2⟩
  · exact ⟨m3, u3⟩

/-- The executable fold data of an `R` node: `T` is the local rotation `nodeRot`, `Q` the child
slots `rQ`. -/
theorem rfold_R (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hrep : t.toSpqrTree.Represents g) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hR : t.toSpqrTree.type i = .R)
    (hloc : IsPlanarEmbedding (t.localSkeleton i) (t.toSpqrTree.nVerts i) (t.nodeRot i))
    (hclosed : t.NodeRotClosed) (hcor : t.NodeCorners) (s : EmbedState)
    (h : t.GluedFaces g (i + 1) s) :
    t.RFold i (t.toSpqrTree.neRange i).1 (t.toSpqrTree.nEdges i) s (t.nodeRot i).rot (t.rQ i s) := by
  have hne := hrep.ne
  have hsz := t.nodeRot_size_le hwf hi hR hloc
  have hs : (t.nodeRot i).size = 4 * t.toSpqrTree.nEdges i := by
    rw [hloc.size, t.localSkeleton_length]
  have hen := t.neEn_R hwf hi
  have h6 := (t.shape_R hwf hi hR).2.1
  have dec : ∀ l, 4 ≤ l → l < 4 * t.toSpqrTree.nEdges i →
      ∃ j r, j < t.toSpqrTree.nEdges i - 1 ∧ r < 4 ∧ l = 4 * (j + 1) + r :=
    fun l hl4 hl => ⟨l / 4 - 1, l % 4, by omega, by omega, by omega⟩
  have hdisj : ∀ a b, a ∈ t.children i → b ∈ t.children i → a ≠ b → ∀ q,
      (t.pieceBelow g a).Mem q → (t.pieceBelow g b).Mem q → False :=
    fun a b ha hb hab q hqa hqb => t.maximal_pieces_disjoint hwf hsh
      (t.child_maximal hwf hi ha) (t.child_maximal hwf hi hb) hab hqa hqb
  refine
    { rot := fun l hl => (t.rotR hwf hi hR hloc hclosed hl).1
      T_lt := fun l hl => (t.rotR hwf hi hR hloc hclosed hl).2
      T_inv := fun l hl => RotationSystem.rot_rot hloc.total hloc.involution (by rw [hs]; exact hl)
      T_edge := ?_, qe := ?_, Q_lt := ?_, Q_inj := ?_, v_ne := ?_
      v_lt := fun l hc => (t.corner_nv_R hwf hi hR hloc hclosed hcor s hc).2.2
      v_some1 := ?_, w_ne := ?_, w_lt := ?_, wQ := ?_, ww := ?_, row := ?_
      nE_pos := by omega }
  · intro l hl
    have he := t.localSkeleton_getElem? hwf hi (k := l / 4) (by omega)
    have hnvs := t.nvs_R hwf hi hR (k := l / 4) (by omega)
    have huv : (t.nodeEdges[(t.toSpqrTree.neRange i).1 + l / 4]!).nvs.1 - (t.toSpqrTree.nvRange i).1 ≠
        (t.nodeEdges[(t.toSpqrTree.neRange i).1 + l / 4]!).nvs.2 - (t.toSpqrTree.nvRange i).1 := by omega
    have h3 : SpqrTree.ThreeConnected (t.toSpqrTree.nVerts i) (t.localSkeleton i) :=
      hrep.r_three_connected i hi hR
    have hconn := SpqrTree.ThreeConnected.edgesConn_eraseIdx h3 hloc.verts he huv
    have hget : (t.nodeRot i).get (4 * (l / 4) + l % 4) = some ((t.nodeRot i).rot l) := by
      rw [show 4 * (l / 4) + l % 4 = l by omega]
      exact RotationSystem.get_eq_rot hloc.total (by rw [hs]; exact hl)
    exact deg_of_edgesConn hloc he huv hconn (l % 4) (by omega) _ hget
  · intro l hl4 hl s' hs'
    obtain ⟨j, r, hj, hr, rfl⟩ := dec l hl4 hl
    obtain ⟨hcm, hic, -, hslot, -⟩ := t.child_cert_R hwf hne hsep hi hR hj s h
    have hc := t.sC_getElem?_R hwf hi hR hj
    have htq := (t.child_R (g := g) hwf hsep hi hR hj hc).2.2.2.2.2.2.2 r hr s'
    rw [show 4 * (t.toSpqrTree.neRange i).1 + (4 * (j + 1) + r) =
      4 * ((t.toSpqrTree.neRange i).1 + (j + 1)) + r by omega, htq,
      getElem!_eq_of_getElem?_eq (hs' _ (Nat.ne_of_gt hic)), hslot r hr]
  · intro l hl4 hl
    obtain ⟨j, r, hj, hr, rfl⟩ := dec l hl4 hl
    exact (Array.getElem?_eq_some_iff.1 (t.slot_R hwf hne hsep hi hR s h hj hr).2).1
  · intro l l' hl4 hl hl4' hl' heq
    obtain ⟨j, r, hj, hr, rfl⟩ := dec l hl4 hl
    obtain ⟨j', r', hj', hr', rfl⟩ := dec l' hl4' hl'
    by_cases hjj : j = j'
    · subst hjj
      obtain ⟨-, -, -, -, ρ, hc⟩ := t.child_cert_R hwf hne hsep hi hR hj s h
      have hd := hc.slots_distinct (t.sV_ne_R hwf hsep hi hR (k := j + 1) (by omega))
      interval_cases r <;> interval_cases r' <;> (try simp only [Nat.add_zero] at heq hd) <;> omega
    · exfalso
      have hcc : (t.sC i)[j]! ≠ (t.sC i)[j']! := fun e => hjj (t.sC_inj_R hwf hi hR hj hj' e)
      have m1 := (t.slot_R hwf hne hsep hi hR s h hj hr).1
      have m2 := (t.slot_R hwf hne hsep hi hR s h hj' hr').1
      rw [heq] at m1
      exact hdisj _ _ (t.child_cert_R hwf hne hsep hi hR hj s h).1
        (t.child_cert_R hwf hne hsep hi hR hj' s h).1 hcc _ m1 m2
  · intro l hl4 hl heq
    obtain ⟨hv, hnv⟩ := t.vOf_R hwf hi hR hl
    obtain ⟨d, hd, -⟩ := t.nodeVerts_nv hwf hi hnv
    have hV := hwf.own.nv_vert _ (Array.getElem?_eq_some_iff.1 hd).1
    rw [getElem!_of_getElem? hd] at hV
    rw [hv] at heq
    unfold sW at heq
    rw [getElem!_of_getElem? hd] at heq
    rw [heq, hR] at hV
    cases hV
  · intro l hc
    obtain ⟨-, -, a0, a1, ρ, h0, h1, -⟩ := t.corner_V_R hwf hne hsep hi hR hloc hclosed hcor s h hc
    rw [h1]; rfl
  · intro l hc
    obtain ⟨-, -, a0, a1, ρ, h0, h1, hop⟩ := t.corner_V_R hwf hne hsep hi hR hloc hclosed hcor s h hc
    unfold w0 w1
    rw [h0, h1]
    simp only [Option.getD_some]
    have := hop.left_dir; have := hop.right_dir; omega
  · intro l hc
    obtain ⟨-, -, a0, a1, ρ, h0, h1, hop⟩ := t.corner_V_R hwf hne hsep hi hR hloc hclosed hcor s h hc
    obtain ⟨la, lb, hla, hlb, -, -⟩ := hop.boundary
    have u0 := (hop.unset a0 (Piece.mem_of_loc hla)).2 (Or.inl rfl)
    have u1 := (hop.unset a1 (Piece.mem_of_loc hlb)).2 (Or.inr rfl)
    unfold w0 w1
    rw [h0, h1]
    exact ⟨(Array.getElem?_eq_some_iff.1 u0).1, (Array.getElem?_eq_some_iff.1 u1).1⟩
  · intro l hc l' hl4' hl'
    obtain ⟨hw, htw, a0, a1, ρ, h0, h1, hop⟩ := t.corner_V_R hwf hne hsep hi hR hloc hclosed hcor s h hc
    obtain ⟨la, lb, hla, hlb, -, -⟩ := hop.boundary
    obtain ⟨j, r, hj, hr, rfl⟩ := dec l' hl4' hl'
    have hcm := (t.child_cert_R hwf hne hsep hi hR hj s h).1
    have hcV := (t.child_R (g := g) hwf hsep hi hR hj (t.sC_getElem?_R hwf hi hR hj)).2.1
    have hne' : t.vOf (4 * (t.toSpqrTree.neRange i).1 + l) ≠ (t.sC i)[j]! := fun e => hcV (e ▸ htw)
    have hm := (t.slot_R hwf hne hsep hi hR s h hj hr).1
    unfold w0 w1
    rw [h0, h1]
    simp only [Option.getD_some]
    exact ⟨fun e => hdisj _ _ hw hcm hne' _ (by rw [e]; exact Piece.mem_of_loc hla) hm,
      fun e => hdisj _ _ hw hcm hne' _ (by rw [e]; exact Piece.mem_of_loc hlb) hm⟩
  · intro l l' hc hc' hll
    obtain ⟨hw, htw, a0, a1, ρ, h0, h1, hop⟩ := t.corner_V_R hwf hne hsep hi hR hloc hclosed hcor s h hc
    obtain ⟨hw', htw', b0, b1, ρ', h0', h1', hop'⟩ :=
      t.corner_V_R hwf hne hsep hi hR hloc hclosed hcor s h hc'
    obtain ⟨hnv, hnce, -⟩ := t.corner_nv_R hwf hi hR hloc hclosed hcor s hc
    obtain ⟨hnv', -, -⟩ := t.corner_nv_R hwf hi hR hloc hclosed hcor s hc'
    have hca := t.corner_R hwf hi hR hloc hclosed s hc
    have hca' := t.corner_R hwf hi hR hloc hclosed s hc'
    have hnn : (t.nodeEdges[(t.toSpqrTree.neRange i).1 + l / 4]!).nvs.2 ≠
        (t.nodeEdges[(t.toSpqrTree.neRange i).1 + l' / 4]!).nvs.2 := by
      intro e
      obtain ⟨ta, -, -, huniq⟩ := hcor.inner i _ hi (Or.inr (Or.inr hR)) hnv hnce
      have e1 := huniq _ hca
      rw [← e] at hca'
      have e2 := huniq _ hca'
      omega
    have hvv : t.vOf (4 * (t.toSpqrTree.neRange i).1 + l) ≠
        t.vOf (4 * (t.toSpqrTree.neRange i).1 + l') := by
      intro e
      rw [(t.vOf_R hwf hi hR hc.2.1).1, (t.vOf_R hwf hi hR hc'.2.1).1] at e
      obtain ⟨d, hd, -⟩ := t.nodeVerts_nv hwf hi hnv
      obtain ⟨d', hd', -⟩ := t.nodeVerts_nv hwf hi hnv'
      unfold sW at e
      rw [getElem!_of_getElem? hd, getElem!_of_getElem? hd'] at e
      exact hnn (hsep.nv_vert_inj i _ _ d d' hi (Or.inr (Or.inr hR)) hnv hnv' hd hd' e)
    obtain ⟨la, lb, hla, hlb, -, -⟩ := hop.boundary
    obtain ⟨la', lb', hla', hlb', -, -⟩ := hop'.boundary
    unfold w0 w1
    rw [h0, h1, h0', h1']
    simp only [Option.getD_some]
    exact ⟨fun e => hdisj _ _ hw hw' hvv _ (Piece.mem_of_loc hla) (by rw [e]; exact Piece.mem_of_loc hla'),
      fun e => hdisj _ _ hw hw' hvv _ (Piece.mem_of_loc hla) (by rw [e]; exact Piece.mem_of_loc hlb'),
      fun e => hdisj _ _ hw hw' hvv _ (Piece.mem_of_loc hlb) (by rw [e]; exact Piece.mem_of_loc hla'),
      fun e => hdisj _ _ hw hw' hvv _ (Piece.mem_of_loc hlb) (by rw [e]; exact Piece.mem_of_loc hlb')⟩
  · have hsz' : i < s.outerE.size := by rw [h.outer_size]; exact hi
    refine ⟨s.outerE[i], Array.getElem?_eq_getElem hsz', ?_⟩
    have := h.outer_row_size i hi
    rwa [getElem!_pos s.outerE i hsz'] at this

end Spqr.PlanarSpqrTree
