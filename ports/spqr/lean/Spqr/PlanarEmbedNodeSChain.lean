import Spqr.PlanarEmbedNodeS
import Spqr.Proofs.PieceJoin
import Spqr.Proofs.PiecePerm

/-!
# The `S` node: the chain of glued children

`sL i m` lists the children glued after path edges `1..m` of the `S` node `i` (the first non-`V`
child, then alternately the `V` item of the next node-vertex and the next non-`V` child);
`chain g i m` is the corresponding piece. `sV i j` is the original vertex of node-vertex `s + j`.
-/

namespace Spqr.PlanarSpqrTree

variable (t : PlanarSpqrTree)

/-- The non-`V` children of `i`, in order. -/
def sC (i : Nat) : List Nat := (t.children i).filter fun c => t.toSpqrTree.type c ≠ .V

/-- The `V` item of node-vertex `nv`. -/
def sW (nv : Nat) : Nat := (t.nodeVerts[nv]!).vert

/-- The children glued after path edges `1..m`. -/
def sL (i m : Nat) : List Nat :=
  (t.sC i)[0]! :: (List.range m).flatMap fun j =>
    [t.sW ((t.toSpqrTree.nvRange i).1 + (j + 1)), (t.sC i)[j + 1]!]

/-- The piece glued after path edges `1..m`. -/
def chain (g : Graph) (i m : Nat) : Piece :=
  ⟨(t.sL i m).flatMap t.edgesBelow, fun e => g.edges[e]!, g.nv, 0, 0⟩

/-- Original vertex of node-vertex `s + j` of node `i`. -/
def sV (i j : Nat) : Nat := (t.toSpqrTree.nvOrig ((t.toSpqrTree.nvRange i).1 + j)).getD 0

theorem sL_succ (i m : Nat) :
    t.sL i (m + 1) = t.sL i m ++ [t.sW ((t.toSpqrTree.nvRange i).1 + (m + 1)), (t.sC i)[m + 1]!] := by
  simp [sL, List.range_succ, List.flatMap_append]

theorem chain_succ_ves (g : Graph) (i m : Nat) :
    (t.chain g i (m + 1)).ves = (t.chain g i m).ves ++
      (t.edgesBelow (t.sW ((t.toSpqrTree.nvRange i).1 + (m + 1))) ++ t.edgesBelow ((t.sC i)[m + 1]!)) := by
  simp [chain, sL_succ, List.flatMap_append]

theorem chain_zero_ves (g : Graph) (i : Nat) : (t.chain g i 0).ves = t.edgesBelow ((t.sC i)[0]!) := by
  simp [chain, sL]

theorem chain_eq_pieceBelow_with (g : Graph) (i m c : Nat) :
    t.chain g i m = {t.pieceBelow g c with ves := (t.chain g i m).ves} := rfl

theorem neOrig_iff (T : SpqrTree) (ne : Nat) (p : Nat × Nat) :
    T.neOrig ne = some p ↔
      ∃ nvs, T.nvsOf ne = some nvs ∧ T.nvOrig nvs.1 = some p.1 ∧ T.nvOrig nvs.2 = some p.2 := by
  show (T.nvsOf ne).bind (fun nvs => (T.nvOrig nvs.1).bind fun a =>
    (T.nvOrig nvs.2).bind fun b => some (a, b)) = some p ↔ _
  simp only [Option.bind_eq_some_iff, Option.some.injEq]
  constructor
  · rintro ⟨nvs, h1, a, h2, b, h3, rfl⟩
    exact ⟨nvs, h1, h2, h3⟩
  · rintro ⟨nvs, h1, h2, h3⟩
    exact ⟨nvs, h1, p.1, h2, p.2, h3, rfl⟩

section S

variable {i : Nat}

theorem nvsOf_S (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hS : t.toSpqrTree.type i = .S)
    {k : Nat} (hk : k < t.toSpqrTree.nVerts i) :
    t.toSpqrTree.nvsOf ((t.toSpqrTree.neRange i).1 + k) =
      some (if k = 0 then ((t.toSpqrTree.nvRange i).1, (t.toSpqrTree.nvRange i).1 + t.toSpqrTree.nVerts i - 1)
      else ((t.toSpqrTree.nvRange i).1 + k - 1, (t.toSpqrTree.nvRange i).1 + k)) := by
  have hn := t.nEdges_S hwf hi hS
  have hsz := PlanarRot.neRange_le_size t.toSpqrTree hwf i hi
  rw [← t.nvs_S hwf hi hS hk]
  unfold SpqrTree.nvsOf
  rw [Array.getElem?_eq_getElem (show (t.toSpqrTree.neRange i).1 + k < t.nodeEdges.size by omega),
    getElem!_pos t.nodeEdges ((t.toSpqrTree.neRange i).1 + k) (by omega)]
  rfl

theorem sC_length (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hS : t.toSpqrTree.type i = .S) :
    (t.sC i).length = t.toSpqrTree.nVerts i - 1 := (t.noncap_S hwf hi hS).1

theorem sC_getElem? (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hS : t.toSpqrTree.type i = .S)
    {j : Nat} (hj : j < t.toSpqrTree.nVerts i - 1) : (t.sC i)[j]? = some (t.sC i)[j]! := by
  have := t.sC_length hwf hi hS
  rw [getElem!_pos (t.sC i) j (by omega), List.getElem?_eq_getElem (by omega)]

/-- The original endpoints of path edge `k` of an `S` node are `sV (k - 1)`, `sV k`. -/
theorem neOrig_S (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hS : t.toSpqrTree.type i = .S) {k : Nat} (hk1 : 1 ≤ k) (hk : k < t.toSpqrTree.nVerts i) :
    t.toSpqrTree.neOrig ((t.toSpqrTree.neRange i).1 + k) = some (t.sV i (k - 1), t.sV i k) := by
  obtain ⟨hc, -, hcap, -, -, -, horig, -⟩ :=
    t.child_S (g := g) hwf hsep hi hS (j := k - 1) (c := (t.sC i)[k - 1]!) (by omega)
      (t.sC_getElem? hwf hi hS (by omega))
  have hcne : t.toSpqrTree.capNe (t.sC i)[k - 1]! = some (t.toSpqrTree.neRange (t.sC i)[k - 1]!).1 := by
    simp only [SpqrTree.capNe, hcap, ↓reduceIte]
  obtain ⟨p, hp⟩ := hsep.cap_orig _ _ (t.child_data hwf hi hc).1 hcne
  rw [show k - 1 + 1 = k by omega] at horig
  have hp' : t.toSpqrTree.neOrig ((t.toSpqrTree.neRange i).1 + k) = some p := horig.trans hp
  rw [hp']
  rw [neOrig_iff] at hp'
  obtain ⟨nvs, hnvs, h1, h2⟩ := hp'
  rw [t.nvsOf_S hwf hi hS hk, if_neg (by omega)] at hnvs
  cases hnvs
  simp only at h1 h2
  have e1 : t.sV i (k - 1) = p.1 := by
    unfold sV
    rw [show (t.toSpqrTree.nvRange i).1 + (k - 1) = (t.toSpqrTree.nvRange i).1 + k - 1 by omega, h1]; rfl
  have e2 : t.sV i k = p.2 := by unfold sV; rw [h2]; rfl
  rw [e1, e2]

theorem nvOrig_S (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hS : t.toSpqrTree.type i = .S) {j : Nat} (hj : j < t.toSpqrTree.nVerts i) :
    t.toSpqrTree.nvOrig ((t.toSpqrTree.nvRange i).1 + j) = some (t.sV i j) := by
  have h3 := (t.shape_S hwf hi hS).1
  by_cases hj' : j + 1 < t.toSpqrTree.nVerts i
  · have h := t.neOrig_S (g := g) hwf hsep hi hS (k := j + 1) (by omega) hj'
    rw [neOrig_iff] at h
    obtain ⟨nvs, hnvs, h1, -⟩ := h
    rw [t.nvsOf_S hwf hi hS hj', if_neg (by omega)] at hnvs
    cases hnvs
    simpa using h1
  · have h := t.neOrig_S (g := g) hwf hsep hi hS (k := j) (by omega) hj
    rw [neOrig_iff] at h
    obtain ⟨nvs, hnvs, -, h2⟩ := h
    rw [t.nvsOf_S hwf hi hS hj, if_neg (by omega)] at hnvs
    cases hnvs
    simpa using h2

/-- The cap of an `S` node has original endpoints `sV 0`, `sV (n - 1)`. -/
theorem capOrig_S (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hS : t.toSpqrTree.type i = .S) {p : Nat × Nat}
    (hp : t.toSpqrTree.neOrig (t.toSpqrTree.neRange i).1 = some p) :
    p = (t.sV i 0, t.sV i (t.toSpqrTree.nVerts i - 1)) := by
  have h3 := (t.shape_S hwf hi hS).1
  rw [neOrig_iff] at hp
  obtain ⟨nvs, hnvs, h1, h2⟩ := hp
  have := t.nvsOf_S hwf hi hS (k := 0) (by omega)
  rw [Nat.add_zero] at this
  rw [this, if_pos rfl] at hnvs
  cases hnvs
  simp only at h1 h2
  have h0 := t.nvOrig_S (g := g) hwf hsep hi hS (j := 0) (by omega)
  rw [Nat.add_zero] at h0
  rw [h0] at h1
  rw [show (t.toSpqrTree.nvRange i).1 + t.toSpqrTree.nVerts i - 1 =
    (t.toSpqrTree.nvRange i).1 + (t.toSpqrTree.nVerts i - 1) by omega,
    t.nvOrig_S (g := g) hwf hsep hi hS (by omega)] at h2
  rw [Option.some.injEq] at h1 h2
  exact Prod.ext h1.symm h2.symm

/-- The `j`-th non-`V` child of an `S` node carries a `Capped` certificate at `sV j`, `sV (j + 1)`. -/
theorem child_cert_S (hwf : t.toSpqrTree.WF) (hne : t.ne = g.ne) (hsep : t.toSpqrTree.PieceSep g)
    (hi : i < t.size) (hS : t.toSpqrTree.type i = .S) {j : Nat} (hj : j < t.toSpqrTree.nVerts i - 1)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    (t.sC i)[j]! ∈ t.children i ∧ t.toSpqrTree.type (t.sC i)[j]! ≠ .V ∧
    ∃ a0 a1 a2 a3 ρ, s.outerE[(t.sC i)[j]!]![0]! = some a0 ∧ s.outerE[(t.sC i)[j]!]![1]! = some a1 ∧
      s.outerE[(t.sC i)[j]!]![2]! = some a2 ∧ s.outerE[(t.sC i)[j]!]![3]! = some a3 ∧
      (t.pieceBelow g (t.sC i)[j]!).Capped s.rotAdj ρ a0 a1 a2 a3 (t.sV i j) (t.sV i (j + 1)) := by
  obtain ⟨hc, hcV, hcap, hcI, hcO, -, horig, -⟩ :=
    t.child_S (g := g) hwf hsep hi hS hj (t.sC_getElem? hwf hi hS hj)
  refine ⟨hc, hcV, ?_⟩
  have hp : t.toSpqrTree.neOrig (t.toSpqrTree.neRange (t.sC i)[j]!).1 = some (t.sV i j, t.sV i (j + 1)) := by
    rw [← horig, t.neOrig_S (g := g) hwf hsep hi hS (k := j + 1) (by omega) (by omega)]
    simp
  exact t.capped_child g hwf hne hsep hi hc hcI hcO hcap hp s h


/-! ### Node-vertices and `V` items of an `S` node -/

theorem nodeVerts_nv (hwf : t.toSpqrTree.WF) (hi : i < t.size) {nv : Nat}
    (h : t.toSpqrTree.NvOf i nv) : ∃ d, t.nodeVerts[nv]? = some d ∧ d.node = i := by
  have := hwf.own.nv_node i nv hi h.1 h.2
  unfold SpqrTree.nodeOfNv at this
  exact Option.map_eq_some_iff.1 this

theorem not_capEnd_S (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hS : t.toSpqrTree.type i = .S)
    {j : Nat} (hj1 : 1 ≤ j) (hj : j < t.toSpqrTree.nVerts i - 1) :
    ¬ t.toSpqrTree.CapEnd i ((t.toSpqrTree.nvRange i).1 + j) := by
  rintro ⟨ne, d, hne, hd, hend⟩
  rw [t.capNe_S hS] at hne
  cases hne
  have h0 := t.nvs_S hwf hi hS (k := 0) (by omega)
  rw [Nat.add_zero, getElem!_of_getElem? hd, if_pos rfl] at h0
  rw [h0] at hend
  simp only at hend
  omega

/-- The `V` item of an interior node-vertex `s + j` (`1 ≤ j ≤ n - 2`) of an `S` node. -/
theorem sW_S (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hS : t.toSpqrTree.type i = .S) {j : Nat} (hj1 : 1 ≤ j) (hj : j < t.toSpqrTree.nVerts i - 1) :
    (∃ d, t.nodeVerts[(t.toSpqrTree.nvRange i).1 + j]? = some d ∧
      d.vert = t.sW ((t.toSpqrTree.nvRange i).1 + j)) ∧
    t.sW ((t.toSpqrTree.nvRange i).1 + j) ∈ t.children i ∧
    t.toSpqrTree.type (t.sW ((t.toSpqrTree.nvRange i).1 + j)) = .V ∧
    t.origId[t.sW ((t.toSpqrTree.nvRange i).1 + j)]! = some (t.sV i j) := by
  have hnv : t.toSpqrTree.NvOf i ((t.toSpqrTree.nvRange i).1 + j) := by
    constructor <;> unfold SpqrTree.nVerts at hj <;> omega
  obtain ⟨d, hd, -⟩ := t.nodeVerts_nv hwf hi hnv
  have hdv : d.vert = t.sW ((t.toSpqrTree.nvRange i).1 + j) := by
    unfold sW; rw [getElem!_of_getElem? hd]
  have hpar := (hsep.nv_child i _ d hi (Or.inl hS) hnv hd).2 (t.not_capEnd_S hwf hi hS hj1 hj)
  rw [hdv] at hpar
  have hws : t.sW ((t.toSpqrTree.nvRange i).1 + j) < t.size := by
    unfold SpqrTree.parent at hpar
    cases hp : t.par[t.sW ((t.toSpqrTree.nvRange i).1 + j)]? with
    | none => rw [hp] at hpar; simp at hpar
    | some _ => exact lt_of_lt_of_eq (Array.getElem?_eq_some_iff.1 hp).1 hwf.sizes.par
  have hbd : (t.toSpqrTree.nvRange i).1 + j < t.nodeVerts.size := (Array.getElem?_eq_some_iff.1 hd).1
  refine ⟨⟨d, hd, hdv⟩, t.mem_children_of_parent hwf hi hws hpar, hwf.own.nv_vert _ hbd, ?_⟩
  have ho := t.nvOrig_S (g := g) hwf hsep hi hS (j := j) (by omega)
  unfold SpqrTree.nvOrig at ho
  rw [hd] at ho
  have ho2 : t.origId[d.vert]?.getD none = some (t.sV i j) := ho
  rw [hdv] at ho2
  cases ho' : t.origId[t.sW ((t.toSpqrTree.nvRange i).1 + j)]? with
  | none => rw [ho'] at ho2; simp at ho2
  | some o =>
    rw [ho'] at ho2
    simp only [Option.getD_some] at ho2
    rw [getElem!_of_getElem? ho', ho2]

/-- The lower boundary of the `V` item at interior node-vertex `s + j`. -/
theorem sW_boundary (hwf : t.toSpqrTree.WF) (hne : t.ne = g.ne) (hsep : t.toSpqrTree.PieceSep g)
    (hi : i < t.size) (hS : t.toSpqrTree.type i = .S) {j : Nat} (hj1 : 1 ≤ j)
    (hj : j < t.toSpqrTree.nVerts i - 1) (s : EmbedState) (h : t.GluedUpTo g (i + 1) s) :
    (s.outerE[t.sW ((t.toSpqrTree.nvRange i).1 + j)]![0]! = none ∧
      t.edgesBelow (t.sW ((t.toSpqrTree.nvRange i).1 + j)) = []) ∨
    ∃ w0 w1 ρ, s.outerE[t.sW ((t.toSpqrTree.nvRange i).1 + j)]![0]! = some w0 ∧
      s.outerE[t.sW ((t.toSpqrTree.nvRange i).1 + j)]![1]! = some w1 ∧
      (t.pieceBelow g (t.sW ((t.toSpqrTree.nvRange i).1 + j))).OpenEmbedding s.rotAdj ρ w0 w1 (t.sV i j) := by
  obtain ⟨-, hw, htw, ho⟩ := t.sW_S (g := g) hwf hsep hi hS hj1 hj
  exact t.q_lower_boundary g hwf hne hi hw htw ho s h

/-! ### Attachment of the children at node-vertices -/

theorem nvInc_C (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hS : t.toSpqrTree.type i = .S) {j c nv : Nat} (hj : j < t.toSpqrTree.nVerts i - 1)
    (hc : (t.sC i)[j]? = some c) (h : t.toSpqrTree.NvInc i c nv) :
    nv = (t.toSpqrTree.nvRange i).1 + j ∨ nv = (t.toSpqrTree.nvRange i).1 + j + 1 := by
  obtain ⟨-, hcV, hcap, -, -, htw, -, -⟩ := t.child_S (g := g) hwf hsep hi hS hj hc
  rcases h with ⟨d, hd, hdv⟩ | ⟨ne, tw, d, -, htw', hcne, hd, hend⟩
  · exfalso
    have := hwf.own.nv_vert nv (Array.getElem?_eq_some_iff.1 hd).1
    rw [getElem!_of_getElem? hd, hdv] at this
    exact hcV this
  · have hcne' : t.toSpqrTree.capNe c = some (t.toSpqrTree.neRange c).1 := by
      simp only [SpqrTree.capNe, hcap, ↓reduceIte]
    rw [hcne'] at hcne
    cases hcne
    have h1 := hwf.twins.twin_invol _ _ htw'
    have h2 := hwf.twins.twin_invol _ _ htw
    rw [h1] at h2
    cases h2
    have hn := t.nvs_S hwf hi hS (k := j + 1) (by omega)
    rw [getElem!_of_getElem? hd, if_neg (Nat.succ_ne_zero j)] at hn
    rw [hn] at hend
    simp only at hend
    omega

theorem nvInc_W (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hS : t.toSpqrTree.type i = .S) {j nv : Nat} (hj1 : 1 ≤ j) (hj : j < t.toSpqrTree.nVerts i - 1)
    (hnv : t.toSpqrTree.NvOf i nv) (h : t.toSpqrTree.NvInc i (t.sW ((t.toSpqrTree.nvRange i).1 + j)) nv) :
    nv = (t.toSpqrTree.nvRange i).1 + j := by
  obtain ⟨⟨d', hd', hdv'⟩, -, htw, -⟩ := t.sW_S (g := g) hwf hsep hi hS hj1 hj
  rcases h with ⟨d, hd, hdv⟩ | ⟨ne, tw, d, -, -, hcne, -, -⟩
  · have hnv' : t.toSpqrTree.NvOf i ((t.toSpqrTree.nvRange i).1 + j) := by
      constructor <;> unfold SpqrTree.nVerts at hj <;> omega
    exact hsep.nv_vert_inj i nv _ d d' hi (Or.inl hS) hnv hnv' hd hd' (hdv.trans hdv'.symm)
  · exfalso
    unfold SpqrTree.capNe SpqrTree.hasCap at hcne
    rw [htw] at hcne
    simp [NodeType.isNode] at hcne

theorem mem_sL (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hS : t.toSpqrTree.type i = .S)
    {m c' : Nat} (hm : m ≤ t.toSpqrTree.nVerts i - 2) (h : c' ∈ t.sL i m) :
    (∃ j, j ≤ m ∧ (t.sC i)[j]? = some c') ∨
    (∃ j, 1 ≤ j ∧ j ≤ m ∧ c' = t.sW ((t.toSpqrTree.nvRange i).1 + j)) := by
  have h3 := (t.shape_S hwf hi hS).1
  unfold sL at h
  simp only [List.mem_cons, List.mem_flatMap, List.mem_range, List.not_mem_nil, or_false] at h
  rcases h with rfl | ⟨j, hj, rfl | rfl⟩
  · exact Or.inl ⟨0, Nat.zero_le _, t.sC_getElem? hwf hi hS (j := 0) (by omega)⟩
  · exact Or.inr ⟨j + 1, by omega, by omega, rfl⟩
  · exact Or.inl ⟨j + 1, by omega, t.sC_getElem? hwf hi hS (j := j + 1) (by omega)⟩

theorem mem_sL_children (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hS : t.toSpqrTree.type i = .S) {m c' : Nat} (hm : m ≤ t.toSpqrTree.nVerts i - 2)
    (h : c' ∈ t.sL i m) : c' ∈ t.children i := by
  have h3 := (t.shape_S hwf hi hS).1
  rcases t.mem_sL hwf hi hS hm h with ⟨j, hj, hc⟩ | ⟨j, hj1, hj, rfl⟩
  · exact (t.child_S (g := g) hwf hsep hi hS (j := j) (by omega) hc).1
  · exact (t.sW_S (g := g) hwf hsep hi hS hj1 (by omega)).2.1

theorem nvInc_sL (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hS : t.toSpqrTree.type i = .S) {m c' nv : Nat} (hm : m ≤ t.toSpqrTree.nVerts i - 2)
    (hc' : c' ∈ t.sL i m) (hnv : t.toSpqrTree.NvOf i nv) (h : t.toSpqrTree.NvInc i c' nv) :
    nv ≤ (t.toSpqrTree.nvRange i).1 + m + 1 := by
  have h3 := (t.shape_S hwf hi hS).1
  rcases t.mem_sL hwf hi hS hm hc' with ⟨j, hj, hc⟩ | ⟨j, hj1, hj, rfl⟩
  · have := t.nvInc_C (g := g) hwf hsep hi hS (j := j) (by omega) hc h; omega
  · have := t.nvInc_W (g := g) hwf hsep hi hS hj1 (by omega) hnv h; omega

theorem sC_nodup (hwf : t.toSpqrTree.WF) (hi : i < t.size) : (t.sC i).Nodup :=
  (t.children_nodup hwf hi).filter _

theorem ne_sL_C (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hS : t.toSpqrTree.type i = .S) {m c' : Nat} (hm : m + 1 ≤ t.toSpqrTree.nVerts i - 2)
    (hc' : c' ∈ t.sL i m) : c' ≠ (t.sC i)[m + 1]! := by
  intro heq
  rcases t.mem_sL hwf hi hS (by omega) hc' with ⟨j, hj, hc⟩ | ⟨j, hj1, hj, rfl⟩
  · have hc2 := t.sC_getElem? hwf hi hS (j := m + 1) (by omega)
    rw [← heq, ← hc] at hc2
    have hl := t.sC_length hwf hi hS
    have := (List.Nodup.getElem?_inj (by omega) (t.sC_nodup hwf hi)).1 hc2
    omega
  · have h1 := (t.sW_S (g := g) hwf hsep hi hS hj1 (by omega)).2.2.1
    have h2 := (t.child_S (g := g) hwf hsep hi hS (j := m + 1) (by omega)
      (t.sC_getElem? hwf hi hS (by omega))).2.1
    rw [heq] at h1
    exact h2 h1

theorem ne_sL_W (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hS : t.toSpqrTree.type i = .S) {m c' : Nat} (hm : m + 1 ≤ t.toSpqrTree.nVerts i - 2)
    (hc' : c' ∈ t.sL i m) : c' ≠ t.sW ((t.toSpqrTree.nvRange i).1 + (m + 1)) := by
  intro heq
  obtain ⟨⟨d', hd', hdv'⟩, -, htw', -⟩ := t.sW_S (g := g) hwf hsep hi hS (j := m + 1) (by omega) (by omega)
  rcases t.mem_sL hwf hi hS (by omega) hc' with ⟨j, hj, hc⟩ | ⟨j, hj1, hj, rfl⟩
  · have h2 := (t.child_S (g := g) hwf hsep hi hS (j := j) (by omega) hc).2.1
    rw [heq] at h2
    exact h2 htw'
  · obtain ⟨⟨d, hd, hdv⟩, -, -, -⟩ := t.sW_S (g := g) hwf hsep hi hS hj1 (by omega)
    have hnv : t.toSpqrTree.NvOf i ((t.toSpqrTree.nvRange i).1 + j) := by
      constructor <;> unfold SpqrTree.nVerts at hm <;> omega
    have hnv' : t.toSpqrTree.NvOf i ((t.toSpqrTree.nvRange i).1 + (m + 1)) := by
      constructor <;> unfold SpqrTree.nVerts at hm <;> omega
    have := hsep.nv_vert_inj i _ _ d d' hi (Or.inl hS) hnv hnv' hd hd' (by rw [hdv, hdv', heq])
    omega

theorem W_ne_C (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hS : t.toSpqrTree.type i = .S) {j : Nat} (hj1 : 1 ≤ j) (hj : j < t.toSpqrTree.nVerts i - 1) :
    t.sW ((t.toSpqrTree.nvRange i).1 + j) ≠ (t.sC i)[j]! := by
  intro heq
  have h1 := (t.sW_S (g := g) hwf hsep hi hS hj1 hj).2.2.1
  have h2 := (t.child_S (g := g) hwf hsep hi hS hj (t.sC_getElem? hwf hi hS hj)).2.1
  rw [heq] at h1
  exact h2 h1

/-- Two distinct children of an `S` node touching the same vertex `x` are both attached at a
node-vertex with original `x`. -/
theorem attach_S (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size) (hS : t.toSpqrTree.type i = .S)
    {a b x : Nat} (ha : a ∈ t.children i) (hb : b ∈ t.children i) (hab : a ≠ b)
    (hta : t.toSpqrTree.Touches g a x) (htb : t.toSpqrTree.Touches g b x) :
    ∃ nv, t.toSpqrTree.NvOf i nv ∧ t.toSpqrTree.nvOrig nv = some x ∧
      t.toSpqrTree.NvInc i a nv ∧ t.toSpqrTree.NvInc i b nv := by
  rw [← t.children_eq] at ha hb
  obtain ⟨nv, hnv, ho⟩ := hsep.node_attach i a b x hi (Or.inl hS) ha hb hab hta htb
  exact ⟨nv, hnv, ho, hsep.node_touch i a x nv hi (Or.inl hS) ha hta hnv ho,
    hsep.node_touch i b x nv hi (Or.inl hS) hb htb hnv ho⟩

theorem hasEdge_chain {m x : Nat} (h : HasEdge (t.chain g i m).es x) :
    ∃ c' ∈ t.sL i m, HasEdge (t.pieceBelow g c').es x := by
  obtain ⟨p, hp, hx⟩ := h
  simp only [Piece.es, chain, List.map_flatMap, List.mem_flatMap] at hp
  obtain ⟨c', hc', hp⟩ := hp
  exact ⟨c', hc', p, hp, hx⟩

theorem sep_sL (hwf : t.toSpqrTree.WF) (hne : t.ne = g.ne) (hsep : t.toSpqrTree.PieceSep g)
    (hi : i < t.size) (hS : t.toSpqrTree.type i = .S) {m c x : Nat}
    (hm : m + 1 ≤ t.toSpqrTree.nVerts i - 2) (hc : c ∈ t.children i)
    (hcL : ∀ c' ∈ t.sL i m, c' ≠ c)
    (hcnv : ∀ nv, t.toSpqrTree.NvOf i nv → t.toSpqrTree.NvInc i c nv →
      (t.toSpqrTree.nvRange i).1 + (m + 1) ≤ nv)
    (hx : HasEdge (t.chain g i m).es x) (hxc : HasEdge (t.pieceBelow g c).es x) :
    x = t.sV i (m + 1) := by
  obtain ⟨c', hc', hx'⟩ := t.hasEdge_chain hx
  obtain ⟨nv, hnv, ho, h1, h2⟩ := t.attach_S hsep hi hS
    (t.mem_sL_children (g := g) hwf hsep hi hS (by omega) hc') hc (hcL c' hc')
    (t.touches_of_hasEdge hwf g hne hx') (t.touches_of_hasEdge hwf g hne hxc)
  have := t.nvInc_sL (g := g) hwf hsep hi hS (by omega) hc' hnv h1
  have := hcnv nv hnv h2
  have heq : nv = (t.toSpqrTree.nvRange i).1 + (m + 1) := by omega
  rw [heq, t.nvOrig_S (g := g) hwf hsep hi hS (by omega)] at ho
  exact (Option.some.inj ho).symm

theorem sep_chain_C (hwf : t.toSpqrTree.WF) (hne : t.ne = g.ne) (hsep : t.toSpqrTree.PieceSep g)
    (hi : i < t.size) (hS : t.toSpqrTree.type i = .S) {m : Nat}
    (hm : m + 1 ≤ t.toSpqrTree.nVerts i - 2) :
    ∀ x, HasEdge (t.chain g i m).es x → HasEdge (t.pieceBelow g (t.sC i)[m + 1]!).es x →
      x = t.sV i (m + 1) := by
  intro x hx hxc
  have hc := t.sC_getElem? hwf hi hS (j := m + 1) (by omega)
  refine t.sep_sL hwf hne hsep hi hS hm (t.child_S (g := g) hwf hsep hi hS (j := m + 1) (by omega) hc).1
    (fun c' hc' => t.ne_sL_C (g := g) hwf hsep hi hS hm hc') ?_ hx hxc
  intro nv _ h
  have := t.nvInc_C (g := g) hwf hsep hi hS (j := m + 1) (by omega) hc h
  omega

theorem sep_chain_W (hwf : t.toSpqrTree.WF) (hne : t.ne = g.ne) (hsep : t.toSpqrTree.PieceSep g)
    (hi : i < t.size) (hS : t.toSpqrTree.type i = .S) {m : Nat}
    (hm : m + 1 ≤ t.toSpqrTree.nVerts i - 2) :
    ∀ x, HasEdge (t.chain g i m).es x →
      HasEdge (t.pieceBelow g (t.sW ((t.toSpqrTree.nvRange i).1 + (m + 1)))).es x →
      x = t.sV i (m + 1) := by
  intro x hx hxc
  refine t.sep_sL hwf hne hsep hi hS hm (t.sW_S (g := g) hwf hsep hi hS (by omega) (by omega)).2.1
    (fun c' hc' => t.ne_sL_W (g := g) hwf hsep hi hS hm hc') ?_ hx hxc
  intro nv hnv h
  have := t.nvInc_W (g := g) hwf hsep hi hS (by omega) (by omega) hnv h
  omega

theorem sep_W_C (hwf : t.toSpqrTree.WF) (hne : t.ne = g.ne) (hsep : t.toSpqrTree.PieceSep g)
    (hi : i < t.size) (hS : t.toSpqrTree.type i = .S) {j : Nat} (hj1 : 1 ≤ j)
    (hj : j < t.toSpqrTree.nVerts i - 1) :
    ∀ x, HasEdge (t.pieceBelow g (t.sW ((t.toSpqrTree.nvRange i).1 + j))).es x →
      HasEdge (t.pieceBelow g (t.sC i)[j]!).es x → x = t.sV i j := by
  intro x hxw hxc
  have hc := t.sC_getElem? hwf hi hS hj
  obtain ⟨nv, hnv, ho, h1, h2⟩ := t.attach_S hsep hi hS
    (t.sW_S (g := g) hwf hsep hi hS hj1 hj).2.1 (t.child_S (g := g) hwf hsep hi hS hj hc).1
    (t.W_ne_C (g := g) hwf hsep hi hS hj1 hj)
    (t.touches_of_hasEdge hwf g hne hxw) (t.touches_of_hasEdge hwf g hne hxc)
  have := t.nvInc_W (g := g) hwf hsep hi hS hj1 hj hnv h1
  rw [this, t.nvOrig_S (g := g) hwf hsep hi hS (by omega)] at ho
  exact (Option.some.inj ho).symm

/-! ### Disjointness of the glued pieces -/

theorem disjoint_sL (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size) (hS : t.toSpqrTree.type i = .S)
    {m c : Nat} (hm : m ≤ t.toSpqrTree.nVerts i - 2) (hc : c ∈ t.children i)
    (hcL : ∀ c' ∈ t.sL i m, c' ≠ c) : List.Disjoint (t.chain g i m).ves (t.edgesBelow c) := by
  rw [List.disjoint_left]
  intro e he hec
  simp only [chain, List.mem_flatMap] at he
  obtain ⟨c', hc', he⟩ := he
  exact List.disjoint_left.1 (t.maximal_pieces_disjoint hwf hsh
    (t.child_maximal hwf hi (t.mem_sL_children (g := g) hwf hsep hi hS hm hc'))
    (t.child_maximal hwf hi hc) (hcL c' hc')) he hec

theorem disjoint_chain_C (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size) (hS : t.toSpqrTree.type i = .S)
    {m : Nat} (hm : m + 1 ≤ t.toSpqrTree.nVerts i - 2) :
    List.Disjoint (t.chain g i m).ves (t.edgesBelow (t.sC i)[m + 1]!) :=
  t.disjoint_sL (g := g) hwf hsh hsep hi hS (by omega)
    (t.child_S (g := g) hwf hsep hi hS (by omega) (t.sC_getElem? hwf hi hS (by omega))).1
    (fun c' hc' => t.ne_sL_C (g := g) hwf hsep hi hS hm hc')

theorem disjoint_chain_W (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size) (hS : t.toSpqrTree.type i = .S)
    {m : Nat} (hm : m + 1 ≤ t.toSpqrTree.nVerts i - 2) :
    List.Disjoint (t.chain g i m).ves (t.edgesBelow (t.sW ((t.toSpqrTree.nvRange i).1 + (m + 1)))) :=
  t.disjoint_sL (g := g) hwf hsh hsep hi hS (by omega)
    (t.sW_S (g := g) hwf hsep hi hS (by omega) (by omega)).2.1
    (fun c' hc' => t.ne_sL_W (g := g) hwf hsep hi hS hm hc')

theorem disjoint_W_C (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size) (hS : t.toSpqrTree.type i = .S)
    {j : Nat} (hj1 : 1 ≤ j) (hj : j < t.toSpqrTree.nVerts i - 1) :
    List.Disjoint (t.edgesBelow (t.sW ((t.toSpqrTree.nvRange i).1 + j))) (t.edgesBelow (t.sC i)[j]!) :=
  t.maximal_pieces_disjoint hwf hsh
    (t.child_maximal hwf hi (t.sW_S (g := g) hwf hsep hi hS hj1 hj).2.1)
    (t.child_maximal hwf hi (t.child_S (g := g) hwf hsep hi hS hj (t.sC_getElem? hwf hi hS hj)).1)
    (t.W_ne_C (g := g) hwf hsep hi hS hj1 hj)

end S

end Spqr.PlanarSpqrTree
