import Spqr.PlanarEmbedNodeSChain

/-!
# The `S` node: the node fold block by block

The node loop of an `S` node with `n` node-vertices runs over `n` blocks of four quarter-edges
(`blockFold`): block `0` (the cap) exposes the outer slots of the first and last non-`V` children,
block `k` (`1 ≤ k ≤ n - 2`) glues the next `V` item and non-`V` child at node-vertex `s + k`
(`Capped.attachOpen` + `Capped.join`), and block `n - 1` is a no-op.
-/

namespace Spqr.PlanarSpqrTree

open EmbedM

variable (t : PlanarSpqrTree)

theorem getElem!_eq_of_getElem?_eq {α : Type} [Inhabited α] {a b : Array α} {j : Nat}
    (h : a[j]? = b[j]?) : a[j]! = b[j]! := by
  by_cases ha : j < a.size
  · rw [Array.getElem?_eq_getElem ha] at h
    obtain ⟨hb, hv⟩ := Array.getElem?_eq_some_iff.1 h.symm
    rw [getElem!_pos a j ha, getElem!_pos b j hb, hv]
  · rw [Array.getElem?_eq_none_iff.2 (by omega)] at h
    have hb : ¬ j < b.size := by have := Array.getElem?_eq_none_iff.1 h.symm; omega
    rw [getElem!_neg a j ha, getElem!_neg b j hb]

theorem setOuter_outerE_eq' (i k : Nat) (q : Option Nat) (s : EmbedState) {r : Array (Option Nat)}
    (h : s.outerE[i]? = some r) : ((setOuter i k q).run s).2.outerE[i]? = some (r.set! k q) := by
  obtain ⟨hi, hv⟩ := Array.getElem?_eq_some_iff.1 h
  rw [setOuter_outerE_eq _ _ _ _ hi, hv]

theorem setOuter_row_ne (i k : Nat) (q : Option Nat) (s : EmbedState) {j : Nat} (hj : j ≠ i) (r : Nat) :
    ((setOuter i k q).run s).2.outerE[j]![r]! = s.outerE[j]![r]! := by
  rw [getElem!_eq_of_getElem?_eq (setOuter_outerE_ne i k q s j hj)]

theorem ite_perm4 {α : Type} (q a b c d : Nat) (x y z w f : α) (hab : a ≠ b) (hcd : c ≠ d)
    (hac : a ≠ c) (had : a ≠ d) (hbc : b ≠ c) (hbd : b ≠ d) :
    (if q = c then z else if q = d then w else if q = a then x else if q = b then y else f) =
      (if q = a then x else if q = b then y else if q = d then w else if q = c then z else f) := by
  split_ifs <;> first | rfl | omega

theorem ite_swap2 {α : Type} (q a b : Nat) (x y f : α) (hab : a ≠ b) :
    (if q = a then x else if q = b then y else f) = (if q = b then y else if q = a then x else f) := by
  split_ifs <;> first | rfl | omega

theorem array_four {α : Type} (r : Array α) (h : r.size = 4) :
    r = #[r[0], r[1], r[2], r[3]] := by
  apply Array.ext
  · simp [h]
  · intro i h1 h2
    simp at h2
    interval_cases i <;> rfl

theorem nodeStep_outerE_ge (i neSt : Nat) (s : EmbedState) {ta : Nat} (hta : 4 * (neSt + 1) ≤ ta) :
    (t.nodeStep i neSt s ta).outerE = s.outerE := by
  unfold nodeStep
  cases t.neRotAdj[ta]! with
  | none => rfl
  | some tb =>
    by_cases h1 : tb < ta
    · simp [h1]
    · have h2 : ¬ ta < 4 * (neSt + 1) := by omega
      simp only [h1, h2, if_false]
      split
      · split <;> simp [link_outerE]
      · simp [link_outerE]

/-- The state after blocks `0..m` of the node loop. -/
def blockFold (i neSt : Nat) (s : EmbedState) (m : Nat) : EmbedState :=
  (List.range' (4 * neSt) (4 * (m + 1))).foldl (t.nodeStep i neSt) s

theorem range'_four (a : Nat) : List.range' a 4 = [a, a + 1, a + 2, a + 3] := by
  simp [List.range'_succ]

theorem blockFold_succ (i neSt : Nat) (s : EmbedState) (m : Nat) :
    t.blockFold i neSt s (m + 1) =
      (List.range' (4 * (neSt + (m + 1))) 4).foldl (t.nodeStep i neSt) (t.blockFold i neSt s m) := by
  unfold blockFold
  rw [← List.foldl_append, show 4 * (neSt + (m + 1)) = 4 * neSt + 1 * (4 * (m + 1)) by omega,
    List.range'_append, show 4 * (m + 1) + 4 = 4 * (m + 1 + 1) by omega]

theorem blockFold_outerE_succ (i neSt : Nat) (s : EmbedState) (m : Nat) :
    (t.blockFold i neSt s (m + 1)).outerE = (t.blockFold i neSt s m).outerE := by
  rw [blockFold_succ, range'_four]
  simp only [List.foldl_cons, List.foldl_nil]
  rw [t.nodeStep_outerE_ge _ _ _ (by omega), t.nodeStep_outerE_ge _ _ _ (by omega),
    t.nodeStep_outerE_ge _ _ _ (by omega), t.nodeStep_outerE_ge _ _ _ (by omega)]

theorem blockFold_outerE (i neSt : Nat) (s : EmbedState) (m : Nat) :
    (t.blockFold i neSt s m).outerE = (t.blockFold i neSt s 0).outerE := by
  induction m with
  | zero => rfl
  | succ m ih => rw [blockFold_outerE_succ, ih]

theorem blockFold_outerE_ne (i neSt : Nat) (s : EmbedState) (m : Nat) {j : Nat} (hj : j ≠ i) :
    (t.blockFold i neSt s m).outerE[j]? = s.outerE[j]? :=
  t.nodeFold_outerE_ne i neSt _ s j hj

theorem blockFold_rotAdj_size (i neSt : Nat) (s : EmbedState) (m : Nat) :
    (t.blockFold i neSt s m).rotAdj.size = s.rotAdj.size :=
  t.nodeFold_rotAdj_size i neSt _ s

section S

variable {i : Nat}

theorem rotS_S (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hS : t.toSpqrTree.type i = .S)
    (hlay : ∃ ev mr cv, t.neRotAdj.extract (4 * (t.toSpqrTree.neRange i).1) (4 * (t.toSpqrTree.neRange i).2) =
      layoutRot (t.toSpqrTree.type i) (t.toSpqrTree.nVerts i) (t.toSpqrTree.neRange i).1
        (t.toSpqrTree.neRange i).2 ev mr cv)
    {k r : Nat} (hk : k < t.toSpqrTree.nVerts i) (hr : r < 4) {ta tb : Nat}
    (hta : ta = 4 * ((t.toSpqrTree.neRange i).1 + k) + r)
    (htb : tb = rotS (t.toSpqrTree.nVerts i) (t.toSpqrTree.neRange i).1 k r) :
    t.neRotAdj[ta]! = some tb := by
  subst hta htb; exact t.neRotAdj_S hwf hi hS hlay hk hr

/-- Invariant of the `S` fold after blocks `0..m` (`m ≤ n - 2`): rows of other items unchanged,
`rotAdj` unchanged outside `chain m`, and `chain m` is `Capped` at `sV 0`, `sV (m + 1)` with the
first child's slots `0, 1` and the `m`-th child's slots `2, 3`. -/
structure InvS (g : Graph) (i : Nat) (s s' : EmbedState) (m : Nat) : Prop where
  outer_ne : ∀ j, j ≠ i → s'.outerE[j]? = s.outerE[j]?
  rot_size : s'.rotAdj.size = s.rotAdj.size
  frame : ∀ q, ¬ (t.chain g i m).Mem q → s'.rotAdj[q]? = s.rotAdj[q]?
  capped : ∃ a0 a1 a2 a3 ρ, s.outerE[(t.sC i)[0]!]![0]! = some a0 ∧
    s.outerE[(t.sC i)[0]!]![1]! = some a1 ∧
    s.outerE[(t.sC i)[m]!]![2]! = some a2 ∧ s.outerE[(t.sC i)[m]!]![3]! = some a3 ∧
    (t.chain g i m).Capped s'.rotAdj ρ a0 a1 a2 a3 (t.sV i 0) (t.sV i (m + 1))

/-- Block `0` (the cap edge): four `setOuter`s exposing the first child's slots `0, 1` and the last
child's slots `2, 3`; `rotAdj` untouched. -/
theorem cap_block_S (hwf : t.toSpqrTree.WF) (hne : t.ne = g.ne) (hsep : t.toSpqrTree.PieceSep g)
    (hi : i < t.size) (hS : t.toSpqrTree.type i = .S)
    (hlay : ∃ ev mr cv, t.neRotAdj.extract (4 * (t.toSpqrTree.neRange i).1) (4 * (t.toSpqrTree.neRange i).2) =
      layoutRot (t.toSpqrTree.type i) (t.toSpqrTree.nVerts i) (t.toSpqrTree.neRange i).1
        (t.toSpqrTree.neRange i).2 ev mr cv)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    (t.blockFold i (t.toSpqrTree.neRange i).1 s 0).rotAdj = s.rotAdj ∧
    ∃ a0 a1 e2 e3, s.outerE[(t.sC i)[0]!]![0]! = some a0 ∧ s.outerE[(t.sC i)[0]!]![1]! = some a1 ∧
      s.outerE[(t.sC i)[t.toSpqrTree.nVerts i - 2]!]![2]! = some e2 ∧
      s.outerE[(t.sC i)[t.toSpqrTree.nVerts i - 2]!]![3]! = some e3 ∧
      (t.blockFold i (t.toSpqrTree.neRange i).1 s 0).outerE[i]? =
        some #[some a0, some a1, some e2, some e3] := by
  have h3 := (t.shape_S hwf hi hS).1
  obtain ⟨hc0, -, a0, a1, -, -, -, ha0, ha1, -, -, -⟩ :=
    t.child_cert_S hwf hne hsep hi hS (j := 0) (by omega) s h
  obtain ⟨hcl, -, -, -, e2, e3, -, -, -, he2, he3, -⟩ :=
    t.child_cert_S hwf hne hsep hi hS (j := t.toSpqrTree.nVerts i - 2) (by omega) s h
  have hc0i : (t.sC i)[0]! ≠ i := fun e => by have := (t.child_data hwf hi hc0).2.2; omega
  have hcli : (t.sC i)[t.toSpqrTree.nVerts i - 2]! ≠ i := fun e => by
    have := (t.child_data hwf hi hcl).2.2; omega
  obtain ⟨-, -, -, -, -, -, -, hq0⟩ :=
    t.child_S (g := g) hwf hsep hi hS (j := 0) (by omega) (t.sC_getElem? hwf hi hS (by omega))
  obtain ⟨-, -, -, -, -, -, -, hql⟩ :=
    t.child_S (g := g) hwf hsep hi hS (j := t.toSpqrTree.nVerts i - 2) (by omega)
      (t.sC_getElem? hwf hi hS (by omega))
  have hr0 := t.rotS_S hwf hi hS hlay (k := 0) (r := 0) (by omega) (by omega)
    (ta := 4 * (t.toSpqrTree.neRange i).1) (tb := 4 * ((t.toSpqrTree.neRange i).1 + (0 + 1)) + 1)
    (by omega) (by simp [rotS] <;> omega)
  have hr1 := t.rotS_S hwf hi hS hlay (k := 0) (r := 1) (by omega) (by omega)
    (ta := 4 * (t.toSpqrTree.neRange i).1 + 1) (tb := 4 * ((t.toSpqrTree.neRange i).1 + (0 + 1)) + 0)
    (by omega) (by simp [rotS] <;> omega)
  have hr2 := t.rotS_S hwf hi hS hlay (k := 0) (r := 2) (by omega) (by omega)
    (ta := 4 * (t.toSpqrTree.neRange i).1 + 2)
    (tb := 4 * ((t.toSpqrTree.neRange i).1 + (t.toSpqrTree.nVerts i - 2 + 1)) + 3)
    (by omega) (by simp [rotS] <;> omega)
  have hr3 := t.rotS_S hwf hi hS hlay (k := 0) (r := 3) (by omega) (by omega)
    (ta := 4 * (t.toSpqrTree.neRange i).1 + 3)
    (tb := 4 * ((t.toSpqrTree.neRange i).1 + (t.toSpqrTree.nVerts i - 2 + 1)) + 2)
    (by omega) (by simp [rotS] <;> omega)
  have hsl0 : 2 * QE.side (4 * (t.toSpqrTree.neRange i).1) + (1 - QE.dir (4 * (t.toSpqrTree.neRange i).1)) = 1 := by
    unfold QE.side QE.dir; omega
  have hsl1 : 2 * QE.side (4 * (t.toSpqrTree.neRange i).1 + 1) + (1 - QE.dir (4 * (t.toSpqrTree.neRange i).1 + 1)) = 0 := by
    unfold QE.side QE.dir; omega
  have hsl2 : 2 * QE.side (4 * (t.toSpqrTree.neRange i).1 + 2) + (1 - QE.dir (4 * (t.toSpqrTree.neRange i).1 + 2)) = 3 := by
    unfold QE.side QE.dir; omega
  have hsl3 : 2 * QE.side (4 * (t.toSpqrTree.neRange i).1 + 3) + (1 - QE.dir (4 * (t.toSpqrTree.neRange i).1 + 3)) = 2 := by
    unfold QE.side QE.dir; omega
  have hsz : i < s.outerE.size := by rw [h.outer_size]; exact hi
  obtain ⟨r, hrow, hr4⟩ : ∃ r, s.outerE[i]? = some r ∧ r.size = 4 :=
    ⟨_, Array.getElem?_eq_getElem hsz, by rw [← getElem!_pos s.outerE i hsz]; exact h.outer_row_size i hi⟩
  unfold blockFold
  rw [show 4 * (0 + 1) = 4 by rfl, range'_four]
  simp only [List.foldl_cons, List.foldl_nil]
  rw [t.nodeStep_cap i _ s hr0 (by omega) (by omega), hsl0, hq0 1 (by omega)]
  rw [t.nodeStep_cap i _ _ hr1 (by omega) (by omega), hsl1, hq0 0 (by omega), setOuter_row_ne _ _ _ _ hc0i]
  rw [t.nodeStep_cap i _ _ hr2 (by omega) (by omega), hsl2, hql 3 (by omega),
    setOuter_row_ne _ _ _ _ hcli, setOuter_row_ne _ _ _ _ hcli]
  rw [t.nodeStep_cap i _ _ hr3 (by omega) (by omega), hsl3, hql 2 (by omega),
    setOuter_row_ne _ _ _ _ hcli, setOuter_row_ne _ _ _ _ hcli, setOuter_row_ne _ _ _ _ hcli]
  rw [ha0, ha1, he2, he3]
  refine ⟨by simp only [setOuter_rotAdj], a0, a1, e2, e3, rfl, rfl, rfl, rfl, ?_⟩
  rw [setOuter_outerE_eq' _ _ _ _ (setOuter_outerE_eq' _ _ _ _ (setOuter_outerE_eq' _ _ _ _
    (setOuter_outerE_eq' _ _ _ _ hrow)))]
  rw [array_four r hr4]
  rfl


theorem treeQe_congr {s s' : EmbedState} (h : s'.outerE = s.outerE) (q : Nat) :
    t.treeQe s' q = t.treeQe s q := by
  unfold treeQe; rw [h]

theorem mem_chain_iff (g : Graph) (m q : Nat) :
    (t.chain g i m).Mem q ↔ ∃ c' ∈ t.sL i m, (t.pieceBelow g c').Mem q := by
  simp only [Piece.Mem, chain, pieceBelow, List.mem_flatMap]

theorem mem_chain_bound (hwf : t.toSpqrTree.WF) (g : Graph) {m q : Nat} (hq : (t.chain g i m).Mem q) :
    q < 4 * t.ne := by
  obtain ⟨c', -, hq⟩ := (t.mem_chain_iff g m q).1 hq
  exact t.mem_pieceBelow_bound hwf g hq

theorem sV_inj (hwf : t.toSpqrTree.WF) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hS : t.toSpqrTree.type i = .S) {j j' : Nat} (hj : j < t.toSpqrTree.nVerts i)
    (hj' : j' < t.toSpqrTree.nVerts i) (h : t.sV i j = t.sV i j') : j = j' := by
  have h1 := t.nvOrig_S (g := g) hwf hsep hi hS hj
  have h2 := t.nvOrig_S (g := g) hwf hsep hi hS hj'
  rw [h, ← h2] at h1
  have hnv : t.toSpqrTree.NvOf i ((t.toSpqrTree.nvRange i).1 + j) := by
    constructor <;> unfold SpqrTree.nVerts at hj <;> omega
  have hnv' : t.toSpqrTree.NvOf i ((t.toSpqrTree.nvRange i).1 + j') := by
    constructor <;> unfold SpqrTree.nVerts at hj' <;> omega
  have := hsep.nv_orig_inj i _ _ hi (Or.inl hS) hnv hnv' h1
  omega

theorem chain_zero_eq (g : Graph) : t.chain g i 0 = t.pieceBelow g (t.sC i)[0]! := by
  simp [chain, sL, pieceBelow]

theorem chain_succ_eq (g : Graph) (m : Nat) :
    t.chain g i (m + 1) =
      {t.chain g i m with ves := ((t.chain g i m).ves ++
        t.edgesBelow (t.sW ((t.toSpqrTree.nvRange i).1 + (m + 1)))) ++ t.edgesBelow (t.sC i)[m + 1]!} := by
  have := t.chain_succ_ves g i m
  rw [← List.append_assoc] at this
  exact congrArg (fun v => ({t.chain g i m with ves := v} : Piece)) this

theorem invS_zero (hwf : t.toSpqrTree.WF) (hne : t.ne = g.ne) (hsep : t.toSpqrTree.PieceSep g)
    (hi : i < t.size) (hS : t.toSpqrTree.type i = .S)
    (hlay : ∃ ev mr cv, t.neRotAdj.extract (4 * (t.toSpqrTree.neRange i).1) (4 * (t.toSpqrTree.neRange i).2) =
      layoutRot (t.toSpqrTree.type i) (t.toSpqrTree.nVerts i) (t.toSpqrTree.neRange i).1
        (t.toSpqrTree.neRange i).2 ev mr cv)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    t.InvS g i s (t.blockFold i (t.toSpqrTree.neRange i).1 s 0) 0 := by
  have h3 := (t.shape_S hwf hi hS).1
  obtain ⟨hrot, a0, a1, e2, e3, ha0, ha1, -, -, -⟩ := t.cap_block_S hwf hne hsep hi hS hlay s h
  obtain ⟨-, -, b0, b1, b2, b3, ρ, hb0, hb1, hb2, hb3, hcap⟩ :=
    t.child_cert_S hwf hne hsep hi hS (j := 0) (by omega) s h
  rw [ha0] at hb0; rw [ha1] at hb1; cases hb0; cases hb1
  refine ⟨fun j hj => t.blockFold_outerE_ne _ _ _ _ hj, by rw [hrot], fun q _ => by rw [hrot],
    a0, a1, b2, b3, ρ, ha0, ha1, hb2, hb3, ?_⟩
  rw [hrot, chain_zero_eq]
  exact hcap

set_option maxHeartbeats 1000000 in
/-- Block `m + 1` (`m + 1 ≤ n - 2`): the corner at node-vertex `s + (m + 1)` attaches the `V`
item (if it has edges) and joins the next child. -/
theorem invS_step (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape) (hne : t.ne = g.ne)
    (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size) (hS : t.toSpqrTree.type i = .S)
    (hlay : ∃ ev mr cv, t.neRotAdj.extract (4 * (t.toSpqrTree.neRange i).1) (4 * (t.toSpqrTree.neRange i).2) =
      layoutRot (t.toSpqrTree.type i) (t.toSpqrTree.nVerts i) (t.toSpqrTree.neRange i).1
        (t.toSpqrTree.neRange i).2 ev mr cv)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) {m : Nat} (hm : m + 1 ≤ t.toSpqrTree.nVerts i - 2)
    (hinv : t.InvS g i s (t.blockFold i (t.toSpqrTree.neRange i).1 s m) m) :
    t.InvS g i s (t.blockFold i (t.toSpqrTree.neRange i).1 s (m + 1)) (m + 1) := by
  have h3 := (t.shape_S hwf hi hS).1
  rw [blockFold_succ, range'_four]
  simp only [List.foldl_cons, List.foldl_nil]
  generalize t.blockFold i (t.toSpqrTree.neRange i).1 s m = s' at hinv
  obtain ⟨houter, hrsz, hframe, a0, a1, a2, a3, ρ, ha0, ha1, ha2, ha3, hcap⟩ := hinv
  have hcm := t.sC_getElem? hwf hi hS (j := m) (by omega)
  have hcm1 := t.sC_getElem? hwf hi hS (j := m + 1) (by omega)
  obtain ⟨hcm_mem, -, -, -, -, -, -, hqm⟩ := t.child_S (g := g) hwf hsep hi hS (j := m) (by omega) hcm
  obtain ⟨hcm1_mem, -, -, -, -, -, -, hqm1⟩ :=
    t.child_S (g := g) hwf hsep hi hS (j := m + 1) (by omega) hcm1
  obtain ⟨-, -, b0, b1, b2, b3, σ, hb0, hb1, hb2, hb3, hcapb⟩ :=
    t.child_cert_S hwf hne hsep hi hS (j := m + 1) (by omega) s h
  obtain ⟨-, hw_mem, -, -⟩ := t.sW_S (g := g) hwf hsep hi hS (j := m + 1) (by omega) (by omega)
  have hwi : t.sW ((t.toSpqrTree.nvRange i).1 + (m + 1)) ≠ i := fun e => by
    have := (t.child_data hwf hi hw_mem).2.2; omega
  have hcmi : (t.sC i)[m]! ≠ i := fun e => by have := (t.child_data hwf hi hcm_mem).2.2; omega
  have hcm1i : (t.sC i)[m + 1]! ≠ i := fun e => by have := (t.child_data hwf hi hcm1_mem).2.2; omega
  have hr0 := t.rotS_S hwf hi hS hlay (k := m + 1) (r := 0) (by omega) (by omega)
    (ta := 4 * ((t.toSpqrTree.neRange i).1 + (m + 1))) rfl rfl
  have hr1 := t.rotS_S hwf hi hS hlay (k := m + 1) (r := 1) (by omega) (by omega)
    (ta := 4 * ((t.toSpqrTree.neRange i).1 + (m + 1)) + 1) rfl rfl
  have hr2 := t.rotS_S hwf hi hS hlay (k := m + 1) (r := 2) (by omega) (by omega)
    (ta := 4 * ((t.toSpqrTree.neRange i).1 + (m + 1)) + 2)
    (tb := 4 * ((t.toSpqrTree.neRange i).1 + (m + 1 + 1)) + 1) rfl
    (by simp [rotS, show ¬ (m + 1 + 1 = t.toSpqrTree.nVerts i) by omega] <;> omega)
  have hr3 := t.rotS_S hwf hi hS hlay (k := m + 1) (r := 3) (by omega) (by omega)
    (ta := 4 * ((t.toSpqrTree.neRange i).1 + (m + 1)) + 3)
    (tb := 4 * ((t.toSpqrTree.neRange i).1 + (m + 1 + 1)) + 0) rfl
    (by simp [rotS, show ¬ (m + 1 + 1 = t.toSpqrTree.nVerts i) by omega] <;> omega)
  have hlt0 : rotS (t.toSpqrTree.nVerts i) (t.toSpqrTree.neRange i).1 (m + 1) 0 <
      4 * ((t.toSpqrTree.neRange i).1 + (m + 1)) := by
    simp only [rotS, Nat.succ_ne_zero, ite_false]; split_ifs <;> omega
  have hlt1 : rotS (t.toSpqrTree.nVerts i) (t.toSpqrTree.neRange i).1 (m + 1) 1 <
      4 * ((t.toSpqrTree.neRange i).1 + (m + 1)) + 1 := by
    simp only [rotS, Nat.succ_ne_zero, ite_false]; split_ifs <;> omega
  rw [t.nodeStep_skip _ _ _ hr0 hlt0, t.nodeStep_skip _ _ _ hr1 hlt1,
    t.nodeStep_corner _ _ _ hr2 (by omega) (by omega) (by omega) (by omega)]
  have hv : t.nodeVerts[(t.nodeEdges[QE.edge (4 * ((t.toSpqrTree.neRange i).1 + (m + 1)) + 2)]!).nvs.2]!.vert =
      t.sW ((t.toSpqrTree.nvRange i).1 + (m + 1)) := by
    have : QE.edge (4 * ((t.toSpqrTree.neRange i).1 + (m + 1)) + 2) = (t.toSpqrTree.neRange i).1 + (m + 1) := by
      unfold QE.edge; omega
    rw [this, t.nvs_S hwf hi hS (k := m + 1) (by omega), if_neg (Nat.succ_ne_zero m)]
    rfl
  simp only [hv]
  have hq2 : t.treeQe s' (4 * ((t.toSpqrTree.neRange i).1 + (m + 1)) + 2) = some a2 := by
    rw [hqm 2 (by omega), getElem!_eq_of_getElem?_eq (houter _ hcmi), ha2]
  have hq3 : t.treeQe s' (4 * ((t.toSpqrTree.neRange i).1 + (m + 1)) + 3) = some a3 := by
    rw [hqm 3 (by omega), getElem!_eq_of_getElem?_eq (houter _ hcmi), ha3]
  have hqb1 : t.treeQe s' (4 * ((t.toSpqrTree.neRange i).1 + (m + 1 + 1)) + 1) = some b1 := by
    rw [hqm1 1 (by omega), getElem!_eq_of_getElem?_eq (houter _ hcm1i), hb1]
  have hqb0 : t.treeQe s' (4 * ((t.toSpqrTree.neRange i).1 + (m + 1 + 1)) + 0) = some b0 := by
    rw [hqm1 0 (by omega), getElem!_eq_of_getElem?_eq (houter _ hcm1i), hb0]
  have hw0 : s'.outerE[t.sW ((t.toSpqrTree.nvRange i).1 + (m + 1))]![0]! =
      s.outerE[t.sW ((t.toSpqrTree.nvRange i).1 + (m + 1))]![0]! := by
    rw [getElem!_eq_of_getElem?_eq (houter _ hwi)]
  have hw1 : s'.outerE[t.sW ((t.toSpqrTree.nvRange i).1 + (m + 1))]![1]! =
      s.outerE[t.sW ((t.toSpqrTree.nvRange i).1 + (m + 1))]![1]! := by
    rw [getElem!_eq_of_getElem?_eq (houter _ hwi)]
  rw [hq2, hqb1, hw0, hw1]
  -- membership, bounds and distinctness of the linked slots
  have hmemA2 : (t.chain g i m).Mem a2 := by
    obtain ⟨l2, l3, hl2, -, -, -⟩ := hcap.pair2; exact Piece.mem_of_loc hl2
  have hmemA3 : (t.chain g i m).Mem a3 := by
    obtain ⟨l2, l3, -, hl3, -, -⟩ := hcap.pair2; exact Piece.mem_of_loc hl3
  have hmemB0 : (t.pieceBelow g (t.sC i)[m + 1]!).Mem b0 := by
    obtain ⟨l0, l1, hl0, -, -, -⟩ := hcapb.pair0; exact Piece.mem_of_loc hl0
  have hmemB1 : (t.pieceBelow g (t.sC i)[m + 1]!).Mem b1 := by
    obtain ⟨l0, l1, -, hl1, -, -⟩ := hcapb.pair0; exact Piece.mem_of_loc hl1
  have hbdA : ∀ q, (t.chain g i m).Mem q → q < s'.rotAdj.size := by
    intro q hq; rw [hrsz, h.rot_size]; exact t.mem_chain_bound hwf g hq
  have hbdC : ∀ q, (t.pieceBelow g (t.sC i)[m + 1]!).Mem q → q < s'.rotAdj.size := by
    intro q hq; rw [hrsz, h.rot_size]; exact t.mem_pieceBelow_bound hwf g hq
  have hbdW : ∀ q, (t.pieceBelow g (t.sW ((t.toSpqrTree.nvRange i).1 + (m + 1)))).Mem q →
      q < s'.rotAdj.size := by
    intro q hq; rw [hrsz, h.rot_size]; exact t.mem_pieceBelow_bound hwf g hq
  have hdisAC := t.disjoint_chain_C (g := g) hwf hsh hsep hi hS hm
  have hdisAW := t.disjoint_chain_W (g := g) hwf hsh hsep hi hS hm
  have hdisWC := t.disjoint_W_C (g := g) hwf hsh hsep hi hS (j := m + 1) (by omega) (by omega)
  have hAC : ∀ q, (t.chain g i m).Mem q → (t.pieceBelow g (t.sC i)[m + 1]!).Mem q → False :=
    fun q h1 h2 => List.disjoint_left.1 hdisAC h1 h2
  have hAW : ∀ q, (t.chain g i m).Mem q →
      (t.pieceBelow g (t.sW ((t.toSpqrTree.nvRange i).1 + (m + 1)))).Mem q → False :=
    fun q h1 h2 => List.disjoint_left.1 hdisAW h1 h2
  have hWC : ∀ q, (t.pieceBelow g (t.sW ((t.toSpqrTree.nvRange i).1 + (m + 1)))).Mem q →
      (t.pieceBelow g (t.sC i)[m + 1]!).Mem q → False :=
    fun q h1 h2 => List.disjoint_left.1 hdisWC h1 h2
  have ha23 : a2 ≠ a3 := by have := hcap.dir2; have := hcap.dir3; omega
  have hb01 : b0 ≠ b1 := by have := hcapb.dir0; have := hcapb.dir1; omega
  have ha2b0 : a2 ≠ b0 := fun e => hAC a2 hmemA2 (by rw [e]; exact hmemB0)
  have ha2b1 : a2 ≠ b1 := fun e => hAC a2 hmemA2 (by rw [e]; exact hmemB1)
  have ha3b0 : a3 ≠ b0 := fun e => hAC a3 hmemA3 (by rw [e]; exact hmemB0)
  have ha3b1 : a3 ≠ b1 := fun e => hAC a3 hmemA3 (by rw [e]; exact hmemB1)
  have huv : t.sV i 0 ≠ t.sV i (m + 1) := fun e => by
    have := t.sV_inj (g := g) hwf hsep hi hS (j := 0) (j' := m + 1) (by omega) (by omega) e; omega
  have hvw : t.sV i (m + 1) ≠ t.sV i (m + 1 + 1) := fun e => by
    have := t.sV_inj (g := g) hwf hsep hi hS (j := m + 1) (j' := m + 1 + 1) (by omega) (by omega) e
    omega
  have hsepAC := t.sep_chain_C (g := g) hwf hne hsep hi hS hm
  have hsepAW := t.sep_chain_W (g := g) hwf hne hsep hi hS hm
  have hsepWC := t.sep_W_C (g := g) hwf hne hsep hi hS (j := m + 1) (by omega) (by omega)
  have hmem1 : ∀ q, (t.chain g i (m + 1)).Mem q ↔ ((t.chain g i m).Mem q ∨
      (t.pieceBelow g (t.sW ((t.toSpqrTree.nvRange i).1 + (m + 1)))).Mem q ∨
      (t.pieceBelow g (t.sC i)[m + 1]!).Mem q) := by
    intro q
    show QE.edge q ∈ (t.chain g i (m + 1)).ves ↔ QE.edge q ∈ (t.chain g i m).ves ∨
      QE.edge q ∈ t.edgesBelow (t.sW ((t.toSpqrTree.nvRange i).1 + (m + 1))) ∨
      QE.edge q ∈ t.edgesBelow (t.sC i)[m + 1]!
    rw [chain_succ_ves, List.mem_append, List.mem_append]
  rcases t.sW_boundary hwf hne hsep hi hS (j := m + 1) (by omega) (by omega) s h.toGluedUpTo with
    ⟨hnone, hnil⟩ | ⟨w0, w1, τ, hw0', hw1', hopen⟩
  · -- no `V` piece at this corner: a single link `a2 ↔ b1`, then `a3 ↔ b0`
    split
    next hcond => rw [hnone] at hcond; exact absurd hcond (by simp)
    rw [t.nodeStep_link _ _ _ hr3 (by omega) (by omega) (by omega),
      t.treeQe_congr (link_outerE _ _ _), t.treeQe_congr (link_outerE _ _ _), hq3, hqb0]
    have hnoW : ∀ q, ¬ (t.pieceBelow g (t.sW ((t.toSpqrTree.nvRange i).1 + (m + 1)))).Mem q := by
      intro q hq; simp only [Piece.Mem, pieceBelow, hnil, List.not_mem_nil] at hq
    have hsz1 : ((link (some a2) (some b1)).run s').2.rotAdj.size = s'.rotAdj.size :=
      link_rotAdj_size _ _ _
    refine ⟨?_, ?_, ?_, a0, a1, b2, b3, ?_⟩
    · intro j hj; rw [link_outerE, link_outerE]; exact houter j hj
    · rw [link_rotAdj_size, link_rotAdj_size]; exact hrsz
    · intro q hq
      rw [hmem1] at hq
      push_neg at hq
      obtain ⟨hqA, -, hqC⟩ := hq
      rw [link_rotAdj_get_other _ _ _ (fun e => hqA (by rw [e]; exact hmemA3)) (fun e => hqC (by rw [e]; exact hmemB0)),
        link_rotAdj_get_other _ _ _ (fun e => hqA (by rw [e]; exact hmemA2)) (fun e => hqC (by rw [e]; exact hmemB1))]
      exact hframe q hqA
    have hE : ({t.chain g i m with ves := t.edgesBelow (t.sC i)[m + 1]!} : Piece).Capped s'.rotAdj σ
        b0 b1 b2 b3 (t.sV i (m + 1)) (t.sV i (m + 1 + 1)) :=
      hcapb.frame fun q hq => hframe q fun hq' => hAC q hq' hq
    obtain ⟨ρ', hρ'⟩ := Piece.Capped.join hcap hE huv hvw hdisAC hsepAC (A' :=
      ((link (some a3) (some b0)).run ((link (some a2) (some b1)).run s').2).2.rotAdj) (by
      intro q
      rw [link_rotAdj_get _ _ _ (by rw [hsz1]; exact hbdA a3 hmemA3) (by rw [hsz1]; exact hbdC b0 hmemB0)
        ha3b0, link_rotAdj_get _ _ _ (hbdA a2 hmemA2) (hbdC b1 hmemB1) ha2b1,
        ite_perm4 q a2 b1 a3 b0 _ _ _ _ _ ha2b1 ha3b0 ha23 ha2b0 ha3b1.symm hb01.symm])
    refine ⟨ρ', ha0, ha1, hb2, hb3, ?_⟩
    rw [chain_succ_eq, hnil, List.append_nil]
    exact hρ'
  · -- the `V` piece is attached by `a2 ↔ w1`, `b1 ↔ w0`, then the child by `a3 ↔ b0`
    split
    swap
    next hcond => exact absurd (by rw [hw0']; rfl) hcond
    rw [hw0', hw1']
    have hs2 : ((link (some b1) (some w0)).run ((link (some a2) (some w1)).run s').2).2.outerE =
        s'.outerE := by rw [link_outerE, link_outerE]
    rw [t.nodeStep_link _ _ _ hr3 (by omega) (by omega) (by omega),
      t.treeQe_congr hs2, t.treeQe_congr hs2, hq3, hqb0]
    have hmemW0 : (t.pieceBelow g (t.sW ((t.toSpqrTree.nvRange i).1 + (m + 1)))).Mem w0 := by
      obtain ⟨l0, l1, hl0, -, -, -⟩ := hopen.boundary; exact Piece.mem_of_loc hl0
    have hmemW1 : (t.pieceBelow g (t.sW ((t.toSpqrTree.nvRange i).1 + (m + 1)))).Mem w1 := by
      obtain ⟨l0, l1, -, hl1, -, -⟩ := hopen.boundary; exact Piece.mem_of_loc hl1
    have hw01 : w0 ≠ w1 := by have := hopen.left_dir; have := hopen.right_dir; omega
    have ha2w0 : a2 ≠ w0 := fun e => hAW a2 hmemA2 (by rw [e]; exact hmemW0)
    have ha2w1 : a2 ≠ w1 := fun e => hAW a2 hmemA2 (by rw [e]; exact hmemW1)
    have ha3w0 : a3 ≠ w0 := fun e => hAW a3 hmemA3 (by rw [e]; exact hmemW0)
    have ha3w1 : a3 ≠ w1 := fun e => hAW a3 hmemA3 (by rw [e]; exact hmemW1)
    have hw0b0 : w0 ≠ b0 := fun e => hWC w0 hmemW0 (by rw [e]; exact hmemB0)
    have hw0b1 : w0 ≠ b1 := fun e => hWC w0 hmemW0 (by rw [e]; exact hmemB1)
    have hw1b0 : w1 ≠ b0 := fun e => hWC w1 hmemW1 (by rw [e]; exact hmemB0)
    have hw1b1 : w1 ≠ b1 := fun e => hWC w1 hmemW1 (by rw [e]; exact hmemB1)
    have hsz1 : ((link (some a2) (some w1)).run s').2.rotAdj.size = s'.rotAdj.size :=
      link_rotAdj_size _ _ _
    have hsz2 : ((link (some b1) (some w0)).run ((link (some a2) (some w1)).run s').2).2.rotAdj.size =
        s'.rotAdj.size := by rw [link_rotAdj_size, hsz1]
    refine ⟨?_, ?_, ?_, a0, a1, b2, b3, ?_⟩
    · intro j hj; rw [link_outerE, link_outerE, link_outerE]; exact houter j hj
    · rw [link_rotAdj_size, link_rotAdj_size, link_rotAdj_size]; exact hrsz
    · intro q hq
      rw [hmem1] at hq
      push_neg at hq
      obtain ⟨hqA, hqW, hqC⟩ := hq
      rw [link_rotAdj_get_other _ _ _ (fun e => hqA (by rw [e]; exact hmemA3)) (fun e => hqC (by rw [e]; exact hmemB0)),
        link_rotAdj_get_other _ _ _ (fun e => hqC (by rw [e]; exact hmemB1)) (fun e => hqW (by rw [e]; exact hmemW0)),
        link_rotAdj_get_other _ _ _ (fun e => hqA (by rw [e]; exact hmemA2)) (fun e => hqW (by rw [e]; exact hmemW1))]
      exact hframe q hqA
    have hW : ({t.chain g i m with ves := t.edgesBelow (t.sW ((t.toSpqrTree.nvRange i).1 + (m + 1)))} :
        Piece).OpenEmbedding s'.rotAdj τ w0 w1 (t.sV i (m + 1)) :=
      hopen.frame fun q hq => hframe q fun hq' => hAW q hq' hq
    obtain ⟨ρ1, hρ1⟩ := Piece.Capped.attachOpen hcap hW huv hdisAW hsepAW
      (A' := ((link (some a2) (some w1)).run s').2.rotAdj)
      (fun q => link_rotAdj_get _ _ _ (hbdA a2 hmemA2) (hbdW w1 hmemW1) ha2w1 q)
    have hE : ({t.chain g i m with ves := t.edgesBelow (t.sC i)[m + 1]!} : Piece).Capped
        ((link (some a2) (some w1)).run s').2.rotAdj σ b0 b1 b2 b3 (t.sV i (m + 1)) (t.sV i (m + 1 + 1)) := by
      refine hcapb.frame fun q hq => ?_
      rw [link_rotAdj_get_other _ _ _ (fun e => hAC q (by rw [e]; exact hmemA2) hq) (fun e => hWC q (by rw [e]; exact hmemW1) hq)]
      exact hframe q fun hq' => hAC q hq' hq
    have hsepU : ∀ x,
        HasEdge ({t.chain g i m with ves := (t.chain g i m).ves ++ t.edgesBelow (t.sW ((t.toSpqrTree.nvRange i).1 + (m + 1)))} : Piece).es x →
        HasEdge ({t.chain g i m with ves := t.edgesBelow (t.sC i)[m + 1]!} : Piece).es x →
        x = t.sV i (m + 1) := by
      intro x hx hxc
      obtain ⟨p, hp, hpx⟩ := hx
      simp only [Piece.es, List.map_append, List.mem_append] at hp
      rcases hp with hp | hp
      · exact hsepAC x ⟨p, hp, hpx⟩ hxc
      · exact hsepWC x ⟨p, hp, hpx⟩ hxc
    obtain ⟨ρ', hρ'⟩ := Piece.Capped.join hρ1 hE huv hvw
      (List.disjoint_append_left.2 ⟨hdisAC, hdisWC⟩) hsepU (A' :=
      ((link (some a3) (some b0)).run ((link (some b1) (some w0)).run
        ((link (some a2) (some w1)).run s').2).2).2.rotAdj) (by
      intro q
      rw [link_rotAdj_get _ _ _ (by rw [hsz2]; exact hbdA a3 hmemA3) (by rw [hsz2]; exact hbdC b0 hmemB0)
        ha3b0, link_rotAdj_get _ _ _ (by rw [hsz1]; exact hbdC b1 hmemB1) (by rw [hsz1]; exact hbdW w0 hmemW0)
        hw0b1.symm, ite_perm4 q b1 w0 a3 b0 _ _ _ _ _ hw0b1.symm ha3b0 ha3b1.symm hb01.symm ha3w0.symm hw0b0,
        ite_swap2 q b1 w0 _ _ _ hw0b1.symm])
    refine ⟨ρ', ha0, ha1, hb2, hb3, ?_⟩
    rw [chain_succ_eq]
    exact hρ'

theorem invS_all (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape) (hne : t.ne = g.ne)
    (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size) (hS : t.toSpqrTree.type i = .S)
    (hlay : ∃ ev mr cv, t.neRotAdj.extract (4 * (t.toSpqrTree.neRange i).1) (4 * (t.toSpqrTree.neRange i).2) =
      layoutRot (t.toSpqrTree.type i) (t.toSpqrTree.nVerts i) (t.toSpqrTree.neRange i).1
        (t.toSpqrTree.neRange i).2 ev mr cv)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) (m : Nat) (hm : m ≤ t.toSpqrTree.nVerts i - 2) :
    t.InvS g i s (t.blockFold i (t.toSpqrTree.neRange i).1 s m) m := by
  induction m with
  | zero => exact t.invS_zero hwf hne hsep hi hS hlay s h
  | succ m ih => exact t.invS_step hwf hsh hne hsep hi hS hlay s h hm (ih (by omega))

/-- The last block (the last path edge) only meets smaller quarter-edges: a no-op. -/
theorem last_block_S (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hS : t.toSpqrTree.type i = .S)
    (hlay : ∃ ev mr cv, t.neRotAdj.extract (4 * (t.toSpqrTree.neRange i).1) (4 * (t.toSpqrTree.neRange i).2) =
      layoutRot (t.toSpqrTree.type i) (t.toSpqrTree.nVerts i) (t.toSpqrTree.neRange i).1
        (t.toSpqrTree.neRange i).2 ev mr cv)
    (s : EmbedState) :
    t.blockFold i (t.toSpqrTree.neRange i).1 s (t.toSpqrTree.nVerts i - 2 + 1) =
      t.blockFold i (t.toSpqrTree.neRange i).1 s (t.toSpqrTree.nVerts i - 2) := by
  have h3 := (t.shape_S hwf hi hS).1
  rw [blockFold_succ, range'_four]
  simp only [List.foldl_cons, List.foldl_nil]
  have hk : t.toSpqrTree.nVerts i - 2 + 1 < t.toSpqrTree.nVerts i := by omega
  have hr := fun r (hr : r < 4) => t.rotS_S hwf hi hS hlay (k := t.toSpqrTree.nVerts i - 2 + 1) (r := r) hk hr
    (ta := 4 * ((t.toSpqrTree.neRange i).1 + (t.toSpqrTree.nVerts i - 2 + 1)) + r) rfl rfl
  have hr0 := t.rotS_S hwf hi hS hlay (k := t.toSpqrTree.nVerts i - 2 + 1) (r := 0) hk (by omega)
    (ta := 4 * ((t.toSpqrTree.neRange i).1 + (t.toSpqrTree.nVerts i - 2 + 1))) (by omega) rfl
  have hne1 : t.toSpqrTree.nVerts i - 2 + 1 ≠ 1 := by omega
  have hlt : ∀ r, r < 4 → rotS (t.toSpqrTree.nVerts i) (t.toSpqrTree.neRange i).1
      (t.toSpqrTree.nVerts i - 2 + 1) r < 4 * ((t.toSpqrTree.neRange i).1 + (t.toSpqrTree.nVerts i - 2 + 1)) + r := by
    intro r hr
    simp only [rotS, Nat.succ_ne_zero, ite_false, hne1, show t.toSpqrTree.nVerts i - 2 + 1 + 1 = t.toSpqrTree.nVerts i by omega,
      ite_true]
    split_ifs <;> omega
  rw [t.nodeStep_skip _ _ _ hr0 (by have := hlt 0 (by omega); omega),
    t.nodeStep_skip _ _ _ (hr 1 (by omega)) (by have := hlt 1 (by omega); omega),
    t.nodeStep_skip _ _ _ (hr 2 (by omega)) (by have := hlt 2 (by omega); omega),
    t.nodeStep_skip _ _ _ (hr 3 (by omega)) (by have := hlt 3 (by omega); omega)]

end S

end Spqr.PlanarSpqrTree
