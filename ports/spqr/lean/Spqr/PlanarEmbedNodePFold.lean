import Spqr.PlanarEmbedNodeP
import Spqr.Proofs.PieceParJoin

/-!
# The `P` node fold: `InvP` and `nodeFold_capped_P`

`pChain g i m` is the piece of the first `m + 1` children of the `P` node `i` (in bond order);
`InvP` is the invariant after blocks `0..m` of the node loop: the chain is `Capped` at the cap's
original endpoints `sV 0`, `sV 1` with pairs `(c0, c3)` on the first child and `(c1, c2)` on child
`m`. Block `m + 1` is `Capped.parJoin`; the last block is a no-op; the full chain is
`pieceBelow g i` itself (no `V` children).
-/

namespace Spqr.PlanarSpqrTree

variable (t : PlanarSpqrTree)

open EmbedM

section P

variable {i : Nat}

/-- The first `m + 1` children of a `P` node, in bond order. -/
def pL (i m : Nat) : List Nat := (List.range (m + 1)).map fun j => (t.sC i)[j]!

/-- The piece of the first `m + 1` children. -/
def pChain (g : Graph) (i m : Nat) : Piece :=
  {t.pieceBelow g i with ves := (t.pL i m).flatMap t.edgesBelow}

theorem pL_succ (i m : Nat) : t.pL i (m + 1) = t.pL i m ++ [(t.sC i)[m + 1]!] := by
  unfold pL; rw [List.range_succ, List.map_append]; rfl

theorem pChain_zero_eq (g : Graph) : t.pChain g i 0 = t.pieceBelow g (t.sC i)[0]! := by
  simp [pChain, pL, pieceBelow]

theorem pChain_succ_eq (g : Graph) (m : Nat) :
    t.pChain g i (m + 1) =
      {t.pChain g i m with ves := (t.pChain g i m).ves ++ t.edgesBelow (t.sC i)[m + 1]!} := by
  simp [pChain, pL_succ, List.flatMap_append]

theorem mem_pL {m c' : Nat} (h : c' ∈ t.pL i m) : ∃ j, j ≤ m ∧ (t.sC i)[j]! = c' := by
  unfold pL at h
  simp only [List.mem_map, List.mem_range] at h
  obtain ⟨j, hj, rfl⟩ := h
  exact ⟨j, by omega, rfl⟩

theorem mem_pChain_iff (g : Graph) (m q : Nat) :
    (t.pChain g i m).Mem q ↔ ∃ c' ∈ t.pL i m, (t.pieceBelow g c').Mem q := by
  simp only [Piece.Mem, pChain, pieceBelow, List.mem_flatMap]

theorem mem_pChain_bound (hwf : t.toSpqrTree.WF) (g : Graph) {m q : Nat}
    (hq : (t.pChain g i m).Mem q) : q < 4 * t.ne := by
  obtain ⟨c', -, hq⟩ := (t.mem_pChain_iff g m q).1 hq
  exact t.mem_pieceBelow_bound hwf g hq

theorem sC_ne_P (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hP : t.toSpqrTree.type i = .P)
    {j j' : Nat} (hj : j < t.toSpqrTree.nEdges i - 1) (hj' : j' < t.toSpqrTree.nEdges i - 1)
    (hne : j ≠ j') : (t.sC i)[j]! ≠ (t.sC i)[j']! := by
  intro e
  have h1 := t.sC_getElem?_P hwf hi hP hj
  have h2 := t.sC_getElem?_P hwf hi hP hj'
  rw [e, ← h2] at h1
  have hl := t.sC_length_P hwf hi hP
  exact hne ((List.Nodup.getElem?_inj (by omega) (t.sC_nodup hwf hi)).1 h1)

theorem disjoint_pChain_C (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hi : i < t.size) (hP : t.toSpqrTree.type i = .P) {m : Nat}
    (hm : m + 1 < t.toSpqrTree.nEdges i - 1) :
    List.Disjoint (t.pChain g i m).ves (t.edgesBelow (t.sC i)[m + 1]!) := by
  rw [List.disjoint_left]
  intro e he he'
  simp only [pChain, List.mem_flatMap] at he
  obtain ⟨c', hc', he⟩ := he
  obtain ⟨j, hj, rfl⟩ := t.mem_pL hc'
  have hcj := List.mem_of_getElem? (t.sC_getElem?_P hwf hi hP (j := j) (by omega))
  have hcm := List.mem_of_getElem? (t.sC_getElem?_P hwf hi hP (j := m + 1) hm)
  rw [sC, List.mem_filter] at hcj hcm
  exact List.disjoint_left.1 (t.maximal_pieces_disjoint hwf hsh (t.child_maximal hwf hi hcj.1)
    (t.child_maximal hwf hi hcm.1) (t.sC_ne_P hwf hi hP (by omega) hm (by omega))) he he'

theorem sep_pChain_C (hwf : t.toSpqrTree.WF) (hne : t.ne = g.ne) (hsep : t.toSpqrTree.PieceSep g)
    (hi : i < t.size) (hP : t.toSpqrTree.type i = .P) {m : Nat}
    (hm : m + 1 < t.toSpqrTree.nEdges i - 1) :
    ∀ x, HasEdge (t.pChain g i m).es x → HasEdge (t.pieceBelow g (t.sC i)[m + 1]!).es x →
      x = t.sV i 0 ∨ x = t.sV i 1 := by
  intro x hx hxc
  obtain ⟨p, hp, hpx⟩ := hx
  simp only [Piece.es, pChain, List.mem_map, List.mem_flatMap] at hp
  obtain ⟨e, ⟨c', hc', he⟩, rfl⟩ := hp
  obtain ⟨j, hj, rfl⟩ := t.mem_pL hc'
  have hcj := List.mem_of_getElem? (t.sC_getElem?_P hwf hi hP (j := j) (by omega))
  have hcm := List.mem_of_getElem? (t.sC_getElem?_P hwf hi hP (j := m + 1) hm)
  rw [sC, List.mem_filter] at hcj hcm
  have hta : t.toSpqrTree.Touches g (t.sC i)[j]! x :=
    t.touches_of_hasEdge hwf g hne ⟨_, List.mem_map.2 ⟨e, he, rfl⟩, hpx⟩
  have htb : t.toSpqrTree.Touches g (t.sC i)[m + 1]! x := t.touches_of_hasEdge hwf g hne hxc
  exact t.attach_P (g := g) hwf hsep hi hP hcj.1 hcm.1 (t.sC_ne_P hwf hi hP (by omega) hm (by omega))
    hta htb

/-- Invariant after blocks `0..m` of a `P` node's loop. -/
structure InvP (g : Graph) (i : Nat) (s s' : EmbedState) (m : Nat) : Prop where
  outer_ne : ∀ j, j ≠ i → s'.outerE[j]? = s.outerE[j]?
  rot_size : s'.rotAdj.size = s.rotAdj.size
  frame : ∀ q, ¬ (t.pChain g i m).Mem q → s'.rotAdj[q]? = s.rotAdj[q]?
  capped : ∃ a0 a1 a2 a3 ρ, s.outerE[(t.sC i)[0]!]![0]! = some a0 ∧
    s.outerE[(t.sC i)[m]!]![1]! = some a1 ∧
    s.outerE[(t.sC i)[m]!]![2]! = some a2 ∧ s.outerE[(t.sC i)[0]!]![3]! = some a3 ∧
    (t.pChain g i m).Capped s'.rotAdj ρ a0 a1 a2 a3 (t.sV i 0) (t.sV i 1)

/-- Block `0` (the cap edge): four `setOuter`s exposing slots `0, 3` of the first child and `1, 2`
of the last child; `rotAdj` untouched. -/
theorem cap_block_P (hwf : t.toSpqrTree.WF) (hne : t.ne = g.ne) (hsep : t.toSpqrTree.PieceSep g)
    (hi : i < t.size) (hP : t.toSpqrTree.type i = .P) (hlay : t.LayoutAt i)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    (t.blockFold i (t.toSpqrTree.neRange i).1 s 0).rotAdj = s.rotAdj ∧
    ∃ a0 e1 e2 a3, s.outerE[(t.sC i)[0]!]![0]! = some a0 ∧
      s.outerE[(t.sC i)[t.toSpqrTree.nEdges i - 2]!]![1]! = some e1 ∧
      s.outerE[(t.sC i)[t.toSpqrTree.nEdges i - 2]!]![2]! = some e2 ∧
      s.outerE[(t.sC i)[0]!]![3]! = some a3 ∧
      (t.blockFold i (t.toSpqrTree.neRange i).1 s 0).outerE[i]? =
        some #[some a0, some e1, some e2, some a3] := by
  have h3 := (t.shape_P hwf hi hP).2.1
  obtain ⟨hc0, a0, -, -, a3, -, ha0, -, -, ha3, -⟩ :=
    t.child_cert_P hwf hne hsep hi hP (j := 0) (by omega) s h
  obtain ⟨hcl, -, e1, e2, -, -, -, he1, he2, -, -⟩ :=
    t.child_cert_P hwf hne hsep hi hP (j := t.toSpqrTree.nEdges i - 2) (by omega) s h
  have hc0i : (t.sC i)[0]! ≠ i := fun e => by have := (t.child_data hwf hi hc0).2.2; omega
  have hcli : (t.sC i)[t.toSpqrTree.nEdges i - 2]! ≠ i := fun e => by
    have := (t.child_data hwf hi hcl).2.2; omega
  obtain ⟨-, -, -, -, -, -, -, hq0⟩ :=
    t.child_P (g := g) hwf hsep hi hP (j := 0) (by omega) (t.sC_getElem?_P hwf hi hP (by omega))
  obtain ⟨-, -, -, -, -, -, -, hql⟩ :=
    t.child_P (g := g) hwf hsep hi hP (j := t.toSpqrTree.nEdges i - 2) (by omega)
      (t.sC_getElem?_P hwf hi hP (by omega))
  have hr0 := t.rotP_P hwf hi hP hlay (k := 0) (r := 0) (by omega) (by omega)
    (ta := 4 * (t.toSpqrTree.neRange i).1)
    (tb := 4 * ((t.toSpqrTree.neRange i).1 + (t.toSpqrTree.nEdges i - 2 + 1)) + 1)
    (by omega) (by simp [rotP]; omega)
  have hr1 := t.rotP_P hwf hi hP hlay (k := 0) (r := 1) (by omega) (by omega)
    (ta := 4 * (t.toSpqrTree.neRange i).1 + 1) (tb := 4 * ((t.toSpqrTree.neRange i).1 + (0 + 1)) + 0)
    (by omega) (by simp [rotP]; omega)
  have hr2 := t.rotP_P hwf hi hP hlay (k := 0) (r := 2) (by omega) (by omega)
    (ta := 4 * (t.toSpqrTree.neRange i).1 + 2) (tb := 4 * ((t.toSpqrTree.neRange i).1 + (0 + 1)) + 3)
    (by omega) (by simp [rotP]; omega)
  have hr3 := t.rotP_P hwf hi hP hlay (k := 0) (r := 3) (by omega) (by omega)
    (ta := 4 * (t.toSpqrTree.neRange i).1 + 3)
    (tb := 4 * ((t.toSpqrTree.neRange i).1 + (t.toSpqrTree.nEdges i - 2 + 1)) + 2)
    (by omega) (by simp [rotP]; omega)
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
  rw [t.nodeStep_cap i _ s hr0 (by omega) (by omega), hsl0, hql 1 (by omega)]
  rw [t.nodeStep_cap i _ _ hr1 (by omega) (by omega), hsl1, hq0 0 (by omega), setOuter_row_ne _ _ _ _ hc0i]
  rw [t.nodeStep_cap i _ _ hr2 (by omega) (by omega), hsl2, hq0 3 (by omega),
    setOuter_row_ne _ _ _ _ hc0i, setOuter_row_ne _ _ _ _ hc0i]
  rw [t.nodeStep_cap i _ _ hr3 (by omega) (by omega), hsl3, hql 2 (by omega),
    setOuter_row_ne _ _ _ _ hcli, setOuter_row_ne _ _ _ _ hcli, setOuter_row_ne _ _ _ _ hcli]
  rw [ha0, ha3, he1, he2]
  refine ⟨by simp only [setOuter_rotAdj], a0, e1, e2, a3, rfl, rfl, rfl, rfl, ?_⟩
  rw [setOuter_outerE_eq' _ _ _ _ (setOuter_outerE_eq' _ _ _ _ (setOuter_outerE_eq' _ _ _ _
    (setOuter_outerE_eq' _ _ _ _ hrow)))]
  rw [array_four r hr4]
  rfl

theorem invP_zero (hwf : t.toSpqrTree.WF) (hne : t.ne = g.ne) (hsep : t.toSpqrTree.PieceSep g)
    (hi : i < t.size) (hP : t.toSpqrTree.type i = .P) (hlay : t.LayoutAt i)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    t.InvP g i s (t.blockFold i (t.toSpqrTree.neRange i).1 s 0) 0 := by
  have h3 := (t.shape_P hwf hi hP).2.1
  obtain ⟨hrot, -⟩ := t.cap_block_P hwf hne hsep hi hP hlay s h
  obtain ⟨-, b0, b1, b2, b3, ρ, hb0, hb1, hb2, hb3, hcap⟩ :=
    t.child_cert_P hwf hne hsep hi hP (j := 0) (by omega) s h
  refine ⟨fun j hj => t.blockFold_outerE_ne _ _ _ _ hj, by rw [hrot], fun q _ => by rw [hrot],
    b0, b1, b2, b3, ρ, hb0, hb1, hb2, hb3, ?_⟩
  rw [hrot, pChain_zero_eq]
  exact hcap

set_option maxHeartbeats 1000000 in
/-- Block `m + 1` (`m + 1 ≤ k - 2`): links `a1 ↔ b0`, `a2 ↔ b3` between child `m` and child
`m + 1` = `Capped.parJoin`. -/
theorem invP_step (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape) (hne : t.ne = g.ne)
    (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size) (hP : t.toSpqrTree.type i = .P)
    (hlay : t.LayoutAt i)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) {m : Nat} (hm : m + 1 ≤ t.toSpqrTree.nEdges i - 2)
    (hinv : t.InvP g i s (t.blockFold i (t.toSpqrTree.neRange i).1 s m) m) :
    t.InvP g i s (t.blockFold i (t.toSpqrTree.neRange i).1 s (m + 1)) (m + 1) := by
  have h3 := (t.shape_P hwf hi hP).2.1
  rw [blockFold_succ, range'_four]
  simp only [List.foldl_cons, List.foldl_nil]
  generalize t.blockFold i (t.toSpqrTree.neRange i).1 s m = s' at hinv
  obtain ⟨houter, hrsz, hframe, a0, a1, a2, a3, ρ, ha0, ha1, ha2, ha3, hcap⟩ := hinv
  have hcm := t.sC_getElem?_P hwf hi hP (j := m) (by omega)
  have hcm1 := t.sC_getElem?_P hwf hi hP (j := m + 1) (by omega)
  obtain ⟨hcm_mem, -, -, -, -, -, -, hqm⟩ := t.child_P (g := g) hwf hsep hi hP (j := m) (by omega) hcm
  obtain ⟨hcm1_mem, -, -, -, -, -, -, hqm1⟩ :=
    t.child_P (g := g) hwf hsep hi hP (j := m + 1) (by omega) hcm1
  obtain ⟨-, b0, b1, b2, b3, σ, hb0, hb1, hb2, hb3, hcapb⟩ :=
    t.child_cert_P hwf hne hsep hi hP (j := m + 1) (by omega) s h
  have hcmi : (t.sC i)[m]! ≠ i := fun e => by have := (t.child_data hwf hi hcm_mem).2.2; omega
  have hcm1i : (t.sC i)[m + 1]! ≠ i := fun e => by have := (t.child_data hwf hi hcm1_mem).2.2; omega
  have hr0 := t.rotP_P hwf hi hP hlay (k := m + 1) (r := 0) (by omega) (by omega)
    (ta := 4 * ((t.toSpqrTree.neRange i).1 + (m + 1))) rfl rfl
  have hr3 := t.rotP_P hwf hi hP hlay (k := m + 1) (r := 3) (by omega) (by omega)
    (ta := 4 * ((t.toSpqrTree.neRange i).1 + (m + 1)) + 3) rfl rfl
  have hr1 := t.rotP_P hwf hi hP hlay (k := m + 1) (r := 1) (by omega) (by omega)
    (ta := 4 * ((t.toSpqrTree.neRange i).1 + (m + 1)) + 1)
    (tb := 4 * ((t.toSpqrTree.neRange i).1 + (m + 1 + 1)) + 0) rfl
    (by simp [rotP, show ¬ (m + 1 + 1 = t.toSpqrTree.nEdges i) by omega])
  have hr2 := t.rotP_P hwf hi hP hlay (k := m + 1) (r := 2) (by omega) (by omega)
    (ta := 4 * ((t.toSpqrTree.neRange i).1 + (m + 1)) + 2)
    (tb := 4 * ((t.toSpqrTree.neRange i).1 + (m + 1 + 1)) + 3) rfl
    (by simp [rotP, show ¬ (m + 1 + 1 = t.toSpqrTree.nEdges i) by omega])
  have hlt0 : rotP (t.toSpqrTree.nEdges i) (t.toSpqrTree.neRange i).1 (m + 1) 0 <
      4 * ((t.toSpqrTree.neRange i).1 + (m + 1)) := by
    simp only [rotP, Nat.succ_ne_zero, ite_false]; split_ifs <;> omega
  have hlt3 : rotP (t.toSpqrTree.nEdges i) (t.toSpqrTree.neRange i).1 (m + 1) 3 <
      4 * ((t.toSpqrTree.neRange i).1 + (m + 1)) + 3 := by
    simp only [rotP, Nat.succ_ne_zero, ite_false]; split_ifs <;> omega
  rw [t.nodeStep_skip _ _ _ hr0 hlt0,
    t.nodeStep_link _ _ _ hr1 (by omega) (by omega) (by omega)]
  have hq1 : t.treeQe s' (4 * ((t.toSpqrTree.neRange i).1 + (m + 1)) + 1) = some a1 := by
    rw [hqm 1 (by omega), getElem!_eq_of_getElem?_eq (houter _ hcmi), ha1]
  have hq2 : t.treeQe s' (4 * ((t.toSpqrTree.neRange i).1 + (m + 1)) + 2) = some a2 := by
    rw [hqm 2 (by omega), getElem!_eq_of_getElem?_eq (houter _ hcmi), ha2]
  have hqb0 : t.treeQe s' (4 * ((t.toSpqrTree.neRange i).1 + (m + 1 + 1)) + 0) = some b0 := by
    rw [hqm1 0 (by omega), getElem!_eq_of_getElem?_eq (houter _ hcm1i), hb0]
  have hqb3 : t.treeQe s' (4 * ((t.toSpqrTree.neRange i).1 + (m + 1 + 1)) + 3) = some b3 := by
    rw [hqm1 3 (by omega), getElem!_eq_of_getElem?_eq (houter _ hcm1i), hb3]
  rw [hq1, hqb0, t.nodeStep_link _ _ _ hr2 (by omega) (by omega) (by omega),
    t.treeQe_congr (link_outerE _ _ _), t.treeQe_congr (link_outerE _ _ _), hq2, hqb3,
    t.nodeStep_skip _ _ _ hr3 hlt3]
  -- membership, bounds and distinctness of the linked slots
  have hmemA1 : (t.pChain g i m).Mem a1 := by
    obtain ⟨l0, l1, -, hl1, -, -⟩ := hcap.pair0; exact Piece.mem_of_loc hl1
  have hmemA2 : (t.pChain g i m).Mem a2 := by
    obtain ⟨l2, l3, hl2, -, -, -⟩ := hcap.pair2; exact Piece.mem_of_loc hl2
  have hmemB0 : (t.pieceBelow g (t.sC i)[m + 1]!).Mem b0 := by
    obtain ⟨l0, l1, hl0, -, -, -⟩ := hcapb.pair0; exact Piece.mem_of_loc hl0
  have hmemB3 : (t.pieceBelow g (t.sC i)[m + 1]!).Mem b3 := by
    obtain ⟨l2, l3, -, hl3, -, -⟩ := hcapb.pair2; exact Piece.mem_of_loc hl3
  have hbdA : ∀ q, (t.pChain g i m).Mem q → q < s'.rotAdj.size := by
    intro q hq; rw [hrsz, h.rot_size]; exact t.mem_pChain_bound hwf g hq
  have hbdC : ∀ q, (t.pieceBelow g (t.sC i)[m + 1]!).Mem q → q < s'.rotAdj.size := by
    intro q hq; rw [hrsz, h.rot_size]; exact t.mem_pieceBelow_bound hwf g hq
  have hdisAC := t.disjoint_pChain_C (g := g) hwf hsh hi hP (m := m) (by omega)
  have hAC : ∀ q, (t.pChain g i m).Mem q → (t.pieceBelow g (t.sC i)[m + 1]!).Mem q → False :=
    fun q h1 h2 => List.disjoint_left.1 hdisAC h1 h2
  have ha12 : a1 ≠ a2 := by have := hcap.dir1; have := hcap.dir2; omega
  have hb03 : b0 ≠ b3 := by have := hcapb.dir0; have := hcapb.dir3; omega
  have ha1b0 : a1 ≠ b0 := fun e => hAC a1 hmemA1 (by rw [e]; exact hmemB0)
  have ha1b3 : a1 ≠ b3 := fun e => hAC a1 hmemA1 (by rw [e]; exact hmemB3)
  have ha2b0 : a2 ≠ b0 := fun e => hAC a2 hmemA2 (by rw [e]; exact hmemB0)
  have ha2b3 : a2 ≠ b3 := fun e => hAC a2 hmemA2 (by rw [e]; exact hmemB3)
  have huv := t.sV01_ne_P (g := g) hwf hsep hi hP
  have hsepAC := t.sep_pChain_C (g := g) hwf hne hsep hi hP (m := m) (by omega)
  have hmem1 : ∀ q, (t.pChain g i (m + 1)).Mem q ↔
      ((t.pChain g i m).Mem q ∨ (t.pieceBelow g (t.sC i)[m + 1]!).Mem q) := by
    intro q
    show QE.edge q ∈ (t.pChain g i (m + 1)).ves ↔ QE.edge q ∈ (t.pChain g i m).ves ∨
      QE.edge q ∈ t.edgesBelow (t.sC i)[m + 1]!
    rw [pChain_succ_eq, List.mem_append]
  have hsz1 : ((link (some a1) (some b0)).run s').2.rotAdj.size = s'.rotAdj.size :=
    link_rotAdj_size _ _ _
  refine ⟨?_, ?_, ?_, a0, b1, b2, a3, ?_⟩
  · intro j hj; rw [link_outerE, link_outerE]; exact houter j hj
  · rw [link_rotAdj_size, link_rotAdj_size]; exact hrsz
  · intro q hq
    rw [hmem1, not_or] at hq
    obtain ⟨hqA, hqC⟩ := hq
    rw [link_rotAdj_get_other _ _ _ (fun e => hqA (by rw [e]; exact hmemA2)) (fun e => hqC (by rw [e]; exact hmemB3)),
      link_rotAdj_get_other _ _ _ (fun e => hqA (by rw [e]; exact hmemA1)) (fun e => hqC (by rw [e]; exact hmemB0))]
    exact hframe q hqA
  have hE : ({t.pChain g i m with ves := t.edgesBelow (t.sC i)[m + 1]!} : Piece).Capped s'.rotAdj σ
      b0 b1 b2 b3 (t.sV i 0) (t.sV i 1) :=
    hcapb.frame fun q hq => hframe q fun hq' => hAC q hq' hq
  obtain ⟨ρ', hρ'⟩ := Piece.Capped.parJoin hcap hE huv hdisAC hsepAC (A' :=
    ((link (some a2) (some b3)).run ((link (some a1) (some b0)).run s').2).2.rotAdj) (by
    intro q
    rw [link_rotAdj_get _ _ _ (by rw [hsz1]; exact hbdA a2 hmemA2) (by rw [hsz1]; exact hbdC b3 hmemB3)
      ha2b3, link_rotAdj_get _ _ _ (hbdA a1 hmemA1) (hbdC b0 hmemB0) ha1b0]
    split_ifs <;> first | rfl | omega)
  refine ⟨ρ', ha0, hb1, hb2, ha3, ?_⟩
  rw [pChain_succ_eq]
  exact hρ'

theorem invP_all (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape) (hne : t.ne = g.ne)
    (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size) (hP : t.toSpqrTree.type i = .P)
    (hlay : t.LayoutAt i) (s : EmbedState) (h : t.GluedFaces g (i + 1) s) (m : Nat)
    (hm : m ≤ t.toSpqrTree.nEdges i - 2) :
    t.InvP g i s (t.blockFold i (t.toSpqrTree.neRange i).1 s m) m := by
  induction m with
  | zero => exact t.invP_zero hwf hne hsep hi hP hlay s h
  | succ m ih => exact t.invP_step hwf hsh hne hsep hi hP hlay s h hm (ih (by omega))

/-- The last block only meets smaller quarter-edges: a no-op. -/
theorem last_block_P (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hP : t.toSpqrTree.type i = .P)
    (hlay : t.LayoutAt i) (s : EmbedState) :
    t.blockFold i (t.toSpqrTree.neRange i).1 s (t.toSpqrTree.nEdges i - 2 + 1) =
      t.blockFold i (t.toSpqrTree.neRange i).1 s (t.toSpqrTree.nEdges i - 2) := by
  have h3 := (t.shape_P hwf hi hP).2.1
  rw [blockFold_succ, range'_four]
  simp only [List.foldl_cons, List.foldl_nil]
  have hk : t.toSpqrTree.nEdges i - 2 + 1 < t.toSpqrTree.nEdges i := by omega
  have hr := fun r (hr : r < 4) => t.rotP_P hwf hi hP hlay (k := t.toSpqrTree.nEdges i - 2 + 1) (r := r) hk hr
    (ta := 4 * ((t.toSpqrTree.neRange i).1 + (t.toSpqrTree.nEdges i - 2 + 1)) + r) rfl rfl
  have hr0 := t.rotP_P hwf hi hP hlay (k := t.toSpqrTree.nEdges i - 2 + 1) (r := 0) hk (by omega)
    (ta := 4 * ((t.toSpqrTree.neRange i).1 + (t.toSpqrTree.nEdges i - 2 + 1))) (by omega) rfl
  have hlt : ∀ r, r < 4 → rotP (t.toSpqrTree.nEdges i) (t.toSpqrTree.neRange i).1
      (t.toSpqrTree.nEdges i - 2 + 1) r <
        4 * ((t.toSpqrTree.neRange i).1 + (t.toSpqrTree.nEdges i - 2 + 1)) + r := by
    intro r hr
    simp only [rotP, Nat.succ_ne_zero, ite_false,
      show t.toSpqrTree.nEdges i - 2 + 1 + 1 = t.toSpqrTree.nEdges i by omega, ite_true]
    split_ifs <;> omega
  rw [t.nodeStep_skip _ _ _ hr0 (by have := hlt 0 (by omega); omega),
    t.nodeStep_skip _ _ _ (hr 1 (by omega)) (by have := hlt 1 (by omega); omega),
    t.nodeStep_skip _ _ _ (hr 2 (by omega)) (by have := hlt 2 (by omega); omega),
    t.nodeStep_skip _ _ _ (hr 3 (by omega)) (by have := hlt 3 (by omega); omega)]

theorem fold_eq_blockFold_P (hwf : t.toSpqrTree.WF) (hi : i < t.size) (hP : t.toSpqrTree.type i = .P)
    (s : EmbedState) :
    (List.range' (4 * t.neBounds[i]!) (4 * t.neBounds[i + 1]! - 4 * t.neBounds[i]!)).foldl
      (t.nodeStep i t.neBounds[i]!) s =
      t.blockFold i (t.toSpqrTree.neRange i).1 s (t.toSpqrTree.nEdges i - 2 + 1) := by
  have h3 := (t.shape_P hwf hi hP).2.1
  unfold SpqrTree.nEdges at h3 ⊢
  rw [PlanarRot.neRange_eq] at h3 ⊢
  unfold blockFold
  dsimp only at h3 ⊢
  rw [show 4 * t.neBounds[i + 1]! - 4 * t.neBounds[i]! =
    4 * (t.neBounds[i + 1]! - t.neBounds[i]! - 2 + 1 + 1) by omega]

/-- The chain of all children of a `P` node is its piece. -/
theorem pChain_eq_pieceBelow (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size) (hP : t.toSpqrTree.type i = .P) :
    t.pChain g i (t.toSpqrTree.nEdges i - 2) = t.pieceBelow g i := by
  have h3 := (t.shape_P hwf hi hP).2.1
  have hl := t.sC_length_P hwf hi hP
  have hpL : t.pL i (t.toSpqrTree.nEdges i - 2) = t.sC i := by
    apply List.ext_getElem?
    intro j
    unfold pL
    by_cases hj : j < t.toSpqrTree.nEdges i - 2 + 1
    · rw [List.getElem?_map, List.getElem?_eq_getElem (by simpa using hj), List.getElem_range,
        Option.map_some, t.sC_getElem?_P hwf hi hP (by omega)]
    · rw [List.getElem?_map, List.getElem?_eq_none (by simp; omega), List.getElem?_eq_none (by omega)]
      rfl
  unfold pChain
  rw [hpL, t.sC_eq_children_P (g := g) hwf hsep hi hP]
  have he : t.edgesBelow i = (t.children i).flatMap t.edgesBelow := by
    rw [t.edgesBelow_eq hwf hsh i hi, ← t.type_eq_of_lt i hi, hP]; rfl
  unfold pieceBelow
  rw [he]

/-- `nodeFold_capped` for `P` nodes: iterated `Capped.parJoin` along the bond layout. -/
theorem nodeFold_capped_P (g : Graph) (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hrep : t.toSpqrTree.Represents g) (hsep : t.toSpqrTree.PieceSep g) (hi : i < t.size)
    (hP : t.toSpqrTree.type i = .P) (hlay : t.LayoutAt i)
    {ne : Nat} (hcne : t.toSpqrTree.capNe i = some ne) {p : Nat × Nat}
    (hp : t.toSpqrTree.neOrig ne = some p)
    (s : EmbedState) (h : t.GluedFaces g (i + 1) s) :
    let s' := (List.range' (4 * t.neBounds[i]!) (4 * t.neBounds[i + 1]! - 4 * t.neBounds[i]!)).foldl
      (t.nodeStep i t.neBounds[i]!) s
    (∀ q, ¬(t.pieceBelow g i).Mem q → s'.rotAdj[q]? = s.rotAdj[q]?) ∧
    ∃ a0 a1 a2 a3 ρ, s'.outerE[i]? = some #[some a0, some a1, some a2, some a3] ∧
      (t.pieceBelow g i).Capped s'.rotAdj ρ a0 a1 a2 a3 p.1 p.2 := by
  dsimp only
  have hne := hrep.ne
  have h3 := (t.shape_P hwf hi hP).2.1
  rw [t.fold_eq_blockFold_P hwf hi hP, t.last_block_P hwf hi hP hlay]
  obtain ⟨houter, hrsz, hframe, a0, a1, a2, a3, ρ, ha0, ha1, ha2, ha3, hcap⟩ :=
    t.invP_all hwf hsh hne hsep hi hP hlay s h (t.toSpqrTree.nEdges i - 2) le_rfl
  obtain ⟨-, b0, e1, e2, b3, hb0, he1, he2, hb3, hrow⟩ := t.cap_block_P hwf hne hsep hi hP hlay s h
  rw [ha0] at hb0; rw [ha1] at he1; rw [ha2] at he2; rw [ha3] at hb3
  cases hb0; cases he1; cases he2; cases hb3
  rw [t.capNe_P hP] at hcne
  cases hcne
  have hp' := t.capOrig_P (g := g) hwf hsep hi hP hp
  subst hp'
  rw [t.pChain_eq_pieceBelow (g := g) hwf hsh hsep hi hP] at hframe hcap
  refine ⟨hframe, a0, a1, a2, a3, ρ, ?_, hcap⟩
  rw [t.blockFold_outerE, hrow]

end P

end Spqr.PlanarSpqrTree
