import Mathlib.Tactic.IntervalCases
import Spqr.Proofs.Planar
import Spqr.PlanarRelabel

/-!
# The S and P local rotation systems are planar embeddings

Closed forms of `layoutRot .S` (cycle) and `layoutRot .P` (bond), and the proof that each is an
`IsPlanarEmbedding` of the corresponding skeleton: the vertex / face steps are explicit, the
vertex step of the cycle and the face step of the bond are fixed-point-free involutions, and the
other step has exactly four orbits (one per direction and side of the cap).
-/

namespace Spqr

/-! ### Folds appending fixed-size blocks -/

theorem foldl_append_blocks_size {α} (f : Nat → Array α) (hf : ∀ i, (f i).size = 4) (init : Array α) :
    ∀ m, ((List.range m).foldl (fun acc i => acc ++ f i) init).size = init.size + 4 * m := by
  intro m
  induction m with
  | zero => simp
  | succ m ih => rw [List.range_succ, List.foldl_append, List.foldl_cons, List.foldl_nil, Array.size_append, ih, hf]; omega

theorem foldl_append_blocks_get {α} (f : Nat → Array α) (hf : ∀ i, (f i).size = 4) (init : Array α) :
    ∀ m q, ((List.range m).foldl (fun acc i => acc ++ f i) init)[q]? =
      if q < init.size then init[q]?
      else if q < init.size + 4 * m then (f ((q - init.size) / 4))[(q - init.size) % 4]? else none := by
  intro m
  induction m with
  | zero =>
    intro q
    simp only [List.range_zero, List.foldl_nil, Nat.mul_zero, Nat.add_zero]
    by_cases h : q < init.size
    · simp only [h, ↓reduceIte]
    · simp only [h, ↓reduceIte]
      exact Array.getElem?_eq_none (by omega)
  | succ m ih =>
    intro q
    rw [List.range_succ, List.foldl_append, List.foldl_cons, List.foldl_nil, Array.getElem?_append,
      foldl_append_blocks_size f hf, ih]
    by_cases h1 : q < init.size
    · simp only [h1, show q < init.size + 4 * m by omega, ↓reduceIte]
    · by_cases h2 : q < init.size + 4 * m
      · simp only [h1, h2, show q < init.size + 4 * (m + 1) by omega, ↓reduceIte]
      · simp only [h1, h2, ↓reduceIte]
        by_cases h3 : q < init.size + 4 * (m + 1)
        · simp only [h3, ↓reduceIte]
          have e1 : (q - init.size) / 4 = m := by omega
          have e2 : (q - init.size) % 4 = q - (init.size + 4 * m) := by omega
          rw [e1, e2]
        · simp only [h3, ↓reduceIte]
          exact Array.getElem?_eq_none (by rw [hf]; omega)

/-! ### Closed forms -/

/-- The four matches of edge `ne` (local index `1 ≤ ne < n`) of an `S` layout. -/
def blockS (neSt neEn ne : Nat) : Array (Option Nat) :=
  let (a, b) := if ne - 1 == neSt then (4 * neSt + 1, 4 * neSt + 0) else (4 * (ne - 1) + 3, 4 * (ne - 1) + 2)
  let (c, d) := if ne + 1 == neEn then (4 * neSt + 3, 4 * neSt + 2) else (4 * (ne + 1) + 1, 4 * (ne + 1) + 0)
  #[some a, some b, some c, some d]

/-- The four matches of edge `neSt + k` of a `P` layout. -/
def blockP (neSt neEn k : Nat) : Array (Option Nat) :=
  let ne := neSt + k
  let prv := (if ne == neSt then neEn else ne) - 1
  let nxt := if ne + 1 == neEn then neSt else ne + 1
  #[some (4 * prv + 1), some (4 * nxt + 0), some (4 * nxt + 3), some (4 * prv + 2)]

theorem blockS_size (neSt neEn ne : Nat) : (blockS neSt neEn ne).size = 4 := by
  unfold blockS; split <;> split <;> rfl

theorem blockP_size (neSt neEn k : Nat) : (blockP neSt neEn k).size = 4 := rfl

theorem layoutRot_S_eq (n neSt : Nat) (hn : n ≠ 1) (ev : List Nat) (mr : Nat → Array (Option Nat)) (cv : Nat) :
    layoutRot .S n neSt (neSt + n) ev mr cv =
      (List.range (n - 1)).foldl (fun acc i => acc ++ blockS neSt (neSt + n) (neSt + 1 + i))
        #[some (4 * (neSt + 1) + 1), some (4 * (neSt + 1) + 0), some (4 * (neSt + n - 1) + 3), some (4 * (neSt + n - 1) + 2)] := by
  unfold layoutRot
  simp only [beq_iff_eq, hn, ite_false, show (NodeType.S == NodeType.Q) = false from rfl,
    show (NodeType.S == NodeType.I) = false from rfl, show (NodeType.S == NodeType.P) = false from rfl,
    Bool.false_or, Bool.false_eq_true, show (NodeType.S == NodeType.S) = true from rfl, ite_true,
    Nat.add_sub_cancel_left]
  congr 1
  funext acc i
  simp only [blockS, beq_iff_eq]

theorem layoutRot_P_eq (k neSt : Nat) (ev : List Nat) (mr : Nat → Array (Option Nat)) (cv : Nat) :
    layoutRot .P 2 neSt (neSt + k) ev mr cv =
      (List.range k).foldl (fun acc i => acc ++ blockP neSt (neSt + k) i) #[] := by
  unfold layoutRot
  simp only [beq_iff_eq, show (2 : Nat) ≠ 1 by omega, ite_false, show (NodeType.P == NodeType.Q) = false from rfl,
    show (NodeType.P == NodeType.I) = false from rfl, show (NodeType.P == NodeType.P) = true from rfl,
    Bool.false_or, Bool.false_eq_true, ite_true, Nat.add_sub_cancel_left]
  congr 1
  funext acc i
  simp only [blockP, beq_iff_eq]

/-- The match of quarter-edge `4 * ne + r` in the `S` layout of `n` edges starting at `neSt`. -/
def rotS (n neSt ne r : Nat) : Nat :=
  if ne = 0 then
    if r = 0 then 4 * (neSt + 1) + 1
    else if r = 1 then 4 * (neSt + 1)
    else if r = 2 then 4 * (neSt + n - 1) + 3
    else 4 * (neSt + n - 1) + 2
  else
    if r = 0 then (if ne = 1 then 4 * neSt + 1 else 4 * (neSt + ne - 1) + 3)
    else if r = 1 then (if ne = 1 then 4 * neSt else 4 * (neSt + ne - 1) + 2)
    else if r = 2 then (if ne + 1 = n then 4 * neSt + 3 else 4 * (neSt + ne + 1) + 1)
    else (if ne + 1 = n then 4 * neSt + 2 else 4 * (neSt + ne + 1))

/-- The match of quarter-edge `4 * ne + r` in the `P` layout of `k` edges starting at `neSt`. -/
def rotP (k neSt ne r : Nat) : Nat :=
  let prv := if ne = 0 then k - 1 else ne - 1
  let nxt := if ne + 1 = k then 0 else ne + 1
  if r = 0 then 4 * (neSt + prv) + 1
  else if r = 1 then 4 * (neSt + nxt)
  else if r = 2 then 4 * (neSt + nxt) + 3
  else 4 * (neSt + prv) + 2

theorem layoutRot_S_size (n neSt : Nat) (hn : 2 ≤ n) (ev mr cv) :
    (layoutRot .S n neSt (neSt + n) ev mr cv).size = 4 * n := by
  rw [layoutRot_S_eq n neSt (by omega),
    foldl_append_blocks_size (fun i => blockS neSt (neSt + n) (neSt + 1 + i)) (fun _ => blockS_size _ _ _)]
  simp; omega

theorem layoutRot_S_get (n neSt : Nat) (hn : 2 ≤ n) (ev mr cv) (ne r : Nat) (hne : ne < n) (hr : r < 4) :
    (layoutRot .S n neSt (neSt + n) ev mr cv)[4 * ne + r]? = some (some (rotS n neSt ne r)) := by
  rw [layoutRot_S_eq n neSt (by omega),
    foldl_append_blocks_get (fun i => blockS neSt (neSt + n) (neSt + 1 + i)) (fun _ => blockS_size _ _ _)]
  simp only [List.size_toArray, List.length_cons, List.length_nil]
  by_cases h0 : ne = 0
  · subst h0
    rw [if_pos (by omega)]
    unfold rotS
    interval_cases r <;> simp
  · rw [if_neg (by omega), if_pos (by omega)]
    have e1 : (4 * ne + r - 4) / 4 = ne - 1 := by omega
    have e2 : (4 * ne + r - 4) % 4 = r := by omega
    rw [e1, e2]
    unfold blockS rotS
    simp only [beq_iff_eq, if_neg h0]
    have c1 : (neSt + 1 + (ne - 1) - 1 = neSt) ↔ ne = 1 := by omega
    have c2 : (neSt + 1 + (ne - 1) + 1 = neSt + n) ↔ ne + 1 = n := by omega
    simp only [c1, c2]
    have e3 : 4 * (neSt + 1 + (ne - 1) - 1) = 4 * (neSt + ne - 1) := by omega
    have e4 : 4 * (neSt + 1 + (ne - 1) + 1) = 4 * (neSt + ne + 1) := by omega
    interval_cases r <;> simp [e3, e4] <;> (try split_ifs) <;> omega

theorem layoutRot_P_size (k neSt : Nat) (ev mr cv) :
    (layoutRot .P 2 neSt (neSt + k) ev mr cv).size = 4 * k := by
  rw [layoutRot_P_eq, foldl_append_blocks_size (fun i => blockP neSt (neSt + k) i) (fun _ => blockP_size _ _ _)]
  simp

theorem layoutRot_P_get (k neSt : Nat) (ev mr cv) (ne r : Nat) (hne : ne < k) (hr : r < 4) :
    (layoutRot .P 2 neSt (neSt + k) ev mr cv)[4 * ne + r]? = some (some (rotP k neSt ne r)) := by
  rw [layoutRot_P_eq, foldl_append_blocks_get (fun i => blockP neSt (neSt + k) i) (fun _ => blockP_size _ _ _)]
  simp only [Array.size_empty, Nat.zero_add, Nat.not_lt_zero, ite_false, Nat.sub_zero]
  rw [if_pos (by omega)]
  have e1 : (4 * ne + r) / 4 = ne := by omega
  have e2 : (4 * ne + r) % 4 = r := by omega
  rw [e1, e2]
  unfold blockP rotP
  have c0 : (neSt + ne = neSt) ↔ ne = 0 := by omega
  have c1 : (neSt + ne + 1 = neSt + k) ↔ ne + 1 = k := by omega
  simp only [beq_iff_eq, c0, c1]
  interval_cases r <;> simp <;> (try split_ifs) <;> omega

/-! ### Edge / side decomposition of the closed forms -/

/-- Edge of `rotS n 0 ne r`. -/
def rotSE (n ne r : Nat) : Nat :=
  if ne = 0 then (if r < 2 then 1 else n - 1)
  else if r < 2 then ne - 1 else if ne + 1 = n then 0 else ne + 1

/-- Residue of `rotS n 0 ne r`. -/
def rotSR (n ne r : Nat) : Nat :=
  if ne = 0 then (if r < 2 then 1 - r else 5 - r)
  else if r < 2 then (if ne = 1 then 1 - r else 3 - r)
  else if ne + 1 = n then 5 - r else 3 - r

/-- Edge of `rotP k 0 ne r`. -/
def rotPE (k ne r : Nat) : Nat :=
  if r = 0 ∨ r = 3 then (if ne = 0 then k - 1 else ne - 1) else (if ne + 1 = k then 0 else ne + 1)

/-- Residue of `rotP k 0 ne r`: the other side's quarter-edge with the opposite direction. -/
def rotR (r : Nat) : Nat := if r < 2 then 1 - r else 5 - r

theorem rotS_eq (n ne r : Nat) (hr : r < 4) : rotS n 0 ne r = 4 * rotSE n ne r + rotSR n ne r := by
  unfold rotS rotSE rotSR
  interval_cases r <;> split_ifs <;> first | contradiction | omega

theorem rotP_eq (k ne r : Nat) (hr : r < 4) : rotP k 0 ne r = 4 * rotPE k ne r + rotR r := by
  unfold rotP rotPE rotR
  dsimp only
  interval_cases r <;> split_ifs <;> first | contradiction | omega

theorem rotSE_lt (n ne r : Nat) (hn : 2 ≤ n) (hne : ne < n) : rotSE n ne r < n := by
  unfold rotSE; split_ifs <;> omega
theorem rotSR_lt (n ne r : Nat) (hr : r < 4) : rotSR n ne r < 4 := by
  unfold rotSR; split_ifs <;> omega
theorem rotPE_lt (k ne r : Nat) (hne : ne < k) : rotPE k ne r < k := by
  unfold rotPE; split_ifs <;> omega
theorem rotR_lt (r : Nat) (hr : r < 4) : rotR r < 4 := by
  unfold rotR; split_ifs <;> omega

theorem exists_decomp (q : Nat) : ∃ ne r, r < 4 ∧ q = 4 * ne + r :=
  ⟨q / 4, q % 4, Nat.mod_lt _ (by omega), (Nat.div_add_mod q 4).symm⟩

/-! ### The skeletons and local rotation systems -/

/-- The skeleton of an `S` node on `n` vertices: the cap `(0, n-1)` followed by the path edges. -/
def cycleEdges (n : Nat) : List (Nat × Nat) := (0, n - 1) :: (List.range (n - 1)).map fun i => (i, i + 1)

/-- The skeleton of a `P` node with `k` edges. -/
def bondEdges (k : Nat) : List (Nat × Nat) := List.replicate k (0, 1)

/-- `layoutRot .S` with node-edges renumbered from `0`. -/
def cycleRot (n : Nat) : RotationSystem := ⟨layoutRot .S n 0 n [] (fun _ => #[]) 0⟩

/-- `layoutRot .P` with node-edges renumbered from `0`. -/
def bondRot (k : Nat) : RotationSystem := ⟨layoutRot .P 2 0 k [] (fun _ => #[]) 0⟩

theorem cycleEdges_length (n : Nat) : (cycleEdges n).length = n - 1 + 1 := by simp [cycleEdges]
theorem bondEdges_length (k : Nat) : (bondEdges k).length = k := by simp [bondEdges]

theorem RotationSystem.get_eq_none_of_ge (rs : RotationSystem) (q : Nat) (h : rs.size ≤ q) : rs.get q = none := by
  unfold RotationSystem.get; rw [Array.getElem?_eq_none h]; rfl

theorem RotationSystem.stepFn_vertexStep_of_ge (rs : RotationSystem) (q : Nat) (h : rs.size ≤ QE.flipDir q) :
    stepFn rs.vertexStep q = q := by
  unfold stepFn RotationSystem.vertexStep; rw [rs.get_eq_none_of_ge _ h]; rfl

theorem RotationSystem.stepFn_faceStep_of_ge (rs : RotationSystem) (q : Nat) (h : rs.size ≤ QE.across q) :
    stepFn rs.faceStep q = q := by
  unfold stepFn RotationSystem.faceStep; rw [rs.get_eq_none_of_ge _ h]; rfl

theorem cycleRot_size (n : Nat) (hn : 2 ≤ n) : (cycleRot n).size = 4 * n := by
  unfold cycleRot RotationSystem.size
  have := layoutRot_S_size n 0 hn [] (fun _ => #[]) 0
  rwa [Nat.zero_add] at this

theorem cycleRot_get (n : Nat) (hn : 2 ≤ n) (ne r : Nat) (hne : ne < n) (hr : r < 4) :
    (cycleRot n).get (4 * ne + r) = some (rotS n 0 ne r) := by
  unfold cycleRot RotationSystem.get
  have := layoutRot_S_get n 0 hn [] (fun _ => #[]) 0 ne r hne hr
  rw [Nat.zero_add] at this
  simp only [this]; rfl

theorem bondRot_size (k : Nat) : (bondRot k).size = 4 * k := by
  unfold bondRot RotationSystem.size
  have := layoutRot_P_size k 0 [] (fun _ => #[]) 0
  rwa [Nat.zero_add] at this

theorem bondRot_get (k ne r : Nat) (hne : ne < k) (hr : r < 4) :
    (bondRot k).get (4 * ne + r) = some (rotP k 0 ne r) := by
  unfold bondRot RotationSystem.get
  have := layoutRot_P_get k 0 [] (fun _ => #[]) 0 ne r hne hr
  rw [Nat.zero_add] at this
  simp only [this]; rfl

/-- `flipDir` on a residue: `0 ↔ 1`, `2 ↔ 3`. -/
def flipR (r : Nat) : Nat := if r % 2 = 0 then r + 1 else r - 1

theorem flipDir_mk (ne r : Nat) (hr : r < 4) : QE.flipDir (4 * ne + r) = 4 * ne + flipR r := by
  rw [QE.flipDir_eq]; unfold flipR; split_ifs <;> omega

theorem across_mk (ne r : Nat) (hr : r < 4) : QE.across (4 * ne + r) = 4 * ne + (3 - r) := by
  rw [QE.across_eq]; omega

theorem cycleRot_vstep (n : Nat) (hn : 2 ≤ n) (ne r : Nat) (hne : ne < n) (hr : r < 4) :
    stepFn (cycleRot n).vertexStep (4 * ne + r) = 4 * rotSE n ne (flipR r) + rotSR n ne (flipR r) := by
  unfold stepFn RotationSystem.vertexStep
  rw [flipDir_mk ne r hr, cycleRot_get n hn ne _ hne (by unfold flipR; split_ifs <;> omega),
    rotS_eq n ne _ (by unfold flipR; split_ifs <;> omega)]
  rfl

theorem cycleRot_fstep (n : Nat) (hn : 2 ≤ n) (ne r : Nat) (hne : ne < n) (hr : r < 4) :
    stepFn (cycleRot n).faceStep (4 * ne + r) = 4 * rotSE n ne (3 - r) + rotSR n ne (3 - r) := by
  unfold stepFn RotationSystem.faceStep
  rw [across_mk ne r hr, cycleRot_get n hn ne _ hne (by omega), rotS_eq n ne _ (by omega)]
  rfl

theorem bondRot_vstep (k ne r : Nat) (hne : ne < k) (hr : r < 4) :
    stepFn (bondRot k).vertexStep (4 * ne + r) = 4 * rotPE k ne (flipR r) + rotR (flipR r) := by
  unfold stepFn RotationSystem.vertexStep
  rw [flipDir_mk ne r hr, bondRot_get k ne _ hne (by unfold flipR; split_ifs <;> omega),
    rotP_eq k ne _ (by unfold flipR; split_ifs <;> omega)]
  rfl

theorem bondRot_fstep (k ne r : Nat) (hne : ne < k) (hr : r < 4) :
    stepFn (bondRot k).faceStep (4 * ne + r) = 4 * rotPE k ne (3 - r) + rotR (3 - r) := by
  unfold stepFn RotationSystem.faceStep
  rw [across_mk ne r hr, bondRot_get k ne _ hne (by omega), rotP_eq k ne _ (by omega)]
  rfl

/-! ### Residue arithmetic -/

theorem flipR_lt (r : Nat) (hr : r < 4) : flipR r < 4 := by unfold flipR; split_ifs <;> omega
theorem flipR_flipR (r : Nat) (hr : r < 4) : flipR (flipR r) = r := by unfold flipR; split_ifs <;> omega

/-- `rotS` is an involution. -/
theorem cycle_inv (n ne r : Nat) (hn : 2 ≤ n) (hne : ne < n) (hr : r < 4) :
    rotSE n (rotSE n ne r) (rotSR n ne r) = ne ∧ rotSR n (rotSE n ne r) (rotSR n ne r) = r := by
  unfold rotSE rotSR
  interval_cases r <;> split_ifs <;> first | contradiction | omega

theorem rotSE_flipR (n ne r : Nat) (hr : r < 4) : rotSE n ne (flipR r) = rotSE n ne r := by
  unfold rotSE flipR
  interval_cases r <;> split_ifs <;> first | contradiction | omega

theorem rotSR_flipR (n ne r : Nat) (hr : r < 4) : rotSR n ne (flipR r) = flipR (rotSR n ne r) := by
  unfold rotSR flipR
  interval_cases r <;> split_ifs <;> first | contradiction | omega

/-- The vertex step of the cycle is an involution. -/
theorem cycle_vinv (n ne r : Nat) (hn : 2 ≤ n) (hne : ne < n) (hr : r < 4) :
    rotSE n (rotSE n ne r) (flipR (rotSR n ne r)) = ne ∧
    rotSR n (rotSE n ne r) (flipR (rotSR n ne r)) = flipR r := by
  rw [rotSE_flipR _ _ _ (rotSR_lt n ne r hr), rotSR_flipR _ _ _ (rotSR_lt n ne r hr),
    (cycle_inv n ne r hn hne hr).1, (cycle_inv n ne r hn hne hr).2]
  exact ⟨rfl, rfl⟩

/-- `rotP` is an involution. -/
theorem bond_inv (k ne r : Nat) (hne : ne < k) (hr : r < 4) :
    rotPE k (rotPE k ne r) (rotR r) = ne ∧ rotR (rotR r) = r := by
  unfold rotPE rotR
  interval_cases r <;> split_ifs <;> first | contradiction | omega

/-- The face step of the bond is an involution. -/
theorem bond_finv (k ne r : Nat) (hne : ne < k) (hr : r < 4) :
    rotPE k (rotPE k ne r) (3 - rotR r) = ne ∧ rotR (3 - rotR r) = 3 - r := by
  unfold rotPE rotR
  interval_cases r <;> split_ifs <;> first | contradiction | omega

/-! ### Cycle: vertex step is a fixed-point-free involution, face step has four orbits -/

section Cycle
variable (n : Nat) (hn : 2 ≤ n)
include hn

theorem cycleRot_vstep_lt (q : Nat) (hq : q < 4 * n) : stepFn (cycleRot n).vertexStep q < 4 * n := by
  obtain ⟨ne, r, hr, rfl⟩ := exists_decomp q
  rw [cycleRot_vstep n hn ne r (by omega) hr]
  have := rotSE_lt n ne (flipR r) hn (by omega)
  have := rotSR_lt n ne (flipR r) (flipR_lt r hr)
  omega

theorem cycleRot_vstep_invol (q : Nat) (hq : q < 4 * n) :
    stepFn (cycleRot n).vertexStep (stepFn (cycleRot n).vertexStep q) = q := by
  obtain ⟨ne, r, hr, rfl⟩ := exists_decomp q
  have hne : ne < n := by omega
  rw [cycleRot_vstep n hn ne r hne hr,
    cycleRot_vstep n hn _ _ (rotSE_lt n ne _ hn hne) (rotSR_lt n ne _ (flipR_lt r hr)),
    (cycle_vinv n ne (flipR r) hn hne (flipR_lt r hr)).1,
    (cycle_vinv n ne (flipR r) hn hne (flipR_lt r hr)).2, flipR_flipR r hr]

theorem cycleRot_vstep_ne (q : Nat) (hq : q < 4 * n) : stepFn (cycleRot n).vertexStep q ≠ q := by
  obtain ⟨ne, r, hr, rfl⟩ := exists_decomp q
  rw [cycleRot_vstep n hn ne r (by omega) hr]
  unfold rotSE rotSR flipR
  interval_cases r <;> split_ifs <;> first | contradiction | omega

theorem cycleRot_vertexOrbits : (cycleRot n).numVertexOrbits = 2 * n := by
  unfold RotationSystem.numVertexOrbits
  rw [cycleRot_size n hn, numOrbits_involution _ _ (cycleRot_vstep_lt n hn) (cycleRot_vstep_invol n hn)
    (cycleRot_vstep_ne n hn)]
  omega

/-- The face of a cycle quarter-edge: the cap's residues are swapped (`r ↦ r + 2`). -/
def cycleFaceLabel (q : Nat) : Nat := if q / 4 = 0 then (q % 4 + 2) % 4 else q % 4

theorem cycleRot_fstep_label (x : Nat) :
    cycleFaceLabel (stepFn (cycleRot n).faceStep x) = cycleFaceLabel x := by
  by_cases hx : x < 4 * n
  · obtain ⟨ne, r, hr, rfl⟩ := exists_decomp x
    rw [cycleRot_fstep n hn ne r (by omega) hr]
    unfold cycleFaceLabel rotSE rotSR
    interval_cases r <;> split_ifs <;> first | contradiction | omega
  · rw [(cycleRot n).stepFn_faceStep_of_ge x (by rw [cycleRot_size n hn, QE.across_eq]; omega)]

theorem cycleRot_fstep_lt (q : Nat) (hq : q < 4 * n) : stepFn (cycleRot n).faceStep q < 4 * n := by
  obtain ⟨ne, r, hr, rfl⟩ := exists_decomp q
  rw [cycleRot_fstep n hn ne r (by omega) hr]
  have := rotSE_lt n ne (3 - r) hn (by omega)
  have := rotSR_lt n ne (3 - r) (by omega)
  omega

theorem cycleRot_fstep_eq (ne r : Nat) (hne : ne < n) (hr : r < 4) :
    stepFn (cycleRot n).faceStep (4 * ne + r) =
      if ne = 0 then (if r < 2 then 4 * (n - 1) + (r + 2) else 4 + (r - 2))
      else if r < 2 then (if ne + 1 = n then r + 2 else 4 * (ne + 1) + r)
      else (if ne = 1 then r - 2 else 4 * (ne - 1) + r) := by
  rw [cycleRot_fstep n hn ne r hne hr]
  unfold rotSE rotSR
  interval_cases r <;> split_ifs <;> first | contradiction | omega

theorem cycleRot_fiter_lo (ne r : Nat) (hne : 1 ≤ ne) (hr : r < 2) :
    ∀ i, ne + i + 1 ≤ n → (stepFn (cycleRot n).faceStep)^[i] (4 * ne + r) = 4 * (ne + i) + r := by
  intro i
  induction i with
  | zero => intro _; simp
  | succ i ih =>
    intro h
    rw [Function.iterate_succ_apply', ih (by omega), cycleRot_fstep_eq n hn _ _ (by omega) (by omega)]
    split_ifs <;> omega

theorem cycleRot_fiter_hi (ne r : Nat) (hr : 2 ≤ r) (hr4 : r < 4) :
    ∀ i, i + 1 ≤ ne → ne < n → (stepFn (cycleRot n).faceStep)^[i] (4 * ne + r) = 4 * (ne - i) + r := by
  intro i
  induction i with
  | zero => intro _ _; simp
  | succ i ih =>
    intro h hne
    rw [Function.iterate_succ_apply', ih (by omega) hne, cycleRot_fstep_eq n hn _ _ (by omega) (by omega)]
    split_ifs <;> omega

theorem cycleRot_isOrbitMin_iff (q : Nat) (hq : q < 4 * n) :
    isOrbitMin (cycleRot n).faceStep (4 * n) q = true ↔ q < 4 := by
  constructor
  · intro h
    by_contra h4
    obtain ⟨ne, r, hr, rfl⟩ := exists_decomp q
    have hne : 1 ≤ ne := by omega
    have hne' : ne < n := by omega
    suffices hs : isOrbitMin (cycleRot n).faceStep (4 * n) (4 * ne + r) = false by rw [h] at hs; cases hs
    by_cases hr2 : r < 2
    · apply not_isOrbitMin_of_iter _ _ _ (n - ne) (by omega)
      rw [show n - ne = (n - 1 - ne) + 1 by omega, Function.iterate_succ_apply',
        cycleRot_fiter_lo n hn ne r hne hr2 _ (by omega), cycleRot_fstep_eq n hn _ _ (by omega) (by omega)]
      split_ifs <;> omega
    · obtain ⟨m, rfl⟩ : ∃ m, ne = m + 1 := ⟨ne - 1, by omega⟩
      apply not_isOrbitMin_of_iter _ _ _ (m + 1) (by omega)
      rw [Function.iterate_succ_apply', cycleRot_fiter_hi n hn (m + 1) r (by omega) hr m (by omega) hne',
        show m + 1 - m = 1 by omega, cycleRot_fstep_eq n hn 1 r (by omega) hr]
      split_ifs <;> omega
  · intro h4
    apply isOrbitMin_of_label _ _ _ cycleFaceLabel (cycleRot_fstep_label n hn)
    intro x hx
    unfold cycleFaceLabel at hx
    split_ifs at hx <;> omega

theorem cycleRot_faceOrbits : (cycleRot n).numFaceOrbits = 4 := by
  unfold RotationSystem.numFaceOrbits
  rw [cycleRot_size n hn, numOrbits_eq_card _ _ (fun q => q < 4) (cycleRot_isOrbitMin_iff n hn)]
  have : (Finset.range (4 * n)).filter (fun q => q < 4) = Finset.range 4 := by
    ext q; simp only [Finset.mem_filter, Finset.mem_range]; omega
  rw [this, Finset.card_range]

end Cycle

/-! ### Bond: face step is a fixed-point-free involution, vertex step has four orbits -/

section Bond
variable (k : Nat) (hk : 1 ≤ k)
include hk

theorem bondRot_fstep_lt (q : Nat) (hq : q < 4 * k) : stepFn (bondRot k).faceStep q < 4 * k := by
  obtain ⟨ne, r, hr, rfl⟩ := exists_decomp q
  rw [bondRot_fstep k ne r (by omega) hr]
  have := rotPE_lt k ne (3 - r) (by omega)
  have := rotR_lt (3 - r) (by omega)
  omega

theorem bondRot_fstep_invol (q : Nat) (hq : q < 4 * k) :
    stepFn (bondRot k).faceStep (stepFn (bondRot k).faceStep q) = q := by
  obtain ⟨ne, r, hr, rfl⟩ := exists_decomp q
  have hne : ne < k := by omega
  rw [bondRot_fstep k ne r hne hr, bondRot_fstep k _ _ (rotPE_lt k ne _ hne) (rotR_lt _ (by omega)),
    (bond_finv k ne (3 - r) hne (by omega)).1, (bond_finv k ne (3 - r) hne (by omega)).2]
  omega

theorem bondRot_fstep_ne (q : Nat) (hq : q < 4 * k) : stepFn (bondRot k).faceStep q ≠ q := by
  obtain ⟨ne, r, hr, rfl⟩ := exists_decomp q
  rw [bondRot_fstep k ne r (by omega) hr]
  unfold rotPE rotR
  interval_cases r <;> split_ifs <;> first | contradiction | omega

theorem bondRot_faceOrbits : (bondRot k).numFaceOrbits = 2 * k := by
  unfold RotationSystem.numFaceOrbits
  rw [bondRot_size k, numOrbits_involution _ _ (bondRot_fstep_lt k hk) (bondRot_fstep_invol k hk)
    (bondRot_fstep_ne k hk)]
  omega

theorem bondRot_vstep_label (x : Nat) : stepFn (bondRot k).vertexStep x % 4 = x % 4 := by
  by_cases hx : x < 4 * k
  · obtain ⟨ne, r, hr, rfl⟩ := exists_decomp x
    rw [bondRot_vstep k ne r (by omega) hr]
    unfold rotPE rotR flipR
    interval_cases r <;> split_ifs <;> first | contradiction | omega
  · rw [(bondRot k).stepFn_vertexStep_of_ge x (by rw [bondRot_size k, QE.flipDir_eq]; omega)]

theorem bondRot_vstep_eq (ne r : Nat) (hne : ne < k) (hr : r < 4) :
    stepFn (bondRot k).vertexStep (4 * ne + r) =
      if r = 0 ∨ r = 3 then (if ne + 1 = k then r else 4 * (ne + 1) + r)
      else (if ne = 0 then 4 * (k - 1) + r else 4 * (ne - 1) + r) := by
  rw [bondRot_vstep k ne r hne hr]
  unfold rotPE rotR flipR
  interval_cases r <;> split_ifs <;> first | contradiction | omega

theorem bondRot_viter_up (ne r : Nat) (hr : r = 0 ∨ r = 3) :
    ∀ i, ne + i + 1 ≤ k → (stepFn (bondRot k).vertexStep)^[i] (4 * ne + r) = 4 * (ne + i) + r := by
  intro i
  induction i with
  | zero => intro _; simp
  | succ i ih =>
    intro h
    rw [Function.iterate_succ_apply', ih (by omega), bondRot_vstep_eq k hk _ _ (by omega) (by omega)]
    split_ifs <;> omega

theorem bondRot_viter_down (ne r : Nat) (hr : r = 1 ∨ r = 2) (hne : ne < k) :
    ∀ i, i ≤ ne → (stepFn (bondRot k).vertexStep)^[i] (4 * ne + r) = 4 * (ne - i) + r := by
  intro i
  induction i with
  | zero => intro _; simp
  | succ i ih =>
    intro h
    rw [Function.iterate_succ_apply', ih (by omega), bondRot_vstep_eq k hk _ _ (by omega) (by omega)]
    split_ifs <;> omega

theorem bondRot_isOrbitMin_iff (q : Nat) (hq : q < 4 * k) :
    isOrbitMin (bondRot k).vertexStep (4 * k) q = true ↔ q < 4 := by
  constructor
  · intro h
    by_contra h4
    obtain ⟨ne, r, hr, rfl⟩ := exists_decomp q
    have hne : 1 ≤ ne := by omega
    have hne' : ne < k := by omega
    suffices hs : isOrbitMin (bondRot k).vertexStep (4 * k) (4 * ne + r) = false by rw [h] at hs; cases hs
    by_cases hr2 : r = 0 ∨ r = 3
    · apply not_isOrbitMin_of_iter _ _ _ (k - ne) (by omega)
      rw [show k - ne = (k - 1 - ne) + 1 by omega, Function.iterate_succ_apply',
        bondRot_viter_up k hk ne r hr2 _ (by omega), bondRot_vstep_eq k hk _ _ (by omega) (by omega)]
      split_ifs <;> omega
    · apply not_isOrbitMin_of_iter _ _ _ ne (by omega)
      rw [bondRot_viter_down k hk ne r (by omega) hne' ne le_rfl]
      omega
  · intro h4
    apply isOrbitMin_of_label _ _ _ (· % 4) (bondRot_vstep_label k hk)
    intro x hx
    omega

theorem bondRot_vertexOrbits : (bondRot k).numVertexOrbits = 4 := by
  unfold RotationSystem.numVertexOrbits
  rw [bondRot_size k, numOrbits_eq_card _ _ (fun q => q < 4) (bondRot_isOrbitMin_iff k hk)]
  have : (Finset.range (4 * k)).filter (fun q => q < 4) = Finset.range 4 := by
    ext q; simp only [Finset.mem_filter, Finset.mem_range]; omega
  rw [this, Finset.card_range]

end Bond

/-! ### Non-isolated vertices and components of the cycle and the bond -/

theorem nonIsolated_cycle (n : Nat) (hn : 2 ≤ n) (v : Nat) (hv : v < n) : nonIsolated (cycleEdges n) v = true := by
  unfold nonIsolated
  rw [List.any_eq_true]
  by_cases h : v = n - 1
  · exact ⟨(0, n - 1), List.mem_cons_self, by simp [h]⟩
  · exact ⟨(v, v + 1), List.mem_cons_of_mem _ (List.mem_map.mpr ⟨v, List.mem_range.mpr (by omega), rfl⟩), by simp⟩

theorem nonIsolated_bond (k : Nat) (hk : 1 ≤ k) (v : Nat) (hv : v < 2) : nonIsolated (bondEdges k) v = true := by
  unfold nonIsolated
  rw [List.any_eq_true]
  refine ⟨(0, 1), List.mem_replicate.mpr ⟨by omega, rfl⟩, ?_⟩
  interval_cases v <;> simp

theorem numNonIsolated_cycle (n : Nat) (hn : 2 ≤ n) : numNonIsolated (cycleEdges n) n = n := by
  unfold numNonIsolated
  rw [List.filter_eq_self.mpr (fun v hv => nonIsolated_cycle n hn v (List.mem_range.mp hv)), List.length_range]

theorem numNonIsolated_bond (k : Nat) (hk : 1 ≤ k) : numNonIsolated (bondEdges k) 2 = 2 := by
  unfold numNonIsolated
  rw [List.filter_eq_self.mpr (fun v hv => nonIsolated_bond k hk v (List.mem_range.mp hv)), List.length_range]

/-- One relaxation step along one edge. -/
def relaxEdge (l : List Nat) (p : Nat × Nat) : List Nat :=
  let m := min (l[p.1]?.getD p.1) (l[p.2]?.getD p.2)
  (l.set p.1 m).set p.2 m

theorem relaxLabels_eq (es : List (Nat × Nat)) (l : List Nat) : relaxLabels es l = es.foldl relaxEdge l := rfl

/-- Labels of length `n` that are `0` on every vertex `v < n` with `v ≤ j ∨ v = n - 1`. -/
def ZeroUpTo (n j : Nat) (l : List Nat) : Prop := l.length = n ∧ ∀ v, v < n → (v ≤ j ∨ v = n - 1) → l[v]? = some 0

/-- All-zero labels. -/
def AllZero (n : Nat) (l : List Nat) : Prop := l.length = n ∧ ∀ v, v < n → l[v]? = some 0

theorem relaxEdge_set_zero (n : Nat) (l : List Nat) (hl : l.length = n) (a b : Nat) (ha : a < n) (hb : b < n)
    (h0 : l[a]? = some 0 ∨ l[b]? = some 0) (v : Nat) (hv : v < n) :
    (relaxEdge l (a, b))[v]? = if v = a ∨ v = b then some 0 else l[v]? := by
  have hm : min (l[a]?.getD a) (l[b]?.getD b) = 0 := by
    rcases h0 with h | h <;> rw [h] <;> simp
  simp only [relaxEdge, hm]
  by_cases h2 : b = v
  · subst h2; rw [List.getElem?_set_self (by simp [hl]; omega)]; simp
  · rw [List.getElem?_set_ne h2]
    by_cases h1 : a = v
    · subst h1; rw [List.getElem?_set_self (by omega)]; simp
    · rw [List.getElem?_set_ne h1, if_neg (by omega)]

theorem relaxEdge_length (l : List Nat) (p : Nat × Nat) : (relaxEdge l p).length = l.length := by
  simp [relaxEdge]

theorem relaxEdge_allZero (n : Nat) (l : List Nat) (p : Nat × Nat) (h : AllZero n l) (hp : p.1 < n ∧ p.2 < n) :
    AllZero n (relaxEdge l p) := by
  obtain ⟨hl, hz⟩ := h
  refine ⟨by rw [relaxEdge_length, hl], fun v hv => ?_⟩
  rw [relaxEdge_set_zero n l hl p.1 p.2 hp.1 hp.2 (Or.inl (hz _ hp.1)) v hv]
  split_ifs
  · rfl
  · exact hz v hv

theorem foldl_relaxEdge_allZero (n : Nat) : ∀ (es : List (Nat × Nat)), (∀ p ∈ es, p.1 < n ∧ p.2 < n) →
    ∀ l, AllZero n l → AllZero n (es.foldl relaxEdge l) := by
  intro es
  induction es with
  | nil => intro _ l h; exact h
  | cons p es ih =>
    intro hes l h
    rw [List.foldl_cons]
    exact ih (fun q hq => hes q (List.mem_cons_of_mem _ hq)) _ (relaxEdge_allZero n l p h (hes p List.mem_cons_self))

theorem iterate_relaxLabels_allZero (n : Nat) (es : List (Nat × Nat)) (hes : ∀ p ∈ es, p.1 < n ∧ p.2 < n) :
    ∀ (k : Nat) (l : List Nat), AllZero n l → AllZero n ((relaxLabels es)^[k] l) := by
  intro k
  induction k with
  | zero => intro l h; exact h
  | succ k ih =>
    intro l h
    rw [Function.iterate_succ_apply]
    exact ih _ (foldl_relaxEdge_allZero n es hes l h)

theorem numComponents_eq_one (es : List (Nat × Nat)) (nVerts : Nat) (h0 : 0 < nVerts)
    (hz : AllZero nVerts (compLabels es nVerts)) (hni : ∀ v, v < nVerts → nonIsolated es v = true) :
    numComponents es nVerts = 1 := by
  unfold numComponents
  simp only
  rw [length_filter_range_eq_card]
  have : (Finset.range nVerts).filter
      (fun v => (nonIsolated es v && (compLabels es nVerts)[v]?.getD v == v) = true) =
      (Finset.range nVerts).filter (fun v => v = 0) := by
    apply Finset.filter_congr
    intro v hv
    rw [Finset.mem_range] at hv
    rw [hni v hv, hz.2 v hv]
    simp only [Option.getD_some, Bool.true_and, beq_iff_eq]
    exact eq_comm
  rw [this, Finset.filter_eq', if_pos (Finset.mem_range.mpr h0), Finset.card_singleton]

theorem cycleEdges_verts (n : Nat) (hn : 2 ≤ n) : ∀ p ∈ cycleEdges n, p.1 < n ∧ p.2 < n := by
  intro p hp
  simp only [cycleEdges, List.mem_cons, List.mem_map, List.mem_range] at hp
  rcases hp with rfl | ⟨i, hi, rfl⟩ <;> simp <;> omega

theorem bondEdges_verts (k : Nat) : ∀ p ∈ bondEdges k, p.1 < 2 ∧ p.2 < 2 := by
  intro p hp
  rw [bondEdges, List.mem_replicate] at hp
  rw [hp.2]; simp

theorem cycle_path_fold (n : Nat) (hn : 2 ≤ n) : ∀ j, j ≤ n - 1 → ∀ l, ZeroUpTo n 0 l →
    ZeroUpTo n j (((List.range j).map fun i => (i, i + 1)).foldl relaxEdge l) := by
  intro j
  induction j with
  | zero => intro _ l h; simpa using h
  | succ j ih =>
    intro hj l h
    rw [List.range_succ, List.map_append, List.foldl_append, List.map_singleton, List.foldl_cons, List.foldl_nil]
    obtain ⟨hl, hz⟩ := ih (by omega) l h
    refine ⟨by rw [relaxEdge_length, hl], fun v hv hv' => ?_⟩
    rw [relaxEdge_set_zero n _ hl j (j + 1) (by omega) (by omega) (Or.inl (hz j (by omega) (Or.inl le_rfl))) v hv]
    split_ifs with h
    · rfl
    · exact hz v hv (by omega)

theorem compLabels_cycle (n : Nat) (hn : 2 ≤ n) : AllZero n (compLabels (cycleEdges n) n) := by
  unfold compLabels
  rw [show n = (n - 1) + 1 by omega, Function.iterate_succ_apply, show n - 1 + 1 = n by omega]
  apply iterate_relaxLabels_allZero n _ (cycleEdges_verts n hn)
  rw [relaxLabels_eq, cycleEdges, List.foldl_cons]
  have h0 : ZeroUpTo n 0 (relaxEdge (List.range n) (0, n - 1)) := by
    refine ⟨by rw [relaxEdge_length, List.length_range], fun v hv hv' => ?_⟩
    rw [relaxEdge_set_zero n _ (List.length_range) 0 (n - 1) (by omega) (by omega)
      (Or.inl (by rw [List.getElem?_range (by omega)])) v hv]
    rw [if_pos (by omega)]
  obtain ⟨hl, hz⟩ := cycle_path_fold n hn (n - 1) le_rfl _ h0
  exact ⟨hl, fun v hv => hz v hv (Or.inl (by omega))⟩

theorem compLabels_bond (k : Nat) (hk : 1 ≤ k) : AllZero 2 (compLabels (bondEdges k) 2) := by
  unfold compLabels
  rw [Function.iterate_succ_apply]
  apply iterate_relaxLabels_allZero 2 _ (bondEdges_verts k)
  rw [relaxLabels_eq, bondEdges, show k = (k - 1) + 1 by omega, List.replicate_succ, List.foldl_cons]
  apply foldl_relaxEdge_allZero 2 _ (bondEdges_verts (k - 1))
  exact ⟨rfl, by decide⟩

theorem numComponents_cycle (n : Nat) (hn : 2 ≤ n) : numComponents (cycleEdges n) n = 1 :=
  numComponents_eq_one _ _ (by omega) (compLabels_cycle n hn) (nonIsolated_cycle n hn)

theorem numComponents_bond (k : Nat) (hk : 1 ≤ k) : numComponents (bondEdges k) 2 = 1 :=
  numComponents_eq_one _ _ (by omega) (compLabels_bond k hk) (nonIsolated_bond k hk)

/-! ### Vertices of quarter-edges -/

theorem cycleEdges_get (n ne : Nat) (hne : ne < n) :
    (cycleEdges n)[ne]? = some (if ne = 0 then (0, n - 1) else (ne - 1, ne)) := by
  cases ne with
  | zero => rfl
  | succ m => simp [cycleEdges, List.getElem?_range (show m < n - 1 by omega)]

theorem bondEdges_get (k ne : Nat) (hne : ne < k) : (bondEdges k)[ne]? = some (0, 1) :=
  List.getElem?_replicate_of_lt hne

/-- The vertex of quarter-edge `4 * ne + r` of the cycle. -/
def cycVert (n ne r : Nat) : Nat :=
  if ne = 0 then (if r / 2 = 0 then 0 else n - 1) else (if r / 2 = 0 then ne - 1 else ne)

theorem vert_cycle (n ne r : Nat) (hne : ne < n) (hr : r < 4) :
    QE.vert (cycleEdges n) (4 * ne + r) = some (cycVert n ne r) := by
  unfold QE.vert
  rw [QE.edge_mk _ _ hr, QE.side_mk _ _ hr, cycleEdges_get n ne hne, Option.map_some]
  unfold cycVert
  split_ifs <;> rfl

theorem vert_bond (k ne r : Nat) (hne : ne < k) (hr : r < 4) :
    QE.vert (bondEdges k) (4 * ne + r) = some (if r / 2 = 0 then 0 else 1) := by
  unfold QE.vert
  rw [QE.edge_mk _ _ hr, QE.side_mk _ _ hr, bondEdges_get k ne hne, Option.map_some]

theorem cycVert_rot (n ne r : Nat) (hn : 2 ≤ n) (hne : ne < n) (hr : r < 4) :
    cycVert n (rotSE n ne r) (rotSR n ne r) = cycVert n ne r := by
  unfold cycVert rotSE rotSR
  interval_cases r <;> split_ifs <;> first | contradiction | omega

theorem rotSR_parity (n ne r : Nat) (hr : r < 4) : rotSR n ne r % 2 ≠ r % 2 := by
  unfold rotSR
  interval_cases r <;> split_ifs <;> first | contradiction | omega

theorem rotR_parity (r : Nat) (hr : r < 4) : rotR r % 2 ≠ r % 2 := by
  unfold rotR
  interval_cases r <;> split_ifs <;> first | contradiction | omega

theorem rotR_side (r : Nat) (hr : r < 4) : rotR r / 2 = r / 2 := by
  unfold rotR
  interval_cases r <;> split_ifs <;> first | contradiction | omega

/-! ### The planar embeddings -/

theorem cycleRot_isPlanarEmbedding (n : Nat) (hn : 2 ≤ n) : IsPlanarEmbedding (cycleEdges n) n (cycleRot n) where
  size := by rw [cycleRot_size n hn, cycleEdges_length]; omega
  verts := cycleEdges_verts n hn
  total := by
    intro q hq
    rw [cycleRot_size n hn] at hq
    obtain ⟨ne, r, hr, rfl⟩ := exists_decomp q
    rw [cycleRot_get n hn ne r (by omega) hr]; rfl
  involution := by
    intro q hq r' hr'
    rw [cycleRot_size n hn] at hq ⊢
    obtain ⟨ne, r, hr, rfl⟩ := exists_decomp q
    have hne : ne < n := by omega
    rw [cycleRot_get n hn ne r hne hr, rotS_eq n ne r hr, Option.mem_def, Option.some.injEq] at hr'
    subst hr'
    have h1 := rotSE_lt n ne r hn hne
    have h2 := rotSR_lt n ne r hr
    refine ⟨by omega, ?_⟩
    rw [cycleRot_get n hn _ _ h1 h2, rotS_eq n _ _ h2, (cycle_inv n ne r hn hne hr).1,
      (cycle_inv n ne r hn hne hr).2]
  opposite_dir := by
    intro q hq r' hr'
    rw [cycleRot_size n hn] at hq
    obtain ⟨ne, r, hr, rfl⟩ := exists_decomp q
    rw [cycleRot_get n hn ne r (by omega) hr, rotS_eq n ne r hr, Option.mem_def, Option.some.injEq] at hr'
    subst hr'
    rw [QE.dir_mk _ _ hr, QE.dir_mk _ _ (rotSR_lt n ne r hr)]
    exact rotSR_parity n ne r hr
  same_vertex := by
    intro q hq r' hr'
    rw [cycleRot_size n hn] at hq
    obtain ⟨ne, r, hr, rfl⟩ := exists_decomp q
    have hne : ne < n := by omega
    rw [cycleRot_get n hn ne r hne hr, rotS_eq n ne r hr, Option.mem_def, Option.some.injEq] at hr'
    subst hr'
    rw [vert_cycle n ne r hne hr, vert_cycle n _ _ (rotSE_lt n ne r hn hne) (rotSR_lt n ne r hr),
      cycVert_rot n ne r hn hne hr]
  vertex_orbits := by rw [cycleRot_vertexOrbits n hn, numNonIsolated_cycle n hn]
  euler := by
    unfold EulerFormula
    rw [cycleRot_faceOrbits n hn, numNonIsolated_cycle n hn, numComponents_cycle n hn, cycleEdges_length]
    omega

theorem bondRot_isPlanarEmbedding (k : Nat) (hk : 1 ≤ k) : IsPlanarEmbedding (bondEdges k) 2 (bondRot k) where
  size := by rw [bondRot_size k, bondEdges_length]
  verts := bondEdges_verts k
  total := by
    intro q hq
    rw [bondRot_size k] at hq
    obtain ⟨ne, r, hr, rfl⟩ := exists_decomp q
    rw [bondRot_get k ne r (by omega) hr]; rfl
  involution := by
    intro q hq r' hr'
    rw [bondRot_size k] at hq ⊢
    obtain ⟨ne, r, hr, rfl⟩ := exists_decomp q
    have hne : ne < k := by omega
    rw [bondRot_get k ne r hne hr, rotP_eq k ne r hr, Option.mem_def, Option.some.injEq] at hr'
    subst hr'
    have h1 := rotPE_lt k ne r hne
    have h2 := rotR_lt r hr
    refine ⟨by omega, ?_⟩
    rw [bondRot_get k _ _ h1 h2, rotP_eq k _ _ h2, (bond_inv k ne r hne hr).1, (bond_inv k ne r hne hr).2]
  opposite_dir := by
    intro q hq r' hr'
    rw [bondRot_size k] at hq
    obtain ⟨ne, r, hr, rfl⟩ := exists_decomp q
    rw [bondRot_get k ne r (by omega) hr, rotP_eq k ne r hr, Option.mem_def, Option.some.injEq] at hr'
    subst hr'
    rw [QE.dir_mk _ _ hr, QE.dir_mk _ _ (rotR_lt r hr)]
    exact rotR_parity r hr
  same_vertex := by
    intro q hq r' hr'
    rw [bondRot_size k] at hq
    obtain ⟨ne, r, hr, rfl⟩ := exists_decomp q
    have hne : ne < k := by omega
    rw [bondRot_get k ne r hne hr, rotP_eq k ne r hr, Option.mem_def, Option.some.injEq] at hr'
    subst hr'
    rw [vert_bond k ne r hne hr, vert_bond k _ _ (rotPE_lt k ne r hne) (rotR_lt r hr), rotR_side r hr]
  vertex_orbits := by rw [bondRot_vertexOrbits k hk, numNonIsolated_bond k hk]
  euler := by
    unfold EulerFormula
    rw [bondRot_faceOrbits k hk, numNonIsolated_bond k hk, numComponents_bond k hk, bondEdges_length]
    omega

end Spqr
