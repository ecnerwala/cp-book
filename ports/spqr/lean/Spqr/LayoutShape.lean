import Spqr.StLayout

/-!
# Per-type skeleton layouts

Exact contents of `layoutNode type …` for every node type: the `edges` array, the local
`adjBounds`, the `adjDat` slots and the rows `LayoutR.row` (`2 nv` = lower neighbours of node-vert
`nv`, `2 nv + 1` = higher neighbours). `Layout.Shape` mirrors `SpqrTree.Shape` on one node's
layout so that the per-node transport is a rewrite.
-/

namespace Spqr

namespace LayoutShape

open LayoutR (row rowBound forIn_id_yield getElem!_replicate')

/-! ### Array / `setNe` lemmas -/

theorem getElem!_set!_eq {α : Type} [Inhabited α] (xs : Array α) (i j : Nat) (a : α)
    (hj : j < xs.size) :
    (xs.set! i a)[j]! = if j = i then a else xs[j]! := by
  split
  · subst ‹j = i›; exact Array.getElem!_set!_self _ _ _ hj
  · exact Array.getElem!_set!_ne _ _ _ _ (Ne.symm ‹_›)

theorem getElem!_modify_eq {α : Type} [Inhabited α] (xs : Array α) (i j : Nat) (f : α → α)
    (hj : j < xs.size) :
    (xs.modify i f)[j]! = if j = i then f xs[j]! else xs[j]! := by
  split
  · subst ‹j = i›; exact LayoutR.getElem!_modify_self' _ _ _ hj
  · exact LayoutR.getElem!_modify_ne' _ _ _ _ (Ne.symm ‹_›)

theorem ite_of_pos {α : Type} {c : Prop} [Decidable c] (h : c) (a b : α) :
    (if c then a else b) = a := by simp [h]
theorem ite_of_neg {α : Type} {c : Prop} [Decidable c] (h : ¬ c) (a b : α) :
    (if c then a else b) = b := by simp [h]

theorem empty_edges_size (n m : Nat) : (Layout.empty n m).edges.size = m := by simp [Layout.empty]
theorem empty_adjBounds_size (n m : Nat) : (Layout.empty n m).adjBounds.size = 2 * n + 1 := by
  simp [Layout.empty]
theorem empty_adjDat_size (n m : Nat) : (Layout.empty n m).adjDat.size = 2 * m := by
  simp [Layout.empty]
theorem empty_edges_get (n m k : Nat) : (Layout.empty n m).edges[k]! = default := by
  simp [Layout.empty, getElem!_replicate']
theorem empty_adjBounds_get (n m k : Nat) : (Layout.empty n m).adjBounds[k]! = 0 := by
  simp [Layout.empty, getElem!_replicate']
theorem empty_adjDat_get (n m k : Nat) : (Layout.empty n m).adjDat[k]! = default := by
  simp [Layout.empty, getElem!_replicate']

theorem setNe_adjBounds (l : Layout) (neSt node ne : Nat) (nvs nds : Nat × Nat) :
    (l.setNe neSt node ne nvs nds).adjBounds = l.adjBounds := rfl
theorem setNe_edges_size (l : Layout) (neSt node ne : Nat) (nvs nds : Nat × Nat) :
    (l.setNe neSt node ne nvs nds).edges.size = l.edges.size := by simp [Layout.setNe]
theorem setNe_adjDat_size (l : Layout) (neSt node ne : Nat) (nvs nds : Nat × Nat) :
    (l.setNe neSt node ne nvs nds).adjDat.size = l.adjDat.size := by simp [Layout.setNe]

theorem setNe_edges_get (l : Layout) (neSt node ne : Nat) (nvs nds : Nat × Nat) (j : Nat)
    (hj : j < l.edges.size) :
    (l.setNe neSt node ne nvs nds).edges[j]! =
      if j = ne - neSt then ⟨node, none, nvs⟩ else l.edges[j]! := by
  simp only [Layout.setNe]; exact getElem!_set!_eq _ _ _ _ hj

theorem setNe_adjDat_get (l : Layout) (neSt node ne : Nat) (nvs nds : Nat × Nat) (j : Nat)
    (hj : j < l.adjDat.size) :
    (l.setNe neSt node ne nvs nds).adjDat[j]! =
      if j = nds.2 - 2 * neSt then ⟨ne, nvs.1⟩
      else if j = nds.1 - 2 * neSt then ⟨ne, nvs.2⟩ else l.adjDat[j]! := by
  simp only [Layout.setNe]
  rw [getElem!_set!_eq _ _ _ _ (by simpa using hj)]
  split
  · rfl
  · exact getElem!_set!_eq _ _ _ _ hj

/-! ### V and F -/

theorem layoutNode_V (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat)) :
    layoutNode .V node nvSt nvEn neSt neEn E = Layout.empty (nvEn - nvSt) (neEn - neSt) := rfl

/-- A fold that only rewrites `adjBounds` is a fold on `adjBounds`. -/
theorem foldl_adjBounds {α : Type} (L : List α) (f : Array Nat → α → Array Nat) (l : Layout) :
    L.foldl (fun l i => { l with adjBounds := f l.adjBounds i }) l =
      { l with adjBounds := L.foldl f l.adjBounds } := by
  induction L generalizing l with
  | nil => rfl
  | cons a L ih => simp only [List.foldl_cons]; exact ih _

/-- `layoutNode .F`: all local bounds `2 neSt`, nothing else. -/
def runF (nvSt nvEn neSt neEn : Nat) : Layout :=
  { Layout.empty (nvEn - nvSt) (neEn - neSt) with
    adjBounds := (List.range' (2 * nvSt + 1) (2 * nvEn + 1 - (2 * nvSt + 1))).foldl
      (fun arr i => arr.set! (i - 2 * nvSt) (2 * neSt))
      (Layout.empty (nvEn - nvSt) (neEn - neSt)).adjBounds }

theorem layoutNode_F_eq (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat)) :
    layoutNode .F node nvSt nvEn neSt neEn E = runF nvSt nvEn neSt neEn := by
  simp only [layoutNode, Id.run, Std.Legacy.Range.forIn_eq_forIn_range', Std.Legacy.Range.size,
    bind, pure, forIn_id_yield, Nat.add_sub_cancel, Nat.div_one]
  exact foldl_adjBounds _ (fun arr i => arr.set! (i - 2 * nvSt) (2 * neSt)) _

/-- Writing `g i` at `i - off` for `i ∈ [a, a + n)`, with `off ≤ a`. -/
theorem foldl_set!_range' (arr : Array Nat) (a n off : Nat) (g : Nat → Nat) (hoff : off ≤ a)
    (j : Nat) (hj : j < arr.size) :
    ((List.range' a n).foldl (fun arr i => arr.set! (i - off) (g i)) arr)[j]! =
      if a ≤ j + off ∧ j + off < a + n then g (j + off) else arr[j]! := by
  induction n generalizing a arr with
  | zero => simp only [List.range'_zero, List.foldl_nil]; split_ifs <;> omega
  | succ n ih =>
    rw [List.range'_succ, List.foldl_cons, ih _ _ (by omega) (by simpa using hj)]
    rw [getElem!_set!_eq _ _ _ _ hj]
    split_ifs <;> first | rfl | omega | (congr 1; omega)

theorem foldl_set!_range'_size (arr : Array Nat) (a n off : Nat) (g : Nat → Nat) :
    ((List.range' a n).foldl (fun arr i => arr.set! (i - off) (g i)) arr).size = arr.size := by
  induction n generalizing a arr with
  | zero => rfl
  | succ n ih => rw [List.range'_succ, List.foldl_cons, ih]; simp

theorem runF_edges (nvSt nvEn neSt neEn : Nat) :
    (runF nvSt nvEn neSt neEn).edges = Array.replicate (neEn - neSt) default := rfl
theorem runF_adjDat (nvSt nvEn neSt neEn : Nat) :
    (runF nvSt nvEn neSt neEn).adjDat = Array.replicate (2 * (neEn - neSt)) default := rfl
theorem runF_adjBounds_size (nvSt nvEn neSt neEn : Nat) :
    (runF nvSt nvEn neSt neEn).adjBounds.size = 2 * (nvEn - nvSt) + 1 := by
  simp only [runF]
  rw [foldl_set!_range'_size (g := fun _ => 2 * neSt), empty_adjBounds_size]
theorem runF_adjBounds_get (nvSt nvEn neSt neEn : Nat) (j : Nat) (hj : j < 2 * (nvEn - nvSt) + 1) :
    (runF nvSt nvEn neSt neEn).adjBounds[j]! = if 1 ≤ j then 2 * neSt else 0 := by
  simp only [runF]
  rw [foldl_set!_range' _ _ _ _ (fun _ => 2 * neSt) (by omega) j (by rw [empty_adjBounds_size]; omega),
    empty_adjBounds_get]
  split_ifs <;> omega

theorem runF_rowBound (nvSt nvEn neSt neEn : Nat) (r : Nat) (hr1 : 2 * nvSt ≤ r) (hr2 : r ≤ 2 * nvEn) :
    rowBound nvSt neSt (runF nvSt nvEn neSt neEn) r = 2 * neSt := by
  unfold rowBound
  split
  · rfl
  · rw [runF_adjBounds_get _ _ _ _ _ (by omega)]; split_ifs <;> omega

theorem runF_row (nvSt nvEn neSt neEn : Nat) (r : Nat) (hr1 : 2 * nvSt ≤ r) (hr2 : r < 2 * nvEn) :
    row nvSt neSt (runF nvSt nvEn neSt neEn) r = [] := by
  unfold row
  rw [runF_rowBound _ _ _ _ _ hr1 (by omega), runF_rowBound _ _ _ _ _ (by omega) (by omega)]
  simp

/-! ### Single-vertex nodes (Q self-loop, O) -/

/-- `layoutNode` of a non-F/V type with one node-vertex: one loop edge `(nvSt, nvSt)`. -/
def runLoop (node nvSt nvEn neSt neEn : Nat) : Layout :=
  ({ Layout.empty (nvEn - nvSt) (neEn - neSt) with
      adjBounds := ((Layout.empty (nvEn - nvSt) (neEn - neSt)).adjBounds.set! 1 (2 * neSt + 1)).set! 2
        (2 * neSt + 2) }).setNe neSt node neSt (nvSt, nvSt) (2 * neSt + 1, 2 * neSt)

theorem layoutNode_loop_eq (type : NodeType) (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat))
    (ht : type ≠ .F ∧ type ≠ .V) (hv : nvEn - nvSt = 1) :
    layoutNode type node nvSt nvEn neSt neEn E = runLoop node nvSt nvEn neSt neEn := by
  cases type <;> simp only [ne_eq, not_true_eq_false, false_and, and_false] at ht <;>
    simp only [layoutNode, Id.run, pure, hv, beq_self_eq_true, ↓reduceIte, runLoop]

theorem runLoop_edges_size (node nvSt nvEn neSt neEn : Nat) :
    (runLoop node nvSt nvEn neSt neEn).edges.size = neEn - neSt := by
  simp [runLoop, Layout.setNe, Layout.empty]
theorem runLoop_adjBounds_size (node nvSt nvEn neSt neEn : Nat) :
    (runLoop node nvSt nvEn neSt neEn).adjBounds.size = 2 * (nvEn - nvSt) + 1 := by
  simp [runLoop, Layout.setNe, Layout.empty]
theorem runLoop_adjDat_size (node nvSt nvEn neSt neEn : Nat) :
    (runLoop node nvSt nvEn neSt neEn).adjDat.size = 2 * (neEn - neSt) := by
  simp [runLoop, Layout.setNe, Layout.empty]

theorem runLoop_edges_get (node nvSt nvEn neSt neEn : Nat) (he : neEn - neSt = 1) :
    (runLoop node nvSt nvEn neSt neEn).edges[0]! = ⟨node, none, (nvSt, nvSt)⟩ := by
  rw [runLoop, setNe_edges_get _ _ _ _ _ _ _ (by simp [Layout.empty]; omega)]
  simp

theorem runLoop_adjBounds_get (node nvSt nvEn neSt neEn : Nat) (hv : nvEn - nvSt = 1) (j : Nat)
    (hj : j ≤ 2) :
    (runLoop node nvSt nvEn neSt neEn).adjBounds[j]! =
      if j = 0 then 0 else if j = 1 then 2 * neSt + 1 else 2 * neSt + 2 := by
  rw [runLoop, setNe_adjBounds]
  simp only
  rw [getElem!_set!_eq _ _ _ _ (by simp [Layout.empty]; omega),
    getElem!_set!_eq _ _ _ _ (by simp [Layout.empty]; omega), empty_adjBounds_get]
  split_ifs <;> omega

theorem runLoop_adjDat_get (node nvSt nvEn neSt neEn : Nat) (he : neEn - neSt = 1) (j : Nat)
    (hj : j < 2) :
    (runLoop node nvSt nvEn neSt neEn).adjDat[j]! = ⟨neSt, nvSt⟩ := by
  rw [runLoop, setNe_adjDat_get _ _ _ _ _ _ _ (by simp [Layout.empty]; omega)]
  simp only [Nat.add_sub_cancel_left, Nat.sub_self]
  split_ifs <;> first | rfl | omega

theorem runLoop_rowBound (node nvSt nvEn neSt neEn : Nat) (hv : nvEn - nvSt = 1) (r : Nat)
    (hr1 : 2 * nvSt ≤ r) (hr2 : r ≤ 2 * nvEn) :
    rowBound nvSt neSt (runLoop node nvSt nvEn neSt neEn) r = 2 * neSt + (r - 2 * nvSt) := by
  unfold rowBound
  split
  · omega
  · rw [runLoop_adjBounds_get _ _ _ _ _ hv _ (by omega)]; split_ifs <;> omega

theorem runLoop_row (node nvSt nvEn neSt neEn : Nat) (hv : nvEn - nvSt = 1) (he : neEn - neSt = 1)
    (r : Nat) (hr1 : 2 * nvSt ≤ r) (hr2 : r < 2 * nvEn) :
    row nvSt neSt (runLoop node nvSt nvEn neSt neEn) r = [⟨neSt, nvSt⟩] := by
  unfold row
  rw [runLoop_rowBound _ _ _ _ _ hv _ hr1 (by omega), runLoop_rowBound _ _ _ _ _ hv _ (by omega) (by omega),
    show 2 * neSt + (r + 1 - 2 * nvSt) - (2 * neSt + (r - 2 * nvSt)) = 1 by omega]
  simp only [List.range_one, List.map_cons, List.map_nil, Nat.add_zero, List.cons.injEq, and_true]
  rw [runLoop_adjDat_get _ _ _ _ _ he _ (by omega)]

/-! ### Q / I with two node-vertices -/

/-- `layoutNode .Q`/`.I` with two node-vertices: one edge `(nvSt, nvSt + 1)`. -/
def runQI (node nvSt nvEn neSt neEn : Nat) : Layout :=
  ({ Layout.empty (nvEn - nvSt) (neEn - neSt) with
      adjBounds := (((((Layout.empty (nvEn - nvSt) (neEn - neSt)).adjBounds.set! 1 (2 * neSt)).set! 2
        (2 * neSt + 1)).set! 3 (2 * neSt + 2)).set! 4 (2 * neSt + 2)) }).setNe neSt node neSt
    (nvSt, nvSt + 1) (2 * neSt, 2 * neSt + 1)

theorem layoutNode_QI_eq (type : NodeType) (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat))
    (ht : type = .Q ∨ type = .I) (hv : nvEn - nvSt ≠ 1) :
    layoutNode type node nvSt nvEn neSt neEn E = runQI node nvSt nvEn neSt neEn := by
  have hv' : (nvEn - nvSt == 1) = false := by simpa using hv
  rcases ht with rfl | rfl <;>
    simp only [layoutNode, Id.run, pure, hv', Bool.false_eq_true, ↓reduceIte,
      beq_self_eq_true, Bool.true_or, Bool.or_true] <;> rfl

theorem runQI_edges_size (node nvSt nvEn neSt neEn : Nat) :
    (runQI node nvSt nvEn neSt neEn).edges.size = neEn - neSt := by
  simp [runQI, Layout.setNe, Layout.empty]
theorem runQI_adjBounds_size (node nvSt nvEn neSt neEn : Nat) :
    (runQI node nvSt nvEn neSt neEn).adjBounds.size = 2 * (nvEn - nvSt) + 1 := by
  simp [runQI, Layout.setNe, Layout.empty]
theorem runQI_adjDat_size (node nvSt nvEn neSt neEn : Nat) :
    (runQI node nvSt nvEn neSt neEn).adjDat.size = 2 * (neEn - neSt) := by
  simp [runQI, Layout.setNe, Layout.empty]

theorem runQI_edges_get (node nvSt nvEn neSt neEn : Nat) (he : neEn - neSt = 1) :
    (runQI node nvSt nvEn neSt neEn).edges[0]! = ⟨node, none, (nvSt, nvSt + 1)⟩ := by
  rw [runQI, setNe_edges_get _ _ _ _ _ _ _ (by simp [Layout.empty]; omega)]
  simp

theorem runQI_adjBounds_get (node nvSt nvEn neSt neEn : Nat) (hv : nvEn - nvSt = 2) (j : Nat)
    (hj : j ≤ 4) :
    (runQI node nvSt nvEn neSt neEn).adjBounds[j]! =
      if j = 0 then 0 else if j = 1 then 2 * neSt else if j = 2 then 2 * neSt + 1
      else 2 * neSt + 2 := by
  rw [runQI, setNe_adjBounds]
  simp only
  have hs : (Layout.empty (nvEn - nvSt) (neEn - neSt)).adjBounds.size = 5 := by
    simp [Layout.empty]; omega
  rw [getElem!_set!_eq _ _ _ _ (by simp [Layout.empty]; omega),
    getElem!_set!_eq _ _ _ _ (by simp [Layout.empty]; omega),
    getElem!_set!_eq _ _ _ _ (by simp [Layout.empty]; omega),
    getElem!_set!_eq _ _ _ _ (by simp [Layout.empty]; omega), empty_adjBounds_get]
  split_ifs <;> omega

theorem runQI_adjDat_get (node nvSt nvEn neSt neEn : Nat) (he : neEn - neSt = 1) (j : Nat)
    (hj : j < 2) :
    (runQI node nvSt nvEn neSt neEn).adjDat[j]! =
      if j = 0 then ⟨neSt, nvSt + 1⟩ else ⟨neSt, nvSt⟩ := by
  rw [runQI, setNe_adjDat_get _ _ _ _ _ _ _ (by simp [Layout.empty]; omega)]
  simp only [Nat.add_sub_cancel_left, Nat.sub_self]
  split_ifs <;> first | rfl | omega

theorem runQI_rowBound (node nvSt nvEn neSt neEn : Nat) (hv : nvEn - nvSt = 2) (r : Nat)
    (hr1 : 2 * nvSt ≤ r) (hr2 : r ≤ 2 * nvEn) :
    rowBound nvSt neSt (runQI node nvSt nvEn neSt neEn) r =
      if r ≤ 2 * nvSt + 1 then 2 * neSt else if r = 2 * nvSt + 2 then 2 * neSt + 1
      else 2 * neSt + 2 := by
  unfold rowBound
  split
  · split_ifs <;> omega
  · rw [runQI_adjBounds_get _ _ _ _ _ hv _ (by omega)]; split_ifs <;> omega

theorem runQI_row (node nvSt nvEn neSt neEn : Nat) (hv : nvEn - nvSt = 2) (he : neEn - neSt = 1)
    (r : Nat) (hr1 : 2 * nvSt ≤ r) (hr2 : r < 2 * nvEn) :
    row nvSt neSt (runQI node nvSt nvEn neSt neEn) r =
      if r = 2 * nvSt + 1 then [⟨neSt, nvSt + 1⟩]
      else if r = 2 * nvSt + 2 then [⟨neSt, nvSt⟩] else [] := by
  unfold row
  rw [runQI_rowBound _ _ _ _ _ hv _ hr1 (by omega), runQI_rowBound _ _ _ _ _ hv _ (by omega) (by omega)]
  by_cases h1 : r = 2 * nvSt + 1
  · subst h1
    simp only [Nat.le_refl, ↓reduceIte, Nat.add_sub_cancel_left,
      List.range_one, List.map_cons, List.map_nil, List.cons.injEq, and_true,
      show ¬ (2 * nvSt + 1 + 1 ≤ 2 * nvSt + 1) by omega, show 2 * nvSt + 1 + 1 = 2 * nvSt + 2 by omega]
    rw [runQI_adjDat_get _ _ _ _ _ he _ (by omega)]; simp
  by_cases h2 : r = 2 * nvSt + 2
  · subst h2
    simp only [show ¬ (2 * nvSt + 2 ≤ 2 * nvSt + 1) by omega, show ¬ (2 * nvSt + 2 + 1 ≤ 2 * nvSt + 1) by omega,
      show 2 * nvSt + 2 + 1 ≠ 2 * nvSt + 2 by omega, ↓reduceIte, h1,
      show 2 * neSt + 2 - (2 * neSt + 1) = 1 by omega, List.range_one, List.map_cons, List.map_nil,
      Nat.add_zero, List.cons.injEq, and_true]
    rw [runQI_adjDat_get _ _ _ _ _ he _ (by omega)]; simp
  · simp only [h1, h2, ↓reduceIte]
    split_ifs <;> simp <;> omega

/-! ### P -/

/-- Initial P layout: bounds `[_, 2 neSt, 2 neSt + nEdges, 2 neSt + 2 nEdges, 2 neSt + 2 nEdges]`. -/
def initP (nvSt nvEn neSt neEn : Nat) : Layout :=
  { Layout.empty (nvEn - nvSt) (neEn - neSt) with
    adjBounds := (((((Layout.empty (nvEn - nvSt) (neEn - neSt)).adjBounds.set! 1 (2 * neSt)).set! 2
      (2 * neSt + (neEn - neSt))).set! 3 (2 * neSt + 2 * (neEn - neSt))).set! 4
      (2 * neSt + 2 * (neEn - neSt))) }

/-- Step `k` of the P fill: edge `neSt + k` at slots `k` (row `2 nvSt + 1`) and `2 nEdges - 1 - k`
(row `2 nvSt + 2`, reversed). -/
def stepP (node nvSt neSt neEn : Nat) (l : Layout) (k : Nat) : Layout :=
  l.setNe neSt node (neSt + k) (nvSt, nvSt + 1) (2 * neSt + k, 2 * neEn - 1 - k)

def runP (node nvSt nvEn neSt neEn : Nat) : Layout :=
  (List.range' 0 (neEn - neSt)).foldl (stepP node nvSt neSt neEn) (initP nvSt nvEn neSt neEn)

theorem layoutNode_P_eq (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat))
    (hv : nvEn - nvSt ≠ 1) :
    layoutNode .P node nvSt nvEn neSt neEn E = runP node nvSt nvEn neSt neEn := by
  have hv' : (nvEn - nvSt == 1) = false := by simpa using hv
  simp only [layoutNode, Id.run, pure, hv', Bool.false_eq_true, ↓reduceIte, beq_self_eq_true,
    Std.Legacy.Range.forIn_eq_forIn_range', Std.Legacy.Range.size, forIn_id_yield,
    Nat.add_sub_cancel, Nat.div_one, Nat.sub_zero]
  rfl

theorem initP_adjBounds_size (nvSt nvEn neSt neEn : Nat) :
    (initP nvSt nvEn neSt neEn).adjBounds.size = 2 * (nvEn - nvSt) + 1 := by
  simp [initP, Layout.empty]

theorem initP_adjBounds_get (nvSt nvEn neSt neEn : Nat) (hv : nvEn - nvSt = 2) (j : Nat) (hj : j ≤ 4) :
    (initP nvSt nvEn neSt neEn).adjBounds[j]! =
      if j = 0 then 0 else if j = 1 then 2 * neSt else if j = 2 then 2 * neSt + (neEn - neSt)
      else 2 * neSt + 2 * (neEn - neSt) := by
  simp only [initP]
  rw [getElem!_set!_eq _ _ _ _ (by simp [Layout.empty]; omega),
    getElem!_set!_eq _ _ _ _ (by simp [Layout.empty]; omega),
    getElem!_set!_eq _ _ _ _ (by simp [Layout.empty]; omega),
    getElem!_set!_eq _ _ _ _ (by simp [Layout.empty]; omega), empty_adjBounds_get]
  split_ifs <;> omega

/-- The P fill after `m` steps: edges `0 … m - 1` written, slots `k < m` hold `⟨neSt + k, nvSt + 1⟩`,
slots `2 nE - 1 - k` hold `⟨neSt + k, nvSt⟩`, everything else untouched. -/
theorem foldl_stepP (node nvSt neSt neEn : Nat) (l : Layout) (nE : Nat) (hnE : neEn = neSt + nE)
    (hs : l.edges.size = nE) (hd : l.adjDat.size = 2 * nE) (m : Nat) (hm : m ≤ nE) :
    let l' := (List.range' 0 m).foldl (stepP node nvSt neSt neEn) l
    l'.adjBounds = l.adjBounds ∧ l'.edges.size = nE ∧ l'.adjDat.size = 2 * nE ∧
    (∀ j, j < nE → l'.edges[j]! = if j < m then ⟨node, none, (nvSt, nvSt + 1)⟩ else l.edges[j]!) ∧
    (∀ j, j < 2 * nE → l'.adjDat[j]! =
      if j < m then ⟨neSt + j, nvSt + 1⟩
      else if 2 * nE - m ≤ j then ⟨neSt + (2 * nE - 1 - j), nvSt⟩ else l.adjDat[j]!) := by
  intro l'
  induction m with
  | zero =>
    simp only [l', List.range'_zero, List.foldl_nil]
    refine ⟨trivial, hs, hd, ?_, ?_⟩
    · intro j _; simp
    · intro j hj
      rw [ite_of_neg (by omega), ite_of_neg (by omega)]
  | succ m ih =>
    obtain ⟨b, se, sd, ge, gd⟩ := ih (by omega)
    simp only [l', List.range'_1_concat, List.foldl_append, List.foldl_cons, List.foldl_nil, Nat.zero_add]
    refine ⟨b, ?_, ?_, ?_, ?_⟩
    · rw [stepP, setNe_edges_size, se]
    · rw [stepP, setNe_adjDat_size, sd]
    · intro j hj
      rw [stepP, setNe_edges_get _ _ _ _ _ _ _ (by omega), ge j hj, Nat.add_sub_cancel_left]
      split_ifs <;> first | rfl | omega
    · intro j hj
      rw [stepP, setNe_adjDat_get _ _ _ _ _ _ _ (by omega), gd j hj]
      simp only [Nat.add_sub_cancel_left, hnE]
      split_ifs <;> first | rfl | omega | (congr 1; omega)

theorem runP_spec (node nvSt nvEn neSt neEn : Nat) (he : neSt ≤ neEn) :
    let l' := runP node nvSt nvEn neSt neEn
    l'.adjBounds = (initP nvSt nvEn neSt neEn).adjBounds ∧ l'.edges.size = neEn - neSt ∧
    l'.adjDat.size = 2 * (neEn - neSt) ∧
    (∀ j, j < neEn - neSt → l'.edges[j]! = ⟨node, none, (nvSt, nvSt + 1)⟩) ∧
    (∀ j, j < 2 * (neEn - neSt) → l'.adjDat[j]! =
      if j < neEn - neSt then ⟨neSt + j, nvSt + 1⟩ else ⟨neEn - 1 - (j - (neEn - neSt)), nvSt⟩) := by
  intro l'
  obtain ⟨b, se, sd, ge, gd⟩ := foldl_stepP node nvSt neSt neEn (initP nvSt nvEn neSt neEn)
    (neEn - neSt) (by omega) (by simp [initP, Layout.empty]) (by simp [initP, Layout.empty]) _
    (Nat.le_refl _)
  simp only [l', runP]
  refine ⟨b, se, sd, ?_, ?_⟩
  · intro j hj; rw [ge j hj, ite_of_pos hj]
  · intro j hj
    rw [gd j hj]
    split_ifs <;> first | rfl | omega | (congr 1; omega)

theorem runP_rowBound (node nvSt nvEn neSt neEn : Nat) (hv : nvEn - nvSt = 2) (he : neSt ≤ neEn)
    (r : Nat) (hr1 : 2 * nvSt ≤ r) (hr2 : r ≤ 2 * nvEn) :
    rowBound nvSt neSt (runP node nvSt nvEn neSt neEn) r =
      if r ≤ 2 * nvSt + 1 then 2 * neSt else if r = 2 * nvSt + 2 then 2 * neSt + (neEn - neSt)
      else 2 * neSt + 2 * (neEn - neSt) := by
  unfold rowBound
  split
  · split_ifs <;> omega
  · rw [(runP_spec node nvSt nvEn neSt neEn he).1, initP_adjBounds_get _ _ _ _ hv _ (by omega)]
    split_ifs <;> omega

/-- Rows of a P node: `2 nvSt + 1` lists the parallel edges in increasing order (towards
`nvSt + 1`), `2 nvSt + 2` lists them in decreasing order (towards `nvSt`); the other two are
empty. -/
theorem runP_row (node nvSt nvEn neSt neEn : Nat) (hv : nvEn - nvSt = 2) (he : neSt ≤ neEn) (r : Nat)
    (hr1 : 2 * nvSt ≤ r) (hr2 : r < 2 * nvEn) :
    row nvSt neSt (runP node nvSt nvEn neSt neEn) r =
      if r = 2 * nvSt + 1 then (List.range (neEn - neSt)).map fun k => ⟨neSt + k, nvSt + 1⟩
      else if r = 2 * nvSt + 2 then (List.range (neEn - neSt)).map fun k => ⟨neEn - 1 - k, nvSt⟩
      else [] := by
  obtain ⟨-, -, -, -, gd⟩ := runP_spec node nvSt nvEn neSt neEn he
  have hB := runP_rowBound node nvSt nvEn neSt neEn hv he
  by_cases h1 : r = 2 * nvSt + 1
  · subst h1
    have b0 := hB (2 * nvSt + 1) (by omega) (by omega)
    rw [ite_of_pos (by omega)] at b0
    have b1 := hB (2 * nvSt + 1 + 1) (by omega) (by omega)
    rw [ite_of_neg (by omega), ite_of_pos (by omega)] at b1
    unfold row
    rw [b0, b1, ite_of_pos rfl, Nat.add_sub_cancel_left]
    apply List.map_congr_left
    intro k hk
    rw [List.mem_range] at hk
    rw [show 2 * neSt + k - 2 * neSt = k by omega, gd _ (by omega), ite_of_pos hk]
  by_cases h2 : r = 2 * nvSt + 2
  · subst h2
    have b0 := hB (2 * nvSt + 2) (by omega) (by omega)
    rw [ite_of_neg (by omega), ite_of_pos rfl] at b0
    have b1 := hB (2 * nvSt + 2 + 1) (by omega) (by omega)
    rw [ite_of_neg (by omega), ite_of_neg (by omega)] at b1
    unfold row
    rw [b0, b1, ite_of_neg h1, ite_of_pos rfl,
      show 2 * neSt + 2 * (neEn - neSt) - (2 * neSt + (neEn - neSt)) = neEn - neSt by omega]
    apply List.map_congr_left
    intro k hk
    rw [List.mem_range] at hk
    rw [show 2 * neSt + (neEn - neSt) + k - 2 * neSt = neEn - neSt + k by omega,
      gd _ (by omega), ite_of_neg (by omega), Nat.add_sub_cancel_left]
  · rw [ite_of_neg h1, ite_of_neg h2]
    have b0 := hB r hr1 (by omega)
    have b1 := hB (r + 1) (by omega) (by omega)
    rcases (show r = 2 * nvSt ∨ r = 2 * nvSt + 3 by omega) with rfl | rfl
    · rw [ite_of_pos (by omega)] at b0; rw [ite_of_pos (by omega)] at b1
      unfold row; rw [b0, b1]; simp
    · rw [ite_of_neg (by omega), ite_of_neg (by omega)] at b0
      rw [ite_of_neg (by omega), ite_of_neg (by omega)] at b1
      unfold row; rw [b0, b1]; simp

/-! ### S -/

/-- Closes `first | rfl | omega` goals between `NodeEdge`/`NodeAdj` records. -/
macro "rec_omega" : tactic =>
  `(tactic| first | rfl | omega |
      (simp only [NodeEdge.mk.injEq, NodeAdj.mk.injEq, Prod.mk.injEq, List.cons.injEq, true_and,
        and_true] <;> omega))

/-- S bounds before the cap: local bound `j ↦ j + 2 neSt` for `1 ≤ j ≤ 2 n`, then row `2 nvSt`
loses its slot to row `2 nvSt + 1` (`modify 1 (· - 1)`) and row `2 nvEn - 1` gives its slot to
row `2 nvEn - 2` (`modify (2 nvEn - 1 - 2 nvSt) (· + 1)`). -/
def boundsS (nvSt nvEn neSt neEn : Nat) : Array Nat :=
  (((List.range' (2 * nvSt + 1) (2 * nvEn + 1 - (2 * nvSt + 1))).foldl
    (fun arr i => arr.set! (i - 2 * nvSt) (i - 2 * nvSt + 2 * neSt))
    (Layout.empty (nvEn - nvSt) (neEn - neSt)).adjBounds).modify 1 (· - 1)).modify
      (2 * nvEn - 1 - 2 * nvSt) (· + 1)

/-- Bounds plus the cap `(nvSt, nvEn - 1)` at `neSt`, slots `0` and `2 nEdges - 1`. -/
def initS (node nvSt nvEn neSt neEn : Nat) : Layout :=
  ({ Layout.empty (nvEn - nvSt) (neEn - neSt) with adjBounds := boundsS nvSt nvEn neSt neEn }).setNe
    neSt node neSt (nvSt, nvEn - 1) (2 * neSt, 2 * neEn - 1)

/-- Path edge `i`: `neSt + i = (nvSt + i - 1, nvSt + i)` at slots `2 i - 1` and `2 i`. -/
def stepS (node nvSt neSt : Nat) (l : Layout) (i : Nat) : Layout :=
  l.setNe neSt node (neSt + i) (nvSt + i - 1, nvSt + i) (2 * (neSt + i) - 1, 2 * (neSt + i))

def runS (node nvSt nvEn neSt neEn : Nat) : Layout :=
  (List.range' 1 (neEn - neSt - 1)).foldl (stepS node nvSt neSt) (initS node nvSt nvEn neSt neEn)

theorem layoutNode_S_eq (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat))
    (hv : nvEn - nvSt ≠ 1) :
    layoutNode .S node nvSt nvEn neSt neEn E = runS node nvSt nvEn neSt neEn := by
  have hv' : (nvEn - nvSt == 1) = false := by simpa using hv
  simp only [layoutNode, Id.run, pure, hv', Bool.false_eq_true, ↓reduceIte, beq_self_eq_true,
    Std.Legacy.Range.forIn_eq_forIn_range', Std.Legacy.Range.size, forIn_id_yield,
    Nat.add_sub_cancel, Nat.div_one]
  rw [foldl_adjBounds _ (fun arr i => arr.set! (i - 2 * nvSt) (i - 2 * nvSt + 2 * neSt))]
  rfl

theorem boundsS_size (nvSt nvEn neSt neEn : Nat) :
    (boundsS nvSt nvEn neSt neEn).size = 2 * (nvEn - nvSt) + 1 := by
  simp only [boundsS, Array.size_modify]
  rw [foldl_set!_range'_size (g := fun i => i - 2 * nvSt + 2 * neSt), empty_adjBounds_size]

theorem boundsS_get (nvSt nvEn neSt neEn : Nat) (hv : 2 ≤ nvEn - nvSt) (j : Nat)
    (hj : j ≤ 2 * (nvEn - nvSt)) :
    (boundsS nvSt nvEn neSt neEn)[j]! =
      if j = 0 then 0 else if j = 1 then 2 * neSt
      else if j = 2 * (nvEn - nvSt) - 1 then 2 * neSt + 2 * (nvEn - nvSt) else 2 * neSt + j := by
  have hsz : ((List.range' (2 * nvSt + 1) (2 * nvEn + 1 - (2 * nvSt + 1))).foldl
      (fun arr i => arr.set! (i - 2 * nvSt) (i - 2 * nvSt + 2 * neSt))
      (Layout.empty (nvEn - nvSt) (neEn - neSt)).adjBounds).size = 2 * (nvEn - nvSt) + 1 := by
    rw [foldl_set!_range'_size (g := fun i => i - 2 * nvSt + 2 * neSt), empty_adjBounds_size]
  simp only [boundsS]
  rw [getElem!_modify_eq _ _ _ _ (by simp only [Array.size_modify]; omega),
    getElem!_modify_eq _ _ _ _ (by omega),
    foldl_set!_range' _ _ _ _ (fun i => i - 2 * nvSt + 2 * neSt) (by omega) j
      (by rw [empty_adjBounds_size]; omega), empty_adjBounds_get]
  split_ifs <;> omega

theorem initS_adjBounds (node nvSt nvEn neSt neEn : Nat) :
    (initS node nvSt nvEn neSt neEn).adjBounds = boundsS nvSt nvEn neSt neEn := rfl
theorem initS_edges_size (node nvSt nvEn neSt neEn : Nat) :
    (initS node nvSt nvEn neSt neEn).edges.size = neEn - neSt := by
  simp [initS, Layout.setNe, Layout.empty]
theorem initS_adjDat_size (node nvSt nvEn neSt neEn : Nat) :
    (initS node nvSt nvEn neSt neEn).adjDat.size = 2 * (neEn - neSt) := by
  simp [initS, Layout.setNe, Layout.empty]

theorem initS_edges_get (node nvSt nvEn neSt neEn : Nat) (j : Nat) (hj : j < neEn - neSt) :
    (initS node nvSt nvEn neSt neEn).edges[j]! =
      if j = 0 then ⟨node, none, (nvSt, nvEn - 1)⟩ else default := by
  rw [initS, setNe_edges_get _ _ _ _ _ _ _ (by simp [Layout.empty]; omega), Nat.sub_self]
  simp only [empty_edges_get]

theorem initS_adjDat_get (node nvSt nvEn neSt neEn : Nat) (j : Nat) (hj : j < 2 * (neEn - neSt)) :
    (initS node nvSt nvEn neSt neEn).adjDat[j]! =
      if j = 2 * (neEn - neSt) - 1 then ⟨neSt, nvSt⟩
      else if j = 0 then ⟨neSt, nvEn - 1⟩ else default := by
  rw [initS, setNe_adjDat_get _ _ _ _ _ _ _ (by simp [Layout.empty]; omega)]
  simp only [empty_adjDat_get, Nat.sub_self]
  split_ifs <;> rec_omega

/-- The S path fill after steps `1 … m`. -/
theorem foldl_stepS (node nvSt neSt : Nat) (l : Layout) (nE : Nat)
    (hs : l.edges.size = nE) (hd : l.adjDat.size = 2 * nE) (m : Nat) (hm : m + 1 ≤ nE) :
    let l' := (List.range' 1 m).foldl (stepS node nvSt neSt) l
    l'.adjBounds = l.adjBounds ∧ l'.edges.size = nE ∧ l'.adjDat.size = 2 * nE ∧
    (∀ j, j < nE → l'.edges[j]! =
      if 1 ≤ j ∧ j ≤ m then ⟨node, none, (nvSt + j - 1, nvSt + j)⟩ else l.edges[j]!) ∧
    (∀ j, j < 2 * nE → l'.adjDat[j]! =
      if 1 ≤ j ∧ j ≤ 2 * m then
        (if j % 2 = 1 then ⟨neSt + (j + 1) / 2, nvSt + (j + 1) / 2⟩
          else ⟨neSt + j / 2, nvSt + j / 2 - 1⟩)
      else l.adjDat[j]!) := by
  intro l'
  induction m with
  | zero =>
    simp only [l', List.range'_zero, List.foldl_nil]
    refine ⟨trivial, hs, hd, ?_, ?_⟩
    · intro j _; rw [ite_of_neg (by omega)]
    · intro j _; rw [ite_of_neg (by omega)]
  | succ m ih =>
    obtain ⟨b, se, sd, ge, gd⟩ := ih (by omega)
    simp only [l', List.range'_1_concat, List.foldl_append, List.foldl_cons, List.foldl_nil]
    refine ⟨b, ?_, ?_, ?_, ?_⟩
    · rw [stepS, setNe_edges_size, se]
    · rw [stepS, setNe_adjDat_size, sd]
    · intro j hj
      rw [stepS, setNe_edges_get _ _ _ _ _ _ _ (by omega), ge j hj, Nat.add_sub_cancel_left]
      split_ifs <;> rec_omega
    · intro j hj
      rw [stepS, setNe_adjDat_get _ _ _ _ _ _ _ (by omega), gd j hj]
      simp only
      split_ifs <;> rec_omega

theorem runS_spec (node nvSt nvEn neSt neEn : Nat) (he : neSt + 1 ≤ neEn) :
    let l' := runS node nvSt nvEn neSt neEn
    l'.adjBounds = boundsS nvSt nvEn neSt neEn ∧ l'.edges.size = neEn - neSt ∧
    l'.adjDat.size = 2 * (neEn - neSt) ∧
    (∀ j, j < neEn - neSt → l'.edges[j]! =
      if j = 0 then ⟨node, none, (nvSt, nvEn - 1)⟩ else ⟨node, none, (nvSt + j - 1, nvSt + j)⟩) ∧
    (∀ j, j < 2 * (neEn - neSt) → l'.adjDat[j]! =
      if j = 0 then ⟨neSt, nvEn - 1⟩
      else if j = 2 * (neEn - neSt) - 1 then ⟨neSt, nvSt⟩
      else if j % 2 = 1 then ⟨neSt + (j + 1) / 2, nvSt + (j + 1) / 2⟩
      else ⟨neSt + j / 2, nvSt + j / 2 - 1⟩) := by
  intro l'
  obtain ⟨b, se, sd, ge, gd⟩ := foldl_stepS node nvSt neSt (initS node nvSt nvEn neSt neEn)
    (neEn - neSt) (initS_edges_size ..) (initS_adjDat_size ..) (neEn - neSt - 1) (by omega)
  simp only [l', runS]
  refine ⟨b, se, sd, ?_, ?_⟩
  · intro j hj
    rw [ge j hj, initS_edges_get _ _ _ _ _ _ hj]
    split_ifs <;> rec_omega
  · intro j hj
    rw [gd j hj, initS_adjDat_get _ _ _ _ _ _ hj]
    split_ifs <;> rec_omega

theorem runS_rowBound (node nvSt nvEn neSt neEn : Nat) (hv : 2 ≤ nvEn - nvSt) (he : neSt + 1 ≤ neEn)
    (r : Nat) (hr1 : 2 * nvSt ≤ r) (hr2 : r ≤ 2 * nvEn) :
    rowBound nvSt neSt (runS node nvSt nvEn neSt neEn) r =
      if r ≤ 2 * nvSt + 1 then 2 * neSt
      else if 2 * nvEn - 1 ≤ r then 2 * neSt + 2 * (nvEn - nvSt) else 2 * neSt + (r - 2 * nvSt) := by
  unfold rowBound
  split
  · split_ifs <;> omega
  · rw [(runS_spec node nvSt nvEn neSt neEn he).1, boundsS_get _ _ _ _ hv _ (by omega)]
    split_ifs <;> omega

/-- Rows of an S node (cycle `nvSt, …, nvEn - 1` closed by the cap at `neSt`): the ends `nvSt` /
`nvEn - 1` have two neighbours in one row (cap first for `nvSt`, cap last for `nvEn - 1`), the
interior vertices one neighbour per row. -/
theorem runS_row (node nvSt nvEn neSt neEn : Nat) (hv : 3 ≤ nvEn - nvSt) (he : neEn - neSt = nvEn - nvSt)
    (r : Nat) (hr1 : 2 * nvSt ≤ r) (hr2 : r < 2 * nvEn) :
    row nvSt neSt (runS node nvSt nvEn neSt neEn) r =
      if r = 2 * nvSt then []
      else if r = 2 * nvSt + 1 then [⟨neSt, nvEn - 1⟩, ⟨neSt + 1, nvSt + 1⟩]
      else if r = 2 * nvEn - 2 then [⟨neEn - 1, nvEn - 2⟩, ⟨neSt, nvSt⟩]
      else if r = 2 * nvEn - 1 then []
      else if r % 2 = 0 then [⟨neSt + (r / 2 - nvSt), r / 2 - 1⟩]
      else [⟨neSt + (r / 2 - nvSt) + 1, r / 2 + 1⟩] := by
  have he' : neSt + 1 ≤ neEn := by omega
  obtain ⟨-, -, -, -, gd⟩ := runS_spec node nvSt nvEn neSt neEn he'
  have hB := runS_rowBound node nvSt nvEn neSt neEn (by omega) he'
  have b0 := hB r hr1 (by omega)
  have b1 := hB (r + 1) (by omega) (by omega)
  have h2 : List.range 2 = [0, 1] := rfl
  unfold row
  by_cases c0 : r = 2 * nvSt
  · rw [ite_of_pos c0]; subst c0
    rw [ite_of_pos (by omega)] at b0; rw [ite_of_pos (by omega)] at b1
    rw [b0, b1]; simp
  by_cases c1 : r = 2 * nvSt + 1
  · rw [ite_of_neg c0, ite_of_pos c1]; subst c1
    rw [ite_of_pos (by omega)] at b0; rw [ite_of_neg (by omega), ite_of_neg (by omega)] at b1
    rw [b0, b1, show 2 * neSt + (2 * nvSt + 1 + 1 - 2 * nvSt) - 2 * neSt = 2 by omega, h2]
    simp only [List.map_cons, List.map_nil, Nat.add_sub_cancel_left]
    rw [gd 0 (by omega), gd 1 (by omega), ite_of_neg (show ¬ (1 = 2 * (neEn - neSt) - 1) by omega)]
    simp
  by_cases c2 : r = 2 * nvEn - 2
  · rw [ite_of_neg c0, ite_of_neg c1, ite_of_pos c2]; subst c2
    rw [ite_of_neg (by omega), ite_of_neg (by omega)] at b0
    rw [ite_of_neg (by omega), ite_of_pos (by omega)] at b1
    rw [b0, b1, show 2 * neSt + 2 * (nvEn - nvSt) - (2 * neSt + (2 * nvEn - 2 - 2 * nvSt)) = 2 by omega,
      h2]
    simp only [List.map_cons, List.map_nil]
    rw [gd _ (by omega), gd _ (by omega)]
    split_ifs <;> rec_omega
  by_cases c3 : r = 2 * nvEn - 1
  · rw [ite_of_neg c0, ite_of_neg c1, ite_of_neg c2, ite_of_pos c3]; subst c3
    rw [ite_of_neg (by omega), ite_of_pos (by omega)] at b0
    rw [ite_of_neg (by omega), ite_of_pos (by omega)] at b1
    rw [b0, b1]; try simp
  rw [ite_of_neg c0, ite_of_neg c1, ite_of_neg c2, ite_of_neg c3]
  rw [ite_of_neg (by omega), ite_of_neg (by omega)] at b0
  rw [ite_of_neg (by omega), ite_of_neg (by omega)] at b1
  rw [b0, b1, show 2 * neSt + (r + 1 - 2 * nvSt) - (2 * neSt + (r - 2 * nvSt)) = 1 by omega]
  simp only [List.range_one, List.map_cons, List.map_nil, Nat.add_zero]
  rw [gd _ (by omega)]
  split_ifs <;> rec_omega

/-! ### R -/

section R
open LayoutR

theorem fillStep_edges_size (nvSt neSt node : Nat) (s : Layout × Nat) (p : Nat × Nat) :
    (fillStep nvSt neSt node s p).1.edges.size = s.1.edges.size := by
  simp [fillStep, Layout.setNe, countStep, inc]

theorem fillStep_snd (nvSt neSt node : Nat) (s : Layout × Nat) (p : Nat × Nat) :
    (fillStep nvSt neSt node s p).2 = s.2 - 1 := rfl

/-- The reverse fill writes `L[i]` as node-edge `c - 1 - i`. -/
theorem foldl_fillStep_edges (nvSt neSt node : Nat) (L : List (Nat × Nat)) :
    ∀ (l : Layout) (c : Nat), neSt + L.length ≤ c → c ≤ neSt + l.edges.size →
      (L.foldl (fillStep nvSt neSt node) (l, c)).1.edges.size = l.edges.size ∧
      ∀ j, j < l.edges.size →
        (L.foldl (fillStep nvSt neSt node) (l, c)).1.edges[j]! =
          if c ≤ j + neSt + L.length ∧ j + neSt < c then
            ⟨node, none, L.getD (c - 1 - neSt - j) default⟩
          else l.edges[j]! := by
  induction L with
  | nil =>
    intro l c _ _
    refine ⟨rfl, fun j _ => ?_⟩
    simp only [List.length_nil, List.foldl_nil]
    rw [ite_of_neg (by omega)]
  | cons p L ih =>
    intro l c h1 h2
    rw [List.foldl_cons]
    have hsz : (fillStep nvSt neSt node (l, c) p).1.edges.size = l.edges.size :=
      fillStep_edges_size ..
    simp only [List.length_cons] at h1
    obtain ⟨sz, g⟩ := ih (fillStep nvSt neSt node (l, c) p).1 (fillStep nvSt neSt node (l, c) p).2
      (by rw [fillStep_snd]; omega) (by rw [fillStep_snd, hsz]; omega)
    refine ⟨sz.trans hsz, fun j hj => ?_⟩
    rw [g j (by omega), fillStep_snd]
    simp only [fillStep, List.length_cons]
    rw [setNe_edges_get _ _ _ _ _ _ _ (by rw [countStep_edges]; omega), countStep_edges]
    by_cases hc : c - 1 ≤ j + neSt + L.length ∧ j + neSt < c - 1
    · rw [ite_of_pos hc, ite_of_pos (by omega)]
      obtain ⟨k, hk⟩ : ∃ k, c - 1 - neSt - j = k + 1 := ⟨c - 1 - 1 - neSt - j, by omega⟩
      rw [hk, List.getD_cons_succ]
      congr 3; omega
    · rw [ite_of_neg hc]
      by_cases hj0 : j = c - 1 - neSt
      · rw [ite_of_pos hj0, ite_of_pos (by omega), show c - 1 - neSt - j = 0 by omega,
          List.getD_cons_zero]
      · rw [ite_of_neg hj0, ite_of_neg (by omega)]

/-- Edges of an R layout: the cap `(nvSt, nvEn - 1)` at local index `0`, then `E` in order. -/
theorem run_edges (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat)) (hv : nvSt + 2 ≤ nvEn)
    (hE : ∀ q ∈ E, nvSt ≤ q.1 ∧ q.1 < q.2 ∧ q.2 < nvEn) (hne : neEn = neSt + E.length + 1) :
    (run node nvSt nvEn neSt neEn E).edges.size = E.length + 1 ∧
    ∀ j, j < E.length + 1 → (run node nvSt nvEn neSt neEn E).edges[j]! =
      if j = 0 then ⟨node, none, (nvSt, nvEn - 1)⟩ else ⟨node, none, E.getD (j - 1) default⟩ := by
  simp only [run]
  set l₀ := inc nvSt (inc nvSt (Layout.empty (nvEn - nvSt) (neEn - neSt)) (2 * nvSt + 2))
    (2 * nvEn - 1) with hl₀
  have hsz0 : l₀.adjBounds.size = 2 * (nvEn - nvSt) + 1 := by simp [hl₀, inc, Layout.empty]
  have he0 : l₀.edges = (Layout.empty (nvEn - nvSt) (neEn - neSt)).edges := rfl
  set l₁ := E.foldl (countStep nvSt) l₀ with hl₁
  obtain ⟨e1, -, sz1, -⟩ := foldl_countStep nvSt nvEn E hE l₀ hsz0
  obtain ⟨e2, -, -, -, -⟩ := foldl_prefixStep nvSt (2 * nvEn + 1 - (2 * nvSt + 1)) (2 * nvSt + 1)
    l₁ (2 * neSt) (Nat.le_refl _) (by rw [sz1, hsz0]; omega)
  set l₂ := ((List.range' (2 * nvSt + 1) (2 * nvEn + 1 - (2 * nvSt + 1))).foldl (prefixStep nvSt)
    (l₁, 2 * neSt)).1 with hl₂
  have he2 : (inc nvSt l₂ (2 * nvSt + 2)).edges = (Layout.empty (nvEn - nvSt) (neEn - neSt)).edges := by
    show l₂.edges = _
    rw [e2, e1, he0]
  have hesz : (inc nvSt l₂ (2 * nvSt + 2)).edges.size = E.length + 1 := by
    rw [he2, empty_edges_size]; omega
  obtain ⟨sz3, g3⟩ := foldl_fillStep_edges nvSt neSt node E.reverse (inc nvSt l₂ (2 * nvSt + 2)) neEn
    (by simp; omega) (by rw [hesz]; omega)
  set l₃ := (E.reverse.foldl (fillStep nvSt neSt node) (inc nvSt l₂ (2 * nvSt + 2), neEn)).1 with hl₃
  have hsz : (inc nvSt (l₃.setNe neSt node neSt (nvSt, nvEn - 1) (2 * neSt, 2 * neEn - 1))
      (2 * nvEn - 1)).edges.size = E.length + 1 := by
    show (l₃.setNe neSt node neSt (nvSt, nvEn - 1) (2 * neSt, 2 * neEn - 1)).edges.size = _
    rw [setNe_edges_size, sz3, hesz]
  refine ⟨hsz, fun j hj => ?_⟩
  show (l₃.setNe neSt node neSt (nvSt, nvEn - 1) (2 * neSt, 2 * neEn - 1)).edges[j]! = _
  rw [setNe_edges_get _ _ _ _ _ _ _ (by rw [sz3, hesz]; exact hj), Nat.sub_self]
  by_cases hj0 : j = 0
  · rw [ite_of_pos hj0, ite_of_pos hj0]
  · rw [ite_of_neg hj0, ite_of_neg hj0, g3 j (by rw [hesz]; exact hj)]
    simp only [List.length_reverse]
    rw [ite_of_pos (by omega), List.getD_eq_getElem?_getD, List.getD_eq_getElem?_getD,
      List.getElem?_reverse (by omega)]
    congr 3; omega

/-- The R skeleton as a list: the cap followed by `E`. -/
theorem run_skeleton (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat)) (hv : nvSt + 2 ≤ nvEn)
    (hE : ∀ q ∈ E, nvSt ≤ q.1 ∧ q.1 < q.2 ∧ q.2 < nvEn) (hne : neEn = neSt + E.length + 1) :
    (run node nvSt nvEn neSt neEn E).edges.toList.map (·.nvs) = (nvSt, nvEn - 1) :: E := by
  obtain ⟨hsz, hget⟩ := run_edges node nvSt nvEn neSt neEn E hv hE hne
  apply List.ext_getElem
  · simp [hsz]
  · intro k h1 h2
    simp only [List.getElem_map, Array.getElem_toList]
    rw [← getElem!_pos _ _ (by simp at h1; omega), hget k (by simp at h1; omega)]
    cases k with
    | zero => rfl
    | succ k =>
      simp only [Nat.succ_ne_zero, ↓reduceIte, Nat.add_sub_cancel, List.getElem_cons_succ]
      rw [List.getD_eq_getElem?_getD, List.getElem?_eq_getElem (by simpa using h2)]
      rfl

/-- Bounds of an R layout, in terms of `LayoutR.start`. -/
theorem run_rowBound (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat)) (hv : nvSt + 2 ≤ nvEn)
    (hE : ∀ q ∈ E, nvSt ≤ q.1 ∧ q.1 < q.2 ∧ q.2 < nvEn) (hne : neEn = neSt + E.length + 1)
    (r : Nat) (hr1 : 2 * nvSt ≤ r) (hr2 : r ≤ 2 * nvEn) :
    rowBound nvSt neSt (run node nvSt nvEn neSt neEn E) r =
      start nvSt neSt (allE nvSt nvEn E) (r + 1) := by
  obtain ⟨hb, -⟩ := run_spec node nvSt nvEn neSt neEn E hv hE hne
  unfold rowBound
  split
  · subst ‹r = 2 * nvSt›; rw [start_low _ _ _ _ (Nat.le_refl _)]
  · rw [hb r (by omega) (by omega)]

theorem run_rowBound_mono (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat)) (hv : nvSt + 2 ≤ nvEn)
    (hE : ∀ q ∈ E, nvSt ≤ q.1 ∧ q.1 < q.2 ∧ q.2 < nvEn) (hne : neEn = neSt + E.length + 1)
    (r : Nat) (hr1 : 2 * nvSt ≤ r) (hr2 : r < 2 * nvEn) :
    rowBound nvSt neSt (run node nvSt nvEn neSt neEn E) r ≤
      rowBound nvSt neSt (run node nvSt nvEn neSt neEn E) (r + 1) := by
  rw [run_rowBound _ _ _ _ _ _ hv hE hne _ hr1 (by omega),
    run_rowBound _ _ _ _ _ _ hv hE hne _ (by omega) (by omega)]
  exact start_le _ _ _ (by omega)

theorem run_rowBound_last (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat)) (hv : nvSt + 2 ≤ nvEn)
    (hE : ∀ q ∈ E, nvSt ≤ q.1 ∧ q.1 < q.2 ∧ q.2 < nvEn) (hne : neEn = neSt + E.length + 1) :
    rowBound nvSt neSt (run node nvSt nvEn neSt neEn E) (2 * nvEn) = 2 * neEn := by
  rw [run_rowBound _ _ _ _ _ _ hv hE hne _ (by omega) (Nat.le_refl _),
    start_top _ _ _ _ (allE_bounds nvSt nvEn E hv hE) (by omega), allE_length]
  omega

/-- Row destinations of an R layout (`LayoutR.row_destNv` restated on `run`). -/
theorem run_row_destNv (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat)) (hv : nvSt + 2 ≤ nvEn)
    (hE : ∀ q ∈ E, nvSt ≤ q.1 ∧ q.1 < q.2 ∧ q.2 < nvEn) (hne : neEn = neSt + E.length + 1)
    (r : Nat) (hr1 : 2 * nvSt ≤ r) (hr2 : r ≤ 2 * nvEn - 1) :
    (row nvSt neSt (run node nvSt nvEn neSt neEn E) r).map (·.destNv) = rowDest nvSt nvEn E (r + 1) := by
  rw [← layoutNode_R_eq _ _ _ _ _ _ hv]
  exact row_destNv node nvSt nvEn neSt neEn E hv hE hne r hr1 hr2

end R

end LayoutShape

end Spqr
