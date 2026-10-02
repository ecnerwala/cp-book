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

end LayoutShape

end Spqr
