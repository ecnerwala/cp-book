import Spqr.PlanarEmbedSteps

/-!
# The edges below an item, by children

Under `Preorder.subtree_eq` and `ChildShape.tile` the subtree range of an item is the item itself
followed by the subtree ranges of its children in order, so `edgesBelow i` is the item's own edge
(for a `Q` item) followed by the `edgesBelow` of its children.
-/

namespace Spqr

namespace PlanarSpqrTree

variable (t : PlanarSpqrTree)

theorem getElem!_nat (a : Array Nat) (j : Nat) : a[j]! = a[j]?.getD 0 := by
  rw [getElem!_def]; cases a[j]? <;> rfl

theorem children_eq (i : Nat) : t.toSpqrTree.children i = t.children i := by
  unfold SpqrTree.children PlanarSpqrTree.children SpqrTree.chRange
  simp only [getElem!_nat]

/-- The subtree sizes of the first `k` children add up to the range from `i + 1` to the `k`-th
child. -/
theorem range'_children (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape) (i : Nat)
    (hi : i < t.size) :
    List.range' (i + 1) (t.subtreeEnd[i]! - (i + 1)) =
      (t.children i).flatMap fun c => List.range' c (t.subtreeEnd[c]! - c) := by
  have hsub := hwf.preorder.subtree_eq i hi
  have htile := hsh.tile i hi
  rw [t.children_eq] at hsub htile
  have key : ∀ k, k ≤ (t.children i).length →
      List.range' (i + 1) ((((t.children i).take k).map fun c => t.subtreeEnd[c]! - c).sum) =
        ((t.children i).take k).flatMap fun c => List.range' c (t.subtreeEnd[c]! - c) := by
    intro k
    induction k with
    | zero => intro _; simp
    | succ k ih =>
      intro hk
      have hk' : k < (t.children i).length := hk
      rw [List.take_add_one, List.getElem?_eq_getElem hk', Option.toList_some, List.map_append,
        List.sum_append, List.map_singleton, List.sum_singleton, List.flatMap_append,
        ← List.range'_append_1, ih (Nat.le_of_succ_le hk), List.flatMap_singleton]
      congr 2
      have := htile k hk'
      rw [getElem!_pos _ k hk'] at this
      omega
  have := key _ le_rfl
  rw [List.take_length] at this
  rw [hsub, Nat.add_sub_cancel_left, this]

/-- `edgesBelow i` is the item's own edge (a `Q` item) followed by its children's. -/
theorem edgesBelow_eq (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape) (i : Nat)
    (hi : i < t.size) :
    t.edgesBelow i =
      (if t.types[i]! == .Q then t.origId[i]! else none).toList ++
        (t.children i).flatMap t.edgesBelow := by
  have hsub := hwf.preorder.subtree_eq i hi
  have hle : i + 1 ≤ t.subtreeEnd[i]! := by rw [hsub]; omega
  unfold edgesBelow
  rw [show t.subtreeEnd[i]! - i = (t.subtreeEnd[i]! - (i + 1)) + 1 by omega, List.range'_succ,
    List.filterMap_cons, t.range'_children hwf hsh i hi, List.filterMap_flatMap]
  split <;> rename_i heq <;> rw [heq] <;> rfl

end PlanarSpqrTree

end Spqr
