import Spqr.PlanarEmbedQExec

/-!
# Unfolding the `S`/`P`/`R` branch of `embedItem`

The node loop over the quarter-edges `[4 * neSt : 4 * neEn]` as a `List.foldl` of the pure step
`nodeStep`: every node quarter-edge `ta` with rotation successor `tb > ta` links the exposed slots of
the children owning the twins of `ta` and `tb` (through the vertex item's exposed pair when the
corner is `(ta % 4 = 2, tb % 4 = 1)`), and the cap edge's corners expose the opposite slots.
-/

namespace Spqr.PlanarSpqrTree

open EmbedM

variable (t : PlanarSpqrTree)

theorem forIn_state_foldl {σ α : Type} (l : List α) (f : α → PUnit → StateM σ (ForInStep PUnit))
    (g : σ → α → σ) (hf : ∀ a s, (f a ()).run s = (ForInStep.yield (), g s a)) (s : σ) :
    (forIn l () f).run s = ((), l.foldl g s) := by
  induction l generalizing s with
  | nil => rfl
  | cons a l ih =>
    rw [List.forIn_cons, StateT.run_bind, hf]
    exact ih _

/-- The exposed slot of the child owning the twin of the node quarter-edge `q`. -/
def treeQe (s : EmbedState) (q : Nat) : Option Nat :=
  s.outerE[t.nodeEdges[(t.nodeEdges[QE.edge q]!.twin).getD 0]!.node]![2 * QE.side q + QE.dir q]!

/-- One iteration of the node loop at quarter-edge `ta`. -/
def nodeStep (i neSt : Nat) (s : EmbedState) (ta : Nat) : EmbedState :=
  match t.neRotAdj[ta]! with
  | none => s
  | some tb =>
    if tb < ta then s
    else
      let qb := t.treeQe s tb
      if ta < 4 * (neSt + 1) then ((setOuter i (2 * QE.side ta + (1 - QE.dir ta)) qb).run s).2
      else
        let qa := t.treeQe s ta
        if ta % 4 == 2 && tb % 4 == 1 then
          let v := t.nodeVerts[(t.nodeEdges[QE.edge ta]!).nvs.2]!.vert
          if s.outerE[v]![0]!.isSome then
            ((link qb s.outerE[v]![0]!).run ((link qa s.outerE[v]![1]!).run s).2).2
          else ((link qa qb).run s).2
        else ((link qa qb).run s).2

theorem embedItem_node (i : Nat) (ht : t.types[i]! = .S ∨ t.types[i]! = .P ∨ t.types[i]! = .R)
    (s : EmbedState) :
    ((t.embedItem i).run s).2 =
      (List.range' (4 * t.neBounds[i]!) (4 * t.neBounds[i + 1]! - 4 * t.neBounds[i]!)).foldl
        (t.nodeStep i t.neBounds[i]!) s := by
  have idbind {α β : Type} (a : Id α) (f : α → Id β) : (a >>= f) = f a := rfl
  unfold embedItem
  rcases ht with ht | ht | ht <;> rw [ht] <;>
  simp only [Std.Legacy.Range.forIn_eq_forIn_range', Std.Legacy.Range.size, Nat.add_sub_cancel,
    Nat.div_one]
  all_goals
    refine (congrArg Prod.snd (forIn_state_foldl _ _ (t.nodeStep i t.neBounds[i]!) ?_ s)).trans rfl
    intro ta s
    show (_ : EmbedM (ForInStep PUnit)).run s = _
    unfold nodeStep
    rcases h : t.neRotAdj[ta]! with _ | tb
    · rfl
    · simp only []
      by_cases h1 : tb < ta
      · simp only [h1, ↓reduceIte]; rfl
      · simp only [h1, ↓reduceIte, StateT.run_bind, outer_run, idbind, treeQe]
        split_ifs <;> simp only [*, ↓reduceIte, StateT.run_bind, outer_run, link_outerE] <;> rfl
theorem nodeStep_outerE_ne (i neSt : Nat) (s : EmbedState) (ta j : Nat) (hj : j ≠ i) :
    (t.nodeStep i neSt s ta).outerE[j]? = s.outerE[j]? := by
  unfold nodeStep
  split
  · rfl
  · dsimp only
    split_ifs <;> simp only [link_outerE, setOuter_outerE_ne _ _ _ _ _ hj]

theorem nodeStep_outerE_size (i neSt : Nat) (s : EmbedState) (ta : Nat) :
    (t.nodeStep i neSt s ta).outerE.size = s.outerE.size := by
  unfold nodeStep
  split
  · rfl
  · dsimp only
    split_ifs <;> simp only [link_outerE, setOuter_outerE_size]

theorem nodeStep_rotAdj_size (i neSt : Nat) (s : EmbedState) (ta : Nat) :
    (t.nodeStep i neSt s ta).rotAdj.size = s.rotAdj.size := by
  unfold nodeStep
  split
  · rfl
  · dsimp only
    split_ifs <;> simp only [link_rotAdj_size, setOuter_rotAdj]

theorem nodeStep_exposedAt_ne (i neSt : Nat) (s : EmbedState) (ta j q : Nat) (hj : j ≠ i) :
    (t.nodeStep i neSt s ta).exposedAt j q ↔ s.exposedAt j q := by
  simp [EmbedState.exposedAt, nodeStep_outerE_ne _ _ _ _ _ _ hj]

theorem nodeFold_outerE_ne (i neSt : Nat) (l : List Nat) (s : EmbedState) (j : Nat) (hj : j ≠ i) :
    (l.foldl (t.nodeStep i neSt) s).outerE[j]? = s.outerE[j]? := by
  induction l generalizing s with
  | nil => rfl
  | cons a l ih => rw [List.foldl_cons, ih, nodeStep_outerE_ne _ _ _ _ _ _ hj]

theorem nodeFold_outerE_size (i neSt : Nat) (l : List Nat) (s : EmbedState) :
    (l.foldl (t.nodeStep i neSt) s).outerE.size = s.outerE.size := by
  induction l generalizing s with
  | nil => rfl
  | cons a l ih => rw [List.foldl_cons, ih, nodeStep_outerE_size]

theorem nodeFold_rotAdj_size (i neSt : Nat) (l : List Nat) (s : EmbedState) :
    (l.foldl (t.nodeStep i neSt) s).rotAdj.size = s.rotAdj.size := by
  induction l generalizing s with
  | nil => rfl
  | cons a l ih => rw [List.foldl_cons, ih, nodeStep_rotAdj_size]

theorem nodeFold_exposedAt_ne (i neSt : Nat) (l : List Nat) (s : EmbedState) (j q : Nat) (hj : j ≠ i) :
    (l.foldl (t.nodeStep i neSt) s).exposedAt j q ↔ s.exposedAt j q := by
  simp [EmbedState.exposedAt, nodeFold_outerE_ne _ _ _ _ _ _ hj]

end Spqr.PlanarSpqrTree
