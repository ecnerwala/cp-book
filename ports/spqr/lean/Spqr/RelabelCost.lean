import Lean
import Spqr.Relabel
import Spqr.WalkCost

/-!
# Step bound for phase 3

`relabel` ticks `1 + [cap edge] + [two own node-verts]` when it numbers a node and `2` per child
when it lays out the child list. Against the potential `4 * types.size + 2 * chDat.size` each
primitive action is non-increasing, so `ticks ≤ 4 * #nodes numbered + 2 * #child slots`. Turning
that into `O(V + E)` needs the (unproved, see `relabelRun_sizes_le`) fact that the item tree is a
tree, so that `relabel` visits each item at most once.
-/

namespace Spqr
open StateRun

def RelabelState.pot (s : RelabelState) : Nat := 4 * s.types.size + 2 * s.chDat.size

/-- `m` ticks at most as much as it raises the potential. -/
def ROk (m : RelabelM α) : Prop := ∀ s, (m.run s).2.ticks + s.pot ≤ s.ticks + (m.run s).2.pot

attribute [run_simp] RelabelState.pot Array.size_push Array.size_append List.size_toArray Array.size_set!
  Array.size_setIfInBounds Array.size_modify

namespace ROk

variable {α β : Type}

theorem pure (a : α) : ROk (pure a : RelabelM α) := fun s => by run_simp; omega
theorem bind {m : RelabelM α} {f : α → RelabelM β} (hm : ROk m) (hf : ∀ a, ROk (f a)) :
    ROk (m >>= f) := fun s => by
  rw [run_bind, runBind_eq]
  have h1 := hm s
  have h2 := hf (m.run s).1 (m.run s).2
  omega
theorem map (g : α → β) {m : RelabelM α} (hm : ROk m) : ROk (g <$> m) := fun s => by
  rw [run_map, runMap_eq]; exact hm s
theorem ite (c : Prop) [Decidable c] {m1 m2 : RelabelM α} (h1 : ROk m1) (h2 : ROk m2) :
    ROk (if c then m1 else m2) := by
  split <;> assumption
theorem dite (c : Prop) [Decidable c] {m1 : c → RelabelM α} {m2 : ¬c → RelabelM α}
    (h1 : ∀ h, ROk (m1 h)) (h2 : ∀ h, ROk (m2 h)) : ROk (dite c m1 m2) := by
  split
  · exact h1 _
  · exact h2 _
theorem forIn_list (l : List α) (init : β) (f : α → β → RelabelM (ForInStep β))
    (hf : ∀ a b, ROk (f a b)) : ROk (forIn l init f) := by
  induction l generalizing init with
  | nil => simp only [List.forIn_nil]; exact pure _
  | cons a l ih =>
    simp only [List.forIn_cons]
    refine bind (hf a init) fun r => ?_
    cases r with
    | done b => exact pure _
    | yield b => exact ih b
theorem get : ROk (get : RelabelM RelabelState) := fun s => by run_simp; omega
theorem modify (f : RelabelState → RelabelState) (h : ∀ s, (f s).ticks + s.pot ≤ s.ticks + (f s).pot) :
    ROk (modify f : RelabelM Unit) := fun s => by
  run_simp; exact h s
theorem item (i : ItemId) : ROk (RelabelM.item i) := fun s => by
  simp only [run_simp, RelabelM.item]; omega
theorem orderedChildren (it : Item) (ch : List ItemId) (n : Nat) :
    ROk (RelabelM.orderedChildren it ch n) := fun s => by
  simp only [run_simp, RelabelM.orderedChildren]
  split <;> (try dsimp only) <;> omega

end ROk

/-! ## `ok_auto` -/

open Lean Elab Tactic Meta in
def okLeaf : Name → Option Name
  | ``Pure.pure => some ``ROk.pure
  | ``MonadState.get => some ``ROk.get
  | ``MonadStateOf.get => some ``ROk.get
  | ``getThe => some ``ROk.get
  | ``RelabelM.item => some ``ROk.item
  | ``RelabelM.orderedChildren => some ``ROk.orderedChildren
  | _ => none

open Lean Elab Tactic Meta in
/-- Solve `ROk m` by recursion on the syntax of `m`. Unknown actions are closed by a hypothesis;
`modify` leaves that `omega` cannot close are left to the caller. -/
def okAuto (leftovers : IO.Ref (Array MVarId)) : Nat → MVarId → TacticM Unit
  | 0, _ => throwError "ok_auto: out of fuel"
  | fuel + 1, g => g.withContext do
    if ← g.isAssigned then return
    let ty := (← instantiateMVars (← g.getType)).consumeMData
    if ty.isForall then
      let (_, g') ← g.intro1
      return ← okAuto leftovers fuel g'
    unless ty.isAppOfArity ``ROk 2 do
      throwError "ok_auto: not an ROk goal{indentExpr ty} ({ty.getAppFn} / {ty.getAppNumArgs})"
    let m := (ty.getArg! 1).headBeta
    let recurse (gs : List MVarId) : TacticM Unit := do
      for g' in gs do
        unless ← g'.isAssigned do
          let t := (← instantiateMVars (← g'.getType)).consumeMData
          if t.isForall || t.isAppOf ``ROk then okAuto leftovers fuel g'
    let applyLemma (n : Name) : TacticM (List MVarId) := do
      g.apply (← mkConstWithFreshMVarLevels n)
    let changeTo (m' : Expr) : TacticM Unit := do
      okAuto leftovers fuel (← g.change (mkApp2 ty.getAppFn (ty.getArg! 0) m'))
    match m with
    | .letE _ _ v b _ => changeTo (b.instantiate1 v)
    | _ =>
    let .const n _ := m.getAppFn | throwError "ok_auto: unexpected action{indentExpr m}"
    if n == ``Bind.bind then recurse (← applyLemma ``ROk.bind)
    else if n == ``Functor.map then recurse (← applyLemma ``ROk.map)
    else if n == ``ite then recurse (← applyLemma ``ROk.ite)
    else if n == ``dite then recurse (← applyLemma ``ROk.dite)
    else if n == ``letFun then changeTo ((m.getArg! 3).beta #[m.getArg! 2])
    else if n == ``ForIn.forIn then recurse (← applyLemma ``ROk.forIn_list)
    else if n == ``modify then
      for g' in ← applyLemma ``ROk.modify do
        unless ← g'.isAssigned do
          if (← instantiateMVars (← g'.getType)).isForall then
            setGoals [g']
            try evalTactic (← `(tactic| intro s; run_simp; (repeat' split) <;> (try dsimp only) <;> omega))
            catch _ => leftovers.modify (·.push g')
    else if let some lem := okLeaf n then recurse (← applyLemma lem)
    else if ← isMatcherApp m then
      if let .reduced m' ← reduceMatcher? m then return ← changeTo m'
      setGoals [g]
      evalTactic (← `(tactic| split))
      recurse (← getUnsolvedGoals)
    else
      for d in ← getLCtx do
        if d.isImplementationDetail then continue
        try
          let gs ← g.apply d.toExpr
          if gs.isEmpty then return
          throwError "partial"
        catch _ => pure ()
      throwError "ok_auto: no rule for{indentExpr m}"

open Lean Elab Tactic in
elab "ok_auto" : tactic => do
  let g ← getMainGoal
  let rest ← getGoals
  let leftovers ← IO.mkRef #[]
  okAuto leftovers 100000 g
  setGoals ((← leftovers.get).toList ++ (← rest.filterM fun g => return !(← g.isAssigned)))

theorem ok_relabel : ∀ fuel cur parent parNv capTwin, ROk (relabel fuel cur parent parNv capTwin)
  | 0, _, _, _, _ => ROk.pure ()
  | fuel + 1, cur, parent, parNv, capTwin => by
    have ih := ok_relabel fuel
    unfold relabel
    ok_auto

/-- Phase 3 ticks at most `4` per node numbered and `2` per child slot written. -/
theorem relabel_ticks_le_sizes (g : Graph) (items : Array Item) :
    (relabelRun g items).ticks ≤
      4 * (relabelRun g items).types.size + 2 * (relabelRun g items).chDat.size := by
  have h := ok_relabel items.size rootItem none none none (RelabelState.init g items)
  unfold relabelRun
  simp only [RelabelState.pot, RelabelState.init, List.size_toArray, List.length_nil] at h ⊢
  omega

/-- `relabel` numbers at most `items.size` nodes and writes at most `items.size` child slots.
MISSING LEMMA (`sorry`): the items produced by `Graph.walk` form a rooted tree under `Item.ch`
(each item is a child of at most one item, none is its own ancestor), so the preorder traversal
visits each item at most once. That walk invariant is not proved here. -/
theorem relabelRun_sizes_le (g : Graph) (tern : Bool) (forest : List DfsTree) :
    (relabelRun g (g.walk tern forest).items).types.size ≤ (g.walk tern forest).items.size ∧
    (relabelRun g (g.walk tern forest).items).chDat.size ≤ (g.walk tern forest).items.size := by
  sorry

/-- Steps per vertex / edge in phase 3: `6` per item, `items ≤ 1 + (V + E) + walk ticks`. -/
abbrev relabelC : Nat := 6 * (walkC + 2)

/-- Phase 3 is `O(V + E)`: `288 * (V + E) + 6` steps. Depends on `relabelRun_sizes_le` and
`dfsForest_size_le`. -/
theorem relabel_ticks_le (g : Graph) (tern : Bool) (vertOrder edgeOrder : List Nat) :
    (relabelRun g (g.walk tern (g.dfsForest vertOrder edgeOrder)).items).ticks ≤
      relabelC * (g.nv + g.ne) + 6 := by
  have h1 := relabel_ticks_le_sizes g (g.walk tern (g.dfsForest vertOrder edgeOrder)).items
  have h2 := relabelRun_sizes_le g tern (g.dfsForest vertOrder edgeOrder)
  have h3 := walk_items_le g tern (g.dfsForest vertOrder edgeOrder)
  have h4 := walk_ticks_le g tern vertOrder edgeOrder
  simp only [relabelC, walkC] at *
  omega

end Spqr
