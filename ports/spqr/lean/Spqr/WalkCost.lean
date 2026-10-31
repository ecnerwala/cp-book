import Lean
import Spqr.WalkFast
import Spqr.Proofs.Dfs
import Spqr.StateRun

/-!
# Step bound for the walk

`WalkState.ticks` counts loop iterations, tstack pushes / pops and item allocations. The bound is
amortized through the potential `pot = ticks + 4 * |tstack|`: every ear pushed prepays the loop
iteration that eventually pops it, so each `loop` of `finishEdge` has non-positive amortized cost
and the walk of a DFS forest costs `O(#nodes + #out-edges) = O(V + E)`.

The proofs go through `Cost m w` ("running `m` raises `pot` by at most `w`, and allocates at most
as many items as it ticks"), a Hoare-style `Spec` for the loop bodies, and a small tactic
`cost_auto` that walks the syntax of a `do` block and adds up the costs of its primitives.
-/

namespace Spqr.Fast
open StateRun

def WalkState.pot (s : WalkState) : Nat := s.ticks + 4 * s.tstack.size

/-- Going from `s` to `s'` allocated at most as many items as it ticked. -/
def WalkState.Alloc (s s' : WalkState) : Prop := s'.items.size + s.ticks ≤ s.items.size + s'.ticks

theorem WalkState.Alloc.trans {s s' s'' : WalkState} (h1 : s.Alloc s') (h2 : s'.Alloc s'') : s.Alloc s'' := by
  unfold WalkState.Alloc at *; omega

namespace WalkM

/-- Running `m` raises the potential by at most `w` (and allocates at most as it ticks). -/
def Cost (m : WalkM α) (w : Nat) : Prop :=
  ∀ s, (m.run s).2.pot ≤ s.pot + w ∧ s.Alloc (m.run s).2

/-- From a state satisfying `P`, `m` ends in a state satisfying `Q` with `pot' + b ≤ pot + a`. -/
def Spec (P : WalkState → Prop) (m : WalkM α) (Q : WalkState → Prop) (a b : Nat) : Prop :=
  ∀ s, P s → Q (m.run s).2 ∧ (m.run s).2.pot + b ≤ s.pot + a ∧ s.Alloc (m.run s).2

attribute [run_simp] WalkState.pot WalkState.Alloc Array.size_push Array.size_pop Array.size_modify Array.size_set!
  Array.size_setIfInBounds

section Rules
variable {m : WalkM α} {f : α → WalkM β}

theorem Cost.mono (h : Cost m w) (hw : w ≤ w') : Cost m w' :=
  fun s => ⟨Nat.le_trans (h s).1 (by omega), (h s).2⟩

theorem cost_pure (a : α) : Cost (pure a : WalkM α) 0 := fun s => by simp only [run_simp]; omega
theorem cost_bind (hm : Cost m w1) (hf : ∀ a, Cost (f a) w2) : Cost (m >>= f) (w1 + w2) := fun s => by
  rw [run_bind, runBind_eq]
  obtain ⟨h1, a1⟩ := hm s
  obtain ⟨h2, a2⟩ := hf (m.run s).1 (m.run s).2
  exact ⟨by omega, a1.trans a2⟩
theorem cost_map (g : α → β) (hm : Cost m w) : Cost (g <$> m) w := fun s => by rw [run_map, runMap_eq]; exact hm s
theorem cost_ite (c : Prop) [Decidable c] {m1 m2 : WalkM α} (h1 : Cost m1 w1) (h2 : Cost m2 w2) :
    Cost (if c then m1 else m2) (max w1 w2) := by
  split
  · exact h1.mono (Nat.le_max_left _ _)
  · exact h2.mono (Nat.le_max_right _ _)
theorem cost_dite (c : Prop) [Decidable c] {m1 : c → WalkM α} {m2 : ¬c → WalkM α}
    (h1 : ∀ h, Cost (m1 h) w1) (h2 : ∀ h, Cost (m2 h) w2) : Cost (dite c m1 m2) (max w1 w2) := by
  split
  · exact (h1 _).mono (Nat.le_max_left _ _)
  · exact (h2 _).mono (Nat.le_max_right _ _)
theorem cost_get : Cost (get : WalkM WalkState) 0 := fun s => by simp only [run_simp]; omega
theorem cost_modify (g : WalkState → WalkState)
    (h : ∀ s, (g s).pot ≤ s.pot ∧ (g s).items.size ≤ s.items.size ∧ s.ticks ≤ (g s).ticks) :
    Cost (modify g : WalkM Unit) 0 := fun s => by
  have := h s; simp only [run_simp] at *; omega

theorem Spec.weaken (h : Spec P m Q a b) (ha : a ≤ a') (hb : b' ≤ b) : Spec P m Q a' b' :=
  fun s hs => ⟨(h s hs).1, by have := (h s hs).2.1; omega, (h s hs).2.2⟩
theorem Spec.conseq (h : Spec P m Q a b) (hP : ∀ s, P' s → P s) (hQ : ∀ s, Q s → Q' s) :
    Spec P' m Q' a b :=
  fun s hs => ⟨hQ _ (h s (hP s hs)).1, (h s (hP s hs)).2⟩
theorem spec_bind (hm : Spec P m Q a1 b1) (hf : ∀ x, Spec Q (f x) R a2 b2) :
    Spec P (m >>= f) R (a1 + a2) (b1 + b2) := fun s hs => by
  rw [run_bind, runBind_eq]
  obtain ⟨hQ, h1, al1⟩ := hm s hs
  obtain ⟨hR, h2, al2⟩ := hf _ _ hQ
  exact ⟨hR, by omega, al1.trans al2⟩
theorem spec_ite (c : Prop) [Decidable c] {m1 m2 : WalkM α} (h1 : Spec P m1 Q a b)
    (h2 : Spec P m2 Q a b) : Spec P (if c then m1 else m2) Q a b := by
  split <;> assumption
theorem spec_pure (x : α) (hPQ : ∀ s, P s → Q s) : Spec P (pure x : WalkM α) Q 0 0 :=
  fun s hs => ⟨hPQ s hs, by simp only [run_simp]; omega, by simp only [run_simp]; omega⟩
/-- A step that never pops keeps any lower bound on the stack size. -/
theorem spec_of_cost (h : Cost m w) (hsz : ∀ s, s.tstack.size ≤ (m.run s).2.tstack.size) :
    Spec (fun s => k ≤ s.tstack.size) m (fun s => k ≤ s.tstack.size) w 0 :=
  fun s hs => ⟨Nat.le_trans hs (hsz s), by have := (h s).1; omega, (h s).2⟩

theorem Spec.pays (h : Spec P m Q a b) (hab : a + 1 ≤ b) :
    ∀ s, P s → (m.run s).2.pot + 1 ≤ s.pot ∧ s.Alloc (m.run s).2 :=
  fun s hs => ⟨by have := (h s hs).2.1; omega, (h s hs).2.2⟩

end Rules

section Leaves

theorem cost_tick : Cost tick 1 := fun s => by simp only [run_simp, tick]; omega
theorem cost_modifyItem (i : ItemId) (g : Item → Item) : Cost (modifyItem i g) 0 := fun s => by
  simp only [run_simp, modifyItem]; omega
theorem cost_getItem (i : ItemId) : Cost (getItem i) 0 := fun s => by simp only [run_simp, getItem]; omega
theorem cost_allocItem (t : NodeType) : Cost (allocItem t) 1 := fun s => by
  simp only [run_simp, allocItem]; omega
theorem cost_stackDir (d : Nat) : Cost (stackDir d) 0 := fun s => by simp only [run_simp, stackDir]; omega
theorem cost_setStackDir (d : Nat) (b : Bool) : Cost (setStackDir d b) 0 := fun s => by
  simp only [run_simp, setStackDir]; omega
theorem cost_makeVs (v d : Nat) : Cost (makeVs v d) 0 := fun s => by simp only [run_simp, makeVs]; omega
theorem cost_cur : Cost cur 0 := fun s => by simp only [run_simp, cur]; omega
theorem cost_nxt : Cost nxt 0 := fun s => by simp only [run_simp, nxt]; omega
theorem cost_tstackSize : Cost tstackSize 0 := fun s => by simp only [run_simp, tstackSize]; omega
theorem cost_modifyCur (g : TEntry → TEntry) : Cost (modifyCur g) 0 := fun s => by
  simp only [run_simp, modifyCur]; omega
theorem cost_modifyNxt (g : TEntry → TEntry) : Cost (modifyNxt g) 0 := fun s => by
  simp only [run_simp, modifyNxt]; split <;> (try simp only [run_simp]) <;> omega
theorem cost_popTstack : Cost popTstack 1 := fun s => by simp only [run_simp, popTstack]; omega
theorem cost_pushTstack (v d i : Nat) : Cost (pushTstack v d i) 5 := fun s => by
  simp only [run_simp, pushTstack]; omega
theorem cost_pushVertTstack (v d : Nat) : Cost (pushVertTstack v d) 5 := cost_pushTstack _ _ _
theorem cost_pushEdgeTstack (v d e : Nat) : Cost (pushEdgeTstack v d e) 5 :=
  (cost_bind cost_get fun _ => cost_pushTstack _ _ _).mono (by omega)
theorem cost_mergeTstackTops : Cost mergeTstackTops 1 :=
  (cost_bind cost_popTstack fun _ => cost_modifyCur _).mono (by omega)

theorem sz_nxt (s : WalkState) : s.tstack.size ≤ (nxt.run s).2.tstack.size := by simp only [run_simp, nxt]; omega
theorem sz_cur (s : WalkState) : s.tstack.size ≤ (cur.run s).2.tstack.size := by simp only [run_simp, cur]; omega
theorem sz_setStackDir (d : Nat) (b : Bool) (s : WalkState) :
    s.tstack.size ≤ ((setStackDir d b).run s).2.tstack.size := by simp only [run_simp, setStackDir]; omega
theorem sz_modifyCur (g : TEntry → TEntry) (s : WalkState) :
    s.tstack.size ≤ ((modifyCur g).run s).2.tstack.size := by simp only [run_simp, modifyCur]; omega

/-- A real pop: from `k + 1` entries, `popTstack` pays 4 for its one tick. -/
theorem spec_popTstack :
    Spec (fun s => k + 1 ≤ s.tstack.size) popTstack (fun s => k ≤ s.tstack.size) 1 4 := fun s hs => by
  simp only [run_simp, popTstack]; omega
theorem spec_mergeTstackTops :
    Spec (fun s => k + 1 ≤ s.tstack.size) mergeTstackTops (fun s => k ≤ s.tstack.size) 1 4 :=
  (spec_bind spec_popTstack fun _ => spec_of_cost (cost_modifyCur _) (sz_modifyCur _)).weaken
    (by omega) (by omega)

end Leaves

section Loops

theorem cost_loop {cond : WalkM Bool} {body : WalkM Unit} {P : WalkState → Prop}
    (hc : ∀ s, (cond.run s).2 = s) (hcP : ∀ s, (cond.run s).1 = true → P s)
    (hb : ∀ s, P s → (body.run s).2.pot + 1 ≤ s.pot ∧ s.Alloc (body.run s).2) :
    ∀ n, Cost (loop n cond body) 0
  | 0 => cost_pure ()
  | n + 1 => fun s => by
    have ih := cost_loop hc hcP hb n
    simp only [loop, run_bind, runBind_eq, hc]
    split
    · rename_i h
      simp only [run_bind, runBind_eq]
      obtain ⟨h1, a1⟩ := hb s (hcP s h)
      obtain ⟨h2, a2⟩ := ih (tick.run (body.run s).2).2
      have h3 : (tick.run (body.run s).2).2.pot = (body.run s).2.pot + 1 ∧
          (body.run s).2.Alloc (tick.run (body.run s).2).2 := by simp only [run_simp, tick]; omega
      exact ⟨by omega, a1.trans (h3.2.trans a2)⟩
    · simp only [run_simp]; omega

theorem closeCond_state (d : Nat) (s : WalkState) : ((closeCond d).run s).2 = s := by
  simp only [run_simp, closeCond, tstackSize, nxt]
theorem closeCond_true (d : Nat) (s : WalkState) (h : ((closeCond d).run s).1 = true) :
    2 ≤ s.tstack.size := by
  simp only [run_simp, closeCond, tstackSize, nxt] at h
  simp at h
  omega

theorem firstIdxCond_state (fo : Nat) (s : WalkState) : ((firstIdxCond fo).run s).2 = s := by
  simp only [run_simp, firstIdxCond, cur]
theorem firstIdxCond_true (fo : Nat) (s : WalkState) (h : ((firstIdxCond fo).run s).1 = true) :
    0 + 1 ≤ s.tstack.size := by
  simp only [run_simp, firstIdxCond, cur] at h
  cases hsz : s.tstack.size with
  | zero =>
    simp [Array.eq_empty_of_size_eq_zero hsz, Array.back!, show (default : TEntry).firstIdx = 0 from rfl] at h
  | succ n => omega

theorem sizeCond_state (n : Nat) (s : WalkState) : ((sizeCond n).run s).2 = s := by
  simp only [run_simp, sizeCond, tstackSize]
theorem sizeCond_true (n : Nat) (s : WalkState) (h : ((sizeCond n).run s).1 = true) :
    0 + 1 ≤ s.tstack.size := by
  simp only [run_simp, sizeCond, tstackSize] at h
  simp at h
  omega

theorem cost_loop_firstIdx (n fo : Nat) : Cost (loop n (firstIdxCond fo) mergeTstackTops) 0 :=
  cost_loop (firstIdxCond_state fo) (firstIdxCond_true fo) ((spec_mergeTstackTops (k := 0)).pays (by omega)) n
theorem cost_loop_size (n k : Nat) : Cost (loop n (sizeCond k) mergeTstackTops) 0 :=
  cost_loop (sizeCond_state k) (sizeCond_true k) ((spec_mergeTstackTops (k := 0)).pays (by omega)) n

end Loops

/-! ## `cost_auto` -/

open Lean Elab Tactic Meta in
/-- Leaf lemma for a primitive action, by its head constant. -/
def costLeaf : Name → Option Name
  | ``tick => some ``cost_tick
  | ``Pure.pure => some ``cost_pure
  | ``MonadState.get => some ``cost_get
  | ``MonadStateOf.get => some ``cost_get
  | ``getThe => some ``cost_get
  | ``getItem => some ``cost_getItem
  | ``allocItem => some ``cost_allocItem
  | ``stackDir => some ``cost_stackDir
  | ``setStackDir => some ``cost_setStackDir
  | ``makeVs => some ``cost_makeVs
  | ``cur => some ``cost_cur
  | ``nxt => some ``cost_nxt
  | ``tstackSize => some ``cost_tstackSize
  | ``modifyItem => some ``cost_modifyItem
  | ``modifyCur => some ``cost_modifyCur
  | ``modifyNxt => some ``cost_modifyNxt
  | ``popTstack => some ``cost_popTstack
  | ``pushTstack => some ``cost_pushTstack
  | ``pushVertTstack => some ``cost_pushVertTstack
  | ``pushEdgeTstack => some ``cost_pushEdgeTstack
  | ``mergeTstackTops => some ``cost_mergeTstackTops
  | ``maybeUnwrapNxt => some `Spqr.Fast.WalkM.cost_maybeUnwrapNxt
  | ``finishTstackTop => some `Spqr.Fast.WalkM.cost_finishTstackTop
  | ``Spqr.Fast.finishEdge => some `Spqr.Fast.cost_finishEdge
  | _ => none

open Lean Elab Tactic Meta in
/-- Solve `Cost m ?w` by recursion on the syntax of `m`, computing `?w`. -/
def costAuto : Nat → MVarId → TacticM Unit
  | 0, _ => throwError "cost_auto: out of fuel"
  | fuel + 1, g => g.withContext do
    if ← g.isAssigned then return
    let ty ← instantiateMVars (← g.getType)
    if ty.isForall then
      let (_, g') ← g.intro1
      return ← costAuto fuel g'
    unless ty.isAppOfArity ``Cost 3 do
      throwError "cost_auto: not a Cost goal{indentExpr ty}"
    let m := (ty.getArg! 1).headBeta
    let recurse (gs : List MVarId) : TacticM Unit := do
      for g' in gs do
        unless ← g'.isAssigned do
          let t ← instantiateMVars (← g'.getType)
          if t.isForall || t.isAppOf ``Cost then costAuto fuel g'
    let applyLemma (n : Name) : TacticM (List MVarId) := do
      g.apply (← mkConstWithFreshMVarLevels n)
    let changeTo (m' : Expr) : TacticM Unit := do
      let g' ← g.change (mkApp3 ty.getAppFn (ty.getArg! 0) m' (ty.getArg! 2))
      costAuto fuel g'
    match m with
    | .letE _ _ v b _ => changeTo (b.instantiate1 v)
    | _ =>
    let .const n _ := m.getAppFn | throwError "cost_auto: unexpected action{indentExpr m}"
    if n == ``Bind.bind then recurse (← applyLemma ``cost_bind)
    else if n == ``Functor.map then recurse (← applyLemma ``cost_map)
    else if n == ``ite then recurse (← applyLemma ``cost_ite)
    else if n == ``dite then recurse (← applyLemma ``cost_dite)
    else if n == ``letFun then changeTo ((m.getArg! 3).beta #[m.getArg! 2])
    else if n == ``modify then
      for g' in ← applyLemma ``cost_modify do
        unless ← g'.isAssigned do
          if (← instantiateMVars (← g'.getType)).isForall then
            setGoals [g']
            evalTactic (← `(tactic| intro s; run_simp; omega))
    else if n == ``loop then
      let lem ← match (m.getArg! 1).getAppFn with
        | .const ``closeCond _ => pure `Spqr.Fast.WalkM.cost_loop_close
        | .const ``firstIdxCond _ => pure ``cost_loop_firstIdx
        | .const ``sizeCond _ => pure ``cost_loop_size
        | _ => throwError "cost_auto: unknown loop condition{indentExpr m}"
      recurse (← applyLemma lem)
    else if let some lem := costLeaf n then recurse (← applyLemma lem)
    else if ← isMatcherApp m then
      if let .reduced m' ← reduceMatcher? m then return ← changeTo m'
      setGoals [g]
      evalTactic (← `(tactic| apply Cost.mono (w := max _ _) _ (Nat.le_refl _)))
      let gs ← getUnsolvedGoals
      let some gc ← gs.findM? (fun g' => do
          return (← instantiateMVars (← g'.getType)).isAppOf ``Cost)
        | throwError "cost_auto: lost the goal"
      setGoals [gc]
      evalTactic (← `(tactic| split))
      let [a1, a2] ← getUnsolvedGoals | throwError "cost_auto: expected a two-armed match"
      setGoals [a1]
      evalTactic (← `(tactic| apply Cost.mono _ (Nat.le_max_left _ _)))
      recurse (← getUnsolvedGoals)
      setGoals [a2]
      evalTactic (← `(tactic| apply Cost.mono _ (Nat.le_max_right _ _)))
      recurse (← getUnsolvedGoals)
    else
      setGoals [g]
      evalTactic (← `(tactic| assumption))

open Lean Elab Tactic in
elab "cost_auto" : tactic => do
  let g ← getMainGoal
  let rest ← getGoals
  costAuto 100000 g
  setGoals (← rest.filterM fun g => return !(← g.isAssigned))

theorem cost_maybeUnwrapNxt (t : NodeType) : Cost (maybeUnwrapNxt t) 1 := by
  apply Cost.mono
  · unfold maybeUnwrapNxt; cost_auto
  · omega

theorem cost_finishTstackTop (i : ItemId) : Cost (finishTstackTop i) 0 := by
  apply Cost.mono
  · unfold finishTstackTop; cost_auto
  · omega

theorem sz_maybeUnwrapNxt (t : NodeType) (s : WalkState) :
    s.tstack.size ≤ ((maybeUnwrapNxt t).run s).2.tstack.size := by
  simp only [maybeUnwrapNxt, allocItem, nxt, stackDir, getItem, modifyNxt]
  run_simp
  (repeat' split) <;> simp

theorem sz_finishTstackTop (i : ItemId) (s : WalkState) :
    s.tstack.size ≤ ((finishTstackTop i).run s).2.tstack.size := by
  simp only [run_simp, finishTstackTop, cur, stackDir, makeVs, modifyItem, modifyCur]
  omega

theorem spec_rest (t : NodeType) :
    Spec (fun s => 1 ≤ s.tstack.size)
      (do let item ← maybeUnwrapNxt t; mergeTstackTops; finishTstackTop item) (fun _ => True) 2 4 := by
  refine (spec_bind (a1 := 1) (b1 := 0) (a2 := 1) (b2 := 4)
    (spec_of_cost (cost_maybeUnwrapNxt _) (sz_maybeUnwrapNxt _)) fun i => ?_).weaken (by omega) (by omega)
  refine (spec_bind (a1 := 1) (b1 := 4) (a2 := 0) (b2 := 0) (spec_mergeTstackTops (k := 0))
    fun _ => ?_).weaken (by omega) (by omega)
  exact (spec_of_cost (k := 0) (cost_finishTstackTop _) (sz_finishTstackTop _)).conseq
    (fun _ _ => Nat.zero_le _) (fun _ _ => trivial)

/-- One iteration of loop 1 pays for itself: 3 ticks against at least one real pop. -/
theorem spec_closeBody (d : Nat) (e : Bool) :
    Spec (fun s => 2 ≤ s.tstack.size) (closeBody d e) (fun _ => True) 3 4 := by
  unfold closeBody
  simp only [pure_bind]
  refine (spec_bind (a1 := 0) (b1 := 0) (a2 := 3) (b2 := 4) (spec_of_cost cost_nxt sz_nxt)
    fun t => ?_).weaken (by omega) (by omega)
  refine spec_ite _ ?_ ?_
  · refine (spec_bind (a1 := 0) (b1 := 0) (a2 := 3) (b2 := 4) (spec_of_cost cost_nxt sz_nxt)
      fun t => ?_).weaken (by omega) (by omega)
    refine (spec_bind (a1 := 0) (b1 := 0) (a2 := 3) (b2 := 4)
      (spec_of_cost (cost_setStackDir _ _) (sz_setStackDir _ _)) fun _ => ?_).weaken (by omega) (by omega)
    exact (spec_bind (a1 := 1) (b1 := 4) (a2 := 2) (b2 := 4) (spec_mergeTstackTops (k := 1))
      fun _ => spec_rest _).weaken (by omega) (by omega)
  · refine (spec_bind (a1 := 0) (b1 := 0) (a2 := 3) (b2 := 4) (spec_of_cost cost_nxt sz_nxt)
      fun t => ?_).weaken (by omega) (by omega)
    refine (spec_bind (a1 := 0) (b1 := 0) (a2 := 3) (b2 := 4) (spec_of_cost cost_cur sz_cur)
      fun c => ?_).weaken (by omega) (by omega)
    refine spec_ite _ ?_ ?_ <;>
      exact ((spec_rest _).conseq (fun _ h => by omega) (fun _ _ => trivial)).weaken (by omega) (by omega)

theorem cost_loop_close (n d : Nat) (e : Bool) : Cost (loop n (closeCond d) (closeBody d e)) 0 :=
  cost_loop (closeCond_state d) (closeCond_true d) ((spec_closeBody d e).pays (by omega)) n

end WalkM

open Spqr.Fast.WalkM

/-- Finishing one out-edge costs a constant, amortized. -/
theorem cost_finishEdge (v d : Nat) (o : DfsOut) (orig : Nat) (hv : Bool) :
    Cost (finishEdge v d o orig hv) 40 := by
  apply Cost.mono
  · unfold finishEdge; cost_auto
  · omega

end Spqr.Fast

namespace Spqr

mutual
def DfsTree.size : DfsTree → Nat
  | .node _ outs => 1 + DfsOut.sizeList outs
def DfsOut.size : DfsOut → Nat
  | .back .. => 1
  | .tree _ _ child => 1 + child.size
def DfsOut.sizeList : List DfsOut → Nat
  | [] => 0
  | o :: rest => o.size + DfsOut.sizeList rest
end

theorem DfsTree.one_le_size (t : DfsTree) : 1 ≤ t.size := by
  cases t; simp only [DfsTree.size]; omega

mutual
theorem DfsTree.size_eq (t : DfsTree) : t.size = t.verts.length + t.edges.length := by
  match t with
  | .node v outs =>
    simp only [DfsTree.size, DfsTree.verts, DfsTree.edges, List.length_cons, DfsOut.sizeList_eq outs]
    omega
theorem DfsOut.sizeList_eq (outs : List DfsOut) :
    DfsOut.sizeList outs = (DfsOut.vertsList outs).length + (DfsOut.edgesList outs).length := by
  match outs with
  | [] => rfl
  | .back e a b :: rest =>
    simp only [DfsOut.sizeList, DfsOut.size, DfsOut.vertsList, DfsOut.edgesList, List.length_cons,
      DfsOut.sizeList_eq rest]
    omega
  | .tree e a child :: rest =>
    simp only [DfsOut.sizeList, DfsOut.size, DfsOut.vertsList, DfsOut.edgesList, List.length_cons,
      List.length_append, DfsOut.sizeList_eq rest, DfsTree.size_eq child]
    omega
end

end Spqr

namespace Spqr.Fast
open Spqr.Fast.WalkM

/-- Steps per DFS node / out-edge. -/
abbrev walkC : Nat := 46

mutual
theorem cost_walkTree (t : DfsTree) (d : Nat) : Cost (walkTree t d) (walkC * t.size) := by
  match t with
  | .node v outs =>
    have ih := cost_walkOuts v d outs false
    apply Cost.mono
    · unfold walkTree; cost_auto
    · simp only [DfsTree.size, walkC]; omega
theorem cost_walkOuts (v d : Nat) (outs : List DfsOut) (hv : Bool) :
    Cost (walkOuts v d outs hv) (walkC * DfsOut.sizeList outs) := by
  match outs with
  | [] => exact (cost_pure _).mono (by simp [DfsOut.sizeList])
  | o :: rest =>
    have ih1 := cost_walkOut v d o hv
    have ih2 : ∀ hv', Cost (walkOuts v d rest hv') (walkC * DfsOut.sizeList rest) :=
      fun hv' => cost_walkOuts v d rest hv'
    apply Cost.mono
    · unfold walkOuts; exact cost_bind ih1 ih2
    · simp only [DfsOut.sizeList, walkC]; omega
theorem cost_walkOut (v d : Nat) (o : DfsOut) (hv : Bool) : Cost (walkOut v d o hv) (walkC * o.size) := by
  cases o with
  | back e dest cls =>
    apply Cost.mono
    · unfold walkOut; cost_auto
    · simp only [DfsOut.size, walkC]; omega
  | tree e cls child =>
    have ih := cost_walkTree child (d + 1)
    apply Cost.mono
    · unfold walkOut; cost_auto
    · simp only [DfsOut.size, walkC]; omega
end

def forestSize (forest : List DfsTree) : Nat := (forest.map DfsTree.size).sum

theorem forestSize_eq (forest : List DfsTree) :
    forestSize forest = (forest.flatMap DfsTree.verts).length + (forest.flatMap DfsTree.edges).length := by
  induction forest with
  | nil => rfl
  | cons t rest ih =>
    simp only [forestSize, List.map_cons, List.sum_cons, List.flatMap_cons, List.length_append,
      DfsTree.size_eq] at *
    omega

theorem cost_walkForest (forest : List DfsTree) : Cost (walkForest forest) ((walkC + 1) * forestSize forest) := by
  induction forest with
  | nil => exact (cost_pure _).mono (by simp [forestSize])
  | cons t rest ih =>
    unfold walkForest at *
    simp only [List.forM]
    have h := cost_walkTree t 0
    have h1 := t.one_le_size
    refine (cost_bind (cost_bind h fun _ => cost_bind cost_popTstack fun _ => cost_modifyItem _ _)
      fun _ => ih).mono ?_
    simp only [forestSize, List.map_cons, List.sum_cons, walkC] at *
    omega

/-- Phase 2 step bound, in terms of the size of the forest it walks. -/
theorem walk_ticks_le_forest (g : Graph) (tern : Bool) (forest : List DfsTree) :
    (g.walkFast tern forest).ticks ≤ (walkC + 1) * forestSize forest := by
  have h := (cost_walkForest forest (WalkState.init g tern)).1
  unfold Graph.walkFast
  simp only [WalkState.pot, WalkState.init, Array.size_empty] at h ⊢
  omega

/-- The walk allocates at most one item per tick. -/
theorem walk_items_le (g : Graph) (tern : Bool) (forest : List DfsTree) :
    (g.walkFast tern forest).items.size ≤ 1 + g.nv + g.ne + (g.walkFast tern forest).ticks := by
  have h := (cost_walkForest forest (WalkState.init g tern)).2
  unfold Graph.walkFast
  simp [WalkState.Alloc, WalkState.init, initialItems] at h ⊢
  omega

/-- The DFS forest has exactly one node per vertex and one out-edge per edge
(`dfsForest_spanning'`), so `forestSize = V + E`. -/
theorem dfsForest_size_le (g : Graph) (hg : g.WF) (vertOrder edgeOrder : List Nat)
    (hvo : OrderOK g.nv vertOrder) (heo : OrderOK g.ne edgeOrder) :
    forestSize (g.dfsForest vertOrder edgeOrder) ≤ g.nv + g.ne := by
  obtain ⟨h1, h2⟩ := dfsForest_spanning' hg hvo heo
  exact Nat.le_of_eq (by rw [forestSize_eq, h1.length_eq, h2.length_eq, List.length_range, List.length_range])

/-- Phase 2 is `O(V + E)`: `47 * (V + E) + 0` steps. -/
theorem walk_ticks_le (g : Graph) (hg : g.WF) (tern : Bool) (vertOrder edgeOrder : List Nat)
    (hvo : OrderOK g.nv vertOrder) (heo : OrderOK g.ne edgeOrder) :
    (g.walkFast tern (g.dfsForest vertOrder edgeOrder)).ticks ≤ (walkC + 1) * (g.nv + g.ne) + 0 :=
  Nat.le_trans (walk_ticks_le_forest ..) (by
    have := dfsForest_size_le g hg vertOrder edgeOrder hvo heo
    simp only [walkC]; omega)

end Spqr.Fast
