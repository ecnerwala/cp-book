import Lean
import Spqr.RelabelFast
import Spqr.WalkCost
import Spqr.ItemTree
import Spqr.Correctness
import Mathlib.Tactic.SplitIfs

/-!
# Step bound for phase 3

`relabel` ticks `1 + [cap edge] + [two own node-verts]` when it numbers a node and `2` per child
when it lays out the child list. Against the potential `4 * types.size + 2 * chDat.size` each
primitive action is non-increasing, so `ticks ≤ 4 * #nodes numbered + 2 * #child slots`. Turning
that into `O(V + E)` needs the (unproved, see `relabelRun_sizes_le`) fact that the item tree is a
tree, so that `relabel` visits each item at most once.
-/

namespace Spqr.Fast
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
theorem orderedChildren (it : Spqr.Item) (n : Nat) :
    ROk (RelabelM.orderedChildren it n) := fun s => by
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
theorem relabel_ticks_le_sizes (g : Graph) (items : Array Spqr.Item) :
    (relabelRun g items).ticks ≤
      4 * (relabelRun g items).types.size + 2 * (relabelRun g items).chDat.size := by
  have h := ok_relabel items.size rootItem none none none (RelabelState.init g items)
  unfold relabelRun
  simp only [RelabelState.pot, RelabelState.init, Spqr.RelabelState.init, List.size_toArray,
    List.length_nil] at h ⊢
  omega

namespace RelabelM

/-- Weakest precondition: `Q` holds of the result and final state of `m` run from `s`. -/
def wp (m : RelabelM α) (Q : α → RelabelState → Prop) (s : RelabelState) : Prop :=
  Q (m.run s).1 (m.run s).2

variable {s : RelabelState}

theorem wp_pure {Q : α → RelabelState → Prop} (a : α) : wp (pure a : RelabelM α) Q s = Q a s := rfl
theorem wp_bind {Q : α → RelabelState → Prop} (m : RelabelM β) (f : β → RelabelM α) :
    wp (m >>= f) Q s = wp m (fun b s' => wp (f b) Q s') s := rfl
theorem wp_map {Q : α → RelabelState → Prop} (f : β → α) (m : RelabelM β) :
    wp (f <$> m) Q s = wp m (fun b s' => Q (f b) s') s := rfl
theorem wp_get {Q : RelabelState → RelabelState → Prop} : wp (get : RelabelM RelabelState) Q s = Q s s := rfl
theorem wp_modify {Q : Unit → RelabelState → Prop} (f : RelabelState → RelabelState) :
    wp (modify f : RelabelM Unit) Q s = Q () (f s) := rfl
theorem wp_ite {Q : α → RelabelState → Prop} (c : Prop) [Decidable c] (m₁ m₂ : RelabelM α) :
    wp (if c then m₁ else m₂) Q s = if c then wp m₁ Q s else wp m₂ Q s := by split <;> rfl
theorem wp_item {Q : Spqr.Item → RelabelState → Prop} (i : ItemId) : wp (item i) Q s = Q s.items[i]! s := rfl

theorem orderedChildren_run_snd (it : Spqr.Item) (n : Nat) (s : RelabelState) :
    ((orderedChildren it n).run s).2 = s := by
  unfold orderedChildren; split <;> rfl

theorem orderedChildren_run_fst_perm (it : Spqr.Item) (n : Nat) (s : RelabelState) :
    ((orderedChildren it n).run s).1.Perm it.ch := by
  unfold orderedChildren; split
  · exact List.Perm.refl _
  · show (List.map (fun x : Nat × ItemId => x.2) (List.mergeSort (it.ch.map _) _)).Perm it.ch
    refine ((List.mergeSort_perm _ _).map _).trans ?_
    simp [List.map_map, Function.comp_def]

theorem wp_orderedChildren {Q : List ItemId → RelabelState → Prop} (it : Spqr.Item) (n : Nat)
    (h : ∀ l, l.Perm it.ch → Q l s) : wp (orderedChildren it n) Q s := by
  unfold wp; rw [orderedChildren_run_snd]; exact h _ (orderedChildren_run_fst_perm it n s)

theorem wp_forIn_inv {α β : Type} (l : List α) (init : β) (f : α → β → RelabelM (ForInStep β))
    (Q : β → RelabelState → Prop) (I : List α → β → RelabelState → Prop) (s : RelabelState)
    (hinit : I l init s)
    (hstep : ∀ a ∈ l, ∀ rest b s', I (a :: rest) b s' →
      wp (f a b) (fun r s'' => ∃ b', r = .yield b' ∧ I rest b' s'') s')
    (hfin : ∀ b s', I [] b s' → Q b s') : wp (forIn l init f) Q s := by
  induction l generalizing init s with
  | nil => exact hfin _ _ hinit
  | cons a l ih =>
    obtain ⟨b', hr, hI⟩ := hstep a (List.mem_cons_self ..) l init s hinit
    have e : (forIn (a :: l) init f).run s = (forIn l b' f).run ((f a init).run s).2 := by
      rw [List.forIn_cons]
      show (match ((f a init).run s).1 with
        | ForInStep.done b => pure b | ForInStep.yield b => forIn l b f).run ((f a init).run s).2 = _
      rw [hr]
    unfold wp; rw [e]
    exact ih b' _ hI fun a' ha' => hstep a' (List.mem_cons_of_mem _ ha')

/-- Package the inductive hypothesis for one child in the shape of the loop invariant. -/
theorem bound_step {items : Spqr.Items} (m : RelabelM Unit) (k : Nat)
    (hm : ∀ s : RelabelState, s.items = items → (m.run s).2.items = items ∧
      (m.run s).2.types.size ≤ s.types.size + k ∧ (m.run s).2.chDat.size + 1 ≤ s.chDat.size + k)
    (s : RelabelState) (hs : s.items = items) (r T C : Nat) (h1 : s.types.size + k + r ≤ T)
    (h2 : s.chDat.size + k + r ≤ C + 1) :
    (m.run s).2.items = items ∧ (m.run s).2.types.size + r ≤ T ∧ (m.run s).2.chDat.size + r ≤ C := by
  obtain ⟨e1, e2, e3⟩ := hm s hs
  exact ⟨e1, by omega, by omega⟩

/-- `s'` agrees with `s` on everything the size bound tracks. -/
abbrev Frame (s s' : RelabelState) : Prop :=
  s'.items = s.items ∧ s'.types.size = s.types.size ∧ s'.chDat.size = s.chDat.size

theorem wp_frame {Q : Unit → RelabelState → Prop} (m : RelabelM Unit)
    (hm : ∀ s, Frame s (m.run s).2) (hk : ∀ s', Frame s s' → Q () s') : wp m Q s :=
  hk _ (hm s)

theorem wp_unit {Q : Unit → RelabelState → Prop} (m : RelabelM Unit)
    (hk : ∀ s', s'.items = (m.run s).2.items → s'.types.size = (m.run s).2.types.size →
      s'.chDat.size = (m.run s).2.chDat.size → Q () s') : wp m Q s :=
  hk _ rfl rfl rfl

theorem wp_jp {Q : Unit → RelabelState → Prop} (m k : RelabelM Unit) (s : RelabelState)
    (h : ∀ s, ∃ s', Frame s s' ∧ (m.run s).2 = (k.run s').2)
    (hk : ∀ s', Frame s s' → wp k Q s') : wp m Q s := by
  obtain ⟨s', hf, he⟩ := h s
  have := hk s' hf
  unfold wp at this ⊢
  exact he ▸ this

/-- A `match`/`if` arm that runs a frame-preserving `modify` before the shared continuation. -/
theorem arm_modify (f : RelabelState → RelabelState) (k : RelabelM Unit) (s : RelabelState)
    (hf : Frame s (f s)) :
    ∃ s', Frame s s' ∧ ((modify f >>= fun _ => k).run s).2 = (k.run s').2 :=
  ⟨f s, hf, rfl⟩

/-- An arm that is just the shared continuation. -/
theorem arm_id (k : RelabelM Unit) (s : RelabelState) :
    ∃ s', Frame s s' ∧ (k.run s).2 = (k.run s').2 :=
  ⟨s, ⟨rfl, rfl, rfl⟩, rfl⟩

end RelabelM

/-- Discharge one arm of a `match`/`if`: a frame-preserving `modify` before the shared
continuation, or the continuation itself. -/
macro "jp_arm" : tactic =>
  `(tactic| first
      | exact RelabelM.arm_modify _ _ _ ⟨rfl, rfl, rfl⟩
      | exact RelabelM.arm_id _ _)

/-- `relabel fuel cur` numbers at most `|desc cur|` items and writes at most `|desc cur| - 1`
child slots, for `items` a rooted tree. States are kept abstract (`wp_frame`/`wp_unit`/`wp_jp`)
so the proof is linear in the size of `relabel`'s body. -/
theorem relabel_sizes {g : Graph} {items : Spqr.Items} (ht : items.Tree g) :
    ∀ (fuel : Nat) (cur : ItemId) (p pn ct : Option Nat) (s : RelabelState), cur < items.size →
      s.items = items →
      ((relabel fuel cur p pn ct).run s).2.items = items ∧
      ((relabel fuel cur p pn ct).run s).2.types.size ≤ s.types.size + (items.desc cur).card ∧
      ((relabel fuel cur p pn ct).run s).2.chDat.size + 1 ≤ s.chDat.size + (items.desc cur).card
  | 0, cur, p, pn, ct, s, hcur, hs => by
    have := Items.one_le_card_desc hcur
    show s.items = items ∧ s.types.size ≤ s.types.size + _ ∧ s.chDat.size + 1 ≤ s.chDat.size + _
    exact ⟨hs, by omega, by omega⟩
  | fuel + 1, cur, p, pn, ct, s, hcur, hs => by
    have ih := relabel_sizes ht fuel
    have hsum := ht.sum_card_desc_children hcur
    show RelabelM.wp (relabel (fuel + 1) cur p pn ct) (fun _ s' => s'.items = items ∧
      s'.types.size ≤ s.types.size + (items.desc cur).card ∧
      s'.chDat.size + 1 ≤ s.chDat.size + (items.desc cur).card) s
    unfold relabel
    simp only [RelabelM.wp_bind, RelabelM.wp_get, RelabelM.wp_item, hs, Items.getElem!_ch hcur]
    -- number `cur`
    refine RelabelM.wp_unit _ fun s₁ hi₁ ht₁ hc₁ => ?_
    simp only [run_simp, Array.size_push] at hi₁ ht₁ hc₁
    try dsimp only
    -- V/Q bookkeeping
    refine RelabelM.wp_jp _ _ _ (fun s => by split <;> jp_arm)
      fun s₂ ⟨hi₂, ht₂, hc₂⟩ => ?_
    simp only [RelabelM.wp_bind, RelabelM.wp_get]
    -- node-verts
    refine RelabelM.wp_frame _ (fun s => ⟨rfl, rfl, rfl⟩) fun s₃ ⟨hi₃, ht₃, hc₃⟩ => ?_
    try dsimp only
    -- R: vertex positions
    refine RelabelM.wp_jp _ _ _ (fun s => by split <;> jp_arm)
      fun s₄ ⟨hi₄, ht₄, hc₄⟩ => ?_
    simp only [RelabelM.wp_bind, RelabelM.wp_get]
    refine RelabelM.wp_orderedChildren _ _ fun children hperm => ?_
    rw [Items.getElem!_ch hcur] at hperm
    have hsum' : (children.map fun c => (items.desc c).card).sum =
        ((items.ch cur).map fun c => (items.desc c).card).sum := (hperm.map _).sum_eq
    -- child slots
    refine RelabelM.wp_unit _ fun s₅ hi₅ ht₅ hc₅ => ?_
    simp only [run_simp, Array.size_append, List.size_toArray] at hi₅ ht₅ hc₅
    try dsimp only
    -- node edges / bounds
    refine RelabelM.wp_frame _ (fun s => ⟨rfl, rfl, rfl⟩) fun s₆ ⟨hi₆, ht₆, hc₆⟩ => ?_
    try dsimp only
    -- cap twin: the continuation differs between the branches (`curNe` starts at `neSt + 1`)
    split
    all_goals
      refine RelabelM.wp_jp _ _ _ (fun s => by jp_arm) fun s₇ ⟨hi₇, ht₇, hc₇⟩ => ?_
      simp only [RelabelM.wp_bind]
      -- children
      refine RelabelM.wp_forIn_inv _ _ _ _ (fun rest b s' => s'.items = items ∧
        s'.types.size + (rest.map fun c => (items.desc c).card).sum ≤
          s.types.size + 1 + (children.map fun c => (items.desc c).card).sum ∧
        s'.chDat.size + (rest.map fun c => (items.desc c).card).sum ≤
          s.chDat.size + (children.map fun c => (items.desc c).card).sum + rest.length) _
        ?_ ?_ ?_
      · exact ⟨by rw [hi₇, hi₆, hi₅, hi₄, hi₃, hi₂, hi₁, hs], by omega, by omega⟩
      · intro c hc rest b s' ⟨hs', h1, h2⟩
        have hclt : c < items.size := ht.ch_lt cur c (hperm.subset hc)
        simp only [RelabelM.wp_bind, RelabelM.wp_get, RelabelM.wp_modify, RelabelM.wp_pure,
          RelabelM.wp_ite]
        split_ifs
        all_goals
          unfold RelabelM.wp
          try dsimp only
          refine ⟨_, rfl, RelabelM.bound_step _ _ (fun t ht' => ih c _ _ _ t hclt ht') _ ?_ _ _ _
            ?_ ?_⟩
          · exact hs'
          all_goals
            try simp only [Array.size_set!]
            simp only [List.map_cons, List.sum_cons, List.length_cons] at h1 h2
            omega
      · intro b s' ⟨hs', h1, h2⟩
        simp only [List.map_nil, List.sum_nil, List.length_nil, Nat.add_zero] at h1 h2
        try dsimp only
        refine RelabelM.wp_frame _ (fun s => ⟨rfl, rfl, rfl⟩) fun s₈ ⟨hi₈, ht₈, hc₈⟩ => ?_
        exact ⟨by rw [hi₈, hs'], by omega, by omega⟩

/-- `relabel` numbers at most `items.size` nodes and writes at most `items.size` child slots when
`items` is a rooted tree (`Items.Tree`). -/
theorem relabelRun_sizes_le (g : Graph) (items : Spqr.Items) (ht : items.Tree g) :
    (relabelRun g items).types.size ≤ items.size ∧ (relabelRun g items).chDat.size ≤ items.size := by
  have h0 : rootItem < items.size := Nat.lt_of_lt_of_le (by unfold rootItem; omega) ht.size
  obtain ⟨-, h1, h2⟩ :=
    relabel_sizes ht items.size rootItem none none none (RelabelState.init g items) h0 rfl
  have hc := @Items.card_desc_le items rootItem
  have e1 : (RelabelState.init g items).types.size = 0 := rfl
  have e2 : (RelabelState.init g items).chDat.size = 0 := rfl
  unfold relabelRun
  omega

/-- Steps per vertex / edge in phase 3: `6` per item, `items ≤ 1 + (V + E) + walk ticks`. -/
abbrev relabelC : Nat := 6 * (walkC + 2)

/-- Phase 3 is `O(V + E)`: `288 * (V + E) + 6` steps, given that the walk's items form a rooted
tree. `ht` is discharged by `walk_items_wf` (Correctness.lean, admitted upstream) via
`Items.WF.tree`, transported along `Graph.walkFast_items` (see `relabel_ticks_le'`);
`walk_typing` alone only gives the typing fields (`size/root/vert/edge/node/v_children`) of
`Items.Tree`, not the acyclicity / unique-parent structure the traversal bound needs. -/
theorem relabel_ticks_le (g : Graph) (hg : g.WF) (tern : Bool) (vertOrder edgeOrder : List Nat)
    (hvo : OrderOK g.nv vertOrder) (heo : OrderOK g.ne edgeOrder)
    (ht : Items.Tree g ((g.walkFast tern (g.dfsForest vertOrder edgeOrder)).items.map Item.toSlow)) :
    (relabelRun g ((g.walkFast tern (g.dfsForest vertOrder edgeOrder)).items.map Item.toSlow)).ticks ≤
      relabelC * (g.nv + g.ne) + 6 := by
  have h1 := relabel_ticks_le_sizes g ((g.walkFast tern (g.dfsForest vertOrder edgeOrder)).items.map Item.toSlow)
  have h2 := relabelRun_sizes_le g _ ht
  have h3 := walk_items_le g tern (g.dfsForest vertOrder edgeOrder)
  have h4 := walk_ticks_le g hg tern vertOrder edgeOrder hvo heo
  simp only [relabelC, walkC, Array.size_map] at *
  omega

/-- `relabel_ticks_le` with `ht` discharged by the admitted `walk_items_wf`; hence `sorryAx`. -/
theorem relabel_ticks_le' (g : Graph) (hg : g.WF) (tern : Bool) (vertOrder edgeOrder : List Nat)
    (hvo : OrderOK g.nv vertOrder) (heo : OrderOK g.ne edgeOrder) :
    (relabelRun g ((g.walkFast tern (g.dfsForest vertOrder edgeOrder)).items.map Item.toSlow)).ticks ≤
      relabelC * (g.nv + g.ne) + 6 :=
  relabel_ticks_le g hg tern vertOrder edgeOrder hvo heo
    (by rw [Graph.walkFast_items]; exact (walk_items_wf g tern vertOrder edgeOrder).tree)

end Spqr.Fast
