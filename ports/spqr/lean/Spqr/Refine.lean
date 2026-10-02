import Lean
import Spqr.StateRun

/-!
# Simulation of one `StateT σ Id` program by another

`Sim proj rel mf m` says that running `mf` from any state `s` and running `m` from `proj s` end in
states related by `proj` and return values related by `rel`. The lemmas below let a simulation be
established by walking the syntax of the two programs in lockstep; `sim_auto` does that walk.
-/

namespace Spqr.Refine

variable {σ τ α β α' β' : Type}

def Sim (proj : σ → τ) (rel : α → β → Prop) (mf : StateT σ Id α) (m : StateT τ Id β) : Prop :=
  ∀ s, proj (mf.run s).2 = (m.run (proj s)).2 ∧ rel (mf.run s).1 (m.run (proj s)).1

namespace Sim
open StateRun

variable {proj : σ → τ} {rel : α → β → Prop}

theorem pure (a : α) (b : β) (h : rel a b) : Sim proj rel (pure a) (pure b) := fun _ =>
  ⟨rfl, h⟩

theorem bind {rel1 : α → β → Prop} {mf : StateT σ Id α} {m : StateT τ Id β}
    {ff : α → StateT σ Id α'} {f : β → StateT τ Id β'} {rel2 : α' → β' → Prop}
    (h1 : Sim proj rel1 mf m) (h2 : ∀ a b, rel1 a b → Sim proj rel2 (ff a) (f b)) :
    Sim proj rel2 (mf >>= ff) (m >>= f) := fun s => by
  rw [run_bind, runBind_eq, run_bind, runBind_eq]
  obtain ⟨hs, hr⟩ := h1 s
  have h := h2 _ _ hr (mf.run s).2
  rw [hs] at h
  exact h

theorem map {mf : StateT σ Id α} {m : StateT τ Id β} {gf : α → α'} {g : β → β'} {rel2 : α' → β' → Prop}
    (h1 : Sim proj rel mf m) (h2 : ∀ a b, rel a b → rel2 (gf a) (g b)) :
    Sim proj rel2 (gf <$> mf) (g <$> m) := fun s => by
  rw [run_map, runMap_eq, run_map, runMap_eq]
  exact ⟨(h1 s).1, h2 _ _ (h1 s).2⟩

theorem ite {c c' : Prop} [Decidable c] [Decidable c'] {tf ef : StateT σ Id α} {t e : StateT τ Id β}
    (hc : c ↔ c') (ht : Sim proj rel tf t) (he : Sim proj rel ef e) :
    Sim proj rel (if c then tf else ef) (if c' then t else e) := by
  by_cases h : c
  · have h' := hc.1 h; simp only [h, h', ↓reduceIte]; exact ht
  · have h' : ¬c' := fun h' => h (hc.2 h'); simp only [h, h', ↓reduceIte]; exact he

theorem dite {c c' : Prop} [Decidable c] [Decidable c'] {tf : c → StateT σ Id α} {ef : ¬c → StateT σ Id α}
    {t : c' → StateT τ Id β} {e : ¬c' → StateT τ Id β}
    (hc : c ↔ c') (ht : ∀ h h', Sim proj rel (tf h) (t h')) (he : ∀ h h', Sim proj rel (ef h) (e h')) :
    Sim proj rel (_root_.dite c tf ef) (_root_.dite c' t e) := by
  by_cases h : c
  · have h' := hc.1 h; simp only [h, h', ↓reduceDIte]; exact ht _ _
  · have h' : ¬c' := fun h' => h (hc.2 h'); simp only [h, h', ↓reduceDIte]; exact he _ _

/-- A fast-side action that does not affect the projection (a tick). -/
theorem skipL {a : StateT σ Id Unit} {ff : Unit → StateT σ Id α} {m : StateT τ Id β}
    (h : ∀ s, proj (a.run s).2 = proj s) (hk : Sim proj rel (ff ()) m) : Sim proj rel (a >>= ff) m := fun s => by
  rw [run_bind, runBind_eq]
  have := hk (a.run s).2
  rw [h] at this
  exact this

theorem get : Sim proj (fun a b => b = proj a) (get : StateT σ Id σ) (get : StateT τ Id τ) := fun _ =>
  ⟨rfl, rfl⟩

theorem modify {ff : σ → σ} {f : τ → τ} (h : ∀ s, proj (ff s) = f (proj s)) :
    Sim proj Eq (modify ff : StateT σ Id Unit) (modify f) := fun s =>
  ⟨h s, rfl⟩

theorem forM {γ : Type} (l : List γ) {ff : γ → StateT σ Id PUnit} {f : γ → StateT τ Id PUnit}
    (h : ∀ x, Sim proj Eq (ff x) (f x)) : Sim proj Eq (l.forM ff) (l.forM f) := by
  induction l with
  | nil => exact pure _ _ rfl
  | cons x l ih => exact bind (h x) fun _ _ _ => ih

theorem forIn {γ : Type} (l : List γ) (init : α) {ff : γ → α → StateT σ Id (ForInStep α)}
    {f : γ → α → StateT τ Id (ForInStep α)} (h : ∀ x b, Sim proj Eq (ff x b) (f x b)) :
    Sim proj Eq (forIn l init ff) (forIn l init f) := by
  induction l generalizing init with
  | nil => simp only [List.forIn_nil]; exact pure _ _ rfl
  | cons x l ih =>
    simp only [List.forIn_cons]
    refine bind (h x init) fun r _ hr => ?_
    subst hr
    cases r with
    | done b => exact pure _ _ rfl
    | yield b => exact ih b

end Sim

/-! ## `sim_auto` -/

/-- Leaf rule: closes `Sim proj rel mf m` for a pair of primitive actions. Extended per module with
`macro_rules`. -/
syntax "sim_leaf" : tactic
macro_rules | `(tactic| sim_leaf) => `(tactic| fail "sim_leaf: no rule")

/-- Discharges side conditions (`c ↔ c'`, `rel a b`, `proj (ff s) = f (proj s)`). -/
syntax "sim_side" : tactic
macro_rules
  | `(tactic| sim_side) =>
    `(tactic| (intros; (try subst_vars); first | rfl | exact Iff.rfl | trivial | (simp [sim_simp]; done)))

open Lean Elab Tactic Meta in
/-- Fast-side actions skipped with `Sim.skipL`. Extended per module. -/
initialize simSilentExt : SimplePersistentEnvExtension Name NameSet ←
  registerSimplePersistentEnvExtension {
    addEntryFn := fun s n => s.insert n
    addImportedFn := fun as => as.foldl (fun s a => a.foldl (·.insert ·) s) {} }

open Lean Elab Command in
elab "sim_silent " n:ident : command => do
  let n ← liftCoreM <| realizeGlobalConstNoOverloadWithInfo n
  modifyEnv fun env => simSilentExt.addEntry env n

register_option sim_auto.trace : Bool := { defValue := false, descr := "trace sim_auto steps" }

open Lean Elab Tactic Meta in
/-- Walks two `StateT σ Id` programs in lockstep, applying the `Sim` lemmas; leaves go to
`sim_leaf`, side conditions to `sim_side`. -/
def simAuto : Nat → MVarId → TacticM Unit
  | 0, _ => throwError "sim_auto: out of fuel"
  | fuel + 1, g => g.withContext do
    if ← g.isAssigned then return
    let ty := (← instantiateMVars (← g.getType)).consumeMData
    if ty.isForall then
      let (h, g') ← g.intro1
      let g' ← g'.withContext do
        let hty := (← instantiateMVars (← h.getType)).headBeta
        if hty.isAppOfArity ``Eq 3 then
          try pure (← subst g' h) catch _ => pure g'
        else pure g'
      let g' ← g'.withContext do
        try
          let gs ← Lean.Elab.Tactic.run g' (evalTactic (← `(tactic| simp only [sim_simp])))
          match gs with | [g''] => pure g'' | _ => pure g'
        catch _ => pure g'
      return ← simAuto fuel g'
    unless ty.isAppOfArity ``Sim 8 do
      throwError "sim_auto: not a Sim goal{indentExpr ty}"
    let isSimGoal (g' : MVarId) : TacticM Bool := do
      let t := (← instantiateMVars (← g'.getType)).consumeMData
      return t.getForallBody.isAppOf ``Sim
    let recurse (gs : List MVarId) : TacticM Unit := do
      for g' in gs do
        unless ← g'.isAssigned do
          if ← isSimGoal g' then simAuto fuel g'
          else
            setGoals [g']
            evalTactic (← `(tactic| sim_side))
    let tryTac (tac : TacticM Unit) : TacticM Bool := do
      let s ← saveState
      try
        tac
        return true
      catch _ =>
        s.restore
        return false
    -- local hypotheses (induction hypotheses, `Sim` assumptions)
    for d in (← getLCtx) do
      if d.isImplementationDetail then continue
      if (← instantiateMVars d.type).getForallBody.isAppOf ``Sim then
        if ← tryTac (do recurse (← withTransparency .instances <| g.apply (mkFVar d.fvarId))) then return
    let mf := (ty.getArg! 6).headBeta
    let m := (ty.getArg! 7).headBeta
    let changeTo (mf' m' : Expr) : TacticM Unit := do
      let args := ty.getAppArgs
      let ty' := mkAppN ty.getAppFn (args.set! 6 mf' |>.set! 7 m')
      simAuto fuel (← g.change ty')
    let applyLemma (n : Name) : TacticM Unit := do
      recurse (← withTransparency .instances <| g.apply (← mkConstWithFreshMVarLevels n))
    match mf, m with
    | .letE _ _ v b _, _ => return ← changeTo (b.instantiate1 v) m
    | _, .letE _ _ v b _ => return ← changeTo mf (b.instantiate1 v)
    | _, _ => pure ()
    let nf? := mf.getAppFn.constName?
    let ns? := m.getAppFn.constName?
    if sim_auto.trace.get (← getOptions) then logInfo m!"sim_auto: {nf?} vs {ns?}"
    if nf? == some ``letFun then return ← changeTo ((mf.getArg! 3).beta #[mf.getArg! 2]) m
    if ns? == some ``letFun then return ← changeTo mf ((m.getArg! 3).beta #[m.getArg! 2])
    if ← isMatcherApp mf then
      if let .reduced mf' ← reduceMatcher? mf then return ← changeTo mf' m
    if ← isMatcherApp m then
      if let .reduced m' ← reduceMatcher? m then return ← changeTo mf m'
    let tryLeaf : TacticM Bool := tryTac do
      setGoals [g]
      evalTactic (← `(tactic| sim_leaf))
      recurse (← getUnsolvedGoals)
    let isGet (n? : Option Name) : Bool :=
      n? == some ``MonadState.get || n? == some ``MonadStateOf.get || n? == some ``getThe
    if isGet nf? && isGet ns? then return ← applyLemma ``Sim.get
    -- an unconstrained result relation: a leaf may fix it, otherwise it is `Eq`
    let rel ← instantiateMVars (ty.getArg! 5)
    if rel.getAppFn.isMVar then
      if ← tryLeaf then return
      let α := ty.getArg! 2
      unless ← isDefEq rel (mkApp (mkConst ``Eq [Level.one]) α) do
        throwError "sim_auto: cannot default the result relation to Eq in{indentExpr ty}"
    if nf? == some ``Bind.bind then
      if let some na := (mf.getArg! 4).getAppFn.constName? then
        if (simSilentExt.getState (← getEnv)).contains na then
          return ← applyLemma ``Sim.skipL
    let both (n : Name) : Bool := nf? == some n && ns? == some n
    if both ``Bind.bind then return ← applyLemma ``Sim.bind
    if both ``Functor.map then return ← applyLemma ``Sim.map
    if both ``ite then return ← applyLemma ``Sim.ite
    if both ``dite then return ← applyLemma ``Sim.dite
    if both ``Pure.pure then return ← applyLemma ``Sim.pure
    if both ``modify then return ← applyLemma ``Sim.modify
    if both ``ForIn.forIn then return ← applyLemma ``Sim.forIn
    if both ``List.forM then return ← applyLemma ``Sim.forM
    -- a `match` on a variable: case on it, so that both sides reduce
    let casesOn (e : Expr) : TacticM Bool := do
      let some n := e.getAppFn.constName? | return false
      let some info ← getMatcherInfo? n | return false
      let args := e.getAppArgs
      for i in [0:info.numDiscrs] do
        let d := args[info.getFirstDiscrPos + i]!
        if let .fvar id := d then
          let subgoals ← g.cases id
          recurse (subgoals.map (·.mvarId)).toList
          return true
        if !d.hasLooseBVars && !d.hasMVar then
          if ← tryTac (do
              let (fvars, g') ← g.generalize #[{ expr := d }]
              let subgoals ← g'.cases fvars[0]!
              recurse (subgoals.map (·.mvarId)).toList) then
            return true
      return false
    if ← isMatcherApp mf then
      if ← casesOn mf then return
    if ← isMatcherApp m then
      if ← casesOn m then return
    if (← isMatcherApp mf) || (← isMatcherApp m) then
      setGoals [g]
      evalTactic (← `(tactic| split))
      return ← recurse (← getUnsolvedGoals)
    if ← tryLeaf then return
    if let some mf' ← delta? mf then return ← changeTo mf'.headBeta m
    if let some m' ← delta? m then return ← changeTo mf m'.headBeta
    throwError "sim_auto: no rule for{indentExpr mf}\nvs{indentExpr m}"

open Lean Elab Tactic in
elab "sim_auto" : tactic => do
  let g ← getMainGoal
  let rest ← getGoals
  simAuto 100000 g
  setGoals (← rest.filterM fun g => return !(← g.isAssigned))

end Spqr.Refine
