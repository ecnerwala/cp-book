import Lean
import Spqr.RunSimp

/-! Lemmas pushing `StateT.run` through the `StateT σ Id` combinators, and a tactic that
normalizes a computation applied to a concrete state with them. -/

namespace Spqr.StateRun

variable {σ α β : Type}

/-- `(m >>= f).run s` after `m` has run. Kept as a separate function so that `simp` normalizes
the pair before substituting it into `f` (substituting first would leave stale `Decidable`
instances inside `if`s, which `simp` does not rewrite). -/
def runBind (p : α × σ) (f : α → StateT σ Id β) : β × σ := (f p.1).run p.2
def runMap (p : α × σ) (f : α → β) : β × σ := (f p.1, p.2)

@[run_simp] theorem pure_id (a : α) : (pure a : Id α) = a := rfl
@[run_simp] theorem run_pure (a : α) (s : σ) : (pure a : StateT σ Id α).run s = (a, s) := rfl
@[run_simp] theorem run_bind (m : StateT σ Id α) (f : α → StateT σ Id β) (s : σ) :
    (m >>= f).run s = runBind (m.run s) f := rfl
@[run_simp] theorem run_map (f : α → β) (m : StateT σ Id α) (s : σ) :
    (f <$> m).run s = runMap (m.run s) f := rfl
@[run_simp] theorem runBind_mk (a : α) (s : σ) (f : α → StateT σ Id β) :
    runBind (a, s) f = (f a).run s := rfl
@[run_simp] theorem runMap_mk (a : α) (s : σ) (f : α → β) : runMap (a, s) f = (f a, s) := rfl
@[run_simp low] theorem runBind_eq (p : α × σ) (f : α → StateT σ Id β) :
    runBind p f = (f p.1).run p.2 := rfl
@[run_simp low] theorem runMap_eq (p : α × σ) (f : α → β) : runMap p f = (f p.1, p.2) := rfl
@[run_simp] theorem run_get (s : σ) : (get : StateT σ Id σ).run s = (s, s) := rfl
@[run_simp] theorem run_modify (f : σ → σ) (s : σ) : (modify f : StateT σ Id Unit).run s = ((), f s) := rfl
@[run_simp] theorem run_modifyGet (f : σ → α × σ) (s : σ) : (modifyGet f : StateT σ Id α).run s = f s := rfl
@[run_simp] theorem run_ite (c : Prop) [Decidable c] (m1 m2 : StateT σ Id α) (s : σ) :
    (if c then m1 else m2).run s = if c then m1.run s else m2.run s := by split <;> rfl
@[run_simp] theorem run_dite (c : Prop) [Decidable c] (m1 : c → StateT σ Id α) (m2 : ¬c → StateT σ Id α)
    (s : σ) : (dite c m1 m2).run s = dite c (fun h => (m1 h).run s) (fun h => (m2 h).run s) := by
  split <;> rfl

end Spqr.StateRun

/-- Normalize `(m : StateT σ Id α).run s`. -/
macro "run_simp" : tactic => `(tactic| simp only [run_simp])
