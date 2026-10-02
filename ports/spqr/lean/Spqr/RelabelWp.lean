import Spqr.RelabelGhost
import Spqr.ItemTree

/-!
# A `wp` calculus for the ghost `relabel`

Same shape as `Spqr.Fast.RelabelM.wp` (`RelabelCost.lean`): predicate transformers for the
`StateM` primitives, `wp_jp` to pass a `match` / `if` whose arms all continue into the same
(inlined) continuation, and `wp_forIn_inv` for the children loop.
-/

namespace Spqr.Ghost.RelabelM

variable {α β : Type}

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
theorem wp_item {Q : Item → RelabelState → Prop} (i : ItemId) : wp (item i) Q s = Q s.items[i]! s := rfl

theorem wp_mono {Q Q' : α → RelabelState → Prop} {m : RelabelM α} (h : wp m Q s)
    (himp : ∀ a s', Q a s' → Q' a s') : wp m Q' s := himp _ _ h

theorem orderedChildren_run_snd (it : Item) (n : Nat) (s : RelabelState) :
    ((orderedChildren it n).run s).2 = s := by
  unfold orderedChildren; split <;> rfl

/-- The children returned by `orderedChildren`, as a pure function of the state. -/
def orderedList (it : Item) (n : Nat) (s : RelabelState) : List ItemId :=
  if it.type ≠ .R then it.ch
  else it.ch.mergeSort fun a b =>
    (if a < 1 + s.g.nv then 2 * (s.vertPos[a - 1]! - n)
     else (s.vertPos[s.items[a]!.vs.1.getD 0]! - n) + (s.vertPos[s.items[a]!.vs.2.getD 0]! - n)) ≤
    (if b < 1 + s.g.nv then 2 * (s.vertPos[b - 1]! - n)
     else (s.vertPos[s.items[b]!.vs.1.getD 0]! - n) + (s.vertPos[s.items[b]!.vs.2.getD 0]! - n))

theorem orderedChildren_run_fst (it : Item) (n : Nat) (s : RelabelState) :
    ((orderedChildren it n).run s).1 = orderedList it n s := by
  unfold orderedChildren orderedList
  by_cases h : it.type = .R
  · simp [h]; rfl
  · simp [h]; rfl

theorem wp_orderedChildren {Q : List ItemId → RelabelState → Prop} (it : Item) (n : Nat)
    (h : Q (orderedList it n s) s) : wp (orderedChildren it n) Q s := by
  unfold wp; rw [orderedChildren_run_snd, orderedChildren_run_fst]; exact h

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

/-- Pass a `match` / `if` whose arms all end in the same continuation `k`: `m` run from `s` is `k`
run from `σ`. -/
theorem wp_jp {Q : Unit → RelabelState → Prop} (m k : RelabelM Unit) (σ s : RelabelState)
    (h : (m.run s).2 = (k.run σ).2) (hk : wp k Q σ) : wp m Q s := by
  unfold wp at *
  rw [h, show (m.run s).1 = (k.run σ).1 from rfl]
  exact hk

theorem arm_modify (f : RelabelState → RelabelState) (k : RelabelM Unit) (s : RelabelState) :
    ((modify f >>= fun _ => k).run s).2 = (k.run (f s)).2 := rfl

theorem arm_pure (k : RelabelM Unit) (s : RelabelState) :
    ((pure () >>= fun _ => k).run s).2 = (k.run s).2 := rfl

theorem arm_id (k : RelabelM Unit) (s : RelabelState) : (k.run s).2 = (k.run s).2 := rfl

end Spqr.Ghost.RelabelM

/-- Discharge one arm of a `match` / `if` against the shared continuation. -/
macro "jp_arm" : tactic =>
  `(tactic| first
      | exact Spqr.Ghost.RelabelM.arm_modify _ _ _
      | exact Spqr.Ghost.RelabelM.arm_pure _ _
      | exact Spqr.Ghost.RelabelM.arm_id _ _)
