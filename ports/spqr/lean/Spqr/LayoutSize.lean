import Lean
import Spqr.Relabel

/-!
# Array sizes of `layoutNode`

`layoutNode type node nvSt nvEn neSt neEn E` has `neEn - neSt` edges, `2 (nvEn - nvSt) + 1`
adjacency bounds and `2 (neEn - neSt)` adjacency slots, whatever the type and `E`: every write is a
`set!` / `modify` of `Layout.empty`.
-/

namespace Spqr

namespace Layout

/-- Array sizes of a `Layout` for `nV` node-verts and `nE` node-edges. -/
def Sz (nV nE : Nat) (l : Layout) : Prop :=
  l.edges.size = nE ∧ l.adjBounds.size = 2 * nV + 1 ∧ l.adjDat.size = 2 * nE

variable {nV nE : Nat} {l : Layout}

theorem Sz.empty (nV nE : Nat) : Sz nV nE (Layout.empty nV nE) := by simp [Sz, Layout.empty]

theorem Sz.setNe (h : Sz nV nE l) (neSt node ne : Nat) (nvs nds : Nat × Nat) :
    Sz nV nE (l.setNe neSt node ne nvs nds) := by
  simpa [Sz, Layout.setNe] using h

theorem Sz.mk_b (h : Sz nV nE l) {b : Array Nat} (hb : b.size = l.adjBounds.size) :
    Sz nV nE ⟨l.edges, b, l.adjDat⟩ := by
  obtain ⟨h1, h2, h3⟩ := h
  exact ⟨h1, hb.trans h2, h3⟩

end Layout

/-- A loop in `Id` preserves `P` if every step does. -/
theorem forIn_id_preserves {α β : Type} (P : β → Prop) (l : List α) (init : β)
    (f : α → β → Id (ForInStep β)) (h0 : P init)
    (hf : ∀ a b, P b → ∀ b', (f a b = ForInStep.yield b' ∨ f a b = ForInStep.done b') → P b') :
    P (forIn (m := Id) l init f) := by
  induction l generalizing init with
  | nil => exact h0
  | cons a l ih =>
    rw [List.forIn_cons]
    show P (match f a init with | .done b => b | .yield b => forIn l b f)
    cases hfa : f a init with
    | done b => exact hf a init h0 b (Or.inr hfa)
    | yield b => exact ih b (hf a init h0 b (Or.inl hfa))

theorem forIn_range_preserves {β : Type} (P : β → Prop) (r : Std.Legacy.Range) (init : β)
    (f : Nat → β → Id (ForInStep β)) (h0 : P init)
    (hf : ∀ a b, P b → ∀ b', (f a b = ForInStep.yield b' ∨ f a b = ForInStep.done b') → P b') :
    P (forIn (m := Id) r init f) := by
  rw [Std.Legacy.Range.forIn_eq_forIn_range']
  exact forIn_id_preserves P _ init f h0 hf

theorem forIn_id_preserves_fst {α β γ : Type} (P : β → Prop) (l : List α) (init : β × γ)
    (f : α → β × γ → Id (ForInStep (β × γ))) (h0 : P init.1)
    (hf : ∀ a b, P b.1 → ∀ b', (f a b = ForInStep.yield b' ∨ f a b = ForInStep.done b') → P b'.1) :
    P (forIn (m := Id) l init f).1 :=
  forIn_id_preserves (fun p => P p.1) l init f h0 hf

theorem forIn_range_preserves_fst {β γ : Type} (P : β → Prop) (r : Std.Legacy.Range) (init : β × γ)
    (f : Nat → β × γ → Id (ForInStep (β × γ))) (h0 : P init.1)
    (hf : ∀ a b, P b.1 → ∀ b', (f a b = ForInStep.yield b' ∨ f a b = ForInStep.done b') → P b'.1) :
    P (forIn (m := Id) r init f).1 :=
  forIn_range_preserves (fun p => P p.1) r init f h0 hf

theorem forIn_id_preserves_mfst {α β γ : Type} (P : β → Prop) (l : List α) (init : MProd β γ)
    (f : α → MProd β γ → Id (ForInStep (MProd β γ))) (h0 : P init.1)
    (hf : ∀ a b, P b.1 → ∀ b', (f a b = ForInStep.yield b' ∨ f a b = ForInStep.done b') → P b'.1) :
    P (forIn (m := Id) l init f).1 :=
  forIn_id_preserves (fun p => P p.1) l init f h0 hf

theorem forIn_range_preserves_mfst {β γ : Type} (P : β → Prop) (r : Std.Legacy.Range) (init : MProd β γ)
    (f : Nat → MProd β γ → Id (ForInStep (MProd β γ))) (h0 : P init.1)
    (hf : ∀ a b, P b.1 → ∀ b', (f a b = ForInStep.yield b' ∨ f a b = ForInStep.done b') → P b'.1) :
    P (forIn (m := Id) r init f).1 :=
  forIn_range_preserves (fun p => P p.1) r init f h0 hf

open Lean Elab Tactic Meta in
/-- One step towards a `Layout.Sz` goal, dispatching on the head symbol of the layout (so that no
unification ever has to evaluate a loop). -/
elab "sz_step" : tactic => do
  let g ← getMainGoal
  let t ← whnfR (← instantiateMVars (← g.getType))
  let some (nV, nE, l) := t.app3? ``Spqr.Layout.Sz | throwError "sz_step: not a Sz goal"
  let l ← whnfCore l
  let g ← g.replaceTargetDefEq (mkApp3 t.getAppFn nV nE l)
  replaceMainGoal [g]
  match l.getAppFn with
  | .fvar _ => evalTactic (← `(tactic| assumption))
  | .const n _ =>
    if n == ``Spqr.Layout.setNe then evalTactic (← `(tactic| apply Spqr.Layout.Sz.setNe))
    else if n == ``Spqr.Layout.mk then evalTactic (← `(tactic| apply Spqr.Layout.Sz.mk_b))
    else if n == ``Spqr.Layout.empty then evalTactic (← `(tactic| exact Spqr.Layout.Sz.empty _ _))
    else if n == ``ForIn.forIn then
      evalTactic (← `(tactic| first
        | refine Spqr.forIn_range_preserves (fun l => Spqr.Layout.Sz _ _ l) _ _ _ ?_ ?_
        | refine Spqr.forIn_id_preserves (fun l => Spqr.Layout.Sz _ _ l) _ _ _ ?_ ?_))
    else if n == ``MProd.fst || n == ``Prod.fst then
      evalTactic (← `(tactic| first
        | refine Spqr.forIn_range_preserves_mfst (fun l => Spqr.Layout.Sz _ _ l) _ _ _ ?_ ?_
        | refine Spqr.forIn_id_preserves_mfst (fun l => Spqr.Layout.Sz _ _ l) _ _ _ ?_ ?_
        | refine Spqr.forIn_range_preserves_fst (fun l => Spqr.Layout.Sz _ _ l) _ _ _ ?_ ?_
        | refine Spqr.forIn_id_preserves_fst (fun l => Spqr.Layout.Sz _ _ l) _ _ _ ?_ ?_))
    else throwError "sz_step: unknown head {n}"
  | _ => throwError "sz_step: unexpected layout term"

/-- Close a `Layout.Sz` goal built from `Layout.empty` by `setNe` / `set!` / `modify` and loops. -/
syntax "sz_auto" : tactic
macro_rules
  | `(tactic| sz_auto) => `(tactic| repeat' first
      | assumption
      | sz_step
      | (intro a b hb b' h; rcases h with h | h <;> first | cases h | (obtain ⟨a1, a2⟩ := a; cases h))
      | simp only [Array.size_set!, Array.size_modify])

theorem layoutNode_sz_F (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat)) :
    Layout.Sz (nvEn - nvSt) (neEn - neSt) (layoutNode .F node nvSt nvEn neSt neEn E) := by
  unfold layoutNode
  simp only [Id.run, bind, pure, Bool.or_eq_true, beq_iff_eq, reduceCtorEq, or_false, false_or, or_self,
    ↓reduceIte, or_true, true_or]
  try dsimp only
  repeat' split
  all_goals sz_auto
theorem layoutNode_sz_V (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat)) :
    Layout.Sz (nvEn - nvSt) (neEn - neSt) (layoutNode .V node nvSt nvEn neSt neEn E) := by
  unfold layoutNode
  simp only [Id.run, bind, pure, Bool.or_eq_true, beq_iff_eq, reduceCtorEq, or_false, false_or, or_self,
    ↓reduceIte, or_true, true_or]
  try dsimp only
  repeat' split
  all_goals sz_auto
theorem layoutNode_sz_Q (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat)) :
    Layout.Sz (nvEn - nvSt) (neEn - neSt) (layoutNode .Q node nvSt nvEn neSt neEn E) := by
  unfold layoutNode
  simp only [Id.run, bind, pure, Bool.or_eq_true, beq_iff_eq, reduceCtorEq, or_false, false_or, or_self,
    ↓reduceIte, or_true, true_or]
  try dsimp only
  repeat' split
  all_goals sz_auto
theorem layoutNode_sz_I (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat)) :
    Layout.Sz (nvEn - nvSt) (neEn - neSt) (layoutNode .I node nvSt nvEn neSt neEn E) := by
  unfold layoutNode
  simp only [Id.run, bind, pure, Bool.or_eq_true, beq_iff_eq, reduceCtorEq, or_false, false_or, or_self,
    ↓reduceIte, or_true, true_or]
  try dsimp only
  repeat' split
  all_goals sz_auto
theorem layoutNode_sz_O (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat)) :
    Layout.Sz (nvEn - nvSt) (neEn - neSt) (layoutNode .O node nvSt nvEn neSt neEn E) := by
  unfold layoutNode
  simp only [Id.run, bind, pure, Bool.or_eq_true, beq_iff_eq, reduceCtorEq, or_false, false_or, or_self,
    ↓reduceIte, or_true, true_or]
  try dsimp only
  repeat' split
  all_goals sz_auto
theorem layoutNode_sz_P (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat)) :
    Layout.Sz (nvEn - nvSt) (neEn - neSt) (layoutNode .P node nvSt nvEn neSt neEn E) := by
  unfold layoutNode
  simp only [Id.run, bind, pure, Bool.or_eq_true, beq_iff_eq, reduceCtorEq, or_false, false_or, or_self,
    ↓reduceIte, or_true, true_or]
  try dsimp only
  repeat' split
  all_goals sz_auto
theorem layoutNode_sz_S (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat)) :
    Layout.Sz (nvEn - nvSt) (neEn - neSt) (layoutNode .S node nvSt nvEn neSt neEn E) := by
  unfold layoutNode
  simp only [Id.run, bind, pure, Bool.or_eq_true, beq_iff_eq, reduceCtorEq, or_false, false_or, or_self,
    ↓reduceIte, or_true, true_or]
  try dsimp only
  repeat' split
  all_goals sz_auto
theorem layoutNode_sz_R (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat)) :
    Layout.Sz (nvEn - nvSt) (neEn - neSt) (layoutNode .R node nvSt nvEn neSt neEn E) := by
  unfold layoutNode
  simp only [Id.run, bind, pure, Bool.or_eq_true, beq_iff_eq, reduceCtorEq, or_false, false_or, or_self,
    ↓reduceIte, or_true, true_or]
  try dsimp only
  repeat' split
  all_goals sz_auto

theorem layoutNode_sz (type : NodeType) (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat)) :
    Layout.Sz (nvEn - nvSt) (neEn - neSt) (layoutNode type node nvSt nvEn neSt neEn E) := by
  cases type
  case F => exact layoutNode_sz_F ..
  case V => exact layoutNode_sz_V ..
  case Q => exact layoutNode_sz_Q ..
  case I => exact layoutNode_sz_I ..
  case O => exact layoutNode_sz_O ..
  case P => exact layoutNode_sz_P ..
  case S => exact layoutNode_sz_S ..
  case R => exact layoutNode_sz_R ..

end Spqr
