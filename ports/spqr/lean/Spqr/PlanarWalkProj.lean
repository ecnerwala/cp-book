import Spqr.PlanarWalk
import Lean

/-!
# The planar walk projects to the ordinary walk

`Lifts f m m'` says that the planar computation `m` is the ordinary computation `m'` on the
base component (returning `f` of its result), and `AuxOnly m` that `m` leaves the base component
alone. Every planar step is either a lifted base step or an auxiliary step, so `planarFinishEdge`
lifts `finishEdge` and `planarWalkForest` lifts `walkForest`.
-/

namespace Spqr

def Lifts (m : PlanarWalkM α) (m' : WalkM α) : Prop :=
  ∀ s, m' s.base = ((m s).1, (m s).2.base)

def AuxOnly (m : PlanarWalkM α) : Prop := ∀ s, (m s).2.base = s.base

namespace Lifts

theorem liftW (m : WalkM α) : Lifts (PlanarWalkM.liftW m) m := fun _ => rfl

theorem pure (a : α) : Lifts (Pure.pure a) (Pure.pure a) := fun _ => rfl

theorem bind {m : PlanarWalkM α} {m' : WalkM α} {k : α → PlanarWalkM γ} {k' : α → WalkM γ}
    (h : Lifts m m') (hk : ∀ a, Lifts (k a) (k' a)) : Lifts (m >>= k) (m' >>= k') := by
  intro s
  show k' (m' s.base).1 (m' s.base).2 = _
  rw [h s]
  exact hk _ _

theorem auxBind {m : PlanarWalkM α} {k : α → PlanarWalkM γ} {m' : WalkM γ}
    (h : AuxOnly m) (hk : ∀ a, Lifts (k a) m') : Lifts (m >>= k) m' := by
  intro s
  show m' s.base = ((k (m s).1 (m s).2).1, (k (m s).1 (m s).2).2.base)
  rw [← h s]
  exact hk _ _

theorem auxThen {m : PlanarWalkM Unit} {k : PlanarWalkM γ} {m' : WalkM γ}
    (h : AuxOnly m) (hk : Lifts k m') : Lifts (m *> k) m' := by
  intro s
  show m' s.base = ((k (m s).2).1, (k (m s).2).2.base)
  rw [← h s]
  exact hk _

theorem ite (c : Prop) [Decidable c] {t e : PlanarWalkM α} {t' e' : WalkM α}
    (ht : c → Lifts t t') (he : ¬c → Lifts e e') :
    Lifts (if c then t else e) (if c then t' else e') := by
  split
  · exact ht ‹_›
  · exact he ‹_›


end Lifts

namespace AuxOnly

theorem pure (a : α) : AuxOnly (Pure.pure a : PlanarWalkM α) := fun _ => rfl

theorem bind {m : PlanarWalkM α} {k : α → PlanarWalkM β}
    (h : AuxOnly m) (hk : ∀ a, AuxOnly (k a)) : AuxOnly (m >>= k) := by
  intro s
  show (k (m s).1 (m s).2).2.base = s.base
  rw [hk _ _, h s]

theorem ite (c : Prop) [Decidable c] {t e : PlanarWalkM α}
    (ht : c → AuxOnly t) (he : ¬c → AuxOnly e) : AuxOnly (if c then t else e) := by
  split
  · exact ht ‹_›
  · exact he ‹_›

theorem modifyAux (f : PlanarAux → PlanarAux) : AuxOnly (PlanarWalkM.modifyAux f) := fun _ => rfl
theorem getAux : AuxOnly PlanarWalkM.getAux := fun _ => rfl
theorem liftQem (m : QemM α) : AuxOnly (PlanarWalkM.liftQem m) := fun _ => rfl
theorem get : AuxOnly (PlanarWalkM.liftW get) := fun _ => rfl
theorem cur : AuxOnly (PlanarWalkM.liftW WalkM.cur) := fun _ => rfl
theorem nxt : AuxOnly (PlanarWalkM.liftW WalkM.nxt) := fun _ => rfl
theorem tstackSize : AuxOnly (PlanarWalkM.liftW WalkM.tstackSize) := fun _ => rfl
theorem stackDir (d : Nat) : AuxOnly (PlanarWalkM.liftW (WalkM.stackDir d)) := fun _ => rfl
theorem getItem (i : ItemId) : AuxOnly (PlanarWalkM.liftW (WalkM.getItem i)) := fun _ => rfl
theorem curPl : AuxOnly PlanarWalkM.curPl := fun _ => rfl
theorem nxtPl : AuxOnly PlanarWalkM.nxtPl := fun _ => rfl
theorem modifyCurPl (f) : AuxOnly (PlanarWalkM.modifyCurPl f) := fun _ => rfl
theorem modifyNxtPl (f) : AuxOnly (PlanarWalkM.modifyNxtPl f) := fun _ => rfl
theorem pushPl (e) : AuxOnly (PlanarWalkM.pushPl e) := fun _ => rfl
theorem popPl : AuxOnly PlanarWalkM.popPl := fun _ => rfl
theorem plAt (i) : AuxOnly (PlanarWalkM.plAt i) := fun _ => rfl
theorem modifyPlAt (i f) : AuxOnly (PlanarWalkM.modifyPlAt i f) := fun _ => rfl
theorem topDepthAt (i) : AuxOnly (PlanarWalkM.topDepthAt i) := fun _ => rfl
theorem nodeIdx (i) : AuxOnly (PlanarWalkM.nodeIdx i) := fun _ => rfl
theorem setNodePlanarity (i p) : AuxOnly (PlanarWalkM.setNodePlanarity i p) := fun _ => rfl
theorem getNodePlanarity (i) : AuxOnly (PlanarWalkM.getNodePlanarity i) := fun _ => rfl
theorem getItemFlips (i) : AuxOnly (PlanarWalkM.getItemFlips i) := fun _ => rfl
theorem modifyItemFlips (i f) : AuxOnly (PlanarWalkM.modifyItemFlips i f) := fun _ => rfl
theorem setItemFlips (i f) : AuxOnly (PlanarWalkM.setItemFlips i f) := fun _ => rfl
theorem makeEdgePlanarity (i t b) : AuxOnly (PlanarWalkM.makeEdgePlanarity i t b) := fun _ => rfl

end AuxOnly

end Spqr

namespace Spqr

theorem AuxOnly.forM (l : List α) (f : α → PlanarWalkM Unit) (h : ∀ x, AuxOnly (f x)) : AuxOnly (l.forM f) := by
  induction l with
  | nil => exact AuxOnly.pure _
  | cons x rest ih => exact AuxOnly.bind (h x) fun _ => ih

theorem Lifts.fmap {m : PlanarWalkM α} {m' : WalkM α} (h : Lifts m m') (g : α → β) :
    Lifts (g <$> m) (g <$> m') := fun s => by
  show (g (m' s.base).1, (m' s.base).2) = (g (m s).1, (m s).2.base)
  rw [h s]

theorem Lifts.forM (l : List α) {f : α → PlanarWalkM Unit} {f' : α → WalkM Unit}
    (h : ∀ x, Lifts (f x) (f' x)) : Lifts (l.forM f) (l.forM f') := by
  induction l with
  | nil => exact Lifts.pure _
  | cons x rest ih => exact Lifts.bind (h x) fun _ => ih

syntax "aux_only0" : tactic
macro_rules
  | `(tactic| aux_only0) => `(tactic| first
    | with_reducible exact AuxOnly.pure _
    | with_reducible first
      | exact AuxOnly.modifyAux _ | exact AuxOnly.liftQem _ | exact AuxOnly.get | exact AuxOnly.cur
      | exact AuxOnly.nxt | exact AuxOnly.tstackSize | exact AuxOnly.stackDir _ | exact AuxOnly.getItem _
      | exact AuxOnly.getAux | exact AuxOnly.curPl | exact AuxOnly.nxtPl | exact AuxOnly.modifyCurPl _
      | exact AuxOnly.modifyNxtPl _ | exact AuxOnly.pushPl _ | exact AuxOnly.popPl | exact AuxOnly.plAt _
      | exact AuxOnly.modifyPlAt _ _ | exact AuxOnly.topDepthAt _ | exact AuxOnly.nodeIdx _
      | exact AuxOnly.setNodePlanarity _ _ | exact AuxOnly.getNodePlanarity _ | exact AuxOnly.getItemFlips _
      | exact AuxOnly.modifyItemFlips _ _ | exact AuxOnly.setItemFlips _ _ | exact AuxOnly.makeEdgePlanarity _ _ _
    | (with_reducible refine AuxOnly.forM _ _ fun _ => ?_; aux_only0)
    | (with_reducible refine AuxOnly.bind ?_ fun _ => ?_ <;> aux_only0)
    | (with_reducible refine AuxOnly.ite _ (fun _ => ?_) (fun _ => ?_) <;> aux_only0)
    | (split <;> aux_only0)
    | (dsimp only; aux_only0))

namespace PlanarWalkM

theorem auxOnly_closeBackedges : AuxOnly closeBackedges := by unfold closeBackedges; aux_only0
theorem auxOnly_flipBeforeMerge (d fo b) : AuxOnly (flipBeforeMerge d fo b) := by unfold flipBeforeMerge; aux_only0
theorem auxOnly_pruneBackedges (d) : AuxOnly (pruneBackedges d) := by unfold pruneBackedges; aux_only0
theorem auxOnly_flipForLowval (l o) : AuxOnly (flipForLowval l o) := by unfold flipForLowval; aux_only0
theorem auxOnly_foldSides (d l) : AuxOnly (foldSides d l) := by unfold foldSides; aux_only0

end PlanarWalkM

syntax "aux_only" : tactic
macro_rules
  | `(tactic| aux_only) => `(tactic| first
    | with_reducible exact AuxOnly.pure _
    | with_reducible first
      | exact AuxOnly.modifyAux _ | exact AuxOnly.liftQem _ | exact AuxOnly.get | exact AuxOnly.cur
      | exact AuxOnly.nxt | exact AuxOnly.tstackSize | exact AuxOnly.stackDir _ | exact AuxOnly.getItem _
      | exact AuxOnly.getAux | exact AuxOnly.curPl | exact AuxOnly.nxtPl | exact AuxOnly.modifyCurPl _
      | exact AuxOnly.modifyNxtPl _ | exact AuxOnly.pushPl _ | exact AuxOnly.popPl | exact AuxOnly.plAt _
      | exact AuxOnly.modifyPlAt _ _ | exact AuxOnly.topDepthAt _ | exact AuxOnly.nodeIdx _
      | exact AuxOnly.setNodePlanarity _ _ | exact AuxOnly.getNodePlanarity _ | exact AuxOnly.getItemFlips _
      | exact AuxOnly.modifyItemFlips _ _ | exact AuxOnly.setItemFlips _ _ | exact AuxOnly.makeEdgePlanarity _ _ _
    | with_reducible first
      | exact PlanarWalkM.auxOnly_closeBackedges | exact PlanarWalkM.auxOnly_flipBeforeMerge _ _ _
      | exact PlanarWalkM.auxOnly_pruneBackedges _ | exact PlanarWalkM.auxOnly_flipForLowval _ _
      | exact PlanarWalkM.auxOnly_foldSides _ _
    | (with_reducible refine AuxOnly.forM _ _ fun _ => ?_; aux_only)
    | (with_reducible refine AuxOnly.bind ?_ fun _ => ?_ <;> aux_only)
    | (with_reducible refine AuxOnly.ite _ (fun _ => ?_) (fun _ => ?_) <;> aux_only)
    | (split <;> aux_only)
    | (dsimp only; aux_only))

end Spqr

namespace Spqr

theorem Lifts.bindAux {m : PlanarWalkM Unit} {m' : WalkM Unit} (h : Lifts m m')
    {k : Unit → PlanarWalkM Unit} (hk : ∀ a, AuxOnly (k a)) : Lifts (m >>= k) m' := by
  intro s
  show m' s.base = ((k (m s).1 (m s).2).1, (k (m s).1 (m s).2).2.base)
  rw [hk _ _, h s]

theorem Lifts.loop (fuel : Nat) {c : PlanarWalkM Bool} {c' : WalkM Bool} {b : PlanarWalkM Unit} {b' : WalkM Unit}
    (hc : Lifts c c') (hb : Lifts b b') :
    Lifts (PlanarWalkM.loop fuel c b) (WalkM.loop fuel c' b') := by
  induction fuel with
  | zero => exact Lifts.pure _
  | succ n ih =>
    unfold PlanarWalkM.loop WalkM.loop
    exact Lifts.bind hc fun a => Lifts.ite _ (fun _ => Lifts.bind hb fun _ => ih) (fun _ => Lifts.pure _)

theorem Lifts.loopFirst (fuel : Nat) {c : PlanarWalkM Bool} {c' : WalkM Bool} {b : Bool → PlanarWalkM Unit}
    {b' : WalkM Unit} (hc : Lifts c c') (hb : ∀ x, Lifts (b x) b') (first : Bool) :
    Lifts (PlanarWalkM.loopFirst fuel c b first) (WalkM.loop fuel c' b') := by
  induction fuel generalizing first with
  | zero => exact Lifts.pure _
  | succ n ih =>
    unfold PlanarWalkM.loopFirst WalkM.loop
    exact Lifts.bind hc fun a => Lifts.ite _ (fun _ => Lifts.bind (hb _) fun _ => ih _) (fun _ => Lifts.pure _)

namespace PlanarWalkM

theorem lifts_allocItem (type) : Lifts (allocItem type) (WalkM.allocItem type) := fun _ => rfl
theorem lifts_pushTstack (v d item pl) : Lifts (pushTstack v d item pl) (WalkM.pushTstack v d item) := fun _ => rfl
theorem lifts_pushVertTstack (v d) : Lifts (pushVertTstack v d) (WalkM.pushVertTstack v d) := fun _ => rfl
theorem lifts_pushEdgeTstack (v d e b) : Lifts (pushEdgeTstack v d e b) (WalkM.pushEdgeTstack v d e) := fun _ => rfl
theorem lifts_mergeTstackTops : Lifts mergeTstackTops WalkM.mergeTstackTops := fun _ => rfl

end PlanarWalkM

syntax "lifts0" : tactic
macro_rules
  | `(tactic| lifts0) => `(tactic| first
    | with_reducible exact Lifts.pure _
    | with_reducible exact Lifts.liftW _
    | with_reducible first
      | exact PlanarWalkM.lifts_allocItem _ | exact PlanarWalkM.lifts_pushVertTstack _ _
      | exact PlanarWalkM.lifts_pushEdgeTstack _ _ _ _ | exact PlanarWalkM.lifts_mergeTstackTops
    | (with_reducible refine Lifts.auxBind ?_ fun _ => ?_
       · aux_only
       · lifts0)
    | (with_reducible refine Lifts.auxThen ?_ ?_
       · aux_only
       · lifts0)
    | (with_reducible refine Lifts.bind ?_ fun _ => ?_ <;> lifts0)
    | (with_reducible refine Lifts.ite _ (fun _ => ?_) (fun _ => ?_) <;> lifts0)
    | (with_reducible refine Lifts.fmap ?_ _; lifts0)
    | (with_reducible refine Lifts.loop _ ?_ ?_ <;> lifts0)
    | (with_reducible refine Lifts.loopFirst _ ?_ (fun _ => ?_) _ <;> lifts0)
    | (with_reducible refine Lifts.bindAux ?_ fun _ => ?_
       · lifts0
       · aux_only)
    | (split <;> lifts0)
    | (dsimp only; lifts0)
    | skip)

namespace PlanarWalkM

theorem lifts_maybeUnwrapNxt (type b) : Lifts (maybeUnwrapNxt type b) (WalkM.maybeUnwrapNxt type) := by
  unfold maybeUnwrapNxt WalkM.maybeUnwrapNxt
  lifts0

theorem lifts_finishTstackTop (item b) : Lifts (finishTstackTop item b) (WalkM.finishTstackTop item) := by
  unfold finishTstackTop; lifts0

end PlanarWalkM

set_option linter.deprecated false in
theorem Lifts.haveVal {β : Sort u} (v : β) {rest : β → PlanarWalkM α} {rest' : β → WalkM α}
    (h : Lifts (rest v) (rest' v)) : Lifts (letFun v rest) (letFun v rest') := h

set_option linter.deprecated false in
theorem Lifts.haveJp1 {β : Sort u} {f : β → PlanarWalkM γ} {f' : β → WalkM γ}
    {rest : (β → PlanarWalkM γ) → PlanarWalkM α} {rest' : (β → WalkM γ) → WalkM α}
    (hf : ∀ x, Lifts (f x) (f' x))
    (h : ∀ jp jp', (∀ x, Lifts (jp x) (jp' x)) → Lifts (rest jp) (rest' jp')) :
    Lifts (letFun f rest) (letFun f' rest') := h f f' hf

set_option linter.deprecated false in
theorem Lifts.haveJp2 {β₁ : Sort u} {β₂ : Sort v} {f : β₁ → β₂ → PlanarWalkM γ} {f' : β₁ → β₂ → WalkM γ}
    {rest : (β₁ → β₂ → PlanarWalkM γ) → PlanarWalkM α} {rest' : (β₁ → β₂ → WalkM γ) → WalkM α}
    (hf : ∀ x y, Lifts (f x y) (f' x y))
    (h : ∀ jp jp', (∀ x y, Lifts (jp x y) (jp' x y)) → Lifts (rest jp) (rest' jp')) :
    Lifts (letFun f rest) (letFun f' rest') := h f f' hf

set_option linter.deprecated false in
theorem Lifts.haveJpOnly1 {β : Sort u} {f : β → PlanarWalkM α} {rest : (β → PlanarWalkM α) → PlanarWalkM α}
    {m' : WalkM α} (hf : ∀ x, Lifts (f x) m')
    (h : ∀ jp, (∀ x, Lifts (jp x) m') → Lifts (rest jp) m') : Lifts (letFun f rest) m' := h f hf

set_option linter.deprecated false in
theorem Lifts.haveJpOnly2 {β₁ : Sort u} {β₂ : Sort v} {f : β₁ → β₂ → PlanarWalkM α}
    {rest : (β₁ → β₂ → PlanarWalkM α) → PlanarWalkM α} {m' : WalkM α} (hf : ∀ x y, Lifts (f x y) m')
    (h : ∀ jp, (∀ x y, Lifts (jp x y) m') → Lifts (rest jp) m') : Lifts (letFun f rest) m' := h f hf

theorem Lifts.iteAux (c : Prop) [Decidable c] {t e : PlanarWalkM α} {m' : WalkM α}
    (ht : c → Lifts t m') (he : ¬c → Lifts e m') : Lifts (if c then t else e) m' := by
  split
  · exact ht ‹_›
  · exact he ‹_›

open Lean Elab Tactic Meta in
/-- `lift_have R jp1 jp2 only1 only2` peels one `have` off a goal `R m m'`. A `have` present on
both sides: a pure value is zeta-substituted, a join point (function-typed `have`) is abstracted via
`jp1` / `jp2` (the `Lifts.haveJp1` / `Lifts.haveJp2` shapes). A `have` present only in `m`: a value
is zeta-substituted, a join point is abstracted via `only1` / `only2` (`Lifts.haveJpOnly1/2`). -/
elab "lift_have" rel:ident jp1:ident jp2:ident only1:ident only2:ident : tactic => do
  let relName ← realizeGlobalConstNoOverloadWithInfo rel
  let jp1Name ← realizeGlobalConstNoOverloadWithInfo jp1
  let jp2Name ← realizeGlobalConstNoOverloadWithInfo jp2
  let only1Name ← realizeGlobalConstNoOverloadWithInfo only1
  let only2Name ← realizeGlobalConstNoOverloadWithInfo only2
  let g ← getMainGoal
  let ty ← instantiateMVars (← g.getType)
  unless ty.isAppOfArity relName 3 do throwError "lift_have: not a {relName} goal"
  let args := ty.getAppArgs
  let α := args[0]!
  let .letE n t v b _ := args[1]! | throwError "lift_have: no have"
  let mkRel (m m' : Expr) : Expr := mkApp3 ty.getAppFn α m m'
  let applyJp (g : MVarId) (lem : Name) : TacticM Unit := do
    let gs ← g.apply (← mkConstWithFreshMVarLevels lem)
    let gs ← gs.filterM fun g => return !(← g.isAssigned)
    let mut out := #[]
    for g in gs do
      let (_, g) ← g.intros
      out := out.push g
    replaceMainGoal out.toList
  match args[2]! with
  | .letE n' t' v' b' _ =>
    if t.isForall then
      let lhs ← mkAppM ``letFun #[v, Expr.lam n t b .default]
      let rhs ← mkAppM ``letFun #[v', Expr.lam n' t' b' .default]
      let g' ← g.replaceTargetDefEq (mkRel lhs rhs)
      applyJp g' (if t.getForallArity == 2 then jp2Name else jp1Name)
    else
      unless ← isDefEq v v' do throwError "lift_have: have values differ"
      replaceMainGoal [← g.replaceTargetDefEq (mkRel (b.instantiate1 v) (b'.instantiate1 v'))]
  | m' =>
    if t.isForall then
      let lhs ← mkAppM ``letFun #[v, Expr.lam n t b .default]
      let g' ← g.replaceTargetDefEq (mkRel lhs m')
      applyJp g' (if t.getForallArity == 2 then only2Name else only1Name)
    else
      replaceMainGoal [← g.replaceTargetDefEq (mkRel (b.instantiate1 v) m')]

syntax "lifts" : tactic
macro_rules
  | `(tactic| lifts) => `(tactic| first
    | with_reducible exact Lifts.pure _
    | with_reducible exact Lifts.liftW _
    | with_reducible first
      | exact PlanarWalkM.lifts_allocItem _ | exact PlanarWalkM.lifts_pushVertTstack _ _
      | exact PlanarWalkM.lifts_pushEdgeTstack _ _ _ _ | exact PlanarWalkM.lifts_mergeTstackTops
      | exact PlanarWalkM.lifts_maybeUnwrapNxt _ _ | exact PlanarWalkM.lifts_finishTstackTop _ _
    | (with_reducible refine Lifts.bind (Lifts.liftW _) fun _ => ?_; lifts)
    | (with_reducible refine Lifts.auxBind ?_ fun _ => ?_
       · aux_only
       · lifts)
    | (with_reducible refine Lifts.auxThen ?_ ?_
       · aux_only
       · lifts)
    | with_reducible solve_by_elim only [*]
    | (lift_have Lifts Lifts.haveJp1 Lifts.haveJp2 Lifts.haveJpOnly1 Lifts.haveJpOnly2 <;> lifts)
    | (with_reducible refine Lifts.bind ?_ fun _ => ?_ <;> lifts)
    | (with_reducible refine Lifts.ite _ (fun _ => ?_) (fun _ => ?_) <;> lifts)
    | (with_reducible refine Lifts.iteAux _ (fun _ => ?_) (fun _ => ?_) <;> lifts)
    | (with_reducible refine Lifts.fmap ?_ _; lifts)
    | (with_reducible refine Lifts.loop _ ?_ ?_ <;> lifts)
    | (with_reducible refine Lifts.loopFirst _ ?_ (fun _ => ?_) _ <;> lifts)
    | (with_reducible refine Lifts.bindAux ?_ fun _ => ?_
       · lifts
       · aux_only)
    | (rename_i h1 _ h2; cases h1.symm.trans h2)
    | (split <;> lifts)
    | (dsimp only; lifts)
    | skip)

end Spqr

namespace Spqr

namespace PlanarWalkM

theorem lifts_planarFinishEdge (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool) :
    Lifts (planarFinishEdge curV d o origTstack hasVert) (finishEdge curV d o origTstack hasVert) := by
  unfold planarFinishEdge finishEdge; lifts

end PlanarWalkM

mutual
theorem lifts_planarWalkTree (t : DfsTree) (d : Nat) : Lifts (planarWalkTree t d) (walkTree t d) := by
  unfold planarWalkTree walkTree
  match t with
  | .node v outs =>
    with_reducible refine Lifts.bind (Lifts.liftW _) fun _ => ?_
    with_reducible refine Lifts.bind (lifts_planarWalkOuts v d outs false) fun _ => ?_
    lifts

theorem lifts_planarWalkOuts (v d : Nat) (outs : List DfsOut) (hasVert : Bool) :
    Lifts (planarWalkOuts v d outs hasVert) (walkOuts v d outs hasVert) := by
  unfold planarWalkOuts walkOuts
  match outs with
  | [] => lifts
  | o :: rest =>
    with_reducible refine Lifts.bind (lifts_planarWalkOut v d o hasVert) fun _ => ?_
    exact lifts_planarWalkOuts v d rest _

theorem lifts_planarWalkOut (v d : Nat) (o : DfsOut) (hasVert : Bool) :
    Lifts (planarWalkOut v d o hasVert) (walkOut v d o hasVert) := by
  unfold planarWalkOut walkOut
  match o with
  | .tree a b child =>
    lifts
    all_goals first
      | exact lifts_planarWalkTree child (d + 1)
      | exact PlanarWalkM.lifts_planarFinishEdge _ _ _ _ _
  | .back a b c =>
    lifts
    all_goals exact PlanarWalkM.lifts_planarFinishEdge _ _ _ _ _
end

theorem lifts_planarWalkForest (forest : List DfsTree) : Lifts (planarWalkForest forest) (walkForest forest) := by
  unfold planarWalkForest walkForest
  refine Lifts.forM _ fun t => ?_
  with_reducible refine Lifts.bind (lifts_planarWalkTree t 0) fun _ => ?_
  lifts

/-- The planar walk is the ordinary walk plus auxiliary state. -/
theorem planarWalk_base (g : Graph) (ternarize : Bool) (forest : List DfsTree) :
    (g.planarWalk ternarize forest).base = g.walk ternarize forest := by
  have h := lifts_planarWalkForest forest (PlanarWalkState.init g ternarize)
  unfold Graph.planarWalk Graph.walk
  exact (congrArg Prod.snd h).symm

theorem planarWalk_proj (g : Graph) (ternarize : Bool) (forest : List DfsTree) :
    (g.planarWalk ternarize forest).base.items = (g.walk ternarize forest).items := by
  rw [planarWalk_base]

end Spqr
