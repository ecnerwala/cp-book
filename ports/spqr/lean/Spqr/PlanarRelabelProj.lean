import Spqr.Build
import Spqr.PlanarRelabel
import Spqr.PlanarWalkProj

/-!
# Projection of the planar relabel onto the ordinary relabel

`LiftsR m m'` says the `PlanarRelabelM` computation `m` acts on the `RelabelState` component
exactly as `m'`; `AuxOnlyR m` says `m` only touches the planarity component. `planarRelabel` is
`relabel` interleaved with `AuxOnlyR` steps, which gives `planarRelabelTree_base` and the
projection theorem `planarRelabel_proj`.
-/

namespace Spqr

def LiftsR (m : PlanarRelabelM α) (m' : RelabelM α) : Prop :=
  ∀ s, m' s.base = ((m s).1, (m s).2.base)

def AuxOnlyR (m : PlanarRelabelM α) : Prop := ∀ s, (m s).2.base = s.base

namespace LiftsR

theorem liftR (m : RelabelM α) : LiftsR (PlanarRelabelM.liftR m) m := fun _ => rfl

theorem pure (a : α) : LiftsR (Pure.pure a) (Pure.pure a) := fun _ => rfl

theorem bind {m : PlanarRelabelM α} {m' : RelabelM α} {k : α → PlanarRelabelM γ} {k' : α → RelabelM γ}
    (h : LiftsR m m') (hk : ∀ a, LiftsR (k a) (k' a)) : LiftsR (m >>= k) (m' >>= k') := by
  intro s
  show k' (m' s.base).1 (m' s.base).2 = _
  rw [h s]
  exact hk _ _

theorem auxBind {m : PlanarRelabelM α} {k : α → PlanarRelabelM γ} {m' : RelabelM γ}
    (h : AuxOnlyR m) (hk : ∀ a, LiftsR (k a) m') : LiftsR (m >>= k) m' := by
  intro s
  show m' s.base = ((k (m s).1 (m s).2).1, (k (m s).1 (m s).2).2.base)
  rw [← h s]
  exact hk _ _

theorem auxThen {m : PlanarRelabelM Unit} {k : PlanarRelabelM γ} {m' : RelabelM γ}
    (h : AuxOnlyR m) (hk : LiftsR k m') : LiftsR (m *> k) m' := by
  intro s
  show m' s.base = ((k (m s).2).1, (k (m s).2).2.base)
  rw [← h s]
  exact hk _

theorem bindAux {m : PlanarRelabelM Unit} {m' : RelabelM Unit} (h : LiftsR m m')
    {k : Unit → PlanarRelabelM Unit} (hk : ∀ a, AuxOnlyR (k a)) : LiftsR (m >>= k) m' := by
  intro s
  show m' s.base = ((k (m s).1 (m s).2).1, (k (m s).1 (m s).2).2.base)
  rw [hk _ _, h s]

theorem ite (c : Prop) [Decidable c] {t e : PlanarRelabelM α} {t' e' : RelabelM α}
    (ht : c → LiftsR t t') (he : ¬c → LiftsR e e') :
    LiftsR (if c then t else e) (if c then t' else e') := by
  split
  · exact ht ‹_›
  · exact he ‹_›

theorem forIn (l : List γ) (init : σ) {f : γ → σ → PlanarRelabelM (ForInStep σ)}
    {f' : γ → σ → RelabelM (ForInStep σ)} (h : ∀ x b, LiftsR (f x b) (f' x b)) :
    LiftsR (forIn l init f) (forIn l init f') := by
  induction l generalizing init with
  | nil => exact LiftsR.pure _
  | cons x rest ih =>
    simp only [List.forIn_cons]
    refine LiftsR.bind (h x init) fun r => ?_
    cases r
    · exact LiftsR.pure _
    · exact ih _

set_option linter.deprecated false in
theorem haveJp1 {β : Sort u} {f : β → PlanarRelabelM γ} {f' : β → RelabelM γ}
    {rest : (β → PlanarRelabelM γ) → PlanarRelabelM α} {rest' : (β → RelabelM γ) → RelabelM α}
    (hf : ∀ x, LiftsR (f x) (f' x))
    (h : ∀ jp jp', (∀ x, LiftsR (jp x) (jp' x)) → LiftsR (rest jp) (rest' jp')) :
    LiftsR (letFun f rest) (letFun f' rest') := h f f' hf

set_option linter.deprecated false in
theorem haveJp2 {β₁ : Sort u} {β₂ : Sort v} {f : β₁ → β₂ → PlanarRelabelM γ} {f' : β₁ → β₂ → RelabelM γ}
    {rest : (β₁ → β₂ → PlanarRelabelM γ) → PlanarRelabelM α} {rest' : (β₁ → β₂ → RelabelM γ) → RelabelM α}
    (hf : ∀ x y, LiftsR (f x y) (f' x y))
    (h : ∀ jp jp', (∀ x y, LiftsR (jp x y) (jp' x y)) → LiftsR (rest jp) (rest' jp')) :
    LiftsR (letFun f rest) (letFun f' rest') := h f f' hf

set_option linter.deprecated false in
theorem haveJpOnly1 {β : Sort u} {f : β → PlanarRelabelM α} {rest : (β → PlanarRelabelM α) → PlanarRelabelM α}
    {m' : RelabelM α} (hf : ∀ x, LiftsR (f x) m')
    (h : ∀ jp, (∀ x, LiftsR (jp x) m') → LiftsR (rest jp) m') : LiftsR (letFun f rest) m' := h f hf

set_option linter.deprecated false in
theorem haveJpOnly2 {β₁ : Sort u} {β₂ : Sort v} {f : β₁ → β₂ → PlanarRelabelM α}
    {rest : (β₁ → β₂ → PlanarRelabelM α) → PlanarRelabelM α} {m' : RelabelM α} (hf : ∀ x y, LiftsR (f x y) m')
    (h : ∀ jp, (∀ x y, LiftsR (jp x y) m') → LiftsR (rest jp) m') : LiftsR (letFun f rest) m' := h f hf

theorem iteAux (c : Prop) [Decidable c] {t e : PlanarRelabelM α} {m' : RelabelM α}
    (ht : c → LiftsR t m') (he : ¬c → LiftsR e m') : LiftsR (if c then t else e) m' := by
  split
  · exact ht ‹_›
  · exact he ‹_›

end LiftsR

namespace AuxOnlyR

theorem pure (a : α) : AuxOnlyR (Pure.pure a : PlanarRelabelM α) := fun _ => rfl

theorem bind {m : PlanarRelabelM α} {k : α → PlanarRelabelM β}
    (h : AuxOnlyR m) (hk : ∀ a, AuxOnlyR (k a)) : AuxOnlyR (m >>= k) := by
  intro s
  show (k (m s).1 (m s).2).2.base = s.base
  rw [hk _ _, h s]

theorem ite (c : Prop) [Decidable c] {t e : PlanarRelabelM α}
    (ht : c → AuxOnlyR t) (he : ¬c → AuxOnlyR e) : AuxOnlyR (if c then t else e) := by
  split
  · exact ht ‹_›
  · exact he ‹_›

theorem getAux : AuxOnlyR PlanarRelabelM.getAux := fun _ => rfl
theorem modifyAux (f) : AuxOnlyR (PlanarRelabelM.modifyAux f) := fun _ => rfl
theorem liftAux (m : StateM PlanarRelabelAux α) : AuxOnlyR (PlanarRelabelM.liftAux m) := fun _ => rfl
theorem setupNode (g type cur) : AuxOnlyR (PlanarRelabelM.setupNode g type cur) := fun _ => rfl
theorem applyFlips (g it flips) : AuxOnlyR (PlanarRelabelM.applyFlips g it flips) := fun _ => rfl

end AuxOnlyR

syntax "aux_only_r" : tactic
macro_rules
  | `(tactic| aux_only_r) => `(tactic| first
    | with_reducible exact AuxOnlyR.pure _
    | with_reducible first
      | exact AuxOnlyR.getAux | exact AuxOnlyR.modifyAux _ | exact AuxOnlyR.liftAux _
      | exact AuxOnlyR.setupNode _ _ _ | exact AuxOnlyR.applyFlips _ _ _
    | (with_reducible refine AuxOnlyR.bind ?_ fun _ => ?_ <;> aux_only_r)
    | (with_reducible refine AuxOnlyR.ite _ (fun _ => ?_) (fun _ => ?_) <;> aux_only_r)
    | (split <;> aux_only_r)
    | (dsimp only; aux_only_r))

syntax "lifts_r" : tactic
macro_rules
  | `(tactic| lifts_r) => `(tactic| first
    | with_reducible exact LiftsR.pure _
    | with_reducible exact LiftsR.liftR _
    | with_reducible solve_by_elim only [*]
    | (with_reducible refine LiftsR.bind (LiftsR.liftR _) fun _ => ?_; lifts_r)
    | (with_reducible refine LiftsR.auxBind ?_ fun _ => ?_
       · aux_only_r
       · lifts_r)
    | (with_reducible refine LiftsR.auxThen ?_ ?_
       · aux_only_r
       · lifts_r)
    | (lift_have LiftsR LiftsR.haveJp1 LiftsR.haveJp2 LiftsR.haveJpOnly1 LiftsR.haveJpOnly2 <;> lifts_r)
    | (with_reducible refine LiftsR.bind ?_ fun _ => ?_ <;> lifts_r)
    | (with_reducible refine LiftsR.ite _ (fun _ => ?_) (fun _ => ?_) <;> lifts_r)
    | (with_reducible refine LiftsR.iteAux _ (fun _ => ?_) (fun _ => ?_) <;> lifts_r)
    | (with_reducible refine LiftsR.forIn _ _ fun _ _ => ?_; lifts_r)
    | (with_reducible refine LiftsR.bindAux ?_ fun _ => ?_
       · lifts_r
       · aux_only_r)
    | (rename_i h1 _ h2; cases h1.symm.trans h2)
    | (split <;> lifts_r)
    | (dsimp only; lifts_r)
    | skip)

set_option maxRecDepth 8000 in
theorem liftsR_planarRelabel (fuel : Nat) (cur : ItemId) (parent parNv capTwin : Option Nat) :
    LiftsR (planarRelabel fuel cur parent parNv capTwin) (relabel fuel cur parent parNv capTwin) := by
  induction fuel generalizing cur parent parNv capTwin with
  | zero => exact LiftsR.pure _
  | succ fuel ih =>
    unfold planarRelabel relabel
    lifts_r

theorem planarRelabelTree_base (g : Graph) (w : PlanarWalkState) :
    (planarRelabelTree g w).toSpqrTree = relabelTree g w.base.items := by
  rw [relabelTree_eq_ofState]
  have h := liftsR_planarRelabel w.base.items.size rootItem none none none (PlanarRelabelState.init g w)
  exact (congrArg (fun r => SpqrTree.ofState g r.2) h).symm

/-- The planar SPQR tree is the ordinary SPQR tree plus planarity data. -/
theorem planarRelabel_proj (g : Graph) (ternarize : Bool) (vertOrder edgeOrder : List Nat) :
    (g.planarSpqrTree ternarize vertOrder edgeOrder).toSpqrTree = g.spqrTree ternarize vertOrder edgeOrder := by
  unfold Graph.planarSpqrTree Graph.spqrTree
  rw [planarRelabelTree_base, planarWalk_base]

end Spqr
