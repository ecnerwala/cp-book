import Spqr.Frame
import Spqr.WalkTyping
import Spqr.StEar

/-!
# `stackDir` frame facts

`walkTree t d` only writes `stackDir[d']` for `d' ≥ d`: `walkOut` writes slot `d`, the recursive
`walkTree child (d + 1)` writes slots `≥ d + 1`, and inside `finishEdge` the only write is
`setStackDir (← nxt).topDepth` under the guard `(← nxt).topDepth > d` (`loop1Type`). Everything
else leaves `stackDir` untouched.
-/

namespace Spqr
open WalkM

/-- `m` never changes `stackDir`. -/
def KeepsAll (m : WalkM α) : Prop := ∀ s, wp m (fun _ s' => s'.stackDir = s.stackDir) s

/-- The `wp` simp set computing straight-line `WalkM` blocks. -/
syntax "wp_simp" : tactic
macro_rules
  | `(tactic| wp_simp) => `(tactic| simp only [wp_bind, wp_ite, wp_pure, wp_map, wp_get, wp_set,
      wp_modify, wp_modifyItem, wp_getItem, wp_allocItem, wp_stackDir, wp_setStackDir, wp_makeVs,
      wp_cur, wp_nxt, wp_tstackSize, wp_modifyCur, wp_modifyNxt, wp_popTstack, wp_pushTstack,
      wp_pushVertTstack, wp_pushEdgeTstack, wp_mergeTstackTops])

/-- `m` leaves `stackDir[d']` unchanged for every `d' < d`. -/
def KeepsBelow (d : Nat) (m : WalkM α) : Prop :=
  ∀ s d', d' < d → wp m (fun _ s' => s'.stackDir[d']! = s.stackDir[d']!) s

namespace KeepsAll

theorem keepsBelow {m : WalkM α} (h : KeepsAll m) (d : Nat) : KeepsBelow d m :=
  fun s d' _ => by
    have e : ((m.run s).2).stackDir = s.stackDir := h s
    exact congrArg (fun a : Array Bool => a[d']!) e

theorem pure (a : α) : KeepsAll (pure a : WalkM α) := fun _ => rfl

theorem bind {m : WalkM β} {f : β → WalkM α} (hm : KeepsAll m) (hf : ∀ b, KeepsAll (f b)) :
    KeepsAll (m >>= f) := by
  intro s
  show (((f (m.run s).1).run (m.run s).2)).2.stackDir = s.stackDir
  rw [hf _ _, hm s]

theorem map {m : WalkM β} (f : β → α) (hm : KeepsAll m) : KeepsAll (f <$> m) := hm

theorem ite (c : Prop) [Decidable c] {a b : WalkM α} (ha : KeepsAll a) (hb : KeepsAll b) :
    KeepsAll (if c then a else b) := by
  split <;> assumption

theorem modify (f : WalkState → WalkState) (hf : ∀ s, (f s).stackDir = s.stackDir) :
    KeepsAll (modify f : WalkM Unit) := hf

theorem get : KeepsAll (get : WalkM WalkState) := fun _ => rfl
theorem stackDir (d : Nat) : KeepsAll (stackDir d) := fun _ => rfl
theorem makeVs (a b : Nat) : KeepsAll (makeVs a b) := fun _ => rfl
theorem cur : KeepsAll cur := fun _ => rfl
theorem nxt : KeepsAll nxt := fun _ => rfl
theorem tstackSize : KeepsAll tstackSize := fun _ => rfl
theorem getItem (i : ItemId) : KeepsAll (getItem i) := fun _ => rfl
theorem modifyItem (i : ItemId) (f : Item → Item) : KeepsAll (modifyItem i f) := fun _ => rfl
theorem allocItem (ty : NodeType) : KeepsAll (allocItem ty) := fun _ => rfl
theorem modifyCur (f : TEntry → TEntry) : KeepsAll (modifyCur f) := fun _ => rfl
theorem modifyNxt (f : TEntry → TEntry) : KeepsAll (modifyNxt f) := fun _ => rfl
theorem popTstack : KeepsAll popTstack := fun _ => rfl
theorem pushTstack (v d : Nat) (i : ItemId) : KeepsAll (pushTstack v d i) := fun _ => rfl
theorem pushVertTstack (v d : Nat) : KeepsAll (pushVertTstack v d) := fun _ => rfl
theorem pushEdgeTstack (v d e : Nat) : KeepsAll (pushEdgeTstack v d e) := fun _ => rfl
theorem mergeTstackTops : KeepsAll mergeTstackTops := fun _ => rfl

theorem loop (n : Nat) {cond : WalkM Bool} {body : WalkM Unit} (hc : KeepsAll cond)
    (hb : KeepsAll body) : KeepsAll (loop n cond body) := by
  induction n with
  | zero => exact pure ()
  | succ n ih =>
    simp only [WalkM.loop]
    exact bind hc fun b => ite _ (bind hb fun _ => ih) (pure ())

theorem maybeUnwrapNxt (ty : NodeType) : KeepsAll (maybeUnwrapNxt ty) := by
  intro s; unfold WalkM.maybeUnwrapNxt; wp_simp
  repeat' split
  all_goals trivial

theorem finishTstackTop (i : ItemId) : KeepsAll (finishTstackTop i) := fun _ => rfl

theorem loop1Cond (d : Nat) : KeepsAll (loop1Cond d) := by
  unfold Spqr.loop1Cond
  exact bind tstackSize fun _ => bind nxt fun _ => pure _

theorem loop2Cond (fo : Nat) : KeepsAll (loop2Cond fo) := by
  unfold Spqr.loop2Cond
  exact bind cur fun _ => pure _

theorem loop3Cond (o : Nat) : KeepsAll (loop3Cond o) := by
  unfold Spqr.loop3Cond
  exact bind tstackSize fun _ => pure _

theorem mergeLate (d : Nat) : KeepsAll (mergeLate d) := by
  unfold Spqr.mergeLate
  refine bind get fun s => bind cur fun t => ite _ ?_ (pure _)
  exact bind tstackSize fun n => bind (loop n (loop2Cond _) mergeTstackTops) fun _ => pure _

theorem condP (curV lowval : Nat) (isType1 : Bool) : KeepsAll (condP curV lowval isType1) := by
  unfold Spqr.condP
  exact bind tstackSize fun _ => bind nxt fun _ => bind nxt fun _ => pure _

theorem finishP (curV lowval : Nat) (isType1 : Bool) : KeepsAll (finishP curV lowval isType1) := by
  intro s; unfold Spqr.finishP Spqr.condP WalkM.maybeUnwrapNxt WalkM.finishTstackTop; wp_simp
  repeat' split
  all_goals trivial

theorem finishTail (curV d : Nat) (hasVert isSingle : Bool) :
    KeepsAll (finishTail curV d hasVert isSingle) := by
  intro s; unfold Spqr.finishTail; wp_simp
  repeat' split
  all_goals trivial

theorem finishRest (curV d lowval : Nat) (isType1 hasVert isSingle : Bool) :
    KeepsAll (finishRest curV d lowval isType1 hasVert isSingle) := by
  unfold Spqr.finishRest
  exact bind (finishP _ _ _) fun _ => finishTail _ _ _ _

theorem closeVertTail (curV : Nat) (edgeDir isSingle : Bool) {k : Bool → WalkM Bool}
    (hk : ∀ b, KeepsAll (k b)) (item : Option ItemId) :
    KeepsAll (closeVertTail curV edgeDir isSingle k item) := by
  unfold Spqr.closeVertTail
  refine bind mergeTstackTops fun _ => bind mergeTstackTops fun _ => bind (modifyCur _) fun _ => ?_
  cases item with
  | some i => exact bind (finishTstackTop i) fun _ => hk _
  | none => exact hk _

theorem closeVert (curV : Nat) (edgeDir isType1 : Bool) (origTstack : Nat) (isSingle : Bool)
    {k : Bool → WalkM Bool} (hk : ∀ b, KeepsAll (k b)) :
    KeepsAll (closeVert curV edgeDir isType1 origTstack isSingle k) := by
  unfold Spqr.closeVert
  refine ite _ ?_ ?_
  · exact bind tstackSize fun n => bind (loop n (loop3Cond _) mergeTstackTops) fun _ =>
      closeVertTail _ _ _ hk _
  · exact bind (map _ (maybeUnwrapNxt _)) fun i => closeVertTail _ _ _ hk _

theorem finishBoundary (curV d : Nat) (o : DfsOut) (qItem : ItemId) (hasVert : Bool) :
    KeepsAll (finishBoundary curV d o qItem hasVert) := by
  intro s; unfold Spqr.finishBoundary; wp_simp
  repeat' split
  all_goals trivial

theorem finishBack (curV d : Nat) (o : DfsOut) (hasVert : Bool) :
    KeepsAll (finishBack curV d o hasVert) := by
  unfold Spqr.finishBack
  exact bind (pushEdgeTstack _ _ _) fun _ => bind (modify _ fun _ => rfl) fun _ => finishRest _ _ _ _ _ _

end KeepsAll

namespace KeepsBelow

theorem mono {m : WalkM α} {d d' : Nat} (h : KeepsBelow d m) (hd : d' ≤ d) : KeepsBelow d' m :=
  fun s i hi => h s i (Nat.lt_of_lt_of_le hi hd)

theorem wp_of {d d' : Nat} {m : WalkM α} (h : KeepsBelow d m) (hd' : d' < d) {s₀ s : WalkState}
    (hs : s.stackDir[d']! = s₀.stackDir[d']!) :
    wp m (fun _ s' => s'.stackDir[d']! = s₀.stackDir[d']!) s :=
  wp_mono _ (h s d' hd') fun _ _ e => e.trans hs

theorem bind {d : Nat} {m : WalkM β} {f : β → WalkM α} (hm : KeepsBelow d m)
    (hf : ∀ b, KeepsBelow d (f b)) : KeepsBelow d (m >>= f) := by
  intro s i hi
  show (((f (m.run s).1).run (m.run s).2)).2.stackDir[i]! = s.stackDir[i]!
  rw [hf _ _ i hi, hm s i hi]

theorem ite {d : Nat} (c : Prop) [Decidable c] {a b : WalkM α} (ha : c → KeepsBelow d a)
    (hb : ¬ c → KeepsBelow d b) : KeepsBelow d (if c then a else b) := by
  split
  · exact ha ‹_›
  · exact hb ‹_›

theorem setStackDir {d i : Nat} (hi : d ≤ i) (b : Bool) : KeepsBelow d (setStackDir i b) := by
  intro s j hj
  exact Array.getElem!_set!_ne _ _ _ _ (by omega)

theorem loop {d : Nat} (n : Nat) {cond : WalkM Bool} {body : WalkM Unit} (hc : KeepsBelow d cond)
    (hb : KeepsBelow d body) : KeepsBelow d (loop n cond body) := by
  induction n with
  | zero => exact (KeepsAll.pure ()).keepsBelow d
  | succ n ih =>
    simp only [WalkM.loop]
    exact bind hc fun b => ite _ (fun _ => bind hb fun _ => ih) fun _ => (KeepsAll.pure ()).keepsBelow d

theorem loop1Type (d : Nat) (edgeDir : Bool) : KeepsBelow d (loop1Type d edgeDir) := by
  intro s d' hd'
  unfold Spqr.loop1Type; wp_simp
  repeat' split
  all_goals first | trivial | rfl | exact Array.getElem!_set!_ne _ _ _ _ (by omega)

theorem loop1Body (d : Nat) (edgeDir : Bool) : KeepsBelow d (loop1Body d edgeDir) := by
  unfold Spqr.loop1Body
  refine bind (loop1Type d edgeDir) fun ty => bind ((KeepsAll.maybeUnwrapNxt ty).keepsBelow d) fun i => ?_
  exact bind (KeepsAll.mergeTstackTops.keepsBelow d) fun _ => (KeepsAll.finishTstackTop i).keepsBelow d

theorem closeEars (nxtV d e : Nat) (edgeDir : Bool) : KeepsBelow d (closeEars nxtV d e edgeDir) := by
  unfold Spqr.closeEars
  refine bind ((KeepsAll.pushEdgeTstack _ _ _).keepsBelow d) fun _ =>
    bind (KeepsAll.tstackSize.keepsBelow d) fun n => ?_
  exact loop n ((KeepsAll.loop1Cond d).keepsBelow d) (loop1Body d edgeDir)

theorem finishTree (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert edgeDir : Bool) :
    KeepsBelow d (finishTree curV d o origTstack hasVert edgeDir) := by
  unfold Spqr.finishTree
  refine bind (closeEars _ _ _ _) fun _ => bind ((KeepsAll.mergeLate d).keepsBelow d) fun isSingle => ?_
  refine ite _ (fun _ => ?_) fun _ => (KeepsAll.finishRest _ _ _ _ _ _).keepsBelow d
  exact (KeepsAll.closeVert _ _ _ _ _ fun b => KeepsAll.finishRest _ _ _ _ _ _).keepsBelow d

theorem finishEdge (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool) :
    KeepsBelow d (finishEdge curV d o origTstack hasVert) := by
  rw [finishEdge_eq]
  unfold finishEdge'
  refine bind (KeepsAll.get.keepsBelow d) fun s => bind ((KeepsAll.stackDir d).keepsBelow d) fun edgeDir => ?_
  refine ite _ (fun _ => (KeepsAll.finishBoundary _ _ _ _ _).keepsBelow d) fun _ => ?_
  refine bind ((KeepsAll.makeVs _ _).keepsBelow d) fun vs => bind ((KeepsAll.modifyItem _ _).keepsBelow d) fun _ => ?_
  exact ite _ (fun _ => finishTree _ _ _ _ _ _) fun _ => (KeepsAll.finishBack _ _ _ _).keepsBelow d

end KeepsBelow

theorem walk_keepsBelow :
    (∀ (t : DfsTree) (d : Nat), KeepsBelow d (walkTree t d)) ∧
    (∀ (v d : Nat) (outs : List DfsOut) (hasVert : Bool), KeepsBelow d (walkOuts v d outs hasVert)) ∧
    (∀ (v d : Nat) (o : DfsOut) (hasVert : Bool), KeepsBelow d (walkOut v d o hasVert)) := by
  refine walkTree.mutual_induct _ _ _ ?_ ?_ ?_ ?_
  · intro d v outs ih
    simp only [walkTree]
    refine KeepsBelow.bind ((KeepsAll.modify _ fun _ => rfl).keepsBelow d) fun _ =>
      KeepsBelow.bind ih fun hasVert => ?_
    split
    · exact (KeepsAll.pure _).keepsBelow d
    · exact KeepsBelow.bind (KeepsBelow.setStackDir (Nat.le_refl d) _) fun _ =>
        (KeepsAll.pushVertTstack _ _).keepsBelow d
  · intro o v d hasVert ih
    cases o with
    | tree e cls child =>
      obtain ⟨cv, couts⟩ := child
      dsimp only at ih
      have hfin : ∀ (hasVert₁ : Bool) (origTstack : Nat),
          KeepsBelow d ((modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }) >>= fun _ =>
            walkTree (.node cv couts) (d + 1) >>= fun _ =>
            finishEdge v d (.tree e cls (.node cv couts)) origTstack hasVert₁) := fun _ _ =>
        KeepsBelow.bind ((KeepsAll.modify _ fun _ => rfl).keepsBelow d) fun _ =>
          KeepsBelow.bind (ih.mono (Nat.le_succ d)) fun _ => KeepsBelow.finishEdge _ _ _ _ _
      intro s d' hd'
      unfold walkOut
      dsimp only
      rw [bind_stackDir]
      refine bind_spec ((KeepsBelow.setStackDir (Nat.le_refl d) _).wp_of hd' rfl) fun _ s₁ h₁ => ?_
      split
      · refine bind_spec (((KeepsAll.pushVertTstack v d).keepsBelow d).wp_of hd' h₁) fun _ s₂ h₂ => ?_
        rw [bind_pure, bind_tstackSize]
        exact (hfin _ _).wp_of hd' h₂
      · rw [bind_pure, bind_tstackSize]
        exact (hfin _ _).wp_of hd' h₁
    | back e dest cls =>
      intro s d' hd'
      unfold walkOut
      dsimp only
      rw [bind_stackDir]
      refine bind_spec ((KeepsBelow.setStackDir (Nat.le_refl d) _).wp_of hd' rfl) fun _ s₁ h₁ => ?_
      split
      · refine bind_spec (((KeepsAll.pushVertTstack v d).keepsBelow d).wp_of hd' h₁) fun _ s₂ h₂ => ?_
        rw [bind_pure, bind_tstackSize]
        exact (KeepsBelow.finishEdge _ _ _ _ _).wp_of hd' h₂
      · rw [bind_pure, bind_tstackSize]
        exact (KeepsBelow.finishEdge _ _ _ _ _).wp_of hd' h₁
  · intro v d hasVert
    exact (KeepsAll.pure _).keepsBelow d
  · intro v d hasVert o rest ih₁ ih₂
    simp only [walkOuts]
    exact KeepsBelow.bind ih₁ fun hasVert' => ih₂ hasVert'

/-- Walking a subtree at depth `d` leaves `stackDir` below `d` unchanged. -/
theorem walkTree_stackDir_below (t : DfsTree) (d : Nat) (s : WalkState) :
    ∀ d', d' < d → ((walkTree t d).run s).2.stackDir[d']! = s.stackDir[d']! :=
  fun d' hd' => walk_keepsBelow.1 t d s d' hd'

end Spqr
