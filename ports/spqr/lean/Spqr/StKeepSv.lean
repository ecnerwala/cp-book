import Spqr.StFrame
/-!
# `stackVerts` frame facts

The only write to `stackVerts` is `walkTree`'s `stackVerts.set! d v`; `walkTree t d` leaves
`stackVerts[d']` unchanged for every `d' < d`.
-/

namespace Spqr
open WalkM

/-- `m` never changes `stackVerts`. -/
def KeepsSv (m : WalkM α) : Prop := ∀ s, wp m (fun _ s' => s'.stackVerts = s.stackVerts) s

/-- `m` leaves `stackVerts[d']` unchanged for every `d' < d`. -/
def KeepsSvBelow (d : Nat) (m : WalkM α) : Prop :=
  ∀ s d', d' < d → wp m (fun _ s' => s'.stackVerts[d']! = s.stackVerts[d']!) s

namespace KeepsSv

theorem keepsBelow {m : WalkM α} (h : KeepsSv m) (d : Nat) : KeepsSvBelow d m :=
  fun s d' _ => by
    have e : ((m.run s).2).stackVerts = s.stackVerts := h s
    exact congrArg (fun a : Array Nat => a[d']!) e

theorem pure (a : α) : KeepsSv (pure a : WalkM α) := fun _ => rfl

theorem bind {m : WalkM β} {f : β → WalkM α} (hm : KeepsSv m) (hf : ∀ b, KeepsSv (f b)) :
    KeepsSv (m >>= f) := by
  intro s
  show (((f (m.run s).1).run (m.run s).2)).2.stackVerts = s.stackVerts
  rw [hf _ _, hm s]

theorem map {m : WalkM β} (f : β → α) (hm : KeepsSv m) : KeepsSv (f <$> m) := hm

theorem ite (c : Prop) [Decidable c] {a b : WalkM α} (ha : KeepsSv a) (hb : KeepsSv b) :
    KeepsSv (if c then a else b) := by
  split <;> assumption

theorem modify (f : WalkState → WalkState) (hf : ∀ s, (f s).stackVerts = s.stackVerts) :
    KeepsSv (modify f : WalkM Unit) := hf

theorem get : KeepsSv (get : WalkM WalkState) := fun _ => rfl
theorem stackDir (d : Nat) : KeepsSv (stackDir d) := fun _ => rfl
theorem setStackDir (d : Nat) (b : Bool) : KeepsSv (setStackDir d b) := fun _ => rfl
theorem makeVs (a b : Nat) : KeepsSv (makeVs a b) := fun _ => rfl
theorem cur : KeepsSv cur := fun _ => rfl
theorem nxt : KeepsSv nxt := fun _ => rfl
theorem tstackSize : KeepsSv tstackSize := fun _ => rfl
theorem getItem (i : ItemId) : KeepsSv (getItem i) := fun _ => rfl
theorem modifyItem (i : ItemId) (f : Item → Item) : KeepsSv (modifyItem i f) := fun _ => rfl
theorem allocItem (ty : NodeType) : KeepsSv (allocItem ty) := fun _ => rfl
theorem modifyCur (f : TEntry → TEntry) : KeepsSv (modifyCur f) := fun _ => rfl
theorem modifyNxt (f : TEntry → TEntry) : KeepsSv (modifyNxt f) := fun _ => rfl
theorem popTstack : KeepsSv popTstack := fun _ => rfl
theorem pushTstack (v d : Nat) (i : ItemId) : KeepsSv (pushTstack v d i) := fun _ => rfl
theorem pushVertTstack (v d : Nat) : KeepsSv (pushVertTstack v d) := fun _ => rfl
theorem pushEdgeTstack (v d e : Nat) : KeepsSv (pushEdgeTstack v d e) := fun _ => rfl
theorem mergeTstackTops : KeepsSv mergeTstackTops := fun _ => rfl

theorem loop (n : Nat) {cond : WalkM Bool} {body : WalkM Unit} (hc : KeepsSv cond)
    (hb : KeepsSv body) : KeepsSv (loop n cond body) := by
  induction n with
  | zero => exact pure ()
  | succ n ih =>
    simp only [WalkM.loop]
    exact bind hc fun b => ite _ (bind hb fun _ => ih) (pure ())

theorem maybeUnwrapNxt (ty : NodeType) : KeepsSv (maybeUnwrapNxt ty) := by
  intro s; unfold WalkM.maybeUnwrapNxt; wp_simp
  repeat' split
  all_goals trivial

theorem finishTstackTop (i : ItemId) : KeepsSv (finishTstackTop i) := fun _ => rfl

theorem loop1Cond (d : Nat) : KeepsSv (loop1Cond d) := by
  unfold Spqr.loop1Cond
  exact bind tstackSize fun _ => bind nxt fun _ => pure _

theorem loop2Cond (fo : Nat) : KeepsSv (loop2Cond fo) := by
  unfold Spqr.loop2Cond
  exact bind cur fun _ => pure _

theorem loop3Cond (o : Nat) : KeepsSv (loop3Cond o) := by
  unfold Spqr.loop3Cond
  exact bind tstackSize fun _ => pure _

theorem mergeLate (d : Nat) : KeepsSv (mergeLate d) := by
  unfold Spqr.mergeLate
  refine bind get fun s => bind cur fun t => ite _ ?_ (pure _)
  exact bind tstackSize fun n => bind (loop n (loop2Cond _) mergeTstackTops) fun _ => pure _

theorem finishP (curV lowval : Nat) (isType1 : Bool) : KeepsSv (finishP curV lowval isType1) := by
  intro s; unfold Spqr.finishP Spqr.condP WalkM.maybeUnwrapNxt WalkM.finishTstackTop; wp_simp
  repeat' split
  all_goals trivial

theorem finishTail (curV d : Nat) (hasVert isSingle : Bool) :
    KeepsSv (finishTail curV d hasVert isSingle) := by
  intro s; unfold Spqr.finishTail; wp_simp
  repeat' split
  all_goals trivial

theorem finishRest (curV d lowval : Nat) (isType1 hasVert isSingle : Bool) :
    KeepsSv (finishRest curV d lowval isType1 hasVert isSingle) := by
  unfold Spqr.finishRest
  exact bind (finishP _ _ _) fun _ => finishTail _ _ _ _

theorem closeVertTail (curV : Nat) (edgeDir isSingle : Bool) {k : Bool → WalkM Bool}
    (hk : ∀ b, KeepsSv (k b)) (item : Option ItemId) :
    KeepsSv (closeVertTail curV edgeDir isSingle k item) := by
  unfold Spqr.closeVertTail
  refine bind mergeTstackTops fun _ => bind mergeTstackTops fun _ => bind (modifyCur _) fun _ => ?_
  cases item with
  | some i => exact bind (finishTstackTop i) fun _ => hk _
  | none => exact hk _

theorem closeVert (curV : Nat) (edgeDir isType1 : Bool) (origTstack : Nat) (isSingle : Bool)
    {k : Bool → WalkM Bool} (hk : ∀ b, KeepsSv (k b)) :
    KeepsSv (closeVert curV edgeDir isType1 origTstack isSingle k) := by
  unfold Spqr.closeVert
  refine ite _ ?_ ?_
  · exact bind tstackSize fun n => bind (loop n (loop3Cond _) mergeTstackTops) fun _ =>
      closeVertTail _ _ _ hk _
  · exact bind (map _ (maybeUnwrapNxt _)) fun i => closeVertTail _ _ _ hk _

theorem finishBoundary (curV d : Nat) (o : DfsOut) (qItem : ItemId) (hasVert : Bool) :
    KeepsSv (finishBoundary curV d o qItem hasVert) := by
  intro s; unfold Spqr.finishBoundary; wp_simp
  repeat' split
  all_goals trivial

theorem finishBack (curV d : Nat) (o : DfsOut) (hasVert : Bool) :
    KeepsSv (finishBack curV d o hasVert) := by
  unfold Spqr.finishBack
  exact bind (pushEdgeTstack _ _ _) fun _ => bind (modify _ fun _ => rfl) fun _ => finishRest _ _ _ _ _ _

theorem loop1Type (d : Nat) (edgeDir : Bool) : KeepsSv (loop1Type d edgeDir) := by
  intro s
  unfold Spqr.loop1Type; wp_simp
  repeat' split
  all_goals first | trivial | rfl

theorem loop1Body (d : Nat) (edgeDir : Bool) : KeepsSv (loop1Body d edgeDir) := by
  unfold Spqr.loop1Body
  refine bind (loop1Type d edgeDir) fun ty => bind (maybeUnwrapNxt ty) fun i => ?_
  exact bind mergeTstackTops fun _ => finishTstackTop i

theorem closeEars (nxtV d e : Nat) (edgeDir : Bool) : KeepsSv (closeEars nxtV d e edgeDir) := by
  unfold Spqr.closeEars
  refine bind (pushEdgeTstack _ _ _) fun _ => bind tstackSize fun n => ?_
  exact loop n (loop1Cond d) (loop1Body d edgeDir)

theorem finishTree (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert edgeDir : Bool) :
    KeepsSv (finishTree curV d o origTstack hasVert edgeDir) := by
  unfold Spqr.finishTree
  refine bind (closeEars _ _ _ _) fun _ => bind (mergeLate d) fun isSingle => ?_
  refine ite _ ?_ (finishRest _ _ _ _ _ _)
  exact closeVert _ _ _ _ _ fun b => finishRest _ _ _ _ _ _

theorem finishEdge (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool) :
    KeepsSv (finishEdge curV d o origTstack hasVert) := by
  rw [finishEdge_eq]
  unfold finishEdge'
  refine bind get fun s => bind (stackDir d) fun edgeDir => ?_
  refine ite _ (finishBoundary _ _ _ _ _) ?_
  refine bind (makeVs _ _) fun vs => bind (modifyItem _ _) fun _ => ?_
  exact ite _ (finishTree _ _ _ _ _ _) (finishBack _ _ _ _)

end KeepsSv

namespace KeepsSvBelow

theorem mono {m : WalkM α} {d d' : Nat} (h : KeepsSvBelow d m) (hd : d' ≤ d) : KeepsSvBelow d' m :=
  fun s i hi => h s i (Nat.lt_of_lt_of_le hi hd)

theorem wp_of {d d' : Nat} {m : WalkM α} (h : KeepsSvBelow d m) (hd' : d' < d) {s₀ s : WalkState}
    (hs : s.stackVerts[d']! = s₀.stackVerts[d']!) :
    wp m (fun _ s' => s'.stackVerts[d']! = s₀.stackVerts[d']!) s := by
  have := h s d' hd'
  show _ = _
  rw [this, hs]

theorem bind {d : Nat} {m : WalkM β} {f : β → WalkM α} (hm : KeepsSvBelow d m)
    (hf : ∀ b, KeepsSvBelow d (f b)) : KeepsSvBelow d (m >>= f) := by
  intro s d' hd'
  show ((f (m.run s).1).run (m.run s).2).2.stackVerts[d']! = s.stackVerts[d']!
  rw [hf _ _ d' hd', hm s d' hd']

theorem setSv {d : Nat} (v : Nat) :
    KeepsSvBelow d (modify fun s => { s with stackVerts := s.stackVerts.set! d v } : WalkM Unit) := by
  intro s j hj
  exact Array.getElem!_set!_ne _ _ _ _ (by omega)

end KeepsSvBelow

theorem walk_keepsSvBelow :
    (∀ (t : DfsTree) (d : Nat), KeepsSvBelow d (walkTree t d)) ∧
    (∀ (v d : Nat) (outs : List DfsOut) (hasVert : Bool), KeepsSvBelow (d + 1) (walkOuts v d outs hasVert)) ∧
    (∀ (v d : Nat) (o : DfsOut) (hasVert : Bool), KeepsSvBelow (d + 1) (walkOut v d o hasVert)) := by
  refine walkTree.mutual_induct _ _ _ ?_ ?_ ?_ ?_
  · intro d v outs ih
    simp only [walkTree]
    refine KeepsSvBelow.bind (KeepsSvBelow.setSv v) fun _ =>
      KeepsSvBelow.bind (ih.mono (Nat.le_succ d)) fun hasVert => ?_
    split
    · exact (KeepsSv.pure _).keepsBelow d
    · exact KeepsSvBelow.bind ((KeepsSv.setStackDir _ _).keepsBelow d) fun _ =>
        (KeepsSv.pushVertTstack _ _).keepsBelow d
  · intro o v d hasVert ih
    cases o with
    | tree e cls child =>
      obtain ⟨cv, couts⟩ := child
      dsimp only at ih
      have hfin : ∀ (hasVert₁ : Bool) (origTstack : Nat),
          KeepsSvBelow (d + 1) ((modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }) >>= fun _ =>
            walkTree (.node cv couts) (d + 1) >>= fun _ =>
            finishEdge v d (.tree e cls (.node cv couts)) origTstack hasVert₁) := fun _ _ =>
        KeepsSvBelow.bind ((KeepsSv.modify _ fun _ => rfl).keepsBelow _) fun _ =>
          KeepsSvBelow.bind ih fun _ => (KeepsSv.finishEdge _ _ _ _ _).keepsBelow _
      intro s d' hd'
      unfold walkOut
      dsimp only
      rw [bind_stackDir]
      refine bind_spec (((KeepsSv.setStackDir _ _).keepsBelow _).wp_of hd' rfl) fun _ s₁ h₁ => ?_
      split
      · refine bind_spec (((KeepsSv.pushVertTstack v d).keepsBelow _).wp_of hd' h₁) fun _ s₂ h₂ => ?_
        rw [bind_pure, bind_tstackSize]
        exact (hfin _ _).wp_of hd' h₂
      · rw [bind_pure, bind_tstackSize]
        exact (hfin _ _).wp_of hd' h₁
    | back e dest cls =>
      intro s d' hd'
      unfold walkOut
      dsimp only
      rw [bind_stackDir]
      refine bind_spec (((KeepsSv.setStackDir _ _).keepsBelow _).wp_of hd' rfl) fun _ s₁ h₁ => ?_
      split
      · refine bind_spec (((KeepsSv.pushVertTstack v d).keepsBelow _).wp_of hd' h₁) fun _ s₂ h₂ => ?_
        rw [bind_pure, bind_tstackSize]
        exact ((KeepsSv.finishEdge _ _ _ _ _).keepsBelow _).wp_of hd' h₂
      · rw [bind_pure, bind_tstackSize]
        exact ((KeepsSv.finishEdge _ _ _ _ _).keepsBelow _).wp_of hd' h₁
  · intro v d hasVert
    exact (KeepsSv.pure _).keepsBelow _
  · intro v d hasVert o rest ih₁ ih₂
    simp only [walkOuts]
    exact KeepsSvBelow.bind ih₁ fun hasVert' => ih₂ hasVert'

/-- Walking a subtree at depth `d` leaves `stackVerts` below `d` unchanged. -/
theorem walkTree_stackVerts_below (t : DfsTree) (d : Nat) (s : WalkState) :
    ∀ d', d' < d → ((walkTree t d).run s).2.stackVerts[d']! = s.stackVerts[d']! :=
  fun d' hd' => walk_keepsSvBelow.1 t d s d' hd'

theorem walkOut_stackVerts_le (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) :
    ∀ d', d' ≤ d → ((walkOut v d o hasVert).run s).2.stackVerts[d']! = s.stackVerts[d']! :=
  fun d' hd' => walk_keepsSvBelow.2.2 v d o hasVert s d' (Nat.lt_succ_of_le hd')

end Spqr
