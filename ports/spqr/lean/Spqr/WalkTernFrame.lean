import Spqr.EarFrame

/-!
# `ternarize` is never written

Frame of the `ternarize` flag through `walkOutPre` and `finishEdge`, by the same traversal as
`EarFrame.lean`'s `keep_*` (no hypotheses: every primitive is a record update that leaves the
flag alone). The backbone (`WalkBackbone.lean`) carries `s.ternarize = G.tern` from it, and
`walk_ternarize` (`WalkItemsWF.lean`) is its public form.
-/

namespace Spqr.WalkState
open WalkM

variable {b : Bool} {s : WalkState}

theorem tern_mergeTstackTops (h₀ : s.ternarize = b) : wp mergeTstackTops (fun _ s' => s'.ternarize = b) s := by
  simp only [wp_mergeTstackTops]; exact h₀

theorem tern_finishTstackTop (item : ItemId) (h₀ : s.ternarize = b) :
    wp (finishTstackTop item) (fun _ s' => s'.ternarize = b) s := by
  unfold finishTstackTop
  simp only [wp_bind, wp_cur, wp_stackDir, wp_makeVs, wp_modifyItem, wp_modifyCur]
  exact h₀

theorem tern_maybeUnwrapNxt (ty : NodeType) (h₀ : s.ternarize = b) :
    wp (maybeUnwrapNxt ty) (fun _ s' => s'.ternarize = b) s := by
  unfold maybeUnwrapNxt
  simp only [wp_bind, wp_ite, wp_get, wp_allocItem, wp_pure, wp_nxt, wp_stackDir, wp_getItem,
    wp_modifyNxt]
  split
  · exact h₀
  · split
    · exact h₀
    · exact h₀

theorem tern_loop1Body (d : Nat) (edgeDir : Bool) (h₀ : s.ternarize = b) :
    wp (loop1Body d edgeDir) (fun _ s' => s'.ternarize = b) s := by
  unfold loop1Body loop1Type
  simp only [wp_bind, wp_ite, wp_nxt, wp_cur, wp_setStackDir, wp_mergeTstackTops, wp_pure]
  split
  · exact wp_mono _ (tern_maybeUnwrapNxt .S h₀) fun item s₁ h₁ => tern_finishTstackTop item h₁
  · split <;> exact wp_mono _ (tern_maybeUnwrapNxt _ h₀) fun item s₁ h₁ => tern_finishTstackTop item h₁

theorem tern_loop (n : Nat) (cond : WalkM Bool) (body : WalkM Unit) (hc : ∀ s, (cond.run s).2 = s)
    (hb : ∀ s, s.ternarize = b → wp body (fun _ s' => s'.ternarize = b) s) (h₀ : s.ternarize = b) :
    wp (loop n cond body) (fun _ s' => s'.ternarize = b) s :=
  wp_loop (fun s' => s'.ternarize = b) n cond body
    (fun s₁ h₁ => by show (cond.run s₁).2.ternarize = b; rw [hc]; exact h₁) hb h₀ fun _ h => h

theorem tern_finishRest (curV d lowval : Nat) (isType1 hasVert isSingle : Bool) (h₀ : s.ternarize = b) :
    wp (finishRest curV d lowval isType1 hasVert isSingle) (fun _ s' => s'.ternarize = b) s := by
  unfold finishRest finishP finishTail condP
  simp only [wp_bind, wp_ite, wp_tstackSize, wp_nxt, wp_pure, wp_pushVertTstack, wp_mergeTstackTops]
  split
  · refine wp_mono _ (tern_maybeUnwrapNxt .P h₀) fun item s₁ h₁ => ?_
    refine wp_mono _ (tern_finishTstackTop item h₁) fun _ s₂ h₂ => ?_
    try simp only [wp_bind, wp_ite, wp_pure, wp_pushVertTstack, wp_mergeTstackTops]
    split <;> (try split) <;> exact h₂
  · split <;> (try split) <;> exact h₀

theorem tern_closeVertTail (curV : Nat) (edgeDir isSingle : Bool) (k : Bool → WalkM Bool)
    (item : Option ItemId) (hk : ∀ b' s, s.ternarize = b → wp (k b') (fun _ s' => s'.ternarize = b) s)
    (h₀ : s.ternarize = b) :
    wp (closeVertTail curV edgeDir isSingle k item) (fun _ s' => s'.ternarize = b) s := by
  unfold closeVertTail
  cases item with
  | none =>
    simp only [wp_bind, wp_mergeTstackTops, wp_modifyCur]
    exact hk _ _ h₀
  | some item =>
    simp only [wp_bind, wp_mergeTstackTops, wp_modifyCur]
    exact wp_mono _ (tern_finishTstackTop item h₀) fun _ s₂ h₂ => hk _ _ h₂

theorem tern_closeVert (curV : Nat) (edgeDir isType1 : Bool) (origTstack : Nat) (isSingle : Bool)
    (k : Bool → WalkM Bool) (hk : ∀ b' s, s.ternarize = b → wp (k b') (fun _ s' => s'.ternarize = b) s)
    (h₀ : s.ternarize = b) :
    wp (closeVert curV edgeDir isType1 origTstack isSingle k) (fun _ s' => s'.ternarize = b) s := by
  unfold closeVert
  simp only [wp_bind, wp_ite, wp_tstackSize, wp_map]
  split
  · refine wp_mono _ (tern_loop _ _ _ (fun s => by rw [run_loop3Cond])
      (fun s h => tern_mergeTstackTops h) h₀) fun _ s₃ h₃ => ?_
    exact tern_closeVertTail _ _ _ _ _ hk h₃
  · refine wp_mono _ (tern_maybeUnwrapNxt _ h₀) fun item s₃ h₃ => ?_
    exact tern_closeVertTail _ _ _ _ _ hk h₃

theorem tern_finishTree (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert edgeDir : Bool)
    (h₀ : s.ternarize = b) :
    wp (finishTree curV d o origTstack hasVert edgeDir) (fun _ s' => s'.ternarize = b) s := by
  unfold finishTree closeEars mergeLate
  simp only [wp_bind, wp_pushEdgeTstack, wp_tstackSize, wp_get, wp_cur, wp_ite, wp_pure]
  refine wp_mono _ (tern_loop _ _ _ (fun s => by rw [run_loop1Cond])
    (fun s h => tern_loop1Body d _ h) h₀) fun _ s₁ h₁ => ?_
  have fin : ∀ (b' : Bool) (s₂ : WalkState), s₂.ternarize = b →
      if hasVert = true then
        wp (closeVert curV edgeDir o.cls.isType1 origTstack b'
          (finishRest curV d (o.cls.lowval d) o.cls.isType1 hasVert)) (fun _ s' => s'.ternarize = b) s₂
      else wp (finishRest curV d (o.cls.lowval d) o.cls.isType1 hasVert b') (fun _ s' => s'.ternarize = b) s₂ := by
    intro b' s₂ h₂
    split
    · exact tern_closeVert _ _ _ _ _ _ (fun b s h => tern_finishRest _ _ _ _ _ _ h) h₂
    · exact tern_finishRest _ _ _ _ _ _ h₂
  split
  · exact wp_mono _ (tern_loop _ _ _ (fun s => by rw [run_loop2Cond])
      (fun s h => tern_mergeTstackTops h) h₁) fun _ s₂ h₂ => fin false s₂ h₂
  · exact fin true s₁ h₁

theorem tern_finishBoundary (curV d : Nat) (o : DfsOut) (qItem : ItemId) (hasVert : Bool)
    (h₀ : s.ternarize = b) :
    wp (finishBoundary curV d o qItem hasVert) (fun _ s' => s'.ternarize = b) s := by
  unfold finishBoundary
  simp only [wp_bind, wp_modifyItem, wp_modify, wp_ite, wp_allocItem, wp_makeVs, wp_popTstack,
    wp_pure]
  split <;> (try split) <;> exact h₀

theorem tern_finishEdge (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool)
    (h₀ : s.ternarize = b) :
    wp (finishEdge curV d o origTstack hasVert) (fun _ s' => s'.ternarize = b) s := by
  by_cases hge : o.cls.lowval d ≥ d
  · simp only [finishEdge_eq, finishEdge', wp_bind, wp_get, wp_stackDir, hge, ↓reduceIte]
    exact tern_finishBoundary _ _ _ _ _ h₀
  · by_cases ht : o.cls.isTree = true
    · simp only [finishEdge_eq, finishEdge', wp_bind, wp_get, wp_stackDir, hge, ht, ↓reduceIte,
        wp_makeVs, wp_modifyItem]
      exact tern_finishTree _ _ _ _ _ _ h₀
    · have ht' : o.cls.isTree = false := Bool.eq_false_iff.2 ht
      simp only [finishEdge_eq, finishEdge', finishBack, wp_bind, wp_get, wp_stackDir, hge, ht',
        Bool.false_eq_true, ↓reduceIte, wp_makeVs, wp_modifyItem, wp_pushEdgeTstack, wp_modify]
      exact tern_finishRest _ _ _ _ _ _ h₀

theorem tern_walkOutPre (v d : Nat) (o : DfsOut) (hasVert : Bool) (h₀ : s.ternarize = b) :
    wp (walkOutPre v d o hasVert) (fun _ s' => s'.ternarize = b) s := by
  unfold walkOutPre
  simp only [wp_bind, wp_stackDir, wp_setStackDir, wp_ite, wp_pushVertTstack, wp_pure]
  split <;> exact h₀

end Spqr.WalkState
