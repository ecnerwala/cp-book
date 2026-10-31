import Spqr.Frame

/-! # `vStart` provenance of `finishEdge` (`PROOF.md` §7)

Every tstack entry after `finishEdge curV d o _ _` has `vStart = curV`, or `vStart = o.dest` (the
tree-edge entry of a returning tree edge), or inherits the `vStart` of an entry that was already on
the stack: the primitives only push entries at `curV`/`o.dest`, merge into the lower entry (keeping
its `vStart`), relabel the merged vertex ear to `curV`, pop, or touch spans only. -/

namespace Spqr
open WalkM

/-- Every tstack entry's `vStart` satisfies `V`. -/
abbrev VsIn (V : Nat → Prop) (s : WalkState) : Prop := ∀ t ∈ s.tstack, V t.vStart

namespace VsIn
variable {V : Nat → Prop}

theorem tail {ts : List TEntry} (h : ∀ t ∈ ts, V t.vStart) : ∀ t ∈ ts.tail, V t.vStart :=
  fun t ht => h t (List.mem_of_mem_tail ht)

theorem cons {e : TEntry} {ts : List TEntry} (hv : V e.vStart) (h : ∀ t ∈ ts, V t.vStart) :
    ∀ t ∈ e :: ts, V t.vStart := by
  intro t ht
  rcases List.mem_cons.mp ht with rfl | ht
  · exact hv
  · exact h t ht

theorem merge : ∀ {ts : List TEntry}, (∀ t ∈ ts, V t.vStart) → ∀ t ∈ mergeTop ts, V t.vStart
  | [], _, t, ht => by simp [mergeTop] at ht
  | [_], _, t, ht => by simp [mergeTop] at ht
  | b :: a :: rest, h, t, ht => by
    simp only [mergeTop, List.mem_cons] at ht
    rcases ht with rfl | ht
    · exact h a (by simp)
    · exact h t (by simp [ht])

theorem modifyCur (f : TEntry → TEntry) {ts : List TEntry} (h : ∀ t ∈ ts, V t.vStart)
    (hf : ∀ t ∈ ts, V (f t).vStart) :
    ∀ t ∈ (match ts with | a :: rest => f a :: rest | [] => []), V t.vStart := by
  cases ts with
  | nil => simp
  | cons a rest =>
    intro t ht
    rcases List.mem_cons.mp ht with rfl | ht
    · exact hf _ (List.mem_cons_self ..)
    · exact h _ (List.mem_cons_of_mem _ ht)

theorem modifyNxt (f : TEntry → TEntry) (hf : ∀ t, (f t).vStart = t.vStart) {ts : List TEntry}
    (h : ∀ t ∈ ts, V t.vStart) :
    ∀ t ∈ (match ts with | a :: b :: rest => a :: f b :: rest | l => l), V t.vStart := by
  match ts with
  | [] => exact h
  | [_] => exact h
  | a :: b :: rest =>
    intro t ht
    rcases List.mem_cons.mp ht with rfl | ht
    · exact h _ (List.mem_cons_self ..)
    rcases List.mem_cons.mp ht with rfl | ht
    · rw [hf]; exact h _ (by simp)
    · exact h _ (by simp [ht])

end VsIn

section
variable {V : Nat → Prop} {s : WalkState}

theorem vsIn_mergeTstackTops (h : VsIn V s) : wp mergeTstackTops (fun _ s' => VsIn V s') s := by
  simp only [wp_mergeTstackTops]; exact VsIn.merge h

theorem vsIn_finishTstackTop (item : ItemId) (h : VsIn V s) :
    wp (finishTstackTop item) (fun _ s' => VsIn V s') s := by
  unfold finishTstackTop
  simp only [wp_bind, wp_cur, wp_stackDir, wp_makeVs, wp_modifyItem, wp_modifyCur]
  exact VsIn.modifyCur _ h (by intro t ht; exact h t ht)

theorem vsIn_maybeUnwrapNxt (ty : NodeType) (h : VsIn V s) :
    wp (maybeUnwrapNxt ty) (fun _ s' => VsIn V s') s := by
  unfold maybeUnwrapNxt
  simp only [wp_bind, wp_ite, wp_get, wp_allocItem, wp_pure, wp_nxt, wp_stackDir, wp_getItem,
    wp_modifyNxt]
  split
  · exact h
  · split
    · exact VsIn.modifyNxt _ (by intro _; rfl) h
    · exact h

theorem vsIn_loop1Body (d : Nat) (edgeDir : Bool) (h : VsIn V s) :
    wp (loop1Body d edgeDir) (fun _ s' => VsIn V s') s := by
  unfold loop1Body loop1Type
  simp only [wp_bind, wp_ite, wp_nxt, wp_cur, wp_setStackDir, wp_mergeTstackTops, wp_pure]
  split
  · exact wp_mono _ (vsIn_maybeUnwrapNxt .S (VsIn.merge h)) fun item s₁ h₁ =>
      vsIn_finishTstackTop item (VsIn.merge h₁)
  · split <;> exact wp_mono _ (vsIn_maybeUnwrapNxt _ h) fun item s₁ h₁ =>
      vsIn_finishTstackTop item (VsIn.merge h₁)

theorem vsIn_loop (n : Nat) (cond : WalkM Bool) (body : WalkM Unit) (hc : ∀ s, (cond.run s).2 = s)
    (hb : ∀ s, VsIn V s → wp body (fun _ s' => VsIn V s') s) (h : VsIn V s) :
    wp (loop n cond body) (fun _ s' => VsIn V s') s :=
  wp_loop (fun s' => VsIn V s') n cond body
    (fun s₁ h₁ => by show VsIn V (cond.run s₁).2; rw [hc]; exact h₁) hb h fun _ h => h

theorem vsIn_finishRest (curV d lowval : Nat) (isType1 hasVert isSingle : Bool) (hV : V curV)
    (h : VsIn V s) :
    wp (finishRest curV d lowval isType1 hasVert isSingle) (fun _ s' => VsIn V s') s := by
  unfold finishRest finishP finishTail condP
  simp only [wp_bind, wp_ite, wp_tstackSize, wp_nxt, wp_pure, wp_pushVertTstack, wp_mergeTstackTops]
  split
  · refine wp_mono _ (vsIn_maybeUnwrapNxt .P h) fun item s₁ h₁ => ?_
    refine wp_mono _ (vsIn_finishTstackTop item (VsIn.merge h₁)) fun _ s₂ h₂ => ?_
    try simp only [wp_bind, wp_ite, wp_pure, wp_pushVertTstack, wp_mergeTstackTops]
    split <;> (try split) <;>
      first | exact h₂ | exact VsIn.merge (VsIn.cons hV h₂) | exact VsIn.cons hV h₂
  · split <;> (try split) <;>
      first | exact h | exact VsIn.merge (VsIn.cons hV h) | exact VsIn.cons hV h

theorem vsIn_closeVertTail (curV : Nat) (edgeDir isSingle : Bool) (k : Bool → WalkM Bool)
    (item : Option ItemId) (hV : V curV)
    (hk : ∀ b s, VsIn V s → wp (k b) (fun _ s' => VsIn V s') s) (h : VsIn V s) :
    wp (closeVertTail curV edgeDir isSingle k item) (fun _ s' => VsIn V s') s := by
  unfold closeVertTail
  cases item with
  | none =>
    simp only [wp_bind, wp_mergeTstackTops, wp_modifyCur]
    exact hk _ _ (VsIn.modifyCur _ (VsIn.merge (VsIn.merge h)) (by intros; exact hV))
  | some item =>
    simp only [wp_bind, wp_mergeTstackTops, wp_modifyCur]
    exact wp_mono _ (vsIn_finishTstackTop item
      (VsIn.modifyCur _ (VsIn.merge (VsIn.merge h)) (by intros; exact hV))) fun _ s₂ h₂ => hk _ _ h₂

theorem vsIn_closeVert (curV : Nat) (edgeDir isType1 : Bool) (origTstack : Nat) (isSingle : Bool)
    (k : Bool → WalkM Bool) (hV : V curV)
    (hk : ∀ b s, VsIn V s → wp (k b) (fun _ s' => VsIn V s') s) (h : VsIn V s) :
    wp (closeVert curV edgeDir isType1 origTstack isSingle k) (fun _ s' => VsIn V s') s := by
  unfold closeVert
  simp only [wp_bind, wp_ite, wp_tstackSize, wp_map]
  split
  · refine wp_mono _ (vsIn_loop _ _ _ (fun s => by rw [run_loop3Cond])
      (fun s h => vsIn_mergeTstackTops h) h) fun _ s₃ h₃ => ?_
    exact vsIn_closeVertTail _ _ _ _ _ hV hk h₃
  · refine wp_mono _ (vsIn_maybeUnwrapNxt _ h) fun item s₃ h₃ => ?_
    exact vsIn_closeVertTail _ _ _ _ _ hV hk h₃

theorem vsIn_finishTree (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert edgeDir : Bool)
    (hV : V curV) (hD : V o.dest) (h : VsIn V s) :
    wp (finishTree curV d o origTstack hasVert edgeDir) (fun _ s' => VsIn V s') s := by
  unfold finishTree closeEars mergeLate
  simp only [wp_bind, wp_pushEdgeTstack, wp_tstackSize, wp_get, wp_cur, wp_ite, wp_pure]
  refine wp_mono _ (vsIn_loop _ _ _ (fun s => by rw [run_loop1Cond])
    (fun s h => vsIn_loop1Body d _ h) (VsIn.cons hD h)) fun _ s₁ h₁ => ?_
  have fin : ∀ (b : Bool) (s₂ : WalkState), VsIn V s₂ →
      if hasVert = true then
        wp (closeVert curV edgeDir o.cls.isType1 origTstack b
          (finishRest curV d (o.cls.lowval d) o.cls.isType1 hasVert)) (fun _ s' => VsIn V s') s₂
      else wp (finishRest curV d (o.cls.lowval d) o.cls.isType1 hasVert b)
        (fun _ s' => VsIn V s') s₂ := by
    intro b s₂ h₂
    split
    · exact vsIn_closeVert _ _ _ _ _ _ hV (fun b s h => vsIn_finishRest _ _ _ _ _ _ hV h) h₂
    · exact vsIn_finishRest _ _ _ _ _ _ hV h₂
  split
  · exact wp_mono _ (vsIn_loop _ _ _ (fun s => by rw [run_loop2Cond])
      (fun s h => vsIn_mergeTstackTops h) h₁) fun _ s₂ h₂ => fin false s₂ h₂
  · exact fin true s₁ h₁

theorem vsIn_finishBoundary (curV d : Nat) (o : DfsOut) (qItem : ItemId) (hasVert : Bool)
    (h : VsIn V s) : wp (finishBoundary curV d o qItem hasVert) (fun _ s' => VsIn V s') s := by
  unfold finishBoundary
  simp only [wp_bind, wp_modifyItem, wp_modify, wp_ite, wp_allocItem, wp_makeVs, wp_popTstack,
    wp_pure]
  split <;> (try split) <;> first | exact h | exact VsIn.tail h | exact VsIn.tail (VsIn.tail h)

theorem vsIn_finishEdge (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool)
    (hV : V curV) (hD : o.cls.isTree = true → o.cls.lowval d < d → V o.dest) (h : VsIn V s) :
    wp (finishEdge curV d o origTstack hasVert) (fun _ s' => VsIn V s') s := by
  by_cases hge : o.cls.lowval d ≥ d
  · simp only [finishEdge_eq, finishEdge', wp_bind, wp_get, wp_stackDir, hge, ↓reduceIte]
    exact vsIn_finishBoundary _ _ _ _ _ h
  · by_cases ht : o.cls.isTree = true
    · simp only [finishEdge_eq, finishEdge', wp_bind, wp_get, wp_stackDir, hge, ht, ↓reduceIte,
        wp_makeVs, wp_modifyItem]
      exact vsIn_finishTree _ _ _ _ _ _ hV (hD ht (Nat.lt_of_not_le hge)) h
    · have ht' : o.cls.isTree = false := Bool.eq_false_iff.2 ht
      simp only [finishEdge_eq, finishEdge', finishBack, wp_bind, wp_get, wp_stackDir, hge, ht',
        Bool.false_eq_true, ↓reduceIte, wp_makeVs, wp_modifyItem, wp_pushEdgeTstack, wp_modify]
      exact vsIn_finishRest _ _ _ _ _ _ hV (VsIn.cons hV h)

end

/-- Every tstack entry after `finishEdge` starts at `curV`, at the child `o.dest` of a returning
tree edge, or at the `vStart` of an entry that was already on the stack. -/
theorem finishEdge_vStart (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool)
    (s : WalkState) :
    ∀ t ∈ ((finishEdge curV d o origTstack hasVert).run s).2.tstack,
      t.vStart = curV ∨ (o.cls.isTree = true ∧ o.cls.lowval d < d ∧ t.vStart = o.dest) ∨
        ∃ t₀ ∈ s.tstack, t.vStart = t₀.vStart :=
  vsIn_finishEdge (V := fun x => x = curV ∨ (o.cls.isTree = true ∧ o.cls.lowval d < d ∧ x = o.dest) ∨
      ∃ t₀ ∈ s.tstack, x = t₀.vStart)
    curV d o origTstack hasVert (Or.inl rfl) (fun ht hlt => Or.inr (Or.inl ⟨ht, hlt, rfl⟩))
    fun t ht => Or.inr (Or.inr ⟨t, ht, rfl⟩)

end Spqr
