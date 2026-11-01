import Spqr.EarSides

/-!
# `RootOK` at the forest level

Root-level facts for `walk_sides`: the walk keeps `g` and `stackDir.size` fixed (`Fr`), every
root out-edge is a boundary that pops exactly the child's entries (`finishEdge_root`), so after
`walkTree t 0` the stack is the root's single vertex entry on side 2 (`walkTree_rootOK`); and the
forest threading `sidesForest_of_roots`/`walk_sides_of_roots`: `SidesForest` from the per-root
inputs `RootsBook` (`GuardsTree`/`BookTree`/`Inv' 0` at every root).
-/

namespace Spqr
open WalkM

/-! ### `RootOK` after `walkTree t 0`

At depth 0 every out-edge is a block boundary (`lowval ≥ 0`): `finishEdge` pops exactly the
entries the child's walk pushed (`bd_bridge`/`bd_comp`) and never pushes, so the stack is empty
and `hasVert = false` when the root's own vertex entry is pushed with `stackDir[0] = true`. -/

namespace WalkState

/-- The graph and the `stackDir` size are untouched by the walk. -/
def Fr (s s' : WalkState) : Prop := s'.g = s.g ∧ s'.stackDir.size = s.stackDir.size

theorem Fr.refl : Fr s s := ⟨rfl, rfl⟩
theorem Fr.trans {s₁ s₂ : WalkState} (h₁ : Fr s s₁) (h₂ : Fr s₁ s₂) : Fr s s₂ :=
  ⟨h₂.1.trans h₁.1, h₂.2.trans h₁.2⟩
theorem Fr.setStackDir (d : Nat) (b : Bool) :
    Fr s { s with stackDir := s.stackDir.set! d b } := ⟨rfl, by simp⟩

theorem fr_maybeUnwrapNxt (ty : NodeType) : wp (maybeUnwrapNxt ty) (fun _ s' => Fr s s') s := by
  unfold maybeUnwrapNxt
  simp only [wp_bind, wp_ite, wp_get, wp_allocItem, wp_pure, wp_nxt, wp_stackDir, wp_getItem,
    wp_modifyNxt]
  split <;> (try split) <;> exact ⟨rfl, rfl⟩

theorem fr_finishTstackTop (item : ItemId) : wp (finishTstackTop item) (fun _ s' => Fr s s') s := by
  unfold finishTstackTop
  simp only [wp_bind, wp_cur, wp_stackDir, wp_makeVs, wp_modifyItem, wp_modifyCur]
  exact ⟨rfl, rfl⟩

theorem fr_mergeTstackTops : wp mergeTstackTops (fun _ s' => Fr s s') s := by
  simp only [wp_mergeTstackTops]; exact ⟨rfl, rfl⟩

theorem fr_loop1Body (d : Nat) (edgeDir : Bool) :
    wp (loop1Body d edgeDir) (fun _ s' => Fr s s') s := by
  unfold loop1Body loop1Type
  simp only [wp_bind, wp_ite, wp_nxt, wp_cur, wp_setStackDir, wp_mergeTstackTops, wp_pure]
  split
  · refine wp_mono _ (fr_maybeUnwrapNxt _) fun item s₁ h₁ => ?_
    exact wp_mono _ (fr_finishTstackTop item) fun _ s₂ h₂ =>
      ((Fr.setStackDir _ _).trans h₁).trans h₂
  · split <;>
    · refine wp_mono _ (fr_maybeUnwrapNxt _) fun item s₁ h₁ => ?_
      exact wp_mono _ (fr_finishTstackTop item) fun _ s₂ h₂ => h₁.trans h₂

theorem fr_loop (n : Nat) (cond : WalkM Bool) (body : WalkM Unit)
    (hc : ∀ s, (cond.run s).2 = s) (hb : ∀ s, wp body (fun _ s' => Fr s s') s) :
    wp (loop n cond body) (fun _ s' => Fr s s') s :=
  wp_loop (fun s' => Fr s s') n cond body
    (fun s₁ h₁ => by show Fr s (cond.run s₁).2; rw [hc]; exact h₁)
    (fun s₁ h₁ => wp_mono _ (hb s₁) fun _ s₂ h₂ => h₁.trans h₂) Fr.refl fun _ h => h

theorem fr_closeVertTail (curV : Nat) (edgeDir isSingle : Bool) (k : Bool → WalkM Bool)
    (item : Option ItemId) (hk : ∀ b s, wp (k b) (fun _ s' => Fr s s') s) :
    wp (closeVertTail curV edgeDir isSingle k item) (fun _ s' => Fr s s') s := by
  unfold closeVertTail
  cases item with
  | none => simp only [wp_bind, wp_mergeTstackTops, wp_modifyCur]; exact hk _ _
  | some item =>
    simp only [wp_bind, wp_mergeTstackTops, wp_modifyCur]
    exact wp_mono _ (fr_finishTstackTop item) fun _ s₂ h₂ => wp_mono _ (hk _ _) fun _ s₃ h₃ => h₂.trans h₃

theorem fr_finishRest (curV d lowval : Nat) (isType1 hasVert isSingle : Bool) :
    wp (finishRest curV d lowval isType1 hasVert isSingle) (fun _ s' => Fr s s') s := by
  unfold finishRest finishP finishTail condP
  simp only [wp_bind, wp_ite, wp_tstackSize, wp_nxt, wp_pure, wp_pushVertTstack, wp_mergeTstackTops]
  split
  · refine wp_mono _ (fr_maybeUnwrapNxt _) fun item s₁ h₁ => ?_
    refine wp_mono _ (fr_finishTstackTop item) fun _ s₂ h₂ => ?_
    try simp only [wp_bind, wp_ite, wp_pure, wp_pushVertTstack, wp_mergeTstackTops]
    split <;> (try split) <;> exact h₁.trans h₂
  · split <;> (try split) <;> exact ⟨rfl, rfl⟩

theorem fr_closeVert (curV : Nat) (edgeDir isType1 : Bool) (origTstack : Nat) (isSingle : Bool)
    (k : Bool → WalkM Bool) (hk : ∀ b s, wp (k b) (fun _ s' => Fr s s') s) :
    wp (closeVert curV edgeDir isType1 origTstack isSingle k) (fun _ s' => Fr s s') s := by
  unfold closeVert
  simp only [wp_bind, wp_ite, wp_tstackSize, wp_map]
  split
  · refine wp_mono _ (fr_loop _ _ _ (fun s => by rw [run_loop3Cond]) fun s => fr_mergeTstackTops)
      fun _ s₃ h₃ => ?_
    exact wp_mono _ (fr_closeVertTail _ _ _ _ _ hk) fun _ s₄ h₄ => h₃.trans h₄
  · refine wp_mono _ (fr_maybeUnwrapNxt _) fun item s₃ h₃ => ?_
    exact wp_mono _ (fr_closeVertTail _ _ _ _ _ hk) fun _ s₄ h₄ => h₃.trans h₄

theorem fr_finishTree (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert edgeDir : Bool) :
    wp (finishTree curV d o origTstack hasVert edgeDir) (fun _ s' => Fr s s') s := by
  unfold finishTree closeEars mergeLate
  simp only [wp_bind, wp_pushEdgeTstack, wp_tstackSize, wp_get, wp_cur, wp_ite, wp_pure]
  refine wp_mono _ (fr_loop _ _ _ (fun s => by rw [run_loop1Cond]) fun s => fr_loop1Body d _)
    fun _ s₁ h₁ => ?_
  have fin : ∀ (b : Bool) (s₂ : WalkState), Fr s s₂ →
      if hasVert = true then
        wp (closeVert curV edgeDir o.cls.isType1 origTstack b
          (finishRest curV d (o.cls.lowval d) o.cls.isType1 hasVert)) (fun _ s' => Fr s s') s₂
      else wp (finishRest curV d (o.cls.lowval d) o.cls.isType1 hasVert b) (fun _ s' => Fr s s') s₂ := by
    intro b s₂ h₂
    split
    · exact wp_mono _ (fr_closeVert _ _ _ _ _ _ fun b s => fr_finishRest ..) fun _ s₃ h₃ => h₂.trans h₃
    · exact wp_mono _ (fr_finishRest ..) fun _ s₃ h₃ => h₂.trans h₃
  split
  · exact wp_mono _ (fr_loop _ _ _ (fun s => by rw [run_loop2Cond]) fun s => fr_mergeTstackTops)
      fun _ s₂ h₂ => fin false s₂ (h₁.trans h₂)
  · exact fin true s₁ h₁

theorem fr_finishBoundary (curV d : Nat) (o : DfsOut) (qItem : ItemId) (hasVert : Bool) :
    wp (finishBoundary curV d o qItem hasVert) (fun _ s' => Fr s s') s := by
  unfold finishBoundary
  simp only [wp_bind, wp_modifyItem, wp_modify, wp_ite, wp_allocItem, wp_makeVs, wp_popTstack,
    wp_pure]
  split <;> (try split) <;> exact ⟨rfl, rfl⟩

theorem fr_finishEdge (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool) :
    wp (finishEdge curV d o origTstack hasVert) (fun _ s' => Fr s s') s := by
  by_cases hge : o.cls.lowval d ≥ d
  · simp only [finishEdge_eq, finishEdge', wp_bind, wp_get, wp_stackDir, hge, ↓reduceIte]
    exact fr_finishBoundary ..
  · by_cases ht : o.cls.isTree = true
    · simp only [finishEdge_eq, finishEdge', wp_bind, wp_get, wp_stackDir, hge, ht, ↓reduceIte,
        wp_makeVs, wp_modifyItem]
      exact wp_mono _ (fr_finishTree ..) fun _ s₂ h₂ => Fr.trans ⟨rfl, rfl⟩ h₂
    · have ht' : o.cls.isTree = false := Bool.eq_false_iff.2 ht
      simp only [finishEdge_eq, finishEdge', finishBack, wp_bind, wp_get, wp_stackDir, hge, ht',
        Bool.false_eq_true, ↓reduceIte, wp_makeVs, wp_modifyItem, wp_pushEdgeTstack, wp_modify]
      exact wp_mono _ (fr_finishRest ..) fun _ s₂ h₂ => Fr.trans ⟨rfl, rfl⟩ h₂

abbrev FrTreeP (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  wp (walkTree t d) (fun _ s' => Fr s s') s
abbrev FrOutsP (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  wp (walkOuts v d outs hasVert) (fun _ s' => Fr s s') s
abbrev FrOutP (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  wp (walkOut v d o hasVert) (fun _ s' => Fr s s') s

mutual
theorem frTreeP : ∀ (t : DfsTree) (d : Nat) (s : WalkState), FrTreeP t d s
  | .node v outs, d, s => by
    unfold FrTreeP walkTree
    simp only [wp_bind, wp_modify]
    refine wp_mono _ (frOutsP v d outs false _) fun hv s' h' => ?_
    cases hv
    · simp only [Bool.false_eq_true, ↓reduceIte, wp_bind, wp_setStackDir, wp_pushVertTstack]
      exact h'.trans (Fr.setStackDir _ _)
    · exact h'

theorem frOutsP : ∀ (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState),
    FrOutsP v d outs hasVert s
  | v, d, [], hasVert, s => by unfold FrOutsP walkOuts; simp only [wp_pure]; exact Fr.refl
  | v, d, o :: rest, hasVert, s => by
    unfold FrOutsP walkOuts
    simp only [wp_bind]
    exact wp_mono _ (frOutP v d o hasVert s) fun hv' s' h' =>
      wp_mono _ (frOutsP v d rest hv' s') fun _ s'' h'' => h'.trans h''

theorem frOutP : ∀ (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState), FrOutP v d o hasVert s
  | v, d, o, hasVert, s => by
    unfold FrOutP
    rw [walkOut_eq, wp_bind]
    have pre : wp (walkOutPre v d o hasVert) (fun _ s' => Fr s s') s := by
      unfold walkOutPre
      simp only [wp_bind, wp_stackDir, wp_setStackDir, wp_ite, wp_pushVertTstack, wp_pure]
      split <;> exact Fr.setStackDir _ _
    refine wp_mono _ pre fun hv' s₁ h₁ => ?_
    unfold walkOutRest
    rw [wp_bind, wp_tstackSize]
    cases o with
    | back e cls dest =>
      exact wp_mono _ (fr_finishEdge ..) fun _ s₂ h₂ => h₁.trans h₂
    | tree e cls child =>
      simp only [wp_bind, wp_modify]
      refine wp_mono _ (frTreeP child (d + 1) _) fun _ s₃ h₃ => ?_
      exact wp_mono _ (fr_finishEdge ..) fun _ s₄ h₄ => (h₁.trans h₃).trans h₄
end

/-- `finishEdge` of a root out-edge: the boundary branch pops exactly the child's entries. -/
theorem finishEdge_root {v : Nat} {o : DfsOut} (hb : FinishBook v 0 o 0 false s)
    (hback : o.cls.isTree = false → s.tstack = []) :
    wp (finishEdge v 0 o 0 false) (fun hv s' => hv = false ∧ s'.tstack = []) s := by
  obtain ⟨sub, base, hlen, hE⟩ := hb.ear
  have hbase : base = [] := List.length_eq_zero_iff.mp hlen
  subst hbase
  have hts : s.tstack = sub := by simpa using hE.tstack
  have hge : o.cls.lowval 0 ≥ 0 := Nat.zero_le _
  rw [finishEdge_eq]
  simp only [finishEdge', wp_bind, wp_get, wp_stackDir, hge, ↓reduceIte]
  unfold finishBoundary
  simp only [wp_bind, wp_modifyItem, wp_modify, wp_ite, wp_allocItem, wp_makeVs, wp_popTstack,
    wp_pure]
  split
  · rename_i hT
    split
    · rename_i hL
      obtain ⟨t, hsub, -, -⟩ := hE.bd_bridge hT (by simpa using hL)
      simp [hts, hsub]
    · rename_i hL
      obtain ⟨t₁, t₂, hsub, -, -, -, -⟩ := hE.bd_comp hT hge (by simpa using hL)
      simp [hts, hsub]
  · rename_i hT
    simp [hback (Bool.eq_false_iff.2 hT)]

theorem rootOuts : ∀ (v : Nat) (outs : List DfsOut) (s : WalkState),
    BookOuts v 0 outs false s → s.tstack = [] →
    wp (walkOuts v 0 outs false) (fun hv s' => hv = false ∧ s'.tstack = []) s
  | v, [], s, _, hts => by unfold walkOuts; simp [wp_pure, hts]
  | v, o :: rest, s, hb, hts => by
    unfold BookOuts at hb
    have hout : wp (walkOut v 0 o false) (fun hv s' => hv = false ∧ s'.tstack = []) s := by
      have hb₀ := hb.1
      unfold BookOut at hb₀
      have hb₁ := hb₀.2
      rw [walkOut_eq, wp_bind]
      have hpre : ∀ Q : Bool → WalkState → Prop,
          wp (walkOutPre v 0 o false) Q s = Q false { s with stackDir := s.stackDir.set! 0 false } := by
        intro Q
        unfold walkOutPre
        (try simp only [wp_bind, wp_stackDir, wp_setStackDir, wp_ite, wp_pure])
        (try rfl)
      rw [hpre] at hb₁ ⊢
      unfold walkOutRest
      rw [wp_bind, wp_tstackSize]
      simp only [hts, List.length_nil] at hb₁ ⊢
      cases o with
      | back e cls dest =>
        exact finishEdge_root hb₁ fun _ => by simp
      | tree e cls child =>
        simp only [wp_bind, wp_modify] at hb₁ ⊢
        exact wp_imp (wp_of_forall fun _ s₃ hb₃ => finishEdge_root hb₃ fun h => by
          rw [hb₃.tree.2 ⟨_, _, _, rfl⟩] at h; cases h) hb₁.2
    unfold walkOuts
    simp only [wp_bind]
    refine wp_imp (wp_imp (wp_of_forall fun hv' s' (h : hv' = false ∧ s'.tstack = [])
      (hb' : BookOuts v 0 rest hv' s') => ?_) hout) hb.2
    obtain ⟨rfl, hts'⟩ := h
    exact rootOuts v rest s' hb' hts'

/-- After `walkTree t 0` from an empty tstack the stack is the root's vertex entry, on side 2
(`stackDir[0] = true`). -/
theorem walkTree_rootOK (t : DfsTree) (s : WalkState)
    (hi : ∀ v outs, t = .node v outs →
      ({ s with stackVerts := s.stackVerts.set! 0 v } : WalkState).Inv' 0)
    (hs : Shape s) (hg : GuardsTree t 0 s) (hb : BookTree t 0 s)
    (hts : s.tstack = []) (hsd : 0 < s.stackDir.size) :
    RootOK ((walkTree t 0).run s).2 := by
  obtain ⟨v, outs⟩ := t
  show wp (walkTree (.node v outs) 0) (fun _ s' => RootOK s') s
  unfold walkTree
  simp only [wp_bind, wp_modify]
  unfold GuardsTree at hg; unfold BookTree at hb
  have h₁ := invOuts v 0 outs false _ (hi v outs rfl) hs.frame' hg hb
  have h₂ := rootOuts v outs _ hb hts
  have h₃ := frOutsP v 0 outs false { s with stackVerts := s.stackVerts.set! 0 v }
  refine wp_mono _ (wp_and h₁ (wp_and h₂ h₃)) fun hv s' h => ?_
  obtain ⟨⟨-, -, hvb⟩, ⟨hv0, hts'⟩, -, hsz⟩ := h
  subst hv0
  obtain ⟨hvlt, -, -⟩ := hvb rfl
  simp only [Bool.false_eq_true, ↓reduceIte, wp_bind, wp_setStackDir, wp_pushVertTstack]
  have hsz' : 0 < s'.stackDir.size := by rw [hsz]; exact hsd
  have hdir : (s'.stackDir.setIfInBounds 0 true)[0]! = true :=
    Array.getElem!_set!_self _ _ _ hsz'
  refine ⟨⟨v, 0, s'.nxtEdgeIdx, setSides (s'.stackDir.set! 0 true)[0]! [vertItem v] []⟩, ?_, ?_,
    fun c hc => ⟨v, hvlt, ?_⟩⟩
  · simp [hts']
  · simp [hsz', setSides]
  · simpa [hdir, hsz', setSides] using hc

/-! ### Forest level -/

/-- The per-root inputs of `walk_sides`: `GuardsTree`/`BookTree` (the admitted `walkTree_guards`/
`walkTree_book`) and `Inv' 0` at the start of every root, threaded through the root pop/append of
`walkForest`. -/
def RootsBook : List DfsTree → WalkState → Prop
  | [], _ => True
  | t :: rest, s => GuardsTree t 0 s ∧ BookTree t 0 s ∧ s.Inv' 0 ∧
      wp (walkTree t 0) (fun _ s₁ =>
        wp (popTstack >>= fun top => modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 })
          (fun _ s₂ => RootsBook rest s₂) s₁) s

theorem Inv'.stackVerts_of_nil {s : WalkState} {D : Nat} (h : s.Inv' D) (hts : s.tstack = []) (sv : Array Nat) :
    ({ s with stackVerts := sv } : WalkState).Inv' D :=
  ⟨fun _ _ _ h' => by simp [hts] at h',
   fun i h1 h2 => ⟨(h.nodes i h1 h2).conn, (h.nodes i h1 h2).attached⟩⟩

/-- `SidesForest` from the per-root inputs: `SidesTree` by `walkTree_sides`, `RootOK` by
`walkTree_rootOK`, and `Shape`/empty stack/`stackDir` size carried through the root pop/append. -/
theorem sidesForest_of_roots : ∀ (forest : List DfsTree) (s : WalkState),
    RootsBook forest s → Shape s → s.tstack = [] → 0 < s.stackDir.size → SidesForest forest s
  | [], _, _, _, _, _ => trivial
  | t :: rest, s, h, hs, hts, hsd => by
    obtain ⟨hg, hb, hi, hrb⟩ := h
    have hi' : ∀ v outs, t = .node v outs →
        ({ s with stackVerts := s.stackVerts.set! 0 v } : WalkState).Inv' 0 :=
      fun _ _ _ => hi.stackVerts_of_nil hts _
    refine ⟨walkTree_sides t 0 s hi' hs hg hb, ?_⟩
    have hrk : wp (walkTree t 0) (fun _ s' => RootOK s') s := walkTree_rootOK t s hi' hs hg hb hts hsd
    have hinv := invTree t 0 s hi' hs hg hb
    have hfr : wp (walkTree t 0) (fun _ s' => Fr s s') s := frTreeP t 0 s
    refine wp_mono _ (wp_and hrk (wp_and hinv (wp_and hfr hrb)))
      fun _ s₁ ⟨hrk₁, ⟨_, hs₁⟩, hfr₁, hrb₁⟩ => ⟨hrk₁, ?_⟩
    obtain ⟨tt, hts₁, -, -⟩ := hrk₁
    simp only [wp_bind, wp_popTstack, wp_modifyItem] at hrb₁ ⊢
    refine sidesForest_of_roots rest _ hrb₁ ?_ (by simp [hts₁]) ?_
    · refine (hs₁.tstack (l := s₁.tstack.tail) fun e he => hs₁.span e (List.mem_of_mem_tail he)).modify
        rootItem _ (fun _ => rfl) fun hj c hc => ?_
      simp only [List.mem_append] at hc
      rcases hc with hc | hc
      · exact hs₁.ch_lt rootItem c (by simpa [Items.ch, Items.IsParent, hj] using hc)
      · exact hs₁.span tt (by simp [hts₁]) c (List.mem_append_right _ (by simpa [hts₁] using hc))
    · show 0 < s₁.stackDir.size
      rw [hfr₁.2]; exact hsd

theorem init_shape (g : Graph) (ternarize : Bool) : Shape (WalkState.init g ternarize) where
  size := by simp [WalkState.init, Items.initialItems_size]
  root := by simp [WalkState.init, Items.initialItems_type, rootItem]
  vert v hv := by simp [WalkState.init, Items.initialItems_type, vertItem] at hv ⊢; omega
  edge e he := by simp [WalkState.init, Items.initialItems_type, edgeItem] at he ⊢; omega
  ch_lt p c hp := by simp [WalkState.init, Items.IsParent, Items.initialItems_ch] at hp
  span t ht := by simp [WalkState.init] at ht

theorem init_inv (g : Graph) (ternarize : Bool) : (WalkState.init g ternarize).Inv' 0 :=
  ⟨fun _ _ _ h => by simp [WalkState.init] at h,
   fun i h1 h2 => by simp [WalkState.init, Items.initialItems_size] at h1 h2; omega⟩

/-- `walk_sides` reduced to the per-root inputs. -/
theorem walk_sides_of_roots (g : Graph) (ternarize : Bool) (forest : List DfsTree) (hf : ForestOK g forest)
    (hrb : RootsBook forest (WalkState.init g ternarize)) :
    SidesForest forest (WalkState.init g ternarize) := by
  cases forest with
  | nil => trivial
  | cons t rest =>
    obtain ⟨v, outs⟩ := t
    have hv : v < g.nv := hf.verts_lt v (by simp [DfsTree.verts])
    exact sidesForest_of_roots _ _ hrb (init_shape g ternarize) rfl (by simp [WalkState.init]; omega)

end WalkState

end Spqr
