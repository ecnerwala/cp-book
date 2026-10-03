import Spqr.WalkInv
import Spqr.WalkPlace

/-!
# The item frame of a subtree walk

`walkTree` writes the child list `ch` of a fixed item (`0 < i < 1 + nv + ne`) only for the `V` item
of the current vertex and the `Q` item of the edge being finished: every other `ch` write goes to a
node allocated by the walk (`allocItem`) or reopened by `maybeUnwrapNxt`, whose type is `S`/`P`/`R`
and hence (the fixed items being typed `F`/`V`/`Q`) is not a fixed item. `Keep j s s'` is the frame
of one step for the fixed item `j` (graph, array sizes, item types, `ch j`), `Types g s` the typing
of the fixed items it needs; `walkTree_frame` is the per-subtree statement `EarWalk.lean` uses.
-/

namespace Spqr
open WalkM

namespace WalkState

/-- The fixed items are typed `F`/`V`/`Q` and allocated. -/
structure Types (g : Graph) (s : WalkState) : Prop where
  g_eq : s.g = g
  size : 1 + g.nv + g.ne ≤ s.items.size
  root : Items.type s.items rootItem = .F
  vert : ∀ v, v < g.nv → Items.type s.items (vertItem v) = .V
  edge : ∀ e, e < g.ne → Items.type s.items (edgeItem g e) = .Q

theorem Place.types {g : Graph} {P X : ItemId → Prop} {s : WalkState} (h : s.Place g P X) : Types g s :=
  ⟨h.g_eq, h.size, h.root_type, h.vert, h.edge⟩

/-- One walk step keeps the graph, the array sizes, the types of the existing items and `ch j`. -/
structure Keep (j : Nat) (s s' : WalkState) : Prop where
  g : s'.g = s.g
  sv : s'.stackVerts.size = s.stackVerts.size
  sd : s'.stackDir.size = s.stackDir.size
  fo : s'.firstOccurrence.size = s.firstOccurrence.size
  size : s.items.size ≤ s'.items.size
  type : ∀ k, k < s.items.size → Items.type s'.items k = Items.type s.items k
  ch : Items.ch s'.items j = Items.ch s.items j

variable {j : Nat} {s : WalkState}

theorem Keep.refl : Keep j s s := ⟨rfl, rfl, rfl, rfl, Nat.le_refl _, fun _ _ => rfl, rfl⟩
theorem Keep.trans {s₁ s₂ : WalkState} (h₁ : Keep j s s₁) (h₂ : Keep j s₁ s₂) : Keep j s s₂ :=
  ⟨h₂.g.trans h₁.g, h₂.sv.trans h₁.sv, h₂.sd.trans h₁.sd, h₂.fo.trans h₁.fo, h₁.size.trans h₂.size,
    fun k hk => (h₂.type k (Nat.lt_of_lt_of_le hk h₁.size)).trans (h₁.type k hk), h₂.ch.trans h₁.ch⟩

theorem Types.of_keep {g : Graph} (hT : Types g s) {s' : WalkState} (h : Keep j s s') : Types g s' where
  g_eq := h.g.trans hT.g_eq
  size := hT.size.trans h.size
  root := by rw [h.type _ (by have := hT.size; show 0 < _; omega)]; exact hT.root
  vert v hv := by rw [h.type _ (by have := hT.size; show 1 + v < _; omega)]; exact hT.vert v hv
  edge e he := by rw [h.type _ (by have := hT.size; show 1 + g.nv + e < _; omega)]; exact hT.edge e he

/-- A fixed item is never a node of type `S`/`P`/`R`. -/
theorem Types.ne_of_type {g : Graph} (hT : Types g s) (hj : j < 1 + g.nv + g.ne)
    {x : Nat} (hx : Items.type s.items x ≠ .F ∧ Items.type s.items x ≠ .V ∧ Items.type s.items x ≠ .Q) :
    x ≠ j := by
  rintro rfl
  have h1 : (x : Nat) < 1 + g.nv + g.ne := hj
  have h0 : 0 < x := by
    rcases Nat.eq_zero_or_pos x with h | h
    · subst h; exact absurd hT.root hx.1
    · exact h
  rcases Nat.lt_or_ge (x : Nat) (1 + g.nv) with h | h
  · have := hT.vert (x - 1) (by omega)
    rw [show vertItem (x - 1) = x from Nat.add_sub_of_le h0] at this
    exact hx.2.1 this
  · have := hT.edge (x - (1 + g.nv)) (Nat.sub_lt_left_of_lt_add h h1)
    rw [show edgeItem g (x - (1 + g.nv)) = x from Nat.add_sub_of_le h] at this
    exact hx.2.2 this

theorem Keep.tstack (ts : List TEntry) : Keep j s { s with tstack := ts } :=
  ⟨rfl, rfl, rfl, rfl, Nat.le_refl _, fun _ _ => rfl, rfl⟩
theorem Keep.setStackDir (d : Nat) (b : Bool) : Keep j s { s with stackDir := s.stackDir.set! d b } :=
  ⟨rfl, rfl, by simp, rfl, Nat.le_refl _, fun _ _ => rfl, rfl⟩
theorem Keep.modifyItem_ne {i : Nat} (f : Item → Item) (hf : ∀ it, (f it).type = it.type) (hij : i ≠ j) :
    Keep j s { s with items := s.items.modify i f } :=
  ⟨rfl, rfl, rfl, rfl, by simp, fun k _ => Items.type_modify _ _ _ _ hf, Items.ch_modify_ne _ _ _ _ hij⟩
theorem Keep.modifyItem_vs (i : Nat) (vs : Option Nat × Option Nat) :
    Keep j s { s with items := s.items.modify i fun it => { it with vs := vs } } :=
  ⟨rfl, rfl, rfl, rfl, by simp, fun k _ => Items.type_modify _ _ _ _ fun _ => rfl,
    Items.ch_modify_of_ch _ _ _ _ fun _ => rfl⟩
theorem Keep.push (ty : NodeType) : Keep j s { s with items := s.items.push ⟨ty, (none, none), []⟩ } :=
  ⟨rfl, rfl, rfl, rfl, by simp, fun k hk => by rw [Items.type_push]; simp [Nat.ne_of_lt hk],
    Items.ch_push_nil _ rfl _⟩

theorem Keep.frame {s' : WalkState} (hg : s'.g = s.g) (hsv : s'.stackVerts.size = s.stackVerts.size)
    (hsd : s'.stackDir.size = s.stackDir.size) (hfo : s'.firstOccurrence.size = s.firstOccurrence.size)
    (hit : s'.items = s.items) : Keep j s s' :=
  ⟨hg, hsv, hsd, hfo, hit ▸ Nat.le_refl _, fun _ _ => by rw [hit], by rw [hit]⟩

theorem type_eq_of_lt {items : Items} {i : Nat} (h : i < items.size) :
    Items.type items i = items[i]!.type := by
  simp [Items.type, Array.getElem?_eq_getElem h, getElem!_pos items i h]

/-! ### `finishEdge`

Every lemma is relative to a base state `s₀`: from `Keep j s₀ s` it concludes `Keep j s₀ s'` for
the state `s'` after the operation (the `Keep` hypothesis comes last so that the current state is
read off the goal). -/

variable {s₀ : WalkState}

theorem keep_mergeTstackTops (h₀ : Keep j s₀ s) :
    wp mergeTstackTops (fun _ s' => Keep j s₀ s') s := by
  simp only [wp_mergeTstackTops]; exact h₀.trans (Keep.frame rfl rfl rfl rfl rfl)

theorem keep_finishTstackTop {item : ItemId} (hij : item ≠ j) (h₀ : Keep j s₀ s) :
    wp (finishTstackTop item) (fun _ s' => Keep j s₀ s') s := by
  unfold finishTstackTop
  simp only [wp_bind, wp_cur, wp_stackDir, wp_makeVs, wp_modifyItem, wp_modifyCur]
  exact h₀.trans ⟨rfl, rfl, rfl, rfl, by simp, fun k _ => Items.type_modify _ _ _ _ fun _ => rfl,
    Items.ch_modify_ne _ _ _ _ hij⟩

/-- `maybeUnwrapNxt` returns a fresh node or a reopened node of type `ty`; with `ty` a node type, it
is not the fixed item `j`. -/
theorem keep_maybeUnwrapNxt {g : Graph} (hT : Types g s₀) (hj : j < 1 + g.nv + g.ne) (ty : NodeType)
    (hty : ty ≠ .F ∧ ty ≠ .V ∧ ty ≠ .Q) (h₀ : Keep j s₀ s) :
    wp (maybeUnwrapNxt ty) (fun item s' => Keep j s₀ s' ∧ item ≠ j) s := by
  have hT := hT.of_keep h₀
  have hsz : s.items.size ≠ j := by have := hT.size; omega
  unfold maybeUnwrapNxt
  simp only [wp_bind, wp_ite, wp_get, wp_allocItem, wp_pure, wp_nxt, wp_stackDir, wp_getItem,
    wp_modifyNxt]
  split
  · exact ⟨h₀.trans (Keep.push _), hsz⟩
  · split
    · rename_i heq
      refine ⟨h₀.trans (Keep.tstack _), ?_⟩
      set x := (getSide s.tstack.tail.head!.spans s.stackDir[s.tstack.tail.head!.topDepth]!).head!
      rcases Nat.lt_or_ge x s.items.size with hx | hx
      · refine hT.ne_of_type hj ?_
        rw [type_eq_of_lt hx, beq_iff_eq.mp heq]; exact hty
      · intro h; rw [h] at hx; exact absurd (Nat.lt_of_lt_of_le hj hT.size) (Nat.not_lt.mpr hx)
    · exact ⟨h₀.trans (Keep.push _), hsz⟩

theorem keep_loop1Body {g : Graph} (hT : Types g s₀) (hj : j < 1 + g.nv + g.ne) (d : Nat)
    (edgeDir : Bool) (h₀ : Keep j s₀ s) : wp (loop1Body d edgeDir) (fun _ s' => Keep j s₀ s') s := by
  unfold loop1Body loop1Type
  simp only [wp_bind, wp_ite, wp_nxt, wp_cur, wp_setStackDir, wp_mergeTstackTops, wp_pure]
  split
  · refine wp_mono _ (keep_maybeUnwrapNxt hT hj .S (by decide)
      (h₀.trans (Keep.frame rfl rfl (by simp) rfl rfl))) fun item s₁ ⟨h₁, hne⟩ => ?_
    exact keep_finishTstackTop hne (h₁.trans (Keep.frame rfl rfl rfl rfl rfl))
  · split <;>
    · refine wp_mono _ (keep_maybeUnwrapNxt hT hj _ (by decide) h₀) fun item s₁ ⟨h₁, hne⟩ => ?_
      exact keep_finishTstackTop hne (h₁.trans (Keep.frame rfl rfl rfl rfl rfl))

theorem keep_loop (n : Nat) (cond : WalkM Bool) (body : WalkM Unit) (hc : ∀ s, (cond.run s).2 = s)
    (hb : ∀ s, Keep j s₀ s → wp body (fun _ s' => Keep j s₀ s') s) (h₀ : Keep j s₀ s) :
    wp (loop n cond body) (fun _ s' => Keep j s₀ s') s :=
  wp_loop (fun s' => Keep j s₀ s') n cond body
    (fun s₁ h₁ => by show Keep j s₀ (cond.run s₁).2; rw [hc]; exact h₁) hb h₀ fun _ h => h

theorem keep_finishRest {g : Graph} (hT : Types g s₀) (hj : j < 1 + g.nv + g.ne)
    (curV d lowval : Nat) (isType1 hasVert isSingle : Bool) (h₀ : Keep j s₀ s) :
    wp (finishRest curV d lowval isType1 hasVert isSingle) (fun _ s' => Keep j s₀ s') s := by
  unfold finishRest finishP finishTail condP
  simp only [wp_bind, wp_ite, wp_tstackSize, wp_nxt, wp_pure, wp_pushVertTstack, wp_mergeTstackTops]
  split
  · refine wp_mono _ (keep_maybeUnwrapNxt hT hj .P (by decide) h₀) fun item s₁ ⟨h₁, hne⟩ => ?_
    refine wp_mono _ (keep_finishTstackTop hne (h₁.trans (Keep.frame rfl rfl rfl rfl rfl)))
      fun _ s₂ h₂ => ?_
    try simp only [wp_bind, wp_ite, wp_pure, wp_pushVertTstack, wp_mergeTstackTops]
    split <;> (try split) <;> exact h₂.trans (Keep.frame rfl rfl rfl rfl rfl)
  · split <;> (try split) <;> exact h₀.trans (Keep.frame rfl rfl rfl rfl rfl)

theorem keep_closeVertTail (curV : Nat) (edgeDir isSingle : Bool) (k : Bool → WalkM Bool)
    (item : Option ItemId) (hitem : ∀ i, item = some i → i ≠ j)
    (hk : ∀ b s, Keep j s₀ s → wp (k b) (fun _ s' => Keep j s₀ s') s) (h₀ : Keep j s₀ s) :
    wp (closeVertTail curV edgeDir isSingle k item) (fun _ s' => Keep j s₀ s') s := by
  unfold closeVertTail
  cases item with
  | none =>
    simp only [wp_bind, wp_mergeTstackTops, wp_modifyCur]
    exact hk _ _ (h₀.trans (Keep.frame rfl rfl rfl rfl rfl))
  | some item =>
    simp only [wp_bind, wp_mergeTstackTops, wp_modifyCur]
    exact wp_mono _ (keep_finishTstackTop (hitem item rfl) (h₀.trans (Keep.frame rfl rfl rfl rfl rfl)))
      fun _ s₂ h₂ => hk _ _ h₂

theorem keep_closeVert {g : Graph} (hT : Types g s₀) (hj : j < 1 + g.nv + g.ne) (curV : Nat)
    (edgeDir isType1 : Bool) (origTstack : Nat) (isSingle : Bool) (k : Bool → WalkM Bool)
    (hk : ∀ b s, Keep j s₀ s → wp (k b) (fun _ s' => Keep j s₀ s') s) (h₀ : Keep j s₀ s) :
    wp (closeVert curV edgeDir isType1 origTstack isSingle k) (fun _ s' => Keep j s₀ s') s := by
  unfold closeVert
  simp only [wp_bind, wp_ite, wp_tstackSize, wp_map]
  split
  · refine wp_mono _ (keep_loop _ _ _ (fun s => by rw [run_loop3Cond])
      (fun s h => keep_mergeTstackTops h) h₀) fun _ s₃ h₃ => ?_
    exact keep_closeVertTail _ _ _ _ _ (fun _ h => nomatch h) hk h₃
  · refine wp_mono _ (keep_maybeUnwrapNxt hT hj _ (by split <;> decide) h₀) fun item s₃ ⟨h₃, hne⟩ => ?_
    exact keep_closeVertTail _ _ _ _ _ (fun i h => Option.some.inj h ▸ hne) hk h₃

theorem keep_finishTree {g : Graph} (hT : Types g s₀) (hj : j < 1 + g.nv + g.ne) (curV d : Nat)
    (o : DfsOut) (origTstack : Nat) (hasVert edgeDir : Bool) (h₀ : Keep j s₀ s) :
    wp (finishTree curV d o origTstack hasVert edgeDir) (fun _ s' => Keep j s₀ s') s := by
  unfold finishTree closeEars mergeLate
  simp only [wp_bind, wp_pushEdgeTstack, wp_tstackSize, wp_get, wp_cur, wp_ite, wp_pure]
  refine wp_mono _ (keep_loop _ _ _ (fun s => by rw [run_loop1Cond])
    (fun s h => keep_loop1Body hT hj d _ h) (h₀.trans (Keep.frame rfl rfl rfl rfl rfl)))
    fun _ s₁ h₁ => ?_
  have fin : ∀ (b : Bool) (s₂ : WalkState), Keep j s₀ s₂ →
      if hasVert = true then
        wp (closeVert curV edgeDir o.cls.isType1 origTstack b
          (finishRest curV d (o.cls.lowval d) o.cls.isType1 hasVert)) (fun _ s' => Keep j s₀ s') s₂
      else wp (finishRest curV d (o.cls.lowval d) o.cls.isType1 hasVert b) (fun _ s' => Keep j s₀ s') s₂ := by
    intro b s₂ h₂
    split
    · exact keep_closeVert hT hj _ _ _ _ _ _ (fun b s h => keep_finishRest hT hj _ _ _ _ _ _ h) h₂
    · exact keep_finishRest hT hj _ _ _ _ _ _ h₂
  split
  · exact wp_mono _ (keep_loop _ _ _ (fun s => by rw [run_loop2Cond])
      (fun s h => keep_mergeTstackTops h) h₁) fun _ s₂ h₂ => fin false s₂ h₂
  · exact fin true s₁ h₁

theorem keep_finishBoundary {g : Graph} (hT : Types g s₀) (hj : j < 1 + g.nv + g.ne) (curV d : Nat)
    (o : DfsOut) (qItem : ItemId) (hasVert : Bool) (hq : qItem ≠ j) (hv : vertItem curV ≠ j)
    (h₀ : Keep j s₀ s) :
    wp (finishBoundary curV d o qItem hasVert) (fun _ s' => Keep j s₀ s') s := by
  have hsz : s.items.size ≠ j := by have := (hT.of_keep h₀).size; omega
  unfold finishBoundary
  simp only [wp_bind, wp_modifyItem, wp_modify, wp_ite, wp_allocItem, wp_makeVs, wp_popTstack,
    wp_pure]
  split <;> (try split) <;>
  · refine h₀.trans ⟨rfl, rfl, rfl, rfl, by simp, fun k hk => ?_, ?_⟩
    · simp [Items.type_modify, Items.type_push, Nat.ne_of_lt hk]
    · simp [Items.ch_modify_ne, Items.ch_push_nil, hv, hq, hsz]

theorem keep_finishEdge {g : Graph} (hT : Types g s₀) (hj : j < 1 + g.nv + g.ne) (curV d : Nat)
    (o : DfsOut) (origTstack : Nat) (hasVert : Bool) (hv : vertItem curV ≠ j)
    (he : edgeItem g o.e ≠ j) (h₀ : Keep j s₀ s) :
    wp (finishEdge curV d o origTstack hasVert) (fun _ s' => Keep j s₀ s') s := by
  have he' : edgeItem s.g o.e ≠ j := by rw [(hT.of_keep h₀).g_eq]; exact he
  by_cases hge : o.cls.lowval d ≥ d
  · simp only [finishEdge_eq, finishEdge', wp_bind, wp_get, wp_stackDir, hge, ↓reduceIte]
    exact keep_finishBoundary hT hj _ _ _ _ _ he' hv h₀
  · by_cases ht : o.cls.isTree = true
    · simp only [finishEdge_eq, finishEdge', wp_bind, wp_get, wp_stackDir, hge, ht, ↓reduceIte,
        wp_makeVs, wp_modifyItem]
      exact keep_finishTree hT hj _ _ _ _ _ _ (h₀.trans ⟨rfl, rfl, rfl, rfl, by simp,
        fun k _ => Items.type_modify _ _ _ _ fun _ => rfl, Items.ch_modify_of_ch _ _ _ _ fun _ => rfl⟩)
    · have ht' : o.cls.isTree = false := Bool.eq_false_iff.2 ht
      simp only [finishEdge_eq, finishEdge', finishBack, wp_bind, wp_get, wp_stackDir, hge, ht',
        Bool.false_eq_true, ↓reduceIte, wp_makeVs, wp_modifyItem, wp_pushEdgeTstack, wp_modify]
      exact keep_finishRest hT hj _ _ _ _ _ _ (h₀.trans ⟨rfl, rfl, rfl, by simp, by simp,
        fun k _ => Items.type_modify _ _ _ _ fun _ => rfl, Items.ch_modify_of_ch _ _ _ _ fun _ => rfl⟩)

/-! ### The walk -/

abbrev KTree (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  ∀ (g : Graph) (j : Nat) (s₀ : WalkState), Types g s₀ → j < 1 + g.nv + g.ne →
    (∀ v ∈ t.verts, vertItem v ≠ j) → (∀ e ∈ t.edges, edgeItem g e ≠ j) → Keep j s₀ s →
    wp (walkTree t d) (fun _ s' => Keep j s₀ s') s
abbrev KOuts (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  ∀ (g : Graph) (j : Nat) (s₀ : WalkState), Types g s₀ → j < 1 + g.nv + g.ne → vertItem v ≠ j →
    (∀ w ∈ DfsOut.vertsList outs, vertItem w ≠ j) → (∀ e ∈ DfsOut.edgesList outs, edgeItem g e ≠ j) →
    Keep j s₀ s → wp (walkOuts v d outs hasVert) (fun _ s' => Keep j s₀ s') s
abbrev KOut (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  ∀ (g : Graph) (j : Nat) (s₀ : WalkState), Types g s₀ → j < 1 + g.nv + g.ne → vertItem v ≠ j →
    (∀ w ∈ o.verts, vertItem w ≠ j) → (∀ e ∈ o.edges, edgeItem g e ≠ j) → Keep j s₀ s →
    wp (walkOut v d o hasVert) (fun _ s' => Keep j s₀ s') s

mutual
theorem kTree : ∀ (t : DfsTree) (d : Nat) (s : WalkState), KTree t d s
  | .node v outs, d, s => by
    intro g j s₀ hT hj hv he h₀
    simp only [DfsTree.verts, DfsTree.edges, List.mem_cons, forall_eq_or_imp] at hv he
    unfold walkTree
    simp only [wp_bind, wp_modify]
    refine wp_mono _ (kOuts v d outs false _ g j s₀ hT hj hv.1 hv.2 he
      (h₀.trans (Keep.frame rfl (by simp) rfl rfl rfl))) fun hv' s' h' => ?_
    cases hv'
    · simp only [Bool.false_eq_true, ↓reduceIte, wp_bind, wp_setStackDir, wp_pushVertTstack]
      exact h'.trans (Keep.frame rfl rfl (by simp) rfl rfl)
    · exact h'

theorem kOuts : ∀ (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState),
    KOuts v d outs hasVert s
  | v, d, [], hasVert, s => by
    intro g j s₀ _ _ _ _ _ h₀
    unfold walkOuts; simp only [wp_pure]; exact h₀
  | v, d, o :: rest, hasVert, s => by
    intro g j s₀ hT hj hv hw he h₀
    rw [DfsOut.vertsList_eq, List.flatMap_cons] at hw
    rw [DfsOut.edgesList_eq, List.flatMap_cons] at he
    simp only [List.mem_append, ← DfsOut.vertsList_eq, ← DfsOut.edgesList_eq] at hw he
    unfold walkOuts
    simp only [wp_bind]
    refine wp_mono _ (kOut v d o hasVert s g j s₀ hT hj hv (fun w h => hw w (.inl h))
      (fun e h => he e (.inl h)) h₀) fun hv' s' h' => ?_
    exact kOuts v d rest hv' s' g j s₀ hT hj hv (fun w h => hw w (.inr h)) (fun e h => he e (.inr h)) h'

theorem kOut : ∀ (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState), KOut v d o hasVert s
  | v, d, o, hasVert, s => by
    intro g j s₀ hT hj hv hw he h₀
    rw [walkOut_eq, wp_bind]
    have pre : wp (walkOutPre v d o hasVert) (fun _ s' => Keep j s₀ s') s := by
      unfold walkOutPre
      simp only [wp_bind, wp_stackDir, wp_setStackDir, wp_ite, wp_pushVertTstack, wp_pure]
      split <;> exact h₀.trans (Keep.frame rfl rfl (by simp) rfl rfl)
    refine wp_mono _ pre fun hv' s₁ h₁ => ?_
    unfold walkOutRest
    rw [wp_bind, wp_tstackSize]
    cases o with
    | back e dest cls =>
      exact keep_finishEdge hT hj _ _ _ _ _ hv (he e (by simp [DfsOut.edges])) h₁
    | tree e cls child =>
      simp only [wp_bind, wp_modify]
      refine wp_mono _ (kTree child (d + 1) _ g j s₀ hT hj
        (fun w h => hw w (by simpa [DfsOut.verts] using h))
        (fun e h => he e (by simp [DfsOut.edges, h]))
        (h₁.trans (Keep.frame rfl rfl rfl (by simp) rfl))) fun _ s₃ h₃ => ?_
      exact keep_finishEdge hT hj _ _ _ _ _ hv (he e (by simp [DfsOut.edges])) h₃
end

/-- The frame of a subtree walk: the graph and the array sizes are untouched, and the child list of
every fixed item other than the `V` items of the tree's vertices and the `Q` items of its edges is
unchanged (the walk only writes `ch` of those and of the `S`/`P`/`R` nodes it allocates or reopens;
the typing `Types g s` is what rules out reopening a fixed item). -/
theorem walkTree_frame (t : DfsTree) (d : Nat) (s : WalkState) {g : Graph} (hT : Types g s) :
    wp (walkTree t d) (fun _ s' => s'.g = s.g ∧ s'.stackVerts.size = s.stackVerts.size ∧
      s'.stackDir.size = s.stackDir.size ∧ s'.firstOccurrence.size = s.firstOccurrence.size ∧
      ∀ i, 0 < i → i < 1 + s.g.nv + s.g.ne → (∀ v ∈ t.verts, i ≠ vertItem v) →
        (∀ e ∈ t.edges, i ≠ edgeItem s.g e) → Items.ch s'.items i = Items.ch s.items i) s := by
  have hg := hT.g_eq
  have h0 := kTree t d s g 0 s hT (by omega) (fun v _ => by simp [vertItem])
    (fun e _ => by simp [edgeItem]) Keep.refl
  refine ⟨h0.g, h0.sv, h0.sd, h0.fo, fun i _ hi hv he => ?_⟩
  rw [hg] at hi he
  exact (kTree t d s g i s hT hi (fun v h => (hv v h).symm) (fun e h => (he e h).symm) Keep.refl).ch

end WalkState
end Spqr
