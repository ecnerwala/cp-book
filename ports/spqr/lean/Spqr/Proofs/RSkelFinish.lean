import Spqr.Proofs.RSkelInv

/-!
# `Items.RSkelInv` across `finishEdge` (PROOF.md §4.5)

`KeepsR m s`: running `m` from `s` keeps `Items.RSkelInv`. Every block of `finishEdge` keeps it by
frame reasoning (childless pushes, parentless modifications) except the two R closes — the `.R`
iterate of Loop 1 and the type-1 vertex close with `isSingle = false` — whose content is taken as a
hypothesis (`hR1`/`hR2` of `keepsR_finishEdge`) and supplied by the R-branch theory.
-/

namespace Spqr
open WalkState WalkM

def KeepsR {α : Type} (m : WalkM α) (s : WalkState) : Prop :=
  Items.RSkelInv s.g s.items → Items.RSkelInv s.g (m.run s).2.items

theorem two_entries {s : WalkState} (h : 2 ≤ s.tstack.length) :
    ∃ a b rest, s.tstack = a :: b :: rest := by
  match hts : s.tstack with
  | [] | [_] => rw [hts] at h; simp at h
  | a :: b :: rest => exact ⟨a, b, rest, rfl⟩

theorem one_entry {s : WalkState} (h : s.tstack ≠ []) : ∃ t rest, s.tstack = t :: rest := by
  match hts : s.tstack with
  | [] => exact absurd hts h
  | t :: rest => exact ⟨t, rest, rfl⟩

/-- Pushing a fresh `.R` item and filling it: the old items are untouched, the new one is given. -/
theorem Items.RSkelInv.push_modify_R {g : Graph} {items : Items} (h : Items.RSkelInv g items)
    (hc : ∀ i, i < items.size → ∀ c ∈ items.ch i, c < items.size)
    (hroot : ∀ p, ¬ items.IsParent p items.size) (f : Item → Item)
    (hnew : Items.RSkel3 g ((items.push ⟨.R, (none, none), []⟩).modify items.size f) items.size) :
    Items.RSkelInv g ((items.push ⟨.R, (none, none), []⟩).modify items.size f) := by
  intro i hi hty
  rw [Array.size_modify, Array.size_push] at hi
  rcases Nat.lt_or_eq_of_le (Nat.le_of_lt_succ hi) with hlt | rfl
  · have hne := Nat.ne_of_lt hlt
    rw [Items.type_modify_of_ne _ _ hne, Items.type_push_of_ne _ hne] at hty
    refine ((h i hlt hty).push_nil hlt (hc i hlt) _ rfl).modify_of_not_below ?_ _
    rw [Items.Below_push_nil _ rfl]
    exact Items.not_below_of_ne_root hne hroot
  · exact hnew

theorem keepsR_of_items_eq {α : Type} {m : WalkM α} {s : WalkState} (hit : (m.run s).2.items = s.items) :
    KeepsR m s := fun h => by rw [hit]; exact h

theorem keepsR_mergeTop (s : WalkState) : KeepsR mergeTstackTops s :=
  keepsR_of_items_eq (by rw [run_mergeTstackTops])

theorem loop_items_eq (cond : WalkM Bool) (body : WalkM Unit) (fuel : Nat)
    (hcond : ∀ s, (cond.run s).2 = s) (hbody : ∀ s, (body.run s).2.items = s.items) (s : WalkState) :
    ((loop fuel cond body).run s).2.items = s.items := by
  induction fuel generalizing s with
  | zero => rfl
  | succ fuel ih =>
    rw [loop_succ_run fuel cond body s (hcond s)]
    split
    · rw [ih, hbody]
    · rfl

theorem mergeLoop_items_eq (cond : WalkM Bool) (fuel : Nat) (hcond : ∀ s, (cond.run s).2 = s)
    (s : WalkState) : ((loop fuel cond mergeTstackTops).run s).2.items = s.items :=
  loop_items_eq cond _ fuel hcond (fun s => by rw [run_mergeTstackTops]) s

theorem Items.type_finishTop (items : Items) (item : ItemId) (vs : Option Nat × Option Nat) (ch : List ItemId)
    (p : ItemId) :
    Items.type (items.modify item fun it => { it with vs := vs, ch := ch }) p = items.type p :=
  Items.type_modify_type_eq item (fun it => { it with vs := vs, ch := ch }) (fun _ => rfl) p

theorem keepsR_finishTop {s : WalkState} {t : TEntry} {rest : List TEntry} (hts : s.tstack = t :: rest)
    (item : ItemId) (hroot : ∀ p, ¬ Items.IsParent s.items p item)
    (hnew : Items.type s.items item = .R →
      Items.RSkel3 s.g ((finishTstackTop item).run s).2.items item) :
    KeepsR (finishTstackTop item) s := by
  intro h
  rw [finishTstackTop_run_eq s item t rest hts] at hnew ⊢
  exact h.modify_root hroot _ fun hR => hnew (by rwa [Items.type_finishTop] at hR)

theorem keepsR_maybeUnwrap {s : WalkState} (hs : Shape s) (ty : NodeType) (hty : ty ≠ .R) {a b : TEntry}
    {rest : List TEntry} (hts : s.tstack = a :: b :: rest) : KeepsR (maybeUnwrapNxt ty) s := by
  intro h
  rw [maybeUnwrapNxt_run_eq ty s a b rest hts _ rfl _ rfl]
  split
  · rw [run_allocItem]; exact h.push_nil (fun i _ c hc => hs.ch_lt i c hc) _ rfl hty
  · split
    · exact h
    · rw [run_allocItem]; exact h.push_nil (fun i _ c hc => hs.ch_lt i c hc) _ rfl hty

theorem maybeUnwrapNxt_result_type {s : WalkState} (hs : Shape s) (ty : NodeType) (hok : UnwrapOk ty s)
    {a b : TEntry} {rest : List TEntry} (hts : s.tstack = a :: b :: rest) :
    Items.type ((maybeUnwrapNxt ty).run s).2.items ((maybeUnwrapNxt ty).run s).1 = ty := by
  have hn : nxtE s = b := by rw [nxtE, hts]; rfl
  have hd : nxtDir s = s.stackDir[b.topDepth]! := by rw [nxtDir, hn]
  have hh : nxtHead s = (getSide b.spans s.stackDir[b.topDepth]!).head! := by rw [nxtHead, hn, hd]
  rw [maybeUnwrapNxt_run_eq ty s a b rest hts _ rfl _ rfl]
  split
  · rw [run_allocItem]; exact Items.type_push_size _
  · rename_i h1
    split
    · rename_i h2
      have hu := hok.unwrap h1 (by rw [hh]; exact h2)
      have hmem : nxtHead s ∈ (nxtE s).spans.1 ++ (nxtE s).spans.2 := by
        have := hu.single
        cases hdir : nxtDir s <;> simp [getSide, hdir] at this <;> simp [this]
      have hlt := hs.span _ (by rw [hn, hts]; simp) _ hmem
      dsimp only
      rw [← hh] at h2 ⊢
      rw [getElem!_pos s.items s.nxtHead hlt] at h2
      simp [Items.type, Array.getElem?_eq_getElem hlt, h2]
    · rw [run_allocItem]; exact Items.type_push_size _

theorem keepsR_closeTwo {D : Nat} {s : WalkState} (hi : s.Inv' D) (hs : Shape s) (item : ItemId)
    (hok : CloseTwoOk D s) (hf : ItemFree s item)
    (hnew : Items.type s.items item = .R →
      Items.RSkel3 s.g ((finishTstackTop item).run (mergeTstackTops.run s).2).2.items item) :
    Items.RSkelInv s.g s.items →
      Items.RSkelInv s.g ((finishTstackTop item).run (mergeTstackTops.run s).2).2.items := by
  intro h
  have st := Step.mergeTop (v := 0) hi hs hok.merge
  have h₁ : Items.RSkelInv s.g (mergeTstackTops.run s).2.items := keepsR_mergeTop s h
  obtain ⟨t, rest, hts⟩ := one_entry hok.finish.nonempty
  have hts' : (mergeTstackTops.run s).2.tstack = t :: rest := hts
  have := keepsR_finishTop hts' item hf.merge.root fun hR => by
    rw [st.g]
    exact hnew (by rw [run_mergeTstackTops] at hR; exact hR)
  unfold KeepsR at this
  rw [st.g] at this
  exact this h₁

/-- A loop whose body is a step keeping the invariant under `Ok` keeps the invariant. -/
theorem keepsR_loop {D v : Nat} (cond : WalkM Bool) (body : WalkM Unit) (Ok : WalkState → Prop) (fuel : Nat)
    (hcond : ∀ s, (cond.run s).2 = s)
    (hbody : ∀ s, v < s.g.nv → s.Inv' D → Shape s → (cond.run s).1 = true → Ok s →
      Step D v s (body.run s).2 ∧ KeepsR body s)
    {s : WalkState} (hv : v < s.g.nv) (hi : s.Inv' D) (hs : Shape s)
    (hok : ∀ k, (∀ j, j ≤ k → (cond.run (iter body j s)).1 = true) → Ok (iter body k s)) :
    KeepsR (loop fuel cond body) s := by
  induction fuel generalizing s with
  | zero => exact fun h => h
  | succ fuel ih =>
    intro h
    rw [loop_succ_run fuel cond body s (hcond s)]
    by_cases hc : (cond.run s).1 = true
    · simp only [hc, ↓reduceIte]
      obtain ⟨st, hk⟩ := hbody s hv hi hs hc (hok 0 fun j hj => by rw [Nat.le_zero.1 hj]; exact hc)
      have h' := hk h
      rw [← st.g] at h'
      have := ih (by rw [st.g]; exact hv) st.inv st.shape (fun k hk => hok (k + 1) fun j hj => by
        cases j with
        | zero => exact hc
        | succ j => exact hk j (Nat.le_of_succ_le_succ hj)) h'
      rwa [st.g] at this
    · simp only [hc, Bool.false_eq_true, ↓reduceIte]; exact h

theorem l1S₁_items (d : Nat) (edgeDir : Bool) (s : WalkState) : (l1S₁ d edgeDir s).items = s.items := by
  unfold l1S₁ after; rw [loop1Type_run]
  split
  · rw [run_mergeTstackTops]
  · split <;> rfl

theorem keepsR_loop1Body {D d : Nat} {edgeDir : Bool} {s : WalkState} (hi : s.Inv' D) (hs : Shape s)
    (hok : Loop1BodyOk D d edgeDir s)
    (hR : l1Ty d edgeDir s = .R → KeepsR (loop1Body d edgeDir) s) :
    KeepsR (loop1Body d edgeDir) s := by
  by_cases hty : l1Ty d edgeDir s = .R
  · exact hR hty
  intro h
  rw [loop1Body_run_eq]; unfold closeAt
  have st₁ : Step D 0 s (l1S₁ d edgeDir s) := Step.loop1Type hi hs hok.mergeS
  have h₁ : Items.RSkelInv s.g (l1S₁ d edgeDir s).items := by rw [l1S₁_items]; exact h
  have r := maybeUnwrapNxt_spec (v := 0) st₁.inv st₁.shape (loop1Type_result d edgeDir s) hok.unwrap
  obtain ⟨a, b, rest, hts⟩ := two_entries hok.unwrap.two
  have h₂ : Items.RSkelInv s.g (l1S₂ d edgeDir s).items := by
    have := keepsR_maybeUnwrap st₁.shape _ hty hts
    unfold KeepsR at this
    rw [st₁.g] at this
    exact this h₁
  have hty₂ := maybeUnwrapNxt_result_type st₁.shape _ hok.unwrap hts
  have := keepsR_closeTwo r.step.inv r.step.shape _ hok.close r.free
    fun hR' => absurd (hty₂.symm.trans hR') hty
  rw [r.step.g, st₁.g] at this
  exact this h₂

theorem keepsR_closeEars {D nxtV d e v : Nat} {edgeDir : Bool} {s : WalkState} (hi : s.Inv' D) (hs : Shape s)
    (hv : v < s.g.nv) (hok : CloseEarsOk D nxtV d e edgeDir s)
    (hR : ∀ k, (∀ j, j ≤ k → result (loop1Cond d) (iter (loop1Body d edgeDir) j (ceS₁ nxtV d e s)) = true) →
      l1Ty d edgeDir (iter (loop1Body d edgeDir) k (ceS₁ nxtV d e s)) = .R →
      KeepsR (loop1Body d edgeDir) (iter (loop1Body d edgeDir) k (ceS₁ nxtV d e s))) :
    KeepsR (closeEars nxtV d e edgeDir) s := by
  intro h
  have st₁ : Step D v s (ceS₁ nxtV d e s) := Step.pushEdge hi hs nxtV d e hok.e_lt hok.q hok.ends hok.d_le
  have h₁ : Items.RSkelInv s.g (ceS₁ nxtV d e s).items := h
  have := keepsR_loop (D := D) (v := v) (loop1Cond d) (loop1Body d edgeDir)
    (fun s' => Loop1BodyOk D d edgeDir s' ∧ (l1Ty d edgeDir s' = .R → KeepsR (loop1Body d edgeDir) s'))
    (ceS₁ nxtV d e s).tstack.length (fun _ => rfl)
    (fun s' hv' hi' hs' _ hok' => ⟨Step.loop1Body hi' hs' hv' hok'.1, keepsR_loop1Body hi' hs' hok'.1 hok'.2⟩)
    (by rw [st₁.g]; exact hv) st₁.inv st₁.shape (fun k hk => ⟨hok.body k hk, hR k hk⟩)
  unfold KeepsR at this
  rw [st₁.g] at this
  exact this h₁

theorem mergeLate_items (d : Nat) (s : WalkState) : (after (mergeLate d) s).items = s.items := by
  unfold after; rw [mergeLate_run]
  split
  · exact mergeLoop_items_eq _ _ (fun _ => rfl) s
  · rfl

theorem keepsR_finishP {D : Nat} {s : WalkState} (hi : s.Inv' D) (hs : Shape s) {curV lowval : Nat} {isType1 : Bool} (hok : FinishPOk D curV lowval isType1 s) :
    KeepsR (finishP curV lowval isType1) s := by
  intro h
  simp only [Spqr.finishP, WalkM.run_bind, run_condP]
  by_cases hc : result (condP curV lowval isType1) s = true
  · have h' : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == curV) &&
        (s.tstack.tail.head!.topDepth == lowval)) = true := hc
    simp only [h', ↓reduceIte, WalkM.run_bind]
    obtain ⟨hu, hcl⟩ := hok.ok hc
    obtain ⟨a, b, rest, hts⟩ := two_entries hu.two
    have r := maybeUnwrapNxt_spec (v := 0) hi hs (by decide) hu
    have h₂ := keepsR_maybeUnwrap hs .P (by decide) hts h
    have hty₂ := maybeUnwrapNxt_result_type hs .P hu hts
    have := keepsR_closeTwo r.step.inv r.step.shape _ hcl r.free
      fun hR' => absurd (hty₂.symm.trans hR') (by decide)
    rw [r.step.g] at this
    exact this h₂
  · have h' : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == curV) &&
        (s.tstack.tail.head!.topDepth == lowval)) = false := Bool.eq_false_iff.2 hc
    simp only [h', Bool.false_eq_true, ↓reduceIte, WalkM.pure_run]
    exact h

theorem finishTail_items (curV d : Nat) (hasVert isSingle : Bool) (s : WalkState) :
    ((finishTail curV d hasVert isSingle).run s).2.items = s.items := by
  cases hasVert <;> cases isSingle <;>
    simp only [finishTail, Bool.not_false, Bool.not_true, Bool.false_eq_true, ↓reduceIte, WalkM.run_bind,
      WalkM.pure_run, run_mergeTstackTops] <;> rfl

theorem keepsR_finishRest {D : Nat} {s : WalkState} (hi : s.Inv' D) (hs : Shape s)
    {curV d lowval : Nat} {isType1 hasVert isSingle : Bool}
    (hok : FinishRestOk D curV d lowval isType1 hasVert isSingle s) :
    KeepsR (finishRest curV d lowval isType1 hasVert isSingle) s := by
  intro h
  have h₁ := keepsR_finishP hi hs hok.p h
  show Items.RSkelInv s.g ((finishTail curV d hasVert isSingle).run (after (finishP curV lowval isType1) s)).2.items
  rw [finishTail_items]
  exact h₁

theorem cvS₅_items (curV : Nat) (edgeDir isType1 : Bool) (origTstack : Nat) (isSingle : Bool) (s : WalkState) :
    (cvS₅ curV edgeDir isType1 origTstack isSingle s).items = (cvS₂ isType1 origTstack isSingle s).items := by
  unfold cvS₅ cvS₄ cvS₃ after retarget
  simp only [run_modifyCur, run_mergeTstackTops]

theorem cvS₁_items (isType1 : Bool) (origTstack : Nat) (isSingle : Bool) (s : WalkState) :
    (cvS₁ isType1 origTstack isSingle s).items = s.items := by
  cases isType1
  · unfold cvS₁ after
    simp only [WalkState.vertPre, Bool.not_false, ↓reduceIte, WalkM.run_bind, run_tstackSize, WalkM.pure_run]
    exact mergeLoop_items_eq _ _ (fun _ => rfl) s
  · rfl

/-- `closeVert'` keeps the invariant; the type-1 R close (`isSingle = false`) is taken as `hR`: the
item `(cvS₁ …).items.size` allocated by `maybeUnwrapNxt .R` is `RSkel3` in the output. -/
theorem keepsR_closeVert' {D v : Nat} {s : WalkState} (hi : s.Inv' D) (hs : Shape s) (hv : v < s.g.nv)
    {curV : Nat} {edgeDir isType1 : Bool} {origTstack : Nat} {isSingle : Bool}
    (hok : CloseVertOk D curV edgeDir isType1 origTstack isSingle s)
    (hR : isType1 = true → isSingle = false →
      Items.RSkel3 s.g ((closeVert' curV edgeDir isType1 origTstack isSingle).run s).2.items
        (cvS₁ isType1 origTstack isSingle s).items.size) :
    KeepsR (closeVert' curV edgeDir isType1 origTstack isSingle) s := by
  intro h
  have st₁ : Step D v s (cvS₁ isType1 origTstack isSingle s) := Step.vertPre hi hs hv hok.loop3
  have h₁ : Items.RSkelInv s.g (cvS₁ isType1 origTstack isSingle s).items := by rw [cvS₁_items]; exact h
  have hv₁ : v < (cvS₁ isType1 origTstack isSingle s).g.nv := by rw [st₁.g]; exact hv
  obtain ⟨st₂, hfree⟩ := Step.vertUnwrap (v := v) st₁.inv st₁.shape (isType1 := isType1)
    (isSingle := cvB₁ isType1 origTstack isSingle s) (fun h => by subst h; exact hok.unwrap rfl)
  have st₂ : Step D v (cvS₁ isType1 origTstack isSingle s) (cvS₂ isType1 origTstack isSingle s) := st₂
  have st₃ : Step D v _ (cvS₃ isType1 origTstack isSingle s) := Step.mergeTop st₂.inv st₂.shape hok.merge₁
  have st₄ : Step D v _ (cvS₄ isType1 origTstack isSingle s) := Step.mergeTop st₃.inv st₃.shape hok.merge₂
  have st₅ : Step D v _ (cvS₅ curV edgeDir isType1 origTstack isSingle s) :=
    Step.retarget st₄.inv st₄.shape curV edgeDir hok.retarget
  have hg₅ : (cvS₅ curV edgeDir isType1 origTstack isSingle s).g = s.g := by
    rw [st₅.g, st₄.g, st₃.g, st₂.g, st₁.g]
  cases isType1
  · show Items.RSkelInv s.g (cvS₅ curV edgeDir false origTstack isSingle s).items
    rw [cvS₅_items]; exact h₁
  · obtain ⟨a, b, rest, hts⟩ := two_entries (hok.unwrap rfl).two
    have hts₁ : (cvS₁ true origTstack isSingle s).tstack = a :: b :: rest := hts
    obtain ⟨t5, rest5, hts5⟩ := one_entry (hok.finish rfl).nonempty
    cases isSingle
    · have key : ((closeVert' curV edgeDir true origTstack false).run s).2 =
          ((finishTstackTop ((maybeUnwrapNxt .R).run (cvS₁ true origTstack false s)).1).run
            (cvS₅ curV edgeDir true origTstack false s)).2 := rfl
      have hR' := hR rfl rfl
      rw [key, finishTstackTop_run_eq _ _ t5 rest5 hts5] at hR' ⊢
      dsimp only at hR' ⊢
      rw [cvS₅_items] at hR' ⊢
      have e2 := maybeUnwrapNxt_run_eq .R (cvS₁ true origTstack false s) a b rest hts₁ _ rfl _ rfl
      simp only [true_or, ↓reduceIte, run_allocItem] at e2
      unfold cvS₂ after at hR' ⊢
      rw [e2] at hR' ⊢
      dsimp only at hR' ⊢
      exact h₁.push_modify_R (fun i _ c hc => st₁.shape.ch_lt i c hc) st₁.shape.no_parent_size _ hR'
    · have key : ((closeVert' curV edgeDir true origTstack true).run s).2 =
          ((finishTstackTop ((maybeUnwrapNxt .S).run (cvS₁ true origTstack true s)).1).run
            (cvS₅ curV edgeDir true origTstack true s)).2 := rfl
      rw [key]
      have h₂ : Items.RSkelInv s.g (cvS₂ true origTstack true s).items := by
        have := keepsR_maybeUnwrap st₁.shape .S (by decide) hts₁
        unfold KeepsR at this
        rw [st₁.g] at this
        exact this h₁
      have hty₂ := maybeUnwrapNxt_result_type st₁.shape .S (hok.unwrap rfl) hts₁
      have hf : ItemFree (cvS₂ true origTstack true s) ((maybeUnwrapNxt .S).run (cvS₁ true origTstack true s)).1 :=
        hfree _ rfl
      have hf₅ := ((hf.merge).merge).retarget curV edgeDir
      have := keepsR_finishTop hts5 _ hf₅.root fun hR' => by
        rw [cvS₅_items] at hR'
        exact absurd (hty₂.symm.trans hR') (by decide)
      unfold KeepsR at this
      rw [hg₅] at this
      exact this (by rw [cvS₅_items]; exact h₂)

/-- `finishEdge` (a returning edge, `lv < d`) keeps `Items.RSkelInv`, given that the pending edge's
Q item is parentless (`hq`) and the two R closes are `RSkel3` (`hR1`: the `.R` iterates of Loop 1;
`hR2`: the type-1 vertex close with `feSingle = false`). -/
theorem keepsR_finishEdge {D : Nat} {s : WalkState} (curV d lv : Nat) (kind : RetKind) (o : DfsOut)
    (origTstack : Nat) (hasVert : Bool) (ho : o.cls = .ret lv kind) (hlow : lv < d) (hv : curV < s.g.nv)
    (hi : s.Inv' D) (hs : Shape s) (hok : FinishOk D curV d lv o origTstack hasVert s)
    (hq : ∀ p, ¬ Items.IsParent s.items p (edgeItem s.g o.e))
    (hR1 : o.cls.isTree = true → ∀ k,
      (∀ j, j ≤ k → result (loop1Cond d) (iter (loop1Body d s.stackDir[d]!) j (ceS₁ o.dest d o.e (feS₀ d o s))) = true) →
      l1Ty d s.stackDir[d]! (iter (loop1Body d s.stackDir[d]!) k (ceS₁ o.dest d o.e (feS₀ d o s))) = .R →
      KeepsR (loop1Body d s.stackDir[d]!) (iter (loop1Body d s.stackDir[d]!) k (ceS₁ o.dest d o.e (feS₀ d o s))))
    (hR2 : o.cls.isTree = true → hasVert = true → o.cls.isType1 = true → feSingle d o s = false →
      Items.RSkel3 s.g (feS₃ curV d o origTstack s).items
        (cvS₁ o.cls.isType1 origTstack (feSingle d o s) (feS₂ d o s)).items.size) :
    KeepsR (finishEdge curV d o origTstack hasVert) s := by
  intro h
  have hlv : o.cls.lowval d = lv := by rw [ho]; rfl
  have hge : ¬ (lv ≥ d) := by omega
  rw [finishEdge_eq]
  simp only [finishEdge', finishTree, finishBack, hlv, hge, ↓reduceIte, WalkM.run_bind, WalkM.get_run,
    run_stackDir, run_makeVs, run_modifyItem]
  have st₀ : Step D curV s (feS₀ d o s) :=
    Step.modifyVs hi hs (edgeItem s.g o.e) _ (by show 1 + s.g.nv + o.e < _; have := hok.e_lt; omega)
  have h₀ : Items.RSkelInv s.g (feS₀ d o s).items := by
    unfold feS₀ after; rw [run_modifyItem]
    exact h.modify_root_of_ne hq _ (fun _ => rfl) (by rw [hs.edge o.e hok.e_lt]; decide)
  have hv₀ : curV < (feS₀ d o s).g.nv := by rw [st₀.g]; exact hv
  by_cases ht : o.cls.isTree = true
  · simp only [ht, ↓reduceIte, closeVert_eq, WalkM.run_bind]
    have st₁ : Step D curV _ (feS₁ d o s) := Step.closeEars st₀.inv st₀.shape hv₀ (hok.ears ht)
    have h₁ : Items.RSkelInv s.g (feS₁ d o s).items := by
      have := keepsR_closeEars st₀.inv st₀.shape hv₀ (hok.ears ht) (hR1 ht)
      unfold KeepsR at this
      rw [st₀.g] at this
      exact this h₀
    have hv₁ : curV < (feS₁ d o s).g.nv := by rw [st₁.g]; exact hv₀
    have st₂ : Step D curV _ (feS₂ d o s) := Step.mergeLate st₁.inv st₁.shape hv₁ (hok.late ht)
    have h₂ : Items.RSkelInv s.g (feS₂ d o s).items := by
      show Items.RSkelInv s.g (after (mergeLate d) (feS₁ d o s)).items
      rw [mergeLate_items]; exact h₁
    have hv₂ : curV < (feS₂ d o s).g.nv := by rw [st₂.g]; exact hv₁
    have hg₂ : (feS₂ d o s).g = s.g := by rw [st₂.g, st₁.g, st₀.g]
    cases hasVert
    · simp only [Bool.false_eq_true, ↓reduceIte]
      have := keepsR_finishRest st₂.inv st₂.shape (hok.rest_tree ht rfl)
      unfold KeepsR at this
      rw [hg₂] at this
      exact this h₂
    · simp only [↓reduceIte, WalkM.run_bind]
      have st₃ : Step D curV _ (feS₃ curV d o origTstack s) :=
        Step.closeVert' st₂.inv st₂.shape hv₂ (hok.vert ht rfl)
      have h₃ : Items.RSkelInv s.g (feS₃ curV d o origTstack s).items := by
        have := keepsR_closeVert' st₂.inv st₂.shape hv₂ (hok.vert ht rfl)
          (fun h1 hB => by rw [hg₂]; exact hR2 ht rfl h1 hB)
        unfold KeepsR at this
        rw [hg₂] at this
        exact this h₂
      have hv₃ : curV < (feS₃ curV d o origTstack s).g.nv := by rw [st₃.g]; exact hv₂
      have := keepsR_finishRest st₃.inv st₃.shape (hok.rest_vert ht rfl)
      unfold KeepsR at this
      rw [st₃.g, hg₂] at this
      exact this h₃
  · have ht' : o.cls.isTree = false := Bool.eq_false_iff.2 ht
    simp only [ht', Bool.false_eq_true, ↓reduceIte, WalkM.run_bind]
    have hq' : Items.ch (feS₀ d o s).items (edgeItem (feS₀ d o s).g o.e) = [] := by
      show Items.ch (s.items.modify _ _) (edgeItem s.g o.e) = []
      rw [Items.ch_modify_ch_eq (edgeItem s.g o.e)
        (fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) (fun _ => rfl)]
      exact hok.q ht'
    have st₁ : Step D curV _ (after (pushEdgeTstack curV lv o.e) (feS₀ d o s)) :=
      Step.pushEdge st₀.inv st₀.shape curV lv o.e hok.e_lt hq' (hok.ends ht') (hok.lv_le ht')
    have st₂ : Step D curV _ (feBack curV lv d o s) := Step.frame st₁.inv st₁.shape rfl rfl rfl rfl
    have hv₂ : curV < (feBack curV lv d o s).g.nv := by rw [st₂.g, st₁.g]; exact hv₀
    have h₂ : Items.RSkelInv s.g (feBack curV lv d o s).items := h₀
    have := keepsR_finishRest st₂.inv st₂.shape (hok.rest_back ht')
    unfold KeepsR at this
    rw [st₂.g, st₁.g, st₀.g] at this
    exact this h₂

end Spqr
