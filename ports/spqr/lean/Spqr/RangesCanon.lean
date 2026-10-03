import Spqr.WalkWF
import Spqr.RangesCloseTree

/-! # Canonicity of the unternarized walk (`CloseCanon`, `CanonInv`)

`Items.Canonical` (no S under S, no P under P) for `ternarize = false`. Only `finishTstackTop x`
writes the children of an S/P node, so canonicity is kept by every walk primitive except that one,
where it needs the per-site fact `FinishCanon x r`: no child on the closing side has `x`'s type
(the C++ guarantees it through `maybeUnwrapNxt`, which reuses a same-type `nxt` node instead of
nesting). `CloseCanon` states that fact at the three `finishTstackTop` sites of one `finishEdge`
call; `finishEdge_canon` (below) is the step lemma the walk induction (`WalkBackbone.lean`) uses,
and the executable checker (`checks/WalkInvCheck/Ranges.lean`, `checkCanon`) evaluates every field
(`canon_*`). -/

namespace Spqr.WalkState
open WalkM

/-- `finishTstackTop x` from `r` keeps canonicity: if not ternarizing and `x` is S or P, no item on
the closing side of the top entry has `x`'s type. -/
def FinishCanon (x : ItemId) (r : WalkState) : Prop :=
  r.ternarize = false → Items.type r.items x = .S ∨ Items.type r.items x = .P →
    ∀ c ∈ getSide (curE r).spans r.stackDir[(curE r).topDepth]!, Items.type r.items c ≠ Items.type r.items x

/-- The P node `finishP` closes from `r` (when `condP` holds), and the state it closes it at. -/
def pItem (r : WalkState) : ItemId := result (maybeUnwrapNxt .P) r
def pPre (r : WalkState) : WalkState := after mergeTstackTops (after (maybeUnwrapNxt .P) r)

/-- `FinishCanon` at the three `finishTstackTop` sites of one `finishEdge curV d o origTstack hasVert`
call from `s`: the P close (`finishP`, whenever `condP` holds at `feRest`), the type-1 vertex close
(`cvS₅`, for the unwrapped node) and every loop-1 iteration whose first `k` conditions held. -/
structure CloseCanon (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool)
    (s : WalkState) : Prop where
  p_site : o.cls.lowval d < d → o.cls.isType1 = true →
    result (condP curV (o.cls.lowval d) true) (feRest curV d o origTstack hasVert s) = true →
    FinishCanon (pItem (feRest curV d o origTstack hasVert s)) (pPre (feRest curV d o origTstack hasVert s))
  v_site : o.cls.isTree = true → o.cls.lowval d < d → hasVert = true → o.cls.isType1 = true →
    FinishCanon ((maybeUnwrapNxt (if feSingle d o s then NodeType.S else .R)).run (feS₂ d o s)).1
      (cvS₅ curV s.stackDir[d]! true origTstack (feSingle d o s) (feS₂ d o s))
  l1_site : o.cls.isTree = true → o.cls.lowval d < d → ∀ k,
    (∀ j, j ≤ k → result (loop1Cond d) (l1Iter d o s j) = true) →
    FinishCanon (l1Node d o s k) (l1Pre d o s k)

/-! ## Preservation -/

end Spqr.WalkState

namespace Spqr.Items

variable {items : Items}

theorem Canonical.of_types_ch {items' : Items} (h : Canonical items)
    (hty : ∀ p, items'.type p = items.type p) (hch : ∀ p, items'.ch p = items.ch p) :
    Canonical items' := by
  intro p c hpc
  rw [IsParent_congr hch] at hpc
  rw [hty, hty]; exact h p c hpc

theorem Canonical.push (h : Canonical items) (hs : ∀ p c, items.IsParent p c → c < items.size)
    (x : Item) (hx : x.ch = []) : Canonical (items.push x) := by
  intro p c hpc
  have hpc' : items.IsParent p c := by rwa [IsParent, ch_push_nil _ hx] at hpc
  rw [type_push_of_ne _ (Nat.ne_of_lt (parent_lt hpc')), type_push_of_ne _ (Nat.ne_of_lt (hs p c hpc'))]
  exact h p c hpc'

theorem Canonical.modify (h : Canonical items) (j : ItemId) (f : Item → Item)
    (hf : ∀ it, (f it).type = it.type)
    (hch : items.type j = .S ∨ items.type j = .P → ∀ hj : j < items.size, ∀ c ∈ (f items[j]).ch,
      items.type c ≠ items.type j) :
    Canonical (items.modify j f) := by
  intro p c hpc
  rw [type_modify_type_eq j f hf, type_modify_type_eq j f hf]
  rcases IsParent_modify hpc with hpc | ⟨rfl, hj, hc⟩
  · exact h p c hpc
  · by_cases hsp : items.type p = .S ∨ items.type p = .P
    · have hne := hch hsp hj c hc
      exact ⟨fun hcS hpS => hne (hcS.trans hpS.symm), fun hcP hpP => hne (hcP.trans hpP.symm)⟩
    · exact ⟨fun _ hpS => hsp (.inl hpS), fun _ hpP => hsp (.inr hpP)⟩

theorem Canonical.modify_nonSP (h : Canonical items) (j : ItemId) (f : Item → Item)
    (hf : ∀ it, (f it).type = it.type) (hS : items.type j ≠ .S) (hP : items.type j ≠ .P) :
    Canonical (items.modify j f) :=
  h.modify j f hf fun hsp => absurd hsp (by rintro (h | h) <;> contradiction)

theorem Canonical.modify_ch_nonSP (h : Canonical items) (j : ItemId) (g : Item → List ItemId)
    (hS : items.type j ≠ .S) (hP : items.type j ≠ .P) :
    Canonical (items.modify j fun it => { it with ch := g it }) :=
  h.modify_nonSP j _ (fun _ => rfl) hS hP

theorem type_modify_ch' (j : ItemId) (g : Item → List ItemId) (p : ItemId) :
    type (items.modify j fun it => { it with ch := g it }) p = type items p :=
  type_modify items j p (fun it => { it with ch := g it }) fun _ => rfl

theorem Canonical.modify_vs (h : Canonical items) (j : ItemId) (vsv : Option Nat × Option Nat) :
    Canonical (items.modify j fun it => { it with vs := vsv }) :=
  h.of_types_ch (type_modify_vs j vsv) (ch_modify_vs j vsv)

end Spqr.Items

namespace Spqr.WalkState
open WalkM

variable {σ : List Nat} {n D : Nat} {s : WalkState}

theorem CanonInv.frame {s' : WalkState} (h : s.CanonInv) (ht : s'.ternarize = s.ternarize)
    (hi : s'.items = s.items) : s'.CanonInv := fun h' => by
  rw [hi]; exact h (by rw [← ht]; exact h')

theorem CanonInv.alloc (h : s.CanonInv) (hs : Shape s) (ty : NodeType) :
    ({ s with items := s.items.push ⟨ty, (none, none), []⟩ } : WalkState).CanonInv :=
  fun ht => (h ht).push hs.ch_lt _ rfl

theorem CanonInv.modifyVs (h : s.CanonInv) (j : ItemId) (vsv : Option Nat × Option Nat) :
    ({ s with items := s.items.modify j fun it => { it with vs := vsv } } : WalkState).CanonInv :=
  fun ht => (h ht).modify_vs j vsv

theorem CanonInv.mergeTop (h : s.CanonInv) : (after mergeTstackTops s).CanonInv := by
  rw [after_mergeTstackTops]; exact h.frame rfl rfl

theorem CanonInv.loop (cond : WalkM Bool) (body : WalkM Unit) (fuel : Nat)
    (hcond : ∀ s, (cond.run s).2 = s)
    (hbody : ∀ k, (∀ j, j ≤ k → (cond.run (iter body j s)).1 = true) →
      (iter body k s).CanonInv → (after body (iter body k s)).CanonInv)
    (h : s.CanonInv) : (after (WalkM.loop fuel cond body) s).CanonInv := by
  induction fuel generalizing s with
  | zero => exact h
  | succ fuel ih =>
    show ((WalkM.loop (fuel + 1) cond body).run s).2.CanonInv
    rw [loop_succ_run fuel cond body s (hcond s)]
    by_cases hc : (cond.run s).1 = true
    · simp only [hc, ↓reduceIte]
      have h₁ := hbody 0 (fun j hj => by rw [Nat.le_zero.1 hj]; exact hc) h
      exact ih (s := (body.run s).2) (fun k hk => hbody (k + 1) fun j hj => by
        cases j with
        | zero => exact hc
        | succ j => exact hk j (Nat.le_of_succ_le_succ hj)) h₁
    · simp only [hc]; exact h

theorem CanonInv.mergeLoop (h : s.CanonInv) (cond : WalkM Bool) (hcond : ∀ s, (cond.run s).2 = s)
    (fuel : Nat) : (after (WalkM.loop fuel cond mergeTstackTops) s).CanonInv :=
  CanonInv.loop cond mergeTstackTops fuel hcond (fun _ _ hk => hk.mergeTop) h

theorem CanonInv.mergeLate (h : s.CanonInv) (d : Nat) : (after (Spqr.mergeLate d) s).CanonInv := by
  show ((Spqr.mergeLate d).run s).2.CanonInv
  rw [mergeLate_run]
  by_cases hc : (curE s).firstIdx > s.firstOccurrence[d]!
  · simp only [hc, ↓reduceIte]
    exact h.mergeLoop _ (fun _ => rfl) _
  · simp only [hc, ↓reduceIte]; exact h

theorem CanonInv.vertPre (h : s.CanonInv) (isType1 : Bool) (origTstack : Nat) (isSingle : Bool) :
    (cvS₁ isType1 origTstack isSingle s).CanonInv := by
  cases isType1
  · show ((WalkState.vertPre false origTstack isSingle).run s).2.CanonInv
    simp only [WalkState.vertPre, Bool.not_false, ↓reduceIte, WalkM.run_bind, run_tstackSize,
      WalkM.pure_run]
    exact h.mergeLoop _ (fun _ => rfl) _
  · exact h

theorem CanonInv.maybeUnwrap (h : s.CanonInv) (hs : Shape s) (ty : NodeType) {a b : TEntry}
    {rest : List TEntry} (hts : s.tstack = a :: b :: rest) :
    (after (maybeUnwrapNxt ty) s).CanonInv := by
  show ((maybeUnwrapNxt ty).run s).2.CanonInv
  rw [maybeUnwrapNxt_run_eq ty s a b rest hts _ rfl _ rfl]
  split
  · rw [run_allocItem]; exact h.alloc hs ty
  · split
    · exact h.frame rfl rfl
    · rw [run_allocItem]; exact h.alloc hs ty

theorem CanonInv.vertUnwrap (h : s.CanonInv) (hs : Shape s) {isType1 isSingle : Bool}
    (hok : isType1 = true → UnwrapOk (if isSingle then .S else .R) s) :
    (after (vertUnwrap isType1 isSingle) s).CanonInv := by
  cases isType1
  · exact h
  · obtain ⟨a, b, rest, hts⟩ := two_entries_of_le (hok rfl).two
    show ((some <$> maybeUnwrapNxt _).run s).2.CanonInv
    rw [WalkM.map_run]
    exact h.maybeUnwrap hs _ hts

theorem CanonInv.retarget (h : s.CanonInv) (curV : Nat) (dir : Bool) :
    (after (WalkState.retarget curV dir) s).CanonInv := by
  show ((WalkState.retarget curV dir).run s).2.CanonInv
  rw [WalkState.retarget, run_modifyCur]; exact h.frame rfl rfl

theorem CanonInv.pushVert (h : s.CanonInv) (v d : Nat) : (after (pushVertTstack v d) s).CanonInv := by
  show ((pushVertTstack v d).run s).2.CanonInv
  rw [pushVertTstack, run_pushTstack]; exact h.frame rfl rfl

theorem CanonInv.finishTail (h : s.CanonInv) (curV d : Nat) (hasVert isSingle : Bool) :
    (after (Spqr.finishTail curV d hasVert isSingle) s).CanonInv := by
  cases hasVert
  · have h₁ := h.pushVert curV d
    cases isSingle
    · show (after mergeTstackTops (after (pushVertTstack curV d) s)).CanonInv
      exact h₁.mergeTop
    · exact h₁
  · exact h

theorem CanonInv.finishTop (h : s.CanonInv) (x : ItemId) (hc : FinishCanon x s) {t : TEntry}
    {rest : List TEntry} (hts : s.tstack = t :: rest) : (after (finishTstackTop x) s).CanonInv := by
  show ((finishTstackTop x).run s).2.CanonInv
  rw [finishTstackTop_run_eq s x t rest hts]
  refine fun ht => (h ht).modify x _ (fun _ => rfl) fun hsp _ c hc' => hc ht hsp c ?_
  simpa [curE, hts] using hc'

theorem CanonInv.pSite (h : s.CanonInv) (hs : Shape s) {curV lv : Nat} {b : Bool}
    (hok : result (condP curV lv b) s = true → FinishCanon (pItem s) (pPre s)) :
    (after (finishP curV lv b) s).CanonInv := by
  show ((finishP curV lv b).run s).2.CanonInv
  simp only [finishP, WalkM.run_bind, run_condP]
  by_cases hc : result (condP curV lv b) s = true
  · have hc' : (b && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == curV) &&
        (s.tstack.tail.head!.topDepth == lv)) = true := hc
    simp only [hc', ↓reduceIte, WalkM.run_bind]
    have h2 : 2 ≤ s.tstack.length :=
      of_decide_eq_true (Bool.and_eq_true_iff.1 (Bool.and_eq_true_iff.1 (Bool.and_eq_true_iff.1 hc').1).1).2
    obtain ⟨a, b', rest, hts⟩ := two_entries_of_le h2
    have h₁ : (after (maybeUnwrapNxt .P) s).CanonInv := h.maybeUnwrap hs .P hts
    have h₂ : (pPre s).CanonInv := h₁.mergeTop
    obtain ⟨y', items, heq, -, -, -⟩ := maybeUnwrapNxt_run .P s hts
    have hts₁ : (after (maybeUnwrapNxt .P) s).tstack = a :: y' :: rest := by
      show ((maybeUnwrapNxt .P).run s).2.tstack = _; rw [heq]
    obtain ⟨y'', heq₂, -, -, -⟩ := mergeTstackTops_run _ hts₁
    have hts₂ : (pPre s).tstack = y'' :: rest := by
      show (mergeTstackTops.run _).2.tstack = _; rw [heq₂]
    exact h₂.finishTop _ (hok hc) hts₂
  · have hc' : (b && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == curV) &&
        (s.tstack.tail.head!.topDepth == lv)) = false := Bool.eq_false_iff.2 hc
    simp only [hc', Bool.false_eq_true, ↓reduceIte]
    exact h

theorem CanonInv.l1Type (h : s.CanonInv) (d : Nat) (dir : Bool) : (l1S₁ d dir s).CanonInv := by
  show ((loop1Type d dir).run s).2.CanonInv
  rw [loop1Type_run]
  split
  · show (after mergeTstackTops { s with stackDir := s.stackDir.set! (nxtE s).topDepth dir }).CanonInv
    exact (h.frame (s' := { s with stackDir := s.stackDir.set! (nxtE s).topDepth dir }) rfl rfl).mergeTop
  · split <;> exact h

/-- One loop-1 iteration from its merged pre-close state. -/
theorem CanonInv.l1Body (h : s.CanonInv) (d : Nat) (dir : Bool) (hs : Shape (l1S₁ d dir s))
    (hok : UnwrapOk (l1Ty d dir s) (l1S₁ d dir s))
    (hc : FinishCanon (result (maybeUnwrapNxt (l1Ty d dir s)) (l1S₁ d dir s))
      (after mergeTstackTops (l1S₂ d dir s))) :
    (after (loop1Body d dir) s).CanonInv := by
  obtain ⟨a, b, rest, hts⟩ := two_entries_of_le hok.two
  have h₁ := (h.l1Type d dir).maybeUnwrap hs (l1Ty d dir s) hts
  have h₂ : (after mergeTstackTops (l1S₂ d dir s)).CanonInv := h₁.mergeTop
  obtain ⟨y', items, heq, -, -, -⟩ := maybeUnwrapNxt_run (l1Ty d dir s) _ hts
  have hts₁ : (l1S₂ d dir s).tstack = a :: y' :: rest := by
    show ((maybeUnwrapNxt _).run _).2.tstack = _; rw [heq]
  obtain ⟨y'', heq₂, -, -, -⟩ := mergeTstackTops_run _ hts₁
  have hts₂ : (after mergeTstackTops (l1S₂ d dir s)).tstack = y'' :: rest := by
    show (mergeTstackTops.run _).2.tstack = _; rw [heq₂]
  exact h₂.finishTop _ hc hts₂

/-- The boundary branch: only `vs` writes, allocations, and `ch` writes to the Q and V items. -/
theorem CanonInv.boundary (h : s.CanonInv) (hs : Shape s) (curV d : Nat) (o : DfsOut)
    (hasVert : Bool) (he : o.e < s.g.ne) (hv : curV < s.g.nv) :
    (after (finishBoundary curV d o (edgeItem s.g o.e) hasVert) s).CanonInv := by
  have hq : Items.type s.items (edgeItem s.g o.e) = .Q := hs.edge o.e he
  have hV : Items.type s.items (vertItem curV) = .V := hs.vert curV hv
  have hqsz : edgeItem s.g o.e < s.items.size := by
    show 1 + s.g.nv + o.e < _; have := hs.size; omega
  have hvsz : vertItem curV < s.items.size := by
    show 1 + curV < _; have := hs.size; omega
  have hqne : edgeItem s.g o.e ≠ s.items.size := Nat.ne_of_lt hqsz
  have hvne : vertItem curV ≠ s.items.size := Nat.ne_of_lt hvsz
  have hpq : ∀ (vs₀ : Option Nat × Option Nat) (x : Item),
      Items.type ((s.items.modify (edgeItem s.g o.e) fun it => { it with vs := vs₀ }).push x)
        (edgeItem s.g o.e) = .Q := by
    intro vs₀ x
    rw [Items.type_push_of_ne _ (by rw [Array.size_modify]; exact hqne), Items.type_modify_vs, hq]
  have hpv : ∀ (vs₀ : Option Nat × Option Nat) (x : Item),
      Items.type ((s.items.modify (edgeItem s.g o.e) fun it => { it with vs := vs₀ }).push x)
        (vertItem curV) = .V := by
    intro vs₀ x
    rw [Items.type_push_of_ne _ (by rw [Array.size_modify]; exact hvne), Items.type_modify_vs, hV]
  show wp (finishBoundary curV d o (edgeItem s.g o.e) hasVert) (fun _ s' => s'.CanonInv) s
  unfold finishBoundary
  simp only [wp_bind, wp_modifyItem, wp_modify, wp_ite, wp_allocItem, wp_makeVs, wp_popTstack,
    wp_pure]
  have h₀ := h.modifyVs (edgeItem s.g o.e) (some curV, none)
  have hs' : ∀ p c, Items.IsParent (s.items.modify (edgeItem s.g o.e) fun it =>
      { it with vs := (some curV, none) }) p c →
      c < (s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) }).size := by
    intro p c hpc
    rw [Array.size_modify]
    exact hs.ch_lt p c ((Items.IsParent_congr (Items.ch_modify_vs _ _)).1 hpc)
  split
  · split
    · intro ht
      dsimp only at ht ⊢
      simp only [Array.size_modify]
      refine (((h₀ ht).push hs' _ rfl).modify_vs _ _).modify_ch_nonSP _ _ ?_ ?_
        |>.modify_ch_nonSP _ _ ?_ ?_
      · rw [Items.type_modify_vs, hpq]; decide
      · rw [Items.type_modify_vs, hpq]; decide
      · rw [Items.type_modify_ch', Items.type_modify_vs, hpv]; decide
      · rw [Items.type_modify_ch', Items.type_modify_vs, hpv]; decide
    · intro ht
      dsimp only at ht ⊢
      refine ((h₀ ht).modify_ch_nonSP _ _ ?_ ?_).modify_ch_nonSP _ _ ?_ ?_
      · rw [Items.type_modify_vs, hq]; decide
      · rw [Items.type_modify_vs, hq]; decide
      · rw [Items.type_modify_ch', Items.type_modify_vs, hV]; decide
      · rw [Items.type_modify_ch', Items.type_modify_vs, hV]; decide
  · intro ht
    dsimp only at ht ⊢
    simp only [Array.size_modify]
    refine (((h₀ ht).push hs' _ rfl).modify_vs _ _).modify_ch_nonSP _ _ ?_ ?_
      |>.modify_ch_nonSP _ _ ?_ ?_
    · rw [Items.type_modify_vs, hpq]; decide
    · rw [Items.type_modify_vs, hpq]; decide
    · rw [Items.type_modify_ch', Items.type_modify_vs, hpv]; decide
    · rw [Items.type_modify_ch', Items.type_modify_vs, hpv]; decide

theorem CanonInv.rest (h : s.CanonInv) (hs : Shape s) {curV d lv : Nat}
    {isType1 hasVert isSingle : Bool}
    (hok : result (condP curV lv isType1) s = true → FinishCanon (pItem s) (pPre s)) :
    (after (finishRest curV d lv isType1 hasVert isSingle) s).CanonInv := by
  show ((finishRest curV d lv isType1 hasVert isSingle).run s).2.CanonInv
  simp only [finishRest, WalkM.run_bind]
  exact (h.pSite hs hok).finishTail curV d hasVert isSingle

/-- `finishEdge` keeps canonicity under `CloseBase` and the three site facts. -/
theorem finishEdge_canon (h : CloseBase σ n D curV d o origTstack hasVert s)
    (hcn : CloseCanon curV d o origTstack hasVert s) (hc : s.CanonInv) :
    (after (finishEdge curV d o origTstack hasVert) s).CanonInv := by
  have hs := h.rgs.2.1
  by_cases hge : d ≤ o.cls.lowval d
  · have hge' : o.cls.lowval d ≥ d := hge
    show ((finishEdge curV d o origTstack hasVert).run s).2.CanonInv
    rw [finishEdge_eq]
    simp only [finishEdge', WalkM.run_bind, WalkM.get_run, run_stackDir, hge', ↓reduceIte]
    exact hc.boundary hs curV d o hasVert h.book.e_lt h.book.v_lt
  have hlow : o.cls.lowval d < d := Nat.lt_of_not_le hge
  obtain ⟨lv, kind, ho, hl⟩ := ret_of_lowval_lt hlow
  have hok := h.finishOk ho hl
  have hlv : o.cls.lowval d = lv := by rw [ho]; rfl
  have hge' : ¬ (lv ≥ d) := by omega
  have hv := h.book.v_lt
  have hj : edgeItem s.g o.e < 1 + s.g.nv + s.g.ne := by
    show 1 + s.g.nv + o.e < _; have := hok.e_lt; omega
  have st₀ : Step D curV s (feS₀ d o s) := Step.modifyVs h.rgs.1.inv hs (edgeItem s.g o.e) _ hj
  have hv₀ : curV < (feS₀ d o s).g.nv := by rw [st₀.g]; exact hv
  have h₀ : (feS₀ d o s).CanonInv := hc.modifyVs _ _
  have pSite : ∀ r, feRest curV d o origTstack hasVert s = r →
      result (condP curV lv o.cls.isType1) r = true → FinishCanon (pItem r) (pPre r) := by
    intro r hr hcond
    have hc' : (o.cls.isType1 && decide (r.tstack.length ≥ 2) && (r.tstack.tail.head!.vStart == curV) &&
        (r.tstack.tail.head!.topDepth == lv)) = true := hcond
    have hb : o.cls.isType1 = true :=
      (Bool.and_eq_true_iff.1 (Bool.and_eq_true_iff.1 (Bool.and_eq_true_iff.1 hc').1).1).1
    rw [hb] at hcond
    have := hcn.p_site hlow hb (by rw [hr, hlv]; exact hcond)
    rwa [hr] at this
  show ((finishEdge curV d o origTstack hasVert).run s).2.CanonInv
  rw [finishEdge_eq]
  simp only [finishEdge', finishTree, finishBack, hlv, hge', ↓reduceIte, WalkM.run_bind, WalkM.get_run,
    run_stackDir, run_makeVs, run_modifyItem]
  by_cases ht : o.cls.isTree = true
  · simp only [ht, ↓reduceIte, closeVert_eq, WalkM.run_bind]
    have st₁ : Step D curV _ (feS₁ d o s) := Step.closeEars st₀.inv st₀.shape hv₀ (hok.ears ht)
    have hv₁ : curV < (feS₁ d o s).g.nv := by rw [st₁.g]; exact hv₀
    have st₂ : Step D curV _ (feS₂ d o s) := Step.mergeLate st₁.inv st₁.shape hv₁ (hok.late ht)
    have hv₂ : curV < (feS₂ d o s).g.nv := by rw [st₂.g]; exact hv₁
    have h₁ : (feS₁ d o s).CanonInv := by
      show (after (loop _ (loop1Cond d) (loop1Body d s.stackDir[d]!))
        (ceS₁ o.dest d o.e (feS₀ d o s))).CanonInv
      refine CanonInv.loop _ _ _ (fun _ => rfl) (fun k hk hk' => ?_) ?_
      · have hok₁ : Loop1BodyOk D d s.stackDir[d]! (l1Iter d o s k) := (hok.ears ht).body k hk
        have hadj : Loop1BodyAdj σ d s.stackDir[d]! (l1Iter d o s k) :=
          ((h.finishR.1 lv kind ho hl).ears ht).body k hk
        have st := rangesInv_l1Iter h ht hlow k (fun j hj => hk j hj.le)
        have st₁ : RgStep σ (n + 1) D curV _ (l1S₁ d s.stackDir[d]! (l1Iter d o s k)) :=
          RgStep.loop1Type st.ranges st.step.shape h.nodup (st.hσ h.rgs.2.2) hok₁.mergeS hadj.mergeS
        exact hk'.l1Body d _ st₁.step.shape hok₁.unwrap (hcn.l1_site ht hlow k hk)
      · show ((pushEdgeTstack o.dest d o.e).run (feS₀ d o s)).2.CanonInv
        rw [run_pushEdgeTstack]; exact h₀.frame rfl rfl
    have h₂ : (feS₂ d o s).CanonInv := h₁.mergeLate d
    cases hasVert
    · simp only [Bool.false_eq_true, ↓reduceIte]
      have hr : feRest curV d o origTstack false s = feS₂ d o s := by simp [feRest, ht]
      exact h₂.rest (st₀.trans (st₁.trans st₂)).shape (pSite _ hr)
    · simp only [↓reduceIte, WalkM.run_bind]
      have hcv := hok.vert ht rfl
      have st₃ : Step D curV _ (cvS₁ o.cls.isType1 origTstack (feSingle d o s) (feS₂ d o s)) :=
        Step.vertPre st₂.inv st₂.shape hv₂ hcv.loop3
      have h₃ := h₂.vertPre o.cls.isType1 origTstack (feSingle d o s)
      have h₄ : (cvS₂ o.cls.isType1 origTstack (feSingle d o s) (feS₂ d o s)).CanonInv :=
        h₃.vertUnwrap st₃.shape (fun h => by rw [h]; exact hcv.unwrap h)
      have h₅ : (cvS₃ o.cls.isType1 origTstack (feSingle d o s) (feS₂ d o s)).CanonInv := by
        show (after mergeTstackTops _).CanonInv; exact h₄.mergeTop
      have h₆ : (cvS₄ o.cls.isType1 origTstack (feSingle d o s) (feS₂ d o s)).CanonInv := by
        show (after mergeTstackTops _).CanonInv; exact h₅.mergeTop
      have h₇ := h₆.retarget curV s.stackDir[d]!
      have h₈ : (feS₃ curV d o origTstack s).CanonInv := by
        cases h1 : o.cls.isType1
        · have hS₃ : feS₃ curV d o origTstack s =
              cvS₅ curV s.stackDir[d]! false origTstack (feSingle d o s) (feS₂ d o s) := by
            simp only [feS₃, h1, closeVert', cvS₅, cvS₄, cvS₃, cvS₂, cvS₁, after, vertPre,
              vertUnwrap, vertFinish, Bool.not_false, Bool.false_eq_true, ↓reduceIte, WalkM.run_bind,
              WalkM.pure_run, run_tstackSize]
          rw [hS₃]; rw [h1] at h₇; exact h₇
        · have hS₃ : feS₃ curV d o origTstack s =
              after (finishTstackTop
                ((maybeUnwrapNxt (if feSingle d o s then .S else .R)).run (feS₂ d o s)).1)
                (cvS₅ curV s.stackDir[d]! true origTstack (feSingle d o s) (feS₂ d o s)) := by
            simp only [feS₃, h1, closeVert', cvS₅, cvS₄, cvS₃, cvS₂, cvS₁, cvB₁, after, result, vertPre,
              vertUnwrap, vertFinish, Bool.not_true, Bool.false_eq_true, ↓reduceIte, WalkM.run_bind,
              WalkM.pure_run, WalkM.map_run]
            rfl
          have hne := (hcv.finish h1).nonempty
          rw [h1] at h₇ hne
          obtain ⟨t, rest, hts⟩ := List.exists_cons_of_ne_nil hne
          rw [hS₃]
          exact h₇.finishTop _ (hcn.v_site ht hlow rfl h1) hts
      have st₈ : Step D curV _ (feS₃ curV d o origTstack s) :=
        Step.closeVert' st₂.inv st₂.shape hv₂ hcv
      have hr : feRest curV d o origTstack true s = feS₃ curV d o origTstack s := by simp [feRest, ht]
      exact h₈.rest (st₀.trans (st₁.trans (st₂.trans st₈))).shape (pSite _ hr)
  · have ht' : o.cls.isTree = false := Bool.eq_false_iff.2 ht
    simp only [ht', Bool.false_eq_true, ↓reduceIte, WalkM.run_bind]
    have hq : Items.ch (feS₀ d o s).items (edgeItem (feS₀ d o s).g o.e) = [] := by
      show Items.ch (s.items.modify _ _) (edgeItem s.g o.e) = []
      rw [Items.ch_modify_ch_eq (edgeItem s.g o.e)
        (fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) })
        (fun _ => rfl)]
      exact hok.q ht'
    have st₁ : Step D curV _ (after (pushEdgeTstack curV lv o.e) (feS₀ d o s)) :=
      Step.pushEdge st₀.inv st₀.shape curV lv o.e hok.e_lt hq (hok.ends ht') (hok.lv_le ht')
    have st₂ : Step D curV _ (feBack curV lv d o s) :=
      Step.frame st₁.inv st₁.shape rfl rfl rfl rfl
    have hB : (feBack curV lv d o s).CanonInv := by
      have h₁ : (after (pushEdgeTstack curV lv o.e) (feS₀ d o s)).CanonInv := by
        show ((pushEdgeTstack curV lv o.e).run (feS₀ d o s)).2.CanonInv
        rw [run_pushEdgeTstack]; exact h₀.frame rfl rfl
      exact h₁.frame rfl rfl
    have hr : feRest curV d o origTstack hasVert s = feBack curV lv d o s := by simp [feRest, ht', hlv]
    exact hB.rest (st₀.trans (st₁.trans st₂)).shape (pSite _ hr)

end Spqr.WalkState
