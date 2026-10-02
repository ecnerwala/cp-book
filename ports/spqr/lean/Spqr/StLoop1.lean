import Spqr.StUnwrap
import Spqr.EarLoop1

/-!
# Loop 1 of `finishEdge` keeps the st-reading

The st-side invariant of loop 1 (`closeEars`): the entries above a fixed suffix `B` of the stack read
as fixed pieces `ps` (up to expansion), the item facts hold for fixed `blocks`, and `stackDir` at
depths `≤ d` is kept. The ear-side shape (`L1Ctx`, `L1Inv` of `EarLoop1.lean`) supplies the stack
shape at every iteration; each iteration is `StSim.unwrapMergeClose`.
-/

namespace Spqr
open WalkM WalkState

theorem mem_of_mem_getSide {p : List ItemId × List ItemId} {dir : Bool} {c : ItemId}
    (h : c ∈ getSide p dir) : c ∈ p.1 ++ p.2 := by
  cases dir <;> simp [getSide] at h <;> simp [h]

theorem StItems.congr {g : Graph} {s s' : WalkState} {blocks : List StBlock} (h : StItems g s blocks)
    (hread : readStack s'.tstack = readStack s.tstack) (hitems : s'.items = s.items) :
    StItems g s' blocks := by
  obtain ⟨roots, nodup, bounded, chLt, chNodup, closed⟩ := h
  exact ⟨by rw [hread, hitems]; exact roots, by rw [hread]; exact nodup,
    by rw [hread, hitems]; exact bounded, by rw [hitems]; exact chLt, by rw [hitems]; exact chNodup,
    by rw [hread, hitems]; exact closed⟩

/-- `L1Unwrap` read before `finishEdge` holds at an iteration state `st = c :: R` whose entries `R` are
kept (`L1Keep`) and whose top items are roots or on the original stack. -/
theorem L1Unwrap.transport {s st : WalkState} {R : List TEntry} {ty : NodeType} {t c : TEntry}
    (hK : L1Keep s st R) (htR : t ∈ R) (hts : st.tstack = c :: R) (hR : ∀ u ∈ R, u ∈ s.tstack)
    (hc : ∀ i ∈ c.spans.1 ++ c.spans.2,
      (∀ p, ¬ Items.IsParent s.items p i) ∨ ∃ u ∈ s.tstack, i ∈ u.spans.1 ++ u.spans.2)
    (hbd : ∀ i ∈ t.spans.1 ++ t.spans.2, i < s.items.size) (hpos : 0 < s.items.size)
    (hU : L1Unwrap s ty t) : L1Unwrap st ty t := by
  intro dir h hh hT
  have hlt : h < s.items.size := by
    rcases hside : getSide t.spans dir with _ | ⟨x, xs⟩
    · rw [hside, List.head!_nil] at hh
      subst hh; exact hpos
    · rw [hside] at hh
      simp only [List.head!_cons] at hh
      subst hh
      exact hbd x (mem_of_mem_getSide (dir := dir) (by rw [hside]; exact List.mem_cons_self))
  have hT' : Items.type s.items h = ty := by rw [← hK.type h hlt]; exact hT
  obtain ⟨hsingle, hroot, hch⟩ := hU dir h hh hT'
  have hmem : h ∈ t.spans.1 ++ t.spans.2 :=
    mem_of_mem_getSide (dir := dir) (by rw [hsingle]; exact List.mem_singleton_self _)
  refine ⟨hsingle, fun p hp => hroot p (hK.parent t htR h hmem p hp), ?_⟩
  intro x hx u hu hxu
  rw [hK.ch t htR h hmem] at hx
  rw [hts] at hu
  rcases List.mem_cons.1 hu with rfl | hu
  · rcases hc x hxu with hr | ⟨u', hu', hxu'⟩
    · exact hr h hx
    · exact hch x hx u' hu' hxu'
  · exact hch x hx u (hR u hu) hxu

/-- A fuelled loop runs some number of body iterations, each after a true condition. -/
theorem loop_run_iter (cond : WalkM Bool) (body : WalkM Unit) (hcond : ∀ s, (cond.run s).2 = s)
    (fuel : Nat) (s : WalkState) :
    ∃ k, ((loop fuel cond body).run s).2 = iter body k s ∧
      ∀ j, j < k → (cond.run (iter body j s)).1 = true := by
  induction fuel generalizing s with
  | zero => exact ⟨0, rfl, fun j hj => absurd hj (Nat.not_lt_zero _)⟩
  | succ fuel ih =>
    rw [loop_succ_run fuel cond body s (hcond s)]
    by_cases hc : (cond.run s).1 = true
    · simp only [hc, ↓reduceIte]
      obtain ⟨k, hk, hj⟩ := ih (body.run s).2
      refine ⟨k + 1, hk, fun j hj' => ?_⟩
      cases j with
      | zero => exact hc
      | succ j => exact hj j (Nat.lt_of_succ_lt_succ hj')
    · simp only [hc]
      exact ⟨0, rfl, fun j hj => absurd hj (Nat.not_lt_zero _)⟩

/-- The st-side invariant of loop 1 at an iteration state `st`, relative to the state `s` before
`finishEdge`: the entries above the fixed suffix `B` read as the fixed pieces `ps`, the item facts
hold for the fixed `blocks`, and the directions at depths `≤ d` are unchanged. -/
structure L1StInv (g : Graph) (s : WalkState) (d : Nat) (ps : List StPiece) (blocks : List StBlock)
    (B : List TEntry) (st : WalkState) : Prop where
  read : ∃ new, st.tstack = new ++ B ∧ StRead st.items new ps
  items : StItems g st blocks
  dirs : ∀ k, k ≤ d → st.stackDir[k]! = s.stackDir[k]!

theorem l1St_step {D d : Nat} {o : DfsOut} {s st : WalkState} {hi lo base : List TEntry}
    {g : Graph} {ps : List StPiece} {blocks : List StBlock} {pre B : List TEntry}
    (hc : L1Ctx D d o s hi lo base) (hI : L1Inv D d o s hi lo base st)
    (hcond : result (loop1Cond d) st = true) (hB : base = pre ++ B)
    (hJ : L1StInv g s d ps blocks B st) :
    L1StInv g s d ps blocks B ((loop1Body d s.stackDir[d]!).run st).2 := by
  obtain ⟨done, rest, c, hreach, hP, hK, hF⟩ := hI
  have hts := hP.tstack
  rw [result_loop1Cond_iff] at hcond
  obtain ⟨t, rest', rfl⟩ : ∃ t rest', rest = t :: rest' := by
    cases rest with
    | cons t rest' => exact ⟨t, rest', rfl⟩
    | nil =>
      exfalso
      obtain ⟨l, lo', rfl⟩ : ∃ l lo', lo = l :: lo' := by
        cases lo with
        | nil => exact absurd rfl hc.lo_ne
        | cons l lo' => exact ⟨l, lo', rfl⟩
      have h1 := hc.lo_top l rfl
      have h2 := hcond.2
      rw [hts] at h2
      simp at h2
      omega
  have hpos : 0 < s.items.size := by have := hc.shape.size; omega
  have hmemS : ∀ u ∈ rest' ++ lo ++ base, u ∈ s.tstack := fun u hu => by
    rw [hc.tstack, hreach.append]
    simp only [List.mem_append, List.mem_cons] at hu ⊢
    tauto
  have htHi : t ∈ hi := by rw [hreach.append]; simp
  have htS : t ∈ s.tstack := by rw [hc.tstack]; simp [htHi]
  obtain ⟨new, hnew, hRd⟩ := hJ.read
  have hnew' : new = c :: t :: (rest' ++ lo ++ pre) := by
    apply List.append_cancel_right (bs := B)
    rw [← hnew, hts, hB]
    simp only [List.cons_append, List.append_assoc]
  subst hnew'
  have hnxt : nxtE st = t := by simp [nxtE, hts]
  have hcur : curE st = c := by simp [curE, hts]
  have hroots : ∀ i ∈ c.spans.1 ++ c.spans.2,
      (∀ p, ¬ Items.IsParent s.items p i) ∨ ∃ u ∈ s.tstack, i ∈ u.spans.1 ++ u.spans.2 := hP.root
  rw [loop1Body_run_eq]
  unfold closeAt l1Ty l1S₁ result after
  rw [loop1Type_run, hnxt, hcur]
  split
  · rename_i hgt
    obtain ⟨-, t', rest'', rfl, -, hU⟩ := (hc.spec done t rest' hreach).2 hgt
    have ht'Hi : t' ∈ hi := by rw [hreach.append]; simp
    have ht'S : t' ∈ s.tstack := by rw [hc.tstack]; simp [ht'Hi]
    have hsplit : rest'' ++ lo ++ base = (rest'' ++ lo ++ pre) ++ B := by
      rw [hB]; simp only [List.append_assoc]
    dsimp only
    set stM : WalkState := { st with stackDir := st.stackDir.set! t.topDepth s.stackDir[d]! } with hstM
    have htsM : stM.tstack = c :: t :: (t' :: (rest'' ++ lo ++ base)) := by rw [hstM]; exact hts
    have hs₁ := run_mergeTstackTops_cons_cons stM c t (t' :: (rest'' ++ lo ++ base)) htsM
    generalize (mergeTstackTops.run stM).2 = st₁ at hs₁ ⊢
    have hts₁' : st₁.tstack = TEntry.mergeInto c t :: (t' :: (rest'' ++ lo ++ base)) := by rw [hs₁]
    have hts₁ : st₁.tstack = TEntry.mergeInto c t :: t' :: ((rest'' ++ lo ++ pre) ++ B) := by
      rw [hts₁', hsplit]
    have hitems₁ : st₁.items = st.items := by rw [hs₁, hstM]
    have hsd₁ : st₁.stackDir = st.stackDir.set! t.topDepth s.stackDir[d]! := by rw [hs₁, hstM]
    have hd₁ : st₁.stackDir[d]! = s.stackDir[d]! := by
      rw [hsd₁, Array.getElem!_set!_ne _ _ _ _ (Nat.ne_of_gt hgt)]; exact hF.dir
    have hmin : min t'.topDepth (TEntry.mergeInto c t).topDepth = d := by
      show min t'.topDepth (min t.topDepth c.topDepth) = d
      rw [hP.top, Nat.min_eq_right (Nat.le_of_lt hgt), Nat.min_eq_right (hc.hi_top t' ht'Hi)]
    have hc₁ : getSide (TEntry.mergeInto c t).spans
        (!st₁.stackDir[min t'.topDepth (TEntry.mergeInto c t).topDepth]!) = [] := by
      rw [hmin, hd₁]; exact TEntry.mergeInto_side_nil _ c t hP.side (hc.side t htHi)
    have ht₁ : getSide t'.spans (!st₁.stackDir[min t'.topDepth (TEntry.mergeInto c t).topDepth]!) = [] := by
      rw [hmin, hd₁]; exact hc.side t' ht'Hi
    have hU₁ : L1Unwrap st₁ .S t' := by
      refine L1Unwrap.transport (R := t' :: (rest'' ++ lo ++ base)) ((hK.items_eq hitems₁).mono ?_)
        List.mem_cons_self hts₁' ?_ ?_ (hc.shape.span t' ht'S) hpos hU
      · intro u hu
        simp only [List.mem_cons, List.mem_append] at hu ⊢
        tauto
      · intro u hu
        rcases List.mem_cons.1 hu with rfl | hu
        · exact ht'S
        · exact hmemS u (by simp only [List.cons_append, List.mem_cons]; exact Or.inr hu)
      · intro i hi
        rcases (TEntry.mem_mergeInto_spans c t i).1 hi with hi | hi
        · exact hroots i hi
        · exact Or.inr ⟨t, htS, hi⟩
    have hR₁ : StRead st₁.items (TEntry.mergeInto c t :: t' :: (rest'' ++ lo ++ pre)) ps := by
      unfold StRead; rw [readStack_mergeInto_cons, hitems₁]; exact hRd
    have hI₁ : StItems g st₁ blocks :=
      StItems.congr hJ.items (by rw [hts₁', hts, readStack_mergeInto_cons]; rfl) hitems₁
    have res := StSim.unwrapMergeClose st₁ .S (TEntry.mergeInto c t) t' (rest'' ++ lo ++ pre) B ps blocks
      hts₁ hc₁ ht₁ (Or.inl rfl) (fun _ => hU₁) hR₁ hI₁
    dsimp only at res
    obtain ⟨hts', hsd', -, -, -, hR', hI'⟩ := res
    refine ⟨⟨_, hts', hR'⟩, hI', fun k hk => ?_⟩
    rw [hsd', hsd₁, Array.getElem!_set!_ne _ _ _ _ (by omega)]
    exact hJ.dirs k hk
  · rename_i hle
    have htd : t.topDepth = d := by
      have h2 := hcond.2
      rw [hts] at h2
      simp at h2
      omega
    have hsplit : rest' ++ lo ++ base = (rest' ++ lo ++ pre) ++ B := by
      rw [hB]; simp only [List.append_assoc]
    have hts₀ : st.tstack = c :: t :: ((rest' ++ lo ++ pre) ++ B) := by
      rw [hts]; simp only [List.cons_append, List.append_assoc, hB]
    have hmin : min t.topDepth c.topDepth = d := by rw [hP.top, htd, Nat.min_self]
    have hc₀ : getSide c.spans (!st.stackDir[min t.topDepth c.topDepth]!) = [] := by
      rw [hmin, hF.dir]; exact hP.side
    have ht₀ : getSide t.spans (!st.stackDir[min t.topDepth c.topDepth]!) = [] := by
      rw [hmin, hF.dir]; exact hc.side t htHi
    have hRd' : StRead st.items (c :: t :: (rest' ++ lo ++ pre)) ps := hRd
    have hfin : ∀ ty, ty = .S ∨ ty = .P ∨ ty = .R →
        (ty ≠ .R → L1Unwrap st ty t) →
        L1StInv g s d ps blocks B
          ((finishTstackTop ((maybeUnwrapNxt ty).run st).1).run
            (mergeTstackTops.run ((maybeUnwrapNxt ty).run st).2).2).2 := by
      intro ty hty hU
      have res := StSim.unwrapMergeClose st ty c t (rest' ++ lo ++ pre) B ps blocks hts₀ hc₀ ht₀ hty
        (fun h => hU h) hRd' hJ.items
      dsimp only at res
      obtain ⟨hts', hsd', -, -, -, hR', hI'⟩ := res
      refine ⟨⟨_, hts', hR'⟩, hI', fun k hk => ?_⟩
      rw [hsd']
      exact hJ.dirs k hk
    split
    · rename_i hbeq
      dsimp only
      refine hfin .P (Or.inr (Or.inl rfl)) fun _ => ?_
      have hbot : t.vStart = l1Bot o done := by
        rw [← hP.bot]; exact Nat.eq_of_beq_eq_true hbeq
      have hU := ((hc.spec done t rest' hreach).1 htd).2 hbot
      exact L1Unwrap.transport (R := t :: (rest' ++ lo ++ base)) hK List.mem_cons_self hts
        (fun u hu => by
          rcases List.mem_cons.1 hu with rfl | hu
          · exact htS
          · exact hmemS u hu)
        hroots (hc.shape.span t htS) hpos hU
    · dsimp only
      exact hfin .R (Or.inr (Or.inr rfl)) fun h => absurd rfl h

theorem l1St_iter {D d : Nat} {o : DfsOut} {s : WalkState} {hi lo base : List TEntry} {v₀ : Nat}
    {g : Graph} {ps : List StPiece} {blocks : List StBlock} {pre B : List TEntry}
    (hc : L1Ctx D d o s hi lo base) {st₀ : WalkState}
    (h0 : L1Inv D d o s hi lo base st₀) (hv : v₀ < s.g.nv) (hB : base = pre ++ B)
    (hJ : L1StInv g s d ps blocks B st₀) (k : Nat)
    (hk : ∀ j, j < k → result (loop1Cond d) (iter (loop1Body d s.stackDir[d]!) j st₀) = true) :
    L1StInv g s d ps blocks B (iter (loop1Body d s.stackDir[d]!) k st₀) := by
  induction k with
  | zero => exact hJ
  | succ k ih =>
    rw [iter_succ']
    have hI := l1_iter hc h0 hv k fun j hj => hk j (Nat.lt_succ_of_lt hj)
    exact l1St_step hc hI (hk k (Nat.lt_succ_self k)) hB (ih fun j hj => hk j (Nat.lt_succ_of_lt hj))

/-- Loop 1 (`loop _ (loop1Cond d) (loop1Body d edgeDir)`, `edgeDir = stackDir[d]`) keeps the
st-reading, given the ear-side context and invariant of `EarLoop1.lean`. -/
theorem l1St_loop {D d : Nat} {o : DfsOut} {s : WalkState} {hi lo base : List TEntry} {v₀ : Nat}
    {g : Graph} {ps : List StPiece} {blocks : List StBlock} {pre B : List TEntry}
    (hc : L1Ctx D d o s hi lo base) {st₀ : WalkState}
    (h0 : L1Inv D d o s hi lo base st₀) (hv : v₀ < s.g.nv) (hB : base = pre ++ B)
    (hJ : L1StInv g s d ps blocks B st₀) (fuel : Nat) :
    L1StInv g s d ps blocks B
      ((loop fuel (loop1Cond d) (loop1Body d s.stackDir[d]!)).run st₀).2 := by
  obtain ⟨k, hk, hj⟩ := loop_run_iter (loop1Cond d) (loop1Body d s.stackDir[d]!) (fun _ => rfl) fuel st₀
  rw [hk]
  exact l1St_iter hc h0 hv hB hJ k fun j hj' => hj j hj'

end Spqr
