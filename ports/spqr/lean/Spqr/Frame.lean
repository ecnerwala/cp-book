import Spqr.Sim
import Mathlib.Data.List.Basic

/-!
`finishEdge` split into named blocks (`finishEdge_eq`), and the simulation lemmas showing each block
commutes with `lift bot` when the entries of `bot` are strictly below the current depth.
-/

namespace Spqr
open WalkM

/-! ### Blocks of `finishEdge` -/

def loop1Cond (d : Nat) : WalkM Bool := do return (← tstackSize) ≥ 2 && (← nxt).topDepth ≥ d

/-- Classify the ear below the top (merging the vertex entry of an S-ear into it). -/
def loop1Type (d : Nat) (edgeDir : Bool) : WalkM NodeType := do
  if (← nxt).topDepth > d then
    setStackDir (← nxt).topDepth edgeDir
    mergeTstackTops
    pure .S
  else if (← nxt).vStart == (← cur).vStart then pure .P
  else pure .R

def loop1Body (d : Nat) (edgeDir : Bool) : WalkM Unit := do
  let item ← maybeUnwrapNxt (← loop1Type d edgeDir)
  mergeTstackTops
  finishTstackTop item

/-- Push the tree edge and close the ears whose top is at or below `d`. -/
def closeEars (nxtV d e : Nat) (edgeDir : Bool) : WalkM Unit := do
  pushEdgeTstack nxtV d e
  loop (← tstackSize) (loop1Cond d) (loop1Body d edgeDir)

def loop2Cond (fo : Nat) : WalkM Bool := do return (← cur).firstIdx > fo

/-- Merge the ears whose first back edge comes after one of `cur`'s. Returns `isSingle`. -/
def mergeLate (d : Nat) : WalkM Bool := do
  let fo := (← get).firstOccurrence[d]!
  if (← cur).firstIdx > fo then
    loop (← tstackSize) (loop2Cond fo) mergeTstackTops
    return false
  return true

def loop3Cond (origTstack : Nat) : WalkM Bool := do return (← tstackSize) > origTstack + 3

/-- The `hasVert` case: wrap `cur` into the pending vertex ear, then continue with `k isSingle`. -/
def closeVert (curV : Nat) (edgeDir isType1 : Bool) (origTstack : Nat) (isSingle : Bool)
    (k : Bool → WalkM Bool) : WalkM Bool := do
  let mut isSingle := isSingle
  if !isType1 then
    loop (← tstackSize) (loop3Cond origTstack) mergeTstackTops
    isSingle := false
  let item ← if isType1 then some <$> maybeUnwrapNxt (if isSingle then .S else .R) else pure none
  mergeTstackTops
  mergeTstackTops
  modifyCur fun t => { t with vStart := curV, spans := setSides (!edgeDir) (t.spans.1 ++ t.spans.2) [] }
  match item with
  | some item => finishTstackTop item; k true
  | none => k isSingle

/-- Is the ear below the top the P-ear of a parallel back edge? -/
def condP (curV lowval : Nat) (isType1 : Bool) : WalkM Bool := do
  return isType1 && (← tstackSize) ≥ 2 && (← nxt).vStart == curV && (← nxt).topDepth == lowval

def finishP (curV lowval : Nat) (isType1 : Bool) : WalkM Unit := do
  if ← condP curV lowval isType1 then
    let item ← maybeUnwrapNxt .P
    mergeTstackTops
    finishTstackTop item

/-- The first-edge vertex push. -/
def finishTail (curV d : Nat) (hasVert isSingle : Bool) : WalkM Bool := do
  if !hasVert then
    pushVertTstack curV d
    if !isSingle then mergeTstackTops
    return true
  return hasVert

/-- The tail of `finishEdge`: the P-check and the first-edge vertex push. -/
def finishRest (curV d lowval : Nat) (isType1 hasVert isSingle : Bool) : WalkM Bool := do
  finishP curV lowval isType1
  finishTail curV d hasVert isSingle

/-- The `lowval ≥ d` case of `finishEdge`. -/
def finishBoundary (curV d : Nat) (o : DfsOut) (qItem : ItemId) (hasVert : Bool) : WalkM Bool := do
  modifyItem qItem fun it => { it with vs := (some curV, none) }
  modify fun s => { s with totBlocks := s.totBlocks + 1 }
  if o.cls.isTree then
    if o.cls.lowval d == d + 1 then
      let item ← allocItem .I
      let vs ← makeVs o.dest d
      modifyItem item fun it => { it with vs := vs }
      let t ← popTstack
      modifyItem qItem fun it => { it with ch := item :: t.spans.2 }
    else
      let backedge ← popTstack
      let t ← popTstack
      modifyItem qItem fun it => { it with ch := backedge.spans.1 ++ t.spans.2 }
  else
    modify fun s => { s with totSelfLoops := s.totSelfLoops + 1 }
    let item ← allocItem .O
    modifyItem item fun it => { it with vs := (some curV, none) }
    modifyItem qItem fun it => { it with ch := [item] }
  modifyItem (vertItem curV) fun it => { it with ch := it.ch ++ [qItem] }
  return hasVert

def finishEdge' (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool) : WalkM Bool := do
  let g := (← get).g
  let lowval := o.cls.lowval d
  let isType1 := o.cls.isType1
  let edgeDir ← stackDir d
  if lowval ≥ d then finishBoundary curV d o (edgeItem g o.e) hasVert
  else
    let vs ← makeVs o.dest d
    modifyItem (edgeItem g o.e) fun it => { it with vs := vs }
    if o.cls.isTree then
      closeEars o.dest d o.e edgeDir
      let isSingle ← mergeLate d
      if hasVert then
        closeVert curV edgeDir isType1 origTstack isSingle (finishRest curV d lowval isType1 hasVert)
      else finishRest curV d lowval isType1 hasVert isSingle
    else
      pushEdgeTstack curV lowval o.e
      modify fun s =>
        { s with firstOccurrence := s.firstOccurrence.modify lowval (min · s.nxtEdgeIdx),
                 nxtEdgeIdx := s.nxtEdgeIdx + 1 }
      finishRest curV d lowval isType1 hasVert true

theorem ite_bind (c : Prop) [Decidable c] (a b : WalkM α) (f : α → WalkM β) :
    (if c then a else b) >>= f = if c then a >>= f else b >>= f := by split <;> rfl

theorem finishEdge_eq (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool) :
    finishEdge curV d o origTstack hasVert = finishEdge' curV d o origTstack hasVert := by
  simp only [finishEdge, finishEdge', finishBoundary, closeEars, loop1Cond, loop1Type, loop1Body, mergeLate,
    loop2Cond, loop3Cond, closeVert, finishRest, condP, finishP, finishTail, bind_assoc, Bool.false_eq_true, Bool.not_false,
    Bool.not_true, ite_true, ite_false, pure_bind, ite_bind]
  rfl

/-! ### Stack shapes -/

/-- The bottom entry (if any) is at or above depth `d`. -/
def Inv1 (d : Nat) (l : List TEntry) : Prop := ∀ e ∈ l.getLast?, e.topDepth ≤ d

/-- The bottom entry (if any) was opened at or before edge index `fo`. -/
def Inv2 (fo : Nat) (l : List TEntry) : Prop := ∀ e ∈ l.getLast?, e.firstIdx ≤ fo

theorem Inv1.of_le {d : Nat} {x y y' : TEntry} {rest : List TEntry} (h : Inv1 d (x :: y :: rest))
    (hy : y'.topDepth ≤ y.topDepth) : Inv1 d (y' :: rest) := by
  intro e he
  cases rest with
  | nil => simp at he; subst he; exact Nat.le_trans hy (h y (by simp))
  | cons z rest => exact h e (by simpa using he)

theorem Inv1.cons_of_le {d : Nat} {x x' : TEntry} {rest : List TEntry} (h : Inv1 d (x :: rest))
    (hx : x'.topDepth ≤ x.topDepth) : Inv1 d (x' :: rest) := by
  intro e he
  cases rest with
  | nil => simp at he; subst he; exact Nat.le_trans hx (h x (by simp))
  | cons z rest => exact h e (by simpa using he)

theorem Inv1.cons_cons_of_le {d : Nat} {x y y' : TEntry} {rest : List TEntry}
    (h : Inv1 d (x :: y :: rest)) (hy : y'.topDepth ≤ y.topDepth) : Inv1 d (x :: y' :: rest) := by
  intro e he
  cases rest with
  | nil => simp at he; subst he; exact Nat.le_trans hy (h y (by simp))
  | cons z rest => exact h e (by simpa using he)

theorem Inv1.push {d : Nat} {x : TEntry} {l : List TEntry} (h : Inv1 d l) (hx : x.topDepth ≤ d) :
    Inv1 d (x :: l) := by
  intro e he
  cases l with
  | nil => simp at he; subst he; exact hx
  | cons z rest => exact h e (by simpa using he)

theorem Inv2.of_eq {fo : Nat} {x y y' : TEntry} {rest : List TEntry} (h : Inv2 fo (x :: y :: rest))
    (hy : y'.firstIdx = y.firstIdx) : Inv2 fo (y' :: rest) := by
  intro e he
  cases rest with
  | nil => simp at he; subst he; exact hy ▸ h y (by simp)
  | cons z rest => exact h e (by simpa using he)

/-! ### Run lemmas for the compound primitives -/

theorem idBind (x : Id α) (f : α → Id β) : x >>= f = f x := rfl

def mergeTops : List TEntry → List TEntry
  | b :: a :: rest =>
    { a with topDepth := min a.topDepth b.topDepth,
             spans := (b.spans.1 ++ a.spans.1, a.spans.2 ++ b.spans.2) } :: rest
  | _ => []

theorem run_mergeTstackTops (s : WalkState) :
    mergeTstackTops.run s = ((), { s with tstack := mergeTops s.tstack }) := by
  rcases s with ⟨g, tern, items, sv, sd, nei, fo, ts, tb, tsl⟩
  match ts with
  | [] | [_] | _ :: _ :: _ => rfl

theorem maybeUnwrapNxt_run (t : NodeType) (s : WalkState) {x y : TEntry} {rest : List TEntry}
    (h : s.tstack = x :: y :: rest) :
    ∃ y' items, ((maybeUnwrapNxt t).run s).2 = { s with tstack := x :: y' :: rest, items := items } ∧
      y'.vStart = y.vStart ∧ y'.topDepth = y.topDepth ∧ y'.firstIdx = y.firstIdx := by
  rcases s with ⟨g, tern, items, sv, sd, nei, fo, ts, tb, tsl⟩
  simp only at h; subst h
  unfold WalkM.maybeUnwrapNxt
  simp only [StateT.run_bind, StateT.run_get, pure_bind]
  split
  · exact ⟨y, _, rfl, rfl, rfl, rfl⟩
  · simp only [StateT.run_bind, run_nxt, run_stackDir, run_getItem, idBind]
    split
    · exact ⟨_, _, rfl, rfl, rfl, rfl⟩
    · exact ⟨y, _, rfl, rfl, rfl, rfl⟩

theorem finishTstackTop_run (i : ItemId) (s : WalkState) {x : TEntry} {rest : List TEntry}
    (h : s.tstack = x :: rest) :
    ∃ x' items, ((finishTstackTop i).run s).2 = { s with tstack := x' :: rest, items := items } ∧
      x'.vStart = x.vStart ∧ x'.topDepth = x.topDepth ∧ x'.firstIdx = x.firstIdx := by
  rcases s with ⟨g, tern, items, sv, sd, nei, fo, ts, tb, tsl⟩
  simp only at h; subst h
  exact ⟨_, _, rfl, rfl, rfl, rfl⟩

theorem mergeTstackTops_run (s : WalkState) {x y : TEntry} {rest : List TEntry}
    (h : s.tstack = x :: y :: rest) :
    ∃ y', (WalkM.mergeTstackTops.run s).2 = { s with tstack := y' :: rest } ∧
      y'.vStart = y.vStart ∧ y'.topDepth ≤ y.topDepth ∧ y'.firstIdx = y.firstIdx := by
  rw [run_mergeTstackTops, h]
  exact ⟨_, rfl, rfl, Nat.min_le_left _ _, rfl⟩

/-! ### Simulation of the blocks -/

variable {bot : List TEntry} {P : WalkState → Prop}

/-- `bot` lies strictly below depth `d` and does not belong to vertex `curV`. -/
def Below (d curV : Nat) (bot : List TEntry) : Prop := ∀ e ∈ bot, e.topDepth < d ∧ e.vStart ≠ curV

theorem Sim.spu {m m' : WalkM Unit} (h : Sim bot P m m' R) :
    Sim bot P m m' fun a a' s' => a = a' ∧ ∃ s, P s ∧ (a', s') = m'.run s :=
  Sim.mono (Sim.sp h) (fun _ => id) fun _ _ _ h => ⟨rfl, h.2⟩

theorem Sim.maybeUnwrapNxt (t : NodeType) (hP : ∀ s, P s → 2 ≤ s.tstack.length) :
    Sim bot P (maybeUnwrapNxt t) (maybeUnwrapNxt t) fun a a' _ => a = a' := by
  unfold WalkM.maybeUnwrapNxt
  refine Sim.bind_get fun a => ?_
  simp only [lift_ternarize]
  refine Sim.iteb _ (fun _ => Sim.allocItem t) fun _ => ?_
  refine Sim.bind_eq (Sim.sp (Sim.nxt fun s h => hP s h.1)) fun x => ?_
  refine Sim.bind_eq (Sim.sp (Sim.stackDir _)) fun dir => ?_
  refine Sim.bind_eq (Sim.sp (Sim.getItem _)) fun it => ?_
  refine Sim.iteb _ (fun _ => ?_) fun _ => ?_
  · refine Sim.seq (Sim.modifyNxt _ ?_) (Sim.pure _ _ fun _ _ => rfl)
    rintro s ⟨s₂, ⟨s₁, ⟨s₀, hs, h0⟩, h1⟩, h2⟩
    simp only [run_stackDir, run_getItem, run_nxt, Prod.mk.injEq] at h0 h1 h2
    rw [h2.2, h1.2, h0.2]; exact hP _ hs.1
  · exact Sim.allocItem t

theorem run_loop1Cond (d : Nat) (s : WalkState) :
    (loop1Cond d).run s = (decide (s.tstack.length ≥ 2) && decide (s.tstack.tail.head!.topDepth ≥ d), s) := rfl

theorem Sim.loop1Cond (d : Nat) (hB : ∀ e ∈ bot, e.topDepth < d) :
    Sim bot P (loop1Cond d) (loop1Cond d) fun b b' s =>
      b = b' ∧ (b' = true → 2 ≤ s.tstack.length ∧ d ≤ s.tstack.tail.head!.topDepth) := by
  rintro ⟨g, tern, items, sv, sd, nei, fo, ts, tb, tsl⟩ hs
  simp only [run_loop1Cond, lift]
  refine ⟨_, rfl, ?_, fun h => by simpa using h⟩
  match ts, bot, hB with
  | [], [], _ | [], [_], _ | [_], [], _ => rfl
  | [], _ :: b :: _, hB => simp [List.head!, Nat.not_le.2 (hB b (by simp))]
  | [_], b :: _, hB => simp [List.head!, Nat.not_le.2 (hB b (by simp))]
  | _ :: _ :: _, _, _ => simp [List.head!]

/-- Loop 1 state: at most `n` entries, bottom entry at depth `≤ d`. -/
def I1 (d n : Nat) (s : WalkState) : Prop := s.tstack.length ≤ n ∧ Inv1 d s.tstack

theorem Sim.loop1Type (d : Nat) (edgeDir : Bool) (n : Nat) :
    Sim bot (fun s => I1 d (n + 1) s ∧ 2 ≤ s.tstack.length ∧ d ≤ s.tstack.tail.head!.topDepth)
      (loop1Type d edgeDir) (loop1Type d edgeDir) fun a a' s => a = a' ∧ 2 ≤ s.tstack.length ∧ I1 d (n + 1) s := by
  unfold Spqr.loop1Type
  refine Sim.bind_eq (Sim.sp (Sim.nxt ?_)) fun t => ?_
  · exact fun s h => h.2.1
  refine Sim.ite _ _ (fun _ _ => Iff.rfl) (fun ht _ => ?_) fun ht _ => ?_
  · refine Sim.bind_eq (Sim.sp (Sim.nxt fun s h => ?_)) fun t₂ => ?_
    · obtain ⟨s₀, hs, h0⟩ := h
      simp only [run_nxt, Prod.mk.injEq] at h0
      rw [h0.2]; exact hs.2.1
    refine Sim.bind_eq (Sim.spu (Sim.setStackDir _ _)) fun _ => ?_
    refine Sim.bind_eq (Sim.spu (Sim.mergeTstackTops fun s h => ?_)) fun _ => ?_
    · obtain ⟨s₂, ⟨s₁, ⟨s₀, hs, h0⟩, h1⟩, h2⟩ := h
      simp only [run_nxt, run_setStackDir, Prod.mk.injEq] at h0 h1 h2
      obtain ⟨x, y, rest, hl⟩ : ∃ x y rest, s₀.tstack = x :: y :: rest := by
        match hl : s₀.tstack, hs.2.1 with
        | x :: y :: l, _ => exact ⟨x, y, l, rfl⟩
      rw [hl] at h0
      simp only [List.tail_cons, List.head!] at h0 ht
      obtain ⟨hty, rfl⟩ := h0
      obtain ⟨-, rfl⟩ := h1
      have hl' : s.tstack = x :: y :: rest := by rw [h2.2]; exact hl
      rw [hl']
      cases rest with
      | nil => exact absurd (hs.1.2 y (by simp [hl])) (Nat.not_le.2 (hty ▸ ht))
      | cons z rest => simp
    refine Sim.pure _ _ fun s h => ⟨rfl, ?_⟩
    obtain ⟨s₃, ⟨s₂, ⟨s₁, ⟨s₀, hs, h0⟩, h1⟩, h2⟩, h3⟩ := h
    simp only [run_nxt, run_setStackDir, Prod.mk.injEq] at h0 h1 h2
    obtain ⟨x, y, rest, hl⟩ : ∃ x y rest, s₀.tstack = x :: y :: rest := by
      match hl : s₀.tstack, hs.2.1 with
      | x :: y :: l, _ => exact ⟨x, y, l, rfl⟩
    rw [hl] at h0
    simp only [List.tail_cons, List.head!] at h0 ht
    obtain ⟨hty, rfl⟩ := h0
    obtain ⟨-, rfl⟩ := h1
    have hl' : s₃.tstack = x :: y :: rest := by rw [h2.2]; exact hl
    have e3 : s = (WalkM.mergeTstackTops.run s₃).2 := congrArg Prod.snd h3
    cases rest with
    | nil => exact absurd (hs.1.2 y (by simp [hl])) (Nat.not_le.2 (hty ▸ ht))
    | cons z rest =>
      obtain ⟨y', hm, -, hy, -⟩ := mergeTstackTops_run s₃ hl'
      rw [e3, hm]
      have hlen := hs.1.1; have hinv := hs.1.2
      rw [hl] at hlen hinv
      refine ⟨by simp, by simp at hlen ⊢; omega, ?_⟩
      simp only
      exact hinv.of_le hy
  · refine Sim.bind_eq (Sim.sp (Sim.nxt fun s h => ?_)) fun a => ?_
    · obtain ⟨s₀, hs, h0⟩ := h
      simp only [run_nxt, Prod.mk.injEq] at h0
      rw [h0.2]; exact hs.2.1
    refine Sim.bind_eq (Sim.sp (Sim.cur fun s h => ?_)) fun b => ?_
    · obtain ⟨s₁, ⟨s₀, hs, h0⟩, h1⟩ := h
      simp only [run_nxt, Prod.mk.injEq] at h0 h1
      rw [h1.2, h0.2]; exact List.ne_nil_of_length_pos (by have := hs.2.1; omega)
    have hQ : ∀ s, (∃ s₁, (∃ s₀, (∃ s', (I1 d (n + 1) s' ∧ 2 ≤ s'.tstack.length ∧
        d ≤ s'.tstack.tail.head!.topDepth) ∧ (t, s₀) = WalkM.nxt.run s') ∧ (a, s₁) = WalkM.nxt.run s₀) ∧
        (b, s) = WalkM.cur.run s₁) → 2 ≤ s.tstack.length ∧ I1 d (n + 1) s := by
      rintro s ⟨s₁, ⟨s₀, ⟨s', hs, h0⟩, h1⟩, h2⟩
      simp only [run_nxt, run_cur, Prod.mk.injEq] at h0 h1 h2
      rw [h2.2, h1.2, h0.2]; exact ⟨hs.2.1, hs.1⟩
    exact Sim.ite _ _ (fun _ _ => Iff.rfl) (fun _ _ => Sim.pure _ _ fun s h => ⟨rfl, hQ s h⟩)
      fun _ _ => Sim.pure _ _ fun s h => ⟨rfl, hQ s h⟩

theorem Sim.loop1Body (d : Nat) (edgeDir : Bool) (n : Nat) :
    Sim bot (fun s => I1 d (n + 1) s ∧ 2 ≤ s.tstack.length ∧ d ≤ s.tstack.tail.head!.topDepth)
      (loop1Body d edgeDir) (loop1Body d edgeDir) fun _ _ s => I1 d n s := by
  unfold Spqr.loop1Body
  refine Sim.bind_eq (Sim.loop1Type d edgeDir n) fun type => ?_
  refine Sim.bind_eq (Sim.sp (Sim.maybeUnwrapNxt type ?_)) fun item => ?_
  · exact fun s h => h.1
  refine Sim.bind_eq (Sim.spu (Sim.mergeTstackTops fun s h => ?_)) fun _ => ?_
  · obtain ⟨s₀, hs, h0⟩ := h
    obtain ⟨x, y, rest, hl⟩ : ∃ x y rest, s₀.tstack = x :: y :: rest := by
      match hl : s₀.tstack, hs.1 with
      | x :: y :: l, _ => exact ⟨x, y, l, rfl⟩
    obtain ⟨y', items, hu, -⟩ := maybeUnwrapNxt_run type s₀ hl
    have : s = ((WalkM.maybeUnwrapNxt type).run s₀).2 := congrArg Prod.snd h0
    rw [this, hu]; simp
  refine Sim.mono (Sim.spu (Sim.finishTstackTop item ?_)) (fun _ h => h) fun _ _ s h => ?_
  · intro s h
    obtain ⟨s₁, ⟨s₀, hs, h0⟩, h1⟩ := h
    obtain ⟨x, y, rest, hl⟩ : ∃ x y rest, s₀.tstack = x :: y :: rest := by
      match hl : s₀.tstack, hs.1 with
      | x :: y :: l, _ => exact ⟨x, y, l, rfl⟩
    obtain ⟨y', items, hu, -⟩ := maybeUnwrapNxt_run type s₀ hl
    have e0 : s₁ = ((WalkM.maybeUnwrapNxt type).run s₀).2 := congrArg Prod.snd h0
    obtain ⟨y'', hm, -⟩ := mergeTstackTops_run s₁ (x := x) (y := y') (rest := rest) (by rw [e0, hu])
    have e1 : s = (WalkM.mergeTstackTops.run s₁).2 := congrArg Prod.snd h1
    rw [e1, hm]; simp
  · obtain ⟨-, s₂, ⟨s₁, ⟨s₀, hs, h0⟩, h1⟩, h2⟩ := h
    obtain ⟨x, y, rest, hl⟩ : ∃ x y rest, s₀.tstack = x :: y :: rest := by
      match hl : s₀.tstack, hs.1 with
      | x :: y :: l, _ => exact ⟨x, y, l, rfl⟩
    obtain ⟨y', items, hu, -, hy', -⟩ := maybeUnwrapNxt_run type s₀ hl
    have e0 : s₁ = ((WalkM.maybeUnwrapNxt type).run s₀).2 := congrArg Prod.snd h0
    obtain ⟨y'', hm, -, hy'', -⟩ := mergeTstackTops_run s₁ (x := x) (y := y') (rest := rest) (by rw [e0, hu])
    have e1 : s₂ = (WalkM.mergeTstackTops.run s₁).2 := congrArg Prod.snd h1
    obtain ⟨y''', items', hf, -, hy''', -⟩ :=
      finishTstackTop_run item s₂ (x := y'') (rest := rest) (by rw [e1, hm])
    have e2 : s = ((WalkM.finishTstackTop item).run s₂).2 := congrArg Prod.snd h2
    rw [e2, hf]
    have hlen := hs.2.1; have hinv := hs.2.2
    rw [hl] at hlen hinv
    refine ⟨by simp at hlen ⊢; omega, ?_⟩
    simp only
    exact (hinv.of_le (y' := y''') (by omega)).cons_of_le (Nat.le_refl _) |>.cons_of_le (Nat.le_refl _)

/-- Loop 1 exit: fewer than two entries, or the entry below the top is above depth `d`. -/
def Q1 (d : Nat) (s : WalkState) : Prop :=
  Inv1 d s.tstack ∧ ¬ (2 ≤ s.tstack.length ∧ d ≤ s.tstack.tail.head!.topDepth)

theorem Sim.loop1 (d : Nat) (edgeDir : Bool) (hB : ∀ e ∈ bot, e.topDepth < d) (n m : Nat)
    (h : n ≤ m) :
    Sim bot (I1 d n) (WalkM.loop m (Spqr.loop1Cond d) (Spqr.loop1Body d edgeDir))
      (WalkM.loop n (Spqr.loop1Cond d) (Spqr.loop1Body d edgeDir)) fun _ _ s => Q1 d s :=
  Sim.loop (I1 d) (fun n s => I1 d n s ∧ 2 ≤ s.tstack.length ∧ d ≤ s.tstack.tail.head!.topDepth)
    (Q1 d) _ _ _ _ (fun _ => rfl)
    (fun n => Sim.mono (Sim.sp (Sim.loop1Cond d hB)) (fun _ => id)
      fun b b' s ⟨⟨hb, hJ⟩, s₀, hs, h0⟩ => by
        simp only [run_loop1Cond, Prod.mk.injEq] at h0
        obtain ⟨rfl, rfl⟩ := h0
        refine ⟨hb, fun ht => ⟨hs, hJ ht⟩, fun hf => ⟨hs.2, fun hc => ?_⟩⟩
        simp at hf; omega)
    (fun s ⟨⟨hl, _⟩, h2, _⟩ => by omega)
    (fun n => Sim.loop1Body d edgeDir n) n m h

theorem Sim.closeEars (nxtV d e : Nat) (edgeDir : Bool) (hB : ∀ e ∈ bot, e.topDepth < d) :
    Sim bot (fun s => Inv1 d s.tstack) (closeEars nxtV d e edgeDir) (closeEars nxtV d e edgeDir)
      fun _ _ s => Q1 d s := by
  unfold Spqr.closeEars
  refine Sim.seq (Q := fun s' => ∃ s : WalkState, Inv1 d s.tstack ∧ ((), s') = (WalkM.pushEdgeTstack nxtV d e).run s)
    (Sim.mono (Sim.sp (Sim.pushEdgeTstack nxtV d e)) (fun _ => id) fun _ _ _ h => h.2) ?_
  refine Sim.bind_tstackSize fun n => ?_
  refine Sim.mono (Sim.loop1 d edgeDir hB n (n + bot.length) (by omega)) ?_ fun _ _ _ h => h
  rintro s ⟨⟨s₀, hs, h0⟩, rfl⟩
  rw [run_pushEdgeTstack, Prod.mk.injEq] at h0
  obtain ⟨-, rfl⟩ := h0
  exact ⟨Nat.le_refl _, hs.push (Nat.le_refl _)⟩

/-! ### Loop 2: merge ears whose first back edge is after `fo` -/

/-- Loop 2 state: at most `n` nonempty entries, bottom entry at depth `≤ d` with `firstIdx ≤ fo`. -/
def I2 (d fo n : Nat) (s : WalkState) : Prop :=
  s.tstack.length ≤ n ∧ s.tstack ≠ [] ∧ Inv1 d s.tstack ∧ Inv2 fo s.tstack

def Q2 (d fo : Nat) (s : WalkState) : Prop :=
  s.tstack ≠ [] ∧ Inv1 d s.tstack ∧ Inv2 fo s.tstack ∧ s.tstack.head!.firstIdx ≤ fo

theorem run_loop2Cond (fo : Nat) (s : WalkState) :
    (loop2Cond fo).run s = (decide (s.tstack.head!.firstIdx > fo), s) := rfl

theorem Sim.loop2Cond (fo : Nat) (hP : ∀ s, P s → s.tstack ≠ []) :
    Sim bot P (loop2Cond fo) (loop2Cond fo) fun b b' s =>
      b = b' ∧ (b' = true → fo < s.tstack.head!.firstIdx) ∧ (b' = false → s.tstack.head!.firstIdx ≤ fo) := by
  rintro ⟨g, tern, items, sv, sd, nei, fo', ts, tb, tsl⟩ hs
  simp only [run_loop2Cond, lift]
  cases ts with
  | nil => exact absurd rfl (hP _ hs)
  | cons x l => exact ⟨_, rfl, rfl, by simp, by simp⟩

theorem Sim.loop2Body (d fo n : Nat) :
    Sim bot (fun s => I2 d fo (n + 1) s ∧ 2 ≤ s.tstack.length) WalkM.mergeTstackTops WalkM.mergeTstackTops
      fun _ _ s => I2 d fo n s := by
  refine Sim.mono (Sim.sp (Sim.mergeTstackTops ?_)) (fun _ => id) ?_
  · exact fun s h => h.2
  rintro _ _ s ⟨-, s₀, hs, h0⟩
  obtain ⟨x, y, rest, hl⟩ : ∃ x y rest, s₀.tstack = x :: y :: rest := by
    match hl : s₀.tstack, hs.2 with
    | x :: y :: l, _ => exact ⟨x, y, l, rfl⟩
  obtain ⟨y', hm, -, hy, hf⟩ := mergeTstackTops_run s₀ hl
  have e : s = (WalkM.mergeTstackTops.run s₀).2 := congrArg Prod.snd h0
  rw [e, hm]
  obtain ⟨⟨hn, -, h1, h2⟩, -⟩ := hs
  rw [hl] at hn h1 h2
  exact ⟨by simp at hn ⊢; omega, by simp, h1.of_le hy, h2.of_eq hf⟩

theorem Sim.loop2 (d fo n m : Nat) (h : n ≤ m) :
    Sim bot (I2 d fo n) (WalkM.loop m (Spqr.loop2Cond fo) WalkM.mergeTstackTops)
      (WalkM.loop n (Spqr.loop2Cond fo) WalkM.mergeTstackTops) fun _ _ s => Q2 d fo s :=
  Sim.loop (I2 d fo) (fun n s => I2 d fo n s ∧ 2 ≤ s.tstack.length) (Q2 d fo) _ _ _ _ (fun _ => rfl)
    (fun n => Sim.mono (Sim.sp (Sim.loop2Cond fo fun s h => h.2.1)) (fun _ => id)
      fun b b' s ⟨⟨hb, ht, hf⟩, s₀, hs, h0⟩ => by
        simp only [run_loop2Cond, Prod.mk.injEq] at h0
        obtain ⟨rfl, rfl⟩ := h0
        refine ⟨hb, fun h => ⟨hs, ?_⟩, fun h => ⟨hs.2.1, hs.2.2.1, hs.2.2.2, hf h⟩⟩
        have := ht h
        match hl : s.tstack, hs.2.1, hs.2.2.2, this with
        | [x], _, h2, h3 => exact absurd (h2 x (by simp)) (by simp [hl, List.head!] at h3; omega)
        | _ :: _ :: _, _, _, _ => simp)
    (fun s ⟨⟨hl, _⟩, h2⟩ => by omega)
    (fun n => Sim.loop2Body d fo n) n m h

theorem Sim.mergeLate (d : Nat) :
    Sim bot (fun s => s.tstack ≠ [] ∧ Inv1 d s.tstack ∧ Inv2 s.firstOccurrence[d]! s.tstack)
      (mergeLate d) (mergeLate d) fun b b' s => b = b' ∧ s.tstack ≠ [] ∧ Inv1 d s.tstack := by
  unfold Spqr.mergeLate
  refine Sim.bind_get fun a => ?_
  simp only [lift_firstOccurrence]
  refine Sim.bind_eq (Sim.sp (Sim.cur fun s h => h.1.1)) fun c => ?_
  refine Sim.ite _ _ (fun _ _ => Iff.rfl) (fun _ _ => ?_) fun _ _ => ?_
  · refine Sim.bind_tstackSize fun n => ?_
    refine Sim.seq (Sim.mono (Sim.loop2 d _ n (n + bot.length) (by omega)) ?_ fun _ _ _ h => h) ?_
    · rintro s ⟨⟨s₀, ⟨hs, rfl⟩, h0⟩, rfl⟩
      simp only [run_cur, Prod.mk.injEq] at h0
      obtain ⟨-, rfl⟩ := h0
      exact ⟨Nat.le_refl _, hs.1, hs.2.1, hs.2.2⟩
    exact Sim.pure _ _ fun s h => ⟨rfl, h.1, h.2.1⟩
  · refine Sim.pure _ _ ?_
    rintro s ⟨s₀, ⟨hs, rfl⟩, h0⟩
    simp only [run_cur, Prod.mk.injEq] at h0
    obtain ⟨-, rfl⟩ := h0
    exact ⟨rfl, hs.1, hs.2.1⟩

/-! ### The P-check and the vertex push -/

theorem run_condP (curV lowval : Nat) (isType1 : Bool) (s : WalkState) :
    (condP curV lowval isType1).run s =
      (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == curV) &&
        (s.tstack.tail.head!.topDepth == lowval), s) := rfl

theorem Sim.condP (curV lowval : Nat) (isType1 : Bool) (hB : ∀ e ∈ bot, e.vStart ≠ curV) :
    Sim bot P (condP curV lowval isType1) (condP curV lowval isType1) fun b b' s =>
      b = b' ∧ (b' = true → 2 ≤ s.tstack.length) := by
  rintro ⟨g, tern, items, sv, sd, nei, fo, ts, tb, tsl⟩ hs
  simp only [run_condP, lift]
  refine ⟨_, rfl, ?_, fun h => by simp at h; exact h.1.1.2⟩
  match ts, bot, hB with
  | [], [], _ | [], [_], _ | [_], [], _ => rfl
  | [], _ :: b :: _, hB => simp [List.head!, hB b (by simp)]
  | [_], b :: _, hB => simp [List.head!, hB b (by simp)]
  | _ :: _ :: _, _, _ => simp [List.head!]

theorem Sim.finishP (curV lowval : Nat) (isType1 : Bool) (hB : ∀ e ∈ bot, e.vStart ≠ curV) :
    Sim bot (fun s => s.tstack ≠ []) (finishP curV lowval isType1) (finishP curV lowval isType1)
      fun _ _ s => s.tstack ≠ [] := by
  unfold Spqr.finishP
  refine Sim.bind_eq (Q := fun c s => (c = true → 2 ≤ s.tstack.length) ∧
      ∃ s₀ : WalkState, s₀.tstack ≠ [] ∧ (c, s) = (Spqr.condP curV lowval isType1).run s₀)
    (Sim.mono (Sim.sp (Sim.condP curV lowval isType1 hB)) (fun _ => id)
      fun _ _ _ h => ⟨h.1.1, h.1.2, h.2⟩) fun c => ?_
  refine Sim.iteb c (fun hc => ?_) fun _ => Sim.pure _ _ ?_
  · have hP : ∀ s, (c = true → 2 ≤ s.tstack.length) ∧
        (∃ s₀ : WalkState, s₀.tstack ≠ [] ∧ (c, s) = (Spqr.condP curV lowval isType1).run s₀) →
        2 ≤ s.tstack.length := fun s h => h.1 hc
    refine Sim.bind_eq (Q := fun _ s => 2 ≤ s.tstack.length)
      (Sim.mono (Sim.sp (Sim.maybeUnwrapNxt .P hP)) (fun _ => id) ?_) fun item => ?_
    · rintro a a' s ⟨rfl, s₀, hs, h0⟩
      obtain ⟨x, y, rest, hl⟩ : ∃ x y rest, s₀.tstack = x :: y :: rest := by
        match hl : s₀.tstack, hP s₀ hs with
        | x :: y :: l, _ => exact ⟨x, y, l, rfl⟩
      obtain ⟨y', items, hu, -⟩ := maybeUnwrapNxt_run .P s₀ hl
      have e : s = ((WalkM.maybeUnwrapNxt .P).run s₀).2 := congrArg Prod.snd h0
      exact ⟨rfl, by rw [e, hu]; simp⟩
    refine Sim.bind (R := fun _ _ s => s.tstack ≠ [])
      (Sim.mono (Sim.sp (Sim.mergeTstackTops fun _ h => h)) (fun _ => id) ?_) fun _ _ => ?_
    · rintro _ _ s ⟨-, s₁, hs, h1⟩
      obtain ⟨x, y, rest, hl⟩ : ∃ x y rest, s₁.tstack = x :: y :: rest := by
        match hl : s₁.tstack, hs with
        | x :: y :: l, _ => exact ⟨x, y, l, rfl⟩
      obtain ⟨y', hm, -⟩ := mergeTstackTops_run s₁ hl
      have e : s = (WalkM.mergeTstackTops.run s₁).2 := congrArg Prod.snd h1
      rw [e, hm]; simp
    refine Sim.mono (Sim.sp (Sim.finishTstackTop item fun _ h => h)) (fun _ => id) ?_
    rintro _ _ s ⟨-, s₂, hs, h2⟩
    obtain ⟨x, l, hl⟩ := List.exists_cons_of_ne_nil hs
    obtain ⟨x', items, hf, -⟩ := finishTstackTop_run item s₂ hl
    have e : s = ((WalkM.finishTstackTop item).run s₂).2 := congrArg Prod.snd h2
    rw [e, hf]; simp
  · rintro s ⟨-, s₀, hs, h0⟩
    simp only [run_condP, Prod.mk.injEq] at h0
    exact h0.2 ▸ hs

theorem Sim.finishTail (curV d : Nat) (hasVert isSingle : Bool) :
    Sim bot (fun s => s.tstack ≠ []) (finishTail curV d hasVert isSingle)
      (finishTail curV d hasVert isSingle) fun a a' _ => a = a' := by
  rintro ⟨g, tern, items, sv, sd, nei, fo, ts, tb, tsl⟩ hs
  obtain ⟨x, l, rfl⟩ := List.exists_cons_of_ne_nil hs
  refine ⟨_, ?_, rfl⟩
  cases hasVert <;> cases isSingle <;> rfl

theorem Sim.finishRest (curV d lowval : Nat) (isType1 hasVert isSingle : Bool)
    (hB : ∀ e ∈ bot, e.vStart ≠ curV) :
    Sim bot (fun s => s.tstack ≠ []) (finishRest curV d lowval isType1 hasVert isSingle)
      (finishRest curV d lowval isType1 hasVert isSingle) fun a a' _ => a = a' :=
  Sim.seq (Sim.finishP curV lowval isType1 hB) (Sim.finishTail curV d hasVert isSingle)

/-! ### The block boundary -/

theorem Sim.finishBoundary (curV d : Nat) (o : DfsOut) (qItem : ItemId) (hasVert : Bool) :
    Sim bot (fun s => o.cls.isTree = true →
        if o.cls.lowval d == d + 1 then s.tstack ≠ [] else 2 ≤ s.tstack.length)
      (finishBoundary curV d o qItem hasVert) (finishBoundary curV d o qItem hasVert)
      fun a a' _ => a = a' := by
  rintro ⟨g, tern, items, sv, sd, nei, fo, ts, tb, tsl⟩ hs
  refine ⟨_, ?_, rfl⟩
  unfold Spqr.finishBoundary
  cases hT : o.cls.isTree
  · rfl
  · have hs := hs hT
    cases hL : (o.cls.lowval d == d + 1)
    · simp only [hL, Bool.false_eq_true, ite_false] at hs
      match ts, hs with
      | _ :: _ :: _, _ => rfl
    · simp only [hL, ite_true] at hs
      match ts, hs with
      | _ :: _, _ => rfl

/-! ### Remaining obligations -/

/-- The `hasVert` block. Not yet closed: `loop3Cond` is already synchronised by the `+ bot.length`
offset of `origTstack`; what remains is the join-point/`let mut` structure and `Sim.map`. -/
theorem Sim.closeVert (curV : Nat) (edgeDir isType1 : Bool) (origTstack : Nat) (isSingle : Bool)
    (k : Bool → WalkM Bool)
    (hk : ∀ b, Sim bot (fun s => s.tstack ≠ []) (k b) (k b) fun a a' _ => a = a') :
    Sim bot (fun s => 3 ≤ s.tstack.length ∧ (isType1 = false → origTstack + 3 ≤ s.tstack.length))
      (closeVert curV edgeDir isType1 (origTstack + bot.length) isSingle k)
      (closeVert curV edgeDir isType1 origTstack isSingle k) fun a a' _ => a = a' := by
  sorry

/-- The walk of a subtree never inspects entries strictly below its depth. This is the walk
invariant: it supplies the block preconditions (`Inv1`, `Inv2`, the `≥ 3` shapes) at each call of
`finishEdge`, which the block lemmas above then discharge. -/
theorem Sim.walkTree (t : DfsTree) (d : Nat)
    (hB : ∀ e ∈ bot, e.topDepth < d ∧ e.vStart ∉ t.verts) :
    Sim bot (fun s => s.tstack = [] ∧ ∀ e ∈ bot, e.firstIdx ≤ s.nxtEdgeIdx)
      (walkTree t d) (walkTree t d) fun _ _ _ => True := by
  sorry

end Spqr
