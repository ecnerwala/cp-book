import Spqr.Sim
import Spqr.WalkTyping
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

/-- Wrap `cur` into the pending vertex ear, then continue with `k`. -/
def closeVertTail (curV : Nat) (edgeDir isSingle : Bool) (k : Bool → WalkM Bool)
    (item : Option ItemId) : WalkM Bool := do
  mergeTstackTops
  mergeTstackTops
  modifyCur fun t => { t with vStart := curV, spans := setSides (!edgeDir) (t.spans.1 ++ t.spans.2) [] }
  match item with
  | some item => finishTstackTop item; k true
  | none => k isSingle

/-- The `hasVert` case: close the sub-ear down to three entries (type 2), or unwrap its top (type 1). -/
def closeVert (curV : Nat) (edgeDir isType1 : Bool) (origTstack : Nat) (isSingle : Bool)
    (k : Bool → WalkM Bool) : WalkM Bool := do
  if !isType1 then
    loop (← tstackSize) (loop3Cond origTstack) mergeTstackTops
    closeVertTail curV edgeDir false k none
  else
    let item ← some <$> maybeUnwrapNxt (if isSingle then .S else .R)
    closeVertTail curV edgeDir isSingle k item

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

/-- Record the endpoints of the returning edge; returns the edge direction. -/
def finishSetup (d : Nat) (o : DfsOut) : WalkM Bool := do
  let g := (← get).g
  let edgeDir ← stackDir d
  let vs ← makeVs o.dest d
  modifyItem (edgeItem g o.e) fun it => { it with vs := vs }
  pure edgeDir

/-- A returning tree edge: close the ears of the subtree, then the vertex ear. -/
def finishTree (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert edgeDir : Bool) :
    WalkM Bool := do
  closeEars o.dest d o.e edgeDir
  let isSingle ← mergeLate d
  if hasVert then
    closeVert curV edgeDir o.cls.isType1 origTstack isSingle
      (finishRest curV d (o.cls.lowval d) o.cls.isType1 hasVert)
  else finishRest curV d (o.cls.lowval d) o.cls.isType1 hasVert isSingle

/-- A returning back edge: push it as a one-edge ear. -/
def finishBack (curV d : Nat) (o : DfsOut) (hasVert : Bool) : WalkM Bool := do
  pushEdgeTstack curV (o.cls.lowval d) o.e
  modify fun s =>
    { s with firstOccurrence := s.firstOccurrence.modify (o.cls.lowval d) (min · s.nxtEdgeIdx),
             nxtEdgeIdx := s.nxtEdgeIdx + 1 }
  finishRest curV d (o.cls.lowval d) o.cls.isType1 hasVert true

def finishEdge' (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool) : WalkM Bool := do
  let g := (← get).g
  let edgeDir ← stackDir d
  if o.cls.lowval d ≥ d then finishBoundary curV d o (edgeItem g o.e) hasVert
  else
    let vs ← makeVs o.dest d
    modifyItem (edgeItem g o.e) fun it => { it with vs := vs }
    if o.cls.isTree then finishTree curV d o origTstack hasVert edgeDir
    else finishBack curV d o hasVert

theorem ite_bind (c : Prop) [Decidable c] (a b : WalkM α) (f : α → WalkM β) :
    (if c then a else b) >>= f = if c then a >>= f else b >>= f := by split <;> rfl

theorem finishEdge_eq (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool) :
    finishEdge curV d o origTstack hasVert = finishEdge' curV d o origTstack hasVert := by
  simp only [finishEdge, finishEdge', finishBoundary, closeEars, loop1Cond, loop1Type, loop1Body, mergeLate,
    loop2Cond, loop3Cond, closeVert, closeVertTail, finishRest, finishTree, finishBack, condP, finishP, finishTail, bind_assoc, Bool.false_eq_true, Bool.not_false,
    Bool.not_true, ite_true, ite_false, pure_bind, ite_bind]
  generalize o.cls.isType1 = b
  cases b <;> rfl

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
def I1 (d n : Nat) (s : WalkState) : Prop := s.tstack.length ≤ n ∧ s.tstack ≠ [] ∧ Inv1 d s.tstack

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
      | nil => exact absurd (hs.1.2.2 y (by simp [hl])) (Nat.not_le.2 (hty ▸ ht))
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
    | nil => exact absurd (hs.1.2.2 y (by simp [hl])) (Nat.not_le.2 (hty ▸ ht))
    | cons z rest =>
      obtain ⟨y', hm, -, hy, -⟩ := mergeTstackTops_run s₃ hl'
      rw [e3, hm]
      have hlen := hs.1.1; have hinv := hs.1.2.2
      rw [hl] at hlen hinv
      refine ⟨by simp, by simp at hlen ⊢; omega, by simp, ?_⟩
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
    have hlen := hs.2.1; have hinv := hs.2.2.2
    rw [hl] at hlen hinv
    refine ⟨by simp at hlen ⊢; omega, by simp, ?_⟩
    simp only
    exact (hinv.of_le (y' := y''') (by omega)).cons_of_le (Nat.le_refl _) |>.cons_of_le (Nat.le_refl _)

/-- Loop 1 exit: fewer than two entries, or the entry below the top is above depth `d`. -/
def Q1 (d : Nat) (s : WalkState) : Prop :=
  s.tstack ≠ [] ∧ Inv1 d s.tstack ∧ ¬ (2 ≤ s.tstack.length ∧ d ≤ s.tstack.tail.head!.topDepth)

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
        refine ⟨hb, fun ht => ⟨hs, hJ ht⟩, fun hf => ⟨hs.2.1, hs.2.2, fun hc => ?_⟩⟩
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
  exact ⟨Nat.le_refl _, by simp, hs.push (Nat.le_refl _)⟩

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
        | [x], _, h2, h3 => exact absurd (h2 x (by simp)) (by simp [List.head!] at h3; omega)
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

/-! ### Length-tracking primitives -/

theorem Sim.mergeTstackTops_len (Pl : Nat → Prop) :
    Sim bot (fun s => 2 ≤ s.tstack.length ∧ Pl s.tstack.length) WalkM.mergeTstackTops
      WalkM.mergeTstackTops fun _ _ s => Pl (s.tstack.length + 1) := by
  refine Sim.mono (Sim.sp (Sim.mergeTstackTops fun _ h => h.1)) (fun _ => id) ?_
  rintro _ _ s ⟨-, s₀, hs, h0⟩
  obtain ⟨x, y, rest, hl⟩ : ∃ x y rest, s₀.tstack = x :: y :: rest := by
    match hl : s₀.tstack, hs.1 with
    | x :: y :: l, _ => exact ⟨x, y, l, rfl⟩
  obtain ⟨y', hm, -⟩ := mergeTstackTops_run s₀ hl
  have e : s = (WalkM.mergeTstackTops.run s₀).2 := congrArg Prod.snd h0
  rw [e, hm]; have := hs.2; rw [hl] at this; simpa using this

theorem Sim.maybeUnwrapNxt_len (t : NodeType) (Pl : Nat → Prop) :
    Sim bot (fun s => 2 ≤ s.tstack.length ∧ Pl s.tstack.length) (WalkM.maybeUnwrapNxt t)
      (WalkM.maybeUnwrapNxt t) fun a a' s => a = a' ∧ Pl s.tstack.length := by
  refine Sim.mono (Sim.sp (Sim.maybeUnwrapNxt t fun _ h => h.1)) (fun _ => id) ?_
  rintro _ _ s ⟨he, s₀, hs, h0⟩
  obtain ⟨x, y, rest, hl⟩ : ∃ x y rest, s₀.tstack = x :: y :: rest := by
    match hl : s₀.tstack, hs.1 with
    | x :: y :: l, _ => exact ⟨x, y, l, rfl⟩
  obtain ⟨y', items, hu, -⟩ := maybeUnwrapNxt_run t s₀ hl
  have e : s = ((WalkM.maybeUnwrapNxt t).run s₀).2 := congrArg Prod.snd h0
  refine ⟨he, ?_⟩
  rw [e, hu]; have := hs.2; rw [hl] at this; simpa using this

theorem Sim.modifyCur_ne (f : TEntry → TEntry) :
    Sim bot (fun s => s.tstack ≠ []) (WalkM.modifyCur f) (WalkM.modifyCur f)
      fun _ _ s => s.tstack ≠ [] := by
  refine Sim.mono (Sim.sp (Sim.modifyCur f fun _ h => h)) (fun _ => id) ?_
  rintro _ _ s ⟨-, s₀, hs, h0⟩
  obtain ⟨x, l, hl⟩ := List.exists_cons_of_ne_nil hs
  have e : s = ((WalkM.modifyCur f).run s₀).2 := congrArg Prod.snd h0
  rw [e, run_modifyCur, hl]; simp

theorem Sim.finishTstackTop_ne (i : ItemId) :
    Sim bot (fun s => s.tstack ≠ []) (WalkM.finishTstackTop i) (WalkM.finishTstackTop i)
      fun _ _ s => s.tstack ≠ [] := by
  refine Sim.mono (Sim.sp (Sim.finishTstackTop i fun _ h => h)) (fun _ => id) ?_
  rintro _ _ s ⟨-, s₀, hs, h0⟩
  obtain ⟨x, l, hl⟩ := List.exists_cons_of_ne_nil hs
  obtain ⟨x', items, hf, -⟩ := finishTstackTop_run i s₀ hl
  have e : s = ((WalkM.finishTstackTop i).run s₀).2 := congrArg Prod.snd h0
  rw [e, hf]; simp

/-! ### Loop 3: merge down to the three entries of the sub-ear -/

theorem run_loop3Cond (o : Nat) (s : WalkState) :
    (loop3Cond o).run s = (decide (s.tstack.length > o + 3), s) := rfl

theorem Sim.loop3Cond (o : Nat) :
    Sim bot P (loop3Cond (o + bot.length)) (loop3Cond o) fun b b' s =>
      b = b' ∧ (b' = true → o + 3 < s.tstack.length) ∧ (b' = false → s.tstack.length ≤ o + 3) := by
  intro s _
  simp only [run_loop3Cond, lift_tstack, List.length_append]
  refine ⟨_, rfl, ?_, fun h => by simpa using h, fun h => by simpa using h⟩
  simp only [decide_eq_decide]; omega

def I3 (o n : Nat) (s : WalkState) : Prop := s.tstack.length ≤ n ∧ o + 3 ≤ s.tstack.length

theorem Sim.loop3Body (o n : Nat) :
    Sim bot (fun s => I3 o (n + 1) s ∧ o + 3 < s.tstack.length) WalkM.mergeTstackTops
      WalkM.mergeTstackTops fun _ _ s => I3 o n s :=
  Sim.mono (Sim.mergeTstackTops_len fun l => l ≤ n + 1 ∧ o + 3 < l)
    (fun s ⟨⟨h1, _⟩, h2⟩ => ⟨by omega, h1, h2⟩) fun _ _ _ ⟨h1, h2⟩ => ⟨by omega, by omega⟩

theorem Sim.loop3 (o n m : Nat) (h : n ≤ m) :
    Sim bot (I3 o n) (WalkM.loop m (Spqr.loop3Cond (o + bot.length)) WalkM.mergeTstackTops)
      (WalkM.loop n (Spqr.loop3Cond o) WalkM.mergeTstackTops) fun _ _ s => s.tstack.length = o + 3 :=
  Sim.loop (I3 o) (fun n s => I3 o n s ∧ o + 3 < s.tstack.length) (fun s => s.tstack.length = o + 3)
    _ _ _ _ (fun _ => rfl)
    (fun n => Sim.mono (Sim.sp (Sim.loop3Cond o)) (fun _ => id)
      fun b b' s ⟨⟨hb, ht, hf⟩, s₀, hs, h0⟩ => by
        simp only [run_loop3Cond, Prod.mk.injEq] at h0
        obtain ⟨rfl, rfl⟩ := h0
        exact ⟨hb, fun h => ⟨hs, ht h⟩, fun h => by have := hf h; have := hs.2; omega⟩)
    (fun s ⟨⟨hl, _⟩, h2⟩ => by omega)
    (fun n => Sim.loop3Body o n) n m h

theorem Sim.closeVertTail (curV : Nat) (edgeDir isSingle : Bool) (k : Bool → WalkM Bool)
    (hk : ∀ b, Sim bot (fun s => s.tstack ≠ []) (k b) (k b) fun a a' _ => a = a') (item : Option ItemId) :
    Sim bot (fun s => 3 ≤ s.tstack.length) (closeVertTail curV edgeDir isSingle k item)
      (closeVertTail curV edgeDir isSingle k item) fun a a' _ => a = a' := by
  unfold Spqr.closeVertTail
  refine Sim.seq (Sim.mono (Sim.mergeTstackTops_len fun l => 3 ≤ l)
    (fun s h => ⟨by omega, h⟩) fun _ _ _ h => h) ?_
  refine Sim.seq (Sim.mono (Sim.mergeTstackTops_len fun l => 2 ≤ l)
    (fun s h => ⟨by omega, by omega⟩) fun _ _ _ h => h) ?_
  refine Sim.seq (Sim.mono (Sim.modifyCur_ne _)
    (fun s h => List.ne_nil_of_length_pos (by omega)) fun _ _ _ h => h) ?_
  cases item with
  | some i => exact Sim.seq (Sim.finishTstackTop_ne i) (hk true)
  | none => exact hk isSingle

/-- The `hasVert` block. -/
theorem Sim.closeVert (curV : Nat) (edgeDir isType1 : Bool) (origTstack : Nat) (isSingle : Bool)
    (k : Bool → WalkM Bool)
    (hk : ∀ b, Sim bot (fun s => s.tstack ≠ []) (k b) (k b) fun a a' _ => a = a') :
    Sim bot (fun s => 3 ≤ s.tstack.length ∧ (isType1 = false → origTstack + 3 ≤ s.tstack.length))
      (closeVert curV edgeDir isType1 (origTstack + bot.length) isSingle k)
      (closeVert curV edgeDir isType1 origTstack isSingle k) fun a a' _ => a = a' := by
  unfold Spqr.closeVert
  cases isType1
  · simp only [Bool.not_false, ite_true]
    refine Sim.bind_tstackSize fun n => ?_
    refine Sim.seq (Sim.mono (Sim.loop3 origTstack n (n + bot.length) (by omega))
      (fun s ⟨⟨_, h⟩, hn⟩ => ⟨by omega, h trivial⟩) fun _ _ _ h => h) ?_
    exact Sim.mono (Sim.closeVertTail curV edgeDir false k hk none) (fun s h => by omega)
      fun _ _ _ h => h
  · simp only [Bool.not_true, Bool.false_eq_true, ite_false]
    refine Sim.bind_eq (Q := fun _ s => 3 ≤ s.tstack.length)
      (Sim.mono (Sim.map some (Sim.maybeUnwrapNxt_len _ fun l => 3 ≤ l))
        (fun s h => ⟨by omega, h.1⟩) fun a a' s ⟨x, x', hx, hx', he, hl⟩ => ⟨hx ▸ hx' ▸ he ▸ rfl, hl⟩)
      fun item => Sim.closeVertTail curV edgeDir isSingle k hk item

/-! ### `finishEdge` -/

/-- The tstack-shape guards `finishEdge` relies on, as facts about the unlifted run: the boundary
pops find their entries; the bottom entry is at depth `≤ d` (loop 1); after loop 1 the bottom entry's
first back edge is not after `firstOccurrence[d]` (loop 2); after loop 2 the vertex ear has its
three entries (loop 3). -/
def FinishGuards (d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool)
    (s : WalkState) : Prop :=
  (d ≤ o.cls.lowval d → o.cls.isTree = true →
    if o.cls.lowval d == d + 1 then s.tstack ≠ [] else 2 ≤ s.tstack.length) ∧
  (o.cls.lowval d < d → o.cls.isTree = true →
    wp (finishSetup d o) (fun edgeDir s₁ => Inv1 d s₁.tstack ∧
      wp (closeEars o.dest d o.e edgeDir) (fun _ s₂ =>
        Inv2 s₂.firstOccurrence[d]! s₂.tstack ∧
        (hasVert = true → wp (mergeLate d) (fun _ s₃ => 3 ≤ s₃.tstack.length ∧
          (o.cls.isType1 = false → origTstack + 3 ≤ s₃.tstack.length)) s₂)) s₁) s)

theorem Sim.finishTree (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert edgeDir : Bool)
    (hB : Below d curV bot) :
    Sim bot (fun s => Inv1 d s.tstack ∧
        wp (Spqr.closeEars o.dest d o.e edgeDir) (fun _ s₂ =>
          Inv2 s₂.firstOccurrence[d]! s₂.tstack ∧
          (hasVert = true → wp (Spqr.mergeLate d) (fun _ s₃ => 3 ≤ s₃.tstack.length ∧
            (o.cls.isType1 = false → origTstack + 3 ≤ s₃.tstack.length)) s₂)) s)
      (finishTree curV d o (origTstack + bot.length) hasVert edgeDir)
      (finishTree curV d o origTstack hasVert edgeDir) fun a a' _ => a = a' := by
  unfold Spqr.finishTree
  refine Sim.seq (Q := fun s₂ => Q1 d s₂ ∧ Inv2 s₂.firstOccurrence[d]! s₂.tstack ∧
      (hasVert = true → wp (Spqr.mergeLate d) (fun _ s₃ => 3 ≤ s₃.tstack.length ∧
        (o.cls.isType1 = false → origTstack + 3 ≤ s₃.tstack.length)) s₂))
    (Sim.mono (Sim.sp (Sim.mono (Sim.closeEars _ _ _ _ fun e he => (hB e he).1) (fun s h => h.1)
      fun _ _ _ h => h)) (fun _ => id) ?_) ?_
  · rintro _ _ s₂ ⟨hq, s₁, ⟨-, hw⟩, h1⟩
    have e : s₂ = ((Spqr.closeEars o.dest d o.e edgeDir).run s₁).2 := congrArg Prod.snd h1
    subst e
    exact ⟨hq, hw⟩
  refine Sim.bind_eq (Q := fun _ s₃ => s₃.tstack ≠ [] ∧ (hasVert = true → 3 ≤ s₃.tstack.length ∧
      (o.cls.isType1 = false → origTstack + 3 ≤ s₃.tstack.length)))
    (Sim.mono (Sim.sp (Sim.mono (Sim.mergeLate d) (fun s h => ⟨h.1.1, h.1.2.1, h.2.1⟩)
      fun _ _ _ h => h)) (fun _ => id) ?_) fun isSingle => ?_
  · rintro b b' s₃ ⟨⟨hb, hne, -⟩, s₂, ⟨-, -, hw⟩, h2⟩
    have e : s₃ = ((Spqr.mergeLate d).run s₂).2 := congrArg Prod.snd h2
    subst e
    exact ⟨hb, hne, hw⟩
  refine Sim.iteb _ (fun hv => ?_) fun hv => ?_
  · exact Sim.mono (Sim.closeVert curV edgeDir _ origTstack isSingle _ fun b =>
      Sim.finishRest curV d _ _ hasVert b fun e he => (hB e he).2) (fun s h => h.2 hv)
      fun _ _ _ h => h
  · exact Sim.mono (Sim.finishRest curV d _ _ hasVert isSingle fun e he => (hB e he).2)
      (fun s h => h.1) fun _ _ _ h => h

theorem Sim.finishBack (curV d : Nat) (o : DfsOut) (hasVert : Bool) (hB : Below d curV bot) :
    Sim bot (fun _ => True) (finishBack curV d o hasVert) (finishBack curV d o hasVert)
      fun a a' _ => a = a' := by
  unfold Spqr.finishBack
  refine Sim.seq (Q := fun s => s.tstack ≠ [])
    (Sim.mono (Sim.sp (Sim.pushEdgeTstack _ _ _)) (fun _ => id) ?_) ?_
  · rintro _ _ s ⟨-, s₀, -, h0⟩
    rw [run_pushEdgeTstack, Prod.mk.injEq] at h0
    obtain ⟨-, rfl⟩ := h0
    simp
  refine Sim.seq (Q := fun s => s.tstack ≠ [])
    (Sim.mono (Sim.sp (Sim.modify_of _ fun _ => rfl)) (fun _ => id) ?_)
    (Sim.finishRest curV d _ _ hasVert true fun e he => (hB e he).2)
  rintro _ _ s ⟨-, s₀, hs, h0⟩
  have e : s = ((modify (fun s : WalkState =>
    { s with firstOccurrence := s.firstOccurrence.modify (o.cls.lowval d) (min · s.nxtEdgeIdx),
             nxtEdgeIdx := s.nxtEdgeIdx + 1 }) : WalkM Unit).run s₀).2 := congrArg Prod.snd h0
  subst e
  exact hs

theorem Sim.finishEdge' (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool)
    (hB : Below d curV bot) :
    Sim bot (FinishGuards d o origTstack hasVert)
      (finishEdge' curV d o (origTstack + bot.length) hasVert)
      (finishEdge' curV d o origTstack hasVert) fun a a' _ => a = a' := by
  unfold Spqr.finishEdge'
  refine Sim.bind_get fun s₀ => ?_
  simp only [lift_g]
  refine Sim.bind_eq (Sim.sp (Sim.stackDir d)) fun edgeDir => ?_
  refine Sim.ite _ _ (fun _ _ => Iff.rfl) (fun hc _ => ?_) fun hc _ => ?_
  · refine Sim.mono (Sim.finishBoundary curV d o _ hasVert) ?_ fun _ _ _ h => h
    rintro s ⟨s₁, ⟨hg, rfl⟩, h1⟩
    rw [run_stackDir, Prod.mk.injEq] at h1
    obtain ⟨-, rfl⟩ := h1
    exact hg.1 hc
  · refine Sim.bind_eq (Sim.sp (Sim.makeVs _ _)) fun vs => ?_
    refine Sim.bind_eq (Sim.mono (Sim.sp (Sim.modifyItem _ _)) (fun _ => id) fun _ _ _ h => ⟨rfl, h.2⟩)
      fun _ => ?_
    refine Sim.iteb _ (fun ht => ?_) fun ht => ?_
    · refine Sim.mono (Sim.finishTree curV d o origTstack hasVert edgeDir hB) ?_ fun _ _ _ h => h
      rintro s ⟨s₃, ⟨s₂, ⟨s₁, ⟨hg, rfl⟩, h1⟩, h2⟩, h3⟩
      rw [run_stackDir, Prod.mk.injEq] at h1
      rw [run_makeVs, Prod.mk.injEq] at h2
      rw [run_modifyItem, Prod.mk.injEq] at h3
      obtain ⟨rfl, rfl⟩ := h1
      obtain ⟨rfl, rfl⟩ := h2
      obtain ⟨-, rfl⟩ := h3
      exact hg.2 (Nat.lt_of_not_le hc) ht
    · exact Sim.mono (Sim.finishBack curV d o hasVert hB) (fun _ _ => trivial) fun _ _ _ h => h

theorem Sim.finishEdge (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool)
    (hB : Below d curV bot) :
    Sim bot (FinishGuards d o origTstack hasVert)
      (finishEdge curV d o (origTstack + bot.length) hasVert)
      (finishEdge curV d o origTstack hasVert) fun a a' _ => a = a' := by
  rw [finishEdge_eq, finishEdge_eq]
  exact Sim.finishEdge' curV d o origTstack hasVert hB

/-! ### The walk -/

/-- The prelude of `walkOut`: set the direction of depth `d` and push the vertex entry if needed. -/
def walkOutPre (v d : Nat) (o : DfsOut) (hasVert : Bool) : WalkM Bool := do
  let lowval := o.cls.lowval d
  let lowDir ← stackDir lowval
  setStackDir d (if lowval ≥ d then false else !lowDir)
  if !hasVert && lowval < d && o.cls.isType1 then
    pushVertTstack v d
    pure true
  else pure hasVert

def walkOutRest (v d : Nat) (o : DfsOut) (hasVert : Bool) : WalkM Bool := do
  let origTstack ← tstackSize
  match o with
  | .tree _ _ child =>
    modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }
    walkTree child (d + 1)
  | .back .. => pure ()
  finishEdge v d o origTstack hasVert

theorem walkOut_eq (v d : Nat) (o : DfsOut) (hasVert : Bool) :
    walkOut v d o hasVert = walkOutPre v d o hasVert >>= walkOutRest v d o := by
  unfold walkOut walkOutPre walkOutRest
  simp only [bind_assoc, ite_bind]
  rfl

theorem Sim.walkOutPre (v d : Nat) (o : DfsOut) (hasVert : Bool) :
    Sim bot P (walkOutPre v d o hasVert) (walkOutPre v d o hasVert) fun a a' _ => a = a' := by
  unfold Spqr.walkOutPre
  refine Sim.bind_eq (Q := fun _ _ => True)
    (Sim.mono (Sim.stackDir _) (fun _ => id) fun _ _ _ h => ⟨h, trivial⟩) fun lowDir => ?_
  refine Sim.seq (Q := fun _ => True) (Sim.setStackDir _ _) ?_
  exact Sim.iteb _ (fun _ => Sim.seq (Q := fun _ => True) (Sim.pushVertTstack _ _)
    (Sim.pure _ _ fun _ _ => rfl)) fun _ => Sim.pure _ _ fun _ _ => rfl

mutual
/-- The tstack guards (`FinishGuards`) hold at every `finishEdge` call of the walk of `t` from `s`. -/
def GuardsTree (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  match t with
  | .node v outs => GuardsOuts v d outs false { s with stackVerts := s.stackVerts.set! d v }

def GuardsOuts (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  match outs with
  | [] => True
  | o :: rest => GuardsOut v d o hasVert s ∧
      wp (walkOut v d o hasVert) (fun hasVert' s' => GuardsOuts v d rest hasVert' s') s

def GuardsOut (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  wp (walkOutPre v d o hasVert) (fun hasVert' s₁ =>
    match o with
    | .tree _ _ child =>
      wp (modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }) (fun _ s₂ =>
        GuardsTree child (d + 1) s₂ ∧
        wp (walkTree child (d + 1)) (fun _ s₃ => FinishGuards d o s₁.tstack.length hasVert' s₃) s₂) s₁
    | .back .. => FinishGuards d o s₁.tstack.length hasVert' s₁) s
end

theorem DfsOut.vertsList_cons (o : DfsOut) (rest : List DfsOut) :
    DfsOut.vertsList (o :: rest) = DfsOut.vertsList [o] ++ DfsOut.vertsList rest := by
  cases o <;> simp [DfsOut.vertsList]

theorem Sim.walk_aux :
    (∀ (t : DfsTree) (d : Nat) (bot : List TEntry),
      (∀ e ∈ bot, e.topDepth < d ∧ e.vStart ∉ t.verts) →
      Sim bot (GuardsTree t d) (walkTree t d) (walkTree t d) fun _ _ _ => True) ∧
    (∀ (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (bot : List TEntry),
      (∀ e ∈ bot, e.topDepth < d ∧ e.vStart ≠ v ∧ e.vStart ∉ DfsOut.vertsList outs) →
      Sim bot (GuardsOuts v d outs hasVert) (walkOuts v d outs hasVert) (walkOuts v d outs hasVert)
        fun a a' _ => a = a') ∧
    (∀ (v d : Nat) (o : DfsOut) (hasVert : Bool) (bot : List TEntry),
      (∀ e ∈ bot, e.topDepth < d ∧ e.vStart ≠ v ∧ e.vStart ∉ DfsOut.vertsList [o]) →
      Sim bot (GuardsOut v d o hasVert) (walkOut v d o hasVert) (walkOut v d o hasVert)
        fun a a' _ => a = a') := by
  refine walkTree.mutual_induct _ _ _ ?_ ?_ ?_ ?_
  · intro d v outs ih bot hB
    simp only [walkTree]
    refine Sim.seq (Q := GuardsOuts v d outs false)
      (Sim.mono (Sim.sp (Sim.modify_of _ fun _ => rfl)) (fun _ => id) ?_) ?_
    · rintro _ _ s ⟨-, s₀, hs, h0⟩
      have e : s = ((modify fun s : WalkState =>
        { s with stackVerts := s.stackVerts.set! d v } : WalkM Unit).run s₀).2 := congrArg Prod.snd h0
      subst e
      exact hs
    refine Sim.bind_eq (Q := fun _ _ => True)
      (Sim.mono (ih bot ?_) (fun _ => id) fun _ _ _ h => ⟨h, trivial⟩) fun hv => ?_
    · intro e he
      obtain ⟨h1, h2⟩ := hB e he
      simp only [DfsTree.verts, List.mem_cons, not_or] at h2
      exact ⟨h1, h2.1, h2.2⟩
    exact Sim.iteb _ (fun _ => Sim.pure _ _ fun _ _ => trivial)
      fun _ => Sim.seq (Q := fun _ => True) (Sim.setStackDir _ _) (Sim.pushVertTstack _ _)
  · intro o v d hasVert ih bot hB
    rw [walkOut_eq]
    refine Sim.bind_eq (Q := fun hv' s₁ => ∃ s, GuardsOut v d o hasVert s ∧
        (hv', s₁) = (Spqr.walkOutPre v d o hasVert).run s)
      (Sim.mono (Sim.sp (Sim.walkOutPre v d o hasVert)) (fun _ => id) fun _ _ _ h => h) fun hv' => ?_
    unfold walkOutRest
    refine Sim.bind_tstackSize fun n => ?_
    have hB' : Below d v bot := fun e he => ⟨(hB e he).1, (hB e he).2.1⟩
    cases o with
    | back e dest cls =>
      simp only
      refine Sim.mono (Sim.finishEdge v d _ n hv' hB') ?_ fun _ _ _ h => h
      rintro s ⟨⟨s₀, hg, h0⟩, rfl⟩
      unfold GuardsOut at hg
      rw [WalkM.wp, ← h0] at hg
      exact hg
    | tree e cls child =>
      dsimp only at ih
      simp only
      refine Sim.seq (Q := fun s₂ => GuardsTree child (d + 1) s₂ ∧
          wp (walkTree child (d + 1)) (fun _ s₃ => FinishGuards d (.tree e cls child) n hv' s₃) s₂)
        (Sim.mono (Sim.sp (Sim.modify_of _ fun _ => rfl)) (fun _ => id) ?_) ?_
      · rintro _ _ s₂ ⟨-, s₁, ⟨⟨s₀, hg, h0⟩, rfl⟩, h1⟩
        have e2 : s₂ = ((modify fun s : WalkState =>
          { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne } : WalkM Unit).run s₁).2 :=
          congrArg Prod.snd h1
        subst e2
        unfold GuardsOut at hg
        rw [WalkM.wp, ← h0] at hg
        exact hg
      refine Sim.seq (Q := fun s₃ => FinishGuards d (.tree e cls child) n hv' s₃)
        (Sim.mono (Sim.sp (Sim.mono (ih bot ?_) (fun s h => h.1) fun _ _ _ h => h)) (fun _ => id) ?_)
        (Sim.finishEdge v d _ n hv' hB')
      · intro e' he'
        obtain ⟨h1, -, h2⟩ := hB e' he'
        refine ⟨by omega, ?_⟩
        simpa [DfsOut.vertsList] using h2
      · rintro _ _ s₃ ⟨-, s₂, ⟨-, hw⟩, h2⟩
        have e3 : s₃ = ((walkTree child (d + 1)).run s₂).2 := congrArg Prod.snd h2
        subst e3
        exact hw
  · intro v d hasVert bot hB
    simp only [walkOuts]
    exact Sim.pure _ _ fun _ _ => rfl
  · intro v d hasVert o rest ih₁ ih₂ bot hB
    simp only [walkOuts]
    rw [DfsOut.vertsList_cons] at hB
    simp only [List.mem_append, not_or] at hB
    refine Sim.bind_eq (Q := fun hv' s' => GuardsOuts v d rest hv' s')
      (Sim.mono (Sim.sp (Sim.mono (ih₁ bot fun e he => ⟨(hB e he).1, (hB e he).2.1, (hB e he).2.2.1⟩)
        (fun s h => h.1) fun _ _ _ h => h)) (fun _ => id) ?_)
      fun hv' => ih₂ hv' bot fun e he => ⟨(hB e he).1, (hB e he).2.1, (hB e he).2.2.2⟩
    rintro a a' s' ⟨hb, s, ⟨-, hw⟩, h0⟩
    refine ⟨hb, ?_⟩
    rw [show a' = ((walkOut v d o hasVert).run s).1 from congrArg Prod.fst h0,
      show s' = ((walkOut v d o hasVert).run s).2 from congrArg Prod.snd h0]
    exact hw

/-- The walk of a subtree never inspects entries strictly below its depth, given the tstack guards
(`GuardsTree`) along the unlifted run. -/
theorem Sim.walkTree (t : DfsTree) (d : Nat)
    (hB : ∀ e ∈ bot, e.topDepth < d ∧ e.vStart ∉ t.verts) :
    Sim bot (GuardsTree t d) (walkTree t d) (walkTree t d) fun _ _ _ => True :=
  Sim.walk_aux.1 t d bot hB

end Spqr
