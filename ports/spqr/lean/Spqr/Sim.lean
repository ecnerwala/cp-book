import Spqr.Walk
namespace Spqr
open WalkM

/-- `s` with `bot` underneath its tstack. -/
def lift (bot : List TEntry) (s : WalkState) : WalkState := { s with tstack := s.tstack ++ bot }

/-- `s` is `s'` with `k` extra entries at the bottom of its tstack, all other fields equal. -/
def Lifts (k : Nat) (s s' : WalkState) : Prop := ∃ bot, bot.length = k ∧ s = lift bot s'

/-- Running `m` with `bot` under the tstack simulates running `m'` without it: the final states
agree up to `bot`, and the results are related by `R` (they differ for `tstackSize`). -/
def Sim (bot : List TEntry) (P : WalkState → Prop) (m m' : WalkM α) (R : α → α → WalkState → Prop) :
    Prop :=
  ∀ s, P s → ∃ a, m.run (lift bot s) = (a, lift bot (m'.run s).2) ∧ R a (m'.run s).1 (m'.run s).2

section runLemmas
variable (s : WalkState)
theorem run_stackDir (d : Nat) : (stackDir d).run s = (s.stackDir[d]!, s) := rfl
theorem run_setStackDir (d : Nat) (b : Bool) :
    (setStackDir d b).run s = ((), { s with stackDir := s.stackDir.set! d b }) := rfl
theorem run_modifyItem (i : ItemId) (f : Item → Item) :
    (modifyItem i f).run s = ((), { s with items := s.items.modify i f }) := rfl
theorem run_getItem (i : ItemId) : (getItem i).run s = (s.items[i]!, s) := rfl
theorem run_allocItem (t : NodeType) :
    (allocItem t).run s = (s.items.size, { s with items := s.items.push ⟨t, (none, none), []⟩ }) := rfl
theorem run_makeVs (v d : Nat) :
    (makeVs v d).run s = (setSides s.stackDir[d]! (some s.stackVerts[d]!) (some v), s) := rfl
theorem run_cur : cur.run s = (s.tstack.head!, s) := rfl
theorem run_nxt : nxt.run s = (s.tstack.tail.head!, s) := rfl
theorem run_tstackSize : tstackSize.run s = (s.tstack.length, s) := rfl
theorem run_modifyCur (f : TEntry → TEntry) : (modifyCur f).run s =
    ((), { s with tstack := match s.tstack with | a :: rest => f a :: rest | [] => [] }) := rfl
theorem run_modifyNxt (f : TEntry → TEntry) : (modifyNxt f).run s =
    ((), { s with tstack := match s.tstack with | a :: b :: rest => a :: f b :: rest | l => l }) := rfl
theorem run_popTstack : popTstack.run s = (s.tstack.head!, { s with tstack := s.tstack.tail }) := rfl
theorem run_pushTstack (v d : Nat) (i : ItemId) : (pushTstack v d i).run s =
    ((), { s with tstack := ⟨v, d, s.nxtEdgeIdx, setSides s.stackDir[d]! [i] []⟩ :: s.tstack }) := rfl
theorem run_pushEdgeTstack (v d e : Nat) : (pushEdgeTstack v d e).run s =
    ((), { s with tstack := ⟨v, d, s.nxtEdgeIdx, setSides s.stackDir[d]! [edgeItem s.g e] []⟩ :: s.tstack }) :=
  rfl
end runLemmas

namespace Sim
variable {bot : List TEntry}

theorem pure (a a' : α) (h : ∀ s, P s → R a a' s) : Sim bot P (pure a) (pure a') R :=
  fun s hs => ⟨a, rfl, h s hs⟩

theorem bind {m m' : WalkM α} {f f' : α → WalkM β} (h₁ : Sim bot P m m' R)
    (h₂ : ∀ a a', Sim bot (R a a') (f a) (f' a') R₂) : Sim bot P (m >>= f) (m' >>= f') R₂ := by
  intro s hs
  obtain ⟨a, h, hr⟩ := h₁ s hs
  obtain ⟨b, hb, hb'⟩ := h₂ a _ _ hr
  exact ⟨b, by rw [StateT.run_bind, StateT.run_bind, h]; exact hb, hb'⟩

theorem mono {m m' : WalkM α} (h : Sim bot P m m' R) (hP : ∀ s, P' s → P s)
    (hR : ∀ a a' s, R a a' s → R' a a' s) : Sim bot P' m m' R' := by
  intro s hs
  obtain ⟨a, h, hr⟩ := h s (hP s hs)
  exact ⟨a, h, hR _ _ _ hr⟩

theorem ite {m₁ m₂ m₁' m₂' : WalkM α} (c c' : Prop) [Decidable c] [Decidable c']
    (hc : ∀ s, P s → (c ↔ c')) (h₁ : c → c' → Sim bot P m₁ m₁' R)
    (h₂ : ¬c → ¬c' → Sim bot P m₂ m₂' R) :
    Sim bot P (if c then m₁ else m₂) (if c' then m₁' else m₂') R := by
  intro s hs
  by_cases h : c
  · have h' := (hc s hs).mp h
    rw [ite_eq_left h, ite_eq_left h']; exact h₁ h h' s hs
  · have h' := fun h' => h ((hc s hs).mpr h')
    rw [ite_eq_right h, ite_eq_right h']; exact h₂ h h' s hs

theorem same (m : WalkM α) (h : ∀ s, P s → ∃ a, m.run (lift bot s) = (a, lift bot (m.run s).2) ∧ R a (m.run s).1 (m.run s).2) : Sim bot P m m R := h

/-- Strengthen the postcondition with the short run's exact result and state. -/
theorem sp {m m' : WalkM α} (h : Sim bot P m m' R) :
    Sim bot P m m' fun a a' s' => R a a' s' ∧ ∃ s, P s ∧ (a', s') = m'.run s := by
  intro s hs
  obtain ⟨a, h, hr⟩ := h s hs
  exact ⟨a, h, hr, s, hs, rfl⟩

theorem bind_eq {m m' : WalkM α} {f f' : α → WalkM β} {Q : α → WalkState → Prop}
    (h₁ : Sim bot P m m' fun a a' s => a = a' ∧ Q a' s)
    (h₂ : ∀ a, Sim bot (Q a) (f a) (f' a) R₂) : Sim bot P (m >>= f) (m' >>= f') R₂ :=
  bind h₁ fun a a' s ⟨hs, hq⟩ => by subst hs; exact h₂ a s hq

theorem bind_get {f f' : WalkState → WalkM β}
    (h : ∀ a, Sim bot (fun s => P s ∧ s = a) (f (lift bot a)) (f' a) R) :
    Sim bot P (get >>= f) (get >>= f') R :=
  bind (m := get) (m' := get) (R := fun a a' s => a = lift bot a' ∧ a' = s ∧ P s)
    (fun s hs => ⟨_, rfl, rfl, rfl, hs⟩) fun a a' s ⟨ha, ha', hs⟩ => ha ▸ h a' s ⟨hs, ha'.symm⟩

theorem bind_tstackSize {f f' : Nat → WalkM β}
    (h : ∀ n, Sim bot (fun s => P s ∧ n = s.tstack.length) (f (n + bot.length)) (f' n) R) :
    Sim bot P (tstackSize >>= f) (tstackSize >>= f') R :=
  bind (m := tstackSize) (m' := tstackSize)
    (R := fun n n' s => n = n' + bot.length ∧ n' = s.tstack.length ∧ P s)
    (fun s hs => ⟨_, rfl, by simp [run_tstackSize, lift], rfl, hs⟩)
    fun a a' s ⟨ha, ha', hs⟩ => ha ▸ h a' s ⟨hs, ha'⟩

theorem map {m m' : WalkM α} (g : α → β) (h : Sim bot P m m' R) :
    Sim bot P (g <$> m) (g <$> m') fun b b' s => ∃ a a', b = g a ∧ b' = g a' ∧ R a a' s := by
  intro s hs
  obtain ⟨a, h, hr⟩ := h s hs
  refine ⟨g a, ?_, a, _, rfl, rfl, hr⟩
  show (g (m.run (lift bot s)).1, (m.run (lift bot s)).2) = (g a, lift bot (m'.run s).2)
  rw [h]

theorem iteb {m₁ m₂ m₁' m₂' : WalkM α} (b : Bool) (h₁ : b = true → Sim bot P m₁ m₁' R)
    (h₂ : b = false → Sim bot P m₂ m₂' R) :
    Sim bot P (if b = true then m₁ else m₂) (if b = true then m₁' else m₂') R := by
  cases b
  · exact h₂ rfl
  · exact h₁ rfl

theorem iteb_eq {m₁ m₂ m₁' m₂' : WalkM α} (b b' : Bool) (hb : ∀ s, P s → b = b')
    (h₁ : b' = true → Sim bot P m₁ m₁' R) (h₂ : b' = false → Sim bot P m₂ m₂' R) :
    Sim bot P (if b = true then m₁ else m₂) (if b' = true then m₁' else m₂') R := by
  intro s hs
  rw [hb s hs]
  cases b'
  · exact h₂ rfl s hs
  · exact h₁ rfl s hs

theorem seq {m m' : WalkM Unit} {f f' : WalkM β} (h₁ : Sim bot P m m' fun _ _ s => Q s)
    (h₂ : Sim bot Q f f' R₂) : Sim bot P (m >>= fun _ => f) (m' >>= fun _ => f') R₂ :=
  bind h₁ fun _ _ => h₂

/-- The precondition only matters on the short side; the results need not be related. -/
theorem drop {m m' : WalkM α} (h : Sim bot P m m' R) : Sim bot P m m' fun _ _ _ => True :=
  mono h (fun _ => id) fun _ _ _ _ => trivial

end Sim

@[simp] theorem lift_tstack (bot : List TEntry) (s : WalkState) : (lift bot s).tstack = s.tstack ++ bot := rfl
@[simp] theorem lift_g (bot : List TEntry) (s : WalkState) : (lift bot s).g = s.g := rfl
@[simp] theorem lift_ternarize (bot : List TEntry) (s : WalkState) : (lift bot s).ternarize = s.ternarize := rfl
@[simp] theorem lift_firstOccurrence (bot : List TEntry) (s : WalkState) :
    (lift bot s).firstOccurrence = s.firstOccurrence := rfl
@[simp] theorem lift_nxtEdgeIdx (bot : List TEntry) (s : WalkState) : (lift bot s).nxtEdgeIdx = s.nxtEdgeIdx := rfl
@[simp] theorem lift_stackDir (bot : List TEntry) (s : WalkState) : (lift bot s).stackDir = s.stackDir := rfl
@[simp] theorem lift_items (bot : List TEntry) (s : WalkState) : (lift bot s).items = s.items := rfl


/-- A state update that does not look at the tstack. -/
theorem Sim.modify_of (f : WalkState → WalkState) (h : ∀ s, f (lift bot s) = lift bot (f s)) :
    Sim bot P (modify f) (modify f) fun _ _ _ => True := by
  intro s _; exact ⟨(), by show ((), f (lift bot s)) = ((), lift bot (f s)); rw [h], trivial⟩

theorem Sim.stackDir (d : Nat) : Sim bot P (stackDir d) (stackDir d) fun a a' _ => a = a' :=
  fun s _ => ⟨_, rfl, rfl⟩

theorem Sim.setStackDir (d : Nat) (b : Bool) :
    Sim bot P (setStackDir d b) (setStackDir d b) fun _ _ _ => True :=
  Sim.modify_of _ fun _ => rfl

theorem Sim.modifyItem (i : ItemId) (f : Item → Item) :
    Sim bot P (modifyItem i f) (modifyItem i f) fun _ _ _ => True :=
  Sim.modify_of _ fun _ => rfl

theorem Sim.getItem (i : ItemId) : Sim bot P (getItem i) (getItem i) fun a a' _ => a = a' :=
  fun s _ => ⟨_, rfl, rfl⟩

theorem Sim.allocItem (t : NodeType) : Sim bot P (allocItem t) (allocItem t) fun a a' _ => a = a' :=
  fun s _ => ⟨_, rfl, rfl⟩

theorem Sim.makeVs (v d : Nat) : Sim bot P (makeVs v d) (makeVs v d) fun a a' _ => a = a' :=
  fun s _ => ⟨_, rfl, rfl⟩

theorem Sim.pushTstack (v d : Nat) (i : ItemId) :
    Sim bot P (pushTstack v d i) (pushTstack v d i) fun _ _ _ => True :=
  Sim.modify_of _ fun _ => rfl

theorem Sim.pushVertTstack (v d : Nat) :
    Sim bot P (pushVertTstack v d) (pushVertTstack v d) fun _ _ _ => True :=
  Sim.pushTstack v d _

theorem Sim.pushEdgeTstack (v d e : Nat) :
    Sim bot P (pushEdgeTstack v d e) (pushEdgeTstack v d e) fun _ _ _ => True :=
  fun s _ => ⟨_, rfl, trivial⟩

theorem Sim.tstackSize : Sim bot P tstackSize tstackSize fun n n' _ => n = n' + bot.length :=
  fun s _ => ⟨_, rfl, by simp [run_tstackSize]⟩

/-- Destructure `s` with a known tstack so both runs compute by `rfl`. -/
syntax "sim_rfl" ident ident : tactic
macro_rules
  | `(tactic| sim_rfl $s $h) =>
    `(tactic| (rcases $s:ident with ⟨g, tern, items, sv, sd, nei, fo, ts, tb, tsl⟩; simp only at $h:ident; subst $h:ident; rfl))

theorem Sim.cur (hP : ∀ s, P s → s.tstack ≠ []) : Sim bot P cur cur fun a a' _ => a = a' := by
  intro s hs
  obtain ⟨x, l, h⟩ := List.exists_cons_of_ne_nil (hP s hs)
  exact ⟨_, rfl, by sim_rfl s h⟩

theorem Sim.nxt (hP : ∀ s, P s → 2 ≤ s.tstack.length) : Sim bot P nxt nxt fun a a' _ => a = a' := by
  intro s hs
  obtain ⟨x, y, l, h⟩ : ∃ x y l, s.tstack = x :: y :: l := by
    match hl : s.tstack, hP s hs with
    | x :: y :: l, _ => exact ⟨x, y, l, rfl⟩
  exact ⟨_, rfl, by sim_rfl s h⟩

theorem Sim.popTstack (hP : ∀ s, P s → s.tstack ≠ []) :
    Sim bot P popTstack popTstack fun a a' _ => a = a' := by
  intro s hs
  obtain ⟨x, l, h⟩ := List.exists_cons_of_ne_nil (hP s hs)
  exact ⟨_, by sim_rfl s h, by sim_rfl s h⟩

theorem Sim.modifyCur (f : TEntry → TEntry) (hP : ∀ s, P s → s.tstack ≠ []) :
    Sim bot P (modifyCur f) (modifyCur f) fun _ _ _ => True := by
  intro s hs
  obtain ⟨x, l, h⟩ := List.exists_cons_of_ne_nil (hP s hs)
  exact ⟨(), by sim_rfl s h, trivial⟩

theorem Sim.modifyNxt (f : TEntry → TEntry) (hP : ∀ s, P s → 2 ≤ s.tstack.length) :
    Sim bot P (modifyNxt f) (modifyNxt f) fun _ _ _ => True := by
  intro s hs
  obtain ⟨x, y, l, h⟩ : ∃ x y l, s.tstack = x :: y :: l := by
    match hl : s.tstack, hP s hs with
    | x :: y :: l, _ => exact ⟨x, y, l, rfl⟩
  exact ⟨(), by sim_rfl s h, trivial⟩

theorem Sim.mergeTstackTops (hP : ∀ s, P s → 2 ≤ s.tstack.length) :
    Sim bot P mergeTstackTops mergeTstackTops fun _ _ _ => True := by
  intro s hs
  obtain ⟨x, y, l, h⟩ : ∃ x y l, s.tstack = x :: y :: l := by
    match hl : s.tstack, hP s hs with
    | x :: y :: l, _ => exact ⟨x, y, l, rfl⟩
  exact ⟨(), by sim_rfl s h, trivial⟩

theorem Sim.finishTstackTop (i : ItemId) (hP : ∀ s, P s → s.tstack ≠ []) :
    Sim bot P (finishTstackTop i) (finishTstackTop i) fun _ _ _ => True := by
  intro s hs
  obtain ⟨x, l, h⟩ := List.exists_cons_of_ne_nil (hP s hs)
  exact ⟨(), by sim_rfl s h, trivial⟩

theorem Sim.get : Sim bot P get get fun a a' s' => a = lift bot a' ∧ a' = s' ∧ P s' :=
  fun s hs => ⟨_, rfl, rfl, rfl, hs⟩

/-- `loop` on both sides, with `I n` meaning at most `n` more iterations remain (so the short
side's fuel suffices) and `J n` the states where the body runs with `n` iterations of fuel left. -/
theorem Sim.loop (I J : Nat → WalkState → Prop) (Q : WalkState → Prop)
    (cond cond' : WalkM Bool) (body body' : WalkM Unit)
    (hpure : ∀ s, (cond'.run s).2 = s)
    (hcond : ∀ n, Sim bot (I n) cond cond' fun b b' s => b = b' ∧ (b' = true → J n s) ∧ (b' = false → Q s))
    (hJ0 : ∀ s, ¬ J 0 s)
    (hbody : ∀ n, Sim bot (J (n + 1)) body body' fun _ _ s => I n s) :
    ∀ n m, n ≤ m → Sim bot (I n) (loop m cond body) (loop n cond' body') fun _ _ s => Q s := by
  intro n
  induction n with
  | zero =>
    intro m _ s hs
    obtain ⟨b, hb, hbb, hJ, hQ⟩ := hcond 0 s hs
    rw [hpure] at hb hJ hQ
    have hb' : (cond'.run s).1 = false := by
      cases h : (cond'.run s).1
      · rfl
      · exact absurd (hJ h) (hJ0 s)
    refine ⟨(), ?_, hQ hb'⟩
    cases m with
    | zero => rfl
    | succ m =>
      show (cond >>= fun c => if c then body >>= fun _ => WalkM.loop m cond body else Pure.pure ()).run
        (lift bot s) = _
      rw [StateT.run_bind, hb, hbb, hb']; rfl
  | succ n ih =>
    intro m hm s hs
    obtain ⟨m, rfl⟩ : ∃ m', m = m' + 1 := ⟨m - 1, by omega⟩
    have hp := hpure s
    obtain ⟨b, hb, hbb, hJ, hQ⟩ := hcond (n + 1) s hs
    have e1 : (WalkM.loop (m + 1) cond body).run (lift bot s) =
        (cond >>= fun c => if c then body >>= fun _ => WalkM.loop m cond body else Pure.pure ()).run
          (lift bot s) := rfl
    have e2 : (WalkM.loop (n + 1) cond' body').run s =
        (cond' >>= fun c => if c then body' >>= fun _ => WalkM.loop n cond' body' else Pure.pure ()).run
          s := rfl
    rw [e1, e2, StateT.run_bind, StateT.run_bind, hb, hbb]
    generalize hc : cond'.run s = r at hp hJ hQ ⊢
    obtain ⟨c, s'⟩ := r
    simp only at hp hJ hQ
    subst hp
    cases c
    · exact ⟨(), rfl, hQ rfl⟩
    · obtain ⟨⟨⟩, hb2, hI⟩ := hbody n s' (hJ rfl)
      show ∃ a, (body >>= fun _ => WalkM.loop m cond body).run (lift bot s') =
          (a, lift bot ((body' >>= fun _ => WalkM.loop n cond' body').run s').2) ∧
        Q ((body' >>= fun _ => WalkM.loop n cond' body').run s').2
      rw [StateT.run_bind, StateT.run_bind, hb2]
      generalize body'.run s' = r at hI ⊢
      obtain ⟨_, s₂⟩ := r
      obtain ⟨⟨⟩, hl, hQ'⟩ := ih m (by omega) s₂ hI
      exact ⟨(), by show (WalkM.loop m cond body).run (lift bot s₂) = _; rw [hl]; rfl, hQ'⟩

end Spqr
