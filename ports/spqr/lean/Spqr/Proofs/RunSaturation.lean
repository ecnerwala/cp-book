import Spqr.Proofs.RInv

namespace Spqr

namespace Pieces

theorem WF.sepClass_of_mem {g : Graph} {P : Pieces} (hP : P.WF g) (h2 : g.TwoConnected)
    {a b i e e' : Nat} (ha : P.Skel g a) (hb : P.Skel g b) (hnt : ¬P.TermPair a b)
    (he : P.Mem i e) (he' : P.Mem i e') : g.SepClass a b e e' := by
  have hi := hP.mem_k he
  have live : (P.x i ≠ a ∧ P.x i ≠ b) ∨ (P.y i ≠ a ∧ P.y i ≠ b) := by
    by_contra hn
    have hx : P.x i = a ∨ P.x i = b := by tauto
    have hy : P.y i = a ∨ P.y i = b := by tauto
    rcases hx with hx | hx <;> rcases hy with hy | hy
    · exact hP.ne i hi (hx.trans hy.symm)
    · exact hnt ⟨i, hi, .inl ⟨hx.symm, hy.symm⟩⟩
    · exact hnt ⟨i, hi, .inr ⟨hy.symm, hx.symm⟩⟩
    · exact hP.ne i hi (hx.trans hy.symm)
  obtain ⟨_, _, _, _, reach⟩ := hP.live_terminal h2 (skel_ok ha hb) hi live
  obtain ⟨v, hv, hr⟩ := reach e he
  obtain ⟨v', hv', hr'⟩ := reach e' he'
  exact .of_reach hv hv' (hr.trans hr'.symm)

theorem WF.sepClass_constant {g : Graph} {P : Pieces} (hP : P.WF g) (h2 : g.TwoConnected)
    {a b i e e' e₀ : Nat} (ha : P.Skel g a) (hb : P.Skel g b) (hnt : ¬P.TermPair a b)
    (he : P.Mem i e) (he' : P.Mem i e') :
    g.SepClass a b e₀ e ↔ g.SepClass a b e₀ e' := by
  have h := hP.sepClass_of_mem h2 ha hb hnt he he'
  exact ⟨fun hc => hc.trans h, fun hc => hc.trans h.symm⟩

end Pieces

def runEdges {α : Type*} (edges : α → Nat → Prop) (L : List α) (e : Nat) : Prop :=
  ∃ i ∈ L, edges i e

/-- Selecting one edge per piece transports an aligned edge interval to a run of pieces. -/
theorem runEdges_of_interval {α : Type*} {edges : α → Nat → Prop} {K : Nat → Prop}
    {order interval : List Nat} {select : Nat → Option α} {L : List α}
    (hin : interval <:+: order) (hL : order.filterMap select = L)
    (hK : ∀ e, K e ↔ e ∈ interval)
    (hconstant : ∀ m ∈ order, ∀ i, select m = some i → ∀ e, edges i e → (K m ↔ K e))
    (hown : ∀ e, K e → ∃ i, edges i e ∧ ∃ m ∈ order, select m = some i) :
    ∃ R, R <:+: L ∧ ∀ e, K e ↔ runEdges edges R e := by
  refine ⟨interval.filterMap select, hL ▸ hin.filterMap select, fun e => ⟨?_, ?_⟩⟩
  · intro he
    obtain ⟨i, hi, m, hm, hmi⟩ := hown e he
    exact ⟨i, List.mem_filterMap.2 ⟨m, (hK m).1 ((hconstant m hm i hmi e hi).2 he), hmi⟩, hi⟩
  · rintro ⟨i, hi, he⟩
    obtain ⟨m, hm, hmi⟩ := List.mem_filterMap.1 hi
    exact (hconstant m (hin.subset hm) i hmi e he).1 ((hK m).2 hm)

/-- No proper run of two or more pieces is 2-attached. -/
def Graph.RunSaturated {α : Type*} (g : Graph) (edges : α → Nat → Prop) (L : List α) : Prop :=
  ∀ R, R <:+: L → 2 ≤ R.length → R.length < L.length →
    ∀ a b, ¬g.TwoAttached (runEdges edges R) a b

namespace Graph

theorem RunSaturated.short_or_eq {α : Type*} {g : Graph} {edges : α → Nat → Prop}
    {L R : List α} (h : g.RunSaturated edges L) (hin : R <:+: L) {a b : Nat}
    (ha : g.TwoAttached (runEdges edges R) a b) : R.length ≤ 1 ∨ R = L := by
  by_cases heq : R.length = L.length
  · exact .inr (hin.eq_of_length heq)
  · left
    have hle := hin.length_le
    have hlt : R.length < L.length := by omega
    by_contra hn
    exact h R hin (by omega) hlt a b ha

theorem RunSaturated.laminar {α : Type*} {g : Graph} {edges : α → Nat → Prop}
    {L R : List α} {K : Nat → Prop} (h : g.RunSaturated edges L) (hin : R <:+: L)
    (hK : ∀ e, K e ↔ runEdges edges R e) {a b : Nat} (ha : g.TwoAttached K a b) :
    (∃ i ∈ L, ∀ e, K e → edges i e) ∨
      (∀ e, K e → ¬runEdges edges L e) ∨ (∀ e, runEdges edges L e → K e) := by
  have ha' : g.TwoAttached (runEdges edges R) a b :=
    (TwoAttached.congr fun e _ => hK e).1 ha
  rcases h.short_or_eq hin ha' with hshort | rfl
  · cases R with
    | nil =>
      exact .inr (.inl fun e he => by simpa [runEdges] using (hK e).1 he)
    | cons i R =>
      have hr : R = [] := List.eq_nil_of_length_eq_zero (by simpa using hshort)
      subst hr
      refine .inl ⟨i, hin.subset (by simp), fun e he => ?_⟩
      simpa [runEdges] using (hK e).1 he
  · exact .inr (.inr fun e he => (hK e).2 he)

end Graph

namespace Pieces

theorem WF.class_laminar_of_interval {g : Graph} {P : Pieces} (hP : P.WF g)
    (h2 : g.TwoConnected) {a b e₀ : Nat} {U : Nat → Prop}
    (ha : P.Skel g a) (hb : P.Skel g b) (hnt : ¬P.TermPair a b)
    {order interval L : List Nat} {select : Nat → Option Nat}
    (hin : interval <:+: order) (hL : order.filterMap select = L)
    (hK : ∀ e, g.SepClass a b e₀ e ↔ e ∈ interval)
    (hselect : ∀ m ∈ order, ∀ i, select m = some i → P.Mem i m)
    (hown : ∀ e, g.SepClass a b e₀ e → ∃ i ∈ L, P.Mem i e)
    (hcover : ∀ e, U e ↔ runEdges P.Mem L e)
    (hsat : g.RunSaturated P.Mem L) : P.LaminarWith U (g.SepClass a b e₀) := by
  have marked : ∀ i ∈ L, ∃ m ∈ order, select m = some i := by
    intro i hi
    rw [← hL] at hi
    exact List.mem_filterMap.1 hi
  obtain ⟨R, hR, hKR⟩ := runEdges_of_interval hin hL hK
    (fun m hm i hmi e he => hP.sepClass_constant h2 ha hb hnt (hselect m hm i hmi) he)
    (fun e he => by
      obtain ⟨i, hi, hie⟩ := hown e he
      exact ⟨i, hie, marked i hi⟩)
  rcases hsat.laminar hR hKR Graph.sepClass_twoAttached with h | h | h
  · obtain ⟨i, hi, hsub⟩ := h
    obtain ⟨m, hm, hmi⟩ := marked i hi
    exact .inl ⟨i, hP.mem_k (hselect m hm i hmi), hsub⟩
  · exact .inr (.inl fun e he hu => h e he ((hcover e).1 hu))
  · exact .inr (.inr fun e he => h e ((hcover e).1 he))

end Pieces

namespace WalkState

def rChildren (s : WalkState) (i : ItemId) : List ItemId :=
  (Items.ch s.items i).filter fun c => decide (Items.type s.items c ≠ .V)

def backEdges (dfs : DfsData) (a b e : Nat) : Prop :=
  ∃ o ∈ dfs.outs b, o.isTree = false ∧ o.dest = a ∧ o.e = e

open Classical in
noncomputable def activeEntries (s : WalkState) (L : List TEntry) : List TEntry :=
  L.filter fun t => decide (∃ e, e < s.g.ne ∧ t.edges s.g s.items e)

/-- Saturation below the schedule frontier and inside closed R items.
The whole list of an R item's children is allowed to attach at its cap. -/
structure Saturated (dfs : DfsData) (origTstack : Nat) (s : WalkState) : Prop where
  stack : ∀ R, R <:+: s.tstack.drop (s.tstack.length - origTstack) →
    2 ≤ (s.activeEntries R).length → ∀ a b,
      ¬s.g.TwoAttached (runEdges (fun t => t.edges s.g s.items) R) a b
  closed : ∀ i, i < s.items.size → Items.type s.items i = .R →
    s.g.RunSaturated (Items.EdgeBelow s.g s.items) (s.rChildren i)
  type1 : ∀ i, i < s.items.size → Items.type s.items i = .R → ∀ c ∈ s.rChildren i, ∀ a b,
    (∀ e, backEdges dfs a b e → Items.EdgeBelow s.g s.items i e) →
    (∃ e, backEdges dfs a b e ∧ ¬Items.EdgeBelow s.g s.items c e) →
    (∃ e, Items.EdgeBelow s.g s.items i e ∧
      ¬Items.EdgeBelow s.g s.items c e ∧ ¬backEdges dfs a b e) →
    ¬s.g.TwoAttached (fun e => Items.EdgeBelow s.g s.items c e ∨ backEdges dfs a b e) a b

theorem Saturated.closed_laminar {dfs : DfsData} {origTstack : Nat} {s : WalkState}
    (h : Saturated dfs origTstack s) {i : ItemId} (hi : i < s.items.size) (ht : Items.type s.items i = .R)
    {R : List ItemId} {K : Nat → Prop} (hin : R <:+: s.rChildren i)
    (hK : ∀ e, K e ↔ runEdges (Items.EdgeBelow s.g s.items) R e)
    (hcover : ∀ e, Items.EdgeBelow s.g s.items i e ↔
      runEdges (Items.EdgeBelow s.g s.items) (s.rChildren i) e)
    {a b : Nat} (ha : s.g.TwoAttached K a b) :
    (∃ c ∈ s.rChildren i, ∀ e, K e → Items.EdgeBelow s.g s.items c e) ∨
      (∀ e, K e → ¬Items.EdgeBelow s.g s.items i e) ∨
      (∀ e, Items.EdgeBelow s.g s.items i e → K e) := by
  rcases (h.closed i hi ht).laminar hin hK ha with h | h | h
  · exact .inl h
  · exact .inr (.inl fun e he ht => h e he ((hcover e).1 ht))
  · exact .inr (.inr fun e he => h e ((hcover e).1 he))

theorem entryLaminar_of_runSaturated {s : WalkState} {t : TEntry} {R : List ItemId}
    {K : Nat → Prop} (h : s.g.RunSaturated (Items.EdgeBelow s.g s.items) (s.entryPieceItems t))
    (hin : R <:+: s.entryPieceItems t) (hK : ∀ e, K e ↔ runEdges (Items.EdgeBelow s.g s.items) R e)
    (hcover : ∀ e, t.edges s.g s.items e ↔
      runEdges (Items.EdgeBelow s.g s.items) (s.entryPieceItems t) e)
    {a b : Nat} (ha : s.g.TwoAttached K a b) : s.EntryLaminar t K := by
  rcases h.laminar hin hK ha with h | h | h
  · exact .inl h
  · exact .inr (.inl fun e he ht => h e he ((hcover e).1 ht))
  · exact .inr (.inr fun e he => h e ((hcover e).1 he))

end WalkState

end Spqr
