import Mathlib.Logic.Relation
import Spqr.ItemSpec

/-!
# Connected and 2-attached edge sets

The pure graph-theory core of the walk invariant (`WalkSpec.lean`): edge sets `E : Nat → Prop` of
`g` that are connected, and whose attachment to the rest of `g` is confined to two terminals.
-/

namespace Spqr

namespace Graph

variable (g : Graph)

/-- Edge `e` is incident to `v`. -/
def Inc (e v : Nat) : Prop := (g.edges[e]!).1 = v ∨ (g.edges[e]!).2 = v

/-- `u` and `w` are joined by an edge of `E`. -/
def AdjIn (E : Nat → Prop) (u w : Nat) : Prop :=
  ∃ e, e < g.ne ∧ E e ∧ Items.PairEq (u, w) g.edges[e]!

/-- The edge set `E` is connected (as a subgraph; isolated vertices are not part of it). -/
def ConnEdges (E : Nat → Prop) : Prop :=
  ∀ e e', e < g.ne → e' < g.ne → E e → E e' →
    Relation.ReflTransGen (g.AdjIn E) (g.edges[e]!).1 (g.edges[e']!).1

/-- Every vertex incident to both an edge in `E` and an edge outside `E` is `a` or `b`. -/
def TwoAttached (E : Nat → Prop) (a b : Nat) : Prop :=
  ∀ v e e', e < g.ne → e' < g.ne → E e → ¬ E e' → g.Inc e v → g.Inc e' v → v = a ∨ v = b

/-- Every vertex incident to both an edge in `E` and an edge outside `E` satisfies `S`. -/
def AttachedIn (E : Nat → Prop) (S : Nat → Prop) : Prop :=
  ∀ v e e', e < g.ne → e' < g.ne → E e → ¬ E e' → g.Inc e v → g.Inc e' v → S v

/-- Every edge at `v` is in `E`. -/
def Interior (E : Nat → Prop) (v : Nat) : Prop := ∀ e, e < g.ne → g.Inc e v → E e

/-- Some edge of `E` is incident to `v`. -/
def Touches (E : Nat → Prop) (v : Nat) : Prop := ∃ e, e < g.ne ∧ E e ∧ g.Inc e v

variable {g}
variable {E E₁ E₂ : Nat → Prop}

theorem AdjIn.symm {u w : Nat} (h : g.AdjIn E u w) : g.AdjIn E w u := by
  obtain ⟨e, he, hE, hp⟩ := h
  refine ⟨e, he, hE, ?_⟩
  rcases hp with hp | hp
  · exact .inr (by rw [Prod.ext_iff] at hp ⊢; exact ⟨hp.2, hp.1⟩)
  · exact .inl (by rw [Prod.ext_iff] at hp ⊢; exact ⟨hp.2, hp.1⟩)

theorem AdjIn.mono (h : ∀ e, e < g.ne → E₁ e → E₂ e) {u w : Nat} (ha : g.AdjIn E₁ u w) :
    g.AdjIn E₂ u w := by
  obtain ⟨e, he, hE, hp⟩ := ha
  exact ⟨e, he, h e he hE, hp⟩

theorem reach_symm {u w : Nat} (h : Relation.ReflTransGen (g.AdjIn E) u w) :
    Relation.ReflTransGen (g.AdjIn E) w u := by
  induction h with
  | refl => exact .refl
  | tail _ h ih => exact .head h.symm ih

theorem reach_mono (h : ∀ e, e < g.ne → E₁ e → E₂ e) {u w : Nat} (hr : Relation.ReflTransGen (g.AdjIn E₁) u w) :
    Relation.ReflTransGen (g.AdjIn E₂) u w := by
  induction hr with
  | refl => exact .refl
  | tail _ h' ih => exact .tail ih (h'.mono h)

/-- An edge of `E` reaches both its endpoints. -/
theorem reach_of_inc {e v : Nat} (he : e < g.ne) (hE : E e) (hv : g.Inc e v) :
    Relation.ReflTransGen (g.AdjIn E) (g.edges[e]!).1 v := by
  rcases hv with hv | hv
  · rw [hv]
  · exact .single ⟨e, he, hE, .inl (by rw [← hv])⟩

theorem ConnEdges.congr (h : ∀ e, e < g.ne → (E₁ e ↔ E₂ e)) : g.ConnEdges E₁ ↔ g.ConnEdges E₂ := by
  constructor <;> intro hc e e' he he' hE hE'
  · exact reach_mono (fun e he => (h e he).1) (hc e e' he he' ((h e he).2 hE) ((h e' he').2 hE'))
  · exact reach_mono (fun e he => (h e he).2) (hc e e' he he' ((h e he).1 hE) ((h e' he').1 hE'))

theorem TwoAttached.congr (h : ∀ e, e < g.ne → (E₁ e ↔ E₂ e)) {a b : Nat} :
    g.TwoAttached E₁ a b ↔ g.TwoAttached E₂ a b := by
  constructor <;> intro ha v e e' he he' hE hE' hv hv'
  · exact ha v e e' he he' ((h e he).2 hE) (fun h' => hE' ((h e' he').1 h')) hv hv'
  · exact ha v e e' he he' ((h e he).1 hE) (fun h' => hE' ((h e' he').2 h')) hv hv'

theorem Interior.congr (h : ∀ e, e < g.ne → (E₁ e ↔ E₂ e)) {v : Nat} :
    g.Interior E₁ v ↔ g.Interior E₂ v := by
  constructor <;> intro hi e he hv
  · exact (h e he).1 (hi e he hv)
  · exact (h e he).2 (hi e he hv)

theorem TwoAttached.comm {a b : Nat} (h : g.TwoAttached E a b) : g.TwoAttached E b a :=
  fun v e e' he he' hE hE' hv hv' => (h v e e' he he' hE hE' hv hv').symm

/-- A terminal with no outside edge can be replaced by anything. -/
theorem TwoAttached.of_interior_left {a b c : Nat} (h : g.TwoAttached E a b) (hi : g.Interior E a) :
    g.TwoAttached E c b := by
  intro v e e' he he' hE hE' hv hv'
  rcases h v e e' he he' hE hE' hv hv' with rfl | rfl
  · exact absurd (hi e' he' hv') hE'
  · exact .inr rfl

theorem ConnEdges.empty (h : ∀ e, e < g.ne → ¬ E e) : g.ConnEdges E :=
  fun e _ he _ hE _ => absurd hE (h e he)

theorem TwoAttached.empty {a b : Nat} (h : ∀ e, e < g.ne → ¬ E e) : g.TwoAttached E a b :=
  fun _ e _ he _ hE _ _ _ => absurd hE (h e he)

/-- A single edge `e₀ = (a, b)` is 2-attached at its endpoints. -/
theorem TwoAttached.single {e₀ : Nat} (h : ∀ e, E e ↔ e = e₀) :
    g.TwoAttached E (g.edges[e₀]!).1 (g.edges[e₀]!).2 := by
  intro v e _ _ _ hE _ hv _
  rw [(h e).1 hE] at hv
  rcases hv with rfl | rfl <;> simp

theorem ConnEdges.single {e₀ : Nat} (h : ∀ e, E e ↔ e = e₀) : g.ConnEdges E := by
  intro e e' _ _ hE hE'
  rw [(h e).1 hE, (h e').1 hE']

/-- Two connected edge sets that share a vertex (unless one is empty) have a connected union. -/
theorem ConnEdges.union (h₁ : g.ConnEdges E₁) (h₂ : g.ConnEdges E₂)
    (hx : (∃ e, e < g.ne ∧ E₁ e) → (∃ e, e < g.ne ∧ E₂ e) → ∃ x, g.Touches E₁ x ∧ g.Touches E₂ x) :
    g.ConnEdges fun e => E₁ e ∨ E₂ e := by
  have cross : ∀ e e', e < g.ne → e' < g.ne → E₁ e → E₂ e' →
      Relation.ReflTransGen (g.AdjIn fun e => E₁ e ∨ E₂ e) (g.edges[e]!).1 (g.edges[e']!).1 := by
    intro e e' he he' hE hE'
    obtain ⟨x, ⟨e₁, he₁, hE₁, hx₁⟩, ⟨e₂, he₂, hE₂, hx₂⟩⟩ := hx ⟨e, he, hE⟩ ⟨e', he', hE'⟩
    have p₁ := reach_mono (E₂ := fun e => E₁ e ∨ E₂ e) (fun _ _ => Or.inl) (h₁ e e₁ he he₁ hE hE₁)
    have p₂ := reach_mono (E₂ := fun e => E₁ e ∨ E₂ e) (fun _ _ => Or.inl) (reach_of_inc he₁ hE₁ hx₁)
    have p₃ := reach_mono (E₂ := fun e => E₁ e ∨ E₂ e) (fun _ _ => Or.inr) (reach_of_inc he₂ hE₂ hx₂)
    have p₄ := reach_mono (E₂ := fun e => E₁ e ∨ E₂ e) (fun _ _ => Or.inr) (h₂ e₂ e' he₂ he' hE₂ hE')
    exact p₁.trans (p₂.trans ((reach_symm p₃).trans p₄))
  intro e e' he he' hE hE'
  rcases hE with hE | hE <;> rcases hE' with hE' | hE'
  · exact reach_mono (fun _ _ => Or.inl) (h₁ e e' he he' hE hE')
  · exact cross e e' he he' hE hE'
  · exact reach_symm (cross e' e he' he hE' hE)
  · exact reach_mono (fun _ _ => Or.inr) (h₂ e e' he he' hE hE')

/-- The union of two 2-attached edge sets is 2-attached at `a, b` once every old terminal is
either a new terminal or interior to the union. -/
theorem TwoAttached.union {a₁ b₁ a₂ b₂ a b : Nat}
    (h₁ : g.TwoAttached E₁ a₁ b₁) (h₂ : g.TwoAttached E₂ a₂ b₂)
    (hv : ∀ v ∈ [a₁, b₁, a₂, b₂], v = a ∨ v = b ∨ g.Interior (fun e => E₁ e ∨ E₂ e) v) :
    g.TwoAttached (fun e => E₁ e ∨ E₂ e) a b := by
  intro v e e' he he' hE hE' hv' hv''
  have key : ∀ w ∈ [a₁, b₁, a₂, b₂], v = w → v = a ∨ v = b := by
    rintro w hw rfl
    rcases hv v hw with h | h | h
    · exact .inl h
    · exact .inr h
    · exact absurd (h e' he' hv'') hE'
  rcases hE with hE | hE
  · rcases h₁ v e e' he he' hE (fun h => hE' (.inl h)) hv' hv'' with h | h
    · exact key _ (by simp) h
    · exact key _ (by simp) h
  · rcases h₂ v e e' he he' hE (fun h => hE' (.inr h)) hv' hv'' with h | h
    · exact key _ (by simp) h
    · exact key _ (by simp) h

theorem twoAttached_iff {a b : Nat} : g.TwoAttached E a b ↔ g.AttachedIn E fun v => v = a ∨ v = b := Iff.rfl

theorem AttachedIn.congr (h : ∀ e, e < g.ne → (E₁ e ↔ E₂ e)) {S : Nat → Prop} :
    g.AttachedIn E₁ S ↔ g.AttachedIn E₂ S := by
  constructor <;> intro ha v e e' he he' hE hE' hv hv'
  · exact ha v e e' he he' ((h e he).2 hE) (fun h' => hE' ((h e' he').1 h')) hv hv'
  · exact ha v e e' he he' ((h e he).1 hE) (fun h' => hE' ((h e' he').2 h')) hv hv'

theorem AttachedIn.mono {S₁ S₂ : Nat → Prop} (h : ∀ v, S₁ v → S₂ v) (ha : g.AttachedIn E S₁) :
    g.AttachedIn E S₂ :=
  fun v e e' he he' hE hE' hv hv' => h v (ha v e e' he he' hE hE' hv hv')

theorem AttachedIn.empty {S : Nat → Prop} (h : ∀ e, e < g.ne → ¬ E e) : g.AttachedIn E S :=
  fun _ e _ he _ hE _ _ _ => absurd hE (h e he)

/-- An attachment vertex of the union is an attachment vertex of a part that is not interior to
the union. -/
theorem AttachedIn.union {S₁ S₂ : Nat → Prop} (h₁ : g.AttachedIn E₁ S₁) (h₂ : g.AttachedIn E₂ S₂) :
    g.AttachedIn (fun e => E₁ e ∨ E₂ e)
      fun v => (S₁ v ∨ S₂ v) ∧ ¬ g.Interior (fun e => E₁ e ∨ E₂ e) v := by
  intro v e e' he he' hE hE' hv hv'
  refine ⟨?_, fun hi => hE' (hi e' he' hv')⟩
  rcases hE with hE | hE
  · exact .inl (h₁ v e e' he he' hE (fun h => hE' (.inl h)) hv hv')
  · exact .inr (h₂ v e e' he he' hE (fun h => hE' (.inr h)) hv hv')

/-- Attachment vertices touch `E` and are not interior to it. -/
theorem AttachedIn.strengthen {S : Nat → Prop} (h : g.AttachedIn E S) :
    g.AttachedIn E fun v => S v ∧ g.Touches E v ∧ ¬ g.Interior E v :=
  fun v e e' he he' hE hE' hv hv' =>
    ⟨h v e e' he he' hE hE' hv hv', ⟨e, he, hE, hv⟩, fun hi => hE' (hi e' he' hv')⟩

end Graph

end Spqr
