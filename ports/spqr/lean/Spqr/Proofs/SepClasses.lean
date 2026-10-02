import Mathlib.Tactic.Push
import Spqr.Proofs.Contract
import Spqr.Proofs.SepPairExhaust

/-!
# Walk pieces are separation classes (PROOF.md §4.2b ↔ §3)

The walk invariant (`WalkSpec.lean`) describes each tstack entry / closed item by an edge set `E`
that is connected and 2-attached at two terminals `{u, v}`. Over a block, such an `E` (nonempty,
not everything) is a union of separation classes of `{u, v}` (`TwoAttached.sepClass_mem`), the
terminals are distinct and both touched by `E`, and `{u, v}` is a separation pair — unless `E` or
its complement is a single `u–v` edge (`TwoAttached.sepPair_or_single`). Through `sepPair_iff'`
the pair is then a type-1 or type-2 pair of the DFS tree (`DfsData.twoAttached_type1_or_type2`).
Connectedness of `E` is not needed for any of this.
-/

namespace Spqr

namespace Graph

variable {g : Graph} {E : Nat → Prop}

theorem TwoAttached.compl {a b : Nat} (h : g.TwoAttached E a b) :
    g.TwoAttached (fun e => ¬E e) a b :=
  fun w e e' he he' hE hE' hv hv' => h w e' e he' he (not_not.1 hE') hE hv' hv

/-- A walk from a vertex of `E` either stays at vertices of `E` or passes a terminal of `E`. -/
theorem TwoAttached.reach_cross {ok : Nat → Prop} {u v x y : Nat} (h : g.TwoAttached E u v)
    (hr : g.Reach ok x y) (hx : g.Touches E x) :
    g.Touches E y ∨ ∃ w, ok w ∧ (w = u ∨ w = v) ∧ g.Touches E w := by
  induction hr with
  | refl _ => exact .inl hx
  | tail hr hadj _ ih =>
    rcases ih with ih | ih
    · obtain ⟨f, hj⟩ := hadj
      by_cases hf : E f
      · exact .inl (Touches.of_isEnd hf hj.symm.isEnd)
      · have ⟨e₁, he₁, hE₁, hy⟩ := ih
        have hf' := isEnd_iff.1 hj.isEnd
        exact .inr ⟨_, hr.ok_right, h _ e₁ f he₁ hf'.1 hE₁ hf hy hf'.2, ih⟩
    · exact .inr ih

/-- In a block, a nonempty proper 2-attached `E` touches a terminal other than any given `z`. -/
theorem TwoAttached.exists_touch_ne (h2 : g.TwoConnected) {u v : Nat} (h : g.TwoAttached E u v)
    {e₁ e₂ : Nat} (he₁ : e₁ < g.ne) (hE₁ : E e₁) (he₂ : e₂ < g.ne) (hE₂ : ¬E e₂) (z : Nat) :
    ∃ w, w ≠ z ∧ (w = u ∨ w = v) ∧ g.Touches E w := by
  rcases h2 z e₁ e₂ he₁ he₂ with rfl | ⟨x, y, hx, hy, hr⟩
  · exact absurd hE₁ hE₂
  · rcases h.reach_cross hr (Touches.of_isEnd hE₁ hx) with ht | ⟨w, hw, huv, ht⟩
    · have ⟨e₃, he₃, hE₃, hy₃⟩ := ht
      have hy' := isEnd_iff.1 hy
      exact ⟨y, hr.ok_right, h _ e₃ e₂ he₃ hy'.1 hE₃ hE₂ hy₃ hy'.2, ht⟩
    · exact ⟨w, hw, huv, ht⟩

theorem TwoAttached.ne (h2 : g.TwoConnected) {u v : Nat} (h : g.TwoAttached E u v)
    {e₁ e₂ : Nat} (he₁ : e₁ < g.ne) (hE₁ : E e₁) (he₂ : e₂ < g.ne) (hE₂ : ¬E e₂) : u ≠ v := by
  obtain ⟨w, hw, huv, -⟩ := h.exists_touch_ne h2 he₁ hE₁ he₂ hE₂ u
  rcases huv with rfl | rfl
  · exact absurd rfl hw
  · exact hw.symm

theorem TwoAttached.touches (h2 : g.TwoConnected) {u v : Nat} (h : g.TwoAttached E u v)
    {e₁ e₂ : Nat} (he₁ : e₁ < g.ne) (hE₁ : E e₁) (he₂ : e₂ < g.ne) (hE₂ : ¬E e₂) :
    g.Touches E u ∧ g.Touches E v := by
  constructor
  · obtain ⟨w, hw, huv, ht⟩ := h.exists_touch_ne h2 he₁ hE₁ he₂ hE₂ v
    rcases huv with rfl | rfl
    · exact ht
    · exact absurd rfl hw
  · obtain ⟨w, hw, huv, ht⟩ := h.exists_touch_ne h2 he₁ hE₁ he₂ hE₂ u
    rcases huv with rfl | rfl
    · exact absurd rfl hw
    · exact ht

/-- In a block, a 2-attached `E = {e₁}` with something outside it is a `u–v` edge. -/
theorem TwoAttached.joins_of_single (h2 : g.TwoConnected) {u v : Nat} (h : g.TwoAttached E u v)
    {e₁ e₂ : Nat} (he₁ : e₁ < g.ne) (hE₁ : E e₁) (hone : ∀ e, e < g.ne → E e → e = e₁)
    (he₂ : e₂ < g.ne) (hE₂ : ¬E e₂) : g.Joins e₁ u v := by
  have hends : ∀ x y, g.Joins e₁ x y → x = u ∨ x = v := by
    intro x y hj
    by_contra hx
    push Not at hx
    have hall : ∀ f, g.IsEnd f x → f = e₁ := by
      intro f hf
      have hf' := isEnd_iff.1 hf
      by_contra hfe
      have hEf : ¬E f := fun hEf => hfe (hone f hf'.1 hEf)
      rcases h x e₁ f he₁ hf'.1 hE₁ hEf (isEnd_iff.1 hj.isEnd).2 hf'.2 with h' | h'
      · exact hx.1 h'
      · exact hx.2 h'
    have hstay : ∀ x' z, g.Reach (· ≠ y) x' z → x' = x → z = x := by
      intro x' z hz
      induction hz with
      | refl _ => exact id
      | tail _ hadj hok ih =>
        intro hx'
        obtain ⟨f, hf⟩ := hadj
        rw [ih hx'] at hf
        have := hall f hf.isEnd
        subst this
        rcases hf.eq_or hj with ⟨-, h'⟩ | ⟨-, h'⟩
        · exact absurd h' hok
        · exact h'
    rcases h2 y e₁ e₂ he₁ he₂ with rfl | ⟨x', y', hx', hy', hr⟩
    · exact hE₂ hE₁
    · have hx'x : x' = x := by
        rcases hx'.eq_or hj with h' | h'
        · exact h'
        · exact absurd h' hr.ok_left
      rw [hstay x' y' hr hx'x] at hy'
      exact hE₂ (hall e₂ hy' ▸ hE₁)
  rcases hpq : g.edges[e₁]! with ⟨p, q⟩
  have hj : g.Joins e₁ p q := joins_iff.2 ⟨he₁, .inl hpq⟩
  have hp := hends p q hj
  have hq := hends q p hj.symm
  have hpq' : p ≠ q := by
    rintro rfl
    rcases h2 p e₁ e₂ he₁ he₂ with rfl | ⟨x', _, hx', _, hr⟩
    · exact hE₂ hE₁
    · rcases hx'.eq_or hj with h' | h' <;> exact hr.ok_left h'
  rcases hp with rfl | rfl <;> rcases hq with rfl | rfl
  · exact absurd rfl hpq'
  · exact hj
  · exact hj.symm
  · exact absurd rfl hpq'

/-- A nonempty proper 2-attached edge set of a block is cut off by a separation pair, unless it or
its complement is a single `u–v` edge. -/
theorem TwoAttached.sepPair_or_single (h2 : g.TwoConnected) {u v : Nat} (h : g.TwoAttached E u v)
    (hne : ∃ e, e < g.ne ∧ E e) (hpr : ∃ e, e < g.ne ∧ ¬E e) :
    g.SeparationPair u v ∨ (∃ e₀, g.Joins e₀ u v ∧ ∀ e, e < g.ne → (E e ↔ e = e₀)) ∨
      (∃ e₀, g.Joins e₀ u v ∧ ∀ e, e < g.ne → (¬E e ↔ e = e₀)) := by
  obtain ⟨e₁, he₁, hE₁⟩ := hne
  obtain ⟨e₂, he₂, hE₂⟩ := hpr
  have huv := h.ne h2 he₁ hE₁ he₂ hE₂
  have cl : ∀ {e e'}, E e → g.SepClass u v e e' → E e' := fun hE hc => (h.sepClass_mem hc).1 hE
  have cl' : ∀ {e e'}, ¬E e → g.SepClass u v e e' → ¬E e' :=
    fun hE hc hE' => hE ((h.sepClass_mem hc).2 hE')
  have n12 : ¬g.SepClass u v e₁ e₂ := fun hc => hE₂ (cl hE₁ hc)
  by_cases hA : ∃ e, e < g.ne ∧ E e ∧ ¬g.SepClass u v e₁ e
  · obtain ⟨e, he, hE, hc⟩ := hA
    exact .inl ⟨huv, .inl ⟨e₁, e, e₂, he₁, he, he₂, hc, n12, fun hc => hE₂ (cl hE hc)⟩⟩
  by_cases hB : ∃ e, e < g.ne ∧ ¬E e ∧ ¬g.SepClass u v e₂ e
  · obtain ⟨e, he, hE, hc⟩ := hB
    exact .inl ⟨huv, .inl ⟨e₂, e, e₁, he₂, he, he₁, hc, fun hc => n12 hc.symm,
      fun hc => cl' hE hc hE₁⟩⟩
  push Not at hA hB
  by_cases hC : ∃ e, e < g.ne ∧ E e ∧ e ≠ e₁
  · obtain ⟨e, he, hE, hee₁⟩ := hC
    by_cases hD : ∃ e, e < g.ne ∧ ¬E e ∧ e ≠ e₂
    · obtain ⟨e', he', hE', hee₂⟩ := hD
      exact .inl ⟨huv, .inr ⟨e₁, e, e₂, e', he₁, he, he₂, he', hee₁.symm, hA e he hE,
        hee₂.symm, hB e' he' hE', n12⟩⟩
    · push Not at hD
      refine .inr (.inr ⟨e₂, h.compl.joins_of_single h2 he₂ hE₂ hD he₁ (not_not.2 hE₁), ?_⟩)
      exact fun e he => ⟨hD e he, fun h' => h' ▸ hE₂⟩
  · push Not at hC
    refine .inr (.inl ⟨e₁, h.joins_of_single h2 he₁ hE₁ hC he₂ hE₂, ?_⟩)
    exact fun e he => ⟨hC e he, fun h' => h' ▸ hE₁⟩

/-- Packaging for the walk invariant: a nonempty proper 2-attached `E` of a block is a union of
separation classes of its terminals `{u, v}`, which are distinct and both touched by `E`, and
`{u, v}` is a separation pair unless `E` or its complement is a single `u–v` edge. -/
theorem twoAttached_union_classes (h2 : g.TwoConnected) {u v : Nat} (h : g.TwoAttached E u v)
    (hne : ∃ e, e < g.ne ∧ E e) (hpr : ∃ e, e < g.ne ∧ ¬E e) :
    (∀ e e', g.SepClass u v e e' → (E e ↔ E e')) ∧ u ≠ v ∧ g.Touches E u ∧ g.Touches E v ∧
      (g.SeparationPair u v ∨ (∃ e₀, g.Joins e₀ u v ∧ ∀ e, e < g.ne → (E e ↔ e = e₀)) ∨
        (∃ e₀, g.Joins e₀ u v ∧ ∀ e, e < g.ne → (¬E e ↔ e = e₀))) := by
  obtain ⟨e₁, he₁, hE₁⟩ := hne
  obtain ⟨e₂, he₂, hE₂⟩ := hpr
  have ht := h.touches h2 he₁ hE₁ he₂ hE₂
  exact ⟨fun _ _ hc => h.sepClass_mem hc, h.ne h2 he₁ hE₁ he₂ hE₂, ht.1, ht.2,
    h.sepPair_or_single h2 ⟨e₁, he₁, hE₁⟩ ⟨e₂, he₂, hE₂⟩⟩

end Graph

namespace DfsData

variable {g : Graph} {d : DfsData} {E : Nat → Prop}

/-- Over a block with its sorted DFS tree, the terminals of a nonempty proper 2-attached edge set
form a type-1 or type-2 pair (in one of the two orders), unless the set or its complement is a
single `u–v` edge. -/
theorem twoAttached_type1_or_type2 (hs : d.Spec g) (hr : d.Rooted g) (h2 : g.TwoConnected)
    {u v : Nat} (h : g.TwoAttached E u v) (hne : ∃ e, e < g.ne ∧ E e)
    (hpr : ∃ e, e < g.ne ∧ ¬E e) :
    ((d.Type1Pair u v g ∨ d.Type2Pair u v) ∨ (d.Type1Pair v u g ∨ d.Type2Pair v u)) ∨
      (∃ e₀, g.Joins e₀ u v ∧ ∀ e, e < g.ne → (E e ↔ e = e₀)) ∨
      (∃ e₀, g.Joins e₀ u v ∧ ∀ e, e < g.ne → (¬E e ↔ e = e₀)) := by
  obtain ⟨-, -, hu, hv, hsp | hsp⟩ := Graph.twoAttached_union_classes h2 h hne hpr
  · obtain ⟨eu, -, hu⟩ := hu.isEnd
    obtain ⟨ev, -, hv⟩ := hv.isEnd
    exact .inl ((sepPair_iff' hs hr h2 (hr _ _ hu) (hr _ _ hv)).mp hsp)
  · exact .inr hsp

end DfsData

end Spqr
