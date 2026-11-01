import Spqr.Proofs.RSkel3

/-!
# The R-skeleton invariant of the walk (PROOF.md §4.5, item level)

`Items.RSkelInv g items`: every R item has a 3-connected skeleton (`Items.RSkel3`). The walk only
changes items by pushing childless items and by modifying parentless ones (the item being closed,
the pending edge's Q item, the current vertex's V item), so the invariant is kept by every step
other than an R close, where the new R item is `RBranch.rSkel3`.
-/

namespace Spqr

def Items.RSkelInv (g : Graph) (items : Items) : Prop :=
  ∀ i, i < items.size → items.type i = .R → Items.RSkel3 g items i

theorem Items.not_below_of_ne_root {items : Items} {i j : ItemId} (hij : i ≠ j)
    (hj : ∀ p, ¬ items.IsParent p j) : ¬ items.Below i j :=
  fun h => hij (Items.Below.eq_of_no_parent hj h)

theorem Items.RSkelInv.push_nil {g : Graph} {items : Items} (h : Items.RSkelInv g items)
    (hc : ∀ i, i < items.size → ∀ c ∈ items.ch i, c < items.size) (x : Item) (hx : x.ch = [])
    (hR : x.type ≠ .R) : Items.RSkelInv g (items.push x) := by
  intro i hi hty
  rw [Array.size_push] at hi
  rcases Nat.lt_or_eq_of_le (Nat.le_of_lt_succ hi) with hlt | rfl
  · rw [Items.type_push_of_ne x (Nat.ne_of_lt hlt)] at hty
    exact (h i hlt hty).push_nil hlt (hc i hlt) x hx
  · rw [Items.type_push_size] at hty
    exact absurd hty hR

theorem Items.RSkelInv.modify_root {g : Graph} {items : Items} (h : Items.RSkelInv g items) {j : ItemId}
    (hj : ∀ p, ¬ items.IsParent p j) (f : Item → Item)
    (hnew : Items.type (items.modify j f) j = .R → Items.RSkel3 g (items.modify j f) j) :
    Items.RSkelInv g (items.modify j f) := by
  intro i hi hty
  rw [Array.size_modify] at hi
  by_cases hij : i = j
  · subst hij; exact hnew hty
  · rw [Items.type_modify_of_ne j f hij] at hty
    exact (h i hi hty).modify_of_not_below (Items.not_below_of_ne_root hij hj) f

/-- A parentless non-R item may be modified freely. -/
theorem Items.RSkelInv.modify_root_of_ne {g : Graph} {items : Items} (h : Items.RSkelInv g items) {j : ItemId}
    (hj : ∀ p, ¬ items.IsParent p j) (f : Item → Item) (hty : ∀ it, (f it).type = it.type)
    (hR : items.type j ≠ .R) : Items.RSkelInv g (items.modify j f) :=
  h.modify_root hj f fun hR' => by
    by_cases hjs : j < items.size
    · simp [Items.type, Array.getElem_modify, hjs, hty] at hR'
      exact (hR (by simp [Items.type, Array.getElem?_eq_getElem hjs, hR'])).elim
    · rw [Items.type, Array.getElem?_modify, Array.getElem?_eq_none (Nat.le_of_not_lt hjs)] at hR'
      simp at hR'

namespace WalkState

variable {s : WalkState} {dfs : DfsData} {d : Nat} {cur nxt : TEntry} {rest : List TEntry}

/-- The R close of Loop 1 keeps the invariant: the new item is `RBranch.rSkel3`, the old ones are
untouched (the fresh index is parentless). -/
theorem rCloseItems_rSkelInv (h : Items.RSkelInv s.g s.items) (hsh : Shape s)
    (hb : s.RBranch d cur nxt rest) (hR : s.RTop dfs cur nxt) (hi : s.Inv' (d + 1))
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g) (dir : Bool)
    (hside : getSide (TEntry.mergeInto cur nxt).spans (!dir) = []) :
    Items.RSkelInv s.g (s.rCloseItems d cur nxt dir) := by
  intro i hi' hty
  have hsz : (s.rCloseItems d cur nxt dir).size = s.items.size + 1 := by
    simp [rCloseItems]
  rw [hsz] at hi'
  rcases Nat.lt_or_eq_of_le (Nat.le_of_lt_succ hi') with hlt | rfl
  · have hne := Nat.ne_of_lt hlt
    unfold rCloseItems at hty ⊢
    rw [Items.type_modify_of_ne _ _ hne, Items.type_push_of_ne _ hne] at hty
    refine ((h i hlt hty).push_nil hlt (fun c hc => hsh.ch_lt i c hc) _ rfl).modify_of_not_below ?_ _
    rw [Items.Below_push_nil _ rfl]
    exact Items.not_below_of_ne_root hne hsh.no_parent_size
  · exact RBranch.rSkel3 hsh hb hR hi h2 hsp hrt dir hside

/-- `loop1Type` answers `.R` without changing the state. -/
theorem l1S₁_of_R {edgeDir : Bool} (hty : l1Ty d edgeDir s = .R) : l1S₁ d edgeDir s = s := by
  unfold l1S₁ after; have h := hty; unfold l1Ty result at h; rw [loop1Type_run] at h ⊢
  by_cases h1 : (nxtE s).topDepth > d
  · simp [h1] at h
  · by_cases h2 : ((nxtE s).vStart == (curE s).vStart) = true
    · simp [h1, h2] at h
    · simp [h1, h2]

/-- The R iterate of Loop 1 as a state: `allocItem .R`, merge, close (`rCloseItems`). -/
theorem loop1Body_run_of_R {edgeDir : Bool} (hts : s.tstack = cur :: nxt :: rest)
    (hcd : cur.topDepth = d) (hnd : nxt.topDepth = d) (hty : l1Ty d edgeDir s = .R) :
    ((loop1Body d edgeDir).run s).2 =
      { s with
        items := s.rCloseItems d cur nxt (s.stackDir[d]!)
        tstack := { (TEntry.mergeInto cur nxt) with
          spans := setSides (s.stackDir[d]!) [s.items.size] [] } :: rest } := by
  rw [loop1Body_run_eq, l1S₁_of_R hty, hty]
  unfold closeAt result after
  rw [maybeUnwrapNxt_run_eq .R s cur nxt rest hts _ rfl _ rfl]
  simp only [true_or, ↓reduceIte, run_allocItem]
  have e1 := mergeTstackTops_run_eq { s with items := s.items.push ⟨.R, (none, none), []⟩ } cur nxt rest hts
  rw [e1]
  have e2 := finishTstackTop_run_eq
    { s with
      items := s.items.push ⟨.R, (none, none), []⟩
      tstack := TEntry.mergeInto cur nxt :: rest } s.items.size (TEntry.mergeInto cur nxt) rest rfl
  rw [e2]
  have htd : (TEntry.mergeInto cur nxt).topDepth = d := by
    simp [TEntry.mergeInto, hcd, hnd]
  simp only [rCloseItems, htd]
  rfl

end WalkState

end Spqr
