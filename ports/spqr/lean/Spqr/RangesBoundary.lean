import Spqr.RangesStep
import Spqr.WalkInv

set_option maxHeartbeats 2000000

namespace Spqr
open WalkM
namespace WalkState
variable {s : WalkState} {σ : List Nat} {n D : Nat}

theorem RangesInv.advance (h : s.RangesInv σ n D) {m : Nat} (hn : n ≤ m) :
    s.RangesInv σ m D :=
  ⟨h.inv, fun t ht e he hte => Nat.lt_of_lt_of_le (h.processed t ht e he hte) hn,
    h.ordered, h.convex, h.closed⟩

theorem RangesInv.modifyVs_leaf (j : ItemId) (vs : Option Nat × Option Nat)
    (h : s.RangesInv σ n D) (hj : Items.ch s.items j = []) :
    RangesInv σ n D { s with items := s.items.modify j fun it => { it with vs := vs } } := by
  have hch := Items.ch_modify_ch_eq j (fun it => { it with vs := vs }) (fun _ => rfl) (items := s.items)
  have hty := Items.type_modify_type_eq j (fun it => { it with vs := vs }) (fun _ => rfl) (items := s.items)
  have hB : ∀ a i, Items.Below (s.items.modify j fun it => { it with vs := vs }) a i ↔ Items.Below s.items a i :=
    fun _ _ => Items.Below_modify_ch_eq j (fun it => { it with vs := vs }) fun _ => rfl
  have hN : ∀ a i, Items.BelowNoV (s.items.modify j fun it => { it with vs := vs }) a i ↔ Items.BelowNoV s.items a i :=
    fun _ _ => Items.BelowNoV_congr hch (fun _ c _ => hty c)
  refine h.items_congr (h.inv.modifyVs_leaf j vs hj)
    (fun t _ e => TEntry.edges_congr (fun i _ _ => hB i _) e)
    (fun t _ e => TEntry.piece_congr (fun i _ => hty i) (fun i _ _ => hN i _) e) ?_
  intro i hi ht a b c hab hbc hc ha hz
  exact (hB i _).2 (h.closed i (by simpa using hi) (by rwa [hty] at ht)
    a b c hab hbc hc ((hN i _).1 ha) ((hN i _).1 hz))

theorem RangesInv.pop' (h : s.RangesInv σ n D)
    (hi : Inv' D { s with tstack := s.tstack.tail }) :
    RangesInv σ n D { s with tstack := s.tstack.tail } := by
  refine ⟨hi, fun t ht => h.processed t (List.mem_of_mem_tail ht), ?_,
    fun t ht => h.convex t (List.mem_of_mem_tail ht), h.closed⟩
  intro above t below hts
  cases hs : s.tstack with
  | nil => simp [hs] at hts
  | cons x rest =>
    have hr : rest = above ++ t :: below := by simpa [hs] using hts
    exact h.ordered (x :: above) t below (by rw [hs, hr]; rfl)

/-- Writing a free root preserves the old ranges; only that root's new range needs checking. -/
theorem RangesInv.modifyCh (h : s.RangesInv σ n D) (j : ItemId) (f : Item → Item)
    (hj : j < 1 + s.g.nv + s.g.ne) (hty : ∀ it, (f it).type = it.type)
    (hroot : ∀ p, ¬ Items.IsParent s.items p j)
    (hfree : ∀ t ∈ s.tstack, j ∉ t.spans.1 ++ t.spans.2)
    (hclosed : Items.type s.items j ∉ [NodeType.F, .V] →
      ∀ a b c, a ≤ b → b ≤ c → c < σ.length →
        Items.PieceEdge s.g (s.items.modify j f) j σ[a]! →
        Items.PieceEdge s.g (s.items.modify j f) j σ[c]! →
        Items.EdgeBelow s.g (s.items.modify j f) j σ[b]!) :
    RangesInv σ n D { s with items := s.items.modify j f } := by
  have hnb : ∀ i, i ≠ j → ¬ Items.Below s.items i j :=
    fun i hi hb => hi (Items.Below.eq_of_no_parent hroot hb)
  have hT := Items.type_modify_type_eq j f hty (items := s.items)
  refine h.items_congr (h.inv.modifyCh j f hj hroot hfree)
    (fun t ht e => TEntry.edges_modify_of_not_mem j f hroot (hfree t ht) e)
    (fun t ht e => TEntry.piece_congr (fun i _ => hT i)
      (fun i hi _ => Items.BelowNoV_modify_of_not_below j f hty
        (hnb i (fun heq => hfree t ht (heq ▸ hi)))) e) ?_
  intro i hi ht a b c hab hbc hc ha hz
  rw [hT] at ht
  by_cases heq : i = j
  · subst i; exact hclosed ht a b c hab hbc hc ha hz
  · exact (Items.Below_modify_of_not_below j f (hnb i heq)).2
      (h.closed i (by simpa using hi) ht a b c hab hbc hc
        ((Items.BelowNoV_modify_of_not_below j f hty (hnb i heq)).1 ha)
        ((Items.BelowNoV_modify_of_not_below j f hty (hnb i heq)).1 hz))

end WalkState

namespace Items
variable {g : Graph} {items : Items} {σ : List Nat} {n e : Nat} {xs : List ItemId}

/-- The children of a boundary cap fill the postorder up to its pending edge. -/
structure CapRange (σ : List Nat) (n : Nat) (g : Graph) (items : Items) (xs : List ItemId) : Prop where
  processed : ∀ i ∈ xs, items.type i ≠ .V → ∀ a, a < σ.length →
    items.PieceEdge g i σ[a]! → a < n
  fill : ∀ a b, a ≤ b → b < n →
    (∃ i ∈ xs, items.type i ≠ .V ∧ items.PieceEdge g i σ[a]!) →
    ∃ i ∈ xs, items.EdgeBelow g i σ[b]!

theorem CapRange.congr {items' : Items} (h : CapRange σ n g items xs)
    (hty : ∀ i ∈ xs, items'.type i = items.type i)
    (hB : ∀ i ∈ xs, ∀ e, items'.EdgeBelow g i e ↔ items.EdgeBelow g i e)
    (hN : ∀ i ∈ xs, ∀ e, items'.PieceEdge g i e ↔ items.PieceEdge g i e) :
    CapRange σ n g items' xs := by
  refine ⟨fun i hi ht a ha hp => h.processed i hi (by rwa [hty i hi] at ht) a ha
    ((hN i hi _).1 hp), ?_⟩
  rintro a b hab hb ⟨i, hi, ht, hp⟩
  obtain ⟨j, hj, he⟩ := h.fill a b hab hb ⟨i, hi, by rwa [hty i hi] at ht, (hN i hi _).1 hp⟩
  exact ⟨j, hj, (hB j hj _).2 he⟩

theorem CapRange.modifyVs (h : CapRange σ n g items xs) (j : ItemId) (vs : Option Nat × Option Nat) :
    CapRange σ n g (items.modify j fun it => { it with vs := vs }) xs :=
  h.congr (fun i _ => type_modify_type_eq j (fun it => { it with vs := vs }) (fun _ => rfl) i)
    (fun _ _ _ => Below_modify_ch_eq j (fun it => { it with vs := vs }) fun _ => rfl)
    (fun _ _ _ => BelowNoV_congr (ch_modify_ch_eq j (fun it => { it with vs := vs }) fun _ => rfl)
      (fun _ c _ => type_modify_type_eq j (fun it => { it with vs := vs }) (fun _ => rfl) c))

theorem CapRange.alloc (h : CapRange σ n g items xs) (ty : NodeType)
    (hxs : ∀ i ∈ xs, i < items.size)
    (hch : ∀ p c, items.IsParent p c → c < items.size) :
    CapRange σ n g (items.push ⟨ty, (none, none), []⟩) xs :=
  h.congr (fun i hi => type_push_of_ne _ (Nat.ne_of_lt (hxs i hi)))
    (fun _ _ _ => Below_push_nil _ rfl)
    (fun _ _ _ => BelowNoV_congr (ch_push_nil _ rfl)
      (fun p c hpc => type_push_of_ne _ (Nat.ne_of_lt (hch p c hpc))))

theorem CapRange.cons_leaf (h : CapRange σ n g items xs) (j : ItemId)
    (hj : g.ne + (1 + g.nv) ≤ j) (hleaf : items.ch j = [])
    (hσ : ∀ e ∈ σ, e < g.ne) (hn : n ≤ σ.length) : CapRange σ n g items (j :: xs) := by
  have hempty : ∀ k, k < σ.length → ¬ items.PieceEdge g j σ[k]! := by
    intro k hk hp
    have he := hσ σ[k]! (by rw [getElem!_pos σ k hk]; exact List.getElem_mem hk)
    rcases hp.below.head_cases with heq | ⟨c, hc, _⟩
    · have heq : j = 1 + g.nv + σ[k]! := heq
      rw [heq] at hj
      omega
    · simp [IsParent, hleaf] at hc
  refine ⟨?_, ?_⟩
  · intro i hi ht a ha hp
    rcases List.mem_cons.1 hi with rfl | hi
    · exact (hempty a ha hp).elim
    · exact h.processed i hi ht a ha hp
  · rintro a b hab hb ⟨i, hi, ht, hp⟩
    rcases List.mem_cons.1 hi with rfl | hi
    · exact (hempty a (by omega) hp).elim
    · obtain ⟨k, hk, he⟩ := h.fill a b hab hb ⟨i, hi, ht, hp⟩
      exact ⟨k, List.mem_cons_of_mem _ hk, he⟩

theorem cap_convex (hr : CapRange σ n g items xs) (hnd : σ.Nodup)
    (hpos : σ[n]? = some e) (hj : edgeItem g e < items.size)
    (hroot : ∀ p, ¬ items.IsParent p (edgeItem g e)) (hfree : edgeItem g e ∉ xs) :
    ∀ a b c, a ≤ b → b ≤ c → c < σ.length →
      PieceEdge g (items.modify (edgeItem g e) fun it => { it with ch := xs }) (edgeItem g e) σ[a]! →
      PieceEdge g (items.modify (edgeItem g e) fun it => { it with ch := xs }) (edgeItem g e) σ[c]! →
      EdgeBelow g (items.modify (edgeItem g e) fun it => { it with ch := xs }) (edgeItem g e) σ[b]! := by
  let j := edgeItem g e
  let f : Item → Item := fun it => { it with ch := xs }
  have hch : Items.ch (items.modify j f) j = xs := by rw [ch_modify_at j f hj]
  have hty := type_modify_type_eq j f (fun _ => rfl) (items := items)
  have hnb : ∀ i ∈ xs, ¬ items.Below i j := fun i hi hb =>
    hfree (Items.Below.eq_of_no_parent hroot hb ▸ hi)
  have hn : n < σ.length := List.getElem?_eq_some_iff.1 hpos |>.1
  have hnval : σ[n]! = e := by rw [getElem!_pos σ n hn]; exact (List.getElem?_eq_some_iff.1 hpos).2
  have piece : ∀ k, k < σ.length → PieceEdge g (items.modify j f) j σ[k]! →
      k = n ∨ ∃ i ∈ xs, items.type i ≠ .V ∧ items.PieceEdge g i σ[k]! := by
    intro k hk hp
    rcases Relation.ReflTransGen.cases_head hp with heq | ⟨i, ⟨hi, ht⟩, hp⟩
    · have he : σ[k]! = e := (edgeItem_inj heq).symm
      have := WalkState.idxOf_getElem! hnd hk
      have := WalkState.idxOf_getElem! hnd hn
      rw [he] at *; rw [hnval] at *; exact .inl (by omega)
    · have hi : i ∈ xs := by simpa [Items.IsParent, hch] using hi
      exact .inr ⟨i, hi, by rwa [hty] at ht,
        (BelowNoV_modify_of_not_below j f (fun _ => rfl) (hnb i hi)).1 hp⟩
  intro a b c hab hbc hc ha hz
  have ha' := piece a (by omega) ha
  have hc' : c ≤ n := by
    rcases piece c hc hz with rfl | ⟨i, hi, ht, hp⟩
    · exact Nat.le_refl _
    · exact Nat.le_of_lt (hr.processed i hi ht c hc hp)
  by_cases hbn : b = n
  · subst b; rw [hnval]; exact .refl
  · have hb : b < n := by omega
    have ha' : ∃ i ∈ xs, items.type i ≠ .V ∧ items.PieceEdge g i σ[a]! := by
      rcases ha' with ha' | ha'
      · omega
      · exact ha'
    obtain ⟨i, hi, hp⟩ := hr.fill a b hab hb ha'
    exact .head (by change i ∈ Items.ch (items.modify j f) j; rw [hch]; exact hi)
      ((Below_modify_of_not_below j f (hnb i hi)).2 hp)

end Items

namespace WalkState

theorem RangesInv.cap_of_entries {s : WalkState} {σ : List Nat} {n D : Nat}
    (h : s.RangesInv σ n D) (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne)
    (ts : List TEntry) (xs : List ItemId) (hsub : ∀ t ∈ ts, t ∈ s.tstack)
    (hxs : ∀ i, i ∈ xs ↔ ∃ t ∈ ts, i ∈ t.spans.1 ++ t.spans.2)
    (hfill : ∀ a b, a ≤ b → b < n → (∃ t ∈ ts, t.piece s.g s.items σ[a]!) →
      ∃ t ∈ ts, t.edges s.g s.items σ[b]!) : Items.CapRange σ n s.g s.items xs := by
  refine ⟨?_, ?_⟩
  · intro i hi ht a ha hp
    obtain ⟨t, htt, hit⟩ := (hxs i).1 hi
    have hh := h.processed t (hsub t htt) σ[a]! (getElem!_lt hσ ha) ⟨i, hit, hp.edgeBelow⟩
    rwa [idxOf_getElem! hnd ha] at hh
  · rintro a b hab hb ⟨i, hi, ht, hp⟩
    obtain ⟨t, htt, hit⟩ := (hxs i).1 hi
    obtain ⟨u, hu, j, hj, he⟩ := hfill a b hab hb ⟨t, htt, i, hit, ht, hp⟩
    exact ⟨j, (hxs j).2 ⟨u, hu, hj⟩, he⟩

end WalkState
end Spqr
