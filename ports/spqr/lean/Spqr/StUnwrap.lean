import Spqr.StClose

/-!
# `StSim` through `maybeUnwrapNxt; mergeTstackTops; finishTstackTop`

The common tail of a loop-1 iteration (`loop1Body`) and of the P-check (`finishP`): the entry
under the top is reopened (when it is a single node of the wanted type) or a fresh node is
allocated, the two top entries merge, and the merged entry closes into the node.
-/

namespace Spqr
open WalkM

/-! ### Reading helpers -/

theorem readL_append (l l' : List TEntry) : readL (l ++ l') = readL l ++ readL l' := by
  induction l with
  | nil => rfl
  | cons t l ih => simp [readL, ih]

theorem readR_append (l l' : List TEntry) : readR (l ++ l') = readR l' ++ readR l := by
  induction l with
  | nil => simp [readR]
  | cons t l ih => simp [readR, ih]

theorem mem_readStack_append_left {x : ItemId} {l l' : List TEntry} (h : x ∈ readStack l) :
    x ∈ readStack (l ++ l') := by
  simp only [readStack, readL_append, readR_append, List.mem_append] at h ⊢
  tauto

/-- Reading `u :: l`, as a permutation. -/
theorem readStack_cons_perm' (u : TEntry) (l : List TEntry) :
    (readStack (u :: l)).Perm (u.spans.1 ++ u.spans.2 ++ readStack l) := by
  rw [List.perm_iff_count]; intro a
  simp [readStack, readL, readR, List.count_append]; omega

theorem mem_readStack_cons {x : ItemId} {u : TEntry} {l : List TEntry} :
    x ∈ readStack (u :: l) ↔ x ∈ u.spans.1 ++ u.spans.2 ∨ x ∈ readStack l := by
  rw [(readStack_cons_perm' u l).mem_iff, List.mem_append]

theorem mem_expandItem {i : ItemId} {ch l : List ItemId} {x : ItemId} :
    x ∈ expandItem i ch l ↔ (x ∈ l ∧ x ≠ i) ∨ (i ∈ l ∧ x ∈ ch) := by
  simp only [expandItem, List.mem_flatMap]
  constructor
  · rintro ⟨a, ha, hx⟩
    by_cases h : a = i
    · subst h; simp at hx; exact Or.inr ⟨ha, hx⟩
    · simp [h] at hx; subst hx; exact Or.inl ⟨ha, h⟩
  · rintro (⟨hx, hne⟩ | ⟨hi, hx⟩)
    · exact ⟨x, hx, by simp [hne]⟩
    · exact ⟨i, hi, by simp [hx]⟩

theorem expandItem_cons (i : ItemId) (ch : List ItemId) (y : ItemId) (l : List ItemId) :
    expandItem i ch (y :: l) = (if y = i then ch else [y]) ++ expandItem i ch l := by
  simp [expandItem]

theorem nodup_expandItem {i : ItemId} {ch l : List ItemId} (hl : l.Nodup) (hch : ch.Nodup)
    (hd : ∀ x ∈ ch, x ∉ l) : (expandItem i ch l).Nodup := by
  induction l with
  | nil => simp [expandItem]
  | cons y l ih =>
    rw [List.nodup_cons] at hl
    have ih := ih hl.2 (fun x hx hxl => hd x hx (List.mem_cons_of_mem _ hxl))
    rw [expandItem_cons, List.nodup_append]
    by_cases hy : y = i
    · subst hy
      rw [if_pos rfl]
      refine ⟨hch, ih, fun x hx x' hx' hxx' => ?_⟩
      subst hxx'
      rcases mem_expandItem.1 hx' with ⟨hxl, _⟩ | ⟨hyl, _⟩
      · exact hd x hx (List.mem_cons_of_mem _ hxl)
      · exact hl.1 hyl
    · rw [ite_eq_right_iff.2 (fun e => absurd e hy)]
      refine ⟨List.nodup_singleton _, ih, fun x hx x' hx' hxx' => ?_⟩
      simp at hx; subst hx; subst hxx'
      rcases mem_expandItem.1 hx' with ⟨hxl, _⟩ | ⟨_, hxc⟩
      · exact hl.1 hxl
      · exact hd _ hxc List.mem_cons_self

/-! ### `mergeTstackTops` -/

theorem mergeTops_cons_cons (c t : TEntry) (R : List TEntry) :
    mergeTops (c :: t :: R) = TEntry.mergeInto c t :: R := rfl

theorem readStack_mergeInto_cons (c t : TEntry) (R : List TEntry) :
    readStack (TEntry.mergeInto c t :: R) = readStack (c :: t :: R) :=
  readStack_mergeTops c t R

theorem readL_mergeInto_cons (c t : TEntry) (R : List TEntry) :
    readL (TEntry.mergeInto c t :: R) = readL (c :: t :: R) := by
  simp [readL, TEntry.mergeInto]

theorem readR_mergeInto_cons (c t : TEntry) (R : List TEntry) :
    readR (TEntry.mergeInto c t :: R) = readR (c :: t :: R) := by
  simp [readR, TEntry.mergeInto]

theorem TEntry.mergeInto_side_nil (dir : Bool) (c t : TEntry) (hc : getSide c.spans dir = [])
    (ht : getSide t.spans dir = []) : getSide (TEntry.mergeInto c t).spans dir = [] := by
  cases dir <;> simp_all [TEntry.mergeInto, getSide]

theorem run_mergeTstackTops_cons_cons (s : WalkState) (c t : TEntry) (R : List TEntry)
    (hts : s.tstack = c :: t :: R) :
    (mergeTstackTops.run s).2 = { s with tstack := TEntry.mergeInto c t :: R } := by
  rw [run_mergeTstackTops, hts, mergeTops_cons_cons]

/-! ### The allocation path -/

/-- `allocItem ty` (`item = size`), `mergeTstackTops`, `finishTstackTop item` on `c :: t :: _`
(both one-sided on `stackDir[t.topDepth]`, which is also the direction at the merged depth). -/
theorem StSim.allocMergeClose {g : Graph} (s : WalkState) (ty : NodeType) (c t : TEntry)
    (new base : List TEntry) (ps : List StPiece) (blocks : List StBlock)
    (hts : s.tstack = c :: t :: (new ++ base))
    (hc : getSide c.spans (!s.stackDir[min t.topDepth c.topDepth]!) = [])
    (ht : getSide t.spans (!s.stackDir[min t.topDepth c.topDepth]!) = [])
    (hty : ty ≠ .V ∧ ty ≠ .Q)
    (hR : StRead s.items (c :: t :: new) ps) (hI : StItems g s blocks) :
    let s' := ((finishTstackTop s.items.size).run
      (mergeTstackTops.run { s with items := s.items.push ⟨ty, (none, none), []⟩ }).2).2
    s'.tstack = { TEntry.mergeInto c t with spans := setSides s.stackDir[min t.topDepth c.topDepth]! [s.items.size] [] } ::
        (new ++ base) ∧
      s'.stackDir = s.stackDir ∧ s'.items.size = s.items.size + 1 ∧
      Items.type s'.items s.items.size = ty ∧
      StRead s'.items ({ TEntry.mergeInto c t with spans := setSides s.stackDir[min t.topDepth c.topDepth]! [s.items.size] [] } :: new) ps ∧
      StItems g s' blocks := by
  dsimp only
  have hs₂ := run_mergeTstackTops_cons_cons { s with items := s.items.push ⟨ty, (none, none), []⟩ } c t
    (new ++ base) hts
  generalize (mergeTstackTops.run { s with items := s.items.push ⟨ty, (none, none), []⟩ }).2 = s₂ at hs₂ ⊢
  have hts₂ : s₂.tstack = TEntry.mergeInto c t :: (new ++ base) := by rw [hs₂]
  have hitems₂ : s₂.items = s.items.push ⟨ty, (none, none), []⟩ := by rw [hs₂]
  have hsd₂ : s₂.stackDir = s.stackDir := by rw [hs₂]
  have hread₂ : readStack s₂.tstack = readStack s.tstack := by
    rw [hts₂, readStack_mergeInto_cons, hts]
  have hsz : s₂.items.size = s.items.size + 1 := by rw [hitems₂, Array.size_push]
  have htop : (TEntry.mergeInto c t).topDepth = min t.topDepth c.topDepth := rfl
  have hside : getSide (TEntry.mergeInto c t).spans (!s₂.stackDir[(TEntry.mergeInto c t).topDepth]!) = [] := by
    rw [hsd₂, htop]; exact TEntry.mergeInto_side_nil _ c t hc ht
  have hdir₂ : s₂.stackDir[(TEntry.mergeInto c t).topDepth]! = s.stackDir[min t.topDepth c.topDepth]! := by
    rw [hsd₂, htop]
  have hlt : ∀ x ∈ readStack s.tstack, x < s.items.size := fun x hx => hI.bounded x hx x .refl
  have hsub : ∀ x ∈ readStack (c :: t :: new), x ∈ readStack s.tstack := fun x hx => by
    rw [hts]; exact mem_readStack_append_left (l' := base) hx
  have hnotin : s.items.size ∉ readStack s.tstack := fun h => Nat.lt_irrefl _ (hlt _ h)
  have hparent₂ : ∀ p x, Items.IsParent s₂.items p x ↔ p ≠ s.items.size ∧ Items.IsParent s.items p x := by
    intro p x
    rw [hitems₂]
    unfold Items.IsParent
    rw [Items.ch_push]
    by_cases hp : p = s.items.size <;> simp [hp]
  have hbelow₂ : ∀ x y, x < s.items.size → Items.Below s₂.items x y → Items.Below s.items x y := by
    intro x y hx h
    rw [hitems₂] at h
    exact h.of_push _ hI.chLt hx
  have hbelow₂' : ∀ x y, x ∈ readStack s.tstack → Items.Below s.items x y → Items.Below s₂.items x y := by
    intro x y hx h
    rw [hitems₂]
    exact h.push _ (hI.bounded x hx)
  have htype₂ : ∀ j, j ≠ s.items.size → Items.type s₂.items j = Items.type s.items j := by
    intro j hj; rw [hitems₂, Items.type_push, ite_eq_right_iff.2 (fun e => absurd e hj)]
  have hR₂ : StRead s₂.items (TEntry.mergeInto c t :: new) ps := by
    unfold StRead at hR ⊢
    rw [readL_mergeInto_cons, readR_mergeInto_cons, hitems₂]
    exact ⟨hR.1.push _ (fun x hx y hy => hI.bounded x (hsub x (mem_readStack_of_readL hx)) y hy),
      hR.2.push _ (fun x hx y hy => hI.bounded x (hsub x (mem_readStack_of_readR hx)) y hy)⟩
  have hnb : ∀ x ∈ readStack (TEntry.mergeInto c t :: new), ¬ Items.Below s₂.items x s.items.size := by
    intro x hx h
    rw [readStack_mergeInto_cons] at hx
    exact Nat.lt_irrefl _ (hI.bounded x (hsub x hx) _ (hbelow₂ x _ (hlt x (hsub x hx)) h))
  have hnew : s.items.size ∉ readStack new := fun h =>
    hnotin (hsub _ (mem_readStack_cons.2 (Or.inr (mem_readStack_cons.2 (Or.inr h)))))
  have hty₂ : ¬ (Items.type s₂.items s.items.size = .V ∨ Items.type s₂.items s.items.size = .Q) := by
    rw [hitems₂, Items.type_push_size]
    rintro (h | h)
    · exact hty.1 h
    · exact hty.2 h
  have hi₂ : s.items.size < s₂.items.size := by rw [hsz]; exact Nat.lt_succ_self _
  obtain ⟨hts', hR'⟩ := StRead.finishTstackTop s₂ s.items.size (TEntry.mergeInto c t) new base ps hts₂
    hside hi₂ hty₂ hnb hnew hR₂
  have hitems' := finishTstackTop_items s₂ s.items.size (TEntry.mergeInto c t) (new ++ base) hts₂
  rw [hdir₂] at hts' hR'
  refine ⟨hts', ?_, ?_, ?_, hR', ?_⟩
  · show ((finishTstackTop s.items.size).run s₂).2.stackDir = s.stackDir
    rw [← hsd₂]
    obtain ⟨x', items', h, -⟩ := finishTstackTop_run s.items.size s₂ hts₂
    rw [h]
  · show ((finishTstackTop s.items.size).run s₂).2.items.size = s.items.size + 1
    rw [hitems', Array.size_modify, hsz]
  · show Items.type ((finishTstackTop s.items.size).run s₂).2.items s.items.size = ty
    rw [hitems', Items.type_modify s₂.items s.items.size s.items.size
      (fun it => { it with
        vs := setSides s₂.stackDir[(TEntry.mergeInto c t).topDepth]!
          (some (TEntry.top s₂ (TEntry.mergeInto c t))) (some (TEntry.mergeInto c t).vStart),
        ch := getSide (TEntry.mergeInto c t).spans s₂.stackDir[(TEntry.mergeInto c t).topDepth]! })
      (fun _ => rfl), hitems₂, Items.type_push_size]
  · refine StItems.close s₂ blocks s.items.size (TEntry.mergeInto c t) (new ++ base) hts₂ hside hi₂
      ?_ ?_ ?_ ?_ ?_ ?_ ?_ ?_
    · intro p hp
      exact Nat.lt_irrefl _ (hI.chLt p _ ((hparent₂ p _).1 hp).2)
    · rw [hread₂]; exact hnotin
    · intro x hx p hpx
      exact hI.roots x (by rw [hts]; exact mem_readStack_cons.2 (Or.inr (mem_readStack_cons.2 (Or.inr hx))))
        p ((hparent₂ p x).1 hpx).2
    · rw [hread₂]; exact hI.nodup
    · intro x hx y hy
      rw [hread₂] at hx
      rw [hsz]
      exact Nat.lt_succ_of_lt (hI.bounded x hx y (hbelow₂ x y (hlt x hx) hy))
    · intro p x hpx
      rw [hsz]
      exact Nat.lt_succ_of_lt (hI.chLt p x ((hparent₂ p x).1 hpx).2)
    · intro p
      rw [hitems₂, Items.ch_push]
      split
      · exact List.nodup_nil
      · exact hI.chNodup p
    · intro j hj hji hty'
      rw [hsz] at hj
      have hj' : j < s.items.size := by omega
      rw [htype₂ j hji] at hty'
      rcases hI.closed j hj' hty' with ⟨x, hx, hxj⟩ | ⟨b, hb, hbj⟩
      · exact Or.inl ⟨x, by rw [hread₂]; exact hx, hbelow₂' x j hx hxj⟩
      · right
        refine ⟨b, hb, ?_⟩
        rw [hitems₂]
        exact hbj.push _ hI.chLt hj'

/-! ### The reopen path -/

/-- Reopening `t = ⟨_, _, _, setSides dir [h] []⟩` (`h` a non-leaf on the stack) to `h`'s children,
`mergeTstackTops`, `finishTstackTop h`. -/
theorem StSim.reopenMergeClose {g : Graph} (s : WalkState) (c t : TEntry) (h : ItemId)
    (new base : List TEntry) (ps : List StPiece) (blocks : List StBlock)
    (hts : s.tstack = c :: t :: (new ++ base))
    (hc : getSide c.spans (!s.stackDir[t.topDepth]!) = [])
    (htsp : t.spans = setSides s.stackDir[t.topDepth]! [h] [])
    (hdir : s.stackDir[min t.topDepth c.topDepth]! = s.stackDir[t.topDepth]!)
    (hty : ¬ (Items.type s.items h = .V ∨ Items.type s.items h = .Q))
    (hR : StRead s.items (c :: t :: new) ps) (hI : StItems g s blocks) :
    let s' := ((finishTstackTop h).run (mergeTstackTops.run { s with tstack :=
      c :: { t with spans := setSides s.stackDir[t.topDepth]! (Items.ch s.items h) [] } :: (new ++ base) }).2).2
    s'.tstack = { TEntry.mergeInto c t with spans := setSides s.stackDir[t.topDepth]! [h] [] } :: (new ++ base) ∧
      s'.stackDir = s.stackDir ∧ s'.items.size = s.items.size ∧
      Items.type s'.items h = Items.type s.items h ∧
      StRead s'.items ({ TEntry.mergeInto c t with spans := setSides s.stackDir[t.topDepth]! [h] [] } :: new) ps ∧
      StItems g s' blocks := by
  dsimp only
  generalize hd : s.stackDir[t.topDepth]! = dir at hc htsp hdir ⊢
  set t' : TEntry := { t with spans := setSides dir (Items.ch s.items h) [] } with ht'
  have hs₂ := run_mergeTstackTops_cons_cons { s with tstack := c :: t' :: (new ++ base) } c t' (new ++ base) rfl
  generalize (mergeTstackTops.run { s with tstack := c :: t' :: (new ++ base) }).2 = s₂ at hs₂ ⊢
  have hts₂ : s₂.tstack = TEntry.mergeInto c t' :: (new ++ base) := by rw [hs₂]
  have hitems₂ : s₂.items = s.items := by rw [hs₂]
  have hsd₂ : s₂.stackDir = s.stackDir := by rw [hs₂]
  have hht : h ∈ t.spans.1 ++ t.spans.2 := by rw [htsp]; cases dir <;> simp [setSides]
  have hhmem : h ∈ readStack s.tstack := by
    rw [hts]; exact mem_readStack_cons.2 (Or.inr (mem_readStack_cons.2 (Or.inl hht)))
  have hh : h < s.items.size := hI.bounded h hhmem h .refl
  have hroot : ∀ p, ¬ Items.IsParent s.items p h := hI.roots h hhmem
  have hchfresh : ∀ x ∈ Items.ch s.items h, x ∉ readStack s.tstack := fun x hx hxs => hI.roots x hxs h hx
  have hhch : h ∉ Items.ch s.items h := fun e => hroot h e
  have hnd : (readStack s.tstack).Nodup := hI.nodup
  have hnd₁ : (c.spans.1 ++ c.spans.2 ++ readStack (t :: (new ++ base))).Nodup := by
    rw [hts] at hnd; exact (readStack_cons_perm' c _).nodup_iff.1 hnd
  have hnd₂ : (t.spans.1 ++ t.spans.2 ++ readStack (new ++ base)).Nodup :=
    (readStack_cons_perm' t _).nodup_iff.1 hnd₁.of_append_right
  have ha : h ∉ c.spans.1 ++ c.spans.2 := fun hc' =>
    List.disjoint_of_nodup_append hnd₁ hc' (mem_readStack_cons.2 (Or.inl hht))
  have hrest : h ∉ readStack (new ++ base) := fun hr => List.disjoint_of_nodup_append hnd₂ hht hr
  have hnew : h ∉ readStack new := fun hn => hrest (mem_readStack_append_left hn)
  have hread₂ : readStack s₂.tstack = expandItem h (Items.ch s.items h) (readStack s.tstack) := by
    rw [hts₂, readStack_mergeInto_cons, hts]
    exact readStack_reopen c t (new ++ base) dir h _ htsp ha hrest
  have hreadseg : readStack (TEntry.mergeInto c t' :: new) =
      expandItem h (Items.ch s.items h) (readStack (c :: t :: new)) := by
    rw [readStack_mergeInto_cons]
    exact readStack_reopen c t new dir h _ htsp ha hnew
  have hreadsegL : readL (TEntry.mergeInto c t' :: new) =
      expandItem h (Items.ch s.items h) (readL (c :: t :: new)) := by
    rw [readL_mergeInto_cons]
    exact readL_reopen c t new dir h _ htsp ha hnew
  have hreadsegR : readR (TEntry.mergeInto c t' :: new) =
      expandItem h (Items.ch s.items h) (readR (c :: t :: new)) := by
    rw [readR_mergeInto_cons]
    exact readR_reopen c t new dir h _ htsp ha hnew
  have hmem₂ : ∀ x, x ∈ readStack s₂.tstack ↔
      (x ∈ readStack s.tstack ∧ x ≠ h) ∨ x ∈ Items.ch s.items h := by
    intro x; rw [hread₂, mem_expandItem]
    constructor
    · rintro (h1 | ⟨_, h2⟩)
      · exact Or.inl h1
      · exact Or.inr h2
    · rintro (h1 | h2)
      · exact Or.inl h1
      · exact Or.inr ⟨hhmem, h2⟩
  have hnotin₂ : h ∉ readStack s₂.tstack := fun hm => by
    rcases (hmem₂ h).1 hm with ⟨_, e⟩ | e
    · exact e rfl
    · exact hhch e
  have htop : (TEntry.mergeInto c t').topDepth = min t.topDepth c.topDepth := rfl
  have hside : getSide (TEntry.mergeInto c t').spans (!s₂.stackDir[(TEntry.mergeInto c t').topDepth]!) = [] := by
    rw [hsd₂, htop, hdir]
    exact TEntry.mergeInto_side_nil _ c t' hc (by rw [ht']; exact getSide_setSides_other _ _ _)
  have hdir₂ : s₂.stackDir[(TEntry.mergeInto c t').topDepth]! = dir := by rw [hsd₂, htop, hdir]
  have hR₂ : StRead s₂.items (TEntry.mergeInto c t' :: new) ps := by
    unfold StRead at hR ⊢
    rw [hitems₂, hreadsegL, hreadsegR]
    exact ⟨(ExpandsList.expandItem_self_iff hty).2 hR.1, (ExpandsList.expandItem_self_iff hty).2 hR.2⟩
  have hnb : ∀ x ∈ readStack (TEntry.mergeInto c t' :: new), ¬ Items.Below s₂.items x h := by
    intro x hx
    rw [hitems₂]
    refine Items.not_below_root_of_ne hroot fun e => ?_
    subst e
    rw [hreadseg, mem_expandItem] at hx
    rcases hx with ⟨_, e⟩ | ⟨_, e⟩
    · exact e rfl
    · exact hhch e
  have hi₂ : h < s₂.items.size := by rw [hitems₂]; exact hh
  have hty₂ : ¬ (Items.type s₂.items h = .V ∨ Items.type s₂.items h = .Q) := by rw [hitems₂]; exact hty
  obtain ⟨hts', hR'⟩ := StRead.finishTstackTop s₂ h (TEntry.mergeInto c t') new base ps hts₂ hside hi₂
    hty₂ hnb hnew hR₂
  have hitems' := finishTstackTop_items s₂ h (TEntry.mergeInto c t') (new ++ base) hts₂
  rw [hdir₂] at hts' hR'
  have hent : ({ TEntry.mergeInto c t' with spans := setSides dir [h] [] } : TEntry) =
      { TEntry.mergeInto c t with spans := setSides dir [h] [] } := rfl
  rw [hent] at hts' hR'
  refine ⟨hts', ?_, ?_, ?_, hR', ?_⟩
  · rw [← hsd₂]
    obtain ⟨x', items', e, -⟩ := finishTstackTop_run h s₂ hts₂
    rw [e]
  · rw [hitems', Array.size_modify, hitems₂]
  · rw [hitems', Items.type_modify s₂.items h h
      (fun it => { it with
        vs := setSides s₂.stackDir[(TEntry.mergeInto c t').topDepth]!
          (some (TEntry.top s₂ (TEntry.mergeInto c t'))) (some (TEntry.mergeInto c t').vStart),
        ch := getSide (TEntry.mergeInto c t').spans s₂.stackDir[(TEntry.mergeInto c t').topDepth]! })
      (fun _ => rfl), hitems₂]
  · refine StItems.close s₂ blocks h (TEntry.mergeInto c t') (new ++ base) hts₂ hside hi₂
      ?_ ?_ ?_ ?_ ?_ ?_ ?_ ?_
    · rw [hitems₂]; exact hroot
    · exact hnotin₂
    · intro x hx p hpx
      rw [hitems₂] at hpx
      exact hI.roots x (by rw [hts]; exact mem_readStack_cons.2 (Or.inr (mem_readStack_cons.2 (Or.inr hx)))) p hpx
    · rw [hread₂]; exact nodup_expandItem hnd (hI.chNodup h) hchfresh
    · intro x hx y hy
      rw [hitems₂] at hy ⊢
      rcases (hmem₂ x).1 hx with ⟨hxs, _⟩ | hxc
      · exact hI.bounded x hxs y hy
      · exact hI.bounded h hhmem y (Relation.ReflTransGen.head (show Items.IsParent s.items h x from hxc) hy)
    · intro p x hpx
      rw [hitems₂] at hpx ⊢
      exact hI.chLt p x hpx
    · intro p; rw [hitems₂]; exact hI.chNodup p
    · intro j hj hjh hty'
      rw [hitems₂] at hj hty' ⊢
      rcases hI.closed j hj hty' with ⟨x, hx, hxj⟩ | ⟨b, hb, hbj⟩
      · left
        by_cases hxh : x = h
        · subst hxh
          rcases Relation.ReflTransGen.cases_head hxj with e | ⟨c', hc', hcj⟩
          · exact absurd e.symm hjh
          · exact ⟨c', (hmem₂ c').2 (Or.inr hc'), hcj⟩
        · exact ⟨x, (hmem₂ x).2 (Or.inl ⟨hx, hxh⟩), hxj⟩
      · exact Or.inr ⟨b, hb, hbj⟩

/-! ### `maybeUnwrapNxt` -/

theorem Items.getElem!_type_of_lt {items : Items} {i : ItemId} (hi : i < items.size) :
    items[i]!.type = Items.type items i := by
  simp [Items.type, getElem!_pos items i hi, Array.getElem?_eq_getElem hi]

theorem Items.getElem!_ch_of_lt {items : Items} {i : ItemId} (hi : i < items.size) :
    items[i]!.ch = Items.ch items i := by
  simp [Items.ch, getElem!_pos items i hi, Array.getElem?_eq_getElem hi]

theorem Items.getElem!_type_of_le {items : Items} {i : ItemId} (hi : items.size ≤ i) :
    items[i]!.type = .F := by
  rw [getElem!_neg items i (Nat.not_lt.2 hi)]; rfl

theorem spans_eq_setSides_of_sides {sp : List ItemId × List ItemId} {dir : Bool} {l : List ItemId}
    (h1 : getSide sp dir = l) (h2 : getSide sp (!dir) = []) : sp = setSides dir l [] := by
  obtain ⟨a, b⟩ := sp; cases dir <;> simp_all [getSide, setSides]

/-- `maybeUnwrapNxt ty; mergeTstackTops; finishTstackTop` on `c :: t :: _` (both one-sided on
`stackDir[t.topDepth]`, which is also the direction at the merged depth), for `ty ∈ {S, P, R}`.
`hU` is `L1Unwrap s ty t` (`Spqr/EarInv.lean`) spelled out. -/
theorem StSim.unwrapMergeClose {g : Graph} (s : WalkState) (ty : NodeType) (c t : TEntry)
    (new base : List TEntry) (ps : List StPiece) (blocks : List StBlock)
    (hts : s.tstack = c :: t :: (new ++ base))
    (hc : getSide c.spans (!s.stackDir[min t.topDepth c.topDepth]!) = [])
    (ht : getSide t.spans (!s.stackDir[min t.topDepth c.topDepth]!) = [])
    (hty : ty = .S ∨ ty = .P ∨ ty = .R)
    (hU : ty ≠ .R → ∀ dir h, (getSide t.spans dir).head! = h → Items.type s.items h = ty →
      getSide t.spans dir = [h] ∧ (∀ p, ¬ Items.IsParent s.items p h) ∧
      ∀ c ∈ Items.ch s.items h, ∀ u ∈ s.tstack, c ∉ u.spans.1 ++ u.spans.2)
    (hR : StRead s.items (c :: t :: new) ps) (hI : StItems g s blocks) :
    let r := (maybeUnwrapNxt ty).run s
    let s' := ((finishTstackTop r.1).run (mergeTstackTops.run r.2).2).2
    s'.tstack = { TEntry.mergeInto c t with
        spans := setSides s.stackDir[min t.topDepth c.topDepth]! [r.1] [] } :: (new ++ base) ∧
      s'.stackDir = s.stackDir ∧ s.items.size ≤ s'.items.size ∧ r.1 < s'.items.size ∧
      Items.type s'.items r.1 = ty ∧
      StRead s'.items ({ TEntry.mergeInto c t with
        spans := setSides s.stackDir[min t.topDepth c.topDepth]! [r.1] [] } :: new) ps ∧
      StItems g s' blocks := by
  dsimp only
  have hty' : ty ≠ .V ∧ ty ≠ .Q := by rcases hty with rfl | rfl | rfl <;> simp
  have alloc := StSim.allocMergeClose s ty c t new base ps blocks hts hc ht hty' hR hI
  dsimp only at alloc
  have halloc : (allocItem ty).run s = (s.items.size, { s with items := s.items.push ⟨ty, (none, none), []⟩ }) := rfl
  have allocGoal : let r := (allocItem ty).run s
      let s' := ((finishTstackTop r.1).run (mergeTstackTops.run r.2).2).2
      s'.tstack = { TEntry.mergeInto c t with
          spans := setSides s.stackDir[min t.topDepth c.topDepth]! [r.1] [] } :: (new ++ base) ∧
        s'.stackDir = s.stackDir ∧ s.items.size ≤ s'.items.size ∧ r.1 < s'.items.size ∧
        Items.type s'.items r.1 = ty ∧
        StRead s'.items ({ TEntry.mergeInto c t with
          spans := setSides s.stackDir[min t.topDepth c.topDepth]! [r.1] [] } :: new) ps ∧
        StItems g s' blocks := by
    dsimp only
    rw [halloc]
    dsimp only
    exact ⟨alloc.1, alloc.2.1, by rw [alloc.2.2.1]; exact Nat.le_succ _, by rw [alloc.2.2.1]; exact Nat.lt_succ_self _,
      alloc.2.2.2.1, alloc.2.2.2.2.1, alloc.2.2.2.2.2⟩
  dsimp only at allocGoal
  rw [maybeUnwrapNxt_run_eq ty s c t (new ++ base) hts _ rfl _ rfl]
  split
  · exact allocGoal
  · split
    · rename_i hA hT
      have hne : ty ≠ .R := fun e => hA (Or.inl e)
      generalize hh : (getSide t.spans s.stackDir[t.topDepth]!).head! = h at hT ⊢
      have hlt : h < s.items.size := by
        by_contra hge
        rw [Items.getElem!_type_of_le (Nat.not_lt.1 hge)] at hT
        rcases hty with rfl | rfl | rfl <;> cases hT
      rw [Items.getElem!_type_of_lt hlt] at hT
      obtain ⟨hsingle, -, -⟩ := hU hne _ h hh hT
      by_cases hd : s.stackDir[t.topDepth]! = s.stackDir[min t.topDepth c.topDepth]!
      · rw [hd] at hsingle hh ⊢
        have hc' : getSide c.spans (!s.stackDir[t.topDepth]!) = [] := by rw [hd]; exact hc
        have ht' : getSide t.spans (!s.stackDir[t.topDepth]!) = [] := by rw [hd]; exact ht
        have htsp : t.spans = setSides s.stackDir[t.topDepth]! [h] [] := by
          rw [hd]; exact spans_eq_setSides_of_sides hsingle ht
        have htyh : ¬ (Items.type s.items h = .V ∨ Items.type s.items h = .Q) := by
          rw [hT]; exact fun e => e.elim hty'.1 hty'.2
        have reopen := StSim.reopenMergeClose s c t h new base ps blocks hts hc' htsp hd.symm htyh hR hI
        dsimp only at reopen
        rw [hd] at reopen
        rw [Items.getElem!_ch_of_lt hlt]
        dsimp only
        exact ⟨reopen.1, reopen.2.1, by rw [reopen.2.2.1], by rw [reopen.2.2.1]; exact hlt,
          by rw [reopen.2.2.2.1, hT], reopen.2.2.2.2.1, reopen.2.2.2.2.2⟩
      · exfalso
        have hd' : s.stackDir[t.topDepth]! = !s.stackDir[min t.topDepth c.topDepth]! := by
          revert hd
          cases s.stackDir[t.topDepth]! <;> cases s.stackDir[min t.topDepth c.topDepth]! <;> simp
        rw [hd', ht] at hsingle
        exact List.cons_ne_nil _ _ hsingle.symm
    · exact allocGoal

end Spqr
