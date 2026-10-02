import Spqr.StMerge
import Spqr.StUnwrap

/-! # The vertex close (`closeVert'`) and the st-reading

`StItemsX`: `StItems` while a fresh / reopened item `item` is pending (between `maybeUnwrapNxt` and
its `finishTstackTop`): the stack may hold `item`'s children, and `item` itself is not yet placed. -/

namespace Spqr
open WalkM WalkState

structure StItemsX (g : Graph) (s : WalkState) (blocks : List StBlock) (item : ItemId) : Prop where
  roots : ∀ x ∈ readStack s.tstack.tail, ∀ p, ¬ Items.IsParent s.items p x
  nodup : (readStack s.tstack).Nodup
  bounded : ∀ x ∈ readStack s.tstack, ∀ y, Items.Below s.items x y → y < s.items.size
  chLt : ∀ p c, Items.IsParent s.items p c → c < s.items.size
  chNodup : ∀ p, (Items.ch s.items p).Nodup
  closed : ∀ j, j < s.items.size → j ≠ item →
    Items.type s.items j = .S ∨ Items.type s.items j = .P ∨ Items.type s.items j = .R →
    (∃ x ∈ readStack s.tstack, Items.Below s.items x j) ∨ ∃ b ∈ blocks, InBlock g s.items b j
  lt : item < s.items.size
  root : ∀ p, ¬ Items.IsParent s.items p item
  notin : item ∉ readStack s.tstack
  notBelow : ∀ x ∈ readStack s.tstack, ¬ Items.Below s.items x item

theorem StItemsX.perm {g : Graph} {s s' : WalkState} {blocks : List StBlock} {item : ItemId}
    (h : StItemsX g s blocks item) (hread : (readStack s'.tstack).Perm (readStack s.tstack))
    (htail : ∀ x ∈ readStack s'.tstack.tail, x ∈ readStack s.tstack.tail)
    (hitems : s'.items = s.items) : StItemsX g s' blocks item := by
  obtain ⟨roots, nodup, bounded, chLt, chNodup, closed, lt, root, notin, notBelow⟩ := h
  refine ⟨fun x hx => by rw [hitems]; exact roots x (htail x hx), hread.nodup_iff.2 nodup,
    fun x hx => by rw [hitems]; exact bounded x (hread.mem_iff.1 hx), by rw [hitems]; exact chLt,
    by rw [hitems]; exact chNodup, fun j hj hji hty => ?_, by rw [hitems]; exact lt,
    by rw [hitems]; exact root, fun h => notin (hread.mem_iff.1 h),
    fun x hx => by rw [hitems]; exact notBelow x (hread.mem_iff.1 hx)⟩
  rw [hitems] at hj hty ⊢
  rcases closed j hj hji hty with ⟨x, hx, hxj⟩ | h
  · exact Or.inl ⟨x, hread.mem_iff.2 hx, hxj⟩
  · exact Or.inr h

theorem StItemsX.close {g : Graph} {s : WalkState} {blocks : List StBlock} {item : ItemId}
    (h : StItemsX g s blocks item) (t : TEntry) (rest : List TEntry) (hts : s.tstack = t :: rest)
    (hside : getSide t.spans (!s.stackDir[t.topDepth]!) = []) :
    StItems g ((WalkM.finishTstackTop item).run s).2 blocks :=
  StItems.close s blocks item t rest hts hside h.lt h.root h.notin
    (fun x hx => h.roots x (by rw [hts]; exact hx)) h.nodup h.bounded h.chLt h.chNodup h.closed

/-! ### The fold (`retarget`) -/

theorem readStack_fold_perm (dir : Bool) (curV : Nat) (t : TEntry) (rest : List TEntry) :
    (readStack ({ t with vStart := curV, spans := setSides dir (t.spans.1 ++ t.spans.2) [] } :: rest)).Perm
      (readStack (t :: rest)) := by
  rw [List.perm_iff_count]; intro a
  cases dir <;> simp [readStack, readL, readR, setSides, List.count_append] <;> omega

theorem StRead.fold {items : Items} (dir : Bool) (curV : Nat) (t : TEntry) (ps : List StPiece)
    (h : StRead items [t] ps) :
    StRead items [{ t with vStart := curV, spans := setSides dir (t.spans.1 ++ t.spans.2) [] }]
      [⟨dir, stNest ps⟩] := by
  have hflat : ExpandsList items (t.spans.1 ++ t.spans.2) (stNest ps) := by
    have := h.1.append h.2
    simpa [readL, readR, stNest] using this
  unfold StRead
  cases dir <;> simp [readL, readR, setSides, stNestL, stNestR]
  · exact ⟨hflat, .nil⟩
  · exact ⟨.nil, hflat⟩

/-! ### Loop 3 -/

theorem mergeTops_length (l : List TEntry) (h : 2 ≤ l.length) : (mergeTops l).length = l.length - 1 := by
  match l, h with
  | b :: a :: rest, _ => simp [mergeTops]

theorem loop3_length (o : Nat) : ∀ (fuel : Nat) (s : WalkState), o + 3 ≤ s.tstack.length →
    s.tstack.length ≤ fuel + o + 3 →
    ((loop fuel (loop3Cond o) mergeTstackTops).run s).2.tstack.length = o + 3 := by
  intro fuel
  induction fuel with
  | zero => intro s h1 h2; show s.tstack.length = o + 3; omega
  | succ fuel ih =>
    intro s h1 h2
    rw [loop_succ_run _ _ _ _ rfl, run_loop3Cond]
    dsimp only
    split
    · rename_i hc
      rw [decide_eq_true_eq] at hc
      rw [run_mergeTstackTops]
      refine ih _ ?_ ?_ <;> dsimp only <;> rw [mergeTops_length _ (by omega)] <;> omega
    · rename_i hc
      rw [decide_eq_true_eq] at hc
      show s.tstack.length = o + 3; omega

/-! ### `maybeUnwrapNxt; mergeTstackTops` -/

/-- `maybeUnwrapNxt ty; mergeTstackTops` on `c :: t :: _` (both one-sided on the direction at the
merged depth), for `ty ∈ {S, P, R}`: the merged top reads as `c :: t`, the returned item is pending. -/
theorem StSim.unwrapMerge {g : Graph} (s : WalkState) (ty : NodeType) (c t : TEntry)
    (new base : List TEntry) (ps : List StPiece) (blocks : List StBlock)
    (hts : s.tstack = c :: t :: (new ++ base))
    (hc : getSide c.spans (!s.stackDir[min t.topDepth c.topDepth]!) = [])
    (ht : getSide t.spans (!s.stackDir[min t.topDepth c.topDepth]!) = [])
    (hty : ty = .S ∨ ty = .P ∨ ty = .R)
    (hU : ty ≠ .R → ∀ h, (getSide t.spans s.stackDir[t.topDepth]!).head! = h →
      Items.type s.items h = ty → getSide t.spans s.stackDir[t.topDepth]! = [h])
    (hR : StRead s.items (c :: t :: new) ps) (hI : StItems g s blocks) :
    let r := (maybeUnwrapNxt ty).run s
    let s' := (mergeTstackTops.run r.2).2
    ∃ m : TEntry, s'.tstack = m :: (new ++ base) ∧ m.vStart = t.vStart ∧
      m.topDepth = min t.topDepth c.topDepth ∧
      getSide m.spans (!s.stackDir[min t.topDepth c.topDepth]!) = [] ∧
      s'.stackDir = s.stackDir ∧ s'.g = s.g ∧ s'.stackVerts = s.stackVerts ∧
      s.items.size ≤ s'.items.size ∧ Items.type s'.items r.1 = ty ∧
      StRead s'.items (m :: new) ps ∧ StItemsX g s' blocks r.1 := by
  dsimp only
  have hty' : ty ≠ .V ∧ ty ≠ .Q := by rcases hty with rfl | rfl | rfl <;> simp
  have hlt : ∀ x ∈ readStack s.tstack, x < s.items.size := fun x hx => hI.bounded x hx x .refl
  have hsub : ∀ x ∈ readStack (c :: t :: new), x ∈ readStack s.tstack := fun x hx => by
    rw [hts]; exact mem_readStack_append_left (l' := base) hx
  have hsubR : ∀ x ∈ readStack (new ++ base), x ∈ readStack s.tstack := fun x hx => by
    rw [hts]; exact mem_readStack_cons.2 (Or.inr (mem_readStack_cons.2 (Or.inr hx)))
  have htop : (TEntry.mergeInto c t).topDepth = min t.topDepth c.topDepth := rfl
  -- the allocation path
  have alloc : let r := (allocItem ty).run s
      let s' := (mergeTstackTops.run r.2).2
      ∃ m : TEntry, s'.tstack = m :: (new ++ base) ∧ m.vStart = t.vStart ∧
        m.topDepth = min t.topDepth c.topDepth ∧
        getSide m.spans (!s.stackDir[min t.topDepth c.topDepth]!) = [] ∧
        s'.stackDir = s.stackDir ∧ s'.g = s.g ∧ s'.stackVerts = s.stackVerts ∧
        s.items.size ≤ s'.items.size ∧ Items.type s'.items r.1 = ty ∧
        StRead s'.items (m :: new) ps ∧ StItemsX g s' blocks r.1 := by
    dsimp only
    have halloc : (allocItem ty).run s = (s.items.size, { s with items := s.items.push ⟨ty, (none, none), []⟩ }) := rfl
    rw [halloc]
    dsimp only
    have hs₂ := run_mergeTstackTops_cons_cons { s with items := s.items.push ⟨ty, (none, none), []⟩ } c t
      (new ++ base) hts
    generalize (mergeTstackTops.run { s with items := s.items.push ⟨ty, (none, none), []⟩ }).2 = s₂ at hs₂ ⊢
    have hts₂ : s₂.tstack = TEntry.mergeInto c t :: (new ++ base) := by rw [hs₂]
    have hitems₂ : s₂.items = s.items.push ⟨ty, (none, none), []⟩ := by rw [hs₂]
    have hread₂ : readStack s₂.tstack = readStack s.tstack := by
      rw [hts₂, readStack_mergeInto_cons, hts]
    have hsz : s₂.items.size = s.items.size + 1 := by rw [hitems₂, Array.size_push]
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
    refine ⟨TEntry.mergeInto c t, hts₂, rfl, rfl, TEntry.mergeInto_side_nil _ c t hc ht, by rw [hs₂],
      by rw [hs₂], by rw [hs₂], by rw [hsz]; exact Nat.le_succ _, by rw [hitems₂, Items.type_push_size], ?_, ?_⟩
    · unfold StRead at hR ⊢
      rw [readL_mergeInto_cons, readR_mergeInto_cons, hitems₂]
      exact ⟨hR.1.push _ (fun x hx y hy => hI.bounded x (hsub x (mem_readStack_of_readL hx)) y hy),
        hR.2.push _ (fun x hx y hy => hI.bounded x (hsub x (mem_readStack_of_readR hx)) y hy)⟩
    refine ⟨?_, by rw [hread₂]; exact hI.nodup, ?_, ?_, ?_, ?_, by rw [hsz]; exact Nat.lt_succ_self _, ?_,
      by rw [hread₂]; exact hnotin, ?_⟩
    · intro x hx p hpx
      rw [hts₂] at hx
      exact hI.roots x (hsubR x hx) p ((hparent₂ p x).1 hpx).2
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
    · intro p hp
      exact Nat.lt_irrefl _ (hI.chLt p _ ((hparent₂ p _).1 hp).2)
    · intro x hx h
      rw [hread₂] at hx
      exact Nat.lt_irrefl _ (hI.bounded x hx _ (hbelow₂ x _ (hlt x hx) h))
  dsimp only at alloc
  rw [maybeUnwrapNxt_run_eq ty s c t (new ++ base) hts _ rfl _ rfl]
  split
  · exact alloc
  split
  · rename_i hA hT
    have hne : ty ≠ .R := fun e => hA (Or.inl e)
    generalize hh : (getSide t.spans s.stackDir[t.topDepth]!).head! = h at hT ⊢
    have hlth : h < s.items.size := by
      by_contra hge
      rw [Items.getElem!_type_of_le (Nat.not_lt.1 hge)] at hT
      rcases hty with rfl | rfl | rfl <;> cases hT
    rw [Items.getElem!_type_of_lt hlth] at hT
    have hsingle := hU hne h hh hT
    have hd : s.stackDir[t.topDepth]! = s.stackDir[min t.topDepth c.topDepth]! := by
      by_contra hd
      have hd' : s.stackDir[t.topDepth]! = !s.stackDir[min t.topDepth c.topDepth]! := by
        revert hd
        cases s.stackDir[t.topDepth]! <;> cases s.stackDir[min t.topDepth c.topDepth]! <;> simp
      rw [hd', ht] at hsingle
      exact List.cons_ne_nil _ _ hsingle.symm
    have htsp : t.spans = setSides s.stackDir[t.topDepth]! [h] [] :=
      spans_eq_setSides_of_sides hsingle (by rw [hd]; exact ht)
    rw [Items.getElem!_ch_of_lt hlth]
    dsimp only
    generalize hdir : s.stackDir[t.topDepth]! = dir at hd htsp ⊢
    rw [← hd] at hc ht ⊢
    set t' : TEntry := { t with spans := setSides dir (Items.ch s.items h) [] } with ht'
    have hs₂ := run_mergeTstackTops_cons_cons { s with tstack := c :: t' :: (new ++ base) } c t' (new ++ base) rfl
    generalize (mergeTstackTops.run { s with tstack := c :: t' :: (new ++ base) }).2 = s₂ at hs₂ ⊢
    have hts₂ : s₂.tstack = TEntry.mergeInto c t' :: (new ++ base) := by rw [hs₂]
    have hitems₂ : s₂.items = s.items := by rw [hs₂]
    have hht : h ∈ t.spans.1 ++ t.spans.2 := by rw [htsp]; cases dir <;> simp [setSides]
    have hhmem : h ∈ readStack s.tstack := by
      rw [hts]; exact mem_readStack_cons.2 (Or.inr (mem_readStack_cons.2 (Or.inl hht)))
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
    have htyh : ¬ (Items.type s.items h = .V ∨ Items.type s.items h = .Q) := by
      rw [hT]; exact fun e => e.elim hty'.1 hty'.2
    refine ⟨TEntry.mergeInto c t', hts₂, rfl, rfl,
      TEntry.mergeInto_side_nil _ c t' hc (by rw [ht']; exact getSide_setSides_other _ _ _),
      by rw [hs₂], by rw [hs₂], by rw [hs₂], by rw [hitems₂], by rw [hitems₂, hT], ?_, ?_⟩
    · unfold StRead at hR ⊢
      rw [hitems₂, readL_mergeInto_cons, readR_mergeInto_cons,
        readL_reopen c t new dir h _ htsp ha hnew, readR_reopen c t new dir h _ htsp ha hnew]
      exact ⟨(ExpandsList.expandItem_self_iff htyh).2 hR.1, (ExpandsList.expandItem_self_iff htyh).2 hR.2⟩
    refine ⟨?_, by rw [hread₂]; exact nodup_expandItem hnd (hI.chNodup h) hchfresh, ?_,
      by rw [hitems₂]; exact hI.chLt, by rw [hitems₂]; exact hI.chNodup, ?_, by rw [hitems₂]; exact hlth,
      by rw [hitems₂]; exact hroot, hnotin₂, ?_⟩
    · intro x hx p hpx
      rw [hts₂] at hx
      rw [hitems₂] at hpx
      exact hI.roots x (hsubR x hx) p hpx
    · intro x hx y hy
      rw [hitems₂] at hy ⊢
      rcases (hmem₂ x).1 hx with ⟨hxs, _⟩ | hxc
      · exact hI.bounded x hxs y hy
      · exact hI.bounded h hhmem y (Relation.ReflTransGen.head (show Items.IsParent s.items h x from hxc) hy)
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
    · intro x hx
      rw [hitems₂]
      exact Items.not_below_root_of_ne hroot fun e => hnotin₂ (e ▸ hx)
  · exact alloc

end Spqr
