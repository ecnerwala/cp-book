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

/-! ### Splitting readings -/

theorem ExpandsList.unique {items : Items} {xs L L' : List ItemId} (h : ExpandsList items xs L)
    (h' : ExpandsList items xs L') : L = L' := by
  induction h generalizing L' with
  | nil => cases h'; rfl
  | leaf hx _ ih =>
    cases h' with
    | leaf _ h'' => rw [ih h'']
    | node hn => exact absurd hx hn
  | node hx _ ih =>
    cases h' with
    | leaf hl => exact absurd hl hx
    | node _ h'' => exact ih h''

theorem ExpandsList.split {items : Items} {a b L : List ItemId} (h : ExpandsList items (a ++ b) L) :
    ∃ A B, L = A ++ B ∧ ExpandsList items a A ∧ ExpandsList items b B := by
  generalize hab : a ++ b = ab at h
  induction h generalizing a with
  | nil =>
    obtain ⟨rfl, rfl⟩ := List.append_eq_nil_iff.1 hab
    exact ⟨[], [], rfl, .nil, .nil⟩
  | leaf hx hxs ih =>
    rename_i x xs L
    cases a with
    | nil => simp at hab; subst hab; exact ⟨[], x :: L, rfl, .nil, .leaf hx hxs⟩
    | cons y a' =>
      simp at hab; obtain ⟨rfl, hab⟩ := hab
      obtain ⟨A, B, rfl, hA, hB⟩ := ih hab
      exact ⟨y :: A, B, rfl, .leaf hx hA, hB⟩
  | node hx hxs ih =>
    rename_i x xs L
    cases a with
    | nil => simp at hab; subst hab; exact ⟨[], L, rfl, .nil, .node hx hxs⟩
    | cons y a' =>
      simp at hab; obtain ⟨rfl, hab⟩ := hab
      obtain ⟨A, B, rfl, hA, hB⟩ := ih (a := Items.ch items y ++ a') (by rw [List.append_assoc, hab])
      exact ⟨A, B, rfl, .node hx hA, hB⟩

/-- Concatenating segments: `hi` sits above `lo`, so its pieces come later. -/
theorem StRead.append {items : Items} {hi lo : List TEntry} {ps qs : List StPiece}
    (h₁ : StRead items hi qs) (h₂ : StRead items lo ps) : StRead items (hi ++ lo) (ps ++ qs) := by
  unfold StRead at *
  rw [readL_append, readR_append, stNestL_append, stNestR_append]
  exact ⟨h₁.1.append h₂.1, h₂.2.append h₁.2⟩

/-- A reading of `hi ++ lo` whose upper part is known splits (leaf expansions are unique). -/
theorem StRead.split {items : Items} {hi lo : List TEntry} {ps qs : List StPiece}
    (h : StRead items (hi ++ lo) (ps ++ qs)) (h₁ : StRead items hi qs) : StRead items lo ps := by
  unfold StRead at *
  rw [readL_append, readR_append, stNestL_append, stNestR_append] at h
  obtain ⟨A, B, hAB, hA, hB⟩ := h.1.split
  obtain ⟨A', B', hAB', hA', hB'⟩ := h.2.split
  rw [hA.unique h₁.1] at hAB
  rw [hB'.unique h₁.2] at hAB'
  exact ⟨List.append_cancel_left hAB ▸ hB, List.append_cancel_right hAB' ▸ hA'⟩

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
    (hd : s.stackDir[t.topDepth]! = s.stackDir[min t.topDepth c.topDepth]!)
    (hty : ty = .S ∨ ty = .P ∨ ty = .R)
    (hU : ty ≠ .R → ∀ h, (getSide t.spans s.stackDir[t.topDepth]!).head! = h →
      Items.type s.items h = ty → getSide t.spans s.stackDir[t.topDepth]! = [h])
    (hU' : ∀ h, (getSide t.spans s.stackDir[t.topDepth]!).head! = h →
      Items.type s.items h = ty → getSide t.spans (!s.stackDir[t.topDepth]!) = [])
    (hR : StRead s.items (c :: t :: new) ps) (hI : StItems g s blocks) :
    let r := (maybeUnwrapNxt ty).run s
    let s' := (mergeTstackTops.run r.2).2
    ∃ m : TEntry, s'.tstack = m :: (new ++ base) ∧ m.vStart = t.vStart ∧
      m.topDepth = min t.topDepth c.topDepth ∧
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
    refine ⟨TEntry.mergeInto c t, hts₂, rfl, rfl, by rw [hs₂],
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
    have htsp : t.spans = setSides s.stackDir[t.topDepth]! [h] [] :=
      spans_eq_setSides_of_sides hsingle (hU' h hh hT)
    rw [Items.getElem!_ch_of_lt hlth]
    dsimp only
    generalize hdir : s.stackDir[t.topDepth]! = dir at hd htsp ⊢
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
    refine ⟨TEntry.mergeInto c t', hts₂, rfl, rfl, by rw [hs₂], by rw [hs₂], by rw [hs₂], by rw [hitems₂], by rw [hitems₂, hT], ?_, ?_⟩
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

theorem StItems.perm {g : Graph} {s s' : WalkState} {blocks : List StBlock} (h : StItems g s blocks)
    (hread : (readStack s'.tstack).Perm (readStack s.tstack)) (hitems : s'.items = s.items) :
    StItems g s' blocks := by
  obtain ⟨roots, nodup, bounded, chLt, chNodup, closed⟩ := h
  refine ⟨fun x hx => by rw [hitems]; exact roots x (hread.mem_iff.1 hx), hread.nodup_iff.2 nodup,
    fun x hx => by rw [hitems]; exact bounded x (hread.mem_iff.1 hx), by rw [hitems]; exact chLt,
    by rw [hitems]; exact chNodup, fun j hj hty => ?_⟩
  rw [hitems] at hj hty ⊢
  rcases closed j hj hty with ⟨x, hx, hxj⟩ | h
  · exact Or.inl ⟨x, hread.mem_iff.2 hx, hxj⟩
  · exact Or.inr h

theorem loop3_state (o fuel : Nat) (st : WalkState) :
    ∃ l, ((loop fuel (loop3Cond o) mergeTstackTops).run st).2 = { st with tstack := l } := by
  obtain ⟨k, hk, -⟩ := loop_run_iter' (loop3Cond o) mergeTstackTops (fun _ => rfl) fuel st
  exact ⟨_, by rw [hk, iter_mergeTstackTops]⟩

/-! ### `closeVert'` -/

/-- The vertex close on `c :: mid ++ [py, vy] ++ base` reading as `ps` (the child's pieces and its
tree edge): afterwards one entry reads as the single piece `⟨!edgeDir, stNest ps⟩`. For type 1,
`mid = []` and `py` is the single-item entry at the lowpoint depth `lv`. -/
theorem closeVert_st {g : Graph} {s st : WalkState} {d lv : Nat} {ps : List StPiece}
    {blocks : List StBlock} {base : List TEntry} (curV origTstack : Nat) (edgeDir isType1 isSingle : Bool)
    (c : TEntry) (mid : List TEntry) (py vy : TEntry)
    {pre B : List TEntry} {qs : List StPiece} (hbase : base = pre ++ B)
    (hJ : L1StInv g s d ps blocks base st) (hJ' : L1StInv g s d (qs ++ ps) blocks B st)
    (hts : st.tstack = c :: mid ++ [py, vy] ++ base)
    (hlv : lv ≤ d) (horig : base.length = origTstack) (hsd : edgeDir = !s.stackDir[lv]!)
    (hcl : lv ≤ c.topDepth) (hvyt : lv ≤ vy.topDepth)
    (hmid : isType1 = true → mid = [])
    (hpy : isType1 = true → py.topDepth = lv ∧ ∃ i, py.spans = setSides s.stackDir[lv]! [i] []) :
    let s' := ((closeVert' curV edgeDir isType1 origTstack isSingle).run st).2
    ∃ m : TEntry, s'.tstack = m :: base ∧ m.vStart = curV ∧ s'.stackDir = st.stackDir ∧ s'.g = st.g ∧
      StRead s'.items (m :: pre) (qs ++ [⟨!edgeDir, stNest ps⟩]) ∧ StItems g s' blocks ∧
      (isType1 = true → m.topDepth = lv ∧ getSide m.spans (!st.stackDir[lv]!) = []) := by
  dsimp only
  have hrun : ((closeVert' curV edgeDir isType1 origTstack isSingle).run st).2 =
      ((vertFinish (result (vertUnwrap isType1 (cvB₁ isType1 origTstack isSingle st))
        (cvS₁ isType1 origTstack isSingle st)) (cvB₁ isType1 origTstack isSingle st)).run
        (cvS₅ curV edgeDir isType1 origTstack isSingle st)).2 := rfl
  rw [hrun]
  have hdirlv : st.stackDir[lv]! = s.stackDir[lv]! := hJ.dirs lv hlv
  cases isType1
  · -- type 2: loop 3, two merges, the fold
    have hS₁ : cvS₁ false origTstack isSingle st =
        ((loop st.tstack.length (loop3Cond origTstack) mergeTstackTops).run st).2 := by
      simp only [cvS₁, after, vertPre, Bool.not_false, ↓reduceIte, WalkM.run_bind, run_tstackSize]
      rfl
    obtain ⟨l₁, hl₁⟩ := loop3_state origTstack st.tstack.length st
    have hlen : (cvS₁ false origTstack isSingle st).tstack.length = origTstack + 3 := by
      rw [hS₁]
      refine loop3_length origTstack _ st ?_ (by omega)
      rw [hts]; simp; omega
    have hJ₁ : L1StInv g s d ps blocks base (cvS₁ false origTstack isSingle st) := by
      rw [hS₁]
      refine L1StInv.mergeLoop _ (fun _ => rfl) _ hJ ?_
      rw [← hS₁, hlen]; omega
    obtain ⟨new₁, hts₁, hR₁⟩ := hJ₁.read
    have hlen₁ : new₁.length = 3 := by
      rw [hts₁, List.length_append] at hlen; omega
    obtain ⟨a, b, c', rfl⟩ := List.length_eq_three.1 hlen₁
    have hJ₁' : L1StInv g s d (qs ++ ps) blocks B (cvS₁ false origTstack isSingle st) := by
      rw [hS₁]
      refine L1StInv.mergeLoop _ (fun _ => rfl) _ hJ' ?_
      rw [← hS₁, hlen, ← horig, hbase, List.length_append]; omega
    obtain ⟨new₁', hts₁', hR₁'⟩ := hJ₁'.read
    have hnew₁' : new₁' = [a, b, c'] ++ pre :=
      List.append_cancel_right (hts₁'.symm.trans (by rw [hts₁, hbase, List.append_assoc]))
    subst hnew₁'
    generalize hS : cvS₁ false origTstack isSingle st = S₁ at hts₁ hR₁ hJ₁ hS₁ hR₁'
    have hS₂ : cvS₂ false origTstack isSingle st = S₁ := by rw [← hS]; rfl
    have hS₃ : cvS₃ false origTstack isSingle st = { S₁ with tstack := TEntry.mergeInto a b :: c' :: base } := by
      show (mergeTstackTops.run (cvS₂ false origTstack isSingle st)).2 = _
      rw [hS₂]; exact run_mergeTstackTops_cons_cons _ a b (c' :: base) hts₁
    have hS₄ : cvS₄ false origTstack isSingle st =
        { S₁ with tstack := TEntry.mergeInto (TEntry.mergeInto a b) c' :: base } := by
      show (mergeTstackTops.run (cvS₃ false origTstack isSingle st)).2 = _
      rw [hS₃]; exact run_mergeTstackTops_cons_cons _ (TEntry.mergeInto a b) c' base rfl
    have hS₅ : cvS₅ curV edgeDir false origTstack isSingle st = ((retarget curV edgeDir).run (cvS₄ false origTstack isSingle st)).2 := rfl
    rw [retarget_run_eq curV edgeDir _ _ base (show (cvS₄ false origTstack isSingle st).tstack = _ by rw [hS₄]), hS₄] at hS₅
    have hfin : ((vertFinish (result (vertUnwrap false (cvB₁ false origTstack isSingle st)) S₁)
        (cvB₁ false origTstack isSingle st)).run
        (cvS₅ curV edgeDir false origTstack isSingle st)).2 = cvS₅ curV edgeDir false origTstack isSingle st := rfl
    rw [hfin, hS₅]
    dsimp only
    have hsd₁ : S₁.stackDir = st.stackDir := by rw [hS₁, hl₁]
    have hg₁ : S₁.g = st.g := by rw [hS₁, hl₁]
    refine ⟨_, rfl, rfl, hsd₁, hg₁, ?_, ?_, fun h => absurd h Bool.false_ne_true⟩
    · have hM : StRead S₁.items [TEntry.mergeInto (TEntry.mergeInto a b) c'] ps := by
        unfold StRead at hR₁ ⊢
        rw [readL_mergeInto_cons, readR_mergeInto_cons, readL_mergeInto_cons, readR_mergeInto_cons]
        exact hR₁
      have hM' : StRead S₁.items ([TEntry.mergeInto (TEntry.mergeInto a b) c'] ++ pre) (qs ++ ps) := by
        unfold StRead at hR₁' ⊢
        rw [List.singleton_append, readL_mergeInto_cons, readR_mergeInto_cons, readL_mergeInto_cons,
          readR_mergeInto_cons]
        exact hR₁'
      exact StRead.append (StRead.fold _ _ _ _ hM) (hM'.split hM)
    · refine StItems.perm hJ₁.items ?_ rfl
      dsimp only
      refine (readStack_fold_perm _ _ _ _).trans ?_
      rw [readStack_mergeInto_cons, readStack_mergeInto_cons, hts₁]
      exact List.Perm.refl _
  · -- type 1: the unwrap, two merges, the fold, the close
    have hmid' := hmid rfl; subst hmid'
    obtain ⟨hpyt, i, hpys⟩ := hpy rfl
    have hts' : st.tstack = c :: py :: ([vy] ++ base) := by rw [hts]; rfl
    have hmin : min py.topDepth c.topDepth = lv := by rw [hpyt]; exact Nat.min_eq_left hcl
    have hd : st.stackDir[py.topDepth]! = st.stackDir[min py.topDepth c.topDepth]! := by rw [hmin, hpyt]
    have hpys' : py.spans = setSides st.stackDir[py.topDepth]! [i] [] := by rw [hpyt, hdirlv]; exact hpys
    have htyS : (if isSingle then NodeType.S else NodeType.R) = .S ∨
        (if isSingle then NodeType.S else NodeType.R) = .P ∨
        (if isSingle then NodeType.S else NodeType.R) = .R := by cases isSingle <;> simp
    obtain ⟨new, hnew, hRn⟩ := hJ.read
    have hnew' : new = [c, py, vy] := List.append_cancel_right (hnew.symm.trans hts')
    subst hnew'
    obtain ⟨m, hts₃, hmv, hmt, hsd₃, hg₃, hsv₃, hsz₃, hty₃, hR₃, hX₃⟩ :=
      StSim.unwrapMerge st (if isSingle then NodeType.S else NodeType.R) c py [vy] base ps blocks hts' hd htyS
        (fun _ h hh _ => by rw [hpys', getSide_setSides] at hh ⊢; rw [← hh]; rfl)
        (fun h _ _ => by rw [hpys', getSide_setSides_other])
        hRn hJ.items
    obtain ⟨new', hnew'', hRn'⟩ := hJ'.read
    have hnew''' : new' = c :: py :: vy :: pre :=
      List.append_cancel_right (hnew''.symm.trans (by rw [hts', hbase]; rfl))
    subst hnew'''
    obtain ⟨m', hts₃', -, -, -, -, -, -, -, hR₃', -⟩ :=
      StSim.unwrapMerge st (if isSingle then NodeType.S else NodeType.R) c py (vy :: pre) B (qs ++ ps) blocks
        (by rw [hts', hbase]; rfl) hd htyS
        (fun _ h hh _ => by rw [hpys', getSide_setSides] at hh ⊢; rw [← hh]; rfl)
        (fun h _ _ => by rw [hpys', getSide_setSides_other])
        hRn' hJ.items
    have hmm : m' = m := (List.cons.inj (hts₃'.symm.trans hts₃)).1
    rw [hmm] at hR₃'
    have hS₃ : cvS₃ true origTstack isSingle st =
        (mergeTstackTops.run ((maybeUnwrapNxt (if isSingle then NodeType.S else NodeType.R)).run st).2).2 := rfl
    rw [← hS₃] at hts₃ hsd₃ hg₃ hsv₃ hsz₃ hty₃ hR₃ hX₃ hR₃'
    generalize hS : cvS₃ true origTstack isSingle st = S₃ at hts₃ hsd₃ hg₃ hsv₃ hsz₃ hty₃ hR₃ hX₃ hR₃'
    have hitem : result (vertUnwrap true (cvB₁ true origTstack isSingle st)) (cvS₁ true origTstack isSingle st) =
        some ((maybeUnwrapNxt (if isSingle then NodeType.S else NodeType.R)).run st).1 := rfl
    rw [hitem]
    generalize hI : ((maybeUnwrapNxt (if isSingle then NodeType.S else NodeType.R)).run st).1 = item at hty₃ hX₃
    have hS₄ : cvS₄ true origTstack isSingle st = { S₃ with tstack := TEntry.mergeInto m vy :: base } := by
      show (mergeTstackTops.run (cvS₃ true origTstack isSingle st)).2 = _
      rw [hS]; exact run_mergeTstackTops_cons_cons _ m vy base hts₃
    have hS₅ : cvS₅ curV edgeDir true origTstack isSingle st = ((retarget curV edgeDir).run (cvS₄ true origTstack isSingle st)).2 := rfl
    rw [retarget_run_eq curV edgeDir _ _ base (show (cvS₄ true origTstack isSingle st).tstack = _ by rw [hS₄]), hS₄] at hS₅
    set M : TEntry := TEntry.mergeInto m vy with hM
    generalize hF : ({ M with vStart := curV, spans := setSides (!edgeDir) (M.spans.1 ++ M.spans.2) [] } : TEntry) = F at hS₅
    have hfin : ((vertFinish (some item) (cvB₁ true origTstack isSingle st)).run
        (cvS₅ curV edgeDir true origTstack isSingle st)).2 =
        ((finishTstackTop item).run (cvS₅ curV edgeDir true origTstack isSingle st)).2 := rfl
    rw [hfin, hS₅]
    generalize hS₅' : ({ S₃ with tstack := F :: base } : WalkState) = S₅
    have hts₅ : S₅.tstack = F :: ([] ++ base) := by rw [← hS₅']; rfl
    have hitems₅ : S₅.items = S₃.items := by rw [← hS₅']
    have hsd₅ : S₅.stackDir = S₃.stackDir := by rw [← hS₅']
    have hFt : F.topDepth = lv := by
      rw [← hF]; show M.topDepth = lv
      rw [hM]; show min vy.topDepth m.topDepth = lv
      rw [hmt, hmin]; exact Nat.min_eq_right hvyt
    have hside : getSide F.spans (!S₅.stackDir[F.topDepth]!) = [] := by
      rw [hsd₅, hsd₃, hFt, hdirlv, ← hF]
      show getSide (setSides (!edgeDir) _ []) (!s.stackDir[lv]!) = []
      rw [hsd, Bool.not_not]; exact getSide_setSides_other _ _ _
    have hX₅ : StItemsX g S₅ blocks item := by
      refine hX₃.perm ?_ ?_ hitems₅
      · rw [hts₅, hts₃, ← hF]
        refine (readStack_fold_perm _ _ _ _).trans ?_
        rw [hM, readStack_mergeInto_cons]; rfl
      · intro x hx
        rw [hts₅] at hx; rw [hts₃]
        exact mem_readStack_cons.2 (Or.inr hx)
    have hMr : StRead S₃.items [M] ps := by
      unfold StRead at hR₃ ⊢
      rw [hM, readL_mergeInto_cons, readR_mergeInto_cons]
      exact hR₃
    have hMr' : StRead S₃.items ([M] ++ pre) (qs ++ ps) := by
      unfold StRead at hR₃' ⊢
      rw [List.singleton_append, hM, readL_mergeInto_cons, readR_mergeInto_cons]
      exact hR₃'
    have hR₅ : StRead S₅.items (F :: pre) (qs ++ [⟨!edgeDir, stNest ps⟩]) := by
      rw [hitems₅, ← hF]
      exact StRead.append (StRead.fold _ _ _ _ hMr) (hMr'.split hMr)
    have htyi : ¬ (Items.type S₅.items item = .V ∨ Items.type S₅.items item = .Q) := by
      rw [hitems₅, hty₃]; cases isSingle <;> simp
    have hts₅' : S₅.tstack = F :: (pre ++ B) := by rw [hts₅, hbase]; rfl
    obtain ⟨hts₆, hR₆⟩ := StRead.finishTstackTop S₅ item F pre B (qs ++ [⟨!edgeDir, stNest ps⟩]) hts₅' hside
      hX₅.lt htyi
      (fun x hx => hX₅.notBelow x (by rw [hts₅']; exact mem_readStack_append_left hx))
      (fun hx => hX₅.notin (by rw [hts₅']; exact mem_readStack_cons.2 (Or.inr (mem_readStack_append_left hx))))
      hR₅
    obtain ⟨x', items', hrun₆, -⟩ := finishTstackTop_run item S₅ hts₅
    refine ⟨{ F with spans := setSides S₅.stackDir[F.topDepth]! [item] [] }, by rw [hts₆, hbase],
      ?_, ?_, ?_, ?_, hX₅.close F base hts₅ hside, fun _ => ⟨hFt, ?_⟩⟩
    · show F.vStart = curV; rw [← hF]
    · rw [hrun₆]; show S₅.stackDir = st.stackDir; rw [hsd₅, hsd₃]
    · rw [hrun₆]; show S₅.g = st.g; rw [← hS₅', hg₃]
    · exact hR₆
    · show getSide (setSides S₅.stackDir[F.topDepth]! [item] []) _ = []
      rw [hsd₅, hsd₃, hFt]; exact getSide_setSides_other _ _ _

/-! ### `finishP` and `finishTail` -/

theorem mergeTstackTops_g (s : WalkState) : (mergeTstackTops.run s).2.g = s.g := by
  rfl

theorem finishTstackTop_g (item : ItemId) (s : WalkState) :
    ((finishTstackTop item).run s).2.g = s.g := by
  rfl

theorem maybeUnwrapNxt_g (ty : NodeType) (s : WalkState) :
    ((maybeUnwrapNxt ty).run s).2.g = s.g := by
  obtain ⟨g, tern, items, sv, sd, nei, fo, tstack, tb, tsl⟩ := s
  simp only [maybeUnwrapNxt, WalkM.run_bind, WalkM.get_run, modifyNxt, Bool.or_eq_true, beq_iff_eq]
  by_cases h : ty = NodeType.R ∨ tern = true
  · simp only [h, ↓reduceIte]; rfl
  · simp only [h, ↓reduceIte, WalkM.run_bind, run_nxt, run_stackDir, run_getItem]
    by_cases h₂ : items[(getSide tstack.tail.head!.spans sd[tstack.tail.head!.topDepth]!).head!]!.type = ty <;>
      simp only [h₂, ↓reduceIte] <;> rfl

/-- The type-1 P-check on `c :: new ++ base` (`c` one-sided on `stackDir[lowval]`, no entry of
`base` starts at `curV`): the P-target `t` (`L1Unwrap s .P t`, one-sided) is closed with `c` into
one P item; otherwise nothing happens. -/
theorem finishP_st {g : Graph} (s : WalkState) (curV lowval : Nat) (isType1 : Bool) (c : TEntry)
    (new base : List TEntry) (ps : List StPiece) (blocks : List StBlock)
    (hts : s.tstack = c :: (new ++ base))
    (hB : ∀ t ∈ base, t.vStart ≠ curV)
    (hP : isType1 = true → ∀ t ∈ new, t.vStart = curV → t.topDepth = lowval →
      getSide c.spans (!s.stackDir[lowval]!) = [] ∧ lowval ≤ c.topDepth ∧
      getSide t.spans (!s.stackDir[lowval]!) = [] ∧ L1Unwrap s .P t)
    (hR : StRead s.items (c :: new) ps) (hI : StItems g s blocks) :
    let s' := ((finishP curV lowval isType1).run s).2
    ∃ new', new' ≠ [] ∧ s'.tstack = new' ++ base ∧ s'.stackDir = s.stackDir ∧ s'.g = s.g ∧
      s.items.size ≤ s'.items.size ∧ StRead s'.items new' ps ∧ StItems g s' blocks := by
  dsimp only
  simp only [Spqr.finishP, WalkM.run_bind, run_condP]
  by_cases h : result (condP curV lowval isType1) s = true
  · have h' : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == curV) &&
        (s.tstack.tail.head!.topDepth == lowval)) = true := h
    simp only [h', ↓reduceIte, WalkM.run_bind]
    simp only [Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq] at h'
    obtain ⟨⟨⟨h1, hlen⟩, hvs⟩, htd⟩ := h'
    rw [hts] at hlen hvs htd
    simp only [List.tail_cons] at hvs htd
    obtain ⟨t, new', rfl⟩ : ∃ t new', new = t :: new' := by
      cases new with
      | nil =>
        exfalso
        cases base with
        | nil => simp at hlen
        | cons b bs => exact hB b List.mem_cons_self (by simpa using hvs)
      | cons t new' => exact ⟨t, new', rfl⟩
    simp only [List.cons_append, List.head!_cons] at hvs htd
    obtain ⟨hc, hct, ht, hU⟩ := hP h1 t List.mem_cons_self hvs htd
    have hmin : min t.topDepth c.topDepth = lowval := by rw [htd]; exact Nat.min_eq_left hct
    have res := StSim.unwrapMergeClose s .P c t new' base ps blocks hts (by rw [hmin]; exact hc)
      (by rw [hmin]; exact ht) (Or.inr (Or.inl rfl)) (fun _ => hU) hR hI
    dsimp only at res
    obtain ⟨hts', hsd', hsz', -, -, hR', hI'⟩ := res
    refine ⟨_, List.cons_ne_nil _ _, hts', hsd', ?_, hsz', hR', hI'⟩
    rw [finishTstackTop_g, mergeTstackTops_g, maybeUnwrapNxt_g]
  · have h' : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == curV) &&
        (s.tstack.tail.head!.topDepth == lowval)) = false := Bool.eq_false_iff.2 h
    simp only [h', Bool.false_eq_true, ↓reduceIte, WalkM.pure_run]
    refine ⟨c :: new, List.cons_ne_nil _ _, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> trivial

/-- The first-edge vertex push on `new ++ base`: the piece `⟨stackDir[d], [V curV]⟩` is appended
(merged into the top entry unless `isSingle`). -/
theorem finishTail_st {g : Graph} (s : WalkState) (curV d : Nat) (hasVert isSingle : Bool)
    (new base : List TEntry) (ps : List StPiece) (blocks : List StBlock)
    (hts : s.tstack = new ++ base)
    (hne : hasVert = false → isSingle = false → new ≠ [])
    (hvert : hasVert = false → vertItem curV < s.items.size ∧ Items.type s.items (vertItem curV) = .V ∧
      (∀ p, ¬ Items.IsParent s.items p (vertItem curV)) ∧ vertItem curV ∉ readStack s.tstack)
    (hR : StRead s.items new ps) (hI : StItems g s blocks) :
    let r := (finishTail curV d hasVert isSingle).run s
    r.1 = true ∧ r.2.stackDir = s.stackDir ∧ r.2.g = s.g ∧ r.2.items = s.items ∧
    ∃ new', r.2.tstack = new' ++ base ∧
      StRead r.2.items new' (if hasVert then ps else ps ++ [⟨s.stackDir[d]!, [vertItem curV]⟩]) ∧
      StItems g r.2 blocks := by
  dsimp only
  cases hasVert
  · obtain ⟨hv, hvty, hvroot, hvfree⟩ := hvert rfl
    set E : TEntry := ⟨curV, d, s.nxtEdgeIdx, setSides s.stackDir[d]! [vertItem curV] []⟩ with hE
    have hR₁ : StRead s.items (E :: new) (ps ++ [⟨s.stackDir[d]!, [vertItem curV]⟩]) :=
      StRead.pushEntry curV d s.nxtEdgeIdx s.stackDir[d]! (vertItem curV) (Or.inl hvty) hR
    have hI₁ : StItems g { s with tstack := E :: s.tstack } blocks :=
      StItems.pushEntry curV d s.nxtEdgeIdx s.stackDir[d]! (vertItem curV) hv hvroot hvfree hI
    cases isSingle
    · have hr : (finishTail curV d false false).run s =
          (true, (mergeTstackTops.run { s with tstack := E :: s.tstack }).2) := rfl
      rw [hr]
      obtain ⟨c, new', rfl⟩ : ∃ c new', new = c :: new' := by
        cases new with
        | nil => exact absurd rfl (hne rfl rfl)
        | cons c new' => exact ⟨c, new', rfl⟩
      have hs₂ := run_mergeTstackTops_cons_cons { s with tstack := E :: s.tstack } E c (new' ++ base)
        (by show E :: s.tstack = _; rw [hts]; rfl)
      rw [hs₂]
      refine ⟨rfl, rfl, rfl, rfl, TEntry.mergeInto E c :: new', rfl, ?_, ?_⟩
      · unfold StRead at hR₁ ⊢
        rw [readL_mergeInto_cons, readR_mergeInto_cons]
        exact hR₁
      · refine StItems.perm hI₁ ?_ rfl
        dsimp only
        rw [readStack_mergeInto_cons, hts]
        exact List.Perm.refl _
    · have hr : (finishTail curV d false true).run s = (true, { s with tstack := E :: s.tstack }) := rfl
      rw [hr]
      exact ⟨rfl, rfl, rfl, rfl, E :: new, by rw [hts]; rfl, hR₁, hI₁⟩
  · exact ⟨rfl, rfl, rfl, rfl, new, hts, hR, hI⟩

end Spqr
