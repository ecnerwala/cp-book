import Spqr.StSim
import Spqr.StWalk

/-!
# `StSim` through `finishTstackTop`

The stack reading (`StRead`) and the item facts (`StItems`) across a close of the one-sided top
entry into `item`.
-/

namespace Spqr
open WalkM

/-! ### `Below` framing -/

theorem Items.Below.modify_of_not_below {items : Items} {i x y : ItemId} (f : Item → Item)
    (hx : ¬ Items.Below items x i) (h : Items.Below items x y) : Items.Below (items.modify i f) x y := by
  induction h with
  | refl => exact .refl
  | @tail b c hxb hbc ih =>
    refine ih.tail ?_
    have hbi : b ≠ i := fun e => hx (e ▸ hxb)
    show c ∈ Items.ch (items.modify i f) b
    rw [Items.ch_modify_ne _ _ _ _ (Ne.symm hbi)]
    exact hbc

theorem Items.Below.of_modify {items : Items} {i x y : ItemId} (f : Item → Item)
    (hx : ¬ Items.Below items x i) (h : Items.Below (items.modify i f) x y) : Items.Below items x y := by
  induction h with
  | refl => exact .refl
  | @tail b c hxb hbc ih =>
    refine ih.tail ?_
    have hbi : b ≠ i := fun e => hx (e ▸ ih)
    have hbc : c ∈ Items.ch (items.modify i f) b := hbc
    rwa [Items.ch_modify_ne _ _ _ _ (Ne.symm hbi)] at hbc

theorem Items.Below.push {items : Items} (it : Item) {x y : ItemId}
    (hb : ∀ z, Items.Below items x z → z < items.size) (h : Items.Below items x y) : Items.Below (items.push it) x y := by
  induction h with
  | refl => exact .refl
  | @tail b c hxb hbc ih =>
    refine ih.tail ?_
    show c ∈ Items.ch (items.push it) b
    rw [Items.ch_push, ite_eq_right_iff.2 (fun e => absurd e (Nat.ne_of_lt (hb b hxb)))]
    exact hbc

theorem Items.Below.of_push {items : Items} (it : Item) {x y : ItemId}
    (hlt : ∀ p c, items.IsParent p c → c < items.size) (hx : x < items.size)
    (h : Items.Below (items.push it) x y) : Items.Below items x y := by
  induction h with
  | refl => exact .refl
  | @tail b c hxb hbc ih =>
    refine ih.tail ?_
    have hb : b < items.size := by
      rcases Relation.ReflTransGen.cases_tail ih with e | ⟨p, _, hp⟩
      · exact e ▸ hx
      · exact hlt p b hp
    have hbc : c ∈ Items.ch (items.push it) b := hbc
    rwa [Items.ch_push, ite_eq_right_iff.2 (fun e => absurd e (Nat.ne_of_lt hb))] at hbc

/-! ### The close -/

theorem run_finishTstackTop_tstack (s : WalkState) (item : ItemId) (t : TEntry) (rest : List TEntry)
    (hts : s.tstack = t :: rest) :
    ((finishTstackTop item).run s).2.tstack =
      { t with spans := setSides s.stackDir[t.topDepth]! [item] [] } :: rest := by
  rcases s with ⟨g, tern, items, sv, sd, nei, fo, ts, tb, tsl⟩
  simp only at hts; subst hts
  rfl

/-- Closing the one-sided top `t` into `item` (not on the stack below): expanding `item` back
recovers the reading. -/
theorem readStack_close (dir : Bool) (item : ItemId) (t : TEntry) (rest : List TEntry)
    (hside : getSide t.spans (!dir) = []) (hnew : item ∉ readStack rest) :
    expandItem item (getSide t.spans dir)
        (readStack ({ t with spans := setSides dir [item] [] } :: rest)) =
      readStack (t :: rest) := by
  rw [readStack_cons_setSides, expandItem_append, expandItem_append, expandItem_of_not_mem _ _ _ hnew]
  cases dir <;> simp [getSide] at hside ⊢ <;>
    simp [readStack, readL, readR, expandItem, stNestL, stNestR, hside]

/-- `finishTstackTop item` on the one-sided top `t` of the segment above `base` keeps the reading. -/
theorem StRead.finishTstackTop (s : WalkState) (item : ItemId) (t : TEntry) (new base : List TEntry)
    (ps : List StPiece) (hts : s.tstack = t :: (new ++ base))
    (hside : getSide t.spans (!s.stackDir[t.topDepth]!) = [])
    (hi : item < s.items.size)
    (hty : ¬ (Items.type s.items item = .V ∨ Items.type s.items item = .Q))
    (hnb : ∀ x ∈ readStack (t :: new), ¬ Items.Below s.items x item)
    (hnew : item ∉ readStack new)
    (h : StRead s.items (t :: new) ps) :
    ((finishTstackTop item).run s).2.tstack =
        { t with spans := setSides s.stackDir[t.topDepth]! [item] [] } :: (new ++ base) ∧
      StRead ((finishTstackTop item).run s).2.items
        ({ t with spans := setSides s.stackDir[t.topDepth]! [item] [] } :: new) ps := by
  refine ⟨run_finishTstackTop_tstack s item t (new ++ base) hts, ?_⟩
  rw [finishTstackTop_items s item t (new ++ base) hts]
  exact ExpandsList.close
    (f := fun it => { it with
      vs := setSides s.stackDir[t.topDepth]! (some (t.top s)) (some t.vStart),
      ch := getSide t.spans s.stackDir[t.topDepth]! })
    (fun _ => rfl) hi hty (hl := readStack_close _ item t new hside hnew) rfl hnb h


/-! ### `StItems` through the close -/

theorem Items.Below.lt_of_chLt {items : Items} (hlt : ∀ p c, Items.IsParent items p c → c < items.size)
    {x y : ItemId} (hx : x < items.size) (h : Items.Below items x y) : y < items.size := by
  rcases Relation.ReflTransGen.cases_tail h with e | ⟨p, _, hp⟩
  · exact e ▸ hx
  · exact hlt p y hp

/-- `VsOrientedAt` only looks at `i`, its children and the descendants of `i`. -/
theorem VsOrientedAt.congr {g : Graph} {items items' : Items} {b : StBlock} {i : ItemId} {L : List ItemId}
    (hvs : Items.vs items' i = Items.vs items i) (hch : Items.ch items' i = Items.ch items i)
    (hc : ∀ c ∈ Items.ch items i, Items.type items' c = Items.type items c ∧ Items.vs items' c = Items.vs items c)
    (hbi : ∀ e, Items.Below items' i e → Items.Below items i e)
    (hbc : ∀ c ∈ Items.ch items i, ∀ e, Items.Below items' c e → Items.Below items i e → Items.Below items c e)
    (h : VsOrientedAt g items b i L) : VsOrientedAt g items' b i L := by
  obtain ⟨h1, h2, h3, h4, h5, h6⟩ := h
  refine ⟨h1, fun e he hb hbe => h2 e he hb (hbi _ hbe), hvs ▸ h3, ?_, ?_, ?_⟩
  · intro s t hst c hc' hV
    exact h4 s t (hvs ▸ hst) c (hch ▸ hc') ((hc c (hch ▸ hc')).1 ▸ hV)
  · intro c hc' hV
    rw [hch] at hc'
    rw [(hc c hc').2]
    exact h5 c hc' ((hc c hc').1 ▸ hV)
  · intro c hc' hV u v huv e he hb hbe y hy
    rw [hch] at hc'
    refine h6 c hc' ((hc c hc').1 ▸ hV) u v ((hc c hc').2 ▸ huv) e he hb ?_ y hy
    have hp : Items.IsParent items' i c := by unfold Items.IsParent; rw [hch]; exact hc'
    exact hbc c hc' (edgeItem g e) hbe (hbi _ (Relation.ReflTransGen.head hp hbe))

theorem InBlock.congr {g : Graph} {items items' : Items} {b : StBlock} {i : ItemId}
    (hE : ∀ y, Items.Below items i y →
      Items.type items' y = Items.type items y ∧ Items.ch items' y = Items.ch items y)
    (hvs : Items.vs items' i = Items.vs items i)
    (hc : ∀ c ∈ Items.ch items i, Items.vs items' c = Items.vs items c)
    (hbi : ∀ e, Items.Below items' i e → Items.Below items i e)
    (hbc : ∀ c ∈ Items.ch items i, ∀ e, Items.Below items' c e → Items.Below items i e → Items.Below items c e)
    (h : InBlock g items b i) : InBlock g items' b i := by
  obtain ⟨L, hL, hseg, hor⟩ := h
  refine ⟨L, hL.congr (fun x hx => by simp at hx; subst hx; exact hE), hseg, hor.congr hvs (hE i .refl).2 ?_ hbi hbc⟩
  intro c hc'
  exact ⟨(hE c (.single hc')).1, hc c hc'⟩

/-- Modifying a root `j ≠ i` keeps `i` finished. -/
theorem InBlock.modify_root {g : Graph} {items : Items} {b : StBlock} {i j : ItemId} (f : Item → Item)
    (hr : ∀ p, ¬ Items.IsParent items p j) (hij : i ≠ j) (h : InBlock g items b i) :
    InBlock g (items.modify j f) b i := by
  have hnb : ∀ x, Items.Below items i x → x ≠ j := fun x hx e =>
    Items.not_below_root_of_ne hr hij (e ▸ hx)
  refine h.congr (fun y hy => ?_) ?_ (fun c hc => ?_) (fun e he => he.of_modify f (Items.not_below_root_of_ne hr hij))
    (fun c hc e he _ => he.of_modify f (Items.not_below_root_of_ne hr (hnb c (.single hc))))
  · have := hnb y hy
    exact ⟨by simp [Items.type, Array.getElem?_modify, this.symm], Items.ch_modify_ne _ _ _ _ this.symm⟩
  · exact Items.vs_modify_ne _ _ _ _ hij.symm
  · exact Items.vs_modify_ne _ _ _ _ (hnb c (.single hc)).symm

/-- Pushing a new item keeps `i` finished. -/
theorem InBlock.push {g : Graph} {items : Items} {b : StBlock} {i : ItemId} (it : Item)
    (hlt : ∀ p c, Items.IsParent items p c → c < items.size) (hi : i < items.size)
    (h : InBlock g items b i) : InBlock g (items.push it) b i := by
  have hlt' : ∀ y, Items.Below items i y → y ≠ items.size := fun y hy =>
    Nat.ne_of_lt (Items.Below.lt_of_chLt hlt hi hy)
  refine h.congr (fun y hy => ?_) ?_ (fun c hc => ?_) (fun e he => he.of_push it hlt hi)
    (fun c hc e he _ => he.of_push it hlt (hlt i c hc))
  · rw [Items.type_push, Items.ch_push, ite_eq_right_iff.2 (fun e => absurd e (hlt' y hy)), ite_eq_right_iff.2 (fun e => absurd e (hlt' y hy))]; exact ⟨rfl, rfl⟩
  · rw [Items.vs_push, ite_eq_right_iff.2 (fun e => absurd e (hlt' i .refl))]
  · rw [Items.vs_push, ite_eq_right_iff.2 (fun e => absurd e (hlt' c (.single hc)))]

/-- Reading the stack with a one-sided entry on top, as a permutation. -/
theorem readStack_cons_perm (dir : Bool) (t : TEntry) (rest : List TEntry)
    (hside : getSide t.spans (!dir) = []) :
    (readStack (t :: rest)).Perm (getSide t.spans dir ++ readStack rest) := by
  rw [List.perm_iff_count]; intro a
  cases dir <;> simp [getSide] at hside ⊢ <;> simp [readStack, readL, readR, hside, List.count_append]; omega

theorem readStack_close_perm (dir : Bool) (item : ItemId) (t : TEntry) (rest : List TEntry) :
    (readStack ({ t with spans := setSides dir [item] [] } :: rest)).Perm (item :: readStack rest) := by
  rw [readStack_cons_setSides, List.perm_iff_count]; intro a
  cases dir <;> simp [stNestL, stNestR, List.count_append, List.count_cons]

/-- The item facts across `finishTstackTop item` on the one-sided top `t`: `item` is a root not on
the stack, the rest of the stack holds roots, and every other S / P / R item is live or finished. -/
theorem StItems.close {g : Graph} (s : WalkState) (blocks : List StBlock) (item : ItemId)
    (t : TEntry) (rest : List TEntry) (hts : s.tstack = t :: rest)
    (hside : getSide t.spans (!s.stackDir[t.topDepth]!) = [])
    (hi : item < s.items.size)
    (hroot : ∀ p, ¬ Items.IsParent s.items p item)
    (hitem : item ∉ readStack s.tstack)
    (roots : ∀ x ∈ readStack rest, ∀ p, ¬ Items.IsParent s.items p x)
    (nodup : (readStack s.tstack).Nodup)
    (bounded : ∀ x ∈ readStack s.tstack, ∀ y, Items.Below s.items x y → y < s.items.size)
    (chLt : ∀ p c, Items.IsParent s.items p c → c < s.items.size)
    (chNodup : ∀ p, (Items.ch s.items p).Nodup)
    (closed : ∀ j, j < s.items.size → j ≠ item →
      Items.type s.items j = .S ∨ Items.type s.items j = .P ∨ Items.type s.items j = .R →
      (∃ x ∈ readStack s.tstack, Items.Below s.items x j) ∨ ∃ b ∈ blocks, InBlock g s.items b j) :
    StItems g ((WalkM.finishTstackTop item).run s).2 blocks := by
  have hts' := run_finishTstackTop_tstack s item t rest hts
  have hitems := finishTstackTop_items s item t rest hts
  have hperm := readStack_cons_perm s.stackDir[t.topDepth]! t rest hside
  have hperm' := readStack_close_perm s.stackDir[t.topDepth]! item t rest
  have hmem : ∀ x, x ∈ readStack s.tstack ↔
      x ∈ getSide t.spans s.stackDir[t.topDepth]! ∨ x ∈ readStack rest := by
    intro x; rw [hts, hperm.mem_iff, List.mem_append]
  have hmem' : ∀ x, x ∈ readStack ((WalkM.finishTstackTop item).run s).2.tstack ↔
      x = item ∨ x ∈ readStack rest := by
    intro x; rw [hts', hperm'.mem_iff, List.mem_cons]
  have hnd : (getSide t.spans s.stackDir[t.topDepth]! ++ readStack rest).Nodup :=
    hperm.nodup_iff.1 (hts ▸ nodup)
  have hCnd := hnd.of_append_left
  have hrestnd := hnd.of_append_right
  have hdisj : ∀ x, x ∈ getSide t.spans s.stackDir[t.topDepth]! → x ∉ readStack rest :=
    fun x hx hx' => List.disjoint_of_nodup_append hnd hx hx'
  have hCsize : ∀ x ∈ getSide t.spans s.stackDir[t.topDepth]!, x < s.items.size :=
    fun x hx => bounded x ((hmem x).2 (Or.inl hx)) x .refl
  have hitemC : item ∉ getSide t.spans s.stackDir[t.topDepth]! := fun h => hitem ((hmem item).2 (Or.inl h))
  have hitemR : item ∉ readStack rest := fun h => hitem ((hmem item).2 (Or.inr h))
  have hnb : ∀ x ∈ readStack s.tstack, ¬ Items.Below s.items x item := fun x hx =>
    Items.not_below_root_of_ne hroot (fun e => hitem (e ▸ hx))
  have hch_item : Items.ch ((WalkM.finishTstackTop item).run s).2.items item =
      getSide t.spans s.stackDir[t.topDepth]! := by
    rw [hitems, Items.ch_modify_self _ _ _ hi]
  have hch_ne : ∀ p, p ≠ item →
      Items.ch ((WalkM.finishTstackTop item).run s).2.items p = Items.ch s.items p := fun p hp => by
    rw [hitems, Items.ch_modify_ne _ _ _ _ (Ne.symm hp)]
  have hsize : ((WalkM.finishTstackTop item).run s).2.items.size = s.items.size := by
    rw [hitems, Array.size_modify]
  have htype : ∀ j, Items.type ((WalkM.finishTstackTop item).run s).2.items j = Items.type s.items j :=
    fun j => by rw [hitems]; exact Items.type_modify _ _ _ _ (fun _ => rfl)
  have hparent : ∀ p c, Items.IsParent ((WalkM.finishTstackTop item).run s).2.items p c ↔
      (p = item ∧ c ∈ getSide t.spans s.stackDir[t.topDepth]!) ∨
        (p ≠ item ∧ Items.IsParent s.items p c) := by
    intro p c
    by_cases hp : p = item
    · subst hp; simp [Items.IsParent, hch_item]
    · simp [Items.IsParent, hch_ne p hp, hp]
  have hbelow : ∀ x y, x ∈ readStack s.tstack →
      Items.Below ((WalkM.finishTstackTop item).run s).2.items x y → Items.Below s.items x y := by
    intro x y hx h
    rw [hitems] at h
    exact h.of_modify _ (hnb x hx)
  have hbelow' : ∀ x y, x ∈ readStack s.tstack →
      Items.Below s.items x y → Items.Below ((WalkM.finishTstackTop item).run s).2.items x y := by
    intro x y hx h
    rw [hitems]
    exact h.modify_of_not_below _ (hnb x hx)
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩
  · intro x hx p hpx
    rcases (hparent p x).1 hpx with ⟨rfl, hxC⟩ | ⟨hp, hpx'⟩
    · rcases (hmem' x).1 hx with rfl | hx
      · exact hitemC hxC
      · exact hdisj x hxC hx
    · rcases (hmem' x).1 hx with rfl | hx
      · exact hroot p hpx'
      · exact roots x hx p hpx'
  · rw [hts', hperm'.nodup_iff]
    exact List.nodup_cons.2 ⟨hitemR, hrestnd⟩
  · intro x hx y hy
    rw [hsize]
    rcases (hmem' x).1 hx with rfl | hx
    · rcases Relation.ReflTransGen.cases_head hy with rfl | ⟨c, hc, hcy⟩
      · exact hi
      · rcases (hparent x c).1 hc with ⟨_, hcC⟩ | ⟨h, _⟩
        · exact bounded c ((hmem c).2 (Or.inl hcC)) y (hbelow c y ((hmem c).2 (Or.inl hcC)) hcy)
        · exact absurd rfl h
    · exact bounded x ((hmem x).2 (Or.inr hx)) y (hbelow x y ((hmem x).2 (Or.inr hx)) hy)
  · intro p c hpc
    rw [hsize]
    rcases (hparent p c).1 hpc with ⟨_, hcC⟩ | ⟨_, hpc⟩
    · exact hCsize c hcC
    · exact chLt p c hpc
  · intro p
    by_cases hp : p = item
    · subst hp; rw [hch_item]; exact hCnd
    · rw [hch_ne p hp]; exact chNodup p
  · intro j hj hty
    rw [hsize] at hj
    rw [htype] at hty
    by_cases hji : j = item
    · subst hji
      exact Or.inl ⟨j, (hmem' j).2 (Or.inl rfl), .refl⟩
    rcases closed j hj hji hty with ⟨x, hx, hxj⟩ | ⟨b, hb, hbj⟩
    · left
      rcases (hmem x).1 hx with hxC | hxR
      · refine ⟨item, (hmem' item).2 (Or.inl rfl), .head ?_ (hbelow' x j hx hxj)⟩
        exact (hparent item x).2 (Or.inl ⟨rfl, hxC⟩)
      · exact ⟨x, (hmem' x).2 (Or.inr hxR), hbelow' x j hx hxj⟩
    · right
      refine ⟨b, hb, ?_⟩
      rw [hitems]
      exact hbj.modify_root _ hroot hji

end Spqr
