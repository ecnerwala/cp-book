import Spqr.StTree
import Spqr.StBdPop
import Spqr.StRetFrame
import Spqr.StVStart
import Spqr.WalkCover
import Spqr.WalkInv

/-!
# The block-completing steps of the simulation (PROOF.md §7.6)

`finishEdge` at a block boundary (`d ≤ lowval`) pops the entries of the finished block and makes
their items children of the boundary edge's Q item; popping the root entry of a tree of the forest
closes the root block. Both steps turn the live S / P / R items of the popped entries into
`InBlock` items of the new block (`StItems.closed`). The bookkeeping (`BdPop`, `StRead.complete`,
`StBdPop.lean`) is proved; what remains admitted is `finishBoundary_stLive`, the open-block invariant
`StLive` at the boundary (the orientation clauses of the new block), whose preservation along the
walk is PROOF.md §7.6, step (c). `Place.fresh` gives the freshness of a not-yet-pushed fixed item from the
coverage invariant.
-/

namespace Spqr

open WalkM WalkState

/-- A fixed item not yet pushed is a root and is off the stack. -/
theorem Place.fresh {g : Graph} {P X : ItemId → Prop} {s : WalkState} (h : s.Place g P X) {i : ItemId}
    (h0 : 0 < i) (hi : i < 1 + g.nv + g.ne) (hP : ¬ P i) :
    (∀ p, ¬ Items.IsParent s.items p i) ∧ i ∉ readStack s.tstack := by
  have hc := h.cnt_eq_zero h0 hi hP
  refine ⟨fun p => noParent_of_cnt_eq_zero hc p, fun hm => ?_⟩
  have := spansCount_pos_of_mem_readStack hm
  simp only [WalkState.cnt] at hc; omega

/-- Admitted: the open-block invariant `StLive` (StBdPop.lean) at a boundary tree edge, on the
pre-state: every S / P / R item below a span item of the popped entries `sub` (not hanging under a
V / Q item, i.e. not in an already completed block) is `InBlock` of the block
`⟨some (curV, o.dest), stNest ps⟩` closed there — its leaves are a segment of `stNest ps` (trivial
from `StRead`), and the orientation clauses `VsOrientedAt` hold. This is what `check_stsim` checks
as m8 (against the blocks of the truncated reference `refBlocks g (prev ++ [truncTree fs t])`,
seeds 0..1000, 0 violations); it is the live hypothesis `hlive` of `finishBoundary_st`, to be
supplied by the backbone induction carrying `StLive` per open segment (PROOF.md §7 "StLive
obligations"). -/
theorem finishBoundary_stLive {D curV d : Nat} {o : DfsOut} {hasVert : Bool} {s : WalkState}
    {sub base : List TEntry} {g : Graph} {ps : List StPiece} {blocks : List StBlock}
    (hE : s.EarFinish curV d o hasVert sub base) (hi : s.Inv' D) (hs : Shape s)
    (hD : D = if o.cls.isTree then d + 1 else d) (hok : BoundaryOk D curV d o s)
    (hb : FinishBook curV d o base.length hasVert s) (hge : d ≤ o.cls.lowval d)
    (ht : o.cls.isTree = true) (hg : s.g = g)
    (hR : StRead s.items sub ps) (hI : StItems g s blocks) :
    StLive g s.items sub ⟨some (curV, o.dest), stNest ps⟩ := by
  sorry

/-- `finishEdge` at a block boundary. The popped entries `sub` read as the pieces `ps` of the
child's subtree (`[]` for a back edge; `[t]` at a bridge, `[t₁, t₂]` at a component edge:
`EarFinish.bd_bridge`/`bd_comp`); the edge's Q item becomes the root of the block
`⟨some (curV, o.dest), stNest ps⟩`, and every S / P / R item of `sub` becomes `InBlock` of it
(`StRead.complete`, with the orientation clauses `finishBoundary_vsOrientedAt`) or is already
finished (`StItems.finished`). The tstack is `base` afterwards, `stackDir` and the returned
`hasVert` are unchanged, and the items below `base` keep type and children (`BdPop`). -/
theorem finishBoundary_st {D curV d : Nat} {o : DfsOut} {hasVert : Bool} {s : WalkState}
    {sub base : List TEntry} {g : Graph} {ps : List StPiece} {blocks : List StBlock}
    (hE : s.EarFinish curV d o hasVert sub base) (hi : s.Inv' D) (hs : Shape s)
    (hD : D = if o.cls.isTree then d + 1 else d) (hok : BoundaryOk D curV d o s)
    (hb : FinishBook curV d o base.length hasVert s) (hge : d ≤ o.cls.lowval d) (hg : s.g = g)
    (hR : StRead s.items sub ps) (hI : StItems g s blocks)
    (hlive : o.cls.isTree = true → StLive g s.items sub ⟨some (curV, o.dest), stNest ps⟩) :
    let r := (finishEdge curV d o base.length hasVert).run s
    r.1 = hasVert ∧ r.2.g = s.g ∧ r.2.stackDir = s.stackDir ∧ r.2.tstack = base ∧
    StItems g r.2 (blocks ++ (if o.cls.isTree then [⟨some (curV, o.dest), stNest ps⟩] else [])) ∧
    ∀ x ∈ readStack base, ∀ y, Items.Below s.items x y →
      Items.type r.2.items y = Items.type s.items y ∧ Items.ch r.2.items y = Items.ch s.items y ∧
      Items.vs r.2.items y = Items.vs s.items y := by
  have hq_root := hok.q_root
  have hv_root := hok.v_root
  have hq_free : edgeItem s.g o.e ∉ readStack s.tstack := fun h => by
    obtain ⟨u, hu, hx⟩ := mem_readStack_exists h; exact hok.q_free u hu hx
  have hv_free : vertItem curV ∉ readStack s.tstack := fun h => by
    obtain ⟨u, hu, hx⟩ := mem_readStack_exists h; exact hok.v_free u hu hx
  have hq_ty : Items.type s.items (edgeItem s.g o.e) = .Q := hs.edge o.e hb.e_lt
  have hv_ty : Items.type s.items (vertItem curV) = .V := hs.vert curV hb.v_lt
  have hq_lt : edgeItem s.g o.e < s.items.size := by
    have h1 := hs.size; have h2 := hb.e_lt; show 1 + s.g.nv + o.e < s.items.size; omega
  have hv_lt : vertItem curV < s.items.size := by
    have h1 := hs.size; have h2 := hb.v_lt; show 1 + curV < s.items.size; omega
  have hqv : edgeItem s.g o.e ≠ vertItem curV := fun e => by rw [e, hv_ty] at hq_ty; cases hq_ty
  have hframe : ∀ {items' : Items}, BdPop s sub (edgeItem s.g o.e) (vertItem curV) items' →
      ∀ x ∈ readStack base, ∀ y, Items.Below s.items x y →
      Items.type items' y = Items.type s.items y ∧ Items.ch items' y = Items.ch s.items y ∧
      Items.vs items' y = Items.vs s.items y := by
    intro items' hP x hx y hy
    have hxs : x ∈ readStack s.tstack := by rw [hE.tstack]; exact mem_readStack_append.2 (.inr hx)
    have hyq : y ≠ edgeItem s.g o.e := fun e =>
      Items.not_below_root_of_ne hq_root (fun e' => hq_free (e' ▸ hxs)) (e ▸ hy)
    have hyv : y ≠ vertItem curV := fun e =>
      Items.not_below_root_of_ne hv_root (fun e' => hv_free (e' ▸ hxs)) (e ▸ hy)
    exact hP.frame y hyq hyv (hI.bounded x hxs y hy)
  have hge' : o.cls.lowval d ≥ d := hge
  show wp (finishEdge curV d o base.length hasVert) (fun r s' => r = hasVert ∧ s'.g = s.g ∧
    s'.stackDir = s.stackDir ∧ s'.tstack = base ∧
    StItems g s' (blocks ++ (if o.cls.isTree then [⟨some (curV, o.dest), stNest ps⟩] else [])) ∧
    ∀ x ∈ readStack base, ∀ y, Items.Below s.items x y →
      Items.type s'.items y = Items.type s.items y ∧ Items.ch s'.items y = Items.ch s.items y ∧
      Items.vs s'.items y = Items.vs s.items y) s
  simp only [finishEdge_eq, finishEdge', wp_bind, wp_get, wp_stackDir, hge', ↓reduceIte]
  unfold finishBoundary
  by_cases ht : o.cls.isTree = true
  · have hbd : ∀ i, (∃ x ∈ readStack sub, Items.Below s.items x i) →
        Items.type s.items i = .S ∨ Items.type s.items i = .P ∨ Items.type s.items i = .R →
        ∃ b ∈ blocks ++ [⟨some (curV, o.dest), stNest ps⟩], InBlock g s.items b i :=
      fun i hx hty => hR.complete rfl hI.finished
        (fun i L hx hty hL => (hlive ht).vsOrientedAt hx hty hL)
        hx hty
    by_cases hl : o.cls.lowval d = d + 1
    · obtain ⟨t, rfl, -, -⟩ := hE.bd_bridge ht hl
      have hts : s.tstack = t :: base := hE.tstack
      have ht_mem : t ∈ s.tstack := by rw [hts]; exact List.mem_cons_self
      simp only [wp_bind, wp_modifyItem, wp_modify, wp_allocItem, wp_makeVs, wp_popTstack, wp_pure, ht,
        hl, beq_self_eq_true, ↓reduceIte]
      rw [hts]
      simp only [List.head!_cons, List.tail_cons, Array.size_modify, true_and]
      refine ⟨StItems.bdPop hI hE.tstack rfl ?hP hq_root hq_free hv_root hv_free hq_ty hv_ty
        hq_lt hv_lt hbd, hframe ?hP⟩
      refine BdPop.mk' hqv hq_lt hv_lt ?_ ?_ ?_ ?_ ?_
      · simp only [Array.size_modify, Array.size_push]; omega
      · intro y hy
        have hne : y ≠ s.items.size := Nat.ne_of_lt hy
        refine ⟨?_, ?_, fun hyq => ?_⟩
        · simp [Items.type, Array.getElem?_modify, Array.getElem?_push, Array.size_modify, hne, hne.symm]
          cases s.items[y]? <;> simp
          split <;> simp
        · simp [Items.ch, Array.getElem?_modify, Array.getElem?_push, Array.size_modify, hne, hne.symm]
          cases s.items[y]? <;> simp
          split <;> simp
        · simp [Items.vs, Array.getElem?_modify, Array.getElem?_push, Array.size_modify, hne, hne.symm,
            hyq.symm]
      · intro y hy
        have hnone : s.items[y]? = none := Array.getElem?_eq_none hy
        by_cases hys : y = s.items.size
        · subst hys
          simp [Items.type, Items.ch, Array.size_modify,
            Array.getElem_modify, Array.getElem_push]
        · simp [Items.type, Items.ch, Array.getElem?_modify, Array.getElem?_push, Array.size_modify, hnone, hys,
            Ne.symm hys]
      · intro c hc
        rcases List.mem_cons.1 hc with rfl | hc
        · exact .inr ⟨Nat.le_refl _, by simp [Array.size_modify, Array.size_push]⟩
        · exact .inl (mem_readStack_of_mem List.mem_cons_self (List.mem_append_right _ hc))
      · refine List.nodup_cons.2 ⟨fun h => ?_, ?_⟩
        · have := hI.bounded _ (mem_readStack_of_mem ht_mem (List.mem_append_right _ h)) _ .refl
          exact Nat.lt_irrefl _ this
        · simpa using readStack_nodup_mix hI.nodup ht_mem ht_mem |>.sublist (List.sublist_append_right _ _)
    · obtain ⟨t₁, t₂, rfl, -, -, -, -⟩ := hE.bd_comp ht hge hl
      have hts : s.tstack = t₁ :: t₂ :: base := hE.tstack
      have hl' : (o.cls.lowval d == d + 1) = false := beq_eq_false_iff_ne.2 hl
      have h₁ : t₁ ∈ s.tstack := by rw [hts]; exact List.mem_cons_self
      have h₂ : t₂ ∈ s.tstack := by rw [hts]; exact List.mem_cons_of_mem _ List.mem_cons_self
      simp only [wp_bind, wp_modifyItem, wp_modify, wp_popTstack, wp_pure, ht,
        hl', Bool.false_eq_true, ↓reduceIte]
      rw [hts]
      simp only [List.head!_cons, List.tail_cons, true_and]
      suffices hP : BdPop s [t₁, t₂] (edgeItem s.g o.e) (vertItem curV)
          (((s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) }).modify
            (edgeItem s.g o.e) fun it => { it with ch := t₁.spans.1 ++ t₂.spans.2 }).modify
            (vertItem curV) fun it => { it with ch := it.ch ++ [edgeItem s.g o.e] }) by
        exact ⟨StItems.bdPop hI hE.tstack rfl hP hq_root hq_free hv_root hv_free hq_ty hv_ty hq_lt hv_lt hbd,
          hframe hP⟩
      refine BdPop.mk' hqv hq_lt hv_lt ?_ ?_ ?_ ?_ ?_
      · simp only [Array.size_modify]; exact Nat.le_refl _
      · intro y hy
        refine ⟨?_, ?_, fun hyq => ?_⟩
        · simp [Items.type, Array.getElem?_modify]
          cases s.items[y]? <;> simp
          split <;> simp
        · simp [Items.ch, Array.getElem?_modify]
          cases s.items[y]? <;> simp
          split <;> simp
        · simp [Items.vs, Array.getElem?_modify, hyq.symm]
      · intro y hy
        have hnone : s.items[y]? = none := Array.getElem?_eq_none hy
        simp [Items.type, Items.ch, Array.getElem?_modify, hnone]
      · intro c hc
        rcases List.mem_append.1 hc with hc | hc
        · exact .inl (mem_readStack_of_mem List.mem_cons_self (List.mem_append_left _ hc))
        · exact .inl (mem_readStack_of_mem (List.mem_cons_of_mem _ List.mem_cons_self) (List.mem_append_right _ hc))
      · exact readStack_nodup_mix hI.nodup h₁ h₂
  · have ht' : o.cls.isTree = false := Bool.eq_false_iff.2 ht
    obtain rfl := hE.back_nil ht'
    have hts : s.tstack = base := hE.tstack
    have hbd : ∀ i, (∃ x ∈ readStack ([] : List TEntry), Items.Below s.items x i) →
        Items.type s.items i = .S ∨ Items.type s.items i = .P ∨ Items.type s.items i = .R →
        ∃ b ∈ blocks ++ [], InBlock g s.items b i :=
      fun _ ⟨_, hx, _⟩ _ => absurd hx (by simp [readStack, readL, readR])
    simp only [wp_bind, wp_modifyItem, wp_modify, wp_allocItem, wp_pure, ht',
      Bool.false_eq_true, ↓reduceIte, Array.size_modify, true_and]
    rw [hts]
    refine ⟨rfl, StItems.bdPop hI hE.tstack rfl ?hPb hq_root hq_free hv_root hv_free hq_ty hv_ty
      hq_lt hv_lt hbd, hframe ?hPb⟩
    refine BdPop.mk' hqv hq_lt hv_lt ?_ ?_ ?_ ?_ ?_
    · simp only [Array.size_modify, Array.size_push]; omega
    · intro y hy
      have hne : y ≠ s.items.size := Nat.ne_of_lt hy
      refine ⟨?_, ?_, fun hyq => ?_⟩
      · simp [Items.type, Array.getElem?_modify, Array.getElem?_push, Array.size_modify, hne, hne.symm]
        cases s.items[y]? <;> simp
        split <;> simp
      · simp [Items.ch, Array.getElem?_modify, Array.getElem?_push, Array.size_modify, hne, hne.symm]
        cases s.items[y]? <;> simp
        split <;> simp
      · simp [Items.vs, Array.getElem?_modify, Array.getElem?_push, Array.size_modify, hne, hne.symm,
          hyq.symm]
    · intro y hy
      have hnone : s.items[y]? = none := Array.getElem?_eq_none hy
      by_cases hys : y = s.items.size
      · subst hys
        simp [Items.type, Items.ch, Array.size_modify,
          Array.getElem_modify, Array.getElem_push]
      · simp [Items.type, Items.ch, Array.getElem?_modify, Array.getElem?_push, Array.size_modify, hnone, hys,
          Ne.symm hys]
    · intro c hc
      rw [List.mem_singleton] at hc; subst hc
      exact .inr ⟨Nat.le_refl _, by simp [Array.size_modify, Array.size_push]⟩
    · exact List.nodup_singleton _

/-- Frame facts of `finishEdge` on a returning edge (`lowval = lv < d`), the hypotheses of
`finishEdge_st`: for `hasVert = false` the vertex item is still a root off the stack after `finishP`
(the `hvf` of `finishTree_st`/`finishBack_st`; from `EarFinish.v_root`/`vert_free` by per-primitive
frames), and the items below the untouched base `B` keep their type and children. -/
theorem finishRet_frame_st {g : Graph} {P X : ItemId → Prop} {blocks : List StBlock} {d lv : Nat}
    {kind : RetKind} {o : DfsOut} {s : WalkState} {curV : Nat} {hasVert : Bool} {sub pre B : List TEntry}
    (hE : s.EarFinish curV d o hasVert sub (pre ++ B)) (hfull : s.Full g P X) (hI : StItems g s blocks)
    (hv : curV < s.g.nv) (ho : o.cls = .ret lv kind) (hlow : lv < d) (hB : ∀ t ∈ B, t.vStart ≠ curV) :
    (hasVert = false →
      (∀ p, ¬ Items.IsParent (fePState curV lv d o s).items p (vertItem curV)) ∧
      vertItem curV ∉ readStack (fePState curV lv d o s).tstack) ∧
    ∀ x ∈ readStack B, ∀ y, Items.Below s.items x y →
      Items.type ((finishEdge curV d o (pre ++ B).length hasVert).run s).2.items y = Items.type s.items y ∧
      Items.ch ((finishEdge curV d o (pre ++ B).length hasVert).run s).2.items y = Items.ch s.items y ∧
      Items.vs ((finishEdge curV d o (pre ++ B).length hasVert).run s).2.items y = Items.vs s.items y :=
  finishRet_frame hE hfull hI hv ho hlow hB

/-- Popping the root entry of a tree of the forest onto `rootItem` closes the root block
`⟨none, stNest ps⟩` (`ps` the pieces of the whole tree). The entry holds only vertex items
(`RootOK`), so every S / P / R item below it is already finished (`StItems.finished`): the new block
has no S / P / R items of its own, and the old blocks are kept by the root's children change. -/
theorem rootPop_st {g : Graph} {s : WalkState} {t : TEntry} {ps : List StPiece} {blocks : List StBlock}
    (hts : s.tstack = [t]) (hg : s.g = g) (hi : s.Inv' 0) (hs : Shape s) (hrk : RootOK s)
    (hroot : ∀ p, ¬ Items.IsParent s.items p rootItem)
    (hR : StRead s.items [t] ps) (hI : StItems g s blocks) :
    StItems g
      { s with items := s.items.modify rootItem fun it => { it with ch := it.ch ++ t.spans.2 },
               tstack := [] }
      (blocks ++ [⟨none, stNest ps⟩]) := by
  obtain ⟨t', ht', hs1, hs2⟩ := hrk
  rw [hts] at ht'
  have e := (List.cons.inj ht').1
  subst e
  set f : Item → Item := fun it => { it with ch := it.ch ++ t.spans.2 } with hf
  have hr0 : rootItem < s.items.size := by have := hs.size; show 0 < s.items.size; omega
  have hread : readStack s.tstack = t.spans.2 := by rw [hts]; simp [readStack, readL, readR, hs1]
  have hty : ∀ j, Items.type (s.items.modify rootItem f) j = Items.type s.items j :=
    fun j => Items.type_modify _ _ _ _ (fun _ => rfl)
  have hch_ne : ∀ p, p ≠ rootItem → Items.ch (s.items.modify rootItem f) p = Items.ch s.items p :=
    fun p hp => Items.ch_modify_ne _ _ _ _ (Ne.symm hp)
  have hch0 : Items.ch (s.items.modify rootItem f) rootItem = Items.ch s.items rootItem ++ t.spans.2 := by
    rw [Items.ch_modify_self _ _ _ hr0]; simp [f, Items.ch, Array.getElem?_eq_getElem hr0]
  have hneF : ∀ x, Items.type s.items x ≠ .F → x ≠ rootItem := fun x hx e => hx (e ▸ hs.root)
  refine ⟨fun x hx => by simp [readStack, readL, readR] at hx, by simp [readStack, readL, readR],
    fun x hx => by simp [readStack, readL, readR] at hx, fun p c hpc => ?_, fun p => ?_,
    fun i hi hsp => ?_, fun x i hVQ hb hsp => ?_⟩
  · change c ∈ Items.ch (s.items.modify rootItem f) p at hpc
    change c < (s.items.modify rootItem f).size
    rw [Array.size_modify]
    by_cases hp : p = rootItem
    · subst hp
      rw [hch0] at hpc
      rcases List.mem_append.1 hpc with h | h
      · exact hI.chLt _ _ h
      · obtain ⟨v, hv, rfl⟩ := hs2 c h
        have := hs.size
        show 1 + v < _
        omega
    · rw [hch_ne p hp] at hpc
      exact hI.chLt p c hpc
  · change (Items.ch (s.items.modify rootItem f) p).Nodup
    by_cases hp : p = rootItem
    · subst hp
      rw [hch0]
      refine List.Nodup.append (hI.chNodup rootItem) (by rw [← hread]; exact hI.nodup) fun c hc hc' => ?_
      exact hI.roots c (by rw [hread]; exact hc') rootItem hc
    · rw [hch_ne p hp]
      exact hI.chNodup p
  · change i < (s.items.modify rootItem f).size at hi
    change Items.type (s.items.modify rootItem f) i = .S ∨ Items.type (s.items.modify rootItem f) i = .P ∨
      Items.type (s.items.modify rootItem f) i = .R at hsp
    rw [Array.size_modify] at hi
    rw [hty] at hsp
    have hi0 : i ≠ rootItem := hneF i fun e => by rw [e] at hsp; rcases hsp with h | h | h <;> cases h
    have hIB : ∃ b ∈ blocks, InBlock g s.items b i := by
      rcases hI.closed i hi hsp with ⟨x, hx, hb⟩ | h
      · rw [hread] at hx
        obtain ⟨v, hv, rfl⟩ := hs2 x hx
        exact hI.finished _ i (Or.inl (hs.vert v hv)) hb hsp
      · exact h
    obtain ⟨b, hb, hB⟩ := hIB
    exact Or.inr ⟨b, List.mem_append_left _ hb, hB.modify_root f hroot hi0⟩
  · change Items.type (s.items.modify rootItem f) x = .V ∨ Items.type (s.items.modify rootItem f) x = .Q at hVQ
    change Items.type (s.items.modify rootItem f) i = .S ∨ Items.type (s.items.modify rootItem f) i = .P ∨
      Items.type (s.items.modify rootItem f) i = .R at hsp
    change Items.Below (s.items.modify rootItem f) x i at hb
    rw [hty] at hVQ hsp
    have hx0 : x ≠ rootItem := hneF x fun e => by rw [e] at hVQ; rcases hVQ with h | h <;> cases h
    have hi0 : i ≠ rootItem := hneF i fun e => by rw [e] at hsp; rcases hsp with h | h | h <;> cases h
    have hb' := hb.of_modify f (Items.not_below_root_of_ne hroot hx0)
    obtain ⟨b, hb, hB⟩ := hI.finished x i hVQ hb' hsp
    exact ⟨b, List.mem_append_left _ hb, hB.modify_root f hroot hi0⟩

end Spqr
