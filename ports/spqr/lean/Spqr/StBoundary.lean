import Spqr.StTree
import Spqr.StRetFrame
import Spqr.StVStart
import Spqr.WalkCover
import Spqr.WalkInv

/-!
# The block-completing steps of the simulation (PROOF.md §7.6)

`finishEdge` at a block boundary (`d ≤ lowval`) pops the entries of the finished block and makes
their items children of the boundary edge's Q item; popping the root entry of a tree of the forest
closes the root block. Both steps turn the live S / P / R items of the popped entries into
`InBlock` items of the new block (`StItems.closed`), which needs the `vs`-position invariant of
the stack spans (PROOF.md §7.6, step (c)); they are the two admissions of the simulation.
`Place.fresh` gives the freshness of a not-yet-pushed fixed item from the coverage invariant.
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

/-- Admitted: `finishEdge` at a block boundary. The popped entries `sub` read as the pieces `ps`
of the child's subtree (`[]` for a back edge); the edge's Q item becomes the root of the block
`⟨some (curV, o.dest), stNest ps⟩` (the block's st-order), and every S / P / R item of `sub`
becomes `InBlock` of it (the orientation clauses of `VsOrientedAt` need the `vs`-position
invariant of PROOF.md §7.6 (c)). The tstack is `base` afterwards, `stackDir` and the returned
`hasVert` are unchanged. Checked at every boundary edge by `check_stsim`. -/
theorem finishBoundary_st {D curV d : Nat} {o : DfsOut} {hasVert : Bool} {s : WalkState}
    {sub base : List TEntry} {g : Graph} {ps : List StPiece} {blocks : List StBlock}
    (hE : s.EarFinish curV d o hasVert sub base) (hi : s.Inv' D) (hs : Shape s)
    (hD : D = if o.cls.isTree then d + 1 else d) (hok : BoundaryOk D curV d o s)
    (hb : FinishBook curV d o base.length hasVert s) (hge : d ≤ o.cls.lowval d) (hg : s.g = g)
    (hR : StRead s.items sub ps) (hI : StItems g s blocks) :
    let r := (finishEdge curV d o base.length hasVert).run s
    r.1 = hasVert ∧ r.2.g = s.g ∧ r.2.stackDir = s.stackDir ∧ r.2.tstack = base ∧
    StItems g r.2 (blocks ++ (if o.cls.isTree then [⟨some (curV, o.dest), stNest ps⟩] else [])) ∧
    ∀ x ∈ readStack base, ∀ y, Items.Below s.items x y →
      Items.type r.2.items y = Items.type s.items y ∧ Items.ch r.2.items y = Items.ch s.items y := by
  sorry

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
      Items.ch ((finishEdge curV d o (pre ++ B).length hasVert).run s).2.items y = Items.ch s.items y :=
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
