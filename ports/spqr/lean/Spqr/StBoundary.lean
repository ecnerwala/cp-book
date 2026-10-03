import Spqr.StTree
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

theorem spansCount_pos_of_mem_readStack : ∀ {ts : List TEntry} {x : ItemId}, x ∈ readStack ts →
    0 < spansCount ts x
  | [], _, h => by simp [readStack, readL, readR] at h
  | t :: ts, x, h => by
    rw [spansCount_cons]
    rcases mem_readStack_cons.1 h with h | h
    · exact Nat.lt_of_lt_of_le (List.count_pos_iff.2 h) (Nat.le_add_right _ _)
    · exact Nat.lt_of_lt_of_le (spansCount_pos_of_mem_readStack h) (Nat.le_add_left _ _)

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
    StItems g r.2 (blocks ++ (if o.cls.isTree then [⟨some (curV, o.dest), stNest ps⟩] else [])) := by
  sorry

/-- Admitted: popping the root entry of a tree of the forest onto `rootItem` closes the root block
`⟨none, stNest ps⟩` (`ps` the pieces of the whole tree): every S / P / R item of the entry becomes
`InBlock` of it. Checked after every tree by `check_stsim`. -/
theorem rootPop_st {g : Graph} {s : WalkState} {t : TEntry} {ps : List StPiece} {blocks : List StBlock}
    (hts : s.tstack = [t]) (hg : s.g = g) (hi : s.Inv' 0) (hs : Shape s) (hrk : RootOK s)
    (hroot : ∀ p, ¬ Items.IsParent s.items p rootItem)
    (hR : StRead s.items [t] ps) (hI : StItems g s blocks) :
    StItems g
      { s with items := s.items.modify rootItem fun it => { it with ch := it.ch ++ t.spans.2 },
               tstack := [] }
      (blocks ++ [⟨none, stNest ps⟩]) := by
  sorry

end Spqr
