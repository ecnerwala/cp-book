import Spqr.WalkWF
import Spqr.PieceSep
import Spqr.Proofs.Dfs
import Spqr.WalkItemsWF
import Spqr.RelabelPieceSep
import Spqr.RangesPiece

namespace Spqr

/-- Admitted (named hypothesis for the walk induction, `WalkBackbone.lean`): the final items satisfy
`PieceFacts` (`RangesPiece.lean`). Within a tree the backbone carries `PieceInv d s` (`PieceFacts` plus
the path clause `root_path`), kept by every primitive (`PieceInv.alloc`/`modify`/`modifyVs`/`setVert`/
`exit`/`frame`/`pushVert`/`mergeTop`/`retarget`/…) and by `finishEdge_piece` from `CloseBase` +
`CloseContent` + `ClosePiece`; a root append keeps `PieceFacts` (`PieceFacts.rootAppend`) and the next
root re-enters `PieceInv 0` (`PieceFacts.root`). Checker `check_walkinv`: `checkPieceInv` at the tree
entry/end, before/after every `finishEdge` and (`PieceFacts`) after every root append, `checkPiece` for
`ClosePiece`, `Ranges.checkFinal` on the final items. -/
theorem walk_pieceInv (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    (g.walk tern (g.dfsForest vo eo)).PieceFacts := by
  sorry

theorem walk_g' (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) : (g.walk tern (g.dfsForest vo eo)).g = g := by
  rcases Nat.eq_zero_or_pos g.nv with hnv | hnv
  · rw [dfsForest_nil_of_nv_zero g hnv hvo]; rfl
  · exact walk_g g tern _ hnv (dfsForest_bounded g hg hvo heo)

theorem walk_q_upper (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    Items.QUpper g (g.walk tern (g.dfsForest vo eo)).items :=
  by have h := (walk_pieceInv g hg tern vo eo hvo heo).q_upper; rwa [walk_g' g hg tern vo eo hvo heo] at h

theorem walk_root_sep (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    Items.RootSep g (g.walk tern (g.dfsForest vo eo)).items :=
  by have h := (walk_pieceInv g hg tern vo eo hvo heo).root_sep; rwa [walk_g' g hg tern vo eo hvo heo] at h

theorem walk_root_v (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    Items.RootV (g.walk tern (g.dfsForest vo eo)).items :=
  (walk_pieceInv g hg tern vo eo hvo heo).root_v

theorem walk_q_child_vs (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    Items.QChildVs g (g.walk tern (g.dfsForest vo eo)).items :=
  by have h := (walk_pieceInv g hg tern vo eo hvo heo).q_child_vs; rwa [walk_g' g hg tern vo eo hvo heo] at h

theorem walk_p_child_vs (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    Items.PChildVs (g.walk tern (g.dfsForest vo eo)).items :=
  (walk_pieceInv g hg tern vo eo hvo heo).p_child_vs

/-- Walk-side separation of the uncapped pieces (F components, V blocks, root Qs). -/
theorem spqrTree_pieceSep (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    (g.spqrTree tern vo eo).PieceSep g := by
  rw [spqrTree_eq]
  obtain ⟨idx, hok⟩ := relabelOK_of_wf g _ (walk_items_wf g hg tern vo eo hvo heo)
  obtain ⟨σ, hr⟩ := walk_ranges' g hg tern vo eo hvo heo
  exact hok.pieceSep hg hr (walk_q_upper g hg tern vo eo hvo heo)
    (walk_root_sep g hg tern vo eo hvo heo) (walk_root_v g hg tern vo eo hvo heo)
    (walk_q_child_vs g hg tern vo eo hvo heo) (walk_p_child_vs g hg tern vo eo hvo heo)

end Spqr
