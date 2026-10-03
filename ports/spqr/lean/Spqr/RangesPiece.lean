import Spqr.RelabelPieceSep
import Spqr.RangesCanon

/-! # The PieceSep item facts as a walk-state invariant (`PieceFacts`, `PieceInv`, `ClosePiece`)

`spqrTree_pieceSep` needs five facts about the final items (`Items.QUpper`/`RootSep`/`RootV`/
`QChildVs`/`PChildVs`, `RelabelPieceSep.lean`). `PieceFacts s` states them for the items of a walk
state; `PieceInv d s` adds the frame clause `root_path` (no vertex item of the open DFS path lies
below a root child), which is what keeps `RootSep` through the boundary writes. They are kept by
every walk primitive given the per-site record `ClosePiece` (the orientation of the closed node at
a boundary and `FinishPiece` at the three `finishTstackTop` sites: a P's new children carry its
`vs`); `finishEdge_piece` is the step lemma the walk induction (`WalkBackbone.lean`) uses, and the
executable checker (`checks/WalkInvCheck/Ranges.lean`, `checkPieceInv`/`checkPiece`) evaluates
every clause (`piece_*`). -/

namespace Spqr.WalkState
open WalkM

/-- The five final-items facts of `spqrTree_pieceSep` at a walk state. -/
structure PieceFacts (s : WalkState) : Prop where
  q_upper : Items.QUpper s.g s.items
  root_sep : Items.RootSep s.g s.items
  root_v : Items.RootV s.items
  q_child_vs : Items.QChildVs s.g s.items
  p_child_vs : Items.PChildVs s.items

/-- `PieceFacts` plus the frame clause: no vertex item of the open path `stackVerts[0..d]` lies below
a root child (root children are the finished DFS trees). -/
structure PieceInv (d : Nat) (s : WalkState) : Prop extends PieceFacts s where
  root_path : ∀ a, Items.IsParent s.items rootItem a → ∀ k, k ≤ d →
    ¬ Items.Below s.items a (vertItem s.stackVerts[k]!)

/-- `finishTstackTop x` from `r` keeps `PChildVs`: if `x` is a P, every non-V item on the closing
side of the top entry `t` already carries `makeVs t.vStart t.topDepth`, the `vs` the close writes
to `x`. -/
def FinishPiece (x : ItemId) (r : WalkState) : Prop :=
  Items.type r.items x = .P →
    ∀ c ∈ getSide (curE r).spans r.stackDir[(curE r).topDepth]!, Items.type r.items c ≠ .V →
      Items.vs r.items c =
        setSides r.stackDir[(curE r).topDepth]! (some r.stackVerts[(curE r).topDepth]!) (some (curE r).vStart)

/-- The PieceSep site facts of one `finishEdge curV d o origTstack hasVert` call from `s`: at a
boundary the new I node (`bd_bridge`) or the closed block node (`bd_node`) is oriented
`(curV, dest)`; `FinishPiece` at the P close (`finishP`, whenever `condP` holds at `feRest`), the
type-1 vertex close (`cvS₅`) and every loop-1 iteration whose first `k` conditions held. -/
structure ClosePiece (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool)
    (s : WalkState) : Prop where
  bd_bridge : o.cls.isTree = true → d ≤ o.cls.lowval d → o.cls.lowval d = d + 1 →
    setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) = (some curV, some o.dest)
  bd_node : o.cls.isTree = true → d ≤ o.cls.lowval d → o.cls.lowval d ≠ d + 1 →
    ∀ b ∈ s.tstack.head?, ∀ c, b.spans.1 = [c] → Items.vs s.items c = (some curV, some o.dest)
  p_site : o.cls.lowval d < d → o.cls.isType1 = true →
    result (condP curV (o.cls.lowval d) true) (feRest curV d o origTstack hasVert s) = true →
    FinishPiece (pItem (feRest curV d o origTstack hasVert s)) (pPre (feRest curV d o origTstack hasVert s))
  v_site : o.cls.isTree = true → o.cls.lowval d < d → hasVert = true → o.cls.isType1 = true →
    FinishPiece ((maybeUnwrapNxt (if feSingle d o s then NodeType.S else .R)).run (feS₂ d o s)).1
      (cvS₅ curV s.stackDir[d]! true origTstack (feSingle d o s) (feS₂ d o s))
  l1_site : o.cls.isTree = true → o.cls.lowval d < d → ∀ k,
    (∀ j, j ≤ k → result (loop1Cond d) (l1Iter d o s j) = true) →
    FinishPiece (l1Node d o s k) (l1Pre d o s k)

end Spqr.WalkState
