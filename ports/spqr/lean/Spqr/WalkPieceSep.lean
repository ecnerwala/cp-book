import Spqr.WalkWF
import Spqr.PieceSep
import Spqr.Proofs.Dfs
import Spqr.WalkItemsWF
import Spqr.RelabelPieceSep

namespace Spqr

/-- Admitted (named hypothesis for the walk induction; checker `check_walkinv` field
`ranges.q_upper` on the final items): a child of a V item has that vertex as first endpoint. -/
theorem walk_q_upper (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    Items.QUpper g (g.walk tern (g.dfsForest vo eo)).items := by
  sorry

/-- Admitted (named hypothesis for the walk induction; checker field `ranges.root_sep` on the
final items): distinct root children share no vertex. -/
theorem walk_root_sep (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    Items.RootSep g (g.walk tern (g.dfsForest vo eo)).items := by
  sorry

/-- Walk-side separation of the uncapped pieces (F components, V blocks, root Qs). -/
theorem spqrTree_pieceSep (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    (g.spqrTree tern vo eo).PieceSep g := by
  sorry

end Spqr
