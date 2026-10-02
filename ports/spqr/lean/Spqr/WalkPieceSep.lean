import Spqr.WalkWF
import Spqr.PieceSep
import Spqr.Proofs.Dfs

namespace Spqr

/-- Walk-side separation of the uncapped pieces (F components, V blocks, root Qs). -/
theorem spqrTree_pieceSep (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    (g.spqrTree tern vo eo).PieceSep g := by
  sorry

end Spqr
