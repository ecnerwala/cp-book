import Spqr.RelabelRep
import Spqr.PieceSep
import Spqr.Ranges

/-!
# Relabel: `SpqrTree.PieceSep` from the final items

`spqrTree_pieceSep` (`WalkPieceSep.lean`) is `Items.WF` + `Items.Ranges` + two walk-side facts about
the final items (`Items.QUpper`, `Items.RootSep`), transported through `RelabelOK` (`RelabelRep.lean`).
-/

namespace Spqr.Items

variable (g : Graph) (items : Items)

/-- A child of a V item records that vertex as its first endpoint (`finishBoundary` writes
`vs := (some curV, none)` and appends the block root to `vertItem curV`). -/
def QUpper : Prop :=
  ∀ v c, v < g.nv → items.IsParent (vertItem v) c → (items.vs c).1 = some v

/-- Distinct children of the root (the DFS components) share no vertex. -/
def RootSep : Prop :=
  ∀ a b, items.IsParent rootItem a → items.IsParent rootItem b → a ≠ b →
    ∀ v e e', e < g.ne → e' < g.ne → g.Inc e v → g.Inc e' v →
      items.EdgeBelow g a e → items.EdgeBelow g b e' → False

end Spqr.Items
