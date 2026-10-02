import Spqr.Spec

/-!
# Children shapes the gluing step relies on

`SpqrTree.WF` constrains the node-vertex layout and the twin edges but not which children an
`O`/`I` leaf may have (`nv_layout` allows a `V` child under an `O` item). `embedItem` relies on the
phase-2 children shapes (`Items.Shapes`); `ChildShape` restates the ones it needs on the output
tree (derived from `Items.Shapes` through the relabel interface in `RelabelChildShape.lean`).
-/

namespace Spqr

namespace SpqrTree

variable (t : SpqrTree)

structure ChildShape : Prop where
  /-- `O` and `I` items are leaves. -/
  leaf : ∀ i, i < t.size → t.type i = .O ∨ t.type i = .I → t.children i = []

end SpqrTree

end Spqr
