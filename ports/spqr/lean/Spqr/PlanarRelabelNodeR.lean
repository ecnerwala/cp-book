import Spqr.PlanarRelabel
import Spqr.RelabelSpec
import Spqr.PlanarNodeSpec
import Spqr.Proofs.PlanarWalkFacts

/-!
# `RelabelNodeR`: what `planarRelabel` records for a planar R node

The relabel-side counterpart of `PlanarFinish`: for an R node `n` of the output tree there is an
item `it` of the walk, recorded `.planar m`, whose children `planarRelabel` laid out in
`Items.ordered` order with vertex positions `pos`; the node's ranges, skeleton and `neRotAdj`
rows are determined by these data and by the walk's `qem` after `setupNode` (`capLinked`) and
`applyFlips` (`flipped`). The rows are stated for the entries whose target edge is one of the
node's own virtual edges (where `rotEdgeNe` is the node's numbering `neSt + idxOf`); the
remaining entries are excluded by `PlanarFinish.closed`.
-/

namespace Spqr

open PlanarRelabelM

/-- The data `planarRelabel` leaves for the R node `n` of `T`, built from the walk state `w`. -/
def RelabelNodeR (g : Graph) (w : PlanarWalkState) (T : PlanarSpqrTree) (n : Nat) : Prop :=
  ∃ (it : ItemId) (m : Array Nat) (children : List ItemId) (pos : Nat → Nat),
    let items := w.base.items
    let nvSt := (T.toSpqrTree.nvRange n).1
    let nvEn := (T.toSpqrTree.nvRange n).2
    let neSt := (T.toSpqrTree.neRange n).1
    let neEn := (T.toSpqrTree.neRange n).2
    let edgeVes := (children.filter (· ≥ 1 + g.nv)).map (· - (1 + g.nv))
    let ves := capVe g :: edgeVes
    let Q := flipped g items[it]!.ch w.aux.itemFlips[it]! (capLinked g w.aux.qem m)
    1 + g.nv + g.ne ≤ it ∧ it < items.size ∧
    items[it]!.type = .R ∧
    w.aux.nodePlanarity[it - (1 + g.nv + g.ne)]? = some (.planar m) ∧
    children = Items.ordered g items it nvSt pos ∧
    Items.PosOK nvSt (Items.nvList g items it) pos ∧
    nvEn = nvSt + (Items.nvList g items it).length ∧
    neEn = neSt + edgeVes.length + 1 ∧
    T.toSpqrTree.skeleton n = (nvSt, nvEn - 1) :: Items.edgeChildren g items pos children ∧
    4 * neEn ≤ T.neRotAdj.size ∧
    ∀ l, l < 4 * ves.length →
      (Q[4 * ves[l / 4]! + l % 4]! = none → T.neRotAdj[4 * neSt + l]! = none) ∧
      ∀ o, Q[4 * ves[l / 4]! + l % 4]! = some o → QE.edge o ∈ ves →
        T.neRotAdj[4 * neSt + l]! =
          some (4 * (neSt + ves.idxOf (QE.edge o)) + (o &&& 2) + (1 - l % 2))

end Spqr
