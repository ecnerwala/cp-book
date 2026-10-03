import Spqr.PlanarInv

/-!
# `PlanarFinish`: the planarity record of a closed node, read at the end of the walk

When the planar walk closes an S/P/R item (`finishTstackTop`) it records
`nodePlanarity[item] = .planar m`, the four exposed ends of the entry's piece (`finishMatches`),
and the flips of the item's children (`itemFlips`). `planarRelabel` later reads the node's local
rotation system off `qem` (`setupNode` links the cap slots `8 ne + s` to `m[s]`, `applyFlips`
swaps the two corners at each endpoint of every flipped child edge, `mapRot` transports the
matches to node quarter-edges). `PlanarFinish` states, on the final walk state alone, that this
reading is a planar embedding of the node's skeleton: it is the walk-side fact behind
`nodePlanar_sound_R` (PROOF.md §8.4), to be carried by the backbone induction
(`Spqr/WalkBackbone.lean`) and checked by `checks/WalkInvCheck/Planar.lean`.
-/

namespace Spqr

/-- The walk-side virtual edge of the cap of every node (`setupNode` uses the slots `8 ne + s`). -/
def capVe (g : Graph) : Nat := 2 * g.ne

/-- The skeleton of the node item `it` as a `Piece` over the walk's virtual edges: the cap first,
then the virtual edges of the edge children in `ch` order, with the items' `vs` as endpoints. -/
def nodeSkel (g : Graph) (items : Array Item) (it : ItemId) : Piece where
  ves := capVe g :: (items[it]!.ch.filter (· ≥ 1 + g.nv)).map (· - (1 + g.nv))
  ends := fun ve =>
    let vs := if ve = capVe g then items[it]!.vs else items[1 + g.nv + ve]!.vs
    (vs.1.getD 0, vs.2.getD 0)
  nVerts := g.nv
  bot := 0
  top := 0

/-- `qem` after `setupNode`: the cap slots linked to the recorded matches `m`. -/
def capLinked (g : Graph) (qem : Qem) (m : Array Nat) : Qem :=
  (List.range 4).foldl (fun q s => (q.set! (4 * capVe g + s) (some m[s]!)).set! m[s]! (some (4 * capVe g + s))) qem

/-- `qem` after `applyFlips`: the two corners at each endpoint of every flipped edge child swapped. -/
def flipped (g : Graph) (ch : List ItemId) (flips : List Bool) (qem : Qem) : Qem :=
  (ch.zip flips).foldl (fun q (c, f) =>
    if c ≥ 1 + g.nv && f then
      let ve := c - (1 + g.nv)
      (q.swapIfInBounds (4 * ve + 0) (4 * ve + 1)).swapIfInBounds (4 * ve + 2) (4 * ve + 3)
    else q) qem

/-- The local rotation system `mapRot` reads for the piece `P` off `qem`: local quarter-edge
`4 k + z` of `P.ves[k]` faces the local quarter-edge of `qem[4 P.ves[k] + z]`, on the same side
and with the opposite direction. -/
def readRot (P : Piece) (qem : Qem) : RotationSystem :=
  ⟨(Array.range (4 * P.ves.length)).map fun l =>
    (qem[4 * P.ves[l / 4]! + l % 4]!).bind fun o =>
      (P.ves.findIdx? (· == QE.edge o)).map fun k => 4 * k + (o &&& 2) + (1 - l % 2)⟩

/-- The rotation system `planarRelabel` reads for the closed node `it` with matches `m`, from the
final walk state `w`: cap links, the children's flips, then `readRot`. -/
def finishRot (g : Graph) (w : PlanarWalkState) (it : ItemId) (m : Array Nat) : RotationSystem :=
  readRot (nodeSkel g w.base.items it)
    (flipped g w.base.items[it]!.ch w.aux.itemFlips[it]! (capLinked g w.aux.qem m))

/-- **Planar finish.** In the final walk state `w` of `g`, every S/P/R item recorded planar with
matches `m` has `m` of size four listing distinct quarter-edges of its edge children, and the
rotation system `planarRelabel` reads for it is a planar embedding of its skeleton. -/
structure PlanarFinish (g : Graph) (w : PlanarWalkState) : Prop where
  match_ends : ∀ k m, w.aux.nodePlanarity[k]? = some (.planar m) →
    let it := 1 + g.nv + g.ne + k
    w.base.items[it]!.type ∈ [NodeType.S, .P, .R] →
    m.size = 4 ∧ (∀ s, s < 4 → 1 + g.nv + QE.edge m[s]! ∈ w.base.items[it]!.ch) ∧
      m.toList.Nodup
  /-- Every child of the item is a V item or a walk virtual edge below the cap (`setupNode` and
  `applyFlips` index `qem`/`rotEdgeNe` by `c - (1 + nv)`). -/
  ch_lt : ∀ k m, w.aux.nodePlanarity[k]? = some (.planar m) →
    let it := 1 + g.nv + g.ne + k
    w.base.items[it]!.type ∈ [NodeType.S, .P, .R] →
    ∀ c ∈ w.base.items[it]!.ch, c < 1 + g.nv + 2 * g.ne
  /-- The `qem` entries `mapRot` reads (cap linked, flips applied) are quarter-edges of the skeleton. -/
  closed : ∀ k m, w.aux.nodePlanarity[k]? = some (.planar m) →
    let it := 1 + g.nv + g.ne + k
    w.base.items[it]!.type ∈ [NodeType.S, .P, .R] →
    ∀ ve ∈ (nodeSkel g w.base.items it).ves, ∀ z, z < 4 →
      ∃ o, (flipped g w.base.items[it]!.ch w.aux.itemFlips[it]! (capLinked g w.aux.qem m))[4 * ve + z]! =
        some o ∧ QE.edge o ∈ (nodeSkel g w.base.items it).ves
  embedded : ∀ k m, w.aux.nodePlanarity[k]? = some (.planar m) →
    let it := 1 + g.nv + g.ne + k
    w.base.items[it]!.type ∈ [NodeType.S, .P, .R] →
    IsPlanarEmbedding (nodeSkel g w.base.items it).es g.nv (finishRot g w it m)

end Spqr
