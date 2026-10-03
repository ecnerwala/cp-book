import WalkInvCheck.Common
import Spqr.Proofs.PlanarWalkFacts
/-!
# `PlanarFinish` (`Spqr/Proofs/PlanarWalkFacts.lean`) on the final planar walk state

Separately reported fields, per S/P/R item recorded `.planar m`:
* `pf.size`: `m.size = 4`; `pf.mem`: each `m[s]` is a quarter-edge of an edge child;
  `pf.nodup`: the four matches are distinct; `pf.ch_lt`: children below `1 + nv + 2 ne`;
* `pf.closed`: the `qem` entries of the skeleton's quarter-edges (cap linked, flips applied) point
  back into the skeleton;
* `pf.planar`: `finishRot` is a planar embedding of the skeleton (`PlanarFinish.embedded`).

The children's flips are part of the fact: without `flipped` (reading `readRot P (capLinked g qem m)`
directly) the embedding test fails on 530 of the 6002 runs over seeds 0..3000, first at seed 6
(`tern = false`).
-/
namespace WalkInvCheck.Planar
open Spqr

def finishB (seed : Nat) (g : Graph) (tern : Bool) (f : List DfsTree) : List Viol := Id.run do
  let w := g.planarWalk tern f
  let items := w.base.items
  let mut vs : List Viol := []
  let bad (field info : String) : Viol := ⟨seed, tern, field, info⟩
  for k in List.range w.aux.nodePlanarity.size do
    let it := 1 + g.nv + g.ne + k
    match w.aux.nodePlanarity[k]! with
    | .planar m =>
      let ty := items[it]!.type
      unless ty == .S || ty == .P || ty == .R do continue
      let site := s!"item {it} {repr ty} m={m}"
      unless m.size == 4 do vs := vs ++ [bad "pf.size" site]
      unless (List.range 4).all fun s => (items[it]!.ch).contains (1 + g.nv + QE.edge m[s]!) do
        vs := vs ++ [bad "pf.mem" site]
      unless m.toList.Nodup do vs := vs ++ [bad "pf.nodup" site]
      unless items[it]!.ch.all (· < 1 + g.nv + 2 * g.ne) do vs := vs ++ [bad "pf.ch_lt" site]
      let P := nodeSkel g items it
      let qf := flipped g items[it]!.ch w.aux.itemFlips[it]! (capLinked g w.aux.qem m)
      unless P.ves.all fun ve => (List.range 4).all fun z =>
          match qf[4 * ve + z]! with
          | some o => P.ves.contains (QE.edge o)
          | none => false do
        vs := vs ++ [bad "pf.closed" s!"{site} ves={P.ves} qem={P.ves.map fun ve => (List.range 4).map fun z => qf[4 * ve + z]!}"]
      unless IsPlanarEmbedding P.es g.nv (finishRot g w it m) do
        vs := vs ++ [bad "pf.planar" s!"{site} es={P.es} rot={(finishRot g w it m).rotAdj}"]
    | _ => pure ()
  return vs

end WalkInvCheck.Planar
