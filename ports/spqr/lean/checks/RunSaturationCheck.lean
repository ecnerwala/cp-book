import Spqr.Walk
import Spqr.SepPair

open Spqr

partial def belowEdges (g : Graph) (items : Array Item) (i : ItemId) : List Nat :=
  (if 1 + g.nv ≤ i ∧ i < 1 + g.nv + g.ne then [i - 1 - g.nv] else []) ++
    items[i]!.ch.flatMap (belowEdges g items)

def boundaryVertices (g : Graph) (es : List Nat) : List Nat :=
  (List.range g.nv).filter fun v =>
    (List.range g.ne).any (fun e => es.contains e && (g.edges[e]!.1 == v || g.edges[e]!.2 == v)) &&
    (List.range g.ne).any (fun e => !es.contains e && (g.edges[e]!.1 == v || g.edges[e]!.2 == v))

#eval do
  let g : Graph := ⟨4, #[(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]⟩
  let forest := g.dfsForest [] []
  let items : Array Item := (g.walk false forest).items
  IO.println s!"postorder={edgePostorderForest forest}"
  for i in List.range items.size do
    if items[i]!.type == NodeType.R then
      let ch : List Nat := items[i]!.ch
      let nc := ch.filter fun (c : Nat) => items[c]!.type != NodeType.V
      let es := nc.flatMap (belowEdges g items)
      IO.println s!"R={i} cap={items[i]!.vs} ch={ch} nonV={nc} edges={es} boundary={boundaryVertices g es}"
      IO.println s!"child caps={nc.map fun (c : Nat) => items[c]!.vs}"
