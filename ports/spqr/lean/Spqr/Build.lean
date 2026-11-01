import Spqr.DfsFast
import Spqr.WalkFast
import Spqr.RelabelFast

namespace Spqr

/-- Build the SPQR tree of `g`. `vertOrder` / `edgeOrder` are (partial) orders in which DFS roots
and adjacency entries are tried; `ternarize` forbids reusing S / P nodes, giving a tree with only
binary/ternary merges. -/
def Graph.spqrTree (g : Graph) (ternarize : Bool := false) (vertOrder edgeOrder : List Nat := []) : SpqrTree :=
  let forest := g.dfsForestFast vertOrder edgeOrder
  let w := g.walkFast ternarize forest
  Fast.relabelTreeFast g (w.items.map Fast.Item.toSlow)

end Spqr
