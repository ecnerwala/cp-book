import Spqr.Dfs
import Spqr.Walk
import Spqr.Relabel

namespace Spqr

/-- Build the SPQR tree of `g`. `vertOrder` / `edgeOrder` are (partial) orders in which DFS roots
and adjacency entries are tried; `ternarize` forbids reusing S / P nodes, giving a tree with only
binary/ternary merges. -/
def Graph.spqrTree (g : Graph) (ternarize : Bool := false) (vertOrder edgeOrder : List Nat := []) : SpqrTree :=
  let forest := g.dfsForest vertOrder edgeOrder
  let w := g.walk ternarize forest
  relabelTree g w.items

end Spqr
