namespace Spqr

/-- An undirected multigraph on vertices `0, …, nv-1`; `edges[e] = (u, v)`. Self-loops and
parallel edges are allowed. -/
structure Graph where
  nv : Nat
  edges : Array (Nat × Nat)
deriving Repr

namespace Graph

def ne (g : Graph) : Nat := g.edges.size

end Graph

/-- The indices listed in `order` (in that order), followed by the remaining indices in `[0, n)`
in increasing order. -/
def inOrder (n : Nat) (order : List Nat) : List Nat :=
  order ++ (List.range n).filter (fun i => !order.contains i)

/-- Adjacency lists: `adj[v]` lists `(dest, e)` for every edge `e` incident to `v`, in `edgeOrder`
order. A self-loop appears once in its vertex's list. -/
def Graph.adjacency (g : Graph) (edgeOrder : List Nat) : Array (List (Nat × Nat)) :=
  let adj := (inOrder g.ne edgeOrder).foldl (init := Array.replicate g.nv [])
    fun adj e =>
      let (u, v) := g.edges[e]!
      let adj := adj.modify u ((v, e) :: ·)
      if u != v then adj.modify v ((u, e) :: ·) else adj
  adj.map List.reverse

end Spqr
