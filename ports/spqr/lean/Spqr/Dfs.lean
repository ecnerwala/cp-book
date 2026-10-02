import Spqr.Graph

/-!
# Phase 1: lowpoint DFS

A DFS forest is built over the input graph. Every out-edge of a vertex at depth `d` (tree edges to
children, and back edges to proper ancestors or to the vertex itself) is classified by where its
subtree returns to, and the out-edges of each vertex are then stably sorted by that class. This
ordering is what makes the phase-2 walk visit each vertex's return edges from the shallowest
ancestor outwards, exactly as in the C++ implementation.
-/

namespace Spqr

/-- For out-edges returning to a proper ancestor at `lowval < d`: whether the edge is a back edge,
a tree child whose second-lowest return is also a proper ancestor (`type2`), or a tree child all
of whose other returns stay at depth `≥ d` (`type1`). -/
inductive RetKind where
  | type1Child
  | backEdge
  | type2Child
deriving DecidableEq, Repr

/-- Classification of an out-edge of a vertex at depth `d`. -/
inductive OutClass where
  /-- A tree edge whose subtree has no return to depth `≤ d`. -/
  | bridge
  /-- A tree edge whose subtree returns to depth exactly `d` and no higher. -/
  | component
  | selfLoop
  /-- An edge returning to depth `lowval < d`. -/
  | ret (lowval : Nat) (kind : RetKind)
deriving DecidableEq, Repr

namespace RetKind

def rank : RetKind → Nat
  | type1Child => 0
  | backEdge => 1
  | type2Child => 2

end RetKind

namespace OutClass

/-- Sort rank; equal to the C++ `key = 3 * (lowval' + 2) + kind`. -/
def rank : OutClass → Nat
  | bridge => 0
  | component => 3
  | selfLoop => 4
  | ret lowval kind => 3 * (lowval + 2) + kind.rank

def isTree : OutClass → Bool
  | selfLoop => false
  | ret _ .backEdge => false
  | _ => true

def isType1 : OutClass → Bool
  | ret _ .type2Child => false
  | _ => true

/-- The depth the edge returns to, where a bridge from depth `d` "returns" to `d + 1` and a
component / self-loop to `d`. -/
def lowval (d : Nat) : OutClass → Nat
  | bridge => d + 1
  | component => d
  | selfLoop => d
  | ret lowval _ => lowval

end OutClass

mutual
/-- A DFS tree rooted at `v`, with its out-edges in phase-2 order. -/
inductive DfsTree where
  | node (v : Nat) (outs : List DfsOut)
/-- An out-edge `e` from the tree's vertex to `dest`. -/
inductive DfsOut where
  | back (e dest : Nat) (cls : OutClass)
  | tree (e : Nat) (cls : OutClass) (child : DfsTree)
end

namespace DfsTree
def v : DfsTree → Nat | node v _ => v
def outs : DfsTree → List DfsOut | node _ outs => outs
end DfsTree

namespace DfsOut
def e : DfsOut → Nat | back e _ _ => e | tree e _ _ => e
def cls : DfsOut → OutClass | back _ _ cls => cls | tree _ cls _ => cls
def dest : DfsOut → Nat | back _ dest _ => dest | tree _ _ child => child.v
end DfsOut

mutual
/-- Vertices of the tree, in preorder. -/
def DfsTree.verts : DfsTree → List Nat
  | .node v outs => v :: DfsOut.vertsList outs
def DfsOut.vertsList : List DfsOut → List Nat
  | [] => []
  | .back .. :: rest => DfsOut.vertsList rest
  | .tree _ _ child :: rest => child.verts ++ DfsOut.vertsList rest
end

mutual
/-- Edges of the tree (tree and back edges), in preorder. -/
def DfsTree.edges : DfsTree → List Nat
  | .node _ outs => DfsOut.edgesList outs
def DfsOut.edgesList : List DfsOut → List Nat
  | [] => []
  | .back e _ _ :: rest => e :: DfsOut.edgesList rest
  | .tree e _ child :: rest => e :: child.edges ++ DfsOut.edgesList rest
end

/-- The two smallest distinct depths reachable from a subtree, with `(d, d)` meaning "none". -/
abbrev Lowvals := Nat × Nat

def mergeLowvals (l n : Lowvals) : Lowvals :=
  if n.1 < l.1 then (n.1, min n.2 l.1)
  else (l.1, min l.2 (if n.1 == l.1 then n.2 else n.1))

def classify (d : Nat) (isTree : Bool) (n : Lowvals) : OutClass :=
  if n.1 ≥ d then
    if n.1 == d then (if isTree then .component else .selfLoop) else .bridge
  else
    .ret n.1 (if isTree then (if n.2 < d then .type2Child else .type1Child) else .backEdge)

/-- One DFS visit of `v` at depth `d`, entered via edge `prvE`. `fuel` bounds the recursion depth;
it is never exhausted when `fuel ≥ nv - d`. Returns the tree, the subtree's lowvals, and the
updated depth array. -/
def dfsVisit (adj : Array (List (Nat × Nat))) :
    Nat → (v d : Nat) → Option Nat → Array (Option Nat) → DfsTree × Lowvals × Array (Option Nat)
  | 0, v, _, _, depth => (.node v [], (0, 0), depth)
  | fuel + 1, v, d, prvE, depth =>
    let depth := depth.set! v (some d)
    let (outs, lv, depth) := adj[v]!.foldl (init := (([] : List DfsOut), ((d, d) : Lowvals), depth))
      fun (outs, lv, depth) (nxt, e) =>
        if some e == prvE || (depth[nxt]!).any (· > d) then (outs, lv, depth)
        else match depth[nxt]! with
          | none =>
            let (child, n, depth) := dfsVisit adj fuel nxt (d + 1) (some e) depth
            (.tree e (classify d true n) child :: outs, mergeLowvals lv n, depth)
          | some nd =>
            (.back e nxt (classify d false (nd, d)) :: outs, mergeLowvals lv (nd, d), depth)
    (.node v (outs.reverse.mergeSort fun a b => a.cls.rank ≤ b.cls.rank), lv, depth)

/-- The DFS forest of `g`, rooted in `vertOrder` order. -/
def Graph.dfsForest (g : Graph) (vertOrder edgeOrder : List Nat) : List DfsTree :=
  let adj := g.adjacency edgeOrder
  let (roots, _) := (inOrder g.nv vertOrder).foldl (init := (([] : List DfsTree), Array.replicate g.nv (none : Option Nat)))
    fun (roots, depth) rt =>
      if depth[rt]!.isSome then (roots, depth)
      else
        let (t, _, depth) := dfsVisit adj g.nv rt 0 none depth
        (t :: roots, depth)
  roots.reverse

end Spqr
