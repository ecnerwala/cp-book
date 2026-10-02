import Spqr.Walk

/-!
# The st-order reference

An Even–Tarjan style description of the order in which the walk (`Spqr.Walk`) lists the children
of its S / P / R items. The DFS tree is walked in the same lowval-sorted order; every edge is
handled at its deeper endpoint, exactly where `walkOut` / `finishEdge` see it, and is placed on
one side of the open DFS path:

* an item pushed at depth `l` sits on side `dirs[l]`, where `dirs[l]` is the direction the vertex
  at depth `l` chose for its current out-edge (`walkOut`'s `setStackDir`: `false` for a block
  boundary edge, otherwise the opposite of the direction at the edge's lowval);
* the pieces of one block nest: later pieces lie outside earlier ones (`stNest`);
* when a returning tree edge is finished at a vertex that already has its vertex item
  (`finishEdge`'s vertex close), the whole sub-ear is folded onto side `dirs[lowval]` as one piece.

`refOrder` is the resulting st-order of the leaf items (vertices and edges) of each block; the
walk's child list of an S / P / R item is the restriction of `refOrder` to the item's subtree.
-/

namespace Spqr

/-- One push onto the ear stack, abstracted: `items` (already in st-order) on side `side`. -/
structure StPiece where
  side : Bool
  items : List ItemId

/-- The st-order of a stack of pieces, bottom first: later pieces lie outside earlier ones. -/
def stNestL : List StPiece → List ItemId
  | [] => []
  | p :: ps => stNestL ps ++ (if p.side then [] else p.items)
def stNestR : List StPiece → List ItemId
  | [] => []
  | p :: ps => (if p.side then p.items else []) ++ stNestR ps
def stNest (ps : List StPiece) : List ItemId := stNestL ps ++ stNestR ps

mutual
/-- The pieces pushed while walking the subtree `t` at depth `d`, with `dirs` the directions
chosen along the path above `t`; also the st-orders of the blocks completed inside `t`. -/
def refTree (g : Graph) (t : DfsTree) (d : Nat) (dirs : List Bool) :
    List StPiece × List (List ItemId) :=
  match t with
  | .node v outs =>
    let (ps, blocks, hasVert) := refOuts g v d dirs outs false
    (if hasVert then ps else ps ++ [StPiece.mk true [vertItem v]], blocks)

def refOuts (g : Graph) (v d : Nat) (dirs : List Bool) :
    List DfsOut → Bool → List StPiece × List (List ItemId) × Bool
  | [], hasVert => ([], [], hasVert)
  | o :: rest, hasVert =>
    let (ps, blocks, hasVert) := refOut g v d dirs o hasVert
    let (ps', blocks', hasVert') := refOuts g v d dirs rest hasVert
    (ps ++ ps', blocks ++ blocks', hasVert')

/-- The out-edge `o` of `v`, handled at `v`. -/
def refOut (g : Graph) (v d : Nat) (dirs : List Bool) (o : DfsOut) (hasVert : Bool) :
    List StPiece × List (List ItemId) × Bool :=
  let l := o.cls.lowval d
  if d ≤ l then
    match o with
    | .tree _ _ child =>
      let (ps, blocks) := refTree g child (d + 1) (dirs ++ [false])
      ([], blocks ++ [stNest ps], hasVert)
    | .back .. => ([], [], hasVert)
  else
    let lowDir := dirs.getD l false
    let sd := !lowDir
    let (pre, hasVert) :=
      if !hasVert && o.cls.isType1 then ([StPiece.mk sd [vertItem v]], true) else ([], hasVert)
    let (mid, blocks) :=
      match o with
      | .tree e _ child =>
        let (ps, blocks) := refTree g child (d + 1) (dirs ++ [sd])
        let sub := ps ++ [StPiece.mk sd [edgeItem g e]]
        (if hasVert then [StPiece.mk lowDir (stNest sub)] else sub, blocks)
      | .back e _ _ => ([StPiece.mk lowDir [edgeItem g e]], [])
    let (post, hasVert) := if hasVert then ([], true) else ([StPiece.mk sd [vertItem v]], true)
    (pre ++ mid ++ post, blocks, hasVert)
end

/-- The st-orders of all blocks of the forest, in completion order. -/
def refBlocks (g : Graph) (forest : List DfsTree) : List (List ItemId) :=
  forest.flatMap fun t =>
    let (ps, blocks) := refTree g t 0 []
    blocks ++ [stNest ps]

/-- The reference st-order of the leaf items (blocks concatenated; no item lies in two blocks). -/
def refOrder (g : Graph) (forest : List DfsTree) : List ItemId :=
  (refBlocks g forest).flatMap id

end Spqr
