import Spqr.Walk

/-!
# Phase 2, organized by ears

`walkEarTree` is the walk of `Spqr.Walk` with the recursion re-associated around ears. An ear
starts at a vertex `v` with its first returning out-edge and follows *type-2 first children*
downwards (the first-child chain); it ends at the chain vertex whose first returning edge is
type-1 (a back edge or a type-1 child, both of which collapse to a single "virtual back edge").
`descend` walks down a chain, handling the block-boundary edges of every chain vertex on the way
and suspending each chain vertex as a `Frame`; `ascend` then finishes the chain bottom-up: for each
frame, `finishEdge` for the chain edge, the vertex's remaining out-edges (each a sub-ear, or a
trivial single-edge ear), and the vertex's own ear if it still has none.

The tstack contributions of one ear are thus delimited by the frames, which is what the walk
invariant (`PROOF.md`, §4) is stated over. `walkEarTree` computes exactly the same state as
`walkTree` (`walkEarTree_eq_walkTree`).
-/

namespace Spqr
open WalkM

/-- A suspended chain vertex `v` at depth `d`: its first returning out-edge `o` is a type-2 child
whose subtree is being walked; `rest` are the out-edges after `o`. -/
structure Frame where
  v : Nat
  d : Nat
  o : DfsOut
  origTstack : Nat
  hasVert : Bool
  rest : List DfsOut

/-- Whether the ear continues through `o`: a tree edge to a type-2 child. -/
def DfsOut.continuesChain (d : Nat) : DfsOut → Bool
  | .tree _ cls _ => cls.lowval d < d && !cls.isType1
  | .back .. => false

/-- `walkOut` with the subtree walk abstracted. -/
def walkOutWith (walkChild : DfsTree → Nat → WalkM Unit) (v d : Nat) (o : DfsOut) (hasVert : Bool) :
    WalkM Bool := do
  let lowval := o.cls.lowval d
  let lowDir ← stackDir lowval
  setStackDir d (if lowval ≥ d then false else !lowDir)
  let hasVert ← do
    if !hasVert && lowval < d && o.cls.isType1 then
      pushVertTstack v d
      pure true
    else pure hasVert
  let origTstack ← tstackSize
  match o with
  | .tree _ _ child =>
    modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }
    walkChild child (d + 1)
  | .back .. => pure ()
  finishEdge v d o origTstack hasVert

def walkOutsWith (walkChild : DfsTree → Nat → WalkM Unit) (v d : Nat) :
    List DfsOut → Bool → WalkM Bool
  | [], hasVert => return hasVert
  | o :: rest, hasVert => do
    walkOutsWith walkChild v d rest (← walkOutWith walkChild v d o hasVert)

/-- A vertex with no returning edge of its own becomes its own ear. -/
def finishVert (v d : Nat) (hasVert : Bool) : WalkM Unit := do
  unless hasVert do
    setStackDir d true
    pushVertTstack v d

mutual
/-- Ear-structured walk of the subtree `t` at depth `d`; `fuel ≥` the number of vertices of `t`. -/
def walkEarTree (fuel : Nat) (t : DfsTree) (d : Nat) : WalkM Unit := do
  ascend fuel (← descend fuel t d [])
termination_by (fuel, 2, 0)

/-- Walk down the first-child chain starting at `t`, suspending the chain vertices in `acc`
(deepest first). The chain's last vertex is completed here. -/
def descend : Nat → DfsTree → Nat → List Frame → WalkM (List Frame)
  | 0, _, _, acc => pure acc
  | fuel + 1, .node v outs, d, acc => do
    modify fun s => { s with stackVerts := s.stackVerts.set! d v }
    let (boundary, rets) := outs.span fun o => o.cls.lowval d ≥ d
    let hasVert ← earOuts (fuel + 1) v d boundary false
    match rets with
    | .tree e cls child :: rest =>
      if cls.lowval d < d && !cls.isType1 then
        setStackDir d (!(← stackDir (cls.lowval d)))
        let origTstack ← tstackSize
        modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }
        descend fuel child (d + 1) (⟨v, d, .tree e cls child, origTstack, hasVert, rest⟩ :: acc)
      else
        finishVert v d (← earOuts (fuel + 1) v d (.tree e cls child :: rest) hasVert)
        pure acc
    | rets =>
      finishVert v d (← earOuts (fuel + 1) v d rets hasVert)
      pure acc
termination_by fuel _ _ _ => (fuel, 1, 0)

/-- Finish the suspended chain vertices, deepest first. -/
def ascend : Nat → List Frame → WalkM Unit
  | 0, _ | _, [] => pure ()
  | fuel + 1, f :: fs => do
    let hasVert ← finishEdge f.v f.d f.o f.origTstack f.hasVert
    finishVert f.v f.d (← earOuts (fuel + 1) f.v f.d f.rest hasVert)
    ascend (fuel + 1) fs
termination_by fuel frames => (fuel, 1, frames.length)

/-- `walkOutWith (walkEarTree fuel)`. -/
def earOut : Nat → Nat → Nat → DfsOut → Bool → WalkM Bool
  | 0, _, _, _, hasVert => pure hasVert
  | fuel + 1, v, d, o, hasVert => do
    let lowval := o.cls.lowval d
    let lowDir ← stackDir lowval
    setStackDir d (if lowval ≥ d then false else !lowDir)
    let hasVert ← do
      if !hasVert && lowval < d && o.cls.isType1 then
        pushVertTstack v d
        pure true
      else pure hasVert
    let origTstack ← tstackSize
    match o with
    | .tree _ _ child =>
      modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }
      walkEarTree fuel child (d + 1)
    | .back .. => pure ()
    finishEdge v d o origTstack hasVert
termination_by fuel _ _ _ _ => (fuel, 0, 0)

def earOuts (fuel : Nat) (v d : Nat) : List DfsOut → Bool → WalkM Bool
  | [], hasVert => return hasVert
  | o :: rest, hasVert => do earOuts fuel v d rest (← earOut fuel v d o hasVert)
termination_by outs _ => (fuel, 0, outs.length + 1)
end

def walkEarForest (fuel : Nat) (forest : List DfsTree) : WalkM Unit :=
  forest.forM fun t => do
    walkEarTree fuel t 0
    let top ← popTstack
    modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 }

/-- Phase 2 entry point, ear-structured. -/
def Graph.walkEar (g : Graph) (ternarize : Bool) (forest : List DfsTree) : WalkState :=
  (walkEarForest g.nv forest).run (WalkState.init g ternarize) |>.2

end Spqr
