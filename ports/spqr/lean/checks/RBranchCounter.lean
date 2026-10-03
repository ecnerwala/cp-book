import Spqr.Proofs.RLoop1
import Mathlib.Tactic.IntervalCases

/-!
# Kernel-checked counterexample: the bottom of Loop 1's R piece is not the child

`RBranch.cur_c` claimed `cur.vStart = stackVerts[d + 1]` (the piece bottom is the child) at every
`.R` iterate of Loop 1. False: an S merge earlier in the loop moves the bottom down the path.
Here (gen.py seed 2895, `ternarize = false`) the tree edge `13 : 4 → 0` is finished at depth
`d = 1`; iterate 0 merges the piece `{13}` with the `(0, 2)`-entry of the child's ear into an S
node whose bottom is `5 = stackVerts[3]`, and iterate 1 is `.R` with head `(5, 1)` while
`stackVerts[2] = 0`. The state is the walk's own prefix (`dfsForest_eq`; `site` is the
`walkOutPre` prefix exactly as `RFinishEdgeCounter`). What holds instead is `RBranch.mid`:
`stackVerts[d + 1]` is the bottom or interior to the piece (`child_edges` below: the child's
three edges all sit in `cur`'s S node `25`), the `TwoAttached.of_term` shape.
-/

namespace Spqr.RBranchCounter

open WalkM WalkState

def g : Graph := ⟨7, #[(0, 5), (6, 4), (1, 2), (5, 2), (2, 5), (4, 6), (3, 1), (5, 3), (3, 6),
  (4, 2), (5, 2), (5, 0), (4, 6), (0, 4)]⟩

def t3 : DfsTree := .node 3 [.back 8 6 (.ret 0 .backEdge), .back 7 5 (.ret 3 .backEdge)]
def t1 : DfsTree := .node 1 [.tree 6 (.ret 0 .type2Child) t3]
def t2 : DfsTree := .node 2 [.tree 2 (.ret 0 .type2Child) t1, .back 9 4 (.ret 1 .backEdge),
  .back 4 5 (.ret 3 .backEdge), .back 3 5 (.ret 3 .backEdge)]
def t5 : DfsTree := .node 5 [.tree 10 (.ret 0 .type2Child) t2, .back 0 0 (.ret 2 .backEdge)]
def t0 : DfsTree := .node 0 [.tree 11 (.ret 0 .type2Child) t5]
def o13 : DfsOut := .tree 13 (.ret 0 .type1Child) t0
def t4 : DfsTree := .node 4 [o13, .back 5 6 (.ret 0 .backEdge), .back 1 6 (.ret 0 .backEdge)]
def rootOut : DfsOut := .tree 12 .component t4

theorem dfsForest_eq :
    g.dfsForest [6, 5, 4, 3, 2, 1, 0] [13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0] =
      [.node 6 [rootOut]] := by cbv

/-- The state at vertex `4` (depth `1`) just before the child `0` is walked. -/
def site : WalkState :=
  let s := WalkState.init g false
  let s := { s with stackVerts := s.stackVerts.set! 0 6 }
  let s := after (walkOutPre 6 0 rootOut false) s
  let s := { s with firstOccurrence := s.firstOccurrence.set! 0 g.ne }
  let s := { s with stackVerts := s.stackVerts.set! 1 4 }
  let s := after (walkOutPre 4 1 o13 false) s
  { s with firstOccurrence := s.firstOccurrence.set! 1 g.ne }

/-- The `finishEdge` site of the tree edge `13`. -/
def returned : WalkState := after (walkTree t0 2) site

/-- Loop 1's second iterate. -/
def it : WalkState := rl1Iter 1 o13 returned 1

def cur : TEntry := ⟨5, 1, 5, ([], [25])⟩
def nxt : TEntry := ⟨3, 1, 1, ([], [15, 22, 3, 17, 23, 6])⟩
def rest : List TEntry := [⟨3, 0, 0, ([16], [])⟩, ⟨3, 6, 0, ([], [4])⟩, ⟨4, 1, 0, ([], [5])⟩]

theorem loop1Cond_iter : ∀ j, j ≤ 1 → result (loop1Cond 1) (rl1Iter 1 o13 returned j) = true := by
  intro j hj
  interval_cases j <;> cbv

theorem it_type : l1Ty 1 returned.stackDir[1]! it = .R := by cbv

theorem it_tstack : it.tstack = cur :: nxt :: rest := by cbv

theorem it_stackVerts : it.stackVerts = #[6, 4, 0, 5, 2, 1, 3] := by cbv

theorem not_cur_c : cur.vStart ≠ it.stackVerts[1 + 1]! := by cbv

/-- The former contract (`cur.vStart = stackVerts[d + 1]` at every `.R` iterate) fails here. -/
theorem counter :
    (∀ j, j ≤ 1 → result (loop1Cond 1) (rl1Iter 1 o13 returned j) = true) ∧
      l1Ty 1 returned.stackDir[1]! it = .R ∧ it.tstack = cur :: nxt :: rest ∧
      cur.vStart ≠ it.stackVerts[1 + 1]! :=
  ⟨loop1Cond_iter, it_type, it_tstack, not_cur_c⟩

/-- The child `0` has exactly the edges `0, 11, 13` (`0 → 5`, `5 → 0` and the tree edge), all in
the S node `25` — `stackVerts[d + 1]` is interior to the piece. -/
theorem child_edges :
    (List.range g.ne).filter (fun e => g.edges[e]!.1 == 0 || g.edges[e]!.2 == 0) = [0, 11, 13] := by
  cbv

#print axioms counter
#print axioms dfsForest_eq

end Spqr.RBranchCounter
