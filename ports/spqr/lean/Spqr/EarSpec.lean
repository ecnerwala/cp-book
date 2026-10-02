import Spqr.Ear
import Spqr.Correctness

/-!
# Ear-level statements for phase 2

The walk proof (`PROOF.md`, §4) is stated on the ear-structured walk `walkEarTree`, which is
transported to `walkTree` by `walkEarTree_eq_walkTree`.
-/

namespace Spqr

/-- `walkEarTree` is the walk, given enough fuel. -/
theorem walkEarTree_eq_walkTree (t : DfsTree) (d fuel : Nat) (h : t.verts.length < fuel) :
    walkEarTree fuel t d = walkTree t d := by
  sorry

theorem walkEar_eq_walk (g : Graph) (tern : Bool) (vo eo : List Nat) :
    g.walkEar tern (g.dfsForest vo eo) = g.walk tern (g.dfsForest vo eo) := by
  sorry

/-- `m` never looks at the tstack it was started on: it computes what it would on an empty tstack,
leaving the original tstack underneath. -/
def TstackLocal (m : WalkM α) (s : WalkState) : Prop :=
  let (a, s') := m.run s
  let (a', s'') := m.run { s with tstack := [] }
  a = a' ∧ s' = { s'' with tstack := s''.tstack ++ s.tstack }

/-- An entry a subtree walk at depth `d` never interacts with: it returns above `d`, is not
started at a vertex of the subtree, and predates every back edge the subtree will see. -/
def TEntry.BelowSubtree (e : TEntry) (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  e.topDepth < d ∧ e.vStart ∉ t.verts ∧ e.firstIdx ≤ s.nxtEdgeIdx

/-- Frame rule: the walk of a subtree is local to the entries it creates. (Hypotheses on `t`
being a well-formed DFS subtree, from `Spqr.Proofs.Dfs`, still to be added.) -/
theorem walkTree_local (t : DfsTree) (d : Nat) (s : WalkState)
    (hS : ∀ e ∈ s.tstack, e.BelowSubtree t d s) : TstackLocal (walkTree t d) s := by
  sorry

/-- PROOF.md Lemma 4.3: a non-first out-edge of `v` (a sub-ear, or a single back edge) nets out to
one entry `(v, lowval)` on top of the stack — a fresh one, or merged (P) into an entry with the
same `(vStart, topDepth)` that was already on top. (DFS well-formedness hypotheses to be added.) -/
theorem earOut_one_entry (fuel v d : Nat) (o : DfsOut) (s : WalkState)
    (hret : o.cls.lowval d < d) (hfuel : fuel > 0) :
    let s' := ((earOut fuel v d o true).run s).2
    ∃ e rest, s'.tstack = e :: rest ∧ e.vStart = v ∧ e.topDepth = o.cls.lowval d ∧
      (rest = s.tstack ∨ ∃ e₀, s.tstack = e₀ :: rest ∧ e₀.vStart = v ∧ e₀.topDepth = e.topDepth) := by
  sorry

/-- PROOF.md Lemma 4.4: once a whole ear has been finished (the chain edge at a frame has been
closed by `finishEdge` and the vertex's remaining out-edges processed), the frame's vertex owns
exactly one entry: the ear as a single two-terminal piece. -/
theorem ascend_frame_one_entry (fuel : Nat) (f : Frame) (fs : List Frame) (s : WalkState)
    (hfuel : fuel > 0) :
    let s₁ := ((finishEdge f.v f.d f.o f.origTstack f.hasVert >>= fun hv =>
      earOuts fuel f.v f.d f.rest hv >>= finishVert f.v f.d).run s).2
    ∃ e, s₁.tstack.length = f.origTstack + 1 ∧ s₁.tstack.head? = some e ∧ e.vStart = f.v := by
  sorry

end Spqr
