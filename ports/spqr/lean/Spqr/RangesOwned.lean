import Spqr.RangesCloseTree

/-!
# Ownership of the processed prefix (`walk_rootsCover`)

`Owned` is the invariant threaded through the walk of a vertex `v` (depth `d`, subtree starting at
position `n₀` of the schedule `σ`, root tree starting at `nR`, `orig` stack entries at its start,
`sts[k]` the start of the path vertex at depth `k`, `P` the placement predicate): every processed
edge of the root tree is held by a stack entry or hangs under the vertex item of a path vertex; the
vertex items of strict ancestors hold only edges before the next path vertex's start; the entries
above the `orig` old ones hold only edges of `v`'s subtree, the old ones do not bottom at `v`; every
entry bottoms at a visited vertex and unvisited vertex items are empty.
(Checker: `checkOwned`, kinds `own_*`, at every `finishEdge` pre-state, P site and vertex end.)
-/

namespace Spqr
open WalkM
namespace WalkState

variable {σ : List Nat} {n D : Nat} {s : WalkState}

structure Owned (σ : List Nat) (nR n₀ n v d orig : Nat) (P : ItemId → Prop) (sts : List Nat)
    (s : WalkState) : Prop where
  sv : s.stackVerts[d]! = v
  len : orig ≤ s.tstack.length
  cover : ∀ b, nR ≤ b → b < n → (∃ t ∈ s.tstack, t.edges s.g s.items σ[b]!) ∨
    ∃ k, k ≤ d ∧ Items.EdgeBelow s.g s.items (vertItem s.stackVerts[k]!) σ[b]!
  anc : ∀ k, k < d → ∀ e, e < s.g.ne → Items.EdgeBelow s.g s.items (vertItem s.stackVerts[k]!) e →
    σ.idxOf e < sts[k + 1]!
  vertHi : ∀ e, e < s.g.ne → Items.EdgeBelow s.g s.items (vertItem v) e → σ.idxOf e < n
  vertLo : ∀ e, e < s.g.ne → Items.EdgeBelow s.g s.items (vertItem v) e → n₀ ≤ σ.idxOf e
  new : ∀ t ∈ s.tstack.take (s.tstack.length - orig), ∀ e, e < s.g.ne → t.edges s.g s.items e →
    n₀ ≤ σ.idxOf e
  old : ∀ t ∈ s.tstack.drop (s.tstack.length - orig), t.vStart ≠ v
  vis : ∀ t ∈ s.tstack, P (vertItem t.vStart) ∨ ∃ k, k ≤ d ∧ t.vStart = s.stackVerts[k]!
  fresh : ∀ w, w < s.g.nv → ¬ P (vertItem w) → (∀ k, k ≤ d → w ≠ s.stackVerts[k]!) →
    Items.ch s.items (vertItem w) = []

section Admissions
variable {curV d : Nat} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}

/-- P-site coverage at a `finishEdge` call: with the prefix owned (`Owned`), the merge base `nxt`
(bottoming at `curV`, hence above the old entries) holds only edges of `curV`'s subtree, and the
subtree's processed edges `n₀..n` are all on the stack, so `MergeBaseCover σ n₀ (n + 1)` holds
(checker: `checkOwned`, kinds `own_*`). -/
theorem finishP_ownership (h : CloseBase σ n D curV d o origTstack hasVert s)
    {nR n₀ orig : Nat} {P : ItemId → Prop} {sts : List Nat}
    (ho : Owned σ nR n₀ n curV d orig P sts s) (horig : orig ≤ origTstack)
    (hsts : ∀ k, k ≤ d → sts[k]! ≤ n₀)
    (hhv : o.cls.lowval d < d → o.cls.isType1 = true → hasVert = true) :
    FinishPOwnership σ n curV d o origTstack hasVert s := by
  sorry

end Admissions

end WalkState
end Spqr
