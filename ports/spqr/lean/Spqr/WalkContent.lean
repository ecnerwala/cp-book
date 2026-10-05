import Spqr.WalkSpec

/-!
# Open-stack content

`EntryContent D s above t`: what the span items of one open entry `t` are, given the entries `above`
it (checker `content.*` in `checks/WalkInvCheck/Content.lean`). The span items are S/P/R/Q/V roots
(Q leaves); a non-V item has two distinct terminals, each touched by the item's edges and either a
terminal of the entry in the sense of `Term'` (`vStart`, the open path above `topDepth`, the
`vStart`s of the entries above) or one of the entry's own vertex items; a vertex item of the entry
that is not a `Term'` vertex is touched by the entry and has all its edges in the entry or the
entries above it; and every vertex interior to the entry but to no single span item is one of its
vertex items. `ContentInv D s` states it for every entry of the stack.
-/

namespace Spqr
namespace WalkState

structure EntryContent (D : Nat) (s : WalkState) (above : List TEntry) (t : TEntry) : Prop where
  kinds : ∀ j ∈ t.spans.1 ++ t.spans.2,
    Items.type s.items j ∈ [NodeType.S, .P, .R, .Q, .V] ∧
    (Items.type s.items j = .Q → Items.ch s.items j = [])
  two : ∀ j ∈ t.spans.1 ++ t.spans.2, Items.type s.items j ≠ .V →
    ∃ a b, Items.vs s.items j = (some a, some b) ∧ a ≠ b
  root : ∀ j ∈ t.spans.1 ++ t.spans.2, ∀ p, ¬ Items.IsParent s.items p j
  ends : ∀ j ∈ t.spans.1 ++ t.spans.2, Items.type s.items j ≠ .V →
    ∀ a b, Items.vs s.items j = (some a, some b) → ∀ w, w = a ∨ w = b →
      s.g.Touches (Items.EdgeBelow s.g s.items j) w ∧
      (t.Term' D s above w ∨ vertItem w ∈ t.spans.1 ++ t.spans.2)
  vtouch : ∀ w, vertItem w ∈ t.spans.1 ++ t.spans.2 → ¬ t.Term' D s above w →
    s.g.Touches (t.edges s.g s.items) w
  vabove : ∀ w, vertItem w ∈ t.spans.1 ++ t.spans.2 → ¬ t.Term' D s above w →
    s.g.Interior (fun e => t.edges s.g s.items e ∨ ∃ t' ∈ above, t'.edges s.g s.items e) w
  vinner : ∀ w, w < s.g.nv →
    s.g.Touches (t.edges s.g s.items) w → s.g.Interior (t.edges s.g s.items) w →
    (∀ j ∈ t.spans.1 ++ t.spans.2, ¬ s.g.Interior (Items.EdgeBelow s.g s.items j) w) →
    vertItem w ∈ t.spans.1 ++ t.spans.2

def ContentInv (D : Nat) (s : WalkState) : Prop :=
  ∀ above t below, s.tstack = above ++ t :: below → s.EntryContent D above t

end WalkState
end Spqr
