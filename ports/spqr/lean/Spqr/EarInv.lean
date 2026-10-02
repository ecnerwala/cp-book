import Spqr.WalkSpec

/-!
# The attachment set of an open entry, corrected

`WalkSpec.EntryInv D` lets an open entry `t` be attached only at `t.vStart` and at the open-path
vertices `stackVerts[k]`, `topDepth ≤ k ≤ D` (`TEntry.Term D`). That is **false** mid-walk, for
every `D`, once a type-2 chain is being ascended (`checks/InvCheck.lean`, `cexBuried`): on the cycle
`0-1-2-3-4-5-6-0` with chords `6-1`, `5-2` and the pendant ear `4-7-3`, after `finishEdge` of the
frame `(4, 4)` (type-2 edge `4→5`) the entry `(vStart 6, topDepth 5, [], [Q(5,6)])` is still open
underneath `(5,2)`, `(5,5)` (loop 1 stops at `(5,2)`, whose `topDepth 2 < 4`); its edge `5-6` is
attached at `5`, which is `stackVerts[5]` for `D ≥ 5` but, after the sibling `4→7` has set
`stackVerts[5] := 7`, is `stackVerts[k]` for no `k`. So `WalkInv.ear_lower` and `walkTree_inv'`'s
`Inv d` post-condition are false, and no re-indexing of `D` repairs `Inv`.

What does hold (0 violations on 9000 random multigraphs with up to 12 vertices / 24 edges, at every
`walkTree` start/end and after every `walkOut`; `checks/InvCheck.lean`): the extra attachment
vertices of a buried entry are `vStart`s of entries *above* it — every finished chain vertex `x`
of an open ear still owns an entry `(x, …)` on the stack above everything created inside `x`'s
subtree — and with that clause the depth index `D = d` (the current depth, also for a tree edge's
subtree) suffices. `Term'`/`EntryInv'`/`Inv'` below state exactly that; `Inv' D` is implied by
`Inv D` and is what the walk induction should carry.
-/

namespace Spqr

namespace TEntry

/-- `Term D` plus the bottoms of the entries above `t` on the stack. -/
def Term' (D : Nat) (s : WalkState) (above : List TEntry) (t : TEntry) (v : Nat) : Prop :=
  t.Term D s v ∨ ∃ t' ∈ above, v = t'.vStart

theorem Term'_of_Term {D : Nat} {s : WalkState} {above : List TEntry} {t : TEntry} {v : Nat}
    (h : t.Term D s v) : t.Term' D s above v := .inl h

end TEntry

namespace WalkState

variable {s : WalkState}

/-- `EntryInv D` with the corrected attachment set; `above` are the entries on top of `t`. -/
structure EntryInv' (D : Nat) (s : WalkState) (above : List TEntry) (t : TEntry) : Prop where
  conn : s.g.ConnEdges (t.edges s.g s.items)
  attached : s.g.AttachedIn (t.edges s.g s.items) (t.Term' D s above)

/-- The corrected walk invariant: every open entry is `EntryInv'` relative to the entries above it,
and every closed item is a connected 2-attached piece. -/
structure Inv' (D : Nat) (s : WalkState) : Prop where
  entries : ∀ above t below, s.tstack = above ++ t :: below → s.EntryInv' D above t
  nodes : ∀ i, 1 + s.g.nv + s.g.ne ≤ i → i < s.items.size → s.ItemInv i

theorem EntryInv'.of_entryInv {D : Nat} {above : List TEntry} {t : TEntry} (h : s.EntryInv D t) :
    s.EntryInv' D above t :=
  ⟨h.conn, fun v e e' he he' hE hE' hv hv' => .inl (h.attached v e e' he he' hE hE' hv hv')⟩

theorem Inv'.of_inv {D : Nat} (h : s.Inv D) : s.Inv' D :=
  ⟨fun _ t _ hts => .of_entryInv (h.entries t (by rw [hts]; simp)), h.nodes⟩

theorem Inv'.mono {D D' : Nat} (h : s.Inv' D) (hD : D ≤ D') : s.Inv' D' :=
  ⟨fun above t below hts =>
    ⟨(h.entries above t below hts).conn, fun v e e' he he' hE hE' hv hv' =>
      match (h.entries above t below hts).attached v e e' he he' hE hE' hv hv' with
      | .inl (.inl hv) => .inl (.inl hv)
      | .inl (.inr ⟨k, hk₁, hk₂, hk₃⟩) => .inl (.inr ⟨k, hk₁, Nat.le_trans hk₂ hD, hk₃⟩)
      | .inr h => .inr h⟩,
    h.nodes⟩

end WalkState
end Spqr
