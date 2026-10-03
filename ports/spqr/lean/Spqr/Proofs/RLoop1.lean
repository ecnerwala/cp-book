import Spqr.Proofs.RInvFrame
import Spqr.EarLoop1

/-!
# The R branch of Loop 1 under the ear context

`loop1_rBranch` (`RInvFrame.lean`) has no access to the ear layer, so even the positional field
`RBranch.cur_top` (every loop-1 iterate's head tops out at `d`) is out of its reach: after an S
merge the head's `topDepth` is `min` of the merged entries', and only the ear layer's `hi`/`lo`
split (`L1Ctx.hi_top`) says the merged entries top out at `≥ d`. Here the R branch is read through
`L1Ctx`/`L1Inv` (`EarLoop1.lean`, `L1Ctx.ofEar` + `l1_init` at the `finishEdge` site): the loop-1
iterate carries an `L1Piece` on top — the head tops out at `d` (`L1Piece.top`), its bottom is
`l1Bot o done`, its edges are `l1Edges o s done` (the tree edge plus the consumed entries') — which
gives `cur_top`, `cur_ne` and the shape outright (`loop1_r_shape_ctx`). The remaining `RBranch`
fields and `RTop` are the named admission `loop1_rBranch_content_ctx`.
-/

namespace Spqr
open WalkM

namespace WalkState
variable {s : WalkState} {dfs : DfsData}

/-- The `k`-th loop-1 iterate of the tree edge `o` finished at depth `d` from the site state `s`. -/
abbrev l1Iter (d : Nat) (o : DfsOut) (s : WalkState) (k : Nat) : WalkState :=
  iter (loop1Body d s.stackDir[d]!) k (ceS₁ o.dest d o.e (feS₀ d o s))

/-- At an `.R` iterate under the ear context the head is the `L1Piece`: it tops out at `d` and
holds the tree edge. -/
theorem loop1_r_shape_ctx {D d : Nat} {o : DfsOut} {hi lo base : List TEntry} {v₀ : Nat}
    (hc : L1Ctx D d o s hi lo base)
    (h0 : L1Inv D d o s hi lo base (ceS₁ o.dest d o.e (feS₀ d o s)))
    (hv : v₀ < s.g.nv) (he : o.e < s.g.ne) (k : Nat)
    (hk : ∀ j, j ≤ k → result (loop1Cond d) (l1Iter d o s j) = true)
    (hty : l1Ty d s.stackDir[d]! (l1Iter d o s k) = .R) :
    ∃ cur nxt rest, (l1Iter d o s k).tstack = cur :: nxt :: rest ∧
      cur.topDepth = d ∧ nxt.topDepth = d ∧ nxt.vStart ≠ cur.vStart ∧
      (∃ e, e < s.g.ne ∧ cur.edges s.g (l1Iter d o s k).items e) := by
  obtain ⟨cur, nxt, rest, hts, hnt, hne⟩ := loop1_r_shape (hk k le_rfl) hty
  obtain ⟨done, rest', c, -, hP, -, -⟩ := l1_iter hc h0 hv k fun j hj => hk j (Nat.le_of_lt hj)
  have hcur : c = cur := by
    have h := hP.tstack
    rw [hts] at h
    exact (List.cons.inj h).1.symm
  subst hcur
  exact ⟨c, nxt, rest, hts, hP.top, hnt, hne, o.e, he, (hP.edges o.e he).2 (.inl rfl)⟩

/-- **Named admission** (Lemma 4.3 at the R branch, ear-context form). Exact obligation: at an
`.R` iterate of Loop 1 whose shape `cur :: nxt :: rest` with `cur.topDepth = nxt.topDepth = d`,
`nxt.vStart ≠ cur.vStart` and `cur` holding an edge is given (`loop1_r_shape_ctx`), the remaining
`RBranch` fields and `RTop`: `cur_c` (`cur.vStart = stackVerts[d+1]`: the consumed entries all
start at the child at an R iterate — `L1Piece.bot` gives `l1Bot o done`), `cur_piece`/`cur_vs`
(`cur`'s closed items are `(child, stackVerts[d])` pieces), `interior` (the child's subtree is
finished, all its edges are in `rU cur nxt`: postorder interval ownership, `RangesInv`/
`RunSaturation` Facts C/D), `nxt_ne`, `proper`, `nxt_touch_top`/`nxt_touch_bot`, `nxt_no_cu`, and
`RTop` (`cur`/`nxt` are `EntryR` and edge-disjoint: `nxt` is a frontier entry topping out at `d`
not covered by `RInvFront.base`, so it needs the settled-frontier content, cf.
`FinishRShape.settled`; `cur` needs saturation of the child's `(child, stackVerts[d])` classes).
Checked at every R iterate of seeds 0..300 × both modes + 6000 random multigraphs
(`checks/RFinishEdgeCheck.lean`, `RTop`/`RBranch` lines; 0 failures). -/
theorem loop1_rBranch_content_ctx {D d : Nat} {o : DfsOut} {hi lo base : List TEntry}
    (hc : L1Ctx D d o s hi lo base)
    (h0 : L1Inv D d o s hi lo base (ceS₁ o.dest d o.e (feS₀ d o s)))
    (hi₀ : (feS₀ d o s).Inv' D) (hs₀ : Shape (feS₀ d o s))
    (hok : CloseEarsOk D o.dest d o.e s.stackDir[d]! (feS₀ d o s))
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hchild : o.dest = s.stackVerts[d + 1]!) {origTstack : Nat}
    (hR : (feS₀ d o s).RInvFront dfs s.stackVerts[d]! d origTstack) (k : Nat)
    (hk : ∀ j, j ≤ k → result (loop1Cond d) (l1Iter d o s j) = true)
    (hty : l1Ty d s.stackDir[d]! (l1Iter d o s k) = .R)
    (cur nxt : TEntry) (rest : List TEntry)
    (hts : (l1Iter d o s k).tstack = cur :: nxt :: rest)
    (hct : cur.topDepth = d) (hnt : nxt.topDepth = d) (hne : nxt.vStart ≠ cur.vStart)
    (hce : ∃ e, e < s.g.ne ∧ cur.edges s.g (l1Iter d o s k).items e) :
    (l1Iter d o s k).RBranch d cur nxt rest ∧ (l1Iter d o s k).RTop dfs cur nxt := by
  sorry

/-- `loop1_rBranch` under the ear context: shape from `loop1_r_shape_ctx`, content from
`loop1_rBranch_content_ctx`. -/
theorem loop1_rBranch_ctx {D d : Nat} {o : DfsOut} {hi lo base : List TEntry} {v₀ : Nat}
    (hc : L1Ctx D d o s hi lo base)
    (h0 : L1Inv D d o s hi lo base (ceS₁ o.dest d o.e (feS₀ d o s)))
    (hv : v₀ < s.g.nv) (he : o.e < s.g.ne)
    (hi₀ : (feS₀ d o s).Inv' D) (hs₀ : Shape (feS₀ d o s))
    (hok : CloseEarsOk D o.dest d o.e s.stackDir[d]! (feS₀ d o s))
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hchild : o.dest = s.stackVerts[d + 1]!) {origTstack : Nat}
    (hR : (feS₀ d o s).RInvFront dfs s.stackVerts[d]! d origTstack) (k : Nat)
    (hk : ∀ j, j ≤ k → result (loop1Cond d) (l1Iter d o s j) = true)
    (hty : l1Ty d s.stackDir[d]! (l1Iter d o s k) = .R) :
    ∃ cur nxt rest, (l1Iter d o s k).RBranch d cur nxt rest ∧ (l1Iter d o s k).RTop dfs cur nxt := by
  obtain ⟨cur, nxt, rest, hts, hct, hnt, hne, hce⟩ := loop1_r_shape_ctx hc h0 hv he k hk hty
  exact ⟨cur, nxt, rest, loop1_rBranch_content_ctx hc h0 hi₀ hs₀ hok h2 hsp hrt hchild hR k hk hty
    cur nxt rest hts hct hnt hne hce⟩

end WalkState
end Spqr
