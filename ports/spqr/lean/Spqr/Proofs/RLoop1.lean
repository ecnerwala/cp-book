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
gives `cur_top`, `cur_ne` and the shape outright (`loop1_r_shape_ctx`), and `L1Close.mid` for the
next entry gives `interior` once `cur_c` is known (`loop1_r_interior_ctx`). `cur_c`, the
interval/saturation fields (`RBranchFields`) and `RTop` are the named admissions
`loop1_rBranch_cur_c_ctx`, `loop1_rBranch_fields_ctx`, `loop1_rTop_ctx`.
-/

namespace Spqr
open WalkM

namespace WalkState
variable {s : WalkState} {dfs : DfsData}

/-- The `k`-th loop-1 iterate of the tree edge `o` finished at depth `d` from the site state `s`. -/
abbrev rl1Iter (d : Nat) (o : DfsOut) (s : WalkState) (k : Nat) : WalkState :=
  iter (loop1Body d s.stackDir[d]!) k (ceS₁ o.dest d o.e (feS₀ d o s))

/-- At an `.R` iterate under the ear context the head is the `L1Piece`: it tops out at `d` and
holds the tree edge. -/
theorem loop1_r_shape_ctx {D d : Nat} {o : DfsOut} {hi lo base : List TEntry} {v₀ : Nat}
    (hc : L1Ctx D d o s hi lo base)
    (h0 : L1Inv D d o s hi lo base (ceS₁ o.dest d o.e (feS₀ d o s)))
    (hv : v₀ < s.g.nv) (he : o.e < s.g.ne) (k : Nat)
    (hk : ∀ j, j ≤ k → result (loop1Cond d) (rl1Iter d o s j) = true)
    (hty : l1Ty d s.stackDir[d]! (rl1Iter d o s k) = .R) :
    ∃ cur nxt rest, (rl1Iter d o s k).tstack = cur :: nxt :: rest ∧
      cur.topDepth = d ∧ nxt.topDepth = d ∧ nxt.vStart ≠ cur.vStart ∧
      (∃ e, e < s.g.ne ∧ cur.edges s.g (rl1Iter d o s k).items e) := by
  obtain ⟨cur, nxt, rest, hts, hnt, hne⟩ := loop1_r_shape (hk k le_rfl) hty
  obtain ⟨done, rest', c, -, hP, -, -⟩ := l1_iter hc h0 hv k fun j hj => hk j (Nat.le_of_lt hj)
  have hcur : c = cur := by
    have h := hP.tstack
    rw [hts] at h
    exact (List.cons.inj h).1.symm
  subst hcur
  exact ⟨c, nxt, rest, hts, hP.top, hnt, hne, o.e, he, (hP.edges o.e he).2 (.inl rfl)⟩

/-- The `RBranch` fields beyond the shape, `cur_c`, `interior` and `cur_ne`. -/
structure RBranchFields (s : WalkState) (d : Nat) (cur nxt : TEntry) : Prop where
  cur_piece : ∃ i, i ∈ s.entryPieceItems cur
  cur_vs : ∀ i ∈ s.entryPieceItems cur,
    Items.vs s.items i = (some cur.vStart, some s.stackVerts[d]!) ∨
      Items.vs s.items i = (some s.stackVerts[d]!, some cur.vStart)
  nxt_ne : ∃ e, e < s.g.ne ∧ nxt.edges s.g s.items e
  proper : ∃ e, e < s.g.ne ∧ ¬s.rU cur nxt e
  nxt_touch_top : s.g.Touches (nxt.edges s.g s.items) s.stackVerts[d]!
  nxt_touch_bot : s.g.Touches (nxt.edges s.g s.items) nxt.vStart
  nxt_no_cu : ∀ e, nxt.edges s.g s.items e → ¬s.g.Joins e cur.vStart s.stackVerts[d]!

/-- At an `.R` iterate whose head starts at the child, the child is interior to `rU cur nxt`:
`Loop1Spec` gives `L1Close.mid` for `nxt` (the next entry of `hi`, at depth `d`), whose first
disjunct is excluded by `nxt.vStart ≠ cur.vStart` and whose third by the tree edge `o.e ∈ cur`. -/
theorem loop1_r_interior_ctx {D d : Nat} {o : DfsOut} {hi lo base : List TEntry} {v₀ : Nat}
    (hc : L1Ctx D d o s hi lo base)
    (h0 : L1Inv D d o s hi lo base (ceS₁ o.dest d o.e (feS₀ d o s))) (hv : v₀ < s.g.nv)
    (hok : CloseEarsOk D o.dest d o.e s.stackDir[d]! (feS₀ d o s))
    (hchild : o.dest = s.stackVerts[d + 1]!) (k : Nat)
    (hk : ∀ j, j ≤ k → result (loop1Cond d) (rl1Iter d o s j) = true)
    (cur nxt : TEntry) (rest : List TEntry)
    (hts : (rl1Iter d o s k).tstack = cur :: nxt :: rest)
    (hnt : nxt.topDepth = d) (hne : nxt.vStart ≠ cur.vStart)
    (hcc : cur.vStart = (rl1Iter d o s k).stackVerts[d + 1]!) :
    (rl1Iter d o s k).g.Interior ((rl1Iter d o s k).rU cur nxt) cur.vStart := by
  obtain ⟨done, rest', c, hreach, hP, hK, hF⟩ := l1_iter hc h0 hv k fun j hj => hk j (Nat.le_of_lt hj)
  have hcur : cur = c := by
    have h := hP.tstack
    rw [hts] at h
    exact (List.cons.inj h).1
  subst hcur
  obtain ⟨rest'', rfl⟩ : ∃ rest'', rest' = nxt :: rest'' := by
    cases rest' with
    | cons t rest'' =>
      have h := hP.tstack
      rw [hts] at h
      simp only [List.cons_append, List.cons.injEq] at h
      exact ⟨rest'', by rw [h.2.1]⟩
    | nil =>
      exfalso
      obtain ⟨l, lo', rfl⟩ : ∃ l lo', lo = l :: lo' := by
        cases lo with
        | nil => exact absurd rfl hc.lo_ne
        | cons l lo' => exact ⟨l, lo', rfl⟩
      have h := hP.tstack
      rw [hts] at h
      simp only [List.nil_append, List.cons_append, List.cons.injEq] at h
      have := hc.lo_top l rfl
      rw [← h.2.1] at this
      omega
  have hclose := ((hc.spec done nxt rest'' hreach).1 hnt).1
  have hsv := hF.sv
  have hg := hF.g
  have hinc : s.g.Inc o.e s.stackVerts[d + 1]! :=
    hchild ▸ (Graph.inc_of_pairEq (hok.ends : Items.PairEq (o.dest, s.stackVerts[d]!) s.g.edges[o.e]!)).1
  have hcc' : cur.vStart = s.stackVerts[d + 1]! := by rw [hcc, hsv]
  intro e he hince
  rw [hg] at he hince
  rw [hcc'] at hince
  show TEntry.edges _ _ cur e ∨ TEntry.edges _ _ nxt e
  rw [hg]
  rcases hclose.mid with h | h | h
  · exact absurd (h.symm.trans hcc'.symm) hne
  · rcases h e he hince with h' | h'
    · exact .inl ((hP.edges e he).2 h')
    · exact .inr ((hK.edges (by simp) e).2 h')
  · exact absurd ⟨o.e, (hok.e_lt : o.e < s.g.ne), (by unfold l1Edges; exact Or.inl (Or.inl rfl)), hinc⟩ h

/-- **Named admission** (Lemma 4.3 at the R branch, `cur_c`). Exact obligation: at an `.R` iterate
of Loop 1 the head `cur` (the `L1Piece`, bottom `l1Bot o done`) starts at the child
`stackVerts[d+1]`: every entry of `hi` consumed so far either started at the child or was S-merged
into a piece that does (ear-layer fact about the range `hi`; not derivable from `Loop1Spec` alone).
Checked at every R iterate of seeds 0..300 × both modes + 6000 random multigraphs
(`checks/RFinishEdgeCheck.lean`, `RBranch` lines; 0 failures). -/
theorem loop1_rBranch_cur_c_ctx {D d : Nat} {o : DfsOut} {hi lo base : List TEntry} {v₀ : Nat}
    (hc : L1Ctx D d o s hi lo base)
    (h0 : L1Inv D d o s hi lo base (ceS₁ o.dest d o.e (feS₀ d o s))) (hv : v₀ < s.g.nv)
    (hi₀ : (feS₀ d o s).Inv' D) (hs₀ : Shape (feS₀ d o s))
    (hok : CloseEarsOk D o.dest d o.e s.stackDir[d]! (feS₀ d o s))
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hchild : o.dest = s.stackVerts[d + 1]!) {origTstack : Nat}
    (hR : (feS₀ d o s).RInvFront dfs s.stackVerts[d]! d origTstack) (k : Nat)
    (hk : ∀ j, j ≤ k → result (loop1Cond d) (rl1Iter d o s j) = true)
    (hty : l1Ty d s.stackDir[d]! (rl1Iter d o s k) = .R)
    (cur nxt : TEntry) (rest : List TEntry)
    (hts : (rl1Iter d o s k).tstack = cur :: nxt :: rest)
    (hct : cur.topDepth = d) (hnt : nxt.topDepth = d) (hne : nxt.vStart ≠ cur.vStart)
    (hce : ∃ e, e < s.g.ne ∧ cur.edges s.g (rl1Iter d o s k).items e) :
    cur.vStart = (rl1Iter d o s k).stackVerts[d + 1]! := by
  sorry

/-- **Named admission** (Lemma 4.3 at the R branch, interval/saturation fields). Exact obligation:
`RBranchFields` at the iterate — `cur_piece`/`cur_vs` (`cur`'s closed items are
`(child, stackVerts[d])` pieces), `nxt_ne`, `proper` (an edge outside `rU cur nxt`, e.g. the
pending parent edge), `nxt_touch_top`/`nxt_touch_bot`, `nxt_no_cu` (edges joining the child to
`stackVerts[d]` are the child's back edges, P-merged into `cur`): `RangesInv`/`RunSaturation`
Facts C/D at the iterate. Checked as for `loop1_rBranch_cur_c_ctx` (0 failures). -/
theorem loop1_rBranch_fields_ctx {D d : Nat} {o : DfsOut} {hi lo base : List TEntry} {v₀ : Nat}
    (hc : L1Ctx D d o s hi lo base)
    (h0 : L1Inv D d o s hi lo base (ceS₁ o.dest d o.e (feS₀ d o s))) (hv : v₀ < s.g.nv)
    (hi₀ : (feS₀ d o s).Inv' D) (hs₀ : Shape (feS₀ d o s))
    (hok : CloseEarsOk D o.dest d o.e s.stackDir[d]! (feS₀ d o s))
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hchild : o.dest = s.stackVerts[d + 1]!) {origTstack : Nat}
    (hR : (feS₀ d o s).RInvFront dfs s.stackVerts[d]! d origTstack) (k : Nat)
    (hk : ∀ j, j ≤ k → result (loop1Cond d) (rl1Iter d o s j) = true)
    (hty : l1Ty d s.stackDir[d]! (rl1Iter d o s k) = .R)
    (cur nxt : TEntry) (rest : List TEntry)
    (hts : (rl1Iter d o s k).tstack = cur :: nxt :: rest)
    (hct : cur.topDepth = d) (hnt : nxt.topDepth = d) (hne : nxt.vStart ≠ cur.vStart)
    (hce : ∃ e, e < s.g.ne ∧ cur.edges s.g (rl1Iter d o s k).items e) :
    (rl1Iter d o s k).RBranchFields d cur nxt := by
  sorry

/-- **Named admission** (Lemma 4.3 at the R branch, `RTop`). Exact obligation: `cur`/`nxt` are
`EntryR` and edge-disjoint at the iterate — `nxt` is a frontier entry topping out at `d` not
covered by `RInvFront.base`, so it needs the settled-frontier content (cf. `FinishRShape.settled`);
`cur` needs saturation of the child's `(child, stackVerts[d])` classes. Checked as for
`loop1_rBranch_cur_c_ctx` (`RTop` lines; 0 failures). -/
theorem loop1_rTop_ctx {D d : Nat} {o : DfsOut} {hi lo base : List TEntry} {v₀ : Nat}
    (hc : L1Ctx D d o s hi lo base)
    (h0 : L1Inv D d o s hi lo base (ceS₁ o.dest d o.e (feS₀ d o s))) (hv : v₀ < s.g.nv)
    (hi₀ : (feS₀ d o s).Inv' D) (hs₀ : Shape (feS₀ d o s))
    (hok : CloseEarsOk D o.dest d o.e s.stackDir[d]! (feS₀ d o s))
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hchild : o.dest = s.stackVerts[d + 1]!) {origTstack : Nat}
    (hR : (feS₀ d o s).RInvFront dfs s.stackVerts[d]! d origTstack) (k : Nat)
    (hk : ∀ j, j ≤ k → result (loop1Cond d) (rl1Iter d o s j) = true)
    (hty : l1Ty d s.stackDir[d]! (rl1Iter d o s k) = .R)
    (cur nxt : TEntry) (rest : List TEntry)
    (hts : (rl1Iter d o s k).tstack = cur :: nxt :: rest)
    (hct : cur.topDepth = d) (hnt : nxt.topDepth = d) (hne : nxt.vStart ≠ cur.vStart)
    (hce : ∃ e, e < s.g.ne ∧ cur.edges s.g (rl1Iter d o s k).items e) :
    (rl1Iter d o s k).RTop dfs cur nxt := by
  sorry

/-- Lemma 4.3 at the R branch, ear-context form: the shape (`loop1_r_shape_ctx`) plus `cur_c`
(`loop1_rBranch_cur_c_ctx`), `interior` (`loop1_r_interior_ctx`), the interval/saturation fields
(`loop1_rBranch_fields_ctx`) and `RTop` (`loop1_rTop_ctx`). -/
theorem loop1_rBranch_content_ctx {D d : Nat} {o : DfsOut} {hi lo base : List TEntry} {v₀ : Nat}
    (hc : L1Ctx D d o s hi lo base)
    (h0 : L1Inv D d o s hi lo base (ceS₁ o.dest d o.e (feS₀ d o s))) (hv : v₀ < s.g.nv)
    (hi₀ : (feS₀ d o s).Inv' D) (hs₀ : Shape (feS₀ d o s))
    (hok : CloseEarsOk D o.dest d o.e s.stackDir[d]! (feS₀ d o s))
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hchild : o.dest = s.stackVerts[d + 1]!) {origTstack : Nat}
    (hR : (feS₀ d o s).RInvFront dfs s.stackVerts[d]! d origTstack) (k : Nat)
    (hk : ∀ j, j ≤ k → result (loop1Cond d) (rl1Iter d o s j) = true)
    (hty : l1Ty d s.stackDir[d]! (rl1Iter d o s k) = .R)
    (cur nxt : TEntry) (rest : List TEntry)
    (hts : (rl1Iter d o s k).tstack = cur :: nxt :: rest)
    (hct : cur.topDepth = d) (hnt : nxt.topDepth = d) (hne : nxt.vStart ≠ cur.vStart)
    (hce : ∃ e, e < s.g.ne ∧ cur.edges s.g (rl1Iter d o s k).items e) :
    (rl1Iter d o s k).RBranch d cur nxt rest ∧ (rl1Iter d o s k).RTop dfs cur nxt := by
  have hcc := loop1_rBranch_cur_c_ctx hc h0 hv hi₀ hs₀ hok h2 hsp hrt hchild hR k hk hty
    cur nxt rest hts hct hnt hne hce
  have hf := loop1_rBranch_fields_ctx hc h0 hv hi₀ hs₀ hok h2 hsp hrt hchild hR k hk hty
    cur nxt rest hts hct hnt hne hce
  obtain ⟨-, -, -, -, -, -, hF⟩ := l1_iter hc h0 hv k fun j hj => hk j (Nat.le_of_lt hj)
  refine ⟨?_, loop1_rTop_ctx hc h0 hv hi₀ hs₀ hok h2 hsp hrt hchild hR k hk hty cur nxt rest hts
    hct hnt hne hce⟩
  have hce' : ∃ e, e < (rl1Iter d o s k).g.ne ∧ cur.edges (rl1Iter d o s k).g (rl1Iter d o s k).items e := by
    rw [hF.g]; exact hce
  exact
    { tstack := hts
      cur_top := hct
      nxt_top := hnt
      ne := hne
      cur_c := hcc
      cur_piece := hf.cur_piece
      cur_vs := hf.cur_vs
      interior := loop1_r_interior_ctx hc h0 hv hok hchild k hk cur nxt rest hts hnt hne hcc
      cur_ne := hce'
      nxt_ne := hf.nxt_ne
      proper := hf.proper
      nxt_touch_top := hf.nxt_touch_top
      nxt_touch_bot := hf.nxt_touch_bot
      nxt_no_cu := hf.nxt_no_cu }

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
    (hk : ∀ j, j ≤ k → result (loop1Cond d) (rl1Iter d o s j) = true)
    (hty : l1Ty d s.stackDir[d]! (rl1Iter d o s k) = .R) :
    ∃ cur nxt rest, (rl1Iter d o s k).RBranch d cur nxt rest ∧ (rl1Iter d o s k).RTop dfs cur nxt := by
  obtain ⟨cur, nxt, rest, hts, hct, hnt, hne, hce⟩ := loop1_r_shape_ctx hc h0 hv he k hk hty
  exact ⟨cur, nxt, rest, loop1_rBranch_content_ctx hc h0 hv hi₀ hs₀ hok h2 hsp hrt hchild hR k hk hty
    cur nxt rest hts hct hnt hne hce⟩

end WalkState
end Spqr
