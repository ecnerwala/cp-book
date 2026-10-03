import Spqr.RangesWalk

/-! # The close records at the six close sites of `finishEdge`

`CloseInv` preservation through one `finishEdge`, split into one named statement per block that
creates or completes an item record (PROOF.md §4.6). Each is stated at the state where its block
runs, under `CloseCtx`: what the walk induction knows at a `finishEdge` call (the range
invariant, the ear book-keeping, the frontier, and the scheduled adjacency `FinishR`). The
blocks between them (`feS₀`, `mergeLate`, the `closeVert'` merges/unwrap/retarget, `finishTail`)
only move spans and are covered by the proved `CloseInv` frame lemmas. The empirical checker
(`checks/RangesInvCheck.lean`, `closeSites`) evaluates every `CloseAt` clause at every one of
these block boundaries, under the statement's name. -/

namespace Spqr.WalkState
open WalkM

/-- The context of a `finishEdge curV d o origTstack hasVert` call of the walk at state `s`. -/
structure CloseCtx (σ : List Nat) (n D curV d : Nat) (o : DfsOut) (origTstack : Nat)
    (hasVert : Bool) (s : WalkState) : Prop where
  hD : D = if o.cls.isTree then d + 1 else d
  nodup : σ.Nodup
  lt : ∀ e ∈ σ, e < s.g.ne
  pos : σ[n]? = some o.e
  block : o.block <:+: σ
  ranges : s.RangesInv σ n D
  shape : Shape s
  guards : FinishGuards d o origTstack hasVert s
  book : FinishBook curV d o origTstack hasVert s
  frontier : Frontier (o := o) d origTstack s
  finishR : FinishR σ n curV d o origTstack hasVert s
  close : s.CloseInv
  /-- DFS facts: the out-edge's endpoints; the open path's tree edges, still pending (after `o.e`
  in `σ`); a returning child has an edge of its own. -/
  ends : o.Ends s.g curV
  path : ∀ k, k < d → ∃ e, e < s.g.ne ∧ n < σ.idxOf e ∧
    Items.PairEq (s.stackVerts[k]!, s.stackVerts[k + 1]!) s.g.edges[e]!
  dest_edge : o.cls.isTree = true → o.cls.lowval d < d →
    ∃ e, e < s.g.ne ∧ e ≠ o.e ∧ s.g.Inc e o.dest

/-- The state `finishRest` (hence `finishP`) runs from: after `closeVert'` for a tree edge with
a vertex ear, after `mergeLate` for a tree edge without one, after the push for a back edge. -/
def feRest (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool) (s : WalkState) :
    WalkState :=
  if o.cls.isTree then
    if hasVert then feS₃ curV d o origTstack s else feS₂ d o s
  else feBack curV (o.cls.lowval d) d o s

variable {σ : List Nat} {n D curV d : Nat} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}
  {s : WalkState}

theorem CloseCtx.idxOf (hc : CloseCtx σ n D curV d o origTstack hasVert s) : σ.idxOf o.e = n := by
  obtain ⟨hn, he⟩ := List.getElem?_eq_some_iff.mp hc.pos
  rw [← he]; exact hc.nodup.idxOf_getElem n hn

theorem CloseCtx.ne_of_path (hc : CloseCtx σ n D curV d o origTstack hasVert s) {e : Nat}
    (h : n < σ.idxOf e) : e ≠ o.e := fun heq => by
  rw [heq, hc.idxOf] at h; exact Nat.lt_irrefl _ h

theorem cnt_eq_zero_of_free {i : ItemId} (hs : ∀ t ∈ s.tstack, i ∉ t.spans.1 ++ t.spans.2)
    (hp : ∀ p, ¬ Items.IsParent s.items p i) : s.cnt i = 0 := by
  have h1 : spansCount s.tstack i = 0 := by
    apply List.sum_eq_zero
    intro x hx
    obtain ⟨t, ht, rfl⟩ := List.mem_map.mp hx
    exact List.count_eq_zero.mpr (hs t ht)
  have h2 : chCount s.items i = 0 := Finset.sum_eq_zero fun j _ => List.count_eq_zero.mpr (hp j)
  dsimp [cnt]; omega

theorem pairEq_swap' {a b : Nat} {q : Nat × Nat} (h : Items.PairEq (a, b) q) :
    Items.PairEq (b, a) q := by
  rcases h with h | h
  · exact Or.inr (by rw [← h])
  · obtain ⟨ha, hb⟩ := Prod.mk.inj h; exact Or.inl (by rw [ha, hb])

/-- `feS₀` writes the terminals of the pending edge into its root, loose Q item. -/
theorem feS₀_closeInv (hc : CloseCtx σ n D curV d o origTstack hasVert s) : (feS₀ d o s).CloseInv := by
  obtain ⟨sub, base, -, hE⟩ := hc.book.ear
  exact hc.close.modifyLoose (x := edgeItem s.g o.e) (by show s.g.nv < 1 + s.g.nv + o.e; omega)
    (cnt_eq_zero_of_free hE.q_free hE.q_root) _ (fun _ => rfl)

/-- Pushing the pending edge after `feS₀`: the leaf-Q record (`CloseInv.pushEdge`) from the
written terminals `(stackVerts[d], o.dest)`, their distinctness, and another edge at each. -/
theorem leafQ_closeInv (hc : CloseCtx σ n D curV d o origTstack hasVert s) {v dd : Nat}
    (hne : s.stackVerts[d]! ≠ o.dest) (hp : Items.PairEq (s.stackVerts[d]!, o.dest) s.g.edges[o.e]!)
    (h₁ : ∃ e, e < s.g.ne ∧ e ≠ o.e ∧ s.g.Inc e s.stackVerts[d]!)
    (h₂ : ∃ e, e < s.g.ne ∧ e ≠ o.e ∧ s.g.Inc e o.dest) :
    (after (pushEdgeTstack v dd o.e) (feS₀ d o s)).CloseInv := by
  have hlt : edgeItem s.g o.e < s.items.size := by
    have := hc.shape.size; have := hc.book.e_lt; show 1 + s.g.nv + o.e < _; omega
  have ht : Items.type (feS₀ d o s).items (edgeItem (feS₀ d o s).g o.e) = .Q := by
    show Items.type (s.items.modify _ _) (edgeItem s.g o.e) = .Q
    rw [Items.type_modify (f := fun it =>
      { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) _ _ _
      (fun _ => rfl)]
    exact hc.shape.edge _ hc.book.e_lt
  have hch : Items.ch (feS₀ d o s).items (edgeItem (feS₀ d o s).g o.e) = [] := by
    show Items.ch (s.items.modify _ _) (edgeItem s.g o.e) = []
    rw [Items.ch_modify_of_ch (f := fun it =>
      { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) _ _ _
      (fun _ => rfl)]
    exact hc.book.q
  have hvs : Items.vs (feS₀ d o s).items (edgeItem (feS₀ d o s).g o.e) =
      setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) := by
    show Items.vs (s.items.modify _ _) (edgeItem s.g o.e) = _
    rw [Items.vs_modify_self _ _ _ hlt]
  have hatt : ∀ w, s.g.Inc o.e w → (∃ e, e < s.g.ne ∧ e ≠ o.e ∧ s.g.Inc e w) →
      Items.Att (feS₀ d o s).g (feS₀ d o s).items (edgeItem (feS₀ d o s).g o.e) w := by
    intro w hw ⟨e', he', hne', hinc⟩
    refine ⟨o.e, e', hc.book.e_lt, he', hw, hinc, Relation.ReflTransGen.refl, fun hb => hne' ?_⟩
    have h := Items.below_eq_of_ch_nil hch hb
    have : 1 + s.g.nv + e' = 1 + s.g.nv + o.e := h
    omega
  have h₀ := feS₀_closeInv hc
  obtain ⟨hi₁, hi₂⟩ := Graph.inc_of_pairEq hp
  cases hdir : s.stackDir[d]! with
  | false =>
    exact h₀.pushEdge v dd o.e _ _ ht hch (by rw [hvs, hdir]; rfl) hne hp (hatt _ hi₁ h₁) (hatt _ hi₂ h₂)
  | true =>
    exact h₀.pushEdge v dd o.e _ _ ht hch (by rw [hvs, hdir]; rfl) (Ne.symm hne) (pairEq_swap' hp)
      (hatt _ hi₂ h₂) (hatt _ hi₁ h₁)

/-- The parent edge of `curV = stackVerts[d]` (`d ≥ 1`): pending, hence not `o.e`. -/
theorem CloseCtx.parent_edge (hc : CloseCtx σ n D curV d o origTstack hasVert s) (hd : 0 < d) :
    ∃ e, e < s.g.ne ∧ e ≠ o.e ∧ s.g.Inc e s.stackVerts[d]! := by
  obtain ⟨e, he, hn, hp⟩ := hc.path (d - 1) (by omega)
  rw [show d - 1 + 1 = d by omega] at hp
  exact ⟨e, he, hc.ne_of_path hn, (Graph.inc_of_pairEq hp).2⟩

/-- Leaf Q at a returning tree edge: `pushEdgeTstack o.dest d o.e` after `feS₀` keeps every
record (`CloseInv.pushEdge`: distinct endpoints and both terminal `Att` witnesses of
`edgeItem o.e`, whose other incident edges are not below it). -/
theorem closeEars_closeAt (hc : CloseCtx σ n D curV d o origTstack hasVert s)
    (ht : o.cls.isTree = true) (hlow : o.cls.lowval d < d) :
    (ceS₁ o.dest d o.e (feS₀ d o s)).CloseInv := by
  obtain ⟨lv, kind, ho, hlv⟩ := ret_of_lowval_lt hlow
  obtain ⟨sub, base, -, hE⟩ := hc.book.ear
  have hends := hc.book.ends lv kind ho
  simp only [ht, ↓reduceIte] at hends
  exact leafQ_closeInv hc (hE.path_child ht d (Nat.le_refl d)) (pairEq_swap' hends)
    (hc.parent_edge (by omega)) (hc.dest_edge ht hlow)

/-- One iteration of loop 1: the S/P/R record created (or reopened) by `maybeUnwrapNxt` and
closed by `finishTstackTop` over the merged entry — attachment/terminal facts, the interior
V-child equivalence, two terminals of every non-V child, and the S/P/R shape clause. -/
theorem loop1Body_closeAt (hc : CloseCtx σ n D curV d o origTstack hasVert s)
    (ht : o.cls.isTree = true) (hlow : o.cls.lowval d < d) (k : Nat)
    (hk : ∀ j, j ≤ k → result (loop1Cond d)
      (iter (loop1Body d s.stackDir[d]!) j (ceS₁ o.dest d o.e (feS₀ d o s))) = true)
    (h : (iter (loop1Body d s.stackDir[d]!) k (ceS₁ o.dest d o.e (feS₀ d o s))).CloseInv) :
    (after (loop1Body d s.stackDir[d]!)
      (iter (loop1Body d s.stackDir[d]!) k (ceS₁ o.dest d o.e (feS₀ d o s)))).CloseInv := by
  sorry

/-- The `some item` arm of `closeVertTail` (`vertFinish`): the S or R record of the vertex ear,
from the state `cvS₅` after the merges and the retarget (the `none` arm is the identity). -/
theorem closeVertTail_closeAt (hc : CloseCtx σ n D curV d o origTstack hasVert s)
    (ht : o.cls.isTree = true) (hlow : o.cls.lowval d < d) (hv : hasVert = true)
    (h : (cvS₅ curV s.stackDir[d]! o.cls.isType1 origTstack (feSingle d o s) (feS₂ d o s)).CloseInv) :
    (feS₃ curV d o origTstack s).CloseInv := by
  sorry

/-- `finishP`: when `condP` holds, the new or reused P record over the merged `(curV,
stackVerts[lowval])` entry — at least two virtual edges, all equal to its terminals, no V child. -/
theorem finishP_closeAt (hc : CloseCtx σ n D curV d o origTstack hasVert s)
    (hlow : o.cls.lowval d < d) (h : (feRest curV d o origTstack hasVert s).CloseInv) :
    (after (finishP curV (o.cls.lowval d) o.cls.isType1)
      (feRest curV d o origTstack hasVert s)).CloseInv := by
  sorry

/-- Leaf Q at a returning back edge: `pushEdgeTstack curV lowval o.e` after `feS₀`
(`CloseInv.pushEdge` as for `closeEars_closeAt`; the counter writes are frames). -/
theorem finishBack_closeAt (hc : CloseCtx σ n D curV d o origTstack hasVert s)
    (ht : o.cls.isTree = false) (hlow : o.cls.lowval d < d) :
    (feBack curV (o.cls.lowval d) d o s).CloseInv := by
  obtain ⟨lv, kind, ho, hlv⟩ := ret_of_lowval_lt hlow
  obtain ⟨sub, base, -, hE⟩ := hc.book.ear
  have hends := hc.book.ends lv kind ho
  simp only [ht, Bool.false_eq_true, ↓reduceIte] at hends
  obtain ⟨e, dest, cls, rfl⟩ : ∃ e dest cls, o = .back e dest cls := by
    cases o with
    | back e dest cls => exact ⟨e, dest, cls, rfl⟩
    | tree e cls child => exact absurd (hc.book.tree.mpr ⟨e, cls, child, rfl⟩) (by simp [ht])
  have hends' : Items.PairEq (curV, dest) s.g.edges[e]! := by
    have h := hc.ends; simp only [DfsOut.Ends] at h; exact h
  have hends : Items.PairEq (curV, s.stackVerts[lv]!) s.g.edges[e]! := hends
  have hlvd : s.stackVerts[lv]! ≠ curV := by rw [← hE.sv_d]; exact hE.path lv d hlv (Nat.le_refl d)
  have hdest : dest = s.stackVerts[lv]! := by
    generalize s.g.edges[e]! = q at hends hends'
    obtain ⟨x, y⟩ := q
    simp only [Items.PairEq, Prod.mk.injEq] at hends hends'
    omega
  obtain ⟨e₂, he₂, hn₂, hp₂⟩ := hc.path lv hlv
  have ho' : cls = .ret lv kind := ho
  have hlv' : (DfsOut.back e dest cls).cls.lowval d = lv := by
    show cls.lowval d = lv; rw [ho']; rfl
  rw [hlv']
  have h := leafQ_closeInv hc (v := curV) (dd := lv) (by rw [hE.sv_d]; show curV ≠ dest; rw [hdest]; exact Ne.symm hlvd)
    (by rw [hE.sv_d]; exact hends') (hc.parent_edge (by omega))
    ⟨e₂, he₂, hc.ne_of_path hn₂, by show s.g.Inc e₂ dest; rw [hdest]; exact (Graph.inc_of_pairEq hp₂).1⟩
  exact h.frame rfl rfl fun _ h => h

/-- `finishBoundary`: the nonempty Q record of `o.e` with children `I :: t.spans.2` (bridge),
`backedge.spans.1 ++ t.spans.2` (completed block), or `[O]` (self-loop), then appended under
`vertItem curV` (`CloseInv.vertex_append`). -/
theorem finishBoundary_closeAt (hc : CloseCtx σ n D curV d o origTstack hasVert s)
    (hge : d ≤ o.cls.lowval d) :
    (after (finishEdge curV d o origTstack hasVert) s).CloseInv := by
  sorry

end Spqr.WalkState
