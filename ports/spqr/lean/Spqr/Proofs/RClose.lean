import Spqr.RClose
import Spqr.Proofs.RMax
import Spqr.Proofs.SepClasses

/-!
# The R case of loop 1 yields an `RCloseShape` (PROOF.md §4.5, walk side)

From `Inv d` (the depth-indexed walk invariant at the current depth `d`) and the stack shape `RStep`, the structural fields of `RCloseShape` follow for
`U = rU`, `P = ofItems (rPieceItems)`, terminals `nxt.vStart`, `stackVerts[d]`; the content fields
are `RContent`. Hence (`RCloseShape.threeConnected`) the R skeleton is 3-connected.
-/

namespace Spqr

namespace Pieces

variable {g : Graph} {items : Items} {L : List ItemId}

theorem ofItems_k : (ofItems g items L).k = L.length := rfl

theorem ofItems_mem_iff
    (hdisj : ∀ i ∈ L, ∀ j ∈ L, i ≠ j → ∀ e, items.EdgeBelow g i e → ¬items.EdgeBelow g j e)
    (hnodup : L.Nodup) {i e : Nat} :
    (ofItems g items L).Mem i e ↔ e < g.ne ∧ ∃ h : i < L.length, items.EdgeBelow g L[i] e := by
  simp only [Mem, ofItems]
  by_cases he : e < g.ne
  · simp only [he, ↓reduceIte, true_and, List.findIdx?_eq_some_iff_getElem, decide_eq_true_eq]
    constructor
    · rintro ⟨h, hb, -⟩; exact ⟨h, hb⟩
    · rintro ⟨h, hb⟩
      refine ⟨h, hb, fun j hj hbj => ?_⟩
      have hne : L[j] ≠ L[i] := fun heq => by
        have := (hnodup.getElem_inj_iff).1 heq; omega
      exact hdisj _ (List.getElem_mem _) _ (List.getElem_mem _) hne e (by simpa using hbj) hb
  · simp [he]

theorem ofItems_x {i : Nat} (hi : i < L.length) {x y : Nat}
    (hvs : Items.vs items L[i] = (some x, some y)) : (ofItems g items L).x i = x := by
  simp [ofItems, List.getElem?_eq_getElem hi, hvs]

theorem ofItems_y {i : Nat} (hi : i < L.length) {x y : Nat}
    (hvs : Items.vs items L[i] = (some x, some y)) : (ofItems g items L).y i = y := by
  simp [ofItems, List.getElem?_eq_getElem hi, hvs]

end Pieces

namespace WalkState

variable {s : WalkState} {d : Nat} {cur nxt : TEntry} {rest : List TEntry}

theorem rU_of_mem_rItems {i e : Nat} (hi : i ∈ rItems cur nxt)
    (hb : Items.EdgeBelow s.g s.items i e) : s.rU cur nxt e := by
  simp only [rItems, List.mem_append] at hi
  rcases hi with (hi | hi) | (hi | hi)
  · exact .inl ⟨i, List.mem_append_left _ hi, hb⟩
  · exact .inr ⟨i, List.mem_append_left _ hi, hb⟩
  · exact .inr ⟨i, List.mem_append_right _ hi, hb⟩
  · exact .inl ⟨i, List.mem_append_right _ hi, hb⟩

/-- The merged items are well-formed pieces inside `U`. -/
theorem PieceItems.wf (hp : s.PieceItems (s.rPieceItems cur nxt)) (h2 : s.g.TwoConnected)
    (hproper : ∃ e, e < s.g.ne ∧ ¬s.rU cur nxt e) :
    (Pieces.ofItems s.g s.items (s.rPieceItems cur nxt)).WF s.g ∧
      ∀ i e, (Pieces.ofItems s.g s.items (s.rPieceItems cur nxt)).Mem i e → s.rU cur nxt e := by
  set L := s.rPieceItems cur nxt with hL
  have mem : ∀ {i e}, (Pieces.ofItems s.g s.items L).Mem i e ↔
      e < s.g.ne ∧ ∃ h : i < L.length, Items.EdgeBelow s.g s.items L[i] e :=
    Pieces.ofItems_mem_iff hp.disj hp.nodup
  have sub : ∀ i e, (Pieces.ofItems s.g s.items L).Mem i e → s.rU cur nxt e := by
    intro i e h
    obtain ⟨-, hi, hb⟩ := mem.1 h
    exact rU_of_mem_rItems (List.mem_of_mem_filter (List.getElem_mem hi)) hb
  have memE : ∀ {i} (hi : i < L.length), ∀ e, e < s.g.ne →
      ((Pieces.ofItems s.g s.items L).Mem i e ↔ Items.EdgeBelow s.g s.items L[i] e) := by
    intro i hi e he
    rw [mem]; exact ⟨fun h => h.2.2, fun h => ⟨he, hi, h⟩⟩
  have hLmem : ∀ {i} (hi : i < L.length), L[i] ∈ L := fun hi => List.getElem_mem hi
  have proper : ∀ i, i < L.length → ∃ e, e < s.g.ne ∧ ¬(Pieces.ofItems s.g s.items L).Mem i e := by
    intro i hi
    obtain ⟨e, he, hU⟩ := hproper
    exact ⟨e, he, fun h => hU (sub i e h)⟩
  have att : ∀ i, i < L.length → s.g.TwoAttached ((Pieces.ofItems s.g s.items L).Mem i)
      ((Pieces.ofItems s.g s.items L).x i) ((Pieces.ofItems s.g s.items L).y i) := by
    intro i hi
    obtain ⟨x, y, hvs⟩ := hp.vs _ (hLmem hi)
    rw [Pieces.ofItems_x hi hvs, Pieces.ofItems_y hi hvs]
    exact (Graph.TwoAttached.congr (memE hi)).2 (hp.attached _ (hLmem hi) x y hvs)
  have ne : ∀ i, i < L.length → ∃ e, e < s.g.ne ∧ (Pieces.ofItems s.g s.items L).Mem i e := by
    intro i hi
    obtain ⟨e, he, hb⟩ := hp.ne _ (hLmem hi)
    exact ⟨e, he, (memE hi e he).2 hb⟩
  refine ⟨⟨?_, ?_, att, ?_, ?_, proper⟩, sub⟩
  · intro i e h
    obtain ⟨he, hi, -⟩ := mem.1 h
    exact ⟨hi, he⟩
  · intro i hi
    exact (Graph.ConnEdges.congr (memE hi)).2 (hp.conn _ (hLmem hi))
  · intro i hi
    obtain ⟨-, -, hu, hv, -⟩ := Graph.twoAttached_union_classes h2 (att i hi) (ne i hi) (proper i hi)
    exact ⟨hu, hv⟩
  · intro i hi
    obtain ⟨-, hne, -⟩ := Graph.twoAttached_union_classes h2 (att i hi) (ne i hi) (proper i hi)
    exact hne

/-- The R case of loop 1 is an `RCloseShape`: structural fields from `Inv d` and the stack shape,
content fields from `RContent`. -/
theorem RStep.rCloseShape (h : s.Inv d) (h2 : s.g.TwoConnected) (hr : s.RStep d cur nxt rest)
    {dfs : DfsData} (hc : s.RContent dfs d cur nxt) :
    RCloseShape s.g dfs (Pieces.ofItems s.g s.items (s.rPieceItems cur nxt)) (s.rU cur nxt)
      nxt.vStart s.stackVerts[d]! := by
  have hcur := h.entries cur (by simp [hr.tstack])
  have hnxt := h.entries nxt (by simp [hr.tstack])
  obtain ⟨e₀, he₀, hU₀⟩ := hr.proper
  have hcur_out : ∃ e, e < s.g.ne ∧ ¬cur.edges s.g s.items e := ⟨e₀, he₀, fun h => hU₀ (.inl h)⟩
  have hnxt_out : ∃ e, e < s.g.ne ∧ ¬nxt.edges s.g s.items e := ⟨e₀, he₀, fun h => hU₀ (.inr h)⟩
  have hcurA : s.g.TwoAttached (cur.edges s.g s.items) cur.vStart s.stackVerts[d]! := by
    have := TwoAttached.of_term hcur.attached fun k h1 h2 => by have := hr.cur_top; omega
    rwa [hr.cur_top] at this
  have hnxtA : s.g.TwoAttached (nxt.edges s.g s.items) nxt.vStart s.stackVerts[d]! := by
    have := TwoAttached.of_term hnxt.attached fun k h1 h2 => by have := hr.nxt_top; omega
    rwa [hr.nxt_top] at this
  obtain ⟨wf, sub⟩ := hr.pieces.wf h2 ⟨e₀, he₀, hU₀⟩
  have conn : s.g.ConnEdges (s.rU cur nxt) := by
    refine Graph.ConnEdges.union hcur.conn hnxt.conn fun _ _ => ⟨s.stackVerts[d]!, ?_, ?_⟩
    · obtain ⟨e₁, he₁, hE₁⟩ := hr.cur_ne
      obtain ⟨e₂, he₂, hE₂⟩ := hcur_out
      exact (hcurA.touches h2 he₁ hE₁ he₂ hE₂).2
    · obtain ⟨e₁, he₁, hE₁⟩ := hr.nxt_ne
      obtain ⟨e₂, he₂, hE₂⟩ := hnxt_out
      exact (hnxtA.touches h2 he₁ hE₁ he₂ hE₂).2
  have att : s.g.TwoAttached (s.rU cur nxt) nxt.vStart s.stackVerts[d]! := by
    refine Graph.TwoAttached.union hcurA hnxtA ?_
    intro v hv
    simp only [List.mem_cons, List.mem_singleton, List.not_mem_nil, or_false] at hv
    rcases hv with rfl | rfl | rfl | rfl
    · exact .inr (.inr hr.interior)
    · exact .inr (.inl rfl)
    · exact .inl rfl
    · exact .inr (.inl rfl)
  obtain ⟨-, hne, hts, htt, -⟩ :=
    Graph.twoAttached_union_classes h2 att
      (hr.cur_ne.imp fun e h => ⟨h.1, .inl h.2⟩) ⟨e₀, he₀, hU₀⟩
  exact ⟨wf, sub, conn, att, hts, htt, hne, ⟨e₀, he₀, hU₀⟩, hc.single, hc.maximal, hc.type1,
    hc.bond, hc.type2⟩

/-- The R skeleton closed by loop 1 is 3-connected. -/
theorem RStep.threeConnected (h : s.Inv d) (h2 : s.g.TwoConnected) (hr : s.RStep d cur nxt rest)
    {dfs : DfsData} (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g) (hc : s.RContent dfs d cur nxt) :
    (((Pieces.ofItems s.g s.items (s.rPieceItems cur nxt)).addParent s.g (s.rU cur nxt)
      nxt.vStart s.stackVerts[d]!).contract s.g).ThreeConnected :=
  (hr.rCloseShape h h2 hc).threeConnected hsp hrt h2

/-- `nxt` touches `cur.vStart`: otherwise `cur.vStart` is interior to `cur` and `cur` is attached at
`stackVerts[d]` alone, contradicting 2-connectivity. -/
theorem RStep.nxt_touches (h2 : s.g.TwoConnected) (hr : s.RStep d cur nxt rest)
    (hcurA : s.g.TwoAttached (cur.edges s.g s.items) cur.vStart s.stackVerts[d]!) :
    s.g.Touches (nxt.edges s.g s.items) cur.vStart := by
  obtain ⟨e₀, he₀, hU₀⟩ := hr.proper
  obtain ⟨e₁, he₁, hE₁⟩ := hr.cur_ne
  by_contra hno
  have hint : s.g.Interior (cur.edges s.g s.items) cur.vStart := by
    intro e he hv
    rcases hr.interior e he hv with h | h
    · exact h
    · exact absurd ⟨e, he, h, hv⟩ hno
  have hA : s.g.TwoAttached (cur.edges s.g s.items) s.stackVerts[d]! s.stackVerts[d]! := by
    intro v e e' he he' hE hE' hv hv'
    rcases hcurA v e e' he he' hE hE' hv hv' with rfl | h
    · exact absurd (hint e' he' hv') hE'
    · exact .inl h
  exact hA.ne h2 he₁ hE₁ he₀ (fun h => hU₀ (.inl h)) rfl

/-- `Inv d` never holds at an R close with edge-disjoint entries: `Inv d` would attach `nxt` at
`{nxt.vStart, stackVerts[d]}` only, but `cur.vStart` (the child `stackVerts[d+1]`) is touched by
both entries (`RStep.nxt_touches`) and is neither (`RStep.ne`, `TwoAttached.ne`). The invariant that
holds at the R branch of `finishEdge` at depth `d` is `Inv (d+1)`; use `RStep.rCloseShape'`. -/
theorem RStep.inv_absurd (h : s.Inv d) (h2 : s.g.TwoConnected) (hr : s.RStep d cur nxt rest)
    (hdisj : ∀ e, cur.edges s.g s.items e → nxt.edges s.g s.items e → False) : False := by
  have hcur := h.entries cur (by simp [hr.tstack])
  have hnxt := h.entries nxt (by simp [hr.tstack])
  obtain ⟨e₀, he₀, hU₀⟩ := hr.proper
  obtain ⟨e₁, he₁, hE₁⟩ := hr.cur_ne
  have hcurA : s.g.TwoAttached (cur.edges s.g s.items) cur.vStart s.stackVerts[d]! := by
    have := TwoAttached.of_term hcur.attached fun k h1 h2 => by have := hr.cur_top; omega
    rwa [hr.cur_top] at this
  have hnxtA : s.g.TwoAttached (nxt.edges s.g s.items) nxt.vStart s.stackVerts[d]! := by
    have := TwoAttached.of_term hnxt.attached fun k h1 h2 => by have := hr.nxt_top; omega
    rwa [hr.nxt_top] at this
  obtain ⟨e₂, he₂, hE₂, hv₂⟩ := hr.nxt_touches h2 hcurA
  obtain ⟨e₃, he₃, hE₃, hv₃⟩ := (hcurA.touches h2 he₁ hE₁ he₀ (fun h => hU₀ (.inl h))).1
  rcases hnxtA cur.vStart e₂ e₃ he₂ he₃ hE₂ (fun h => hdisj e₃ hE₃ h) hv₂ hv₃ with h | h
  · exact hr.ne h.symm
  · exact hcurA.ne h2 he₁ hE₁ he₀ (fun h => hU₀ (.inl h)) h

/-- The R case of loop 1 is an `RCloseShape`, from the invariant `Inv D` at any depth `D ≥ d` whose
intermediate stack vertices `stackVerts[d+1..D]` all equal `cur.vStart` (at the R branch of
`finishEdge` at depth `d`: `D = d+1` and `cur.vStart = stackVerts[d+1]`, the child). `cur` is then
2-attached at `{cur.vStart, stackVerts[d]}`; `nxt` may attach at `cur.vStart`, which is interior to
`U`, so `U` is 2-attached at `{nxt.vStart, stackVerts[d]}`. -/
theorem RStep.rCloseShape' {D : Nat} (h : s.Inv D)
    (hmid : ∀ k, d < k → k ≤ D → s.stackVerts[k]! = cur.vStart)
    (h2 : s.g.TwoConnected) (hr : s.RStep d cur nxt rest)
    {dfs : DfsData} (hc : s.RContent dfs d cur nxt) :
    RCloseShape s.g dfs (Pieces.ofItems s.g s.items (s.rPieceItems cur nxt)) (s.rU cur nxt)
      nxt.vStart s.stackVerts[d]! := by
  have hcur := h.entries cur (by simp [hr.tstack])
  have hnxt := h.entries nxt (by simp [hr.tstack])
  obtain ⟨e₀, he₀, hU₀⟩ := hr.proper
  obtain ⟨e₁, he₁, hE₁⟩ := hr.cur_ne
  have hcurA : s.g.TwoAttached (cur.edges s.g s.items) cur.vStart s.stackVerts[d]! := by
    have := TwoAttached.of_term hcur.attached fun k h1 h2 =>
      .inl (hmid k (by have := hr.cur_top; omega) h2)
    rwa [hr.cur_top] at this
  have hcurT := (hcurA.touches h2 he₁ hE₁ he₀ (fun h => hU₀ (.inl h))).1
  have hnxtC := hr.nxt_touches h2 hcurA
  obtain ⟨wf, sub⟩ := hr.pieces.wf h2 ⟨e₀, he₀, hU₀⟩
  have conn : s.g.ConnEdges (s.rU cur nxt) :=
    Graph.ConnEdges.union hcur.conn hnxt.conn fun _ _ => ⟨cur.vStart, hcurT, hnxtC⟩
  have att : s.g.TwoAttached (s.rU cur nxt) nxt.vStart s.stackVerts[d]! := by
    intro v e e' he he' hE hE' hv hv'
    have hvc : v ≠ cur.vStart := by
      rintro rfl
      exact hE' (hr.interior e' he' hv')
    rcases hE with hE | hE
    · rcases hcurA v e e' he he' hE (fun h => hE' (.inl h)) hv hv' with h | h
      · exact absurd h hvc
      · exact .inr h
    · rcases hnxt.attached v e e' he he' hE (fun h => hE' (.inr h)) hv hv' with h | ⟨k, h1, h2, rfl⟩
      · exact .inl h
      · rcases Nat.eq_or_lt_of_le h1 with h1 | h1
        · exact .inr (by rw [← h1, hr.nxt_top])
        · exact absurd (hmid k (by have := hr.nxt_top; omega) h2) hvc
  obtain ⟨-, hne, hts, htt, -⟩ :=
    Graph.twoAttached_union_classes h2 att ⟨e₁, he₁, .inl hE₁⟩ ⟨e₀, he₀, hU₀⟩
  exact ⟨wf, sub, conn, att, hts, htt, hne, ⟨e₀, he₀, hU₀⟩, hc.single, hc.maximal, hc.type1,
    hc.bond, hc.type2⟩

/-- `RStep.threeConnected` under the invariant that actually holds at the R branch. -/
theorem RStep.threeConnected' {D : Nat} (h : s.Inv D)
    (hmid : ∀ k, d < k → k ≤ D → s.stackVerts[k]! = cur.vStart)
    (h2 : s.g.TwoConnected) (hr : s.RStep d cur nxt rest)
    {dfs : DfsData} (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g) (hc : s.RContent dfs d cur nxt) :
    (((Pieces.ofItems s.g s.items (s.rPieceItems cur nxt)).addParent s.g (s.rU cur nxt)
      nxt.vStart s.stackVerts[d]!).contract s.g).ThreeConnected :=
  (hr.rCloseShape' h hmid h2 hc).threeConnected hsp hrt h2

end WalkState

end Spqr
