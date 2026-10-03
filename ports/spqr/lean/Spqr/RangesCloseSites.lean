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

namespace Spqr.Items

theorem vs_modify_ch {items : Items} (j : ItemId) (cs : List ItemId) (p : ItemId) :
    Items.vs (items.modify j fun it => { it with ch := cs }) p = items.vs p :=
  vs_modify_of_vs items j p (fun it => { it with ch := cs }) fun _ => rfl

theorem type_modify_ch {items : Items} (j : ItemId) (cs : List ItemId) (p : ItemId) :
    Items.type (items.modify j fun it => { it with ch := cs }) p = items.type p :=
  type_modify items j p (fun it => { it with ch := cs }) fun _ => rfl

theorem type_modify_vs {items : Items} (j : ItemId) (vsv : Option Nat × Option Nat) (p : ItemId) :
    Items.type (items.modify j fun it => { it with vs := vsv }) p = items.type p :=
  type_modify items j p (fun it => { it with vs := vsv }) fun _ => rfl

theorem ch_modify_vs {items : Items} (j : ItemId) (vsv : Option Nat × Option Nat) (p : ItemId) :
    Items.ch (items.modify j fun it => { it with vs := vsv }) p = items.ch p :=
  ch_modify_ch_eq j (fun it => { it with vs := vsv }) (fun _ => rfl) p

theorem Below_modify_vs {items : Items} (j : ItemId) (vsv : Option Nat × Option Nat) (a i : ItemId) :
    Items.Below (items.modify j fun it => { it with vs := vsv }) a i ↔ items.Below a i :=
  Below_modify_ch_eq j (fun it => { it with vs := vsv }) fun _ => rfl

/-- The record of a vertex item: no I/O children, every child a nonempty Q. -/
theorem CloseAt.vertex' {g : Graph} {items : Items} {v : Nat} (hv : v < g.nv)
    (ht : items.type (vertItem v) = .V)
    (hio : ∀ c, items.IsParent (vertItem v) c → items.type c ≠ .I ∧ items.type c ≠ .O)
    (hne : ∀ c, items.IsParent (vertItem v) c → items.ch c ≠ []) :
    CloseAt g items (vertItem v) := by
  constructor
  · simp [ht]
  · simp [ht]
  · simp [ht]
  · simp [ht]
  · simp [ht]
  · intro c h hio'; exact (hio'.elim (hio c h).1 (hio c h).2).elim
  · intro e he; exact absurd he (by show 1 + v ≠ 1 + g.nv + e; omega)
  · intro e he; exact absurd he (by show 1 + v ≠ 1 + g.nv + e; omega)
  · intro _ _ _ c h; exact hne c h
  · simp [ht]
  · simp [ht]
  · simp [ht]

/-- The record of a nonempty (root) Q item attached only at `u`. -/
theorem CloseAt.rootQ {g : Graph} {items : Items} {e u c : Nat} (he : e < g.ne)
    (ht : items.type (edgeItem g e) = .Q) (hne : items.ch (edgeItem g e) ≠ [])
    (hvs : items.vs (edgeItem g e) = (some u, none)) (hu : g.Inc e u)
    (hatt : ∀ v, items.Att g (edgeItem g e) v → v = u)
    (hct : items.type c ∉ [NodeType.F, .V]) (hcq : items.type c = .Q → items.ch c = [])
    (hloop : (g.edges[e]!).1 = (g.edges[e]!).2 →
      items.ch (edgeItem g e) = [c] ∧ items.vs c = (some u, none))
    (hnon : (g.edges[e]!).1 ≠ (g.edges[e]!).2 → ∃ w, w < g.nv ∧ PairEq (u, w) g.edges[e]! ∧
      items.ch (edgeItem g e) = [c, vertItem w] ∧
      ∃ a b, items.vs c = (some a, some b) ∧ PairEq (a, b) (u, w)) :
    CloseAt g items (edgeItem g e) := by
  constructor
  · intro _ v hv; left; rw [hvs, hatt v hv]
  · rintro (h | ⟨_, h⟩)
    · simp [ht] at h
    · exact absurd h hne
  · simp [ht]
  · simp [ht]
  · intro _ _ h; simp [ht] at h
  · intro _ _ _; exact ht
  · intro _ _ _ h; exact absurd h hne
  · intro e' he' _ _
    obtain rfl : e = e' := by have : 1 + g.nv + e = 1 + g.nv + e' := he'; omega
    exact ⟨u, c, hvs, hu, hct, hcq, hloop, hnon⟩
  · intro v hv hvn; exact absurd hv (by show 1 + g.nv + e ≠ 1 + v; omega)
  · simp [ht]
  · simp [ht]
  · simp [ht]

end Spqr.Items

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
  dest_lt : o.dest < s.g.nv
  /-- Boundary (`d ≤ lowval`) facts the ear layer does not export: a self-loop returns to `curV`;
  the popped block entry holds exactly the child's vertex item; at a component edge the
  `(o.dest, lowval)` entry above it holds one closed node, or leaf Q, on the block's terminals
  `curV, o.dest`. -/
  bd_loop : o.cls.isTree = false → d ≤ o.cls.lowval d → o.dest = curV
  bd_vert : o.cls.isTree = true → d ≤ o.cls.lowval d →
    ∀ t ∈ (if o.cls.lowval d = d + 1 then s.tstack.head? else s.tstack.tail.head?),
      t.spans.2 = [vertItem o.dest]
  bd_node : o.cls.isTree = true → d ≤ o.cls.lowval d → o.cls.lowval d ≠ d + 1 →
    ∀ b ∈ s.tstack.head?, ∃ c, b.spans.1 = [c] ∧ Items.type s.items c ∉ [NodeType.F, .V] ∧
      (Items.type s.items c = .Q → Items.ch s.items c = []) ∧
      ∃ a b', Items.vs s.items c = (some a, some b') ∧ Items.PairEq (a, b') (curV, o.dest)

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

theorem cnt_lt_size (hs : Shape s) {i : ItemId} (hc : 0 < s.cnt i) : i < s.items.size := by
  by_contra hn
  have h1 : spansCount s.tstack i = 0 := by
    apply List.sum_eq_zero
    intro x hx
    obtain ⟨t, ht, rfl⟩ := List.mem_map.mp hx
    exact List.count_eq_zero.mpr fun hm => hn (hs.span t ht i hm)
  have h2 : chCount s.items i = 0 :=
    Finset.sum_eq_zero fun j _ => List.count_eq_zero.mpr fun hm => hn (hs.ch_lt j i hm)
  have : s.cnt i = 0 := by dsimp [cnt]; omega
  omega

theorem CloseInv.alloc' (h : s.CloseInv) (hs : Shape s) (ty : NodeType) :
    ({ s with items := s.items.push ⟨ty, (none, none), []⟩ } : WalkState).CloseInv := by
  constructor
  intro i _ hc
  have hc' : i ≤ s.g.nv ∨ 0 < s.cnt i := by
    simpa only [cnt, Items.chCount_push s.items ⟨ty, (none, none), []⟩ rfl] using hc
  have hi : i < s.items.size := by
    rcases hc' with hc' | hc'
    · have := hs.size; omega
    · exact cnt_lt_size hs hc'
  apply (h.closed i hi hc').frame
  refine ⟨Items.type_push_of_ne _ (Nat.ne_of_lt hi), Items.vs_push_of_ne _ (Nat.ne_of_lt hi),
    Items.ch_push_nil _ rfl _, ?_, ?_, ?_, ?_⟩
  · intro c hc; exact Items.type_push_of_ne _ (Nat.ne_of_lt (hs.ch_lt i c hc))
  · intro c hc; exact Items.vs_push_of_ne _ (Nat.ne_of_lt (hs.ch_lt i c hc))
  · intro c _; exact Items.ch_push_nil _ rfl _
  · intro j _ e; exact Items.Below_push_nil _ rfl

/-- `CloseInv.vertex_append` from `Shape` alone: the vertex record only needs its children to be
nonempty Qs, and every existing child is a counted item. -/
theorem CloseInv.vertex_append' (h : s.CloseInv) (hs : Shape s) (v e : Nat) (hv : v < s.g.nv)
    (he : e < s.g.ne) (hz : s.cnt (vertItem v) = 0)
    (hq : Items.CloseAt s.g s.items (edgeItem s.g e))
    (hqne : Items.ch s.items (edgeItem s.g e) ≠ []) :
    ({ s with items := s.items.modify (vertItem v) fun it =>
      { it with ch := it.ch ++ [edgeItem s.g e] } } : WalkState).CloseInv := by
  have hvlt : vertItem v < s.items.size := by have := hs.size; show 1 + v < _; omega
  have hn : Items.NoParent s.items (vertItem v) := noParent_of_cnt_eq_zero hz
  have hvert := h.closed (vertItem v) hvlt (Or.inl (by show 1 + v ≤ s.g.nv; omega))
  have htv := hs.vert v hv
  have heq : (s.items.modify (vertItem v) fun it => { it with ch := it.ch ++ [edgeItem s.g e] }) =
      s.items.modify (vertItem v) fun it =>
        { it with ch := Items.ch s.items (vertItem v) ++ [edgeItem s.g e] } := by
    apply Array.ext
    · simp
    · intro i hi₁ hi₂
      simp only [Array.getElem_modify]
      split
      · next hEq => subst i; simp [Items.ch, Array.getElem?_eq_getElem hvlt]
      · rfl
  rw [heq]
  have hch : ∀ c, Items.IsParent (s.items.modify (vertItem v) fun it =>
      { it with ch := Items.ch s.items (vertItem v) ++ [edgeItem s.g e] }) (vertItem v) c ↔
      c ∈ Items.ch s.items (vertItem v) ∨ c = edgeItem s.g e := by
    intro c
    simp only [Items.IsParent, Items.ch_modify_self _ _ _ hvlt, List.mem_append, List.mem_singleton]
  have hcv : ∀ c, c ∈ Items.ch s.items (vertItem v) ∨ c = edgeItem s.g e → vertItem v ≠ c := by
    rintro c (hc | rfl) heq
    · subst heq; exact hn _ hc
    · have : 1 + v = 1 + s.g.nv + e := heq; omega
  apply h.writeChildren hvlt hn
  · apply Items.CloseAt.vertex' hv
    · rw [Items.type_modify_ch]; exact htv
    · intro c hc
      have hc' := (hch c).1 hc
      rw [Items.type_modify_ch]
      rcases hc' with hc' | rfl
      · constructor <;> intro hio <;>
          · have := hvert.io_parent c hc' (by simp [hio]); rw [htv] at this; cases this
      · rw [hs.edge e he]; simp
    · intro c hc
      have hc' := (hch c).1 hc
      rw [Items.ch_modify_ne _ _ _ _ (hcv c hc')]
      rcases hc' with hc' | rfl
      · exact hvert.q_under_v v rfl hv c hc'
      · exact hqne
  · intro c hc
    rcases List.mem_append.mp hc with hc | hc
    · have hpos : 0 < s.cnt c := by
        have := Items.count_le_chCount s.items hvlt c
        have := List.count_pos_iff.mpr hc
        dsimp [cnt]; omega
      exact h.closed c (cnt_lt_size hs hpos) (Or.inr hpos)
    · obtain rfl := List.mem_singleton.mp hc; exact hq

theorem eq_of_inc_pairEq {g : Graph} {e a b v : Nat} (hp : Items.PairEq (a, b) g.edges[e]!)
    (hv : g.Inc e v) : v = a ∨ v = b := by
  have hv' : (g.edges[e]!).1 = v ∨ (g.edges[e]!).2 = v := hv
  generalize g.edges[e]! = q at hp hv'
  obtain ⟨x, y⟩ := q
  simp only [Items.PairEq, Prod.mk.injEq] at hp hv'
  omega

theorem CloseCtx.ends_tree (hc : CloseCtx σ n D curV d o origTstack hasVert s)
    (hT : o.cls.isTree = true) : Items.PairEq (o.dest, curV) s.g.edges[o.e]! := by
  obtain ⟨e, cls, child, rfl⟩ := hc.book.tree.mp hT
  have h := hc.ends
  unfold DfsOut.Ends at h
  exact h.1

theorem CloseCtx.ends_back (hc : CloseCtx σ n D curV d o origTstack hasVert s)
    (hT : o.cls.isTree = false) : Items.PairEq (curV, o.dest) s.g.edges[o.e]! := by
  cases o with
  | back e dest cls => have h := hc.ends; unfold DfsOut.Ends at h; exact h
  | tree e cls child => exact absurd (hc.book.tree.mpr ⟨e, cls, child, rfl⟩) (by simp [hT])

/-- `finishBoundary`: the nonempty Q record of `o.e` with children `I :: t.spans.2` (bridge),
`backedge.spans.1 ++ t.spans.2` (completed block), or `[O]` (self-loop), then appended under
`vertItem curV` (`CloseInv.vertex_append'`). -/
theorem finishBoundary_closeAt (hc : CloseCtx σ n D curV d o origTstack hasVert s)
    (hge : d ≤ o.cls.lowval d) :
    (after (finishEdge curV d o origTstack hasVert) s).CloseInv := by
  obtain ⟨sub, base, hlen, hE⟩ := hc.book.ear
  have hok := ear_boundary hge hc.guards hc.ranges.inv hc.shape hc.book hc.hD
  have hs := hc.shape
  have hi := hc.ranges.inv
  have hge' : o.cls.lowval d ≥ d := hge
  have helt := hc.book.e_lt
  have hvlt := hc.book.v_lt
  have hqlt : edgeItem s.g o.e < 1 + s.g.nv + s.g.ne := by show 1 + s.g.nv + o.e < _; omega
  have hqsz : edgeItem s.g o.e < s.items.size := Nat.lt_of_lt_of_le hqlt hs.size
  have hvsz : vertItem curV < s.items.size := by show 1 + curV < _; have := hs.size; omega
  set s₀ := { s with items := s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) } }
    with hs₀
  have b₀ : BStep D s s₀ := BStep.modifyVs hi hs _ _ hqlt
  have h₀ : s₀.CloseInv := hc.close.modifyLoose (x := edgeItem s.g o.e)
    (by show s.g.nv < 1 + s.g.nv + o.e; omega) (cnt_eq_zero_of_free hE.q_free hE.q_root) _
    (fun _ => rfl)
  set s₁ := { s₀ with totBlocks := s₀.totBlocks + 1 } with hs₁
  have b₁ : BStep D s s₁ := b₀.trans (BStep.frame b₀.inv b₀.shape)
  have h₁ : s₁.CloseInv := h₀.frame rfl rfl fun _ h => h
  have hsz₁ : s₁.items.size = s.items.size := by simp [hs₁, hs₀]
  have hq_root₁ : ∀ p, ¬ Items.IsParent s₁.items p (edgeItem s.g o.e) := fun p h =>
    hok.q_root p ((isParent_modifyVs_iff _ _ _ _ _).1 h)
  have hv_root₁ : ∀ p, ¬ Items.IsParent s₁.items p (vertItem curV) := fun p h =>
    hok.v_root p ((isParent_modifyVs_iff _ _ _ _ _).1 h)
  have hB₁ : ∀ a i, Items.Below s₁.items a i ↔ Items.Below s.items a i := fun a i =>
    Items.Below_modify_vs _ _ _ _
  have hty₁ : ∀ j, Items.type s₁.items j = Items.type s.items j := fun j =>
    Items.type_modify_vs _ _ _
  have hch₁ : ∀ j, Items.ch s₁.items j = Items.ch s.items j := fun j =>
    Items.ch_modify_vs _ _ _
  have hvs₁ : Items.vs s₁.items (edgeItem s.g o.e) = (some curV, none) := by
    show Items.vs (s.items.modify _ _) _ = _
    rw [Items.vs_modify_self _ _ _ hqsz]
  have hvs₁' : ∀ j, j ≠ edgeItem s.g o.e → Items.vs s₁.items j = Items.vs s.items j := fun j hj =>
    Items.vs_modify_ne _ _ _ _ (Ne.symm hj)
  have hsv : s.stackVerts[d]! = curV := hE.sv_d
  show wp (finishEdge curV d o origTstack hasVert) (fun _ s' => s'.CloseInv) s
  rw [finishEdge_eq]
  simp only [finishEdge', wp_bind, wp_get, wp_stackDir, hge', ↓reduceIte]
  unfold finishBoundary
  simp only [wp_bind, wp_modifyItem, wp_modify, wp_ite, wp_allocItem, wp_makeVs, wp_popTstack,
    wp_pure]
  split
  · rename_i hT
    have hends := hc.ends_tree hT
    have hdne : o.dest ≠ curV := fun h =>
      hE.path_child hT d (Nat.le_refl d) (by rw [hsv, h])
    have hedge : (s.g.edges[o.e]!).1 ≠ (s.g.edges[o.e]!).2 := by
      intro h
      generalize s.g.edges[o.e]! = q at hends h
      obtain ⟨x, y⟩ := q
      simp only [Items.PairEq, Prod.mk.injEq] at hends h
      omega
    split
    · -- bridge
      rename_i hL
      have hL' : o.cls.lowval d = d + 1 := by simpa using hL
      obtain ⟨t, hsub, hvs', hdep⟩ := hE.bd_bridge hT hL'
      have hts : s.tstack = t :: base := by rw [hE.tstack, hsub]; rfl
      have ht : t ∈ s.tstack := by rw [hts]; exact List.mem_cons_self ..
      have hside : t.spans.1 = [] := by
        have := hE.bd_side hT hge; rw [if_pos hL', hts] at this; exact this t rfl
      have hvert : t.spans.2 = [vertItem o.dest] := by
        have := hc.bd_vert hT hge; rw [if_pos hL', hts] at this; exact this t rfl
      have hhead : s.tstack.head! = t := by rw [hts]; rfl
      set vsI := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) with hvsI
      set s₂ := { s₁ with items := s₁.items.push ⟨.I, (none, none), []⟩ } with hs₂
      have b₂ : BStep D s s₂ := b₁.trans (BStep.alloc b₁.inv b₁.shape .I)
      have h₂ : s₂.CloseInv := h₁.alloc' b₁.shape .I
      set s₃ := { s₂ with items := s₂.items.modify s₁.items.size fun it => { it with vs := vsI } }
        with hs₃
      have b₃ : BStep D s s₃ := b₂.trans (BStep.modifyVs_leaf b₂.inv b₂.shape s₁.items.size vsI
        (Items.ch_push_size _ rfl))
      have hszne : s₁.items.size ≠ edgeItem s.g o.e := by rw [hsz₁]; exact Nat.ne_of_gt hqsz
      have hszne' : s₁.items.size ≠ vertItem curV := by rw [hsz₁]; exact Nat.ne_of_gt hvsz
      have h₃ : s₃.CloseInv := by
        apply h₂.modifyLoose (x := s₁.items.size)
          (by show s.g.nv < s₁.items.size; rw [hsz₁]; have := hs.size; omega) _
          (fun it => { it with vs := vsI }) (fun _ => rfl)
        apply cnt_eq_zero_of_free
        · intro t' ht' hm
          exact Nat.lt_irrefl _ (hsz₁ ▸ hs.span t' ht' _ hm)
        · intro p h
          have h' := (isParent_push_iff _ _ _ _).1 h
          rw [Items.IsParent, hch₁] at h'
          exact Nat.lt_irrefl _ (hsz₁ ▸ hs.ch_lt p _ h')
      set s₄ := { s₃ with tstack := s₃.tstack.tail } with hs₄
      have b₄ : BStep D s s₄ := b₃.trans
        (BStep.pop' b₃.inv b₃.shape (s₀ := s) t base hts rfl rfl
          (fun a i => (Items.Below_modify_vs _ _ _ _).trans
            ((Items.Below_push_nil _ rfl).trans (hB₁ a i)))
          (hok.gone hT t base hts))
      have h₄ : s₄.CloseInv := h₃.pop
      have hB₄ : ∀ a i, Items.Below s₄.items a i ↔ Items.Below s.items a i := fun a i =>
        (Items.Below_modify_vs _ _ _ _).trans
          ((Items.Below_push_nil _ rfl).trans (hB₁ a i))
      have hn₄ : Items.NoParent s₄.items (edgeItem s.g o.e) := fun p h =>
        hq_root₁ p ((isParent_push_iff _ _ _ _).1 ((isParent_modifyVs_iff _ _ _ _ _).1 h))
      have hvroot₄ : ∀ p, ¬ Items.IsParent s₄.items p (vertItem curV) := fun p h =>
        hv_root₁ p ((isParent_push_iff _ _ _ _).1 ((isParent_modifyVs_iff _ _ _ _ _).1 h))
      have hchI : Items.ch s₄.items s₁.items.size = [] := by
        show Items.ch ((s₁.items.push _).modify _ _) _ = []
        rw [Items.ch_modify_vs]; exact Items.ch_push_size _ rfl
      have htyI : Items.type s₄.items s₁.items.size = .I := by
        show Items.type ((s₁.items.push _).modify _ _) _ = _
        rw [Items.type_modify_vs]; exact Items.type_push_size _
      have hvsI' : Items.vs s₄.items s₁.items.size = vsI := by
        show Items.vs ((s₁.items.push _).modify _ _) _ = _
        rw [Items.vs_modify_self _ _ _ (by simp)]
      have hty₄ : ∀ j, j ≠ s₁.items.size → Items.type s₄.items j = Items.type s.items j := by
        intro j hj
        show Items.type ((s₁.items.push _).modify _ _) _ = _
        rw [Items.type_modify_vs, Items.type_push_of_ne _ hj, hty₁]
      have hvs₄ : ∀ j, j ≠ s₁.items.size → Items.vs s₄.items j = Items.vs s₁.items j := by
        intro j hj
        show Items.vs ((s₁.items.push _).modify _ _) _ = _
        rw [Items.vs_modify_ne _ _ _ _ (Ne.symm hj), Items.vs_push_of_ne _ hj]
      have hch₄ : ∀ j, Items.ch s₄.items j = Items.ch s.items j := by
        intro j
        show Items.ch ((s₁.items.push _).modify _ _) _ = _
        rw [Items.ch_modify_vs, Items.ch_push_nil _ rfl, hch₁]
      have hsz₄ : s₄.items.size = s.items.size + 1 := by simp [hs₄, hs₃, hs₂, hsz₁]
      have hqsz₄ : edgeItem s.g o.e < s₄.items.size := by rw [hsz₄]; exact Nat.lt_succ_of_lt hqsz
      set cs := s₁.items.size :: s.tstack.head!.spans.2 with hcs₀
      have hcs : cs = s₁.items.size :: t.spans.2 := by rw [hcs₀, hhead]
      set items₅ := s₄.items.modify (edgeItem s.g o.e) fun it => { it with ch := cs } with hitems₅
      have hch₅ : Items.ch items₅ (edgeItem s.g o.e) = cs := by
        rw [hitems₅, Items.ch_modify_self _ _ _ hqsz₄]
      have hpar₅ : ∀ c, Items.IsParent items₅ (edgeItem s.g o.e) c ↔ c ∈ cs := fun c => by
        rw [Items.IsParent, hch₅]
      have hbelow₅ : ∀ c ∈ t.spans.2, ∀ j, Items.Below items₅ c j ↔ Items.Below s.items c j := by
        intro c hcm j
        rw [hitems₅, Items.Below_modify_of_not_below (edgeItem s.g o.e) (fun it => { it with ch := cs })
          (fun hb => hok.q_free t ht (List.mem_append_right _ ((hn₄.below_eq hb) ▸ hcm)))]
        exact hB₄ c j
      have hbelowI : ∀ j, Items.Below items₅ s₁.items.size j → j = s₁.items.size := by
        intro j hj
        rcases hj.head_cases with h | ⟨c, hc, _⟩
        · exact h.symm
        · exfalso
          rw [Items.IsParent, hitems₅, Items.ch_modify_ne _ _ _ _ (Ne.symm hszne), hchI] at hc
          exact List.not_mem_nil hc
      have hEB : ∀ e', e' < s.g.ne →
          (Items.EdgeBelow s.g items₅ (edgeItem s.g o.e) e' ↔ e' = o.e ∨ t.edges s.g s.items e') := by
        intro e' he'
        constructor
        · intro hb
          rcases hb.head_cases with h | ⟨c, hc, hb'⟩
          · left; have : 1 + s.g.nv + o.e = 1 + s.g.nv + e' := h; omega
          · right
            rw [hpar₅, hcs, List.mem_cons] at hc
            rcases hc with rfl | hc
            · exfalso
              have h2 : 1 + s.g.nv + e' = s₁.items.size := hbelowI _ hb'
              rw [hsz₁] at h2; have := hs.size; omega
            · exact ⟨c, List.mem_append_right _ hc, (hbelow₅ c hc _).1 hb'⟩
        · rintro (rfl | ⟨c, hc, hb⟩)
          · exact .refl
          · rw [hside, List.nil_append] at hc
            exact .head ((hpar₅ c).2 (by rw [hcs]; exact List.mem_cons_of_mem _ hc))
              ((hbelow₅ c hc _).2 hb)
      have hsubt : ∀ e', e' < s.g.ne → s.g.Inc e' o.dest → e' = o.e ∨ t.edges s.g s.items e' := by
        intro e' he' hinc
        by_cases heq : e' = o.e
        · exact Or.inl heq
        · obtain ⟨t', ht', hte⟩ := hE.sub_cover e' he' heq (hE.dest_edges hT e' he' hinc)
          rw [hsub, List.mem_singleton] at ht'
          exact Or.inr (ht' ▸ hte)
      have hatt : ∀ v, Items.Att s.g items₅ (edgeItem s.g o.e) v → v = curV := by
        rintro v ⟨e₁, e₂, he₁, he₂, hi₁, hi₂, hb₁, hb₂⟩
        rw [hEB e₁ he₁] at hb₁
        rw [hEB e₂ he₂] at hb₂
        have hnd : v ≠ o.dest := fun hv => hb₂ (hsubt e₂ he₂ (hv ▸ hi₂))
        rcases hb₁ with rfl | hb₁
        · rcases eq_of_inc_pairEq hends hi₁ with h | h
          · exact absurd h hnd
          · exact h
        · have hE' := (hi.entries [] t base (by rw [hts]; rfl)).attached v e₁ e₂ he₁ he₂ hb₁
            (fun h => hb₂ (Or.inr h)) hi₁ hi₂
          rw [TEntry.Term'_nil] at hE'
          rcases hE' with h | ⟨k, hk₁, hk₂, hk⟩
          · exact absurd (h.trans hvs') hnd
          · rw [hdep] at hk₁
            rw [hc.hD, if_pos hT] at hk₂
            obtain rfl : k = d + 1 := by omega
            exact absurd (hk.trans (hE.sv_child hT)) hnd
      have hnew : Items.CloseAt s.g items₅ (edgeItem s.g o.e) := by
        refine Items.CloseAt.rootQ (c := s₁.items.size) helt ?_ (by rw [hch₅, hcs]; simp) ?_
          (Graph.inc_of_pairEq hends).2 hatt ?_ ?_ (fun h => absurd h hedge) ?_
        · rw [hitems₅, Items.type_modify_ch _ _ _, hty₄ _ (Ne.symm hszne)]
          exact hs.edge o.e helt
        · rw [hitems₅, Items.vs_modify_ch, hvs₄ _ (Ne.symm hszne)]; exact hvs₁
        · rw [hitems₅, Items.type_modify_ch _ _ _, htyI]; simp
        · rw [hitems₅, Items.type_modify_ch _ _ _, htyI]; intro h; cases h
        · intro _
          refine ⟨o.dest, hc.dest_lt, pairEq_swap' hends, by rw [hch₅, hcs, hvert], ?_⟩
          rw [hitems₅, Items.vs_modify_ch, hvsI', hvsI, hsv]
          cases s.stackDir[d]!
          · exact ⟨_, _, rfl, Or.inl rfl⟩
          · exact ⟨_, _, rfl, Or.inr rfl⟩
      have hcs' : ∀ c ∈ cs, Items.CloseAt s.g s₄.items c := by
        intro c hcm
        rw [hcs, List.mem_cons] at hcm
        rcases hcm with rfl | hcm
        · exact Items.CloseAt.leafIO (by rw [hsz₁]; exact hs.size) (Or.inl htyI) hchI
        · have hpos : 0 < s₃.cnt c := by
            have : 0 < spansCount s₃.tstack c := spansCount_pos_of_mem_head!
              (by show c ∈ s.tstack.head!.spans.1 ++ _; rw [hhead]; exact List.mem_append_right _ hcm)
            dsimp [cnt]; omega
          exact h₃.closed c (cnt_lt_size b₃.shape hpos) (Or.inr hpos)
      have h₅ := h₄.writeChildren hqsz₄ hn₄ cs hnew hcs'
      have b₅ : BStep D s { s₄ with items := items₅ } := b₄.trans
        (BStep.modifyCh b₄.inv b₄.shape (edgeItem s.g o.e) (fun it => { it with ch := cs }) hqlt
          (fun _ => rfl) hn₄ (fun t' ht' => hok.q_free t' (List.mem_of_mem_tail ht'))
          (fun _ c hc => by
            show c < s₄.items.size
            rcases List.mem_cons.1 hc with rfl | hc
            · rw [hsz₄, hsz₁]; exact Nat.lt_succ_self _
            · rw [hhead] at hc
              exact Nat.lt_of_lt_of_le (hs.span t ht c (List.mem_append_right _ hc))
                (by rw [hsz₄]; omega)))
      have hz : ({ s₄ with items := items₅ } : WalkState).cnt (vertItem curV) = 0 := by
        apply cnt_eq_zero_of_free
        · intro t' ht'; exact hok.v_free t' (List.mem_of_mem_tail ht')
        · intro p h
          rcases Items.IsParent_modify h with h | ⟨rfl, _, hcm⟩
          · exact hvroot₄ p h
          · have hcm' : vertItem curV ∈ cs := hcm
            rw [hcs] at hcm'
            rcases List.mem_cons.1 hcm' with h | h
            · exact hszne' h.symm
            · exact hok.v_free t ht (List.mem_append_right _ h)
      have h₆ := h₅.vertex_append' b₅.shape curV o.e hvlt helt hz hnew
        (by rw [hch₅, hcs]; simp)
      exact h₆
    · -- completed block
      rename_i hL
      have hL' : o.cls.lowval d ≠ d + 1 := by simpa using hL
      obtain ⟨t₁, t₂, hsub, hv₁, hd₁, hv₂, hd₂⟩ := hE.bd_comp hT hge hL'
      have hts : s.tstack = t₁ :: t₂ :: base := by rw [hE.tstack, hsub]; rfl
      have ht₁ : t₁ ∈ s.tstack := by rw [hts]; exact List.mem_cons_self ..
      have ht₂ : t₂ ∈ s.tstack := by rw [hts]; exact List.mem_cons_of_mem _ (List.mem_cons_self ..)
      have hsides := hE.bd_side hT hge
      rw [if_neg hL', hts] at hsides
      have hside₁ : t₁.spans.2 = [] := hsides.1 t₁ rfl
      have hside₂ : t₂.spans.1 = [] := hsides.2 t₂ rfl
      have hvert : t₂.spans.2 = [vertItem o.dest] := by
        have := hc.bd_vert hT hge; rw [if_neg hL', hts] at this; exact this t₂ rfl
      obtain ⟨c, hc₁, hct, hcq, a, b', hcvs, hcp⟩ := by
        have := hc.bd_node hT hge hL'; rw [hts] at this; exact this t₁ rfl
      have hcne : c ≠ edgeItem s.g o.e := fun h =>
        hok.q_free t₁ ht₁ (List.mem_append_left _ (by rw [hc₁, ← h]; exact List.mem_singleton_self _))
      have hhead₁ : s.tstack.head! = t₁ := by rw [hts]; rfl
      have hhead₂ : s.tstack.tail.head! = t₂ := by rw [hts]; rfl
      set s₂ := { s₁ with tstack := s₁.tstack.tail } with hs₂
      set s₃ := { s₂ with tstack := s₂.tstack.tail } with hs₃
      have b₂ : BStep D s s₂ := b₁.trans
        (BStep.pop' b₁.inv b₁.shape (s₀ := s) t₁ (t₂ :: base) hts rfl rfl hB₁
          (hok.gone hT t₁ (t₂ :: base) hts))
      have b₃ : BStep D s s₃ := b₂.trans
        (BStep.pop' b₂.inv b₂.shape (s₀ := s) t₂ base
          (by show s.tstack.tail = _; rw [hts, List.tail_cons]) rfl rfl hB₁
          (hok.gone₂ hT hL' t₁ t₂ base hts))
      have h₃ : s₃.CloseInv := h₁.pop.pop
      have hn₃ : Items.NoParent s₃.items (edgeItem s.g o.e) := hq_root₁
      have hqsz₃ : edgeItem s.g o.e < s₃.items.size := by rw [hsz₁]; exact hqsz
      set cs := s.tstack.head!.spans.1 ++ s.tstack.tail.head!.spans.2 with hcs₀
      have hcs : cs = t₁.spans.1 ++ t₂.spans.2 := by rw [hcs₀, hhead₁, hhead₂]
      have hcs₂ : cs = [c, vertItem o.dest] := by rw [hcs, hc₁, hvert]; rfl
      set items₄ := s₃.items.modify (edgeItem s.g o.e) fun it => { it with ch := cs } with hitems₄
      have hch₄ : Items.ch items₄ (edgeItem s.g o.e) = cs := by
        rw [hitems₄, Items.ch_modify_self _ _ _ hqsz₃]
      have hpar₄ : ∀ c', Items.IsParent items₄ (edgeItem s.g o.e) c' ↔ c' ∈ cs := fun c' => by
        rw [Items.IsParent, hch₄]
      have hbelow₄ : ∀ c' ∈ cs, ∀ j, Items.Below items₄ c' j ↔ Items.Below s.items c' j := by
        intro c' hcm j
        rw [hitems₄, Items.Below_modify_of_not_below (edgeItem s.g o.e) (fun it => { it with ch := cs })
          (fun hb => ?_)]
        · exact hB₁ c' j
        · have := hn₃.below_eq hb
          subst this
          rw [hcs, List.mem_append] at hcm
          rcases hcm with hcm | hcm
          · exact hok.q_free t₁ ht₁ (List.mem_append_left _ hcm)
          · exact hok.q_free t₂ ht₂ (List.mem_append_right _ hcm)
      have hEB : ∀ e', e' < s.g.ne → (Items.EdgeBelow s.g items₄ (edgeItem s.g o.e) e' ↔
          e' = o.e ∨ t₁.edges s.g s.items e' ∨ t₂.edges s.g s.items e') := by
        intro e' he'
        constructor
        · intro hb
          rcases hb.head_cases with h | ⟨c', hc', hb'⟩
          · left; have : 1 + s.g.nv + o.e = 1 + s.g.nv + e' := h; omega
          · right
            have hb'' := (hbelow₄ c' ((hpar₄ c').1 hc') _).1 hb'
            rw [hpar₄, hcs, List.mem_append] at hc'
            rcases hc' with hc' | hc'
            · exact Or.inl ⟨c', List.mem_append_left _ hc', hb''⟩
            · exact Or.inr ⟨c', List.mem_append_right _ hc', hb''⟩
        · rintro (rfl | ⟨c', hc', hb⟩ | ⟨c', hc', hb⟩)
          · exact .refl
          · rw [hside₁, List.append_nil] at hc'
            have hm : c' ∈ cs := by rw [hcs]; exact List.mem_append_left _ hc'
            exact .head ((hpar₄ c').2 hm) ((hbelow₄ c' hm _).2 hb)
          · rw [hside₂, List.nil_append] at hc'
            have hm : c' ∈ cs := by rw [hcs]; exact List.mem_append_right _ hc'
            exact .head ((hpar₄ c').2 hm) ((hbelow₄ c' hm _).2 hb)
      have hsubt : ∀ e', e' < s.g.ne → s.g.Inc e' o.dest →
          e' = o.e ∨ t₁.edges s.g s.items e' ∨ t₂.edges s.g s.items e' := by
        intro e' he' hinc
        by_cases heq : e' = o.e
        · exact Or.inl heq
        · obtain ⟨t', ht', hte⟩ := hE.sub_cover e' he' heq (hE.dest_edges hT e' he' hinc)
          rw [hsub, List.mem_cons, List.mem_singleton] at ht'
          rcases ht' with rfl | rfl
          · exact Or.inr (Or.inl hte)
          · exact Or.inr (Or.inr hte)
      have hterm₁ : ∀ v, t₁.Term D s v → v ≠ o.dest → v = curV := by
        intro v h hnd
        rcases h with h | ⟨k, hk₁, hk₂, hk⟩
        · exact absurd (h.trans hv₁) hnd
        · rw [hd₁] at hk₁
          rw [hc.hD, if_pos hT] at hk₂
          rcases Nat.lt_or_ge k (d + 1) with hk' | hk'
          · obtain rfl : k = d := by omega
            exact hk.trans hsv
          · obtain rfl : k = d + 1 := by omega
            exact absurd (hk.trans (hE.sv_child hT)) hnd
      have hterm₂ : ∀ v, t₂.Term D s v → v ≠ o.dest → False := by
        intro v h hnd
        rcases h with h | ⟨k, hk₁, hk₂, hk⟩
        · exact hnd (h.trans hv₂)
        · rw [hd₂] at hk₁
          rw [hc.hD, if_pos hT] at hk₂
          obtain rfl : k = d + 1 := by omega
          exact hnd (hk.trans (hE.sv_child hT))
      have hatt : ∀ v, Items.Att s.g items₄ (edgeItem s.g o.e) v → v = curV := by
        rintro v ⟨e₁, e₂, he₁, he₂, hi₁, hi₂, hb₁, hb₂⟩
        rw [hEB e₁ he₁] at hb₁
        rw [hEB e₂ he₂] at hb₂
        have hnd : v ≠ o.dest := fun hv => hb₂ (hsubt e₂ he₂ (hv ▸ hi₂))
        rcases hb₁ with rfl | hb₁ | hb₁
        · rcases eq_of_inc_pairEq hends hi₁ with h | h
          · exact absurd h hnd
          · exact h
        · have hE' := (hi.entries [] t₁ (t₂ :: base) (by rw [hts]; rfl)).attached v e₁ e₂ he₁ he₂
            hb₁ (fun h => hb₂ (Or.inr (Or.inl h))) hi₁ hi₂
          rw [TEntry.Term'_nil] at hE'
          exact hterm₁ v hE' hnd
        · have hE' := (hi.entries [t₁] t₂ base (by rw [hts]; rfl)).attached v e₁ e₂ he₁ he₂
            hb₁ (fun h => hb₂ (Or.inr (Or.inr h))) hi₁ hi₂
          rcases hE' with h | ⟨t', ht', h⟩
          · exact (hterm₂ v h hnd).elim
          · rw [List.mem_singleton] at ht'
            subst ht'
            exact hterm₁ v h hnd
      have hnew : Items.CloseAt s.g items₄ (edgeItem s.g o.e) := by
        refine Items.CloseAt.rootQ (c := c) helt ?_ (by rw [hch₄, hcs₂]; simp) ?_
          (Graph.inc_of_pairEq hends).2 hatt ?_ ?_ (fun h => absurd h hedge) ?_
        · rw [hitems₄, Items.type_modify_ch _ _ _, hty₁]; exact hs.edge o.e helt
        · rw [hitems₄, Items.vs_modify_ch]; exact hvs₁
        · rw [hitems₄, Items.type_modify_ch _ _ _, hty₁]; exact hct
        · rw [hitems₄, Items.type_modify_ch _ _ _, hty₁,
            Items.ch_modify_ne _ _ _ _ (Ne.symm hcne), hch₁]
          exact hcq
        · intro _
          refine ⟨o.dest, hc.dest_lt, pairEq_swap' hends, by rw [hch₄, hcs₂], a, b', ?_, hcp⟩
          rw [hitems₄, Items.vs_modify_ch, hvs₁' c hcne]; exact hcvs
      have hcs' : ∀ c' ∈ cs, Items.CloseAt s.g s₃.items c' := by
        intro c' hcm
        have hpos : 0 < s₁.cnt c' := by
          rw [hcs, List.mem_append] at hcm
          rcases hcm with hcm | hcm
          · have : 0 < spansCount s₁.tstack c' := spansCount_pos_of_mem_head!
              (by show c' ∈ s.tstack.head!.spans.1 ++ _; rw [hhead₁]; exact List.mem_append_left _ hcm)
            dsimp [cnt]; omega
          · have : 0 < spansCount s₁.tstack c' := spansCount_pos_of_mem_tail_head!
              (by show c' ∈ s.tstack.tail.head!.spans.1 ++ _; rw [hhead₂]
                  exact List.mem_append_right _ hcm)
            dsimp [cnt]; omega
        exact h₁.closed c' (cnt_lt_size b₁.shape hpos) (Or.inr hpos)
      have h₄ := h₃.writeChildren hqsz₃ hn₃ cs hnew hcs'
      have b₄ : BStep D s { s₃ with items := items₄ } := b₃.trans
        (BStep.modifyCh b₃.inv b₃.shape (edgeItem s.g o.e) (fun it => { it with ch := cs }) hqlt
          (fun _ => rfl) hn₃
          (fun t' ht' => hok.q_free t' (List.mem_of_mem_tail (List.mem_of_mem_tail ht')))
          (fun _ c' hc' => by
            show c' < s₁.items.size
            rw [hsz₁]
            have hc'' : c' ∈ cs := hc'
            rw [hcs] at hc''
            rcases List.mem_append.1 hc'' with hc' | hc'
            · exact hs.span t₁ ht₁ c' (List.mem_append_left _ hc')
            · exact hs.span t₂ ht₂ c' (List.mem_append_right _ hc')))
      have hz : ({ s₃ with items := items₄ } : WalkState).cnt (vertItem curV) = 0 := by
        apply cnt_eq_zero_of_free
        · intro t' ht'
          exact hok.v_free t' (List.mem_of_mem_tail (List.mem_of_mem_tail ht'))
        · intro p h
          rcases Items.IsParent_modify h with h | ⟨rfl, _, hcm⟩
          · exact hv_root₁ p h
          · have hcm' : vertItem curV ∈ cs := hcm
            rw [hcs] at hcm'
            rcases List.mem_append.1 hcm' with h | h
            · exact hok.v_free t₁ ht₁ (List.mem_append_left _ h)
            · exact hok.v_free t₂ ht₂ (List.mem_append_right _ h)
      have h₅ := h₄.vertex_append' b₄.shape curV o.e hvlt helt hz hnew
        (by rw [hch₄, hcs₂]; simp)
      exact h₅
  · -- self-loop
    rename_i hT
    have hT' : o.cls.isTree = false := by simpa using hT
    have hends := hc.ends_back hT'
    have hloop := hc.bd_loop hT' hge
    have hedge : (s.g.edges[o.e]!).1 = (s.g.edges[o.e]!).2 := by
      rw [hloop] at hends
      generalize s.g.edges[o.e]! = q at hends
      obtain ⟨x, y⟩ := q
      simp only [Items.PairEq, Prod.mk.injEq] at hends
      omega
    set s₂ := { s₁ with totSelfLoops := s₁.totSelfLoops + 1 } with hs₂
    have b₂ : BStep D s s₂ := b₁.trans (BStep.frame b₁.inv b₁.shape)
    have h₂ : s₂.CloseInv := h₁.frame rfl rfl fun _ h => h
    set s₃ := { s₂ with items := s₂.items.push ⟨.O, (none, none), []⟩ } with hs₃
    have b₃ : BStep D s s₃ := b₂.trans (BStep.alloc b₂.inv b₂.shape .O)
    have h₃ : s₃.CloseInv := h₂.alloc' b₂.shape .O
    set s₄ := { s₃ with items := s₃.items.modify s₂.items.size fun it => { it with vs := (some curV, none) } }
      with hs₄
    have hsz₂ : s₂.items.size = s.items.size := hsz₁
    have b₄ : BStep D s s₄ := b₃.trans (BStep.modifyVs_leaf b₃.inv b₃.shape s₂.items.size
      (some curV, none) (Items.ch_push_size _ rfl))
    have hszne : s₂.items.size ≠ edgeItem s.g o.e := by rw [hsz₂]; exact Nat.ne_of_gt hqsz
    have hszne' : s₂.items.size ≠ vertItem curV := by rw [hsz₂]; exact Nat.ne_of_gt hvsz
    have h₄ : s₄.CloseInv := by
      apply h₃.modifyLoose (x := s₂.items.size)
        (by show s.g.nv < s₂.items.size; rw [hsz₂]; have := hs.size; omega) _
        (fun it => { it with vs := (some curV, none) }) (fun _ => rfl)
      apply cnt_eq_zero_of_free
      · intro t' ht' hm
        exact Nat.lt_irrefl _ (hsz₂ ▸ hs.span t' ht' _ hm)
      · intro p h
        have h' := (isParent_push_iff _ _ _ _).1 h
        rw [Items.IsParent, hch₁] at h'
        exact Nat.lt_irrefl _ (hsz₂ ▸ hs.ch_lt p _ h')
    have hn₄ : Items.NoParent s₄.items (edgeItem s.g o.e) := fun p h =>
      hq_root₁ p ((isParent_push_iff _ _ _ _).1 ((isParent_modifyVs_iff _ _ _ _ _).1 h))
    have hvroot₄ : ∀ p, ¬ Items.IsParent s₄.items p (vertItem curV) := fun p h =>
      hv_root₁ p ((isParent_push_iff _ _ _ _).1 ((isParent_modifyVs_iff _ _ _ _ _).1 h))
    have hchO : Items.ch s₄.items s₂.items.size = [] := by
      show Items.ch ((s₁.items.push _).modify _ _) _ = []
      rw [Items.ch_modify_vs]; exact Items.ch_push_size _ rfl
    have htyO : Items.type s₄.items s₂.items.size = .O := by
      show Items.type ((s₁.items.push _).modify _ _) _ = _
      rw [Items.type_modify_vs]; exact Items.type_push_size _
    have hvsO : Items.vs s₄.items s₂.items.size = (some curV, none) := by
      show Items.vs ((s₁.items.push _).modify _ _) _ = _
      rw [Items.vs_modify_self _ _ _ (by show s₁.items.size < (s₁.items.push _).size; simp)]
    have hty₄ : ∀ j, j ≠ s₂.items.size → Items.type s₄.items j = Items.type s.items j := by
      intro j hj
      show Items.type ((s₁.items.push _).modify _ _) _ = _
      rw [Items.type_modify_vs, Items.type_push_of_ne _ hj, hty₁]
    have hvs₄ : ∀ j, j ≠ s₂.items.size → Items.vs s₄.items j = Items.vs s₁.items j := by
      intro j hj
      show Items.vs ((s₁.items.push _).modify _ _) _ = _
      rw [Items.vs_modify_ne _ _ _ _ (Ne.symm hj), Items.vs_push_of_ne _ hj]
    have hsz₄ : s₄.items.size = s.items.size + 1 := by simp [hs₄, hs₃, hs₂, hsz₁]
    have hqsz₄ : edgeItem s.g o.e < s₄.items.size := by rw [hsz₄]; exact Nat.lt_succ_of_lt hqsz
    set cs := [s₂.items.size] with hcs
    set items₅ := s₄.items.modify (edgeItem s.g o.e) fun it => { it with ch := cs } with hitems₅
    have hch₅ : Items.ch items₅ (edgeItem s.g o.e) = cs := by
      rw [hitems₅, Items.ch_modify_self _ _ _ hqsz₄]
    have hEB : ∀ e', e' < s.g.ne →
        Items.EdgeBelow s.g items₅ (edgeItem s.g o.e) e' → e' = o.e := by
      intro e' he' hb
      rcases hb.head_cases with h | ⟨c, hc, hb'⟩
      · have : 1 + s.g.nv + o.e = 1 + s.g.nv + e' := h; omega
      · exfalso
        rw [Items.IsParent, hch₅, hcs, List.mem_singleton] at hc
        subst hc
        rcases hb'.head_cases with h | ⟨c', hc', _⟩
        · have h2 : 1 + s.g.nv + e' = s₂.items.size := h.symm
          rw [hsz₂] at h2; have := hs.size; omega
        · rw [Items.IsParent, hitems₅, Items.ch_modify_ne _ _ _ _ (Ne.symm hszne), hchO] at hc'
          exact List.not_mem_nil hc'
    have hatt : ∀ v, Items.Att s.g items₅ (edgeItem s.g o.e) v → v = curV := by
      rintro v ⟨e₁, e₂, he₁, he₂, hi₁, hi₂, hb₁, _⟩
      obtain rfl := hEB e₁ he₁ hb₁
      rcases eq_of_inc_pairEq hends hi₁ with h | h
      · exact h
      · exact h.trans hloop
    have hnew : Items.CloseAt s.g items₅ (edgeItem s.g o.e) := by
      refine Items.CloseAt.rootQ (c := s₂.items.size) helt ?_ (by rw [hch₅, hcs]; simp) ?_
        (Graph.inc_of_pairEq hends).1 hatt ?_ ?_ ?_ (fun h => absurd hedge h)
      · rw [hitems₅, Items.type_modify_ch _ _ _, hty₄ _ (Ne.symm hszne)]
        exact hs.edge o.e helt
      · rw [hitems₅, Items.vs_modify_ch, hvs₄ _ (Ne.symm hszne)]; exact hvs₁
      · rw [hitems₅, Items.type_modify_ch _ _ _, htyO]; simp
      · rw [hitems₅, Items.type_modify_ch _ _ _, htyO]; intro h; cases h
      · intro _
        refine ⟨by rw [hch₅, hcs], ?_⟩
        rw [hitems₅, Items.vs_modify_ch]; exact hvsO
    have hcs' : ∀ c ∈ cs, Items.CloseAt s.g s₄.items c := by
      intro c hcm
      rw [hcs, List.mem_singleton] at hcm
      subst hcm
      exact Items.CloseAt.leafIO (by rw [hsz₂]; exact hs.size) (Or.inr htyO) hchO
    have h₅ := h₄.writeChildren hqsz₄ hn₄ cs hnew hcs'
    have b₅ : BStep D s { s₄ with items := items₅ } := b₄.trans
      (BStep.modifyCh b₄.inv b₄.shape (edgeItem s.g o.e) (fun it => { it with ch := cs }) hqlt
        (fun _ => rfl) hn₄ (fun t' ht' => hok.q_free t' ht')
        (fun _ c hc => by
          show c < s₄.items.size
          rw [hsz₄, List.mem_singleton.1 hc, hsz₂]; exact Nat.lt_succ_self _))
    have hz : ({ s₄ with items := items₅ } : WalkState).cnt (vertItem curV) = 0 := by
      apply cnt_eq_zero_of_free
      · intro t' ht'; exact hok.v_free t' ht'
      · intro p h
        rcases Items.IsParent_modify h with h | ⟨rfl, _, hcm⟩
        · exact hvroot₄ p h
        · exact hszne' (List.mem_singleton.1 hcm).symm
    have h₆ := h₅.vertex_append' b₅.shape curV o.e hvlt helt hz hnew
      (by rw [hch₅, hcs]; simp)
    exact h₆

end Spqr.WalkState
