import Spqr.Proofs.RInvFrame
import Spqr.WalkInv

/-!
# The schedule frontier from the ear contract

`finishEdge_frontier` establishes the R session's `WalkState.Frontier` (`Spqr/Proofs/RInvFrame.lean`)
at a `finishEdge` site from `FinishBook.ear` (`EarAt`): `size`/`owns`/`base_disj` are the
`EarFinish` ownership fields; `loop1` follows the loop-1 iterates through `L1Inv` (`l1_init`,
`l1_iter`), `loop2` places the loop-2 iterates between `EarLate` (`feS₁`) and `EarClose` (`feS₂`),
and `loop3` folds the `feS₂` shape. `walkTree_frontiers` gives `FrontiersTree` by the walk
induction, under the hypotheses of `walkTree_inv'`.
-/

namespace Spqr
open WalkM

namespace WalkState

theorem frontierOwns_of_split {origTstack : Nat} {base front : List TEntry} {E : Nat → Prop}
    {s : WalkState} (hts : s.tstack = front ++ base) (hl : base.length = origTstack)
    (h : ∀ e, e < s.g.ne → ((∃ t ∈ front, t.edges s.g s.items e) ↔ E e)) :
    FrontierOwns origTstack base E s := by
  have hlen : s.tstack.length - origTstack = front.length := by
    rw [hts, List.length_append, hl]; omega
  refine ⟨by rw [hts, List.length_append, hl]; omega, ?_, ?_⟩
  · rw [hlen, hts, List.drop_left]
  · rw [hlen, hts, List.take_left]; exact h

/-- Merging the top two entries `k` times inside the front segment keeps `FrontierOwns`. -/
theorem frontierOwns_iter_merge {origTstack : Nat} {base R' : List TEntry} {c₀ : TEntry}
    {E : Nat → Prop} {st : WalkState} (hts : st.tstack = c₀ :: R' ++ base) (hl : base.length = origTstack)
    (h : ∀ e, e < st.g.ne → ((∃ t ∈ c₀ :: R', t.edges st.g st.items e) ↔ E e)) (k : Nat)
    (hk : k ≤ R'.length) :
    FrontierOwns origTstack base E (iter mergeTstackTops k st) ∧
      (iter mergeTstackTops k st).tstack.length = R'.length + 1 - k + origTstack := by
  have hts' : st.tstack = c₀ :: (R' ++ base) := hts
  rw [iter_merge_eq st c₀ (R' ++ base) hts' k (by rw [List.length_append]; omega),
    List.take_append_of_le_length hk, List.drop_append_of_le_length hk]
  refine ⟨frontierOwns_of_split (front := l2Cur c₀ (R'.take k) :: R'.drop k) rfl hl ?_, ?_⟩
  · intro e he
    rw [← h e he]
    constructor
    · rintro ⟨t, ht, hte⟩
      rcases List.mem_cons.1 ht with rfl | ht
      · rcases (l2Cur_edges c₀ (R'.take k) e).1 hte with hte | ⟨u, hu, hue⟩
        · exact ⟨c₀, by simp, hte⟩
        · exact ⟨u, List.mem_cons_of_mem _ (List.mem_of_mem_take hu), hue⟩
      · exact ⟨t, List.mem_cons_of_mem _ (List.mem_of_mem_drop ht), hte⟩
    · rintro ⟨t, ht, hte⟩
      rcases List.mem_cons.1 ht with rfl | ht
      · exact ⟨l2Cur t (R'.take k), by simp, (l2Cur_edges _ _ _).2 (.inl hte)⟩
      · have : t ∈ R'.take k ++ R'.drop k := by rw [List.take_append_drop]; exact ht
        rcases List.mem_append.1 this with h1 | h1
        · exact ⟨l2Cur c₀ (R'.take k), by simp, (l2Cur_edges _ _ _).2 (.inr ⟨t, h1, hte⟩)⟩
        · exact ⟨t, List.mem_cons_of_mem _ h1, hte⟩
  · simp only [List.length_cons, List.length_append, List.length_drop, hl]; omega

/-- The ear contract at a `finishEdge` site gives the R session's schedule frontier. -/
theorem finishEdge_frontier {D curV d : Nat} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}
    {s : WalkState} (hb : FinishBook curV d o origTstack hasVert s) (hi : s.Inv' D) (hs : Shape s)
    (hD : D = if o.cls.isTree then d + 1 else d) : Frontier (o := o) d origTstack s := by
  obtain ⟨sub, base, hl, hE⟩ := hb.ear
  have hdrop : s.tstack.drop (s.tstack.length - origTstack) = base := by
    rw [hE.tstack, List.length_append, hl, Nat.add_sub_cancel, List.drop_left]
  have htake : s.tstack.take (s.tstack.length - origTstack) = sub := by
    rw [hE.tstack, List.length_append, hl, Nat.add_sub_cancel, List.take_left]
  refine ⟨?size, ?owns, ?base_disj, ?loop1, ?loop2, ?loop3⟩
  case size => rw [hE.tstack, List.length_append, hl]; omega
  case owns =>
    intro e he
    rw [htake]
    constructor
    · rintro (rfl | ⟨t, ht, hte⟩)
      · exact Or.inl rfl
      · exact hE.sub_edges t ht e he hte
    · intro hse
      by_cases heq : e = o.e
      · exact Or.inl heq
      · exact Or.inr (hE.sub_cover e he heq hse)
  case base_disj => rw [hdrop]; exact hE.base_disj
  case loop1 =>
    intro ht hlow k hk
    dsimp only
    have hD' : D = d + 1 := by rw [hD, if_pos ht]
    obtain ⟨hi', lo, hsplit, hc⟩ := L1Ctx.ofEar hE hs hD' ht hlow
    have hends : Items.PairEq (o.dest, s.stackVerts[d]!) s.g.edges[o.e]! := by
      obtain ⟨lv, kind, ho⟩ : ∃ lv kind, o.cls = .ret lv kind := by
        cases h : o.cls <;> simp [OutClass.lowval, h] at hlow ⊢; omega
      simpa [ht] using hb.ends lv kind ho
    obtain ⟨done, rest, c, hreach, hP, hK, hF⟩ :=
      l1_iter hc (l1_init hE hi hs hD' hb.e_lt hb.q hends hsplit) hb.v_lt k hk
    have hsub : sub = done ++ rest ++ lo := by rw [hsplit, hreach.append]
    have hmemR : ∀ t, t ∈ rest ++ lo → t ∈ rest ++ lo ++ base := fun t ht =>
      List.mem_append_left _ ht
    have hmemS : ∀ t, t ∈ rest ++ lo → t ∈ sub := fun t ht => by
      rw [hsub, List.append_assoc]; exact List.mem_append_right _ ht
    refine ⟨?_, ?_⟩
    · rw [hdrop]
      refine frontierOwns_of_split (front := c :: (rest ++ lo)) (by rw [hP.tstack]; simp) hl ?_
      intro e he
      rw [hF.g] at he ⊢
      constructor
      · rintro ⟨t, ht', hte⟩
        rcases List.mem_cons.1 ht' with rfl | ht'
        · rcases (hP.edges e he).1 hte with heq | ⟨u, hu, hue⟩
          · exact Or.inl heq
          · exact hE.sub_edges u (by rw [hsub]; simp [hu]) e he hue
        · rw [hK.edges (hmemR t ht')] at hte
          exact hE.sub_edges t (hmemS t ht') e he hte
      · intro hse
        by_cases heq : e = o.e
        · exact ⟨c, by simp, (hP.edges e he).2 (Or.inl heq)⟩
        · obtain ⟨t, ht', hte⟩ := hE.sub_cover e he heq hse
          rw [hsub, List.append_assoc] at ht'
          rcases List.mem_append.1 ht' with hd | hr
          · exact ⟨c, by simp, (hP.edges e he).2 (Or.inr ⟨t, hd, hte⟩)⟩
          · exact ⟨t, List.mem_cons_of_mem _ hr, (hK.edges (hmemR t hr) e).2 hte⟩
    · intro hcond
      rw [result_loop1Cond_iff, hP.tstack] at hcond
      obtain ⟨l0, lo', hlo⟩ : ∃ l0 lo', lo = l0 :: lo' := by
        cases lo with
        | nil => exact absurd rfl hc.lo_ne
        | cons a b => exact ⟨a, b, rfl⟩
      have hl0 := hc.lo_top l0 (by rw [hlo]; rfl)
      cases rest with
      | nil => simp [hlo] at hcond; omega
      | cons t rest' =>
        rw [hP.tstack]
        simp only [List.length_cons, List.length_append, hl, hlo]
        split <;> omega
  case loop2 =>
    intro ht hlow k hk
    dsimp only
    rw [hdrop]
    obtain ⟨c₀, R, hL⟩ := hE.late ht hlow
    obtain ⟨c, mid, py, vy, hC⟩ := hE.close ht hlow
    have hfeS₂ : feS₂ d o s = ((mergeLate d).run (feS₁ d o s)).2 := rfl
    have hrun := mergeLate_run d (feS₁ d o s)
    have hlen₂ : (feS₂ d o s).tstack.length = mid.length + 3 + origTstack := by
      rw [hC.tstack]; simp [hl]; omega
    have hown₂ : ∀ e, e < (feS₂ d o s).g.ne →
        ((∃ t ∈ c :: (mid ++ [py, vy]), t.edges (feS₂ d o s).g (feS₂ d o s).items e) ↔ subEdges o e) := by
      intro e he
      rw [hC.g] at he ⊢
      constructor
      · rintro ⟨t, ht', hte⟩
        exact hC.sub_edges t (by simpa using ht') e he hte
      · intro hse
        obtain ⟨t, ht', hte⟩ := hC.sub_cover e he hse
        exact ⟨t, by simpa using ht', hte⟩
    by_cases hc₀ : (curE (feS₁ d o s)).firstIdx > (feS₁ d o s).firstOccurrence[d]!
    · rw [if_pos hc₀] at hrun
      obtain ⟨k₂, hk₂, heq, hall, hend⟩ := loop_run_iter (feS₁ d o s).tstack.length
        (loop2Cond (feS₁ d o s).firstOccurrence[d]!) mergeTstackTops (feS₁ d o s)
        (fun s => by rw [run_loop2Cond])
      have hS₂ : feS₂ d o s = iter mergeTstackTops k₂ (feS₁ d o s) := by rw [hfeS₂, hrun]; exact heq
      have hk₂R : k₂ ≤ R.length := by
        by_contra hgt
        have h1 : (iter mergeTstackTops R.length (feS₁ d o s)).tstack.length ≤ 1 := by
          rw [iter_merge_eq _ c₀ R hL.tstack R.length (Nat.le_refl _)]; simp
        have h2 := iter_merge_short _ h1 (k₂ - R.length)
        rw [← iter_add, Nat.sub_add_cancel (by omega), ← hS₂, hlen₂] at h2
        omega
      have hiter := iter_merge_eq _ c₀ R hL.tstack k₂ hk₂R
      rw [hS₂, hiter] at hlen₂ hown₂
      have hts₂ : l2Cur c₀ (R.take k₂) :: R.drop k₂ = c :: (mid ++ [py, vy] ++ base) := by
        have h := hC.tstack
        rw [hS₂, hiter] at h
        exact h
      have hc : l2Cur c₀ (R.take k₂) = c := (List.cons_eq_cons.1 hts₂).1
      have hR : R.drop k₂ = mid ++ [py, vy] ++ base := (List.cons_eq_cons.1 hts₂).2
      have hRsplit : R = (R.take k₂ ++ (mid ++ [py, vy])) ++ base := by
        rw [List.append_assoc, ← hR, List.take_append_drop]
      have hts₁ : (feS₁ d o s).tstack = c₀ :: (R.take k₂ ++ (mid ++ [py, vy])) ++ base := by
        rw [hL.tstack]; exact congrArg (List.cons c₀) hRsplit
      have hown₁ : ∀ e, e < (feS₁ d o s).g.ne →
          ((∃ t ∈ c₀ :: (R.take k₂ ++ (mid ++ [py, vy])), t.edges (feS₁ d o s).g (feS₁ d o s).items e) ↔
            subEdges o e) := by
        intro e he
        rw [← hown₂ e he]
        constructor
        · rintro ⟨t, ht', hte⟩
          rcases List.mem_cons.1 ht' with rfl | ht'
          · exact ⟨c, by simp, by rw [← hc]; exact (l2Cur_edges _ _ _).2 (.inl hte)⟩
          · rcases List.mem_append.1 ht' with h1 | h1
            · exact ⟨c, by simp, by rw [← hc]; exact (l2Cur_edges _ _ _).2 (.inr ⟨t, h1, hte⟩)⟩
            · exact ⟨t, List.mem_cons_of_mem _ h1, hte⟩
        · rintro ⟨t, ht', hte⟩
          rcases List.mem_cons.1 ht' with rfl | ht'
          · rw [← hc] at hte
            rcases (l2Cur_edges _ _ _).1 hte with hte | ⟨u, hu, hue⟩
            · exact ⟨c₀, by simp, hte⟩
            · exact ⟨u, by simp [hu], hue⟩
          · exact ⟨t, by simp [ht'], hte⟩
      have hkk : k ≤ k₂ := by
        by_contra hgt
        rcases hend with hfuel | hfalse
        · have := hk₂R
          rw [hfuel, hL.tstack] at this
          simp only [List.length_cons] at this
          omega
        · have := hk k₂ (by omega)
          unfold result at this
          rw [this] at hfalse
          exact Bool.noConfusion hfalse
      obtain ⟨h1, h2⟩ := frontierOwns_iter_merge hts₁ hl hown₁ k
        (by simp only [List.length_append, List.length_take, List.length_cons, List.length_nil]; omega)
      refine ⟨h1, fun _ => ?_⟩
      rw [h2]
      simp only [List.length_append, List.length_cons, List.length_nil, List.length_take,
        Nat.min_eq_left hk₂R]
      omega
    · rw [if_neg hc₀] at hrun
      have hS₂ : feS₂ d o s = feS₁ d o s := by rw [hfeS₂, hrun]
      have hk0 : k = 0 := by
        by_contra hne
        have := hk 0 (Nat.pos_of_ne_zero hne)
        rw [result_loop2Cond_iff] at this
        exact hc₀ this
      subst hk0
      show FrontierOwns origTstack base (subEdges o) (feS₁ d o s) ∧
        (result (loop2Cond (feS₁ d o s).firstOccurrence[d]!) (feS₁ d o s) = true →
          origTstack + 2 ≤ (feS₁ d o s).tstack.length)
      rw [← hS₂]
      refine ⟨frontierOwns_of_split (front := c :: (mid ++ [py, vy])) (by rw [hC.tstack]; simp) hl hown₂,
        fun _ => by rw [hlen₂]; omega⟩
  case loop3 =>
    intro ht hlow _ k hk
    dsimp only
    rw [hdrop]
    obtain ⟨c, mid, py, vy, hC⟩ := hE.close ht hlow
    have hts : (feS₂ d o s).tstack = c :: (mid ++ [py, vy]) ++ base := by rw [hC.tstack]; simp
    have hown : ∀ e, e < (feS₂ d o s).g.ne →
        ((∃ t ∈ c :: (mid ++ [py, vy]), t.edges (feS₂ d o s).g (feS₂ d o s).items e) ↔ subEdges o e) := by
      intro e he
      rw [hC.g] at he ⊢
      constructor
      · rintro ⟨t, ht', hte⟩
        exact hC.sub_edges t (by simpa using ht') e he hte
      · intro hse
        obtain ⟨t, ht', hte⟩ := hC.sub_cover e he hse
        exact ⟨t, by simpa using ht', hte⟩
    have hk' : k ≤ mid.length + 2 := by
      by_contra hgt
      have := hk (mid.length + 2) (by omega)
      unfold result at this
      rw [run_loop3Cond] at this
      simp only [decide_eq_true_eq] at this
      rw [(frontierOwns_iter_merge hts hl hown (mid.length + 2) (by simp)).2] at this
      simp only [List.length_append, List.length_cons, List.length_nil] at this
      omega
    refine ⟨(frontierOwns_iter_merge hts hl hown k
      (by simp only [List.length_append, List.length_cons, List.length_nil]; omega)).1, fun hcond => ?_⟩
    unfold result at hcond
    rw [run_loop3Cond] at hcond
    simp only [decide_eq_true_eq] at hcond
    omega

/-! `FrontiersTree` by the walk induction (the shape of `invTree`/`invOuts`/`invOut`): `Inv'`/`Shape`
are threaded to each `finishEdge` site, where `finishEdge_frontier` applies. -/

abbrev FrTree (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  (∀ v outs, t = .node v outs → ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).Inv' d) →
  Shape s → GuardsTree t d s → BookTree t d s → FrontiersTree t d s

abbrev FrOuts (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  s.Inv' d → Shape s → GuardsOuts v d outs hasVert s → BookOuts v d outs hasVert s →
  FrontiersOuts v d outs hasVert s

abbrev FrOut (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  s.Inv' d → Shape s → GuardsOut v d o hasVert s → BookOut v d o hasVert s →
  FrontiersOut v d o hasVert s

mutual
theorem frTree : ∀ (t : DfsTree) (d : Nat) (s : WalkState), FrTree t d s
  | .node v outs, d, s => fun hi hs hg hb => by
    unfold FrontiersTree; unfold GuardsTree at hg; unfold BookTree at hb
    exact frOuts v d outs false _ (hi v outs rfl) hs.frame' hg hb

theorem frOuts : ∀ (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState),
    FrOuts v d outs hasVert s
  | v, d, [], hasVert, s => fun _ _ _ _ => by
    unfold FrontiersOuts; trivial
  | v, d, o :: rest, hasVert, s => fun hi hs hg hb => by
    unfold FrontiersOuts; unfold GuardsOuts at hg; unfold BookOuts at hb
    refine ⟨frOut v d o hasVert s hi hs hg.1 hb.1, ?_⟩
    exact wp_imp (wp_imp (wp_imp (wp_of_forall fun hv' s' ⟨hi', hs'⟩ hg' hb' =>
      frOuts v d rest hv' s' hi' hs' hg' hb') (invOut v d o hasVert s hi hs hg.1 hb.1)) hg.2) hb.2

theorem frOut : ∀ (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState), FrOut v d o hasVert s
  | v, d, o, hasVert, s => fun hi hs hg hb => by
    unfold FrontiersOut; unfold GuardsOut at hg; unfold BookOut at hb
    refine wp_imp (wp_imp (wp_imp (wp_of_forall fun hv' s₁ ⟨hi₁, hs₁⟩ hg₁ hb₁ => ?_)
      (walkOutPre_inv hi hs hb.1)) hg) hb.2
    cases o with
    | back e cls dest =>
      exact finishEdge_frontier (D := d) hb₁ hi₁ hs₁
        (by rw [if_neg fun h => by obtain ⟨_, _, _, h⟩ := hb₁.tree.1 h; cases h])
    | tree e cls child =>
      try simp only [wp_modify] at hg₁ hb₁
      simp only [wp_modify]
      have pre : ∀ w outs, child = .node w outs →
          ({ { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne } with
            stackVerts := s₁.stackVerts.set! (d + 1) w } : WalkState).Inv' (d + 1) :=
        fun w outs _ => (hi₁.frame' (s' := { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne })).setSv w
      refine ⟨frTree child (d + 1) _ pre hs₁.frame' hg₁.1 hb₁.1, ?_⟩
      exact wp_imp (wp_imp (wp_of_forall fun _ s₃ ⟨hg₃, hb₃⟩ ⟨hi₃, hs₃⟩ =>
        finishEdge_frontier (D := d + 1) hb₃ hi₃ hs₃ (by rw [if_pos (hb₃.tree.2 ⟨_, _, _, rfl⟩)]))
        (wp_and hg₁.2 hb₁.2)) (invTree child (d + 1) _ pre hs₁.frame' hg₁.1 hb₁.1)
end

/-- The schedule frontier holds at every `finishEdge` of the walk of `t`, given the ear guards and
the bookkeeping facts (the hypotheses of `walkTree_inv'`). -/
theorem walkTree_frontiers (t : DfsTree) (d : Nat) (s : WalkState)
    (hi : ∀ v outs, t = .node v outs → ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).Inv' d)
    (hs : Shape s) (hg : GuardsTree t d s) (hb : BookTree t d s) : FrontiersTree t d s :=
  frTree t d s hi hs hg hb

end WalkState
end Spqr
