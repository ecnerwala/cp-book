import Spqr.Proofs.SepPair

/-!
# The class of `a → a'` for a type-2 pair

For `a = anc b l` with `a` not the root and the *above*/*between* parts separated by `{a, b}`, the
separation class of the tree edge `a → a'` toward `b` consists of: that edge, the out-edges of
`T_{a'} − T_b`, and the blocks of the out-edges of `b` sorted after `ret l backEdge`.
-/

namespace Spqr

namespace Graph

variable {g : Graph} {ok : Nat → Prop}

theorem IsEnd.eq_or {e u w y : Nat} (hy : g.IsEnd e y) (hj : g.Joins e u w) : y = u ∨ y = w := by
  obtain ⟨z, hz⟩ := hy
  rcases hz.eq_or hj with ⟨h, -⟩ | ⟨h, -⟩
  · exact .inl h
  · exact .inr h

theorem EdgeConn.trans {e e' e'' : Nat} (h : g.EdgeConn ok e e') (h' : g.EdgeConn ok e' e'') :
    g.EdgeConn ok e e'' := by
  rcases h with rfl | ⟨x, y, hx, hy, hr⟩
  · exact h'
  rcases h' with rfl | ⟨y', z, hy', hz, hr'⟩
  · exact .inr ⟨x, y, hx, hy, hr⟩
  · refine .inr ⟨x, z, hx, hz, hr.trans ?_⟩
    obtain ⟨y₁, hj⟩ := hy
    rcases hy'.eq_or hj with rfl | rfl
    · exact hr'
    · exact (Reach.tail (.refl hr.ok_right) hj.adj hr'.ok_left).trans hr'

theorem SepClass.exists_ok_end {a b e e' : Nat} (h : g.SepClass a b e e') (hne : e ≠ e') :
    ∃ y, g.IsEnd e' y ∧ y ≠ a ∧ y ≠ b := by
  rcases h with h | ⟨x, y, -, hy, hr⟩
  · exact absurd h hne
  · exact ⟨y, hy, hr.ok_right⟩

end Graph

namespace OutClass

theorem AttachesAbove.rank_lt {l : Nat} {c : OutClass} (h : c.AttachesAbove l) :
    c.rank < (OutClass.ret l .backEdge).rank := by
  obtain ⟨l', k, rfl, hl⟩ := h.ret
  have := k.rank_le
  show 3 * (l' + 2) + k.rank < 3 * (l + 2) + 1
  omega

theorem AttachesBetween.rank_gt {l : Nat} {c : OutClass} (h : c.AttachesBetween l) :
    (OutClass.ret l .backEdge).rank < c.rank := by
  obtain ⟨l', k, rfl, hl, hk⟩ := h.ret
  show 3 * (l + 2) + 1 < 3 * (l' + 2) + k.rank
  rcases Nat.lt_or_eq_of_le hl with hlt | rfl
  · omega
  · rw [hk rfl]; show 3 * (l + 2) + 1 < 3 * (l + 2) + 2; omega

end OutClass

namespace DfsData

variable {g : Graph} {d : DfsData}

section Spec

variable (hs : d.Spec g)
include hs

theorem IsType1To.tree_cls {v l : Nat} {o : DfsOut} (ho : o ∈ d.outs v) (ht : o.isTree = true)
    (h : o.cls.IsType1To l) : o.cls = .ret l .type1Child := by
  cases hc : o.cls with
  | ret l' k =>
    rw [hc] at h
    cases k with
    | type1Child => rw [show l' = l from h]
    | backEdge =>
      have := ((hs.cls_backEdge v o ho l').mp hc).1
      rw [ht] at this; cases this
    | type2Child => exact h.elim
  | bridge => rw [hc] at h; exact h.elim
  | component => rw [hc] at h; exact h.elim
  | selfLoop => rw [hc] at h; exact h.elim

theorem subtree_edge_conn {a b c e e' : Nat} (hab : d.Anc a b) (hp : d.IsParent b c)
    (he : d.EndIn c e g) (he' : d.EndIn c e' g) : g.SepClass a b e e' := by
  obtain ⟨x, hx, hcx⟩ := he
  obtain ⟨y, hy, hcy⟩ := he'
  exact .of_reach hx hy (reach_in_subtree hs hcx hcy fun w hw => subtree_ok hs hab hp hw)

theorem between_edge_conn {a a' b e₀ e : Nat} (hpa : d.IsParent a a') (ha'b : d.Anc a' b)
    (hne : a' ≠ b) (he₀ : g.IsEnd e₀ a') (he : d.Between a' b e g) : g.SepClass a b e₀ e := by
  obtain ⟨w, hw, haw, hbw⟩ := he
  exact .of_reach he₀ hw
    (between_reach hs hpa (Anc.refl a') (fun h => hne (ha'b.antisymm hs h)) haw hbw)

/-- A child of `b` whose subtree meets the class of `a → a'` is sorted after `ret (depth a)
backEdge`. -/
theorem type2_child_rank (h2 : g.TwoConnected) {a a' b : Nat} (hpa : d.IsParent a a')
    (ha'b : d.Anc a' b) (hne : a' ≠ b)
    (hsep : ∀ e e', d.Above a' a b e g → d.Between a' b e' g → ¬g.SepClass a b e e')
    {e₀ : Nat} (hj₀ : g.Joins e₀ a a') {o : DfsOut} (ho : o ∈ d.outs b) (ht : o.isTree = true)
    {e' : Nat} (he : d.EndIn o.dest e' g) (h : g.SepClass a b e₀ e') :
    (OutClass.ret (d.depth a) .backEdge).rank < o.cls.rank := by
  have hab : d.Anc a b := hpa.anc.trans ha'b
  have hda' := hs.depth_parent _ _ hpa
  have hda'b := ha'b.depth_le hs
  have hp : d.IsParent b o.dest := isParent_of_tree ho ht
  have hdc := hs.depth_parent _ _ hp
  obtain ⟨r, -, hrb⟩ := hab.cases_tail.resolve_left fun h => by subst h; omega
  rcases tree_out_trichotomy hs h2 (d.depth a) hrb ho ht with hA | hT | hB
  · obtain ⟨l', hl', hret⟩ := AttachesAbove.returns hs ho ht hA
    obtain ⟨eA, heA, hAb⟩ := attaches_above hs hab hpa hret hl'
    exact absurd (h.trans (subtree_edge_conn hs hab hp he heA)).symm
      (hsep eA e₀ hAb ⟨a', hj₀.symm.isEnd, Anc.refl _, fun h => hne (ha'b.antisymm hs h)⟩)
  · have hcls := IsType1To.tree_cls hs ho ht hT
    obtain ⟨y, hy, hcy⟩ := (type1_class hs hab ho hcls he).mp h.symm
    have := hcy.depth_le hs
    rcases hy.eq_or hj₀ with rfl | rfl <;> omega
  · exact hB.rank_gt

/-- Fact D, abstract form: the separation class of `a → a'` for a type-2 pair `{a, b}`. -/
theorem type2_class_iff (h2 : g.TwoConnected) {a a' b : Nat} (hpa : d.IsParent a a')
    (ha'b : d.Anc a' b) (hne : a' ≠ b)
    (hsep : ∀ e e', d.Above a' a b e g → d.Between a' b e' g → ¬g.SepClass a b e e')
    {o₀ : DfsOut} (ho₀ : o₀ ∈ d.outs a) (ht₀ : o₀.isTree = true) (hd₀ : o₀.dest = a') {e' : Nat} :
    g.SepClass a b o₀.e e' ↔
      e' = o₀.e ∨ (∃ v, ∃ o ∈ d.outs v, o.e = e' ∧ d.Anc a' v ∧ ¬d.Anc b v) ∨
        ∃ o ∈ d.outs b, (OutClass.ret (d.depth a) .backEdge).rank < o.cls.rank ∧
          ((o.isTree = true ∧ d.EndIn o.dest e' g) ∨ (o.isTree = false ∧ o.e = e')) := by
  have hab : d.Anc a b := hpa.anc.trans ha'b
  have hda' := hs.depth_parent _ _ hpa
  have hda'b := ha'b.depth_le hs
  have hj₀ : g.Joins o₀.e a a' := hd₀ ▸ hs.joins a o₀ ho₀
  have hba' : ¬d.Anc b a' := fun h => hne (ha'b.antisymm hs h)
  have hBet₀ : d.Between a' b o₀.e g := ⟨a', hj₀.symm.isEnd, Anc.refl _, hba'⟩
  have hnot : ∀ e, d.Above a' a b e g → ¬g.SepClass a b o₀.e e :=
    fun e hA h => hsep e o₀.e hA hBet₀ h.symm
  have hAboveDepth : ∀ e w, g.IsEnd e w → d.depth w < d.depth a → d.Above a' a b e g :=
    fun e w he hw =>
      ⟨w, he, fun h => by rw [h] at hw; omega, fun h => by rw [h] at hw; omega,
        fun h => h.ne_of_depth_lt hs (by omega)⟩
  obtain ⟨r, -, hrb⟩ := hab.cases_tail.resolve_left fun h => by subst h; omega
  constructor
  · intro h
    by_cases he₀ : e' = o₀.e
    · exact .inl he₀
    obtain ⟨y, hy, hya, hyb⟩ := h.exists_ok_end (Ne.symm he₀)
    obtain ⟨z, hjy⟩ := hy
    obtain ⟨v, o, ho, rfl, hvo⟩ := edge_out_cases hs hjy
    have hj := hs.joins v o ho
    have hends : ∀ w, g.IsEnd o.e w → w = v ∨ w = o.dest := fun w hw => hw.eq_or hj
    by_cases hva : v = a
    · subst hva
      exfalso
      cases ht : o.isTree
      · have hdest : d.Anc o.dest v := hs.back_anc v o ho ht
        by_cases hdv : o.dest = v
        · rw [hdv] at hvo
          rcases hvo with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ <;> exact hya rfl
        · exact hnot _ (hAboveDepth _ _ hj.symm.isEnd (hdest.depth_lt hs hdv)) h
      · have hp : d.IsParent v o.dest := isParent_of_tree ho ht
        have hdd := hs.depth_parent _ _ hp
        have hda' : o.dest ≠ a' := fun hd =>
          he₀ (congrArg DfsOut.e (hs.tree_inj v o ho o₀ ho₀ ht ht₀ (hd.trans hd₀.symm)))
        refine hnot _ ⟨o.dest, hj.symm.isEnd, fun h => by rw [h] at hdd; omega,
          fun h => by rw [h] at hdd; omega, fun h => ?_⟩ h
        have := h.depth_lt hs (Ne.symm hda')
        omega
    by_cases hv : d.Anc a' v
    · by_cases hbv : d.Anc b v
      · by_cases hvb : v = b
        · subst hvb
          by_cases hrank : (OutClass.ret (d.depth a) .backEdge).rank < o.cls.rank
          · refine .inr (.inr ⟨o, ho, hrank, ?_⟩)
            cases ht : o.isTree
            · exact .inr ⟨rfl, rfl⟩
            · exact .inl ⟨rfl, o.dest, hj.symm.isEnd, Anc.refl _⟩
          exfalso
          cases ht : o.isTree
          · cases hc : o.cls with
            | bridge =>
              have := ((hs.cls_bridge _ o ho).mp hc).1
              rw [ht] at this; cases this
            | component =>
              have := ((hs.cls_component _ o ho).mp hc).1
              rw [ht] at this; cases this
            | selfLoop =>
              obtain ⟨-, hdest⟩ := (hs.cls_selfLoop _ o ho).mp hc
              rcases hends y hjy.isEnd with rfl | rfl
              · exact hyb rfl
              · exact hyb hdest
            | ret l₀ k₀ =>
              have hb : o.cls = .ret l₀ .backEdge := by
                cases k₀ with
                | backEdge => exact hc
                | type1Child => exact absurd ((hs.cls_type1 _ o ho l₀).mp hc).1 (by rw [ht]; decide)
                | type2Child => exact absurd ((hs.cls_type2 _ o ho l₀).mp hc).1 (by rw [ht]; decide)
              obtain ⟨-, hdne, hdep⟩ := (hs.cls_backEdge _ o ho l₀).mp hb
              have hdest : d.Anc o.dest v := hs.back_anc v o ho ht
              rw [hb] at hrank
              have hr : ¬3 * (d.depth a + 2) + 1 < 3 * (l₀ + 2) + 1 := hrank
              rcases Nat.lt_or_eq_of_le (show l₀ ≤ d.depth a by omega) with hlt | heq
              · exact hnot _ (hAboveDepth _ _ hj.symm.isEnd (by omega)) h
              · have : o.dest = a := hdest.eq_of_depth_eq hs hab (by omega)
                rcases hends y hjy.isEnd with rfl | rfl
                · exact hyb rfl
                · exact hya this
          · exact hrank (type2_child_rank hs h2 hpa ha'b hne hsep hj₀ ho ht
              ⟨o.dest, hj.symm.isEnd, Anc.refl _⟩ h)
        · obtain ⟨c, hbc, hcv⟩ := hbv.child (Ne.symm hvb)
          obtain ⟨oc, hoc, htc, hdc⟩ := hbc
          have he : d.EndIn oc.dest o.e g := ⟨v, hj.isEnd, by rw [hdc]; exact hcv⟩
          exact .inr (.inr ⟨oc, hoc, type2_child_rank hs h2 hpa ha'b hne hsep hj₀ hoc htc he h,
            .inl ⟨htc, he⟩⟩)
      · exact .inr (.inl ⟨v, o, ho, rfl, hv, hbv⟩)
    · exact absurd h (hnot _ ⟨v, hj.isEnd, hva, fun h => hv (h ▸ ha'b), hv⟩)
  · rintro (rfl | ⟨v, o, ho, rfl, hv, hbv⟩ | ⟨o, ho, hrank, ⟨ht, he⟩ | ⟨ht, rfl⟩⟩)
    · exact .refl _
    · exact between_edge_conn hs hpa ha'b hne hj₀.symm.isEnd ⟨v, (hs.joins v o ho).isEnd, hv, hbv⟩
    · rcases tree_out_trichotomy hs h2 (d.depth a) hrb ho ht with hA | hT | hB
      · exact absurd hrank (Nat.not_lt.mpr (Nat.le_of_lt (OutClass.AttachesAbove.rank_lt hA)))
      · rw [IsType1To.tree_cls hs ho ht hT] at hrank
        exact absurd hrank (by show ¬3 * (d.depth a + 2) + 1 < 3 * (d.depth a + 2) + 0; omega)
      · obtain ⟨l', hl, hl', hret⟩ := AttachesBetween.returns hs ho ht hB
        have hp : d.IsParent b o.dest := isParent_of_tree ho ht
        obtain ⟨e'', he'', hBet⟩ := attaches_between hs hab hpa ha'b hp hret hl hl'
        exact (between_edge_conn hs hpa ha'b hne hj₀.symm.isEnd hBet).trans
          (subtree_edge_conn hs hab hp he'' he)
    · have hj := hs.joins b o ho
      cases hc : o.cls with
      | bridge =>
        have := ((hs.cls_bridge _ o ho).mp hc).1
        rw [ht] at this; cases this
      | component =>
        have := ((hs.cls_component _ o ho).mp hc).1
        rw [ht] at this; cases this
      | selfLoop =>
        rw [hc] at hrank
        exact absurd hrank (by show ¬3 * (d.depth a + 2) + 1 < 4; omega)
      | ret l₀ k₀ =>
        have hb : o.cls = .ret l₀ .backEdge := by
          cases k₀ with
          | backEdge => exact hc
          | type1Child => exact absurd ((hs.cls_type1 _ o ho l₀).mp hc).1 (by rw [ht]; decide)
          | type2Child => exact absurd ((hs.cls_type2 _ o ho l₀).mp hc).1 (by rw [ht]; decide)
        obtain ⟨-, hdne, hdep⟩ := (hs.cls_backEdge _ o ho l₀).mp hb
        have hdest : d.Anc o.dest b := hs.back_anc b o ho ht
        rw [hb] at hrank
        have hr : 3 * (d.depth a + 2) + 1 < 3 * (l₀ + 2) + 1 := hrank
        have ha'd : d.Anc a' o.dest := by
          rcases ha'b.comparable hs hdest with h | h
          · exact h
          · have := h.depth_le hs
            rw [h.eq_of_depth_eq hs (Anc.refl a') (by omega)]
            exact Anc.refl _
        exact between_edge_conn hs hpa ha'b hne hj₀.symm.isEnd
          ⟨o.dest, hj.symm.isEnd, ha'd, fun h => hdne (hdest.antisymm hs h)⟩

end Spec

end DfsData

end Spqr
