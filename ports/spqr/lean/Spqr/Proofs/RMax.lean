import Spqr.RMax
import Spqr.Proofs.Contract
import Spqr.Proofs.Blocks

/-!
# R-maximality (PROOF.md §4.5, row `spqrTree_r_three_connected`)

From the walk-side facts `RCloseShape` and the separation-pair layer (`sepPair_iff'`,
`threeConnected_contract_iff`):
* `RCloseShape.not_sepPair`: no skeleton pair of `U` is a separation pair of the block — a
  type-1 or type-2 class at such a pair is 2-attached at the pair, touches both its vertices and
  is proper, and no such class is laminar with `U` (`not_laminarWith`);
* `RCloseShape.threeConnected`: the R skeleton `(P.addParent g U s t).contract g` is 3-connected.
-/

namespace Spqr

namespace Graph

variable {g : Graph} {E K : Nat → Prop}

theorem adjIn_of_joins {f u w : Nat} (hj : g.Joins f u w) (hE : E f) : g.AdjIn E u w := by
  refine ⟨f, hj.lt, hE, ?_⟩
  rcases (joins_iff.1 hj).2 with h | h
  · exact .inl h.symm
  · exact .inr (by rw [h])

/-- A vertex of `E` with an edge outside `E` is a terminal. -/
theorem TwoAttached.boundary {a b v e : Nat} (h : g.TwoAttached E a b) (hv : g.Touches E v)
    (he : g.IsEnd e v) (hE : ¬E e) : v = a ∨ v = b := by
  obtain ⟨e₁, he₁, hE₁, hv₁⟩ := hv
  exact h v e₁ e he₁ (isEnd_iff.1 he).1 hE₁ hE hv₁ (isEnd_iff.1 he).2

/-- A non-terminal vertex of `E` is interior. -/
theorem TwoAttached.interior {a b v : Nat} (h : g.TwoAttached E a b) (hv : g.Touches E v)
    (ha : v ≠ a) (hb : v ≠ b) : g.Interior E v := fun e he hv' => by
  by_contra hE
  rcases h.boundary hv (isEnd_iff.2 ⟨he, hv'⟩) hE with h' | h'
  · exact ha h'
  · exact hb h'

/-- A separation class of `{a, b}` is 2-attached at `{a, b}`. -/
theorem sepClass_twoAttached {a b e₀ : Nat} : g.TwoAttached (g.SepClass a b e₀) a b := by
  intro v e e' he he' hE hE' hv hv'
  by_contra hn
  exact hE' (hE.trans (.of_reach (isEnd_iff.2 ⟨he, hv⟩) (isEnd_iff.2 ⟨he', hv'⟩)
    (.refl ⟨fun h => hn (.inl h), fun h => hn (.inr h)⟩)))

/-- In a block, a nonempty proper edge set 2-attached at `{a, b}` has an outside edge at `b`. -/
theorem TwoAttached.exists_not_mem_end (h2 : g.TwoConnected) {a b : Nat}
    (h : g.TwoAttached K a b) (hne : ∃ e, e < g.ne ∧ K e) (hprop : ∃ e, e < g.ne ∧ ¬K e) :
    ∃ f, f < g.ne ∧ ¬K f ∧ g.IsEnd f b := by
  obtain ⟨eK, heK, hK⟩ := hne
  obtain ⟨eo, heo, ho⟩ := hprop
  rcases h2 a eo eK heo heK with rfl | ⟨x, y, hx, hy, hr⟩
  · exact absurd hK ho
  by_cases hxK : g.Touches K x
  · rcases h.boundary hxK hx ho with h' | h'
    · exact absurd h' hr.ok_left
    · exact ⟨eo, heo, ho, h' ▸ hx⟩
  rcases hr.exit (ok' := fun v => ¬g.Touches K v) hxK with hin | ⟨u, v, hu, ⟨f, hf⟩, hva, hvK⟩
  · exact absurd (Touches.of_isEnd hK hy) hin.ok_right.2
  · have hnf : ¬K f := fun hKf => hu.ok_right.2 (Touches.of_isEnd hKf hf.isEnd)
    rcases h.boundary (not_not.1 hvK) hf.symm.isEnd hnf with h' | h'
    · exact absurd h' hva
    · exact ⟨f, hf.lt, hnf, h' ▸ hf.symm.isEnd⟩

/-- The complement of a 2-attached edge set of a block is connected. -/
theorem TwoAttached.conn_compl (h2 : g.TwoConnected) {s t : Nat} (h : g.TwoAttached E s t) :
    g.ConnEdges fun e => ¬E e := by
  intro e e' he he' hE hE'
  have toFst : ∀ {e y}, e < g.ne → ¬E e → g.IsEnd e y →
      Relation.ReflTransGen (g.AdjIn fun e => ¬E e) (g.edges[e]!).1 y := by
    intro e y he hE ⟨w, hw⟩
    rcases (joins_iff.1 hw).2 with h' | h'
    · rw [h']
    · exact .single ⟨e, he, hE, .inl (by rw [h'])⟩
  rcases h2 t e e' he he' with rfl | ⟨y, y', hy, hy', hr⟩
  · exact .refl
  have key : ∀ {z}, g.Reach (· ≠ t) y z →
      (Relation.ReflTransGen (g.AdjIn fun e => ¬E e) y z ∧ ∃ f, f < g.ne ∧ ¬E f ∧ g.IsEnd f z) ∨
        (Relation.ReflTransGen (g.AdjIn fun e => ¬E e) y s ∧ g.Touches E z) := by
    intro z hz
    induction hz with
    | refl => exact .inl ⟨.refl, e, he, hE, hy⟩
    | @tail z z' hr₀ hadj hok ih =>
      obtain ⟨f, hf⟩ := hadj
      rcases ih with ⟨hR, f₀, -, hE₀, hz⟩ | ⟨hR, hz⟩
      · by_cases hEf : E f
        · rcases h.boundary (Touches.of_isEnd hEf hf.isEnd) hz hE₀ with h' | h'
          · exact .inr ⟨h' ▸ hR, Touches.of_isEnd hEf hf.symm.isEnd⟩
          · exact absurd h' hr₀.ok_right
        · exact .inl ⟨hR.tail (adjIn_of_joins hf hEf), f, hf.lt, hEf, hf.symm.isEnd⟩
      · by_cases hEf : E f
        · exact .inr ⟨hR, Touches.of_isEnd hEf hf.symm.isEnd⟩
        · rcases h.boundary hz hf.isEnd hEf with h' | h'
          · exact .inl ⟨hR.tail (h' ▸ adjIn_of_joins hf hEf), f, hf.lt, hEf, hf.symm.isEnd⟩
          · exact absurd h' hr₀.ok_right
  have hyy' : Relation.ReflTransGen (g.AdjIn fun e => ¬E e) y y' := by
    rcases key hr with ⟨hR, -⟩ | ⟨hR, hz⟩
    · exact hR
    · rcases h.boundary hz hy' hE' with h' | h'
      · exact h' ▸ hR
      · exact absurd h' hr.ok_right
  exact ((toFst he hE hy).trans hyy').trans (Graph.reach_symm (toFst he' hE' hy'))

/-- An isolated vertex is in no separation pair. -/
theorem not_sepPair_of_isolated (h2 : g.TwoConnected) {a b : Nat} (ha : ∀ e, ¬g.IsEnd e a) :
    ¬g.SeparationPair a b := by
  intro h
  obtain ⟨e, e', he, he', hn⟩ := h.exists_not_sepClass
  apply hn
  rcases h2 b e e' he he' with h' | ⟨x, y, hx, hy, hr⟩
  · exact .inl h'
  · refine .inr ⟨x, y, hx, hy, ?_⟩
    exact (hr.and_of_adj (P := (· ≠ a)) (fun h' => ha e (h' ▸ hx))
      (fun _ z ⟨f, hj⟩ h' => ha f (h' ▸ hj.symm.isEnd))).mono fun _ h => ⟨h.2, h.1⟩

end Graph

namespace Pieces

variable {g : Graph} {P : Pieces} {U K : Nat → Prop} {s t : Nat}

theorem SkelPair.symm {a b : Nat} (h : P.SkelPair g U s t a b) : P.SkelPair g U s t b a := by
  obtain ⟨ha, hb, hsa, hsb, hnt, hst⟩ := h
  refine ⟨hb, ha, hsb, hsa, fun ⟨i, hi, h⟩ => hnt ⟨i, hi, ?_⟩, fun h => hst ?_⟩
  · rcases h with ⟨h1, h2⟩ | ⟨h1, h2⟩
    · exact .inr ⟨h2, h1⟩
    · exact .inl ⟨h2, h1⟩
  · rcases h with ⟨h1, h2⟩ | ⟨h1, h2⟩
    · exact .inr ⟨h2, h1⟩
    · exact .inl ⟨h2, h1⟩

/-- No class 2-attached at a skeleton pair `{a, b}` of `U`, touching `a` and `b` and proper, is
laminar with `U`. -/
theorem not_laminarWith (h2 : g.TwoConnected) (hU : g.TwoAttached U s t)
    {a b : Nat} (hab : P.SkelPair g U s t a b) (hne : a ≠ b) (hK : g.TwoAttached K a b)
    (hKa : g.Touches K a) (hKb : g.Touches K b) (hprop : ∃ e, e < g.ne ∧ ¬K e)
    (hl : P.LaminarWith U K) : False := by
  obtain ⟨hUa, hUb, hsa, hsb, hnt, hst⟩ := hab
  obtain ⟨eK, heK, hKe, hvK⟩ := hKa
  rcases hl with ⟨i, hi, hsub⟩ | hdisj | hsup
  · have ha := terminal_of_skel hsa hi ⟨eK, heK, hsub _ hKe, hvK⟩
    obtain ⟨eb, heb, hKb', hvb⟩ := hKb
    have hb := terminal_of_skel hsb hi ⟨eb, heb, hsub _ hKb', hvb⟩
    apply hnt ⟨i, hi, ?_⟩
    rcases ha with ha | ha <;> rcases hb with hb | hb
    · exact absurd (ha.trans hb.symm) hne
    · exact .inl ⟨ha, hb⟩
    · exact .inr ⟨ha, hb⟩
    · exact absurd (ha.trans hb.symm) hne
  · have ha := hU.boundary hUa (Graph.isEnd_iff.2 ⟨heK, hvK⟩) (hdisj _ hKe)
    obtain ⟨eb, heb, hKb', hvb⟩ := hKb
    have hb := hU.boundary hUb (Graph.isEnd_iff.2 ⟨heb, hvb⟩) (hdisj _ hKb')
    apply hst
    rcases ha with ha | ha <;> rcases hb with hb | hb
    · exact absurd (ha.trans hb.symm) hne
    · exact .inl ⟨ha, hb⟩
    · exact .inr ⟨ha, hb⟩
    · exact absurd (ha.trans hb.symm) hne
  · have hKne : ∃ e, e < g.ne ∧ K e := ⟨eK, heK, hKe⟩
    have int : ∀ w, g.Touches U w → w ≠ s → w ≠ t → g.Interior K w := fun w hw hs ht e he hv =>
      hsup e (hU.interior hw hs ht e he hv)
    by_cases has : a = s
    · have hbt : b ≠ t := fun hbt => hst (.inl ⟨has, hbt⟩)
      have hbs : b ≠ s := fun hbs => hne (has.trans hbs.symm)
      obtain ⟨f, hf, hnf, hfb⟩ := hK.exists_not_mem_end h2 hKne hprop
      exact hnf (int b hUb hbs hbt f hf (Graph.isEnd_iff.1 hfb).2)
    by_cases hat : a = t
    · have hbs : b ≠ s := fun hbs => hst (.inr ⟨hat, hbs⟩)
      have hbt : b ≠ t := fun hbt => hne (hat.trans hbt.symm)
      obtain ⟨f, hf, hnf, hfb⟩ := hK.exists_not_mem_end h2 hKne hprop
      exact hnf (int b hUb hbs hbt f hf (Graph.isEnd_iff.1 hfb).2)
    obtain ⟨f, hf, hnf, hfa⟩ := hK.comm.exists_not_mem_end h2 hKne hprop
    exact hnf (int a hUa has hat f hf (Graph.isEnd_iff.1 hfa).2)

/-! ### Terminal pairs of a skeleton -/

/-- A terminal pair whose piece's complement is one class of the pair does not separate the
skeleton. -/
theorem WF.not_sepPair_terminal (hP : P.WF g) (h2 : g.TwoConnected) {i : Nat} (hk : i < P.k)
    (hmax : ∀ e e', e < g.ne → e' < g.ne → ¬P.Mem i e → ¬P.Mem i e' →
      g.SepClass (P.x i) (P.y i) e e') :
    ¬(P.contract g).SeparationPair (P.x i) (P.y i) := by
  obtain ⟨fi, hfi⟩ := orig_inr_iff.2 hk
  have hji : (P.contract g).Joins fi (P.x i) (P.y i) := hfi.joins_inr.2 (.inl ⟨rfl, rfl⟩)
  have key : ∀ f f', f < (P.contract g).ne → f' < (P.contract g).ne → f ≠ fi → f' ≠ fi →
      (P.contract g).SepClass (P.x i) (P.y i) f f' := by
    intro f f' hf hf' hne hne'
    obtain ⟨e, he, hef⟩ := hP.pre_total hf
    obtain ⟨e', he', hef'⟩ := hP.pre_total hf'
    have hm : ¬P.Mem i e := fun hm => hne ((hef.orig_of_mem hm).inj hfi)
    have hm' : ¬P.Mem i e' := fun hm => hne' ((hef'.orig_of_mem hm).inj hfi)
    refine hP.edgeConn_contract h2 (fun v hv => ?_) hef hef' (hmax e e' he he' hm hm')
    by_cases hx : v = P.x i
    · exact hx ▸ hP.skel_x hk
    · have hy : v = P.y i := by by_contra hy; exact hv ⟨hx, hy⟩
      exact hy ▸ hP.skel_y hk
  rintro ⟨-, ⟨f₁, f₂, f₃, h₁, h₂, h₃, h₁₂, h₁₃, h₂₃⟩ |
    ⟨f₁, f₁', f₂, f₂', h₁, h₁', h₂, h₂', hn₁, hc₁, hn₂, hc₂, hn⟩⟩
  · by_cases hf₁ : f₁ = fi
    · subst hf₁
      have hf₂ : f₂ ≠ f₁ := fun h => h₁₂ (h ▸ .inl rfl)
      have hf₃ : f₃ ≠ f₁ := fun h => h₁₃ (h ▸ .inl rfl)
      exact h₂₃ (key f₂ f₃ h₂ h₃ hf₂ hf₃)
    by_cases hf₂ : f₂ = fi
    · subst hf₂
      have hf₃ : f₃ ≠ f₂ := fun h => h₂₃ (h ▸ .inl rfl)
      exact h₁₃ (key f₁ f₃ h₁ h₃ hf₁ hf₃)
    exact h₁₂ (key f₁ f₂ h₁ h₂ hf₁ hf₂)
  · have hf₁ : f₁ ≠ fi := fun h => by subst h; exact hn₁ (Graph.sepClass_eq_of_joins hji hc₁)
    have hf₂ : f₂ ≠ fi := fun h => by subst h; exact hn₂ (Graph.sepClass_eq_of_joins hji hc₂)
    exact hn (key f₁ f₂ h₁ h₂ hf₁ hf₂)

/-! ### `addParent` -/

section AddParent

variable (hsub : ∀ i e, P.Mem i e → U e)
include hsub

theorem addParent_mem_iff {i e : Nat} (hi : i < P.k) :
    (P.addParent g U s t).Mem i e ↔ P.Mem i e := by
  simp only [Mem, addParent]
  split_ifs with hU he
  · exact Iff.rfl
  · exact ⟨fun h => absurd (Option.some.inj h).symm (Nat.ne_of_lt hi),
      fun h => absurd (hsub i e h) hU⟩
  · exact ⟨fun h => h.elim, fun h => absurd (hsub i e h) hU⟩

omit hsub in
theorem addParent_mem_top_iff (hP : P.WF g) {e : Nat} :
    (P.addParent g U s t).Mem P.k e ↔ e < g.ne ∧ ¬U e := by
  simp only [Mem, addParent]
  split_ifs with hU he
  · exact ⟨fun h => absurd (hP.lt _ _ h).1 (Nat.lt_irrefl _), fun h => absurd hU h.2⟩
  · exact ⟨fun _ => ⟨he, hU⟩, fun _ => rfl⟩
  · exact ⟨fun h => h.elim, fun h => absurd h.1 he⟩

omit hsub in
theorem addParent_x_of_lt {i : Nat} (hi : i < P.k) : (P.addParent g U s t).x i = P.x i := by
  simp [addParent, Nat.ne_of_lt hi]

omit hsub in
theorem addParent_y_of_lt {i : Nat} (hi : i < P.k) : (P.addParent g U s t).y i = P.y i := by
  simp [addParent, Nat.ne_of_lt hi]

omit hsub in
theorem addParent_x_top : (P.addParent g U s t).x P.k = s := by simp [addParent]

omit hsub in
theorem addParent_y_top : (P.addParent g U s t).y P.k = t := by simp [addParent]

theorem addParent_touches_iff {i w : Nat} (hi : i < P.k) :
    g.Touches ((P.addParent g U s t).Mem i) w ↔ g.Touches (P.Mem i) w := by
  constructor <;> rintro ⟨e, he, hE, hv⟩
  · exact ⟨e, he, (addParent_mem_iff hsub hi).1 hE, hv⟩
  · exact ⟨e, he, (addParent_mem_iff hsub hi).2 hE, hv⟩

theorem addParent_int_iff {i w : Nat} (hi : i < P.k) :
    (P.addParent g U s t).Int g i w ↔ P.Int g i w := by
  simp only [Int, addParent_touches_iff hsub hi, addParent_x_of_lt hi, addParent_y_of_lt hi]
  exact and_congr_left' ⟨fun _ => hi, fun _ => Nat.lt_succ_of_lt hi⟩

theorem addParent_wf (hP : P.WF g) (h2 : g.TwoConnected) (hU : g.TwoAttached U s t)
    (hts : g.Touches U s) (hst : s ≠ t) (hprop : ∃ e, e < g.ne ∧ ¬U e) :
    (P.addParent g U s t).WF g := by
  have top : ∀ e, e < g.ne → ((P.addParent g U s t).Mem P.k e ↔ ¬U e) := fun e he => by
    rw [addParent_mem_top_iff hP]; exact ⟨fun h => h.2, fun h => ⟨he, h⟩⟩
  have lt_or : ∀ i, i < P.k + 1 → i < P.k ∨ i = P.k := fun i hi => by omega
  refine ⟨fun i e h => ?_, fun i hi => ?_, fun i hi => ?_, fun i hi => ?_, fun i hi => ?_,
    fun i hi => ?_⟩
  · simp only [Mem, addParent] at h
    split_ifs at h with hU he
    · exact ⟨Nat.lt_succ_of_lt (hP.lt i e h).1, (hP.lt i e h).2⟩
    · exact ⟨by rw [← Option.some.inj h]; exact Nat.lt_succ_self _, he⟩
  · rcases lt_or i hi with hi | rfl
    · exact (Graph.ConnEdges.congr fun e _ => addParent_mem_iff hsub hi).2 (hP.conn i hi)
    · exact (Graph.ConnEdges.congr top).2 (hU.conn_compl h2)
  · rcases lt_or i hi with hi | rfl
    · rw [addParent_x_of_lt hi, addParent_y_of_lt hi]
      exact (Graph.TwoAttached.congr fun e _ => addParent_mem_iff hsub hi).2 (hP.attached i hi)
    · rw [addParent_x_top, addParent_y_top]
      refine (Graph.TwoAttached.congr top).2 fun v e e' he he' hE hE' hv hv' => ?_
      exact hU v e' e he' he (not_not.1 hE') hE hv' hv
  · rcases lt_or i hi with hi | rfl
    · rw [addParent_x_of_lt hi, addParent_y_of_lt hi]
      exact ⟨(addParent_touches_iff hsub hi).2 (hP.touch i hi).1,
        (addParent_touches_iff hsub hi).2 (hP.touch i hi).2⟩
    · rw [addParent_x_top, addParent_y_top]
      obtain ⟨e, he, hE⟩ := hts.isEnd
      obtain ⟨f, hf, hnf, hfs⟩ := hU.comm.exists_not_mem_end h2 ⟨e, hE.lt, he⟩ hprop
      obtain ⟨f', hf', hnf', hft⟩ := hU.exists_not_mem_end h2 ⟨e, hE.lt, he⟩ hprop
      exact ⟨⟨f, hf, (top f hf).2 hnf, (Graph.isEnd_iff.1 hfs).2⟩,
        ⟨f', hf', (top f' hf').2 hnf', (Graph.isEnd_iff.1 hft).2⟩⟩
  · rcases lt_or i hi with hi | rfl
    · rw [addParent_x_of_lt hi, addParent_y_of_lt hi]; exact hP.ne i hi
    · rw [addParent_x_top, addParent_y_top]; exact hst
  · rcases lt_or i hi with hi | rfl
    · obtain ⟨e, he, hn⟩ := hP.proper i hi
      exact ⟨e, he, fun h => hn ((addParent_mem_iff hsub hi).1 h)⟩
    · obtain ⟨e, he, hE, -⟩ := hts
      exact ⟨e, he, fun h => (top e he).1 h hE⟩

end AddParent

end Pieces

/-! ### The DFS side -/

namespace DfsData

variable {g : Graph} {d : DfsData} {P : Pieces} {U : Nat → Prop} {s t : Nat}

section Spec

variable (hs : d.Spec g)
include hs

theorem IsParent.joins' {p c : Nat} (h : d.IsParent p c) : ∃ o ∈ d.outs p, g.Joins o.e p c := by
  obtain ⟨o, ho, -, rfl⟩ := h
  exact ⟨o, ho, hs.joins p o ho⟩

end Spec

end DfsData

section Main

open DfsData

variable {g : Graph} {d : DfsData} {P : Pieces} {U : Nat → Prop} {s t : Nat}
variable (hs : d.Spec g)
include hs

/-- A type-1 pair at a skeleton pair of `U` contradicts `RCloseShape.type1` / `bond`. -/
theorem RCloseShape.not_type1 (h : RCloseShape g d P U s t) (h2 : g.TwoConnected) {a b : Nat}
    (hab : P.SkelPair g U s t a b) (ht : d.Type1Pair a b g) : False := by
  obtain ⟨hanc, hne, ⟨o, ho, hcls, e, e', he, he', hee', hne₁, hne₂⟩ | ⟨e₁, e₂, e₃, -, h₁₂, -, -, hj₁, hj₂⟩⟩ := ht
  · obtain ⟨htree, -, ⟨u, o', ho', hcu, hback, hdep⟩, -, -⟩ := (hs.cls_type1 b o ho _).mp hcls
    have hbc : d.IsParent b o.dest := isParent_of_tree ho htree
    have hjo := hs.joins b o ho
    have hKb : g.Touches (d.EndIn o.dest · g) b :=
      ⟨o.e, hjo.lt, ⟨o.dest, hjo.symm.isEnd, .refl _⟩, (Graph.isEnd_iff.1 hjo.isEnd).2⟩
    have hua : o'.dest = a :=
      (hs.back_anc u o' ho' hback).eq_of_depth_eq hs (hanc.trans (hbc.anc.trans hcu)) hdep
    have hjo' := hs.joins u o' ho'
    have hKa : g.Touches (d.EndIn o.dest · g) a :=
      ⟨o'.e, hjo'.lt, ⟨u, hjo'.isEnd, hcu⟩, (Graph.isEnd_iff.1 (hua ▸ hjo'.symm.isEnd)).2⟩
    have hK : g.TwoAttached (d.EndIn o.dest · g) a b :=
      (Graph.TwoAttached.congr fun e _ =>
        type1_class hs hanc ho hcls (e := o.e) ⟨o.dest, hjo.symm.isEnd, .refl _⟩).1
        Graph.sepClass_twoAttached
    exact Pieces.not_laminarWith h2 h.attached hab hne hK hKa hKb ⟨e, he, hne₁⟩
      (h.type1 a b hab hanc o ho hcls)
  · obtain ⟨hUa, hUb, hsa, hsb, hnt, hst⟩ := hab
    by_cases hU : U e₁ ∧ U e₂
    · obtain ⟨i, hi₁, hi₂⟩ := h.bond a b e₁ e₂ h₁₂ hj₁ hj₂ hU.1 hU.2
      have hk := (h.wf.lt i e₁ hi₁).1
      have ha := Pieces.terminal_of_skel hsa hk (Graph.Touches.of_isEnd hi₁ hj₁.isEnd)
      have hb := Pieces.terminal_of_skel hsb hk (Graph.Touches.of_isEnd hi₁ hj₁.symm.isEnd)
      apply hnt ⟨i, hk, ?_⟩
      rcases ha with ha | ha <;> rcases hb with hb | hb
      · exact absurd (ha.trans hb.symm) hne
      · exact .inl ⟨ha, hb⟩
      · exact .inr ⟨ha, hb⟩
      · exact absurd (ha.trans hb.symm) hne
    · have : ∃ e, g.Joins e a b ∧ ¬U e := by
        by_cases hU₁ : U e₁
        · exact ⟨e₂, hj₂, fun hU₂ => hU ⟨hU₁, hU₂⟩⟩
        · exact ⟨e₁, hj₁, hU₁⟩
      obtain ⟨e, hj, hnU⟩ := this
      have ha := h.attached.boundary hUa hj.isEnd hnU
      have hb := h.attached.boundary hUb hj.symm.isEnd hnU
      apply hst
      rcases ha with ha | ha <;> rcases hb with hb | hb
      · exact absurd (ha.trans hb.symm) hne
      · exact .inl ⟨ha, hb⟩
      · exact .inr ⟨ha, hb⟩
      · exact absurd (ha.trans hb.symm) hne

/-- A type-2 pair at a skeleton pair of `U` contradicts `RCloseShape.type2`. -/
theorem RCloseShape.not_type2 (h : RCloseShape g d P U s t) (h2 : g.TwoConnected) {a b : Nat}
    (hab : P.SkelPair g U s t a b) (ht : d.Type2Pair a b) : False := by
  obtain ⟨q, a', hqa, hpa, ha'b, hne, hno, hst⟩ := id ht
  obtain ⟨q₂, a₂, -, hpa₂, ha₂b, -, hsep⟩ := ht.above_between hs
  have ha₂ : a₂ = a' := ha₂b.eq_of_depth_eq hs ha'b
    (by rw [hs.depth_parent _ _ hpa₂, hs.depth_parent _ _ hpa])
  rw [ha₂] at hsep
  obtain ⟨o₀, ho₀, htree₀, hdest₀⟩ := id hpa
  have hj₀ : g.Joins o₀.e a a' := hdest₀ ▸ hs.joins a o₀ ho₀
  have hda' := hs.depth_parent _ _ hpa
  have hdq := hs.depth_parent _ _ hqa
  have hdb := ha'b.depth_le hs
  have hab' : a ≠ b := fun h => by rw [h] at hda'; omega
  have hBet₀ : d.Between a' b o₀.e g :=
    ⟨a', hj₀.symm.isEnd, .refl _, fun h => hne (ha'b.antisymm hs h)⟩
  have hK : g.TwoAttached (g.SepClass a b o₀.e) a b := Graph.sepClass_twoAttached
  have hKa : g.Touches (g.SepClass a b o₀.e) a :=
    Graph.Touches.of_isEnd (Graph.EdgeConn.refl _) hj₀.isEnd
  obtain ⟨p, hpb⟩ : ∃ p, d.Anc a' p ∧ d.IsParent p b := by
    rcases ha'b.cases_tail with h | h
    · exact absurd h.symm hne
    · exact h
  obtain ⟨op, hop, hjp⟩ := hpb.2.joins' hs
  have hKb : g.Touches (g.SepClass a b o₀.e) b :=
    Graph.Touches.of_isEnd
      (between_edge_conn hs hpa ha'b hne hj₀.symm.isEnd ⟨p, hjp.isEnd, hpb.1, hpb.2.not_anc hs⟩)
      hjp.symm.isEnd
  obtain ⟨oq, hoq, hjq⟩ := hqa.joins' hs
  have hprop : ∃ e, e < g.ne ∧ ¬g.SepClass a b o₀.e e := by
    refine ⟨oq.e, hjq.lt, fun hc => hsep oq.e o₀.e ⟨q, hjq.isEnd, hqa.ne hs, ?_, ?_⟩ hBet₀ hc.symm⟩
    · intro h; rw [h] at hdq; omega
    · intro h; have := h.depth_le hs; omega
  exact Pieces.not_laminarWith h2 h.attached hab hab' hK hKa hKb hprop
    (h.type2 a b hab ht o₀ ho₀ htree₀ (hdest₀ ▸ ha'b))

/-- No skeleton pair of `U` is a separation pair of the block. -/
theorem RCloseShape.not_sepPair (h : RCloseShape g d P U s t) (hr : d.Rooted g)
    (h2 : g.TwoConnected) {a b : Nat} (hab : P.SkelPair g U s t a b) :
    ¬g.SeparationPair a b := by
  intro hsep
  obtain ⟨ea, hEa, -, hva⟩ := hab.1
  obtain ⟨eb, hEb, -, hvb⟩ := hab.2.1
  have hra : d.Anc d.root a := hr ea a (Graph.isEnd_iff.2 ⟨hEa, hva⟩)
  have hrb : d.Anc d.root b := hr eb b (Graph.isEnd_iff.2 ⟨hEb, hvb⟩)
  rcases (d.sepPair_iff' hs hr h2 hra hrb).1 hsep with h' | h' <;> rcases h' with h' | h'
  · exact h.not_type1 hs h2 hab h'
  · exact h.not_type2 hs h2 hab h'
  · exact h.not_type1 hs h2 hab.symm h'
  · exact h.not_type2 hs h2 hab.symm h'

/-- The R skeleton is 3-connected. -/
theorem RCloseShape.threeConnected (h : RCloseShape g d P U s t) (hr : d.Rooted g)
    (h2 : g.TwoConnected) : ((P.addParent g U s t).contract g).ThreeConnected := by
  have hP' := Pieces.addParent_wf h.sub h.wf h2 h.attached h.touch_s h.ne h.proper
  rw [hP'.threeConnected_contract_iff h2]
  have top : ∀ e, e < g.ne → ((P.addParent g U s t).Mem P.k e ↔ ¬U e) := fun e he => by
    rw [Pieces.addParent_mem_top_iff h.wf]; exact ⟨fun h => h.2, fun h => ⟨he, h⟩⟩
  have touchU : ∀ {w}, (P.addParent g U s t).Skel g w → (∃ e, g.IsEnd e w) → g.Touches U w := by
    rintro w hw ⟨e, he⟩
    by_cases hU : U e
    · exact Graph.Touches.of_isEnd hU he
    · by_cases hws : w = s
      · exact hws ▸ h.touch_s
      by_cases hwt : w = t
      · exact hwt ▸ h.touch_t
      exact absurd ⟨Nat.lt_succ_self _, Graph.Touches.of_isEnd ((top e he.lt).2 hU) he,
        by rw [Pieces.addParent_x_top]; exact hws, by rw [Pieces.addParent_y_top]; exact hwt⟩
        (hw P.k)
  refine ⟨fun a b ha hb hnt hsep => ?_, fun i hi => ?_⟩
  · by_cases hea : ∃ e, g.IsEnd e a
    · by_cases heb : ∃ e, g.IsEnd e b
      · refine h.not_sepPair hs hr h2 ⟨touchU ha hea, touchU hb heb, fun i hi => ?_,
          fun i hi => ?_, fun hnt' => ?_, fun hst => ?_⟩ hsep
        · exact ha i ((Pieces.addParent_int_iff h.sub hi.1).2 hi)
        · exact hb i ((Pieces.addParent_int_iff h.sub hi.1).2 hi)
        · obtain ⟨i, hi, hi'⟩ := hnt'
          exact hnt ⟨i, Nat.lt_succ_of_lt hi, by
            rwa [Pieces.addParent_x_of_lt hi, Pieces.addParent_y_of_lt hi]⟩
        · exact hnt ⟨P.k, Nat.lt_succ_self _, by
            rwa [Pieces.addParent_x_top, Pieces.addParent_y_top]⟩
      · exact Graph.not_sepPair_of_isolated h2 (fun e he => heb ⟨e, he⟩) hsep.symm
    · exact Graph.not_sepPair_of_isolated h2 (fun e he => hea ⟨e, he⟩) hsep
  · have hi' : i < P.k + 1 := hi
    rcases (by omega : i < P.k ∨ i = P.k) with hi | rfl
    · refine hP'.not_sepPair_terminal h2 (Nat.lt_succ_of_lt hi) fun e e' he he' hn hn' => ?_
      rw [Pieces.addParent_x_of_lt hi, Pieces.addParent_y_of_lt hi]
      exact h.maximal i hi e e' he he' (fun m => hn ((Pieces.addParent_mem_iff h.sub hi).2 m))
        (fun m => hn' ((Pieces.addParent_mem_iff h.sub hi).2 m))
    · refine hP'.not_sepPair_terminal h2 (Nat.lt_succ_self _) fun e e' he he' hn hn' => ?_
      rw [Pieces.addParent_x_top, Pieces.addParent_y_top]
      exact h.single e e' he he' (not_not.1 fun hU => hn ((top e he).2 hU))
        (not_not.1 fun hU => hn' ((top e' he').2 hU))

end Main

end Spqr
