import Mathlib.Tactic.Push
import Spqr.SepPairExhaust
import Spqr.Proofs.SepPair
import Spqr.Proofs.Type2

/-!
# Every separation pair of a block is a type-1 or type-2 pair (PROOF.md §4.5)

Over a block (`TwoConnected`, `Rooted`) with its sorted DFS tree (`Spec`):
* `sepPair_of_type1`, `sepPair_of_type2`: soundness of the two splits;
* `type1_or_type2_of_sepPair`: exhaustiveness — with no type-1 child and no type-2 split at
  `{a, b}` all edges not joining `a` and `b` form one separation class, and at most one edge joins
  them, so `{a, b}` is not a separation pair;
* `sepPair_iff`, `sepPair_iff'`, `three_connected_of_no_split`.
-/

namespace Spqr

namespace Graph

variable {g : Graph}

theorem Reach.and_of_adj {ok P : Nat → Prop} {x y : Nat} (h : g.Reach ok x y) (hx : P x)
    (hP : ∀ y z, g.Adj y z → P z) : g.Reach (fun w => ok w ∧ P w) x y := by
  induction h with
  | refl hok => exact .refl ⟨hok, hx⟩
  | tail _ hadj hok ih => exact ih.tail hadj ⟨hok, hP _ _ hadj⟩

theorem sepClass_comm {a b e e' : Nat} : g.SepClass a b e e' ↔ g.SepClass b a e e' :=
  ⟨fun h => h.mono fun _ h => ⟨h.2, h.1⟩, fun h => h.mono fun _ h => ⟨h.2, h.1⟩⟩

theorem SeparationPair.symm {a b : Nat} (h : g.SeparationPair a b) : g.SeparationPair b a := by
  obtain ⟨hne, h⟩ := h
  have : g.SepClass a b = g.SepClass b a := by
    funext e e'
    exact propext sepClass_comm
  rw [this] at h
  exact ⟨Ne.symm hne, h⟩

/-- An edge joining `a` and `b` is its own separation class. -/
theorem sepClass_eq_of_joins {a b e e' : Nat} (hj : g.Joins e a b) (h : g.SepClass a b e e') :
    e = e' := by
  rcases h with h | ⟨x, y, hx, hy, hw⟩
  · exact h
  · rcases hx.eq_or hj with h | h
    · exact absurd h hw.ok_left.1
    · exact absurd h hw.ok_left.2

theorem SeparationPair.two_le {a b : Nat} (h : g.SeparationPair a b) : 2 ≤ g.ne := by
  obtain ⟨e, e', he, he', hs⟩ := h.exists_not_sepClass
  have : e ≠ e' := fun h => hs (by rw [h]; exact EdgeConn.refl _)
  omega

end Graph

namespace DfsData

variable {g : Graph} {d : DfsData}

/-- Vertices of the *above* side of `{a, b}`: outside `T_{a'} ∪ {a, b}`, or in a child subtree of
`b` returning strictly above `a`. -/
def AboveSide (d : DfsData) (a a' b x : Nat) : Prop :=
  x ≠ a ∧ x ≠ b ∧
    (¬d.Anc a' x ∨ ∃ c, d.IsParent b c ∧ d.Anc c x ∧ ∃ l, l < d.depth a ∧ d.Returns c l)

section Spec

variable (hs : d.Spec g)
include hs

/-- In a block the root has one child. -/
theorem root_one_child (h2 : g.TwoConnected) {c c' : Nat} (hc : d.IsParent d.root c)
    (hc' : d.IsParent d.root c') : c = c' := by
  by_contra hne
  obtain ⟨e, he⟩ := hc.joins hs
  obtain ⟨e', he'⟩ := hc'.joins hs
  rcases h2 d.root e e' he.lt he'.lt with rfl | ⟨x, y, hx, hy, hw⟩
  · rcases he.eq_or he' with ⟨-, h⟩ | ⟨h, h'⟩
    · exact hne h
    · exact hne (h'.trans h)
  · have hxc := (hx.eq_or he).resolve_left hw.ok_left
    have hyc := (hy.eq_or he').resolve_left hw.ok_right
    rw [hxc, hyc] at hw
    have hcc : ¬d.Anc c c' := fun h => by
      have h3 : c = d.root := root_anc_eq hs (h.to_parent hs hne hc')
      rw [h3] at hc
      exact not_isParent_root hs hc
    obtain ⟨u, v, -, -, -, hv, hvr⟩ := exit_subtree hs hc (.refl c) hcc hw fun h => h rfl
    exact hvr (root_anc_eq hs hv)

/-- A vertex strictly above `a` reaches the root avoiding `{a, b}`. -/
theorem reach_root_above {a b w : Nat} (hab : d.Anc a b) (hw : d.Anc d.root w)
    (hwa : d.depth w < d.depth a) : g.Reach (fun x => x ≠ a ∧ x ≠ b) w d.root :=
  reach_root_of_not_anc hs hw (fun h => h.ne_of_depth_lt hs hwa)
    (fun h => (hab.trans h).ne_of_depth_lt hs hwa)

/-- A subtree `T_c` disjoint from `{a, b}` returning strictly above `a` reaches the root. -/
theorem reach_root_of_returns (hr : d.Rooted g) {a b c x l : Nat} (hab : d.Anc a b)
    (hok : ∀ z, d.Anc c z → z ≠ a ∧ z ≠ b) (hx : d.Anc c x) (hret : d.Returns c l)
    (hl : l < d.depth a) : g.Reach (fun z => z ≠ a ∧ z ≠ b) x d.root := by
  obtain ⟨u, o, ho, hcu, hb, rfl⟩ := hret
  have hj := hs.joins u o ho
  have hru : d.Anc d.root u := hr _ _ hj.isEnd
  have hwu := hs.back_anc u o ho hb
  have hrw : d.Anc d.root o.dest := by
    rcases hwu.comparable hs hru with h | h
    · rw [root_anc_eq hs h]; exact .refl _
    · exact h
  have hdb := hab.depth_le hs
  have hwok : o.dest ≠ a ∧ o.dest ≠ b :=
    ⟨fun h => by rw [h] at hl; omega, fun h => by rw [h] at hl; omega⟩
  exact ((reach_in_subtree hs hx hcu hok).tail hj.adj hwok).trans (reach_root_above hs hab hrw hl)

/-- Vertices outside `T_{a'} ∪ {a, b}` reach the root avoiding `{a, b}` (`a` not the root). -/
theorem reach_root_outside (hr : d.Rooted g) (h2 : g.TwoConnected) {q a a' b x : Nat}
    (hqa : d.IsParent q a) (hpa : d.IsParent a a') (ha'b : d.Anc a' b) (hx : d.Anc d.root x)
    (hxa : x ≠ a) (hx' : ¬d.Anc a' x) :
    g.Reach (fun z => z ≠ a ∧ z ≠ b) x d.root := by
  have hab : d.Anc a b := hpa.anc.trans ha'b
  by_cases hax : d.Anc a x
  · obtain ⟨c, hac, hcx⟩ := hax.child (Ne.symm hxa)
    have hne : c ≠ a' := fun h => hx' (h ▸ hcx)
    obtain ⟨l, hl, hret⟩ := child_returns_above hs h2 hqa hac
    refine reach_root_of_returns hs hr hab (fun z hz => ⟨?_, ?_⟩) hcx hret hl
    · intro h
      have h1 := hs.depth_parent _ _ hac
      have h2 := hz.depth_le hs
      rw [h] at h2; omega
    · exact fun h => not_anc_of_sibling hs hpa hac hne ha'b hz (by rw [h]; exact .refl _)
  · exact reach_root_of_not_anc hs hx hax fun h => hax (hab.trans h)

/-- A vertex of depth strictly between `a` and `b` above a vertex of `T_b` lies in
`T_{a'} − T_b`. -/
theorem mem_between {a a' b w u : Nat} (hpa : d.IsParent a a') (ha'b : d.Anc a' b)
    (hwu : d.Anc w u) (hbu : d.Anc b u) (h1 : d.depth a < d.depth w) (h2 : d.depth w < d.depth b) :
    d.Anc a' w ∧ ¬d.Anc b w := by
  have hda' := hs.depth_parent _ _ hpa
  have hwb : d.Anc w b := (hwu.comparable hs hbu).resolve_right fun h => h.ne_of_depth_lt hs h2
  refine ⟨?_, fun h => h.ne_of_depth_lt hs h2⟩
  rcases ha'b.comparable hs hwb with h | h
  · exact h
  · have := h.depth_le hs
    rw [h.eq_of_depth_eq hs (Anc.refl a') (by omega)]
    exact .refl _

theorem aboveSide_closed {a a' b y z : Nat} (hpa : d.IsParent a a')
    (hno : d.NoBothSides a b) (hst : d.BetweenStays a a' b) (hy : d.AboveSide a a' b y)
    (hadj : g.Adj y z) (hza : z ≠ a) (hzb : z ≠ b) : d.AboveSide a a' b z := by
  refine ⟨hza, hzb, ?_⟩
  by_cases hz : d.Anc a' z
  case neg => exact .inl hz
  have hda' := hs.depth_parent _ _ hpa
  obtain ⟨hya, hyb, hy' | ⟨c, hbc, hcy, l, hl, hret⟩⟩ := hy
  · rcases adj_comparable hs hadj with hyz | hzy
    case inr => exact absurd (hz.trans hzy) hy'
    have hya' : d.Anc y a := by
      rcases hyz.comparable hs (hpa.anc.trans hz) with h | h
      · exact h
      · exfalso
        obtain ⟨c₀, hc₀, hc₀y⟩ := h.child (Ne.symm hya)
        have : c₀ = a' := (hc₀y.trans hyz).eq_of_depth_eq hs hz
          (by rw [hs.depth_parent _ _ hc₀, hda'])
        exact hy' (this ▸ hc₀y)
    have hdy : d.depth y < d.depth a := hya'.depth_lt hs hya
    obtain ⟨e, he⟩ := hadj
    have hyz' : y ≠ z := fun h => hy' (h ▸ hz)
    rcases edge_up hs he.symm hyz hyz' with hpar | ⟨o, ho, -, hback, hdest⟩
    · exfalso
      have := IsParent.eq_of_anc hs hz hpar hy'
      rw [this] at hpar
      exact hya (hs.parent_unique _ _ _ hpar hpa)
    · by_cases hbz : d.Anc b z
      · obtain ⟨c, hbc, hcz⟩ := hbz.child (Ne.symm hzb)
        exact .inr ⟨c, hbc, hcz, d.depth y, hdy, z, o, ho, hcz, hback, by rw [hdest]⟩
      · exfalso
        have := hst z hz hbz o ho hback
        rw [hdest] at this; omega
  · rcases adj_comparable hs hadj with hyz | hzy
    · exact .inr ⟨c, hbc, hcy.trans hyz, l, hl, hret⟩
    by_cases hcz : d.Anc c z
    · exact .inr ⟨c, hbc, hcz, l, hl, hret⟩
    exfalso
    have hzc : d.Anc z c := (hzy.comparable hs hcy).resolve_right hcz
    have hzc' : z ≠ c := fun h => hcz (by rw [h]; exact .refl _)
    have hzb' : d.Anc z b := hzc.to_parent hs hzc' hbc
    have hdz1 : d.depth a < d.depth z := by have := hz.depth_le hs; omega
    have hdz2 : d.depth z < d.depth b := hzb'.depth_lt hs hzb
    obtain ⟨e, he⟩ := hadj
    have hzy' : z ≠ y := fun h => hcz (by rw [h]; exact hcy)
    rcases edge_up hs he hzy hzy' with hpar | ⟨o, ho, -, hback, hdest⟩
    · have := IsParent.eq_of_anc hs hcy hpar hcz
      rw [this] at hpar
      exact hzb (hs.parent_unique _ _ _ hpar hbc)
    · exact hno c hbc ⟨l, hl, hret⟩ (d.depth z) hdz1 hdz2 ⟨y, o, ho, hcy, hback, by rw [hdest]⟩

theorem AboveSide.reach {a a' b x y : Nat} (hpa : d.IsParent a a')
    (hno : d.NoBothSides a b) (hst : d.BetweenStays a a' b) (hx : d.AboveSide a a' b x)
    (h : g.Reach (fun z => z ≠ a ∧ z ≠ b) x y) : d.AboveSide a a' b y := by
  induction h with
  | refl => exact hx
  | tail _ hadj hok ih => exact aboveSide_closed hs hpa hno hst ih hadj hok.1 hok.2

/-- Endpoints of an edge with an endpoint on the above side are on the above side. -/
theorem AboveSide.end {a a' b e x w : Nat} (hpa : d.IsParent a a') (hno : d.NoBothSides a b)
    (hst : d.BetweenStays a a' b) (hx : d.AboveSide a a' b x) (he : g.IsEnd e x)
    (hw : g.IsEnd e w) (hwa : w ≠ a) (hwb : w ≠ b) : d.AboveSide a a' b w := by
  obtain ⟨v, hv⟩ := he
  rcases hw.eq_or hv with h | h
  · exact h ▸ hx
  · exact h ▸ aboveSide_closed hs hpa hno hst hx ⟨e, hv⟩ (h ▸ hwa) (h ▸ hwb)

omit hs in
theorem not_aboveSide_between {a a' b y : Nat} (hy : d.Anc a' y) (hby : ¬d.Anc b y) :
    ¬d.AboveSide a a' b y := by
  rintro ⟨-, -, h | ⟨c, hbc, hcy, -⟩⟩
  · exact h hy
  · exact hby (hbc.anc.trans hcy)

/-- A type-2 pair satisfies the above/between separation hypothesis of Facts B/C
(`type2_first_out`). -/
theorem Type2Pair.above_between {a b : Nat} (h : d.Type2Pair a b) :
    ∃ q a', d.IsParent q a ∧ d.IsParent a a' ∧ d.Anc a' b ∧ a' ≠ b ∧
      ∀ e e', d.Above a' a b e g → d.Between a' b e' g → ¬g.SepClass a b e e' := by
  obtain ⟨q, a', hqa, hpa, ha'b, hne, hno, hst⟩ := h
  have hda' := hs.depth_parent _ _ hpa
  refine ⟨q, a', hqa, hpa, ha'b, hne, fun e e' ⟨x, hx, hxa, hxb, hx'⟩ ⟨y, hy, hy', hby⟩ hsep => ?_⟩
  have hxs : d.AboveSide a a' b x := ⟨hxa, hxb, .inl hx'⟩
  have hya : y ≠ a := fun h => by have := hy'.depth_le hs; rw [h] at this; omega
  have hyb : y ≠ b := fun h => hby (by rw [h]; exact .refl _)
  refine not_aboveSide_between (a := a) hy' hby ?_
  rcases hsep with rfl | ⟨x', y', hx', hy'', hw⟩
  · exact hxs.end hs hpa hno hst hx hy hya hyb
  · exact (((hxs.end hs hpa hno hst hx hx' hw.ok_left.1 hw.ok_left.2).reach hs hpa hno hst hw).end
      hs hpa hno hst hy'' hy hya hyb)

/-! ### Soundness -/

theorem sepPair_of_type2 (hr : d.Rooted g) (h2 : g.TwoConnected) {a b : Nat}
    (h : d.Type2Pair a b) : g.SeparationPair a b := by
  obtain ⟨q, a', hqa, hpa, ha'b, hne, hno, hst⟩ := h
  have hab : d.Anc a b := hpa.anc.trans ha'b
  have hdq := hs.depth_parent _ _ hqa
  have hda' := hs.depth_parent _ _ hpa
  have hdb := ha'b.depth_le hs
  obtain ⟨e₁, he₁⟩ := hqa.joins hs
  obtain ⟨l, hl, u, o, ho, hu, hback, hlo⟩ := child_returns_above hs h2 hqa hpa
  have hj := hs.joins u o ho
  have hdu := hu.depth_le hs
  obtain ⟨e₂, he₂⟩ := hpa.joins hs
  obtain ⟨x, hx, hxb⟩ := ha'b.child hne
  obtain ⟨e₂', he₂'⟩ := hx.joins hs
  have hdx := hs.depth_parent _ _ hx
  refine ⟨fun h => by rw [h] at hda'; omega,
    .inr ⟨e₁, o.e, e₂, e₂', he₁.lt, hj.lt, he₂.lt, he₂'.lt, ?_, ?_, ?_, ?_, ?_⟩⟩
  · rintro rfl
    rcases he₁.eq_or hj with ⟨-, h⟩ | ⟨-, h⟩
    · rw [← h] at hlo; omega
    · rw [← h] at hdu; omega
  · have hrq : d.Anc d.root q := hr _ _ he₁.isEnd
    have hrw : d.Anc d.root o.dest := hr _ _ hj.symm.isEnd
    exact .of_reach he₁.isEnd hj.symm.isEnd
      ((reach_root_above hs hab hrq (by omega)).trans
        (reach_root_above hs hab hrw (by omega)).symm)
  · rintro rfl
    rcases he₂.eq_or he₂' with ⟨h, -⟩ | ⟨h, -⟩
    · rw [h] at hda'; omega
    · rw [← h] at hdx; omega
  · exact .of_reach he₂.symm.isEnd he₂'.isEnd (.refl ⟨fun h => by rw [h] at hda'; omega, hne⟩)
  · rintro (h | ⟨y, z, hy, hz, hw⟩)
    · rw [h] at he₁
      rcases he₁.eq_or he₂ with ⟨h, -⟩ | ⟨h, -⟩
      · rw [h] at hdq; omega
      · rw [h] at hdq; omega
    · have hyq : y = q := (hy.eq_or he₁).resolve_right hw.ok_left.1
      have hza' : z = a' := (hz.eq_or he₂).resolve_left hw.ok_right.1
      rw [hyq, hza'] at hw
      have hq : d.AboveSide a a' b q :=
        ⟨fun h => by rw [h] at hdq; omega, fun h => by rw [h] at hdq; omega,
          .inl fun h => by have := h.depth_le hs; omega⟩
      obtain ⟨-, -, h | ⟨c, hbc, hca', -⟩⟩ := hq.reach hs hpa hno hst hw
      · exact h (.refl _)
      · exact hne (ha'b.antisymm hs (hbc.anc.trans hca'))

theorem sepPair_of_type1 {a b : Nat} (h : d.Type1Pair a b g) : g.SeparationPair a b := by
  obtain ⟨hab, hne, ⟨o, ho, hcls, e, e', he, he', hee', hce, hce'⟩ |
    ⟨e₁, e₂, e₃, he₃, h12, h13, h23, hj₁, hj₂⟩⟩ := h
  · obtain ⟨ht, -, hret, -, -⟩ := (hs.cls_type1 b o ho _).mp hcls
    have hj := hs.joins b o ho
    obtain ⟨u, o', ho', hu, hb', -⟩ := hret
    have hj' := hs.joins u o' ho'
    have hC : d.EndIn o.dest o.e g := ⟨o.dest, hj.symm.isEnd, .refl _⟩
    have hC' : d.EndIn o.dest o'.e g := ⟨u, hj'.isEnd, hu⟩
    have hne' : o.e ≠ o'.e := fun h => by
      obtain ⟨-, rfl⟩ := hs.out_inj _ _ o ho o' ho' h
      rw [ht] at hb'; cases hb'
    have hiff := fun {e'} => type1_class hs hab ho hcls hC (e' := e')
    refine ⟨hne, ?_⟩
    by_cases hsep : g.SepClass a b e e'
    · exact .inr ⟨o.e, o'.e, e, e', hj.lt, hj'.lt, he, he', hne', hiff.mpr hC', hee', hsep,
        fun h => hce (hiff.mp h)⟩
    · exact .inl ⟨o.e, e, e', hj.lt, he, he', fun h => hce (hiff.mp h), fun h => hce' (hiff.mp h),
        hsep⟩
  · exact ⟨hne, .inl ⟨e₁, e₂, e₃, hj₁.lt, hj₂.lt, he₃, fun h => h12 (Graph.sepClass_eq_of_joins hj₁ h),
      fun h => h13 (Graph.sepClass_eq_of_joins hj₁ h), fun h => h23 (Graph.sepClass_eq_of_joins hj₂ h)⟩⟩

/-! ### Exhaustiveness -/

/-- An edge not joining `a` and `b` has an endpoint outside `{a, b}`. -/
theorem exists_end_not_ab (h2 : g.TwoConnected) (hne : 2 ≤ g.ne) {a b e : Nat} (he : e < g.ne)
    (hab : ¬g.Joins e a b) : ∃ x, g.IsEnd e x ∧ x ≠ a ∧ x ≠ b := by
  obtain ⟨v, o, ho, rfl⟩ := hs.edge_out e he
  have hj := hs.joins v o ho
  have hvw : v ≠ o.dest := fun h => h2.not_loop hne (by rw [← h] at hj; exact hj)
  by_cases hva : v = a
  · subst hva
    exact ⟨o.dest, hj.symm.isEnd, Ne.symm hvw, fun h => hab (by rw [← h]; exact hj)⟩
  · by_cases hvb : v = b
    · subst hvb
      exact ⟨o.dest, hj.symm.isEnd, fun h => hab (by rw [← h]; exact hj.symm), Ne.symm hvw⟩
    · exact ⟨v, hj.isEnd, hva, hvb⟩

/-- `a` not the root, no type-1 child of `b` to `a`, no type-2 split: every vertex other than
`a`, `b` reaches the root avoiding `{a, b}`. -/
theorem reach_root_of_not_type2 (hr : d.Rooted g) (h2 : g.TwoConnected) {q p a a' b x : Nat}
    (hqa : d.IsParent q a) (hpa : d.IsParent a a') (ha'b : d.Anc a' b) (hpb : d.IsParent p b)
    (ht1 : ∀ o ∈ d.outs b, o.cls ≠ .ret (d.depth a) .type1Child) (hn2 : ¬d.Type2Pair a b)
    (hx : d.Anc d.root x) (hxa : x ≠ a) (hxb : x ≠ b) :
    g.Reach (fun z => z ≠ a ∧ z ≠ b) x d.root := by
  have hab : d.Anc a b := hpa.anc.trans ha'b
  have hda' := hs.depth_parent _ _ hpa
  have hbet : ∀ y, d.Anc a' y → ¬d.Anc b y → g.Reach (fun z => z ≠ a ∧ z ≠ b) y d.root := by
    intro y hy hby
    have hne : a' ≠ b := fun h => hby (h ▸ hy)
    by_cases hst : d.BetweenStays a a' b
    · have hno : ¬d.NoBothSides a b := fun hno => hn2 ⟨q, a', hqa, hpa, ha'b, hne, hno, hst⟩
      unfold NoBothSides at hno
      push Not at hno
      obtain ⟨c, hbc, ⟨l₁, hl₁, hret₁⟩, l₂, hl₂, hl₂', u, o, ho, hcu, hback, hlo⟩ := hno
      have hj := hs.joins u o ho
      obtain ⟨hwa', hbw⟩ := mem_between hs hpa ha'b (hs.back_anc u o ho hback)
        (hbc.anc.trans hcu) (by omega) (by omega)
      exact ((between_reach hs hpa hy hby hwa' hbw).tail hj.adj.symm (subtree_ok hs hab hbc hcu)).trans
        (reach_root_of_returns hs hr hab (fun z hz => subtree_ok hs hab hbc hz) hcu hret₁ hl₁)
    · unfold BetweenStays at hst
      push Not at hst
      obtain ⟨u, hu, hbu, o, ho, hback, hlt⟩ := hst
      have hj := hs.joins u o ho
      have hdb := hab.depth_le hs
      exact ((between_reach hs hpa hy hby hu hbu).tail hj.adj
        ⟨fun h => by rw [h] at hlt; omega, fun h => by rw [h] at hlt; omega⟩).trans
        (reach_root_above hs hab (hr _ _ hj.symm.isEnd) (Nat.lt_of_not_le hlt))
  by_cases hx' : d.Anc a' x
  case neg => exact reach_root_outside hs hr h2 hqa hpa ha'b hx hxa hx'
  by_cases hbx : d.Anc b x
  case neg => exact hbet x hx' hbx
  obtain ⟨c, ⟨o, ho, ht, rfl⟩, hcx⟩ := hbx.child (Ne.symm hxb)
  have hbc : d.IsParent b o.dest := isParent_of_tree ho ht
  have hok : ∀ z, d.Anc o.dest z → z ≠ a ∧ z ≠ b := fun z hz => subtree_ok hs hab hbc hz
  rcases tree_out_trichotomy hs h2 (d.depth a) hpb ho ht with hup | h1 | hdn
  · obtain ⟨l, hl, hret⟩ := AttachesAbove.returns hs ho ht hup
    exact reach_root_of_returns hs hr hab hok hcx hret hl
  · exact absurd (IsType1To.tree_cls hs ho ht h1) (ht1 o ho)
  · obtain ⟨l, hl, hl', u, o', ho', hcu, hback, hlo⟩ := AttachesBetween.returns hs ho ht hdn
    have hj' := hs.joins u o' ho'
    obtain ⟨hwa', hbw⟩ := mem_between hs hpa ha'b (hs.back_anc u o' ho' hback)
      (hbc.anc.trans hcu) (by omega) (by omega)
    exact ((reach_in_subtree hs hcx hcu hok).tail hj'.adj
      ⟨fun h => by rw [h] at hlo; omega, fun h => by rw [h] at hlo; omega⟩).trans (hbet _ hwa' hbw)

/-- `a` the root, no type-1 child of `b` to `a`: every vertex other than `a`, `b` reaches `a`'s
child `a'` avoiding `{a, b}`. -/
theorem reach_child_of_root (h2 : g.TwoConnected) {p a' b x : Nat}
    (hpa : d.IsParent d.root a') (ha'b : d.Anc a' b) (hpb : d.IsParent p b)
    (ht1 : ∀ o ∈ d.outs b, o.cls ≠ .ret (d.depth d.root) .type1Child) (hx : d.Anc d.root x)
    (hxa : x ≠ d.root) (hxb : x ≠ b) : g.Reach (fun z => z ≠ d.root ∧ z ≠ b) x a' := by
  have hab := hpa.anc.trans ha'b
  have hda' := hs.depth_parent _ _ hpa
  have hd0 := hs.depth_root
  obtain ⟨c₀, hc₀, hc₀x⟩ := hx.child (Ne.symm hxa)
  rw [root_one_child hs h2 hc₀ hpa] at hc₀x
  by_cases hbx : d.Anc b x
  case neg => exact between_reach hs hpa hc₀x hbx (.refl _) fun h => hbx (h.trans hc₀x)
  obtain ⟨c, ⟨o, ho, ht, rfl⟩, hcx⟩ := hbx.child (Ne.symm hxb)
  have hbc : d.IsParent b o.dest := isParent_of_tree ho ht
  have hok : ∀ z, d.Anc o.dest z → z ≠ d.root ∧ z ≠ b := fun z hz => subtree_ok hs hab hbc hz
  rcases tree_out_trichotomy hs h2 (d.depth d.root) hpb ho ht with hup | h1 | hdn
  · obtain ⟨l, hl, -⟩ := AttachesAbove.returns hs ho ht hup
    omega
  · exact absurd (IsType1To.tree_cls hs ho ht h1) (ht1 o ho)
  · obtain ⟨l, hl, hl', u, o', ho', hcu, hback, hlo⟩ := AttachesBetween.returns hs ho ht hdn
    have hj' := hs.joins u o' ho'
    obtain ⟨hwa', hbw⟩ := mem_between hs hpa ha'b (hs.back_anc u o' ho' hback)
      (hbc.anc.trans hcu) (by omega) (by omega)
    exact ((reach_in_subtree hs hcx hcu hok).tail hj'.adj
      ⟨fun h => by rw [h] at hlo; omega, fun h => by rw [h] at hlo; omega⟩).trans
      (between_reach hs hpa hwa' hbw (.refl _) fun h => hbw (h.trans hwa'))

/-- With no type-1 child of `b` to `a` and no type-2 split, all edges not joining `a` and `b`
are one separation class. -/
theorem sepClass_of_not_type2 (hr : d.Rooted g) (h2 : g.TwoConnected) {a b : Nat}
    (hab : d.Anc a b) (hne : a ≠ b) (hne2 : 2 ≤ g.ne)
    (ht1 : ∀ o ∈ d.outs b, o.cls ≠ .ret (d.depth a) .type1Child) (hn2 : ¬d.Type2Pair a b)
    {e e' : Nat} (he : e < g.ne) (he' : e' < g.ne) (hje : ¬g.Joins e a b)
    (hje' : ¬g.Joins e' a b) : g.SepClass a b e e' := by
  obtain ⟨a', hpa, ha'b⟩ := hab.child hne
  obtain ⟨p, hpb⟩ : ∃ p, d.IsParent p b := by
    rcases ha'b.cases_tail with h | ⟨p, -, hp⟩
    · exact ⟨a, by rw [h]; exact hpa⟩
    · exact ⟨p, hp⟩
  obtain ⟨x, hx, hxa, hxb⟩ := exists_end_not_ab hs h2 hne2 he hje
  obtain ⟨y, hy, hya, hyb⟩ := exists_end_not_ab hs h2 hne2 he' hje'
  have hrx := hr _ _ hx
  have hry := hr _ _ hy
  by_cases hroot : a = d.root
  · subst hroot
    exact .of_reach hx hy ((reach_child_of_root hs h2 hpa ha'b hpb ht1 hrx hxa hxb).trans
      (reach_child_of_root hs h2 hpa ha'b hpb ht1 hry hya hyb).symm)
  · obtain ⟨e₀, he₀⟩ := hpa.joins hs
    obtain ⟨q, hqa⟩ : ∃ q, d.IsParent q a := by
      rcases (hr _ _ he₀.isEnd).cases_tail with h | ⟨q, -, hq⟩
      · exact absurd h hroot
      · exact ⟨q, hq⟩
    exact .of_reach hx hy
      ((reach_root_of_not_type2 hs hr h2 hqa hpa ha'b hpb ht1 hn2 hrx hxa hxb).trans
        (reach_root_of_not_type2 hs hr h2 hqa hpa ha'b hpb ht1 hn2 hry hya hyb).symm)

theorem type1_or_type2_of_sepPair (hr : d.Rooted g) (h2 : g.TwoConnected) {a b : Nat}
    (hab : d.Anc a b) (h : g.SeparationPair a b) : d.Type1Pair a b g ∨ d.Type2Pair a b := by
  by_contra hn
  push Not at hn
  obtain ⟨hn1, hn2⟩ := hn
  have hne := h.1
  have hne2 := h.two_le
  by_cases ht1 : ∃ o ∈ d.outs b, o.cls = .ret (d.depth a) .type1Child
  · obtain ⟨o, ho, hcls⟩ := ht1
    have hC : ∀ e e', e < g.ne → e' < g.ne → e ≠ e' →
        d.EndIn o.dest e g ∨ d.EndIn o.dest e' g := by
      intro e e' he he' hee'
      by_contra hc
      push Not at hc
      exact hn1 ⟨hab, hne, .inl ⟨o, ho, hcls, e, e', he, he', hee', hc.1, hc.2⟩⟩
    have hcls' : ∀ {e e'}, d.EndIn o.dest e g → d.EndIn o.dest e' g → g.SepClass a b e e' :=
      fun he he' => (type1_class hs hab ho hcls he).mpr he'
    have hcls'' : ∀ {e e'}, g.SepClass a b e e' → d.EndIn o.dest e g → d.EndIn o.dest e' g :=
      fun hs' he => (type1_class hs hab ho hcls he).mp hs'
    rcases h.2 with ⟨e₁, e₂, e₃, h1, h2, h3, n12, n13, n23⟩ |
      ⟨e₁, e₁', e₂, e₂', h1, h1', h2, h2', d1, s1, d2, s2, n12⟩
    · have ne12 : e₁ ≠ e₂ := fun h => n12 (by rw [h]; exact Graph.EdgeConn.refl _)
      have ne13 : e₁ ≠ e₃ := fun h => n13 (by rw [h]; exact Graph.EdgeConn.refl _)
      have ne23 : e₂ ≠ e₃ := fun h => n23 (by rw [h]; exact Graph.EdgeConn.refl _)
      rcases hC e₁ e₂ h1 h2 ne12 with c1 | c2
      · rcases hC e₂ e₃ h2 h3 ne23 with c2 | c3
        · exact n12 (hcls' c1 c2)
        · exact n13 (hcls' c1 c3)
      · rcases hC e₁ e₃ h1 h3 ne13 with c1 | c3
        · exact n12 (hcls' c1 c2)
        · exact n23 (hcls' c2 c3)
    · have ne12 : e₁ ≠ e₂ := fun h => n12 (by rw [h]; exact Graph.EdgeConn.refl _)
      rcases hC e₁ e₂ h1 h2 ne12 with c1 | c2
      · rcases hC e₂ e₂' h2 h2' d2 with c2 | c2'
        · exact n12 (hcls' c1 c2)
        · exact n12 (hcls' c1 (hcls'' s2.symm c2'))
      · rcases hC e₁ e₁' h1 h1' d1 with c1 | c1'
        · exact n12 (hcls' c1 c2)
        · exact n12 (hcls' (hcls'' s1.symm c1') c2)
  · push Not at ht1
    by_cases hbond : ∃ e₁ e₂, e₁ ≠ e₂ ∧ g.Joins e₁ a b ∧ g.Joins e₂ a b
    · obtain ⟨e₁, e₂, h12, hj₁, hj₂⟩ := hbond
      have hall : ∀ e, e < g.ne → e = e₁ ∨ e = e₂ := by
        intro e he
        by_contra hc
        push Not at hc
        exact hn1 ⟨hab, hne, .inr ⟨e₁, e₂, e, he, h12, Ne.symm hc.1, Ne.symm hc.2, hj₁, hj₂⟩⟩
      have hsing : ∀ e e', e < g.ne → g.SepClass a b e e' → e = e' := by
        intro e e' he hs'
        rcases hall e he with rfl | rfl
        · exact Graph.sepClass_eq_of_joins hj₁ hs'
        · exact Graph.sepClass_eq_of_joins hj₂ hs'
      rcases h.2 with ⟨f₁, f₂, f₃, h1, h2, h3, n12, n13, n23⟩ |
        ⟨f₁, f₁', -, -, h1, -, -, -, d1, s1, -, -, -⟩
      · rcases hall f₁ h1 with rfl | rfl <;> rcases hall f₂ h2 with rfl | rfl <;>
          rcases hall f₃ h3 with rfl | rfl
        all_goals first
          | exact n12 (Graph.EdgeConn.refl _)
          | exact n13 (Graph.EdgeConn.refl _)
          | exact n23 (Graph.EdgeConn.refl _)
      · exact d1 (hsing f₁ f₁' h1 s1)
    · push Not at hbond
      have hone : ∀ {e e'}, e < g.ne → e' < g.ne → ¬g.Joins e a b → ¬g.Joins e' a b →
          g.SepClass a b e e' :=
        fun he he' hje hje' => sepClass_of_not_type2 hs hr h2 hab hne hne2 ht1 hn2 he he' hje hje'
      rcases h.2 with ⟨f₁, f₂, f₃, h1, h2, h3, n12, n13, n23⟩ |
        ⟨f₁, f₁', f₂, f₂', h1, -, h2, -, d1, s1, d2, s2, n12⟩
      · have ne12 : f₁ ≠ f₂ := fun h => n12 (by rw [h]; exact Graph.EdgeConn.refl _)
        have ne13 : f₁ ≠ f₃ := fun h => n13 (by rw [h]; exact Graph.EdgeConn.refl _)
        have ne23 : f₂ ≠ f₃ := fun h => n23 (by rw [h]; exact Graph.EdgeConn.refl _)
        by_cases j1 : g.Joins f₁ a b <;> by_cases j2 : g.Joins f₂ a b <;>
          by_cases j3 : g.Joins f₃ a b
        · exact hbond f₁ f₂ ne12 j1 j2
        · exact hbond f₁ f₂ ne12 j1 j2
        · exact hbond f₁ f₃ ne13 j1 j3
        · exact n23 (hone h2 h3 j2 j3)
        · exact hbond f₂ f₃ ne23 j2 j3
        · exact n13 (hone h1 h3 j1 j3)
        · exact n12 (hone h1 h2 j1 j2)
        · exact n12 (hone h1 h2 j1 j2)
      · have j1 : ¬g.Joins f₁ a b := fun j => d1 (Graph.sepClass_eq_of_joins j s1)
        have j2 : ¬g.Joins f₂ a b := fun j => d2 (Graph.sepClass_eq_of_joins j s2)
        exact n12 (hone h1 h2 j1 j2)

/-- Hopcroft–Tarjan: the separation pairs `{a, b}`, `a` an ancestor of `b`, of a block are exactly
the type-1 and type-2 pairs. -/
theorem sepPair_iff (hr : d.Rooted g) (h2 : g.TwoConnected) {a b : Nat} (hab : d.Anc a b) :
    g.SeparationPair a b ↔ d.Type1Pair a b g ∨ d.Type2Pair a b :=
  ⟨type1_or_type2_of_sepPair hs hr h2 hab,
    fun h => h.elim (sepPair_of_type1 hs) (sepPair_of_type2 hs hr h2)⟩

/-- `sepPair_iff` for two vertices of the tree in either order (Fact A). -/
theorem sepPair_iff' (hr : d.Rooted g) (h2 : g.TwoConnected) {a b : Nat} (ha : d.Anc d.root a)
    (hb : d.Anc d.root b) :
    g.SeparationPair a b ↔
      (d.Type1Pair a b g ∨ d.Type2Pair a b) ∨ (d.Type1Pair b a g ∨ d.Type2Pair b a) := by
  constructor
  · intro h
    rcases sepPair_comparable hs hr h2 ha hb h with hab | hba
    · exact .inl ((sepPair_iff hs hr h2 hab).mp h)
    · exact .inr ((sepPair_iff hs hr h2 hba).mp h.symm)
  · rintro (h | h)
    · exact h.elim (sepPair_of_type1 hs) (sepPair_of_type2 hs hr h2)
    · exact (h.elim (sepPair_of_type1 hs) (sepPair_of_type2 hs hr h2)).symm

omit hs in
/-- Removing a non-vertex `a` changes nothing: in a block every two edges stay connected. -/
theorem sepClass_of_not_vertex (hr : d.Rooted g) (h2 : g.TwoConnected) {a b e e' : Nat}
    (ha : ¬d.Anc d.root a) (he : e < g.ne) (he' : e' < g.ne) : g.SepClass a b e e' := by
  rcases h2 b e e' he he' with h | ⟨x, y, hx, hy, hw⟩
  · exact .inl h
  · refine .inr ⟨x, y, hx, hy, ?_⟩
    have hxa : x ≠ a := fun h => ha (by rw [← h]; exact hr _ _ hx)
    exact (hw.and_of_adj hxa fun y z ⟨e, hj⟩ h => ha (by rw [← h]; exact hr _ _ hj.symm.isEnd)).mono
      fun _ h => ⟨h.2, h.1⟩

/-- A block in which no type-1 or type-2 split applies is 3-connected. -/
theorem three_connected_of_no_split (hr : d.Rooted g) (h2 : g.TwoConnected)
    (h : ∀ a b, ¬d.Type1Pair a b g ∧ ¬d.Type2Pair a b) : g.ThreeConnected := by
  intro a b hab
  obtain ⟨e, e', he, he', hsep⟩ := hab.exists_not_sepClass
  by_cases ha : d.Anc d.root a
  · by_cases hb : d.Anc d.root b
    · rcases (sepPair_iff' hs hr h2 ha hb).mp hab with h' | h'
      · exact h'.elim (h a b).1 (h a b).2
      · exact h'.elim (h b a).1 (h b a).2
    · exact hsep (Graph.sepClass_comm.mpr (sepClass_of_not_vertex hr h2 hb he he'))
  · exact hsep (sepClass_of_not_vertex hr h2 ha he he')

end Spec

end DfsData

end Spqr
