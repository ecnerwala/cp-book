import Spqr.SepPair
import Mathlib.Tactic.ByContra

/-!
# Separation pairs: Facts A–C of PROOF.md §3
-/

namespace Spqr

namespace Graph

variable {g : Graph} {ok ok' : Nat → Prop}

theorem Joins.symm {e x y : Nat} (h : g.Joins e x y) : g.Joins e y x := Or.symm h

theorem Joins.lt {e x y : Nat} (h : g.Joins e x y) : e < g.ne := by
  rcases h with h | h <;> exact Array.getElem?_eq_some_iff.mp h |>.1

theorem Joins.eq_or {e x y x' y' : Nat} (h : g.Joins e x y) (h' : g.Joins e x' y') :
    (x = x' ∧ y = y') ∨ (x = y' ∧ y = x') := by
  rcases h with h | h <;> rcases h' with h' | h' <;> simp_all

theorem Joins.isEnd {e x y : Nat} (h : g.Joins e x y) : g.IsEnd e x := ⟨y, h⟩

theorem Joins.adj {e x y : Nat} (h : g.Joins e x y) : g.Adj x y := ⟨e, h⟩

theorem Adj.symm {x y : Nat} (h : g.Adj x y) : g.Adj y x := h.imp fun _ h => h.symm

theorem IsEnd.lt {e x : Nat} (h : g.IsEnd e x) : e < g.ne := h.choose_spec.lt

theorem Reach.ok_left {x y : Nat} (h : g.Reach ok x y) : ok x := by
  induction h with
  | refl h => exact h
  | tail _ _ _ ih => exact ih

theorem Reach.ok_right {x y : Nat} (h : g.Reach ok x y) : ok y := by
  cases h with
  | refl h => exact h
  | tail _ _ h => exact h

theorem Reach.trans {x y z : Nat} (h : g.Reach ok x y) (h' : g.Reach ok y z) : g.Reach ok x z := by
  induction h' with
  | refl _ => exact h
  | tail _ hadj hok ih => exact ih.tail hadj hok

theorem Reach.head {x y z : Nat} (hx : ok x) (hadj : g.Adj x y) (h : g.Reach ok y z) :
    g.Reach ok x z :=
  (Reach.tail (.refl hx) hadj h.ok_left).trans h

theorem Reach.symm {x y : Nat} (h : g.Reach ok x y) : g.Reach ok y x := by
  induction h with
  | refl h => exact .refl h
  | tail _ hadj hok ih => exact ih.head hok hadj.symm

theorem Reach.mono (hok : ∀ x, ok x → ok' x) {x y : Nat} (h : g.Reach ok x y) : g.Reach ok' x y := by
  induction h with
  | refl h => exact .refl (hok _ h)
  | tail _ hadj hok' ih => exact ih.tail hadj (hok _ hok')

/-- A walk from inside `ok'` to outside it either stays inside or has a first exit edge. -/
theorem Reach.exit {x y : Nat} (h : g.Reach ok x y) (hx : ok' x) :
    g.Reach (fun v => ok v ∧ ok' v) x y ∨
      ∃ u v, g.Reach (fun v => ok v ∧ ok' v) x u ∧ g.Adj u v ∧ ok v ∧ ¬ok' v := by
  induction h with
  | refl h => exact .inl (.refl ⟨h, hx⟩)
  | @tail u v _ hadj hok ih =>
    rcases ih with ih | ⟨u', v', h1, h2, h3, h4⟩
    · by_cases hv : ok' v
      · exact .inl (ih.tail hadj ⟨hok, hv⟩)
      · exact .inr ⟨_, _, ih, hadj, hok, hv⟩
    · exact .inr ⟨u', v', h1, h2, h3, h4⟩

theorem EdgeConn.refl (e : Nat) : g.EdgeConn ok e e := .inl rfl

theorem EdgeConn.symm {e e' : Nat} (h : g.EdgeConn ok e e') : g.EdgeConn ok e' e := by
  rcases h with h | ⟨x, y, hx, hy, h⟩
  · exact .inl h.symm
  · exact .inr ⟨y, x, hy, hx, h.symm⟩

theorem EdgeConn.of_reach {e e' x y : Nat} (hx : g.IsEnd e x) (hy : g.IsEnd e' y)
    (h : g.Reach ok x y) : g.EdgeConn ok e e' := .inr ⟨x, y, hx, hy, h⟩

theorem EdgeConn.mono (hok : ∀ x, ok x → ok' x) {e e' : Nat} (h : g.EdgeConn ok e e') :
    g.EdgeConn ok' e e' := by
  rcases h with h | ⟨x, y, hx, hy, h⟩
  · exact .inl h
  · exact .inr ⟨x, y, hx, hy, h.mono hok⟩

theorem SeparationPair.exists_not_sepClass {a b : Nat} (h : g.SeparationPair a b) :
    ∃ e e', e < g.ne ∧ e' < g.ne ∧ ¬g.SepClass a b e e' := by
  rcases h with ⟨_, ⟨e₁, e₂, _, h₁, h₂, _, h, _⟩ | ⟨e₁, _, e₂, _, h₁, _, h₂, _, _, _, _, _, h⟩⟩
  · exact ⟨e₁, e₂, h₁, h₂, h⟩
  · exact ⟨e₁, e₂, h₁, h₂, h⟩

/-- A block with two edges has no loops. -/
theorem TwoConnected.not_loop (h2 : g.TwoConnected) (hne : 2 ≤ g.ne) {e x : Nat}
    (h : g.Joins e x x) : False := by
  obtain ⟨e', he', hne'⟩ : ∃ e', e' < g.ne ∧ e' ≠ e :=
    if h0 : e = 0 then ⟨1, by omega, by omega⟩ else ⟨0, by omega, Ne.symm h0⟩
  rcases h2 x e e' h.lt he' with h' | ⟨y, _, ⟨z, hy⟩, _, hr⟩
  · exact hne' h'.symm
  · obtain ⟨h1, _⟩ | ⟨_, h1⟩ := h.eq_or hy <;> exact hr.ok_left h1.symm

end Graph

namespace DfsData

variable {g : Graph} {d : DfsData}

theorem Anc.refl (x : Nat) : d.Anc x x := Relation.ReflTransGen.refl

theorem Anc.trans {x y z : Nat} (h : d.Anc x y) (h' : d.Anc y z) : d.Anc x z :=
  Relation.ReflTransGen.trans h h'

theorem IsParent.anc {p c : Nat} (h : d.IsParent p c) : d.Anc p c := Relation.ReflTransGen.single h

theorem Anc.tail {x y z : Nat} (h : d.Anc x y) (h' : d.IsParent y z) : d.Anc x z :=
  Relation.ReflTransGen.tail h h'

theorem Anc.head {x y z : Nat} (h : d.IsParent x y) (h' : d.Anc y z) : d.Anc x z :=
  Relation.ReflTransGen.head h h'

theorem Anc.cases_tail {x y : Nat} (h : d.Anc x y) : y = x ∨ ∃ p, d.Anc x p ∧ d.IsParent p y :=
  Relation.ReflTransGen.cases_tail h

/-- The child of `a` toward a proper descendant `x`. -/
theorem Anc.child {a x : Nat} (h : d.Anc a x) (hne : a ≠ x) : ∃ c, d.IsParent a c ∧ d.Anc c x := by
  induction h using Relation.ReflTransGen.head_induction_on with
  | refl => exact absurd rfl hne
  | head hp h _ => exact ⟨_, hp, h⟩

section Spec

variable (hs : d.Spec g)
include hs

theorem IsParent.joins {p c : Nat} (h : d.IsParent p c) : ∃ e, g.Joins e p c := by
  obtain ⟨o, ho, _, rfl⟩ := h
  exact ⟨_, hs.joins p o ho⟩

theorem IsParent.adj {p c : Nat} (h : d.IsParent p c) : g.Adj p c := h.joins hs

theorem Anc.depth_le {x y : Nat} (h : d.Anc x y) : d.depth x ≤ d.depth y := by
  induction h with
  | refl => exact Nat.le_refl _
  | tail _ hp ih => rw [hs.depth_parent _ _ hp]; omega

theorem Anc.depth_lt {x y : Nat} (h : d.Anc x y) (hne : x ≠ y) : d.depth x < d.depth y := by
  rcases h.cases_tail with rfl | ⟨p, hp, hpy⟩
  · exact absurd rfl hne
  · have := hp.depth_le hs; rw [hs.depth_parent _ _ hpy]; omega

theorem Anc.antisymm {x y : Nat} (h : d.Anc x y) (h' : d.Anc y x) : x = y := by
  by_contra hne
  have := h.depth_lt hs hne
  have := h'.depth_lt hs (Ne.symm hne)
  omega

theorem Anc.ne_of_depth_lt {x y : Nat} (h : d.Anc x y) (hlt : d.depth y < d.depth x) : False := by
  have := h.depth_le hs; omega

/-- Ancestors of a vertex are linearly ordered. -/
theorem Anc.comparable {u v x : Nat} (hu : d.Anc u x) (hv : d.Anc v x) : d.Anc u v ∨ d.Anc v u := by
  induction hv with
  | refl => exact .inl hu
  | tail hv hp ih =>
    rename_i y z
    rcases hu.cases_tail with rfl | ⟨p, hp', hpz⟩
    · exact .inr (hv.tail hp)
    · rw [hs.parent_unique _ _ _ hpz hp] at hp'
      exact ih hp'

theorem Anc.eq_of_depth_eq {u v x : Nat} (hu : d.Anc u x) (hv : d.Anc v x)
    (h : d.depth u = d.depth v) : u = v := by
  rcases hu.comparable hs hv with h' | h'
  · by_contra hne; have := h'.depth_lt hs hne; omega
  · by_contra hne; have := h'.depth_lt hs (Ne.symm hne); omega

/-- A proper ancestor of `c` is an ancestor of `c`'s parent. -/
theorem Anc.to_parent {v p c : Nat} (hv : d.Anc v c) (hne : v ≠ c) (hp : d.IsParent p c) :
    d.Anc v p := by
  rcases hv.cases_tail with rfl | ⟨q, hq, hqc⟩
  · exact absurd rfl hne
  · rwa [hs.parent_unique _ _ _ hqc hp] at hq

theorem root_anc_eq {v : Nat} (h : d.Anc v d.root) : v = d.root := by
  by_contra hne
  have := h.depth_lt hs hne
  rw [hs.depth_root] at this
  omega

theorem not_isParent_root {p : Nat} (h : d.IsParent p d.root) : False := by
  have := hs.depth_parent _ _ h
  rw [hs.depth_root] at this
  omega

/-- Adjacent vertices are ancestor/descendant: no cross edges. -/
theorem adj_comparable {x y : Nat} (h : g.Adj x y) : d.Anc x y ∨ d.Anc y x := by
  obtain ⟨e, he⟩ := h
  obtain ⟨v, o, ho, rfl⟩ := hs.edge_out e he.lt
  have hj := hs.joins v o ho
  have hvo : d.Anc v o.dest ∨ d.Anc o.dest v := by
    cases o with
    | back => exact .inr (hs.back_anc v _ ho rfl)
    | tree => exact .inl (IsParent.anc ⟨_, ho, rfl, rfl⟩)
  rcases he.eq_or hj with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩
  · exact hvo
  · exact hvo.symm

/-- The tree path from `x` up to its ancestor `y`, inside `ok`. -/
theorem Anc.reach_up {ok : Nat → Prop} {y x : Nat} (h : d.Anc y x)
    (hok : ∀ w, d.Anc y w → d.Anc w x → ok w) : g.Reach ok x y := by
  induction h with
  | refl => exact .refl (hok _ (.refl _) (.refl _))
  | tail hyz hp ih =>
    rename_i z w
    refine Graph.Reach.head (hok _ (hyz.tail hp) (.refl _)) (hp.adj hs).symm (ih ?_)
    intro v h1 h2
    exact hok v h1 (h2.tail hp)

theorem Anc.reach_down {ok : Nat → Prop} {y x : Nat} (h : d.Anc y x)
    (hok : ∀ w, d.Anc y w → d.Anc w x → ok w) : g.Reach ok y x :=
  (h.reach_up hs hok).symm

/-- Two vertices of `T_c` are joined inside `T_c`. -/
theorem reach_in_subtree {ok : Nat → Prop} {c x y : Nat} (hx : d.Anc c x) (hy : d.Anc c y)
    (hok : ∀ w, d.Anc c w → ok w) : g.Reach ok x y :=
  (hx.reach_up hs fun w h _ => hok w h).trans (hy.reach_down hs fun w h _ => hok w h)

/-- A walk from `T_c` to outside `T_c` that avoids `c`'s parent `a` leaves `T_c` along an edge to
a proper ancestor of `a`. -/
theorem exit_subtree {ok : Nat → Prop} {a c x y : Nat} (hp : d.IsParent a c) (hx : d.Anc c x)
    (hy : ¬d.Anc c y) (h : g.Reach ok x y) (ha : ¬ok a) :
    ∃ u v, g.Reach (fun w => ok w ∧ d.Anc c w) x u ∧ g.Adj u v ∧ ok v ∧ d.Anc v a ∧ v ≠ a := by
  rcases h.exit hx with h | ⟨u, v, h1, h2, h3, h4⟩
  · exact absurd h.ok_right.2 hy
  · have hu : d.Anc c u := h1.ok_right.2
    have hva : v ≠ a := fun h => ha (h ▸ h3)
    refine ⟨u, v, h1, h2, h3, ?_, hva⟩
    rcases adj_comparable hs h2 with huv | hvu
    · exact absurd (hu.trans huv) h4
    · rcases hvu.comparable hs hu with hvc | hcv
      · exact hvc.to_parent hs (fun h => h4 (h ▸ .refl _)) hp
      · exact absurd hcv h4

/-- In a block, a child subtree `T_c` of `a` has an edge to a proper ancestor of `a`, provided
some edge has an endpoint outside `T_c ∪ {a}`. -/
theorem exists_edge_above (h2 : g.TwoConnected) {a c : Nat} (hp : d.IsParent a c) {e' y : Nat}
    (hy : g.IsEnd e' y) (hyc : ¬d.Anc c y) (hya : y ≠ a) :
    ∃ u v, d.Anc c u ∧ g.Adj u v ∧ d.Anc v a ∧ v ≠ a := by
  obtain ⟨e, he⟩ := hp.joins hs
  rcases h2 a e e' he.lt hy.lt with rfl | ⟨x, y', hx, hy', hr⟩
  · obtain ⟨z, hz⟩ := hy
    rcases he.eq_or hz with ⟨h, _⟩ | ⟨_, h⟩
    · exact absurd h.symm hya
    · exact absurd (h ▸ .refl _) hyc
  · have hxa : x ≠ a := hr.ok_left
    obtain ⟨z, hz⟩ := hx
    obtain rfl : c = x := by
      rcases he.eq_or hz with ⟨h, _⟩ | ⟨_, h⟩
      · exact absurd h.symm hxa
      · exact h
    by_cases hcy : d.Anc c y'
    · obtain ⟨z, hz⟩ := hy
      obtain ⟨z', hz'⟩ := hy'
      obtain rfl : z = y' := by
        rcases hz.eq_or hz' with ⟨h, _⟩ | ⟨_, h⟩
        · exact absurd (h ▸ hcy) hyc
        · exact h
      refine ⟨_, y, hcy, hz.adj.symm, ?_, hya⟩
      rcases adj_comparable hs hz.adj with h | h
      · rcases h.comparable hs hcy with h' | h'
        · exact h'.to_parent hs (fun h => hyc (h ▸ .refl _)) hp
        · exact absurd h' hyc
      · exact absurd (hcy.trans h) hyc
    · obtain ⟨u, v, h1, h3, _, h4, h5⟩ := exit_subtree hs hp (.refl _) hcy hr (fun h => h rfl)
      exact ⟨u, v, h1.ok_right.2, h3, h4, h5⟩

/-- Along the tree path from a vertex outside `T_a ∪ T_b` to the root nothing is removed. -/
theorem reach_root_of_not_anc {a b x : Nat} (hx : d.Anc d.root x) (ha : ¬d.Anc a x)
    (hb : ¬d.Anc b x) : g.Reach (fun w => w ≠ a ∧ w ≠ b) x d.root :=
  hx.reach_up hs fun _ _ hw => ⟨fun h => ha (h ▸ hw), fun h => hb (h ▸ hw)⟩

/-- Fact A, the core: with `a`, `b` incomparable, every vertex of `T_a − a` reaches the root
avoiding `{a, b}`. -/
theorem reach_root_of_anc (h2 : g.TwoConnected) {a b x : Nat} (hab : ¬d.Anc a b)
    (hba : ¬d.Anc b a) (hra : d.Anc d.root a) (hrb : d.Anc d.root b) (hx : d.Anc a x)
    (hxa : x ≠ a) :
    g.Reach (fun w => w ≠ a ∧ w ≠ b) x d.root := by
  obtain ⟨c, hac, hcx⟩ := hx.child (Ne.symm hxa)
  have hroot : a ≠ d.root := fun h => hab (h ▸ hrb)
  obtain ⟨p, -, hpa⟩ := hra.cases_tail.resolve_left hroot
  obtain ⟨e', he'⟩ := hpa.joins hs
  have hdc := hs.depth_parent _ _ hac
  have hdp := hs.depth_parent _ _ hpa
  obtain ⟨u, v, hcu, huv, hva, hvne⟩ := exists_edge_above hs h2 hac he'.isEnd
    (fun h => h.ne_of_depth_lt hs (by omega)) (fun h => by subst h; omega)
  have hokT : ∀ w, d.Anc c w → w ≠ a ∧ w ≠ b := fun w hw =>
    ⟨fun h => hw.ne_of_depth_lt hs (by subst h; omega), fun h => hab (h ▸ hac.anc.trans hw)⟩
  have hokv : v ≠ a ∧ v ≠ b := ⟨hvne, fun h => hba (h ▸ hva)⟩
  refine ((reach_in_subtree hs hcx hcu hokT).tail huv hokv).trans ?_
  rcases hra.comparable hs hva with hrv | hvr
  · exact reach_root_of_not_anc hs hrv (fun h => hvne (h.antisymm hs hva).symm)
      (fun h => hba (h.trans hva))
  · exact root_anc_eq hs hvr ▸ .refl hokv

/-- Fact A: in a block, the two vertices of a separation pair are ancestor/descendant. -/
theorem sepPair_comparable (hr : d.Rooted g) (h2 : g.TwoConnected) {a b : Nat}
    (ha : d.Anc d.root a) (hb : d.Anc d.root b) (hab : g.SeparationPair a b) :
    d.Anc a b ∨ d.Anc b a := by
  by_contra hn
  have ⟨hab', hba'⟩ : ¬d.Anc a b ∧ ¬d.Anc b a := ⟨fun h => hn (.inl h), fun h => hn (.inr h)⟩
  obtain ⟨e, e', he, he', hsep⟩ := hab.exists_not_sepClass
  have hne : 2 ≤ g.ne := by
    have : e ≠ e' := fun h => hsep (h ▸ .refl _)
    omega
  -- every endpoint not in `{a, b}` reaches the root
  have key : ∀ x, d.Anc d.root x → x ≠ a → x ≠ b → g.Reach (fun w => w ≠ a ∧ w ≠ b) x d.root := by
    intro x hx hxa hxb
    by_cases hax : d.Anc a x
    · exact reach_root_of_anc hs h2 hab' hba' ha hb hax hxa
    by_cases hbx : d.Anc b x
    · exact (reach_root_of_anc hs h2 hba' hab' hb ha hbx hxb).mono fun _ h => h.symm
    exact reach_root_of_not_anc hs hx hax hbx
  -- every edge has an endpoint not in `{a, b}`
  have ends : ∀ e, e < g.ne → ∃ x, g.IsEnd e x ∧ x ≠ a ∧ x ≠ b := by
    intro e he
    obtain ⟨v, o, ho, rfl⟩ := hs.edge_out e he
    have hj := hs.joins v o ho
    by_cases hv : v ≠ a ∧ v ≠ b
    · exact ⟨v, hj.isEnd, hv⟩
    by_cases hw : o.dest ≠ a ∧ o.dest ≠ b
    · exact ⟨o.dest, hj.symm.isEnd, hw⟩
    exfalso
    have hcomp := adj_comparable hs hj.adj
    have hvw : v ≠ o.dest := fun h => h2.not_loop hne (h ▸ hj)
    have : (v = a ∧ o.dest = b) ∨ (v = b ∧ o.dest = a) := by omega
    rcases this with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩
    · exact hcomp.elim hab' hba'
    · exact hcomp.elim hba' hab'
  obtain ⟨x, hx, hxa, hxb⟩ := ends e he
  obtain ⟨y, hy, hya, hyb⟩ := ends e' he'
  exact hsep (.of_reach hx hy ((key x (hr _ _ hx) hxa hxb).trans (key y (hr _ _ hy) hya hyb).symm))

/-! ### Fact B: the classes at `a = anc b l` -/

theorem edge_out_cases {e y z : Nat} (hj : g.Joins e y z) :
    ∃ v, ∃ o ∈ d.outs v, o.e = e ∧ ((v = y ∧ o.dest = z) ∨ (v = z ∧ o.dest = y)) := by
  obtain ⟨v, o, ho, rfl⟩ := hs.edge_out e hj.lt
  exact ⟨v, o, ho, rfl, (hs.joins v o ho).eq_or hj⟩

omit hs in
theorem isParent_of_tree {v : Nat} {o : DfsOut} (ho : o ∈ d.outs v) (ht : o.isTree = true) :
    d.IsParent v o.dest := ⟨o, ho, ht, rfl⟩

/-- An edge from `u` to a proper ancestor `v` is a back edge out of `u`, unless `v` is `u`'s
parent. -/
theorem edge_up {e u v : Nat} (hj : g.Joins e u v) (hv : d.Anc v u) (hne : v ≠ u) :
    d.IsParent v u ∨ ∃ o ∈ d.outs u, o.e = e ∧ o.isTree = false ∧ o.dest = v := by
  obtain ⟨w, o, ho, rfl, ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩⟩ := edge_out_cases hs hj
  · cases ht : o.isTree
    · exact .inr ⟨o, ho, rfl, ht, rfl⟩
    · exact absurd (hv.antisymm hs (isParent_of_tree ho ht).anc) hne
  · cases ht : o.isTree
    · exact absurd (hv.antisymm hs (hs.back_anc _ o ho ht)) hne
    · exact .inl (isParent_of_tree ho ht)

/-- A tree edge into `T_c` from outside `T_c` enters at `c`. -/
theorem IsParent.eq_of_anc {c z y : Nat} (hcy : d.Anc c y) (hzy : d.IsParent z y)
    (hz : ¬d.Anc c z) : y = c := by
  rcases hcy.cases_tail with h | ⟨q, hq, hqy⟩
  · exact h
  · exact absurd (hs.parent_unique _ _ _ hqy hzy ▸ hq) hz

theorem subtree_ok {a b c w : Nat} (hab : d.Anc a b) (hp : d.IsParent b c) (hw : d.Anc c w) :
    w ≠ a ∧ w ≠ b := by
  have h1 := hs.depth_parent _ _ hp
  have h2 := hw.depth_le hs
  have h3 := hab.depth_le hs
  exact ⟨fun h => by subst h; omega, fun h => by subst h; omega⟩

/-- If `T_c` (child of `b`) returns only to `a` and below `b`, walks avoiding `{a, b}` cannot
leave `T_c`. -/
theorem reach_stays_in_subtree {a b c x y : Nat} (hab : d.Anc a b) (hp : d.IsParent b c)
    (hlo : ∀ l, l < d.depth a → ¬d.Returns c l)
    (hmid : ∀ l, d.depth a < l → l < d.depth b → ¬d.Returns c l)
    (hx : d.Anc c x) (h : g.Reach (fun w => w ≠ a ∧ w ≠ b) x y) : d.Anc c y := by
  induction h with
  | refl => exact hx
  | tail _ hadj hok ih =>
    rename_i y z _
    rcases adj_comparable hs hadj with hyz | hzy
    · exact ih.trans hyz
    by_contra hcz
    have hne : z ≠ y := fun h => hcz (h ▸ ih)
    obtain ⟨e, he⟩ := hadj
    rcases edge_up hs he hzy hne with hpar | ⟨o, ho, -, hback, rfl⟩
    · obtain rfl := IsParent.eq_of_anc hs ih hpar hcz
      exact hok.2 (hs.parent_unique _ _ _ hpar hp)
    · have hret : d.Returns c (d.depth o.dest) := ⟨y, o, ho, ih, hback, rfl⟩
      have hzc : d.Anc o.dest c := (hzy.comparable hs ih).resolve_right hcz
      have hzb : d.Anc o.dest b := hzc.to_parent hs (fun h => hcz (h ▸ .refl _)) hp
      have hd := hzb.depth_le hs
      rcases Nat.lt_trichotomy (d.depth o.dest) (d.depth a) with hl | hl | hl
      · exact hlo _ hl hret
      · exact hok.1 (hzb.eq_of_depth_eq hs hab hl)
      · rcases Nat.lt_or_ge (d.depth o.dest) (d.depth b) with hl' | hl'
        · exact hmid _ hl hl' hret
        · exact hok.2 (hzb.eq_of_depth_eq hs (.refl _) (by omega))

/-- Fact B, type-1 children: if `b → c` is classified `ret (depth a) type1Child` then the edges
with an endpoint in `T_c` (`b–c`, `E(T_c)`, back edges out of `T_c`) form one separation class of
`{a, b}`. -/
theorem type1_class {a b : Nat} {o : DfsOut} (hab : d.Anc a b) (ho : o ∈ d.outs b)
    (hcls : o.cls = .ret (d.depth a) .type1Child) {e e' : Nat} (he : d.EndIn o.dest e g) :
    g.SepClass a b e e' ↔ d.EndIn o.dest e' g := by
  obtain ⟨ht, -, -, hlo, hmid⟩ := (hs.cls_type1 b o ho _).mp hcls
  have hp : d.IsParent b o.dest := isParent_of_tree ho ht
  obtain ⟨x, hx, hcx⟩ := he
  constructor
  · rintro (rfl | ⟨y, z, hy, hz, hr⟩)
    · exact ⟨x, hx, hcx⟩
    · obtain ⟨x', hx'⟩ := hx
      obtain ⟨y', hy'⟩ := hy
      have hcy : d.Anc o.dest y := by
        rcases hx'.eq_or hy' with ⟨h, -⟩ | ⟨-, h⟩
        · exact h ▸ hcx
        · exact reach_stays_in_subtree hs hab hp hlo hmid hcx
            (Graph.Reach.tail (.refl (subtree_ok hs hab hp hcx)) (h ▸ hx'.adj) hr.ok_left)
      exact ⟨z, hz, reach_stays_in_subtree hs hab hp hlo hmid hcy hr⟩
  · rintro ⟨y, hy, hcy⟩
    exact .of_reach hx hy (reach_in_subtree hs hcx hcy fun w hw => subtree_ok hs hab hp hw)

/-- A child `c` of `b` returning above `a` attaches to the *above* part. -/
theorem attaches_above {a a' b c l : Nat} (hab : d.Anc a b) (hpa : d.IsParent a a')
    (h : d.Returns c l) (hl : l < d.depth a) :
    ∃ e, d.EndIn c e g ∧ d.Above a' a b e g := by
  obtain ⟨u, o, ho, hcu, -, rfl⟩ := h
  have hj := hs.joins u o ho
  have hd := hab.depth_le hs
  refine ⟨o.e, ⟨u, hj.isEnd, hcu⟩, o.dest, hj.symm.isEnd, ?_, ?_, ?_⟩
  · exact fun h => by rw [h] at hl; omega
  · exact fun h => by rw [h] at hl; omega
  · exact fun h => h.ne_of_depth_lt hs (by rw [hs.depth_parent _ _ hpa]; omega)

/-- A child `c` of `b` returning strictly between `a` and `b` attaches to the *between* part. -/
theorem attaches_between {a a' b c l : Nat} (hab : d.Anc a b) (hpa : d.IsParent a a')
    (ha'b : d.Anc a' b) (hp : d.IsParent b c) (h : d.Returns c l) (hl : d.depth a < l)
    (hl' : l < d.depth b) : ∃ e, d.EndIn c e g ∧ d.Between a' b e g := by
  obtain ⟨u, o, ho, hcu, hback, rfl⟩ := h
  have hj := hs.joins u o ho
  have hwu : d.Anc o.dest u := hs.back_anc u o ho hback
  have hbu : d.Anc b u := hp.anc.trans hcu
  have hwb : d.Anc o.dest b :=
    (hwu.comparable hs hbu).resolve_right fun h => h.ne_of_depth_lt hs hl'
  have haw : d.Anc a o.dest :=
    (hab.comparable hs hwb).resolve_right fun h => h.ne_of_depth_lt hs hl
  obtain ⟨c', hc', hc'w⟩ := haw.child fun h => by rw [h] at hl; omega
  obtain rfl : c' = a' := (hc'w.trans hwb).eq_of_depth_eq hs ha'b
    (by rw [hs.depth_parent _ _ hc', hs.depth_parent _ _ hpa])
  exact ⟨o.e, ⟨u, hj.isEnd, hcu⟩, o.dest, hj.symm.isEnd, hc'w,
    fun h => h.ne_of_depth_lt hs hl'⟩

/-- The *between* part `T_{a'} − T_b` is connected avoiding `{a, b}`. -/
theorem between_reach {a a' b x y : Nat} (hpa : d.IsParent a a') (hx : d.Anc a' x)
    (hbx : ¬d.Anc b x) (hy : d.Anc a' y) (hby : ¬d.Anc b y) :
    g.Reach (fun w => w ≠ a ∧ w ≠ b) x y := by
  have hok : ∀ z, ¬d.Anc b z → ∀ w, d.Anc a' w → d.Anc w z → w ≠ a ∧ w ≠ b :=
    fun z hz w hw hwz =>
      ⟨fun h => hw.ne_of_depth_lt hs (by subst h; rw [hs.depth_parent _ _ hpa]; omega),
        fun h => hz (h ▸ hwz)⟩
  exact (hx.reach_up hs (hok x hbx)).trans (hy.reach_down hs (hok y hby))

/-- In a block, every child subtree of a non-root vertex `b` returns strictly above `b`
(`lowpt1 c < depth b`). -/
theorem child_returns_above (h2 : g.TwoConnected) {p b c : Nat} (hpb : d.IsParent p b)
    (hp : d.IsParent b c) : ∃ l, l < d.depth b ∧ d.Returns c l := by
  obtain ⟨e', he'⟩ := hpb.joins hs
  have hdb := hs.depth_parent _ _ hpb
  have hdc := hs.depth_parent _ _ hp
  obtain ⟨u, v, hcu, huv, hvb, hvb'⟩ := exists_edge_above hs h2 hp he'.isEnd
    (fun h => h.ne_of_depth_lt hs (by omega)) (fun h => by subst h; omega)
  have hbu : d.Anc b u := hp.anc.trans hcu
  have hvu : d.Anc v u := hvb.trans hbu
  have hne : v ≠ u := fun h => hvb' (hvb.antisymm hs (h ▸ hbu))
  obtain ⟨e, he⟩ := huv
  rcases edge_up hs he hvu hne with hpar | ⟨o, ho, -, hback, rfl⟩
  · have hcv : ¬d.Anc c v := fun h =>
      h.ne_of_depth_lt hs (by have := hvb.depth_lt hs hvb'; omega)
    obtain rfl := IsParent.eq_of_anc hs hcu hpar hcv
    exact absurd (hs.parent_unique _ _ _ hpar hp) hvb'
  · exact ⟨_, hvb.depth_lt hs hvb', u, o, ho, hcu, hback, rfl⟩

/-- In a block, every tree out-edge of a non-root vertex is a `ret` edge. -/
theorem cls_ret_of_tree (h2 : g.TwoConnected) {p b : Nat} {o : DfsOut} (hpb : d.IsParent p b)
    (ho : o ∈ d.outs b) (ht : o.isTree = true) : ∃ l k, o.cls = .ret l k ∧ k ≠ .backEdge := by
  obtain ⟨l, hl, hret⟩ := child_returns_above hs h2 hpb (isParent_of_tree ho ht)
  cases hcls : o.cls with
  | bridge => exact absurd hret (((hs.cls_bridge b o ho).mp hcls).2 l (Nat.le_of_lt hl))
  | component => exact absurd hret (((hs.cls_component b o ho).mp hcls).2.2 l hl)
  | selfLoop => exact absurd ht (by rw [((hs.cls_selfLoop b o ho).mp hcls).1]; decide)
  | ret l' k =>
    refine ⟨l', k, rfl, fun h => ?_⟩
    subst h
    exact absurd ht (by rw [((hs.cls_backEdge b o ho l').mp hcls).1]; decide)

/-- The lowpoint data of a `ret` tree edge. -/
theorem ret_lowpt {v l : Nat} {k : RetKind} {o : DfsOut} (ho : o ∈ d.outs v)
    (ht : o.isTree = true) (hc : o.cls = .ret l k) :
    l < d.depth v ∧ d.Returns o.dest l ∧ ∀ l', l' < l → ¬d.Returns o.dest l' := by
  cases k with
  | type1Child =>
    obtain ⟨-, h1, h2, h3, -⟩ := (hs.cls_type1 v o ho l).mp hc
    exact ⟨h1, h2, h3⟩
  | type2Child =>
    obtain ⟨-, h1, h2, h3, -⟩ := (hs.cls_type2 v o ho l).mp hc
    exact ⟨h1, h2, h3⟩
  | backEdge => exact absurd ht (by rw [((hs.cls_backEdge v o ho l).mp hc).1]; decide)

theorem rank_le_of_sublist {v : Nat} {o o' : DfsOut} (h : [o, o'].Sublist (d.outs v)) :
    o.cls.rank ≤ o'.cls.rank :=
  (List.pairwise_cons.mp ((hs.sorted v).sublist h)).1 o' (List.mem_singleton.mpr rfl)

end Spec

end DfsData

namespace RetKind

theorem rank_le (k : RetKind) : k.rank ≤ 2 := by cases k <;> decide

end RetKind

namespace OutClass

theorem AttachesAbove.ret {l : Nat} {c : OutClass} (h : c.AttachesAbove l) :
    ∃ l' k, c = .ret l' k ∧ l' < l := by
  cases c with
  | ret l' k => exact ⟨l', k, rfl, h⟩
  | bridge => exact h.elim
  | component => exact h.elim
  | selfLoop => exact h.elim

theorem AttachesBetween.ret {l : Nat} {c : OutClass} (h : c.AttachesBetween l) :
    ∃ l' k, c = .ret l' k ∧ l ≤ l' ∧ (l = l' → k = .type2Child) := by
  cases c with
  | ret l' k =>
    cases k with
    | type2Child => exact ⟨l', _, rfl, h, fun _ => rfl⟩
    | type1Child =>
      exact ⟨l', _, rfl, Nat.le_of_lt h, fun h' => by subst h'; exact absurd h (Nat.lt_irrefl _)⟩
    | backEdge =>
      exact ⟨l', _, rfl, Nat.le_of_lt h, fun h' => by subst h'; exact absurd h (Nat.lt_irrefl _)⟩
  | bridge => exact h.elim
  | component => exact h.elim
  | selfLoop => exact h.elim

theorem attachesBetween_of_lt {l l' : Nat} {k : RetKind} (h : l < l') :
    (ret l' k).AttachesBetween l := by
  cases k
  · exact h
  · exact h
  · exact Nat.le_of_lt h

end OutClass

namespace DfsData

variable {g : Graph} {d : DfsData}

section Spec

variable (hs : d.Spec g)
include hs

/-- Fact B, the sorted out-list: before a child attaching *above*, every tree edge attaches
*above* (the *above* children form a prefix). -/
theorem above_prefix (h2 : g.TwoConnected) {p b l : Nat} {o o' : DfsOut} (hpb : d.IsParent p b)
    (h : [o, o'].Sublist (d.outs b)) (ht : o.isTree = true) (h' : o'.cls.AttachesAbove l) :
    o.cls.AttachesAbove l := by
  have hr := rank_le_of_sublist hs h
  obtain ⟨l₁, k₁, hc, -⟩ := cls_ret_of_tree hs h2 hpb (h.subset (by simp)) ht
  obtain ⟨l₂, k₂, hc', hl⟩ := h'.ret
  rw [hc, hc'] at hr
  have hr' : 3 * (l₁ + 2) + k₁.rank ≤ 3 * (l₂ + 2) + k₂.rank := hr
  have := k₁.rank_le
  have := k₂.rank_le
  rw [hc]
  show l₁ < l
  omega

/-- Fact B, the sorted out-list: after a child attaching *between*, every tree edge attaches
*between* (the *between* children form a suffix). -/
theorem between_suffix (h2 : g.TwoConnected) {p b l : Nat} {o o' : DfsOut}
    (hpb : d.IsParent p b) (h : [o, o'].Sublist (d.outs b)) (ht' : o'.isTree = true)
    (h' : o.cls.AttachesBetween l) : o'.cls.AttachesBetween l := by
  have hr := rank_le_of_sublist hs h
  obtain ⟨l₂, k₂, hc', hk⟩ := cls_ret_of_tree hs h2 hpb (h.subset (by simp)) ht'
  obtain ⟨l₁, k₁, hc, hl, hk₁⟩ := h'.ret
  rw [hc, hc'] at hr
  have hr' : 3 * (l₁ + 2) + k₁.rank ≤ 3 * (l₂ + 2) + k₂.rank := hr
  have := k₁.rank_le
  have := k₂.rank_le
  rw [hc']
  rcases Nat.lt_or_ge l l₂ with hlt | hge
  · exact OutClass.attachesBetween_of_lt hlt
  obtain rfl : l₂ = l := by omega
  obtain rfl : l₁ = l₂ := by omega
  obtain rfl := hk₁ rfl
  have e1 : RetKind.rank .type2Child = 2 := rfl
  cases k₂ with
  | type2Child => exact Nat.le_refl _
  | type1Child => have e2 : RetKind.rank .type1Child = 0 := rfl; omega
  | backEdge => exact absurd rfl hk

/-- Every tree out-edge of `b` attaches above, is a type-1 edge to `l`, or attaches between. -/
theorem tree_out_trichotomy (h2 : g.TwoConnected) {p b : Nat} (l : Nat) {o : DfsOut}
    (hpb : d.IsParent p b) (ho : o ∈ d.outs b) (ht : o.isTree = true) :
    o.cls.AttachesAbove l ∨ o.cls.IsType1To l ∨ o.cls.AttachesBetween l := by
  obtain ⟨l', k, hc, hk⟩ := cls_ret_of_tree hs h2 hpb ho ht
  rw [hc]
  rcases Nat.lt_trichotomy l' l with h | rfl | h
  · exact .inl h
  · cases k with
    | type1Child => exact .inr (.inl rfl)
    | type2Child => exact .inr (.inr (Nat.le_refl _))
    | backEdge => exact absurd rfl hk
  · exact .inr (.inr (OutClass.attachesBetween_of_lt h))

/-- A tree edge classified *above* `l` has `lowpt1 < l`. -/
theorem AttachesAbove.returns {b l : Nat} {o : DfsOut} (ho : o ∈ d.outs b) (ht : o.isTree = true)
    (h : o.cls.AttachesAbove l) : ∃ l', l' < l ∧ d.Returns o.dest l' := by
  obtain ⟨l', k, hc, hl⟩ := h.ret
  exact ⟨l', hl, (ret_lowpt hs ho ht hc).2.1⟩

/-- A tree edge classified *between* `l` and `b` has a return in `(l, depth b)`. -/
theorem AttachesBetween.returns {b l : Nat} {o : DfsOut} (ho : o ∈ d.outs b)
    (ht : o.isTree = true) (h : o.cls.AttachesBetween l) :
    ∃ l', l < l' ∧ l' < d.depth b ∧ d.Returns o.dest l' := by
  obtain ⟨l', k, hc, hl, hk⟩ := h.ret
  obtain ⟨hd, hret, -⟩ := ret_lowpt hs ho ht hc
  rcases Nat.lt_or_eq_of_le hl with hlt | rfl
  · exact ⟨l', hlt, hd, hret⟩
  · obtain ⟨-, -, -, -, l'', h1, h2, h3⟩ :=
      (hs.cls_type2 b o ho _).mp (by rw [hk rfl] at hc; exact hc)
    exact ⟨l'', h1, h2, h3⟩

end Spec

end DfsData

end Spqr
