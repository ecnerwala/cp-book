import Spqr.Blocks
import Spqr.Proofs.SepPair

/-!
# Blocks are the `lowval ≥ d` branches (PROOF.md §2, Lemma 2.1)

Over a `DfsData` satisfying the phase-1 `Spec`:
* `blockRoot_cut`: the parent `p` of a `component` / `bridge` child `c` separates the edges with
  deeper endpoint in `T_c` from all other edges; `deepIn_of_sameBlock` restates this as "the edges
  with deeper endpoint in `T_c` are a union of blocks".
* `ret_child_sameBlock`, `ret_back_sameBlock`: for a `ret` child `c` of `v`, the edges with
  deeper endpoint in `T_c` that are not inside a nested block, in particular the back edges
  returning to `lowpt1 c`, lie in the block of the tree edge `v → c`.
* `sameBlock_iff`: `SameBlock e e' ↔ e = e' ∨ ∃ c, InBlock c e ∧ InBlock c e'` — the blocks are
  the singletons `{bridge}`, `{loop}` and, for every block root `c`, the non-bridge non-loop edges
  whose deeper endpoint has block top `c`.
-/

namespace Spqr

namespace DfsOut

theorem deep_tree {v : Nat} {o : DfsOut} (h : o.isTree = true) : o.deep v = o.dest := by
  cases o <;> simp_all [deep, dest, isTree]

theorem deep_back {v : Nat} {o : DfsOut} (h : o.isTree = false) : o.deep v = v := by
  cases o <;> simp_all [deep, isTree]

theorem shallow_tree {v : Nat} {o : DfsOut} (h : o.isTree = true) : o.shallow v = v := by
  cases o <;> simp_all [shallow, isTree]

theorem shallow_back {v : Nat} {o : DfsOut} (h : o.isTree = false) : o.shallow v = o.dest := by
  cases o <;> simp_all [shallow, dest, isTree]

theorem ends_cases (v : Nat) (o : DfsOut) :
    (o.isTree = true ∧ o.deep v = o.dest ∧ o.shallow v = v) ∨
      (o.isTree = false ∧ o.deep v = v ∧ o.shallow v = o.dest) := by
  cases o <;> simp [deep, shallow, dest, isTree]

end DfsOut

namespace Graph

variable {g : Graph}

theorem SameBlock.symm {e e' : Nat} (h : g.SameBlock e e') : g.SameBlock e' e :=
  fun v => (h v).symm

theorem SameBlock.refl (e : Nat) : g.SameBlock e e := fun _ => .inl rfl

theorem twoConnected_iff_sameBlock :
    g.TwoConnected ↔ ∀ e e', e < g.ne → e' < g.ne → g.SameBlock e e' :=
  ⟨fun h e e' he he' v => h v e e' he he', fun h v e e' he he' => h e e' he he' v⟩

end Graph

namespace DfsData

variable {g : Graph} {d : DfsData}

theorem BlockTop.self {c : Nat} (h : d.BlockRoot c) : d.BlockTop c c :=
  ⟨h, .refl _, fun _ _ h => h⟩

theorem BlockTop.of_anc {c w x : Nat} (h : d.BlockTop c x) (hcw : d.Anc c w) (hwx : d.Anc w x) :
    d.BlockTop c w :=
  ⟨h.1, hcw, fun c' hc' hc'w => h.2.2 c' hc' (hc'w.trans hwx)⟩

theorem IsRoot.anc_eq {c x : Nat} (hr : d.IsRoot c) (h : d.Anc x c) : x = c := by
  rcases h.cases_tail with h | ⟨p, -, hp⟩
  · exact h.symm
  · exact absurd hp (hr p)

theorem IsRoot.exists_parent {c : Nat} (h : ¬d.IsRoot c) : ∃ p, d.IsParent p c := by
  simpa [IsRoot] using h

section Spec

variable (hs : d.Spec g)
include hs

theorem IsParent.ne {p c : Nat} (hp : d.IsParent p c) : p ≠ c := fun h => by
  have := hs.depth_parent _ _ hp; rw [h] at this; omega

theorem IsParent.not_anc {p c : Nat} (hp : d.IsParent p c) : ¬d.Anc c p := fun h =>
  hp.ne hs (hp.anc.antisymm hs h)

theorem BlockTop.unique {c c' x : Nat} (h : d.BlockTop c x) (h' : d.BlockTop c' x) : c = c' :=
  (h'.2.2 c h.1 h.2.1).antisymm hs (h.2.2 c' h'.1 h'.2.1)

omit hs in
theorem BlockTop.parent {c x : Nat} (h : d.BlockTop c x) (hne : x ≠ c) :
    ∃ q, d.IsParent q x ∧ d.BlockTop c q := by
  rcases h.2.1.cases_tail with h' | ⟨q, hq, hqx⟩
  · exact absurd h' hne
  · exact ⟨q, hqx, h.1, hq, fun c' hc' hc'q => h.2.2 c' hc' (hc'q.tail hqx)⟩

theorem BlockTop.not_blockRoot {c x : Nat} (h : d.BlockTop c x) (hne : x ≠ c) :
    ¬d.BlockRoot x := fun hx => hne ((h.2.2 x hx (.refl _)).antisymm hs h.2.1)

/-- Every vertex has a block top. -/
theorem exists_blockTop (x : Nat) : ∃ c, d.BlockTop c x := by
  suffices ∀ n x, d.depth x = n → ∃ c, d.BlockTop c x from this _ x rfl
  intro n
  induction n using Nat.strongRecOn with
  | _ n ih =>
    intro x hn
    by_cases hr : d.BlockRoot x
    · exact ⟨x, BlockTop.self hr⟩
    · obtain ⟨q, hq⟩ := IsRoot.exists_parent fun h => hr (.inl h)
      obtain ⟨c, hc⟩ := ih (d.depth q) (by rw [hs.depth_parent _ _ hq] at hn; omega) q rfl
      exact ⟨c, hc.1, hc.2.1.tail hq, fun c' hc' hc'x =>
        hc.2.2 c' hc' (hc'x.to_parent hs (fun h => hr (h ▸ hc')) hq)⟩

/-- The class of a tree out-edge. -/
theorem tree_cls {v : Nat} {o : DfsOut} (ho : o ∈ d.outs v) (ht : o.isTree = true) :
    o.cls = .bridge ∨ o.cls = .component ∨ ∃ l k, o.cls = .ret l k := by
  cases hc : o.cls with
  | bridge => exact .inl rfl
  | component => exact .inr (.inl rfl)
  | selfLoop => exact absurd ht (by rw [((hs.cls_selfLoop v o ho).mp hc).1]; decide)
  | ret l k => exact .inr (.inr ⟨l, k, rfl⟩)

theorem not_bridge_of_back {v : Nat} {o : DfsOut} (ho : o ∈ d.outs v) (hb : o.isTree = false) :
    o.cls ≠ .bridge := fun h => by
  have := ((hs.cls_bridge v o ho).mp h).1; rw [hb] at this; cases this

theorem dest_ne_of_not_selfLoop {v : Nat} {o : DfsOut} (ho : o ∈ d.outs v) (hb : o.isTree = false)
    (h : o.cls ≠ .selfLoop) : o.dest ≠ v := fun h' => h ((hs.cls_selfLoop v o ho).mpr ⟨hb, h'⟩)

/-- A block root entered by a tree edge is a `component` or `bridge` child, so `T_c` has no return
strictly above the parent. -/
theorem BlockRoot.parent_cls {c p : Nat} (hc : d.BlockRoot c) (hp : d.IsParent p c) :
    ∃ o ∈ d.outs p, o.isTree = true ∧ o.dest = c ∧ (o.cls = .component ∨ o.cls = .bridge) := by
  rcases hc with hr | ⟨p', o, ho, ht, hd, hcls⟩
  · exact absurd hp (hr p)
  · obtain rfl := hs.parent_unique _ _ _ ⟨o, ho, ht, hd⟩ hp
    exact ⟨o, ho, ht, hd, hcls⟩

theorem BlockRoot.no_returns {c p : Nat} (hc : d.BlockRoot c) (hp : d.IsParent p c) :
    ∀ l, l < d.depth p → ¬d.Returns c l := by
  obtain ⟨o, ho, -, rfl, h | h⟩ := hc.parent_cls hs hp
  · exact ((hs.cls_component p o ho).mp h).2.2
  · exact fun l hl => ((hs.cls_bridge p o ho).mp h).2 l (Nat.le_of_lt hl)

omit hs in
theorem BlockRoot.of_tree {p : Nat} {o : DfsOut} (ho : o ∈ d.outs p) (ht : o.isTree = true)
    (h : o.cls = .component ∨ o.cls = .bridge) : d.BlockRoot o.dest :=
  .inr ⟨p, o, ho, ht, rfl, h⟩

/-- A `ret` child is not a block root. -/
theorem not_blockRoot_of_ret {v l : Nat} {k : RetKind} {o : DfsOut} (ho : o ∈ d.outs v)
    (ht : o.isTree = true) (hc : o.cls = .ret l k) : ¬d.BlockRoot o.dest := fun h => by
  obtain ⟨o', ho', ht', hd', hcls'⟩ := h.parent_cls hs (isParent_of_tree ho ht)
  obtain rfl := hs.tree_inj v o' ho' o ho ht' ht hd'
  rcases hcls' with h | h <;> rw [hc] at h <;> cases h

/-- The parent edge of a non-block-root is a `ret` edge: its subtree returns above the parent. -/
theorem returns_of_not_blockRoot {q x : Nat} (hq : d.IsParent q x) (hx : ¬d.BlockRoot x) :
    ∃ l, d.Returns x l ∧ l < d.depth q := by
  obtain ⟨o, ho, ht, rfl⟩ := hq
  rcases tree_cls hs ho ht with h | h | ⟨l, k, h⟩
  · exact absurd (BlockRoot.of_tree ho ht (.inr h)) hx
  · exact absurd (BlockRoot.of_tree ho ht (.inl h)) hx
  · obtain ⟨hl, hr, -⟩ := ret_lowpt hs ho ht h
    exact ⟨l, hr, hl⟩

omit hs in
theorem Returns.mono {c c' l : Nat} (h : d.Returns c' l) (hcc' : d.Anc c c') : d.Returns c l := by
  obtain ⟨u, o, ho, hu, hb, hl⟩ := h
  exact ⟨u, o, ho, hcc'.trans hu, hb, hl⟩

/-- A block root with a nontrivial block is a `component` child (in particular not a root). -/
theorem exists_parent_of_blockTop {c x : Nat} (h : d.BlockTop c x) (hne : x ≠ c) :
    ∃ p, d.IsParent p c := by
  obtain ⟨u, hcu, hux⟩ := h.2.1.child (Ne.symm hne)
  have hu : ¬d.BlockRoot u := fun hu =>
    hcu.ne hs ((h.2.2 u hu hux).antisymm hs hcu.anc).symm
  obtain ⟨l, hl, hlt⟩ := returns_of_not_blockRoot hs hcu hu
  obtain ⟨s, o, ho, hus, hb, rfl⟩ := hl
  have hcs : d.Anc c s := hcu.anc.trans hus
  have hwc : d.Anc o.dest c :=
    ((hs.back_anc s o ho hb).comparable hs hcs).resolve_right fun h => h.ne_of_depth_lt hs hlt
  rcases hwc.cases_tail with h | ⟨p, -, hp⟩
  · rw [h] at hlt; omega
  · exact ⟨p, hp⟩

/-- Endpoints of an out-edge. -/
theorem isEnd_iff {v x : Nat} {o : DfsOut} (ho : o ∈ d.outs v) :
    g.IsEnd o.e x ↔ x = o.deep v ∨ x = o.shallow v := by
  have hj := hs.joins v o ho
  have : g.IsEnd o.e x ↔ x = v ∨ x = o.dest := by
    constructor
    · rintro ⟨y, hy⟩
      rcases hy.eq_or hj with ⟨rfl, -⟩ | ⟨rfl, -⟩
      · exact .inl rfl
      · exact .inr rfl
    · rintro (rfl | rfl)
      · exact hj.isEnd
      · exact hj.symm.isEnd
  rw [this]
  rcases DfsOut.ends_cases v o with ⟨-, h1, h2⟩ | ⟨-, h1, h2⟩ <;> rw [h1, h2] <;> exact or_comm

theorem isEnd_deep {v : Nat} {o : DfsOut} (ho : o ∈ d.outs v) : g.IsEnd o.e (o.deep v) :=
  (isEnd_iff hs ho).mpr (.inl rfl)

theorem isEnd_shallow {v : Nat} {o : DfsOut} (ho : o ∈ d.outs v) : g.IsEnd o.e (o.shallow v) :=
  (isEnd_iff hs ho).mpr (.inr rfl)

theorem shallow_anc_deep {v : Nat} {o : DfsOut} (ho : o ∈ d.outs v) :
    d.Anc (o.shallow v) (o.deep v) := by
  rcases DfsOut.ends_cases v o with ⟨ht, h1, h2⟩ | ⟨hb, h1, h2⟩ <;> rw [h1, h2]
  · exact (isParent_of_tree ho ht).anc
  · exact hs.back_anc v o ho hb

theorem deep_ne_shallow {v : Nat} {o : DfsOut} (ho : o ∈ d.outs v) (hn : ¬o.isSingleton) :
    o.deep v ≠ o.shallow v := by
  rcases DfsOut.ends_cases v o with ⟨ht, h1, h2⟩ | ⟨hb, h1, h2⟩ <;> rw [h1, h2]
  · exact Ne.symm ((isParent_of_tree ho ht).ne hs)
  · exact Ne.symm (dest_ne_of_not_selfLoop hs ho hb fun h => hn (.inr h))

theorem adj_deep_shallow {v : Nat} {o : DfsOut} (ho : o ∈ d.outs v) :
    g.Adj (o.deep v) (o.shallow v) := by
  have hj := (hs.joins v o ho).adj
  rcases DfsOut.ends_cases v o with ⟨-, h1, h2⟩ | ⟨-, h1, h2⟩ <;> rw [h1, h2]
  · exact hj.symm
  · exact hj

/-- An edge index determines its out-edge. -/
theorem out_eq_of_e {v v' : Nat} {o o' : DfsOut} (ho : o ∈ d.outs v) (ho' : o' ∈ d.outs v')
    (h : o.e = o'.e) : v = v' ∧ o = o' := hs.out_inj v v' o ho o' ho' h

/-! ### Leaving a subtree -/

/-- A walk leaving `T_c` (child of `p`) does so along the tree edge `c–p` or along a back edge
from `T_c` to a proper ancestor `z` of `c`. -/
theorem exit_edge {ok : Nat → Prop} {p c x y : Nat} (hp : d.IsParent p c) (hx : d.Anc c x)
    (hy : ¬d.Anc c y) (h : g.Reach ok x y) :
    ∃ u z, d.Anc c u ∧ ok u ∧ ok z ∧ d.Anc z p ∧ ((u = c ∧ z = p) ∨ d.Returns c (d.depth z)) := by
  rcases h.exit hx with h | ⟨u, z, h1, h2, h3, h4⟩
  · exact absurd h.ok_right.2 hy
  · have hu : d.Anc c u := h1.ok_right.2
    have hzu : d.Anc z u := by
      rcases adj_comparable hs h2 with h | h
      · exact absurd (hu.trans h) h4
      · exact h
    have hzu' : z ≠ u := fun h => h4 (h ▸ hu)
    have hzc : d.Anc z c := (hzu.comparable hs hu).resolve_right h4
    have hzc' : z ≠ c := fun h => h4 (h ▸ .refl _)
    refine ⟨u, z, hu, h1.ok_right.1, h3, hzc.to_parent hs hzc' hp, ?_⟩
    obtain ⟨e, he⟩ := h2
    rcases edge_up hs he hzu hzu' with hpar | ⟨o, ho, -, hback, rfl⟩
    · have huc := IsParent.eq_of_anc hs hu hpar h4
      rw [huc] at hpar
      exact .inl ⟨huc, hs.parent_unique _ _ _ hpar hp⟩
    · exact .inr ⟨u, o, ho, hu, hback, rfl⟩

/-- Removing the parent `p` of a block root `c` disconnects `T_c` from the outside. -/
theorem BlockRoot.sep {ok : Nat → Prop} {c p x y : Nat} (hc : d.BlockRoot c) (hp : d.IsParent p c)
    (hx : d.Anc c x) (hy : ¬d.Anc c y) (h : g.Reach ok x y) : ok p := by
  obtain ⟨u, z, hu, hou, hoz, hzp, hcase⟩ := exit_edge hs hp hx hy h
  rcases hcase with ⟨-, rfl⟩ | hret
  · exact hoz
  · by_cases hz : z = p
    · exact hz ▸ hoz
    · exact absurd hret (hc.no_returns hs hp _ (hzp.depth_lt hs hz))

/-- Removing the child end `c` of a bridge disconnects `T_c − c` from the outside. -/
theorem bridge_sep {ok : Nat → Prop} {p x y : Nat} {o : DfsOut} (ho : o ∈ d.outs p)
    (hb : o.cls = .bridge) (hx : d.Anc o.dest x) (hy : ¬d.Anc o.dest y) (h : g.Reach ok x y) :
    ok o.dest := by
  obtain ⟨ht, hno⟩ := (hs.cls_bridge p o ho).mp hb
  obtain ⟨u, z, hu, hou, hoz, hzp, hcase⟩ := exit_edge hs (isParent_of_tree ho ht) hx hy h
  rcases hcase with ⟨rfl, -⟩ | hret
  · exact hou
  · exact absurd hret (hno _ (hzp.depth_le hs))

/-- Nothing leaves the tree of a root. -/
theorem IsRoot.reach {ok : Nat → Prop} {c x y : Nat} (hr : d.IsRoot c) (hx : d.Anc c x)
    (h : g.Reach ok x y) : d.Anc c y := by
  induction h with
  | refl => exact hx
  | tail _ hadj _ ih =>
    rcases adj_comparable hs hadj with h | h
    · exact ih.trans h
    · rcases h.comparable hs ih with h' | h'
      · rw [hr.anc_eq h']; exact .refl _
      · exact h'

/-! ### Where the endpoints of an edge lie -/

/-- If the deeper endpoint of an edge is in `T_c` for a block root `c` with parent `p`, its
shallower endpoint is in `T_c` or is `p`. -/
theorem shallow_mem {c p v : Nat} {o : DfsOut} (hc : d.BlockRoot c) (hp : d.IsParent p c)
    (ho : o ∈ d.outs v) (hd : d.Anc c (o.deep v)) : d.Anc c (o.shallow v) ∨ o.shallow v = p := by
  rcases DfsOut.ends_cases v o with ⟨ht, h1, h2⟩ | ⟨hb, h1, h2⟩ <;> rw [h1] at hd <;> rw [h2]
  · by_cases hv : d.Anc c v
    · exact .inl hv
    · have hpar := isParent_of_tree ho ht
      have := IsParent.eq_of_anc hs hd hpar hv
      rw [this] at hpar
      exact .inr (hs.parent_unique _ _ _ hpar hp)
  · have hwv := hs.back_anc v o ho hb
    rcases hwv.comparable hs hd with hwc | hcw
    · by_cases hwc' : o.dest = c
      · exact .inl (hwc' ▸ .refl _)
      · have hwp := hwc.to_parent hs hwc' hp
        have hret : d.Returns c (d.depth o.dest) := ⟨v, o, ho, hd, hb, rfl⟩
        by_contra hne
        exact hc.no_returns hs hp _ (hwp.depth_lt hs fun h => hne (.inr h)) hret
    · exact .inl hcw

theorem shallow_mem_root {c v : Nat} {o : DfsOut} (hr : d.IsRoot c) (ho : o ∈ d.outs v)
    (hd : d.Anc c (o.deep v)) : d.Anc c (o.shallow v) := by
  rcases (shallow_anc_deep hs ho).comparable hs hd with h | h
  · rw [hr.anc_eq h]; exact .refl _
  · exact h

theorem shallow_blockTop {c p v : Nat} {o : DfsOut} (hp : d.IsParent p c) (ho : o ∈ d.outs v)
    (hd : d.BlockTop c (o.deep v)) : d.BlockTop c (o.shallow v) ∨ o.shallow v = p :=
  (shallow_mem hs hd.1 hp ho hd.2.1).imp_left fun h => hd.of_anc h (shallow_anc_deep hs ho)

theorem ends_in {c p v x : Nat} {o : DfsOut} (hc : d.BlockRoot c) (hp : d.IsParent p c)
    (ho : o ∈ d.outs v) (hd : d.Anc c (o.deep v)) (hx : g.IsEnd o.e x) (hxp : x ≠ p) :
    d.Anc c x := by
  rcases (isEnd_iff hs ho).mp hx with rfl | rfl
  · exact hd
  · exact (shallow_mem hs hc hp ho hd).resolve_right hxp

theorem ends_in_root {c v x : Nat} {o : DfsOut} (hr : d.IsRoot c) (ho : o ∈ d.outs v)
    (hd : d.Anc c (o.deep v)) (hx : g.IsEnd o.e x) : d.Anc c x := by
  rcases (isEnd_iff hs ho).mp hx with rfl | rfl
  · exact hd
  · exact shallow_mem_root hs hr ho hd

theorem ends_out {c v y : Nat} {o : DfsOut} (ho : o ∈ d.outs v) (hd : ¬d.Anc c (o.deep v))
    (hy : g.IsEnd o.e y) : ¬d.Anc c y := by
  rcases (isEnd_iff hs ho).mp hy with rfl | rfl
  · exact hd
  · exact fun h => hd (h.trans (shallow_anc_deep hs ho))

/-! ### (a) `component` / `bridge` children are block boundaries -/

/-- The parent `p` of a block root `c` separates every edge with deeper endpoint in `T_c` from
every edge whose deeper endpoint is not. -/
theorem blockRoot_cut {c p e e' : Nat} (hc : d.BlockRoot c) (hp : d.IsParent p c)
    (he : d.DeepIn c e) (he' : e' < g.ne) (he'' : ¬d.DeepIn c e') :
    ¬g.EdgeConn (· ≠ p) e e' := by
  obtain ⟨v, o, ho, rfl, hd⟩ := he
  obtain ⟨v', o', ho', rfl⟩ := hs.edge_out e' he'
  have hd' : ¬d.Anc c (o'.deep v') := fun h => he'' ⟨v', o', ho', rfl, h⟩
  rintro (h | ⟨x, y, hx, hy, hr⟩)
  · exact he'' ⟨v, o, ho, h, hd⟩
  · exact hc.sep hs hp (ends_in hs hc hp ho hd hx hr.ok_left) (ends_out hs ho' hd' hy) hr rfl

theorem root_cut {c e e' : Nat} {ok : Nat → Prop} (hr : d.IsRoot c) (he : d.DeepIn c e)
    (he' : e' < g.ne) (he'' : ¬d.DeepIn c e') : ¬g.EdgeConn ok e e' := by
  obtain ⟨v, o, ho, rfl, hd⟩ := he
  obtain ⟨v', o', ho', rfl⟩ := hs.edge_out e' he'
  have hd' : ¬d.Anc c (o'.deep v') := fun h => he'' ⟨v', o', ho', rfl, h⟩
  rintro (h | ⟨x, y, hx, hy, hr'⟩)
  · exact he'' ⟨v, o, ho, h, hd⟩
  · exact ends_out hs ho' hd' hy (hr.reach hs (ends_in_root hs hr ho hd hx) hr')

/-- The edges with deeper endpoint in `T_c`, `c` a block root, are a union of blocks. -/
theorem deepIn_of_sameBlock {c e e' : Nat} (hc : d.BlockRoot c) (he : d.DeepIn c e)
    (he' : e' < g.ne) (h : g.SameBlock e e') : d.DeepIn c e' := by
  by_contra he''
  rcases hc with hr | ⟨p, o, ho, ht, rfl, hcls⟩
  · exact root_cut hs hr he he' he'' (h 0)
  · exact blockRoot_cut hs (.inr ⟨p, o, ho, ht, rfl, hcls⟩) (isParent_of_tree ho ht) he he' he''
      (h p)

/-! ### (c) two edges of one block are not separated by any vertex -/

/-- For a non-block-root `x` with block top `c` and parent `q`, `T_x` has a back edge to a proper
ancestor `w` of `q`, and `w` is in the block (`BlockTop c w`) or is the parent of `c`. -/
theorem ret_up {c q x : Nat} (hx : d.BlockTop c x) (hne : x ≠ c) (hq : d.IsParent q x) :
    ∃ s w, d.Anc x s ∧ g.Adj s w ∧ d.Anc w q ∧ w ≠ q ∧ (d.BlockTop c w ∨ d.IsParent w c) := by
  obtain ⟨l, hl, hlt⟩ := returns_of_not_blockRoot hs hq (hx.not_blockRoot hs hne)
  obtain ⟨s, o, ho, hxs, hb, rfl⟩ := hl
  have hws := hs.back_anc s o ho hb
  have hwq : d.Anc o.dest q :=
    (hws.comparable hs (hq.anc.trans hxs)).resolve_right fun h => h.ne_of_depth_lt hs hlt
  have hwq' : o.dest ≠ q := fun h => by rw [h] at hlt; omega
  refine ⟨s, o.dest, hxs, (hs.joins s o ho).adj, hwq, hwq', ?_⟩
  have hcs : d.Anc c s := hx.2.1.trans hxs
  rcases hws.comparable hs hcs with hwc | hcw
  · by_cases hwc' : o.dest = c
    · exact .inl (by rw [hwc']; exact BlockTop.self hx.1)
    · rcases hwc.cases_tail with h | ⟨p, hwp, hpc⟩
      · exact absurd h.symm hwc'
      · have hret : d.Returns c (d.depth o.dest) := ⟨s, o, ho, hcs, hb, rfl⟩
        have hle := hwp.depth_le hs
        have hlt' : ¬d.depth o.dest < d.depth p := fun h => hx.1.no_returns hs hpc _ h hret
        exact .inr ((hwp.eq_of_depth_eq hs (.refl _) (by omega)) ▸ hpc)
  · exact .inl (hx.of_anc hcw (hwq.trans hq.anc))

/-- In the block topped by `c` (child of `p`), every vertex other than `v` reaches `c` avoiding
`v`, or `p` when `v = c`. -/
theorem reach_top {c p : Nat} (hp : d.IsParent p c) (v : Nat) :
    ∀ x, (d.BlockTop c x ∨ x = p) → x ≠ v → g.Reach (· ≠ v) x (if v = c then p else c) := by
  suffices ∀ n x, d.depth x = n → (d.BlockTop c x ∨ x = p) → x ≠ v →
      g.Reach (· ≠ v) x (if v = c then p else c) from fun x => this _ x rfl
  intro n
  induction n using Nat.strongRecOn with
  | _ n ih =>
    intro x hn hx hxv
    rcases hx with hx | rfl
    · by_cases hxc : x = c
      · subst hxc
        have hvc : ¬v = x := fun h => hxv h.symm
        simp only [hvc, ↓reduceIte]
        exact .refl hxv
      · obtain ⟨q, hqx, hq⟩ := hx.parent hxc
        have hdq : d.depth q < n := by rw [← hn, hs.depth_parent _ _ hqx]; omega
        by_cases hqv : q = v
        · subst hqv
          obtain ⟨s, w, hxs, hsw, hwq, hwq', hw⟩ := ret_up hs hx hxc hqx
          have hdw : d.depth w < n := by have := hwq.depth_le hs; omega
          have hw' : d.BlockTop c w ∨ w = p := hw.imp_right fun h => hs.parent_unique _ _ _ h hp
          have h1 : g.Reach (· ≠ q) x s := reach_in_subtree hs (.refl _) hxs fun z hz h =>
            hqx.ne hs (hqx.anc.antisymm hs (h ▸ hz))
          exact (h1.tail hsw hwq').trans (ih _ hdw w rfl hw' hwq')
        · exact (Graph.Reach.tail (.refl hxv : g.Reach (· ≠ v) x x) (hqx.adj hs).symm hqv).trans
            (ih _ hdq q rfl (.inl hq) hqv)
    · by_cases hvc : v = c
      · subst hvc; simp only [↓reduceIte]; exact .refl hxv
      · simp only [hvc, ↓reduceIte]
        exact Graph.Reach.tail (.refl hxv : g.Reach (· ≠ v) x x) (hp.adj hs) (Ne.symm hvc)

/-- A nontrivial edge of the block topped by `c` has an endpoint `≠ v` in the block. -/
theorem exists_end_ne {c p v w : Nat} {o : DfsOut} (hp : d.IsParent p c) (ho : o ∈ d.outs w)
    (hn : ¬o.isSingleton) (hd : d.BlockTop c (o.deep w)) :
    ∃ x, g.IsEnd o.e x ∧ x ≠ v ∧ (d.BlockTop c x ∨ x = p) := by
  by_cases h : o.deep w = v
  · refine ⟨o.shallow w, isEnd_shallow hs ho, fun h' => deep_ne_shallow hs ho hn (h.trans h'.symm),
      shallow_blockTop hs hp ho hd⟩
  · exact ⟨o.deep w, isEnd_deep hs ho, h, .inl hd⟩

theorem sameBlock_of_inBlock {c e e' : Nat} (h : d.InBlock c e) (h' : d.InBlock c e') :
    g.SameBlock e e' := by
  obtain ⟨a, o, ho, rfl, hn, hd⟩ := h
  obtain ⟨a', o', ho', rfl, hn', hd'⟩ := h'
  obtain ⟨p, hp⟩ : ∃ p, d.IsParent p c := by
    by_cases hc : o.deep a = c
    · rcases DfsOut.ends_cases a o with ⟨ht, h1, -⟩ | ⟨hb, h1, -⟩
      · exact ⟨a, hc ▸ h1 ▸ isParent_of_tree ho ht⟩
      · have hne := dest_ne_of_not_selfLoop hs ho hb fun h => hn (.inr h)
        rcases (hs.back_anc a o ho hb).cases_tail with h | ⟨q, -, hq⟩
        · exact absurd h.symm hne
        · exact ⟨q, (hc ▸ h1) ▸ hq⟩
    · exact exists_parent_of_blockTop hs hd hc
  intro v
  obtain ⟨x, hx, hxv, hxb⟩ := exists_end_ne hs hp (v := v) ho hn hd
  obtain ⟨x', hx', hxv', hxb'⟩ := exists_end_ne hs hp (v := v) ho' hn' hd'
  exact .of_reach hx hx' ((reach_top hs hp v x hxb hxv).trans (reach_top hs hp v x' hxb' hxv').symm)

/-! ### (c) edges of different blocks are separated -/

theorem not_sameBlock_of_outside {c a e' : Nat} {o : DfsOut} (ho : o ∈ d.outs a)
    (hd : d.BlockTop c (o.deep a)) (hout : ∀ y, g.IsEnd e' y → ¬d.Anc c y) :
    ¬g.SameBlock o.e e' := by
  intro h
  by_cases hr : d.IsRoot c
  · rcases h 0 with heq | ⟨x, y, hx, hy, hw⟩
    · exact hout _ (heq ▸ isEnd_deep hs ho) hd.2.1
    · exact hout y hy (hr.reach hs (ends_in_root hs hr ho hd.2.1 hx) hw)
  · obtain ⟨p, hp⟩ := IsRoot.exists_parent hr
    rcases h p with heq | ⟨x, y, hx, hy, hw⟩
    · exact hout _ (heq ▸ isEnd_deep hs ho) hd.2.1
    · exact hd.1.sep hs hp (ends_in hs hd.1 hp ho hd.2.1 hx hw.ok_left) (hout y hy) hw rfl

/-- A bridge or a loop is separated from every other edge. -/
theorem not_sameBlock_of_singleton {a e' : Nat} {o : DfsOut} (ho : o ∈ d.outs a)
    (hn : o.isSingleton) (he' : e' < g.ne) (hne : o.e ≠ e') : ¬g.SameBlock o.e e' := by
  obtain ⟨a', o', ho', rfl⟩ := hs.edge_out e' he'
  intro h
  rcases hn with hb | hl
  · obtain ⟨ht, hno⟩ := (hs.cls_bridge a o ho).mp hb
    have hp := isParent_of_tree ho ht
    have hc : d.BlockRoot o.dest := BlockRoot.of_tree ho ht (.inr hb)
    by_cases hin : d.Anc o.dest (o'.deep a')
    · have hall : ∀ y, g.IsEnd o'.e y → d.Anc o.dest y := by
        intro y hy
        rcases (isEnd_iff hs ho').mp hy with rfl | rfl
        · exact hin
        · rcases DfsOut.ends_cases a' o' with ⟨ht', h1, h2⟩ | ⟨hb', h1, h2⟩ <;>
            rw [h1] at hin <;> rw [h2]
          · by_contra ha'
            have hpar := isParent_of_tree ho' ht'
            have := IsParent.eq_of_anc hs hin hpar ha'
            rw [this] at hpar
            obtain rfl := hs.parent_unique _ _ _ hpar hp
            exact hne (congrArg DfsOut.e (hs.tree_inj _ o ho o' ho' ht ht' this.symm))
          · rcases (hs.back_anc a' o' ho' hb').comparable hs hin with hwc | hcw
            · by_contra hne'
              have hwp := hwc.to_parent hs (fun h => hne' (h ▸ .refl _)) hp
              exact hno _ (hwp.depth_le hs) ⟨a', o', ho', hin, hb', rfl⟩
            · exact hcw
      rcases h o.dest with heq | ⟨x, y, hx, hy, hw⟩
      · exact hne heq
      · have hx' : x = a := by
          rcases (isEnd_iff hs ho).mp hx with rfl | rfl
          · exact absurd (DfsOut.deep_tree ht) hw.ok_left
          · exact DfsOut.shallow_tree ht
        subst hx'
        exact bridge_sep hs ho hb (hall y hy) (hp.not_anc hs) hw.symm rfl
    · rcases h a with heq | ⟨x, y, hx, hy, hw⟩
      · exact hne heq
      · have hx' : d.Anc o.dest x := by
          rcases (isEnd_iff hs ho).mp hx with rfl | rfl
          · rw [DfsOut.deep_tree ht]; exact .refl _
          · exact absurd (DfsOut.shallow_tree ht) hw.ok_left
        exact hc.sep hs hp hx' (ends_out hs ho' hin hy) hw rfl
  · obtain ⟨hb, hd⟩ := (hs.cls_selfLoop a o ho).mp hl
    rcases h a with heq | ⟨x, y, hx, hy, hw⟩
    · exact hne heq
    · rcases (isEnd_iff hs ho).mp hx with rfl | rfl
      · exact hw.ok_left (DfsOut.deep_back hb)
      · exact hw.ok_left ((DfsOut.shallow_back hb).trans hd)

theorem inBlock_of_sameBlock {e e' : Nat} (he : e < g.ne) (he' : e' < g.ne) (hne : e ≠ e')
    (h : g.SameBlock e e') : ∃ c, d.InBlock c e ∧ d.InBlock c e' := by
  obtain ⟨a, o, ho, rfl⟩ := hs.edge_out e he
  obtain ⟨a', o', ho', rfl⟩ := hs.edge_out e' he'
  have hn : ¬o.isSingleton := fun hn => not_sameBlock_of_singleton hs ho hn he' hne h
  have hn' : ¬o'.isSingleton := fun hn =>
    not_sameBlock_of_singleton hs ho' hn he (Ne.symm hne) h.symm
  obtain ⟨c, hc⟩ := exists_blockTop hs (o.deep a)
  obtain ⟨c', hc'⟩ := exists_blockTop hs (o'.deep a')
  suffices c = c' from ⟨c, ⟨a, o, ho, rfl, hn, hc⟩, ⟨a', o', ho', rfl, hn', this ▸ hc'⟩⟩
  by_contra hcc
  by_cases hanc : d.Anc c c'
  · have hout : ¬d.Anc c' (o.deep a) := fun h => hcc ((hc.2.2 c' hc'.1 h).antisymm hs hanc).symm
    exact not_sameBlock_of_outside hs ho' hc' (fun y hy => ends_out hs ho hout hy) h.symm
  · have hout : ¬d.Anc c (o'.deep a') := fun h => hanc (hc'.2.2 c hc.1 h)
    exact not_sameBlock_of_outside hs ho hc (fun y hy => ends_out hs ho' hout hy) h

/-- Lemma 2.1: the blocks are the singletons `{bridge}`, `{loop}`, and for every block root `c`
the set of non-bridge non-loop edges whose deeper endpoint has block top `c`. -/
theorem sameBlock_iff {e e' : Nat} (he : e < g.ne) (he' : e' < g.ne) :
    g.SameBlock e e' ↔ e = e' ∨ ∃ c, d.InBlock c e ∧ d.InBlock c e' := by
  constructor
  · intro h
    by_cases hne : e = e'
    · exact .inl hne
    · exact .inr (inBlock_of_sameBlock hs he he' hne h)
  · rintro (rfl | ⟨c, h, h'⟩)
    · exact .refl e
    · exact sameBlock_of_inBlock hs h h'

theorem InBlock.unique {c c' e : Nat} (h : d.InBlock c e) (h' : d.InBlock c' e) : c = c' := by
  obtain ⟨a, o, ho, rfl, -, hc⟩ := h
  obtain ⟨a', o', ho', he, -, hc'⟩ := h'
  obtain ⟨rfl, rfl⟩ := out_eq_of_e hs ho' ho he
  exact hc.unique hs hc'

/-! ### (b) `ret` children stay in their parent's block -/

/-- The tree edge `v → c` of a `ret` child lies in the block of `v`'s block top. -/
theorem ret_child_inBlock {v l : Nat} {k : RetKind} {o : DfsOut} (ho : o ∈ d.outs v)
    (ht : o.isTree = true) (hc : o.cls = .ret l k) {c₀ : Nat} (hv : d.BlockTop c₀ v) :
    d.InBlock c₀ o.e := by
  have hp := isParent_of_tree ho ht
  refine ⟨v, o, ho, rfl, ?_, ?_⟩
  · rintro (h | h) <;> rw [hc] at h <;> cases h
  · rw [DfsOut.deep_tree ht]
    refine ⟨hv.1, hv.2.1.tail hp, fun c' hc' h => ?_⟩
    have hne : c' ≠ o.dest := fun h' => not_blockRoot_of_ret hs ho ht hc (h' ▸ hc')
    exact hv.2.2 c' hc' (h.to_parent hs hne hp)

/-- For a `ret` child `c` of `v`: every non-bridge non-loop edge with deeper endpoint in `T_c`
not inside a nested block is in the block of the tree edge `v → c`. -/
theorem ret_child_sameBlock {v l : Nat} {k : RetKind} {o : DfsOut} (ho : o ∈ d.outs v)
    (ht : o.isTree = true) (hc : o.cls = .ret l k) {a' : Nat} {o' : DfsOut} (ho' : o' ∈ d.outs a')
    (hn' : ¬o'.isSingleton) (hd' : d.Anc o.dest (o'.deep a'))
    (hnest : ∀ c', d.BlockRoot c' → d.Anc c' (o'.deep a') → d.Anc c' o.dest) :
    g.SameBlock o.e o'.e := by
  obtain ⟨c₀, hv⟩ := exists_blockTop hs v
  have hp := isParent_of_tree ho ht
  refine sameBlock_of_inBlock hs (ret_child_inBlock hs ho ht hc hv) ⟨a', o', ho', rfl, hn', ?_⟩
  refine ⟨hv.1, (hv.2.1.tail hp).trans hd', fun c' hc' h => ?_⟩
  have hne : c' ≠ o.dest := fun h' => not_blockRoot_of_ret hs ho ht hc (h' ▸ hc')
  exact hv.2.2 c' hc' ((hnest c' hc' h).to_parent hs hne hp)

/-- For a `ret l _` child `c` of `v`: the back edges out of `T_c` returning to depth `l`
(`lowpt1 c`) are in the block of the tree edge `v → c`. -/
theorem ret_back_sameBlock {v l : Nat} {k : RetKind} {o : DfsOut} (ho : o ∈ d.outs v)
    (ht : o.isTree = true) (hc : o.cls = .ret l k) {u : Nat} {o' : DfsOut} (ho' : o' ∈ d.outs u)
    (hb' : o'.isTree = false) (hu : d.Anc o.dest u) (hl : d.depth o'.dest = l) :
    g.SameBlock o.e o'.e := by
  obtain ⟨hlt, -, -⟩ := ret_lowpt hs ho ht hc
  have hp := isParent_of_tree ho ht
  have hdu : d.depth v < d.depth u := by
    have := hu.depth_le hs; rw [hs.depth_parent _ _ hp] at this; omega
  refine ret_child_sameBlock hs ho ht hc ho' ?_ ?_ ?_
  · rintro (h | h)
    · exact not_bridge_of_back hs ho' hb' h
    · have := ((hs.cls_selfLoop u o' ho').mp h).2; rw [this] at hl; omega
  · rwa [DfsOut.deep_back hb']
  · rw [DfsOut.deep_back hb']
    intro c' hc' hc'u
    rcases hc'u.comparable hs hu with h | h
    · exact h
    · by_contra hne
      rcases h.cases_tail with h' | ⟨q, hq, hqc'⟩
      · exact hne (h' ▸ .refl _)
      · have hret : d.Returns c' l := ⟨u, o', ho', hc'u, hb', hl⟩
        have := hq.depth_le hs
        have := hs.depth_parent _ _ hp
        exact hc'.no_returns hs hqc' l (by omega) hret

end Spec

end DfsData

end Spqr
