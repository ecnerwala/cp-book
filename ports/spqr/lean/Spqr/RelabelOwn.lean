import Mathlib.Data.List.Sort
import Mathlib.Data.List.Nodup
import Spqr.RelabelSpec
import Spqr.ItemTree
import Spqr.LayoutShape

/-!
# Phase 3b: `Bijections`, `Ownership`, `Twins` of `relabelTree` from the per-node interface

Everything here is derived from `relabel_node_spec` (taken as a hypothesis through its `∃ idx`)
and `Items.WF`, plus the extra item-level hypotheses collected in `Items.OwnExtra` that
`Items.WF` does not promise (see its docstring). The one `Layout`-local fact needed, the edge
records of `layoutNode` (`layoutNode_edges`), is read off the per-type forms of `LayoutShape.lean`;
only the R edges are re-derived here without the endpoint hypothesis of `LayoutShape.run_edges`,
so that `Twins` does not depend on `ROriented`.
-/

namespace Spqr

/-- `omega` after unfolding the `ItemId` abbreviation (`omega` does not see through it). -/
macro "iomega" : tactic => `(tactic| ((try unfold ItemId at *); (try delta ItemId at *); omega))

/-! ### Generic array / list helpers -/

theorem Array.getElem!_eq_getD_getElem? {α : Type} [Inhabited α] (a : Array α) (i : Nat) :
    a[i]! = a[i]?.getD default := by
  rw [Array.getElem!_eq_getD, Array.getD_eq_getD_getElem?]

theorem Array.getElem!_of_lt {α : Type} [Inhabited α] (a : Array α) (i : Nat) (h : i < a.size) :
    a[i]? = some a[i]! := by
  rw [Array.getElem!_eq_getD_getElem?, Array.getElem?_eq_getElem h]; rfl

/-- A monotone chain of bounds. -/
theorem bounds_le_of_le (b : Array Nat)
    (hmono : ∀ i, i + 1 < b.size → b[i]! ≤ b[i + 1]!) :
    ∀ i j, i ≤ j → j < b.size → b[i]! ≤ b[j]! := by
  intro i j hij hj
  induction j, hij using Nat.le_induction with
  | base => exact le_rfl
  | succ j _ ih => exact (ih (by iomega)).trans (hmono j hj)

/-- Locate the interval `[b[i], b[i+1])` containing `x`. -/
theorem bounds_locate (b : Array Nat) (h0 : b[0]! = 0) :
    ∀ n, n < b.size → ∀ x, x < b[n]! → ∃ i, i < n ∧ b[i]! ≤ x ∧ x < b[i + 1]! := by
  intro n
  induction n with
  | zero => intro _ x hx; iomega
  | succ n ih =>
    intro hn x hx
    by_cases h : x < b[n]!
    · obtain ⟨i, hi, h1, h2⟩ := ih (by iomega) x h
      exact ⟨i, by iomega, h1, h2⟩
    · exact ⟨n, by iomega, by iomega, hx⟩

theorem Option.toList_length_le {α : Type} (o : Option α) : o.toList.length ≤ 1 := by
  cases o <;> simp

theorem List.getElem?_getD_eq_getElem! {α : Type} [Inhabited α] (a : Array α) (i : Nat) :
    a[i]?.getD default = a[i]! := (Array.getElem!_eq_getD_getElem? a i).symm

/-! ### Item typing from the id layout -/

namespace Items

variable {g : Graph} {items : Items}

theorem Tree.type_V_iff (ht : items.Tree g) {i : ItemId} (hi : i < items.size) :
    items.type i = .V ↔ 1 ≤ i ∧ i < 1 + g.nv := by
  have hroot : items.type 0 = .F := ht.root
  constructor
  · intro h
    rcases Nat.lt_or_ge i 1 with h1 | h1
    · have : i = 0 := by iomega
      subst this; rw [hroot] at h; cases h
    rcases Nat.lt_or_ge i (1 + g.nv) with h2 | h2
    · exact ⟨h1, h2⟩
    rcases Nat.lt_or_ge i (1 + g.nv + g.ne) with h3 | h3
    · have := ht.edge (i - 1 - g.nv) (by iomega)
      have he : edgeItem g (i - 1 - g.nv) = i := by simp only [edgeItem]; iomega
      rw [he, h] at this; cases this
    · have := ht.node i h3 hi
      simp [h] at this
  · rintro ⟨h1, h2⟩
    have := ht.vert (i - 1) (by iomega)
    have hv : vertItem (i - 1) = i := by simp only [vertItem]; iomega
    rwa [hv] at this

theorem Tree.type_Q_iff (ht : items.Tree g) {i : ItemId} (hi : i < items.size) :
    items.type i = .Q ↔ 1 + g.nv ≤ i ∧ i < 1 + g.nv + g.ne := by
  have hroot : items.type 0 = .F := ht.root
  constructor
  · intro h
    rcases Nat.lt_or_ge i 1 with h1 | h1
    · have : i = 0 := by iomega
      subst this; rw [hroot] at h; cases h
    rcases Nat.lt_or_ge i (1 + g.nv) with h2 | h2
    · have := ht.vert (i - 1) (by iomega)
      have hv : vertItem (i - 1) = i := by simp only [vertItem]; iomega
      rw [hv, h] at this; cases this
    rcases Nat.lt_or_ge i (1 + g.nv + g.ne) with h3 | h3
    · exact ⟨h2, h3⟩
    · have := ht.node i h3 hi
      simp [h] at this
  · rintro ⟨h1, h2⟩
    have := ht.edge (i - 1 - g.nv) (by iomega)
    have he : edgeItem g (i - 1 - g.nv) = i := by simp only [edgeItem]; iomega
    rwa [he] at this

theorem Tree.type_F_iff (ht : items.Tree g) {i : ItemId} (hi : i < items.size) :
    items.type i = .F ↔ i = 0 := by
  constructor
  · intro h
    rcases Nat.lt_or_ge i 1 with h1 | h1
    · iomega
    rcases Nat.lt_or_ge i (1 + g.nv) with h2 | h2
    · have := (ht.type_V_iff hi).2 ⟨by iomega, h2⟩
      rw [h] at this; cases this
    rcases Nat.lt_or_ge i (1 + g.nv + g.ne) with h3 | h3
    · have := (ht.type_Q_iff hi).2 ⟨h2, h3⟩
      rw [h] at this; cases this
    · have := ht.node i h3 hi
      simp [h] at this
  · rintro rfl; exact ht.root

theorem Tree.isNode_iff (ht : items.Tree g) {i : ItemId} (hi : i < items.size) :
    (items.type i).isNode = true ↔ 1 + g.nv ≤ i := by
  have hV := ht.type_V_iff hi
  have hF := ht.type_F_iff hi
  cases h : items.type i <;> simp [NodeType.isNode, h] at hV hF ⊢ <;> iomega

theorem Tree.child_lt (ht : items.Tree g) {p c : ItemId} (hc : c ∈ items.ch p) : c < items.size :=
  ht.ch_lt p c hc

theorem Tree.child_pos (ht : items.Tree g) {p c : ItemId} (hc : c ∈ items.ch p) : 0 < c := by
  rcases Nat.eq_zero_or_pos c with h | h
  · subst h; exact absurd hc (ht.root_no_parent p)
  · exact h

theorem Tree.child_ne (ht : items.Tree g) {p c : ItemId} (hc : c ∈ items.ch p) : c ≠ p := by
  intro h; subst h; exact ht.acyclic c (Relation.TransGen.single hc)

/-- The parent of a non-root item. -/
theorem Tree.parent_exists (ht : items.Tree g) {c : ItemId} (hc : c < items.size) (h0 : c ≠ 0) :
    ∃ p, c ∈ items.ch p := by
  obtain ⟨p, hp, _⟩ := ht.unique_parent c (by iomega) hc
  exact ⟨p, hp⟩

theorem Tree.parent_lt (_ht : items.Tree g) {p c : ItemId} (hc : c ∈ items.ch p) : p < items.size := by
  by_contra h
  have : items.ch p = [] := by
    simp [Items.ch, Array.getElem?_eq_none_iff.2 (show items.size ≤ p by iomega)]
  rw [this] at hc; cases hc

/-! ### Node-vert lists -/

theorem nvList_length (i : ItemId) :
    (items.nvList g i).length =
      (items.vs i).1.toList.length + ((items.ch i).filter (· < 1 + g.nv)).length +
        (items.vs i).2.toList.length := by
  simp [Items.nvList, Nat.add_assoc]

theorem nvList_mem_lt (ht : items.Tree g) (he : items.Endpoints g) {i : ItemId} {v : Nat}
    (hv : v ∈ items.nvList g i) : v < g.nv := by
  simp only [Items.nvList, List.mem_append, List.mem_map, List.mem_filter] at hv
  rcases hv with (hv | ⟨c, ⟨hc, hlt⟩, rfl⟩) | hv
  · exact he.vs_lt i v (Or.inl (by cases h : (items.vs i).1 <;> simp [h] at hv ⊢; exact hv.symm))
  · have h1 : 0 < c := ht.child_pos hc
    have h2 : c < 1 + g.nv := by simpa using hlt
    iomega
  · exact he.vs_lt i v (Or.inr (by cases h : (items.vs i).2 <;> simp [h] at hv ⊢; exact hv.symm))


/-! ### Extra item-level hypotheses -/

/-- What `Ownership`/`Twins` need beyond `Items.WF` (all true of the walk's output; none is
promised by `Items.WF`):
* `nv_nodup`: a node's node-vert list has no repeats (`WF` allows `I`/`P` nodes with `u = v`);
* `q_leaf_of_node`: a Q child of a node is a leaf (block-root Qs hang under V items);
* `r_edges_in_nv`: both endpoints of an R node's virtual edges are node-verts (`WF` only
  constrains `some` endpoints, so a loop-Q/O child would contribute `getD 0`);
* `q_vert_child`: the `vertItem v` child of a block-root Q is a real V item (`Shapes.q_children`
  does not bound `v`). -/
structure OwnExtra (g : Graph) (items : Items) : Prop where
  nv_nodup : ∀ i, i < items.size → (items.nvList g i).Nodup
  q_leaf_of_node : ∀ p c, items.IsParent p c → (items.type p).isNode = true →
    items.type c = .Q → items.ch c = []
  r_edges_in_nv : ∀ i, i < items.size → items.type i = .R →
    ∀ p ∈ items.virtualEdges i, p.1 ∈ items.nvList g i ∧ p.2 ∈ items.nvList g i
  q_vert_child : ∀ i, i < items.size → items.type i = .Q →
    ∀ v, vertItem v ∈ items.ch i → v < g.nv

/-! ### Children filters -/

theorem Tree.filter_V_eq (ht : items.Tree g) (i : ItemId) :
    ((items.ch i).filter fun c => items.type c = .V) = (items.ch i).filter (· < 1 + g.nv) := by
  apply List.filter_congr
  intro c hc
  have h1 := ht.child_pos hc
  rw [decide_eq_decide, ht.type_V_iff (ht.child_lt hc)]
  constructor
  · exact fun h => h.2
  · intro h; constructor
    · iomega
    · exact h

theorem Tree.filter_nonV_eq (ht : items.Tree g) (i : ItemId) :
    ((items.ch i).filter fun c => items.type c ≠ .V) = (items.ch i).filter (· ≥ 1 + g.nv) := by
  apply List.filter_congr
  intro c hc
  have h1 := ht.child_pos hc
  rw [decide_eq_decide, ne_eq, ht.type_V_iff (ht.child_lt hc)]
  constructor
  · intro h; by_contra h'; apply h; constructor
    · iomega
    · iomega
  · intro h ⟨_, h2⟩; iomega

theorem Tree.countP_nonV (ht : items.Tree g) (i : ItemId) :
    (items.ch i).countP (· ≥ 1 + g.nv) = (items.virtualEdges i).length := by
  rw [List.countP_eq_length_filter, Items.virtualEdges, List.length_map, ht.filter_nonV_eq]

theorem Tree.child_ge_of_ne_V (ht : items.Tree g) {p c : ItemId} (hc : c ∈ items.ch p)
    (h : items.type c ≠ .V) : 1 + g.nv ≤ c := by
  have h1 := ht.child_pos hc
  by_contra h'
  apply h
  rw [ht.type_V_iff (ht.child_lt hc)]
  iomega

theorem Tree.child_lt_of_V (ht : items.Tree g) {p c : ItemId} (hc : c ∈ items.ch p)
    (h : items.type c = .V) : c < 1 + g.nv := ((ht.type_V_iff (ht.child_lt hc)).1 h).2

theorem Tree.type_ne_F_of_child (ht : items.Tree g) {p c : ItemId} (hc : c ∈ items.ch p) :
    items.type c ≠ .F := by
  intro h
  have := (ht.type_F_iff (ht.child_lt hc)).1 h
  have := ht.child_pos hc
  iomega

theorem Tree.isNode_of_ge (ht : items.Tree g) {c : ItemId} (hc : c < items.size)
    (h : 1 + g.nv ≤ c) : (items.type c).isNode = true := (ht.isNode_iff hc).2 h

/-! ### `ordered` -/

theorem ordered_perm (i : ItemId) (nvSt : Nat) (pos : Nat → Nat) :
    (items.ordered g i nvSt pos).Perm (items.ch i) := by
  unfold Items.ordered
  split
  · exact List.Perm.refl _
  · exact List.mergeSort_perm _ _

theorem ordered_eq_of_ne_R {i : ItemId} (h : items.type i ≠ .R) (nvSt : Nat) (pos : Nat → Nat) :
    items.ordered g i nvSt pos = items.ch i := by
  unfold Items.ordered; simp [h]

theorem mem_ordered {i c : ItemId} {nvSt : Nat} {pos : Nat → Nat} :
    c ∈ items.ordered g i nvSt pos ↔ c ∈ items.ch i := (ordered_perm i nvSt pos).mem_iff

/-! ### Per-type counting facts -/

theorem nEdges_eq {i : ItemId} (hn : (items.type i).isNode = true) :
    items.nEdges g i = (items.ch i).countP (· ≥ 1 + g.nv) + items.capCount i := by
  simp [Items.nEdges, hn]

theorem nEdges_eq_zero {i : ItemId} (hn : (items.type i).isNode = false) :
    items.nEdges g i = 0 := by
  simp [Items.nEdges, Items.capCount, Items.hasCap, hn]

theorem hasCap_of_ne_Q {i : ItemId} (hn : (items.type i).isNode = true) (hq : items.type i ≠ .Q) :
    items.hasCap i = true := by
  simp [Items.hasCap, hn, hq]

theorem hasCap_Q {i : ItemId} (hq : items.type i = .Q) :
    items.hasCap i = (items.ch i).isEmpty := by
  simp [Items.hasCap, hq, NodeType.isNode]

theorem capCount_le (i : ItemId) : items.capCount i ≤ 1 := by
  unfold Items.capCount; split <;> omega

theorem WF.vs_fst (hw : items.WF g) {i : ItemId} (hi : i < items.size)
    (hn : (items.type i).isNode = true) : ∃ u, (items.vs i).1 = some u := by
  have h := hw.endpoints.vs_shape i hi
  cases ht : items.type i <;> simp only [ht] at h hn <;> simp [NodeType.isNode] at hn
  · rcases h with ⟨v, hv⟩ | ⟨u, v, hv⟩ <;> simp [hv]
  · obtain ⟨u, v, hv⟩ := h; simp [hv]
  · obtain ⟨v, hv⟩ := h; simp [hv]
  · obtain ⟨u, v, hv⟩ := h; simp [hv]
  · obtain ⟨u, v, hv⟩ := h; simp [hv]
  · obtain ⟨u, v, hv⟩ := h; simp [hv]

theorem WF.vs_two (hw : items.WF g) {i : ItemId} (hi : i < items.size)
    (hn : items.type i = .I ∨ items.type i = .S ∨ items.type i = .P ∨ items.type i = .R) :
    ∃ u v, items.vs i = (some u, some v) := by
  have h := hw.endpoints.vs_shape i hi
  rcases hn with ht | ht | ht | ht <;> simp only [ht] at h <;> exact h

theorem WF.q_cases (hw : items.WF g) {i : ItemId} (hi : i < items.size) (hQ : items.type i = .Q) :
    items.ch i = [] ∨ ∃ c, items.type c ∉ [NodeType.F, .V, .Q] ∧
      (items.ch i = [c] ∨ ∃ v, items.ch i = [c, vertItem v]) := by
  have hr := (hw.tree.type_Q_iff hi).1 hQ
  have he : edgeItem g (i - 1 - g.nv) = i := by simp only [edgeItem]; iomega
  have := hw.shapes.q_children (i - 1 - g.nv) (by iomega)
  rwa [he] at this

theorem WF.one_le_nvList (hw : items.WF g) {i : ItemId} (hi : i < items.size)
    (hn : (items.type i).isNode = true) : 1 ≤ (items.nvList g i).length := by
  obtain ⟨u, hu⟩ := hw.vs_fst hi hn
  rw [nvList_length, hu]; simp only [Option.toList_some, List.length_singleton]; omega

/-- The hypotheses of `layoutNode_edges` and the two-vertex facts, per node type. -/
theorem WF.layout_hyps (hw : items.WF g) (hx : items.OwnExtra g) {i : ItemId} (hi : i < items.size)
    (hn : (items.type i).isNode = true) :
    ((items.nvList g i).length = 1 → items.nEdges g i ≤ 1) ∧
    (items.type i = .Q ∨ items.type i = .I → items.nEdges g i ≤ 1) ∧
    (items.type i = .S → items.nEdges g i = (items.nvList g i).length) ∧
    (items.type i = .R → 4 ≤ (items.nvList g i).length) ∧
    (items.hasCap i = true → items.type i = .Q ∨ items.type i = .I ∨ items.type i = .P →
      (items.nvList g i).length ≤ 2) ∧
    (items.type i = .O → (items.nvList g i).length = 1) := by
  have ht := hw.tree
  have hlen := nvList_length (g := g) (items := items) i
  have hne := nEdges_eq (g := g) hn
  have hcap := capCount_le (items := items) i
  cases hty : items.type i
  all_goals simp only [hty] at hn hne ⊢
  all_goals simp [NodeType.isNode] at hn
  · -- Q
    rcases hw.q_cases hi hty with h0 | ⟨c, hc, h1 | ⟨v, h1⟩⟩
    · simp only [h0, List.countP_nil, List.filter_nil, List.length_nil] at hne hlen
      simp only [Items.hasCap, hty] at hcap ⊢
      refine ⟨fun _ => by omega, fun _ => by omega, by simp, by simp, fun _ _ => ?_, by simp⟩
      have := Option.toList_length_le (items.vs i).1
      have := Option.toList_length_le (items.vs i).2
      omega
    · have hcge : 1 + g.nv ≤ c := ht.child_ge_of_ne_V (by rw [h1]; simp) (by simp at hc; tauto)
      have hcapv : items.capCount i = 0 := by
        simp [Items.capCount, Items.hasCap, hty, h1]
      simp only [h1, List.countP_cons, List.countP_nil] at hne
      simp [hcge, hcapv] at hne
      refine ⟨fun _ => by omega, fun _ => by omega, by simp, by simp, ?_, by simp⟩
      intro hcp; simp [Items.hasCap, hty, h1] at hcp
    · have hcge : 1 + g.nv ≤ c := ht.child_ge_of_ne_V (by rw [h1]; simp) (by simp at hc; tauto)
      have hcapv : items.capCount i = 0 := by
        simp [Items.capCount, Items.hasCap, hty, h1]
      have hv : v < g.nv := hx.q_vert_child i hi hty v (by rw [h1]; simp)
      have hv' : ¬ g.nv ≤ v := by omega
      simp only [h1, List.countP_cons, List.countP_nil] at hne
      simp [hcge, hcapv, vertItem, hv'] at hne
      simp only [h1] at hlen
      simp [vertItem, hcge] at hlen
      have := Option.toList_length_le (items.vs i).1
      have := hw.one_le_nvList hi (by simp [hty, NodeType.isNode])
      refine ⟨fun h => by omega, fun _ => by omega, by simp, by simp, ?_, by simp⟩
      intro hcp; simp [Items.hasCap, hty, h1] at hcp
  · -- I
    have h0 := hw.shapes.i_o_leaf i hi (Or.inl hty)
    obtain ⟨u, v, huv⟩ := hw.vs_two hi (Or.inl hty)
    simp only [h0, List.countP_nil, List.filter_nil, List.length_nil, huv, Option.toList_some,
      List.length_singleton] at hne hlen
    refine ⟨fun _ => by omega, fun _ => by omega, by simp, by simp, fun _ _ => by omega, by simp⟩
  · -- O
    have h0 := hw.shapes.i_o_leaf i hi (Or.inr hty)
    have hvs := hw.endpoints.vs_shape i hi
    simp only [hty] at hvs
    obtain ⟨v, hv⟩ := hvs
    simp only [h0, List.countP_nil, List.filter_nil, List.length_nil, hv, Option.toList_some,
      Option.toList_none, List.length_singleton] at hne hlen
    refine ⟨fun _ => by omega, by simp, by simp, by simp, by simp, fun _ => by omega⟩
  · -- S
    obtain ⟨u, v, xs, huv, hxs, hk, hperm⟩ := hw.shapes.s_shape i hi hty
    have hcnt := ht.countP_nonV (g := g) i
    rw [hperm.length_eq, List.length_zip] at hcnt
    simp only [List.length_cons, List.length_append] at hcnt
    have hV : ((items.ch i).filter (· < 1 + g.nv)).length = xs.length := by
      rw [← ht.filter_V_eq, hxs, List.length_map]
    have hc1 : items.capCount i = 1 := by
      simp [Items.capCount, Items.hasCap, hty, NodeType.isNode]
    rw [huv] at hlen; simp only [Option.toList_some, List.length_singleton] at hlen
    refine ⟨fun _ => by omega, by simp, fun _ => by omega, by simp, by simp, by simp⟩
  · -- P
    obtain ⟨_, hnoV, _⟩ := hw.shapes.p_shape i hi hty
    obtain ⟨u, v, huv⟩ := hw.vs_two hi (Or.inr (Or.inr (Or.inl hty)))
    have hV : ((items.ch i).filter (· < 1 + g.nv)) = [] := by
      rw [← ht.filter_V_eq, List.filter_eq_nil_iff]
      intro c hc; simpa using hnoV c hc
    rw [huv, hV] at hlen; simp only [Option.toList_some, List.length_singleton, List.length_nil] at hlen
    refine ⟨fun _ => by omega, by simp, by simp, by simp, fun _ _ => by omega, by simp⟩
  · -- R
    obtain ⟨h2, _, _, _⟩ := hw.shapes.r_shape i hi hty
    obtain ⟨u, v, huv⟩ := hw.vs_two hi (Or.inr (Or.inr (Or.inr hty)))
    rw [ht.filter_V_eq] at h2
    rw [huv] at hlen; simp only [Option.toList_some, List.length_singleton] at hlen
    refine ⟨fun _ => by omega, by simp, by simp, fun _ => by omega, by simp, by simp⟩

/-! ### Positions of V children; `ordered` keeps the V children in `ch` order -/

theorem nvList_getElem?_mid (i : ItemId) {p : Nat}
    (hp : p < ((items.ch i).filter (· < 1 + g.nv)).length) :
    (items.nvList g i)[(items.vs i).1.toList.length + p]? =
      some (((items.ch i).filter (· < 1 + g.nv))[p] - 1) := by
  unfold Items.nvList
  rw [List.getElem?_append_left (by simp; omega), List.getElem?_append_right (by simp),
    Nat.add_sub_cancel_left, List.getElem?_map, List.getElem?_eq_getElem hp, Option.map_some]

theorem idxOf_eq_of_getElem? {l : List Nat} (hnd : l.Nodup) {k u : Nat} (h : l[k]? = some u) :
    l.idxOf u = k := by
  obtain ⟨hk, hu⟩ := List.getElem?_eq_some_iff.1 h
  have hmem : u ∈ l := hu ▸ List.getElem_mem hk
  have hlt : l.idxOf u < l.length := List.idxOf_lt_length_iff.2 hmem
  exact (hnd.getElem_inj_iff (hi := hlt) (hj := hk)).1 ((List.getElem_idxOf hlt).trans hu.symm)

theorem PosOK.pos_sub (hnd : (items.nvList g i).Nodup)
    (hpos : Items.PosOK nvSt (items.nvList g i) pos) {v : Nat} (hv : v ∈ items.nvList g i) :
    pos v - nvSt = (items.nvList g i).idxOf v :=
  (idxOf_eq_of_getElem? hnd (hpos v hv).2).symm

theorem PosOK.pos_lt (hpos : Items.PosOK nvSt (items.nvList g i) pos) {v : Nat}
    (hv : v ∈ items.nvList g i) : pos v - nvSt < (items.nvList g i).length :=
  (List.getElem?_eq_some_iff.1 (hpos v hv).2).1

/-- Under `PosOK`, the `p`-th V child sits at node-vert `nvSt + |vs.1| + p`. -/
theorem pos_V_child (hnd : (items.nvList g i).Nodup)
    (hpos : Items.PosOK nvSt (items.nvList g i) pos) {p : Nat}
    (hp : p < ((items.ch i).filter (· < 1 + g.nv)).length) :
    pos (((items.ch i).filter (· < 1 + g.nv))[p] - 1) - nvSt = (items.vs i).1.toList.length + p := by
  have h := nvList_getElem?_mid (g := g) (items := items) i hp
  have hmem : ((items.ch i).filter (· < 1 + g.nv))[p] - 1 ∈ items.nvList g i :=
    List.mem_iff_getElem?.2 ⟨_, h⟩
  rw [PosOK.pos_sub hnd hpos hmem, idxOf_eq_of_getElem? hnd h]

theorem loc_V {c : ItemId} (hc : c < 1 + g.nv) (nvSt : Nat) (pos : Nat → Nat) :
    items.loc g nvSt pos c = 2 * (pos (c - 1) - nvSt) := by
  simp [Items.loc, hc]

theorem ordered_filter_V (ht : items.Tree g) {i : ItemId} (hnd : (items.nvList g i).Nodup)
    {nvSt : Nat} {pos : Nat → Nat}
    (hpos : items.type i = .R → Items.PosOK nvSt (items.nvList g i) pos) :
    (items.ordered g i nvSt pos).filter (· < 1 + g.nv) = (items.ch i).filter (· < 1 + g.nv) := by
  by_cases hR : items.type i = .R
  swap
  · rw [ordered_eq_of_ne_R hR]
  have hpos := hpos hR
  have hperm : ((items.ordered g i nvSt pos).filter (· < 1 + g.nv)).Perm
      ((items.ch i).filter (· < 1 + g.nv)) := (ordered_perm i nvSt pos).filter _
  have hLpw : ((items.ch i).filter (· < 1 + g.nv)).Pairwise
      fun a b => items.loc g nvSt pos a < items.loc g nvSt pos b := by
    rw [List.pairwise_iff_getElem]
    intro p q hp hq hpq
    have ha := (List.mem_filter.1 (List.getElem_mem hp)).2
    have hb := (List.mem_filter.1 (List.getElem_mem hq)).2
    rw [loc_V (by simpa using ha), loc_V (by simpa using hb), pos_V_child hnd hpos hp,
      pos_V_child hnd hpos hq]
    omega
  have hkey : ∀ a ∈ (items.ch i).filter (· < 1 + g.nv), ∀ b ∈ (items.ch i).filter (· < 1 + g.nv),
      a ≠ b → items.loc g nvSt pos a ≠ items.loc g nvSt pos b := by
    intro a ha b hb hab
    obtain ⟨p, hp, rfl⟩ := List.mem_iff_getElem.1 ha
    obtain ⟨q, hq, rfl⟩ := List.mem_iff_getElem.1 hb
    rcases lt_trichotomy p q with h | h | h
    · exact Nat.ne_of_lt (List.pairwise_iff_getElem.1 hLpw p q hp hq h)
    · subst h; exact absurd rfl hab
    · exact Nat.ne_of_gt (List.pairwise_iff_getElem.1 hLpw q p hq hp h)
  have hMnd : ((items.ordered g i nvSt pos).filter (· < 1 + g.nv)).Nodup :=
    List.filter_sublist.nodup ((ordered_perm i nvSt pos).nodup_iff.2 (ht.ch_nodup i))
  have hMle : ((items.ordered g i nvSt pos).filter (· < 1 + g.nv)).Pairwise
      fun a b => items.loc g nvSt pos a ≤ items.loc g nvSt pos b := by
    have h2 : items.ordered g i nvSt pos =
        (items.ch i).mergeSort fun a b => decide (items.loc g nvSt pos a ≤ items.loc g nvSt pos b) := by
      unfold Items.ordered; simp [hR]
    rw [h2]
    refine ((List.pairwise_mergeSort (fun a b c => ?_) (fun a b => ?_) (items.ch i)).filter _).imp ?_
    · simp only [decide_eq_true_eq]; omega
    · simp only [Bool.or_eq_true, decide_eq_true_eq]; omega
    · simp only [decide_eq_true_eq]; exact id
  have hMpw : ((items.ordered g i nvSt pos).filter (· < 1 + g.nv)).Pairwise
      fun a b => items.loc g nvSt pos a < items.loc g nvSt pos b := by
    refine (hMle.and hMnd).imp_of_mem ?_
    intro a b ha hb ⟨hle, hne⟩
    exact lt_of_le_of_ne hle (hkey a (hperm.mem_iff.1 ha) b (hperm.mem_iff.1 hb) hne)
  exact List.Perm.eq_of_pairwise (fun a b _ _ h1 h2 => absurd h2 (not_lt.2 h1.le)) hMpw hLpw hperm

theorem edgeChildren_length (i : ItemId) (nvSt : Nat) (pos : Nat → Nat) :
    (items.edgeChildren g pos (items.ordered g i nvSt pos)).length =
      (items.ch i).countP (· ≥ 1 + g.nv) := by
  unfold Items.edgeChildren
  rw [List.length_map, ← List.countP_eq_length_filter, (ordered_perm i nvSt pos).countP_eq]

end Items

/-! ### `layoutNode` edges (a `Layout`-local fact, in the scope of `LayoutShape`) -/

/-- Endpoints (node-vert ids) of edge `k` of `layoutNode`. -/
def layoutNvs (ty : NodeType) (nvSt nvEn : Nat) (ec : List (Nat × Nat)) (k : Nat) : Nat × Nat :=
  if nvEn - nvSt = 1 then (nvSt, nvSt)
  else if ty = .Q ∨ ty = .I ∨ ty = .P then (nvSt, nvSt + 1)
  else if ty = .S then (if k = 0 then (nvSt, nvEn - 1) else (nvSt + k - 1, nvSt + k))
  else (if k = 0 then (nvSt, nvEn - 1) else ec[k - 1]?.getD (0, 0))

namespace LayoutEdges

open LayoutR LayoutShape

theorem inc_edges (nvSt : Nat) (l : Layout) (i : Nat) : (inc nvSt l i).edges = l.edges := rfl

theorem foldl_countStep_edges (nvSt : Nat) (P : List (Nat × Nat)) :
    ∀ l : Layout, (P.foldl (countStep nvSt) l).edges = l.edges := by
  induction P with
  | nil => intro l; rfl
  | cons p P ih => intro l; rw [List.foldl_cons, ih, countStep_edges]

theorem foldl_prefixStep_edges (nvSt : Nat) (L : List Nat) :
    ∀ s : Layout × Nat, (L.foldl (prefixStep nvSt) s).1.edges = s.1.edges := by
  induction L with
  | nil => intro s; rfl
  | cons i L ih => intro s; rw [List.foldl_cons, ih]; rfl

/-- The layout after the counting and prefix-sum passes of `run`. -/
def prefixed (nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat)) : Layout :=
  ((List.range' (2 * nvSt + 1) (2 * nvEn + 1 - (2 * nvSt + 1))).foldl (prefixStep nvSt)
    (E.foldl (countStep nvSt) (inc nvSt (inc nvSt (Layout.empty (nvEn - nvSt) (neEn - neSt))
      (2 * nvSt + 2)) (2 * nvEn - 1)), 2 * neSt)).1

theorem run_eq (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat)) :
    run node nvSt nvEn neSt neEn E =
      inc nvSt ((E.reverse.foldl (fillStep nvSt neSt node)
        (inc nvSt (prefixed nvSt nvEn neSt neEn E) (2 * nvSt + 2), neEn)).1.setNe
          neSt node neSt (nvSt, nvEn - 1) (2 * neSt, 2 * neEn - 1)) (2 * nvEn - 1) := rfl

theorem prefixed_edges_size (nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat)) :
    (prefixed nvSt nvEn neSt neEn E).edges.size = neEn - neSt := by
  unfold prefixed
  rw [foldl_prefixStep_edges, foldl_countStep_edges, inc_edges, inc_edges, empty_edges_size]

/-- Edges of the R layout: the cap at `0`, then the children in order. Unlike
`LayoutShape.run_edges` this needs nothing about the endpoints in `E`. -/
theorem run_edges (node nvSt nvEn neSt neEn : Nat) (E : List (Nat × Nat))
    (hne : neEn = neSt + E.length + 1) (k : Nat) (hk : k < neEn - neSt) :
    (run node nvSt nvEn neSt neEn E).edges[k]! =
      ⟨node, none, if k = 0 then (nvSt, nvEn - 1) else E[k - 1]?.getD (0, 0)⟩ := by
  rw [run_eq, inc_edges]
  have hl₂ : (inc nvSt (prefixed nvSt nvEn neSt neEn E) (2 * nvSt + 2)).edges.size =
      neEn - neSt := by
    rw [inc_edges, prefixed_edges_size]
  obtain ⟨hsz, hfill⟩ := foldl_fillStep_edges nvSt neSt node E.reverse
    (inc nvSt (prefixed nvSt nvEn neSt neEn E) (2 * nvSt + 2)) neEn
    (by rw [List.length_reverse]; omega) (by rw [hl₂]; omega)
  rw [setNe_edges_get _ _ _ _ _ _ _ (by rw [hsz, hl₂]; exact hk), Nat.sub_self]
  by_cases hk0 : k = 0
  · rw [ite_of_pos hk0, ite_of_pos hk0]
  · rw [ite_of_neg hk0, ite_of_neg hk0, hfill k (by rw [hl₂]; exact hk), List.length_reverse,
      ite_of_pos ⟨by omega, by omega⟩, List.getD_eq_getElem?_getD,
      List.getElem?_reverse (by omega)]
    congr 3; omega

end LayoutEdges

open LayoutEdges LayoutShape in
/-- Edge `k` of `layoutNode` for a node type: owned by `node`, no twin, endpoints `layoutNvs`. -/
theorem layoutNode_edges (ty : NodeType) (node nvSt nvEn neSt neEn : Nat) (ec : List (Nat × Nat))
    (hty : ty.isNode = true) (h1 : nvEn - nvSt = 1 → neEn - neSt ≤ 1)
    (hQI : ty = .Q ∨ ty = .I → neEn - neSt ≤ 1)
    (hR : ty = .R → nvSt + 2 ≤ nvEn ∧ neEn = neSt + ec.length + 1)
    (hO : ty = .O → nvEn - nvSt = 1)
    (k : Nat) (hk : k < neEn - neSt) :
    (layoutNode ty node nvSt nvEn neSt neEn ec).edges[k]! =
      ⟨node, none, layoutNvs ty nvSt nvEn ec k⟩ := by
  unfold layoutNvs
  by_cases hv1 : nvEn - nvSt = 1
  · have hk0 : k = 0 := by have := h1 hv1; omega
    subst hk0
    rw [ite_eq_left hv1, layoutNode_loop_eq _ _ _ _ _ _ _
      ⟨by rintro rfl; simp [NodeType.isNode] at hty, by rintro rfl; simp [NodeType.isNode] at hty⟩
      hv1]
    exact runLoop_edges_get _ _ _ _ _ (by have := h1 hv1; omega)
  · rw [ite_eq_right hv1]
    cases ty <;> simp [NodeType.isNode] at hty
    case O => exact absurd (hO rfl) hv1
    case R =>
      obtain ⟨hv, hne⟩ := hR rfl
      rw [LayoutR.layoutNode_R_eq _ _ _ _ _ _ hv, run_edges _ _ _ _ _ _ hne k hk]
      simp
    case Q =>
      have hk0 : k = 0 := by have := hQI (Or.inl rfl); omega
      subst hk0
      rw [layoutNode_QI_eq _ _ _ _ _ _ _ (Or.inl rfl) hv1]
      exact runQI_edges_get _ _ _ _ _ (by have := hQI (Or.inl rfl); omega)
    case I =>
      have hk0 : k = 0 := by have := hQI (Or.inr rfl); omega
      subst hk0
      rw [layoutNode_QI_eq _ _ _ _ _ _ _ (Or.inr rfl) hv1]
      exact runQI_edges_get _ _ _ _ _ (by have := hQI (Or.inr rfl); omega)
    case P =>
      rw [layoutNode_P_eq _ _ _ _ _ _ hv1]
      obtain ⟨-, -, -, ge, -⟩ := runP_spec node nvSt nvEn neSt neEn (by omega)
      exact ge k hk
    case S =>
      rw [layoutNode_S_eq _ _ _ _ _ _ hv1]
      obtain ⟨-, -, -, ge, -⟩ := runS_spec node nvSt nvEn neSt neEn (by omega)
      rw [ge k hk]
      split <;> rfl

/-! ### Range endpoints as `[·]!` -/

namespace SpqrTree
variable (t : SpqrTree)
theorem nvRange_fst (n : Nat) : (t.nvRange n).1 = t.nvBounds[n]! :=
  (Array.getElem!_eq_getD_getElem? _ _).symm
theorem nvRange_snd (n : Nat) : (t.nvRange n).2 = t.nvBounds[n + 1]! :=
  (Array.getElem!_eq_getD_getElem? _ _).symm
theorem neRange_fst (n : Nat) : (t.neRange n).1 = t.neBounds[n]! :=
  (Array.getElem!_eq_getD_getElem? _ _).symm
theorem neRange_snd (n : Nat) : (t.neRange n).2 = t.neBounds[n + 1]! :=
  (Array.getElem!_eq_getD_getElem? _ _).symm
theorem chRange_fst (n : Nat) : (t.chRange n).1 = t.chBounds[n]! :=
  (Array.getElem!_eq_getD_getElem? _ _).symm
theorem chRange_snd (n : Nat) : (t.chRange n).2 = t.chBounds[n + 1]! :=
  (Array.getElem!_eq_getD_getElem? _ _).symm
end SpqrTree

/-! ### The relabel interface as one hypothesis -/

/-- `relabel_node_spec` unpacked: well-formed items, the global record and every per-item record. -/
structure RelabelAll (g : Graph) (items : Items) (t : SpqrTree) (idx : ItemId → Nat) : Prop where
  wf : items.WF g
  gl : RelabelIdx g items t idx
  node : ∀ i, i < items.size → RelabelNode g items t idx i

namespace RelabelAll

variable {g : Graph} {items : Items} {t : SpqrTree} {idx : ItemId → Nat}
variable (H : RelabelAll g items t idx)
include H

theorem size : t.size = items.size := H.gl.size

theorem idx_lt {i : ItemId} (hi : i < items.size) : idx i < t.size := (H.node i hi).idx_lt

theorem idx_inj {i j : ItemId} (hi : i < items.size) (hj : j < items.size) (h : idx i = idx j) :
    i = j := H.gl.inj i j hi hj h

theorem idx_surj {n : Nat} (hn : n < t.size) : ∃ i, i < items.size ∧ idx i = n := by
  have hs := H.size
  have := Finset.surj_on_of_inj_on_of_card_le (s := Finset.range items.size)
    (t := Finset.range t.size) (fun i _ => idx i)
    (fun i hi => Finset.mem_range.2 (H.idx_lt (Finset.mem_range.1 hi)))
    (fun i j hi hj h => H.idx_inj (Finset.mem_range.1 hi) (Finset.mem_range.1 hj) h)
    (by simp [hs]) n (Finset.mem_range.2 hn)
  obtain ⟨i, hi, rfl⟩ := this
  exact ⟨i, Finset.mem_range.1 hi, rfl⟩

theorem type {i : ItemId} (hi : i < items.size) : t.type (idx i) = items.type i := (H.node i hi).type

theorem tree : items.Tree g := H.wf.tree

/-! #### Bounds -/


theorem nv_mono : ∀ n, n + 1 < t.nvBounds.size → t.nvBounds[n]! ≤ t.nvBounds[n + 1]! := by
  intro n hn
  rw [H.gl.sizes.nvBounds] at hn
  obtain ⟨i, hi, rfl⟩ := H.idx_surj (by omega : n < t.size)
  rw [← t.nvRange_fst, ← t.nvRange_snd, (H.node i hi).nv_range]; omega

theorem ne_mono : ∀ n, n + 1 < t.neBounds.size → t.neBounds[n]! ≤ t.neBounds[n + 1]! := by
  intro n hn
  rw [H.gl.sizes.neBounds] at hn
  obtain ⟨i, hi, rfl⟩ := H.idx_surj (by omega : n < t.size)
  rw [← t.neRange_fst, ← t.neRange_snd, (H.node i hi).ne_range]; omega

theorem ch_mono : ∀ n, n + 1 < t.chBounds.size → t.chBounds[n]! ≤ t.chBounds[n + 1]! := by
  intro n hn
  rw [H.gl.sizes.chBounds] at hn
  obtain ⟨i, hi, rfl⟩ := H.idx_surj (by omega : n < t.size)
  rw [← t.chRange_fst, ← t.chRange_snd, (H.node i hi).ch_range]; omega

/-- `nvEn a ≤ nvSt b` for `a < b`. -/
theorem nvEn_le_nvSt {a b : Nat} (hab : a < b) (hb : b < t.size) :
    (t.nvRange a).2 ≤ (t.nvRange b).1 := by
  rw [t.nvRange_snd, t.nvRange_fst]
  exact bounds_le_of_le _ H.nv_mono _ _ hab (by rw [H.gl.sizes.nvBounds]; omega)

theorem neEn_le_neSt {a b : Nat} (hab : a < b) (hb : b < t.size) :
    (t.neRange a).2 ≤ (t.neRange b).1 := by
  rw [t.neRange_snd, t.neRange_fst]
  exact bounds_le_of_le _ H.ne_mono _ _ hab (by rw [H.gl.sizes.neBounds]; omega)

theorem nvEn_le_size {a : Nat} (ha : a < t.size) : (t.nvRange a).2 ≤ t.nodeVerts.size := by
  rw [t.nvRange_snd, ← H.gl.nv_last]
  exact bounds_le_of_le _ H.nv_mono _ _ ha (by rw [H.gl.sizes.nvBounds]; omega)

theorem neEn_le_size {a : Nat} (ha : a < t.size) : (t.neRange a).2 ≤ t.nodeEdges.size := by
  rw [t.neRange_snd, ← H.gl.ne_last]
  exact bounds_le_of_le _ H.ne_mono _ _ ha (by rw [H.gl.sizes.neBounds]; omega)

theorem nvSt_le_nvEn {a : Nat} (ha : a < t.size) : (t.nvRange a).1 ≤ (t.nvRange a).2 := by
  rw [t.nvRange_snd, t.nvRange_fst]
  exact H.nv_mono a (by rw [H.gl.sizes.nvBounds]; omega)

theorem neSt_le_neEn {a : Nat} (ha : a < t.size) : (t.neRange a).1 ≤ (t.neRange a).2 := by
  rw [t.neRange_snd, t.neRange_fst]
  exact H.ne_mono a (by rw [H.gl.sizes.neBounds]; omega)

/-- Node-vert ranges of distinct nodes are disjoint. -/
theorem nv_disjoint {a b x : Nat} (ha : a < t.size) (hb : b < t.size) (hab : a ≠ b)
    (hxa : (t.nvRange a).1 ≤ x ∧ x < (t.nvRange a).2)
    (hxb : (t.nvRange b).1 ≤ x ∧ x < (t.nvRange b).2) : False := by
  rcases Nat.lt_or_gt_of_ne hab with h | h
  · have := H.nvEn_le_nvSt h hb; omega
  · have := H.nvEn_le_nvSt h ha; omega

theorem ne_disjoint {a b x : Nat} (ha : a < t.size) (hb : b < t.size) (hab : a ≠ b)
    (hxa : (t.neRange a).1 ≤ x ∧ x < (t.neRange a).2)
    (hxb : (t.neRange b).1 ≤ x ∧ x < (t.neRange b).2) : False := by
  rcases Nat.lt_or_gt_of_ne hab with h | h
  · have := H.neEn_le_neSt h hb; omega
  · have := H.neEn_le_neSt h ha; omega

/-- Every node-vert lies in the range of some item. -/
theorem nv_locate {nv : Nat} (h : nv < t.nodeVerts.size) :
    ∃ i, i < items.size ∧ (t.nvRange (idx i)).1 ≤ nv ∧ nv < (t.nvRange (idx i)).2 := by
  have hsz := H.gl.sizes.nvBounds
  obtain ⟨n, hn, h1, h2⟩ := bounds_locate t.nvBounds H.gl.nv_zero t.size
    (by omega) nv (by rw [H.gl.nv_last]; exact h)
  obtain ⟨i, hi, rfl⟩ := H.idx_surj hn
  exact ⟨i, hi, by rw [t.nvRange_fst]; exact h1, by rw [t.nvRange_snd]; exact h2⟩

theorem ne_locate {ne : Nat} (h : ne < t.nodeEdges.size) :
    ∃ i, i < items.size ∧ (t.neRange (idx i)).1 ≤ ne ∧ ne < (t.neRange (idx i)).2 := by
  have hsz := H.gl.sizes.neBounds
  obtain ⟨n, hn, h1, h2⟩ := bounds_locate t.neBounds H.gl.ne_zero t.size
    (by omega) ne (by rw [H.gl.ne_last]; exact h)
  obtain ⟨i, hi, rfl⟩ := H.idx_surj hn
  exact ⟨i, hi, by rw [t.neRange_fst]; exact h1, by rw [t.neRange_snd]; exact h2⟩

/-! #### Bijections -/

theorem bijections : t.Bijections := by
  have ht := H.tree
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · intro v hv
    rw [H.gl.nv] at hv
    have hi : vertItem v < items.size := by have := ht.size; simp only [vertItem]; iomega
    refine ⟨idx (vertItem v), H.gl.vert_index v hv, ?_, ?_⟩
    · rw [H.type hi, ht.vert v hv]
    · rw [(H.node _ hi).orig, Items.origOf, ht.vert v hv]; simp [vertItem]
  · intro n hn hV
    obtain ⟨i, hi, rfl⟩ := H.idx_surj hn
    rw [H.type hi] at hV
    have hr := (ht.type_V_iff hi).1 hV
    refine ⟨i - 1, ?_, (H.node i hi).vert_index hV⟩
    rw [(H.node i hi).orig, Items.origOf, hV]
  · intro e he
    rw [H.gl.ne] at he
    have hi : edgeItem g e < items.size := by have := ht.size; simp only [edgeItem]; iomega
    refine ⟨idx (edgeItem g e), H.gl.edge_index e he, ?_, ?_⟩
    · rw [H.type hi, ht.edge e he]
    · rw [(H.node _ hi).orig, Items.origOf, ht.edge e he]; simp only [edgeItem]; congr 1; omega
  · intro n hn hQ
    obtain ⟨i, hi, rfl⟩ := H.idx_surj hn
    rw [H.type hi] at hQ
    refine ⟨i - 1 - g.nv, ?_, ((H.node i hi).edge_index hQ).1⟩
    rw [(H.node i hi).orig, Items.origOf, hQ]
  · intro n hn hV hQ
    obtain ⟨i, hi, rfl⟩ := H.idx_surj hn
    rw [H.type hi] at hV hQ
    rw [(H.node i hi).orig, Items.origOf]
    cases h : items.type i <;> simp_all

/-! #### Ownership: node-verts -/

theorem type_child {p c : ItemId} (hc : c ∈ items.ch p) : t.type (idx c) = items.type c :=
  H.type (H.tree.child_lt hc)

theorem children_perm {i : ItemId} (hi : i < items.size) :
    (t.children (idx i)).Perm ((items.ch i).map idx) := by
  obtain ⟨pos, hl⟩ := (H.node i hi).layout
  rw [hl.children]
  exact (Items.ordered_perm i _ pos).map idx

theorem hasCap_eq {i : ItemId} (hi : i < items.size) : t.hasCap (idx i) = items.hasCap i := by
  unfold SpqrTree.hasCap Items.hasCap
  rw [H.type hi]
  by_cases hQ : items.type i = .Q
  · have : ((t.children (idx i)).any fun c => t.type c != .V) = !(items.ch i).isEmpty := by
      rcases H.wf.q_cases hi hQ with h0 | ⟨c, hc, h1⟩
      · have h2 : t.children (idx i) = [] := by
          have := H.children_perm hi; rw [h0] at this; exact this.eq_nil
        simp [h2, h0]
      · have hcm : c ∈ items.ch i := by rcases h1 with h1 | ⟨v, h1⟩ <;> simp [h1]
        have hne : (items.ch i).isEmpty = false := by rcases h1 with h1 | ⟨v, h1⟩ <;> simp [h1]
        rw [hne, Bool.not_false, List.any_eq_true]
        refine ⟨idx c, (H.children_perm hi).mem_iff.2 (List.mem_map_of_mem hcm), ?_⟩
        rw [H.type_child hcm, bne_iff_ne]
        simp at hc; exact hc.2.1
    rw [this]
  · have h := beq_eq_false_iff_ne.2 hQ
    simp only [h, Bool.false_and]

theorem nodeVert_at {i : ItemId} (hi : i < items.size) {k : Nat}
    (hk : k < (items.nvList g i).length) :
    t.nodeVerts[(t.nvRange (idx i)).1 + k]! = ⟨idx i, idx (vertItem (items.nvList g i)[k])⟩ := by
  rw [(H.node i hi).node_verts k hk,
    H.gl.vert_index _ (Items.nvList_mem_lt H.tree H.wf.endpoints (List.getElem_mem hk))]
  rfl

theorem nodeVertsOf_eq {i : ItemId} (hi : i < items.size) :
    t.nodeVertsOf (idx i) = (items.nvList g i).map fun v => ⟨idx i, idx (vertItem v)⟩ := by
  unfold SpqrTree.nodeVertsOf
  rw [(H.node i hi).nv_range, Nat.add_sub_cancel_left]
  apply List.ext_getElem (by simp)
  intro k h1 h2
  simp only [List.getElem_map, List.getElem_range]
  rw [← Array.getElem!_eq_getD_getElem?, H.nodeVert_at hi (by simpa using h2)]

theorem vertItem_lt {v : Nat} (hv : v < g.nv) : vertItem v < items.size := by
  have := H.tree.size; simp only [vertItem]; iomega

/-! #### Ownership: node-edges -/

theorem nodeEdge_at (hx : items.OwnExtra g) {i : ItemId} (hi : i < items.size) {k : Nat}
    (hk : k < items.nEdges g i) :
    ∃ pos, RelabelLayout g items t idx i pos ∧
      (t.nodeEdges[(t.neRange (idx i)).1 + k]!).node = idx i ∧
      (t.nodeEdges[(t.neRange (idx i)).1 + k]!).nvs =
        layoutNvs (items.type i) (t.nvRange (idx i)).1 (t.nvRange (idx i)).2
          (items.edgeChildren g pos (items.ordered g i (t.nvRange (idx i)).1 pos)) k := by
  obtain ⟨pos, hl⟩ := (H.node i hi).layout
  have hn : (items.type i).isNode = true := by
    by_contra h
    rw [Items.nEdges_eq_zero (Bool.eq_false_iff.2 h)] at hk; omega
  have hy := H.wf.layout_hyps hx hi hn
  have hnv := (H.node i hi).nv_range
  have hne := (H.node i hi).ne_range
  have hec := Items.edgeChildren_length (g := g) (items := items) i (t.nvRange (idx i)).1 pos
  have hL := layoutNode_edges (items.type i) (idx i) (t.nvRange (idx i)).1 (t.nvRange (idx i)).2
    (t.neRange (idx i)).1 (t.neRange (idx i)).2
    (items.edgeChildren g pos (items.ordered g i (t.nvRange (idx i)).1 pos)) hn
    (fun h => by have := hy.1 (by omega); omega)
    (fun h => by have := hy.2.1 h; omega)
    (fun hR => by
      have := hy.2.2.2.1 hR
      have hc := Items.hasCap_of_ne_Q hn (by rw [hR]; decide)
      refine ⟨by omega, ?_⟩
      rw [hec, hne, Items.nEdges_eq hn]; simp only [Items.capCount, hc, ↓reduceIte]; omega)
    (fun hO => by have := hy.2.2.2.2.2 hO; omega)
    k (by omega)
  refine ⟨pos, hl, ?_, ?_⟩
  · rw [hl.edge_node k hk]; unfold nodeLayout; rw [hL]
  · rw [hl.edge_nvs k hk]; unfold nodeLayout; rw [hL]

omit H in
theorem isNode_of_nEdges_pos {i : ItemId} (h : 0 < items.nEdges g i) :
    (items.type i).isNode = true := by
  by_contra h'
  rw [Items.nEdges_eq_zero (Bool.eq_false_iff.2 h')] at h; omega

/-- The `k`-th non-V child of `i` in output order has both endpoints among `i`'s node-verts
(R nodes), and they are the pair `edgeChildren` lists at `k`. -/
theorem edgeChild_at (hx : items.OwnExtra g) {i : ItemId} (hi : i < items.size)
    (hR : items.type i = .R) {pos : Nat → Nat} {nvSt : Nat} {k : Nat}
    (hk : k < (items.edgeChildren g pos (items.ordered g i nvSt pos)).length) :
    ∃ u v, u ∈ items.nvList g i ∧ v ∈ items.nvList g i ∧ (u, v) ∈ items.virtualEdges i ∧
      (items.edgeChildren g pos (items.ordered g i nvSt pos))[k]? = some (pos u, pos v) := by
  have hk' : k < ((items.ordered g i nvSt pos).filter (· ≥ 1 + g.nv)).length := by
    simpa [Items.edgeChildren] using hk
  unfold Items.edgeChildren
  rw [List.getElem?_map, List.getElem?_eq_getElem hk', Option.map_some]
  have hmem := List.mem_filter.1 (List.getElem_mem hk')
  have hc : ((items.ordered g i nvSt pos).filter (· ≥ 1 + g.nv))[k]'hk' ∈ items.ch i :=
    Items.mem_ordered.1 hmem.1
  have hge : 1 + g.nv ≤ ((items.ordered g i nvSt pos).filter (· ≥ 1 + g.nv))[k]'hk' := by
    simpa using hmem.2
  have hV : items.type (((items.ordered g i nvSt pos).filter (· ≥ 1 + g.nv))[k]'hk') ≠ .V := by
    intro h; have := ((H.tree.type_V_iff (H.tree.child_lt hc)).1 h).2; iomega
  have hve : ((items.vs _).1.getD 0, (items.vs _).2.getD 0) ∈ items.virtualEdges i :=
    List.mem_map_of_mem (List.mem_filter.2 ⟨hc, by simpa using hV⟩)
  obtain ⟨h1, h2⟩ := hx.r_edges_in_nv i hi hR _ hve
  exact ⟨_, _, h1, h2, hve, rfl⟩

theorem layoutNvs_bounds (hx : items.OwnExtra g) (hor : items.ROriented g) {i : ItemId}
    (hi : i < items.size) {pos : Nat → Nat} (hl : RelabelLayout g items t idx i pos) {k : Nat}
    (hk : k < items.nEdges g i) :
    (t.nvRange (idx i)).1 ≤ (layoutNvs (items.type i) (t.nvRange (idx i)).1 (t.nvRange (idx i)).2
        (items.edgeChildren g pos (items.ordered g i (t.nvRange (idx i)).1 pos)) k).1 ∧
    (layoutNvs (items.type i) (t.nvRange (idx i)).1 (t.nvRange (idx i)).2
        (items.edgeChildren g pos (items.ordered g i (t.nvRange (idx i)).1 pos)) k).1 ≤
      (layoutNvs (items.type i) (t.nvRange (idx i)).1 (t.nvRange (idx i)).2
        (items.edgeChildren g pos (items.ordered g i (t.nvRange (idx i)).1 pos)) k).2 ∧
    (layoutNvs (items.type i) (t.nvRange (idx i)).1 (t.nvRange (idx i)).2
        (items.edgeChildren g pos (items.ordered g i (t.nvRange (idx i)).1 pos)) k).2 <
      (t.nvRange (idx i)).2 := by
  have hn := isNode_of_nEdges_pos (g := g) (items := items) (i := i) (by omega)
  have hy := H.wf.layout_hyps hx hi hn
  have hnv := (H.node i hi).nv_range
  have h1 := H.wf.one_le_nvList hi hn
  have hec := Items.edgeChildren_length (g := g) (items := items) i (t.nvRange (idx i)).1 pos
  unfold layoutNvs
  split_ifs with hA hB hC hk0 hk0
  · exact ⟨le_rfl, le_rfl, by omega⟩
  · exact ⟨le_rfl, by omega, by omega⟩
  · exact ⟨le_rfl, by omega, by omega⟩
  · have := hy.2.2.1 hC; refine ⟨by omega, by omega, by omega⟩
  · exact ⟨le_rfl, by omega, by omega⟩
  · have hR : items.type i = .R := by
      cases hty : items.type i <;> simp [hty, NodeType.isNode] at hn hA hB hC ⊢
      have := hy.2.2.2.2.2 hty; omega
    have hcap := Items.hasCap_of_ne_Q hn (by rw [hR]; decide)
    have hlen : k - 1 < (items.edgeChildren g pos (items.ordered g i (t.nvRange (idx i)).1 pos)).length := by
      rw [hec]; have := Items.nEdges_eq (g := g) hn
      simp only [Items.capCount, hcap, ↓reduceIte] at this; omega
    obtain ⟨u, v, hu, hv, huv, heq⟩ := H.edgeChild_at hx hi hR hlen
    rw [heq, Option.getD_some]
    have hpos := hl.pos_ok hR
    have hu1 := (hpos u hu).1
    have hv1 := (hpos v hv).1
    have hu2 := Items.PosOK.pos_lt hpos hu
    have hv2 := Items.PosOK.pos_lt hpos hv
    obtain ⟨hnd, hord⟩ := hor i hi hR
    have h3 := hord _ huv
    rw [← Items.PosOK.pos_sub hnd hpos hu, ← Items.PosOK.pos_sub hnd hpos hv] at h3
    exact ⟨by omega, by omega, by omega⟩

theorem layoutNvs_cap (hx : items.OwnExtra g) {i : ItemId} (hi : i < items.size)
    (hcap : items.hasCap i = true) (ec : List (Nat × Nat)) :
    layoutNvs (items.type i) (t.nvRange (idx i)).1 (t.nvRange (idx i)).2 ec 0 =
      ((t.nvRange (idx i)).1, (t.nvRange (idx i)).2 - 1) := by
  have hn : (items.type i).isNode = true := by
    unfold Items.hasCap at hcap; exact (Bool.and_eq_true_iff.1 hcap).1
  have hy := H.wf.layout_hyps hx hi hn
  have hnv := (H.node i hi).nv_range
  have h1 := H.wf.one_le_nvList hi hn
  unfold layoutNvs
  simp only [↓reduceIte]
  split_ifs with hA hB hC
  · rw [Prod.mk.injEq]; omega
  · have := hy.2.2.2.2.1 hcap hB; rw [Prod.mk.injEq]; omega
  · rfl
  · rfl

/-! #### Ownership -/

theorem ownership (hx : items.OwnExtra g) (hor : items.ROriented g) : t.Ownership := by
  have ht := H.tree
  refine ⟨H.nv_mono, H.ne_mono, H.gl.nv_last, H.gl.ne_last, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · -- nv_node
    intro n nv hn h1 h2
    obtain ⟨i, hi, rfl⟩ := H.idx_surj hn
    have hk : nv - (t.nvRange (idx i)).1 < (items.nvList g i).length := by
      have := (H.node i hi).nv_range; omega
    have hsz : nv < t.nodeVerts.size := lt_of_lt_of_le h2 (H.nvEn_le_size hn)
    unfold SpqrTree.nodeOfNv
    rw [Array.getElem!_of_lt _ _ hsz]
    have := H.nodeVert_at hi hk
    rw [Nat.add_sub_cancel' h1] at this
    rw [this]; rfl
  · -- nv_vert
    intro nv hnv
    obtain ⟨i, hi, h1, h2⟩ := H.nv_locate hnv
    have hk : nv - (t.nvRange (idx i)).1 < (items.nvList g i).length := by
      have := (H.node i hi).nv_range; omega
    have := H.nodeVert_at hi hk
    rw [Nat.add_sub_cancel' h1] at this
    rw [this]
    show t.type (idx (vertItem _)) = .V
    have hv := Items.nvList_mem_lt ht H.wf.endpoints (List.getElem_mem hk)
    rw [H.type (H.vertItem_lt hv)]
    exact ht.vert _ hv
  · -- nv_distinct
    intro n hn
    obtain ⟨i, hi, rfl⟩ := H.idx_surj hn
    rw [H.nodeVertsOf_eq hi, List.map_map]
    have hf : ((fun x : NodeVert => x.vert) ∘ fun v => (⟨idx i, idx (vertItem v)⟩ : NodeVert)) =
        fun v => idx (vertItem v) := rfl
    rw [hf]
    refine (hx.nv_nodup i hi).map_on ?_
    intro x hx' y hy hxy
    have hx1 := Items.nvList_mem_lt ht H.wf.endpoints hx'
    have hy1 := Items.nvList_mem_lt ht H.wf.endpoints hy
    have := H.idx_inj (H.vertItem_lt hx1) (H.vertItem_lt hy1) hxy
    simp only [vertItem] at this; iomega
  · -- ne_node
    intro n ne hn h1 h2
    obtain ⟨i, hi, rfl⟩ := H.idx_surj hn
    have hk : ne - (t.neRange (idx i)).1 < items.nEdges g i := by
      have := (H.node i hi).ne_range; omega
    have hsz : ne < t.nodeEdges.size := lt_of_lt_of_le h2 (H.neEn_le_size hn)
    obtain ⟨pos, hl, hnode, -⟩ := H.nodeEdge_at hx hi hk
    rw [Nat.add_sub_cancel' h1] at hnode
    unfold SpqrTree.nodeOfNe
    rw [Array.getElem!_of_lt _ _ hsz, Option.map_some, hnode]
  · -- ne_nvs
    intro ne hne n hn
    obtain ⟨i, hi, h1, h2⟩ := H.ne_locate hne
    have hk : ne - (t.neRange (idx i)).1 < items.nEdges g i := by
      have := (H.node i hi).ne_range; omega
    obtain ⟨pos, hl, hnode, hnvs⟩ := H.nodeEdge_at hx hi hk
    rw [Nat.add_sub_cancel' h1] at hnode hnvs
    have hn' : n = idx i := by
      unfold SpqrTree.nodeOfNe at hn
      rw [Array.getElem!_of_lt _ _ hne, Option.map_some, hnode] at hn
      exact (Option.some.inj hn).symm
    subst hn'
    rw [hnvs]
    exact H.layoutNvs_bounds hx hor hi hl hk
  · -- vert_par_nv
    intro n p hn hV hp
    obtain ⟨j, hj, rfl⟩ := H.idx_surj hn
    rw [H.type hj] at hV
    have hj0 : j ≠ 0 := fun h => by rw [(ht.type_F_iff hj).2 h] at hV; cases hV
    obtain ⟨q, hq⟩ := ht.parent_exists hj hj0
    have hqs := ht.parent_lt hq
    have hpar := (H.node q hqs).child_par j hq
    rw [hp] at hpar
    obtain rfl := Option.some.inj hpar
    obtain ⟨pos, hl⟩ := (H.node q hqs).layout
    have hjV : j < 1 + g.nv := ht.child_lt_of_V hq hV
    have hfil := Items.ordered_filter_V ht (hx.nv_nodup q hqs) hl.pos_ok
    have hmem : j ∈ (items.ordered g q (t.nvRange (idx q)).1 pos).filter (· < 1 + g.nv) := by
      rw [hfil]; exact List.mem_filter.2 ⟨hq, by simpa using hjV⟩
    obtain ⟨k, hk, hjk⟩ := List.mem_iff_getElem.1 hmem
    have hvp := hl.vert_par_nv k hk
    rw [hjk] at hvp
    have hk' : k < ((items.ch q).filter (· < 1 + g.nv)).length := hfil ▸ hk
    have hlen := Items.nvList_length (g := g) (items := items) q
    have hnv := (H.node q hqs).nv_range
    refine ⟨_, hvp, by omega, by omega, ?_⟩
    rw [Nat.add_assoc, H.nodeVert_at hqs (by omega)]
    show idx (vertItem _) = idx j
    congr 1
    have hmid := Items.nvList_getElem?_mid (g := g) (items := items) q hk'
    rw [List.getElem?_eq_getElem (by omega), Option.some.injEq] at hmid
    rw [hmid, ← List.getElem_of_eq hfil hk, hjk]
    have := ht.child_pos hq
    simp only [vertItem]; iomega
  · -- vert_par_nv_none
    intro n hn hV
    obtain ⟨j, hj, rfl⟩ := H.idx_surj hn
    rw [H.type hj] at hV
    by_cases hj0 : j = 0
    · subst hj0
      show t.vertParNv[idx rootItem]! = none
      rw [H.gl.root]; exact H.gl.root_par_nv
    · obtain ⟨q, hq⟩ := ht.parent_exists hj hj0
      exact (H.node q (ht.parent_lt hq)).child_par_nv_none j hq (ht.child_ge_of_ne_V hq hV)
  · -- nv_layout
    intro n hn
    obtain ⟨i, hi, rfl⟩ := H.idx_surj hn
    obtain ⟨pos, hl⟩ := (H.node i hi).layout
    refine ⟨(items.vs i).1.toList.map (fun v => ⟨idx i, idx (vertItem v)⟩),
      (items.vs i).2.toList.map (fun v => ⟨idx i, idx (vertItem v)⟩), ?_,
      by rw [List.length_map]; exact Option.toList_length_le _,
      by rw [List.length_map]; exact Option.toList_length_le _, ?_⟩
    · rw [H.nodeVertsOf_eq hi]
      unfold Items.nvList
      rw [List.map_append, List.map_append, List.map_map]
      congr 2
      rw [hl.children, List.filter_map, List.map_map]
      have hfil : (items.ordered g i (t.nvRange (idx i)).1 pos).filter
          ((fun c => decide (t.type c = .V)) ∘ idx) =
          (items.ordered g i (t.nvRange (idx i)).1 pos).filter (· < 1 + g.nv) := by
        apply List.filter_congr
        intro c hc
        have hc' := Items.mem_ordered.1 hc
        simp only [Function.comp, H.type_child hc', decide_eq_decide]
        rw [ht.type_V_iff (ht.child_lt hc')]
        have := ht.child_pos hc'
        constructor
        · exact fun h => h.2
        · intro h; exact ⟨by iomega, h⟩
      rw [hfil, Items.ordered_filter_V ht (hx.nv_nodup i hi) hl.pos_ok]
      apply List.map_congr_left
      intro c hc
      have := ht.child_pos (List.mem_filter.1 hc).1
      simp only [Function.comp]
      congr 2
      simp only [vertItem]; iomega
    · intro hcap nvs hnvs
      rw [H.hasCap_eq hi] at hcap
      have hn1 : 0 < items.nEdges g i := by
        have hn := (Bool.and_eq_true_iff.1 (by unfold Items.hasCap at hcap; exact hcap)).1
        rw [Items.nEdges_eq hn]; simp [Items.capCount, hcap]
      obtain ⟨pos', hl', -, hnvs'⟩ := H.nodeEdge_at hx hi hn1
      rw [Nat.add_zero] at hnvs'
      have hsz : (t.neRange (idx i)).1 < t.nodeEdges.size := by
        have := (H.node i hi).ne_range; have := H.neEn_le_size (H.idx_lt hi); omega
      unfold SpqrTree.nvsOf at hnvs
      rw [Array.getElem!_of_lt _ _ hsz, Option.map_some, Option.some.injEq] at hnvs
      rw [← hnvs, hnvs', H.layoutNvs_cap hx hi hcap]
      exact ⟨rfl, rfl⟩

/-! #### Twins -/

omit H in
theorem filter_nonV_length {i : ItemId} (hn : (items.type i).isNode = true) (nvSt : Nat)
    (pos : Nat → Nat) :
    ((items.ordered g i nvSt pos).filter (· ≥ 1 + g.nv)).length =
      items.nEdges g i - items.capCount i := by
  rw [← List.countP_eq_length_filter, (Items.ordered_perm i nvSt pos).countP_eq,
    Items.nEdges_eq hn, Nat.add_sub_cancel]

/-- A non-V child of a node has a cap (a Q child of a node is a leaf, by `OwnExtra`). -/
theorem child_hasCap (hx : items.OwnExtra g) {p c : ItemId} (hp : (items.type p).isNode = true)
    (hc : c ∈ items.ch p) (hge : 1 + g.nv ≤ c) : items.hasCap c = true := by
  have hn := H.tree.isNode_of_ge (H.tree.child_lt hc) hge
  by_cases hQ : items.type c = .Q
  · rw [Items.hasCap_Q hQ, hx.q_leaf_of_node p c hc hp hQ]; rfl
  · exact Items.hasCap_of_ne_Q hn hQ

theorem child_nEdges_pos (hx : items.OwnExtra g) {p c : ItemId} (hp : (items.type p).isNode = true)
    (hc : c ∈ items.ch p) (hge : 1 + g.nv ≤ c) : 0 < items.nEdges g c := by
  have hcap := H.child_hasCap hx hp hc hge
  rw [Items.nEdges_eq (H.tree.isNode_of_ge (H.tree.child_lt hc) hge)]
  simp [Items.capCount, hcap]

omit H in
theorem twin_none_of_ge {ne : Nat} (h : t.nodeEdges.size ≤ ne) : t.twin ne = none := by
  unfold SpqrTree.twin; rw [Array.getElem?_eq_none h]; rfl

omit H in
theorem twin_of_lt {ne : Nat} (h : ne < t.nodeEdges.size) :
    t.twin ne = (t.nodeEdges[ne]?.getD default).twin := by
  unfold SpqrTree.twin; rw [Array.getElem?_eq_getElem h]; rfl

theorem ne_decomp {ne : Nat} (h : ne < t.nodeEdges.size) :
    ∃ j k, j < items.size ∧ k < items.nEdges g j ∧ ne = (t.neRange (idx j)).1 + k := by
  obtain ⟨j, hj, h1, h2⟩ := H.ne_locate h
  have := (H.node j hj).ne_range
  exact ⟨j, ne - (t.neRange (idx j)).1, hj, by omega, by omega⟩

theorem ne_lt_size {j : ItemId} (hj : j < items.size) {k : Nat} (hk : k < items.nEdges g j) :
    (t.neRange (idx j)).1 + k < t.nodeEdges.size := by
  have := (H.node j hj).ne_range; have := H.neEn_le_size (H.idx_lt hj); omega

theorem nodeOfNe_at (hx : items.OwnExtra g) {j : ItemId} (hj : j < items.size) {k : Nat}
    (hk : k < items.nEdges g j) : t.nodeOfNe ((t.neRange (idx j)).1 + k) = some (idx j) := by
  obtain ⟨pos, hl, hnode, -⟩ := H.nodeEdge_at hx hj hk
  unfold SpqrTree.nodeOfNe
  rw [Array.getElem!_of_lt _ _ (H.ne_lt_size hj hk), Option.map_some, hnode]

/-- The cap of a node `j` is twinned with edge `capCount p + k'` of its parent `p` when `p` is a
node (where `j` is the `k'`-th non-V child of `p`), and has no twin otherwise. -/
theorem twin_cap {j : ItemId} (hj : j < items.size) (hcap : items.hasCap j = true) :
    (∃ p, j ∈ items.ch p ∧ (items.type p).isNode = true ∧
      ∃ pos, RelabelLayout g items t idx p pos ∧
      ∃ k', ∃ hk' : k' < ((items.ordered g p (t.nvRange (idx p)).1 pos).filter (· ≥ 1 + g.nv)).length,
        ((items.ordered g p (t.nvRange (idx p)).1 pos).filter (· ≥ 1 + g.nv))[k'] = j ∧
        t.twin (t.neRange (idx j)).1 = some ((t.neRange (idx p)).1 + items.capCount p + k') ∧
        t.twin ((t.neRange (idx p)).1 + items.capCount p + k') = some (t.neRange (idx j)).1) ∨
    (∃ p, j ∈ items.ch p ∧ (items.type p).isNode = false ∧ t.twin (t.neRange (idx j)).1 = none) := by
  have ht := H.tree
  have hn : (items.type j).isNode = true := by
    unfold Items.hasCap at hcap; exact (Bool.and_eq_true_iff.1 hcap).1
  have hge : 1 + g.nv ≤ j := (ht.isNode_iff hj).1 hn
  obtain ⟨p, hp⟩ := ht.parent_exists hj (by iomega)
  have hps := ht.parent_lt hp
  by_cases hpn : (items.type p).isNode = true
  · left
    obtain ⟨pos, hl⟩ := (H.node p hps).layout
    have hmem : j ∈ (items.ordered g p (t.nvRange (idx p)).1 pos).filter (· ≥ 1 + g.nv) :=
      List.mem_filter.2 ⟨Items.mem_ordered.2 hp, by simpa using hge⟩
    obtain ⟨k', hk', hjk⟩ := List.mem_iff_getElem.1 hmem
    obtain ⟨h1, h2⟩ := hl.twin hpn k' hk'
    rw [hjk] at h1 h2
    exact ⟨p, hp, hpn, pos, hl, k', hk', hjk, h2, h1⟩
  · right
    exact ⟨p, hp, Bool.eq_false_iff.2 hpn,
      (H.node p hps).child_cap_twin_none hpn j hp hcap⟩

theorem twin_invol : ∀ ne ne', t.twin ne = some ne' → t.twin ne' = some ne := by
  intro ne ne' h
  by_cases hne : ne < t.nodeEdges.size
  swap
  · rw [twin_none_of_ge (Nat.le_of_not_lt hne)] at h; cases h
  obtain ⟨j, k, hj, hk, rfl⟩ := H.ne_decomp hne
  have hn := isNode_of_nEdges_pos (g := g) (items := items) (i := j) (by omega)
  by_cases hc : items.capCount j ≤ k
  · obtain ⟨pos, hl⟩ := (H.node j hj).layout
    have hk' : k - items.capCount j <
        ((items.ordered g j (t.nvRange (idx j)).1 pos).filter (· ≥ 1 + g.nv)).length := by
      rw [filter_nonV_length hn]; omega
    obtain ⟨h1, h2⟩ := hl.twin hn (k - items.capCount j) hk'
    rw [show (t.neRange (idx j)).1 + items.capCount j + (k - items.capCount j) =
      (t.neRange (idx j)).1 + k by omega] at h1 h2
    rw [h1] at h
    obtain rfl := Option.some.inj h
    exact h2
  · have hk0 : k = 0 := by have := Items.capCount_le (items := items) j; omega
    have hcap : items.hasCap j = true := by
      unfold Items.capCount at hc; by_contra h'; simp [h'] at hc
    subst hk0
    rw [Nat.add_zero] at h ⊢
    rcases H.twin_cap hj hcap with ⟨p, -, -, pos, -, k', -, -, h1, h2⟩ | ⟨p, -, -, h1⟩
    · rw [h1] at h; obtain rfl := Option.some.inj h; exact h2
    · rw [h1] at h; cases h

theorem twin_ne (hx : items.OwnExtra g) : ∀ ne ne', t.twin ne = some ne' → ne ≠ ne' := by
  intro ne ne' h
  by_cases hne : ne < t.nodeEdges.size
  swap
  · rw [twin_none_of_ge (Nat.le_of_not_lt hne)] at h; cases h
  obtain ⟨j, k, hj, hk, rfl⟩ := H.ne_decomp hne
  have hn := isNode_of_nEdges_pos (g := g) (items := items) (i := j) (by omega)
  have hnode := H.nodeOfNe_at hx hj hk
  by_cases hc : items.capCount j ≤ k
  · obtain ⟨pos, hl⟩ := (H.node j hj).layout
    have hk' : k - items.capCount j <
        ((items.ordered g j (t.nvRange (idx j)).1 pos).filter (· ≥ 1 + g.nv)).length := by
      rw [filter_nonV_length hn]; omega
    obtain ⟨h1, -⟩ := hl.twin hn (k - items.capCount j) hk'
    rw [show (t.neRange (idx j)).1 + items.capCount j + (k - items.capCount j) =
      (t.neRange (idx j)).1 + k by omega, h] at h1
    obtain rfl := Option.some.inj h1
    have hmem := List.mem_filter.1 (List.getElem_mem hk')
    have hcm := Items.mem_ordered.1 hmem.1
    have hge : 1 + g.nv ≤ ((items.ordered g j (t.nvRange (idx j)).1 pos).filter (· ≥ 1 + g.nv))[k - items.capCount j] := by
      simpa using hmem.2
    have hpos := H.child_nEdges_pos hx hn hcm hge
    have hnode' := H.nodeOfNe_at hx (H.tree.child_lt hcm) hpos
    rw [Nat.add_zero] at hnode'
    intro heq
    rw [heq, hnode'] at hnode
    exact H.tree.child_ne hcm (H.idx_inj (H.tree.child_lt hcm) hj (Option.some.inj hnode))
  · have hk0 : k = 0 := by have := Items.capCount_le (items := items) j; omega
    have hcap : items.hasCap j = true := by
      unfold Items.capCount at hc; by_contra h'; simp [h'] at hc
    subst hk0
    rw [Nat.add_zero] at h hnode ⊢
    rcases H.twin_cap hj hcap with ⟨p, hp, hpn, pos, hl, k', hk', -, h1, -⟩ | ⟨p, -, -, h1⟩
    · rw [h1] at h; obtain rfl := Option.some.inj h
      have hps := H.tree.parent_lt hp
      have hkp : items.capCount p + k' < items.nEdges g p := by
        rw [filter_nonV_length hpn] at hk'; omega
      have hnode' := H.nodeOfNe_at hx hps hkp
      rw [← Nat.add_assoc] at hnode'
      intro heq
      rw [heq, hnode'] at hnode
      exact H.tree.child_ne hp (H.idx_inj hj hps (Option.some.inj hnode).symm)
    · rw [h1] at h; cases h

theorem twins (hx : items.OwnExtra g) : t.Twins := by
  have ht := H.tree
  refine ⟨H.twin_invol, H.twin_ne hx, ?_, ?_, ?_⟩
  · -- twin_parent
    intro n ne hn hcap
    obtain ⟨j, hj, rfl⟩ := H.idx_surj hn
    unfold SpqrTree.capNe at hcap
    rw [H.hasCap_eq hj] at hcap
    split_ifs at hcap with hc
    obtain rfl := Option.some.inj hcap
    rcases H.twin_cap hj hc with ⟨p, hp, hpn, pos, hl, k', hk', -, h1, -⟩ | ⟨p, hp, hpn, h1⟩
    · left
      have hps := ht.parent_lt hp
      refine ⟨idx p, (H.node p hps).child_par j hp, by rw [H.type hps]; exact hpn, _, h1, ?_⟩
      have hkp : items.capCount p + k' < items.nEdges g p := by
        rw [filter_nonV_length hpn] at hk'; omega
      have := H.nodeOfNe_at hx hps hkp
      rwa [← Nat.add_assoc] at this
    · right
      have hps := ht.parent_lt hp
      exact ⟨h1, idx p, (H.node p hps).child_par j hp, by rw [H.type hps]; simp [hpn]⟩
  · -- noncap_children
    intro n hn hnode
    obtain ⟨j, hj, rfl⟩ := H.idx_surj hn
    rw [H.type hj] at hnode
    obtain ⟨pos, hl⟩ := (H.node j hj).layout
    have hcc : (if t.hasCap (idx j) = true then 1 else 0) = items.capCount j := by
      rw [H.hasCap_eq hj]; rfl
    rw [hcc, hl.children, List.filter_map, List.map_map]
    have hfil : (items.ordered g j (t.nvRange (idx j)).1 pos).filter
        ((fun c => decide (t.type c ≠ .V)) ∘ idx) =
        (items.ordered g j (t.nvRange (idx j)).1 pos).filter (· ≥ 1 + g.nv) := by
      apply List.filter_congr
      intro c hc
      have hc' := Items.mem_ordered.1 hc
      simp only [Function.comp, H.type_child hc', decide_eq_decide]
      constructor
      · exact fun h => ht.child_ge_of_ne_V hc' h
      · intro h hV; have := ht.child_lt_of_V hc' hV; iomega
    rw [hfil]
    have hne := (H.node j hj).ne_range
    have hlen := filter_nonV_length (g := g) (items := items) hnode (t.nvRange (idx j)).1 pos
    apply List.ext_getElem
    · simp [SpqrTree.nodeEdgesOf, hne, hlen]
    intro k' h1 h2
    have hk' : k' < ((items.ordered g j (t.nvRange (idx j)).1 pos).filter (· ≥ 1 + g.nv)).length := by
      simpa using h2
    have hkn : items.capCount j + k' < items.nEdges g j := by rw [hlen] at hk'; omega
    obtain ⟨htw, -⟩ := hl.twin hnode k' hk'
    simp only [List.getElem_map, List.getElem_drop, SpqrTree.nodeEdgesOf, List.getElem_range,
      Function.comp]
    rw [← twin_of_lt (H.ne_lt_size hj hkn), ← Nat.add_assoc, htw]
    have hmem := List.mem_filter.1 (List.getElem_mem hk')
    have hcm := Items.mem_ordered.1 hmem.1
    have hge : 1 + g.nv ≤ ((items.ordered g j (t.nvRange (idx j)).1 pos).filter (· ≥ 1 + g.nv))[k'] := by
      simpa using hmem.2
    unfold SpqrTree.capNe
    rw [H.hasCap_eq (ht.child_lt hcm), H.child_hasCap hx hnode hcm hge]
    rfl
  · -- cap_none
    intro n hn hnode
    obtain ⟨j, hj, rfl⟩ := H.idx_surj hn
    rw [H.type hj] at hnode
    unfold SpqrTree.nEdges
    rw [(H.node j hj).ne_range, Items.nEdges_eq_zero (Bool.eq_false_iff.2 hnode)]
    omega

end RelabelAll

/-! ### Assembly -/

/-- Phase 3b: the output tree's bijections, ownership, and twin pairing, from `relabel_node_spec`.
`Bijections` needs only `Items.WF`; `Ownership` and `Twins` need `Items.OwnExtra` (and
`Ownership.ne_nvs` for R nodes needs `Items.ROriented`). -/
theorem relabelTree_own (g : Graph) (items : Items) (h : items.WF g) (hx : items.OwnExtra g)
    (hor : items.ROriented g) :
    (relabelTree g items).Bijections ∧ (relabelTree g items).Ownership ∧
      (relabelTree g items).Twins := by
  obtain ⟨idx, -, hidx, hnode⟩ := relabel_node_spec g items h
  have H : RelabelAll g items (relabelTree g items) idx := ⟨h, hidx, hnode⟩
  exact ⟨H.bijections, H.ownership hx hor, H.twins hx⟩

theorem relabelTree_bijections (g : Graph) (items : Items) (h : items.WF g) :
    (relabelTree g items).Bijections := by
  obtain ⟨idx, -, hidx, hnode⟩ := relabel_node_spec g items h
  exact (RelabelAll.mk h hidx hnode).bijections

theorem relabelTree_twins (g : Graph) (items : Items) (h : items.WF g) (hx : items.OwnExtra g) :
    (relabelTree g items).Twins := by
  obtain ⟨idx, -, hidx, hnode⟩ := relabel_node_spec g items h
  exact (RelabelAll.mk h hidx hnode).twins hx

end Spqr
