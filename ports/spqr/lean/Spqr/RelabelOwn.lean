import Mathlib.Data.List.Sort
import Mathlib.Data.List.Nodup
import Spqr.RelabelSpec
import Spqr.ItemTree

/-!
# Phase 3b: `Bijections`, `Ownership`, `Twins` of `relabelTree` from the per-node interface

Everything here is derived from `relabel_node_spec` (taken as a hypothesis through its `∃ idx`)
and `Items.WF`, plus the extra item-level hypotheses collected in `Items.OwnExtra` that
`Items.WF` does not promise (see its docstring), and the `layoutNode` edge facts
`layoutNode_edges` (a `Layout`-local statement in the scope of the `LayoutShape` work).
-/

namespace Spqr

/-- `iomega` after unfolding the `ItemId` abbreviation (`iomega` does not see through it). -/
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
      (items.nvList g i).length ≤ 2) := by
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
      refine ⟨fun _ => by omega, fun _ => by omega, by simp, by simp, fun _ _ => ?_⟩
      have := Option.toList_length_le (items.vs i).1
      have := Option.toList_length_le (items.vs i).2
      omega
    · have hcge : 1 + g.nv ≤ c := ht.child_ge_of_ne_V (by rw [h1]; simp) (by simp at hc; tauto)
      have hcapv : items.capCount i = 0 := by
        simp [Items.capCount, Items.hasCap, hty, h1]
      simp only [h1, List.countP_cons, List.countP_nil] at hne
      simp [hcge, hcapv] at hne
      refine ⟨fun _ => by omega, fun _ => by omega, by simp, by simp, ?_⟩
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
      refine ⟨fun h => by omega, fun _ => by omega, by simp, by simp, ?_⟩
      intro hcp; simp [Items.hasCap, hty, h1] at hcp
  · -- I
    have h0 := hw.shapes.i_o_leaf i hi (Or.inl hty)
    obtain ⟨u, v, huv⟩ := hw.vs_two hi (Or.inl hty)
    simp only [h0, List.countP_nil, List.filter_nil, List.length_nil, huv, Option.toList_some,
      List.length_singleton] at hne hlen
    refine ⟨fun _ => by omega, fun _ => by omega, by simp, by simp, fun _ _ => by omega⟩
  · -- O
    have h0 := hw.shapes.i_o_leaf i hi (Or.inr hty)
    simp only [h0, List.countP_nil] at hne
    refine ⟨fun _ => by omega, by simp, by simp, by simp, by simp⟩
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
    refine ⟨fun _ => by omega, by simp, fun _ => by omega, by simp, by simp⟩
  · -- P
    obtain ⟨_, hnoV, _⟩ := hw.shapes.p_shape i hi hty
    obtain ⟨u, v, huv⟩ := hw.vs_two hi (Or.inr (Or.inr (Or.inl hty)))
    have hV : ((items.ch i).filter (· < 1 + g.nv)) = [] := by
      rw [← ht.filter_V_eq, List.filter_eq_nil_iff]
      intro c hc; simpa using hnoV c hc
    rw [huv, hV] at hlen; simp only [Option.toList_some, List.length_singleton, List.length_nil] at hlen
    refine ⟨fun _ => by omega, by simp, by simp, by simp, fun _ _ => by omega⟩
  · -- R
    obtain ⟨h2, _, _, _⟩ := hw.shapes.r_shape i hi hty
    obtain ⟨u, v, huv⟩ := hw.vs_two hi (Or.inr (Or.inr (Or.inr hty)))
    rw [ht.filter_V_eq] at h2
    rw [huv] at hlen; simp only [Option.toList_some, List.length_singleton] at hlen
    refine ⟨fun _ => by omega, by simp, by simp, fun _ => by omega, by simp⟩

end Items

/-! ### `layoutNode` edges (a `Layout`-local fact, in the scope of `LayoutShape`) -/

/-- Endpoints (node-vert ids) of edge `k` of `layoutNode`. -/
def layoutNvs (ty : NodeType) (nvSt nvEn : Nat) (ec : List (Nat × Nat)) (k : Nat) : Nat × Nat :=
  if nvEn - nvSt = 1 then (nvSt, nvSt)
  else if ty = .Q ∨ ty = .I ∨ ty = .P then (nvSt, nvSt + 1)
  else if ty = .S then (if k = 0 then (nvSt, nvEn - 1) else (nvSt + k - 1, nvSt + k))
  else (if k = 0 then (nvSt, nvEn - 1) else ec[k - 1]?.getD (0, 0))

/-- Edge `k` of `layoutNode` for a node type: owned by `node`, no twin, endpoints `layoutNvs`.
Stated `Layout`-locally (session C's territory); admitted here. -/
theorem layoutNode_edges (ty : NodeType) (node nvSt nvEn neSt neEn : Nat) (ec : List (Nat × Nat))
    (hty : ty.isNode = true) (h1 : nvEn - nvSt = 1 → neEn - neSt ≤ 1)
    (hQI : ty = .Q ∨ ty = .I → neEn - neSt ≤ 1)
    (hR : ty = .R → nvSt + 2 ≤ nvEn ∧ neEn = neSt + ec.length + 1)
    (k : Nat) (hk : k < neEn - neSt) :
    (layoutNode ty node nvSt nvEn neSt neEn ec).edges[k]! =
      ⟨node, none, layoutNvs ty nvSt nvEn ec k⟩ := by
  sorry

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

end RelabelAll

end Spqr
