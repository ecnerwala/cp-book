import Spqr.Ranges
import Spqr.ItemTree

/-!
# `Items.Endpoints` and `Items.Shapes` from `Items.Ranges`

Pure item-level reasoning: given the item tree (`Items.Tree`), the typing facts the walk supplies
independently of any span reasoning (`TypingFacts`, a subset of `WalkTyping`) and the range/attachment
structure `Items.Ranges`, every clause of `Items.Endpoints` and `Items.Shapes` follows, except
`Shapes.canonical`, which is false for `ternarize = true` and is taken as a hypothesis
(PROOF.md §4.6).
-/

namespace Spqr
namespace Items

variable {g : Graph} {items : Items}

set_option linter.unusedSectionVars false

macro "romega" : tactic => `(tactic| ((try unfold ItemId at *); (try delta ItemId at *); (try simp only [vertItem, edgeItem, rootItem] at *); omega))

/-- The typing facts the range layer takes from `WalkTyping`. -/
structure TypingFacts (g : Graph) (items : Items) : Prop where
  i_o_leaf : ∀ i, i < items.size → items.type i = .I ∨ items.type i = .O → items.ch i = []
  vs_shape : ∀ i, i < items.size →
    match items.type i with
    | .F | .V => items.vs i = (none, none)
    | .O => ∃ v, items.vs i = (some v, none)
    | .Q => (∃ v, items.vs i = (some v, none)) ∨ (∃ u v, items.vs i = (some u, some v))
    | _ => ∃ u v, items.vs i = (some u, some v)
  vs_lt : ∀ i u, (items.vs i).1 = some u ∨ (items.vs i).2 = some u → u < g.nv

/-! ### Small facts -/

theorem type_cases (t : NodeType) :
    t = .F ∨ t = .V ∨ t = .Q ∨ t = .I ∨ t = .O ∨ t = .S ∨ t = .P ∨ t = .R := by
  cases t <;> simp

theorem pairEq_loop {a b : Nat} {q : Nat × Nat} (h : PairEq (a, b) q) (hl : q.1 = q.2) : a = b := by
  rcases h with h | h <;> simp [Prod.ext_iff] at h <;> romega

theorem pairEq_ne {a b : Nat} {q : Nat × Nat} (h : PairEq (a, b) q) (hl : q.1 ≠ q.2) : a ≠ b :=
  fun hab => hl (by rcases h with h | h <;> simp [Prod.ext_iff] at h <;> romega)

theorem pairEq_mem {a b : Nat} {q : Nat × Nat} (h : PairEq (a, b) q) :
    (a = q.1 ∨ a = q.2) ∧ (b = q.1 ∨ b = q.2) := by
  rcases h with h | h <;> simp [Prod.ext_iff] at h <;> romega

theorem parent_lt {p c : ItemId} (h : items.IsParent p c) : p < items.size := by
  by_contra hp
  simp only [IsParent, ch] at h
  rw [Array.getElem?_eq_none_iff.2 (by romega)] at h
  simp at h

theorem classify {i : ItemId} :
    i = rootItem ∨ (∃ v, v < g.nv ∧ i = vertItem v) ∨ (∃ e, e < g.ne ∧ i = edgeItem g e) ∨
      1 + g.nv + g.ne ≤ i := by
  by_cases h0 : i = 0
  · exact Or.inl h0
  by_cases h1 : i < 1 + g.nv
  · exact Or.inr (Or.inl ⟨i - 1, by romega, by romega⟩)
  by_cases h2 : i < 1 + g.nv + g.ne
  · exact Or.inr (Or.inr (Or.inl ⟨i - 1 - g.nv, by romega, by romega⟩))
  exact Or.inr (Or.inr (Or.inr (by romega)))

theorem Tree.type_V_iff' (ht : items.Tree g) {i : ItemId} (hi : i < items.size) :
    items.type i = .V ↔ ∃ v, v < g.nv ∧ i = vertItem v := by
  constructor
  · intro h
    rcases classify (g := g) (i := i) with rfl | ⟨v, hv, rfl⟩ | ⟨e, he, rfl⟩ | hn
    · rw [ht.root] at h; cases h
    · exact ⟨v, hv, rfl⟩
    · rw [ht.edge e he] at h; cases h
    · have := ht.node i hn hi; rw [h] at this; simp at this
  · rintro ⟨v, hv, rfl⟩; exact ht.vert v hv

theorem Tree.type_Q_iff' (ht : items.Tree g) {i : ItemId} (hi : i < items.size) :
    items.type i = .Q ↔ ∃ e, e < g.ne ∧ i = edgeItem g e := by
  constructor
  · intro h
    rcases classify (g := g) (i := i) with rfl | ⟨v, hv, rfl⟩ | ⟨e, he, rfl⟩ | hn
    · rw [ht.root] at h; cases h
    · rw [ht.vert v hv] at h; cases h
    · exact ⟨e, he, rfl⟩
    · have := ht.node i hn hi; rw [h] at this; simp at this
  · rintro ⟨e, he, rfl⟩; exact ht.edge e he

theorem Tree.type_ne_F' (ht : items.Tree g) {p c : ItemId} (h : items.IsParent p c) :
    items.type c ≠ .F := by
  intro hF
  have hc := ht.ch_lt p c h
  rcases classify (g := g) (i := c) with rfl | ⟨v, hv, rfl⟩ | ⟨e, he, rfl⟩ | hn
  · exact ht.root_no_parent p h
  · rw [ht.vert v hv] at hF; cases hF
  · rw [ht.edge e he] at hF; cases hF
  · have := ht.node c hn hc; rw [hF] at this; simp at this

theorem Tree.has_parent (ht : items.Tree g) {c : ItemId} (hc : c < items.size)
    (h0 : items.type c ≠ .F) : ∃ p, items.IsParent p c := by
  have : c ≠ rootItem := fun h => h0 (h ▸ ht.root)
  obtain ⟨p, hp, -⟩ := ht.unique_parent c (by unfold rootItem at this; romega) hc
  exact ⟨p, hp⟩

theorem Tree.edgeItem_lt' (ht : items.Tree g) {e : Nat} (he : e < g.ne) : edgeItem g e < items.size := by
  have := ht.size; unfold edgeItem; romega

theorem Tree.mem_vchildren (ht : items.Tree g) {i x : ItemId}
    (hx : x ∈ (items.ch i).filter fun c => items.type c = NodeType.V) :
    ∃ v, v < g.nv ∧ x = vertItem v ∧ items.IsParent i (vertItem v) := by
  have h := List.mem_filter.1 hx
  obtain ⟨v, hv, rfl⟩ := (ht.type_V_iff' (ht.ch_lt i x h.1)).1 (by simpa using h.2)
  exact ⟨v, hv, rfl, h.1⟩

theorem Tree.vchildren_nodup (ht : items.Tree g) (i : ItemId) :
    (((items.ch i).filter fun c => items.type c = NodeType.V).map (· - 1)).Nodup := by
  refine List.Nodup.map_on ?_ ((ht.ch_nodup i).filter _)
  intro x hx y hy hxy
  obtain ⟨a, -, rfl, -⟩ := ht.mem_vchildren hx
  obtain ⟨b, -, rfl, -⟩ := ht.mem_vchildren hy
  (try simp only at hxy); romega

theorem Tree.vchildren_mem (ht : items.Tree g) {i : ItemId} {x : Nat}
    (hx : x ∈ ((items.ch i).filter fun c => items.type c = NodeType.V).map (· - 1)) :
    x < g.nv ∧ items.IsParent i (vertItem x) := by
  obtain ⟨c, hc, rfl⟩ := List.mem_map.1 hx
  obtain ⟨v, hv, rfl, hp⟩ := ht.mem_vchildren hc
  unfold vertItem at *
  simpa using ⟨hv, hp⟩

theorem ch_ne_nil_of_parent {p c : ItemId} (h : items.IsParent p c) : items.ch p ≠ [] := by
  intro h0; rw [IsParent, h0] at h; simp at h

/-! ### `Endpoints` -/

section
variable (ht : items.Tree g) (hf : TypingFacts g items) {σ : List Nat} (hr : items.Ranges g σ)
include ht hf hr

theorem vs_of_type {i : ItemId} (hi : i < items.size) (h : items.type i ∈ [NodeType.S, .P, .R]) :
    ∃ u v, items.vs i = (some u, some v) ∧ u ≠ v := by
  have hs := hf.vs_shape i hi
  simp only [List.mem_cons, List.not_mem_nil, or_false] at h
  rcases h with h | h | h <;> rw [h] at hs <;>
    obtain ⟨u, v, huv⟩ := hs <;> exact ⟨u, v, huv, hr.vs_ne i hi (by simp [h]) u v huv⟩

/-- An endpoint of an S/P/R node differs from each of its V children. -/
theorem vs_ne_vchild {i : ItemId} (hi : i < items.size) (hSPR : items.type i ∈ [NodeType.S, .P, .R])
    {u x : Nat} (hu : items.IsVs i u) (hx : x < g.nv) (hpx : items.IsParent i (vertItem x)) : u ≠ x := by
  rintro rfl
  obtain ⟨e, e', -, he', -, hinc', -, hnb⟩ := hr.vs_att i hi (Or.inl hSPR) u hu
  exact hnb (((hr.interior i u hi hx hSPR).1 hpx).1.2 e' he' hinc')

theorem q_vs : ∀ e, e < g.ne → ∀ u, (items.vs (edgeItem g e)).1 = some u →
    (u = (g.edges[e]!).1 ∨ u = (g.edges[e]!).2) ∧
    ((items.vs (edgeItem g e)).2 = none ↔ items.ch (edgeItem g e) ≠ []) ∧
    ((g.edges[e]!).1 = (g.edges[e]!).2 → items.ch (edgeItem g e) ≠ []) ∧
    ∀ v, (items.vs (edgeItem g e)).2 = some v → PairEq (u, v) g.edges[e]! := by
  intro e he u hu
  by_cases h0 : items.ch (edgeItem g e) = []
  · obtain ⟨a, b, hvs, hab, hpe⟩ := hr.q_leaf e he h0
    rw [hvs] at hu ⊢
    simp only [Option.some.injEq] at hu
    subst hu
    refine ⟨(pairEq_mem hpe).1, by simp [h0], fun hl => absurd (pairEq_loop hpe hl) hab, ?_⟩
    intro v hv
    simp only [Option.some.injEq] at hv
    subst hv; exact hpe
  · obtain ⟨u', c, hvs, hinc, -⟩ := hr.q_root e he h0
    rw [hvs] at hu ⊢
    simp only [Option.some.injEq] at hu
    subst hu
    exact ⟨hinc.imp Eq.symm Eq.symm, by simp [h0], fun _ => h0, fun v hv => by simp at hv⟩

theorem q_root' : ∀ e, e < g.ne →
    (∀ c, items.ch (edgeItem g e) = [c] → (g.edges[e]!).1 = (g.edges[e]!).2) ∧
    (∀ c w u, items.ch (edgeItem g e) = [c, vertItem w] → (items.vs (edgeItem g e)).1 = some u →
      PairEq (u, w) g.edges[e]!) := by
  intro e he
  constructor
  · intro c hch
    by_contra hl
    obtain ⟨u, c', -, -, -, -, -, hnl⟩ := hr.q_root e he (by rw [hch]; simp)
    obtain ⟨w, -, -, hch', -⟩ := hnl hl
    rw [hch] at hch'; simp at hch'
  · intro c w u hch hu
    obtain ⟨u', c', hvs, -, -, -, hloop, hnl⟩ := hr.q_root e he (by rw [hch]; simp)
    rw [hvs] at hu
    simp only [Option.some.injEq] at hu
    subst hu
    by_cases hl : (g.edges[e]!).1 = (g.edges[e]!).2
    · obtain ⟨hch', -⟩ := hloop hl
      rw [hch] at hch'; simp at hch'
    · obtain ⟨w', -, hpe, hch', -⟩ := hnl hl
      rw [hch] at hch'
      simp only [List.cons.injEq, and_true] at hch'
      obtain ⟨-, hw⟩ := hch'
      have : w = w' := by romega
      subst this; exact hpe

theorem child_vs_in_parent : ∀ p c, items.IsParent p c → items.type p ∉ [NodeType.F, .V] →
    items.type c ≠ .V → ∀ u, (items.vs c).1 = some u ∨ (items.vs c).2 = some u →
      (items.vs p).1 = some u ∨ (items.vs p).2 = some u ∨ items.IsParent p (vertItem u) := by
  intro p c hpc hpF hcV u hu
  have hp := parent_lt hpc
  have hc := ht.ch_lt p c hpc
  rcases type_cases (items.type p) with h | h | h | h | h | h | h | h
  · simp [h] at hpF
  · simp [h] at hpF
  · -- block-root Q
    obtain ⟨e, he, rfl⟩ := (ht.type_Q_iff' hp).1 h
    obtain ⟨u', c', hvs, -, -, -, hloop, hnl⟩ := hr.q_root e he (ch_ne_nil_of_parent hpc)
    by_cases hl : (g.edges[e]!).1 = (g.edges[e]!).2
    · obtain ⟨hch, hvc⟩ := hloop hl
      rw [IsParent, hch] at hpc
      simp only [List.mem_singleton] at hpc
      subst hpc
      rw [hvc] at hu
      simp only [Option.some.injEq, reduceCtorEq, or_false] at hu
      subst hu; left; rw [hvs]
    · obtain ⟨w, hw, -, hch, a, b, hvc, hab⟩ := hnl hl
      rw [IsParent, hch] at hpc
      simp only [List.mem_cons, List.mem_singleton, List.not_mem_nil, or_false] at hpc
      rcases hpc with rfl | rfl
      · rw [hvc] at hu
        simp only [Option.some.injEq] at hu
        simp only [PairEq, Prod.mk.injEq] at hab
        rcases hab with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ <;> rcases hu with rfl | rfl <;>
          first
          | (left; simp [hvs]; done)
          | (right; right; simp [IsParent, hch]; done)
      · exact absurd (ht.vert w hw) hcV
  · have := hf.i_o_leaf p hp (Or.inl h); exact absurd hpc (by rw [IsParent, this]; simp)
  · have := hf.i_o_leaf p hp (Or.inr h); exact absurd hpc (by rw [IsParent, this]; simp)
  all_goals
  have hSPR : items.type p ∈ [NodeType.S, .P, .R] := by simp [h]
  have hatt : items.Att g c u := by
    rcases type_cases (items.type c) with hc' | hc' | hc' | hc' | hc' | hc' | hc' | hc'
    · exact absurd hc' (ht.type_ne_F' hpc)
    · exact absurd hc' hcV
    · obtain ⟨e, he, rfl⟩ := (ht.type_Q_iff' hc).1 hc'
      by_cases h0 : items.ch (edgeItem g e) = []
      · exact hr.vs_att _ hc (Or.inr ⟨hc', h0⟩) u hu
      · exfalso
        obtain ⟨u', c', hvs, -⟩ := hr.q_root e he h0
        obtain ⟨a, b, hab, -⟩ := hr.child_two p _ hpc hSPR hcV
        rw [hvs] at hab; cases hab
    · have := hr.io_parent p c hpc (Or.inl hc'); rw [this] at hSPR; simp at hSPR
    · have := hr.io_parent p c hpc (Or.inr hc'); rw [this] at hSPR; simp at hSPR
    · exact hr.vs_att c hc (Or.inl (by simp [hc'])) u hu
    · exact hr.vs_att c hc (Or.inl (by simp [hc'])) u hu
    · exact hr.vs_att c hc (Or.inl (by simp [hc'])) u hu
  obtain ⟨e, e', he, he', hinc, hinc', hbel, hnbel⟩ := hatt
  by_cases hout : ∃ e'', e'' < g.ne ∧ g.Inc e'' u ∧ ¬ items.EdgeBelow g p e''
  · obtain ⟨e'', he'', hinc'', hnb⟩ := hout
    have hpatt : items.Att g p u :=
      ⟨e, e'', he, he'', hinc, hinc'', Relation.ReflTransGen.head hpc hbel, hnb⟩
    rcases hr.att_vs p hp hpF u hpatt with h1 | h1
    · exact Or.inl h1
    · exact Or.inr (Or.inl h1)
  · push Not at hout
    right; right
    have hu' : u < g.nv := hf.vs_lt c u hu
    refine (hr.interior p u hp hu' hSPR).2 ⟨⟨⟨e, he, hinc⟩, fun e'' he'' hi => hout e'' he'' hi⟩, ?_⟩
    intro c' hc' hall
    by_cases hcc : c' = c
    · subst hcc; exact hnbel (hall e' he' hinc')
    · have h1 : edgeItem g e ∈ items.desc c := mem_desc.2 ⟨ht.edgeItem_lt' he, hbel⟩
      have h2 : edgeItem g e ∈ items.desc c' := mem_desc.2 ⟨ht.edgeItem_lt' he, hall e he hinc⟩
      exact Finset.disjoint_left.1 (ht.desc_disjoint hc' hpc hcc) h2 h1

theorem nv_nodup : ∀ i, i < items.size →
    ((items.vs i).1.toList ++ ((items.ch i).filter fun c => items.type c = NodeType.V).map (· - 1) ++
      (items.vs i).2.toList).Nodup := by
  intro i hi
  have hxs := ht.vchildren_nodup i
  rcases type_cases (items.type i) with h | h | h | h | h | h | h | h
  · have hvs : items.vs i = (none, none) := by have := hf.vs_shape i hi; rw [h] at this; exact this
    rw [hvs]; simpa using hxs
  · have hvs : items.vs i = (none, none) := by have := hf.vs_shape i hi; rw [h] at this; exact this
    rw [hvs]; simpa using hxs
  · obtain ⟨e, he, rfl⟩ := (ht.type_Q_iff' hi).1 h
    by_cases h0 : items.ch (edgeItem g e) = []
    · obtain ⟨a, b, hvs, hab, -⟩ := hr.q_leaf e he h0
      rw [hvs, h0]; simp [hab]
    · obtain ⟨u, c, hvs, -, hcFV, -, hloop, hnl⟩ := hr.q_root e he h0
      have hcV : items.type c ≠ .V := fun hV => hcFV (by simp [hV])
      by_cases hl : (g.edges[e]!).1 = (g.edges[e]!).2
      · obtain ⟨hch, -⟩ := hloop hl
        rw [hvs, hch]; simp [hcV]
      · obtain ⟨w, hw, hpe, hch, -⟩ := hnl hl
        have huw : u ≠ w := pairEq_ne hpe hl
        have hvw : items.type (vertItem w) = .V := ht.vert w hw
        have hw1 : vertItem w - 1 = w := by romega
        rw [hvs, hch]; simp [hcV, hvw, hw1, huw, huw.symm]
  · -- I: a bridge under a block-root Q
    have hch := hf.i_o_leaf i hi (Or.inl h)
    have hvs := hf.vs_shape i hi
    rw [h] at hvs
    obtain ⟨a, b, hvs⟩ := hvs
    rw [hvs, hch]
    suffices a ≠ b by simp [this]
    obtain ⟨p, hp⟩ := ht.has_parent hi (by rw [h]; decide)
    have hpQ := hr.io_parent p i hp (Or.inl h)
    obtain ⟨e, he, rfl⟩ := (ht.type_Q_iff' (parent_lt hp)).1 hpQ
    obtain ⟨u, c, -, -, -, -, hloop, hnl⟩ := hr.q_root e he (ch_ne_nil_of_parent hp)
    by_cases hl : (g.edges[e]!).1 = (g.edges[e]!).2
    · obtain ⟨hch', hvc⟩ := hloop hl
      rw [IsParent, hch'] at hp
      simp only [List.mem_singleton] at hp
      subst hp; rw [hvs] at hvc; cases hvc
    · obtain ⟨w, hw, hpe, hch', a', b', hvc, hab⟩ := hnl hl
      rw [IsParent, hch'] at hp
      simp only [List.mem_cons, List.mem_singleton, List.not_mem_nil, or_false] at hp
      rcases hp with rfl | rfl
      · rw [hvs] at hvc
        simp only [Prod.mk.injEq, Option.some.injEq] at hvc
        obtain ⟨rfl, rfl⟩ := hvc
        exact pairEq_ne hab (pairEq_ne hpe hl)
      · rw [ht.vert w hw] at h; cases h
  · have hch := hf.i_o_leaf i hi (Or.inr h)
    have hvs := hf.vs_shape i hi
    rw [h] at hvs
    obtain ⟨v, hvs⟩ := hvs
    rw [hvs, hch]; simp
  all_goals
  have hSPR : items.type i ∈ [NodeType.S, .P, .R] := by simp [h]
  obtain ⟨u, v, hvs, huv⟩ := vs_of_type ht hf hr hi hSPR
  have hu : items.IsVs i u := Or.inl (by rw [hvs])
  have hv : items.IsVs i v := Or.inr (by rw [hvs])
  rw [hvs]
  simp only [Option.toList_some, List.singleton_append]
  refine List.Nodup.cons ?_ (List.Nodup.append hxs (List.nodup_singleton v) ?_)
  · intro hmem
    rcases List.mem_append.1 hmem with hmem | hmem
    · obtain ⟨hx, hpx⟩ := ht.vchildren_mem hmem
      exact vs_ne_vchild ht hf hr hi hSPR hu hx hpx rfl
    · exact huv (List.mem_singleton.1 hmem)
  · intro x hx hx'
    obtain rfl := List.mem_singleton.1 hx'
    obtain ⟨hx, hpx⟩ := ht.vchildren_mem hx
    exact vs_ne_vchild ht hf hr hi hSPR hv hx hpx rfl

theorem separation : ∀ i, i < items.size → items.type i ∉ [NodeType.F, .V] →
    ∀ v e e', v < g.nv → e < g.ne → e' < g.ne →
      ((g.edges[e]!).1 = v ∨ (g.edges[e]!).2 = v) → ((g.edges[e']!).1 = v ∨ (g.edges[e']!).2 = v) →
      items.EdgeBelow g i e → ¬ items.EdgeBelow g i e' →
      (items.vs i).1 = some v ∨ (items.vs i).2 = some v := by
  intro i hi hFV v e e' _ he he' hinc hinc' hb hnb
  exact hr.att_vs i hi hFV v ⟨e, e', he, he', hinc, hinc', hb, hnb⟩

theorem interior' : ∀ i v, i < items.size → v < g.nv → items.type i ∈ [NodeType.S, .P, .R] →
    (items.IsParent i (vertItem v) ↔
      (∃ e, e < g.ne ∧ ((g.edges[e]!).1 = v ∨ (g.edges[e]!).2 = v)) ∧
      (∀ e, e < g.ne → ((g.edges[e]!).1 = v ∨ (g.edges[e]!).2 = v) → items.EdgeBelow g i e) ∧
      ∀ c, items.IsParent i c → ¬ ∀ e, e < g.ne → ((g.edges[e]!).1 = v ∨ (g.edges[e]!).2 = v) →
        items.EdgeBelow g c e) := by
  intro i v hi hv hSPR
  rw [hr.interior i v hi hv hSPR]
  constructor
  · rintro ⟨⟨h1, h2⟩, h3⟩; exact ⟨h1, h2, h3⟩
  · rintro ⟨h1, h2, h3⟩; exact ⟨⟨h1, h2⟩, h3⟩

theorem endpoints_of_ranges : items.Endpoints g where
  vs_shape := hf.vs_shape
  vs_lt := hf.vs_lt
  q_vs := q_vs ht hf hr
  q_root := q_root' ht hf hr
  child_vs_in_parent := child_vs_in_parent ht hf hr
  nv_nodup := nv_nodup ht hf hr
  separation := separation ht hf hr
  interior := interior' ht hf hr

/-! ### `Shapes` -/

theorem q_children : ∀ e, e < g.ne → items.ch (edgeItem g e) = [] ∨
    ∃ c, items.type c ∉ [NodeType.F, .V] ∧ (items.type c = .Q → items.ch c = []) ∧
      (items.ch (edgeItem g e) = [c] ∨ ∃ v, v < g.nv ∧ items.ch (edgeItem g e) = [c, vertItem v]) := by
  intro e he
  by_cases h0 : items.ch (edgeItem g e) = []
  · exact Or.inl h0
  · obtain ⟨u, c, -, -, hcFV, hcQ, hloop, hnl⟩ := hr.q_root e he h0
    refine Or.inr ⟨c, hcFV, hcQ, ?_⟩
    by_cases hl : (g.edges[e]!).1 = (g.edges[e]!).2
    · exact Or.inl (hloop hl).1
    · obtain ⟨w, hw, -, hch, -⟩ := hnl hl
      exact Or.inr ⟨w, hw, hch⟩

theorem o_parent : ∀ p c, items.IsParent p c → items.type c = .O →
    items.type p = .Q ∧ items.ch p = [c] := by
  intro p c hpc hO
  have hc := ht.ch_lt p c hpc
  have hpQ := hr.io_parent p c hpc (Or.inr hO)
  refine ⟨hpQ, ?_⟩
  obtain ⟨e, he, rfl⟩ := (ht.type_Q_iff' (parent_lt hpc)).1 hpQ
  obtain ⟨u, c', -, -, -, -, hloop, hnl⟩ := hr.q_root e he (ch_ne_nil_of_parent hpc)
  by_cases hl : (g.edges[e]!).1 = (g.edges[e]!).2
  · obtain ⟨hch, -⟩ := hloop hl
    rw [IsParent, hch] at hpc
    simp only [List.mem_singleton] at hpc
    subst hpc; exact hch
  · obtain ⟨w, hw, -, hch, a, b, hvc, -⟩ := hnl hl
    rw [IsParent, hch] at hpc
    simp only [List.mem_cons, List.mem_singleton, List.not_mem_nil, or_false] at hpc
    rcases hpc with rfl | rfl
    · have := hf.vs_shape c hc
      rw [hO] at this
      obtain ⟨v, hv⟩ := this
      rw [hvc] at hv; cases hv
    · rw [ht.vert w hw] at hO; cases hO

theorem q_leaf_of_node : ∀ p c, items.IsParent p c → items.type p ∉ [NodeType.F, .V] →
    items.type c = .Q → items.ch c = [] := by
  intro p c hpc hpFV hcQ
  have hp := parent_lt hpc
  have hc := ht.ch_lt p c hpc
  by_contra hne
  rcases type_cases (items.type p) with h | h | h | h | h | h | h | h
  · simp [h] at hpFV
  · simp [h] at hpFV
  · obtain ⟨e, he, rfl⟩ := (ht.type_Q_iff' hp).1 h
    obtain ⟨u, c', -, -, -, hcQ', hloop, hnl⟩ := hr.q_root e he (ch_ne_nil_of_parent hpc)
    have hcc : c = c' := by
      by_cases hl : (g.edges[e]!).1 = (g.edges[e]!).2
      · obtain ⟨hch, -⟩ := hloop hl
        rw [IsParent, hch] at hpc; simpa using hpc
      · obtain ⟨w, hw, -, hch, -⟩ := hnl hl
        rw [IsParent, hch] at hpc
        simp only [List.mem_cons, List.mem_singleton, List.not_mem_nil, or_false] at hpc
        rcases hpc with h1 | h1
        · exact h1
        · subst h1; rw [ht.vert w hw] at hcQ; cases hcQ
    subst hcc; exact hne (hcQ' hcQ)
  · have := hf.i_o_leaf p hp (Or.inl h); exact absurd hpc (by rw [IsParent, this]; simp)
  · have := hf.i_o_leaf p hp (Or.inr h); exact absurd hpc (by rw [IsParent, this]; simp)
  all_goals
  obtain ⟨e, he, rfl⟩ := (ht.type_Q_iff' hc).1 hcQ
  obtain ⟨u, c', hvs, -⟩ := hr.q_root e he hne
  obtain ⟨a, b, hab, -⟩ := hr.child_two p _ hpc (by simp [h]) (by rw [hcQ]; decide)
  rw [hvs] at hab; cases hab

theorem s_shape : ∀ i, i < items.size → items.type i = .S → ∃ u v xs,
    items.vs i = (some u, some v) ∧
    ((items.ch i).filter fun c => items.type c = NodeType.V) = xs.map vertItem ∧
    1 ≤ xs.length ∧ (items.virtualEdges i).Perm (List.zip (u :: xs) (xs ++ [v])) := by
  intro i hi hS
  obtain ⟨u, v, xs, h1, h2, h3, h4⟩ := hr.s_order i hi hS
  exact ⟨u, v, xs, h1, h2, h3, h4 ▸ List.Perm.refl _⟩

theorem r_shape : ∀ i, i < items.size → items.type i = .R →
    2 ≤ ((items.ch i).filter fun c => items.type c = NodeType.V).length ∧
    5 ≤ (items.virtualEdges i).length ∧
    ((items.virtualEdges i).map fun q => min q.1 q.2 + g.nv * max q.1 q.2).Nodup ∧
    (∀ q ∈ items.virtualEdges i, q.1 ≠ q.2) ∧
    (∀ u v, items.vs i = (some u, some v) → ∀ q ∈ items.virtualEdges i, ¬ PairEq q (u, v)) ∧
    ∀ c, items.IsParent i c → items.type c ≠ .V → ∃ u v, items.vs c = (some u, some v) := by
  intro i hi hR
  obtain ⟨h1, h2, h3, h4⟩ := hr.r_shape i hi hR
  have hSPR : items.type i ∈ [NodeType.S, .P, .R] := by simp [hR]
  have two : ∀ c, items.IsParent i c → items.type c ≠ .V →
      ∃ a b, items.vs c = (some a, some b) ∧ a ≠ b :=
    fun c hc hV => hr.child_two i c hc hSPR hV
  refine ⟨h1, h2, h3, ?_, h4, fun c hc hV => (two c hc hV).imp fun _ h => h.imp fun _ h => h.1⟩
  intro q hq
  obtain ⟨c, hc, rfl⟩ := List.mem_map.1 hq
  have hc' := List.mem_filter.1 hc
  obtain ⟨a, b, hab, hne⟩ := two c hc'.1 (by simpa using hc'.2)
  simp [hab, hne]

theorem shapes_of_ranges
    (hcanon : ∀ p c, items.IsParent p c →
      (items.type c = .S → items.type p ≠ .S) ∧ (items.type c = .P → items.type p ≠ .P)) :
    items.Shapes g where
  i_o_leaf := hf.i_o_leaf
  q_children := q_children ht hf hr
  o_parent := o_parent ht hf hr
  q_leaf_of_node := q_leaf_of_node ht hf hr
  p_shape := hr.p_shape
  s_shape := s_shape ht hf hr
  s_order := fun i hi hS => by
    obtain ⟨u, v, xs, h1, h2, -, h4⟩ := hr.s_order i hi hS
    exact ⟨u, v, xs, h1, h2, h4⟩
  r_shape := r_shape ht hf hr
  canonical := hcanon

/-- `Items.WF` from the item tree, the typing facts, the ranges, and canonicity (the only clause
not implied by `Ranges`: it fails for `ternarize = true`). -/
theorem wf_of_ranges
    (hcanon : ∀ p c, items.IsParent p c →
      (items.type c = .S → items.type p ≠ .S) ∧ (items.type c = .P → items.type p ≠ .P)) :
    items.WF g :=
  ⟨ht, endpoints_of_ranges ht hf hr, shapes_of_ranges ht hf hr hcanon⟩

end

end Items
end Spqr
