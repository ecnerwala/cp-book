import Spqr.RelabelRep
import Spqr.PieceSep
import Spqr.Ranges

/-!
# Relabel: `SpqrTree.PieceSep` from the final items

`spqrTree_pieceSep` (`WalkPieceSep.lean`) is `Items.WF` + `Items.Ranges` + two walk-side facts about
the final items (`Items.QUpper`, `Items.RootSep`), transported through `RelabelOK` (`RelabelRep.lean`).
-/

namespace Spqr.Items

variable (g : Graph) (items : Items)

/-- A child of a V item records that vertex as its first endpoint (`finishBoundary` writes
`vs := (some curV, none)` and appends the block root to `vertItem curV`). -/
def QUpper : Prop :=
  ∀ v c, v < g.nv → items.IsParent (vertItem v) c → (items.vs c).1 = some v

/-- Distinct children of the root (the DFS components) share no vertex. -/
def RootSep : Prop :=
  ∀ a b, items.IsParent rootItem a → items.IsParent rootItem b → a ≠ b →
    ∀ v e e', e < g.ne → e' < g.ne → g.Inc e v → g.Inc e' v →
      items.EdgeBelow g a e → items.EdgeBelow g b e' → False

/-- Every child of the root is a V item (a DFS root's edges all lie in its blocks). -/
def RootV : Prop := ∀ c, items.IsParent rootItem c → items.type c = .V

/-- The non-V child of a block-root Q is oriented like the edge: `vs c = (u, some w)` for the
children `[c, vertItem w]`, `vs c = (u, none)` for a loop's `[c]`. -/
def QChildVs : Prop :=
  ∀ e, e < g.ne →
    (∀ c w, items.ch (edgeItem g e) = [c, vertItem w] →
      items.vs c = ((items.vs (edgeItem g e)).1, some w)) ∧
    ∀ c, items.ch (edgeItem g e) = [c] → items.vs c = ((items.vs (edgeItem g e)).1, none)

/-- The non-V children of a P item all carry the item's own `vs`, in the same order. -/
def PChildVs : Prop :=
  ∀ i c, i < items.size → items.type i = .P → items.IsParent i c → items.type c ≠ .V →
    items.vs c = items.vs i

end Spqr.Items

namespace Spqr.RelabelOK

open Items

variable {g : Graph} {items : Items} {t : SpqrTree} {idx : ItemId → Nat}
  (h : RelabelOK g items t idx)
include h

theorem children_eq {i : ItemId} (hi : i < items.size) (hR : items.type i ≠ .R) :
    t.children (idx i) = (items.ch i).map idx := by
  obtain ⟨pos, hl⟩ := (h.node i hi).layout
  rw [hl.children, Items.ordered, if_pos hR]

theorem origId_vertItem {v : Nat} (hv : v < g.nv) : t.origId[idx (vertItem v)]! = some v := by
  rw [(h.node _ (h.vertItem_lt hv)).orig, Items.origOf, h.type_vertItem hv]
  simp [vertItem]

theorem origId_edgeItem {e : Nat} (he : e < g.ne) : t.origId[idx (edgeItem g e)]! = some e := by
  rw [(h.node _ (h.edgeItem_lt he)).orig, Items.origOf, h.type_edgeItem he]
  dsimp only
  congr 1
  riomega

theorem edgeFlipped_eq {e : Nat} (he : e < g.ne) :
    t.edgeFlipped[e]! = ((items.vs (edgeItem g e)).1 != some (g.edges[e]!).1) := by
  have := ((h.node _ (h.edgeItem_lt he)).edge_index (h.type_edgeItem he)).2
  have hidx : edgeItem g e - 1 - g.nv = e := by riomega
  rw [hidx] at this
  exact this

theorem vertItem_of_V {i : ItemId} (hV : items.type i = .V) :
    ∃ v, v < g.nv ∧ i = vertItem v := by
  obtain ⟨h1, h2⟩ := h.tree.V_range hV
  exact ⟨i - 1, by omega, by riomega⟩

theorem edgeItem_of_Q {i : ItemId} (hi : i < items.size) (hQ : items.type i = .Q) :
    ∃ e, e < g.ne ∧ i = edgeItem g e := h.type_Q_eq hi hQ

theorem touches_iff {a : ItemId} (ha : a < items.size) (v : Nat) :
    t.Touches g (idx a) v ↔ ∃ e, e < g.ne ∧ g.Inc e v ∧ items.EdgeBelow g a e := by
  unfold SpqrTree.Touches SpqrTree.Graph.Incident Graph.Inc
  constructor
  · rintro ⟨e, ⟨he, hinc⟩, hin⟩
    exact ⟨e, he, hinc, (h.edgeIn_iff ha he).1 hin⟩
  · rintro ⟨e, he, hinc, hb⟩
    exact ⟨e, ⟨he, hinc⟩, (h.edgeIn_iff ha he).2 hb⟩

theorem parent_some {c : ItemId} (hc : c < items.size) {p : Nat} (hp : t.parent (idx c) = some p) :
    ∃ p', p' < items.size ∧ idx p' = p ∧ items.IsParent p' c := by
  obtain ⟨c', p', hc', hp', hcc, rfl, hpar⟩ := h.parent_cases hp
  rw [h.idx_inj hc' hc hcc] at hpar
  exact ⟨p', hp', rfl, hpar⟩

theorem mem_children {i c : ItemId} (hi : i < items.size) (hc : c ∈ items.ch i) :
    idx c ∈ t.children (idx i) := (h.mem_children_iff hi _).2 ⟨c, hc, rfl⟩

theorem type_ne_of_mem_children {i : ItemId} (hi : i < items.size) {n : Nat}
    (hn : n ∈ t.children (idx i)) : ∃ c ∈ items.ch i, idx c = n := (h.mem_children_iff hi n).1 hn

omit h in
theorem map_eq_pair {α β : Type} {f : α → β} {l : List α} {a b : β} (hl : l.map f = [a, b]) :
    ∃ a' b', l = [a', b'] ∧ f a' = a ∧ f b' = b := by
  cases l with
  | nil => cases hl
  | cons a' l =>
    cases l with
    | nil => cases hl
    | cons b' l =>
      cases l with
      | nil =>
        simp only [List.map_cons, List.map_nil, List.cons.injEq, and_true] at hl
        exact ⟨a', b', rfl, hl.1, hl.2⟩
      | cons _ _ => cases hl

omit h in
theorem map_eq_single {α β : Type} {f : α → β} {l : List α} {a : β} (hl : l.map f = [a]) :
    ∃ a', l = [a'] ∧ f a' = a := by
  cases l with
  | nil => cases hl
  | cons a' l =>
    cases l with
    | nil =>
      simp only [List.map_cons, List.map_nil, List.cons.injEq, and_true] at hl
      exact ⟨a', rfl, hl⟩
    | cons _ _ => cases hl

theorem root_of_F {i : ItemId} (hi : i < items.size) (hF : items.type i = .F) : i = rootItem := by
  by_contra h0
  by_cases h1 : i < 1 + g.nv
  · have : i = vertItem (i - 1) := by riomega
    rw [this, h.type_vertItem (by riomega)] at hF; cases hF
  · by_cases h2 : i < 1 + g.nv + g.ne
    · have : i = edgeItem g (i - 1 - g.nv) := by riomega
      rw [this, h.type_edgeItem (by riomega)] at hF; cases hF
    · have := h.tree.node i (by riomega) hi
      simp [hF] at this

/-- A non-F/V item with `vs = (some u, none)` that is a leaf if Q is an `O` leaf. -/
theorem type_O_of_vs_none {c : ItemId} (hcs : c < items.size) (hcF : items.type c ≠ .F)
    (hcV : items.type c ≠ .V) (hcq : items.type c = .Q → items.ch c = []) {u : Nat}
    (hvs : items.vs c = (some u, none)) : items.type c = .O := by
  have hs := h.endpoints.vs_shape c hcs
  cases hty : items.type c with
  | F => exact absurd hty hcF
  | V => exact absurd hty hcV
  | O => rfl
  | Q =>
    obtain ⟨e, he, rfl⟩ := h.type_Q_eq hcs hty
    have := ((h.endpoints.q_vs e he u (by rw [hvs])).2.1).1 (by rw [hvs])
    exact absurd (hcq hty) this
  | I => rw [hty, hvs] at hs; obtain ⟨_, _, hs⟩ := hs; cases hs
  | S => rw [hty, hvs] at hs; obtain ⟨_, _, hs⟩ := hs; cases hs
  | P => rw [hty, hvs] at hs; obtain ⟨_, _, hs⟩ := hs; cases hs
  | R => rw [hty, hvs] at hs; obtain ⟨_, _, hs⟩ := hs; cases hs

omit h in
theorem hasCap_of_node {c : ItemId} (hcF : items.type c ≠ .F) (hcV : items.type c ≠ .V)
    (hcq : items.type c = .Q → items.ch c = []) : items.hasCap c = true := by
  cases hty : items.type c with
  | F => exact absurd hty hcF
  | V => exact absurd hty hcV
  | Q => simp [Items.hasCap, hty, hcq hty, NodeType.isNode]
  | _ => simp [Items.hasCap, hty, NodeType.isNode]

omit h in
theorem isNode_of_hasCap {c : ItemId} (hc : items.hasCap c = true) : (items.type c).isNode = true := by
  simp only [Items.hasCap, Bool.and_eq_true] at hc; exact hc.1

omit h in
theorem type_ne_V_of_isNode {c : ItemId} (hn : (items.type c).isNode = true) : items.type c ≠ .V := by
  intro hV; rw [hV] at hn; cases hn

omit h in
theorem type_ne_F_of_isNode {c : ItemId} (hn : (items.type c).isNode = true) : items.type c ≠ .F := by
  intro hV; rw [hV] at hn; cases hn

theorem capNe_some {c : ItemId} (hcs : c < items.size) {ne : Nat} (hcap : t.capNe (idx c) = some ne) :
    items.hasCap c = true ∧ ne = (t.neRange (idx c)).1 := by
  unfold SpqrTree.capNe at hcap
  split at hcap
  · rename_i hc; rw [h.hasCap_eq hcs] at hc; cases hcap; exact ⟨hc, rfl⟩
  · cases hcap

theorem q_not_R {e : Nat} (he : e < g.ne) : items.type (edgeItem g e) ≠ .R := by
  rw [h.type_edgeItem he]; decide

/-! ### The `Q` fields -/

section Q
variable (hr : Items.Ranges g items σ)
include hr

theorem pieceSep_q_shape : ∀ i, i < t.size → t.type i = .Q →
    t.children i = [] ∨ (∃ c, t.type c = .O ∧ t.children i = [c]) ∨
      ∃ c w, t.type c ≠ .V ∧ t.type c ≠ .O ∧ t.hasCap c = true ∧
        t.type w = .V ∧ t.children i = [c, w] := by
  intro n hn hQ
  obtain ⟨i, hi, rfl⟩ := h.idx_surj hn
  rw [h.type_eq hi] at hQ
  obtain ⟨e, he, rfl⟩ := h.type_Q_eq hi hQ
  rw [h.children_eq hi (h.q_not_R he)]
  rcases h.shapes.q_children e he with h0 | ⟨c, hc, hcq, hch⟩
  · left; rw [h0]; rfl
  have hcm : c ∈ items.ch (edgeItem g e) := by rcases hch with h1 | ⟨v, -, h1⟩ <;> simp [h1]
  have hcs := h.ch_lt hi hcm
  have hcV : items.type c ≠ .V := fun hV => hc (by simp [hV])
  have hcF : items.type c ≠ .F := fun hV => hc (by simp [hV])
  rcases hch with h1 | ⟨v, hv, h1⟩
  · right; left
    refine ⟨idx c, ?_, by rw [h1]; rfl⟩
    rw [h.type_eq hcs]
    have hloop := (h.endpoints.q_root e he).1 c h1
    obtain ⟨u, c', -, -, -, -, hl, -⟩ := hr.q_root e he (by rw [h1]; simp)
    obtain ⟨hch', hvsc⟩ := hl hloop
    rw [h1] at hch'
    simp only [List.cons.injEq, and_true] at hch'
    subst hch'
    exact h.type_O_of_vs_none hcs hcF hcV hcq hvsc
  · right; right
    refine ⟨idx c, idx (vertItem v), ?_, ?_, ?_, ?_, by rw [h1]; rfl⟩
    · rw [h.type_eq hcs]; exact hcV
    · rw [h.type_eq hcs]
      intro hO
      have := (h.shapes.o_parent _ _ hcm hO).2
      rw [h1] at this; cases this
    · rw [h.hasCap_eq hcs]; exact hasCap_of_node hcF hcV hcq
    · rw [h.type_eq (h.vertItem_lt hv)]; exact h.type_vertItem hv

theorem pieceSep_q_leaf_parent (hrv : items.RootV) :
    ∀ i p, i < t.size → t.type i = .Q → t.children i = [] → t.parent i = some p →
      t.type p ≠ .F ∧ t.type p ≠ .V := by
  intro n p hn hQ hch hp
  obtain ⟨i, hi, rfl⟩ := h.idx_surj hn
  rw [h.type_eq hi] at hQ
  obtain ⟨e, he, rfl⟩ := h.type_Q_eq hi hQ
  rw [h.children_eq hi (h.q_not_R he), List.map_eq_nil_iff] at hch
  obtain ⟨p', hp's, rfl, hpar⟩ := h.parent_some hi hp
  rw [h.type_eq hp's]
  constructor
  · intro hF
    have := h.root_of_F hp's hF
    subst this
    have := hrv _ hpar
    rw [hQ] at this; cases this
  · intro hV
    obtain ⟨v, hv, rfl⟩ := h.vertItem_of_V hV
    exact hr.q_under_v v _ hv hpar hch

theorem pieceSep_q_loop : ∀ i e, i < t.size → t.type i = .Q → t.origId[i]! = some e →
    ((g.edges[e]!).1 = (g.edges[e]!).2 ↔ ∃ c ∈ t.children i, t.type c = .O) := by
  intro n e hn hQ horig
  obtain ⟨i, hi, rfl⟩ := h.idx_surj hn
  rw [h.type_eq hi] at hQ
  obtain ⟨e', he', rfl⟩ := h.type_Q_eq hi hQ
  rw [h.origId_edgeItem he'] at horig
  cases horig
  rw [h.children_eq hi (h.q_not_R he')]
  constructor
  · intro hl
    obtain ⟨u, -, hu, -⟩ := h.q_pair he'
    have hne := (h.endpoints.q_vs e he' u hu).2.2.1 hl
    obtain ⟨u', c, -, -, hcF, hcq, hloop, -⟩ := hr.q_root e he' hne
    obtain ⟨hch, hvsc⟩ := hloop hl
    have hcm : c ∈ items.ch (edgeItem g e) := by rw [hch]; simp
    refine ⟨idx c, List.mem_map_of_mem hcm, ?_⟩
    rw [h.type_eq (h.ch_lt hi hcm)]
    exact h.type_O_of_vs_none (h.ch_lt hi hcm) (fun hF => hcF (by simp [hF]))
      (fun hV => hcF (by simp [hV])) hcq hvsc
  · rintro ⟨c, hcm, hO⟩
    obtain ⟨c', hc', rfl⟩ := List.mem_map.1 hcm
    rw [h.type_eq (h.ch_lt hi hc')] at hO
    exact (h.endpoints.q_root e he').1 c' (h.shapes.o_parent _ _ hc' hO).2

theorem pieceSep_q_upper (hqu : items.QUpper g) :
    ∀ i p e, i < t.size → t.type i = .Q → t.parent i = some p → t.type p = .V →
      t.origId[i]! = some e →
      t.origId[p]! = some (if t.edgeFlipped[e]! then (g.edges[e]!).2 else (g.edges[e]!).1) := by
  intro n p e hn hQ hp hpV horig
  obtain ⟨i, hi, rfl⟩ := h.idx_surj hn
  rw [h.type_eq hi] at hQ
  obtain ⟨e', he', rfl⟩ := h.type_Q_eq hi hQ
  rw [h.origId_edgeItem he'] at horig
  cases horig
  obtain ⟨p', hp's, rfl, hpar⟩ := h.parent_some hi hp
  rw [h.type_eq hp's] at hpV
  obtain ⟨v, hv, rfl⟩ := h.vertItem_of_V hpV
  rw [h.origId_vertItem hv, h.edgeFlipped_eq he']
  have hu := hqu v _ hv hpar
  rw [hu]
  rcases (h.endpoints.q_vs e he' v hu).1 with h1 | h1
  · rw [h1]; simp
  · by_cases h12 : (g.edges[e]!).1 = (g.edges[e]!).2
    · rw [h1, ← h12]; simp
    · rw [h1]
      have : (some (g.edges[e]!).2 != some (g.edges[e]!).1) = true := by
        simp [bne_iff_ne, Ne.symm h12]
      rw [if_pos this]

/-- The children of a block-root `Q` with two children: `[c, vertItem w]`, `{u, w}` the edge. -/
theorem q_two_children {e : Nat} (he : e < g.ne) {c w : ItemId}
    (hch : items.ch (edgeItem g e) = [c, w]) :
    ∃ u w', w = vertItem w' ∧ w' < g.nv ∧ items.vs (edgeItem g e) = (some u, none) ∧
      (g.edges[e]!).1 ≠ (g.edges[e]!).2 ∧ Items.PairEq (u, w') g.edges[e]! ∧
      items.type c ∉ [NodeType.F, .V] ∧ (items.type c = .Q → items.ch c = []) := by
  obtain ⟨u, c0, hvs, -, hcF, hcq, hloop, hnl⟩ := hr.q_root e he (by rw [hch]; simp)
  by_cases h12 : (g.edges[e]!).1 = (g.edges[e]!).2
  · obtain ⟨h1, -⟩ := hloop h12; rw [hch] at h1; cases h1
  · obtain ⟨w', hw', hpq, h1, -⟩ := hnl h12
    rw [hch] at h1
    simp only [List.cons.injEq, and_true] at h1
    obtain ⟨rfl, rfl⟩ := h1
    exact ⟨u, w', rfl, hw', hvs, h12, hpq, hcF, hcq⟩

theorem pieceSep_q_lower : ∀ i c w e, i < t.size → t.type i = .Q → t.children i = [c, w] →
    t.origId[i]! = some e →
    t.origId[w]! = some (if t.edgeFlipped[e]! then (g.edges[e]!).1 else (g.edges[e]!).2) := by
  intro n c w e hn hQ hch horig
  obtain ⟨i, hi, rfl⟩ := h.idx_surj hn
  rw [h.type_eq hi] at hQ
  obtain ⟨e', he', rfl⟩ := h.type_Q_eq hi hQ
  rw [h.origId_edgeItem he'] at horig
  cases horig
  rw [h.children_eq hi (h.q_not_R he')] at hch
  obtain ⟨c', w', hch', rfl, rfl⟩ := map_eq_pair hch
  obtain ⟨u, w, rfl, hw, hvs, h12, hpq, -, -⟩ := h.q_two_children hr he' hch'
  rw [h.origId_vertItem hw, h.edgeFlipped_eq he', hvs]
  rcases hpq with h2 | h2
  · rw [← h2]; simp
  · simp only [Prod.mk.injEq] at h2
    obtain ⟨rfl, rfl⟩ := h2
    have : (some (g.edges[e]!).2 != some (g.edges[e]!).1) = true := by
      simp [bne_iff_ne, Ne.symm h12]
    rw [if_pos this]

theorem pieceSep_q_cap_orient (hqc : items.QChildVs g) :
    ∀ i c e ne p, i < t.size → t.type i = .Q → c ∈ t.children i →
      t.origId[i]! = some e → t.capNe c = some ne → t.neOrig ne = some p →
      p = if t.edgeFlipped[e]! then ((g.edges[e]!).2, (g.edges[e]!).1) else g.edges[e]! := by
  intro n c e ne p hn hQ hcm horig hcap hor
  obtain ⟨i, hi, rfl⟩ := h.idx_surj hn
  rw [h.type_eq hi] at hQ
  obtain ⟨e', he', rfl⟩ := h.type_Q_eq hi hQ
  rw [h.origId_edgeItem he'] at horig
  cases horig
  rw [h.children_eq hi (h.q_not_R he')] at hcm
  obtain ⟨c', hc', rfl⟩ := List.mem_map.1 hcm
  have hcs := h.ch_lt hi hc'
  obtain ⟨hcap', rfl⟩ := h.capNe_some hcs hcap
  have hn' := isNode_of_hasCap hcap'
  have hcV := type_ne_V_of_isNode hn'
  obtain ⟨x, y, hx, hy, hcq, -⟩ := h.child_cap_orig hi (by rw [hQ]; rfl) hc' hcV
  rw [hcq] at hor
  cases hor
  obtain ⟨u, c0, hvs, hinc, -, -, hloop, hnl⟩ := hr.q_root e he' (by intro h0; rw [h0] at hc'; cases hc')
  rw [h.edgeFlipped_eq he', hvs]
  by_cases h12 : (g.edges[e]!).1 = (g.edges[e]!).2
  · obtain ⟨h1, -⟩ := hloop h12
    have hvc := (hqc e he').2 _ h1
    rw [h1] at hc'
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hc'
    subst hc'
    rw [hvs] at hvc
    rw [hvc] at hx hy
    simp only [Option.some.injEq, Option.getD_none] at hx hy
    subst hx; subst hy
    have hu : u = (g.edges[e]!).1 := by rcases hinc with h3 | h3 <;> omega
    subst hu
    generalize g.edges[e]! = pr at *
    obtain ⟨p1, p2⟩ := pr
    simp only at h12
    subst h12
    split <;> rfl
  · obtain ⟨w, hw, hpq, h1, -⟩ := hnl h12
    have hvc := (hqc e he').1 _ _ h1
    rw [h1] at hc'
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hc'
    rcases hc' with rfl | rfl
    · rw [hvs] at hvc
      rw [hvc] at hx hy
      simp only [Option.some.injEq, Option.getD_some] at hx hy
      subst hx; subst hy
      rcases hpq with h2 | h2
      · rw [← h2]; simp
      · simp only [Prod.mk.injEq] at h2
        obtain ⟨rfl, rfl⟩ := h2
        have : (some (g.edges[e]!).2 != some (g.edges[e]!).1) = true := by
          simp [bne_iff_ne, Ne.symm h12]
        rw [if_pos this]
    · exact absurd (h.type_vertItem hw) hcV

end Q

theorem pieceSep_cap_orig : ∀ c ne, c < t.size → t.capNe c = some ne → ∃ p, t.neOrig ne = some p := by
  intro n ne hn hcap
  obtain ⟨c, hc, rfl⟩ := h.idx_surj hn
  obtain ⟨hcap', rfl⟩ := h.capNe_some hc hcap
  obtain ⟨x, -, hx⟩ := h.cap_orig hc (isNode_of_hasCap hcap') hcap'
  exact ⟨_, hx⟩

end Spqr.RelabelOK
