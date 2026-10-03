import Spqr.RelabelRep
import Spqr.RelabelOwn
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
  rw [hl.children, Items.ordered, ite_eq_left hR]

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
      rw [ite_eq_left this]

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
    rw [ite_eq_left this]

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
        rw [ite_eq_left this]
    · exact absurd (h.type_vertItem hw) hcV

end Q

theorem pieceSep_cap_orig : ∀ c ne, c < t.size → t.capNe c = some ne → ∃ p, t.neOrig ne = some p := by
  intro n ne hn hcap
  obtain ⟨c, hc, rfl⟩ := h.idx_surj hn
  obtain ⟨hcap', rfl⟩ := h.capNe_some hc hcap
  obtain ⟨x, -, hx⟩ := h.cap_orig hc (isNode_of_hasCap hcap') hcap'
  exact ⟨_, hx⟩

/-! ### Node fields -/

omit h in
theorem isNode_of_SPR {ty : NodeType} (hT : ty = .S ∨ ty = .P ∨ ty = .R) : ty.isNode = true := by
  rcases hT with rfl | rfl | rfl <;> rfl

omit h in
theorem SPR_mem {ty : NodeType} (hT : ty = .S ∨ ty = .P ∨ ty = .R) : ty ∈ [NodeType.S, .P, .R] := by
  rcases hT with rfl | rfl | rfl <;> simp

theorem not_root_of_child {p c : ItemId} (hc : items.IsParent p c) : c ≠ rootItem := by
  rintro rfl
  exact h.tree.root_no_parent p hc

theorem pieceSep_node_parent (hrv : items.RootV) :
    ∀ i p, i < t.size → (t.type i = .S ∨ t.type i = .P ∨ t.type i = .R) →
      t.parent i = some p → t.type p ≠ .F ∧ t.type p ≠ .V := by
  intro n p hn hT hp
  obtain ⟨i, hi, rfl⟩ := h.idx_surj hn
  rw [h.type_eq hi] at hT
  obtain ⟨p', hp's, rfl, hpar⟩ := h.parent_some hi hp
  rw [h.type_eq hp's]
  constructor
  · intro hF
    have := h.root_of_F hp's hF
    subst this
    have := hrv _ hpar
    rcases hT with hT | hT | hT <;> rw [hT] at this <;> cases this
  · intro hV
    obtain ⟨v, hv, rfl⟩ := h.vertItem_of_V hV
    have := h.tree.v_children v _ hv hpar
    rcases hT with hT | hT | hT <;> rw [hT] at this <;> cases this

theorem pieceSep_node_child_cap (hr : Items.Ranges g items σ) :
    ∀ i c, i < t.size → (t.type i = .S ∨ t.type i = .P ∨ t.type i = .R) →
      c ∈ t.children i → t.type c ≠ .V → t.hasCap c = true ∧ t.type c ≠ .I ∧ t.type c ≠ .O := by
  intro n c hn hT hcm hcV
  obtain ⟨i, hi, rfl⟩ := h.idx_surj hn
  rw [h.type_eq hi] at hT
  obtain ⟨c', hc', rfl⟩ := h.type_ne_of_mem_children hi hcm
  have hcs := h.ch_lt hi hc'
  rw [h.type_eq hcs] at hcV ⊢
  rw [h.hasCap_eq hcs]
  have hIO : ∀ ty, items.type c' = ty → ty = .I ∨ ty = .O → False := by
    intro ty hty hio
    have := hr.io_parent i c' hc' (by rw [hty]; exact hio)
    rcases hT with hT | hT | hT <;> rw [hT] at this <;> cases this
  refine ⟨?_, fun hI => hIO _ hI (Or.inl rfl), fun hO => hIO _ hO (Or.inr rfl)⟩
  refine hasCap_of_node ?_ hcV fun hQ => h.shapes.q_leaf_of_node i c' hc' ?_ hQ
  · intro hF
    exact h.not_root_of_child hc' (h.root_of_F hcs hF)
  · rcases hT with hT | hT | hT <;> rw [hT] <;> decide

theorem edgeIn_self {e : Nat} (he : e < g.ne) : t.EdgeIn (idx (edgeItem g e)) e :=
  (h.edgeIn_iff (h.edgeItem_lt he) he).2 Relation.ReflTransGen.refl

theorem pieceSep_cap_nonempty (hr : Items.Ranges g items σ) :
    ∀ c, c < t.size → t.hasCap c = true → t.type c ≠ .I → t.type c ≠ .O → ∃ e, t.EdgeIn c e := by
  intro n hn hcap hI hO
  obtain ⟨c, hc, rfl⟩ := h.idx_surj hn
  rw [h.hasCap_eq hc] at hcap
  rw [h.type_eq hc] at hI hO
  have hn' := isNode_of_hasCap hcap
  cases hty : items.type c with
  | F => rw [hty] at hn'; cases hn'
  | V => rw [hty] at hn'; cases hn'
  | I => exact absurd hty hI
  | O => exact absurd hty hO
  | Q =>
    obtain ⟨e, he, rfl⟩ := h.type_Q_eq hc hty
    exact ⟨e, h.edgeIn_self he⟩
  | S | P | R =>
    obtain ⟨u, hu⟩ := h.vs_fst hc hn'
    obtain ⟨e, -, he, -, -, -, hb, -⟩ := hr.vs_att c hc (Or.inl (by rw [hty]; simp)) u (Or.inl hu)
    exact ⟨e, (h.edgeIn_iff hc he).2 hb⟩

theorem pieceSep_v_nonempty : ∀ i, i < t.size → t.type i = .V →
    ∀ c ∈ t.children i, ∃ e, t.EdgeIn c e := by
  intro n hn hV c hcm
  obtain ⟨i, hi, rfl⟩ := h.idx_surj hn
  rw [h.type_eq hi] at hV
  obtain ⟨v, hv, rfl⟩ := h.vertItem_of_V hV
  obtain ⟨c', hc', rfl⟩ := h.type_ne_of_mem_children hi hcm
  have hQ := h.tree.v_children v c' hv hc'
  obtain ⟨e, he, rfl⟩ := h.type_Q_eq (h.ch_lt hi hc') hQ
  exact ⟨e, h.edgeIn_self he⟩

/-! ### Node-vertices -/

theorem nvOf_iff {i : ItemId} (hi : i < items.size) (nv : Nat) :
    t.NvOf (idx i) nv ↔ ∃ k, k < (items.nvList g i).length ∧ nv = (t.nvRange (idx i)).1 + k := by
  have := (h.node i hi).nv_range
  unfold SpqrTree.NvOf
  constructor
  · rintro ⟨h1, h2⟩; exact ⟨nv - (t.nvRange (idx i)).1, by omega, by omega⟩
  · rintro ⟨k, hk, rfl⟩; omega

theorem vertIndex_eq {v : Nat} (hv : v < g.nv) : t.vertIndex[v]! = some (idx (vertItem v)) := by
  have := (h.node _ (h.vertItem_lt hv)).vert_index (h.type_vertItem hv)
  rwa [show vertItem v - 1 = v by riomega] at this

theorem nodeVerts_get' {i : ItemId} (hi : i < items.size) {k : Nat}
    (hk : k < (items.nvList g i).length) :
    t.nodeVerts[(t.nvRange (idx i)).1 + k]? =
      some ⟨idx i, idx (vertItem (items.nvList g i)[k])⟩ := by
  rw [(h.nodeVerts_get hi hk).2, h.vertIndex_eq (h.nvList_lt hi (List.getElem_mem hk)), Option.getD_some]

theorem SPR_vs {i : ItemId} (hi : i < items.size) (hT : items.type i = .S ∨ items.type i = .P ∨ items.type i = .R) :
    ∃ u v, items.vs i = (some u, some v) := by
  rcases hT with hT | hT | hT
  · obtain ⟨u, v, -, hvs, -⟩ := h.nvList_S hi hT; exact ⟨u, v, hvs⟩
  · obtain ⟨u, v, hvs, -⟩ := h.nvList_P hi hT; exact ⟨u, v, hvs⟩
  · obtain ⟨u, v, hvs, -⟩ := h.nvList_R hi hT; exact ⟨u, v, hvs⟩

theorem mem_mid_iff {i : ItemId} (hi : i < items.size) {x : Nat} (hx : x < g.nv) :
    x ∈ ((items.ch i).filter (· < 1 + g.nv)).map (· - 1) ↔ items.IsParent i (vertItem x) := by
  rw [List.mem_map]
  constructor
  · rintro ⟨c, hc, rfl⟩
    obtain ⟨hcm, hlt⟩ := List.mem_filter.1 hc
    have hV := (h.type_V_iff hi hcm).2 (by simpa using hlt)
    obtain ⟨w, -, rfl⟩ := h.vertItem_of_V hV
    rwa [show vertItem w - 1 = w by riomega]
  · intro hp
    exact ⟨vertItem x, List.mem_filter.2 ⟨hp, by rw [decide_eq_true_eq]; riomega⟩, by riomega⟩

/-- The vertices of an S/P/R node: its two endpoints and its V children. -/
theorem mem_nvList_iff {i : ItemId} (hi : i < items.size)
    (hT : items.type i = .S ∨ items.type i = .P ∨ items.type i = .R) {w : Nat} (hw : w < g.nv) :
    w ∈ items.nvList g i ↔
      (items.vs i).1 = some w ∨ (items.vs i).2 = some w ∨ items.IsParent i (vertItem w) := by
  obtain ⟨u, v, hvs⟩ := h.SPR_vs hi hT
  rw [nvList_eq, hvs, ← h.mem_mid_iff hi hw]
  simp only [Option.toList_some, List.mem_append, List.mem_singleton, Option.some.injEq]
  constructor
  · rintro ((rfl | hm) | rfl)
    · exact Or.inl rfl
    · exact Or.inr (Or.inr hm)
    · exact Or.inr (Or.inl rfl)
  · rintro (h1 | h2 | hm)
    · exact Or.inl (Or.inl h1.symm)
    · exact Or.inr h2.symm
    · exact Or.inl (Or.inr hm)

theorem nvOrig_of_mem {i : ItemId} (hi : i < items.size) {w : Nat} (hw : w ∈ items.nvList g i) :
    ∃ nv, t.NvOf (idx i) nv ∧ t.nvOrig nv = some w := by
  obtain ⟨k, hk, hkw⟩ := List.getElem_of_mem hw
  exact ⟨(t.nvRange (idx i)).1 + k, (h.nvOf_iff hi _).2 ⟨k, hk, rfl⟩, by rw [h.nvOrig_get hi hk, hkw]⟩

/-- The middle of an S/P/R node's `nvList` lists exactly its V children. -/
theorem nvList_mid_iff {i : ItemId} (hi : i < items.size)
    (hT : items.type i = .S ∨ items.type i = .P ∨ items.type i = .R) {k : Nat}
    (hk : k < (items.nvList g i).length) :
    items.IsParent i (vertItem (items.nvList g i)[k]) ↔ k ≠ 0 ∧ k ≠ (items.nvList g i).length - 1 := by
  obtain ⟨u, v, hvs⟩ := h.SPR_vs hi hT
  have hnd := h.nv_nodup hi
  set mid := ((items.ch i).filter (· < 1 + g.nv)).map (· - 1) with hmid
  have hnl : items.nvList g i = u :: mid ++ [v] := by rw [nvList_eq, hvs]; rfl
  have hmem : ∀ x, x < g.nv → (x ∈ mid ↔ items.IsParent i (vertItem x)) := fun x hx => h.mem_mid_iff hi hx
  have hlen : (items.nvList g i).length = mid.length + 2 := by rw [hnl]; simp
  have hget : (items.nvList g i)[k] = (u :: (mid ++ [v]))[k]'(by rw [← List.cons_append, ← hnl]; exact hk) :=
    List.getElem_of_eq hnl hk
  have hx := h.nvList_lt hi (List.getElem_mem hk)
  rw [hnl, List.cons_append] at hnd
  have hnd' := hnd
  rw [List.nodup_cons, List.mem_append] at hnd'
  obtain ⟨hu, hnd2⟩ := hnd'
  rw [List.nodup_append] at hnd2
  obtain ⟨-, -, hdisj⟩ := hnd2
  rw [← hmem _ hx, hget]
  simp only [hlen]
  constructor
  · intro hmx
    constructor
    · rintro rfl; exact hu (Or.inl hmx)
    · intro hlast
      have : (u :: (mid ++ [v]))[k]'(by simp; omega) = v := by
        subst hlast
        simp [List.getElem_cons_succ, List.getElem_append_right]
      rw [this] at hmx
      exact hdisj _ hmx v (List.mem_singleton_self v) rfl
  · rintro ⟨h0, hl⟩
    obtain ⟨k, rfl⟩ : ∃ k', k = k' + 1 := ⟨k - 1, by omega⟩
    rw [List.getElem_cons_succ, List.getElem_append_left (by rw [hlen] at hk; omega)]
    exact List.getElem_mem _

theorem capEnd_iff {i : ItemId} (hi : i < items.size)
    (hT : items.type i = .S ∨ items.type i = .P ∨ items.type i = .R) (nv : Nat) :
    t.CapEnd (idx i) nv ↔
      nv = (t.nvRange (idx i)).1 ∨ nv = (t.nvRange (idx i)).1 + ((items.nvList g i).length - 1) := by
  have hn := isNode_of_SPR hT
  have hc : items.hasCap i = true := by
    rcases hT with hT | hT | hT <;> simp [Items.hasCap, hT, NodeType.isNode]
  have hcap : t.capNe (idx i) = some (t.neRange (idx i)).1 := by
    unfold SpqrTree.capNe; rw [h.hasCap_eq hi, ite_eq_left hc]
  have h1 : 1 ≤ items.nEdges g i := by
    rw [h.nEdges_node hi hn]; simp [Items.capCount, hc]
  obtain ⟨-, hget⟩ := h.nodeEdges_get hi (k := 0) h1
  rw [Nat.add_zero] at hget
  have hnvs := h.edge0_nvs hi hn h1
  unfold SpqrTree.CapEnd
  constructor
  · rintro ⟨ne, d, hne, hd, hor⟩
    rw [hcap] at hne; cases hne
    rw [hget] at hd; cases hd
    rw [hnvs] at hor; exact hor.imp Eq.symm Eq.symm
  · intro hor
    exact ⟨_, _, hcap, hget, by rw [hnvs]; exact hor.imp Eq.symm Eq.symm⟩

theorem pieceSep_v_child_nv : ∀ i c, i < t.size → (t.type i = .S ∨ t.type i = .P ∨ t.type i = .R) →
    c ∈ t.children i → t.type c = .V → ∃ nv d, t.NvOf i nv ∧ t.nodeVerts[nv]? = some d ∧ d.vert = c := by
  intro n c hn hT hcm hcV
  obtain ⟨i, hi, rfl⟩ := h.idx_surj hn
  rw [h.type_eq hi] at hT
  obtain ⟨c', hc', rfl⟩ := h.type_ne_of_mem_children hi hcm
  rw [h.type_eq (h.ch_lt hi hc')] at hcV
  obtain ⟨x, hx, rfl⟩ := h.vertItem_of_V hcV
  obtain ⟨u, v, hvs⟩ := h.SPR_vs hi hT
  have hxm : x ∈ items.nvList g i := by
    rw [nvList_eq, hvs]
    simp only [Option.toList_some, List.mem_append, List.mem_map, List.mem_filter, List.mem_singleton]
    exact Or.inl (Or.inr ⟨vertItem x, ⟨hc', by rw [decide_eq_true_eq]; riomega⟩, by riomega⟩)
  obtain ⟨k, hk, hkx⟩ := List.getElem_of_mem hxm
  refine ⟨(t.nvRange (idx i)).1 + k, _, (h.nvOf_iff hi _).2 ⟨k, hk, rfl⟩, h.nodeVerts_get' hi hk, ?_⟩
  simp [hkx]

theorem pieceSep_nv_vert_inj : ∀ i nv nv' d d', i < t.size →
    (t.type i = .S ∨ t.type i = .P ∨ t.type i = .R) →
    t.NvOf i nv → t.NvOf i nv' → t.nodeVerts[nv]? = some d → t.nodeVerts[nv']? = some d' →
    d.vert = d'.vert → nv = nv' := by
  intro n nv nv' d d' hn hT h1 h2 hd hd' hvv
  obtain ⟨i, hi, rfl⟩ := h.idx_surj hn
  obtain ⟨k, hk, rfl⟩ := (h.nvOf_iff hi _).1 h1
  obtain ⟨k', hk', rfl⟩ := (h.nvOf_iff hi _).1 h2
  rw [h.nodeVerts_get' hi hk] at hd
  rw [h.nodeVerts_get' hi hk'] at hd'
  cases hd; cases hd'
  simp only at hvv
  have hv := h.nvList_lt hi (List.getElem_mem hk)
  have hv' := h.nvList_lt hi (List.getElem_mem hk')
  have := h.idx_inj (h.vertItem_lt hv) (h.vertItem_lt hv') hvv
  have hkk : k = k' := (h.nv_nodup hi).getElem_inj_iff.1 (by riomega)
  rw [hkk]

theorem pieceSep_nv_orig_inj : ∀ i nv nv', i < t.size →
    (t.type i = .S ∨ t.type i = .P ∨ t.type i = .R) →
    t.NvOf i nv → t.NvOf i nv' → t.nvOrig nv = t.nvOrig nv' → nv = nv' := by
  intro n nv nv' hn hT h1 h2 ho
  obtain ⟨i, hi, rfl⟩ := h.idx_surj hn
  obtain ⟨k, hk, rfl⟩ := (h.nvOf_iff hi _).1 h1
  obtain ⟨k', hk', rfl⟩ := (h.nvOf_iff hi _).1 h2
  rw [h.nvOrig_get hi hk, h.nvOrig_get hi hk', Option.some.injEq] at ho
  have hkk : k = k' := (h.nv_nodup hi).getElem_inj_iff.1 ho
  rw [hkk]

theorem pieceSep_nv_child : ∀ i nv d, i < t.size → (t.type i = .S ∨ t.type i = .P ∨ t.type i = .R) →
    t.NvOf i nv → t.nodeVerts[nv]? = some d → (t.parent d.vert = some i ↔ ¬ t.CapEnd i nv) := by
  intro n nv d hn hT hnv hd
  obtain ⟨i, hi, rfl⟩ := h.idx_surj hn
  rw [h.type_eq hi] at hT
  obtain ⟨k, hk, rfl⟩ := (h.nvOf_iff hi _).1 hnv
  rw [h.nodeVerts_get' hi hk] at hd
  cases hd
  have hv := h.nvList_lt hi (List.getElem_mem hk)
  rw [h.capEnd_iff hi hT, h.parent_eq_iff hi (h.vertItem_lt hv), h.nvList_mid_iff hi hT hk]
  omega

/-! ### Virtual edges, exactly oriented -/

/-- Virtual edge `j` of node `a` (`capCount a + j` in its node-edge range) joins the node-vertices
holding the cap endpoints `(x, y)` of the `j`-th non-V ordered child, in that order; so its
`neOrig` equals the child's cap `neOrig` exactly (not just up to swapping). -/
theorem virt_edge (hqc : items.QChildVs g) (hpc : items.PChildVs) {a : ItemId} (ha : a < items.size)
    (hn : (items.type a).isNode = true) {pos : Nat → Nat} (hl : RelabelLayout g items t idx a pos)
    {j : Nat} (hj : j < ((items.ordered g a (t.nvRange (idx a)).1 pos).filter (· ≥ 1 + g.nv)).length) :
    ∃ x y ka kb, ka < (items.nvList g a).length ∧ kb < (items.nvList g a).length ∧
      (items.nvList g a)[ka]? = some x ∧ (items.nvList g a)[kb]? = some y ∧
      (items.vs ((items.ordered g a (t.nvRange (idx a)).1 pos).filter (· ≥ 1 + g.nv))[j]).1 = some x ∧
      (items.vs ((items.ordered g a (t.nvRange (idx a)).1 pos).filter (· ≥ 1 + g.nv))[j]).2.getD x = y ∧
      (t.nodeEdges[(t.neRange (idx a)).1 + (items.capCount a + j)]!).nvs =
        ((t.nvRange (idx a)).1 + ka, (t.nvRange (idx a)).1 + kb) ∧
      t.neOrig ((t.neRange (idx a)).1 + (items.capCount a + j)) = some (x, y) ∧
      t.neOrig (t.neRange (idx ((items.ordered g a (t.nvRange (idx a)).1 pos).filter
        (· ≥ 1 + g.nv))[j])).1 = some (x, y) := by
  set F := (items.ordered g a (t.nvRange (idx a)).1 pos).filter (· ≥ 1 + g.nv) with hF
  have hperm : F.Perm ((items.ch a).filter (· ≥ 1 + g.nv)) := ordered_filter_perm a _ pos _
  have hFlen : F.length = (items.virtualEdges a).length := by
    rw [hperm.length_eq, h.virt_length ha, List.countP_eq_length_filter]
  have hcF : F[j] ∈ F := List.getElem_mem hj
  have hcm : F[j] ∈ items.ch a := List.mem_of_mem_filter (hperm.mem_iff.1 hcF)
  have hcge : F[j] ≥ 1 + g.nv := by simpa using List.of_mem_filter (hperm.mem_iff.1 hcF)
  have hcV : items.type F[j] ≠ .V := fun hV =>
    absurd ((h.type_V_iff ha hcm).1 hV) (not_lt.2 hcge)
  obtain ⟨x, y, hx, hy, hcq, hO⟩ := h.child_cap_orig ha hn hcm hcV
  have hxm : x ∈ items.nvList g a := h.mem_nvList_of_child hn hcm hcV (Or.inl hx)
  have hcs := h.ch_lt ha hcm
  have hym : y ∈ items.nvList g a := by
    by_cases hO' : items.type F[j] = .O
    · obtain ⟨v, hv, -⟩ := h.nvList_O hcs hO'
      rw [hv] at hy; simp at hy; subst hy; exact hxm
    · exact h.mem_nvList_of_child hn hcm hcV (Or.inr (hO hO'))
  have hnE := h.nEdges_node ha hn
  have hnv := (h.node a ha).nv_range
  have hne := (h.node a ha).ne_range
  have hk : items.capCount a + j < items.nEdges g a := by omega
  have hlay := hl.edge_nvs _ hk
  rw [nodeLayout, hnv, hne, Items.edgeChildren, ← hF] at hlay
  have hnotO : items.type a ≠ .Q → items.type F[j] ≠ .O :=
    fun hQ hO' => hQ (h.wf.shapes.o_parent a _ hcm hO').1
  set nvSt := (t.nvRange (idx a)).1 with hnvSt
  set neSt := (t.neRange (idx a)).1 with hneSt
  set len := (items.nvList g a).length with hlen0
  cases hty : items.type a
  · rw [hty] at hn; cases hn
  · rw [hty] at hn; cases hn
  · -- Q: a block root with one non-V child
    obtain ⟨e, he, rfl⟩ := h.type_Q_eq ha hty
    obtain ⟨u, v, hu, -, hcase, -, hne1⟩ := h.q_pair he
    rcases hcase with ⟨h0, -, -⟩ | ⟨hne0, hvs, c, hc⟩
    · rw [h0] at hcm; cases hcm
    have hcap0 : items.capCount (edgeItem g e) = 0 := by
      simp [Items.capCount, Items.hasCap, hty, hne0]
    have hj0 : j = 0 := by omega
    subst hj0
    simp only [hcap0, Nat.zero_add, Nat.add_zero] at hlay hk ⊢
    have hpos := nvList_pos (g := g) hu
    have hnvs := h.edge0_nvs ha hn (by omega)
    have hp := h.neOrig_eq ha (k := 0) (a := 0) (b := len - 1) (by omega) (by omega) (by omega)
      (by rw [Nat.add_zero]; exact hnvs)
    rw [Nat.add_zero] at hp
    rcases hc with ⟨hch, hnl⟩ | ⟨hch, hnl, -⟩
    · have hvc := (hqc e he).1 _ _ hch
      rw [hvs] at hvc
      have hF0 : F[0] = c := by
        rw [hch] at hcm
        simp only [List.mem_cons, List.not_mem_nil, or_false] at hcm
        rcases hcm with hcm | hcm
        · exact hcm
        · exfalso
          have hv' := h.nvList_lt ha (v := v) (by rw [hnl]; simp)
          rw [hcm] at hcge; riomega
      have hx' : x = u := by have := hx; rw [hF0, hvc] at this; simpa using this.symm
      have hy' : y = v := by have := hy; rw [hF0, hvc] at this; simpa using this.symm
      refine ⟨x, y, 0, len - 1, by omega, by omega, ?_, ?_, hx, hy, by simpa using hnvs, ?_, hcq⟩
      · simp [hnl, hx']
      · simp [hlen0, hnl, hy']
      · rw [hp]; simp [hlen0, hnl, hx', hy']
    · have hvc := (hqc e he).2 _ hch
      rw [hvs] at hvc
      have hF0 : F[0] = c := by
        rw [hch] at hcm
        simpa using hcm
      have hx' : x = u := by have := hx; rw [hF0, hvc] at this; simpa using this.symm
      have hy' : y = x := by have := hy; rw [hF0, hvc] at this; simpa using this.symm
      refine ⟨x, y, 0, len - 1, by omega, by omega, ?_, ?_, hx, hy, by simpa using hnvs, ?_, hcq⟩
      · simp [hnl, hx']
      · simp [hlen0, hnl, hy', hx']
      · rw [hp]; simp [hlen0, hnl, hx', hy']
  · rw [h.shapes.i_o_leaf a ha (Or.inl hty)] at hcm; cases hcm
  · rw [h.shapes.i_o_leaf a ha (Or.inr hty)] at hcm; cases hcm
  · -- S
    rw [hty] at hlay
    obtain ⟨u, v, xs, hvs, hnl, hxs, hvirt⟩ := h.nvList_S ha hty
    have hcap1 : items.capCount a = 1 := by
      simp [Items.capCount, Items.hasCap, hty, NodeType.isNode]
    have hFeq : F = (items.ch a).filter fun c => items.type c ≠ .V := by
      rw [hF, Items.ordered, ite_eq_left (by rw [hty]; decide), h.filter_ge_eq ha]
    have hvirt' : items.virtualEdges a =
        F.map fun c => ((items.vs c).1.getD 0, (items.vs c).2.getD 0) := by
      rw [Items.virtualEdges, hFeq]
    have hlen : len = xs.length + 2 := by rw [hlen0, hnl]; simp
    have hvl : (items.virtualEdges a).length = xs.length + 1 := by
      rw [hvirt, List.length_zip]; simp
    rw [hcap1] at hlay hk ⊢
    rw [hvl, hcap1] at hnE
    have hs := LayoutFacts.s_edge (idx a) nvSt (nvSt + len) neSt (neSt + items.nEdges g a)
      (F.map fun c => (pos ((items.vs c).1.getD 0), pos ((items.vs c).2.getD 0)))
      (by omega) (by omega) (1 + j) (by omega)
    rw [ite_eq_right (by omega)] at hs
    rw [hs] at hlay
    have hnvs : (t.nodeEdges[neSt + (1 + j)]!).nvs = (nvSt + j, nvSt + (j + 1)) := by
      rw [hlay]; exact Prod.ext (by dsimp only; omega) (by dsimp only; omega)
    have hp := h.neOrig_eq ha (k := 1 + j) (a := j) (b := j + 1) hk (by omega) (by omega) hnvs
    have hy' := hO (hnotO (by rw [hty]; decide))
    have hvj : (items.virtualEdges a)[j]'(by omega) = (x, y) := by
      rw [List.getElem_of_eq hvirt', List.getElem_map, hx, hy']; rfl
    rw [List.getElem_of_eq hvirt, List.getElem_zip] at hvj
    have hj1 : j < (u :: xs).length := by simp; omega
    have hj2 : j < (xs ++ [v]).length := by simp; omega
    have e1 : (items.nvList g a)[j]? = some ((u :: xs)[j]'hj1) := by
      rw [hnl, List.getElem?_append_left hj1, List.getElem?_eq_getElem]
    have e2 : (items.nvList g a)[j + 1]? = some ((xs ++ [v])[j]'hj2) := by
      rw [hnl, List.cons_append, List.getElem?_cons_succ, List.getElem?_eq_getElem]
    simp only [Prod.mk.injEq] at hvj
    rw [hvj.1] at e1
    rw [hvj.2] at e2
    refine ⟨x, y, j, j + 1, by omega, by omega, e1, e2, hx, hy, hnvs, ?_, hcq⟩
    rw [hp, (List.getElem_eq_iff _).2 e1, (List.getElem_eq_iff _).2 e2]
  · -- P
    rw [hty] at hlay
    obtain ⟨u, v, hvs, hnl⟩ := h.nvList_P ha hty
    have hcap1 : items.capCount a = 1 := by
      simp [Items.capCount, Items.hasCap, hty, NodeType.isNode]
    have hlen : len = 2 := by rw [hlen0, hnl]; rfl
    rw [hcap1] at hlay hk ⊢
    rw [hcap1] at hnE
    rw [hlen] at hlay
    have hs := LayoutFacts.p_edge (idx a) nvSt neSt (neSt + items.nEdges g a)
      (F.map fun c => (pos ((items.vs c).1.getD 0), pos ((items.vs c).2.getD 0))) (1 + j) (by omega)
    rw [hs] at hlay
    have hnvs : (t.nodeEdges[neSt + (1 + j)]!).nvs = (nvSt + 0, nvSt + 1) := by rw [hlay]; rfl
    have hp := h.neOrig_eq ha (k := 1 + j) (a := 0) (b := 1) hk (by omega) (by omega) hnvs
    have hy' := hO (hnotO (by rw [hty]; decide))
    have hvc := hpc a _ ha hty hcm hcV
    rw [hvs] at hvc
    have hx' : x = u := by have := hx; rw [hvc] at this; simpa using this.symm
    have hy'' : y = v := by have := hy'; rw [hvc] at this; simpa using this.symm
    refine ⟨x, y, 0, 1, by omega, by omega, by simp [hnl, hx'], by simp [hnl, hy''], hx, hy, hnvs, ?_, hcq⟩
    rw [hp]
    simp [hnl, hx', hy'']
  · -- R
    rw [hty] at hlay
    obtain ⟨u, v, hvs, hlen4⟩ := h.nvList_R ha hty
    have hcap1 : items.capCount a = 1 := by
      simp [Items.capCount, Items.hasCap, hty, NodeType.isNode]
    rw [hcap1] at hlay hk ⊢
    rw [hcap1] at hnE
    have hy' := hO (hnotO (by rw [hty]; decide))
    have hpos := hl.pos_ok hty
    obtain ⟨hxle, hxp⟩ := hpos x hxm
    obtain ⟨hyle, hyp⟩ := hpos y hym
    have hxa := (List.getElem?_eq_some_iff.1 hxp).1
    have hyb := (List.getElem?_eq_some_iff.1 hyp).1
    set E := F.map fun c => (pos ((items.vs c).1.getD 0), pos ((items.vs c).2.getD 0)) with hE
    have hElen : E.length = F.length := List.length_map _
    have hs := LayoutFacts.r_edge (idx a) nvSt (nvSt + len) neSt (neSt + items.nEdges g a) E
      (by omega) (by omega) (1 + j) (by omega)
    rw [ite_eq_right (by omega), show 1 + j - 1 = j by omega, getElem!_pos E j (by omega)] at hs
    have hEj : E[j]'(by omega) = (pos x, pos y) := by
      simp only [hE, List.getElem_map, hx, hy', Option.getD_some]
    rw [hEj] at hs
    rw [hs] at hlay
    have hnvs : (t.nodeEdges[neSt + (1 + j)]!).nvs = (nvSt + (pos x - nvSt), nvSt + (pos y - nvSt)) := by
      rw [hlay]; exact Prod.ext (by dsimp only; omega) (by dsimp only; omega)
    have hp := h.neOrig_eq ha (k := 1 + j) (a := pos x - nvSt) (b := pos y - nvSt) hk hxa hyb hnvs
    refine ⟨x, y, pos x - nvSt, pos y - nvSt, hxa, hyb, hxp, hyp, hx, hy, hnvs, ?_, hcq⟩
    rw [hp, (List.getElem_eq_iff _).2 hxp, (List.getElem_eq_iff _).2 hyp]

theorem pieceSep_twin_orient (hqc : items.QChildVs g) (hpc : items.PChildVs) :
    ∀ ne tw, t.twin ne = some tw → t.neOrig ne = t.neOrig tw := by
  intro ne ne' ht
  obtain ⟨a, ha, h1, h2⟩ := h.ne_mem_range (twin_some_lt ht)
  obtain ⟨pos, hl⟩ := (h.node a ha).layout
  have hner := (h.node a ha).ne_range
  rw [hner] at h2
  by_cases hn : (items.type a).isNode = true
  swap
  · exfalso
    have := nEdges_not_node (g := g) (i := a) (by simpa using hn)
    omega
  have hnE := h.nEdges_node ha hn
  have hFlen : ((items.ordered g a (t.nvRange (idx a)).1 pos).filter (· ≥ 1 + g.nv)).length =
      (items.virtualEdges a).length := by
    rw [(ordered_filter_perm a _ pos _).length_eq, h.virt_length ha, List.countP_eq_length_filter]
  by_cases hk : (t.neRange (idx a)).1 + items.capCount a ≤ ne
  · -- `ne` is a virtual edge of `a`
    have hj : ne - ((t.neRange (idx a)).1 + items.capCount a) <
        ((items.ordered g a (t.nvRange (idx a)).1 pos).filter (· ≥ 1 + g.nv)).length := by omega
    obtain ⟨x, y, -, -, -, -, -, -, -, -, -, hp', hq'⟩ := h.virt_edge hqc hpc ha hn hl hj
    have hne : (t.neRange (idx a)).1 + (items.capCount a +
        (ne - ((t.neRange (idx a)).1 + items.capCount a))) = ne := by omega
    rw [hne] at hp'
    have htw := (hl.twin hn _ hj).1
    rw [Nat.add_assoc, hne, ht] at htw
    cases htw
    rw [hp', hq']
  · -- `ne` is the cap of `a`: twinned with a virtual edge of its parent
    have hc : items.hasCap a = true := by
      by_contra hc
      have : items.capCount a = 0 := by simp [Items.capCount, hc]
      omega
    have hne : ne = (t.neRange (idx a)).1 := by
      have : items.capCount a = 1 := by simp [Items.capCount, hc]
      omega
    subst hne
    have ha0 : a ≠ rootItem := by
      rintro rfl; rw [h.tree.root] at hn; cases hn
    obtain ⟨p0, hp0, -⟩ := h.tree.unique_parent a (Nat.pos_of_ne_zero ha0) ha
    have hp0s := parent_lt hp0
    by_cases hn0 : (items.type p0).isNode = true
    · obtain ⟨pos0, hl0⟩ := (h.node p0 hp0s).layout
      have hperm := Items.ordered_perm (items := items) (g := g) (i := p0) (nvSt := (t.nvRange (idx p0)).1) (pos := pos0)
      have haV : items.type a ≠ .V := by intro hV; rw [hV] at hn; cases hn
      have hage : ¬ a < 1 + g.nv := fun hlt => haV ((h.type_V_iff hp0s hp0).2 hlt)
      have hmem : a ∈ (items.ordered g p0 (t.nvRange (idx p0)).1 pos0).filter (· ≥ 1 + g.nv) :=
        List.mem_filter.2 ⟨hperm.mem_iff.2 hp0, by simpa using hage⟩
      obtain ⟨j, hj, hja⟩ := List.getElem_of_mem hmem
      obtain ⟨x, y, -, -, -, -, -, -, -, -, -, hp', hq'⟩ := h.virt_edge hqc hpc hp0s hn0 hl0 hj
      have htw := (hl0.twin hn0 j hj).2
      rw [hja] at htw hq'
      have htw := htw hc
      rw [ht] at htw
      cases htw
      rw [hq', Nat.add_assoc, hp']
    · exfalso
      have := (h.node p0 hp0s).child_cap_twin_none (by simpa using hn0) a hp0 hc
      rw [this] at ht; cases ht

/-! ### Tree order -/

theorem below_lt {a b : ItemId} (ha : a < items.size) (hab : items.Below a b) : b < items.size := by
  induction hab with
  | refl => exact ha
  | tail _ hbc _ => exact h.tree.ch_lt _ _ hbc

theorem below_antisymm {a b : ItemId} (hab : items.Below a b) (hba : items.Below b a) : a = b := by
  rcases Relation.reflTransGen_iff_eq_or_transGen.1 hab with rfl | h1
  · rfl
  rcases Relation.reflTransGen_iff_eq_or_transGen.1 hba with rfl | h2
  · rfl
  exact absurd (h1.trans h2) (h.tree.acyclic a)

theorem not_below_parent {p c : ItemId} (hpc : items.IsParent p c) : ¬ items.Below c p :=
  fun hcp => h.tree.acyclic p (Relation.TransGen.head' hpc hcp)

omit h in
/-- An edge below a V item lies below one of its (block-root Q) children. -/
theorem v_unwind {v e : Nat} (hv : v < g.nv) (hb : items.EdgeBelow g (vertItem v) e) :
    ∃ c, items.IsParent (vertItem v) c ∧ items.EdgeBelow g c e := by
  rcases Relation.ReflTransGen.cases_head hb with heq | ⟨c, hc, hce⟩
  · exfalso; riomega
  · exact ⟨c, hc, hce⟩

omit h in
theorem root_unwind {e : Nat} (hb : items.EdgeBelow g rootItem e) :
    ∃ c, items.IsParent rootItem c ∧ items.EdgeBelow g c e := by
  rcases Relation.ReflTransGen.cases_head hb with heq | ⟨c, hc, hce⟩
  · exfalso; riomega
  · exact ⟨c, hc, hce⟩

/-- A child of a V item is a block-root Q with `vs = (some v, none)`. -/
theorem q_block_vs (hr : Items.Ranges g items σ) (hqu : items.QUpper g) {v c : Nat} (hv : v < g.nv)
    (hc : items.IsParent (vertItem v) c) : items.type c = .Q ∧ items.vs c = (some v, none) := by
  have hQ := h.tree.v_children v c hv hc
  have hcs := h.tree.ch_lt _ _ hc
  obtain ⟨e, he, rfl⟩ := h.type_Q_eq hcs hQ
  obtain ⟨u, -, hvs, -⟩ := hr.q_root e he (hr.q_under_v v _ hv hc)
  have hu := hqu v _ hv hc
  rw [hvs] at hu
  simp only [Option.some.injEq] at hu
  subst hu
  exact ⟨hQ, hvs⟩

theorem not_below_of_sibling {p a b : ItemId} (ha : items.IsParent p a) (hb : items.IsParent p b)
    (hne : a ≠ b) {e : Nat} (he : items.EdgeBelow g a e) : ¬ items.EdgeBelow g b e := by
  intro he'
  have hd := h.tree.desc_disjoint ha hb hne
  have hes : edgeItem g e < items.size := h.below_lt (h.tree.ch_lt _ _ ha) he
  exact Finset.disjoint_left.1 hd (mem_desc.2 ⟨hes, he⟩) (mem_desc.2 ⟨hes, he'⟩)

theorem pieceSep_root_disjoint (hrs : items.RootSep g) :
    ∀ a ∈ t.children 0, ∀ b ∈ t.children 0, a ≠ b →
      ∀ v, t.Touches g a v → t.Touches g b v → False := by
  intro a ha b hb hne v hta htb
  have hrs' : rootItem < items.size := by have := h.tree.size; riomega
  rw [← h.ridx.root] at ha hb
  obtain ⟨a', ha', rfl⟩ := h.type_ne_of_mem_children hrs' ha
  obtain ⟨b', hb', rfl⟩ := h.type_ne_of_mem_children hrs' hb
  obtain ⟨e, he, hinc, hbe⟩ := (h.touches_iff (h.tree.ch_lt _ _ ha') v).1 hta
  obtain ⟨e', he', hinc', hbe'⟩ := (h.touches_iff (h.tree.ch_lt _ _ hb') v).1 htb
  exact hrs a' b' ha' hb' (fun heq => hne (by rw [heq])) v e e' he he' hinc hinc' hbe hbe'

/-! ### Attachment -/

omit h in
theorem inc_lt (hg : g.WF) {e w : Nat} (he : e < g.ne) (hw : g.Inc e w) : w < g.nv := by
  unfold Graph.Inc at hw
  rw [getElem!_pos g.edges e he] at hw
  have := hg _ (Array.getElem_mem he)
  rcases hw with rfl | rfl
  · exact this.1
  · exact this.2

theorem sep_vs (hg : g.WF) {i : ItemId} (hi : i < items.size) (hFV : items.type i ∉ [NodeType.F, .V])
    {w e e' : Nat} (he : e < g.ne) (he' : e' < g.ne) (hw : g.Inc e w) (hw' : g.Inc e' w)
    (hb : items.EdgeBelow g i e) (hnb : ¬ items.EdgeBelow g i e') :
    (items.vs i).1 = some w ∨ (items.vs i).2 = some w :=
  h.endpoints.separation i hi hFV w e e' (inc_lt hg he hw) he he' hw hw' hb hnb

/-- A vertex touched below a V item `v` and by an edge not below `vertItem v` is `v` itself. -/
theorem v_touch_eq (hg : g.WF) (hr : Items.Ranges g items σ) (hqu : items.QUpper g)
    {v w e e' : Nat} (hv : v < g.nv) (he : e < g.ne) (he' : e' < g.ne)
    (hb : items.EdgeBelow g (vertItem v) e) (hw : g.Inc e w) (hw' : g.Inc e' w)
    (hnb : ¬ items.EdgeBelow g (vertItem v) e') : w = v := by
  obtain ⟨c, hc, hce⟩ := v_unwind hv hb
  have hcs := h.tree.ch_lt _ _ hc
  obtain ⟨hQ, hvs⟩ := h.q_block_vs hr hqu hv hc
  have hnb' : ¬ items.EdgeBelow g c e' := fun hce' => hnb (Relation.ReflTransGen.head hc hce')
  have := h.sep_vs hg hcs (by rw [hQ]; decide) he he' hw hw' hce hnb'
  rw [hvs] at this
  simp only [Option.some.injEq, reduceCtorEq, or_false] at this
  exact this.symm

theorem child_hasCap {i c : ItemId} (hi : i < items.size)
    (hT : items.type i = .S ∨ items.type i = .P ∨ items.type i = .R) (hc : c ∈ items.ch i)
    (hcV : items.type c ≠ .V) : items.hasCap c = true := by
  have hcs := h.ch_lt hi hc
  refine hasCap_of_node ?_ hcV fun hQ => h.shapes.q_leaf_of_node i c hc ?_ hQ
  · intro hF; exact h.not_root_of_child hc (h.root_of_F hcs hF)
  · rcases hT with hT | hT | hT <;> rw [hT] <;> decide

/-- At a vertex `w` of an S/P/R node some edge lies outside any given child. -/
theorem node_vertex_outside (hr : Items.Ranges g items σ) {i : ItemId} (hi : i < items.size)
    (hT : items.type i = .S ∨ items.type i = .P ∨ items.type i = .R) {w : Nat} (hw : w < g.nv)
    (hwm : w ∈ items.nvList g i) {c : ItemId} (hc : items.IsParent i c) :
    ∃ e', e' < g.ne ∧ g.Inc e' w ∧ ¬ items.EdgeBelow g c e' := by
  rcases (h.mem_nvList_iff hi hT hw).1 hwm with h1 | h2 | hp
  · obtain ⟨e, e', -, he', -, hinc', -, hnb⟩ := hr.vs_att i hi (Or.inl (SPR_mem hT)) w (Or.inl h1)
    exact ⟨e', he', hinc', fun hb => hnb (Relation.ReflTransGen.head hc hb)⟩
  · obtain ⟨e, e', -, he', -, hinc', -, hnb⟩ := hr.vs_att i hi (Or.inl (SPR_mem hT)) w (Or.inr h2)
    exact ⟨e', he', hinc', fun hb => hnb (Relation.ReflTransGen.head hc hb)⟩
  · obtain ⟨-, -, hch⟩ := (h.endpoints.interior i w hi hw (SPR_mem hT)).1 hp
    have := hch c hc
    push_neg at this
    obtain ⟨e', he', hinc', hnb⟩ := this
    exact ⟨e', he', hinc', hnb⟩

theorem pieceSep_v_attach (hg : g.WF) (hr : Items.Ranges g items σ) (hqu : items.QUpper g) :
    ∀ i v, i < t.size → t.type i = .V → t.origId[i]! = some v →
      ∀ a ∈ t.children i, ∀ b ∈ t.children i, a ≠ b →
        ∀ w, t.Touches g a w → t.Touches g b w → w = v := by
  intro n v hn hV horig a ha b hb hne w hta htb
  obtain ⟨i, hi, rfl⟩ := h.idx_surj hn
  rw [h.type_eq hi] at hV
  obtain ⟨v', hv', rfl⟩ := h.vertItem_of_V hV
  rw [h.origId_vertItem hv'] at horig
  cases horig
  obtain ⟨a', ha', rfl⟩ := h.type_ne_of_mem_children hi ha
  obtain ⟨b', hb', rfl⟩ := h.type_ne_of_mem_children hi hb
  have hne' : a' ≠ b' := fun heq => hne (by rw [heq])
  obtain ⟨e, he, hinc, hbe⟩ := (h.touches_iff (h.tree.ch_lt _ _ ha') w).1 hta
  obtain ⟨e', he', hinc', hbe'⟩ := (h.touches_iff (h.tree.ch_lt _ _ hb') w).1 htb
  obtain ⟨hQ, hvs⟩ := h.q_block_vs hr hqu hv' ha'
  have hnb := h.not_below_of_sibling hb' ha' hne'.symm hbe'
  have := h.sep_vs hg (h.tree.ch_lt _ _ ha') (by rw [hQ]; decide) he he' hinc hinc' hbe hnb
  rw [hvs] at this
  simp only [Option.some.injEq, reduceCtorEq, or_false] at this
  exact this.symm

theorem pieceSep_q_root_attach (hg : g.WF) (hr : Items.Ranges g items σ) :
    ∀ i c w e, i < t.size → t.type i = .Q →
      t.children i = [c, w] → t.type w = .V → t.origId[i]! = some e →
        ∀ v e', t.Touches g i v → SpqrTree.Graph.Incident g v e' → ¬ t.EdgeIn i e' →
          v = (g.edges[e]!).1 ∨ v = (g.edges[e]!).2 := by
  intro n c w e hn hQ hch hwV horig v e' htv hinc' hnin
  obtain ⟨i, hi, rfl⟩ := h.idx_surj hn
  rw [h.type_eq hi] at hQ
  obtain ⟨e₀, he₀, rfl⟩ := h.type_Q_eq hi hQ
  rw [h.origId_edgeItem he₀] at horig
  cases horig
  rw [h.children_eq hi (h.q_not_R he₀)] at hch
  obtain ⟨c', w', hch', -, -⟩ := map_eq_pair hch
  obtain ⟨u, c0, hvs, hinc, -⟩ := hr.q_root e he₀ (by rw [hch']; simp)
  obtain ⟨e₁, he₁, hinc₁, hb₁⟩ := (h.touches_iff hi v).1 htv
  obtain ⟨he', hinc'⟩ := hinc'
  have hnb : ¬ items.EdgeBelow g (edgeItem g e) e' := fun hb => hnin ((h.edgeIn_iff hi he').2 hb)
  have := h.sep_vs hg hi (by rw [hQ]; decide) he₁ he' hinc₁ hinc' hb₁ hnb
  rw [hvs] at this
  simp only [Option.some.injEq, reduceCtorEq, or_false] at this
  subst this
  exact Or.imp Eq.symm Eq.symm hinc

theorem pieceSep_q_lower_attach (hr : Items.Ranges g items σ) (hg : g.WF) (hqu : items.QUpper g) :
    ∀ i c w e, i < t.size → t.type i = .Q → t.children i = [c, w] →
      t.origId[i]! = some e → ∀ v, (t.Touches g c v ∨ SpqrTree.Graph.Incident g v e) → t.Touches g w v →
        v = if t.edgeFlipped[e]! then (g.edges[e]!).1 else (g.edges[e]!).2 := by
  intro n c w e hn hQ hch horig v hcv hwv
  have hlow := h.pieceSep_q_lower hr n c w e hn hQ hch horig
  obtain ⟨i, hi, rfl⟩ := h.idx_surj hn
  rw [h.type_eq hi] at hQ
  obtain ⟨e₀, he₀, rfl⟩ := h.type_Q_eq hi hQ
  rw [h.origId_edgeItem he₀] at horig
  cases horig
  rw [h.children_eq hi (h.q_not_R he₀)] at hch
  obtain ⟨c', w', hch', rfl, rfl⟩ := map_eq_pair hch
  obtain ⟨u, x, rfl, hx, -, -, -, -, -⟩ := h.q_two_children hr he₀ hch'
  rw [h.origId_vertItem hx, Option.some.injEq] at hlow
  rw [← hlow]
  have hcm : c' ∈ items.ch (edgeItem g e) := by rw [hch']; simp
  have hxm : vertItem x ∈ items.ch (edgeItem g e) := by rw [hch']; simp
  have hcx : c' ≠ vertItem x := by
    intro heq
    rw [heq] at hch'
    exact absurd (h.tree.ch_nodup (edgeItem g e)) (by rw [hch']; simp)
  obtain ⟨e₂, he₂, hinc₂, hb₂⟩ := (h.touches_iff (h.vertItem_lt hx) v).1 hwv
  rcases hcv with hcv | ⟨he₁, hinc₁⟩
  · obtain ⟨e₁, he₁, hinc₁, hb₁⟩ := (h.touches_iff (h.tree.ch_lt _ _ hcm) v).1 hcv
    exact h.v_touch_eq hg hr hqu hx he₂ he₁ hb₂ hinc₂ hinc₁ (h.not_below_of_sibling hcm hxm hcx hb₁)
  · exact h.v_touch_eq hg hr hqu hx he₂ he₁ hb₂ hinc₂ hinc₁ (h.not_below_parent hxm)

theorem pieceSep_node_attach (hg : g.WF) (hr : Items.Ranges g items σ) (hqu : items.QUpper g) :
    ∀ i a b w, i < t.size → (t.type i = .S ∨ t.type i = .P ∨ t.type i = .R) →
      a ∈ t.children i → b ∈ t.children i → a ≠ b → t.Touches g a w → t.Touches g b w →
      ∃ nv, t.NvOf i nv ∧ t.nvOrig nv = some w := by
  intro n a b w hn hT ha hb hne hta htb
  obtain ⟨i, hi, rfl⟩ := h.idx_surj hn
  rw [h.type_eq hi] at hT
  obtain ⟨a', ha', rfl⟩ := h.type_ne_of_mem_children hi ha
  obtain ⟨b', hb', rfl⟩ := h.type_ne_of_mem_children hi hb
  have hne' : a' ≠ b' := fun heq => hne (by rw [heq])
  have hbs := h.tree.ch_lt _ _ hb'
  obtain ⟨e, he, hinc, hbe⟩ := (h.touches_iff (h.tree.ch_lt _ _ ha') w).1 hta
  obtain ⟨e', he', hinc', hbe'⟩ := (h.touches_iff hbs w).1 htb
  have hnb := h.not_below_of_sibling ha' hb' hne' hbe
  have hw := inc_lt hg he hinc
  refine h.nvOrig_of_mem hi ((h.mem_nvList_iff hi hT hw).2 ?_)
  by_cases hbV : items.type b' = .V
  · obtain ⟨x, hx, rfl⟩ := h.vertItem_of_V hbV
    have := h.v_touch_eq hg hr hqu hx he' he hbe' hinc' hinc hnb
    subst this
    exact Or.inr (Or.inr hb')
  · have hbF : items.type b' ≠ .F := fun hF => h.not_root_of_child hb' (h.root_of_F hbs hF)
    have hFV : items.type b' ∉ [NodeType.F, .V] := by simp [hbF, hbV]
    have hvb := h.sep_vs hg hbs hFV he' he hinc' hinc hbe' hnb
    exact h.endpoints.child_vs_in_parent i b' hb'
      (by rcases hT with hT | hT | hT <;> simp [hT]) hbV w hvb

theorem nodeEdge_node {i : ItemId} (hi : i < items.size) {k : Nat} (hk : k < items.nEdges g i) :
    t.nodeOfNe ((t.neRange (idx i)).1 + k) = some (idx i) := by
  obtain ⟨-, -, hnode, -⟩ := RelabelAll.nodeEdge_at ⟨h.wf, h.ridx, h.node⟩ hi hk
  unfold SpqrTree.nodeOfNe
  rw [(h.nodeEdges_get hi hk).2]
  exact congrArg some hnode

theorem pieceSep_node_touch (hg : g.WF) (hr : Items.Ranges g items σ) (hqu : items.QUpper g)
    (hqc : items.QChildVs g) (hpc : items.PChildVs) :
    ∀ i c w nv, i < t.size → (t.type i = .S ∨ t.type i = .P ∨ t.type i = .R) →
      c ∈ t.children i → t.Touches g c w → t.NvOf i nv → t.nvOrig nv = some w →
      t.NvInc i c nv := by
  intro n c w nv hn hT hc htc hnv horig
  obtain ⟨i, hi, rfl⟩ := h.idx_surj hn
  rw [h.type_eq hi] at hT
  obtain ⟨c', hc', rfl⟩ := h.type_ne_of_mem_children hi hc
  have hcs := h.tree.ch_lt _ _ hc'
  obtain ⟨e, he, hinc, hbe⟩ := (h.touches_iff hcs w).1 htc
  obtain ⟨k, hk, rfl⟩ := (h.nvOf_iff hi _).1 hnv
  rw [h.nvOrig_get hi hk, Option.some.injEq] at horig
  have hw := inc_lt hg he hinc
  have hwm : w ∈ items.nvList g i := horig ▸ List.getElem_mem hk
  obtain ⟨e', he', hinc', hnb⟩ := h.node_vertex_outside hr hi hT hw hwm hc'
  by_cases hcV : items.type c' = .V
  · obtain ⟨x, hx, rfl⟩ := h.vertItem_of_V hcV
    have hwx := h.v_touch_eq hg hr hqu hx he he' hbe hinc hinc' hnb
    subst hwx
    exact Or.inl ⟨_, h.nodeVerts_get' hi hk, by rw [horig]⟩
  · have hn := isNode_of_SPR hT
    obtain ⟨pos, hl⟩ := (h.node i hi).layout
    set F := (items.ordered g i (t.nvRange (idx i)).1 pos).filter (· ≥ 1 + g.nv) with hF
    have hperm : F.Perm ((items.ch i).filter (· ≥ 1 + g.nv)) := ordered_filter_perm i _ pos _
    have hage : ¬ c' < 1 + g.nv := fun hlt => hcV ((h.type_V_iff hi hc').2 hlt)
    have hmem : c' ∈ F := hperm.mem_iff.2 (List.mem_filter.2 ⟨hc', by simpa using hage⟩)
    obtain ⟨j, hj, hjc⟩ := List.getElem_of_mem hmem
    obtain ⟨x, y, ka, kb, hka, hkb, hxa, hyb, hx, hy, hnvs, -, -⟩ := h.virt_edge hqc hpc hi hn hl hj
    rw [hjc] at hx hy
    have htw := (hl.twin hn j hj).1
    rw [hjc] at htw
    have hcap : items.hasCap c' = true := h.child_hasCap hi hT hc' hcV
    have hcapNe : t.capNe (idx c') = some (t.neRange (idx c')).1 := by
      unfold SpqrTree.capNe; rw [h.hasCap_eq hcs, ite_eq_left hcap]
    have hFlen : F.length = (items.virtualEdges i).length := by
      rw [hperm.length_eq, h.virt_length hi, List.countP_eq_length_filter]
    have hklt : items.capCount i + j < items.nEdges g i := by rw [h.nEdges_node hi hn]; omega
    obtain ⟨-, hget⟩ := h.nodeEdges_get hi hklt
    have hnode := h.nodeEdge_node hi hklt
    have hcF : items.type c' ≠ .F := fun hF => h.not_root_of_child hc' (h.root_of_F hcs hF)
    have hFV : items.type c' ∉ [NodeType.F, .V] := by simp [hcF, hcV]
    have hvs := h.sep_vs hg hcs hFV he he' hinc hinc' hbe hnb
    have hwxy : w = x ∨ w = y := by
      rcases hvs with h1 | h2
      · rw [hx, Option.some.injEq] at h1; exact Or.inl h1.symm
      · rw [h2, Option.getD_some] at hy; exact Or.inr hy
    have hnd := h.nv_nodup hi
    rw [List.getElem?_eq_getElem hka, Option.some.injEq] at hxa
    rw [List.getElem?_eq_getElem hkb, Option.some.injEq] at hyb
    refine Or.inr ⟨_, _, _, hnode, by rw [← Nat.add_assoc]; exact htw, hcapNe, hget, ?_⟩
    rw [hnvs]
    rcases hwxy with rfl | rfl
    · have : k = ka := hnd.getElem_inj_iff.1 (horig.trans hxa.symm)
      subst this
      exact Or.inl rfl
    · have : k = kb := hnd.getElem_inj_iff.1 (horig.trans hyb.symm)
      subst this
      exact Or.inr rfl

/-- `SpqrTree.PieceSep` of the relabelled tree from the item-level facts. -/
theorem pieceSep (hg : g.WF) (hr : Items.Ranges g items σ) (hqu : items.QUpper g)
    (hrs : items.RootSep g) (hrv : items.RootV) (hqc : items.QChildVs g) (hpc : items.PChildVs) :
    t.PieceSep g where
  root_disjoint := h.pieceSep_root_disjoint hrs
  v_attach := h.pieceSep_v_attach hg hr hqu
  v_nonempty := h.pieceSep_v_nonempty
  cap_orig := h.pieceSep_cap_orig
  twin_orient := h.pieceSep_twin_orient hqc hpc
  cap_nonempty := h.pieceSep_cap_nonempty hr
  q_root_attach := h.pieceSep_q_root_attach hg hr
  q_shape := h.pieceSep_q_shape hr
  q_leaf_parent := h.pieceSep_q_leaf_parent hr hrv
  q_loop := h.pieceSep_q_loop hr
  q_upper := h.pieceSep_q_upper hr hqu
  q_lower := h.pieceSep_q_lower hr
  q_cap_orient := h.pieceSep_q_cap_orient hr hqc
  q_lower_attach := h.pieceSep_q_lower_attach hr hg hqu
  nv_child := h.pieceSep_nv_child
  v_child_nv := h.pieceSep_v_child_nv
  node_parent := h.pieceSep_node_parent hrv
  nv_vert_inj := h.pieceSep_nv_vert_inj
  node_child_cap := h.pieceSep_node_child_cap hr
  nv_orig_inj := h.pieceSep_nv_orig_inj
  node_touch := h.pieceSep_node_touch hg hr hqu hqc hpc
  node_attach := h.pieceSep_node_attach hg hr hqu

end Spqr.RelabelOK
