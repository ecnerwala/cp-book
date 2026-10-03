import Spqr.RangesFinal
import Spqr.WalkPlace
import Spqr.ItemAcyc

namespace Spqr
namespace Items

/-- The attachment and shape record of one closed item. -/
structure CloseAt (g : Graph) (items : Items) (i : ItemId) : Prop where
  att_vs : items.type i ∈ [NodeType.S, .P, .R, .Q] →
    ∀ v, items.Att g i v → items.IsVs i v
  vs_att : items.type i ∈ [NodeType.S, .P, .R] ∨ (items.type i = .Q ∧ items.ch i = []) →
    ∀ v, items.IsVs i v → items.Att g i v
  vs_ne : items.type i ∈ [NodeType.S, .P, .R] →
    ∀ u v, items.vs i = (some u, some v) → u ≠ v
  interior : items.type i ∈ [NodeType.S, .P, .R] → ∀ v, v < g.nv →
    (items.IsParent i (vertItem v) ↔
      items.Inner g i v ∧ ∀ c, items.IsParent i c →
        ¬ ∀ e, e < g.ne → g.Inc e v → items.EdgeBelow g c e)
  child_two : ∀ c, items.IsParent i c → items.type i ∈ [NodeType.S, .P, .R] →
    items.type c ≠ .V → ∃ a b, items.vs c = (some a, some b) ∧ a ≠ b
  io_parent : ∀ c, items.IsParent i c → items.type c = .I ∨ items.type c = .O → items.type i = .Q
  q_leaf : ∀ e, i = edgeItem g e → e < g.ne → items.ch i = [] →
    ∃ a b, items.vs i = (some a, some b) ∧ a ≠ b ∧ PairEq (a, b) g.edges[e]!
  q_root : ∀ e, i = edgeItem g e → e < g.ne → items.ch i ≠ [] →
    ∃ u c, items.vs i = (some u, none) ∧ g.Inc e u ∧ items.type c ∉ [NodeType.F, .V] ∧
      (items.type c = .Q → items.ch c = []) ∧
      ((g.edges[e]!).1 = (g.edges[e]!).2 → items.ch i = [c] ∧ items.vs c = (some u, none)) ∧
      ((g.edges[e]!).1 ≠ (g.edges[e]!).2 → ∃ w, w < g.nv ∧ PairEq (u, w) g.edges[e]! ∧
        items.ch i = [c, vertItem w] ∧
        ∃ a b, items.vs c = (some a, some b) ∧ PairEq (a, b) (u, w))
  q_under_v : ∀ v, i = vertItem v → v < g.nv → ∀ c, items.IsParent i c → items.ch c ≠ []
  p_shape : items.type i = .P →
    2 ≤ (items.virtualEdges i).length ∧ (∀ c, items.IsParent i c → items.type c ≠ .V) ∧
    ∀ q ∈ items.virtualEdges i, ∃ u v, items.vs i = (some u, some v) ∧ PairEq q (u, v)
  s_order : items.type i = .S → ∃ u v xs,
    items.vs i = (some u, some v) ∧ ((items.ch i).filter fun c => items.type c = .V) = xs.map vertItem ∧
    1 ≤ xs.length ∧ items.virtualEdges i = List.zip (u :: xs) (xs ++ [v])
  r_shape : items.type i = .R →
    2 ≤ ((items.ch i).filter fun c => items.type c = .V).length ∧ 5 ≤ (items.virtualEdges i).length ∧
    ((items.virtualEdges i).map fun q => min q.1 q.2 + g.nv * max q.1 q.2).Nodup ∧
    ∀ u v, items.vs i = (some u, some v) → ∀ q ∈ items.virtualEdges i, ¬ PairEq q (u, v)

theorem CloseAt.root {g : Graph} {items : Items} (ht : items.type rootItem = .F)
    (hn : items.ch rootItem = []) : CloseAt g items rootItem := by
  constructor
  · simp [ht]
  · simp [ht]
  · simp [ht]
  · simp [ht]
  · simp [ht]
  · simp [hn, IsParent]
  · intro e he; have he' : 0 = 1 + g.nv + e := he; omega
  · intro e he; have he' : 0 = 1 + g.nv + e := he; omega
  · intro v hv; have hv' : 0 = 1 + v := hv; omega
  · simp [ht]
  · simp [ht]
  · simp [ht]

theorem CloseAt.root_of_children {g : Graph} {items : Items} (ht : items.type rootItem = .F)
    (hc : ∀ c, items.IsParent rootItem c → items.type c ≠ .I ∧ items.type c ≠ .O) :
    CloseAt g items rootItem := by
  constructor
  · simp [ht]
  · simp [ht]
  · simp [ht]
  · simp [ht]
  · simp [ht]
  · intro c h hio; exact (hio.elim (hc c h).1 (hc c h).2).elim
  · intro e he; have he' : 0 = 1 + g.nv + e := he; omega
  · intro e he; have he' : 0 = 1 + g.nv + e := he; omega
  · intro v hv; have hv' : 0 = 1 + v := hv; omega
  · simp [ht]
  · simp [ht]
  · simp [ht]

theorem CloseAt.vertex {g : Graph} {items : Items} {v : Nat} (hv : v < g.nv)
    (ht : items.type (vertItem v) = .V)
    (hc : ∀ c, items.IsParent (vertItem v) c → items.type c = .Q ∧ items.ch c ≠ []) :
    CloseAt g items (vertItem v) := by
  constructor
  · simp [ht]
  · simp [ht]
  · simp [ht]
  · simp [ht]
  · simp [ht]
  · intro c h hio; simp [(hc c h).1] at hio
  · intro e he; have he' : 1 + v = 1 + g.nv + e := he; omega
  · intro e he; have he' : 1 + v = 1 + g.nv + e := he; omega
  · intro w hw _ c h; exact (hc c h).2
  · simp [ht]
  · simp [ht]
  · simp [ht]

theorem CloseAt.leafQ {g : Graph} {items : Items} {e a b : Nat}
    (ht : items.type (edgeItem g e) = .Q) (hc : items.ch (edgeItem g e) = [])
    (hvs : items.vs (edgeItem g e) = (some a, some b)) (hne : a ≠ b)
    (hpair : PairEq (a, b) g.edges[e]!)
    (ha : items.Att g (edgeItem g e) a) (hb : items.Att g (edgeItem g e) b) :
    CloseAt g items (edgeItem g e) := by
  constructor
  · intro _ v hv
    obtain ⟨e', _, _, _, hinc, _, hbelow, _⟩ := hv
    have heq := below_eq_of_ch_nil hc hbelow
    have heq' : e' = e := by have : 1 + g.nv + e' = 1 + g.nv + e := heq; omega
    subst e'
    have hinc' : v = a ∨ v = b := by
      rcases hpair with h | h <;> rcases hinc with h' | h'
      · exact Or.inl (h'.symm.trans (congrArg Prod.fst h).symm)
      · exact Or.inr (h'.symm.trans (congrArg Prod.snd h).symm)
      · exact Or.inr (h'.symm.trans (congrArg Prod.snd h).symm)
      · exact Or.inl (h'.symm.trans (congrArg Prod.fst h).symm)
    simpa [IsVs, hvs, eq_comm] using hinc'
  · intro _ v hv
    have hv' : a = v ∨ b = v := by simpa [IsVs, hvs] using hv
    exact hv'.elim (fun h => h ▸ ha) (fun h => h ▸ hb)
  · simp [ht]
  · simp [ht]
  · simp [ht]
  · simp [IsParent, hc]
  · intro e' heq _ _
    have heq' : e = e' := by have : 1 + g.nv + e = 1 + g.nv + e' := heq; omega
    subst e'
    exact ⟨a, b, hvs, hne, hpair⟩
  · simp [hc]
  · intro v h hv; have h' : 1 + g.nv + e = 1 + v := h; omega
  · simp [ht]
  · simp [ht]
  · simp [ht]

theorem CloseAt.leafIO {g : Graph} {items : Items} {i : ItemId}
    (hi : 1 + g.nv + g.ne ≤ i) (ht : items.type i = .I ∨ items.type i = .O)
    (hc : items.ch i = []) : CloseAt g items i := by
  constructor
  · rcases ht with ht | ht <;> simp [ht]
  · rcases ht with ht | ht <;> simp [ht]
  · rcases ht with ht | ht <;> simp [ht]
  · rcases ht with ht | ht <;> simp [ht]
  · rcases ht with ht | ht <;> simp [ht]
  · simp [IsParent, hc]
  · intro e he hel
    have he' : i = 1 + g.nv + e := he
    rw [he'] at hi
    omega
  · simp [hc]
  · intro v he hvl
    have he' : i = 1 + v := he
    rw [he'] at hi
    omega
  · rcases ht with ht | ht <;> simp [ht]
  · rcases ht with ht | ht <;> simp [ht]
  · rcases ht with ht | ht <;> simp [ht]

theorem parent_lt {items : Items} {p c : ItemId} (h : items.IsParent p c) : p < items.size := by
  by_contra hn
  simp [IsParent, ch_of_le _ _ (Nat.le_of_not_lt hn)] at h

/-- The parts of an item and its children observable by its close record. -/
structure CloseFrame (g : Graph) (items items' : Items) (i : ItemId) : Prop where
  type : items'.type i = items.type i
  vs : items'.vs i = items.vs i
  ch : items'.ch i = items.ch i
  child_type : ∀ c, items.IsParent i c → items'.type c = items.type c
  child_vs : ∀ c, items.IsParent i c → items'.vs c = items.vs c
  child_ch : ∀ c, items.IsParent i c → items'.ch c = items.ch c
  edges : ∀ j, j = i ∨ items.IsParent i j → ∀ e,
    items'.EdgeBelow g j e ↔ items.EdgeBelow g j e

theorem CloseFrame.vchildren {g : Graph} {items items' : Items} {i : ItemId}
    (h : CloseFrame g items items' i) :
    (items'.ch i).filter (fun c => items'.type c = .V) =
      (items.ch i).filter (fun c => items.type c = .V) := by
  rw [h.ch]
  apply List.filter_congr
  intro c hc
  rw [h.child_type c hc]

theorem CloseFrame.virtualEdges {g : Graph} {items items' : Items} {i : ItemId}
    (h : CloseFrame g items items' i) : items'.virtualEdges i = items.virtualEdges i := by
  have hc : (items'.ch i).filter (fun c => items'.type c ≠ .V) =
      (items.ch i).filter (fun c => items.type c ≠ .V) := by
    rw [h.ch]
    apply List.filter_congr
    intro c hc
    rw [h.child_type c hc]
  unfold Items.virtualEdges
  rw [hc]
  apply List.map_congr_left
  intro c hc
  rw [h.child_vs c (List.mem_filter.mp hc).1]

theorem CloseAt.frame {g : Graph} {items items' : Items} {i : ItemId}
    (h : CloseAt g items i) (f : CloseFrame g items items' i) : CloseAt g items' i := by
  have hp : ∀ c, items'.IsParent i c ↔ items.IsParent i c := fun c => by simp [IsParent, f.ch]
  have hv : ∀ v, items'.IsVs i v ↔ items.IsVs i v := fun v => by simp [IsVs, f.vs]
  have ha : ∀ v, items'.Att g i v ↔ items.Att g i v := fun v => by
    simp only [Att, f.edges i (Or.inl rfl)]
  have hn : ∀ v, items'.Inner g i v ↔ items.Inner g i v := fun v => by
    simp only [Inner, f.edges i (Or.inl rfl)]
  constructor
  · intro ht v hav
    exact (hv v).2 (h.att_vs (by rwa [f.type] at ht) v ((ha v).1 hav))
  · intro ht v hvv
    exact (ha v).2 (h.vs_att (by simpa only [f.type, f.ch] using ht) v ((hv v).1 hvv))
  · simpa only [f.type, f.vs] using h.vs_ne
  · intro ht v hvg
    rw [hp, hn]
    have hall : (∀ c, items'.IsParent i c → ¬ ∀ e, e < g.ne → g.Inc e v → items'.EdgeBelow g c e) ↔
        (∀ c, items.IsParent i c → ¬ ∀ e, e < g.ne → g.Inc e v → items.EdgeBelow g c e) := by
      simp only [hp]
      apply forall_congr'; intro c
      apply imp_congr_right; intro hc
      simp only [f.edges c (Or.inr hc)]
    rw [hall]
    exact h.interior (by rwa [f.type] at ht) v hvg
  · intro c hc ht hcv
    have hc' := (hp c).1 hc
    simpa only [f.child_vs c hc'] using
      h.child_two c hc' (by rwa [f.type] at ht) (by rwa [f.child_type c hc'] at hcv)
  · intro c hc ht
    have hc' := (hp c).1 hc
    rw [f.type]
    exact h.io_parent c hc' (by rwa [f.child_type c hc'] at ht)
  · simpa only [f.vs, f.ch] using h.q_leaf
  · intro e hei he hch
    obtain ⟨u, c, hu, heuv, hct, hcq, hloop, hnon⟩ := h.q_root e hei he (by rwa [f.ch] at hch)
    have hc : items.IsParent i c := by
      by_cases hl : (g.edges[e]!).1 = (g.edges[e]!).2
      · simp [IsParent, (hloop hl).1]
      · obtain ⟨w, -, -, hw, -⟩ := hnon hl
        simp [IsParent, hw]
    refine ⟨u, c, by rwa [f.vs], heuv, by rwa [f.child_type c hc], ?_, ?_, ?_⟩
    · simpa only [f.child_type c hc, f.child_ch c hc] using hcq
    · simpa only [f.ch, f.child_vs c hc] using hloop
    · simpa only [f.ch, f.child_vs c hc] using hnon
  · intro v hiv hvg c hc
    have hc' := (hp c).1 hc
    rw [f.child_ch c hc']
    exact h.q_under_v v hiv hvg c hc'
  · intro ht
    obtain ⟨he, hcv, hvs⟩ := h.p_shape (by rwa [f.type] at ht)
    rw [f.virtualEdges]
    refine ⟨he, fun c hc => ?_, ?_⟩
    · have hc' := (hp c).1 hc
      rw [f.child_type c hc']; exact hcv c hc'
    · simpa only [f.vs] using hvs
  · simpa only [f.type, f.vs, f.vchildren, f.virtualEdges] using h.s_order
  · simpa only [f.type, f.vs, f.vchildren, f.virtualEdges] using h.r_shape

theorem CloseAt.modify_of_not_below {g : Graph} {items : Items} {i j : ItemId}
    (h : CloseAt g items i) (f : Item → Item) (hj : ¬ items.Below i j) :
    CloseAt g (items.modify j f) i := by
  have hne : i ≠ j := fun hij => hj (by subst j; exact .refl)
  have hn : ∀ c, items.IsParent i c → ¬ items.Below c j :=
    fun c hc hb => hj ((Relation.ReflTransGen.single hc).trans hb)
  have hcn : ∀ c, items.IsParent i c → c ≠ j := fun c hc hcj => hn c hc (by subst j; exact .refl)
  apply h.frame
  refine ⟨?_, vs_modify_of_ne _ _ hne, ch_modify_of_ne _ _ hne,
    ?_, (fun c hc => vs_modify_of_ne _ _ (hcn c hc)),
    (fun c hc => ch_modify_of_ne _ _ (hcn c hc)), ?_⟩
  · simp [type, Array.getElem?_modify, Ne.symm hne]
  · intro c hc; simp [type, Array.getElem?_modify, Ne.symm (hcn c hc)]
  · intro k hk e
    rcases hk with rfl | hk
    · exact Below_modify_of_not_below _ _ hj
    · exact Below_modify_of_not_below _ _ (hn k hk)

end Items

namespace WalkState
open WalkM

/-- Closed records are required for the root, all vertices, and items occurring in a span or child list.
Freshly allocated and reopened items are loose until their next close. -/
structure CloseInv (s : WalkState) : Prop where
  closed : ∀ i, i < s.items.size → i ≤ s.g.nv ∨ 0 < s.cnt i → Items.CloseAt s.g s.items i

theorem init_closeInv (g : Graph) (tern : Bool) : (init g tern).CloseInv := by
  constructor
  intro i hi hc
  rcases hc with hc | hc
  · change i ≤ g.nv at hc
    by_cases hz : i = 0
    · subst i
      exact Items.CloseAt.root (by simp [init, Items.initialItems_type, rootItem])
        (by simp [init, Items.initialItems_ch])
    · have heq : i = vertItem (i - 1) := by dsimp [vertItem]; omega
      rw [heq]
      apply Items.CloseAt.vertex (by dsimp [init] at hc ⊢; omega)
      · simp [init, Items.initialItems_type, vertItem]; omega
      · intro c hh; simp [Items.IsParent, init, Items.initialItems_ch] at hh
  · rw [init_cnt] at hc; omega

theorem CloseInv.frame {s s' : WalkState} (h : s.CloseInv) (hg : s'.g = s.g)
    (hi : s'.items = s.items) (hc : ∀ i, 0 < s'.cnt i → 0 < s.cnt i) : s'.CloseInv := by
  constructor
  intro i hil hl
  rw [hg, hi]
  have hl' : i ≤ s.g.nv ∨ 0 < s'.cnt i := by simpa only [hg] using hl
  exact h.closed i (by rwa [hi] at hil) (hl'.imp_right (hc i))

theorem CloseInv.pop {s : WalkState} (h : s.CloseInv) :
    ({ s with tstack := s.tstack.tail } : WalkState).CloseInv := by
  apply h.frame (s' := { s with tstack := s.tstack.tail }) rfl rfl
  intro i hi
  have := spansCount_tail_le s.tstack i
  dsimp [cnt] at hi ⊢; omega

theorem CloseInv.mergeTop {s : WalkState} (h : s.CloseInv) :
    ({ s with tstack := WalkM.mergeTop s.tstack } : WalkState).CloseInv := by
  apply h.frame (s' := { s with tstack := WalkM.mergeTop s.tstack }) rfl rfl
  intro i hi
  have := spansCount_mergeTop_le s.tstack i
  dsimp [cnt] at hi ⊢; omega

theorem CloseInv.push {s : WalkState} (h : s.CloseInv) (v d i : Nat)
    (hi : Items.CloseAt s.g s.items i) : (after (pushTstack v d i) s).CloseInv := by
  constructor
  intro j hj hl
  by_cases hji : j = i
  · subst hji; exact hi
  · apply h.closed j hj
    rcases hl with hr | hl
    · exact Or.inl hr
    · apply Or.inr
      change 0 < spansCount (_ :: s.tstack) j + chCount s.items j at hl
      rw [spansCount_cons, count_setSides] at hl
      simpa [List.count_singleton, Ne.symm hji, cnt] using hl

theorem CloseInv.pushVert {s : WalkState} (h : s.CloseInv) (v d : Nat)
    (hv : v < s.g.nv) (hs : s.g.nv < s.items.size) :
    (after (pushVertTstack v d) s).CloseInv := by
  apply h.push v d
  apply h.closed (vertItem v) (by dsimp [vertItem]; omega)
  exact Or.inl (by dsimp [vertItem]; omega)

theorem CloseInv.pushEdge {s : WalkState} (h : s.CloseInv) (v d e a b : Nat)
    (ht : Items.type s.items (edgeItem s.g e) = .Q) (hc : Items.ch s.items (edgeItem s.g e) = [])
    (hvs : Items.vs s.items (edgeItem s.g e) = (some a, some b)) (hne : a ≠ b)
    (hp : Items.PairEq (a, b) s.g.edges[e]!)
    (ha : Items.Att s.g s.items (edgeItem s.g e) a) (hb : Items.Att s.g s.items (edgeItem s.g e) b) :
    (after (pushEdgeTstack v d e) s).CloseInv := by
  exact h.push v d _ (Items.CloseAt.leafQ ht hc hvs hne hp ha hb)

theorem CloseInv.alloc {s : WalkState} {P X : ItemId → Prop} (h : s.CloseInv)
    (hp : s.Place s.g P X) (ty : NodeType) : (after (allocItem ty) s).CloseInv := by
  change ({ s with items := s.items.push { type := ty } } : WalkState).CloseInv
  constructor
  intro i hi hc
  have hc' : i ≤ s.g.nv ∨ 0 < s.cnt i := by
    simpa only [cnt, Items.chCount_push s.items { type := ty } rfl] using hc
  have hi' : i < s.items.size := by
    rcases hc' with hc' | hc'
    · have := hp.size; omega
    · exact hp.alloc i hc'
  have hchild : ∀ c, Items.IsParent s.items i c → c < s.items.size := by
    intro c hc
    have hpos := List.count_pos_iff.mpr hc
    have hle := Items.count_le_chCount s.items hi' c
    apply hp.alloc c
    dsimp [cnt]; omega
  apply (h.closed i hi' hc').frame
  refine ⟨Items.type_push_of_ne _ (Nat.ne_of_lt hi'), Items.vs_push_of_ne _ (Nat.ne_of_lt hi'),
    Items.ch_push_nil _ rfl _, ?_, ?_, ?_, ?_⟩
  · intro c hc; exact Items.type_push_of_ne _ (Nat.ne_of_lt (hchild c hc))
  · intro c hc; exact Items.vs_push_of_ne _ (Nat.ne_of_lt (hchild c hc))
  · intro c _; exact Items.ch_push_nil _ rfl _
  · intro j _ e; exact Items.Below_push_nil _ rfl

theorem CloseInv.finishTop {s : WalkState} (h : s.CloseInv) {x : ItemId}
    (hx : x < s.items.size) (hz : s.cnt x = 0) (dir : Bool) (vs : Option Nat × Option Nat)
    (hnew : Items.CloseAt s.g
      (s.items.modify x fun it => { it with vs := vs, ch := getSide s.tstack.head!.spans dir }) x) :
    ({ s with
        items := s.items.modify x fun it => { it with vs := vs, ch := getSide s.tstack.head!.spans dir },
        tstack := match s.tstack with
          | a :: rest => { a with spans := setSides dir [x] [] } :: rest
          | [] => [] } : WalkState).CloseInv := by
  have hn : ∀ p, ¬ Items.IsParent s.items p x := by
    intro p hp
    have hpos := List.count_pos_iff.mpr hp
    have hle := Items.count_le_chCount s.items (Items.parent_lt hp) x
    dsimp [cnt] at hz; omega
  constructor
  intro i hi hc
  by_cases hix : i = x
  · subst i; exact hnew
  · have hle : spansCount (match s.tstack with
          | a :: rest => { a with spans := setSides dir [x] [] } :: rest
          | [] => []) i +
        chCount (s.items.modify x fun it => { it with vs := vs, ch := getSide s.tstack.head!.spans dir }) i
          ≤ s.cnt i := by
      have h1 := Items.chCount_modify s.items hx
        (fun it => { it with vs := vs, ch := getSide s.tstack.head!.spans dir }) i
      have h2 := count_getSide_le s.tstack.head!.spans dir i
      dsimp [cnt]
      cases ht : s.tstack with
      | nil =>
        have hs : getSide ([] : List TEntry).head!.spans dir = [] := by cases dir <;> rfl
        simp only [ht, hs, spansCount_nil, List.count_nil] at h1 ⊢
        omega
      | cons a rest =>
        simp only [ht, head!_cons, spansCount_cons, count_setSides, List.append_nil,
          List.count_singleton, beq_iff_eq, Ne.symm hix, ite_false, Nat.zero_add] at h1 h2 ⊢
        omega
    have hc' : i ≤ s.g.nv ∨ 0 < s.cnt i := hc.imp_right fun hpos => lt_of_lt_of_le hpos hle
    exact (h.closed i (by simpa using hi) hc').modify_of_not_below
      (fun it => { it with vs := vs, ch := getSide s.tstack.head!.spans dir })
      fun hb => hix (Items.Below.eq_of_no_parent hn hb)

theorem CloseInv.unwrap {s : WalkState} (h : s.CloseInv) {a t : TEntry} {rest : List TEntry}
    (hts : s.tstack = a :: t :: rest) {x : ItemId} (hx : x < s.items.size) (dir : Bool) :
    ({ s with tstack := a :: { t with spans := setSides dir (Items.ch s.items x) [] } :: rest } : WalkState).CloseInv := by
  apply h.frame (s' := { s with tstack := a :: { t with spans := setSides dir (Items.ch s.items x) [] } :: rest }) rfl rfl
  intro i hi
  have hc := Items.count_le_chCount s.items hx i
  simp only [cnt, hts, spansCount_cons, count_setSides, List.append_nil] at hi ⊢
  omega

theorem CloseInv.modifyLoose {s : WalkState} (h : s.CloseInv) {x : ItemId}
    (hx : s.g.nv < x) (hz : s.cnt x = 0) (f : Item → Item) (hf : ∀ it, (f it).ch = it.ch) :
    ({ s with items := s.items.modify x f } : WalkState).CloseInv := by
  have hn : ∀ p, ¬ Items.IsParent s.items p x := by
    intro p hp
    have hpos := List.count_pos_iff.mpr hp
    have hle := Items.count_le_chCount s.items (Items.parent_lt hp) x
    dsimp [cnt] at hz; omega
  have hcnt : ∀ i, ({ s with items := s.items.modify x f } : WalkState).cnt i = s.cnt i := by
    intro i
    simp only [cnt, Items.chCount_modify_of_ch s.items x f hf]
  constructor
  intro i hi hc
  rw [hcnt] at hc
  have hix : i ≠ x := by
    intro heq; subst i
    rcases hc with hc | hc
    · exact (Nat.not_le_of_lt hx hc).elim
    · omega
  exact (h.closed i (by simpa using hi) hc).modify_of_not_below f
    fun hb => hix (Items.Below.eq_of_no_parent hn hb)

theorem CloseInv.writeChildren {s : WalkState} (h : s.CloseInv) {x : ItemId}
    (hx : x < s.items.size) (hn : Items.NoParent s.items x) (cs : List ItemId)
    (hnew : Items.CloseAt s.g (s.items.modify x fun it => { it with ch := cs }) x)
    (hcs : ∀ c ∈ cs, Items.CloseAt s.g s.items c) :
    ({ s with items := s.items.modify x fun it => { it with ch := cs } } : WalkState).CloseInv := by
  constructor
  intro i hi hc
  by_cases hix : i = x
  · subst i; exact hnew
  apply Items.CloseAt.modify_of_not_below (f := fun it => { it with ch := cs })
    (hj := fun hb => hix (hn.below_eq hb))
  by_cases hic : i ∈ cs
  · exact hcs i hic
  · apply h.closed i (by simpa using hi)
    rcases hc with hc | hc
    · exact Or.inl hc
    · apply Or.inr
      have hh := Items.chCount_modify s.items hx (fun it => { it with ch := cs }) i
      have hz : cs.count i = 0 := List.count_eq_zero.mpr hic
      dsimp only at hh
      rw [hz, Nat.add_zero] at hh
      dsimp only [cnt] at hc ⊢
      omega

theorem CloseInv.root_append {s : WalkState} {P X : ItemId → Prop} (h : s.CloseInv)
    (hp : s.Place s.g P X) (cs : List ItemId)
    (hcs : ∀ c ∈ cs, ∃ v, v < s.g.nv ∧ c = vertItem v) :
    ({ s with items := s.items.modify rootItem fun it => { it with ch := it.ch ++ cs } } : WalkState).CloseInv := by
  have hr : rootItem < s.items.size := Nat.lt_of_lt_of_le (by show 0 < 1 + s.g.nv + s.g.ne; omega) hp.size
  have hn : Items.NoParent s.items rootItem := by
    intro p hc
    have := List.count_pos_iff.mpr hc
    have hpc := Items.count_le_chCount s.items (Items.parent_lt hc) rootItem
    have hz := hp.root
    dsimp [cnt] at hz
    omega
  have hroot := h.closed rootItem hr (Or.inl (Nat.zero_le _))
  have heq : (s.items.modify rootItem fun it => { it with ch := it.ch ++ cs }) =
      (s.items.modify rootItem fun it => { it with ch := Items.ch s.items rootItem ++ cs }) := by
    apply Array.ext
    · simp
    · intro i hi₁ hi₂
      simp only [Array.getElem_modify]
      split
      · next hEq => subst i; simp [Items.ch, Array.getElem?_eq_getElem hr]
      · rfl
  rw [heq]
  apply h.writeChildren hr hn
  · apply Items.CloseAt.root_of_children
    · rw [Items.type_modify s.items rootItem _ (fun it => { it with ch := Items.ch s.items rootItem ++ cs }) (fun _ => rfl)]
      exact hp.root_type
    · intro c hc
      have hc' : c ∈ Items.ch s.items rootItem ++ cs := by
        simpa only [Items.IsParent, Items.ch_modify_self _ _ _ hr] using hc
      simp only [Items.type_modify s.items rootItem _ (fun it => { it with ch := Items.ch s.items rootItem ++ cs }) (fun _ => rfl)]
      rcases List.mem_append.mp hc' with hc' | hc'
      · constructor <;> intro ht
        · have := hroot.io_parent c hc' (Or.inl ht); simp [hp.root_type] at this
        · have := hroot.io_parent c hc' (Or.inr ht); simp [hp.root_type] at this
      · obtain ⟨v, hv, rfl⟩ := hcs c hc'
        simp [hp.vert v hv]
  · intro c hc
    rcases List.mem_append.mp hc with hc | hc
    · apply h.closed c
      · apply hp.alloc c
        have := Items.count_le_chCount s.items hr c
        have := List.count_pos_iff.mpr hc
        dsimp [cnt]; omega
      · apply Or.inr
        have := Items.count_le_chCount s.items hr c
        have := List.count_pos_iff.mpr hc
        dsimp [cnt]; omega
    · obtain ⟨v, hv, rfl⟩ := hcs c hc
      apply h.closed
      · have := hp.size; dsimp [vertItem]; omega
      · exact Or.inl (by dsimp [vertItem]; omega)

theorem CloseInv.of_tree {s : WalkState} (h : s.CloseInv) (ht : Items.Tree s.g s.items)
    (hty : WalkTyping s.g s.items) : Items.CloseFacts s.g s.items := by
  have hall : ∀ i, i < s.items.size → Items.CloseAt s.g s.items i := by
    intro i hi
    apply h.closed i hi
    by_cases hz : i = rootItem
    · subst i; exact Or.inl (Nat.zero_le _)
    · obtain ⟨p, hp, -⟩ := ht.unique_parent i (Nat.pos_of_ne_zero hz) hi
      have hpos := List.count_pos_iff.mpr hp
      have hle := Items.count_le_chCount s.items (Items.parent_lt hp) i
      exact Or.inr (by dsimp [cnt]; omega)
  have he : ∀ e, e < s.g.ne → edgeItem s.g e < s.items.size := by
    intro e he; have := hty.size; change 1 + s.g.nv + e < s.items.size; omega
  refine ⟨?_, fun i hi => (hall i hi).vs_att, fun i hi => (hall i hi).vs_ne,
    fun i v hi hv ht => (hall i hi).interior ht v hv,
    fun p c hp => (hall p (Items.parent_lt hp)).child_two c hp,
    fun p c hp => (hall p (Items.parent_lt hp)).io_parent c hp,
    fun e he' => (hall _ (he e he')).q_leaf e rfl he',
    fun e he' => (hall _ (he e he')).q_root e rfl he', ?_,
    fun i hi => (hall i hi).p_shape, fun i hi => (hall i hi).s_order,
    fun i hi => (hall i hi).r_shape⟩
  · intro e he' v
    exact (hall _ (he e he')).att_vs (by rw [hty.edge e he']; simp) v
  · intro v c hv hp
    exact (hall _ (Items.parent_lt hp)).q_under_v v rfl hv c hp

end WalkState
end Spqr
