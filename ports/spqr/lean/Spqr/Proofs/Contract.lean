import Spqr.Contract
import Spqr.Proofs.SepPairExhaust
import Mathlib.Data.List.Nodup

/-!
# Separation pairs of a skeleton

`Pieces.WF.sepPair_contract_lift`, `Pieces.WF.sepPair_contract_of`,
`Pieces.WF.threeConnected_contract_iff` / `threeConnected_contract_iff_dfs` (PROOF.md §4.5/§5), and
`Graph.TwoAttached.sepClass_mem` (separation classes respect a 2-attached edge set).
-/

namespace Spqr

namespace Graph

variable {g : Graph}

theorem joins_iff {e x y : Nat} :
    g.Joins e x y ↔ e < g.ne ∧ (g.edges[e]! = (x, y) ∨ g.edges[e]! = (y, x)) := by
  unfold Joins
  constructor
  · rintro (h | h) <;> obtain ⟨hlt, h⟩ := Array.getElem?_eq_some_iff.1 h
    · exact ⟨hlt, .inl (by rw [getElem!_pos g.edges e hlt, h])⟩
    · exact ⟨hlt, .inr (by rw [getElem!_pos g.edges e hlt, h])⟩
  · rintro ⟨hlt, h | h⟩
    · exact .inl (Array.getElem?_eq_some_iff.2
        ⟨hlt, by rw [getElem!_pos g.edges e (show e < g.edges.size from hlt)] at h; exact h⟩)
    · exact .inr (Array.getElem?_eq_some_iff.2
        ⟨hlt, by rw [getElem!_pos g.edges e (show e < g.edges.size from hlt)] at h; exact h⟩)

theorem isEnd_iff {e v : Nat} : g.IsEnd e v ↔ e < g.ne ∧ g.Inc e v := by
  constructor
  · rintro ⟨w, h⟩
    obtain ⟨hlt, h | h⟩ := joins_iff.1 h
    · exact ⟨hlt, .inl (by rw [h])⟩
    · exact ⟨hlt, .inr (by rw [h])⟩
  · rintro ⟨hlt, h | h⟩
    · exact ⟨(g.edges[e]!).2, joins_iff.2 ⟨hlt, .inl (by rw [← h])⟩⟩
    · exact ⟨(g.edges[e]!).1, joins_iff.2 ⟨hlt, .inr (by rw [← h])⟩⟩

theorem joins_fst {e : Nat} (he : e < g.ne) : g.Joins e (g.edges[e]!).1 (g.edges[e]!).2 :=
  joins_iff.2 ⟨he, .inl rfl⟩

theorem Touches.isEnd {E : Nat → Prop} {v : Nat} (h : g.Touches E v) :
    ∃ e, E e ∧ g.IsEnd e v := by
  obtain ⟨e, he, hE, hv⟩ := h
  exact ⟨e, hE, isEnd_iff.2 ⟨he, hv⟩⟩

theorem Touches.of_isEnd {E : Nat → Prop} {e v : Nat} (hE : E e) (h : g.IsEnd e v) :
    g.Touches E v :=
  ⟨e, (isEnd_iff.1 h).1, hE, (isEnd_iff.1 h).2⟩

theorem AdjIn.adj {E : Nat → Prop} {u w : Nat} (h : g.AdjIn E u w) :
    ∃ e, E e ∧ g.Joins e u w := by
  obtain ⟨e, he, hE, hp | hp⟩ := h
  · exact ⟨e, hE, joins_iff.2 ⟨he, .inl hp.symm⟩⟩
  · exact ⟨e, hE, joins_iff.2 ⟨he, .inr (by rw [Prod.ext_iff] at hp; exact Prod.ext hp.2.symm hp.1.symm)⟩⟩

/-- A walk along edges of `E` is a walk in `g` whose vertices all touch `E`. -/
theorem reach_of_adjIn {E ok : Nat → Prop} {u w : Nat}
    (h : Relation.ReflTransGen (g.AdjIn E) u w) (hok : ∀ v, g.Touches E v → ok v)
    (hu : g.Touches E u) : g.Reach ok u w := by
  induction h with
  | refl => exact .refl (hok _ hu)
  | tail _ hadj ih =>
    obtain ⟨e, hE, hj⟩ := hadj.adj
    exact ih.tail hj.adj (hok _ (Touches.of_isEnd hE hj.symm.isEnd))

theorem ConnEdges.reach {E ok : Nat → Prop} (hc : g.ConnEdges E) {u w : Nat}
    (hu : g.Touches E u) (hw : g.Touches E w) (hok : ∀ v, g.Touches E v → ok v) :
    g.Reach ok u w := by
  obtain ⟨e, he, hE, hv⟩ := hu
  obtain ⟨e', he', hE', hv'⟩ := hw
  have h1 := reach_of_adjIn (reach_of_inc he hE hv) hok ⟨e, he, hE, .inl rfl⟩
  have h2 := reach_of_adjIn (hc e e' he he' hE hE') hok ⟨e, he, hE, .inl rfl⟩
  have h3 := reach_of_adjIn (reach_of_inc he' hE' hv') hok ⟨e', he', hE', .inl rfl⟩
  exact h1.symm.trans (h2.trans h3)

/-- A walk avoiding the terminals of a 2-attached `E` that starts at a vertex of `E` stays at
vertices of `E`. -/
theorem TwoAttached.reach_touches {E ok : Nat → Prop} {a b u w : Nat} (h : g.TwoAttached E a b)
    (ha : ¬ok a) (hb : ¬ok b) (hr : g.Reach ok u w) (hu : g.Touches E u) : g.Touches E w := by
  induction hr with
  | refl => exact hu
  | tail hr hadj _ ih =>
    obtain ⟨e, hj⟩ := hadj
    obtain ⟨e₁, he₁, hE₁, hv₁⟩ := ih
    have he := (isEnd_iff.1 hj.isEnd)
    have hok := hr.ok_right
    by_cases hE : E e
    · exact Touches.of_isEnd hE hj.symm.isEnd
    · rcases h _ e₁ e he₁ he.1 hE₁ hE hv₁ he.2 with h' | h'
      · rw [h'] at hok; exact absurd hok ha
      · rw [h'] at hok; exact absurd hok hb

/-- Classes avoiding the terminals of a 2-attached `E` are inside `E` or disjoint from it. -/
theorem TwoAttached.edgeConn_mem {E ok : Nat → Prop} {a b e e' : Nat} (h : g.TwoAttached E a b)
    (ha : ¬ok a) (hb : ¬ok b) (he : E e) (hc : g.EdgeConn ok e e') : E e' := by
  rcases hc with rfl | ⟨u, w, hu, hw, hr⟩
  · exact he
  · obtain ⟨e₁, he₁, hE₁, hv₁⟩ := h.reach_touches ha hb hr (Touches.of_isEnd he hu)
    by_cases hE : E e'
    · exact hE
    · have hw' := isEnd_iff.1 hw
      have hok := hr.ok_right
      rcases h _ e₁ e' he₁ hw'.1 hE₁ hE hv₁ hw'.2 with h' | h'
      · rw [h'] at hok; exact absurd hok ha
      · rw [h'] at hok; exact absurd hok hb

/-- If `E` is 2-attached at `{a, b}`, every separation class of `{a, b}` is inside `E` or disjoint
from it. -/
theorem TwoAttached.sepClass_mem {E : Nat → Prop} {a b e e' : Nat} (h : g.TwoAttached E a b)
    (hc : g.SepClass a b e e') : (E e ↔ E e') :=
  ⟨fun he => h.edgeConn_mem (fun h => h.1 rfl) (fun h => h.2 rfl) he hc,
    fun he => h.edgeConn_mem (fun h => h.1 rfl) (fun h => h.2 rfl) he hc.symm⟩

end Graph

namespace Pieces

variable {g : Graph} {P : Pieces}

theorem ne_contract : (P.contract g).ne = (P.origins g).length := by
  simp [contract, Graph.ne]

theorem origins_nodup : (P.origins g).Nodup := by
  refine List.nodup_append.2 ⟨?_, ?_, ?_⟩
  · exact List.Nodup.map (fun _ _ h => Sum.inl.inj h) (List.nodup_range.filter _)
  · exact List.Nodup.map (fun _ _ h => Sum.inr.inj h) List.nodup_range
  · intro s hs s' hs'
    obtain ⟨_, -, rfl⟩ := List.mem_map.1 hs
    obtain ⟨_, -, rfl⟩ := List.mem_map.1 hs'
    exact Sum.inl_ne_inr

theorem mem_origins_inl {e : Nat} :
    Sum.inl e ∈ P.origins g ↔ e < g.ne ∧ P.piece e = none := by
  simp only [origins, List.mem_append, List.mem_map, List.mem_filter, List.mem_range,
    Sum.inl.injEq, exists_eq_right, Option.isNone_iff_eq_none, reduceCtorEq, and_false,
    exists_false, or_false]

theorem mem_origins_inr {i : Nat} : Sum.inr i ∈ P.origins g ↔ i < P.k := by
  simp only [origins, List.mem_append, List.mem_map, List.mem_filter, List.mem_range,
    Sum.inr.injEq, exists_eq_right, reduceCtorEq, and_false, exists_false, false_or]

theorem Orig.lt {f : Nat} {s : Nat ⊕ Nat} (h : P.Orig g f s) : f < (P.contract g).ne := by
  rw [ne_contract]; exact (List.getElem?_eq_some_iff.1 h).1

theorem Orig.inj {f f' : Nat} {s : Nat ⊕ Nat} (h : P.Orig g f s) (h' : P.Orig g f' s) : f = f' :=
  (origins_nodup.getElem?_inj (List.getElem?_eq_some_iff.1 h).1).1 (h.trans h'.symm)

theorem Orig.fn {f : Nat} {s s' : Nat ⊕ Nat} (h : P.Orig g f s) (h' : P.Orig g f s') : s = s' :=
  Option.some.inj (h.symm.trans h')

theorem orig_total {f : Nat} (hf : f < (P.contract g).ne) : ∃ s, P.Orig g f s := by
  rw [ne_contract] at hf
  exact ⟨_, List.getElem?_eq_getElem hf⟩

theorem orig_inl_iff {e : Nat} :
    (∃ f, P.Orig g f (.inl e)) ↔ e < g.ne ∧ P.piece e = none := by
  rw [← mem_origins_inl, List.mem_iff_getElem?]; rfl

theorem orig_inr_iff {i : Nat} : (∃ f, P.Orig g f (.inr i)) ↔ i < P.k := by
  rw [← mem_origins_inr, List.mem_iff_getElem?]; rfl

theorem Orig.inl_lt {f e : Nat} (h : P.Orig g f (.inl e)) : e < g.ne ∧ P.piece e = none :=
  orig_inl_iff.1 ⟨f, h⟩

theorem Orig.inr_lt {f i : Nat} (h : P.Orig g f (.inr i)) : i < P.k :=
  orig_inr_iff.1 ⟨f, h⟩

theorem contract_edges_getElem? {f : Nat} :
    (P.contract g).edges[f]? = ((P.origins g)[f]?).map (P.originEnds g) := by
  simp [contract]

theorem joins_contract_iff {f u v : Nat} :
    (P.contract g).Joins f u v ↔
      ∃ s, P.Orig g f s ∧ (P.originEnds g s = (u, v) ∨ P.originEnds g s = (v, u)) := by
  unfold Graph.Joins Orig
  rw [contract_edges_getElem?]
  constructor
  · rintro (h | h) <;> obtain ⟨s, hs, h⟩ := Option.map_eq_some_iff.1 h
    · exact ⟨s, hs, .inl h⟩
    · exact ⟨s, hs, .inr h⟩
  · rintro ⟨s, hs, h | h⟩
    · exact .inl (by rw [hs, Option.map_some, h])
    · exact .inr (by rw [hs, Option.map_some, h])

theorem Orig.joins_inl {f e u v : Nat} (h : P.Orig g f (.inl e)) :
    (P.contract g).Joins f u v ↔ g.Joins e u v := by
  rw [joins_contract_iff, Graph.joins_iff]
  constructor
  · rintro ⟨s, hs, hv⟩
    obtain rfl := h.fn hs
    exact ⟨h.inl_lt.1, hv⟩
  · rintro ⟨-, hv⟩
    exact ⟨_, h, hv⟩

theorem Orig.joins_inr {f i u v : Nat} (h : P.Orig g f (.inr i)) :
    (P.contract g).Joins f u v ↔ (u = P.x i ∧ v = P.y i) ∨ (u = P.y i ∧ v = P.x i) := by
  rw [joins_contract_iff]
  constructor
  · rintro ⟨s, hs, hv⟩
    obtain rfl := h.fn hs
    simp only [originEnds, Prod.ext_iff] at hv
    rcases hv with ⟨h1, h2⟩ | ⟨h1, h2⟩
    · exact .inl ⟨h1.symm, h2.symm⟩
    · exact .inr ⟨h2.symm, h1.symm⟩
  · rintro (⟨rfl, rfl⟩ | ⟨rfl, rfl⟩)
    · exact ⟨_, h, .inl rfl⟩
    · exact ⟨_, h, .inr rfl⟩

theorem Orig.isEnd_inr {f i : Nat} (h : P.Orig g f (.inr i)) :
    (P.contract g).IsEnd f (P.x i) ∧ (P.contract g).IsEnd f (P.y i) :=
  ⟨⟨_, h.joins_inr.2 (.inl ⟨rfl, rfl⟩)⟩, ⟨_, h.joins_inr.2 (.inr ⟨rfl, rfl⟩)⟩⟩

theorem Orig.isEnd_inl {f e u : Nat} (h : P.Orig g f (.inl e)) :
    (P.contract g).IsEnd f u ↔ g.IsEnd e u :=
  ⟨fun ⟨v, hv⟩ => ⟨v, h.joins_inl.1 hv⟩, fun ⟨v, hv⟩ => ⟨v, h.joins_inl.2 hv⟩⟩

theorem Orig.isEnd_inr_iff {f i u : Nat} (h : P.Orig g f (.inr i)) :
    (P.contract g).IsEnd f u ↔ u = P.x i ∨ u = P.y i := by
  constructor
  · rintro ⟨v, hv⟩
    rcases h.joins_inr.1 hv with ⟨h, -⟩ | ⟨h, -⟩
    · exact .inl h
    · exact .inr h
  · rintro (rfl | rfl)
    · exact h.isEnd_inr.1
    · exact h.isEnd_inr.2

theorem Mem.fn {i j e : Nat} (h : P.Mem i e) (h' : P.Mem j e) : i = j :=
  Option.some.inj (h.symm.trans h')

theorem Mem.ne_none {i e : Nat} (h : P.Mem i e) : P.piece e ≠ none := by
  rw [Mem] at h; rw [h]; exact Option.some_ne_none _

theorem not_mem_of_none {e : Nat} (h : P.piece e = none) (i : Nat) : ¬P.Mem i e := by
  intro h'; rw [Mem, h] at h'; exact Option.some_ne_none _ h'.symm

theorem exists_mem_of_ne_none {e : Nat} (h : P.piece e ≠ none) : ∃ i, P.Mem i e :=
  Option.ne_none_iff_exists'.1 h

theorem Img.lt {e f : Nat} (h : P.Img g e f) : f < (P.contract g).ne :=
  h.elim (fun h => h.2.lt) fun ⟨_, _, h⟩ => h.lt

theorem Img.fn {e f f' : Nat} (h : P.Img g e f) (h' : P.Img g e f') : f = f' := by
  rcases h with ⟨hn, h⟩ | ⟨i, hi, h⟩ <;> rcases h' with ⟨hn', h'⟩ | ⟨j, hj, h'⟩
  · exact h.inj h'
  · exact absurd hn hj.ne_none
  · exact absurd hn' hi.ne_none
  · obtain rfl := hi.fn hj
    exact h.inj h'

theorem Img.inj {e e' f : Nat} (h : P.Img g e f) (h' : P.Img g e' f) :
    e = e' ∨ ∃ i, P.Mem i e ∧ P.Mem i e' := by
  rcases h with ⟨hn, h⟩ | ⟨i, hi, h⟩ <;> rcases h' with ⟨hn', h'⟩ | ⟨j, hj, h'⟩
  · exact .inl (Sum.inl.inj (h.fn h'))
  · exact absurd (h.fn h') Sum.inl_ne_inr
  · exact absurd (h.fn h') Sum.inr_ne_inl
  · obtain rfl := Sum.inr.inj (h.fn h')
    exact .inr ⟨i, hi, hj⟩

theorem Img.orig_of_mem {e f i : Nat} (h : P.Img g e f) (hi : P.Mem i e) :
    P.Orig g f (.inr i) := by
  rcases h with ⟨hn, -⟩ | ⟨j, hj, h⟩
  · exact absurd hn hi.ne_none
  · obtain rfl := hi.fn hj
    exact h

theorem Img.isEnd_of_none {e f u : Nat} (h : P.Img g e f) (hn : P.piece e = none) :
    (P.contract g).IsEnd f u ↔ g.IsEnd e u := by
  rcases h with ⟨-, h⟩ | ⟨i, hi, -⟩
  · exact h.isEnd_inl
  · exact absurd hn hi.ne_none

theorem not_skel {w : Nat} (h : ¬P.Skel g w) : ∃ i, P.Int g i w := by
  by_contra h'
  exact h fun i hi => h' ⟨i, hi⟩

/-- A skeleton vertex touching piece `i` is one of its terminals. -/
theorem terminal_of_skel {i w : Nat} (hs : P.Skel g w) (hk : i < P.k)
    (ht : g.Touches (P.Mem i) w) : w = P.x i ∨ w = P.y i := by
  by_contra h
  exact hs i ⟨hk, ht, fun h' => h (.inl h'), fun h' => h (.inr h')⟩

theorem terminal_other {i t₀ t₁ v : Nat}
    (ht : (t₀ = P.x i ∧ t₁ = P.y i) ∨ (t₀ = P.y i ∧ t₁ = P.x i))
    (hv : v = P.x i ∨ v = P.y i) (hne : v ≠ t₀) : v = t₁ := by
  rcases ht with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ <;> rcases hv with h | h <;>
    first | exact h | exact absurd h hne

section WF

variable (hP : P.WF g)
include hP

theorem WF.mem_lt {i e : Nat} (h : P.Mem i e) : e < g.ne := (hP.lt i e h).2
theorem WF.mem_k {i e : Nat} (h : P.Mem i e) : i < P.k := (hP.lt i e h).1

/-- Edges at an interior vertex of piece `i` lie in `i`. -/
theorem WF.int_mem {i w e : Nat} (hi : P.Int g i w) (he : g.IsEnd e w) : P.Mem i e := by
  obtain ⟨hk, ⟨e₁, he₁, hE₁, hv₁⟩, hx, hy⟩ := hi
  by_contra hE
  have he' := Graph.isEnd_iff.1 he
  rcases hP.attached i hk w e₁ e he₁ he'.1 hE₁ hE hv₁ he'.2 with h | h
  · exact hx h
  · exact hy h

/-- A vertex incident to an edge outside piece `i` is a skeleton vertex as far as `i` is
concerned; so endpoints of edges in no piece, and terminals, are skeleton vertices. -/
theorem WF.skel_of_isEnd {w e : Nat} (he : g.IsEnd e w)
    (hE : ∀ i, P.Mem i e → w = P.x i ∨ w = P.y i) : P.Skel g w := by
  rintro i ⟨hk, ⟨e₁, he₁, hE₁, hv₁⟩, hx, hy⟩
  have he' := Graph.isEnd_iff.1 he
  by_cases hi : P.Mem i e
  · rcases hE i hi with h | h
    · exact hx h
    · exact hy h
  · rcases hP.attached i hk w e₁ e he₁ he'.1 hE₁ hi hv₁ he'.2 with h | h
    · exact hx h
    · exact hy h

theorem WF.skel_of_none {w e : Nat} (he : g.IsEnd e w) (hE : P.piece e = none) : P.Skel g w :=
  hP.skel_of_isEnd he fun i h => absurd h (not_mem_of_none hE i)

theorem WF.skel_x {i : Nat} (hk : i < P.k) : P.Skel g (P.x i) := by
  obtain ⟨e, hE, he⟩ := (hP.touch i hk).1.isEnd
  exact hP.skel_of_isEnd he fun j hj => by obtain rfl := hE.fn hj; exact .inl rfl

theorem WF.skel_y {i : Nat} (hk : i < P.k) : P.Skel g (P.y i) := by
  obtain ⟨e, hE, he⟩ := (hP.touch i hk).2.isEnd
  exact hP.skel_of_isEnd he fun j hj => by obtain rfl := hE.fn hj; exact .inr rfl

theorem WF.skel_end {f u : Nat} (h : (P.contract g).IsEnd f u) : P.Skel g u := by
  obtain ⟨s, hs⟩ := orig_total h.lt
  rcases s with e | i
  · exact hP.skel_of_none (hs.isEnd_inl.1 h) hs.inl_lt.2
  · rcases hs.isEnd_inr_iff.1 h with rfl | rfl
    · exact hP.skel_x hs.inr_lt
    · exact hP.skel_y hs.inr_lt

/-- The image of `e` is incident to every skeleton endpoint of `e`. -/
theorem WF.img_isEnd {e f u : Nat} (h : P.Img g e f) (hu : g.IsEnd e u) (hs : P.Skel g u) :
    (P.contract g).IsEnd f u := by
  rcases h with ⟨hn, h⟩ | ⟨i, hi, h⟩
  · exact h.isEnd_inl.2 hu
  · exact h.isEnd_inr_iff.2
      (terminal_of_skel hs (hP.mem_k hi) (Graph.Touches.of_isEnd hi hu))

theorem WF.img_total {e : Nat} (he : e < g.ne) : ∃ f, P.Img g e f := by
  by_cases hn : P.piece e = none
  · obtain ⟨f, hf⟩ := orig_inl_iff.2 ⟨he, hn⟩
    exact ⟨f, .inl ⟨hn, hf⟩⟩
  · obtain ⟨i, hi⟩ := exists_mem_of_ne_none hn
    obtain ⟨f, hf⟩ := orig_inr_iff.2 (hP.mem_k hi)
    exact ⟨f, .inr ⟨i, hi, hf⟩⟩

theorem WF.pre_total {f : Nat} (hf : f < (P.contract g).ne) : ∃ e, e < g.ne ∧ P.Img g e f := by
  obtain ⟨s, hs⟩ := orig_total hf
  rcases s with e | i
  · exact ⟨e, hs.inl_lt.1, .inl ⟨hs.inl_lt.2, hs⟩⟩
  · obtain ⟨e, he, hE, -⟩ := (hP.touch i hs.inr_lt).1
    exact ⟨e, he, .inr ⟨i, hE, hs⟩⟩

/-- A walk of `g` from a skeleton vertex lifts to the skeleton: to its end if that is a skeleton
vertex, else to a live terminal of the piece it is interior to. -/
theorem WF.reach_contract_aux {ok : Nat → Prop} {u w : Nat} (hr : g.Reach ok u w)
    (hu : P.Skel g u) :
    (P.Skel g w → (P.contract g).Reach ok u w) ∧
      ∀ i, P.Int g i w → ∃ t, (t = P.x i ∨ t = P.y i) ∧ ok t ∧ (P.contract g).Reach ok u t := by
  induction hr with
  | refl hok => exact ⟨fun _ => .refl hok, fun i hi => absurd hi (hu i)⟩
  | @tail w z hr hadj hok ih =>
    obtain ⟨e, hj⟩ := hadj
    have hw := hr.ok_right
    have step : ∀ {i t}, i < P.k → (t = P.x i ∨ t = P.y i) → (z = P.x i ∨ z = P.y i) →
        (P.contract g).Reach ok u t → (P.contract g).Reach ok u z := by
      intro i t hk ht hz hrt
      by_cases htz : t = z
      · subst htz; exact hrt
      · obtain ⟨f, hf⟩ := orig_inr_iff.2 hk
        refine hrt.tail ⟨f, hf.joins_inr.2 ?_⟩ hok
        rcases ht with h | h <;> rcases hz with h' | h'
        · exact absurd (h.trans h'.symm) htz
        · exact .inl ⟨h, h'⟩
        · exact .inr ⟨h, h'⟩
        · exact absurd (h.trans h'.symm) htz
    constructor
    · intro hz
      by_cases hws : P.Skel g w
      · have h1 := ih.1 hws
        by_cases hE : P.piece e = none
        · obtain ⟨f, hf⟩ := orig_inl_iff.2 ⟨hj.lt, hE⟩
          exact h1.tail ⟨f, hf.joins_inl.2 hj⟩ hok
        · obtain ⟨i, hi⟩ := exists_mem_of_ne_none hE
          have hk := hP.mem_k hi
          exact step hk (terminal_of_skel hws hk (Graph.Touches.of_isEnd hi hj.isEnd))
            (terminal_of_skel hz hk (Graph.Touches.of_isEnd hi hj.symm.isEnd)) h1
      · obtain ⟨i, hi⟩ := not_skel hws
        obtain ⟨t, ht, -, hrt⟩ := ih.2 i hi
        have hE := hP.int_mem hi hj.isEnd
        exact step hi.1 ht (terminal_of_skel hz hi.1 (Graph.Touches.of_isEnd hE hj.symm.isEnd))
          hrt
    · intro j hj'
      have hE := hP.int_mem hj' hj.symm.isEnd
      by_cases hws : P.Skel g w
      · exact ⟨w, terminal_of_skel hws hj'.1 (Graph.Touches.of_isEnd hE hj.isEnd), hw, ih.1 hws⟩
      · obtain ⟨i, hi⟩ := not_skel hws
        obtain rfl := (hP.int_mem hi hj.isEnd).fn hE
        exact ih.2 i hi

theorem WF.reach_contract {ok : Nat → Prop} {u w : Nat} (hr : g.Reach ok u w)
    (hu : P.Skel g u) (hw : P.Skel g w) : (P.contract g).Reach ok u w :=
  (hP.reach_contract_aux hr hu).1 hw

/-- Terminals of a piece are joined inside it when every removed vertex is a skeleton vertex and
both terminals are alive. -/
theorem WF.reach_terminals {ok : Nat → Prop} (hok : ∀ v, ¬ok v → P.Skel g v) {i : Nat}
    (hk : i < P.k) (hx : ok (P.x i)) (hy : ok (P.y i)) : g.Reach ok (P.x i) (P.y i) :=
  (hP.conn i hk).reach (hP.touch i hk).1 (hP.touch i hk).2 fun v hv => by
    by_contra h
    rcases terminal_of_skel (hok v h) hk hv with rfl | rfl
    · exact h hx
    · exact h hy

/-- A walk of the skeleton is a walk of `g`. -/
theorem WF.reach_of_contract {ok : Nat → Prop} (hok : ∀ v, ¬ok v → P.Skel g v) {u w : Nat}
    (hr : (P.contract g).Reach ok u w) : g.Reach ok u w := by
  induction hr with
  | refl h => exact .refl h
  | tail hr hadj hokz ih =>
    obtain ⟨f, hj⟩ := hadj
    obtain ⟨s, hs⟩ := orig_total hj.lt
    rcases s with e | i
    · exact ih.tail ⟨e, hs.joins_inl.1 hj⟩ hokz
    · rcases hs.joins_inr.1 hj with ⟨h1, h2⟩ | ⟨h1, h2⟩
      · rw [h2] at hokz ⊢; rw [h1] at ih hr
        exact ih.trans (hP.reach_terminals hok hs.inr_lt hr.ok_right hokz)
      · rw [h2] at hokz ⊢; rw [h1] at ih hr
        exact ih.trans (hP.reach_terminals hok hs.inr_lt hokz hr.ok_right).symm

/-- In a block, every edge of piece `i` reaches the terminal `t₁` inside the piece avoiding the
other terminal `t₀`, and `t₁` has an edge outside the piece. -/
theorem WF.exit_edge (h2 : g.TwoConnected) {i t₀ t₁ : Nat} (hk : i < P.k)
    (ht : (t₀ = P.x i ∧ t₁ = P.y i) ∨ (t₀ = P.y i ∧ t₁ = P.x i)) {e : Nat} (he : P.Mem i e) :
    (∃ v, g.IsEnd e v ∧ g.Reach (fun w => w ≠ t₀ ∧ g.Touches (P.Mem i) w) v t₁) ∧
      ∃ e', e' < g.ne ∧ ¬P.Mem i e' ∧ g.IsEnd e' t₁ := by
  obtain ⟨e', he', hE'⟩ := hP.proper i hk
  rcases h2 t₀ e e' (hP.mem_lt he) he' with rfl | ⟨v, v', hv, hv', hr⟩
  · exact absurd he hE'
  rcases hr.exit (Graph.Touches.of_isEnd he hv) with hin | ⟨u, z, hin, hadj, -, hnz⟩
  · have hv'' := Graph.isEnd_iff.1 hv'
    obtain ⟨e₁, he₁, hE₁, hv₁⟩ := hin.ok_right.2
    have hvt := terminal_other ht (hP.attached i hk v' e₁ e' he₁ hv''.1 hE₁ hE' hv₁ hv''.2)
      hin.ok_right.1
    rw [hvt] at hin hv'
    exact ⟨⟨v, hv, hin⟩, e', he', hE', hv'⟩
  · obtain ⟨e₀, hj⟩ := hadj
    have hE₀ : ¬P.Mem i e₀ := fun h => hnz (Graph.Touches.of_isEnd h hj.symm.isEnd)
    obtain ⟨e₁, he₁, hE₁, hv₁⟩ := hin.ok_right.2
    have hj' := Graph.isEnd_iff.1 hj.isEnd
    have hut := terminal_other ht (hP.attached i hk u e₁ e₀ he₁ hj'.1 hE₁ hE₀ hv₁ hj'.2)
      hin.ok_right.1
    have hj0 := hj.isEnd
    rw [hut] at hin hj0
    exact ⟨⟨v, hv, hin⟩, e₀, hj'.1, hE₀, hj0⟩

/-- A piece with a live terminal has a live terminal `t`, incident to an edge outside the piece,
reached from each of its edges by live vertices. -/
theorem WF.live_terminal (h2 : g.TwoConnected) {ok : Nat → Prop}
    (hok : ∀ v, ¬ok v → P.Skel g v) {i : Nat} (hk : i < P.k)
    (hlive : ok (P.x i) ∨ ok (P.y i)) :
    ∃ t, (t = P.x i ∨ t = P.y i) ∧ ok t ∧ (∃ e', e' < g.ne ∧ ¬P.Mem i e' ∧ g.IsEnd e' t) ∧
      ∀ e, P.Mem i e → ∃ v, g.IsEnd e v ∧ g.Reach ok v t := by
  have key : ∀ t₀ t₁, ((t₀ = P.x i ∧ t₁ = P.y i) ∨ (t₀ = P.y i ∧ t₁ = P.x i)) → ok t₁ →
      (∀ v, v ≠ t₀ → g.Touches (P.Mem i) v → ok v) →
      ∃ t, (t = P.x i ∨ t = P.y i) ∧ ok t ∧ (∃ e', e' < g.ne ∧ ¬P.Mem i e' ∧ g.IsEnd e' t) ∧
        ∀ e, P.Mem i e → ∃ v, g.IsEnd e v ∧ g.Reach ok v t := by
    intro t₀ t₁ ht h1 hall
    refine ⟨t₁, ht.elim (fun h => .inr h.2) (fun h => .inl h.2), h1, ?_, fun e he => ?_⟩
    · obtain ⟨e, -, he, -⟩ := (hP.touch i hk).1
      exact (hP.exit_edge h2 hk ht he).2
    · obtain ⟨v, hv, hr⟩ := (hP.exit_edge h2 hk ht he).1
      exact ⟨v, hv, hr.mono fun w hw => hall w hw.1 hw.2⟩
  by_cases hx : ok (P.x i)
  · refine key (P.y i) (P.x i) (.inr ⟨rfl, rfl⟩) hx fun v hv htv => ?_
    by_contra h
    rcases terminal_of_skel (hok v h) hk htv with rfl | rfl
    · exact h hx
    · exact hv rfl
  · have hy := hlive.resolve_left hx
    refine key (P.x i) (P.y i) (.inl ⟨rfl, rfl⟩) hy fun v hv htv => ?_
    by_contra h
    rcases terminal_of_skel (hok v h) hk htv with rfl | rfl
    · exact hv rfl
    · exact h hy

/-- The virtual edge of a piece with a live terminal is in a nontrivial class. -/
theorem WF.virtual_nontrivial (h2 : g.TwoConnected) {ok : Nat → Prop}
    (hok : ∀ v, ¬ok v → P.Skel g v) {i f : Nat} (hf : P.Orig g f (.inr i))
    (hlive : ok (P.x i) ∨ ok (P.y i)) :
    ∃ f', f' < (P.contract g).ne ∧ f' ≠ f ∧ (P.contract g).EdgeConn ok f f' := by
  have hk := hf.inr_lt
  obtain ⟨t, ht, hokt, ⟨e', he', hE', hte⟩, -⟩ := hP.live_terminal h2 hok hk hlive
  obtain ⟨f', hf'⟩ := hP.img_total he'
  have hts : P.Skel g t := ht.elim (fun h => h ▸ hP.skel_x hk) (fun h => h ▸ hP.skel_y hk)
  refine ⟨f', hf'.lt, fun h => ?_, .of_reach (hf.isEnd_inr_iff.2 ht) (hP.img_isEnd hf' hte hts)
    (.refl hokt)⟩
  subst h
  rcases hf' with ⟨-, h⟩ | ⟨j, hj, h⟩
  · exact Sum.inl_ne_inr (h.fn hf)
  · obtain rfl := Sum.inr.inj (h.fn hf)
    exact hE' hj

/-- Edge classes of `g` map to edge classes of the skeleton. -/
theorem WF.edgeConn_contract (h2 : g.TwoConnected) {ok : Nat → Prop}
    (hok : ∀ v, ¬ok v → P.Skel g v) {e e' f f' : Nat} (he : P.Img g e f) (he' : P.Img g e' f')
    (hc : g.EdgeConn ok e e') : (P.contract g).EdgeConn ok f f' := by
  have aux : ∀ {e f v}, P.Img g e f → g.IsEnd e v → ok v →
      (∃ u, P.Skel g u ∧ (P.contract g).IsEnd f u ∧ g.Reach ok u v) ∨
        ∃ i, P.Mem i e ∧ ¬ok (P.x i) ∧ ¬ok (P.y i) := by
    intro e f v hi hv hokv
    rcases hi with ⟨hn, h⟩ | ⟨i, hi, h⟩
    · exact .inl ⟨v, hP.skel_of_none hv hn, h.isEnd_inl.2 hv, .refl hokv⟩
    · by_cases hlive : ok (P.x i) ∨ ok (P.y i)
      · have hk := hP.mem_k hi
        obtain ⟨t, ht, hokt, -, hr⟩ := hP.live_terminal h2 hok hk hlive
        obtain ⟨v₀, hv₀, hr⟩ := hr e hi
        refine .inl ⟨t, ht.elim (fun h => h ▸ hP.skel_x hk) (fun h => h ▸ hP.skel_y hk),
          h.isEnd_inr_iff.2 ht, hr.symm.trans ?_⟩
        by_cases hv₀v : v₀ = v
        · subst hv₀v; exact .refl hokv
        · obtain ⟨w, hw⟩ := hv₀
          obtain ⟨w', hw'⟩ := hv
          rcases hw.eq_or hw' with ⟨h1, -⟩ | ⟨h1, h2⟩
          · exact absurd h1 hv₀v
          · exact .tail (.refl hr.ok_left) ⟨e, h2 ▸ hw⟩ hokv
      · exact .inr ⟨i, hi, fun h => hlive (.inl h), fun h => hlive (.inr h)⟩
  have dead : ∀ {e e' f f' i}, P.Img g e f → P.Img g e' f' → g.EdgeConn ok e e' → P.Mem i e →
      ¬ok (P.x i) → ¬ok (P.y i) → f = f' := by
    intro e e' f f' i he he' hc hi hx hy
    have hi' := (hP.attached i (hP.mem_k hi)).edgeConn_mem hx hy hi hc
    exact (he.orig_of_mem hi).inj (he'.orig_of_mem hi')
  rcases hc with rfl | ⟨v, v', hv, hv', hr⟩
  · rw [he.fn he']; exact .refl _
  have hc : g.EdgeConn ok e e' := .inr ⟨v, v', hv, hv', hr⟩
  rcases aux he hv hr.ok_left with ⟨u, hus, hu, hr₁⟩ | ⟨i, hi, hx, hy⟩
  · rcases aux he' hv' hr.ok_right with ⟨u', hus', hu', hr₂⟩ | ⟨i, hi, hx, hy⟩
    · exact .of_reach hu hu' (hP.reach_contract (hr₁.trans (hr.trans hr₂.symm)) hus hus')
    · rw [dead he' he hc.symm hi hx hy]; exact .refl _
  · rw [dead he he' hc hi hx hy]; exact .refl _

/-- Edge classes of the skeleton map back to edge classes of `g`, for edges whose piece (if any)
has a live terminal. -/
theorem WF.edgeConn_of_contract (h2 : g.TwoConnected) {ok : Nat → Prop}
    (hok : ∀ v, ¬ok v → P.Skel g v) {e e' f f' : Nat} (he : P.Img g e f) (he' : P.Img g e' f')
    (hl : ∀ i, P.Mem i e → ok (P.x i) ∨ ok (P.y i))
    (hl' : ∀ i, P.Mem i e' → ok (P.x i) ∨ ok (P.y i))
    (hc : (P.contract g).EdgeConn ok f f') : g.EdgeConn ok e e' := by
  have aux : ∀ {e f u}, P.Img g e f → (∀ i, P.Mem i e → ok (P.x i) ∨ ok (P.y i)) →
      (P.contract g).IsEnd f u → ok u → ∃ v, g.IsEnd e v ∧ g.Reach ok v u := by
    intro e f u hi hl hu hoku
    rcases hi with ⟨hn, h⟩ | ⟨i, hi, h⟩
    · exact ⟨u, h.isEnd_inl.1 hu, .refl hoku⟩
    · have hk := hP.mem_k hi
      obtain ⟨t, ht, hokt, -, hr⟩ := hP.live_terminal h2 hok hk (hl i hi)
      obtain ⟨v, hv, hr⟩ := hr e hi
      refine ⟨v, hv, hr.trans ?_⟩
      rcases h.isEnd_inr_iff.1 hu with rfl | rfl <;> rcases ht with rfl | rfl
      · exact .refl hokt
      · exact (hP.reach_terminals hok hk hoku hokt).symm
      · exact hP.reach_terminals hok hk hokt hoku
      · exact .refl hokt
  rcases hc with rfl | ⟨u, u', hu, hu', hr⟩
  · rcases he.inj he' with rfl | ⟨i, hi, hi'⟩
    · exact .refl _
    · obtain ⟨t, -, -, -, ht⟩ := hP.live_terminal h2 hok (hP.mem_k hi) (hl i hi)
      obtain ⟨v, hv, hr⟩ := ht e hi
      obtain ⟨v', hv', hr'⟩ := ht e' hi'
      exact .of_reach hv hv' (hr.trans hr'.symm)
  · obtain ⟨v, hv, hvu⟩ := aux he hl hu hr.ok_left
    obtain ⟨v', hv', hvu'⟩ := aux he' hl' hu' hr.ok_right
    exact .of_reach hv hv' (hvu.trans ((hP.reach_of_contract hok hr).trans hvu'.symm))

end WF

theorem skel_ok {a b : Nat} (ha : P.Skel g a) (hb : P.Skel g b) :
    ∀ v, ¬(v ≠ a ∧ v ≠ b) → P.Skel g v := by
  intro v hv
  by_cases h : v = a
  · subst h; exact ha
  by_cases h' : v = b
  · subst h'; exact hb
  exact absurd ⟨h, h'⟩ hv

section Main

variable (hP : P.WF g) (h2 : g.TwoConnected)
include hP h2

/-- (1) A separation pair of the skeleton whose vertices are skeleton vertices is a separation
pair of `g`. -/
theorem WF.sepPair_contract_lift {a b : Nat} (ha : P.Skel g a) (hb : P.Skel g b)
    (h : (P.contract g).SeparationPair a b) : g.SeparationPair a b := by
  have hok := skel_ok ha hb
  obtain ⟨hab, h⟩ := h
  refine ⟨hab, ?_⟩
  have key : ∀ {f f' e e'}, P.Img g e f → P.Img g e' f' → ¬(P.contract g).SepClass a b f f' →
      ¬g.SepClass a b e e' :=
    fun he he' hn hc => hn (hP.edgeConn_contract h2 hok he he' hc)
  have key2 : ∀ {f f' e e'}, P.Img g e f → P.Img g e' f' → f ≠ f' →
      (P.contract g).SepClass a b f f' → e ≠ e' ∧ g.SepClass a b e e' := by
    intro f f' e e' he he' hne hc
    have live : ∀ {f e}, P.Img g e f → (∃ u, (P.contract g).IsEnd f u ∧ u ≠ a ∧ u ≠ b) →
        ∀ i, P.Mem i e → (P.x i ≠ a ∧ P.x i ≠ b) ∨ (P.y i ≠ a ∧ P.y i ≠ b) := by
      intro f e he hu i hi
      obtain ⟨u, hu, hoku⟩ := hu
      rcases (he.orig_of_mem hi).isEnd_inr_iff.1 hu with rfl | rfl
      · exact .inl hoku
      · exact .inr hoku
    rcases hc with h | ⟨u, u', hu, hu', hr⟩
    · exact absurd h hne
    refine ⟨fun h => hne (by subst h; exact he.fn he'), ?_⟩
    exact hP.edgeConn_of_contract h2 hok he he' (live he ⟨u, hu, hr.ok_left⟩)
      (live he' ⟨u', hu', hr.ok_right⟩) (.inr ⟨u, u', hu, hu', hr⟩)
  rcases h with ⟨f₁, f₂, f₃, l1, l2, l3, h12, h13, h23⟩ |
    ⟨f₁, f₁', f₂, f₂', l1, l1', l2, l2', hne₁, hc₁, hne₂, hc₂, h12⟩
  · obtain ⟨e₁, he₁, hi₁⟩ := hP.pre_total l1
    obtain ⟨e₂, he₂, hi₂⟩ := hP.pre_total l2
    obtain ⟨e₃, he₃, hi₃⟩ := hP.pre_total l3
    exact .inl ⟨e₁, e₂, e₃, he₁, he₂, he₃, key hi₁ hi₂ h12, key hi₁ hi₃ h13, key hi₂ hi₃ h23⟩
  · obtain ⟨e₁, he₁, hi₁⟩ := hP.pre_total l1
    obtain ⟨e₁', he₁', hi₁'⟩ := hP.pre_total l1'
    obtain ⟨e₂, he₂, hi₂⟩ := hP.pre_total l2
    obtain ⟨e₂', he₂', hi₂'⟩ := hP.pre_total l2'
    obtain ⟨n1, c1⟩ := key2 hi₁ hi₁' hne₁ hc₁
    obtain ⟨n2, c2⟩ := key2 hi₂ hi₂' hne₂ hc₂
    exact .inr ⟨e₁, e₁', e₂, e₂', he₁, he₁', he₂, he₂', n1, c1, n2, c2, key hi₁ hi₂ h12⟩

/-- (2) A separation pair of `g` whose vertices are skeleton vertices and which is not the
terminal pair of a piece is a separation pair of the skeleton. -/
theorem WF.sepPair_contract_of {a b : Nat} (ha : P.Skel g a) (hb : P.Skel g b)
    (hnt : ¬P.TermPair a b) (h : g.SeparationPair a b) : (P.contract g).SeparationPair a b := by
  have hok := skel_ok ha hb
  obtain ⟨hab, h⟩ := h
  refine ⟨hab, ?_⟩
  have live : ∀ i, i < P.k → (P.x i ≠ a ∧ P.x i ≠ b) ∨ (P.y i ≠ a ∧ P.y i ≠ b) := by
    intro i hk
    by_contra hc
    have hx : P.x i = a ∨ P.x i = b := by
      by_contra h'
      exact hc (.inl ⟨fun h => h' (.inl h), fun h => h' (.inr h)⟩)
    have hy : P.y i = a ∨ P.y i = b := by
      by_contra h'
      exact hc (.inr ⟨fun h => h' (.inl h), fun h => h' (.inr h)⟩)
    refine hnt ⟨i, hk, ?_⟩
    rcases hx with hx | hx <;> rcases hy with hy | hy
    · exact absurd (hx.trans hy.symm) (hP.ne i hk)
    · exact .inl ⟨hx.symm, hy.symm⟩
    · exact .inr ⟨hy.symm, hx.symm⟩
    · exact absurd (hx.trans hy.symm) (hP.ne i hk)
  have live' : ∀ e i, P.Mem i e → (P.x i ≠ a ∧ P.x i ≠ b) ∨ (P.y i ≠ a ∧ P.y i ≠ b) :=
    fun _ i hi => live i (hP.mem_k hi)
  have key : ∀ {e e' f f'}, P.Img g e f → P.Img g e' f' → ¬g.SepClass a b e e' →
      ¬(P.contract g).SepClass a b f f' :=
    fun he he' hn hc => hn (hP.edgeConn_of_contract h2 hok he he' (live' _) (live' _) hc)
  have key2 : ∀ {e e' f}, P.Img g e f → e ≠ e' → e' < g.ne → g.SepClass a b e e' →
      ∃ f', f' < (P.contract g).ne ∧ f ≠ f' ∧ (P.contract g).SepClass a b f f' := by
    intro e e' f he hne he' hc
    obtain ⟨f', hf'⟩ := hP.img_total he'
    by_cases hff : f = f'
    · subst hff
      rcases he.inj hf' with h | ⟨i, hi, -⟩
      · exact absurd h hne
      · obtain ⟨f'', l, n, c⟩ :=
          hP.virtual_nontrivial h2 hok (he.orig_of_mem hi) (live i (hP.mem_k hi))
        exact ⟨f'', l, Ne.symm n, c⟩
    · exact ⟨f', hf'.lt, hff, hP.edgeConn_contract h2 hok he hf' hc⟩
  rcases h with ⟨e₁, e₂, e₃, l1, l2, l3, h12, h13, h23⟩ |
    ⟨e₁, e₁', e₂, e₂', l1, l1', l2, l2', hne₁, hc₁, hne₂, hc₂, h12⟩
  · obtain ⟨f₁, hf₁⟩ := hP.img_total l1
    obtain ⟨f₂, hf₂⟩ := hP.img_total l2
    obtain ⟨f₃, hf₃⟩ := hP.img_total l3
    exact .inl ⟨f₁, f₂, f₃, hf₁.lt, hf₂.lt, hf₃.lt, key hf₁ hf₂ h12, key hf₁ hf₃ h13,
      key hf₂ hf₃ h23⟩
  · obtain ⟨f₁, hf₁⟩ := hP.img_total l1
    obtain ⟨f₂, hf₂⟩ := hP.img_total l2
    obtain ⟨f₁', m1, n1, c1⟩ := key2 hf₁ hne₁ l1' hc₁
    obtain ⟨f₂', m2, n2, c2⟩ := key2 hf₂ hne₂ l2' hc₂
    exact .inr ⟨f₁, f₁', f₂, f₂', hf₁.lt, m1, hf₂.lt, m2, n1, c1, n2, c2, key hf₁ hf₂ h12⟩

/-- The skeleton of a block is a block. -/
theorem WF.twoConnected_contract : (P.contract g).TwoConnected := by
  intro v f f' hf hf'
  obtain ⟨e, he, hi⟩ := hP.pre_total hf
  obtain ⟨e', he', hi'⟩ := hP.pre_total hf'
  have hc := h2 v e e' he he'
  by_cases hv : P.Skel g v
  · exact hP.edgeConn_contract h2 (fun w hw => by obtain rfl := not_not.1 hw; exact hv) hi hi' hc
  · have hc' := hP.edgeConn_contract h2 (ok := fun _ => True) (fun _ h => absurd trivial h) hi hi'
      (hc.mono fun _ _ => trivial)
    rcases hc' with h | ⟨u, u', hu, hu', hr⟩
    · exact .inl h
    · have hr' : (P.contract g).Reach (fun w => True ∧ w ≠ v) u u' :=
        hr.and_of_adj (fun h => hv (h ▸ hP.skel_end hu))
          (fun _ z hadj h => by obtain ⟨f₀, hj⟩ := hadj; exact hv (h ▸ hP.skel_end hj.symm.isEnd))
      exact .inr ⟨u, u', hu, hu', hr'.mono fun _ h => h.2⟩

/-- A vertex interior to a piece is not on the skeleton, so is in no separation pair of it. -/
theorem WF.not_sepPair_of_not_skel {a b : Nat} (ha : ¬P.Skel g a) :
    ¬(P.contract g).SeparationPair a b := by
  intro h
  obtain ⟨f, f', hf, hf', hn⟩ := h.exists_not_sepClass
  apply hn
  rcases hP.twoConnected_contract h2 b f f' hf hf' with h | ⟨u, u', hu, hu', hr⟩
  · exact .inl h
  · have hr' : (P.contract g).Reach (fun w => w ≠ b ∧ w ≠ a) u u' :=
      hr.and_of_adj (fun h => ha (h ▸ hP.skel_end hu))
        (fun _ z hadj h => by obtain ⟨f₀, hj⟩ := hadj; exact ha (h ▸ hP.skel_end hj.symm.isEnd))
    exact .inr ⟨u, u', hu, hu', hr'.mono fun _ h => ⟨h.2, h.1⟩⟩

/-- (3) The skeleton is 3-connected iff no separation pair of `g` has both vertices on the
skeleton other than the terminal pairs of its pieces, and no terminal pair separates the
skeleton. -/
theorem WF.threeConnected_contract_iff :
    (P.contract g).ThreeConnected ↔
      (∀ a b, P.Skel g a → P.Skel g b → ¬P.TermPair a b → ¬g.SeparationPair a b) ∧
        ∀ i, i < P.k → ¬(P.contract g).SeparationPair (P.x i) (P.y i) := by
  constructor
  · intro h
    exact ⟨fun a b ha hb hnt hs => h a b (hP.sepPair_contract_of h2 ha hb hnt hs),
      fun i _ => h _ _⟩
  · rintro ⟨h, ht⟩ a b hs
    by_cases ha : P.Skel g a
    · by_cases hb : P.Skel g b
      · by_cases hnt : P.TermPair a b
        · obtain ⟨i, hk, ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩⟩ := hnt
          · exact ht i hk hs
          · exact ht i hk hs.symm
        · exact h a b ha hb hnt (hP.sepPair_contract_lift h2 ha hb hs)
      · exact hP.not_sepPair_of_not_skel h2 hb hs.symm
    · exact hP.not_sepPair_of_not_skel h2 ha hs

/-- With `sepPair_iff'`: an R skeleton over a lowpoint-sorted DFS tree is 3-connected iff no
type-1/type-2 pair of the block has both vertices on the skeleton (terminal pairs of pieces
excepted) and no terminal pair separates the skeleton. -/
theorem WF.threeConnected_contract_iff_dfs {d : DfsData} (hs : d.Spec g) (hr : d.Rooted g) :
    (P.contract g).ThreeConnected ↔
      (∀ a b, P.Skel g a → P.Skel g b → ¬P.TermPair a b →
          d.Anc d.root a → d.Anc d.root b →
          ¬(d.Type1Pair a b g ∨ d.Type2Pair a b) ∧ ¬(d.Type1Pair b a g ∨ d.Type2Pair b a)) ∧
        ∀ i, i < P.k → ¬(P.contract g).SeparationPair (P.x i) (P.y i) := by
  rw [hP.threeConnected_contract_iff h2]
  refine and_congr_left' ⟨fun h a b ha hb hnt hra hrb => ?_, fun h a b ha hb hnt hsep => ?_⟩
  · exact not_or.1 fun h' => h a b ha hb hnt ((d.sepPair_iff' hs hr h2 hra hrb).2 h')
  · obtain ⟨e, e', he, he', hn⟩ := hsep.exists_not_sepClass
    by_cases hra : d.Anc d.root a
    · by_cases hrb : d.Anc d.root b
      · exact not_or.2 (h a b ha hb hnt hra hrb) ((d.sepPair_iff' hs hr h2 hra hrb).1 hsep)
      · exact hn (Graph.sepClass_comm.2 (d.sepClass_of_not_vertex hr h2 hrb he he'))
    · exact hn (d.sepClass_of_not_vertex hr h2 hra he he')

end Main

end Pieces

end Spqr
