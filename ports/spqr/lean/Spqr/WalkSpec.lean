import Mathlib.Tactic.Set
import Spqr.GraphLemmas
import Spqr.WalkTyping

/-!
# The walk invariant: soundness and completeness per step

Stated on `walkTree` directly (no ear structure). Every tstack entry `t` describes a connected set
of original edges `t.edges` whose attachment to the rest of the block is exactly the two terminals
`t.vStart` and `stackVerts[t.topDepth]`. *Soundness*: `mergeTstackTops` joins two entries sharing a
terminal, so the union is again connected and 2-attached. *Completeness*: every item `finishEdge`
closes is 2-attached and connected, i.e. cut off by a genuine 2-separation, and the closed items
partition the edges. The stack-shape facts these steps rely on (`size ≥ origTstack + 3`, `nxt`
being the chain's entry) are proved on the ear layer (`EarSpec.lean`).
-/

namespace Spqr

namespace Items

variable {items items' : Items}

theorem IsParent_congr (h : ∀ p, items'.ch p = items.ch p) {p c : ItemId} :
    items'.IsParent p c ↔ items.IsParent p c := by
  simp [IsParent, h]

theorem Below_congr (h : ∀ p, items'.ch p = items.ch p) {a i : ItemId} :
    items'.Below a i ↔ items.Below a i := by
  have : items'.IsParent = items.IsParent := funext fun _ => funext fun _ => propext (IsParent_congr h)
  simp [Below, this]

theorem ch_push_nil (x : Item) (hx : x.ch = []) (p : ItemId) : Items.ch (items.push x) p = items.ch p := by
  simp only [ch, Array.getElem?_push]
  split
  · subst p; simp [hx]
  · rfl

theorem Below_push_nil (x : Item) (hx : x.ch = []) {a i : ItemId} :
    Items.Below (items.push x) a i ↔ items.Below a i :=
  Below_congr (ch_push_nil x hx)

theorem ch_modify_of_ne (j : ItemId) (f : Item → Item) {p : ItemId} (h : p ≠ j) :
    Items.ch (items.modify j f) p = items.ch p := by
  simp [ch, Array.getElem?_modify, h.symm]

theorem ch_modify_at (j : ItemId) (f : Item → Item) (hj : j < items.size) :
    Items.ch (items.modify j f) j = (f items[j]).ch := by
  simp [ch, Array.getElem?_modify, Array.getElem?_eq_getElem hj]

theorem vs_modify_at (j : ItemId) (f : Item → Item) (hj : j < items.size) :
    Items.vs (items.modify j f) j = (f items[j]).vs := by
  simp [vs, Array.getElem?_modify, Array.getElem?_eq_getElem hj]

theorem vs_modify_of_ne (j : ItemId) (f : Item → Item) {p : ItemId} (h : p ≠ j) :
    Items.vs (items.modify j f) p = items.vs p := by
  simp [vs, Array.getElem?_modify, h.symm]

theorem Below_modify_ch_eq (j : ItemId) (f : Item → Item) (hf : ∀ it, (f it).ch = it.ch) {a i : ItemId} :
    Items.Below (items.modify j f) a i ↔ items.Below a i := by
  refine Below_congr fun p => ?_
  by_cases h : p = j
  · subst h
    by_cases hj : p < items.size
    · rw [ch_modify_at p f hj, hf]; simp [ch, Array.getElem?_eq_getElem hj]
    · simp [ch, Array.getElem?_modify, Array.getElem?_eq_none (Nat.le_of_not_lt hj)]
  · exact ch_modify_of_ne j f h

/-- Modifying an item outside the subtree of `a` does not change that subtree. -/
theorem Below_modify_of_not_below (j : ItemId) (f : Item → Item) {a i : ItemId}
    (hj : ¬ items.Below a j) : Items.Below (items.modify j f) a i ↔ items.Below a i := by
  constructor <;> intro h
  · induction h with
    | refl => exact .refl
    | @tail k l _ hkl ih =>
      have hk : k ≠ j := fun hkj => hj (hkj ▸ ih)
      exact ih.tail (by simpa [IsParent, ch_modify_of_ne j f hk] using hkl)
  · induction h with
    | refl => exact .refl
    | @tail k l hak hkl ih =>
      have hk : k ≠ j := fun hkj => hj (hkj ▸ hak)
      exact ih.tail (by simpa [IsParent, ch_modify_of_ne j f hk] using hkl)

theorem Below.eq_of_no_parent {a j : ItemId} (hj : ∀ p, ¬ items.IsParent p j) (h : items.Below a j) :
    a = j := by
  rcases h.cases_tail with h | ⟨c, _, hc⟩
  · exact h.symm
  · exact absurd hc (hj c)

theorem Below.head_cases {a i : ItemId} (h : items.Below a i) :
    a = i ∨ ∃ c, items.IsParent a c ∧ items.Below c i := by
  rcases h.cases_head with h | ⟨c, hc, h⟩
  · exact .inl h
  · exact .inr ⟨c, hc, h⟩

end Items

theorem mem_setSides (dir : Bool) (a : List α) (i : α) :
    i ∈ (setSides dir a []).1 ++ (setSides dir a []).2 ↔ i ∈ a := by
  cases dir <;> simp [setSides]

theorem mem_of_getSide_nil (dir : Bool) (p : List α × List α) (h : getSide p (!dir) = []) (i : α) :
    i ∈ p.1 ++ p.2 ↔ i ∈ getSide p dir := by
  cases dir <;> simp_all [getSide]

theorem setSides_eq (dir : Bool) (a b u v : α) (h : setSides dir a b = (u, v)) :
    (u = a ∧ v = b) ∨ (u = b ∧ v = a) := by
  cases dir <;> simp_all [setSides]

namespace TEntry

/-- The original edges covered by an entry: everything below either side's items. -/
def edges (g : Graph) (items : Items) (t : TEntry) (e : Nat) : Prop :=
  ∃ i ∈ t.spans.1 ++ t.spans.2, items.EdgeBelow g i e

/-- Second terminal, read off the vertex stack. -/
def top (s : WalkState) (t : TEntry) : Nat := s.stackVerts[t.topDepth]!

variable {g : Graph} {items items' : Items} {t t' : TEntry}

theorem edges_of_spans_eq (h : t.spans = t'.spans) : t.edges g items = t'.edges g items := by
  funext e; unfold edges; rw [h]

theorem edges_congr (h : ∀ i ∈ t.spans.1 ++ t.spans.2, ∀ e, items'.EdgeBelow g i e ↔ items.EdgeBelow g i e)
    (e : Nat) : t.edges g items' e ↔ t.edges g items e := by
  simp only [edges]
  constructor <;> rintro ⟨i, hi, hb⟩
  · exact ⟨i, hi, (h i hi e).1 hb⟩
  · exact ⟨i, hi, (h i hi e).2 hb⟩

theorem edges_modify_of_not_mem (j : ItemId) (f : Item → Item) (hj : ∀ p, ¬ items.IsParent p j)
    (hmem : j ∉ t.spans.1 ++ t.spans.2) (e : Nat) :
    t.edges g (items.modify j f) e ↔ t.edges g items e :=
  edges_congr (fun _ hi _ => Items.Below_modify_of_not_below j f fun hb =>
    hmem ((Items.Below.eq_of_no_parent hj hb) ▸ hi)) e

theorem top_congr {s s' : WalkState} (h : s'.stackVerts = s.stackVerts) : t.top s' = t.top s := by
  unfold TEntry.top; rw [h]

end TEntry

namespace WalkState

/-- Invariant on one tstack entry. -/
structure EntryInv (s : WalkState) (t : TEntry) : Prop where
  conn : s.g.ConnEdges (t.edges s.g s.items)
  attached : s.g.TwoAttached (t.edges s.g s.items) t.vStart (t.top s)

/-- Invariant on a closed item (its subtree is a 2-attached connected piece). -/
structure ItemInv (s : WalkState) (i : ItemId) : Prop where
  conn : s.g.ConnEdges (Items.EdgeBelow s.g s.items i)
  attached : ∀ u v, Items.vs s.items i = (some u, some v) →
    s.g.TwoAttached (Items.EdgeBelow s.g s.items i) u v

/-- The walk invariant: every open entry and every allocated node is a 2-attached connected piece. -/
structure Inv (s : WalkState) : Prop where
  entries : ∀ t ∈ s.tstack, s.EntryInv t
  nodes : ∀ i, 1 + s.g.nv + s.g.ne ≤ i → i < s.items.size → s.ItemInv i

theorem EntryInv.congr {s s' : WalkState} {t : TEntry} (hg : s'.g = s.g)
    (hsv : s'.stackVerts = s.stackVerts)
    (hE : ∀ e, e < s.g.ne → (t.edges s.g s'.items e ↔ t.edges s.g s.items e)) (h : s.EntryInv t) :
    s'.EntryInv t := by
  obtain ⟨conn, att⟩ := h
  refine ⟨?_, ?_⟩ <;> rw [hg]
  · exact (Graph.ConnEdges.congr hE).2 conn
  · rw [TEntry.top_congr hsv]; exact (Graph.TwoAttached.congr hE).2 att

theorem ItemInv.congr {s s' : WalkState} {i : ItemId} (hg : s'.g = s.g)
    (hvs : Items.vs s'.items i = Items.vs s.items i)
    (hE : ∀ e, e < s.g.ne → (Items.EdgeBelow s.g s'.items i e ↔ Items.EdgeBelow s.g s.items i e))
    (h : s.ItemInv i) : s'.ItemInv i := by
  obtain ⟨conn, att⟩ := h
  refine ⟨?_, ?_⟩ <;> rw [hg]
  · exact (Graph.ConnEdges.congr hE).2 conn
  · intro u v huv; rw [hvs] at huv; exact (Graph.TwoAttached.congr hE).2 (att u v huv)

end WalkState

open WalkM

theorem mergeTstackTops_run_eq (s : WalkState) (cur nxt : TEntry) (rest : List TEntry)
    (hs : s.tstack = cur :: nxt :: rest) :
    mergeTstackTops.run s = ((), { s with tstack :=
      { nxt with topDepth := min nxt.topDepth cur.topDepth,
                 spans := (cur.spans.1 ++ nxt.spans.1, nxt.spans.2 ++ cur.spans.2) } :: rest }) := by
  obtain ⟨g, tern, items, sv, sd, nei, fo, tstack, tb, tsl⟩ := s
  subst hs; rfl

/-- Soundness of a merge. `cur` (the top) is merged into `nxt`; the result has terminals
`nxt.vStart` and `stackVerts[min nxt.topDepth cur.topDepth]`. Hypotheses: the two edge sets share a
vertex (unless one is empty), and every old terminal is a new terminal or becomes interior to the
union. The three uses in `finishEdge` are: series (`nxt.top = cur.vStart`, interior), parallel
(same two terminals), rigid (`nxt.top = cur.top`, `cur.vStart` interior). -/
theorem mergeTstackTops_sound (s : WalkState) (cur nxt : TEntry) (rest : List TEntry)
    (hs : s.tstack = cur :: nxt :: rest) (h : ∀ t ∈ s.tstack, s.EntryInv t)
    (hx : (∃ e, e < s.g.ne ∧ cur.edges s.g s.items e) → (∃ e, e < s.g.ne ∧ nxt.edges s.g s.items e) →
      ∃ x, s.g.Touches (cur.edges s.g s.items) x ∧ s.g.Touches (nxt.edges s.g s.items) x)
    (hterm : ∀ v ∈ [cur.vStart, cur.top s, nxt.vStart, nxt.top s],
      v = nxt.vStart ∨ v = s.stackVerts[min nxt.topDepth cur.topDepth]! ∨
        s.g.Interior (fun e => cur.edges s.g s.items e ∨ nxt.edges s.g s.items e) v) :
    let s' := (mergeTstackTops.run s).2
    ∀ t ∈ s'.tstack, s'.EntryInv t := by
  intro s' t ht
  have hs' : s' = { s with tstack :=
      { nxt with topDepth := min nxt.topDepth cur.topDepth,
                 spans := (cur.spans.1 ++ nxt.spans.1, nxt.spans.2 ++ cur.spans.2) } :: rest } := by
    simp only [s', mergeTstackTops_run_eq s cur nxt rest hs]
  rw [hs'] at ht ⊢
  simp only [List.mem_cons] at ht
  have hcur := h cur (by simp [hs])
  have hnxt := h nxt (by simp [hs])
  rcases ht with rfl | ht
  · have hm : ∀ e, e < s.g.ne → (TEntry.edges s.g s.items
        { nxt with topDepth := min nxt.topDepth cur.topDepth,
                   spans := (cur.spans.1 ++ nxt.spans.1, nxt.spans.2 ++ cur.spans.2) } e ↔
        (cur.edges s.g s.items e ∨ nxt.edges s.g s.items e)) := by
      intro e _
      simp only [TEntry.edges, List.mem_append]
      constructor
      · rintro ⟨i, hi, hb⟩
        rcases hi with (hi | hi) | (hi | hi)
        · exact .inl ⟨i, .inl hi, hb⟩
        · exact .inr ⟨i, .inl hi, hb⟩
        · exact .inr ⟨i, .inr hi, hb⟩
        · exact .inl ⟨i, .inr hi, hb⟩
      · rintro (⟨i, hi | hi, hb⟩ | ⟨i, hi | hi, hb⟩)
        · exact ⟨i, .inl (.inl hi), hb⟩
        · exact ⟨i, .inr (.inr hi), hb⟩
        · exact ⟨i, .inl (.inr hi), hb⟩
        · exact ⟨i, .inr (.inl hi), hb⟩
    refine ⟨(Graph.ConnEdges.congr hm).2 (Graph.ConnEdges.union hcur.conn hnxt.conn hx), ?_⟩
    exact (Graph.TwoAttached.congr hm).2 (Graph.TwoAttached.union hcur.attached hnxt.attached hterm)
  · exact WalkState.EntryInv.congr (s := s) rfl rfl (fun _ _ => Iff.rfl) (h t (by simp [hs, ht]))

theorem makeVs_run (s : WalkState) (vStart topDepth : Nat) :
    (makeVs vStart topDepth).run s =
      (setSides s.stackDir[topDepth]! (some s.stackVerts[topDepth]!) (some vStart), s) := rfl

theorem finishTstackTop_run_eq (s : WalkState) (item : ItemId) (t : TEntry) (rest : List TEntry)
    (hs : s.tstack = t :: rest) :
    (finishTstackTop item).run s = ((), { s with
      items := s.items.modify item fun it =>
        { it with vs := setSides s.stackDir[t.topDepth]! (some s.stackVerts[t.topDepth]!) (some t.vStart),
                  ch := getSide t.spans s.stackDir[t.topDepth]! },
      tstack := { t with spans := setSides s.stackDir[t.topDepth]! [item] [] } :: rest }) := by
  obtain ⟨g, tern, items, sv, sd, nei, fo, tstack, tb, tsl⟩ := s
  subst hs; rfl

/-- Closing the top entry `t` into `item` records a 2-attached connected piece with terminals
`makeVs t.vStart t.topDepth`. `item` is an allocated node that is not yet attached anywhere (no
parent, in no entry's spans), and `t`'s items all lie on the side `stackDir[t.topDepth]`. -/
theorem finishTstackTop_complete (s : WalkState) (item : ItemId) (t : TEntry) (rest : List TEntry)
    (hs : s.tstack = t :: rest) (h : s.Inv)
    (hitem : item < s.items.size) (hnode : 1 + s.g.nv + s.g.ne ≤ item)
    (hroot : ∀ p, ¬ Items.IsParent s.items p item)
    (hfree : ∀ t' ∈ s.tstack, item ∉ t'.spans.1 ++ t'.spans.2)
    (hside : getSide t.spans (!s.stackDir[t.topDepth]!) = []) :
    ((finishTstackTop item).run s).2.Inv := by
  rw [finishTstackTop_run_eq s item t rest hs]
  set dir := s.stackDir[t.topDepth]!
  set f : Item → Item := fun it =>
    { it with vs := setSides dir (some s.stackVerts[t.topDepth]!) (some t.vStart),
              ch := getSide t.spans dir }
  have ht := h.entries t (by simp [hs])
  have hitem_edge : ∀ e, e < s.g.ne → edgeItem s.g e ≠ item := by
    intro e he; show 1 + s.g.nv + e ≠ item; omega
  have hnb : ∀ i, i ≠ item → ¬ Items.Below s.items i item := fun i hi hb =>
    hi (Items.Below.eq_of_no_parent hroot hb)
  -- the closed item's subtree is exactly `t`'s edge set
  have hsub : ∀ e, e < s.g.ne →
      (Items.EdgeBelow s.g (s.items.modify item f) item e ↔ t.edges s.g s.items e) := by
    intro e he
    have hch : Items.ch (s.items.modify item f) item = getSide t.spans dir := by
      rw [Items.ch_modify_at item f hitem]
    have hne : ∀ c ∈ getSide t.spans dir, c ≠ item := fun c hc hce => by
      subst hce; exact hfree t (by simp [hs]) ((mem_of_getSide_nil dir t.spans hside c).2 hc)
    simp only [TEntry.edges, Items.EdgeBelow, mem_of_getSide_nil dir t.spans hside]
    constructor
    · intro hb
      rcases hb.head_cases with heq | ⟨c, hc, hb⟩
      · exact absurd heq.symm (hitem_edge e he)
      · simp only [Items.IsParent, hch] at hc
        exact ⟨c, hc, (Items.Below_modify_of_not_below item f (hnb c (hne c hc))).1 hb⟩
    · rintro ⟨c, hc, hb⟩
      exact .head (by simpa [Items.IsParent, hch] using hc)
        ((Items.Below_modify_of_not_below item f (hnb c (hne c hc))).2 hb)
  refine ⟨?_, ?_⟩
  · intro t' ht'
    simp only [List.mem_cons] at ht'
    rcases ht' with rfl | ht'
    · have hE : ∀ e, e < s.g.ne →
          (TEntry.edges s.g (s.items.modify item f) { t with spans := setSides dir [item] [] } e ↔
            t.edges s.g s.items e) := by
        intro e he
        rw [← hsub e he]
        simp only [TEntry.edges]
        constructor
        · rintro ⟨i, hi, hb⟩
          rwa [List.mem_singleton.1 ((mem_setSides dir [item] i).1 hi)] at hb
        · exact fun hb => ⟨item, (mem_setSides dir [item] item).2 (List.mem_singleton.2 rfl), hb⟩
      exact ⟨(Graph.ConnEdges.congr hE).2 ht.conn, (Graph.TwoAttached.congr hE).2 ht.attached⟩
    · exact WalkState.EntryInv.congr (s := s) rfl rfl
        (fun e _ => TEntry.edges_modify_of_not_mem item f hroot (hfree t' (by simp [hs, ht'])) e)
        (h.entries t' (by simp [hs, ht']))
  · intro i hi hsize
    simp only [Array.size_modify] at hsize
    by_cases hi' : i = item
    · subst hi'
      refine ⟨(Graph.ConnEdges.congr hsub).2 ht.conn, ?_⟩
      intro u v huv
      rw [Items.vs_modify_at i f hitem] at huv
      rcases setSides_eq _ _ _ _ _ huv with ⟨hu, hv⟩ | ⟨hu, hv⟩ <;>
        simp only [Option.some.injEq] at hu hv <;> subst hu hv
      · exact (Graph.TwoAttached.congr hsub).2 ht.attached.comm
      · exact (Graph.TwoAttached.congr hsub).2 ht.attached
    · exact WalkState.ItemInv.congr (s := s) rfl (Items.vs_modify_of_ne item f hi')
        (fun _ _ => Items.Below_modify_of_not_below item f (hnb i hi')) (h.nodes i hi hsize)

namespace Graph
variable {g : Graph} {E : Nat → Prop}

theorem inc_of_pairEq {e a b : Nat} (hab : Items.PairEq (a, b) g.edges[e]!) : g.Inc e a ∧ g.Inc e b := by
  rcases hab with hab | hab
  · exact ⟨.inl (congrArg Prod.fst hab).symm, .inr (congrArg Prod.snd hab).symm⟩
  · exact ⟨.inr (congrArg Prod.fst hab).symm, .inl (congrArg Prod.snd hab).symm⟩

theorem TwoAttached.single_of_pairEq {e₀ a b : Nat} (h : ∀ e, E e ↔ e = e₀)
    (hab : Items.PairEq (a, b) g.edges[e₀]!) : g.TwoAttached E a b := by
  intro v e _ _ _ hE _ hv _
  rw [(h e).1 hE] at hv
  rcases hab with hab | hab
  · have h1 : a = (g.edges[e₀]!).1 := congrArg Prod.fst hab
    have h2 : b = (g.edges[e₀]!).2 := congrArg Prod.snd hab
    subst h1 h2
    rcases hv with rfl | rfl <;> simp
  · have h1 : a = (g.edges[e₀]!).2 := congrArg Prod.fst hab
    have h2 : b = (g.edges[e₀]!).1 := congrArg Prod.snd hab
    subst h1 h2
    rcases hv with rfl | rfl <;> simp
end Graph

theorem edgeItem_inj {g : Graph} {e e' : Nat} (h : edgeItem g e = edgeItem g e') : e = e' := by
  have : 1 + g.nv + e = 1 + g.nv + e' := h
  omega

theorem getSide_setSides_not {α : Type} (dir : Bool) (a : List α) :
    getSide (setSides dir a []) (!dir) = [] := by cases dir <;> rfl

theorem getSide_merge_nil {α : Type} (dir : Bool) (p q : List α × List α)
    (hp : getSide p dir = []) (hq : getSide q dir = []) :
    getSide (p.1 ++ q.1, q.2 ++ p.2) dir = [] := by cases dir <;> simp_all [getSide]

theorem mem_merge {α : Type} (p q : List α × List α) (i : α) :
    i ∈ (p.1 ++ q.1, q.2 ++ p.2).1 ++ (p.1 ++ q.1, q.2 ++ p.2).2 ↔ i ∈ p.1 ++ p.2 ∨ i ∈ q.1 ++ q.2 := by
  simp only [List.mem_append]
  constructor <;> intro h <;> rcases h with (h | h) | h | h <;> simp [h]

namespace Items
variable {items : Items}

theorem ch_modify_ch_eq (j : ItemId) (f : Item → Item) (hf : ∀ it, (f it).ch = it.ch) (p : ItemId) :
    Items.ch (items.modify j f) p = items.ch p := by
  by_cases h : p = j
  · subst h
    by_cases hj : p < items.size
    · rw [ch_modify_at p f hj, hf]; simp [ch, Array.getElem?_eq_getElem hj]
    · simp [ch, Array.getElem?_modify, Array.getElem?_eq_none (Nat.le_of_not_lt hj)]
  · exact ch_modify_of_ne j f h

theorem getElem!_modify_type (j : ItemId) (f : Item → Item) (hf : ∀ it, (f it).type = it.type)
    {p : ItemId} (hp : p < items.size) : (items.modify j f)[p]!.type = items.type p := by
  have hp' : p < (items.modify j f).size := by simpa using hp
  rw [getElem!_pos (items.modify j f) p hp', Array.getElem_modify]
  simp only [type, Array.getElem?_eq_getElem hp, Option.map_some, Option.getD_some]
  split <;> simp [hf]

theorem getElem!_modify_ch (j : ItemId) (f : Item → Item) (hf : ∀ it, (f it).ch = it.ch)
    {p : ItemId} (hp : p < items.size) : (items.modify j f)[p]!.ch = items.ch p := by
  have hp' : p < (items.modify j f).size := by simpa using hp
  rw [getElem!_pos (items.modify j f) p hp', Array.getElem_modify]
  simp only [ch, Array.getElem?_eq_getElem hp, Option.map_some, Option.getD_some]
  split <;> simp [hf]

theorem vs_push_of_ne (x : Item) {p : ItemId} (h : p ≠ items.size) :
    Items.vs (items.push x) p = items.vs p := by
  simp [vs, Array.getElem?_push, h]

theorem ch_push_size (x : Item) (hx : x.ch = []) : Items.ch (items.push x) items.size = [] := by
  rw [ch_push_nil x hx]; simp [ch]

theorem type_modify_type_eq (j : ItemId) (f : Item → Item) (hf : ∀ it, (f it).type = it.type) (p : ItemId) :
    Items.type (items.modify j f) p = items.type p := by
  simp only [type, Array.getElem?_modify]
  split
  · cases items[p]? <;> simp [hf]
  · rfl

theorem type_eq_getElem! {i : ItemId} (hi : i < items.size) : items.type i = items[i]!.type := by
  simp [type, getElem!_pos items i hi, Array.getElem?_eq_getElem hi]

theorem ch_eq_getElem! {i : ItemId} (hi : i < items.size) : items.ch i = items[i]!.ch := by
  simp [ch, getElem!_pos items i hi, Array.getElem?_eq_getElem hi]

theorem Tree.modify {g : Graph} (ht : Items.Tree g items) (j : ItemId) (f : Item → Item)
    (hty : ∀ it, (f it).type = it.type) (hch : ∀ it, (f it).ch = it.ch) : Items.Tree g (items.modify j f) := by
  have hC : ∀ p, Items.ch (items.modify j f) p = items.ch p := ch_modify_ch_eq j f hch
  have hT : ∀ p, Items.type (items.modify j f) p = items.type p := type_modify_type_eq j f hty
  have hP : Items.IsParent (items.modify j f) = items.IsParent :=
    funext fun _ => funext fun _ => propext (IsParent_congr hC)
  have hB : Items.Below (items.modify j f) = items.Below :=
    funext fun _ => funext fun _ => propext (Below_congr hC)
  obtain ⟨size, root, vert, edge, node, ch_lt, uniq, rnp, nodup, reach, vch, rch⟩ := ht
  exact ⟨by simpa using size, by rw [hT]; exact root, fun v hv => by rw [hT]; exact vert v hv,
    fun e he => by rw [hT]; exact edge e he, fun i hi hs => by rw [hT]; exact node i hi (by simpa using hs),
    by rw [hP]; simpa using ch_lt, by rw [hP]; simpa using uniq, by rw [hP]; exact rnp,
    fun p => by rw [hC]; exact nodup p, fun i hi => by rw [hB]; exact reach i (by simpa using hi),
    fun v c h => by rw [hT]; exact vch v c (by rwa [hP] at h),
    fun c h => by rw [hT]; exact rch c (by rwa [hP] at h)⟩

theorem node_of_type_P {g : Graph} (ht : Items.Tree g items) {i : ItemId} (hi : i < items.size)
    (hP : items.type i = .P) : 1 + g.nv + g.ne ≤ i := by
  refine Decidable.byContradiction fun hlt => ?_
  rcases Nat.lt_or_ge i (1 + g.nv) with h1 | h1
  · rcases Nat.eq_zero_or_pos i with rfl | h0
    · exact absurd (ht.root.symm.trans hP) (by decide)
    · have hv := ht.vert (i - 1) (by omega)
      have : vertItem (i - 1) = i := show (1 + (i - 1) : Nat) = i by omega
      rw [this, hP] at hv; exact absurd hv (by decide)
  · have hv := ht.edge (i - (1 + g.nv)) (by omega)
    have : edgeItem g (i - (1 + g.nv)) = i := show (1 + g.nv + (i - (1 + g.nv)) : Nat) = i by omega
    rw [this, hP] at hv; exact absurd hv (by decide)

end Items

namespace TEntry
variable {g : Graph} {items : Items}

theorem edges_edgeEntry (dir : Bool) (vStart topDepth fi e : Nat)
    (hq : Items.ch items (edgeItem g e) = []) (e' : Nat) :
    TEntry.edges g items ⟨vStart, topDepth, fi, setSides dir [edgeItem g e] []⟩ e' ↔ e' = e := by
  simp only [edges, mem_setSides, List.mem_singleton, exists_eq_left, Items.EdgeBelow]
  constructor
  · intro hb
    rcases hb.head_cases with heq | ⟨c, hc, _⟩
    · exact (edgeItem_inj heq).symm
    · simp [Items.IsParent, hq] at hc
  · rintro rfl; exact .refl

theorem edges_single {t : TEntry} (dir : Bool) (i : ItemId) (hside : getSide t.spans (!dir) = [])
    (hi : getSide t.spans dir = [i]) (e : Nat) : t.edges g items e ↔ items.EdgeBelow g i e := by
  simp only [edges, mem_of_getSide_nil dir t.spans hside, hi, List.mem_singleton, exists_eq_left]

theorem edges_unwrap (dir : Bool) (vStart topDepth fi : Nat) (i : ItemId)
    (hi : ∀ e, e < g.ne → edgeItem g e ≠ i) {e : Nat} (he : e < g.ne) :
    TEntry.edges g items ⟨vStart, topDepth, fi, setSides dir (items.ch i) []⟩ e ↔ items.EdgeBelow g i e := by
  simp only [edges, mem_setSides, Items.EdgeBelow]
  constructor
  · rintro ⟨c, hc, hb⟩; exact .head hc hb
  · intro hb
    rcases hb.head_cases with heq | ⟨c, hc, hb⟩
    · exact absurd heq.symm (hi e he)
    · exact ⟨c, hc, hb⟩

end TEntry

section Run
variable {α β : Type}

theorem WalkM.run_bind (x : WalkM α) (f : α → WalkM β) (s : WalkState) :
    (x >>= f).run s = (f (x.run s).1).run (x.run s).2 := rfl
theorem List.head!_cons' {γ : Type} [Inhabited γ] (a : γ) (l : List γ) : (a :: l).head! = a := rfl
theorem WalkM.get_run (s : WalkState) : (get : WalkM WalkState).run s = (s, s) := rfl
theorem WalkM.pure_run (a : α) (s : WalkState) : (pure a : WalkM α).run s = (a, s) := rfl
theorem WalkM.modify_run (f : WalkState → WalkState) (s : WalkState) :
    (modify f : WalkM Unit).run s = ((), f s) := rfl
theorem WalkM.stackDir_run (d : Nat) (s : WalkState) : (stackDir d).run s = (s.stackDir[d]!, s) := rfl
theorem WalkM.nxt_run (s : WalkState) : nxt.run s = (s.tstack.tail.head!, s) := rfl
theorem WalkM.tstackSize_run (s : WalkState) : tstackSize.run s = (s.tstack.length, s) := rfl
theorem WalkM.getItem_run (i : ItemId) (s : WalkState) : (getItem i).run s = (s.items[i]!, s) := rfl
theorem WalkM.modifyItem_run (i : ItemId) (f : Item → Item) (s : WalkState) :
    (modifyItem i f).run s = ((), { s with items := s.items.modify i f }) := rfl
theorem WalkM.allocItem_run (ty : NodeType) (s : WalkState) :
    (allocItem ty).run s = (s.items.size, { s with items := s.items.push ⟨ty, (none, none), []⟩ }) := rfl
theorem WalkM.pushTstack_run (vStart topDepth : Nat) (item : ItemId) (s : WalkState) :
    (pushTstack vStart topDepth item).run s = ((), { s with
      tstack := ⟨vStart, topDepth, s.nxtEdgeIdx, setSides s.stackDir[topDepth]! [item] []⟩ :: s.tstack }) := rfl
theorem WalkM.pushEdgeTstack_run (vStart topDepth e : Nat) (s : WalkState) :
    (pushEdgeTstack vStart topDepth e).run s = ((), { s with
      tstack := ⟨vStart, topDepth, s.nxtEdgeIdx, setSides s.stackDir[topDepth]! [edgeItem s.g e] []⟩ :: s.tstack }) := rfl

theorem maybeUnwrapNxt_P_run (s : WalkState) (a b : TEntry) (rest : List TEntry)
    (hs : s.tstack = a :: b :: rest) :
    (maybeUnwrapNxt .P).run s =
      if s.ternarize then (allocItem .P).run s
      else
        let topDir := s.stackDir[b.topDepth]!
        let item := (getSide b.spans topDir).head!
        if s.items[item]!.type == .P then
          (item, { s with tstack := a :: { b with spans := setSides topDir s.items[item]!.ch [] } :: rest })
        else (allocItem .P).run s := by
  obtain ⟨g, tern, items, sv, sd, nei, fo, tstack, tb, tsl⟩ := s
  subst hs
  cases tern
  · simp only [maybeUnwrapNxt, WalkM.run_bind, WalkM.get_run, WalkM.nxt_run, WalkM.stackDir_run,
      WalkM.getItem_run, modifyNxt, Bool.or_false, show (NodeType.P == NodeType.R) = false from rfl,
      Bool.false_eq_true, ↓reduceIte, List.tail_cons, List.head!_cons']
    by_cases hh : (items[(getSide b.spans sd[b.topDepth]!).head!]!.type == NodeType.P) = true
    · simp only [hh, ↓reduceIte]; rfl
    · simp only [hh, Bool.false_eq_true, ↓reduceIte]
  · simp only [maybeUnwrapNxt, WalkM.run_bind, WalkM.get_run, Bool.or_true, ↓reduceIte]

end Run

namespace WalkState

/-- Allocating a fresh node item preserves the invariant: it has no edges below it yet. -/
theorem Inv.alloc (s : WalkState) (ty : NodeType) (h : s.Inv) (hsize : 1 + s.g.nv + s.g.ne ≤ s.items.size) :
    WalkState.Inv { s with items := s.items.push ⟨ty, (none, none), []⟩ } := by
  have hB : ∀ a i, Items.Below (s.items.push ⟨ty, (none, none), []⟩) a i ↔ Items.Below s.items a i :=
    fun _ _ => Items.Below_push_nil _ rfl
  refine ⟨fun t ht => EntryInv.congr (s := s) rfl rfl (fun e _ => TEntry.edges_congr (fun i _ e => hB i _) e)
    (h.entries t ht), fun i hi hsz => ?_⟩
  have hsz' : i < s.items.size + 1 := by simpa using hsz
  by_cases hi' : i = s.items.size
  · subst hi'
    have hE : ∀ e, e < s.g.ne → ¬ Items.EdgeBelow s.g (s.items.push ⟨ty, (none, none), []⟩) s.items.size e := by
      intro e he hb
      rcases hb.head_cases with heq | ⟨c, hc, _⟩
      · have : s.items.size = 1 + s.g.nv + e := heq
        omega
      · simp [Items.IsParent, Items.ch_push_size ⟨ty, (none, none), []⟩ rfl] at hc
    exact ⟨Graph.ConnEdges.empty hE, fun _ _ _ => Graph.TwoAttached.empty hE⟩
  · exact ItemInv.congr (s := s) rfl (Items.vs_push_of_ne _ hi') (fun e _ => hB i _)
      (h.nodes i hi (by omega))

/-- Pushing the vertex item of `v` as a new entry at depth `d` with `stackVerts[d] = v`. -/
theorem Inv.pushVert (s : WalkState) (v d : Nat) (h : s.Inv) (hd : s.stackVerts[d]! = v)
    (hc : s.g.ConnEdges (Items.EdgeBelow s.g s.items (vertItem v)))
    (ha : s.g.TwoAttached (Items.EdgeBelow s.g s.items (vertItem v)) v v) :
    ((pushVertTstack v d).run s).2.Inv := by
  rw [pushVertTstack, WalkM.pushTstack_run]
  refine ⟨fun t ht => ?_, fun i hi hsz => ItemInv.congr (s := s) rfl rfl (fun _ _ => Iff.rfl) (h.nodes i hi hsz)⟩
  simp only [List.mem_cons] at ht
  rcases ht with rfl | ht
  · have hE : ∀ e, TEntry.edges s.g s.items ⟨v, d, s.nxtEdgeIdx, setSides s.stackDir[d]! [vertItem v] []⟩ e ↔
        Items.EdgeBelow s.g s.items (vertItem v) e := by
      intro e; simp only [TEntry.edges, mem_setSides, List.mem_singleton, exists_eq_left]
    refine ⟨(Graph.ConnEdges.congr fun e _ => hE e).2 hc, ?_⟩
    show s.g.TwoAttached _ v s.stackVerts[d]!
    rw [hd]; exact (Graph.TwoAttached.congr fun e _ => hE e).2 ha
  · exact EntryInv.congr (s := s) rfl rfl (fun _ _ => Iff.rfl) (h.entries t ht)

/-- Pushing the back edge `e = {curV, stackVerts[lowval]}` as a new entry, after recording its `vs`. -/
theorem Inv.pushEdge {s s' : WalkState} (curV lowval e : Nat) (vsv : Option Nat × Option Nat) (h : s.Inv)
    (hs' : s' = { s with
      items := s.items.modify (edgeItem s.g e) fun it => { it with vs := vsv },
      nxtEdgeIdx := s.nxtEdgeIdx + 1,
      firstOccurrence := s.firstOccurrence.modify lowval (min · s.nxtEdgeIdx),
      tstack := ⟨curV, lowval, s.nxtEdgeIdx, setSides s.stackDir[lowval]! [edgeItem s.g e] []⟩ :: s.tstack })
    (he : e < s.g.ne) (hq : Items.ch s.items (edgeItem s.g e) = [])
    (hend : Items.PairEq (curV, s.stackVerts[lowval]!) s.g.edges[e]!) : s'.Inv := by
  subst hs'
  have hB : ∀ a i, Items.Below (s.items.modify (edgeItem s.g e) fun it => { it with vs := vsv }) a i ↔
      Items.Below s.items a i :=
    fun _ _ => Items.Below_modify_ch_eq (edgeItem s.g e) (fun it => { it with vs := vsv }) fun _ => rfl
  have hq' : Items.ch (s.items.modify (edgeItem s.g e) fun it => { it with vs := vsv }) (edgeItem s.g e) = [] := by
    rw [Items.ch_modify_ch_eq (edgeItem s.g e) (fun it => { it with vs := vsv }) (fun _ => rfl), hq]
  refine ⟨fun t ht => ?_, fun i hi hsz => ?_⟩
  · simp only [List.mem_cons] at ht
    rcases ht with rfl | ht
    · have hE : ∀ e', TEntry.edges s.g (s.items.modify (edgeItem s.g e) fun it => { it with vs := vsv })
          ⟨curV, lowval, s.nxtEdgeIdx, setSides s.stackDir[lowval]! [edgeItem s.g e] []⟩ e' ↔ e' = e :=
        TEntry.edges_edgeEntry s.stackDir[lowval]! curV lowval s.nxtEdgeIdx e hq'
      exact ⟨Graph.ConnEdges.single hE, Graph.TwoAttached.single_of_pairEq hE hend⟩
    · exact EntryInv.congr (s := s) rfl rfl (fun e _ => TEntry.edges_congr (fun i _ e => hB i _) e) (h.entries t ht)
  · have hi' : 1 + s.g.nv + s.g.ne ≤ i := hi
    have hsz' : i < s.items.size := by simpa using hsz
    have hne : i ≠ edgeItem s.g e := by show i ≠ 1 + s.g.nv + e; omega
    exact ItemInv.congr (s := s) rfl (Items.vs_modify_of_ne _ _ hne) (fun e _ => hB i _) (h.nodes i hi' hsz')

/-- One walk step that keeps the invariant and does not disturb the graph, the vertex stack,
or the item tree shape. -/
structure Step (v : Nat) (s s' : WalkState) : Prop where
  inv : s'.Inv
  g : s'.g = s.g
  sv : s'.stackVerts = s.stackVerts
  below : ∀ i, Items.Below s'.items (vertItem v) i ↔ Items.Below s.items (vertItem v) i

theorem Step.trans {v : Nat} {s₁ s₂ s₃ : WalkState} (h₁ : Step v s₁ s₂) (h₂ : Step v s₂ s₃) : Step v s₁ s₃ :=
  ⟨h₂.inv, h₂.g.trans h₁.g, h₂.sv.trans h₁.sv, fun i => (h₂.below i).trans (h₁.below i)⟩

theorem Inv.pushVert_of_step {s S : WalkState} (curV d : Nat) (st : Step curV s S) (hd : s.stackVerts[d]! = curV)
    (hvc : s.g.ConnEdges (Items.EdgeBelow s.g s.items (vertItem curV)))
    (hva : s.g.TwoAttached (Items.EdgeBelow s.g s.items (vertItem curV)) curV curV) :
    ((pushVertTstack curV d).run S).2.Inv := by
  refine Inv.pushVert S curV d st.inv (by rw [st.sv, hd]) ?_ ?_
  · rw [st.g]; exact (Graph.ConnEdges.congr fun e _ => st.below _).2 hvc
  · rw [st.g]; exact (Graph.TwoAttached.congr fun e _ => st.below _).2 hva

/-- Shape of the entry `b` under the freshly pushed back edge `{curV, stackVerts[lv]}` when the
type-1 P-check fires: `b` is a single root item touching `curV`, on the side `stackDir[lv]`, and
shares no item with the entries below it. -/
structure BackCheckOk (s : WalkState) (curV lv : Nat) (b : TEntry) (rest : List TEntry) : Prop where
  touch : s.g.Touches (b.edges s.g s.items) curV
  side : getSide b.spans (!s.stackDir[lv]!) = []
  single : ∃ i, getSide b.spans s.stackDir[lv]! = [i]
  root : ∀ i ∈ b.spans.1 ++ b.spans.2, ∀ p, ¬ Items.IsParent s.items p i
  disj : ∀ i ∈ b.spans.1 ++ b.spans.2, ∀ t ∈ rest, i ∉ t.spans.1 ++ t.spans.2

/-- `mergeTstackTops` followed by `finishTstackTop item` when `cur` is the single back edge `e` and
`nxt = b` has the same terminals. -/
theorem Step.mergeFinish (s : WalkState) (curV lv fi e item : Nat) (b : TEntry) (rest : List TEntry)
    (hs : s.tstack = ⟨curV, lv, fi, setSides s.stackDir[lv]! [edgeItem s.g e] []⟩ :: b :: rest)
    (h : s.Inv) (hcur : curV < s.g.nv) (he : e < s.g.ne) (hq : Items.ch s.items (edgeItem s.g e) = [])
    (hend : Items.PairEq (curV, s.stackVerts[lv]!) s.g.edges[e]!)
    (hb1 : b.vStart = curV) (hb2 : b.topDepth = lv)
    (htouch : s.g.Touches (b.edges s.g s.items) curV)
    (hside : getSide b.spans (!s.stackDir[lv]!) = [])
    (hitem : item < s.items.size) (hnode : 1 + s.g.nv + s.g.ne ≤ item)
    (hroot : ∀ p, ¬ Items.IsParent s.items p item)
    (hfree : ∀ t ∈ s.tstack, item ∉ t.spans.1 ++ t.spans.2) :
    Step curV s ((finishTstackTop item).run (mergeTstackTops.run s).2).2 := by
  set q : TEntry := ⟨curV, lv, fi, setSides s.stackDir[lv]! [edgeItem s.g e] []⟩ with hqdef
  have hE : ∀ e', q.edges s.g s.items e' ↔ e' = e := TEntry.edges_edgeEntry _ _ _ _ _ hq
  obtain ⟨b', hb'⟩ : ∃ b' : TEntry,
    b' = ⟨b.vStart, min b.topDepth q.topDepth, b.firstIdx, (q.spans.1 ++ b.spans.1, b.spans.2 ++ q.spans.2)⟩ :=
    ⟨_, rfl⟩
  obtain ⟨s₂, hs₂⟩ : ∃ s₂ : WalkState, s₂ = ({ s with tstack := b' :: rest } : WalkState) := ⟨_, rfl⟩
  have hrun : mergeTstackTops.run s = ((), s₂) := by
    rw [mergeTstackTops_run_eq s q b rest hs, hs₂, hb']
  rw [hrun]
  have hinv₂ : s₂.Inv := by
    subst hs₂
    refine ⟨?_, fun i hi hsz => ItemInv.congr (s := s) rfl rfl (fun _ _ => Iff.rfl) (h.nodes i hi hsz)⟩
    have := mergeTstackTops_sound s q b rest hs h.entries ?_ ?_
    · rwa [mergeTstackTops_run_eq s q b rest hs, ← hb'] at this
    · intro _ _
      exact ⟨curV, ⟨e, he, (hE e).2 rfl, (Graph.inc_of_pairEq hend).1⟩, htouch⟩
    · simp only [hqdef, TEntry.top, hb1, hb2, Nat.min_self, List.mem_cons,
        List.not_mem_nil, or_false]
      rintro v (rfl | rfl | rfl | rfl) <;> simp
  have hside₂ : getSide b'.spans (!s₂.stackDir[b'.topDepth]!) = [] := by
    rw [hs₂, hb']
    show getSide (q.spans.1 ++ b.spans.1, b.spans.2 ++ q.spans.2) (!s.stackDir[min b.topDepth lv]!) = []
    rw [hb2, Nat.min_self]
    exact getSide_merge_nil _ _ _ (getSide_setSides_not _ _) hside
  have hfree₂ : ∀ t ∈ s₂.tstack, item ∉ t.spans.1 ++ t.spans.2 := by
    intro t ht hmem
    rw [hs₂] at ht
    simp only [List.mem_cons] at ht
    rcases ht with rfl | ht
    · rw [hb'] at hmem
      rcases (mem_merge q.spans b.spans item).1 hmem with hm | hm
      · exact hfree q (by simp [hs]) hm
      · exact hfree b (by simp [hs]) hm
    · exact hfree t (by simp [hs, ht]) hmem
  have hst : s₂.tstack = b' :: rest := by rw [hs₂]
  refine ⟨finishTstackTop_complete s₂ item b' rest hst hinv₂ (by rw [hs₂]; exact hitem) (by rw [hs₂]; exact hnode)
    (by rw [hs₂]; exact hroot) hfree₂ hside₂, ?_, ?_, ?_⟩
  · rw [finishTstackTop_run_eq s₂ item b' rest hst, hs₂]
  · rw [finishTstackTop_run_eq s₂ item b' rest hst, hs₂]
  · intro i
    rw [finishTstackTop_run_eq s₂ item b' rest hst, hs₂]
    refine Items.Below_modify_of_not_below _ _ fun hb => ?_
    have := hb.eq_of_no_parent hroot
    show False
    have : 1 + curV = item := this
    omega

theorem Step.backCheck (s : WalkState) (curV lv fi e : Nat) (b : TEntry) (rest : List TEntry)
    (hs : s.tstack = ⟨curV, lv, fi, setSides s.stackDir[lv]! [edgeItem s.g e] []⟩ :: b :: rest)
    (h : s.Inv) (ht : Items.Tree s.g s.items) (hcur : curV < s.g.nv) (he : e < s.g.ne)
    (hq : Items.ch s.items (edgeItem s.g e) = [])
    (hend : Items.PairEq (curV, s.stackVerts[lv]!) s.g.edges[e]!)
    (hb1 : b.vStart = curV) (hb2 : b.topDepth = lv)
    (hspan : ∀ t ∈ s.tstack, ∀ i ∈ t.spans.1 ++ t.spans.2, i < s.items.size)
    (hok : BackCheckOk s curV lv b rest) :
    Step curV s ((finishTstackTop ((maybeUnwrapNxt .P).run s).1).run
      (mergeTstackTops.run ((maybeUnwrapNxt .P).run s).2).2).2 := by
  set q : TEntry := ⟨curV, lv, fi, setSides s.stackDir[lv]! [edgeItem s.g e] []⟩ with hqdef
  obtain ⟨i, hi⟩ := hok.single
  have hib : i ∈ b.spans.1 ++ b.spans.2 := by
    rw [mem_of_getSide_nil _ _ hok.side, hi]; exact List.mem_singleton_self i
  have hilt : i < s.items.size := hspan b (by simp [hs]) i hib
  have hbE : ∀ e', b.edges s.g s.items e' ↔ Items.EdgeBelow s.g s.items i e' :=
    TEntry.edges_single _ i hok.side hi
  have alloc : Step curV s ((finishTstackTop ((allocItem .P).run s).1).run
      (mergeTstackTops.run ((allocItem .P).run s).2).2).2 := by
    rw [WalkM.allocItem_run]
    obtain ⟨s₁, hs₁⟩ : ∃ s₁ : WalkState,
      s₁ = ({ s with items := s.items.push ⟨.P, (none, none), []⟩ } : WalkState) := ⟨_, rfl⟩
    rw [← hs₁]
    have hB : ∀ a j, Items.Below s₁.items a j ↔ Items.Below s.items a j :=
      fun _ _ => by rw [hs₁]; exact Items.Below_push_nil _ rfl
    have hch : ∀ p, Items.ch s₁.items p = Items.ch s.items p := fun p => by
      rw [hs₁]; exact Items.ch_push_nil ⟨.P, (none, none), []⟩ rfl p
    have st : Step curV s s₁ := ⟨by rw [hs₁]; exact Inv.alloc s .P h ht.size, by rw [hs₁], by rw [hs₁], hB _⟩
    have hsz : s₁.items.size = s.items.size + 1 := by rw [hs₁]; simp
    refine st.trans (Step.mergeFinish s₁ curV lv fi e s.items.size b rest (by rw [hs₁]; exact hs) st.inv
      (by rw [hs₁]; exact hcur) (by rw [hs₁]; exact he) (by rw [hs₁]; exact (Items.ch_push_nil ⟨.P, (none, none), []⟩ rfl _).trans hq)
      (by rw [hs₁]; exact hend) hb1 hb2 ?_
      (by rw [hs₁]; exact hok.side) (by omega) (by rw [hs₁]; exact ht.size) ?_ ?_)
    · obtain ⟨e', he', hb, hinc⟩ := hok.touch
      rw [hs₁]
      exact ⟨e', he', (TEntry.edges_congr (fun j _ e => Items.Below_push_nil _ rfl) e').2 hb, hinc⟩
    · intro p hp
      exact Nat.lt_irrefl _ (ht.ch_lt p s.items.size (by rwa [Items.IsParent, hch] at hp))
    · intro t ht' hmem
      rw [hs₁] at ht'
      exact Nat.lt_irrefl _ (hspan t ht' _ hmem)
  rw [maybeUnwrapNxt_P_run s q b rest hs]
  by_cases htern : s.ternarize = true
  · rw [ite_eq_left htern]; exact alloc
  · rw [ite_eq_right htern]
    dsimp only
    have hi' : getSide b.spans s.stackDir[b.topDepth]! = [i] := by rw [hb2]; exact hi
    rw [hi', List.head!_cons']
    by_cases hP : (s.items[i]!.type == NodeType.P) = true
    · rw [ite_eq_left hP]
      dsimp only
      have hiP : Items.type s.items i = .P := by
        rw [Items.type_eq_getElem! hilt]; simpa using hP
      have hinode := Items.node_of_type_P ht hilt hiP
      have hie : ∀ e', e' < s.g.ne → edgeItem s.g e' ≠ i := fun e' he' heq => by
        have := ht.edge e' he'; rw [heq, hiP] at this; exact absurd this (by decide)
      obtain ⟨b', hb'⟩ : ∃ b' : TEntry,
        b' = ({ b with spans := setSides s.stackDir[b.topDepth]! s.items[i]!.ch [] } : TEntry) := ⟨_, rfl⟩
      rw [← hb']
      obtain ⟨s₁, hs₁⟩ : ∃ s₁ : WalkState, s₁ = ({ s with tstack := q :: b' :: rest } : WalkState) := ⟨_, rfl⟩
      rw [← hs₁]
      have hb'E : ∀ e', e' < s.g.ne → (b'.edges s.g s.items e' ↔ b.edges s.g s.items e') := by
        intro e' he'
        rw [hbE, hb', ← Items.ch_eq_getElem! hilt]
        exact TEntry.edges_unwrap _ _ _ _ i hie he'
      have hinv₁ : s₁.Inv := by
        subst hs₁
        refine ⟨fun t ht' => ?_, fun j hj hsz => ItemInv.congr (s := s) rfl rfl (fun _ _ => Iff.rfl) (h.nodes j hj hsz)⟩
        simp only [List.mem_cons] at ht'
        rcases ht' with ht' | ht' | ht'
        · rw [ht']; exact EntryInv.congr (s := s) rfl rfl (fun _ _ => Iff.rfl) (h.entries _ (by simp [hs]))
        · rw [ht']
          obtain ⟨conn, att⟩ := h.entries b (by simp [hs])
          refine ⟨(Graph.ConnEdges.congr hb'E).2 conn, ?_⟩
          have h1 : b'.vStart = b.vStart := by rw [hb']
          have h2 : b'.topDepth = b.topDepth := by rw [hb']
          show s.g.TwoAttached (b'.edges s.g s.items) b'.vStart s.stackVerts[b'.topDepth]!
          rw [h1, h2]
          exact (Graph.TwoAttached.congr hb'E).2 att
        · exact EntryInv.congr (s := s) rfl rfl (fun _ _ => Iff.rfl) (h.entries t (by simp [hs, ht']))
      have st : Step curV s s₁ := ⟨hinv₁, by rw [hs₁], by rw [hs₁], fun _ => by rw [hs₁]⟩
      have hb'1 : b'.vStart = curV := by rw [hb']; exact hb1
      have hb'2 : b'.topDepth = lv := by rw [hb']; exact hb2
      have hside' : getSide b'.spans (!s₁.stackDir[lv]!) = [] := by
        rw [hb', hs₁]; show getSide (setSides s.stackDir[b.topDepth]! _ []) (!s.stackDir[lv]!) = []
        rw [hb2]; exact getSide_setSides_not _ _
      refine st.trans (Step.mergeFinish s₁ curV lv fi e i b' rest (by rw [hs₁]) hinv₁ (by rw [hs₁]; exact hcur)
        (by rw [hs₁]; exact he)
        (by rw [hs₁]; exact hq) (by rw [hs₁]; exact hend) hb'1 hb'2 ?_ hside' (by rw [hs₁]; exact hilt)
        (by rw [hs₁]; exact hinode) (by rw [hs₁]; exact hok.root i hib) ?_)
      · obtain ⟨e', he', hb, hinc⟩ := hok.touch
        rw [hs₁]
        exact ⟨e', he', (hb'E e' he').2 hb, hinc⟩
      · intro t ht' hmem
        rw [hs₁] at ht'
        simp only [List.mem_cons] at ht'
        rcases ht' with rfl | rfl | ht'
        · rw [mem_setSides, List.mem_singleton] at hmem
          exact hie e he hmem.symm
        · rw [hb', mem_setSides, ← Items.ch_eq_getElem! hilt] at hmem
          exact hok.root i hib i hmem
        · exact hok.disj i hib t ht' hmem
    · rw [ite_eq_right hP]; exact alloc

/-- `finishEdge` on a back edge `o = (curV, stackVerts[lv])` with `lv < d` preserves the invariant,
given the stack-shape facts the type-1 P-check relies on. -/
theorem finishEdge_back_inv (curV d lv : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool)
    (s : WalkState) (ho : o.cls = .ret lv .backEdge) (hlow : lv < d)
    (h : s.Inv) (ht : Items.Tree s.g s.items) (hcur : curV < s.g.nv) (he : o.e < s.g.ne)
    (hq : Items.ch s.items (edgeItem s.g o.e) = [])
    (hend : Items.PairEq (curV, s.stackVerts[lv]!) s.g.edges[o.e]!)
    (hd : s.stackVerts[d]! = curV)
    (hvc : s.g.ConnEdges (Items.EdgeBelow s.g s.items (vertItem curV)))
    (hva : s.g.TwoAttached (Items.EdgeBelow s.g s.items (vertItem curV)) curV curV)
    (hspan : ∀ t ∈ s.tstack, ∀ i ∈ t.spans.1 ++ t.spans.2, i < s.items.size)
    (hok : ∀ b rest, s.tstack = b :: rest → b.vStart = curV → b.topDepth = lv →
      BackCheckOk s curV lv b rest) :
    ((finishEdge curV d o origTstack hasVert).run s).2.Inv := by
  have hsz : edgeItem s.g o.e < s.items.size := by
    have := ht.size; show 1 + s.g.nv + o.e < s.items.size; omega
  have hB : ∀ a i, Items.Below (s.items.modify (edgeItem s.g o.e) fun it =>
      { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) a i ↔
      Items.Below s.items a i := fun _ _ => Items.Below_modify_ch_eq (edgeItem s.g o.e)
        (fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) fun _ => rfl
  have hC : ∀ p, Items.ch (s.items.modify (edgeItem s.g o.e) fun it =>
      { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) p =
      Items.ch s.items p := Items.ch_modify_ch_eq (edgeItem s.g o.e)
        (fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) fun _ => rfl
  simp only [finishEdge, ho, OutClass.lowval, OutClass.isTree, OutClass.isType1, WalkM.run_bind,
    WalkM.get_run, WalkM.stackDir_run, Bool.false_eq_true, ↓reduceIte, Nat.not_le.mpr hlow,
    Bool.true_and, Bool.not_true, makeVs_run, WalkM.modifyItem_run, WalkM.pushEdgeTstack_run,
    WalkM.modify_run, WalkM.tstackSize_run, WalkM.nxt_run, List.tail_cons, List.length_cons,
    Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq]
  have st₁ := Inv.pushEdge (s := s) curV lv o.e (setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest))
    h rfl he hq hend
  have step₁ : Step curV s _ := ⟨st₁, rfl, rfl, hB _⟩
  by_cases hc : (s.tstack.length + 1 ≥ 2 ∧ s.tstack.head!.vStart = curV) ∧ s.tstack.head!.topDepth = lv
  · simp only [hc, and_self, ite_true]
    obtain ⟨b, rest, hst⟩ : ∃ b rest, s.tstack = b :: rest := by
      cases h' : s.tstack with
      | nil => rw [h'] at hc; simp at hc
      | cons b rest => exact ⟨b, rest, rfl⟩
    rw [hst] at hc
    simp only [List.head!_cons'] at hc
    obtain ⟨⟨-, hb1⟩, hb2⟩ := hc
    have hok' := hok b rest hst hb1 hb2
    have step₂ := Step.backCheck _ curV lv s.nxtEdgeIdx o.e b rest (by simp only [hst]) st₁
      (Items.Tree.modify ht (edgeItem s.g o.e)
        (fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) })
        (fun _ => rfl) (fun _ => rfl)) hcur he (by rw [hC]; exact hq) hend hb1 hb2 ?_ ?_
    · cases hasVert
      · simp only [Bool.not_false, ↓reduceIte, WalkM.run_bind, WalkM.pure_run]
        exact Inv.pushVert_of_step curV d (step₁.trans step₂) hd hvc hva
      · simpa only [Bool.not_true, Bool.false_eq_true, ↓reduceIte, WalkM.run_bind, WalkM.pure_run]
          using step₂.inv
    · intro t ht' i hi
      rw [hst] at ht'
      simp only [List.mem_cons] at ht'
      simp only [Array.size_modify]
      rcases ht' with rfl | ht'
      · rw [mem_setSides, List.mem_singleton] at hi; rw [hi]; exact hsz
      · exact hspan t (by rw [hst]; exact List.mem_cons.2 ht') i hi
    · obtain ⟨touch, side, single, root, disj⟩ := hok'
      refine ⟨?_, side, single, fun i hi p hp => root i hi p ?_, disj⟩
      · obtain ⟨e', he', hb, hinc⟩ := touch
        exact ⟨e', he', (TEntry.edges_congr (fun j _ e => hB j _) e').2 hb, hinc⟩
      · rwa [Items.IsParent, hC] at hp
  · simp only [hc, ite_false]
    cases hasVert
    · simp only [Bool.not_false, ↓reduceIte, WalkM.run_bind, WalkM.pure_run]
      exact Inv.pushVert_of_step curV d step₁ hd hvc hva
    · simpa only [Bool.not_true, Bool.false_eq_true, ↓reduceIte, WalkM.pure_run] using st₁

end WalkState

/-- `finishEdge` preserves the invariant. -/
theorem finishEdge_inv (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool)
    (s : WalkState) (h : s.Inv) :
    ((finishEdge curV d o origTstack hasVert).run s).2.Inv := by
  sorry

/-- The whole walk preserves the invariant. -/
theorem walkTree_inv (t : DfsTree) (d : Nat) (s : WalkState) (h : s.Inv) :
    ((walkTree t d).run s).2.Inv := by
  sorry

/-- Completeness: after the walk, the allocated nodes partition the edges (every edge is below
exactly one child chain from the root) — the `Items.Tree` content of `Items.WF`. -/
theorem walk_nodes_partition (g : Graph) (tern : Bool) (forest : List DfsTree) :
    let s := g.walk tern forest
    ∀ e, e < g.ne → ∃ p, Items.IsParent s.items p (edgeItem g e) ∧
      ∀ p', Items.IsParent s.items p' (edgeItem g e) → p' = p := by
  sorry

end Spqr
