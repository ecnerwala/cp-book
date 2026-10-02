import Mathlib.Tactic.Set
import Spqr.GraphLemmas

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

theorem ch_modify_self (j : ItemId) (f : Item → Item) (hj : j < items.size) :
    Items.ch (items.modify j f) j = (f items[j]).ch := by
  simp [ch, Array.getElem?_modify, Array.getElem?_eq_getElem hj]

theorem vs_modify_of_ne (j : ItemId) (f : Item → Item) {p : ItemId} (h : p ≠ j) :
    Items.vs (items.modify j f) p = items.vs p := by
  simp [vs, Array.getElem?_modify, h.symm]

theorem vs_modify_self (j : ItemId) (f : Item → Item) (hj : j < items.size) :
    Items.vs (items.modify j f) j = (f items[j]).vs := by
  simp [vs, Array.getElem?_modify, Array.getElem?_eq_getElem hj]

theorem Below_modify_ch_eq (j : ItemId) (f : Item → Item) (hf : ∀ it, (f it).ch = it.ch) {a i : ItemId} :
    Items.Below (items.modify j f) a i ↔ items.Below a i := by
  refine Below_congr fun p => ?_
  by_cases h : p = j
  · subst h
    by_cases hj : p < items.size
    · rw [ch_modify_self p f hj, hf]; simp [ch, Array.getElem?_eq_getElem hj]
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

theorem mergeTstackTops_run (s : WalkState) (cur nxt : TEntry) (rest : List TEntry)
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
    simp only [s', mergeTstackTops_run s cur nxt rest hs]
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

theorem finishTstackTop_run (s : WalkState) (item : ItemId) (t : TEntry) (rest : List TEntry)
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
  rw [finishTstackTop_run s item t rest hs]
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
      rw [Items.ch_modify_self item f hitem]
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
      rw [Items.vs_modify_self i f hitem] at huv
      rcases setSides_eq _ _ _ _ _ huv with ⟨hu, hv⟩ | ⟨hu, hv⟩ <;>
        simp only [Option.some.injEq] at hu hv <;> subst hu hv
      · exact (Graph.TwoAttached.congr hsub).2 ht.attached.comm
      · exact (Graph.TwoAttached.congr hsub).2 ht.attached
    · exact WalkState.ItemInv.congr (s := s) rfl (Items.vs_modify_of_ne item f hi')
        (fun _ _ => Items.Below_modify_of_not_below item f (hnb i hi')) (h.nodes i hi hsize)

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
