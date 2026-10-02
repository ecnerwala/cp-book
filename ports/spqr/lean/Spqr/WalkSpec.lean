import Mathlib.Tactic.Set
import Spqr.GraphLemmas
import Spqr.WalkTyping
import Spqr.Frame

/-!
# The walk invariant: soundness and completeness per step

Stated on `walkTree` directly (no ear structure). Every tstack entry `t` describes a connected set
of original edges `t.edges` whose attachment to the rest of the block is confined to its terminals
`t.Term D`: the bottom `t.vStart` and the path vertices `stackVerts[k]` for `t.topDepth ≤ k ≤ D`,
where `D` is the depth the walk is currently at (a type-2 entry keeps attachments at intermediate
chain vertices until they are closed). Every closed item is 2-attached at its `vs`. *Soundness*:
`mergeTstackTops` joins two entries sharing a vertex, so the union is again connected and attached
within the merged terminals. *Completeness*: every item `finishEdge` closes is 2-attached and
connected, i.e. cut off by a genuine 2-separation. The stack-shape facts these steps rely on are
collected in `FinishOk` and are to be discharged on the ear layer (`EarSpec.lean`).
-/

namespace Spqr
open WalkM

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

theorem ch_modify_ch_eq (j : ItemId) (f : Item → Item) (hf : ∀ it, (f it).ch = it.ch) (p : ItemId) :
    Items.ch (items.modify j f) p = items.ch p := by
  by_cases h : p = j
  · subst h
    by_cases hj : p < items.size
    · rw [ch_modify_at p f hj, hf]; simp [ch, Array.getElem?_eq_getElem hj]
    · simp [ch, Array.getElem?_modify, Array.getElem?_eq_none (Nat.le_of_not_lt hj)]
  · exact ch_modify_of_ne j f h

/-- Modifying `j` changes `ch` at `j` only. -/
theorem IsParent_modify {j : ItemId} {f : Item → Item} {p c : ItemId} (h : Items.IsParent (items.modify j f) p c) :
    items.IsParent p c ∨ (p = j ∧ ∃ hj : j < items.size, c ∈ (f items[j]).ch) := by
  by_cases hp : p = j
  · subst hp
    by_cases hj : p < items.size
    · exact .inr ⟨rfl, hj, by rwa [IsParent, ch_modify_at p f hj] at h⟩
    · left
      rw [IsParent, ch, Array.getElem?_modify, Array.getElem?_eq_none (Nat.le_of_not_lt hj)] at h
      simp at h
  · left; rwa [IsParent, ch_modify_of_ne j f hp] at h

theorem vs_push_of_ne (x : Item) {p : ItemId} (h : p ≠ items.size) :
    Items.vs (items.push x) p = items.vs p := by
  simp [vs, Array.getElem?_push, h]

theorem ch_push_size (x : Item) (hx : x.ch = []) : Items.ch (items.push x) items.size = [] := by
  rw [ch_push_nil x hx]; simp [ch]

theorem type_push_of_ne (x : Item) {p : ItemId} (h : p ≠ items.size) :
    Items.type (items.push x) p = items.type p := by
  simp [type, Array.getElem?_push, h]

theorem type_push_size (x : Item) : Items.type (items.push x) items.size = x.type := by
  simp [type]

theorem type_modify_type_eq (j : ItemId) (f : Item → Item) (hf : ∀ it, (f it).type = it.type) (p : ItemId) :
    Items.type (items.modify j f) p = items.type p := by
  simp only [type, Array.getElem?_modify]
  split
  · cases items[p]? <;> simp [hf]
  · rfl

theorem type_eq_getElem {i : ItemId} (hi : i < items.size) : items.type i = items[i].type := by
  simp [type, Array.getElem?_eq_getElem hi]

theorem ch_eq_getElem {i : ItemId} (hi : i < items.size) : items.ch i = items[i].ch := by
  simp [ch, Array.getElem?_eq_getElem hi]

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

theorem getSide_setSides_not {α : Type} (dir : Bool) (a : List α) :
    getSide (setSides dir a []) (!dir) = [] := by cases dir <;> rfl

theorem getSide_merge_nil {α : Type} (dir : Bool) (p q : List α × List α)
    (hp : getSide p dir = []) (hq : getSide q dir = []) :
    getSide (p.1 ++ q.1, q.2 ++ p.2) dir = [] := by cases dir <;> simp_all [getSide]

theorem mem_merge {α : Type} (p q : List α × List α) (i : α) :
    i ∈ (p.1 ++ q.1, q.2 ++ p.2).1 ++ (p.1 ++ q.1, q.2 ++ p.2).2 ↔ i ∈ p.1 ++ p.2 ∨ i ∈ q.1 ++ q.2 := by
  simp only [List.mem_append]
  constructor <;> intro h <;> rcases h with (h | h) | h | h <;> simp [h]

theorem edgeItem_inj {g : Graph} {e e' : Nat} (h : edgeItem g e = edgeItem g e') : e = e' := by
  have : 1 + g.nv + e = 1 + g.nv + e' := h
  omega

namespace TEntry

/-- The upper terminal of an entry. -/
def top (s : WalkState) (t : TEntry) : Nat := s.stackVerts[t.topDepth]!

/-- The original edges covered by an entry: everything below either side's items. -/
def edges (g : Graph) (items : Items) (t : TEntry) (e : Nat) : Prop :=
  ∃ i ∈ t.spans.1 ++ t.spans.2, items.EdgeBelow g i e

/-- The vertices at which an entry may be attached while the walk is at depth `D`: its bottom
`vStart` and the path vertices at depths `topDepth..D`. -/
def Term (D : Nat) (s : WalkState) (t : TEntry) (v : Nat) : Prop :=
  v = t.vStart ∨ ∃ k, t.topDepth ≤ k ∧ k ≤ D ∧ v = s.stackVerts[k]!

/-- The result of `mergeTstackTops` on `cur :: nxt :: _`. -/
def mergeInto (cur nxt : TEntry) : TEntry :=
  { nxt with topDepth := min nxt.topDepth cur.topDepth,
             spans := (cur.spans.1 ++ nxt.spans.1, nxt.spans.2 ++ cur.spans.2) }

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

theorem edges_mergeInto (cur nxt : TEntry) (e : Nat) :
    (mergeInto cur nxt).edges g items e ↔ cur.edges g items e ∨ nxt.edges g items e := by
  simp only [edges, mergeInto, List.mem_append]
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

theorem mem_mergeInto (cur nxt : TEntry) (i : ItemId) :
    i ∈ (mergeInto cur nxt).spans.1 ++ (mergeInto cur nxt).spans.2 ↔
      i ∈ cur.spans.1 ++ cur.spans.2 ∨ i ∈ nxt.spans.1 ++ nxt.spans.2 :=
  mem_merge cur.spans nxt.spans i

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

theorem Term_congr {D : Nat} {s s' : WalkState} (hsv : s'.stackVerts = s.stackVerts)
    (hv : t'.vStart = t.vStart) (hd : t'.topDepth = t.topDepth) : t'.Term D s' = t.Term D s := by
  funext v; simp only [Term, hsv, hv, hd]

end TEntry

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

/-! ### The invariant -/

namespace WalkState

/-- Invariant on one tstack entry at walk depth `D`. -/
structure EntryInv (D : Nat) (s : WalkState) (t : TEntry) : Prop where
  conn : s.g.ConnEdges (t.edges s.g s.items)
  attached : s.g.AttachedIn (t.edges s.g s.items) (t.Term D s)

/-- Invariant on a closed item (its subtree is a 2-attached connected piece). -/
structure ItemInv (s : WalkState) (i : ItemId) : Prop where
  conn : s.g.ConnEdges (Items.EdgeBelow s.g s.items i)
  attached : ∀ u v, Items.vs s.items i = (some u, some v) →
    s.g.TwoAttached (Items.EdgeBelow s.g s.items i) u v

/-- The walk invariant at depth `D`: every open entry is a connected piece attached within its
terminals, every allocated node is a 2-attached connected piece. -/
structure Inv (D : Nat) (s : WalkState) : Prop where
  entries : ∀ t ∈ s.tstack, s.EntryInv D t
  nodes : ∀ i, 1 + s.g.nv + s.g.ne ≤ i → i < s.items.size → s.ItemInv i

/-- The item facts the closing steps rely on: ids of vertices/edges have their types, children
and span items are allocated. -/
structure Shape (s : WalkState) : Prop where
  size : 1 + s.g.nv + s.g.ne ≤ s.items.size
  root : Items.type s.items rootItem = .F
  vert : ∀ v, v < s.g.nv → Items.type s.items (vertItem v) = .V
  edge : ∀ e, e < s.g.ne → Items.type s.items (edgeItem s.g e) = .Q
  ch_lt : ∀ p c, Items.IsParent s.items p c → c < s.items.size
  span : ∀ t ∈ s.tstack, ∀ i ∈ t.spans.1 ++ t.spans.2, i < s.items.size

variable {D : Nat} {s s' : WalkState}

theorem EntryInv.congr {t t' : TEntry} (hg : s'.g = s.g) (hsv : s'.stackVerts = s.stackVerts)
    (hv : t'.vStart = t.vStart) (hd : t'.topDepth = t.topDepth)
    (hE : ∀ e, e < s.g.ne → (t'.edges s.g s'.items e ↔ t.edges s.g s.items e)) (h : s.EntryInv D t) :
    s'.EntryInv D t' := by
  obtain ⟨conn, att⟩ := h
  refine ⟨?_, ?_⟩ <;> rw [hg]
  · exact (Graph.ConnEdges.congr hE).2 conn
  · rw [TEntry.Term_congr hsv hv hd]; exact (Graph.AttachedIn.congr hE).2 att

theorem EntryInv.mono {D' : Nat} {t : TEntry} (hD : D ≤ D') (h : s.EntryInv D t) : s.EntryInv D' t :=
  ⟨h.conn, h.attached.mono fun v hv => by
    rcases hv with hv | ⟨k, h1, h2, h3⟩
    · exact .inl hv
    · exact .inr ⟨k, h1, Nat.le_trans h2 hD, h3⟩⟩

theorem ItemInv.congr {i : ItemId} (hg : s'.g = s.g)
    (hvs : Items.vs s'.items i = Items.vs s.items i)
    (hE : ∀ e, e < s.g.ne → (Items.EdgeBelow s.g s'.items i e ↔ Items.EdgeBelow s.g s.items i e))
    (h : s.ItemInv i) : s'.ItemInv i := by
  obtain ⟨conn, att⟩ := h
  refine ⟨?_, ?_⟩ <;> rw [hg]
  · exact (Graph.ConnEdges.congr hE).2 conn
  · intro u v huv; rw [hvs] at huv; exact (Graph.TwoAttached.congr hE).2 (att u v huv)

theorem Inv.mono {D' : Nat} (hD : D ≤ D') (h : s.Inv D) : s.Inv D' :=
  ⟨fun t ht => (h.entries t ht).mono hD, h.nodes⟩

/-- Changing fields other than `g`, `stackVerts`, `items`, `tstack` keeps the invariant. -/
theorem Inv.frame (hg : s'.g = s.g) (hsv : s'.stackVerts = s.stackVerts) (hi : s'.items = s.items)
    (hts : s'.tstack = s.tstack) (h : s.Inv D) : s'.Inv D := by
  refine ⟨fun t ht => EntryInv.congr hg hsv rfl rfl (fun _ _ => by rw [hi]) (h.entries t (hts ▸ ht)),
    fun i hi' hsz => ItemInv.congr hg (by rw [hi]) (fun _ _ => by rw [hi]) (h.nodes i ?_ ?_)⟩
  · rwa [hg] at hi'
  · rwa [hi] at hsz

theorem Shape.frame (hg : s'.g = s.g) (hi : s'.items = s.items) (hts : s'.tstack = s.tstack)
    (h : Shape s) : Shape s' := by
  obtain ⟨size, root, vert, edge, ch_lt, span⟩ := h
  exact ⟨by rw [hg, hi]; exact size, by rw [hi]; exact root, by rw [hg, hi]; exact vert,
    by rw [hg, hi]; exact edge, by rw [hi]; exact ch_lt, by rw [hi, hts]; exact span⟩

/-- The attachment of an entry whose terminals above `topDepth` are its bottom, interior, or
untouched: it is 2-attached at `vStart` and `stackVerts[topDepth]`. -/
theorem TwoAttached.of_term {E : Nat → Prop} {t : TEntry} (h : s.g.AttachedIn E (t.Term D s))
    (hmid : ∀ k, t.topDepth < k → k ≤ D → s.stackVerts[k]! = t.vStart ∨
      s.g.Interior E s.stackVerts[k]! ∨ ¬ s.g.Touches E s.stackVerts[k]!) :
    s.g.TwoAttached E t.vStart s.stackVerts[t.topDepth]! := by
  refine h.strengthen.mono fun v ⟨hv, htouch, hint⟩ => ?_
  rcases hv with hv | ⟨k, h1, h2, rfl⟩
  · exact .inl hv
  · rcases Nat.eq_or_lt_of_le h1 with rfl | h1
    · exact .inr rfl
    · rcases hmid k h1 h2 with h | h | h
      · exact .inl h
      · exact absurd h hint
      · exact absurd htouch h

end WalkState

/-! ### Running the primitives -/

section Run
variable {α β : Type}

theorem WalkM.run_bind (x : WalkM α) (f : α → WalkM β) (s : WalkState) :
    (x >>= f).run s = (f (x.run s).1).run (x.run s).2 := rfl
theorem List.head!_cons' {γ : Type} [Inhabited γ] (a : γ) (l : List γ) : (a :: l).head! = a := rfl
theorem WalkM.get_run (s : WalkState) : (get : WalkM WalkState).run s = (s, s) := rfl
theorem WalkM.pure_run (a : α) (s : WalkState) : (pure a : WalkM α).run s = (a, s) := rfl
theorem WalkM.modify_run (f : WalkState → WalkState) (s : WalkState) :
    (modify f : WalkM Unit).run s = ((), f s) := rfl
theorem WalkM.map_run (f : α → β) (x : WalkM α) (s : WalkState) :
    (f <$> x).run s = (f (x.run s).1, (x.run s).2) := rfl

theorem mergeTstackTops_run_eq (s : WalkState) (cur nxt : TEntry) (rest : List TEntry)
    (hs : s.tstack = cur :: nxt :: rest) :
    mergeTstackTops.run s = ((), { s with tstack := TEntry.mergeInto cur nxt :: rest }) := by
  obtain ⟨g, tern, items, sv, sd, nei, fo, tstack, tb, tsl⟩ := s
  subst hs; rfl

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

theorem maybeUnwrapNxt_run_eq (ty : NodeType) (s : WalkState) (a b : TEntry) (rest : List TEntry)
    (hs : s.tstack = a :: b :: rest) (dir : Bool) (hdir : dir = s.stackDir[b.topDepth]!) (i : ItemId)
    (hi : i = (getSide b.spans dir).head!) :
    (maybeUnwrapNxt ty).run s =
      if ty = .R ∨ s.ternarize = true then (allocItem ty).run s
      else if s.items[i]!.type = ty then
        (i, { s with tstack := a :: { b with spans := setSides dir s.items[i]!.ch [] } :: rest })
      else (allocItem ty).run s := by
  obtain ⟨g, tern, items, sv, sd, nei, fo, tstack, tb, tsl⟩ := s
  subst hs hdir hi
  simp only [maybeUnwrapNxt, WalkM.run_bind, WalkM.get_run, modifyNxt, Bool.or_eq_true, beq_iff_eq]
  by_cases h : ty = NodeType.R ∨ tern = true
  · simp only [h, ↓reduceIte]
  · simp only [h, ↓reduceIte]
    split <;> rename_i h' <;>
      simp only [WalkM.run_bind, run_nxt, run_stackDir, run_getItem,
        List.tail_cons, List.head!_cons', h', ↓reduceIte]
    rfl

theorem pure_tstack_eq (s : WalkState) : ({ s with tstack := s.tstack } : WalkState) = s := rfl

end Run

/-! ### Shape preservation -/

namespace WalkState
variable {D : Nat} {s s' : WalkState}

/-- Items of type S/P/R are allocated nodes. -/
theorem Shape.node_of_type (h : Shape s) {i : ItemId} (hi : i < s.items.size) {ty : NodeType}
    (hty : ty ∉ [NodeType.F, .V, .Q]) (hP : Items.type s.items i = ty) : 1 + s.g.nv + s.g.ne ≤ i := by
  refine Decidable.byContradiction fun hlt => ?_
  rcases Nat.lt_or_ge i (1 + s.g.nv) with h1 | h1
  · rcases Nat.eq_zero_or_pos i with rfl | h0
    · have hr : Items.type s.items 0 = .F := h.root
      exact hty (by rw [← hP]; exact (congrArg (fun x => x ∈ [NodeType.F, .V, .Q]) hr).mpr (by simp))
    · have hv := h.vert (i - 1) (by omega)
      have : vertItem (i - 1) = i := show (1 + (i - 1) : Nat) = i by omega
      rw [this, hP] at hv; exact hty (by rw [hv]; simp)
  · have hv := h.edge (i - (1 + s.g.nv)) (by omega)
    have : edgeItem s.g (i - (1 + s.g.nv)) = i := show (1 + s.g.nv + (i - (1 + s.g.nv)) : Nat) = i by omega
    rw [this, hP] at hv; exact hty (by rw [hv]; simp)

theorem Shape.edgeItem_ne (_h : Shape s) {e i : Nat} (_he : e < s.g.ne) (hi : 1 + s.g.nv + s.g.ne ≤ i) :
    edgeItem s.g e ≠ i := by
  show 1 + s.g.nv + e ≠ i; omega

theorem Shape.no_parent_size (h : Shape s) (p : ItemId) : ¬ Items.IsParent s.items p s.items.size :=
  fun hp => Nat.lt_irrefl _ (h.ch_lt p _ hp)

theorem Shape.tstack (h : Shape s) {l : List TEntry}
    (hl : ∀ t ∈ l, ∀ i ∈ t.spans.1 ++ t.spans.2, i < s.items.size) : Shape { s with tstack := l } :=
  ⟨h.size, h.root, h.vert, h.edge, h.ch_lt, hl⟩

theorem Shape.stackDir (h : Shape s) (sd : Array Bool) : Shape { s with stackDir := sd } :=
  ⟨h.size, h.root, h.vert, h.edge, h.ch_lt, h.span⟩

theorem Shape.push (h : Shape s) (ty : NodeType) :
    Shape { s with items := s.items.push ⟨ty, (none, none), []⟩ } := by
  have hch := Items.ch_push_nil (items := s.items) ⟨ty, (none, none), []⟩ rfl
  have hsz := h.size
  have h0 : rootItem ≠ s.items.size := by show (0 : Nat) ≠ _; omega
  have hv : ∀ v, vertItem v ≠ s.items.size ∨ s.g.nv ≤ v := fun v => by
    show (1 + v : Nat) ≠ _ ∨ _; omega
  have he : ∀ e, edgeItem s.g e ≠ s.items.size ∨ s.g.ne ≤ e := fun e => by
    show (1 + s.g.nv + e : Nat) ≠ _ ∨ _; omega
  refine ⟨?_, ?_, fun v hv' => ?_, fun e he' => ?_, fun p c hp => ?_, fun t ht i hi => ?_⟩
  · simp only [Array.size_push]; omega
  · rw [Items.type_push_of_ne _ h0]; exact h.root
  · rw [Items.type_push_of_ne _ ((hv v).resolve_right (Nat.not_le.2 hv'))]; exact h.vert v hv'
  · rw [Items.type_push_of_ne _ ((he e).resolve_right (Nat.not_le.2 he'))]; exact h.edge e he'
  · rw [Items.IsParent, hch] at hp; simp only [Array.size_push]; exact Nat.lt_succ_of_lt (h.ch_lt p c hp)
  · simp only [Array.size_push]; exact Nat.lt_succ_of_lt (h.span t ht i hi)

theorem Shape.modify (h : Shape s) (j : ItemId) (f : Item → Item) (hty : ∀ it, (f it).type = it.type)
    (hch : ∀ hj : j < s.items.size, ∀ c ∈ (f s.items[j]).ch, c < s.items.size) :
    Shape { s with items := s.items.modify j f } := by
  have hT := Items.type_modify_type_eq (items := s.items) j f hty
  refine ⟨by simpa using h.size, by rw [hT]; exact h.root, fun v hv => by rw [hT]; exact h.vert v hv,
    fun e he => by rw [hT]; exact h.edge e he, fun p c hp => ?_, fun t ht i hi => by simpa using h.span t ht i hi⟩
  simp only [Array.size_modify]
  rcases Items.IsParent_modify hp with hp | ⟨rfl, hj, hc⟩
  · exact h.ch_lt p c hp
  · exact hch hj c hc

end WalkState

/-! ### Invariant preservation by the primitives -/

namespace WalkState
variable {D : Nat} {s s' : WalkState}

/-- Allocating a fresh node item preserves the invariant: it has no edges below it yet. -/
theorem Inv.alloc (ty : NodeType) (h : s.Inv D) (hsize : 1 + s.g.nv + s.g.ne ≤ s.items.size) :
    Inv D { s with items := s.items.push ⟨ty, (none, none), []⟩ } := by
  have hB : ∀ a i, Items.Below (s.items.push ⟨ty, (none, none), []⟩) a i ↔ Items.Below s.items a i :=
    fun _ _ => Items.Below_push_nil _ rfl
  refine ⟨fun t ht => EntryInv.congr (s := s) rfl rfl rfl rfl (fun e _ => TEntry.edges_congr (fun i _ e => hB i _) e)
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
  · exact ItemInv.congr (s := s) rfl (Items.vs_push_of_ne _ hi') (fun e _ => hB i _) (h.nodes i hi (by omega))

/-- Recording `vs` on a non-node item (an edge item) preserves the invariant. -/
theorem Inv.modifyVs (j : ItemId) (vsv : Option Nat × Option Nat) (h : s.Inv D) (hj : j < 1 + s.g.nv + s.g.ne) :
    Inv D { s with items := s.items.modify j fun it => { it with vs := vsv } } := by
  have hB : ∀ a i, Items.Below (s.items.modify j fun it => { it with vs := vsv }) a i ↔ Items.Below s.items a i :=
    fun _ _ => Items.Below_modify_ch_eq j (fun it => { it with vs := vsv }) fun _ => rfl
  refine ⟨fun t ht => EntryInv.congr (s := s) rfl rfl rfl rfl (fun e _ => TEntry.edges_congr (fun i _ e => hB i _) e)
    (h.entries t ht), fun i hi hsz => ?_⟩
  have hsz' : i < s.items.size := by simpa using hsz
  have hne : i ≠ j := by intro h; subst h; exact absurd hi (Nat.not_le.2 hj)
  exact ItemInv.congr (s := s) rfl (Items.vs_modify_of_ne _ _ hne) (fun e _ => hB i _) (h.nodes i hi hsz')

/-- Pushing the vertex item of `v` as a new entry. -/
theorem Inv.pushVert (v d : Nat) (h : s.Inv D)
    (hc : s.g.ConnEdges (Items.EdgeBelow s.g s.items (vertItem v)))
    (ha : s.g.TwoAttached (Items.EdgeBelow s.g s.items (vertItem v)) v v) :
    ((pushVertTstack v d).run s).2.Inv D := by
  rw [pushVertTstack, run_pushTstack]
  refine ⟨fun t ht => ?_, fun i hi hsz => ItemInv.congr (s := s) rfl rfl (fun _ _ => Iff.rfl) (h.nodes i hi hsz)⟩
  simp only [List.mem_cons] at ht
  rcases ht with rfl | ht
  · have hE : ∀ e, TEntry.edges s.g s.items ⟨v, d, s.nxtEdgeIdx, setSides s.stackDir[d]! [vertItem v] []⟩ e ↔
        Items.EdgeBelow s.g s.items (vertItem v) e := by
      intro e; simp only [TEntry.edges, mem_setSides, List.mem_singleton, exists_eq_left]
    refine ⟨(Graph.ConnEdges.congr fun e _ => hE e).2 hc, ?_⟩
    refine ((Graph.AttachedIn.congr fun e _ => hE e).2 (Graph.twoAttached_iff.1 ha)).mono fun x hx => ?_
    rcases hx with rfl | rfl <;> exact .inl rfl
  · exact EntryInv.congr (s := s) rfl rfl rfl rfl (fun _ _ => Iff.rfl) (h.entries t ht)

/-- Pushing the edge `e = {vStart, stackVerts[topDepth]}` as a new entry at `topDepth ≤ D`. -/
theorem Inv.pushEdge (vStart topDepth e : Nat) (h : s.Inv D) (_he : e < s.g.ne)
    (hq : Items.ch s.items (edgeItem s.g e) = [])
    (hend : Items.PairEq (vStart, s.stackVerts[topDepth]!) s.g.edges[e]!) (hD : topDepth ≤ D) :
    ((pushEdgeTstack vStart topDepth e).run s).2.Inv D := by
  rw [run_pushEdgeTstack]
  refine ⟨fun t ht => ?_, fun i hi hsz => ItemInv.congr (s := s) rfl rfl (fun _ _ => Iff.rfl) (h.nodes i hi hsz)⟩
  simp only [List.mem_cons] at ht
  rcases ht with rfl | ht
  · have hE := TEntry.edges_edgeEntry (items := s.items) s.stackDir[topDepth]! vStart topDepth s.nxtEdgeIdx e hq
    refine ⟨Graph.ConnEdges.single hE, (Graph.twoAttached_iff.1 (Graph.TwoAttached.single_of_pairEq hE hend)).mono fun x hx => ?_⟩
    rcases hx with rfl | rfl
    · exact .inl rfl
    · exact .inr ⟨topDepth, Nat.le_refl _, hD, rfl⟩
  · exact EntryInv.congr (s := s) rfl rfl rfl rfl (fun _ _ => Iff.rfl) (h.entries t ht)

/-- Hypotheses under which merging `cur` into `nxt` keeps the invariant: the two edge sets share a
vertex (unless one is empty), and `cur`'s bottom is a terminal of the merged entry or interior to
the union. -/
structure MergeOk (D : Nat) (s : WalkState) (cur nxt : TEntry) : Prop where
  share : (∃ e, e < s.g.ne ∧ cur.edges s.g s.items e) → (∃ e, e < s.g.ne ∧ nxt.edges s.g s.items e) →
    ∃ x, s.g.Touches (cur.edges s.g s.items) x ∧ s.g.Touches (nxt.edges s.g s.items) x
  bottom : (TEntry.mergeInto cur nxt).Term D s cur.vStart ∨
    s.g.Interior (fun e => cur.edges s.g s.items e ∨ nxt.edges s.g s.items e) cur.vStart

theorem MergeOk.congr {cur nxt : TEntry} (hg : s'.g = s.g) (hsv : s'.stackVerts = s.stackVerts)
    (hc : ∀ e, e < s.g.ne → (cur.edges s.g s'.items e ↔ cur.edges s.g s.items e))
    (hn : ∀ e, e < s.g.ne → (nxt.edges s.g s'.items e ↔ nxt.edges s.g s.items e))
    (h : MergeOk D s cur nxt) : MergeOk D s' cur nxt := by
  obtain ⟨share, bottom⟩ := h
  refine ⟨?_, ?_⟩ <;> rw [hg]
  · rintro ⟨e, he, hc'⟩ ⟨e', he', hn'⟩
    obtain ⟨x, ⟨e₁, he₁, h₁, hi₁⟩, ⟨e₂, he₂, h₂, hi₂⟩⟩ := share ⟨e, he, (hc e he).1 hc'⟩ ⟨e', he', (hn e' he').1 hn'⟩
    exact ⟨x, ⟨e₁, he₁, (hc e₁ he₁).2 h₁, hi₁⟩, ⟨e₂, he₂, (hn e₂ he₂).2 h₂, hi₂⟩⟩
  · rw [TEntry.Term_congr hsv rfl rfl]
    rcases bottom with h | h
    · exact .inl h
    · refine .inr fun e he hi => ?_
      have := h e he hi
      show TEntry.edges s.g s'.items cur e ∨ TEntry.edges s.g s'.items nxt e
      rw [hc e he, hn e he]; exact this

/-- The merged entry satisfies the invariant. -/
theorem EntryInv.merge {cur nxt : TEntry} (hc : s.EntryInv D cur) (hn : s.EntryInv D nxt)
    (hok : MergeOk D s cur nxt) : s.EntryInv D (TEntry.mergeInto cur nxt) := by
  have hm := TEntry.edges_mergeInto (g := s.g) (items := s.items) cur nxt
  refine ⟨(Graph.ConnEdges.congr fun e _ => hm e).2 (Graph.ConnEdges.union hc.conn hn.conn hok.share), ?_⟩
  refine (Graph.AttachedIn.congr fun e _ => hm e).2 ((Graph.AttachedIn.union hc.attached hn.attached).mono ?_)
  rintro v ⟨hv, hint⟩
  rcases hv with (rfl | ⟨k, h1, h2, h3⟩) | (hv | ⟨k, h1, h2, h3⟩)
  · rcases hok.bottom with h | h
    · exact h
    · exact absurd h hint
  · exact .inr ⟨k, Nat.le_trans (Nat.min_le_right _ _) h1, h2, h3⟩
  · exact .inl hv
  · exact .inr ⟨k, Nat.le_trans (Nat.min_le_left _ _) h1, h2, h3⟩

/-- Soundness of `mergeTstackTops`. -/
theorem mergeTstackTops_sound (cur nxt : TEntry) (rest : List TEntry) (hs : s.tstack = cur :: nxt :: rest)
    (h : s.Inv D) (hok : MergeOk D s cur nxt) : (mergeTstackTops.run s).2.Inv D := by
  rw [mergeTstackTops_run_eq s cur nxt rest hs]
  refine ⟨fun t ht => ?_, fun i hi hsz => ItemInv.congr (s := s) rfl rfl (fun _ _ => Iff.rfl) (h.nodes i hi hsz)⟩
  simp only [List.mem_cons] at ht
  rcases ht with rfl | ht
  · exact EntryInv.congr (s := s) rfl rfl rfl rfl (fun _ _ => Iff.rfl)
      (EntryInv.merge (h.entries cur (by simp [hs])) (h.entries nxt (by simp [hs])) hok)
  · exact EntryInv.congr (s := s) rfl rfl rfl rfl (fun _ _ => Iff.rfl) (h.entries t (by simp [hs, ht]))

/-- Closing the top entry `t` into `item` records a 2-attached connected piece with terminals
`makeVs t.vStart t.topDepth`. `item` is an allocated node that is not yet attached anywhere (no
parent, in no entry's spans), `t`'s items all lie on the side `stackDir[t.topDepth]`, and `t` has no
attachment at the intermediate depths `(t.topDepth, D]`. -/
theorem finishTstackTop_complete (item : ItemId) (t : TEntry) (rest : List TEntry)
    (hs : s.tstack = t :: rest) (h : s.Inv D)
    (hitem : item < s.items.size) (hnode : 1 + s.g.nv + s.g.ne ≤ item)
    (hroot : ∀ p, ¬ Items.IsParent s.items p item)
    (hfree : ∀ t' ∈ s.tstack, item ∉ t'.spans.1 ++ t'.spans.2)
    (hside : getSide t.spans (!s.stackDir[t.topDepth]!) = [])
    (hmid : ∀ k, t.topDepth < k → k ≤ D → s.stackVerts[k]! = t.vStart ∨
      s.g.Interior (t.edges s.g s.items) s.stackVerts[k]! ∨ ¬ s.g.Touches (t.edges s.g s.items) s.stackVerts[k]!) :
    ((finishTstackTop item).run s).2.Inv D := by
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
    · refine EntryInv.congr (s := s) (t := t) rfl rfl rfl rfl (fun e he => ?_) ht
      show TEntry.edges s.g (s.items.modify item f) _ e ↔ _
      rw [← hsub e he]
      simp only [TEntry.edges]
      constructor
      · rintro ⟨i, hi, hb⟩
        rwa [List.mem_singleton.1 ((mem_setSides dir [item] i).1 hi)] at hb
      · exact fun hb => ⟨item, (mem_setSides dir [item] item).2 (List.mem_singleton.2 rfl), hb⟩
    · exact EntryInv.congr (s := s) rfl rfl rfl rfl
        (fun e _ => TEntry.edges_modify_of_not_mem item f hroot (hfree t' (by simp [hs, ht'])) e)
        (h.entries t' (by simp [hs, ht']))
  · intro i hi hsize
    simp only [Array.size_modify] at hsize
    by_cases hi' : i = item
    · subst hi'
      refine ⟨(Graph.ConnEdges.congr hsub).2 ht.conn, ?_⟩
      intro u v huv
      rw [Items.vs_modify_at i f hitem] at huv
      have h2 := (Graph.TwoAttached.congr hsub).2 (TwoAttached.of_term ht.attached hmid)
      rcases setSides_eq _ _ _ _ _ huv with ⟨hu, hv⟩ | ⟨hu, hv⟩ <;>
        simp only [Option.some.injEq] at hu hv <;> subst hu hv
      · exact h2.comm
      · exact h2
    · exact ItemInv.congr (s := s) rfl (Items.vs_modify_of_ne item f hi')
        (fun _ _ => Items.Below_modify_of_not_below item f (hnb i hi')) (h.nodes i hi hsize)

/-! ### Steps -/

/-- One walk step that keeps the invariant and the shape and does not disturb the graph, the
vertex stack, or the subtree of `vertItem v`. -/
structure Step (D v : Nat) (s s' : WalkState) : Prop where
  inv : s'.Inv D
  shape : Shape s'
  g : s'.g = s.g
  sv : s'.stackVerts = s.stackVerts
  below : ∀ i, Items.Below s'.items (vertItem v) i ↔ Items.Below s.items (vertItem v) i

variable {v : Nat}

theorem Step.refl (hi : s.Inv D) (hs : Shape s) : Step D v s s := ⟨hi, hs, rfl, rfl, fun _ => Iff.rfl⟩

theorem Step.trans {s₁ s₂ s₃ : WalkState} (h₁ : Step D v s₁ s₂) (h₂ : Step D v s₂ s₃) : Step D v s₁ s₃ :=
  ⟨h₂.inv, h₂.shape, h₂.g.trans h₁.g, h₂.sv.trans h₁.sv, fun i => (h₂.below i).trans (h₁.below i)⟩

theorem Step.of_eq (h : Step D v s s') {s'' : WalkState} (he : s'' = s') : Step D v s s'' := he ▸ h

theorem Step.frame (hi : s.Inv D) (hs : Shape s) (hg : s'.g = s.g) (hsv : s'.stackVerts = s.stackVerts)
    (hitems : s'.items = s.items) (hts : s'.tstack = s.tstack) : Step D v s s' :=
  ⟨hi.frame hg hsv hitems hts, hs.frame hg hitems hts, hg, hsv, fun _ => by rw [hitems]⟩

theorem Step.alloc (hi : s.Inv D) (hs : Shape s) (ty : NodeType) :
    Step D v s { s with items := s.items.push ⟨ty, (none, none), []⟩ } :=
  ⟨hi.alloc ty hs.size, hs.push ty, rfl, rfl, fun _ => Items.Below_push_nil _ rfl⟩

theorem Step.modifyVs (hi : s.Inv D) (hs : Shape s) (j : ItemId) (vsv : Option Nat × Option Nat)
    (hj : j < 1 + s.g.nv + s.g.ne) :
    Step D v s { s with items := s.items.modify j fun it => { it with vs := vsv } } := by
  refine ⟨hi.modifyVs j vsv hj, ?_, rfl, rfl, fun _ => Items.Below_modify_ch_eq j (fun it => { it with vs := vsv }) fun _ => rfl⟩
  refine hs.modify j (fun it => { it with vs := vsv }) (fun _ => rfl) fun hj' c hc => hs.ch_lt j c ?_
  rw [Items.IsParent, Items.ch_eq_getElem hj']
  exact hc

theorem Step.pushVert (hi : s.Inv D) (hs : Shape s) (d : Nat) (hv : v < s.g.nv)
    (hc : s.g.ConnEdges (Items.EdgeBelow s.g s.items (vertItem v)))
    (ha : s.g.TwoAttached (Items.EdgeBelow s.g s.items (vertItem v)) v v) :
    Step D v s ((pushVertTstack v d).run s).2 := by
  refine ⟨hi.pushVert v d hc ha, ?_, ?_, ?_, ?_⟩ <;> rw [pushVertTstack, run_pushTstack]
  · refine hs.tstack fun t ht i hit => ?_
    simp only [List.mem_cons] at ht
    rcases ht with rfl | ht
    · rw [mem_setSides, List.mem_singleton] at hit; subst hit
      have := hs.size; show 1 + v < _; omega
    · exact hs.span t ht i hit
  all_goals first | rfl | exact fun _ => Iff.rfl

theorem Step.pushEdge (hi : s.Inv D) (hs : Shape s) (vStart topDepth e : Nat) (he : e < s.g.ne)
    (hq : Items.ch s.items (edgeItem s.g e) = [])
    (hend : Items.PairEq (vStart, s.stackVerts[topDepth]!) s.g.edges[e]!) (hD : topDepth ≤ D) :
    Step D v s ((pushEdgeTstack vStart topDepth e).run s).2 := by
  refine ⟨hi.pushEdge vStart topDepth e he hq hend hD, ?_, ?_, ?_, ?_⟩ <;> rw [run_pushEdgeTstack]
  · refine hs.tstack fun t ht i hit => ?_
    simp only [List.mem_cons] at ht
    rcases ht with rfl | ht
    · rw [mem_setSides, List.mem_singleton] at hit; subst hit
      have := hs.size; show 1 + s.g.nv + e < _; omega
    · exact hs.span t ht i hit
  all_goals first | rfl | exact fun _ => Iff.rfl

/-- `mergeTstackTops` keeps the invariant when the top two entries satisfy `MergeOk` (on a stack
with fewer than two entries it empties the stack). -/
def MergeTopOk (D : Nat) (s : WalkState) : Prop :=
  ∀ cur nxt rest, s.tstack = cur :: nxt :: rest → MergeOk D s cur nxt

theorem Step.mergeTop (hi : s.Inv D) (hs : Shape s) (hok : MergeTopOk D s) :
    Step D v s (mergeTstackTops.run s).2 := by
  match hts : s.tstack with
  | [] | [_] =>
    rw [run_mergeTstackTops, hts]
    exact ⟨⟨fun t ht => by simp [mergeTops] at ht, fun i hi' hsz => ItemInv.congr (s := s) rfl rfl (fun _ _ => Iff.rfl) (hi.nodes i hi' hsz)⟩,
      hs.tstack (by simp [mergeTops]), rfl, rfl, fun _ => Iff.rfl⟩
  | cur :: nxt :: rest =>
    refine ⟨mergeTstackTops_sound cur nxt rest hts hi (hok cur nxt rest hts), ?_, ?_, ?_, ?_⟩ <;>
      rw [mergeTstackTops_run_eq s cur nxt rest hts]
    · refine hs.tstack fun t ht i hit => ?_
      simp only [List.mem_cons] at ht
      rcases ht with rfl | ht
      · rcases (TEntry.mem_mergeInto cur nxt i).1 hit with h | h
        · exact hs.span cur (by simp [hts]) i h
        · exact hs.span nxt (by simp [hts]) i h
      · exact hs.span t (by simp [hts, ht]) i hit
    all_goals first | rfl | exact fun _ => Iff.rfl

/-- The state after `k` iterations of `body`. -/
def iter (body : WalkM Unit) : Nat → WalkState → WalkState
  | 0, s => s
  | k + 1, s => iter body k (body.run s).2

theorem loop_succ_run (fuel : Nat) (cond : WalkM Bool) (body : WalkM Unit) (s : WalkState)
    (hcond : (cond.run s).2 = s) :
    (loop (fuel + 1) cond body).run s =
      if (cond.run s).1 = true then (loop fuel cond body).run (body.run s).2 else ((), s) := by
  simp only [loop, WalkM.run_bind, hcond]
  split <;> rfl

/-- A loop whose body is a step whenever its (state-independent) condition holds and the iteration
hypothesis `Ok` holds is a step. -/
theorem Step.loop (cond : WalkM Bool) (body : WalkM Unit) (Ok : WalkState → Prop) (fuel : Nat)
    (hcond : ∀ s, (cond.run s).2 = s)
    (hbody : ∀ s, v < s.g.nv → s.Inv D → Shape s → (cond.run s).1 = true → Ok s → Step D v s (body.run s).2)
    (hv : v < s.g.nv) (hi : s.Inv D) (hs : Shape s)
    (hok : ∀ k, (∀ j, j ≤ k → (cond.run (iter body j s)).1 = true) → Ok (iter body k s)) :
    Step D v s ((loop fuel cond body).run s).2 := by
  induction fuel generalizing s with
  | zero => exact Step.refl hi hs
  | succ fuel ih =>
    rw [loop_succ_run fuel cond body s (hcond s)]
    by_cases hc : (cond.run s).1 = true
    · simp only [hc, ↓reduceIte]
      have st := hbody s hv hi hs hc (hok 0 fun j hj => by rw [Nat.le_zero.1 hj]; exact hc)
      refine st.trans (ih (by rw [st.g]; exact hv) st.inv st.shape fun k hk => hok (k + 1) fun j hj => ?_)
      cases j with
      | zero => exact hc
      | succ j => exact hk j (Nat.le_of_succ_le_succ hj)
    · simp only [hc]; exact Step.refl hi hs


/-! ### Atomic steps with stack-shape hypotheses -/

/-- The state after running `m` from `s`. -/
def after {α : Type} (m : WalkM α) (s : WalkState) : WalkState := (m.run s).2
/-- The result of running `m` from `s`. -/
def result {α : Type} (m : WalkM α) (s : WalkState) : α := (m.run s).1
/-- The top entry (`cur`). -/
def curE (s : WalkState) : TEntry := s.tstack.head!
/-- The entry below the top (`nxt`). -/
def nxtE (s : WalkState) : TEntry := s.tstack.tail.head!
/-- The side of `nxt`. -/
def nxtDir (s : WalkState) : Bool := s.stackDir[(nxtE s).topDepth]!
/-- The item `maybeUnwrapNxt` inspects. -/
def nxtHead (s : WalkState) : ItemId := (getSide (nxtE s).spans (nxtDir s)).head!

/-- An allocated node that is attached nowhere yet: no parent and in no entry's spans. -/
structure ItemFree (s : WalkState) (item : ItemId) : Prop where
  lt : item < s.items.size
  node : 1 + s.g.nv + s.g.ne ≤ item
  root : ∀ p, ¬ Items.IsParent s.items p item
  free : ∀ t ∈ s.tstack, item ∉ t.spans.1 ++ t.spans.2

theorem ItemFree.frame {item : ItemId} (hg : s'.g = s.g) (hitems : s'.items = s.items)
    (hsp : ∀ t ∈ s'.tstack, ∀ i ∈ t.spans.1 ++ t.spans.2, ∃ t' ∈ s.tstack, i ∈ t'.spans.1 ++ t'.spans.2)
    (h : ItemFree s item) : ItemFree s' item :=
  ⟨by rw [hitems]; exact h.lt, by rw [hg]; exact h.node, by rw [hitems]; exact h.root,
   fun t ht hi => let ⟨t', ht', hi'⟩ := hsp t ht item hi; h.free t' ht' hi'⟩

theorem ItemFree.merge {item : ItemId} (h : ItemFree s item) : ItemFree (mergeTstackTops.run s).2 item := by
  match hts : s.tstack with
  | [] | [_] =>
    rw [run_mergeTstackTops, hts]
    exact ItemFree.frame (s := s) rfl rfl (fun t ht => by simp [mergeTops] at ht) h
  | cur :: nxt :: rest =>
    rw [mergeTstackTops_run_eq s cur nxt rest hts]
    refine ItemFree.frame (s := s) rfl rfl (fun t ht i hi => ?_) h
    simp only [List.mem_cons] at ht
    rcases ht with rfl | ht
    · rcases (TEntry.mem_mergeInto cur nxt i).1 hi with h | h
      · exact ⟨cur, by simp [hts], h⟩
      · exact ⟨nxt, by simp [hts], h⟩
    · exact ⟨t, by simp [hts, ht], hi⟩

/-- `finishTstackTop` on `s`: the top entry is one-sided at `stackDir[topDepth]` and has no attachment
at the intermediate depths `(topDepth, D]`. -/
structure FinishTopOk (D : Nat) (s : WalkState) : Prop where
  nonempty : s.tstack ≠ []
  side : getSide (curE s).spans (!s.stackDir[(curE s).topDepth]!) = []
  mid : ∀ k, (curE s).topDepth < k → k ≤ D → s.stackVerts[k]! = (curE s).vStart ∨
    s.g.Interior ((curE s).edges s.g s.items) s.stackVerts[k]! ∨
    ¬ s.g.Touches ((curE s).edges s.g s.items) s.stackVerts[k]!

theorem Step.finishTop (hi : s.Inv D) (hs : Shape s) (hv : v < s.g.nv) (item : ItemId)
    (hok : FinishTopOk D s) (hf : ItemFree s item) : Step D v s ((finishTstackTop item).run s).2 := by
  match hts : s.tstack with
  | [] => exact absurd hts hok.nonempty
  | t :: rest =>
    have hc : curE s = t := by rw [curE, hts, List.head!_cons']
    have hside := hok.side
    have hmid := hok.mid
    rw [hc] at hside hmid
    have hnb : ¬ Items.Below s.items (vertItem v) item := fun hb => by
      have h1 : 1 + v = item := Items.Below.eq_of_no_parent hf.root hb
      have h2 := hf.node; omega
    refine ⟨finishTstackTop_complete item t rest hts hi hf.lt hf.node hf.root hf.free hside hmid, ?_, ?_, ?_, ?_⟩ <;>
      rw [finishTstackTop_run_eq s item t rest hts]
    · refine (hs.modify item (fun it =>
          { it with vs := setSides s.stackDir[t.topDepth]! (some s.stackVerts[t.topDepth]!) (some t.vStart),
                    ch := getSide t.spans s.stackDir[t.topDepth]! }) (fun _ => rfl) fun _ c hc' =>
          hs.span t (by simp [hts]) c ((mem_of_getSide_nil _ t.spans hside c).2 hc')).tstack fun t' ht' i hit => ?_
      simp only [Array.size_modify]
      simp only [List.mem_cons] at ht'
      rcases ht' with rfl | ht'
      · rw [mem_setSides, List.mem_singleton] at hit; subst hit; exact hf.lt
      · exact hs.span t' (by simp [hts, ht']) i hit
    all_goals first | rfl | exact fun _ => Items.Below_modify_of_not_below item _ hnb

/-- The shape `maybeUnwrapNxt` relies on when it reuses the item `nxtHead s` heading `nxt`: `nxt` is
exactly that item on its side, and the item is a root lying in no other entry. -/
structure UnwrapAt (s : WalkState) : Prop where
  side : getSide (nxtE s).spans (!nxtDir s) = []
  single : getSide (nxtE s).spans (nxtDir s) = [nxtHead s]
  root : ∀ p, ¬ Items.IsParent s.items p (nxtHead s)
  fresh : ∀ t ∈ curE s :: s.tstack.tail.tail, nxtHead s ∉ t.spans.1 ++ t.spans.2

/-- `maybeUnwrapNxt ty` on `s`: at least two entries, and `UnwrapAt` whenever it actually unwraps. -/
structure UnwrapOk (ty : NodeType) (s : WalkState) : Prop where
  two : 2 ≤ s.tstack.length
  unwrap : ¬ (ty = .R ∨ s.ternarize = true) → s.items[nxtHead s]!.type = ty → UnwrapAt s

/-- What `maybeUnwrapNxt` guarantees: a step, and the returned item is free. -/
structure UnwrapRes (D v : Nat) (s : WalkState) (r : ItemId × WalkState) : Prop where
  step : Step D v s r.2
  free : ItemFree r.2 r.1

theorem Step.allocRes (hi : s.Inv D) (hs : Shape s) (ty : NodeType) :
    UnwrapRes D v s ((allocItem ty).run s) := by
  rw [run_allocItem]
  refine ⟨Step.alloc hi hs ty, by simp, hs.size, fun p hp => ?_, fun t ht hmem => ?_⟩
  · rw [Items.IsParent, Items.ch_push_nil _ rfl] at hp
    exact hs.no_parent_size p hp
  · exact Nat.lt_irrefl _ (hs.span t ht _ hmem)

theorem maybeUnwrapNxt_spec {ty : NodeType} (hi : s.Inv D) (hs : Shape s) (hty : ty ∉ [NodeType.F, .V, .Q])
    (hok : UnwrapOk ty s) : UnwrapRes D v s ((maybeUnwrapNxt ty).run s) := by
  match hts : s.tstack with
  | [] | [_] => have := hok.two; rw [hts] at this; simp at this
  | a :: b :: rest =>
    have hn : nxtE s = b := by rw [nxtE, hts]; rfl
    have hc : curE s = a := by rw [curE, hts]; rfl
    have hd : nxtDir s = s.stackDir[b.topDepth]! := by rw [nxtDir, hn]
    have hh : nxtHead s = (getSide b.spans s.stackDir[b.topDepth]!).head! := by rw [nxtHead, hn, hd]
    rw [maybeUnwrapNxt_run_eq ty s a b rest hts _ rfl _ rfl]
    by_cases h1 : ty = .R ∨ s.ternarize = true
    · simp only [h1, ↓reduceIte]; exact Step.allocRes hi hs ty
    simp only [h1, ↓reduceIte]
    by_cases h2 : s.items[(getSide b.spans s.stackDir[b.topDepth]!).head!]!.type = ty
    · simp only [h2, ↓reduceIte]
      obtain ⟨side, single, root, fresh⟩ := hok.unwrap h1 (by rw [hh]; exact h2)
      rw [hn, hd] at side single
      rw [hh] at single root fresh
      rw [hc, hts] at fresh
      simp only [List.tail_cons] at fresh
      set dir := s.stackDir[b.topDepth]!
      set i := (getSide b.spans dir).head!
      have hib : i ∈ b.spans.1 ++ b.spans.2 := by
        rw [mem_of_getSide_nil dir b.spans side, single]; exact List.mem_singleton_self i
      have hilt : i < s.items.size := hs.span b (by simp [hts]) i hib
      have hget : s.items[i]! = s.items[i] := getElem!_pos s.items i hilt
      have hity : Items.type s.items i = ty := by rw [Items.type_eq_getElem hilt, ← hget]; exact h2
      have hinode := hs.node_of_type hilt hty hity
      have hie : ∀ e, e < s.g.ne → edgeItem s.g e ≠ i := fun e he => hs.edgeItem_ne he hinode
      have hch : s.items[i]!.ch = Items.ch s.items i := by rw [Items.ch_eq_getElem hilt, hget]
      set b' : TEntry := { b with spans := setSides dir s.items[i]!.ch [] } with hb'
      have hb'E : ∀ e, e < s.g.ne → (b'.edges s.g s.items e ↔ b.edges s.g s.items e) := by
        intro e he
        rw [TEntry.edges_single dir i side single, hb', hch]
        exact TEntry.edges_unwrap dir b.vStart b.topDepth b.firstIdx i hie he
      have hinv : Inv D { s with tstack := a :: b' :: rest } := by
        refine ⟨fun t ht => ?_, fun j hj hsz => ItemInv.congr (s := s) rfl rfl (fun _ _ => Iff.rfl) (hi.nodes j hj hsz)⟩
        simp only [List.mem_cons] at ht
        rcases ht with ht | ht | ht
        · rw [ht]; exact EntryInv.congr (s := s) rfl rfl rfl rfl (fun _ _ => Iff.rfl) (hi.entries a (by simp [hts]))
        · rw [ht]; exact EntryInv.congr (s := s) (t := b) rfl rfl rfl rfl hb'E (hi.entries b (by simp [hts]))
        · exact EntryInv.congr (s := s) rfl rfl rfl rfl (fun _ _ => Iff.rfl) (hi.entries t (by simp [hts, ht]))
      have hshape : Shape { s with tstack := a :: b' :: rest } := by
        refine hs.tstack fun t ht j hj => ?_
        simp only [List.mem_cons] at ht
        rcases ht with ht | ht | ht
        · rw [ht] at hj; exact hs.span a (by simp [hts]) j hj
        · rw [ht, hb', mem_setSides, hch] at hj
          exact hs.ch_lt i j hj
        · exact hs.span t (by simp [hts, ht]) j hj
      refine ⟨⟨hinv, hshape, rfl, rfl, fun _ => Iff.rfl⟩, hilt, hinode, root, fun t ht hmem => ?_⟩
      simp only [List.mem_cons] at ht
      rcases ht with ht | ht | ht
      · rw [ht] at hmem; exact fresh a (by simp) hmem
      · rw [ht, hb', mem_setSides, hch] at hmem
        exact root i hmem
      · exact fresh t (by simp [ht]) hmem
    · simp only [h2, ↓reduceIte]; exact Step.allocRes hi hs ty

/-- `mergeTstackTops; finishTstackTop item`. -/
structure CloseTwoOk (D : Nat) (s : WalkState) : Prop where
  merge : MergeTopOk D s
  finish : FinishTopOk D (after mergeTstackTops s)

theorem Step.closeTwo (hi : s.Inv D) (hs : Shape s) (hv : v < s.g.nv) (item : ItemId)
    (hok : CloseTwoOk D s) (hf : ItemFree s item) :
    Step D v s ((finishTstackTop item).run (mergeTstackTops.run s).2).2 := by
  have st := Step.mergeTop (v := v) hi hs hok.merge
  exact st.trans (Step.finishTop st.inv st.shape (by rw [st.g]; exact hv) item hok.finish hf.merge)

/-- Re-targeting the top entry's bottom to `curV` (the `modifyCur` of the vertex close): the old
bottom is `curV`, interior, or a path vertex. -/
structure RetargetOk (D curV : Nat) (s : WalkState) : Prop where
  nonempty : s.tstack ≠ []
  old : (curE s).vStart = curV ∨ s.g.Interior ((curE s).edges s.g s.items) (curE s).vStart ∨
    ∃ k, (curE s).topDepth ≤ k ∧ k ≤ D ∧ (curE s).vStart = s.stackVerts[k]!

def retarget (curV : Nat) (edgeDir : Bool) : WalkM Unit :=
  modifyCur fun t => { t with vStart := curV, spans := setSides (!edgeDir) (t.spans.1 ++ t.spans.2) [] }

theorem retarget_run_eq (curV : Nat) (edgeDir : Bool) (s : WalkState) (t : TEntry) (rest : List TEntry)
    (hs : s.tstack = t :: rest) :
    (retarget curV edgeDir).run s =
      ((), { s with tstack := { t with vStart := curV, spans := setSides (!edgeDir) (t.spans.1 ++ t.spans.2) [] } :: rest }) := by
  rw [retarget, run_modifyCur, hs]

theorem Step.retarget (hi : s.Inv D) (hs : Shape s) (curV : Nat) (edgeDir : Bool) (hok : RetargetOk D curV s) :
    Step D v s ((retarget curV edgeDir).run s).2 := by
  match hts : s.tstack with
  | [] => exact absurd hts hok.nonempty
  | t :: rest =>
    have hc : curE s = t := by rw [curE, hts, List.head!_cons']
    have old := hok.old
    rw [hc] at old
    rw [retarget_run_eq curV edgeDir s t rest hts]
    set t' : TEntry := { t with vStart := curV, spans := setSides (!edgeDir) (t.spans.1 ++ t.spans.2) [] } with ht'
    have hE : ∀ e, t'.edges s.g s.items e ↔ t.edges s.g s.items e := by
      intro e; simp only [TEntry.edges, ht', mem_setSides]
    have ht := hi.entries t (by simp [hts])
    have hinv : Inv D { s with tstack := t' :: rest } := by
      refine ⟨fun u hu => ?_, fun j hj hsz => ItemInv.congr (s := s) rfl rfl (fun _ _ => Iff.rfl) (hi.nodes j hj hsz)⟩
      simp only [List.mem_cons] at hu
      rcases hu with rfl | hu
      · refine ⟨(Graph.ConnEdges.congr fun e _ => hE e).2 ht.conn, ?_⟩
        refine (Graph.AttachedIn.congr fun e _ => hE e).2 (ht.attached.strengthen.mono ?_)
        rintro x ⟨hx, htouch, hint⟩
        show x = curV ∨ ∃ k, t.topDepth ≤ k ∧ k ≤ D ∧ x = s.stackVerts[k]!
        rcases hx with rfl | ⟨k, h1, h2, h3⟩
        · rcases old with h | h | ⟨k, h1, h2, h3⟩
          · exact .inl h
          · exact absurd h hint
          · exact .inr ⟨k, h1, h2, h3⟩
        · exact .inr ⟨k, h1, h2, h3⟩
      · exact EntryInv.congr (s := s) rfl rfl rfl rfl (fun _ _ => Iff.rfl) (hi.entries u (by simp [hts, hu]))
    refine ⟨hinv, hs.tstack fun u hu i hi' => ?_, rfl, rfl, fun _ => Iff.rfl⟩
    simp only [List.mem_cons] at hu
    rcases hu with rfl | hu
    · rw [ht', mem_setSides] at hi'
      exact hs.span t (by simp [hts]) i hi'
    · exact hs.span u (by simp [hts, hu]) i hi'

theorem ItemFree.retarget {item : ItemId} (curV : Nat) (edgeDir : Bool) (h : ItemFree s item) :
    ItemFree ((retarget curV edgeDir).run s).2 item := by
  match hts : s.tstack with
  | [] =>
    rw [WalkState.retarget, run_modifyCur, hts]
    exact ItemFree.frame (s := s) rfl rfl (fun t ht => by simp at ht) h
  | t :: rest =>
    rw [retarget_run_eq curV edgeDir s t rest hts]
    refine ItemFree.frame (s := s) rfl rfl (fun u hu i hi => ?_) h
    simp only [List.mem_cons] at hu
    rcases hu with rfl | hu
    · rw [mem_setSides] at hi; exact ⟨t, by simp [hts], hi⟩
    · exact ⟨u, by simp [hts, hu], hi⟩

/-! ### The blocks of `finishEdge` -/

theorem loop1Type_run (d : Nat) (edgeDir : Bool) (s : WalkState) : (loop1Type d edgeDir).run s =
    if (nxtE s).topDepth > d then
      (.S, (mergeTstackTops.run { s with stackDir := s.stackDir.set! (nxtE s).topDepth edgeDir }).2)
    else if (nxtE s).vStart == (curE s).vStart then (.P, s) else (.R, s) := by
  simp only [loop1Type, WalkM.run_bind, run_nxt, nxtE, curE]
  by_cases h : s.tstack.tail.head!.topDepth > d
  · simp only [h, ↓reduceIte, WalkM.run_bind, run_nxt, run_setStackDir, WalkM.pure_run]
  · simp only [h, ↓reduceIte, WalkM.run_bind, run_nxt, run_cur]
    by_cases h' : (s.tstack.tail.head!.vStart == s.tstack.head!.vStart) = true
    · simp only [h', ↓reduceIte, WalkM.pure_run]
    · simp only [h', Bool.false_eq_true, ↓reduceIte, WalkM.pure_run]

theorem loop1Type_result (d : Nat) (edgeDir : Bool) (s : WalkState) :
    result (loop1Type d edgeDir) s ∉ [NodeType.F, .V, .Q] := by
  unfold result; rw [loop1Type_run]
  split
  · simp
  · split <;> simp

/-- The S case of loop 1 merges the entry at depth `> d` into the tree edge. -/
theorem Step.loop1Type (hi : s.Inv D) (hs : Shape s) {d : Nat} {edgeDir : Bool}
    (hok : (nxtE s).topDepth > d → MergeTopOk D { s with stackDir := s.stackDir.set! (nxtE s).topDepth edgeDir }) :
    Step D v s (after (loop1Type d edgeDir) s) := by
  unfold after; rw [loop1Type_run]
  by_cases h : (nxtE s).topDepth > d
  · simp only [h, ↓reduceIte]
    have st : Step D v s { s with stackDir := s.stackDir.set! (nxtE s).topDepth edgeDir } :=
      Step.frame hi hs rfl rfl rfl rfl
    exact st.trans (Step.mergeTop st.inv st.shape (hok h))
  · simp only [h, ↓reduceIte]; split <;> exact Step.refl hi hs

def l1Ty (d : Nat) (edgeDir : Bool) (s : WalkState) : NodeType := result (loop1Type d edgeDir) s
def l1S₁ (d : Nat) (edgeDir : Bool) (s : WalkState) : WalkState := after (loop1Type d edgeDir) s
def l1S₂ (d : Nat) (edgeDir : Bool) (s : WalkState) : WalkState :=
  after (maybeUnwrapNxt (l1Ty d edgeDir s)) (l1S₁ d edgeDir s)

/-- Loop-1 body: classify (merging in the S case), unwrap, merge, close. -/
structure Loop1BodyOk (D d : Nat) (edgeDir : Bool) (s : WalkState) : Prop where
  mergeS : (nxtE s).topDepth > d → MergeTopOk D { s with stackDir := s.stackDir.set! (nxtE s).topDepth edgeDir }
  unwrap : UnwrapOk (l1Ty d edgeDir s) (l1S₁ d edgeDir s)
  close : CloseTwoOk D (l1S₂ d edgeDir s)

theorem Step.loop1Body (hi : s.Inv D) (hs : Shape s) (hv : v < s.g.nv) {d : Nat} {edgeDir : Bool}
    (hok : Loop1BodyOk D d edgeDir s) : Step D v s ((loop1Body d edgeDir).run s).2 := by
  have st₁ : Step D v s (l1S₁ d edgeDir s) := Step.loop1Type hi hs hok.mergeS
  have hv₁ : v < (l1S₁ d edgeDir s).g.nv := by rw [st₁.g]; exact hv
  have r := maybeUnwrapNxt_spec (v := v) st₁.inv st₁.shape (loop1Type_result d edgeDir s) hok.unwrap
  have st₃ := Step.closeTwo r.step.inv r.step.shape (by rw [r.step.g]; exact hv₁) _ hok.close r.free
  exact st₁.trans (r.step.trans st₃)

def ceS₁ (nxtV d e : Nat) (s : WalkState) : WalkState := after (pushEdgeTstack nxtV d e) s

/-- `closeEars`: push the tree edge `e = {nxtV, stackVerts[d]}` and run loop 1. -/
structure CloseEarsOk (D nxtV d e : Nat) (edgeDir : Bool) (s : WalkState) : Prop where
  e_lt : e < s.g.ne
  q : Items.ch s.items (edgeItem s.g e) = []
  ends : Items.PairEq (nxtV, s.stackVerts[d]!) s.g.edges[e]!
  d_le : d ≤ D
  body : ∀ k, (∀ j, j ≤ k → result (loop1Cond d) (iter (loop1Body d edgeDir) j (ceS₁ nxtV d e s)) = true) →
    Loop1BodyOk D d edgeDir (iter (loop1Body d edgeDir) k (ceS₁ nxtV d e s))

theorem Step.closeEars (hi : s.Inv D) (hs : Shape s) (hv : v < s.g.nv) {nxtV d e : Nat} {edgeDir : Bool}
    (hok : CloseEarsOk D nxtV d e edgeDir s) : Step D v s ((closeEars nxtV d e edgeDir).run s).2 := by
  have st₁ : Step D v s (ceS₁ nxtV d e s) := Step.pushEdge hi hs nxtV d e hok.e_lt hok.q hok.ends hok.d_le
  exact st₁.trans (Step.loop (loop1Cond d) (Spqr.loop1Body d edgeDir) (Loop1BodyOk D d edgeDir) _ (fun _ => rfl)
    (fun _ hv hi hs _ hok => Step.loop1Body hi hs hv hok) (by rw [st₁.g]; exact hv) st₁.inv st₁.shape hok.body)

theorem mergeLate_run (d : Nat) (s : WalkState) : (mergeLate d).run s =
    if (curE s).firstIdx > s.firstOccurrence[d]! then
      (false, ((loop s.tstack.length (loop2Cond s.firstOccurrence[d]!) mergeTstackTops).run s).2)
    else (true, s) := by
  simp only [mergeLate, WalkM.run_bind, run_cur, WalkM.get_run, curE]
  by_cases h : s.tstack.head!.firstIdx > s.firstOccurrence[d]!
  · simp only [h, ↓reduceIte, WalkM.run_bind, run_tstackSize, WalkM.pure_run]
  · simp only [h, ↓reduceIte, WalkM.pure_run]

/-- `mergeLate`: every late merge satisfies `MergeTopOk`. -/
structure MergeLateOk (D d : Nat) (s : WalkState) : Prop where
  body : (curE s).firstIdx > s.firstOccurrence[d]! → ∀ k,
    (∀ j, j ≤ k → result (loop2Cond s.firstOccurrence[d]!) (iter mergeTstackTops j s) = true) →
    MergeTopOk D (iter mergeTstackTops k s)

theorem Step.mergeLate (hi : s.Inv D) (hs : Shape s) (hv : v < s.g.nv) {d : Nat} (hok : MergeLateOk D d s) :
    Step D v s ((mergeLate d).run s).2 := by
  rw [mergeLate_run]
  by_cases h : (curE s).firstIdx > s.firstOccurrence[d]!
  · simp only [h, ↓reduceIte]
    exact Step.loop _ mergeTstackTops (MergeTopOk D) _ (fun _ => rfl)
      (fun _ _ hi hs _ hok => Step.mergeTop hi hs hok) hv hi hs (hok.body h)
  · simp only [h, ↓reduceIte]; exact Step.refl hi hs

def vertPre (isType1 : Bool) (origTstack : Nat) (isSingle : Bool) : WalkM Bool :=
  if !isType1 then do loop (← tstackSize) (loop3Cond origTstack) mergeTstackTops; pure false else pure isSingle

def vertUnwrap (isType1 isSingle : Bool) : WalkM (Option ItemId) :=
  if isType1 then some <$> maybeUnwrapNxt (if isSingle then .S else .R) else pure none

def vertFinish (item : Option ItemId) (isSingle : Bool) : WalkM Bool :=
  match item with
  | some item => do finishTstackTop item; pure true
  | none => pure isSingle

/-- `closeVert` without its continuation. -/
def closeVert' (curV : Nat) (edgeDir isType1 : Bool) (origTstack : Nat) (isSingle : Bool) : WalkM Bool := do
  let isSingle ← vertPre isType1 origTstack isSingle
  let item ← vertUnwrap isType1 isSingle
  mergeTstackTops
  mergeTstackTops
  retarget curV edgeDir
  vertFinish item isSingle

theorem closeVert_eq (curV : Nat) (edgeDir isType1 : Bool) (origTstack : Nat) (isSingle : Bool)
    (k : Bool → WalkM Bool) :
    closeVert curV edgeDir isType1 origTstack isSingle k =
      closeVert' curV edgeDir isType1 origTstack isSingle >>= k := by
  cases isType1 <;> simp only [closeVert, closeVertTail, closeVert', vertPre, vertUnwrap, retarget, vertFinish, Bool.not_false,
    Bool.not_true, Bool.false_eq_true, ↓reduceIte, bind_assoc, pure_bind, map_eq_pure_bind]

def cvB₁ (isType1 : Bool) (origTstack : Nat) (isSingle : Bool) (s : WalkState) : Bool :=
  result (vertPre isType1 origTstack isSingle) s
def cvS₁ (isType1 : Bool) (origTstack : Nat) (isSingle : Bool) (s : WalkState) : WalkState :=
  after (vertPre isType1 origTstack isSingle) s
def cvS₂ (isType1 : Bool) (origTstack : Nat) (isSingle : Bool) (s : WalkState) : WalkState :=
  after (vertUnwrap isType1 (cvB₁ isType1 origTstack isSingle s)) (cvS₁ isType1 origTstack isSingle s)
def cvS₃ (isType1 : Bool) (origTstack : Nat) (isSingle : Bool) (s : WalkState) : WalkState :=
  after mergeTstackTops (cvS₂ isType1 origTstack isSingle s)
def cvS₄ (isType1 : Bool) (origTstack : Nat) (isSingle : Bool) (s : WalkState) : WalkState :=
  after mergeTstackTops (cvS₃ isType1 origTstack isSingle s)
def cvS₅ (curV : Nat) (edgeDir isType1 : Bool) (origTstack : Nat) (isSingle : Bool) (s : WalkState) : WalkState :=
  after (retarget curV edgeDir) (cvS₄ isType1 origTstack isSingle s)

/-- `closeVert`: loop 3 (type 2), the unwrap (type 1), the two merges into the vertex entry, the
re-targeting to `curV`, and the close (type 1). -/
structure CloseVertOk (D curV : Nat) (edgeDir isType1 : Bool) (origTstack : Nat) (isSingle : Bool)
    (s : WalkState) : Prop where
  loop3 : isType1 = false → ∀ k,
    (∀ j, j ≤ k → result (loop3Cond origTstack) (iter mergeTstackTops j s) = true) →
    MergeTopOk D (iter mergeTstackTops k s)
  unwrap : isType1 = true → UnwrapOk (if isSingle then .S else .R) s
  merge₁ : MergeTopOk D (cvS₂ isType1 origTstack isSingle s)
  merge₂ : MergeTopOk D (cvS₃ isType1 origTstack isSingle s)
  retarget : RetargetOk D curV (cvS₄ isType1 origTstack isSingle s)
  finish : isType1 = true → FinishTopOk D (cvS₅ curV edgeDir isType1 origTstack isSingle s)

theorem Step.vertPre (hi : s.Inv D) (hs : Shape s) (hv : v < s.g.nv) {isType1 : Bool} {origTstack : Nat}
    {isSingle : Bool}
    (hok : isType1 = false → ∀ k,
      (∀ j, j ≤ k → result (loop3Cond origTstack) (iter mergeTstackTops j s) = true) →
      MergeTopOk D (iter mergeTstackTops k s)) :
    Step D v s ((vertPre isType1 origTstack isSingle).run s).2 := by
  cases isType1
  · simp only [WalkState.vertPre, Bool.not_false, ↓reduceIte, WalkM.run_bind, run_tstackSize, WalkM.pure_run]
    exact Step.loop _ mergeTstackTops (MergeTopOk D) _ (fun _ => rfl)
      (fun _ _ hi hs _ hok => Step.mergeTop hi hs hok) hv hi hs (hok rfl)
  · exact Step.refl hi hs

theorem Step.vertUnwrap (hi : s.Inv D) (hs : Shape s) {isType1 isSingle : Bool}
    (hok : isType1 = true → UnwrapOk (if isSingle then .S else .R) s) :
    Step D v s ((vertUnwrap isType1 isSingle).run s).2 ∧
      ∀ i, ((vertUnwrap isType1 isSingle).run s).1 = some i → ItemFree ((vertUnwrap isType1 isSingle).run s).2 i := by
  cases isType1
  · exact ⟨Step.refl hi hs, fun i h => by simp [WalkState.vertUnwrap] at h⟩
  · have r := maybeUnwrapNxt_spec (v := v) hi hs (ty := if isSingle then .S else .R) (by cases isSingle <;> decide) (hok rfl)
    simp only [WalkState.vertUnwrap, ↓reduceIte, WalkM.map_run, Option.some.injEq]
    exact ⟨r.step, fun i h => h ▸ r.free⟩

theorem Step.closeVert' (hi : s.Inv D) (hs : Shape s) (hv : v < s.g.nv) {curV : Nat} {edgeDir isType1 : Bool}
    {origTstack : Nat} {isSingle : Bool} (hok : CloseVertOk D curV edgeDir isType1 origTstack isSingle s) :
    Step D v s ((closeVert' curV edgeDir isType1 origTstack isSingle).run s).2 := by
  have st₁ : Step D v s (cvS₁ isType1 origTstack isSingle s) := Step.vertPre hi hs hv hok.loop3
  have hv₁ : v < (cvS₁ isType1 origTstack isSingle s).g.nv := by rw [st₁.g]; exact hv
  obtain ⟨st₂, hfree⟩ := Step.vertUnwrap (v := v) st₁.inv st₁.shape (isType1 := isType1)
    (isSingle := cvB₁ isType1 origTstack isSingle s)
    (fun h => by subst h; exact hok.unwrap rfl)
  have st₂ : Step D v (cvS₁ isType1 origTstack isSingle s) (cvS₂ isType1 origTstack isSingle s) := st₂
  have st₃ : Step D v _ (cvS₃ isType1 origTstack isSingle s) := Step.mergeTop st₂.inv st₂.shape hok.merge₁
  have st₄ : Step D v _ (cvS₄ isType1 origTstack isSingle s) := Step.mergeTop st₃.inv st₃.shape hok.merge₂
  have st₅ : Step D v _ (cvS₅ curV edgeDir isType1 origTstack isSingle s) :=
    Step.retarget st₄.inv st₄.shape curV edgeDir hok.retarget
  have hv₅ : v < (cvS₅ curV edgeDir isType1 origTstack isSingle s).g.nv := by
    rw [st₅.g, st₄.g, st₃.g, st₂.g]; exact hv₁
  have st := st₁.trans (st₂.trans (st₃.trans (st₄.trans st₅)))
  cases isType1
  · exact st
  · have hf : ItemFree (cvS₂ true origTstack isSingle s) ((maybeUnwrapNxt (if isSingle then .S else .R)).run s).1 :=
      hfree _ rfl
    have hf₅ := ((hf.merge).merge).retarget curV edgeDir
    exact st.trans (Step.finishTop st₅.inv st₅.shape hv₅ _ (hok.finish rfl) hf₅)

/-- `finishP`: the type-1 P-check closes `nxt` as a P-node with the back edge. -/
structure FinishPOk (D curV lowval : Nat) (isType1 : Bool) (s : WalkState) : Prop where
  ok : result (condP curV lowval isType1) s = true →
    UnwrapOk .P s ∧ CloseTwoOk D (after (maybeUnwrapNxt .P) s)

theorem Step.finishP (hi : s.Inv D) (hs : Shape s) (hv : v < s.g.nv) {curV lowval : Nat} {isType1 : Bool}
    (hok : FinishPOk D curV lowval isType1 s) : Step D v s ((finishP curV lowval isType1).run s).2 := by
  simp only [Spqr.finishP, WalkM.run_bind, run_condP]
  by_cases h : result (condP curV lowval isType1) s = true
  · have h' : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == curV) &&
        (s.tstack.tail.head!.topDepth == lowval)) = true := h
    simp only [h', ↓reduceIte, WalkM.run_bind]
    obtain ⟨hu, hc⟩ := hok.ok h
    have r := maybeUnwrapNxt_spec (v := v) hi hs (by decide) hu
    exact r.step.trans (Step.closeTwo r.step.inv r.step.shape (by rw [r.step.g]; exact hv) _ hc r.free)
  · have h' : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == curV) &&
        (s.tstack.tail.head!.topDepth == lowval)) = false := Bool.eq_false_iff.2 h
    simp only [h', Bool.false_eq_true, ↓reduceIte]
    exact Step.refl hi hs

/-- `finishTail`: the first-edge vertex push needs the blocks already hanging at `curV` to be a
connected piece attached only at `curV`; the merge into the first ear needs `MergeTopOk`. -/
structure FinishTailOk (D curV d : Nat) (hasVert isSingle : Bool) (s : WalkState) : Prop where
  conn : hasVert = false → s.g.ConnEdges (Items.EdgeBelow s.g s.items (vertItem curV))
  att : hasVert = false → s.g.TwoAttached (Items.EdgeBelow s.g s.items (vertItem curV)) curV curV
  merge : hasVert = false → isSingle = false → MergeTopOk D (after (pushVertTstack curV d) s)

theorem Step.finishTail (hi : s.Inv D) (hs : Shape s) {curV d : Nat} (hv : curV < s.g.nv) {hasVert isSingle : Bool}
    (hok : FinishTailOk D curV d hasVert isSingle s) :
    Step D curV s ((finishTail curV d hasVert isSingle).run s).2 := by
  cases hasVert
  · have st₁ := Step.pushVert hi hs d hv (hok.conn rfl) (hok.att rfl)
    cases isSingle
    · exact st₁.trans (Step.mergeTop st₁.inv st₁.shape (hok.merge rfl rfl))
    · exact st₁
  · exact Step.refl hi hs

/-- `finishRest = finishP; finishTail`. -/
structure FinishRestOk (D curV d lowval : Nat) (isType1 hasVert isSingle : Bool) (s : WalkState) : Prop where
  p : FinishPOk D curV lowval isType1 s
  tail : FinishTailOk D curV d hasVert isSingle (after (finishP curV lowval isType1) s)

theorem Step.finishRest (hi : s.Inv D) (hs : Shape s) {curV d lowval : Nat} (hv : curV < s.g.nv)
    {isType1 hasVert isSingle : Bool} (hok : FinishRestOk D curV d lowval isType1 hasVert isSingle s) :
    Step D curV s ((finishRest curV d lowval isType1 hasVert isSingle).run s).2 := by
  have st₁ := Step.finishP (v := curV) hi hs hv hok.p
  exact st₁.trans (Step.finishTail st₁.inv st₁.shape (by rw [st₁.g]; exact hv) hok.tail)

/-! ### `finishEdge` -/

def feS₀ (d : Nat) (o : DfsOut) (s : WalkState) : WalkState :=
  after (modifyItem (edgeItem s.g o.e) fun it =>
    { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) s
def feS₁ (d : Nat) (o : DfsOut) (s : WalkState) : WalkState :=
  after (closeEars o.dest d o.e s.stackDir[d]!) (feS₀ d o s)
def feSingle (d : Nat) (o : DfsOut) (s : WalkState) : Bool := result (mergeLate d) (feS₁ d o s)
def feS₂ (d : Nat) (o : DfsOut) (s : WalkState) : WalkState := after (mergeLate d) (feS₁ d o s)
def feS₃ (curV d : Nat) (o : DfsOut) (origTstack : Nat) (s : WalkState) : WalkState :=
  after (closeVert' curV s.stackDir[d]! o.cls.isType1 origTstack (feSingle d o s)) (feS₂ d o s)
def feB₃ (curV d : Nat) (o : DfsOut) (origTstack : Nat) (s : WalkState) : Bool :=
  result (closeVert' curV s.stackDir[d]! o.cls.isType1 origTstack (feSingle d o s)) (feS₂ d o s)
def feBack (curV lv d : Nat) (o : DfsOut) (s : WalkState) : WalkState :=
  after (modify fun s => { s with firstOccurrence := s.firstOccurrence.modify lv (min · s.nxtEdgeIdx),
                                  nxtEdgeIdx := s.nxtEdgeIdx + 1 })
    (after (pushEdgeTstack curV lv o.e) (feS₀ d o s))

/-- The stack-shape hypotheses of `finishEdge` for an out-edge `o.cls = .ret lv kind`, `lv < d`, block
by block (each at the state where the block runs). -/
structure FinishOk (D curV d lv : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool) (s : WalkState) : Prop where
  e_lt : o.e < s.g.ne
  ears : o.cls.isTree = true → CloseEarsOk D o.dest d o.e s.stackDir[d]! (feS₀ d o s)
  late : o.cls.isTree = true → MergeLateOk D d (feS₁ d o s)
  vert : o.cls.isTree = true → hasVert = true →
    CloseVertOk D curV s.stackDir[d]! o.cls.isType1 origTstack (feSingle d o s) (feS₂ d o s)
  rest_vert : o.cls.isTree = true → hasVert = true →
    FinishRestOk D curV d lv o.cls.isType1 hasVert (feB₃ curV d o origTstack s) (feS₃ curV d o origTstack s)
  rest_tree : o.cls.isTree = true → hasVert = false →
    FinishRestOk D curV d lv o.cls.isType1 hasVert (feSingle d o s) (feS₂ d o s)
  q : o.cls.isTree = false → Items.ch s.items (edgeItem s.g o.e) = []
  ends : o.cls.isTree = false → Items.PairEq (curV, s.stackVerts[lv]!) s.g.edges[o.e]!
  lv_le : o.cls.isTree = false → lv ≤ D
  rest_back : o.cls.isTree = false →
    FinishRestOk D curV d lv o.cls.isType1 hasVert true (feBack curV lv d o s)

/-- `finishEdge` preserves the invariant on a returning edge, under `FinishOk`. -/
theorem finishEdge_inv (curV d lv : Nat) (kind : RetKind) (o : DfsOut) (origTstack : Nat) (hasVert : Bool)
    (ho : o.cls = .ret lv kind) (hlow : lv < d) (hv : curV < s.g.nv) (hi : s.Inv D) (hs : Shape s)
    (hok : FinishOk D curV d lv o origTstack hasVert s) :
    Step D curV s ((finishEdge curV d o origTstack hasVert).run s).2 := by
  have hlv : o.cls.lowval d = lv := by rw [ho]; rfl
  have hge : ¬ (lv ≥ d) := by omega
  rw [finishEdge_eq]
  simp only [finishEdge', finishTree, finishBack, hlv, hge, ↓reduceIte, WalkM.run_bind, WalkM.get_run,
    run_stackDir, run_makeVs, run_modifyItem]
  have st₀ : Step D curV s (feS₀ d o s) :=
    Step.modifyVs hi hs (edgeItem s.g o.e) _ (by show 1 + s.g.nv + o.e < _; have := hok.e_lt; omega)
  have hv₀ : curV < (feS₀ d o s).g.nv := by rw [st₀.g]; exact hv
  by_cases ht : o.cls.isTree = true
  · simp only [ht, ↓reduceIte, closeVert_eq, WalkM.run_bind]
    have st₁ : Step D curV _ (feS₁ d o s) := Step.closeEars st₀.inv st₀.shape hv₀ (hok.ears ht)
    have hv₁ : curV < (feS₁ d o s).g.nv := by rw [st₁.g]; exact hv₀
    have st₂ : Step D curV _ (feS₂ d o s) := Step.mergeLate st₁.inv st₁.shape hv₁ (hok.late ht)
    have hv₂ : curV < (feS₂ d o s).g.nv := by rw [st₂.g]; exact hv₁
    cases hasVert
    · simp only [Bool.false_eq_true, ↓reduceIte]
      exact st₀.trans (st₁.trans (st₂.trans (Step.finishRest st₂.inv st₂.shape hv₂ (hok.rest_tree ht rfl))))
    · simp only [↓reduceIte, WalkM.run_bind]
      have st₃ : Step D curV _ (feS₃ curV d o origTstack s) :=
        Step.closeVert' st₂.inv st₂.shape hv₂ (hok.vert ht rfl)
      have hv₃ : curV < (feS₃ curV d o origTstack s).g.nv := by rw [st₃.g]; exact hv₂
      exact st₀.trans (st₁.trans (st₂.trans (st₃.trans
        (Step.finishRest st₃.inv st₃.shape hv₃ (hok.rest_vert ht rfl)))))
  · have ht' : o.cls.isTree = false := Bool.eq_false_iff.2 ht
    simp only [ht', Bool.false_eq_true, ↓reduceIte, WalkM.run_bind]
    have hq : Items.ch (feS₀ d o s).items (edgeItem (feS₀ d o s).g o.e) = [] := by
      show Items.ch (s.items.modify _ _) (edgeItem s.g o.e) = []
      rw [Items.ch_modify_ch_eq (edgeItem s.g o.e)
        (fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) (fun _ => rfl)]
      exact hok.q ht'
    have st₁ : Step D curV _ (after (pushEdgeTstack curV lv o.e) (feS₀ d o s)) :=
      Step.pushEdge st₀.inv st₀.shape curV lv o.e hok.e_lt hq (hok.ends ht') (hok.lv_le ht')
    have st₂ : Step D curV _ (feBack curV lv d o s) := Step.frame st₁.inv st₁.shape rfl rfl rfl rfl
    have hv₂ : curV < (feBack curV lv d o s).g.nv := by rw [st₂.g, st₁.g]; exact hv₀
    exact st₀.trans (st₁.trans (st₂.trans (Step.finishRest st₂.inv st₂.shape hv₂ (hok.rest_back ht'))))

/-- The back-edge case: push the back edge `o.e = {curV, stackVerts[lv]}` at depth `lv`, then the
type-1 P-check and the vertex push. -/
theorem finishEdge_back_inv (curV d lv : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool)
    (ho : o.cls = .ret lv .backEdge) (hlow : lv < d) (hv : curV < s.g.nv) (hi : s.Inv D) (hs : Shape s)
    (he : o.e < s.g.ne) (hq : Items.ch s.items (edgeItem s.g o.e) = [])
    (hend : Items.PairEq (curV, s.stackVerts[lv]!) s.g.edges[o.e]!) (hD : lv ≤ D)
    (hrest : FinishRestOk D curV d lv true hasVert true (feBack curV lv d o s)) :
    ((finishEdge curV d o origTstack hasVert).run s).2.Inv D := by
  have ht : o.cls.isTree = false := by rw [ho]; rfl
  have h1 : o.cls.isType1 = true := by rw [ho]; rfl
  refine (finishEdge_inv curV d lv .backEdge o origTstack hasVert ho hlow hv hi hs
    ⟨he, ?_, ?_, ?_, ?_, ?_, fun _ => hq, fun _ => hend, fun _ => hD, fun _ => ?_⟩).inv
  all_goals first | exact fun h => absurd h (by rw [ht]; decide) | (rw [h1]; exact hrest)

end WalkState

/-- The whole walk preserves the invariant (the ear frame rule fixes the depth bound). -/
theorem walkTree_inv {D : Nat} (t : DfsTree) (d : Nat) (s : WalkState) (h : s.Inv D) :
    ((walkTree t d).run s).2.Inv D := by
  sorry

/-- Completeness: after the walk, the allocated nodes partition the edges (every edge is below
exactly one child chain from the root) — the `Items.Tree` content of `Items.WF`. -/
theorem walk_nodes_partition (g : Graph) (tern : Bool) (forest : List DfsTree) :
    let s := g.walk tern forest
    ∀ e, e < g.ne → ∃ p, Items.IsParent s.items p (edgeItem g e) ∧
      ∀ p', Items.IsParent s.items p' (edgeItem g e) → p' = p := by
  sorry

end Spqr
