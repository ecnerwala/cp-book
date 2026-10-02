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


end WalkState

/-- `finishEdge` preserves the invariant. -/
theorem finishEdge_inv {D : Nat} (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool)
    (s : WalkState) (h : s.Inv D) :
    ((finishEdge curV d o origTstack hasVert).run s).2.Inv D := by
  sorry

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
