import Spqr.WalkSpec

/-!
# Corrected attachment invariant for open entries (PROOF.md §4.2b correction)

`WalkSpec.Inv D` attaches every open entry within its own terminals `TEntry.Term D`. That is false
mid-walk: after `finishEdge` of a type-2 chain frame an entry can stay open below the entries of
the chain and remain attached at a vertex that is the `vStart` of an entry *above* it, and that
vertex later drops out of `stackVerts` (cycle `0..6`, chords `6-1`, `5-2`, ear `4-7-3`: the entry
`(6,5,[],[Q(5,6)])` stays attached at `5` while `stackVerts[5]` becomes `7`).

The corrected attachment set `TEntry.Term'` adds the terminals of the entries above; `Stack` reads
the stack top-down carrying that context, and `Inv'` is `Inv` with `Stack` in place of the
per-entry `EntryInv`. The context is empty at the top, so the top entry satisfies `EntryInv D`
exactly as before and every close (`finishTstackTop` acts on the top) still records a strictly
2-attached item.

This file proves the primitive steps for `Inv'` (`Step'`): the item-only steps need nothing new;
merging and re-targeting need, in addition to `MergeOk`/`RetargetOk`, that the entries below are
edge-disjoint from the touched ones (the `Interior` case of `MergeOk.bottom`/`RetargetOk.old` is
otherwise not transferable to an entry attached there); popping needs that no remaining entry is
attached at a lost terminal unless it is its own.
-/

namespace Spqr
open WalkM

namespace TEntry

variable {D D' : Nat} {s : WalkState} {t : TEntry} {v : Nat}

theorem Term.mono (hD : D ≤ D') (h : t.Term D s v) : t.Term D' s v := by
  rcases h with hv | ⟨k, h1, h2, h3⟩
  · exact .inl hv
  · exact .inr ⟨k, h1, Nat.le_trans h2 hD, h3⟩

/-- `v` is a terminal of `t` or of an entry above `t`. -/
def Term' (D : Nat) (s : WalkState) (above : List TEntry) (t : TEntry) (v : Nat) : Prop :=
  t.Term D s v ∨ ∃ t' ∈ above, t'.Term D s v

theorem Term'_of_Term {above : List TEntry} (h : t.Term D s v) : t.Term' D s above v := .inl h

theorem Term'_nil : t.Term' D s [] v ↔ t.Term D s v := by simp [Term']

end TEntry

namespace WalkState

variable {D D' : Nat} {s s' : WalkState}

/-- `EntryInv` with the attachment set widened by the terminals of the entries `above`. -/
structure EntryInv' (D : Nat) (s : WalkState) (above : List TEntry) (t : TEntry) : Prop where
  conn : s.g.ConnEdges (t.edges s.g s.items)
  attached : s.g.AttachedIn (t.edges s.g s.items) (t.Term' D s above)

/-- The entries `l` of the stack, read top-down, below the entries `above` (nearest first). -/
def Stack (D : Nat) (s : WalkState) : List TEntry → List TEntry → Prop
  | _, [] => True
  | above, t :: rest => s.EntryInv' D above t ∧ s.Stack D (t :: above) rest

structure Inv' (D : Nat) (s : WalkState) : Prop where
  entries : s.Stack D [] s.tstack
  nodes : ∀ i, 1 + s.g.nv + s.g.ne ≤ i → i < s.items.size → s.ItemInv i

theorem Stack_nil {A : List TEntry} : s.Stack D A [] := trivial

theorem Stack_cons {A rest : List TEntry} {t : TEntry} :
    s.Stack D A (t :: rest) ↔ s.EntryInv' D A t ∧ s.Stack D (t :: A) rest := Iff.rfl

theorem EntryInv'.of_entryInv {above : List TEntry} {t : TEntry} (h : s.EntryInv D t) :
    s.EntryInv' D above t :=
  ⟨h.conn, h.attached.mono fun _ hv => .inl hv⟩

theorem EntryInv'.toEntryInv {t : TEntry} (h : s.EntryInv' D [] t) : s.EntryInv D t :=
  ⟨h.conn, h.attached.mono fun v hv => by
    rcases hv with hv | ⟨_, h, _⟩
    · exact hv
    · simp at h⟩

theorem EntryInv'.congr {A : List TEntry} {t : TEntry} (hg : s'.g = s.g) (hsv : s'.stackVerts = s.stackVerts)
    (hE : ∀ e, e < s.g.ne → (t.edges s.g s'.items e ↔ t.edges s.g s.items e))
    (h : s.EntryInv' D A t) : s'.EntryInv' D A t := by
  refine ⟨?_, ?_⟩ <;> rw [hg]
  · exact (Graph.ConnEdges.congr hE).2 h.conn
  · refine (Graph.AttachedIn.congr hE).2 (h.attached.mono fun v hv => ?_)
    simp only [TEntry.Term', TEntry.Term, hsv] at hv ⊢
    exact hv

theorem EntryInv'.mono {A : List TEntry} {t : TEntry} (hD : D ≤ D') (h : s.EntryInv' D A t) :
    s.EntryInv' D' A t :=
  ⟨h.conn, h.attached.mono fun _ hv => hv.elim (fun h => .inl (h.mono hD))
    fun ⟨t', h1, h2⟩ => .inr ⟨t', h1, h2.mono hD⟩⟩

/-- Replacing the context `A` by `A'`: every terminal of `A` at which an entry of `l` is attached
is a terminal of that entry or of `A'`. -/
theorem Stack.transfer {A A' l : List TEntry}
    (h : ∀ t ∈ l, ∀ v, (∃ t' ∈ A, t'.Term D s v) → s.g.Touches (t.edges s.g s.items) v →
      t.Term D s v ∨ ∃ t' ∈ A', t'.Term D s v)
    (hs : s.Stack D A l) : s.Stack D A' l := by
  induction l generalizing A A' with
  | nil => exact Stack_nil
  | cons t rest ih =>
    obtain ⟨ht, hrest⟩ := Stack_cons.1 hs
    refine Stack_cons.2 ⟨⟨ht.conn, ht.attached.strengthen.mono ?_⟩, ih (A := t :: A) (A' := t :: A') ?_ hrest⟩
    · rintro v ⟨hv, htouch, -⟩
      rcases hv with hv | hv
      · exact .inl hv
      · exact h t (by simp) v hv htouch
    · intro u hu v hv htouch
      obtain ⟨t', ht', hT⟩ := hv
      simp only [List.mem_cons] at ht'
      rcases ht' with rfl | ht'
      · exact .inr ⟨t', by simp, hT⟩
      · rcases h u (by simp [hu]) v ⟨t', ht', hT⟩ htouch with h' | ⟨t'', h1, h2⟩
        · exact .inl h'
        · exact .inr ⟨t'', by simp [h1], h2⟩

theorem Stack.mono_above {A A' l : List TEntry} (h : ∀ t' ∈ A, t' ∈ A') (hs : s.Stack D A l) :
    s.Stack D A' l :=
  hs.transfer fun _ _ _ ⟨t', ht', hT⟩ _ => .inr ⟨t', h t' ht', hT⟩

theorem Stack.monoD {A l : List TEntry} (hD : D ≤ D') (hs : s.Stack D A l) : s.Stack D' A l := by
  induction l generalizing A with
  | nil => exact Stack_nil
  | cons t rest ih =>
    obtain ⟨ht, hrest⟩ := Stack_cons.1 hs
    exact Stack_cons.2 ⟨ht.mono hD, ih hrest⟩

theorem Stack.congr {A l : List TEntry} (hg : s'.g = s.g) (hsv : s'.stackVerts = s.stackVerts)
    (hE : ∀ t ∈ l, ∀ e, e < s.g.ne → (t.edges s.g s'.items e ↔ t.edges s.g s.items e))
    (hs : s.Stack D A l) : s'.Stack D A l := by
  induction l generalizing A with
  | nil => exact Stack_nil
  | cons t rest ih =>
    obtain ⟨ht, hrest⟩ := Stack_cons.1 hs
    exact Stack_cons.2 ⟨ht.congr hg hsv (hE t (by simp)), ih (fun u hu => hE u (by simp [hu])) hrest⟩

theorem Stack.of_entryInv {A : List TEntry} : ∀ {l : List TEntry}, (∀ t ∈ l, s.EntryInv D t) → s.Stack D A l
  | [], _ => Stack_nil
  | t :: rest, h => Stack_cons.2 ⟨.of_entryInv (h t (by simp)), Stack.of_entryInv fun u hu => h u (by simp [hu])⟩

theorem Inv'.of_inv (h : s.Inv D) : s.Inv' D := ⟨Stack.of_entryInv h.entries, h.nodes⟩

theorem Inv'.mono (h : s.Inv' D) (hD : D ≤ D') : s.Inv' D' := ⟨h.entries.monoD hD, h.nodes⟩

/-- The top entry satisfies the original `EntryInv`. -/
theorem Inv'.top {t : TEntry} {rest : List TEntry} (h : s.Inv' D) (hs : s.tstack = t :: rest) :
    s.EntryInv D t := by
  have := h.entries; rw [hs] at this
  exact (Stack_cons.1 this).1.toEntryInv

theorem Inv'.frame (hg : s'.g = s.g) (hsv : s'.stackVerts = s.stackVerts) (hi : s'.items = s.items)
    (hts : s'.tstack = s.tstack) (h : s.Inv' D) : s'.Inv' D := by
  refine ⟨by rw [hts]; exact h.entries.congr hg hsv (fun _ _ _ _ => by rw [hi]),
    fun i hi' hsz => ItemInv.congr hg (by rw [hi]) (fun _ _ => by rw [hi]) (h.nodes i ?_ ?_)⟩
  · rwa [hg] at hi'
  · rwa [hi] at hsz

/-! ### Item-only steps -/

theorem Inv'.alloc (ty : NodeType) (h : s.Inv' D) (hsize : 1 + s.g.nv + s.g.ne ≤ s.items.size) :
    Inv' D { s with items := s.items.push ⟨ty, (none, none), []⟩ } := by
  have hB : ∀ a i, Items.Below (s.items.push ⟨ty, (none, none), []⟩) a i ↔ Items.Below s.items a i :=
    fun _ _ => Items.Below_push_nil _ rfl
  refine ⟨h.entries.congr (s := s) rfl rfl (fun _ _ e _ => TEntry.edges_congr (fun i _ e => hB i _) e),
    fun i hi hsz => ?_⟩
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

theorem Inv'.modifyVs (j : ItemId) (vsv : Option Nat × Option Nat) (h : s.Inv' D) (hj : j < 1 + s.g.nv + s.g.ne) :
    Inv' D { s with items := s.items.modify j fun it => { it with vs := vsv } } := by
  have hB : ∀ a i, Items.Below (s.items.modify j fun it => { it with vs := vsv }) a i ↔ Items.Below s.items a i :=
    fun _ _ => Items.Below_modify_ch_eq j (fun it => { it with vs := vsv }) fun _ => rfl
  refine ⟨h.entries.congr (s := s) rfl rfl (fun _ _ e _ => TEntry.edges_congr (fun i _ e => hB i _) e),
    fun i hi hsz => ?_⟩
  have hsz' : i < s.items.size := by simpa using hsz
  have hne : i ≠ j := by intro h; subst h; exact absurd hi (Nat.not_le.2 hj)
  exact ItemInv.congr (s := s) rfl (Items.vs_modify_of_ne _ _ hne) (fun e _ => hB i _) (h.nodes i hi hsz')

/-! ### Stack steps -/

theorem Inv'.pushVert (v d : Nat) (h : s.Inv' D)
    (hc : s.g.ConnEdges (Items.EdgeBelow s.g s.items (vertItem v)))
    (ha : s.g.TwoAttached (Items.EdgeBelow s.g s.items (vertItem v)) v v) :
    ((pushVertTstack v d).run s).2.Inv' D := by
  rw [pushVertTstack, run_pushTstack]
  refine ⟨Stack_cons.2 ⟨?_, (h.entries.mono_above (by simp)).congr (s := s) rfl rfl (fun _ _ _ _ => Iff.rfl)⟩,
    fun i hi hsz => ItemInv.congr (s := s) rfl rfl (fun _ _ => Iff.rfl) (h.nodes i hi hsz)⟩
  have hE : ∀ e, TEntry.edges s.g s.items ⟨v, d, s.nxtEdgeIdx, setSides s.stackDir[d]! [vertItem v] []⟩ e ↔
      Items.EdgeBelow s.g s.items (vertItem v) e := by
    intro e; simp only [TEntry.edges, mem_setSides, List.mem_singleton, exists_eq_left]
  refine ⟨(Graph.ConnEdges.congr fun e _ => hE e).2 hc, ?_⟩
  refine ((Graph.AttachedIn.congr fun e _ => hE e).2 (Graph.twoAttached_iff.1 ha)).mono fun x hx => ?_
  rcases hx with rfl | rfl <;> exact .inl (.inl rfl)

theorem Inv'.pushEdge (vStart topDepth e : Nat) (h : s.Inv' D) (_he : e < s.g.ne)
    (hq : Items.ch s.items (edgeItem s.g e) = [])
    (hend : Items.PairEq (vStart, s.stackVerts[topDepth]!) s.g.edges[e]!) (hD : topDepth ≤ D) :
    ((pushEdgeTstack vStart topDepth e).run s).2.Inv' D := by
  rw [run_pushEdgeTstack]
  refine ⟨Stack_cons.2 ⟨?_, (h.entries.mono_above (by simp)).congr (s := s) rfl rfl (fun _ _ _ _ => Iff.rfl)⟩,
    fun i hi hsz => ItemInv.congr (s := s) rfl rfl (fun _ _ => Iff.rfl) (h.nodes i hi hsz)⟩
  have hE := TEntry.edges_edgeEntry (items := s.items) s.stackDir[topDepth]! vStart topDepth s.nxtEdgeIdx e hq
  refine ⟨Graph.ConnEdges.single hE,
    (Graph.twoAttached_iff.1 (Graph.TwoAttached.single_of_pairEq hE hend)).mono fun x hx => .inl ?_⟩
  rcases hx with rfl | rfl
  · exact .inl rfl
  · exact .inr ⟨topDepth, Nat.le_refl _, hD, rfl⟩

/-- `mergeTstackTops` under `MergeOk`, when the entries below are edge-disjoint from the two merged
ones (so none of them is attached at a vertex interior to the union). -/
theorem Inv'.mergeTop (cur nxt : TEntry) (rest : List TEntry) (hs : s.tstack = cur :: nxt :: rest)
    (h : s.Inv' D) (hok : MergeOk D s cur nxt)
    (hdisj : ∀ t ∈ rest, ∀ e, e < s.g.ne → t.edges s.g s.items e →
      ¬ (cur.edges s.g s.items e ∨ nxt.edges s.g s.items e)) :
    (mergeTstackTops.run s).2.Inv' D := by
  rw [mergeTstackTops_run_eq s cur nxt rest hs]
  have hst := h.entries; rw [hs] at hst
  obtain ⟨hc, hn, hrest⟩ := (Stack_cons.1 hst).imp_right Stack_cons.1
  have hc' : s.EntryInv D cur := hc.toEntryInv
  have hcur : ∀ v, ¬ s.g.Interior (fun e => cur.edges s.g s.items e ∨ nxt.edges s.g s.items e) v →
      cur.Term D s v → (TEntry.mergeInto cur nxt).Term D s v := by
    rintro v hint (rfl | ⟨k, h1, h2, h3⟩)
    · rcases hok.bottom with h | h
      · exact h
      · exact absurd h hint
    · exact .inr ⟨k, Nat.le_trans (Nat.min_le_right _ _) h1, h2, h3⟩
  have hnxt : ∀ v, nxt.Term D s v → (TEntry.mergeInto cur nxt).Term D s v := by
    rintro v (hv | ⟨k, h1, h2, h3⟩)
    · exact .inl hv
    · exact .inr ⟨k, Nat.le_trans (Nat.min_le_left _ _) h1, h2, h3⟩
  refine ⟨Stack_cons.2 ⟨?_, ?_⟩, fun i hi hsz => ItemInv.congr (s := s) rfl rfl (fun _ _ => Iff.rfl) (h.nodes i hi hsz)⟩
  · have hm := TEntry.edges_mergeInto (g := s.g) (items := s.items) cur nxt
    refine EntryInv'.congr (s := s) rfl rfl (fun _ _ => Iff.rfl) ⟨?_, ?_⟩
    · exact (Graph.ConnEdges.congr fun e _ => hm e).2 (Graph.ConnEdges.union hc'.conn hn.conn hok.share)
    · refine (Graph.AttachedIn.congr fun e _ => hm e).2 ((Graph.AttachedIn.union hc'.attached hn.attached).mono ?_)
      rintro v ⟨hv, hint⟩
      refine .inl ?_
      rcases hv with hv | hv | ⟨t', ht', hT⟩
      · exact hcur v hint hv
      · exact hnxt v hv
      · simp only [List.mem_singleton] at ht'; subst ht'
        exact hcur v hint hT
  · refine (hrest.transfer ?_).congr (s := s) rfl rfl (fun _ _ _ _ => Iff.rfl)
    intro t ht v hv htouch
    refine .inr ⟨TEntry.mergeInto cur nxt, by simp, ?_⟩
    obtain ⟨t', ht', hT⟩ := hv
    simp only [List.mem_cons, List.not_mem_nil, or_false] at ht'
    rcases ht' with rfl | rfl
    · exact hnxt v hT
    · refine hcur v (fun hint => ?_) hT
      obtain ⟨e, he, hte, hinc⟩ := htouch
      exact hdisj t ht e he hte (hint e he hinc)

/-- `retarget` under `RetargetOk.old` for the top entry `t`, when the entries below are
edge-disjoint from `t`. -/
theorem Inv'.retarget (curV : Nat) (edgeDir : Bool) (t : TEntry) (rest : List TEntry)
    (hs : s.tstack = t :: rest) (h : s.Inv' D)
    (old : t.vStart = curV ∨ s.g.Interior (t.edges s.g s.items) t.vStart ∨
      ∃ k, t.topDepth ≤ k ∧ k ≤ D ∧ t.vStart = s.stackVerts[k]!)
    (hdisj : ∀ u ∈ rest, ∀ e, e < s.g.ne → u.edges s.g s.items e → ¬ t.edges s.g s.items e) :
    ((retarget curV edgeDir).run s).2.Inv' D := by
  rw [retarget_run_eq curV edgeDir s t rest hs]
  set t' : TEntry := { t with vStart := curV, spans := setSides (!edgeDir) (t.spans.1 ++ t.spans.2) [] } with ht'
  have hE : ∀ e, t'.edges s.g s.items e ↔ t.edges s.g s.items e := by
    intro e; simp only [TEntry.edges, ht', mem_setSides]
  have hst := h.entries; rw [hs] at hst
  obtain ⟨ht, hrest⟩ := Stack_cons.1 hst
  have hT : ∀ v, ¬ s.g.Interior (t.edges s.g s.items) v → t.Term D s v → t'.Term D s v := by
    rintro v hint (rfl | ⟨k, h1, h2, h3⟩)
    · rcases old with h | h | ⟨k, h1, h2, h3⟩
      · exact .inl h
      · exact absurd h hint
      · exact .inr ⟨k, h1, h2, h3⟩
    · exact .inr ⟨k, h1, h2, h3⟩
  refine ⟨Stack_cons.2 ⟨?_, ?_⟩, fun i hi hsz => ItemInv.congr (s := s) rfl rfl (fun _ _ => Iff.rfl) (h.nodes i hi hsz)⟩
  · refine EntryInv'.congr (s := s) rfl rfl (fun _ _ => Iff.rfl) ⟨(Graph.ConnEdges.congr fun e _ => hE e).2 ht.conn, ?_⟩
    refine (Graph.AttachedIn.congr fun e _ => hE e).2 (ht.attached.strengthen.mono ?_)
    rintro v ⟨hv, -, hint⟩
    rcases hv with hv | ⟨_, h, _⟩
    · exact .inl (hT v hint hv)
    · simp at h
  · refine (hrest.transfer ?_).congr (s := s) rfl rfl (fun _ _ _ _ => Iff.rfl)
    intro u hu v hv htouch
    obtain ⟨t'', ht'', hT'⟩ := hv
    simp only [List.mem_singleton] at ht''; subst ht''
    refine .inr ⟨t', by simp, hT v (fun hint => ?_) hT'⟩
    obtain ⟨e, he, hue, hinc⟩ := htouch
    exact hdisj u hu e he hue (hint e he hinc)

/-- Popping the top entry `t`: every remaining entry attached at a terminal of `t` has it as its
own terminal. -/
theorem Inv'.pop (t : TEntry) (rest : List TEntry) (hs : s.tstack = t :: rest) (h : s.Inv' D)
    (hgone : ∀ u ∈ rest, ∀ v, t.Term D s v → s.g.Touches (u.edges s.g s.items) v → u.Term D s v) :
    Inv' D { s with tstack := rest } := by
  have hst := h.entries; rw [hs] at hst
  obtain ⟨-, hrest⟩ := Stack_cons.1 hst
  refine ⟨(hrest.transfer ?_).congr (s := s) rfl rfl (fun _ _ _ _ => Iff.rfl),
    fun i hi hsz => ItemInv.congr (s := s) rfl rfl (fun _ _ => Iff.rfl) (h.nodes i hi hsz)⟩
  intro u hu v hv htouch
  obtain ⟨t'', ht'', hT⟩ := hv
  simp only [List.mem_singleton] at ht''; subst ht''
  exact .inl (hgone u hu v hT htouch)

/-- `finishTstackTop` (the hypotheses of `finishTstackTop_complete`). -/
theorem Inv'.finishTop (item : ItemId) (t : TEntry) (rest : List TEntry)
    (hs : s.tstack = t :: rest) (h : s.Inv' D)
    (hitem : item < s.items.size) (hnode : 1 + s.g.nv + s.g.ne ≤ item)
    (hroot : ∀ p, ¬ Items.IsParent s.items p item)
    (hfree : ∀ t' ∈ s.tstack, item ∉ t'.spans.1 ++ t'.spans.2)
    (hside : getSide t.spans (!s.stackDir[t.topDepth]!) = [])
    (hmid : ∀ k, t.topDepth < k → k ≤ D → s.stackVerts[k]! = t.vStart ∨
      s.g.Interior (t.edges s.g s.items) s.stackVerts[k]! ∨ ¬ s.g.Touches (t.edges s.g s.items) s.stackVerts[k]!) :
    ((finishTstackTop item).run s).2.Inv' D := by
  rw [finishTstackTop_run_eq s item t rest hs]
  set dir := s.stackDir[t.topDepth]!
  set f : Item → Item := fun it =>
    { it with vs := setSides dir (some s.stackVerts[t.topDepth]!) (some t.vStart),
              ch := getSide t.spans dir }
  have hst := h.entries; rw [hs] at hst
  obtain ⟨ht', hrest⟩ := Stack_cons.1 hst
  have ht := ht'.toEntryInv
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
  have hE' : ∀ e, e < s.g.ne →
      (TEntry.edges s.g (s.items.modify item f) { t with spans := setSides dir [item] [] } e ↔
        t.edges s.g s.items e) := by
    intro e he
    rw [← hsub e he]
    simp only [TEntry.edges]
    constructor
    · rintro ⟨i, hi, hb⟩
      rwa [List.mem_singleton.1 ((mem_setSides dir [item] i).1 hi)] at hb
    · exact fun hb => ⟨item, (mem_setSides dir [item] item).2 (List.mem_singleton.2 rfl), hb⟩
  refine ⟨Stack_cons.2 ⟨?_, ?_⟩, ?_⟩
  · exact ⟨(Graph.ConnEdges.congr hE').2 ht.conn,
      (Graph.AttachedIn.congr hE').2 (ht.attached.mono fun _ hv => .inl hv)⟩
  · refine (hrest.transfer ?_).congr (s := s) rfl rfl
      (fun u hu e _ => TEntry.edges_modify_of_not_mem item f hroot (hfree u (by simp [hs, hu])) e)
    intro u _ v hv _
    obtain ⟨t'', ht'', hT⟩ := hv
    rw [List.mem_singleton] at ht''; rw [ht''] at hT
    exact .inr ⟨{ t with spans := setSides dir [item] [] }, by simp, hT⟩
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

end WalkState

end Spqr
