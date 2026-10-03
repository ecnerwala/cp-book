import Spqr.RangesWalk

/-! # The close records at the six close sites of `finishEdge`

`CloseInv` preservation through one `finishEdge`, split into one named statement per block that
creates or completes an item record (PROOF.md §4.6). Each is stated at the state where its block
runs, under `CloseCtx`: what the walk induction knows at a `finishEdge` call (the range
invariant, the ear book-keeping, the frontier, and the scheduled adjacency `FinishR`). The
blocks between them (`feS₀`, `mergeLate`, the `closeVert'` merges/unwrap/retarget, `finishTail`)
only move spans and are covered by the proved `CloseInv` frame lemmas. The empirical checker
(`checks/RangesInvCheck.lean`, `closeSites`) evaluates every `CloseAt` clause at every one of
these block boundaries, under the statement's name. -/

namespace Spqr.Items

theorem vs_modify_ch {items : Items} (j : ItemId) (cs : List ItemId) (p : ItemId) :
    Items.vs (items.modify j fun it => { it with ch := cs }) p = items.vs p :=
  vs_modify_of_vs items j p (fun it => { it with ch := cs }) fun _ => rfl

theorem type_modify_ch {items : Items} (j : ItemId) (cs : List ItemId) (p : ItemId) :
    Items.type (items.modify j fun it => { it with ch := cs }) p = items.type p :=
  type_modify items j p (fun it => { it with ch := cs }) fun _ => rfl

theorem type_modify_vs {items : Items} (j : ItemId) (vsv : Option Nat × Option Nat) (p : ItemId) :
    Items.type (items.modify j fun it => { it with vs := vsv }) p = items.type p :=
  type_modify items j p (fun it => { it with vs := vsv }) fun _ => rfl

theorem ch_modify_vs {items : Items} (j : ItemId) (vsv : Option Nat × Option Nat) (p : ItemId) :
    Items.ch (items.modify j fun it => { it with vs := vsv }) p = items.ch p :=
  ch_modify_ch_eq j (fun it => { it with vs := vsv }) (fun _ => rfl) p

theorem Below_modify_vs {items : Items} (j : ItemId) (vsv : Option Nat × Option Nat) (a i : ItemId) :
    Items.Below (items.modify j fun it => { it with vs := vsv }) a i ↔ items.Below a i :=
  Below_modify_ch_eq j (fun it => { it with vs := vsv }) fun _ => rfl

/-- The record of a vertex item: no I/O children, every child a nonempty Q. -/
theorem CloseAt.vertex' {g : Graph} {items : Items} {v : Nat} (hv : v < g.nv)
    (ht : items.type (vertItem v) = .V)
    (hio : ∀ c, items.IsParent (vertItem v) c → items.type c ≠ .I ∧ items.type c ≠ .O)
    (hne : ∀ c, items.IsParent (vertItem v) c → items.ch c ≠ []) :
    CloseAt g items (vertItem v) := by
  constructor
  · simp [ht]
  · simp [ht]
  · simp [ht]
  · simp [ht]
  · simp [ht]
  · intro c h hio'; exact (hio'.elim (hio c h).1 (hio c h).2).elim
  · intro e he; exact absurd he (by show 1 + v ≠ 1 + g.nv + e; omega)
  · intro e he; exact absurd he (by show 1 + v ≠ 1 + g.nv + e; omega)
  · intro _ _ _ c h; exact hne c h
  · simp [ht]
  · simp [ht]
  · simp [ht]

/-- The record of a nonempty (root) Q item attached only at `u`. -/
theorem CloseAt.rootQ {g : Graph} {items : Items} {e u c : Nat} (he : e < g.ne)
    (ht : items.type (edgeItem g e) = .Q) (hne : items.ch (edgeItem g e) ≠ [])
    (hvs : items.vs (edgeItem g e) = (some u, none)) (hu : g.Inc e u)
    (hatt : ∀ v, items.Att g (edgeItem g e) v → v = u)
    (hct : items.type c ∉ [NodeType.F, .V]) (hcq : items.type c = .Q → items.ch c = [])
    (hloop : (g.edges[e]!).1 = (g.edges[e]!).2 →
      items.ch (edgeItem g e) = [c] ∧ items.vs c = (some u, none))
    (hnon : (g.edges[e]!).1 ≠ (g.edges[e]!).2 → ∃ w, w < g.nv ∧ PairEq (u, w) g.edges[e]! ∧
      items.ch (edgeItem g e) = [c, vertItem w] ∧
      ∃ a b, items.vs c = (some a, some b) ∧ PairEq (a, b) (u, w)) :
    CloseAt g items (edgeItem g e) := by
  constructor
  · intro _ v hv; left; rw [hvs, hatt v hv]
  · rintro (h | ⟨_, h⟩)
    · simp [ht] at h
    · exact absurd h hne
  · simp [ht]
  · simp [ht]
  · intro _ _ h; simp [ht] at h
  · intro _ _ _; exact ht
  · intro _ _ _ h; exact absurd h hne
  · intro e' he' _ _
    obtain rfl : e = e' := by have : 1 + g.nv + e = 1 + g.nv + e' := he'; omega
    exact ⟨u, c, hvs, hu, hct, hcq, hloop, hnon⟩
  · intro v hv hvn; exact absurd hv (by show 1 + g.nv + e ≠ 1 + v; omega)
  · simp [ht]
  · simp [ht]
  · simp [ht]

theorem PairEq.trans' {p q r : Nat × Nat} (h₁ : PairEq p q) (h₂ : PairEq q r) : PairEq p r := by
  obtain ⟨a, b⟩ := p; obtain ⟨c, d⟩ := q; obtain ⟨e, f⟩ := r
  simp [PairEq] at *
  omega

theorem PairEq.ne_of_ne {a b u v : Nat} (h : PairEq (a, b) (u, v)) (huv : u ≠ v) : a ≠ b := by
  simp [PairEq] at h
  omega

theorem mem_SPRQ_ne {t : NodeType} (h : t ∈ [NodeType.S, .P, .R, .Q]) :
    t ≠ .F ∧ t ≠ .V ∧ t ≠ .I ∧ t ≠ .O := by
  simp only [List.mem_cons, List.not_mem_nil, or_false] at h
  rcases h with rfl | rfl | rfl | rfl <;> decide

/-- The close record of a P node `x` with children `cs`, each a closed S/P/R or leaf Q on the
terminals `{u, v}` and attached nowhere else; both terminals touch `x` and keep an edge outside it. -/
theorem CloseAt.pNode {g : Graph} {items : Items} {x : ItemId} {u v : Nat} {cs : List ItemId}
    (hx : 1 + g.nv + g.ne ≤ x) (ht : items.type x = .P) (hch : items.ch x = cs)
    (hvs : items.vs x = (some u, some v)) (huv : u ≠ v) (hlen : 2 ≤ cs.length)
    (hV : ∀ w, w < g.nv → items.type (vertItem w) = .V)
    (hty : ∀ c ∈ cs, items.type c ∈ [NodeType.S, .P, .R, .Q])
    (hcvs : ∀ c ∈ cs, ∃ a b, items.vs c = (some a, some b) ∧ PairEq (a, b) (u, v))
    (hatt : ∀ c ∈ cs, ∀ w, items.Att g c w → w = u ∨ w = v)
    (htouch : ∀ w, w = u ∨ w = v → ∃ e, e < g.ne ∧ g.Inc e w ∧ items.EdgeBelow g x e)
    (hpend : ∀ w, w = u ∨ w = v → ∃ e, e < g.ne ∧ g.Inc e w ∧ ¬ items.EdgeBelow g x e) :
    CloseAt g items x := by
  have hpar : ∀ c, items.IsParent x c ↔ c ∈ cs := fun c => by simp [IsParent, hch]
  have hsplit : ∀ e, e < g.ne → items.EdgeBelow g x e → ∃ c ∈ cs, items.EdgeBelow g c e := by
    intro e he hb
    rcases hb.head_cases with heq | ⟨c, hc, hb⟩
    · have : x = 1 + g.nv + e := heq
      rw [this] at hx
      omega
    · exact ⟨c, (hpar c).1 hc, hb⟩
  have hjoin : ∀ c ∈ cs, ∀ e, items.EdgeBelow g c e → items.EdgeBelow g x e :=
    fun c hc e hb => Relation.ReflTransGen.head ((hpar c).2 hc) hb
  have hIsVs : ∀ w, items.IsVs x w ↔ w = u ∨ w = v := by
    intro w
    simp only [IsVs, hvs]
    constructor
    · rintro (h | h)
      · exact .inl (Option.some.inj h).symm
      · exact .inr (Option.some.inj h).symm
    · rintro (rfl | rfl) <;> simp
  have hattx : ∀ w, items.Att g x w → w = u ∨ w = v := by
    rintro w ⟨e, e', he, he', hi, hi', hb, hnb⟩
    obtain ⟨c, hc, hbc⟩ := hsplit e he hb
    exact hatt c hc w ⟨e, e', he, he', hi, hi', hbc, fun h => hnb (hjoin c hc e' h)⟩
  have hne : ∀ c ∈ cs, items.type c ≠ .V := fun c hc => (mem_SPRQ_ne (hty c hc)).2.1
  constructor
  · intro _ w hw; exact (hIsVs w).2 (hattx w hw)
  · intro _ w hw
    obtain ⟨e, he, hi, hb⟩ := htouch w ((hIsVs w).1 hw)
    obtain ⟨e', he', hi', hnb⟩ := hpend w ((hIsVs w).1 hw)
    exact ⟨e, e', he, he', hi, hi', hb, hnb⟩
  · intro _ a b hab
    rw [hvs] at hab
    simp only [Prod.mk.injEq, Option.some.injEq] at hab
    omega
  · intro _ w hw
    constructor
    · intro hp
      exact absurd (hV w hw) (hne _ ((hpar _).1 hp))
    · rintro ⟨⟨⟨e₀, he₀, hi₀⟩, hall⟩, hno⟩
      obtain ⟨c, hc, hbc⟩ := hsplit e₀ he₀ (hall e₀ he₀ hi₀)
      have hnot := hno c ((hpar c).2 hc)
      simp only [not_forall, Classical.not_imp] at hnot
      obtain ⟨e₁, he₁, hi₁, hnb₁⟩ := hnot
      obtain ⟨e₂, he₂, hi₂, hnb₂⟩ :=
        hpend w (hatt c hc w ⟨e₀, e₁, he₀, he₁, hi₀, hi₁, hbc, hnb₁⟩)
      exact (hnb₂ (hall e₂ he₂ hi₂)).elim
  · intro c hc _ _
    obtain ⟨a, b, hab, hp⟩ := hcvs c ((hpar c).1 hc)
    exact ⟨a, b, hab, hp.ne_of_ne huv⟩
  · intro c hc hio
    have := mem_SPRQ_ne (hty c ((hpar c).1 hc))
    rcases hio with hio | hio
    · exact absurd hio this.2.2.1
    · exact absurd hio this.2.2.2
  · intro e he hel _
    have : x = 1 + g.nv + e := he
    rw [this] at hx
    omega
  · intro e he hel _
    have : x = 1 + g.nv + e := he
    rw [this] at hx
    omega
  · intro w hw hwl
    have : x = 1 + w := hw
    rw [this] at hx
    omega
  · intro _
    have hfilt : (items.ch x).filter (fun c => items.type c ≠ .V) = cs := by
      rw [hch]
      exact List.filter_eq_self.mpr fun c hc => decide_eq_true (hne c hc)
    refine ⟨?_, fun c hc => hne c ((hpar c).1 hc), ?_⟩
    · unfold virtualEdges
      rw [hfilt, List.length_map]
      exact hlen
    · intro q hq
      unfold virtualEdges at hq
      rw [hfilt, List.mem_map] at hq
      obtain ⟨c, hc, rfl⟩ := hq
      obtain ⟨a, b, hab, hp⟩ := hcvs c hc
      exact ⟨u, v, hvs, by rw [hab]; exact hp⟩
  · intro h; rw [ht] at h; cases h
  · intro h; rw [ht] at h; cases h

/-- `CloseAt.pNode` with the terminals written by `makeVs`. -/
theorem CloseAt.pNode' {g : Graph} {items : Items} {x : ItemId} {u v : Nat} {cs : List ItemId}
    (dir : Bool) (hx : 1 + g.nv + g.ne ≤ x) (ht : items.type x = .P) (hch : items.ch x = cs)
    (hvs : items.vs x = setSides dir (some u) (some v)) (huv : u ≠ v) (hlen : 2 ≤ cs.length)
    (hV : ∀ w, w < g.nv → items.type (vertItem w) = .V)
    (hty : ∀ c ∈ cs, items.type c ∈ [NodeType.S, .P, .R, .Q])
    (hcvs : ∀ c ∈ cs, ∃ a b, items.vs c = (some a, some b) ∧ PairEq (a, b) (u, v))
    (hatt : ∀ c ∈ cs, ∀ w, items.Att g c w → w = u ∨ w = v)
    (htouch : ∀ w, w = u ∨ w = v → ∃ e, e < g.ne ∧ g.Inc e w ∧ items.EdgeBelow g x e)
    (hpend : ∀ w, w = u ∨ w = v → ∃ e, e < g.ne ∧ g.Inc e w ∧ ¬ items.EdgeBelow g x e) :
    CloseAt g items x := by
  cases dir
  · exact pNode hx ht hch hvs huv hlen hV hty hcvs hatt htouch hpend
  · exact pNode hx ht hch hvs huv.symm hlen hV hty
      (fun c hc => let ⟨a, b, hab, hp⟩ := hcvs c hc; ⟨a, b, hab, Or.symm hp⟩)
      (fun c hc w hw => (hatt c hc w hw).symm) (fun w hw => htouch w hw.symm)
      (fun w hw => hpend w hw.symm)

/-- The record of a freshly closed S/R node `x` with children `cs` and terminals `(u, v)`. -/
theorem CloseAt.node {g : Graph} {items : Items} {x : ItemId} {u v : Nat} {cs : List ItemId}
    (hx : 1 + g.nv + g.ne ≤ x) (ht : items.type x = .S ∨ items.type x = .R) (hch : items.ch x = cs)
    (hvs : items.vs x = (some u, some v)) (huv : u ≠ v)
    (hty : ∀ c ∈ cs, items.type c ∈ [NodeType.S, .P, .R, .Q, .V])
    (htwo : ∀ c ∈ cs, items.type c ≠ .V → ∃ a b, items.vs c = (some a, some b) ∧ a ≠ b)
    (hatt : ∀ w, items.Att g x w → w = u ∨ w = v)
    (htouch : ∀ w, w = u ∨ w = v → ∃ e, e < g.ne ∧ g.Inc e w ∧ items.EdgeBelow g x e)
    (hpend : ∀ w, w = u ∨ w = v → ∃ e, e < g.ne ∧ g.Inc e w ∧ ¬ items.EdgeBelow g x e)
    (hint : ∀ w, w < g.nv → (vertItem w ∈ cs ↔ items.Inner g x w ∧
      ∀ c ∈ cs, ¬ ∀ e, e < g.ne → g.Inc e w → items.EdgeBelow g c e))
    (hs : items.type x = .S → ∃ xs, (cs.filter fun c => items.type c = .V) = xs.map vertItem ∧
      1 ≤ xs.length ∧ items.virtualEdges x = List.zip (u :: xs) (xs ++ [v]))
    (hr : items.type x = .R → 2 ≤ (cs.filter fun c => items.type c = .V).length ∧
      5 ≤ (items.virtualEdges x).length ∧
      ((items.virtualEdges x).map fun q => min q.1 q.2 + g.nv * max q.1 q.2).Nodup ∧
      ∀ q ∈ items.virtualEdges x, ¬ PairEq q (u, v)) :
    CloseAt g items x := by
  have hpar : ∀ c, items.IsParent x c ↔ c ∈ cs := fun c => by simp [IsParent, hch]
  have hIsVs : ∀ w, items.IsVs x w ↔ w = u ∨ w = v := by
    intro w
    simp only [IsVs, hvs]
    constructor
    · rintro (h | h)
      · exact .inl (Option.some.inj h).symm
      · exact .inr (Option.some.inj h).symm
    · rintro (rfl | rfl) <;> simp
  constructor
  · intro _ w hw; exact (hIsVs w).2 (hatt w hw)
  · intro _ w hw
    obtain ⟨e, he, hi, hb⟩ := htouch w ((hIsVs w).1 hw)
    obtain ⟨e', he', hi', hnb⟩ := hpend w ((hIsVs w).1 hw)
    exact ⟨e, e', he, he', hi, hi', hb, hnb⟩
  · intro _ a b hab
    rw [hvs] at hab
    simp only [Prod.mk.injEq, Option.some.injEq] at hab
    omega
  · intro _ w hw
    rw [hpar, hint w hw]
    constructor
    · rintro ⟨h1, h2⟩; exact ⟨h1, fun c hc => h2 c ((hpar c).1 hc)⟩
    · rintro ⟨h1, h2⟩; exact ⟨h1, fun c hc => h2 c ((hpar c).2 hc)⟩
  · intro c hc _ hcV; exact htwo c ((hpar c).1 hc) hcV
  · intro c hc hio
    have := hty c ((hpar c).1 hc)
    rcases hio with h | h <;> (rw [h] at this; simp at this)
  · intro e he hel _
    have : x = 1 + g.nv + e := he
    rw [this] at hx
    omega
  · intro e he hel _
    have : x = 1 + g.nv + e := he
    rw [this] at hx
    omega
  · intro w hw hwl
    have : x = 1 + w := hw
    rw [this] at hx
    omega
  · intro h; rcases ht with h' | h' <;> rw [h'] at h <;> cases h
  · intro h
    obtain ⟨xs, hf, hl, hve⟩ := hs h
    exact ⟨u, v, xs, hvs, by rw [hch]; exact hf, hl, hve⟩
  · intro h
    obtain ⟨h1, h2, h3, h4⟩ := hr h
    refine ⟨by rw [hch]; exact h1, h2, h3, fun a b hab => ?_⟩
    rw [hvs] at hab
    simp only [Prod.mk.injEq, Option.some.injEq] at hab
    obtain ⟨rfl, rfl⟩ := hab
    exact h4

end Spqr.Items

namespace Spqr.WalkState
open WalkM

/-- The state `finishRest` (hence `finishP`) runs from: after `closeVert'` for a tree edge with
a vertex ear, after `mergeLate` for a tree edge without one, after the push for a back edge. -/
def feRest (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool) (s : WalkState) :
    WalkState :=
  if o.cls.isTree then
    if hasVert then feS₃ curV d o origTstack s else feS₂ d o s
  else feBack curV (o.cls.lowval d) d o s

/-- The top two entries at a type-1 P close (`condP` holds at `r`, so `nxt` is the enclosing
`(curV, lv)` entry; `cur` holds the closed child): each is a single item on the side
`stackDir[lv]` whose terminals are `{stackVerts[lv], curV}`, touches both and is attached nowhere
else; the items, and the children of a P item among them, are S/P/R or leaf Q; each terminal keeps
an edge outside both entries; each item of the two entries is a root listed once on the stack. -/
structure PSite (curV lv : Nat) (r : WalkState) : Prop where
  shape : Shape r
  stack : ∃ cur rest, r.tstack = cur :: rest ∧ cur.vStart = curV ∧ cur.topDepth = lv
  ne : curV ≠ r.stackVerts[lv]!
  single : ∀ t ∈ r.tstack.take 2, ∃ j, t.spans = setSides r.stackDir[lv]! [j] [] ∧
    ∃ a b, Items.vs r.items j = (some a, some b) ∧ Items.PairEq (a, b) (r.stackVerts[lv]!, curV)
  att : ∀ t ∈ r.tstack.take 2, ∀ w, r.g.Touches (t.edges r.g r.items) w →
    w = curV ∨ w = r.stackVerts[lv]! ∨ r.g.Interior (t.edges r.g r.items) w
  touch : ∀ t ∈ r.tstack.take 2, r.g.Touches (t.edges r.g r.items) curV ∧
    r.g.Touches (t.edges r.g r.items) r.stackVerts[lv]!
  pend : ∀ w, w = curV ∨ w = r.stackVerts[lv]! → ∃ e, e < r.g.ne ∧ r.g.Inc e w ∧
    ∀ t ∈ r.tstack.take 2, ¬ t.edges r.g r.items e
  kinds : ∀ t ∈ r.tstack.take 2, ∀ j ∈ t.spans.1 ++ t.spans.2, ∀ j',
    j' = j ∨ (Items.type r.items j = .P ∧ Items.IsParent r.items j j') →
    Items.type r.items j' ∈ [NodeType.S, .P, .R, .Q] ∧
      (Items.type r.items j' = .Q → Items.ch r.items j' = [])
  once : ∀ t ∈ r.tstack.take 2, ∀ j ∈ t.spans.1 ++ t.spans.2,
    (∀ p, ¬ Items.IsParent r.items p j) ∧ spansCount r.tstack j = 1

/-- The children side, terminals and virtual edges of the node closed from the top entry `t`. -/
def vKids (t : TEntry) (r : WalkState) : List ItemId := getSide t.spans r.stackDir[t.topDepth]!
def vTerms (curV : Nat) (t : TEntry) (r : WalkState) : Nat × Nat :=
  setSides r.stackDir[t.topDepth]! r.stackVerts[t.topDepth]! curV
def vVirt (t : TEntry) (r : WalkState) : List (Nat × Nat) :=
  ((vKids t r).filter fun c => Items.type r.items c ≠ .V).map fun c =>
    ((Items.vs r.items c).1.getD 0, (Items.vs r.items c).2.getD 0)

/-- The vertex close of a type-1 tree edge (`closeVertTail_closeAt`), at the state `cvS₅` the
`finishTstackTop x` of the `some item` arm runs from: `t` is the top entry, re-targeted to `curV`,
with everything on its `stackDir[topDepth]` side; `x` is the free S/R node about to receive that
side as children; the terminals `stackVerts[topDepth] ≠ curV` are touched and are the only
attachments of the side, each with an edge outside it; children are S/P/R/Q/V, non-V ones with two
distinct terminals; a vertex item is a child iff all its edges lie in the side and in no single
child; and the S path / R shape clauses hold for the side's virtual edges. -/
structure VSite (curV : Nat) (x : ItemId) (t : TEntry) (r : WalkState) : Prop where
  shape : Shape r
  stack : ∃ rest, r.tstack = t :: rest
  vstart : t.vStart = curV
  side : getSide t.spans (!r.stackDir[t.topDepth]!) = []
  free : ItemFree r x
  ty : Items.type r.items x = .S ∨ Items.type r.items x = .R
  ne : curV ≠ r.stackVerts[t.topDepth]!
  kinds : ∀ c ∈ vKids t r, Items.type r.items c ∈ [NodeType.S, .P, .R, .Q, .V]
  two : ∀ c ∈ vKids t r, Items.type r.items c ≠ .V →
    ∃ a b, Items.vs r.items c = (some a, some b) ∧ a ≠ b
  att : ∀ w, r.g.Touches (t.edges r.g r.items) w →
    w = curV ∨ w = r.stackVerts[t.topDepth]! ∨ r.g.Interior (t.edges r.g r.items) w
  touch : r.g.Touches (t.edges r.g r.items) curV ∧
    r.g.Touches (t.edges r.g r.items) r.stackVerts[t.topDepth]!
  pend : ∀ w, w = curV ∨ w = r.stackVerts[t.topDepth]! → ∃ e, e < r.g.ne ∧ r.g.Inc e w ∧
    ¬ t.edges r.g r.items e
  inner : ∀ w, w < r.g.nv → (vertItem w ∈ vKids t r ↔
    r.g.Touches (t.edges r.g r.items) w ∧ r.g.Interior (t.edges r.g r.items) w ∧
    ∀ c ∈ vKids t r, ¬ r.g.Interior (Items.EdgeBelow r.g r.items c) w)
  s_order : Items.type r.items x = .S → ∃ xs,
    ((vKids t r).filter fun c => Items.type r.items c = .V) = xs.map vertItem ∧ 1 ≤ xs.length ∧
    vVirt t r = List.zip ((vTerms curV t r).1 :: xs) (xs ++ [(vTerms curV t r).2])
  r_shape : Items.type r.items x = .R →
    2 ≤ ((vKids t r).filter fun c => Items.type r.items c = .V).length ∧ 5 ≤ (vVirt t r).length ∧
    ((vVirt t r).map fun q => min q.1 q.2 + r.g.nv * max q.1 q.2).Nodup ∧
    ∀ q ∈ vVirt t r, ¬ Items.PairEq q (vTerms curV t r)

/-- The context of a `finishEdge curV d o origTstack hasVert` call of the walk at state `s`. -/
structure CloseCtx (σ : List Nat) (n D curV d : Nat) (o : DfsOut) (origTstack : Nat)
    (hasVert : Bool) (s : WalkState) : Prop where
  hD : D = if o.cls.isTree then d + 1 else d
  nodup : σ.Nodup
  lt : ∀ e ∈ σ, e < s.g.ne
  pos : σ[n]? = some o.e
  block : o.block <:+: σ
  ranges : s.RangesInv σ n D
  shape : Shape s
  guards : FinishGuards d o origTstack hasVert s
  book : FinishBook curV d o origTstack hasVert s
  frontier : Frontier (o := o) d origTstack s
  finishR : FinishR σ n curV d o origTstack hasVert s
  close : s.CloseInv
  /-- DFS facts: the out-edge's endpoints; the open path's tree edges, still pending (after `o.e`
  in `σ`); a returning child has an edge of its own. -/
  ends : o.Ends s.g curV
  path : ∀ k, k < d → ∃ e, e < s.g.ne ∧ n < σ.idxOf e ∧
    Items.PairEq (s.stackVerts[k]!, s.stackVerts[k + 1]!) s.g.edges[e]!
  dest_edge : o.cls.isTree = true → o.cls.lowval d < d →
    ∃ e, e < s.g.ne ∧ e ≠ o.e ∧ s.g.Inc e o.dest
  dest_lt : o.dest < s.g.nv
  /-- Boundary (`d ≤ lowval`) facts the ear layer does not export: a self-loop returns to `curV`;
  the popped block entry holds exactly the child's vertex item; at a component edge the
  `(o.dest, lowval)` entry above it holds one closed node, or leaf Q, on the block's terminals
  `curV, o.dest`. -/
  bd_loop : o.cls.isTree = false → d ≤ o.cls.lowval d → o.dest = curV
  bd_vert : o.cls.isTree = true → d ≤ o.cls.lowval d →
    ∀ t ∈ (if o.cls.lowval d = d + 1 then s.tstack.head? else s.tstack.tail.head?),
      t.spans.2 = [vertItem o.dest]
  bd_node : o.cls.isTree = true → d ≤ o.cls.lowval d → o.cls.lowval d ≠ d + 1 →
    ∀ b ∈ s.tstack.head?, ∃ c, b.spans.1 = [c] ∧ Items.type s.items c ∉ [NodeType.F, .V] ∧
      (Items.type s.items c = .Q → Items.ch s.items c = []) ∧
      ∃ a b', Items.vs s.items c = (some a, some b') ∧ Items.PairEq (a, b') (curV, o.dest)
  /-- The P site (`finishP_closeAt`): `PSite` at the state `finishP` runs from, whenever `condP`
  holds there. -/
  p_site : o.cls.lowval d < d → o.cls.isType1 = true →
    result (condP curV (o.cls.lowval d) true) (feRest curV d o origTstack hasVert s) = true →
    PSite curV (o.cls.lowval d) (feRest curV d o origTstack hasVert s)
  /-- The vertex-close site (`closeVertTail_closeAt`): `VSite` for the unwrapped item at `cvS₅`. -/
  v_site : o.cls.isTree = true → o.cls.lowval d < d → hasVert = true → o.cls.isType1 = true →
    ∃ t, VSite curV ((maybeUnwrapNxt (if feSingle d o s then NodeType.S else .R)).run (feS₂ d o s)).1 t
      (cvS₅ curV s.stackDir[d]! true origTstack (feSingle d o s) (feS₂ d o s))

variable {σ : List Nat} {n D curV d : Nat} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}
  {s : WalkState}

theorem CloseCtx.idxOf (hc : CloseCtx σ n D curV d o origTstack hasVert s) : σ.idxOf o.e = n := by
  obtain ⟨hn, he⟩ := List.getElem?_eq_some_iff.mp hc.pos
  rw [← he]; exact hc.nodup.idxOf_getElem n hn

theorem CloseCtx.ne_of_path (hc : CloseCtx σ n D curV d o origTstack hasVert s) {e : Nat}
    (h : n < σ.idxOf e) : e ≠ o.e := fun heq => by
  rw [heq, hc.idxOf] at h; exact Nat.lt_irrefl _ h

theorem cnt_eq_zero_of_free {i : ItemId} (hs : ∀ t ∈ s.tstack, i ∉ t.spans.1 ++ t.spans.2)
    (hp : ∀ p, ¬ Items.IsParent s.items p i) : s.cnt i = 0 := by
  have h1 : spansCount s.tstack i = 0 := by
    apply List.sum_eq_zero
    intro x hx
    obtain ⟨t, ht, rfl⟩ := List.mem_map.mp hx
    exact List.count_eq_zero.mpr (hs t ht)
  have h2 : chCount s.items i = 0 := Finset.sum_eq_zero fun j _ => List.count_eq_zero.mpr (hp j)
  dsimp [cnt]; omega

theorem pairEq_swap' {a b : Nat} {q : Nat × Nat} (h : Items.PairEq (a, b) q) :
    Items.PairEq (b, a) q := by
  rcases h with h | h
  · exact Or.inr (by rw [← h])
  · obtain ⟨ha, hb⟩ := Prod.mk.inj h; exact Or.inl (by rw [ha, hb])

/-- `feS₀` writes the terminals of the pending edge into its root, loose Q item. -/
theorem feS₀_closeInv (hc : CloseCtx σ n D curV d o origTstack hasVert s) : (feS₀ d o s).CloseInv := by
  obtain ⟨sub, base, -, hE⟩ := hc.book.ear
  exact hc.close.modifyLoose (x := edgeItem s.g o.e) (by show s.g.nv < 1 + s.g.nv + o.e; omega)
    (cnt_eq_zero_of_free hE.q_free hE.q_root) _ (fun _ => rfl)

/-- Pushing the pending edge after `feS₀`: the leaf-Q record (`CloseInv.pushEdge`) from the
written terminals `(stackVerts[d], o.dest)`, their distinctness, and another edge at each. -/
theorem leafQ_closeInv (hc : CloseCtx σ n D curV d o origTstack hasVert s) {v dd : Nat}
    (hne : s.stackVerts[d]! ≠ o.dest) (hp : Items.PairEq (s.stackVerts[d]!, o.dest) s.g.edges[o.e]!)
    (h₁ : ∃ e, e < s.g.ne ∧ e ≠ o.e ∧ s.g.Inc e s.stackVerts[d]!)
    (h₂ : ∃ e, e < s.g.ne ∧ e ≠ o.e ∧ s.g.Inc e o.dest) :
    (after (pushEdgeTstack v dd o.e) (feS₀ d o s)).CloseInv := by
  have hlt : edgeItem s.g o.e < s.items.size := by
    have := hc.shape.size; have := hc.book.e_lt; show 1 + s.g.nv + o.e < _; omega
  have ht : Items.type (feS₀ d o s).items (edgeItem (feS₀ d o s).g o.e) = .Q := by
    show Items.type (s.items.modify _ _) (edgeItem s.g o.e) = .Q
    rw [Items.type_modify (f := fun it =>
      { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) _ _ _
      (fun _ => rfl)]
    exact hc.shape.edge _ hc.book.e_lt
  have hch : Items.ch (feS₀ d o s).items (edgeItem (feS₀ d o s).g o.e) = [] := by
    show Items.ch (s.items.modify _ _) (edgeItem s.g o.e) = []
    rw [Items.ch_modify_of_ch (f := fun it =>
      { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) _ _ _
      (fun _ => rfl)]
    exact hc.book.q
  have hvs : Items.vs (feS₀ d o s).items (edgeItem (feS₀ d o s).g o.e) =
      setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) := by
    show Items.vs (s.items.modify _ _) (edgeItem s.g o.e) = _
    rw [Items.vs_modify_self _ _ _ hlt]
  have hatt : ∀ w, s.g.Inc o.e w → (∃ e, e < s.g.ne ∧ e ≠ o.e ∧ s.g.Inc e w) →
      Items.Att (feS₀ d o s).g (feS₀ d o s).items (edgeItem (feS₀ d o s).g o.e) w := by
    intro w hw ⟨e', he', hne', hinc⟩
    refine ⟨o.e, e', hc.book.e_lt, he', hw, hinc, Relation.ReflTransGen.refl, fun hb => hne' ?_⟩
    have h := Items.below_eq_of_ch_nil hch hb
    have : 1 + s.g.nv + e' = 1 + s.g.nv + o.e := h
    omega
  have h₀ := feS₀_closeInv hc
  obtain ⟨hi₁, hi₂⟩ := Graph.inc_of_pairEq hp
  cases hdir : s.stackDir[d]! with
  | false =>
    exact h₀.pushEdge v dd o.e _ _ ht hch (by rw [hvs, hdir]; rfl) hne hp (hatt _ hi₁ h₁) (hatt _ hi₂ h₂)
  | true =>
    exact h₀.pushEdge v dd o.e _ _ ht hch (by rw [hvs, hdir]; rfl) (Ne.symm hne) (pairEq_swap' hp)
      (hatt _ hi₂ h₂) (hatt _ hi₁ h₁)

/-- The parent edge of `curV = stackVerts[d]` (`d ≥ 1`): pending, hence not `o.e`. -/
theorem CloseCtx.parent_edge (hc : CloseCtx σ n D curV d o origTstack hasVert s) (hd : 0 < d) :
    ∃ e, e < s.g.ne ∧ e ≠ o.e ∧ s.g.Inc e s.stackVerts[d]! := by
  obtain ⟨e, he, hn, hp⟩ := hc.path (d - 1) (by omega)
  rw [show d - 1 + 1 = d by omega] at hp
  exact ⟨e, he, hc.ne_of_path hn, (Graph.inc_of_pairEq hp).2⟩

/-- Leaf Q at a returning tree edge: `pushEdgeTstack o.dest d o.e` after `feS₀` keeps every
record (`CloseInv.pushEdge`: distinct endpoints and both terminal `Att` witnesses of
`edgeItem o.e`, whose other incident edges are not below it). -/
theorem closeEars_closeAt (hc : CloseCtx σ n D curV d o origTstack hasVert s)
    (ht : o.cls.isTree = true) (hlow : o.cls.lowval d < d) :
    (ceS₁ o.dest d o.e (feS₀ d o s)).CloseInv := by
  obtain ⟨lv, kind, ho, hlv⟩ := ret_of_lowval_lt hlow
  obtain ⟨sub, base, -, hE⟩ := hc.book.ear
  have hends := hc.book.ends lv kind ho
  simp only [ht, ↓reduceIte] at hends
  exact leafQ_closeInv hc (hE.path_child ht d (Nat.le_refl d)) (pairEq_swap' hends)
    (hc.parent_edge (by omega)) (hc.dest_edge ht hlow)

/-- One iteration of loop 1: the S/P/R record created (or reopened) by `maybeUnwrapNxt` and
closed by `finishTstackTop` over the merged entry — attachment/terminal facts, the interior
V-child equivalence, two terminals of every non-V child, and the S/P/R shape clause. -/
theorem loop1Body_closeAt (hc : CloseCtx σ n D curV d o origTstack hasVert s)
    (ht : o.cls.isTree = true) (hlow : o.cls.lowval d < d) (k : Nat)
    (hk : ∀ j, j ≤ k → result (loop1Cond d)
      (iter (loop1Body d s.stackDir[d]!) j (ceS₁ o.dest d o.e (feS₀ d o s))) = true)
    (h : (iter (loop1Body d s.stackDir[d]!) k (ceS₁ o.dest d o.e (feS₀ d o s))).CloseInv) :
    (after (loop1Body d s.stackDir[d]!)
      (iter (loop1Body d s.stackDir[d]!) k (ceS₁ o.dest d o.e (feS₀ d o s)))).CloseInv := by
  sorry

/-- `finishP`: when `condP` holds, the new or reused P record over the merged `(curV,
stackVerts[lowval])` entry — at least two virtual edges, all equal to its terminals, no V child. -/
theorem getSide_setSides' {α : Type} (dir : Bool) (a b : α) : getSide (setSides dir a b) dir = a := by
  cases dir <;> rfl

theorem setSides_some (dir : Bool) (a b : Nat) :
    setSides dir (some a) (some b) = (some (setSides dir a b).1, some (setSides dir a b).2) := by
  cases dir <;> rfl

theorem setSides_ne {a b : Nat} (dir : Bool) (h : a ≠ b) :
    (setSides dir a b).1 ≠ (setSides dir a b).2 := by
  cases dir <;> simp [setSides] <;> omega

theorem eq_setSides_or (dir : Bool) (a b w : Nat) :
    (w = (setSides dir a b).1 ∨ w = (setSides dir a b).2) ↔ (w = a ∨ w = b) := by
  cases dir <;> simp [setSides, or_comm]

theorem mem_of_mem_getSide {p : List ItemId × List ItemId} {dir : Bool} {c : ItemId}
    (h : c ∈ getSide p dir) : c ∈ p.1 ++ p.2 := by
  cases dir <;> simp_all [getSide]

theorem mem_sides_iff {p : List ItemId × List ItemId} {dir : Bool} {c : ItemId}
    (hside : getSide p (!dir) = []) : c ∈ p.1 ++ p.2 ↔ c ∈ getSide p dir := by
  cases dir <;> simp_all [getSide]

theorem getSide_mergeInto {dir : Bool} {cur nxt : TEntry} {a b : List ItemId}
    (hc : cur.spans = setSides dir a []) (hn : nxt.spans = setSides dir b []) :
    getSide (TEntry.mergeInto cur nxt).spans dir = if dir then b ++ a else a ++ b := by
  cases dir <;> simp [TEntry.mergeInto, getSide, setSides, hc, hn]

theorem CloseInv.mergeTop' (h : s.CloseInv) {cur nxt : TEntry} {rest : List TEntry}
    (hts : s.tstack = cur :: nxt :: rest) :
    ({ s with tstack := TEntry.mergeInto cur nxt :: rest } : WalkState).CloseInv := by
  have := h.mergeTop
  rwa [hts] at this

theorem CloseInv.finishTop' (h : s.CloseInv) {t : TEntry} {rest : List TEntry}
    (hts : s.tstack = t :: rest) {x : ItemId} (hx : x < s.items.size) (hz : s.cnt x = 0) (dir : Bool)
    (vs : Option Nat × Option Nat)
    (hnew : Items.CloseAt s.g
      (s.items.modify x fun it => { it with vs := vs, ch := getSide t.spans dir }) x) :
    ({ s with
        items := s.items.modify x fun it => { it with vs := vs, ch := getSide t.spans dir },
        tstack := { t with spans := setSides dir [x] [] } :: rest } : WalkState).CloseInv := by
  have := h.finishTop hx hz dir vs (by rw [hts]; exact hnew)
  rw [hts] at this
  exact this

/-- Leaf Q at a returning back edge: `pushEdgeTstack curV lowval o.e` after `feS₀`
(`CloseInv.pushEdge` as for `closeEars_closeAt`; the counter writes are frames). -/
theorem finishBack_closeAt (hc : CloseCtx σ n D curV d o origTstack hasVert s)
    (ht : o.cls.isTree = false) (hlow : o.cls.lowval d < d) :
    (feBack curV (o.cls.lowval d) d o s).CloseInv := by
  obtain ⟨lv, kind, ho, hlv⟩ := ret_of_lowval_lt hlow
  obtain ⟨sub, base, -, hE⟩ := hc.book.ear
  have hends := hc.book.ends lv kind ho
  simp only [ht, Bool.false_eq_true, ↓reduceIte] at hends
  obtain ⟨e, dest, cls, rfl⟩ : ∃ e dest cls, o = .back e dest cls := by
    cases o with
    | back e dest cls => exact ⟨e, dest, cls, rfl⟩
    | tree e cls child => exact absurd (hc.book.tree.mpr ⟨e, cls, child, rfl⟩) (by simp [ht])
  have hends' : Items.PairEq (curV, dest) s.g.edges[e]! := by
    have h := hc.ends; simp only [DfsOut.Ends] at h; exact h
  have hends : Items.PairEq (curV, s.stackVerts[lv]!) s.g.edges[e]! := hends
  have hlvd : s.stackVerts[lv]! ≠ curV := by rw [← hE.sv_d]; exact hE.path lv d hlv (Nat.le_refl d)
  have hdest : dest = s.stackVerts[lv]! := by
    generalize s.g.edges[e]! = q at hends hends'
    obtain ⟨x, y⟩ := q
    simp only [Items.PairEq, Prod.mk.injEq] at hends hends'
    omega
  obtain ⟨e₂, he₂, hn₂, hp₂⟩ := hc.path lv hlv
  have ho' : cls = .ret lv kind := ho
  have hlv' : (DfsOut.back e dest cls).cls.lowval d = lv := by
    show cls.lowval d = lv; rw [ho']; rfl
  rw [hlv']
  have h := leafQ_closeInv hc (v := curV) (dd := lv) (by rw [hE.sv_d]; show curV ≠ dest; rw [hdest]; exact Ne.symm hlvd)
    (by rw [hE.sv_d]; exact hends') (hc.parent_edge (by omega))
    ⟨e₂, he₂, hc.ne_of_path hn₂, by show s.g.Inc e₂ dest; rw [hdest]; exact (Graph.inc_of_pairEq hp₂).1⟩
  exact h.frame rfl rfl fun _ h => h

theorem cnt_lt_size (hs : Shape s) {i : ItemId} (hc : 0 < s.cnt i) : i < s.items.size := by
  by_contra hn
  have h1 : spansCount s.tstack i = 0 := by
    apply List.sum_eq_zero
    intro x hx
    obtain ⟨t, ht, rfl⟩ := List.mem_map.mp hx
    exact List.count_eq_zero.mpr fun hm => hn (hs.span t ht i hm)
  have h2 : chCount s.items i = 0 :=
    Finset.sum_eq_zero fun j _ => List.count_eq_zero.mpr fun hm => hn (hs.ch_lt j i hm)
  have : s.cnt i = 0 := by dsimp [cnt]; omega
  omega

theorem CloseInv.alloc' (h : s.CloseInv) (hs : Shape s) (ty : NodeType) :
    ({ s with items := s.items.push ⟨ty, (none, none), []⟩ } : WalkState).CloseInv := by
  constructor
  intro i _ hc
  have hc' : i ≤ s.g.nv ∨ 0 < s.cnt i := by
    simpa only [cnt, Items.chCount_push s.items ⟨ty, (none, none), []⟩ rfl] using hc
  have hi : i < s.items.size := by
    rcases hc' with hc' | hc'
    · have := hs.size; omega
    · exact cnt_lt_size hs hc'
  apply (h.closed i hi hc').frame
  refine ⟨Items.type_push_of_ne _ (Nat.ne_of_lt hi), Items.vs_push_of_ne _ (Nat.ne_of_lt hi),
    Items.ch_push_nil _ rfl _, ?_, ?_, ?_, ?_⟩
  · intro c hc; exact Items.type_push_of_ne _ (Nat.ne_of_lt (hs.ch_lt i c hc))
  · intro c hc; exact Items.vs_push_of_ne _ (Nat.ne_of_lt (hs.ch_lt i c hc))
  · intro c _; exact Items.ch_push_nil _ rfl _
  · intro j _ e; exact Items.Below_push_nil _ rfl

/-- `finishP` preserves `CloseInv` at a type-1 P site. -/
theorem PSite.closeInv {curV lv : Nat} {r : WalkState} (hp : PSite curV lv r) (h : r.CloseInv)
    (hc : result (condP curV lv true) r = true) :
    (after (finishP curV lv true) r).CloseInv := by
  have h' : (true && decide (r.tstack.length ≥ 2) && (r.tstack.tail.head!.vStart == curV) &&
      (r.tstack.tail.head!.topDepth == lv)) = true := hc
  obtain ⟨cur, rest₀, hts, hcv, hct⟩ := hp.stack
  obtain ⟨nxt, rest, rfl⟩ : ∃ nxt rest, rest₀ = nxt :: rest := by
    cases rest₀ with
    | nil => rw [hts] at h'; simp at h'
    | cons nxt rest => exact ⟨nxt, rest, rfl⟩
  have h'' := h'
  rw [hts] at h''
  have hnv : nxt.vStart = curV :=
    beq_iff_eq.1 (Bool.and_eq_true_iff.1 (Bool.and_eq_true_iff.1 h'').1).2
  have hnt : nxt.topDepth = lv := beq_iff_eq.1 (Bool.and_eq_true_iff.1 h'').2
  have hcur : cur ∈ r.tstack.take 2 := by rw [hts]; simp
  have hnxt : nxt ∈ r.tstack.take 2 := by rw [hts]; simp
  have hcur' : cur ∈ r.tstack := by rw [hts]; simp
  have hnxt' : nxt ∈ r.tstack := by rw [hts]; simp
  obtain ⟨c, hcsp, a₁, b₁, hcvs, hcpe⟩ := hp.single cur hcur
  obtain ⟨i, hisp, a₂, b₂, hivs, hipe⟩ := hp.single nxt hnxt
  have hcmem : c ∈ cur.spans.1 ++ cur.spans.2 := by
    rw [hcsp]; exact (mem_setSides _ _ _).2 (List.mem_singleton.2 rfl)
  have himem : i ∈ nxt.spans.1 ++ nxt.spans.2 := by
    rw [hisp]; exact (mem_setSides _ _ _).2 (List.mem_singleton.2 rfl)
  have hclt : c < r.items.size := hp.shape.span cur hcur' c hcmem
  have hilt : i < r.items.size := hp.shape.span nxt hnxt' i himem
  have hcE : ∀ e, cur.edges r.g r.items e ↔ Items.EdgeBelow r.g r.items c e := fun e =>
    TEntry.edges_single r.stackDir[lv]! c (by rw [hcsp]; exact getSide_setSides_not _ _)
      (by rw [hcsp]; exact getSide_setSides' _ _ _) e
  have hiE : ∀ e, nxt.edges r.g r.items e ↔ Items.EdgeBelow r.g r.items i e := fun e =>
    TEntry.edges_single r.stackDir[lv]! i (by rw [hisp]; exact getSide_setSides_not _ _)
      (by rw [hisp]; exact getSide_setSides' _ _ _) e
  have hattc : ∀ w, Items.Att r.g r.items c w → w = r.stackVerts[lv]! ∨ w = curV := by
    rintro w ⟨e, e', he, he', hi, hi', hb, hnb⟩
    rcases hp.att cur hcur w ⟨e, he, (hcE e).2 hb, hi⟩ with hw | hw | hw
    · exact .inr hw
    · exact .inl hw
    · exact absurd ((hcE e').1 (hw e' he' hi')) hnb
  have hatti : ∀ w, Items.Att r.g r.items i w → w = r.stackVerts[lv]! ∨ w = curV := by
    rintro w ⟨e, e', he, he', hi, hi', hb, hnb⟩
    rcases hp.att nxt hnxt w ⟨e, he, (hiE e).2 hb, hi⟩ with hw | hw | hw
    · exact .inr hw
    · exact .inl hw
    · exact absurd ((hiE e').1 (hw e' he' hi')) hnb
  have hkc := hp.kinds cur hcur c hcmem c (Or.inl rfl)
  have hki := hp.kinds nxt hnxt i himem i (Or.inl rfl)
  have hne : r.stackVerts[lv]! ≠ curV := fun h => hp.ne h.symm
  have htd : (TEntry.mergeInto cur nxt).topDepth = lv := by
    show min nxt.topDepth cur.topDepth = lv
    rw [hnt, hct, Nat.min_self]
  have htv : (TEntry.mergeInto cur nxt).vStart = curV := hnv
  have hi : i = (getSide nxt.spans r.stackDir[lv]!).head! := by
    rw [hisp, getSide_setSides']; rfl
  have hdir : r.stackDir[lv]! = r.stackDir[nxt.topDepth]! := by rw [hnt]
  have htc : ∀ w, w = r.stackVerts[lv]! ∨ w = curV →
      ∃ e, e < r.g.ne ∧ r.g.Inc e w ∧ Items.EdgeBelow r.g r.items c e := by
    intro w hw
    obtain ⟨hcv', hu'⟩ := hp.touch cur hcur
    rcases hw with rfl | rfl
    · obtain ⟨e, he, hE, hi⟩ := hu'; exact ⟨e, he, hi, (hcE e).1 hE⟩
    · obtain ⟨e, he, hE, hi⟩ := hcv'; exact ⟨e, he, hi, (hcE e).1 hE⟩
  have hpd : ∀ w, w = r.stackVerts[lv]! ∨ w = curV → ∃ e, e < r.g.ne ∧ r.g.Inc e w ∧
      ¬ Items.EdgeBelow r.g r.items c e ∧ ¬ Items.EdgeBelow r.g r.items i e := by
    intro w hw
    obtain ⟨e, he, hi, hn⟩ := hp.pend w hw.symm
    exact ⟨e, he, hi, fun hb => hn cur hcur ((hcE e).2 hb), fun hb => hn nxt hnxt ((hiE e).2 hb)⟩
  have hmemL : ∀ (L : List ItemId) j,
      j ∈ (if r.stackDir[lv]! then L ++ [c] else [c] ++ L) ↔ j = c ∨ j ∈ L := by
    intro L j
    cases r.stackDir[lv]! <;> simp [or_comm]
  have hlenL : ∀ L : List ItemId, (if r.stackDir[lv]! then L ++ [c] else [c] ++ L).length =
      L.length + 1 := by
    intro L
    cases r.stackDir[lv]! <;> simp
  show ((finishP curV lv true).run r).2.CloseInv
  simp only [Spqr.finishP, WalkM.run_bind, run_condP]
  simp only [h', ↓reduceIte, WalkM.run_bind]
  rw [maybeUnwrapNxt_run_eq .P r cur nxt rest hts r.stackDir[lv]! hdir i hi]
  have alloc_case : ((finishTstackTop ((allocItem .P).run r).1).run
      (mergeTstackTops.run ((allocItem .P).run r).2).2).2.CloseInv := by
    rw [run_allocItem]
    dsimp only
    rw [mergeTstackTops_run_eq { r with items := r.items.push ⟨.P, (none, none), []⟩ } cur nxt rest
      hts]
    dsimp only
    rw [finishTstackTop_run_eq
      { r with
          items := r.items.push ⟨.P, (none, none), []⟩,
          tstack := TEntry.mergeInto cur nxt :: rest }
      r.items.size (TEntry.mergeInto cur nxt) rest rfl]
    dsimp only
    have hroot₁ : ∀ p, ¬ Items.IsParent (r.items.push ⟨.P, (none, none), []⟩) p r.items.size := by
      intro p hpar
      have hpar' : r.items.size ∈ Items.ch (r.items.push ⟨.P, (none, none), []⟩) p := hpar
      rw [Items.ch_push_nil _ rfl] at hpar'
      exact absurd (hp.shape.ch_lt p _ hpar') (lt_irrefl _)
    have hx₁ : r.items.size < (r.items.push ⟨.P, (none, none), []⟩).size := by simp
    have h₂ : ({ r with
        items := r.items.push ⟨.P, (none, none), []⟩,
        tstack := TEntry.mergeInto cur nxt :: rest } : WalkState).CloseInv :=
      (h.alloc' hp.shape .P).mergeTop' hts
    refine CloseInv.finishTop' h₂ rfl hx₁ ?_ _ _ ?_
    · refine cnt_eq_zero_of_free ?_ hroot₁
      intro t ht hmem
      have ht' : t ∈ TEntry.mergeInto cur nxt :: rest := ht
      rcases List.mem_cons.1 ht' with rfl | ht'
      · rcases (TEntry.mem_mergeInto cur nxt _).1 hmem with hm | hm
        · exact absurd (hp.shape.span cur hcur' _ hm) (lt_irrefl _)
        · exact absurd (hp.shape.span nxt hnxt' _ hm) (lt_irrefl _)
      · have ht'' : t ∈ r.tstack := by
          rw [hts]; exact List.mem_cons_of_mem _ (List.mem_cons_of_mem _ ht')
        exact absurd (hp.shape.span t ht'' _ hmem) (lt_irrefl _)
    · dsimp only
      rw [htd, htv, getSide_mergeInto hcsp hisp]
      have hEB : ∀ j, j ≠ r.items.size → ∀ e,
          Items.EdgeBelow r.g ((r.items.push ⟨.P, (none, none), []⟩).modify r.items.size fun it =>
            { it with vs := setSides r.stackDir[lv]! (some r.stackVerts[lv]!) (some curV),
                      ch := if r.stackDir[lv]! then [i] ++ [c] else [c] ++ [i] }) j e ↔
          Items.EdgeBelow r.g r.items j e := fun j hj e =>
        (Items.Below_modify_of_not_below _ _ fun hb => hj (Items.Below.eq_of_no_parent hroot₁ hb)).trans
          (Items.Below_push_nil _ rfl)
      have hAtt : ∀ j, j ≠ r.items.size → ∀ w,
          Items.Att r.g ((r.items.push ⟨.P, (none, none), []⟩).modify r.items.size fun it =>
            { it with vs := setSides r.stackDir[lv]! (some r.stackVerts[lv]!) (some curV),
                      ch := if r.stackDir[lv]! then [i] ++ [c] else [c] ++ [i] }) j w →
          Items.Att r.g r.items j w := by
        rintro j hj w ⟨e, e', he, he', hi, hi', hb, hnb⟩
        exact ⟨e, e', he, he', hi, hi', (hEB j hj e).1 hb, fun hb' => hnb ((hEB j hj e').2 hb')⟩
      have hch : Items.ch ((r.items.push ⟨.P, (none, none), []⟩).modify r.items.size fun it =>
            { it with vs := setSides r.stackDir[lv]! (some r.stackVerts[lv]!) (some curV),
                      ch := if r.stackDir[lv]! then [i] ++ [c] else [c] ++ [i] }) r.items.size =
          if r.stackDir[lv]! then [i] ++ [c] else [c] ++ [i] := Items.ch_modify_at _ _ hx₁
      have hf : ∀ it : Item,
          ({ it with
              vs := setSides r.stackDir[lv]! (some r.stackVerts[lv]!) (some curV),
              ch := if r.stackDir[lv]! then [i] ++ [c] else [c] ++ [i] } : Item).type = it.type :=
        fun _ => rfl
      have hpar : ∀ j, j = c ∨ j = i → Items.IsParent
          ((r.items.push ⟨.P, (none, none), []⟩).modify r.items.size fun it =>
            { it with vs := setSides r.stackDir[lv]! (some r.stackVerts[lv]!) (some curV),
                      ch := if r.stackDir[lv]! then [i] ++ [c] else [c] ++ [i] }) r.items.size j := by
        intro j hj
        show j ∈ Items.ch _ _
        rw [hch]
        exact (hmemL [i] j).2 (by simpa using hj)
      have hcne : c ≠ r.items.size := Nat.ne_of_lt hclt
      have hine : i ≠ r.items.size := Nat.ne_of_lt hilt
      refine Items.CloseAt.pNode' (u := r.stackVerts[lv]!) (v := curV) r.stackDir[lv]! hp.shape.size
        ?_ hch (Items.vs_modify_at _ _ hx₁) hne ?_ ?_ ?_ ?_ ?_ ?_ ?_
      · rw [Items.type_modify_type_eq _ _ hf, Items.type_push_size]
      · rw [hlenL]; simp
      · intro w hw
        rw [Items.type_modify_type_eq _ _ hf,
          Items.type_push_of_ne _ (Nat.ne_of_lt (by have := hp.shape.size; show 1 + w < _; omega))]
        exact hp.shape.vert w hw
      · intro j hj
        rw [Items.type_modify_type_eq _ _ hf]
        rcases (hmemL [i] j).1 hj with rfl | hj
        · rw [Items.type_push_of_ne _ hcne]; exact hkc.1
        · rw [List.mem_singleton.1 hj, Items.type_push_of_ne _ hine]; exact hki.1
      · intro j hj
        rcases (hmemL [i] j).1 hj with rfl | hj
        · rw [Items.vs_modify_of_ne _ _ hcne, Items.vs_push_of_ne _ hcne]; exact ⟨a₁, b₁, hcvs, hcpe⟩
        · rw [List.mem_singleton.1 hj, Items.vs_modify_of_ne _ _ hine, Items.vs_push_of_ne _ hine]
          exact ⟨a₂, b₂, hivs, hipe⟩
      · intro j hj w hw
        rcases (hmemL [i] j).1 hj with rfl | hj
        · exact hattc w (hAtt _ hcne w hw)
        · rw [List.mem_singleton.1 hj] at hw; exact hatti w (hAtt _ hine w hw)
      · intro w hw
        obtain ⟨e, he, hi, hb⟩ := htc w hw
        exact ⟨e, he, hi, Relation.ReflTransGen.head (hpar c (Or.inl rfl)) ((hEB c hcne e).2 hb)⟩
      · intro w hw
        obtain ⟨e, he, hi, hnc, hni⟩ := hpd w hw
        refine ⟨e, he, hi, fun hb => ?_⟩
        rcases hb.head_cases with heq | ⟨j, hj, hb'⟩
        · have : r.items.size = 1 + r.g.nv + e := heq
          have := hp.shape.size
          omega
        · have hj' : j ∈ Items.ch _ _ := hj
          rw [hch] at hj'
          rcases (hmemL [i] j).1 hj' with rfl | hj'
          · exact hnc ((hEB _ hcne e).1 hb')
          · rw [List.mem_singleton.1 hj'] at hb'
            exact hni ((hEB _ hine e).1 hb')
  by_cases htern : r.ternarize = true
  · rw [if_pos (Or.inr htern)]; exact alloc_case
  rw [if_neg (by simp [htern])]
  by_cases hty : r.items[i]!.type = .P
  swap
  · rw [if_neg hty]; exact alloc_case
  rw [if_pos hty]
  dsimp only
  obtain ⟨nxt', hn⟩ : ∃ nxt' : TEntry,
      nxt' = { nxt with spans := setSides r.stackDir[lv]! r.items[i]!.ch [] } := ⟨_, rfl⟩
  rw [← hn]
  have hchi : r.items[i]!.ch = Items.ch r.items i := by simp [Items.ch_eq_getElem hilt, hilt]
  have htyi : Items.type r.items i = .P := by rw [Items.type_eq_getElem hilt]; simpa [hilt] using hty
  obtain ⟨hroot, hsc⟩ := hp.once nxt hnxt i himem
  have hfresh : ∀ t ∈ cur :: rest, i ∉ t.spans.1 ++ t.spans.2 := by
    rw [hts, spansCount_cons, spansCount_cons] at hsc
    have hpos := List.count_pos_iff.2 himem
    intro t ht hmem
    have hpos' := List.count_pos_iff.2 hmem
    rcases List.mem_cons.1 ht with rfl | ht
    · omega
    · have : (t.spans.1 ++ t.spans.2).count i ≤ spansCount rest i :=
        List.single_le_sum (fun _ _ => Nat.zero_le _) _
          (List.mem_map_of_mem (f := fun t : TEntry => (t.spans.1 ++ t.spans.2).count i) ht)
      omega
  have hcni : c ≠ i := fun h => hfresh cur (List.mem_cons_self ..) (h ▸ hcmem)
  have hchni : ∀ j ∈ Items.ch r.items i, j ≠ i := fun j hj h => hroot i (h ▸ hj)
  have hPi : Items.CloseAt r.g r.items i := by
    refine h.closed i hilt (Or.inr ?_)
    have := spansCount_pos_of_mem_tail_head! (ts := r.tstack) (c := i)
      (by show i ∈ r.tstack.tail.head!.spans.1 ++ _; rw [hts]; exact himem)
    dsimp [cnt]; omega
  have hps := hPi.p_shape htyi
  have hchlen : 2 ≤ (Items.ch r.items i).length := le_trans hps.1 (by
    unfold Items.virtualEdges; rw [List.length_map]; exact List.length_filter_le _ _)
  have hchild : ∀ j ∈ Items.ch r.items i, Items.type r.items j ∈ [NodeType.S, .P, .R, .Q] ∧
      (∃ a b, Items.vs r.items j = (some a, some b) ∧
        Items.PairEq (a, b) (r.stackVerts[lv]!, curV)) ∧
      ∀ w, Items.Att r.g r.items j w → w = r.stackVerts[lv]! ∨ w = curV := by
    intro j hj
    have hk := hp.kinds nxt hnxt i himem j (Or.inr ⟨htyi, hj⟩)
    obtain ⟨a, b, hab, -⟩ := hPi.child_two j hj (by simp [htyi]) (hps.2.1 j hj)
    have hq : (a, b) ∈ Items.virtualEdges r.items i :=
      List.mem_map.2 ⟨j, List.mem_filter.2 ⟨hj, decide_eq_true (hps.2.1 j hj)⟩, by rw [hab]; rfl⟩
    obtain ⟨u', v', hu'v', hpe⟩ := hps.2.2 (a, b) hq
    rw [hivs] at hu'v'
    simp only [Prod.mk.injEq, Option.some.injEq] at hu'v'
    obtain ⟨rfl, rfl⟩ := hu'v'
    have hpe' := hpe.trans' hipe
    refine ⟨hk.1, ⟨a, b, hab, hpe'⟩, fun w hw => ?_⟩
    have hjlt := hp.shape.ch_lt i j hj
    have hpos : 0 < r.cnt j := by
      have := Items.count_le_chCount r.items hilt j
      have := List.count_pos_iff.mpr hj
      dsimp [cnt]; omega
    have hv := (h.closed j hjlt (Or.inr hpos)).att_vs hk.1 w hw
    simp only [Items.IsVs, hab, Option.some.injEq] at hv
    simp [Items.PairEq] at hpe'
    omega
  have hn' : nxt'.spans = setSides r.stackDir[lv]! (Items.ch r.items i) [] := by
    rw [hn, ← hchi]
  have htd' : (TEntry.mergeInto cur nxt').topDepth = lv := by
    show min nxt'.topDepth cur.topDepth = lv
    rw [hn]; exact htd
  have htv' : (TEntry.mergeInto cur nxt').vStart = curV := by
    show nxt'.vStart = curV
    rw [hn]; exact hnv
  rw [mergeTstackTops_run_eq { r with tstack := cur :: nxt' :: rest } cur nxt' rest rfl]
  dsimp only
  rw [finishTstackTop_run_eq { r with tstack := TEntry.mergeInto cur nxt' :: rest } i
    (TEntry.mergeInto cur nxt') rest rfl]
  dsimp only
  have h₁ : ({ r with tstack := cur :: nxt' :: rest } : WalkState).CloseInv := by
    have := h.unwrap hts hilt r.stackDir[lv]!
    rw [← hchi] at this
    rw [hn]; exact this
  have h₂ : ({ r with tstack := TEntry.mergeInto cur nxt' :: rest } : WalkState).CloseInv :=
    h₁.mergeTop' rfl
  refine CloseInv.finishTop' h₂ rfl hilt ?_ _ _ ?_
  · refine cnt_eq_zero_of_free ?_ hroot
    intro t ht hmem
    have ht' : t ∈ TEntry.mergeInto cur nxt' :: rest := ht
    rcases List.mem_cons.1 ht' with rfl | ht'
    · rcases (TEntry.mem_mergeInto _ _ _).1 hmem with hm | hm
      · exact hfresh cur (List.mem_cons_self ..) hm
      · rw [hn', mem_setSides] at hm
        exact hroot i hm
    · exact hfresh t (List.mem_cons_of_mem _ ht') hmem
  · dsimp only
    rw [htd', htv', getSide_mergeInto hcsp hn']
    have hEB : ∀ j, j ≠ i → ∀ e,
        Items.EdgeBelow r.g (r.items.modify i fun it =>
          { it with vs := setSides r.stackDir[lv]! (some r.stackVerts[lv]!) (some curV),
                    ch := if r.stackDir[lv]! then Items.ch r.items i ++ [c]
                      else [c] ++ Items.ch r.items i }) j e ↔
        Items.EdgeBelow r.g r.items j e := fun j hj e =>
      Items.Below_modify_of_not_below _ _ fun hb => hj (Items.Below.eq_of_no_parent hroot hb)
    have hAtt : ∀ j, j ≠ i → ∀ w,
        Items.Att r.g (r.items.modify i fun it =>
          { it with vs := setSides r.stackDir[lv]! (some r.stackVerts[lv]!) (some curV),
                    ch := if r.stackDir[lv]! then Items.ch r.items i ++ [c]
                      else [c] ++ Items.ch r.items i }) j w →
        Items.Att r.g r.items j w := by
      rintro j hj w ⟨e, e', he, he', hi, hi', hb, hnb⟩
      exact ⟨e, e', he, he', hi, hi', (hEB j hj e).1 hb, fun hb' => hnb ((hEB j hj e').2 hb')⟩
    have hch : Items.ch (r.items.modify i fun it =>
          { it with vs := setSides r.stackDir[lv]! (some r.stackVerts[lv]!) (some curV),
                    ch := if r.stackDir[lv]! then Items.ch r.items i ++ [c]
                      else [c] ++ Items.ch r.items i }) i =
        if r.stackDir[lv]! then Items.ch r.items i ++ [c] else [c] ++ Items.ch r.items i :=
      Items.ch_modify_at _ _ hilt
    have hf : ∀ it : Item,
        ({ it with
            vs := setSides r.stackDir[lv]! (some r.stackVerts[lv]!) (some curV),
            ch := if r.stackDir[lv]! then Items.ch r.items i ++ [c]
              else [c] ++ Items.ch r.items i } : Item).type = it.type :=
      fun _ => rfl
    have hpar : ∀ j, j = c ∨ j ∈ Items.ch r.items i → Items.IsParent
        (r.items.modify i fun it =>
          { it with vs := setSides r.stackDir[lv]! (some r.stackVerts[lv]!) (some curV),
                    ch := if r.stackDir[lv]! then Items.ch r.items i ++ [c]
                      else [c] ++ Items.ch r.items i }) i j := by
      intro j hj
      show j ∈ Items.ch _ _
      rw [hch]
      exact (hmemL _ j).2 hj
    refine Items.CloseAt.pNode' (u := r.stackVerts[lv]!) (v := curV) r.stackDir[lv]!
      (hp.shape.node_of_type hilt (by simp) htyi) ?_ hch (Items.vs_modify_at _ _ hilt) hne ?_ ?_ ?_ ?_ ?_ ?_ ?_
    · rw [Items.type_modify_type_eq _ _ hf]; exact htyi
    · rw [hlenL]; omega
    · intro w hw
      rw [Items.type_modify_type_eq _ _ hf]
      exact hp.shape.vert w hw
    · intro j hj
      rw [Items.type_modify_type_eq _ _ hf]
      rcases (hmemL _ j).1 hj with rfl | hj
      · exact hkc.1
      · exact (hchild j hj).1
    · intro j hj
      rcases (hmemL _ j).1 hj with rfl | hj
      · rw [Items.vs_modify_of_ne _ _ hcni]; exact ⟨a₁, b₁, hcvs, hcpe⟩
      · rw [Items.vs_modify_of_ne _ _ (hchni j hj)]; exact (hchild j hj).2.1
    · intro j hj w hw
      rcases (hmemL _ j).1 hj with rfl | hj
      · exact hattc w (hAtt _ hcni w hw)
      · exact (hchild j hj).2.2 w (hAtt _ (hchni j hj) w hw)
    · intro w hw
      obtain ⟨e, he, hi, hb⟩ := htc w hw
      exact ⟨e, he, hi, Relation.ReflTransGen.head (hpar c (Or.inl rfl)) ((hEB c hcni e).2 hb)⟩
    · intro w hw
      obtain ⟨e, he, hi, hnc, hni⟩ := hpd w hw
      refine ⟨e, he, hi, fun hb => ?_⟩
      rcases hb.head_cases with heq | ⟨j, hj, hb'⟩
      · have h3 : i = 1 + r.g.nv + e := heq
        have h4 := hp.shape.node_of_type hilt (by simp) htyi
        rw [h3] at h4
        omega
      · have hj' : j ∈ Items.ch _ _ := hj
        rw [hch] at hj'
        rcases (hmemL _ j).1 hj' with rfl | hj'
        · exact hnc ((hEB _ hcni e).1 hb')
        · exact hni (Relation.ReflTransGen.head hj' ((hEB _ (hchni j hj') e).1 hb'))

/-- The `some item` arm of `closeVertTail` (`vertFinish`): the S or R record of the vertex ear,
from the state `cvS₅` after the merges and the retarget (the `none` arm is the identity). -/
theorem VSite.closeInv {curV : Nat} {x : ItemId} {t : TEntry} {r : WalkState}
    (hv : VSite curV x t r) (h : r.CloseInv) : (after (finishTstackTop x) r).CloseInv := by
  obtain ⟨rest, hts⟩ := hv.stack
  show ((finishTstackTop x).run r).2.CloseInv
  rw [finishTstackTop_run_eq r x t rest hts]
  dsimp only
  refine CloseInv.finishTop' h hts hv.free.lt (cnt_eq_zero_of_free hv.free.free hv.free.root) _ _ ?_
  obtain ⟨items', hI⟩ : ∃ I : Items, I = r.items.modify x fun it =>
      { it with vs := setSides r.stackDir[t.topDepth]! (some r.stackVerts[t.topDepth]!) (some t.vStart),
                ch := getSide t.spans r.stackDir[t.topDepth]! } := ⟨_, rfl⟩
  rw [← hI]
  have hf : ∀ it : Item,
      ({ it with
          vs := setSides r.stackDir[t.topDepth]! (some r.stackVerts[t.topDepth]!) (some t.vStart),
          ch := getSide t.spans r.stackDir[t.topDepth]! } : Item).type = it.type := fun _ => rfl
  have htc : ∀ c, Items.type items' c = Items.type r.items c := fun c => by
    rw [hI]; exact Items.type_modify_type_eq _ _ hf c
  have hch : Items.ch items' x = vKids t r := by rw [hI]; exact Items.ch_modify_at _ _ hv.free.lt
  have hvs : Items.vs items' x =
      (some (vTerms curV t r).1, some (vTerms curV t r).2) := by
    rw [hI, Items.vs_modify_at _ _ hv.free.lt]
    show setSides _ (some _) (some t.vStart) = _
    rw [hv.vstart]; exact setSides_some _ _ _
  have hxs : x ∉ t.spans.1 ++ t.spans.2 := hv.free.free t (by rw [hts]; simp)
  have hxne : ∀ c ∈ vKids t r, c ≠ x := fun c hc heq => hxs (heq ▸ mem_of_mem_getSide hc)
  have hvsc : ∀ c ∈ vKids t r, Items.vs items' c = Items.vs r.items c := fun c hc => by
    rw [hI]; exact Items.vs_modify_of_ne _ _ (hxne c hc)
  have hEB : ∀ j, j ≠ x → ∀ e, Items.EdgeBelow r.g items' j e ↔ Items.EdgeBelow r.g r.items j e :=
    fun j hj e => by
      rw [hI]
      exact Items.Below_modify_of_not_below _ _ fun hb => hj (Items.Below.eq_of_no_parent hv.free.root hb)
  have hpar : ∀ c, Items.IsParent items' x c ↔ c ∈ vKids t r := fun c => by
    show c ∈ Items.ch items' x ↔ _
    rw [hch]
  have hEx : ∀ e, e < r.g.ne → (Items.EdgeBelow r.g items' x e ↔ t.edges r.g r.items e) := by
    intro e he
    constructor
    · intro hb
      rcases hb.head_cases with heq | ⟨c, hc, hb'⟩
      · have : x = 1 + r.g.nv + e := heq
        have h4 := hv.free.node
        rw [this] at h4
        omega
      · have hc' := (hpar c).1 hc
        exact ⟨c, mem_of_mem_getSide hc', (hEB c (hxne c hc') e).1 hb'⟩
    · rintro ⟨i, hi, hb⟩
      have hi' := (mem_sides_iff hv.side).1 hi
      exact Relation.ReflTransGen.head ((hpar i).2 hi') ((hEB i (hxne i hi') e).2 hb)
  have hor : ∀ w, (w = (vTerms curV t r).1 ∨ w = (vTerms curV t r).2) ↔
      (w = curV ∨ w = r.stackVerts[t.topDepth]!) := fun w => by
    rw [vTerms, eq_setSides_or, or_comm]
  have hfilt : ((vKids t r).filter fun c => Items.type items' c = .V) =
      (vKids t r).filter fun c => Items.type r.items c = .V :=
    List.filter_congr fun c _ => by rw [htc c]
  have hVE : Items.virtualEdges items' x = vVirt t r := by
    unfold Items.virtualEdges vVirt
    rw [hch, List.filter_congr (fun c _ => by rw [htc c])]
    exact List.map_congr_left fun c hc => by rw [hvsc c (List.mem_of_mem_filter hc)]
  refine Items.CloseAt.node hv.free.node (by rw [htc]; exact hv.ty) hch hvs
    (setSides_ne _ hv.ne.symm) (fun c hc => by rw [htc]; exact hv.kinds c hc)
    (fun c hc hne => by rw [htc] at hne; rw [hvsc c hc]; exact hv.two c hc hne) ?_ ?_ ?_ ?_ ?_ ?_
  · rintro w ⟨e, e', he, he', hi, hi', hb, hnb⟩
    refine (hor w).2 ?_
    rcases hv.att w ⟨e, he, (hEx e he).1 hb, hi⟩ with h1 | h1 | h1
    · exact .inl h1
    · exact .inr h1
    · exact (hnb ((hEx e' he').2 (h1 e' he' hi'))).elim
  · intro w hw
    rcases (hor w).1 hw with rfl | rfl
    · obtain ⟨e, he, hE, hi⟩ := hv.touch.1; exact ⟨e, he, hi, (hEx e he).2 hE⟩
    · obtain ⟨e, he, hE, hi⟩ := hv.touch.2; exact ⟨e, he, hi, (hEx e he).2 hE⟩
  · intro w hw
    obtain ⟨e, he, hi, hne⟩ := hv.pend w ((hor w).1 hw)
    exact ⟨e, he, hi, fun hb => hne ((hEx e he).1 hb)⟩
  · intro w hw
    rw [hv.inner w hw]
    constructor
    · rintro ⟨⟨e, he, hE, hi⟩, hIn, hno⟩
      refine ⟨⟨⟨e, he, hi⟩, fun e' he' hi' => (hEx e' he').2 (hIn e' he' hi')⟩, fun c hc hall =>
        hno c hc fun e' he' hi' => (hEB c (hxne c hc) e').1 (hall e' he' hi')⟩
    · rintro ⟨⟨⟨e, he, hi⟩, hall⟩, hno⟩
      refine ⟨⟨e, he, (hEx e he).1 (hall e he hi), hi⟩, fun e' he' hi' => (hEx e' he').1 (hall e' he' hi'),
        fun c hc hIn => hno c hc fun e' he' hi' => (hEB c (hxne c hc) e').2 (hIn e' he' hi')⟩
  · intro hS
    rw [htc] at hS
    obtain ⟨xs, hf, hl, hve⟩ := hv.s_order hS
    exact ⟨xs, by rw [hfilt]; exact hf, hl, by rw [hVE]; exact hve⟩
  · intro hR
    rw [htc] at hR
    obtain ⟨h1, h2, h3, h4⟩ := hv.r_shape hR
    rw [hfilt, hVE]
    exact ⟨h1, h2, h3, h4⟩

theorem closeVertTail_closeAt (hc : CloseCtx σ n D curV d o origTstack hasVert s)
    (ht : o.cls.isTree = true) (hlow : o.cls.lowval d < d) (hv : hasVert = true)
    (h : (cvS₅ curV s.stackDir[d]! o.cls.isType1 origTstack (feSingle d o s) (feS₂ d o s)).CloseInv) :
    (feS₃ curV d o origTstack s).CloseInv := by
  cases h1 : o.cls.isType1
  · have hS₃ : feS₃ curV d o origTstack s =
        cvS₅ curV s.stackDir[d]! false origTstack (feSingle d o s) (feS₂ d o s) := by
      simp only [feS₃, h1, closeVert', cvS₅, cvS₄, cvS₃, cvS₂, cvS₁, after, vertPre,
        vertUnwrap, vertFinish, Bool.not_false, Bool.false_eq_true, ↓reduceIte, WalkM.run_bind,
        WalkM.pure_run, run_tstackSize]
    rw [hS₃]; rw [h1] at h; exact h
  · have hS₃ : feS₃ curV d o origTstack s =
        after (finishTstackTop ((maybeUnwrapNxt (if feSingle d o s then .S else .R)).run (feS₂ d o s)).1)
          (cvS₅ curV s.stackDir[d]! true origTstack (feSingle d o s) (feS₂ d o s)) := by
      simp only [feS₃, h1, closeVert', cvS₅, cvS₄, cvS₃, cvS₂, cvS₁, cvB₁, after, result, vertPre,
        vertUnwrap, vertFinish, Bool.not_true, Bool.false_eq_true, ↓reduceIte, WalkM.run_bind,
        WalkM.pure_run, WalkM.map_run]
      rfl
    rw [hS₃]; rw [h1] at h
    obtain ⟨t, hvs⟩ := hc.v_site ht hlow hv h1
    exact hvs.closeInv h

theorem finishP_closeInv_of_not {curV lv : Nat} {b : Bool} (h : s.CloseInv)
    (hc : result (condP curV lv b) s = false) : (after (finishP curV lv b) s).CloseInv := by
  have hc' : (b && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == curV) &&
      (s.tstack.tail.head!.topDepth == lv)) = false := hc
  show ((finishP curV lv b).run s).2.CloseInv
  simp only [Spqr.finishP, WalkM.run_bind, run_condP]
  simp only [hc', Bool.false_eq_true, ↓reduceIte]
  exact h

theorem finishP_closeAt (hc : CloseCtx σ n D curV d o origTstack hasVert s)
    (hlow : o.cls.lowval d < d) (h : (feRest curV d o origTstack hasVert s).CloseInv) :
    (after (finishP curV (o.cls.lowval d) o.cls.isType1)
      (feRest curV d o origTstack hasVert s)).CloseInv := by
  by_cases hcond : result (condP curV (o.cls.lowval d) o.cls.isType1)
      (feRest curV d o origTstack hasVert s) = true
  · have ht1 : o.cls.isType1 = true := by
      have h5 : (o.cls.isType1 && decide ((feRest curV d o origTstack hasVert s).tstack.length ≥ 2) &&
          ((feRest curV d o origTstack hasVert s).tstack.tail.head!.vStart == curV) &&
          ((feRest curV d o origTstack hasVert s).tstack.tail.head!.topDepth == o.cls.lowval d)) =
          true := hcond
      simp only [Bool.and_eq_true] at h5
      exact h5.1.1.1
    rw [ht1] at hcond ⊢
    exact (hc.p_site hlow ht1 hcond).closeInv h hcond
  · exact finishP_closeInv_of_not h (Bool.eq_false_iff.2 hcond)

/-- `CloseInv.vertex_append` from `Shape` alone: the vertex record only needs its children to be
nonempty Qs, and every existing child is a counted item. -/
theorem CloseInv.vertex_append' (h : s.CloseInv) (hs : Shape s) (v e : Nat) (hv : v < s.g.nv)
    (he : e < s.g.ne) (hz : s.cnt (vertItem v) = 0)
    (hq : Items.CloseAt s.g s.items (edgeItem s.g e))
    (hqne : Items.ch s.items (edgeItem s.g e) ≠ []) :
    ({ s with items := s.items.modify (vertItem v) fun it =>
      { it with ch := it.ch ++ [edgeItem s.g e] } } : WalkState).CloseInv := by
  have hvlt : vertItem v < s.items.size := by have := hs.size; show 1 + v < _; omega
  have hn : Items.NoParent s.items (vertItem v) := noParent_of_cnt_eq_zero hz
  have hvert := h.closed (vertItem v) hvlt (Or.inl (by show 1 + v ≤ s.g.nv; omega))
  have htv := hs.vert v hv
  have heq : (s.items.modify (vertItem v) fun it => { it with ch := it.ch ++ [edgeItem s.g e] }) =
      s.items.modify (vertItem v) fun it =>
        { it with ch := Items.ch s.items (vertItem v) ++ [edgeItem s.g e] } := by
    apply Array.ext
    · simp
    · intro i hi₁ hi₂
      simp only [Array.getElem_modify]
      split
      · next hEq => subst i; simp [Items.ch, Array.getElem?_eq_getElem hvlt]
      · rfl
  rw [heq]
  have hch : ∀ c, Items.IsParent (s.items.modify (vertItem v) fun it =>
      { it with ch := Items.ch s.items (vertItem v) ++ [edgeItem s.g e] }) (vertItem v) c ↔
      c ∈ Items.ch s.items (vertItem v) ∨ c = edgeItem s.g e := by
    intro c
    simp only [Items.IsParent, Items.ch_modify_self _ _ _ hvlt, List.mem_append, List.mem_singleton]
  have hcv : ∀ c, c ∈ Items.ch s.items (vertItem v) ∨ c = edgeItem s.g e → vertItem v ≠ c := by
    rintro c (hc | rfl) heq
    · subst heq; exact hn _ hc
    · have : 1 + v = 1 + s.g.nv + e := heq; omega
  apply h.writeChildren hvlt hn
  · apply Items.CloseAt.vertex' hv
    · rw [Items.type_modify_ch]; exact htv
    · intro c hc
      have hc' := (hch c).1 hc
      rw [Items.type_modify_ch]
      rcases hc' with hc' | rfl
      · constructor <;> intro hio <;>
          · have := hvert.io_parent c hc' (by simp [hio]); rw [htv] at this; cases this
      · rw [hs.edge e he]; simp
    · intro c hc
      have hc' := (hch c).1 hc
      rw [Items.ch_modify_ne _ _ _ _ (hcv c hc')]
      rcases hc' with hc' | rfl
      · exact hvert.q_under_v v rfl hv c hc'
      · exact hqne
  · intro c hc
    rcases List.mem_append.mp hc with hc | hc
    · have hpos : 0 < s.cnt c := by
        have := Items.count_le_chCount s.items hvlt c
        have := List.count_pos_iff.mpr hc
        dsimp [cnt]; omega
      exact h.closed c (cnt_lt_size hs hpos) (Or.inr hpos)
    · obtain rfl := List.mem_singleton.mp hc; exact hq

theorem eq_of_inc_pairEq {g : Graph} {e a b v : Nat} (hp : Items.PairEq (a, b) g.edges[e]!)
    (hv : g.Inc e v) : v = a ∨ v = b := by
  have hv' : (g.edges[e]!).1 = v ∨ (g.edges[e]!).2 = v := hv
  generalize g.edges[e]! = q at hp hv'
  obtain ⟨x, y⟩ := q
  simp only [Items.PairEq, Prod.mk.injEq] at hp hv'
  omega

theorem CloseCtx.ends_tree (hc : CloseCtx σ n D curV d o origTstack hasVert s)
    (hT : o.cls.isTree = true) : Items.PairEq (o.dest, curV) s.g.edges[o.e]! := by
  obtain ⟨e, cls, child, rfl⟩ := hc.book.tree.mp hT
  have h := hc.ends
  unfold DfsOut.Ends at h
  exact h.1

theorem CloseCtx.ends_back (hc : CloseCtx σ n D curV d o origTstack hasVert s)
    (hT : o.cls.isTree = false) : Items.PairEq (curV, o.dest) s.g.edges[o.e]! := by
  cases o with
  | back e dest cls => have h := hc.ends; unfold DfsOut.Ends at h; exact h
  | tree e cls child => exact absurd (hc.book.tree.mpr ⟨e, cls, child, rfl⟩) (by simp [hT])

/-- `finishBoundary`: the nonempty Q record of `o.e` with children `I :: t.spans.2` (bridge),
`backedge.spans.1 ++ t.spans.2` (completed block), or `[O]` (self-loop), then appended under
`vertItem curV` (`CloseInv.vertex_append'`). -/
theorem finishBoundary_closeAt (hc : CloseCtx σ n D curV d o origTstack hasVert s)
    (hge : d ≤ o.cls.lowval d) :
    (after (finishEdge curV d o origTstack hasVert) s).CloseInv := by
  obtain ⟨sub, base, hlen, hE⟩ := hc.book.ear
  have hok := ear_boundary hge hc.guards hc.ranges.inv hc.shape hc.book hc.hD
  have hs := hc.shape
  have hi := hc.ranges.inv
  have hge' : o.cls.lowval d ≥ d := hge
  have helt := hc.book.e_lt
  have hvlt := hc.book.v_lt
  have hqlt : edgeItem s.g o.e < 1 + s.g.nv + s.g.ne := by show 1 + s.g.nv + o.e < _; omega
  have hqsz : edgeItem s.g o.e < s.items.size := Nat.lt_of_lt_of_le hqlt hs.size
  have hvsz : vertItem curV < s.items.size := by show 1 + curV < _; have := hs.size; omega
  set s₀ := { s with items := s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) } }
    with hs₀
  have b₀ : BStep D s s₀ := BStep.modifyVs hi hs _ _ hqlt
  have h₀ : s₀.CloseInv := hc.close.modifyLoose (x := edgeItem s.g o.e)
    (by show s.g.nv < 1 + s.g.nv + o.e; omega) (cnt_eq_zero_of_free hE.q_free hE.q_root) _
    (fun _ => rfl)
  set s₁ := { s₀ with totBlocks := s₀.totBlocks + 1 } with hs₁
  have b₁ : BStep D s s₁ := b₀.trans (BStep.frame b₀.inv b₀.shape)
  have h₁ : s₁.CloseInv := h₀.frame rfl rfl fun _ h => h
  have hsz₁ : s₁.items.size = s.items.size := by simp [hs₁, hs₀]
  have hq_root₁ : ∀ p, ¬ Items.IsParent s₁.items p (edgeItem s.g o.e) := fun p h =>
    hok.q_root p ((isParent_modifyVs_iff _ _ _ _ _).1 h)
  have hv_root₁ : ∀ p, ¬ Items.IsParent s₁.items p (vertItem curV) := fun p h =>
    hok.v_root p ((isParent_modifyVs_iff _ _ _ _ _).1 h)
  have hB₁ : ∀ a i, Items.Below s₁.items a i ↔ Items.Below s.items a i := fun a i =>
    Items.Below_modify_vs _ _ _ _
  have hty₁ : ∀ j, Items.type s₁.items j = Items.type s.items j := fun j =>
    Items.type_modify_vs _ _ _
  have hch₁ : ∀ j, Items.ch s₁.items j = Items.ch s.items j := fun j =>
    Items.ch_modify_vs _ _ _
  have hvs₁ : Items.vs s₁.items (edgeItem s.g o.e) = (some curV, none) := by
    show Items.vs (s.items.modify _ _) _ = _
    rw [Items.vs_modify_self _ _ _ hqsz]
  have hvs₁' : ∀ j, j ≠ edgeItem s.g o.e → Items.vs s₁.items j = Items.vs s.items j := fun j hj =>
    Items.vs_modify_ne _ _ _ _ (Ne.symm hj)
  have hsv : s.stackVerts[d]! = curV := hE.sv_d
  show wp (finishEdge curV d o origTstack hasVert) (fun _ s' => s'.CloseInv) s
  rw [finishEdge_eq]
  simp only [finishEdge', wp_bind, wp_get, wp_stackDir, hge', ↓reduceIte]
  unfold finishBoundary
  simp only [wp_bind, wp_modifyItem, wp_modify, wp_ite, wp_allocItem, wp_makeVs, wp_popTstack,
    wp_pure]
  split
  · rename_i hT
    have hends := hc.ends_tree hT
    have hdne : o.dest ≠ curV := fun h =>
      hE.path_child hT d (Nat.le_refl d) (by rw [hsv, h])
    have hedge : (s.g.edges[o.e]!).1 ≠ (s.g.edges[o.e]!).2 := by
      intro h
      generalize s.g.edges[o.e]! = q at hends h
      obtain ⟨x, y⟩ := q
      simp only [Items.PairEq, Prod.mk.injEq] at hends h
      omega
    split
    · -- bridge
      rename_i hL
      have hL' : o.cls.lowval d = d + 1 := by simpa using hL
      obtain ⟨t, hsub, hvs', hdep⟩ := hE.bd_bridge hT hL'
      have hts : s.tstack = t :: base := by rw [hE.tstack, hsub]; rfl
      have ht : t ∈ s.tstack := by rw [hts]; exact List.mem_cons_self ..
      have hside : t.spans.1 = [] := by
        have := hE.bd_side hT hge; rw [if_pos hL', hts] at this; exact this t rfl
      have hvert : t.spans.2 = [vertItem o.dest] := by
        have := hc.bd_vert hT hge; rw [if_pos hL', hts] at this; exact this t rfl
      have hhead : s.tstack.head! = t := by rw [hts]; rfl
      set vsI := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) with hvsI
      set s₂ := { s₁ with items := s₁.items.push ⟨.I, (none, none), []⟩ } with hs₂
      have b₂ : BStep D s s₂ := b₁.trans (BStep.alloc b₁.inv b₁.shape .I)
      have h₂ : s₂.CloseInv := h₁.alloc' b₁.shape .I
      set s₃ := { s₂ with items := s₂.items.modify s₁.items.size fun it => { it with vs := vsI } }
        with hs₃
      have b₃ : BStep D s s₃ := b₂.trans (BStep.modifyVs_leaf b₂.inv b₂.shape s₁.items.size vsI
        (Items.ch_push_size _ rfl))
      have hszne : s₁.items.size ≠ edgeItem s.g o.e := by rw [hsz₁]; exact Nat.ne_of_gt hqsz
      have hszne' : s₁.items.size ≠ vertItem curV := by rw [hsz₁]; exact Nat.ne_of_gt hvsz
      have h₃ : s₃.CloseInv := by
        apply h₂.modifyLoose (x := s₁.items.size)
          (by show s.g.nv < s₁.items.size; rw [hsz₁]; have := hs.size; omega) _
          (fun it => { it with vs := vsI }) (fun _ => rfl)
        apply cnt_eq_zero_of_free
        · intro t' ht' hm
          exact Nat.lt_irrefl _ (hsz₁ ▸ hs.span t' ht' _ hm)
        · intro p h
          have h' := (isParent_push_iff _ _ _ _).1 h
          rw [Items.IsParent, hch₁] at h'
          exact Nat.lt_irrefl _ (hsz₁ ▸ hs.ch_lt p _ h')
      set s₄ := { s₃ with tstack := s₃.tstack.tail } with hs₄
      have b₄ : BStep D s s₄ := b₃.trans
        (BStep.pop' b₃.inv b₃.shape (s₀ := s) t base hts rfl rfl
          (fun a i => (Items.Below_modify_vs _ _ _ _).trans
            ((Items.Below_push_nil _ rfl).trans (hB₁ a i)))
          (hok.gone hT t base hts))
      have h₄ : s₄.CloseInv := h₃.pop
      have hB₄ : ∀ a i, Items.Below s₄.items a i ↔ Items.Below s.items a i := fun a i =>
        (Items.Below_modify_vs _ _ _ _).trans
          ((Items.Below_push_nil _ rfl).trans (hB₁ a i))
      have hn₄ : Items.NoParent s₄.items (edgeItem s.g o.e) := fun p h =>
        hq_root₁ p ((isParent_push_iff _ _ _ _).1 ((isParent_modifyVs_iff _ _ _ _ _).1 h))
      have hvroot₄ : ∀ p, ¬ Items.IsParent s₄.items p (vertItem curV) := fun p h =>
        hv_root₁ p ((isParent_push_iff _ _ _ _).1 ((isParent_modifyVs_iff _ _ _ _ _).1 h))
      have hchI : Items.ch s₄.items s₁.items.size = [] := by
        show Items.ch ((s₁.items.push _).modify _ _) _ = []
        rw [Items.ch_modify_vs]; exact Items.ch_push_size _ rfl
      have htyI : Items.type s₄.items s₁.items.size = .I := by
        show Items.type ((s₁.items.push _).modify _ _) _ = _
        rw [Items.type_modify_vs]; exact Items.type_push_size _
      have hvsI' : Items.vs s₄.items s₁.items.size = vsI := by
        show Items.vs ((s₁.items.push _).modify _ _) _ = _
        rw [Items.vs_modify_self _ _ _ (by simp)]
      have hty₄ : ∀ j, j ≠ s₁.items.size → Items.type s₄.items j = Items.type s.items j := by
        intro j hj
        show Items.type ((s₁.items.push _).modify _ _) _ = _
        rw [Items.type_modify_vs, Items.type_push_of_ne _ hj, hty₁]
      have hvs₄ : ∀ j, j ≠ s₁.items.size → Items.vs s₄.items j = Items.vs s₁.items j := by
        intro j hj
        show Items.vs ((s₁.items.push _).modify _ _) _ = _
        rw [Items.vs_modify_ne _ _ _ _ (Ne.symm hj), Items.vs_push_of_ne _ hj]
      have hch₄ : ∀ j, Items.ch s₄.items j = Items.ch s.items j := by
        intro j
        show Items.ch ((s₁.items.push _).modify _ _) _ = _
        rw [Items.ch_modify_vs, Items.ch_push_nil _ rfl, hch₁]
      have hsz₄ : s₄.items.size = s.items.size + 1 := by simp [hs₄, hs₃, hs₂, hsz₁]
      have hqsz₄ : edgeItem s.g o.e < s₄.items.size := by rw [hsz₄]; exact Nat.lt_succ_of_lt hqsz
      set cs := s₁.items.size :: s.tstack.head!.spans.2 with hcs₀
      have hcs : cs = s₁.items.size :: t.spans.2 := by rw [hcs₀, hhead]
      set items₅ := s₄.items.modify (edgeItem s.g o.e) fun it => { it with ch := cs } with hitems₅
      have hch₅ : Items.ch items₅ (edgeItem s.g o.e) = cs := by
        rw [hitems₅, Items.ch_modify_self _ _ _ hqsz₄]
      have hpar₅ : ∀ c, Items.IsParent items₅ (edgeItem s.g o.e) c ↔ c ∈ cs := fun c => by
        rw [Items.IsParent, hch₅]
      have hbelow₅ : ∀ c ∈ t.spans.2, ∀ j, Items.Below items₅ c j ↔ Items.Below s.items c j := by
        intro c hcm j
        rw [hitems₅, Items.Below_modify_of_not_below (edgeItem s.g o.e) (fun it => { it with ch := cs })
          (fun hb => hok.q_free t ht (List.mem_append_right _ ((hn₄.below_eq hb) ▸ hcm)))]
        exact hB₄ c j
      have hbelowI : ∀ j, Items.Below items₅ s₁.items.size j → j = s₁.items.size := by
        intro j hj
        rcases hj.head_cases with h | ⟨c, hc, _⟩
        · exact h.symm
        · exfalso
          rw [Items.IsParent, hitems₅, Items.ch_modify_ne _ _ _ _ (Ne.symm hszne), hchI] at hc
          exact List.not_mem_nil hc
      have hEB : ∀ e', e' < s.g.ne →
          (Items.EdgeBelow s.g items₅ (edgeItem s.g o.e) e' ↔ e' = o.e ∨ t.edges s.g s.items e') := by
        intro e' he'
        constructor
        · intro hb
          rcases hb.head_cases with h | ⟨c, hc, hb'⟩
          · left; have : 1 + s.g.nv + o.e = 1 + s.g.nv + e' := h; omega
          · right
            rw [hpar₅, hcs, List.mem_cons] at hc
            rcases hc with rfl | hc
            · exfalso
              have h2 : 1 + s.g.nv + e' = s₁.items.size := hbelowI _ hb'
              rw [hsz₁] at h2; have := hs.size; omega
            · exact ⟨c, List.mem_append_right _ hc, (hbelow₅ c hc _).1 hb'⟩
        · rintro (rfl | ⟨c, hc, hb⟩)
          · exact .refl
          · rw [hside, List.nil_append] at hc
            exact .head ((hpar₅ c).2 (by rw [hcs]; exact List.mem_cons_of_mem _ hc))
              ((hbelow₅ c hc _).2 hb)
      have hsubt : ∀ e', e' < s.g.ne → s.g.Inc e' o.dest → e' = o.e ∨ t.edges s.g s.items e' := by
        intro e' he' hinc
        by_cases heq : e' = o.e
        · exact Or.inl heq
        · obtain ⟨t', ht', hte⟩ := hE.sub_cover e' he' heq (hE.dest_edges hT e' he' hinc)
          rw [hsub, List.mem_singleton] at ht'
          exact Or.inr (ht' ▸ hte)
      have hatt : ∀ v, Items.Att s.g items₅ (edgeItem s.g o.e) v → v = curV := by
        rintro v ⟨e₁, e₂, he₁, he₂, hi₁, hi₂, hb₁, hb₂⟩
        rw [hEB e₁ he₁] at hb₁
        rw [hEB e₂ he₂] at hb₂
        have hnd : v ≠ o.dest := fun hv => hb₂ (hsubt e₂ he₂ (hv ▸ hi₂))
        rcases hb₁ with rfl | hb₁
        · rcases eq_of_inc_pairEq hends hi₁ with h | h
          · exact absurd h hnd
          · exact h
        · have hE' := (hi.entries [] t base (by rw [hts]; rfl)).attached v e₁ e₂ he₁ he₂ hb₁
            (fun h => hb₂ (Or.inr h)) hi₁ hi₂
          rw [TEntry.Term'_nil] at hE'
          rcases hE' with h | ⟨k, hk₁, hk₂, hk⟩
          · exact absurd (h.trans hvs') hnd
          · rw [hdep] at hk₁
            rw [hc.hD, if_pos hT] at hk₂
            obtain rfl : k = d + 1 := by omega
            exact absurd (hk.trans (hE.sv_child hT)) hnd
      have hnew : Items.CloseAt s.g items₅ (edgeItem s.g o.e) := by
        refine Items.CloseAt.rootQ (c := s₁.items.size) helt ?_ (by rw [hch₅, hcs]; simp) ?_
          (Graph.inc_of_pairEq hends).2 hatt ?_ ?_ (fun h => absurd h hedge) ?_
        · rw [hitems₅, Items.type_modify_ch _ _ _, hty₄ _ (Ne.symm hszne)]
          exact hs.edge o.e helt
        · rw [hitems₅, Items.vs_modify_ch, hvs₄ _ (Ne.symm hszne)]; exact hvs₁
        · rw [hitems₅, Items.type_modify_ch _ _ _, htyI]; simp
        · rw [hitems₅, Items.type_modify_ch _ _ _, htyI]; intro h; cases h
        · intro _
          refine ⟨o.dest, hc.dest_lt, pairEq_swap' hends, by rw [hch₅, hcs, hvert], ?_⟩
          rw [hitems₅, Items.vs_modify_ch, hvsI', hvsI, hsv]
          cases s.stackDir[d]!
          · exact ⟨_, _, rfl, Or.inl rfl⟩
          · exact ⟨_, _, rfl, Or.inr rfl⟩
      have hcs' : ∀ c ∈ cs, Items.CloseAt s.g s₄.items c := by
        intro c hcm
        rw [hcs, List.mem_cons] at hcm
        rcases hcm with rfl | hcm
        · exact Items.CloseAt.leafIO (by rw [hsz₁]; exact hs.size) (Or.inl htyI) hchI
        · have hpos : 0 < s₃.cnt c := by
            have : 0 < spansCount s₃.tstack c := spansCount_pos_of_mem_head!
              (by show c ∈ s.tstack.head!.spans.1 ++ _; rw [hhead]; exact List.mem_append_right _ hcm)
            dsimp [cnt]; omega
          exact h₃.closed c (cnt_lt_size b₃.shape hpos) (Or.inr hpos)
      have h₅ := h₄.writeChildren hqsz₄ hn₄ cs hnew hcs'
      have b₅ : BStep D s { s₄ with items := items₅ } := b₄.trans
        (BStep.modifyCh b₄.inv b₄.shape (edgeItem s.g o.e) (fun it => { it with ch := cs }) hqlt
          (fun _ => rfl) hn₄ (fun t' ht' => hok.q_free t' (List.mem_of_mem_tail ht'))
          (fun _ c hc => by
            show c < s₄.items.size
            rcases List.mem_cons.1 hc with rfl | hc
            · rw [hsz₄, hsz₁]; exact Nat.lt_succ_self _
            · rw [hhead] at hc
              exact Nat.lt_of_lt_of_le (hs.span t ht c (List.mem_append_right _ hc))
                (by rw [hsz₄]; omega)))
      have hz : ({ s₄ with items := items₅ } : WalkState).cnt (vertItem curV) = 0 := by
        apply cnt_eq_zero_of_free
        · intro t' ht'; exact hok.v_free t' (List.mem_of_mem_tail ht')
        · intro p h
          rcases Items.IsParent_modify h with h | ⟨rfl, _, hcm⟩
          · exact hvroot₄ p h
          · have hcm' : vertItem curV ∈ cs := hcm
            rw [hcs] at hcm'
            rcases List.mem_cons.1 hcm' with h | h
            · exact hszne' h.symm
            · exact hok.v_free t ht (List.mem_append_right _ h)
      have h₆ := h₅.vertex_append' b₅.shape curV o.e hvlt helt hz hnew
        (by rw [hch₅, hcs]; simp)
      exact h₆
    · -- completed block
      rename_i hL
      have hL' : o.cls.lowval d ≠ d + 1 := by simpa using hL
      obtain ⟨t₁, t₂, hsub, hv₁, hd₁, hv₂, hd₂⟩ := hE.bd_comp hT hge hL'
      have hts : s.tstack = t₁ :: t₂ :: base := by rw [hE.tstack, hsub]; rfl
      have ht₁ : t₁ ∈ s.tstack := by rw [hts]; exact List.mem_cons_self ..
      have ht₂ : t₂ ∈ s.tstack := by rw [hts]; exact List.mem_cons_of_mem _ (List.mem_cons_self ..)
      have hsides := hE.bd_side hT hge
      rw [if_neg hL', hts] at hsides
      have hside₁ : t₁.spans.2 = [] := hsides.1 t₁ rfl
      have hside₂ : t₂.spans.1 = [] := hsides.2 t₂ rfl
      have hvert : t₂.spans.2 = [vertItem o.dest] := by
        have := hc.bd_vert hT hge; rw [if_neg hL', hts] at this; exact this t₂ rfl
      obtain ⟨c, hc₁, hct, hcq, a, b', hcvs, hcp⟩ := by
        have := hc.bd_node hT hge hL'; rw [hts] at this; exact this t₁ rfl
      have hcne : c ≠ edgeItem s.g o.e := fun h =>
        hok.q_free t₁ ht₁ (List.mem_append_left _ (by rw [hc₁, ← h]; exact List.mem_singleton_self _))
      have hhead₁ : s.tstack.head! = t₁ := by rw [hts]; rfl
      have hhead₂ : s.tstack.tail.head! = t₂ := by rw [hts]; rfl
      set s₂ := { s₁ with tstack := s₁.tstack.tail } with hs₂
      set s₃ := { s₂ with tstack := s₂.tstack.tail } with hs₃
      have b₂ : BStep D s s₂ := b₁.trans
        (BStep.pop' b₁.inv b₁.shape (s₀ := s) t₁ (t₂ :: base) hts rfl rfl hB₁
          (hok.gone hT t₁ (t₂ :: base) hts))
      have b₃ : BStep D s s₃ := b₂.trans
        (BStep.pop' b₂.inv b₂.shape (s₀ := s) t₂ base
          (by show s.tstack.tail = _; rw [hts, List.tail_cons]) rfl rfl hB₁
          (hok.gone₂ hT hL' t₁ t₂ base hts))
      have h₃ : s₃.CloseInv := h₁.pop.pop
      have hn₃ : Items.NoParent s₃.items (edgeItem s.g o.e) := hq_root₁
      have hqsz₃ : edgeItem s.g o.e < s₃.items.size := by rw [hsz₁]; exact hqsz
      set cs := s.tstack.head!.spans.1 ++ s.tstack.tail.head!.spans.2 with hcs₀
      have hcs : cs = t₁.spans.1 ++ t₂.spans.2 := by rw [hcs₀, hhead₁, hhead₂]
      have hcs₂ : cs = [c, vertItem o.dest] := by rw [hcs, hc₁, hvert]; rfl
      set items₄ := s₃.items.modify (edgeItem s.g o.e) fun it => { it with ch := cs } with hitems₄
      have hch₄ : Items.ch items₄ (edgeItem s.g o.e) = cs := by
        rw [hitems₄, Items.ch_modify_self _ _ _ hqsz₃]
      have hpar₄ : ∀ c', Items.IsParent items₄ (edgeItem s.g o.e) c' ↔ c' ∈ cs := fun c' => by
        rw [Items.IsParent, hch₄]
      have hbelow₄ : ∀ c' ∈ cs, ∀ j, Items.Below items₄ c' j ↔ Items.Below s.items c' j := by
        intro c' hcm j
        rw [hitems₄, Items.Below_modify_of_not_below (edgeItem s.g o.e) (fun it => { it with ch := cs })
          (fun hb => ?_)]
        · exact hB₁ c' j
        · have := hn₃.below_eq hb
          subst this
          rw [hcs, List.mem_append] at hcm
          rcases hcm with hcm | hcm
          · exact hok.q_free t₁ ht₁ (List.mem_append_left _ hcm)
          · exact hok.q_free t₂ ht₂ (List.mem_append_right _ hcm)
      have hEB : ∀ e', e' < s.g.ne → (Items.EdgeBelow s.g items₄ (edgeItem s.g o.e) e' ↔
          e' = o.e ∨ t₁.edges s.g s.items e' ∨ t₂.edges s.g s.items e') := by
        intro e' he'
        constructor
        · intro hb
          rcases hb.head_cases with h | ⟨c', hc', hb'⟩
          · left; have : 1 + s.g.nv + o.e = 1 + s.g.nv + e' := h; omega
          · right
            have hb'' := (hbelow₄ c' ((hpar₄ c').1 hc') _).1 hb'
            rw [hpar₄, hcs, List.mem_append] at hc'
            rcases hc' with hc' | hc'
            · exact Or.inl ⟨c', List.mem_append_left _ hc', hb''⟩
            · exact Or.inr ⟨c', List.mem_append_right _ hc', hb''⟩
        · rintro (rfl | ⟨c', hc', hb⟩ | ⟨c', hc', hb⟩)
          · exact .refl
          · rw [hside₁, List.append_nil] at hc'
            have hm : c' ∈ cs := by rw [hcs]; exact List.mem_append_left _ hc'
            exact .head ((hpar₄ c').2 hm) ((hbelow₄ c' hm _).2 hb)
          · rw [hside₂, List.nil_append] at hc'
            have hm : c' ∈ cs := by rw [hcs]; exact List.mem_append_right _ hc'
            exact .head ((hpar₄ c').2 hm) ((hbelow₄ c' hm _).2 hb)
      have hsubt : ∀ e', e' < s.g.ne → s.g.Inc e' o.dest →
          e' = o.e ∨ t₁.edges s.g s.items e' ∨ t₂.edges s.g s.items e' := by
        intro e' he' hinc
        by_cases heq : e' = o.e
        · exact Or.inl heq
        · obtain ⟨t', ht', hte⟩ := hE.sub_cover e' he' heq (hE.dest_edges hT e' he' hinc)
          rw [hsub, List.mem_cons, List.mem_singleton] at ht'
          rcases ht' with rfl | rfl
          · exact Or.inr (Or.inl hte)
          · exact Or.inr (Or.inr hte)
      have hterm₁ : ∀ v, t₁.Term D s v → v ≠ o.dest → v = curV := by
        intro v h hnd
        rcases h with h | ⟨k, hk₁, hk₂, hk⟩
        · exact absurd (h.trans hv₁) hnd
        · rw [hd₁] at hk₁
          rw [hc.hD, if_pos hT] at hk₂
          rcases Nat.lt_or_ge k (d + 1) with hk' | hk'
          · obtain rfl : k = d := by omega
            exact hk.trans hsv
          · obtain rfl : k = d + 1 := by omega
            exact absurd (hk.trans (hE.sv_child hT)) hnd
      have hterm₂ : ∀ v, t₂.Term D s v → v ≠ o.dest → False := by
        intro v h hnd
        rcases h with h | ⟨k, hk₁, hk₂, hk⟩
        · exact hnd (h.trans hv₂)
        · rw [hd₂] at hk₁
          rw [hc.hD, if_pos hT] at hk₂
          obtain rfl : k = d + 1 := by omega
          exact hnd (hk.trans (hE.sv_child hT))
      have hatt : ∀ v, Items.Att s.g items₄ (edgeItem s.g o.e) v → v = curV := by
        rintro v ⟨e₁, e₂, he₁, he₂, hi₁, hi₂, hb₁, hb₂⟩
        rw [hEB e₁ he₁] at hb₁
        rw [hEB e₂ he₂] at hb₂
        have hnd : v ≠ o.dest := fun hv => hb₂ (hsubt e₂ he₂ (hv ▸ hi₂))
        rcases hb₁ with rfl | hb₁ | hb₁
        · rcases eq_of_inc_pairEq hends hi₁ with h | h
          · exact absurd h hnd
          · exact h
        · have hE' := (hi.entries [] t₁ (t₂ :: base) (by rw [hts]; rfl)).attached v e₁ e₂ he₁ he₂
            hb₁ (fun h => hb₂ (Or.inr (Or.inl h))) hi₁ hi₂
          rw [TEntry.Term'_nil] at hE'
          exact hterm₁ v hE' hnd
        · have hE' := (hi.entries [t₁] t₂ base (by rw [hts]; rfl)).attached v e₁ e₂ he₁ he₂
            hb₁ (fun h => hb₂ (Or.inr (Or.inr h))) hi₁ hi₂
          rcases hE' with h | ⟨t', ht', h⟩
          · exact (hterm₂ v h hnd).elim
          · rw [List.mem_singleton] at ht'
            subst ht'
            exact hterm₁ v h hnd
      have hnew : Items.CloseAt s.g items₄ (edgeItem s.g o.e) := by
        refine Items.CloseAt.rootQ (c := c) helt ?_ (by rw [hch₄, hcs₂]; simp) ?_
          (Graph.inc_of_pairEq hends).2 hatt ?_ ?_ (fun h => absurd h hedge) ?_
        · rw [hitems₄, Items.type_modify_ch _ _ _, hty₁]; exact hs.edge o.e helt
        · rw [hitems₄, Items.vs_modify_ch]; exact hvs₁
        · rw [hitems₄, Items.type_modify_ch _ _ _, hty₁]; exact hct
        · rw [hitems₄, Items.type_modify_ch _ _ _, hty₁,
            Items.ch_modify_ne _ _ _ _ (Ne.symm hcne), hch₁]
          exact hcq
        · intro _
          refine ⟨o.dest, hc.dest_lt, pairEq_swap' hends, by rw [hch₄, hcs₂], a, b', ?_, hcp⟩
          rw [hitems₄, Items.vs_modify_ch, hvs₁' c hcne]; exact hcvs
      have hcs' : ∀ c' ∈ cs, Items.CloseAt s.g s₃.items c' := by
        intro c' hcm
        have hpos : 0 < s₁.cnt c' := by
          rw [hcs, List.mem_append] at hcm
          rcases hcm with hcm | hcm
          · have : 0 < spansCount s₁.tstack c' := spansCount_pos_of_mem_head!
              (by show c' ∈ s.tstack.head!.spans.1 ++ _; rw [hhead₁]; exact List.mem_append_left _ hcm)
            dsimp [cnt]; omega
          · have : 0 < spansCount s₁.tstack c' := spansCount_pos_of_mem_tail_head!
              (by show c' ∈ s.tstack.tail.head!.spans.1 ++ _; rw [hhead₂]
                  exact List.mem_append_right _ hcm)
            dsimp [cnt]; omega
        exact h₁.closed c' (cnt_lt_size b₁.shape hpos) (Or.inr hpos)
      have h₄ := h₃.writeChildren hqsz₃ hn₃ cs hnew hcs'
      have b₄ : BStep D s { s₃ with items := items₄ } := b₃.trans
        (BStep.modifyCh b₃.inv b₃.shape (edgeItem s.g o.e) (fun it => { it with ch := cs }) hqlt
          (fun _ => rfl) hn₃
          (fun t' ht' => hok.q_free t' (List.mem_of_mem_tail (List.mem_of_mem_tail ht')))
          (fun _ c' hc' => by
            show c' < s₁.items.size
            rw [hsz₁]
            have hc'' : c' ∈ cs := hc'
            rw [hcs] at hc''
            rcases List.mem_append.1 hc'' with hc' | hc'
            · exact hs.span t₁ ht₁ c' (List.mem_append_left _ hc')
            · exact hs.span t₂ ht₂ c' (List.mem_append_right _ hc')))
      have hz : ({ s₃ with items := items₄ } : WalkState).cnt (vertItem curV) = 0 := by
        apply cnt_eq_zero_of_free
        · intro t' ht'
          exact hok.v_free t' (List.mem_of_mem_tail (List.mem_of_mem_tail ht'))
        · intro p h
          rcases Items.IsParent_modify h with h | ⟨rfl, _, hcm⟩
          · exact hv_root₁ p h
          · have hcm' : vertItem curV ∈ cs := hcm
            rw [hcs] at hcm'
            rcases List.mem_append.1 hcm' with h | h
            · exact hok.v_free t₁ ht₁ (List.mem_append_left _ h)
            · exact hok.v_free t₂ ht₂ (List.mem_append_right _ h)
      have h₅ := h₄.vertex_append' b₄.shape curV o.e hvlt helt hz hnew
        (by rw [hch₄, hcs₂]; simp)
      exact h₅
  · -- self-loop
    rename_i hT
    have hT' : o.cls.isTree = false := by simpa using hT
    have hends := hc.ends_back hT'
    have hloop := hc.bd_loop hT' hge
    have hedge : (s.g.edges[o.e]!).1 = (s.g.edges[o.e]!).2 := by
      rw [hloop] at hends
      generalize s.g.edges[o.e]! = q at hends
      obtain ⟨x, y⟩ := q
      simp only [Items.PairEq, Prod.mk.injEq] at hends
      omega
    set s₂ := { s₁ with totSelfLoops := s₁.totSelfLoops + 1 } with hs₂
    have b₂ : BStep D s s₂ := b₁.trans (BStep.frame b₁.inv b₁.shape)
    have h₂ : s₂.CloseInv := h₁.frame rfl rfl fun _ h => h
    set s₃ := { s₂ with items := s₂.items.push ⟨.O, (none, none), []⟩ } with hs₃
    have b₃ : BStep D s s₃ := b₂.trans (BStep.alloc b₂.inv b₂.shape .O)
    have h₃ : s₃.CloseInv := h₂.alloc' b₂.shape .O
    set s₄ := { s₃ with items := s₃.items.modify s₂.items.size fun it => { it with vs := (some curV, none) } }
      with hs₄
    have hsz₂ : s₂.items.size = s.items.size := hsz₁
    have b₄ : BStep D s s₄ := b₃.trans (BStep.modifyVs_leaf b₃.inv b₃.shape s₂.items.size
      (some curV, none) (Items.ch_push_size _ rfl))
    have hszne : s₂.items.size ≠ edgeItem s.g o.e := by rw [hsz₂]; exact Nat.ne_of_gt hqsz
    have hszne' : s₂.items.size ≠ vertItem curV := by rw [hsz₂]; exact Nat.ne_of_gt hvsz
    have h₄ : s₄.CloseInv := by
      apply h₃.modifyLoose (x := s₂.items.size)
        (by show s.g.nv < s₂.items.size; rw [hsz₂]; have := hs.size; omega) _
        (fun it => { it with vs := (some curV, none) }) (fun _ => rfl)
      apply cnt_eq_zero_of_free
      · intro t' ht' hm
        exact Nat.lt_irrefl _ (hsz₂ ▸ hs.span t' ht' _ hm)
      · intro p h
        have h' := (isParent_push_iff _ _ _ _).1 h
        rw [Items.IsParent, hch₁] at h'
        exact Nat.lt_irrefl _ (hsz₂ ▸ hs.ch_lt p _ h')
    have hn₄ : Items.NoParent s₄.items (edgeItem s.g o.e) := fun p h =>
      hq_root₁ p ((isParent_push_iff _ _ _ _).1 ((isParent_modifyVs_iff _ _ _ _ _).1 h))
    have hvroot₄ : ∀ p, ¬ Items.IsParent s₄.items p (vertItem curV) := fun p h =>
      hv_root₁ p ((isParent_push_iff _ _ _ _).1 ((isParent_modifyVs_iff _ _ _ _ _).1 h))
    have hchO : Items.ch s₄.items s₂.items.size = [] := by
      show Items.ch ((s₁.items.push _).modify _ _) _ = []
      rw [Items.ch_modify_vs]; exact Items.ch_push_size _ rfl
    have htyO : Items.type s₄.items s₂.items.size = .O := by
      show Items.type ((s₁.items.push _).modify _ _) _ = _
      rw [Items.type_modify_vs]; exact Items.type_push_size _
    have hvsO : Items.vs s₄.items s₂.items.size = (some curV, none) := by
      show Items.vs ((s₁.items.push _).modify _ _) _ = _
      rw [Items.vs_modify_self _ _ _ (by show s₁.items.size < (s₁.items.push _).size; simp)]
    have hty₄ : ∀ j, j ≠ s₂.items.size → Items.type s₄.items j = Items.type s.items j := by
      intro j hj
      show Items.type ((s₁.items.push _).modify _ _) _ = _
      rw [Items.type_modify_vs, Items.type_push_of_ne _ hj, hty₁]
    have hvs₄ : ∀ j, j ≠ s₂.items.size → Items.vs s₄.items j = Items.vs s₁.items j := by
      intro j hj
      show Items.vs ((s₁.items.push _).modify _ _) _ = _
      rw [Items.vs_modify_ne _ _ _ _ (Ne.symm hj), Items.vs_push_of_ne _ hj]
    have hsz₄ : s₄.items.size = s.items.size + 1 := by simp [hs₄, hs₃, hs₂, hsz₁]
    have hqsz₄ : edgeItem s.g o.e < s₄.items.size := by rw [hsz₄]; exact Nat.lt_succ_of_lt hqsz
    set cs := [s₂.items.size] with hcs
    set items₅ := s₄.items.modify (edgeItem s.g o.e) fun it => { it with ch := cs } with hitems₅
    have hch₅ : Items.ch items₅ (edgeItem s.g o.e) = cs := by
      rw [hitems₅, Items.ch_modify_self _ _ _ hqsz₄]
    have hEB : ∀ e', e' < s.g.ne →
        Items.EdgeBelow s.g items₅ (edgeItem s.g o.e) e' → e' = o.e := by
      intro e' he' hb
      rcases hb.head_cases with h | ⟨c, hc, hb'⟩
      · have : 1 + s.g.nv + o.e = 1 + s.g.nv + e' := h; omega
      · exfalso
        rw [Items.IsParent, hch₅, hcs, List.mem_singleton] at hc
        subst hc
        rcases hb'.head_cases with h | ⟨c', hc', _⟩
        · have h2 : 1 + s.g.nv + e' = s₂.items.size := h.symm
          rw [hsz₂] at h2; have := hs.size; omega
        · rw [Items.IsParent, hitems₅, Items.ch_modify_ne _ _ _ _ (Ne.symm hszne), hchO] at hc'
          exact List.not_mem_nil hc'
    have hatt : ∀ v, Items.Att s.g items₅ (edgeItem s.g o.e) v → v = curV := by
      rintro v ⟨e₁, e₂, he₁, he₂, hi₁, hi₂, hb₁, _⟩
      obtain rfl := hEB e₁ he₁ hb₁
      rcases eq_of_inc_pairEq hends hi₁ with h | h
      · exact h
      · exact h.trans hloop
    have hnew : Items.CloseAt s.g items₅ (edgeItem s.g o.e) := by
      refine Items.CloseAt.rootQ (c := s₂.items.size) helt ?_ (by rw [hch₅, hcs]; simp) ?_
        (Graph.inc_of_pairEq hends).1 hatt ?_ ?_ ?_ (fun h => absurd hedge h)
      · rw [hitems₅, Items.type_modify_ch _ _ _, hty₄ _ (Ne.symm hszne)]
        exact hs.edge o.e helt
      · rw [hitems₅, Items.vs_modify_ch, hvs₄ _ (Ne.symm hszne)]; exact hvs₁
      · rw [hitems₅, Items.type_modify_ch _ _ _, htyO]; simp
      · rw [hitems₅, Items.type_modify_ch _ _ _, htyO]; intro h; cases h
      · intro _
        refine ⟨by rw [hch₅, hcs], ?_⟩
        rw [hitems₅, Items.vs_modify_ch]; exact hvsO
    have hcs' : ∀ c ∈ cs, Items.CloseAt s.g s₄.items c := by
      intro c hcm
      rw [hcs, List.mem_singleton] at hcm
      subst hcm
      exact Items.CloseAt.leafIO (by rw [hsz₂]; exact hs.size) (Or.inr htyO) hchO
    have h₅ := h₄.writeChildren hqsz₄ hn₄ cs hnew hcs'
    have b₅ : BStep D s { s₄ with items := items₅ } := b₄.trans
      (BStep.modifyCh b₄.inv b₄.shape (edgeItem s.g o.e) (fun it => { it with ch := cs }) hqlt
        (fun _ => rfl) hn₄ (fun t' ht' => hok.q_free t' ht')
        (fun _ c hc => by
          show c < s₄.items.size
          rw [hsz₄, List.mem_singleton.1 hc, hsz₂]; exact Nat.lt_succ_self _))
    have hz : ({ s₄ with items := items₅ } : WalkState).cnt (vertItem curV) = 0 := by
      apply cnt_eq_zero_of_free
      · intro t' ht'; exact hok.v_free t' ht'
      · intro p h
        rcases Items.IsParent_modify h with h | ⟨rfl, _, hcm⟩
        · exact hvroot₄ p h
        · exact hszne' (List.mem_singleton.1 hcm).symm
    have h₆ := h₅.vertex_append' b₅.shape curV o.e hvlt helt hz hnew
      (by rw [hch₅, hcs]; simp)
    exact h₆

end Spqr.WalkState
