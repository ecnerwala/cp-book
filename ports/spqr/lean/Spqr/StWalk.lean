import Spqr.WalkSpec
import Spqr.StSpec
import Spqr.RelabelSt
import Spqr.WalkWF
import Spqr.StRef
import Spqr.StRefEt
import Spqr.StRestrict

/-!
# Walk-level st-ordering invariant

The walk builds the Even–Tarjan st-numbering of every block by ear insertion: each tstack entry
holds a list of the eventual st-ordering with a hole in the middle for the part of the DFS path
that is still open.  `spans.1` is the finished part to the left of the hole, `spans.2` the part to
the right, and `stackDir`/`edgeDir` select the side a new piece is attached to (`setSides`).

Everything here is additive to `Walk.lean`/`WalkSpec.lean`.  The data-level facts about
`setSides`, `pushTstack`, `mergeTstackTops` and the child fold are proved; the preservation of the
invariant through `finishEdge` and its consequence `walk_st` are admitted.
-/

namespace Spqr

/-- Vertices of the V items of a list of items, in order. -/
def Items.vertsOf (items : Items) (l : List ItemId) : List Nat :=
  (l.filter fun c => items.type c = .V).map (· - 1)

namespace TEntry

/-- The items of the entry, read around the deeper part `inner`. -/
def wrap (t : TEntry) (inner : List ItemId) : List ItemId :=
  t.spans.1 ++ inner ++ t.spans.2

/-- One of the two sides is empty: every piece of the entry was attached on the same side. -/
def OneSided (t : TEntry) : Prop := t.spans.1 = [] ∨ t.spans.2 = []

/-- All items of the entry are on the side selected by `dir`. -/
def OnSide (t : TEntry) (dir : Bool) : Prop := getSide t.spans (!dir) = []

theorem OnSide.oneSided {t : TEntry} {dir : Bool} (h : t.OnSide dir) : t.OneSided := by
  unfold OnSide getSide at h
  cases dir <;> simp_all [OneSided]

end TEntry

/-- The items of a tstack segment (top first) in nesting order: each entry's `spans.1` to the left
of everything below it and its `spans.2` to the right. -/
def nest : List TEntry → List ItemId
  | [] => []
  | t :: rest => t.wrap (nest rest)

theorem setSides_onSide (dir : Bool) (l : List ItemId) (v td fi : Nat) :
    (⟨v, td, fi, setSides dir l []⟩ : TEntry).OnSide dir := by
  unfold TEntry.OnSide getSide setSides
  cases dir <;> rfl

theorem getSide_setSides (dir : Bool) (a b : α) : getSide (setSides dir a b) dir = a := by
  unfold getSide setSides; cases dir <;> rfl

theorem getSide_setSides_other (dir : Bool) (a b : α) : getSide (setSides dir a b) (!dir) = b := by
  unfold getSide setSides; cases dir <;> rfl

/-- `pushTstack` creates an entry whose single item sits on the side of `stackDir[topDepth]`. -/
theorem pushTstack_onSide (s : WalkState) (vStart topDepth item : Nat) :
    ∀ t, ((WalkM.pushTstack vStart topDepth item).run s).2.tstack.head? = some t →
      t.OnSide s.stackDir[topDepth]! := by
  intro t ht
  have : ((WalkM.pushTstack vStart topDepth item).run s).2.tstack.head? =
      some ⟨vStart, topDepth, s.nxtEdgeIdx, setSides s.stackDir[topDepth]! [item] []⟩ := rfl
  rw [this, Option.some.injEq] at ht
  subst ht
  exact setSides_onSide _ _ _ _ _

/-- Merging two entries whose items are on the same side keeps them on that side
(`ear_uniform_side`, data part: pieces attached along one ear share a side and stay one-sided). -/
theorem merge_onSide (a b : TEntry) (dir : Bool) (ha : a.OnSide dir) (hb : b.OnSide dir) :
    (⟨a.vStart, min a.topDepth b.topDepth, a.firstIdx,
      (b.spans.1 ++ a.spans.1, a.spans.2 ++ b.spans.2)⟩ : TEntry).OnSide dir := by
  unfold TEntry.OnSide getSide at *
  cases dir <;> simp_all

/-- The child fold when leaving a type-2 child puts the whole entry on the side of `!edgeDir`
(`type2_entry_one_sided`, data part). -/
theorem fold_onSide (t : TEntry) (edgeDir : Bool) (curV : Nat) :
    ({ t with vStart := curV, spans := setSides (!edgeDir) (t.spans.1 ++ t.spans.2) [] } : TEntry).OnSide
      (!edgeDir) :=
  setSides_onSide _ _ _ _ _

/-- The vertex list of an entry read as a closed item: terminals from `makeVs`, interior from the
side selected by `stackDir[topDepth]`. -/
def WalkState.entryVertList (s : WalkState) (t : TEntry) : List Nat :=
  let dir := s.stackDir[t.topDepth]!
  let vs := setSides dir (some (t.top s)) (some t.vStart)
  vs.1.toList ++ Items.vertsOf s.items (getSide t.spans dir) ++ vs.2.toList

/-- Skeleton edges of an entry read as a closed item: the `vs` of its non-V items plus the
terminal pair. -/
def WalkState.entryEdges (s : WalkState) (t : TEntry) : List (Nat × Nat) :=
  let dir := s.stackDir[t.topDepth]!
  let vs := setSides dir (some (t.top s)) (some t.vStart)
  (vs.1.getD 0, vs.2.getD 0) ::
    ((getSide t.spans dir).filter fun c => Items.type s.items c ≠ .V).map fun c =>
      ((Items.vs s.items c).1.getD 0, (Items.vs s.items c).2.getD 0)

/-- The top entry can be closed into an st-numbered item. -/
def WalkState.TopClosable (s : WalkState) : Prop :=
  ∀ t, s.tstack.head? = some t →
    Items.StList (s.entryVertList t) (s.entryEdges t) ∧
    ∀ p ∈ (s.entryEdges t).tail,
      (s.entryVertList t).idxOf p.1 < (s.entryVertList t).idxOf p.2

/-- Order facts about the two sides of a live entry relative to a fixed numbering `ord` of the
block (the eventual st-numbering): the V items are distinct and increasing, and lie on the
`stackDir[topDepth]` side of the entry's upper terminal. -/
structure WalkState.StSides (s : WalkState) (ord : Nat → Nat) (t : TEntry) : Prop where
  nodup : (Items.vertsOf s.items (t.spans.1 ++ t.spans.2)).Nodup
  sorted : (Items.vertsOf s.items (t.spans.1 ++ t.spans.2)).Pairwise fun x y => ord x < ord y

/-- The top entry is in addition oriented: read along `stackDir[topDepth]`, its vertices lie
strictly between the upper terminal `top` and the lower terminal `vStart`.  A vertex entry
(`vStart = stackVerts[topDepth]`, its own V item as the only item) is exempt where the terminals
coincide with the vertex. -/
structure WalkState.StEntry (s : WalkState) (ord : Nat → Nat) (t : TEntry) : Prop
    extends s.StSides ord t where
  oriented : ∀ x ∈ Items.vertsOf s.items (t.spans.1 ++ t.spans.2), x ≠ t.top s →
    (s.stackDir[t.topDepth]! = false → ord (t.top s) < ord x) ∧
    (s.stackDir[t.topDepth]! = true → ord x < ord (t.top s))
  bottom : ∀ x ∈ Items.vertsOf s.items (t.spans.1 ++ t.spans.2), x ≠ t.vStart →
    (s.stackDir[t.topDepth]! = false → ord x < ord t.vStart) ∧
    (s.stackDir[t.topDepth]! = true → ord t.vStart < ord x)
  ends : t.vStart ≠ t.top s →
    (s.stackDir[t.topDepth]! = false → ord (t.top s) < ord t.vStart) ∧
    (s.stackDir[t.topDepth]! = true → ord t.vStart < ord (t.top s))

/-- The vertices of an entry: its V items, both sides. -/
def WalkState.entryVerts (s : WalkState) (t : TEntry) : List Nat :=
  Items.vertsOf s.items (t.spans.1 ++ t.spans.2)

/-- The skeleton edges carried by the non-V items of `l` (their `vs`). -/
def WalkState.itemEdges (s : WalkState) (l : List ItemId) : List (Nat × Nat) :=
  (l.filter fun c => Items.type s.items c ≠ .V).map fun c =>
    ((Items.vs s.items c).1.getD 0, (Items.vs s.items c).2.getD 0)

/-- A vertex an entry's edges may reach at depth `d`: one of its own vertices, a terminal, or an
open-path vertex in its hole `(topDepth, d]`. -/
def WalkState.EntryReach (s : WalkState) (d : Nat) (t : TEntry) (x : Nat) : Prop :=
  x ∈ s.entryVerts t ∨ x = t.top s ∨ x = t.vStart ∨
    ∃ k, t.topDepth < k ∧ k ≤ d ∧ s.stackVerts[k]! = x

/-- Hole closure of a live entry relative to the numbering `ord`: its skeleton edges are oriented
along `ord` and stay within reach, and every vertex of the entry other than its lower terminal
has a lower and a higher neighbour along them (the hole, i.e. the open path through
`(topDepth, d]`, may supply the endpoint).  The lower terminal is exempt because a vertex entry
`[V v]` with `vStart = v` carries no edge until the vertex is finished. -/
structure WalkState.StHole (s : WalkState) (ord : Nat → Nat) (d : Nat) (t : TEntry) : Prop where
  edges : ∀ p ∈ s.itemEdges (t.spans.1 ++ t.spans.2),
    ord p.1 < ord p.2 ∧ s.EntryReach d t p.1 ∧ s.EntryReach d t p.2
  lower : ∀ x ∈ s.entryVerts t, x ≠ t.vStart → ∃ p ∈ s.itemEdges (t.spans.1 ++ t.spans.2),
    (p.1 = x ∧ ord p.2 < ord x) ∨ (p.2 = x ∧ ord p.1 < ord x)
  upper : ∀ x ∈ s.entryVerts t, x ≠ t.vStart → ∃ p ∈ s.itemEdges (t.spans.1 ++ t.spans.2),
    (p.1 = x ∧ ord x < ord p.2) ∨ (p.2 = x ∧ ord x < ord p.1)

/-- Close-site condition on the hole of `t` (the st-side of `WalkSpec.FinishTopOk.mid`): every
open-path vertex in `(topDepth, d]` is the lower terminal, one of the entry's own vertices, or
not an endpoint of any of its edges. -/
def WalkState.HoleClosed (s : WalkState) (d : Nat) (t : TEntry) : Prop :=
  ∀ k, t.topDepth < k → k ≤ d →
    s.stackVerts[k]! = t.vStart ∨ s.stackVerts[k]! ∈ s.entryVerts t ∨
    ∀ p ∈ s.itemEdges (t.spans.1 ++ t.spans.2), p.1 ≠ s.stackVerts[k]! ∧ p.2 ≠ s.stackVerts[k]!

/-- The st-invariant of a walk state at depth `d`: every entry attached at or above the current
depth has sorted, distinct sides; the top entry is in addition oriented towards its upper
terminal; and whenever `finishEdge` closes the top entry it is st-numbered. -/
structure WalkState.StInv (s : WalkState) (d : Nat) (ord : Nat → Nat) : Prop where
  sides : ∀ t ∈ s.tstack, t.topDepth ≤ d → s.StSides ord t
  top : ∀ t, s.tstack.head? = some t → t.topDepth ≤ d → s.StEntry ord t
  hole : ∀ t ∈ s.tstack, t.topDepth ≤ d → s.StHole ord d t
  items : ∀ i, i < s.items.size →
    Items.type s.items i = .S ∨ Items.type s.items i = .P ∨ Items.type s.items i = .R →
    Items.StItem s.items i
  /-- Every entry finished by an open-path vertex strictly below its upper terminal (`vStart`
  is `stackVerts[j]` for some `topDepth < j ≤ d`: back edges, type-1/type-2 folds and P
  pieces of the vertices on the path) lies on the side `stackDir[topDepth]`. Entries whose
  `vStart` has been popped (the eager vertex merge after a type-2 first edge) and entries at
  `topDepth = d` from earlier out-edges of the current vertex (`stackDir[d]` has been reset
  since) can be two-sided, see `PROOF.md` §7.4. -/
  onSide : ∀ t ∈ s.tstack, t.topDepth < d →
    (∃ j, t.topDepth < j ∧ j ≤ d ∧ s.stackVerts[j]! = t.vStart) →
    t.OnSide s.stackDir[t.topDepth]!

/-- Along the first-child chain of an ear every vertex has the ear's lowval, so `stackDir` is
constant along it (`hchain`, from `first_ret_lowval`, `chain_stackDir_step` and
`walkTree_stackDir_below` in `Spqr.StFrame`), and an entry finished on the chain
(`StInv.onSide`) is on the side `stackDir[l + 1]` of the whole chain.
The earlier statement with the hypothesis `∀ d' ∈ (l, d], stackVerts[d'] ≠ t.vStart` and
`t.topDepth ≤ d` is false: `gen.py` seed 6, depth `d = 14`, `l = 13`, has the entry
`vStart = 10, topDepth = 14` (back edge `10 → 2` pushed while `stackDir[14] = 1`, with
`stackDir[14] = 0` after the next out-edge of `2`) on `spans.2`, see `PROOF.md` §7.4. -/
theorem chain_stackDir_const (s : WalkState) (l d : Nat) (ord : Nat → Nat) (hinv : s.StInv d ord)
    (hchain : ∀ d', l < d' → d' ≤ d → s.stackDir[d']! = s.stackDir[l + 1]!) :
    ∀ t ∈ s.tstack, l < t.topDepth → t.topDepth < d →
      (∃ j, t.topDepth < j ∧ j ≤ d ∧ s.stackVerts[j]! = t.vStart) →
      t.OnSide s.stackDir[l + 1]! := by
  intro t ht hl hd hj
  have h := hinv.onSide t ht hd hj
  rwa [hchain t.topDepth hl (Nat.le_of_lt hd)] at h

theorem idxOf_lt_idxOf_iff (ord : Nat → Nat) (xs : List Nat)
    (hs : xs.Pairwise fun x y => ord x < ord y) {a b : Nat} (ha : a ∈ xs) (hb : b ∈ xs) :
    xs.idxOf a < xs.idxOf b ↔ ord a < ord b := by
  have ha' := List.idxOf_lt_length_iff.2 ha
  have hb' := List.idxOf_lt_length_iff.2 hb
  have key : ∀ {x y : Nat}, x ∈ xs → y ∈ xs → xs.idxOf x < xs.idxOf y → ord x < ord y := by
    intro x y hx hy hlt
    have hx' := List.idxOf_lt_length_iff.2 hx
    have hy' := List.idxOf_lt_length_iff.2 hy
    have := List.pairwise_iff_getElem.1 hs _ _ hx' hy' hlt
    rwa [List.getElem_idxOf, List.getElem_idxOf] at this
  refine ⟨key ha hb, fun h => ?_⟩
  rcases Nat.lt_trichotomy (xs.idxOf a) (xs.idxOf b) with hlt | heq | hgt
  · exact hlt
  · have : a = b := by
      rw [← List.getElem_idxOf ha', ← List.getElem_idxOf hb']
      simp only [heq]
    subst this; exact absurd h (Nat.lt_irrefl _)
  · exact absurd (key hb ha hgt) (Nat.not_lt.2 (Nat.le_of_lt h))

/-- An `ord`-increasing vertex list `a :: verts ++ [b]` whose edges are oriented along `ord`,
stay inside the list and give every interior vertex a lower and a higher neighbour is an
`StList`, with every edge oriented along it. -/
theorem stList_of_sorted (ord : Nat → Nat) (a b : Nat) (verts : List Nat) (es : List (Nat × Nat))
    (hab : ord a < ord b) (ha : ∀ x ∈ verts, ord a < ord x) (hb : ∀ x ∈ verts, ord x < ord b)
    (hsorted : verts.Pairwise fun x y => ord x < ord y)
    (hes : ∀ p ∈ es, ord p.1 < ord p.2 ∧ p.1 ∈ a :: (verts ++ [b]) ∧ p.2 ∈ a :: (verts ++ [b]))
    (hlow : ∀ x ∈ verts, ∃ p ∈ es, (p.1 = x ∧ ord p.2 < ord x) ∨ (p.2 = x ∧ ord p.1 < ord x))
    (hup : ∀ x ∈ verts, ∃ p ∈ es, (p.1 = x ∧ ord x < ord p.2) ∨ (p.2 = x ∧ ord x < ord p.1)) :
    Items.StList (a :: (verts ++ [b])) ((a, b) :: es) ∧
    ∀ p ∈ es, (a :: (verts ++ [b])).idxOf p.1 < (a :: (verts ++ [b])).idxOf p.2 := by
  set xs := a :: (verts ++ [b]) with hxs
  have hpw : xs.Pairwise fun x y => ord x < ord y := by
    rw [hxs, List.pairwise_cons, List.pairwise_append]
    refine ⟨fun x hx => ?_, hsorted, by simp, fun x hx y hy => ?_⟩
    · simp only [List.mem_append, List.mem_singleton] at hx
      rcases hx with hx | rfl
      · exact ha x hx
      · exact hab
    · simp only [List.mem_singleton] at hy; subst hy; exact hb x hx
  have hmem : ∀ x, x ∈ xs ↔ x = a ∨ x ∈ verts ∨ x = b := by
    intro x; simp [hxs]
  have hidx : ∀ {x y : Nat}, x ∈ xs → y ∈ xs → (xs.idxOf x < xs.idxOf y ↔ ord x < ord y) :=
    fun hx hy => idxOf_lt_idxOf_iff ord xs hpw hx hy
  have hirr : Std.Irrefl fun x y : Nat => ord x < ord y := ⟨fun x => Nat.lt_irrefl (ord x)⟩
  have hnodup : xs.Nodup := hpw.nodup
  have hesx : ∀ p ∈ (a, b) :: es, ord p.1 < ord p.2 ∧ p.1 ∈ xs ∧ p.2 ∈ xs := by
    intro p hp
    simp only [List.mem_cons] at hp
    rcases hp with rfl | hp
    · exact ⟨hab, by simp [hxs], by simp [hxs]⟩
    · exact hes p hp
  refine ⟨⟨hnodup, fun p hp => ?_, fun x hx hhd hlast => ?_⟩, fun p hp => ?_⟩
  · obtain ⟨h1, h2, h3⟩ := hesx p hp
    exact ⟨h2, h3, fun h => by rw [h] at h1; exact Nat.lt_irrefl _ h1⟩
  · have hxv : x ∈ verts := by
      rcases (hmem x).1 hx with rfl | hxv | rfl
      · exact absurd rfl hhd
      · exact hxv
      · exact absurd (by rw [hxs, ← List.cons_append, List.getLast?_concat]) hlast
    refine ⟨?_, ?_⟩
    · obtain ⟨p, hp, hcase⟩ := hlow x hxv
      obtain ⟨-, h1, h2⟩ := hesx p (List.mem_cons_of_mem _ hp)
      refine ⟨p, List.mem_cons_of_mem _ hp, ?_⟩
      rcases hcase with ⟨rfl, hlt⟩ | ⟨rfl, hlt⟩
      · exact Or.inl ⟨rfl, (hidx h2 h1).2 hlt⟩
      · exact Or.inr ⟨rfl, (hidx h1 h2).2 hlt⟩
    · obtain ⟨p, hp, hcase⟩ := hup x hxv
      obtain ⟨-, h1, h2⟩ := hesx p (List.mem_cons_of_mem _ hp)
      refine ⟨p, List.mem_cons_of_mem _ hp, ?_⟩
      rcases hcase with ⟨rfl, hlt⟩ | ⟨rfl, hlt⟩
      · exact Or.inl ⟨rfl, (hidx h1 h2).2 hlt⟩
      · exact Or.inr ⟨rfl, (hidx h2 h1).2 hlt⟩
  · obtain ⟨h1, h2, h3⟩ := hes p hp
    exact (hidx h2 h3).2 h1

theorem getSide_eq_append_of_onSide {t : TEntry} {dir : Bool} (h : t.OnSide dir) :
    getSide t.spans dir = t.spans.1 ++ t.spans.2 := by
  unfold TEntry.OnSide getSide at *
  cases dir <;> simp_all

theorem entryVertList_eq (s : WalkState) (t : TEntry) :
    s.entryVertList t =
      (setSides s.stackDir[t.topDepth]! (some (t.top s)) (some t.vStart)).1.toList ++
        Items.vertsOf s.items (getSide t.spans s.stackDir[t.topDepth]!) ++
        (setSides s.stackDir[t.topDepth]! (some (t.top s)) (some t.vStart)).2.toList := rfl

theorem entryEdges_eq (s : WalkState) (t : TEntry) :
    s.entryEdges t =
      ((setSides s.stackDir[t.topDepth]! (some (t.top s)) (some t.vStart)).1.getD 0,
        (setSides s.stackDir[t.topDepth]! (some (t.top s)) (some t.vStart)).2.getD 0) ::
        s.itemEdges (getSide t.spans s.stackDir[t.topDepth]!) := rfl

/-- Closability of a one-sided live entry from the invariant: `StEntry` gives the `ord`-increasing
vertex list between the terminals, `StHole` the neighbours, and the close-site condition
`HoleClosed` keeps the edges inside the list; the terminals are distinct and not among the
entry's vertices.  (The earlier statement `StInv d ord → TopClosable` is not provable: the top
entry's edges may reach into its hole before the hole has been closed, see `PROOF.md` §7.4.) -/
theorem stInv_topClosable (s : WalkState) (d : Nat) (ord : Nat → Nat) (hinv : s.StInv d ord)
    (t : TEntry) (ht : s.tstack.head? = some t) (htd : t.topDepth ≤ d)
    (hside : t.OnSide s.stackDir[t.topDepth]!) (hclosed : s.HoleClosed d t)
    (hne : t.vStart ≠ t.top s) (hvv : t.vStart ∉ s.entryVerts t) (htop : t.top s ∉ s.entryVerts t) :
    Items.StList (s.entryVertList t) (s.entryEdges t) ∧
    ∀ p ∈ (s.entryEdges t).tail, (s.entryVertList t).idxOf p.1 < (s.entryVertList t).idxOf p.2 := by
  have htm : t ∈ s.tstack := List.mem_of_mem_head? ht
  have hent := hinv.top t ht htd
  have hhole := hinv.hole t htm htd
  have hsideq := getSide_eq_append_of_onSide hside
  have hreach : ∀ x, s.EntryReach d t x →
      (∃ p ∈ s.itemEdges (t.spans.1 ++ t.spans.2), p.1 = x ∨ p.2 = x) →
      x ∈ s.entryVerts t ∨ x = t.top s ∨ x = t.vStart := by
    intro x hr ⟨p, hp, hpx⟩
    rcases hr with h | h | h | ⟨k, hk1, hk2, hk⟩
    · exact Or.inl h
    · exact Or.inr (Or.inl h)
    · exact Or.inr (Or.inr h)
    · rcases hclosed k hk1 hk2 with h | h | h
      · exact Or.inr (Or.inr (hk ▸ h))
      · exact Or.inl (hk ▸ h)
      · exact absurd hpx (by have := h p hp; rw [hk] at this; tauto)
  have hmem : ∀ a b x, x ∈ s.entryVerts t ∨ x = a ∨ x = b → x ∈ a :: (s.entryVerts t ++ [b]) := by
    intro a b x h; simp only [List.mem_cons, List.mem_append]; tauto
  have hmem' : ∀ a b x, x ∈ s.entryVerts t ∨ x = b ∨ x = a → x ∈ a :: (s.entryVerts t ++ [b]) := by
    intro a b x h; simp only [List.mem_cons, List.mem_append]; tauto
  have hes : ∀ (a b : Nat),
      (∀ x, x ∈ s.entryVerts t ∨ x = t.top s ∨ x = t.vStart → x ∈ a :: (s.entryVerts t ++ [b])) →
      ∀ p ∈ s.itemEdges (t.spans.1 ++ t.spans.2),
        ord p.1 < ord p.2 ∧ p.1 ∈ a :: (s.entryVerts t ++ [b]) ∧ p.2 ∈ a :: (s.entryVerts t ++ [b]) := by
    intro a b hm p hp
    obtain ⟨h1, h2, h3⟩ := hhole.edges p hp
    exact ⟨h1, hm _ (hreach p.1 h2 ⟨p, hp, Or.inl rfl⟩), hm _ (hreach p.2 h3 ⟨p, hp, Or.inr rfl⟩)⟩
  have hsorted : (s.entryVerts t).Pairwise fun x y => ord x < ord y := hent.sorted
  have hor : ∀ x ∈ s.entryVerts t, _ := fun x hx => hent.oriented x hx (ne_of_mem_of_not_mem hx htop)
  have hbot : ∀ x ∈ s.entryVerts t, _ := fun x hx => hent.bottom x hx (ne_of_mem_of_not_mem hx hvv)
  have hends := hent.ends hne
  have hlow : ∀ x ∈ s.entryVerts t, _ := fun x hx => hhole.lower x hx (ne_of_mem_of_not_mem hx hvv)
  have hup : ∀ x ∈ s.entryVerts t, _ := fun x hx => hhole.upper x hx (ne_of_mem_of_not_mem hx hvv)
  rw [entryVertList_eq, entryEdges_eq, hsideq]
  generalize hdir : s.stackDir[t.topDepth]! = dir at hor hbot hends ⊢
  cases dir
  · simp only [setSides, Bool.false_eq_true, ↓reduceIte, Option.toList,
      List.cons_append, Option.getD_some, List.tail_cons]
    exact stList_of_sorted ord (t.top s) t.vStart (s.entryVerts t) _ (hends.1 rfl)
      (fun x hx => (hor x hx).1 rfl) (fun x hx => (hbot x hx).1 rfl) hsorted
      (hes _ _ (hmem _ _)) hlow hup
  · simp only [setSides, ↓reduceIte, Option.toList,
      List.cons_append, Option.getD_some, List.tail_cons]
    exact stList_of_sorted ord t.vStart (t.top s) (s.entryVerts t) _ (hends.2 rfl)
      (fun x hx => (hbot x hx).2 rfl) (fun x hx => (hor x hx).2 rfl) hsorted
      (hes _ _ (hmem' _ _)) hlow hup

theorem finishTstackTop_items (s : WalkState) (item : ItemId) (t : TEntry) (rest : List TEntry)
    (hts : s.tstack = t :: rest) :
    ((WalkM.finishTstackTop item).run s).2.items =
      s.items.modify item fun it => { it with
        vs := setSides s.stackDir[t.topDepth]! (some (t.top s)) (some t.vStart),
        ch := getSide t.spans s.stackDir[t.topDepth]! } := by
  rcases s with ⟨g, tern, items, sv, sd, nei, fo, ts, tb, tsl⟩
  simp only at hts; subst hts
  rfl

/-- Closing the closable top entry `t` into `item` (not itself one of the entry's items) yields an
st-numbered item: `finishTstackTop` writes `entryVertList`/`entryEdges` into `item`. -/
theorem finishTstackTop_stItem (s : WalkState) (item : ItemId) (h : s.TopClosable)
    (t : TEntry) (ht : s.tstack.head? = some t) (hitem : item < s.items.size)
    (hnot : item ∉ getSide t.spans s.stackDir[t.topDepth]!) :
    Items.StItem ((WalkM.finishTstackTop item).run s).2.items item := by
  obtain ⟨rest, hts⟩ : ∃ rest, s.tstack = t :: rest := by
    cases hs : s.tstack with
    | nil => simp [hs] at ht
    | cons a r => simp [hs] at ht; exact ⟨r, by rw [ht]⟩
  obtain ⟨hlist, hor⟩ := h t ht
  have heq := finishTstackTop_items s item t rest hts
  generalize hdir : s.stackDir[t.topDepth]! = dir at heq hnot hlist hor
  generalize hvs : setSides dir (some (t.top s)) (some t.vStart) = vs at heq
  generalize hitems' : ((WalkM.finishTstackTop item).run s).2.items = items' at heq
  have hvs' : Items.vs items' item = vs := by
    simp [heq, Items.vs, Array.getElem?_modify, Array.getElem?_eq_getElem hitem]
  have hch' : Items.ch items' item = getSide t.spans dir := by
    simp [heq, Items.ch, Array.getElem?_modify, Array.getElem?_eq_getElem hitem]
  have htype' : ∀ c, Items.type items' c = Items.type s.items c := by
    intro c
    by_cases hc : item = c
    · subst hc; simp [heq, Items.type, Array.getElem?_modify, Array.getElem?_eq_getElem hitem]
    · simp [heq, Items.type, Array.getElem?_modify, hc]
  have hvsc : ∀ c, c ≠ item → Items.vs items' c = Items.vs s.items c := by
    intro c hc
    simp [heq, Items.vs, Array.getElem?_modify, Ne.symm hc]
  have hfilt : ∀ l : List ItemId,
      (l.filter fun c => Items.type items' c = .V) = l.filter fun c => Items.type s.items c = .V :=
    fun l => List.filter_congr fun c _ => by simp [htype']
  have hfilt' : ∀ l : List ItemId,
      (l.filter fun c => Items.type items' c ≠ .V) = l.filter fun c => Items.type s.items c ≠ .V :=
    fun l => List.filter_congr fun c _ => by simp [htype']
  have hvert : Items.vertList items' item = s.entryVertList t := by
    simp only [Items.vertList, WalkState.entryVertList, hvs', hch', Items.vertsOf, hdir, hvs, hfilt]
  have hedges : Items.virtualEdges items' item =
      ((getSide t.spans dir).filter fun c => Items.type s.items c ≠ .V).map fun c =>
        ((Items.vs s.items c).1.getD 0, (Items.vs s.items c).2.getD 0) := by
    simp only [Items.virtualEdges, hch', hfilt']
    exact List.map_congr_left fun c hc => by
      rw [hvsc c fun hci => hnot (hci ▸ (List.mem_filter.1 hc).1)]
  have hvs_shape : ∃ a b, vs = (some a, some b) := by
    rw [← hvs]; simp only [setSides]; cases dir <;> exact ⟨_, _, rfl⟩
  obtain ⟨a, b, hab⟩ := hvs_shape
  have hent : s.entryEdges t = (a, b) :: Items.virtualEdges items' item := by
    simp only [WalkState.entryEdges, hedges, hdir, hvs, hab]; rfl
  refine ⟨a, b, by rw [hvs', hab], ?_, ?_⟩
  · rw [hvert, ← hent]; exact hlist
  · intro p hp
    rw [hvert]
    exact hor p (by rw [hent, List.tail_cons]; exact hp)

/-- Close site: under `StInv` the top entry, once one-sided with its hole closed and its terminals
distinct and not among its vertices, closes into an st-numbered item (`stInv_topClosable` fed to
`finishTstackTop_stItem`). -/
theorem stInv_finishTstackTop_stItem (s : WalkState) (d : Nat) (ord : Nat → Nat)
    (hinv : s.StInv d ord) (item : ItemId) (t : TEntry) (ht : s.tstack.head? = some t)
    (htd : t.topDepth ≤ d) (hside : t.OnSide s.stackDir[t.topDepth]!) (hclosed : s.HoleClosed d t)
    (hne : t.vStart ≠ t.top s) (hvv : t.vStart ∉ s.entryVerts t) (htop : t.top s ∉ s.entryVerts t)
    (hitem : item < s.items.size) (hnot : item ∉ getSide t.spans s.stackDir[t.topDepth]!) :
    Items.StItem ((WalkM.finishTstackTop item).run s).2.items item :=
  finishTstackTop_stItem s item
    (fun t' ht' => by
      rw [ht] at ht'; cases ht'
      exact stInv_topClosable s d ord hinv t ht htd hside hclosed hne hvv htop)
    t ht hitem hnot

end Spqr
