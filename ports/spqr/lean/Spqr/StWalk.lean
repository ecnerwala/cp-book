import Spqr.WalkSpec
import Spqr.StSpec

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

structure WalkState.StEntry (s : WalkState) (ord : Nat → Nat) (t : TEntry) : Prop
    extends s.StSides ord t where
  oriented : ∀ x ∈ Items.vertsOf s.items (t.spans.1 ++ t.spans.2), x ≠ t.top s →
    (s.stackDir[t.topDepth]! = false → ord (t.top s) < ord x) ∧
    (s.stackDir[t.topDepth]! = true → ord x < ord (t.top s))

/-- The st-invariant of a walk state at depth `d`: every entry attached at or above the current
depth has sorted, distinct sides; the top entry is in addition oriented towards its upper
terminal; and whenever `finishEdge` closes the top entry it is st-numbered. -/
structure WalkState.StInv (s : WalkState) (d : Nat) (ord : Nat → Nat) : Prop where
  sides : ∀ t ∈ s.tstack, t.topDepth ≤ d → s.StSides ord t
  top : ∀ t, s.tstack.head? = some t → t.topDepth ≤ d → s.StEntry ord t
  items : ∀ i, i < s.items.size →
    Items.type s.items i = .S ∨ Items.type s.items i = .P ∨ Items.type s.items i = .R →
    Items.StItem s.items i

/-- Admitted: along the first-child chain of an ear every vertex has the ear's lowval, so
`stackDir` is constant along it and every piece attached along the ear lands on the same side
(the semantic half of `ear_uniform_side`; the data half is `merge_onSide`). -/
theorem chain_stackDir_const (s : WalkState) (l d : Nat) (hl : l < d)
    (hchain : ∀ d', l < d' → d' ≤ d → s.stackDir[d']! = s.stackDir[l + 1]!) :
    ∀ t ∈ s.tstack, l < t.topDepth → t.topDepth ≤ d →
      (∀ d', l < d' → d' ≤ d → s.stackVerts[d']! ≠ t.vStart) →
      t.OnSide s.stackDir[l + 1]! := by
  sorry

/-- Admitted: the top entry is closable whenever `finishEdge` closes it
(the `finishTstackTop` calls in the type-1, type-2 and block-boundary branches). -/
theorem finishEdge_topClosable (s : WalkState) (d : Nat) (ord : Nat → Nat) (hinv : s.StInv d ord) :
    s.TopClosable := by
  sorry

/-- Admitted: `finishEdge` preserves the st-invariant. -/
theorem finishEdge_stInv (s : WalkState) (curV d : Nat) (o : DfsOut) (origTstack : Nat)
    (hasVert : Bool) (ord : Nat → Nat) (hinv : s.StInv d ord) :
    ((finishEdge curV d o origTstack hasVert).run s).2.StInv d ord := by
  sorry

/-- Admitted: closing the top entry into `item` yields an st-numbered item. -/
theorem finishTstackTop_stItem (s : WalkState) (item : ItemId) (h : s.TopClosable)
    (hitem : item < s.items.size) :
    Items.StItem ((WalkM.finishTstackTop item).run s).2.items item := by
  sorry

end Spqr
