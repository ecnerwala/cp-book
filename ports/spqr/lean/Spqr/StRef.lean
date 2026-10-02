import Spqr.ItemSpec
import Spqr.Frame
import Spqr.StSpec

/-!
# The st-order reference

An Even–Tarjan style description of the order in which the walk (`Spqr.Walk`) lists the children
of its S / P / R items. The DFS tree is walked in the same lowval-sorted order; every edge is
handled at its deeper endpoint, exactly where `walkOut` / `finishEdge` see it, and is placed on
one side of the open DFS path:

* an item pushed at depth `l` sits on side `dirs[l]`, where `dirs[l]` is the direction the vertex
  at depth `l` chose for its current out-edge (`walkOut`'s `setStackDir`: `false` for a block
  boundary edge, otherwise the opposite of the direction at the edge's lowval);
* the pieces of one block nest: later pieces lie outside earlier ones (`stNest`);
* when a returning tree edge is finished at a vertex that already has its vertex item
  (`finishEdge`'s vertex close), the whole sub-ear is folded onto side `dirs[lowval]` as one piece.

`refOrder` is the resulting st-order of the leaf items (vertices and edges) of each block; the
walk's child list of an S / P / R item is the restriction of `refOrder` to the item's subtree.
-/

namespace Spqr
open WalkM

/-- One push onto the ear stack, abstracted: `items` (already in st-order) on side `side`. -/
structure StPiece where
  side : Bool
  items : List ItemId

/-- The st-order of a stack of pieces, bottom first: later pieces lie outside earlier ones. -/
def stNestL : List StPiece → List ItemId
  | [] => []
  | p :: ps => stNestL ps ++ (if p.side then [] else p.items)
def stNestR : List StPiece → List ItemId
  | [] => []
  | p :: ps => (if p.side then p.items else []) ++ stNestR ps
def stNest (ps : List StPiece) : List ItemId := stNestL ps ++ stNestR ps

/-- A completed block: `root = some (v, w)` is the block-boundary tree edge it hangs from (`v` the
cut vertex, `w` its DFS child; the block's node is oriented `(v, w)`), `none` for the one-vertex
block of a DFS root; `items` is its st-order (vertex and edge items). -/
structure StBlock where
  root : Option (Nat × Nat)
  items : List ItemId

mutual
/-- The pieces pushed while walking the subtree `t` at depth `d`, with `dirs` the directions
chosen along the path above `t`; also the st-orders of the blocks completed inside `t`. -/
def refTree (g : Graph) (t : DfsTree) (d : Nat) (dirs : List Bool) :
    List StPiece × List StBlock :=
  match t with
  | .node v outs =>
    let (ps, blocks, hasVert) := refOuts g v d dirs outs false
    (if hasVert then ps else ps ++ [StPiece.mk true [vertItem v]], blocks)

def refOuts (g : Graph) (v d : Nat) (dirs : List Bool) :
    List DfsOut → Bool → List StPiece × List StBlock × Bool
  | [], hasVert => ([], [], hasVert)
  | o :: rest, hasVert =>
    let (ps, blocks, hasVert) := refOut g v d dirs o hasVert
    let (ps', blocks', hasVert') := refOuts g v d dirs rest hasVert
    (ps ++ ps', blocks ++ blocks', hasVert')

/-- The out-edge `o` of `v`, handled at `v`. -/
def refOut (g : Graph) (v d : Nat) (dirs : List Bool) (o : DfsOut) (hasVert : Bool) :
    List StPiece × List StBlock × Bool :=
  let l := o.cls.lowval d
  if d ≤ l then
    match o with
    | .tree _ _ child =>
      let (ps, blocks) := refTree g child (d + 1) (dirs ++ [false])
      ([], blocks ++ [⟨some (v, o.dest), stNest ps⟩], hasVert)
    | .back .. => ([], [], hasVert)
  else
    let lowDir := dirs.getD l false
    let sd := !lowDir
    let (pre, hasVert) :=
      if !hasVert && o.cls.isType1 then ([StPiece.mk sd [vertItem v]], true) else ([], hasVert)
    let (mid, blocks) :=
      match o with
      | .tree e _ child =>
        let (ps, blocks) := refTree g child (d + 1) (dirs ++ [sd])
        let sub := ps ++ [StPiece.mk sd [edgeItem g e]]
        (if hasVert then [StPiece.mk lowDir (stNest sub)] else sub, blocks)
      | .back e _ _ => ([StPiece.mk lowDir [edgeItem g e]], [])
    let (post, hasVert) := if hasVert then ([], true) else ([StPiece.mk sd [vertItem v]], true)
    (pre ++ mid ++ post, blocks, hasVert)
end

/-- The st-orders of all blocks of the forest, in completion order. -/
def refBlocks (g : Graph) (forest : List DfsTree) : List StBlock :=
  forest.flatMap fun t =>
    let (ps, blocks) := refTree g t 0 []
    blocks ++ [⟨none, stNest ps⟩]

/-- The reference st-order of the leaf items (blocks concatenated; no item lies in two blocks). -/
def refOrder (g : Graph) (forest : List DfsTree) : List ItemId :=
  (refBlocks g forest).flatMap StBlock.items

/-! ### Orientation: the block's vertex sequence and `vs` -/

/-- The vertices of a block in st-order: the cut vertex it hangs from, then the vertices of its
vertex items. -/
def StBlock.seq (g : Graph) (b : StBlock) : List Nat :=
  (b.root.map (·.1)).toList ++ b.items.filterMap fun x =>
    if 1 ≤ x ∧ x < 1 + g.nv then some (x - 1) else none

/-- The edges of a block: the boundary tree edge it hangs from and its edge items. -/
def StBlock.edges (g : Graph) (b : StBlock) : List (Nat × Nat) :=
  b.root.toList ++ b.items.filterMap fun x =>
    if 1 + g.nv ≤ x ∧ x < 1 + g.nv + g.ne then some g.edges[x - (1 + g.nv)]! else none

/-- Even–Tarjan on a block: its vertex sequence is an st-numbering of its edges. -/
def StBlock.St (g : Graph) (b : StBlock) : Prop := Items.StList (b.seq g) (b.edges g)

/-- `a` comes strictly before `b` in `xs`. -/
def Precedes (xs : List Nat) (a b : Nat) : Prop := a ∈ xs ∧ b ∈ xs ∧ xs.idxOf a < xs.idxOf b

instance (xs : List Nat) (a b : Nat) : Decidable (Precedes xs a b) := by unfold Precedes; infer_instance

/-- A `vs` pair oriented along `xs`. -/
def Oriented (xs : List Nat) : Option Nat × Option Nat → Prop
  | (some a, some b) => Precedes xs a b
  | _ => False

instance (xs : List Nat) (p : Option Nat × Option Nat) : Decidable (Oriented xs p) := by
  rcases p with ⟨_ | a, _ | b⟩ <;> simp only [Oriented] <;> infer_instance

instance (g : Graph) (b : StBlock) : Decidable (b.St g) := by unfold StBlock.St Items.StList; infer_instance

/-! ### Reading a tstack as pieces (the simulation relation for `walk_st'`) -/

/-- The st-order read off a tstack (top first): entries nest, higher entries outside lower ones
(`mergeTstackTops` puts the top entry's spans outside the next one's). -/
def readL : List TEntry → List ItemId
  | [] => []
  | t :: below => t.spans.1 ++ readL below
def readR : List TEntry → List ItemId
  | [] => []
  | t :: below => readR below ++ t.spans.2
def readStack (ts : List TEntry) : List ItemId := readL ts ++ readR ts

theorem stNestL_append (ps qs : List StPiece) : stNestL (ps ++ qs) = stNestL qs ++ stNestL ps := by
  induction ps with
  | nil => simp [stNestL]
  | cons p ps ih => simp [stNestL, ih]

theorem stNestR_append (ps qs : List StPiece) : stNestR (ps ++ qs) = stNestR ps ++ stNestR qs := by
  induction ps with
  | nil => simp [stNestR]
  | cons p ps ih => simp [stNestR, ih]

/-- Pushing the pieces `qs` on top of `ps` puts them outside. -/
theorem stNest_append (ps qs : List StPiece) :
    stNest (ps ++ qs) = stNestL qs ++ stNest ps ++ stNestR qs := by
  simp [stNest, stNestL_append, stNestR_append]

theorem stNest_single (side : Bool) (items : List ItemId) :
    stNest [⟨side, items⟩] = items := by
  cases side <;> simp [stNest, stNestL, stNestR]

/-- Reading the stack with one more entry on top: the entry is the piece `⟨dir, items⟩` if its
spans are `setSides dir items []`. -/
theorem readStack_cons_setSides (v d idx : Nat) (dir : Bool) (items : List ItemId)
    (ts : List TEntry) :
    readStack (⟨v, d, idx, setSides dir items []⟩ :: ts) =
      stNestL [⟨dir, items⟩] ++ readStack ts ++ stNestR [⟨dir, items⟩] := by
  cases dir <;> simp [readStack, readL, readR, setSides, stNestL, stNestR]

/-- `pushTstack` / `pushEdgeTstack` push the piece `⟨stackDir[d], [i]⟩`. -/
theorem readStack_pushTstack (s : WalkState) (v d : Nat) (i : ItemId) :
    readStack ((pushTstack v d i).run s).2.tstack =
      stNestL [⟨s.stackDir[d]!, [i]⟩] ++ readStack s.tstack ++ stNestR [⟨s.stackDir[d]!, [i]⟩] := by
  rw [run_pushTstack]
  exact readStack_cons_setSides _ _ _ _ _ _

/-- `mergeTstackTops` does not change the reading. -/
theorem readStack_mergeTops (b a : TEntry) (rest : List TEntry) :
    readStack (mergeTops (b :: a :: rest)) = readStack (b :: a :: rest) := by
  simp [mergeTops, readStack, readL, readR]

theorem readStack_mergeTstackTops (s : WalkState) (h : 2 ≤ s.tstack.length) :
    readStack (mergeTstackTops.run s).2.tstack = readStack s.tstack := by
  rw [run_mergeTstackTops]
  match s.tstack, h with
  | b :: a :: rest, _ => exact readStack_mergeTops b a rest

/-- The vertex-close fold (`closeVertTail`) turns the top entry into the single piece
`⟨dir, spans.1 ++ spans.2⟩`. -/
theorem readStack_fold (dir : Bool) (curV : Nat) (t : TEntry) (rest : List TEntry) :
    readStack ({ t with vStart := curV, spans := setSides dir (t.spans.1 ++ t.spans.2) [] } :: rest) =
      stNestL [⟨dir, t.spans.1 ++ t.spans.2⟩] ++ readStack rest ++
        stNestR [⟨dir, t.spans.1 ++ t.spans.2⟩] :=
  readStack_cons_setSides _ _ _ _ _ _

/-- Replace the item `i` by `ch` (the children it was closed with). -/
def expandItem (i : ItemId) (ch : List ItemId) (l : List ItemId) : List ItemId :=
  l.flatMap fun x => if x = i then ch else [x]

theorem expandItem_append (i : ItemId) (ch l₁ l₂ : List ItemId) :
    expandItem i ch (l₁ ++ l₂) = expandItem i ch l₁ ++ expandItem i ch l₂ := by
  simp [expandItem]

theorem expandItem_of_not_mem (i : ItemId) (ch l : List ItemId) (h : i ∉ l) :
    expandItem i ch l = l := by
  induction l with
  | nil => rfl
  | cons x l ih =>
    simp only [List.mem_cons, not_or] at h
    have := ih h.2
    simp only [expandItem, List.flatMap_cons, Ne.symm h.1, ite_false, List.singleton_append] at this ⊢
    rw [this]

/-- `finishTstackTop item` on a one-sided top entry: expanding `item` back to the entry's items
recovers the reading (the item is new, so it occurs nowhere else). -/
theorem readStack_finishTstackTop (s : WalkState) (item : ItemId) (t : TEntry) (rest : List TEntry)
    (ht : s.tstack = t :: rest)
    (hside : getSide t.spans (!s.stackDir[t.topDepth]!) = [])
    (hnew : item ∉ readStack rest) :
    expandItem item (getSide t.spans s.stackDir[t.topDepth]!)
        (readStack ((finishTstackTop item).run s).2.tstack) =
      readStack s.tstack := by
  rcases s with ⟨g, tern, items, sv, sd, nei, fo, ts, tb, tsl⟩
  simp only at ht; subst ht
  show expandItem _ _ (readStack ({ t with spans := setSides sd[t.topDepth]! [item] [] } :: rest)) = _
  rw [readStack_cons_setSides, expandItem_append, expandItem_append, expandItem_of_not_mem _ _ _ hnew]
  simp only at hside
  cases hd : sd[t.topDepth]! <;> simp [hd, getSide] at hside ⊢ <;>
    simp [readStack, readL, readR, expandItem, stNestL, stNestR, hside]

/-- `maybeUnwrapNxt`'s reopen: the entry under the top holding the single item `i` (a closed S / P
node) is reopened to `i`'s children; expanding `i` to `ch` recovers the reading (`i` occurs nowhere
else on the stack). -/
theorem readStack_reopen (a b : TEntry) (rest : List TEntry) (dir : Bool) (i : ItemId)
    (ch : List ItemId) (hb : b.spans = setSides dir [i] []) (ha : i ∉ a.spans.1 ++ a.spans.2)
    (hrest : i ∉ readStack rest) :
    readStack (a :: { b with spans := setSides dir ch [] } :: rest) =
      expandItem i ch (readStack (a :: b :: rest)) := by
  have h1 : i ∉ a.spans.1 := fun h => ha (List.mem_append_left _ h)
  have h2 : i ∉ a.spans.2 := fun h => ha (List.mem_append_right _ h)
  have h3 : i ∉ readL rest := fun h => hrest (List.mem_append_left _ h)
  have h4 : i ∉ readR rest := fun h => hrest (List.mem_append_right _ h)
  have h7 : ∀ l, expandItem i ch (i :: l) = ch ++ expandItem i ch l := by
    intro l; simp [expandItem]
  cases dir <;> simp [readStack, readL, readR, setSides, hb, expandItem_append,
    expandItem_of_not_mem _ _ _ h1, expandItem_of_not_mem _ _ _ h2,
    expandItem_of_not_mem _ _ _ h3, expandItem_of_not_mem _ _ _ h4, h7]

theorem readStack_modifyNxt_reopen (s : WalkState) (a b : TEntry) (rest : List TEntry)
    (dir : Bool) (i : ItemId) (ch : List ItemId) (ht : s.tstack = a :: b :: rest)
    (hb : b.spans = setSides dir [i] []) (ha : i ∉ a.spans.1 ++ a.spans.2)
    (hrest : i ∉ readStack rest) :
    readStack ((modifyNxt fun t => { t with spans := setSides dir ch [] }).run s).2.tstack =
      expandItem i ch (readStack s.tstack) := by
  rw [run_modifyNxt, ht]
  exact readStack_reopen a b rest dir i ch hb ha hrest

/-! ### The statement: the walk's child lists are restrictions of `refOrder` -/

/-- The leaf items (vertices and edges) under `i`, in `ch` order; `fuel ≥` the depth of `i`. -/
def Items.leaves (items : Items) : Nat → ItemId → List ItemId
  | 0, i => [i]
  | fuel + 1, i =>
    match Items.type items i with
    | .V | .Q => [i]
    | _ => (Items.ch items i).flatMap (Items.leaves items fuel)

/-- Collapse runs of equal elements. -/
def collapseRuns : List ItemId → List ItemId
  | [] => []
  | [x] => [x]
  | x :: y :: xs => if x = y then collapseRuns (y :: xs) else x :: collapseRuns (y :: xs)

/-- The restriction of `order` to the children of `i`: each leaf is replaced by the child of `i`
it lies under, and runs are collapsed. -/
def restrictCh (items : Items) (fuel : Nat) (order : List ItemId) (i : ItemId) : List ItemId :=
  collapseRuns (order.filterMap fun x =>
    (Items.ch items i).find? fun c => x ∈ Items.leaves items fuel c)

/-- The orientation facts the walk provides for every S / P / R item `i`: `i` lies in one block
`b` of `blocks` (all its leaves are items of `b`, and an edge of `b` — an edge item of `b` or an
edge on `b.root`'s endpoints — that is below `i` is a leaf of `i`, i.e. not in a block hanging
off a V child of `i`), the `vs` of `i` and of its non-V children are oriented along `b.seq`
(`makeVs` orients by the same `stackDir` bit as the splice side), the V children of `i` lie
strictly between its endpoints in `b.seq`, and the block edges below a non-V child `c` have their
endpoints between `c`'s endpoints (the items' vertex sets are nested sub-ears). Differentially
tested by `check_stref`. -/
def VsOriented (g : Graph) (items : Items) (blocks : List StBlock) : Prop :=
  ∀ i, i < items.size →
    Items.type items i = .S ∨ Items.type items i = .P ∨ Items.type items i = .R →
    ∃ b ∈ blocks, (∀ x ∈ Items.leaves items items.size i, x ∈ b.items) ∧
      (∀ e, e < g.ne → (edgeItem g e ∈ b.items ∨ ∃ r ∈ b.root, Items.PairEq g.edges[e]! r) →
        Items.Below items i (edgeItem g e) → edgeItem g e ∈ Items.leaves items items.size i) ∧
      Oriented (b.seq g) (Items.vs items i) ∧
      (∀ s t, Items.vs items i = (some s, some t) → ∀ c ∈ Items.ch items i, Items.type items c = .V →
        Precedes (b.seq g) s (c - 1) ∧ Precedes (b.seq g) (c - 1) t) ∧
      (∀ c ∈ Items.ch items i, Items.type items c ≠ .V → Oriented (b.seq g) (Items.vs items c)) ∧
      ∀ c ∈ Items.ch items i, Items.type items c ≠ .V → ∀ u v, Items.vs items c = (some u, some v) →
        ∀ e, e < g.ne → (edgeItem g e ∈ b.items ∨ ∃ r ∈ b.root, Items.PairEq g.edges[e]! r) →
        Items.Below items c (edgeItem g e) →
        ∀ y, (g.edges[e]!).1 = y ∨ (g.edges[e]!).2 = y →
          (u = y ∨ Precedes (b.seq g) u y) ∧ (y = v ∨ Precedes (b.seq g) y v)

/-- The walk lists the children of every S / P / R item in the reference st-order
(differentially tested by `check_stref`). -/
theorem walk_st' (g : Graph) (tern : Bool) (vo eo : List Nat) (i : ItemId)
    (hi : i < (g.walk tern (g.dfsForest vo eo)).items.size)
    (ht : Items.type (g.walk tern (g.dfsForest vo eo)).items i = .S ∨
      Items.type (g.walk tern (g.dfsForest vo eo)).items i = .P ∨
      Items.type (g.walk tern (g.dfsForest vo eo)).items i = .R) :
    Items.ch (g.walk tern (g.dfsForest vo eo)).items i =
      restrictCh (g.walk tern (g.dfsForest vo eo)).items
        (g.walk tern (g.dfsForest vo eo)).items.size
        (refOrder g (g.dfsForest vo eo)) i := by
  sorry

/-- The walk's `vs` are oriented along the reference blocks (proved in the same simulation as
`walk_st'`; differentially tested by `check_stref`). -/
theorem walk_vsOriented (g : Graph) (tern : Bool) (vo eo : List Nat) :
    VsOriented g (g.walk tern (g.dfsForest vo eo)).items (refBlocks g (g.dfsForest vo eo)) := by
  sorry


end Spqr
