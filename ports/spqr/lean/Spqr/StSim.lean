import Spqr.StRef

/-!
# The simulation relation for `walk_st'` / `walk_vsOriented` (PROOF.md §7.6)

The walk's tstack above the base of the current frame reads (`readStack`) as `stNest` of the
reference's pieces for that frame, up to expanding the S / P / R items closed inside the ear back to
their leaves (`Expands`); every S / P / R item is either live (below an item of the open stack) or
finished, i.e. a contiguous segment of a reference block with the `VsOriented` clauses already
satisfied; the blocks used are those of the *truncated* DFS tree (`truncTree`: along the open path
every vertex keeps its finished out-edges and the current tree edge), whose st-orders only grow
outwards as the walk proceeds. Checked at every `finishEdge` boundary by `check_stsim`.
-/

namespace Spqr

/-! ### Leaf expansion -/

mutual
/-- `Expands items x L`: the leaves of `x` (V / Q items are their own leaf). -/
inductive Expands (items : Items) : ItemId → List ItemId → Prop
  | leaf (x : ItemId) (h : Items.type items x = .V ∨ Items.type items x = .Q) : Expands items x [x]
  | node (x : ItemId) (L : List ItemId) (h : ¬ (Items.type items x = .V ∨ Items.type items x = .Q))
      (hL : ExpandsList items (Items.ch items x) L) : Expands items x L
/-- Concatenated expansion of a list of items. -/
inductive ExpandsList (items : Items) : List ItemId → List ItemId → Prop
  | nil : ExpandsList items [] []
  | cons {x : ItemId} {L xs Ls : List ItemId} (hx : Expands items x L) (hxs : ExpandsList items xs Ls) :
      ExpandsList items (x :: xs) (L ++ Ls)
end

/-- The reference's `dirs` are the walk's `stackDir` below `d`. -/
def DirsOf (s : WalkState) (d : Nat) : List Bool := (List.range d).map fun k => s.stackDir[k]!

/-- The stack segment `new` reads, up to expansion, as the pieces `ps`. -/
def StRead (items : Items) (new : List TEntry) (ps : List StPiece) : Prop :=
  ExpandsList items (readStack new) (stNest ps)

/-! ### Items against blocks -/

/-- The body of `VsOriented` for one item `i` and one block `b`. -/
def VsOrientedAt (g : Graph) (items : Items) (b : StBlock) (i : ItemId) : Prop :=
  (∀ x ∈ Items.leaves items items.size i, x ∈ b.items) ∧
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

theorem vsOriented_iff (g : Graph) (items : Items) (blocks : List StBlock) :
    VsOriented g items blocks ↔
      ∀ i, i < items.size → Items.type items i = .S ∨ Items.type items i = .P ∨ Items.type items i = .R →
        ∃ b ∈ blocks, VsOrientedAt g items b i := Iff.rfl

/-- `i` is finished in the block `b`: its leaves are a contiguous segment of `b.items` and the
`VsOriented` clauses hold. -/
structure InBlock (g : Graph) (items : Items) (b : StBlock) (i : ItemId) : Prop where
  seg : ∃ A L B, Expands items i L ∧ b.items = A ++ L ++ B
  oriented : VsOrientedAt g items b i

/-- The items relative to the stack and the blocks: stack items are roots and distinct; every
S / P / R item is live (below a stack item) or finished in one of `blocks`. -/
structure StItems (g : Graph) (s : WalkState) (blocks : List StBlock) : Prop where
  roots : ∀ x ∈ readStack s.tstack, ∀ p, ¬ Items.IsParent s.items p x
  nodup : (readStack s.tstack).Nodup
  closed : ∀ i, i < s.items.size →
    Items.type s.items i = .S ∨ Items.type s.items i = .P ∨ Items.type s.items i = .R →
    (∃ x ∈ readStack s.tstack, Items.Below s.items x i) ∨ ∃ b ∈ blocks, InBlock g s.items b i

/-! ### The truncated tree and the relation at the boundaries -/

/-- A frame of the open path: vertex `v`, its finished out-edges `done`, the current tree edge
`o` (whose child is the next frame). -/
structure PathFrame where
  v : Nat
  done : List DfsOut
  o : DfsOut

/-- The DFS tree truncated along the open path, with the subtree `t` at its bottom. -/
def truncTree : List PathFrame → DfsTree → DfsTree
  | [], t => t
  | f :: fs, t => .node f.v (f.done ++ [.tree f.o.e f.o.cls (truncTree fs t)])

/-- The relation at the end of `walkTree t d` (`d = fs.length`) under the path `fs`, with `prev`
the finished trees of the forest and `base` the tstack when the walk of `t` started. -/
structure StSim (g : Graph) (prev : List DfsTree) (fs : List PathFrame) (t : DfsTree)
    (base : List TEntry) (s : WalkState) : Prop where
  read : ∃ new, s.tstack = new ++ base ∧ StRead s.items new (refTree g t fs.length (DirsOf s fs.length)).1
  items : StItems g s (refBlocks g (prev ++ [truncTree fs t]))

/-- The relation at the start of an out-edge of `v` (`d = fs.length`): the out-edges `done` are
finished, `hasVert` is the walk's flag. -/
structure StSimOuts (g : Graph) (prev : List DfsTree) (fs : List PathFrame) (v : Nat)
    (done : List DfsOut) (hasVert : Bool) (base : List TEntry) (s : WalkState) : Prop where
  read : ∃ new, s.tstack = new ++ base ∧
    StRead s.items new (refOuts g v fs.length (DirsOf s fs.length) done false).1
  hasVert : (refOuts g v fs.length (DirsOf s fs.length) done false).2.2 = hasVert
  items : StItems g s (refBlocks g (prev ++ [truncTree fs (.node v done)]))

end Spqr
