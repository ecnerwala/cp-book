import Mathlib.Logic.Relation
import Spqr.Dfs

/-!
# Separation pairs relative to a lowpoint-sorted DFS tree

Pure graph theory over the data produced by phase 1 (PROOF.md §3). Paths are defined directly on
the multigraph (`Graph.Reach`), the DFS forest is abstracted to per-vertex sorted out-lists and
depths (`DfsData`), and the phase-1 facts another module proves are collected in `DfsData.Spec`
(`DfsForestSpec g forest` for the forest actually computed). The facts themselves are proven in
`Spqr.Proofs.SepPair`.
-/

namespace Spqr

namespace Graph

variable (g : Graph)

/-- Edge `e` joins `x` and `y` (in either orientation). -/
def Joins (e x y : Nat) : Prop := g.edges[e]? = some (x, y) ∨ g.edges[e]? = some (y, x)

/-- `x` is an endpoint of `e`. -/
def IsEnd (e x : Nat) : Prop := ∃ y, g.Joins e x y

def Adj (x y : Nat) : Prop := ∃ e, g.Joins e x y

/-- Walks in `g` all of whose vertices (endpoints included) satisfy `ok`. -/
inductive Reach (ok : Nat → Prop) : Nat → Nat → Prop
  | refl {x} : ok x → Reach ok x x
  | tail {x y z} : Reach ok x y → g.Adj y z → ok z → Reach ok x z

/-- Edges `e`, `e'` are joined by a walk whose vertices satisfy `ok` (so edges with no `ok`
endpoint are only joined to themselves). With `ok = (· ∉ {a, b})` these are the separation classes
of `{a, b}`. -/
def EdgeConn (ok : Nat → Prop) (e e' : Nat) : Prop :=
  e = e' ∨ ∃ x y, g.IsEnd e x ∧ g.IsEnd e' y ∧ g.Reach ok x y

/-- No vertex separates two edges: the edge set is a block. -/
def TwoConnected : Prop := ∀ v e e', e < g.ne → e' < g.ne → g.EdgeConn (· ≠ v) e e'

/-- Separation classes of `E` with `{a, b}` removed. -/
def SepClass (a b : Nat) : Nat → Nat → Prop := g.EdgeConn fun x => x ≠ a ∧ x ≠ b

/-- Standard SPQR convention: `{a, b}` is a separation pair iff the separation classes number at
least two with at least two edges each, or at least three. -/
def SeparationPair (a b : Nat) : Prop :=
  a ≠ b ∧
    ((∃ e₁ e₂ e₃, e₁ < g.ne ∧ e₂ < g.ne ∧ e₃ < g.ne ∧
        ¬g.SepClass a b e₁ e₂ ∧ ¬g.SepClass a b e₁ e₃ ∧ ¬g.SepClass a b e₂ e₃) ∨
      (∃ e₁ e₁' e₂ e₂', e₁ < g.ne ∧ e₁' < g.ne ∧ e₂ < g.ne ∧ e₂' < g.ne ∧
        e₁ ≠ e₁' ∧ g.SepClass a b e₁ e₁' ∧ e₂ ≠ e₂' ∧ g.SepClass a b e₂ e₂' ∧ ¬g.SepClass a b e₁ e₂))

end Graph

def DfsOut.isTree : DfsOut → Bool
  | .back .. => false
  | .tree .. => true

/-- A DFS forest seen abstractly: the sorted out-list and depth of every vertex, and the root of
the tree the facts below are about. -/
structure DfsData where
  root : Nat
  depth : Nat → Nat
  outs : Nat → List DfsOut

namespace DfsData

variable (d : DfsData)

/-- `p → c` is a tree edge. -/
def IsParent (p c : Nat) : Prop := ∃ o ∈ d.outs p, o.isTree = true ∧ o.dest = c

/-- `a` is an ancestor-or-self of `x`; `T_a = {x | d.Anc a x}`. -/
def Anc : Nat → Nat → Prop := Relation.ReflTransGen d.IsParent

/-- Some back edge out of (or inside) `T_c` lands at depth `l`. -/
def Returns (c l : Nat) : Prop :=
  ∃ u, ∃ o ∈ d.outs u, d.Anc c u ∧ o.isTree = false ∧ d.depth o.dest = l

/-- The phase-1 facts (PROOF.md §1, Lemmas 1.1–1.2 and "no cross edges") about `g`'s forest. -/
structure Spec (g : Graph) : Prop where
  /-- Every out-edge of `v` is an edge of `g` from `v` to its `dest`. -/
  joins : ∀ v, ∀ o ∈ d.outs v, g.Joins o.e v o.dest
  /-- Every edge appears as an out-edge ... -/
  edge_out : ∀ e, e < g.ne → ∃ v, ∃ o ∈ d.outs v, o.e = e
  /-- ... exactly once. -/
  out_inj : ∀ v v', ∀ o ∈ d.outs v, ∀ o' ∈ d.outs v', o.e = o'.e → v = v' ∧ o = o'
  nodup : ∀ v, (d.outs v).Nodup
  /-- A child is entered by one tree edge. -/
  tree_inj : ∀ v, ∀ o ∈ d.outs v, ∀ o' ∈ d.outs v,
    o.isTree = true → o'.isTree = true → o.dest = o'.dest → o = o'
  /-- Back edges go to an ancestor-or-self: no cross edges. -/
  back_anc : ∀ v, ∀ o ∈ d.outs v, o.isTree = false → d.Anc o.dest v
  parent_unique : ∀ p p' c, d.IsParent p c → d.IsParent p' c → p = p'
  depth_parent : ∀ p c, d.IsParent p c → d.depth c = d.depth p + 1
  depth_root : d.depth d.root = 0
  /-- Out-lists are sorted by `OutClass.rank`. -/
  sorted : ∀ v, (d.outs v).Pairwise fun o o' => o.cls.rank ≤ o'.cls.rank
  /-- `classify` agrees with the lowpoint table. -/
  cls_bridge : ∀ v, ∀ o ∈ d.outs v,
    o.cls = .bridge ↔ o.isTree = true ∧ ∀ l, l ≤ d.depth v → ¬d.Returns o.dest l
  cls_component : ∀ v, ∀ o ∈ d.outs v,
    o.cls = .component ↔
      o.isTree = true ∧ d.Returns o.dest (d.depth v) ∧ ∀ l, l < d.depth v → ¬d.Returns o.dest l
  cls_selfLoop : ∀ v, ∀ o ∈ d.outs v, o.cls = .selfLoop ↔ o.isTree = false ∧ o.dest = v
  cls_backEdge : ∀ v, ∀ o ∈ d.outs v, ∀ l,
    o.cls = .ret l .backEdge ↔ o.isTree = false ∧ o.dest ≠ v ∧ d.depth o.dest = l
  cls_type1 : ∀ v, ∀ o ∈ d.outs v, ∀ l,
    o.cls = .ret l .type1Child ↔
      o.isTree = true ∧ l < d.depth v ∧ d.Returns o.dest l ∧ (∀ l', l' < l → ¬d.Returns o.dest l') ∧
        ∀ l', l < l' → l' < d.depth v → ¬d.Returns o.dest l'
  cls_type2 : ∀ v, ∀ o ∈ d.outs v, ∀ l,
    o.cls = .ret l .type2Child ↔
      o.isTree = true ∧ l < d.depth v ∧ d.Returns o.dest l ∧ (∀ l', l' < l → ¬d.Returns o.dest l') ∧
        ∃ l', l < l' ∧ l' < d.depth v ∧ d.Returns o.dest l'

/-- Every endpoint of an edge of `g` lies in the tree of `root`: the facts of §3 are about one
block, whose DFS forest is a single tree. -/
def Rooted (g : Graph) : Prop := ∀ e x, g.IsEnd e x → d.Anc d.root x

/-- `e` has an endpoint in `T_c`. -/
def EndIn (c e : Nat) (g : Graph) : Prop := ∃ x, g.IsEnd e x ∧ d.Anc c x

/-- Edges of the *above* part of `E − {a, b}` for `a = anc b l` with child `a'` toward `b`: those
with an endpoint outside `T_{a'} ∪ {a}` (`b`'s own back edges above `a` included). -/
def Above (a' a b e : Nat) (g : Graph) : Prop :=
  ∃ x, g.IsEnd e x ∧ x ≠ a ∧ x ≠ b ∧ ¬d.Anc a' x

/-- Edges of the *between* part: those with an endpoint in `T_{a'} − T_b`. -/
def Between (a' b e : Nat) (g : Graph) : Prop := ∃ x, g.IsEnd e x ∧ d.Anc a' x ∧ ¬d.Anc b x

end DfsData

/-! ### Where an out-edge of `b` attaches, relative to `a = anc b l` -/

namespace OutClass

/-- Attaches to the *above* part: returns strictly above `l`. -/
def AttachesAbove (l : Nat) : OutClass → Prop
  | ret l' _ => l' < l
  | _ => False

/-- A type-1 edge to `l`: a type-1 child or a back edge with `lowval = l`. -/
def IsType1To (l : Nat) : OutClass → Prop
  | ret l' .type1Child => l' = l
  | ret l' .backEdge => l' = l
  | _ => False

/-- A child attaching to the *between* part: a type-2 child to `l`, or anything returning strictly
between `l` and `b`. -/
def AttachesBetween (l : Nat) : OutClass → Prop
  | ret l' .type2Child => l ≤ l'
  | ret l' _ => l < l'
  | _ => False

end OutClass

/-! ### The concrete forest -/

mutual
/-- The out-list of vertex `v` in the tree (empty if absent). -/
def DfsTree.outsAt (v : Nat) : DfsTree → List DfsOut
  | .node w outs => (if w = v then outs else []) ++ DfsOut.outsAtList v outs
def DfsOut.outsAtList (v : Nat) : List DfsOut → List DfsOut
  | [] => []
  | .back .. :: rest => DfsOut.outsAtList v rest
  | .tree _ _ child :: rest => child.outsAt v ++ DfsOut.outsAtList v rest
end

mutual
/-- The depth of `v` in the tree rooted at depth `d₀` (`none` if absent). -/
def DfsTree.depthAt (v d₀ : Nat) : DfsTree → Option Nat
  | .node w outs => if w = v then some d₀ else DfsOut.depthAtList v (d₀ + 1) outs
def DfsOut.depthAtList (v d₀ : Nat) : List DfsOut → Option Nat
  | [] => none
  | .back .. :: rest => DfsOut.depthAtList v d₀ rest
  | .tree _ _ child :: rest => (child.depthAt v d₀).orElse fun _ => DfsOut.depthAtList v d₀ rest
end

/-! ### The walk's edge order -/

mutual
/-- The edges below a tree in the order `finishEdge` reaches them: a vertex's out-edges in sorted
order, each tree edge right after its subtree's edges (PROOF.md §3, Fact D). -/
def DfsTree.edgePostorder : DfsTree → List Nat
  | .node _ outs => DfsOut.edgePostorderList outs
def DfsOut.edgePostorderList : List DfsOut → List Nat
  | [] => []
  | .back e _ _ :: rest => e :: DfsOut.edgePostorderList rest
  | .tree e _ child :: rest => child.edgePostorder ++ e :: DfsOut.edgePostorderList rest
end

/-- The edges of an out-edge's block: its subtree's postorder followed by the edge itself. -/
def DfsOut.block : DfsOut → List Nat
  | .back e _ _ => [e]
  | .tree e _ child => child.edgePostorder ++ [e]

mutual
/-- The blocks of all tree out-edges below a tree. -/
def DfsTree.blocks : DfsTree → List (List Nat)
  | .node _ outs => DfsOut.blocksList outs
def DfsOut.blocksList : List DfsOut → List (List Nat)
  | [] => []
  | .back .. :: rest => DfsOut.blocksList rest
  | .tree e _ child :: rest =>
    (child.edgePostorder ++ [e]) :: (child.blocks ++ DfsOut.blocksList rest)
end

/-- The edge order of the whole forest (`nxtEdgeIdx` counts the back edges in it). -/
def edgePostorderForest (forest : List DfsTree) : List Nat := forest.flatMap DfsTree.edgePostorder

/-- `s` is a subtree of `t`. -/
inductive DfsTree.Sub : DfsTree → DfsTree → Prop
  | refl (t) : Sub t t
  | step {s v outs e cls child} : Sub s child → DfsOut.tree e cls child ∈ outs → Sub s (.node v outs)

/-- Two edge intervals are nested or disjoint. -/
def Laminar (B₁ B₂ : List Nat) : Prop := B₁ ⊆ B₂ ∨ B₂ ⊆ B₁ ∨ ∀ x, x ∈ B₁ → x ∈ B₂ → False

/-- The abstract view of a phase-1 forest; `root` is the first tree's root. -/
def DfsData.ofForest (forest : List DfsTree) : DfsData where
  root := (forest.head?.map DfsTree.v).getD 0
  depth v := (forest.findSome? (·.depthAt v 0)).getD 0
  outs v := forest.flatMap (·.outsAt v)

/-- The phase-1 hypotheses, for the forest actually computed. -/
abbrev DfsForestSpec (g : Graph) (forest : List DfsTree) : Prop := (DfsData.ofForest forest).Spec g

end Spqr
