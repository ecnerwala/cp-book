# Why the SPQR walk is correct — a top-down proof sketch

This is the natural-language proof that the Lean implementation in `lean/Spqr/` (and hence the
C++ it mirrors) computes the canonical SPQR tree. It is written to be formalized: every statement
carries a tag — **[def]** already a Lean definition, **[lemma]** stated/provable with the current
definitions, **[hard]** the genuinely difficult steps. Section 6 maps the statements to Lean files
and child-session work packages.

Notation. `G = (V, E)` is a multigraph with loops and parallel edges, `nv = |V|`, `ne = |E|`.
Phase 1 (`Spqr/Dfs.lean`) produces a DFS forest; `depth v`, `parent`, `anc v l` (the ancestor of
`v` at depth `l`), `T_c` (the subtree of `c`, as a vertex set) are the usual notions. Every
non-tree edge is a *back edge* between a vertex and a proper ancestor (or itself: a loop) —
there are no cross edges **[lemma, DFS]**. For a vertex `c` at depth `d`,
`lowpt1 c < lowpt2 c` are the two smallest depths `≤ d` reached by a back edge from `T_c`
(`Lowvals`, `mergeLowvals` **[def]**; "none" is encoded as `d`).

## 1. The sorted DFS forest

`classify d isTree (lowpt1, lowpt2)` **[def]** labels the out-edge `v → c` (or back edge `v → a`)
of `v` at depth `d`:

| class | meaning |
|---|---|
| `bridge` | tree edge, `lowpt1 c = d+1` — `T_c` never returns |
| `component` | tree edge, `lowpt1 c = d` — `T_c` returns to `v` only: `v` is a cut vertex, `T_c ∪ v` opens a new block |
| `selfLoop` | back edge `v → v` |
| `ret l type1Child` | tree edge, `lowpt1 c = l < d`, `lowpt2 c ≥ d` — the subtree attaches to the rest **only through `v` and `anc v l`** |
| `ret l backEdge` | back edge to `anc v l`, `l < d` |
| `ret l type2Child` | tree edge, `lowpt1 c = l < d`, `lowpt2 c < d` — the subtree also attaches strictly between |

*Terminology.* "type-1 child" (subtree hangs off `{v, anc v l}` only) is the Hopcroft–Tarjan type-1 split; "type-2 child" here just means a returning child that is not type-1. "Type-2 *pair*/split" (§3) is the HT separation-pair notion, which is about first-child chains, not about a single child. Only the relative order "type-1 edges to `l` before type-2 children to `l`" matters for correctness; the type-1-child vs back-edge tiebreak only affects output order.

`OutClass.rank` **[def]** orders them `bridge < component < selfLoop <` all `ret`, and the `ret`
edges by `(l, kind)` with `type1Child < backEdge < type2Child`; `dfsVisit` stably sorts each
vertex's out-list by rank. We call the `ret` edges the *returning* edges of `v`, and a type-1 child
or a back edge a *type-1 edge*: both are, for everything above `v`, just an edge `v — anc v l`
(the type-1 child will be collapsed into exactly that — Lemma 4.3). An ear only "continues"
through type-2 children.

**Lemma 1.1 (spanning) [lemma, DFS].** `dfsForest` visits every vertex exactly once, and every
edge appears exactly once as a `DfsOut` of its shallower endpoint (loops: of their vertex). Hence
`Σ |DfsTree.verts| = nv`, `Σ |DfsTree.edges| = ne`, both without repetition.

**Lemma 1.2 (lowpoints) [lemma, DFS].** The `Lowvals` returned by `dfsVisit` for `c` are the two
smallest distinct depths reached by back edges out of `T_c` (or `(d, d)`); hence `classify`
agrees with the table above.

## 2. Blocks (the `lowval ≥ d` cases)

**Lemma 2.1 [lemma].** Vertex `v` is a cut vertex of its connected component iff it has a
`component` out-edge or is a root with ≥ 2 tree children; the blocks (2-connected components of
the *edge set*; a bridge and a loop are blocks on their own) are: each `bridge` edge, each loop,
and for each `component` edge `v → c` the edge set
`B(v,c) = {v–c} ∪ edges(T_c) ∪ {back edges from T_c}` minus the blocks nested inside `T_c`.

In `finishEdge` these are exactly the `lowval ≥ d` branches: a `bridge` becomes a `Q → I` item
with the child's vertex as the I's single child; a loop a `Q → O`; a `component` edge a
`Q` whose children are the block's top-level ear (everything on the tstack above `origTstack`,
which by Lemma 4.3 is one entry) and the child's vertex item. In all three cases the `Q` is made a
child of `V v`, so the forest/block-cut tree structure (`F → V → Q → … → V → Q → …`) follows from
Lemma 2.1 directly. Nothing else in the walk touches these items.

So from here on fix a block `B` with top vertex `r` (the `v` of its `component` edge, or a root)
and treat it as a 2-connected graph with ≥ 2 edges; the DFS tree restricted to `B` is a DFS tree
of `B` rooted at `r` with the same lowpoints (the deeper blocks have been cut off, and their only
contact with `B` is a cut vertex, which never changes a lowpoint of `B`'s vertices).

## 3. Where the 2-cuts are (the graph theory)

A *separation pair* of `B` is `{a, b}` such that the *separation classes* — equivalence classes of
`E(B)` under "joined by a path avoiding `a` and `b` internally"; every edge `a–b` is its own class
— number at least 2 after merging… (we use the standard SPQR convention: ≥ 2 classes with ≥ 2
edges, or ≥ 3 classes). All statements in this section are about the fixed sorted DFS tree of `B`.

**Fact A (ancestor/descendant) [lemma, graph theory].** If `{a, b}` is a separation pair of `B`
then one of `a, b` is an ancestor of the other.
*Proof.* If not, `T_a` and `T_b` are disjoint and `R := B − T_a − T_b` is connected by tree edges
(it contains the root and both parents). Each child subtree of `a` has `lowpt1 < depth a` (else
`a` would be a cut vertex of `B`), so it has a back edge to a proper ancestor of `a`, which lies
in `R` (it is not `b`, and not in `T_b`, else `a ∈ T_b`). Same for `b`. There is no edge `a–b`
(adjacent vertices are ancestor/descendant). So all of `E − {a,b}` is one class.

Fix a separation pair `a = anc b l`, `depth b = d > l`. Let `a'` be the child of `a` toward `b`,
`P = ` the tree path strictly between `a` and `b` (depths `l+1 … d−1`, possibly empty).
`B − {a, b}` decomposes into: **above** = `B − T_{a'}` (plus `a`'s other children: all together
one class, since `B` is 2-connected and `a` is not a cut vertex — the children of `a` other than
`a'` return above `a`), **between** = `T_{a'} − T_b` with its hanging subtrees, and for each
child `c` of `b` the subtree `T_c`.

**Fact B (classes at an ancestor/descendant pair) [lemma, graph theory].** Each child `c` of `b`
has `lowpt1 c < d` and attaches to: *above* iff `lowpt1 c < l`; *between* iff some back edge of
`T_c` lands in `(l, d)`, i.e. `lowpt1 c > l` or (`lowpt1 c = l` and `lowpt2 c < d`); and only to
`a` iff `lowpt1 c = l` and `lowpt2 c ≥ d` — i.e. iff `b → c` is a **type-1 child** of `b` with
`lowval l`. In the sorted out-list of `b`, the children attaching above form a **prefix**
(`lowval < l`), then come the type-1 edges with `lowval = l`, then the children attaching between
(a **suffix**: `lowval = l` type-2, then `lowval > l`). `b`'s own back edges do not merge classes
(they are incident to the removed `b`).

Consequently a separation pair `{a, b}` is of one (or both) of two kinds:

* **type 1**: `b` has a type-1 edge `b → c` with `lowval l` (so `E(T_c) ∪ {b–c} ∪ back edges of
  T_c` is a class, and the rest is non-empty); or
* **type 2**: no child of `b` attaches both above and between, *between* does not attach above
  (no back edge from `P ∪ (hanging subtrees of P) ∪ b` lands at depth `< l` — `b`'s own back
  edges above `a` are allowed, they are their own classes/`above`), and both sides are non-trivial.
  Then `above ∪ prefix` and `between ∪ suffix` are unions of classes separated by `{a, b}`.

**Fact C (type-2 pairs live on first-child chains) [lemma, graph theory; Lean:
`type2_first_out`].** Let `{a, b}` be a pair whose *above* and *between* parts are separated by
`{a, b}` (no edge of *above* is in the class of an edge of *between* — this is what "type-2"
contributes), and let `a` have a parent (`a` is not the root). Then `b` is reached from `a` by one
tree edge to some child `a'` followed by a chain of **first** out-edges in sorted order: for each
`x` on `P ∪ {b}` with parent `p ≠ a`, `p → x` is the *first* out-edge of `p` — no back edge, no
type-1 edge and no other child of `p` precedes it.
*Proof.* `a` is not a cut vertex, so `T_{a'}` has a back edge to a proper ancestor of `a`. If it
left from `T_{a'} − T_x` it would join *between* to *above*; hence it leaves from `T_x`, and
`lowpt1 x < l`. An out-edge of `p` sorted before `p → x` has `lowval ≤ lowpt1 x < l`, so it (a
back edge of `p`, or the subtree of another child of `p`) reaches a vertex above `a` from `p ∈ P`,
again joining *between* to *above*: contradiction.

*Caveat (root case).* The hypothesis that `a` is not the root is necessary. If `a` is the root
(`l = 0`) the *above* part is empty, `{a, b}` separates only through parallel `a–b` edges (a
bond, ≥ 3 classes), and the first-child claim fails: `T_x` returns to `a` only (`lowpt1 x = 0`),
so an earlier out-edge of `p` with `lowval 0` — a back edge to `a`, a type-1 child returning only
to `a`, or a type-2 child with `lowpt1 = 0` that ties with `p → x` in rank — is possible. (The
"type-1 edges returning exactly to `l` may precede `p → x`" exception of an earlier version of
this note only arises there; for a non-root `a` nothing precedes `p → x`.) The walk must handle
the root/bond case separately — via `firstOccurrence`/`firstIdx`, see 4.4.

Similarly, for `a` itself: `a'` need not be `a`'s first child — a type-2 pair is attached to `a`
from whichever child.

**Fact D (laminarity) [lemma, graph theory; Lean: `Spqr/Proofs/{Postorder,Interval}.lean`].**
Order the edges of `B` by *postorder* σ: `dfsVisit`'s finishing order — a vertex's out-edges in
sorted order, each tree edge placed right after its subtree's edges. (This is the order in which
`finishEdge` is called, and `nxtEdgeIdx`/`firstIdx` count the back edges in it.) In Lean σ is
`DfsTree.edgePostorder` / `edgePostorderForest`, defined purely from the `DfsTree`; the *block* of
a tree edge `p → c` is `σ(T_c) ++ [p → c]` (`DfsOut.block`), and `DfsForestSpec` additionally
assumes the forest's vertex list and edge list are duplicate-free. Then every separation class is
a contiguous interval of σ **up to S/P multiplicity**:

* *Type-1 classes* (`type1_class_interval`): for `a = anc b l` and a child `c` of `b` with
  `b → c` classified `ret l type1Child`, the class of `b → c` is exactly the block of `b → c`,
  an interval of σ. Blocks of tree edges are pairwise nested or disjoint
  (`type1_classes_laminar`, structural: `blocks_laminar_forest`).
* *Type-2 classes* (`type2_class_interval`): under the hypotheses of Fact C (`a` not the root,
  `a'` its child toward `b`, *above* and *between* separated by `{a, b}`), the class of `a → a'`
  is `DfsTree.type2Block l (a → a') T_{a'} T_b` = `σ(dropWhile (rank ≤ rank(ret l backEdge))
  outs(b)) ++ (σ(T_{a'}) − σ(T_b)) ++ [a → a']`, a suffix of the block of `a → a'` and hence an
  interval of σ. The proof uses Fact C: every tree edge on the path `a' ⇝ b` is first in its
  out-list, so `σ(T_b)` is a prefix of `σ(T_{a'})` (`type2_chain_prefix`).
* *Nesting, type-2 vs type-1* (`type2Block_laminar_block`): the type-2 class of `{a, b}` and the
  block of any tree edge `p → w` are nested or disjoint **unless** `w` lies strictly below `a'` on
  the tree path to `b`. That exception is exactly the S caveat: a cycle (S) or bond (P) can be
  split at any point, so the maximal pieces of the canonical decomposition, not the individual
  classes, are what is unique; e.g. in a cycle `x – a – v – w – b – x` the type-2 class
  `{a–v, v–w, w–b}` of `{a, b}` and the type-1 class `{v–w, w–b, b–x}` of `{x, v}` overlap.
* *Nesting, type-2 vs type-2*: **not proved in Lean.** Structurally the two `type2Block`s need
  not be nested (e.g. same `a, a'`, `b₂` below `b₁`, with out-edges of `b₁` sorted before
  `ret l backEdge`, or `b₂` in a child of `b₁` sorted after it while `b₂` itself has out-edges
  sorted before); ruling these out needs the semantic facts of §1 (in a block no child of `b₂`
  below a non-returning child of `b₁` can return to depth `≤ l`), which were not formalised.

This is the statement that makes a *stack* the right data structure: the open intervals at any
time are nested.

## 4. The walk: ears and the tstack

### 4.1 Ears

Define, in the sorted DFS tree of `B`, the *first-child chain* from a vertex `c`: `c = x_1`, and
`x_{i+1}` is the target of the first `ret` out-edge of `x_i` as long as that edge is a **type-2
child**; the chain ends at `x_k` whose first `ret` edge is a type-1 edge (back edge or type-1
child) to depth `l`. Note `lowval` is constant along the chain: `lowpt1 x_1 = … = lowpt1 x_k = l`
(the first edge has the minimum `lowval` and it is attained by the chain).

An **ear** is: a non-first `ret` edge `u → c` of some vertex `u` (its *base*), or the `component`
edge `r → c` opening the block, together with the first-child chain from `c` and the type-1 edge
closing it. Its terminals are `u` and `anc u l` (for the opening ear, `r` and `r` itself: the
first cycle). Every vertex except `r` is interior to exactly one ear (the one whose chain contains
it), and every edge belongs to exactly one ear (chain edges and the closing edge to the chain's
ear; a non-first back edge is its own trivial ear; a non-first tree edge starts a new ear). This is
the open ear decomposition of `B` induced by the DFS **[lemma]**: each ear is a path between two
vertices already present, because the chain's lowpoint is an ancestor of the base.

### 4.2 What a tstack entry is

A `TEntry (vStart, topDepth, firstIdx, spans)` **[def]** denotes the edge set
`E(t) = ⋃ { subtree edges of item | item ∈ spans.1 ++ spans.2 }` (a `Q` item is its edge, a `V`
item is empty, a node item is the union over its children) together with the vertex set touched.
Its *top* is `anc cur topDepth` (`stackVerts[topDepth]`), its *bottom* is `vStart`.
`firstIdx` is the postorder index (back-edge count) at which the entry was created.

**Invariant W (per ear) [hard — the central lemma].** Consider the moment the walk is at vertex `x`
of depth `d`, about to process its out-edge `o` (in `walkOut`). Let `x` be interior to the ear
`ε` with base `u`, chain `x_1 … x_k`, say `x = x_i`. The tstack is
```
  (tstack as it was when ε's base edge was entered)        -- entries of enclosing ears
  ++ Stack(ε, i)                                             -- this ear's open entries
```
where `Stack(ε, i)` consists of, from bottom to top:
1. for `j = i, i−1, …` (going *up* the chain, i.e. the entries pushed deeper are on top — see
   below), the entries contributed by the already-finished non-first edges of `x_j`, each a single
   entry `(x_j, l_j, …)` with `E = ` the whole sub-ear (Lemma 4.3), **merged** with its neighbours
   whenever the merge rules 4.4 fired — merging never crosses from `Stack(ε, i)` into the
   enclosing part, except for the case documented in 4.4(c);
2. the vertex item `V x_j` of each chain vertex whose first `ret` edge has been *finished* (pushed
   by `walkOut` before a type-1 first edge, or by the `!hasVert` branch of `finishEdge` after a
   type-2 first edge);
3. every entry `t` satisfies: all vertices of `E(t)` other than `anc topDepth`, `vStart`, and the
   chain vertices `x_i, …` currently on the DFS path have **all** their incident edges in `E(t)`
   (they are *finished*: interior), and `E(t)` is connected; entries are pairwise edge-disjoint and
   together with the entries of the enclosing ears partition the edges processed so far in `B`,
   minus those already placed in finished items;
4. `topDepth` is weakly *increasing* from bottom to top within the group pushed by one vertex (the
   out-edges are sorted by `lowval`), and the `firstIdx` of entries are increasing from bottom to
   top.

Statement 3 is the "every closed item is a separation class" direction; statement 1 is what makes
the induction go through: an ear's walk is self-contained.

### 4.2b Proof architecture: soundness and completeness per step

The decomposition theorems do not need the full interval structure. The per-step invariant on the
global tstack has two local parts per entry `E = (vStart, topDepth, …)` with edge set `edges(E)`:

* **connected**: `edges(E)` is connected, and so is each item already closed;
* **2-attached**: every edge of the block not in `edges(E)` meets `edges(E)` only at
  `{vStart, stackVerts[topDepth]}`.

*Soundness*: every `mergeTstackTops` joins two entries sharing a terminal, so the union is again
connected and 2-attached (with the terminals computed by the `min topDepth` rule). *Completeness*:
every item that `finishEdge` closes (S/P/R, or the single entry handed up as a virtual edge) is a
2-attached connected edge set, i.e. cut off by a genuine 2-vertex separation of the block, and the
closed items partition the edges. This gives `Items.Tree`, `Items.Endpoints` and the S/P/Q/I/O shapes
without ears. The ear view (§4.1–4.3) is used only for the *stack-shape* facts these steps rely on
(`size ≥ origTstack + 3`, `nxt` is the chain/vertex entry, the merge loops never cross a frame) and
for maximality (§4.5), where Facts C–D are needed.

### 4.3 Lemma (ear collapse) [hard, the inductive step]

After `finishEdge u d (u → c)` for a **non-first** `ret` edge `u → c` with `lowval l`, the tstack
is exactly what it was before `walkOut` for that edge, plus **one** entry
`(vStart = u, topDepth = l, E = {u–c} ∪ edges(T_c) ∪ back edges out of T_c)`.
If `u → c` is type-1 the entry's spans are a single item `S`/`R`/`P` with `vs = {u, anc u l}`
(the "virtual back edge"); if type-2 it is an open entry whose spans still list `V` items of
chain vertices of the sub-ear that attach strictly between `l` and `d` (statement 3 holds: those
are the only unfinished vertices).

*Proof shape.* Induction over the sub-ear's chain `c = x_1, …, x_k` and Invariant W. When
`finishEdge` runs for `u → c`:
* `pushEdgeTstack c d e` pushes the tree edge as `(c, d, {u–c})`.
* **Loop 1** (`nxt.topDepth ≥ d`): every entry of `Stack(ε', ·)` whose top is `u` or deeper gets
  closed: `topDepth > d` means the entry is a piece hanging below `u` on the chain (its top is a
  chain vertex `x_j`): it is **S**-composed (series) with the current piece — the two pieces share
  the single vertex `x_j`, which is now finished, so their union is again a 2-terminal piece with
  `vStart` unchanged and top `u` (`setStackDir` records the side). `topDepth = d` means it returns
  to `u` itself: if its bottom is the same `vStart` it is **P**arallel to the current piece
  (same two terminals), else the two pieces share the top `u` and have different bottoms `y ≠ y'`,
  and — this is the key case — `y` is the bottom of the *current* piece, which (by statement 3,
  with `y` no longer on the path) now has all its edges inside the union: so the union is a
  2-terminal piece `{y', u}` whose skeleton is not a cycle or bond: **R**. Each closure runs
  `finishTstackTop`, creating the item with `vs = {u, bottom}` and the entry's items as children.
  Loop 1 ends with `cur` = one piece with top `u`, bottom = some chain vertex `y`, and everything
  below has `topDepth < d`.
* **Loop 2** (`firstIdx > firstOccurrence[d]`): `firstOccurrence[d]` is the postorder index of the
  first back edge *into* `u` from `T_c` (reset when the tree edge was entered). Entries created
  after that back edge, i.e. the sub-ear pieces lying between the piece that first returned to `u`
  and `cur`, cannot be separated from `cur` by any pair `{anc, u}`: they are glued into `cur`
  (this is the Fact C situation — a type-2 chain that is interrupted; the pieces returning to `u`
  from inside force the merge). After Loop 2 `cur` is the maximal piece with top `≤ d` built from
  the sub-ear, with `isSingle = false` iff anything was merged.
* **Type-1 (`hasVert`, `isType1`)**: by Fact B everything in `T_c` returns to `l` or `≥ d`; the
  stack above `origTstack` is now `[V y, (y, l, back-edge piece), cur]` — exactly three entries
  (the assertion in the C++). `maybeUnwrapNxt (S or R)`, merge the back-edge piece, merge `V y`,
  set `vStart := u`, `finishTstackTop` → the single item with `vs = {u, anc u l}`. `isSingle`
  decides S vs R: a single chain with one return is a cycle through `u, …, y, anc l`; anything
  merged in Loop 2 adds a chord → R (and Loop 1's R/P closures are already items).
* **Type-2 (`hasVert`, `!isType1`)**: merge everything above `origTstack + 3`, then the two more
  merges, set `vStart := u`: one open entry with `topDepth = l`. It is not closed because its
  piece attaches at `u`, `anc l` *and* at the depths in `(l, d)` where `T_c` returns.
* **First edge (`!hasVert`)**: the P-check / `V u` push: `u` becomes interior to the enclosing ear
  (its `V` item goes onto the stack), the sub-ear's open entries stay open (they are now the
  enclosing ear's `Stack(ε, i−1)` contributions). Nothing is merged across.
* **P-check** (any type-1 edge): if `nxt` is `(u, l)` as well — a previous type-1 edge of `u` to
  the same ancestor — the two are parallel: `maybeUnwrapNxt .P` reuses an existing `P` item
  (canonicality: no `P` under `P`) or allocates one, and merges.

The statements to prove for each bullet are of the form "the union of the merged entries is a
piece with the stated terminals" (soundness, gives `Items.Endpoints.separation`) and "nothing that
could still be separated was merged" (maximality, see 4.5).

### 4.4 Side bookkeeping is ordering only

`stackDir`, `edgeDir`, `setSides`, the choice of `spans.1` vs `spans.2`, and the final fold
`setSides (!edgeDir) (spans.1 ++ spans.2) []` implement the st-ordering layered on the
decomposition (children of a node come out in s–t order; this is what the planar variant and the
"center edge" adjacency need). They never affect *which* items are merged, only the order of `ch`
lists. The item-level spec (`ItemSpec.lean`) therefore quantifies over `ch` as a list up to
permutation, and the st-order is a separate theorem to be added later.

What the side bookkeeping computes (the invariant for the st theorem): it is exactly the
Even–Tarjan st-numbering of the block, read off the lowval-sorted DFS tree. In ear terms (4.1)
the bottom of an ear is innermost; as the walk returns upwards every further piece (type-1
child, back edge, finished sub-ear) is either *prepended* or *appended* to the current
sequence, the side being chosen by `stackDir`/`edgeDir` via `setSides`. So a tstack entry
stands for a list of the final ordering **with a hole in the middle** for the deeper, not yet
finished part: `spans.1` is what lies left of the hole, `spans.2` what lies right of it, and the
two terminals of the entry are the two ends. When the hole closes (`mergeTstackTops`,
`finishTstackTop`) the sides concatenate around the inner piece and become the `ch` order of the
item; every interior vertex then has a neighbour on each side, which is the st property
(`Items.StNumbered`, checked empirically by `cpp/check_st_planar.cpp`).

Two facts make most entries one-sided (they are what the one-sided asserts in the C++ rely on):
along an ear every vertex has the same lowval (the ear's return depth), so `stackDir` is constant
along the first-child chain and every piece attached along the ear goes to the *same* side — the
hole is two-sided only at ear boundaries, where a finished inner ear gets wrapped from both
sides; and type-2 separation pairs always lie along a single ear (Fact C), so the piece a type-2
split cuts off is uniformly on one side and its entry has the empty-side shape. Only type-1
closes / ear boundaries produce genuinely two-sided entries.

### 4.4b Planar variant: what a tstack entry stores

In the planar variant a tstack entry represents its piece as an interior spine along the ear
plus the two outer-face boundary walks along the outside of the piece; only the exposed ends of
those two walks (quarter-edges) are stored, together with, on each side, a linked list of the far
ends of the back edges leaving the piece on that side. Conceptually the whole upper DFS stack
(the back-edge targets, i.e. the ancestors above the piece) is contracted into a single vertex at
which the exposed ends also live. Embedding a new piece or back edge is choosing which boundary
walk it goes on (consistent with its st side); the piece is non-planar when both sides are
already blocked. The invariant to prove is therefore: the piece has a planar embedding with both
terminals on the outer face, and the exposed ends + side lists describe exactly its outer-face
boundary split at the terminals; merges glue embeddings along the shared terminal.

### 4.5 Maximality and R skeletons [hard]

Soundness (each item is a separation class, Fact B direction "⇐") comes from statement 3 of
Invariant W. For the canonical decomposition we also need: **no S/P/R item could be split
further**, equivalently every separation pair of `B` appears as a pair of vertices of some `S`
cycle, the two vertices of a `P`, or the endpoints of a virtual edge — and hence
**every R skeleton is 3-connected**. By Facts A–C every separation pair `{a, b}` is a type-1 or
type-2 pair with `b` a chain vertex and `a = anc b l`; the walk pushes `(b, l)`-entries for exactly
these candidates (type-1 edges of `b`; the first-child chain through `b`), and the only events that
*destroy* a candidate (Loop 1's R case, Loop 2) happen exactly when Fact B/C say the pair is not a
separation pair (a child attaching both above and between; between attaching above). So the pairs
that survive to be closed are all separation pairs and every separation pair is closed. The
3-connectivity of an R skeleton then follows from: a separation pair of the skeleton would be a
separation pair of `B` (virtual edges stand for connected pieces attached at their two ends, so a
cut of the skeleton lifts to a cut of `B`), which would have been closed as an S or P.

The graph-theoretic half is done: `Proofs/SepPairExhaust.lean` proves `sepPair_iff` — over a block
with its sorted DFS tree, `{a, b}` with `a` an ancestor of `b` is a separation pair iff it is a
`Type1Pair` (a type-1 child of `b` to `depth a` plus two further edges, or a bond: two parallel
`a–b` edges plus a third edge) or a `Type2Pair` (`a` not the root, `b` strictly below `a`'s child
`a'`, no child of `b` returning both above `a` and into `(a, b)`, no back edge from `T_{a'} − T_b`
above `a`) — and the corollary `three_connected_of_no_split`: a block with no type-1 and no type-2
pair has no separation pair at all. What remains is the walk-side half below.

Formalizing 4.5 directly on the walk is the research-scale part. The planned route avoids a
separate "all separation pairs are found" theorem by proving the stronger local statement in 4.3:
each merge rule fires iff the corresponding candidate is not a separation pair, by Fact B/C applied
to the lowpoints that the `OutClass`es carry.

## 5. Phase 3: relabel

`relabelTree` **[def]** takes the item array and produces `SpqrTree`. It is a plain preorder walk:
items are numbered in DFS order from `rootItem`, `par`/`subtreeEnd` are the standard preorder
facts, `chBounds`/`ch` is the CSR of children, `nvBounds`/`nodeVerts` lists for each node its
skeleton vertices (`vs` endpoints then the `V` children's vertices), `neBounds`/`nodeEdges` its
skeleton edges (one per non-`V` child, as a *virtual* edge with the child's `vs`; plus the node's
own cap edge if it has a parent), `twin` pairs each child's virtual edge in the parent with the
cap edge of the child. Everything in `Spec.lean`'s `WF` is a statement about this encoding of
`Items.Tree`, and `Represents.twin_glue`/`q_endpoints` are `Items.Endpoints` transported along the
relabeling **[lemma, mechanical but large]**; `r_three_connected` and `canonical` are
`Items.Shapes` transported.

## 6. Lean plan (what is proved where)

| statement | file | status |
|---|---|---|
| graph / DFS / items / walk / relabel definitions | `Graph.lean … Build.lean` | def |
| output spec `SpqrTree.WF`, `Represents` | `Spec.lean` | def |
| phase-2 contract `Items.WF` | `ItemSpec.lean` | def |
| 1.1, 1.2 DFS spanning + lowpoints | `Proofs/Dfs.lean` (`dfsForest_spanning`, `dfsForest_wf`, `classify_child_*`) | proved |
| 2.1 blocks ↔ `lowval ≥ d` branches (`sameBlock_iff`, `blockRoot_cut`, `ret_child_sameBlock`) | `Blocks.lean`, `Proofs/Blocks.lean` | proved |
| Facts A–C (`sepPair_comparable`, type-1 class / above-between / sorted prefix-suffix, `type2_first_out`) | `SepPair.lean`, `Proofs/SepPair.lean` | proved |
| Fact D (laminar intervals): `edgePostorder`, type-1 classes are intervals + laminar, `type2_class_interval`, `type2Block_laminar_block` | `Proofs/{Postorder,Interval,Type2}.lean` | proved (type-2 vs type-2 nesting open, see §3) |
| 4.5 exhaustiveness: separation pairs of a block = type-1 ∪ type-2 pairs (`sepPair_iff`, `three_connected_of_no_split`) | `SepPairExhaust.lean`, `Proofs/SepPairExhaust.lean` | proved |
| ear-structured walk (`descend`/`ascend` over chain `Frame`s) | `Ear.lean` | def |
| `walkEarTree = walkTree` (`walkEarTree_eq_walkTree`, `walkEar_eq_walk`) | `EarSpec.lean` | proved |
| frame rule `walkTree_local` via `Lifts`/`Sim` simulation (`Sim.closeEars`, `Sim.mergeLate`, `Sim.finishRest`, `Sim.finishBoundary` proved) | `Sim.lean`, `Frame.lean`, `EarSpec.lean` | `Sim.closeVert` (loop 3 mechanics) and `Sim.walkTree` (invariant threading) admitted; rest proved |
| typing/allocation part of `Items.WF` (`Items.Tree` sizes/types, I/O leaves, `vs_shape`, `vs_lt`): `walk_typing` | `WalkTyping.lean` | proved (`walk_q_children` sorry: needs span shape) |
| §4.2b walk invariant (`EntryInv`, `Inv`): closure lemmas (`GraphLemmas.lean`), `mergeTstackTops_sound`, `finishTstackTop_complete`, `finishEdge_back_inv` (under explicit stack-shape hyps `BackCheckOk`) | `GraphLemmas.lean`, `WalkSpec.lean` | proved; `finishEdge_inv` (remaining branches), `walkTree_inv`, `walk_nodes_partition` sorry |
| Invariant W, Lemmas 4.3/4.4 (`earOut_one_entry`, `ascend_frame_one_entry`) | `EarSpec.lean` | sorry / hard |
| 4.5 maximality / R 3-connected | `spqrTree_r_three_connected` | hard |
| 5 relabel: `Items.WF → WF ∧ Represents` | `relabelTree_wf`, `relabelTree_represents` | sorry |
| 7 st-order spec `StOrder`, `Items.StNumbered`, split `spqrTree_st = relabel_st ∘ walk_st` | `StSpec.lean` | def / proved split |
| 7 relabel-side: `vchildren_nv_increasing`, `orderedChildren_sorted`, `edgeChildren_dominance` | `StSpec.lean` | proved |
| 7 relabel-side: `layoutNode_r_bracket`, `relabel_st` | `StSpec.lean` | sorry |
| 7 walk-side: `WalkState.StInv`, data lemmas `pushTstack_onSide`, `merge_onSide`, `fold_onSide` | `StWalk.lean` | def / proved |
| 7 walk-side: `chain_stackDir_const`, `finishEdge_topClosable`, `finishEdge_stInv`, `finishTstackTop_stItem`, `walk_st` | `StWalk.lean`, `StSpec.lean` | sorry |

Work packages for child sessions, in dependency order:
* **DFS**: 1.1, 1.2, no cross edges, `lowpt` characterization of `OutClass`.
* **Relabel**: 5, independent of the walk (assumes `Items.WF`).
* **Ear reformulation**: the iterative `walkEar` intermediate and its equality to `walkTree`.
* **Graph theory**: Facts A–D on an abstract sorted DFS tree.
* **Walk invariant**: Invariant W + Lemma 4.3 on `walkEar`, giving `Items.Tree`, `Endpoints`.
* **Shapes / maximality**: `Items.Shapes`, R 3-connectivity via 4.5.
* **st-order**: §7, `StSpec.lean` / `StWalk.lean`; relabel side first, then the walk invariant.
* **planarity**: §8, `Planar.lean` / `PlanarWalk.lean` / `PlanarRelabel.lean` / `PlanarEmbed.lean` / `PlanarSpec.lean`; projection theorems proved, embedding soundness/completeness admitted.

## 7. st-ordering

The side bookkeeping of §4.4 computes, for every block, the Even–Tarjan st-numbering of the
lowval-sorted DFS tree: the bottom of an ear is innermost, and as the walk returns upwards each
other subtree or back-edge piece is either prepended or appended to the sequence built so far.

### 7.1 The output statement (`StSpec.lean`)

For a node `i` with node-verts `[s, e)` and skeleton `skeleton i` (the `nvs` of its node-edges):

* `SpqrTree.StNumbered i`: every skeleton edge has `nvs.1 < nvs.2`, and every interior node-vert
  `s < nv < e-1` has a skeleton edge to a smaller and one to a larger node-vert.
  So the node-vert order is an st-numbering of the skeleton with `s` first and `e-1` last.
* `SpqrTree.EdgeDominance i`: among the non-cap node-edges (`ownEdges`, i.e. `nodeEdgesOf i`
  without position 0 when `hasCap i`), a later edge with distinct `nvs` never has both endpoints
  `≤` those of an earlier one.
  The cap `(s, e-1)` is excluded: it dominates every other edge but sits at position 0.
* `SpqrTree.AdjBracket nv`: adjacency row `2 nv` holds the incidences with `destNv < nv`, row
  `2 nv + 1` those with `destNv > nv`, each with non-increasing `destNv` (ties only for parallel
  edges in P nodes).
  Rows `2nv ++ 2nv+1` are the "center is longest" bracket order of the header comment.
* `SpqrTree.StOrder` bundles the three (st-numbering for S/P/R nodes, dominance for every node,
  brackets for every node-vert of a node with more than one vertex).
* `vchildren_nv_increasing` **[proved]**: V children appear in `children i` in increasing
  node-vert order; this is already contained in `WF.nv_layout` (the V children are the middle of
  the node-vert list, in child order) and is just read off it.

`spqrTree_st : (g.spqrTree tern vo eo).StOrder` **[proved from the two admitted phase theorems]**
follows the same route as `spqrTree_wf`: `spqrTree_eq`, then `relabel_st` applied to `walk_st`
and `walk_items_wf`.

### 7.2 The item-level contract `Items.StNumbered`

For an S/P/R item with `vs = (some s, some t)` let `vertList = s :: (vertices of the V children,
in `ch` order) ++ [t]` and let the skeleton edges be `(s, t)` plus `virtualEdges` (the `vs` of
the non-V children).
`Items.StItem i` says `vertList` is an `Items.StList` of these edges (distinct vertices, every
edge inside the list and non-loop, every interior vertex with an edge to an earlier and to a later
position) and every virtual edge is oriented from the earlier to the later position.
`Items.StNumbered` quantifies this over all S/P/R items.
The split is `walk_st : Items.StNumbered (walk …).items` (§7.4) and
`relabel_st : Items.StNumbered items → Items.WF g items → (relabelTree g items).StOrder` (§7.3).

### 7.3 Relabel: from `ch` order to the output order

* S and P nodes keep `ch` (`orderedChildren_eq_of_ne_R` **[proved]**), so `StNumbered` is the
  item fact transported along `nv_layout`, and dominance is vacuous for S (the edges are a path
  `s → … → t`) and a matter of equal endpoints for P.
* R nodes: `orderedChildren` stably sorts `ch` by `loc`, the sum of the two endpoint positions of
  an edge child (`2 · position` for a V child).
  `orderedChildren_sorted` **[proved]**: the result is a permutation of `ch` sorted by `loc`
  (`List.mergeSort_perm`, `List.pairwise_mergeSort`).
  `edgeChildren_dominance` **[proved]**: if `p` precedes `q` in a list sorted by endpoint sum and
  `q ≠ p` then `q` does not dominate `p` componentwise — equal sums with componentwise `≤` force
  equality.
  Positions are those of `vertList`, so by `Items.StItem` every edge is oriented low → high and
  the V children of the sorted list are still in increasing position (`loc` is monotone in the
  position of a V child, and the sort is stable), which is what `nv_layout` then records as the
  node-vert order.
* Adjacency rows (`layoutNode_r_bracket` **[sorry]**): R nodes count each edge `(a, c)`, `a < c`,
  once for row `2 a + 1` (higher neighbours of `a`) and once for row `2 c` (lower neighbours of
  `c`); the counts are kept one slot up (`2 a + 2`, `2 c + 1`) so that the prefix sum turns slot
  `r + 1` into a cursor that starts at the beginning of row `r` and, after the fill, ends at its
  end — the final `adjBounds`.
  The fill walks `edgeChildren.reverse` and advances the two cursors of each edge, so the rows are
  filled front to back in *reverse* dominance order: row `2 nv + 1` receives the edges `(nv, c)`
  with `c` decreasing (among edges with the same first endpoint, dominance order is increasing
  `c`), and row `2 nv` the edges `(a, nv)` with `a` decreasing.
  The cap is written last, at node-edge position `neSt` (position 0, `capNe`), with incidences at
  the first slot of row `2 s + 1` and the last slot of row `2 (e-1)`; it is the only edge outside
  the dominance order.
  The remaining work is mechanical: the layout is a loop over `Layout` with `inc`/`setNe`, and the
  statement has to be read through `adjBounds`/`adjDat` offsets.

### 7.4 Walk: the hole invariant (`StWalk.lean`)

A tstack entry `t` represents a contiguous piece of the eventual st-ordering of its block *with a
hole*: `spans.1` is the finished part to the left of the hole, `spans.2` the finished part to the
right, and the hole holds the deeper, not yet finished part — the open DFS path below the entry
and whatever is still on the tstack above it.
Reading a stack segment top-first, `nest (t :: rest) = t.spans.1 ++ nest rest ++ t.spans.2`
(`TEntry.wrap`) is the piece with the hole filled in by the entries above.

The side a new piece goes to is decided by `setSides dir a b = if dir then (b, a) else (a, b)`:

* `pushTstack` puts its single item on the `stackDir[topDepth]` side (`pushTstack_onSide`
  **[proved]**),
* `mergeTstackTops` wraps `a` around `b`: `(b.spans.1 ++ a.spans.1, a.spans.2 ++ b.spans.2)` —
  the newer entry `a` is outside, the older `b` inside, on both sides,
* leaving a type-2 child folds the entry to one side of `!edgeDir`:
  `setSides (!edgeDir) (spans.1 ++ spans.2) []` (`fold_onSide` **[proved]**),
* closing an entry (`finishTstackTop`) reads `ch := getSide spans topDir` and
  `vs := makeVs = setSides topDir (top, vStart)`: so `vs = (top, vStart)` with `ch = spans.1`
  when `stackDir[topDepth] = false`, and `vs = (vStart, top)` with `ch = spans.2` otherwise.
  Either way the children list reads from the upper terminal `top = stackVerts[topDepth]` towards
  `vStart`, and `vertList = s :: V children ++ [t]` is the entry's sequence between its two
  terminals.

**Ear uniform side (`ear_uniform_side`).** Along the first-child chain of an ear every vertex has
the ear's lowval `l` (§4.1), so `stackDir d' = !stackDir l` is constant along the chain and every
piece attached along the ear is pushed on the same side.
The data half is `merge_onSide` **[proved]**: merging two entries that are on the same side stays
on that side.
The semantic half, `chain_stackDir_const` **[sorry]**: an entry whose pieces were all attached
along a chain of constant `stackDir` is `OnSide` that direction.
Consequently the spans of an entry are one-sided inside an ear; genuinely two-sided entries arise
only at ear boundaries, when a finished inner ear (`spans.1`, `spans.2` both non-empty after the
wrap) is enclosed by pieces from both sides.

**Type-2 entries are one-sided (`type2_entry_one_sided`).** A type-2 separation pair `(a, b)` lies
along a single ear (`b` is reached from `a`'s child by first children, Fact C), so the piece a
type-2 split cuts off was attached uniformly on one side; the fold `fold_onSide` makes this
literal in the data.
Only type-1 closes and ear boundaries produce two-sided entries.

**The invariant.** `WalkState.StInv s d ord`, for the eventual numbering `ord` of the block:
every entry with `topDepth ≤ d` has distinct, `ord`-increasing V items on `spans.1 ++ spans.2`
(`StSides`); the top entry is in addition oriented — its V items lie after `stackVerts[topDepth]`
when `stackDir[topDepth] = false` and before it otherwise (`StEntry`); and every closed S/P/R item
is `StItem`.
`WalkState.TopClosable` is the intrinsic st-property of the top entry at a close: `entryVertList`
(terminals from `makeVs` around the V items of `getSide spans topDir`) is an `StList` of
`entryEdges`, with every virtual edge oriented.
`finishTstackTop_stItem` **[sorry, mechanical]** turns `TopClosable` into `Items.StItem` for the
closed item; `finishEdge_topClosable` **[sorry]** says the top entry is closable whenever
`finishEdge` closes it, and `finishEdge_stInv` **[sorry]** is the preservation of `StInv` through
`finishEdge`; `walk_st` follows from `StInv.items` at the end of the walk.

**What was validated empirically** (trace of the C++ walk on `gen.py` seeds 0..149, every
tstack snapshot at the start of an out-edge; the Lean walk is byte-identical on these inputs):
the V items of every entry are distinct and increasing in the final order of their parent item
(`StSides`, 0 violations); the top entry is oriented towards `stackVerts[topDepth]` as in
`StEntry` (215 live cases, 0 violations; non-top entries whose `topDepth` slot has since been
reused by a sibling can violate it, which is why `StEntry` is stated for the top entry only);
every closed S/P/R item satisfies `Items.StItem` (986 items, 0 violations).
The naive strengthening "every entry's `nest` reading is an st-list with the hole at the far end
from the terminal" is **false** when the open path changes direction (an entry attached across a
dir-0 sub-chain and a dir-1 chain has path vertices on both sides of its items, e.g. seed 20 of
`gen.py`), so the hole must be modelled as the path itself, ordered by `stackDir` per depth;
that refinement is left to the preservation proof.

### 7.5 Work packages

| lemma | file | status |
|---|---|---|
| `SpqrTree.StNumbered`, `EdgeDominance`, `AdjBracket`, `StOrder` | `StSpec.lean` | def |
| `vchildren_nv_increasing` | `StSpec.lean` | proved |
| `Items.StList`, `Items.StItem`, `Items.StNumbered` | `StSpec.lean` | def |
| `spqrTree_st` from `walk_st`, `relabel_st`, `walk_items_wf` | `StSpec.lean` | proved (modulo the three) |
| `orderedChildren_eq_of_ne_R`, `orderedChildren_sorted` | `StSpec.lean` | proved |
| `pairwise_dominance_of_sorted_sum`, `edgeChildren_dominance` | `StSpec.lean` | proved |
| `layoutNode_r_bracket` | `StSpec.lean` | sorry (mechanical) |
| `relabel_st` | `StSpec.lean` | sorry (7.3) |
| `TEntry.wrap`, `OneSided`, `OnSide`, `nest` | `StWalk.lean` | def |
| `pushTstack_onSide`, `merge_onSide`, `fold_onSide`, `getSide_setSides` | `StWalk.lean` | proved |
| `WalkState.StSides`, `StEntry`, `StInv`, `TopClosable`, `entryVertList`, `entryEdges` | `StWalk.lean` | def |
| `chain_stackDir_const` (`ear_uniform_side`, semantic half) | `StWalk.lean` | sorry |
| `finishTstackTop_stItem` | `StWalk.lean` | sorry (mechanical) |
| `finishEdge_topClosable`, `finishEdge_stInv` | `StWalk.lean` | sorry (hard) |
| `walk_st` | `StSpec.lean` | sorry (from `finishEdge_stInv`) |

## 8. Planarity

The planar variant (`planar_spqr_tree`, `with_planarity`) runs the same walk and relabel with extra
state; the Lean port (`PlanarWalk.lean`, `PlanarRelabel.lean`, `PlanarEmbed.lean`) is a
conservative extension of the ordinary one, and the pure specification is `Planar.lean`.

### 8.1 Quarter-edges and rotation systems

A quarter-edge is `q = 4 e + 2 side + dir` **[def `QE`]**: `side` picks an endpoint of edge `e`
(`side = 0` is `edges[e].1`), `dir` one of the two corners at that endpoint. `q ^^^ 1` is the other
corner at the same endpoint, `q ^^^ 3` the corner across the edge. A `RotationSystem` **[def]** is
`rotAdj : Array (Option Nat)`; `rotAdj[q]` is the corner facing `q` around their common vertex.
Well-formedness **[def `IsEmbedding`]**: total on the `4 |E|` quarter-edges, an involution, pairing
corners at the same vertex with opposite `dir`, and with exactly `2 · #(non-isolated vertices)`
orbits of `q ↦ rotAdj[q ^^^ 1]` (every vertex gives two orbits, one per `dir`, so the pairing is a
genuine cyclic order at each vertex). Faces are the orbits of `q ↦ rotAdj[q ^^^ 3]`, again counted
twice. `IsPlanarEmbedding` **[def]** adds Euler's formula per component,
`#faces + 2 V = 2 (2 C + E)` with `C` the number of components with an edge and `V` the non-isolated
vertices; `Planar es n` **[def]** is the existence of such a rotation system. All of these are
decidable, which is what `CheckPlanarLean.lean` executes.

### 8.2 What a planar tstack entry is (Invariant P)

Fix a tstack entry `t` with edge set `E(t)` (§4.2), bottom terminal `vStart` and top `anc topDepth`.
Think of the whole *upper* DFS stack — every ancestor at depth `≤ topDepth`, i.e. every vertex a
back edge of `E(t)` can return to — contracted into a single vertex `T`. The *piece* of `t` is
`E(t)` with `T` as its top terminal: an interior spine along the ear (the tree path from `vStart`
down), the sub-ears already merged into it, and the back edges from the spine to `T`.

`tstack_planarity_t` **[def `Planarity`]** stores, per side (left / right of the spine), the
exposed ends of the piece's two outer-face boundary walks:

* `bot_ends` **[`PlSide.bot`]**: the outer and inner exposed quarter-edge of the walk along the
  tree part of the boundary, attached to the bottom-most / top-most vertex of the spine;
* `top_ends`, `top_depths` **[`PlSide.top`]**: the outer-most and inner-most exposed back edge
  leaving the piece on that side towards `T`, and the depths they return to; the back edges
  in between are chained through `quarter_edge_matches` (**[`Qem`]**, `link`) — this is the
  "linked list of the far ends of the back edges on that side". The ends are stored at `T`
  too: `rotAdj` entries at quarter-edges of back edges whose top is an ancestor are filled in
  while the ancestor is still on the DFS stack.

Side `0` holds a minimal return (`sides[0].top_depths[0] = topDepth` whenever any back edge is
open), and within a side depths increase inwards (`top_depths[0] ≤ top_depths[1]`).

**Invariant P (per entry).** Whenever `t` is on the tstack with `planarity = some p`:

1. *(embedded)* the piece of `t` has a planar embedding `ρ_t` (a rotation system on its
   quarter-edges, with `T` contracted) in which both terminals `vStart` and `T` lie on the outer
   face, and `qem` restricted to the quarter-edges of `E(t)` agrees with `ρ_t` on every pair of
   corners that are *not* exposed;
2. *(boundary)* the two outer-face boundary walks of `ρ_t`, split at the terminals, are exactly
   described by the two sides: on each side the walk runs `bot_ends[0]` (at `vStart`) … along the
   spine and the merged sub-ears … `bot_ends[1]`, then along the chain of back-edge ends
   `top_ends[0]` … `top_ends[1]` at `T`, the ends still unmatched in `qem` being precisely the four
   `bot_ends` / `top_ends` entries; the back edges on a side are nested, the outer one returning
   shallower (`top_depths` increasing inwards);
3. the span items of each side (`spans`, with their flip bits `PlEntry.flips`) are the pieces lying
   along that boundary walk, in order, and a flipped item has its own two sides swapped.

`planarity = none` **[`PlEntry.pl = none`]** records that the piece (with its terminals forced onto
the outer face) is nonplanar.

### 8.3 Why the merge test is exactly the obstruction

`merge_tstack_tops` **[`mergeTstackTops`]** glues the top entry `b` (the deeper part of the ear,
pushed later) under the entry `a` below it along their shared terminal (`a`'s bottom spine vertex =
`b`'s top, both terminals of the union being `a.vStart` and `T`). Embedding the union means
choosing, on each side, which boundary walk `b`'s side goes on; `flip_tstack_planarity`
**[`flipEntry`, `flipBeforeMerge`]** has already made that choice consistently with the st-side
the spans were assigned (`b`'s minimal return goes on side `0`). `merge_planarity` then joins side
by side **[`mergeSide`]**: `a.bot_ends[1]` is linked to `b.bot_ends[0]` (the spine continues), and
`a.top_ends[1]` to `b.top_ends[0]` (the chains of back-edge ends concatenate), unless
`a.top_depths[1] > b.top_depths[0]`.

Why this is the obstruction: two back edges `x → anc p` and `y → anc q` drawn on the same side of
the spine, with `y` above `x`, are two arcs from the spine to `T` that must *nest*: the outer arc
starts lower (`x`) and, once `T` is expanded back into the ancestor path, ends shallower
(`p ≤ q`); otherwise their endpoints interleave along the cycle spine + ancestor path and the arcs
cross. `b`'s back edges start lower than all of `a`'s and are placed inside `a`'s chain, so they
must end deeper than all of `a`'s: `b.top_depths[0] ≥ a.top_depths[1]`. If `a`'s innermost return
is deeper than `b`'s outermost, that pair crosses.
Since the flip already put `b`'s shallowest return on the side where the deepest open return is
side `0`'s minimum, failure on either side means both sides are blocked by a previous back edge
returning deeper — the obstruction Andrew's description names — and `a.planarity := none`. This is
the LR-planarity conflict for same-side return edges, specialised to the two pieces being merged.

The other planarity steps keep Invariant P:

* `make_edge_planarity` **[`makeEdgePlanarity`]**: a single virtual edge is a planar piece; a
  tree edge exposes both of its corners on each side (`bot_ends`), a back edge exposes its lower
  corners as `bot_ends` and its upper corners as the (one-element) side-`0` chain of returns.
* `finish_tstack_top` / `node_planarity` **[`finishMatches`]**: closing an S/P/R item replaces the
  piece by its cap edge; the four exposed ends of the piece are recorded as the cap's matches, so
  the cap is a planar piece with the same boundary. `maybe_unwrap_nxt` **[`unwrapPlanarity`]** is
  the inverse (the recorded matches become the ends again).
* closing the back edges when the walk returns to their target depth **[`closeSide`,
  `pruneSide`]**: a return to the current vertex `v` is no longer a back edge to `T` but part of
  the spine at `v`; its ends move from `top_ends` to `bot_ends` (`link bot_ends[1] top_ends[1]`,
  pop the chain, `edge_top_depths` recovers the next return depth);
* leaving a child **[`foldPlanarity`]**: when the ear is left through its base, side `1` holds
  only returns to `lowval`, and is folded around onto side `0`, so the finished sub-ear is a piece
  with its two terminals `lowval`-ancestor and base on the outer face.

### 8.4 Why each S/P/R local embedding is planar

Relabel **[`planarRelabel`, `layoutRot`, `setupNode`, `applyFlips`]** maps the recorded matches
into `ne_embedding.rot_adj` over the node-edge quarter-edges `4 ne + 2 side + dir`.

* **S**: `layoutRot .S` is the cycle `nvSt, nvSt+1, …` closed by the cap, embedded as a polygon:
  two faces, `2 V` vertex orbits, `2 · (2 - V + V) = 4` face orbits.
* **P**: the bond of `k` parallel edges in the order the children were listed; `k` faces,
  Euler gives `2 · (2 - 2 + k) = 2k` face orbits.
* **R**: the rotation is `node_planarity`'s matches at the moment the R item was finished, mapped
  through `mapRot` (the cap's quarter-edges renumbered to the node edge, children's virtual edges
  flipped as `applyFlips` says). By Invariant P at that moment the piece was embedded with its
  terminals on the outer face, and the cap closes the outer face, so the restriction is a planar
  embedding of the skeleton. Nonplanar R nodes get `node_planar = false` and all `rot_adj` entries
  unset; S and P nodes are always planar.

### 8.5 Why gluing through twins preserves planarity

`planar_embed` **[`planarEmbed`, `embedItem`]** processes items bottom-up over the preorder.
Every item reports the exposed ends of its cap (`outerE`); an S/P/R node links its children's
exposed ends according to its local rotation, through the twin of each non-cap node-edge, i.e. a
2-sum along the virtual edge: the child's embedding (with its cap on the outer face — the cap is
an edge of the skeleton, so it lies on two faces, and we can choose either) is inserted into the
face of the parent's embedding on the side of the virtual edge. The 2-sum of two planar
embeddings along an edge is planar: the two faces incident to the virtual edge are merged into
one, vertices `V_1 + V_2 - 2`, edges `E_1 + E_2 - 2`, faces `F_1 + F_2 - 2`, and Euler's formula is
preserved. Vertex items splice the blocks hanging off a vertex into its rotation (one new face
merge per block, each component is embedded separately, which is the per-component Euler
formula in `IsPlanarEmbedding`). If some node is nonplanar the result is `none`
(`planarEmbed_isSome_iff`, proved).

### 8.6 Lean plan (planarity)

| statement | file | status |
|---|---|---|
| quarter-edges, `RotationSystem`, `IsEmbedding`, `IsPlanarEmbedding`, `Planar` | `Planar.lean` | def |
| planar walk (`planarWalk`), relabel (`planarSpqrTree`), gluing (`planarEmbed`) | `PlanarWalk.lean`, `PlanarRelabel.lean`, `PlanarEmbed.lean` | def |
| `planarWalk_base`, `planarWalk_proj` (planar walk = ordinary walk + aux) | `PlanarWalkProj.lean` | **proved** (`propext`, `Quot.sound`) |
| `planarRelabelTree_base`, `planarRelabel_proj` | `PlanarRelabelProj.lean` | **proved** (`propext`, `Quot.sound`) |
| `planarEmbed_isSome_iff` | `PlanarSpec.lean` | **proved** (`propext`, `Classical.choice`, `Quot.sound`) |
| Invariant P (§8.2) as a Lean predicate on `PlanarWalkState` | — | to state |
| `nodePlanar_sound` (S, P cases: direct from `layoutRot`; R case: Invariant P at finish) | `PlanarSpec.lean` | sorry |
| `nodePlanar_complete` (Kuratowski-style certificate from the §8.3 crossing) | `PlanarSpec.lean` | sorry, hard |
| `planarEmbed_sound` (2-sum along twins, per component) | `PlanarSpec.lean` | sorry |
| `spqrTree_planar` (`→` from `planarEmbed_sound`; `←` needs completeness + skeletons are minors of `g`) | `PlanarSpec.lean` | sorry |

Admitted cases, precisely: `nodePlanar_sound` is admitted for all three node types (the S and P
cases are routine counting over `layoutRot`; the R case needs Invariant P); `nodePlanar_complete`
entirely; `planarEmbed_sound` entirely (its content is the 2-sum lemma plus the `V`/`Q`/`F` item
splicing); `spqrTree_planar` entirely. Everything else in this section is proved.

Work packages:
* **Invariant P**: state §8.2 on `PlanarWalkState` (per `plStack` entry, over `qem`), prove it
  for `makeEdgePlanarity`, `mergeSide`, `closeSide`/`pruneSide`, `foldPlanarity`,
  `finishMatches`/`unwrapPlanarity`; the merge-side crossing argument gives the `none` case.
* **Local embeddings**: `nodePlanar_sound` S and P by computation on `layoutRot`; R from
  Invariant P via `planarRelabel`'s `mapRot`.
* **Gluing**: a 2-sum lemma on `IsPlanarEmbedding` (face/vertex orbit counting under `link`
  of two exposed corners) and the bottom-up induction over `embedItem`.
* **Completeness**: the crossing of §8.3 as a `K₅`/`K₃,₃` subdivision — the hard, optional one.
