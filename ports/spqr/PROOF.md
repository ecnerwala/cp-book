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

**Fact C (type-2 pairs live on first-child chains) [lemma, graph theory].** If `{a, b}` is a
type-2 pair then `b` is reached from `a` by one tree edge to some child `a'` followed by a chain
of **first** children in sorted order: for each `x` on `P ∪ {b}` with parent `p ≠ a`, `p → x` is
the first `ret` edge of `p`.
*Proof.* Let `p ∈ P` have a child `x` toward `b` and suppose an earlier `ret` edge `p → y` (or
back edge `p → anc`) exists. Then `lowpt1 y ≤ lowpt1 x ≤ l` (`T_x ⊇ T_b` reaches `a`), so `T_y`
(or the back edge) reaches depth `≤ l`; it is part of *between* yet attaches to `above ∪ {a}`. If
it reaches `< l`, between and above are joined: contradiction. If it reaches exactly `l` only…
then `p → y` is type-1 for `p` with the same `l` — it is a class of its own hanging at `{p, a}`,
and `p` is interior to the between piece; the pair `{a,b}` still separates. *(So the precise
statement: non-first edges on the chain are allowed only if they are type-1 edges returning
exactly to `a`; these are "virtual back edges" `p — a` and get glued to the between piece as
leaves. The algorithm handles exactly this via `firstOccurrence`/`firstIdx`, see 4.4.)*

Similarly, for `a` itself: `a'` need not be `a`'s first child — a type-2 pair is attached to `a`
from whichever child.

**Fact D (laminarity) [lemma, graph theory; to be stated precisely].** Order the edges of `B` by
*postorder* σ: `dfsVisit`'s finishing order — a vertex's out-edges in sorted order, each tree edge
placed right after its subtree's edges. (This is the order in which `finishEdge` is called, and
`nxtEdgeIdx`/`firstIdx` count the back edges in it.) Then every separation class is a contiguous
interval of σ **up to S/P multiplicity**: the classes of a type-1 pair are intervals; the classes
`between ∪ suffix` of a type-2 pair are an interval; and the family of all such intervals is
laminar. The S/P caveat: a cycle (S) or bond (P) can be split at any point, so the maximal pieces
of the canonical decomposition, not the individual classes, are what is unique. This is the
statement that makes a *stack* the right data structure: the open intervals at any time are nested.

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
| 1.1, 1.2 DFS spanning + lowpoints | `Correctness.lean` (`dfsForest_spanning`) | sorry |
| 2.1 blocks ↔ `lowval ≥ d` branches | — | to state |
| Facts A–D | — (pure graph theory over the sorted DFS tree) | to state |
| ear-structured walk (`descend`/`ascend` over chain `Frame`s) | `Ear.lean` | def |
| `walkEarTree = walkTree`, frame rule `walkTree_local` | `EarSpec.lean` | sorry |
| Invariant W, Lemmas 4.3/4.4 (`earOut_one_entry`, `ascend_frame_one_entry`) | `EarSpec.lean` | sorry / hard |
| 4.5 maximality / R 3-connected | `spqrTree_r_three_connected` | hard |
| 5 relabel: `Items.WF → WF ∧ Represents` | `relabelTree_wf`, `relabelTree_represents` | sorry |

Work packages for child sessions, in dependency order:
* **DFS**: 1.1, 1.2, no cross edges, `lowpt` characterization of `OutClass`.
* **Relabel**: 5, independent of the walk (assumes `Items.WF`).
* **Ear reformulation**: the iterative `walkEar` intermediate and its equality to `walkTree`.
* **Graph theory**: Facts A–D on an abstract sorted DFS tree.
* **Walk invariant**: Invariant W + Lemma 4.3 on `walkEar`, giving `Items.Tree`, `Endpoints`.
* **Shapes / maximality**: `Items.Shapes`, R 3-connectivity via 4.5.
