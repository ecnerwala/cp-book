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
* *Nesting, type-2 vs type-2* (`type2Block_laminar_type2Block`): the type-2 classes of two
  pairs `{a₁, b₁}`, `{a₂, b₂}` (each under the hypotheses of `type2_class_interval`) are nested
  or disjoint **unless** one pair's `a'` lies strictly below the other's `a'` on the other's tree
  path `a' ⇝ b` (inclusive of `b`) — again the S caveat, and the exception is necessary: in
  `r → a₁ → a₂ → b₁ → c → b₂` with back edges `b₂ → r`, `b₂ → b₁` (2-connected), `{a₁, b₁}` and
  `{a₂, b₂}` are both type-2 pairs (`a₂' = b₁`) cutting the cycle `r a₁ a₂ b₁ b₂`, with classes
  `{a₁–a₂, a₂–b₁}` and `{a₂–b₁, b₁–c, c–b₂, b₂–b₁}` that overlap without nesting. The proof
  reduces to the structural cases of `type2Block_laminar_block` when `a₁' ≠ a₂'`; for `a₁' = a₂'`
  (so `a₁ = a₂`) it uses §1: `b₁, b₂` are comparable (otherwise `T_{b₁}` returns above `a` inside
  the *between* part of `{a, b₂}`, `subtree_returns_above`), and if `b₂` is below `b₁` the child
  of `b₁` containing `b₂` has `lowpt1 < l`, so it is sorted before `ret l backEdge` and the class
  of `{a, b₁}` is contained in that of `{a, b₂}` (`type2Block_subset_type2Block`).

What §4.5 needs is only this weaker form: two classes that are *not* cuts of one S cycle (one
pair's path `a' ⇝ b` is not a proper continuation of the other's) are nested or disjoint, so the
candidates open at any time form a stack, and the classes of one S cycle are handled by the
cycle's own frame.

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
global tstack (`WalkSpec.Inv D`, with `D` the current depth bound: `d + 1` while a tree edge's
subtree is being walked, `d` for a back edge) has two local parts per entry `t = (vStart, topDepth, …)`
with edge set `edges(t)` (`WalkSpec.EntryInv D`):

* **connected**: `edges(t)` is connected (`Graph.ConnEdges`), and so is each item already closed;
* **attached at the open path**: every vertex incident to an edge of `edges(t)` and to an edge
  outside it is `vStart` or `stackVerts[k]` for some `topDepth ≤ k ≤ D` (`Graph.AttachedIn` at
  `TEntry.Term D`).

The second part is deliberately weaker than "2-attached at `{vStart, stackVerts[topDepth]}`":
strict 2-attachment is *false* mid-walk for a type-2 entry. After the type-2 first edge `u→c` of a
vertex `u` at depth `d` with lowpoint `l` has been finished (loop 1 closed the ears of the subtree
and loop 2 merged the late entries), the top entry `(u, l, …)` still contains the chain vertices at
depths in `(l, d)`, whose remaining out-edges have not been processed, so it is attached at
`u`, at `stackVerts[l]` *and* at those intermediate chain vertices (this is exactly the attachment
set §4.3 describes: `u`, `anc l`, and the depths in `(l, d)`). Only once the vertex close
(`closeVert`/the `!hasVert` push) has absorbed the chain does the entry become two-terminal. Closed
items, by contrast, are strictly connected and 2-attached (`WalkSpec.ItemInv`): `finishTstackTop`
closes an entry `t` only under the extra hypothesis that `edges(t)` has **no attachment at
`stackVerts[k]` for `topDepth < k ≤ D`** other than `vStart` (`FinishTopOk.mid`: every such
vertex is `vStart`, interior to `edges(t)`, or untouched), which collapses `Term D` to the two
terminals `{vStart, stackVerts[topDepth]}` (`TwoAttached.of_term`). This hypothesis is vacuous for
the loop-1 closures (`topDepth ≥ d`, the sub-ear is finished) and is the type-1 condition for the
vertex close and the P-check.

**Correction (ear session 3, `EarInv.lean`, `checks/InvCheck.lean`).** The index `D = d + 1`
"while a tree edge's subtree is being walked" is **not** enough, and no choice of `D` is: the
sub-ears that loop 1 closes are only those *above* the first entry with `topDepth < d`; an entry
of a deeper chain vertex that returns above `d` stays open underneath it. Cycle `0-1-2-3-4-5-6-0`
with chords `6-1`, `5-2` and the pendant ear `4-7-3`: after `finishEdge` of the type-2 frame
`(4, 4)` the stack is `(4,4) (5,4) (5,2) (5,5) (6,5)=[Q(5,6)] (6,1) (6,0) (6,6)`; `(6,5)` is attached
at `5`, so `Inv 4` fails (`WalkInv.ear_lower` is false and `walkTree_inv'`'s `Inv d` post-condition
with it), and as soon as the sibling `4→7` is entered (`stackVerts[5] := 7`) vertex `5` is
`stackVerts[k]` for no `k`, so `Inv D` fails for every `D`. The attachment set has to include the
bottoms of the entries above: `TEntry.Term' D s above t v := Term D s t v ∨ ∃ t' ∈ above, v =
t'.vStart` (`EarInv.EntryInv'`/`Inv'`, implied by `Inv`). With that clause `D = d` (the current
depth, under tree edges too) is enough — 0 violations on 9000 random multigraphs (≤ 12 vertices,
≤ 24 edges) at every `walkTree` start/end and after every `walkOut`, while `Term d` fails 1827
times and `Term (d+1)` 728 times on the same snapshots; connectivity of every open entry also had
0 violations. The per-block lemmas of `WalkSpec.lean` are unaffected (they are stated for `Inv D`
with `D` free and transfer to `Inv'`, whose extra clause is preserved by pushes/merges that only add
entries above or merge adjacent ones); the walk induction `WalkInv.walkTree_inv'` has to be
restated for `Inv' d`.

*Soundness*: every `mergeTstackTops` joins two entries sharing a terminal, so the union is again
connected and attached inside the union of the two `Term` sets, which the `min topDepth` rule
re-expresses as `Term` of the merged entry (`mergeTstackTops_sound`, hypotheses `MergeOk`).
*Completeness*: every item that `finishEdge` closes (S/P/R, or the single entry handed up as a
virtual edge) is a 2-attached connected edge set, i.e. cut off by a genuine 2-vertex separation of
the block (`finishTstackTop_complete`), and the closed items partition the edges. This gives
`Items.Tree`, `Items.Endpoints` and the S/P/Q/I/O shapes without ears. `finishEdge_inv` proves
the step for all returning out-edges (type-1 child with/without the vertex entry, type-2 child
with its three loops, back edge) under the per-block stack-shape hypotheses `FinishOk` (which
entries are on top, their terminals, `topDepth` relations, the one-sided spans the C++ asserts).
The ear view (§4.1–4.3) is used only to discharge those (`walkTree_guards`/`EarShape`: `size ≥
origTstack + 3`, `nxt` is the chain/vertex entry, the merge loops never cross a frame) and for
maximality (§4.5), where Facts C–D are needed.

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

**Lean target (`WalkCover.lean`).** The "nothing is dropped" half of `Items.Tree` is proved from a
single admitted predicate, `walk_sides : SidesForest forest (WalkState.init g tern)` (hypotheses:
`ForestOK`, `∀ t ∈ forest, t.WF []`). `SidesForest`/`SidesTree`/`SidesOuts`/`SidesOut` mirror the
walk in the style of `GuardsTree` and assert, at each discard site, one of:

- `MergeOK s` (`2 ≤ s.tstack.length`) before every `mergeTstackTops`, and the `origTstack + 3`
  shape of the block branch — these are exactly `EarSpec.walkTree_guards` (`FinishGuards`/
  `GuardsTree`: the popped entries exist);
- `CloseOK s` before every `finishTstackTop`: the top entry `t` satisfies
  `t.OnSide s.stackDir[t.topDepth]!`, i.e. the side `getSide t.spans (!topDir)` it discards is `[]`;
- `UnwrapOK s ty` before `maybeUnwrapNxt ty` (when it may reuse): if the kept side of the entry
  below the top starts with an item of type `ty`, that entry is `[x]` on its own side;
- `BoundaryOK curV d o s` in the block branch: `spans.1 = []` of the top (type-2 close, lowval
  `= d+1`) resp. `spans.2 = []` of the back-edge entry and `spans.1 = []` of the vertex entry below
  it; and on each kept side no item `c` has `Items.Below s.items c (vertItem curV)` (the Q item
  takes the kept sides as its children and is then put under `vertItem curV`; this clause is what
  keeps `ch` acyclic, see below);
- `RootOK s` after `walkTree t 0`: one entry is left, `spans.1 = []`, and `spans.2` consists of
  `vertItem`s of real vertices (it is `[vertItem t.v]`).

Which ear fact discharges each: an entry closed by `finishTstackTop` is a whole ear (or a type-2
piece of one, Fact C), so by the paragraph above all its pieces were attached on the side
`stackDir[topDepth]` — `StWalk.chain_stackDir_const` (one lowval along the ear ⇒ one `stackDir`)
gives `TEntry.OnSide t stackDir[t.topDepth]!`, which is literally `CloseOK`; `setSides_onSide`,
`merge_onSide`, `fold_onSide` show `pushVertTstack`/`pushEdgeTstack`/`mergeTstackTops`/the final
fold preserve `OnSide` for a fixed `dir`, so the per-ear `EarShape` invariant of the ear session
should carry `OnSide (stackDir[topDepth])` for every entry above the ear's boundary.
`BoundaryOK` is the same statement at the ear boundary, in the orientation `finishEdge` fixes
there (`setStackDir d false` for a tree child with lowval `≥ d`, so the vertex ear's entry has
`spans.1 = []`); `UnwrapOK` is the sub-ear case of `CloseOK` (a finished sub-ear sits alone in
its entry, on the side it was finished on); `RootOK` is `BoundaryOK` at depth 0 after the root's
`pushVertTstack` with `stackDir[0] = false`. The extra `¬ Below c (vertItem curV)` clause of
`BoundaryOK` is the ear fact that the popped spans, and everything already hanging below them,
were built while walking the subtree under `o`, which never places `vertItem curV`.

**Acyclicity (`ItemAcyc.lean`, `WalkState.Full.acyc`).** Exact placement plus coverage say every
non-root item has exactly one parent but not that `Items.IsParent` is well-founded, so `Full`
also carries

```
Items.NoParent items r := ∀ p, ¬ items.IsParent p r
Items.Acyc items      := ∀ i < items.size, ∃ r, items.NoParent r ∧ items.Below r i
```

(a forest: every item sits below a parentless one). The tempting "parent id > child id" is false
for the walk: in the block branch the Q item `edgeItem o.e` (allocated at the start, so older than
every node item) receives the popped spans as children and then goes under `vertItem curV`.
`Acyc` is preserved by `Acyc.push` (a fresh childless item), `Acyc.append` (append `L` to `ch a`
when no `c ∈ L` is above `a` — `pushEdgeTstack`'s self-loop O item, `root_append`, and the Q
re-wiring `acyc_q_write`, which is where `BoundaryOK`'s new clause is used), `Acyc.dropCh`
(clearing a `ch` list whose members have no other parent, by `Place.le` — the `maybeUnwrapNxt`
reuse), and `Acyc.of_ch_eq` (everything that leaves `ch` alone: `finishTstackTop` only moves items
from spans to a fresh node, `mergeTstackTops` only moves between spans). `Acyc.reach` then gives
`walk_reach`: by `walk_covered` every item `0 < i < size` has a parent, so the parentless root
that `Acyc` provides is `rootItem`.

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

**R-maximality, formal statement (`RMax.lean`, `Proofs/RMax.lean`).** An R close of `finishEdge`
(Loop 1, entries with equal `topDepth` and distinct `vStart`s) is described by the walk-free
predicate `RCloseShape g d P U s t`: `U` is the union of the closed entries' edges, `s` the top
entry's `vStart`, `t = stackVerts[topDepth]`, and `P : Pieces` the merged entries' items (the
maximal sub-pieces of `U`; `P.addParent g U s t` adds the complement of `U` as the parent piece with
terminals `s, t`, so the R skeleton is `(P.addParent g U s t).contract g`). Its fields:

* graph structure of the close (statement 3 of Invariant W): `wf : P.WF g`, `sub` (sub-pieces are
  inside `U`), `conn`, `attached : TwoAttached U s t`, `touch_s`, `touch_t`, `ne : s ≠ t`, `proper`;
* `single`: `U` is one separation class of `{s, t}` (Loop 1 chose R, not P);
* `maximal`: each sub-piece's complement is one class of its terminal pair (the sub-piece was closed
  as a maximal item);
* `bond`: two parallel edges of `U` lie in a common sub-piece (the P-check merged them);
* `type1`: for a skeleton pair `{a, b}` of `U` (`SkelPair`: vertices of `U`, interior to no
  sub-piece, not a sub-piece's terminal pair, not `{s, t}`) with `a` an ancestor of `b`, the class
  `T_c ∪ {b–c}` of a type-1 child `c` of `b` returning to `depth a` is `LaminarWith U`: inside one
  sub-piece, disjoint from `U`, or containing `U` — the type-1 close / P-check at `b` ran before `U`
  formed;
* `type2`: for a skeleton pair `{a, b}` that is a `Type2Pair`, the class of the tree edge `a → a'`
  towards `b` is `LaminarWith U` — the `firstIdx > firstOccurrence[d]` test fired a
  `finishTstackTop` before `U` formed.

Proved from these (pure graph theory + §3 Facts, `Proofs/RMax.lean`): `RCloseShape.not_sepPair`
— no `SkelPair` of `U` is a separation pair of the block (`sepPair_iff'` reduces to a type-1 /
type-2 pair; its class is 2-attached at `{a, b}`, touches both and is proper, and
`Pieces.not_laminarWith` shows no such class can be laminar with `U`: inside a sub-piece makes
`{a, b}` that piece's terminals, disjoint from `U` makes it `{s, t}`, containing `U` contradicts
`exists_not_mem_end`) — and `RCloseShape.threeConnected :
((P.addParent g U s t).contract g).ThreeConnected` via `threeConnected_contract_iff`, with the
terminal pairs handled by `Pieces.WF.not_sepPair_terminal` from `maximal` / `single`
(`addParent_wf` uses `TwoAttached.conn_compl`: the complement of a 2-attached set in a block is
connected).

Walk side (`RClose.lean`, `Proofs/RClose.lean`). At an R step of Loop 1 (`loop1Type` returns
`.R`: `tstack = cur :: nxt :: rest`, `cur.topDepth = nxt.topDepth = d`, `nxt.vStart ≠ cur.vStart`)
the closed set is `U = rU = edges cur ∪ edges nxt`, the terminals are `nxt.vStart` and
`stackVerts[d]`, and the sub-pieces are `Pieces.ofItems` of the merged non-`V` items
(`rPieceItems`). `WalkState.RStep` records this shape plus what the ear structure says about it:
`cur.vStart` is interior to `U`, both entries are nonempty, `U` is proper, and the merged items
are pairwise edge-disjoint 2-terminal pieces (`PieceItems`). The invariant hypothesis is not
`Inv d` (at the R branch `cur` is attached at `stackVerts[d+1] = cur.vStart`, so `Inv d` is
contradictory there) but `EarInv.Inv' D` (§4.2b) with `stackVerts[k] = cur.vStart` for `d < k ≤ D`;
at the branch `D = d+1`. `cur` is the top entry, so its `Term'` is `Term`; `nxt`'s extra `Term'`
vertex is `cur.vStart`, interior to `U`, so `union_attached` needs no further shape fact. Proved:
`RStep.rCloseShape' : Inv' D → (∀ k, d < k → k ≤ D → stackVerts[k] = cur.vStart)
→ TwoConnected → RStep → RContent → RCloseShape …` — the structural fields `wf`, `sub`, `conn`,
`attached`, `touch_s`, `touch_t`, `ne`, `proper` come from `Inv' D` (`EntryInv'` via
`TwoAttached.of_term`, `ConnEdges.union`/`TwoAttached.union`, `twoAttached_union_classes`) and
`PieceItems` (`ofItems_mem_iff`) — and `RStep.threeConnected'`.

Content from the history of the walk (`RInv.lean`, `Proofs/RInv.lean`). What the five content
fields of `RContent` and the non-`Inv` fields of `RStep` need is a per-entry predicate
`EntryR dfs t`: the closed non-`V` items of `t` are pairwise edge-disjoint 2-terminal pieces
(`pieces`), each one's complement is a single separation class of its terminals (`maximal`: every
P-check at its bottom vertex has fired), parallel edges of `t` lie in one closed item (`bond`),
if `t` is not 2-attached at `{vStart, stackVerts[topDepth]}` its edges are one class of that pair
(`single`), and the type-1 / type-2 classes at skeleton pairs of `t` are laminar with `t`
(`type1`, `type2`: every type-1 close below and every `firstIdx > firstOccurrence` merge at or
below the entry has fired). `RTop dfs cur nxt` is `EntryR` of the two top entries plus their
edge-disjointness; `RBranch d cur nxt rest` the branch shape (`loop1Type = .R`, `cur.vStart =
stackVerts[d+1]`, `cur` has a closed item, `nxt` touches both terminals and has no `cur.vStart`–
`stackVerts[d]` edge). Proved: `RBranch.rStep : RTop → RBranch → RStep`, `RBranch.rContent :
Inv' (d+1) → TwoConnected → dfs.Spec → RTop → RBranch → RContent` (all five fields), and
`RBranch.threeConnected` (the contracted R skeleton is 3-connected) from `Inv' (d+1) ∧ RTop ∧
RBranch` alone. `RInv dfs s` (all entries `EntryR`, pairwise disjoint) is not an invariant of every
intermediate state: the single-edge entry a back edge `v → stackVerts[l]` pushes is not `maximal`
until the P-check with the next `(v, l)` class, and the tree-edge entry of `pushEdgeTstack` is not
either. The settled form is `RInvAt dfs v` (entries with bottom `v` exempt). `Proofs/RInvFrame.lean`
proves `EntryR`/`RInvAt` depend on the state only through `g`, `stackVerts`, the tstack, and the
types/terminals/edge sets of the entries' items (`EntryR.congr`, `RInvAt.congr`), and transports
them across the bookkeeping steps of `finishEdge` (`modifyItem` of a free item, `setStackDir`,
`modify` of `firstOccurrence`, `pushVertTstack`: `EntryR.vert`). Admitted (named):
`finishEdge_rInvAt` (the content blocks of `finishEdge` keep `RInvAt curV`), `walkTree_rInvAt`
(finishing a child settles its entries), `loop1_rBranch` (every `.R` iterate of Loop 1 is
`RBranch ∧ RTop`); `loop1_r_threeConnected` combines the last with `RBranch.threeConnected`
(`closeEars_iter_step` gives `Inv (d+1)` at the iterate via `Step`; to be restated for `Inv'`).
`loop1_rBranch`'s hypotheses (`Inv D`, `Shape`, `CloseEarsOk`, `RInvAt`) say nothing about which
edges the Loop-1 entries hold, so `interior`, `proper`, `nxt_ne` and `cur_c` (after an S merge the
head's `vStart` is the S entry's) are not derivable from them: its proof needs the Loop-1 ear
content (entries with `topDepth > d` at `closeEars` have `vStart = nxtV`; the entries with
`topDepth ≥ d` hold exactly the tree edge and the child's subtree edges) as an extra hypothesis
or from `EarShape`. The shape fields (`tstack`, `cur_top`, `nxt_top`, `ne`) follow from
`run_loop1Cond`, `loop1Type_run` and the head-`topDepth` induction over the iterates.

Item level (`Proofs/RItems.lean`). `Items.RSkel3 g items i`: the skeleton of the R item `i` (its
non-`V` children as pieces, the complement of its edge set as the parent piece at `i`'s terminals,
contracted) is 3-connected. `RBranch.rSkel3` proves it for the item Loop 1's R branch closes
(`rCloseItems`: `allocItem .R`, merge, `finishTstackTop`), from `RBranch.threeConnected`
(under `Inv' (d+1)`) and `Shape`. Admitted: `items_r_three_connected` (every R item of `g.walk` on a block is `RSkel3`;
needs `loop1_rBranch`, the same argument for the type-1 R close of `finishEdge`, and that closed
items are never modified again). For a graph with several blocks the parent piece must be the
complement within the item's block; `RSkel3` with the whole graph is the block statement only.

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

Per-node layout (`LayoutShape.lean`) **[proved]**: for every type the `Layout` that `layoutNode`
returns for one node is pinned down locally, so the transport above is a rewrite once the per-node
interface (`RelabelNode`, `RelabelSpec.lean`) identifies node `i`'s `nodeEdges`/`adjBounds`/`adjDat`
segment with `layoutNode it.type … edgeChildren`. `Layout.Shape type nvSt nvEn` is `SpqrTree.Shape`
verbatim on a single node's `edges`; `Layout.Local node nvSt nvEn neSt neEn` collects the local
forms of `WF.adj_bounds_mono`/`adj_last`/`adj_incident`/`adj_dest` and `Ownership.ne_nvs` (sizes `neEn - neSt`, `2(nvEn - nvSt) + 1`, `2(neEn - neSt)`; every edge has
`node`, `twin := none`, endpoints `nvSt ≤ a ≤ c < nvEn`; bounds monotone and ending at `2·neEn`;
row `2nv` is a permutation of the edges with `nvs.2 = nv`, row `2nv + 1` of those with
`nvs.1 = nv`; every entry's `destNv` is an endpoint of its edge). Closed forms of the rows:
`runF_row` (all bounds `2 neSt`), `runLoop_row`, `runQI_row`, `runP_row` (parallel edges; the row of
`nvSt + 1` in reverse edge order), `runS_row` (closing edge `(nvSt, nvEn-1)` at `neSt`, path edge `k` at
`neSt + 1 + k`), and for R `run_entries` on top of `StLayout.run_spec` (slot `0` and the last slot
are the cap, fill slots carry child edges in strictly decreasing `ne` order). Hypotheses are those
of `Items.Shapes` read through `vertPos`: Q/I `nvEn - nvSt ∈ {1, 2}`, `neEn - neSt = 1`; P
`nvEn - nvSt = 2`, `3 ≤ neEn - neSt` (`p_shape`: ≥ 2 virtual edges + the cap); S
`3 ≤ nvEn - nvSt = neEn - neSt` (`s_shape`); R `nvSt + 2 ≤ nvEn`, `4 ≤ nvEn - nvSt` (`r_shape`:
≥ 2 V children), `6 ≤ E.length + 1`, `neEn = neSt + E.length + 1`, every child edge
`nvSt ≤ a < c < nvEn` (`r_shape`'s `q.1 ≠ q.2` plus the orientation `edgeChildren_dominance`
gives), and `((nvSt, nvEn-1) :: E).Nodup`. The last one needs `r_shape`'s clause that no child
virtual edge is parallel to the node's own `vs` (the cap edge; the key-`Nodup` clause alone only
makes the children pairwise distinct): `r_skeleton_nodup` derives it from the two clauses through
the position map `vertPos` (injective on the node's vertices).
`PlanarSpec.neRotAdj_segment` reads no adjacency row (it is about `layoutRot`, which does not look
at the `Layout`); the layout-side facts its transport needs are the `Local` sizes and
`bound_last` (the `adjBounds.extract 1 …` pushed by `planarRelabel` has `2(nvEn - nvSt)` entries
and ends at `2·neEn`).

The structural half of this is done through the per-node interface of `RelabelSpec.lean`
(`RelabelIdx`: the preorder index `idx` is a bijection with root `0`, agrees with
`vertIndex`/`edgeIndex` on `V`/`Q` items, CSR zero/last facts; `RelabelNode`: each item's slot,
ranges, `node_verts`, `child_par`, and a pos-dependent `RelabelLayout` whose edge/adjacency rows
are `layoutNode`'s). `RelabelOwn.lean` takes `relabel_node_spec` (**[sorry]**, the glue
`∃ idx, …` for the recursive fold) as a hypothesis and proves `relabelTree_own :
Items.WF g items → Items.OwnExtra g items → Items.ROriented g items →
(relabelTree g items).Bijections ∧ .Ownership ∧ .Twins` **[proved]**:
* `Bijections` from `Items.Tree` only (types of `1+v`, `1+nv+e`, `RelabelIdx.vert_index/edge_index`,
  `RelabelNode.orig`);
* `Ownership` from `Items.Tree`, `Endpoints.vs_shape`, `Shapes.*_shape`/`i_o_leaf`/`q_children`
  (the per-type `layoutNode` hypotheses), the `layoutNode` edge records (`layoutNode_edges`, read
  off `LayoutShape.lean`'s `layoutNode_*_eq`/`run*_edges_get`/`run*_spec`; the R edges are
  re-derived there without `LayoutShape.run_edges`'s endpoint hypothesis so that `Twins` does not
  need `ROriented`), and the two extra hypotheses: `Items.ROriented`
  (`ne_nvs` for R nodes needs the ordered endpoints to lie in the node-vert list, which is the
  st-ordering fact of §7, not part of `Items.WF`) and `Items.OwnExtra` — four facts about the walk's
  items that `Items.WF` does not state: node-vert lists are nodup (`WF` admits `I`/`P` with `u = v`,
  which would make `nv_distinct` false), a `Q` child of a node is a leaf, both endpoints of an R
  node's virtual edges are node-verts (`WF` only constrains `some` endpoints), and the `vertItem v`
  child of a block-root `Q` has `v < nv`. These are to be discharged by the walk (`WalkTyping.lean`
  / `StSpec.lean`), not by relabel;
* `Twins` from `Items.Tree` and `OwnExtra.q_leaf_of_node` (a capped child has at least one edge),
  via `RelabelLayout.twin` both ways, `child_cap_twin_none`, injectivity of `idx` and disjointness
  of the `neRange`s.

## 6. Lean plan (what is proved where)

| statement | file | status |
|---|---|---|
| graph / DFS / items / walk / relabel definitions | `Graph.lean … Build.lean` | def |
| output spec `SpqrTree.WF`, `Represents` | `Spec.lean` | def |
| phase-2 contract `Items.WF` | `ItemSpec.lean` | def |
| 1.1, 1.2 DFS spanning + lowpoints | `Proofs/Dfs.lean` (`dfsForest_spanning`, `dfsForest_wf`, `classify_child_*`) | proved |
| 2.1 blocks ↔ `lowval ≥ d` branches (`sameBlock_iff`, `blockRoot_cut`, `ret_child_sameBlock`) | `Blocks.lean`, `Proofs/Blocks.lean` | proved |
| Facts A–C (`sepPair_comparable`, type-1 class / above-between / sorted prefix-suffix, `type2_first_out`) | `SepPair.lean`, `Proofs/SepPair.lean` | proved |
| Fact D (laminar intervals: `type1_class_interval`, `type2_class_interval`, `type2Block_laminar_block`, `type2Block_laminar_type2Block`) | `Proofs/Postorder.lean`, `Proofs/Interval.lean` | proved (S-caveat exceptions) |
| 4.5 exhaustiveness: separation pairs of a block = type-1 ∪ type-2 pairs (`sepPair_iff`, `three_connected_of_no_split`) | `SepPairExhaust.lean`, `Proofs/SepPairExhaust.lean` | proved |
| 5 skeleton/contraction: `Graph.contract` of a laminar family of 2-attached pieces; `sepPair_contract_lift`, `sepPair_contract_of` (non-terminal pairs), `threeConnected_contract_iff_dfs` (R skeleton 3-connected ↔ no type-1/type-2 pair among non-terminal skeleton vertices), `TwoAttached.sepClass_mem` | `Contract.lean`, `Proofs/Contract.lean` | proved |
| walk pieces ↔ separation classes: a nonempty proper `TwoAttached` set of a block is a union of `{u,v}`-classes with `{u,v}` a separation pair, or it / its complement is a single `u–v` edge (`twoAttached_union_classes`, `twoAttached_type1_or_type2`) | `Proofs/SepClasses.lean` | proved |
| ear-structured walk (`descend`/`ascend` over chain `Frame`s) | `Ear.lean` | def |
| `walkEarTree = walkTree` (`walkEarTree_eq_walkTree`, `walkEar_eq_walk`) | `EarSpec.lean` | proved |
| span discipline / placement (`Place`: each item id placed ≤ 1 time over spans + ch lists; `walk_place`, `walk_ch_nodup`, `walk_parent_unique`, `root_no_parent`) | `WalkPlace.lean` | proved; `walk_tstack_nil`, `walk_covered`, `walk_root_children`, `walk_reach` moved to `WalkCover.lean` |
| coverage + reachability half of `Items.Tree` (`WalkState.Full` = exact placement + `Acyc`; `walk_tstack_nil`, `walk_covered`, `walk_root_children`, `walk_reach`, `walk_tree`) | `WalkCover.lean`, `ItemAcyc.lean` | proved from `walk_sides : SidesForest forest (WalkState.init g tern)` (§4.4; the only `sorry` in the file, now including the `¬ Below c (vertItem curV)` clause of `BoundaryOK`) |
| linear phase 2/3 refinements `walkFast` / `relabelTreeFast` (`CatList` spans, `Array` tstack, ticks): `walkFast_items`, `relabelTreeFast_eq`; `spqrTree_eq` routed through them | `CatList.lean`, `Refine.lean`, `WalkFast.lean`, `RelabelFast.lean`, `Correctness.lean` | proved |
| step bounds: `walk_ticks_le : (g.walkFast tern (g.dfsForest vo eo)).ticks ≤ 47·(nv+ne)` (via `walk_ticks_le_forest`, `dfsForest_size_le` from `dfsForest_spanning'`); `relabelRun_sizes_le` (`Items.desc` cardinalities), `relabel_ticks_le : ticks ≤ 288·(nv+ne) + 6` | `WalkCost.lean`, `ItemTree.lean`, `RelabelCost.lean` | proved; the relabel bounds take `Items.Tree` as hypothesis (`relabel_ticks_le'` discharges it with the admitted `walk_items_wf`) |
| frame rule `walkTree_local` via `Lifts`/`Sim` simulation (`Sim.closeEars`, `Sim.mergeLate`, `Sim.finishRest`, `Sim.finishBoundary` proved) | `Sim.lean`, `Frame.lean`, `EarSpec.lean` | `Sim.closeVert`, `Sim.finishEdge`, `Sim.walkTree` proved; `walkTree_local` reduces to the stack-shape invariant `walkTree_guards` (admitted, with `earOut_one_entry` / `ascend_frame_one_entry`) |
| typing/allocation part of `Items.WF` (`Items.Tree` sizes/types, I/O leaves, `vs_shape`, `vs_lt`): `walk_typing` | `WalkTyping.lean` | proved (`walk_q_children` sorry: needs span shape) |
| §4.2b walk invariant `Inv D` (`EntryInv D`: connected + attached at `vStart`/`stackVerts[topDepth..D]`; closed items 2-attached): closure lemmas (`GraphLemmas.lean`: `AttachedIn`, `twoAttached_iff`), `mergeTstackTops_sound`, `finishTstackTop_complete`, `Shape`/`Step` infrastructure, per-block lemmas `Step.closeEars`/`mergeLate`/`closeVert'`/`finishRest` under `CloseEarsOk`/`MergeLateOk`/`CloseVertOk`/`FinishRestOk` | `GraphLemmas.lean`, `WalkSpec.lean` | proved |
| §4.2b `finishEdge_inv` (type-1 ± vertex entry, type-2 three loops, back edge) under `FinishOk`; `finishEdge_back_inv` corollary | `WalkSpec.lean` | proved; `FinishOk ← FinishGuards`/`EarShape`, `walkTree_inv` sorry — and **not provable as stated**: `ear_lower` is false and `Inv D` fails for every `D` under a type-2 chain with a sibling subtree (§4.2b correction, `EarInv.lean`); restate for `EarInv.Inv' d`; `walk_nodes_partition` proved in `WalkPlace.lean` under `ForestOK` + coverage |
| Invariant W, Lemma 4.3 (`earOut_one_entry`) | `EarSpec.lean` | sorry / hard |
| Lemma 4.4 (`ascend_frame_one_entry`: a finished frame's vertex owns one entry) | `EarSpec.lean` | **false** for chain frames (cycle `0..5` + chord `5-1`, frame `(4,4)`: five entries); removed, the collapse holds only at the ear's top (= `earOut_one_entry`) |
| boundary branch of `finishEdge` keeps `Inv D ∧ Shape` (`finishBoundary_inv`, via `BStep`, under `BoundaryOk`: popped entries exist, Q/V items are roots not on any span) | `WalkInv.lean` | proved (`BoundaryOk ← ear_boundary` sorry); the former `Step` form is false (the Q item goes under `vertItem curV`), as is `VertBook`'s `hasVert = false → ch (vertItem v) = []` (bridge `1-2` before back edge `1→0`) — replaced by connectivity + `TwoAttached v v` of the vertex item |
| corrected attachment set `TEntry.Term'`, `EntryInv'`, `Inv'` (`Inv'.of_inv`, `Inv'.mono`); empirical check `checks/InvCheck.lean` | `EarInv.lean` | def + proved; the walk induction for `Inv'` is open |
| 4.5 maximality: `RCloseShape` ⇒ no skeleton pair separates (`RCloseShape.not_sepPair`), R skeleton 3-connected (`RCloseShape.threeConnected`) | `RMax.lean`, `Proofs/RMax.lean` | proved; `RStep.rCloseShape'`/`RStep.threeConnected'` (`Proofs/RClose.lean`) give it for Loop 1's R step from `EarInv.Inv' D` + `stackVerts[d+1..D] = cur.vStart` + `RStep` + `RContent` (`Inv d` is contradictory there) |
| 4.5 walk side: `EntryR`/`RTop`/`RBranch`/`RInvAt`; `RBranch.rStep`, `RBranch.rContent` (all five content fields), `RBranch.threeConnected` (from `Inv' (d+1)`); `EntryR.congr`/`RInvAt.congr` + bookkeeping frames; `Items.RSkel3`, `RBranch.rSkel3` | `RInv.lean`, `Proofs/RInv.lean`, `Proofs/RInvFrame.lean`, `Proofs/RItems.lean` | proved; admitted: `finishEdge_rInvAt`, `walkTree_rInvAt`, `loop1_rBranch` (history preservation), `items_r_three_connected` (all R items of the walk on a block); `spqrTree_r_three_connected` (relabel transport): hard |
| 5 relabel: `Items.WF → WF ∧ Represents` | `relabelTree_wf`, `relabelTree_represents` | sorry |
| 5 relabel, per-node layout: `Layout.Shape`/`Layout.Local` for F, V, Q-loop/O, Q/I, P, S, R (`shape_*`, `local_*`), exact rows (`runF_row`, `runLoop_row`, `runQI_row`, `runP_row`, `runS_row`, `run_entries`) | `LayoutShape.lean` | proved (standard axioms); `r_skeleton_nodup` discharges the R `Nodup` hypothesis from `r_shape` |
| 5 relabel, structural part: `relabelTree_own : Items.WF → Items.OwnExtra → Items.ROriented → Bijections ∧ Ownership ∧ Twins` (also `relabelTree_bijections` from `WF` alone, `relabelTree_twins` from `WF` + `OwnExtra`); `layoutNode_edges` | `RelabelOwn.lean` | proved modulo `relabel_node_spec` (`RelabelSpec.lean`, sorry); `OwnExtra` is a stated hypothesis to be discharged by the walk |
| 7 st-order spec `StOrder`, `Items.StNumbered`, split `spqrTree_st = relabel_st ∘ walk_st` | `StSpec.lean`, `StWalk.lean` | def / proved split |
| 7 relabel-side: `vchildren_nv_increasing`, `orderedChildren_sorted`, `edgeChildren_dominance`, `layoutNode_r_bracket` | `StSpec.lean`, `StLayout.lean` | proved |
| 7 relabel-side: `relabel_st` | `StSpec.lean` | sorry |
| 7 walk-side: `WalkState.StInv`, data lemmas `pushTstack_onSide`, `merge_onSide`, `fold_onSide` | `StWalk.lean` | def / proved |
| 7 walk-side: `finishTstackTop_stItem`; ear lowvals `first_ret_lowval`, `chain_stackDir_step` | `StWalk.lean`, `StEar.lean` | proved |
| 7 walk-side: `StInv.onSide` field, `chain_stackDir_const` (corrected statement, see 7.4) | `StWalk.lean` | def / proved |
| 7 walk-side: `StInv.hole` (`StHole`/`HoleClosed`), `stInv_topClosable`, `stInv_finishTstackTop_stItem` (close site, see 7.4) | `StWalk.lean` | def / proved |
| 7 walk-side: `finishEdge_stInv` (under `FinishGuards`/`EarsOnSide`), `walk_stInv`; `walk_st`, `spqrTree_st` from `walk_stInv`; route changed to the `StRef.lean` reference order (§7.6) | `StWalk.lean`, `StRef.lean` | sorry / proved (`walkTree_stackDir_below` in `StFrame.lean` proved); reference tested, equality proof not started |

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
* Adjacency rows (`layoutNode_r_bracket` **[proved]**, `StLayout.lean`: `layoutNode .R` is
  rewritten as a fold `LayoutR.run` — count, prefix sum, reverse fill, cap — whose row contents are
  computed exactly by `LayoutR.run_spec`): R nodes count each edge `(a, c)`, `a < c`,
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
* What blocks `relabel_st` **[sorry]** itself is not a layout fact but the glue: a per-node
  characterization of the recursive `relabel` fold — for every node `i` of `relabelTree g items`,
  `nvRange i`, `neRange i`, `nodeEdgesOf i` and the rows `adjRow (2 nv)`, `adjRow (2 nv + 1)` are
  those of `layoutNode it.type … edgeChildren` for the item `it` numbered `i`, with `edgeChildren`
  the `vertPos` images of the non-V children of `orderedChildren it nvSt`.
  No such lemma exists yet (`relabelTree_wf` is admitted for the same reason), so it is the exact
  statement to add, as `relabel_node_layout`; with it `relabel_st` is `layoutNode_r_bracket` plus
  the S/P row shapes read off `layoutNode` and `Items.StNumbered` transported along `vertPos`.

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
The semantic half, `chain_stackDir_const` **[proved]** from the `StInv.onSide` field: an entry
finished by an open-path vertex on a chain of constant `stackDir` is `OnSide` that direction.
Its hypothesis `hchain` (constant `stackDir` along the chain) is the lowval fact, proved in
`StEar.lean`: `DfsOut.lowval_eq_lmin` (the lowval of a well-formed out-edge is the minimum of its
return depths, from `DfsOut.WF` = `dfsVisit_spec`'s per-edge classification), `first_ret_lowval`
(the first returning out-edge of a child classified `ret l _` — the first child continuing the
ear — returns to exactly `l`, using `classify_eq_ret_iff_tree` for `l = lmin` and the
`OutClass.rank` sortedness of `DfsTree.WF` to put the minimal edge first), and `chain_stackDir_step`
(`walkOut` at depth `d + 1` sets `stackDir[d+1] = !stackDir[l] = stackDir[d]`).
Threading these along the chain needs the frame fact `walkTree_stackDir_below` **[proved]**
(`StFrame.lean`: walking a subtree at depth `d` leaves `stackDir` below `d` unchanged), by the
mutual induction over `walkTree`/`walkOuts`/`walkOut` with the block decomposition `finishEdge_eq`:
the only `stackDir` write inside `finishEdge` is `loop1Type`'s, at the depth of a stack entry
strictly above `d`.
The conclusion `t.OnSide dir` is a history property of how the entry `t` was built (every
`pushTstack`/`mergeTstackTops` along the chain used the same `edgeDir`), not determined by the
current state plus `hchain`; it is the field `StInv.onSide`: every entry with `topDepth < d`
whose `vStart` is an open-path vertex strictly below its top (`stackVerts[j]` for some
`topDepth < j ≤ d`) is `OnSide stackDir[topDepth]`.  `chain_stackDir_const` is that field
rewritten along `hchain`.
Two earlier formulations are **false** (traces of the C++ walk, `gen.py` seeds 0..499):

* `chain_stackDir_const` with the hypothesis "`t.vStart ∉ stackVerts(l, d]`" and
  `t.topDepth ≤ d` (the entry was *not* finished on the chain) — seed 6, current depth
  `d = 14`, open path `0, 1, 18, 4, 19, 8, 15, 11, 13, 6, 3, 9, 12, 16, 2`,
  `stackDir = 0, 1, …, 1, 0`, `l = 13`: the entry `vStart = 10, topDepth = 14` is the back edge
  `10 → 2` pushed while `stackDir[14] = 1` (first out-edge of `2`, lowval 0), so it sits on
  `spans.2`; the next out-edge of `2` (lowval 4) reset `stackDir[14] = 0`.  Entries at
  `topDepth = d` from earlier out-edges of the current vertex are never closed individually
  (a later entry with `topDepth < d` always sits above them, or they are flattened by the
  parent's fold), so no clause of `StInv` speaks about them.
* "every live entry with `topDepth < d` not containing the V item of an open-path vertex is
  `OnSide stackDir[topDepth]`" — seed 6, `d = 11`, open path
  `0, 1, 18, 4, 19, 8, 15, 11, 13, 6, 3, 9`: the entry `vStart = 16, topDepth = 0` is two-sided
  (`spans.1 = Q, V10, S, Q, V2, Q, P, Q, Q`, `spans.2 = Q, Q, P, Q, V12`).  It is the type-2 fold
  at vertex `16` (depth 13) that then absorbed the eagerly merged vertex entries of `12`
  (depth 12, first out-edge type 2, `stackDir[12] = 1` vs. the fold side `!edgeDir = 0`); the
  vertices `16, 2, 10, 12` have since been popped.  Such entries stay two-sided until the first
  ancestor with `hasVert` folds them (type-2 branch) and are never closed before that.

Entries from the current subtree are one-sided in the way the closes need: at the start of
`finishEdge` for a returning tree edge (`lowval < d`), every entry above `origTstack` with
`topDepth ≥ d` is `OnSide stackDir[d]` (3686 checks, 0 violations; the chain fact for the
`topDepth > d` entries is `chain_stackDir_step`), which is the `CloseOK` of the `loop1`
closes; the only exception is the component-boundary case `lowval ≥ d`, where the child's
vertex entry (`topDepth = d + 1`, side `!stackDir[lowval]`) may be on the other side.  This
segment property is indexed by `origTstack` and belongs to the `FinishGuards`/`EarShape`
hypotheses of `finishEdge_stInv`, not to the depth-indexed `StInv`.
Genuinely two-sided entries thus arise only from the eager vertex merge, when a finished inner
ear is enclosed by pieces from both sides.

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
`finishTstackTop_stItem` **[proved]** turns `TopClosable` into `Items.StItem` for the closed item
(`finishTstackTop` writes exactly `entryVertList`/`entryEdges` into it; the item must not be one of
the entry's own items, which `allocItem`/`maybeUnwrapNxt` guarantee); `finishEdge_stInv`
**[sorry]** is the preservation of `StInv` through `finishEdge` (stated under the structural
`FinishGuards` and the named ear fact `EarsOnSide`, below); `walk_st` **[proved]** follows from
`StInv.items` at the end of the walk via the admitted walk-level induction `walk_stInv`
**[sorry]** (`∃ ord, (g.walk …).StInv 0 ord`, the `walkTree_inv'` induction).

**Closability.** `finishEdge_topClosable` as originally stated (`StInv d ord → TopClosable`) is
not provable: `StInv` recorded only the order of the vertices on the two sides, and the top
entry's skeleton edges may reach into its hole before the hole is closed.  The hole-closure
content is now the field `StInv.hole` (`StHole`, for every entry with `topDepth ≤ d`): every
skeleton edge `(p₁, p₂)` of the entry (`itemEdges`, the non-V items' `vs`) is oriented
`ord p₁ < ord p₂` with both ends *in reach* (`EntryReach`: a V vertex of the entry, a terminal, or
an open path vertex `stackVerts[k]`, `topDepth < k ≤ d`), and every V vertex other than the lower
terminal `vStart` has a lower and a higher neighbour among these edges.  The close-site condition
`HoleClosed d t` says every open path vertex of `(topDepth, d]` is `vStart`, a V vertex of the
entry, or touched by none of its edges.  `stInv_topClosable` **[proved]**: `StInv d ord`, top entry
`t` with `topDepth ≤ d`, `OnSide stackDir[topDepth]`, `HoleClosed`, and terminals distinct and not
among the entry's vertices give `TopClosable` for `t` (via the generic `stList_of_sorted`: a
strictly `ord`-sorted vertex list between two terminals with oriented, in-list edges and a
lower/upper neighbour per interior vertex is an `StList`); `stInv_finishTstackTop_stItem`
**[proved]** chains it into `finishTstackTop_stItem`.  The side/hole hypotheses are exactly what
the close sites of `finishEdge` provide: `OnSide` from `EarsOnSide` (loop1) / `fold_onSide`
(type-1 fold, `stackDir[lowval] = !edgeDir`), `HoleClosed` because the closed entry's hole is the
finished child path (type-1: `(lowval, d]`, all its vertices are V items of the merged entries or
`curV = vStart`), and the terminal facts from `makeVs` being a separation pair.

**Vertex entries.** The first versions of `StEntry.oriented`/`bottom`/`ends` and `StHole.lower`/
`upper` were stated for every vertex; they are **false** for a vertex entry: `gen.py` seed 5, the
snapshot at the start of out-edge 72 (depth 3, path `4, 0, 7, 14`) has the live entry
`vStart = 7, topDepth = 2, spans = ([], [V7])` — the vertex entry of `7`, pushed by `walkOut` and
sitting *below* `origTstack`, so that it is not merged by `7`'s own type-1 closes (those merge the
child's vertex entry: the three entries above `origTstack` at the merge are the child's vertex
entry, the lowval back edge and the tree edge) and carries no edge until `7` itself is finished;
its `V7 = vStart = top`.  Hence the clauses exempt `x = vStart` (and `x = top`), `ends` is
conditional on `vStart ≠ top`, and the close sites supply `vStart ∉ entryVerts`, `top ∉ entryVerts`,
`vStart ≠ top`.

**What was validated empirically** (trace of the C++ walk on `gen.py` seeds 0..149, every
tstack snapshot at the start of an out-edge; the Lean walk is byte-identical on these inputs):
the V items of every entry are distinct and increasing in the final order of their parent item
(`StSides`, 0 violations); the top entry is oriented towards `stackVerts[topDepth]` as in
`StEntry` (215 live cases, 0 violations; non-top entries whose `topDepth` slot has since been
reused by a sibling can violate it, which is why `StEntry` is stated for the top entry only);
every closed S/P/R item satisfies `Items.StItem` (986 items, 0 violations).
The naive strengthening "every entry's `nest` reading is an st-list with the hole at the far end
from the terminal" is **false** when the open path changes direction.
Concrete failing configuration (`gen.py` seed 20, the snapshot at the start of out-edge 72 of
the walk, current depth 11): the open DFS path is `39, 7, 8, 20, 36, 6, 4, 0, 3, 10, 27, 21`
(depths 0..11) with `stackDir = 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0` — depths 1..6 form the ear
with lowval 1 (`stackDir = 1`), the ear with lowval 3 starts at depth 7 (vertex 0) where
`stackDir` flips to 0, and depth 11 (vertex 21) has lowval 5.
The tstack entry in question is `vStart = 21, topDepth = 1` (upper terminal `stackVerts[1] = 7`),
`spans.1 = []`,
`spans.2 = [Q(21,5), V5, Q(5,7), Q(5,33), Q(3,33), V33, S(33,7), Q(5,6), Q(5,16), S(5,12), Q(27,12),
V12, Q(12,16), Q(10,16), V16, Q(16,6)]`.
Read as a closed item with the hole at the far end (`vs = (21, 7)`, children `spans.2`) its vertex
list is `21, 5, 33, 12, 16, 7`, and vertex `16` has neighbours `5` and `12` before it but none
after it: its other neighbours `6` (depth 5, on the `stackDir = 1` segment) and `10` (depth 9,
on the `stackDir = 0` segment) are open path vertices, i.e. they lie *in the hole*, on different
sides of the entry's items in the final order.
So the entry's items are not an st-list on their own; they only become one once the path
segments are inserted according to their per-depth `stackDir` — which is exactly what the
corrected invariant (`StSides` for every entry, orientation only for the top entry) does not
claim, and what the `StItem` check at every close (986 items, 0 violations) confirms.
The hole must be modelled as the path itself, ordered by `stackDir` per depth; `StHole`/
`EntryReach`/`HoleClosed` above do this by allowing edge ends on the open path and closing only
when the path segment of the hole has been absorbed.

### 7.6 The st-order reference (`StRef.lean`)

Instead of carrying `StInv` through `finishEdge`, the walk-side proof goes through a *reference
order*: an Even–Tarjan style description of the order in which the walk lists the children of its
S / P / R items, stated directly on the lowval-sorted DFS tree, without the tstack. Item
boundaries (which edges form which S / P / R item, `maybeUnwrapNxt` reuse) are taken from the
walk's `Items`; the reference recomputes only the order.

Every edge is handled at its deeper endpoint, in that vertex's sorted out-edge order (exactly where
`walkOut` / `finishEdge` see it). `refTree g t d dirs` walks the subtree `t` at depth `d` with
`dirs : List Bool` the directions chosen along the path above (`walkOut`'s `setStackDir`: `false`
for a block-boundary edge, otherwise `!dirs[lowval]`) and a `hasVert` flag, and produces the list
of *pieces* `⟨side, items⟩` in push order:

* the vertex item `V v` is a piece on side `!dirs[l]` at `v`'s first returning out-edge of lowval
  `l` (before the sub-walk if it is type 1, after the child and the tree edge if it is a type-2
  child — `finishTail`), on side `true` at the end of `v` if `v` has no returning edge;
* a back edge to depth `l` is a piece on side `dirs[l]`;
* a returning tree edge `(v, c)` of lowval `l` contributes the child's pieces followed by the tree
  edge `Q e` on side `!dirs[l]`; if `v` already has its vertex item (the `closeVert` case of
  `finishEdge`) the whole sub-ear is *folded* into the single piece `⟨dirs[l], stNest sub⟩`;
* an edge with `lowval ≥ d` starts a new block (its pieces are read separately).

`stNest` reads a block's pieces bottom-up as `L ++ R`: a piece on side `false` is prepended to `L`,
one on side `true` appended to `R`, so later pieces lie outside earlier ones (the walk's
`mergeTstackTops` puts the top entry's spans outside the next entry's: `(b.1 ++ a.1, a.2 ++ b.2)`).
The folds are the only non-flat ingredient: without them (pure nesting) the order is wrong on
143 of the seeds 0..300 (the sub-ear's right pieces would end up outside the parent's items);
with them `refOrder` restricted to an item's subtree equals the walk's `ch` on every S / P / R
item for seeds 0..300 (`check_stref`, `compare_stref.sh`) and, for the same rule prototyped in
Python, on seeds 301..1000. (`finishTstackTop`'s re-siding of a closed entry, `maybeUnwrapNxt`'s
unwrapping and the P-closes do not move items relative to each other, which is why they do not
appear in the reference.)

`walk_st` is now derived from two admissions: `walk_st'` — for every S / P / R item, `ch i` is
the restriction of `refOrder` — to be proved by simulation in the style of `Sim.lean`, relating
the reading `readStack` of the live tstack above the current block's base to `stNest` of the
pieces pushed since (`readStack_pushTstack`, `readStack_mergeTstackTops`, `readStack_fold`,
`readStack_finishTstackTop` are the per-primitive steps; the closes are up to `expandItem`); and
`stItem_of_refOrder` — an item listed in the reference order is in s-t order — by induction over
`refTree` without the tstack (every piece spliced at depth `l` joins the open path on side
`dirs[l]`, so every vertex other than a block's endpoints has a neighbour on each side; the
restriction to an item keeps this since the item's virtual edges are the ends of the sub-ears).
The one walk-side ingredient of the latter is the orientation of the children's `vs`: `makeVs`
uses the same `stackDir[d]` as the splice side, which is why `stItem_of_refOrder` is stated for
the walk's items rather than for arbitrary `WF` items with `ch` in reference order (for those the
`vs` could be flipped). The §7.4 `StInv` route (`finishEdge_stInv`, `walk_stInv`,
`walk_st_of_stInv`) is kept as a documented alternative and no longer feeds `walk_st`.

### 7.5 Work packages

| lemma | file | status |
|---|---|---|
| `SpqrTree.StNumbered`, `EdgeDominance`, `AdjBracket`, `StOrder` | `StSpec.lean` | def |
| `vchildren_nv_increasing` | `StSpec.lean` | proved |
| `Items.StList`, `Items.StItem`, `Items.StNumbered` | `StSpec.lean` | def |
| `orderedChildren_eq_of_ne_R`, `orderedChildren_sorted` | `StSpec.lean` | proved |
| `pairwise_dominance_of_sorted_sum`, `edgeChildren_dominance` | `StSpec.lean` | proved |
| `LayoutR.run`, `layoutNode_R_eq`, `run_spec`, `layoutNode_R_bracket` | `StLayout.lean` | proved |
| `layoutNode_r_bracket` | `StSpec.lean` | proved |
| `relabel_st` | `StSpec.lean` | sorry (7.3: needs the relabel-fold glue `relabel_node_layout`) |
| `TEntry.wrap`, `OneSided`, `OnSide`, `nest` | `StWalk.lean` | def |
| `pushTstack_onSide`, `merge_onSide`, `fold_onSide`, `getSide_setSides` | `StWalk.lean` | proved |
| `WalkState.StSides`, `StEntry`, `StInv`, `TopClosable`, `entryVertList`, `entryEdges` | `StWalk.lean` | def |
| `DfsOut.lowval_eq_lmin`, `first_ret_lowval`, `chain_stackDir_step` | `StEar.lean` | proved |
| `walkTree_stackDir_below` (frame: `stackDir` below `d` unchanged by `walkTree _ d`) | `StFrame.lean` | proved |
| `StInv.onSide`, `chain_stackDir_const` (`ear_uniform_side`, semantic half; corrected statement) | `StWalk.lean` | def / proved |
| `finishTstackTop_items`, `finishTstackTop_stItem` | `StWalk.lean` | proved |
| `StInv.hole` (`StHole`, `EntryReach`, `HoleClosed`), `idxOf_lt_idxOf_iff`, `stList_of_sorted`, `stInv_topClosable`, `stInv_finishTstackTop_stItem` | `StWalk.lean` | def / proved (replaces `finishEdge_topClosable`, see 7.4) |
| `EarsOnSide` (named ear-shape hypothesis), `finishEdge_stInv` | `StWalk.lean` | def / sorry (hard; push/merge/fold/close blocks via `Step`/`Sim`) |
| `walk_stInv` | `StWalk.lean` | sorry (the `walkTree_inv'`-shaped induction; needs `StInv (d+1) → StInv d` at returns) |
| `walk_st`, `spqrTree_st` | `StWalk.lean` | proved (from `walk_st'`, `stItem_of_refOrder`, `walk_items_wf`, `relabel_st`); `walk_st_of_stInv` is the same from the alternative `walk_stInv` route |
| `refTree`/`refOrder`, `restrictCh`, `check_stref` differential test (§7.6) | `StRef.lean`, `CheckStRef.lean` | def / tested seeds 0..300 (0 mismatches) |
| `walk_st'` (`ch i = restrictCh … (refOrder …) i` for S/P/R items) | `StRef.lean` | sorry (the simulation) |
| `stItem_of_refOrder` (an item listed in the reference order is in s-t order) | `StRef.lean` | sorry (Even–Tarjan on the reference, induction over `refTree`; needs the `makeVs` orientation of the children's `vs`) |
| reading a tstack as pieces: `readStack`, `stNest_append`, `readStack_pushTstack`, `readStack_mergeTstackTops`, `readStack_fold`, `readStack_finishTstackTop` (per-primitive steps of the simulation relation `readStack stack = stNest pieces`, up to `expandItem` at closes) | `StRef.lean` | proved |

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
| orbit counting: `numOrbits` = orbit minima, `numOrbits_involution`; `numComponents` by label relaxation | `Planar.lean`, `Proofs/Planar.lean` | def / **proved** |
| closed forms of `layoutRot .S` / `.P`; `cycleRot_isPlanarEmbedding`, `bondRot_isPlanarEmbedding` | `PlanarLayout.lean` | **proved** (`propext`, `Classical.choice`, `Quot.sound`) |
| `layoutRot_S_shift`, `layoutRot_P_shift` (renumbering a node's edges from 0) | `PlanarInv.lean` | **proved** |
| Invariant P (§8.2) as a Lean structure: `Piece`, `SideWalk`, `InvariantP`, `StackInv` | `PlanarInv.lean` | def (the deliverable is the statement) |
| `planarWalkOut_stackInv` (the walk preserves Invariant P) | `PlanarInv.lean` | sorry |
| `neRotAdj_segment` (relabel bookkeeping: node `i`'s `neRotAdj` segment is its `layoutRot`) | `PlanarSpec.lean` | sorry |
| `nodePlanar_sound_S`, `nodePlanar_sound_P` | `PlanarSpec.lean` | proved modulo `neRotAdj_segment` and `Shape` (via `spqrTree_wf`, itself admitted in `relabelTree_wf`) |
| `nodePlanar_sound_R` (Invariant P at finish, mapped by `mapRot`) | `PlanarSpec.lean` | sorry |
| `nodePlanar_sound` = S ∨ P ∨ R cases | `PlanarSpec.lean` | proved from the three |
| `nodePlanar_complete` (Kuratowski-style certificate from the §8.3 crossing) | `PlanarSpec.lean` | sorry, hard |
| `twoSum_planar`, `oneSum_planar`, `disjointUnion_planar` (gluing on `Planar`) | `PlanarInv.lean` | sorry (explicit splice with `f₁ + f₂ − 2` faces: `PlanarGlue`) |
| `WalkInv` / `Preserves` (Hoare triple on `PlanarWalkM`), `Preserves.frame`, `Preserves.popPair`, `planarFinishEdge_inv`, `planarWalkOut_inv` (+ `Tree`/`Outs`), `planarWalkOut_stackInv` | `PlanarInvSteps.lean` | def / **proved** modulo the per-step lemmas |
| `pushVertTstack_inv` (a fresh vertex entry has no exposed ends, so `StackInv` exempts it) | `PlanarInvSteps.lean` | **proved** |
| per-step lemmas `pushEdgeTstack_inv`, `mergeTstackTops_inv`, `maybeUnwrapNxt_inv`, `finishTstackTop_inv`, `closeBackedges_inv`, `flipBeforeMerge_inv`, `pruneBackedges_inv`, `flipForLowval_inv`, `modifyCur_foldSides_inv` (one `planarFinishEdge` step each preserves Invariant P) | `PlanarInvSteps.lean` | sorry |
| gluing invariant `GluedUpTo` (per maximal processed item: restricted rotation agrees with a planar embedding of the piece below it, exposed ends = unset ends, on one face), `gluedUpTo_init`, fold `forM_reverse_range_inv`, `gluedUpTo_planarEmbed` | `PlanarEmbedSteps.lean` | def / **proved** modulo the steps |
| `embedItem_step_F` / `_V` / `_Q` / `_leaf` / `_node` (one `embedItem` preserves `GluedUpTo`) | `PlanarEmbedSteps.lean` | sorry (F: `disjointUnion_planar`; V, Q: `oneSum_planar`; node: `twoSum_planar` + `nodePlanar_sound`) |
| `glued_root` (`GluedUpTo 0` at the root `F` item is `IsPlanarEmbedding` of `g`) | `PlanarSpec.lean` | sorry |
| `planarEmbed_sound` | `PlanarSpec.lean` | proved from `gluedUpTo_planarEmbed` + `glued_root` |
| `spqrTree_planar` (`→` from `planarEmbed_sound`; `←` needs completeness + skeletons are minors of `g`) | `PlanarSpec.lean` | sorry |

Admitted, precisely (`#print axioms` reports `sorryAx` for each): `neRotAdj_segment`;
`nodePlanar_sound_R`; the nine per-step lemmas `*_inv` of `PlanarInvSteps.lean` (hence
`planarWalkOut_stackInv`); `nodePlanar_complete`; `twoSum_planar`,
`oneSum_planar`, `disjointUnion_planar`; `embedItem_step_F`, `embedItem_step_V`,
`embedItem_step_Q`, `embedItem_step_leaf`, `embedItem_step_node`, `glued_root` (hence
`planarEmbed_sound`); `spqrTree_planar`. The S and P
cases of `nodePlanar_sound` are proved except for `neRotAdj_segment` and the `Shape` of the
skeleton (`spqrTree_wf`, which inherits `relabelTree_wf`'s `sorry`); the local counting itself
(`cycleRot_isPlanarEmbedding`, `bondRot_isPlanarEmbedding`) is fully proved.

Work packages:
* **Invariant P**: `InvariantP` / `StackInv` (`PlanarInv.lean`) are stated; prove
  `planarWalkOut_stackInv` step by step for `makeEdgePlanarity`, `mergeSide`,
  `closeSide`/`pruneSide`, `foldPlanarity`, `finishMatches`/`unwrapPlanarity`; the merge-side
  crossing argument gives the `none` case.
* **Relabel transport**: `neRotAdj_segment` (induction over `planarRelabel`: `neRotAdj` grows by
  exactly the node's `layoutRot` when `neBounds` is pushed).
* **Local embeddings**: S and P done; R (`nodePlanar_sound_R`) from Invariant P at the finish of
  the R item via `planarRelabel`'s `mapRot`.
* **Gluing**: `twoSum_planar` (explicit splice, `PlanarGlue`), `oneSum_planar`,
  `disjointUnion_planar`; the bottom-up induction over `embedItem` is done
  (`forM_reverse_range_inv`), what remains are the per-item steps `embedItem_step_*` (each one
  `link` = one splice of two pieces' exposed ends) and `glued_root`.
* **Completeness**: the crossing of §8.3 as a `K₅`/`K₃,₃` subdivision — the hard, optional one.
