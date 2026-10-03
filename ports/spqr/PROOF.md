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
terminals of the entries above: `TEntry.Term' D s above t v := Term D s t v ∨ ∃ t' ∈ above, t'.Term
D s v` (`WalkState.EntryInv'`/`Inv'` in `WalkSpec.lean`, implied by `Inv`). With that clause `D = d` (the current
depth, under tree edges too) is enough — 0 violations on 9000 random multigraphs (≤ 12 vertices,
≤ 24 edges) at every `walkTree` start/end and after every `walkOut`, while `Term d` fails 1827
times and `Term (d+1)` 728 times on the same snapshots; connectivity of every open entry also had
0 violations. `Step`, `BStep` and the per-block lemmas of `WalkSpec.lean` are now stated for
`Inv' D`; the extra clause costs three ear facts the depth-indexed version did not need: entries below
a merge/re-target are edge-disjoint from the touched entries (`MergeTopOk`, `RetargetOk.disj`; the
context of a popped entry is otherwise lost), a popped block's terminals are touched by no entry below
(`BoundaryOk.gone`/`gone₂`, the articulation vertex `curV` separates it), and lowering `Inv' (d+1)`
to `Inv' d` after a tree edge (`WalkInv.ear_lower'`, derived from `EarFinish.lower`: Loop 1 leaves every entry still attached at
`stackVerts[d+1]` as a terminal or under an entry with that `vStart`). Closed items stay strictly
2-attached: `finishTstackTop` acts on the top entry, whose context is empty (`Inv'.finishTop`), and at
Loop 1's R close the only context terminal of `nxt` is `cur`'s, interior to the union (`RStep.union_attached`).
The walk induction `WalkInv.walkTree_inv'` is restated for `Inv' d` and proved modulo the `ear_*` admissions.

**Guards from the ear content (ear session 5, `EarShape.lean`).** `EarShape.finishGuards`
(`FinishGuards` from the predecessor's depth/first-index shape `EarShape` + `E.bot = d+1`,
`E.top ≤ d`, `lowval = E.lam`, type-2 tree edge into the head's `vStart`) is false: the one-entry
stack `[(0, d+1, 0, _)]` at `d = 1` satisfies every field of `EarShape ⟨0,1,2,0,0,0⟩` but `Inv1 1`
needs an entry of `topDepth < 1` (`EarShape.finishGuards_false`, kernel-checked). The shape the
walk actually carries is `EarFinish` (via `BookTree`), so `EarShape.finishGuards`/`finishEdge`
are replaced by `finishGuards_of_ear : EarAt → FinishGuards` (boundary cases from
`bd_bridge`/`bd_comp`, `Inv1` from `loop1`'s `Units` range + `bottom`, the vertex-ear lengths from
`loops`) and the new field `late_fo` (after loop 1 the bottom entry predates `firstOccurrence[d]`
— `Inv2` of `FinishGuards`; 0 violations on 3000 seeds). `GuardsTree` follows from `BookTree` by
the walk induction (`guards_of_book`/`walkTree_guards'`), so `walkTree_guards` is exactly the
ear contract `walkTree_book`.

**Ear content, restated (ear session 4, `EarInv.lean`, `checks/EarCheck.lean`).** The first
version of `WalkState.EarFinish` had fields that are false for the real walk (3000 random
multigraphs, every `finishEdge`): `base_top` (`∀ t ∈ base, topDepth ≤ d` — a vertex entry `V y` of a
finished deeper vertex stays buried under a lower piece: seed 1, `finishEdge` of `3→2` at `d = 1`
with `(1,4)=[V 1]` on the stack), `loop1_bot` (`topDepth > d ⇒ vStart = child` in loop 1's range —
an S-merged entry keeps a deeper chain vertex as `vStart`), `touch_top` for all entries of
`topDepth ≤ d+1` (buried `V` entries hold only block `Q` items), `vert`/`vert_free` with
`vStart = curV` (after a late merge the entry holding `V curV` has a chain vertex as `vStart`),
`loop1_side` without `lowval < d`. They were dropped or weakened (`loop1_touch` only in loop 1's
range, `vert` only `topDepth ≤ d`) and the following fields added, all with 0 violations: `sub_bot`,
`path` (distinct path vertices, for `FinishTopOk.mid`), `span_disj`, `p_entry` (the P-merge target
`(curV, lowval)` in `base` is a single root item on side `stackDir[lowval]` attached only at `curV`,
`stackVerts[lowval]`), and the author's two structural facts: `bottom` — for a returning tree edge
the child's entries end, bottom-up, with the vertex entry of the chain bottom `y` and the `(y,
lowval)` piece (a single root item on side `stackDir[lowval]`, touching `y` and `stackVerts[lowval]`;
`EarBottom`), and both are byte-identical to when `y`'s walk ended at every inner `finishEdge` and
after loops 1–2 (checked on 6945 returning tree edges); `loops` — after loops 1–2 the child's entries
are `[c, py, vy]` for a type-1 edge (`c` reaches `curV` and `o.dest`; the three cover the sub-ear
edges) and `c :: mid ++ [py, vy]` for type 2, with `lowval ≤ c.topDepth
≤ d`. Not true (so not stated): `c.topDepth = d`, `c.vStart = o.dest`, `c` one-sided on
`stackDir[d]`, `c` not touching `stackVerts[lowval]` (loop 2 can merge a lowval piece into `c`, making it two-sided with `topDepth =
lowval`: seed 2, `finishEdge` of `4→1` at `d = 1`). `FinishBook` now carries `EarAt` (`origTstack`
indexed), `finishOk_of_guards` and the `ear_*` admissions take `EarFinish` as hypothesis.

**Loop 2 and the vertex close (ear session 5, `EarInv.lean`: `EarLate`, `EarClose`, `FoldSpec`).**
`ear_mergeLate` needs `MergeOk` at every loop-2 merge; `ear_closeVert`/`ear_finishP_vert`/
`ear_tail_tree` need the state after loops 1–2 beyond the shape `loops` gives. Both are now fields of
`EarFinish`, read on the intermediate states the block lemmas are stated at (loops 1–2 only merge or
close entries, so the facts are stated where they are consumed rather than transported): `late` —
at `feS₁` the stack is `c₀ :: R`, pairwise edge-disjoint, and `l2Cur c₀ done = foldl mergeInto`
merged with the next entry is `MergeOk (d+1)` at every split reached while `firstIdx >
firstOccurrence[d]` held (`loop2Cond`; `mergeInto` keeps `nxt.firstIdx`, so the condition at
iterate `k` is about the `k`-th merged entry); `close` — at `feS₂` the child's entries are
`c :: mid ++ [py, vy]` over `base` (whose edge sets are unchanged), pairwise edge/span-disjoint, hold
exactly the sub-ear edges (`sub_edges`/`sub_cover`), every edge at the chain bottom `y =
vy.vStart` is a sub-ear edge (`y_edges`: `y` is in the child's subtree, so `y` is *interior* to the
fold — this is how `RetargetOk.old` and `MergeOk.bottom` hold although `vy.vStart ≠ o.dest` for a
type-2 chain), `py` is still the single root item of `EarBottom` on side `stackDir[lowval]`, folding
`mid ++ [py, vy]` into `c` top-down is `MergeOk (d+1)` at every step (`FoldSpec`: loop 3, then the two
merges of the vertex close; the `!hasVert` merge of `V curV` into `c` uses `vert_disj`/`vert_touch`:
the vertex item's edges are on no entry and touch `curV`), and for type 1 `mid = []` and the three
entries touch only `curV`, `stackVerts[lowval]` and interior vertices (`FinishTopOk.mid` of the
close and of the P merge). Path bookkeeping (`sv_d`, `sv_child`, `path_child`, `dir_d`) is
recorded as fields too; `vy_top` (`d < vy.topDepth`, the `EarBottom` fact read at `feS₂`) gives
`topDepth` of the type-1 close result `= lowval`; `base_root` (a root span item of a `base` entry at
`s` is still a root at `feS₂`: loops 1–2 never make a `base` item a child) keeps `p_entry`'s root
fact alive for the P-check after the vertex close; `c_edge` (the top piece `c` owns the tree edge
`o.e` at `feS₂`) gives the shared vertex of the `!hasVert` merges of `V curV` into `c`. All 0 violations on 3000 random multigraphs
(`closeCheck`: 1864 fold merges and 301 late merges on the first 300 seeds alone); nothing had to be
restated. Frame facts (`g`/`stackVerts`/`stackDir[lowval]` at `feS₂`) are part of `EarClose` for
now; they are derivable from the block `Step`s and move out of the contract once the preservation
induction carries it.

**Boundary and lowering fields (ear session 5).** `ear_boundary` and `ear_lower'` are derived from
five more `EarFinish` fields. `dest_edges` (every edge at `o.dest` is a `subEdges` edge after a tree
edge, so by `base_disj` no `base` entry touches `o.dest`), `bd_noVert` (`d ≤ lowval → hasVert = false`:
the vertex entry of `curV` is pushed only after its boundary edges; gives `BoundaryOk.v_free`),
`bd_bridge` (`lowval = d + 1`: `sub = [(o.dest, d + 1)]`), `bd_comp` (`d ≤ lowval ≠ d + 1`:
`sub = [(o.dest, lowval), (o.dest, d + 1)]` — the block piece sits on the child's vertex entry),
`bd_term` (an entry of `tstack.tail` touching `curV` has `vStart = curV` or `topDepth ≤ d`, the
terminal `BoundaryOk.gone`/`gone₂` need), and `lower` (the `ear_lower'` fact at the final state, with
the `g`/`stackVerts` frame until the preservation induction supplies it). `bd_side` (ear session 6)
is the orientation half of `BoundaryOK` (§4.4) at the boundary: the bridge entry has `spans.1 = []`;
at a component edge the `(o.dest, lowval)` entry has `spans.2 = []` and the `(o.dest, d + 1)` entry
below it `spans.1 = []` (0 violations on 3000 seeds). The first boundary version
was **false**: `bd_shape` (`sub = [t]` with `t.topDepth = lowval`) and `bd_vert` (the entry under it
is `(curV, d)`) fail on almost every component edge (2942 of 3000 seeds; seed 0, `finishEdge` of
`2→1` at `d = 2`, `lowval = 2`: `tstack = [(1, 2), (1, 3)]`, the second popped entry is the child's
vertex entry `V 1`, not `(curV, d)`). `BoundaryOk.gone`/`gone₂` are now stated under
`o.cls.isTree = true` (and `lowval ≠ d + 1` for `gone₂`), exactly where `finishBoundary_inv` uses
them: for a back edge at the boundary (a self-loop, `lowval = d`) the top entry is an arbitrary
`base` entry and nothing separates it from the ones below. All 0 violations on 3000 random multigraphs.

**Export to the range layer: `WalkState.Frontier` (`Proofs/RInvFrame.lean`, stated by the R
session; derived in `EarFrontier.lean`).** The one positional fact crossing into the range/saturation
layer, stated with no reference to `EarFinish`: for `finishEdge _ d o origTstack _` at `s`, the
entries strictly above the `origTstack` enclosing ones (`tstack.take (length - origTstack)`) plus
the pending `o.e` own exactly `subEdges o` (`owns`), the enclosing entries own none of them
(`base_disj`), and at every iterate of loops 1–3 reached by the loop condition `FrontierOwns` holds
(the same `base` is the bottom, the entries above it own exactly `subEdges o`) with enough entries
strictly above the boundary for the next body (`loop1`/`loop2`/`loop3`). `finishEdge_frontier`
derives it from `FinishBook.ear` (`EarAt`) + `Inv'`/`Shape`: `owns`/`base_disj` from
`sub_edges`/`sub_cover`/`base_disj`; `loop1` by following the loop-1 iterates through `L1Inv`
(`L1Ctx.ofEar`, `l1_init`, `l1_iter`: the kept/closed pieces own `subEdges o`, the `d`-or-deeper
part of `sub` bounds the iterate count); `loop2` by placing the loop-2 iterates between `EarLate`
(`feS₁`, via `iter_merge_eq`) and `EarClose` (`feS₂`), merges preserving ownership (`l2Cur_edges`);
`loop3` by folding the `feS₂` shape `c :: mid ++ [py, vy] ++ base`. `walkTree_frontiers` gives
`FrontiersTree` by the walk induction (same hypotheses as `walkTree_inv'`: `Inv'`, `Shape`,
`GuardsTree`, `BookTree`). Standard axioms.

*Derivations (`EarLoop2.lean`)*: `ear_mergeLate` from `late` (`mergeLateOk_of_late`: `iter
mergeTstackTops k` on `c₀ :: R` is `l2Cur c₀ (R.take k) :: R.drop k`, the loop-2 condition at every
`j ≤ k` gives the `firstIdx` facts, `MergeTopOk` from `merge` + `disj`). `ear_closeVert` from
`close` (`closeVertOk_of_close`): the bounded loop 3 runs `k` iterations (`loop_run_iter`) with each
`MergeTopOk` from `FoldSpec` + pairwise disjointness (`mergeTopOk_iter_fold`); type 2 stops at
`[c, py, vy] ++ base` by `loop3Cond` (`origTstack = base.length`, now a hypothesis `hlen` of
`finishOk_of_guards` from `EarAt`), the two merges are `fold`'s last two steps, the retarget's `old`
is `y`-interior (`y_edges`/`sub_cover`) and `topDepth ≤ d` (`c_top`); type 1 unwraps `py` (its single
root item, `py_item`; `maybeUnwrapNxt_edges` keeps edge sets, so `MergeOk` transports along
`MergeOk.congr_entry`), merges, retargets (`old` by `type1`'s touch set with `path`/`sv_child`/
`path_child` separating `curV`, `stackVerts[lowval]`, `o.dest`), and `vertFinish` closes a
`(curV, lowval)` entry on `!stackDir[lowval] = stackDir[d]` (`dir_d`, `py_top`, `vy_top`), with
`FinishTopOk.mid` again from `type1`. The call site hoists the `Step`s to `feS₁`/`feS₂` (shape at
`feS₂`), as `hD`/`hlen`/`hs₂` extra hypotheses. `ear_finishP_vert` from `close`
(`finishPOk_type1_of_close`): type 2 has `condP` false; for type 1 the close is replayed as in
`closeVertOk_type1` (`cvS₂`..`cvS₅`, `finishTstackTop` modifies only the fresh `item`), so `feS₃` is
`r' :: base` with `r'` a `(curV, lowval)` single-item entry whose edges are `c ∪ py ∪ vy`
(`maybeUnwrapNxt_edges`, now also returning that the unwrapped `nxt`'s spans are old spans or
children of the unwrapped item; `maybeUnwrapNxt_ch`: the unwrap keeps every `ch`). If `condP`
holds, `nxt = b ∈ base` is `p_entry`'s `(curV, lowval)` root piece (`base_root` for rootness at
`feS₂`, `IsParent_modify` for the close), and the rest is `ear_finishP_back`'s argument: `UnwrapOk .P`
from `span_disj`/`base_root`, `MergeOk` via `touch_bot`/`c_edge`+`sub_cover`/`hends`, `FinishTopOk`
with `mid` from `type1`'s touch set and `p_entry`'s `att`, separated from `stackVerts[lowval]` by
`path`/`path_child`. Extra call-site hypotheses: `hi₂`/`hs₂` (`Step.mergeLate`), `hs₃`
(`Step.closeVert'` from `ear_closeVert`), `he`, `hends` (tree form).
`ear_tail_tree` (tree edge, no vertex entry, not `feSingle`): `ear_condP_tree` factors
`ear_finishP_tree`'s argument, so `finishP` is the identity on `feS₂`; `pushVertTstack` pushes the
`V curV` entry `ve` (`run_pushTstack`) on top of `c :: mid ++ [py, vy] ++ base` (`close`), and
`MergeTopOk` of `ve`/`c` is `vert_touch`+`c_edge`+`hends` (share `curV`), `vert_disj` (disjointness
from every entry), `sv_d`/`c_top` (bottom). Extra call-site hypotheses: `hsingle`, `hD`, `he`, `hends`.

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

**The `walkTree_ear` induction (`WalkInv.lean`, `cTree`/`cOuts`/`cOut`).** Mutual structural
induction over `walkTree`/`walkOuts`/`walkOut`, mirroring the frame `kTree`/`kOuts`/`kOut`
(`EarFrame.lean`), under `Types g s`, the DFS side conditions (`WF`, `Ends`, nodup of
`anc ++ v :: vertsList`, bounds, rank-sortedness of the outs) and the array sizes.
`CTree t d s`: `EarTree t d s` and, after the walk, `TreeEnd` — the end-of-outs context
`EarCtx t.v d done' [] hv' base bE sv sd sE` with `done'.map (·.1) = t.outs`, the end push
`L'`/`push'` (`push' ↔ hv' = false`) and the final `stackDir` write `D₃` (`pushEnd sE D₃ L'`),
exactly the data `TreeSite` consumes at the parent's site. `COuts`: from `EarCtx v d done rest
hasVert base bE sv sd s` at the state before `walkOut v d o` (`o :: rest`), `EarOuts` and the
end-of-outs `EarCtx`. `COut`: `EarOut` and `EarCtx` with `(o, hasVert')` appended to `done` after
`walkOut`. `walkOut` is split by `walkOut_eq` into `walkOutPre` (`wp_walkOutPre`: the `L`/`push`
stack shape) and `walkOutRest`; `EarOut`'s per-out body is `OutMid`. At a back out the site is
`earAt_back_of_ctx` and the step `ctx_step_back`. At a tree out `(e, cls, node y outs')`:
`ctx_init_child` gives the child's entry context over `L ++ s.tstack`; the IH `cTree` gives
`EarTree` of the child, its `TreeEnd` and its item frame (`KeptRel` + freshness, below);
`kTree`/`walkTree_stackDir_below`/the item frame
(hence `walkTree_qroot_kept`/`walkTree_vroot_kept`)/`walkTree_below_kept` give the frame below
`d + 1`, the untouched items and the two foreign roots `edgeItem e`/`vertItem v`; `endsOut_wf` (every edge below a WF out ends at `v`, in its subtree,
or at an ancestor) separates the done outs' edges from the child's vertices (`done_nc`) and
`ends_of_wf_boundary` gives `TreeSite.ends`; this assembles `TreeSite`, hence `EarAt`
(`earAt_tree_of_ctx`) and the step (`ctx_step_tree`). `FinishEar.vert` at the site is `vert_book`
transported along `walkTree_below_kept`. The root (`walkTree_ear`) is `ctx_init_root` + `cTree`
at `d = 0` with `base = []`.

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

**Lean status (`EarSides.lean`).** `finishEdge_sides` proves `FinishSides curV d o origTstack hasVert s`
(every `MergeOK`/`CloseOK`/`UnwrapOK`/`BoundaryOK` site inside one `finishEdge`) from
`FinishGuards`, `FinishBook` and `Inv'`/`Shape`, with no new side invariant: `CloseOK` is
`FinishTopOk.side`, `UnwrapOK` is `UnwrapOk.unwrap` (`UnwrapAt.single`/`side`), `MergeOK` at the
loop-1 sites is `loop1Cond`, at the loop-2/3 sites the `origTstack + 2` bound of
`WalkState.Frontier.loop2`/`loop3` (`finishEdge_frontier`), after the vertex unwrap the `3 ≤`
clause of `FinishGuards`, and at the vertex push the stack is nonempty (`loops`) and untouched by
`finishP` (`ear_condP_tree`); `BoundaryOK` is `bd_side` plus `v_root`/`vert_free`/`bd_noVert`
(no span item is below the fresh `vertItem curV`). `walkTree_sides` threads it through the walk
(`SidesTree`, the `walkTree_frontiers` induction). `EarRoot.lean` does the forest level:
`walkTree_rootOK` (after `walkTree t 0` from an empty stack the tstack is the root's single vertex
entry on side 2 — every depth-0 out-edge is a boundary, `finishEdge_root` pops exactly the child's
entries by `bd_bridge`/`bd_comp`, with the frame `Fr` = `g`/`stackDir.size` fixed), `init_shape`/
`init_inv`, and `sidesForest_of_roots`/`walk_sides_of_roots`: `SidesForest forest (init g tern)`
from `RootsBook forest (init g tern)` = `GuardsTree`/`BookTree`/`Inv' 0` at the start of every root,
threaded through the root pop/append (`Shape`, empty stack and `stackDir.size` are carried; `Inv' 0`
at a later root is an input because the root append needs `rootItem` parentless, a placement fact).
`RootsBook forest (init)` is derived in `EarWalk.lean` (`rootsBook_of_state`, threading `RootState` =
`Place` of the roots walked so far, `g`/array sizes, `Shape`, empty stack, `Inv' 0`, and the `ch` of
every fixed item not yet walked empty): `BookTree` at a root is the admitted `walkTree_book` (restated:
`t.WF []`, the endpoint fact `t.Ends g`, bounds/nodup, arrays sized `g.nv`, empty stack, `Inv' 0`,
`Shape`, and freshness — childless and parentless — of the tree's *own* `V`/`Q` items only; the old
`hfresh`, every fixed item childless, holds only at the first root), `GuardsTree` is `walkTree_guards'`
(from `BookTree`), `Inv' 0` at a later root is `Inv' 0` after the previous root walk (`invTree`) through
the pop (the stack is `[tt]` by `RootOK`, so `Stack` is vacuous) and the append onto `rootItem`
(`Inv'.modifyCh`, `rootItem` parentless by `Place.root`), and the later root's freshness is `Place`
(parentless: `cnt = 0` for an unpushed fixed item) plus the frame `walkTree_frame` (`EarFrame.lean`,
proved: `g`, the array sizes and the `ch` of fixed items outside the subtree are untouched). The frame
needed a typing hypothesis `Types g s` (`g`, `1 + nv + ne ≤ items.size`, `rootItem : F`, `vertItem : V`,
`edgeItem : Q`; `Place.types` supplies it): without it the statement does not follow, since
`maybeUnwrapNxt ty` reopens an existing item of type `ty` and writes its `ch`, so a fixed item typed
`S`/`P`/`R` could be rewritten. The proof is a `Keep j s₀ s` relation (fixed base state `s₀`: `g`, the
array sizes, `items.size ≤`, the types of the original items and `ch j` are kept) threaded through
`keep_maybeUnwrapNxt` … `keep_finishEdge` and the mutual `kTree`/`kOuts`/`kOut`, for every `j` that is
not a `vertItem`/`edgeItem` of the subtree; `Types.of_keep` transports the typing to the current state
and `Types.ne_of_type` rules out `j` as the unwrapped item. `walk_sides`/`walk_full` and their consumers (`walk_tstack_nil` … `walk_tree`,
`WalkTree`) moved from `WalkCover.lean` to `EarWalk.lean` (below `EarRoot`/`EarShape`), statements
unchanged except for the added hypothesis `hends : ∀ t ∈ forest, t.Ends g` (`DfsTree.Ends`,
`WalkInv.lean`: a tree edge joins the child's vertex and `v`, a back edge joins `v` and its destination
— `FinishBook.ends` needs it and `t.WF []` does not imply it); `walk_sides` is now
`walk_sides_of_roots` + `rootsBook_of_state`, so its only admission is `walkTree_book`. What is left
in `walkTree_book` is the ear preservation induction itself: `FinishBook.ear` is the whole `EarFinish`
contract, and the public `walkTree_guards` (guards only come from the ear, `finishGuards_of_ear`),
the tree-edge case of `earOut_one_entry` and the restated `finishGuards`/`finishEdge` all reduce to
that one induction. The split (statement unchanged): `EarTree` (= `BookTree` with the `ear` and
`vert` fields only, `FinishEar`) is the single admission `walkTree_ear`, and `walkTree_book` is proved
from it by the mutual `bTree`/`bOuts`/`bOut` induction (`WalkInv.lean`, standard axioms) over the
ancestor path `anc` (`d = anc.length`, `anc[k] = stackVerts[k]` for `k < d`, `stackVerts[d] = v`):
`ends` is `DfsTree.Ends` + `WF` (`ret lv` with `anc[lv]` the back edge's destination; the tree edge's
other end is `stackVerts[d]` after the child walk) through the frame `Keep D j` of `EarFrame.lean`,
now also keeping `stackVerts[k]` for `k < D` (`D ≤ d` for `walkTree t d`, `D ≤ d + 1` for
`walkOut`/`walkOuts`), and `q` is the freshness of the tree's `Q` items transported by `Keep` with
the exclusion set the edges of the *remaining* outs (a `Q` item is written only at its own
`finishEdge`; `walkOutPre` touches no items, `keep_walkOutPre`). The public `walkTree_guards`
(`EarShape.lean`) is now stated at a root with exactly `walkTree_book`'s hypotheses and is
`walkTree_guards' ∘ walkTree_book` (the admitted `EarSpec` form, with only `WF`/empty
tstack/`nxtEdgeIdx`/`firstOccurrence` size, is gone; `walkTree_local` takes `GuardsTree` as a
hypothesis). `earOut_one_entry` (`EarSpec`, over the fuel-indexed `earOut` with no walk-state
hypotheses) is replaced by `finishEdge_one_entry` (`EarShape.lean`, standard axioms): at a
`finishEdge` site with `EarAt curV d o origTstack true`, a returning **type-1** edge (`lowval < d`;
a tree edge closing the chain, or a back edge) leaves exactly one entry `(curV, lowval)` above the
`base` (`tstack.drop (length - origTstack)`), merged into the entry below exactly when that is the
P-ear `(curV, lowval)`. Two statement changes: the site-level form (the `EarAt` of the site is what
`EarTree` supplies, and `base` replaces the start tstack — the child's walk leaves the entries below
untouched only through the frame, not through the contract), and the `isType1` hypothesis: for a
type-2 edge the collapsed entry's `topDepth` is the minimum over the child's `mid` entries, and no
contract field bounds those below `lowval` (true for the DFS's `low`, but not derivable from
`EarFinish`). The proof is the run of `closeVert'`/`finishRest` on the `EarClose` shape
`c :: [py, vy] ++ base` (`mid = []` for type 1, `lowval ≤ c.topDepth`, `py.topDepth = lowval`,
`d < vy.topDepth`), via `maybeUnwrapNxt_run`/`merge_run_shape`/`retarget_run_eq`/
`finishTstackTop_run` and `finishRest_one_entry` (the P-check).

**The `EarTree` induction step (next, not started).** `walkTree_ear : EarTree t 0 s` asks for
`EarAt curV d o origTstack hasVert` at *every* `finishEdge` site, so the induction over `outs`
needs a between-edges invariant `EarCtx v d done base s` (the state after `walkOut v d oᵢ` for the
finished outs `done`, before `walkOutPre` of the next) from which the next site's `EarAt` follows,
and which `walkOut` re-establishes. From the consumers above and the proved pieces its content
is: (i) the stack is `top ++ base` with `base` the entries at entry to `v` (`FinishBook.origTstack`
frame: `walkTree_frame`/`Sim.walkTree` give that `base` is untouched by the child's walk, and the
`Keep` frame gives `stackVerts[k]`, `k ≤ d`, `firstOccurrence` below `d`); (ii) `top` is one entry
per finished out of `v` with `lowval < d` (`finishEdge_one_entry` for type 1 — the type-2 case
needs the DFS `low` fact `lowval ≤ topDepth` of every entry the child leaves, which is `Inv1`'s
`topDepth` bound through loops 1–2 and should be added to `EarClose` as `mid_top : ∀ t ∈ mid,
lowval ≤ t.topDepth`, 0 violations expected), each `(v, lowvalᵢ)`, pairwise edge/span-disjoint,
holding exactly that out's `subEdges` (this is `sub_edges`/`sub_cover`/`disj`/`span_disj` of the
*parent's* `EarAt` restricted to the entries above the parent's chain), plus, once `hasVert`, the
vertex entry `V v` at the bottom of `top` (`vert`/`vert_free`/`v_root`, `vy_top`); (iii) the child
site's `base`-side fields are the parent's: `p_entry` (the P-ear `(v, lowval)` below the top is the
single-piece entry left by the *first* out with that lowval — `bottom.py` for the first child,
`finishEdge_one_entry` after), `base_root`/`q_root`/`q_free` (freshness of the tree's own items,
`walkTree_book`'s `hefresh`/`hvfresh` through `Keep`), `touch_bot` (every entry of `top` touches
`v`), `base_disj` (nothing in `base` holds a `subEdges` edge: `t.edges.Nodup` + ownership); (iv)
the loop-1/loop-2 contracts `loop1`/`late` at the next tree edge are about the *child's* entries
after the child's walk, i.e. the child's own `EarCtx` at the end of its outs plus the first-child
chain shape (`bottom`: the chain bottom's `V y`/`(y, lowval)` pair is never merged until the ear's
top), which is the part with no proved piece yet; `EarCheck`'s `specWalk`/`lateWalk`/`closeCheck`
check exactly these at the site, so the first step is to add `EarCtx` to `checks/EarCheck.lean`
between the outs (one check per `walkOut` return) and run seeds 0..3000 before stating it in Lean.
The leaf of the induction (back edge, `sub = []`) needs only (i)–(iii).

**`EarCtx` checked (ear session 7, `EarCtx.lean`, `checks/EarCheck.lean` `ctxCheck`).** The
between-edges invariant is now stated as `WalkState.EarCtx v d done rest hasVert base bE sv s`
(statement only, no proof) and every clause is checked at the start of `walkOuts` and after every
`walkOut` return: 0 violations on seeds 0..3000, `ternarize = false` and `true`. The dump corrected
(ii) above in three ways. (a) The per-out entries are the entries *above* the one holding `V v`,
not all of `top`: the first out of `v` is the chain child, whose pieces stay as separate open
entries with `vStart ≠ v` (seed 1, `v = 5`, `d = 5`: `[V 5, (6,5), (6,3), (6,1), V 6] ++ base`).
(b) `V v` is a single entry at depth `d` only until a type-2 out closes with the vertex entry; after
that the entry holding `vertItem v` is the merged one with a chain vertex as `vStart` and
`topDepth ≤ d` (seed 1, `v = 4`, `d = 3`: `(1,1,[29,7,15,23,30,5])` holds `V 4`), and `hasVert`
can also be set by `finishEdge` itself (seed 1, `v = 5`: `V 5` pushed by a boundary finish), so
`done` carries per out whether the vertex entry existed at its `finishEdge`. (c) Only type-1 outs
P-merge: a type-2 out at the same lowval leaves a second, multi-item `(v, l)` entry (seed 1,
`v = 1`, `d = 4`: `(1,1,[29,7,15,23])` above `(1,1,[16])`), so "one entry per lowval, single root
item" holds only at lowvals all of whose post-push outs are type 1 (`allType1`/`CtxSingle`);
in general the entries above `V v` have non-increasing `topDepth` and strictly decreasing
`firstIdx` top-down, touch `v` and `stackVerts[topDepth]`, lie on side `stackDir[topDepth]`, are
attached only at `v` and `stackVerts[k]`, `topDepth ≤ k ≤ d`, and hold exactly the post-push outs'
`subEdges` (up to the `V v` entry). The header of `EarCtx.lean` lists which `EarFinish` field each
clause is meant to re-establish; the fields it does not cover (`loop1*`, `late*`, `loops`, `close`,
`bottom`, `boundary`/`bd_*`, `lower`, `dir_d`) need the child's end-of-outs `EarCtx` plus the chain
anchor.

**`EarAt` from `EarCtx` (ear session 8, `EarCtxAt.lean`).** `EarCtx` gained the clauses the site
proofs need, each dump-checked 0..3000 both modes before being stated: `CtxTop` (the stack above
`base` split at the entry holding `V v` when `hasVert`, type-1 ownership above it, `cover` of the
done outs' `subEdges`, `bot` freshness below it, `vitems`/`qitems` item coverage), the direction
list `sd`, `hv_ret`, `vert_book`/`vert_disj` (with `hasVert = false` the edges below `vertItem v`
are the returning done outs' edges, connected, two-attached at `v`, and disjoint from every
entry), `touch_bot`, `span_root`, and `v_fresh` extended to the stack arrays and `Touches`. One
first version of `hv_ret` (`hasVert = true → ∃ returning done out`) was reported with 11,119
violations and weakened by `rest ≠ []`; the violations were an artefact of the instrumentation,
which checks the state after `walkTree`'s own end-of-outs vertex push as `hasVert = true` (no
`finishEdge` ran there). With that site marked (`endPush`) the unconditional `hv_ret` has 0
violations and is the clause used (the bridge case of `ctx_step_tree` needs it: the child's
end-of-outs `hasVert` is `false` when all its outs are boundary). `CtxTop.ret` (`top ≠ [] →`
some done out returns; `ctx_top_ret`, 0 violations with the same `endPush` guard) gives the bridge
case the child's empty `top`; the stronger "every `top` entry without `V v` has an edge" is false
(`ctx_top_edge`: a returning child's own vertex entry `V y` can survive in the parent's `top` with
no edge). `TreeSite` carries, besides the child's end-of-outs `EarCtx` and frame, `items_kept`,
`hdir'` (the end push has `dir' = true`), `bridge_bd`/`comp_ret` (the child's done outs are all
boundary for a bridge, and return exactly to `d` for a component; `EarDfs.bridge_outs_bd`/
`comp_outs_ret` from `WF`), `e_ends`, `verts_lt`, `hsz`, `rest_nd`/`rest_nv`, `c_reach`/`edges_lt`
(the child's edges are connected to `y` inside any edge set containing them; `out_reach`/
`outs_reach` from `Ends`). `ctx_step_tree` is a proved case split (`H.tree`/`H.cls_ret`) over
`ctx_step_tree_bridge` (proved), `ctx_step_tree_comp` (proved) and `ctx_step_tree_ret` (proved
up to the named admission `tree_ret_shape`, below). `tree_comp_shape` (dump-checked
`compEndCheck`, 0 violations 0..3000 both modes; proved, standard axioms) is the child's
end-of-outs stack at a component edge: `hasVert = true` (no end push) and exactly
`t₁ :: ⟨y, d+1, _, ([], [V y])⟩` above the parent's `tstack`, with `t₁.vStart = y`,
`t₁.topDepth = d`, `t₁.spans.2 = []`. It is derived (`tree_comp_shape_of_shape`) from the child's
end-of-outs `EarCtx` (`CtxTop.split`) and a second between-edges invariant `CtxShape v d done
hasVert base s` (`EarCtxAt.lean`), carried by the same `cTree`/`cOuts`/`cOut` induction as an
additive conclusion (`COuts`/`COut` yield it after every out, `CTree` yields `TreeEndS`, i.e.
`TreeEnd` plus the child's end-of-outs `CtxShape`; `TreeSite.shape'`, `TreeSite.comp_ex`/`comp_t1`
from `WF`): `ret_hv` (a returning done out forces `hasVert = true`), `t1_flag` (a type-1 returning
done out was finished with `hasVert = true`), `t1` (if `hasVert` and every returning done out is
type 1, the stack is `above ++ vt :: base` with `vt` the `V v` entry and `above` strictly
decreasing in `topDepth`), `vt_side` (an entry of `v` holding `V v` has spans
`setSides (!stackDir[lowval]) [V v] []` for some returning done out). `CtxShape` is re-established
by `ctxShape_init` (empty `done`), `ctxShape_same` (boundary outs: stack unchanged),
`ctxShape_pushBack` (returning back edge: `Q e` push, optional `V v` push), `ctxShape_mergeP` (the
P-merge replaces the `(v, lowval)` entry by a single-item entry) and `ctxShape_ret` (returning tree
edge: `RetTop` + the optional `V v` push/merge), all standard axioms; the proved `ctx_step_*`
lemmas return it alongside `EarCtx`/`OutFrame`. From the shape
`ctx_step_tree_comp` is `earCtx_comp` (standard axioms) over the concrete post-state: the two
entries are popped, `Q e` gains `t₁.spans.1 ++ [V y]` as children and goes under `V v`; the
child's edges are exactly those below `t₁.spans.1` or `V y` (`CtxTop.edges`/`cover` of the child's
context, whose `top` is forced to be `[t₁, V y]`, plus `vert_edges`), `v_root`/`q_fresh`/`v_fresh`
use `CtxTop.vitems`/`qitems` of `t₁`, and `vert_book`'s connectivity of the child's edges to `y`
is `c_reach`. `ctx_step_tree_bridge` (standard axioms) is
`earCtx_bridge` over the concrete `finishBoundary` post-state: `bridge_bd` makes the child's
`hasVert` false and (via `CtxTop.ret`) its `top` empty, so `sE.tstack = s.tstack` and the end push
`V y` is the only entry popped; the items gain `Q e` under `V v` and the fresh `I` item and `V y`
under `Q e`. `below_kept` (no axioms) transports `Items.Below` from any item untouched by the
child's walk (`items_kept`: not a child vertex/edge item, allocated before) — every span item of
the parent's `tstack`, `V v` and `Q e` qualify — so the parent's edge sets, `touch_bot`, `disj`,
`q_fresh`/`v_fresh` are unchanged; `vert_book` is rebuilt from the parent's book, the child's book
(`ctx'.vert_book`, whose edges are exactly the child's edges by `bridge_bd`) and the bridge `e`
(`e_ends`), with `ends`/`comp` for `TwoAttached`. `ctx_step_tree_ret` (returning child,
`lowval < d`) is reduced operationally: `wp_finishEdge_tree` (proved) rewrites `finishEdge` at a
tree edge to `finishRest` (the P-check `finishP`, then `finishTail`'s first-edge vertex push) over
the loop-1–3 state `feS₃`/`feS₂` (`modifyVs`, `closeEars`, `mergeLate`, `closeVert'`). The stack
after loops 1–3 is the named admission `tree_ret_shape`, stated as the structure `RetTop v d s e
cls y outs L hv sX R` (`EarCtxAt.lean`; dump-checked literally by `retCheck` in `EarCheck.lean`,
0 violations 0..3000 both modes): relative to the parent's pre-push state `s` the post-loop state
`sX` has `tstack = R ++ L ++ s.tstack`, `g`/`stackVerts ≤ d`/`stackDir < d` unchanged, items only
grown, every item that is not a child vertex/edge item or `Q e` keeping its children and parents
(`kept`), `ch_lt`, and the retained entries `R`: parentless span items (`span_root`) that are fresh
or child/`Q e` items (`span_new`), `touch_bot`, edges exactly `subEdges o` (`edges`), touching
only `v`, the child's vertices or `stackVerts[lowval..d]` (`touch`), pairwise edge/span disjoint,
`vStart ∈ {v} ∪ child.verts`; with `hv = true` exactly one entry `⟨v, lowval, …⟩` carrying `V v`
as its single parentless item on the active side (`vert`); with `hv = false` no entry starts at
`v` and the entry holding `e` is the head (`noVert`). An earlier `retCheck` asserting a single
retained entry in the `hv = false` case was false (type-2 children keep several) and was dropped.
From `RetTop`, the `finishRest` cases are: `hv = false` — `cls.isType1 = false` (`hpush`), so
`condP` splits on `result`: P-merge of the fresh `V v` entry into the head of `R`
(`TEntry.mergeInto`) or a plain push of `V v`; `hv = true` — `result condP` false leaves `R`, true
is `ctx_step_tree_ret_P` (proved, standard axioms): `condP` forces `push = false` (a pushed `V v`
entry would sit at depth `d ≠ lowval`), so `hasVert = true` and `L = []`, the retained `(v, lowval)`
entry carries one item (`RetTop.vert`, type 1), and the entry below it is the enclosing parallel
ear; `earCtx_ret` with `T = R` gives the context at `sX`, `mergeP_single` the merged side and
`earCtx_mergeP` the post-state — `earCtx_mergeP` was generalised for this from the `Q e` entry to
any single-item `(v, lowval)` entry of `above` (its touch facts come from the entry's `CtxEntry`/
`CtxSingle`; the `Q e`-specific hypotheses `dest_edges`-style `hend`/`hdest`/`hqch` are gone).
The other cases end in `earCtx_ret` (proved, standard axioms): from `TreeSite`,
`RetTop` and the final top `T` (`R`, `V v :: R`, or `mergeInto (V v) c :: R'`) it rebuilds the
parent's `EarCtx` with `hasVert = true` and `o` appended to `done` — the old entries keep their
edge sets via `below_kept` over `RetTop.kept` (every span item of `L ++ s.tstack`, `V v`, and the
remaining outs' items qualify), `CtxTop.split` puts the retained vertex entry (or the fresh/merged
`V v` entry) as the new `vt` with the old `above` (ordered by `lowval_le_of_rank` + `vfirst`), the
old `V v` ownership (`vert_edges`) is unchanged because a returning out adds nothing below `V v`,
and `v_fresh`'s touch clause for the remaining subtrees uses the new `TreeSite.done_rest` (the
done outs' edges avoid the remaining subtrees' vertices; from `endsOut_wf` + nodup in `cOut`, like
`ctx_init_child`'s `hnc`). One
`EarFinish` statement change: `dest_edges` (every edge at the child's root is in `subEdges o`) is not
derivable from `walkTree_ear`'s hypotheses (graph edges outside the DFS tree are allowed) and was
only used for its touch consequence; it is replaced by `base_touch` (no enclosing `base` entry
touches `o.dest`), which `ear_boundary` consumes directly. `earAt_back_of_ctx` derives every
`EarFinish` field at a back-edge site from `EarCtx` of the parent plus the WF side conditions
(rank order of `done`, `Inc` of the done edges, nodup, the back class' lowval `< d`) and the exact
`walkOutPre` stack shape (`L`/`push`); no admissions (standard axioms). `earAt_tree_of_ctx` does the
same at a tree-edge site from `TreeSite`: the parent `EarCtx` before the site, the child's
end-of-outs `EarCtx y (d+1) done' [] hv' (L ++ s.tstack) …` at `sE`, the frame of the child's walk
(`g`, `stackVerts`/`stackDir` below `d + 1`, the two foreign roots `edgeItem e`/`vertItem v`), and
the child's end push (`L'`/`push'`/`D₃`, `pushEnd`). Derived: `tstack`, `back_nil`, `base_bot`,
`sub_bot`, `path`, `disj`, `span_disj`, `sub_edges`, `base_disj`, `sub_cover`, `touch_bot`, `vert`,
`vert_free`, `q_free`, `q_root`, `v_root`, `p_entry`, `sv_d`, `sv_child`, `path_child`, `dir_d`,
`boundary`, `base_touch`, `bd_noVert`. Named admissions (each the field verbatim over the site's
`sub`, `EarCtxAt.lean`): `earAt_tree_loop1_side`, `earAt_tree_loop1_touch`, `earAt_tree_loop1`,
`earAt_tree_bottom`, `earAt_tree_loops`, `earAt_tree_late`, `earAt_tree_late_fo`,
`earAt_tree_close`, `earAt_tree_bd_bridge`, `earAt_tree_bd_comp`, `earAt_tree_bd_term`,
`earAt_tree_bd_side`, `earAt_tree_lower` — the chain-anchor clauses they need are not yet in
`EarCtx`. Two `TreeSite` statement corrections while assembling: `ends` is guarded by
`d ≤ cls.lowval d` (only a boundary out has every subtree edge's endpoints in `v :: child.verts`;
for a returning out a back edge leaves the subtree — `ends` is consumed only by `boundary`), and
`hpush` records the actual push condition of `walkOutPre`, `hasVert = false ∧ lowval < d ∧
isType1` (the type-2 case pushes no vertex entry). `walkTree_ear` is now a **proved assembly**
(`WalkInv.lean`, `cTree`/`cOuts`/`cOut`, §4.3) over `earAt_back_of_ctx`/`earAt_tree_of_ctx`,
the proved entry contexts `ctx_init_root` (root: empty stack, `done = []`) and `ctx_init_child`
(child: base `L ++ s.tstack`, from the parent's `EarCtx`; its extra hypotheses `hinc`/`hnc`/`hndc`
— the done outs' edges are incident to `v` and avoid the child's vertices, the child's vertices are
nodup — come from `endsOut_wf` and the tree's nodup hypotheses in `cOut`), the back-edge step
`ctx_step_back` (proved, standard axioms), the tree step `ctx_step_tree` (proved over
`tree_ret_shape`; `tree_comp_shape` is proved via `CtxShape`), and the named admission
`tree_ret_shape` (dump-checked 0..3000 both modes: `ctxCheck` at every `walkOuts` step,
`retCheck`). The item frame of a child walk — formerly the admission
`walkTree_items_kept`, dump-checked `kept_items_ch`/`kept_items_par` — is a conclusion of the same
induction: `CTree` additionally yields `KeptRel g s t.verts t.edges s'` (a child walk changes
neither the children nor the parents of any allocated item other than its own vertices'/edges'
items) together with the freshness of the entries it leaves above `base` (every item they span is
a `V`/`Q` item of the tree or allocated during the walk; dump-checked `freshCheck`); `COut` yields
`OutFrame v o base s s'` per out (`ctx_step_back_frame`; `ctx_step_tree_frame` =
`outFrame_bridge`/`outFrame_comp`/`outFrame_retTop` + `outFrame_mergeP` over
`TreeSite.frame_child`, which combines the child's `KeptRel`/freshness — the additive `TreeSite`
fields `items_kept`/`child_fresh` — with the `L`/`L'` pushes via `outFrame_child`; the untouched
items of an out are those outside `V v`, the out's items and the inherited top entries' spans,
`Untouched`); `COuts` composes the frames (`outFrame_step`: an out's frame preserves
`KeptRel`/`Fresh` relative to the state at entry of the enclosing `walkTree`, since an untouched
item is neither a tree item nor spanned by a fresh entry; dump-checked `siteKeptCheck`).
`walkTree_qroot_kept`/`walkTree_vroot_kept` are corollaries of the kept frame, and
`walkTree_below_kept` (the edge set below `V v`, `v` outside the subtree, is unchanged;
dump-checked `kept_v_below`) follows by `below_kept` — the items below `V v` before the walk are
allocated (`ch_lt`) and are not the tree's items, which are roots
(`EarCtx.v_fresh`/`q_fresh` in `cOut`). `ends_of_wf_boundary` (dump-checked `kept_ends`) is proved
(standard axioms): `bdOut_wf`/`bdOuts_wf` — below the ancestor path `anc ++ pre`, if every return
depth of an out is `≥ anc.length`, its edges end in `pre`, at the vertex or in the out's vertices
(a back edge's destination is `(anc ++ pre ++ [w])[lowval]`, `WF`); at a boundary tree edge the
lowval is the least return depth (`lowval_eq_lmin`), so `pre = []`. `ctx_step_back` splits on `d ≤ lowval`: the boundary (self-loop) case
`ctx_step_back_boundary` is proved (standard axioms) — `hasVert = false` there (`hv_ret` + rank
order), nothing is pushed, `WF` forces `dest = v` (`classify` at index `d`), and `finishBoundary`
only relinks items, which `earCtx_selfLoop` states over an abstract post-state (the items gain
`edgeItem e` under `vertItem v` and the fresh `O` item `s.items.size` under `edgeItem e`; every
other field is transported through `IsParent`/`Below` characterisations); the returning case
(`finishBack`: `Q` push, P-check, first-edge vertex push) is the named admission
`ctx_step_back_ret`. For this, `EarCtx` gained the two allocation clauses `span_lt` (span items are
`< items.size`) and `ch_lt` (children are `< items.size`; `Shape.ch_lt` at the root), both
dump-checked first, and `ctx_step_back` takes the WF side conditions `e < ne`, `1 + nv + ne ≤
items.size`, `PairEq (v, dest) edges[e]`, `d ≤ lowval → dest = v`, and the bounds/freshness of the
remaining outs' edges and vertices (`hrest_e`/`hrest_v`), all discharged in `cOut`.

The returning case is now split at the P-check. `earCtx_pushBack` (standard axioms) is the
abstract post-state lemma for the `Q e` push: from `EarCtx` before the out and the equations of the
state after `pushEdgeTstack`/`firstOccurrence.modify`/`nxtEdgeIdx + 1` (graph, `stackVerts`,
`stackDir` below `d`, `items.size`, `Items.ch` unchanged; `tstack = Q e :: L ++ tstack` with `L`
the first-edge `V v` push iff `hasVert = false`) it re-establishes
`EarCtx v d (done ++ [(o, true)]) rest true`. For it `EarCtx` gained two more dump-checked clauses:
`noVert_after` (`hasVert = false → afterVert done = []`) and `vfirst` (an entry started at `v` that
does not carry `V v` has `firstIdx < nxtEdgeIdx`; the unrestricted `firstIdx < nxtEdgeIdx` is
false, 24462 violations). `ctx_step_back_ret` is proved from it: `finishEdge` unfolds to
`finishRest` over that state (`ctx_step_back_ret_rest`), `by_cases` on `result condP`; when the
P-check is false `finishP` is the identity and `finishTail` returns `hasVert || push = true`
(`isType1 = true` for every returning back edge), so the goal is `earCtx_pushBack`. When `condP`
holds (the entry below `Q e` is the parallel ear `(v, lowval)`), `ctx_step_back_ret_P` (standard
axioms) takes the `EarCtx` established by `earCtx_pushBack` to the context after the P-merge:
`condP` forces `push = false` (a pushed `V v` entry would sit below `Q e` at depth `d ≠ lowval`),
`mergeP_single` shows the `(v, lowval)` entry's active side is a single parentless item `i`
(`CtxTop`: a type-1 ear's entry carries one item; `allType1` over the done outs via the rank
order), and `earCtx_mergeP` is the abstract post-state lemma for the merge itself: `Q e` and the
`(v, lowval)` entry are replaced by one entry `(v, lowval, firstIdx_a, [p])` where `p` is either
the fresh `items.size` (ternarized, or `i` not a `P` item; children `[Q e, i]`) or the reused `P`
item `i` (children `Q e :: ch i`); every clause is transported through the `IsParent`
characterisation `IsParent s' p' c' ↔ IsParent s p' c' ∨ (p' = p ∧ c' ∈ S)` and
`below_addChildren` (`Below s' x i ↔ Below s x i ∨ (Below s x p ∧ ∃ r ∈ S, Below s r i)`). The
three `maybeUnwrapNxt` branches are the two instances of `earCtx_mergeP`. `ctx_step_back` takes
two more side conditions (`stackVerts[lowval] = dest` for returning outs, via `DfsOut.WF` +
`classify_eq_ret_iff` + the ancestor/`stackVerts` correspondence; done edges' endpoints avoid the
remaining tree outs' subtrees, via `endsOut_wf` + the tree's nodup vertices), both discharged in
`cOut`.

**`walkTree_ear`/`walkTree_book` are false without edge completeness (`EarFalse.lean`).** Trying
`ctx_step_tree` for a boundary tree out exposed that the parent's `vert_book` cannot be
re-established: the hypotheses of `walkTree_ear` (`t.WF []`, `t.Ends g`, bounds, nodup, the start
state, freshness) say nothing about graph edges that are incident to a vertex of `t` but are not
edges of `t` (`WF` only classifies the tree's own edges), and `VertBook`'s `TwoAttached` is then
false. Kernel-checked counterexample `walkTree_ear_false`/`walkTree_book_false` (standard axioms;
`decide +kernel` runs the walk on the concrete state): the path `0 - 1 - 2` (two bridges) with an
extra edge `2 - 3` outside the tree, `nv = 4`; after the bridge `1 - 2` is finished the edge set
below `V 1` is `{1 - 2}`, attached at vertex `2` by the edge `2 - 3`, which is neither terminal
(`1`, `1`). The same gap is why `EarFinish.dest_edges` had to become `base_touch`. Statement
correction (applied): the per-tree edge completeness hypothesis `hcomp : ∀ e < ne, ∀ x, Inc e x →
x ∈ t.verts → e ∈ t.edges` on `walkTree_ear`, `walkTree_book`, `walkTree_guards`,
`RootState.book`/`RootState.step`/`rootsBook_of_state` (per tree of the forest), and the forest
form `hecov : ∀ e < ne, e ∈ forest.flatMap edges` (which `walk_tree` already took) on
`walk_sides`/`walk_full`/`walk_tstack_nil`/`walk_root_children` and, outside the ear files,
`RangesWalk.forest_ranges_of_cover`/`walk_rangesInv_of_cover` (the only other `RootState.book`
consumer; `WalkItemsWF.walk_rangesInv` supplies `hecov` from its edge permutation). `hecov` gives
`hcomp` for every tree of the forest (`comp_of_forest`, `EarWalk.lean`, standard axioms: an edge of
another tree has both endpoints in that tree by `Ends`/`WF` (`tree_inc_verts`), disjoint from
`t.verts` by `ForestOK.verts_nodup`). Inside the induction `CTree`/`COuts`/`COut` carry the
relative form "every edge incident to the subtree is a subtree edge or a `pe` edge, and `pe` edges
join ancestors" (`pe` = the parent tree edge for a child, `False` at the root), which `cOut` passes
to the child (its own edge `e` is the `pe` edge; done/rest subtrees' edges avoid the child by
`endsOut_wf` + nodup) and records in `TreeSite.comp` for `ctx_step_tree`. `walk_tree`'s statement
is unchanged.

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
(finishing a child settles its entries), `loop1_rBranch_content` (every `.R` iterate of Loop 1 is
`RBranch ∧ RTop`; its stack shape `tstack`/`nxt_top`/`ne` is the proved `loop1_r_shape`, from
`run_loop1Cond` and `loop1Type_run`, and `loop1_rBranch` is the composite); `loop1_r_threeConnected`
combines the last with `RBranch.threeConnected`
(`closeEars_iter_step` gives `Inv (d+1)` at the iterate via `Step`; to be restated for `Inv'`).
`loop1_rBranch`'s hypotheses (`Inv D`, `Shape`, `CloseEarsOk`, `RInvAt`) say nothing about which
edges the Loop-1 entries hold, so `interior`, `proper`, `nxt_ne` and `mid` (after an S merge the
head's `vStart` is the S entry's, so the child is only interior to the head) are not derivable from them: its proof needs the Loop-1 ear
content (entries with `topDepth > d` at `closeEars` have `vStart = nxtV`; the entries with
`topDepth ≥ d` hold exactly the tree edge and the child's subtree edges) as an extra hypothesis
or from `EarShape`. The shape fields (`tstack`, `cur_top`, `nxt_top`, `ne`) follow from
`run_loop1Cond`, `loop1Type_run` and the head-`topDepth` induction over the iterates.
R-4 (`Proofs/RLoop1.lean`): the R-skeleton pipeline now reads the R branch through the ear
context. At the `finishEdge` site `FinishRShape.ear` (= `FinishBook.ear`) gives `EarFinish`, so
`L1Ctx.ofEar` + `l1_init` (`EarLoop1.lean`) give `L1Ctx`/`L1Inv` at `ceS₁`, and `l1_iter` carries
the `L1Piece` to every iterate `l1Iter d o s k`: `loop1_r_shape_ctx` (proved, standard axioms)
yields the shape `cur :: nxt :: rest` with `cur.topDepth = d` (`L1Piece.top` — the head-`topDepth`
induction is the ear layer's), `nxt.topDepth = d`, `nxt.vStart ≠ cur.vStart` and `cur` holding the
tree edge (`L1Piece.edges`, `l1Edges`). R-5: **`RBranch.cur_c` (`cur.vStart = stackVerts[d+1]`,
the piece bottom is the child) was false.** Kernel-checked counterexample
`checks/RBranchCounter.lean` (gen.py seed 2895, `ternarize = false`, the walk's own prefix as in
`RFinishEdgeCounter`): finishing the tree edge `13 : 4 → 0` at depth `1`, iterate 0 of Loop 1
S-merges the piece `{13}` with the child's `(0, 2)` entry, moving the bottom to `5 =
stackVerts[3]`, and iterate 1 is `.R` with head `(5, 1)` while `stackVerts[2] = 0`
(`counter : ¬ it.RBranch 1 cur nxt rest`, standard axioms). Restated minimally as the
`FinishTopOk.mid`/`TwoAttached.of_term` shape `mid : stackVerts[d+1] = cur.vStart ∨ Interior
(cur.edges) stackVerts[d+1]` (the child's edges are all in the piece once the bottom has moved);
`RBranch.hmid` and the `hmid` hypotheses of `RStep.cur_attached`/`union_attached`/`rCloseShape'`/
`threeConnected'` take the disjunctive form (an interior intermediate vertex cannot be an
attachment of `cur` or of `rU cur nxt`), nothing else consumed `cur_c`. Checked as `rbranch mid`
on seeds 0..400 × both modes + 6000 random multigraphs, 0 failures (the union form and `interior`
too). `interior` is proved from the ear context (`loop1_r_interior_ctx`, standard axioms):
`Loop1Spec` gives `L1Close` for `nxt` (the next entry of `hi`, at depth `d`; `rest = []` is
excluded by `lo_top`/`lo_ne`); its `bottom` puts the piece bottom `l1Bot o done` at `nxt.vStart`
(excluded by `ne`), at `stackVerts[d]` (excluded: `hi` entries and the child do not start at
`curV` — `EarFinish.sub_bot`/`path_child`, passed down as the hypotheses `hhb`/`hdb` of
`loop1_r_interior_ctx`/`loop1_rBranch_content_ctx`/`loop1_rBranch_ctx`/`loop1_r_keepsR'`), at the
child (then `L1Close.mid`, whose first disjunct is `nxt.vStart = cur.vStart` and whose third is
refuted by the tree edge `o.e ∈ cur`) or interior to the union. `loop1_rBranch_content_ctx` is
proved from three narrower named admissions: `loop1_rBranch_mid_ctx` (`mid`: `L1Close.mid` gives
the child's interiority only for `rU cur nxt`; the `cur`-alone form is R content — the child's
back edges to `stackVerts[d]` are P-merged into the piece), `loop1_rBranch_fields_ctx` (`RBranchFields`: `cur_piece`, `cur_vs`, `nxt_ne`, `proper`,
`nxt_touch_*`, `nxt_no_cu` — `RangesInv`/`RunSaturation` content at the iterate) and
`loop1_rTop_ctx` (`RTop`: the settled-frontier/saturation content). `loop1_rBranch_content_ctx`
also gained the hypothesis `hv : v₀ < s.g.nv` (needed by `l1_iter`; available at its only call
site `loop1_rBranch_ctx`). `loop1_rBranch_ctx` composes them and `loop1_r_keepsR'`
(`Proofs/RSkelWalk.lean`) replaces `loop1_r_keepsR` in `keepsR_finishEdge_site`. The
ear-context-free `loop1_rBranch_content`/`loop1_rBranch` keep their statements and are no longer on
the pipeline.

Item level (`Proofs/RItems.lean`). `Items.RSkel3 g items i`: the skeleton of the R item `i` (its
non-`V` children as pieces, the complement of its edge set as the parent piece at `i`'s terminals,
contracted) is 3-connected. `RBranch.rSkel3` proves it for the item Loop 1's R branch closes
(`rCloseItems`: `allocItem .R`, merge, `finishTstackTop`), from `RBranch.threeConnected`
(under `Inv' (d+1)`) and `Shape`. Admitted: `items_r_three_connected` (every R item of `g.walk` on a block is `RSkel3`;
needs `loop1_rBranch`, the same argument for the type-1 R close of `finishEdge`, and that closed
items are never modified again). For a graph with several blocks the parent piece must be the
complement within the item's block; `RSkel3` with the whole graph is the block statement only.

Skeleton persistence frames are proved in `Proofs/RItems.lean`: `Pieces.contract_congr` reads
only valid edge/piece indices; `Items.RSkel3.congr` needs only the item's children, their types
and non-V terminals/edge sets, and the item's own terminals/edge set.
`Items.RSkel3.modify_of_not_below` preserves a closed skeleton when another subtree is modified;
`Items.RSkel3.push_nil` preserves it when a childless item is allocated.

R-5 (`Proofs/RVert.lean`): the type-1 vertex close with `feSingle = false` builds an R item
(`maybeUnwrapNxt .R`, two merges, `retarget`, `finishTstackTop`). On the `EarClose` stack
`c :: py :: vy :: base` at `feS₂` its children are `vItems c py vy` (the concatenated spans), its
edge set `vU c py vy`, its terminals `curV` and `stackVerts[min(c, py, vy).topDepth] =
stackVerts[lv]` (`c_top`/`py_top`/`vy_top`). `vClose_rSkel3` (proved, standard axioms) turns an
`RCloseShape` for these data into `Items.RSkel3` of the built item, mirroring `RBranch.rSkel3`
(`cvS₅_type1` computes the state, `getSide_vMerged` reads the children from `FinishTopOk.side`);
`closeVert_type1_rSkel3` (`Proofs/RSkelWalk.lean`) is proved from it. The remaining content is the
named admission `closeVert_type1_rCloseShape`: the `RCloseShape` itself (pieces pairwise
edge-disjoint and 2-attached from `EntryR` at `feS₂`, union = the sub-ear at `(curV,
stackVerts[lv])` from `EarClose.sub_edges`/`sub_cover`/`type1`, and the HT fields `single`/
`maximal`/`type1`/`bond`/`type2` from run saturation and `RangesInv`).
All four have only `propext`, `Classical.choice`, `Quot.sound` as axioms.
The walk-level ownership argument that all later writes satisfy these frames remains open.

**`dfsForest_rooted` (R-5, restated and proved).** As first stated —
`(DfsData.ofForest (g.dfsForest vo eo)).Rooted g` for 2-connected `g` — it is false:
`ofForest` takes `root` from the first tree and `dfsForest` roots a tree at every vertex in
`inOrder` order, so an isolated vertex ordered first is `root` with no out-edges
(`checks/RRootedCounter.lean`: `g = ⟨3, #[(1,2),(1,2),(1,2)]⟩`, `vo = eo = []`, kernel-checked,
standard axioms). The statement now is `∃ r, dfs.depth r = 0 ∧ ({ dfs with root := r }).Rooted g`
(re-rooted at the root of the tree holding the edges), proved from `Spec.edge_out`/`joins`/
`back_anc` (the endpoints of every edge are `Anc`-comparable), `Anc.comparable` (a depth-0
ancestor of one endpoint is an ancestor of the other) and 2-connectivity (every edge is joined to
edge `0` by a walk). `walk_rSkelInv` runs the `rk*` induction on the re-rooted data (`Spec` with
`depth_root := hr0`; all other fields are root-free).

R completeness must come from postorder interval ownership and saturation, not `EarFinish`
or Loop-1 ear content. The pending `RangesInv` interface must connect entries to the actual
processed prefix, DFS path, and completed vertices. At the `ceS₁` boundary it must account
for `mid` (the returned child is the bottom or interior to the piece), interior ownership of the returned child, disjointness, nonemptiness, and a
block-local complement. An additional saturation argument must derive `EntryR`'s `single`,
`maximal`, `bond`, `type1`, and `type2` fields; these are consumed by `RBranch.rContent`,
not consequences of edge-disjointness alone. Closed-piece maximality must survive later writes.

The necessary adjacent-merge witness lemma needs the path/vertex ownership information.
Edge intervals and pairwise adjacent saturation alone are insufficient: take the union of
the cliques on `{0,1,2,3}` and `{0,1,4,5}`, sharing edge `{0,1}`. Give each edge a singleton
interval. Every singleton is 2-attached, but every pair of distinct edges has at least three
boundary vertices (all vertex degrees are at least three), so every ordering is saturated
under a two-entry/2-attachment check. Nevertheless `{0,1}` separates the graph into the
classes through `{2,3}`, through `{4,5}`, and the singleton shared edge.
This is a counterexample to that abstract sufficiency claim, not a reachable walk state.
No new ear-content assumption is used to fill the gap; the interface restatement remains open.

Run saturation is now stated in `Proofs/RunSaturation.lean` as
`WalkState.Saturated dfs origTstack s`: no 2-attached run with at least two edge-bearing entries
below the schedule frontier; no proper 2-attached run of at least two non-V children of a
closed R item; and no proper 2-attached union of an R child with additional back edges at a
vertex. It is a hypothesis, not yet an invariant theorem. Empty entries do not count toward
the stack-run length. The closed-item clause excludes the whole child list, which is
2-attached at the cap even for K4. `checks/RunSaturationCheck.lean` exhibits the actual K4
output: the R item's five non-V children own `[2,4,5,1,3]` and attach at `{0,1}`. The DFS
postorder is `[1,2,4,5,3,0]`, so child order is not itself the postorder.

Proved, with standard axioms only: `Pieces.WF.sepClass_of_mem` and `.sepClass_constant`
(a skeleton pair other than a piece's terminal pair cannot split that piece),
`runEdges_of_interval` (selecting one marker per piece transports an aligned edge interval
to a run of pieces), `Graph.RunSaturated.short_or_eq` and `.laminar`,
`Pieces.WF.class_laminar_of_interval`, `Saturated.closed_laminar`, and
`entryLaminar_of_runSaturated`. These are conditional derivations, not proofs that the
walk maintains saturation. The class lemma requires the class to be covered by the listed
pieces and their marker order to be the corresponding postorder subsequence; neither that
alignment nor the crossing-class cases follows from the current range interface.
`EntryR.maximal`, `.bond`, and `.single` still need derivations, as does saturation inside
open entries (the stack-run and closed-R clauses alone do not state it).

Coverage and laminar containment range over valid graph edges. `Items.EdgeBelow` itself is
an ancestry predicate and does not enforce `e < g.ne`: an allocated item `i` contains the
out-of-range index `i - 1 - g.nv` naming itself. On valid K4, R item 11 contains index 6,
whereas none of its non-V children contains that index. `checks/REdgeDomainCheck.lean`
proves this failure of unrestricted child coverage, and proves coverage for all valid edges,
using only standard axioms. The run-to-item bridges now take bounded coverage and identify
classes with the valid-edge restriction of each run. `EntryLaminar` restricts its containment
clause; `RContent` and `RCloseShape` use `Pieces.LaminarWith` on `e < g.ne ∧ U e`.
The type-1 saturation clause's properness witness must also be a valid edge.
`runEdges_restrict` and `Graph.RunSaturated.restrict` justify the restriction. The affected
R-maximality, R-branch, persistence, and HT-to-cut transports remain standard-axiom proofs.

The child-return settling claim also fails on a reachable state.
`checks/RInvReturnCheck.lean` proves `dfsForest_eq` and `returned_not_rInvAt` with only
`propext`, `Classical.choice`, and `Quot.sound`. For the graph
`[(0,1), (1,2), (1,2), (1,2), (0,2)]`, its actual DFS visits `0 → 1 → 2`.
The check executes the two `walkOutPre` prefixes and then `walkTree` at vertex 2.
The returned stack contains `⟨2, 1, 1, ([], [9])⟩`; item 9 is a P item with
terminals `(2,1)` and children Q6/Q7 (edges 2 and 3).
Edges 1 and 0 are outside this item and are not in the same `{2,1}` separation class:
edge 1 has both endpoints deleted, whereas edge 0 belongs to the path through vertex 0.
Thus `returned.RInvAt dfs 1` is false for every `dfs`, before the pending tree edge 1
is processed. This checks the formal `Graph.SepClass`, not a Boolean approximation.
It refutes the intended reachable-state conclusion; the check does not package all the
antecedents of the admitted `walkTree_rInvAt` into a formal negation of that implication.
The return interface must distinguish the preserved base from the provisional frontier,
and settle the latter only after incorporating the pending tree edge. In particular,
`walkTree_rInvAt` cannot currently supply the parent-settled hypothesis consumed by
`finishEdge_rInvAt`; adding interval/run saturation facts cannot repair this timing issue.
The corrected preservation contracts remain open rather than assuming provisional
frontier pieces are maximal.

The schedule-specific input is `WalkState.Frontier (o := o) d origTstack s`
(`Proofs/RInvFrame.lean`), separate from interval ownership and saturation. Since the stack is
top-first and `origTstack` counts the preserved bottom entries, the split is at
`s.tstack.length - origTstack`, using `take` for the frontier and `drop` for the base.
Before `finishEdge`, frontier edges together with the pending `o.e` are exactly `subEdges o`;
the base owns none of these edges. At every reached iterate of loops 1–3, `FrontierOwns`
requires that same base suffix and exact subtree-plus-tree-edge ownership above it.
The loop fields also require enough entries above the boundary for every enabled merge
(three entries for Loop 1's S case, two otherwise). The saturation boundary is this
`origTstack` boundary, never a `topDepth ≤ d` cut.
`finishEdge_rInvAt` now takes this single hypothesis; `walkTree_rInvAt` takes `FrontiersTree`,
which threads it through the actual nested walk calls, like `GuardsTree`. Both proofs remain
admitted. Exporting the schedule fact from `EarFinish.bottom`/`sub_cover` belongs to the ear
layer; it has not been re-derived from stack shape here. `RInvAt`, `RTop`, and `loop1_rBranch`
still await the interval/saturation interface rather than hiding their content in `Frontier`.

The graph-theoretic HT-to-cut implication is proved in `Proofs/ThreeConnected.lean`.
`Graph.ThreeConnected.relabel` takes `TwoConnected`, HT `ThreeConnected`, and a nodup list
of at least four active vertices containing every edge endpoint, and proves the cut-based
`SpqrTree.ThreeConnected` after `idxOf` relabelling. It permits isolated original vertex IDs
outside that list and does not require simplicity. `TwoConnected.two_incident_of_vertices`
gives two distinct incident edges at every active vertex; disconnected surviving vertices
would therefore supply two separation classes of size at least two. All five supporting
theorems and the relabelling theorem have only standard axioms.
`Pieces.ofItems_addParent_edges` identifies the contracted edge list with the item endpoint
pairs followed by the parent pair, provided every edge of `U` belongs to an item in `L`.
`Items.rSkeleton_perm_contract` then identifies its relabelling with `Items.rSkeleton` up to
permutation. Its coverage hypothesis is precisely:
`∀ e, e < g.ne → items.EdgeBelow g i e → ∃ c ∈ items.ch i, items.type c ≠ .V ∧ items.EdgeBelow g c e`.
These two theorems also have only standard axioms.
`Graph.ThreeConnected.twoConnected` proves that HT 3-connectivity implies edge-based
2-connectivity when there are at least three edges. First, deleting the endpoints of a third
edge rules out disconnected edge components. If deleting one vertex separates two edges,
each resulting edge class must be a singleton: otherwise isolate an edge incident to the
deleted vertex and obtain three separation classes. A third edge then gives a contradiction.
The three-edge lower bound is necessary (two disjoint edges have no HT separation pair).
`Items.RSkel3.rThreeConnected` now gives the item-specific bridge from `RSkel3` to the
cut-based item skeleton specification, assuming `Items.WF` and the coverage clause above.
The skeleton's edge count supplies its 2-connectivity via `ThreeConnected.twoConnected`;
neither block 2-connectivity nor child-plus-parent piece WF is needed for this bridge.
It derives the active vertex list: endpoints are on the cap or a child edge; a V child has
an incident edge covered by a non-V child, and the `interior`/`separation` clauses force that
vertex to be a terminal of such a child. The bridge handles either cap orientation using
`SpqrTree.ThreeConnected_congr_undirected`; it does not require a simplicity hypothesis.
These lemmas have only standard axioms.

`Items.WF.q_root_covers_block` proves that, in a 2-connected graph, every Q item with
nonempty children contains every physical edge. Its single recorded endpoint and
`Endpoints.separation` make its edge set 1-attached; edge connectivity after deleting that
endpoint propagates membership from its own edge to every edge.
`vertex_child_leaf_of_block` then proves that V children of S/P/R nodes are leaves,
provided every Q child of such a V item has nonempty children: otherwise that Q, and hence
the V child, contains every edge, contradicting `Endpoints.interior` for its parent.
`nonV_child_cover_of_block` derives the required coverage from this leaf property and
item ancestry. `RSkel3.rThreeConnected_of_block` applies the connectivity bridge with this
derived coverage. All four lemmas have only standard axioms.
The remaining walk placement hypothesis is exactly
`∀ v, v < g.nv → ∀ c, items.IsParent (vertItem v) c → items.ch c ≠ []`.
It is not a field of `Items.WF`: `Shapes.q_leaf_of_node` only forces Q children of
non-F/V parents to be leaves. Exporting the placement fact from the walk, and the
block-local version for an arbitrary input graph, remain open.

`Ranges.q_root` and `Ranges.q_leaf` condition on nonempty/empty Q children without
constraining their parent; `Items.Typing.v_children` in `WalkTyping.lean` only supplies the Q type.
Missing range/typing clause (proposal, not an admission): `v_child_root : ∀ v c, v < g.nv → items.IsParent (vertItem v) c → items.ch c ≠ []`.

`WalkState.RReturn before after dfs parent` is the provisional child-return postcondition:
the bottom `before.tstack.length` entries equal the original stack; only entries in that
preserved base, except those starting at `parent`, must satisfy `EntryR`; the entire stack
remains edge-disjoint. The prefix above the base is provisional until the parent's pending
tree edge is incorporated by `finishEdge`. `WalkTreeRReturnSpec` names the preservation
obligation as a `Prop` definition, not a theorem/admission. It starts settled at the parent
and relates the walked out-list to `dfs.outs c`. The old `walkTree_rInvAt` statement and
`returned_not_rInvAt` counterexample are retained; no preservation proof is claimed.
**Depth bound on settling (`RInvTop`, `RInvFront`; `checks/RFinishEdgeCounter.lean`).**
The parent exemption alone is still too strong, in both directions. `Deep` (4 vertices,
`0-1-2-3`, `3→1` doubled, `3→0`): after the parent `2` (depth 2) processes the tree edge `2→3`,
the entry `(3, 1)` (bottom 3, top at depth 1) holds the P piece `{3,1}` of the two back edges
while the path class `1-2-3` is still a second `{3,1}` class — `finished_not_rInvAt` refutes the
settled-at-`curV` output of the old `finishEdge_rInvAt` although its base is settled
(`base_entryR`). `Base` (seed 3757 of the extra corpus, 5 used vertices): before the return from
`2` to its parent `3`, the grandparent's base entry `(1, 0)` holds the closed S piece `{6,1}` of
the first ear `1-4-6` while the tree edge `1→3` and the back edge `1→6` are outside it, so
`RReturn.entries` with only `vStart ≠ parent` exempt is false (`not_base_entries`). The
invariant the walk actually keeps is bounded by the top depth: an entry is settled only once the
walk is back at its top, i.e. `RInvTop s dfs v d` (= `RInvAt` restricted to `d ≤ t.topDepth`),
`RInvFront s dfs v d origTstack` (base below the split, same bound, plus whole-stack
disjointness) is the input of `finishEdge_rInvTop` (back-edge branch proved; tree-edge branch proved up to the admitted `finishEdge_tree_top_settled_first`; replacing `finishEdge_rInvAt`;
`walkTree_rInvAt` is deleted), `RReturn before after dfs parent d` carries the same bound, and
`WalkTreeRReturnSpec` starts from `RInvTop` at the parent. `loop1_rBranch` takes `RInvFront`.
`checks/RFinishEdgeCheck.lean` evaluates both contracts at every `finishEdge` site on blocks
(all six `EntryR` fields via an executable `SepClass`/`Type2Pair`): on seeds 0..300, the fixed
regression and 6000 extra random multigraphs, both ternarize modes (303 block graphs, 6414
sites), the old contract A fails on 146 runs and contract B, whole-stack disjointness and the
R-branch `RTop` check (34 R closes) fail on none. The R-side frame lemmas (`congr`, `of_eq`,
`setStackDir`, `modify`, `modifyItem_free`, `pushVert`, `toRTop`) are restated for `RInvTop`.
`checks/RInvReturnCheck.lean` kernel-checks the fixed example's preserved base and its
settled-entry clause. Its executable checks the base split, disjointness (including item IDs
beyond physical edges), and preservation of every state observation read by `EntryR` on
non-parent base entries, on both ternarize modes of seeds 0..400 from `gen.py`.
Base shape/disjointness and instrumentation agreement pass all 802 runs; the content-frame
check applies to the 91 graphs satisfying edge-based `TwoConnected` (182 runs), as required
by the contract. The fixed regression also passes. This checks conditional preservation,
not all `EntryR` premises at arbitrary walk states. Empty vertex entries have vacuous
`EntryR` content and need no terminal-value frame: stale depth slots may be overwritten
while walking another subtree. An initial stronger terminal-frame check flagged 78
seed/mode runs: 38 involved empty entries, and the other 40 were outside `TwoConnected`.
R-5: the terminal-value frame is observed only for base entries with `topDepth ≤ d` (the
parent's depth). Seed 390 (both modes) has non-empty buried base entries with `topDepth > d`
(single-`Q` entries `(3, 14)`/`(11, 13)` of finished deeper vertices, the `base_top`-false
shape) whose slot the child's walk rewrites; `RInvG`/`EntryR` still holds for them (their
fields are vacuous for one edge) — `checks/RFinishEdgeCheck.lean`'s contract-B lines check
that exact statement at every site of seeds 0..400 with 0 failures, so the positional
terminal frame is a sufficient condition that is false in general, not a contract failure.
Reproduce with `lake env lean --run checks/RInvReturnCheck.lean`, supplying the concatenated
outputs of `python3 ../gen.py 0` through `python3 ../gen.py 400`, prefixed by `401`.

**Exports needed from ear/Ranges (R-5).** The remaining R admissions take facts of the ear and
range layers as explicit hypotheses; this lists, per site, the exact statement each export must
have so ear-5 (`walkTree_ear`) and Ranges-3 can discharge them by composition. Sites are the
states of the `rsTree`/`rsOuts`/`rsOut` induction (`RSide.lean`), i.e. where the `GuardsTree`/
`BookTree` conjuncts are read: a *child-entry site* is the `walkOut v d o hasVert` call state `s`
(before `walkOutPre`) for a tree out `o = .tree e cls (.node c couts)`; a *vertex site* is a
`walkOut` call / end-of-`walkOuts` state with `hasVert = false`; a *finish site* is the
`finishEdge v d o origTstack hasVert` call state `s₃` (after the child's `walkTree`) with
`o.cls.lowval d < d` and `origTstack = B` the stack length before the child. `rl1Iter d o s k`
(`RLoop1.lean`) and `l1Iter d o s k` (`RangesCloseSites.lean`) have the same body.

Ear (ear-5) — as a `BookTree`-style predicate (`CtxTree`) carried by `walkTree_ear`:
- E1 (child-entry and vertex sites): `EarCtx v d done rest hasVert base bE sv s` for some
  `done rest base bE sv`, `rest = o :: _` (resp. `rest = []`). Fields read: `split`
  (`CtxEntry.depth`: `above` tops out `< d`; `vt.topDepth ≤ d`; `below` has `vStart ≠ v`),
  `base_bot`, `vert_edges`, `vert_free`, `v_root`, `disj`, `span_disj`, `afterVert_ret`. Derived
  from it by R-5: `∀ t ∈ s.tstack, t.vStart = v → t.topDepth ≤ d` (`rSide_entry_site`'s
  parent-start bound).
- E2 (new field candidate `bd_free`; for `rSide_vertFree_site` and `FinishRShape.vert_own`): at
  a vertex site, and at a finish site with `hasVert = false`,
  `∀ t ∈ s.tstack, ∀ e, e < s.g.ne → t.edges s.g s.items e → ¬ Items.EdgeBelow s.g s.items (vertItem v) e`
  (the blocks already closed at `v` are owned by no open entry; with `vert_edges` this is
  `∀ t ∈ s.tstack, ∀ e < ne, t.edges e → ∀ o ∈ done, d ≤ o.1.cls.lowval d → ¬ subEdges o.1 e`).
  The `e ≥ s.g.ne` residue of `VertFree` (item ids aliasing `edgeItem e`) is E4.
- E3 (field `stab_single`; for `rSide_entry_site`'s `stab`): at the child-entry site
  `walkTree c (d + 1)`, every open `t` with `t.topDepth = d + 1` (the only entries whose
  `EntryR.single` reads the updated `stackVerts[d + 1]`, `EntryR.congr_top`) satisfies
  `s.g.TwoAttached (t.edges s.g s.items) t.vStart c ∨ ∀ e e', e < s.g.ne → e' < s.g.ne → t.edges s.g s.items e → t.edges s.g s.items e' → s.g.SepClass t.vStart c e e'`,
  i.e. `EntryR.single` holds in the updated state outright. *Refuted original* (`buried_vacuous`):
  "every open `t` with `d < t.topDepth` is a single edge or `TwoAttached … vStart vStart`" (so that
  `single` is vacuous). `checks/WalkInvCheck.lean` (seeds 0..3000 × both modes) found 6 violations
  (seeds 666, 989, 1966), all buried P-bonds attached at `{vStart, stackVerts[topDepth]}` — two
  stale vertices of a finished sibling subtree (`RFinishEdgeCheck`'s 70 buried entries on seeds
  0..400 + 6000 random were all single-edge by luck). Kernel-checked counterexample
  `WalkInvCheck.E3False.e3_buried_vacuous_false` (`checks/WalkInvCheck/E3False.lean`, axioms
  `[propext, Quot.sound]`): the edge-minimal block of seed 1966 (8 vertices, 11 edges), the walk run
  by `decide +kernel` to the child-entry site of `7` under `1` (depth 2) after the subtree
  `1 → 4 → 6`; the open entry `⟨6, 3, 2, ([], [21])⟩` (`topDepth = 3 = d + 1`) owns the P item of the
  two parallel edges `4 – 6`, attached at `4 = stackVerts[3]` (stale) and `6 = vStart`. The restated
  field holds there because both edges share `4 ∉ {6, 7}`; it is what `stab` actually consumes
  (`EntryR.set_stackVerts` for `topDepth ≠ d + 1`, `stab_single` for `topDepth = d + 1`). Checked:
  `e3.stab_single` in `checks/WalkInvCheck.lean`, 0 violations on all graphs (not only blocks).
- E4 (item-forest facts at every site; ear or Ranges, whichever carries them — `ItemFree`/`cnt`
  are Ranges-3's): `∀ t ∈ s.tstack, ∀ i ∈ t.spans.1 ++ t.spans.2, ∀ p, ¬ Items.IsParent s.items p i`
  (span items are roots) and `∀ c p p', Items.IsParent s.items p c → Items.IsParent s.items p' c → p = p'`
  (unique parent). With `v_root`/`vert_free` these give the `e ≥ s.g.ne` part of `VertFree`/
  `vert_own` (an id below both a span item and `vertItem v` forces the span item under
  `vertItem v`, i.e. equal to it).
- E5 (new `Loop1Spec`/`L1Close` field candidate `mid_cur`; for `loop1_rBranch_mid_ctx`): at a
  reached split `L1Reach d hi done (t :: rest)` of a tree edge `o` with `lowval < d`, with `t`
  closing at depth `d`:
  `l1Bot o done = s.stackVerts[d + 1]! ∨ s.g.Interior (l1Edges o s done) s.stackVerts[d + 1]!`
  (the child is the piece's bottom or interior to the piece *alone*; `L1Close.mid` only says it
  for `l1Edges ∪ t.edges`). Already checked: `RFinishEdgeCheck` `rbranch mid` lines, 6518 sites,
  0 failures.
- No new ear field for `loop1_rTop_ctx`/`feS₂_top_entryR`/`settled`: they read `L1Ctx.hi_top`,
  `L1Inv` (`L1Piece`, `L1Keep.edges`), `EarFinish.loops`/`late`/`close`, all exported already,
  on top of R1/R2.

Ranges (Ranges-3) — as an `RgTree`-style predicate carried by the walk, or per-site records:
- R1 (finish site): `CloseCtx σ n (d + 1) v d o B hasVert s₃` for some `σ n` (whole record; fields
  read: `ranges`, `nodup`, `lt`, `pos`, `block`, `l1_site`, `v_site`, `p_site`, `finishR`, `close`,
  `dest_edge`). Consumers: `rSide_finish_content_site` (`settled`), `feS₂_top_entryR`,
  `loop1_rTop_ctx`, `loop1_rBranch_fields_ctx` (`cur_piece`/`cur_vs`/`proper`/`nxt_touch_*` from
  `VSite` at `l1_site`), `closeVert_type1_rCloseShape` (`v_site`).
- R2: `RangesInv σ n (d + 1)` at the states R-5 reasons at — `l1Iter d o s k` for every `k` with
  `∀ j ≤ k, result (loop1Cond d) (l1Iter d o s j) = true`, `feS₁ d o s`, `feS₂ d o s`,
  `feP v d o s` — as lemmas from `CloseCtx` (if `RangesInv.iter_mergeAdj`/`mergeAdj_of_frontier`/
  `RgStep.*` already give these, name the lemma; else export).
- R3 (child-entry and vertex sites): `RangesInv σ n d` at `s` (`processed`-style ownership for
  `rSide_vertFree_site`; E4 if span-rootness lives here).

Pure R content with no export (R-5 proves it from the above): the `entry_cur` saturation at an
R iterate, `feS₂_top_entryR`'s class saturation, `closeVert_type1_rCloseShape`'s HT fields,
`settled`.

**R site obligations (R-5 handoff to the backbone induction).** The remaining named R admissions,
each with the site it sits at and the `EarCtx`/`CtxEntry`/`EarFinish`/`RangesInv`/`CloseInv`
fields that discharge it once one induction carries them together with `RInvTop`/`RInvFront`/
`RSkelInv`. Sites as above (child-entry site = `walkOut v d o hasVert` call state; vertex site =
`walkOut` call / end of `walkOuts` with `hasVert = false`; finish site = `finishEdge v d o B hasVert`
call state after the child walk, `lowval < d`, `B = origTstack`). Every statement below was
dump-checked on seeds 0..400 × both modes + 6000 random multigraphs (`checks/RFinishEdgeCheck.lean`,
0 failures).

1. `feS₂_top_entryR` (`Proofs/RInvTree.lean`): `(feS₂ d o s).EntryR dfs c` for the top `c` of
   `feS₂ d o s` with `c.topDepth = d`, `c.vStart ≠ curV`, at a first-edge (`hasVert = false`)
   tree-edge finish site (`o.cls = .ret lv kind`, `kind ≠ .backEdge`, `lv < d`; hypotheses
   `Inv' D`, `Shape`, `FinishOk`, `FinishGuards`, `Frontier`, `FinishRShape`, `RInvFront dfs curV d
   origTstack`, chain/`Spec`/`Rooted`/`TwoConnected`). `c` is loop 1's output (loop 2 did not fire,
   `finishEdge_tree_top_settled_first`). Discharge: `pieces`/`maximal` from `RangesInv σ n (d + 1)`
   at the site transported to `feS₂` (`convex`/`closed`: each closed item owns a σ-interval, so its
   complement is one class) with `CloseInv.closed` at `feS₂` (`CloseAt.att_vs`/`vs_att`: the item's
   `vs` are its attachments) and `CloseCtx.l1_site k`'s `VSite` (`kinds`/`two`/`r_shape`/`p_shape`/
   `s_order`) for the piece shapes; `bond` from `EntryR.bond` of the consumed range entries
   (`RInvFront.base` for those above `origTstack` with `topDepth ≥ d`) plus the P-merge of the tree
   edge (`EarFinish.p_entry`); `single` from `EarFinish.loop1` (`Loop1Spec`: `L1Close.mid`, bottom
   at `nxt.vStart`) and `RangesInv.ordered` (every `(nxt.vStart, curV)` class edge lies in the
   consumed range); `type1`/`type2` from the consumed entries' `EntryR.type1`/`type2`
   (`EntryLaminar` is closed under the union of pieces sharing the terminal pair) — R content on top
   of those fields.
2. `loop1_rTop_ctx` (`Proofs/RLoop1.lean`): `(rl1Iter d o s k).RTop dfs cur nxt` at an `.R`
   iterate with `tstack = cur :: nxt :: rest`, `cur.topDepth = nxt.topDepth = d`,
   `nxt.vStart ≠ cur.vStart`, from `L1Ctx`/`L1Inv` (`ceS₁ o.dest d o.e (feS₀ d o s)`),
   `CloseEarsOk`, `RInvFront dfs stackVerts[d] d origTstack` at `feS₀`. `entry_nxt`: `nxt` is the
   next unconsumed entry of the range `hi` (`L1Inv`'s `L1Reach`, edges unchanged by `L1Keep.edges`);
   it tops out at `d` and starts below the child, so it is not in `RInvFront.base` — its `EntryR` is
   the child-return settledness of the frontier (`walkTree_rSide`'s `settled`, item 7), which the
   backbone must carry as `EntryR` for every entry with `topDepth = d`, `vStart ≠ curV` at
   `feS₀`. `disj`: `RInvFront.disj` at `feS₀` transported across the iterates (merges keep pairwise
   disjointness; the `e ≥ ne` part needs `L1Ctx.disj`/`span_disj`/`q_root`/`q_free` at item level).
   `entry_cur`: the merged piece's saturation — R content on `RangesInv.convex`/`ordered` at the
   iterate (R2) and `CloseCtx.l1_site k`'s `VSite`.
3. `loop1_rBranch_fields_ctx` (`Proofs/RLoop1.lean`): `(rl1Iter d o s k).RBranchFields d cur nxt`
   (same site/hypotheses as 2). `cur_piece` from `VSite.kinds`/`two` at `l1_site k` (a non-`V`
   kid exists; `L1Piece` gives the tree edge's `Q`); `cur_vs` from `CloseInv.closed`
   (`CloseAt.att_vs`/`vs_att`) + `VSite.att` (every kid attaches inside `{cur.vStart,
   stackVerts[d]}`); `nxt_ne`/`nxt_touch_top`/`nxt_touch_bot` from `EarFinish.loop1_touch` (range
   entries touch both terminals; nonempty); `proper` from `VSite.pend` at the post-merge state
   (an unprocessed edge at `curV` outside the top entry) with `RangesInv.processed`; `nxt_no_cu`
   from the consumed entry's `EntryR.bond` (`RInvFront.base`) — all `(cur.vStart, curV)` edges were
   bonded into the item loop 1 merged into `cur`, so none is in `nxt` (`RInvFront.disj`).
4. `loop1_rBranch_mid_ctx` (`Proofs/RLoop1.lean`):
   `stackVerts[d + 1] = cur.vStart ∨ Interior (cur.edges) stackVerts[d + 1]` (same site). `L1Close.mid`
   only gives the `l1Edges ∪ t.edges` form; discharge = the `Loop1Spec` field `mid_cur` (E5 above):
   at a reached split closing `t` at `d`, `l1Bot o done = stackVerts[d + 1] ∨ Interior (l1Edges o s
   done) stackVerts[d + 1]`.
5. `rSide_entry_site` (`Proofs/RSide.lean`): at a child-entry site with `v = stackVerts[d]`,
   `IsParent v c`, chain, `RInvTop dfs v d`: (a) `∀ t ∈ tstack, EntryR t → EntryR t` after
   `stackVerts.set! (d + 1) c`; (b) `∀ t ∈ tstack, t.vStart = v → t.topDepth ≤ d`. (b) =
   `EarCtx.split` at the call state (`CtxEntry.depth` for `above`, `vt.topDepth ≤ d`, `below` has
   `vStart ≠ v`) + `EarCtx.base_bot`. (a) = `EntryR.set_stackVerts` for `topDepth ≠ d + 1`, and
   for `topDepth = d + 1` the field E3 (`stab_single`: `EntryR.single` in the updated state — the
   only field reading `stackVerts[topDepth]`; the original `buried_vacuous` candidate is refuted,
   see E3 above).
6. `rSide_vertFree_site` (`Proofs/RSide.lean`): `VertFree v s` at a vertex site with `VertBook v
   false`, `RWalk`. `e < ne`: `EarCtx.vert_edges` (edges under `vertItem v` = `subEdges` of the
   `done` outs with `d ≤ lowval`) + field E2 `bd_free` (no open entry owns one). `e ≥ ne`: E4
   (span items are roots; unique parent) with `EarCtx.v_root`/`vert_free` — not in `RangesInv`/
   `CloseInv` either (`ItemFree`-style facts).
7. `rSide_finish_content_site` (`Proofs/RSide.lean`): at a tree-edge finish site (child walked:
   `RWalk dfs ((v, d, B) :: F) c (d + 1) s`, `FinishGuards`/`FinishBook`/`Frontier`, chain):
   `settled` (every entry loop 1 leaves above the base with `topDepth ≥ d`, `vStart ≠ v` is
   `EntryR` at `feS₁`) — the frontier entries' `EntryR` (item 2's `entry_nxt`) and the merged
   piece (item 1/2's `entry_cur`), R content on R1/R2 + `EarFinish.loop1`; `unwrap` (`hasVert`,
   type 1: `Exempt v d (nxtE (feS₂))`) — provable now from `FinishBook.ear`: `EarFinish.close`
   gives `feS₂`'s stack `c :: py :: vy :: base` with `py.topDepth = lv < d`; `vert_own`
   (`hasVert = false`: `VertFree v (feP v d o s)`) — E2 at the site + transport through loops 1–2
   and the P-check (`L1Keep.edges`, `maybeUnwrapNxt_edges`, merge lemmas) + E4 for `e ≥ ne`.
8. `closeVert_type1_rCloseShape` (`Proofs/RVert.lean`): `RCloseShape g dfs (Pieces.ofItems
   (vPieceItems c py vy)) (vU c py vy) curV stackVerts[lv]` at `feS₂ d o s` of a type-1 tree-edge
   site with `hasVert = true`, `feSingle = false`, `EarClose` stack `c :: py :: vy :: base`,
   `base.length = origTstack`. `wf`/`sub`: `CloseCtx.v_site`'s `VSite` (`kinds`, `two`, `free`) +
   `CloseInv.closed` (`child_two`, `att_vs`); `conn`/`attached`/`touch_s`/`touch_t`/`ne`:
   `VSite.att`/`touch`/`ne` + `EarClose`'s `vert_touch`/`c_edge`; `proper`: `CloseCtx.path lv`
   (an unprocessed edge at depth `lv`) with `RangesInv.processed`; `maximal`: `RangesInv.closed`
   at `feS₂` + `EntryR.maximal` of `c` (item 1) — `py`'s items (`topDepth = lv < d`) are in no
   `RInv*` record, so their maximality is R content on `closed`; `bond`: `EntryR.bond` of `c` +
   `VSite.p_shape`; `single`: `RangesInv.ordered` + `EarClose` (the vertex entry holds the whole
   child subtree); `type1`/`type2`: `EntryR.type1`/`type2` of `c` + laminarity closure (R content).
9. `walk_items_rThreeConnected` (`Proofs/RWalkItems.lean`): `Items.RThreeConnected g (g.walk tern
   (g.dfsForest vo eo)).items` for general `g.WF`. No site; block-local restatement of
   `walk_items_rThreeConnected_of_twoConnected` via `Blocks.lean` (`q_root_covers_block`,
   `nonV_child_cover_of_block`, per-block `Represents`, `RSkel3.rThreeConnected_of_block`); no
   ear/Ranges field involved.
10. `walkTree_rSide_spec` (`Proofs/RInvWalk.lean`): the pre-R-4 statement (yields `BookTree ∧
    RSideTree` from `GuardsTree`/`FrontiersTree` alone), kept only for `walkTreeRReturnSpec`; not
    provable as stated (`BookTree` at a non-root entry is root-threaded ear content) — drop with
    its consumer or restate as `walkTree_rSide`.

Cross-cutting: every unbounded-`e` conclusion (`RTop.disj`, `VertFree`, `vert_own`, `RInvFront.disj`
transport) needs the item-forest facts E4 at the site in addition to the bounded ear fields.

### 4.6 Ranges: `Endpoints`/`Shapes` without the tstack (`Ranges.lean`, `RangesWF.lean`)

Almost all of `Items.Endpoints`/`Items.Shapes` is a consequence of the *final* item tree alone, read
through the DFS edge postorder `σ = edgePostorderForest forest` (each edge listed at its deeper
endpoint, in `walkOut` order — the order `StRef.refTree` uses) plus attachment-point facts; no
ear/tstack reasoning is needed. `Items.Ranges g items σ` (`Spqr/Ranges.lean`) is that final-state
property. Notation: `EdgeBelow i e` = `edgeItem e` is a descendant of `i`; `Att i v` = `v` has an
incident edge below `i` and one not below `i`; `Inner i v` = `v` has an edge and all its incident
edges are below `i`; `IsVs i v` = `v` is one of `vs i`; `PieceEdge i e` = `edgeItem e` is reached
from `i` without passing through a `V` item (the edges of `i`'s own piece, excluding blocks hanging
at its interior vertices).

Fields (all checked by `lake build check_ranges`, 0 violations on seeds 0..400 of `gen.py` and on
a single edge, a self-loop, isolated vertices, a star, a triangle, a digon, two blocks at a cut
vertex, each with `tern ∈ {0,1}`):
* `convex` (R1): for a node item `i` (type ∉ {F, V}) and piece edges `σ[a]`, `σ[c]` of `i` with
  `a ≤ b ≤ c`, `σ[b]` is below `i` — the piece spans an interval of `σ` whose holes are blocks hanging
  at `i`'s interior vertices (which are also below `i`). Laminarity along `IsParent` is
  `ItemTree.Tree.desc_disjoint`.
* attachments (R2): `att_vs` (`Att i v → IsVs i v` for every node item), `vs_att` (`IsVs i v → Att i v`
  at S/P/R and at leaf Q), `vs_ne` (S/P/R endpoints distinct), `interior`
  (`IsParent i (vertItem v) ↔ Inner i v ∧ no child owns all of v's edges`, at S/P/R),
  `child_two` (a non-V child of an S/P/R node has two distinct endpoints), `io_parent`
  (an I/O item's parent is a Q), `q_leaf` (a leaf Q's `vs` is its edge), `q_root` (a Q with children:
  `vs = (u, none)`, `u` on the edge; loop ⇒ `ch = [c]` with `vs c = (u, none)`; non-loop ⇒
  `ch = [c, vertItem w]`, `w < nv`, `{u, w}` = the edge, `vs c = (a, b)` a permutation of `(u, w)`;
  `c ∉ {F, V}`, and `c` a leaf if it is a Q), `q_under_v` (a child of a V item is a block root,
  `ch c ≠ []`; consumed by the R coverage bridge), `p_shape`, `s_order`, `r_shape` as in `Shapes`.

False candidates (recorded, restated):
* *All* edges below a node form an interval of `σ` — **false**: seed 348, S node `{e2, e0}` at
  σ-positions `{4, 5}` has a digon hanging at its interior vertex 1 at positions `{0, 1}`, with another
  piece's back edges at `{2, 3}` between. Restated on piece edges with holes allowed (`convex` above);
  the strict form "piece edges are consecutive" is also false (seed 0, Q item 9: the hanging block
  sits inside the piece's interval).
* `Shapes.s_shape` with `2 ≤ xs.length` — **false**: a triangle gives S with `vs = (0,1)`, one V child
  `[2]`, `virt = [(0,2), (2,1)]` (the third edge is the parent cap). Corrected to `1 ≤ xs.length`
  (`ItemSpec.lean`, `RelabelRep.lean`, `RelabelSt.lean`).
* `Shapes.r_shape` with `6 ≤ virtualEdges.length` — **false**: K4 gives R with 5 virtual edges
  (+ parent cap). Corrected to `5 ≤`.
* `Shapes.q_children` with `type c ∉ [F, V, Q]` — **false**: a digon (tree edge + a parallel back
  edge; seed 380, `e = 1`) gives the block root `Q → [Q leaf, V w]`. Corrected to
  `type c ∉ [F, V] ∧ (type c = Q → ch c = [])` (also in `WalkTyping.walk_q_children`,
  `RelabelOwn`, `RelabelRep`, `RelabelAdj`).
* `Shapes.canonical` — **false for `tern = true`**: P under P in ~55% of seeds (e.g. seed 0, items
  62 → 61), S under S (seed 386). Canonicity is therefore no longer a clause of `Items.Shapes` /
  `Represents`: it is the separate predicates `Items.Canonical` / `SpqrTree.Canonical`, with
  `spqrTree_canonical : (g.spqrTree false vo eo).Canonical` (via `RelabelOK.canonical` /
  `relabelTree_canonical`) and the walk-side `walk_canonical` (`WalkItemsWF.lean`) for
  `tern = false` only, a proved assembly over the named hypotheses `walk_canonInv`/`walk_ternarize`
  (below).

`RangesWF.lean`: `wf_of_ranges : Items.Tree → Items.TypingFacts → Items.Ranges → Items.WF`
(`endpoints_of_ranges`, `shapes_of_ranges`; standard axioms). `TypingFacts` = the `i_o_leaf`,
`vs_shape`, `vs_lt` fields of `WalkTyping`. Every clause of `Endpoints` and `Shapes` is derived
from `Ranges` (the Q clauses from `q_leaf`/`q_root`; `child_vs_in_parent`
from `vs_att` at the child + `att_vs`/`interior` at the parent + `desc_disjoint`; `nv_nodup` from
`vs_ne`, `vs_att` vs `interior` and the injectivity of `vertItem`; `separation` = `att_vs`;
`interior` = `interior`). `convex` is not needed for `WF` at all — it is the field the walk-time
invariant (`RangesInv`, step 4) and the R/st layers consume. `walk_items_wf` (`WalkItemsWF.lean`, under
`g.WF`/`OrderOK`) is `wf_of_ranges walk_tree.toTree walk_typing.toTypingFacts walk_ranges`: `Tree`
comes from `EarWalk.walk_tree` (which, since `walk_sides` moved under the ear layer, also takes
`∀ t ∈ forest, t.Ends g`, supplied by `dfsForest_ends` (`WalkInv.lean`, proved: the `adjacency`-level
endpoint fact `adjacency_ends` — `(y, e) ∈ adj[x]` means `g.edges[e]` is `(x, y)` or `(y, x)`, from
`AdjInv` — carried through the `dfsStep`/`forestStep` folds by `dfsVisit_ends`), the `Endpoints`/`Shapes` half
rests on the named admission `walk_ranges` instead of the ear layer. `WalkItemsWF` sits above `WalkCover`/`WalkTyping`
(the DFS prerequisites `Bounded`/`ForestOK`/coverage come from `dfsForest_spanning'`/`dfsForest_wf`;
the empty graph, where `walk_tree`'s `0 < g.nv` fails, is the initial items), so the st layer
(`StWalk.walk_st`/`spqrTree_st`) takes `Items.WF` as a hypothesis and `StOriented`/`Correctness`/
`RelabelChildShape` supply it; `spqrTree_represents`/`spqrTree_canonical`/`spqrTree_childShape`
therefore carry `g.WF`/`OrderOK` like `spqrTree_wf'`.

`RangesInv.lean` (step 4, the walk-time invariant behind `walk_ranges`): `WalkState.RangesInv σ n D s`
after the first `n` edges of `σ` at depth `D`: `inv : Inv' D` (the attachment half, `WalkSpec.lean`;
there is no "≤ 2 terminals" claim — an open chain is attached at every path vertex between its
`topDepth` and `D`, and at `vStart`s of entries above it, exactly `Term'`), `processed` (entries own
processed edges only), `ordered` (each entry's *piece* edges lie strictly after all edges of lower
entries in `σ`; `TEntry.piece` = below a non-V span item without passing a V
item), `convex` (an entry's piece is an interval of `σ` whose holes are *edges of the entry* — blocks
hanging at its interior vertices), `closed` (`Ranges.convex` for every allocated node item). It is
schedule-agnostic (no ear, no merge order). Preservation, all standard axioms:
`RangesInv.alloc`/`pushVert`/`pushEdge` (the new edge is `σ[n]`, so `n+1`), `mergeTop` (under
`MergeOk` + edge-disjointness as `Inv'.mergeTop`, plus the *local* range condition `hadj`: the two
merged entries are adjacent in `σ` — every edge strictly between a piece edge of `nxt` and one of
`cur` is an edge of one of them), `finishTop` (as `Inv'.finishTop`: the closed item takes over the
entry's edges and a subset of its pieces, so `closed` for it is the entry's `convex`). The
`finishEdge`/`walkTree` preservation is proved under its range-side site hypotheses;
their instantiation for the real walk and `walk_ranges` from the final state stay admitted.
Checked at every `finishEdge` by `checks/RangesInvCheck.lean` (`lake env lean`; reimplements
`walkTree` around `finishEdge`, asserts the instrumented items equal the library walk's; seeds
0..400 × tern, tiny graphs): `processed`, `ordered`, `convex`, `closed` 0 violations.
`RangesFrontier.lean`: `boundaryAdj_of_book` discharges `BoundaryAdj` from `FinishBook.ear`,
the pending edge's postorder position and `o.block <:+: σ` (`subEdges_interval`,
`subEdges_iff_mem_block`, `infix_interval`). `RangesInv.mergeAdj_of_frontier` derives
`MergeAdj` when the top two entries lie above `origTstack`, `FrontierOwns` covers an interval,
and the state has `RangesInv`: the interval supplies a stack owner of each gap, and the
stronger `ordered` excludes lower owners (`mergeAdj_of_cover`). `iter_merge_ranges` /
`iter_mergeAdj` mutually thread that argument with range preservation through guarded
merge iterates; `mergeLateAdj_of_frontier` and `loop3Adj_of_frontier` discharge loops 2 and 3
from the corresponding `Frontier` fields and `MergeLateOk` / `CloseVertOk` guards.
`FrontierOwns.replaceNxt` / `mergeTop` / `congr` / `unwrap` preserve ownership through the
local stack changes. `loop1BodyAdj_of_frontier` discharges both sites of a loop-1 body,
including its optional S merge and following unwrap, from the frontier's 3/2-entry bound.
`RgStep.iter_loop1_ranges` threads preservation through these bodies, and
`closeEarsAdj_of_frontier` supplies all of `CloseEarsAdj` from the loop-1 frontier export.
`finishTailAdj_of_vert` discharges the tail merge because the pushed vertex has no piece
edges; the `PushVertR` processed-edge bound remains a separate obligation.
`Frontier.owns_late` obtains ownership at the exit of loop 2.
`closeVertAdj_of_frontier` supplies all three vertex-close adjacency fields; its three-entry
frontier bound follows from `FinishBook.late_length`.
`finishPAdj_of_frontier` handles the P unwrap/merge when an enclosing interval frontier owns
both entries; supplying that enclosing frontier for the actual P merge remains open.
`MergeBaseCover` instead splits ownership at the child's starting postorder position:
the settled base covers the preceding suffix and the frontier covers the new child block.
`RangesInv.mergeAdj_of_baseCover` and `finishPAdj_of_baseCover` derive adjacency from this
explicit induction hypothesis (`FinishPCover`), with no new admission. Both ownership clauses
at actual P sites pass seeds 0–400 × both modes in `RangesInvCheck`.
`RangesSchedule.finishR_of_cover` assembles the entire `FinishR` bundle from the ear guards,
book and frontier, the postorder position, and `FinishCover`. The latter retains only the
P-site ownership and the unpushed vertex's processed-edge bound; all other adjacency is derived.
`Step.pushVertR` transports that bound through the returning-edge primitives.
`RangesWalk.scheduleTree` / `scheduleOuts` / `scheduleOut` prove the mutual instantiation of
`RgTree` / `RgOuts` / `RgOut` from `Guards*`, `Book*`, `Frontiers*`, and the reduced `Cover*`
ownership obligations. `PostAt` derives every edge's position and subtree interval from the
DFS edge-postorder concatenations. These theorems use only standard axioms; instantiating
the remaining `Cover*` obligations on actual DFS walks is still required.
`init_rangesInv` and `RangesInv.root_append` handle the initial state and attaching each root.
`forest_ranges_of_cover` and `walk_rangesInv_of_cover` assemble the forest using `RootState`,
`walkTree_frontiers` and the postorder permutation. Their only nonstandard axiom dependency is
the ear layer's existing `walkTree_book` admission (through `RootState.book` / `RootState.step`);
the remaining range hypothesis is explicitly `RootsCover`, not an admitted local proof.
The actual `walk_rangesInv` assembly now lives in `WalkItemsWF` (avoiding the ear import cycle)
and threads `g.WF` and both `OrderOK` hypotheses through `walk_ranges`. Its residual range
admission is named `walk_rootsCover`; the statement specifies precisely the P-site coverage
and unpushed-vertex bounds still to establish. `Place.pushVertR` proves the latter bound from
the existing placement invariant once its pushed-edge predicate is bounded by the prefix.
`PostAt.idx_bounds` and `pushed_past` prove the prefix arithmetic; `walkTree_past`,
`walkOutPre_place` / `walkOutPre_past`, and `finishEdge_past` transport it through the child,
pre-push, and edge finish. `RootState.pushVertR` supplies it at each root from the already
processed forest prefix. `walkOut_past` carries the bound through a complete out-edge call.
`FinishPOwnership` isolates the two P-site fields from the vertex bound, and
`coverOut_back` / `coverOut_tree` construct `CoverOut` using placement and postorder bounds,
with only that P ownership and the recursive child's `CoverTree` left as hypotheses.
All these lemmas use only standard axioms. The enclosing `CoverOuts` / `CoverTree` mutual
induction and forest assembly of `RootsCover` remain unfinished. The first-tree-edge P check with
`hasVert = false` is false by `ear_condP_tree`; `FinishCover` therefore requires P ownership
only for the `p_vert` and `p_back` sites, not `p_tree`.
`FinishCover` asks for P coverage only on the return branch (`lowval < d`), where P finishing
actually executes. The empirical checker now also checks the unpushed vertex bound and every
`CloseFacts` clause on live items (`cnt > 0`) before and after each edge finish. Seeds 0–400,
both modes, and the tiny examples have no violations. Allocated but loose items are excluded:
freshly allocated nodes and reopened nodes are not closed records until they are reattached.
`RangesClose.CloseAt` records all the checked attachment, child, Q/I/O, and P/S/R shape
clauses for one item. `CloseInv` requires this record for the root, every vertex item (even
before it is pushed), and items with positive span/child occurrence count. The stronger
all-vertices clause has zero required failures on seeds 0–400 × both modes; it is checked
by `checkClose`. Its initial-state proof is complete. `CloseInv.of_tree` recovers
`Items.CloseFacts` when every non-root item has a parent; this uses the existing `Items.Tree`
contract, not another reachability admission. `CloseAt.frame` transports a record from its
own type/terminals/children, its children's fields, and their edge-below predicates.
`CloseInv.alloc`, `pop`, `mergeTop`, `unwrap`, `modifyLoose`, and conditional `push`/`finishTop`
are proved; closing a
zero-count item preserves every other record via `CloseAt.modify_of_not_below`.
`CloseAt.root_of_children`, `vertex`, `leafQ`, and `leafIO` construct complete per-type records.
`CloseInv.pushVert` now preserves the invariant without a new-record hypothesis. `pushEdge`
constructs a leaf-Q record from its distinct endpoints, graph-edge equality, and both terminal
attachments; those actual-site attachment hypotheses still need the ear/DFS context.
`CloseInv.writeChildren` preserves all other records when writing a no-parent item, and
`CloseInv.root_append` constructs the root record and preserves all records when appending
vertex items. These constructors and the strengthened preservation lemmas have only standard
axioms. `modifyLoose` now requires the modified item to lie beyond the fixed vertex interval.
`walkOutPre_closeInv` handles the real pre-push code, and `rootAppend_closeInv` handles the
actual pop/append action from `RootOK` and `Place`. `CloseInv.vertex_append` constructs the
vertex's updated record at a boundary from the already-closed, nonempty Q child, the old
V-to-Q child typing, and the vertex's zero count. It remains to construct that Q record at
the loop/bridge/block boundary, rather than assuming it.
Node-finish and boundary new records and edge terminal attachments are proved per site
(`RangesCloseSites.lean`) and threaded through the walk (`RangesCloseTree.lean`, DFS-side site
facts proved by `dsTree`); the five per-site theorems
`closeCtx_bd_vert`/`bd_node`/`p_site`/`v_site`/`l1_site` (`RangesCloseTree.lean`) are proved from
`CloseBase` and the block-content record `CloseContent` (`RangesCloseContent.lean`), which the
backbone induction supplies at every `finishEdge` (named hypothesis `closeBase_content`).
Canonicity (`Items.Canonical`, `tern = false`: no S under S / P under P) is the walk-state
predicate `WalkState.CanonInv s := s.ternarize = false → Items.Canonical s.items` (`WalkWF.lean`),
kept by every `finishEdge` primitive: the only item writes are `allocItem` (a new leaf),
`finishTstackTop x` (writes the active side of the top entry as `x`'s children) and the Q/V `ch`
writes of `finishBoundary`, so `finishEdge_canon : CloseBase → CloseCanon → CanonInv s →
CanonInv (finishEdge … s)` (`RangesCanon.lean`, standard axioms) needs exactly the site record
`CloseCanon`: at the three `finishTstackTop` sites (P merge `p_site`, type-1 vertex close `v_site`,
loop-1 iterate `l1_site`) no child on the finished side has the type of the node being finished
(`FinishCanon x s`; checker `checkCanon`, fields `canon_p`/`canon_v`/`canon_l1`, 0 violations
0..3000 × both modes). Threading `CanonInv` to the final state is the backbone's: `walk_canonInv`
(`CanonInv` of `g.walk tern (g.dfsForest vo eo)`) and the frame fact `walk_ternarize` are the
named hypotheses, and `walk_canonical` is their assembly (`WalkItemsWF.lean`). `walk_q_children`
is `walk_items_wf.shapes.q_children`, restated over `g.dfsForest vo eo` under `g.WF`/`OrderOK`
(its former arbitrary-forest form had no callers).

The five final-items facts of `spqrTree_pieceSep` are handled the same way (`RangesPiece.lean`):
`PieceFacts s` (`Items.QUpper`/`RootSep`/`RootV`/`QChildVs`/`PChildVs` of `s.items`) and
`PieceInv d s` (`PieceFacts` plus `root_path`: no vertex item of the open path `stackVerts[0..d]`
lies below a root child — the clause that keeps `RootSep` through `finishBoundary`'s append under
`vertItem curV`). Preservation is proved for every primitive (`PieceInv.frame`/`exit`/`setVert`/
`init`/`alloc`/`modify`/`modifyVs`/`vAppend`/`mergeTop`/`mergeLoop`/`mergeLate`/`vertPre`/`retarget`/
`pushVert`/`finishTail`/`l1Type`/`maybeUnwrap`/`finishTop`/`vertUnwrap`/`pSite`/`qClose`/`boundary`/
`l1Body`/`rest`) and `finishEdge_piece : CloseBase → CloseContent → ClosePiece → PieceInv d s →
PieceInv d (finishEdge … s)` (standard axioms), where the site record `ClosePiece` holds the
boundary orientations (`bd_bridge`: the bridge's `makeVs` is `(curV, dest)`; `bd_node`: the closed
node's `vs` is `(curV, dest)`) and `FinishPiece x s` at the three `finishTstackTop` sites (an
`x` of type P only receives children already carrying its `vs`). Between trees only `PieceFacts`
holds (`PieceFacts.rootAppend`; `root_path` is a within-tree clause — at a root end the finished
root itself is a root child), and the next root re-enters `PieceInv 0` (`PieceFacts.root`).
Checker (`checks/WalkInvCheck/Ranges.lean`): `checkPieceInv` evaluates every clause at the tree
entry/end, before/after every `finishEdge` and (`PieceFacts` only) after every root append;
`checkPiece` evaluates `ClosePiece` and the frame hypotheses of `finishEdge_piece`; fields
`piece_*`, 0 violations 0..3000 × both modes. Threading to the final state is the backbone's:
`walk_pieceInv : (g.walk tern (g.dfsForest vo eo)).PieceFacts` (`WalkPieceSep.lean`) is the named
hypothesis and `walk_q_upper`/`walk_root_sep`/`walk_root_v`/`walk_q_child_vs`/`walk_p_child_vs`
are its projections.
`walk_closeFacts` is the proved tree/typing assembly from `walk_closeInv`, not an independent
close-facts admission.

The former unrestricted `walk_closeFacts` statement is false for malformed graphs:
`nv = 0`, `edges = #[(0, 0)]`, `vo = eo = []`, `tern = false` leaves the Q item's `vs` as
`(none, none)`, contradicting `q_leaf`. `closeFacts_needs_graph_wf` in `checks/RangesInvCheck.lean`
is a kernel-checked counterexample (no `native_decide`). The wrapper now takes `g.WF`, `OrderOK`
for both orders, and `0 < g.nv`; `walk_items_wf` already handles the empty valid graph separately.
No protected proof-target definition was changed.
All twenty-two lemmas have only standard axioms. Full `FinishAdj` still needs the P merge
into the base; the child's frontier lemma alone applies only
to merges above the split. The tail merge against a newly pushed V entry has no piece
on that entry, so its `MergeAdj` obligation is vacuous.
The final assembly must avoid importing `RangesFrontier` into `WalkWF`: the current
`EarFrontier → RInvFrame → EarShape → EarSpec → WalkWF` path would make a cycle.
`WalkItemsWF`, already above the ear layer, can host that assembly.
False candidates (recorded): *every processed edge is owned by some entry or the root* — false
(star: after the bridge `0–1` closes its Q hangs under `vertItem 0`, which is on no entry yet);
*≤ 2 attachments per entry* — false (K4, entry `vStart 2, topDepth 0` attached at `0, 1, 2`);
*attachments ⊆ `vStart` ∪ path* — false (seed 1: attached at the `vStart` 4 of the entry above),
i.e. `Term'` is tight. **Saturation** (R layer): "below the frontier, no two adjacent entries (no run
of consecutive entries) have attachment union ≤ 2" is false for the attachment-set measure, pairwise
and run-wise: a vertex entry with no hanging edges has an empty attachment set (triangle + triangle
at vertex 2, back edge `2→0`: entries `vertItem 2`, `vertItem 1`, union attachments `{2}`), and a
piece whose `vStart` is a path vertex `w` sits next to the (lazily merged) `vertItem w` entry with
union attachments = its own two (K4 at `curV 3`: `[(2,0,{e(0,2)}), (2,2,[v2]), (1,1,[v1])]`, union
`{0, 2}`). Saturation must be stated in terms of the entries' terminals (`vStart`, `stackVerts[topDepth]`)
and the depth separation of vertex entries, which is the merge test the code applies; it is left to
the R layer as a hypothesis over `tstack.drop (tstack.length - origTstack)` (`Frontier`,
`Proofs/RInvFrame.lean`) and is not a field of `RangesInv`.

The item shape counts exclude the parent cap. `Shapes.r_shape` therefore requires five
non-V children, giving six skeleton edges with the cap. `Shapes.s_shape` requires one
V child, giving a path of two child edges and a three-edge cycle with the cap.
`checks/SkeletonShapeCheck.lean` kernel-checks both sharp cases on actual walk outputs:
K4 has R item 11 with five non-V children and two V children; the triangle has S item 7
with two non-V children and one V child. The former bounds of six non-V R children and
two V S children were false. `RelabelOK.nvList_S` and `.nvList_S'` use the corrected S bound;
the public minimum skeleton sizes remain six edges for R and three vertices/edges for S.

The public R-connectivity statement also needs valid input orders. On the valid K4 graph
`[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]`, `spqrTree false [] [6]` has an R node at index 3
with only one vertex. Edge ID 6 is outside the graph's edge range. The corresponding walk
item 11 still has type R but has terminals `(some 0, none)`, so it cannot satisfy `RSkel3`.
`checks/RInvalidOrderCheck.lean` kernel-checks the graph WF, invalid order, output node,
item terminals, and failures of both connectivity specifications using standard axioms.
`items_r_three_connected`, `spqrTree_r_three_connected`, and `spqrTree_represents` now take
the standard `g.WF`, `OrderOK g.nv vo`, and `OrderOK g.ne eo` hypotheses used by the DFS
and `spqrTree_wf'` interfaces. The planar embedding soundness caller already supplies these.
The item and public R-connectivity theorems remain admitted; correcting their domain does
not discharge the walk invariants or the block-local reduction.

*Adjacency of the merged entries (`RangesInv.mergeTop`'s `hadj`) is a schedule fact, not a range
fact.* `checks/RangesInvCheck.lean` candidates at every `finishEdge`: `cand_gapfree` (adjacent
entries are `σ`-adjacent: the positions strictly between `nxt`'s piece and `cur`'s piece are edges
of one of them) has no violations on seeds 0..400 × tern + tiny graphs, but it does not survive
through entries with an empty piece (vertex entries), and the two generalisations that would make
it preservable fail: `cand_runs` (runs of consecutive entries are `σ`-convex, 1570 observations)
and `cand_reach_all` (every position from an entry's first piece edge up to the edge being pushed
is owned by an entry at or above it, 448). First counterexample, seed 1, `tern = false`, at
`curV = 6`, `d = 6`: the self-loop `e13 = (5,5)` is a hanging block under `vertItem 5`, whose
vertex entry is not pushed yet (vertex 6 is a type-2 child of 5), and `σ`-position 4 of `e13` lies
between the pieces of `(1,1,1,[Q8])` and `(6,1,2,[Q15])`/`(6,6,2,[vert 6])`. Holes between
non-adjacent pieces are therefore edges below the *unpushed* vertex items of the DFS path, and the
fact that the code never merges across such a hole (the vertex entry is pushed and merged first)
is `hasVert` bookkeeping, i.e. schedule-specific. Consequently `finishEdge_rangesInv` takes `hadj`
per merge site as a hypothesis bundle (`FinishAdj`, `RangesStep.lean`, mirroring `FinishOk`'s
`MergeTopOk` sites) instead of deriving it from `ordered` + stack adjacency.

*Induction (`RangesStep.lean`, `RangesTree.lean`, `RangesFinal.lean`).* `finishEdge_rangesInv`
(standard axioms): `RangesInv σ n D` is preserved across the returning-edge branches of `finishEdge`
under `FinishOk` + `FinishAdj` (`RgStep.loop`/`loop1Body`/`closeEars`/`mergeLate`/`closeVert'`/
`finishP`/`finishTail`/`finishRest`, each `Step`-shaped lemma paired with the range invariant).
`walkTree_rangesInv` (`rgTree`/`rgOuts`/`rgOut`, mutual like `invTree`/`GuardsTree`): under
`GuardsTree`/`BookTree` and the range-side hypotheses `RgTree` (`σ[n]? = some o.e` at every edge,
`FinishR` = `FinishAdj` at returning edges / `BoundaryAdj` at block boundaries) the state after
`walkTree t d` satisfies `RangesInv σ (n + t.edgePostorder.length) d`. The block-boundary case
`finishBoundary_rangesInv` is proved under `BoundaryOk` and the unchanged `BoundaryAdj`:
`RangesBoundary.lean` transports ranges through leaf endpoint writes, stack pops, and free-root
child writes; `Items.CapRange`/`cap_convex` fill the popped spans up to the pending edge.
Both it and `walkTree_rangesInv` print only `[propext, Classical.choice, Quot.sound]`.
The stronger `ordered` clause (pieces versus all lower edges, not just lower pieces) is preserved
by all the same local lemmas and passes seeds 0..400 inclusive × both modes. The checker also
checks `MergeAdj` at the actual loop-1/2/3, vertex, P, and tail merge sites, with zero violations.
`ranges_of_rangesInv` (standard axioms):
`RangesInv σ n D s → WalkTyping s.g s.items → Items.CloseFacts s.g s.items → Items.Ranges s.g s.items σ`.
Exactly two `Ranges` clauses come from the invariant: `convex` (= `closed`) and `att_vs` for the
*allocated node items* (`Inv'.nodes`: `TwoAttached` at their `vs`; `I`/`O` leaves own no edge).
`Items.CloseFacts` is the remainder, which `RangesInv` does not carry because it is about an item
at its close (its `vs` just written from the entries' terminals, its children's shape): `att_vs`
for the Q items (not in `Inv'.nodes`), `vs_att`, `vs_ne`, `interior`, `child_two`, `io_parent`,
`q_leaf`, `q_root`, `q_under_v`, `p_shape`, `s_order`, `r_shape`. The remaining range/close
hypothesis is `closeBase_content` (`CloseContent` — the boundary-entry and `PSite`/`VSite` block
contents at a `finishEdge`, stated over the exported `CloseBase`; discharged by the backbone
induction from the primitives).
`walk_rangesInv` assembles the forest from the ear exports plus `RootsCover`; `walk_closeFacts`
assembles the final close records using the item tree and typing. The mutual scheduler and
local close-record proofs use only standard axioms. The forest assembly also consumes the
ear layer's `walkTree_ear` admission; the final wrappers retain `sorryAx` until all of these
obligations are discharged. `walk_ranges` (`WalkItemsWF.lean`, under `g.WF`, both `OrderOK`,
`0 < g.nv`, and a bounded/covering forest) is `ranges_of_rangesInv walk_rangesInv walk_typing
walk_closeFacts` with `walk_g : (g.walk tern forest).g = g` from `walkForest_typing`.

Remaining close-site obligations, precisely:
* `closeEars` / `finishBack`: apply `CloseInv.pushEdge` after `finishSetup` by deriving distinct
  endpoints and both terminal `Att` witnesses for the leaf Q from the returning-edge contract.
  The Q constructor itself and preservation of the other records are proved.
* `loop1Body`: proved (`loop1Body_closeAt`), see below.
* `closeVertTail`: proved (`closeVertTail_closeAt`), see below; the `none` arm only
  merges/modifies the stack and creates no item record.
* `finishP`: the analogous P record after unwrap/merge, including the minimum virtual-edge
  count and absence of V children; both newly allocated and reused item paths remain.
* `finishBoundary`: proved (`finishBoundary_closeAt`), see below.
* Actual vertex pre-push and root pop/append are proved (`walkOutPre_closeInv`,
  `rootAppend_closeInv`). The end-of-tree and first-edge vertex pushes use `CloseInv.pushVert`.
  These are threaded through the mutual walk/forest induction, with placement and `RootOK`,
  in `RangesCloseTree.lean` (`ccTree`/`ccOuts`/`ccOut`, `forest_closeInv`).

These sites are now separate named admissions in `RangesCloseSites.lean`, each stated as
`CloseInv` preservation of its block at the state where the block runs, under `CloseCtx` (the
`finishEdge` call's `RangesInv`, `Shape`, `FinishGuards`, `FinishBook`, `Frontier`, `FinishR`,
`CloseInv`, plus the DFS facts `ends : o.Ends`, `path` — the open path's tree edges
`(stackVerts[k], stackVerts[k+1])`, `k < d`, each still pending, i.e. after `o.e` in `σ` — and
`dest_edge` — a returning child has an edge `≠ o.e` at `o.dest`): `closeEars_closeAt` (`ceS₁`
after `feS₀`) and `finishBack_closeAt` (`feBack`) are proved (`leafQ_closeInv`: the leaf-Q
record from the written terminals `setSides dir (stackVerts[d]) o.dest`, distinct by
`EarFinish.path`/`path_child`, with the parent edge / `dest_edge` / the path edge at
`stackVerts[lv]` as the second `Att` witness; `feS₀` itself is `CloseInv.modifyLoose` since
the pending Q is in no span and has no parent). `finishBoundary_closeAt` (the whole boundary
`finishEdge`) is proved, mirroring `finishBoundary_rangesInv` step by step: the pending Q's
terminal `(some curV, none)` is a loose write (`CloseInv.modifyLoose`); the I/O leaf is
allocated (`CloseInv.alloc'`, the `Place`-free variant from `Shape`) and its terminals written
loose; the one or two popped entries transfer through `CloseInv.pop` with
`BoundaryOk.gone`/`gone₂`; the root-Q record is `CloseAt.rootQ` (children `[I, vertItem dest]`
at a bridge, `[node, vertItem dest]` at a completed block, `[O]` at a self-loop; the Q's `Att`
is `curV` only, from the popped entries' terminals via `EarFinish.sub_cover`/`dest_edges`/
`sv_child`/`path_child`, `Inv'.entries`, and the open-path edges of `CloseCtx.path`) written by
`CloseInv.writeChildren`, after which `CloseInv.vertex_append'` (the `Shape` variant) hangs the
Q under `vertItem curV` (`cnt (vertItem curV) = 0` from `BoundaryOk.v_root`/`v_free`). It
consumes four new `CloseCtx` fields, facts about the pre-state that the walk induction must
supply and that no current `EarFinish`/`RangesInv` field exports: `dest_lt : o.dest < nv`;
`bd_loop` (a non-tree boundary edge is a self-loop, `o.dest = curV`); `bd_vert` (the block
entry `(o.dest, d+1)` — the top entry at a bridge, the second at a completed block — has
`spans.2 = [vertItem o.dest]`, i.e. nothing was merged into the child vertex's entry);
`bd_node` (at a completed block the top entry's `spans.1` is a single non-V/F node, a leaf if
Q, whose terminals are `{curV, o.dest}`). The checker evaluates these four at the pre-state
of every block boundary (`checkCtx`, kinds `ctx_*`). `finishP_closeAt` (from `feRest`, the
state `finishRest` runs from) is proved: when `condP` fails, `finishP` is the identity
(`finishP_closeInv_of_not`); otherwise `PSite.closeInv` runs the three primitives by their run
equations (`maybeUnwrapNxt_run_eq`/`mergeTstackTops_run_eq`/`finishTstackTop_run_eq`): the P
item is either allocated (`CloseInv.alloc'`) or the second entry's own P reused (its children
become that entry's active side; `CloseInv.unwrap`), the two top entries merge
(`CloseInv.mergeTop'`), and the merged active side becomes the children of the P record
(`CloseInv.finishTop'`), whose clauses are `CloseAt.pNode'`: distinct terminals
`{curV, stackVerts[lv]}`, every child S/P/R/leaf-Q with those terminals, `Att` of the P equal
to the two terminals (a child's attachment is a terminal or interior; each terminal has an edge
below and one outside), at least two children, no V child. The facts about the two entries come
from one new `CloseCtx` field `p_site : lowval < d → isType1 → condP → PSite curV lowval feRest`,
where `PSite` (`RangesCloseSites.lean`) records: stack head `(curV, lv)`;
`curV ≠ stackVerts[lv]`; each of the two top entries has a single S/P/R/Q item on its active
side with terminals `{stackVerts[lv], curV}` and an empty passive side; its attachments are
terminals or interior; it touches both terminals; each terminal has an edge outside both
entries; S/P/R/leaf-Q kinds (also for a P item's children); every spanned item is a root
listed once on the stack. The checker evaluates `PSite` at every `feRest` state where `condP`
holds (`checkP`, kinds `psite_*`): 0 failures. `closeVertTail_closeAt` (`cvS₅ → feS₃`) is
proved: for a type-2 edge `feS₃` is `cvS₅` itself; for a type-1 edge it is
`finishTstackTop x` at `cvS₅` for the item `x` returned by `maybeUnwrapNxt`, and
`VSite.closeInv` runs it by `finishTstackTop_run_eq` and `CloseInv.finishTop'`, with the
S/R record `CloseAt.node` (distinct terminals; `Att` of the node equal to the terminals;
each terminal touched and with an edge outside; children S/P/R/Q/V, non-V ones with two
distinct terminals; a vertex item is a child iff all its edges lie below the node and in no
single child; the S path and R shape clauses on the node's virtual edges). Its facts come from
one new `CloseCtx` field `v_site : isTree → lowval < d → hasVert → isType1 → ∃ t, VSite curV x t
cvS₅`, where `VSite` (`RangesCloseSites.lean`) records at `cvS₅`: top entry `t` with
`vStart = curV` and an empty passive side; `x` free (`ItemFree`) of type S or R;
`curV ≠ stackVerts[t.topDepth]`; the active side's kinds and two-terminal children; its
attachments are terminals or interior, both terminals touched and with an edge outside the
side; the interior-vertex equivalence; and the S/R clauses on the side's virtual edges
(`vKids`/`vTerms`/`vVirt`). The checker evaluates `VSite` at every type-1 vertex close
(`checkV`, kinds `vsite_*`): 0 failures. `loop1Body_closeAt` (one loop-1 iteration under the
loop condition) is proved by the same record: `loop1Type` only re-sets a stack direction and
merges (`CloseInv.l1S₁`, via `CloseInv.frame`/`mergeTop`), `maybeUnwrapNxt` allocates or reuses
(`CloseInv.maybeUnwrap`: `CloseInv.alloc'`/`unwrap` from `Shape` and a two-entry stack), the
merge is `CloseInv.mergeTop` (`after_mergeTstackTops`), and the close is `VSite.closeInv` with
`curV := t.vStart` of the merged entry; `VSite` and `CloseAt.node` carry the P clause too
(`p_shape`: at least two virtual edges, all equal to the terminals, no V child). The facts come
from the new `CloseCtx` field `l1_site : isTree → lowval < d → ∀ k, (loop condition up to k) →
Shape (l1S₁ (l1Iter k)) ∧ (two entries) ∧ ∃ t, VSite t.vStart x t (after mergeTstackTops (l1S₂
(l1Iter k)))`. The checker evaluates it at every loop-1 iteration (`checkL1`, kinds `l1site_*`):
0 failures. All six close sites are now proved. The checker evaluates
every `CloseAt` clause at every one of these block boundaries (`closeSites`, with the theorem
name as the site; the frame blocks `feS₀`/`mergeLate`/`closeVert'.pre`/`finishTail` as well):
0 failures on seeds 0..400 × both modes + tiny graphs. `walk_closeInv` is now the proved
assembly over these six names (`RangesCloseTree.lean`): `finishEdge_closeInv` runs the six
sites in `finishEdge`'s block order from a `CloseCtx` (`CloseInv.loop`/`mergeLate`/`vertPre`/
`vertUnwrap`/`retarget`/`finishTail` for the frame blocks, `Step.*` for the intermediate
`Shape`s); `ccTree`/`ccOuts`/`ccOut` mirror `rgTree` with `CloseInv` threaded (`walkOutPre_closeInv'`,
`CloseInv.frame`, `CloseInv.pushVert`), building each `CloseCtx` (`CloseCtx.of_exports`) from
`CloseBase`: the range/ear/frontier/R-side exports (`RgS`, `FinishGuards`, `FinishBook`, `Frontier`,
`FinishR`), the prior `CloseInv`, and `DfsSite` — the DFS-side site facts `pos`/`block` (postorder
position and block), `ends`/`path`/`dest_edge`/`dest_lt`/`bd_loop` — carried in the `wp` shape of
`RgTree` (`CsTree`/`CsOuts`/`CsOut`) and proved on the real walk by the mutual `dsTree`/`dsOuts`/`dsOut`
from `DfsTree.WF`/`Ends`, the ancestor path (`AncPath`, joined by edges scheduled after the current
block, from `PostAt`), the `Keep` frame and `classify` (`classify_back`). `forest_closeInv`/`walk_closeInv'`
mirror `forest_ranges_of_cover` with `rootAppend_closeInv`. The per-site theorems, stated over
`CloseBase` at a `finishEdge` and checked 0..400 × both modes (0 failures) under the same checker
labels: `closeCtx_bd_vert`/`closeCtx_bd_node` (boundary entry shapes, `checkCtx`), `closeCtx_p_site`
(`PSite`, `checkP`), `closeCtx_v_site` (`VSite` at the type-1 vertex close, `checkV`), `closeCtx_l1_site`
(`VSite` at each loop-1 iteration, `checkL1`), are proved from `CloseBase` and `CloseContent`
(`RangesCloseContent.lean`): the item-level facts no export gives — `bd_vert` (block entry
`spans.2 = [vertItem dest]`), `bd_node` (completed block: one non-F/V node, leaf if Q, terminals
`{curV, dest}`), `p_site : PContent` (= `PSite` minus `shape`/`ne`/`pend`) at `feRest`, `v_site :
VContent` (= `VSite` minus `shape`/`ne`/`att`/`pend`) at `cvS₅`, `l1_site : ne ∧ pend ∧ VContent` at
each loop-1 iterate (`l1Pre`/`l1Node`). The dropped clauses come from `CloseBase`: `shape` via
`RgStep` (`CloseBase.rgFeRest`, `RgStep.cvS₅`, `rangesInv_l1Iter`), `ne` via `EarFinish.path`/`sv_d`
(the retargeted entry's `topDepth` is `lowval`: `cvS₅_topDepth`), `att` via `Inv'.entries` +
`FinishTopOk.mid` (`att_of_finishTop`), `pend` via the open path and `RangesInv.processed`
(`CloseBase.path_pend`). `CloseContent` is checked field by field in `checks/WalkInvCheck/Ranges.lean`
(`checkContent`, labels `ranges.content_*`): 0 violations over seeds 0..3000 × both modes; its only
consumer-side hypothesis is `closeBase_content : CloseBase → CloseContent`, which the backbone
induction (`WalkBackbone.lean`) discharges at the `finishEdge` sites. For `walk_rootsCover`, the P-site coverage is the theorem `finishP_ownership`
(`RangesOwned.lean`), stated at a `finishEdge` over `CloseBase` and the prefix-ownership invariant
`Owned σ nR n₀ n v d orig P sts s` (every processed edge of the root tree since `nR` is on the stack or
under the vertex item of a path vertex; strict ancestors' vertex items hold only edges before the next
path start `sts[k+1]`; `vertItem v` holds exactly edges of `n₀..n`; the entries above the `orig` old ones
hold only edges `≥ n₀` and the old ones do not bottom at `v`; entries bottom at visited vertices and
unvisited vertex items are empty), with `nR ≤ n₀` and `sts[k] ≤ n₀` on the path; the checker `checkOwned`
(kinds `own_*`) evaluates all clauses at every `finishEdge` pre-state, P site and vertex end, 0..400 × both
modes, 0 failures. Proof (standard axioms): the merge base `nxt` bottoms at `curV` (`condP`), so by `old`
it is a new entry and by `new` holds only edges `≥ n₀`, which makes `MergeBaseCover.base` vacuous for
`lo := min n₀ (n + 1)`; `child` is `cover` transported through the three states — `closeVert'_edges`
(type-1 vertex close: one entry holding the `c`/`py`/`vy` edge sets above an unchanged `base`, via
`unwrap_edges`, `Step.mergeTop`/`retarget` and `finishTop_edges`), the P `unwrap_edges`, and
`EarClose.base_edges`/`sub_cover` — with the `anc`/`sts` bounds excluding strict-ancestor vertex items
and `EarFinish.vert` placing `vertItem curV` on a `base` entry; the returning edge `σ[n] = o.e` is a
sub-edge (tree) or on the pushed edge entry (back). The enclosing induction is set up over the
depth-indexed `OwnedD σ sts origs P d n s` (`Owned` at every path vertex at once: `sts[k]`/`origs[k]`
are the schedule position and stack length at the entry of the path vertex of depth `k`; `cover` over
the root tree `sts[0]..n`; `hi`/`lo` the current vertex item holding exactly `sts[d]..n`; `anc`/`lo`
strict ancestors holding `sts[k]..sts[k+1]`; `new`/`old` for every `k ≤ d`; `vis`/`fresh` as before).
Its single step is the named theorem `finishEdge_ownedD` (`OwnSite` = the `CloseBase` exports
`hD`/`nodup`/`rgs`/`guards`/`book`/`pos`; `OwnedD` at `n` with `sts`/`origs` monotone along the path,
`sts[d] ≤ n`, `origs[d] ≤ origTstack`, the path below `d` avoiding `curV`, the path vertices in range
(`stackVerts[k] < nv`, `k ≤ d`) and `P (vertItem o.dest)` for a returning tree edge (the two
hypotheses added when it was proved; `cvOut` supplies them from `Keep`/`hanc` and the child's `Place`)
⊢ `OwnedD` at `n + 1` after `finishEdge`), proved (standard axioms). Returning edges: `EarFinish`
splits the pre-stack as `sub ++ base`, `Dec` (threaded through `closeEars`/loop 1/`mergeLate`/
`closeVert'`/`finishRest`) describes the post-stack as `N ++ Bt` with `Bt` an untouched suffix of
`base`, every edge of a new entry a sub-edge, `o.e`, an edge of an absorbed base entry (all bottoming at
`curV`, so above the `origs[k]` old ones) or below `vertItem curV`, and every such edge held by `N`;
`OFrame` keeps items below the old size unchanged; `OwnedD.of_split` rebuilds the nine clauses.
Boundary returns: `OwnedD.of_boundary`/`of_boundary'` — the popped entries' edges and `o.e` move under
`vertItem curV` (`Below_modify_ch_append`/`_set` for the `ch` writes of `finishBoundary`), other vertex
items and the base entries are unchanged. `checkOwned` evaluates the `OwnedD` clauses at every vertex
entry (`entry`), `finishEdge` pre-/post-state (`pre`/`post`, the latter at `n + 1`), P site and vertex
end, 0..400 × both modes, 0 failures. The induction itself is proved (`RangesCoverTree.lean`): the mutual
`cvTree`/`cvOuts`/`cvOut` (`RgTree` style) carry `Place`, `OwnedD`, the processed-prefix bound on the pushed
edges, the schedule/stack origins `sts`/`origs` and the open path through `walkOutPre` (`walkOutPre_ownedD`,
`walkOutPre_place`, `walkOutPre_hv`), the child walk and `finishEdge` (`finishEdge_ownedD`,
`finishEdge_place`, `finishEdge_vertCover`), discharge `FinishCover` at every `finishEdge` from
`finishP_ownership` (`OwnedD.toOwned` at the active depth) and `Place.pushVertR`, and build the
`CoverTree`/`CoverOuts`/`CoverOut` schedule (`coverOut_tree`/`coverOut_back`); the vertex end re-pushes
ownership (`OwnedD.pushVert`/`exit`), and `OwnedD.root`/`forest_rootsCover`/`walk_rootsCover'` run the
DFS forest from `RootState` (fresh vertex items childless, `RangesInv.root_append`). `walk_rootsCover`
(`WalkItemsWF.lean`) is that assembly; both of its steps are proved. The other step, `finishEdge_vertCover` (`hasVert` after
`finishEdge` ⊢ `VertCover curV`: every edge below `vertItem curV` is on some stack entry; checker
`checkVCover`, `own_vcover`, at every `finishEdge` post-state and vertex end with the vertex pushed,
0..400 × both modes, 0 failures), is proved (standard axioms) through the stronger `VSpan v s`
(`vertItem v` is below an item spanned by some stack entry), which every primitive of `finishEdge`
keeps: `VSpan.alloc`/`modifyVs`/`push*`/`retarget` (frame), `VSpan.mergeTop` (two entries, derived
from a nonempty result: `VSpan.mergeTop'`, `two_of_mergeTop_ne_nil`), `VSpan.maybeUnwrap` (the
unwrapped `nxt` item is not `vertItem v`: `Shape.vert` and `ty ≠ V`, so its children keep the witness),
`VSpan.finishTop` (`ItemFree`: the witness is not the closed item, which gains it as a child), the
merge loops (`VSpan.mergeLoop`: a merge loop ending nonempty had two entries at every merge,
`mergeLoop_nil`), `VSpan.loop1Body`/`closeEars` (`loop_pred`, `Step.iter'`, `loop1Cond` ⊢ two
entries), `VSpan.mergeLate`/`finishP`/`closeVert'` (`CloseVertOk.retarget.nonempty` ⊢ the two merges
had entries) and `VSpan.finishRest` (`hasVert = false`: the pushed vertex entry, merged onto a nonempty
stack, `finishP_ne_nil`); at the pre-state `hasVert` gives the witness from `EarFinish.vert`, the
boundary branch returns `false`, and `EarFinish.loops` gives the nonempty `feS₂`.

Public site exports (`RangesSites.lean`, importable by any file; nothing in it is admitted beyond the
inductions above): `CbTree`/`CbOuts`/`CbOut` (the `CsTree` shape with `CloseBase` as the conclusion)
state `CloseBase σ n D curV d o origTstack hasVert s` at the pre-state of every `finishEdge` call of
`walkTree`, proved by the mutual `cbTree`/`cbOuts`/`cbOut` from the same exports as `ccTree` (which it
reuses for the `CloseInv` component); `CbForest`/`forest_closeBase`/`walk_closeBase` run the DFS forest
(`walk_closeBase g tern vo eo` for `g.dfsForest vo eo`, from `walk_rootsCover'`). From a `CloseBase` at a
returning tree edge (`isTree`, `lowval d < d`), the range invariant is carried into `finishEdge`:
`rangesInv_l1Iter` gives `RgStep σ (n + 1) D curV s (l1Iter d o s k)` (so `RangesInv σ (n + 1) D` and
`Shape` at the `k`-th loop-1 iterate, when the first `k` loop conditions held; `RgStep.iter`,
`closeEarsAdj_of_frontier`) and `rangesInv_feS₂` the same at `feS₂ d o s` (`mergeLateAdj_of_frontier`);
`CloseBase.finishOk` derives `FinishOk` from the guards/book fields. Both are standard-axiom; the checker
evaluates the `RangesInv` clauses at exactly these states (`checkRI`, labels `rangesInv_l1Iter`/
`rangesInv_feS₂`, 0..400 × both modes, 0 failures). `CloseBase`/`CloseCtx`/`OwnedD` are unchanged.
No new admission was added to the completed local lemmas.

### 4.7 Backbone: one invariant, one induction (`checks/WalkInvCheck.lean`)

The per-layer inductions (ear `cTree`/`cOuts`/`cOut`, Ranges `dsTree`/`cvTree`/`cbTree`, R
`rrTree`/`rkTree`, ST `stWalk`/`stForest`) are to be replaced by ONE conjunction `WalkInv` carried
by ONE induction over `walkTree`/`walkOuts`/`walkOut`/`finishEdge`, every component lemma's
hypothesis taken from the conjunction at the same site. Before stating it, the conjunction is
checked executably: `lake build check_walkinv && ./.lake/build/bin/check_walkinv [lo hi]`
(default seeds 0..3000, both `ternarize` modes, `randGraph` of `checks/EarCheck.lean`) re-runs the
walk by one instrumented traversal around the library `finishEdge` (checked equal to `Graph.walk` on
every run: items and tstack length) and evaluates, per field with the smallest violating seed:

* `inv.*` (`Inv' d`, `WalkSpec.lean`) at every between-edge site;
* `ear.*` (`EarCtx.check` of `checks/EarCheck.lean`: `EarCtx` between edges, `EarFinish`/
  `Loop1Spec`/`EarLate`/`EarClose`/bottoms/`Frontier` at `finishEdge`, `compEnd`/`left` at the end
  of `walkTree`);
* `ranges.*` (`checks/RangesInvCheck.lean`'s `RangesInv`/`FinishAdj`/`CloseInv`/`CloseCtx` and the
  P / vertex / loop-1 sub-sites, `OwnedD`/vertex coverage; its `cand_*` observations excluded);
* `r.*`/`pieces.*`/`rskel.*` (`checks/RFinishEdgeCheck.lean`'s contract-B `EntryR` lines `preB`/
  `postB`/`stab`/`anc`/`disj`, the loop-1 emulation `RTop`/`RBranch`/`rclose`, `FinishRShape`,
  `RInvG` base frames, on blocks; `RSkel3` of every R item at every site);
* `st.*` (`CheckStSim.lean`'s m8–m11 reads/items/truncations at every site, m13–m17 contexts, the
  base/root frames, `finishBoundary`, and the per-segment `StLive ↔ openBlock` pairing §7 listed as
  unchecked: `st.live_cur`/`st.live_pre`/`st.live_end`/`st.live_lower`, `st.seg_lower`);
* the handoff candidates E1–E5 of §4.5 (`e1.entry`, `e2.bd_free`, `e3.stab_single`, `e4.span_root`/
  `e4.unique_parent`, `e5.*`).

Result (seeds 0..3000 × both modes, 6002 runs, 166 block runs): instrumentation matches, 0
violations in every field — after one correction: the original E3 `buried_vacuous` is false (6
violations, kernel counterexample `checks/WalkInvCheck/E3False.lean`), restated as `stab_single`
(§4.5 E3). `check_walkinv minimize <field> <seed>` edge-minimises a violating graph and
`check_walkinv small <field> <n>` searches small random graphs; both were used for the E3
counterexample.

**Stage 2 — the statement (`Spqr/WalkBackbone.lean`).** `WalkInv G t d s` is the conjunction at the
entry of `walkTree t d` (before `stackVerts[d] := t.v`), over one ghost record `TreeGhost`
(`g anc base bE sv sd pe` of the ear frame, the schedule `σ n`, the R data `dfs F`, the coverage
prefix `sts origs P X`, the st path `prev fs segs`): the DFS/frame facts (`Full`, `WF`, `Ends`,
nodup/bounds, edge completeness `comp`/`pe_anc`, `stackVerts`/`stackDir` sizes, `anc_sv`), the ear
context `EarCtx t.v d [] t.outs false …` with `Inv' d` and `Shape`, `RangesInv σ n d` with
`PostAt`/`AncPath`/`CloseInv`, `OwnedD` with its prefix bounds and `P`-freshness, the R context
`RCtx` (on a block, at a non-root entry `d = dp + 1`: `dfs.Spec`/`Rooted`, `Inv' dp` of the parent
state, `IsParent`, the out-lists of the subtree, the ancestor chain, `RInvTop` at the parent, the
frames `F` as `RInvG`, `RSkelInv`), and the st hypotheses of `StTreeP` (`segsStack segs`,
`SegRead`, `StItems (simBlocks …)`, the height bound). `WalkInvEnd G t d s s'` is the conjunction at
the exit: `TreeEnd`, `Inv'`/`Shape`, `RgS σ (n + |edgePostorder|)`, `CloseInv`, `Place (Pushed …)`/
`OwnedD`/`VertCover`, (block) `RWalk dfs F t.v d` + `BotKeep` + `RSkelInv`, and `StSim`.

* `WalkInv.sites` is the cross-layer plumbing: from the conjunction alone it produces every site
  hypothesis the component inductions were stated against — `EarTree` (`cTree`), `BookTree`
  (`bTree`), `GuardsTree` (`gbTree`), `FrontiersTree` (`frTree`), `CsTree` (`dsTree`),
  `CoverTree` (`cvTree`), `RgTree` (`scheduleTree`), `CbTree` (`cbTree`), `RSideTree` (`rsTree`).
  No layer takes an input that another induction carries.
* `walkTree_inv : WalkInv G t d s → wp (walkTree t d) (fun _ s' => WalkInvEnd G t d s s') s`.
  In this stage its body is the conjunction of the existing per-layer tree theorems (`cTree`,
  `invTree`, `rgTree`, `ccTree`, `cvTree`, `rrTree`/`rkTree`, `stWalk`) applied to `sites`; the
  second stage replaces them by one `walkTree.mutual_induct` whose `walkOut` step is assembled
  from the `finishEdge` primitives (`finishEdge_rangesInv`, `finishEdge_closeInv`/`*_closeAt`,
  `keepsR_finishEdge`/`finishEdge_rInvTop`/`finishEdge_rInvG_base`, `finishEdge_st`/
  `finishRet_frame_st`/`finishBoundary_st`, `ctx_step_*`), after which the per-layer inductions
  are deletable.
* **Stage 2b — the between-edge level.** `WalkInvOut G B v d outs₀ done rest hasVert n P s` is the
  conjunction at a between-edge site of `walkOuts v d` (after `done`, before `rest`; `n`/`P` the
  current schedule position and pushed set, `B` the R bottom): the frame facts of `COut`
  (`split`/`sorted`/`WF`/`Ends`/nodup/bounds/`comp`/`pe_anc`/sizes/`anc_sv`/`sv_d`), `EarCtx v d
  done rest hasVert …` + `Inv' d` + `Shape`, `RangesInv σ n d` + `PostAt σ n (edgePostorderList
  rest)` + `AncPath` + `CloseInv`, `Full g P X` + `P_past` + `OwnedD σ sts origs P d n` with the
  prefix bounds and `P`-freshness of `rest`/`v` (`Pv`/`Pe`/`Pcur`/`Pcur'`) + `VertCover`, (block)
  `ROutCtx` (`dfs.Spec`/`Rooted`, `dfs.outs v = outs₀` and the sub-out-lists, `AncChain`, `RWalk
  dfs F v d`, the frame bounds against `B`, `RSkelInv`), and the st data (`StPre … (done.map (·.1))
  hasVert`, the height bound, `base_out`, `segs.length = fs.length`) **plus the §7 pairing as
  fields**: `live_cur` (`∀ b, StLive g items new (openBlock g fs (DirsOf s d) (refOuts … done).1 ++
  [⟨b, [vertItem v]⟩ unless hasVert])` for the current segment `new`, quantified over the side `b`
  the vertex piece will get at the push — the piece is a single V item, so its side is irrelevant
  to `StLive`; checker mirrors both `b`) and `live_lower` (segment `segs[j]` of
  frame `k = fs.length - 1 - j` is `StLive` in `openBlock g (fs.take k) (DirsOf s k) segs[j].2`).
  `WalkInvOut.sites` produces `EarOut`/`BookOut`/`GuardsOut`/`FrontiersOut`/`CsOut`/`CoverOut`/
  `RgOut`/`CbOut`/`RSideOut` from the conjunction (`bOut`, `gbOut`, `frOut`, `dsOut`, `cvOut`,
  `scheduleOut`, `cbOut`, `rsOut`; the per-out edge completeness `comp_out` is derived from the
  list-level `comp` by `endsOut_wf` and the vertex-disjointness of distinct outs).
  `walkOut_inv : WalkInvOut … done (o :: rest) hasVert n P s → wp (walkOut v d o hasVert)
  (fun hv' s' => (hasVert → hv') ∧ ∃ hvF, WalkInvOut … (done ++ [(o, hvF)]) rest hv' (n +
  |o.block|) (Pushed g (P ∨ hv' ∧ · = vertItem v) o.verts o.edges) s')` is the step (in this stage
  still the conjunction of the out-level component theorems `cOut`, `kOut`, `invOut`, `rgOut`,
  `ccOut`, `cvOut`, `walk_full_aux`, `stOut_step`, `rrOut`/`rkOut` applied to `sites`; the next
  stage opens `walkOut` and assembles it from the `finishEdge` primitives), and `walkOuts_inv`
  threads it through `walkOuts` (post: `∃ done'`, `WalkInvOut … done' [] hv' (n + |edgePostorderList
  rest|) (Pushed g (P ∨ …) (vertsList rest) (edgesList rest))`, via `WalkInvOut.congr` and
  `Pushed.append`).
* **New named admission (backbone-specific).** `walkOut_stLive` (`Spqr/WalkBackbone.lean`): from
  `WalkInvOut … (o :: rest) hasVert n P s`, `wp (walkOut v d o hasVert)` of `live_cur` for
  `done ++ [o]`/`hv'` and `live_lower` unchanged — i.e. one `walkOut` keeps the per-segment
  `StLive ↔ openBlock` pairing. This is the §7 "StLive closes" content stated at the only site the
  backbone needs it (checker fields `st.live_cur`/`st.live_lower`, 0 violations on 6002 runs); the
  tree-level `live_end` (what `finishBoundary_stLive` needs at the child's boundary) is the `[]`
  instance of `live_cur` after `walkOuts_inv`, to be plumbed in the next stage.
* **Stage 2c — the tree-level glue.** `walkTree_inv` now runs through `walkOuts_inv`: its body
  is `unfold walkTree`, `WalkInv.toOut`, `walkOuts_inv`, then `WalkInvOut.exit_true`/`exit_false`
  (the tree-level composition `cTree ∧ invTree ∧ … ∧ stWalk` is gone). `WalkInv.toOut :
  WalkInv G (.node v outs) d s → WalkInvOut {G with sts := sts ++ [n], origs := origs ++
  [|tstack|]} |tstack| v d outs [] outs false n P {s with stackVerts[d] := v}` is the entry glue
  (the R part builds `ROutCtx` from `RCtx` via `AncChain`/`EntryR.single` of the frames, taken for
  now from `WalkInv.sites.rside`, i.e. the `rsTree` composition — to be replaced by direct `RCtx`
  reasoning when the mutual induction lands); `WalkInvOut.exit_true` (no vertex entry pushed) and
  `WalkInvOut.exit_false` (`setStackDir d true; pushVertTstack v d`, per primitive: `Step.pushVert`,
  `RgStep.pushVert`, `CloseInv.pushVert`, `Place.cons_fixed`/`Place.fresh`, `OwnedD.pushVert`,
  `vertCover_pushVert`, `RWalk.pushVert` under `rSide_vertFree_site`, `botKeep_cons`,
  `StRead.pushEntry`/`StItems.pushEntry`/`DirsOf_push`, and the new `StLive.cons_vert`: pushing a
  `V` entry keeps a segment live since everything strictly below a `V` item hangs under it) are the
  exit glue. To make the glue close, `WalkInv` gained `segs_len` and the lower-segment pairing
  `live_lower`, and `WalkInvEnd` gained `live` (the current segment is `StLive` in `openBlock g fs
  (DirsOf s' d) (refTree g t d …).1` — the `live_end` field `finishBoundary_stLive` needs at the
  child's boundary) and `live_lower`. `#print axioms`: `StLive.cons_vert`, `WalkInv.d_lt`,
  `WalkInvOut.exit_true` standard; `WalkInv.toOut`, `WalkInvOut.exit_false`, `walkTree_inv` reach
  `sorryAx` (`exit_false` only through `rSide_vertFree_site`; `toOut` through `sites.rside`).
* **Stage 2d — the induction skeleton.** `backbone : (∀ t d, BbTree t d) ∧ (∀ v d outs hv, BbOuts
  v d outs hv) ∧ (∀ v d o hv, BbOut v d o hv)` is one `walkTree.mutual_induct` with the three
  motives `BbTree t d` (`∀ G s, WalkInv G t d s → wp (walkTree t d) (WalkInvEnd G t d s)`),
  `BbOuts v d rest hasVert` (the `walkOuts` statement over `WalkInvOut`), `BbOut v d o hasVert`
  (the `walkOut` statement); the cases are the named site lemmas `bbTree_node` (= `toOut`, the
  out-list hypothesis, `exit_true`/`exit_false`), `bbOuts_nil`, `bbOuts_cons` (standard axioms),
  and `bbOut_step` — whose body in this stage is still the out-level composition `walkOut_inv`
  (the child hypothesis `BbTree child (d + 1)` is **not consumed yet**). `walkTree_inv`/
  `walkOuts_inv` are now `backbone.1`/`backbone.2.1`. Next: open `walkOut` in `bbOut_step` —
  `walkOutPre` site (`wp_walkOutPre` + the `walkOutPre_*` primitives), the child entry
  `WalkInvOut … (.tree e cls child :: rest) → WalkInv G' child (d + 1)` (the per-layer child-entry
  derivations currently inside `invOut`/`dsOut`/`cvOut`/`scheduleOut`/`rsOut`/`stOut_step`/
  `rrOut`), then `finishEdge` from its primitives consuming `WalkInvEnd` of the child.
* Admissions reachable from `walkTree_inv` at stage 2d (historical — ear-5/Ranges-3 have since
  proved several of these; the current set is in stage 2g below; 36 = the 35 below +
  `walkOut_stLive`):
  ear `ctx_step_tree_ret`, `walkTree_below_kept`, `walkTree_items_kept`, `ends_of_wf_boundary`,
  `earAt_tree_{bd_bridge,bd_comp,bd_side,bd_term,bottom,close,late,late_fo,loop1,loop1_side,
  loop1_touch,loops,lower,of_ctx}`, `tree_comp_shape`; Ranges `closeCtx_{bd_node,bd_vert,l1_site,
  p_site,v_site}`, `finishEdge_ownedD`, `finishEdge_vertCover`; R `rSide_{entry,vertFree,
  finish_content}_site`, `loop1_rBranch_fields_ctx`, `loop1_rBranch_mid_ctx`, `loop1_rTop_ctx`,
  `feS₂_top_entryR`, `closeVert_type1_rCloseShape`; ST `finishBoundary_stLive`. Not reachable:
  `loop1_rBranch_content`, `walkTree_rSide_spec`, `walk_items_rThreeConnected` (R-layer glue
  outside the tree theorems), `spqrTree_wf`, the final-items hypotheses of `spqrTree_pieceSep`
  (`walk_pieceInv`, §4.6), `walk_canonInv`/`walk_ternarize`
  (canonicity, §4.6), the planarity admissions.
* **Stage 2e — R context uniform in `d`.** Preparing the child-entry lemma (`WalkInvOut … (.tree e
  cls child :: rest) → WalkInv G child (d + 1)`): the R fields `WalkInv.r`/`WalkInvOut.r` are no
  longer gated on `d = dp + 1`; `RCtx` now carries `par : ∀ dp, d = dp + 1 → RPar dfs t dp s`
  (the former parent-side fields `inv_par`/`parent`/`chain`/`top`) and `root : d = 0 → (∀ v outs,
  t = .node v outs → dfs.depth v = 0) ∧ s.tstack = [] ∧ F = []`, and `ROutCtx` gains `root : d = 0
  → s.tstack = [] ∧ hasVert = false ∧ F = []` (true by `rootOut`: every root between-edge site has
  an empty stack and no vertex entry). `WalkInv.toOut` builds the root `ROutCtx` from `RCtx.root`
  and `walkOut_inv_r` keeps it through a root `walkOut` by `rootOut` + `rkRootOut` (`RSkelInv`) +
  `kOut`. Checker: the `r.anc`/`r.stab` mirrors are no longer skipped at `d = 0` and a new
  `r.root` field checks `tstack = [] ∧ hasVert = false` at every root out-site; still
  `violations=0`.
* **Stage 2f — opening `walkOut`: the `walkOutPre` site and the child entry.** `PreOut G B v d
  outs₀ done o rest hasVert n P s push hv₁ s₁` is the conjunction after `walkOutPre v d o hasVert`
  from a `WalkInvOut` at head `o`: the push flag (`push ↔ ¬hasVert ∧ lowval < d ∧ isType1`),
  `hv₁ = hasVert || push`, the exact state (`stackDir[d]` set, the optional `TEntry` for
  `vertItem v` pushed), and every layer's post-`walkOutPre` contract (`Inv'`/`Shape`,
  `RangesInv`, `CloseInv`, `Full` with `P ∨ (hv₁ ∧ · = vertItem v)`, `OwnedD` + `VertCover`,
  `RWalk`/`BotKeep`/frame bound under 2-connectivity, the ST `PreSpec`, `Keep (d+1) 0`).
  `WalkInvOut.walkOutPre_out : WalkInvOut … (o :: rest) … → wp (walkOutPre v d o hasVert) (fun hv₁
  s₁ => ∃ push, PreOut … push hv₁ s₁)` is the conjunction of `wp_walkOutPre`, `walkOutPre_inv`,
  `walkOutPre_ranges`, `walkOutPre_closeInv'`, `walkOutPre_full`, `walkOutPre_ownedD`,
  `walkOutPre_r` (its `VertFree` hypothesis from `rSide_vertFree_site`, the only admission it
  reaches), `walkOutPre_pre`, `keep_walkOutPre`. **Child entry** `PreOut.child`: for a tree edge
  `.tree e cls (.node c couts)`, after `firstOccurrence[d] := ne`, `WalkInv G' (.node c couts) (d +
  1) s₂` with the child ghost `G' = {G with anc := G.anc ++ [v], base/bE := the new tstack, pe := (·
  = e), n, F := (v, d, |tstack|) :: F, P := P ∨ vertItem v, fs := fs ++ [frame v done o], segs :=
  (new₁, refOuts … done ++ vertex piece) :: segs}` — ear by `ctx_init_child`, `Inv'`/`RangesInv` by
  `frame'`+`setSv`, `PostAt`/`AncPath` from `post_o`/`path_o` (`e` sits at `σ[n + |child
  postorder|]`), `OwnedD.entry`, `RCtx` from `RWalk.top` (`RPar` at `dp = d`, new frame `(v, d,
  |tstack|)` by `RInvTop.toG`), `StPre` from `PreSpec.read`/`tstack`, and the per-segment
  `StLive`: the new current segment is `live_lower 0`, obtained from `live_cur` at side
  `s₁.stackDir[d]!` + `StLive.cons_vert` when the vertex piece was pushed (hence the `∀ b` in
  `live_cur`). `#print axioms`: `PreOut.child` standard; `walkOutPre_out` reaches `sorryAx` only
  through `rSide_vertFree_site`. Not yet wired into `bbOut_step` (next: the `.back` finish and the
  child `WalkInvEnd` → `finishEdge` composition).

* **Stage 2g — `bbOut_step` from the `finishEdge` primitives (`Spqr/WalkBackboneStep.lean`).**
  `bbOut_step` no longer delegates to the out-level composition `walkOut_inv`: it opens `walkOut`
  (`walkOut_eq`, `WalkInvOut.walkOutPre_out` → `PreOut`, `walkOutRest`), and at the `finishEdge`
  site assembles the post-state `WalkInvOut` from the per-primitive lemmas. Shared helper
  `finish_core` (`CloseBase` + `Inv'` + `Full` + `OwnedD` + `Keep` in, one `wp (finishEdge …)` out):
  `finishEdge_step` (`Inv'`/`Shape`), `finishEdge_ranges`, `finishEdge_closeInv` with the `CloseCtx`
  from `CloseCtx.of_exports` + `closeBase_content`, `finishEdge_ownedD`, `finishEdge_vertCover`,
  `finishEdge_full`, `keep_finishEdge`. **Back edge**: `PreOut` normalised, `EarFinish` from
  `BookOut`, then `d ≤ lowval` (boundary: `ear_boundary`/`finishBoundary_st`/`keepsR_finishBoundary`,
  root frame `ROutCtx.root`) vs `lowval < d` (`ret_of_lowval_lt`, `finishOk_of_guards`,
  `finishEdge_frontier`, `earFinish_close_len`, `stRet_finish`, `finishEdge_rInvG_base` per frame,
  `finishEdge_rInvTop`, `finishEdge_bot_keep`, `keepsR_finishEdge_site`). **Tree edge**: child entry
  `PreOut.child`, the induction hypothesis `BbTree child (d+1)` gives the child's `WalkInvEnd`
  (`inv`/`shape`/`ranges`/`owned.exit`/`full`/`vertCover`/`r`/`st`/`segRead`/`vStart`/`live`/`keep`),
  `walk_keepsBelow` for `DirsOf`, the per-segment `live` at the boundary gives `finishBoundary_st`'s
  `StLive` (`openBlockP_snoc_bd`), and the same finish-site lemmas as the back edge. New frame
  field `WalkInvEnd.keep : Keep d 0 s s'` (and `Keep (d+1) 0` in the `BbOuts`/`BbOut` posts), used
  for all `stackVerts`/sizes/item-type facts below the site. `pushed_back_iff`/`pushed_tree_iff`
  re-shape `Full`'s predicate to `Pushed … o.verts o.edges`. `backbone`/`walkTree_inv`/`walkOuts_inv`
  moved to the new file; the out-level `walkOut_inv` stays (unused by the backbone) until stage 3.
  Still composed, not primitive: the ear conjunct (`earOut` = ear-5's `cOut`), the R export
  `WalkInvOut.rshape` (`FinishRShape` at the site, via `rrOut`), and `sites.rside` inside
  `WalkInv.toOut`; `walkOut_stLive` is still the `live_cur`/`live_lower` admission. Not yet carried:
  `CanonInv`/`CloseCanon`, `PieceInv`/`ClosePiece`, `CloseContent` as a field (it is produced inside
  `finish_core` through `closeBase_content`), the `ternarize` frame, ear-5's
  `OutFrame`/`CtxShape`/`TreeEndS`. `#print axioms`: `pushed_back_iff`/`pushed_tree_iff`
  `[propext]`; `PreOut.child` standard; `finish_core` reaches only `closeBase_content`;
  `walkOutPre_out` only `rSide_vertFree_site`; `bbOut_step`/`backbone` reach exactly 25
  admissions: ear `earAt_tree_{bd_bridge,bd_comp,bd_side,bd_term,bottom,close,late,late_fo,loop1,
  loop1_side,loop1_touch,loops,lower,of_ctx}`, `tree_ret_shape`; Ranges `closeBase_content`; R
  `rSide_{entry,finish_content,vertFree}_site`, `loop1_rBranch_{fields,mid}_ctx`, `loop1_rTop_ctx`,
  `feS₂_top_entryR`, `closeVert_type1_rCloseShape`; ST `walkOut_stLive`.

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


**Back-edge branch of `finishEdge_rInvTop` (`Proofs/RInvBack.lean`).** `finishEdge_rInvTop`
is now proved by dispatch on the edge kind: the back-edge branch is `finishEdge_back_rInvTop`
(standard axioms); the tree-edge branch `finishEdge_tree_rInvTop` (same hypotheses, plus
`kind ≠ .backEdge`) is reduced below to the admitted `finishEdge_tree_top_settled_first`. The dispatcher takes one extra call-site fact,
`hback : kind = .backEdge → hasVert = true ∧ s.tstack.length ≤ origTstack` (a back edge is type 1
with `lv < d`, so `walkOutPre` has pushed the vertex entry, and `walkOutRest` reads `origTstack`
with nothing pushed since), under which `RInvFront` is `RInvTop` of the whole stack and
`Frontier.base_disj` says no open entry owns the unprocessed back edge `o.e`. The branch itself
is per-primitive: `RInvTop.modifyVs_free` (`feS₀` rewrites the terminals of the unspanned Q item),
`RInvTop.pushEdge` (the new `(curV, lv)` entry owns exactly `o.e`, `edges_edgeEntry`),
`RInvTop.modify` (bookkeeping), and `RInvTop.finishP` for the type-1 P-check — when `condP`
fires the entry below the top starts at `curV`, so `RInvTop.unwrapNxt_exempt`,
`RInvTop.mergeTop_exempt` and `RInvTop.finishTop_exempt` only have to track whole-stack
disjointness (the unwrapped/merged/finished entries' edge sets are subsets or unions of the
replaced ones; `maybeUnwrapNxt_tstack` keeps the start of the entry below the top), and
`finishTail` is the identity for `hasVert = true`. No ownership fact is re-derived: the one
needed (`o.e` unowned) is read off `Frontier`; `RangesInv.processed` with `FinishAdj.pos` would
give the same fact schedule-agnostically.
**Per-node characterization (`RelabelSpec.lean`, proved in `RelabelMain.lean`) [proved].**
`relabel_node_spec : items.WF g → ∃ idx, idx rootItem = 0 ∧ RelabelIdx g items t idx ∧
∀ i < items.size, RelabelNode g items t idx i` (`t = relabelTree g items`) is the shared
foundation of everything below: `idx` is the preorder numbering (a bijection onto `[0, t.size)`,
`RelabelIdx`), and `RelabelNode` says, for item `i` numbered `idx i`, that `type`/`par`/`origId`
/`vertIndex`/`edgeIndex` are the item's, that `chRange`/`nvRange`/`neRange` have the item's
lengths (`|ch|`, `|nvList|`, `nEdges`), that `nodeVerts[nvSt + k]` is `nvList[k]` re-indexed
through `vertIndex`, and (`RelabelLayout`, for *some* `pos`, the `vertPos` scratch at the time
the node was laid out) that the children are `(ordered g i nvSt pos).map idx`, the node-edges and
the two adjacency rows of each node-vert are those of
`layoutNode (type i) (idx i) nvSt nvEn neSt neEn (edgeChildren g pos (ordered …))`, the k-th
V child gets `vertParNv = some (nvSt + |vs.1| + k)`, the parent's k-th virtual edge has
`twin = some (neSt of the k-th non-V child)`, the child's cap points back **when
`items.hasCap child`**, and children occupy contiguous preorder ranges
(`child_idx`, `subtree_end`).
The proof (`RelabelGhost.lean` … `RelabelProof.lean`) is a `wp` induction over `relabel`
on a ghost copy carrying the numbering `order`, with the per-call contract
`CallPre`/`CallPost` (append-only prefix agreement `Agree`, `Consistent` sizes, and
`NodeS`/`LowB` for every descendant), the children loop by `wp_forIn_inv` with `LoopInv`,
and one abstract `Step*` lemma per primitive of the node body (`entry_of_steps`).

**Tree-edge branch of `finishEdge_rInvTop` (`Proofs/RInvTree.lean`).** The branch is proved
modulo one admission, the `EntryR` of the single entry it builds on top of the stack
(`finishEdge_tree_top_settled`: after the P-check, in the first-edge case, and in the output).
Its `hasVert = true` half is proved: `closeVert` re-targets the top to `curV` (`retarget`), and
`vertFinish`/`finishP`/`finishTail` keep a top starting at `curV` (when the P-check fires the
merged top is `nxt`, whose `vStart = curV` is the P condition), so the output top is exempt
(`TopStart`, `finishEdge_tree_vert_topStart`, standard axioms). The first-edge case
`finishEdge_tree_top_settled_first` (`hasVert = false`) is proved (R-4) modulo the single named
admission `feS₂_top_entryR`: the top `c` of `feS₂` is `EntryR` when `c.topDepth = d` and
`c.vStart ≠ curV`. Everything else is bookkeeping, from the ear contract now carried by
`FinishRShape.ear` (= `FinishBook.ear`'s `EarAt`): `EarFinish.loops` gives
`lowval ≤ c.topDepth ≤ d`, so the `feS₂` top is settled; the P-check keeps a settled top
(`settledTop_finishP`: when `condP` fires the merged top starts at `curV`; `finishP_top`); the
first-edge push puts the exempt vertex entry on top (`topStart_pushVert`), and when loop 2 fired
(`feSingle = false`) the merged output top has `topDepth = min c.topDepth d < d`: loop 2's first
merge pulls `c` down to the entry below loop 1's output (`l2Cur_topDepth_le`, `iter_merge_eq`),
which tops out below `d` by loop 1's exit condition (`RInvG.closeEars`, `result_loop1Cond_iff`).
So the content proper is exactly `feS₂_top_entryR`, i.e. `c.topDepth = d` forces loop 2 not to
have fired and `c` is loop 1's output (the tree edge plus the S/P/R items closed by loop 1).
Checker: 278 of the 6414 sites have a non-exempt first-edge top (`stat:ptop-nonexempt`), none a
non-exempt output top (`stat:post-top-nonexempt`); all pass contract B.
Two invariants carry everything else. `RInvG s dfs v d n` is `RInvFront` positionally: the bottom
`n` entries are settled (depth-bounded, `vStart ≠ v` exempt) and the whole stack is edge-disjoint.
Loop 1 preserves it with `n = origTstack` because, by `Frontier.loop1`, each iteration touches
only the two or three entries above the base (`RInvG.loop1Body`, `RInvG.loop1`,
`RInvG.closeEars`; the primitives `RInvG.pushEdge/mergeTop/finishTop/unwrapNxt/retarget/pushVert`
take their exemption hypotheses only when the touched entry lies in the bottom `n`, so they are
vacuous here). `RInvH s dfs v d` (= `RInvG` at `n = length − 1`) says everything *below the top*
is settled. `RInvG.widen` passes from the first to the second after loop 1 given that the frontier
entries below the top are settled (`FinishRShape.settled`), and loop 2, `closeVert` (loop 3, the
two merges, `retarget`, the type-1 close) and the P-check preserve `RInvH` with no exemption
hypothesis at all — they modify the top only (`RInvH.mergeTop/finishTop/retarget` via
`RInvG.mono`), except the unwrap of the entry below the top (`FinishRShape.unwrap`: the type-1
`closeVert` unwraps an exempt entry; the P-check's target starts at `curV` by `condP`) and the
first-edge vertex push (`RInvH.finishTail`), which moves the old top below the new vertex entry
and is where the admission's part (a) is consumed; part (b) closes the output with
`RInvH.toTop`. The previous candidate hypothesis — "after loop 1 every frontier entry below the
top tops out above `d` and every `closeVert` merge lands in an exempt entry" — is false
(`RFinishEdgeCheck`, 602 sites: the child's vertex entry `(c, d+1)` with `c ≠ curV` sits below the
top after loop 1 and is the target of the second `closeVert` merge; it is settled, `EntryR.vert`,
not exempt), which is what forced the top-hole formulation. `FinishRShape dfs curV d o origTstack
hasVert s` (`pend`, `settled`, `unwrap`, `vert_own`) and the admission's part (a) are checked on
the 6414 finishEdge sites of `checks/RFinishEdgeCheck.lean` (`shape` lines, 0 failures);
`finishEdge_tree_rInvTop_of_top`, `RInvG.*`, `RInvH.*` have standard axioms;
`finishEdge_tree_rInvTop` and the dispatcher depend on `sorryAx` only through
`feS₂_top_entryR` (via `finishEdge_tree_top_settled_first`).

**`WalkTreeRReturnSpec`: induction design (frame proved; induction proved up to the side facts).** The walk induction
(`rgTree`/`rgOuts`/`rgOut` style, `RangesTree.lean`) must carry two facts through the child's walk
from `s` with `n₀ = s.tstack.length` and parent `p = stackVerts[d]`: the parent's base positionally
settled, `RInvG s' dfs p d n₀` (bottom `n₀` entries unchanged as `TEntry`s and `EntryR` at
`(p, d)`, whole stack disjoint), and the child's own `RInvTop s' dfs c (d+1)` (vacuous on the
base, whose entries top out at depth `≤ d`). At the child's `finishEdge` sites the second fact is
`finishEdge_rInvTop` (with `RReturn` of the grandchild giving `RInvFront`); the first is the frame
statement `finishEdge_rInvG_base` (`Proofs/RInvBase.lean`, proved, standard axioms): a
`finishEdge` at `(curV, d)` with `n₀ ≤ origTstack` preserves `RInvG dfs p dp n₀` for any
`(p, dp)`. Its content is purely positional — the `RInvG.*` primitive lemmas of
`Proofs/RInvTree.lean` take their exemption hypotheses only when the touched entry lies in the
bottom `n₀`, so what is needed per primitive is a stack-length bound (`n₀ + 2 ≤ length` for
merges, `n₀ + 1` for top-only steps): loops 1–2 from `Frontier.loop1/2` (`origTstack ≥ n₀`),
loop 3 from its own guard (it fires only above `origTstack + 3`, `loop3_iter_len`), the
`closeVert` unwrap/merges/retarget/close and the P-check from `hclose : origTstack + 3 ≤ length`
after loop 2 (the child's `c :: mid ++ [py, vy]`; `FinishGuards` gives this only for type-2
edges with `hasVert`, so it is a separate hypothesis of the frame, discharged at the call site
from the ear layer) together with `hasVert = true → n₀ + 1 ≤ origTstack` (the vertex entry
`(c, d+1)` was pushed after the descent); on the back-edge branch `hback` (the vertex entry is on
top of the whole base) puts the push and the P-check above `n₀ + 1`. The positional primitives
(`RInvG.pushVert_pos`, `mergeTop_pos`, `retarget_pos`, `mergeLoop`, `mergeLate_pos`,
`closeVert_pos`, `finishP_pos`, `finishTail_pos`, `finishRest_pos`) are the `RInvH.*` lemmas
with the exemptions replaced by length bounds. The contract was dump-checked before the proof
(`checks/RFinishEdgeCheck.lean`, `base` lines: at every `finishEdge` site, for every ancestor
frame `(p, dp, n₀)`, bottom `n₀` unchanged, non-exempt `EntryR` before and after, and the three
side conditions; seeds 0..300 × both modes + 6000 random multigraphs, 19058 frame checks,
0 failures). `Sim.walkTree` cannot be used for the base directly: its precondition `GuardsTree`
is on the unlifted run. Two side facts enter as hypotheses, not admissions: every out-edge inside
a non-root subtree of a block returns (`∀ u ≠ dfs.root, ∀ o ∈ dfs.outs u, o.cls.lowval
(dfs.depth u) < dfs.depth u`; `finishEdge_rInvTop` covers returning edges only, the boundary
branch `finishBoundary` runs at the root), and the ancestor chain `hanc` is re-established at
each `walkTree` entry from `stackVerts.set! (d+1) c` and `dfs.Spec.depth_parent`.

**`WalkTreeRReturnSpec`: the induction (`Proofs/RInvWalk.lean`).** `rrTree`/`rrOuts`/`rrOut`
(mutual, `invTree`/`invOuts`/`invOut` style, consuming `walkOutPre_inv`/`invTree`/`finishOk_of_guards`/
`finishEdge_frontier`/`finishEdge_inv` and the ear bookkeeping `BookTree`) carry `RWalk dfs F v d`
— a list `F` of ancestor frames `(p, dp, n₀)` each held positionally as `RInvG dfs p dp n₀`, plus
`RInvTop dfs v d` for the vertex being walked — through the walk of `v`. At a `finishEdge` site
every frame is preserved by `finishEdge_rInvG_base` (the side conditions `n₀ ≤ origTstack`,
`hasVert → n₀ + 1 ≤ origTstack` come from `BotKeep` — the bottom `n₀` entries of the stack kept as a
list, `finishEdge_bot_keep` — and the walk-level length bounds; `hclose` from `EarFinish.close`,
`earFinish_close_len`), and `RInvTop` at `(v, d)` by `finishEdge_rInvTop` with `RInvFront` taken
from the child's frame `(v, d, origTstack)` pushed onto `F` for the child's walk. Entering a child
`(c, d+1)` from a parent settled on the whole stack (`RInvTop (p, d)`), the new `RInvTop (c, d+1)`
is the parent's `RInvTop` restricted to entries not starting at `p` (those top out at `≤ d`, a
side fact) and transported across `stackVerts.set! (d+1) c` (a second side fact).
`walkTree_rReturn` instantiates `F = [(stackVerts[d], d, s.tstack.length)]` and yields `RReturn`
(proved; its only non-standard axiom use is through `finishEdge_rInvTop`'s admitted
`finishEdge_tree_top_settled`). The side facts are bundled as `RSideTree`/`RSideOuts`/`RSideOut`
(hypotheses of the induction, threaded like `GuardsTree`): the ancestor chain `AncChain` at every
entry, `EntryR` stability under `stackVerts.set! d v`, no entry starting at the parent above
`d - 1`, every out-edge of a non-root vertex returns (`lowval d < d`; its DFS-layer content is
proved as `DfsData.Spec.outs_lowval_lt`, R-5, standard axioms; likewise the entry chain is
`ancChain_child` given `d + 1 < stackVerts.size`, and `EntryR` stability is `EntryR.set_stackVerts`
given `topDepth ≠ d + 1`, both R-5), `VertFree` (no entry owns an
edge below `vertItem v`) before a vertex push, and `FinishRShape` at tree-edge sites. R-5:
`walkTree_rSide` is proved (`Proofs/RSide.lean`) by the mutual induction `rsTree`/`rsOuts`/`rsOut`
mirroring `rrTree` (it runs `rrOut`/`rrTree` alongside for the `RWalk` at the next site, `invOut`/
`invTree` for `Inv'`/`Shape`, `frOut` for the site's `Frontier`, `kOut` for the layout): the chain
is `ancChain_child` with `d + 1 < stackVerts.size` from `ancChain_lt_size` (the `d + 2` distinct
vertices `stackVerts[0..d], c` of `g`), `ret` is `DfsData.Spec.outs_lowval_lt` with the parent of
`v` from `AncChain.parent`, the children's `dfs.outs` come from `DfsTree.Sub` (`ofForest_outs`).
The R content is isolated as three named site admissions: `rSide_entry_site` (`EntryR` stability
after `stackVerts.set! (d+1) c` and the parent-start top bound `≤ d`, at a child entry with the
parent's `RInvTop`), `rSide_vertFree_site` (`VertFree v` under `VertBook v false` + `RWalk`), and
`rSide_finish_content_site` (`FinishRShape`'s `settled`/`unwrap`/`vert_own` at a tree-edge site
given the child's `RWalk` with the parent frame, `FinishGuards`/`FinishBook`/`Frontier` and the
parent chain; `rSide_finish_site` proves `pend` from `EarFinish.q_root`/`q_free` via
`FinishBook.ear` and takes `ear` from it). All six side facts are
dump-checked at every entry/site of seeds 0..400 × both modes + 6000 random multigraphs
(`checks/RFinishEdgeCheck.lean`: `anc`, `stab`, `entry`, `ret`, `vertown`, `FinishRShape` lines;
0 failures). Statement change (R-5, non-proof-target): `walkTree_rSide` no longer takes
`FrontiersTree` (derived by `frOut`) and takes `g.WF`, `stackVerts.size = g.nv` and
`∀ t', t'.Sub (.node c outs) → dfs.outs t'.v = t'.outs` (in place of `outs = dfs.outs c`), all
root-threaded through `rkRootOut`/`rkRootOuts`/`rkRootTree`/`rkForest` (`RootState.sv`,
`ofForest_outs`). Statement correction (R-4): the earlier form
derived `BookTree ∧ RSideTree` from `Inv'`/`Shape`/`GuardsTree`/`FrontiersTree`/DFS facts alone,
which is under-hypothesized — `BookTree` at a non-root entry is ear content that only the
root-threaded `bTree`/`walkTree_book` produces (`gbTree` goes the other way), and `AncChain` for
`k < d` constrains `stackVerts` below `d`, which none of those hypotheses mention. `walkTree_rSide`
now takes `BookTree (.node c outs) (d+1) s` and the parent's chain (`∀ k ≤ d, Anc stackVerts[k]
stackVerts[d] ∧ depth stackVerts[k] = k`) as hypotheses and yields `RSideTree` only; the root glue
(`rkRootOut`, `Proofs/RSkelRoot.lean`) supplies both (`BookOut`'s `BookTree`, and `depth root = 0`
threaded as `hd0` from `DfsData.ofForest_depth_root`). The standalone form is kept as
`walkTree_rSide_spec` (admitted, consumed only by `walkTreeRReturnSpec`, whose statement
`WalkTreeRReturnSpec` is fixed; nothing else uses it). Of `RSideTree`'s parts, the chain at
`(c, d+1)` (`hchain` + `hp` + `Spec.depth_parent`), `FinishRShape.pend` (`EarFinish` `q_root`/
`q_free`) and `FinishRShape.ear` (= `FinishBook.ear`) are bookkeeping; `ret` is a DFS-layer fact
(`Spec.cls_*` + 2-connectivity: no bridge/component out-edge at a non-root vertex); `stab`, the
parent-top bound, `VertFree`, and `FinishRShape.settled`/`unwrap`/`vert_own` are content. Statement correction: `WalkTreeRReturnSpec`
quantified `s.Inv' D` over an unconstrained `D`, which nothing can consume (the walk at `d + 1`
needs `Inv' d`); it now takes `s.Inv' d`.

Two facts about the inputs that `Items.WF` does *not* give, found while stating this:

* **Orientation of R children (`Items.ROriented`).** `layoutNode .R` counts an edge child
  `(a, c)` at bounds slots `2 a + 2` and `2 c + 1` (and the cap at `2 nvSt + 2`, `2 nvEn - 1`),
  so each node's rows start at `2 neSt` and the last bound equals `adjDat.size` only if every
  edge child has `pos a < pos c`, i.e. is oriented along the node-vert order.
  `Items.WF` only has `q.1 ≠ q.2`, so `relabel_adj_spec` (and `SpqrTree.WF`'s `adj_bounds`
  clauses / `Shape .R`'s `p.1 < p.2`) take `items.ROriented g` as an explicit hypothesis;
  it is discharged by the ST layer's `Items.StNumbered` (`StSpec.lean`, every edge oriented
  low → high in `vertList` order).
  The C++ (`spqr_tree.hpp`, the `nvs[0]`/`nvs[1]` counting loop) relies on exactly the same
  orientation: it sets `nvs = {vert_pos[vs[0]], vert_pos[vs[1]]}`, `assert(nvs[0] < nvs[1])`,
  and increments `bounds[2 * nvs[0] + 2]` / `bounds[2 * nvs[1] + 1]`.
* **Q block-roots as children of S/P/R.** `Items.WF` allows a Q item with children
  (`Shapes.q_children`) whose parent is an S/P/R node; then `hasCap Q = false`, the parent's
  virtual edge gets `twin = some neSt_Q`, but `nodeEdges[neSt_Q]` is the Q's own first child
  edge, so the back-pointer (and `SpqrTree.WF.twin_invol`) fails.
  `RelabelLayout.twin` therefore guards the child side by `hasCap`; `relabelTree_wf` as stated
  needs `Shapes` to exclude this (a Q with children has a V/F parent) or the same hypothesis.

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
`PlanarSpec.neRotAdj_segment` (proved in `PlanarRotSpec.lean` from the fold characterization
`planarRelabel_rot_spec`, itself proved by the `RotInv` induction of `PlanarRotFold.lean`) reads no adjacency row (it is about `layoutRot`, which does not look
at the `Layout`); the layout-side facts its transport needs are the `Local` sizes and
`bound_last` (the `adjBounds.extract 1 …` pushed by `planarRelabel` has `2(nvEn - nvSt)` entries
and ends at `2·neEn`).

The structural half of this is done through the per-node interface of `RelabelSpec.lean`
(`RelabelIdx`: the preorder index `idx` is a bijection with root `0`, agrees with
`vertIndex`/`edgeIndex` on `V`/`Q` items, CSR zero/last facts; `RelabelNode`: each item's slot,
ranges, `node_verts`, `child_par`, and a pos-dependent `RelabelLayout` whose edge/adjacency rows
are `layoutNode`'s). `RelabelOwn.lean` takes `relabel_node_spec` (**[proved]** in `RelabelMain.lean`, the glue
`∃ idx, …` for the recursive fold) and proves `relabelTree_own :
Items.WF g items → Items.ROriented g items →
(relabelTree g items).Bijections ∧ .Ownership ∧ .Twins` **[proved]**:
* `Bijections` from `Items.Tree` only (types of `1+v`, `1+nv+e`, `RelabelIdx.vert_index/edge_index`,
  `RelabelNode.orig`);
* `Ownership` from `Items.Tree`, `Endpoints.vs_shape`, `Shapes.*_shape`/`i_o_leaf`/`q_children`
  (the per-type `layoutNode` hypotheses), the `layoutNode` edge records (`layoutNode_edges`, read
  off `LayoutShape.lean`'s `layoutNode_*_eq`/`run*_edges_get`/`run*_spec`; the R edges are
  re-derived there without `LayoutShape.run_edges`'s endpoint hypothesis so that `Twins` does not
  need `ROriented`), and the one extra hypothesis `Items.ROriented` (`ne_nvs` for R nodes needs
  the ordered endpoints to lie in the node-vert list, which is the st-ordering fact of §7, not part
  of `Items.WF`). Four item-level facts this needs live in `Items.WF` for the walk to discharge:
  `Endpoints.nv_nodup` (node-vert lists have no repeats; `I`/`P` nodes have `u ≠ v`),
  `Shapes.q_leaf_of_node` (a `Q` child of a node is a leaf: block-root Qs hang under `F`/`V`),
  `r_shape`'s last clause (non-V children of an R node have two endpoints, so with
  `child_vs_in_parent` both ends of its virtual edges are node-verts), and `q_children`'s `v < nv`
  for the `[c, vertItem v]` block root;
* `Twins` from `Items.Tree` and `Shapes.q_leaf_of_node` (a capped child has at least one edge),
  via `RelabelLayout.twin` both ways, `child_cap_twin_none`, injectivity of `idx` and disjointness
  of the `neRange`s.
`RelabelAdj.lean` adds the CSR facts of `relabel_adj_spec` on the same hypothesis:
`relabelTree_adj : Items.WF g items → Items.ROriented g items →
(∀ n < size, adjBounds[2 nvSt n] = 2 neSt n) ∧ adjBounds[2 |nodeVerts|] = |adjDat|` **[proved]**.
Every item's `nodeLayout` satisfies `Layout.Local` (`layout_local`: the `local_*` instances of
`LayoutShape.lean`, with their hypotheses read off `Items.WF` through `WF.layout_hyps`, `WF.nvList_V`
— V items have no node-verts — `WF.nvList_Q`, `WF.nEdges_pos`, and for R the orientation
`edgeChildren_bounds` from `ROriented`), so `RelabelLayout.adj_bounds` at `j = 2 nVerts` plus
`Local.bound_last` give `adjBounds[2 nvEn] = 2 neEn` whenever `adjBounds[2 nvSt] = 2 neSt`
(`adj_bound_step`; a node without node-verts has no node-edges). Since `nvSt (n+1) = nvEn n`,
induction over the node index gives `adjBounds[2 nvBounds[n]] = 2 neBounds[n]` for all `n ≤ size`
(`adj_bounds_start`), which is the first clause for `n < size` and the second at `n = size`
(`RelabelIdx.nv_last`/`ne_last`/`adj_dat_size`). The Q case needs `Endpoints.q_vs`'s corrected
clause: a block-root Q (a Q with children) records only its upper endpoint, `vs.2 = none` — the
walk sets this for every block boundary, not only self-loops, and the statement `vs.2 = none ↔ loop`
was false for a bridge block (`[I, vertItem v]` under a Q with `vs = (some u, none)`).
`RelabelWF.lean` assembles `SpqrTree.WF` on the same hypothesis, field by field
(`RelabelAll.wf_tree`), closing `relabelTree_wf : Items.WF → Items.ROriented →
(relabelTree g items).WF` **[proved]**. Two `Spec.lean` clauses had to
be corrected (both were false of the output, i.e. statement bugs):
* `Preorder.only_root_F` now has the bound `i < size`: `type` reads `nodeTypes[i]!`, which is
  `default = .F` out of range. Proved as `only_root_F` (from `Items.Tree.type_F_iff` through
  `RelabelNode.type` and the `idx` bijection).
* `WF.adj_incident` listed each node-edge incident to `nv` once, but a loop `(nv, nv)` (a Q node of
  a self-loop) occupies both rows `2 nv` and `2 nv + 1` of the CSR. It now states the rows with
  multiplicity — the row pair of `nv` is a permutation of `filter (nvs.2 = nv) ++ filter (nvs.1 = nv)`
  over all node-edges (`adj_incident'`) — which coincides with the old form whenever no node-edge is
  a loop.
`ROriented` is needed for R nodes (`Ownership.ne_nvs`, `Shape`); `StOriented.lean` derives it from
the st-ordering, `Items.rOriented_of_stNumbered : StNumbered → WF → ROriented` (`vertList =
nvList` under `Items.Tree`), so `walk_items_rOriented' : ROriented (walk …).items` follows from
`walk_st`; `spqrTree_wf'` uses it directly (`spqrTree_eq` lives in `WalkWF.lean` below the st
layer; `walk_st` takes `Items.WF` as a hypothesis, supplied by `WalkItemsWF.walk_items_wf`).
The other fields: `Sizes` is `RelabelIdx.sizes`; `Preorder` from `RelabelLayout.child_idx`/
`subtree_end` (children are numbered `idx i + 1`, then after the previous child's subtree —
`chainEnd`; by induction on `|Items.desc|` every subtree is non-empty and contained in the
parent's, `subtree_props`), giving `children_pairwise`, `subtree_eq` (`children_sum`), and
`parent_eq_iff` (`t.parent n = some (idx i)` iff `n` is a child's index); `Shape` from
`LayoutShape.shape_*` transported along `skeleton_eq` (the node's global `skeleton` is its
`nodeLayout`'s, by `RelabelLayout.edge_nvs`), with the R `Nodup`/non-parallel hypotheses from
`r_shape` and injectivity of `pos` on `nvList` (`PosOK`); `adj_bounds_mono`, `adj_dest`,
`adj_incident'` by locating the node that owns the row/entry (`nv_locate`/`ne_locate`), identifying
the global row with the local one (`global_bound`: `adjBounds[r] = rowBound r` for
`2 nvSt ≤ r ≤ 2 nvEn`, from `relabelTree_adj`'s first clause and `RelabelLayout.adj_bounds`;
`global_row`), and using `Layout.Local`'s `bound_mono`/`adj_dest`/`adj_incident_lo/hi`; node-edges
of other nodes have no endpoint in this node's `nvRange` (`foreign_ne`, from `Local.ne_nvs` and
disjointness), which restricts the global filter to the node's own segment (`row_filter`).

`RelabelRep.lean` transports `Represents` along the same interface **[proved]**: `relabelOK_of_wf` packages `relabel_node_spec`'s witness as `RelabelOK`
(`ridx`, `node`, the `Items.WF` fields), and the fields are
* `nv`/`ne` (`RelabelIdx.Sizes`), `interior`, `canonical` (`Items.Endpoints.interior`,
  `Items.Shapes.canonical` along `idx`/`vertIndex`), `nv_orig_inj` (`Endpoints.nv_nodup` through
  `node_verts`, `nvList` injective on original ids) — from `Items.WF` alone;
* `q_endpoints`, `separation`, `twin_glue`, and the output-level `r_three_connected` additionally
  use three `Items.WF` clauses added for this transport (formerly a separate `Items.RepOK`):
  `Endpoints.q_root` (a block-root Q with children `[c]` is a self-loop and with `[c, vertItem w]`
  has edge `{u, w}` — needed to read the Q's node-verts as its edge; in the walk, `[c]` is written
  only by the self-loop branch of `finishEdge` and `[c, vertItem w]` by the bridge/block branches
  with `w` the tree child, so this is the same tstack-span fact as `walk_q_children`),
  `Shapes.o_parent` (an `O` item is allocated only in the self-loop branch, as the sole child of
  its Q — its cap is a loop `(x, x)`, used by `twin_glue`/`separation` and the R virtual-edge
  endpoints), and `Shapes.s_order` (the positional form of `s_shape`: the S children in `ch`
  order are the path `u → xs → v`, so that `skeleton` of an S node is a path and its virtual
  edges are the consecutive pairs; `s_shape` only gives a permutation; `finishTstackTop` writes
  `ch` as the ear's span in walk direction, and `walk_st'`'s reference order fixes it). The three
  clauses are checked on the walk's output by `lake build check_repok` over `gen.py` seeds
  (0..400, 0 violations). `twin_glue` is
  `RelabelLayout.twin` + `edge_nvs` on both sides plus `Endpoints.child_endpoints`
  (`virt_glue`/`cap_orig`); `r_three_connected` is a pure transport
  `Items.RThreeConnected g items → (relabelTree g items).r_three_connected`, where
  `Items.RThreeConnected` (R items' `rSkeleton` — cap plus child virtual edges, read through the
  node-vert positions — is `ThreeConnected`) is the item-level statement of §4.5
  (`items_r_three_connected`, admitted there), not proved here.
`relabelTree_represents' : Items.WF → Items.RThreeConnected → Represents` and
`relabelTree_represents_of_r` (same with the output-level R clause, the hypothesis
`spqrTree_r_three_connected` provides). `Correctness.relabelTree_represents` now takes
`Items.RThreeConnected` as its explicit hypothesis and is `relabelTree_represents'`;
`spqrTree_represents` (under `g.WF`/`OrderOK`) is `relabelTree_represents_of_r` on `walk_items_wf`
and `spqrTree_r_three_connected`, so its admissions are those of `walk_items_wf` (`walk_ranges`,
`walk_sides`) and `spqrTree_r_three_connected`.

`RelabelSt.lean` proves `relabel_st : Items.StNumbered → Items.WF g → (relabelTree g items).StOrder`
**[proved]** from the same package, see §7.3.

## 6. Lean plan (what is proved where)

| statement | file | status |
|---|---|---|
| graph / DFS / items / walk / relabel definitions | `Graph.lean … Build.lean` | def |
| output spec `SpqrTree.WF`, `Represents` | `Spec.lean` | def |
| phase-2 contract `Items.WF` | `ItemSpec.lean` | def |
| 1.1, 1.2 DFS spanning + lowpoints | `Proofs/Dfs.lean` (`dfsForest_spanning`, `dfsForest_wf`, `classify_child_*`) | proved |
| DFS endpoints `∀ t ∈ g.dfsForest vo eo, t.Ends g` (`dfsForest_ends`: tree edges join child and parent, back edges the vertex and `dest`) | `WalkInv.lean` (`dfsVisit_v`, `dfsVisit_ends`; `adjacency_ends` in `Proofs/Dfs.lean`) | proved (standard axioms) |
| 2.1 blocks ↔ `lowval ≥ d` branches (`sameBlock_iff`, `blockRoot_cut`, `ret_child_sameBlock`) | `Blocks.lean`, `Proofs/Blocks.lean` | proved |
| Facts A–C (`sepPair_comparable`, type-1 class / above-between / sorted prefix-suffix, `type2_first_out`) | `SepPair.lean`, `Proofs/SepPair.lean` | proved |
| Fact D (laminar intervals: `type1_class_interval`, `type2_class_interval`, `type2Block_laminar_block`, `type2Block_laminar_type2Block`) | `Proofs/Postorder.lean`, `Proofs/Interval.lean` | proved (S-caveat exceptions) |
| 4.5 exhaustiveness: separation pairs of a block = type-1 ∪ type-2 pairs (`sepPair_iff`, `three_connected_of_no_split`) | `SepPairExhaust.lean`, `Proofs/SepPairExhaust.lean` | proved |
| 5 skeleton/contraction: `Graph.contract` of a laminar family of 2-attached pieces; `sepPair_contract_lift`, `sepPair_contract_of` (non-terminal pairs), `threeConnected_contract_iff_dfs` (R skeleton 3-connected ↔ no type-1/type-2 pair among non-terminal skeleton vertices), `TwoAttached.sepClass_mem` | `Contract.lean`, `Proofs/Contract.lean` | proved |
| walk pieces ↔ separation classes: a nonempty proper `TwoAttached` set of a block is a union of `{u,v}`-classes with `{u,v}` a separation pair, or it / its complement is a single `u–v` edge (`twoAttached_union_classes`, `twoAttached_type1_or_type2`) | `Proofs/SepClasses.lean` | proved |
| ear-structured walk (`descend`/`ascend` over chain `Frame`s) | `Ear.lean` | def |
| `walkEarTree = walkTree` (`walkEarTree_eq_walkTree`, `walkEar_eq_walk`) | `EarSpec.lean` | proved |
| span discipline / placement (`Place`: each item id placed ≤ 1 time over spans + ch lists; `walk_place`, `walk_ch_nodup`, `walk_parent_unique`, `root_no_parent`) | `WalkPlace.lean` | proved; `walk_tstack_nil`, `walk_covered`, `walk_root_children`, `walk_reach` moved to `WalkCover.lean` |
| coverage + reachability half of `Items.Tree` (`WalkState.Full` = exact placement + `Acyc`; `walk_tstack_nil`, `walk_covered`, `walk_root_children`, `walk_reach`, `walk_tree`) | `WalkCover.lean`, `ItemAcyc.lean` | proved from `walk_sides : SidesForest forest (WalkState.init g tern)` (§4.4; `walk_sides`/`walk_full`/the consumers now live in `EarWalk.lean`, `walk_sides` proved from `RootsBook forest (init)` = `rootsBook_of_state`, whose only admission is `walkTree_book` (restated, true at every root; the frame `walkTree_frame` is proved in `EarFrame.lean` under `Types g s`); new hypotheses `hends : ∀ t ∈ forest, t.Ends g` and `hecov : ∀ e < ne, e ∈ forest.flatMap edges` (forced by `walkTree_ear_false`, §4.2b; `comp_of_forest`)); its per-`finishEdge` content `FinishSides` is derived (`finishEdge_sides`, `walkTree_sides` in `EarSides.lean`, standard axioms), `RootOK` after each root and the forest threading are derived (`walkTree_rootOK`, `sidesForest_of_roots`, `walk_sides_of_roots` in `EarRoot.lean`, standard axioms), leaving exactly `RootsBook forest (init)`: `GuardsTree`/`BookTree`/`Inv' 0` at each root, derived in `EarWalk.lean` (`RootState.book`/`RootState.step`/`rootsBook_of_state`) from `walkTree_book` (admitted) + `walkTree_frame` (proved, `EarFrame.lean`) |
| linear phase 2/3 refinements `walkFast` / `relabelTreeFast` (`CatList` spans, `Array` tstack, ticks): `walkFast_items`, `relabelTreeFast_eq`; `spqrTree_eq` routed through them | `CatList.lean`, `Refine.lean`, `WalkFast.lean`, `RelabelFast.lean`, `WalkWF.lean` | proved |
| step bounds: `walk_ticks_le : (g.walkFast tern (g.dfsForest vo eo)).ticks ≤ 47·(nv+ne)` (via `walk_ticks_le_forest`, `dfsForest_size_le` from `dfsForest_spanning'`); `relabelRun_sizes_le` (`Items.desc` cardinalities), `relabel_ticks_le : ticks ≤ 288·(nv+ne) + 6` | `WalkCost.lean`, `ItemTree.lean`, `RelabelCost.lean` | proved; the relabel bounds take `Items.Tree` as hypothesis (`relabel_ticks_le'` discharges it with the admitted `walk_items_wf`) |
| frame rule `walkTree_local` via `Lifts`/`Sim` simulation (`Sim.closeEars`, `Sim.mergeLate`, `Sim.finishRest`, `Sim.finishBoundary` proved) | `Sim.lean`, `Frame.lean`, `EarSpec.lean` | `Sim.closeVert`, `Sim.finishEdge`, `Sim.walkTree` proved; `walkTree_local` takes `GuardsTree` as hypothesis; the public `walkTree_guards` (`EarShape.lean`) is stated at a root with `walkTree_book`'s hypotheses (incl. `hcomp`, §4.2b) and proved from it; `GuardsTree` is derived from `BookTree` (`EarShape.lean`: `finishGuards_of_ear`, `walkTree_guards'`, standard axioms), so the remaining content is `walkTree_book` |
| typing/allocation part of `Items.WF` (`Items.Tree` sizes/types, I/O leaves, `vs_shape`, `vs_lt`): `walk_typing` | `WalkTyping.lean` | proved (`walk_q_children` = `walk_items_wf.shapes.q_children`, restated over `g.dfsForest vo eo`, `WalkItemsWF.lean`) |
| §4.2b walk invariant `Inv' D` (`EntryInv' D above`: connected + attached at `Term'`; closed items 2-attached): closure lemmas (`GraphLemmas.lean`: `AttachedIn`, `twoAttached_iff`), primitives `Inv'.alloc`/`modifyVs`/`pushVert`/`pushEdge`/`mergeTop`/`retarget`/`pop`/`finishTop`, `Shape`/`Step` infrastructure, per-block lemmas `Step.closeEars`/`mergeLate`/`closeVert'`/`finishRest` under `CloseEarsOk`/`MergeLateOk`/`CloseVertOk`/`FinishRestOk` (`MergeTopOk`/`RetargetOk.disj`: entries below edge-disjoint from the touched ones) | `GraphLemmas.lean`, `WalkSpec.lean` | proved (the depth-indexed `Inv D` versions `mergeTstackTops_sound`/`finishTstackTop_complete` are kept; `Inv D` itself is false mid-walk) |
| §4.2b `finishEdge_inv` (type-1 ± vertex entry, type-2 three loops, back edge) under `FinishOk`; `finishEdge_back_inv` corollary | `WalkSpec.lean` | proved for `Inv' D`; `FinishOk ← FinishGuards`/`EarFinish` (`ear_*` sorries, `ear_finishP_back` derived); `WalkInv.walkTree_inv'` (`Inv' d ∧ Shape`) proved modulo them — the `Inv d` form was **false** (`ear_lower` false, `Inv D` fails for every `D` under a type-2 chain with a sibling subtree; §4.2b correction), `ear_lower'` (`Inv' (d+1) → Inv' d` after a tree edge) derived from `EarFinish.lower` in its place; `walk_nodes_partition` proved in `WalkPlace.lean` under `ForestOK` + coverage |
| Invariant W, Lemma 4.3 (`earOut_one_entry`) | `EarShape.lean` | restated as `finishEdge_one_entry` (at the `finishEdge` site, from `EarAt`, type-1 edges; §4.2b) and proved (standard axioms); the type-2 form needs `lowval ≤ topDepth` of the child's entries, not in the contract |
| between-edges invariant `EarCtx` of `walkOuts` (induction hypothesis for `walkTree_ear`; §4.2b) | `EarCtx.lean`, `EarCtxAt.lean`, `checks/EarCheck.lean` (`ctxCheck`) | stated, dump-checked 0..3000 both modes (0 violations); `EarAt` at a back-edge site derived from it (`earAt_back_of_ctx`, standard axioms); `EarAt` at a tree-edge site derived from it + the child's end-of-outs `EarCtx` + frame (`TreeSite`, `earAt_tree_of_ctx`) up to the 13 named admissions `earAt_tree_{loop1_side,loop1_touch,loop1,bottom,loops,late,late_fo,close,bd_bridge,bd_comp,bd_term,bd_side,lower}`; `walkTree_ear` (with the edge-completeness hypothesis `hcomp`, forced by the kernel-checked `walkTree_ear_false`/`walkTree_book_false` in `EarFalse.lean`) proved as an assembly (`cTree`/`cOuts`/`cOut`, `WalkInv.lean`, §4.3) over these, the proved `ctx_init_root`/`ctx_init_child`/`ctx_step_back_boundary` (standard axioms; `earCtx_selfLoop`), the proved `ctx_step_back_ret` (standard axioms: `earCtx_pushBack` for the `Q` push, `ctx_step_back_ret_P` = `mergeP_single` + `earCtx_mergeP` for the P-merge), the proved `ctx_step_tree_bridge` (standard axioms: `earCtx_bridge`, `below_kept`), the proved `ctx_step_tree_comp` (`earCtx_comp`, standard axioms; `tree_comp_shape` = the child's end-of-outs stack shape at a component edge, dump-checked `compEndCheck`, proved from the second between-edges invariant `CtxShape`, a conclusion of the same induction), the proved `ctx_step_tree_ret` (`wp_finishEdge_tree`, `earCtx_ret`, standard axioms, over the named admission `tree_ret_shape` = `RetTop` after loops 1–3, dump-checked `retCheck`; `ctx_step_tree_ret_P`, the P-merge with the vertex entry, is proved; `ctx_step_tree` is the case split on the class), and the named admission `tree_ret_shape` (dump-checked `retCheck`); the item frame of a child walk (`KeptRel` + freshness, formerly the admission `walkTree_items_kept`; `OutFrame` per out, `outFrame_*`; dump-checked `kept_*`/`freshCheck`/`siteKeptCheck`) is a conclusion of the induction and `walkTree_{qroot,vroot,below}_kept` are corollaries; `endsOut_wf`/`endsOuts_wf`/`ends_of_wf_boundary` (`bdOut_wf`) proved |
| Lemma 4.4 (`ascend_frame_one_entry`: a finished frame's vertex owns one entry) | `EarSpec.lean` | **false** for chain frames (cycle `0..5` + chord `5-1`, frame `(4,4)`: five entries); removed, the collapse holds only at the ear's top (= `earOut_one_entry`) |
| boundary branch of `finishEdge` keeps `Inv' D ∧ Shape` (`finishBoundary_inv`, via `BStep`, under `BoundaryOk`: popped entries exist, Q/V items are roots not on any span, popped blocks' terminals touched by no entry below) | `WalkInv.lean` | proved (`BoundaryOk ← ear_boundary`, derived from `dest_edges`/`bd_noVert`/`bd_bridge`/`bd_comp`/`bd_term`; `gone`/`gone₂` stated under `isTree`); the former `Step` form is false (the Q item goes under `vertItem curV`), as is `VertBook`'s `hasVert = false → ch (vertItem v) = []` (bridge `1-2` before back edge `1→0`) — replaced by connectivity + `TwoAttached v v` of the vertex item |
| corrected attachment set `TEntry.Term'`, `EntryInv'`, `Stack`, `Inv'` (`Inv'.of_inv`, `Inv'.mono`, `Stack_iff`, `Inv'.setSv`); empirical check `checks/InvCheck.lean` | `WalkSpec.lean`, `WalkInv.lean` | def + proved; `Step`/`BStep`/`walkTree_inv'` stated through it |
| 4.1–4.3 ear content at a `finishEdge` (`WalkState.EarFinish`/`EarAt`/`EarBottom`: `sub ++ base` split, pairwise edge/span disjointness, span ownership `subEdges`, distinct path vertices, loop-1 range side `stackDir[d]` and touching, P-merge target `p_entry`, ear-bottom anchor `bottom`, post-loop shape `loops`, V/Q item freshness, block-boundary separation) | `EarInv.lean` | **stated**, every field 0 violations on 3000 random multigraphs (`checks/EarCheck.lean`); `FinishBook.ear` carries it, `finishOk_of_guards`/`ear_*` take it as hypothesis; false first versions (`base_top`, `loop1_bot`, `touch_top`, `vert.vStart`) dropped with counterexamples (§4.2b); carried by the walk: `walkTree_ear` (`EarTree`, the ear/vert fields of `BookTree`) is a proved assembly (`cTree`/`cOuts`/`cOut`, §4.3) over `EarCtx` modulo the named admissions listed in §4.2b; `walkTree_book` proved from it, `bTree`; `EarShape.finishGuards` was false, see §4.2b, replaced by `finishGuards_of_ear`); derived so far: `ear_finishP_back` (the back-edge P merge, from `p_entry`/`touch_bot`/`path`/`disj`/`span_disj`/`q_free`/`q_root`, via `maybeUnwrapNxt_edges`: the unwrap of a one-sided single-item `nxt` keeps every entry's edge set), `ear_loop1` (`Loop1BodyOk` at every loop-1 iterate, `EarLoop1.lean`: invariant `L1Inv` = reached split `L1Reach` of the range + top piece `L1Piece` (bottom `l1Bot`, edges `l1Edges`, one-sided, fresh) + untouched rest `L1Keep` + frame `L1Frame`; `l1_step` discharges `MergeTopOk`/`UnwrapOk`/`CloseTwoOk` from `EarFinish.loop1` (`Loop1Spec`: `L1Merge`/`L1Close`/`L1Unwrap` per reached split) + `loop1_side`/`disj`/`span_disj`/`q_root`/`q_free`, `L1Ctx.ofEar` takes the range as the longest `≥ d` prefix of `sub`, proper by `bottom`), `ear_finishP_tree` (no vertex entry: after loops 1–2 `nxt` is the `(y, lowval)` piece `py` (`loops`, type 1 ⇒ `mid = []`), whose bottom is not `curV` (`sub_bot`), so `condP` is false); loop-2/close contract stated and checked (`late`: `EarLate` at `feS₁`, `close`: `EarClose`/`FoldSpec` at `feS₂`, path fields `sv_d`/`sv_child`/`path_child`/`dir_d`; 0 violations); `ear_mergeLate` derived from `late` (`EarLoop2.lean`: `mergeLateOk_of_late`, standard axioms), `ear_closeVert` derived from `close` (`closeVertOk_of_close`: loop 3 via `loop_run_iter`/`mergeTopOk_iter_fold`, type-1 unwrap via `MergeOk.congr_entry`; standard axioms), `ear_finishP_vert` derived from `close` + `p_entry`/`base_root` (`finishPOk_type1_of_close`; standard axioms); `ear_tail_tree` derived from `close` (`ear_condP_tree`, `vert_touch`/`vert_disj`/`c_edge`; standard axioms); `ear_lower'` derived from `lower`/`sv_child`, `ear_boundary` from `base_touch`/`bd_*`/`sv_d`/`sv_child` (`dest_edges` replaced, §4.2b) (standard axioms); every `ear_*` of `finishEdge_step` is now derived from `EarFinish`; `Frontier d o origTstack s` (standalone frontier/ownership fact for the range layer) derived by `finishEdge_frontier`; `bd_side` (boundary orientation, 0 violations) added for `BoundaryOK`; `finishEdge_sides`/`walkTree_sides` (`EarSides.lean`) derive `FinishSides`/`SidesTree` (standard axioms) |
| 4.5 maximality: `RCloseShape` ⇒ no skeleton pair separates (`RCloseShape.not_sepPair`), R skeleton 3-connected (`RCloseShape.threeConnected`) | `RMax.lean`, `Proofs/RMax.lean` | proved; `RStep.rCloseShape'`/`RStep.threeConnected'` (`Proofs/RClose.lean`) give it for Loop 1's R step from `Inv' D` + `stackVerts[d+1..D] = cur.vStart` + `RStep` + `RContent` (`Inv d` is contradictory there) |
| 4.5 walk side: `EntryR`/`RTop`/`RBranch`/`RInvAt`; `RBranch.rStep`, `RBranch.rContent` (all five content fields), `RBranch.threeConnected` (from `Inv' (d+1)`); `EntryR.congr`/`RInvAt.congr` + bookkeeping frames; `Items.RSkel3`, `RBranch.rSkel3` | `RInv.lean`, `Proofs/RInv.lean`, `Proofs/RInvFrame.lean`, `Proofs/RItems.lean`, `checks/RBranchCounter.lean` | proved; admitted: `finishEdge_rInvAt`, `walkTree_rInvAt`, `loop1_rBranch_mid_ctx`/`loop1_rBranch_fields_ctx`/`loop1_rTop_ctx` (`Proofs/RLoop1.lean`, the ear-context form on the pipeline; `RBranch.cur_c` refuted by `checks/RBranchCounter.lean` (kernel-checked, R-5) and restated as `mid`; shape + `cur_top` + `cur_ne` proved as `loop1_r_shape_ctx` from `L1Piece`, `interior` as `loop1_r_interior_ctx` from `L1Close.bottom`/`mid`, `loop1_rBranch_content_ctx` proved from these; the ear-context-free `loop1_rBranch_content` kept for `loop1_rBranch`'s fixed statement, shape part `loop1_r_shape` proved); `items_r_three_connected` and `spqrTree_r_three_connected` proved modulo the named admissions of §4.5 (`walk_rSkelInv` threading, R-4); `closeVert_type1_rSkel3` proved from `vClose_rSkel3` (`Proofs/RVert.lean`, R-5), admitted: `closeVert_type1_rCloseShape` (the `RCloseShape` of the type-1 vertex close); `dfsForest_rooted` restated (re-rooted `ofForest`, `checks/RRootedCounter.lean` refutes the first form) and proved, standard axioms |
| 4.5 R skeleton persistence: `Pieces.contract_congr`, `Items.RSkel3.congr`, `.modify_of_not_below`, `.push_nil` | `Proofs/RItems.lean` | proved (standard axioms); walk-level ownership of later writes remains open |
| 4.5 HT-to-cut transport: `Graph.ThreeConnected.relabel` | `Proofs/ThreeConnected.lean` | proved (standard axioms), assuming a block and its active vertex list; item-level hypotheses remain open |
| 4.5 Item/contract edge correspondence: `Pieces.ofItems_addParent_edges`, `Items.rSkeleton_perm_contract` | `Proofs/RItems.lean` | proved (standard axioms), under non-V-child edge coverage; deriving coverage for completed R items remains open |
| 4.5 HT 3-connectivity implies 2-connectivity with ≥ 3 edges: `Graph.ThreeConnected.twoConnected` | `Proofs/ThreeConnected.lean` | proved (standard axioms); no block or simplicity assumption |
| 4.5 Item HT-to-cut bridge: `Items.RSkel3.rThreeConnected` | `Proofs/RItems.lean`, `Proofs/ThreeConnected.lean` | proved (standard axioms), including active vertices and either cap orientation; takes item WF and non-V-child edge coverage |
| 4.5 Coverage in a 2-connected block: `WF.q_root_covers_block`, `WF.vertex_child_leaf_of_block`, `WF.nonV_child_cover_of_block`, `RSkel3.rThreeConnected_of_block` | `Proofs/RItems.lean` | proved (standard axioms); requires Q children of V items to have nonempty children, a walk placement fact not yet exported |
| Item S/R shape minimum counts | `ItemSpec.lean`, `checks/SkeletonShapeCheck.lean` | corrected: ≥ 1 V child for S, ≥ 5 non-V children for R; kernel-checked sharp triangle/K4 outputs (standard axioms) |
| R correctness input domain | `Correctness.lean`, `Proofs/RItems.lean`, `checks/RInvalidOrderCheck.lean` | corrected to `g.WF` + `OrderOK` for both orders; K4 with invalid edge order `[6]` kernel-checks failure of the former public target (standard axioms) |
| R coverage and laminarity edge domain | `RInv.lean`, `RClose.lean`, `RMax.lean`, `Proofs/RunSaturation.lean`, `checks/REdgeDomainCheck.lean` | restricted containment and coverage to `e < g.ne`; kernel-checked K4 failure of unrestricted coverage and success of bounded coverage; affected transport proofs audited (standard axioms) |
| 4.5 Run saturation and interval-to-run laminarity | `Proofs/RunSaturation.lean` | `Saturated` stated; eight conditional lemmas proved, standard axioms only; walk preservation and marker alignment remain open |
| 4.5 Depth-bounded settling `RInvTop`/`RInvFront`; `finishEdge_rInvTop` (replaces `finishEdge_rInvAt`): back-edge branch proved (`finishEdge_back_rInvTop`), tree-edge branch proved up to the top entry's `EntryR` in the first-edge case (`finishEdge_tree_rInvTop_of_top`; `hasVert = true` top exempt by `finishEdge_tree_vert_topStart`; `hasVert = false` proved by `finishEdge_tree_top_settled_first` modulo the admitted `feS₂_top_entryR`: `EntryR` of the `feS₂` top when it tops out at `d` and does not start at `curV`) | `RInv.lean`, `Proofs/RInvFrame.lean`, `Proofs/RInvBack.lean`, `Proofs/RInvTree.lean`, `checks/RFinishEdgeCounter.lean`, `checks/RFinishEdgeCheck.lean` | two kernel-checked counterexamples to the parent-only exemption (standard axioms); contract B checked on 6414 sites, 0 failures; frame lemmas, the back-edge branch (`RInvTop.pushEdge/finishP/unwrapNxt_exempt/mergeTop_exempt/finishTop_exempt`) and the tree-edge branch below the top (`RInvG.closeEars`, `RInvH.mergeLate/closeVert'/finishRest`) proved, standard axioms; `FinishRShape` checked on 6414 sites |
| 4.5 Parent-base frame `finishEdge_rInvG_base` (any `finishEdge` with `n₀ ≤ origTstack` preserves `RInvG dfs p dp n₀`; positional `RInvG.*_pos` primitives, `loop3_iter_len`) | `Proofs/RInvBase.lean`, `checks/RFinishEdgeCheck.lean` (`base` lines) | proved (standard axioms); contract dump-checked on every ancestor frame of every `finishEdge` site, seeds 0..300 × both modes + 6000 random, 0 failures |
| 4.5 Child-return induction `rrTree`/`rrOuts`/`rrOut` (`RWalk`: ancestor frames `RInvG` + `RInvTop` of the walked vertex; `BotKeep`/`finishEdge_bot_keep`), `walkTree_rReturn`, `walkTreeRReturnSpec`; side facts `RSideTree`/`RSideOuts`/`RSideOut` (`AncChain`, `EntryR` stability under `stackVerts.set!`, parent entries top out `≤ d-1`, non-root out-edges return, `VertFree`, `FinishRShape`) | `Proofs/RInvWalk.lean`, `checks/RFinishEdgeCheck.lean` (`anc`/`stab`/`entry`/`ret`/`vertown` lines) | induction and `walkTree_rReturn` proved from `finishEdge_rInvTop` + `finishEdge_rInvG_base` (non-standard axiom only via the admitted `finishEdge_tree_top_settled`); side facts dump-checked at every entry/site (seeds 0..400 × both modes + 6000 random, 0 failures); `walkTree_rSide` (R-4 statement correction: takes `BookTree` and the parent's ancestor chain, both root-threaded and supplied by `rkRootOut`; standalone form `walkTree_rSide_spec` kept for `WalkTreeRReturnSpec` only) proved in R-5 by `rsTree`/`rsOuts`/`rsOut` (`Proofs/RSide.lean`) modulo the three site admissions `rSide_entry_site`/`rSide_vertFree_site`/`rSide_finish_content_site` (R content; `FinishRShape.pend`/`ear` proved) |
| 4.5 Child-return settling diagnostic and provisional contract | `Proofs/RInvFrame.lean`, `checks/RInvReturnCheck.lean` | `RReturn`/`WalkTreeRReturnSpec` stated without an admission; legacy conclusion still refuted; fixed base/settled-entry clauses kernel-checked (standard axioms); seeds 0..400 × both modes pass shape/disjointness, with content frames checked on the 91 block inputs (terminal-value frame observed for `topDepth ≤ d` only — seed 390's buried single-`Q` entries, R-5); preservation proof remains open |
| 4.5 Schedule frontier: `Frontier`, `FrontiersTree` | `Proofs/RInvFrame.lean`, `EarFrontier.lean` | stated and threaded into `finishEdge_rInvTop`; ear export proved: `finishEdge_frontier` (from `FinishBook.ear` + `Inv'`/`Shape`), `walkTree_frontiers` (`FrontiersTree` under the `walkTree_inv'` hypotheses), standard axioms; R interval/saturation preservation remains open |
| 4.6 walk-time range invariant `WalkState.RangesInv σ n D` (`Inv' D` + `processed`/`ordered`/`convex`/`closed`; `TEntry.piece`, `Items.BelowNoV_congr`/`_modify_of_not_below`): `RangesInv.alloc`/`pushVert`/`pushEdge`/`mergeTop` (local adjacency `hadj`)/`finishTop` | `RangesInv.lean`, `checks/RangesInvCheck.lean` | proved (standard axioms); 0 violations at every `finishEdge` (seeds 0..400 × tern + tiny graphs); `finishEdge_rangesInv` (`RangesStep.lean`, under `FinishAdj`) and `walkTree_rangesInv` (`RangesTree.lean`, under `GuardsTree`/`BookTree`/`RgTree`) proved; `finishBoundary_rangesInv` proved (`RangesTree.lean`); `finishEdge_ownedD` proved (`RangesOwned.lean`, over `OwnSite`: `OwnedD.of_split`/`of_boundary`, `Dec`, `OFrame`); `closeCtx_bd_vert`/`bd_node`/`p_site`/`v_site`/`l1_site` (`RangesCloseTree.lean`) proved from `CloseBase` and `CloseContent` (`RangesCloseContent.lean`; `checkContent` 0 violations 0..3000 × both modes), the named hypothesis `closeBase_content : CloseBase → CloseContent` left to the backbone induction; `finishEdge_vertCover` proved (`VSpan`, `RangesOwned.lean`); `walk_rootsCover` proved as the assembly over the former (`RangesCoverTree.lean`: `cvTree`/`cvOuts`/`cvOut`, `OwnedD.root`, `forest_rootsCover`, `walk_rootsCover'`) and `walk_closeInv` over the latter (`RangesCloseTree.lean`: `finishEdge_closeInv`, `ccTree`/`ccOuts`/`ccOut`, `CsTree` with `DfsSite` proved by `dsTree`, `CloseCtx.of_exports`, `forest_closeInv`, `walk_closeInv'`) and six checked per-site statements, `closeEars_closeAt`/`finishBack_closeAt` (leaf Q, `leafQ_closeInv`) and `finishBoundary_closeAt` (root Q, `CloseAt.rootQ`; new `CloseCtx` fields `dest_lt`/`bd_loop`/`bd_vert`/`bd_node`, checker `checkCtx`) and `finishP_closeAt` (P record `CloseAt.pNode'` via `PSite.closeInv`; new `CloseCtx` field `p_site : … → PSite`, checker `checkP`) and `closeVertTail_closeAt` (S/R record `CloseAt.node` via `VSite.closeInv`; new `CloseCtx` field `v_site : … → ∃ t, VSite`, checker `checkV`) and `loop1Body_closeAt` (same record via `CloseInv.l1S₁`/`maybeUnwrap`/`mergeTop` + `VSite.closeInv`; new `CloseCtx` field `l1_site`, checker `checkL1`) proved (`RangesCloseSites.lean`, 0 violations at every block boundary, seeds 0..400 × tern + tiny graphs); `ranges_of_rangesInv` (`RangesFinal.lean`) derives `convex` + node `att_vs` from `RangesInv`, the rest is `Items.CloseFacts`; saturation not a field (attachment-count forms false, §4.6); `Owned`/`finishP_ownership` (`RangesOwned.lean`, with `unwrap_edges`/`finishTop_edges`/`closeVert'_edges`): P-site coverage of `walk_rootsCover` proved (standard axioms), checker `checkOwned` (`own_*`) 0 failures 0..400 × both modes; `OwnedD`/`OwnSite`/`finishEdge_ownedD` (proved, standard axioms; `checkOwned` at `entry`/`pre`/`post`/`P`/`end`, 0 failures) and `VertCover`/`finishEdge_vertCover` (proved, standard axioms, via `VSpan` preserved by every `finishEdge` primitive; `checkVCover` `own_vcover` at `post`/`end`, 0 failures) the two steps of the proved `walk_rootsCover` induction (`RangesCoverTree.lean`); public site exports `cbTree`/`forest_closeBase`/`walk_closeBase` (`CloseBase` at every `finishEdge` pre-state) and `rangesInv_l1Iter`/`rangesInv_feS₂` (`RgStep` into loop-1 iterates / `feS₂` from `CloseBase`, standard axioms; checker `checkRI`) in `RangesSites.lean`; canonicity `WalkState.CanonInv`/`CloseCanon`/`FinishCanon`/`finishEdge_canon` (`RangesCanon.lean`, standard axioms; checker `checkCanon` `canon_p`/`canon_v`/`canon_l1` 0 violations 0..3000 × both modes), named hypotheses `walk_canonInv`/`walk_ternarize` for the backbone |
| 5 relabel: `Items.WF → Items.ROriented → WF` | `relabelTree_wf` (`Correctness.lean`, = `RelabelAll.wf_tree`) | proved (`RelabelWF.lean`) |
| 5 relabel: `relabelTree_represents : Items.WF → Items.RThreeConnected → Represents` (`Correctness.lean`, = `relabelTree_represents'`), `relabelTree_represents_of_r` (output-level R clause, used by `spqrTree_represents`); per field `RelabelOK.q_endpoints/twin_glue/nv_orig_inj/separation/interior/canonical/r_three_connected` | `RelabelRep.lean` | proved (every `RelabelOK.*` field is standard-axioms only); needs the `Items.WF` clauses `Endpoints.q_root`, `Shapes.o_parent`, `Shapes.s_order` (§5; checked by `check_repok`); `Items.RThreeConnected` is the item-level R statement (§4.5, `items_r_three_connected`), transported not proved |
| 5 relabel, per-node layout: `Layout.Shape`/`Layout.Local` for F, V, Q-loop/O, Q/I, P, S, R (`shape_*`, `local_*`), exact rows (`runF_row`, `runLoop_row`, `runQI_row`, `runP_row`, `runS_row`, `run_entries`) | `LayoutShape.lean` | proved (standard axioms); `r_skeleton_nodup` discharges the R `Nodup` hypothesis from `r_shape` |
| 5 relabel, structural part: `relabelTree_own : Items.WF → Items.ROriented → Bijections ∧ Ownership ∧ Twins` (also `relabelTree_bijections`, `relabelTree_twins` from `WF` alone); `layoutNode_edges` | `RelabelOwn.lean` | proved (`RelabelSpec.lean`, sorry); the item-level facts it needs are clauses of `Items.WF` (`nv_nodup`, `q_leaf_of_node`, `q_children`'s `v < nv`, `r_shape`'s child endpoints) |
| 5 relabel, CSR bounds: `relabelTree_adj : Items.WF → Items.ROriented → (∀ n, adjBounds[2 nvSt n] = 2 neSt n) ∧ adjBounds[2 |nodeVerts|] = |adjDat|` (the statement of `relabel_adj_spec`); `layout_local` (`Layout.Local` for every item's `nodeLayout`) | `RelabelAdj.lean` | proved; `relabel_adj_spec` itself stays admitted in `RelabelSpec.lean` only because that file cannot import its proof |
| 5 relabel, `WF` assembly: `RelabelAll.wf_tree`, `preorder` (`child_idx`/`subtree_end` chain, `subtree_props`, `parent_eq_iff`), `only_root_F`, `shape` (`skeleton_eq` + `LayoutShape.shape_*`), `adj_bounds_mono`, `adj_dest`, `adj_incident'` (`global_bound`/`global_row`: global CSR rows = `Layout.Local` rows; `foreign_ne`, `row_filter`) | `RelabelWF.lean` | proved (standard axioms) |
| 2/7 `Items.ROriented` of the walk output: `rOriented_of_stNumbered`, `walk_items_rOriented'` | `StOriented.lean` | proved from `walk_st` (+ `walk_items_wf`), under `g.WF`/`OrderOK` like `walk_st`; used by `spqrTree_wf'` (`Correctness.lean`); the hypothesis-free `spqrTree_wf` the planar layer uses is a named admission (§7.6) |
| 2 walk→relabel interface `walk_items_wf : g.WF → OrderOK g.nv vo → OrderOK g.ne eo → Items.WF g (g.walk tern (g.dfsForest vo eo)).items` (`WalkItemsWF.lean`, above `WalkCover`/`WalkTyping`), and `spqrTree_eq` (`WalkWF.lean`, below the st layer) | `walk_items_wf = wf_of_ranges walk_tree.toTree walk_typing.toTypingFacts walk_ranges` (empty graph: `wf_initialItems`); admitted: `walk_rangesInv`/`walk_closeFacts` (§4.6; `walk_ranges` is derived), `walk_canonical` (`tern = false`, used by `spqrTree_canonical`) is the proved assembly over the named hypotheses `walk_canonInv`/`walk_ternarize` (step `finishEdge_canon` proved, `RangesCanon.lean`, §4.6); `spqrTree_eq` proved |
| §4.6 final-state ranges `Items.Ranges g items σ` (piece-edge convexity in the DFS edge postorder + attachment facts), `wf_of_ranges : Tree → TypingFacts → Ranges → canonical → Items.WF` (`endpoints_of_ranges`, `shapes_of_ranges`); checker `check_ranges` (`CheckRanges.lean`) | `Ranges.lean`, `RangesWF.lean` | def / proved (standard axioms); every field 0 violations on seeds 0..400 + tiny graphs; corrected `s_shape` (`1 ≤`), `r_shape` (`5 ≤`), `q_children` (Q leaf child allowed) in `ItemSpec.lean` with counterexamples |
| 7 st-order spec `StOrder`, `Items.StNumbered`, split `spqrTree_st = relabel_st ∘ walk_st` | `StSpec.lean`, `StWalk.lean` | def / proved split |
| 5 relabel per-node interface `RelabelNode`/`RelabelLayout`/`RelabelIdx`, `Items.nvList`/`ordered`/`edgeChildren`/`PosOK`/`hasCap`/`nEdges`/`ROriented` | `RelabelSpec.lean` | def; `relabel_node_spec` proved in `RelabelMain.lean` (`Ghost.relabel_node_spec_proved`) |
| 5 relabel ghost (`order` field) + refinement `relabelTree_eq : relabelTree g items = ofRelabelState (relabelRun g items)` | `RelabelGhost.lean` | proved (`sim_relabel` at `maxHeartbeats 4000000`) |
| 5 relabel state invariants `Consistent`, `Agree` (append-only prefixes), `NodeS`, `LowB`, transport `NodeS.mono`/`LowB.agree`, `RelabelLayout.congr` | `RelabelInv.lean`, `RelabelMono.lean` | proved |
| 5 relabel per-call contract `CallPre`/`CallPost`, node-body `Step*` lemmas, `entry_of_steps`, children loop `LoopInv.pre`/`loop_init`/`LoopInv.step`/`LoopInv.fin` | `RelabelLoop.lean`, `RelabelWp.lean` | proved |
| 5 relabel `relabel_spec` (wp induction over `relabel`), `relabel_node_spec_proved` | `RelabelProof.lean`, `RelabelMain.lean` | proved, standard axioms |
| 5 relabel adjacency CSR `relabel_adj_spec` (needs `ROriented`) | `RelabelSpec.lean` | proved (`RelabelAdj.relabelTree_adj`, re-exported as `relabel_adj_spec`) |
| 7 relabel-side: `vchildren_nv_increasing`, `orderedChildren_sorted`, `edgeChildren_dominance`, `layoutNode_r_bracket` | `StSpec.lean`, `StLayout.lean` | proved |
| 7 relabel-side: `relabel_st : Items.StNumbered → Items.WF → StOrder` (`RelabelOK.stOrder`: `st_S/st_P/st_R`, `dom_node`, `adj_node` via `adjRow_eq_layout` + `LayoutShape` rows + `layoutNode_R_bracket`) | `RelabelSt.lean` | proved (through `relabelOK_of_wf` and `relabelTree_adj`) |
| 7 walk-side: `WalkState.StInv`, data lemmas `pushTstack_onSide`, `merge_onSide`, `fold_onSide` | `StWalk.lean` | def / proved |
| 7 walk-side: `finishTstackTop_stItem`; ear lowvals `first_ret_lowval`, `chain_stackDir_step` | `StWalk.lean`, `StEar.lean` | proved |
| 7 walk-side: `StInv.onSide` field, `chain_stackDir_const` (corrected statement, see 7.4) | `StWalk.lean` | def / proved |
| 7 walk-side: `StInv.hole` (`StHole`/`HoleClosed`), `stInv_topClosable`, `stInv_finishTstackTop_stItem` (close site, see 7.4) | `StWalk.lean` | def / proved |
| 7 walk-side: `finishEdge_stInv` (under `FinishGuards`/`EarsOnSide`), `walk_stInv`; `walk_st`, `spqrTree_st` from `walk_stInv`; route changed to the `StRef.lean` reference order (§7.6) | `StWalk.lean`, `StRef.lean` | sorry / proved (`walkTree_stackDir_below` in `StFrame.lean` proved); reference + `VsOriented` + `StBlock.St` tested seeds 0..1000; `refBlocks_st` proved (`StRefEt.lean`, standard axioms); `stItem_of_refOrder` proved (`StRestrict.lean`, standard axioms); `walk_st'`, `walk_vsOriented` sorry |

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
* `relabel_st` **[proved]** (`RelabelSt.lean`): the glue is the per-node
  interface `RelabelSpec.lean` (`RelabelNode.nv_range/ne_range/node_verts`, `RelabelLayout`'s
  `edge_nvs/adj_bounds/adj_dat/pos_ok`), taken through `RelabelOK` (`RelabelRep.lean`).
  - `StOrder.st`: `skeleton_S` (path `(nvSt, nvEn-1) :: (nvSt+k, nvSt+k+1)`), `skeleton_P`
    (`nEdges` copies of `(nvSt, nvSt+1)`), `skeleton_R` (cap plus `edgeChildren`, every child edge
    `nvSt ≤ a < c < nvEn` by `pos_ok` and `Items.StItem`'s orientation — `r_facts`); then every
    interior node-vert has an edge in and an edge out: S/P by construction, R because
    `Items.StList` gives it on `vertList` and `pos_eq` (`pos x = nvSt + idxOf x`, from `pos_ok` on
    a `Nodup` list) moves it to positions.
  - `StOrder.dom`: F/V/Q/I/O have at most one own edge (`dom_small`, `nEdges_Q = 1` from
    `q_children`); S and P by the explicit skeletons; R by `pairwise_dominance_of_sorted_sum` on
    the `loc`-sorted `edgeChildren` (`Items.ordered_pairwise_loc`).
  - `StOrder.adj`: `nodeOfNv_range` locates `nv` in a node's `nvRange` (`nvBounds_cover`, the
    `RelabelIdx` CSR facts); `adjRow_eq_layout` rewrites `adjRow r` as `LayoutR.row` of the node's
    `layoutNode` given `adjBounds[2 nvSt] = 2 neSt` (`relabelTree_adj`), using `adj_bounds`,
    `adj_dat` and `Layout.Local`'s bound monotonicity (`rowBound_ge/le` keep the row inside the
    node's `adjDat` segment); the per-type rows are `LayoutShape.runQI_row/runP_row/runS_row` (F/V:
    no node-verts with `1 < nVerts`), and R is `layoutNode_R_bracket` with `edgeChildren`'s
    dominance, `Nodup` (`r_skeleton_nodup`) and `neEn = neSt + |E| + 1`.
  `Items.StNumbered` is needed only for `Items.ROriented` (`rOriented`) and the R skeleton facts;
  everything else is `Items.WF`.

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
`stItem_of_refOrder` — an item listed in the reference order is in s-t order — without the tstack:
`refBlocks_st` by induction over `refTree` (every piece spliced at depth `l` joins the open path on
side `dirs[l]`, so every vertex other than a block's endpoints has a neighbour on each side) and
the restriction `stItem_of_block` (`StRestrict.lean`; the item's virtual edges are the ends of the
sub-ears, so the restriction to an item keeps the neighbours).
The one walk-side ingredient of the latter is the orientation of the children's `vs`: `makeVs`
uses the same `stackDir[d]` as the splice side, which is why `stItem_of_refOrder` is stated for
the walk's items rather than for arbitrary `WF` items with `ch` in reference order (for those the
`vs` could be flipped). The §7.4 `StInv` route (`finishEdge_stInv`, `walk_stInv`,
`walk_st_of_stInv`) has been removed; `StInv` and its proved close-site lemmas remain as documentation of the invariant.

**Blocks and orientation.** The reference now also records, per completed block (`StBlock`), the
block-boundary tree edge `(v, w)` it hangs from (`root = none` for the one-vertex block of a DFS
root): the block's vertex sequence `StBlock.seq` is the cut vertex `v` followed by the vertices of
its items, and its edge list `StBlock.edges` is `(v, w)` followed by its edge items. The
orientation statement is `VsOriented g items blocks`: every S / P / R item `i` lies in one block
`b` (all its leaves are items of `b`) and its `vs` and the `vs` of its non-V children are
`Oriented` along `b.seq` — `(some a, some c)` with `a` strictly before `c` (`Precedes`). The
differential test `check_stref` checks, besides `ch = restrictCh …`, `VsOriented` for the walk's
items and `StBlock.St` (Even–Tarjan: `Items.StList (b.seq g) (b.edges g)`) for every reference
block; all three hold on seeds 0..1000. Accordingly `stItem_of_refOrder` takes `VsOriented` as a
hypothesis and splits into `refBlocks_st` (every reference block is st-numbered; the `refTree`
induction) and the restriction argument, and the walk-side admissions are `walk_st'` and
`walk_vsOriented` (both to be proved by the same simulation).

**Even–Tarjan on the reference (`StRefEt.lean`).** `refBlocks_st` is proved by the mutual
induction `refTree_inv`/`refOuts_inv` over `refTree`/`refOuts`. The invariant
`VInv g anc dirs v vs rets hasVert L` on the item list `L = stNest pieces` of the open ear at
vertex `v` (ancestor path `anc`, splice sides `dirs`, return depths `rets` so far): `L` consists
of vertex items (`v` and `vs`, no repetition) and edge items of `g`; every edge item joins a
vertex item of `L` to another one or to `anc[l]` for a return depth `l` (`edge_mem`); every
vertex item other than `v` has a neighbour on each side (`Nb`: an edge item of `L` to a vertex
item before / after it in `L` — `Before`, `Side` — or a back edge to `anc[l]` counted on side
`dirs[l]`); every vertex other than `v` lies on side `dirs[l]` of `v` for one of the return
depths (`root_side`); and `v` has a neighbour on the side `dirs[l]` of its lowest return
(`root_nb`). `VInv.step` inserts one returning piece on side `dirs[l]` of its return depth
(`VInv.step_back` for a back edge; `VInv.step_tree` for a returning subtree, whose own invariant
has `v` as the path vertex at depth `anc.length` — `Nb.popV`/`NbV.toNb` transport it when `v`
becomes an explicit endpoint); boundary edges only add return depths `≥ anc.length`
(`VInv.extend`); `block_st` turns the invariant of the subtree below a boundary tree edge `(v, w)`
into `StBlock.St` (sequence `v :: vertices of the pieces`, the returns to `v` being the
neighbours of the first vertex on the path side), and `root_block_st` is the one-vertex block of
a DFS root. The DFS facts come from `DfsTree.WF` (`dfsForest_wf`), `DfsForestSpec.joins` /
`verts_nodup` (`dfsForestSpec_of_dfsForest`) and the vertex bound (`dfsForest_spanning'`), so
`refBlocks_st` is stated under `g.WF`, `OrderOK g.nv vo`, `OrderOK g.ne eo` — a correction of the
hypothesis-free admission: without `g.WF` it is false, e.g. for `g = ⟨1, #[(0, 1)]⟩` the DFS
reaches the non-vertex `1` and `StBlock.seq`/`edges` read `vertItem 1 = 1 + g.nv + 0` as the edge
item `0`, so the block's edge `(0, 1)` has an endpoint outside `seq` (traced by hand; `#eval`
aborts on out-of-range endpoints). Accordingly `stItem_of_refOrder` takes the conclusion of
`refBlocks_st` as the hypothesis `hbl`, and `walk_st`/`spqrTree_st` take `g.WF` and the order
hypotheses like `dfsForest_spanning`; so do `walk_items_rOriented'` (`StOriented.lean`) and `spqrTree_wf'`
(`Correctness.lean`). The hypothesis-free `spqrTree_wf` that the planar layer consumes
(`planarTree_shape`, `planarRelabel_rot_spec`, `neRotAdj_segment'`) is kept as a named admission: it is
`spqrTree_wf'` under `g.WF`/`OrderOK`, and whether its `ROriented` part holds for malformed graphs is
open. `VsOriented` also records that a block edge (an edge item of `b`, or an edge on `b.root`'s
endpoints) below `i` is a leaf of `i` — i.e. not in a block hanging off a V child of `i` — which the
restriction argument needs to know that a block neighbour of a V child of `i` comes from one of
`i`'s own children; that the V children of `i` lie strictly between `i`'s endpoints in `b.seq`; and
that a block edge below a non-V child `c` has both endpoints between `c`'s endpoints (the items'
vertex sets are nested sub-ears) — all checked by `check_stref`, seeds 0..1000. `refBlocks_root_none`
(the `root = none` block is a DFS root's `[vertItem r]`) and `refBlocks_root_edge` (a block's boundary
edge is a tree edge, hence `PairEq` to some `g.edges[e]!`; from `DfsForestSpec.joins`) are the two
structural facts about blocks. `#print axioms refBlocks_st`: propext, Classical.choice, Quot.sound.

**The simulation (`walk_st'`, `walk_vsOriented`; open).** The relation is `readStack (tstack above
the base) = stNest pieces` with `pieces` the reference's pieces for the open ear, *up to expanding
the items closed inside the ear back to the children they were closed with* (`expandItem`): a
closed item occurs on the stack as the single item of a one-sided entry and the reference lists
its leaves instead. Per-primitive steps, all proved: `readStack_pushTstack` (`pushTstack v d i`
pushes the piece `⟨stackDir[d], [i]⟩`, matching the reference's `dirs[d]`), `readStack_mergeTstackTops`
(no change), `readStack_fold` (the vertex close folds the top entry to one piece on side `!edgeDir`,
matching the reference's `StPiece.mk lowDir (stNest sub)`), `readStack_finishTstackTop` (a close
replaces the entry's items by the new item: `expandItem item … = …`) and `readStack_reopen` /
`readStack_modifyNxt_reopen` (`maybeUnwrapNxt`'s reopen of a closed S / P node under the top is the
inverse expansion).

*The relation (`StSim.lean`, checked by `check_stsim` / `compare_stsim.sh`, seeds 0..1000, 0 violations
at every boundary).* Expansion is the fuel-free inductive `ExpandsList items xs L` (`L` = the
concatenated leaves of `xs`; a V / Q item is its own leaf, any other item is replaced by its children),
`Expands items x L := ExpandsList items [x] L`; `StRead items new ps := ExpandsList items (readL new) (stNestL ps) ∧ ExpandsList items (readR new)
(stNestR ps)` (per side, so that readings of adjacent stack segments concatenate); `DirsOf s d := [stackDir[0], …, stackDir[d-1]]` is the reference's `dirs`. `StSim g prev
fs t base s` is the postcondition of `walkTree t d` (`d = fs.length`, `base` the tstack before, `prev`
the finished trees of the forest): `s.tstack = new ++ base` with `StRead s.items new (refTree g t d
(DirsOf s d)).1`, and `StItems g s (refBlocks g (prev ++ [truncTree fs t]))`; `StSimOuts g prev fs v
done hasVert base s` is the same at the start of an out-edge of `v` with `refOuts g v d dirs done
false` (pieces, and its `hasVert` equals the walk's). `StItems`: the items on the stack are roots
(`¬ IsParent p x`), `readStack tstack` is `Nodup` (what `readStack_reopen` needs — the reopened
children are not elsewhere on the stack), the stack's subtrees and all children are in range (`bounded`, `chLt`, so `allocItem` does not
touch them), child lists have no repeats (`chNodup`, what reopening needs), every S / P / R item strictly below a
leaf-type (V / Q) item is finished (`finished`, sim-4; what the root pop and the boundary completion consume), and every S / P / R item is *live* (`Below x i` for some
stack item `x`) or *finished*: `InBlock g items b i` for a block `b`, i.e. `b.items = A ++ L ++ B`
with `Expands items i L` and `VsOrientedAt g items b i L` (the body of `VsOriented` for one block
with the leaf list `L` in place of `leaves items items.size i`, `vsOriented_iff`; the final step
identifies the two via `Tree.acyclic`/`ch_lt`). The blocks are those of the **truncated tree** `truncTree fs t`: along the open
path (`PathFrame` = vertex, finished out-edges `done`, current tree edge `o`) every vertex keeps
`done ++ [o]` with the next frame as `o`'s child, and `t` at the bottom. Its reference order is the
walk's *prediction* of the block order: it already accounts for the pending `closeVert` folds of
the open frames (the only reordering primitive — a fold moves the whole sub-ear reading `L ++ R` to
one side of the base, so an item closed at depth `d` on side `edgeDir` with `vs = (vStart, curV)`
is on the wrong side of `V curV` before the fold and on the right side after it; a prediction from
the current tstack alone is false there). As the walk proceeds the truncated reference only grows
outwards (`refOuts_append`/`stNest_append`: later pieces wrap the earlier ones; a fold or a
`hasVert` vertex piece wraps the sub-ear as a unit), so `Precedes`/`Oriented` facts of live items
persist, and a finished block of the truncated tree is a final block (its subtree is complete).
The checker mirrors `walkTree`/`walkOuts`/`walkOut` step by step calling the real `finishEdge`,
checks `StSimOuts` at the start of every out-edge, `StSim` of the child plus `StRead (pre ++ new)`
before every `finishEdge` (`sub` above `origTstack` reads as the child's `refTree` pieces) and
`StSim` at the end of every `walkTree`, and finally that the mirrored items equal `g.walk`'s.

`StUnwrap.lean` composes `StClose.lean` into the tail shared by `loop1Body` and `finishP`:
`StSim.unwrapMergeClose` (`maybeUnwrapNxt ty; mergeTstackTops; finishTstackTop`) keeps `StRead` above
`base` and `StItems`, in both the allocation path (`ExpandsList.push`/`InBlock.push`, the new node is
out of every stack subtree by `bounded`) and the reopen path (`readStack_reopen` + `expandItem_self_iff`;
the reopened node's children are fresh on the stack because they have a parent and stack items are
roots); `L1Unwrap` is taken as the hypothesis `hU`, only its `getSide t.spans dir = [h]` clause is used.
`StLoop1.lean` threads it through loop 1: `l1St_loop` takes the ear side's `L1Ctx` (the facts of
`EarFinish` loop 1 consumes) and `L1Inv` (its per-iteration stack shape: the piece `c` at depth `d`,
one-sided, over the untouched rest) as hypotheses and keeps `L1StInv` — the reading above the block's
base and `StItems` are unchanged up to expansion, `stackDir` at depths `≤ d` is unchanged (the type-2
S-merge writes `stackDir[t.topDepth]` with `t.topDepth > d`), so the reference's `dirs` are fixed.
`maybeUnwrapNxt` reads the side `stackDir[t.topDepth]`, which need not be the ear's side: if it is
not, the side is empty and `L1Unwrap` (hence the head is no S / P item) forces the allocation path.
`StEars.lean` composes it into `closeEars_st`: the whole of `feS₁` (the tree edge's `vs`, the
pushed `Q` piece, loop 1) under exactly `l1_init`'s hypotheses.
`StMerge.lean` (`L1StInv.mergeLate`), `StVert.lean` (`closeVert_st`: the `closeVert` unwrap /
reopen / merge / fold of the type-1 close, carrying the lower entries `pre` reading as `qs`;
`finishP_st`: the P merge when a `(curV, lowval)` single-piece entry is on the stack, its side
conditions only under `isType1 = true`; `finishTail_st`: the vertex push, needing `vertItem curV`
root, typed `V`, off the stack only under `hasVert = false`) and `StTree.lean` compose the whole of
`finishTree` (`finishTree_st`, tree edge with `lowval < d`, `hasVert` either way: result pieces
`qs ++ [⟨!stackDir[d], stNest (ps ++ [Q e])⟩]` resp. `qs ++ ps ++ [Q e, V curV]`, the lower
stack `B` untouched, `stackDir[k]` for `k ≤ d` kept) and `finishBack` (`finishBack_st`, back edge
with `lowval < d`: `qs ++ [⟨stackDir[lv], [Q e]⟩]` then the optional `V curV`). Both take the
ear facts as `EarFinish`/`FinishOk` (`loops`, `bottom`, `p_entry`, `q_root`/`q_free`, `rest_*`) and
one extra named hypothesis `hvf`: for `hasVert = false`, `vertItem curV` is a root and off the stack
*after* `finishP` (it is so before `finishEdge` by `EarFinish.v_root`/`vert_free`; carrying it
through `closeEars`/`mergeLate`/`closeVert`/`finishP` is a `Fresh`-style frame lemma per primitive,
not yet written).
**Status (sim-3).** (2)–(4) below are done: `stOut_step`/`stOuts_*`/`stTree_node`/`stWalk`/`stForest`
(`StInduct.lean`) and the final step `walk_sim`/`walk_inBlock`/`restrictCh_eq` (`StFinal.lean`) make
`walk_st'`/`walk_vsOriented` proved modulo the named admission `finishBoundary_st` (item (1),
`StBoundary.lean`);
`hvf` is discharged inside `stOut_step` from `Full`. **Sim-4:** `finishEdge_vStart` is proved
(`StVStart.lean`, a `vStart`-provenance invariant `VsIn` through every `finishEdge` primitive; its
statement gained the disjunct "`vStart = o.dest` for a returning tree edge", which the original
statement lacked — `closeEars` pushes the child's ear with `vStart = o.dest` and nothing below it
need absorb it — and `stRet_finish` consumes the extra case via `o.dest ∈ child.verts`), and
`refOrder_nodup` is proved (`StNodup.lean`: `refTree_items`/`refOuts_items` — the pieces and blocks
of a subtree hold only its vertex/edge items, each once — then `ForestOK`), and `walk_i_parent`
is proved (`StIParent.lean`: the invariant `IOk` — vertex items are `V`, edge items are `Q`, child
ids are in range, no stack span item is an `I`, and every `I` child has a `Q` parent — through every
walk primitive; the only `I` allocation is the bridge branch of `finishBoundary`, whose parent is
`edgeItem g o.e`; `walk_i_parent` now takes `g.WF`/`OrderOK`, needed for `dfsForest_bounded`), and
`finishRet_frame_st` is proved (`StRetFrame.lean`: a protected set `U` of item ids with the frame
`Fr` — `U`-items keep type/children, no non-`U` item ever gets a `U` child — and `Out` — every
tstack entry above the untouched base `B` is `U`-free — carried through every primitive of a
returning `finishEdge`; instantiated with `U` = "below a span item of `B`" for the base frame and
with `U = {V curV}` for the `hvf` facts; the hypotheses are `EarFinish`, `Full` (`cnt ≤ 1` gives
unique parents and a parentless, off-stack `rootItem`), `StItems` and `curV < g.nv`, so
`stRet_finish` now also takes `Full`, which the tree case obtains from `walk_full_aux`), and
`rootPop_st` is proved: `StItems` gained the clause `finished` — every S/P/R item strictly below a
leaf-type (V / Q) item is `InBlock` of some block (block roots hang under the `Q` of the boundary
edge / the `V` of the root entry) — added to `check_stsim` first (`stItemsB`, `m7`; seeds 0..1000,
0 violations) and then carried through every `StItems` step (`congr`/`perm`/`modify_root`/
`pushEntry`/`close`, the fresh-item and reopen sites of `StUnwrap`/`StVert`); since the root entry
holds only vertex items (`RootOK`), every S/P/R item below it is already finished, so the root pop
only has to frame the old blocks past the root's children change (`InBlock.modify_root`) and the
new block `⟨none, stNest ps⟩` is appended with nothing to show. For `finishBoundary_st` the same
clause means: the block completed under `Q e` must contain every S/P/R item below it, which is
exactly the `InBlock` obligation for the popped entries. The
boundary statement is exactly what the induction consumes and is dump-checked by `check_stsim`
(seeds 0..1000). `finishBoundary_st` itself is now proved (`StBdPop.lean` + `StBoundary.lean`)
modulo the single named admission `finishBoundary_stLive` (sim-4; see the end of this block): the popped-ear shape is already
an `EarFinish` field (`back_nil` for a back edge, `bd_bridge` — `sub = [t]` — for `lowval = d + 1`,
`bd_comp` — `sub = [t₁, t₂]` — otherwise; no new ear field was needed); the item-array change of
every branch is abstracted as `BdPop` (frame off `Q e`/`V curV`, fresh items are not S/P/R and
childless, `ch (Q e)` ⊆ `readStack sub` ∪ fresh, `ch (V curV)` gains `Q e`) and `StItems.bdPop`
carries `StItems` through it; `StRead.complete` turns `StRead s.items sub ps` into `InBlock` for
every S/P/R item below a span item of `sub` — an item below a member of an expanded list either
owns a segment of the expansion (`Items.Below.expands_cases`/`ExpandsList.segment_of_mem`) or lies
below a leaf-type item and is already finished (`StItems.finished`). What `StRead` does not carry
is the orientation core: the live hypothesis `hlive` of `finishBoundary_st` (sim-5; `StLive` of the
popped segment `sub` against the closing block, supplied by the named admission
`finishBoundary_stLive`; `StLive.vsOrientedAt` turns it into the
`VsOrientedAt g s.items ⟨some (curV, o.dest), stNest ps⟩ i L` clauses for every S/P/R item `i`
below a span item of `sub` with leaves `L`).

*The open-block invariant (sim-4, `StLive`, `StBdPop.lean`).* `StLive g items new b := ∀ x ∈ readStack
new, ∀ i, Below x i → i is S/P/R → InBlock g items b i`: every live S/P/R item under a stack segment
is already `InBlock` of the block the segment belongs to — for the walk, the block of the
**truncated** reference `refBlocks g (prev ++ [truncTree fs t])` (the design `fe86c46` moved away
from in favour of the completed-only `simBlocks`, which is why `StItems.closed` only says "live or
finished"). `finishBoundary_stLive` is `StLive s.items sub ⟨some (curV, o.dest), stNest ps⟩` at a
boundary tree edge, and is exactly what `finishBoundary_st` needs (`StLive.vsOrientedAt`, by
`ExpandsList.unique`). Checked first (`check_stsim`, `truncItemsB`, at every out-edge, before every
`finishEdge` and at the end of every `walkTree`, against `refBlocks g (prev ++ [truncTree fs ·])`
with `·` = `.node v done` / `.node v (done ++ [o])` / `t`; seeds 0..1000, 0 violations): m8 — every
S/P/R item is `InBlock` of some truncated block (= `StLive` for every open segment, given the
blocks' disjointness); m10 — every live non-V item with two endpoints is `Oriented` in a truncated
block (what the new item's clause 5 needs at a close); m11 — for every live S/P/R item, the block
edges below it lie weakly between its endpoints (clause 6 at a close; note the restriction to the
block's own edges — `Below` also reaches the hanging sub-blocks through `V → Q`). What a
preservation proof of `StLive` through `finishTstackTop` still lacks is a *per-entry* fact for
clauses 3–4 of the new item (`Oriented (setSides dir (stackVerts[t.topDepth]) t.vStart)`, V span
items strictly between); the obvious candidates are FALSE and are kept in `CheckStSim.lean`
(`truncEntryChecks := false`) with the counterexamples: `stackVerts[t.topDepth]` is stale for the
buried entries of finished sibling subtrees (seeds 554, 6: an entry with `topDepth = d` whose
vertex at `d` has since changed), a two-sided merged entry is oriented only by the side
`stackDir[t.topDepth]` picks (seed 0), and the type-2 fold retargets `vStart := curV` so the V span
items of a to-be-folded entry need not lie between its *current* terminals (seed 9: entry `(3, 0)`
holding `V 4` with the final order `0 … 3, 4, 2`). The right per-entry fact must therefore name the
entry's top vertex explicitly (not via `stackVerts`) and be stated for the one-sided entries that
a close actually consumes (loop 1: `loop1_side`; the P entry: `p_entry`; the fold: against
`(stackVerts[c.topDepth], curV)` after the fold, from `FoldSpec`/`EarClose`), i.e. it is an ear-side
invariant in the sense of §7.6's `EarFinish` fields rather than a `StRead` clause.

*The site-level terminal facts (sim-5, `check_stsim` `truncCtxB`/`truncSiteB`, seeds 0..1000, 0
violations).* Stated only for the entries whose terminals are current, against the truncated
blocks: (m16/m17, at every out-edge, before every `finishEdge`, at the end of every `walkTree` of
`v` at depth `d`) every entry `(v, l)` with `l < d` — `EarCtx.above` — has `(stackVerts[l], v)`
`Oriented` by `stackDir[l]` (i.e. `setSides stackDir[l] stackVerts[l] v`, the `Q` convention) and its
V span items strictly between; (m13/m14, before `finishEdge o` of `v` with `lowval < d`) the ear's
terminals `(stackVerts[lowval], v)` are `Oriented` by `stackDir[lowval]` and every V span item of
the child's leftover `sub` (other than `vertItem v`) is strictly between — the fold / tail close;
(m15, same site) for every loop-1 entry `t'` of `sub` (the `topDepth ≥ d` prefix) with
`t'.topDepth = d`, `(v, t'.vStart)` is `Oriented` by `stackDir[d]` and the V span items of the
prefix down to `t'` are strictly between — the loop-1 closes, whose bottom is the consumed entry's
`vStart`. Kept as counterexamples next to the statements: m9/m12 (any entry, via
`stackVerts[t.topDepth]`) fail on seeds 554/6 (stale depth of a buried entry), 0 (two-sided entry)
and 9 (fold retarget); m18 ("loop 1's range has `topDepth ≤ d + 1`") fails on seeds 154, 195, 442 —
the buried `V y`/`(y, l)` pair of a chain bottom sits in loop 1's range with its deeper `topDepth`
(`EarBottom.vy_top`), so the loop-1 statement must be per consumed entry, not per depth.

**StLive obligations (sim-5 handoff for the backbone induction).** The separate `stWalk`
re-threading was stopped (regroup: one backbone induction carries `EarCtx`, `RangesInv`/`CloseInv`
and the St relations together); the components below are what it should reuse.

*Definition to carry (`StBdPop.lean`).* `HangingUnder items x i := ∃ y, Below x y ∧ y ≠ i ∧ y is
V/Q ∧ Below y i` and
`StLive g items new b := ∀ x ∈ readStack new, ∀ i, Below x i → i is S/P/R → ¬ HangingUnder x i →
InBlock g items b i`. The `HangingUnder` filter is necessary: an S/P/R item below a V/Q item of a
live span item belongs to an already completed block hanging off that vertex/edge, not to the open
block (the unfiltered form is checker-false; those items are `InBlock` of a completed block by
`StItems.finished`). The block of an open segment is `openBlock` (`StOpenBlock.lean`):
`openBlock g fs dirs ps` is the block of the truncated reference containing the bottom pieces `ps`
under the path frames `fs` with directions `dirs` — `openBlock_snoc_bd` (the innermost frame is a
boundary edge: `⟨some (f.v, f.o.dest), stNest ps⟩`), `openBlock_snoc_ret` (a returning frame:
the parent's block with `(refOuts g f.v k dirs f.done false).1 ++ retPieces … ps` as bottom
pieces), `openBlock_ctx` (under `TreeFrames fs`: `openBlock g fs dirs ps = ⟨r, A ++ stNest ps ++ B⟩`
with `r`, `A`, `B` independent of `ps` — the context into which the open segment's items are
spliced; `TreeFrames.of_ok` from the per-frame bookkeeping `FrameOK`). The per-segment invariant to
carry alongside `StSim`/`StItems` (`StPre` in `StInduct.lean`, with `segs` the lower segments
`(new_k, qs_k)` of the frames `fs` and `new` the current segment, `tstack = new ++ segsStack segs`):

* current segment at `v`, depth `d = fs.length`, out-edges `done` done, flag `hasVert`:
  `StLive g s.items new (openBlock g fs (DirsOf s d) ((refOuts g v d (DirsOf s d) done false).1 ++
  (if hasVert then [] else [⟨true, [vertItem v]⟩])))`;
* lower segment `k < d`: `StLive g s.items new_k (openBlock g (fs.take k) (DirsOf s k) qs_k)`;
* at the end of `walkTree t d` (what the parent's returning / boundary step consumes):
  `StLive g s'.items new' (openBlock g fs (DirsOf s' d) (refTree g t d (DirsOf s' d)).1)`
  (`refTree_node` relates it to the first form with `hasVert := (refOuts …).2.2`).

Not yet added to `check_stsim` in this per-segment form: m8 checks the union form (every S/P/R
item in *some* truncated block); the pairing segment ↔ `openBlock` is the only new content and
should be added (`openBlock` is computable; sites: `chkOuts` with `above s base`, end of
`chkTree`) before the backbone relies on it.

*Transport lemmas that exist (`StOpen.lean`, all proved, axioms ⊆ {propext, Classical.choice,
Quot.sound}).* `StLive.nil`, `StLive.of_subset` (any segment whose span items are among the
old ones — covers `mergeTstackTops`, pops, and permutations of entries), `StLive.append`,
`StLive.cons`, `StLive.of_leafTypes` (a segment whose span items are all V/Q, e.g. the pushed
`⟨v, d, _, setSides dir [vertItem v] []⟩` / `[edgeItem e]` entries), `StLive.frame` (under
`ItemsFrame items items' new`: type, children and `vs` of everything below the segment's span items
unchanged — the lower segments across any sub-walk / `finishEdge`), `StLive.insert` (fresh items
`P`, `Q` spliced around `M = stNest ps` inside the block `⟨r, A ++ M ++ B⟩`: given `StRead items
new ps`, `P ++ Q` disjoint from the block and from `vertItem r.1`, and nothing below the segment
in `P ++ Q`, the segment stays live in `⟨r, A ++ P ++ M ++ Q ++ B⟩` — the block growing by the
pieces of later edges; built on `InBlock.insert`, `Precedes.insert`, `Oriented.insert`,
`segment_insert`), `StLive.vsOrientedAt` (what `finishBoundary_st` consumes),
`HangingUnder.iff_of_frame`, `Items.Below.{of_ch_frame, to_of_ch_frame}`. `StRetFrame.lean`'s
`Fr` now also keeps `vs` of the protected items, so `finishRet_frame_st` and the frame conclusion
of `finishBoundary_st` deliver `ItemsFrame` for every base segment (lower segments stay live by
`StLive.frame` across every `finishEdge`). `StKeepSv.lean`: `KeepsSv`/`KeepsSvBelow` per primitive
and `walk_keepsSvBelow` → `walkTree_stackVerts_below` (`walkTree t d` fixes `stackVerts[k]` for
`k < d`), `walkOut_stackVerts_le` (`walkOut v d` fixes `k ≤ d`): the path identification
`stackVerts[k]! = fs[k].v` that `openBlock`'s roots need.

*Per-primitive status.* Covered by the lemmas above (composition not written): `walkOutPre`
(`setStackDir` — `DirsOf s d` unchanged; type-1 vertex push: `StLive.cons`/`of_leafTypes` +
`StLive.insert` with `openBlock_ctx` for the new piece `stNest_snoc`), the pushes of `finishBack`
/ `finishBoundary` (`StRead.pushEntry` + `StLive.insert`), all merges (`StLive.of_subset`), the
boundary step for the lower segments (`StLive.frame` with `finishBoundary_st`'s frame conclusion),
the returning step for the lower segments (`finishRet_frame_st`), the root pop (nothing live
remains). **Missing** — the closes, i.e. every allocation of an S/P/R item whose children are the
span items of a consumed one-sided entry: `finishTstackTop` in loop 1 (`loop1Type`), the
`closeVert` fold, `finishP` and `maybeUnwrapNxt`'s reopen (children appended to an existing S/P
item). For the new item `i` the six `InBlock` clauses split as: segment of leaves — from `StRead`
(`StRead.complete` / `Items.Below.expands_cases`); non-V children oriented and block edges below
a non-V child between its endpoints — from the invariant itself (the child is live, hence
`InBlock`, clauses 1 and 6 of the child); `vs i` oriented and the V children strictly between —
exactly the site-level facts m13/m14 (fold and tail), m15 (loop-1 entries with `topDepth = d`),
m16/m17 (`EarCtx.above` entries `(v, l)`, `l < d`, at every site) of `truncCtxB`/`truncSiteB`
(seeds 0..1000). These are statements about the walk state against the truncated reference's
order, so they cannot be `EarFinish` fields; they are additional invariant clauses for the
backbone (`StLiveCtx`: m16/m17 at every site; m13/m15 at the `finishEdge` site). Their expected
proofs are reference-side (`VInv.root_side`/`VInv.step_tree` of `StRefEt.lean`: a piece returning
to depth `l` lies on side `dirs[l]` of the path vertex at `l`, nested inside the pieces of lower
return depth), given the identification of the walk's terminals with the path: `EarFinish.sv_d`
(`stackVerts[d] = curV`), `sv_child`, `path_child`, `path`, `dir_d` (`stackDir[d] =
!stackDir[lowval]`), `EarCtx.sv`/`sv_d`/`path`, and of the consumed entries' sides:
`EarFinish.loop1_side` (loop-1 entries are one-sided on `!stackDir[d]`), `p_entry` (the P entry),
`close`/`bottom`/`late` (`FoldSpec`/`EarClose` for the fold's entry, whose `vStart` becomes
`curV`), `EarCtx.above` (which entries are `(v, l)`). The false per-entry forms m9/m12/m18 and
their seeds are listed above; do not restate them with `stackVerts[t.topDepth]`.

*What `finishBoundary_stLive` needs.* With the invariant carried, it is an instance of the
end-of-`walkTree` form for the child: `openBlock g (fs ++ [⟨curV, done, o⟩]) … ps =
⟨some (curV, o.dest), stNest ps⟩` by `openBlock_snoc_bd` (from `hge : d ≤ lowval`), the popped
segment `sub` being the child's `new'` (`EarFinish.tstack`, `bd_bridge`/`bd_comp`), and the
block root `(curV, o.dest)` from `sv_d`/`sv_child`; no further ear fact is needed at the
boundary site itself — the whole obligation lives in the closes listed above. Without the
invariant it is unprovable from `EarFinish`/`EarCtx` alone (they do not mention the reference
order).

The original plan, for reference:
What remains: (1) `finishBoundary` — the only place a block completes: popping the child's ear
must give `InBlock ⟨some (v, o.dest), stNest ps⟩` for every S/P/R
item below the popped entries, i.e. the segment clause from `StRead sub ps` (done) *and* the
`VsOrientedAt` clauses (`Oriented (b.seq g) (vs i)`, V children between the endpoints, …) — these
need an invariant on the `vs` of every stack span item relative to the entry's `vStart`/`topDepth`
and the st-order, which `StRead` does not carry (the vs-orientation core of `walk_vsOriented`;
now `finishBoundary_stLive`, see above);
(2) the `walkOut` step for a returning edge: `walkOutPre_st`, the child's `StSim` (its `sub`
reads as `refTree child (d+1) (DirsOf s₁ (d+1))`, `DirsOf s₁ (d+1) = DirsOf s d ++ [sd]`), then
`finishEdge_st`, matched against `refOut` (`DirsOf_getD`) and `refOuts (done ++ [o])` (a snoc
lemma for `refOuts`, `truncTree (fs ++ [f]) t = truncTree fs (.node f.v (f.done ++ [.tree … t]))`);
`StItems` must also be transported from the blocks of `truncTree fs (.node v done)` to those of
`truncTree fs (.node v (done ++ [o]))` — only the truncated tree's root block changes, by pieces
prepended on the L side / appended on the R side of `stNest`, so segments stay segments and the
monotone `VsOrientedAt` clauses survive; the block-edge clause needs the new `Q`/`V` items to be
below no finished item (they are roots); (3) the `walkOuts`/`walkTree` induction (`StSimOuts` →
`StSim`, consuming `walkOut_eq`, `finishOk_of_guards`, `finishEdge_step`, `EarAt` from
`FinishBook`); (4) `walk_st'` from the final
`InBlock`: `Expands i L` with `L` a segment of the `Nodup` `refOrder` and the children's expansions
nonempty gives `restrictCh` = `ch` once `Items.leaves items items.size` is shown to agree with
`Expands` (height < size in a tree).
`walk_vsOriented` is `InBlock.oriented` at the end of the walk (everything is finished once the
stack is empty) transported from the truncated to the final blocks (monotone clauses, plus the
block-edge clause from disjointness of the blocks' items).

**The restriction argument (`StRestrict.lean`).** `stItem_of_block` proves `Items.StItem items i`
for one block `b` with `b.St` from `Items.WF`, `ch i = restrictCh … order i` where `order` contains
`b.items` as a segment, and the `VsOriented` clauses for `i`. With `vs i = (some s, some t)` the
item's vertex list is `s :: xs ++ [t]` (`xs` the V children's vertices); the order of two V children
in `ch i` is their order in `b.items` (`before_ch_of_before_order`: `Before b.items → Before order →
Before (filterMap find?) → Before (collapseRuns …)`, since `find?` picks the child itself for a V
child), so `Precedes b.seq` on the item's vertices implies `Precedes (vertList i)` (the between
clause handles `s`/`t`). For an interior `x ∈ xs`, `b.St` gives a block edge `(x, z)` with `z` before
(resp. after) `x` in `b.seq`; it is below `i` (`interior`), hence a leaf of `i` (the block-edge clause),
hence below a non-V child `c` with `vs c = (some u, some v)`; `separation` forces `x ∈ {u, v}`, the span
clause excludes `x = u` (resp. `x = v`), and `(u, v) ∈ virtualEdges i` is the required neighbour. The
third `StItem` clause is the children's orientation transported by the same order lemma.
`stItem_of_refOrder` is `stItem_of_block` on the block supplied by `VsOriented`, with `refBlocks_st`
and `refBlocks_root_edge`. `#print axioms stItem_of_block`, `stItem_of_refOrder`: propext,
Classical.choice, Quot.sound.

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
| `relabel_st` | `RelabelSt.lean` | proved (7.3) |
| `TEntry.wrap`, `OneSided`, `OnSide`, `nest` | `StWalk.lean` | def |
| `pushTstack_onSide`, `merge_onSide`, `fold_onSide`, `getSide_setSides` | `StWalk.lean` | proved |
| `WalkState.StSides`, `StEntry`, `StInv`, `TopClosable`, `entryVertList`, `entryEdges` | `StWalk.lean` | def |
| `DfsOut.lowval_eq_lmin`, `first_ret_lowval`, `chain_stackDir_step` | `StEar.lean` | proved |
| `walkTree_stackDir_below` (frame: `stackDir` below `d` unchanged by `walkTree _ d`) | `StFrame.lean` | proved |
| `StInv.onSide`, `chain_stackDir_const` (`ear_uniform_side`, semantic half; corrected statement) | `StWalk.lean` | def / proved |
| `finishTstackTop_items`, `finishTstackTop_stItem` | `StWalk.lean` | proved |
| `StInv.hole` (`StHole`, `EntryReach`, `HoleClosed`), `idxOf_lt_idxOf_iff`, `stList_of_sorted`, `stInv_topClosable`, `stInv_finishTstackTop_stItem` | `StWalk.lean` | def / proved (replaces `finishEdge_topClosable`, see 7.4) |
| `EarsOnSide`, `finishEdge_stInv`, `walk_stInv`, `walk_st_of_stInv` | `StWalk.lean` | removed (the `StInv` preservation route; superseded by `walk_st'`, nothing depended on them) |
| `walk_st`, `spqrTree_st` (under `g.WF`, `OrderOK g.nv vo`, `OrderOK g.ne eo`, like `dfsForest_spanning`) | `StFinal.lean` (moved from `StWalk.lean`: `walk_st'` is proved below `StWalk`) | proved (from `walk_st'`, `stItem_of_refOrder`, `refBlocks_st`, `walk_items_wf`, `relabel_st`); `walk_st_of_stInv` is the same from the alternative `walk_stInv` route |
| `refTree`/`refOrder`, `restrictCh`, `check_stref` differential test (§7.6) | `StRef.lean`, `CheckStRef.lean` | def / tested seeds 0..300 (0 mismatches) |
| `walk_st'` (`ch i = restrictCh … (refOrder …) i` for S/P/R items; statement change: under `g.WF`, `OrderOK g.nv vo`, `OrderOK g.ne eo`, like `walk_st`) | `StFinal.lean` | proved modulo the named admissions `finishBoundary_st`, `rootPop_st` (`StBoundary.lean`) and the ear-side inputs (`walkTree_ear`, `dfsForest_ends`) |
| `StBlock`, `StBlock.seq`/`edges`/`St`, `Precedes`, `Oriented`, `VsOriented`; `check_stref` checks `VsOriented` and `StBlock.St` too (§7.6) | `StRef.lean`, `CheckStRef.lean` | def / tested seeds 0..1000 (0 mismatches) |
| `walk_vsOriented` (`VsOriented` for the walk's items and `refBlocks`: block membership of the leaves, block edges below `i` are leaves of `i`, `vs` orientation of `i` and its non-V children, V children between `i`'s endpoints, block edges below a non-V child between its endpoints; same statement change as `walk_st'`) | `StFinal.lean` | proved modulo the same admissions as `walk_st'` |
| `refBlocks_root_none` (a block without boundary edge is a DFS root's one-vertex block) | `StRefEt.lean` | proved (`refTree_roots`; axioms propext, Quot.sound) |
| `refBlocks_root_edge` (a block's boundary edge is `PairEq` to some `g.edges[e]!`, `e < g.ne`; under `g.WF`, `OrderOK`) | `StRefEt.lean` | proved (`refTree_root_edges`/`refOuts_root_edges`/`refOut_root_edges` from `TreeJoins`; standard axioms) |
| `refBlocks_st` (every reference block is `StBlock.St`; corrected statement: under `g.WF`, `OrderOK g.nv vo`, `OrderOK g.ne eo`) | `StRefEt.lean` | proved (`refTree_inv`/`refOuts_inv`; axioms propext, Classical.choice, Quot.sound) |
| `Before`, `Side`, `Nb`, `NbV`, `VInv` (`nil`, `single`, `extend`, `step`, `step_back`, `step_tree`), `block_st`, `root_block_st`, `refOut_boundary_back/tree`, `refOut_ret_back/tree`, `refOuts_nil/cons/zero`, `refTree_node` | `StRefEt.lean` | def / proved |
| `Before.*` (`or_of_mem_mem`, `filter`, `filter_of`, `map`, `cons`, `append_left/right`, `filterMap`, …), `collapseRuns_sublist`, `before_collapseRuns`, `Precedes.of_before/before/ne/asymm/trans`, `Items.below_of_mem_leaves`, `leaves_of_type`, `leaves_succ`, `Tree.child_below_child`, `Tree.child_class`, `Tree.v_child_eq`, `Tree.type_ne_F`, `find?_leaves_vertItem`, `before_ch_of_before_order`, `edgeOf_eq_some` | `StRestrict.lean` | proved |
| `stItem_of_block` (the restriction argument for one block, §7.6), `stItem_of_refOrder` (on the block from `VsOriented`, with `refBlocks_st` as `hbl` and `refBlocks_root_edge` as `hre`) | `StRestrict.lean` | proved (axioms propext, Classical.choice, Quot.sound) |
| `Expands`/`ExpandsList`, `DirsOf`, `StRead`, `VsOrientedAt` (`vsOriented_iff`), `InBlock`, `StItems`, `PathFrame`/`truncTree`, `StSim`, `StSimOuts` (the simulation relation, §7.6); `check_stsim` / `compare_stsim.sh` | `StSim.lean`, `CheckStSim.lean` | def / tested seeds 0..1000 at every `finishEdge` boundary (0 violations) |
| `ExpandsList.{append, append_inv, cons_iff, det, congr, modify_of_not_below, modify_root, push, expandItem_self_iff, close}` (expansion algebra: framing under `modify`/`push`, `expandItem` of a node is invisible, closing `i` with children `ch` keeps the reading) | `StSim.lean` | proved |
| `StClose.lean`: `Items.Below.{modify_of_not_below, of_modify, push, of_push, lt_of_chLt}`, `VsOrientedAt.congr`, `InBlock.{congr, modify_root, push}` (framing of the finished-item facts under `modify` of a root / `push`), `readStack_{cons_perm, close, close_perm}`, `StRead.finishTstackTop` (the reading above `base` survives the close of a one-sided top), `StItems.close` (the item facts survive it: `item` is a root not on the stack, the stack below holds roots; `closed` is assumed with `item` exempt) | `StClose.lean` | proved |
| `StUnwrap.lean`: `readL_append`/`readR_append`, `readStack_cons_perm'`, `mem_readStack_cons`, `mem_expandItem`, `nodup_expandItem`, `mergeTops_cons_cons`, `TEntry.mergeInto_side_nil`, `StSim.allocMergeClose` (`allocItem; mergeTstackTops; finishTstackTop size`), `StSim.reopenMergeClose` (`modifyNxt` reopen of a single node `h`; merge; `finishTstackTop h`), `StSim.unwrapMergeClose` (`maybeUnwrapNxt ty; mergeTstackTops; finishTstackTop` on `c :: t :: _`, both one-sided on `stackDir[t.topDepth]` = the direction at the merged depth, `ty ∈ {S, P, R}`, `L1Unwrap s ty t` as the named hypothesis `hU`): the reading above `base` and `StItems` survive, the new top is `{mergeInto c t with spans := setSides dir [item] []}` with `type item = ty` — the common tail of `loop1Body` and `finishP` | `StUnwrap.lean` | proved (axioms propext, Classical.choice, Quot.sound) |
| `StLoop1.lean`: `L1StInv` (the st-side loop-1 invariant: the entries above a fixed suffix `B` read as fixed pieces `ps`, `StItems` for fixed `blocks`, `stackDir[k]` for `k ≤ d` kept), `L1Unwrap.transport` (`L1Unwrap` read before `finishEdge` holds at an iteration state whose kept entries are `L1Keep` and whose top items are roots or on the original stack), `loop_run_iter` (`loop fuel cond body` is `iter body k` with the condition true before each iteration), `l1St_step` (one `loop1Body` iteration under `L1Ctx`/`L1Inv` of `EarLoop1.lean`: the type-2 S-merge sets `stackDir` above `d` only, then `StSim.unwrapMergeClose`), `l1St_iter`, `l1St_loop` (the whole loop 1, with `l1_iter` supplying `L1Inv` at every iteration) | `StLoop1.lean` | proved (axioms propext, Classical.choice, Quot.sound) |
| `StEars.lean`: `ExpandsList.leaves`, `mem_readStack_exists`, `mem_readStack_push`, `nodup_readStack_push`, `StRead.pushEntry` (pushing the one-item entry `⟨dir, [q]⟩`, `q` a leaf, appends the piece), `StItems.modify_root` (a `vs`-only change of a root non-S/P/R item keeps `StItems`), `StItems.pushEntry` (pushing a root childless item absent from the stack keeps `StItems`), `closeEars_st` (`feS₁ d o s` for a tree edge with `lowval < d`: `feS₀` sets `vs` of `Q e`, `pushEdgeTstack` pushes `⟨stackDir[d], [Q e]⟩`, then `l1St_loop` via `L1Ctx.ofEar`/`l1_init`; result `L1StInv` with pieces `ps ++ [⟨stackDir[d], [Q e]⟩]` above `base`) | `StEars.lean` | proved (axioms propext, Classical.choice, Quot.sound) |
| `StMerge.lean`: `mergeTopsN`, `iter_mergeTstackTops`, `mergeTopsN_length`, `mergeTopsN_above` (`k` merges above a suffix `B` keep `readL`/`readR` of the part above `B`), `L1StInv.mergeTopsN`, `L1StInv.mergeLoop` (any `loop _ cond mergeTstackTops` that ends with an entry above `B` keeps `L1StInv`: merges only shorten the stack, so all of them happened above `B`), `L1StInv.mergeLate` (loop 2) | `StMerge.lean` | proved (axioms propext, Classical.choice, Quot.sound) |
| `StVert.lean`: `ExpandsList.{unique, split}`, `StRead.{append, split, fold}`, `StItems.perm`, `closeVert_st` (the type-1 `closeVert` unwrap/reopen/merge/fold, lower entries `pre`/`qs` carried), `finishP_st` (P merge; `hP` side conditions only for `isType1 = true`; result `new' ≠ []`), `finishTail_st` (vertex push / merge; vertex-item conditions only for `hasVert = false`) | `StVert.lean` | proved (axioms propext, Classical.choice, Quot.sound) |
| `StTree.lean`: `mem_readStack_of_mem`, `mem_spans_setSides_single`, `finishTree_st` (`finishTree` for a tree edge with `lowval < d`: `closeEars_st`, `L1StInv.mergeLate`, `closeVert_st`, `finishP_st`, `finishTail_st` composed; `hvf` named hypothesis), `finishBack_st` (`finishBack` for a back edge with `lowval < d`), `walkOutPre_st` (`setStackDir d` keeps `DirsOf s d`; the type-1 vertex push appends the `pre` piece), `DirsOf_getD`, `fePState`, `finishEdge_st` (`finishEdge` for a returning edge via `finishEdge_eq`: `qs ++ mid ++ post` with `lowDir = !stackDir[d]`, `sd = stackDir[d]`, above the untouched `B`; `stackDir[k]`, `k ≤ d`, kept) | `StTree.lean` | proved (axioms propext, Classical.choice, Quot.sound) |
| `StBoundary.lean`: `Place.fresh`, `finishRet_frame_st` (frame facts of a returning `finishEdge`; from `finishRet_frame`, `StRetFrame.lean`) (proved); `finishBoundary_st` (block completion at a boundary edge `lowval ≥ d`: the entries above `base` read as the child's pieces and become the block `⟨some (curV, o.dest), stNest ps⟩`, every S/P/R item of them `InBlock`, `tstack = base`, `stackDir`/`hasVert` kept), `rootPop_st` (the root pop after a tree; proved from `StItems.finished` + `RootOK`), `finishBoundary_stLive` (admitted: `StLive s.items sub ⟨some (curV, o.dest), stNest ps⟩` on the pre-state; sim-5: `finishBoundary_st` takes it as the hypothesis `hlive`, and its frame conclusion now also keeps `vs` of the items below the base) | `StBoundary.lean`, `StBdPop.lean`, `CheckStSim.lean` | `finishBoundary_stLive` **admitted** (`Admitted:` docstring; dump-checked by `check_stsim` m8 on seeds 0..1000, 0 violations); `finishBoundary_st` proved from it (axioms propext, sorryAx, Classical.choice, Quot.sound through the admission only); `rootPop_st` proved (axioms propext, Quot.sound) |
| `StBdPop.lean`: `ExpandsList.nil_inv`/`Expands.leaf_inv`/`Expands.node_inv` (moved from `StFinal.lean`), `BdPop` (the item-array change of a boundary step), `BdPop.{below_old, below_new, inBlock, mk'}`, `StItems.bdPop` (`StItems` through a boundary pop given `InBlock` for the popped S/P/R items), `ExpandsList.segment_of_mem`, `Items.Below.expands_cases`, `StRead.complete`, `readStack_nodup_mix`, `StLive` (the open-block invariant), `StLive.vsOrientedAt` | `StBdPop.lean` | proved (axioms propext, Classical.choice, Quot.sound) |
| `StVStart.lean`: `VsIn` (every stack entry's `vStart` satisfies `V`), `VsIn.{tail, cons, merge, modifyCur, modifyNxt}`, `vsIn_{mergeTstackTops, finishTstackTop, maybeUnwrapNxt, loop1Body, loop, finishRest, closeVertTail, closeVert, finishTree, finishBoundary, finishEdge}`, `finishEdge_vStart` (entries after `finishEdge` start at `curV`, at `o.dest` for a returning tree edge, or at an existing start) | `StVStart.lean` | proved (axioms propext, Quot.sound) |
| `StNodup.lean`: `pieceItems`, `stNest_perm` (`stNest ps ~ pieceItems ps`), `Prov`/`Prov.ne`, `refOut_ret_items`, `refOut_items_of`, `refTree_items`/`refOuts_items` (the pieces and blocks of `refTree g t d dirs` are `Nodup` and hold only vertex items of `t.verts` / edge items of `t.edges`), `refOrder_nodup_of_forestOK`, `refOrder_nodup_of_perm`, `refOrder_nodup` (`StFinal.lean`) | `StNodup.lean`, `StFinal.lean` | proved (axioms propext, Classical.choice, Quot.sound) |
| `StSimLemmas.lean`: `refOuts_snoc`, `refOut_hv_false`/`refOuts_hv_false`, `DirsOf_{length, succ, take, congr}`, `frameBlocks_{append, congr}`, `simBlocks_{snoc, congr, DirsOf}`, `DfsOut.{vertsList, edgesList}_append`, `segsStack`/`SegRead`, `mem_readStack_{append, segsStack}`, `SegRead.{congr, cons}`, `StRead.nil` | `StSimLemmas.lean` | proved |
| `StInduct.lean`: `KeepsSize`, `stOut_step` (one `walkOut`: back / boundary tree / returning tree edge via `walkOutPre_st`, `finishBoundary_st`, `finishEdge_st`), `stOuts_nil`/`stOuts_cons`, `stTree_node`, `stWalk` (the `walkOuts`/`walkTree` mutual induction, `StSimOuts` → `StSim`), `DirsOf_zero`, `refBlocks_snoc`, `stForest` (the forest: `stWalk` per tree, `rootPop_st`, `StItems` over `refBlocks g (pre ++ forest)`) | `StInduct.lean` | proved modulo the `StBoundary.lean` admissions and the ear inputs |
| `StFinal.lean`: `no_cycle`, `chain_lt_size` (parent chains in a tree with unique parents are shorter than `items.size`), `ExpandsList.nil_inv`, `Expands.{leaf_inv, node_inv}`, `leaves_eq_of_expands`/`leaves_size_eq` (`Expands i L` ⇒ `L = Items.leaves items items.size i`), `ExpBelow`, `expandsList_ne_nil`, `spr_of_expBelow`, `spr_ch_ne_nil`, `expands_ne_nil` (an S/P/R item's expansion is nonempty), `collapseRuns_{cons_ne, replicate_append}`, `ExpandsList.expands_of_mem`, `restrict_mid`, `restrictCh_eq` (`restrictCh` of a `Nodup` order containing the expansion of `ch i` as a segment is `ch i`), `Items.WF.no_cycle`, `stItems_init`, `walk_sim` (`stForest` at `WalkState.init`), `walk_inBlock` (empty final stack ⇒ every S/P/R item is `InBlock` of a reference block) | `StFinal.lean` | proved (axioms propext, Classical.choice, Quot.sound for the generic lemmas; `walk_sim`/`walk_inBlock` modulo the admissions) |
| `StIParent.lean`: `SpanAll`, `IOk` (vertex items `V`, edge items `Q`, child ids in range, stack span items are not `I`, every `I` child has a `Q` parent), `IOk.{tail, pop, merge, modifyCur, modifyNxt, push, push_vert, push_edge, modify_vs, modify_ch, modify_edge_ch, …}`, `iOk_{mergeTstackTops, finishTstackTop, maybeUnwrapNxt, loop1Body, loop, finishTail, finishRest, closeVertTail, closeVert, finishTree, finishBoundary, finishEdge, walk_aux, walkForest, init}`, `walk_i_parent'`, `walk_i_parent` (`StFinal.lean`; statement change: under `g.WF`, `OrderOK g.nv vo`, `OrderOK g.ne eo`) | `StIParent.lean`, `StFinal.lean` | proved (axioms propext, Classical.choice, Quot.sound) |
| `StRetFrame.lean`: `spansCount_pos_of_mem_readStack`, `lt_of_mem_ch`, `chCount_pos_of_mem_ch`, `mem_getSide'`, `UOk`/`UOk.{head, head_getSide}`, `Fr` (`U`-items keep type/children, non-`U` items get no `U` child), `Out` (entries above `B` are `U`-free), `FrO`, `Fr.{refl, of_items, push, modify_keep, modify_ch}`, `Out.{of_tstack, length, cons, merge, modifyCur, modifyNxt, head, nxt}`, `length_{mergeTop, modifyCur, modifyNxt, maybeUnwrapNxt, finishTstackTop, loop1Type, loop1Body, mergeLate}`, `wp_cond_result`, `wp_loop_cond`, `wp_loop_or`, `loop3Cond_result`, `after_finishTail_true`, `FrO.{of_eq, merge, merge', cons, unwrapNxt, finishTop, l1Type, l1Body, ears, late, finP, tail, vPre, vUnwrap, closeV, finishEdge_ret}`, `spansCount_append_eq`, `below_cases`, `finishRet_frame` | `StRetFrame.lean`, `StBoundary.lean`, `StInduct.lean` (`stRet_finish` takes `Full`) | proved (axioms propext, Classical.choice, Quot.sound) |
| the open-block invariant `StLive` along the walk (the orientation core of `finishBoundary_st`, i.e. `finishBoundary_stLive`): `check_stsim` m8/m10/m11 (seeds 0..1000, 0 violations); site-level terminal facts m13–m17 (`truncCtxB`/`truncSiteB`, 0..1000) | `CheckStSim.lean` (`truncItemsB`, `truncCtxB`, `truncSiteB`) | checked; preservation open — to be carried by the backbone induction (§7.6 "StLive obligations"); the per-entry terminal facts m9/m12/m18 are false as stated (`truncEntryChecks`; seeds 554/6, 0, 9, 154/195/442); the popped-ear shape turned out to be the existing `EarFinish` fields `back_nil`/`bd_bridge`/`bd_comp` |
| `StOpen.lean` (sim-5): `segment_insert`, `Precedes.insert`, `Oriented.insert`, `seqOf`, `InBlock.insert` (`InBlock` is kept when fresh items are spliced around the item's block context), `ItemsFrame`, `HangingUnder.iff_of_frame`, `StLive.{frame, nil, of_subset, append, cons, insert, of_leafTypes}` | `StOpen.lean` | proved (axioms ⊆ propext, Classical.choice, Quot.sound) |
| `StOpenBlock.lean` (sim-5): `retPieces`/`childPieces`/`refOut_ret_pieces`, `openBlockP`/`openBlock` (the truncated reference's block of an open segment), `openBlock_snoc_bd`, `openBlock_snoc_ret`, `frame_ctx`, `TreeFrames`, `openBlockP_ctx`/`openBlock_ctx` (`openBlock g fs dirs ps = ⟨r, A ++ stNest ps ++ B⟩`), `stNest_snoc`, `FrameOK`/`FrameOK.mono`/`TreeFrames.of_ok` | `StOpenBlock.lean` | proved (axioms propext) |
| `StKeepSv.lean` (sim-5): `KeepsSv`/`KeepsSvBelow` through every walk primitive, `walk_keepsSvBelow`, `walkTree_stackVerts_below`, `walkOut_stackVerts_le` | `StKeepSv.lean` | proved |
| reading a tstack as pieces: `readStack`, `stNest_append`, `readStack_pushTstack`, `readStack_mergeTstackTops`, `readStack_fold`, `readStack_finishTstackTop`, `readStack_reopen`/`readStack_modifyNxt_reopen` (per-primitive steps of the simulation relation `readStack stack = stNest pieces`, up to `expandItem` at closes and reopens) | `StRef.lean` | proved |

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

The second gluing pass adds `SpqrTree.PieceSep` (`PieceSep.lean`): `Touches i v` means an
original edge below `i` is incident to `v`; distinct root children have disjoint `Touches`,
distinct children of `V(v)` touch only at `v`, and a block-root `Q(e)` with children `[c,w]`
(`w` of type V) touches an outside edge only at an endpoint of `e`.
`spqrTree_pieceSep` (`WalkPieceSep.lean`), under graph well-formedness and valid vertex/edge
orders, is threaded through the gluing steps and fold; it is now an assembly
(`RelabelOK.pieceSep`, `RelabelPieceSep.lean`) over the item-level facts plus five named
final-items hypotheses (below).
`lake build check_piece_sep` builds the executable checker. Seeds 0..300 from `gen.py`,
plus a single edge, a three-leaf star, a self-loop and four isolated vertices, give zero
violations of all three fields. These checks are empirical evidence, not a proof of the admission.

`IsPlanarEmbedding.append` (`Proofs/PlanarAppend.lean`) proves the same-vertex-numbering
union needed by F: from embeddings of `es₁` and `es₂` and disjoint incident vertex sets,
`rs₁.union rs₂` embeds `es₁ ++ es₂`. It transports `IsPlanarEmbedding.union` with
`f x = if x < n then x else x - n`; the disjointness hypothesis proves injectivity only on
incident vertices, so isolated vertices need no special case.

The F step was also false with the original `GluedUpTo`, independently of separation.
On the actual single-edge tree for `nv=2`, `edges=[(0,1)]`, the items are
`F, V(0), Q(0), I, V(1)`. Set `rotAdj=[some 1,some 0,none,none]`,
`outerE[1]=[none,none,some 2,some 3]`, and all other rows to four `none`s.
The old invariant holds at index 1, witnessed by `[1,0,3,2]`, but F only closes slots 0/1;
the state is unchanged and the root cannot account for its unmatched quarter-edges 2/3.
`Proofs/PlanarEmbedCounterexample.lean` proves both facts without `sorryAx`.
The old invariant is retained as `GluedPieces`; `GluedSlots` extends it with
four-slot row sizes and exposed-slot support: slots are below 4, F exposes nothing,
and an item whose parent is F or V exposes only slots 0/1.
The initialization and leaf proofs preserve these fields; `badState_excluded` proves
the counterexample is excluded. `CheckPieceSep` checks these fields after every item
of every successful planar fold. Seeds 0..300 and the four tiny cases again give zero violations.
`GluedAttachments` extends `GluedSlots` with `outer_at_vertex`: exposed ends of a child
of a V item for original vertex `v` must be at `v` in `g.edges`.
The Q/node steps still need endpoint / cofacial audits; these corrections are not
a claim that their current statements are sufficient.

The attachment condition is necessary already on the two-edge star
`nv=3, edges=[(0,1),(0,2)]`. Its items are `F,V(0),Q(0),I,V(1),Q(1),I,V(2)`.
Set `rotAdj=[1,0,-,-,5,4,-,-]`, expose `[2,3,-,-]` at item 2 and `[6,7,-,-]`
at item 5, and leave all other outer rows empty. `GluedSlots 2` holds, but processing
V item 1 links quarter-edges 3 and 6, at original vertices 1 and 2 respectively.
No `SameVertex` rotation can agree with that link, so `GluedSlots 1` fails.
`Proofs/PlanarEmbedVCounterexample.lean` proves `badVState_before`, `badVState_after`,
and `badVState_excluded` with only standard axioms. Its literal base tree was also
compared by execution with `starGraph.planarSpqrTree false [0,1,2] [0,1]` and matched;
this comparison is empirical, not a kernel proof of tree well-formedness.
`outer_at_vertex` excludes that state and is preserved by the initialization, leaf,
and F proofs. The checker validates it after every item as well: seeds 0..300 and
the four tiny graphs all give zero violations.

The attachment condition alone still permits a fully closed pair of child pieces.
On the same star, set `rotAdj=[1,0,3,2,5,4,7,6]` and leave every outer row empty.
Each child has a planar one-edge rotation, so `GluedAttachments 2` holds, but V does
no linking and leaves eight vertex orbits instead of `2 * 3 = 6`.
`closedVState_before` and `closedVState_after` prove this failure with only standard
axioms. `GluedOriented` extends `GluedAttachments` with `outer_present`: a processed
nonempty piece whose parent is V has an exposed end. It also records `outer_dir`,
the equality of the quarter-edge direction and slot parity needed by the specified
1-sum transposition. Initialization, leaf, and F preservation are proved; the
closed counterexample is excluded by `closedVState_excluded`.

V's executable loop cannot ignore an empty child after starting its chain: such a
child clears the last endpoint. The structural `PieceSep.v_nonempty` clause says
every child of V contains an original edge; together with `outer_present`, the
paired-end invariant, and slot support this supplies slots 0/1 for each child.
The new structural clause and both boundary conditions were checked before use,
on seeds 0..300 and the single-edge/star/self-loop/isolated-vertices cases, with
zero violations. `spqrTree_pieceSep` is discharged (`RelabelPieceSep.lean`) from `Items.WF` +
`Items.Ranges` through `RelabelOK` — every field as a theorem `RelabelOK.pieceSep_<field>`
(standard axioms), assembled by `RelabelOK.pieceSep`; `walk_ranges'` (`WalkItemsWF.lean`) supplies
`Items.Ranges` on the final items for every graph (`ranges_initialItems` for `g.nv = 0`). Node-side
transport: `nvOf_iff`/`nodeVerts_get'`/`nvList_mid_iff`/`mem_nvList_iff` (node-verts are the
`nvList` positions, V children its interior), `capEnd_iff`, `virt_edge` (exact `nvs`/`neOrig` of
each virtual edge from `LayoutFacts.s_edge`/`p_edge`/`r_edge`), `nodeEdge_node`
(`RelabelAll.nodeEdge_at`), `v_touch_eq`/`sep_vs`/`node_vertex_outside` (`Endpoints.separation`,
`interior`, `Ranges.vs_att`), `below_antisymm`/`not_below_parent` (`Tree.acyclic`). Five walk-side
facts about the final items are consumed: `Items.QUpper` (a child of `vertItem v` has
`(vs).1 = some v`), `Items.RootSep` (distinct root children share no vertex), `Items.RootV` (root
children are V items), `Items.QChildVs` (a block-root Q's non-V child has `vs = (vs_Q.1, some w)`,
or `(vs_Q.1, none)` for a loop), `Items.PChildVs` (a P's non-V children carry its `vs`). They are
the fields of `PieceFacts` (`RangesPiece.lean`, §4.6), whose walk-state invariant `PieceInv d s` is
kept by every primitive and by `finishEdge_piece`; the final-state form is the single named
hypothesis `walk_pieceInv` (`WalkPieceSep.lean`), with `walk_q_upper`/`walk_root_sep`/`walk_root_v`/
`walk_q_child_vs`/`walk_p_child_vs` its projections. Checker: `Ranges.checkFinal` on the final
items (`ranges.q_upper`/`root_sep`/`root_v`/`q_child_vs`/`p_child_vs`) and `checkPieceInv`/`checkPiece`
at every walk site (`piece_*`), 0 violations on seeds 0..3000 × both modes. `#print axioms
spqrTree_pieceSep` reaches `sorryAx` only through `walk_pieceInv` and `walk_items_wf`/`walk_ranges`'
existing admissions.

`embedItem_step_F` is proved in `PlanarEmbedF.lean` with the strengthened invariant.
`only_root_F` reduces the item to index 0; `closeList_get` identifies each
child's final rotation entries with its single-pair closure, with all other children
framed out by `root_disjoint`. `Piece.union_embeddings` iterates the same-numbering
append theorem; `edgesBelow_eq` identifies its concatenation with the root piece.
The root outer row stays unset, and the row-size / slot-support / attachment fields are preserved.
`PlanarEmbedFold.lean` imports this proof instead of the former F admission.
The F theorem and every helper in this group were audited with `#print axioms`:
only `propext`, `Classical.choice`, and `Quot.sound` occur, never `sorryAx`.

`embedItem_step_V` is proved in `PlanarEmbedV.lean`. `embedItem_V` identifies the
executable loop with `vLoop`; `vLoop_open_spec` iterates `Piece.OpenEmbedding.splice`,
framing later children and all other maximal pieces by edge disjointness.
`vLoop_children` handles the empty child list with the empty embedding and the
nonempty case with `v_child_boundary`. `setOuterPair_lookup` exposes only slots 0/1,
and `v_parent_not_v` discharges the attachment/presence conditions for the new V row.
The V theorem and every loop/frame helper were audited with `#print axioms`:
only standard axioms occur. Q and node remain admitted.

Q also needs the boundary of a V item itself to be at that item's original vertex.
On `nv=3, edges=[(0,1),(1,2)]`, the actual tree is the chain
`F,V(0),Q(0),I,V(1),Q(1),I,V(2)` (Q(1)'s parent is V(1)).
Set `rotAdj=[-,-,-,-,5,4,-,-]`, expose `[6,7,-,-]` at V item 4,
`[4,5,-,-]` at Q item 5, and leave all other rows empty. `GluedOriented 3`
holds, but Q item 2 links quarter-edges 3 and 6, at original vertices 1 and 2.
`badQState_before`, `badQState_after`, and `badQState_excluded` in
`Proofs/PlanarEmbedQCounterexample.lean` are kernel proofs with standard axioms only.
The literal tree's projection matched `pathGraph.planarSpqrTree false [0,1,2] [0,1]`
by execution; this comparison is empirical, not a kernel proof of well-formedness.

`GluedVertex` extends `GluedOriented` with `outer_vertex`: a V(v) item's exposed
ends lie at v, occupy only slots 0/1, and exist whenever the processed piece is
nonempty. The latter clauses let Q consume the V child's entire open boundary
and prevent a nonempty closed child from being silently skipped.
Initialization, leaf, F, and V preservation are proved and re-audited with standard
axioms only. For V, the fold's `OpenEmbedding.boundary` localizes both outer ends
at v, and its empty-output case has an empty edge list. All three clauses pass
the extended `check_piece_sep` on seeds 0..300 and the four tiny cases.
The cap endpoint condition is independently necessary on two parallel edges
`nv=2, edges=[(0,1),(0,1)]`, whose tree is `F,V(0),Q(0),Q(1),V(1)`.
Set every rotation entry to `none`, expose `[6,7,4,5]` at capped Q item 3,
and leave every other outer row empty. `GluedVertex 3` holds (the child has the
one-edge rotation), but Q item 2 links quarter-edges 1 and 6 at vertices 0 and 1.
`badCapState_before`, `badCapState_after`, and `badCapState_excluded` in
`Proofs/PlanarEmbedCapCounterexample.lean` prove this with standard axioms only.
The literal tree was also compared with the executable output.

`GluedUpTo` extends `GluedVertex` with `outer_cap`: slots 0/1 are at the first
original cap endpoint, slots 2/3 at the second, and a processed nonempty capped
piece exposes all four slots. Initialization, leaf, F, and V preserve the field
with only standard axioms; F/V have no cap themselves and frame all older rows.
The checker validates cap endpoints and slot presence on seeds 0..300 and the
single-edge/star/self-loop/isolated-vertices cases, with zero violations.
The capped-node-child case of Q still needs its cofacial contract audited;
these endpoint corrections do not claim that Q's current statement is sufficient.

**Cofacial gap (kernel-checked).** `PlanarEmbedFaceCounterexample.lean` now
exhibits the remaining obstruction on three parallel edges. The literal tree's
printed representation matches the executable output. Its P child has two real
edges and the planar rotation `[5,4,7,6,1,0,3,2]`. The bad exposed row is
`[8,5,10,7]`: both pairs have the correct endpoint and direction and are facing
pairs of that rotation, but the local slot-0 and slot-2 ends are `0` and `2`, in
different face orbits. `badFaceState_before` proves the entire current `GluedUpTo`
precondition; `badFaceState_after` proves failure even of `GluedPieces` after Q.
The forced closed rotation then has two face orbits instead of the six required
by Euler's formula. `doubleRot_not_cofacial` proves the boundary obstruction.
All four lemmas have standard axioms only.

The needed restatement is: for each **maximal capped piece**, its planar witness
`ρ` in `GluedUpTo.piece` must additionally satisfy
`ρ.SameFaceOrbit la lc`, where `loc a = some la`, `loc c = some lc`, and `a`, `c`
are slots 0 and 2. This must refer to the same `ρ` that agrees with `rotAdj` and
the two facing pairs, not a separately existential embedding. The extended
`check_piece_sep` closes both boundary pairs, localizes the rotation, and checks
this condition on each maximal capped piece; seeds 0–300 and the four tiny cases
pass. The restatement is `GluedFaces` (`PlanarEmbedFaces.lean`): it extends `GluedUpTo`
with `cap_face`, which for every maximal capped item `j` gives one `ρ` satisfying both
`PieceWitness g s j ρ` (planar embedding of `pieceBelow j`, agreement with `rotAdj`,
unset ↔ exposed, both exposed pairs facing pairs of `ρ`) and `CapFace g s j ρ`
(`ρ.SameFaceOrbit la lc` for the localized slots 0 and 2). The old `GluedUpTo` is kept
unchanged, so `badFaceState_*` remain theorems about the insufficient precondition.
`gluedFaces_init` and the leaf/F/V steps (`PlanarEmbedFacesSteps.lean`) reuse the
`GluedUpTo` steps and frame the untouched maximal pieces (`pieceWitness_frame`,
`gluedFaces_of_frame`); F's root and V's new piece carry no cap.

**Q step (proved).** `PieceSep` gained `q_shape` (a `Q` item's children are `[]`, `[O]`,
or `[c, V w]` with `c` a capped non-`V`/non-`O` item — the earlier `I/S/P/R` form was false:
random trees have capped `Q` children), `cap_orig` (cap edges have original endpoints),
`cap_nonempty` and (for the node step) `twin_orient` (twin node-edges list their original
endpoints in the same order); all checked on seeds 0..300 and the tiny cases with zero violations.
`Proofs/OrbitSplit.lean` proves `conj_numFaceOrbits_split`: conjugating by two cofacial
pairs adds exactly two face orbits. `Proofs/PlanarInsert.lean` defines
`RotationSystem.insert ρ x a0 a2` (`single.union ρ` followed by the two `conj`s the
executable `link` sequence performs) and proves `IsPlanarEmbedding.insert`: hanging an
edge between two exposed ends on a common face of `ρ` is planar. `Proofs/PieceInsert.lean`
lifts this to pieces: `Capped` (same-witness certificate of a capped child) and `Open2`
(two exposed facing pairs), with `Capped.insert`, `Open2.close`, `Open2.splice`,
`single_open2`. `PlanarEmbedQExec.lean` unfolds the executable `Q` branch exactly
(`embedItem_Q_nil/_O/_cons`); `PlanarEmbedQ.lean` proves the three cases —
`q_leaf_glued` (no child: `setOuter4` exposes the four quarter-edges, the cap face is
the single-edge rotation), `q_loop_glued` (`O` child: the edge is a loop,
`RotationSystem.loop`), `q_block_glued` (`[c, V w]`: `Capped.insert` on the child's two
cofacial pairs, then `Open2.close`/`Open2.splice` with the lower `V` piece by
`q_lower_boundary`) — and `embedItem_step_Q_faces` dispatches on `q_shape`. Every
theorem of this group has only `propext`, `Classical.choice`, `Quot.sound`.

**Fold.** `PlanarEmbedFacesFold.lean` dispatches `embedItem_step_faces` (leaf, F, V, Q
proved) and folds `gluedFaces_planarEmbed`; `planarEmbed_sound` now uses it (through
`GluedFaces.toGluedUpTo` into `glued_root`). The S/P/R step `embedItem_step_node_faces`
(`PlanarEmbedNode.lean`) is proved from `embedItem_node` and two lemmas: `node_capped_glued`
(proved: `GluedFaces g i` from a `Piece.Capped` certificate of the node's piece whose four exposed
ends fill the node's row, `rotAdj` framed outside the piece, other rows unchanged — the node's
`PieceWitness`/`CapFace` are read off the certificate, the row's vertices from `pair0`/`pair2` and
`SameVertex`, the parent-type side conditions from `PieceSep.node_parent`) and `nodeFold_capped`
(the one remaining gluing admission: the fold of `nodeStep` over the node's quarter-edges frames
`rotAdj` outside the piece and produces such a certificate at the cap's endpoints). Its executable side is done: `PlanarEmbedNodeExec.lean`
unfolds the node branch exactly (`embedItem_node`: the loop over `[4·neSt : 4·neEn]` is
`List.foldl (nodeStep i neSt)` over `List.range'`, where `nodeStep` links the exposed slots
`treeQe` of the children owning the twins of `ta` and its rotation successor `tb > ta`, through
the vertex item's exposed pair at a `(ta % 4 = 2, tb % 4 = 1)` corner, and exposes the slot
opposite each cap corner), with frame lemmas `nodeStep_*`/`nodeFold_*` (`outerE` of other items,
`outerE`/`rotAdj` sizes, `exposedAt` of other items unchanged). What remains is the semantic
part, in three pieces: (N1) a `TwoSum`-style splice of a planar rotation system *containing* the
virtual edge `e` (the node's `nodeRot`, hypothesis `hloc`) with a *capped* piece (cap absent, its
four exposed slots in two facing pairs, slots 0/2 cofacial — `Piece.Capped`), concretely equal to
the `link`s `nodeStep` performs and planar: this is `TwoSum.splice_isPlanarEmbedding` after
re-inserting the cap with `IsPlanarEmbedding.insert`, or a direct Euler count as in
`Proofs/PlanarGlue.lean`; (N2) iterating (N1) over the node's non-cap edges in `nodeRot` order
with the `Piece.loc` bookkeeping (`ves` of the node = concatenation of the children's `edgesBelow`),
where the twin ↔ exposed-slot correspondence is `Twins.noncap_children` (the twin of a non-cap
node-edge is the cap of the corresponding child) plus `PieceSep.twin_orient` (`neOrig ne = neOrig
tw`: the node's quarter-edge `(side, dir)` of `ne` meets the child's exposed slot `2·side+dir`,
checked in `CheckPieceSep.lean` on seeds 0..300 and the tiny cases); (N3) the inner `V` items at
corners (1-sum through `OpenEmbedding`, as in `Open2.splice`), then the cap's four neighbours
become the node's exposed slots and are cofacial because the cap edge's two sides were distinct
faces of `nodeRot` (2-connected skeleton) — the `CapFace` certificate of the node. The local embedding is the hypothesis `hloc` of the node
step/dispatcher/fold (`IsPlanarEmbedding (localSkeleton i) (nVerts i) (nodeRot i)`, defs moved
to `PlanarNodeSpec.lean`); `planarEmbed_sound` discharges it with `nodePlanar_sound` and
`isPlanar_of_all` from `t.nodePlanar.all id = true`. That needs `i < nodePlanar.size`, i.e.
`planarTree_nodePlanar_size : nodePlanar.size = size` for `g.planarTree …` — a fact
`planarEmbed_sound` genuinely depends on (with a short flag array `hall` would be vacuous); it is
`RotInv.np_size` (`PlanarRotInv.lean`), threaded through `Frame`/`RotInv.step` like `ne_size`.
Facts recorded for (N2)/(N3) in `PieceSep` (checked by `CheckPieceSep.lean`, seeds 0..300 + tiny
cases, 0 violations): `node_parent` (an S/P/R node hangs under a Q or a node), `nv_child` (a
node-vertex's V item is a child of the node iff the node-vertex is not a cap endpoint),
`v_child_nv`, `nv_vert_inj`. `PlanarNodeSpec.NodeCorners` states the corner structure the node
step relies on (`CornerAt i nv ta`: the `(side 1, dir 0)` quarter-edge of a non-cap edge at `nv`
whose rotation successor is a `(side 0, dir 1)` quarter-edge; cap endpoints have none, every inner
node-vertex exactly one with `ta < tb`); `planarTree_nodeCorners` is a named admission (not
reached by `planarEmbed_sound` yet; it will be a hypothesis of `nodeFold_capped` once used).
`Proofs/PlanarDegree.lean` (`lbl_stepC`, `sameOrbit_of_lbl`, `rot_div_ne`) and
`OrbitSplit.conj_not_sameOrbit_split`/`insert_not_sameFaceOrbit` are proved building blocks for
(N1)/(N3). `Proofs/PieceJoin.lean` has the S-node building block `Capped.join` (standard axioms):
two `Capped` pieces meeting exactly in the vertex `v` (caps `u v` and `v w`), glued by the node
step's two links `c2 ↔ d1`, `d0 ↔ c3`, form a `Capped` piece with ends `c0 c1` at `u`, `d2 d3` at
`w`; the witness is `(ρ.union σ).conj l2 (ρ.size + m0)` (as in `IsPlanarEmbedding.splice`), and
`c0`, `d2` stay cofacial because both image swaps of `conj_stepC` merge orbits of the two
components (the outer faces of the two pieces become one face).
`Capped.attachOpen` (same file, standard axioms) is the S-corner variant: a `V` item's open piece
(`OpenEmbedding`, boundary pair `w0 w1` at `v`) hung on the `v`-side of a capped piece by the
node step's corner link `c2 ↔ w1` gives a capped piece with `v`-side pair `w0 c3`, `c0`/`w0`
cofacial; the executable corner (`link qa outerE[v][1]`, then `link qb outerE[v][0]`, then the
`slot 3` link) is `attachOpen` followed by `join`. `PlanarEmbedNode.capped_child` (standard
axioms) generalizes `q_capped_child` to any capped non-`I`/`O` child of any item: in a
`GluedFaces g (i + 1)` state its row is filled and its piece is `Capped` at its cap's original
endpoints (`PieceSep.cap_orig`).
The P-node building blocks (all standard axioms): `IsPlanarEmbedding.ident_conj`
(`Proofs/PlanarIdent.lean`) identifies two already-connected vertices of one embedding by a
`conj` of two cofacial pairs (Euler: the face count grows by two, the vertex count drops by
one, the edge count is unchanged); `IsPlanarEmbedding.splice2` (`Proofs/PlanarSplice2.lean`)
is the 2-sum on the common vertex space: embeddings of `es₁`, `es₂` sharing exactly the vertices
`u ≠ v`, both connecting `u` to `v`, glued by `(rs₁.union rs₂).conj a (rs₁.size + b)` at `u` and
then `conj c (rs₁.size + d)` at `v`, embed `es₁ ++ es₂` provided the two `v`-pairs are cofacial
after the first conjugation (`oneSum_conj`, then `ident_conj` on the shifted copy, then the
relabelling `x ↦ if x < n then x else x - n`, injective on incident vertices by the separation
hypothesis). `Proofs/OrbitSegment.lean` has the orbit bookkeeping for sequential
conjugations: `iterate_eq_of_eqOn`, `IsPermOn.exists_first_hit` (the first iterate reaching
`b`, with no earlier return to `a`), `RotationSystem.sameOrbit_stepC3_rot` (the face orbit of
`rot a`, `rot b` mirrors that of `a`, `b`). `Capped.parJoin` (`Proofs/PieceParJoin.lean`) is
the P-node step: two capped pieces with the same cap endpoints `u, v`, meeting exactly in `u`
and `v`, glued by the node step's two links `c1 ↔ d0` (at `u`) and `c2 ↔ d3` (at `v`), form a
capped piece with pairs `c0 d1` at `u` and `d2 c3` at `v`; `c0` and `d2` stay cofacial because
the face walk from `l0` reaches `l3 ^^^ 3` before any point moved by either conjugation
(`exists_first_hit` + `iterate_eq_of_eqOn`), and the second conjugation sends `l3 ^^^ 3` to
`d2`'s slot.
Route for `nodeFold_capped` (not yet formalised): the fold's links on the children's exposed
slots are the local rotation `nodeRot` transported along `treeQe` (the skeleton's quarter-edges
never enter `rotAdj`). S nodes iterate `Capped.join`/`attachOpen` along the cycle, P nodes
iterate `Capped.parJoin`; the cap edge is never inserted, its two corners become the node's four
exposed slots. For R nodes the general step is `TwoSum.splice_isPlanarEmbedding` with `es₁` the
hybrid skeleton (processed children substituted for their virtual edges), `es₂` the child's piece
with its cap re-inserted by `RotationSystem.insert` between its two cofacial exposed slots
(`WF.face` from `insert_not_sameFaceOrbit`, `WF.conn` from `edgesConn_of_sameOrbit`); the final
`Capped` certificate then needs the inverse of `insert` (delete the cap edge, pair its four
neighbours, the two sides merge into one face) — the dual lemma `IsPlanarEmbedding.uninsert`
below. Its hypothesis "the cap's two sides are distinct faces" is the classical edge-on-a-cycle
fact `IsPlanarEmbedding.cap_sides_distinct`, which needs the genus inequality
`numFaceOrbits + 2·numNonIsolated ≤ 2·(2·numComponents + es.length)` for arbitrary valid rotation
systems (`twoSum_planar` assumes `face ∨ conn`). The S/P routes via `join`/`parJoin` do not
need it.
Genus inequality (`RotationSystem.genus_le`, `Proofs/PlanarGenus.lean`, standard axioms): for every
`IsEmbedding es n rs`, `rs.numFaceOrbits + 2·numNonIsolated es n ≤ 2·(2·numComponents es n +
es.length)`. Induction on `es`, deleting the first edge `p` of `p :: es`: each end of `p` whose
rotation leaves the edge is detached by one `conj` (`DetachAt`: pairs `p`'s two quarter-edges at
that end with each other and its two former neighbours with each other; `+2` vertex orbits via
`conj_orbitCount_split`, `±2` face orbits via `conj_numFaceOrbits_split`/`_merge`), after which the
system is literally `single.union ρ` or `loop.union ρ` (`eq_union_shrink`) with `ρ = shrink _` an
`IsEmbedding` of `es` (`isEmbedding_of_union`); the inductive bound for `ρ` transfers
(`genus_cons_bound`). Graph side (`Proofs/PlanarCompCons.lean`, via the component minima of
`ccCount`): `numNonIsolated_cons_eq` (exact count of endpoints isolated in `es`),
`numComponents_cons_le_succ` (`C(es) ≤ C(p :: es) + 1`), `numComponents_cons_succ_le` (both
endpoints isolated: `+1`), `numComponents_cons_le_of_isolated`/`_of_loop` (no drop), and the
bridge fact `edgesConn_of_not_cofacial`: if the two sides of `p` are not cofacial, the face walk
along one side connects `p.1` to `p.2` inside `es` (`edgesConn_tail_across`), so
`numComponents_cons` applies. Case split: an end whose rotation stays inside `p` is isolated in
`es` (`not_hasEdge_of_pair`, from the vertex-orbit labelling), both ends detached gives the
`detach_both` dichotomy "`F + 4` or (`F` and sides not cofacial)". Checked beforehand on 4000
random embeddings by `checks/PlanarGenusCheck.lean`.
Edge deletion (`Proofs/PlanarUninsert.lean`, standard axioms): `detach_both` additionally returns
the rotation of the detached system above the edge (`rot 0 ↔ rot 1`, `rot 2 ↔ rot 3`, identity
elsewhere) and, when the sides are distinct faces, the cofaciality of `rot 1` and `rot 3` in it
(`sameOrbit_detach_merge`: the four `swapImg`s of the two `conj_stepC`s merge the two side faces
and split off the isolated edge's own face). `IsPlanarEmbedding.uninsert`: if `p :: es` has a
planar embedding `rs`, `p.1 ≠ p.2` and `EdgesConn es p.1 p.2` (not a bridge), then the sides `2`,
`0` of `p` are distinct faces, all four neighbours `rot 0..3` lie above the edge (a neighbour on
the edge itself would make an endpoint isolated in `es`, `not_hasEdge_of_pair`, or force
`p.1 = p.2`), and `σ = shrink H` is a planar embedding of `es` (Euler: `genus_le` for `σ` gives
`C(es) ≥ C(p :: es)`, `numComponents_cons_le_of_conn` gives `≤`; the cofacial case of
`detach_both` would give `F + 2` for `σ`, contradicting `genus_le`) with the explicit rotation
`σ.get q = some (…)` above and `rot 1 - 4`, `rot 3 - 4` cofacial (transported through
`union_stepC`/`sameOrbit_unionStep_right`/`sameOrbit_shift_sub`). `cap_sides_distinct` is the
first conjunct restated as `¬rs.SameFaceOrbit 0 2`. Checked beforehand (`delFirst` in
`checks/PlanarGenusCheck.lean`: embedding, Euler, components, neighbour pairing, cofaciality of the
neighbours; 410 deletions, 0 failures — gating on connectivity in `es`, not on non-isolation,
is essential: bridges fail).
Node-fold bookkeeping (towards `nodeFold_capped_S/P/R`): `PieceSep` gains the checked structural
facts `node_child_cap` (a non-`V` child of an S/P/R node is a capped non-`I`/`O` item),
`nv_orig_inj` (distinct node-vertices of an S/P/R node have distinct original vertices),
`node_touch` (a child touches the original vertex of a node-vertex only where it is attached,
`NvInc`: it is that node-vertex's `V` item, or its cap's twin is a node-edge at that node-vertex)
and `node_attach` (two distinct children meet only at original vertices of node-vertices); all
checked by `check_piece_sep` on seeds 0..300 and the tiny cases, 0 violations.
`Proofs/PiecePerm.lean` (standard axioms): `Piece.Capped.perm` transports a `Capped` certificate
across a permutation of `ves` (the children are glued in skeleton order, `edgesBelow` lists them in
children order), via `glob`/`loc` inverse maps, `IsPlanarEmbedding.reindex` and
`sameFaceOrbit_of_map` (face orbits pull back along a quarter-edge map commuting with `get` and
`^^^ 3`).
`nodeFold_capped_S` (the `S` case of the node fold, standard axioms; `PlanarEmbedNodeS.lean`,
`PlanarEmbedNodeSChain.lean`, `PlanarEmbedNodeSFold.lean`, `PlanarEmbedNodeSMain.lean`): the
executable fold is `blockFold` over the `n` four-quarter-edge blocks of the cycle layout
(`layoutRot .S`, obtained from the new hypothesis `hlay : LayoutAt i` = `neRotAdj_segment`,
threaded through `embedItem_step_node_faces`/`embedItem_step_faces`/`gluedFaces_planarEmbed` and
discharged in `planarEmbed_sound`). Block `0` (the cap edge) is four `setOuter`s exposing slots
`0, 1` of the first non-`V` child and `2, 3` of the last one (`cap_block_S`); block `k`
(`1 ≤ k ≤ n - 2`) is two skips, the corner link at node-vertex `s + k` (`Capped.attachOpen` when
the node-vertex's `V` item is present, read off `nodeStep`'s own `outerE[v][0].isSome` test —
`NodeCorners` is not needed) and the child-to-child link (`Capped.join`), with the exposed slots
transported along `treeQe` (`nodeStep_skip`/`_cap`/`_link`/`_corner`); block `n - 1` is a no-op
(`last_block_S`: all its successors are earlier). The induction invariant `InvS` carries the
`Capped` certificate of the growing chain `chain g i m` (first child, then alternately the `V` item
and the non-`V` child of each interior node-vertex, `sL`), `rotAdj` framed outside the chain and
the other rows unchanged (`invS_zero`, `invS_step`, `invS_all`). The chain's edge list is a
permutation of `edgesBelow i` (every child is a non-`V` child or the `V` item of an interior
node-vertex: `PieceSep.v_child_nv`/`nv_child` with `capEnd_S`; `children_perm_sL`), so
`Capped.perm` gives the certificate of `pieceBelow g i`; the cap's `neOrig` endpoints are the
first and last node-vertex's original vertices (`capOrig_S`). The S-node data layer (`shape_S`,
`nEdges_S`, `capNe_S`, `nvs_S`, `noncap_S`, `neRotAdj_S`, `child_S`, `child_cert_S`,
`sW_boundary`, `nvInc_*`, `sep_*`/`disjoint_*`) is derived from `ChildShape`, `Twins` and
`PieceSep` (`node_child_cap`, `nv_orig_inj`, `node_touch`, `node_attach`).
`nodeFold_capped_P` (the `P` case, standard axioms; `PlanarEmbedNodeP.lean`,
`PlanarEmbedNodePFold.lean`) follows the same template over the bond layout `layoutRot .P` =
`rotP`: the `k - 1` non-cap node-edges are, in order, the twins of the caps of the `k - 1`
children (`noncap_P`, `sC_length_P`, `sC_getElem?_P`; a `P` node has no `V` child, `noV_P`, so
`sC i = children i`), all with the same original endpoints `sV i 0 ≠ sV i 1` (`nvOrig_P`,
`neOrig_P`, `capOrig_P`, `sV01_ne_P`), and two children meet only there (`attach_P` from
`PieceSep.node_attach`). Block `0` exposes slots `0, 3` of the first child and `1, 2` of the last
one (`cap_block_P`); block `m + 1` (`m + 1 ≤ k - 2`) is two skips and the links `a1 ↔ b0`,
`a2 ↔ b3` between children `m` and `m + 1`, i.e. `Capped.parJoin` (`invP_step`, invariant `InvP`
on the chain `pChain g i m` of the first `m + 1` children); block `k - 1` is a no-op
(`last_block_P`). `pChain g i (k - 2) = pieceBelow g i` literally (`pChain_eq_pieceBelow`), so no
permutation is needed. `nodeFold_capped` (`PlanarEmbedNodeFold.lean`) dispatches on the type; the
`R` case is `nodeFold_capped_R` (proved below via `main_R`; it takes the extra hypotheses
`hclosed : NodeRotClosed` and `hcor : NodeCorners`).
`#print axioms planarEmbed_sound` reaches `sorryAx` only through the walk/R-side admissions,
`spqrTree_pieceSep`'s final-items hypothesis `walk_pieceInv`, `planarWalk_planarFinish`,
`planarRelabelTree_relabelNodeR`, `planarWalk_items_wf` (the three ingredients of `nodePlanar_sound_R`),
`planarTree_nodeRotClosed` and
`planarTree_nodeCorners`.

**R node.** Route: the hybrid system. `hloc` is the local skeleton `localSkeleton i` with rotation
`nodeRot i`; `rotR` (`PlanarEmbedNodeR.lean`) identifies the executable `neRotAdj` entries of the
node's segment with `nodeRot i` shifted by `4 * neSt`, which needs `NodeRotClosed`
(`PlanarNodeSpec.lean`: the entries of an `R` node stay inside its own segment — `nodeRot`
subtracts `4 * neSt`, so `hloc` alone cannot say so; admitted as `planarTree_nodeRotClosed`,
checked by `check_piece_sep`'s `rot_closed` on seeds 0..300 and the tiny cases). `child_R` is the
`treeQe` transport of the node quarter-edge `4 (neSt + j + 1) + r` to slot `r` of the `j`-th
non-`V` child (whose cap is the twin of node-edge `j + 1`), `neOrig_R` gives the original
endpoints of every node-edge as `sV` of its node-vertices. The generic side is complete and
independent of the tree (`Proofs/PlanarRGlue.lean`, standard axioms): the abstract model `RG`
(skeleton `cap :: sk` in the original vertex space with rotation `ρ₀`, `nV` corner pieces hung by
`vstep` = one `conj`, `m` capped children substituted by `subst_step` = `RotationSystem.insert` of
the child's cap followed by `IsPlanarEmbedding.subst`), the invariants `InvV`/`InvC` carry the
explicit rotation of every slot through the induction (`csys`), and `RG.Hyp.glue` removes the cap
by `IsPlanarEmbedding.uninsert` (its non-bridge hypothesis from the children's cofacial pairs lifting
`EdgesConn sk` — supplied for the tree by `ThreeConnected.edgesConn_eraseIdx`): the final planar
system on `vflat ++ cflat` has rotation `fin ∘ target` on every piece slot (`fin` pairs the two
slots that faced each cap side), the two pairs `skt 0 ↔ skt 1`, `skt 2 ↔ skt 3`, and `skt 1`,
`skt 3` cofacial. The executable correspondence (`PlanarEmbedNodeRExec.lean`,
`PlanarEmbedNodeRCert.lean`, `PlanarEmbedNodeRGlue.lean`, `PlanarEmbedNodeRMain.lean`, standard
axioms) closes `nodeFold_capped_R`: `RFold` abstracts the node's data (rotation `T = (nodeRot i).rot`
on `[4, 4 nE)`, slot map `Q = rQ i s` = `treeQe` into the children's exposed slots, corner quarter-edges
`Corner` with their `V` item's pair `w0`/`w1`), `InvR`/`nodeFoldR` characterise the fold pointwise:
`rotAdj[Q l] = partner l` (`Q (T l)`, or `w1`/`w0` at a corner), the four cap-facing slots `Q (T c)`,
`c < 4`, keep their old (unset) entry, the corner pairs get `Q l`/`Q (T l)`, everything else is
framed, and the node's `outerE` row is `#[Q (T 1), Q (T 0), Q (T 3), Q (T 2)]`. `rfold_R` supplies
`RFold` from the tree (`slot_R`, `child_cert_R`, `corner_R`/`corner_V_R` from `NodeCorners`,
`capNvs_R`, `capEnd_R`, `sW_R`, `vOf_R`), `rgR` instantiates `RG` with the node's data and `hyp_R`
proves `RG.Hyp` (`rgR_m`, `rgR_Sk`, `rgR_skE`, `rgR_cap`, `rgR_cρ`, `rgR_cl`, `rgR_vρ`, `rgR_vw`).
`pieceR` is the node's piece in gluing order (`vesV ++ vesC`, `pieceR_es`), `locV_sub`/`locC_sub`/
`locQ_R`/`locW_R` transport `Piece.loc` of each child slot to the index `RG.Hyp.glue` assigns to it
(`vIdx`/`cIdx`/`res`/`skt` shifted by the deleted cap), `agrees_V_R`/`agrees_C_R` identify `glue`'s
rotation with the fold's `rotAdj` on every piece slot (`tgtV_ta_full`/`tgtV_tb_full`/`tgtV_other`
for the corner redirections), `unset_R` identifies the unset slots with the four cap-facing ones, and
`main_R` packages planarity, `pair0`/`pair2` (`vert_slot_R` + `same_vertex` of `hloc`), `face`
(`sameFaceOrbit_iff`) and the directions (`dir_slot_R`) into `Piece.Capped` of `pieceR`, moved to
`pieceBelow g i` by `Capped.perm` along `pieceR_perm`; the frame outside the piece is `nodeFoldR`'s
with `memQ_R`/`memW_R`.

`glued_root` is proved. `edgesBelow 0` is a permutation of `range t.ne`, rather
than the identity enumeration claimed by its old docstring. Each item's piece is
included in its parent's by `edgesBelow_eq`; induction along the decreasing parent
chain covers every Q item at the root. Bijections and nodup then give the edge
permutation. With the F row unexposed, the actual rotation is total and transports
along `Piece.loc`. `IsPlanarEmbedding.reindex` preserves the vertex/face orbit
counts by `orbitCount_congr` and the graph counts by equal edge membership.
The generic `glued_root_of` also handles a zero-size tree (its graph is edgeless).
`glued_root_of` and its transport helpers have standard axioms only. The named
`glued_root` obtains `ChildShape` from `spqrTree_childShape`, whose sole admission
is the accepted `walk_items_wf`; its projection and edge-count transport add none.

| statement | file | status |
|---|---|---|
| quarter-edges, `RotationSystem`, `IsEmbedding`, `IsPlanarEmbedding`, `Planar` | `Planar.lean` | def |
| `SpqrTree.PieceSep`, `spqrTree_pieceSep` | `PieceSep.lean`, `WalkPieceSep.lean`, `RelabelPieceSep.lean` | def / **assembled** (`RelabelOK.pieceSep`, every field `RelabelOK.pieceSep_<field>` standard-axiom) over the final-items hypothesis `walk_pieceInv` (`PieceFacts` = `Items.QUpper`/`RootSep`/`RootV`/`QChildVs`/`PChildVs`, kept as `PieceInv` by every primitive and `finishEdge_piece`, §4.6; `check_walkinv` `ranges.q_upper`/`root_sep`/`root_v`/`q_child_vs`/`p_child_vs`, 0 violations 0..3000 × both modes); `sorryAx` only through those and `walk_items_wf`'s admissions |
| `IsPlanarEmbedding.append` (vertex-disjoint edge lists in the same vertex numbering) | `Proofs/PlanarAppend.lean` | **proved** (`propext`, `Classical.choice`, `Quot.sound`) |
| `IsPlanarEmbedding.oneSum_conj`, `IsPlanarEmbedding.splice` | `Proofs/PlanarOneSum.lean`, `Proofs/PlanarSplice.lean` | **proved** (standard axioms); the specified transposition gives the 1-sum rotation, and `splice` keeps the original vertex numbering when the pieces meet only at the attachment vertex |
| `Piece.OpenEmbedding.splice`, `OpenEmbedding.frame`, `v_child_boundary` | `Proofs/PieceSplice.lean`, `PlanarEmbedVBoundary.lean` | **proved** (standard axioms); linking the inner ends of two open pieces agrees with the specified conjugated rotation and leaves exactly the two outer ends open; V children supply these open pieces |
| `embedItem_step_V`, `vLoop_children`, `vLoop_piece`, `vLoop_open_spec`, `embedItem_V`, `vLoop_*` / `setOuterPair_*` frame lemmas | `PlanarEmbedV.lean`, `PlanarEmbedVLoop.lean` | **proved** (standard axioms); all fields of the strengthened `GluedUpTo` are preserved |
| `pathEdgeRot_planar`, `badQState_before`, `badQState_after`, `badQState_excluded`; `GluedUpTo.outer_vertex` | `Proofs/PlanarEmbedQCounterexample.lean`, `PlanarEmbedSteps.lean` | **proved** (standard axioms) / boundary contract; initialization, leaf, F, V preserve the strengthened invariant |
| `badCapState_before`, `badCapState_after`, `badCapState_excluded`; `GluedUpTo.outer_cap` | `Proofs/PlanarEmbedCapCounterexample.lean`, `PlanarEmbedSteps.lean` | **proved** (standard axioms) / cap endpoint and presence contract; initialization, leaf, F, V preserve it |
| `maximal_not_inSubtree`, `maximal_pieces_disjoint`, `v_children_nonempty`, `v_parent_not_v`, `v_pieces_meet` and tree/localization helpers | `PlanarEmbedTree.lean`, `Proofs/PieceLoc.lean` | **proved** (standard axioms); disjointness of maximal pieces follows from the preorder child partition, and V hypotheses give the shared attachment vertex |
| `Piece.loc_data`, `mem_of_loc`, `loc_exists`, `loc_lt`, `loc_injective`, `loc_append_left/right`, `agrees_append` | `Proofs/PieceLoc.lean`, `Proofs/PieceAppend.lean` | **proved** |
| `mem_edgesBelow_data`, `mem_edgesBelow_lt`, `edgeIn_of_mem_edgesBelow`, `edgesBelow_nodup` | `PlanarEmbedEdges.lean` | **proved** |
| `closeOuter_spec`: close a single exposed facing pair, preserve its embedding, and frame all other edges | `PlanarEmbedClose.lean` | **proved** |
| `child_data`, `children_nodup`, `parent_some_iff`, `child_maximal_one`, `maximal_zero_eq`, `mem_pieceBelow_bound`, `touches_of_hasEdge`, `root_pieces_disjoint` | `PlanarEmbedEdges.lean` | **proved** (standard axioms) |
| `closeList_nil/cons`, `closeList_outerE`, `closeList_rotAdj_size`, `closeOuter_get_congr`, `closeList_frame/get`, `forM_closeOuter_run`, `embedItem_F` | `PlanarEmbedCloseList.lean` | **proved** (`propext`, `Quot.sound`) |
| `Piece.union_embeddings` | `Proofs/PieceUnion.lean` | **proved** (standard axioms) |
| `root_child_exposed_mem`, `embedItem_step_F` | `PlanarEmbedF.lean` | **proved** (standard axioms); F admission removed |
| `starEdgeRot_planar`, `badVState_before`, `badVState_after`, `badVState_excluded` | `Proofs/PlanarEmbedVCounterexample.lean` | **proved** (standard axioms); attachment-vertex counterexample and exclusion |
| planar walk (`planarWalk`), relabel (`planarSpqrTree`), gluing (`planarEmbed`) | `PlanarWalk.lean`, `PlanarRelabel.lean`, `PlanarEmbed.lean` | def |
| `planarWalk_base`, `planarWalk_proj` (planar walk = ordinary walk + aux) | `PlanarWalkProj.lean` | **proved** (`propext`, `Quot.sound`) |
| `planarRelabelTree_base`, `planarRelabel_proj` | `PlanarRelabelProj.lean` | **proved** (`propext`, `Quot.sound`) |
| `planarEmbed_isSome_iff` | `PlanarSpec.lean` | **proved** (`propext`, `Classical.choice`, `Quot.sound`) |
| orbit counting: `numOrbits` = orbit minima, `numOrbits_involution`; `numComponents` by label relaxation | `Planar.lean`, `Proofs/Planar.lean` | def / **proved** |
| `numOrbits_eq_orbitCount` (executable orbit minima count = `orbitCount` for a permutation of `range n`) | `Proofs/OrbitCount.lean` | **proved** (`propext`, `Classical.choice`, `Quot.sound`) |
| `numComponents_eq_ccCount` (label relaxation converges in `n` rounds to component minima; `numComponents` = number of `EdgesConn` classes with an edge) | `Proofs/CompCount.lean` | **proved** (standard axioms) |
| `ccCount_union`, `ccCount_ident_of_conn` / `_of_not_conn` (component count under disjoint union and identification of two non-isolated vertices) | `Proofs/CompGlue.lean` | **proved** (standard axioms) |
| `TwoSum.numNonIsolated_edges` (`V + 2 = V₁ + V₂`), `TwoSum.numComponents_edges` (`C + 1 = C₁ + C₂`, via `WF.conn`: the 2-sum is `G₁ ⊔ (G₂ + n₁)` with `u`, `v` identified) | `Proofs/TwoSumCount.lean`, `Proofs/TwoSumComp.lean` | **proved** (standard axioms) |
| `TwoSum.splice_isPlanarEmbedding`, `TwoSum.planar` (explicit splice of two planar embeddings along the virtual edge, `f = f₁ + f₂ − 2`) | `PlanarGlue.lean`, `Proofs/PlanarGlue.lean` | **proved** (`propext`, `Classical.choice`, `Quot.sound`) |
| closed forms of `layoutRot .S` / `.P`; `cycleRot_isPlanarEmbedding`, `bondRot_isPlanarEmbedding` | `PlanarLayout.lean` | **proved** (`propext`, `Classical.choice`, `Quot.sound`) |
| `layoutRot_S_shift`, `layoutRot_P_shift` (renumbering a node's edges from 0) | `PlanarInv.lean` | **proved** |
| Invariant P (§8.2) as a Lean structure: `Piece`, `SideWalk`, `InvariantP`, `StackInv` | `PlanarInv.lean` | def (the deliverable is the statement) |
| `planarWalkOut_stackInv` (the walk preserves Invariant P) | `PlanarInv.lean` | sorry |
| `RotInv` (`neRotAdj` = concatenation of one `layoutRot` block per numbered node, in lockstep with `neBounds`/`nvBounds`), planar-state `wp`/`Frame` calculus, `planarRelabel_rotInv` (the `planarRelabel` fold preserves `RotInv`; one abstract `Frame` per primitive, `match`/`if` arms framed by `wp_jp_match1`/`wp_jp_ite`/`wp_jp_ite3`) | `PlanarRotInv.lean`, `PlanarRotFold.lean` | **proved** (`propext`, `Classical.choice`, `Quot.sound`) |
| `planarRelabel_rot_spec` (`neRotAdj` has four entries per node-edge; node `n`'s entries `4 neSt + j` are those of `layoutRot (type n) (nVerts n) neSt neEn edgeVes mapRot (2 g.ne)`, `mapRot` four entries each, `|edgeVes| + 1 = nEdges n` for R) | `PlanarRotSpec.lean` | **proved** from `planarRelabel_rotInv` (`rot_spec_of_inv`); the size clause uses `layoutRot_size`, hence `WF.shape` via `spqrTree_wf` (admitted in `relabelTree_wf`) |
| `layoutRot_size` (every type's `layoutRot` has `4 · nEdges` entries, from `WF.shape`/`Twins.cap_none`), `neRotAdj_segment` (relabel bookkeeping: node `i`'s `neRotAdj` segment is its `layoutRot`) | `PlanarRotSpec.lean`, `PlanarSpec.lean` | **proved** modulo `spqrTree_wf` (through `planarRelabel_proj`) |
| `nodePlanar_sound_S`, `nodePlanar_sound_P` | `PlanarSpec.lean` | proved modulo `Shape` (via `spqrTree_wf`, itself admitted in `relabelTree_wf`) |
| `PlanarFinish` (walk-side record, kept apart from `WalkInv`: a node recorded `.planar m` has `m` = four exposed ends of its edge children (`match_ends`), children in range (`ch_lt`), the `qem` after `capLinked` (= `setupNode`'s cap links) and `flipped` (= `applyFlips`) closed on the node's virtual edges (`closed`), and `readRot` of it is a planar embedding of `nodeSkel` (`embedded`)), `planarWalk_planarFinish` | `Proofs/PlanarWalkFacts.lean` | def / sorry (checked as separate fields by `checks/WalkInvCheck/Planar.lean`, seeds 0..3000 both modes, 0 violations; the variant without the flips fails at seed 6) |
| `RelabelNodeR` (relabel-side record for an R node `n`: its item `it`, `children = Items.ordered`, `PosOK`, `nvEn = nvSt + |nvList|`, `neEn = neSt + |edgeVes| + 1`, skeleton `(nvSt, nvEn - 1) :: edgeChildren`, rows `neRotAdj[4 neSt + l] = some (4 (neSt + idxOf) + (o &&& 2) + (1 - l % 2))` read off `Q = flipped (capLinked qem m)`), `planarRelabelTree_relabelNodeR` | `PlanarRelabelNodeR.lean` | def / sorry (the `mapRot`/`setupNode → nodeRot` bridge, via `planarRelabel_rotInv`/`relabel_node_spec`; in progress) |
| `readRot_perm` (`readRot` along a permutation of `ves`, by `IsPlanarEmbedding.reindex`), `nodePlanar_sound_R_of` (`PlanarFinish` + `RelabelNodeR` + `Items.WF` ⇒ `nodeRot n` is a planar embedding of `localSkeleton n`: reorder to `Items.ordered`, relabel vertices by `pos - nvSt` with `IsPlanarEmbedding.map`, injective by `PosOK`/`nv_nodup`), `planarWalk_items_wf'` (from `walk_items_wf`), `planarWalk_items_wf` (hypothesis-free form) | `PlanarSoundR.lean` | **proved** (`propext`, `Classical.choice`, `Quot.sound`) / proved / sorry (cf. `spqrTree_wf`) |
| `nodePlanar_sound_R` | `PlanarSpec.lean` | proved from `nodePlanar_sound_R_of` modulo the three admissions above |
| `nodePlanar_sound` = S ∨ P ∨ R cases | `PlanarSpec.lean` | proved from the three |
| `nodePlanar_complete` (Kuratowski-style certificate from the §8.3 crossing) | `PlanarSpec.lean` | sorry, hard |
| `IsPlanarEmbedding.map`, `Planar.map` (relabelling the vertices by a map injective on the non-isolated ones; `numNonIsolated`/`ccCount` are transported) | `Proofs/PlanarMap.lean` | **proved** (`propext`, `Classical.choice`, `Quot.sound`) |
| `RotationSystem.union` (second system's quarter-edges shifted by `rs₁.size`), `union_stepC` = `unionStep`, `union_numFaceOrbits`/`union_numVertexOrbits` (orbits add), `IsPlanarEmbedding.union`, `Planar.union` (`es₁ ++ shiftEdges n₁ es₂` on `n₁ + n₂` vertices) | `Proofs/PlanarUnion.lean` | **proved** (`propext`, `Classical.choice`, `Quot.sound`) |
| `RotationSystem.conj a b` (rotation conjugated by the transposition of two quarter-edges `a`, `b` of the same direction: the rotations at their vertices are spliced into one), `conj_stepC` = two `swapImg`s, `conj_numFaceOrbits`/`conj_numVertexOrbits` (`+ 2 =`, when `a`, `b` lie in different `stepC`-invariant halves), `IsPlanarEmbedding.oneSum` (identify `v₁` and `n₁ + v₂`, both with an edge, on `rs₁.union rs₂`: `F + 2 = F₁ + F₂`, `V + 2 = V₁ + V₂`, `C + 1 = C₁ + C₂`) | `Proofs/PlanarOneSum.lean` | **proved** (`propext`, `Classical.choice`, `Quot.sound`) |
| `disjointUnion_planar` (= `Planar.union` after `disjointUnionEdges_eq`), `oneSum_planar` (`oneSumEdges_eq`: `oneSumEdges` = `collapse (n₁ + v₂) ∘ ident v₁ (n₁ + v₂)` on the union; `IsPlanarEmbedding.oneSum` when both vertices have an edge, otherwise `Planar.map` alone), `twoSum_planar` (`twoSumEdges_eq`: `twoSumEdges` = `collapse₂ (n₁ + u₂) (n₁ + v₂)` on `TwoSum.edges`; `Planar.map` of `TwoSum.planar`) | `PlanarInv.lean` | **proved** (`propext`, `Classical.choice`, `Quot.sound`). `twoSum_planar` was **restated**: its former hypotheses (`he₁`, `he₂`, `u₁ ≠ v₁`, `u₂ ≠ v₂`) do not give `TwoSum.WF`'s `deg₁`/`deg₂`/`face`/`conn`, which the explicit splice needs (with a virtual edge whose end has degree 1 the splice rewires the edge to itself; with a bridge on both sides the 2-sum is two 1-sums, not a splice). It now takes `hdeg₁`, `hdeg₂` (`rs.get (4 e + k)` is never a quarter-edge of `e`), `hface` (one virtual edge separates two face orbits) and `hconn` (one virtual edge is not a bridge of its side) verbatim from `TwoSum.WF`. The unconditional statement is true (degenerate cases are 1-sums / relabellings of the pieces) but not needed: at the call site `embedItem_step_node` the skeleton side is a cycle (`S`, ≥ 3 edges), a bond (`P`, ≥ 3 edges) or a 3-connected `R` skeleton, so `hconn` holds on the skeleton side and `hdeg₁`/`hdeg₂` hold because every skeleton/piece vertex has degree ≥ 2 at the virtual edge; `hface` on the skeleton side follows for `S`/`P` from the closed forms of `layoutRot` (`PlanarLayout.lean`) and for `R` needs the face structure of `nodeRot` — supplying these in `embedItem_step_node` is the open part of step 3 |
| `WalkInv` / `Preserves` (Hoare triple on `PlanarWalkM`), `Preserves.frame`, `Preserves.popPair`, `planarFinishEdge_inv`, `planarWalkOut_inv` (+ `Tree`/`Outs`), `planarWalkOut_stackInv` | `PlanarInvSteps.lean` | def / **proved** modulo the per-step lemmas |
| `pushVertTstack_inv` (a fresh vertex entry has no exposed ends, so `StackInv` exempts it) | `PlanarInvSteps.lean` | **proved** |
| per-step lemmas `pushEdgeTstack_inv`, `mergeTstackTops_inv`, `maybeUnwrapNxt_inv`, `finishTstackTop_inv`, `closeBackedges_inv`, `flipBeforeMerge_inv`, `pruneBackedges_inv`, `flipForLowval_inv`, `modifyCur_foldSides_inv` (one `planarFinishEdge` step each preserves Invariant P) | `PlanarInvSteps.lean` | sorry |
| gluing invariant `GluedUpTo` (per maximal processed item: restricted rotation agrees with a planar embedding of the piece below it, exposed ends = unset ends, on one face; `outer_unprocessed`: unprocessed items `j < i` have no exposed ends), `gluedUpTo_init`, fold `forM_reverse_range_inv` | `PlanarEmbedSteps.lean` | def / **proved** (`propext`, `Classical.choice`, `Quot.sound`) |
| `SpqrTree.ChildShape` (children shapes `embedItem` relies on but `SpqrTree.WF` does not imply — `nv_layout` allows a `V` child under an `O`/`I` leaf; currently `leaf`: `O`/`I` items are childless), `RelabelOK.childShape` (from `Items.Shapes.i_o_leaf` through `RelabelOK.mem_children_iff`), `relabelTree_childShape`, `spqrTree_childShape` | `PlanarShape.lean`, `RelabelChildShape.lean` | def / **proved** (`propext`, `Classical.choice`, `Quot.sound`); `spqrTree_childShape` inherits `walk_items_wf`'s `sorry` like `spqrTree_wf` |
| `embedItem_step_Q` / `_node` (`GluedUpTo`-level Q / node steps) | `PlanarEmbedSteps.lean` | removed: superseded by `embedItem_step_Q_faces` / `embedItem_step_node_faces` over `GluedFaces` (the `GluedUpTo` statement of the Q step is refuted by `badFaceState_*`) |
| `EmbedM.link`/`outer`/`setOuter` frame lemmas on an abstract `EmbedState` (`link_some_run`, `link_outerE`, `link_rotAdj_size`, `link_rotAdj_get`, `link_exposedAt`, `setOuter_run`, `setOuter_rotAdj`, `setOuter_outerE_*`, `setOuter_exposedAt_ne`) | `PlanarEmbedLink.lean` | **proved** (`propext`, `Classical.choice`, `Quot.sound`) |
| `ChildShape.tile` (the `k`-th child of `i` is `i + 1 +` the subtree sizes of the earlier children; from `RelabelAll.child_idx_eq`/`children_sum`), `children_eq` (`toSpqrTree.children` = `PlanarSpqrTree.children`), `range'_children`, `edgesBelow_eq` (`edgesBelow i` = the item's own `Q` edge followed by the children's `edgesBelow`, under `Preorder.subtree_eq` + `ChildShape.tile`) | `PlanarShape.lean`, `RelabelChildShape.lean`, `PlanarEmbedBelow.lean` | **proved** (`propext`, `Classical.choice`, `Quot.sound`; `spqrTree_childShape` inherits `walk_items_wf`'s `sorryAx` like `spqrTree_wf'`) |
| `embedItem_step_leaf_of` (the `O`/`I` step from `GluedUpTo (i+1)` + `subtreeEnd[i] = i+1`), `embedItem_leaf` (`embedItem` is the identity on leaves), `edgesBelow_leaf`, `isPlanarEmbedding_nil`; `subtreeEnd_leaf` (`WF.preorder.subtree_eq` + `ChildShape.leaf`), `embedItem_step_leaf`, dispatcher `embedItem_step`, `gluedUpTo_planarEmbed` | `PlanarEmbedLeaf.lean`, `PlanarEmbedFold.lean` | **proved** (`propext`, `Classical.choice`, `Quot.sound`); `embedItem_step`/`gluedUpTo_planarEmbed` modulo the two admitted steps Q/node. The former unconditional `embedItem_step_*` (no `WF`) were false for an arbitrary `PlanarSpqrTree` (an `O` item with a `Q` item in its `subtreeEnd` range), and `WF` alone is still too weak for the leaf step: `WF` admits an `O` item `2` with a `V` child `3` carrying a block-root `Q` child `4` (`nv_layout` with `a = b = []`, `twin_parent` with the non-node parent `V`); after items `4`, `3` the two ends of `4`'s edge exposed at `3` are `some none` in `rotAdj` but not `exposedAt 2`, so `GluedUpTo 2` fails although `embedItem 2` is the identity — hence `ChildShape.leaf` |
| `glued_root` | `PlanarSpec.lean` | **proved**, `sorryAx` only through accepted `walk_items_wf` via `spqrTree_childShape` |
| `glued_root_of`, `edgesBelow_subset_root`, `mem_edgesBelow_root_iff`, `edgesBelow_root_perm`, `pieceBelow_root_perm` | `PlanarEmbedRoot.lean` | **proved** (standard axioms) |
| `IsPlanarEmbedding.reindex`, `graphCounts_congr`, `Piece.loc_xor` | `Proofs/PlanarReindex.lean`, `Proofs/PieceLoc.lean` | **proved** (standard axioms) |
| `doubleRot_planar`, `badFaceState_before`, `badFaceState_after`, `doubleRot_not_cofacial` | `Proofs/PlanarEmbedFaceCounterexample.lean` | **proved** (standard axioms); the cofacial correction is `GluedFaces.cap_face` |
| `GluedFaces`, `PieceWitness`, `CapFace`, `gluedFaces_init`, `pieceWitness_frame`, `gluedFaces_of_frame`, `embedItem_step_leaf_faces`, `embedItem_step_F_faces`, `embedItem_step_V_faces` | `PlanarEmbedFaces.lean`, `PlanarEmbedFacesSteps.lean` | **proved** (standard axioms) |
| `orbitCount_swapImg_same`, `RotationSystem.conj_numFaceOrbits_split` | `Proofs/OrbitSplit.lean` | **proved** (standard axioms) |
| `RotationSystem.single`, `RotationSystem.loop`, `RotationSystem.insert`, `IsPlanarEmbedding.insert` | `Proofs/PlanarInsert.lean` | **proved** (standard axioms) |
| `Piece.Capped`, `Piece.Open2`, `Capped.insert`, `Open2.close`, `Open2.splice`, `single_open2` | `Proofs/PieceInsert.lean` | **proved** (standard axioms) |
| `Capped.join` (1-sum of two capped pieces at a shared cap endpoint, node-step link order) | `Proofs/PieceJoin.lean` | **proved** (standard axioms) |
| `Capped.attachOpen` (open `V` piece hung at a capped end by the corner link), `capped_child` (child certificate from `GluedFaces`) | `Proofs/PieceJoin.lean`, `PlanarEmbedNode.lean` | **proved** (standard axioms) |
| `IsPlanarEmbedding.ident_conj` (identify two connected vertices by a `conj` of two cofacial pairs), `IsPlanarEmbedding.splice2` (2-sum on the common vertex space by two `conj`s), `iterate_eq_of_eqOn`, `IsPermOn.exists_first_hit`, `RotationSystem.sameOrbit_stepC3_rot`, `Capped.parJoin` (P-node step: two capped pieces with the same cap endpoints, node-step link order, cofaciality of `c0`/`d2` kept) | `Proofs/PlanarIdent.lean`, `Proofs/PlanarSplice2.lean`, `Proofs/OrbitSegment.lean`, `Proofs/PieceParJoin.lean` | **proved** (standard axioms) |
| `RotationSystem.genus_le` (genus inequality `F + 2V ≤ 2(2C + E)` for every `IsEmbedding`, by first-edge deletion: `DetachAt` conjugations, `detach_one`/`detach_both`, `eq_union_shrink`/`isEmbedding_of_union`, `genus_cons_bound`), `edgesConn_of_not_cofacial` (non-cofacial sides ⇒ endpoints connected in the tail), `numNonIsolated_cons_eq`, `numComponents_cons_le_succ`/`_succ_le`/`_le_of_isolated`/`_le_of_loop` | `Proofs/PlanarGenus.lean`, `Proofs/PlanarCompCons.lean` | **proved** (standard axioms) |
| `IsPlanarEmbedding.uninsert` (delete a non-bridge first edge of a planar embedding: sides distinct faces, neighbours above the edge, planar `σ` with explicit rotation and cofacial `rot 1 - 4`/`rot 3 - 4`), `IsPlanarEmbedding.cap_sides_distinct`, `sameOrbit_detach_merge`, `detach_both` (extended with the rotation and cofaciality), `numComponents_cons_le_of_conn` | `Proofs/PlanarUninsert.lean`, `Proofs/PlanarGenus.lean`, `Proofs/PlanarCompCons.lean` | **proved** (standard axioms) |
| `embedItem_Q_nil`, `embedItem_Q_O`, `embedItem_Q_cons`, `qUpper`, `qLower`, `setOuter4` | `PlanarEmbedQExec.lean` | **proved** (exact unfolding of the executable `Q` branch) |
| `q_edge_not_below`, `q_fresh`, `q_lower_boundary`, `q_capped_child`, `q_open_glued`, `q_leaf_glued`, `q_loop_glued`, `q_block_glued`, `embedItem_step_Q_faces` | `PlanarEmbedQ.lean` | **proved** (`propext`, `Classical.choice`, `Quot.sound`) |
| `treeQe`, `nodeStep`, `embedItem_node`, `forIn_state_foldl`, `nodeStep_outerE_ne`/`_outerE_size`/`_rotAdj_size`/`_exposedAt_ne`, `nodeFold_*` | `PlanarEmbedNodeExec.lean` | **proved** (exact unfolding of the executable `S`/`P`/`R` branch + frame lemmas) |
| `node_capped_glued`, `capped_child` | `PlanarEmbedNode.lean` | **proved** (standard axioms) |
| S-node data layer (`shape_S`, `nEdges_S`, `capNe_S`, `nvs_S`, `noncap_S`, `neRotAdj_S`, `nodeStep_skip`/`_cap`/`_link`/`_corner`, `child_S`), chain layer (`sC`, `sW`, `sL`, `chain`, `child_cert_S`, `sW_boundary`, `nvInc_*`, `sep_*`, `disjoint_*`), fold (`blockFold`, `InvS`, `cap_block_S`, `invS_zero`, `invS_step`, `invS_all`, `last_block_S`), `capEnd_S`, `children_perm_sL`, `chain_perm_edgesBelow`, `nodeFold_capped_S` | `PlanarEmbedNodeS.lean`, `PlanarEmbedNodeSChain.lean`, `PlanarEmbedNodeSFold.lean`, `PlanarEmbedNodeSMain.lean` | **proved** (standard axioms) |
| P-node data layer (`shape_P`, `hasCap_P`, `capNe_P`, `nvs_P`, `noncap_P`, `neRotAdj_P`, `rotP_P`, `sC_length_P`, `sC_getElem?_P`, `child_P`, `nvOrig_P`, `neOrig_P`, `capOrig_P`, `sV01_ne_P`, `child_cert_P`, `noV_P`, `sC_eq_children_P`, `attach_P`), fold (`pL`, `pChain`, `InvP`, `cap_block_P`, `invP_zero`, `invP_step`, `invP_all`, `last_block_P`, `fold_eq_blockFold_P`, `pChain_eq_pieceBelow`), `nodeFold_capped_P` | `PlanarEmbedNodeP.lean`, `PlanarEmbedNodePFold.lean` | **proved** (standard axioms) |
| `LayoutAt`, `nodeFold_capped` (dispatch), `nodeFold_capped_R`, `embedItem_step_node_faces` | `PlanarEmbedNodeFold.lean` | **proved** (standard axioms) |
| `IsPlanarEmbedding.subst`, `substSum`, `edgesConn_subst`, `deg_of_edgesConn`; `ThreeConnected.edgesConn_eraseIdx`; `vstep`; `subst_step`, `sh1`, `pick4`; `RG`, `InvV`, `InvC`, `csys`, `invV_all`, `invC_all`, `RG.Hyp.glue`, `cap_conn_children` | `Proofs/PlanarSubst.lean`, `Proofs/ThreeConn.lean`, `Proofs/PlanarVStep.lean`, `Proofs/PlanarSubstStep.lean`, `Proofs/PlanarRGlue.lean` | **proved** (standard axioms) |
| R-node data layer (`shape_R`, `hasCap_R`, `capNe_R`, `neEn_R`, `noncap_R`, `sC_length_R`, `sC_getElem?_R`, `child_R`, `nvs_R`, `nvsOf_R`, `localSkeleton_length`, `localSkeleton_getElem?`, `nodeRot_get`, `nodeRot_size_le`, `rotR`, `neOrig_R`) | `PlanarEmbedNodeR.lean` | **proved** (standard axioms) |
| `NodeRotClosed`; `planarTree_nodeRotClosed` | `PlanarNodeSpec.lean`, `PlanarSpec.lean` | def / sorry (checked: `check_piece_sep` `rot_closed`) |
| `RFold`, `InvR`, `partner`, `nodeFoldR`, `RFold.Q_inj`, `corner_iff` | `PlanarEmbedNodeRExec.lean` | **proved** (standard axioms) |
| `rfold_R`, `slot_R`, `child_cert_R`, `corner_V_R`, `corner_R`, `capNvs_R`, `capEnd_R`, `sW_R`, `vOf_R`, `capT_R`, `T_facts` | `PlanarEmbedNodeRCert.lean` | **proved** (standard axioms) |
| `rgR`, `rgR_m`, `rgR_Sk`, `rgR_skE`, `rgR_cap`, `rgR_cρ`, `rgR_cl`, `rgR_vρ`, `rgR_vw`, `hyp_R` | `PlanarEmbedNodeRGlue.lean` | **proved** (standard axioms) |
| `vesV`, `vesC`, `pieceR`, `pieceR_es`, `pieceR_perm`, `locV_sub`, `locC_sub`, `locQ_R`, `locW_R`, `vert_slot_R`, `dir_slot_R`, `memQ_R`, `memW_R`, `mem_pieceR_iff`, `GlueOut`, `RG.Hyp.glue'`, `tgtV_ta_full`, `tgtV_tb_full`, `agrees_V_R`, `agrees_C_R`, `unset_R`, `fold_eq_R`, `main_R` | `PlanarEmbedNodeRMain.lean` | **proved** (standard axioms) |
| `embedItem_step_faces`, `gluedFaces_planarEmbed` | `PlanarEmbedFacesFold.lean` | **proved** (standard axioms modulo the fold hypotheses) |
| `CornerAt`, `NodeCorners`; `planarTree_nodeCorners` | `PlanarNodeSpec.lean`, `PlanarSpec.lean` | def / sorry (hypothesis of the node fold, reached by `planarEmbed_sound`) |
| `conj_not_sameOrbit_split`, `insert_not_sameFaceOrbit`, `Proofs/PlanarDegree.lean` | `Proofs/OrbitSplit.lean`, `Proofs/PlanarInsert.lean`, `Proofs/PlanarDegree.lean` | **proved** (standard axioms) |
| `localSkeleton`, `nodeRot`, `isPlanar`, `isPlanar_of_all`; `RotInv.np_size`, `planarTree_nodePlanar_size` | `PlanarNodeSpec.lean`, `PlanarRotInv.lean`, `PlanarSpec.lean` | def / **proved** |
| `planarEmbed_sound` | `PlanarSpec.lean` | proved from `gluedFaces_planarEmbed` (with `WF`/`ChildShape` of `g.planarTree …` from `spqrTree_wf'`/`spqrTree_childShape` via `planarRelabel_proj`, hence the `g.WF`/`OrderOK` hypotheses) + `glued_root` |
| `spqrTree_planar` (`→` from `planarEmbed_sound`; `←` needs completeness + skeletons are minors of `g`) | `PlanarSpec.lean` | sorry |

Admitted, precisely (`#print axioms` reports `sorryAx` for each): `spqrTree_wf` (inherited by `planarRelabel_rot_spec`, `neRotAdj_segment`);
`planarWalk_planarFinish`, `planarRelabelTree_relabelNodeR`, `planarWalk_items_wf` (hence
`nodePlanar_sound_R`); the nine per-step lemmas `*_inv` of `PlanarInvSteps.lean` (hence
`planarWalkOut_stackInv`); `nodePlanar_complete`; `spqrTree_pieceSep`'s final-items hypothesis
`walk_pieceInv` (§4.6), `planarTree_nodeRotClosed`, `planarTree_nodeCorners` (hence `planarEmbed_sound`; the old `GluedUpTo`-level
`embedItem_step_Q`/`embedItem_step_node` and the unattempted `TwoSum.planar_left` converse have
been removed); `spqrTree_planar`. The S and P
cases of `nodePlanar_sound` are proved except for the `Shape` of the
skeleton (`spqrTree_wf`, which inherits `relabelTree_wf`'s `sorry`); the local counting itself
(`cycleRot_isPlanarEmbedding`, `bondRot_isPlanarEmbedding`) is fully proved.

Work packages:
* **Invariant P**: `InvariantP` / `StackInv` (`PlanarInv.lean`) are stated; prove
  `planarWalkOut_stackInv` step by step for `makeEdgePlanarity`, `mergeSide`,
  `closeSide`/`pruneSide`, `foldPlanarity`, `finishMatches`/`unwrapPlanarity`; the merge-side
  crossing argument gives the `none` case.
* **Relabel transport**: done. `planarRelabel_rotInv` (`PlanarRotFold.lean`) is the induction over
  `planarRelabel`: `RotInv` (`PlanarRotInv.lean`) says `neRotAdj` is the concatenation of the
  `layoutRot` blocks of the nodes numbered so far; each call appends exactly one block while pushing
  `neBounds`/`nvBounds` (`RotInv.step`), and the recursive calls on the children preserve it. The
  proof follows `relabelRun_sizes_le`: a planar-state `wp`, one `Frame` lemma per primitive
  (`setupNode`, `applyFlips`, `modifyAux`, `liftR (modify _)`), and the `match`/`if` arms framed
  abstractly (`wp_jp_match1`, `wp_jp_ite`, `wp_jp_ite3`) so the continuation is elaborated once.
  `planarRelabel_rot_spec` is `rot_spec_of_inv` at the final state; only its `neRotAdj.size`
  clause and the per-node `j`-range go through `layoutRot_size`, i.e. `WF.shape` (`spqrTree_wf`).
  `neRotAdj_segment` is derived from it (`PlanarRotSpec.lean`).
* **Local embeddings**: S and P done; R reduced (`nodePlanar_sound_R_of`) to the walk-side
  record `PlanarFinish` (checked; its proof is `finishTstackTop` + Invariant P, i.e. the nine
  `*_inv` lemmas) and the relabel-side record `RelabelNodeR` (to be proved from
  `planarRelabel_rotInv` / `relabel_node_spec`).
* **Gluing**: `twoSum_planar` (explicit splice, `PlanarGlue`, under `TwoSum.WF`'s
  `deg`/`face`/`conn`), `oneSum_planar` (`RotationSystem.conj` on `RotationSystem.union`),
  `disjointUnion_planar` (`RotationSystem.union`) are proved; the bottom-up induction over
  `embedItem` is done (`forM_reverse_range_inv`), what remains are the per-item steps
  `embedItem_step_*` (each one `link` = one splice of two pieces' exposed ends; the node step must
  supply `twoSum_planar`'s `hdeg`/`hface`/`hconn` from the skeleton's shape).
* **Completeness**: the crossing of §8.3 as a `K₅`/`K₃,₃` subdivision — the hard, optional one.
