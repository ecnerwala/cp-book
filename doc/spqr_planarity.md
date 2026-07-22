# Per-node planarity testing + skeleton embedding in `spqr_tree`

This note explains how `spqr.hpp` computes, for every SPQR node as it is
sealed, a planarity verdict and a combinatorial embedding of that node's
virtual-edge skeleton — what the SPQR machinery already provides, what had to
be added, and why the split is exactly the way it is.

## What the SPQR construction already provides

The construction shares its preprocessing with left-right (LR) planarity
(de Fraysseix–Rosenstiehl; Brandes' writeup is the reference used here):

- `dfs_lowval` is LR's phase 1 (DFS orientation): it computes `depth`,
  per-edge `lowpt`/`lowpt2`, and orients every edge (tree edges parent→child,
  back edges deep→shallow).
- The bucket key `1 + 2*lowpt + (lowpt2 < depth)` used by `build_sorted_adj`
  is exactly LR's nesting depth, so `dfs_spqr` traverses in LR's required
  order, and the estack entries of a sealed node's range appear in DFS order
  of the node's skeleton: every tree vedge appears *after* all entries of its
  subtree, chords appear at their (deeper) source.
- Crucially, every estack entry / vedge is a *collapsed 2-attachment
  component*: its interlacement behavior within its node's skeleton depends
  only on its two endpoints. So each vedge acts as a single tree or back edge
  of the skeleton, and the LR test factors completely over skeletons — the
  conflict component of a piece of graph closes exactly when its attachment
  set shrinks to a split pair, i.e. exactly when `pop_estack_range` seals it
  into a node.

## What the storage order does *not* provide

The one thing the sealed range's order does **not** preserve is the *nesting
order between siblings*: collapsing a subtree into a single virtual tree edge
can change its lowpoints and its "has a second attachment" (chordality) bit
relative to the original graph, so the original bucket-sorted sibling order
is no longer a valid LR nesting order for the skeleton. (This is observable:
random tests produce sealed R ranges where a skeleton-chordal virtual tree
edge precedes a chord with strictly smaller skeleton nesting depth.) Hence
`emit_rigid` recomputes skeleton lowpoints (one linear scan, using the DFS
postorder that the storage order *does* guarantee) and re-sorts each skeleton
vertex's outgoing vedges by skeleton nesting depth before running the test.

## Per-node emission (`finalize_node` → `emit_node_embedding`)

Every node type gets its skeleton rotation the moment it is created:

- **Q/I**: one edge; the identity rotation.
- **O/S**: the vedge range is already in cycle order (the estack path entries
  plus the cap), so the rotation is read off directly (`emit_cycle`).
- **P**: parallel edges don't interlace; storage order is used
  (`emit_parallel`) — this is where the `(k-1)!/2` embedding freedom lives.
- **R** (`emit_rigid`): the real work, in four small passes over the range:
  1. index the skeleton vertices;
  2. compute skeleton lowpoints bottom-up and re-sort each vertex's outgoing
     vedges by skeleton nesting depth;
  3. run the LR test on the skeleton (Brandes' algorithms 3–5: chords push
     singleton conflict pairs, finished tree edges merge their subtree's
     constraints into the parent's, backtracking trims return edges), with
     the two-coloring kept implicit via `lr_ref`/`lr_side`;
  4. resolve the signs and emit the rotation (Brandes' embedding phase:
     left edges reversed, then right edges; incoming chords are spliced next
     to the tree edge leading to their subtree).

A non-planar R skeleton only clears its own `node_t::planar` flag (its
rotation falls back to a structurally-valid but meaningless grouping);
`block_t::planar` and `is_planar` are the conjunctions. Since each skeleton
has O(its size) work and skeleton sizes sum to O(NVE), the whole thing adds
linear time up to the per-vertex re-sort (a `std::stable_sort` of each
skeleton's events).

## Output format

The output is a rotation system over vedge half-edges, node by node:
half-edge `2*ve+k` is the endpoint of `vedges[ve]` at `vedges[ve].vs[k]`, and
`embed_next[h]` is the next half-edge around the same skeleton vertex
(circular within the node). Rotations of different nodes can be glued along
twin virtual edges (`vedge_t::o_ve`), flipping the child's reflection as
desired — R rotations are unique up to reflection, P rotations may be
permuted, S/Q/I/O are rigid.

## Correctness checks

- `spqr.test.cpp` verifies on ~5000 random multigraphs (plus fixtures: K5,
  K3,3, Petersen, cube, grids, multigraphs with self loops/bridges, and a
  K5+K4 sharing an edge — one block, two R nodes, only one non-planar) that
  every node's rotation is a circular permutation of exactly its incident
  half-edges, that non-R nodes are always planar, and that every planar
  node's skeleton satisfies Euler's formula `V - E + F = 2` — which holds iff
  the rotation is a genus-0 (planar) embedding of the skeleton.
- The same per-node checks pass on 4000 random multigraphs with up to 60
  vertices and 140 edges under ASan/UBSan.
- `is_planar` was cross-checked against networkx's planarity test on 60,000
  random multigraphs (0 mismatches).
