# Unique planar embeddings, triconnectivity, and the LR conflict graph

Goal: prove that (essentially) *a graph has a unique planar embedding iff it is
triconnected*, in a constructive/effective way, and then relate the SPQR tree
to the left–right (LR) planarity test precisely: **if the SPQR algorithm
outputs a single (R) node, is the LR conflict graph connected?** Short answer:
*not literally* (K4 is a counterexample), but a precise corrected statement
holds, and the correction is exactly what you want for producing the conflict
graph online inside the SPQR DFS.

Throughout, $G$ is a biconnected (multi)graph. "Embedding" means a
*combinatorial embedding on the sphere*: a rotation system (cyclic order of
edge-ends around each vertex) whose induced face count satisfies Euler's
formula $V - E + F = 2$. Two embeddings are *equivalent* if they are equal or
mirror images (reverse every rotation). "Unique embedding" means unique up to
this equivalence.

---

## Part I. The abstract theorem

### 1. Faces of an embedding, and the gluing lemma

For a rotation system, faces are the orbits of the *next-dart* permutation:
after entering $v$ along edge $e$, leave along the successor of $e$ in the
rotation at $v$. In a 2-connected plane graph every face boundary is a cycle
(standard; a cut vertex is produced otherwise).

**Lemma 1 (faces determine the embedding).** Let $G$ be 2-connected and let
$\mathcal{F}$ be a set of cycles such that some embedding of $G$ has exactly
$\mathcal{F}$ as its face boundaries. Then that embedding is unique up to
mirror.

*Proof (constructive).* Each edge $e$ lies on exactly two cycles of
$\mathcal{F}$ (the two faces flanking it). Fix a vertex $v$: the faces at $v$
and the edges at $v$ form a bipartite "corner" structure — each face
containing $v$ contributes corners (consecutive edge pairs). Since each edge
at $v$ is in exactly two faces, the corner structure is a disjoint union of
cycles alternating edge/face; a valid rotation at $v$ is exactly a single such
cycle, and it is determined up to reversal. Choosing the direction at one
vertex propagates across each edge (the two faces flanking $e$ must be on
opposite sides of $e$ at both endpoints), and $G$ connected forces a globally
consistent choice — two choices total, mirror images. $\blacksquare$

So uniqueness of the embedding reduces to uniqueness of the *face set*.

### 2. Whitney's theorem, via non-separating induced cycles

**Definition.** A cycle $C \subseteq G$ is *non-separating* if $G - V(C)$ is
connected, and *induced* if it has no chords.

**Lemma 2 (easy direction; only 2-connectivity needed).** If $C$ is an induced
non-separating cycle of a planar $G$ with $G \ne C$, then $C$ bounds a face in
*every* embedding of $G$.

*Proof.* In any embedding, $C$ is a closed curve with two open disks. All of
$G - V(C)$ is connected, so it lies in one disk. Chords don't exist ($C$
induced), so the other disk contains nothing: it is a face. $\blacksquare$

**Lemma 3 (hard direction; needs 3-connectivity).** If $G$ is 3-connected and
planar, then every face boundary of every embedding is an induced
non-separating cycle.

*Proof (constructive: a violation yields a 2-cut).* Let $C$ bound face $f$ in
some embedding.
- If $C$ has a chord $xy$, then $x, y$ split $C$ into two arcs $P_1, P_2$,
  each with an interior vertex (else $xy$ is parallel to an edge of $C$, and
  then $\{x,y\}$ still works as below in a multigraph, or the graph is a
  triangle-with-multiedge, excluded by 3-connectivity of simple graphs; for
  multigraphs run the same argument on the sphere). Any path from
  $\mathrm{int}(P_1)$ to $\mathrm{int}(P_2)$ avoiding $\{x, y\}$ must cross
  the closed curve $C' = xy + P_1$ — impossible since $f$ is a face on one
  side and the chord blocks the other. So $\{x,y\}$ is a 2-cut.
- If $G - V(C)$ is disconnected with components $H_1, H_2, \dots$, each $H_i$
  attaches to $C$ in a set $A_i \subseteq V(C)$; since $f$ is a face, all
  $H_i$ live in the other disk. Take the embedding restricted to
  $C \cup H_1$: $H_1$ occupies a region whose boundary meets $C$ in an arc
  between two attachment vertices $x, y$ that separate $H_1$'s attachments
  from the rest — a Jordan curve through $x, y$ around $H_1$ shows $\{x,y\}$
  is a 2-cut (3-connectivity forbids $|A_1| \le 2$ directly too). $\blacksquare$

**Theorem 4 (Whitney).** A 3-connected planar graph has a unique embedding.

*Proof.* By Lemmas 2 and 3, the face set of *any* embedding equals exactly
$\{$induced non-separating cycles of $G\}$ — an abstract, embedding-free set.
By Lemma 1 the embedding is determined up to mirror. $\blacksquare$

Note the effectivity: the canonical face set is *checkable* (each candidate
cycle: chordlessness and connectivity of the complement), and Lemma 1's proof
is an algorithm reconstructing the rotation from the face set.

### 3. The converse, and the exact classification

The literal converse "not 3-connected $\Rightarrow$ multiple embeddings" is
**false**: cycles ($S$), a triangle, and theta graphs (one $P$ node with 3
members) are 2- but not 3-connected and still have unique embeddings. The
correct statement is in terms of the SPQR tree. Recall the embedding count
(di Battista–Tamassia): for biconnected planar $G$ with SPQR tree $T$,

$$\#\{\text{rotation systems}\} \;=\; 2^{r} \cdot \prod_{P \text{ nodes}} (k_P - 1)!$$

where $r$ = number of R nodes and $k_P$ = number of virtual+real edges of a P
node. (S/Q/I/O nodes contribute factor 1.)

*Lower bound, constructively.* Every SPQR tree edge is a split pair
$\{u, v\}$. Given an embedding, the **flip** at $\{u,v\}$ — reverse the
rotation of every vertex strictly inside one split component and reverse the
positions of that component's edge-ends at $u$ and $v$ — is again an
embedding (mirror a sub-disk of the sphere: faces are preserved piecewise;
the two faces incident to the virtual edge's "slot" get re-glued). Flips at R
nodes change the embedding (the R skeleton has, by Whitney, exactly 2
rotation systems and the flip toggles them); at a P node with $k$ members,
the members can be permuted into any of the $(k-1)!$ cyclic orders around the
poles. These choices are independent across $T$ (they act on disjoint
skeletons), giving the count as a lower bound.

*Upper bound.* Induction on $|T|$: pick a leaf node $N$ with virtual edge to
the rest; any embedding of $G$ restricts to an embedding of the leaf's
skeleton (with the virtual edge realized as a path/region through the rest)
and an embedding of $G / N$ (the rest, with $N$ shrunk to its virtual edge).
The skeleton embedding count is: 1 for S/Q/I/O (a cycle or bond of $\le 2$ is
rigid up to mirror, and the mirror is absorbed by the flip already counted),
2 for R (Whitney), $(k-1)!$ for P. Multiply.

**Corollary 5 (exact uniqueness classification).** A biconnected planar
(multi)graph has a unique embedding iff
$2^{r}\prod_P (k_P-1)! \le 2$, i.e. iff its SPQR tree contains **at most one
"non-rigid" contribution**: either (a) no R and no P node with $k \ge 3$
(pure S/degenerate: cycles, etc.), or (b) exactly one R node and nothing
else non-rigid, or (c) exactly one P node with $k_P = 3$ and no R.
In particular for **simple 3-connected graphs** ($T$ = single R node) the
embedding is unique — this is the headline statement, with (a),(c) as the
degenerate exceptions ("unique iff triconnected, up to cycles and thetas").

---

## Part II. Making it effective: LR planarity and the conflict graph

### 4. The LR framework, self-contained

Fix a DFS of $G$ from root $r$; orient tree edges downward, back edges
upward. For a back edge $f$, $\mathrm{low}(f)$ is its upper endpoint's depth;
for a branch (an outgoing edge $b$ at a vertex $u$), $\mathrm{lowpt}(b)$ is
the minimum return depth in its subtree. Assume (as in `spqr.hpp`'s bucket
sort) branches at every vertex are **sorted by lowpt**.

A planar embedding of $G$, restricted to the moment the DFS sits at vertex
$u$ with active root path $\pi = r \to u$, places each back edge on the
**left or the right** of $\pi$. The LR characterization (de
Fraysseix–Rosenstiehl; Brandes' formulation) is:

**Fork Lemma 6.** Let $u$ be a vertex with branches $b_1, b_2$ (in DFS order,
so $l_1 := \mathrm{lowpt}(b_1) \le l_2 := \mathrm{lowpt}(b_2)$). In any planar
embedding:
1. *(same-constraints)* all return edges of $b_1$'s subtree with return depth
   $> l_2$ lie on one common side, and all return edges of $b_2$'s subtree
   with return depth $> l_1$ lie on one common side;
2. *(different-constraint)* those two groups lie on **opposite** sides,
   whenever both are nonempty.

*Proof.* When $b_2$'s subtree is drawn, the curve
$C = (\text{tree path } l_1 \leadsto u) + b_1\text{'s lowest return}$
is a closed cycle (fundamental cycle of $b_1$'s lowest back edge) and $b_2$'s
subtree sits inside one of the two faces flanking $\pi$ at $u$. Every return
edge of $b_2$ landing strictly above $l_1$ lands on $C$'s path part and must
stay in that same region — one side (this is 1 for $b_2$; symmetrically for
$b_1$ with the roles of the reference cycle exchanged, using $l_2$). For 2:
a $b_1$-return at depth $d_1 \in (l_2, u)$ and a $b_2$-return at
$d_2 \in (l_1, u)$ on the same side would force the two subtrees' curves to
cross (Jordan): $b_2$'s subtree hangs from $u$ and reaches down to $l_2$;
a same-side $b_1$-return into the open interval $(l_2, u)$ is separated from
$\pi$ by that hanging curve. $\blacksquare$

**Definition (conflict graph $\Gamma$).** Vertices = back edges of $G$.
Generate constraints by applying Lemma 6 at every fork and every branch pair:
SAME edges within each qualifying group, DIFF edges across groups. The
*conflict graph* is this constraint graph; a **component** is a connected
component of $\Gamma$ under SAME $\cup$ DIFF.

**Theorem 7 (LR test, statement).** $G$ is planar iff constraints are
satisfiable (each component 2-colorable: contract SAME, check DIFF
bipartite). Moreover every admissible assignment, fed to the deterministic
insertion rule (at each vertex: order branches by lowpt, place each back edge
on its assigned side, nest same-side edges by return depth), produces a valid
planar embedding.

The proof of sufficiency is exactly the correctness of the LR embedding
phase; we take it as given (Brandes, *The left-right planarity test*).

Each component of $\Gamma$ has exactly 2 colorings, so admissible assignments
number $2^{c}$, $c$ = number of components. Flipping **all** components
simultaneously mirrors the embedding. So one is tempted to conjecture:

> unique embedding $\iff$ $c = 1$ $\iff$ SPQR = single node?

### 5. The conjecture is false as stated: K4

Take $K_4$ on $\{1,2,3,4\}$, DFS path $1{-}2{-}3{-}4$, back edges
$(4{\to}1), (4{\to}2), (3{\to}1)$.

- $(4{\to}2)$ vs $(3{\to}1)$: spans $(2,4)$ and $(1,3)$ strictly interlace —
  DIFF. (Via Lemma 6 at fork $3$: branch $b_1$ = tree edge to 4 with
  $l_1 = 1$, branch $b_2 = (3{\to}1)$ with $l_2 = 1$; $b_1$'s returns
  $> l_2 = 1$ is $\{(4{\to}2)\}$… and the different-constraint fires.)
- $(4{\to}1)$: returns to depth $1 = $ global low. It is **never** in a
  qualifying group (its return depth is never $>$ any lowpt) and it strictly
  interlaces nothing (its span contains everything). **Isolated vertex of
  $\Gamma$.**

So $c = 2$ although $K_4$ is 3-connected (SPQR = single R node). Yet $K_4$
has exactly **2 planar rotation systems** (verified by brute force over all
$2^4$ rotation systems, checking Euler's formula) — a mirror pair, i.e. a
unique embedding. Indeed all $2^c = 4$ admissible assignments produce the
same embedding up to mirror: flipping the isolated component $\{(4{\to}1)\}$
alone yields the *same face set* $\{123, 124, 134, 234\}$ (check by tracing),
i.e. it acts as (mirror) $\circ$ (flip of the other component). The
assignment $\to$ embedding map is 2-to-1 here **before** quotienting by
mirror.

**Moral: components of $\Gamma$ overcount the embedding freedom.** Two
distinct corrections are needed, and both are meaningful for the SPQR
correspondence:

1. Some component flips are *embedding-trivial* (act as identity or as global
   mirror). These are "fake" freedoms.
2. P-node freedom ($(k-1)!$ orderings) is **not** a side choice at all: it
   lives in the *tie-breaking order* of parallel members with equal
   lowpt/return depth, which the LR insertion rule fixes arbitrarily.
   $\Gamma$ never sees it.

### 6. The corrected theorem

**Definition.** A component $K$ of $\Gamma$ is *essential* if flipping $K$
(keeping all other components fixed) changes the output embedding up to
mirror.

**Theorem 8 (trichotomy).** Let $G$ be biconnected and planar, $K$ a
component of $\Gamma$. Exactly one of:
(i) flipping $K$ is the identity on embeddings;
(ii) flipping $K$ equals the global mirror;
(iii) flipping $K$ realizes a flip at a split pair — and then a specific
2-cut of $G$ is extracted from $K$'s attachments.

*Proof.* Let $E(K)$ be $K$'s back edges together with the tree paths of their
fundamental cycles, $\mathrm{hi}(K)$ the deepest vertex dominating all
sources, and $\mathrm{lo}(K)$ the minimum return depth in $K$. **Span
lemma:** no back edge $f \notin K$ has a constraint-qualifying attachment
strictly inside the open interval
$(\mathrm{lo}(K), \mathrm{hi}(K))$ *on the part of the tree spanned by $K$* —
otherwise Lemma 6 (at the fork joining $f$'s branch with a $K$-branch) would
put $f$ into $K$. Consequently every connection from the region
$R(K)$ (vertices/edges whose escape routes run through $K$'s span) to the
rest of $G$ passes through the two *poles* $p_{\mathrm{lo}}$ (tree vertex at
depth $\mathrm{lo}(K)$) and $p_{\mathrm{hi}}$:

- If $R(K)$ is a proper region (contains a vertex or a second edge beyond a
  single back edge) and the poles are not the poles of the whole block:
  $\{p_{\mathrm{lo}}, p_{\mathrm{hi}}\}$ is a **2-cut**, and flipping $K$ is
  exactly the split-pair flip of §3 — case (iii).
- If $R(K)$ attaches only at the poles of the entire block (as
  $(4{\to}1)$ in $K_4$: return at the root, span covering the whole path):
  flipping $K$ re-chooses which of the two faces flanking the reference cycle
  receives $R(K)$; when everything else is on one side this is a mirror of
  the sphere — case (ii); when composed effects cancel, case (i). $\blacksquare$

**Corollary 9 (corrected conjecture).** If the SPQR algorithm outputs a
single node (in particular a single R node — $G$ triconnected), then **every
component of $\Gamma$ is inessential**: there is at most one "real"
component, and every other component consists of back edges returning to the
global low (depth of the block's root pole) with full span — pole-attached
singletons whose flip is absorbed by the mirror. Quantitatively: the LR
conflict graph of a triconnected graph is **connected after adding the
relation "$f \sim$ everything" for each back edge $f$ with
$\mathrm{low}(f) = $ root and no strict interlacer**. $K_4$ shows this
correction is necessary (one such singleton) and it is the *only* kind of
overcount: by Theorem 8, any component that is neither pole-attached nor
merged would witness a 2-cut, contradicting the single-node output.

Conversely (from §3's lower bound): every SPQR tree edge (2-cut) produces at
least one essential binary freedom, which by Theorem 8(iii) must show up as a
separate component of $\Gamma$ (or as a tie-order at a P node). So:

$$\{\text{essential components of } \Gamma\} \;\hookleftarrow\hspace{-1.2em}\rightarrow\; \{\text{R-side flips at SPQR tree edges}\}, \qquad \{\text{P freedoms}\} \leftrightarrow \{\text{ties in the sort}\}.$$

---

## Part III. Reading this off the SPQR algorithm (online conflict graph)

The correspondence maps directly onto the machinery in `src/spqr.hpp`:

1. **Sorted adjacency = Brandes' preprocessing.** `build_sorted_adj` buckets
   branches by lowpt (and by "has nontrivial lowval2"), which is exactly the
   branch ordering Lemma 6 assumes. The same DFS drives both algorithms.

2. **`tstack` candidates = component boundaries.** A live `tstack` entry is a
   candidate split pair $\{$`vstart`$,$ top of range$\}$ — precisely a
   potential pole pair $\{p_{\mathrm{lo}}, p_{\mathrm{hi}}\}$ of a
   yet-unmerged conflict component. While the candidate survives, the back
   edges attached strictly inside its interval form a union of $\Gamma$
   components not yet linked to the outside.

3. **`tstack` pops = constraint edges.** Each time the SPQR DFS pops tstack
   entries because a new edge's lowpt reaches below the candidate
   (`while (tstack.back().top_depth > t_ins.top_depth)` for back edges, and
   the pops in the type-1/`first_occurrence` branches for tree edges), that
   is exactly a Fork-Lemma constraint (a SAME/DIFF edge of $\Gamma$) merging
   the candidate's component with an outer one. **To build $\Gamma$ online,
   maintain a union–find over back edges keyed to estack ranges: union
   whenever a tstack pop / `pop_estack_range` merge occurs, and record the
   DIFF parity (one bit per union, from which side of the reference cycle
   triggered the pop).** Since every union event coincides with an existing
   $O(1)$-amortized stack operation, this stays linear (with UF's
   $\alpha$).

4. **`push_estack_p` (P nodes) = tie freedom, no $\Gamma$ edge.** Parallel
   members have equal spans; they generate no side constraints, matching
   §5(2): their $(k-1)!$ freedom is an ordering, not a coloring.

5. **The residual singletons are known in advance.** By Corollary 9 the only
   $\Gamma$-components not merged into the spine when a single R node is
   emitted are back edges into the current block's root edge poles (the
   `start_spqr(cur, nxt, e)` pair) with full span — the algorithm can detect
   them as the back edges whose insertion never pushes a surviving tstack
   candidate (`estack.back().vs == e_ins.vs` at the root, or return depth
   $=$ root with empty qualifying groups) and pre-merge them with the spine,
   after which: **single SPQR node output $\Rightarrow$ single conflict
   component**, on the nose.

### Summary of the quantification

- The naive statement "single SPQR component $\Rightarrow$ connected conflict
  graph" fails, with $K_4$ as the minimal counterexample; the failure mode is
  *only* pole-attached full-span singletons (return depth $=$ block root, no
  strict interlacer), whose flips are mirror/identity on embeddings
  (Theorem 8 (i)/(ii)).
- After pre-merging those (a local, online-detectable condition), the
  statement is true, and moreover *every* surviving component boundary
  certifies a 2-cut (Theorem 8 (iii)) — i.e. the conflict graph components
  and the SPQR tree edges are two views of the same object, computed by the
  same stack events.

### Loose ends / caveats

- Theorem 7's sufficiency direction (admissible assignment $\Rightarrow$
  planar) is imported from Brandes; everything else above is self-contained.
- Theorem 8's span lemma is stated for the constraint set generated by
  Lemma 6 exactly as formulated; if one uses a weaker constraint set (pure
  interlacement only, no fork same-constraints), more inessential components
  appear and the pre-merge rule must be extended accordingly.
- The counting map $\{$assignments$\} \to \{$embeddings$\}$ is $2^{(\#$
  inessential components$)}$-to-1 onto its image before the mirror quotient
  ($K_4$: 4 assignments $\to$ 2 rotation systems); ties at P nodes multiply
  the image independently.
