#pragma once

#include <iostream>
#include <algorithm>
#include <array>
#include <limits>
#include <vector>
#include <cassert>
#include <utility>

/** Block-Cut/SPQR Tree
 *
 *  A class to compute and represent the block-cut and SPQR trees of a static
 *  undirected (multi)graph.
 *
 *  We have several node types:
 *  - Q: a single real edge and a single virtual edge
 *  - I: a bridge edge of the original graph
 *  - O: exactly 1 self loop or 2 parallel edges (a degenerate case of both S and P)
 *  - S: at least 3 series edges (in a cycle)
 *  - P: at least 3 parallel edges
 *  - R: a rigid component
 *  Note that 2 parallel Q nodes are never glued directly together; they are always
 *  glued to an O or P node.
 *
 *  Every real edge is stored in a Q node. The Q node for edge e is numbered e
 *  and contains virtual edge e. You can test that virtual edge corresponds
 *  to a single real edge using `tree.vedges[ve].o_type == node_type::Q` or
 *  less directly `tree.vedges[ve].o_ve < NE`.
 *
 *  We have several named objects: Vertices, Components (CCs), Blocks (BCCs),
 *  Nodes (TCCs), and VirtualEdges (which include edges as described above).
 *  The various objects have the following relationships:
 *   * Components contain Vertices
 *   * Components contain Blocks contain Nodes contain VirtualEdges
 *   * Vertices are adjacent to Blocks, Nodes, and VirtualEdges (many-to-many)
 *
 *  The tree is specified unrooted, but is implicitly rooted: everything is
 *  numbered in postorder (except for Q nodes/their virtual edges). In
 *  particular, Nodes contain intervals of VirtualEdges, Blocks contain
 *  intervals of Nodes, and so on. Note that the Block/Component intervals do
 *  not include their Q Nodes/VEdges!
 *
 *  More specifically, vectors of adjacent/contained objects are all
 *  represented as range_t's, representing subintervals of the corresponding
 *  array. Note that, vertices do not have 1-to-1 relationships, so for those
 *  relations, we need to look up the corresponding array.
 */
struct spqr_tree {
	static constexpr bool PRINT_DEBUG = false;

	// Number of vertices
	int NV = 0;
	// Number of real edges
	int NE = 0;

	// Number of connected components
	int NC = 0;
	// Number of blocks
	int NB = 0;
	// Number of SPQR nodes
	int NN = 0;
	// Number of virtual edges
	int NVE = 0;

	// TODO: make range_t's easier to use, particularly the ones requiring a second indirection.
	struct range_t {
		int st = -1, en = -1;
		int size() const { return en - st; }
		bool contains(int v) const { return st <= v && v < en; }
	};
	struct iterable_range_t : public range_t {
		struct iterator {
			int v;
			int operator * () const { return v; }
			iterator& operator ++ () { ++v; return *this; }
			friend bool operator != (iterator a, iterator b) { return a.v != b.v; }
		};
		iterator begin() const { return iterator{st}; }
		iterator end() const { return iterator{en}; }
	};
	template <typename T> struct bound_array_range_t {
		typename std::vector<T>::const_iterator st, en;
		int size() const { return int(en - st); }
		auto begin() const { return st; }
		auto end() const { return en; }
	};
	template <typename T, std::vector<T> spqr_tree::* array> struct array_range_t : public range_t {
		bound_array_range_t<T> bind(const spqr_tree& tree) const {
			return {(tree.*array).begin() + st, (tree.*array).begin() + en};
		}
	};
	template <typename T, std::vector<T> spqr_tree::* array> struct enumerable_array_range_t : public iterable_range_t {
		bound_array_range_t<T> bind(const spqr_tree& tree) const {
			return {(tree.*array).begin() + st, (tree.*array).begin() + en};
		}
	};

	struct vertex_t;
	std::vector<vertex_t> vertices;
	std::vector<int> vertex_blocks;

	struct component_t;
	std::vector<component_t> components;
	std::vector<int> component_vertices;

	struct block_t;
	std::vector<block_t> blocks;
	std::vector<int> block_vertices;

	struct node_t;
	std::vector<node_t> nodes;
	std::vector<int> node_vertices;

	struct vedge_t;
	std::vector<vedge_t> vedges;

	struct vertex_t {
		int component = -1;
		array_range_t<int, &spqr_tree::vertex_blocks> vertex_blocks;
		// TODO: maybe we want vertex_nodes or even vertex_block_nodes or something
	};
	struct component_t {
		enumerable_array_range_t<block_t, &spqr_tree::blocks> blocks;
		enumerable_array_range_t<node_t, &spqr_tree::nodes> nodes;
		enumerable_array_range_t<vedge_t, &spqr_tree::vedges> vedges;
		array_range_t<int, &spqr_tree::component_vertices> component_vertices;
	};

	struct block_t {
		int component = -1;
		bool planar = true;
		enumerable_array_range_t<node_t, &spqr_tree::nodes> nodes;
		enumerable_array_range_t<vedge_t, &spqr_tree::vedges> vedges;
		array_range_t<int, &spqr_tree::block_vertices> block_vertices;
	};

	enum class node_type : char {
		Q = 'Q', I = 'I', O = 'O', S = 'S', P = 'P', R = 'R'
	};

	struct node_t {
		node_type type;
		int component = -1;
		int block = -1;
		bool planar = true;
		enumerable_array_range_t<vedge_t, &spqr_tree::vedges> vedges;
		array_range_t<int, &spqr_tree::node_vertices> node_vertices;
	};

	struct vedge_t {
		std::array<int, 2> vs;
		int component = -1;
		int block = -1;
		int node = -1;
		bool is_tree = false;

		// Information on the opposite virtual edge and its node
		int o_ve = -1;
		int o_node = -1;
		node_type o_type;
	};

	std::vector<int> depth;

	/** Planarity / combinatorial embedding
	 *
	 *  Planarity is computed per SPQR node as the skeletons are sealed
	 *  (`nodes[n].planar`); only R skeletons can be non-planar, and the graph
	 *  is planar iff every node is (`blocks[b].planar` and `is_planar` are the
	 *  conjunctions).
	 *
	 *  A combinatorial embedding is emitted per skeleton, over the virtual
	 *  edges: half-edge 2*ve+k is the endpoint of vedge ve at vertex
	 *  vedges[ve].vs[k], and embed_next[h] is the next half-edge around that
	 *  vertex within the same node's skeleton (a circular rotation). R
	 *  skeletons are embedded planarly (unique up to reflection) whenever
	 *  planar; S/P/O/I/Q rotations are the trivial ones. Whole-graph
	 *  embeddings can be obtained by gluing the skeleton rotations at twin
	 *  virtual edges (vedges[ve].o_ve), reflecting child skeletons as needed.
	 *  For a non-planar R node the rotation is an arbitrary (but still
	 *  circular) order.
	 */
	bool is_planar = true;
	std::vector<int> embed_next;

private:
	// Pairs of nxt, e
	// Negative e means new block; otherwise, it's 2*e + (has_nontrivial_lowval2)
	std::vector<std::pair<int, int>> adj_lists;
	std::vector<array_range_t<std::pair<int, int>, &spqr_tree::adj_lists>> adj;

	struct bucket_edge_t {
		int v; int e; int cur, nxt;
	};
	std::vector<bucket_edge_t> bucket_edges;

	std::pair<int, int> dfs_lowval(int cur, int d, int prvE) {
		depth[cur] = d;
		int v1 = d, v2 = d;
		for (auto [nxt, e] : adj[cur].bind(*this)) {
			if (e == prvE) continue;
			if (depth[nxt] == -1) {
				const auto [n1, n2] = dfs_lowval(nxt, d+1, e);

				if (n1 >= d) {
					bucket_edges.push_back({0, ~e, cur, nxt});
				} else {
					bucket_edges.push_back({1 + n1 * 2 + (n2 < d), 2*e + (n2 < d), cur, nxt});
				}

				if (n1 < v1) {
					v2 = std::min(v1, n2);
					v1 = n1;
				} else if (n1 == v1) {
					v2 = std::min(v2, n2);
				} else {
					v2 = std::min(v2, n1);
				}
			} else if (depth[nxt] <= d) {
				int nd = depth[nxt];
				if (nd == d) {
					// Self-loop
					bucket_edges.push_back({0, ~e, cur, nxt});
				} else {
					bucket_edges.push_back({1 + nd * 2, 2*e, cur, nxt});
				}

				if (nd < v1) {
					v2 = v1;
					v1 = nd;
				} else if (v1 < nd) {
					v2 = std::min(v2, nd);
				}
			} else {
				// reverse of a back-edge (downwards)
				// skip it
			}
		}
		return {v1, v2};
	}

	void build_sorted_adj(std::vector<std::array<int, 2>> edges, int root) {
		{
			std::vector<int> deg(NV, 0);
			for (auto e : edges) {
				for (int v : e) {
					deg[v]++;
				}
			}
			adj.resize(NV);
			int cur = 0;
			for (int i = 0; i < NV; i++) {
				adj[i].st = adj[i].en = cur;
				cur += deg[i];
			}
			adj_lists.resize(cur);
		}
		for (int e = 0; e < int(edges.size()); e++) {
			auto [u, v] = edges[e];
			adj_lists[adj[u].en++] = {v, e};
			if (u != v) {
				adj_lists[adj[v].en++] = {u, e};
			}
		}

		depth = std::vector<int>(NV, -1);

		// Bucketed so that bridges come first, then BCCs, then children from shallowest to deepest lowval
		bucket_edges.reserve(NE);

		for (int rt = 0; rt < NV; rt++) {
			if (root != -1 && rt != root) continue;
			if (depth[rt] != -1) continue;
			dfs_lowval(rt, 0, -1);
		}
		assert(int(bucket_edges.size()) == NE);

		for (auto& v : adj) {
			v.en = v.st;
		}

		std::vector<int> bucket_st(1 + NV * 2);
		for (int i = 0; i < int(bucket_edges.size()); i++) {
			++bucket_st[bucket_edges[i].v];
		}
		{
			int cur = 0;
			for (int i = 0; i < int(bucket_st.size()); i++) {
				bucket_st[i] = std::exchange(cur, cur + bucket_st[i]);
			}
			assert(cur == int(bucket_edges.size()));
		}
		std::vector<bucket_edge_t> bucket_edges_2(bucket_edges.size());
		for (int i = 0; i < int(bucket_edges.size()); i++) {
			bucket_edges_2[bucket_st[bucket_edges[i].v]++] = bucket_edges[i];
		}
		bucket_edges = {};
		for (auto [_, e, cur, nxt] : bucket_edges_2) {
			adj_lists[adj[cur].en++] = {nxt, e};
		}
	}

	// Tree edges have v = {child, parent}
	// Backedges edges have v = {top, bottom}

	struct vestack_t {
		std::array<int, 2> vs;
		bool is_tree = false;
		int o_ve = -1;
		int o_node = -1;
		node_type o_type;
	};
	std::vector<vestack_t> vestack;

	using vestack_range_t = enumerable_array_range_t<vestack_t, &spqr_tree::vestack>;

	struct estack_t {
		std::array<int, 2> vs;
		vestack_range_t vestack_range;
		node_type type;
		bool is_tree;
	};
	std::vector<estack_t> estack;

	struct tstack_t {
		int idx; // idx in estack
		int vstart;
		int top_depth;
	};
	std::vector<tstack_t> tstack;

	std::vector<int> first_occurrence;

	int cur_component;
	int cur_block;

	// ==== Per-skeleton planarity testing and embedding ====
	//
	// Each sealed node's skeleton is tested and embedded as it is finalized,
	// from the vedge range alone. The storage order (the order the
	// corresponding estack entries were pushed) is a DFS of the skeleton:
	// every tree vedge appears after its whole subtree, so one scan computes
	// skeleton lowpoints bottom-up. Every vedge is a collapsed 2-attachment
	// component whose interlacement behavior depends only on its endpoints,
	// so the left-right planarity test (Brandes' presentation) runs on the
	// skeleton directly, with each vedge acting as a single (tree or back)
	// edge; the two-coloring is kept implicit in lr_ref/lr_side and resolved
	// at the end to emit the (unique up to reflection) rotation.
	struct lr_interval {
		int hi = -1, lo = -1;
		bool empty() const { return lo == -1; }
	};
	struct lr_pair_t {
		lr_interval L, R;
	};
	struct lr_frame_t {
		int sv;         // skeleton vertex index
		int pe;         // the inbound tree edge id (-1 for the root)
		int sb;         // lr_conflicts size when entered
		int lowpt_pe;   // lowpt of the first returning outgoing edge
		int lowpt_edge; // its lowpt edge (a chord id)
	};
	// scratch buffers, reused across nodes
	std::vector<lr_pair_t> lr_conflicts;
	std::vector<lr_frame_t> lr_frames;
	std::vector<int> lr_lowpt, lr_low2, lr_lowedge, lr_ref, lr_chain;
	std::vector<signed char> lr_side;
	std::vector<std::array<int, 2>> lr_out; // (skeleton vertex index, id)
	std::vector<int> v_stamp, v_map;        // vertex-indexed scratch
	int v_stamp_cnt = 0;
	std::vector<int> sv_verts, sv_acc1, sv_acc2;
	std::vector<int> sv_out_st, sv_out, sv_tmp_l, sv_tmp_r;
	std::vector<int> sv_head, sv_lref, sv_rref, sv_suf;
	std::vector<int> rot_next, rot_prev;
	std::vector<std::array<int, 2>> gen_verts;

	// Fallback rotation (used for non-planar R skeletons): group the
	// half-edges of each skeleton vertex in storage order.
	void emit_generic(const node_t& n) {
		++v_stamp_cnt;
		gen_verts.clear();
		for (int ve = n.vedges.st; ve < n.vedges.en; ve++) {
			for (int k = 0; k < 2; k++) {
				int v = vedges[ve].vs[k];
				int h = 2 * ve + k;
				if (v_stamp[v] != v_stamp_cnt) {
					v_stamp[v] = v_stamp_cnt;
					v_map[v] = int(gen_verts.size());
					gen_verts.push_back({h, h});
				} else {
					auto& g = gen_verts[v_map[v]];
					embed_next[g[1]] = h;
					g[1] = h;
				}
			}
		}
		for (auto [first, last] : gen_verts) {
			embed_next[last] = first;
		}
	}

	// S/O skeletons are a single cycle, already in cycle order.
	void emit_cycle(const node_t& n) {
		int st = n.vedges.st, k = n.vedges.en - n.vedges.st;
		for (int z = 0; z < k; z++) {
			int zn = (z + 1 == k) ? 0 : z + 1;
			// consecutive vedges share a vertex: vs[1] of z == vs[0] of zn
			assert(vedges[st + z].vs[1] == vedges[st + zn].vs[0]);
			int a = 2 * (st + z) + 1, b = 2 * (st + zn) + 0;
			embed_next[a] = b;
			embed_next[b] = a;
		}
	}

	// P skeletons: storage order at one endpoint, reversed at the other.
	void emit_parallel(const node_t& n) {
		int st = n.vedges.st, k = n.vedges.en - n.vedges.st;
		int a = vedges[st].vs[0], b = vedges[st].vs[1];
		auto half_at = [&](int z, int v) {
			return 2 * (st + z) + (vedges[st + z].vs[0] == v ? 0 : 1);
		};
		for (int z = 0; z < k; z++) {
			int zn = (z + 1 == k) ? 0 : z + 1;
			embed_next[half_at(z, a)] = half_at(zn, a);
			embed_next[half_at(zn, b)] = half_at(z, b);
		}
	}

	void emit_rigid(node_t& n) {
		const int st = n.vedges.st;
		const int m = n.vedges.en - n.vedges.st; // ids 0..m-1; the cap is id m-1
		const int cap_id = m - 1;
		constexpr int INF = std::numeric_limits<int>::max();

		// Each vedge id is a skeleton tree edge (is_tree, oriented deep->shallow
		// child->parent) or a chord (oriented shallow->deep target->source).
		auto id_is_tree = [&](int id) { return vedges[st + id].is_tree; };
		auto id_child = [&](int id) { // tree: the child (deeper) endpoint
			const auto& ve = vedges[st + id];
			assert(depth[ve.vs[0]] > depth[ve.vs[1]]);
			return ve.vs[0];
		};
		auto id_parent = [&](int id) {
			const auto& ve = vedges[st + id];
			assert(depth[ve.vs[0]] > depth[ve.vs[1]]);
			return ve.vs[1];
		};
		auto id_src = [&](int id) { // chord: the source (deeper) endpoint
			const auto& ve = vedges[st + id];
			assert(depth[ve.vs[1]] > depth[ve.vs[0]]);
			return ve.vs[1];
		};
		auto id_tgt = [&](int id) {
			const auto& ve = vedges[st + id];
			assert(depth[ve.vs[1]] > depth[ve.vs[0]]);
			return ve.vs[0];
		};
		// the event vertex: the skeleton-DFS vertex the edge is outgoing from
		auto id_from = [&](int id) { return id_is_tree(id) ? id_parent(id) : id_src(id); };

		// ==== 1. Index the skeleton vertices ====
		++v_stamp_cnt;
		sv_verts.clear();
		for (int id = 0; id < m; id++) {
			for (int v : vedges[st + id].vs) {
				if (v_stamp[v] != v_stamp_cnt) {
					v_stamp[v] = v_stamp_cnt;
					v_map[v] = int(sv_verts.size());
					sv_verts.push_back(v);
				}
			}
		}
		int nsv = int(sv_verts.size());
		// the skeleton root: the shallower split vertex
		const int rt = std::min(vedges[st + cap_id].vs[0], vedges[st + cap_id].vs[1],
				[&](int a, int b) { return depth[a] < depth[b]; });

		// ==== 2. Skeleton lowpoints and nesting order ====
		//
		// The storage order is a DFS of the skeleton (each tree vedge appears
		// after its whole subtree; a chord cap is the root's earliest-nested
		// outgoing edge and is scanned first), so one scan computes each
		// edge's two lowest return depths bottom-up. The nesting order is NOT
		// simply the storage order: collapsing a subtree into a virtual tree
		// edge can change its lowpoints relative to the original graph's, so
		// each vertex's outgoing edges are re-sorted by the skeleton nesting
		// key (as in build_sorted_adj: 2*lowpt + has_nontrivial_lowpt2).
		sv_acc1.assign(nsv, INF);
		sv_acc2.assign(nsv, INF);
		lr_lowpt.assign(m, -1);
		lr_low2.assign(m, -1);
		lr_lowedge.assign(m, -1);
		lr_out.clear();
		auto scan_edge = [&](int id) {
			int n1, n2, from;
			if (!id_is_tree(id)) {
				from = id_src(id);
				n1 = depth[id_tgt(id)];
				n2 = depth[from];
			} else {
				int c = id_child(id);
				from = id_parent(id);
				int i = v_map[c];
				n1 = std::min(depth[c], sv_acc1[i]);
				n2 = std::min(std::max(depth[c], sv_acc1[i]), sv_acc2[i]);
			}
			lr_lowpt[id] = n1;
			lr_low2[id] = n2;
			int i = v_map[from];
			if (n1 < sv_acc1[i]) {
				sv_acc2[i] = std::min(sv_acc1[i], n2);
				sv_acc1[i] = n1;
			} else {
				sv_acc2[i] = std::min(sv_acc2[i], std::max(n1, std::min(sv_acc1[i], n2)));
			}
			lr_out.push_back({i, id});
		};
		const bool cap_is_tree = id_is_tree(cap_id);
		if (!cap_is_tree) scan_edge(cap_id);
		for (int id = 0; id < cap_id; id++) scan_edge(id);
		if (cap_is_tree) scan_edge(cap_id);

		auto nesting = [&](int id) {
			return 2 * lr_lowpt[id] + (lr_low2[id] < depth[id_from(id)]);
		};
		std::stable_sort(lr_out.begin(), lr_out.end(), [&](std::array<int, 2> a, std::array<int, 2> b) {
			return a[0] != b[0] ? a[0] < b[0] : nesting(a[1]) < nesting(b[1]);
		});
		sv_out_st.assign(nsv + 1, 0);
		sv_out.resize(m);
		for (int z = 0; z < m; z++) {
			sv_out[z] = lr_out[z][1];
			sv_out_st[lr_out[z][0] + 1]++;
		}
		for (int i = 0; i < nsv; i++) {
			sv_out_st[i + 1] += sv_out_st[i];
		}

		// ==== 3. Left-right planarity test (Brandes' presentation) ====
		//
		// A DFS of the skeleton in nesting order: chords push singleton
		// conflict pairs, returning tree edges merge their subtree's
		// constraints into their parent's, backtracking trims return edges.
		// The two-coloring is kept implicit in lr_ref/lr_side and resolved
		// afterwards.
		lr_conflicts.clear();
		lr_frames.clear();
		lr_ref.assign(m, -1);
		lr_side.assign(m, 1);

		bool ok = true;

		auto lowest = [&](const lr_pair_t& P) {
			if (P.L.empty()) return lr_lowpt[P.R.lo];
			if (P.R.empty()) return lr_lowpt[P.L.lo];
			return std::min(lr_lowpt[P.L.lo], lr_lowpt[P.R.lo]);
		};

		// Trim back edges returning to the vertex at depth ud (Brandes,
		// Algorithm 5).
		auto trim = [&](int ud) {
			while (!lr_conflicts.empty() && lowest(lr_conflicts.back()) == ud) {
				lr_pair_t P = lr_conflicts.back();
				lr_conflicts.pop_back();
				if (P.L.lo != -1) lr_side[P.L.lo] = -1;
			}
			if (!lr_conflicts.empty()) {
				lr_pair_t P = lr_conflicts.back();
				lr_conflicts.pop_back();
				while (P.L.hi != -1 && lr_lowpt[P.L.hi] == ud) P.L.hi = lr_ref[P.L.hi];
				if (P.L.hi == -1 && P.L.lo != -1) {
					lr_ref[P.L.lo] = P.R.lo;
					lr_side[P.L.lo] = -1;
					P.L.lo = -1;
				}
				while (P.R.hi != -1 && lr_lowpt[P.R.hi] == ud) P.R.hi = lr_ref[P.R.hi];
				if (P.R.hi == -1 && P.R.lo != -1) {
					lr_ref[P.R.lo] = P.L.lo;
					lr_side[P.R.lo] = -1;
					P.R.lo = -1;
				}
				if (!(P.L.empty() && P.R.empty())) lr_conflicts.push_back(P);
			}
		};

		// Merge the conflict pairs generated by edge ei into the constraints
		// of its parent edge (Brandes, Algorithm 4).
		auto add_constraints = [&](int ei, int sb, int lowpt_pe, int lowpt_edge_pe) -> bool {
			lr_pair_t P;
			// merge return edges of ei into P.R
			do {
				lr_pair_t Q = lr_conflicts.back();
				lr_conflicts.pop_back();
				if (!Q.L.empty()) std::swap(Q.L, Q.R);
				if (!Q.L.empty()) return false;
				if (lr_lowpt[Q.R.lo] > lowpt_pe) {
					// merge intervals
					if (P.R.empty()) {
						P.R.hi = Q.R.hi;
					} else {
						lr_ref[P.R.lo] = Q.R.hi;
					}
					P.R.lo = Q.R.lo;
				} else {
					// align with the parent's lowpoint edge
					lr_ref[Q.R.lo] = lowpt_edge_pe;
				}
			} while (int(lr_conflicts.size()) != sb);
			// merge conflicting return edges of earlier siblings into P.L
			auto conflicting = [&](const lr_interval& I) {
				return !I.empty() && lr_lowpt[I.hi] > lr_lowpt[ei];
			};
			while (!lr_conflicts.empty() && (conflicting(lr_conflicts.back().L) || conflicting(lr_conflicts.back().R))) {
				lr_pair_t Q = lr_conflicts.back();
				lr_conflicts.pop_back();
				if (conflicting(Q.R)) std::swap(Q.L, Q.R);
				if (conflicting(Q.R)) return false;
				// merge interval below lowpt(ei) into P.R
				if (!Q.R.empty()) {
					if (P.R.empty()) {
						P.R.hi = Q.R.hi;
					} else {
						lr_ref[P.R.lo] = Q.R.hi;
					}
					P.R.lo = Q.R.lo;
				}
				if (P.L.empty()) {
					P.L.hi = Q.L.hi;
				} else {
					lr_ref[P.L.lo] = Q.L.hi;
				}
				P.L.lo = Q.L.lo;
			}
			if (!(P.L.empty() && P.R.empty())) lr_conflicts.push_back(P);
			return true;
		};

		// Integrate the just-traversed outgoing edge id (with conflict pairs
		// above sb) into the constraints of frame f's inbound tree edge.
		auto integrate = [&](lr_frame_t& f, int id, int sb) {
			if (lr_lowpt[id] < depth[sv_verts[f.sv]]) {
				// id has a return edge
				if (f.lowpt_edge == -1) {
					f.lowpt_edge = lr_lowedge[id];
					f.lowpt_pe = lr_lowpt[id];
				} else if (!add_constraints(id, sb, f.lowpt_pe, f.lowpt_edge)) {
					ok = false;
				}
			}
		};

		sv_suf.assign(sv_out_st.begin(), sv_out_st.end() - 1); // as cursors
		lr_frames.push_back({v_map[rt], -1, 0, -1, -1});
		while (ok && !lr_frames.empty()) {
			lr_frame_t& f = lr_frames.back();
			int i = f.sv;
			if (sv_suf[i] < sv_out_st[i + 1]) {
				int id = sv_out[sv_suf[i]++];
				if (id_is_tree(id)) {
					lr_frames.push_back({v_map[id_child(id)], id, int(lr_conflicts.size()), -1, -1});
				} else {
					lr_lowedge[id] = id;
					int sb = int(lr_conflicts.size());
					lr_conflicts.push_back({{}, {id, id}});
					integrate(f, id, sb);
				}
			} else {
				// all outgoing edges done: backtrack over the inbound edge
				int pe = f.pe;
				lr_frame_t fin = f;
				lr_frames.pop_back();
				if (pe == -1) break; // the root is done
				lr_lowedge[pe] = fin.lowpt_edge;
				int ud = depth[id_parent(pe)];
				trim(ud);
				if (lr_lowpt[pe] < ud) {
					// the side of pe is determined by the highest return edge
					assert(!lr_conflicts.empty());
					int hL = lr_conflicts.back().L.hi;
					int hR = lr_conflicts.back().R.hi;
					if (hL != -1 && (hR == -1 || lr_lowpt[hL] > lr_lowpt[hR])) {
						lr_ref[pe] = hL;
					} else {
						lr_ref[pe] = hR;
					}
				}
				integrate(lr_frames.back(), pe, fin.sb);
			}
		}
		assert(!ok || lr_conflicts.empty());

		n.planar = ok;
		if (!ok) {
			emit_generic(n);
			return;
		}

		// ==== 4. Emit the rotation (Brandes' embedding phase) ====

		// Resolve the implicit two-coloring: lr_side relative to lr_ref
		// chains becomes an absolute sign.
		lr_chain.clear();
		for (int id0 = 0; id0 < m; id0++) {
			int x = id0;
			while (lr_ref[x] != -1) {
				lr_chain.push_back(x);
				x = lr_ref[x];
			}
			signed char s = lr_side[x];
			while (!lr_chain.empty()) {
				int y = lr_chain.back();
				lr_chain.pop_back();
				s = (signed char)(lr_side[y] * s);
				lr_side[y] = s;
				lr_ref[y] = -1;
			}
		}

		// Reorder each vertex's outgoing edges by signed nesting depth: the
		// left edges reversed, then the right edges.
		for (int i = 0; i < nsv; i++) {
			sv_tmp_l.clear();
			sv_tmp_r.clear();
			for (int z = sv_out_st[i]; z < sv_out_st[i + 1]; z++) {
				(lr_side[sv_out[z]] == 1 ? sv_tmp_r : sv_tmp_l).push_back(sv_out[z]);
			}
			int z = sv_out_st[i];
			for (int j = int(sv_tmp_l.size()) - 1; j >= 0; j--) sv_out[z++] = sv_tmp_l[j];
			for (int id : sv_tmp_r) sv_out[z++] = id;
		}

		// Initial rotations: the signed-order outgoing half-edges.
		rot_next.assign(2 * m, -1);
		rot_prev.assign(2 * m, -1);
		sv_head.assign(nsv, -1);
		sv_lref.assign(nsv, -1);
		sv_rref.assign(nsv, -1);

		auto half_at = [&](int id, int v) {
			assert(vedges[st + id].vs[0] == v || vedges[st + id].vs[1] == v);
			return 2 * id + (vedges[st + id].vs[0] == v ? 0 : 1);
		};
		auto make_first = [&](int i, int h) {
			int a = sv_head[i];
			rot_prev[h] = -1;
			rot_next[h] = a;
			if (a != -1) rot_prev[a] = h;
			sv_head[i] = h;
		};
		auto insert_after = [&](int a, int h) {
			assert(a != -1);
			int b = rot_next[a];
			rot_next[a] = h;
			rot_prev[h] = a;
			rot_next[h] = b;
			if (b != -1) rot_prev[b] = h;
		};
		auto insert_before = [&](int i, int a, int h) {
			assert(a != -1);
			int b = rot_prev[a];
			rot_prev[a] = h;
			rot_next[h] = a;
			rot_prev[h] = b;
			if (b != -1) rot_next[b] = h;
			else sv_head[i] = h;
		};
		auto place_chord = [&](int id) {
			int w = id_tgt(id);
			int j = v_map[w];
			int h = half_at(id, w);
			if (lr_side[id] == 1) {
				insert_after(sv_rref[j], h);
			} else {
				insert_before(j, sv_lref[j], h);
				sv_lref[j] = h;
			}
		};

		for (int i = 0; i < nsv; i++) {
			int v = sv_verts[i];
			int tail = -1;
			for (int z = sv_out_st[i]; z < sv_out_st[i + 1]; z++) {
				int h = half_at(sv_out[z], v);
				if (tail == -1) sv_head[i] = h;
				else {
					rot_next[tail] = h;
					rot_prev[h] = tail;
				}
				tail = h;
			}
		}
		// DFS over the skeleton tree in the signed order: descending a tree
		// edge anchors it at the parent, puts its half-edge first in the
		// child's rotation, and recurses; chords are placed at their target's
		// anchors. sv_suf holds each vertex's cursor, sv_tmp_l is the stack.
		sv_suf.assign(sv_out_st.begin(), sv_out_st.end() - 1);
		sv_tmp_l.assign(1, v_map[rt]);
		while (!sv_tmp_l.empty()) {
			int i = sv_tmp_l.back();
			bool descended = false;
			while (sv_suf[i] < sv_out_st[i + 1]) {
				int id = sv_out[sv_suf[i]++];
				if (id_is_tree(id)) {
					int cv = id_child(id);
					int ci = v_map[cv];
					sv_lref[i] = sv_rref[i] = half_at(id, sv_verts[i]);
					make_first(ci, half_at(id, cv));
					sv_tmp_l.push_back(ci);
					descended = true;
					break;
				}
				place_chord(id);
			}
			if (!descended) sv_tmp_l.pop_back();
		}

		// Close the rotations into circular lists, in global half-edge ids.
		for (int i = 0; i < nsv; i++) {
			int h0 = sv_head[i];
			assert(h0 != -1);
			for (int h = h0;;) {
				int nh = rot_next[h];
				embed_next[2 * st + h] = 2 * st + (nh == -1 ? h0 : nh);
				if (nh == -1) break;
				h = nh;
			}
		}
	}

	void emit_node_embedding(node_t& n) {
		switch (n.type) {
		case node_type::I:
			// a single vedge; the identity rotations are already correct
			break;
		case node_type::O:
		case node_type::S:
			emit_cycle(n);
			break;
		case node_type::P:
			emit_parallel(n);
			break;
		case node_type::R:
			emit_rigid(n);
			break;
		default:
			assert(false);
		}
	}

	void finalize_node(estack_t es, vestack_t cap) {
		assert(es.type != node_type::Q);
		int node = int(nodes.size());

		node_t n;
		n.type = es.type;
		n.component = cur_component;
		n.block = cur_block;

		n.vedges.st = int(vedges.size());
		n.node_vertices.st = int(node_vertices.size());

		auto push_vestack_t = [&](const vestack_t& ves) {
			int ve = int(vedges.size());
			embed_next.push_back(2 * ve);
			embed_next.push_back(2 * ve + 1);
			vedges.push_back(vedge_t{
				ves.vs,
				cur_component,
				cur_block,
				node,
				ves.is_tree,
				ves.o_ve,
				ves.o_node,
				ves.o_type,
			});
			int o_ve = ves.o_ve;
			if (o_ve != -1) {
				vedges[o_ve].o_ve = ve;
				vedges[o_ve].o_node = node;
				vedges[o_ve].o_type = n.type;
			}
		};

		for (const auto& ves : es.vestack_range.bind(*this)) {
			push_vestack_t(ves);
			if (ves.is_tree) {
				int cnd = ves.vs[0];
				if (cnd != cap.vs[0] && cnd != cap.vs[1]) {
					node_vertices.push_back(cnd);
				}
			}
		}
		assert(cap.vs[0] == es.vs[1]);
		assert(cap.vs[1] == es.vs[0]);
		assert(cap.is_tree == !es.is_tree);
		push_vestack_t(cap);
		if (cap.vs[0] == cap.vs[1]) {
			node_vertices.push_back(cap.vs[0]);
		} else {
			for (int v : cap.vs) {
				node_vertices.push_back(v);
			}
		}

		NVE = int(vedges.size());

		n.vedges.en = int(vedges.size());
		n.node_vertices.en = int(node_vertices.size());

		emit_node_embedding(n);
		blocks[cur_block].planar = blocks[cur_block].planar && n.planar;
		is_planar = is_planar && n.planar;

		NN++;
		nodes.push_back(n);
	}

	vestack_t finalize_estack(estack_t es) {
		vestack_t cap;
		cap.vs[0] = es.vs[1];
		cap.vs[1] = es.vs[0];
		cap.is_tree = !es.is_tree;

		int o_ve;
		int o_node;
		if (es.type == node_type::Q) {
			int e = vestack[es.vestack_range.st].o_ve;
			o_ve = o_node = e;
			assert(cap.vs == vedges[e].vs && cap.is_tree == vedges[e].is_tree);
		} else {
			finalize_node(es, cap);
			o_ve = int(vedges.size()) - 1;
			o_node = int(nodes.size()) - 1;
		}

		vestack_t ve;
		ve.vs = es.vs;
		ve.is_tree = es.is_tree;
		ve.o_ve = o_ve;
		ve.o_node = o_node;
		ve.o_type = es.type;
		return ve;
	}

	void push_estack(estack_t e_ins) {
		if constexpr (PRINT_DEBUG) std::cerr << "push_estack " << e_ins.vs[0] << '-' << e_ins.vs[1] << ' ' << e_ins.is_tree << '\n';
		estack.push_back(e_ins);
		int v = e_ins.vs[0];
		if (first_occurrence[v] == -1) {
			first_occurrence[v] = int(estack.size()) - 1;
		}
	}

	void push_estack_p(estack_t e_ins) {
		if constexpr (PRINT_DEBUG) std::cerr << "push_estack_p " << e_ins.vs[0] << '-' << e_ins.vs[1] << ' ' << e_ins.is_tree << '\n';
		if (estack.back().type == node_type::P) {
			int st = e_ins.vestack_range.st;
			vestack[st] = finalize_estack(e_ins);
			vestack.resize(st+1);
			assert(estack.back().vestack_range.en == st);
			estack.back().vestack_range.en ++;
		} else {
			int st = estack.back().vestack_range.st;
			vestack[st] = finalize_estack(estack.back());
			vestack[st+1] = finalize_estack(e_ins);
			vestack.resize(st+2);
			estack.back().vestack_range = {st, st+2};
			estack.back().type = node_type::P;
		}
	}

	void prepare_pop_estack(int z) {
		int v = estack[z].vs[0];
		if (first_occurrence[v] == z) {
			first_occurrence[v] = -1;
		}
	}
	void pop_estack() {
		prepare_pop_estack(int(estack.size()) - 1);
		estack.pop_back();
	}

	// push_estack_p is a lot like push_estack + pop_estack_range, but their
	// handling of first_occurrence would differ, so we'll leave them separate
	std::pair<vestack_range_t, node_type> pop_estack_range(int idx) {
		assert(int(estack.size()) - idx > 1);
		for (int z = idx; z < int(estack.size()); z++) {
			prepare_pop_estack(z);
		}

		bool is_S = (int(estack.size()) - idx == 2);
		if (is_S) {
			assert(estack.back().type != node_type::S);
		}
		bool should_reuse = (is_S && estack[idx].type == node_type::S);
		int st = estack[idx].vestack_range.st;

		int sidx = idx + should_reuse;
		int en = estack[sidx].vestack_range.st;
		for (int z = sidx; z < int(estack.size()); z++) {
			vestack[en++] = finalize_estack(estack[z]);
		}

		vestack.resize(en);
		estack.resize(idx);
		return {{st, en}, is_S ? node_type::S : node_type::R};
	}

	estack_t make_q_node(int e, int cur, int nxt, bool is_tree) {
		vedges[e].vs = {cur, nxt};
		vedges[e].component = cur_component;
		vedges[e].block = cur_block;
		vedges[e].node = e;
		vedges[e].is_tree = !is_tree;
		if (cur == nxt) {
			embed_next[2*e] = 2*e + 1;
			embed_next[2*e + 1] = 2*e;
		} else {
			embed_next[2*e] = 2*e;
			embed_next[2*e + 1] = 2*e + 1;
		}
		// We're gonna put a dummy entry on vestack for 2 reasons:
		// 1. This lets the pop_estack_range logic has 1 scratch space per estack
		// 2. We need to lookup the edge for finalize_estack, which we'll smuggle in as o_ve
		vestack.emplace_back();
		vestack.back().o_ve = e;

		nodes[e].type = node_type::Q;
		nodes[e].component = cur_component;
		nodes[e].block = cur_block;
		nodes[e].vedges = {e, e+1};
		nodes[e].node_vertices = {2*e, 2*e + (cur == nxt ? 1 : 2)};
		node_vertices[2*e] = cur;
		node_vertices[2*e+1] = nxt;

		return estack_t{{nxt, cur}, {int(vestack.size()) - 1, int(vestack.size())}, node_type::Q, is_tree};
	}

	void dfs_spqr(int cur, int cur_low) {
		int cur_depth = depth[cur];
		for (auto [nxt, e] : adj[cur].bind(*this)) {
			if (e < 0) continue;
			bool is_type_1 = !(e & 1);
			e >>= 1;

			int orig_size = int(estack.size());

			bool is_tree = (depth[nxt] > cur_depth);
			if (is_tree) {
				dfs_spqr(nxt, cur_low);

				// Before we insert our edge, we may have to fix the single-tree-edge tstack
				assert(!tstack.empty());
				if (tstack.back().top_depth == cur_depth + 1) {
					assert(tstack.back().idx == int(estack.size())-1);
					tstack.back().top_depth = cur_depth;
				}
			}

			// Deal with the current vedge
			estack_t e_ins = make_q_node(e, cur, nxt, is_tree);
			if (is_tree) {
				while (true) {
					if (estack.back().vs == std::array<int, 2>{cur, nxt}) {
						push_estack_p(e_ins);
						e_ins.vestack_range = estack.back().vestack_range;
						assert(estack.back().type == node_type::P);
						e_ins.type = node_type::P;

						assert(!tstack.empty());
						if (tstack.back().idx == int(estack.size()) - 1) {
							tstack.pop_back();
						}
						// Pop it all the way off since we need to
						// change the backedge to a tree edge
						pop_estack();
					}

					push_estack(e_ins);

					assert(!tstack.empty());
					if (tstack.back().top_depth != cur_depth) break;
					assert(tstack.back().idx > orig_size);
					if (estack[tstack.back().idx].is_tree) {
						assert(estack[tstack.back().idx].vs[0] == tstack.back().vstart);
						if (int(estack.size()) - tstack.back().idx > 2) {
							tstack.push_back({
								tstack.back().idx + 1,
								estack[tstack.back().idx].vs[1],
								tstack.back().top_depth
							});
						}
					}

					int idx = tstack.back().idx;
					assert(idx > orig_size);
					nxt = tstack.back().vstart;
					tstack.pop_back();
					auto [vestack_range, type] = pop_estack_range(idx);
					e_ins = estack_t{{nxt, cur}, vestack_range, type, true};
				}

				if (is_type_1) {
					// Handle type 1 split

					nxt = estack[orig_size].vs[0];

					// It's possible there are candidates before orig_size
					// on the stack, but those are doomed to fail, so
					// it's fine to pop them prematurely.
					while (tstack.back().idx > orig_size) tstack.pop_back();
					assert(!tstack.empty());
					if (tstack.back().idx == orig_size) {
						assert(tstack.back().vstart == cur_low && tstack.back().top_depth == depth[nxt]);
					}

					int idx = orig_size;
					auto [vestack_range, type] = pop_estack_range(idx);
					e_ins = estack_t{{nxt, cur}, vestack_range, type, false};
					if (!estack.empty() && estack.back().vs == e_ins.vs) {
						push_estack_p(e_ins);
						// Pop our tstack and use the previous value
						if (tstack.back().idx == orig_size) tstack.pop_back();
					} else {
						push_estack(e_ins);
						// We don't need to change tstack at all
					}
				} else if (cur_low == cur) {
					// Pop until we get to the entire range; that's
					// the only permissible one, since we previously use something
					while (tstack.back().idx > orig_size) tstack.pop_back();
					assert(!tstack.empty());
				} else if (first_occurrence[cur] != -1) {
					// We're the first node, so cur can be interior, but we must include all backedges to cur
					while (tstack.back().idx > first_occurrence[cur]) tstack.pop_back();
					assert(!tstack.empty());
				} else {
					// This is the only case where we can have the tree-edge start the cut.
					// Notably, this is the only tstack entry with top == cur
					// (otherwise it's always strictly shallower)
					// By the time it's resolved, top will be
					// shallower than cur, but that depends on the
					// next edge, so we'll fix it then.
					// We'll also omit the following (equivalent) entry entirely and
					// handle it when we pop it out.
					tstack.push_back({int(estack.size())-1, nxt, cur_depth});
				}
			} else {
				if (!estack.empty() && estack.back().vs == e_ins.vs) {
					push_estack_p(e_ins);
					// Don't change tstack, keep its old low value
				} else {
					push_estack(e_ins);

					// nxt should be shallower
					assert(depth[cur_low] > depth[nxt]);
					tstack_t t_ins{int(estack.size())-1, cur_low, depth[nxt]};
					assert(tstack.empty() || tstack.back().top_depth <= depth[cur_low]);
					while (!tstack.empty() && tstack.back().top_depth > t_ins.top_depth) {
						t_ins.idx = tstack.back().idx;
						t_ins.vstart = tstack.back().vstart;
						tstack.pop_back();
					}
					tstack.push_back(t_ins);
				}
			}

			cur_low = cur;
		}

		block_vertices.push_back(cur);
		vertex_blocks[vertices[cur].vertex_blocks.en++] = cur_block;
	}

	int start_spqr(int cur, int nxt, int e) {
		int block = NB++;
		blocks.emplace_back();
		blocks[block].component = cur_component;
		blocks[block].nodes.st = int(nodes.size());
		blocks[block].vedges.st = int(vedges.size());
		blocks[block].block_vertices.st = int(block_vertices.size());
		cur_block = block;

		if (nxt != cur) {
			dfs_spqr(nxt, nxt);
			if (estack.empty()) {
				estack.push_back(estack_t{{cur, nxt}, {int(vestack.size()), int(vestack.size())}, node_type::I, false});
			} else {
				assert(int(estack.size()) == 1);
				assert((estack.back().vs == std::array<int, 2>{cur, nxt}));
				assert(!estack.back().is_tree);

				if (estack.back().type == node_type::Q) {
					vestack.back() = finalize_estack(estack.back());
					// Don't have to change anything else, since Q nodes always have size 1
					estack.back().type = node_type::O;
				}

				assert(!tstack.empty() && tstack.back().idx == 0);
				tstack.pop_back();
			}
		} else {
			estack.push_back(estack_t{{cur, cur}, {int(vestack.size()), int(vestack.size())}, node_type::O, true});
		}

		// Make the q node and immediately overwrite its spot with our desired one
		auto q_estack = make_q_node(e, cur, nxt, cur != nxt);
		vestack.back() = finalize_estack(q_estack);

		auto es = estack.back();
		pop_estack();
		finalize_node(es, vestack.back());
		vestack.resize(es.vestack_range.st);

		assert(tstack.empty());
		assert(estack.empty());
		assert(vestack.empty());

		block_vertices.push_back(cur);

		cur_block = -1;
		blocks[block].nodes.en = int(nodes.size());
		blocks[block].vedges.en = int(vedges.size());
		blocks[block].block_vertices.en = int(block_vertices.size());
		return block;
	}

	std::vector<int> vertex_blocks_buf;
	void dfs_block(int cur) {
		int buf_st = int(vertex_blocks_buf.size());
		for (auto [nxt, e] : adj[cur].bind(*this)) {
			assert(depth[nxt] <= depth[cur] + 1);
			if (depth[nxt] < depth[cur]) continue;
			if (nxt != cur) {
				dfs_block(nxt);
			}
			if (e < 0) {
				e = ~e;
				// Root of a block
				int block = start_spqr(cur, nxt, e);

				vertex_blocks_buf.push_back(block);
			}
		}

		vertices[cur].vertex_blocks.st = int(vertex_blocks.size());
		vertex_blocks.insert(vertex_blocks.end(), vertex_blocks_buf.begin() + buf_st, vertex_blocks_buf.end());
		vertex_blocks_buf.resize(buf_st);
		vertices[cur].vertex_blocks.en = int(vertex_blocks.size());

		// Placeholder for our parent
		if (depth[cur] > 0) {
			vertex_blocks.push_back(-1);
		}

		component_vertices.push_back(cur);
		vertices[cur].component = cur_component;
	}

	void build_spqr() {
		vertices.resize(NV);
		vertex_blocks.reserve(NV + NE);

		components.reserve(NV);
		component_vertices.reserve(NV);

		blocks.reserve(NE);
		block_vertices.reserve(NV + NE);

		// # nodes <= 2 * # Q nodes
		nodes.reserve(2 * NE);
		NN = NE;
		nodes.resize(NN);
		// # node_vertices <= # block_vertices + 2 * # nodes
		node_vertices.reserve(NV + 5 * NE);
		node_vertices.resize(2 * NE);

		// # vedges ~ 2 * (# nodes - # blocks)
		vedges.reserve(2 * (2 * NE));
		NVE = NE;
		vedges.resize(NVE);

		embed_next.reserve(2 * 2 * (2 * NE));
		embed_next.resize(2 * NE);
		v_stamp.assign(NV, 0);
		v_map.assign(NV, -1);
		v_stamp_cnt = 0;

		vertex_blocks_buf.reserve(NE);
		vestack.reserve(NE);
		estack.reserve(NE);
		tstack.reserve(NE);
		first_occurrence.assign(NV, -1);

		for (int rt = 0; rt < NV; rt++) {
			if (depth[rt] != 0) continue;
			int component = NC++;
			components.emplace_back();
			component_t& c = components[component];
			c.blocks.st = int(blocks.size());
			c.nodes.st = int(nodes.size());
			c.vedges.st = int(vedges.size());
			c.component_vertices.st = int(component_vertices.size());

			cur_component = component;
			dfs_block(rt);
			cur_component = -1;

			c.blocks.en = int(blocks.size());
			c.nodes.en = int(nodes.size());
			c.vedges.en = int(vedges.size());
			c.component_vertices.en = int(component_vertices.size());
		}

		vertex_blocks_buf = {};
		vestack = {};
		estack = {};
		tstack = {};
		first_occurrence = {};

		lr_conflicts = {};
		lr_frames = {};
		lr_lowpt = lr_low2 = lr_lowedge = lr_ref = lr_chain = {};
		lr_side = {};
		lr_out = {};
		v_stamp = v_map = {};
		sv_verts = sv_acc1 = sv_acc2 = {};
		sv_out_st = sv_out = sv_tmp_l = sv_tmp_r = {};
		sv_head = sv_lref = sv_rref = sv_suf = {};
		rot_next = rot_prev = {};
		gen_verts = {};
	}

public:
	spqr_tree() = default;
	explicit spqr_tree(int NV_, std::vector<std::array<int, 2>> edges, int root = -1) {
		NV = NV_;
		if (NV == 0) return;

		if (root != -1) {
			assert(0 <= root && root < NV);
		}

		NE = int(edges.size());

		build_sorted_adj(std::move(edges), root);

		build_spqr();

		adj = {};
		adj_lists = {};
		// Leave depth since it's sometimes useful
		//depth = {};
	}
};
