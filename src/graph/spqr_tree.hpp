#pragma once

#include <algorithm>
#include <vector>
#include <array>
#include <span>
#include <utility>
#include <cassert>
#include <ranges>
#include <ostream>
#include <expected>
#include <type_traits>
#include <variant>

namespace wala {

struct csr_index {
	std::vector<int> bounds;
	std::ranges::iota_view<int, int> indices(int i) const { return std::views::iota(bounds[i], bounds[i+1]); }
	template <std::ranges::contiguous_range R> auto slice(int i, R&& base) const {
		return std::span(base).subspan(bounds[i], bounds[i+1] - bounds[i]);
	}
	int num_rows() const { return bounds.empty() ? 0 : int(bounds.size()) - 1; }
	int num_entries() const { return bounds.empty() ? 0 : bounds.back(); }
};

template <typename T> struct csr : csr_index {
	std::vector<T> dat;
	std::span<T> operator [](int i) { return slice(i, dat); }
	std::span<const T> operator [](int i) const { return slice(i, dat); }
};

struct csr_index_builder {
	std::vector<int> bounds;
	csr_index_builder() = default;
	explicit csr_index_builder(int N) : bounds(N+1) {}
	void count(int k) { bounds[k+1]++; }
	csr_index finalize() && {
		for (int i = 1; i < int(bounds.size()); i++) {
			bounds[i] += bounds[i-1];
		}
		return {std::move(bounds)};
	}
};

template <typename T> struct csr_builder {
	csr_index idx;
	std::vector<T> dat;
	csr_builder() = default;
	explicit csr_builder(csr_index idx_, std::vector<T>&& dat_buf = {}) : idx(std::move(idx_)), dat(std::move(dat_buf)) {
		dat.resize(idx.num_entries());
		if (!idx.bounds.empty()) {
			idx.bounds.pop_back();
			idx.bounds.insert(idx.bounds.begin(), 0);
		}
	}
	explicit csr_builder(csr_index_builder&& idx_builder, std::vector<T>&& dat_buf = {}) : idx{std::move(idx_builder.bounds)}, dat(std::move(dat_buf)) {
		int l = 0;
		for (int i = 1; i < int(idx.bounds.size()); i++) {
			idx.bounds[i] = std::exchange(l, l + idx.bounds[i]);
		}
		dat.resize(l);
	}
	[[nodiscard]] T& push(int k) { return dat[idx.bounds[k+1]++]; }
	[[nodiscard]] csr<T> finalize() && { return { std::move(idx), std::move(dat) }; }
};

struct planar_spqr_tree;

struct spqr_tree {
	// The SPQR tree of a graph is a canonical/"maximal" decomposition of the graph by 2-vertex cuts.
	// The tree consists of nodes which are graphs of virtual edges (vedges), corresponding to nontrivial 2-vertex cuts.
	// Virtual edges are paired, and we can reassemble the graph by gluing nodes at their matching vedges (and removing the vedge).
	// Real edges are represented as special Q nodes which each contain exactly 1 real edge and exactly 1 vedge.
	//
	// Traditionally, the SPQR tree is defined for each biconnected component,
	// but we will embed the SPQR decompositions inside the block-cut tree to get a (rooted) decomposition of the entire graph.
	//
	// As such, we will have a tree of "items", which consist of SPQR nodes, real vertices, and a special "forest root" item:
	//  - Each vertex will be a child of the topmost node which contains it (or the forest root).
	//  - Each block will be a subtree of nodes rooted at a Q edge, which is the child of one of its vertices.
	//
	// Item types:
	//  F - forest - a root node corresponding to the whole forest.
	//  V - vertex - not really a node, just there because they're mixed into the tree a la block/cut tree
	//  Q - real edge - has exactly 1 vedge and 1 real edge
	//  I - bridge - has exactly 1 vedge connecting to a bridge Q node
	//  O - self-loop - has exactly 1 vedge connecting to a self-loop Q node
	//  S - series - a cycle of >= 3 vedges; note that any 2 vertices of the cycle form a cut
	//  P - parallel - a parallel group of >= 3 vedges with the same endpoints
	//  R - rigid - a 3-vertex-connected component
	//
	// Q nodes occur in 2 places: block roots and block leaves.
	// Block leaf Q's simply have no children.
	// Block root Q's have 2 children: their vedge, and their deeper vertex (unless it's a self-loop).
	//
	// Degenerate blocks:
	//  - a block consisting of a self-loop is a Q node connected to an O node.
	//  - a block consisting of a bridge is a Q node connected to an I node.
	//  - a block consisting of exactly 2 parallel edges is represented by 2 glued Q nodes.
	//
	// We have several id spaces:
	//  - items are in preorder
	//  - node_verts (nv's) are each node's vertices, given in node order then s-t order.
	//  - node_edges (ne's) are each node's vedges, given in node order then a s-t order.
	//  - node_adj is each node_vert's incident vedges, given as 2 lists per nv: left/rightwards edges each in reverse s-t order.
	//  - original verts and original edges can be converted to items as vert_item / edge_item
	//
	// Children of a node will be sorted in s-t order.
	// Specifically vertices are sorted, and edges are guaranteed to satisfy the strong "dominance" partial order:
	// if a.nvs[0] <= b.nvs[0] and a.nvs[1] <= b.nvs[1], then a <= b. (In practice, we'll sort by midpoint.)
	// Adjacency lists are sorted as "center-is-longest", which helps make laminar/bracket cases clean.
	//   (5->4) (5->3) (5->2) (5->1) *vertex 5* (5->9) (5->8) (5->7) (5->6)
	// More specifically, node_adj contains two lists per vertex: 2*nv+0 is leftwards and 2*nv+1 is rightwards.
	//
	// All id's are item indices unless clearly nv/ne id's.
	//
	// In general, there are 2 ways to use the SPQR tree: the rooted view and the unrooted view.
	//  - The rooted view uses par / ch walks, and either treats the tree as 1 top-down big decomposition, or walks in paths up/down the tree with LCA-like queries.
	//  - The unrooted view mostly uses nv/ne/nd lists and works locally within a node/sometimes jumps between them.

	enum class node_type : char {
		F = 'F', V = 'V', Q = 'Q', I = 'I', O = 'O', S = 'S', P = 'P', R = 'R'
	};
	friend std::ostream& operator<<(std::ostream& o, node_type t) { return o << char(t); }

	std::vector<int> vert_index;
	std::vector<int> edge_index;
	std::vector<bool> edge_flipped;

	std::vector<int> par;
	std::vector<int> subtree_end;
	std::vector<node_type> types;
	std::vector<int> orig_id;

	csr<int> ch;
	struct node_vert_t {
		int node;
		int vert;
	};
	std::vector<node_vert_t> node_verts;
	csr_index node_nvs;
	// The nv index of a vertex within its parent node
	std::vector<int> vert_par_nv;
	// TODO: Should we store a vert_nodes CSR?

	struct node_edge_t {
		int node;
		int twin_ne;
		// TODO: Should we store the twin node, the twin node type, and/or twin node type == Q?

		std::array<int, 2> nvs;
	};
	std::vector<node_edge_t> node_edges;
	csr_index node_nes;

	struct node_adj_t {
		int ne;
		int dest_nv;
	};
	csr<node_adj_t> node_adj;

	int size() const { return int(par.size()); }

	// vert_order and edge_order are (prefixes of) permutations of vertex / edge ids;
	// listed ids are visited first in the given order, then the rest in id order.
	// Roots are the first unvisited vertices, and DFS children are explored in edge order.
	// Use planar_spqr_tree::build to also compute the planar embeddings.
	static spqr_tree build(
		int NV,
		const std::vector<std::array<int, 2>>& edges,
		bool ternarize = false,
		std::span<const int> vert_order = {},
		std::span<const int> edge_order = {}
	) {
		return build_impl<false>(NV, edges, ternarize, vert_order, edge_order);
	}

protected:
	template <bool with_planarity>
	static std::conditional_t<with_planarity, planar_spqr_tree, spqr_tree> build_impl(
		int NV,
		const std::vector<std::array<int, 2>>& edges,
		bool ternarize,
		std::span<const int> vert_order,
		std::span<const int> edge_order
	);
};

struct planar_embedding {
	// We'll split our edges up into "quarter-edges", indexed according to
	//   4 * edge + 2 * side + dir, where side is v0 vs v1, and dir is cw vs ccw
	//
	//       1     2
	//    v0 ---e--- v1
	//       0     3
	//
	// We can think of a planar embedding as a collection of 3 involutions on quarter-edges:
	// * qe <-> qe ^ 1 maps quarter edges to their opposite side around the endpoint vertex.
	// * qe <-> qe ^ 3 maps quarter edges to their opposite side along the edge (around the face).
	// * qe <-> rot_adj[qe] maps quarter edges to their facing pair.
	// Walking around a vertex is alternating qe ^ 1 and rot_adj[qe], and walking around a face is qe ^ 3 and rot_adj[qe].
	//
	// Partial embeddings are represented with -1's in the rot_adj array.
	// NB: Helpers do not support -1's. It is up to the user to not access these entries!
	std::vector<int> rot_adj;
};

struct planar_spqr_tree : spqr_tree {
	std::vector<bool> node_planar;
	// Planarity adjacencies: ne_rot_adj is an involution of facing quarter-edges, indexed according to:
	// ne_rot_adj[4 * node_edge + 2 * side + dir]
	// Nonplanar nodes have all entries -1.
	planar_embedding ne_embedding;

	static planar_spqr_tree build(
		int NV,
		const std::vector<std::array<int, 2>>& edges,
		bool ternarize = false,
		std::span<const int> vert_order = {},
		std::span<const int> edge_order = {}
	) {
		return build_impl<true>(NV, edges, ternarize, vert_order, edge_order);
	}
};

// Phase 1: build a DFS skeleton with outedges sorted by lowval
struct lowval_storted_skeleton_t {
	std::vector<int> roots;
	struct key_t { int lowval; bool is_tree; bool is_type_1; };
	struct packed_key_t {
		int v;
		friend auto operator <=> (packed_key_t a, packed_key_t b) = default;
		[[nodiscard]] bool is_new_block() const { return v < 6; }
		[[nodiscard]] bool is_type_2() const { return v % 3 == 2; }
		[[nodiscard]] key_t unpack(int cur_depth) const {
			int lowval = v / 3 - 2; if (lowval < 0) lowval = cur_depth + ~lowval;
			int kind = v % 3;
			bool is_tree = kind != 1;
			bool is_type_1 = kind <= 1;
			return {lowval, is_tree, is_type_1 };
		}
	};
	struct outedge_t { int src, dest; int e_side; packed_key_t key; };
	csr<outedge_t> outedges;

	static lowval_storted_skeleton_t build(
		int NV,
		const std::vector<std::array<int, 2>>& edges,
		std::span<const int> vert_order,
		std::span<const int> edge_order
	) {
		// std::min is by reference, which breaks some optimizations
		auto min = [](auto a, auto b) { return a < b ? a : b; };

		int NE = int(edges.size());
		assert(int(vert_order.size()) <= NV);
		assert(int(edge_order.size()) <= NE);

		// Calls f(i) for i in order, then for the remaining i in [0, n) in increasing order.
		auto for_each_in_order = [] [[gnu::always_inline]] (int n, std::span<const int> order, auto f) -> void {
			for (int i : order) f(i);
			if (int(order.size()) == n) return;
			if (order.empty()) {
				for (int i = 0; i < n; i++) f(i);
			} else if (order.size() == 1) {
				for (int i = 0; i < n; i++) {
					if (i != order[0]) f(i);
				}
			} else {
				std::vector<bool> listed(n);
				for (int i : order) listed[i] = true;
				for (int i = 0; i < n; i++) {
					if (!listed[i]) f(i);
				}
			}
		};

		std::vector<int> roots; roots.reserve(NV);
		csr<outedge_t> outedges;
		{
			std::vector<int> depth(NV, -1);
			// 1a: build a normal adjacency list for the initial lowval dfs
			struct edge_t { int dest; int e; };
			csr_index_builder adj_idx_builder(NV);
			for (auto [u, v] : edges) {
				adj_idx_builder.count(u);
				if (u != v) adj_idx_builder.count(v);
			}
			csr_builder<edge_t> adj_builder(std::move(adj_idx_builder));
			for_each_in_order(NE, edge_order, [&] [[gnu::always_inline]] (int e) -> void {
				auto [u, v] = edges[e];
				adj_builder.push(u) = {v, 2 * e + 0};
				if (u != v) adj_builder.push(v) = {u, 2 * e + 1};
			});
			auto adj = std::move(adj_builder).finalize();

			std::vector<outedge_t> all_outedges; all_outedges.reserve(NE);
			// Return the 2 lowvals from this subtree
			struct dfs_stack_t {
				int cur;
				int prv_e;
				std::array<int, 2> lowvals;
				int ch_idx;
				int ch_end;
			};
			std::vector<dfs_stack_t> stk; stk.reserve(NV);
			auto push_vert = [&] [[gnu::always_inline]] (int cur, int prv_e) -> void {
				int d = int(stk.size());
				depth[cur] = d;
				stk.push_back({cur, prv_e, {d, d}, adj.bounds[cur], adj.bounds[cur+1]});
			};
			auto finish_edge = [&] [[gnu::always_inline]] (bool is_tree, std::array<int, 2> n_lowvals) -> void {
				int d = int(stk.size()) - 1;
				auto& s = stk.back();
				int cur = s.cur;
				assert(s.ch_idx < s.ch_end);
				auto [nxt, e] = adj.dat[s.ch_idx];
				auto& lowvals = s.lowvals;
				s.ch_idx++;

				{
					// Extra bit is 0 for type-1 children, 1 for backedges, 2 for children with lowval2
					// Bridges have lowval -2 (kind 0), and components loops have lowval -1 (components are kind 0, loops are kind 1)
					// We group all backedges together/last to avoid breaking a straight-line graph embedding
					int lowval = n_lowvals[0];
					if (lowval >= d) lowval = ~(lowval - d);
					int kind = 2 * (n_lowvals[1] < d) + !is_tree;
					all_outedges.push_back({cur, nxt, e, packed_key_t{3 * (lowval + 2) + kind}});
				}

				// Keep the 2 distinct mins
				if (n_lowvals[0] < lowvals[0]) lowvals = {n_lowvals[0], min(n_lowvals[1], lowvals[0])};
				else lowvals[1] = min(lowvals[1], n_lowvals[0] == lowvals[0] ? n_lowvals[1] : n_lowvals[0]);
			};
			auto start_edge = [&] [[gnu::always_inline]] () -> void {
				int d = int(stk.size()) - 1;
				auto& s = stk.back();
				assert(s.ch_idx < s.ch_end);
				auto [nxt, e] = adj.dat[s.ch_idx];

				if ((e ^ 1) == s.prv_e || depth[nxt] > d) {
					// skip the edge
					s.ch_idx++; return;
				}

				bool is_tree = depth[nxt] == -1;
				if (is_tree) {
					push_vert(nxt, e);
				} else {
					finish_edge(false, {depth[nxt], d});
				}
			};
			auto pop_vert = [&] [[gnu::always_inline]] () -> std::array<int, 2> {
				auto lowvals = stk.back().lowvals;
				stk.pop_back();
				return lowvals;
			};
			for_each_in_order(NV, vert_order, [&] [[gnu::always_inline]] (int rt) -> void {
				if (depth[rt] == -1) {
					roots.push_back(rt);
					push_vert(rt, -1);
					while (true) {
						if (stk.back().ch_idx == stk.back().ch_end) {
							auto lowvals = pop_vert();
							if (stk.empty()) break;
							finish_edge(true, lowvals);
						} else {
							start_edge();
						}
					}
				}
			});

			csr_index_builder by_key_idx_builder(3*NV+6);
			for (auto edge : all_outedges) by_key_idx_builder.count(edge.key.v);
			csr_builder<outedge_t> by_key_builder(std::move(by_key_idx_builder));
			for (auto edge : all_outedges) by_key_builder.push(edge.key.v) = edge;
			csr<outedge_t> by_key = std::move(by_key_builder).finalize();

			csr_index_builder by_src_idx_builder(NV);
			for (auto edge : by_key.dat) by_src_idx_builder.count(edge.src);
			csr_builder<outedge_t> by_src_builder(std::move(by_src_idx_builder), std::move(all_outedges));
			for (auto edge : by_key.dat) by_src_builder.push(edge.src) = edge;
			outedges = std::move(by_src_builder).finalize();
		}

		return {std::move(roots), std::move(outedges)};
	}
};

// Phase 2: the big ear-decomposition-like walk (spqr_walk), shared by spqr_tree::build, planar_embed and can_planar_embed.
// The walk owns the traversal and the nonplanarity detection; each caller supplies a per-entry payload and the hooks
// that build its own output (SPQR items, quarter edge matchings, ...) off it.
// A return of an ear, as seen by the walk: the exposed end (a quarter edge or edge id, owned by the payload's mode) and its depth.
struct spqr_walk_top_t {
	int end = -1;
	int depth = -1;
};

// Quarter edge bookkeeping shared by the embedding modes of spqr_walk (planar_spqr_tree::build and planar_embed).
// Each (v)edge owns 4 quarter edges; matches pairs them up into the rotation system as the walk closes off pieces of ears.
struct spqr_walk_quarter_edges_t {
	using top_t = spqr_walk_top_t;
	// bot_ends[z] are the outer/innermost exposed quarter edges of the walk down side z of an ear
	// (they're connected to the bottommost/topmost vertices of the tree path)
	struct bot_ends_t {
		std::array<std::array<int, 2>, 2> v{{{-1, -1}, {-1, -1}}};
		std::array<int, 2>& operator [] (int z) { return v[z]; }
		const std::array<int, 2>& operator [] (int z) const { return v[z]; }
	};

	std::vector<int>& matches;
	// top_depth of the entry that pushed each (v)edge, indexed by quarter >> 2
	std::vector<int>& edge_top_depths;

	[[gnu::always_inline]] void link(int a, int b) const {
		matches[a] = b;
		matches[b] = a;
	}
	// Append the inner ends b to the outer ends a
	[[gnu::always_inline]] void merge(bot_ends_t& a, const bot_ends_t& b) const {
		for (int z = 0; z < 2; z++) {
			// If there's no bottom edges, then we must be an isolated vertex, so we can end early.
			if (b[z][0] == -1) continue;
			if (a[z][0] == -1) {
				a[z] = b[z];
				continue;
			}
			link(a[z][1], b[z][0]);
			a[z][1] = b[z][1];
		}
	}
	[[gnu::always_inline]] void flip(bot_ends_t& a) const {
		std::swap(a[0], a[1]);
	}
	// Close off the returns tops of one side into its bottom
	[[gnu::always_inline]] void join_backedges(std::array<int, 2>& be, std::array<top_t, 2> tops) const {
		assert(be[1] != -1);
		link(be[1], tops[1].end);
		be[1] = tops[0].end;
	}
	// Fold side 1 onto side 0; side 1 returns only to lowval
	[[gnu::always_inline]] void join_bottoms(bot_ends_t& be, std::array<std::array<top_t, 2>, 2>& tops, [[maybe_unused]] int lowval) const {
		link(be[0][0], be[1][0]);
		be[0][0] = be[1][1];
		if (tops[1][0].end != -1) {
			// Caller must have checked that we're planar
			assert(tops[1][1].depth == lowval);
			// This is always true
			assert(tops[1][0].depth == lowval);
			link(tops[0][0].end, tops[1][0].end);
			tops[0][0].end = tops[1][1].end;
			// Already true since the backedge was on side 0
			assert(tops[0][0].depth == lowval);
		}
		be[1] = {-1, -1};
	}
	// The innermost return top of one side is finished: close it off into the bottom and return the next return outward
	[[gnu::always_inline]] top_t pop_top(std::array<int, 2>& be, top_t top) const {
		assert(be[1] != -1);
		link(be[1], top.end);
		be[1] = top.end ^ 1;
		top_t nxt;
		nxt.end = std::exchange(matches[be[1]], -1);
		if (nxt.end != -1) {
			matches[nxt.end] = -1;
			nxt.depth = edge_top_depths[nxt.end >> 2];
		}
		return nxt;
	}
};

// A (merged) ear on the tstack. The walk owns these fields; the mode's Payload rides along.
template <typename Payload, bool CHECK, bool STOP>
struct spqr_walk_entry_t {
	using top_t = spqr_walk_top_t;
	// Depth of the shallowest return
	int top_depth = -1;
	// Index of the first edge pushed inside this entry
	int first_idx = -1;
	// Bottom vertex
	int v_start = -1;
	// CHECK: per side, the outer/innermost exposed returns, depths increasing going inwards.
	// The convention is that tops[0][0].depth == top_depth, i.e. at least one minimal return lives on side 0.
	[[no_unique_address]] std::conditional_t<CHECK, std::array<std::array<top_t, 2>, 2>, std::monostate> tops{};
	// With STOP every entry is planar (we would have stopped otherwise)
	[[no_unique_address]] std::conditional_t<CHECK && !STOP, bool, std::monostate> nonplanar{};
	[[no_unique_address]] Payload payload{};

	[[nodiscard]] bool is_planar() const {
		if constexpr (CHECK && !STOP) return !nonplanar;
		else return true;
	}
	void mark_nonplanar() {
		if constexpr (CHECK && !STOP) nonplanar = true;
		else assert(false);
	}
	// Forget all returns; the entry is planar again
	void reset_returns() {
		tops = {};
		nonplanar = {};
	}
};

template <
	typename EnterVert, typename StartEdge, typename VertEntry, typename EdgeEntry,
	typename Merge, typename LinkTops, typename Flip, typename JoinBackedges, typename JoinBottoms, typename PopTop,
	typename AbsorbVert, typename BeginNode, typename FinishNode, typename FoldChild, typename FinishBlock, typename FinishRoot
>
struct spqr_walk_hooks_t {
	// enter_vert(v, depth): v is the dfs vertex at depth (called before its entry is pushed)
	EnterVert enter_vert;
	// start_edge(cur_depth, lowval, is_tree): the vertex at cur_depth is about to process its next outedge
	StartEdge start_edge;
	// vert_entry(t, v): initialize the payload of the fresh entry t for vertex v
	VertEntry vert_entry;
	// edge_entry(t, e_side, is_tree, is_block_root): initialize the payload (and for a backedge, t.tops) of the fresh entry t
	EdgeEntry edge_entry;
	// merge(a, b): the inner entry b is being merged into the outer entry a (the walk merges tops)
	Merge merge;
	// CHECK: link_tops(outer, inner): the return inner now sits directly inside the return outer
	LinkTops link_tops;
	// flip(t, end_idx): the sides of t are being swapped; [t.first_idx, end_idx) are the edge indices inside t
	Flip flip;
	// CHECK: join_backedges(t, z): the returns t.tops[z], all to the current depth, are closed off into the bottom of side z
	JoinBackedges join_backedges;
	// CHECK: join_bottoms(t, lowval): fold side 1 of t (returning only to lowval) onto side 0; the walk clears t.tops[1] afterwards
	JoinBottoms join_bottoms;
	// CHECK: pop_top(t, z) -> top_t: t.tops[z][1] is finished; return the next return outward of it (end == -1 if none)
	PopTop pop_top;
	// absorb_vert(vert_depth, cur_depth): the vertex entry at vert_depth joins the S chain through the current edge
	AbsorbVert absorb_vert;
	// begin_node(nxt, type, is_tree) -> item: nxt and the entry above it are about to become one node
	BeginNode begin_node;
	// finish_node(t, item, is_tree): t is the finished node item; make it the entry of the node's cap virtual edge
	FinishNode finish_node;
	// fold_child(t, cur_depth): t is the finished child subtree, about to be merged into the current vertex's entry
	FoldChild fold_child;
	// finish_block(t, e, lowval, is_tree): t is the finished block rooted at edge e, about to be merged into the current vertex's entry
	FinishBlock finish_block;
	// finish_root(t): t is the last entry of a dfs tree
	FinishRoot finish_root;
};

// Phase 2 of the SPQR build / planarity test: the lowval-ordered dfs over lowval_storted_skeleton_t, maintaining the tstack of ears.
// The walk owns the dfs stack, the tstack skeleton and (CHECK) the per-side returns that detect nonplanarity;
// everything mode-specific lives in each entry's Payload and in the hooks.
// CHECK: detect nonplanarity. STOP: return false at the first nonplanarity instead of marking the entry and continuing.
template <typename Payload, bool CHECK, bool STOP, typename Hooks>
bool spqr_walk(
	int NV,
	const std::vector<std::array<int, 2>>& edges,
	std::span<const int> vert_order,
	std::span<const int> edge_order,
	const Hooks& hooks
) {
	static_assert(CHECK || !STOP);
	using entry_t = spqr_walk_entry_t<Payload, CHECK, STOP>;
	using node_type = spqr_tree::node_type;

	// std::min is by reference, which breaks some optimizations
	auto setmin = [](auto& a, auto b) { if (b < a) a = b; };

	int NE = int(edges.size());
	assert(int(vert_order.size()) <= NV);
	assert(int(edge_order.size()) <= NE);

	auto [roots, outedges] = lowval_storted_skeleton_t::build(NV, edges, vert_order, edge_order);

	int nxt_edge_idx = 0; // Counts edges pushed onto the tstack
	std::vector<int> first_occurrence(NV); // First backedge to this depth

	std::vector<entry_t> tstack; tstack.reserve(NV + NE);
	auto cur_tstack = [&] [[gnu::always_inline]] () -> entry_t& { return tstack.end()[-1]; };
	auto nxt_tstack = [&] [[gnu::always_inline]] () -> entry_t& { return tstack.end()[-2]; };

	auto push_tstack = [&] [[gnu::always_inline]] (int v_start, int top_depth) -> entry_t& {
		entry_t t;
		t.top_depth = top_depth;
		t.first_idx = nxt_edge_idx;
		t.v_start = v_start;
		tstack.push_back(t);
		return tstack.back();
	};
	auto push_vert_tstack = [&] [[gnu::always_inline]] (int v, int depth) -> void {
		hooks.vert_entry(push_tstack(v, depth), v);
	};
	auto push_edge_tstack = [&] [[gnu::always_inline]] (int v_start, int top_depth, int e_side, bool is_tree, bool is_block_root) -> int {
		hooks.edge_entry(push_tstack(v_start, top_depth), e_side, is_tree, is_block_root);
		return nxt_edge_idx++;
	};
	auto flip_tstack = [&] [[gnu::always_inline]] (int i) -> void {
		if constexpr (CHECK) {
			entry_t& t = tstack[i];
			hooks.flip(t, i + 1 == int(tstack.size()) ? nxt_edge_idx : tstack[i+1].first_idx);
			if (t.is_planar()) std::swap(t.tops[0], t.tops[1]);
		}
	};
	auto merge_tstack_tops = [&] [[gnu::always_inline]] () -> void {
		entry_t& a = nxt_tstack();
		const entry_t& b = cur_tstack();
		setmin(a.top_depth, b.top_depth);
		hooks.merge(a, b);
		if constexpr (CHECK) {
			if (!a.is_planar()) {
				// Stays nonplanar
			} else if (!b.is_planar()) {
				a.mark_nonplanar();
			} else {
				for (int z = 0; z < 2; z++) {
					auto& at = a.tops[z];
					const auto& bt = b.tops[z];
					if (bt[0].end == -1) {
						// Do nothing
					} else if (at[0].end == -1) {
						at = bt;
					} else {
						// Caller must check that we're planar
						assert(at[1].depth <= bt[0].depth);
						hooks.link_tops(at[1], bt[0]);
						at[1] = bt[1];
					}
				}
			}
		}
		tstack.pop_back();
	};
	// Merge all backedges returning to cur_depth into the component on top of the tstack
	auto join_backedges_to_top = [&] [[gnu::always_inline]] (int cur_depth) -> void {
		if constexpr (CHECK) {
			entry_t& t = cur_tstack();
			if (!t.is_planar()) return;
			for (int z = 0; z < 2; z++) {
				// in the self-loop case, side 1 is empty; otherwise, both sides have returns
				if (t.tops[z][1].end == -1) continue;
				assert(t.tops[z][0].depth == cur_depth);
				assert(t.tops[z][1].depth == cur_depth);
				hooks.join_backedges(t, z);
				t.tops[z] = {};
			}
		}
	};
	// Fold side 1 of the top of the tstack onto side 0; precondition: side 1 should be the lowval-only side
	auto join_bottoms_to_empty = [&] [[gnu::always_inline]] (int lowval) -> void {
		if constexpr (CHECK) {
			entry_t& t = cur_tstack();
			if (!t.is_planar()) return;
			hooks.join_bottoms(t, lowval);
			t.tops[1] = {};
		}
	};
	// Prune the finished returns to cur_depth off the top of the tstack
	auto pop_tops = [&] [[gnu::always_inline]] (int cur_depth) -> void {
		if constexpr (CHECK) {
			entry_t& t = cur_tstack();
			for (int z = 0; z < 2; z++) {
				auto& tops = t.tops[z];
				while (tops[1].depth == cur_depth) {
					tops[1] = hooks.pop_top(t, z);
					if (tops[1].end == -1) tops = {};
				}
			}
		}
	};
	// The single place nonplanarity is handled: with STOP the caller unwinds, otherwise the top of the tstack is marked.
	auto fail = [&] [[gnu::always_inline]] () -> bool {
		if constexpr (STOP) {
			return false;
		} else {
			cur_tstack().mark_nonplanar();
			return true;
		}
	};

	struct dfs_stack_t {
		int v;
		bool has_vert_tstack;
		int ch_idx;
		int ch_end;
		int orig_tstack;
	};
	std::vector<dfs_stack_t> stk; stk.reserve(NV);
	for (auto rt : roots) {
		auto push_vert = [&] [[gnu::always_inline]] (int cur) -> void {
			int cur_depth = int(stk.size());
			hooks.enter_vert(cur, cur_depth);

			int lo = outedges.bounds[cur];
			int hi = outedges.bounds[cur+1];
			bool has_vert_tstack;
			{
				// Find the first same-BCC edge, and check it's type 2 (has lowval2), if so it's the ear tstack and we defer pushing ourselves.
				int first_edge = lo;
				while (first_edge < hi && outedges.dat[first_edge].key.is_new_block()) first_edge++;
				if (first_edge < hi && outedges.dat[first_edge].key.is_type_2()) {
					// Move first_edge to the beginning
					auto e = outedges.dat[first_edge];
					std::move_backward(outedges.dat.begin() + lo, outedges.dat.begin() + first_edge, outedges.dat.begin() + first_edge + 1);
					outedges.dat[lo] = e;
					has_vert_tstack = false;
				} else {
					push_vert_tstack(cur, cur_depth);
					has_vert_tstack = true;
				}
			}
			stk.push_back({cur, has_vert_tstack, lo, hi, -1});
		};
		// Returns the child to descend into for a tree edge
		auto start_edge = [&] [[gnu::always_inline]] () -> std::optional<int> {
			int cur_depth = int(stk.size()) - 1;
			auto& s = stk.back();
			assert(s.ch_idx < s.ch_end);
			auto [_, nxt, e_side, key] = outedges.dat[s.ch_idx];
			auto [lowval, is_tree, is_type_1] = key.unpack(cur_depth);

			if (lowval >= cur_depth || is_type_1) assert(s.has_vert_tstack);

			hooks.start_edge(cur_depth, lowval, is_tree);

			s.orig_tstack = int(tstack.size());
			if (is_tree) {
				first_occurrence[cur_depth] = NE;
				return nxt;
			} else {
				return std::nullopt;
			}
		};
		// Returns false iff we stopped at a nonplanarity (only with STOP)
		auto finish_edge = [&] [[gnu::always_inline]] () -> bool {
			int cur_depth = int(stk.size()) - 1;
			auto& s = stk.back();
			int cur = s.v;
			assert(s.ch_idx < s.ch_end);

			auto [_, nxt, e_side, key] = outedges.dat[s.ch_idx];
			int e = e_side >> 1;
			s.ch_idx++;

			auto [lowval, is_tree, is_type_1] = key.unpack(cur_depth);

			const int orig_tstack = s.orig_tstack;

			if (lowval >= cur_depth) {
				// A new block hanging off cur, rooted at this edge.
				if (is_tree) {
					push_edge_tstack(nxt, cur_depth, e_side, true, true);
					if (lowval == cur_depth) {
						// Merge the backedge
						merge_tstack_tops();
						join_backedges_to_top(cur_depth);
					}
					// Merge the vertex
					merge_tstack_tops();
					join_bottoms_to_empty(lowval);
				} else {
					assert(nxt == cur);
					push_edge_tstack(cur, lowval, e_side, false, true);
					join_backedges_to_top(cur_depth);
				}
				hooks.finish_block(cur_tstack(), e, lowval, is_tree);
				assert(s.has_vert_tstack);
				// Merge into the vertex tstack
				merge_tstack_tops();
				return true;
			}
			assert(lowval < cur_depth);

			if (is_tree) {
				push_edge_tstack(nxt, cur_depth, e_side, true, false);
				while (nxt_tstack().top_depth >= cur_depth) {
					node_type type = node_type::R;
					if (nxt_tstack().top_depth > cur_depth) {
						// This is a vertex in the tstack, followed by either an S edge, possibly merged with other things
						if (tstack.end()[-3].top_depth < cur_depth) {
							// Not actually a good return, just stop
							break;
						}
						hooks.absorb_vert(nxt_tstack().top_depth, cur_depth);
						// Merge the vertex in
						merge_tstack_tops();
						type = nxt_tstack().top_depth > cur_depth ? node_type::S : node_type::R;
					} else {
						// Same endpoints: this will be a P node
						if (nxt_tstack().v_start == cur_tstack().v_start) type = node_type::P;
					}
					int item = hooks.begin_node(nxt_tstack(), type, type == node_type::S);
					merge_tstack_tops();
					join_backedges_to_top(cur_depth);
					hooks.finish_node(cur_tstack(), item, true);
				}

				if (cur_tstack().first_idx > first_occurrence[cur_depth]) {
					// There are returns to cur_depth under us: check they're laminar / one-sided and merge them in.
					bool ok = [&] [[gnu::always_inline]] () -> bool {
						if constexpr (!CHECK) {
							return true;
						} else {
							int source = int(tstack.size()) - 1;
							do {
								--source;
								if (!tstack[source].is_planar()) {
									cur_tstack().mark_nonplanar();
									return true;
								}
							} while (tstack[source].first_idx > first_occurrence[cur_depth]);

							// From planarity's perspective, we can view each tstack as one of 2 shapes:
							// * tstack[i] can be a single "atom" branching off tstack[i].v_start. It can be:
							//   * A single backedge (type 1)
							//   * A subtree (possibly not biconnected) with at least 2 different-depth backedges on its outside (type 2)
							// * tstack[i] can be a "chunk". Chunks contain:
							//   * A core spanning from tstack[i].v_start (bot[0]) to tstack[i+1].v_start (bot[1]).
							//     * The core is a cyclic outer face: it has 2 *disjoint* paths from bot[0] to bot[1]
							//     * Each path can have backedges from its interior (not bot[0] or bot[1])
							//     * side[0]'s path has a backedge to tstack[i].lowval
							//   * Extra atoms on side 1, anchored at tstack[i].v_start
							//     * these must have lowval < tstack[i].lowval, or can have lowval == tstack[i].lowval and be type 2
							//     * these extra atoms are an entire suffix of v_start's: a chunk will always eat them all
							//   * Most of the time, we can treat the side[1] core and the atoms all as separate backedges from tstack[i].v_start.
							//     * The exception is when side[1].tops[1].depth == lowval: then it's forced to be a type 2 atom or part of the core, which matters.
							//       * TODO: Can we easily distinguish the 2 cases?
							// * tstack[i] can also be a tree vertex or a tree edge (trivial cases)
							//
							// Note that all atoms (including the chunk-extras) at one v_start must be sorted by (lowval, type).
							// However, a chunk can occur later (closer to the top) than its (lowval, type) sort at its v_start.

							// last_top == cur_tstack().top_depth
							int last_top = cur_depth;
							while (int(tstack.size()) > source + 2) {
								if (nxt_tstack().top_depth > cur_depth) {
									// Vertex or tree edge, no conditions
								} else if (nxt_tstack().top_depth == cur_depth) {
									if (nxt_tstack().tops[1][0].depth != -1) {
										// Double-sided to cur_depth, conflicts with cur_tstack()
										assert(last_top < cur_depth);
										// nxt_tstack() is a chunk and both backedges are on the core
										// K33 is:
										// * cur_tstack().tops[0]
										// * nxt_tstack().sides[0].tops[0].base
										// * nxt_tstack().sides[1].tops[0].base
										// + cur
										// + nxt_tstack().v_start
										// + cur_tstack().v_start
										//
										// nxt_tstack().v_start -> cur_tstack().tops[0] is the ear lowval loop
										// cur_tstack().v_start -> cur_tstack().tops[0] is just along cur
										// cur -> nxt_tstack().sides[*].tops[0].base is just the backedge
										// The rest is the outer face of nxt_tstack()
										return fail();
									}
									// We will put cur_depth on side 1 until the bottom
									flip_tstack(int(tstack.size()) - 2);
								} else {
									const auto& nt = nxt_tstack().tops;
									if (nt[1][0].depth != -1 && nt[1][0].depth != cur_depth) {
										// Non-empty on both sides, conflicts with source
										// nxt_stack() is a chunk
										if (nt[1][0].depth == nxt_tstack().top_depth) {
											// if it's core + type-2-atom
											// by the atom ordering, we're guaranteed source->cur isn't from v_start
											// * nxt_tstack().sides[0].tops[0].base
											// * nxt_tstack().sides[1].tops[0].base_fork
											// * source.base
											// + cur
											// + nxt_tstack().v_start
											// + nxt_tstack().top_depth
											//
											// cut nxt_tstack().core.side[1]
											// use nxt_tstack().sides[1].tops[0].prev to get from the fork to above cur down to cur
											//
											// otherwise it's double core
											// * cur
											// * nxt_tstack().sides[0].tops[0].base
											// * nxt_tstack().sides[1].tops[0].base
											// + cur_tstack().v_start
											// + nxt_tstack().v_start
											// + nxt_tstack().top_depth
											// (cut the lowval ear edge)
											//
											return fail();
										} else {
											// If it's an atom
											// by atom ordering, we're guaranteed source->cur isn't from v_start
											// * nxt_tstack().sides[0].tops[0].base
											// * nxt_tstack().sides[1].tops[0].end (go up/down to cur/top_depth)
											// * source.base
											// + cur
											// + nxt_tstack().v_start
											// + nxt_tstack().top_depth
											//
											// otherwise it's double core
											// * cur
											// * nxt_tstack().sides[0].tops[0].base
											// * nxt_tstack().sides[1].tops[0].base
											// + cur_tstack().v_start
											// + nxt_tstack().v_start
											// + nxt_tstack().sides[1].tops[0].end (side 0 gets there from above)
											// (cut the lowval ear edge)
											return fail();
										}
									}
									// Implicitly excludes -1
									if (nt[0][1].depth > last_top) {
										assert(last_top < cur_depth);
										// 3 conflicting edges with nxt_tstack(), cur_tstack(), and source
										return fail();
									}
									last_top = nxt_tstack().top_depth;
								}
								merge_tstack_tops();
							}

							int t0 = nxt_tstack().tops[0][1].depth;
							int t1 = nxt_tstack().tops[1][1].depth;
							assert(t0 == cur_depth || t1 == cur_depth);
							if (std::min(t0, t1) > last_top) {
								assert(last_top < cur_depth);
								return fail();
							}
							if (t0 == cur_depth) {
								// We need to flip cur_tstack and nxt_tstack relative to each other.
								// Flip the one with worse top_depth.
								flip_tstack(cur_tstack().top_depth < nxt_tstack().top_depth ? int(tstack.size()) - 2 : int(tstack.size()) - 1);
							}
							merge_tstack_tops();

							pop_tops(cur_depth);
							return true;
						}
					}();
					if (!ok) return false;
					// Without STOP, catch up on the merges skipped by a nonplanarity
					while (cur_tstack().first_idx > first_occurrence[cur_depth]) {
						merge_tstack_tops();
					}
				}

				if (is_type_1) assert(s.has_vert_tstack);
				if (s.has_vert_tstack) {
					// NB: tstack[orig_size] is the vertex and tstack[orig_size+1] is the backedge; maybe we should reverse them?
					assert(int(tstack.size()) >= orig_tstack + 3);

					if (!is_type_1) {
						bool ok = [&] [[gnu::always_inline]] () -> bool {
							if constexpr (!CHECK) {
								return true;
							} else {
								// The lowval side should be side 1, everything else goes on side 0.
								// The exception is tstack[orig_tstack + 2], which could be == lowval on one/both sides,
								// but is guaranteed to have *something* > lowval by non-type-1-ness
								auto& t = tstack[orig_tstack + 2];
								if (!t.is_planar()) {
									cur_tstack().mark_nonplanar();
									return true;
								}
								{
									const auto& tt = t.tops;
									assert(tt[0][0].depth == t.top_depth);
									assert(tt[0][1].depth != -1);
									if (tt[0][1].depth == lowval) {
										flip_tstack(orig_tstack + 2);
									} else if (tt[1][1].depth != -1 && tt[1][1].depth != lowval) {
										return fail();
									}
									assert(tt[0][1].depth > lowval);
								}
								int last_top = t.tops[0][1].depth;
								for (int i = orig_tstack + 3; i < int(tstack.size()); i++) {
									if (!tstack[i].is_planar()) {
										cur_tstack().mark_nonplanar();
										return true;
									}
									if (tstack[i].top_depth == lowval) {
										flip_tstack(i);
									}
									const auto& it = tstack[i].tops;
									if (it[1][1].depth != -1 && it[1][1].depth != lowval) {
										return fail();
									}
									int next_top = it[0][0].depth;
									if (next_top != -1) {
										if (last_top > next_top) {
											return fail();
										}
										last_top = it[0][1].depth;
									}
								}
								return true;
							}
						}();
						if (!ok) return false;
						while (int(tstack.size()) > orig_tstack + 3) {
							merge_tstack_tops();
						}
					}

					assert(int(tstack.size()) == orig_tstack + 3);
					int item = -1;
					if (is_type_1) {
						item = hooks.begin_node(nxt_tstack(), cur_tstack().top_depth == cur_depth ? node_type::S : node_type::R, false);
					}
					// Merge with the backedge
					merge_tstack_tops();
					// Merge with the vertex
					merge_tstack_tops();

					cur_tstack().v_start = cur;
					assert(cur_tstack().top_depth == lowval);

					// Now that we're leaving the child, fold everything to the correct side.
					hooks.fold_child(cur_tstack(), cur_depth);
					join_bottoms_to_empty(lowval);

					if (is_type_1) {
						hooks.finish_node(cur_tstack(), item, false);
					}
				}
			} else {
				assert(is_type_1);
				int idx = push_edge_tstack(cur, lowval, e_side, false, false);
				setmin(first_occurrence[lowval], idx);
			}

			assert(int(tstack.size()) >= orig_tstack + 1);

			// If is_type_1, the last entry on the tstack is either the vert_tstack, or the previous child as a unit
			if (is_type_1 && nxt_tstack().top_depth == lowval) {
				assert(s.has_vert_tstack);
				// This will be a P node
				int item = hooks.begin_node(nxt_tstack(), node_type::P, false);
				merge_tstack_tops();
				hooks.finish_node(cur_tstack(), item, false);
			}

			if (!s.has_vert_tstack) {
				assert(!is_type_1);
				// Throw cur_vert_node onto the tstack so it'll get interleaved correctly
				push_vert_tstack(cur, cur_depth);
				s.has_vert_tstack = true;
			}
			return true;
		};
		auto pop_vert = [&] [[gnu::always_inline]] () -> void {
			auto& s = stk.back();
			assert(s.ch_idx == s.ch_end);
			assert(s.has_vert_tstack);
			stk.pop_back();
		};

		push_vert(rt);
		while (true) {
			if (stk.back().ch_idx == stk.back().ch_end) {
				pop_vert();
				if (stk.empty()) break;
				if (!finish_edge()) return false;
			} else if (std::optional<int> nxt = start_edge(); nxt) {
				push_vert(*nxt);
			} else {
				// Backedges can't fail
				[[maybe_unused]] bool ok = finish_edge();
				assert(ok);
			}
		}
		assert(int(tstack.size()) == 1);
		hooks.finish_root(tstack.back());
		tstack.pop_back();
	}
	// Every edge gets a tstack entry
	assert(nxt_edge_idx == NE);
	return true;
}


template <bool with_planarity>
std::conditional_t<with_planarity, planar_spqr_tree, spqr_tree> spqr_tree::build_impl(
	int NV,
	const std::vector<std::array<int, 2>>& edges,
	bool ternarize,
	std::span<const int> vert_order,
	std::span<const int> edge_order
) {
	int NE = int(edges.size());
	assert(int(vert_order.size()) <= NV);
	assert(int(edge_order.size()) <= NE);

	// We're going to build a tree of all SPQR *nodes* + all original *vertices* (collectively *items*).
	// Vertices will hang off the first SPQR node containing them, and blocks will be rooted at a topmost Q node for the top edge.
	// Items are 1 + v for vertices, 1 + NV + e for edges, and then the allocated SPQR nodes.
	constexpr int ROOT_ITEM = 0;
	auto vert_item = [&] [[gnu::always_inline]] (int v) -> int { return 1 + v; };
	auto edge_item = [&] [[gnu::always_inline]] (int e) -> int { return 1 + NV + e; };

	// As we build, we will represent the children of our nodes/vertices as linked lists.
	struct item_list {
		// Items are actually 2 * item + planarity_flip (always 0 without planarity)
		std::array<int, 2> v{-1, -1};

		[[nodiscard]] bool empty() const { return v[0] < 0; }
	};
	std::vector<int> ch_nxt; ch_nxt.reserve(1 + NV + NE + NE); ch_nxt.assign(1 + NV + NE, -1);
	std::vector<std::array<int, 2>> item_vs; item_vs.reserve(1 + NV + 2 * NE); item_vs.resize(1 + NV + NE, {-1, -1});
	std::vector<item_list> item_ch; item_ch.reserve(1 + NV + 2 * NE); item_ch.resize(1 + NV + NE, item_list{});
	std::vector<node_type> item_types; item_types.reserve(1 + NV + 2 * NE);
	item_types.resize(1, node_type::F);
	item_types.resize(1 + NV, node_type::V);
	item_types.resize(1 + NV + NE, node_type::Q);
	auto concat = [&] [[gnu::always_inline]] (item_list a, item_list b) -> item_list {
		if (b.empty()) return a;
		if (a.empty()) return b;
		ch_nxt[a.v[1] >> 1] = b.v[0] ^ (a.v[1] & 1);
		return {{a.v[0], b.v[1]}};
	};
	auto unit_list = [&] [[gnu::always_inline]] (int item) -> item_list {
		return {{item << 1, item << 1}};
	};

	struct nonplanarity_certificate_t {
		// TODO: What's the nonplanarity certificate look like?
	};
	// with_planarity: per-node quarter edge matchings of the cap vedge
	std::vector<std::expected<std::array<int, 4>, nonplanarity_certificate_t>> node_planarity;
	if constexpr (with_planarity) node_planarity.reserve(NE);
	// with_planarity: each vedge has 4 quarter edges by 4 * vedge_id + 2 * source_vert + is_cw (is_cw is arbitrary);
	// vedges are identified with what item they cap, numbered by (item - 1 - NV), plus a scratch vedge 2 * NE for phase 3.
	std::vector<int> quarter_edge_matches(with_planarity ? 8 * NE + 4 : 0, -1);
	int tot_blocks = 0;
	int tot_self_loops = 0;

	// Phase 2
	{
		// Helpers for working with std::array<T, 2> - these compile to cmov's better than direct index access.

		// return arr[dir] == a, arr[!dir] == b
		auto set_sides = []<typename T>(bool dir, T a, T b) -> std::array<T, 2> {
			return dir ? std::array<T, 2>{b, a} : std::array<T, 2>{a, b};
		};
		auto get_side = []<typename T>(std::array<T, 2> a, bool dir) -> T {
			return dir ? a[1] : a[0];
		};

		auto alloc_item = [&] [[gnu::always_inline]] (node_type type) -> int {
			int item = int(item_vs.size());
			item_vs.push_back({});
			item_ch.push_back({});
			item_types.push_back(type);
			ch_nxt.push_back(-1);
			if constexpr (with_planarity) node_planarity.emplace_back();
			return item;
		};

		// Most of our code will be in terms of v_start / top_depth, so we'll want to read these out
		std::vector<int> stack_verts(NV);
		std::vector<int8_t> stack_dir(NV); // really bool, but I don't want vector<bool>
		auto make_vs = [&] [[gnu::always_inline]] (int v_start, int top_depth) -> std::array<int, 2> {
			return set_sides(bool(stack_dir[top_depth]), stack_verts[top_depth], v_start);
		};

		using top_t = spqr_walk_top_t;
		using qe_t = spqr_walk_quarter_edges_t;
		std::vector<int> edge_top_depths(with_planarity ? 2 * NE : 0, -1);
		[[maybe_unused]] qe_t qe{quarter_edge_matches, edge_top_depths};

		struct payload_t {
			std::array<item_list, 2> spans;
			[[no_unique_address]] std::conditional_t<with_planarity, qe_t::bot_ends_t, std::monostate> bot_ends;
		};

		// Quarter edges of the (v)edge item, oriented by stack_dir[top_depth]
		auto vedge_planarity = [&] [[gnu::always_inline]] (auto& t, int item, bool is_tree) -> void {
			if constexpr (with_planarity) {
				int ve = item - (1 + NV);
				bool top_dir = stack_dir[t.top_depth];
				edge_top_depths[ve] = t.top_depth;
				auto& be = t.payload.bot_ends;
				if (is_tree) {
					be[0] = {4 * ve + 2 * !top_dir + 0, 4 * ve + 2 * top_dir + 1};
					be[1] = {4 * ve + 2 * !top_dir + 1, 4 * ve + 2 * top_dir + 0};
				} else {
					be[0] = {4 * ve + 2 * !top_dir + 0, 4 * ve + 2 * !top_dir + 1};
					t.tops[0] = {{{4 * ve + 2 * top_dir + 1, t.top_depth}, {4 * ve + 2 * top_dir + 0, t.top_depth}}};
				}
			}
		};
		// t becomes the entry of the (v)edge item
		auto vedge_entry = [&] [[gnu::always_inline]] (auto& t, int item, bool is_tree) -> void {
			t.payload.spans = set_sides(bool(stack_dir[t.top_depth]), unit_list(item), {});
			vedge_planarity(t, item, is_tree);
		};

		spqr_walk_hooks_t hooks{
			.enter_vert = [&] [[gnu::always_inline]] (int v, int depth) -> void {
				stack_verts[depth] = v;
				// Set something arbitrary, this is the normal convention for block-roots
				if (depth == 0) stack_dir[0] = true;
			},
			.start_edge = [&] [[gnu::always_inline]] (int cur_depth, int lowval, bool is_tree) -> void {
				// edge_dir convention: false is forwards, true is backwards.
				// That means that cur is on the edge_dir side and nxt is on the !edge_dir side.
				stack_dir[cur_depth] = (lowval >= cur_depth ? false : !stack_dir[lowval]);
				if (is_tree) stack_dir[cur_depth+1] = !stack_dir[std::min(lowval, cur_depth)];
			},
			.vert_entry = [&] [[gnu::always_inline]] (auto& t, int v) -> void {
				t.payload.spans = set_sides(bool(stack_dir[t.top_depth]), unit_list(vert_item(v)), {});
			},
			.edge_entry = [&] [[gnu::always_inline]] (auto& t, int e_side, bool is_tree, bool is_block_root) -> void {
				int item = edge_item(e_side >> 1);
				if (is_block_root) {
					// The block root is a Q item hanging off the vertex; its entry only collects the block's children.
					// Its quarter edges are scratch: nothing reads a block root's matching.
					vedge_planarity(t, item, is_tree);
				} else {
					item_vs[item] = make_vs(t.v_start, t.top_depth);
					vedge_entry(t, item, is_tree);
				}
			},
			.merge = [&] [[gnu::always_inline]] (auto& a, const auto& b) -> void {
				a.payload.spans[0] = concat(b.payload.spans[0], a.payload.spans[0]);
				a.payload.spans[1] = concat(a.payload.spans[1], b.payload.spans[1]);
				if constexpr (with_planarity) {
					if (a.is_planar() && b.is_planar()) qe.merge(a.payload.bot_ends, b.payload.bot_ends);
				}
			},
			.link_tops = [&] [[gnu::always_inline]] (top_t outer, top_t inner) -> void {
				qe.link(outer.end, inner.end);
			},
			.flip = [&] [[gnu::always_inline]] (auto& t, int) -> void {
				for (auto& span : t.payload.spans) {
					span.v[0] ^= 1;
					span.v[1] ^= 1;
				}
				if (t.is_planar()) qe.flip(t.payload.bot_ends);
			},
			.join_backedges = [&] [[gnu::always_inline]] (auto& t, int z) -> void {
				qe.join_backedges(t.payload.bot_ends[z], t.tops[z]);
			},
			.join_bottoms = [&] [[gnu::always_inline]] (auto& t, int lowval) -> void {
				qe.join_bottoms(t.payload.bot_ends, t.tops, lowval);
			},
			.pop_top = [&] [[gnu::always_inline]] (auto& t, int z) -> top_t {
				return qe.pop_top(t.payload.bot_ends[z], t.tops[z][1]);
			},
			.absorb_vert = [&] [[gnu::always_inline]] (int vert_depth, int cur_depth) -> void {
				// Just backfill this for begin_node
				stack_dir[vert_depth] = stack_dir[cur_depth];
			},
			.begin_node = [&] [[gnu::always_inline]] (auto& t, node_type type, bool is_tree) -> int {
				if (type == node_type::R) return alloc_item(type);

				assert(type == node_type::P || type == node_type::S);

				// If we want to ternarize, never reuse.
				if (ternarize) return alloc_item(type);

				bool top_dir = stack_dir[t.top_depth];
				assert(get_side(t.payload.spans, !top_dir).empty());
				int item = get_side(t.payload.spans, top_dir).v[0] >> 1;
				assert(item == (get_side(t.payload.spans, top_dir).v[1] >> 1));
				if (item_types[item] != type) return alloc_item(type);

				t.payload.spans = set_sides(top_dir, item_ch[item], {});
				if constexpr (with_planarity) {
					// Unwrap the planarity data
					// We don't really need to maintain this at all because S/P nodes are known to be trivially planar
					// The current state is just vedge_planarity(item), which means that it has the right shape, just needs to be relabelled.
					assert(node_planarity[item - (1 + NV + NE)]);
					const auto& matches = *node_planarity[item - (1 + NV + NE)];
					assert(t.is_planar());
					auto& be = t.payload.bot_ends;
					if (is_tree) {
						be[0][0] = matches[2 * !top_dir + 1];
						be[0][1] = matches[2 * top_dir + 0];
						be[1][0] = matches[2 * !top_dir + 0];
						be[1][1] = matches[2 * top_dir + 1];
					} else {
						be[0][0] = matches[2 * !top_dir + 1];
						be[0][1] = matches[2 * !top_dir + 0];
						t.tops[0][0].end = matches[2 * top_dir + 0];
						t.tops[0][1].end = matches[2 * top_dir + 1];
					}
				}
				return item;
			},
			.finish_node = [&] [[gnu::always_inline]] (auto& t, int item, bool is_tree) -> void {
				bool top_dir = stack_dir[t.top_depth];
				assert(get_side(t.payload.spans, !top_dir).empty());

				if constexpr (with_planarity) {
					if (t.is_planar()) {
						const auto& be = t.payload.bot_ends;
						std::array<int, 4> matches{};
						if (is_tree) {
							matches[2 * !top_dir + 1] = be[0][0];
							matches[2 * top_dir + 0] = be[0][1];
							matches[2 * !top_dir + 0] = be[1][0];
							matches[2 * top_dir + 1] = be[1][1];
						} else {
							matches[2 * !top_dir + 1] = be[0][0];
							matches[2 * !top_dir + 0] = be[0][1];
							matches[2 * top_dir + 0] = t.tops[0][0].end;
							matches[2 * top_dir + 1] = t.tops[0][1].end;
						}
						node_planarity[item - (1 + NV + NE)] = matches;
					} else {
						assert(item_types[item] == node_type::R);
						node_planarity[item - (1 + NV + NE)] = std::unexpected(nonplanarity_certificate_t{});
					}
				}
				item_vs[item] = make_vs(t.v_start, t.top_depth);
				item_ch[item] = get_side(t.payload.spans, top_dir);

				// t becomes the node's cap vedge
				t.reset_returns();
				t.payload = {};
				vedge_entry(t, item, is_tree);
			},
			.fold_child = [&] [[gnu::always_inline]] (auto& t, int cur_depth) -> void {
				// The entire subtree should go to the !edge_dir side.
				bool edge_dir = stack_dir[cur_depth];
				t.payload.spans = set_sides(!edge_dir, concat(t.payload.spans[0], t.payload.spans[1]), {});
			},
			.finish_block = [&] [[gnu::always_inline]] (auto& t, int e, int lowval, bool is_tree) -> void {
				int cur_depth = t.top_depth;
				int cur = stack_verts[cur_depth];
				// There's no planarity handling for this because it's just a Q node. I/O nodes also don't need any tracking.
				item_vs[edge_item(e)] = {cur, -1};
				tot_blocks++;
				item_list ch = concat(t.payload.spans[0], t.payload.spans[1]);
				if (is_tree) {
					if (lowval == cur_depth + 1) {
						// Bridge: prepend the I node
						int item = alloc_item(node_type::I);
						item_vs[item] = make_vs(t.v_start, cur_depth);
						ch = concat(unit_list(item), ch);
					}
				} else {
					// Self loop
					assert(ch.empty());
					tot_self_loops++;
					int item = alloc_item(node_type::O);
					// Make sure the nxt is -1 as well
					item_vs[item] = {cur, -1};
					ch = unit_list(item);
				}
				item_ch[edge_item(e)] = ch;
				item_ch[vert_item(cur)] = concat(item_ch[vert_item(cur)], unit_list(edge_item(e)));

				// The block is done; nothing of it merges into the vertex's entry
				t.reset_returns();
				t.payload = {};
			},
			.finish_root = [&] [[gnu::always_inline]] (auto& t) -> void {
				item_ch[ROOT_ITEM] = concat(item_ch[ROOT_ITEM], t.payload.spans[1]);
			},
		};
		spqr_walk<payload_t, with_planarity, false>(NV, edges, vert_order, edge_order, hooks);
	}

	// Phase 3: relabel the full tree in preorder
	int tot_items = int(item_types.size());
	{
		std::vector<int> vert_index(NV, -1);
		std::vector<int> edge_index(NE, -1);
		std::vector<bool> edge_flipped(NE);

		std::vector<int> par(tot_items, -1);
		std::vector<int> subtree_end(tot_items, -1);
		std::vector<node_type> types(tot_items, node_type::F);
		std::vector<int> orig_id(tot_items, -1);

		csr<int> ch;
		ch.bounds.resize(tot_items + 1, 0);
		ch.dat.resize(tot_items - 1);

		// Each node is a child, and additionally most non-block node has 2 cap verts; blocks have 1, and O nodes have 1
		int tot_node_verts = NV + (tot_items - 1 - NV) * 2 - tot_blocks - tot_self_loops;
		std::vector<node_vert_t> node_verts(tot_node_verts);
		csr_index node_nvs; node_nvs.bounds.resize(tot_items + 1);
		std::vector<int> vert_par_nv(tot_items, -1);

		int tot_node_edges = (tot_items - 1 - NV - tot_blocks) * 2;
		std::vector<node_edge_t> node_edges(tot_node_edges);
		csr_index node_nes; node_nes.bounds.resize(tot_items + 1);

		csr<node_adj_t> node_adj;
		node_adj.bounds.resize(tot_node_verts * 2 + 1);
		node_adj.dat.resize(tot_node_edges * 2);

		std::vector<bool> node_planar(with_planarity ? tot_items : 0);
		std::vector<int> ne_rot_adj(with_planarity ? 4 * tot_node_edges : 0, -1);

		std::vector<int> vert_pos_buf(NV, -1);
		std::vector<int> cnts_buf(2 * NV, -1);
		struct ch_buf_t {
			int loc;
			int item_id;
		};
		std::vector<ch_buf_t> ch_buf(tot_items);
		std::vector<int> rot_edge_ne(with_planarity ? 2 * NE + 1 : 0);

		int nxt_unassigned_idx = 0;

		struct dfs_stack_t {
			int cur_idx;
			int ch_idx;
			int ch_end;
			int cur_nv;
			int cur_ne;
		};
		std::vector<dfs_stack_t> stk; stk.reserve(tot_items);
		auto push_item = [&] [[gnu::always_inline]] (int cur_item) -> void {
			int cur_idx = nxt_unassigned_idx++;
			node_type cur_type = types[cur_idx] = item_types[cur_item];
			bool planar = true;
			if (cur_type == node_type::F) {
				assert(cur_item == 0);
			} else if (cur_type == node_type::V) {
				assert(1 <= cur_item && cur_item < 1 + NV);
				int orig_vert = cur_item - 1;
				orig_id[cur_idx] = orig_vert;
				vert_index[orig_vert] = cur_idx;
			} else if (cur_type == node_type::Q) {
				assert(1 + NV <= cur_item && cur_item < 1 + NV + NE);
				int orig_edge = cur_item - 1 - NV;
				orig_id[cur_idx] = orig_edge;
				edge_index[orig_edge] = cur_idx;
				assert(item_vs[cur_item][0] != -1);
				edge_flipped[orig_edge] = item_vs[cur_item][0] != edges[orig_edge][0];
			} else {
				assert(1 + NV + NE <= cur_item);
				if constexpr (with_planarity) {
					if (cur_type == node_type::O || cur_type == node_type::I) {
						// No planarity data was set up
					} else if (cur_type == node_type::S || cur_type == node_type::P || cur_type == node_type::R) {
						const auto& p = node_planarity[cur_item - (1 + NV + NE)];
						if (p) {
							// Make sure this runs before our planarity_flip checks
							for (int s = 0; s < 4; s++) {
								int a = 8 * NE + s, b = (*p)[s];
								quarter_edge_matches[a] = b;
								quarter_edge_matches[b] = a;
							}
						} else {
							// TODO: Any certificate stuff
							planar = false;
						}
					} else assert(false);
				}
			}
			if constexpr (with_planarity) node_planar[cur_idx] = planar;

			// HACK: Fill ch and vert_items in with orig items / orig verts for now,
			// because we don't have the final item id's yet.
			int ch_st = ch.bounds[cur_idx];
			int ch_en = ch_st;
			int nv_st = node_nvs.bounds[cur_idx];
			int nv_en = nv_st;
			int n_edges = 0;
			if (item_vs[cur_item][0] != -1) {
				node_verts[nv_en++] = {cur_idx, item_vs[cur_item][0]};
			}
			if (!item_ch[cur_item].empty()) {
				bool planarity_flip = item_ch[cur_item].v[0] & 1;
				for (int ch_item = item_ch[cur_item].v[0] >> 1; true; planarity_flip ^= (ch_nxt[ch_item] & 1), ch_item = ch_nxt[ch_item] >> 1) {
					ch.dat[ch_en++] = ch_item;
					assert(ch_item >= 1);
					if (ch_item < 1 + NV) {
						node_verts[nv_en++] = {cur_idx, ch_item - 1};
					} else {
						if constexpr (with_planarity) {
							if (cur_type != node_type::R) {
								assert(!planarity_flip);
							} else {
								// Fix the planarity direction right here: reverse quarter_edge_matches upfront;
								// this breaks the involution property, but from here on we'll never read the low bits anyways.
								int ve = ch_item - (1 + NV);
								if (planarity_flip) {
									std::swap(quarter_edge_matches[4 * ve + 0], quarter_edge_matches[4 * ve + 1]);
									std::swap(quarter_edge_matches[4 * ve + 2], quarter_edge_matches[4 * ve + 3]);
								}
							}
						}
						n_edges++;
					}
					if (ch_item == (item_ch[cur_item].v[1] >> 1)) {
						assert(ch_nxt[ch_item] == -1);
						break;
					}
				}
				planarity_flip ^= item_ch[cur_item].v[1] & 1;
				assert(!planarity_flip);
			}
			if (item_vs[cur_item][1] != -1) {
				node_verts[nv_en++] = {cur_idx, item_vs[cur_item][1]};
			}
			ch.bounds[cur_idx+1] = ch_en;
			node_nvs.bounds[cur_idx+1] = nv_en;

			int n_verts = nv_en - nv_st;

			bool is_node = cur_type != node_type::F && cur_type != node_type::V;
			bool has_cap = is_node && !(cur_type == node_type::Q && ch_en - ch_st > 0);

			if (!is_node) n_edges = 0;
			if (has_cap) n_edges++;

			int ne_st = node_nes.bounds[cur_idx];
			int ne_en = node_nes.bounds[cur_idx+1] = ne_st + n_edges;

			auto set_ne = [&] [[gnu::always_inline]] (int ne, std::array<int, 2> nvs, std::array<int, 2> nds, std::array<int, 4> rot_adjs) -> void {
				node_edges[ne].node = cur_idx;
				node_edges[ne].nvs = nvs;
				node_adj.dat[nds[0]] = {ne, nvs[1]};
				node_adj.dat[nds[1]] = {ne, nvs[0]};
				if constexpr (with_planarity) {
					for (int z = 0; z < 4; z++) ne_rot_adj[4 * ne + z] = rot_adjs[z];
				}
			};
			if (cur_type == node_type::F) {
				// Just set node_adj bounds and we're good
				for (int i = 2 * nv_st+1; i <= 2 * nv_en; i++) {
					node_adj.bounds[i] = 2 * ne_st;
				}
			} else if (cur_type == node_type::V) {
				// Nothing to do
			} else if (n_verts == 1) {
				assert(cur_type == node_type::Q || cur_type == node_type::O);
				assert(n_edges == 1);
				node_adj.bounds[2 * nv_st + 1] = 2 * ne_st + 1 * n_edges;
				node_adj.bounds[2 * nv_st + 2] = 2 * ne_st + 2 * n_edges;
				set_ne(ne_st, {nv_st, nv_st}, {2 * ne_st + 1, 2 * ne_st}, {4 * ne_st + 3, 4 * ne_st + 2, 4 * ne_st + 1, 4 * ne_st + 0});
			} else if (cur_type == node_type::Q || cur_type == node_type::I) {
				assert(n_verts == 2);
				assert(n_edges == 1);
				node_adj.bounds[2 * nv_st + 1] = 2 * ne_st + 0 * n_edges;
				node_adj.bounds[2 * nv_st + 2] = 2 * ne_st + 1 * n_edges;
				node_adj.bounds[2 * nv_st + 3] = 2 * ne_st + 2 * n_edges;
				node_adj.bounds[2 * nv_st + 4] = 2 * ne_st + 2 * n_edges;
				set_ne(ne_st, {nv_st, nv_st + 1}, {2 * ne_st, 2 * ne_st + 1}, {4 * ne_st + 1, 4 * ne_st + 0, 4 * ne_st + 3, 4 * ne_st + 2});
			} else if (cur_type == node_type::P) {
				// Special case: tiebreak the parallel edges so they're reversed
				assert(n_verts == 2);
				assert(n_edges >= 3);
				node_adj.bounds[2 * nv_st + 1] = 2 * ne_st + 0 * n_edges;
				node_adj.bounds[2 * nv_st + 2] = 2 * ne_st + 1 * n_edges;
				node_adj.bounds[2 * nv_st + 3] = 2 * ne_st + 2 * n_edges;
				node_adj.bounds[2 * nv_st + 4] = 2 * ne_st + 2 * n_edges;
				for (int ne = ne_st; ne < ne_en; ne++) {
					int ne_prv = (ne == ne_st ? ne_en : ne) - 1;
					int ne_nxt = (ne+1 == ne_en ? ne_st : ne+1);
					std::array<int, 4> rot_adjs{4 * ne_prv + 1, 4 * ne_nxt + 0, 4 * ne_nxt + 3, 4 * ne_prv + 2};
					set_ne(ne, {nv_st, nv_st + 1}, {2 * ne_st + (ne - ne_st), 2 * ne_en - 1 - (ne - ne_st)}, rot_adjs);
				}
			} else if (cur_type == node_type::S) {
				assert(n_verts == n_edges);
				assert(n_verts >= 3);
				for (int i = 2 * nv_st + 1; i <= 2 * nv_en; i++) {
					node_adj.bounds[i] = i + 2 * (ne_st - nv_st);
				}
				// Fix bounds for the cap
				node_adj.bounds[2 * nv_st + 1]--;
				node_adj.bounds[2 * nv_en - 1]++;
				set_ne(ne_st, {nv_st, nv_en - 1}, {2 * ne_st, 2 * ne_en - 1}, {4 * (ne_st+1) + 1, 4 * (ne_st+1) + 0, 4 * (ne_en-1) + 3, 4 * (ne_en-1) + 2});
				for (int i = 1; i < n_edges; i++) {
					int ne = ne_st + i;
					std::array<int, 4> rot_adjs{4 * (ne-1) + 3, 4 * (ne-1) + 2, 4 * (ne+1) + 1, 4 * (ne+1) + 0};
					if (ne-1 == ne_st) { rot_adjs[0] = 4 * ne_st + 1, rot_adjs[1] = 4 * ne_st + 0; }
					if (ne+1 == ne_en) { rot_adjs[2] = 4 * ne_st + 3, rot_adjs[3] = 4 * ne_st + 2; }
					set_ne(ne, {nv_st + i - 1, nv_st + i}, {2 * ne - 1, 2 * ne}, rot_adjs);
				}
			} else if (cur_type == node_type::R) {
				// Bucketsort the children by the midpoint
				for (int nv = nv_st; nv < nv_en; nv++) {
					vert_pos_buf[node_verts[nv].vert] = nv;
				}
				cnts_buf.assign(n_verts * 2 - 1, 0);
				ch_buf.clear();

				assert(has_cap);

				// Cap node_adj bounds
				node_adj.bounds[2 * nv_st + 2]++;
				node_adj.bounds[2 * nv_en - 1]++;

				for (int i = ch_st; i < ch_en; i++) {
					int item = ch.dat[i];
					assert(item >= 1);
					std::array<int, 2> nvs;
					if (item < 1 + NV) {
						nvs = {vert_pos_buf[item-1], vert_pos_buf[item-1]};
					} else {
						nvs = {vert_pos_buf[item_vs[item][0]], vert_pos_buf[item_vs[item][1]]};
						assert(nvs[0] < nvs[1]);
						node_adj.bounds[2 * nvs[0] + 2]++;
						node_adj.bounds[2 * nvs[1] + 1]++;
					}
					int loc = (nvs[0] - nv_st) + (nvs[1] - nv_st);
					ch_buf.emplace_back(loc, item);
					cnts_buf[loc]++;
				}
				int offset = ch_st;
				for (auto& cnt : cnts_buf) {
					offset += cnt;
					cnt = offset;
				}
				for (auto [loc, n] : std::views::reverse(ch_buf)) {
					ch.dat[--cnts_buf[loc]] = n;
				}

				if constexpr (with_planarity) {
					// Set up the reverse mapping for ourselves
					int nxt_ne = ne_en;
					for (int i = ch_en - 1; i >= ch_st; i--) {
						int item = ch.dat[i];
						assert(item >= 1);
						if (item < 1 + NV) continue;
						nxt_ne--;
						rot_edge_ne[item - (1 + NV)] = nxt_ne;
					}
					assert(nxt_ne == ne_st + 1);
					rot_edge_ne[2 * NE] = ne_st;
				}
				auto map_rot_edge = [&] [[gnu::always_inline]] (int ve) -> std::array<int, 4> {
					if constexpr (!with_planarity) return {-1, -1, -1, -1};
					if (!planar) return {-1, -1, -1, -1};
					std::array<int, 4> res{};
					for (int z = 0; z < 4; z++) {
						int o = quarter_edge_matches[4 * ve + z];
						assert(o != -1);
						res[z] = (rot_edge_ne[o >> 2] << 2) + (o & 2) + !(z & 1);
					}
					return res;
				};

				{
					int off = 2 * ne_st;
					for (int i = 2 * nv_st + 1; i <= 2 * nv_en; i++) {
						off += std::exchange(node_adj.bounds[i], off);
					}
					assert(off == 2 * ne_en);
				}

				// Fill in node_edges and node_adj.
				// Reverse order to get the adj in bracket ordering.
				{
					// Handle cap as special: it's first in the node_edges, which means it's in the wrong place for the left endpoint.
					node_adj.bounds[2 * nv_st + 2]++;

					int nxt_ne = ne_en;
					for (int i = ch_en - 1; i >= ch_st; i--) {
						int item = ch.dat[i];
						assert(item >= 1);
						if (item < 1 + NV) continue;
						nxt_ne--;
						auto [v0, v1] = item_vs[item];
						// TODO: Reuse this from the ch pass?
						std::array<int, 2> nvs = {vert_pos_buf[v0], vert_pos_buf[v1]};
						set_ne(nxt_ne, nvs, {
							node_adj.bounds[2 * nvs[0] + 2]++,
							node_adj.bounds[2 * nvs[1] + 1]++,
						}, map_rot_edge(item - (1 + NV)));
					}
					assert(nxt_ne == ne_st + 1);

					// Insert the cap / bump its bound
					set_ne(ne_st, {nv_st, nv_en - 1}, {2 * ne_st, 2 * ne_en - 1}, map_rot_edge(2 * NE));
					node_adj.bounds[2 * nv_en - 1]++;
				}
			} else assert(false);

			int cur_nv = nv_st + (item_vs[cur_item][0] != -1);
			int cur_ne = ne_st + has_cap;
			stk.push_back({cur_idx, ch_st, ch_en, cur_nv, cur_ne});
		};

		auto start_child = [&] [[gnu::always_inline]] () -> int {
			auto& [cur_idx, ch_idx, ch_en, cur_nv, cur_ne] = stk.back();
			assert(ch_idx < ch_en);
			int nxt_item = ch.dat[ch_idx];
			int nxt_idx = nxt_unassigned_idx;
			ch.dat[ch_idx] = nxt_idx;
			par[nxt_idx] = cur_idx;
			int nxt_ne = node_nes.bounds[nxt_idx];
			if (nxt_item < 1 + NV) {
				vert_par_nv[nxt_idx] = cur_nv++;
			} else if (types[cur_idx] != node_type::F && types[cur_idx] != node_type::V) {
				node_edges[cur_ne].twin_ne = nxt_ne;
				node_edges[nxt_ne].twin_ne = cur_ne;
				cur_ne++;
			}

			ch_idx++;
			return nxt_item;
		};

		auto pop_item = [&] [[gnu::always_inline]] () -> void {
			auto [cur_idx, ch_idx, ch_en, cur_nv, cur_ne] = stk.back(); stk.pop_back();
			assert(ch_idx == ch_en);
			subtree_end[cur_idx] = nxt_unassigned_idx;
		};

		par[nxt_unassigned_idx] = -1;
		push_item(ROOT_ITEM);
		while (true) {
			if (stk.back().ch_idx == stk.back().ch_end) {
				pop_item();
				if (stk.empty()) break;
			} else {
				push_item(start_child());
			}
		}

		assert(nxt_unassigned_idx == tot_items);
		assert(ch.bounds.back() == int(ch.dat.size()));
		assert(node_nvs.bounds.back() == int(node_verts.size()));
		assert(node_nes.bounds.back() == int(node_edges.size()));
		assert(node_adj.bounds.back() == int(node_adj.dat.size()));

		// Rewrite node_vertices to the correct index
		for (auto& v : node_verts) {
			v.vert = vert_index[v.vert];
		}

		spqr_tree res{
			std::move(vert_index),
			std::move(edge_index),
			std::move(edge_flipped),
			std::move(par),
			std::move(subtree_end),
			std::move(types),
			std::move(orig_id),
			std::move(ch),
			std::move(node_verts),
			std::move(node_nvs),
			std::move(vert_par_nv),
			std::move(node_edges),
			std::move(node_nes),
			std::move(node_adj),
		};
		if constexpr (with_planarity) {
			return planar_spqr_tree{std::move(res), std::move(node_planar), {std::move(ne_rot_adj)}};
		} else {
			return res;
		}
	}
}

inline std::optional<planar_embedding> planar_embed(const planar_spqr_tree& tree) {
	using node_type = planar_spqr_tree::node_type;

	if (!std::ranges::all_of(tree.node_planar, std::identity{})) {
		return std::nullopt;
	}

	int NE = int(tree.edge_index.size());
	std::vector<int> rot_adj(4 * NE, -1);
	auto link = [&] [[gnu::always_inline]] (int a, int b) -> void {
		assert(a != -1 && b != -1);
		assert(rot_adj[a] == -1 && rot_adj[b] == -1);
		assert((a & 1) != (b & 1));
		rot_adj[a] = b;
		rot_adj[b] = a;
	};

	std::vector<std::array<std::array<int, 2>, 2>> outer_e(tree.size(), {{{-1, -1}, {-1, -1}}});
	for (int i = tree.size() - 1; i >= 0; i--) {
		auto type = tree.types[i];
		if (type == node_type::F) {
			for (int j : tree.ch[i]) {
				assert(tree.types[j] == node_type::V);
				auto [a, b] = outer_e[j][0];
				if (a != -1) {
					link(a, b);
				}
			}
		} else if (type == node_type::V) {
			std::array<int, 2> qes{-1, -1};
			for (int j : tree.ch[i]) {
				assert(tree.types[j] == node_type::Q);
				if (qes[0] == -1) {
					qes = outer_e[j][0];
				} else {
					link(qes[1], outer_e[j][0][0]);
					qes[1] = outer_e[j][0][1];
				}
			}
			outer_e[i][0] = qes;
		} else if (type == node_type::Q) {
			int e = tree.orig_id[i];
			bool flip = tree.edge_flipped[e];
			std::array<std::array<int, 2>, 2> qes = {{{4 * e + 2 * flip + 0, 4 * e + 2 * flip + 1}, {4 * e + 2 * !flip + 0, 4 * e + 2 * !flip + 1}}};
			if (tree.ch[i].empty()) {
				// Just return ourselves
				outer_e[i] = qes;
			} else {
				int j = tree.ch[i][0];
				if (tree.types[j] == node_type::O) {
					link(qes[0][1], qes[1][0]);
					outer_e[i][0] = {qes[0][0], qes[1][1]};
				} else {
					if (tree.types[j] != node_type::I) {
						link(qes[0][1], outer_e[j][0][0]);
						qes[0][1] = outer_e[j][0][1];
						link(qes[1][0], outer_e[j][1][1]);
						qes[1][0] = outer_e[j][1][0];
					}
					{
						int k = tree.ch[i][1];
						if (outer_e[k][0][0] != -1) {
							link(qes[1][1], outer_e[k][0][0]);
							link(qes[1][0], outer_e[k][0][1]);
						} else {
							link(qes[1][1], qes[1][0]);
						}
					}
					outer_e[i][0] = qes[0];
				}
			}
		} else if (type == node_type::O || type == node_type::I) {
			// Do nothing, the Q node handles it
		} else if (type == node_type::P || type == node_type::S || type == node_type::R) {
			// Just merge things according to the ne_embedding
			for (int ta = 4 * tree.node_nes.bounds[i]; ta < 4 * tree.node_nes.bounds[i+1]; ta++) {
				int tb = tree.ne_embedding.rot_adj[ta];
				if (tb < ta) continue;

				auto tree_qe_to_qe = [&] [[gnu::always_inline]] (int t) -> int {
					return outer_e[tree.node_edges[tree.node_edges[t>>2].twin_ne].node][(t >> 1) & 1][t & 1];
				};
				int qb = tree_qe_to_qe(tb);
				if (ta < 4 * (tree.node_nes.bounds[i] + 1)) {
					outer_e[i][(ta >> 1) & 1][!(ta & 1)] = qb;
				} else {
					int qa = tree_qe_to_qe(ta);
					if ((ta & 3) == 2 && (tb & 3) == 1) {
						// We're the transition between left and right of a vertex, splice it in here.
						int v = tree.node_verts[tree.node_edges[ta >> 2].nvs[1]].vert;
						if (outer_e[v][0][0] != -1) {
							link(qa, outer_e[v][0][1]);
							link(qb, outer_e[v][0][0]);
						} else {
							link(qa, qb);
						}
					} else {
						link(qa, qb);
					}
				}
			}
		} else assert(false);
	}

	return planar_embedding{std::move(rot_adj)};
}

inline std::optional<planar_embedding> planar_embed(
	int NV,
	const std::vector<std::array<int, 2>>& edges,
	std::span<const int> vert_order,
	std::span<const int> edge_order
) {
	int NE = int(edges.size());
	// Each edge has 4 quarter edges by 4 * edge_id + 2 * side + is_cw
	std::vector<int> quarter_edge_matches(4 * NE, -1);
	// Edges in tstack push order, and side-flip toggles by tstack index (prefix-xor gives each edge's flip)
	std::vector<int> postorder_edges; postorder_edges.reserve(NE);
	std::vector<bool> postorder_flip(NE + 1, false);
	{
		using top_t = spqr_walk_top_t;
		using qe_t = spqr_walk_quarter_edges_t;
		std::vector<int> edge_top_depths(NE, -1);
		qe_t qe{quarter_edge_matches, edge_top_depths};

		spqr_walk_hooks_t hooks{
			.enter_vert = [] [[gnu::always_inline]] (int, int) -> void {},
			.start_edge = [] [[gnu::always_inline]] (int, int, bool) -> void {},
			.vert_entry = [] [[gnu::always_inline]] (auto&, int) -> void {},
			.edge_entry = [&] [[gnu::always_inline]] (auto& t, int e_side, bool is_tree, bool) -> void {
				edge_top_depths[e_side >> 1] = t.top_depth;
				auto& be = t.payload;
				if (is_tree) {
					be[0] = {2 * (e_side ^ 1) + 0, 2 * e_side + 1};
					be[1] = {2 * (e_side ^ 1) + 1, 2 * e_side + 0};
				} else {
					be[0] = {2 * e_side + 0, 2 * e_side + 1};
					t.tops[0] = {{{2 * (e_side ^ 1) + 1, t.top_depth}, {2 * (e_side ^ 1) + 0, t.top_depth}}};
				}
				postorder_edges.push_back(e_side >> 1);
			},
			.merge = [&] [[gnu::always_inline]] (auto& a, const auto& b) -> void {
				qe.merge(a.payload, b.payload);
			},
			.link_tops = [&] [[gnu::always_inline]] (top_t outer, top_t inner) -> void {
				qe.link(outer.end, inner.end);
			},
			.flip = [&] [[gnu::always_inline]] (auto& t, int end_idx) -> void {
				postorder_flip[t.first_idx].flip();
				postorder_flip[end_idx].flip();
				qe.flip(t.payload);
			},
			.join_backedges = [&] [[gnu::always_inline]] (auto& t, int z) -> void {
				qe.join_backedges(t.payload[z], t.tops[z]);
			},
			.join_bottoms = [&] [[gnu::always_inline]] (auto& t, int lowval) -> void {
				qe.join_bottoms(t.payload, t.tops, lowval);
			},
			.pop_top = [&] [[gnu::always_inline]] (auto& t, int z) -> top_t {
				return qe.pop_top(t.payload[z], t.tops[z][1]);
			},
			.absorb_vert = [] [[gnu::always_inline]] (int, int) -> void {},
			.begin_node = [] [[gnu::always_inline]] (auto&, spqr_tree::node_type, bool) -> int { return -1; },
			.finish_node = [] [[gnu::always_inline]] (auto&, int, bool) -> void {},
			.fold_child = [] [[gnu::always_inline]] (auto&, int) -> void {},
			.finish_block = [] [[gnu::always_inline]] (auto&, int, int, bool) -> void {},
			.finish_root = [&] [[gnu::always_inline]] (auto& t) -> void {
				auto [a, b] = t.payload[0];
				if (a != -1) qe.link(a, b);
			},
		};
		if (!spqr_walk<qe_t::bot_ends_t, true, true>(NV, edges, vert_order, edge_order, hooks)) return std::nullopt;
	}
	{
		std::vector<bool> edge_flip(NE);
		{
			bool planarity_flip = false;
			for (int e = 0; e < NE; e++) {
				planarity_flip ^= postorder_flip[e];
				edge_flip[postorder_edges[e]] = planarity_flip;
			}
			planarity_flip ^= postorder_flip[NE];
			assert(!planarity_flip);
		}
		for (int e = 0; e < NE; e++) {
			if (edge_flip[e]) {
				std::swap(quarter_edge_matches[4*e + 0], quarter_edge_matches[4*e + 1]);
				std::swap(quarter_edge_matches[4*e + 2], quarter_edge_matches[4*e + 3]);
			}
			for (int z = 0; z < 4; z++) {
				quarter_edge_matches[4*e + z] = (quarter_edge_matches[4*e+z] >> 1 << 1) | !(z & 1);
			}
		}
	}
	return planar_embedding{std::move(quarter_edge_matches)};
}

inline bool can_planar_embed(
	int NV,
	const std::vector<std::array<int, 2>>& edges,
	std::span<const int> vert_order,
	std::span<const int> edge_order
) {
	int NE = int(edges.size());
	using top_t = spqr_walk_top_t;
	// prev_edge[e] is the next return outwards of e on its side
	std::vector<top_t> prev_edge(NE);

	spqr_walk_hooks_t hooks{
		.enter_vert = [] [[gnu::always_inline]] (int, int) -> void {},
		.start_edge = [] [[gnu::always_inline]] (int, int, bool) -> void {},
		.vert_entry = [] [[gnu::always_inline]] (auto&, int) -> void {},
		.edge_entry = [] [[gnu::always_inline]] (auto& t, int e_side, bool is_tree, bool) -> void {
			if (is_tree) return;
			int e = e_side >> 1;
			t.tops[0] = {{{e, t.top_depth}, {e, t.top_depth}}};
		},
		.merge = [] [[gnu::always_inline]] (auto&, const auto&) -> void {},
		.link_tops = [&] [[gnu::always_inline]] (top_t outer, top_t inner) -> void {
			prev_edge[inner.end] = outer;
		},
		.flip = [] [[gnu::always_inline]] (auto&, int) -> void {},
		.join_backedges = [] [[gnu::always_inline]] (auto&, int) -> void {},
		.join_bottoms = [] [[gnu::always_inline]] ([[maybe_unused]] auto& t, [[maybe_unused]] int lowval) -> void {
			// The lowval-only returns of side 1 carry no information:
			// side 0's outermost return is already at lowval, and they would be the first ones pruned at depth lowval.
			assert(t.tops[1][0].end == -1 || (t.tops[1][0].depth == lowval && t.tops[1][1].depth == lowval && t.tops[0][0].depth == lowval));
		},
		.pop_top = [&] [[gnu::always_inline]] (auto& t, int z) -> top_t {
			return prev_edge[t.tops[z][1].end];
		},
		.absorb_vert = [] [[gnu::always_inline]] (int, int) -> void {},
		.begin_node = [] [[gnu::always_inline]] (auto&, spqr_tree::node_type, bool) -> int { return -1; },
		.finish_node = [] [[gnu::always_inline]] (auto&, int, bool) -> void {},
		.fold_child = [] [[gnu::always_inline]] (auto&, int) -> void {},
		.finish_block = [] [[gnu::always_inline]] (auto&, int, int, bool) -> void {},
		.finish_root = [] [[gnu::always_inline]] (auto&) -> void {},
	};
	return spqr_walk<std::monostate, true, true>(NV, edges, vert_order, edge_order, hooks);
}

} // namespace wala
