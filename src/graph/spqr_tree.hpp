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

namespace wala {

template <typename T> struct csr {
	std::vector<int> bounds;
	std::vector<T> dat;
	std::span<T> operator [](int i) { return std::span<T>(dat.begin() + bounds[i], dat.begin() + bounds[i + 1]); }
	std::span<const T> operator [](int i) const { return std::span<const T>(dat.begin() + bounds[i], dat.begin() + bounds[i + 1]); }
	std::ranges::iota_view<int, int> indices(int i) const { return std::views::iota(bounds[i], bounds[i+1]); }
	int size() const { return int(bounds.size()) - 1; }
};

template <typename T> struct csr_builder {
	std::vector<int> bounds;
	std::vector<T> dat;
	csr_builder() = default;
	explicit csr_builder(int N) : bounds(N+1) {}

	void count(int k) { bounds[k+1]++; }
	void allocate() {
		int l = 0;
		for (int i = 1; i < int(bounds.size()); i++) {
			bounds[i] = std::exchange(l, l + bounds[i]);
		}
		dat.resize(l);
	}
	[[nodiscard]] T& push(int k) { return dat[bounds[k+1]++]; }
	[[nodiscard]] csr<T> finalize() && { return { std::move(bounds), std::move(dat) }; }
};

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

	std::vector<int> par;
	std::vector<int> subtree_end;
	std::vector<node_type> types;
	std::vector<int> orig_id;

	csr<int> ch;
	struct node_vert_t {
		int node;
		int vert;
	};
	csr<node_vert_t> node_verts;
	// The nv index of a vertex within its parent node
	std::vector<int> vert_par_nv;
	// TODO: Should we store a vert_nodes CSR?

	struct node_edge_t {
		int node;
		int twin_ne;
		// TODO: Should we store the twin node, the twin node type, and/or twin node type == Q?

		std::array<int, 2> nvs;
	};
	csr<node_edge_t> node_edges;

	struct node_adj_t {
		int ne;
		int dest_nv;
	};
	csr<node_adj_t> node_adj;

	// Planarity data, only filled in if built with with_planarity = true (otherwise both are empty).
	std::vector<bool> node_planar;
	// Planarity adjacencies: ne_rot_adj is an involution of facing quarter-edges, indexed according to:
	// ne_rot_adj[4 * node_edge + 2 * side + dir]
	std::vector<int> ne_rot_adj;

	int size() const { return int(par.size()); }

	template <bool with_planarity = true>
	static spqr_tree build(int NV, const std::vector<std::array<int, 2>>& edges, bool ternarize = false) {
		// TODO: Figure out the best way to specify roots; maybe accept a permutation of "root priority"?

		// std::min is by reference, which breaks some optimizations
		auto min = [](auto a, auto b) { return a < b ? a : b; };
		auto setmin = [](auto& a, auto b) { if (b < a) a = b; };

		int NE = int(edges.size());

		std::vector<int> roots; roots.reserve(NV);
		struct outedge_t { int src, dest; int e; int key; };
		csr<outedge_t> outedges;

		// Phase 1: build a sorted skeleton
		{
			std::vector<int> depth(NV, -1);
			// 1a: build a normal adjacency list for the initial lowval dfs
			struct edge_t { int dest; int e; };
			csr_builder<edge_t> adj_builder(NV);
			for (auto [u, v] : edges) {
				adj_builder.count(u);
				if (u != v) adj_builder.count(v);
			}
			adj_builder.allocate();
			for (int e = 0; e < NE; e++) {
				auto [u, v] = edges[e];
				adj_builder.push(u) = {v, e};
				if (u != v) adj_builder.push(v) = {u, e};
			}
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
			auto push_vert = [&](int cur, int prv_e) -> void {
				int d = int(stk.size());
				depth[cur] = d;
				stk.push_back({cur, prv_e, {d, d}, adj.bounds[cur], adj.bounds[cur+1]});
			};
			auto finish_edge = [&](bool is_tree, std::array<int, 2> n_lowvals) -> void {
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
					// We don't really need to distinguish backedges vs type-1 children, but do it just for fun?
					int lowval = n_lowvals[0];
					if (lowval >= d) lowval = ~(lowval - d);
					int kind = 2 * (n_lowvals[1] < d) + !is_tree;
					all_outedges.push_back({cur, nxt, e, 3 * (lowval + 2) + kind});
				}

				// Keep the 2 distinct mins
				if (n_lowvals[0] < lowvals[0]) lowvals = {n_lowvals[0], min(n_lowvals[1], lowvals[0])};
				else lowvals[1] = min(lowvals[1], n_lowvals[0] == lowvals[0] ? n_lowvals[1] : n_lowvals[0]);
			};
			auto start_edge = [&]() -> void {
				int d = int(stk.size()) - 1;
				auto& s = stk.back();
				assert(s.ch_idx < s.ch_end);
				auto [nxt, e] = adj.dat[s.ch_idx];

				if (e == s.prv_e || depth[nxt] > d) {
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
			auto pop_vert = [&]() -> std::array<int, 2> {
				auto lowvals = stk.back().lowvals;
				stk.pop_back();
				return lowvals;
			};
			for (int rt = 0; rt < NV; rt++) {
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
			}

			csr_builder<outedge_t> by_key_builder(3*NV+6);
			for (auto edge : all_outedges) by_key_builder.count(edge.key);
			by_key_builder.allocate();
			for (auto edge : all_outedges) by_key_builder.push(edge.key) = edge;
			csr<outedge_t> by_key = std::move(by_key_builder).finalize();

			csr_builder<outedge_t> by_src_builder(NV);
			// Hack to reuse memory
			by_src_builder.dat = std::move(all_outedges);
			for (auto edge : by_key.dat) by_src_builder.count(edge.src);
			by_src_builder.allocate();
			for (auto edge : by_key.dat) by_src_builder.push(edge.src) = edge;
			outedges = std::move(by_src_builder).finalize();
		}

		// Phase 2: do the big ear-decomposition-like walk

		// We're going to build a tree of all SPQR *nodes* + all original *vertices* (collectively *items*).
		// Vertices will hang off the first SPQR node containing them, and blocks will be rooted at a topmost Q node for the top edge.

		// As we build, we will represent the children of our nodes/vertices as linked lists.
		constexpr int ROOT_ITEM = 0;
		auto vert_item = [&](int v) -> int { return 1 + v; };
		auto edge_item = [&](int e) -> int { return 1 + NV + e; };

		// Helpers for working with std::array<T, 2> - these compile to cmov's better than direct index access.

		// return arr[dir] == a, arr[!dir] == b
		auto set_sides = []<typename T>(bool dir, T a, T b) -> std::array<T, 2> {
			return dir ? std::array<T, 2>{b, a} : std::array<T, 2>{a, b};
		};
		auto get_side = []<typename T>(std::array<T, 2> a, bool dir) -> T {
			return dir ? a[1] : a[0];
		};

		// Planarity flip bits ride along the child linked lists; without planarity they collapse to empty structs.
		struct no_planarity_flip_t {};
		struct item_list {
			std::array<int, 2> v{-1, -1};
			[[no_unique_address]] std::conditional_t<with_planarity, std::array<bool, 2>, no_planarity_flip_t> planarity_flip{};

			[[nodiscard]] bool empty() const { return v[0] == -1; }
		};
		struct ch_nxt_t {
			int nxt = -1;
			[[no_unique_address]] std::conditional_t<with_planarity, bool, no_planarity_flip_t> planarity_flip{};
		};
		std::vector<ch_nxt_t> ch_nxt; ch_nxt.reserve(1 + NV + NE + NE); ch_nxt.assign(1 + NV + NE, {});
		auto concat = [&](item_list a, item_list b) -> item_list {
			if (b.empty()) return a;
			if (a.empty()) return b;
			if constexpr (with_planarity) {
				ch_nxt[a.v[1]] = {b.v[0], a.planarity_flip[1] != b.planarity_flip[0]};
				return {{a.v[0], b.v[1]}, {a.planarity_flip[0], b.planarity_flip[1]}};
			} else {
				ch_nxt[a.v[1]] = {b.v[0]};
				return {{a.v[0], b.v[1]}};
			}
		};
		auto unit_list = [&](int item) -> item_list {
			return {{item, item}};
		};

		std::vector<std::array<int, 2>> item_vs; item_vs.reserve(1 + NV + 2 * NE); item_vs.resize(1 + NV + NE, {-1, -1});
		std::vector<item_list> item_ch; item_ch.reserve(1 + NV + 2 * NE); item_ch.resize(1 + NV + NE, item_list{});
		std::vector<node_type> item_types; item_types.reserve(1 + NV + 2 * NE);
		item_types.resize(1, node_type::F);
		item_types.resize(1 + NV, node_type::V);
		item_types.resize(1 + NV + NE, node_type::Q);

		// Quarter edges for planar embedding building.
		// Each vedge has 4 entries by 4 * vedge_id + 2 * source_vert + is_cw (is_cw is arbitrary)
		// vedges are identified with what item they cap, numbered by (item - 1 - NV)
		std::vector<int> quarter_edge_matches(with_planarity ? 8 * NE + 4 : 0, -1);
		struct nonplanarity_certficate_t {};
		std::vector<std::expected<std::array<int, 4>, nonplanarity_certficate_t>> node_planarity;
		if constexpr (with_planarity) node_planarity.reserve(NE);

		int tot_blocks = 0;
		int tot_self_loops = 0;

		{
			auto alloc_item = [&](node_type type) -> int {
				int item = int(item_vs.size());
				item_vs.push_back({});
				item_ch.push_back({});
				item_types.push_back(type);
				ch_nxt.push_back({});
				if constexpr (with_planarity) node_planarity.emplace_back();
				return item;
			};

			// Declare these here: most of our code will be in terms of v_start / top_depth, so we'll want to read these out
			std::vector<int> stack_verts(NV);
			std::vector<int8_t> stack_dir(NV); // really bool, but I don't want vector<bool>

			auto make_vs = [&](int v_start, int top_depth) -> std::array<int, 2> {
				return set_sides(stack_dir[top_depth], stack_verts[top_depth], v_start);
			};

			int nxt_edge_idx = 0; // Counts backedges only
			std::vector<int> first_occurrence(NV); // First backedge to this depth

			std::vector<int> edge_top_depths(with_planarity ? 2 * NE : 0, -1);

			struct tstack_planarity_side_t {
				// For each side, store pointers to the "linked lists" of the edges inside.
				// v[0] is the outer / longer edges and v[1] is the inner / shorter edges, matching the outside-in sort order.

				// bot_ends are the outer/innermost exposed pieces of the walk down the ear in the tree (they're connected to the bottommost/topmost vertices of the tree path)
				std::array<int, 2> bot_ends{-1, -1};
				// top_ends are the outer/innermost exposed backedges
				std::array<int, 2> top_ends{-1, -1};
				// depths should be increasing going inwards
				std::array<int, 2> top_depths{-1, -1};
			};
			struct tstack_planarity_t {
				// The convention is that sides[0].top_depths[0] == top_depth, i.e. at least one minimal return lives on side 0
				std::array<tstack_planarity_side_t, 2> sides;
			};
			struct tstack_nonplanarity_t {
				// TODO: What's the nonplanarity certificate look like?
			};
			struct tstack_no_planarity_t {};
			// With planarity disabled, the planarity state collapses to an empty struct and every planarity step is skipped.
			using tstack_edge_planarity_t = std::conditional_t<with_planarity, tstack_planarity_t, tstack_no_planarity_t>;
			using tstack_maybe_planarity_t = std::conditional_t<with_planarity, std::expected<tstack_planarity_t, tstack_nonplanarity_t>, tstack_no_planarity_t>;
			auto merge_planarity = [&](tstack_maybe_planarity_t& a, const tstack_maybe_planarity_t& b) -> void {
				if constexpr (!with_planarity) return;
				else {
				if (!a) return;
				if (!b) { a = b; return; }
				tstack_planarity_t res;
				for (int z = 0; z < 2; z++) {
					auto& as = a->sides[z];
					const auto& bs = b->sides[z];
					// If there's no bottom edges, then we must be an isolated vertex, so we can end early.
					if (bs.bot_ends[0] == -1) {
						// Do nothing
					} else if (as.bot_ends[0] == -1) {
						as = bs;
					} else {
						quarter_edge_matches[as.bot_ends[1]] = bs.bot_ends[0];
						quarter_edge_matches[bs.bot_ends[0]] = as.bot_ends[1];
						as.bot_ends[1] = bs.bot_ends[1];

						if (bs.top_ends[0] == -1) {
							// Do nothing
						} else if (as.top_ends[0] == -1) {
							as.top_ends = bs.top_ends;
							as.top_depths = bs.top_depths;
						} else if (as.top_depths[1] > bs.top_depths[0]) {
							// TODO: Certificate
							a = std::unexpected(tstack_nonplanarity_t{});
							return;
						} else {
							quarter_edge_matches[as.top_ends[1]] = bs.top_ends[0];
							quarter_edge_matches[bs.top_ends[0]] = as.top_ends[1];
							as.top_ends[1] = bs.top_ends[1];
							as.top_depths[1] = bs.top_depths[1];
						}
					}
				}
				}
			};
			auto make_edge_planarity = [&](int item, int top_depth, bool is_tree) -> tstack_edge_planarity_t {
				if constexpr (!with_planarity) return {};
				else {
				assert(item >= 1 + NV);
				int ve = item - (1 + NV);
				bool top_dir = stack_dir[top_depth];
				edge_top_depths[ve] = top_depth;
				tstack_planarity_t p;
				if (is_tree) {
					p.sides[0].bot_ends = {4 * ve + 2 * !top_dir + 0, 4 * ve + 2 * top_dir + 1};
					p.sides[1].bot_ends = {4 * ve + 2 * !top_dir + 1, 4 * ve + 2 * top_dir + 0};
				} else {
					p.sides[0].bot_ends = {4 * ve + 2 * !top_dir + 0, 4 * ve + 2 * !top_dir + 1};
					p.sides[0].top_ends = {4 * ve + 2 * top_dir + 1, 4 * ve + 2 * top_dir + 0};
					p.sides[0].top_depths = {top_depth, top_depth};
				}
				return p;
				}
			};
			struct tstack_t {
				int v_start = -1;
				int top_depth = -1;
				int first_idx = -1;
				std::array<item_list, 2> spans;
				[[no_unique_address]] tstack_maybe_planarity_t planarity;
			};
			int tstack_size = 0;
			std::vector<tstack_t> tstack(NV + NE);
			auto cur_tstack = [&]() -> tstack_t& { return tstack[tstack_size-1]; };
			auto nxt_tstack = [&]() -> tstack_t& { return tstack[tstack_size-2]; };

			auto push_tstack = [&](int v_start, int top_depth, int item, tstack_edge_planarity_t planarity) -> void {
				tstack[tstack_size++] = { v_start, top_depth, nxt_edge_idx, set_sides(stack_dir[top_depth], unit_list(item), {}), planarity };
			};
			auto flip_tstack_planarity = [&](tstack_t& a) -> void {
				if constexpr (with_planarity) {
					a.spans[0].planarity_flip[0] ^= 1;
					a.spans[0].planarity_flip[1] ^= 1;
					a.spans[1].planarity_flip[0] ^= 1;
					a.spans[1].planarity_flip[1] ^= 1;
					if (a.planarity) {
						std::swap(a.planarity->sides[0], a.planarity->sides[1]);
					}
				}
			};
			auto merge_tstack_tops = [&]() -> void {
				tstack_t& a = nxt_tstack();
				const tstack_t& b = cur_tstack();
				setmin(a.top_depth, b.top_depth);
				a.spans[0] = concat(b.spans[0], a.spans[0]);
				a.spans[1] = concat(a.spans[1], b.spans[1]);
				merge_planarity(a.planarity, b.planarity);
				tstack_size--;
			};

			auto maybe_unwrap_nxt = [&](node_type type, bool is_tree) -> int {
				tstack_t& t = nxt_tstack();

				if (type == node_type::R) return alloc_item(type);

				assert(type == node_type::P || type == node_type::S);

				// If we want to ternarize, never reuse.
				if (ternarize) return alloc_item(type);

				bool top_dir = stack_dir[t.top_depth];
				assert(get_side(t.spans, !top_dir).empty());
				int item = get_side(t.spans, top_dir).v[0];
				assert(item == get_side(t.spans, top_dir).v[1]);
				if (item_types[item] == type) {
					t.spans = set_sides(top_dir, item_ch[item], {});
					if constexpr (with_planarity) {
						// Unwrap the planarity data
						// We don't really need to maintain this at all because S/P nodes are known to be trivially planar
						// The current state is just make_edge_planarity(wrapped), which means that it has the right shape, just needs to be relabelled.
						assert(node_planarity[item - (1 + NV + NE)]);
						const auto& matches = *node_planarity[item - (1 + NV + NE)];
						assert(t.planarity);
						auto& p = *t.planarity;
						if (is_tree) {
							p.sides[0].bot_ends[0] = matches[2 * !top_dir + 1];
							p.sides[0].bot_ends[1] = matches[2 * top_dir + 0];
							p.sides[1].bot_ends[0] = matches[2 * !top_dir + 0];
							p.sides[1].bot_ends[1] = matches[2 * top_dir + 1];
						} else {
							p.sides[0].bot_ends[0] = matches[2 * !top_dir + 1];
							p.sides[0].bot_ends[1] = matches[2 * !top_dir + 0];
							p.sides[0].top_ends[0] = matches[2 * top_dir + 0];
							p.sides[0].top_ends[1] = matches[2 * top_dir + 1];
						}
					}
					return item;
				} else {
					return alloc_item(type);
				}
			};

			auto finish_tstack_top = [&](int item, bool is_tree) -> void {
				tstack_t& t = cur_tstack();
				bool top_dir = stack_dir[t.top_depth];
				assert(get_side(t.spans, !top_dir).empty());

				if constexpr (!with_planarity) {
				} else if (t.planarity) {
					const auto& p = *t.planarity;
					std::array<int, 4> matches{};
					if (is_tree) {
						matches[2 * !top_dir + 1] = p.sides[0].bot_ends[0];
						matches[2 * top_dir + 0] = p.sides[0].bot_ends[1];
						matches[2 * !top_dir + 0] = p.sides[1].bot_ends[0];
						matches[2 * top_dir + 1] = p.sides[1].bot_ends[1];
					} else {
						matches[2 * !top_dir + 1] = p.sides[0].bot_ends[0];
						matches[2 * !top_dir + 0] = p.sides[0].bot_ends[1];
						matches[2 * top_dir + 0] = p.sides[0].top_ends[0];
						matches[2 * top_dir + 1] = p.sides[0].top_ends[1];
					}
					node_planarity[item - (1 + NV + NE)] = matches;
				} else {
					assert(item_types[item] == node_type::R);
					node_planarity[item - (1 + NV + NE)] = std::unexpected(nonplanarity_certficate_t{});
				}
				item_vs[item] = make_vs(t.v_start, t.top_depth);
				item_ch[item] = get_side(t.spans, top_dir);

				t.spans = set_sides(top_dir, unit_list(item), {});
				t.planarity = make_edge_planarity(item, t.top_depth, is_tree);
			};

			struct dfs_stack_t {
				bool pushed_vert;
				int ch_idx;
				int ch_end;
				int orig_tstack;
			};
			std::vector<dfs_stack_t> stk; stk.reserve(NV);
			for (auto rt : roots) {
				auto push_vert = [&](int cur) -> void {
					stk.push_back({false, outedges.bounds[cur], outedges.bounds[cur+1], -1});
					int cur_depth = int(stk.size()) - 1;
					stack_verts[cur_depth] = cur;
				};
				struct key_t { int lowval; bool is_tree; bool is_type_1; };
				auto decode_key = [&](int cur_depth, int key) -> key_t {
					int lowval = key / 3 - 2; if (lowval < 0) lowval = cur_depth + ~lowval;
					int kind = key % 3;
					bool is_tree = kind != 1;
					bool is_type_1 = kind <= 1;
					return {lowval, is_tree, is_type_1 };
				};
				// return true means jump to start_edge, return false means jump to finish_edge
				auto start_edge = [&]() -> std::optional<int> {
					int cur_depth = int(stk.size()) - 1;
					auto& s = stk.back();
					int cur = stack_verts[cur_depth];
					assert(s.ch_idx < s.ch_end);
					auto [_, nxt, e, key] = outedges.dat[s.ch_idx];
					auto [lowval, is_tree, is_type_1] = decode_key(cur_depth, key);

					// edge_dir convention: false is forwards, true is backwards.
					// That means that cur is on the edge_dir side and nxt is on the !edge_dir side.
					stack_dir[cur_depth] = (lowval >= cur_depth ? false : !stack_dir[lowval]);

					if (!s.pushed_vert && lowval < cur_depth && is_type_1) {
						// Do this with the correct stack_dir set
						push_tstack(cur, cur_depth, vert_item(cur), {});
						s.pushed_vert = true;
					}

					s.orig_tstack = tstack_size;
					if (is_tree) {
						first_occurrence[cur_depth] = NE;
						return nxt;
					} else {
						return std::nullopt;
					}
				};
				auto finish_edge = [&]() -> void {
					int cur_depth = int(stk.size()) - 1;
					auto& s = stk.back();
					int cur = stack_verts[cur_depth];
					assert(s.ch_idx < s.ch_end);

					auto [_, nxt, e, key] = outedges.dat[s.ch_idx];
					s.ch_idx++;

					auto [lowval, is_tree, is_type_1] = decode_key(cur_depth, key);

					const int orig_tstack = s.orig_tstack;
					const bool edge_dir = stack_dir[cur_depth];

					if (lowval >= cur_depth) {
						// There's no planarity handling for this because it's just a Q node. I/O nodes also don't need any tracking.
						item_vs[edge_item(e)] = {cur, -1};
						tot_blocks++;
						if (is_tree) {
							// Bridges and components
							if (lowval == cur_depth + 1) {
								// tstack[tstack_size-1] is currently just smuggling out the child vertex, prepend the bridge component
								// This is just a shortcut for allocating a full I-type tstack
								int item = alloc_item(node_type::I);
								item_vs[item] = make_vs(nxt, cur_depth);
								item_ch[edge_item(e)] = concat(unit_list(item), tstack[--tstack_size].spans[1]);
							} else {
								// tstack[tstack_size-2] is the vertex and tstack[tstack_size-1] is the backedge
								auto backedge = tstack[--tstack_size].spans[0];
								item_ch[edge_item(e)] = concat(backedge, tstack[--tstack_size].spans[1]);
							}
						} else {
							// self loops
							assert(nxt == cur);
							tot_self_loops++;
							int item = alloc_item(node_type::O);
							// Make sure the nxt is -1 as well
							item_vs[item] = {cur, -1};
							item_ch[edge_item(e)] = unit_list(item);
						}
						item_ch[vert_item(cur)] = concat(item_ch[vert_item(cur)], unit_list(edge_item(e)));
						return;
					}

					item_vs[edge_item(e)] = make_vs(nxt, cur_depth);

					assert(lowval < cur_depth);

					bool is_single = true;

					// make_q_node
					if (is_tree) {
						// The span lives on side edge_dir
						push_tstack(nxt, cur_depth, edge_item(e), make_edge_planarity(edge_item(e), cur_depth, true));
						while (tstack_size >= 2 && nxt_tstack().top_depth >= cur_depth) {
							node_type type;
							if (nxt_tstack().top_depth > cur_depth) {
								// Just backfill this for maybe_unwrap
								stack_dir[nxt_tstack().top_depth] = edge_dir;

								// The tstack currently contains a tree-edge followed by a vertex; merge the vertex first
								merge_tstack_tops();

								type = node_type::S;
							} else if (nxt_tstack().v_start == cur_tstack().v_start) {
								// This will be a P node
								type = node_type::P;
							} else {
								type = node_type::R;
							}
							int item = maybe_unwrap_nxt(type, type == node_type::S);
							merge_tstack_tops();
							if constexpr (!with_planarity) {
							} else if (cur_tstack().planarity) {
								// Merge all backedges into the component
								for (auto& side : cur_tstack().planarity->sides) {
									assert(side.bot_ends[1] != -1);
									if (side.top_ends[1] == -1) continue;
									assert(side.top_depths[0] == cur_depth);
									assert(side.top_depths[1] == cur_depth);
									quarter_edge_matches[side.bot_ends[1]] = side.top_ends[1];
									quarter_edge_matches[side.top_ends[1]] = side.bot_ends[1];
									side.bot_ends[1] = side.top_ends[0];
									side.top_depths = {-1, -1};
									side.top_ends = {-1, -1};
								}
							}
							finish_tstack_top(item, true);
						}

						if (cur_tstack().first_idx > first_occurrence[cur_depth]) {
							while (cur_tstack().first_idx > first_occurrence[cur_depth]) {
								if (nxt_tstack().first_idx > first_occurrence[cur_depth]) {
									// We will put cur_depth on side 1 until the bottom
									if (nxt_tstack().top_depth == cur_depth) {
										flip_tstack_planarity(nxt_tstack());
									}
								} else if (!is_single) {
									assert(cur_tstack().top_depth < cur_depth);
									if constexpr (!with_planarity) {
									} else if (nxt_tstack().planarity) {
										if (nxt_tstack().planarity->sides[0].top_depths[1] == cur_depth) {
											// We need to flip cur_tstack and nxt_tstack relative to each other.
											// Flip the one with worse top_depth.
											flip_tstack_planarity(cur_tstack().top_depth < nxt_tstack().top_depth ? nxt_tstack() : cur_tstack());
										} else {
											assert(nxt_tstack().planarity->sides[1].top_depths[1] == cur_depth);
										}
									}
								}
								merge_tstack_tops();
								is_single = false;
							}
							if constexpr (!with_planarity) {
							} else if (cur_tstack().planarity) {
								// Prune off finished cur-side things
								for (auto& side : cur_tstack().planarity->sides) {
									assert(side.bot_ends[1] != -1);
									while (side.top_depths[1] == cur_depth) {
										{
											// Link these to bot_ends[1]
											quarter_edge_matches[side.bot_ends[1]] = side.top_ends[1];
											quarter_edge_matches[side.top_ends[1]] = side.bot_ends[1];
											side.bot_ends[1] = side.top_ends[1] ^ 1;
										}
										side.top_ends[1] = std::exchange(quarter_edge_matches[side.bot_ends[1]], -1);
										if (side.top_ends[1] != -1) {
											quarter_edge_matches[side.top_ends[1]] = -1;
											side.top_depths[1] = edge_top_depths[side.top_ends[1] >> 2];
										} else {
											side.top_depths = {-1, -1};
											side.top_ends = {-1, -1};
										}
									}
								}
							}
						}

						if (s.pushed_vert) {
							// NB: tstack[orig_size] is the vertex and tstack[orig_size+1] is the backedge; maybe we should reverse them?
							assert(tstack_size >= orig_tstack + 3);

							if (!is_type_1) {
								// The lowval side should be side 1, everything else goes on side 0.
								// The exception is tstack[orig_tstack + 2], which could be == lowval on one/both sides,
								// but is guaranteed to have *something* > lowval by non-type-1-ness
								{
									auto& t = tstack[orig_tstack + 2];
									if constexpr (!with_planarity) {
									} else if (t.planarity) {
										assert(t.planarity->sides[0].top_depths[0] == t.top_depth);
										if (t.planarity->sides[0].top_depths[1] == lowval) {
											flip_tstack_planarity(t);
										}
										assert(t.planarity->sides[0].top_depths[1] != -1);
										assert(t.planarity->sides[0].top_depths[1] > lowval);
									}
								}
								for (int i = orig_tstack + 3; i < tstack_size; i++) {
									if (tstack[i].top_depth == lowval) {
										flip_tstack_planarity(tstack[i]);
									}
								}
								while (tstack_size > orig_tstack + 3) {
									merge_tstack_tops();
									is_single = false;
								}
								assert(!is_single);
							}

							assert(tstack_size == orig_tstack + 3);
							int item;
							if (is_type_1) {
								item = maybe_unwrap_nxt(is_single ? node_type::S : node_type::R, false);
							} else {
								// Just for the type checker
								item = -1;
							}
							// Merge with the backedge
							merge_tstack_tops();
							// Merge with the vertex
							merge_tstack_tops();

							cur_tstack().v_start = cur;
							assert(cur_tstack().top_depth == lowval);

							// Fold everything to the correct side now that we're leaving the child.
							// The entire subtree should go to the !edge_dir side.
							cur_tstack().spans = set_sides(!edge_dir, concat(cur_tstack().spans[0], cur_tstack().spans[1]), {});

							[&]() -> void {
								if constexpr (!with_planarity) {
								} else if (cur_tstack().planarity) {
									// precondition: side 1 should be the lowval only side
									auto& sides = cur_tstack().planarity->sides;
									auto& s0 = sides[0];
									auto& s1 = sides[1];
									quarter_edge_matches[s0.bot_ends[0]] = s1.bot_ends[0];
									quarter_edge_matches[s1.bot_ends[0]] = s0.bot_ends[0];
									s0.bot_ends[0] = s1.bot_ends[1];
									if (s1.top_ends[0] != -1) {
										if (s1.top_depths[1] != lowval) {
											assert(!is_type_1);
											cur_tstack().planarity = std::unexpected(tstack_nonplanarity_t{});
											return;
										}
										assert(s1.top_depths[0] == lowval);
										quarter_edge_matches[s0.top_ends[0]] = s1.top_ends[0];
										quarter_edge_matches[s1.top_ends[0]] = s0.top_ends[0];
										s0.top_ends[0] = s1.top_ends[1];
										// Already true since the backedge was on side 0
										assert(s0.top_depths[0] == lowval);
									}
									s1 = tstack_planarity_side_t{};
								}
							}();

							if (is_type_1) {
								finish_tstack_top(item, false);
								is_single = true;
							}
						}
					} else {
						assert(is_type_1);
						// The span lives on side !edge_dir
						push_tstack(cur, lowval, edge_item(e), make_edge_planarity(edge_item(e), lowval, false));
						setmin(first_occurrence[lowval], nxt_edge_idx++);
					}

					// NB: We can do this check in lots of ways, maybe there's a cleaner check
					if (is_type_1 && tstack_size >= 2 && nxt_tstack().v_start == cur && nxt_tstack().top_depth == lowval) {
						// This will be a P node
						int item = maybe_unwrap_nxt(node_type::P, false);
						merge_tstack_tops();
						finish_tstack_top(item, false);
					}

					if (!s.pushed_vert) {
						// Throw cur_vert_node onto the tstack so it'll get interleaved correctly
						push_tstack(cur, cur_depth, vert_item(cur), {});
						s.pushed_vert = true;
						assert(!is_type_1);
						if (!is_single) {
							// Just eagerly merge the vertex into the R to avoid a later spurious finish_tstack
							merge_tstack_tops();
						}
					}
				};
				auto pop_vert = [&]() -> void {
					int cur_depth = int(stk.size()) - 1;
					auto& s = stk.back();
					int cur = stack_verts[cur_depth];
					assert(s.ch_idx == s.ch_end);
					if (!s.pushed_vert) {
						// Either our parent is a bridge edge, or we're just a root.
						// We'll just leave it on tstack for future cleanup, it'll just get popped of immediately.
						// edge_dir == !stack_dir[lowval == cur_depth - 1] == true
						stack_dir[cur_depth] = true;
						push_tstack(cur, cur_depth, vert_item(cur), {});
						s.pushed_vert = true;
					}
					stk.pop_back();
				};

				push_vert(rt);
				while (true) {
					if (stk.back().ch_idx == stk.back().ch_end) {
						pop_vert();
						if (stk.empty()) break;
						finish_edge();
					} else if (std::optional<int> nxt = start_edge(); nxt) {
						push_vert(*nxt);
					} else {
						finish_edge();
					}
				}
				item_ch[ROOT_ITEM] = concat(item_ch[ROOT_ITEM], tstack[--tstack_size].spans[1]);
			}
		}

		// Phase 3: relabel the full tree in preorder
		int tot_items = int(item_types.size());
		{
			std::vector<int> vert_index(NV, -1);
			std::vector<int> edge_index(NE, -1);

			std::vector<int> par(tot_items, -1);
			std::vector<int> subtree_end(tot_items, -1);
			std::vector<node_type> types(tot_items, node_type::F);
			std::vector<int> orig_id(tot_items, -1);

			csr<int> ch;
			ch.bounds.resize(tot_items + 1, 0);
			ch.dat.resize(tot_items - 1);

			// Each node is a child, and additionally most non-block node has 2 cap verts; blocks have 1, and O nodes have 1
			int tot_node_verts = NV + (tot_items - 1 - NV) * 2 - tot_blocks - tot_self_loops;
			csr<node_vert_t> node_verts;
			node_verts.bounds.resize(tot_items + 1);
			node_verts.dat.resize(tot_node_verts);
			std::vector<int> vert_par_nv(tot_items, -1);

			int tot_node_edges = (tot_items - 1 - NV - tot_blocks) * 2;
			csr<node_edge_t> node_edges;
			node_edges.bounds.resize(tot_items + 1);
			node_edges.dat.resize(tot_node_edges);

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
			auto push_item = [&](int cur_item) -> void {
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
				} else {
					assert(1 + NV + NE <= cur_item);
					if constexpr (!with_planarity) {
					} else if (cur_type == node_type::O || cur_type == node_type::I) {
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
				if constexpr (with_planarity) node_planar[cur_idx] = planar;

				// HACK: Fill ch and vert_items in with orig items / orig verts for now,
				// because we don't have the final item id's yet.
				int ch_st = ch.bounds[cur_idx];
				int ch_en = ch_st;
				int nv_st = node_verts.bounds[cur_idx];
				int nv_en = nv_st;
				int n_edges = 0;
				if (item_vs[cur_item][0] != -1) {
					node_verts.dat[nv_en++] = {cur_idx, item_vs[cur_item][0]};
				}
				bool planarity_flip = false;
				if constexpr (with_planarity) planarity_flip = item_ch[cur_item].planarity_flip[0];
				for (int nxt_item = item_ch[cur_item].v[0]; nxt_item != -1; nxt_item = ch_nxt[nxt_item].nxt) {
					ch.dat[ch_en++] = nxt_item;
					assert(nxt_item >= 1);
					if (nxt_item < 1 + NV) {
						node_verts.dat[nv_en++] = {cur_idx, nxt_item - 1};
					} else {
						if constexpr (!with_planarity) {
						} else if (cur_type != node_type::R) {
							assert(!planarity_flip);
						} else {
							// Fix the planarity direction right here: reverse quarter_edge_matches upfront;
							// this breaks the involution property, but from here on we'll never read the low bits anyways.
							int ve = nxt_item - (1 + NV);
							if (planarity_flip) {
								std::swap(quarter_edge_matches[4 * ve + 0], quarter_edge_matches[4 * ve + 1]);
								std::swap(quarter_edge_matches[4 * ve + 2], quarter_edge_matches[4 * ve + 3]);
							}
						}
						n_edges++;
					}
					if (nxt_item == item_ch[cur_item].v[1]) {
						assert(ch_nxt[nxt_item].nxt == -1);
					}
					if constexpr (with_planarity) planarity_flip ^= ch_nxt[nxt_item].planarity_flip;
				}
				if constexpr (with_planarity) {
					planarity_flip ^= item_ch[cur_item].planarity_flip[1];
					assert(!planarity_flip);
				}
				if (item_vs[cur_item][1] != -1) {
					node_verts.dat[nv_en++] = {cur_idx, item_vs[cur_item][1]};
				}
				ch.bounds[cur_idx+1] = ch_en;
				node_verts.bounds[cur_idx+1] = nv_en;

				int n_verts = nv_en - nv_st;

				bool is_node = cur_type != node_type::F && cur_type != node_type::V;
				bool has_cap = is_node && !(cur_type == node_type::Q && ch_en - ch_st > 0);

				if (!is_node) n_edges = 0;
				if (has_cap) n_edges++;

				int ne_st = node_edges.bounds[cur_idx];
				int ne_en = node_edges.bounds[cur_idx+1] = ne_st + n_edges;

				auto set_ne = [&](int ne, std::array<int, 2> nvs, std::array<int, 2> nds, std::array<int, 4> rot_adjs) -> void {
					node_edges.dat[ne].node = cur_idx;
					node_edges.dat[ne].nvs = nvs;
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
						vert_pos_buf[node_verts.dat[nv].vert] = nv;
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
					auto map_rot_edge = [&](int ve) -> std::array<int, 4> {
						if (!with_planarity || !planar) return {-1, -1, -1, -1};
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

			auto start_child = [&]() -> int {
				auto& [cur_idx, ch_idx, ch_en, cur_nv, cur_ne] = stk.back();
				assert(ch_idx < ch_en);
				int nxt_item = ch.dat[ch_idx];
				int nxt_idx = nxt_unassigned_idx;
				ch.dat[ch_idx] = nxt_idx;
				par[nxt_idx] = cur_idx;
				int nxt_ne = node_edges.bounds[nxt_idx];
				if (nxt_item < 1 + NV) {
					vert_par_nv[nxt_idx] = cur_nv++;
				} else if (types[cur_idx] != node_type::F && types[cur_idx] != node_type::V) {
					node_edges.dat[cur_ne].twin_ne = nxt_ne;
					node_edges.dat[nxt_ne].twin_ne = cur_ne;
					cur_ne++;
				}

				ch_idx++;
				return nxt_item;
			};

			auto pop_item = [&]() -> void {
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
			assert(node_verts.bounds.back() == int(node_verts.dat.size()));
			assert(node_edges.bounds.back() == int(node_edges.dat.size()));
			assert(node_adj.bounds.back() == int(node_adj.dat.size()));

			// Rewrite node_vertices to the correct index
			for (auto& v : node_verts.dat) {
				v.vert = vert_index[v.vert];
			}

			return spqr_tree{
				std::move(vert_index),
				std::move(edge_index),
				std::move(par),
				std::move(subtree_end),
				std::move(types),
				std::move(orig_id),
				std::move(ch),
				std::move(node_verts),
				std::move(vert_par_nv),
				std::move(node_edges),
				std::move(node_adj),
				std::move(node_planar),
				std::move(ne_rot_adj),
			};
		}
	}
};

} // namespace wala
