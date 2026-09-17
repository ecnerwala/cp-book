#pragma once

#include <algorithm>
#include <vector>
#include <array>
#include <span>
#include <utility>
#include <cassert>
#include <ranges>

namespace wala {

template <typename T> struct csr {
	std::vector<int> bounds;
	std::vector<T> dat;
	std::span<T> operator [](int i) { return std::span<T>(dat.begin() + bounds[i], dat.begin() + bounds[i + 1]); }
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
	// Real edges are represented as a special Q nodes which each contain exactly 1 real edge and exactly 1 vedge.
	//
	// Traditionally, the SPQR tree is defined for each biconnected component,
	// but we will embed the SPQR decompositions inside the block-cut tree to get a (rooted) decomposition of the entire graph.
	//
	// As such, we will have a tree of "items", which consist of SPQR nodes, real vertices, and a special "forest root" item:
	//  - Each vertex will be a child of the topmost node which contains it (or the forest root).
	//  - Each block will be a subtree of nodes rooted at a Q edge, which is the child of one of its vertices.
	//
	// Item types:
	//   F - forest - a root node corresponding to the whole forst.
	//   V - vertex - not really a node, just there because they're mixed into the tree a la block/cut tree
	//   Q - real edge - has exactly 1 vedge and 1 real edg
	//   I - bridge - has exactly 1 vedge connecting to a bridge Q node
	//   O - self-loop - has exactly 1 vedge connecting to a self-loop Q node
	//   S - series - a cycle of >= 3 vedges; note that any 2 vertices of the cycle form a cut
	//   P - parallel - a parallel group of >= 3 vedges with the same endpoints
	//   R - rigid - a 3-vertex-connected component
	//
	// Degenerate blocks:
	//  - a block consisting of a self-loop is a Q node connected to an O node.
	//  - a block consisting of a bridge is a Q node connected to an I node.
	//  - a block consisting of exactly 2 parallel edges is represented by 2 glued Q nodes.
	enum class node_type : char {
		F = 'F', V = 'V', Q = 'Q', I = 'I', O = 'O', S = 'S', P = 'P', R = 'R'
	};

	std::vector<int> vert_index;
	std::vector<int> edge_index;

	std::vector<int> par;
	std::vector<int> subtree_end;
	std::vector<node_type> types;
	std::vector<int> orig_id;

	csr<int> ch;
	struct nv_t {
		int node_item;
		int vert_item;
	};
	csr<nv_t> node_verts;
	std::vector<int> item_nv;

	struct ne_t {
		int node_item;
		int twin_ne;
		// TODO: Should we store the twin node, the twin node type, and/or twin node type == Q?
		std::array<int, 2> nvs;
		std::array<int, 2> nds;
	};
	csr<ne_t> node_edges;

	struct nd_t {
		int node;
		int src_nv;
		int dest_nv;
		int ne;
		int twin_nd;
	};
	csr<nd_t> node_darts;

	static spqr_tree build(int NV, const std::vector<std::array<int, 2>>& edges) {
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
			auto dfs = [&](this auto&& self, int cur, int d, int prv_e) -> std::array<int, 2> {
				depth[cur] = d;
				std::array<int, 2> lowvals{d, d};
				for (auto [nxt, e] : adj[cur]) {
					if (e == prv_e) continue;
					if (depth[nxt] > d) continue;

					bool is_tree = depth[nxt] == -1;

					auto n_lowvals = is_tree ? self(nxt, d+1, e) : std::array<int, 2>{depth[nxt], d};

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
				}
				return lowvals;
			};
			for (int rt = 0; rt < NV; rt++) {
				if (depth[rt] == -1) {
					roots.push_back(rt);
					dfs(rt, 0, -1);
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
		auto get_side = []<typename T>(std::array<T, 2> a, bool dir) -> T{
			return dir ? a[1] : a[0];
		};

		struct item_list {
			std::array<int, 2> v{-1, -1};
			[[nodiscard]] bool empty() const { return v[0] == -1; }
		};
		std::vector<int> ch_nxt; ch_nxt.reserve(1 + NV + NE + NE); ch_nxt.assign(1 + NV + NE, -1);
		auto concat = [&](item_list a, item_list b) -> item_list {
			if (b.empty()) return a;
			if (a.empty()) return b;
			ch_nxt[a.v[1]] = b.v[0];
			return {{a.v[0], b.v[1]}};
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

		int tot_blocks = 0;
		int tot_self_loops = 0;

		{
			auto alloc_item = [&](node_type type) -> int {
				int item = int(item_vs.size());
				item_vs.push_back({});
				item_ch.push_back({});
				item_types.push_back(type);
				ch_nxt.push_back(-1);
				return item;
			};

			// Declare these here: most of our code will be in terms of v_start / top_depth, so we'll want to read these out
			std::vector<int> stack_verts(NV);
			std::vector<int> stack_dir(NV);

			auto make_vs = [&](int v_start, int top_depth) -> std::array<int, 2> {
				return set_sides(stack_dir[top_depth], stack_verts[top_depth], v_start);
			};

			int nxt_edge_idx = 0; // Counts backedges only
			std::vector<int> first_occurrence(NV); // First backedge to this depth

			struct tstack_t {
				int v_start;
				int top_depth;
				int first_idx;
				std::array<item_list, 2> spans;
			};
			auto make_tstack = [&](int v_start, int top_depth, int item) -> tstack_t {
				return { v_start, top_depth, nxt_edge_idx, set_sides(stack_dir[top_depth], unit_list(item), {}) };
			};
			auto merge_tstack = [&](tstack_t a, tstack_t b) -> tstack_t {
				return {
					a.v_start,
					min(a.top_depth, b.top_depth),
					a.first_idx,
					{concat(b.spans[0], a.spans[0]), concat(a.spans[1], b.spans[1])}
				};
			};

			auto maybe_unwrap = [&](tstack_t& t, node_type type) -> int {
				if (type == node_type::R) return alloc_item(type);

				assert(type == node_type::P || type == node_type::S);
				// TODO: This is the wrong dir for is-tree S-type checks, should we just dir in?
				//bool dir = stack_dir[t.top_depth];

				bool dir = t.spans[0].empty();
				assert(get_side(t.spans, !dir).empty());
				int item = get_side(t.spans, dir).v[0];
				assert(item == get_side(t.spans, dir).v[1]);
				if (item_types[item] == type) {
					t.spans = set_sides(dir, item_ch[item], {});
					return item;
				} else {
					return alloc_item(type);
				}
			};

			auto finish_tstack = [&](tstack_t& t, int item) {
				bool dir = stack_dir[t.top_depth];
				assert(get_side(t.spans, !dir).empty());

				item_vs[item] = make_vs(t.v_start, t.top_depth);
				item_ch[item] = get_side(t.spans, dir);
				t.spans = set_sides(dir, unit_list(item), {});
			};

			std::vector<tstack_t> tstack; tstack.reserve(NV + NE);
			auto pop_tstack = [&]() -> tstack_t {
				tstack_t res = tstack.back();
				tstack.pop_back();
				return res;
			};

			for (auto rt : roots) {
				[&](this auto&& self, int cur, int cur_depth) -> void {
					stack_verts[cur_depth] = cur;
					bool has_return_edge = false;

					for (auto [_, nxt, e, key] : outedges[cur]) {
						int lowval = key / 3 - 2; if (lowval < 0) lowval = cur_depth + ~lowval;
						int kind = key % 3;
						bool is_tree = kind != 1;
						bool is_type_1 = kind <= 1;

						// edge_dir convention: false is forwards, true is backwards.
						// That means that cur is on the edge_dir side and nxt is on the !edge_dir side.
						bool edge_dir = (lowval >= cur_depth ? false : !stack_dir[lowval]);
						stack_dir[cur_depth] = edge_dir;

						int orig_tstack = int(tstack.size());
						if (is_tree) {
							first_occurrence[cur_depth] = NE;
							self(nxt, cur_depth + 1);
						}

						if (lowval >= cur_depth) {
							item_vs[edge_item(e)] = {cur, -1};
							tot_blocks++;
							if (is_tree) {
								// Bridges and components
								if (lowval == cur_depth + 1) {
									// tstack.back() is currently just smuggling out the child vertex, prepend the bridge component
									// This is just a shortcut for allocating a full I-type tstack
									int item = alloc_item(node_type::I);
									item_vs[item] = make_vs(nxt, cur_depth);
									item_ch[edge_item(e)] = concat(unit_list(item), pop_tstack().spans[1]);
								} else {
									// tstack.end()[-2] is the vertex and tstack.end()[-1] is the backedge
									auto backedge = pop_tstack().spans[0];
									item_ch[edge_item(e)] = concat(backedge, pop_tstack().spans[1]);
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
							continue;
						}

						item_vs[edge_item(e)] = make_vs(nxt, cur_depth);

						assert(lowval < cur_depth);

						bool is_single = true;

						// make_q_node
						tstack_t cur_tstack;
						if (is_tree) {
							// The span lives on side edge_dir
							cur_tstack = make_tstack(nxt, cur_depth, edge_item(e));
							while (!tstack.empty() && tstack.back().top_depth >= cur_depth) {
								node_type type;
								if (tstack.back().top_depth > cur_depth) {
									// This will be an S node: merge the vertex in first
									cur_tstack = merge_tstack(pop_tstack(), cur_tstack);

									type = node_type::S;
								} else if (tstack.back().v_start == cur_tstack.v_start) {
									// This will be a P node
									type = node_type::P;
								} else {
									type = node_type::R;
								}
								auto nxt_tstack = pop_tstack();
								int item = maybe_unwrap(nxt_tstack, type);
								cur_tstack = merge_tstack(nxt_tstack, cur_tstack);
								finish_tstack(cur_tstack, item);
							}
							while (cur_tstack.first_idx > first_occurrence[cur_depth]) {
								is_single = false;
								cur_tstack = merge_tstack(pop_tstack(), cur_tstack);
							}

							if (has_return_edge || is_type_1) {
								// NB: tstack[orig_size] is the vertex and tstack[orig_size+1] is the backedge; maybe we should reverse them?
								assert(int(tstack.size()) >= orig_tstack + 2);

								int item;
								if (is_type_1) {
									assert(int(tstack.size()) == orig_tstack + 2);

									auto nxt_tstack = pop_tstack();
									item = maybe_unwrap(nxt_tstack, is_single ? node_type::S : node_type::R);
									cur_tstack = merge_tstack(nxt_tstack, cur_tstack);

									cur_tstack = merge_tstack(pop_tstack(), cur_tstack);
								} else {
									// Just to silence a warning
									item = -1;

									while (int(tstack.size()) > orig_tstack) {
										cur_tstack = merge_tstack(pop_tstack(), cur_tstack);
									}
								}

								cur_tstack.v_start = cur;
								assert(cur_tstack.top_depth == lowval);

								// Fold everything to the correct side now that we're leaving the child.
								// The entire subtree should go to the !edge_dir side.
								cur_tstack.spans = set_sides(!edge_dir, concat(cur_tstack.spans[0], cur_tstack.spans[1]), {});

								// TODO: There's some planarity folding to do here

								if (is_type_1) {
									finish_tstack(cur_tstack, item);
									is_single = true;
								}
							}
						} else {
							assert(is_type_1);
							// The span lives on side !edge_dir
							cur_tstack = make_tstack(cur, lowval, edge_item(e));
							nxt_edge_idx++;
							setmin(first_occurrence[lowval], cur_tstack.first_idx);
						}

						// NB: We can do this check in lots of ways, maybe there's a cleaner check
						if (is_type_1 && has_return_edge && tstack.back().v_start == cur && tstack.back().top_depth == lowval) {
							// This will be a P node
							auto nxt_tstack = pop_tstack();
							int item = maybe_unwrap(nxt_tstack, node_type::P);
							cur_tstack = merge_tstack(nxt_tstack, cur_tstack);
							finish_tstack(cur_tstack, item);
						}

						tstack.push_back(cur_tstack);

						if (!has_return_edge) {
							// Throw cur_vert_node onto the tstack so it'll get interleaved correctly
							tstack_t cur_vert_node = make_tstack(cur, cur_depth, vert_item(cur));

							assert(!tstack.empty());
							if (is_type_1) {
								// Insert it underneath the backedge
								cur_vert_node.first_idx = tstack.back().first_idx;
								std::swap(tstack.back(), cur_vert_node);
								tstack.push_back(cur_vert_node);
							} else if (is_single) {
								// This could be an S node, so leave it separate
								tstack.push_back(cur_vert_node);
							} else {
								// Just eagerly merge it to avoid a later bad finish_tstack
								tstack.back() = merge_tstack(tstack.back(), cur_vert_node);
							}
							has_return_edge = true;
						}
					}
					if (!has_return_edge) {
						// Either our parent is a bridge edge, or we're just a root.
						// We'll just leave it on tstack for future cleanup, it'll just get popped of immediately.
						// edge_dir == !stack_dir[lowval == cur_depth - 1] == true
						stack_dir[cur_depth] = true;
						tstack.push_back(make_tstack(cur, cur_depth, vert_item(cur)));
					}
				}(rt, 0);

				auto component_val = pop_tstack();
				item_ch[ROOT_ITEM] = concat(item_ch[ROOT_ITEM], component_val.spans[1]);
			}
		}

		// Phase 3: relabel the full tree in preorder
		// We have several id spaces:
		// * items are in preorder
		// * node_verts (nv's) are nodes-vertex pairs, given in node-item order then bipolar orientation vertex order
		// * node_edges (ne's) are node-vedge pairs, given in node-item order then a bipolar orientation vedge order
		// * node_darts (nd's) are node-dart (half-vedge) pairs, given in attached nv order
		//
		// Children of a node will be sorted in s-t order.
		// Specifically vertices are sorted, and edges are guaranteed to satisfy the strong "dominance" partial order:
		// if a.nvs[0] <= b.nvs[0] and a.nvs[1] <= b.nvs[1], then a <= b. (In practice, we'll sort by midpoint.)
		//
		// All id's are item indices unless clearly nv/ne/nd id's.
		//
		// * original verts and original edges can be converted to items as vert_item / edge_item
		int tot_items = int(item_types.size());
		{
			std::vector<int> vert_index(NV, -1);
			std::vector<int> edge_index(NE, -1);

			std::vector<int> par; par.reserve(tot_items);
			std::vector<int> subtree_end; subtree_end.reserve(tot_items);
			std::vector<node_type> types; types.reserve(tot_items);
			std::vector<int> orig_id; orig_id.reserve(tot_items);

			csr<int> ch;
			ch.bounds.reserve(tot_items + 1);
			ch.bounds.push_back(0);
			ch.dat.reserve(tot_items - 1);

			// Each node is a child, and additionally most non-block node has 2 cap verts; blocks have 1, and O nodes have 1
			int tot_node_verts = NV + (tot_items - 1 - NV) * 2 - tot_blocks - tot_self_loops;
			csr<nv_t> node_verts;
			node_verts.bounds.reserve(tot_items + 1);
			node_verts.bounds.push_back(0);
			node_verts.dat.reserve(tot_node_verts);
			std::vector<int> item_nv(tot_items, -1);

			int tot_node_edges = (tot_items - 1 - NV - tot_blocks) * 2;
			csr<ne_t> node_edges;
			node_edges.bounds.reserve(tot_items + 1);
			node_edges.bounds.push_back(0);
			node_edges.dat.reserve(tot_node_edges);

			csr<nd_t> node_darts;
			node_darts.bounds.reserve(tot_node_verts * 2 + 1);
			node_darts.bounds.push_back(0);
			node_darts.dat.reserve(tot_node_edges * 2);

			std::vector<int> vert_pos_buf(NV, -1);
			std::vector<int> cnts_buf(2 * NV, -1);
			std::vector<std::pair<int, int>> ch_buf(tot_items);

			[&](this auto&& self, int cur_item, int par_idx) -> void {
				int cur_idx = int(par.size());
				par.push_back(par_idx);
				subtree_end.push_back(-1);
				node_type cur_type = item_types[cur_item];
				types.push_back(cur_type);
				if (cur_type == node_type::F) {
					assert(cur_item == 0);
					orig_id.push_back(-1);
				} else if (cur_type == node_type::V) {
					assert(1 <= cur_item && cur_item < 1 + NV);
					int orig_vert = cur_item - 1;
					orig_id.push_back(orig_vert);
					vert_index[orig_vert] = cur_idx;
				} else if (cur_type == node_type::Q) {
					assert(1 + NV <= cur_item && cur_item < 1 + NV + NE);
					int orig_edge = cur_item - 1 - NV;
					orig_id.push_back(orig_edge);
					edge_index[orig_edge] = cur_idx;
				} else {
					assert(1 + NV + NE <= cur_item);
					orig_id.push_back(-1);
				}

				// HACK: Fill ch and vert_items in with orig items / orig verts for now,
				// because we don't have the final item id's yet.
				assert(ch.bounds.back() == int(ch.dat.size()));
				int ch_st = int(ch.dat.size());
				int nv_st = int(node_verts.dat.size());
				if (item_vs[cur_item][0] != -1) {
					node_verts.dat.push_back({cur_idx, item_vs[cur_item][0]});
				}
				for (int nxt_item = item_ch[cur_item].v[0]; nxt_item != -1; nxt_item = ch_nxt[nxt_item]) {
					ch.dat.push_back(nxt_item);
					assert(nxt_item >= 1);
					if (nxt_item < 1 + NV) {
						node_verts.dat.push_back({cur_idx, nxt_item - 1});
					}
					if (nxt_item == item_ch[cur_item].v[1]) {
						assert(ch_nxt[nxt_item] == -1);
					}
				}
				if (item_vs[cur_item][1] != -1) {
					node_verts.dat.push_back({cur_idx, item_vs[cur_item][1]});
				}
				int ch_en = int(ch.dat.size());
				ch.bounds.push_back(ch_en);
				int nv_en = int(node_verts.dat.size());
				node_verts.bounds.push_back(nv_en);

				int n_verts = nv_en - nv_st;

				bool is_node = cur_type != node_type::F && cur_type != node_type::V;
				bool has_cap = is_node && !(cur_type == node_type::Q && ch_en - ch_st > 0);

				for (int i = 0; i < n_verts; i++) {
					vert_pos_buf[node_verts.dat[i + nv_st].vert_item] = i;
				}
				if (cur_type == node_type::R) {
					// Bucketsort the children by the midpoint if necessary
					cnts_buf.assign(n_verts * 2 - 1, 0);
					ch_buf.clear();

					for (int i = ch_st; i < ch_en; i++) {
						int n = ch.dat[i];
						assert(n >= 1);
						int loc;
						if (n < 1 + NV) {
							loc = 2 * vert_pos_buf[n-1];
						} else {
							loc = vert_pos_buf[item_vs[n][0]] + vert_pos_buf[item_vs[n][1]];
						}
						ch_buf.emplace_back(loc, n);
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
				}

				// Fill in node_edges and node_darts
				int ne_st = int(node_edges.dat.size());
				if (is_node) {
					auto insert_ne = [&](int item) -> void {
						auto [v0, v1] = item_vs[item];
						if (v1 == -1) v1 = v0;
						node_edges.dat.push_back({
							cur_idx,
							// NB: Fill this in later
							-1,
							{nv_st + vert_pos_buf[v0], nv_st + vert_pos_buf[v1]},
							{-1, -1},
						});
					};
					if (has_cap) {
						insert_ne(cur_item);
					}
					for (int i = ch_st; i < ch_en; i++) {
						int n = ch.dat[i];
						assert(n >= 1);
						if (n < 1 + NV) continue;
						insert_ne(n);
					}
				}
				int ne_en = int(node_edges.dat.size());
				node_edges.bounds.push_back(ne_en);

				int n_edges = ne_en - ne_st;

				// Fill in nd
				assert(int(node_darts.bounds.size()) == 2 * nv_st + 1);
				node_darts.bounds.resize(2 * nv_en + 1, 0);
				assert(int(node_darts.dat.size()) == 2 * ne_st);
				node_darts.dat.resize(2 * ne_en);
				if (cur_type == node_type::P) {
					// Special case: tiebreak the parallel edges so they're reversed
					assert(nv_en - nv_st == 2);
					node_darts.bounds[2 * nv_st + 1] = 2 * ne_st + 0 * n_edges;
					node_darts.bounds[2 * nv_st + 2] = 2 * ne_st + 1 * n_edges;
					node_darts.bounds[2 * nv_st + 3] = 2 * ne_st + 2 * n_edges;
					node_darts.bounds[2 * nv_st + 4] = 2 * ne_st + 2 * n_edges;
					for (int i = ne_st; i < ne_en; i++) {
						int nd0 = 2 * ne_st + (i - ne_st);
						int nd1 = 2 * ne_en - 1 - (i - ne_st);
						assert(node_edges.dat[i].nvs[0] == nv_st);
						assert(node_edges.dat[i].nvs[1] == nv_st+1);
						node_edges.dat[i].nds = {nd0, nd1};
						node_darts.dat[nd0] = {
							cur_idx,
							nv_st,
							nv_st+1,
							i,
							nd1,
						};
						node_darts.dat[nd1] = {
							cur_idx,
							nv_st+1,
							nv_st,
							i,
							nd0,
						};
					}
				} else {
					// TODO: Could do other faster things for S nodes / Q nodes / whatever?
					// I'm pretty sure this logic is fine for self-loops (if a little nonsensical).
					for (int i = ne_st; i < ne_en; i++) {
						auto nvs = node_edges.dat[i].nvs;
						// Add to the counts
						node_darts.bounds[2 * nvs[0] + 2]++;
						node_darts.bounds[2 * nvs[1] + 1]++;
					}
					{
						int off = 2 * ne_st;
						for (int i = 2 * nv_st + 1; i <= 2 * nv_en; i++) {
							off += std::exchange(node_darts.bounds[i], off);
						}
						assert(off == 2 * ne_en);
					}
					for (int i = ne_st; i < ne_en; i++) {
						auto nvs = node_edges.dat[i].nvs;
						int nd0 = node_darts.bounds[2 * nvs[0] + 2]++;
						int nd1 = node_darts.bounds[2 * nvs[1] + 1]++;
						node_edges.dat[i].nds = {nd0, nd1};
						node_darts.dat[nd0] = {
							cur_idx,
							nvs[0],
							nvs[1],
							i,
							nd1,
						};
						node_darts.dat[nd1] = {
							cur_idx,
							nvs[1],
							nvs[0],
							i,
							nd0,
						};
					}
				}

				// Safe because ch is reserved to the right size
				{
					int cur_nv = nv_st + (item_vs[cur_item][0] != -1);
					int cur_ne = ne_st + has_cap;
					for (int i = ch_st; i < ch_en; i++) {
						// The index of the next node_edge if it exists
						int nxt_item = ch.dat[i];
						int nxt_idx = int(par.size());
						ch.dat[i] = nxt_idx;
						int nxt_ne = int(node_edges.dat.size());
						self(nxt_item, cur_idx);
						if (nxt_item < 1 + NV) {
							item_nv[nxt_item] = cur_nv++;
						} else if (is_node) {
							node_edges.dat[cur_ne].twin_ne = nxt_ne;
							node_edges.dat[nxt_ne].twin_ne = cur_ne;
							cur_ne++;
						}
					}
					cur_nv += (item_vs[cur_item][1] != -1);
					assert(cur_nv == nv_en);
					assert(cur_ne == ne_en);
					subtree_end[cur_idx] = int(subtree_end.size());
				}
			}(ROOT_ITEM, -1);

			assert(int(par.size()) == tot_items);

			// Rewrite node_vertices to the correct index
			for (auto& v : node_verts.dat) {
				v.vert_item = vert_index[v.vert_item];
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
				std::move(item_nv),
				std::move(node_edges),
				std::move(node_darts),
			};
		}
	}
};
} // namespace wala
