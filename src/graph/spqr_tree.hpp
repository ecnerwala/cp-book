#pragma once

#include <vector>
#include <array>
#include <span>
#include <utility>
#include <cassert>

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
	enum class node_type : char {
		Q = 'Q', I = 'I', O = 'O', S = 'S', P = 'P', R = 'R'
	};

	static spqr_tree build(int NV, const std::vector<std::array<int, 2>>& edges) {
		// TODO: Figure out the best way to specify roots; maybe accept a permutation of "root priority"?

		spqr_tree tree;

		int NE = int(edges.size());

		std::vector<int> roots; roots.reserve(NV);
		struct outedge_t { int src, dest; int e; int key; };
		csr<outedge_t> ch;

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

			std::vector<outedge_t> outedges; outedges.reserve(NE);
			// Return the 2 lowvals from this subtree
			auto dfs = [&](this auto&& self, int cur, int d, int prv_e) -> std::array<int, 2> {
				depth[cur] = d;
				std::array<int, 2> lowvals{d, d};
				for (int z = 0; z < int(adj[cur].size()); z++) {
					auto [nxt, e] = adj[cur][z];
					if (e == prv_e) continue;
					if (depth[nxt] > depth[cur]) continue;

					auto n_lowvals = depth[nxt] == -1 ? self(nxt, d+1, e) : std::array<int, 2>{depth[nxt], d};

					{
						// Extra bit is 0 for type-1 children, 1 for backedges, 2 for children with lowval2
						// Bridges have lowval -2 (kind 0), and components loops have lowval -1 (components are kind 0, loops are kind 1)
						// We don't really need to distinguish backedges vs type-1 children, but do it just for fun?
						int lowval = n_lowvals[0];
						if (lowval >= d) lowval = ~(lowval - d);
						int kind = 2 * (n_lowvals[1] < d) + (depth[nxt] <= depth[cur]);
						outedges.push_back({cur, nxt, e, 3 * (lowval + 2) + kind});
					}

					// Keep the 2 distinct mins
					if (n_lowvals[0] < lowvals[0]) lowvals = {n_lowvals[0], std::min(n_lowvals[1], lowvals[0])};
					else lowvals[1] = std::min(lowvals[1], n_lowvals[0] == lowvals[0] ? n_lowvals[1] : n_lowvals[0]);
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
			for (auto edge : outedges) by_key_builder.count(edge.key);
			by_key_builder.allocate();
			for (auto edge : outedges) by_key_builder.push(edge.key) = edge;
			csr<outedge_t> by_key = std::move(by_key_builder).finalize();

			csr_builder<outedge_t> by_src_builder(NV);
			// Hack to reuse memory
			by_src_builder.dat = std::move(outedges);
			for (auto edge : by_key.dat) by_src_builder.count(edge.src);
			by_src_builder.allocate();
			for (auto edge : by_key.dat) by_src_builder.push(edge.src) = edge;
			ch = std::move(by_src_builder).finalize();
		}


		// Phase 2: do the big ear-decomposition-like walk

		// We're going to build a tree of all SPQR *nodes* + all original *vertices* (collectively *items*).
		// Vertices will hang off the first SPQR node containing them, and blocks will be rooted at a topmost Q node for the top edge.

		// As we build, we will represent the Euler tour of nodes/vertices as a linked list.
		constexpr int ROOT_ITEM = 0;
		auto vert_item = [&](int v) -> int { return 1 + v; };
		auto edge_item = [&](int e) -> int { return 1 + NV + e; };

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

		struct item_t {
			std::array<int, 2> vs;
			item_list ch;
			node_type type;
		};
		std::vector<item_t> items; items.reserve(1 + NV + 2 * NE); items.resize(1 + NV + NE);
		auto alloc_item = [&](item_t dat) -> int {
			int item = int(items.size());
			items.push_back(dat);
			ch_nxt.push_back(-1);
			return item;
		};

		// Declare these here: most of our code will be in terms of v_start / top_depth, so we'll want to read these out
		std::vector<int> stack_verts(NV);
		std::vector<int> stack_dir(NV);

		auto make_node = [&](int v_start, int top_depth, node_type type, item_list inner = item_list{}) -> item_t {
			bool dir = stack_dir[top_depth];
			std::array<int, 2> vs{};
			vs[dir] = stack_verts[top_depth];
			vs[!dir] = v_start;
			return {vs, inner, type};
		};

		int nxt_edge_idx = 0; // Counts backedges only
		std::vector<int> first_occurrence(NV); // First backedge to this depth

		struct tstack_t {
			int v_start;
			int top_depth;
			int first_idx;
			std::array<item_list, 2> spans;
		};
		std::vector<tstack_t> tstack; tstack.reserve(NV + NE);
		auto make_tstack = [&](int v_start, int top_depth, int item) -> tstack_t {
			tstack_t t{ v_start, top_depth, nxt_edge_idx, {} };
			t.spans[stack_dir[top_depth]] = unit_list(item);
			return t;
		};
		auto merge_tstack = [&](tstack_t a, tstack_t b) -> tstack_t {
			return {
				a.v_start,
				std::min(a.top_depth, b.top_depth),
				a.first_idx,
				{concat(b.spans[0], a.spans[0]), concat(a.spans[1], b.spans[1])}
			};
		};

		auto maybe_unwrap = [&](tstack_t& t, node_type type) -> int {
			assert(type == node_type::P || type == node_type::S);
			// TODO: This is the wrong dir for is-tree S-type checks, should we just dir in?
			//bool dir = stack_dir[t.top_depth];

			bool dir = t.spans[0].empty();
			assert(t.spans[!dir].empty());
			int item = t.spans[dir].v[0];
			assert(item == t.spans[dir].v[1]);
			if (items[item].type == type) {
				t.spans[dir] = items[item].ch;
				return item;
			} else {
				return -1;
			}
		};

		auto finish_tstack = [&](tstack_t t, node_type type, int item) -> tstack_t {
			bool dir = stack_dir[t.top_depth];
			assert(t.spans[!dir].empty());

			auto item_dat = make_node(t.v_start, t.top_depth, type, t.spans[dir]);
			if (item == -1) item = alloc_item(item_dat);
			else items[item] = item_dat;
			t.spans[dir] = {item, item};
			return t;
		};

		// TODO: What node type
		items[ROOT_ITEM] = {{-1, -1}, item_list{}, node_type::O};
		for (auto rt : roots) {
			[&](this auto&& self, int cur, int cur_depth) -> void {
				stack_verts[cur_depth] = cur;
				bool has_return_edge = false;

				// TODO: node_type::V?
				items[vert_item(cur)] = {{cur, cur}, item_list{}, node_type::O};

				for (auto [_, nxt, e, key] : ch[cur]) {
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

					items[edge_item(e)] = make_node(nxt, cur_depth, node_type::Q);

					if (lowval >= cur_depth) {
						if (is_tree) {
							// Bridges and components
							if (lowval == cur_depth + 1) {
								// tstack.back() is currently just smuggling out the child vertex, prepend the bridge component
								// This is just a shortcut for allocating a full I-type tstack
								int item = alloc_item(make_node(nxt, cur_depth, node_type::I));
								items[edge_item(e)].ch = concat(unit_list(item), tstack.back().spans[1]);
								tstack.pop_back();
							} else {
								// tstack.end()[-2] is the vertex and tstack.end()[-1] is the backedge
								items[edge_item(e)].ch = concat(tstack.end()[-1].spans[0], tstack.end()[-2].spans[1]);
								tstack.pop_back();
								tstack.pop_back();
							}
						} else {
							// self loops
							assert(nxt == cur);
							int item = alloc_item(make_node(nxt, cur_depth, node_type::O));
							items[edge_item(e)].ch = unit_list(item);
						}
						items[vert_item(cur)].ch = concat(items[vert_item(cur)].ch, unit_list(edge_item(e)));
						continue;
					}

					assert(lowval < cur_depth);

					bool is_single = true;

					// make_q_node
					tstack_t cur_tstack;
					if (is_tree) {
						// The span lives on side edge_dir
						cur_tstack = make_tstack(nxt, cur_depth, edge_item(e));
						while (!tstack.empty() && tstack.back().top_depth >= cur_depth) {
							if (tstack.back().top_depth > cur_depth) {
								// This will be an S node: merge the vertex in first
								cur_tstack = merge_tstack(tstack.back(), cur_tstack);
								tstack.pop_back();

								int item = maybe_unwrap(tstack.back(), node_type::S);
								cur_tstack = finish_tstack(merge_tstack(tstack.back(), cur_tstack), node_type::S, item);
								tstack.pop_back();
							} else if (tstack.back().v_start == cur_tstack.v_start) {
								// This will be a P node
								int item = maybe_unwrap(tstack.back(), node_type::P);
								cur_tstack = finish_tstack(merge_tstack(tstack.back(), cur_tstack), node_type::P, item);
								tstack.pop_back();
							} else {
								cur_tstack = finish_tstack(merge_tstack(tstack.back(), cur_tstack), node_type::R, -1);
								tstack.pop_back();
							}
						}
						while (cur_tstack.first_idx > first_occurrence[cur_depth]) {
							is_single = false;
							cur_tstack = merge_tstack(tstack.back(), cur_tstack);
							tstack.pop_back();
						}

						if (has_return_edge || is_type_1) {
							// NB: tstack[orig_size] is the vertex and tstack[orig_size+1] is the backedge; maybe we should reverse them?
							assert(int(tstack.size()) >= orig_tstack + 2);
							if (is_type_1) assert(int(tstack.size()) == orig_tstack + 2);

							int item;
							if (is_type_1 && is_single) {
								item = maybe_unwrap(tstack.back(), node_type::S);
							} else {
								item = -1;
							}

							while (int(tstack.size()) > orig_tstack) {
								cur_tstack = merge_tstack(tstack.back(), cur_tstack);
								tstack.pop_back();
							}

							cur_tstack.v_start = cur;
							assert(cur_tstack.top_depth == lowval);

							// Fold everything to the correct side now that we're leaving the child.
							// The entire subtree should go to the !edge_dir side.
							cur_tstack.spans[!edge_dir] = concat(cur_tstack.spans[0], cur_tstack.spans[1]);
							cur_tstack.spans[edge_dir] = item_list{};

							// TODO: There's some planarity folding to do here

							if (is_type_1) {
								cur_tstack = finish_tstack(cur_tstack, is_single ? node_type::S : node_type::R, item);
								is_single = true;
							}
						}
					} else {
						assert(is_type_1);
						// The span lives on side !edge_dir
						cur_tstack = make_tstack(cur, lowval, edge_item(e));
						nxt_edge_idx++;
						first_occurrence[lowval] = std::min(first_occurrence[lowval], cur_tstack.first_idx);
					}

					// NB: We can do this check in lots of ways, maybe there's a cleaner check
					if (is_type_1 && has_return_edge && tstack.back().v_start == cur && tstack.back().top_depth == lowval) {
						// This will be a P node
						int item = maybe_unwrap(tstack.back(), node_type::P);
						cur_tstack = finish_tstack(merge_tstack(tstack.back(), cur_tstack), node_type::P, item);
						tstack.pop_back();
					}

					tstack.push_back(cur_tstack);

					if (!has_return_edge) {
						// Throw cur_vert_node onto the tstack so it'll get interleaved correctly
						tstack_t cur_vert_node = make_tstack(cur, cur_depth, vert_item(cur));

						assert(!tstack.empty());
						if (is_type_1) {
							// Insert it underneath the backedge
							cur_vert_node.first_idx = tstack.back().first_idx;
							tstack.insert(tstack.end() - 1, cur_vert_node);
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

			auto component_val = tstack.back();
			tstack.pop_back();
			items[ROOT_ITEM].ch = concat(items[ROOT_ITEM].ch, component_val.spans[1]);
		}

		return tree;
	}
};
} // namespace wala
