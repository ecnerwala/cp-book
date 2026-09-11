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
		struct outedge_t { int src, dest; int e; int val; };
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
			auto dfs = [&](this auto&& self, int cur, int d, int prvE) -> std::array<int, 2> {
				depth[cur] = d;
				std::array<int, 2> val{d, d};
				for (int z = 0; z < int(adj[cur].size()); z++) {
					auto [nxt, e] = adj[cur][z];
					if (e == prvE) continue;
					if (depth[nxt] > depth[cur]) continue;

					auto nval = depth[nxt] == -1 ? self(nxt, d+1, e) : std::array<int, 2>{depth[nxt], d};
					// Extra bit is 0 for type-1 children, 1 for backedges, 2 for children with lowval2
					// Bridges have lowval -2 (val -6), and components loops have lowval -1 (components are val -3, loops are val -2)
					// We don't really need to distinguish backedges vs type-1 children, but do it just for fun?

					int n_lowval = nval[0];
					if (n_lowval >= d) n_lowval = ~(n_lowval - d);

					int n_type = 2 * (nval[1] < d) + (depth[nxt] <= depth[cur]);

					outedges.push_back({cur, nxt, e, 3 * n_lowval + n_type});

					if (nval[0] < val[0]) val = {nval[0], std::min(nval[1], val[0])};
					else val[1] = std::min(val[1], nval[0] == val[0] ? nval[1] : nval[0]);
				}
				return val;
			};
			for (int rt = 0; rt < NV; rt++) {
				if (depth[rt] == -1) {
					roots.push_back(rt);
					dfs(rt, 0, -1);
				}
			}

			csr_builder<outedge_t> depth_builder(3*NV+6);
			for (auto edge : outedges) depth_builder.count(edge.val+6);
			depth_builder.allocate();
			for (auto edge : outedges) depth_builder.push(edge.val+6) = edge;
			csr<outedge_t> by_depth = std::move(depth_builder).finalize();

			csr_builder<outedge_t> ch_builder(NV);
			// Hack to reuse memory
			ch_builder.dat = std::move(outedges);
			for (auto edge : by_depth.dat) ch_builder.count(edge.src);
			ch_builder.allocate();
			for (auto edge : by_depth.dat) ch_builder.push(edge.src) = edge;
			ch = std::move(ch_builder).finalize();
		}


		// Phase 2: do the big ear-decomposition-like walk
		std::vector<int> st_nxt; st_nxt.reserve(2 * (1 + NV + 2 * NE)); st_nxt.assign(2 * (1 + NV + NE), -1);
		struct node_t {
			std::array<int, 2> vs;
			node_type type;
		};
		std::vector<node_t> nodes; nodes.reserve(2 * NE); nodes.resize(NE);
		struct st_list { std::array<int, 2> v{-1, -1}; };
		auto concat = [&](st_list a, st_list b) -> st_list {
			if (b.v[0] == -1) return a;
			if (a.v[0] == -1) return b;
			st_nxt[a.v[1]] = b.v[0];
			return {{a.v[0], b.v[1]}};
		};

		auto wrap_st_list = [&](int n, st_list a) -> st_list {
			if (a.v[0] == -1) st_nxt[2*n] = 2*n+1;
			else st_nxt[2*n] = a.v[0], st_nxt[a.v[1]] = 2*n+1;
			return {{2*n, 2*n+1}};
		};

		// edge_dir == false means cur then nxt
		auto alloc_node = [&](int cur, int nxt, bool edge_dir, node_type type, st_list contents) -> st_list {
			std::array<int, 2> vs;
			vs[edge_dir] = cur;
			vs[!edge_dir] = nxt;
			nodes.push_back({vs, type});
			int n = int(st_nxt.size())/2;
			st_nxt.push_back(-1);
			st_nxt.push_back(-1);
			return wrap_st_list(n, contents);
		};
		struct tstack_t {
			int v_start;
			int top_depth;
			int first_idx;
			int num_edges;
			std::array<st_list, 2> lst;
		};
		auto merge_tstack = [&](tstack_t a, tstack_t b) -> tstack_t {
			return {
				a.v_start,
				std::min(a.top_depth, b.top_depth),
				a.first_idx,
				a.num_edges + b.num_edges,
				{concat(b.lst[0], a.lst[0]), concat(a.lst[1], b.lst[1])}
			};
		};
		std::vector<tstack_t> tstack; tstack.reserve(std::max(1, NE));
		std::vector<int> cur_path(NV);
		std::vector<int> stack_dir(NV);
		std::vector<int> first_occurrence(NV);
		int nxt_edge_idx = 0;
		st_list all_comps{{0, 0}};
		for (auto rt : roots) {
			[&](this auto&& self, int cur, int cur_depth, int cur_lowval) -> void {
				cur_path[cur_depth] = cur;
				bool has_return_edge = false;
				st_list cur_subtree{};
				for (auto [_, nxt, e, val] : ch[cur]) {
					int lowval = (val + 6) / 3 - 2;
					bool is_tree = (val + 6) % 3 != 1;
					bool is_type_1 = ((val + 6) % 3) < 2;

					if (lowval < 0) lowval = cur_depth + ~lowval;

					// edge_dir convention: false is forwards, true is backwards.
					// That means that cur is on the edge_dir side and nxt is on the !edge_dir side.
					bool edge_dir = (val < 0 ? false : !stack_dir[lowval]);
					stack_dir[cur_depth] = edge_dir;

					auto make_node = [&](tstack_t t, bool dir) -> tstack_t {
						assert(t.num_edges >= 2);
						node_type type;
						if (t.num_edges == 2) {
							// TODO: figure this out, including reuse?
						} else {
							type = node_type::R;
						}

						assert(dir == stack_dir[t.top_depth]);
						assert(t.lst[!dir].v[0] == -1);
						t.lst[dir] = alloc_node(t.v_start, cur_path[t.top_depth], !dir, type, t.lst[dir]);
						t.num_edges = 1;
						return t;
					};

					int orig_tstack = int(tstack.size());
					if (is_tree) {
						first_occurrence[cur_depth] = NE;
						self(nxt, cur_depth + 1, lowval);
					}

					int e_n = 1 + NV + e;
					nodes[e].vs[edge_dir] = cur;
					nodes[e].vs[!edge_dir] = nxt;
					nodes[e].type = node_type::Q;

					if (val < 0) {
						st_list block_list;
						if (is_tree) {
							// Bridges and components
							if (lowval == cur_depth + 1) {
								// tstack.back() is currently just smuggling out the child vertex, prepend the bridge component
								block_list = concat(alloc_node(cur, nxt, edge_dir, node_type::I, st_list{}), tstack.back().lst[1]);
								tstack.pop_back();
							} else {
								// tstack.end()[-2] is the vertex and tstack.end()[-1] is the backedge
								block_list = concat(tstack.end()[-1].lst[0], tstack.end()[-2].lst[1]);
								tstack.pop_back();
								tstack.pop_back();
							}
						} else {
							// self loops
							assert(nxt == cur);
							block_list = alloc_node(cur, nxt, edge_dir, node_type::O, st_list{});
						}
						cur_subtree = concat(cur_subtree, wrap_st_list(e_n, block_list));
						continue;
					}

					assert(val >= 0);
					assert(lowval < cur_depth);

					// make_q_node
					tstack_t cur_tstack;
					if (is_tree) {
						// TODO: The is_tree construction here really matches what a backedge from nxt -> cur would look like, maybe we
						// should push it down?
						cur_tstack = {
							nxt,
							cur_depth,
							nxt_edge_idx++,
							1,
							{st_list{}, st_list{}}
						};
						cur_tstack.lst[edge_dir] = wrap_st_list(e_n, st_list{});
						bool must_merge_all = has_return_edge || is_type_1;
						while (!tstack.empty() && tstack.back().top_depth >= cur_depth) {
							cur_tstack = make_node(merge_tstack(tstack.back(), cur_tstack), edge_dir);
							tstack.pop_back();
						}
						while (cur_tstack.first_idx > first_occurrence[cur_depth]) {
							cur_tstack = merge_tstack(tstack.back(), cur_tstack);
							tstack.pop_back();
						}

						if (must_merge_all) {
							// Merge the rest
							while (int(tstack.size()) > orig_tstack) {
								cur_tstack = merge_tstack(tstack.back(), cur_tstack);
								tstack.pop_back();
							}

							assert(cur_tstack.top_depth == lowval);

							// Fold everything to the correct side now that we're leaving the child
							// The entire subtree should go to the !edge_dir side
							cur_tstack.lst[!edge_dir] = concat(cur_tstack.lst[0], cur_tstack.lst[1]);
							cur_tstack.lst[edge_dir] = st_list{};
							cur_tstack.v_start = cur;

							// TODO: Planarity has some logic here
							if (is_type_1) {
								// merge it into a single edge
								cur_tstack = make_node(cur_tstack, !edge_dir);
							}
						}
					} else {
						assert(is_type_1);
						cur_tstack = {
							cur,
							lowval,
							nxt_edge_idx++,
							1,
							{st_list{}, st_list{}}
						};
						cur_tstack.lst[!edge_dir] = wrap_st_list(e_n, st_list{});
						first_occurrence[lowval] = std::min(first_occurrence[lowval], cur_tstack.first_idx);
					}

					// NB: We can do this check in lots of ways, maybe there's a cleaner check
					if (is_type_1 && has_return_edge && tstack.back().num_edges == 1 && tstack.back().v_start == cur && tstack.back().top_depth == lowval) {
						tstack.back() = make_node(merge_tstack(tstack.back(), cur_tstack), !edge_dir);
					} else {
						tstack.push_back(cur_tstack);
					}

					if (!has_return_edge) {
						// Throw cur_vert_node onto the tstack so it'll get interleaved correctly
						tstack_t cur_vert_node{
							cur,
							cur_depth,
							// Copy the back's idx
							tstack.back().first_idx,
							0,
							{st_list{}, st_list{}}
						};
						cur_vert_node.lst[edge_dir] = wrap_st_list(1 + cur, cur_subtree);

						assert(!tstack.empty());
						if (is_type_1) {
							// Insert it underneath the backedge
							tstack.insert(tstack.end() - 1, cur_vert_node);
						} else {
							// Could insert it over, but if we merge this way we prevent some spurious merges
							tstack.back() = merge_tstack(tstack.back(), cur_vert_node);
						}
						has_return_edge = true;
					}
				}
				if (!has_return_edge) {
					// Either our parent is a bridge edge, or we're just a root; either way, we'll just leave it on tstack for future cleanup, it'll just get popped of immediately
					tstack_t cur_vert_node{
						cur,
						cur_depth,
						nxt_edge_idx,
						0,
						{st_list{}, st_list{}}
					};
					cur_vert_node.lst[1] = wrap_st_list(1 + cur, cur_subtree);
					tstack.push_back(cur_vert_node);
				}
			}(rt, 0, 0);
			auto component_val = tstack.back();
			tstack.pop_back();
			all_comps = concat(all_comps, component_val.lst[1]);
		}

		// Cap it off
		st_nxt[all_comps.v[1]] = 1;
		all_comps.v[1] = 1;

		return tree;
	}
};
} // namespace wala
