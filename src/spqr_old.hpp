#pragma once

#include <iostream>
#include <array>
#include <vector>
#include <cassert>

/** Block-Cut/SPQR Tree
 *
 *  A class to compute and represent the block-cut and SPQR trees of a static
 *  undirected (multi)graph.
 *
 *  In this tree, virtual edges (vedges) are numbered 2*k and 2*k+1, with 2*k+1
 *  in the "shallower" node.
 *
 *  Every real edge is stored in a Q node. The Q node for edge e is numbered e
 *  and contains virtual edge 2*e; virtual edge 2*e+1 is the "mirror" of edge e
 *  in the rest of the tree. Thus, you can test that virtual edge corresponds
 *  to a single real edge using `v < 2*NE` (and `(v & 1)`).
 *
 *  The tree is specified unrooted, but is implicitly rooted: vectors are
 *  ordered so that parents always come last. Blocks are numbered in preorder,
 *  and virtual edges are numbered in postorder (besides Q node virtual edges,
 *  see above).
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
 *  The types have the following relationships:
 *   * Components contain Vertices
 *   * Components contain Blocks contain Nodes contain VirtualEdges
 *   * Vertices are adjacent to Blocks, Nodes, and VirtualEdges (many-to-many)
 *  Some of these relationships are denormalized.
 */
struct spqr_tree_old {
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

	struct vertex_t {
		int component = -1;
		std::vector<int> blocks;
		std::vector<int> nodes;
	};
	std::vector<vertex_t> vertices;

	// Just the endpoints
	std::vector<std::array<int, 2>> edges;

	struct vedge_t {
		std::array<int, 2> v;
		int block = -1;
		int node = -1;
		bool is_tree = false;
	};
	std::vector<vedge_t> vedges;

	struct component_t {
		std::vector<int> vertices;
		std::vector<int> blocks;
	};
	std::vector<component_t> components;

	struct block_t {
		int component = -1;
		std::vector<int> vertices;
		std::vector<int> edges;
		std::vector<int> vedges;
		std::vector<int> nodes;
	};
	std::vector<block_t> blocks;

	enum class node_type : char {
		Q = 'Q', I = 'I', O = 'O', S = 'S', P = 'P', R = 'R'
	};
	struct node_t {
		node_type type;
		int block = -1;
		std::vector<int> vertices;
		std::vector<int> vedges; // nonnegative is virtual, negative is real
	};
	std::vector<node_t> nodes;

	// Internal vectors, you can use them if you want
	std::vector<std::vector<std::pair<int, int>>> adj;
	std::vector<int> depth;
	std::vector<int> lowval;
	std::vector<int> lowval2;

private:
	std::vector<std::vector<int>> lowval_buckets;

	std::pair<int, int> dfs_lowval(int cur, int d, int prvE) {
		depth[cur] = d;
		int v1 = d, v2 = d;
		int num_blocks = 0;
		for (auto [nxt, e] : adj[cur]) {
			if (e == prvE) continue;
			if (depth[nxt] == -1) {
				edges[e] = {cur, nxt};

				const auto [n1, n2] = dfs_lowval(nxt, d+1, e);

				if (n1 == d + 1) {
					num_blocks++;
					lowval_buckets[0].push_back(e);
				} else if (n1 == d) {
					num_blocks++;
					lowval_buckets[1].push_back(e);
				} else {
					lowval_buckets[2 + n1 * 2 + (n2 < d)].push_back(e);
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
				edges[e] = {cur, nxt};

				int nd = depth[nxt];
				if (nd == d) {
					// Self-loop
					num_blocks++;
					lowval_buckets[1].push_back(e);
				} else {
					lowval_buckets[2 + nd * 2].push_back(e);
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
		vertices[cur].blocks.reserve(num_blocks + (prvE != -1));
		return {lowval[cur] = v1, lowval2[cur] = v2};
	}

	void build_sorted_adj(int root) {
		{
			std::vector<int> deg(NV, 0);
			for (auto e : edges) {
				for (int v : e) {
					deg[v]++;
				}
			}
			adj = std::vector<std::vector<std::pair<int, int>>>(NV);
			for (int i = 0; i < NV; i++) {
				adj[i].reserve(deg[i]);
			}
		}
		for (int e = 0; e < int(edges.size()); e++) {
			auto [u, v] = edges[e];
			adj[u].push_back({v, e});
			if (u != v) {
				adj[v].push_back({u, e});
			}
		}

		depth = std::vector<int>(NV, -1);
		lowval = std::vector<int>(NV);
		lowval2 = std::vector<int>(NV);

		// Bucketed so that bridges come first, then BCCs, then children from shallowest to deepest lowval
		lowval_buckets = std::vector<std::vector<int>>(2 + NV * 2);

		for (int rt = 0; rt < NV; rt++) {
			if (root != -1 && rt != root) continue;
			if (depth[rt] != -1) continue;
			dfs_lowval(rt, 0, -1);
		}

		for (auto& v : adj) {
			v.clear();
		}
		for (const auto& bucket : lowval_buckets) {
			for (auto e : bucket) {
				adj[edges[e][0]].push_back({edges[e][1], e});
			}
		}

		// Free this; literally no one cares about this
		lowval_buckets = {};
	}

	void cleanup_internals() {
		adj = {};
		depth = {};
		lowval = {};
		lowval2 = {};
	}

	struct estack_t {
		int node;
		std::array<int, 2> vs;
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

	int make_node(node_type type) {
		int node = NN++;
		nodes.emplace_back();
		nodes[node].type = type;
		return node;
	}

	void finalize_node(int node, int ve0, int block) {
		vedges[ve0].node = node;
		nodes[node].vedges.push_back(ve0);

		// Often, this capacity is already there.
		nodes[node].vertices.reserve(nodes[node].vertices.size() + 2);
		for (auto v : vedges[ve0].v) {
			nodes[node].vertices.push_back(v);
			vertices[v].nodes.push_back(node);
		}

		nodes[node].block = block;
		blocks[block].nodes.push_back(node);
	}

	int finalize_estack(estack_t stk, int block) {
		int node = stk.node;
		auto vs = stk.vs;
		bool is_tree = stk.is_tree;

		std::array<int, 2> ves;
		if (node < NE) {
			assert(nodes[node].type == node_type::Q);
			ves[0] = 2*node;
			ves[1] = 2*node + 1;
		} else {
			ves[0] = NVE++;
			ves[1] = NVE++;
			vedges.resize(vedges.size() + 2);
		}
		auto [ve0, ve1] = ves;

		for (auto ve : ves) {
			vedges[ve].block = block;
			blocks[block].vedges.push_back(ve);
		}

		vedges[ve0].v = {vs[1], vs[0]};
		vedges[ve1].v = {vs[0], vs[1]};

		vedges[ve0].is_tree = !is_tree;
		vedges[ve1].is_tree = is_tree;

		assert(is_tree == (depth[vs[0]] > depth[vs[1]]));

		finalize_node(node, ve0, block);

		return ve1;
	}

	void push_estack(estack_t e_ins) {
		if constexpr (PRINT_DEBUG) std::cerr << "push_estack " << e_ins.node << ' ' << e_ins.vs[0] << '-' << e_ins.vs[1] << ' ' << e_ins.is_tree << '\n';
		estack.push_back(e_ins);
		int v = e_ins.vs[0];
		if (first_occurrence[v] == -1) {
			first_occurrence[v] = int(estack.size()) - 1;
		}
	}

	int push_estack_p(estack_t e_ins, int block) {
		if constexpr (PRINT_DEBUG) std::cerr << "push_estack_p " << e_ins.node << ' ' << e_ins.vs[0] << '-' << e_ins.vs[1] << ' ' << e_ins.is_tree << '\n';
		int node = estack.back().node;
		if (estack.back().type != node_type::P) {
			node = make_node(node_type::P);
			nodes[node].vedges.reserve(3);
			int ve = finalize_estack(estack.back(), block);

			vedges[ve].node = node;
			nodes[node].vedges.push_back(ve);

			estack.back().node = node;
			estack.back().type = node_type::P;
		}

		int ve = finalize_estack(e_ins, block);
		vedges[ve].node = node;
		nodes[node].vedges.push_back(ve);

		return node;
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
	
	std::pair<int, node_type> pop_estack_into_node(int idx, int exclude_vertex, int block) {
		assert(int(estack.size()) - idx > 1);
		bool is_S = (int(estack.size()) - idx == 2);
		if (is_S) {
			assert(nodes[estack.back().node].type != node_type::S);
		}
		bool should_reuse = (is_S && estack[idx].type == node_type::S);
		int node;
		if (should_reuse) {
			node = estack[idx].node;
		} else {
			node = make_node(is_S ? node_type::S : node_type::R);
			nodes[node].vedges.reserve(int(estack.size()) - idx + 1);
			// This reserves too much, but better too much than too little!
			nodes[node].vertices.reserve(int(estack.size()) - idx + 1);
		}
		for (int z = idx; z < int(estack.size()); z++) {
			prepare_pop_estack(z);
			if (should_reuse && z == idx) continue;
			int ve = finalize_estack(estack[z], block);
			vedges[ve].node = node;
			nodes[node].vedges.push_back(ve);

			if (estack[z].is_tree) {
				int v = estack[z].vs[0];
				if (v != exclude_vertex) {
					nodes[node].vertices.push_back(v);
					vertices[v].nodes.push_back(node);
				}
			}
		}
		estack.resize(idx);
		return {node, is_S ? node_type::S : node_type::R};
	}

	void dfs_spqr(int cur, int cur_block, int cur_low, int component) {
		for (auto [nxt, e] : adj[cur]) {
			if constexpr (PRINT_DEBUG) {
				std::cerr << "cur " << cur << '\n';
				std::cerr << "cur_block " << cur_block << '\n';
				std::cerr << "cur_low " << cur_low << '\n';
				std::cerr << "nxt " << nxt << '\n';
				std::cerr << "e " << e << '\n';
				std::cerr << '\n';
			}
			int orig_size = int(estack.size());
			assert(depth[nxt] <= depth[cur] + 1);
			if (depth[nxt] == depth[cur] || (depth[nxt] > depth[cur] && lowval[nxt] >= depth[cur])) {
				// Root of a new block

				int block = NB++;
				blocks.emplace_back();
				blocks[block].component = component;
				components[component].blocks.push_back(block);

				int node;
				if (nxt != cur) {
					// set cur_low to nxt
					dfs_spqr(nxt, block, nxt, component);

					// finalize the block
					if (int(estack.size()) == orig_size) {
						// bridge
						assert(lowval[nxt] == depth[nxt]);

						node = make_node(node_type::I);
					} else {
						assert(lowval[nxt] == depth[cur]);
						assert(int(estack.size()) == orig_size + 1);
						assert((estack.back().vs == std::array<int, 2>{cur, nxt}));
						assert(!estack.back().is_tree);

						node = estack.back().node;
						if (estack.back().type == node_type::Q) {
							node = make_node(node_type::O);
							nodes[node].vedges.reserve(2);
							int ve = finalize_estack(estack.back(), block);
							vedges[ve].node = node;
							nodes[node].vedges.push_back(ve);
						}
						pop_estack();

						assert(!tstack.empty() && tstack.back().idx == orig_size);
						tstack.pop_back();
					}
				} else {
					// self-loop
					node = make_node(node_type::O);
				}

				nodes[e].type = node_type::Q;
				blocks[block].edges.push_back(e);

				int ve = finalize_estack(estack_t{e, {nxt, cur}, node_type::Q, cur != nxt}, block);
				assert(ve == 2*e+1);
				if (nxt == cur) {
					// We pushed our current vertex/node pair twice, fix it here
					nodes[e].vertices.pop_back();
					vertices[cur].nodes.pop_back();
				}

				// Now cap off our node with this virtual edge
				finalize_node(node, ve, block);
				if (nxt == cur) {
					// We pushed our current vertex/node pair twice, fix it here
					nodes[node].vertices.pop_back();
					vertices[cur].nodes.pop_back();
				}

				vertices[cur].blocks.push_back(block);
				blocks[block].vertices.push_back(cur);
			} else {
				// Normal SPQR algorithm

				// We're not the root, though we may be the child of the root
				assert(lowval[cur] < depth[cur]);
				assert(cur_block != -1);
				assert(cur_low != -1);

				if (depth[nxt] > depth[cur]) {
					// child vertex
					dfs_spqr(nxt, cur_block, cur_low, component);

					if constexpr (PRINT_DEBUG) {
						std::cerr << "pre push down edge" << '\n';
						for (auto es : estack) {
							std::cerr << es.node << ' ' << es.vs[0] << '-' << es.vs[1] << ' ' << es.is_tree << '\n';
						}
						for (auto ts : tstack) {
							std::cerr << ts.idx << ' ' << ts.vstart << '-' << ts.top_depth << '\n';
						}
						std::cerr << '\n';
					}

					nodes[e].type = node_type::Q;
					blocks[cur_block].edges.push_back(e);

					// Before we insert this, we may have to fix the single-tree-edge tstack
					assert(!tstack.empty());
					int cur_depth = depth[cur];
					if (tstack.back().top_depth == cur_depth + 1) {
						assert(tstack.back().idx == int(estack.size())-1);
						tstack.back().top_depth = cur_depth;
					}

					// Insert this edge and check for type-2 splits
					for (estack_t e_ins{e, {nxt, cur}, node_type::Q, true}; true; ) {
						assert(int(estack.size()) > orig_size);
						if (estack.back().vs == std::array<int, 2>{cur, nxt}) {
							e_ins.node = push_estack_p(e_ins, cur_block);
							e_ins.type = node_type::P;

							assert(!tstack.empty());
							if (tstack.back().idx == int(estack.size())-1) {
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
						auto [node, type] = pop_estack_into_node(idx, nxt, cur_block);
						e_ins = estack_t{node, {nxt, cur}, type, true};
					}

					if (lowval2[nxt] >= depth[cur]) {
						// Handle type 1 split

						{

							int nnxt = estack[orig_size].vs[0];
							assert(depth[nnxt] == lowval[nxt]);
							nxt = nnxt;
						}

						// It's possible there are candidates before orig_size
						// on the stack, but those are doomed to fail, so
						// it's fine to pop them prematurely.
						while (tstack.back().idx > orig_size) tstack.pop_back();
						assert(!tstack.empty());
						if (tstack.back().idx == orig_size) {
							assert(tstack.back().vstart == cur_low && tstack.back().top_depth == depth[nxt]);
						}

						int idx = orig_size;
						auto [node, type] = pop_estack_into_node(idx, -1, cur_block);
						estack_t e_ins = estack_t{node, {nxt, cur}, type, false};
						if (!estack.empty() && estack.back().vs == e_ins.vs) {
							push_estack_p(e_ins, cur_block);
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
					// backedge
					nodes[e].type = node_type::Q;
					blocks[cur_block].edges.push_back(e);

					if constexpr (PRINT_DEBUG) {
						std::cerr << "pre push backedge" << '\n';
						for (auto es : estack) {
							std::cerr << es.node << ' ' << es.vs[0] << '-' << es.vs[1] << ' ' << es.is_tree << '\n';
						}
						for (auto ts : tstack) {
							std::cerr << ts.idx << ' ' << ts.vstart << '-' << ts.top_depth << '\n';
						}
						std::cerr << '\n';
					}

					estack_t e_ins{e, {nxt, cur}, node_type::Q, false};
					if (!estack.empty() && estack.back().vs == e_ins.vs) {
						push_estack_p(e_ins, cur_block);
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

				if constexpr (PRINT_DEBUG) {
					std::cerr << "post push" << '\n';
					for (auto es : estack) {
						std::cerr << es.node << ' ' << es.vs[0] << '-' << es.vs[1] << ' ' << es.is_tree << '\n';
					}
					for (auto ts : tstack) {
						std::cerr << ts.idx << ' ' << ts.vstart << '-' << ts.top_depth << '\n';
					}
					std::cerr << '\n';
				}
			}
		}

		if (cur_block != -1) {
			vertices[cur].blocks.push_back(cur_block);
			blocks[cur_block].vertices.push_back(cur);
		}
		vertices[cur].component = component;
		components[component].vertices.push_back(cur);
	}

	void build_spqr() {
		// Tree edges have v = {child, parent}
		// Backedges edges have v = {top, bottom}

		estack.reserve(NE);
		tstack.reserve(NE);
		first_occurrence = std::vector<int>(NV, -1);

		for (int rt = 0; rt < NV; rt++) {
			if (depth[rt] != 0) continue;
			int component = NC++;
			components.emplace_back();
			dfs_spqr(rt, -1, -1, component);
			assert(estack.empty());
			assert(tstack.empty());
		}
		estack = {};
		tstack = {};
		first_occurrence = {};
	}

public:
	spqr_tree_old() = default;
	// Pass in a list of (undirected) edges
	// Edges may get flipped in the output
	explicit spqr_tree_old(int NV_, std::vector<std::array<int, 2>> edges_, int root = -1) {
		NV = NV_;
		if (NV == 0) return;

		if (root != -1) {
			assert(0 <= root && root < NV);
		}

		edges = std::move(edges_);
		NE = int(edges.size());

		vertices.resize(NV);

		components.reserve(NV);
		blocks.reserve(NE);

		nodes.reserve(2 * NE);
		NN = NE;
		nodes.resize(NN);

		vedges.reserve(2 * (2 * NE));
		NVE = NE * 2;
		vedges.resize(NVE);

		build_sorted_adj(root);
		build_spqr();
		// cleanup_interals();
	}
};
