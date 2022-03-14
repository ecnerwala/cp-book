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
 *  indirectly `tree.vedges[ve].o_ve < NE`.
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
 *  intervals of Nodes, and so on.
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
		struct iterator {
			int v;
			int operator * () const { return v; }
			iterator& operator ++ () { ++v; return *this; }
			friend bool operator != (iterator a, iterator b) { return a.v != b.v; }
		};
		iterator begin() const { return iterator{st}; }
		iterator end() const { return iterator{en}; }
	};
	template <typename T, std::vector<T> spqr_tree::* array> struct array_range_t : public range_t {
		struct bound_array_range_t {
			typename std::vector<T>::const_iterator st, en;
			int size() const { return en - st; }
			auto begin() const { return st; }
			auto end() const { return en; }
		};
		bound_array_range_t bind(const spqr_tree& tree) const {
			return bound_array_range_t{(tree.*array).begin() + st, (tree.*array).begin() + en};
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
		array_range_t<block_t, &spqr_tree::blocks> blocks;
		array_range_t<node_t, &spqr_tree::nodes> nodes;
		array_range_t<vedge_t, &spqr_tree::vedges> vedges;
		array_range_t<int, &spqr_tree::component_vertices> component_vertices;
	};

	struct block_t {
		int component = -1;
		array_range_t<node_t, &spqr_tree::nodes> nodes;
		array_range_t<vedge_t, &spqr_tree::vedges> vedges;
		array_range_t<int, &spqr_tree::block_vertices> block_vertices;
	};

	enum class node_type : char {
		Q = 'Q', I = 'I', O = 'O', S = 'S', P = 'P', R = 'R'
	};

	struct node_t {
		node_type type;
		int component = -1;
		int block = -1;
		array_range_t<vedge_t, &spqr_tree::vedges> vedges;
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

private:
	// Negative e means new block; otherwise, it's 2*e + (has_nontrivial_lowval2)
	std::vector<std::vector<std::pair<int, int>>> adj;

	struct bucket_edge_t {
		int e; int cur, nxt;
	};
	std::vector<std::vector<bucket_edge_t>> lowval_buckets;

	std::pair<int, int> dfs_lowval(int cur, int d, int prvE) {
		depth[cur] = d;
		int v1 = d, v2 = d;
		for (auto [nxt, e] : adj[cur]) {
			if (e == prvE) continue;
			if (depth[nxt] == -1) {
				const auto [n1, n2] = dfs_lowval(nxt, d+1, e);

				if (n1 >= d) {
					lowval_buckets[0].push_back({~e, cur, nxt});
				} else {
					lowval_buckets[1 + n1 * 2 + (n2 < d)].push_back({2*e + (n2 < d), cur, nxt});
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
					lowval_buckets[0].push_back({~e, cur, nxt});
				} else {
					lowval_buckets[1 + nd * 2].push_back({2*e, cur, nxt});
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

		// Bucketed so that bridges come first, then BCCs, then children from shallowest to deepest lowval
		lowval_buckets = std::vector<std::vector<bucket_edge_t>>(1 + NV * 2);

		for (int rt = 0; rt < NV; rt++) {
			if (root != -1 && rt != root) continue;
			if (depth[rt] != -1) continue;
			dfs_lowval(rt, 0, -1);
		}

		for (auto& v : adj) {
			v.clear();
		}
		for (const auto& bucket : lowval_buckets) {
			for (auto [e, cur, nxt] : bucket) {
				adj[cur].push_back({nxt, e});
			}
		}

		// Free this; literally no one cares about this
		lowval_buckets = {};
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

	using vestack_range_t = array_range_t<vestack_t, &spqr_tree::vestack>;

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
			int ve = int(vedges.size()) - 1;
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
		for (auto [nxt, e] : adj[cur]) {
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

		component_vertices.push_back(cur);
		vertices[cur].component = cur_component;

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
		for (auto [nxt, e] : adj[cur]) {
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
	}

	void build_spqr() {
		vertices.resize(NV);
		vertex_blocks.reserve(NV + NE);
		vertex_blocks_buf.reserve(NV + NE);

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
			component_vertices.push_back(rt);
			vertices[rt].component = cur_component;
			cur_component = -1;

			c.blocks.en = int(blocks.size());
			c.nodes.en = int(nodes.size());
			c.vedges.en = int(vedges.size());
			c.component_vertices.en = int(component_vertices.size());
		}

		vestack = {};
		estack = {};
		tstack = {};
		first_occurrence = {};
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
		// Leave depth since it's sometimes useful
		//depth = {};
	}
};
