#include "graph/spqr_tree.hpp"

#include <random>
#include <algorithm>
#include <numeric>
#include <tuple>

#include <catch2/catch_test_macros.hpp>

#define REQUIRE_FAST(...) do { if (!(__VA_ARGS__)) REQUIRE(__VA_ARGS__); } while (0)

TEST_CASE("SPQR Tree", "[spqr_tree]") {
	int NV = 10;
	for (int NE = 0; NE <= NV * NV; NE++) {
		for (int seed = 0; seed < 50; seed++) {
			std::seed_seq seq{NE, seed};
			std::mt19937 mt(seq);
			std::vector<std::array<int, 2>> edges(NE);
			for (auto& e : edges) {
				for (auto& v : e) {
					v = std::uniform_int_distribution<int>(0, NV-1)(mt);
				}
			}

			// order_mode 0: default order, 1: random full permutations, 2: random prefixes of permutations
			int order_mode = seed % 3;
			std::vector<int> vert_order(NV), edge_order(NE);
			std::ranges::iota(vert_order, 0);
			std::ranges::iota(edge_order, 0);
			if (order_mode == 0) {
				vert_order.clear();
				edge_order.clear();
			} else {
				std::ranges::shuffle(vert_order, mt);
				std::ranges::shuffle(edge_order, mt);
				if (order_mode == 2) {
					vert_order.resize(std::uniform_int_distribution<int>(0, NV)(mt));
					edge_order.resize(std::uniform_int_distribution<int>(0, NE)(mt));
				}
			}
			// The effective orders: listed ids first, then the rest in id order
			auto complete_order = [](std::vector<int> order, int n) -> std::vector<int> {
				std::vector<bool> listed(n);
				for (int i : order) listed[i] = true;
				for (int i = 0; i < n; i++) {
					if (!listed[i]) order.push_back(i);
				}
				return order;
			};
			std::vector<int> full_vert_order = complete_order(vert_order, NV);
			std::vector<int> full_edge_order = complete_order(edge_order, NE);

			CAPTURE(NV, NE, seed, edges, vert_order, edge_order);

			for (bool ternarize : {false, true}) {
				CAPTURE(ternarize);

				using wala::spqr_tree;
				using wala::planar_spqr_tree;
				using node_type = spqr_tree::node_type;
				auto spqr = planar_spqr_tree::build(NV, edges, ternarize, vert_order, edge_order);

				{
					// Building without planarity should give the same tree
					auto spqr_np = spqr_tree::build(NV, edges, ternarize, vert_order, edge_order);
					REQUIRE_FAST(spqr_np.vert_index == spqr.vert_index);
					REQUIRE_FAST(spqr_np.edge_index == spqr.edge_index);
					REQUIRE_FAST(spqr_np.par == spqr.par);
					REQUIRE_FAST(spqr_np.subtree_end == spqr.subtree_end);
					REQUIRE_FAST(spqr_np.types == spqr.types);
					REQUIRE_FAST(spqr_np.orig_id == spqr.orig_id);
					REQUIRE_FAST(spqr_np.ch.bounds == spqr.ch.bounds);
					REQUIRE_FAST(spqr_np.ch.dat == spqr.ch.dat);
					auto check_csr_equal = [] <typename T> (const wala::csr<T>& a, const wala::csr<T>& b, auto proj) -> void {
						REQUIRE_FAST(a.bounds == b.bounds);
						REQUIRE_FAST(std::ranges::equal(a.dat, b.dat, {}, proj, proj));
					};
					check_csr_equal(spqr_np.node_verts, spqr.node_verts, [](const spqr_tree::node_vert_t& x) { return std::tuple(x.node, x.vert); });
					REQUIRE_FAST(spqr_np.vert_par_nv == spqr.vert_par_nv);
					check_csr_equal(spqr_np.node_edges, spqr.node_edges, [](const spqr_tree::node_edge_t& x) { return std::tuple(x.node, x.twin_ne, x.nvs); });
					check_csr_equal(spqr_np.node_adj, spqr.node_adj, [](const spqr_tree::node_adj_t& x) { return std::tuple(x.ne, x.dest_nv); });
				}

				// Basic bounds checks
				int num_items = int(spqr.par.size());

				REQUIRE_FAST(int(spqr.vert_index.size()) == NV);
				REQUIRE_FAST(int(spqr.edge_index.size()) == NE);
				REQUIRE_FAST(int(spqr.par.size()) == num_items);
				REQUIRE_FAST(int(spqr.subtree_end.size()) == num_items);
				REQUIRE_FAST(int(spqr.types.size()) == num_items);
				REQUIRE_FAST(int(spqr.orig_id.size()) == num_items);
				REQUIRE_FAST(int(spqr.ch.size()) == num_items);
				REQUIRE_FAST(int(spqr.node_verts.size()) == num_items);
				REQUIRE_FAST(int(spqr.vert_par_nv.size()) == num_items);
				REQUIRE_FAST(int(spqr.node_edges.size()) == num_items);
				REQUIRE_FAST(int(spqr.node_adj.size()) == 2 * int(spqr.node_verts.dat.size()));

				auto check_csr_bounds = [] <typename T> (wala::csr<T> c) -> void {
					REQUIRE_FAST(!c.bounds.empty());
					REQUIRE_FAST(c.bounds.front() == 0);
					REQUIRE_FAST(c.bounds.back() == int(c.dat.size()));
					for (int i = 0; i+1 < int(c.bounds.size()); i++) {
						REQUIRE_FAST(c.bounds[i] <= c.bounds[i+1]);
					}
				};
				check_csr_bounds(spqr.ch);
				check_csr_bounds(spqr.node_verts);
				check_csr_bounds(spqr.node_edges);
				check_csr_bounds(spqr.node_adj);

				// Check tree shape / preorder consistency
				REQUIRE_FAST(num_items >= 1);
				for (int i = 0; i < num_items; i++) {
					if (i > 0) {
						REQUIRE_FAST(spqr.par[i] >= 0);
						REQUIRE_FAST(spqr.par[i] < i);
					} else {
						REQUIRE_FAST(spqr.par[i] == -1);
					}
					int cur_end = i+1;
					for (int ch : spqr.ch[i]) {
						REQUIRE_FAST(ch == cur_end);
						REQUIRE_FAST(spqr.par[ch] == i);
						REQUIRE_FAST(spqr.subtree_end[ch] > ch);
						cur_end = spqr.subtree_end[ch];
					}
					REQUIRE_FAST(spqr.subtree_end[i] == cur_end);
				}
				REQUIRE_FAST(spqr.subtree_end[0] == num_items);

				// Check that all verts/edges are present exactly once
				for (int v = 0; v < NV; v++) {
					int i = spqr.vert_index[v];
					REQUIRE_FAST(0 <= i);
					REQUIRE_FAST(i < num_items);
					REQUIRE_FAST(spqr.types[i] == node_type::V);
					REQUIRE_FAST(spqr.orig_id[i] == v);
				}
				for (int e = 0; e < NE; e++) {
					int i = spqr.edge_index[e];
					REQUIRE_FAST(0 <= i);
					REQUIRE_FAST(i < num_items);
					REQUIRE_FAST(spqr.types[i] == node_type::Q);
					REQUIRE_FAST(spqr.orig_id[i] == e);
				}
				for (int i = 0; i < num_items; i++) {
					node_type i_type = spqr.types[i];

					if (i_type == node_type::V) {
						int v = spqr.orig_id[i];
						REQUIRE_FAST(0 <= v);
						REQUIRE_FAST(v < NV);
						REQUIRE_FAST(spqr.vert_index[v] == i);
					} else if (i_type == node_type::Q) {
						int e = spqr.orig_id[i];
						REQUIRE_FAST(0 <= e);
						REQUIRE_FAST(e < NE);
						REQUIRE_FAST(spqr.edge_index[e] == i);
					} else {
						REQUIRE_FAST(spqr.orig_id[i] == -1);
					}
				}

				// Now, we're guaranteed that edges/vertices are 1-to-1 with Q/V nodes.
				// Check the endpoints match the input
				for (int e = 0; e < NE; e++) {
					auto nvs = spqr.node_verts[spqr.edge_index[e]];
					std::array<int, 2> given_ends{spqr.vert_index[edges[e][0]], spqr.vert_index[edges[e][1]]};
					std::ranges::sort(given_ends);
					if (given_ends[0] == given_ends[1]) {
						REQUIRE_FAST(nvs.size() == 1);
						REQUIRE_FAST(given_ends[0] == nvs[0].vert);
					} else {
						std::array<int, 2> spqr_ends{nvs[0].vert, nvs[1].vert};
						REQUIRE_FAST(given_ends == spqr_ends);
					}
				}

				// Check node shapes/consistency
				for (int i = 0; i < num_items; i++) {
					node_type i_type = spqr.types[i];
					CAPTURE(i);
					CAPTURE(i_type);
					int p = spqr.par[i];
					CAPTURE(p);
					node_type p_type = p == -1 ? node_type::F : spqr.types[p];
					CAPTURE(p_type);
					auto ch = spqr.ch[i];
					auto nvs = spqr.node_verts[i];
					auto nes = spqr.node_edges[i];
					int nv_off = spqr.node_verts.bounds[i];

					for (const auto& nv : nvs) REQUIRE_FAST(nv.node == i);
					for (const auto& ne : nes) REQUIRE_FAST(ne.node == i);

					if (i_type != node_type::V) REQUIRE_FAST(spqr.vert_par_nv[i] == -1);
					else REQUIRE_FAST(spqr.vert_par_nv[i] >= 0);

					if (i == 0) {
						REQUIRE_FAST(p == -1);
						REQUIRE_FAST(i_type == node_type::F);
						REQUIRE_FAST(spqr.node_edges[i].empty());
						REQUIRE_FAST(int(ch.size()) == int(nvs.size()));
						for (int z = 0; z < int(ch.size()); z++) {
							REQUIRE_FAST(spqr.types[ch[z]] == node_type::V);
							REQUIRE_FAST(nvs[z].vert == ch[z]);
							REQUIRE_FAST(spqr.vert_par_nv[ch[z]] == nv_off + z);
							REQUIRE_FAST(spqr.node_adj[2 * (nv_off + z) + 0].empty());
							REQUIRE_FAST(spqr.node_adj[2 * (nv_off + z) + 1].empty());
						}
						// Roots follow vert_order, and each root's first incident edge in edge_order roots its block
						{
							std::vector<bool> seen(NV);
							int z = 0;
							for (int r : full_vert_order) {
								if (seen[r]) continue;
								int rt = spqr.vert_index[r];
								REQUIRE_FAST(z < int(ch.size()));
								REQUIRE_FAST(ch[z] == rt);
								z++;
								for (int j = rt; j < spqr.subtree_end[rt]; j++) {
									if (spqr.types[j] == node_type::V) seen[spqr.orig_id[j]] = true;
								}
								auto it = std::ranges::find_if(full_edge_order, [&](int e) { return edges[e][0] == r || edges[e][1] == r; });
								if (it != full_edge_order.end()) {
									REQUIRE_FAST(spqr.par[spqr.edge_index[*it]] == rt);
								}
							}
							REQUIRE_FAST(z == int(ch.size()));
						}
					} else {
						REQUIRE_FAST(p != -1);
						REQUIRE_FAST(i_type != node_type::F);

						if (i_type == node_type::V) {
							REQUIRE_FAST(nvs.empty());
							REQUIRE_FAST(nes.empty());

							for (int z = 0; z < int(ch.size()); z++) {
								REQUIRE_FAST(spqr.types[ch[z]] == node_type::Q);
							}
						} else if (i_type == node_type::Q && p_type == node_type::V) {
							REQUIRE_FAST(nes.size() == 1);
							REQUIRE_FAST(nvs[0].vert == p);
							REQUIRE_FAST(spqr.types[ch[0]] != node_type::V);
							REQUIRE_FAST(nes[0].twin_ne == spqr.node_edges.bounds[ch[0]]);

							if (edges[spqr.orig_id[i]][0] == edges[spqr.orig_id[i]][1]) {
								// Self-loop Q node
								REQUIRE_FAST(ch.size() == 1);
								REQUIRE_FAST(nvs.size() == 1);
								REQUIRE_FAST((nes[0].nvs == std::array<int, 2>{nv_off + 0, nv_off + 0}));
								REQUIRE_FAST(spqr.types[ch[0]] == node_type::O);
							} else {
								REQUIRE_FAST(ch.size() == 2);
								REQUIRE_FAST(nvs.size() == 2);
								REQUIRE_FAST(spqr.types[ch[1]] == node_type::V);
								REQUIRE_FAST(nvs[1].vert == ch[1]);
								REQUIRE_FAST((nes[0].nvs == std::array<int, 2>{nv_off + 0, nv_off + 1}));
							}
						} else if (i_type == node_type::Q || i_type == node_type::S || i_type == node_type::P || i_type == node_type::R || i_type == node_type::I || i_type == node_type::O) {
							REQUIRE_FAST((p_type == node_type::Q || p_type == node_type::S || p_type == node_type::P || p_type == node_type::R));

							REQUIRE_FAST(!nvs.empty());
							REQUIRE_FAST(!nes.empty());
							REQUIRE_FAST(nes[0].nvs == std::array<int, 2>{nv_off, nv_off + int(nvs.size()) - 1});

							int nxt_nv = 1, nxt_ne = 1;
							int last_loc = 0;
							for (auto j : ch) {
								int loc;
								if (spqr.types[j] == node_type::V) {
									REQUIRE_FAST(nvs[nxt_nv].vert == j);
									REQUIRE_FAST(spqr.vert_par_nv[j] == nv_off + nxt_nv);
									loc = 2 * (nv_off + nxt_nv);
									nxt_nv++;
								} else {
									REQUIRE_FAST(nes[nxt_ne].twin_ne == spqr.node_edges.bounds[j]);
									REQUIRE_FAST(spqr.node_verts.bounds[i] <= nes[nxt_ne].nvs[0]);
									REQUIRE_FAST(nes[nxt_ne].nvs[0] < nes[nxt_ne].nvs[1]);
									REQUIRE_FAST(nes[nxt_ne].nvs[1] < spqr.node_verts.bounds[i+1]);
									loc = nes[nxt_ne].nvs[0] + nes[nxt_ne].nvs[1];
									nxt_ne++;
								}
								REQUIRE_FAST(loc >= last_loc);
								last_loc = loc;
							}
							if (i_type != node_type::O) nxt_nv++;
							REQUIRE_FAST(nxt_nv == int(nvs.size()));
							REQUIRE_FAST(nxt_ne == int(nes.size()));

							if (i_type == node_type::O) {
								REQUIRE_FAST(p_type == node_type::Q);
								REQUIRE_FAST(ch.empty());
							} else if (i_type == node_type::I) {
								REQUIRE_FAST(p_type == node_type::Q);
								REQUIRE_FAST(ch.empty());
							} else if (i_type == node_type::Q) {
								REQUIRE_FAST(ch.empty());
							} else if (i_type == node_type::S) {
								if (ternarize) {
									REQUIRE_FAST(nes.size() == 3);
								} else {
									REQUIRE_FAST(p_type != node_type::S);
								}
								REQUIRE_FAST(nes.size() == nvs.size());
								REQUIRE_FAST(nes.size() >= 3);
								for (int z = 0; z < int(nes.size()); z++) {
									REQUIRE_FAST(nes[z].nvs[0] == (z ? nv_off + z-1 : nv_off));
									REQUIRE_FAST(nes[z].nvs[1] == (z ? nv_off + z-0 : nv_off + int(nvs.size()) - 1));
								}
							} else if (i_type == node_type::P) {
								if (ternarize) {
									REQUIRE_FAST(nes.size() == 3);
								} else {
									REQUIRE_FAST(p_type != node_type::P);
								}
								REQUIRE_FAST(nvs.size() == 2);
								REQUIRE_FAST(nes.size() >= 3);
								for (auto ne : nes) {
									REQUIRE_FAST(ne.nvs[0] == nv_off + 0);
									REQUIRE_FAST(ne.nvs[1] == nv_off + 1);
								}
							} else if (i_type == node_type::R) {
								REQUIRE_FAST(nvs.size() >= 4);
								REQUIRE_FAST(nes.size() >= 6);
								// TODO: What else should we check
							} else REQUIRE_FAST(false);
						} else REQUIRE_FAST(false);
					}

					for (int ne : spqr.node_edges.indices(i)) {
						// Check twins have matching vertices
						int twin_ne = spqr.node_edges.dat[ne].twin_ne;
						REQUIRE_FAST(spqr.node_edges.dat[twin_ne].twin_ne == ne);
						for (int z = 0; z < 2; z++) {
							REQUIRE_FAST(
								spqr.node_verts.dat[spqr.node_edges.dat[ne].nvs[z]].vert ==
								spqr.node_verts.dat[spqr.node_edges.dat[twin_ne].nvs[z]].vert
							);
						}
					}

					// Check node_adj
					// Check the total counts are correct
					REQUIRE_FAST(spqr.node_adj.bounds[2 * spqr.node_verts.bounds[i]] == 2 * spqr.node_edges.bounds[i]);
					REQUIRE_FAST(spqr.node_adj.bounds[2 * spqr.node_verts.bounds[i+1]] == 2 * spqr.node_edges.bounds[i+1]);
					for (int nv : spqr.node_verts.indices(i)) {
						for (auto [ne, dest] : spqr.node_adj[2 * nv + 0]) {
							REQUIRE_FAST(spqr.node_edges.dat[ne].nvs[1] == nv);
							REQUIRE_FAST(spqr.node_edges.dat[ne].nvs[0] == dest);
							REQUIRE_FAST(dest <= nv);
						}
						for (auto [ne, dest] : spqr.node_adj[2 * nv + 1]) {
							REQUIRE_FAST(spqr.node_edges.dat[ne].nvs[0] == nv);
							REQUIRE_FAST(spqr.node_edges.dat[ne].nvs[1] == dest);
							REQUIRE_FAST(dest >= nv);
						}
						if (i_type == node_type::P) {
							// Check that edge ids are strictly decreasing on the left, strictly increasing on the right
							auto adj0 = spqr.node_adj[2 * nv + 0];
							REQUIRE_FAST(std::ranges::adjacent_find(adj0, std::ranges::less_equal{}, &spqr_tree::node_adj_t::ne) == adj0.end());
							auto adj1 = spqr.node_adj[2 * nv + 1];
							REQUIRE_FAST(std::ranges::adjacent_find(adj1, std::ranges::greater_equal{}, &spqr_tree::node_adj_t::ne) == adj1.end());
						} else {
							for (int z = 0; z < 2; z++) {
								// Check that destinations are strictly decreasing
								auto adj = spqr.node_adj[2 * nv + z];
								REQUIRE_FAST(std::ranges::adjacent_find(adj, std::ranges::less_equal{}, &spqr_tree::node_adj_t::dest_nv) == adj.end());
							}
						}
					}
					// Because all adjacency lists are now guaranteed distinct by strict ordering, correct counts imply completeness.
				}

				// Now, check planarity guarantees.
				std::vector<bool> face_vis(spqr.ne_rot_adj.size());
				std::vector<bool> vert_vis(spqr.ne_rot_adj.size());
				for (int i = 0; i < num_items; i++) {
					int rot_st = 4 * spqr.node_edges.bounds[i];
					int rot_en = 4 * spqr.node_edges.bounds[i+1];
					if (spqr.node_planar[i]) {
						for (int a = rot_st; a < rot_en; a++) {
							int b = spqr.ne_rot_adj[a];
							REQUIRE_FAST(b >= rot_st);
							REQUIRE_FAST(b < rot_en);
							// Make sure it's actually an involution/doubly-linked
							REQUIRE_FAST(spqr.ne_rot_adj[b] == a);
							REQUIRE_FAST((b & 1) != (a & 1));
							// Make sure the 2 endpoints have the same vertex
							REQUIRE_FAST(spqr.node_edges.dat[a >> 2].nvs[(a & 2) >> 1] == spqr.node_edges.dat[b >> 2].nvs[(b & 2) >> 1]);
						}

						if (spqr.types[i] == node_type::F || spqr.types[i] == node_type::V) {
							continue;
						}

						int num_verts = int(spqr.node_verts[i].size());
						int num_edges = int(spqr.node_edges[i].size());

						// Verify all edges around a vertex form a single cycle
						int num_vert_cycles = 0;
						for (int a = rot_st; a < rot_en; a++) {
							if (vert_vis[a]) continue;
							num_vert_cycles++;
							int cur = a;
							do {
								vert_vis[cur] = true;
								cur ^= 1;
								vert_vis[cur] = true;
								cur = spqr.ne_rot_adj[cur];
							} while (cur != a);
						}
						REQUIRE_FAST(num_vert_cycles == num_verts);

						// Verify the euler characteristic
						int num_face_cycles = 0;
						for (int a = rot_st; a < rot_en; a++) {
							if (face_vis[a]) continue;
							num_face_cycles++;
							int cur = a;
							do {
								face_vis[cur] = true;
								cur ^= 3;
								face_vis[cur] = true;
								cur = spqr.ne_rot_adj[cur];
							} while (cur != a);
						}
						REQUIRE_FAST(num_face_cycles == num_edges - num_verts + 2);
					} else {
						// Make sure everything's 0-ed out
						for (int b = rot_st; b < rot_en; b++) {
							REQUIRE_FAST(spqr.ne_rot_adj[b] == -1);
						}
					}
				}
			}
		}
	}
}
