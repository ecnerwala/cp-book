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
					auto check_vectors_equal = [] <typename T> (const std::vector<T>& a, const std::vector<T>& b, auto proj) -> void {
						REQUIRE_FAST(std::ranges::equal(a, b, {}, proj, proj));
					};
					check_vectors_equal(spqr_np.node_verts, spqr.node_verts, [](const spqr_tree::node_vert_t& x) { return std::tuple(x.node, x.vert); });
					REQUIRE_FAST(spqr_np.node_nvs.bounds == spqr.node_nvs.bounds);
					REQUIRE_FAST(spqr_np.vert_par_nv == spqr.vert_par_nv);
					check_vectors_equal(spqr_np.node_edges, spqr.node_edges, [](const spqr_tree::node_edge_t& x) { return std::tuple(x.node, x.twin_ne, x.nvs); });
					REQUIRE_FAST(spqr_np.node_nes.bounds == spqr.node_nes.bounds);
					REQUIRE_FAST(spqr_np.node_adj.bounds == spqr.node_adj.bounds);
					check_vectors_equal(spqr_np.node_adj.dat, spqr.node_adj.dat, [](const spqr_tree::node_adj_t& x) { return std::tuple(x.ne, x.dest_nv); });
				}

				// Basic bounds checks
				int num_items = int(spqr.par.size());

				REQUIRE_FAST(int(spqr.vert_index.size()) == NV);
				REQUIRE_FAST(int(spqr.edge_index.size()) == NE);
				REQUIRE_FAST(int(spqr.par.size()) == num_items);
				REQUIRE_FAST(int(spqr.subtree_end.size()) == num_items);
				REQUIRE_FAST(int(spqr.types.size()) == num_items);
				REQUIRE_FAST(int(spqr.orig_id.size()) == num_items);
				REQUIRE_FAST(spqr.ch.num_rows() == num_items);
				REQUIRE_FAST(spqr.node_nvs.num_rows() == num_items);
				REQUIRE_FAST(spqr.node_nvs.num_entries() == int(spqr.node_verts.size()));
				REQUIRE_FAST(int(spqr.vert_par_nv.size()) == num_items);
				REQUIRE_FAST(spqr.node_nes.num_rows() == num_items);
				REQUIRE_FAST(spqr.node_nes.num_entries() == int(spqr.node_edges.size()));
				REQUIRE_FAST(spqr.node_adj.num_rows() == 2 * int(spqr.node_verts.size()));

				auto check_csr_index = [] (const wala::csr_index& c) -> void {
					REQUIRE_FAST(!c.bounds.empty());
					REQUIRE_FAST(c.bounds.front() == 0);
					for (int i = 0; i+1 < int(c.bounds.size()); i++) {
						REQUIRE_FAST(c.bounds[i] <= c.bounds[i+1]);
					}
				};
				auto check_csr = [check_csr_index] <typename T> (const wala::csr<T>& c) -> void {
					check_csr_index(c);
					REQUIRE_FAST(c.num_entries() == int(c.dat.size()));
				};
				check_csr(spqr.ch);
				check_csr_index(spqr.node_nvs);
				check_csr_index(spqr.node_nes);
				check_csr(spqr.node_adj);

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
					auto nvs = spqr.node_nvs.slice(spqr.edge_index[e], spqr.node_verts);
					std::array<int, 2> given_ends{spqr.vert_index[edges[e][0]], spqr.vert_index[edges[e][1]]};
					if (given_ends[0] == given_ends[1]) {
						REQUIRE_FAST(nvs.size() == 1);
						REQUIRE_FAST(given_ends[0] == nvs[0].vert);
						REQUIRE_FAST(!spqr.edge_flipped[e]);
					} else {
						std::array<int, 2> spqr_ends{nvs[0].vert, nvs[1].vert};
						if (spqr.edge_flipped[e]) {
							std::swap(spqr_ends[0], spqr_ends[1]);
						}
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
					auto nvs = spqr.node_nvs.slice(i, spqr.node_verts);
					auto nes = spqr.node_nes.slice(i, spqr.node_edges);
					int nv_off = spqr.node_nvs.bounds[i];

					for (const auto& nv : nvs) REQUIRE_FAST(nv.node == i);
					for (const auto& ne : nes) REQUIRE_FAST(ne.node == i);

					if (i_type != node_type::V) REQUIRE_FAST(spqr.vert_par_nv[i] == -1);
					else REQUIRE_FAST(spqr.vert_par_nv[i] >= 0);

					if (i == 0) {
						REQUIRE_FAST(p == -1);
						REQUIRE_FAST(i_type == node_type::F);
						REQUIRE_FAST(nes.empty());
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
							REQUIRE_FAST(nes[0].twin_ne == spqr.node_nes.bounds[ch[0]]);

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
									REQUIRE_FAST(nes[nxt_ne].twin_ne == spqr.node_nes.bounds[j]);
									REQUIRE_FAST(spqr.node_nvs.bounds[i] <= nes[nxt_ne].nvs[0]);
									REQUIRE_FAST(nes[nxt_ne].nvs[0] < nes[nxt_ne].nvs[1]);
									REQUIRE_FAST(nes[nxt_ne].nvs[1] < spqr.node_nvs.bounds[i+1]);
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

					for (int ne : spqr.node_nes.indices(i)) {
						// Check twins have matching vertices
						int twin_ne = spqr.node_edges[ne].twin_ne;
						REQUIRE_FAST(spqr.node_edges[twin_ne].twin_ne == ne);
						for (int z = 0; z < 2; z++) {
							REQUIRE_FAST(
								spqr.node_verts[spqr.node_edges[ne].nvs[z]].vert ==
								spqr.node_verts[spqr.node_edges[twin_ne].nvs[z]].vert
							);
						}
					}

					// Check node_adj
					// Check the total counts are correct
					REQUIRE_FAST(spqr.node_adj.bounds[2 * spqr.node_nvs.bounds[i]] == 2 * spqr.node_nes.bounds[i]);
					REQUIRE_FAST(spqr.node_adj.bounds[2 * spqr.node_nvs.bounds[i+1]] == 2 * spqr.node_nes.bounds[i+1]);
					for (int nv : spqr.node_nvs.indices(i)) {
						for (auto [ne, dest] : spqr.node_adj[2 * nv + 0]) {
							REQUIRE_FAST(spqr.node_edges[ne].nvs[1] == nv);
							REQUIRE_FAST(spqr.node_edges[ne].nvs[0] == dest);
							REQUIRE_FAST(dest <= nv);
						}
						for (auto [ne, dest] : spqr.node_adj[2 * nv + 1]) {
							REQUIRE_FAST(spqr.node_edges[ne].nvs[0] == nv);
							REQUIRE_FAST(spqr.node_edges[ne].nvs[1] == dest);
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
				{
					auto check_planar_embedding = [](const wala::planar_embedding& pe, int V, const std::vector<std::array<int, 2>>& ends, const std::vector<bool>& is_embedded) {
						CAPTURE(pe.rot_adj);
						int E = int(ends.size());
						assert(int(is_embedded.size()) == E);
						REQUIRE_FAST(int(pe.rot_adj.size()) == 4 * E);

						// Check the involution basics
						for (int a = 0; a < 4 * E; a++) {
							if (!is_embedded[a >> 2]) {
								REQUIRE(pe.rot_adj[a] == -1);
								continue;
							}
							int b = pe.rot_adj[a];
							REQUIRE_FAST(b >= 0);
							REQUIRE_FAST(b < 4 * E);
							// Make sure it's actually an involution/doubly-linked
							REQUIRE_FAST(pe.rot_adj[b] == a);
							REQUIRE_FAST((b & 1) != (a & 1));
							// Make sure the 2 endpoints have the same vertex
							REQUIRE_FAST(ends[a >> 2][(a & 2) >> 1] == ends[b >> 2][(b & 2) >> 1]);
						}

						int expected_vert_cycles = 0;
						int expected_face_cycles = 0;
						{
							// Union find to identify components
							std::vector<bool> has_edge(V, false);
							std::vector<int> par(V, -1);
							auto get_par = [&](int a) -> int {
								while (par[a] >= 0) {
									if (par[par[a]] >= 0) par[a] = par[par[a]];
									a = par[a];
								}
								return a;
							};
							auto merge = [&](int a, int b) -> bool {
								a = get_par(a), b = get_par(b);
								if (a == b) return false;
								if (par[a] > par[b]) std::swap(a, b);
								par[a] += par[b];
								par[b] = a;
								return true;
							};

							for (int e = 0; e < E; e++) {
								if (!is_embedded[e]) continue;
								for (auto u : ends[e]) {
									assert(0 <= u && u < V);
									if (!has_edge[u]) {
										has_edge[u] = true;
										expected_vert_cycles++;
										expected_face_cycles++;
									}
								}
								// F = 2C + E - V
								expected_face_cycles += 1 - 2 * merge(ends[e][0], ends[e][1]);
							}
						}

						// Verify all edges around a vertex form a single cycle
						int num_vert_cycles = 0;
						{
							std::vector<bool> vert_vis(pe.rot_adj.size());
							for (int a = 0; a < 4 * E; a++) {
								if (!is_embedded[a >> 2]) continue;
								if (vert_vis[a]) continue;
								num_vert_cycles++;
								int cur = a;
								do {
									vert_vis[cur] = true;
									cur ^= 1;
									vert_vis[cur] = true;
									cur = pe.rot_adj[cur];
								} while (cur != a);
							}
						}
						// This must be true since all paired quarter-edges share a vertex, so no cycle jumps vertices.
						assert(num_vert_cycles >= expected_vert_cycles);
						REQUIRE_FAST(num_vert_cycles == expected_vert_cycles);

						// Verify the euler characteristic
						int num_face_cycles = 0;
						{
							std::vector<bool> face_vis(pe.rot_adj.size());
							for (int a = 0; a < 4 * E; a++) {
								if (!is_embedded[a >> 2]) continue;
								if (face_vis[a]) continue;
								num_face_cycles++;
								int cur = a;
								do {
									face_vis[cur] = true;
									cur ^= 3;
									face_vis[cur] = true;
									cur = pe.rot_adj[cur];
								} while (cur != a);
							}
						}
						// Worse embeddings can only have larger Euler characteristic
						assert(num_face_cycles >= expected_face_cycles);
						REQUIRE_FAST(num_face_cycles == expected_face_cycles);

						// Final thing: check to make sure that all parallel edges are grouped together correctly;
						// even though it's combinatorially valid, it would be impossible to form a straight-line drawing.
						std::vector<std::pair<std::array<int, 2>, int>> darts; darts.reserve(2 * E);
						for (int i = 0; i < E; i++) {
							for (int z = 0; z < 2; z++) {
								darts.push_back({{ends[i][z], ends[i][!z]}, 2 * i + z});
							}
						}
						std::sort(darts.begin(), darts.end());
						for (int i = 0, j = 0; i < int(darts.size()); i = j) {
							while (j < int(darts.size()) && darts[j].first == darts[i].first) j++;
							// Exclude self-loops
							if (darts[i].first[0] == darts[i].first[1]) continue;
							// Quick optimization
							if (j - i == 1) continue;
							CAPTURE(darts[i].first);
							int num_cuts = 0;
							for (int k = i; k < j; k++) {
								int a = 2 * darts[k].second + 1;
								// See if a is part of a 2-gon
								if (pe.rot_adj[a^3] != (pe.rot_adj[a]^3)) {
									num_cuts++;
								}
							}
							REQUIRE_FAST(num_cuts <= 1);
						}
					};
					{
						INFO("Checking partial embeddings");

						std::vector<std::array<int, 2>> ends(spqr.node_edges.size());
						std::vector<bool> is_embedded(spqr.node_edges.size());
						for (int i = 0; i < num_items; i++) {
							if (spqr.node_planar[i]) {
								for (int a : spqr.node_nes.indices(i)) {
									ends[a] = spqr.node_edges[a].nvs;
									is_embedded[a] = true;
								}
							} else {
								for (int a : spqr.node_nes.indices(i)) {
									ends[a] = {-1, -1};
									is_embedded[a] = false;
								}
							}
						}

						check_planar_embedding(spqr.ne_embedding, int(spqr.node_verts.size()), ends, is_embedded);
					}
					{
						INFO("Checking full embedding");
						auto full_embedding = planar_embed(spqr);
						REQUIRE_FAST(bool(full_embedding) == std::ranges::all_of(spqr.node_planar, std::identity{}));
						if (full_embedding) {
							std::vector<bool> is_embedded(edges.size(), true);
							check_planar_embedding(*full_embedding, NV, edges, is_embedded);
						}
					}
					{
						INFO("Checking direct full embedding");
						auto fast_embedding = wala::planar_embed(NV, edges, vert_order, edge_order);
						REQUIRE_FAST(bool(fast_embedding) == std::ranges::all_of(spqr.node_planar, std::identity{}));
						if (fast_embedding) {
							std::vector<bool> is_embedded(edges.size(), true);
							check_planar_embedding(*fast_embedding, NV, edges, is_embedded);
						}
					}
					{
						bool has_embedding = wala::can_planar_embed(NV, edges, vert_order, edge_order);
						REQUIRE(has_embedding == std::ranges::all_of(spqr.node_planar, std::identity{}));
					}
				}
			}
		}
	}
}
