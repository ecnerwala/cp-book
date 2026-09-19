#include "graph/spqr_tree.hpp"

#include <random>
#include <algorithm>

#include <catch2/catch_test_macros.hpp>

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

			INFO("NV = " << NV);
			INFO("NE = " << NE);
			INFO("seed_seq = {" << NE << "," << seed << "}");

			auto UNSCOPED_INFO_graph = [&]() -> void {
				UNSCOPED_INFO("Graph: " << NE << " edges");
				for (int e = 0; e < NE; e++) {
					UNSCOPED_INFO("Edge " << e << ": " << edges[e][0] << "-" << edges[e][1]);
				}
			};

			// Use like this; it affects the next REQUIRE only
			//UNSCOPED_INFO_graph();

			using wala::spqr_tree;
			using node_type = spqr_tree::node_type;
			auto spqr = spqr_tree::build(NV, edges);

			// Basic bounds checks
			int num_items = int(spqr.par.size());

			REQUIRE(int(spqr.vert_index.size()) == NV);
			REQUIRE(int(spqr.edge_index.size()) == NE);
			REQUIRE(int(spqr.par.size()) == num_items);
			REQUIRE(int(spqr.subtree_end.size()) == num_items);
			REQUIRE(int(spqr.types.size()) == num_items);
			REQUIRE(int(spqr.orig_id.size()) == num_items);
			REQUIRE(int(spqr.ch.size()) == num_items);
			REQUIRE(int(spqr.node_verts.size()) == num_items);
			REQUIRE(int(spqr.vert_par_nv.size()) == num_items);
			REQUIRE(int(spqr.node_edges.size()) == num_items);
			REQUIRE(int(spqr.node_adj.size()) == 2 * int(spqr.node_verts.dat.size()));

			auto check_csr_bounds = [] <typename T> (wala::csr<T> c) -> void {
				REQUIRE(!c.bounds.empty());
				REQUIRE(c.bounds.front() == 0);
				REQUIRE(c.bounds.back() == int(c.dat.size()));
				for (int i = 0; i+1 < int(c.bounds.size()); i++) {
					REQUIRE(c.bounds[i] <= c.bounds[i+1]);
				}
			};
			check_csr_bounds(spqr.ch);
			check_csr_bounds(spqr.node_verts);
			check_csr_bounds(spqr.node_edges);
			check_csr_bounds(spqr.node_adj);

			// Check tree shape / preorder consistency
			REQUIRE(num_items >= 1);
			for (int i = 0; i < num_items; i++) {
				if (i > 0) {
					REQUIRE(spqr.par[i] >= 0);
					REQUIRE(spqr.par[i] < i);
				} else {
					REQUIRE(spqr.par[i] == -1);
				}
				int cur_end = i+1;
				for (int ch : spqr.ch[i]) {
					REQUIRE(ch == cur_end);
					REQUIRE(spqr.par[ch] == i);
					REQUIRE(spqr.subtree_end[ch] > ch);
					cur_end = spqr.subtree_end[ch];
				}
				REQUIRE(spqr.subtree_end[i] == cur_end);
			}
			REQUIRE(spqr.subtree_end[0] == num_items);

			// Check that all verts/edges are present exactly once
			for (int v = 0; v < NV; v++) {
				int i = spqr.vert_index[v];
				REQUIRE(0 <= i);
				REQUIRE(i < num_items);
				REQUIRE(spqr.types[i] == node_type::V);
				REQUIRE(spqr.orig_id[i] == v);
			}
			for (int e = 0; e < NE; e++) {
				int i = spqr.edge_index[e];
				REQUIRE(0 <= i);
				REQUIRE(i < num_items);
				REQUIRE(spqr.types[i] == node_type::Q);
				REQUIRE(spqr.orig_id[i] == e);
			}
			for (int i = 0; i < num_items; i++) {
				node_type i_type = spqr.types[i];

				if (i_type == node_type::V) {
					int v = spqr.orig_id[i];
					REQUIRE(0 <= v);
					REQUIRE(v < NV);
					REQUIRE(spqr.vert_index[v] == i);
				} else if (i_type == node_type::Q) {
					int e = spqr.orig_id[i];
					REQUIRE(0 <= e);
					REQUIRE(e < NE);
					REQUIRE(spqr.edge_index[e] == i);
				} else {
					REQUIRE(spqr.orig_id[i] == -1);
				}
			}

			// Now, we're guaranteed that edges/vertices are 1-to-1 with Q/V nodes.
			// Check the endpoints match the input
			for (int e = 0; e < NE; e++) {
				auto nvs = spqr.node_verts[spqr.edge_index[e]];
				std::array<int, 2> given_ends{spqr.vert_index[edges[e][0]], spqr.vert_index[edges[e][1]]};
				std::ranges::sort(given_ends);
				if (given_ends[0] == given_ends[1]) {
					REQUIRE(nvs.size() == 1);
					REQUIRE(given_ends[0] == nvs[0].vert);
				} else {
					std::array<int, 2> spqr_ends{nvs[0].vert, nvs[1].vert};
					REQUIRE(given_ends == spqr_ends);
				}
			}

			// Check node shapes/consistency
			for (int i = 0; i < num_items; i++) {
				node_type i_type = spqr.types[i];
				INFO("i = " << i);
				INFO("i_type = " << char(i_type));
				int p = spqr.par[i];
				INFO("p = " << p);
				node_type p_type = p == -1 ? node_type::F : spqr.types[p];
				INFO("p_type = " << char(p_type));
				auto ch = spqr.ch[i];
				auto nvs = spqr.node_verts[i];
				auto nes = spqr.node_edges[i];
				int nv_off = spqr.node_verts.bounds[i];

				for (const auto& nv : nvs) REQUIRE(nv.node == i);
				for (const auto& ne : nes) REQUIRE(ne.node == i);

				if (i_type != node_type::V) REQUIRE(spqr.vert_par_nv[i] == -1);
				else REQUIRE(spqr.vert_par_nv[i] >= 0);

				if (i == 0) {
					REQUIRE(p == -1);
					REQUIRE(i_type == node_type::F);
					REQUIRE(spqr.node_edges[i].empty());
					REQUIRE(int(ch.size()) == int(nvs.size()));
					for (int z = 0; z < int(ch.size()); z++) {
						REQUIRE(spqr.types[ch[z]] == node_type::V);
						REQUIRE(nvs[z].vert == ch[z]);
						REQUIRE(spqr.vert_par_nv[ch[z]] == nv_off + z);
						REQUIRE(spqr.node_adj[2 * (nv_off + z) + 0].empty());
						REQUIRE(spqr.node_adj[2 * (nv_off + z) + 1].empty());
					}
				} else {
					REQUIRE(p != -1);
					REQUIRE(i_type != node_type::F);

					if (i_type == node_type::V) {
						REQUIRE(nvs.empty());
						REQUIRE(nes.empty());

						for (int z = 0; z < int(ch.size()); z++) {
							REQUIRE(spqr.types[ch[z]] == node_type::Q);
						}
					} else if (i_type == node_type::Q && p_type == node_type::V) {
						REQUIRE(nes.size() == 1);
						REQUIRE(nvs[0].vert == p);
						REQUIRE(spqr.types[ch[0]] != node_type::V);
						REQUIRE(nes[0].twin_ne == spqr.node_edges.bounds[ch[0]]);

						if (edges[spqr.orig_id[i]][0] == edges[spqr.orig_id[i]][1]) {
							// Self-loop Q node
							REQUIRE(ch.size() == 1);
							REQUIRE(nvs.size() == 1);
							REQUIRE((nes[0].nvs == std::array<int, 2>{nv_off + 0, nv_off + 0}));
							REQUIRE(spqr.types[ch[0]] == node_type::O);
						} else {
							REQUIRE(ch.size() == 2);
							REQUIRE(nvs.size() == 2);
							REQUIRE(spqr.types[ch[1]] == node_type::V);
							REQUIRE(nvs[1].vert == ch[1]);
							REQUIRE((nes[0].nvs == std::array<int, 2>{nv_off + 0, nv_off + 1}));
						}
					} else if (i_type == node_type::Q || i_type == node_type::S || i_type == node_type::P || i_type == node_type::R || i_type == node_type::I || i_type == node_type::O) {
						REQUIRE((p_type == node_type::Q || p_type == node_type::S || p_type == node_type::P || p_type == node_type::R));
						REQUIRE(!(p_type == node_type::S && i_type == node_type::S));
						REQUIRE(!(p_type == node_type::P && i_type == node_type::P));

						REQUIRE(!nvs.empty());
						REQUIRE(!nes.empty());
						REQUIRE(nes[0].nvs == std::array<int, 2>{nv_off, nv_off + int(nvs.size()) - 1});

						int nxt_nv = 1, nxt_ne = 1;
						int last_loc = 0;
						for (auto j : ch) {
							int loc;
							if (spqr.types[j] == node_type::V) {
								REQUIRE(nvs[nxt_nv].vert == j);
								REQUIRE(spqr.vert_par_nv[j] == nv_off + nxt_nv);
								loc = 2 * (nv_off + nxt_nv);
								nxt_nv++;
							} else {
								REQUIRE(nes[nxt_ne].twin_ne == spqr.node_edges.bounds[j]);
								REQUIRE(spqr.node_verts.bounds[i] <= nes[nxt_ne].nvs[0]);
								REQUIRE(nes[nxt_ne].nvs[0] < nes[nxt_ne].nvs[1]);
								REQUIRE(nes[nxt_ne].nvs[1] < spqr.node_verts.bounds[i+1]);
								loc = nes[nxt_ne].nvs[0] + nes[nxt_ne].nvs[1];
								nxt_ne++;
							}
							REQUIRE(loc >= last_loc);
							last_loc = loc;
						}
						if (i_type != node_type::O) nxt_nv++;
						REQUIRE(nxt_nv == int(nvs.size()));
						REQUIRE(nxt_ne == int(nes.size()));

						if (i_type == node_type::O) {
							REQUIRE(p_type == node_type::V);
							REQUIRE(ch.empty());
						} else if (i_type == node_type::I) {
							REQUIRE(p_type == node_type::V);
							REQUIRE(ch.empty());
						} else if (i_type == node_type::Q) {
							REQUIRE(ch.empty());
						} else if (i_type == node_type::S) {
							REQUIRE(nes.size() == nvs.size());
							REQUIRE(nes.size() >= 3);
							for (int z = 0; z < int(nes.size()); z++) {
								REQUIRE(nes[z].nvs[0] == (z ? nv_off + z-1 : nv_off));
								REQUIRE(nes[z].nvs[1] == (z ? nv_off + z-0 : nv_off + int(nvs.size()) - 1));
							}
						} else if (i_type == node_type::P) {
							REQUIRE(nvs.size() == 2);
							REQUIRE(nes.size() >= 3);
							for (auto ne : nes) {
								REQUIRE(ne.nvs[0] == nv_off + 0);
								REQUIRE(ne.nvs[1] == nv_off + 1);
							}
						} else if (i_type == node_type::R) {
							REQUIRE(nvs.size() >= 4);
							REQUIRE(nes.size() >= 6);
							// TODO: What else should we check
						} else REQUIRE(false);
					} else REQUIRE(false);
				}

				for (int ne : spqr.node_edges.indices(i)) {
					// Check twins have matching vertices
					int twin_ne = spqr.node_edges.dat[ne].twin_ne;
					REQUIRE(spqr.node_edges.dat[twin_ne].twin_ne == ne);
					for (int z = 0; z < 2; z++) {
						REQUIRE(
							spqr.node_verts.dat[spqr.node_edges.dat[ne].nvs[z]].vert ==
							spqr.node_verts.dat[spqr.node_edges.dat[twin_ne].nvs[z]].vert
						);
					}
				}

				// Check node_adj
				// Check the total counts are correct
				REQUIRE(spqr.node_adj.bounds[2 * spqr.node_verts.bounds[i]] == 2 * spqr.node_edges.bounds[i]);
				REQUIRE(spqr.node_adj.bounds[2 * spqr.node_verts.bounds[i+1]] == 2 * spqr.node_edges.bounds[i+1]);
				for (int nv : spqr.node_verts.indices(i)) {
					for (auto [ne, dest] : spqr.node_adj[2 * nv + 0]) {
						REQUIRE(spqr.node_edges.dat[ne].nvs[1] == nv);
						REQUIRE(spqr.node_edges.dat[ne].nvs[0] == dest);
						REQUIRE(dest <= nv);
					}
					for (auto [ne, dest] : spqr.node_adj[2 * nv + 1]) {
						REQUIRE(spqr.node_edges.dat[ne].nvs[0] == nv);
						REQUIRE(spqr.node_edges.dat[ne].nvs[1] == dest);
						REQUIRE(dest >= nv);
					}
					if (i_type == node_type::P) {
						// Check that edge ids are strictly decreasing on the left, strictly increasing on the right
						auto adj0 = spqr.node_adj[2 * nv + 0];
						REQUIRE(std::ranges::adjacent_find(adj0, std::ranges::less_equal{}, &spqr_tree::node_adj_t::ne) == adj0.end());
						auto adj1 = spqr.node_adj[2 * nv + 1];
						REQUIRE(std::ranges::adjacent_find(adj1, std::ranges::greater_equal{}, &spqr_tree::node_adj_t::ne) == adj1.end());
					} else {
						for (int z = 0; z < 2; z++) {
							// Check that destinations are strictly decreasing
							auto adj = spqr.node_adj[2 * nv + z];
							REQUIRE(std::ranges::adjacent_find(adj, std::ranges::less_equal{}, &spqr_tree::node_adj_t::dest_nv) == adj.end());
						}
					}
				}
				// Because all adjacency lists are now guaranteed distinct by strict ordering, correct counts imply completeness.
			}
		}
	}
}
