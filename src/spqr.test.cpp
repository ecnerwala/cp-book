#include "spqr.hpp"

#include <map>
#include <random>
#include <algorithm>

#include <catch2/catch_test_macros.hpp>

TEST_CASE("SQPR Tree", "[spqr]") {
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
			UNSCOPED_INFO_graph();

			spqr_tree tree(NV, edges);

			// ==== CHECK SIZES ====
			REQUIRE(tree.NV == NV);
			REQUIRE(tree.NE == NE);
			REQUIRE(tree.NC == int(tree.components.size()));
			REQUIRE(tree.NB == int(tree.blocks.size()));
			REQUIRE(tree.NN == int(tree.nodes.size()));
			REQUIRE(tree.NVE == int(tree.vedges.size()));

			REQUIRE(int(tree.vertices.size()) == NV);
			REQUIRE(int(tree.vertex_blocks.size()) <= NV + NE);
			REQUIRE(int(tree.components.size()) <= NV);
			REQUIRE(int(tree.component_vertices.size()) == NV);
			REQUIRE(int(tree.blocks.size()) <= NE);
			REQUIRE(int(tree.block_vertices.size()) <= NV + NE);
			REQUIRE(int(tree.vertex_blocks.size()) == int(tree.block_vertices.size()));
			REQUIRE(int(tree.nodes.size()) <= 2 * NE);
			REQUIRE(int(tree.node_vertices.size()) <= NV + 5 * NE);
			REQUIRE(int(tree.vedges.size()) <= 4 * NE);

			// ==== CHECK VEDGES ====
			for (int ve = 0; ve < tree.NVE; ve++) {
				INFO("ve = " << ve);
				int node = tree.vedges[ve].node;
				int block = tree.vedges[ve].block;
				int component = tree.nodes[node].component;
				REQUIRE(node != -1);
				REQUIRE(block != -1);
				REQUIRE(component != -1);

				REQUIRE(tree.nodes[node].block == block);
				REQUIRE(tree.blocks[block].component == component);

				REQUIRE(tree.nodes[node].vedges.contains(ve));
				if (ve < tree.NE) {
					REQUIRE(!tree.blocks[block].vedges.contains(ve));
					REQUIRE(!tree.components[component].vedges.contains(ve));
				} else {
					REQUIRE(tree.blocks[block].vedges.contains(ve));
					REQUIRE(tree.components[component].vedges.contains(ve));
				}

				int o_ve = tree.vedges[ve].o_ve;
				REQUIRE(tree.vedges[ve].vs[0] == tree.vedges[o_ve].vs[1]);
				REQUIRE(tree.vedges[ve].vs[1] == tree.vedges[o_ve].vs[0]);
				REQUIRE(tree.vedges[ve].is_tree == !tree.vedges[o_ve].is_tree);
				REQUIRE(tree.vedges[ve].block == tree.vedges[o_ve].block);
				REQUIRE(tree.vedges[ve].node != tree.vedges[o_ve].node);

				REQUIRE(tree.vedges[o_ve].o_ve == ve);
				REQUIRE(tree.vedges[ve].o_node == tree.vedges[o_ve].node);
				REQUIRE(tree.vedges[ve].o_type == tree.nodes[tree.vedges[o_ve].node].type);
			}

			using node_type = spqr_tree::node_type;
			for (int e = 0; e < tree.NE; e++) {
				REQUIRE(
					std::min(edges[e][0], edges[e][1])
					== std::min(tree.vedges[e].vs[0], tree.vedges[e].vs[1])
				);
				REQUIRE(
					std::max(edges[e][0], edges[e][1])
					== std::max(tree.vedges[e].vs[0], tree.vedges[e].vs[1])
				);

				REQUIRE(tree.vedges[e].node == e);
				REQUIRE(tree.nodes[e].type == node_type::Q);
			}

			auto check_vertex_set = [&](auto vedges_range, auto verts_range) {
				auto vedges = vedges_range.bind(tree);
				std::vector<int> ve_verts; ve_verts.reserve(vedges.size()*2);
				for (const auto& vedge : vedges) {
					ve_verts.push_back(vedge.vs[0]);
					ve_verts.push_back(vedge.vs[1]);
				}
				std::sort(ve_verts.begin(), ve_verts.end());
				ve_verts.resize(std::unique(ve_verts.begin(), ve_verts.end()) - ve_verts.begin());

				auto verts = verts_range.bind(tree);
				std::vector<int> sorted_verts(verts.begin(), verts.end());
				std::sort(sorted_verts.begin(), sorted_verts.end());
				REQUIRE(ve_verts == sorted_verts);
			};

			// ==== CHECK NODES ====
			for (int node = 0; node < tree.NN; node++) {
				INFO("node = " << node);
				const auto& node_struct = tree.nodes[node];
				node_type type = node_struct.type;
				REQUIRE((type == node_type::Q) == (node < tree.NE));
				REQUIRE(node_struct.vedges.size() >= 1);

				REQUIRE(node_struct.block != -1);
				REQUIRE(node_struct.component != -1);
				REQUIRE(tree.blocks[node_struct.block].component == node_struct.component);
				if (node < tree.NE) {
					REQUIRE(!tree.blocks[node_struct.block].nodes.contains(node));
					REQUIRE(!tree.components[node_struct.component].nodes.contains(node));
				} else {
					REQUIRE(tree.blocks[node_struct.block].nodes.contains(node));
					REQUIRE(tree.components[node_struct.component].nodes.contains(node));
				}

				// Vertices are unique and are the correct set
				check_vertex_set(node_struct.vedges, node_struct.node_vertices);
				if (type == node_type::Q) {
					REQUIRE(node_struct.vedges.size() == 1);
				} else if (type == node_type::I) {
					REQUIRE(node_struct.vedges.size() == 1);
					REQUIRE(node_struct.node_vertices.size() == 2);
				} else if (type == node_type::O) {
					REQUIRE(node_struct.vedges.size() == node_struct.node_vertices.size());
					REQUIRE(node_struct.vedges.size() <= 2);
					int num_non_tree = 0;
					for (int z = 0; z < int(node_struct.vedges.size()); z++) {
						int ve = node_struct.vedges.st + z;
						auto vs = tree.vedges[ve].vs;
						REQUIRE(vs[0] == tree.node_vertices[z ? node_struct.node_vertices.st + z - 1 : node_struct.node_vertices.en - 1]);
						REQUIRE(vs[1] == tree.node_vertices[node_struct.node_vertices.st + z]);

						num_non_tree += !tree.vedges[ve].is_tree;
					}
					REQUIRE(num_non_tree == 1);
				} else if (type == node_type::S) {
					REQUIRE(node_struct.vedges.size() == node_struct.node_vertices.size());
					REQUIRE(node_struct.vedges.size() >= 3);
					int num_non_tree = 0;
					for (int z = 0; z < int(node_struct.vedges.size()); z++) {
						int ve = node_struct.vedges.st + z;
						auto vs = tree.vedges[ve].vs;
						REQUIRE(vs[0] == tree.node_vertices[z ? node_struct.node_vertices.st + z - 1 : node_struct.node_vertices.en - 1]);
						REQUIRE(vs[1] == tree.node_vertices[node_struct.node_vertices.st + z]);

						num_non_tree += !tree.vedges[ve].is_tree;
					}
					REQUIRE(num_non_tree == 1);
				} else if (type == node_type::P) {
					REQUIRE(node_struct.node_vertices.size() == 2);
					REQUIRE(node_struct.vedges.size() >= 3);
					int num_tree = 0;
					for (int z = 0; z < int(node_struct.vedges.size()); z++) {
						int ve = node_struct.vedges.st + z;
						auto vs = tree.vedges[ve].vs;
						REQUIRE(vs[0] != vs[1]);

						num_tree += tree.vedges[ve].is_tree;
					}
					REQUIRE(num_tree == 1);
				} else if (type == node_type::R) {
					REQUIRE(node_struct.node_vertices.size() >= 4);
					// Check for a few sanity things
					// No self-loop or duplicate edges (no trivial P)
					{
						std::vector<std::array<int, 2>> vedges(node_struct.vedges.size());
						for (int z = 0; z < int(node_struct.vedges.size()); z++) {
							int ve = node_struct.vedges.st + z;
							auto vs = tree.vedges[ve].vs;
							REQUIRE(vs[0] != vs[1]);
							if (vs[0] > vs[1]) std::swap(vs[0], vs[1]);
							vedges[z] = vs;
						}
						std::sort(vedges.begin(), vedges.end());
						REQUIRE(std::unique(vedges.begin(), vedges.end()) == vedges.end());
					}
					// All vertices have degree at least 3 (no trivial S)
					{
						std::map<int, int> degs;
						for (int z = 0; z < int(node_struct.vedges.size()); z++) {
							int ve = node_struct.vedges.st + z;
							auto vs = tree.vedges[ve].vs;
							for (int v : vs) degs[v]++;
						}
						for (auto [v, d] : degs) {
							REQUIRE(d > 2);
						}
					}
				} else REQUIRE(false);
			}

			for (int ve = 0; ve < tree.NVE; ve++) {
				int o_ve = tree.vedges[ve].o_ve;
				node_type t0 = tree.nodes[tree.vedges[ve].node].type;
				node_type t1 = tree.nodes[tree.vedges[o_ve].node].type;
				// 2 Q/S/P nodes cannot be glued together
				if (t0 == node_type::Q || t0 == node_type::S || t0 == node_type::P) {
					REQUIRE(t0 != t1);
				}
				// I and O can only be glued to Q
				if (t0 == node_type::I || t0 == node_type::O) {
					REQUIRE(t1 == node_type::Q);
				}
				if (t1 == node_type::I || t1 == node_type::O) {
					REQUIRE(t0 == node_type::Q);
				}
			}

			{
				int cur = 0;
				for (int node = 0; node < tree.NN; node++) {
					REQUIRE(tree.nodes[node].node_vertices.st == cur);
					int sz = tree.nodes[node].node_vertices.size();
					REQUIRE(sz > 0);
					cur += sz;
					REQUIRE(tree.nodes[node].node_vertices.en == cur);
					if (node < NE) {
						REQUIRE(sz <= 2);
						cur += (2 - sz);
					}
				}
				REQUIRE(cur == int(tree.node_vertices.size()));
			}
			{
				int cur = 0;
				for (int node = 0; node < tree.NN; node++) {
					REQUIRE(tree.nodes[node].vedges.st == cur);
					REQUIRE(tree.nodes[node].vedges.size() >= 1);
					cur += tree.nodes[node].vedges.size();
					REQUIRE(tree.nodes[node].vedges.en == cur);
				}
				REQUIRE(cur == int(tree.vedges.size()));
			}

			// ==== CHECK BLOCKS ====
			for (int block = 0; block < tree.NB; block++) {
				const auto& block_struct = tree.blocks[block];
				INFO("block = " << block);
				REQUIRE(block_struct.nodes.size() >= 1);
				REQUIRE(block_struct.component != -1);
				REQUIRE(tree.components[block_struct.component].blocks.contains(block));

				check_vertex_set(block_struct.vedges, block_struct.block_vertices);
			}

			// ==== CHECK COMPONENTS ====
			for (int component = 0; component < tree.NC; component++) {
				const auto& component_struct = tree.components[component];
				INFO("component = " << component);
				if (component_struct.blocks.size() > 0) {
					check_vertex_set(component_struct.vedges, component_struct.component_vertices);
				} else {
					REQUIRE(component_struct.component_vertices.size() == 1);
				}

				for (int v : component_struct.component_vertices.bind(tree)) {
					REQUIRE(tree.vertices[v].component == component);
				}
			}

			// ==== CHECK VERTEX_BLOCKS ===
			{
				// block / vertex pair
				std::vector<std::pair<int, int>> vb_pairs;
				for (int v = 0; v < tree.NV; v++) {
					for (int b : tree.vertices[v].vertex_blocks.bind(tree)) {
						vb_pairs.emplace_back(b, v);
					}
				}
				REQUIRE(vb_pairs.size() == tree.vertex_blocks.size());
				std::sort(vb_pairs.begin(), vb_pairs.end());

				std::vector<std::pair<int, int>> bv_pairs;
				for (int b = 0; b < tree.NB; b++) {
					for (int v : tree.blocks[b].block_vertices.bind(tree)) {
						bv_pairs.emplace_back(b, v);
					}
				}
				REQUIRE(bv_pairs.size() == tree.block_vertices.size());
				std::sort(bv_pairs.begin(), bv_pairs.end());

				REQUIRE(vb_pairs == bv_pairs);
			}

			// ==== CHECK PLANARITY / EMBEDDING ====
			{
				bool all_nodes_planar = true;
				std::vector<char> block_nodes_planar(tree.NB, 1);
				for (int node = 0; node < tree.NN; node++) {
					if (tree.nodes[node].type != node_type::R) {
						REQUIRE(tree.nodes[node].planar);
					}
					all_nodes_planar = all_nodes_planar && tree.nodes[node].planar;
					if (!tree.nodes[node].planar) {
						block_nodes_planar[tree.nodes[node].block] = 0;
					}
				}
				bool all_blocks_planar = true;
				for (int b = 0; b < tree.NB; b++) {
					REQUIRE(tree.blocks[b].planar == bool(block_nodes_planar[b]));
					all_blocks_planar = all_blocks_planar && tree.blocks[b].planar;
				}
				REQUIRE(tree.is_planar == all_nodes_planar);
				REQUIRE(tree.is_planar == all_blocks_planar);

				REQUIRE(int(tree.embed_next.size()) == 2 * tree.NVE);

				for (int node = 0; node < tree.NN; node++) {
					INFO("node = " << node);
					const auto& node_struct = tree.nodes[node];

					// The rotation at each skeleton vertex is a circular list
					// over exactly its incident half-edges within this node
					std::map<int, std::vector<int>> vertex_halves;
					for (int ve = node_struct.vedges.st; ve < node_struct.vedges.en; ve++) {
						vertex_halves[tree.vedges[ve].vs[0]].push_back(2*ve);
						vertex_halves[tree.vedges[ve].vs[1]].push_back(2*ve+1);
					}
					for (const auto& [v, halves] : vertex_halves) {
						INFO("v = " << v);
						std::vector<char> seen(halves.size(), 0);
						int h = halves[0];
						for (int cnt = 0; cnt < int(halves.size()); cnt++) {
							auto it = std::lower_bound(halves.begin(), halves.end(), h);
							REQUIRE(it != halves.end());
							REQUIRE(*it == h);
							REQUIRE(!seen[it - halves.begin()]);
							seen[it - halves.begin()] = 1;
							h = tree.embed_next[h];
						}
						REQUIRE(h == halves[0]);
					}

					// For planar nodes, the skeleton rotation must satisfy
					// Euler's formula: V - E + F == 2
					if (node_struct.planar) {
						int V = node_struct.node_vertices.size();
						int E = node_struct.vedges.size();
						int F = 0;
						std::map<int, char> face_seen;
						for (int ve = node_struct.vedges.st; ve < node_struct.vedges.en; ve++) {
							for (int k = 0; k < 2; k++) {
								int h0 = 2*ve + k;
								if (face_seen[h0]) continue;
								int h = h0;
								do {
									face_seen[h] = 1;
									h = tree.embed_next[h ^ 1];
								} while (h != h0);
								F++;
							}
						}
						REQUIRE(V - E + F == 2);
					}
				}
			}
		}
	}
}

TEST_CASE("SPQR Tree planarity fixtures", "[spqr]") {
	auto planarity = [](int NV, std::vector<std::array<int, 2>> edges) {
		return spqr_tree(NV, std::move(edges)).is_planar;
	};

	// K4
	REQUIRE(planarity(4, {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}}));
	// K5
	REQUIRE(!planarity(5, {{0,1},{0,2},{0,3},{0,4},{1,2},{1,3},{1,4},{2,3},{2,4},{3,4}}));
	// K5 minus an edge
	REQUIRE(planarity(5, {{0,1},{0,2},{0,3},{0,4},{1,2},{1,3},{1,4},{2,3},{2,4}}));
	// K3,3
	REQUIRE(!planarity(6, {{0,3},{0,4},{0,5},{1,3},{1,4},{1,5},{2,3},{2,4},{2,5}}));
	// Petersen graph
	{
		std::vector<std::array<int, 2>> edges;
		for (int i = 0; i < 5; i++) {
			edges.push_back({i, (i+1) % 5});
			edges.push_back({i, i+5});
			edges.push_back({5 + i, 5 + (i+2) % 5});
		}
		REQUIRE(!planarity(10, edges));
	}
	// Cube graph
	REQUIRE(planarity(8, {{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}}));
	// Multigraph with self loops, parallel edges, and a bridge
	REQUIRE(planarity(5, {{0,0},{0,1},{0,1},{0,1},{1,2},{2,3},{3,4},{4,2},{2,2}}));

	// K5 hanging off a square through a cutvertex: only the K5 block is non-planar
	{
		std::vector<std::array<int, 2>> edges = {{0,1},{1,2},{2,3},{3,0}};
		for (int i = 0; i < 5; i++) {
			for (int j = i+1; j < 5; j++) {
				edges.push_back({i == 0 ? 0 : i + 3, j + 3});
			}
		}
		spqr_tree tree(8, edges);
		REQUIRE(!tree.is_planar);
		int num_planar_blocks = 0;
		for (int b = 0; b < tree.NB; b++) {
			num_planar_blocks += tree.blocks[b].planar;
		}
		REQUIRE(tree.NB == 2);
		REQUIRE(num_planar_blocks == 1);
	}

	// K5 and K4 sharing the edge {0,1}: one block, two R nodes, and only the
	// K5's is non-planar
	{
		std::vector<std::array<int, 2>> edges;
		for (int i = 0; i < 5; i++) {
			for (int j = i+1; j < 5; j++) {
				edges.push_back({i, j});
			}
		}
		std::array<int, 4> k4 = {0, 1, 5, 6};
		for (int i = 0; i < 4; i++) {
			for (int j = i+1; j < 4; j++) {
				if (i == 0 && j == 1) continue;
				edges.push_back({k4[i], k4[j]});
			}
		}
		spqr_tree tree(7, edges);
		REQUIRE(!tree.is_planar);
		REQUIRE(tree.NB == 1);
		REQUIRE(!tree.blocks[0].planar);
		int num_R = 0, num_nonplanar = 0;
		for (int node = 0; node < tree.NN; node++) {
			num_R += tree.nodes[node].type == spqr_tree::node_type::R;
			num_nonplanar += !tree.nodes[node].planar;
		}
		REQUIRE(num_R == 2);
		REQUIRE(num_nonplanar == 1);
	}
}
