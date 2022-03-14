#include "catch.hpp"

#include "spqr.hpp"

#include <map>

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
		}
	}
}
