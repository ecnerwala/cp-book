#include "catch.hpp"

#include "spqr.hpp"
#include "spqr_old.hpp"

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

			REQUIRE(tree.NV == NV);
			REQUIRE(tree.NE == NE);
			for (int ve = 0; ve < tree.NVE; ve++) {
				INFO("ve = " << ve);
				REQUIRE(tree.vedges[ve].block != -1);
				REQUIRE(tree.vedges[ve].node != -1);
				REQUIRE(tree.nodes[tree.vedges[ve].node].block == tree.vedges[ve].block);

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

			for (int node = 0; node < tree.NN; node++) {
				const auto& node_struct = tree.nodes[node];
				node_type type = node_struct.type;
				REQUIRE((type == node_type::Q) == (node < tree.NE));
				REQUIRE(node_struct.vedges.size() >= 1);

				// Vertices are unique and are the correct set
				{
					std::vector<int> ve_verts; ve_verts.reserve(node_struct.vedges.size()*2);
					for (int ve = node_struct.vedges.st; ve < node_struct.vedges.en; ve++) {
						ve_verts.push_back(tree.vedges[ve].vs[0]);
						ve_verts.push_back(tree.vedges[ve].vs[1]);
					}
					std::sort(ve_verts.begin(), ve_verts.end());
					ve_verts.resize(std::unique(ve_verts.begin(), ve_verts.end()) - ve_verts.begin());
					std::vector<int> sorted_verts(
						tree.node_vertices.begin() + node_struct.node_vertices.st,
						tree.node_vertices.begin() + node_struct.node_vertices.en
					);
					std::sort(sorted_verts.begin(), sorted_verts.end());
					REQUIRE(ve_verts == sorted_verts);
				}
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
		}
	}
}

TEST_CASE("SQPR Tree Old", "[spqr_old]") {
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

			spqr_tree_old tree(NV, edges);

			REQUIRE(tree.NV == NV);
			REQUIRE(tree.NE == NE);
			for (int ve = 0; ve < tree.NVE; ve += 2) {
				REQUIRE(tree.vedges[ve].v[0] == tree.vedges[ve+1].v[1]);
				REQUIRE(tree.vedges[ve].v[1] == tree.vedges[ve+1].v[0]);
				REQUIRE(tree.vedges[ve].is_tree == !tree.vedges[ve+1].is_tree);
				REQUIRE(tree.vedges[ve].block == tree.vedges[ve+1].block);
				REQUIRE(tree.vedges[ve].block != -1);
				REQUIRE(tree.vedges[ve+1].block != -1);
				REQUIRE(tree.vedges[ve].node != tree.vedges[ve+1].node);
				REQUIRE(tree.vedges[ve].node != -1);
				REQUIRE(tree.vedges[ve+1].node != -1);
				REQUIRE(tree.nodes[tree.vedges[ve].node].block == tree.vedges[ve].block);
				REQUIRE(tree.nodes[tree.vedges[ve+1].node].block == tree.vedges[ve+1].block);
			}
			using node_type = spqr_tree_old::node_type;
			for (int e = 0; e < tree.NE; e++) {
				REQUIRE(
					std::min(edges[e][0], edges[e][1])
					== std::min(tree.vedges[2*e].v[0], tree.vedges[2*e].v[1])
				);
				REQUIRE(
					std::max(edges[e][0], edges[e][1])
					== std::max(tree.vedges[2*e].v[0], tree.vedges[2*e].v[1])
				);

				REQUIRE(tree.vedges[2*e].node == e);
				REQUIRE(tree.nodes[e].type == node_type::Q);
			}

			for (int node = 0; node < tree.NN; node++) {
				const auto& node_struct = tree.nodes[node];
				node_type type = node_struct.type;
				REQUIRE((type == node_type::Q) == (node < tree.NE));
				REQUIRE(node_struct.vedges.size() >= 1);

				// Vertices are unique and are the correct set
				{
					std::vector<int> ve_verts; ve_verts.reserve(node_struct.vedges.size()*2);
					for (auto ve : node_struct.vedges) {
						ve_verts.push_back(tree.vedges[ve].v[0]);
						ve_verts.push_back(tree.vedges[ve].v[1]);
					}
					std::sort(ve_verts.begin(), ve_verts.end());
					ve_verts.resize(std::unique(ve_verts.begin(), ve_verts.end()) - ve_verts.begin());
					std::vector<int> sorted_verts = node_struct.vertices;
					std::sort(sorted_verts.begin(), sorted_verts.end());
					REQUIRE(ve_verts == sorted_verts);
				}
				if (type == node_type::Q) {
					REQUIRE(node_struct.vedges.size() == 1);
				} else if (type == node_type::I) {
					REQUIRE(node_struct.vedges.size() == 1);
					REQUIRE(node_struct.vertices.size() == 2);
				} else if (type == node_type::O) {
					REQUIRE(node_struct.vedges.size() == node_struct.vertices.size());
					REQUIRE(node_struct.vedges.size() <= 2);
					int num_non_tree = 0;
					for (int z = 0; z < int(node_struct.vedges.size()); z++) {
						int ve = node_struct.vedges[z];
						auto vs = tree.vedges[ve].v;
						REQUIRE(vs[0] == (z ? node_struct.vertices[z-1] : node_struct.vertices.back()));
						REQUIRE(vs[1] == node_struct.vertices[z]);

						num_non_tree += !tree.vedges[ve].is_tree;
					}
					REQUIRE(num_non_tree == 1);
				} else if (type == node_type::S) {
					REQUIRE(node_struct.vedges.size() == node_struct.vertices.size());
					REQUIRE(node_struct.vedges.size() >= 3);
					int num_non_tree = 0;
					for (int z = 0; z < int(node_struct.vedges.size()); z++) {
						int ve = node_struct.vedges[z];
						auto vs = tree.vedges[ve].v;
						REQUIRE(vs[0] == (z ? node_struct.vertices[z-1] : node_struct.vertices.back()));
						REQUIRE(vs[1] == node_struct.vertices[z]);

						num_non_tree += !tree.vedges[ve].is_tree;
					}
					REQUIRE(num_non_tree == 1);
				} else if (type == node_type::P) {
					REQUIRE(node_struct.vertices.size() == 2);
					REQUIRE(node_struct.vedges.size() >= 3);
					int num_tree = 0;
					for (int z = 0; z < int(node_struct.vedges.size()); z++) {
						int ve = node_struct.vedges[z];
						auto vs = tree.vedges[ve].v;
						REQUIRE(vs[0] != vs[1]);

						num_tree += tree.vedges[ve].is_tree;
					}
					REQUIRE(num_tree == 1);
				} else if (type == node_type::R) {
					REQUIRE(node_struct.vertices.size() >= 4);
					// Check for a few sanity things
					// No self-loop or duplicate edges (no trivial P)
					{
						std::vector<std::array<int, 2>> vedges(node_struct.vedges.size());
						for (int z = 0; z < int(node_struct.vedges.size()); z++) {
							int ve = node_struct.vedges[z];
							auto vs = tree.vedges[ve].v;
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
							int ve = node_struct.vedges[z];
							auto vs = tree.vedges[ve].v;
							for (int v : vs) degs[v]++;
						}
						for (auto [v, d] : degs) {
							REQUIRE(d > 2);
						}
					}
				} else REQUIRE(false);
			}
			for (int ve = 0; ve < tree.NVE; ve += 2) {
				node_type t0 = tree.nodes[tree.vedges[ve].node].type;
				node_type t1 = tree.nodes[tree.vedges[ve+1].node].type;
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
		}
	}
}
