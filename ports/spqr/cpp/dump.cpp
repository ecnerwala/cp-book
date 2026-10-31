#include "graph/spqr_tree.hpp"
#include <iostream>

template <typename T> void dump_vec(const char* name, const std::vector<T>& v) {
	std::cout << name << ':';
	for (auto x : v) std::cout << ' ' << x;
	std::cout << '\n';
}

int main() {
	int NV, NE; int ternarize;
	std::cin >> NV >> NE >> ternarize;
	std::vector<std::array<int, 2>> edges(NE);
	for (auto& [u, v] : edges) std::cin >> u >> v;
	int K; std::cin >> K; std::vector<int> vert_order(K); for (auto& x : vert_order) std::cin >> x;
	int L; std::cin >> L; std::vector<int> edge_order(L); for (auto& x : edge_order) std::cin >> x;

	auto t = wala::planar_spqr_tree::build(NV, edges, ternarize, vert_order, edge_order);
	dump_vec("vert_index", t.vert_index);
	dump_vec("edge_index", t.edge_index);
	std::cout << "edge_flipped:"; for (bool x : t.edge_flipped) std::cout << ' ' << int(x); std::cout << '\n';
	dump_vec("par", t.par);
	dump_vec("subtree_end", t.subtree_end);
	std::cout << "types:"; for (auto x : t.types) std::cout << ' ' << x; std::cout << '\n';
	dump_vec("orig_id", t.orig_id);
	dump_vec("ch.bounds", t.ch.bounds);
	dump_vec("ch.dat", t.ch.dat);
	dump_vec("node_verts.bounds", t.node_nvs.bounds);
	std::cout << "node_verts.dat:"; for (auto x : t.node_verts) std::cout << ' ' << x.node << ',' << x.vert; std::cout << '\n';
	dump_vec("vert_par_nv", t.vert_par_nv);
	dump_vec("node_edges.bounds", t.node_nes.bounds);
	std::cout << "node_edges.dat:"; for (auto x : t.node_edges) std::cout << ' ' << x.node << ',' << x.twin_ne << ',' << x.nvs[0] << ',' << x.nvs[1]; std::cout << '\n';
	dump_vec("node_adj.bounds", t.node_adj.bounds);
	std::cout << "node_adj.dat:"; for (auto x : t.node_adj.dat) std::cout << ' ' << x.ne << ',' << x.dest_nv; std::cout << '\n';
	std::cout << "node_planar:"; for (bool x : t.node_planar) std::cout << ' ' << int(x); std::cout << '\n';
	dump_vec("ne_rot_adj", t.ne_embedding.rot_adj);

	// Also the non-planar build must agree on the shared fields
	auto s = wala::spqr_tree::build(NV, edges, ternarize, vert_order, edge_order);
	bool same = s.vert_index == t.vert_index && s.edge_index == t.edge_index && s.edge_flipped == t.edge_flipped && s.par == t.par && s.subtree_end == t.subtree_end && s.types == t.types && s.orig_id == t.orig_id && s.ch.bounds == t.ch.bounds && s.ch.dat == t.ch.dat && s.node_nvs.bounds == t.node_nvs.bounds && s.vert_par_nv == t.vert_par_nv && s.node_nes.bounds == t.node_nes.bounds && s.node_adj.bounds == t.node_adj.bounds;
	for (size_t i = 0; i < s.node_verts.size(); i++) same &= s.node_verts[i].node == t.node_verts[i].node && s.node_verts[i].vert == t.node_verts[i].vert;
	for (size_t i = 0; i < s.node_edges.size(); i++) same &= s.node_edges[i].node == t.node_edges[i].node && s.node_edges[i].twin_ne == t.node_edges[i].twin_ne && s.node_edges[i].nvs == t.node_edges[i].nvs;
	for (size_t i = 0; i < s.node_adj.dat.size(); i++) same &= s.node_adj.dat[i].ne == t.node_adj.dat[i].ne && s.node_adj.dat[i].dest_nv == t.node_adj.dat[i].dest_nv;
	std::cout << "nonplanar_build_same: " << int(same) << '\n';
}
