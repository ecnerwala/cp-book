// Non-planar dump of the reference implementation, used to differential-test the Lean port.
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

	auto t = wala::spqr_tree::build(NV, edges, ternarize, vert_order, edge_order);
	dump_vec("vert_index", t.vert_index);
	dump_vec("edge_index", t.edge_index);
	std::cout << "edge_flipped:"; for (bool x : t.edge_flipped) std::cout << ' ' << int(x); std::cout << '\n';
	dump_vec("par", t.par);
	dump_vec("subtree_end", t.subtree_end);
	std::cout << "types:"; for (auto x : t.types) std::cout << ' ' << x; std::cout << '\n';
	dump_vec("orig_id", t.orig_id);
	dump_vec("ch.bounds", t.ch.bounds);
	dump_vec("ch.dat", t.ch.dat);
	dump_vec("node_nvs.bounds", t.node_nvs.bounds);
	std::cout << "node_verts:"; for (auto x : t.node_verts) std::cout << ' ' << x.node << ',' << x.vert; std::cout << '\n';
	dump_vec("vert_par_nv", t.vert_par_nv);
	dump_vec("node_nes.bounds", t.node_nes.bounds);
	std::cout << "node_edges:"; for (auto x : t.node_edges) std::cout << ' ' << x.node << ',' << x.twin_ne << ',' << x.nvs[0] << ',' << x.nvs[1]; std::cout << '\n';
	dump_vec("node_adj.bounds", t.node_adj.bounds);
	std::cout << "node_adj.dat:"; for (auto x : t.node_adj.dat) std::cout << ' ' << x.ne << ',' << x.dest_nv; std::cout << '\n';
}
