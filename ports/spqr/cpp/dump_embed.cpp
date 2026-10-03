#include "graph/spqr_tree.hpp"
#include <iostream>

// dump.cpp plus the glued embedding of planar_embed ("-" if some node is nonplanar).
int main() {
	int NV, NE; int ternarize;
	std::cin >> NV >> NE >> ternarize;
	std::vector<std::array<int, 2>> edges(NE);
	for (auto& [u, v] : edges) std::cin >> u >> v;
	int K; std::cin >> K; std::vector<int> vert_order(K); for (auto& x : vert_order) std::cin >> x;
	int L; std::cin >> L; std::vector<int> edge_order(L); for (auto& x : edge_order) std::cin >> x;

	auto t = wala::planar_spqr_tree::build(NV, edges, ternarize, vert_order, edge_order);
	auto emb = wala::planar_embed(t);
	std::cout << "planar_embed:";
	if (!emb) std::cout << " -";
	else for (int x : emb->rot_adj) std::cout << ' ' << x;
	std::cout << '\n';
}
