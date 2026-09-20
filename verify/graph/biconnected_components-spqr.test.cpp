// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/biconnected_components

#include <bits/stdc++.h>
#include <cassert>

#include "graph/spqr_tree.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	int N, M; std::cin >> N >> M;
	std::vector<std::array<int, 2>> edges(M);
	for (auto& [x, y] : edges) std::cin >> x >> y;

	auto spqr = wala::spqr_tree::build(N, edges);
	using node_type = wala::spqr_tree::node_type;
	std::vector<std::pair<int, int>> stk; stk.reserve(N);

	std::vector<int> comp_verts; comp_verts.reserve(2 * N);
	std::vector<int> comp_bounds; comp_bounds.reserve(N+1); comp_bounds.push_back(0);

	for (int i = int(spqr.size()) - 1; i >= 0; i--) {
		if (spqr.types[i] == node_type::V) {
			stk.push_back({i, spqr.orig_id[i]});
			if (spqr.par[i] == 0) {
				if (spqr.subtree_end[i] == i+1) {
					// Isolated vertex, print it as a special case
					comp_verts.push_back(stk.back().second);
					comp_bounds.push_back(int(comp_verts.size()));
				}
				stk.pop_back();
			}
		} else if (spqr.types[i] == node_type::Q && spqr.subtree_end[i] > i+1) {
			int p = spqr.par[i];
			assert(spqr.types[p] == node_type::V);
			assert(spqr.types[i+1] != node_type::O);

			comp_verts.push_back(spqr.orig_id[p]);
			int comp_end = spqr.subtree_end[i];
			while (!stk.empty() && stk.back().first < comp_end) {
				comp_verts.push_back(stk.back().second);
				stk.pop_back();
			}
			comp_bounds.push_back(int(comp_verts.size()));
		}
	}
	assert(stk.empty());

	int K = int(comp_bounds.size()) - 1;
	std::cout << K << '\n';
	for (int i = 0; i < K; i++) {
		std::cout << comp_bounds[i+1] - comp_bounds[i];
		for (int j = comp_bounds[i]; j < comp_bounds[i+1]; j++) {
			std::cout << ' ' << comp_verts[j];
		}
		std::cout << '\n';
	}

	return 0;
}
