// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/two_edge_connected_components

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

	std::vector<int> comp_verts(N, -1);
	std::vector<int> comp_bounds; comp_bounds.reserve(N+1); comp_bounds.push_back(0);

	for (int i = int(spqr.size()) - 1; i >= 0; i--) {
		int comp_end = -1;
		if (spqr.types[i] == node_type::V) {
			stk.push_back({i, spqr.orig_id[i]});
			if (spqr.par[i] == 0) {
				comp_end = spqr.subtree_end[i];
			}
		} else if (spqr.types[i] == node_type::I) {
			int p = spqr.par[i];
			assert(p == i-1 && spqr.types[p] == node_type::Q);
			comp_end = spqr.subtree_end[p];
		}
		if (comp_end != -1) {
			// Take off
			int idx = comp_bounds.back();
			while (!stk.empty() && stk.back().first < comp_end) {
				comp_verts[idx++] = stk.back().second;
				stk.pop_back();
			}
			comp_bounds.push_back(idx);
		}
	}
	assert(stk.empty());
	assert(comp_bounds.back() == N);

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
