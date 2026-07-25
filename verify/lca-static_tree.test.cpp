// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/lca

#include <bits/stdc++.h>

#include "static_tree.hpp"

int main() {
	std::ios_base::sync_with_stdio(false); std::cin.tie(nullptr);

	int N, Q; std::cin >> N >> Q;
	std::vector<int> P(N);
	P[0] = -1;
	for (auto& v : P | std::views::drop(1)) std::cin >> v;

	std::vector<std::vector<int>> adj(N);
	for (int i = 1; i < N; i++) {
		adj[P[i]].push_back(i);
		adj[i].push_back(P[i]);
	}

	static_forest_t tree(adj, {0});

	for (int q = 0; q < Q; q++) {
		int u, v; std::cin >> u >> v;
		std::cout << tree.lca(u, v) << '\n';
	}

	return 0;
}
