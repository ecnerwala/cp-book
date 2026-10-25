// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/lca

#include <bits/stdc++.h>
#include <cassert>

#include "containers/vec.hpp"
#include "tree/static_tree.hpp"

int main() {
	std::ios_base::sync_with_stdio(false); std::cin.tie(nullptr);

	int N, Q; std::cin >> N >> Q;
	wala::vec<int> P(N);
	P[0] = -1;
	for (auto& v : P | std::views::drop(1)) std::cin >> v;

	wala::csr_index_builder adj_idx_builder(N);
	for (int i = 1; i < N; i++) {
		adj_idx_builder.count(P[i]);
		adj_idx_builder.count(i);
	}
	wala::csr_builder<int> adj_builder(std::move(adj_idx_builder));
	for (int i = 1; i < N; i++) {
		adj_builder.push(P[i]) = i;
		adj_builder.push(i) = P[i];
	}
	wala::csr<int> adj = std::move(adj_builder).finalize();

	wala::static_forest_t tree(adj, std::array{0});

	for (int q = 0; q < Q; q++) {
		int u, v; std::cin >> u >> v;
		std::cout << tree.lca(u, v) << '\n';
	}

	return 0;
}
