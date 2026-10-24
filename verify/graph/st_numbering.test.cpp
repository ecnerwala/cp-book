// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/st_numbering

#include <bits/stdc++.h>
#include <cassert>

#include "graph/make_st_dag.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	int N_CASES; std::cin >> N_CASES;
	while (N_CASES--) {
		int N, M, S, T; std::cin >> N >> M >> S >> T;
		std::vector<std::array<int, 2>> edges(M);
		wala::csr_index_builder adj_idx_builder(N);
		for (auto& [u, v] : edges) {
			std::cin >> u >> v;
			adj_idx_builder.count(u);
			adj_idx_builder.count(v);
		}
		wala::csr_builder<int> adj_builder(std::move(adj_idx_builder));
		for (auto [u, v] : edges) {
			adj_builder.push(u) = v;
			adj_builder.push(v) = u;
		}
		wala::csr<int> adj = std::move(adj_builder).finalize();

		auto res = wala::make_st_dag(adj, S, T);
		if (int(res.size()) == N) {
			std::cout << "Yes" << '\n';

			std::vector<int> pos(N);
			for (int i = 0; i < N; i++) {
				pos[res[i]] = i;
			}
			for (int i = 0; i < N; i++) {
				std::cout << pos[i] << " \n"[i+1==N];
			}
		} else {
			std::cout << "No" << '\n';
		}
	}

	return 0;
}
