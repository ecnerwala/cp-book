// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/st_numbering

#include <bits/stdc++.h>
#include <cassert>

#include "graph/make_st_dag.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	int N_CASES; std::cin >> N_CASES;
	while (N_CASES--) {
		int N, M, S, T; std::cin >> N >> M >> S >> T;
		std::vector<std::vector<int>> adj(N);
		for (int e = 0; e < M; e++) {
			int u, v; std::cin >> u >> v;
			adj[u].push_back(v);
			adj[v].push_back(u);
		}

		auto res = make_st_dag(adj, S, T);
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
