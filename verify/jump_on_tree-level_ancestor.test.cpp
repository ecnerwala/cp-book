// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/jump_on_tree

#include <bits/stdc++.h>
#include <cassert>

#include "level_ancestor.hpp"

int main() {
	std::ios_base::sync_with_stdio(false); std::cin.tie(nullptr);

	int N, Q; std::cin >> N >> Q;
	std::vector<std::vector<int>> adj(N);
	for (int e = 0; e < N-1; e++) {
		int u, v; std::cin >> u >> v;
		adj[u].push_back(v);
		adj[v].push_back(u);
	}
	std::vector<int> par(N);
	std::vector<int> depth(N);
	[&](this auto&& self, int cur, int prv, int d) -> void {
		par[cur] = prv;
		depth[cur] = d;
		for (int nxt : adj[cur]) {
			if (nxt == prv) continue;
			self(nxt, cur, d+1);
		}
	}(0, -1, 0);
	ecnerwala::level_ancestor la(par);


	for (int q = 0; q < Q; q++) {
		int s, t, i; std::cin >> s >> t >> i;
		std::cout << [&]() -> int {
			int l = la.lca(s, t);
			if (i <= depth[s] - depth[l]) {
				return la.get_ancestor(s, i);
			}
			i -= depth[s] - depth[l];
			if (i <= depth[t] - depth[l]) {
				return la.get_ancestor(t, depth[t] - depth[l] - i);
			}
			return -1;
		}() << '\n';
	}

	return 0;
}
