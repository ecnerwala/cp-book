// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/jump_on_tree

#include <bits/stdc++.h>
#include <cassert>

#include "containers/csr.hpp"
#include "containers/vec.hpp"
#include "tree/level_ancestor.hpp"

int main() {
	std::ios_base::sync_with_stdio(false); std::cin.tie(nullptr);

	int N, Q; std::cin >> N >> Q;
	wala::vec<std::array<int, 2>> edges(N-1);
	for (auto& [u, v] : edges) std::cin >> u >> v;
	wala::csr_index_builder adj_idx_builder(N);
	for (auto [u, v] : edges) {
		adj_idx_builder.count(u);
		adj_idx_builder.count(v);
	}
	wala::csr_builder<int> adj_builder(std::move(adj_idx_builder));
	for (auto [u, v] : edges) {
		adj_builder.push(u) = v;
		adj_builder.push(v) = u;
	}
	wala::csr<int> adj = std::move(adj_builder).finalize();
	wala::vec<int> par(N);
	wala::vec<int> depth(N);
	[&](this auto&& self, int cur, int prv, int d) -> void {
		par[cur] = prv;
		depth[cur] = d;
		for (int nxt : adj[cur]) {
			if (nxt == prv) continue;
			self(nxt, cur, d+1);
		}
	}(0, -1, 0);
	wala::level_ancestor la(par);

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
