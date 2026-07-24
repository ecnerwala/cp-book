// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/common_interval_decomposition_tree

#include <bits/stdc++.h>

#include "perm_tree.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	int N; std::cin >> N;
	std::vector<int> P(N); for (auto& p : P) std::cin >> p;
	PermTree perm(P);
	int nxt_idx = 0;
	std::vector<std::tuple<int, std::array<int, 2>, bool>> res; res.reserve(2*N);
	[&](this auto&& self, int cur, int prv_idx, PermTree::NodeType par_type) -> void {
		if (perm[cur].type == PermTree::NodeType::PARTIAL || perm[cur].type == par_type) {
			self(perm[cur].c[0], prv_idx, par_type);
			self(perm[cur].c[1], prv_idx, par_type);
			return;
		}
		int cur_idx = nxt_idx++;
		res.push_back({prv_idx, {perm[cur].l, perm[cur].r}, perm[cur].type == PermTree::NodeType::FULL});
		if (perm[cur].type == PermTree::NodeType::LEAF) return;
		if (perm[cur].type == PermTree::NodeType::FULL) par_type = PermTree::NodeType::PARTIAL;
		else par_type = perm[cur].type;
		self(perm[cur].c[0], cur_idx, par_type);
		self(perm[cur].c[1], cur_idx, par_type);
	}(perm.root, -1, PermTree::NodeType::PARTIAL);
	std::cout << res.size() << '\n';
	for (auto [par, bounds, is_prime] : res) {
		std::cout << par << ' ' << bounds[0] << ' ' << bounds[1] << ' ' << (is_prime ? "prime" : "linear") << '\n';
	}

	return 0;
}
