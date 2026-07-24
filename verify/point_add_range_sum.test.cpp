// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/point_add_range_sum

#include <bits/stdc++.h>

#include "seg_tree.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	int N, Q; std::cin >> N >> Q;
	seg_tree::in_order_layout layout(N);
	std::vector<int64_t> tree(2*N);
	for (int i = 0; i < N; i++) {
		int64_t x; std::cin >> x;
		tree[layout.get_point(i)] = x;
	}
	for (int i = N-1; i >= 1; i--) {
		tree[i] = tree[2*i] + tree[2*i+1];
	}
	for (int q = 0; q < Q; q++) {
		int t; std::cin >> t;
		if (t == 0) {
			int p; int64_t x; std::cin >> p >> x;
			layout.get_point(p).for_each_up([&](seg_tree::point a) {
				tree[a] += x;
			});
		} else {
			int l, r; std::cin >> l >> r;
			int64_t ans = 0;
			layout.get_range(l, r).for_each([&](seg_tree::point a) {
				ans += tree[a];
			});
			std::cout << ans << '\n';
		}
	}

	return 0;
}
