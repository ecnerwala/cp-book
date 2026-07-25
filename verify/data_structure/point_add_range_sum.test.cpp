// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/point_add_range_sum

#include <bits/stdc++.h>

#include "seg_tree.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	int N, Q; std::cin >> N >> Q;
	seg_tree::in_order_layout layout(N);
	std::vector<int64_t> seg(2*N);
	for (int i = 0; i < N; i++) {
		int64_t x; std::cin >> x;
		seg[layout.get_point(i)] = x;
	}
	for (seg_tree::point a(N-1); a >= 1; a--) {
		seg[a] = seg[a.c(0)] + seg[a.c(1)];
	}
	for (int q = 0; q < Q; q++) {
		int op; std::cin >> op;
		if (op == 0) {
			int p; int64_t x; std::cin >> p >> x;
			layout.get_point(p).for_each_up([&](seg_tree::point a) {
				seg[a] += x;
			});
		} else if (op == 1) {
			int l, r; std::cin >> l >> r;
			int64_t ans = 0;
			layout.get_range(l, r).for_each([&](seg_tree::point a) {
				ans += seg[a];
			});
			std::cout << ans << '\n';
		} else assert(false);
	}

	return 0;
}
