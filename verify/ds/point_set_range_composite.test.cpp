// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/point_set_range_composite

#include <bits/stdc++.h>
#include <cassert>

#include "ds/seg_tree.hpp"
#include "num/linear_fn.hpp"
#include "num/modnum.hpp"

int main() {
	std::ios_base::sync_with_stdio(false); std::cin.tie(nullptr);

	using num = wala::modnum<998244353>;

	int N, Q; std::cin >> N >> Q;
	using linear_fn = wala::linear_fn<num>;
	wala::seg_tree::in_order_layout layout(N);
	std::vector<linear_fn> seg(2*N);
	auto update_node = [&](wala::seg_tree::point a) -> void {
		seg[a] = seg[a.c(1)] * seg[a.c(0)];
	};
	for (int i = 0; i < N; i++) {
		auto& v = seg[layout.get_point(i)];
		std::cin >> v.a >> v.b;
	}
	for (wala::seg_tree::point a(N-1); a > 0; a--) update_node(a);

	for (int q = 0; q < Q; q++) {
		int op; std::cin >> op;
		if (op == 0) {
			int p; std::cin >> p;
			wala::seg_tree::point a = layout.get_point(p);
			std::cin >> seg[a].a >> seg[a].b;
			a.for_parents_up(update_node);
		} else if (op == 1) {
			int l, r; std::cin >> l >> r;
			num x; std::cin >> x;
			layout.get_range(l, r).for_each_l_to_r([&](wala::seg_tree::point a) -> void {
				x = seg[a](x);
			});
			std::cout << x << '\n';
		} else assert(false);
	}

	return 0;
}
