// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/point_set_range_composite

#include <bits/stdc++.h>
#include <cassert>

#include "seg_tree.hpp"
#include "modnum.hpp"

int main() {
	std::ios_base::sync_with_stdio(false); std::cin.tie(nullptr);

	using num = modnum<998244353>;

	int N, Q; std::cin >> N >> Q;
	struct linear_fn {
		std::array<num, 2> v;
	};
	// compute b(a)
	auto compose = [&](linear_fn a, linear_fn b) -> linear_fn {
		return {b.v[0] + a.v[0] * b.v[1], a.v[1] * b.v[1]};
	};
	seg_tree::in_order_layout layout(N);
	std::vector<linear_fn> seg(2*N);
	auto update_node = [&](seg_tree::point a) -> void {
		seg[a] = compose(seg[a.c(0)], seg[a.c(1)]);
	};
	for (int i = 0; i < N; i++) {
		auto& v = seg[layout.get_point(i)];
		std::cin >> v.v[1] >> v.v[0];
	}
	for (seg_tree::point a(N-1); a > 0; a--) update_node(a);

	for (int q = 0; q < Q; q++) {
		int op; std::cin >> op;
		if (op == 0) {
			int p; std::cin >> p;
			seg_tree::point a = layout.get_point(p);
			std::cin >> seg[a].v[1] >> seg[a].v[0];
			a.for_parents_up(update_node);
		} else if (op == 1) {
			int l, r; std::cin >> l >> r;
			num x; std::cin >> x;
			linear_fn res{0, 1};
			layout.get_range(l, r).for_each_l_to_r([&](seg_tree::point a) -> void {
				res = compose(res, seg[a]);
			});
			std::cout << res.v[0] + res.v[1] * x << '\n';
		} else assert(false);
	}

	return 0;
}
