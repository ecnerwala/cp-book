// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/range_affine_range_sum

#include <bits/stdc++.h>
#include <cassert>

#include "ds/seg_tree.hpp"
#include "num/linear_fn.hpp"
#include "num/modnum.hpp"

int main() {
	std::ios_base::sync_with_stdio(false); std::cin.tie(nullptr);

	using num = wala::modnum<998244353>;

	int N, Q; std::cin >> N >> Q;
	std::vector<num> A(N); for (auto& a : A) std::cin >> a;
	using linear_fn = wala::linear_fn<num>;
	struct seg_node {
		linear_fn lazy;
		num tot;
		num cnt;
	};
	wala::seg_tree::in_order_layout layout(N);
	std::vector<seg_node> seg(2*layout.N);
	auto apply_lazy = [&](wala::seg_tree::point a, linear_fn f) -> void {
		seg[a].lazy = f * seg[a].lazy;
		seg[a].tot = seg[a].tot * f.a + seg[a].cnt * f.b;
	};
	auto downdate_node = [&](wala::seg_tree::point a) -> void {
		apply_lazy(a.c(0), seg[a].lazy);
		apply_lazy(a.c(1), seg[a].lazy);
		seg[a].lazy = linear_fn{};
	};
	auto update_node = [&](wala::seg_tree::point a) -> void {
		seg[a].tot = seg[a.c(0)].tot + seg[a.c(1)].tot;
		seg[a].cnt = seg[a.c(0)].cnt + seg[a.c(1)].cnt;
	};
	for (int i = 0; i < N; i++) {
		auto a = layout.get_point(i);
		seg[a].lazy = linear_fn{};
		seg[a].tot = A[i];
		seg[a].cnt = 1;
	}
	for (wala::seg_tree::point a(N-1); a > 0; a--) {
		seg[a].lazy = linear_fn{};
		update_node(a);
	}

	for (int q = 0; q < Q; q++) {
		int op; std::cin >> op;
		if (op == 0) {
			int l, r; std::cin >> l >> r;
			linear_fn f; std::cin >> f.a >> f.b;
			auto rng = layout.get_range(l, r);
			rng.for_parents_down(downdate_node);
			rng.for_each([&](wala::seg_tree::point a) -> void {
				apply_lazy(a, f);
			});
			rng.for_parents_up(update_node);
		} else if (op == 1) {
			int l, r; std::cin >> l >> r;
			auto rng = layout.get_range(l, r);
			rng.for_parents_down(downdate_node);
			num ans = 0;
			rng.for_each([&](wala::seg_tree::point a) -> void {
				ans += seg[a].tot;
			});
			std::cout << ans << '\n';
		} else assert(false);
	}

	return 0;
}
