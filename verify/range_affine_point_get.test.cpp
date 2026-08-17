// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/range_affine_point_get

#include <bits/stdc++.h>
#include <cassert>

#include "ds/seg_tree.hpp"
#include "num/modnum.hpp"

int main() {
	std::ios_base::sync_with_stdio(false); std::cin.tie(nullptr);

	using num = modnum<998244353>;

	int N, Q; std::cin >> N >> Q;
	std::vector<num> A(N); for (auto& a : A) std::cin >> a;
	struct linear_fn {
		std::array<num, 2> v;
	};
	auto apply = [&](linear_fn a, linear_fn b) -> linear_fn {
		return {a.v[0] + a.v[1] * b.v[0], a.v[1] * b.v[1]};
	};
	wala::seg_tree::in_order_layout layout(N);
	std::vector<linear_fn> seg(2*layout.N);
	auto apply_lazy = [&](wala::seg_tree::point a, linear_fn f) -> void {
		seg[a] = apply(f, seg[a]);
	};
	auto downdate_node = [&](wala::seg_tree::point a) -> void {
		apply_lazy(a.c(0), seg[a]);
		apply_lazy(a.c(1), seg[a]);
		seg[a] = {num(0), num(1)};
	};
	for (int i = 0; i < N; i++) {
		auto a = layout.get_point(i);
		seg[a] = {A[i], 0};
	}
	for (wala::seg_tree::point a(N-1); a > 0; a--) {
		seg[a] = {0, 1};
	}

	for (int q = 0; q < Q; q++) {
		int op; std::cin >> op;
		if (op == 0) {
			int l, r; std::cin >> l >> r;
			linear_fn f; std::cin >> f.v[1] >> f.v[0];
			auto rng = layout.get_range(l, r);
			rng.for_parents_down(downdate_node);
			rng.for_each([&](wala::seg_tree::point a) -> void {
				apply_lazy(a, f);
			});
		} else if (op == 1) {
			int i; std::cin >> i;
			auto a = layout.get_point(i);
			a.for_parents_down(downdate_node);
			std::cout << seg[a].v[0] << '\n';
		} else assert(false);
	}

	return 0;
}
