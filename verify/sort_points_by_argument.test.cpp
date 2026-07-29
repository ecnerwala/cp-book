// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/sort_points_by_argument

#include <bits/stdc++.h>
#include <cassert>

#include "geometry/point.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	using pt_t = Point<int, int64_t>;

	int N;
	std::cin >> N;
	std::vector<pt_t> P(N);
	for (auto& p : P) std::cin >> p;
	auto fix_0 = [&](pt_t p) -> pt_t { return p == pt_t(0, 0) ? pt_t(1, 0) : p; };
	auto cmp = angle_cmp_upto(pt_t(-1, 0));
	std::ranges::sort(P, [&](auto a, auto b) -> bool { return cmp(fix_0(a), fix_0(b)); });
	for (auto p : P) {
		std::cout << p.x << ' ' << p.y << '\n';
	}

	return 0;
}
