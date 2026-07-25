// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/sort_points_by_argument

#include <bits/stdc++.h>
#include <cassert>

#include "geometry/point.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	using pt_t = Point<int, int64_t>;

	int N; std::cin >> N;
	std::vector<pt_t> P(N);
	for (auto& p : P) std::cin >> p;
	std::ranges::sort(P, angle_cmp_upto(pt_t(-1, 0)));

	return 0;
}
