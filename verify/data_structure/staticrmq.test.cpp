// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/staticrmq

#include <bits/stdc++.h>

#include "rmq.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	int N, Q; std::cin >> N >> Q;
	std::vector<int> A(N); for (auto& x : A) std::cin >> x;
	RangeMinQuery<int> rmq(A);
	for (int q = 0; q < Q; q++) {
		int l, r; std::cin >> l >> r;
		std::cout << rmq.query(l, r-1) << '\n';
	}

	return 0;
}
