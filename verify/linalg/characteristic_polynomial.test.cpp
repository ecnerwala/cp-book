// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/characteristic_polynomial

#include <bits/stdc++.h>
#include <cassert>

#include "linalg/char_poly.hpp"
#include "num/modnum.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	int N; std::cin >> N;
	using num = wala::modnum<998244353>;
	wala::tensor<num, 2> M({N, N});
	for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) std::cin >> M[{i, j}];
	auto res = wala::charPoly(std::move(M));
	for (int i = 0; i <= N; i++) {
		std::cout << res[i] << " \n"[i==N];
	}

	return 0;
}
