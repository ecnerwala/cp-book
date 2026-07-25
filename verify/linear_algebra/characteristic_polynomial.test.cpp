// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/characteristic_polynomial

#include <bits/stdc++.h>
#include <cassert>

#include "char_poly.hpp"
#include "modnum.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	int N; std::cin >> N;
	using num = modnum<998244353>;
	std::vector<std::vector<num>> M(N, std::vector<num>(N));
	for (auto& r : M) for (auto& v : r) std::cin >> v;
	auto res = charPoly(M);
	for (int i = 0; i <= N; i++) {
		std::cout << res[i] << " \n"[i==N];
	}

	return 0;
}
