// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/polynomial_interpolation

#include <bits/stdc++.h>

#include "fft.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	using num = modnum<998244353>;
	using E = ecnerwala::fft::fft_engine<num>;

	int N; std::cin >> N;
	std::vector<num> X(N); for (auto& x : X) std::cin >> x;
	std::vector<num> Y(N); for (auto& y : Y) std::cin >> y;
	auto res = ecnerwala::poly_interpolate<E>(X, Y);
	for (int i = 0; i < N; i++) {
		std::cout << res[i] << " \n"[i+1==N];
	}

	return 0;
}
