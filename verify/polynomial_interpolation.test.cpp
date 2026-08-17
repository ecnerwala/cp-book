// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/polynomial_interpolation

#include <bits/stdc++.h>
#include <cassert>

#include "num/modnum.hpp"
#include "fft/engines/ntt.hpp"
#include "fft/poly.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	using num = wala::modnum<998244353>;
	using E = wala::fft::engines::ntt<num>;

	int N; std::cin >> N;
	std::vector<num> X(N); for (auto& x : X) std::cin >> x;
	std::vector<num> Y(N); for (auto& y : Y) std::cin >> y;
	auto res = wala::poly::interpolate<E>(X, Y);
	for (int i = 0; i < N; i++) {
		std::cout << res[i] << " \n"[i+1==N];
	}

	return 0;
}
