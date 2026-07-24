// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/composition_of_formal_power_series_large

#include <bits/stdc++.h>

#include "fft.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	using num = modnum<998244353>;
	using E = ecnerwala::fft::fft_engine<num>;
	using ps = ecnerwala::power_series_trunc<E>;

	int N; std::cin >> N;
	ps A(N); for (auto& a : A) std::cin >> a;
	ps B(N); for (auto& b : B) std::cin >> b;
	auto res = poly_compose(A, B);
	for (int i = 0; i < N; i++) {
		std::cout << res[i] << " \n"[i+1==N];
	}

	return 0;
}
