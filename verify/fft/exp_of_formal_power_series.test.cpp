// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/exp_of_formal_power_series

#include <bits/stdc++.h>
#include <cassert>

#include "num/modnum.hpp"
#include "fft/engines/ntt.hpp"
#include "fft/series.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	using num = wala::modnum<998244353>;
	using E = wala::fft::engines::ntt<num>;
	using ps = wala::series::trunc<E>;

	int N; std::cin >> N;
	ps A(N); for (auto& a : A) std::cin >> a;
	ps B = ps_exp(A);
	for (int i = 0; i < N; i++) {
		std::cout << B[i] << " \n"[i+1==N];
	}


	return 0;
}
