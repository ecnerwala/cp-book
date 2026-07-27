// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/exp_of_formal_power_series

#include <bits/stdc++.h>
#include <cassert>

#include "fft/engines/ntt.hpp"
#include "fft/series.hpp"
#include "modnum.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	using num = modnum<998244353>;
	using E = ecnerwala::fft::engines::ntt<num>;
	using ps = ecnerwala::series::trunc<E>;

	int N; std::cin >> N;
	ps A(N); for (auto& a : A) std::cin >> a;
	ps B = ps_exp(A);
	for (int i = 0; i < N; i++) {
		std::cout << B[i] << " \n"[i+1==N];
	}


	return 0;
}
