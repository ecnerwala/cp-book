// competitive-verifier: IGNORE (the testcases are 5GB, too big for CI)
// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/convolution_mod_large

#include <bits/stdc++.h>
#include <cassert>

#include "fft/engines/crt.hpp"
#include "fft/series_core.hpp"
#include "modnum.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	using num = modnum<998244353>;
	using E = ecnerwala::fft::engines::crt<num>;
	using pse = ecnerwala::series::exact<E>;

	int N, M; std::cin >> N >> M;
	pse A(N); for (auto& x : A) std::cin >> x;
	pse B(M); for (auto& x : B) std::cin >> x;
	pse C = A * B;
	for (int i = 0; i < N+M-1; i++) {
		std::cout << C[i] << " \n"[i+1==N+M-1];
	}

	return 0;
}
