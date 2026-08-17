// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/convolution_mod_1000000007

#include <bits/stdc++.h>
#include <cassert>

#include "num/wala::modnum.hpp"
#include "fft/engines/split.hpp"
#include "fft/series_core.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	using num = wala::modnum<int(1e9)+7>;
	using E = wala::fft::engines::split<num>;
	using pse = wala::series::exact<E>;

	int N, M; std::cin >> N >> M;
	pse A(N); for (auto& x : A) std::cin >> x;
	pse B(M); for (auto& x : B) std::cin >> x;
	pse C = A * B;
	for (int i = 0; i < N+M-1; i++) {
		std::cout << C[i] << " \n"[i+1==N+M-1];
	}

	return 0;
}
