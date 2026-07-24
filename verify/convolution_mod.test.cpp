// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/convolution_mod

#include <bits/stdc++.h>

#include "fft.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	using num = modnum<998244353>;
	using E = ecnerwala::fft::fft_engine<num>;

	int N, M; std::cin >> N >> M;
	std::vector<num> a(N); for (auto& x : a) std::cin >> x;
	std::vector<num> b(M); for (auto& x : b) std::cin >> x;
	std::vector<num> c = ecnerwala::fft::multiply<E>(a, b);
	for (int i = 0; i < N+M-1; i++) {
		std::cout << c[i] << " \n"[i+1==N+M-1];
	}

	return 0;
}
