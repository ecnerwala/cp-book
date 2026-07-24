// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/convolution_mod

#include <bits/stdc++.h>

#include "fft.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	using num = modnum<998244353>;
	using E = ecnerwala::fft::fft_engine<num>;

	int N, M; std::cin >> N >> M;
	std::vector<num> A(N); for (auto& x : A) std::cin >> x;
	std::vector<num> B(M); for (auto& x : B) std::cin >> x;
	std::vector<num> C = ecnerwala::fft::multiply<E>(A, B);
	for (int i = 0; i < N+M-1; i++) {
		std::cout << C[i] << " \n"[i+1==N+M-1];
	}

	return 0;
}
