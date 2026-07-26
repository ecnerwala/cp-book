// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/multipoint_evaluation

#include <bits/stdc++.h>
#include <cassert>

#include "fft.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	using num = modnum<998244353>;
	using E = ecnerwala::fft::engines::ntt<num>;

	int N, M; std::cin >> N >> M;
	ecnerwala::poly::vec<E> F(N); for (auto& f : F) std::cin >> f;
	std::vector<num> P(M); for (auto& p : P) std::cin >> p;
	auto res = ecnerwala::poly::multipoint<E>(F, P);
	for (int i = 0; i < M; i++) {
		std::cout << res[i] << " \n"[i+1==M];
	}

	return 0;
}
