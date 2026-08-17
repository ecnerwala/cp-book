// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/multipoint_evaluation

#include <bits/stdc++.h>
#include <cassert>

#include "num/modnum.hpp"
#include "fft/engines/ntt.hpp"
#include "fft/poly.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	using num = wala::modnum<998244353>;
	using E = wala::fft::engines::ntt<num>;

	int N, M; std::cin >> N >> M;
	wala::poly::vec<E> F(N); for (auto& f : F) std::cin >> f;
	std::vector<num> P(M); for (auto& p : P) std::cin >> p;
	auto res = wala::poly::multipoint<E>(F, P);
	for (int i = 0; i < M; i++) {
		std::cout << res[i] << " \n"[i+1==M];
	}

	return 0;
}
