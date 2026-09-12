// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/polynomial_taylor_shift

#include <bits/stdc++.h>
#include <cassert>

#include "num/modnum.hpp"
#include "fft/engines/ntt.hpp"
#include "fft/poly.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	using num = wala::modnum<998244353>;
	using E = wala::fft::engines::ntt<num>;

	int N; num C; std::cin >> N >> C;
	wala::poly::vec<E> A(N); for (auto& a : A) std::cin >> a;
	wala::poly::vec<E> res = wala::poly::taylor_shift(A, C);
	for (int i = 0; i < N; i++) {
		std::cout << res[i] << " \n"[i+1==N];
	}

	return 0;
}
