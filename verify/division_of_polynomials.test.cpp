// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/division_of_polynomials

#include <bits/stdc++.h>
#include <cassert>

#include "modnum.hpp"
#include "fft/engines/ntt.hpp"
#include "fft/poly.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	using num = modnum<998244353>;
	using E = ecnerwala::fft::engines::ntt<num>;

	int N, M; std::cin >> N >> M;
	ecnerwala::poly::vec<E> F(N); for (auto& f : F) std::cin >> f;
	ecnerwala::poly::vec<E> G(M); for (auto& g : G) std::cin >> g;

	auto [Q, R] = ecnerwala::poly::divmod(F, G);
	int rlen = R.len();
	while (rlen > 0 && R[rlen-1] == num(0)) rlen--;
	R.resize(rlen);

	std::cout << Q.len() << ' ' << R.len() << '\n';
	for (int i = 0; i < Q.len(); i++) std::cout << Q[i] << " \n"[i+1==Q.len()];
	if (Q.len() == 0) std::cout << '\n';
	for (int i = 0; i < R.len(); i++) std::cout << R[i] << " \n"[i+1==R.len()];
	if (R.len() == 0) std::cout << '\n';

	return 0;
}
