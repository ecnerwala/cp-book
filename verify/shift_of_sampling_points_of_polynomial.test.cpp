// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/shift_of_sampling_points_of_polynomial

#include <bits/stdc++.h>
#include <cassert>

#include "num/modnum.hpp"
#include "fft/engines/ntt.hpp"
#include "fft/ap_sampled_poly.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	using num = wala::modnum<998244353>;
	using E = wala::fft::engines::ntt<num>;
	using ap_vals = wala::ap_sampled_poly<E>;

	int N, M; num C; std::cin >> N >> M >> C;
	ap_vals F(N);
	for (auto& v : F) std::cin >> v;
	ap_vals res = F.eval_range(C, M);
	for (int i = 0; i < M; i++) {
		std::cout << res[i] << " \n"[i+1==M];
	}

	return 0;
}
