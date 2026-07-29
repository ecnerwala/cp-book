// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/kth_term_of_linearly_recurrent_sequence

#include <bits/stdc++.h>
#include <cassert>

#include "modnum.hpp"
#include "fft/series.hpp"
#include "fft/engines/ntt.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	using num = modnum<998244353>;
	using E = ecnerwala::fft::engines::ntt<num>;
	using ps = ecnerwala::series::trunc<E>;
	using pse = ecnerwala::series::exact<E>;

	int D; int64_t K; std::cin >> D >> K;
	ps S(D); for (auto& v : S) std::cin >> v;
	pse Q(D+1);
	Q[0] = 1;
	for (auto& v : std::span<num>(Q).subspan(1)) {
		std::cin >> v; v = -v;
	}

	std::cout << kth_term_of_linear_recurrence(S, Q, K) << '\n';

	return 0;
}
