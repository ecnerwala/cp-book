// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/dirichlet_convolution_and_prefix_sums

#include <bits/stdc++.h>
#include <cassert>

#include "nt/dirichlet_series.hpp"
#include "num/modnum.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	int T; std::cin >> T;
	while (T--) {
		int64_t N; std::cin >> N;
		static wala::dirichlet_series::div_vector_layout layout;
		layout = N;
		using num = wala::modnum<998244353>;
		using ds_prefix = wala::dirichlet_series::prefix<layout, num>;
		ds_prefix F;
		for (int i = 1; i < layout.len; i++) std::cin >> F.st[i];
		ds_prefix G;
		for (int i = 1; i < layout.len; i++) std::cin >> G.st[i];
		ds_prefix H = F * G;
		for (int i = 1; i < layout.len; i++) std::cout << H.st[i] << " \n"[i+1==layout.len];
	}

	return 0;
}
