// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/sum_of_totient_function

#include <bits/stdc++.h>
#include <cassert>

#include "nt/dirichlet_series.hpp"
#include "num/modnum.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	int64_t N; std::cin >> N;
	static wala::dirichlet_series::div_vector_layout layout;
	layout = N;
	using num = wala::modnum<998244353>;
	using ds_prefix = wala::dirichlet_series::prefix<layout, num>;
	std::cout << (ds_prefix([&](int64_t x) { return num(x) * num(x+1) / num(2); }) / ds_prefix([&](int64_t x) { return num(x); }))[N] << '\n';

	return 0;
}
