// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/number_of_substrings

#include <bits/stdc++.h>
#include <cassert>

#include "suffix_array.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	std::string S; std::cin >> S;
	int N = int(S.size());
	auto sa = SuffixArrayLCP::shift_and_construct(S);
	long long ans = (long long)(N) * (N+1) / 2;
	for (int i = 0; i < N; i++) ans -= sa.lcp[i];
	std::cout << ans << '\n';

	return 0;
}
