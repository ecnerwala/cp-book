// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/suffixarray

#include <bits/stdc++.h>
#include <cassert>

#include "seq/suffix_array.hpp"
#include "fast_io.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	std::string S; std::cin >> S;
	auto sa = wala::SuffixArray::shift_and_construct(S);
	wala::FastWriter out;
	for (int i = 1; i <= int(S.size()); i++) {
		out << sa.sa[i] << " \n"[i==int(S.size())];
	}

	return 0;
}
