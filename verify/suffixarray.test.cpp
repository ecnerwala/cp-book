// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/suffixarray

#include <bits/stdc++.h>

#include "suffix_array.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	std::string S; std::cin >> S;
	auto sa = SuffixArray::shift_and_construct(S);
	for (int i = 1; i <= int(S.size()); i++) {
		std::cout << sa.sa[i] << " \n"[i==int(S.size())];
	}

	return 0;
}
