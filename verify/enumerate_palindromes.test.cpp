// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/enumerate_palindromes

#include <bits/stdc++.h>
#include <cassert>

#include "manacher.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	std::string S;
	std::cin >> S;
	int N = int(S.size());
	auto res = manacher(S);
	for (int i = 1; i <= 2 * N - 1; i++) {
		std::cout << res[i] << " \n"[i == 2 * N - 1];
	}

	return 0;
}
