// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/nim_product_64

#include <bits/stdc++.h>

#include "nim_prod.hpp"

constexpr nim_prod_t nimProd;

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	int T; std::cin >> T;
	while (T--) {
		uint64_t A, B; std::cin >> A >> B;
		std::cout << nimProd(A, B) << '\n';
	}

	return 0;
}
