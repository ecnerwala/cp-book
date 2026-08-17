// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/cartesian_tree

#include <bits/stdc++.h>
#include <cassert>

#include "seq/cartesian_tree.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	int N; std::cin >> N;
	std::vector<int> A(N); for (auto& a : A) std::cin >> a;

	auto ct = CartesianTree::build_min_tree(A);
	for (int i = 0; i < N; i++) {
		int p = ct[2*i+1].p;
		if (p == -1) {
			p = i;
		} else {
			assert(p & 1);
			p /= 2;
		}
		std::cout << p << " \n"[i+1==N];
	}

	return 0;
}
