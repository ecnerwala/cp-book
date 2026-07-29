// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/min_plus_convolution_convex_arbitrary

#include <bits/stdc++.h>
#include <cassert>

#include "smawk.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	int N, M;
	std::cin >> N >> M;
	std::vector<int> A(N);
	for (auto& a : A) std::cin >> a;
	std::vector<int> B(M);
	for (auto& b : B) std::cin >> b;
	const int INF = *std::max_element(A.begin(), A.end()) + *std::max_element(B.begin(), B.end()) + 1;
	auto res = smawk::smawk(N + M - 1, M, [&](int row, int col) -> int {
		int i = row - col;
		return i < 0 ? INF + M + (~i) : i >= N ? INF + (i-N) : (A[i] + B[col]); }, [&]([[maybe_unused]] int row, smawk::value_t<int> col1, smawk::value_t<int> col2) -> bool { return col1.v < col2.v ? 0 : 1; });
	for (int i = 0; i < N + M - 1; i++) {
		std::cout << res[i].v << " \n"[i + 1 == N + M - 1];
	}

	return 0;
}
