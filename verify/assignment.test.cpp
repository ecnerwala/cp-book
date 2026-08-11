// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/assignment

#include <bits/stdc++.h>
#include <cassert>

#include "mcmf.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	int N; std::cin >> N;
	// Costs may be negative; shift them all by OFS so they are nonnegative,
	// then subtract N * OFS from the total.
	const int64_t OFS = 1000000000;
	MCMF_SSPA<int, int64_t> m(2 * N + 2);
	int S = 2 * N, T = 2 * N + 1;
	std::vector<std::vector<int>> eid(N, std::vector<int>(N));
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			int64_t A; std::cin >> A;
			eid[i][j] = m.add_edge(i, N + j, 1, A + OFS);
		}
	}
	for (int i = 0; i < N; i++) m.add_edge(S, i, 1, 0);
	for (int j = 0; j < N; j++) m.add_edge(N + j, T, 1, 0);

	auto [flow, cost] = m.max_flow(S, T);
	assert(flow == N);
	std::cout << cost - int64_t(N) * OFS << '\n';
	std::vector<int> match(N, -1);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (m.flow(eid[i][j]) > 0) match[i] = j;
		}
	}
	for (int i = 0; i < N; i++) std::cout << match[i] << " \n"[i+1==N];

	return 0;
}
