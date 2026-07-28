// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/bipartitematching

#include <bits/stdc++.h>
#include <cassert>

#include "mcmf.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	int L, R, M; std::cin >> L >> R >> M;
	Dinic<int> d(L + R + 2);
	int S = L + R, T = L + R + 1;
	std::vector<int> eid(M);
	std::vector<std::pair<int, int>> es(M);
	for (int i = 0; i < M; i++) {
		int A, B; std::cin >> A >> B;
		es[i] = {A, B};
		eid[i] = d.add_edge(A, L + B, 1);
	}
	for (int a = 0; a < L; a++) d.add_edge(S, a, 1);
	for (int b = 0; b < R; b++) d.add_edge(L + b, T, 1);

	int f = d.max_flow(S, T);
	std::cout << f << '\n';
	for (int i = 0; i < M; i++) {
		if (d.flow(eid[i]) > 0) {
			std::cout << es[i].first << ' ' << es[i].second << '\n';
		}
	}

	return 0;
}
