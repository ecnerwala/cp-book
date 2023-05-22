#include "pq_tree.hpp"

#include "catch.hpp"
#include <bits/stdc++.h>

void check_tree(std::vector<int> A) {
	int N = int(A.size());
	std::vector<std::pair<std::array<int, 2>, std::array<int, 2>>> actual_ranges;
	for (int i = 0; i < N; i++) {
		int lo = A[i], hi = A[i];
		for (int j = i; j < N; j++) {
			lo = std::min(lo, A[j]);
			hi = std::max(hi, A[j]);
			assert(hi - lo >= j - i);
			if (hi - lo == j - i) {
				actual_ranges.push_back({{i, j}, {lo, hi}});
			}
		}
	}

	PQTree pq_tree(A);
	std::vector<std::pair<std::array<int, 2>, std::array<int, 2>>> computed_ranges;
	for (int n = 0; n < pq_tree.size(); n++) {
		const auto& node = pq_tree[n];
		if (node.type != PQTree::NodeType::PARTIAL) {
			computed_ranges.push_back({{node.l, node.r}, {node.lo, node.hi}});
		}
		if (node.type == PQTree::NodeType::LEAF) {
			REQUIRE(node.c[0] == -1);
			REQUIRE(node.c[1] == -1);
			REQUIRE(node.l == n/2);
			REQUIRE(node.r == n/2);
			REQUIRE(node.lo == A[n/2]);
			REQUIRE(node.hi == A[n/2]);
			continue;
		}
		REQUIRE(node.c[0] != -1);
		REQUIRE(node.c[1] != -1);
		REQUIRE(node.l == pq_tree[node.c[0]].l);
		REQUIRE(node.r == pq_tree[node.c[1]].r);
		REQUIRE(pq_tree[node.c[0]].r + 1 == pq_tree[node.c[1]].l);
		REQUIRE(node.lo == std::min(pq_tree[node.c[0]].lo, pq_tree[node.c[1]].lo));
		REQUIRE(node.hi == std::max(pq_tree[node.c[0]].hi, pq_tree[node.c[1]].hi));
		if (node.type == PQTree::NodeType::FULL) {
			// There should be at least 3 pieces
			REQUIRE((
				pq_tree[node.c[0]].type == PQTree::NodeType::PARTIAL
				|| pq_tree[node.c[1]].type == PQTree::NodeType::PARTIAL
			));
		}
		if (node.type == PQTree::NodeType::INCR) {
			REQUIRE(pq_tree[node.c[0]].hi + 1 == pq_tree[node.c[1]].lo);

			REQUIRE(pq_tree[node.c[1]].type != PQTree::NodeType::INCR);
			for (int cur = node.c[0]; pq_tree[cur].type == PQTree::NodeType::INCR; cur = pq_tree[cur].c[0]) {
				int ch = pq_tree[cur].c[1];
				computed_ranges.push_back({{pq_tree[ch].l, node.r}, {pq_tree[ch].lo, node.hi}});
			}
		}
		if (node.type == PQTree::NodeType::DECR) {
			REQUIRE(pq_tree[node.c[0]].lo - 1 == pq_tree[node.c[1]].hi);

			REQUIRE(pq_tree[node.c[1]].type != PQTree::NodeType::DECR);
			for (int cur = node.c[0]; pq_tree[cur].type == PQTree::NodeType::DECR; cur = pq_tree[cur].c[0]) {
				int ch = pq_tree[cur].c[1];
				computed_ranges.push_back({{pq_tree[ch].l, node.r}, {node.lo, pq_tree[ch].hi}});
			}
		}
	}
	std::sort(computed_ranges.begin(), computed_ranges.end());
	REQUIRE(actual_ranges == computed_ranges);
}

TEST_CASE("PQ Tree", "[pq_tree]") {
	for (int N = 1; N <= 7; N++) {
		std::vector<int> A(N);
		std::iota(A.begin(), A.end(), 0);
		do {
			check_tree(A);
		} while (next_permutation(A.begin(), A.end()));
	}
}
