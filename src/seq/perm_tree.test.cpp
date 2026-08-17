#include "perm_tree.hpp"

#include <catch2/catch_test_macros.hpp>
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

	PermTree tree(A);
	std::vector<std::pair<std::array<int, 2>, std::array<int, 2>>> computed_ranges;
	for (int n = 0; n < tree.size(); n++) {
		const auto& node = tree[n];
		if (node.type != PermTree::NodeType::PARTIAL) {
			computed_ranges.push_back({{node.l, node.r}, {node.lo, node.hi}});
		}
		if (node.type == PermTree::NodeType::LEAF) {
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
		REQUIRE(node.l == tree[node.c[0]].l);
		REQUIRE(node.r == tree[node.c[1]].r);
		REQUIRE(tree[node.c[0]].r + 1 == tree[node.c[1]].l);
		REQUIRE(node.lo == std::min(tree[node.c[0]].lo, tree[node.c[1]].lo));
		REQUIRE(node.hi == std::max(tree[node.c[0]].hi, tree[node.c[1]].hi));
		if (node.type == PermTree::NodeType::FULL) {
			// There should be at least 3 pieces
			REQUIRE((
				tree[node.c[0]].type == PermTree::NodeType::PARTIAL
				|| tree[node.c[1]].type == PermTree::NodeType::PARTIAL
			));
		}
		if (node.type == PermTree::NodeType::INCR) {
			REQUIRE(tree[node.c[0]].hi + 1 == tree[node.c[1]].lo);

			REQUIRE(tree[node.c[1]].type != PermTree::NodeType::INCR);
			for (int cur = node.c[0]; tree[cur].type == PermTree::NodeType::INCR; cur = tree[cur].c[0]) {
				int ch = tree[cur].c[1];
				computed_ranges.push_back({{tree[ch].l, node.r}, {tree[ch].lo, node.hi}});
			}
		}
		if (node.type == PermTree::NodeType::DECR) {
			REQUIRE(tree[node.c[0]].lo - 1 == tree[node.c[1]].hi);

			REQUIRE(tree[node.c[1]].type != PermTree::NodeType::DECR);
			for (int cur = node.c[0]; tree[cur].type == PermTree::NodeType::DECR; cur = tree[cur].c[0]) {
				int ch = tree[cur].c[1];
				computed_ranges.push_back({{tree[ch].l, node.r}, {node.lo, tree[ch].hi}});
			}
		}
	}
	std::sort(computed_ranges.begin(), computed_ranges.end());
	REQUIRE(actual_ranges == computed_ranges);
}

TEST_CASE("Permutation Tree", "[perm_tree]") {
	for (int N = 1; N <= 7; N++) {
		std::vector<int> A(N);
		std::iota(A.begin(), A.end(), 0);
		do {
			check_tree(A);
		} while (next_permutation(A.begin(), A.end()));
	}
}
