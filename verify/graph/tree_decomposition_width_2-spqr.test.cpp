// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/tree_decomposition_width_2

#include <bits/stdc++.h>
#include <cassert>

#include "graph/spqr_tree.hpp"

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	std::string lit_p, lit_tw;
	int N, M; std::cin >> lit_p >> lit_tw >> N >> M;
	assert(lit_p == "p");
	assert(lit_tw == "tw");
	std::vector<std::array<int, 2>> edges(M);
	for (auto& [x, y] : edges) { std::cin >> x >> y; x--, y--; }

	auto spqr = wala::spqr_tree::build(N, edges, true);
	using node_type = wala::spqr_tree::node_type;
	int tree_width = 0;
	for (int i = 0; i < int(spqr.size()); i++) {
		auto t = spqr.types[i];
		if (t == node_type::F) {
			tree_width = std::max(tree_width, 0);
		} else if (t == node_type::V || t == node_type::O) {
			tree_width = std::max(tree_width, 0);
		} else if (t == node_type::Q) {
			// This depends if it's a self-loop or not
			tree_width = std::max(tree_width, int(spqr.node_verts[i].size()) - 1);
		} else if (t == node_type::I || t == node_type::P) {
			tree_width = std::max(tree_width, 1);
		} else if (t == node_type::S) {
			tree_width = std::max(tree_width, 2);
		} else if (t == node_type::R) {
			tree_width = std::max(tree_width, 3);
		} else assert(false);
	}

	if (tree_width >= 3) {
		std::cout << -1 << '\n';
	} else {
		std::cout << "s" << ' ' << "td" << ' ' << int(spqr.size()) << ' ' << tree_width << ' ' << N << '\n';
		for (int i = 0; i < int(spqr.size()); i++) {
			std::cout << "b" << ' ' << i+1;
			auto print_vert = [&](int v) -> void {
				std::cout << ' ' << spqr.orig_id[v] + 1;
			};
			auto t = spqr.types[i];
			if (t == node_type::F) {
				// Don't print anything
			} else if (t == node_type::V) {
				print_vert(i);
			} else {
				for (auto nv : spqr.node_verts[i]) {
					print_vert(nv.vert);
				}
			}
			std::cout << '\n';
		}
		for (int i = 1; i < int(spqr.size()); i++) {
			std::cout << i+1 << ' ' << spqr.par[i]+1 << '\n';
		}
	}
	return 0;
}
