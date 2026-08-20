// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/incremental_minimum_spanning_forest

#include <bits/stdc++.h>
#include <cassert>

#include "tree/top_tree.hpp"

namespace {

struct mst_top_tree_node : public wala::top_tree_node_base<mst_top_tree_node> {
	bool lazy_flip_path = false;
	int val = 0;
	int max_val = 0;

	void do_flip_path() {
		assert(is_path);
		std::swap(c[0], c[1]);
		lazy_flip_path ^= 1;
	}

	void downdate() {
		if (lazy_flip_path) {
			assert(is_path);
			if (!is_vert) {
				c[0]->do_flip_path();
				c[1]->do_flip_path();
			}
			lazy_flip_path = false;
		}
	}

	// NOTE: You may assume downdate() has been called on the current node, but
	// it may not have been called on the children! In particular, be careful
	// when accessing grandchildren information.
	void update() {
		if (is_vert) {
			max_val = 0;
		} else if (is_path) {
			max_val = std::max(val, std::max(c[0]->max_val, c[1]->max_val));
		} else {
			max_val = 0;
		}
	}
};

}

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	int N, M; std::cin >> N >> M;
	std::vector<mst_top_tree_node> verts(N);
	for (int i = 0; i < N; i++) {
		verts[i].is_vert = verts[i].is_path = true;
		verts[i].update();
	}
	// NB: Important to reserve() here so we pointers stay stable.
	std::vector<mst_top_tree_node> edges; edges.reserve(N-1);
	std::vector<int> edge_ids; edge_ids.reserve(N-1);
	std::vector<int> ans; ans.reserve(M);
	for (int q = 0; q < M; q++) {
		int u, v, w; std::cin >> u >> v >> w;
		if (u == v) {
			ans.push_back(q);
			continue;
		}
		auto a = &verts[u], b = &verts[v];
		a->make_root();
		assert(!a->p);
		b->expose();
		int e_idx;
		if (a->p) {
			// A moved, so we're looking for the edge to cut
			auto pth = b->p;
			if (w < pth->max_val) {
				// NB: We're playing with fire by not downdating, but it's okay for this particular tree
				while (pth->max_val != pth->val) {
					pth = pth->c[pth->c[1]->max_val > pth->c[0]->max_val];
				}
				e_idx = int(pth - edges.data());
				ans.push_back(edge_ids[e_idx]);
				cut(pth);
			} else {
				e_idx = -1;
				ans.push_back(q);
			}
		} else {
			e_idx = int(edges.size());
			ans.push_back(-1);
			edges.emplace_back();
			edge_ids.emplace_back();
		}
		if (e_idx != -1) {
			edges[e_idx].val = w;
			link(&edges[e_idx], a, b);
			edge_ids[e_idx] = q;
		}
	}
	for (int q = 0; q < M; q++) {
		std::cout << ans[q] << " \n"[q+1==M];
	}

	return 0;
}
