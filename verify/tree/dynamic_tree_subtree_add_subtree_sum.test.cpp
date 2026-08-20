// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/dynamic_tree_subtree_add_subtree_sum

#include <bits/stdc++.h>
#include <cassert>

#include "tree/top_tree.hpp"

struct my_top_tree_node : public wala::top_tree_node_base<my_top_tree_node> {
	bool lazy_flip_path = false;
	int64_t lazy_add = 0;
	int64_t tot = 0;
	int64_t val = 0;
	int cnt = 0;

	void do_flip_path() {
		assert(is_path);
		std::swap(c[0], c[1]);
		lazy_flip_path ^= 1;
	}

	void do_lazy_add(int64_t add) {
		tot += cnt * add;
		val += add;
		lazy_add += add;
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
		if (lazy_add) {
			for (auto v : c) if (v) v->do_lazy_add(lazy_add);
			lazy_add = 0;
		}
	}

	// NOTE: You may assume downdate() has been called on the current node, but
	// it may not have been called on the children! In particular, be careful
	// when accessing grandchildren information.
	void update() {
		assert(!lazy_add);
		if (is_vert) {
			cnt = 1 + (c[0] ? c[0]->cnt : 0) + (c[1] ? c[1]->cnt : 0);
			tot = val + (c[0] ? c[0]->tot : 0) + (c[1] ? c[1]->tot : 0);
		} else if (is_path) {
			cnt = c[0]->cnt + c[1]->cnt;
			tot = c[0]->tot + c[1]->tot;
		} else {
			cnt = c[2]->cnt + (c[0] ? c[0]->cnt : 0) + (c[1] ? c[1]->cnt : 0);
			tot = c[2]->tot + (c[0] ? c[0]->tot : 0) + (c[1] ? c[1]->tot : 0);
		}
	}
};

int main() {
	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

	int N, Q; std::cin >> N >> Q;
	std::vector<my_top_tree_node> verts(N);
	for (int i = 0; i < N; i++) {
		std::cin >> verts[i].val;
		verts[i].is_vert = verts[i].is_path = true;
		verts[i].update();
	}
	std::vector<my_top_tree_node> edges(N-1);
	for (int e = 0; e < N-1; e++) {
		int u, v; std::cin >> u >> v;
		link(&edges[e], &verts[u], &verts[v]);
	}

	for (int q = 0; q < Q; q++) {
		int op; std::cin >> op;
		if (op == 0) {
			int u, v, w, x; std::cin >> u >> v >> w >> x;
			auto e = get_path(&verts[u], &verts[v]);
			cut(e);
			link(e, &verts[w], &verts[x]);
		} else if (op == 1) {
			int v, p, x; std::cin >> v >> p >> x;
			auto sub = get_subtree(&verts[p], &verts[v]);
			sub->do_lazy_add(x);
			sub->downdate();
			sub->update_all();
		} else if (op == 2) {
			int v, p; std::cin >> v >> p;
			auto sub = get_subtree(&verts[p], &verts[v]);
			std::cout << sub->tot << '\n';
		} else assert(false);
	}

	return 0;
}
