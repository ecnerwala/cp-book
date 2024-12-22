#pragma once

#include "yc.hpp"
#include "rmq.hpp"

struct static_forest_t {
	int N;
	// original label to preorder
	std::vector<int> idx;

	// all keys/values are by preorder relabelling
	std::vector<int> preorder;
	std::vector<int> depth;
	std::vector<int> par;
	std::vector<int> sz;

	std::vector<int> depth_val_to_idx;
	RangeMinQuery<int> depth_val_rmq;

	static_forest_t() : N(0) {}
	static_forest_t(std::vector<std::vector<int>> adj, std::vector<int> roots = {}) :
		N(int(adj.size())),
		idx(N, -1),
		preorder(N, -1),
		depth(N, -1),
		par(N, -1),
		sz(N, -1),
		depth_val_to_idx(N, -1)
	{
		{
			int nxt_idx = 0;
			std::vector<int> depth_freq(N);
			std::vector<int> depth_val(N);
			auto build_one_tree = [&](int rt) -> void {
				assert(idx[rt] == -1);
				std::y_combinator([&](auto&& self, int cur, int prv, int par_idx, int d) -> int {
					int cur_idx = idx[cur] = nxt_idx++;
					preorder[cur_idx] = cur;
					par[cur_idx] = par_idx;
					sz[cur_idx] = 1;
					depth[cur_idx] = d;
					depth_val[cur_idx] = ++depth_freq[d];
					for (int nxt : adj[cur]) {
						if (nxt == prv) continue;
						sz[cur_idx] += self(nxt, cur, cur_idx, d+1);
					}
					return sz[cur_idx];
				})(rt, -1, -1, 0);
			};
			if (!roots.empty()) {
				for (int r : roots) build_one_tree(r);
				for (int i = 0; i < N; i++) {
					assert(idx[i] != -1);
				}
			} else {
				for (int rt = 0; rt < N; rt++) {
					if (idx[rt] == -1) {
						build_one_tree(rt);
						roots.push_back(rt);
					}
				}
			}
			for (int i = 1; i < N; i++) {
				depth_freq[i] += depth_freq[i-1];
			}
			for (int i = 0; i < N; i++) {
				depth_val[i] = depth_freq[depth[i]] - depth_val[i];
				assert(depth_val_to_idx[depth_val[i]] == -1);
				depth_val_to_idx[depth_val[i]] = i;
			}
			depth_val_rmq = RangeMinQuery<int>(depth_val);
		}
	}

	int dist(int a, int b) {
		if (a == b) return 0;
		a = idx[a], b = idx[b];
		if (a > b) std::swap(a, b);
		int o = depth_val_to_idx[depth_val_rmq.query(a+1, b)];
		return depth[a] + depth[b] - 2 * (depth[o] - 1);
	}

	int lca(int a, int b) {
		if (a == b) return a;
		a = idx[a], b = idx[b];
		if (a > b) std::swap(a, b);
		int o = depth_val_to_idx[depth_val_rmq.query(a+1, b)];
		return preorder[par[o]];
	}

	// query which subtree of a contains b; b must be inside a 
	// returns -1 if a == b
	int get_subtree(int a, int b) {
		if (a == b) return -1;
		a = idx[a], b = idx[b];
		assert(a < b && b < a + sz[a]);
		return preorder[depth_val_to_idx[depth_val_rmq.query(a+1, b)]];
	}

	// next from a to b, a and b must be in the same tree
	int get_next(int a, int b) {
		if (a == b) return -1;
		a = idx[a], b = idx[b];
		if (a < b && b < a + sz[a]) {
			return preorder[depth_val_to_idx[depth_val_rmq.query(a+1, b)]];
		} else {
			return preorder[par[a]];
		}
	}
};
