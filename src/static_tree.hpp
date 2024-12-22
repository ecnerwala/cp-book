#pragma once

#include "yc.hpp"
#include "rmq.hpp"

struct static_forest_t {
	int N;

private:
	// original label to preorder
	std::vector<int> idx;

	// all keys/values are by preorder relabelling
	std::vector<int> preorder;
	std::vector<int> depth;
	std::vector<int> par;
	std::vector<int> sz;
	std::vector<int> heavy_par;
	std::vector<int> heavy_dist;

	std::vector<int> depth_val_to_idx;
	RangeMinQuery<int> depth_val_rmq;

public:

	static_forest_t() : N(0) {}
	static_forest_t(const std::vector<std::vector<int>>& adj, const std::vector<int>& roots = {}) :
		N(int(adj.size())),
		idx(N, -1),
		preorder(N, -1),
		depth(N, -1),
		par(N, -1),
		sz(N, -1),
		heavy_par(N, -1),
		heavy_dist(N, -1),
		depth_val_to_idx(N, -1)
	{
		{
			int nxt_idx = 0;
			std::vector<int> depth_freq(N, 0);
			std::vector<int> depth_val(N, -1);
			std::vector<int> heavy_child(N, -1);
			auto build_one_tree = [&](int rt) -> void {
				std::y_combinator([&](auto self, int cur, int prv) -> int {
					int cur_sz = 1;
					int cur_heavy = -1;
					int cur_heavy_weight = 0;
					for (int nxt : adj[cur]) {
						if (nxt == prv) continue;
						int n_sz = self(nxt, cur);
						if (n_sz > cur_heavy_weight) {
							cur_heavy = nxt;
							cur_heavy_weight = n_sz;
						}
						cur_sz += n_sz;
					}
					heavy_child[cur] = cur_heavy;
					return cur_sz;
				})(rt, -1);
				assert(idx[rt] == -1);
				std::y_combinator([&](auto&& self, int cur, int prv, int par_idx, int d, bool is_heavy_root) -> void {
					int cur_idx = idx[cur] = nxt_idx++;
					preorder[cur_idx] = cur;
					par[cur_idx] = par_idx;
					depth[cur_idx] = d;
					depth_val[cur_idx] = ++depth_freq[d];
					assert(is_heavy_root == (par_idx == -1 || cur_idx != par_idx + 1));
					if (is_heavy_root) {
						heavy_par[cur_idx] = par_idx;
						heavy_dist[cur_idx] = 1;
					} else {
						assert(par_idx == cur_idx - 1);
						heavy_par[cur_idx] = heavy_par[cur_idx - 1];
						heavy_dist[cur_idx] = heavy_dist[cur_idx - 1] + 1;
					}
					if (heavy_child[cur] != -1) {
						int nxt = heavy_child[cur];
						self(nxt, cur, cur_idx, d+1, false);
					}
					for (int nxt : adj[cur]) {
						if (nxt == prv) continue;
						if (nxt == heavy_child[cur]) continue;
						self(nxt, cur, cur_idx, d+1, true);
					}
					sz[cur_idx] = nxt_idx - cur_idx;
				})(rt, -1, -1, 0, true);
			};
			if (!roots.empty()) {
				for (int r : roots) build_one_tree(r);
			} else {
				for (int rt = 0; rt < N; rt++) {
					if (idx[rt] == -1) {
						build_one_tree(rt);
					}
				}
			}
			for (int i = 0; i < N; i++) {
				assert(idx[i] != -1);
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

	int dist(int a, int b) const {
		if (a == b) return 0;
		a = idx[a], b = idx[b];
		if (a > b) std::swap(a, b);
		int o = depth_val_to_idx[depth_val_rmq.query(a+1, b)];
		return depth[a] + depth[b] - 2 * (depth[o] - 1);
	}

	int lca(int a, int b) const {
		if (a == b) return a;
		a = idx[a], b = idx[b];
		if (a > b) std::swap(a, b);
		int o = depth_val_to_idx[depth_val_rmq.query(a+1, b)];
		return preorder[par[o]];
	}

	// query which subtree of a contains b; b must be inside a
	// returns -1 if a == b
	int get_subtree(int a, int b) const {
		if (a == b) return -1;
		a = idx[a], b = idx[b];
		assert(a < b && b < a + sz[a]);
		return preorder[depth_val_to_idx[depth_val_rmq.query(a+1, b)]];
	}

	// next from a to b, a and b must be in the same tree
	int get_next(int a, int b) const {
		if (a == b) return -1;
		a = idx[a], b = idx[b];
		if (a < b && b < a + sz[a]) {
			return preorder[depth_val_to_idx[depth_val_rmq.query(a+1, b)]];
		} else {
			return preorder[par[a]];
		}
	}

	int get_ancestor(int a, int k) const {
		assert(k >= 0);
		a = idx[a];
		if (k > depth[a]) return -1;
		while (a != -1 && k > 0) {
			if (k >= heavy_dist[a]) {
				k -= heavy_dist[a];
				assert(heavy_par[a] <= a - heavy_dist[a]);
				a = heavy_par[a];
			} else {
				a -= k;
				k = 0;
			}
		}
		return preorder[a];
	}

	int get_depth(int a) const { return depth[idx[a]]; }
	int get_sz(int a) const { return sz[idx[a]]; }
	std::array<int, 2> get_range(int a) const { return {idx[a], idx[a] + sz[idx[a]]}; }
	bool is_ancestor(int a, int b) { return idx[a] <= idx[b] && idx[b] < idx[a] + sz[idx[a]]; }
};
