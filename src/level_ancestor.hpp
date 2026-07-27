#pragma once

#include <algorithm>
#include <utility>
#include <vector>
#include <cassert>

#include "newtype.hpp"
#include "yc.hpp"

namespace ecnerwala {

using std::swap;

// Node is the external node-id space (a newtype or plain int); preorder
// indices live in their own internal newtype space pre_t.
template <IntKey Node = int> struct level_ancestor {
	struct pre_tag_;
	using pre_t = int_newtype<pre_tag_>;

	pre_t N;
	vec<pre_t, Node> preorder;
	vec<Node, pre_t> idx;
	vec<pre_t, std::pair<pre_t, pre_t>> heavyPar; // heavy parent, distance
	level_ancestor() : N(0) {}

	explicit level_ancestor(const vec<Node, Node>& par) : N(pre_t(par.size())), preorder(N), idx(par.size()), heavyPar(N) {
		vec<Node, std::vector<Node>> ch(par.size());
		for (Node i = 0; i < par.size(); ++i) {
			if (par[i] != -1) ch[par[i]].push_back(i);
		}
		vec<Node, int> sz(par.size());
		pre_t nxt_idx = 0;
		for (Node i = 0; i < par.size(); ++i) {
			if (par[i] == -1) {
				std::y_combinator([&](auto self, Node cur) -> void {
					sz[cur] = 1;
					for (Node nxt : ch[cur]) {
						self(nxt);
						sz[cur] += sz[nxt];
					}
					if (!ch[cur].empty()) {
						auto mit = std::max_element(ch[cur].begin(), ch[cur].end(), [&](Node a, Node b) { return sz[a] < sz[b]; });
						swap(*ch[cur].begin(), *mit);
					}
				})(i);
				std::y_combinator([&](auto self, Node cur, bool isRoot = true) -> void {
					preorder[idx[cur] = nxt_idx++] = cur;
					if (isRoot) {
						heavyPar[idx[cur]] = {par[cur] == -1 ? pre_t(-1) : idx[par[cur]], 1};
					} else {
						assert(idx[par[cur]] == idx[cur] - 1);
						heavyPar[idx[cur]] = heavyPar[idx[cur] - 1];
						heavyPar[idx[cur]].second++;
					}
					bool chRoot = false;
					for (Node nxt : ch[cur]) {
						self(nxt, chRoot);
						chRoot = true;
					}
				})(i);
			}
		}
	}
	explicit level_ancestor(std::vector<Node> par) : level_ancestor(vec<Node, Node>(std::move(par))) {}

	Node get_ancestor(Node n, pre_t k) const {
		assert(k >= 0);
		pre_t a = idx[n];
		while (a != -1 && k != 0) {
			if (k >= heavyPar[a].second) {
				k -= heavyPar[a].second;
				assert(heavyPar[a].first <= a - heavyPar[a].second);
				a = heavyPar[a].first;
			} else {
				a -= k;
				k = 0;
			}
		}
		if (a == -1) return Node(-1);
		else return preorder[a];
	}

	Node lca(Node a_, Node b_) const {
		pre_t a = idx[a_], b = idx[b_];
		while (true) {
			if (a > b) swap(a, b);
			assert(a <= b);
			if (a > b - heavyPar[b].second) {
				return preorder[a];
			}
			b = heavyPar[b].first;
			if (b == -1) return Node(-1);
		}
	}

	pre_t dist(Node a_, Node b_) const {
		pre_t a = idx[a_], b = idx[b_];
		pre_t res = 0;
		while (true) {
			if (a > b) swap(a, b);
			assert(a <= b);
			if (a > b - heavyPar[b].second) {
				res += b - a;
				break;
			}
			res += heavyPar[b].second;
			b = heavyPar[b].first;
			if (b == -1) return -1;
		}
		return res;
	}
};

} // namespace ecnerwala
