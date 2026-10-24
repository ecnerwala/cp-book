#pragma once

#include <ranges>
#include <array>

#include "containers/vec.hpp"
#include "reverse_comparator.hpp"

namespace wala {

class CartesianTree {
public:
	struct Node {
		int l, m, r; // inclusive ranges
		std::array<int, 2> c;
		int p;
	};
	vec<Node> nodes;
	int root = -1;

	CartesianTree() {}

	Node& operator [] (int idx) { return nodes[idx]; }
	const Node& operator [] (int idx) const { return nodes[idx]; }

	int size() const { return nodes.size(); }

private:
	CartesianTree(vec<Node>&& nodes_, int root_) : nodes(std::move(nodes_)), root(root_) {}

public:
	// min-cartesian-tree, with earlier cells tiebroken earlier
	template <std::ranges::random_access_range R, typename Comp = std::ranges::less>
	static CartesianTree build_min_tree(const R& v, Comp comp = Comp()) {
		int N = int(std::ranges::size(v));
		vec<Node> nodes(2*N+1);
		bounded_vec<int> stk(with_capacity{N});
		int root = -1;
		for (int i = 0; i <= N; i++) {
			int cur = 2*i;
			nodes[cur].l = i;
			nodes[cur].r = i-1;
			nodes[cur].m = i-1;
			nodes[cur].c = {-1, -1};
			while (!stk.empty() && (i == N || comp(v[i], v[nodes[stk.back()].m]))) {
				int nxt = stk.back(); stk.pop_back();
				nodes[cur].p = nxt;
				nodes[nxt].c[1] = cur;
				nodes[nxt].r = nodes[cur].r;
				cur = nxt;
			}
			if (i == N) {
				root = cur;
				break;
			}
			nodes[2*i+1].l = nodes[cur].l;
			nodes[2*i+1].m = i;
			nodes[cur].p = 2*i+1;
			nodes[2*i+1].c[0] = cur;
			stk.push_back(2*i+1);
		}
		nodes[root].p = -1;
		return {std::move(nodes), root};
	}

	// max-cartesian-tree, with earlier cells tiebroken earlier
	template <std::ranges::random_access_range R, typename Comp = std::ranges::less>
	static CartesianTree build_max_tree(const R& v, Comp comp = Comp()) {
		return build_min_tree(v, reverse_comparator(comp));
	}
};

} // namespace wala
