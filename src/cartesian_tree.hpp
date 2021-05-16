#pragma once

#include <vector>
#include <array>

#include "reverse_comparator.hpp"

class CartesianTree {
public:
	struct Node {
		int l, m, r; // inclusive ranges
		std::array<Node*, 2> c;
	};
	std::vector<Node> nodes;
	Node* root = nullptr;

	CartesianTree() {}

	// We prohibit copies because naively copying would invalidate pointers.
	CartesianTree(const CartesianTree& o) = delete;
	CartesianTree& operator= (const CartesianTree& o) = delete;

	// Note that the moved-from tree may have an invalid root.
	CartesianTree(CartesianTree&& o) noexcept = default;
	CartesianTree& operator= (CartesianTree&& o) noexcept = default;

private:
	CartesianTree(std::vector<Node>&& nodes_, Node* root_) : nodes(std::move(nodes_)), root(root_) {}

public:

	// min-cartesian-tree, with earlier cells tiebroken earlier
	template <typename T, typename Comp = std::less<T>>
	static CartesianTree build_min_tree(const std::vector<T>& v, Comp comp = Comp()) {
		std::vector<Node> nodes(v.size()*2+1);
		std::vector<Node*> stk; stk.reserve(v.size());
		Node* root = nullptr;
		for (int i = 0; i <= int(v.size()); i++) {
			nodes[2*i].l = i;
			nodes[2*i].r = i-1;
			nodes[2*i].m = i-1;
			nodes[2*i].c = {nullptr, nullptr};
			Node* cur = &nodes[2*i];
			while (!stk.empty() && (i == int(v.size()) || comp(v[i], v[stk.back()->m]))) {
				Node* nxt = stk.back(); stk.pop_back();
				nxt->c[1] = cur;
				nxt->r = cur->r;
				cur = nxt;
			}
			if (i == int(v.size())) {
				root = cur;
				break;
			}
			nodes[2*i+1].l = cur->l;
			nodes[2*i+1].m = i;
			nodes[2*i+1].c[0] = cur;
			stk.push_back(&nodes[2*i+1]);
		}
		return {std::move(nodes), root};
	}

	// max-cartesian-tree, with earlier cells tiebroken earlier
	template <typename T, typename Comp = std::less<T>>
	static CartesianTree build_max_tree(const std::vector<T>& v, Comp comp = Comp()) {
		return build_min_tree(v, reverse_comparator(comp));
	}
};
