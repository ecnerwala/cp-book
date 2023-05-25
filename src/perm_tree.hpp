#pragma once

#include <vector>
#include <array>
#include <cassert>

class PermTree {
public:
	enum class NodeType {
		LEAF,
		INCR,
		DECR,
		FULL,
		PARTIAL,
	};

	struct Node {
		std::array<int, 2> c;
		NodeType type;
		int l, r, lo, hi;
	};

	std::vector<Node> nodes;
	int root = -1;

	PermTree() {}
	Node& operator [] (int idx) { return nodes[idx]; }
	const Node& operator [] (int idx) const { return nodes[idx]; }

	int size() const { return int(nodes.size()); }

	PermTree(const std::vector<int>& A) : nodes(int(A.size())*2-1) {
		int N = int(A.size());
		std::vector<int> nxt_earlier(N);
		std::vector<int> prv_earlier(N);
		for (int i = 0; i < N; i++) {
			nxt_earlier[i] = i+1;
			prv_earlier[i] = i-1;
		}
		for (int i = N-1; i >= 0; i--) {
			int a = A[i];
			int p = prv_earlier[a];
			int n = nxt_earlier[a];
			if (p != -1) nxt_earlier[p] = n;
			if (n != N) prv_earlier[n] = p;
		}

		struct cnd_t {
			int left;
			int lo;
			int lo_gap;
			int hi;
			int hi_gap;
			int node;
		};

		std::vector<cnd_t> stk; stk.reserve(N);

		for (int i = 0; i < N; i++) {
			int a = A[i];
			while (true) {
				if (!stk.empty() && (a < stk.back().lo_gap || a > stk.back().hi_gap)) {
					assert(stk.size() >= 2);
					stk.end()[-2].lo = std::min(stk.end()[-2].lo, stk.back().lo);
					stk.end()[-2].hi = std::max(stk.end()[-2].hi, stk.back().hi);

					int n = 2 * stk.back().left - 1;
					nodes[n].c = {stk.end()[-2].node, stk.end()[-1].node};
					nodes[n].type = NodeType::PARTIAL;
					nodes[n].l = stk.end()[-2].left;
					nodes[n].r = i-1;
					nodes[n].lo = stk.end()[-2].lo;
					nodes[n].hi = stk.end()[-2].hi;

					stk.pop_back();

					stk.back().node = n;
				} else {
					break;
				}
			}

			stk.push_back({i, a, prv_earlier[a]+1, a, nxt_earlier[a]-1, 2*i});
			nodes[2*i].type = NodeType::LEAF;
			nodes[2*i].c = {-1, -1};
			nodes[2*i].l = nodes[2*i].r = i;
			nodes[2*i].lo = nodes[2*i].hi = a;

			while (stk.size() >= 2 && std::max(stk.back().hi, stk.end()[-2].hi) - std::min(stk.back().lo, stk.end()[-2].lo) == i - stk.end()[-2].left) {
				// merge these two nodes into one
				stk.end()[-2].lo = std::min(stk.end()[-2].lo, stk.back().lo);
				stk.end()[-2].hi = std::max(stk.end()[-2].hi, stk.back().hi);

				int n = 2 * stk.back().left - 1;
				nodes[n].c = {stk.end()[-2].node, stk.end()[-1].node};
				if (stk.end()[-2].lo == stk.end()[-1].lo) {
					nodes[n].type = NodeType::DECR;
				} else if (stk.end()[-2].hi == stk.end()[-1].hi) {
					nodes[n].type = NodeType::INCR;
				} else {
					nodes[n].type = NodeType::FULL;
				}
				nodes[n].l = stk.end()[-2].left;
				nodes[n].r = i;
				nodes[n].lo = stk.end()[-2].lo;
				nodes[n].hi = stk.end()[-2].hi;

				stk.pop_back();
				stk.back().node = n;
			}
		}

		assert(stk.size() == 1);
		root = stk.back().node;
	}
};
