#pragma once

#include <vector>
#include <array>
#include <cassert>

// Finds an optimal alphabetic (binary) Huffman code, i.e. one that preserves the ordering of the original weights
// Implements the Garsia-Wachs algorithm: https://en.wikipedia.org/wiki/Garsia%E2%80%93Wachs_algorithm
// Returns the code specified as a sequence of depths for each input weight
template <typename T, typename T_sum = T> std::vector<int> alphabetic_huffman_code(std::vector<T> weights) {
	int N = int(weights.size());
	if (N == 0) return {};
	std::vector<std::array<int, 2>> ch; ch.reserve(N-1);

	{
		struct splay_node {
			mutable splay_node* p = nullptr;
			std::array<splay_node*, 2> c{nullptr, nullptr};
			int d() const { return this == p->c[1]; }

			T_sum value;
			T_sum max_value;
			int idx;

			void update() {
				max_value = value;
				for (auto ch : c) {
					if (ch && max_value < ch->max_value) max_value = ch->max_value;
				}
			}

			void rot() {
				assert(p);

				int x = d();
				splay_node* pa = p;
				splay_node* ch = c[!x];

				if (ch) ch->p = pa;
				pa->c[x] = ch;

				if (pa->p) pa->p->c[pa->d()] = this;
				this->p = pa->p;

				this->c[!x] = pa;
				pa->p = this;

				pa->update();
			}

			void splay_no_update(splay_node* top) {
				while (p != top) {
					if (p->p != top) {
						if (p->d() == d()) p->rot();
						else rot();
					}
					rot();
				}
			}
		};
		std::vector<splay_node> nodes(N+1);
		for (int i = 0; i < N; i++) {
			nodes[i].p = &nodes[i+1];
			nodes[i+1].c[0] = &nodes[i];
			nodes[i].value = T_sum(weights[i]);
			nodes[i].idx = i;
		}
		nodes[0].update();
		splay_node* cur = &nodes[1];

		// We'll store our current state as the left spine of some splay tree.
		// All vertices from cur to the root are precisely the vertices that may satisfy w[n-2] <= w[n]
		// (all others provably satisfy w[x-2] > w[x] at all times),
		// so cur is exactly the leftmost vertex that might satisfy w[n-2] <= w[n].
		//
		// We then check this condition, and if it does have w[n-2] <= w[n],
		// we merge w[n-2] and w[n-1] and reinsert somewhere according to Garsia-Wachs,
		// i.e. right after the last element of w[0:n-1] greater than or equal to it.
		// Then, the newly inserted node is added to the candidate chain
		// (exercise: prove that all other positions still satisfy w[x-2] > w[x]).

		while (cur) {
			// Note: cur is not necessarily updated

			// First, grab the 2nd child of the left side of cur
			splay_node* a = cur->c[0];
			assert(a);
			while (a->c[1]) a = a->c[1];
			if (a->c[0]) {
				a = a->c[0];
				while (a->c[1]) a = a->c[1];
			} else {
				a = a->p;
			}
			if (a == cur) {
				// size one, so we're done
				cur->update();
				cur = cur->p;
				continue;
			}
			a->splay_no_update(cur);
			assert(a == cur->c[0]);
			assert(a->c[1] && !a->c[1]->c[0] && !a->c[1]->c[1]);
			if (cur->p && cur->value < a->value) {
				// no merging, so we're done
				a->update();
				cur->update();
				cur = cur->p;
				continue;
			}

			// Otherwise, merge a and a->c[1]
			{
				int n_idx = N + int(ch.size());
				ch.push_back({a->idx, a->c[1]->idx});
				a->idx = n_idx;
			}
			a->value += a->c[1]->value;
			a->c[1]->p = nullptr;
			a->c[1] = nullptr;

			// Now, insert a right after the first guy b which is b.v >= a.v
			if (!a->c[0] || a->c[0]->max_value < a->value) {
				a->c[1] = a->c[0];
				a->c[0] = nullptr;
				a->update();
				// Don't recurse on a, since it has no left child
				continue;
			}

			splay_node* b = a->c[0];
			while (true) {
				assert(b);
				assert(!(b->max_value < a->value));
				if (!b->c[1] || b->c[1]->max_value < a->value) {
					if (b->value < a->value) {
						assert(b->c[0]);
						b = b->c[0];
					} else {
						break;
					}
				} else {
					b = b->c[1];
				}
			}
			b->splay_no_update(a);
			assert(b == a->c[0]);
			if (b->c[1]) b->c[1]->p = a;
			a->c[1] = b->c[1];
			b->c[1] = nullptr;
			b->update();
			cur = a;
			continue;
		}
	}

	// Reconstruct depths
	assert(int(ch.size()) == N-1);
	std::vector<int> res(2*N-1, -1);
	res[2*N-2] = 0;
	for (int i = 2*N-2; i >= N; i--) {
		assert(res[i] != -1);
		res[ch[i-N][0]] = res[i] + 1;
		res[ch[i-N][1]] = res[i] + 1;
	}
	res.resize(N);
	return res;
}

// Returns the lca array of length N - 1, suitable for building a Cartesian tree
inline std::vector<int> binary_code_depths_to_lca_depths(std::vector<int> depths) {
	int N = int(depths.size());
	if (N == 0) return {};
	std::vector<int> res; res.reserve(N-1);
	std::vector<int> stk; stk.reserve(N);
	for (int v : depths) {
		while (!stk.empty() && stk.back() == v) {
			stk.pop_back();
			v--;
		}
		assert(stk.empty() || stk.back() < v);
		if (v != 0) res.push_back(v-1);
		stk.push_back(v);
	}
	assert(int(stk.size()) == 1 && stk.back() == 0);
	return res;
}
