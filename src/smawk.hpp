#pragma once

#include <vector>
#include <cassert>
#include <optional>
#if __cpp_concepts >= 202002
#include <concepts>
#endif

namespace smawk {

template <typename T> struct value_t {
	T v;
	int col;
};

// Get(int row, int col) -> T
// Select(int row, const value_t<T>& opt_0, const value_t<T>& opt_1) returns 0 or 1 for which is better
#if __cpp_concepts >= 202002
template <typename T, typename Get, typename Select> concept totally_monotone_matrix_oracle =
	std::default_initializable<T> && std::movable<T>
	&& std::invocable<Get, int, int> && std::convertible_to<std::invoke_result_t<Get, int, int>, T>
	&& std::predicate<Select, int, const value_t<T>&, const value_t<T>&>;
#endif


template <typename Get, typename Select, typename T = std::invoke_result_t<Get, int, int>>
#if __cpp_concepts >= 202002
requires totally_monotone_matrix_oracle<T, Get, Select>
#endif
class LARSCH {
public:
	int N;
	Get get;
	Select select;
	int L;
	int num_rows;

	std::vector<std::vector<value_t<T>>> stk;
	std::vector<std::pair<value_t<T>, int>> bests;
	LARSCH() {}
	LARSCH(int N_, Get&& get_, Select&& select_) : N(N_), get(std::forward<Get>(get_)), select(std::forward<Select>(select_)) {
		L = N ? 31 - __builtin_clz(N) : 0;
		stk.resize(L);
		bests.resize(L);
		// N >> L == 1, unless N == 0
		for (int i = 0; i < L; i++) {
			stk[i].reserve(N >> (i+1));
		}
		num_rows = 0;
	}

	value_t<T> push_and_query_next() {
		assert(num_rows < N);
		int inp_row = num_rows++;

		int l = 0;
		value_t<T> nbest;
		while (true) {
			int r = inp_row >> l;
			int col = l == 0 ? inp_row : stk[l-1][r].col;
			if (r & 1) {
				int row = ((r+1) << l) - 1;
				value_t<T> prv_col_top;
				if (l == 0) prv_col_top = {get(row, col), col};
				else prv_col_top = {std::move(stk[l-1][r].v), stk[l-1][r].col};

				// just check this guy at this row, and then push it into the next layer, but don't query yet
				if (select(row, bests[l].first, prv_col_top)) {
					// prv_col_top is better here
					bests[l].first = std::move(prv_col_top);
					bests[l].second = (r+1)/2;
					// optimization: since we're the global best, we know we'll kill the entire rest of the stack, so just do it here
					assert(int(stk[l].size()) >= (r+1)/2);
					stk[l].resize((r+1)/2);
				}
			}
			if (l < L) {
				std::optional<value_t<T>> to_push;
				while (int(stk[l].size()) > (r+1)/2) {
					int row = (int(stk[l].size()) << (l+1)) - 1;
					value_t<T> nv{get(row, col), col};
					if (select(row, stk[l].back(), nv)) {
						stk[l].pop_back();
						to_push = std::move(nv);
					} else {
						break;
					}
				}
				if (to_push) {
					stk[l].emplace_back(std::move(*to_push));
				} else {
					int row = (int(stk[l].size()+1) << (l+1)) - 1;
					if (row < N) stk[l].emplace_back(get(row, col), col);
				}
			}
			if (r & 1) {
				// go return
				nbest = std::move(bests[l].first);
				l--;
				break;
			} else if (l == L) {
				// special case: just go down 1 level already
				int row = ((r+1) << l) - 1;
				if (l == 0) nbest = {get(row, col), col};
				else nbest = {std::move(stk[l-1][r].v), stk[l-1][r].col};
				l--;
				break;
			} else if (((r+2) << l) - 1 >= N) {
				// go return
				nbest.col = col;
				break;
			} else {
				l++;
				continue;
			}
			assert(false);
		}
		for (; l >= 0; l--) {
			int r = inp_row >> l;
			assert(!(r & 1));
			int row = ((r+1) << l) - 1;
			bests[l].first = std::move(nbest);
			bool did_set = false;
			while (true) {
				int idx = bests[l].second;
				int col = (l == 0 ? idx : stk[l-1][idx].col);
				value_t<T> cnd;
				if (l > 0 && idx == r) cnd = {std::move(stk[l-1][r].v), col};
				else cnd = {get(row, col), col};
				if (!did_set || select(row, nbest, cnd)) {
					did_set = true;
					nbest = std::move(cnd);
				}
				if (col == bests[l].first.col) break;
				bests[l].second++;
			}
		}
		assert(l == -1);
		return nbest;
	}
};

template <typename Get, typename Select, typename T = std::invoke_result_t<Get&&, int, int>>
#if __cpp_concepts >= 202002
requires totally_monotone_matrix_oracle<T, Get&&, Select&&>
#endif
std::vector<value_t<T>> smawk(int N, int M, Get&& get, Select&& select) {
	// TODO: If M >> N, then we should do an extra layer of column filter on the outside. The cutoff should be M > 2N or so.
	std::vector<value_t<T>> res(N);
	for (int i = 0; i < N; i++) res[i].col = -1;
	std::vector<int> stks(N);
	int L = N ? 31 - __builtin_clz(N) : 0;
	std::vector<int> stk_ends(L+1);
	stk_ends[0] = 0;
	for (int l = 0; l < L; l++) {
		int sz = 0;
		auto check_col = [&](int col, int min_sz) -> void {
			while (sz > min_sz) {
				int row = (sz << (l+1)) - 1;
				value_t<T> cnd(get(row, col), col);
				if (select(row, res[row], cnd)) {
					// we prefer cnd, save this
					res[row] = std::move(cnd);
					sz--;
				} else {
					break;
				}
			}

			if (sz < (N >> (l+1))) {
				int row = ((sz+1) << (l+1)) - 1;
				if (res[row].col == col) {
					stks[stk_ends[l] + sz] = col;
					sz++;
				} else {
					value_t<T> cnd(get(row, col), col);
					// This is a legal optimization, but I'm not sure it buys anything real, so just stub it out with true ||
					if (true || res[row].col == -1 || res[row].col < col || !select(row, cnd, res[row])) {
						res[row] = std::move(cnd);
						stks[stk_ends[l] + sz] = col;
						sz++;
					}
				}
			}
		};
		if (l == 0) {
			for (int col = 0; col < M; col++) {
				check_col(col, 0);
			}
		} else {
			for (int z = stk_ends[l-1]; z < stk_ends[l]; z++) {
				check_col(stks[z], (z - stk_ends[l-1]) / 2);
			}
		}
		assert(sz <= (N >> (l+1)));
		stk_ends[l+1] = stk_ends[l] + sz;
	}
	for (int l = L; l >= 0; l--) {
		int z = l == 0 ? 0 : stk_ends[l-1];
		for (int r = 0; r < (N >> l); r += 2) {
			int row = ((r+1) << l) - 1;
			// TODO: You could not reset this? Not sure if it buys anything real.
			res[row].col = -1;
			for (; z < (l == 0 ? M : stk_ends[l]); z++) {
				int col = l == 0 ? z : stks[z];
				value_t<T> cnd = {get(row, col), col};
				if (res[row].col == -1 || select(row, res[row], cnd)) {
					res[row] = std::move(cnd);
				}
				if ((r+1) < (N >> l) && col == res[((r+2) << l) - 1].col) break;
			}
			assert(res[row].col != -1);
		}
	}
	return res;
}

// namespace smawk
}
