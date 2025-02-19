#pragma once

#include <vector>
#include <cassert>

template <typename T> struct smawk_value_t {
	T v;
	int col;
};

// cmp returns if the first is better
template <typename T, typename Get, typename Select> class LARSCH {
public:
	int N;
	Get get;
	Select select;
	int L;
	int num_rows;

	struct stack_entry_t {
		smawk_value_t<T> v;
		int prv_idx;
	};
	std::vector<std::vector<stack_entry_t>> stk;
	LARSCH() {}
	LARSCH(int N_, Get&& get_, Select&& select_) : N(N_), get(std::forward<Get>(get_)), select(std::forward<Select>(select_)) {
		L = N ? 31 - __builtin_clz(N) : 0;
		stk.resize(L);
		for (int i = 0; i < L; i++) {
			stk[i].reserve(N >> (i+1));
		}
		num_rows = 0;
	}

	smawk_value_t<T> push_and_query_next() {
		assert(num_rows < N);
		int inp_row = num_rows++;
		int inp_col = inp_row;

		for (int l = 0, r = inp_row, col = inp_col; true; l++, r >>= 1) {
			if (r & 1) {
				int row = ((r+1) << l) - 1;
				// just check this guy at this row, and then push it into the next layer, but don't query yet
				smawk_value_t<T> v = l > 0 ? stk[l-1][p_idx].v : smawk_value_t<T>{get(row, inp_col), inp_col};
				if (!select(inp_row, stk[l][r/2].v, v)) {
					stk[l][r/2].v = std::move(v);
					stk[l][r/2].prv_idx = p_idx;
				}
				// TODO: push into next layer
			} else {
			}
		}
	}
};

template <typename Get, typename Select>
LARSCH(int, Get&&, Select&&) -> LARSCH<std::invoke_result_t<Get, int, int>, Get, Select>;

template <typename Get, typename Select> std::vector<smawk_value_t<std::invoke_result_t<Get&&, int, int>>> smawk(int N, int M, Get&& get, Select&& select) {
	// TODO: If M >> N, then we should do an extra layer of column filter on the outside. The cutoff should be M > 2N or so.
	using result_t = std::invoke_result_t<Get&&, int, int>;
	std::vector<smawk_value_t<result_t>> res(N, {result_t(), -1});
	std::vector<int> stks(N);
	int L = N ? 31 - __builtin_clz(N) : 0;
	std::vector<int> stk_ends(L+1);
	stk_ends[0] = 0;
	for (int l = 0; l < L; l++) {
		int sz = 0;
		auto check_col = [&](int col, int min_sz) -> void {
			while (sz > min_sz) {
				int row = (sz << (l+1)) - 1;
				smawk_value_t<result_t> cnd(get(row, col), col);
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
					smawk_value_t<result_t> cnd(get(row, col), col);
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
				smawk_value_t<result_t> cnd = {get(row, col), col};
				if (res[row].col == -1 || select(row, res[row], cnd)) {
					res[row] = cnd;
				}
				if ((r+1) < (N >> l) && col == res[((r+2) << l) - 1].col) break;
			}
			assert(res[row].col != -1);
		}
	}
	return res;
}
