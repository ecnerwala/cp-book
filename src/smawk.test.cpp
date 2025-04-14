#include "smawk.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_get_random_seed.hpp>

#include <random>

std::vector<std::vector<int>> generate_totally_monotone(int N, int M, auto&& rng) {
	std::vector<int> cur_order(M);
	std::iota(cur_order.begin(), cur_order.end(), 0);
	std::vector<int> cur_vals(M);
	std::iota(cur_vals.begin(), cur_vals.end(), 0);
	std::vector<std::vector<int>> perms; perms.reserve(M * (M-1) / 2 + 1);
	perms.push_back(cur_vals);
	{
		std::vector<int> cnds; cnds.reserve(M);
		for (int z = 0; z < M * (M-1) / 2; z++) {
			cnds.clear();
			for (int i = 0; i+1 < M; i++) {
				if (cur_order[i] < cur_order[i+1]) {
					cnds.push_back(i);
				}
			}
			assert(!cnds.empty());
			int i = cnds[std::uniform_int_distribution<int>(0, int(cnds.size()) - 1)(rng)];
			std::swap(cur_order[i], cur_order[i+1]);
			cur_vals[cur_order[i]] = i;
			cur_vals[cur_order[i+1]] = i+1;
			perms.push_back(cur_vals);
		}
	}
	std::vector<std::vector<int>> output; output.reserve(N);
	std::vector<int> stars_bars(M * (M-1) / 2 + N);
	std::fill(stars_bars.begin(), stars_bars.begin() + N, 1);
	std::shuffle(stars_bars.begin(), stars_bars.end(), rng);
	{
		int perm_idx = 0;
		for (auto op : stars_bars) {
			if (op == 0) {
				perm_idx++;
			} else {
				output.push_back(perms[perm_idx]);
			}
		}
	}
	assert(int(output.size()) == N);
	return output;
}

void check_smawk(int N, int M, std::vector<std::vector<int>> mat) {
	auto result = smawk::smawk(N, M, [&](int row, int col) -> int {
		REQUIRE(0 <= row); REQUIRE(row < N);
		REQUIRE(0 <= col); REQUIRE(col < M);
		return mat[row][col];
	}, [&](int r, smawk::value_t<int> cnd1, smawk::value_t<int> cnd2) -> bool {
		REQUIRE(0 <= r); REQUIRE(r < N);
		REQUIRE(0 <= cnd1.col); REQUIRE(cnd1.col < M);
		REQUIRE(0 <= cnd2.col); REQUIRE(cnd2.col < M);
		REQUIRE(cnd1.col < cnd2.col);
		REQUIRE(cnd1.v == mat[r][cnd1.col]);
		REQUIRE(cnd2.v == mat[r][cnd2.col]);
		// cnd2 is better when it's strictly smaller
		return cnd2.v < cnd1.v;
	});
	REQUIRE(int(result.size()) == N);
	for (int i = 0; i < N; i++) {
		int j = int(std::min_element(mat[i].begin(), mat[i].end()) - mat[i].begin());
		REQUIRE(result[i].col == j);
		REQUIRE(result[i].v == mat[i][j]);
	}
}

TEST_CASE("SMAWK", "[smawk]") {
	std::mt19937 mt(Catch::getSeed());
	for (int N : {0, 1, 2, 3, 5, 8, 13}) {
		for (int M : {0, 1, 2, 3, 5, 8, 13}) {
			if (N > 0 && M == 0) continue;
			auto inp = generate_totally_monotone(N, M, mt);
			CAPTURE(N, M, inp);
			check_smawk(N, M, inp);
		}
	}
}

struct move_only_t {
	int v;
	move_only_t() : v(-1) {}
	explicit move_only_t(int v_) : v(v_) {
		assert(v_ != -1);
	}
	move_only_t(move_only_t&& o) {
		v = o.v;
		o.v = -1;
	}
	move_only_t& operator = (move_only_t&& o) {
		v = o.v;
		o.v = -1;
		return *this;
	}
	move_only_t(const move_only_t& o) = delete;
	move_only_t& operator = (const move_only_t& o) = delete;
};

void check_larsch(int N, std::vector<std::vector<int>> mat) {
	const int M = N;
	smawk::LARSCH l(N, [&](int row, int col) -> move_only_t {
		REQUIRE(0 <= row); REQUIRE(row < N);
		REQUIRE(0 <= col); REQUIRE(col < M);
		REQUIRE(col <= row);
		return move_only_t{mat[row][col]};
	}, [&](int r, const smawk::value_t<move_only_t>& cnd1, const smawk::value_t<move_only_t>& cnd2) -> bool {
		REQUIRE(0 <= r); REQUIRE(r < N);
		REQUIRE(0 <= cnd1.col); REQUIRE(cnd1.col < M);
		REQUIRE(0 <= cnd2.col); REQUIRE(cnd2.col < M);
		REQUIRE(cnd1.col < cnd2.col);
		REQUIRE(cnd1.col <= r);
		REQUIRE(cnd2.col <= r);
		REQUIRE(cnd1.v.v == mat[r][cnd1.col]);
		REQUIRE(cnd2.v.v == mat[r][cnd2.col]);
		// cnd2 is better when it's strictly smaller
		return cnd2.v.v < cnd1.v.v;
	});
	for (int i = 0; i < N; i++) {
		auto res = l.push_and_query_next();
		int j = int(std::min_element(mat[i].begin(), mat[i].begin() + i + 1) - mat[i].begin());
		REQUIRE(res.col == j);
		REQUIRE(res.v.v == mat[i][j]);
	}
}

TEST_CASE("LARSCH", "[smawk]") {
	std::mt19937 mt(Catch::getSeed());
	for (int N : {0, 1, 2, 3, 5, 8, 13}) {
		auto inp = generate_totally_monotone(N, N, mt);
		for (int i = 0; i < N; i++) {
			for (int j = i+1; j < N; j++) {
				inp[i][j] = -1;
			}
		}
		CAPTURE(N, inp);
		check_larsch(N, inp);
	}
}
