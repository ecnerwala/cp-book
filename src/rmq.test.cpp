#include <catch2/catch_test_macros.hpp>

#include "rmq.hpp"

#include <vector>
#include <random>

TEST_CASE("RangeMinQuery", "[rmq]") {
	std::mt19937 mt(48);
	for (int N : {1, 2, 3, 5, 10, 20, 33, 48, 100, 163, 512}) {
		std::vector<std::pair<int, int>> data(N);
		for (int i = 0; i < N; i++) {
			data[i] = {mt(), i};
		}

		RangeMinQuery<std::pair<int, int>> minQ(data);
		RangeMaxQuery<std::pair<int, int>> maxQ(data);

		for (int l = 0; l < N; l++) {
			std::pair<int, int> cur_min = data[l];
			std::pair<int, int> cur_max = data[l];
			for (int r = l; r < N; r++) {
				cur_min = min(cur_min, data[r]);
				REQUIRE(minQ.query(l, r) == cur_min);
				cur_max = max(cur_max, data[r]);
				REQUIRE(maxQ.query(l, r) == cur_max);
			}
		}
	}
}
