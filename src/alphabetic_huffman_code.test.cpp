#include "alphabetic_huffman_code.hpp"

#include <catch2/catch_test_macros.hpp>
#include <cassert>
#include <random>

template <typename T> T alphabetic_huffman_code_naive(std::vector<T> weights) {
	int N = int(weights.size());
	if (N == 0) return 0;
	assert(N > 0);
	std::vector<T> dp(N * N);
	for (int i = 0; i < N; i++) {
		dp[i * N + i] = 0;
		T pref = weights[i];
		for (int j = i-1; j >= 0; j--) {
			T v = dp[i * N + i] + dp[(i-1) * N + j];
			for (int k = i-1; k >= j+1; k--) {
				v = std::min(v, dp[i * N + k] + dp[(k-1) * N + j]);
			}
			pref += weights[j];
			v += pref;
			dp[i * N + j] = v;
		}
	}
	return dp[(N-1) * N + 0];
}

TEST_CASE("Alphabetic Huffman Code", "[alphabetic_huffman_code]") {
	std::mt19937 mt(48);
	for (int z = 0; z <= 100000; z++) {
		int N = std::uniform_int_distribution(1, 15)(mt);
		std::vector<int> weights(N);
		for (auto& w : weights) w = std::uniform_int_distribution(0, 15)(mt);
		auto naive_tot = alphabetic_huffman_code_naive(weights);

		auto code_depths = alphabetic_huffman_code(weights);
		auto code_lcp = binary_code_depths_to_lca_depths(code_depths);
		int tot = 0;
		for (int i = 0; i < N; i++) {
			tot += weights[i] * code_depths[i];
		}
		REQUIRE(naive_tot == tot);
	}
}
