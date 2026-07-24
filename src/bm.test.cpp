#include <catch2/catch_test_macros.hpp>

#include "bm.hpp"
#include "modnum.hpp"

using namespace std;

TEST_CASE("Berlekamp Massey", "[bm]") {
	using num = modnum<int(1e9)+7>;
	vector<num> S({0, 1, 1, 2, 3, 5, 8, 13});
	vector<num> tr = BerlekampMassey(S);
	REQUIRE(tr == vector<num>({num(1), num(1)}));
	num res = linearRec(S, tr, 1000);
	REQUIRE(res == num(517691607));
}
