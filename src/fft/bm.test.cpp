#include <catch2/catch_test_macros.hpp>

#include "fft/bm.hpp"
#include "num/modnum.hpp"

using namespace std;

TEST_CASE("Berlekamp Massey", "[bm]") {
	using num = wala::modnum<int(1e9)+7>;
	vector<num> S({0, 1, 1, 2, 3, 5, 8, 13});
	wala::vec<num> tr = wala::BerlekampMassey(S);
	REQUIRE(tr == wala::vec<num>({num(1), num(1)}));
	num res = wala::linearRec(S, tr, 1000);
	REQUIRE(res == num(517691607));
}
