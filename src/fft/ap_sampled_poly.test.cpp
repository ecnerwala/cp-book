#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_get_random_seed.hpp>

#include "fft/ap_sampled_poly.hpp"
#include "fft/engines/ntt.hpp"
#include "modnum.hpp"

namespace ecnerwala {
namespace fft {

using namespace std;

TEST_CASE("ap_sampled_poly eval", "[fft,ap_sampled_poly]") {
	using num = modnum<998244353>;
	using poly_vals = ap_sampled_poly<engines::ntt<num>>;
	mt19937 mt(Catch::getSeed());
	for (int len : {0, 1, 2, 3, 5, 8, 13, 21}) {
		INFO("len = " << len);
		std::vector<num> coeffs(len);
		for (int i = 0; i < len; i++) coeffs[i] = num(mt());
		auto eval_at = [&](num v) {
			num r = 0;
			for (int i = len-1; i >= 0; i--) r = r * v + coeffs[i];
			return r;
		};
		poly_vals v(len);
		for (int i = 0; i < len; i++) v[i] = eval_at(i);
		for (int i = -2 * len; i <= 2 * len; i++) {
			REQUIRE(v.eval_at(i) == eval_at(i));
		}
		auto eval_range = [&](num k, int osz) {
			poly_vals r(osz);
			for (int i = 0; i < osz; i++) {
				r[i] = eval_at(k + num(i));
			}
			return r;
		};
		num k = 1023895;
		for (int osz : {0, 1, 2, 3, 5, 8, 13, 21}) {
			INFO("osz = " << osz);
			REQUIRE(v.eval_range(k, osz) == eval_range(k, osz));
		}
	}
}

}} // namespace ecnerwala::fft
