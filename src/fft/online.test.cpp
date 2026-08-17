#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_get_random_seed.hpp>

#include "fft/online.hpp"
#include "fft/test_util.test.hpp"
#include "num/modnum.hpp"

namespace ecnerwala {
namespace fft {

using namespace std;

TEMPLATE_TEST_CASE("online multiplier", "[fft]", MOD_ENGINES) {
	using E = TestType;
	using num = typename E::value_type;
	mt19937 mt(Catch::getSeed());
	for (int n : {1, 2, 3, 17, 64, 100}) {
		INFO("n = " << n);
		vector<num> f(n), g(n);
		fill_rnd(f, mt);
		fill_rnd(g, mt);
		auto slow = multiply_slow(f, g);
		slow.resize(2*n, num(0));
		online_multiplier<E> om(n);
		for (int i = 0; i < n; i++) {
			om.push(f[i], g[i]);
			REQUIRE(om.back() == slow[i]);
		}
		if ((n & (n-1)) == 0) {
			// the final push of a power-of-two N completes all 2N terms
			for (int i = n; i < 2*n; i++) {
				REQUIRE(om.res[i] == slow[i]);
			}
		}
	}
}

TEMPLATE_TEST_CASE("online squarer", "[fft]", MOD_ENGINES) {
	using E = TestType;
	using num = typename E::value_type;
	mt19937 mt(Catch::getSeed());
	for (int n : {1, 2, 3, 17, 64, 100}) {
		INFO("n = " << n);
		vector<num> f(n);
		fill_rnd(f, mt);
		auto slow = multiply_slow(f, f);
		slow.resize(2*n, num(0));
		online_squarer<E> os(n);
		for (int i = 0; i < n; i++) {
			os.push(f[i]);
			REQUIRE(os.back() == slow[i]);
		}
		if ((n & (n-1)) == 0) {
			for (int i = n; i < 2*n; i++) {
				REQUIRE(os.res[i] == slow[i]);
			}
		}
	}
}

// marker engine: same ring, but declared non-commutative so the squarer
// must take the multiplier path (no cross-term doubling)
struct nc_engine : engines::ntt<modnum<998244353>> {
	static constexpr bool commutative = false;
};
static_assert(fft::engine<nc_engine>);

TEST_CASE("online squarer non-commutative fallback", "[fft]") {
	using num = modnum<998244353>;
	mt19937 mt(Catch::getSeed());
	for (int n : {1, 2, 3, 17, 64, 100}) {
		INFO("n = " << n);
		vector<num> f(n);
		fill_rnd(f, mt);
		auto slow = multiply_slow(f, f);
		slow.resize(2*n, num(0));
		online_squarer<nc_engine> os(n);
		for (int i = 0; i < n; i++) {
			os.push(f[i]);
			REQUIRE(os.back() == slow[i]);
		}
	}
}

}} // namespace ecnerwala::fft
