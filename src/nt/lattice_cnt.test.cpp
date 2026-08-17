#include <catch2/catch_test_macros.hpp>

#include "lattice_cnt.hpp"

using namespace std;

long long lattice_cnt_slow(long long A, long long B, long long C) {
	using ll = long long;
	ll ans = 0;
	for (ll x = 0; A * x <= C; x++) {
		for (ll y = 0; A * x + B * y <= C; y++) {
			ans++;
		}
	}
	return ans;
}

long long mod_count_range_slow(long long a, long long m, long long clo, long long chi, long long nlo, long long nhi) {
	assert(nlo <= nhi);
	assert(clo <= chi);
	long long ans = 0;
	for (long long i = nlo; i < nhi; i++) {
		for (long long j = clo; j < chi; j++) {
			ans += (((a * i - j) % m) == 0);
		}
	}
	return ans;
}

TEST_CASE("Lattice Count", "[lattice_cnt]") {
	for (int a = 0; a <= 50; a++) {
		for (int b = 0; b <= 10; b++) {
			for (int c = -1; c <= 100; c++) {
				if ((a == 0 || b == 0) && c >= 0) continue;
				INFO("a = " << a);
				INFO("b = " << b);
				INFO("c = " << c);
				REQUIRE(lattice_cnt(a, b, c) == lattice_cnt_slow(a, b, c));
			}
		}
	}
}

TEST_CASE("Mod Count (positive)", "[lattice_cnt]") {
	for (int m = 1; m <= 25; m++) {
		for (int a = 0; a <= m+10; a++) {
			for (int c = 0; c <= m; c++) {
				INFO("a = " << a);
				INFO("m = " << m);
				INFO("c = " << c);
				int trueAns = 0;
				for (int n = 1; n <= m+10; n++) {
					INFO("n = " << n);

					trueAns += (a * (n-1) % m) < c;
					REQUIRE(mod_count(a, m, c, n) == trueAns);
				}
			}
		}
	}
}

TEST_CASE("Mod Count (negatives)", "[lattice_cnt]") {
	for (int m : {1, 2, 3, 5, 8, 13, 21}) {
		for (int a : {-10, 0, 1, 2, 3, 5, m, m+5}) {
			auto cnds = {-37, -2*m-1, -m, -m+1, -m/2, -1, 0, 1, m/2, m+1, 2*m-1, 34};
			INFO("a = " << a);
			INFO("m = " << m);
			for (int clo : cnds) {
				for (int nlo : cnds) {
					INFO("clo = " << clo);
					INFO("nlo = " << nlo);
					REQUIRE(mod_count_range(a, m, clo, 47, nlo, 49) == mod_count_range_slow(a, m, clo, 47, nlo, 49));
				}
			}

			for (int chi : cnds) {
				for (int nhi : cnds) {
					INFO("chi = " << chi);
					INFO("nhi = " << nhi);
					REQUIRE(mod_count_range(a, m, -55, chi, -57, nhi) == mod_count_range_slow(a, m, -55, chi, -57, nhi));
				}
			}
		}
	}
}
