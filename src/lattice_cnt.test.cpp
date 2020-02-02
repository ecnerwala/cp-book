#include "catch.hpp"

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

long long mod_count_slow(long long a, long long m, long long c, long long n) {
	long long ans = 0;
	for (int i = 0; i < n; i++) {
		ans += (a * i % m) < c;
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

TEST_CASE("Mod Count", "[lattice_cnt]") {
	for (int m = 1; m <= 25; m++) {
		for (int a = 0; a <= m+10; a++) {
			for (int c = 0; c <= m+1; c++) {
				INFO("a = " << a);
				INFO("m = " << m);
				INFO("c = " << c);
				REQUIRE(mod_count(a, m, c, 100) == mod_count_slow(a, m, c, 100));
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
