#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_get_random_seed.hpp>

#include "fft/core.hpp"
#include "modnum.hpp"

namespace ecnerwala {
namespace fft {

using namespace std;

TEST_CASE("FFT core roundtrip", "[fft]") {
	using num = modnum<998244353>;
	mt19937 mt(Catch::getSeed());
	for (int n : {1, 2, 4, 8, 64, 256}) {
		vector<num> a(n);
		for (num& x : a) { x = num(mt()); }
		vector<num> t = a;
		fft_core<num>::forward(span<num>(t));
		fft_core<num>::inverse(span<num>(t));
		num d = inv(num(n));
		for (int i = 0; i < n; i++) REQUIRE(t[i] * d == a[i]);
	}
}

TEST_CASE("FFT core extend matches direct transform", "[fft]") {
	using num = modnum<998244353>;
	mt19937 mt(Catch::getSeed());
	for (int n : {1, 2, 8, 64}) {
		vector<num> a(n);
		for (num& x : a) { x = num(mt()); }
		// direct size-4n transform of a (zero-padded)
		vector<num> direct(4 * n, num(0));
		copy(a.begin(), a.end(), direct.begin());
		fft_core<num>::forward(span<num>(direct));
		// grown by doubling
		vector<num> t = a;
		fft_core<num>::forward(span<num>(t));
		t.resize(2 * n);
		fft_core<num>::extend(span<num>(t), span<const num>(a));
		t.resize(4 * n);
		fft_core<num>::extend(span<num>(t), span<const num>(a));
		REQUIRE(t == direct);
	}
}

TEST_CASE("FFT core even/odd halves", "[fft]") {
	using num = modnum<998244353>;
	mt19937 mt(Catch::getSeed());
	for (int n : {2, 8, 64}) {
		vector<num> a(2 * n);
		for (num& x : a) { x = num(mt()); }
		vector<num> evens(n), odds(n);
		for (int i = 0; i < n; i++) {
			evens[i] = a[2 * i];
			odds[i] = a[2 * i + 1];
		}
		vector<num> t = a;
		fft_core<num>::forward(span<num>(t));
		vector<num> te(n), to(n);
		fft_core<num>::even_half(span<const num>(t), span<num>(te));
		fft_core<num>::odd_half(span<const num>(t), span<num>(to));
		fft_core<num>::forward(span<num>(evens));
		fft_core<num>::forward(span<num>(odds));
		REQUIRE(te == evens);
		REQUIRE(to == odds);
	}
}

TEST_CASE("FFT transform prefix is transform mod x^n - 1", "[fft]") {
	using num = modnum<998244353>;
	mt19937 mt(Catch::getSeed());
	int n = 16;
	vector<num> a(2 * n);
	for (num& x : a) { x = num(mt()); }
	vector<num> t = a;
	fft_core<num>::forward(span<num>(t));
	vector<num> folded(n);
	for (int i = 0; i < n; i++) folded[i] = a[i] + a[i + n];
	fft_core<num>::forward(span<num>(folded));
	for (int i = 0; i < n; i++) REQUIRE(t[i] == folded[i]);
}

}
} // namespace ecnerwala::fft
