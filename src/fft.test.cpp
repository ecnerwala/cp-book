#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "fft.hpp"

namespace ecnerwala {
namespace fft {

using namespace std;

template <typename T> vector<T> multiply_slow(const vector<T>& a, const vector<T>& b) {
	if (a.empty() || b.empty()) return {};
	vector<T> res(a.size() + b.size() - 1);
	for (int i = 0; i < int(a.size()); i++) {
		for (int j = 0; j < int(b.size()); j++) {
			res[i+j] += a[i] * b[j];
		}
	}
	return res;
}

TEST_CASE("FFT Multiply Mod", "[fft]") {
	using num = modnum<int(1e9)+7>;
	mt19937 mt(48);
	vector<num> a(100);
	vector<num> b(168);
	for (num& x : a) { x = num(mt()); }
	for (num& x : b) { x = num(mt()); }
	REQUIRE(multiply<fft_mod_multiplier<num>>(a,b) == multiply_slow(a, b));
}

TEST_CASE("FFT Inverse", "[fft]") {
	using num = modnum<998244353>;
	mt19937 mt(48);
	vector<num> a(298);
	for (num& x : a) { x = num(mt()); }
	auto i = inverse<multiply_inverser<fft_multiplier<num>, num>>(a);
	auto r = multiply<fft_multiplier<num>>(a, i);
	REQUIRE(r == multiply_slow(a, i));

	r.resize(a.size());
	vector<num> tgt(a.size());
	tgt[0] = 1;
	REQUIRE(r == tgt);
}

}} // namespace ecnerwala::fft
