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

TEST_CASE("poly_ap_values eval", "[fft,poly_ap_values]") {
	using num = modnum<998244353>;
	using poly_vals = poly_ap_values_fft<num>;
	mt19937 mt(48);
	for (int len : {0, 1, 2, 3, 5, 8, 13, 21}) {
		INFO("len = " << len);
		std::vector<num> poly(len);
		for (int i = 0; i < len; i++) poly[i] = num(mt());
		auto eval_at = [&](num v) {
			num r = 0;
			for (int i = len-1; i >= 0; i--) r = r * v + poly[i];
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
