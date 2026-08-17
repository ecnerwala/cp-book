#pragma once

#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "fft/engines/ntt.hpp"
#include "fft/engines/real.hpp"
#include "fft/engines/split.hpp"
#include "fft/engines/crt.hpp"
#include "num/modnum.hpp"

// Shared helpers for the fft/ unit tests.

namespace wala {
namespace fft {

template <typename T> std::vector<T> multiply_slow(const std::vector<T>& a, const std::vector<T>& b) {
	if (a.empty() || b.empty()) return {};
	std::vector<T> res(a.size() + b.size() - 1);
	for (int i = 0; i < int(a.size()); i++) {
		for (int j = 0; j < int(b.size()); j++) {
			res[i+j] += a[i] * b[j];
		}
	}
	return res;
}

// Small values for doubles so products are exactly representable; full range otherwise.
template <typename T> T rnd_val(std::mt19937& mt) {
	if constexpr (std::is_floating_point_v<T>) return T(int(mt() % 1024));
	else return T(mt());
}
template <typename T> void fill_rnd(std::vector<T>& v, std::mt19937& mt) {
	for (T& x : v) x = rnd_val<T>(mt);
}
template <typename T> void check_eq(const std::vector<T>& got, const std::vector<T>& want) {
	REQUIRE(got.size() == want.size());
	for (int i = 0; i < int(got.size()); i++) {
		INFO("i = " << i);
		if constexpr (std::is_floating_point_v<T>) REQUIRE(llround(got[i]) == llround(want[i]));
		else REQUIRE(got[i] == want[i]);
	}
}

#define ALL_ENGINES \
		engines::ntt<modnum<998244353>>, engines::ntt<mod_goldilocks>, engines::real<double>, \
		engines::split<modnum<int(1e9)+7>>, engines::crt<modnum<int(1e9)+7>>
#define MOD_ENGINES \
		engines::ntt<modnum<998244353>>, engines::ntt<mod_goldilocks>, \
		engines::split<modnum<int(1e9)+7>>, engines::crt<modnum<int(1e9)+7>>

}} // namespace wala::fft
