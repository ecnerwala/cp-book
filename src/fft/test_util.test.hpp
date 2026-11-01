#pragma once

#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "fft/engines/ntt.hpp"
#include "fft/engines/real.hpp"
#include "fft/engines/split.hpp"
#include "fft/engines/crt.hpp"
#include "num/modnum.hpp"

// Shared helpers for the fft/ unit tests.

namespace wala::fft {

template <std::ranges::random_access_range A, std::ranges::random_access_range B>
	requires std::same_as<std::ranges::range_value_t<A>, std::ranges::range_value_t<B>>
vec<std::ranges::range_value_t<A>> multiply_slow(const A& a, const B& b) {
	using T = std::ranges::range_value_t<A>;
	if (std::ranges::empty(a) || std::ranges::empty(b)) return {};
	vec<T> res(int(std::ranges::size(a) + std::ranges::size(b)) - 1);
	for (int i = 0; i < int(std::ranges::size(a)); i++) {
		for (int j = 0; j < int(std::ranges::size(b)); j++) {
			res[i+j] += a[i] * b[j];
		}
	}
	return res;
}

// r's first n elements, zero-extended (with fill) or truncated to length n.
template <std::ranges::random_access_range R>
vec<std::ranges::range_value_t<R>> resized(const R& r, int n, std::ranges::range_value_t<R> fill) {
	vec<std::ranges::range_value_t<R>> res(n, fill);
	std::copy(std::ranges::begin(r), std::ranges::begin(r) + std::min(n, int(std::ranges::size(r))), res.begin());
	return res;
}

// Small values for doubles so products are exactly representable; full range otherwise.
template <typename T> T rnd_val(std::mt19937& mt) {
	if constexpr (std::is_floating_point_v<T>) return T(int(mt() % 1024));
	else return T(mt());
}
template <std::ranges::range V> void fill_rnd(V& v, std::mt19937& mt) {
	for (auto& x : v) x = rnd_val<std::ranges::range_value_t<V>>(mt);
}
template <std::ranges::random_access_range A, std::ranges::random_access_range B>
	requires std::same_as<std::ranges::range_value_t<A>, std::ranges::range_value_t<B>>
void check_eq(const A& got, const B& want) {
	using T = std::ranges::range_value_t<A>;
	if constexpr (std::is_floating_point_v<T>) {
		REQUIRE(int(std::ranges::size(got)) == int(std::ranges::size(want)));
		for (int i = 0; i < int(std::ranges::size(got)); i++) {
			CAPTURE(i);
			REQUIRE(llround(got[i]) == llround(want[i]));
		}
	} else {
		REQUIRE(std::ranges::equal(got, want));
	}
}

#define ALL_ENGINES \
		engines::ntt<modnum<998244353>>, engines::ntt<mod_goldilocks>, engines::real<double>, \
		engines::split<modnum<int(1e9)+7>>, engines::crt<modnum<int(1e9)+7>>
#define MOD_ENGINES \
		engines::ntt<modnum<998244353>>, engines::ntt<mod_goldilocks>, \
		engines::split<modnum<int(1e9)+7>>, engines::crt<modnum<int(1e9)+7>>

} // namespace wala::fft
