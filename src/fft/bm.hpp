#pragma once

#include <cassert>
#include <cstdint>
#include <ranges>

#include "containers/vec.hpp"

namespace wala {

template <std::ranges::random_access_range R> requires std::ranges::sized_range<R>
vec<std::ranges::range_value_t<R>> BerlekampMassey(const R& s) {
	using num = std::ranges::range_value_t<R>;
	int n = int(std::ranges::size(s)), L = 0, m = 0;
	vec<num> C(n), B(n), T;
	C[0] = B[0] = 1;

	num b = 1;
	for (int i = 0; i < n; i++) { ++m;
		num d = s[i];
		for (int j = 1; j <= L; j++) d += C[j] * s[i - j];
		if (d == 0) continue;
		T = C; num coef = d / b;
		for (int j = m; j < n; j++) C[j] -= coef * B[j - m];
		if (2 * L > i) continue;
		L = i + 1 - L; B = T; b = d; m = 0;
	}

	vec<num> res(L);
	for (int i = 0; i < L; i++) res[i] = -C[i + 1];
	return res;
}

template <std::ranges::random_access_range R, std::ranges::random_access_range Tr>
	requires std::ranges::sized_range<R> && std::ranges::sized_range<Tr>
		&& std::same_as<std::ranges::range_value_t<R>, std::ranges::range_value_t<Tr>>
std::ranges::range_value_t<Tr> linearRec(const R& S, const Tr& tr, int64_t k) {
	using num = std::ranges::range_value_t<Tr>;
	int n = int(std::ranges::size(tr));
	assert(int(std::ranges::size(S)) >= n);

	auto combine = [&](vec<num> a, vec<num> b, bool e = false) {
		// multiply a * b * x^e
		vec<num> res(a.size() + b.size());
		for (int i = 0; i < a.size(); i++) {
			for (int j = 0; j < b.size(); j++) {
				res[i + j + e] += a[i] * b[j];
			}
		}
		for (int i = res.size() - 1; i >= n; --i) {
			for (int j = 0; j < n; j++) {
				res[i - 1 - j] += res[i] * tr[j];
			}
		}
		vec<num> out(n);
		std::copy(res.begin(), res.begin() + n, out.begin());
		return out;
	};

	vec<num> pol(n);
	if (n > 0) pol[0] = num(1);

	assert(k >= 0);
	for (int i = 64 - 1 - (k == 0 ? 64 : __builtin_clzll(k)); i >= 0; i--) {
		pol = combine(pol, pol, (k >> i) & 1);
	}

	num res = 0;
	for (int i = 0; i < n; i++) res += pol[i] * S[i];
	return res;
}

} // namespace wala
