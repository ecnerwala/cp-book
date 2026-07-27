#pragma once

#include "fft/core.hpp"

namespace ecnerwala {
namespace fft {
namespace engines {

template <typename num> struct ntt {
	using value_type = num;
	static constexpr bool commutative = true;
	using core = fft_core<num>;
	struct transformed {
		vector<num> v;
		int size() const { return sz(v); }
	};
	using product = transformed;
	static constexpr int unit_scale = 0;
	template <int A = 0> using transformed_t = transformed;
	template <int K = 0> using product_t = product;

	static transformed transform(std::span<const num> a, int n) {
		assert(sz(a) <= 2 * n);
		transformed r;
		r.v.assign(n, num(0));
		int lo = min(sz(a), n);
		std::copy(a.begin(), a.begin() + lo, r.v.begin());
		for (int i = n; i < sz(a); i++) r.v[i - n] += a[i];
		core::forward(std::span<num>(r.v));
		return r;
	}
	static void extend_to(transformed& t, int m, std::span<const num> coeffs) {
		assert(!(m & (m-1)) && sz(coeffs) <= 2 * m);
		if (t.size() >= m) return;
		if (t.size() == 0) { t = transform(coeffs, m); return; }
		assert(sz(coeffs) <= 2 * t.size());
		while (t.size() < m) {
			t.v.resize(2 * t.size());
			core::extend(std::span<num>(t.v), coeffs);
		}
	}
	static transformed downsample(const transformed& t, int n, bool odd) {
		transformed r; r.v.resize(n);
		if (odd) core::odd_half(std::span<const num>(t.v), std::span<num>(r.v));
		else core::even_half(std::span<const num>(t.v), std::span<num>(r.v));
		return r;
	}
	static transformed negate_arg(const transformed& t, int n) {
		assert(n >= 2 && t.size() >= n);
		transformed r; r.v.resize(n);
		for (int j = 0; j < n; j++) r.v[j] = t.v[j ^ 1];
		return r;
	}
	static product mul(const transformed& a, const transformed& b, int n) {
		assert(a.size() >= n && b.size() >= n);
		product p; p.v.resize(n);
		for (int i = 0; i < n; i++) p.v[i] = a.v[i] * b.v[i];
		return p;
	}
	static product sq(const transformed& a, int n) { return mul(a, a, n); }
	static product mul2(
		const transformed& a1, const transformed& b1,
		const transformed& a2, const transformed& b2,
		int n
	) {
		assert(a1.size() >= n && b1.size() >= n && a2.size() >= n && b2.size() >= n);
		product p; p.v.resize(n);
		for (int i = 0; i < n; i++) p.v[i] = a1.v[i] * b1.v[i] + a2.v[i] * b2.v[i];
		return p;
	}
	static product add(product&& a, const product& b) {
		assert(a.size() == b.size());
		for (int i = 0; i < a.size(); i++) a.v[i] += b.v[i];
		return std::move(a);
	}
	template <typename Op = assign_op> static void finish(product&& p, std::span<num> out, Op op = {}) {
		int n = p.size();
		assert(sz(out) <= n);
		core::inverse(std::span<num>(p.v));
		num d = inv(num(n));
		for (int i = 0; i < sz(out); i++) op(out[i], p.v[i] * d);
	}
};

static_assert(engine<ntt<modnum<998244353>>>);
static_assert(engine<ntt<mod_goldilocks>>);

/* namespace engines */ }
/* namespace fft */ }
/* namespace ecnerwala */ }
