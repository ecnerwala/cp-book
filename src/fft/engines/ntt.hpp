#pragma once

#include <algorithm>
#include <cassert>
#include <span>
#include <utility>
#include <vector>

#include "fft/core.hpp"
#include "fft/engine.hpp"

namespace ecnerwala::fft::engines {

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

	// Out-parameter forms overwrite out, reusing its capacity.
	// out may alias an input to mul/downsample (they read forward of every write);
	// it must not alias the input to transform/negate_arg.
	static void transform(std::span<const num> a, int n, transformed& out) {
		assert(sz(a) <= 2 * n);
		out.v.assign(n, num(0));
		int lo = min(sz(a), n);
		std::copy(a.begin(), a.begin() + lo, out.v.begin());
		for (int i = n; i < sz(a); i++) out.v[i - n] += a[i];
		core::forward(std::span<num>(out.v));
	}
	static transformed transform(std::span<const num> a, int n) {
		transformed r;
		transform(a, n, r);
		return r;
	}
	static void extend_to(transformed& t, int m, std::span<const num> coeffs) {
		assert(!(m & (m-1)) && sz(coeffs) <= 2 * m);
		if (t.size() >= m) return;
		if (t.size() == 0) { transform(coeffs, m, t); return; }
		while (t.size() < m) {
			int s = t.size();
			t.v.resize(2 * s);
			// coeffs past 2s are zero: they didn't fit in the transform we're a prefix of
			core::extend(std::span<num>(t.v), coeffs.first(size_t(min(sz(coeffs), 2 * s))));
		}
	}
	static void downsample(const transformed& t, int n, bool odd, transformed& out) {
		auto in = std::span<const num>(t.v);
		if (&out != &t) out.v.resize(n);
		auto o = std::span<num>(out.v).first(size_t(n));
		if (odd) core::odd_half(in, o);
		else core::even_half(in, o);
		if (&out == &t) out.v.resize(n);
	}
	static transformed downsample(const transformed& t, int n, bool odd) {
		transformed r;
		downsample(t, n, odd, r);
		return r;
	}
	static void negate_arg(const transformed& t, int n, transformed& out) {
		assert(n >= 2 && t.size() >= n && &out != &t);
		out.v.resize(n);
		for (int j = 0; j < n; j++) out.v[j] = t.v[j ^ 1];
	}
	static transformed negate_arg(const transformed& t, int n) {
		transformed r;
		negate_arg(t, n, r);
		return r;
	}
	static void mul(const transformed& a, const transformed& b, int n, product& out) {
		assert(a.size() >= n && b.size() >= n);
		out.v.resize(n);
		for (int i = 0; i < n; i++) out.v[i] = a.v[i] * b.v[i];
	}
	static product mul(const transformed& a, const transformed& b, int n) {
		product p;
		mul(a, b, n, p);
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

/* namespace ecnerwala::fft::engines */ }
