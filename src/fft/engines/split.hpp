#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

#include "fft/core.hpp"
#include "fft/engine.hpp"

namespace ecnerwala::fft::engines {

// Multiplies mod `mnum` by splitting values into balanced 15-bit halves (each limb in
// [-2^14, 2^14], from the balanced representative |v| <= MOD/2) packed into one complex
// transform per operand.
template <typename mnum> struct split {
	static_assert(sizeof(decltype(mnum::MOD)) <= 4, "limbs must fit 15 bits");
	using value_type = mnum;
	static constexpr bool commutative = true;
	static constexpr int unit_scale = 1;
	using cnum = cplx<double>;
	using core = fft_core<cnum>;
	template <int A = 1> struct transformed_t {
		vector<cnum> v;
		int size() const { return sz(v); }
		transformed_t() = default;
		explicit transformed_t(vector<cnum>&& v_) : v(std::move(v_)) {}
		template <int A2> requires (A2 != A) explicit(A2 > A) transformed_t(transformed_t<A2>&& o)
			: v(std::move(o.v)) {}
	};
	using transformed = transformed_t<1>;
	template <int K> struct product_t {
		// After finish's inverse transforms: lo = (lo*lo, hi*lo), hi = (lo*hi, hi*hi).
		vector<cnum> lo, hi;
		int size() const { return sz(lo); }
		product_t() = default;
		product_t(vector<cnum>&& lo_, vector<cnum>&& hi_) : lo(std::move(lo_)), hi(std::move(hi_)) {}
		template <int K2> requires (K2 != K) explicit(K2 > K) product_t(product_t<K2>&& o)
			: lo(std::move(o.lo)), hi(std::move(o.hi)) {}
	};
	using product = product_t<1>;

	static cnum pack(mnum x) {
		int64_t v = x.balanced();
		int64_t hi = (v + (1 << 14)) >> 15;
		return cnum(double(v - (hi << 15)), double(hi));
	}

	static transformed transform(std::span<const mnum> a, int n) {
		assert(sz(a) <= 2 * n);
		transformed r;
		r.v.assign(n, cnum(0));
		for (int i = 0; i < sz(a); i++) {
			int j = i < n ? i : i - n;
			r.v[j] = r.v[j] + pack(a[i]);
		}
		core::forward(std::span<cnum>(r.v));
		return r;
	}
	static void extend_to(transformed& t, int m, std::span<const mnum> coeffs) {
		assert(!(m & (m-1)) && sz(coeffs) <= 2 * m);
		if (t.size() >= m) return;
		if (t.size() == 0) { t = transform(coeffs, m); return; }
		auto buf = buffer_pool<cnum>::get(sz(coeffs));
		for (int i = 0; i < sz(coeffs); i++) buf[i] = pack(coeffs[i]);
		while (t.size() < m) {
			int s = t.size();
			t.v.resize(2 * s);
			// coeffs past 2s are zero: they didn't fit in the transform we're a prefix of
			core::extend(
				std::span<cnum>(t.v),
				std::span<const cnum>(buf.span()).first(size_t(min(sz(coeffs), 2 * s)))
			);
		}
	}
	template <int A> static transformed_t<A> downsample(const transformed_t<A>& t, int n, bool odd) {
		transformed_t<A> r; r.v.resize(n);
		core::downsample(std::span<const cnum>(t.v), std::span<cnum>(r.v), odd);
		return r;
	}
	template <int K> static product_t<K> downsample(const product_t<K>& p, int n, bool odd) {
		product_t<K> r; r.lo.resize(n); r.hi.resize(n);
		core::downsample(std::span<const cnum>(p.lo), std::span<cnum>(r.lo), odd);
		core::downsample(std::span<const cnum>(p.hi), std::span<cnum>(r.hi), odd);
		return r;
	}
	template <int A> static transformed_t<A> upsample(const transformed_t<A>& t, int n, bool odd) {
		transformed_t<A> r; r.v.resize(n);
		core::upsample(std::span<const cnum>(t.v), std::span<cnum>(r.v), odd);
		return r;
	}
	template <int K> static product_t<K> upsample(const product_t<K>& p, int n, bool odd) {
		product_t<K> r; r.lo.resize(n); r.hi.resize(n);
		core::upsample(std::span<const cnum>(p.lo), std::span<cnum>(r.lo), odd);
		core::upsample(std::span<const cnum>(p.hi), std::span<cnum>(r.hi), odd);
		return r;
	}
	template <int A> static transformed_t<A> negate_arg(const transformed_t<A>& t, int n) {
		assert(n >= 2 && t.size() >= n);
		transformed_t<A> r; r.v.resize(n);
		for (int j = 0; j < n; j++) r.v[j] = t.v[j ^ 1];
		return r;
	}
	template <int A, int B> static transformed_t<A + B> add(transformed_t<A>&& a, const transformed_t<B>& b) {
		transformed_t<A + B> r{std::move(a.v)};
		add_into(r.v, b.v);
		return r;
	}
	// Unpacks b's transform into transforms of its low/high halves via conjugate
	// symmetry, then multiplies both against a's (still packed) transform. The scale
	// parameter only affects the bookkeeping, so the body is a shared untyped impl.
	static void mul_impl(const vector<cnum>& a, const vector<cnum>& b, vector<cnum>& lo, vector<cnum>& hi, int n, bool acc = false) {
		core::init(n);
		lo.resize(n); hi.resize(n);
		for (int i = 0; i < n; i++) {
			int ci = core::conj_index(i);
			cnum g0 = (b[i] + conj(b[ci])) * cnum(0.5);
			cnum t = (b[i] - conj(b[ci])) * cnum(0.5);
			cnum g1 = cnum(t.y, -t.x);
			if (acc) {
				lo[i] = lo[i] + a[i] * g0;
				hi[i] = hi[i] + a[i] * g1;
			} else {
				lo[i] = a[i] * g0;
				hi[i] = a[i] * g1;
			}
		}
	}
	template <int A, int B> static product_t<A * B> mul(const transformed_t<A>& a, const transformed_t<B>& b, int n) {
		assert(a.size() >= n && b.size() >= n);
		product_t<A * B> p;
		mul_impl(a.v, b.v, p.lo, p.hi, n);
		return p;
	}
	template <int A> static product_t<A * A> sq(const transformed_t<A>& a, int n) { return mul(a, a, n); }
	template <int A1, int B1, int A2, int B2>
	static product_t<A1 * B1 + A2 * B2> mul2(
		const transformed_t<A1>& a1, const transformed_t<B1>& b1,
		const transformed_t<A2>& a2, const transformed_t<B2>& b2,
		int n
	) {
		assert(a1.size() >= n && b1.size() >= n && a2.size() >= n && b2.size() >= n);
		product_t<A1 * B1 + A2 * B2> p;
		mul_impl(a1.v, b1.v, p.lo, p.hi, n);
		mul_impl(a2.v, b2.v, p.lo, p.hi, n, true);
		return p;
	}
	static void add_into(vector<cnum>& a, const vector<cnum>& b) {
		assert(sz(a) == sz(b));
		for (int i = 0; i < sz(a); i++) a[i] = a[i] + b[i];
	}
	template <int K1, int K2> static product_t<K1 + K2> add(product_t<K1>&& a, product_t<K2>&& b) {
		product_t<K1 + K2> r{std::move(a.lo), std::move(a.hi)};
		add_into(r.lo, b.lo);
		add_into(r.hi, b.hi);
		return r;
	}
	template <int K = 1, typename Op = assign_op> static void finish(product_t<K>&& p, std::span<mnum> out, Op op = {}) {
		// The fp error budget is divided by the accumulated scale; K <= 2 is very
		// conservative (balanced limbs already left ~2x headroom at max lengths).
		static_assert(K <= 2, "split: accumulated scale too large");
		int n = p.size();
		assert(sz(out) <= n);
		core::inverse(std::span<cnum>(p.lo));
		core::inverse(std::span<cnum>(p.hi));
		const int64_t m = mnum::MOD;
		double d = 1.0 / double(n);
		// llround + a final wrap so negative half-products (e.g. from negate_arg'd
		// transforms) reconstruct correctly.
		for (int i = 0; i < sz(out); i++) {
			int64_t v = (llround(p.lo[i].x * d)
					+ (llround(p.lo[i].y * d) % m << 15)
					+ (llround(p.hi[i].x * d) % m << 15)
					+ (llround(p.hi[i].y * d) % m << 30)) % m;
			if (v < 0) v += m;
			op(out[i], mnum(v));
		}
	}
};

/* namespace ecnerwala::fft::engines */ }
