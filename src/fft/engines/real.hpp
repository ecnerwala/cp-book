#pragma once

#include <algorithm>
#include <cassert>
#include <span>
#include <utility>
#include <vector>

#include "fft/core.hpp"
#include "fft/engine.hpp"

namespace wala::fft::engines {

// Convolve real (floating point) values by packing into complex numbers with
//   a'[t] = a[2t] + i * a[2t+1]
// We use conjugate symmetry to untangle/retangle the two.
// TODO: Add type bounds?
template <typename dbl = double> struct real {
	using value_type = dbl;
	static constexpr bool commutative = true;
	using cnum = cplx<dbl>;
	using core = fft_core<cnum>;
	struct transformed {
		vector<cnum> v;
		int size() const { return 2 * sz(v); }
	};
	using product = transformed;
	// Precision is caller-managed for this engine (see add), so scale is untracked.
	static constexpr int unit_scale = 0;
	template <int A = 0> using transformed_t = transformed;
	template <int K = 0> using product_t = product;

	static int packed_size(int n) { return std::max(n / 2, 1); }
	static void pack(std::span<const dbl> a, std::span<cnum> c) {
		for (int i = 0; i < sz(a); i++) (i & 1 ? c[i/2].y : c[i/2].x) = a[i];
	}
	// Spectrum of the real (odd = false) or imaginary (odd = true) part of the packed
	// sequence at bitrev entry t, by conjugate symmetry with the entry of w^{-k}.
	static cnum part(const transformed& f, int t, bool odd) {
		cnum g = conj(f.v[core::conj_index(t)]);
		return odd ? (f.v[t] - g) * cnum(0, dbl(-0.5)) : (f.v[t] + g) * cnum(dbl(0.5));
	}
	// Given the spectra (s0, s1) of a real sequence x at w_{2mo}^q and w_{2mo}^{q+mo},
	// the packed-transform entry of x at packed size mo: the even/odd interleaves of x
	// have spectra (s0 +- s1)/2 (the odd one twisted by w_{2mo}^{-q}).
	static cnum retangle(cnum s0, cnum s1, int mo, int q) {
		cnum s = (s0 + s1) * cnum(dbl(0.5));
		cnum d = (s0 - s1) * cnum(dbl(0.5)) * core::inv_rt[mo + q];
		return s + cnum(-d.y, d.x);
	}

	static transformed transform(std::span<const dbl> a, int n) {
		assert(sz(a) <= 2 * n);
		transformed r;
		r.v.assign(packed_size(n), cnum(0));
		for (int i = 0; i < sz(a); i++) {
			int j = i < n ? i : i - n;
			((j & 1) ? r.v[j/2].y : r.v[j/2].x) += a[i];
		}
		core::forward(std::span<cnum>(r.v));
		return r;
	}
	static void extend_to(transformed& t, int m, std::span<const dbl> coeffs) {
		assert(!(m & (m-1)) && sz(coeffs) <= 2 * m);
		if (t.size() >= m) return;
		if (t.size() == 0) { t = transform(coeffs, m); return; }
		auto buf = buffer_pool<cnum>::get((sz(coeffs) + 1) / 2);
		std::fill(buf.span().begin(), buf.span().end(), cnum(0));
		pack(coeffs, buf.span());
		while (t.size() < m) {
			int s = sz(t.v);
			t.v.resize(2 * s);
			// packed coeffs past 2s are zero: they didn't fit in the transform we're a prefix of
			core::extend(
				std::span<cnum>(t.v),
				std::span<const cnum>(buf.span()).first(size_t(min(sz(buf.span()), 2 * s)))
			);
		}
	}
	// The complex packing pairs adjacent slots, so the packed buffer keeps the row
	// structure at stride C/2; C == 2 packs a whole row into one slot, so fall back.
	static transformed transform_padded(std::span<const dbl> a, int n, int C) {
		assert(sz(a) == n);
		if (C < 4) return transform(a, n);
		transformed r;
		r.v.assign(packed_size(n), cnum(0));
		pack(a, std::span<cnum>(r.v));
		core::forward_padded(std::span<cnum>(r.v), C / 2);
		return r;
	}
	static transformed downsample(const transformed& t, int n, bool odd) { return half(t, n, odd); }
	static transformed upsample(const transformed& t, int n, bool odd) { return unhalf(t, n, odd); }
	// A(-x) negates the odd (imaginary-slot) coefficients, i.e. conjugates the packed
	// sequence; the transform of a conjugated sequence is the conjugate at w^(-k).
	static transformed negate_arg(const transformed& t, int n) {
		int m = packed_size(n);
		assert(n >= 2 && sz(t.v) >= m);
		transformed r; r.v.resize(m);
		for (int j = 0; j < m; j++) r.v[j] = conj(t.v[core::conj_index(j)]);
		return r;
	}
	static transformed half(const transformed& f, int n, bool odd) {
		assert(n >= 2 && f.size() >= 2 * n);
		int mo = n / 2;
		core::init(2 * mo);
		transformed r; r.v.resize(mo);
		for (int u = 0; u < mo; u++) {
			r.v[u] = retangle(part(f, 2*u, odd), part(f, 2*u+1, odd), mo, core::brev(u, mo));
		}
		return r;
	}
	// Inverse of half: the size-n packed transform of the input spread as evens or odds.
	static transformed unhalf(const transformed& f, int n, bool odd) {
		assert(n >= 2 && 2 * f.size() >= n);
		int mo = packed_size(n);
		transformed r; r.v.resize(mo);
		if (mo == 1) {
			r.v[0] = part(f, 0, false);
		} else {
			// The spread-as-evens sequence packs to the plain (purely real) input sequence,
			// so untangle its packed transform into the full spectrum:
			// A(w) = E(w^2) + w O(w^2) at w = +-w_{2mi}^brev(u, mi).
			int mi = mo / 2;
			core::init(2 * mi);
			for (int u = 0; u < mi; u++) {
				cnum e = part(f, u, false);
				cnum o = part(f, u, true) * core::rt[mi + core::brev(u, mi)];
				r.v[2*u+0] = e + o;
				r.v[2*u+1] = e - o;
			}
		}
		// spreading as odds instead packs into the imaginary slots, scaling the spectrum by i
		if (odd) for (auto& z : r.v) z = cnum(-z.y, z.x);
		return r;
	}
	// One Graeffe step: the size-n/2 packed product of the even part of A(x) * A(-x).
	// The product takes the same value A(w) A(-w) at each pair +-w, so each value is
	// already a spectrum value of the halved (real) result; repack adjacent pairs.
	static product graeffe(const transformed& t, int n) {
		int m = packed_size(n);
		assert(n >= 2 && t.size() >= n);
		core::init(2 * m);
		auto val = [&](int j) {
			cnum x = part(t, j, false), y = part(t, j, true) * core::rt[m + core::brev(j, m)];
			return (x + y) * (x - y);
		};
		product p;
		if (m == 1) {
			p.v = {val(0)};
		} else {
			p.v.resize(m / 2);
			for (int u = 0; u < m/2; u++) {
				p.v[u] = retangle(val(2*u), val(2*u+1), m/2, core::brev(u, m/2));
			}
		}
		return p;
	}
	static product mul(const transformed& a, const transformed& b, int n) {
		int m = packed_size(n);
		assert(a.size() >= n && b.size() >= n);
		core::init(2 * m);
		product p; p.v.resize(m);
		for (int t = 0; t < m; t++) {
			int k = core::brev(t, m);
			cnum w = core::rt[m + k];
			cnum xa = part(a, t, false), ya = part(a, t, true);
			cnum xb = part(b, t, false), yb = part(b, t, true);
			// full spectra at w_{2m}^k and w_{2m}^{k+m} = -w_{2m}^k
			cnum p0 = (xa + w * ya) * (xb + w * yb);
			cnum p1 = (xa - w * ya) * (xb - w * yb);
			p.v[t] = retangle(p0, p1, m, k);
		}
		return p;
	}
	static product sq(const transformed& a, int n) { return mul(a, a, n); }
	static product mul2(
		const transformed& a1, const transformed& b1,
		const transformed& a2, const transformed& b2,
		int n
	) {
		int m = packed_size(n);
		assert(a1.size() >= n && b1.size() >= n && a2.size() >= n && b2.size() >= n);
		core::init(2 * m);
		product p; p.v.resize(m);
		for (int t = 0; t < m; t++) {
			int k = core::brev(t, m);
			cnum w = core::rt[m + k];
			cnum xa1 = part(a1, t, false), ya1 = part(a1, t, true);
			cnum xb1 = part(b1, t, false), yb1 = part(b1, t, true);
			cnum xa2 = part(a2, t, false), ya2 = part(a2, t, true);
			cnum xb2 = part(b2, t, false), yb2 = part(b2, t, true);
			cnum p0 = (xa1 + w * ya1) * (xb1 + w * yb1) + (xa2 + w * ya2) * (xb2 + w * yb2);
			cnum p1 = (xa1 - w * ya1) * (xb1 - w * yb1) + (xa2 - w * ya2) * (xb2 - w * yb2);
			p.v[t] = retangle(p0, p1, m, k);
		}
		return p;
	}
	static product add(product&& a, const product& b) {
		assert(a.size() == b.size());
		for (int i = 0; i < sz(a.v); i++) a.v[i] = a.v[i] + b.v[i];
		return std::move(a);
	}
	template <typename Op = assign_op> static void finish(product&& p, std::span<dbl> out, Op op = {}) {
		int m = sz(p.v);
		assert(sz(out) <= 2 * m);
		core::inverse(std::span<cnum>(p.v));
		dbl d = dbl(1) / dbl(m);
		for (int i = 0; i < sz(out); i++) op(out[i], (i & 1 ? p.v[i/2].y : p.v[i/2].x) * d);
	}
	template <typename Op = assign_op> static void finish_padded(product&& p, std::span<dbl> out, int C, Op op = {}) {
		int m = sz(p.v);
		assert(sz(out) == 2 * m);
		if (C < 4) { finish(std::move(p), out, op); return; }
		core::inverse_padded(std::span<cnum>(p.v), C / 2);
		dbl d = dbl(1) / dbl(m);
		for (int i = 0; i < 2 * m; i += C) {
			for (int j = 0; j < C/2; j++) op(out[i+j], ((i+j) & 1 ? p.v[(i+j)/2].y : p.v[(i+j)/2].x) * d);
		}
	}
};

/* namespace wala::fft::engines */ }
