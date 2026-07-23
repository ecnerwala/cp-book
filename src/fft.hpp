#pragma once

#include <vector>
#include <cmath>
#include <utility>
#include <algorithm>
#include <span>
#include <cstdint>
#include <cassert>
#include <concepts>
#include <memory>
#include <new>

#include "modnum.hpp"

/**
 * Author: Andrew He
 * Source: http://neerc.ifmo.ru/trains/toulouse/2017/fft2.pdf
 * Papers about accuracy: http://www.daemonology.net/papers/fft.pdf, http://www.cs.berkeley.edu/~fateman/papers/fftvsothers.pdf
 * For integers rounding works if $(|a| + |b|)\max(a, b) < \mathtt{\sim} 10^9$, or in theory maybe $10^6$.
 *
 * Layering:
 *   fft_core<num>   root tables + in-place DIF/DIT transforms in bit-reversed order,
 *                   plus transform doubling (extend) and even/odd downsampling.
 *   conv engines    one type per (ring, strategy): fft_engine, fft_real_engine,
 *                   fft_split_engine, crt_engine. Each exposes a transform-domain
 *                   `transformed` type and a product-domain accumulator `product`,
 *                   so caching, doubling, and += accumulation work uniformly.
 *   conv layer      span-based multiply / square / multiply_circular / middle_product,
 *                   with the 2^k+1 wraparound optimization, plus fft_cache: the
 *                   engine-level cached operand (transform + logical length), grown
 *                   by half-cost doubling only when the owner of the coefficients
 *                   feeds them back through extend_to. Cached entry points take each
 *                   operand as a (coefficients, fft_cache&) pair.
 *   value types     power_series<E, exact> (exact = known in full vs truncated mod
 *                   x^len, natural coefficient order; mixed-exactness operators;
 *                   power_series_exact / power_series_trunc aliases), poly<E>
 *                   (exact, stored reversed as the power_series_exact rev(p);
 *                   rev_series / from_rev_series), linear_form<E> (same reversed
 *                   convention), the coefficient-owning cached wrappers
 *                   (prefix_cached_power_series, cached_power_series_exact),
 *                   poly_ap_values<E>, algorithms (inverse, log, exp, pow, compose,
 *                   multipoint evaluate/interpolate), and online (relaxed)
 *                   multiplication.
 *
 * Transform convention: forward is DIF (natural input -> bit-reversed output), inverse
 * is DIT (bit-reversed input -> natural output), unscaled; entry j of a size-n
 * transform is P evaluated at w_n^brev(j). Consequences used throughout:
 *   - The first n entries of a size-2n transform of P are the size-n transform of
 *     P mod (x^n - 1); for deg P < n this is the size-n transform of P, so transforms
 *     can be grown by doubling at half cost (extend) and any power-of-two prefix of a
 *     cached transform is directly usable.
 *   - Entries (2j, 2j+1) are evaluations at (w, -w), so downsampling to the transform
 *     of the even/odd part is a linear pass (even_half / odd_half).
 */

namespace ecnerwala {

template<class T> int sz(T&& arg) { using std::size; return int(size(std::forward<T>(arg))); }
inline int nextPow2(int s) { return 1 << (s > 1 ? 32 - __builtin_clz(s-1) : 0); }

namespace fft {

using std::swap;
using std::vector;
using std::min;
using std::max;

// Complex
template <typename dbl> struct cplx { /// start-hash
	dbl x, y;
	cplx(dbl x_ = 0, dbl y_ = 0) : x(x_), y(y_) { }
	friend cplx operator+(cplx a, cplx b) { return cplx(a.x + b.x, a.y + b.y); }
	friend cplx operator-(cplx a, cplx b) { return cplx(a.x - b.x, a.y - b.y); }
	friend cplx operator*(cplx a, cplx b) { return cplx(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }
	friend cplx conj(cplx a) { return cplx(a.x, -a.y); }
	friend cplx inv(cplx a) { dbl n = (a.x*a.x+a.y*a.y); return cplx(a.x/n,-a.y/n); }
};

// getRoot implementations
template <typename num> struct getRoot {
	static num f(int k) = delete;
};
template <typename dbl> struct getRoot<cplx<dbl>> {
	static cplx<dbl> f(int k) {
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
		dbl a=2*M_PI/k;
		return cplx<dbl>(cos(a),sin(a));
	}
};
template <int MOD> struct primitive_root {
	static const int value;
};
// 998244353 = (119 << 23) + 1 = 2^30 - 2^26 - 2^23 + 1
template <> struct primitive_root<998244353> {
	static const int value = 3;
};
// babybear prime
template <> struct primitive_root<(15 << 27) + 1> {
	static const int value = 31;
};
// koalabear prime
template <> struct primitive_root<(127 << 24) + 1> {
	static const int value = 3;
};
template <> struct primitive_root<(7 << 26) + 1> {
	static const int value = 3;
};
template <> struct primitive_root<(5 << 25) + 1> {
	static const int value = 3;
};
template <int MOD> struct getRoot<modnum<MOD>> {
	static modnum<MOD> f(int k) {
		assert((MOD-1)%k == 0);
		return power(modnum<MOD>(primitive_root<MOD>::value), (MOD-1)/k);
	}
};
template <> struct getRoot<mod_goldilocks> {
	static mod_goldilocks f(int k) {
		assert((mod_goldilocks::MOD-1)%k == 0);
		return power(mod_goldilocks(mod_goldilocks::PRIMITIVE_ROOT), (mod_goldilocks::MOD-1)/k);
	}
};

// Reusable scratch buffers. Not thread-safe by default: this is deliberately plain
// static storage so single-threaded programs pay no TLS indirection; define
// ECNERWALA_FFT_POOL_STORAGE to `thread_local` for multithreaded use.
#ifndef ECNERWALA_FFT_POOL_STORAGE
#define ECNERWALA_FFT_POOL_STORAGE
#endif
template <typename T> struct buffer_pool {
	static inline ECNERWALA_FFT_POOL_STORAGE std::vector<std::vector<T>> free_list;
	struct handle {
		std::vector<T> v;
		explicit handle(int n) {
			if (!free_list.empty()) {
				v = std::move(free_list.back());
				free_list.pop_back();
			}
			v.assign(n, T());
		}
		handle(const handle&) = delete;
		handle& operator=(const handle&) = delete;
		handle(handle&& o) noexcept : v(std::move(o.v)) {}
		~handle() { if (v.capacity()) free_list.push_back(std::move(v)); }
		T& operator[](int i) { return v[i]; }
		operator std::span<T>() { return std::span<T>(v); }
		std::span<T> span() { return std::span<T>(v); }
	};
	static handle get(int n) { return handle(n); }
};

template <typename num> struct fft_core {
	static inline vector<int> rev;
	static inline vector<num> rt, inv_rt;

	static void init(int n) {
		if (n <= sz(rt)) return;
		rev.resize(n);
		for (int i = 0; i < n; i++) {
			rev[i] = (rev[i>>1] | ((i&1)*n)) >> 1;
		}
		rt.reserve(n); inv_rt.reserve(n);
		while (sz(rt) < 2 && sz(rt) < n) { rt.push_back(num(1)); inv_rt.push_back(num(1)); }
		for (int k = sz(rt); k < n; k *= 2) {
			rt.resize(2*k); inv_rt.resize(2*k);
			num z = getRoot<num>::f(2*k);
			num iz = inv(z);
			for (int i = k/2; i < k; i++) {
				rt[2*i] = rt[i], rt[2*i+1] = rt[i]*z;
				inv_rt[2*i] = inv_rt[i], inv_rt[2*i+1] = inv_rt[i]*iz;
			}
		}
	}

	// Bit-reverse of i within a size-n transform (init(n) must have been called).
	static int brev(int i, int n) {
		int s = __builtin_ctz(unsigned(sz(rev)/n));
		return rev[i] >> s;
	}
	// Index of the evaluation point w^{-brev(j)}, i.e. the conjugate point of entry j.
	// In bitrev order, block [h, 2h) holds the odd multiples of w_{2h} with reversed
	// low bits, so negating the exponent reverses each block: XOR with h - 1.
	static int conj_index(int j) {
		return j == 0 ? 0 : j ^ ((1 << (31 - __builtin_clz(unsigned(j)))) - 1);
	}

	// Natural input -> bit-reversed output (DIF).
	static void forward(std::span<num> a) {
		int n = sz(a);
		if (n <= 1) return;
		init(n);
		for (int k = n/2; k >= 1; k /= 2) {
			for (int i = 0; i < n; i += 2*k) {
				for (int j = 0; j < k; j++) {
					num u = a[i+j], v = a[i+j+k];
					a[i+j] = u + v;
					a[i+j+k] = (u - v) * rt[j+k];
				}
			}
		}
	}

	// Bit-reversed input -> natural output, unscaled (result is n times the
	// inverse DFT).
	static void inverse(std::span<num> a) {
		int n = sz(a);
		if (n <= 1) return;
		init(n);
		for (int k = 1; k < n; k *= 2) {
			for (int i = 0; i < n; i += 2*k) {
				for (int j = 0; j < k; j++) {
					num t = inv_rt[j+k] * a[i+j+k];
					a[i+j+k] = a[i+j] - t;
					a[i+j] = a[i+j] + t;
				}
			}
		}
	}

	// t has size 2n; t[0..n) is the size-n transform of coeffs mod x^n - 1
	// (sz(coeffs) <= 2n). Fills t[n..2n) so t becomes the size-2n transform, at the
	// cost of one size-n forward transform (half the cost of recomputing from
	// scratch): the top half is the transform of the twisted fold (c_i - c_{i+n}) *
	// w_{2n}^i, since w_{2n}^n = -1.
	static void extend(std::span<num> t, std::span<const num> coeffs) {
		int n = sz(t) / 2;
		assert(sz(coeffs) <= 2 * n);
		init(sz(t));
		auto b = t.subspan(n, n);
		int lo = min(sz(coeffs), n);
		for (int i = 0; i < lo; i++) {
			// rt[n + i] = w_{2n}^i
			b[i] = coeffs[i] * rt[n + i];
		}
		std::fill(b.begin() + lo, b.end(), num(0));
		for (int i = n; i < sz(coeffs); i++) {
			b[i - n] = b[i - n] - coeffs[i] * rt[i];
		}
		forward(b);
	}

	// From the size-2n transform of P, the size-n transform of the even part
	// (E where P(x) = E(x^2) + x*O(x^2)): a linear pass, no FFT.
	static void even_half(std::span<const num> t, std::span<num> out) {
		int n = sz(out);
		assert(sz(t) >= 2*n);
		num half = inv(num(2));
		for (int j = 0; j < n; j++) out[j] = (t[2*j] + t[2*j+1]) * half;
	}
	// Size-n transform of the odd part O.
	static void odd_half(std::span<const num> t, std::span<num> out) {
		int n = sz(out);
		assert(sz(t) >= 2*n);
		init(2*n);
		num half = inv(num(2));
		for (int j = 0; j < n; j++) {
			// entry j of the size-2n transform pairs (w, -w) with w = w_{2n}^{brev(j, n)}
			out[j] = (t[2*j] - t[2*j+1]) * half * inv_rt[n + brev(j, n)];
		}
	}
};

// Output operations for the finish step: the inverse transform necessarily ends with a
// scaling pass, so arbitrary elementwise output ops fuse in for free.
struct assign_op { template <typename T> void operator()(T& d, T v) const { d = v; } };
struct add_op { template <typename T> void operator()(T& d, T v) const { d += v; } };
struct sub_op { template <typename T> void operator()(T& d, T v) const { d -= v; } };
struct add_twice_op { template <typename T> void operator()(T& d, T v) const { d += v + v; } };

// Conv engine contract. Semantics that can't be expressed in the concept:
//   transformed    operand-domain transform buffer; any power-of-two prefix of a
//                  size-m transform is a valid size-n transform (n <= m) of the
//                  coefficients mod (x^n - 1)
//   product        product-domain accumulator (may differ from transformed)
//   transform      size-n transform of a mod x^n - 1: zero-padded, and a tail beyond
//                  n is folded circularly while packing (no extra pass); requires
//                  sz(a) <= 2n. Lets 2^k+1-length operands transform at size 2^k
//                  (the 2^k+1 cut) without any caller-side folding
//   extend_to      grow a transform to size n by half-cost doubling; requires the
//                  original coefficients and sz(coeffs) <= 2 * current size (a
//                  current transform of the coefficients mod x^size - 1 is a valid
//                  starting prefix; the doubling folds the tail)
//   even/odd_half  size-n transform of the even/odd part of a size-2n transform,
//                  in linear time (no FFT)
//   negate_arg     size-n transform of A(-x) from a size-n transform of A, in linear
//                  time (adjacent bitrev entries are evaluations at (w, -w)); n >= 2
//   mul / sq       pointwise product of size-n prefixes
//   add            pointwise sum of two products, so k products share one inverse
//                  transform. Soundness bounds shrink per addend (the crt
//                  reconstruction range and the split engine's fp error budget are
//                  divided by the addend count), so the products of the inexact
//                  engines carry their accumulated term count as a compile-time
//                  parameter (product_t<K>; mul/sq give K = 1, add gives K1 + K2)
//                  and finish static_asserts a conservative K <= 2; exact engines
//                  are unbounded. Zero runtime cost. `product` is the K = 1 type.
//   finish         inverse transform + scale, then out[i] op= result[i] for
//                  i < sz(out); requires sz(out) <= size of the product
template <typename E>
concept conv_engine = requires(
		std::span<const typename E::value_type> in,
		std::span<typename E::value_type> out,
		typename E::transformed& t,
		const typename E::transformed& ct,
		typename E::product& p,
		int n) {
	typename E::value_type;
	{ E::transform(in, n) } -> std::same_as<typename E::transformed>;
	{ ct.size() } -> std::same_as<int>;
	E::extend_to(t, n, in);
	{ E::even_half(ct, n) } -> std::same_as<typename E::transformed>;
	{ E::odd_half(ct, n) } -> std::same_as<typename E::transformed>;
	{ E::negate_arg(ct, n) } -> std::same_as<typename E::transformed>;
	{ E::mul(ct, ct, n) } -> std::same_as<typename E::product>;
	{ E::sq(ct, n) } -> std::same_as<typename E::product>;
	E::finish(std::move(p), out);
	E::finish(std::move(p), out, add_op{});
	E::finish(E::add(std::move(p), std::move(p)), out);
};

template <typename num> struct fft_engine {
	using value_type = num;
	using core = fft_core<num>;
	struct transformed {
		vector<num> v;
		int size() const { return sz(v); }
	};
	// A pointwise product is itself a valid transform (of a*b mod x^n - 1), so
	// products can be fed back into mul/even_half/etc.
	using product = transformed;

	// Folds a mod x^n - 1 while copying (see the concept preamble).
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
	static void extend_to(transformed& t, int n, std::span<const num> coeffs) {
		assert(sz(coeffs) <= 2 * t.size());
		while (t.size() < n) {
			t.v.resize(2 * t.size());
			core::extend(std::span<num>(t.v), coeffs);
		}
	}
	static transformed even_half(const transformed& t, int n) {
		transformed r; r.v.resize(n);
		core::even_half(std::span<const num>(t.v), std::span<num>(r.v));
		return r;
	}
	static transformed odd_half(const transformed& t, int n) {
		transformed r; r.v.resize(n);
		core::odd_half(std::span<const num>(t.v), std::span<num>(r.v));
		return r;
	}
	// Transform of A(-x): adjacent bitrev entries are evaluations at (w, -w), so swap.
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
	// Exact ring: any number of addends is sound, and the sum is again a transform.
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

// Real convolution with a packed transform: coefficients a[2t], a[2t+1] share one
// complex point, so the logical size-n transform is a size-n/2 complex FFT (halving
// transform cost and cache memory). mul/sq untangle the two real spectra on the fly
// via conjugate symmetry and re-tangle the (real) product's spectrum, all pointwise.
// transformed::size() reports the logical (real) size.
template <typename dbl = double> struct fft_real_engine {
	using value_type = dbl;
	using cnum = cplx<dbl>;
	using core = fft_core<cnum>;
	struct transformed {
		vector<cnum> v;
		int size() const { return 2 * sz(v); }
	};
	using product = transformed;

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

	// Folds a mod x^n - 1 while packing (see the concept preamble): n is even for
	// n >= 2, so wrapped indices keep their parity and land in the matching slot;
	// for n == 1 everything folds into the real slot of the single point.
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
	static void extend_to(transformed& t, int n, std::span<const dbl> coeffs) {
		assert(sz(coeffs) <= 2 * t.size());
		auto buf = buffer_pool<cnum>::get((sz(coeffs) + 1) / 2);
		std::fill(buf.span().begin(), buf.span().end(), cnum(0));
		pack(coeffs, buf.span());
		while (t.size() < n) {
			t.v.resize(2 * sz(t.v));
			core::extend(std::span<cnum>(t.v), std::span<const cnum>(buf.span()));
		}
	}
	static transformed even_half(const transformed& t, int n) { return half(t, n, false); }
	static transformed odd_half(const transformed& t, int n) { return half(t, n, true); }
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
	// Pointwise sum of (re-tangled, real) product spectra. Precision is caller-managed
	// for this engine anyway; each addend adds its magnitude to the fp error budget.
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
};

// Multiplies mod `mnum` by splitting values into balanced 15-bit halves (each limb in
// [-2^14, 2^14], from the balanced representative |v| <= MOD/2) packed into one complex
// transform per operand (so per-operand transforms remain cacheable). Balancing halves
// each limb's magnitude, doubling the fp-precision headroom per operand.
template <typename mnum> struct fft_split_engine {
	using value_type = mnum;
	using cnum = cplx<double>;
	using core = fft_core<cnum>;
	struct transformed {
		vector<cnum> v;
		int size() const { return sz(v); }
	};
	// K = number of accumulated term products (see the conv_engine preamble).
	template <int K> struct product_t {
		// After finish's inverse transforms: lo = (lo*lo, hi*lo), hi = (lo*hi, hi*hi).
		vector<cnum> lo, hi;
		int size() const { return sz(lo); }
	};
	using product = product_t<1>;

	static cnum pack(mnum x) {
		int64_t v = int64_t(int(x));
		if (2 * v > int64_t(mnum::MOD)) v -= mnum::MOD;
		int64_t hi = (v + (1 << 14)) >> 15;
		return cnum(double(v - (hi << 15)), double(hi));
	}

	// Folds a mod x^n - 1 while packing (see the concept preamble); a folded
	// coefficient's limbs just add, staying comfortably in the fp error budget.
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
	static void extend_to(transformed& t, int n, std::span<const mnum> coeffs) {
		assert(sz(coeffs) <= 2 * t.size());
		auto buf = buffer_pool<cnum>::get(sz(coeffs));
		for (int i = 0; i < sz(coeffs); i++) buf[i] = pack(coeffs[i]);
		while (t.size() < n) {
			t.v.resize(2 * t.size());
			core::extend(std::span<cnum>(t.v), std::span<const cnum>(buf.span()));
		}
	}
	static transformed even_half(const transformed& t, int n) {
		transformed r; r.v.resize(n);
		core::even_half(std::span<const cnum>(t.v), std::span<cnum>(r.v));
		return r;
	}
	static transformed odd_half(const transformed& t, int n) {
		transformed r; r.v.resize(n);
		core::odd_half(std::span<const cnum>(t.v), std::span<cnum>(r.v));
		return r;
	}
	// The packed complex sequence's halves stay real (some entries just go negative,
	// which finish's signed reconstruction handles), so the plain complex-transform
	// identities apply: A(-x) swaps the (w, -w) bitrev pairs.
	static transformed negate_arg(const transformed& t, int n) {
		assert(n >= 2 && t.size() >= n);
		transformed r; r.v.resize(n);
		for (int j = 0; j < n; j++) r.v[j] = t.v[j ^ 1];
		return r;
	}
	// Unpacks b's transform into transforms of its low/high halves via conjugate
	// symmetry, then multiplies both against a's (still packed) transform.
	static product mul(const transformed& a, const transformed& b, int n) {
		assert(a.size() >= n && b.size() >= n);
		core::init(n);
		product p;
		p.lo.resize(n); p.hi.resize(n);
		for (int i = 0; i < n; i++) {
			int ci = core::conj_index(i);
			cnum g0 = (b.v[i] + conj(b.v[ci])) * cnum(0.5);
			cnum t = (b.v[i] - conj(b.v[ci])) * cnum(0.5);
			cnum g1 = cnum(t.y, -t.x);
			p.lo[i] = a.v[i] * g0;
			p.hi[i] = a.v[i] * g1;
		}
		return p;
	}
	static product sq(const transformed& a, int n) { return mul(a, a, n); }
	static void add_into(vector<cnum>& a, const vector<cnum>& b) {
		assert(sz(a) == sz(b));
		for (int i = 0; i < sz(a); i++) a[i] = a[i] + b[i];
	}
	template <int K1, int K2> static product_t<K1+K2> add(product_t<K1>&& a, product_t<K2>&& b) {
		product_t<K1+K2> r{std::move(a.lo), std::move(a.hi)};
		add_into(r.lo, b.lo);
		add_into(r.hi, b.hi);
		return r;
	}
	template <int K = 1, typename Op = assign_op> static void finish(product_t<K>&& p, std::span<mnum> out, Op op = {}) {
		// The fp error budget is divided by the addend count; K <= 2 is very
		// conservative (balanced limbs already left ~2x headroom at max lengths).
		static_assert(K <= 2, "fft_split_engine: too many accumulated products");
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

// Multiplies mod `mnum` by running NTTs modulo two FFT-friendly primes and CRT'ing.
// Inputs use balanced representatives (|v| <= MOD/2), so the true integer coefficients
// are bounded by n (MOD/2)^2 -- a 4x larger safe range than unsigned representatives.
template <typename mnum, typename num1 = mod_goldilocks, typename num2 = modnum<(15 << 27) + 1>>
struct crt_engine {
	using value_type = mnum;
	using E1 = fft_engine<num1>;
	using E2 = fft_engine<num2>;
	struct transformed {
		typename E1::transformed t1;
		typename E2::transformed t2;
		int size() const { return t1.size(); }
	};
	// K = number of accumulated term products (see the conv_engine preamble).
	template <int K> struct product_t {
		typename E1::product p1;
		typename E2::product p2;
		int size() const { return sz(p1); }
	};
	using product = product_t<1>;

	static int64_t balanced(mnum x) {
		int64_t v = int64_t(int(x));
		return 2 * v > int64_t(mnum::MOD) ? v - mnum::MOD : v;
	}

	// Folding (sz(a) <= 2n) is inherited from the inner engines' transforms.
	static transformed transform(std::span<const mnum> a, int n) {
		assert(sz(a) <= 2 * n);
		auto b1 = buffer_pool<num1>::get(sz(a));
		auto b2 = buffer_pool<num2>::get(sz(a));
		for (int i = 0; i < sz(a); i++) { int64_t v = balanced(a[i]); b1[i] = num1(v); b2[i] = num2(v); }
		return transformed{
			E1::transform(std::span<const num1>(b1.span()), n),
			E2::transform(std::span<const num2>(b2.span()), n),
		};
	}
	static void extend_to(transformed& t, int n, std::span<const mnum> coeffs) {
		auto b1 = buffer_pool<num1>::get(sz(coeffs));
		auto b2 = buffer_pool<num2>::get(sz(coeffs));
		for (int i = 0; i < sz(coeffs); i++) { int64_t v = balanced(coeffs[i]); b1[i] = num1(v); b2[i] = num2(v); }
		E1::extend_to(t.t1, n, std::span<const num1>(b1.span()));
		E2::extend_to(t.t2, n, std::span<const num2>(b2.span()));
	}
	static transformed even_half(const transformed& t, int n) {
		return transformed{E1::even_half(t.t1, n), E2::even_half(t.t2, n)};
	}
	static transformed odd_half(const transformed& t, int n) {
		return transformed{E1::odd_half(t.t1, n), E2::odd_half(t.t2, n)};
	}
	static transformed negate_arg(const transformed& t, int n) {
		return transformed{E1::negate_arg(t.t1, n), E2::negate_arg(t.t2, n)};
	}
	static product mul(const transformed& a, const transformed& b, int n) {
		return product{E1::mul(a.t1, b.t1, n), E2::mul(a.t2, b.t2, n)};
	}
	static product sq(const transformed& a, int n) { return mul(a, a, n); }
	template <int K1, int K2> static product_t<K1+K2> add(product_t<K1>&& a, product_t<K2>&& b) {
		return product_t<K1+K2>{E1::add(std::move(a.p1), b.p1), E2::add(std::move(a.p2), b.p2)};
	}
	template <int K = 1, typename Op = assign_op> static void finish(product_t<K>&& p, std::span<mnum> out, Op op = {}) {
		// The reconstruction needs |c| < whole/2; balanced inputs bound each addend's
		// true coefficients by n (MOD/2)^2, so the safe length is divided by the
		// addend count. K <= 2 is very conservative (~2^35 even for MOD ~ 2^30).
		static_assert(K <= 2, "crt_engine: too many accumulated products");
		int n = p.size();
		assert(sz(out) <= n);
		auto o1 = buffer_pool<num1>::get(sz(out));
		auto o2 = buffer_pool<num2>::get(sz(out));
		E1::finish(std::move(p.p1), o1.span());
		E2::finish(std::move(p.p2), o2.span());

		// TODO: Could hardcode these
		num1 inv_n2 = inv(num1(num2::MOD));
		num2 inv_n1 = inv(num2(num1::MOD));
		__int128_t whole = __int128_t(num1::MOD) * __int128_t(num2::MOD);

		mnum m1_mod = mnum(num1::MOD);
		mnum m2_mod = mnum(num2::MOD);
		mnum whole_mod = m1_mod * m2_mod;
		for (int i = 0; i < sz(out); i++) {
			num1 v1 = o1[i] * inv_n2;
			num2 v2 = o2[i] * inv_n1;
			mnum o_mod = mnum(uint64_t(v1)) * m2_mod + mnum(int(v2)) * m1_mod;
			__int128_t o_exact = __int128_t(uint64_t(v1)) * __int128_t(num2::MOD) + __int128_t(int(v2)) * __int128_t(num1::MOD);
			if (o_exact >= whole) { o_exact -= whole; o_mod -= whole_mod; }
			// Balanced representative: coefficients with negative true values (e.g.
			// from negate_arg'd transforms) reconstruct as whole + c, so values are
			// only required to have |c| < whole/2 rather than 0 <= c < whole.
			if (o_exact > whole / 2) o_mod -= whole_mod;
			op(out[i], o_mod);
		}
	}
};

static_assert(conv_engine<fft_engine<modnum<998244353>>>);
static_assert(conv_engine<fft_engine<mod_goldilocks>>);
static_assert(conv_engine<fft_real_engine<double>>);
static_assert(conv_engine<fft_split_engine<modnum<int(1e9)+7>>>);
static_assert(conv_engine<crt_engine<modnum<int(1e9)+7>>>);

// The engine-level cached operand: a transform of a coefficient sequence plus its
// logical length (which drives product sizes and the 2^k+1 cut). It does not own
// coefficients: it is built at a fixed size, and grows only when the owner of the
// coefficients feeds them back through extend_to (each doubling costs half a forward
// transform). Any power-of-two prefix of the transform is usable directly (see the
// fft_core preamble), so one transform serves every size up to its own.
template <conv_engine E> struct fft_cache {
	using T = typename E::value_type;
	typename E::transformed t;
	int n = 0;

	fft_cache() = default;
	// transform of coeffs at size `size` (a power of two; sz(coeffs) <= 2 * size)
	fft_cache(std::span<const T> coeffs, int size)
			: t(E::transform(coeffs, size)), n(sz(coeffs)) {}
	// adopt an existing transform (e.g. a pointwise product) of the sequence coeffs
	fft_cache(typename E::transformed t_, std::span<const T> coeffs)
			: t(std::move(t_)), n(sz(coeffs)) {}

	int len() const { return n; }
	int size() const { return t.size(); }
	const typename E::transformed& at_size(int m) const {
		assert(!(m & (m-1)) && m <= t.size());
		return t;
	}
	// Build (first call) or grow the transform to size m; the coefficient owner must
	// feed the same sequence every time. First build happens at
	// max(m, nextPow2(len - 1)): a 2^k+1-length operand builds at 2^k (E::transform
	// folds the top coefficient circularly -- what the 2^k+1 cut consumes, and a
	// valid seed for later extension).
	void extend_to(std::span<const T> coeffs, int m) {
		assert(!(m & (m-1)));
		if (t.size() == 0) {
			n = sz(coeffs);
			t = E::transform(coeffs, max(m, nextPow2(max(n - 1, 1))));
		} else if (t.size() < m) {
			assert(sz(coeffs) == n);
			E::extend_to(t, m, coeffs);
		}
	}
};

// Cyclic convolution of length n (a power of two); sz(a), sz(b) <= 2n (tails fold
// circularly). Writes coefficients [0, sz(out)) of the cyclic product through op;
// requires sz(out) <= n. Safe for out to alias a or b.
template <conv_engine E, typename Op = assign_op>
void multiply_circular(std::span<const typename E::value_type> a, std::span<const typename E::value_type> b,
		std::span<typename E::value_type> out, int n, Op op = {}) {
	assert(!(n & (n-1)));
	auto ta = E::transform(a, n);
	auto tb = E::transform(b, n);
	E::finish(E::mul(ta, tb, n), out, op);
}

// Same out contract as multiply_circular.
template <conv_engine E, typename Op = assign_op>
void square_circular(std::span<const typename E::value_type> a, std::span<typename E::value_type> out, int n, Op op = {}) {
	assert(!(n & (n-1)));
	auto ta = E::transform(a, n);
	E::finish(E::sq(ta, n), out, op);
}

namespace detail {
// Circular size for a product of linear length s: n = nextPow2(s), except when
// s == 2^k + 1 the product runs at 2^k (the cut; the one wrapped coefficient is
// corrected by emit_linear / emit_middle).
struct conv_size { int n; bool cut; };
inline conv_size conv_size_for(int s) {
	int n = nextPow2(s);
	bool cut = (n == 2 * (s - 1));
	return {cut ? n / 2 : n, cut};
}

// Shared tail for linear multiplication: applies the 2^k+1 wraparound correction (when
// s - 1 == n, the only wrapped coefficient is c_n, recoverable since c_0 = a_0 b_0),
// then writes through op.
template <typename T, typename Op>
void emit_linear(std::span<T> buf, int n, int s, bool cut, T c0, std::span<T> out, Op op) {
	T cn{};
	if (cut) {
		cn = buf[0] - c0;
		buf[0] = c0;
	}
	int lim = min(sz(out), min(s, n));
	for (int i = 0; i < lim; i++) op(out[i], buf[i]);
	if (cut && sz(out) >= s) op(out[s-1], cn);
}
}

// Linear multiplication: writes coefficients [0, min(sz(out), sz(a)+sz(b)-1)) of a*b
// through op; a shorter out truncates the product, a longer one leaves the tail
// untouched. Safe for out to alias a or b.
template <conv_engine E, typename Op = assign_op>
void multiply(std::span<const typename E::value_type> a, std::span<const typename E::value_type> b,
		std::span<typename E::value_type> out, Op op = {}) {
	using T = typename E::value_type;
	if (sz(a) == 0 || sz(b) == 0) return;
	int s = sz(a) + sz(b) - 1;
	auto [n, cut] = detail::conv_size_for(s);
	T c0 = a[0] * b[0];
	auto buf = buffer_pool<T>::get(n);
	multiply_circular<E>(a, b, buf.span(), n);
	detail::emit_linear<T>(buf.span(), n, s, cut, c0, out, op);
}

// Cached-operand form: each operand is its coefficients plus the fft_cache the
// owner keeps for them, (re)built or extended in place to the product size. Same
// out contract as the span multiply.
template <conv_engine E, typename Op = assign_op>
void multiply(std::span<const typename E::value_type> a, fft_cache<E>& ta,
		std::span<const typename E::value_type> b, fft_cache<E>& tb,
		std::span<typename E::value_type> out, Op op = {}) {
	using T = typename E::value_type;
	if (sz(a) == 0 || sz(b) == 0) return;
	int s = sz(a) + sz(b) - 1;
	auto [n, cut] = detail::conv_size_for(s);
	T c0 = a[0] * b[0];
	ta.extend_to(a, n);
	tb.extend_to(b, n);
	auto buf = buffer_pool<T>::get(n);
	E::finish(E::mul(ta.at_size(n), tb.at_size(n), n), buf.span());
	detail::emit_linear<T>(buf.span(), n, s, cut, c0, out, op);
}

// out op= a1*b1 + a2*b2, summing the two products with E::add so only one inverse
// transform is paid. Requires the two products to have the same linear length
// s = sz(a) + sz(b) - 1 (so they share a circular size); the 2^k+1 cut still
// applies, since the wraparound correction is linear in the summed product. Writes
// coefficients [0, min(sz(out), s)) of the sum through op.
template <conv_engine E, typename Op = assign_op>
void multiply_add2(std::span<const typename E::value_type> a1, fft_cache<E>& ta1,
		std::span<const typename E::value_type> b1, fft_cache<E>& tb1,
		std::span<const typename E::value_type> a2, fft_cache<E>& ta2,
		std::span<const typename E::value_type> b2, fft_cache<E>& tb2,
		std::span<typename E::value_type> out, Op op = {}) {
	using T = typename E::value_type;
	assert(sz(a1) > 0 && sz(b1) > 0 && sz(a2) > 0 && sz(b2) > 0);
	int s = sz(a1) + sz(b1) - 1;
	assert(sz(a2) + sz(b2) - 1 == s);
	auto [n, cut] = detail::conv_size_for(s);
	T c0 = a1[0] * b1[0] + a2[0] * b2[0];
	ta1.extend_to(a1, n); tb1.extend_to(b1, n);
	ta2.extend_to(a2, n); tb2.extend_to(b2, n);
	auto p = E::add(E::mul(ta1.at_size(n), tb1.at_size(n), n), E::mul(ta2.at_size(n), tb2.at_size(n), n));
	auto buf = buffer_pool<T>::get(n);
	E::finish(std::move(p), buf.span());
	detail::emit_linear<T>(buf.span(), n, s, cut, c0, out, op);
}

// Multiply two cached operands into coefficients plus (for engines where a product
// is itself a transform, product == transformed) the pointwise product transform:
// the product's forward FFT comes free (only the inverse is paid) and larger sizes
// extend from it -- worth it whenever the product's transform is consumed later
// (subproduct-tree build, repeated multiplication). This composes with the 2^k+1
// cut: the size-(s-1) product transform is the transform of the coefficients mod
// x^(s-1) - 1, a valid extension seed. Other engines get an empty (lazily built)
// transform. The result's coefficients live in `coeffs`; `t` is the fft_cache over
// them (owned by whoever owns coeffs).
template <conv_engine E>
void multiply_cached(std::span<const typename E::value_type> a, fft_cache<E>& ta,
		std::span<const typename E::value_type> b, fft_cache<E>& tb,
		std::vector<typename E::value_type>& coeffs, fft_cache<E>& t) {
	using T = typename E::value_type;
	coeffs.assign(size_t(sz(a) && sz(b) ? sz(a) + sz(b) - 1 : 0), T{});
	t = fft_cache<E>();
	if (coeffs.empty()) return;
	int s = sz(coeffs);
	if constexpr (std::same_as<typename E::product, typename E::transformed>) {
		auto [n, cut] = detail::conv_size_for(s);
		T c0 = a[0] * b[0];
		ta.extend_to(a, n);
		tb.extend_to(b, n);
		auto p = E::mul(ta.at_size(n), tb.at_size(n), n);
		auto tp = p;
		auto buf = buffer_pool<T>::get(n);
		E::finish(std::move(p), buf.span());
		detail::emit_linear<T>(buf.span(), n, s, cut, c0, std::span<T>(coeffs), assign_op{});
		t = fft_cache<E>(std::move(tp), std::span<const T>(coeffs));
	} else {
		multiply<E>(a, ta, b, tb, std::span<T>(coeffs));
	}
}

// Same out contract as multiply: writes coefficients [0, min(sz(out), 2 sz(a) - 1)).
template <conv_engine E, typename Op = assign_op>
void square(std::span<const typename E::value_type> a, std::span<typename E::value_type> out, Op op = {}) {
	using T = typename E::value_type;
	if (sz(a) == 0) return;
	int s = 2 * sz(a) - 1;
	auto [n, cut] = detail::conv_size_for(s);
	T c0 = a[0] * a[0];
	auto buf = buffer_pool<T>::get(n);
	square_circular<E>(a, buf.span(), n);
	detail::emit_linear<T>(buf.span(), n, s, cut, c0, out, op);
}

// Cached-operand form; same out contract as the span square.
template <conv_engine E, typename Op = assign_op>
void square(std::span<const typename E::value_type> a, fft_cache<E>& ta,
		std::span<typename E::value_type> out, Op op = {}) {
	using T = typename E::value_type;
	if (sz(a) == 0) return;
	int s = 2 * sz(a) - 1;
	auto [n, cut] = detail::conv_size_for(s);
	T c0 = a[0] * a[0];
	ta.extend_to(a, n);
	auto buf = buffer_pool<T>::get(n);
	E::finish(E::sq(ta.at_size(n), n), buf.span());
	detail::emit_linear<T>(buf.span(), n, s, cut, c0, out, op);
}

template <conv_engine E> vector<typename E::value_type> multiply(
		const vector<typename E::value_type>& a, const vector<typename E::value_type>& b) {
	using T = typename E::value_type;
	if (sz(a) == 0 || sz(b) == 0) return {};
	vector<T> r(sz(a) + sz(b) - 1);
	multiply<E>(std::span<const T>(a), std::span<const T>(b), std::span<T>(r));
	return r;
}

template <conv_engine E> vector<typename E::value_type> square(const vector<typename E::value_type>& a) {
	using T = typename E::value_type;
	if (sz(a) == 0) return {};
	vector<T> r(2 * sz(a) - 1);
	square<E>(std::span<const T>(a), std::span<T>(r));
	return r;
}

namespace detail {
// Shared tail for the middle product: buf holds the size-n circular product of a and
// b; writes coefficients [lb-1, lb-1+sz(out)) of a*b (clamped to la) through op. When
// la == n + 1 (the 2^k+1 cut) two slots carry wraparound: slot lb-1 also holds the
// product's top term ctop = a_top b_top, and c_n is recovered from slot 0 since
// c_0 = a_0 b_0 = c0.
template <typename T, typename Op>
void emit_middle(std::span<T> buf, bool cut, int la, int lb, T c0, T ctop, std::span<T> out, Op op) {
	int m = la - lb + 1;
	T cn{};
	if (cut) {
		cn = buf[0] - c0; // for lb == 1 these coincide: slot 0 = c_0 + c_n and ctop = c_n
		buf[lb - 1] -= ctop;
	}
	int lim = min(sz(out), cut ? m - 1 : m);
	for (int t = 0; t < lim; t++) op(out[t], buf[lb - 1 + t]);
	if (cut && sz(out) >= m) op(out[m-1], cn);
}
}

// Middle product (the transposed multiplication): writes coefficients
// [sz(b)-1, sz(b)-1 + min(sz(out), sz(a)-sz(b)+1)) of a*b through op -- the
// applications of the length-sz(b) sliding dot product of rev(b) against a -- into
// out[0..). Requires sz(a) >= sz(b) > 0. Uses a circular convolution of size ~sz(a),
// exploiting the case sz(a) == 2^k + 1 to stay at size 2^k.
template <conv_engine E, typename Op = assign_op>
void middle_product(std::span<const typename E::value_type> a, std::span<const typename E::value_type> b,
		std::span<typename E::value_type> out, Op op = {}) {
	using T = typename E::value_type;
	if (sz(a) == 0 || sz(b) == 0) return;
	assert(sz(a) >= sz(b));
	if (sz(a) == sz(b)) {
		T r{};
		for (int i = 0; i < sz(a); i++) {
			r += a[i] * b[sz(b) - 1 - i];
		}
		if (sz(out) > 0) op(out[0], r);
		return;
	}
	auto [n, cut] = detail::conv_size_for(sz(a));
	auto buf = buffer_pool<T>::get(n);
	multiply_circular<E>(a, b, buf.span(), n);
	detail::emit_middle<T>(buf.span(), cut, sz(a), sz(b),
			a[0] * b[0], a[sz(a) - 1] * b[sz(b) - 1], out, op);
}

template <conv_engine E> vector<typename E::value_type> middle_product(
		std::span<const typename E::value_type> a, std::span<const typename E::value_type> b) {
	using T = typename E::value_type;
	if (sz(a) == 0 || sz(b) == 0) return {};
	assert(sz(a) >= sz(b));
	vector<T> r(size_t(sz(a) - sz(b) + 1));
	middle_product<E>(a, b, std::span<T>(r));
	return r;
}

// Cached-transform middle product; see the span overload above for the contract.
// The coefficient owners pass coefficients alongside the transforms: the transforms
// are (re)built or extended in place to the product size (the 2^k+1 cut reads a
// prefix -- a transform prefix is the transform of the sequence mod x^n - 1, which
// is exactly the folded operand), and the equal-length dot product plus the cut
// correction read the coefficients directly.
template <conv_engine E, typename Op = assign_op>
void middle_product(std::span<const typename E::value_type> a, fft_cache<E>& ta,
		std::span<const typename E::value_type> b, fft_cache<E>& tb,
		std::span<typename E::value_type> out, Op op = {}) {
	using T = typename E::value_type;
	if (sz(a) == 0 || sz(b) == 0) return;
	assert(sz(a) >= sz(b));
	if (sz(a) == sz(b)) {
		T r{};
		for (int i = 0; i < sz(a); i++) {
			r += a[i] * b[sz(b) - 1 - i];
		}
		if (sz(out) > 0) op(out[0], r);
		return;
	}
	auto [n, cut] = detail::conv_size_for(sz(a));
	ta.extend_to(a, n);
	tb.extend_to(b, n);
	auto buf = buffer_pool<T>::get(n);
	E::finish(E::mul(ta.at_size(n), tb.at_size(n), n), buf.span());
	detail::emit_middle<T>(buf.span(), cut, sz(a), sz(b),
			a[0] * b[0], a[sz(a) - 1] * b[sz(b) - 1], out, op);
}

template <conv_engine E>
vector<typename E::value_type> middle_product(std::span<const typename E::value_type> a, fft_cache<E>& ta,
		std::span<const typename E::value_type> b, fft_cache<E>& tb) {
	using T = typename E::value_type;
	if (sz(a) == 0 || sz(b) == 0) return {};
	assert(sz(a) >= sz(b));
	vector<T> r(size_t(sz(a) - sz(b) + 1));
	middle_product<E>(a, ta, b, tb, std::span<T>(r));
	return r;
}

// Newton inversion: out[0..sz(a)) = coefficients of 1/a mod x^sz(a); requires
// sz(out) >= sz(a) (exactly sz(a) coefficients are written). Generic over any
// engine; per doubling step n -> m = 2n this is 5 transforms of size m, reusing b's
// transform for both circular products; in each product the wraparound only
// contaminates coefficients [0, n) which are already known. (The classic exact
// b*(2 - a*b) step via a pointwise triple product needs 3 transforms of size 2m --
// 6m vs 5m butterfly units -- and benchmarks ~20% slower.)
template <conv_engine E> void inverse(std::span<const typename E::value_type> a, std::span<typename E::value_type> out) {
	using T = typename E::value_type;
	if (sz(a) == 0) return;
	int s = nextPow2(sz(a));
	vector<T> b(s, T{});
	b[0] = inv(a[0]);
	for (int n = 1; n < sz(a); n *= 2) {
		int m = 2 * n;
		auto ta = E::transform(a.first(min(sz(a), m)), m);
		auto tb = E::transform(std::span<const T>(b).first(n), m);
		// e = a*b mod x^m; only e[n..m) is needed (and is wraparound-free).
		auto e = buffer_pool<T>::get(m);
		E::finish(E::mul(ta, tb, m), e.span());
		for (int i = 0; i < n; i++) e[i] = T{};
		auto te = E::transform(std::span<const T>(e.span()), m);
		auto c = buffer_pool<T>::get(m);
		// b' = 2b - b*(a*b): keep b on the left of e = a*b
		E::finish(E::mul(tb, te, m), c.span());
		for (int i = n; i < min(m, sz(a)); i++) b[i] = -c[i];
	}
	std::copy(b.begin(), b.begin() + sz(a), out.begin());
}

template <conv_engine E> vector<typename E::value_type> inverse(const vector<typename E::value_type>& a) {
	using T = typename E::value_type;
	vector<T> r(a.size());
	inverse<E>(std::span<const T>(a), std::span<T>(r));
	return r;
}

/* namespace fft */ }

// Power series; these are assumed to be the min of the length
template <typename E>
struct power_series : public std::vector<typename E::value_type> {
	using T = typename E::value_type;
	using std::vector<T>::vector;

	int ssize() const {
		return int(this->size());
	}
	int len() const {
		return ssize();
	}
	int degree() const {
		return len() - 1;
	}
	void extend(int sz) {
		assert(sz >= ssize());
		this->resize(sz);
	}
	void shrink(int sz) {
		assert(sz <= ssize());
		this->resize(sz);
	}
	// multiply by x^n
	void shift(int n = 1) {
		assert(n >= 0 && n <= ssize());
		std::rotate(this->begin(), this->end()-n, this->end());
		std::fill(this->begin(), this->begin()+n, T(0));
	}
	// divide by x^n and 0-pad
	void unshift(int n = 1) {
		assert(n >= 0 && n <= ssize());
		std::fill(this->begin(), this->begin()+n, T(0));
		std::rotate(this->begin(), this->begin()+n, this->end());
	}
	power_series& operator += (const power_series& o) {
		assert(len() == o.len());
		for (int i = 0; i < int(o.size()); i++) {
			(*this)[i] += o[i];
		}
		return *this;
	}
	friend power_series operator + (const power_series& a, const power_series& b) {
		power_series r(std::min(a.size(), b.size()));
		for (int i = 0; i < r.len(); i++) {
			r[i] = a[i] + b[i];
		}
		return r;
	}
	power_series& operator -= (const power_series& o) {
		assert(len() == o.len());
		for (int i = 0; i < int(o.size()); i++) {
			(*this)[i] -= o[i];
		}
		return *this;
	}
	friend power_series operator - (const power_series& a, const power_series& b) {
		power_series r(std::min(a.size(), b.size()));
		for (int i = 0; i < r.len(); i++) {
			r[i] = a[i] - b[i];
		}
		return r;
	}

	power_series& operator *= (const T& n) {
		for (auto& v : *this) v *= n;
		return *this;
	}
	friend power_series operator * (const power_series& a, const T& n) {
		power_series r(a.size());
		for (int i = 0; i < a.len(); i++) {
			r[i] = a[i] * n;
		}
		return r;
	}
	friend power_series operator * (const T& n, const power_series& a) {
		power_series r(a.size());
		for (int i = 0; i < a.len(); i++) {
			r[i] = n * a[i];
		}
		return r;
	}

	friend power_series operator * (const power_series& a, const power_series& b) {
		if (sz(a) == 0 || sz(b) == 0) return {};
		power_series r(std::min(a.size(), b.size()));
		fft::multiply<E>(std::span<const T>(a), std::span<const T>(b), std::span<T>(r));
		return r;
	}
	power_series& operator *= (const power_series& o) {
		return *this = (*this) * o;
	}
	friend power_series square(const power_series& a) {
		if (sz(a) == 0) return {};
		power_series r(a.size());
		fft::square<E>(std::span<const T>(a), std::span<T>(r));
		return r;
	}

	friend power_series inverse(const power_series& a) {
		power_series r(a.size());
		fft::inverse<E>(std::span<const T>(a), std::span<T>(r));
		return r;
	}
	// TODO: operator / can be done slightly faster than inverse:
	// we only need the n/2 terms of inverse(), and can do the last Newton step directly on the quotient

	friend power_series stretch(const power_series& a, int n) {
		power_series r(a.size());
		for (int i = 0; i*n < int(a.size()); i++) {
			r[i*n] = a[i];
		}
		return r;
	}
	friend power_series deriv_shift(power_series a) {
		for (int i = 0; i < a.len(); i++) {
			a[i] *= i;
		}
		return a;
	}
	friend power_series integ_shift(power_series a) {
		assert(a[0] == 0);
		T f = 1;
		for (int i = 1; i < int(a.size()); i++) {
			a[i] *= f;
			f *= i;
		}
		f = inv(f);
		for (int i = int(a.size()) - 1; i > 0; i--) {
			a[i] *= f;
			f *= i;
		}
		return a;
	}
	friend power_series integ_shift_offset(power_series a, int offset) {
		T f = 1;
		for (int i = 0; i < int(a.size()); i++) {
			a[i] *= f;
			f *= i + offset;
		}
		assert(f != 0);
		f = inv(f);
		for (int i = int(a.size()) - 1; i >= 0; i--) {
			a[i] *= f;
			f *= i + offset;
		}
		return a;
	}
	friend power_series deriv_shift_log(power_series a) {
		auto r = deriv_shift(a);
		return r * inverse(a);
	}
	friend power_series poly_log(power_series a) {
		assert(a[0] == 1);
		return integ_shift(deriv_shift_log(std::move(a)));
	}
	friend power_series poly_exp(power_series a) {
		// See https://mathexp.eu/bostan/publications/BoSc09a.pdf for details
		assert(a.size() >= 1);
		assert(a[0] == 0);
		power_series r(1, T(1)); r.reserve(a.size());
		power_series invR(1, T(1)); invR.reserve(a.size());
		while (r.size() < a.size()) {
			int o_sz = int(r.size());
			int n_sz = std::min(o_sz * 2, int(a.size()));
			power_series t = deriv_shift(power_series(a.begin(), a.begin() + o_sz));
			fft::multiply_circular<E>(std::span<const T>(t), std::span<const T>(r).first(o_sz), std::span<T>(t), o_sz);
			t = deriv_shift(r) - t;
			t *= invR;
			t.resize(n_sz - o_sz);
			power_series v(a.begin() + o_sz, a.begin() + n_sz);
			v -= integ_shift_offset(t, o_sz);
			v *= r;
			r.resize(n_sz);
			std::copy(v.begin(), v.end(), r.begin() + o_sz);
			if (r.size() < a.size()) {
				// double invR via a Newton step
				assert(int(r.size()) == 2 * int(invR.size()));
				int n = int(invR.size());
				int nn = int(r.size());
				power_series tmp(4 * n);
				fft::square<E>(std::span<const T>(invR).first(n), std::span<T>(tmp));
				fft::multiply<E>(std::span<const T>(tmp).first(nn), std::span<const T>(r).first(nn), std::span<T>(tmp));
				invR.resize(nn);
				for (int i = n; i < nn; i++) invR[i] = -tmp[i];
			}
		}
		return r;
	}
	friend power_series poly_pow_monic(power_series a, T k) {
		if (a.empty()) return a;
		assert(a.size() >= 1);
		assert(a[0] == 1);
		a = poly_log(a);
		a *= k;
		return poly_exp(a);
	}
	friend power_series poly_pow(power_series a, int64_t k) {
		assert(k >= 0);
		if (k == 0) {
			power_series r(a.len(), T(0));
			if (r.len() > 0) r[0] = T(1);
			return r;
		}

		int st = 0;
		while (st < a.len() && a[st] == 0) st++;

		if (st > 0 && k > (a.len() - 1) / st) {
			return power_series(a.len(), T(0));
		}

		power_series r(a.begin() + st, a.end() - (st * (k-1)));
		T leading_coeff = r[0];
		r *= inv(leading_coeff);
		r = poly_pow_monic(r, T(k));
		r *= power(leading_coeff, k);
		r.insert(r.begin(), st * k, T(0));
		assert(r.len() == a.len());
		return r;
	}

	friend power_series to_newton_sums(const power_series& a, int deg) {
		auto r = deriv_shift_log(a);
		r[0] = deg;
		for (int i = 1; i < int(r.size()); i++) r[i] = -r[i];
		return r;
	}
	friend power_series from_newton_sums(power_series S, int deg) {
		assert(S[0] == int(deg));
		S[0] = 0;
		for (int i = 1; i < int(S.size()); i++) S[i] = -S[i];
		return poly_exp(integ_shift(std::move(S)));
	}

	// Calculates prod 1/(1-x^i)^{a[i]}
	friend power_series euler_transform(const power_series& a) {
		power_series r = deriv_shift(a);
		std::vector<bool> is_prime(a.size(), true);
		for (int p = 2; p < int(a.size()); p++) {
			if (!is_prime[p]) continue;
			for (int i = 1; i*p < int(a.size()); i++) {
				r[i*p] += r[i];
				is_prime[i*p] = false;
			}
		}
		return poly_exp(integ_shift(r));
	}
	friend power_series inverse_euler_transform(const power_series& a) {
		power_series r = deriv_shift(poly_log(a));
		std::vector<bool> is_prime(a.size(), true);
		for (int p = 2; p < int(a.size()); p++) {
			if (!is_prime[p]) continue;
			for (int i = (int(a.size())-1)/p; i >= 1; i--) {
				r[i*p] -= r[i];
				is_prime[i*p] = false;
			}
		}
		return integ_shift(r);
	}

	// Calculates f(g(x)) mod x^n where deg(g) == n
	friend power_series poly_compose(const power_series& f, const power_series& g) {
		if (sz(g) == 0) return {};

		int m = int(f.size());
		int n = int(g.size());

		// https://arxiv.org/pdf/2404.05177
		// Consider P(y) = f(1/y) has terms from y^{-(m-1)}...y^0 (Laurent series)
		// We want [y^0] P(y) / (1 - y g(x))
		// Let Q_0 = 1 - yg(x)
		// Q_{i+1}(x^2, y) = Q_i(x, y) * Q_i(-x, y) mod x^{ceil(n / 2^i)}
		// deg_y(Q_i) = 2^i, deg_x(Q_i) = ceil(n / 2^i) - 1
		//
		// [y^0] P(y) / Q_l(x^2^l, y) * Q_{l-1}(-x^2^{l-1}, y) * Q_{l-2}(-x^2^{l-2}, y) * ... * Q_0(-x, y)
		// The total y deg of Q_{k-1} ... Q_0 is 2^k-1
		int L = __builtin_ctz(unsigned(nextPow2(n)));
		std::vector<power_series> Q(L+1);
		Q[0] = power_series(4 << L);
		Q[0][0] = 1;
		for (int i = 0; i < n; i++) Q[0][(2 << L) + i] = -g[i];
		for (int l = 1; l <= L; l++) {
			auto a = Q[l-1];
			// negate in place
			for (int i = 1; i < (4 << L); i += 2) Q[l-1][i] = -Q[l-1][i];
			Q[l] = power_series(4 << L);
			// TODO: Could be much more efficient:
			// We only need to do 1 forward FFT and then reflect it, and the backwards can be half the size and compactify at the same time.
			// We could also cache the forward FFT for the backwards pass
			fft::multiply_circular<E>(std::span<const T>(a), std::span<const T>(Q[l-1]), std::span<T>(Q[l]), 4 << L);
			// Compactify
			for (int i = 1; i < (2 << L); i++) {
				Q[l][i] = Q[l][2*i];
			}
			// Undo the circularity since we know it's monic
			for (int i = 0; i < (2 << (L - l)); i++) {
				Q[l][(2 << L) + i] = Q[l][i];
				Q[l][i] = 0;
			}
			Q[l][(2 << L)] -= T(1);
			Q[l][0] = T(1);
			// Zero out xs which are too big
			std::fill(Q[l].begin() + (2 << L) + (1 << (L-l)), Q[l].end(), T(0));
			for (int i = 0; i < (2 << L); i += 2 << (L-l)) {
				for (int j = 0; j < (1 << (L-l)); j++) {
					Q[l][i + (1 << (L-l)) + j] = 0;
				}
			}
		}
		power_series P;
		{
			P = f;
			std::reverse(P.begin(), P.end());
			power_series QL((1 << L) + 1);
			for (int i = 0; i <= (1 << L); i++) {
				QL[i] = Q[L][2 * i];
			}
			QL.resize(m, T(0));
			P *= inverse(QL);
			std::reverse(P.begin(), P.end());
			P.resize(1 << L, T(0));
			std::reverse(P.begin(), P.end());
			P.resize(4 << L, T(0));
			for (int i = (1 << L) - 1; i > 0; i--) {
				P[2*i] = P[i];
				P[i] = T(0);
			}
		}
		for (int l = L-1; l >= 0; l--) {
			// Spread it out, clear the high terms
			for (int i = (2 << L) - 1; i > 0; i--) {
				T v = P[i];
				P[2*i] = ((2*i) & (1 << (L-l))) ? T(0) : v;
				P[i] = T(0);
			}
			fft::multiply_circular<E>(std::span<const T>(Q[l]), std::span<const T>(P), std::span<T>(P), 4 << L);
			for (int i = 0; i < (2 << L); i++) {
				P[i] = P[(2 << L) + i];
				P[(2 << L) + i] = T(0);
			}
		}
		return power_series(P.begin(), P.begin() + n);
	}
};


template <typename E>
std::vector<typename E::value_type> poly_evaluate(
		const std::vector<typename E::value_type>& poly, const std::vector<typename E::value_type>& pts) {
	using T = typename E::value_type;
	if (pts.empty()) return {};
	using ps = power_series<E>;
	std::vector<std::vector<T>> series(pts.size() * 2);
	for (int i = 0; i < int(pts.size()); i++) {
		series[pts.size() + i] = {T(1), -pts[i]};
	}
	for (int i = int(pts.size()) - 1; i > 0; i--) {
		series[i] = fft::multiply<E>(series[2*i], series[2*i+1]);
	}
	{
		ps root(poly.rbegin(), poly.rend());
		assert(int(series[1].size()) == int(pts.size()) + 1);
		ps pts_series(series[1].begin(), series[1].end());
		pts_series.resize(root.size());
		ps top = inverse(pts_series) * root;

		series[1] = top;

		// front-pad it back to pts.size()
		std::reverse(series[1].begin(), series[1].end());
		series[1].resize(pts.size(), T(0));
		std::reverse(series[1].begin(), series[1].end());
	}
	for (int i = 1; i < int(pts.size()); i++) {
		series[2*i+0] = fft::middle_product<E>(series[i], series[2*i+0]);
		series[2*i+1] = fft::middle_product<E>(series[i], series[2*i+1]);
		std::swap(series[2*i+0], series[2*i+1]);
	}
	std::vector<T> ans(pts.size());
	for (int i = 0; i < int(pts.size()); i++) {
		ans[i] = series[pts.size() + i][0];
	}
	return ans;
}

template <typename E>
std::vector<typename E::value_type> poly_interpolate(
		const std::vector<typename E::value_type>& pts, const std::vector<typename E::value_type>& vals) {
	using T = typename E::value_type;
	if (pts.empty()) return {};
	using ps = power_series<E>;
	std::vector<std::vector<T>> series(pts.size() * 2);
	for (int i = 0; i < int(pts.size()); i++) {
		series[pts.size() + i] = {T(1), -pts[i]};
	}
	for (int i = int(pts.size()) - 1; i > 0; i--) {
		series[i] = fft::multiply<E>(series[2*i], series[2*i+1]);
	}
	std::vector<std::vector<T>> series_down(pts.size() * 2);
	{
		assert(int(series[1].size()) == int(pts.size()) + 1);
		ps root(series[1].begin(), series[1].end() - 1);
		ps deriv_root = root;
		for (int i = 0; i < int(pts.size()); i++) {
			deriv_root[i] *= T(int(pts.size()) - i);
		}
		series_down[1] = inverse(root) * deriv_root;
	}
	for (int i = 1; i < int(pts.size()); i++) {
		series_down[2*i+0] = fft::middle_product<E>(series_down[i], series[2*i+1]);
		series_down[2*i+1] = fft::middle_product<E>(series_down[i], series[2*i+0]);
	}
	for (int i = 0; i < int(pts.size()); i++) {
		auto& s = series_down[pts.size() + i];
		assert(int(s.size()) == 1);
		s[0] = vals[i] / s[0];
	}
	for (int i = int(pts.size()) - 1; i > 0; i--) {
		auto a = fft::multiply<E>(series_down[2*i+0], series[2*i+1]);
		auto b = fft::multiply<E>(series_down[2*i+1], series[2*i+0]);
		assert(int(a.size()) == int(series_down[i].size()));
		assert(int(b.size()) == int(series_down[i].size()));
		for (int z = 0; z < int(series_down[i].size()); z++) series_down[i][z] = a[z] + b[z];
	}
	std::vector<T> top = series_down[1];
	std::reverse(top.begin(), top.end());
	return top;
}

// Online (relaxed) multiplication: computes the first 2N terms of f*g given the terms
// one at a time. Each completed power-of-two block's transform is cached and reused for
// all later block products against it.
template <fft::conv_engine E> struct online_multiplier {
	using T = typename E::value_type;
	int N; int i;
	std::vector<T> f, g;
	std::vector<T> res;
	std::vector<fft::fft_cache<E>> f_blocks, g_blocks; // level k: block [2^k, 2^{k+1})

	online_multiplier(int N_) : N(N_), i(0), f(N, T{}), g(N, T{}), res(2*N+1, T(0)) {}

	T peek() {
		return res[i];
	}

	void push(T v_f, T v_g) {
		assert(i < N);
		f[i] = v_f;
		g[i] = v_g;
		if (i == 0) {
			res[0] += v_f * v_g;
		} else {
			res[i] += v_f * g[0];
			res[i] += f[0] * v_g;
			for (int p = 1, k = 0; (i & (p-1)) == (p-1); p <<= 1, k++) {
				int lo1 = p;
				int lo2 = i + 1 - p;
				int s = 2*p - 1;
				auto fb = std::span<const T>(f).subspan(p, p);
				auto gb = std::span<const T>(g).subspan(p, p);
				auto out = std::span<T>(res).subspan(lo1 + lo2, s);
				if (i == 2*p-1) {
					f_blocks.emplace_back();
					g_blocks.emplace_back();
					fft::multiply<E>(fb, f_blocks[k], gb, g_blocks[k], out, fft::add_op{});
					break;
				}
				// both products keep f on the left: f_hi * g_lo + f_lo * g_hi
				fft::fft_cache<E> cf, cg;
				fft::multiply_add2<E>(
						fb, f_blocks[k], std::span<const T>(g).subspan(lo2, p), cg,
						std::span<const T>(f).subspan(lo2, p), cf, gb, g_blocks[k],
						out, fft::add_op{});
			}
		}
		i++;
	}

	T back() {
		return res[i-1];
	}
};

template <fft::conv_engine E> struct online_squarer {
	using T = typename E::value_type;
	int N; int i;
	std::vector<T> f;
	std::vector<T> res;
	std::vector<fft::fft_cache<E>> f_blocks;

	online_squarer(int N_) : N(N_), i(0), f(N, T{}), res(2*N+1, T(0)) {}

	T peek() {
		return res[i];
	}

	void push(T v_f) {
		assert(i < N);
		f[i] = v_f;
		if (i == 0) {
			res[0] += v_f * v_f;
		} else {
			res[i] += (v_f + v_f) * f[0];
			for (int p = 1, k = 0; (i & (p-1)) == (p-1); p <<= 1, k++) {
				int lo1 = p;
				int lo2 = i + 1 - p;
				int s = 2*p - 1;
				auto fb = std::span<const T>(f).subspan(p, p);
				auto out = std::span<T>(res).subspan(lo1 + lo2, s);
				if (i == 2*p-1) {
					f_blocks.emplace_back();
					fft::square<E>(fb, f_blocks[k], out, fft::add_op{});
					break;
				}
				fft::fft_cache<E> cw;
				fft::multiply<E>(fb, f_blocks[k], std::span<const T>(f).subspan(lo2, p), cw,
						out, fft::add_twice_op{});
			}
		}
		i++;
	}

	T back() {
		return res[i-1];
	}
};

// A polynomial represented by its values evaluated at an Arithmetic Progression (AP).
// TODO: The AP is always assumed to be 0..length-1; store an explicit offset/gap instead?
// Maybe not, this is just more convenient.
template <typename E>
struct poly_ap_values : public std::vector<typename E::value_type> {
	using T = typename E::value_type;
	using std::vector<T>::vector;

	int ssize() const {
		return int(this->size());
	}
	int len() const {
		return ssize();
	}
	int degree() const {
		return len() - 1;
	}

	poly_ap_values& operator += (const poly_ap_values& o) {
		assert(len() == o.len());
		for (int i = 0; i < int(o.size()); i++) {
			(*this)[i] += o[i];
		}
		return *this;
	}
	friend poly_ap_values operator + (const poly_ap_values& a, const poly_ap_values& b) {
		assert(a.size() == b.size());
		poly_ap_values r(a.size());
		for (int i = 0; i < r.len(); i++) {
			r[i] = a[i] + b[i];
		}
		return r;
	}
	poly_ap_values& operator -= (const poly_ap_values& o) {
		assert(len() == o.len());
		for (int i = 0; i < int(o.size()); i++) {
			(*this)[i] -= o[i];
		}
		return *this;
	}
	friend poly_ap_values operator - (const poly_ap_values& a, const poly_ap_values& b) {
		assert(a.size() == b.size());
		poly_ap_values r(a.size());
		for (int i = 0; i < r.len(); i++) {
			r[i] = a[i] - b[i];
		}
		return r;
	}

	T eval_at(T k) {
		if (0 <= int(k) && int(k) < len()) {
			return (*this)[int(k)];
		} else {
			// Just do the lagrange interpolation
			std::vector<T> terms(*this);
			{
				// Inverse factorial terms
				T v = 1;
				for (int i = 1; i <= len(); i++) v *= T(i);
				v = inv(v);
				for (int i = len()-1; i >= 0; i--) {
					v *= T(i+1);
					terms[i] *= v;
					terms[len()-1-i] *= (i & 1) ? -v : v;
				}
			}
			{
				// Prefix terms
				T v = 1;
				for (int i = 0; i < len(); i++) {
					terms[i] *= v;
					v *= T(k - i);
				}
			}
			{
				// Suffix terms
				T v = 1;
				for (int i = len() - 1; i >= 0; i--) {
					terms[i] *= v;
					v *= T(k - i);
				}
			}
			T res = 0;
			for (int i = 0; i < len(); i++) res += terms[i];
			return res;
		}
	}

	poly_ap_values eval_range(T k, int osz) {
		if (osz == 0) {
			return poly_ap_values(osz);
		}
		if (len() == 0) {
			return poly_ap_values(osz, T(0));
		}

		// just assume there's no overlap; TODO: Fix this

		std::vector<T> inps(*this);
		{
			// Inverse factorial terms
			T v = 1;
			for (int i = 1; i <= len(); i++) v *= T(i);
			v = inv(v);
			for (int i = len()-1; i >= 0; i--) {
				v *= T(i+1);
				inps[i] *= v;
				inps[len()-1-i] *= (i & 1) ? -v : v;
			}
		}
		std::vector<T> inv_offsets(len() + osz - 1);
		poly_ap_values results(osz);
		{
			T v = 1;
			for (int i = - (len() - 1); i <= osz - 1; i++) {
				inv_offsets[i + (len() - 1)] = v;
				v *= k + i;
				if (i >= 0) results[i] = v;
			}
			// Assert there's no overlap
			assert(v != T(0));
			v = inv(v);
			for (int i = osz - 1; i >= -(len() - 1); i--) {
				inv_offsets[i + (len() - 1)] *= v;
				v *= k + i;
				if (i + (len() - 1) <= osz - 1) {
					results[i + (len() - 1)] *= v;
				}
			}
		}
		std::vector<T> prod = fft::middle_product<E>(inv_offsets, inps);
		assert(int(prod.size()) == osz);
		for (int i = 0; i < osz; i++) results[i] *= prod[i];
		return results;
	}

	void extend_right() {
		this->push_back(eval_at(T(len())));
	}
	void extend_left() {
		this->insert(this->begin(), eval_at(T(-1)));
	}

	[[nodiscard]] poly_ap_values prefix_sum_inclusive() const {
		poly_ap_values r = *this;
		r.extend_right();
		for (int i = 1; i < r.len(); i++) r[i] += r[i-1];
		return r;
	}
};


/* namespace ecnerwala */ }
