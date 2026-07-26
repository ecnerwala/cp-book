#pragma once

#include <vector>
#include <cmath>
#include <utility>
#include <algorithm>
#include <span>
#include <cstdint>
#include <cassert>
#include <concepts>
#include <functional>

#include "modnum.hpp"

/**
 * Author: Andrew He
 * Source: http://neerc.ifmo.ru/trains/toulouse/2017/fft2.pdf
 * Papers about accuracy: http://www.daemonology.net/papers/fft.pdf, http://www.cs.berkeley.edu/~fateman/papers/fftvsothers.pdf
 * For integers rounding works if $(|a| + |b|)\max(a, b) < \mathtt{\sim} 10^9$, or in theory maybe $10^6$.
 *
 * Abstraction layers:
 *   fft_core<num>   FFT itself and other ops on rings with 2^k-th roots of unity. We use bit-reversed indexing in the frequency domain.
 *
 *   engines         Engines for packing/unpacking arbitrary rings for convolution (the `engine` concept).
 *                   Still expose (opaque) transform-domain objects for caching/fusion.
 *
 *   multiply layer  Wrappers for convolving bounded sequences: track length/truncation.
 *
 *   value types     series::vec<E> - R[[x]]
 *                   series::exact<E> - exact (finite-support) power series
 *                   series::trunc<E> - truncated prefix of an (infinite) power series
 *
 *                   polynomials - R[x]. Under x -> 1/x a polynomial becomes a Laurent polynomial in 1/x;
 *                                 shifting by x^{deg P} (reversal) lands it in R[[x]], and we store that exact series.
 *                   poly::vec<E> - polynomial type, supporting natural indexing
 *                   linear_form<E> - finite-support linear forms, via the pairing <P, S> = [x^0] P(1/x) S(x)
 *                                    a linear form is one side of this pairing, applied to the other
 *
 *                   online_multiplier<E> - online (relaxed) multiplication of 2 sequences in n log^2 n time
 *                   poly_ap_values<E> - a polynomial stored as its evaluations on an arithmetic progression
 */

namespace ecnerwala {

template<class T> int sz(T&& arg) { using std::size; return int(size(std::forward<T>(arg))); }
inline int nextPow2(int s) { return 1 << (s > 1 ? 32 - __builtin_clz(s-1) : 0); }

namespace fft {

using std::swap;
using std::vector;
using std::min;
using std::max;

// ==== core: roots, buffers, raw transforms ====

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

// We take the bit-reverse convention: the coefficient of a[i] -> b[j] is omega^{i * bit_reverse(j)}.
// This means that the size 2^{k-1} transform is the prefix of the size 2^k transform (wrapping the input).
//
// We mostly work with spans here:
//   spans of transforms are expected to have length exactly 2^k
//   spans of inputs/outputs are expected to have length [0, 2^{k+1})
//
// Inputs/outputs are treated mod x^{2^k} - 1.
// Their length is allowed to be bigger than 2^k mostly to perform ops on sequences of size n+1 with only transforms of size n.
// The upper bound of 2^{k+1} is arbitrary: we could tighten it to 2^k + 1 or loosen it to infinity, this is just a "defensive" choice.
template <typename num> struct fft_core {

	static inline vector<int> rev;
	// rt[2^k + i] = 1^{i / 2^(k+1)}
	// TODO: can we get rid of inv_rt; alternatively, should we store inv_rt in bit-reverse order?
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

	// bit-reversal of i as a log2(n)-bit number
	static int brev(int i, int n) {
		int s = __builtin_ctz(unsigned(sz(rev)/n));
		return rev[i] >> s;
	}
	// index of the conjugate evaluation point, in the returned bit-reversed order
	static int conj_index(int j) {
		return j == 0 ? 0 : j ^ ((1 << (31 - __builtin_clz(unsigned(j)))) - 1);
	}

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

	// Extend a size 2^{k-1} transform to size 2^k; we need the coefficients.
	// t must have size 2^k, and coeffs must have size at most 2^{k+1}.
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

	// Consider t = transform(P) and P(x) = E(x^2) + O(x^2) * x
	// `even_half` and `odd_half` extract a size 2^{k-1} transform of E/O, respectively.
	static void even_half(std::span<const num> t, std::span<num> out) {
		int n = sz(out);
		assert(sz(t) >= 2*n);
		num half = inv(num(2));
		for (int j = 0; j < n; j++) out[j] = (t[2*j] + t[2*j+1]) * half;
	}
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

// ==== engine concept ====

// Output operations for the finish step to express arbitrary fusion into the output buffer.
struct assign_op { template <typename T> void operator()(T& d, T v) const { d = v; } };
struct add_op { template <typename T> void operator()(T& d, T v) const { d += v; } };
struct sub_op { template <typename T> void operator()(T& d, T v) const { d -= v; } };
struct add_twice_op { template <typename T> void operator()(T& d, T v) const { d += v + v; } };

// `engine` contract
//   engine represents a way of packing/unpacking sequences over an arbitrary ring into FFT-style transforms.
//   We expect transforms/products of transforms to be linear but potentially lossy/imprecise, so we'll track precision
//   at compile-time as a template parameter.
//
//   E::value_type   The ring we operate over
//   E::unit_scale   0 or 1 depending on whether there's error that can accumulate
//   E::commutative  A marker for whether the ring is commutative
//
//   transformed_t<A>  The transform of a sequence. This object owns its data buffer.
//   product_t<A>      The product of 2 transforms. May equal transformed_t, particularly when unit_scale = 0.
//
//   transformed       alias for transformed_t<unit_scale>
//   product           alias for product_t<unit_scale>
//
//   The basic multiplication API is
//      transform(span<const value_type> in, int n) -> transformed_t<unit_scale>
//      mul(transformed_t<A>, transformed_t<B>, int n) -> product_t<A*B>
//      mul2(a1, b1, a2, b2, int n) -> product_t<A1*B1 + A2*B2>, computing a1*b1 + a2*b2 in one pass
//      finish(product_t<A>, span<value_type>& out, Op) -> void
//
//   Input span can be length up to 2n.
//   Output spans can be length up to n; only the prefix that exists is filled.
//   finish applies Op exactly once per out element, in index order, to
//   value_type targets (so ops may be stateful).
//   Transforms can be longer than necessary, and only the relevant prefix is used.
//
//   For non-exact engines, there's some subtlety in whether we wrap before or after packing.
//   We will choose to wrap *after* packing, which hurts error bounds but makes the prefix condition more uniform.
//
//   Additionally, we have APIs to take advantage of linearity in both transformed and product space:
//      add(transformed_t<A>, transformed_t<B>) -> transformed_t<A+B>
//      add(product_t<K1>, product_t<K2>) -> product_t<K1+K2>
//
//   Finally, we expose some additional fast-transform optimization paths.
//   extend_to only operates on transformed_t<unit_scale>; the others are scale-generic.
//   downsample is also defined on product_t (halving before finish saves inverse-transform work).
//      extend_to       build (if empty) or grow a transform to size m by repeated doubling; feed the same coefficients (sz <= 2m) every time
//      downsample      compute the half-sized transform/product of just the even (odd = false) or odd terms of the input
//      negate_arg      size n transform of A(-x)
template <typename E>
concept engine = requires(
	std::span<const typename E::value_type> in,
	std::span<typename E::value_type> out,
	typename E::transformed& t,
	const typename E::transformed& ct,
	typename E::product& p,
	const typename E::product& cp,
	int n
) {
	typename E::value_type;
	{ E::transform(in, n) } -> std::same_as<typename E::transformed>;
	{ ct.size() } -> std::same_as<int>;
	E::extend_to(t, n, in);
	{ E::downsample(ct, n, false) } -> std::same_as<typename E::transformed>;
	{ E::downsample(cp, n, false) } -> std::same_as<typename E::product>;
	{ E::negate_arg(ct, n) } -> std::same_as<typename E::transformed>;
	{ E::mul(ct, ct, n) } -> std::same_as<typename E::product>;
	{ E::sq(ct, n) } -> std::same_as<typename E::product>;
	{ E::mul2(ct, ct, ct, ct, n) } -> std::same_as<typename E::template product_t<2 * E::unit_scale>>;
	E::finish(std::move(p), out);
	E::finish(std::move(p), out, add_op{});
	E::finish(E::add(std::move(p), std::move(p)), out);
	{ E::add(E::transform(in, n), ct) } -> std::same_as<typename E::template transformed_t<2 * E::unit_scale>>;
	{ E::add(std::move(p), std::move(p)) } -> std::same_as<typename E::template product_t<2 * E::unit_scale>>;
	requires std::same_as<std::remove_cvref_t<decltype(E::commutative)>, bool>;
	requires std::same_as<std::remove_cvref_t<decltype(E::unit_scale)>, int>;
};

// ==== scalar engines ====

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
		assert(sz(coeffs) <= 2 * t.size());
		auto buf = buffer_pool<cnum>::get((sz(coeffs) + 1) / 2);
		std::fill(buf.span().begin(), buf.span().end(), cnum(0));
		pack(coeffs, buf.span());
		while (t.size() < m) {
			t.v.resize(2 * sz(t.v));
			core::extend(std::span<cnum>(t.v), std::span<const cnum>(buf.span()));
		}
	}
	static transformed downsample(const transformed& t, int n, bool odd) { return half(t, n, odd); }
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
};

// Multiplies mod `mnum` by splitting values into balanced 15-bit halves (each limb in
// [-2^14, 2^14], from the balanced representative |v| <= MOD/2) packed into one complex
// transform per operand.
// TODO: Add type bounds?
template <typename mnum> struct split {
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
		int64_t v = int64_t(int(x));
		if (2 * v > int64_t(mnum::MOD)) v -= mnum::MOD;
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
		assert(sz(coeffs) <= 2 * t.size());
		auto buf = buffer_pool<cnum>::get(sz(coeffs));
		for (int i = 0; i < sz(coeffs); i++) buf[i] = pack(coeffs[i]);
		while (t.size() < m) {
			t.v.resize(2 * t.size());
			core::extend(std::span<cnum>(t.v), std::span<const cnum>(buf.span()));
		}
	}
	static void downsample_core(std::span<const cnum> in, std::span<cnum> out, bool odd) {
		if (odd) core::odd_half(in, out);
		else core::even_half(in, out);
	}
	template <int A> static transformed_t<A> downsample(const transformed_t<A>& t, int n, bool odd) {
		transformed_t<A> r; r.v.resize(n);
		downsample_core(std::span<const cnum>(t.v), std::span<cnum>(r.v), odd);
		return r;
	}
	template <int K> static product_t<K> downsample(const product_t<K>& p, int n, bool odd) {
		product_t<K> r; r.lo.resize(n); r.hi.resize(n);
		downsample_core(std::span<const cnum>(p.lo), std::span<cnum>(r.lo), odd);
		downsample_core(std::span<const cnum>(p.hi), std::span<cnum>(r.hi), odd);
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

// Multiplies mod `mnum` by running NTTs modulo two FFT-friendly primes and CRT'ing.
// Inputs use balanced representatives (|v| <= MOD/2), so the true integer coefficients
// are bounded by n (MOD/2)^2.
template <typename mnum, typename num1 = mod_goldilocks, typename num2 = modnum<(15 << 27) + 1>>
struct crt {
	using value_type = mnum;
	static constexpr bool commutative = true;
	static constexpr int unit_scale = 1;
	using E1 = ntt<num1>;
	using E2 = ntt<num2>;
	template <int A = 1> struct transformed_t {
		typename E1::transformed t1;
		typename E2::transformed t2;
		int size() const { return t1.size(); }
		transformed_t() = default;
		transformed_t(typename E1::transformed&& t1_, typename E2::transformed&& t2_)
			: t1(std::move(t1_)), t2(std::move(t2_)) {}
		template <int A2> requires (A2 != A) explicit(A2 > A) transformed_t(transformed_t<A2>&& o)
			: t1(std::move(o.t1)), t2(std::move(o.t2)) {}
	};
	using transformed = transformed_t<1>;
	template <int K> struct product_t {
		typename E1::product p1;
		typename E2::product p2;
		int size() const { return sz(p1); }
		product_t() = default;
		product_t(typename E1::product&& p1_, typename E2::product&& p2_)
			: p1(std::move(p1_)), p2(std::move(p2_)) {}
		template <int K2> requires (K2 != K) explicit(K2 > K) product_t(product_t<K2>&& o)
			: p1(std::move(o.p1)), p2(std::move(o.p2)) {}
	};
	using product = product_t<1>;

	static int64_t balanced(mnum x) {
		int64_t v = int64_t(int(x));
		return 2 * v > int64_t(mnum::MOD) ? v - mnum::MOD : v;
	}

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
	static void extend_to(transformed& t, int m, std::span<const mnum> coeffs) {
		if (t.size() >= m) return;
		auto b1 = buffer_pool<num1>::get(sz(coeffs));
		auto b2 = buffer_pool<num2>::get(sz(coeffs));
		for (int i = 0; i < sz(coeffs); i++) { int64_t v = balanced(coeffs[i]); b1[i] = num1(v); b2[i] = num2(v); }
		E1::extend_to(t.t1, m, std::span<const num1>(b1.span()));
		E2::extend_to(t.t2, m, std::span<const num2>(b2.span()));
	}
	template <int A> static transformed_t<A> downsample(const transformed_t<A>& t, int n, bool odd) {
		return transformed_t<A>{E1::downsample(t.t1, n, odd), E2::downsample(t.t2, n, odd)};
	}
	template <int K> static product_t<K> downsample(const product_t<K>& p, int n, bool odd) {
		return product_t<K>{E1::downsample(p.p1, n, odd), E2::downsample(p.p2, n, odd)};
	}
	template <int A> static transformed_t<A> negate_arg(const transformed_t<A>& t, int n) {
		return transformed_t<A>{E1::negate_arg(t.t1, n), E2::negate_arg(t.t2, n)};
	}
	// Exact per prime; the scale tracks the true (integer) coefficient growth.
	template <int A, int B> static transformed_t<A + B> add(transformed_t<A>&& a, const transformed_t<B>& b) {
		return transformed_t<A + B>{E1::add(std::move(a.t1), b.t1), E2::add(std::move(a.t2), b.t2)};
	}
	template <int A, int B> static product_t<A * B> mul(const transformed_t<A>& a, const transformed_t<B>& b, int n) {
		return product_t<A * B>{E1::mul(a.t1, b.t1, n), E2::mul(a.t2, b.t2, n)};
	}
	template <int A> static product_t<A * A> sq(const transformed_t<A>& a, int n) { return mul(a, a, n); }
	template <int A1, int B1, int A2, int B2>
	static product_t<A1 * B1 + A2 * B2> mul2(
		const transformed_t<A1>& a1, const transformed_t<B1>& b1,
		const transformed_t<A2>& a2, const transformed_t<B2>& b2,
		int n
	) {
		return product_t<A1 * B1 + A2 * B2>{
			E1::mul2(a1.t1, b1.t1, a2.t1, b2.t1, n),
			E2::mul2(a1.t2, b1.t2, a2.t2, b2.t2, n),
		};
	}
	template <int K1, int K2> static product_t<K1 + K2> add(product_t<K1>&& a, product_t<K2>&& b) {
		return product_t<K1 + K2>{E1::add(std::move(a.p1), b.p1), E2::add(std::move(a.p2), b.p2)};
	}
	template <int K = 1, typename Op = assign_op> static void finish(product_t<K>&& p, std::span<mnum> out, Op op = {}) {
		// The reconstruction needs |c| < whole/2; balanced inputs bound each addend's
		// true coefficients by n (MOD/2)^2, so the safe length is divided by the
		// accumulated scale. K <= 2 is very conservative (~2^35 even for MOD ~ 2^30).
		static_assert(K <= 2, "crt: accumulated scale too large");
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
			// Balanced representatives: |o| <= whole/2
			if (o_exact > whole / 2) o_mod -= whole_mod;
			op(out[i], o_mod);
		}
	}
};

// Small NxN matrix over num, row-major
// ==== wrapper engines ====

template <typename num, int N> struct mat {
	std::array<num, size_t(N) * N> a{};
	num& operator[](std::array<int, 2> rc) { return a[size_t(rc[0]) * N + rc[1]]; }
	const num& operator[](std::array<int, 2> rc) const { return a[size_t(rc[0]) * N + rc[1]]; }
	num* data() { return a.data(); }
	const num* data() const { return a.data(); }
	mat& operator+=(const mat& o) { for (int i = 0; i < N*N; i++) a[i] += o.a[i]; return *this; }
	friend mat operator+(mat x, const mat& y) { x += y; return x; }
	mat& operator-=(const mat& o) { for (int i = 0; i < N*N; i++) a[i] -= o.a[i]; return *this; }
	friend mat operator-(mat x, const mat& y) { x -= y; return x; }
	friend mat operator*(const mat& x, const mat& y) {
		mat r;
		for (int i = 0; i < N; i++) for (int k = 0; k < N; k++) for (int j = 0; j < N; j++)
			r[{i, j}] += x[{i, k}] * y[{k, j}];
		return r;
	}
	mat& operator*=(const mat& o) { return *this = *this * o; }
	friend bool operator==(const mat&, const mat&) = default;
};

// Truncated polynomial mod x^N over num
template <typename num, int N> struct trunc_series {
	std::array<num, size_t(N)> a{};
	num& operator[](int i) { return a[size_t(i)]; }
	const num& operator[](int i) const { return a[size_t(i)]; }
	num* data() { return a.data(); }
	const num* data() const { return a.data(); }
	trunc_series& operator+=(const trunc_series& o) { for (int i = 0; i < N; i++) a[i] += o.a[i]; return *this; }
	friend trunc_series operator+(trunc_series x, const trunc_series& y) { x += y; return x; }
	trunc_series& operator-=(const trunc_series& o) { for (int i = 0; i < N; i++) a[i] -= o.a[i]; return *this; }
	friend trunc_series operator-(trunc_series x, const trunc_series& y) { x -= y; return x; }
	friend trunc_series operator*(const trunc_series& x, const trunc_series& y) {
		trunc_series r;
		for (int i = 0; i < N; i++) for (int j = 0; j < N - i; j++) r[i + j] += x[i] * y[j];
		return r;
	}
	trunc_series& operator*=(const trunc_series& o) { return *this = *this * o; }
	friend bool operator==(const trunc_series&, const trunc_series&) = default;
};

// componentwise

// These are "componentwise engines" which model free modules/algebras over the underlying ring.
// Matrices are the canonical example: we can take entry-wise transforms, then multiply/add in transformed space.
//
// We start with a shared `componentwise` base class which handles all linear ops, i.e. not mul.
//
// If the underlying has unit_scale = 1, we may need to avoid accumulation of sums in transformed space;
// then, the product and the transformed data may have different dimensions.
// We'll represent this by an array Ofs of prefix offsets mapping each input/transform-space dimension to a range of product-space dimensions.
// Specifically, out[c] = sum prod[Ofs[c]:Ofs[c+1]]

template <int L> constexpr std::array<int, size_t(L) + 1> componentwise_iota = [] {
	std::array<int, size_t(L) + 1> r{};
	for (int i = 0; i <= L; i++) r[size_t(i)] = i;
	return r;
}();

template <engine E, typename V, int L, std::array<int, size_t(L) + 1> Ofs = componentwise_iota<L>>
struct componentwise {
	using S = typename E::value_type;
	using value_type = V;
	static constexpr int P = Ofs[size_t(L)];  // total product components
	static constexpr int unit_scale = E::unit_scale;
	template <int A = unit_scale> struct transformed_t {
		std::array<typename E::template transformed_t<A>, size_t(L)> t;
		int size() const { return t[0].size(); }
		transformed_t() = default;
		template <int A2> requires (A2 != A) explicit(A2 > A) transformed_t(transformed_t<A2>&& o) {
			for (int c = 0; c < L; c++)
				t[c] = typename E::template transformed_t<A>(std::move(o.t[c]));
		}
	};
	using transformed = transformed_t<>;
	// TODO: if E::product_t == E::transformed_t, mirror that here
	template <int K> struct product_t {
		std::array<typename E::template product_t<K>, size_t(P)> t;
		int size() const { return t[0].size(); }
		product_t() = default;
		template <int K2> requires (K2 != K) explicit(K2 > K) product_t(product_t<K2>&& o) {
			for (int c = 0; c < P; c++)
				t[c] = typename E::template product_t<K>(std::move(o.t[c]));
		}
	};

	static transformed transform(std::span<const V> a, int n) {
		transformed r;
		auto buf = buffer_pool<S>::get(sz(a));
		for (int c = 0; c < L; c++) {
			for (int i = 0; i < sz(a); i++) buf[i] = a[i].data()[c];
			r.t[c] = E::transform(std::span<const S>(buf.span()), n);
		}
		return r;
	}
	static void extend_to(transformed& t, int m, std::span<const V> coeffs) {
		if (t.size() >= m) return;
		auto buf = buffer_pool<S>::get(sz(coeffs));
		for (int c = 0; c < L; c++) {
			for (int i = 0; i < sz(coeffs); i++) buf[i] = coeffs[i].data()[c];
			E::extend_to(t.t[c], m, std::span<const S>(buf.span()));
		}
	}
	template <int A> static transformed_t<A> downsample(const transformed_t<A>& t, int n, bool odd) {
		transformed_t<A> r;
		for (int c = 0; c < L; c++) r.t[c] = E::downsample(t.t[c], n, odd);
		return r;
	}
	template <int K> static product_t<K> downsample(const product_t<K>& p, int n, bool odd) {
		product_t<K> r;
		for (int c = 0; c < P; c++) r.t[c] = E::downsample(p.t[c], n, odd);
		return r;
	}
	template <int A> static transformed_t<A> negate_arg(const transformed_t<A>& t, int n) {
		transformed_t<A> r;
		for (int c = 0; c < L; c++) r.t[c] = E::negate_arg(t.t[c], n);
		return r;
	}
	template <int A, int B> static transformed_t<A + B> add(transformed_t<A>&& a, const transformed_t<B>& b) {
		transformed_t<A + B> r;
		for (int c = 0; c < L; c++) r.t[c] = E::add(std::move(a.t[c]), b.t[c]);
		return r;
	}
	template <int K1, int K2> static product_t<K1 + K2> add(product_t<K1>&& a, product_t<K2>&& b) {
		product_t<K1 + K2> r;
		for (int c = 0; c < P; c++) r.t[c] = E::add(std::move(a.t[c]), std::move(b.t[c]));
		return r;
	}
	template <int K, typename Op = assign_op> static void finish(product_t<K>&& p, std::span<V> out, Op op = {}) {
		auto buf = buffer_pool<S>::get(sz(out));
		auto emit = [&](std::span<V> dst) {
			for (int c = 0; c < L; c++) {
				E::finish(std::move(p.t[Ofs[size_t(c)]]), buf.span());
				for (int j = Ofs[size_t(c)] + 1; j < Ofs[size_t(c) + 1]; j++)
					E::finish(std::move(p.t[j]), buf.span(), add_op{});
				for (int i = 0; i < sz(dst); i++) dst[i].data()[c] = buf[i];
			}
		};
		// Op must see each out element whole, exactly once, so compose
		// non-assign ops through an element buffer.
		if constexpr (std::same_as<Op, assign_op>) {
			emit(out);
		} else {
			auto vbuf = buffer_pool<V>::get(sz(out));
			emit(vbuf.span());
			for (int i = 0; i < sz(out); i++) op(out[i], vbuf.span()[i]);
		}
	}
};

// Convolve mat<N> (NxN matrices), with accumulation in product space
template <engine E, int N>
struct matrix : componentwise<E, mat<typename E::value_type, N>, N * N> {
	using base = componentwise<E, mat<typename E::value_type, N>, N * N>;
	static constexpr bool commutative = false;
	static constexpr int unit_scale = base::unit_scale;
	template <int A = unit_scale> using transformed_t = typename base::template transformed_t<A>;
	template <int K> using product_t = typename base::template product_t<K * N>;
	using transformed = typename base::transformed;
	using product = product_t<unit_scale * unit_scale>;

	// right fold over k so a tracked inner engine's per-addend types line up
	template <int A, int B, int k = 0>
	static auto entry(const transformed_t<A>& a, const transformed_t<B>& b, int r, int c, int n) {
		auto e = E::mul(a.t[size_t(r) * N + k], b.t[size_t(k) * N + c], n);
		if constexpr (k + 1 == N) return e;
		else return E::add(std::move(e), entry<A, B, k + 1>(a, b, r, c, n));
	}
	template <int A, int B>
	static product_t<A * B> mul(const transformed_t<A>& a, const transformed_t<B>& b, int n) {
		product_t<A * B> p;
		for (int r = 0; r < N; r++) for (int c = 0; c < N; c++)
			p.t[size_t(r) * N + c] = entry<A, B>(a, b, r, c, n);
		return p;
	}
	template <int A> static auto sq(const transformed_t<A>& a, int n) { return mul(a, a, n); }
	template <int A1, int B1, int A2, int B2, int k = 0>
	static auto entry2(
		const transformed_t<A1>& a1, const transformed_t<B1>& b1,
		const transformed_t<A2>& a2, const transformed_t<B2>& b2,
		int r, int c, int n
	) {
		auto e = E::mul2(
			a1.t[size_t(r) * N + k], b1.t[size_t(k) * N + c],
			a2.t[size_t(r) * N + k], b2.t[size_t(k) * N + c],
			n
		);
		if constexpr (k + 1 == N) return e;
		else return E::add(std::move(e), entry2<A1, B1, A2, B2, k + 1>(a1, b1, a2, b2, r, c, n));
	}
	template <int A1, int B1, int A2, int B2>
	static product_t<A1 * B1 + A2 * B2> mul2(
		const transformed_t<A1>& a1, const transformed_t<B1>& b1,
		const transformed_t<A2>& a2, const transformed_t<B2>& b2,
		int n
	) {
		product_t<A1 * B1 + A2 * B2> p;
		for (int r = 0; r < N; r++) for (int c = 0; c < N; c++)
			p.t[size_t(r) * N + c] = entry2<A1, B1, A2, B2>(a1, b1, a2, b2, r, c, n);
		return p;
	}
};

// Convolve trunc_series<num, N> (power series truncated at N), with accumulation in product space
template <engine E, int N>
struct trunc : componentwise<E, trunc_series<typename E::value_type, N>, N> {
	using base = componentwise<E, trunc_series<typename E::value_type, N>, N>;
	static constexpr bool commutative = E::commutative;
	static constexpr int unit_scale = base::unit_scale;
	template <int A = unit_scale> using transformed_t = typename base::template transformed_t<A>;
	template <int K> using product_t = typename base::template product_t<K * N>;
	using transformed = typename base::transformed;
	using product = product_t<unit_scale * unit_scale>;

	template <int A, int B, int s, int i = 0>
	static auto entry(const transformed_t<A>& a, const transformed_t<B>& b, int n) {
		auto e = E::mul(a.t[size_t(i)], b.t[size_t(s - i)], n);
		if constexpr (i == s) return e;
		else return E::add(std::move(e), entry<A, B, s, i + 1>(a, b, n));
	}
	template <int A, int B>
	static product_t<A * B> mul(const transformed_t<A>& a, const transformed_t<B>& b, int n) {
		product_t<A * B> p;
		[&]<size_t... s_>(std::index_sequence<s_...>) {
			((p.t[s_] = entry<A, B, int(s_)>(a, b, n)), ...);
		}(std::make_index_sequence<size_t(N)>{});
		return p;
	}
	template <int A> static auto sq(const transformed_t<A>& a, int n) { return mul(a, a, n); }
	template <int A1, int B1, int A2, int B2, int s, int i = 0>
	static auto entry2(
		const transformed_t<A1>& a1, const transformed_t<B1>& b1,
		const transformed_t<A2>& a2, const transformed_t<B2>& b2,
		int n
	) {
		auto e = E::mul2(a1.t[size_t(i)], b1.t[size_t(s - i)], a2.t[size_t(i)], b2.t[size_t(s - i)], n);
		if constexpr (i == s) return e;
		else return E::add(std::move(e), entry2<A1, B1, A2, B2, s, i + 1>(a1, b1, a2, b2, n));
	}
	template <int A1, int B1, int A2, int B2>
	static product_t<A1 * B1 + A2 * B2> mul2(
		const transformed_t<A1>& a1, const transformed_t<B1>& b1,
		const transformed_t<A2>& a2, const transformed_t<B2>& b2,
		int n
	) {
		product_t<A1 * B1 + A2 * B2> p;
		[&]<size_t... s_>(std::index_sequence<s_...>) {
			((p.t[s_] = entry2<A1, B1, A2, B2, int(s_)>(a1, b1, a2, b2, n)), ...);
		}(std::make_index_sequence<size_t(N)>{});
		return p;
	}
};

// Stable variants of the wrapper engines: do not accumulate in product space.
// This costs an extra log factor.

template <int N> constexpr std::array<int, size_t(N) * N + 1> matrix_stable_ofs = [] {
	std::array<int, size_t(N) * N + 1> r{};
	for (int i = 0; i <= N * N; i++) r[size_t(i)] = i * N;
	return r;
}();

template <engine E, int N>
struct matrix_stable
		: componentwise<E, mat<typename E::value_type, N>, N * N, matrix_stable_ofs<N>> {
	using base = componentwise<E, mat<typename E::value_type, N>, N * N, matrix_stable_ofs<N>>;
	static constexpr bool commutative = false;
	static constexpr int unit_scale = base::unit_scale;
	template <int A = unit_scale> using transformed_t = typename base::template transformed_t<A>;
	template <int K> using product_t = typename base::template product_t<K>;
	using transformed = typename base::transformed;
	using product = product_t<unit_scale * unit_scale>;

	template <int A, int B>
	static product_t<A * B> mul(const transformed_t<A>& a, const transformed_t<B>& b, int n) {
		product_t<A * B> p;
		// entry (r, c)'s k-th addend a(r,k)*b(k,c), grouped per the offsets
		for (int r = 0; r < N; r++) for (int c = 0; c < N; c++) for (int k = 0; k < N; k++)
			p.t[(size_t(r) * N + c) * N + k] = E::mul(a.t[size_t(r) * N + k], b.t[size_t(k) * N + c], n);
		return p;
	}
	template <int A> static auto sq(const transformed_t<A>& a, int n) { return mul(a, a, n); }
	template <int A1, int B1, int A2, int B2>
	static product_t<A1 * B1 + A2 * B2> mul2(
		const transformed_t<A1>& a1, const transformed_t<B1>& b1,
		const transformed_t<A2>& a2, const transformed_t<B2>& b2,
		int n
	) {
		product_t<A1 * B1 + A2 * B2> p;
		for (int r = 0; r < N; r++) for (int c = 0; c < N; c++) for (int k = 0; k < N; k++)
			p.t[(size_t(r) * N + c) * N + k] = E::mul2(
				a1.t[size_t(r) * N + k], b1.t[size_t(k) * N + c],
				a2.t[size_t(r) * N + k], b2.t[size_t(k) * N + c],
				n
			);
		return p;
	}
};

template <int N> constexpr std::array<int, size_t(N) + 1> trunc_series_stable_ofs = [] {
	std::array<int, size_t(N) + 1> r{};
	for (int i = 0; i <= N; i++) r[size_t(i)] = i * (i + 1) / 2;
	return r;
}();

template <engine E, int N>
struct trunc_stable
		: componentwise<E, trunc_series<typename E::value_type, N>, N, trunc_series_stable_ofs<N>> {
	using base = componentwise<E, trunc_series<typename E::value_type, N>, N, trunc_series_stable_ofs<N>>;
	static constexpr bool commutative = E::commutative;
	static constexpr int unit_scale = base::unit_scale;
	template <int A = unit_scale> using transformed_t = typename base::template transformed_t<A>;
	template <int K> using product_t = typename base::template product_t<K>;
	using transformed = typename base::transformed;
	using product = product_t<unit_scale * unit_scale>;

	template <int A, int B>
	static product_t<A * B> mul(const transformed_t<A>& a, const transformed_t<B>& b, int n) {
		product_t<A * B> p;
		for (int s = 0; s < N; s++) for (int i = 0; i <= s; i++)
			p.t[size_t(trunc_series_stable_ofs<N>[size_t(s)] + i)] = E::mul(a.t[size_t(i)], b.t[size_t(s - i)], n);
		return p;
	}
	template <int A> static auto sq(const transformed_t<A>& a, int n) { return mul(a, a, n); }
	template <int A1, int B1, int A2, int B2>
	static product_t<A1 * B1 + A2 * B2> mul2(
		const transformed_t<A1>& a1, const transformed_t<B1>& b1,
		const transformed_t<A2>& a2, const transformed_t<B2>& b2,
		int n
	) {
		product_t<A1 * B1 + A2 * B2> p;
		for (int s = 0; s < N; s++) for (int i = 0; i <= s; i++)
			p.t[size_t(trunc_series_stable_ofs<N>[size_t(s)] + i)] = E::mul2(
				a1.t[size_t(i)], b1.t[size_t(s - i)],
				a2.t[size_t(i)], b2.t[size_t(s - i)],
				n
			);
		return p;
	}
};

static_assert(engine<ntt<modnum<998244353>>>);
static_assert(engine<ntt<mod_goldilocks>>);
static_assert(engine<real<double>>);
static_assert(engine<split<modnum<int(1e9)+7>>>);
static_assert(engine<crt<modnum<int(1e9)+7>>>);
static_assert(engine<matrix<ntt<modnum<998244353>>, 2>>);
static_assert(engine<trunc<ntt<modnum<998244353>>, 3>>);
// tracked inner engines work when the accumulated scale fits the budget (N <= 2)
static_assert(engine<matrix<split<modnum<int(1e9)+7>>, 2>>);
static_assert(engine<trunc<crt<modnum<int(1e9)+7>>, 2>>);
// the stable variants keep tracked inner engines sound at any N
static_assert(engine<matrix_stable<split<modnum<int(1e9)+7>>, 3>>);
static_assert(engine<trunc_stable<crt<modnum<int(1e9)+7>>, 3>>);

/* namespace engines */ }

// short spelling for E::transformed at use sites
template <engine E> using transformed = typename E::transformed;

// ==== multiply layer ====
// These are free functions to convolve spans.
//
// The interfaces will typically take input spans, an output span, and an Op representing how to fold the result into the output.
// Output spans may alias one of the input spans.
// Output spans may be shorter than expected; the output will just be truncated.
//
// Some functions may also take E::transformed& objects associated with the input
// spans. These will be lazily filled (see E::extend_to) and used if available.

// Circular convolution mod n (power of 2)
template <engine E, typename Op = assign_op>
void multiply_circular(std::span<const typename E::value_type> a, std::span<const typename E::value_type> b,
		std::span<typename E::value_type> out, int n, Op op = {}) {
	assert(!(n & (n-1)));
	auto ta = E::transform(a, n);
	auto tb = E::transform(b, n);
	E::finish(E::mul(ta, tb, n), out, op);
}

template <engine E, typename Op = assign_op>
void square_circular(std::span<const typename E::value_type> a, std::span<typename E::value_type> out, int n, Op op = {}) {
	assert(!(n & (n-1)));
	auto ta = E::transform(a, n);
	E::finish(E::sq(ta, n), out, op);
}

namespace detail {
// Arrays of length 2^k + 1 are somewhat common, so we will optimize them by
// multiplying mod 2^k, and fixing up the leading coefficient.

// Helpers to detect and perform this optimization.
struct conv_size { int n; bool cut; };
inline conv_size conv_size_for(int s) {
	int n = nextPow2(s);
	bool cut = (n == 2 * (s - 1));
	return {cut ? n / 2 : n, cut};
}

// Call op while lazily applying the correction if necessary
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

// Applies op, diverting the wrapped leading coefficient of a cut product:
// out[0] receives c0 and the wraparound term is captured into cn for the
// caller to emit at out[s-1].
template <typename T, typename Op>
struct cut_op {
	Op op;
	T* out0;
	T c0;
	T& cn;
	void operator()(T& x, T v) const {
		if (&x == out0) { cn = v - c0; v = c0; }
		op(x, v);
	}
};

// finish + emit_linear fused: write the finished product directly into out,
// applying the cut correction in place.
template <engine E, typename P, typename Op = assign_op>
void finish_linear(
	P&& p, int n, int s, bool cut,
	typename E::value_type c0, std::span<typename E::value_type> out, Op op = {}
) {
	using T = typename E::value_type;
	if (sz(out) == 0) return;
	int lim = min(sz(out), min(s, n));
	if (!cut) {
		E::finish(std::move(p), out.subspan(0, lim), op);
	} else {
		T cn{};
		E::finish(std::move(p), out.subspan(0, lim), cut_op<T, Op>{op, &out[0], c0, cn});
		if (sz(out) >= s) op(out[s-1], cn);
	}
}

}

template <engine E, typename Op = assign_op>
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

template <engine E, typename Op = assign_op>
void multiply(std::span<const typename E::value_type> a, transformed<E>& ta,
		std::span<const typename E::value_type> b, transformed<E>& tb,
		std::span<typename E::value_type> out, Op op = {}) {
	using T = typename E::value_type;
	if (sz(a) == 0 || sz(b) == 0) return;
	int s = sz(a) + sz(b) - 1;
	auto [n, cut] = detail::conv_size_for(s);
	T c0 = a[0] * b[0];
	E::extend_to(ta, n, a);
	E::extend_to(tb, n, b);
	detail::finish_linear<E>(E::mul(ta, tb, n), n, s, cut, c0, out, op);
}

template <engine E, typename Op = assign_op>
void multiply_add2(std::span<const typename E::value_type> a1, transformed<E>& ta1,
		std::span<const typename E::value_type> b1, transformed<E>& tb1,
		std::span<const typename E::value_type> a2, transformed<E>& ta2,
		std::span<const typename E::value_type> b2, transformed<E>& tb2,
		std::span<typename E::value_type> out, Op op = {}) {
	using T = typename E::value_type;
	assert(sz(a1) > 0 && sz(b1) > 0 && sz(a2) > 0 && sz(b2) > 0);
	int s = sz(a1) + sz(b1) - 1;
	assert(sz(a2) + sz(b2) - 1 == s);
	auto [n, cut] = detail::conv_size_for(s);
	T c0 = a1[0] * b1[0] + a2[0] * b2[0];
	E::extend_to(ta1, n, a1); E::extend_to(tb1, n, b1);
	E::extend_to(ta2, n, a2); E::extend_to(tb2, n, b2);
	detail::finish_linear<E>(E::mul2(ta1, tb1, ta2, tb2, n), n, s, cut, c0, out, op);
}

// As multiply_add2, but also outputs the summed pointwise product as a reusable
// transform of the (full-length) result, like multiply_cached.
template <engine E>
void multiply_add2_cached(
		std::span<const typename E::value_type> a1, transformed<E>& ta1,
		std::span<const typename E::value_type> b1, transformed<E>& tb1,
		std::span<const typename E::value_type> a2, transformed<E>& ta2,
		std::span<const typename E::value_type> b2, transformed<E>& tb2,
		std::vector<typename E::value_type>& coeffs, transformed<E>& t) {
	using T = typename E::value_type;
	assert(sz(a1) > 0 && sz(b1) > 0 && sz(a2) > 0 && sz(b2) > 0);
	int s = sz(a1) + sz(b1) - 1;
	assert(sz(a2) + sz(b2) - 1 == s);
	coeffs.assign(size_t(s), T{});
	t = transformed<E>{};
	if constexpr (std::same_as<typename E::product, transformed<E>>) {
		auto [n, cut] = detail::conv_size_for(s);
		T c0 = a1[0] * b1[0] + a2[0] * b2[0];
		E::extend_to(ta1, n, a1); E::extend_to(tb1, n, b1);
		E::extend_to(ta2, n, a2); E::extend_to(tb2, n, b2);
		auto p = E::mul2(ta1, tb1, ta2, tb2, n);
		auto tp = p;
		detail::finish_linear<E>(std::move(p), n, s, cut, c0, std::span<T>(coeffs));
		t = std::move(tp);
	} else {
		multiply_add2<E>(a1, ta1, b1, tb1, a2, ta2, b2, tb2, std::span<T>(coeffs));
	}
}

// This helper also accepts an output transform which will be populated if it is cheap to do so
template <engine E>
void multiply_cached(std::span<const typename E::value_type> a, transformed<E>& ta,
		std::span<const typename E::value_type> b, transformed<E>& tb,
		std::vector<typename E::value_type>& coeffs, transformed<E>& t) {
	using T = typename E::value_type;
	coeffs.assign(size_t(sz(a) && sz(b) ? sz(a) + sz(b) - 1 : 0), T{});
	t = transformed<E>{};
	if (coeffs.empty()) return;
	int s = sz(coeffs);
	if constexpr (std::same_as<typename E::product, transformed<E>>) {
		auto [n, cut] = detail::conv_size_for(s);
		T c0 = a[0] * b[0];
		E::extend_to(ta, n, a);
		E::extend_to(tb, n, b);
		auto p = E::mul(ta, tb, n);
		auto tp = p;
		detail::finish_linear<E>(std::move(p), n, s, cut, c0, std::span<T>(coeffs));
		t = std::move(tp);
	} else {
		multiply<E>(a, ta, b, tb, std::span<T>(coeffs));
	}
}

template <engine E, typename Op = assign_op>
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

template <engine E, typename Op = assign_op>
void square(std::span<const typename E::value_type> a, transformed<E>& ta,
		std::span<typename E::value_type> out, Op op = {}) {
	using T = typename E::value_type;
	if (sz(a) == 0) return;
	int s = 2 * sz(a) - 1;
	auto [n, cut] = detail::conv_size_for(s);
	T c0 = a[0] * a[0];
	E::extend_to(ta, n, a);
	detail::finish_linear<E>(E::sq(ta, n), n, s, cut, c0, out, op);
}

// As square, but also outputs the pointwise product as a reusable transform of
// the result (empty when the engine's product isn't a transform).
template <engine E>
void square_cached(std::span<const typename E::value_type> a, transformed<E>& ta,
		std::vector<typename E::value_type>& coeffs, transformed<E>& t) {
	using T = typename E::value_type;
	coeffs.assign(size_t(sz(a) ? 2 * sz(a) - 1 : 0), T{});
	t = transformed<E>{};
	if (coeffs.empty()) return;
	int s = sz(coeffs);
	if constexpr (std::same_as<typename E::product, transformed<E>>) {
		auto [n, cut] = detail::conv_size_for(s);
		T c0 = a[0] * a[0];
		E::extend_to(ta, n, a);
		auto p = E::sq(ta, n);
		auto tp = p;
		detail::finish_linear<E>(std::move(p), n, s, cut, c0, std::span<T>(coeffs));
		t = std::move(tp);
	} else {
		square<E>(a, ta, std::span<T>(coeffs));
	}
}

template <engine E> vector<typename E::value_type> multiply(
		const vector<typename E::value_type>& a, const vector<typename E::value_type>& b) {
	using T = typename E::value_type;
	if (sz(a) == 0 || sz(b) == 0) return {};
	vector<T> r(sz(a) + sz(b) - 1);
	multiply<E>(std::span<const T>(a), std::span<const T>(b), std::span<T>(r));
	return r;
}

template <engine E> vector<typename E::value_type> square(const vector<typename E::value_type>& a) {
	using T = typename E::value_type;
	if (sz(a) == 0) return {};
	vector<T> r(2 * sz(a) - 1);
	square<E>(std::span<const T>(a), std::span<T>(r));
	return r;
}

namespace detail {
// emit_linear but for middle_product
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

// Middle product (the transposed multiplication): takes only coefficients of a * b which include terms from all of b.
// Must have len(a) >= len(b)
template <engine E, typename Op = assign_op>
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

// TODO: Let's decide whether to keep vector<> returning forms or not; this
// largely depends on whether we think these functions are a public interface or
// merely convenience for value type implementors.
template <engine E> vector<typename E::value_type> middle_product(
		std::span<const typename E::value_type> a, std::span<const typename E::value_type> b) {
	using T = typename E::value_type;
	if (sz(a) == 0 || sz(b) == 0) return {};
	assert(sz(a) >= sz(b));
	vector<T> r(size_t(sz(a) - sz(b) + 1));
	middle_product<E>(a, b, std::span<T>(r));
	return r;
}

template <engine E, typename Op = assign_op>
void middle_product(std::span<const typename E::value_type> a, transformed<E>& ta,
		std::span<const typename E::value_type> b, transformed<E>& tb,
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
	E::extend_to(ta, n, a);
	E::extend_to(tb, n, b);
	auto buf = buffer_pool<T>::get(n);
	E::finish(E::mul(ta, tb, n), buf.span());
	detail::emit_middle<T>(buf.span(), cut, sz(a), sz(b),
			a[0] * b[0], a[sz(a) - 1] * b[sz(b) - 1], out, op);
}

template <engine E>
vector<typename E::value_type> middle_product(std::span<const typename E::value_type> a, transformed<E>& ta,
		std::span<const typename E::value_type> b, transformed<E>& tb) {
	using T = typename E::value_type;
	if (sz(a) == 0 || sz(b) == 0) return {};
	assert(sz(a) >= sz(b));
	vector<T> r(size_t(sz(a) - sz(b) + 1));
	middle_product<E>(a, ta, b, tb, std::span<T>(r));
	return r;
}

/* namespace fft */ }

// ==== value types ====

// Helper packed bivariate buffer for Kinoshita-Li composition (arXiv:2404.05177).
//
// The motivation is performing Bostan-Mori (Graeffe root-squaring) to compute
// something like [x^n] P / Q_0(x, y) with deg_y(Q_0) = 1 and deg_x(Q_0) = n.
//
// In each step, we want to compute Q_{i+1}(x^2, y) = Q_i(x, y) * Q_i(-x, y).
// This doubles the degree of y and also lets us truncate x at half the previous
// degree, leaving the total size invariant.
//
// We will store Q as a packed buffer with x as the inner dimension to facilitate easy Q(-x) substitution.
// The inner span will be 2*deg(x), and the outer span will be 2*deg(y).
// As we advance, we will also return the cached transform of Q_i(-x, y) for the caller to use in the numerator.
template <fft::engine E> struct packed_bivariate {
	using T = typename E::value_type;
	int L, l;
	std::vector<T> c;

	// Q_0 = 1 - y g(x), deg g < n <= 2^L
	packed_bivariate(int L_, std::span<const T> g) : L(L_), l(0), c(size_t(4) << L) {
		c[0] = T(1);
		for (int i = 0; i < sz(g); i++) c[(2 << L) + i] = -g[i];
	}

	fft::transformed<E> advance() {
		int B = 4 << L;
		auto tq = E::transform(std::span<const T>(c), B);
		auto tn = E::negate_arg(tq, B);
		E::finish(E::mul(tq, tn, B), std::span<T>(c));
		l++;
		// compactify x^2 -> x
		for (int i = 1; i < (2 << L); i++) c[i] = c[2*i];
		// undo the circular wraparound using monicity in y
		for (int i = 0; i < (2 << (L - l)); i++) {
			c[(2 << L) + i] = c[i];
			c[i] = T(0);
		}
		c[2 << L] -= T(1);
		c[0] = T(1);
		// zero x coefficients beyond the level's truncation mod x^(2^(L-l))
		std::fill(c.begin() + (2 << L) + (1 << (L - l)), c.end(), T(0));
		for (int i = 0; i < (2 << L); i += 2 << (L - l)) {
			for (int j = 0; j < (1 << (L - l)); j++) {
				c[i + (1 << (L - l)) + j] = T(0);
			}
		}
		return tn;
	}
};

template <typename A, typename B>
concept same_engine = std::same_as<typename A::engine_t, typename B::engine_t>;

namespace series {

// Non-owning view of power series coefficients: the span pattern (contiguous
// window + series semantics), borrowed from an owning series-like type.
template <fft::engine E, bool exact = false>
struct span {
	using T = typename E::value_type;
	using engine_t = E;
	static constexpr bool exact_v = exact;

	span() = default;
	explicit span(std::span<const T> s_) : s(s_) {}

	int len() const { return sz(s); }
	const T& operator[](int i) const { return s[size_t(i)]; }
	auto begin() const { return s.begin(); }
	auto end() const { return s.end(); }
	std::span<const T> coeffs() const { return s; }
	span underlying() const { return *this; }
	span first(int n) const { return span(s.first(size_t(n))); }

private:
	std::span<const T> s;
};

// `vec` represents both exact (finite) power series (R[x]) and prefixes of infinite power series (R[[x]]), depending on the flag.
// `exact` and `trunc` are aliases.
//
// Operators here are typically permissive: they will accept combinations of unequal types and lengths.
template <fft::engine E, bool exact = false>
struct vec : public std::vector<typename E::value_type> {
	using T = typename E::value_type;
	using engine_t = E;
	static constexpr bool exact_v = exact;
	using std::vector<T>::vector;

	const vec& underlying() const { return *this; }
	// a free const borrow of the coefficients: implicit
	operator span<E, exact>() const {
		return span<E, exact>(std::span<const T>(*this));
	}

	// exact -> trunc is implicit, trunc -> exact is explicit
	template <bool oe> requires (oe && !exact)
	vec(const vec<E, oe>& p) : std::vector<T>(p) {}
	template <bool oe> requires (oe && !exact)
	vec(vec<E, oe>&& p) : std::vector<T>(std::move(p)) {}
	template <bool oe> requires (!oe && exact)
	explicit vec(const vec<E, oe>& p) : std::vector<T>(p) {}
	template <bool oe> requires (!oe && exact)
	explicit vec(vec<E, oe>&& p) : std::vector<T>(std::move(p)) {}

	// adopt a plain coefficient vector
	explicit vec(std::vector<T> v) : std::vector<T>(std::move(v)) {}

	int len() const {
		return int(this->size());
	}
	int degree() const requires (exact) {
		return len() - 1;
	}
	void extend(int sz) {
		assert(sz >= len());
		this->resize(sz);
	}
	void shrink(int sz) {
		assert(sz <= len());
		this->resize(sz);
	}
	// multiply by x^n within the fixed precision window
	void shift_trunc(int n = 1) requires (!exact) {
		assert(n >= 0 && n <= len());
		std::rotate(this->begin(), this->end()-n, this->end());
		std::fill(this->begin(), this->begin()+n, T(0));
	}
	// divide by x^n and 0-pad within the fixed precision window
	void unshift_trunc(int n = 1) requires (!exact) {
		assert(n >= 0 && n <= len());
		std::fill(this->begin(), this->begin()+n, T(0));
		std::rotate(this->begin(), this->begin()+n, this->end());
	}

	// in-place forms require that the result's exactness/length must equal this operand's
	template <bool oe> requires (oe || !exact)
	vec& operator += (const vec<E, oe>& o) {
		if constexpr (exact) { if (o.len() > len()) this->resize(o.len()); }
		else if constexpr (!oe) { if (o.len() < len()) this->resize(o.len()); }
		for (int i = 0; i < std::min(len(), o.len()); i++) {
			(*this)[i] += o[i];
		}
		return *this;
	}
	template <bool oe> requires (oe || !exact)
	vec& operator -= (const vec<E, oe>& o) {
		if constexpr (exact) { if (o.len() > len()) this->resize(o.len()); }
		else if constexpr (!oe) { if (o.len() < len()) this->resize(o.len()); }
		for (int i = 0; i < std::min(len(), o.len()); i++) {
			(*this)[i] -= o[i];
		}
		return *this;
	}

	vec& operator *= (const T& n) {
		for (auto& v : *this) v *= n;
		return *this;
	}
	friend vec operator * (const vec& a, const T& n) {
		vec r(a.size());
		for (int i = 0; i < a.len(); i++) {
			r[i] = a[i] * n;
		}
		return r;
	}
	friend vec operator * (const T& n, const vec& a) {
		vec r(a.size());
		for (int i = 0; i < a.len(); i++) {
			r[i] = n * a[i];
		}
		return r;
	}


	vec& operator *= (const vec& o) {
		return *this = (*this) * o;
	}

	// Newton inversion: 1/a mod x^a.len(). Generic over any engine; per doubling step
	// n -> m = 2n this is 5 transforms of size m, reusing b's transform for both circular
	// products; in each product the wraparound only contaminates coefficients [0, n)
	// which are already known.
	//
	// This is correct for non-commutative rings.
	friend vec inverse(const vec& a) requires (!exact) {
		vec r(a.size());
		if (a.len() == 0) return r;
		int s = nextPow2(a.len());
		std::vector<T> b(size_t(s), T{});
		b[0] = inv(a[0]);
		for (int n = 1; n < a.len(); n *= 2) {
			int m = 2 * n;
			auto ta = E::transform(std::span<const T>(a).first(std::min(a.len(), m)), m);
			auto tb = E::transform(std::span<const T>(b).first(n), m);
			// e = a*b mod x^m; only e[n..m) is needed (and is wraparound-free).
			auto e = fft::buffer_pool<T>::get(m);
			E::finish(E::mul(ta, tb, m), e.span());
			for (int i = 0; i < n; i++) e[i] = T{};
			auto te = E::transform(std::span<const T>(e.span()), m);
			auto c = fft::buffer_pool<T>::get(m);
			// b' = 2b - b*(a*b): keep b on the left of e = a*b
			E::finish(E::mul(tb, te, m), c.span());
			for (int i = n; i < std::min(m, a.len()); i++) b[i] = -c[i];
		}
		std::copy(b.begin(), b.begin() + a.len(), r.begin());
		return r;
	}
	// TODO: operator / can be done slightly faster than inverse:
	// we only need the n/2 terms of inverse(), and can do the last Newton step directly on the quotient

	friend vec stretch(const vec& a, int n) {
		vec r(a.size());
		for (int i = 0; i*n < int(a.size()); i++) {
			r[i*n] = a[i];
		}
		return r;
	}
	friend vec deriv_shift(vec a) {
		for (int i = 0; i < a.len(); i++) {
			a[i] *= i;
		}
		return a;
	}
	friend vec integ_shift(vec a) {
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
	friend vec integ_shift_offset(vec a, int offset) {
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
	friend vec deriv_shift_log(vec a) requires (!exact) {
		auto r = deriv_shift(a);
		return r * inverse(a);
	}
	friend vec poly_log(vec a) requires (!exact) {
		assert(a[0] == 1);
		return integ_shift(deriv_shift_log(std::move(a)));
	}
	friend vec poly_exp(vec a) requires (!exact) {
		// See https://mathexp.eu/bostan/publications/BoSc09a.pdf for details
		assert(a.size() >= 1);
		assert(a[0] == 0);
		vec r(1, T(1)); r.reserve(a.size());
		vec invR(1, T(1)); invR.reserve(a.size());
		while (r.size() < a.size()) {
			int o_sz = int(r.size());
			int n_sz = std::min(o_sz * 2, int(a.size()));
			vec t = deriv_shift(vec(a.begin(), a.begin() + o_sz));
			fft::multiply_circular<E>(std::span<const T>(t), std::span<const T>(r).first(o_sz), std::span<T>(t), o_sz);
			t = deriv_shift(r) - t;
			t *= invR;
			t.resize(n_sz - o_sz);
			vec v(a.begin() + o_sz, a.begin() + n_sz);
			v -= integ_shift_offset(t, o_sz);
			v *= r;
			r.resize(n_sz);
			std::copy(v.begin(), v.end(), r.begin() + o_sz);
			if (r.size() < a.size()) {
				// double invR via a Newton step
				assert(int(r.size()) == 2 * int(invR.size()));
				int n = int(invR.size());
				int nn = int(r.size());
				vec tmp(4 * n);
				fft::square<E>(std::span<const T>(invR).first(n), std::span<T>(tmp));
				fft::multiply<E>(std::span<const T>(tmp).first(nn), std::span<const T>(r).first(nn), std::span<T>(tmp));
				invR.resize(nn);
				for (int i = n; i < nn; i++) invR[i] = -tmp[i];
			}
		}
		return r;
	}
	friend vec poly_pow_monic(vec a, T k) requires (!exact) {
		if (a.empty()) return a;
		assert(a.size() >= 1);
		assert(a[0] == 1);
		a = poly_log(a);
		a *= k;
		return poly_exp(a);
	}
	friend vec poly_pow(vec a, int64_t k) requires (!exact) {
		assert(k >= 0);
		if (k == 0) {
			vec r(a.len(), T(0));
			if (r.len() > 0) r[0] = T(1);
			return r;
		}

		int st = 0;
		while (st < a.len() && a[st] == 0) st++;

		if (st > 0 && k > (a.len() - 1) / st) {
			return vec(a.len(), T(0));
		}

		vec r(a.begin() + st, a.end() - (st * (k-1)));
		T leading_coeff = r[0];
		r *= inv(leading_coeff);
		r = poly_pow_monic(r, T(k));
		r *= power(leading_coeff, k);
		r.insert(r.begin(), st * k, T(0));
		assert(r.len() == a.len());
		return r;
	}

	friend vec to_newton_sums(const vec& a, int deg) requires (!exact) {
		auto r = deriv_shift_log(a);
		r[0] = deg;
		for (int i = 1; i < int(r.size()); i++) r[i] = -r[i];
		return r;
	}
	friend vec from_newton_sums(vec S, int deg) requires (!exact) {
		assert(S[0] == int(deg));
		S[0] = 0;
		for (int i = 1; i < int(S.size()); i++) S[i] = -S[i];
		return poly_exp(integ_shift(std::move(S)));
	}

	// Calculates prod 1/(1-x^i)^{a[i]}
	friend vec euler_transform(const vec& a) requires (!exact) {
		vec r = deriv_shift(a);
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
	friend vec inverse_euler_transform(const vec& a) requires (!exact) {
		vec r = deriv_shift(poly_log(a));
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
	friend vec poly_compose(const vec& f, const vec& g) requires (!exact) {
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
		int B = 4 << L;
		packed_bivariate<E> Q(L, std::span<const T>(g));
		// tneg[l] is the transform of Q_l(-x, y), reused by the pushdown pass below
		std::vector<fft::transformed<E>> tneg;
		tneg.reserve(L);
		for (int l = 1; l <= L; l++) tneg.push_back(Q.advance());
		vec P;
		{
			P = f;
			std::reverse(P.begin(), P.end());
			vec QL((1 << L) + 1);
			for (int i = 0; i <= (1 << L); i++) {
				QL[i] = Q.c[2 * i];
			}
			QL.resize(m, T(0));
			P *= inverse(QL);
			std::reverse(P.begin(), P.end());
			P.resize(1 << L, T(0));
			std::reverse(P.begin(), P.end());
			P.resize(B, T(0));
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
			auto tp = E::transform(std::span<const T>(P), B);
			E::finish(E::mul(tneg[l], tp, B), std::span<T>(P));
			for (int i = 0; i < (2 << L); i++) {
				P[i] = P[(2 << L) + i];
				P[(2 << L) + i] = T(0);
			}
		}
		return vec(P.begin(), P.begin() + n);
	}
};

template <fft::engine E> using exact = vec<E, true>;
template <fft::engine E> using trunc = vec<E, false>;

// Series-like concepts: the binary operators below are written once as constrained
// templates and dispatch on which memoized transforms an operand carries.
// A series-like type exposes its engine/exactness and its coefficients as a
// span borrow; cached wrappers additionally expose their transform
// caches (filling them is logically const).
template <typename S>
concept like = fft::engine<typename S::engine_t> && requires(const S& s) {
	{ S::exact_v } -> std::convertible_to<bool>;
	{ s.len() } -> std::same_as<int>;
	{ s.underlying() } -> std::convertible_to<span<typename S::engine_t, S::exact_v>>;
};
// carries one extendable transform of the whole coefficient sequence
template <typename S>
concept has_cache = like<S> && requires(const S& s) {
	{ s.cache() } -> std::same_as<fft::transformed<typename S::engine_t>&>;
};

// A borrowed series paired with the transform serving it: the
// normalized operand form fed to the cached fft:: entry points. Models has_cache.
template <fft::engine E, bool exact>
struct cached_span {
	using engine_t = E;
	static constexpr bool exact_v = exact;

	span<E, exact> s;
	std::reference_wrapper<fft::transformed<E>> f;

	int len() const { return s.len(); }
	span<E, exact> underlying() const { return s; }
	fft::transformed<E>& cache() const { return f; }
};

// carries transforms of power-of-two prefixes, usable at any precision:
// prefix(n) borrows the length-min(n, len) prefix with its cache.
// Trunc-only: an exact operand participates whole, so has_cache covers it.
template <typename S>
concept has_prefix_cache = like<S> && !S::exact_v && requires(const S& s, int n) {
	{ s.prefix(n) } -> has_cache;
};

// Wrapper around vec which caches the transform of the whole series.
// Ops exploit the cache whenever the whole span participates; a trunc series'
// whole-sequence transform is still useful for middle products and repeated
// full-precision use.
template <fft::engine E, bool exact = true>
struct cached {
	using T = typename E::value_type;
	using engine_t = E;
	static constexpr bool exact_v = exact;

	cached() = default;
	// moving coefficients in or out is free: implicit on rvalues, explicit copy otherwise
	cached(vec<E, exact>&& s_) : s(std::move(s_)) {}
	explicit cached(const vec<E, exact>& s_) : s(s_) {}
	operator vec<E, exact>() && { return std::move(s); }

	int len() const { return s.len(); }
	const vec<E, exact>& underlying() const { return s; }
	const T& operator[](int i) const { return s[size_t(i)]; }
	auto begin() const { return s.cbegin(); }
	auto end() const { return s.cend(); }
	// the transform of underlying(), fed to the cached fft:: entry points alongside it
	fft::transformed<E>& cache() const { return f; }

	template <like S>
	friend bool operator==(const cached& a, const S& b) {
		span<E, S::exact_v> bs = b.underlying();
		return a.len() == bs.len() && std::equal(a.s.begin(), a.s.end(), bs.begin());
	}

private:
	vec<E, exact> s;
	mutable fft::transformed<E> f; // memoized transform: filling it is logically const
};

namespace detail {
// the operand's whole cache if it carries one, else the caller's throwaway cache
template <like S>
fft::transformed<typename S::engine_t>& whole_cache_or(const S& s, fft::transformed<typename S::engine_t>& tmp) {
	if constexpr (has_cache<S>) return s.cache(); else return tmp;
}
/* namespace detail */ }

// Both consume whole-sequence transforms by nature (the full span always
// participates), so only whole caches apply, never prefix caches.
template <like A>
auto square(const A& a) {
	using E = typename A::engine_t;
	using T = typename E::value_type;
	span<E, A::exact_v> av = a.underlying();
	fft::transformed<E> ta_;
	if constexpr (A::exact_v) {
		// like operator*, an exact square returns has_cache, adopting the
		// pointwise product as the result's transform when the engine supports it
		std::vector<T> coeffs;
		fft::transformed<E> f;
		fft::square_cached<E>(av.coeffs(), detail::whole_cache_or(a, ta_), coeffs, f);
		cached<E, true> w(exact<E>(std::move(coeffs)));
		w.cache() = std::move(f);
		return w;
	} else {
		vec<E, false> r(size_t(a.len()), T{});
		fft::square<E>(av.coeffs(), detail::whole_cache_or(a, ta_), std::span<T>(r));
		return r;
	}
}

// a*b + c*d, all exact; returns has_cache, adopting the summed pointwise
// product as the result's transform when the engine supports it. Reuses each
// operand's whole cache. Requires a*b and c*d to have equal length.
template <like A, like B, like C, like D>
	requires same_engine<A, B> && same_engine<A, C> && same_engine<A, D>
		&& A::exact_v && B::exact_v && C::exact_v && D::exact_v
cached<typename A::engine_t> multiply_add2(
		const A& a, const B& b, const C& c, const D& d) {
	using E = typename A::engine_t;
	using T = typename E::value_type;
	span<E, true> av = a.underlying(), bv = b.underlying();
	span<E, true> cv = c.underlying(), dv = d.underlying();
	fft::transformed<E> ta_, tb_, tc_, td_;
	std::vector<T> coeffs;
	fft::transformed<E> f;
	fft::multiply_add2_cached<E>(
		av.coeffs(), detail::whole_cache_or(a, ta_),
		bv.coeffs(), detail::whole_cache_or(b, tb_),
		cv.coeffs(), detail::whole_cache_or(c, tc_),
		dv.coeffs(), detail::whole_cache_or(d, td_),
		coeffs, f
	);
	cached<E> w(exact<E>(std::move(coeffs)));
	w.cache() = std::move(f);
	return w;
}

// coefficients [b.len()-1, a.len()) of a*b; requires a.len() >= b.len() > 0
template <like A, like B> requires same_engine<A, B>
std::vector<typename A::engine_t::value_type> middle_product(const A& a, const B& b) {
	using E = typename A::engine_t;
	span<E, A::exact_v> av = a.underlying();
	span<E, B::exact_v> bv = b.underlying();
	fft::transformed<E> ta_, tb_;
	return fft::middle_product<E>(
		av.coeffs(), detail::whole_cache_or(a, ta_),
		bv.coeffs(), detail::whole_cache_or(b, tb_)
	);
}

namespace detail {
template <bool ea, bool eb> int product_prec(int la, int lb) {
	if constexpr (ea && eb) return la > 0 && lb > 0 ? la + lb - 1 : 0;
	else return ea ? lb : eb ? la : std::min(la, lb);
}

// Normalize a product operand at the given precision to a borrowed series + the
// whole cache serving it: a prefix cache at scale nextPow2(prec), or the whole
// span with the operand's own cache, or a truncated span with the caller's
// throwaway cache.
// A whole cache of an over-length operand (len > prec, which pins the other,
// necessarily trunc, operand's span at exactly prec) is only worth using when
// the untruncated span doesn't grow the transform size: a 2x'd inverse
// transform costs more than the saved forward transform.
template <like S>
auto product_operand(const S& s, int prec, fft::transformed<typename S::engine_t>& tmp) {
	using E = typename S::engine_t;
	if constexpr (has_prefix_cache<S>) {
		return s.prefix(nextPow2(prec));
	} else {
		span<E, S::exact_v> v = s.underlying();
		int used = std::min(s.len(), prec);
		if constexpr (has_cache<S>) {
			if (s.len() <= prec || fft::detail::conv_size_for(s.len() + prec - 1).n
					== fft::detail::conv_size_for(2 * prec - 1).n) {
				return cached_span<E, S::exact_v>{v, s.cache()};
			}
		}
		return cached_span<E, S::exact_v>{v.first(used), tmp};
	}
}
/* namespace detail */ }

template <like A, like B> requires same_engine<A, B>
vec<typename A::engine_t, A::exact_v && B::exact_v> operator + (const A& a, const B& b) {
	using T = typename A::engine_t::value_type;
	span<typename A::engine_t, A::exact_v> av = a.underlying();
	span<typename A::engine_t, B::exact_v> bv = b.underlying();
	int n = (A::exact_v && B::exact_v) ? std::max(a.len(), b.len())
		: A::exact_v ? b.len() : B::exact_v ? a.len() : std::min(a.len(), b.len());
	vec<typename A::engine_t, A::exact_v && B::exact_v> r(size_t(n), T(0));
	for (int i = 0; i < n; i++) {
		r[i] = (i < av.len() ? av[i] : T(0)) + (i < bv.len() ? bv[i] : T(0));
	}
	return r;
}
template <like A, like B> requires same_engine<A, B>
vec<typename A::engine_t, A::exact_v && B::exact_v> operator - (const A& a, const B& b) {
	using T = typename A::engine_t::value_type;
	span<typename A::engine_t, A::exact_v> av = a.underlying();
	span<typename A::engine_t, B::exact_v> bv = b.underlying();
	int n = (A::exact_v && B::exact_v) ? std::max(a.len(), b.len())
		: A::exact_v ? b.len() : B::exact_v ? a.len() : std::min(a.len(), b.len());
	vec<typename A::engine_t, A::exact_v && B::exact_v> r(size_t(n), T(0));
	for (int i = 0; i < n; i++) {
		r[i] = (i < av.len() ? av[i] : T(0)) - (i < bv.len() ? bv[i] : T(0));
	}
	return r;
}

// The single multiplication operator: each operand is normalized to a borrowed
// series + whole cache (see detail::product_operand), then multiplied once.
// An exact x exact product returns a has_cache result, going through
// fft::multiply_cached so the pointwise product is adopted as the result's
// transform whenever the engine supports it.
template <like A, like B> requires same_engine<A, B>
auto operator * (const A& a, const B& b) {
	using E = typename A::engine_t;
	using T = typename E::value_type;
	constexpr bool ea = A::exact_v, eb = B::exact_v;
	int prec = detail::product_prec<ea, eb>(a.len(), b.len());
	if (prec == 0 || a.len() == 0 || b.len() == 0) {
		if constexpr (ea && eb) return cached<E, true>{};
		else return vec<E, false>(size_t(prec), T(0));
	}
	fft::transformed<E> ta_, tb_;
	auto va = detail::product_operand(a, prec, ta_);
	auto vb = detail::product_operand(b, prec, tb_);
	if constexpr (ea && eb) {
		std::vector<T> coeffs;
		fft::transformed<E> f;
		fft::multiply_cached<E>(
			va.underlying().coeffs(), va.cache(),
			vb.underlying().coeffs(), vb.cache(),
			coeffs, f
		);
		cached<E, true> w(exact<E>(std::move(coeffs)));
		w.cache() = std::move(f);
		return w;
	} else {
		vec<E, false> r(size_t(prec), T(0));
		fft::multiply<E>(
			va.underlying().coeffs(), va.cache(),
			vb.underlying().coeffs(), vb.cache(),
			std::span<T>(r)
		);
		return r;
	}
}

// [x^k] p(x)/q(x) (Bostan-Mori) for an exact rational function. Requires q[0] != 0 and
// p.len() < q.len(). Each level uses p(x) q(-x) (keeping the parity-of-k half) and
// q(x) q(-x) (even, giving the next q in x^2); q(-x)'s transform is negate_arg of q's,
// so a level costs 2 forward and 2 inverse transforms.
// TODO: downsample optimization
// TODO: support the kth_term_of_linear_recurrence(trunc, exact) form
template <fft::engine E>
typename E::value_type kth_term_of_rational_function(
	exact<E> p,
	exact<E> q,
	uint64_t k
) {
	using T = typename E::value_type;
	assert(q.len() > 0 && q[0] != T(0));
	assert(p.len() < q.len());
	int d = q.len();
	if (d == 1) return T(0);
	p.resize(d - 1);
	while (k > 0) {
		int n = nextPow2(2 * d - 1);
		auto tq = E::transform(std::span<const T>(q), n);
		auto tnq = E::negate_arg(tq, n);
		auto buf = fft::buffer_pool<T>::get(n);
		auto tp = E::transform(std::span<const T>(p), n);
		E::finish(E::mul(tp, tnq, n), buf.span());
		// deg(p * q(-x)) <= 2d-3 < n: wraparound-free
		for (int j = 0; j < d - 1; j++) p[j] = buf[2*j + int(k & 1)];
		E::finish(E::mul(tq, tnq, n), buf.span());
		// q(x) q(-x) is even with degree <= 2d-2
		for (int j = 0; j < d; j++) q[j] = buf[2*j];
		k >>= 1;
	}
	return p[0] * inv(q[0]);
}

// Wrapper around trunc which caches transform(s[:2^k]) for all k,
// matching the doubling shape of inverse/exp so they can populate the caches.
// TODO: make inverse/exp populate these
template <fft::engine E>
struct prefix_cached {
	using T = typename E::value_type;

	using engine_t = E;
	static constexpr bool exact_v = false;

	prefix_cached() = default;
	// moving coefficients in or out is free: implicit on rvalues, explicit copy otherwise
	prefix_cached(vec<E, false>&& s_) : s(std::move(s_)) {}
	explicit prefix_cached(const vec<E, false>& s_) : s(s_) {}
	operator vec<E, false>() && { return std::move(s); }

	int len() const { return s.len(); }
	const vec<E, false>& underlying() const { return s; }
	const T& operator[](int i) const { return s[size_t(i)]; }
	auto begin() const { return s.cbegin(); }
	auto end() const { return s.cend(); }

	// extend precision: appends coefficients, keeping all covering caches valid
	void append(std::span<const T> tail) {
		s.insert(s.end(), tail.begin(), tail.end());
	}

	// the length-min(n, len) prefix borrowed together with its cache
	cached_span<E, false> prefix(int n) const {
		return {
			span<E, false>(std::span<const T>(s).first(std::min(n, len()))),
			prefix_cache(n)
		};
	}
	// cache over the prefix of length min(n, len()); n a power of two
	fft::transformed<E>& prefix_cache(int n) const {
		assert(n > 0 && !(n & (n-1)));
		int k = __builtin_ctz(unsigned(n));
		if (k >= sz(caches)) caches.resize(size_t(k) + 1);
		auto& c = caches[k];
		int e = std::min(n, len());
		if (c.len != e) {
			c.t = E::transform(std::span<const T>(s).first(e), 2 * n);
			c.len = e;
		}
		return c.t;
	}

private:
	vec<E, false> s;
	// memoized transforms: logically const; len tracks how much of s each covers
	struct entry { fft::transformed<E> t; int len = 0; };
	mutable std::vector<entry> caches;
};

/* namespace series */ }

namespace poly {

// polynomial class
// As above, we represent polynomials by a series::exact containing the coefficients in reverse order.
// This representation should be internal-only:
// all accesses/constructors use the logical order though: P[k] = [x^k] P.
// To use the representation, use rev_series() / from_rev_series()
template <fft::engine E> struct vec {
	using T = typename E::value_type;
	using engine_t = E;
	series::exact<E> c;

	vec() = default;
	// zero polynomial with `len` coefficient slots
	explicit vec(int len) : c(size_t(len), T{}) {}
	// coefficient (x^0-first) order
	vec(std::initializer_list<T> coeffs) : c(std::rbegin(coeffs), std::rend(coeffs)) {}
	explicit vec(std::span<const T> coeffs) : c(coeffs.rbegin(), coeffs.rend()) {}

	const series::exact<E>& rev_series() const { return c; }
	static vec from_rev_series(series::exact<E> s) {
		vec r;
		r.c = std::move(s);
		return r;
	}

	// This should rarely be used
	series::vec<E> unrev_series(int n) const {
		series::vec<E> r(size_t(n), T{});
		std::copy(begin(), begin() + std::min(n, len()), r.begin());
		return r;
	}

	// logical (coefficient) order
	auto begin() { return c.rbegin(); }
	auto end() { return c.rend(); }
	auto begin() const { return c.rbegin(); }
	auto end() const { return c.rend(); }

	int len() const { return c.len(); }
	int degree() const { return len() - 1; }
	T& operator[](int i) { return c[len() - 1 - i]; }
	const T& operator[](int i) const { return c[len() - 1 - i]; }
	T leading() const { return c.front(); }
	// multiply by x^k: appends the new zero constant terms to the storage
	void shift(int k = 1) {
		if (len() > 0) c.insert(c.end(), size_t(k), T(0));
	}
	// grow (zero-filled leading coefficients) or shrink to n coefficients
	void resize(int n) {
		if (n >= len()) c.insert(c.begin(), size_t(n - len()), T(0));
		else c.erase(c.begin(), c.begin() + (len() - n));
	}

	T operator()(const T& x) const {
		T r{};
		for (const T& v : c) r = r * x + v;
		return r;
	}

	vec& operator+=(const vec& o) {
		if (o.len() > len()) resize(o.len());
		for (int i = 0; i < o.len(); i++) (*this)[i] += o[i];
		return *this;
	}
	friend vec operator+(vec a, const vec& b) { a += b; return a; }
	vec& operator-=(const vec& o) {
		if (o.len() > len()) resize(o.len());
		for (int i = 0; i < o.len(); i++) (*this)[i] -= o[i];
		return *this;
	}
	friend vec operator-(vec a, const vec& b) { a -= b; return a; }
	friend bool operator==(const vec& a, const vec& b) { return a.c == b.c; }

	vec& operator*=(const T& n) { for (T& v : c) v *= n; return *this; }
	friend vec operator*(vec a, const T& n) { a *= n; return a; }
	friend vec operator*(const T& n, vec a) { a *= n; return a; }

	vec& operator*=(const vec& o) { return *this = (*this) * o; }
};

// any polynomial representation exposing its reversed coefficient series
template <typename P>
concept like = requires(const P& p) {
	typename P::engine_t;
	{ p.len() } -> std::same_as<int>;
	p.rev_series();
	requires series::like<std::remove_cvref_t<decltype(p.rev_series())>>;
	requires std::remove_cvref_t<decltype(p.rev_series())>::exact_v;
};

// immutable polynomial carrying the whole-sequence transform of its rev_series
template <fft::engine E>
struct cached {
	using T = typename E::value_type;
	using engine_t = E;

	cached() = default;
	// moving coefficients in or out is free: implicit on rvalues, explicit copy otherwise
	cached(vec<E>&& p) : c(std::move(p.c)) {}
	explicit cached(const vec<E>& p) : c(p.c) {}
	operator vec<E>() && { return vec<E>::from_rev_series(std::move(c)); }

	const series::cached<E>& rev_series() const { return c; }
	static cached from_rev_series(series::cached<E> s) {
		cached r;
		r.c = std::move(s);
		return r;
	}

	int len() const { return c.len(); }
	int degree() const { return len() - 1; }
	const T& operator[](int i) const { return c[len() - 1 - i]; }
	T leading() const { return c[0]; }

	T operator()(const T& x) const {
		T r{};
		for (const T& v : c) r = r * x + v;
		return r;
	}

private:
	series::cached<E> c;
};

// rev(a*b) = rev(a)*rev(b); the series product reuses/adopts transforms
template <like A, like B> requires same_engine<A, B>
cached<typename A::engine_t> operator*(const A& a, const B& b) {
	return cached<typename A::engine_t>::from_rev_series(a.rev_series() * b.rev_series());
}
template <like A>
cached<typename A::engine_t> square(const A& a) {
	return cached<typename A::engine_t>::from_rev_series(square(a.rev_series()));
}
// rev(a*b + c*d) = rev(a)*rev(b) + rev(c)*rev(d)
template <like A, like B, like C, like D>
	requires same_engine<A, B> && same_engine<A, C> && same_engine<A, D>
cached<typename A::engine_t> multiply_add2(
		const A& a, const B& b, const C& c, const D& d) {
	return cached<typename A::engine_t>::from_rev_series(
			multiply_add2(a.rev_series(), b.rev_series(), c.rev_series(), d.rev_series()));
}
template <like A, like B> requires same_engine<A, B>
bool operator==(const A& a, const B& b) {
	return a.rev_series() == b.rev_series();
}

/* namespace poly */ }

// finite-support linear form
// These are one side of the pairing <poly::vec P, series::vec S> = [x^0] P(1/x) S(x).
// (Strictly speaking, this is actually <>_d where we take polynomials of degree < d.)
// The main point of this wrapper is that if we have <*, S> and want <P *, S>, that's a middle product by P.
//
// TODO: Should we split it apart into <*, S> and <P, *>?
//
// Some use cases of this pairing:
// <P, 1/(1-ax)> = P(a)
// if we represent P as a "polynomial" in the differential operator D (x^k = k! D^k):
// <P, e^{aD}> = P(a)
template <fft::engine E>
struct linear_form {
	using T = typename E::value_type;
	// coeffs of S in <*, S>; always whole-cached: the kernel transform is
	// what repeated middle products against the same form reuse
	series::cached<E> c;

	linear_form() = default;
	explicit linear_form(int len) : c(series::exact<E>(size_t(len), T{})) {}
	// We don't provide coefficient-list constructors, to avoid ordering confusion.

	const series::cached<E>& rev_series() const { return c; }
	static linear_form from_rev_series(series::cached<E> s) {
		linear_form r;
		r.c = std::move(s);
		return r;
	}
	static linear_form from_poly(const poly::vec<E>& p) { return from_rev_series(series::cached<E>(p.rev_series())); }

	int len() const { return c.len(); }

	// Restrict the form's domain: only valid against exact series of length n
	linear_form for_length(int n) const {
		series::exact<E> r(c.underlying());
		if (n >= len()) r.insert(r.begin(), size_t(n - len()), T(0));
		else r.erase(r.begin(), r.begin() + (len() - n));
		return from_rev_series(std::move(r));
	}

	// the functional p -> p(z) on polynomials of length up to len (weight z^i on [x^i])
	static linear_form polynomial_evaluation(T z, int len) {
		series::exact<E> k(size_t(len), T{});
		T p = T(1);
		for (int i = 0; i < len; i++) { k[i] = p; p *= z; }
		return from_rev_series(std::move(k));
	}

	template <poly::like P>
	T operator()(const P& p) const {
		assert(p.len() <= len());
		T r{};
		for (int i = 0; i < p.len(); i++) r += c[i] * p[i]; // weights multiply from the left
		return r;
	}

	// <*, S> -> <q x *, S>
	template <poly::like P>
	linear_form composed_with(const P& q) const {
		assert(q.len() > 0 && q.len() <= len());
		return from_rev_series(series::exact<E>(middle_product(c, q.rev_series())));
	}

	// <P, *> -> <P, s x *>
	template <series::like S> requires std::same_as<typename S::engine_t, E>
	linear_form composed_with(const S& s) const {
		if constexpr (!S::exact_v) assert(s.len() >= len());
		series::vec<E, S::exact_v> r = c * s;
		r.resize(size_t(len()));
		return from_rev_series(series::exact<E>(std::move(r)));
	}
};

// ==== multipoint evaluation / interpolation ====

// Subproduct tree over points a[0:N]
// BFS-order tree, each node holds prod (x - a[i]) as a cached poly::vec.
template <fft::engine E>
struct subproduct_tree {
	using T = typename E::value_type;
	int N;
	std::vector<poly::cached<E>> nodes;

	explicit subproduct_tree(std::span<const T> pts) : N(sz(pts)), nodes(size_t(2) * N) {
		assert(N > 0);
		for (int i = 0; i < N; i++) {
			nodes[N + i] = poly::vec<E>{-pts[i], T(1)};
		}
		for (int i = N - 1; i > 0; i--) {
			nodes[i] = nodes[2*i] * nodes[2*i+1];
		}
	}

	// number of points under node i
	int size(int i) const { return nodes[i].len() - 1; }
	// rev(prod (x - z_j)) over node i's leaves; length size(i) + 1
	const series::exact<E>& rev_prod(int i) const { return nodes[i].rev_series().underlying(); }

	// Computes, for each i, f(product_{j != i} (1 - a[j] x)). Requires f.len() == N.
	std::vector<T> pushdown(linear_form<E> f) const {
		assert(f.len() == N);
		std::vector<linear_form<E>> down(size_t(2) * N);
		down[1] = std::move(f);
		for (int i = 1; i < N; i++) {
			// the form's kernel transform serves both children's middle products
			down[2*i+0] = down[i].composed_with(nodes[2*i+1]);
			down[2*i+1] = down[i].composed_with(nodes[2*i+0]);
			down[i] = linear_form<E>{}; // done with the parent; free it early
		}
		std::vector<T> out(size_t(N), T{});
		for (int i = 0; i < N; i++) out[i] = down[N + i].rev_series()[0];
		return out;
	}

	// Compute sum_i leaf_vals[i] prod_{j!=i} (x - a[j]) (transpose of pushdown)
	poly::cached<E> combine_up(std::span<const T> leaf_vals) const {
		assert(sz(leaf_vals) == N);
		std::vector<poly::cached<E>> up(size_t(2) * N);
		for (int i = 0; i < N; i++) {
			up[N + i] = poly::vec<E>{leaf_vals[i]};
		}
		for (int i = N - 1; i > 0; i--) {
			up[i] = multiply_add2(up[2*i+0], nodes[2*i+1], up[2*i+1], nodes[2*i+0]);
			up[2*i+0] = poly::cached<E>{};
			up[2*i+1] = poly::cached<E>{};
		}
		return std::move(up[1]);
	}
};

template <fft::engine E>
std::vector<typename E::value_type> poly_evaluate(
	const poly::vec<E>& p,
	std::span<const typename E::value_type> pts
) {
	if (pts.empty()) return {};
	int N = sz(pts);
	subproduct_tree<E> tree{pts};
	series::vec<E> q = tree.rev_prod(1);
	q.resize(p.len()); // inverse precision must cover the form's window
	linear_form<E> f = linear_form<E>::from_poly(p).composed_with(inverse(q));
	return tree.pushdown(f.for_length(N));
}

template <fft::engine E>
poly::vec<E> poly_interpolate(
	std::span<const typename E::value_type> pts,
	std::span<const typename E::value_type> vals
) {
	using T = typename E::value_type;
	assert(sz(pts) == sz(vals));
	if (pts.empty()) return {};
	int N = sz(pts);
	using ps = series::vec<E>;
	subproduct_tree<E> tree{pts};
	ps root = tree.rev_prod(1);
	root.shrink(N);

	// We need to evaluate the derivative of the root at each point
	ps deriv_root = root;
	for (int i = 0; i < N; i++) {
		deriv_root[i] *= T(N - i);
	}
	std::vector<T> denoms = tree.pushdown(
		linear_form<E>::from_rev_series(series::exact<E>(inverse(root) * deriv_root))
	);

	std::vector<T> leaf_vals(size_t(N), T{});
	for (int i = 0; i < N; i++) leaf_vals[i] = vals[i] / denoms[i];
	return tree.combine_up(std::span<const T>(leaf_vals));
}

// ==== online multiplication ====

// Online (relaxed) multiplication: computes the first N terms of f*g given the terms one at a time.
template <fft::engine E> struct online_multiplier {
	using T = typename E::value_type;
	int N; int i;
	std::vector<T> f, g;
	std::vector<T> res;
	std::vector<fft::transformed<E>> f_blocks, g_blocks; // level k: block [2^k, 2^{k+1})

	online_multiplier(int N_) : N(N_), i(0), f(N, T{}), g(N, T{}), res(2*N+1, T{}) {}

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
				fft::transformed<E> cf, cg;
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

template <fft::engine E> struct online_squarer {
	using T = typename E::value_type;
	int N; int i;
	std::vector<T> f;
	std::vector<T> res;
	std::vector<fft::transformed<E>> f_blocks;

	online_squarer(int N_) : N(N_), i(0), f(N, T{}), res(2*N+1, T{}) {}

	T peek() {
		return res[i];
	}

	void push(T v_f) {
		assert(i < N);
		f[i] = v_f;
		if (i == 0) {
			res[0] += v_f * v_f;
		} else {
			if constexpr (E::commutative) res[i] += (v_f + v_f) * f[0];
			else res[i] += v_f * f[0] + f[0] * v_f;
			for (int p = 1, k = 0; (i & (p-1)) == (p-1); p <<= 1, k++) {
				int lo1 = p;
				int lo2 = i + 1 - p;
				int s = 2*p - 1;
				auto fb = std::span<const T>(f).subspan(p, p);
				auto fw = std::span<const T>(f).subspan(lo2, p);
				auto out = std::span<T>(res).subspan(lo1 + lo2, s);
				if (i == 2*p-1) {
					f_blocks.emplace_back();
					fft::square<E>(fb, f_blocks[k], out, fft::add_op{});
					break;
				}
				fft::transformed<E> cw;
				if constexpr (E::commutative) {
					fft::multiply<E>(fb, f_blocks[k], fw, cw, out, fft::add_twice_op{});
				} else {
					// f_hi * f_lo + f_lo * f_hi from the same two transforms
					fft::multiply_add2<E>(fb, f_blocks[k], fw, cw,
							fw, cw, fb, f_blocks[k], out, fft::add_op{});
				}
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
template <fft::engine E>
struct poly_ap_values : public std::vector<typename E::value_type> {
	using T = typename E::value_type;
	using std::vector<T>::vector;

	int len() const {
		return int(this->size());
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

		// Check for overlaps. We're checking in linear time to avoid unpacking T, but it should be plenty fast.
		// If the field is very very small and we wrap around several times, our runtime can be bad...
		// but then something has already gone wrong, why are you evaluating so many points???
		for (int i = -(len() - 1); i <= osz - 1; i++) {
			if (k+i == T(0)) {
				// everything from [i,i+len()-1) of the output is a match
				poly_ap_values res; res.reserve(osz);
				int lo = std::max(0, i);
				int hi = std::min(i+len(), osz);
				{
					auto pref = eval_range(k, lo);
					res.insert(res.end(), pref.begin(), pref.end());
				}
				res.insert(res.end(), this->begin() + (lo - i), this->begin() + (hi - i));
				{
					auto suff = eval_range(k + hi, osz - hi);
					res.insert(res.end(), suff.begin(), suff.end());
				}
				return res;
			}
		}

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
