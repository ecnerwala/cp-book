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
 *   conv_engines    Engines for packing/unpacking arbitrary rings for convolution.
 *                   Still expose (opaque) transform-domain objects for caching/fusion.
 *
 *   multiply layer  Wrappers for convolving bounded sequences: track length/truncation.
 *
 *   value types     power_series - R[[x]]
 *                   power_series_exact<E> - exact (finite-support) power series
 *                   power_series_trunc<E> - truncated prefix of an (infinite) power series
 *
 *                   polynomials - R[x]. Under x -> 1/x a polynomial becomes a Laurent polynomial in 1/x;
 *                                 shifting by x^{deg P} (reversal) lands it in R[[x]], and we store that exact series.
 *                   poly<E> - polynomial type, supporting natural indexing
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

// `conv_engine` contract
//   conv_engine represents a way of packing/unpacking sequences over an arbitrary ring into FFT-style transforms.
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
//      finish(product_t<A>, span<value_type>& out, Op) -> void
//
//   Input span can be length up to 2n.
//   Output spans can be length up to n; only the prefix that exists is filled.
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
//      extend_to       grow a transform by repeated doubling using the base engine's extend_to()
//      downsample      compute the half-sized transform/product of just the even (odd = false) or odd terms of the input
//      negate_arg      size n transform of A(-x)
template <typename E>
concept conv_engine = requires(
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
	E::finish(std::move(p), out);
	E::finish(std::move(p), out, add_op{});
	E::finish(E::add(std::move(p), std::move(p)), out);
	{ E::add(E::transform(in, n), ct) } -> std::same_as<typename E::template transformed_t<2 * E::unit_scale>>;
	{ E::add(std::move(p), std::move(p)) } -> std::same_as<typename E::template product_t<2 * E::unit_scale>>;
	requires std::same_as<std::remove_cvref_t<decltype(E::commutative)>, bool>;
	requires std::same_as<std::remove_cvref_t<decltype(E::unit_scale)>, int>;
};

// ==== scalar engines ====

template <typename num> struct fft_engine {
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
	static void extend_to(transformed& t, int n, std::span<const num> coeffs) {
		assert(sz(coeffs) <= 2 * t.size());
		while (t.size() < n) {
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
template <typename dbl = double> struct fft_real_engine {
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
template <typename mnum> struct fft_split_engine {
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
	static void extend_to(transformed& t, int n, std::span<const mnum> coeffs) {
		assert(sz(coeffs) <= 2 * t.size());
		auto buf = buffer_pool<cnum>::get(sz(coeffs));
		for (int i = 0; i < sz(coeffs); i++) buf[i] = pack(coeffs[i]);
		while (t.size() < n) {
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
	static void mul_impl(const vector<cnum>& a, const vector<cnum>& b, vector<cnum>& lo, vector<cnum>& hi, int n) {
		core::init(n);
		lo.resize(n); hi.resize(n);
		for (int i = 0; i < n; i++) {
			int ci = core::conj_index(i);
			cnum g0 = (b[i] + conj(b[ci])) * cnum(0.5);
			cnum t = (b[i] - conj(b[ci])) * cnum(0.5);
			cnum g1 = cnum(t.y, -t.x);
			lo[i] = a[i] * g0;
			hi[i] = a[i] * g1;
		}
	}
	template <int A, int B> static product_t<A * B> mul(const transformed_t<A>& a, const transformed_t<B>& b, int n) {
		assert(a.size() >= n && b.size() >= n);
		product_t<A * B> p;
		mul_impl(a.v, b.v, p.lo, p.hi, n);
		return p;
	}
	template <int A> static product_t<A * A> sq(const transformed_t<A>& a, int n) { return mul(a, a, n); }
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
		static_assert(K <= 2, "fft_split_engine: accumulated scale too large");
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
struct crt_engine {
	using value_type = mnum;
	static constexpr bool commutative = true;
	static constexpr int unit_scale = 1;
	using E1 = fft_engine<num1>;
	using E2 = fft_engine<num2>;
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
	static void extend_to(transformed& t, int n, std::span<const mnum> coeffs) {
		auto b1 = buffer_pool<num1>::get(sz(coeffs));
		auto b2 = buffer_pool<num2>::get(sz(coeffs));
		for (int i = 0; i < sz(coeffs); i++) { int64_t v = balanced(coeffs[i]); b1[i] = num1(v); b2[i] = num2(v); }
		E1::extend_to(t.t1, n, std::span<const num1>(b1.span()));
		E2::extend_to(t.t2, n, std::span<const num2>(b2.span()));
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
	template <int K1, int K2> static product_t<K1 + K2> add(product_t<K1>&& a, product_t<K2>&& b) {
		return product_t<K1 + K2>{E1::add(std::move(a.p1), b.p1), E2::add(std::move(a.p2), b.p2)};
	}
	template <int K = 1, typename Op = assign_op> static void finish(product_t<K>&& p, std::span<mnum> out, Op op = {}) {
		// The reconstruction needs |c| < whole/2; balanced inputs bound each addend's
		// true coefficients by n (MOD/2)^2, so the safe length is divided by the
		// accumulated scale. K <= 2 is very conservative (~2^35 even for MOD ~ 2^30).
		static_assert(K <= 2, "crt_engine: accumulated scale too large");
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

// componentwise_engine

// These are "componentwise engines" which model free modules/algebras over the underlying ring.
// Matrices are the canonical example: we can take entry-wise transforms, then multiply/add in transformed space.
//
// We start with a shared `componentwise_engine` base class which handles all linear ops, i.e. not mul.
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

template <conv_engine E, typename V, int L, std::array<int, size_t(L) + 1> Ofs = componentwise_iota<L>>
struct componentwise_engine {
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
	static void extend_to(transformed& t, int n, std::span<const V> coeffs) {
		auto buf = buffer_pool<S>::get(sz(coeffs));
		for (int c = 0; c < L; c++) {
			for (int i = 0; i < sz(coeffs); i++) buf[i] = coeffs[i].data()[c];
			E::extend_to(t.t[c], n, std::span<const S>(buf.span()));
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
		for (int c = 0; c < L; c++) {
			E::finish(std::move(p.t[Ofs[size_t(c)]]), buf.span());
			for (int j = Ofs[size_t(c)] + 1; j < Ofs[size_t(c) + 1]; j++)
				E::finish(std::move(p.t[j]), buf.span(), add_op{});
			for (int i = 0; i < sz(out); i++) op(out[i].data()[c], buf[i]);
		}
	}
};

// Convolve mat<N> (NxN matrices), with accumulation in product space
template <conv_engine E, int N>
struct matrix_engine : componentwise_engine<E, mat<typename E::value_type, N>, N * N> {
	using base = componentwise_engine<E, mat<typename E::value_type, N>, N * N>;
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
};

// Convolve trunc_series<num, N> (power series truncated at N), with accumulation in product space
template <conv_engine E, int N>
struct trunc_series_engine : componentwise_engine<E, trunc_series<typename E::value_type, N>, N> {
	using base = componentwise_engine<E, trunc_series<typename E::value_type, N>, N>;
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
};

// Stable variants of the wrapper engines: do not accumulate in product space.
// This costs an extra log factor.

template <int N> constexpr std::array<int, size_t(N) * N + 1> matrix_stable_ofs = [] {
	std::array<int, size_t(N) * N + 1> r{};
	for (int i = 0; i <= N * N; i++) r[size_t(i)] = i * N;
	return r;
}();

template <conv_engine E, int N>
struct matrix_engine_stable
		: componentwise_engine<E, mat<typename E::value_type, N>, N * N, matrix_stable_ofs<N>> {
	using base = componentwise_engine<E, mat<typename E::value_type, N>, N * N, matrix_stable_ofs<N>>;
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
};

template <int N> constexpr std::array<int, size_t(N) + 1> trunc_series_stable_ofs = [] {
	std::array<int, size_t(N) + 1> r{};
	for (int i = 0; i <= N; i++) r[size_t(i)] = i * (i + 1) / 2;
	return r;
}();

template <conv_engine E, int N>
struct trunc_series_engine_stable
		: componentwise_engine<E, trunc_series<typename E::value_type, N>, N, trunc_series_stable_ofs<N>> {
	using base = componentwise_engine<E, trunc_series<typename E::value_type, N>, N, trunc_series_stable_ofs<N>>;
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
};

static_assert(conv_engine<fft_engine<modnum<998244353>>>);
static_assert(conv_engine<fft_engine<mod_goldilocks>>);
static_assert(conv_engine<fft_real_engine<double>>);
static_assert(conv_engine<fft_split_engine<modnum<int(1e9)+7>>>);
static_assert(conv_engine<crt_engine<modnum<int(1e9)+7>>>);
static_assert(conv_engine<matrix_engine<fft_engine<modnum<998244353>>, 2>>);
static_assert(conv_engine<trunc_series_engine<fft_engine<modnum<998244353>>, 3>>);
// tracked inner engines work when the accumulated scale fits the budget (N <= 2)
static_assert(conv_engine<matrix_engine<fft_split_engine<modnum<int(1e9)+7>>, 2>>);
static_assert(conv_engine<trunc_series_engine<crt_engine<modnum<int(1e9)+7>>, 2>>);
// the stable variants keep tracked inner engines sound at any N
static_assert(conv_engine<matrix_engine_stable<fft_split_engine<modnum<int(1e9)+7>>, 3>>);
static_assert(conv_engine<trunc_series_engine_stable<crt_engine<modnum<int(1e9)+7>>, 3>>);

// A thin wrapper around E::transformed which manages the input length, as well as transformed::extend_to.
// TODO: Maybe this is useless? E::transformed presents essentially the same API.
// The best thing here is that we document a contract for E::transformed.
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

// ==== multiply layer ====
// These are free functions to convolve spans.
//
// The interfaces will typically take input spans, an output span, and an Op representing how to fold the result into the output.
// Output spans may alias one of the input spans.
// Output spans may be shorter than expected; the output will just be truncated.
//
// Some functions may also take fft_cache& objects associated with the input spans.
// These will be lazily filled and used if available.

// Circular convolution mod n (power of 2)
template <conv_engine E, typename Op = assign_op>
void multiply_circular(std::span<const typename E::value_type> a, std::span<const typename E::value_type> b,
		std::span<typename E::value_type> out, int n, Op op = {}) {
	assert(!(n & (n-1)));
	auto ta = E::transform(a, n);
	auto tb = E::transform(b, n);
	E::finish(E::mul(ta, tb, n), out, op);
}

template <conv_engine E, typename Op = assign_op>
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
}

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

// This helper also accepts an output-span fft_cache which will be populated if it is cheap to do so
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

// TODO: Let's decide whether to keep vector<> returning forms or not; this
// largely depends on whether we think these functions are a public interface or
// merely convenience for value type implementors.
template <conv_engine E> vector<typename E::value_type> middle_product(
		std::span<const typename E::value_type> a, std::span<const typename E::value_type> b) {
	using T = typename E::value_type;
	if (sz(a) == 0 || sz(b) == 0) return {};
	assert(sz(a) >= sz(b));
	vector<T> r(size_t(sz(a) - sz(b) + 1));
	middle_product<E>(a, b, std::span<T>(r));
	return r;
}

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
template <fft::conv_engine E> struct packed_bivariate {
	using T = typename E::value_type;
	int L, l;
	std::vector<T> c;

	// Q_0 = 1 - y g(x), deg g < n <= 2^L
	packed_bivariate(int L_, std::span<const T> g) : L(L_), l(0), c(size_t(4) << L) {
		c[0] = T(1);
		for (int i = 0; i < sz(g); i++) c[(2 << L) + i] = -g[i];
	}

	typename E::transformed advance() {
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

// Non-owning view of power series coefficients: the span pattern (contiguous
// window + series semantics), borrowed from an owning series-like type.
template <fft::conv_engine E, bool exact = false>
struct power_series_span {
	using T = typename E::value_type;
	using engine_t = E;
	static constexpr bool exact_v = exact;

	power_series_span() = default;
	explicit power_series_span(std::span<const T> s_) : s(s_) {}

	int len() const { return sz(s); }
	const T& operator[](int i) const { return s[size_t(i)]; }
	auto begin() const { return s.begin(); }
	auto end() const { return s.end(); }
	std::span<const T> coeffs() const { return s; }
	power_series_span underlying() const { return *this; }
	power_series_span first(int n) const { return power_series_span(s.first(size_t(n))); }

private:
	std::span<const T> s;
};

// `power_series` represents both exact (finite) power series (R[x]) and prefixes of infinite power series (R[[x]]), depending on the flag.
// `power_series_exact` and `power_series_trunc` are aliases.
//
// Operators here are typically permissive: they will accept combinations of unequal types and lengths.
template <fft::conv_engine E, bool exact = false>
struct power_series : public std::vector<typename E::value_type> {
	using T = typename E::value_type;
	using engine_t = E;
	static constexpr bool exact_v = exact;
	using std::vector<T>::vector;

	const power_series& underlying() const { return *this; }
	// a free const borrow of the coefficients: implicit
	operator power_series_span<E, exact>() const {
		return power_series_span<E, exact>(std::span<const T>(*this));
	}

	// exact -> trunc is implicit, trunc -> exact is explicit
	template <bool oe> requires (oe && !exact)
	power_series(const power_series<E, oe>& p) : std::vector<T>(p) {}
	template <bool oe> requires (oe && !exact)
	power_series(power_series<E, oe>&& p) : std::vector<T>(std::move(p)) {}
	template <bool oe> requires (!oe && exact)
	explicit power_series(const power_series<E, oe>& p) : std::vector<T>(p) {}
	template <bool oe> requires (!oe && exact)
	explicit power_series(power_series<E, oe>&& p) : std::vector<T>(std::move(p)) {}

	// adopt a plain coefficient vector
	explicit power_series(std::vector<T> v) : std::vector<T>(std::move(v)) {}

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
	power_series& operator += (const power_series<E, oe>& o) {
		if constexpr (exact) { if (o.len() > len()) this->resize(o.len()); }
		else if constexpr (!oe) { if (o.len() < len()) this->resize(o.len()); }
		for (int i = 0; i < std::min(len(), o.len()); i++) {
			(*this)[i] += o[i];
		}
		return *this;
	}
	template <bool oe> requires (oe || !exact)
	power_series& operator -= (const power_series<E, oe>& o) {
		if constexpr (exact) { if (o.len() > len()) this->resize(o.len()); }
		else if constexpr (!oe) { if (o.len() < len()) this->resize(o.len()); }
		for (int i = 0; i < std::min(len(), o.len()); i++) {
			(*this)[i] -= o[i];
		}
		return *this;
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

	// mixed-exactness binary operators (hidden friends: found by ADL on either operand)
	template <bool oe>
	friend power_series<E, exact && oe> operator + (const power_series& a, const power_series<E, oe>& b) {
		int n = (exact && oe) ? std::max(a.len(), b.len())
			: exact ? b.len() : oe ? a.len() : std::min(a.len(), b.len());
		power_series<E, exact && oe> r(size_t(n), T(0));
		for (int i = 0; i < n; i++) {
			r[i] = (i < a.len() ? a[i] : T(0)) + (i < b.len() ? b[i] : T(0));
		}
		return r;
	}
	template <bool oe>
	friend power_series<E, exact && oe> operator - (const power_series& a, const power_series<E, oe>& b) {
		int n = (exact && oe) ? std::max(a.len(), b.len())
			: exact ? b.len() : oe ? a.len() : std::min(a.len(), b.len());
		power_series<E, exact && oe> r(size_t(n), T(0));
		for (int i = 0; i < n; i++) {
			r[i] = (i < a.len() ? a[i] : T(0)) - (i < b.len() ? b[i] : T(0));
		}
		return r;
	}
	template <bool oe>
	friend power_series<E, exact && oe> operator * (const power_series& a, const power_series<E, oe>& b) {
		if constexpr (exact && oe) {
			if (a.len() == 0 || b.len() == 0) return {};
			power_series<E, true> r(size_t(a.len() + b.len() - 1), T(0));
			fft::multiply<E>(std::span<const T>(a), std::span<const T>(b), std::span<T>(r));
			return r;
		} else {
			int prec = exact ? b.len() : oe ? a.len() : std::min(a.len(), b.len());
			power_series<E, false> r(size_t(prec), T(0));
			if (prec == 0 || a.len() == 0 || b.len() == 0) return r;
			fft::multiply<E>(
				std::span<const T>(a).first(std::min(a.len(), prec)),
				std::span<const T>(b).first(std::min(b.len(), prec)),
				std::span<T>(r)
			);
			return r;
		}
	}

	power_series& operator *= (const power_series& o) {
		return *this = (*this) * o;
	}
	friend power_series square(const power_series& a) {
		if (sz(a) == 0) return {};
		power_series r(size_t(exact ? 2 * a.len() - 1 : a.len()));
		fft::square<E>(std::span<const T>(a), std::span<T>(r));
		return r;
	}

	// Newton inversion: 1/a mod x^a.len(). Generic over any engine; per doubling step
	// n -> m = 2n this is 5 transforms of size m, reusing b's transform for both circular
	// products; in each product the wraparound only contaminates coefficients [0, n)
	// which are already known.
	//
	// This is correct for non-commutative rings.
	friend power_series inverse(const power_series& a) requires (!exact) {
		power_series r(a.size());
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
	friend power_series deriv_shift_log(power_series a) requires (!exact) {
		auto r = deriv_shift(a);
		return r * inverse(a);
	}
	friend power_series poly_log(power_series a) requires (!exact) {
		assert(a[0] == 1);
		return integ_shift(deriv_shift_log(std::move(a)));
	}
	friend power_series poly_exp(power_series a) requires (!exact) {
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
	friend power_series poly_pow_monic(power_series a, T k) requires (!exact) {
		if (a.empty()) return a;
		assert(a.size() >= 1);
		assert(a[0] == 1);
		a = poly_log(a);
		a *= k;
		return poly_exp(a);
	}
	friend power_series poly_pow(power_series a, int64_t k) requires (!exact) {
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

	friend power_series to_newton_sums(const power_series& a, int deg) requires (!exact) {
		auto r = deriv_shift_log(a);
		r[0] = deg;
		for (int i = 1; i < int(r.size()); i++) r[i] = -r[i];
		return r;
	}
	friend power_series from_newton_sums(power_series S, int deg) requires (!exact) {
		assert(S[0] == int(deg));
		S[0] = 0;
		for (int i = 1; i < int(S.size()); i++) S[i] = -S[i];
		return poly_exp(integ_shift(std::move(S)));
	}

	// Calculates prod 1/(1-x^i)^{a[i]}
	friend power_series euler_transform(const power_series& a) requires (!exact) {
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
	friend power_series inverse_euler_transform(const power_series& a) requires (!exact) {
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
	friend power_series poly_compose(const power_series& f, const power_series& g) requires (!exact) {
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
		std::vector<typename E::transformed> tneg;
		tneg.reserve(L);
		for (int l = 1; l <= L; l++) tneg.push_back(Q.advance());
		power_series P;
		{
			P = f;
			std::reverse(P.begin(), P.end());
			power_series QL((1 << L) + 1);
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
		return power_series(P.begin(), P.begin() + n);
	}
};

template <fft::conv_engine E> using power_series_exact = power_series<E, true>;
template <fft::conv_engine E> using power_series_trunc = power_series<E, false>;

// Series-like concepts: the binary operators below are written once as constrained
// templates and dispatch on which memoized transforms an operand carries.
// A series-like type exposes its engine/exactness and its coefficients as a
// power_series_span borrow; cached wrappers additionally expose their transform
// caches (filling them is logically const).
template <typename S>
concept series_like = fft::conv_engine<typename S::engine_t> && requires(const S& s) {
	{ S::exact_v } -> std::convertible_to<bool>;
	{ s.len() } -> std::same_as<int>;
	{ s.underlying() } -> std::convertible_to<power_series_span<typename S::engine_t, S::exact_v>>;
};
// carries one extendable transform of the whole coefficient sequence
template <typename S>
concept whole_cached = series_like<S> && requires(const S& s) {
	{ s.cache() } -> std::same_as<fft::fft_cache<typename S::engine_t>&>;
};
template <typename A, typename B>
concept same_engine = std::same_as<typename A::engine_t, typename B::engine_t>;

// A borrowed series paired with the fft_cache serving its transforms: the
// normalized operand form fed to the cached fft:: entry points. Models whole_cached.
template <fft::conv_engine E, bool exact>
struct cached_power_series_span {
	using engine_t = E;
	static constexpr bool exact_v = exact;

	power_series_span<E, exact> s;
	std::reference_wrapper<fft::fft_cache<E>> f;

	int len() const { return s.len(); }
	power_series_span<E, exact> underlying() const { return s; }
	fft::fft_cache<E>& cache() const { return f; }
};

// carries transforms of power-of-two prefixes, usable at any precision:
// prefix(n) borrows the length-min(n+1, len) prefix with its cache
template <typename S>
concept prefix_cached = series_like<S> && requires(const S& s, int n) {
	{ s.prefix(n) } -> whole_cached;
};

// Wrapper around power_series which caches the transform of the whole series.
// Ops only exploit the cache at full (exact x exact) precision; a trunc series'
// whole-sequence transform is still useful for middle products and repeated
// full-precision use.
template <fft::conv_engine E, bool exact = true>
struct whole_cached_power_series {
	using T = typename E::value_type;
	using engine_t = E;
	static constexpr bool exact_v = exact;

	whole_cached_power_series() = default;
	// moving coefficients in or out is free: implicit on rvalues, explicit copy otherwise
	whole_cached_power_series(power_series<E, exact>&& s_) : s(std::move(s_)) {}
	explicit whole_cached_power_series(const power_series<E, exact>& s_) : s(s_) {}
	operator power_series<E, exact>() && { return std::move(s); }

	int len() const { return s.len(); }
	const power_series<E, exact>& underlying() const { return s; }
	const T& operator[](int i) const { return s[size_t(i)]; }
	auto begin() const { return s.cbegin(); }
	auto end() const { return s.cend(); }
	// the fft_cache over underlying(), fed to the cached fft:: entry points alongside it
	fft::fft_cache<E>& cache() const { return f; }

	template <series_like S>
	friend bool operator==(const whole_cached_power_series& a, const S& b) {
		power_series_span<E, S::exact_v> bs = b.underlying();
		return a.len() == bs.len() && std::equal(a.s.begin(), a.s.end(), bs.begin());
	}

private:
	power_series<E, exact> s;
	mutable fft::fft_cache<E> f; // memoized transform: filling it is logically const
};

template <whole_cached A>
power_series<typename A::engine_t, A::exact_v> square(const A& a) {
	using T = typename A::engine_t::value_type;
	if (a.len() == 0) return {};
	power_series_span<typename A::engine_t, A::exact_v> av = a.underlying();
	power_series<typename A::engine_t, A::exact_v> r(size_t(A::exact_v ? 2 * a.len() - 1 : a.len()), T{});
	if constexpr (A::exact_v) {
		fft::square<typename A::engine_t>(av.coeffs(), a.cache(), std::span<T>(r));
	} else {
		fft::square<typename A::engine_t>(av.coeffs(), std::span<T>(r));
	}
	return r;
}

// coefficients [b.len()-1, a.len()) of a*b; requires a.len() >= b.len() > 0
template <whole_cached A, whole_cached B> requires same_engine<A, B>
std::vector<typename A::engine_t::value_type> middle_product(const A& a, const B& b) {
	power_series_span<typename A::engine_t, A::exact_v> av = a.underlying();
	power_series_span<typename A::engine_t, B::exact_v> bv = b.underlying();
	return fft::middle_product<typename A::engine_t>(av.coeffs(), a.cache(), bv.coeffs(), b.cache());
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
template <series_like S>
auto product_operand(const S& s, int prec, fft::fft_cache<typename S::engine_t>& tmp) {
	using E = typename S::engine_t;
	if constexpr (prefix_cached<S>) {
		return s.prefix(nextPow2(prec));
	} else {
		power_series_span<E, S::exact_v> v = s.underlying();
		int used = std::min(s.len(), prec);
		if constexpr (whole_cached<S>) {
			if (s.len() <= prec || fft::detail::conv_size_for(s.len() + prec - 1).n
					== fft::detail::conv_size_for(2 * prec - 1).n) {
				return cached_power_series_span<E, S::exact_v>{v, s.cache()};
			}
		}
		return cached_power_series_span<E, S::exact_v>{v.first(used), tmp};
	}
}
/* namespace detail */ }

template <series_like A, series_like B> requires same_engine<A, B>
power_series<typename A::engine_t, A::exact_v && B::exact_v> operator + (const A& a, const B& b) {
	using T = typename A::engine_t::value_type;
	power_series_span<typename A::engine_t, A::exact_v> av = a.underlying();
	power_series_span<typename A::engine_t, B::exact_v> bv = b.underlying();
	int n = (A::exact_v && B::exact_v) ? std::max(a.len(), b.len())
		: A::exact_v ? b.len() : B::exact_v ? a.len() : std::min(a.len(), b.len());
	power_series<typename A::engine_t, A::exact_v && B::exact_v> r(size_t(n), T(0));
	for (int i = 0; i < n; i++) {
		r[i] = (i < av.len() ? av[i] : T(0)) + (i < bv.len() ? bv[i] : T(0));
	}
	return r;
}
template <series_like A, series_like B> requires same_engine<A, B>
power_series<typename A::engine_t, A::exact_v && B::exact_v> operator - (const A& a, const B& b) {
	using T = typename A::engine_t::value_type;
	power_series_span<typename A::engine_t, A::exact_v> av = a.underlying();
	power_series_span<typename A::engine_t, B::exact_v> bv = b.underlying();
	int n = (A::exact_v && B::exact_v) ? std::max(a.len(), b.len())
		: A::exact_v ? b.len() : B::exact_v ? a.len() : std::min(a.len(), b.len());
	power_series<typename A::engine_t, A::exact_v && B::exact_v> r(size_t(n), T(0));
	for (int i = 0; i < n; i++) {
		r[i] = (i < av.len() ? av[i] : T(0)) - (i < bv.len() ? bv[i] : T(0));
	}
	return r;
}

// The single multiplication operator: each operand is normalized to a borrowed
// series + whole cache (see detail::product_operand), then multiplied once.
// An exact x exact product returns a whole_cached result, going through
// fft::multiply_cached so the pointwise product is adopted as the result's
// transform whenever the engine supports it.
template <series_like A, series_like B> requires same_engine<A, B>
auto operator * (const A& a, const B& b) {
	using E = typename A::engine_t;
	using T = typename E::value_type;
	constexpr bool ea = A::exact_v, eb = B::exact_v;
	int prec = detail::product_prec<ea, eb>(a.len(), b.len());
	if (prec == 0 || a.len() == 0 || b.len() == 0) {
		if constexpr (ea && eb) return whole_cached_power_series<E, true>{};
		else return power_series<E, false>(size_t(prec), T(0));
	}
	fft::fft_cache<E> ta_, tb_;
	auto va = detail::product_operand(a, prec, ta_);
	auto vb = detail::product_operand(b, prec, tb_);
	if constexpr (ea && eb) {
		std::vector<T> coeffs;
		fft::fft_cache<E> f;
		fft::multiply_cached<E>(
			va.underlying().coeffs(), va.cache(),
			vb.underlying().coeffs(), vb.cache(),
			coeffs, f
		);
		whole_cached_power_series<E, true> w(power_series_exact<E>(std::move(coeffs)));
		w.cache() = std::move(f);
		return w;
	} else {
		power_series<E, false> r(size_t(prec), T(0));
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
// TODO: support the kth_term_of_linear_recurrence(power_series_trunc, power_series_exact) form
template <fft::conv_engine E>
typename E::value_type kth_term_of_rational_function(
	power_series_exact<E> p,
	power_series_exact<E> q,
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

// Wrapper around power_series_trunc which caches transform(s[:2^k]) for all k,
// matching the doubling shape of inverse/exp so they can populate the caches.
// TODO: make inverse/exp populate these
template <fft::conv_engine E, bool exact = false>
struct prefix_cached_power_series {
	using T = typename E::value_type;

	using engine_t = E;
	static constexpr bool exact_v = exact;

	prefix_cached_power_series() = default;
	// moving coefficients in or out is free: implicit on rvalues, explicit copy otherwise
	prefix_cached_power_series(power_series<E, exact>&& s_) : s(std::move(s_)) {}
	explicit prefix_cached_power_series(const power_series<E, exact>& s_) : s(s_) {}
	operator power_series<E, exact>() && { return std::move(s); }

	int len() const { return s.len(); }
	const power_series<E, exact>& underlying() const { return s; }
	const T& operator[](int i) const { return s[size_t(i)]; }
	auto begin() const { return s.cbegin(); }
	auto end() const { return s.cend(); }

	// extend precision: appends coefficients, keeping all covering caches valid
	void append(std::span<const T> tail) {
		s.insert(s.end(), tail.begin(), tail.end());
	}

	// the length-min(n, len) prefix borrowed together with its cache
	cached_power_series_span<E, exact> prefix(int n) const {
		return {
			power_series_span<E, exact>(std::span<const T>(s).first(std::min(n, len()))),
			prefix_cache(n)
		};
	}
	// cache over the prefix of length min(n, len()); n a power of two
	fft::fft_cache<E>& prefix_cache(int n) const {
		assert(n > 0 && !(n & (n-1)));
		int k = __builtin_ctz(unsigned(n));
		if (k >= sz(caches)) caches.resize(size_t(k) + 1);
		auto& c = caches[k];
		int e = std::min(n, len());
		if (c.len() != e) c = fft::fft_cache<E>(std::span<const T>(s).first(e), 2 * n);
		return c;
	}

private:
	power_series<E, exact> s;
	mutable std::vector<fft::fft_cache<E>> caches; // memoized transforms: logically const
};



// polynomial class
// As above, we represent polynomials by a power_series_exact containing the coefficients in reverse order.
// This representation should be internal-only:
// all accesses/constructors use the logical order though: P[k] = [x^k] P.
// To use the representation, use rev_series() / from_rev_series()
// TODO: This guy also needs caching integration
template <fft::conv_engine E> struct poly {
	using T = typename E::value_type;
	power_series_exact<E> c;

	poly() = default;
	// zero polynomial with `len` coefficient slots
	explicit poly(int len) : c(size_t(len), T{}) {}
	// coefficient (x^0-first) order
	poly(std::initializer_list<T> coeffs) : c(std::rbegin(coeffs), std::rend(coeffs)) {}
	explicit poly(std::span<const T> coeffs) : c(coeffs.rbegin(), coeffs.rend()) {}

	const power_series_exact<E>& rev_series() const { return c; }
	static poly from_rev_series(power_series_exact<E> s) {
		poly r;
		r.c = std::move(s);
		return r;
	}

	// This should rarely be used
	power_series<E> unrev_series(int n) const {
		power_series<E> r(size_t(n), T{});
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

	poly& operator+=(const poly& o) {
		if (o.len() > len()) resize(o.len());
		for (int i = 0; i < o.len(); i++) (*this)[i] += o[i];
		return *this;
	}
	friend poly operator+(poly a, const poly& b) { a += b; return a; }
	poly& operator-=(const poly& o) {
		if (o.len() > len()) resize(o.len());
		for (int i = 0; i < o.len(); i++) (*this)[i] -= o[i];
		return *this;
	}
	friend poly operator-(poly a, const poly& b) { a -= b; return a; }
	friend bool operator==(const poly& a, const poly& b) { return a.c == b.c; }

	poly& operator*=(const T& n) { for (T& v : c) v *= n; return *this; }
	friend poly operator*(poly a, const T& n) { a *= n; return a; }
	friend poly operator*(const T& n, poly a) { a *= n; return a; }

	friend poly operator*(const poly& a, const poly& b) {
		if (a.len() == 0 || b.len() == 0) return {};
		poly r(a.len() + b.len() - 1);
		fft::multiply<E>(std::span<const T>(a.c), std::span<const T>(b.c), std::span<T>(r.c));
		return r;
	}
	poly& operator*=(const poly& o) { return *this = (*this) * o; }
	friend poly square(const poly& a) {
		if (a.len() == 0) return {};
		poly r(2 * a.len() - 1);
		fft::square<E>(std::span<const T>(a.c), std::span<T>(r.c));
		return r;
	}
};

// finite-support linear form
// These are one side of the pairing <poly P, power_series S> = [x^0] P(1/x) S(x).
// (Strictly speaking, this is actually <>_d where we take polynomials of degree < d.)
// The main point of this wrapper is that if we have <*, S> and want <P *, S>, that's a middle product by P.
//
// TODO: Should we split it apart into <*, S> and <P, *>?
//
// Some use cases of this pairing:
// <P, 1/(1-ax)> = P(a)
// if we represent P as a "polynomial" in the differential operator D (x^k = k! D^k):
// <P, e^{aD}> = P(a)
template <fft::conv_engine E>
struct linear_form {
	using T = typename E::value_type;
	power_series_exact<E> c; // coeffs of S in <*, S>

	linear_form() = default;
	explicit linear_form(int len) : c(size_t(len), T{}) {}
	// We don't provide coefficient-list constructors, to avoid ordering confusion.

	const power_series_exact<E>& rev_series() const { return c; }
	static linear_form from_rev_series(power_series_exact<E> s) {
		linear_form r;
		r.c = std::move(s);
		return r;
	}
	static linear_form from_poly(const poly<E>& p) { return from_rev_series(p.rev_series()); }

	int len() const { return c.len(); }

	// Restrict the form's domain: only valid against exact series of length n
	linear_form for_length(int n) const {
		linear_form r = *this;
		if (n >= r.len()) r.c.insert(r.c.begin(), size_t(n - r.len()), T(0));
		else r.c.erase(r.c.begin(), r.c.begin() + (r.len() - n));
		return r;
	}

	// the functional p -> p(z) on polynomials of length up to len (weight z^i on [x^i])
	static linear_form polynomial_evaluation(T z, int len) {
		power_series_exact<E> k(size_t(len), T{});
		T p = T(1);
		for (int i = 0; i < len; i++) { k[i] = p; p *= z; }
		return from_rev_series(std::move(k));
	}

	T operator()(const poly<E>& p) const {
		assert(p.len() <= len());
		T r{};
		for (int i = 0; i < p.len(); i++) r += c[i] * p[i]; // weights multiply from the left
		return r;
	}

	// <*, S> -> <q x *, S>
	linear_form composed_with(const poly<E>& q) const {
		assert(q.len() > 0 && q.len() <= len());
		return from_rev_series(power_series_exact<E>(fft::middle_product<E>(
				std::span<const T>(c), std::span<const T>(q.rev_series()))));
	}

	// <P, *> -> <P, s x *>
	template <bool eb>
	linear_form composed_with(const power_series<E, eb>& s) const {
		if constexpr (!eb) assert(s.len() >= len());
		power_series<E, eb> r = c * s;
		r.resize(size_t(len()));
		return from_rev_series(power_series_exact<E>(std::move(r)));
	}
};

// ==== multipoint evaluation / interpolation ====

// Subproduct tree over points a[0:N]
// BFS-order tree, each node holds prod (1 - a[i] x) with caches.
template <fft::conv_engine E>
struct subproduct_tree {
	using T = typename E::value_type;
	int N;
	std::vector<whole_cached_power_series<E>> nodes;

	explicit subproduct_tree(std::span<const T> pts) : N(sz(pts)), nodes(size_t(2) * N) {
		assert(N > 0);
		for (int i = 0; i < N; i++) {
			nodes[N + i] = whole_cached_power_series<E>(power_series_exact<E>{T(1), -pts[i]});
		}
		for (int i = N - 1; i > 0; i--) {
			nodes[i] = nodes[2*i] * nodes[2*i+1];
		}
	}

	// number of points under node i
	int size(int i) const { return nodes[i].len() - 1; }
	// rev(prod (x - z_j)) over node i's leaves; length size(i) + 1
	const power_series_exact<E>& rev_prod(int i) const { return nodes[i].underlying(); }

	// Computes, for each i, f(product_{j != i} (1 - a[j] x)). Requires f.len() == N.
	std::vector<T> pushdown(linear_form<E> f) const {
		assert(f.len() == N);
		std::vector<linear_form<E>> down(size_t(2) * N);
		down[1] = std::move(f);
		for (int i = 1; i < N; i++) {
			// one transform of the kernel serves both children's middle products
			std::span<const T> k(down[i].rev_series());
			fft::fft_cache<E> ck;
			down[2*i+0] = linear_form<E>::from_rev_series(power_series_exact<E>(fft::middle_product<E>(
					k, ck, std::span<const T>(nodes[2*i+1].underlying()), nodes[2*i+1].cache())));
			down[2*i+1] = linear_form<E>::from_rev_series(power_series_exact<E>(fft::middle_product<E>(
					k, ck, std::span<const T>(nodes[2*i+0].underlying()), nodes[2*i+0].cache())));
		}
		std::vector<T> out(size_t(N), T{});
		for (int i = 0; i < N; i++) out[i] = down[N + i].rev_series()[0];
		return out;
	}

	// Compute sum_i leaf_vals[i] prod_{j!=i} (1-a[j] x) (transpose of pushdown)
	power_series_exact<E> combine_up(std::span<const T> leaf_vals) const {
		assert(sz(leaf_vals) == N);
		std::vector<power_series_exact<E>> up(size_t(2) * N);
		for (int i = 0; i < N; i++) {
			up[N + i] = power_series_exact<E>{leaf_vals[i]};
		}
		for (int i = N - 1; i > 0; i--) {
			power_series_exact<E> r(size_t(size(i)), T{});
			fft::fft_cache<E> cl, cr;
			fft::multiply_add2<E>(
					std::span<const T>(up[2*i+0]), cl,
					std::span<const T>(nodes[2*i+1].underlying()), nodes[2*i+1].cache(),
					std::span<const T>(up[2*i+1]), cr,
					std::span<const T>(nodes[2*i+0].underlying()), nodes[2*i+0].cache(),
					std::span<T>(r));
			up[i] = std::move(r);
		}
		return std::move(up[1]);
	}
};

template <fft::conv_engine E>
std::vector<typename E::value_type> poly_evaluate(
	const poly<E>& p,
	std::span<const typename E::value_type> pts
) {
	if (pts.empty()) return {};
	int N = sz(pts);
	subproduct_tree<E> tree{pts};
	power_series<E> q = tree.rev_prod(1);
	q.resize(p.len()); // inverse precision must cover the form's window
	linear_form<E> f = linear_form<E>::from_poly(p).composed_with(inverse(q));
	return tree.pushdown(f.for_length(N));
}

template <fft::conv_engine E>
poly<E> poly_interpolate(
	std::span<const typename E::value_type> pts,
	std::span<const typename E::value_type> vals
) {
	using T = typename E::value_type;
	assert(sz(pts) == sz(vals));
	if (pts.empty()) return {};
	int N = sz(pts);
	using ps = power_series<E>;
	subproduct_tree<E> tree{pts};
	ps root = tree.rev_prod(1);
	root.shrink(N);

	// We need to evaluate the derivative of the root at each point
	ps deriv_root = root;
	for (int i = 0; i < N; i++) {
		deriv_root[i] *= T(N - i);
	}
	std::vector<T> denoms = tree.pushdown(
		linear_form<E>::from_rev_series(power_series_exact<E>(inverse(root) * deriv_root))
	);

	std::vector<T> leaf_vals(size_t(N), T{});
	for (int i = 0; i < N; i++) leaf_vals[i] = vals[i] / denoms[i];
	return poly<E>::from_rev_series(tree.combine_up(std::span<const T>(leaf_vals)));
}

// ==== online multiplication ====

// Online (relaxed) multiplication: computes the first N terms of f*g given the terms one at a time.
template <fft::conv_engine E> struct online_multiplier {
	using T = typename E::value_type;
	int N; int i;
	std::vector<T> f, g;
	std::vector<T> res;
	std::vector<fft::fft_cache<E>> f_blocks, g_blocks; // level k: block [2^k, 2^{k+1})

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
				fft::fft_cache<E> cw;
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
template <fft::conv_engine E>
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
