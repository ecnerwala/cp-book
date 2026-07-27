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
 *                   poly::form<E> - finite-support linear forms, via the pairing <P, S> = [x^0] P(1/x) S(x)
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


// Constrains two engine-parameterized value types to share the same engine.
template <typename A, typename B>
concept same_engine = std::same_as<typename A::engine_t, typename B::engine_t>;

// short spelling for E::transformed at use sites
template <engine E> using transformed = typename E::transformed;

/* namespace fft */ }

/* namespace ecnerwala */ }
