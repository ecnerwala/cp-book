#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <span>
#include <vector>

#include "fft/engine.hpp"
#include "modnum.hpp"

namespace ecnerwala::fft {

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

/* namespace ecnerwala::fft */ }
