#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

#include "fft/common.hpp"
#include "fft/engine.hpp"
#include "fft/multiply.hpp"
#include "fft/series_core.hpp"

// ==== analytic ops ====
// Free functions over series-like operands; each borrows the operand's span
// and writes a fresh result.
// TODO: reuse/populate the operands' whole/prefix transform caches

namespace ecnerwala::series {

template <like S>
vec<typename S::engine_t, S::exact_v> stretch(const S& a_, int n) {
	using E = typename S::engine_t;
	span<E, S::exact_v> a = a_.underlying();
	vec<E, S::exact_v> r(size_t(a.len()));
	for (int i = 0; i*n < a.len(); i++) {
		r[i*n] = a[i];
	}
	return r;
}
template <like S>
vec<typename S::engine_t, S::exact_v> deriv_shift(const S& a_) {
	using E = typename S::engine_t;
	span<E, S::exact_v> a = a_.underlying();
	vec<E, S::exact_v> r(a.begin(), a.end());
	for (int i = 0; i < r.len(); i++) {
		r[i] *= i;
	}
	return r;
}
template <like S>
vec<typename S::engine_t, S::exact_v> integ_shift(const S& a_) {
	using E = typename S::engine_t;
	using T = typename E::value_type;
	span<E, S::exact_v> a = a_.underlying();
	assert(a[0] == 0);
	vec<E, S::exact_v> r(a.begin(), a.end());
	T f = 1;
	for (int i = 1; i < r.len(); i++) {
		r[i] *= f;
		f *= i;
	}
	f = inv(f);
	for (int i = r.len() - 1; i > 0; i--) {
		r[i] *= f;
		f *= i;
	}
	return r;
}
template <like S>
vec<typename S::engine_t, S::exact_v> integ_shift_offset(const S& a_, int offset) {
	using E = typename S::engine_t;
	using T = typename E::value_type;
	span<E, S::exact_v> a = a_.underlying();
	vec<E, S::exact_v> r(a.begin(), a.end());
	T f = 1;
	for (int i = 0; i < r.len(); i++) {
		r[i] *= f;
		f *= i + offset;
	}
	assert(f != 0);
	f = inv(f);
	for (int i = r.len() - 1; i >= 0; i--) {
		r[i] *= f;
		f *= i + offset;
	}
	return r;
}
template <trunc_like S>
trunc<typename S::engine_t> deriv_shift_log(const S& a) {
	return deriv_shift(a) * ps_inv(a);
}
template <trunc_like S>
trunc<typename S::engine_t> ps_log(const S& a) {
	assert(a.underlying()[0] == 1);
	return integ_shift(deriv_shift_log(a));
}
template <trunc_like S>
trunc<typename S::engine_t> ps_exp(const S& a_) {
	// See https://mathexp.eu/bostan/publications/BoSc09a.pdf for details
	using E = typename S::engine_t;
	using T = typename E::value_type;
	span<E, false> a = a_.underlying();
	assert(a.len() >= 1);
	assert(a[0] == 0);
	trunc<E> r(1, T(1)); r.reserve(size_t(a.len()));
	trunc<E> invR(1, T(1)); invR.reserve(size_t(a.len()));
	while (r.len() < a.len()) {
		int o_sz = r.len();
		int n_sz = std::min(o_sz * 2, a.len());
		trunc<E> t = deriv_shift(trunc<E>(a.begin(), a.begin() + o_sz));
		fft::multiply_circular<E>(std::span<const T>(t), std::span<const T>(r).first(o_sz), std::span<T>(t), o_sz);
		t = deriv_shift(r) - t;
		t *= invR;
		t.resize(size_t(n_sz - o_sz));
		trunc<E> v(a.begin() + o_sz, a.begin() + n_sz);
		v -= integ_shift_offset(t, o_sz);
		v *= r;
		r.resize(size_t(n_sz));
		std::copy(v.begin(), v.end(), r.begin() + o_sz);
		if (r.len() < a.len()) {
			// double invR via a Newton step
			assert(r.len() == 2 * invR.len());
			int n = invR.len();
			int nn = r.len();
			trunc<E> tmp(size_t(4) * n);
			fft::square<E>(std::span<const T>(invR).first(n), std::span<T>(tmp));
			fft::multiply<E>(std::span<const T>(tmp).first(nn), std::span<const T>(r).first(nn), std::span<T>(tmp));
			invR.resize(size_t(nn));
			for (int i = n; i < nn; i++) invR[i] = -tmp[i];
		}
	}
	return r;
}
template <trunc_like S>
trunc<typename S::engine_t> ps_pow_monic(const S& a_, typename S::engine_t::value_type k) {
	using E = typename S::engine_t;
	span<E, false> a = a_.underlying();
	if (a.len() == 0) return {};
	assert(a[0] == 1);
	trunc<E> l = ps_log(a_);
	l *= k;
	return ps_exp(l);
}
template <trunc_like S>
trunc<typename S::engine_t> ps_pow(const S& a_, int64_t k) {
	using E = typename S::engine_t;
	using T = typename E::value_type;
	span<E, false> a = a_.underlying();
	assert(k >= 0);
	if (k == 0) {
		trunc<E> r(size_t(a.len()), T(0));
		if (r.len() > 0) r[0] = T(1);
		return r;
	}

	int st = 0;
	while (st < a.len() && a[st] == 0) st++;

	if (st > 0 && k > (a.len() - 1) / st) {
		return trunc<E>(size_t(a.len()), T(0));
	}

	trunc<E> r(a.begin() + st, a.end() - (st * (k-1)));
	T leading_coeff = r[0];
	r *= inv(leading_coeff);
	r = ps_pow_monic(r, T(k));
	r *= power(leading_coeff, k);
	r.insert(r.begin(), size_t(st * k), T(0));
	assert(r.len() == a.len());
	return r;
}

template <trunc_like S>
trunc<typename S::engine_t> to_newton_sums(const S& a, int deg) {
	auto r = deriv_shift_log(a);
	r[0] = deg;
	for (int i = 1; i < r.len(); i++) r[i] = -r[i];
	return r;
}
template <trunc_like S>
trunc<typename S::engine_t> from_newton_sums(const S& s_, int deg) {
	using E = typename S::engine_t;
	span<E, false> s = s_.underlying();
	assert(s[0] == deg);
	trunc<E> r(s.begin(), s.end());
	r[0] = 0;
	for (int i = 1; i < r.len(); i++) r[i] = -r[i];
	return ps_exp(integ_shift(std::move(r)));
}

// Calculates prod 1/(1-x^i)^{a[i]}
template <trunc_like S>
trunc<typename S::engine_t> euler_transform(const S& a) {
	using E = typename S::engine_t;
	trunc<E> r = deriv_shift(a);
	std::vector<bool> is_prime(size_t(r.len()), true);
	for (int p = 2; p < r.len(); p++) {
		if (!is_prime[p]) continue;
		for (int i = 1; i*p < r.len(); i++) {
			r[i*p] += r[i];
			is_prime[i*p] = false;
		}
	}
	return ps_exp(integ_shift(r));
}
template <trunc_like S>
trunc<typename S::engine_t> inverse_euler_transform(const S& a) {
	using E = typename S::engine_t;
	trunc<E> r = deriv_shift(ps_log(a));
	std::vector<bool> is_prime(size_t(r.len()), true);
	for (int p = 2; p < r.len(); p++) {
		if (!is_prime[p]) continue;
		for (int i = (r.len()-1)/p; i >= 1; i--) {
			r[i*p] -= r[i];
			is_prime[i*p] = false;
		}
	}
	return integ_shift(r);
}

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

// Calculates f(g(x)) mod x^n where deg(g) == n
template <trunc_like SF, trunc_like SG> requires fft::same_engine<SF, SG>
trunc<typename SF::engine_t> ps_compose(const SF& f_, const SG& g_) {
	using E = typename SF::engine_t;
	using T = typename E::value_type;
	span<E, false> f = f_.underlying();
	span<E, false> g = g_.underlying();
	if (g.len() == 0) return {};

	int m = f.len();
	int n = g.len();

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
	packed_bivariate<E> Q(L, g.coeffs());
	// tneg[l] is the transform of Q_l(-x, y), reused by the pushdown pass below
	std::vector<fft::transformed<E>> tneg;
	tneg.reserve(L);
	for (int l = 1; l <= L; l++) tneg.push_back(Q.advance());
	trunc<E> P;
	{
		P = trunc<E>(f.begin(), f.end());
		std::reverse(P.begin(), P.end());
		trunc<E> QL((1 << L) + 1);
		for (int i = 0; i <= (1 << L); i++) {
			QL[i] = Q.c[2 * i];
		}
		QL.resize(size_t(m), T(0));
		P *= ps_inv(QL);
		std::reverse(P.begin(), P.end());
		P.resize(size_t(1) << L, T(0));
		std::reverse(P.begin(), P.end());
		P.resize(size_t(B), T(0));
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
	return trunc<E>(P.begin(), P.begin() + n);
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

/* namespace ecnerwala::series */ }
