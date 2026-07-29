#pragma once

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <initializer_list>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

#include "fft/series_core.hpp"

namespace ecnerwala::poly {

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
	series::trunc<E> unrev_series(int n) const {
		series::trunc<E> r(size_t(n), T{});
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

	const series::cached_exact<E>& rev_series() const { return c; }
	static cached from_rev_series(series::cached_exact<E> s) {
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
	series::cached_exact<E> c;
};

// rev(a*b) = rev(a)*rev(b); the series product reuses/adopts transforms
template <like A, like B> requires fft::same_engine<A, B>
cached<typename A::engine_t> operator*(const A& a, const B& b) {
	return cached<typename A::engine_t>::from_rev_series(a.rev_series() * b.rev_series());
}
template <like A>
cached<typename A::engine_t> square(const A& a) {
	return cached<typename A::engine_t>::from_rev_series(square(a.rev_series()));
}
// rev(a*b + c*d) = rev(a)*rev(b) + rev(c)*rev(d)
template <like A, like B, like C, like D>
	requires fft::same_engine<A, B> && fft::same_engine<A, C> && fft::same_engine<A, D>
cached<typename A::engine_t> multiply_add2(
		const A& a, const B& b, const C& c, const D& d) {
	return cached<typename A::engine_t>::from_rev_series(
			multiply_add2(a.rev_series(), b.rev_series(), c.rev_series(), d.rev_series()));
}
template <like A, like B> requires fft::same_engine<A, B>
bool operator==(const A& a, const B& b) {
	return a.rev_series() == b.rev_series();
}

// finite-support linear form
// These are one side of the pairing <vec P, series::vec S> = [x^0] P(1/x) S(x).
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
struct form {
	using T = typename E::value_type;
	// coeffs of S in <*, S>; always whole-cached: the kernel transform is
	// what repeated middle products against the same form reuse
	series::cached_exact<E> c;

	form() = default;
	explicit form(int len) : c(series::exact<E>(size_t(len), T{})) {}
	// We don't provide coefficient-list constructors, to avoid ordering confusion.

	const series::cached_exact<E>& rev_series() const { return c; }
	static form from_rev_series(series::cached_exact<E> s) {
		form r;
		r.c = std::move(s);
		return r;
	}
	static form from_poly(const vec<E>& p) { return from_rev_series(series::cached_exact<E>(p.rev_series())); }

	int len() const { return c.len(); }

	// Restrict the form's domain: only valid against exact series of length n
	form for_length(int n) const {
		series::exact<E> r(c.begin(), c.end());
		if (n >= len()) r.insert(r.begin(), size_t(n - len()), T(0));
		else r.erase(r.begin(), r.begin() + (len() - n));
		return from_rev_series(std::move(r));
	}

	// the functional p -> p(z) on polynomials of length up to len (weight z^i on [x^i])
	static form polynomial_evaluation(T z, int len) {
		series::exact<E> k(size_t(len), T{});
		T p = T(1);
		for (int i = 0; i < len; i++) { k[i] = p; p *= z; }
		return from_rev_series(std::move(k));
	}

	template <like P>
	T operator()(const P& p) const {
		assert(p.len() <= len());
		T r{};
		for (int i = 0; i < p.len(); i++) r += c[i] * p[i]; // weights multiply from the left
		return r;
	}

	// <*, S> -> <q x *, S>
	template <like P>
	form composed_with(const P& q) const {
		assert(q.len() > 0 && q.len() <= len());
		return from_rev_series(middle_product(c, q.rev_series()));
	}

	// <P, *> -> <P, s x *>
	template <series::like S> requires std::same_as<typename S::engine_t, E>
	form composed_with(const S& s) const {
		if constexpr (!S::exact_v) assert(s.len() >= len());
		series::vec<E, S::exact_v> r = c * s;
		r.resize(size_t(len()));
		return from_rev_series(series::exact<E>(std::move(r)));
	}
};

// ==== multipoint evaluation / interpolation ====

// Subproduct tree over points a[0:N]
// BFS-order tree, each node holds prod (x - a[i]) as a cached vec.
template <fft::engine E>
struct subproduct_tree {
	using T = typename E::value_type;
	int N;
	std::vector<cached<E>> nodes;

	explicit subproduct_tree(std::span<const T> pts) : N(sz(pts)), nodes(size_t(2) * N) {
		assert(N > 0);
		for (int i = 0; i < N; i++) {
			nodes[N + i] = vec<E>{-pts[i], T(1)};
		}
		for (int i = N - 1; i > 0; i--) {
			nodes[i] = nodes[2*i] * nodes[2*i+1];
		}
	}

	// number of points under node i
	int size(int i) const { return nodes[i].len() - 1; }
	// rev(prod (x - z_j)) over node i's leaves; length size(i) + 1
	series::exact<E> rev_prod(int i) const {
		const auto& c = nodes[i].rev_series();
		return series::exact<E>(c.begin(), c.end());
	}

	// Computes, for each i, f(product_{j != i} (1 - a[j] x)). Requires f.len() == N.
	std::vector<T> pushdown(form<E> f) const {
		assert(f.len() == N);
		std::vector<form<E>> down(size_t(2) * N);
		down[1] = std::move(f);
		for (int i = 1; i < N; i++) {
			// the form's kernel transform serves both children's middle products
			down[2*i+0] = down[i].composed_with(nodes[2*i+1]);
			down[2*i+1] = down[i].composed_with(nodes[2*i+0]);
			down[i] = form<E>{}; // done with the parent; free it early
		}
		std::vector<T> out(size_t(N), T{});
		for (int i = 0; i < N; i++) out[i] = down[N + i].rev_series()[0];
		return out;
	}

	// Compute sum_i leaf_vals[i] prod_{j!=i} (x - a[j]) (transpose of pushdown)
	cached<E> combine_up(std::span<const T> leaf_vals) const {
		assert(sz(leaf_vals) == N);
		std::vector<cached<E>> up(size_t(2) * N);
		for (int i = 0; i < N; i++) {
			up[N + i] = vec<E>{leaf_vals[i]};
		}
		for (int i = N - 1; i > 0; i--) {
			up[i] = multiply_add2(up[2*i+0], nodes[2*i+1], up[2*i+1], nodes[2*i+0]);
			up[2*i+0] = cached<E>{};
			up[2*i+1] = cached<E>{};
		}
		return std::move(up[1]);
	}
};

template <fft::engine E>
std::vector<typename E::value_type> multipoint(
	const vec<E>& p,
	std::span<const typename E::value_type> pts
) {
	if (pts.empty()) return {};
	int N = sz(pts);
	subproduct_tree<E> tree{pts};
	series::trunc<E> q = tree.rev_prod(1);
	q.resize(p.len()); // inverse precision must cover the form's window
	form<E> f = form<E>::from_poly(p).composed_with(ps_inv(q));
	return tree.pushdown(f.for_length(N));
}

template <fft::engine E>
vec<E> interpolate(
	std::span<const typename E::value_type> pts,
	std::span<const typename E::value_type> vals
) {
	using T = typename E::value_type;
	assert(sz(pts) == sz(vals));
	if (pts.empty()) return {};
	int N = sz(pts);
	using ps = series::trunc<E>;
	subproduct_tree<E> tree{pts};
	ps root = tree.rev_prod(1);
	root.shrink(N);

	// We need to evaluate the derivative of the root at each point
	ps deriv_root = root;
	for (int i = 0; i < N; i++) {
		deriv_root[i] *= T(N - i);
	}
	std::vector<T> denoms = tree.pushdown(
		form<E>::from_rev_series(series::exact<E>(ps_inv(root) * deriv_root))
	);

	std::vector<T> leaf_vals(size_t(N), T{});
	for (int i = 0; i < N; i++) leaf_vals[i] = vals[i] / denoms[i];
	return tree.combine_up(std::span<const T>(leaf_vals));
}

/* namespace ecnerwala::poly */ }
