#pragma once

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

#include "fft/multiply.hpp"

// ==== value types ====

namespace wala::series {

// A series is either exact (a finite series, R[x] sitting inside R[[x]]: the
// length is just the support bound) or trunc (a known prefix of an infinite
// series: the length is the precision, and products truncate to it).
//
// Owning types: vec<E, exact_> (a std::vector of coefficients) and
// cached<E, exact_> (coefficients paired with their transform, grown on demand).
// Every operation takes its operands through the one non-owning view,
// operand<E, exact_>, which any of the owning types converts to implicitly, so
// cached and plain series mix freely.
// Operations return plain vec results; a trailing result request changes where
// the result goes:
//   multiply(a, b)            fresh vec
//   multiply(a, b, keep)      fresh cached, seeded with the product's transform when the engine allows
//   multiply(a, b, into(x))   into x: a std::span<T> (prefix that fits), a vec& or cached& (resized)
// The into forms are the implementations; the others are thin wrappers.
// Inputs are transformed before the output is sized, so into(x) may alias an operand.
//
// Conversions: exact widens to trunc implicitly; trunc narrows to exact
// explicitly; owning -> operand is implicit; operand/std::span -> owning is an
// explicit copy (a raw span must state its exactness); a vec enters cached
// explicitly (cached(v)) and leaves via coeffs() (const& or moved out, which
// drops the transform).

using fft::keep_t;
using fft::keep;

template <fft::engine E, bool exact_> struct vec;
template <fft::engine E, bool exact_> struct cached;
template <fft::engine E, bool exact_> struct operand;

template <typename S> using engine_of = typename std::remove_cvref_t<S>::engine_t;
template <typename S> constexpr bool exact_of = std::remove_cvref_t<S>::exact_v;

// Anything an operation accepts as a series operand (any cv/ref qualification).
// Operations are templates over the operand type so that a growable cached&
// stays distinguishable from a read-only const cached&.
template <typename S>
concept like = requires {
	typename std::remove_cvref_t<S>::engine_t;
	{ std::remove_cvref_t<S>::exact_v } -> std::convertible_to<bool>;
} && std::convertible_to<S, operand<engine_of<S>, exact_of<S>>>;
template <typename S> concept exact_like = like<S> && exact_of<S>;
template <typename S> concept trunc_like = like<S> && !exact_of<S>;

// Owning series: a std::vector with series semantics.
template <fft::engine E, bool exact_>
struct vec : public std::vector<typename E::value_type> {
	using T = typename E::value_type;
	using engine_t = E;
	static constexpr bool exact_v = exact_;
	using std::vector<T>::vector;

	// exact -> trunc widens implicitly; trunc -> exact is explicit
	template <bool oe> requires (oe != exact_)
	explicit(oe < exact_) vec(const vec<E, oe>& p) : std::vector<T>(p) {}
	template <bool oe> requires (oe != exact_)
	explicit(oe < exact_) vec(vec<E, oe>&& p) : std::vector<T>(std::move(p)) {}

	explicit vec(std::vector<T> v) : std::vector<T>(std::move(v)) {}

	// copy out of a view, of either exactness (owning types go through their coeffs)
	template <bool oe>
	explicit vec(operand<E, oe> o) : std::vector<T>(o.begin(), o.end()) {}

	int len() const {
		return int(this->size());
	}
	int degree() const requires (exact_) {
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
	void shift_trunc(int n = 1) requires (!exact_) {
		assert(n >= 0 && n <= len());
		std::rotate(this->begin(), this->end()-n, this->end());
		std::fill(this->begin(), this->begin()+n, T(0));
	}
	// divide by x^n and 0-pad within the fixed precision window
	void unshift_trunc(int n = 1) requires (!exact_) {
		assert(n >= 0 && n <= len());
		std::fill(this->begin(), this->begin()+n, T(0));
		std::rotate(this->begin(), this->begin()+n, this->end());
	}

	// in-place forms require that the result's exactness/length must equal this operand's
	template <like S> requires (exact_ <= exact_of<S>) && fft::same_engine<vec, S>
	vec& operator += (const S& o) {
		if constexpr (exact_) { if (o.len() > len()) this->resize(o.len()); }
		else if constexpr (!exact_of<S>) { if (o.len() < len()) this->resize(o.len()); }
		for (int i = 0; i < std::min(len(), o.len()); i++) {
			(*this)[i] += o[i];
		}
		return *this;
	}
	template <like S> requires (exact_ <= exact_of<S>) && fft::same_engine<vec, S>
	vec& operator -= (const S& o) {
		if constexpr (exact_) { if (o.len() > len()) this->resize(o.len()); }
		else if constexpr (!exact_of<S>) { if (o.len() < len()) this->resize(o.len()); }
		for (int i = 0; i < std::min(len(), o.len()); i++) {
			(*this)[i] -= o[i];
		}
		return *this;
	}
	template <like S> requires (exact_ <= exact_of<S>) && fft::same_engine<vec, S>
	vec& operator *= (const S& o) {
		multiply(*this, o, into(*this));
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
};

template <fft::engine E> using exact = vec<E, true>;
template <fft::engine E> using trunc = vec<E, false>;

// A series paired with its transform.
// The coefficients are immutable while cached; the transform is built and
// grown (doubling, via E::extend_to) on non-const use, so const access reuses
// it only when it is already large enough.
template <fft::engine E, bool exact_>
struct cached {
	using T = typename E::value_type;
	using engine_t = E;
	static constexpr bool exact_v = exact_;

	cached() = default;
	explicit cached(vec<E, exact_> v) : s(std::move(v)) {}
	// with the transform built at size n (a power of two with len() <= 2n)
	cached(vec<E, exact_> v, int n) : s(std::move(v)) { spectrum(n); }

	const vec<E, exact_>& coeffs() const& { return s; }
	vec<E, exact_> coeffs() && { return std::move(s); }

	int len() const { return s.len(); }
	int degree() const requires (exact_) { return s.degree(); }
	const T& operator[](int i) const { return s[i]; }
	auto begin() const { return s.begin(); }
	auto end() const { return s.end(); }

	// the transform, grown to size >= n (a power of two with len() <= 2n)
	fft::transformed<E>& spectrum(int n) {
		E::extend_to(f, n, std::span<const T>(s));
		return f;
	}
	const fft::transformed<E>& spectrum() const { return f; }

	friend bool operator==(const cached& a, const cached& b) { return a.s == b.s; }
	friend bool operator==(const cached& a, const vec<E, exact_>& b) { return a.s == b; }

private:
	vec<E, exact_> s;
	fft::transformed<E> f;

	friend struct operand<E, true>;
	friend struct operand<E, false>;
	template <fft::engine, bool> friend struct cached_sink;
};

template <fft::engine E> using cached_exact = cached<E, true>;
template <fft::engine E> using cached_trunc = cached<E, false>;

// Non-owning view of a series: its coefficients, plus the transform of an
// underlying cached (growable when borrowed from a non-const one).
// first(n) keeps the transform only when it still views the whole series.
template <fft::engine E, bool exact_>
struct operand {
	using T = typename E::value_type;
	using engine_t = E;
	static constexpr bool exact_v = exact_;

	operand() = default;
	explicit operand(std::span<const T> c) : s(c) {}
	operand(const vec<E, exact_>& v) : s(v) {}
	operand(const vec<E, true>& v) requires (!exact_) : s(v) {}
	template <bool oe> requires (oe > exact_)
	operand(operand<E, oe> o) : s(o.s), f(o.f), g(o.g) {}

	template <bool oe> requires (oe >= exact_)
	operand(const cached<E, oe>& c) : s(c.s), f(&c.f) {}
	template <bool oe> requires (oe >= exact_)
	operand(cached<E, oe>& c) : s(c.s), f(&c.f), g(&c.f) {}

	// a contiguous range, so it converts to std::span<const T>
	std::span<const T> coeffs() const { return s; }
	int len() const { return int(s.size()); }
	size_t size() const { return s.size(); }
	const T* data() const { return s.data(); }
	const T& operator[](int i) const { return s[i]; }
	auto begin() const { return s.begin(); }
	auto end() const { return s.end(); }

	operand first(int n) const {
		assert(0 <= n && n <= len());
		operand r(s.first(n));
		if (n == len()) {
			r.f = f;
			r.g = g;
		}
		return r;
	}

	// This view using t (growable) as its transform if it has none of its own,
	// so one lent transform can serve several operations on a plain series.
	operand with_scratch(fft::transformed<E>& t) const {
		operand r = *this;
		if (!r.f) {
			r.f = &t;
			r.g = &t;
		}
		return r;
	}

	// Whether the view carries a transform usable at size n: growable, or already that large.
	bool has_spectrum(int n) const {
		return g || (f && int(f->size()) >= n);
	}
	// The transform for a size-n product over all len() coefficients:
	// the view's own (grown to n) when usable, else scratch, built there.
	const fft::transformed<E>& spectrum(int n, fft::transformed<E>& scratch) const {
		if (g) {
			E::extend_to(*g, n, s);
			return *g;
		}
		if (f && int(f->size()) >= n) return *f;
		E::extend_to(scratch, n, s);
		return scratch;
	}

private:
	std::span<const T> s;
	const fft::transformed<E>* f = nullptr; // the cached's transform, when viewing one
	fft::transformed<E>* g = nullptr;       // the same, when it may be grown (non-const cached)

	friend struct operand<E, false>;
};

// ==== result sinks ====

// A sink answers prepare(n) with the writable span for an n-coefficient result
// (the prefix that fits, for a std::span), and cached sinks expose spectrum()
// for the operation to seed or clear.
template <typename T>
struct span_sink {
	static constexpr bool keeps = false;
	static constexpr bool exact_v = false;
	std::span<T> out;
	std::span<T> prepare(int n) { return out.first(std::min(n, int(out.size()))); }
};
template <fft::engine E, bool exact_>
struct vec_sink {
	static constexpr bool keeps = false;
	static constexpr bool exact_v = exact_;
	vec<E, exact_>& out;
	std::span<typename E::value_type> prepare(int n) {
		out.resize(n);
		return out;
	}
};
template <fft::engine E, bool exact_>
struct cached_sink {
	static constexpr bool keeps = true;
	static constexpr bool exact_v = exact_;
	cached<E, exact_>& out;
	std::span<typename E::value_type> prepare(int n) {
		out.s.resize(n);
		return out.s;
	}
	fft::transformed<E>& spectrum() { return out.f; }
};

template <typename T> span_sink<T> into(std::span<T> out) { return {out}; }
template <fft::engine E, bool exact_> vec_sink<E, exact_> into(vec<E, exact_>& out) { return {out}; }
template <fft::engine E, bool exact_> cached_sink<E, exact_> into(cached<E, exact_>& out) { return {out}; }

template <typename S, typename E>
concept sink = requires(S s, int n) {
	{ s.prepare(n) } -> std::same_as<std::span<typename E::value_type>>;
	{ S::keeps } -> std::convertible_to<bool>;
	{ S::exact_v } -> std::convertible_to<bool>;
};
// an exact result may land in any sink; a trunc result only in a trunc or span sink
template <typename S, typename E, bool exact_>
concept sink_for = sink<S, E> && (exact_ || !S::exact_v);

namespace detail {

// zero-length results: n zeros, and no transform
template <fft::engine E, typename S>
void write_zero(S& out, int n) {
	using T = typename E::value_type;
	auto o = out.prepare(n);
	std::fill(o.begin(), o.end(), T{});
	if constexpr (S::keeps) out.spectrum() = fft::transformed<E>{};
}

template <fft::engine E, typename S>
void clear_spectrum(S& out) {
	if constexpr (S::keeps) out.spectrum() = fft::transformed<E>{};
}

template <bool ea, bool eb> int product_prec(int la, int lb) {
	if constexpr (ea && eb) return la > 0 && lb > 0 ? la + lb - 1 : 0;
	else return ea ? lb : eb ? la : std::min(la, lb);
}

// The coefficients an operand contributes to a product at precision prec:
// its prefix of length prec, or the whole operand when that keeps a usable
// transform and doesn't grow the transform size (an over-length operand pins
// the other, necessarily trunc, operand at exactly prec; a 2x'd inverse
// transform costs more than the saved forward transform).
template <fft::engine E, bool ex>
operand<E, ex> product_operand(operand<E, ex> o, int prec) {
	if (o.len() <= prec) return o;
	int n = fft::conv_size_for(o.len() + prec - 1).n;
	if (o.has_spectrum(n) && n == fft::conv_size_for(2 * prec - 1).n) return o;
	return o.first(prec);
}

/* namespace detail */ }

// ==== products ====

// a * b, at precision product_prec (exact x exact: full; else truncated to the
// trunc operand's precision).
template <like A, like B, sink_for<engine_of<A>, exact_of<A> && exact_of<B>> S>
	requires fft::same_engine<A, B>
void multiply(A&& a_, B&& b_, S out) {
	using E = engine_of<A>;
	using T = typename E::value_type;
	constexpr bool ea = exact_of<A>, eb = exact_of<B>;
	operand<E, ea> a = a_;
	operand<E, eb> b = b_;
	int prec = detail::product_prec<ea, eb>(a.len(), b.len());
	if (prec == 0 || a.len() == 0 || b.len() == 0) return detail::write_zero<E>(out, prec);
	a = detail::product_operand(a, prec);
	b = detail::product_operand(b, prec);
	int s = a.len() + b.len() - 1;
	auto [n, cut] = fft::conv_size_for(s);
	T c0 = a[0] * b[0];
	fft::transformed<E> sa, sb;
	auto p = E::mul(a.spectrum(n, sa), b.spectrum(n, sb), n);
	std::span<T> o = out.prepare(prec);
	if constexpr (S::keeps && ea && eb) {
		fft::finish_linear<E>(std::move(p), n, s, cut, c0, o, out.spectrum());
	} else {
		detail::clear_spectrum<E>(out);
		fft::finish_linear<E>(std::move(p), n, s, cut, c0, o);
	}
}

template <like A, like B> requires fft::same_engine<A, B>
vec<engine_of<A>, exact_of<A> && exact_of<B>> multiply(A&& a, B&& b) {
	vec<engine_of<A>, exact_of<A> && exact_of<B>> r;
	multiply(a, b, into(r));
	return r;
}
template <like A, like B> requires fft::same_engine<A, B>
cached<engine_of<A>, exact_of<A> && exact_of<B>> multiply(A&& a, B&& b, keep_t) {
	cached<engine_of<A>, exact_of<A> && exact_of<B>> r;
	multiply(a, b, into(r));
	return r;
}
template <like A, like B> requires fft::same_engine<A, B>
vec<engine_of<A>, exact_of<A> && exact_of<B>> operator * (A&& a, B&& b) {
	return multiply(a, b);
}

// a^2 (at a's precision when trunc)
template <like A, sink_for<engine_of<A>, exact_of<A>> S>
void square(A&& a_, S out) {
	using E = engine_of<A>;
	using T = typename E::value_type;
	constexpr bool ea = exact_of<A>;
	operand<E, ea> a = a_;
	int prec = detail::product_prec<ea, ea>(a.len(), a.len());
	if (prec == 0) return detail::write_zero<E>(out, prec);
	int s = 2 * a.len() - 1;
	auto [n, cut] = fft::conv_size_for(s);
	T c0 = a[0] * a[0];
	fft::transformed<E> sa;
	auto p = E::sq(a.spectrum(n, sa), n);
	std::span<T> o = out.prepare(prec);
	if constexpr (S::keeps && ea) {
		fft::finish_linear<E>(std::move(p), n, s, cut, c0, o, out.spectrum());
	} else {
		detail::clear_spectrum<E>(out);
		fft::finish_linear<E>(std::move(p), n, s, cut, c0, o);
	}
}
template <like A>
vec<engine_of<A>, exact_of<A>> square(A&& a) {
	vec<engine_of<A>, exact_of<A>> r;
	square(a, into(r));
	return r;
}
template <like A>
cached<engine_of<A>, exact_of<A>> square(A&& a, keep_t) {
	cached<engine_of<A>, exact_of<A>> r;
	square(a, into(r));
	return r;
}

// a*b + c*d, all exact, with a*b and c*d of equal nonzero length.
template <exact_like A, exact_like B, exact_like C, exact_like D, sink_for<engine_of<A>, true> S>
	requires fft::same_engine<A, B> && fft::same_engine<A, C> && fft::same_engine<A, D>
void multiply_add2(A&& a_, B&& b_, C&& c_, D&& d_, S out) {
	using E = engine_of<A>;
	using T = typename E::value_type;
	operand<E, true> a = a_, b = b_, c = c_, d = d_;
	assert(a.len() > 0 && b.len() > 0 && c.len() > 0 && d.len() > 0);
	int s = a.len() + b.len() - 1;
	assert(c.len() + d.len() - 1 == s);
	auto [n, cut] = fft::conv_size_for(s);
	T c0 = a[0] * b[0] + c[0] * d[0];
	fft::transformed<E> sa, sb, sc, sd;
	auto p = E::mul2(a.spectrum(n, sa), b.spectrum(n, sb), c.spectrum(n, sc), d.spectrum(n, sd), n);
	std::span<T> o = out.prepare(s);
	if constexpr (S::keeps) {
		fft::finish_linear<E>(std::move(p), n, s, cut, c0, o, out.spectrum());
	} else {
		fft::finish_linear<E>(std::move(p), n, s, cut, c0, o);
	}
}
template <exact_like A, exact_like B, exact_like C, exact_like D>
	requires fft::same_engine<A, B> && fft::same_engine<A, C> && fft::same_engine<A, D>
exact<engine_of<A>> multiply_add2(A&& a, B&& b, C&& c, D&& d) {
	exact<engine_of<A>> r;
	multiply_add2(a, b, c, d, into(r));
	return r;
}
template <exact_like A, exact_like B, exact_like C, exact_like D>
	requires fft::same_engine<A, B> && fft::same_engine<A, C> && fft::same_engine<A, D>
cached_exact<engine_of<A>> multiply_add2(A&& a, B&& b, C&& c, D&& d, keep_t) {
	cached_exact<engine_of<A>> r;
	multiply_add2(a, b, c, d, into(r));
	return r;
}

// coefficients [b.len()-1, a.len()) of a*b; requires a.len() >= b.len() > 0.
// The kernel b participates whole, so it must be exact; the result mirrors a's kind.
// The result's transform is never seeded (the middle slice has no pointwise product).
template <like A, exact_like B, sink_for<engine_of<A>, exact_of<A>> S>
	requires fft::same_engine<A, B>
void middle_product(A&& a_, B&& b_, S out) {
	using E = engine_of<A>;
	using T = typename E::value_type;
	operand<E, exact_of<A>> a = a_;
	operand<E, true> b = b_;
	assert(a.len() >= b.len() && b.len() > 0);
	int m = a.len() - b.len() + 1;
	detail::clear_spectrum<E>(out);
	if (a.len() == b.len()) {
		T r{};
		fft::middle_product_dot<T>(a, b, std::span<T>(&r, 1));
		std::span<T> o = out.prepare(1);
		if (o.size() > 0) o[0] = r;
		return;
	}
	auto [n, cut] = fft::conv_size_for(a.len());
	T c0 = a[0] * b[0], ctop = a[a.len() - 1] * b[b.len() - 1];
	fft::transformed<E> sa, sb;
	auto p = E::mul(a.spectrum(n, sa), b.spectrum(n, sb), n);
	std::span<T> o = out.prepare(m);
	fft::finish_middle<E>(std::move(p), n, cut, a.len(), b.len(), c0, ctop, o);
}
template <like A, exact_like B> requires fft::same_engine<A, B>
vec<engine_of<A>, exact_of<A>> middle_product(A&& a, B&& b) {
	vec<engine_of<A>, exact_of<A>> r;
	middle_product(a, b, into(r));
	return r;
}
template <like A, exact_like B> requires fft::same_engine<A, B>
cached<engine_of<A>, exact_of<A>> middle_product(A&& a, B&& b, keep_t) {
	cached<engine_of<A>, exact_of<A>> r;
	middle_product(a, b, into(r));
	return r;
}

// ==== elementwise ====

template <like A, like B> requires fft::same_engine<A, B>
vec<engine_of<A>, exact_of<A> && exact_of<B>> operator + (const A& a, const B& b) {
	using T = typename engine_of<A>::value_type;
	constexpr bool ea = exact_of<A>, eb = exact_of<B>;
	int n = (ea && eb) ? std::max(a.len(), b.len()) : ea ? b.len() : eb ? a.len() : std::min(a.len(), b.len());
	vec<engine_of<A>, ea && eb> r(size_t(n), T(0));
	for (int i = 0; i < n; i++) {
		r[i] = (i < a.len() ? a[i] : T(0)) + (i < b.len() ? b[i] : T(0));
	}
	return r;
}
template <like A, like B> requires fft::same_engine<A, B>
vec<engine_of<A>, exact_of<A> && exact_of<B>> operator - (const A& a, const B& b) {
	using T = typename engine_of<A>::value_type;
	constexpr bool ea = exact_of<A>, eb = exact_of<B>;
	int n = (ea && eb) ? std::max(a.len(), b.len()) : ea ? b.len() : eb ? a.len() : std::min(a.len(), b.len());
	vec<engine_of<A>, ea && eb> r(size_t(n), T(0));
	for (int i = 0; i < n; i++) {
		r[i] = (i < a.len() ? a[i] : T(0)) - (i < b.len() ? b[i] : T(0));
	}
	return r;
}

// ==== inversion ====

// Newton inversion: 1/a mod x^a.len(). Generic over any engine; per doubling step
// n -> m = 2n this is 5 transforms of size m, reusing b's transform for both circular
// products; in each product the wraparound only contaminates coefficients [0, n)
// which are already known.
//
// This is correct for non-commutative rings.
template <trunc_like S>
trunc<engine_of<S>> ps_inv(const S& a_) {
	using E = engine_of<S>;
	using T = typename E::value_type;
	operand<E, false> a = a_;
	int N = a.len();
	trunc<E> r(size_t(N), T{});
	if (N == 0) return r;
	int s = nextPow2(N);
	std::vector<T> b(size_t(s), T{});
	b[0] = inv(a[0]);
	for (int n = 1; n < N; n *= 2) {
		int m = 2 * n;
		auto ta = E::transform(a.first(std::min(N, m)).coeffs(), m);
		auto tb = E::transform(std::span<const T>(b).first(n), m);
		// e = a*b mod x^m; only e[n..m) is needed (and is wraparound-free).
		auto e = fft::buffer_pool<T>::get(m);
		E::finish(E::mul(ta, tb, m), e.span());
		for (int i = 0; i < n; i++) e[i] = T{};
		auto te = E::transform(std::span<const T>(e.span()), m);
		auto c = fft::buffer_pool<T>::get(m);
		// b' = 2b - b*(a*b): keep b on the left of e = a*b
		E::finish(E::mul(tb, te, m), c.span());
		for (int i = n; i < std::min(m, N); i++) b[i] = -c[i];
	}
	std::copy(b.begin(), b.begin() + N, r.begin());
	return r;
}
// TODO: operator / can be done slightly faster than ps_inv:
// we only need the n/2 terms of ps_inv(), and can do the last Newton step directly on the quotient

/* namespace wala::series */ }
