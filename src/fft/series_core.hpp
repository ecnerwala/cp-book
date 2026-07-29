#pragma once

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

#include "fft/multiply.hpp"

// ==== value types ====

namespace ecnerwala::series {

// A series is either exact (a finite series, R[x] sitting inside R[[x]]: the
// length is just the support bound) or trunc (a known prefix of an infinite
// series: the length is the precision, and products truncate to it).

// Non-owning view of power series coefficients: the span pattern (contiguous
// window + series semantics), borrowed from an owning series-like type.
template <fft::engine E, bool exact_>
struct span {
	using T = typename E::value_type;
	using engine_t = E;
	static constexpr bool exact_v = exact_;

	span() = default;
	explicit span(std::span<const T> s_) : s(s_) {}

	int len() const { return sz(s); }
	const T& operator[](int i) const { return s[size_t(i)]; }
	auto begin() const { return s.begin(); }
	auto end() const { return s.end(); }
	std::span<const T> coeffs() const { return s; }
	// a series-like operand borrows straight into the engine primitives
	operator std::span<const T>() const { return s; }
	span underlying() const { return *this; }
	span first(int n) const { return span(s.first(size_t(n))); }

private:
	std::span<const T> s;
};

// `vec` represents both exact (finite) power series (R[x]) and prefixes of infinite power series (R[[x]]), depending on the flag.
// `exact` and `trunc` are aliases.
//
// Operators here are typically permissive: they will accept combinations of unequal types and lengths.
template <fft::engine E, bool exact_>
struct vec : public std::vector<typename E::value_type> {
	using T = typename E::value_type;
	using engine_t = E;
	static constexpr bool exact_v = exact_;
	using std::vector<T>::vector;

	const vec& underlying() const { return *this; }
	// a free const borrow of the coefficients: implicit
	operator span<E, exact_>() const {
		return span<E, exact_>(std::span<const T>(*this));
	}

	// exact -> trunc is implicit, trunc -> exact is explicit
	template <bool oe> requires (oe != exact_)
	explicit(oe < exact_) vec(const vec<E, oe>& p) : std::vector<T>(p) {}
	template <bool oe> requires (oe != exact_)
	explicit(oe < exact_) vec(vec<E, oe>&& p) : std::vector<T>(std::move(p)) {}

	// adopt a plain coefficient vector
	explicit vec(std::vector<T> v) : std::vector<T>(std::move(v)) {}

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
	template <bool oe> requires (exact_ <= oe)
	vec& operator += (const vec<E, oe>& o) {
		if constexpr (exact_) { if (o.len() > len()) this->resize(o.len()); }
		else if constexpr (!oe) { if (o.len() < len()) this->resize(o.len()); }
		for (int i = 0; i < std::min(len(), o.len()); i++) {
			(*this)[i] += o[i];
		}
		return *this;
	}
	template <bool oe> requires (exact_ <= oe)
	vec& operator -= (const vec<E, oe>& o) {
		if constexpr (exact_) { if (o.len() > len()) this->resize(o.len()); }
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

};

template <fft::engine E> using exact = vec<E, true>;
template <fft::engine E> using trunc = vec<E, false>;

// Series-like concepts: the binary operators below are written once as constrained
// templates and dispatch on which memoized transforms an operand carries.
// A series-like type exposes its engine/exactness and its coefficients as a
// span borrow; cached wrappers additionally expose their transform
// caches (filling them is logically const).
template <typename S>
concept like = fft::engine<typename S::engine_t> && requires(const S& s, int i) {
	{ S::exact_v } -> std::convertible_to<bool>;
	{ s.len() } -> std::same_as<int>;
	{ s[i] } -> std::convertible_to<const typename S::engine_t::value_type&>;
	// borrows straight into the engine primitives
	{ std::span<const typename S::engine_t::value_type>(s) };
	{ s.underlying() } -> std::convertible_to<span<typename S::engine_t, S::exact_v>>;
};
template <typename S>
concept exact_like = like<S> && S::exact_v;
template <typename S>
concept trunc_like = like<S> && !S::exact_v;

// carries one extendable transform of the whole coefficient sequence
template <typename S>
concept has_cache = like<S> && requires(const S& s) {
	{ s.cache() } -> std::same_as<fft::transformed<typename S::engine_t>&>;
};

// A borrowed series paired with the transform serving it: the
// normalized operand form fed to the cached fft:: entry points. Models has_cache.
template <fft::engine E, bool exact_>
struct cached_span {
	using engine_t = E;
	static constexpr bool exact_v = exact_;

	span<E, exact_> s;
	std::reference_wrapper<fft::transformed<E>> f;

	int len() const { return s.len(); }
	const typename E::value_type& operator[](int i) const { return s[i]; }
	operator std::span<const typename E::value_type>() const { return s.coeffs(); }
	span<E, exact_> underlying() const { return s; }
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
template <fft::engine E, bool exact_>
struct cached {
	using T = typename E::value_type;
	using engine_t = E;
	static constexpr bool exact_v = exact_;

	cached() = default;
	// moving coefficients in or out is free: implicit on rvalues, explicit copy otherwise
	cached(vec<E, exact_>&& s_) : s(std::move(s_)) {}
	explicit cached(const vec<E, exact_>& s_) : s(s_) {}
	operator vec<E, exact_>() && { return std::move(s); }

	int len() const { return s.len(); }
	const vec<E, exact_>& underlying() const { return s; }
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
	vec<E, exact_> s;
	mutable fft::transformed<E> f; // memoized transform: filling it is logically const
};

template <fft::engine E> using cached_exact = cached<E, true>;
template <fft::engine E> using cached_trunc = cached<E, false>;

namespace detail {
// the operand's whole cache if it carries one, else the caller's throwaway cache
template <like S>
fft::transformed<typename S::engine_t>& whole_cache_or(const S& s, fft::transformed<typename S::engine_t>& tmp) {
	if constexpr (has_cache<S>) return s.cache(); else return tmp;
}
// Normalize a whole-span operand to a cached_span: the coefficients borrowed
// together with the cache serving them (the operand's own, or tmp otherwise).
// The whole-span multiply/square/middle_product paths run entirely on this form.
template <like S>
cached_span<typename S::engine_t, S::exact_v> whole_operand(const S& s, fft::transformed<typename S::engine_t>& tmp) {
	return {span<typename S::engine_t, S::exact_v>(s.underlying()), whole_cache_or(s, tmp)};
}
/* namespace detail */ }

// Newton inversion: 1/a mod x^a.len(). Generic over any engine; per doubling step
// n -> m = 2n this is 5 transforms of size m, reusing b's transform for both circular
// products; in each product the wraparound only contaminates coefficients [0, n)
// which are already known.
//
// This is correct for non-commutative rings.
// TODO: reuse/populate the operand's whole/prefix transform caches
template <trunc_like S>
trunc<typename S::engine_t> ps_inv(const S& a_) {
	using E = typename S::engine_t;
	using T = typename E::value_type;
	span<E, false> a = a_.underlying();
	trunc<E> r(size_t(a.len()), T{});
	if (a.len() == 0) return r;
	int s = nextPow2(a.len());
	std::vector<T> b(size_t(s), T{});
	b[0] = inv(a[0]);
	for (int n = 1; n < a.len(); n *= 2) {
		int m = 2 * n;
		auto ta = E::transform(a.first(std::min(a.len(), m)), m);
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
// TODO: operator / can be done slightly faster than ps_inv:
// we only need the n/2 terms of ps_inv(), and can do the last Newton step directly on the quotient

// Both consume whole-sequence transforms by nature (the full span always
// participates), so only whole caches apply, never prefix caches.
template <like A>
auto square(const A& a) {
	using E = typename A::engine_t;
	using T = typename E::value_type;
	fft::transformed<E> ta_;
	auto av = detail::whole_operand(a, ta_);
	if constexpr (A::exact_v) {
		// like operator*, an exact square returns has_cache, adopting the
		// pointwise product as the result's transform when the engine supports it
		std::vector<T> coeffs;
		fft::transformed<E> f;
		fft::square_cached<E>(av, av.cache(), coeffs, f);
		cached<E, true> w(exact<E>(std::move(coeffs)));
		w.cache() = std::move(f);
		return w;
	} else {
		trunc<E> r(size_t(a.len()), T{});
		fft::square<E>(av, av.cache(), std::span<T>(r));
		return r;
	}
}

// a*b + c*d, all exact; returns has_cache, adopting the summed pointwise
// product as the result's transform when the engine supports it. Reuses each
// operand's whole cache. Requires a*b and c*d to have equal length.
template <like A, like B, like C, like D>
	requires fft::same_engine<A, B> && fft::same_engine<A, C> && fft::same_engine<A, D>
		&& A::exact_v && B::exact_v && C::exact_v && D::exact_v
cached<typename A::engine_t, true> multiply_add2(
		const A& a, const B& b, const C& c, const D& d) {
	using E = typename A::engine_t;
	using T = typename E::value_type;
	fft::transformed<E> ta_, tb_, tc_, td_;
	auto av = detail::whole_operand(a, ta_), bv = detail::whole_operand(b, tb_);
	auto cv = detail::whole_operand(c, tc_), dv = detail::whole_operand(d, td_);
	std::vector<T> coeffs;
	fft::transformed<E> f;
	fft::multiply_add2_cached<E>(
		av, av.cache(),
		bv, bv.cache(),
		cv, cv.cache(),
		dv, dv.cache(),
		coeffs, f
	);
	cached<E, true> w(exact<E>(std::move(coeffs)));
	w.cache() = std::move(f);
	return w;
}

// coefficients [b.len()-1, a.len()) of a*b; requires a.len() >= b.len() > 0.
// The kernel b participates whole, so it must be exact; the result mirrors a's kind.
template <like A, exact_like B> requires fft::same_engine<A, B>
vec<typename A::engine_t, A::exact_v> middle_product(const A& a, const B& b) {
	using E = typename A::engine_t;
	fft::transformed<E> ta_, tb_;
	auto av = detail::whole_operand(a, ta_);
	auto bv = detail::whole_operand(b, tb_);
	return vec<E, A::exact_v>(fft::middle_product<E>(
		av, av.cache(),
		bv, bv.cache()
	));
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

template <like A, like B> requires fft::same_engine<A, B>
vec<typename A::engine_t, A::exact_v && B::exact_v> operator + (const A& a, const B& b) {
	using T = typename A::engine_t::value_type;
	int n = (A::exact_v && B::exact_v) ? std::max(a.len(), b.len())
		: A::exact_v ? b.len() : B::exact_v ? a.len() : std::min(a.len(), b.len());
	vec<typename A::engine_t, A::exact_v && B::exact_v> r(size_t(n), T(0));
	for (int i = 0; i < n; i++) {
		r[i] = (i < a.len() ? a[i] : T(0)) + (i < b.len() ? b[i] : T(0));
	}
	return r;
}
template <like A, like B> requires fft::same_engine<A, B>
vec<typename A::engine_t, A::exact_v && B::exact_v> operator - (const A& a, const B& b) {
	using T = typename A::engine_t::value_type;
	int n = (A::exact_v && B::exact_v) ? std::max(a.len(), b.len())
		: A::exact_v ? b.len() : B::exact_v ? a.len() : std::min(a.len(), b.len());
	vec<typename A::engine_t, A::exact_v && B::exact_v> r(size_t(n), T(0));
	for (int i = 0; i < n; i++) {
		r[i] = (i < a.len() ? a[i] : T(0)) - (i < b.len() ? b[i] : T(0));
	}
	return r;
}

// The single multiplication operator: each operand is normalized to a borrowed
// series + whole cache (see detail::product_operand), then multiplied once.
// An exact x exact product returns a has_cache result, going through
// fft::multiply_cached so the pointwise product is adopted as the result's
// transform whenever the engine supports it.
template <like A, like B> requires fft::same_engine<A, B>
auto operator * (const A& a, const B& b) {
	using E = typename A::engine_t;
	using T = typename E::value_type;
	constexpr bool ea = A::exact_v, eb = B::exact_v;
	int prec = detail::product_prec<ea, eb>(a.len(), b.len());
	if (prec == 0 || a.len() == 0 || b.len() == 0) {
		if constexpr (ea && eb) return cached<E, true>{};
		else return trunc<E>(size_t(prec), T(0));
	}
	fft::transformed<E> ta_, tb_;
	auto va = detail::product_operand(a, prec, ta_);
	auto vb = detail::product_operand(b, prec, tb_);
	if constexpr (ea && eb) {
		std::vector<T> coeffs;
		fft::transformed<E> f;
		fft::multiply_cached<E>(
			va, va.cache(),
			vb, vb.cache(),
			coeffs, f
		);
		cached<E, true> w(exact<E>(std::move(coeffs)));
		w.cache() = std::move(f);
		return w;
	} else {
		trunc<E> r(size_t(prec), T(0));
		fft::multiply<E>(
			va, va.cache(),
			vb, vb.cache(),
			std::span<T>(r)
		);
		return r;
	}
}

// Wrapper around trunc which caches transform(s[:2^k]) for all k,
// matching the doubling shape of ps_inv/exp so they can populate the caches.
// TODO: make ps_inv/exp populate these
template <fft::engine E>
struct prefix_cached {
	using T = typename E::value_type;

	using engine_t = E;
	static constexpr bool exact_v = false;

	prefix_cached() = default;
	// moving coefficients in or out is free: implicit on rvalues, explicit copy otherwise
	prefix_cached(trunc<E>&& s_) : s(std::move(s_)) {}
	explicit prefix_cached(const trunc<E>& s_) : s(s_) {}
	operator trunc<E>() && { return std::move(s); }

	int len() const { return s.len(); }
	const trunc<E>& underlying() const { return s; }
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
	trunc<E> s;
	// memoized transforms: logically const; len tracks how much of s each covers
	struct entry { fft::transformed<E> t; int len = 0; };
	mutable std::vector<entry> caches;
};

/* namespace ecnerwala::series */ }
