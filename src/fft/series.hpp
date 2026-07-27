#pragma once

#include "fft/multiply.hpp"

// ==== value types ====

namespace ecnerwala::series {

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
	friend vec ps_inv(const vec& a) requires (!exact) {
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
	// TODO: operator / can be done slightly faster than ps_inv:
	// we only need the n/2 terms of ps_inv(), and can do the last Newton step directly on the quotient

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
		return r * ps_inv(a);
	}
	friend vec ps_log(vec a) requires (!exact) {
		assert(a[0] == 1);
		return integ_shift(deriv_shift_log(std::move(a)));
	}
	friend vec ps_exp(vec a) requires (!exact) {
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
	friend vec ps_pow_monic(vec a, T k) requires (!exact) {
		if (a.empty()) return a;
		assert(a.size() >= 1);
		assert(a[0] == 1);
		a = ps_log(a);
		a *= k;
		return ps_exp(a);
	}
	friend vec ps_pow(vec a, int64_t k) requires (!exact) {
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
		r = ps_pow_monic(r, T(k));
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
		return ps_exp(integ_shift(std::move(S)));
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
		return ps_exp(integ_shift(r));
	}
	friend vec inverse_euler_transform(const vec& a) requires (!exact) {
		vec r = deriv_shift(ps_log(a));
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
	friend vec ps_compose(const vec& f, const vec& g) requires (!exact) {
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
			P *= ps_inv(QL);
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
	requires fft::same_engine<A, B> && fft::same_engine<A, C> && fft::same_engine<A, D>
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
template <like A, like B> requires fft::same_engine<A, B>
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

template <like A, like B> requires fft::same_engine<A, B>
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
template <like A, like B> requires fft::same_engine<A, B>
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
template <like A, like B> requires fft::same_engine<A, B>
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
// matching the doubling shape of ps_inv/exp so they can populate the caches.
// TODO: make ps_inv/exp populate these
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

/* namespace ecnerwala::series */ }
