#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_get_random_seed.hpp>

#include "fft/series.hpp"
#include "fft/poly.hpp"
#include "fft/test_util.test.hpp"
#include "num/modnum.hpp"

namespace wala {
namespace fft {

using namespace std;

// series::like and poly::like are disjoint: neither family's types satisfy the other's concept
namespace {
using CE = engines::ntt<modnum<998244353>>;
static_assert(!series::like<poly::vec<CE>>);
static_assert(!series::like<poly::cached<CE>>);
static_assert(!series::like<poly::form<CE>>);
static_assert(!poly::like<series::exact<CE>>);
static_assert(!poly::like<series::trunc<CE>>);
static_assert(!poly::like<series::cached_exact<CE>>);
static_assert(!poly::like<series::cached_trunc<CE>>);
static_assert(!poly::like<series::operand<CE, true>>);
// every owning type and view is a series operand, cv/ref-qualified or not
static_assert(series::like<series::exact<CE>> && series::like<const series::trunc<CE>&>);
static_assert(series::like<series::cached_exact<CE>&> && series::like<const series::cached_trunc<CE>&>);
static_assert(series::like<series::operand<CE, true>> && series::like<series::operand<CE, false>>);
// exact widens to trunc, cached to operand, never the other way implicitly
static_assert(std::is_convertible_v<series::cached_exact<CE>&, series::operand<CE, false>>);
static_assert(!std::is_convertible_v<series::cached_trunc<CE>&, series::operand<CE, true>>);
static_assert(!std::is_convertible_v<series::cached_exact<CE>, series::exact<CE>>);
static_assert(!std::is_convertible_v<series::exact<CE>, series::cached_exact<CE>>);
static_assert(!std::is_convertible_v<std::span<const CE::value_type>, series::operand<CE, true>>);
static_assert(std::is_constructible_v<series::operand<CE, true>, std::span<const CE::value_type>>);
static_assert(!std::is_convertible_v<series::operand<CE, true>, series::exact<CE>>);
static_assert(std::is_constructible_v<series::exact<CE>, series::operand<CE, false>>);
}

// Archetype exposing exactly the series::like contract, nothing more.
// Instantiating the generic algorithms against it proves they only use
// contract expressions (concepts can't enforce that on function bodies).
namespace {
template <bool exact_>
struct archetype {
	using engine_t = CE;
	static constexpr bool exact_v = exact_;
	series::vec<CE, exact_> v;
	int len() const { return v.len(); }
	const typename CE::value_type& operator[](int i) const { return v[i]; }
	operator series::operand<CE, exact_>() const { return v; }
};
static_assert(series::like<archetype<true>> && series::like<archetype<false>>);

[[maybe_unused]] void archetype_instantiations(
	const archetype<true>& e, const archetype<false>& t, uint64_t k
) {
	series::stretch(e, 2); series::stretch(t, 2);
	series::deriv_shift(t); series::integ_shift(t); series::integ_shift_offset(t, 1);
	series::ogf_to_egf(e); series::egf_to_ogf(t);
	series::deriv_shift_log(t); series::ps_log(t); series::ps_exp(t);
	series::ps_pow_monic(t, {}); series::ps_pow(t, int64_t(k)); series::ps_inv(t);
	series::to_newton_sums(t, 1); series::from_newton_sums(t, 1);
	series::euler_transform(t); series::inverse_euler_transform(t);
	series::ps_compose(t, t);
	series::square(e); series::square(t);
	series::multiply_add2(e, e, e, e);
	series::middle_product(t, e); series::middle_product(e, e);
	series::operator*(e, e); series::operator*(e, t); series::operator*(t, t);
	series::operator+(e, t); series::operator-(t, t);
	series::kth_term_of_rational_function(e, e, k);
	series::kth_term_of_linear_recurrence(t, e, k);
	series::multiply(e, t, series::keep); series::square(t, series::keep);
	series::cached_exact<CE> ce;
	series::multiply(e, e, series::into(ce));
}
}

TEMPLATE_TEST_CASE("FFT Inverse", "[fft]", MOD_ENGINES) {
	using E = TestType;
	using num = typename E::value_type;
	mt19937 mt(Catch::getSeed());
	series::trunc<E> a(size_t(298));
	fill_rnd(a, mt);
	if (a[0] == 0) a[0] = 1;
	auto i = ps_inv(a);
	auto r = multiply_slow<num>(a, i);
	r.resize(a.size());
	vector<num> tgt(a.size());
	tgt[0] = 1;
	REQUIRE(r == tgt);
}

TEMPLATE_TEST_CASE("Bostan-Mori kth_term_of_rational_function", "[fft]", MOD_ENGINES) {
	using E = TestType;
	using num = typename E::value_type;
	mt19937 mt(Catch::getSeed());
	for (int d : {1, 2, 3, 8, 20}) {
		vector<num> p(d - 1), q(d);
		fill_rnd(p, mt);
		fill_rnd(q, mt);
		if (q[0] == 0) q[0] = 1;
		// reference: power series division to many terms
		int terms = 300;
		vector<num> ser(terms);
		num iq0 = inv(q[0]);
		for (int i = 0; i < terms; i++) {
			num v = i < int(p.size()) ? p[i] : num(0);
			for (int j = 1; j <= min<int>(i, d - 1); j++) v -= q[j] * ser[i-j];
			ser[i] = v * iq0;
		}
		series::exact<E> xp(p.begin(), p.end()), xq(q.begin(), q.end());
		for (uint64_t k : {uint64_t(0), uint64_t(1), uint64_t(7), uint64_t(100), uint64_t(299)}) {
			INFO("d = " << d << ", k = " << k);
			REQUIRE(kth_term_of_rational_function(xp, xq, k) == ser[k]);
		}
	}
}

TEMPLATE_TEST_CASE("Bostan-Mori kth_term_of_linear_recurrence", "[fft]", MOD_ENGINES) {
	using E = TestType;
	using num = typename E::value_type;
	mt19937 mt(Catch::getSeed());
	for (int d : {1, 2, 3, 8, 20}) {
		vector<num> p(d - 1), q(d);
		fill_rnd(p, mt);
		fill_rnd(q, mt);
		if (q[0] == 0) q[0] = 1;

		int terms = 300;
		vector<num> ser(terms);
		num iq0 = inv(q[0]);
		for (int i = 0; i < terms; i++) {
			if (i < int(p.size())) {
				ser[i] = p[i];
			} else {
				num v = 0;
				for (int j = 1; j <= min<int>(i, d - 1); j++) v -= q[j] * ser[i-j];
				ser[i] = v * iq0;
			}
		}
		series::trunc<E> xp(p.begin(), p.end()); series::exact<E> xq(q.begin(), q.end());
		for (uint64_t k : {uint64_t(0), uint64_t(1), uint64_t(7), uint64_t(100), uint64_t(299)}) {
			INFO("d = " << d << ", k = " << k);
			REQUIRE(kth_term_of_linear_recurrence(xp, xq, k) == ser[k]);
		}
	}
}

TEMPLATE_TEST_CASE("series result placement: plain, keep, into", "[fft]", MOD_ENGINES) {
	using E = TestType;
	using num = typename E::value_type;
	mt19937 mt(Catch::getSeed());
	series::exact<E> a(37), b(21), c(5);
	fill_rnd(a, mt);
	fill_rnd(b, mt);
	fill_rnd(c, mt);
	series::trunc<E> t(30);
	fill_rnd(t, mt);
	// ordinary products are plain vecs of the natural exactness
	auto p = a * b;
	static_assert(std::is_same_v<decltype(p), series::exact<E>>);
	check_eq(std::vector<num>(p), multiply_slow(a, b));
	static_assert(std::is_same_v<decltype(a * t), series::trunc<E>>);
	static_assert(std::is_same_v<decltype(series::square(t)), series::trunc<E>>);
	static_assert(std::is_same_v<decltype(series::middle_product(t, b)), series::trunc<E>>);
	// keep returns a cached whose transform (when the engine seeds it) is directly usable
	auto k = series::multiply(a, b, series::keep);
	static_assert(std::is_same_v<decltype(k), series::cached_exact<E>>);
	REQUIRE(k == p);
	if constexpr (std::same_as<typename E::product, fft::transformed<E>>) {
		REQUIRE(k.spectrum().size() > 0);
	}
	REQUIRE(k * c == p * c);
	REQUIRE(series::square(k, series::keep) == p * p);
	REQUIRE(series::multiply_add2(a, b, b, a, series::keep) == p + p);
	static_assert(std::is_same_v<decltype(series::multiply(a, t, series::keep)), series::cached_trunc<E>>);
	REQUIRE(series::multiply(a, t, series::keep) == a * t);
	// into: caller-owned outputs are resized (vec, cached) or filled up to their size (span)
	series::exact<E> out;
	series::multiply(a, b, series::into(out));
	REQUIRE(out == p);
	series::multiply(a, c, series::into(out));
	REQUIRE(out == a * c);
	series::cached_exact<E> kout;
	series::multiply(a, b, series::into(kout));
	REQUIRE(kout == p);
	REQUIRE(kout * c == p * c);
	std::vector<num> buf(10);
	series::multiply(a, b, series::into(std::span<num>(buf)));
	REQUIRE(equal(buf.begin(), buf.end(), p.begin()));
	series::trunc<E> tout;
	series::multiply(t, a, series::into(tout));
	REQUIRE(tout == t * a);
	// aliasing an operand is fine: the inputs are transformed before the output is written
	series::exact<E> a2 = a;
	series::multiply(a2, b, series::into(a2));
	REQUIRE(a2 == p);
	a2 = a;
	a2 *= b;
	REQUIRE(a2 == p);
	series::cached_exact<E> ca2(a);
	series::multiply(ca2, b, series::into(ca2));
	REQUIRE(ca2 == p);
}

TEMPLATE_TEST_CASE("series::cached operands: growth and const reuse", "[fft]", MOD_ENGINES) {
	using E = TestType;
	mt19937 mt(Catch::getSeed());
	series::exact<E> a(37), b(21), big(300);
	fill_rnd(a, mt);
	fill_rnd(b, mt);
	fill_rnd(big, mt);
	// a cached built from coefficients has no transform until a non-const use grows it
	series::cached_exact<E> ca(a);
	REQUIRE(ca.spectrum().size() == 0);
	REQUIRE(ca * b == a * b);
	int n1 = ca.spectrum().size();
	REQUIRE(n1 == fft::conv_size_for(a.len() + b.len() - 1).n);
	REQUIRE(ca * big == a * big);
	REQUIRE(ca.spectrum().size() == fft::conv_size_for(a.len() + big.len() - 1).n);
	// a const cached is reused when large enough, otherwise left alone
	const series::cached_exact<E> cb(b, 64);
	REQUIRE(cb.spectrum().size() == 64);
	REQUIRE(cb * a == a * b);
	REQUIRE(cb * big == b * big);
	REQUIRE(cb.spectrum().size() == 64);
	REQUIRE(series::middle_product(big, cb) == series::middle_product(big, b));
	// cached and plain operands are interchangeable at every mixed shape
	series::trunc<E> pa(40), pb(25), small(3);
	fill_rnd(pa, mt);
	fill_rnd(pb, mt);
	fill_rnd(small, mt);
	series::cached_trunc<E> qa(pa);
	const series::cached_trunc<E> qb(pb, 32);
	REQUIRE(qa * qb == pa * pb);
	REQUIRE(qa * pb == pa * pb);
	REQUIRE(pa * qb == pa * pb);
	REQUIRE(qa * small == pa * small);
	REQUIRE(small * qa == small * pa);
	REQUIRE(qa * a == pa * a);
	REQUIRE(ca * pb == a * pb);
	// leaving cached: coeffs() by const& or by move
	const series::exact<E>& ref = ca.coeffs();
	REQUIRE(ref == a);
	series::exact<E> moved = std::move(ca).coeffs();
	REQUIRE(moved == a);
	// views: first() keeps the transform only for the whole series
	series::operand<E, true> ob = cb;
	REQUIRE(ob.has_spectrum(64));
	REQUIRE(!ob.has_spectrum(128));
	REQUIRE(!ob.first(10).has_spectrum(1));
	REQUIRE(ob.first(b.len()).has_spectrum(64));
	REQUIRE(series::exact<E>(ob.first(10)) == series::exact<E>(b.begin(), b.begin() + 10));
}

TEST_CASE("series::vec mixed exactness operators", "[fft]") {
	using num = modnum<998244353>;
	using E = engines::ntt<num>;
	using xps = series::exact<E>;
	using ps = series::trunc<E>;
	mt19937 mt(Catch::getSeed());
	xps a(37), b(23);
	fill_rnd(a, mt);
	fill_rnd(b, mt);
	// exact * exact is the full product
	xps p = a * b;
	check_eq(std::vector<num>(p), multiply_slow(a, b));
	REQUIRE(square(a) == a * a);
	// exact +/- extend to the max length
	xps s = a + b, d = a - b;
	REQUIRE(s.len() == 37);
	for (int i = 0; i < 37; i++) REQUIRE(s[i] == a[i] + (i < 23 ? b[i] : num(0)));
	for (int i = 0; i < 37; i++) REQUIRE(d[i] == a[i] - (i < 23 ? b[i] : num(0)));
	// an exact operand doesn't lower a truncated result's precision
	ps t(b.begin(), b.end());
	ps m = a * t;
	REQUIRE(m.len() == 23);
	for (int i = 0; i < 23; i++) REQUIRE(m[i] == p[i]);
	REQUIRE(t * a == m);
	ps st = a + t, dt = t - a;
	REQUIRE(st.len() == 23);
	for (int i = 0; i < 23; i++) REQUIRE(st[i] == a[i] + t[i]);
	for (int i = 0; i < 23; i++) REQUIRE(dt[i] == t[i] - a[i]);
	// truncated * truncated is the min precision
	ps u(a.begin(), a.end());
	REQUIRE((u * t).len() == 23);
	REQUIRE(u * t == m);
	REQUIRE((u + t).len() == 23);
	// square of a truncated series keeps its precision
	REQUIRE(square(t).len() == 23);
	for (int i = 0; i < 23; i++) REQUIRE(square(t)[i] == (b * b)[i]);
	// exact -> truncated is implicit (forgetting exactness); the reverse is explicit
	ps forgot = p;
	REQUIRE(forgot.len() == p.len());
	REQUIRE(equal(forgot.begin(), forgot.end(), p.begin()));
	xps back(forgot);
	REQUIRE(back == p);
	static_assert(std::is_convertible_v<xps, ps>);
	static_assert(!std::is_convertible_v<ps, xps>);
	static_assert(std::is_constructible_v<xps, ps>);
}

TEST_CASE("series::vec log/exp/pow", "[fft]") {
	using num = modnum<998244353>;
	using ps = series::trunc<engines::ntt<num>>;
	mt19937 mt(Catch::getSeed());
	for (int len : {1, 2, 3, 17, 100}) {
		INFO("len = " << len);
		ps a(len);
		for (num& x : a) { x = num(mt()); }
		a[0] = 1;
		auto l = ps_log(a);
		auto e = ps_exp(l);
		REQUIRE(e == a);

		// ps_pow vs repeated multiplication
		ps p3 = a * a;
		p3 *= a;
		ps q = a;
		q[0] = 1;
		REQUIRE(ps_pow(q, 3) == p3);
	}
	{
		// pow with valuation
		ps a(20, num(0));
		for (int i = 3; i < 20; i++) a[i] = num(mt());
		if (a[3] == 0) a[3] = 1;
		ps p2 = a * a;
		REQUIRE(ps_pow(a, 2) == p2);
		// valuation * exponent overflowing the length gives 0
		REQUIRE(ps_pow(a, 100) == ps(20, num(0)));
	}
}

TEST_CASE("series::vec ps_inv", "[fft]") {
	using num = modnum<998244353>;
	using ps = series::trunc<engines::ntt<num>>;
	mt19937 mt(Catch::getSeed());
	ps a(100);
	for (num& x : a) { x = num(mt()); }
	if (a[0] == 0) a[0] = 1;
	ps i = ps_inv(a);
	ps prod = a * i;
	ps tgt(a.size(), num(0));
	tgt[0] = 1;
	REQUIRE(prod == tgt);
}

TEST_CASE("series::vec compose", "[fft]") {
	using num = modnum<998244353>;
	using ps = series::trunc<engines::ntt<num>>;
	mt19937 mt(Catch::getSeed());
	for (int n : {1, 2, 3, 8, 20, 33}) {
		INFO("n = " << n);
		int m = n + 2;
		ps f(m), g(n);
		for (num& x : f) { x = num(mt()); }
		for (num& x : g) { x = num(mt()); }
		g[0] = 0;
		// naive composition mod x^n
		ps expected(n, num(0));
		ps gp(n, num(0));
		gp[0] = 1;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) expected[j] += f[i] * gp[j];
			gp *= g;
		}
		REQUIRE(ps_compose(f, g) == expected);
	}
}

}} // namespace wala::fft
