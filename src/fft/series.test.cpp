#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_get_random_seed.hpp>

#include "fft/series.hpp"
#include "fft/poly.hpp"
#include "fft/test_util.test.hpp"
#include "num/modnum.hpp"

namespace ecnerwala {
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
static_assert(!poly::like<series::prefix_cached<CE>>);
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
	operator std::span<const typename CE::value_type>() const { return std::span<const typename CE::value_type>(v); }
	operator series::span<CE, exact_>() const { return v; }
	series::span<CE, exact_> first(int n) const { return v.first(n); }
};
static_assert(series::like<archetype<true>> && series::like<archetype<false>>);
static_assert(!series::has_cache<archetype<false>> && !series::has_cache_opt<archetype<false>>);
static_assert(!series::has_prefix_cache<archetype<false>>);

[[maybe_unused]] void archetype_instantiations(
	const archetype<true>& e, const archetype<false>& t, uint64_t k
) {
	series::stretch(e, 2); series::stretch(t, 2);
	series::deriv_shift(t); series::integ_shift(t); series::integ_shift_offset(t, 1);
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
	series::with_len(e, 4); series::with_len(t, 2);
	series::maybe_cached<CE, true>{e};
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

TEST_CASE("series::vec cached wrappers", "[fft]") {
	using num = modnum<998244353>;
	using E = engines::ntt<num>;
	mt19937 mt(Catch::getSeed());
	// series::cached works with the cached fft:: entry points
	series::exact<E> a(37), b(21);
	fill_rnd(a, mt);
	fill_rnd(b, mt);
	series::cached_exact<E> ca(a), cb(b);
	series::exact<E> got(size_t(a.len() + b.len() - 1));
	fft::multiply<E>(span<const num>(ca), ca.cache(),
			span<const num>(cb), cb.cache(), span<num>(got));
	REQUIRE(got == a * b);
	REQUIRE((ca * cb) == (a * b));
	REQUIRE(middle_product(ca, cb) == fft::middle_product<E>(a, b));
	REQUIRE(square(ca) == square(a));
	// the same transform serves multiply and square of the same coefficients
	fft::transformed<E> fa, fb;
	series::exact<E> got2(size_t(a.len() + b.len() - 1));
	fft::multiply<E>(span<const num>(a), fa, span<const num>(b), fb, span<num>(got2));
	REQUIRE(got2 == a * b);
	fft::square<E>(span<const num>(a), fa, span<num>(got2));
	series::exact<E> asq = square(a);
	REQUIRE(equal(got2.begin(), got2.end(), asq.begin()));
	// cached_power_series products match plain products at all mixed shapes (see the
	// templated multiply_cached test for the transform-seeding path on all engines)
	series::trunc<E> pa(40), pb(25);
	fill_rnd(pa, mt);
	fill_rnd(pb, mt);
	series::prefix_cached<E> qa(pa), qb(pb);
	REQUIRE((qa * qb) == (pa * pb));
	REQUIRE((qa * pb) == (pa * pb));
	REQUIRE((pa * qb) == (pa * pb));
	// prefix caches survive precision extension: covered prefixes are reused, the
	// clamped full cache is rebuilt, and results still match
	series::trunc<E> tail(15);
	fill_rnd(tail, mt);
	qa.append(span<const num>(tail));
	series::trunc<E> pa2 = pa;
	pa2.insert(pa2.end(), tail.begin(), tail.end());
	for (int p : {8, 16, 32}) {
		auto pv = qa.first(min(p, qa.len()));
		REQUIRE(pv.len() == min(p, qa.len()));
		REQUIRE(pv[0] == pa2[0]);
	}
	REQUIRE((qa * qb) == (pa2 * pb));
	// products against many smaller operands reuse per-scale prefix caches
	for (int k : {1, 2, 3, 5, 17, 33, 100}) {
		series::trunc<E> small(size_t(k), num{});
		fill_rnd(small, mt);
		REQUIRE((qa * small) == (pa2 * small));
		REQUIRE((small * qa) == (small * pa2));
	}
	// a whole cache also serves truncated products when the operand fits under
	// the precision; oversized operands fall back to a truncated span
	series::cached_trunc<E> wt{series::trunc<E>(pa)};
	series::trunc<E> big(100);
	fill_rnd(big, mt);
	REQUIRE((wt * big) == (pa * big));
	REQUIRE((big * wt) == (big * pa));
	REQUIRE((wt * pb) == (pa * pb));
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

}} // namespace ecnerwala::fft
