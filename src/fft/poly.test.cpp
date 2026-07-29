#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_get_random_seed.hpp>

#include "fft/poly.hpp"
#include "fft/series.hpp"
#include "fft/test_util.test.hpp"
#include "modnum.hpp"

namespace ecnerwala {
namespace fft {

using namespace std;

TEST_CASE("poly::form evaluation and transposed multiplication", "[fft]") {
	using num = modnum<998244353>;
	using E = engines::ntt<num>;
	mt19937 mt(Catch::getSeed());
	int n = 40;
	num z = num(mt());
	auto f = poly::form<E>::polynomial_evaluation(z, n);
	vector<num> sv(30);
	fill_rnd(sv, mt);
	poly::vec<E> s((span<const num>(sv)));
	REQUIRE(f(s) == s(z));
	// f(S * q) == composed_with(q)(S)
	vector<num> qv(11), s2v(n - 11);
	fill_rnd(qv, mt);
	fill_rnd(s2v, mt);
	poly::vec<E> q((span<const num>(qv))), s2((span<const num>(s2v)));
	auto fq = f.composed_with(q);
	REQUIRE(fq(s2) == f(s2 * q));
	// evaluation functional composed with q evaluates S * q at z
	REQUIRE(fq(s2) == (s2 * q)(z));
	// composing with a power series (living in 1/x) multiplies the storages,
	// prefix-truncated back to length n
	series::trunc<E> t(size_t(n), num{});
	fill_rnd(t, mt);
	auto ft = f.composed_with(t);
	REQUIRE(ft.len() == n);
	for (int j = 0; j < n; j++) {
		num want{};
		for (int d = 0; d <= j; d++) want += t[d] * f.rev_series()[j - d];
		REQUIRE(ft.rev_series()[j] == want);
	}
	// an exact series may be any length; the tail beyond it is zero
	series::exact<E> e(t.begin(), t.begin() + 11);
	auto fe = f.composed_with(e);
	REQUIRE(fe.len() == n);
	for (int j = 0; j < n; j++) {
		num want{};
		for (int d = 0; d <= j && d < 11; d++) want += t[d] * f.rev_series()[j - d];
		REQUIRE(fe.rev_series()[j] == want);
	}
}

TEST_CASE("poly::vec reversed storage and series interop", "[fft]") {
	using num = modnum<998244353>;
	using E = engines::ntt<num>;
	mt19937 mt(Catch::getSeed());
	vector<num> pa(37), pb(23);
	fill_rnd(pa, mt);
	fill_rnd(pb, mt);
	poly::vec<E> a((span<const num>(pa))), b((span<const num>(pb)));
	// indexing is coefficient order, storage is reversed
	REQUIRE(a[0] == pa[0]);
	REQUIRE(a.leading() == pa[36]);
	REQUIRE(a.rev_series()[0] == pa[36]);
	REQUIRE(a.rev_series()[36] == pa[0]);
	// products convolve the storage directly
	poly::vec<E> p = a * b;
	check_eq(vector<num>(p.begin(), p.end()), multiply_slow(pa, pb));
	REQUIRE(square(a) == a * a);
	num x = num(mt());
	REQUIRE(p(x) == a(x) * b(x));
	// +/- align at x^0 (the shared storage tail)
	poly::vec<E> s = a + b, d = b - a;
	REQUIRE(s.len() == 37);
	for (int i = 0; i < 37; i++) REQUIRE(s[i] == pa[i] + (i < 23 ? pb[i] : num(0)));
	for (int i = 0; i < 37; i++) REQUIRE(d[i] == (i < 23 ? pb[i] : num(0)) - pa[i]);
	// multiplying by x^k appends zeros to the storage; coefficients shift up
	poly::vec<E> g = a;
	g.shift(2);
	REQUIRE(g.len() == 39);
	REQUIRE(g[0] == num(0));
	REQUIRE(g[1] == num(0));
	for (int i = 0; i < 37; i++) REQUIRE(g[i + 2] == pa[i]);
	REQUIRE(g.rev_series().data()[0] == pa[36]);
	// named conversions use the reversed convention and round-trip freely
	const series::exact<E>& ra = a.rev_series();
	REQUIRE(ra.len() == 37);
	for (int i = 0; i < 37; i++) REQUIRE(ra[i] == pa[36 - i]);
	REQUIRE(poly::vec<E>::from_rev_series(ra) == a);
	// a poly::vec's natural-order coefficients as an exact series (for series products)
	series::exact<E> xa(a.begin(), a.end());
	REQUIRE(equal(xa.begin(), xa.end(), pa.begin(), pa.end()));
	REQUIRE(a.unrev_series(10) == series::trunc<E>(pa.begin(), pa.begin() + 10));
	// the storage transform serves transposed products: middle product against rev_series()
	std::vector<num> vals(60);
	fill_rnd(vals, mt);
	series::cached_exact<E> cv(series::exact<E>(vals.begin(), vals.end()));
	series::cached_exact<E> ca(a.rev_series());
	auto mp = middle_product(cv, ca);
	auto naive = [&](int j) {
		num r{};
		for (int t = 0; t < 37; t++) r += pa[t] * vals[j + t];
		return r;
	};
	for (int j = 0; j < sz(mp); j++) REQUIRE(mp[size_t(j)] == naive(j));
}

TEST_CASE("poly::cached products", "[fft]") {
	using num = modnum<998244353>;
	using E = engines::ntt<num>;
	mt19937 mt(Catch::getSeed());
	vector<num> pa(37), pb(23);
	fill_rnd(pa, mt);
	fill_rnd(pb, mt);
	poly::vec<E> a((span<const num>(pa))), b((span<const num>(pb)));
	// poly::vec products return poly::cached, adopting the product transform
	auto p = a * b;
	static_assert(std::is_same_v<decltype(p), poly::cached<E>>);
	REQUIRE(p.rev_series().cache().size() > 0);
	poly::vec<E> pp = a * b; // naming the plain type moves out and drops the transform
	check_eq(vector<num>(pp.begin(), pp.end()), multiply_slow(pa, pb));
	// cached operands reuse and chain; results compare across representations
	poly::cached<E> ca(a), cb(b);
	REQUIRE(ca == a);
	REQUIRE(ca * cb == p);
	REQUIRE(ca * b == p);
	auto sq = square(ca);
	static_assert(std::is_same_v<decltype(sq), poly::cached<E>>);
	REQUIRE(sq == a * a);
	num x = num(mt());
	REQUIRE(p(x) == ca(x) * cb(x));
	// moving out drops down to a plain mutable poly::vec
	poly::vec<E> q = std::move(ca);
	REQUIRE(q == a);
}

TEST_CASE("poly::multipoint and poly::interpolate", "[fft]") {
	using num = modnum<998244353>;
	mt19937 mt(Catch::getSeed());
	for (int n : {1, 2, 3, 8, 17, 40}) {
		INFO("n = " << n);
		vector<num> coeffs(n);
		for (num& x : coeffs) { x = num(mt()); }
		poly::vec<engines::ntt<num>> p((span<const num>(coeffs)));
		vector<num> pts(n);
		for (int i = 0; i < n; i++) pts[i] = num(1000 + i);
		auto vals = poly::multipoint<engines::ntt<num>>(p, pts);
		for (int i = 0; i < n; i++) {
			REQUIRE(vals[i] == p(pts[i]));
		}
		auto rec = poly::interpolate<engines::ntt<num>>(pts, vals);
		REQUIRE(rec == p);
	}
}

}
} // namespace ecnerwala::fft
