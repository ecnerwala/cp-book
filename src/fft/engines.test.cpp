#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_get_random_seed.hpp>

#include "fft/engines/ntt.hpp"
#include "fft/engines/real.hpp"
#include "fft/engines/split.hpp"
#include "fft/engines/crt.hpp"
#include "fft/multiply.hpp"
#include "fft/test_util.test.hpp"
#include "num/modnum.hpp"

namespace wala {
namespace fft {

using namespace std;

// engine concept sanity checks (the base engines; the wrapping algebras are
// checked in engines/algebras.test.cpp)
static_assert(engine<engines::ntt<modnum<998244353>>>);
static_assert(engine<engines::ntt<mod_goldilocks>>);
static_assert(engine<engines::real<double>>);
static_assert(engine<engines::split<modnum<int(1e9)+7>>>);
static_assert(engine<engines::crt<modnum<int(1e9)+7>>>);

TEMPLATE_TEST_CASE("FFT multiply sizes", "[fft]", ALL_ENGINES) {
	using E = TestType;
	using num = typename E::value_type;
	mt19937 mt(Catch::getSeed());
	for (int la : {1, 2, 3, 17, 100}) {
		for (int lb : {1, 2, 16, 17, 168}) {
			vector<num> a(la), b(lb);
			fill_rnd(a, mt);
			fill_rnd(b, mt);
			INFO("la = " << la << ", lb = " << lb);
			check_eq(multiply<E>(a, b), multiply_slow(a, b));
		}
	}
}

TEMPLATE_TEST_CASE("transform add", "[fft]", MOD_ENGINES) {
	// (a+b)*c via a transform-domain sum; for the inexact engines this exercises the
	// scale-tracked transformed_t (add gives A=2, mul gives product<2>)
	using E = TestType;
	using num = typename E::value_type;
	mt19937 mt(Catch::getSeed());
	for (int n : {2, 16, 64}) {
		vector<num> a(n), b(n), c(n);
		fill_rnd(a, mt); fill_rnd(b, mt); fill_rnd(c, mt);
		auto tab = E::add(E::transform(span<const num>(a), n), E::transform(span<const num>(b), n));
		auto tc = E::transform(span<const num>(c), n);
		vector<num> got(n);
		E::finish(E::mul(tab, tc, n), span<num>(got));
		vector<num> ab(n);
		for (int i = 0; i < n; i++) ab[i] = a[i] + b[i];
		auto want = multiply_slow(ab, c);
		// compare mod x^n - 1 (circular)
		for (int i = n; i < sz(want); i++) want[i - n] += want[i];
		want.resize(n);
		check_eq(got, want);
	}
	// scale conversions: widening is implicit, narrowing is an explicit downcast
	if constexpr (E::unit_scale != 0) {
		using T1 = typename E::template transformed_t<1>;
		using T2 = typename E::template transformed_t<2>;
		STATIC_REQUIRE(std::is_convertible_v<T1&&, T2>);
		STATIC_REQUIRE(!std::is_convertible_v<T2&&, T1>);
		STATIC_REQUIRE(std::is_constructible_v<T1, T2&&>);
		using P1 = typename E::template product_t<1>;
		using P2 = typename E::template product_t<2>;
		STATIC_REQUIRE(std::is_convertible_v<P1&&, P2>);
		STATIC_REQUIRE(!std::is_convertible_v<P2&&, P1>);
		STATIC_REQUIRE(std::is_constructible_v<P1, P2&&>);
	}
}

TEST_CASE("FFT double engine downsample", "[fft]") {
	using E = engines::real<double>;
	mt19937 mt(Catch::getSeed());
	int n = 64;
	vector<double> a(2*n), evens(n), odds(n);
	for (int i = 0; i < 2*n; i++) {
		a[i] = rnd_val<double>(mt);
		(i & 1 ? odds : evens)[i/2] = a[i];
	}
	auto t = E::transform(span<const double>(a), 2*n);
	vector<double> b(n);
	fill_rnd(b, mt);
	auto tb = E::transform(span<const double>(b), n);
	for (bool odd : {false, true}) {
		auto th = E::downsample(t, n, odd);
		auto& h = odd ? odds : evens;
		vector<double> got(n), want(n);
		E::finish(E::mul(th, tb, n), span<double>(got));
		multiply_circular<E>(span<const double>(h), span<const double>(b), span<double>(want), n);
		check_eq(got, want);
	}
}

TEMPLATE_TEST_CASE("product downsample", "[fft]", ALL_ENGINES) {
	using E = TestType;
	using num = typename E::value_type;
	mt19937 mt(Catch::getSeed());
	for (int n : {2, 16, 64}) {
		vector<num> a(2*n), b(2*n);
		fill_rnd(a, mt);
		fill_rnd(b, mt);
		auto p = E::mul(E::transform(span<const num>(a), 2*n), E::transform(span<const num>(b), 2*n), 2*n);
		vector<num> full(2*n);
		E::finish(auto(p), span<num>(full));
		for (bool odd : {false, true}) {
			auto ph = E::downsample(p, n, odd);
			vector<num> got(n);
			E::finish(std::move(ph), span<num>(got));
			vector<num> want(n);
			for (int i = 0; i < n; i++) want[i] = full[2*i + odd];
			check_eq(got, want);
		}
	}
}

TEMPLATE_TEST_CASE("transform upsample", "[fft]", ALL_ENGINES) {
	using E = TestType;
	using num = typename E::value_type;
	mt19937 mt(Catch::getSeed());
	for (int n : {2, 16, 64}) {
		vector<num> a(n), b(2*n);
		fill_rnd(a, mt);
		fill_rnd(b, mt);
		auto ta = E::transform(span<const num>(a), n);
		auto tb = E::transform(span<const num>(b), 2*n);
		for (bool odd : {false, true}) {
			// multiplying the upsampled transform recovers the product with the spread input
			auto tu = E::upsample(ta, 2*n, odd);
			vector<num> spread(2*n);
			for (int i = 0; i < n; i++) spread[2*i + odd] = a[i];
			vector<num> got(2*n), want(2*n);
			E::finish(E::mul(tu, tb, 2*n), span<num>(got));
			multiply_circular<E>(span<const num>(spread), span<const num>(b), span<num>(want), 2*n);
			check_eq(got, want);
			// round trip: downsample recovers the input transform, the other parity is zero
			// (tb's size-n prefix is the transform of b wrapped mod x^n - 1)
			vector<num> bw(n);
			for (int i = 0; i < n; i++) bw[i] = b[i] + b[n + i];
			vector<num> rt_got(n), rt_want(n), zero(n);
			E::finish(E::mul(E::downsample(tu, n, odd), tb, n), span<num>(rt_got));
			multiply_circular<E>(span<const num>(a), span<const num>(bw), span<num>(rt_want), n);
			check_eq(rt_got, rt_want);
			E::finish(E::mul(E::downsample(tu, n, !odd), tb, n), span<num>(rt_got));
			check_eq(rt_got, zero);
		}
	}
}

TEMPLATE_TEST_CASE("product upsample", "[fft]", ALL_ENGINES) {
	using E = TestType;
	using num = typename E::value_type;
	mt19937 mt(Catch::getSeed());
	for (int n : {2, 16, 64}) {
		vector<num> a(n), b(n);
		fill_rnd(a, mt);
		fill_rnd(b, mt);
		auto p = E::mul(E::transform(span<const num>(a), n), E::transform(span<const num>(b), n), n);
		vector<num> full(n);
		E::finish(auto(p), span<num>(full));
		for (bool odd : {false, true}) {
			auto pu = E::upsample(p, 2*n, odd);
			vector<num> got(2*n);
			E::finish(std::move(pu), span<num>(got));
			vector<num> want(2*n);
			for (int i = 0; i < n; i++) want[2*i + odd] = full[i];
			check_eq(got, want);
		}
	}
}

TEMPLATE_TEST_CASE("graeffe", "[fft]", ALL_ENGINES) {
	using E = TestType;
	using num = typename E::value_type;
	mt19937 mt(Catch::getSeed());
	for (int n : {2, 16, 64}) {
		vector<num> a(n);
		fill_rnd(a, mt);
		auto ta = E::transform(span<const num>(a), n);
		vector<num> got(n/2);
		E::finish(E::graeffe(ta, n), span<num>(got));
		// even part of a(x) * a(-x) mod x^n - 1
		vector<num> b = a;
		for (int i = 1; i < n; i += 2) b[i] = -b[i];
		auto full = multiply_slow(a, b);
		full.resize(size_t(2*n), num(0));
		vector<num> want(n/2);
		for (int i = 0; i < n/2; i++) want[i] = full[2*i] + full[2*i + n];
		check_eq(got, want);
	}
}

TEMPLATE_TEST_CASE("padded transform and finish", "[fft]", ALL_ENGINES) {
	using E = TestType;
	using num = typename E::value_type;
	mt19937 mt(Catch::getSeed());
	for (int n : {4, 16, 64}) {
		for (int C : {2, 4, n}) {
			INFO("n = " << n << ", C = " << C);
			vector<num> a(n), b(n);
			fill_rnd(b, mt);
			for (int i = 0; i < n; i++) {
				if (i % C < C/2) a[i] = rnd_val<num>(mt);
			}
			auto ta = E::transform_padded(span<const num>(a), n, C);
			auto tb = E::transform(span<const num>(b), n);
			vector<num> want(n);
			multiply_circular<E>(span<const num>(a), span<const num>(b), span<num>(want), n);
			// transform_padded matches the plain transform
			vector<num> got(n);
			E::finish(E::mul(ta, tb, n), span<num>(got));
			check_eq(got, want);
			// finish_padded writes the live positions and leaves the rest unspecified
			vector<num> got_padded(n);
			E::finish_padded(E::mul(ta, tb, n), span<num>(got_padded), C);
			for (int i = 0; i < n; i++) {
				if (i % C >= C/2) { got_padded[i] = num(0); want[i] = num(0); }
			}
			check_eq(got_padded, want);
			// and composes with a non-assign op on the live positions
			vector<num> acc(n), expected(n);
			fill_rnd(acc, mt);
			for (int i = 0; i < n; i++) {
				expected[i] = i % C < C/2 ? acc[i] + want[i] : acc[i];
			}
			E::finish_padded(E::mul(ta, tb, n), span<num>(acc), C, add_op{});
			for (int i = 0; i < n; i++) {
				if (i % C >= C/2) acc[i] = expected[i];
			}
			check_eq(acc, expected);
		}
	}
}

TEMPLATE_TEST_CASE("FFT Square", "[fft]", ALL_ENGINES) {
	using E = TestType;
	using num = typename E::value_type;
	mt19937 mt(Catch::getSeed());
	for (int la : {1, 3, 17, 100}) {
		INFO("la = " << la);
		vector<num> a(la);
		fill_rnd(a, mt);
		check_eq(square<E>(a), multiply_slow(a, a));
	}
}

TEMPLATE_TEST_CASE("FFT multiply with add-into op", "[fft]", ALL_ENGINES) {
	using E = TestType;
	using num = typename E::value_type;
	mt19937 mt(Catch::getSeed());
	vector<num> a(37), b(23);
	fill_rnd(a, mt);
	fill_rnd(b, mt);
	vector<num> out(a.size() + b.size() - 1);
	fill_rnd(out, mt);
	vector<num> expected = out;
	auto slow = multiply_slow(a, b);
	for (int i = 0; i < int(slow.size()); i++) expected[i] += slow[i];
	multiply<E>(span<const num>(a), span<const num>(b), span<num>(out), add_op{});
	check_eq(out, expected);
}

TEMPLATE_TEST_CASE("FFT cached multiply", "[fft]", ALL_ENGINES) {
	using E = TestType;
	using num = typename E::value_type;
	mt19937 mt(Catch::getSeed());
	vector<num> a(37), b(23), c(100);
	fill_rnd(a, mt);
	fill_rnd(b, mt);
	fill_rnd(c, mt);
	// caller-owned coefficients + lazily built transforms
	fft::transformed<E> ca, cb, cc;
	{
		vector<num> out(a.size() + b.size() - 1);
		multiply<E>(span<const num>(a), ca, span<const num>(b), cb, span<num>(out));
		check_eq(out, multiply_slow(a, b));
	}
	{
		// reuse of a's transform at a larger size exercises extend
		vector<num> out(a.size() + c.size() - 1);
		multiply<E>(span<const num>(a), ca, span<const num>(c), cc, span<num>(out));
		check_eq(out, multiply_slow(a, c));
	}
	{
		// and shrink back down (prefix use)
		vector<num> out(a.size() + b.size() - 1);
		multiply<E>(span<const num>(a), ca, span<const num>(b), cb, span<num>(out));
		check_eq(out, multiply_slow(a, b));
	}
	{
		vector<num> out(2 * a.size() - 1);
		square<E>(span<const num>(a), ca, span<num>(out));
		check_eq(out, multiply_slow(a, a));
	}
	{
		// multiply_cached: coefficients match, and the seeded (or lazily built)
		// transform is directly usable in further products, including after extend
		vector<num> ab;
		fft::transformed<E> cab;
		multiply_cached<E>(span<const num>(a), ca, span<const num>(b), cb, ab, cab);
		check_eq(ab, multiply_slow(a, b));
		vector<num> out(ab.size() + c.size() - 1);
		multiply<E>(span<const num>(ab), cab, span<const num>(c), cc, span<num>(out));
		check_eq(out, multiply_slow(ab, c));
	}
	{
		// 2^k+1 by 1: the cut fires with the long operand's folded size-2^k transform
		vector<num> d(33), e(1);
		fill_rnd(d, mt);
		fill_rnd(e, mt);
		fft::transformed<E> cd, ce;
		vector<num> out(33);
		multiply<E>(span<const num>(d), cd, span<const num>(e), ce, span<num>(out));
		check_eq(out, multiply_slow(d, e));
		// the folded transform must still extend correctly for a bigger product
		vector<num> out2(d.size() + c.size() - 1);
		multiply<E>(span<const num>(d), cd, span<const num>(c), cc, span<num>(out2));
		check_eq(out2, multiply_slow(d, c));
	}
	{
		// 2^k+1 shapes hit the cut, seeding a folded transform; extending it must fold
		vector<num> d(33), e(33);
		fill_rnd(d, mt);
		fill_rnd(e, mt);
		fft::transformed<E> cd, ce;
		vector<num> de;
		fft::transformed<E> cde;
		multiply_cached<E>(span<const num>(d), cd, span<const num>(e), ce, de, cde);
		check_eq(de, multiply_slow(d, e));
		vector<num> out(de.size() + c.size() - 1);
		multiply<E>(span<const num>(de), cde, span<const num>(c), cc, span<num>(out));
		check_eq(out, multiply_slow(de, c));
	}
}

TEMPLATE_TEST_CASE("FFT multiply_add2", "[fft]", ALL_ENGINES) {
	using E = TestType;
	using num = typename E::value_type;
	mt19937 mt(Catch::getSeed());
	// pairs share a linear length; includes 2^k+1 cut shapes and length-1 operands
	vector<array<int, 4>> cases = {
		{1, 2, 2, 1}, {2, 2, 2, 2}, {3, 3, 3, 3}, {1, 5, 3, 3},
		{4, 5, 8, 1}, {17, 16, 32, 1}, {65, 65, 129, 1}, {100, 29, 64, 65},
	};
	for (auto [la1, lb1, la2, lb2] : cases) {
		vector<num> a1(la1), b1(lb1), a2(la2), b2(lb2);
		fill_rnd(a1, mt); fill_rnd(b1, mt); fill_rnd(a2, mt); fill_rnd(b2, mt);
		vector<num> want = multiply_slow(a1, b1);
		vector<num> p2 = multiply_slow(a2, b2);
		for (int i = 0; i < sz(p2); i++) want[i] += p2[i];
		fft::transformed<E> ca1, cb1, ca2, cb2;
		vector<num> got(want.size());
		multiply_add2<E>(span<const num>(a1), ca1, span<const num>(b1), cb1,
				span<const num>(a2), ca2, span<const num>(b2), cb2, span<num>(got));
		INFO("la1 = " << la1 << ", lb1 = " << lb1 << ", la2 = " << la2 << ", lb2 = " << lb2);
		check_eq(got, want);
	}
}

TEMPLATE_TEST_CASE("FFT middle product", "[fft]", ALL_ENGINES) {
	using E = TestType;
	using num = typename E::value_type;
	mt19937 mt(Catch::getSeed());
	for (int la : {5, 16, 17, 33, 100}) {
		for (int lb : {1, 3, 5}) {
			if (lb > la) continue;
			vector<num> a(la), b(lb);
			fill_rnd(a, mt);
			fill_rnd(b, mt);
			auto full = multiply_slow(a, b);
			vector<num> expected(full.begin() + (lb - 1), full.begin() + la);
			INFO("la = " << la << ", lb = " << lb);
			check_eq(middle_product<E>(a, b), expected);
			// out-span form accumulates through the op
			vector<num> acc(expected.size(), num(1));
			middle_product<E>(span<const num>(a), span<const num>(b), span<num>(acc), fft::add_op{});
			for (auto& v : acc) v -= num(1);
			check_eq(acc, expected);
		}
	}
	// equal lengths: a single dot product
	vector<num> a(9), b(9);
	fill_rnd(a, mt);
	fill_rnd(b, mt);
	auto full = multiply_slow(a, b);
	check_eq(middle_product<E>(a, b), vector<num>(1, full[8]));
}

TEMPLATE_TEST_CASE("negate_arg transforms", "[fft]", ALL_ENGINES) {
	using E = TestType;
	using num = typename E::value_type;
	mt19937 mt(Catch::getSeed());
	vector<num> delta(1, num(1));
	for (int n : {2, 4, 16, 64}) {
		vector<num> a(n);
		fill_rnd(a, mt);
		auto ta = E::transform(span<const num>(a), n);
		auto td = E::transform(span<const num>(delta), n);
		// negate_arg: multiplying by delta's transform recovers coefficients of a(-x)
		vector<num> got(n), want = a;
		E::finish(E::mul(E::negate_arg(ta, n), td, n), span<num>(got));
		for (int i = 1; i < n; i += 2) want[i] = -want[i];
		check_eq(got, want);
	}
}

}} // namespace wala::fft
