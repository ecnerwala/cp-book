#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_get_random_seed.hpp>

#include "fft.hpp"

namespace ecnerwala {
namespace fft {

using namespace std;

template <typename T> vector<T> multiply_slow(const vector<T>& a, const vector<T>& b) {
	if (a.empty() || b.empty()) return {};
	vector<T> res(a.size() + b.size() - 1);
	for (int i = 0; i < int(a.size()); i++) {
		for (int j = 0; j < int(b.size()); j++) {
			res[i+j] += a[i] * b[j];
		}
	}
	return res;
}

// Small values for doubles so products are exactly representable; full range otherwise.
template <typename T> T rnd_val(mt19937& mt) {
	if constexpr (std::is_floating_point_v<T>) return T(int(mt() % 1024));
	else return T(mt());
}
template <typename T> void fill_rnd(vector<T>& v, mt19937& mt) {
	for (T& x : v) x = rnd_val<T>(mt);
}
template <typename T> void check_eq(const vector<T>& got, const vector<T>& want) {
	REQUIRE(got.size() == want.size());
	for (int i = 0; i < int(got.size()); i++) {
		INFO("i = " << i);
		if constexpr (std::is_floating_point_v<T>) REQUIRE(llround(got[i]) == llround(want[i]));
		else REQUIRE(got[i] == want[i]);
	}
}

#define ALL_ENGINES \
		fft_engine<modnum<998244353>>, fft_engine<mod_goldilocks>, fft_real_engine<double>, \
		fft_split_engine<modnum<int(1e9)+7>>, crt_engine<modnum<int(1e9)+7>>
#define MOD_ENGINES \
		fft_engine<modnum<998244353>>, fft_engine<mod_goldilocks>, \
		fft_split_engine<modnum<int(1e9)+7>>, crt_engine<modnum<int(1e9)+7>>

TEST_CASE("FFT core roundtrip", "[fft]") {
	using num = modnum<998244353>;
	mt19937 mt(Catch::getSeed());
	for (int n : {1, 2, 4, 8, 64, 256}) {
		vector<num> a(n);
		for (num& x : a) { x = num(mt()); }
		vector<num> t = a;
		fft_core<num>::forward(span<num>(t));
		fft_core<num>::inverse(span<num>(t));
		num d = inv(num(n));
		for (int i = 0; i < n; i++) REQUIRE(t[i] * d == a[i]);
	}
}

TEST_CASE("FFT core extend matches direct transform", "[fft]") {
	using num = modnum<998244353>;
	mt19937 mt(Catch::getSeed());
	for (int n : {1, 2, 8, 64}) {
		vector<num> a(n);
		for (num& x : a) { x = num(mt()); }
		// direct size-4n transform of a (zero-padded)
		vector<num> direct(4*n, num(0));
		copy(a.begin(), a.end(), direct.begin());
		fft_core<num>::forward(span<num>(direct));
		// grown by doubling
		vector<num> t = a;
		fft_core<num>::forward(span<num>(t));
		t.resize(2*n);
		fft_core<num>::extend(span<num>(t), span<const num>(a));
		t.resize(4*n);
		fft_core<num>::extend(span<num>(t), span<const num>(a));
		REQUIRE(t == direct);
	}
}

TEST_CASE("FFT core even/odd halves", "[fft]") {
	using num = modnum<998244353>;
	mt19937 mt(Catch::getSeed());
	for (int n : {2, 8, 64}) {
		vector<num> a(2*n);
		for (num& x : a) { x = num(mt()); }
		vector<num> evens(n), odds(n);
		for (int i = 0; i < n; i++) { evens[i] = a[2*i]; odds[i] = a[2*i+1]; }
		vector<num> t = a;
		fft_core<num>::forward(span<num>(t));
		vector<num> te(n), to(n);
		fft_core<num>::even_half(span<const num>(t), span<num>(te));
		fft_core<num>::odd_half(span<const num>(t), span<num>(to));
		fft_core<num>::forward(span<num>(evens));
		fft_core<num>::forward(span<num>(odds));
		REQUIRE(te == evens);
		REQUIRE(to == odds);
	}
}

TEST_CASE("FFT transform prefix is transform mod x^n - 1", "[fft]") {
	using num = modnum<998244353>;
	mt19937 mt(Catch::getSeed());
	int n = 16;
	vector<num> a(2*n);
	for (num& x : a) { x = num(mt()); }
	vector<num> t = a;
	fft_core<num>::forward(span<num>(t));
	vector<num> folded(n);
	for (int i = 0; i < n; i++) folded[i] = a[i] + a[i+n];
	fft_core<num>::forward(span<num>(folded));
	for (int i = 0; i < n; i++) REQUIRE(t[i] == folded[i]);
}

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

TEST_CASE("FFT double engine even/odd half", "[fft]") {
	using E = fft_real_engine<double>;
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
		auto th = odd ? E::odd_half(t, n) : E::even_half(t, n);
		auto& h = odd ? odds : evens;
		vector<double> got(n), want(n);
		E::finish(E::mul(th, tb, n), span<double>(got));
		multiply_circular<E>(span<const double>(h), span<const double>(b), span<double>(want), n);
		check_eq(got, want);
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
	// caller-owned coefficients + lazily built fft_cache transforms
	fft_cache<E> ca, cb, cc;
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
		fft_cache<E> cab;
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
		fft_cache<E> cd, ce;
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
		fft_cache<E> cd, ce;
		vector<num> de;
		fft_cache<E> cde;
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
		fft_cache<E> ca1, cb1, ca2, cb2;
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

TEMPLATE_TEST_CASE("FFT Inverse", "[fft]", MOD_ENGINES) {
	using E = TestType;
	using num = typename E::value_type;
	mt19937 mt(Catch::getSeed());
	vector<num> a(298);
	fill_rnd(a, mt);
	if (a[0] == 0) a[0] = 1;
	auto i = inverse<E>(a);
	auto r = multiply_slow(a, i);
	r.resize(a.size());
	vector<num> tgt(a.size());
	tgt[0] = 1;
	REQUIRE(r == tgt);
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

TEST_CASE("power_series log/exp/pow", "[fft]") {
	using num = modnum<998244353>;
	using ps = power_series<fft_engine<num>>;
	mt19937 mt(Catch::getSeed());
	for (int len : {1, 2, 3, 17, 100}) {
		INFO("len = " << len);
		ps a(len);
		for (num& x : a) { x = num(mt()); }
		a[0] = 1;
		auto l = poly_log(a);
		auto e = poly_exp(l);
		REQUIRE(e == a);

		// poly_pow vs repeated multiplication
		ps p3 = a * a;
		p3 *= a;
		ps q = a;
		q[0] = 1;
		REQUIRE(poly_pow(q, 3) == p3);
	}
	{
		// pow with valuation
		ps a(20, num(0));
		for (int i = 3; i < 20; i++) a[i] = num(mt());
		if (a[3] == 0) a[3] = 1;
		ps p2 = a * a;
		REQUIRE(poly_pow(a, 2) == p2);
		// valuation * exponent overflowing the length gives 0
		REQUIRE(poly_pow(a, 100) == ps(20, num(0)));
	}
}

TEST_CASE("power_series inverse", "[fft]") {
	using num = modnum<998244353>;
	using ps = power_series<fft_engine<num>>;
	mt19937 mt(Catch::getSeed());
	ps a(100);
	for (num& x : a) { x = num(mt()); }
	if (a[0] == 0) a[0] = 1;
	ps i = inverse(a);
	ps prod = a * i;
	ps tgt(a.size(), num(0));
	tgt[0] = 1;
	REQUIRE(prod == tgt);
}

TEST_CASE("power_series compose", "[fft]") {
	using num = modnum<998244353>;
	using ps = power_series<fft_engine<num>>;
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
		REQUIRE(poly_compose(f, g) == expected);
	}
}

TEST_CASE("poly_evaluate and poly_interpolate", "[fft]") {
	using num = modnum<998244353>;
	mt19937 mt(Catch::getSeed());
	for (int n : {1, 2, 3, 8, 17, 40}) {
		INFO("n = " << n);
		vector<num> poly(n);
		for (num& x : poly) { x = num(mt()); }
		vector<num> pts(n);
		for (int i = 0; i < n; i++) pts[i] = num(1000 + i);
		auto vals = poly_evaluate<fft_engine<num>>(poly, pts);
		for (int i = 0; i < n; i++) {
			num expect = 0;
			for (int j = n - 1; j >= 0; j--) expect = expect * pts[i] + poly[j];
			REQUIRE(vals[i] == expect);
		}
		auto rec = poly_interpolate<fft_engine<num>>(pts, vals);
		REQUIRE(rec == poly);
	}
}

TEMPLATE_TEST_CASE("online multiplier", "[fft]", MOD_ENGINES) {
	using E = TestType;
	using num = typename E::value_type;
	mt19937 mt(Catch::getSeed());
	for (int n : {1, 2, 3, 17, 64, 100}) {
		INFO("n = " << n);
		vector<num> f(n), g(n);
		fill_rnd(f, mt);
		fill_rnd(g, mt);
		auto slow = multiply_slow(f, g);
		slow.resize(2*n, num(0));
		online_multiplier<E> om(n);
		for (int i = 0; i < n; i++) {
			om.push(f[i], g[i]);
			REQUIRE(om.back() == slow[i]);
		}
		if ((n & (n-1)) == 0) {
			// the final push of a power-of-two N completes all 2N terms
			for (int i = n; i < 2*n; i++) {
				REQUIRE(om.res[i] == slow[i]);
			}
		}
	}
}

TEMPLATE_TEST_CASE("online squarer", "[fft]", MOD_ENGINES) {
	using E = TestType;
	using num = typename E::value_type;
	mt19937 mt(Catch::getSeed());
	for (int n : {1, 2, 3, 17, 64, 100}) {
		INFO("n = " << n);
		vector<num> f(n);
		fill_rnd(f, mt);
		auto slow = multiply_slow(f, f);
		slow.resize(2*n, num(0));
		online_squarer<E> os(n);
		for (int i = 0; i < n; i++) {
			os.push(f[i]);
			REQUIRE(os.back() == slow[i]);
		}
		if ((n & (n-1)) == 0) {
			for (int i = n; i < 2*n; i++) {
				REQUIRE(os.res[i] == slow[i]);
			}
		}
	}
}

TEST_CASE("poly_ap_values eval", "[fft,poly_ap_values]") {
	using num = modnum<998244353>;
	using poly_vals = poly_ap_values<fft_engine<num>>;
	mt19937 mt(Catch::getSeed());
	for (int len : {0, 1, 2, 3, 5, 8, 13, 21}) {
		INFO("len = " << len);
		std::vector<num> poly(len);
		for (int i = 0; i < len; i++) poly[i] = num(mt());
		auto eval_at = [&](num v) {
			num r = 0;
			for (int i = len-1; i >= 0; i--) r = r * v + poly[i];
			return r;
		};
		poly_vals v(len);
		for (int i = 0; i < len; i++) v[i] = eval_at(i);
		for (int i = -2 * len; i <= 2 * len; i++) {
			REQUIRE(v.eval_at(i) == eval_at(i));
		}
		auto eval_range = [&](num k, int osz) {
			poly_vals r(osz);
			for (int i = 0; i < osz; i++) {
				r[i] = eval_at(k + num(i));
			}
			return r;
		};
		num k = 1023895;
		for (int osz : {0, 1, 2, 3, 5, 8, 13, 21}) {
			INFO("osz = " << osz);
			REQUIRE(v.eval_range(k, osz) == eval_range(k, osz));
		}
	}
}

}} // namespace ecnerwala::fft
