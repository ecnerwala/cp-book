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

template <typename E, bool online, int N>
void test_matrix_engine(mt19937& mt) {
	using M = typename E::value_type;
	using num = std::remove_reference_t<decltype(std::declval<M&>()[{0, 0}])>;
	auto rnd_mat = [&]() {
		M m;
		for (int r = 0; r < N; r++) for (int c = 0; c < N; c++) m[{r, c}] = rnd_val<num>(mt);
		return m;
	};
	for (int la : {1, 2, 3, 17, 33}) {
		for (int lb : {1, 2, 16, 17}) {
			vector<M> a(la), b(lb);
			for (M& m : a) m = rnd_mat();
			for (M& m : b) m = rnd_mat();
			INFO("la = " << la << ", lb = " << lb);
			REQUIRE(multiply<E>(a, b) == multiply_slow(a, b));
		}
	}
	// square of a matrix sequence must keep both cross orders
	int n = 33;
	vector<M> f(n);
	for (M& m : f) m = rnd_mat();
	auto slow = multiply_slow(f, f);
	slow.resize(2*n, M{});
	vector<M> got(2*n - 1);
	square<E>(span<const M>(f), span<M>(got));
	REQUIRE(got == vector<M>(slow.begin(), slow.begin() + 2*n - 1));
	if constexpr (online) {
		online_squarer<E> os(n);
		for (int i = 0; i < n; i++) {
			os.push(f[i]);
			REQUIRE(os.back() == slow[i]);
		}
	}
}

TEMPLATE_TEST_CASE("matrix engine", "[fft]",
		fft_engine<modnum<998244353>>,
		fft_split_engine<modnum<int(1e9)+7>>,
		crt_engine<modnum<int(1e9)+7>>) {
	using IE = TestType;
	// the tracked engines' scale budget admits N = 2 (entries are N-addend sums), and
	// the non-commutative online squarer accumulates two N-addend products per window
	// (scale 2N), exceeding it
	constexpr int N = IE::unit_scale == 0 ? 3 : 2;
	mt19937 mt(Catch::getSeed());
	test_matrix_engine<matrix_engine<IE, N>, IE::unit_scale == 0, N>(mt);
	// the stable variant works at any N (and its online squarer stays at scale 2)
	test_matrix_engine<matrix_engine_stable<IE, 3>, true, 3>(mt);
}

template <typename E, typename num, int N>
void test_trunc_series_engine(mt19937& mt) {
	using P = typename E::value_type;
	auto rnd_p = [&]() {
		P p;
		for (int i = 0; i < N; i++) p[i] = rnd_val<num>(mt);
		return p;
	};
	for (int la : {1, 2, 3, 17, 33}) {
		for (int lb : {1, 2, 16, 17}) {
			vector<P> a(la), b(lb);
			for (P& p : a) p = rnd_p();
			for (P& p : b) p = rnd_p();
			INFO("la = " << la << ", lb = " << lb);
			REQUIRE(multiply<E>(a, b) == multiply_slow(a, b));
		}
	}
}

TEMPLATE_TEST_CASE("trunc_series engine", "[fft]",
		fft_engine<modnum<998244353>>,
		fft_split_engine<modnum<int(1e9)+7>>,
		crt_engine<modnum<int(1e9)+7>>) {
	using IE = TestType;
	using num = typename IE::value_type;
	constexpr int N = IE::unit_scale == 0 ? 3 : 2;
	mt19937 mt(Catch::getSeed());
	test_trunc_series_engine<trunc_series_engine<IE, N>, num, N>(mt);
	test_trunc_series_engine<trunc_series_engine_stable<IE, 3>, num, 3>(mt);
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

TEMPLATE_TEST_CASE("Bostan-Mori kth_term", "[fft]", MOD_ENGINES) {
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
		power_series_exact<E> xp(p.begin(), p.end()), xq(q.begin(), q.end());
		for (uint64_t k : {uint64_t(0), uint64_t(1), uint64_t(7), uint64_t(100), uint64_t(299)}) {
			INFO("d = " << d << ", k = " << k);
			REQUIRE(kth_term<E>(xp, xq, k) == ser[k]);
		}
	}
}

TEST_CASE("power_series cached wrappers", "[fft]") {
	using num = modnum<998244353>;
	using E = fft_engine<num>;
	mt19937 mt(Catch::getSeed());
	// cached_power_series_exact works with the cached fft:: entry points
	power_series_exact<E> a(37), b(21);
	fill_rnd(a, mt);
	fill_rnd(b, mt);
	cached_power_series_exact<E> ca(a), cb(b);
	power_series_exact<E> got(size_t(a.len() + b.len() - 1));
	fft::multiply<E>(span<const num>(ca.underlying()), ca.cache(),
			span<const num>(cb.underlying()), cb.cache(), span<num>(got));
	REQUIRE(got == a * b);
	REQUIRE((ca * cb) == (a * b));
	REQUIRE(middle_product(ca, cb) == fft::middle_product<E>(a, b));
	REQUIRE(square(ca) == square(a));
	// the same fft_cache serves multiply and square of the same coefficients
	fft::fft_cache<E> fa, fb;
	power_series_exact<E> got2(size_t(a.len() + b.len() - 1));
	fft::multiply<E>(span<const num>(a), fa, span<const num>(b), fb, span<num>(got2));
	REQUIRE(got2 == a * b);
	fft::square<E>(span<const num>(a), fa, span<num>(got2));
	power_series_exact<E> asq = square(a);
	REQUIRE(equal(got2.begin(), got2.end(), asq.begin()));
	// cached_power_series products match plain products at all mixed shapes (see the
	// templated multiply_cached test for the transform-seeding path on all engines)
	power_series<E> pa(40), pb(25);
	fill_rnd(pa, mt);
	fill_rnd(pb, mt);
	prefix_cached_power_series<E> qa(pa), qb(pb);
	REQUIRE((qa * qb) == (pa * pb));
	REQUIRE((qa * pb) == (pa * pb));
	REQUIRE((pa * qb) == (pa * pb));
	// prefix caches survive precision extension: covered prefixes are reused, the
	// clamped full cache is rebuilt, and results still match
	power_series<E> tail(15);
	fill_rnd(tail, mt);
	qa.append(span<const num>(tail));
	power_series<E> pa2 = pa;
	pa2.insert(pa2.end(), tail.begin(), tail.end());
	for (int p : {8, 16, 32, 64}) {
		auto& pc = qa.prefix_cache(p);
		REQUIRE(pc.len() == min(p + 1, qa.len()));
		REQUIRE(qa.prefix(p)[0] == pa2[0]);
	}
	REQUIRE((qa * qb) == (pa2 * pb));
	// products against many smaller operands reuse per-scale prefix caches
	for (int k : {1, 2, 3, 5, 17, 33, 100}) {
		power_series<E> small(size_t(k), num{});
		fill_rnd(small, mt);
		REQUIRE((qa * small) == (pa2 * small));
		REQUIRE((small * qa) == (small * pa2));
	}
	// exact cached series give full products, mixing with truncated operands
	prefix_cached_power_series<E, true> xqa(a), xqb(b);
	REQUIRE((xqa * xqb) == (a * b));
	REQUIRE((xqa * b) == (a * b));
	REQUIRE((xqa * qb) == (a * pb));
	REQUIRE((pb * xqa) == (pb * a));
}

TEST_CASE("linear_form evaluation and transposed multiplication", "[fft]") {
	using num = modnum<998244353>;
	using E = fft_engine<num>;
	mt19937 mt(Catch::getSeed());
	int n = 40;
	num z = num(mt());
	auto f = linear_form<E>::polynomial_evaluation(z, n);
	vector<num> sv(30);
	fill_rnd(sv, mt);
	poly<E> s((span<const num>(sv)));
	REQUIRE(f(s) == s(z));
	// f(S * q) == composed_with(q)(S)
	vector<num> qv(11), s2v(n - 11);
	fill_rnd(qv, mt);
	fill_rnd(s2v, mt);
	poly<E> q((span<const num>(qv))), s2((span<const num>(s2v)));
	auto fq = f.composed_with(q);
	REQUIRE(fq(s2) == f(s2 * q));
	// evaluation functional composed with q evaluates S * q at z
	REQUIRE(fq(s2) == (s2 * q)(z));
	// composing with a power series (living in 1/x) multiplies the storages,
	// prefix-truncated back to length n
	power_series<E> t(size_t(n), num{});
	fill_rnd(t, mt);
	auto ft = f.composed_with(t);
	REQUIRE(ft.len() == n);
	for (int j = 0; j < n; j++) {
		num want{};
		for (int d = 0; d <= j; d++) want += t[d] * f.rev_series()[j - d];
		REQUIRE(ft.rev_series()[j] == want);
	}
	// an exact series may be any length; the tail beyond it is zero
	power_series_exact<E> e(t.begin(), t.begin() + 11);
	auto fe = f.composed_with(e);
	REQUIRE(fe.len() == n);
	for (int j = 0; j < n; j++) {
		num want{};
		for (int d = 0; d <= j && d < 11; d++) want += t[d] * f.rev_series()[j - d];
		REQUIRE(fe.rev_series()[j] == want);
	}
}

TEST_CASE("power_series mixed exactness operators", "[fft]") {
	using num = modnum<998244353>;
	using E = fft_engine<num>;
	using xps = power_series_exact<E>;
	using ps = power_series<E>;
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

TEST_CASE("poly reversed storage and series interop", "[fft]") {
	using num = modnum<998244353>;
	using E = fft_engine<num>;
	mt19937 mt(Catch::getSeed());
	vector<num> pa(37), pb(23);
	fill_rnd(pa, mt);
	fill_rnd(pb, mt);
	poly<E> a((span<const num>(pa))), b((span<const num>(pb)));
	// indexing is coefficient order, storage is reversed
	REQUIRE(a[0] == pa[0]);
	REQUIRE(a.leading() == pa[36]);
	REQUIRE(a.rev_series()[0] == pa[36]);
	REQUIRE(a.rev_series()[36] == pa[0]);
	// products convolve the storage directly
	poly<E> p = a * b;
	check_eq(vector<num>(p.begin(), p.end()), multiply_slow(pa, pb));
	REQUIRE(square(a) == a * a);
	num x = num(mt());
	REQUIRE(p(x) == a(x) * b(x));
	// +/- align at x^0 (the shared storage tail)
	poly<E> s = a + b, d = b - a;
	REQUIRE(s.len() == 37);
	for (int i = 0; i < 37; i++) REQUIRE(s[i] == pa[i] + (i < 23 ? pb[i] : num(0)));
	for (int i = 0; i < 37; i++) REQUIRE(d[i] == (i < 23 ? pb[i] : num(0)) - pa[i]);
	// multiplying by x^k appends zeros to the storage; coefficients shift up
	poly<E> g = a;
	g.shift(2);
	REQUIRE(g.len() == 39);
	REQUIRE(g[0] == num(0));
	REQUIRE(g[1] == num(0));
	for (int i = 0; i < 37; i++) REQUIRE(g[i + 2] == pa[i]);
	REQUIRE(g.rev_series().data()[0] == pa[36]);
	// named conversions use the reversed convention and round-trip freely
	const power_series_exact<E>& ra = a.rev_series();
	REQUIRE(ra.len() == 37);
	for (int i = 0; i < 37; i++) REQUIRE(ra[i] == pa[36 - i]);
	REQUIRE(poly<E>::from_rev_series(ra) == a);
	// a poly's natural-order coefficients as an exact series (for series products)
	power_series_exact<E> xa(a.begin(), a.end());
	REQUIRE(equal(xa.begin(), xa.end(), pa.begin(), pa.end()));
	REQUIRE(a.unrev_series(10) == power_series<E>(pa.begin(), pa.begin() + 10));
	// the storage transform serves transposed products: middle product against rev_series()
	std::vector<num> vals(60);
	fill_rnd(vals, mt);
	cached_power_series_exact<E> cv(power_series_exact<E>(vals.begin(), vals.end()));
	cached_power_series_exact<E> ca(a.rev_series());
	auto mp = middle_product(cv, ca);
	auto naive = [&](int j) {
		num r{};
		for (int t = 0; t < 37; t++) r += pa[t] * vals[j + t];
		return r;
	};
	for (int j = 0; j < sz(mp); j++) REQUIRE(mp[size_t(j)] == naive(j));
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
		vector<num> coeffs(n);
		for (num& x : coeffs) { x = num(mt()); }
		poly<fft_engine<num>> p((span<const num>(coeffs)));
		vector<num> pts(n);
		for (int i = 0; i < n; i++) pts[i] = num(1000 + i);
		auto vals = poly_evaluate<fft_engine<num>>(p, pts);
		for (int i = 0; i < n; i++) {
			REQUIRE(vals[i] == p(pts[i]));
		}
		auto rec = poly_interpolate<fft_engine<num>>(pts, vals);
		REQUIRE(rec == p);
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

// marker engine: same ring, but declared non-commutative so the squarer
// must take the multiplier path (no cross-term doubling)
struct nc_engine : fft_engine<modnum<998244353>> {
	static constexpr bool commutative = false;
};
static_assert(fft::conv_engine<nc_engine>);

TEST_CASE("online squarer non-commutative fallback", "[fft]") {
	using num = modnum<998244353>;
	mt19937 mt(Catch::getSeed());
	for (int n : {1, 2, 3, 17, 64, 100}) {
		INFO("n = " << n);
		vector<num> f(n);
		fill_rnd(f, mt);
		auto slow = multiply_slow(f, f);
		slow.resize(2*n, num(0));
		online_squarer<nc_engine> os(n);
		for (int i = 0; i < n; i++) {
			os.push(f[i]);
			REQUIRE(os.back() == slow[i]);
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
