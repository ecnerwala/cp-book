#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_get_random_seed.hpp>

#include "fft/engines/algebras.hpp"
#include "fft/engines/ntt.hpp"
#include "fft/engines/split.hpp"
#include "fft/engines/crt.hpp"
#include "fft/multiply.hpp"
#include "fft/online.hpp"
#include "fft/test_util.test.hpp"
#include "modnum.hpp"

namespace ecnerwala {
namespace fft {

using namespace std;

// engine concept sanity checks for the wrapping algebras
static_assert(engine<engines::matrix<engines::ntt<modnum<998244353>>, 2>>);
static_assert(engine<engines::trunc<engines::ntt<modnum<998244353>>, 3>>);
// tracked inner engines work when the accumulated scale fits the budget (N <= 2)
static_assert(engine<engines::matrix<engines::split<modnum<int(1e9) + 7>>, 2>>);
static_assert(engine<engines::trunc<engines::crt<modnum<int(1e9) + 7>>, 2>>);
// the stable variants keep tracked inner engines sound at any N
static_assert(engine<engines::matrix_stable<engines::split<modnum<int(1e9) + 7>>, 3>>);
static_assert(engine<engines::trunc_stable<engines::crt<modnum<int(1e9) + 7>>, 3>>);

template <typename E, bool online, int N>
void test_matrix_engine(mt19937& mt) {
	using M = typename E::value_type;
	using num = std::remove_reference_t<decltype(std::declval<M&>()[{0, 0}])>;
	auto rnd_mat = [&]() {
		M m;
		for (int r = 0; r < N; r++)
			for (int c = 0; c < N; c++) m[{r, c}] = rnd_val<num>(mt);
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
	slow.resize(2 * n, M{});
	vector<M> got(2 * n - 1);
	square<E>(span<const M>(f), span<M>(got));
	REQUIRE(got == vector<M>(slow.begin(), slow.begin() + 2 * n - 1));
	if constexpr (online) {
		online_squarer<E> os(n);
		for (int i = 0; i < n; i++) {
			os.push(f[i]);
			REQUIRE(os.back() == slow[i]);
		}
	}
}

TEMPLATE_TEST_CASE("matrix engine", "[fft]",
	engines::ntt<modnum<998244353>>,
	engines::split<modnum<int(1e9) + 7>>,
	engines::crt<modnum<int(1e9) + 7>>) {
	using IE = TestType;
	// the tracked engines' scale budget admits N = 2 (entries are N-addend sums), and
	// the non-commutative online squarer accumulates two N-addend products per window
	// (scale 2N), exceeding it
	constexpr int N = IE::unit_scale == 0 ? 3 : 2;
	mt19937 mt(Catch::getSeed());
	test_matrix_engine<engines::matrix<IE, N>, IE::unit_scale == 0, N>(mt);
	// the stable variant works at any N (and its online squarer stays at scale 2)
	test_matrix_engine<engines::matrix_stable<IE, 3>, true, 3>(mt);
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
	engines::ntt<modnum<998244353>>,
	engines::split<modnum<int(1e9) + 7>>,
	engines::crt<modnum<int(1e9) + 7>>) {
	using IE = TestType;
	using num = typename IE::value_type;
	constexpr int N = IE::unit_scale == 0 ? 3 : 2;
	mt19937 mt(Catch::getSeed());
	test_trunc_series_engine<engines::trunc<IE, N>, num, N>(mt);
	test_trunc_series_engine<engines::trunc_stable<IE, 3>, num, 3>(mt);
}

}
} // namespace ecnerwala::fft
