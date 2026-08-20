#include <bits/stdc++.h>
#include <catch2/catch_test_macros.hpp>

#include "num/linear_fn.hpp"
#include "num/modnum.hpp"

namespace wala {

using namespace std;

TEST_CASE("linear_fn composition", "[num]") {
	using num = modnum<998244353>;
	using fn = linear_fn<num>;

	fn f{num(2), num(3)}, g{num(5), num(7)};
	num x(11);
	REQUIRE(f(x) == num(2) * x + num(3));
	REQUIRE((f * g)(x) == f(g(x)));
	REQUIRE((g * f)(x) == g(f(x)));

	fn id;
	REQUIRE(id(x) == x);
	REQUIRE(f * id == f);
	REQUIRE(id * f == f);

	fn h = f;
	h *= g;
	REQUIRE(h == f * g);
}

} // namespace wala
