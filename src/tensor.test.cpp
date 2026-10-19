#include <catch2/catch_test_macros.hpp>

#include "tensor.hpp"

#include <string>

TEST_CASE("Tensor", "[wala::tensor]") {
	using ten = wala::tensor<std::string, 2>;
	ten a({2, 3});
	a[{0,0}] = "0";
	a[{0,1}] = "1";
	a[{0,2}] = "2";
	a[{1,0}] = "3";
	a[{1,1}] = "4";
	a[{1,2}] = "5";

	const ten const_a = a;
	ten b = a;
	REQUIRE(b[{0,0}] == "0");
	REQUIRE(b[{0,1}] == "1");
	REQUIRE(b[{0,2}] == "2");
	REQUIRE(b[{1,0}] == "3");
	REQUIRE(b[{1,1}] == "4");
	REQUIRE(b[{1,2}] == "5");

	// Bounds checked
	REQUIRE(b.at({0,0}) == "0");
	REQUIRE(b.at({0,1}) == "1");
	REQUIRE(b.at({0,2}) == "2");
	REQUIRE(b.at({1,0}) == "3");
	REQUIRE(b.at({1,1}) == "4");
	REQUIRE(b.at({1,2}) == "5");

	REQUIRE(*b[0][0] == "0");
	REQUIRE(*b[0][1] == "1");
	REQUIRE(*b[0][2] == "2");
	REQUIRE(*b[1][0] == "3");
	REQUIRE(*b[1][1] == "4");
	REQUIRE(*b[1][2] == "5");

	REQUIRE(*const_a[0][0] == "0");
	REQUIRE(*const_a[0][1] == "1");
	REQUIRE(*const_a[0][2] == "2");
	REQUIRE(*const_a[1][0] == "3");
	REQUIRE(*const_a[1][1] == "4");
	REQUIRE(*const_a[1][2] == "5");
}

TEST_CASE("Tensor copy and move", "[wala::tensor]") {
	using ten = wala::tensor<std::string, 2>;
	ten a({2, 3}, "x");
	a[{1, 2}] = "y";

	ten b = a;
	REQUIRE(b[{1, 2}] == "y");
	b[{0, 0}] = "z";
	REQUIRE(a[{0, 0}] == "x");

	ten c = std::move(a);
	REQUIRE(c[{0, 0}] == "x");
	REQUIRE(c[{1, 2}] == "y");

	ten d;
	d = std::move(c);
	REQUIRE(d[{1, 2}] == "y");
	REQUIRE(d[{0, 1}] == "x");

	d = b;
	REQUIRE(d[{0, 0}] == "z");
	REQUIRE(d[{1, 2}] == "y");
	d[{0, 0}] = "w";
	REQUIRE(b[{0, 0}] == "z");

	swap(b, d);
	REQUIRE(b[{0, 0}] == "w");
	REQUIRE(d[{0, 0}] == "z");

	// Moved-from is a valid tensor and can be reassigned
	a = ten({1, 1}, "q");
	REQUIRE(a[{0, 0}] == "q");

	wala::tensor<std::unique_ptr<int>, 1> p({3});
	*p[1] = std::make_unique<int>(5);
	wala::tensor<std::unique_ptr<int>, 1> q = std::move(p);
	REQUIRE(**q[1] == 5);
	REQUIRE(*q[0] == nullptr);
	wala::tensor<std::unique_ptr<int>, 1> r;
	r = std::move(q);
	REQUIRE(**r[1] == 5);
}
