#include <catch2/catch_test_macros.hpp>

#include "tensor.hpp"

#include <string>

TEST_CASE("Tensor", "[tensor]") {
	using ten = tensor<std::string, 2>;
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
