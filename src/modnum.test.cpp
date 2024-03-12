#include "modnum.hpp"
#include <catch2/catch_test_macros.hpp>
#include <numeric> // Include for std::lcm and std::gcd

TEST_CASE("Mod Constraint Regression Test", "[mod_constraint]") {
    for (int a_mod = 1; a_mod <= 10; ++a_mod) {
        for (int a_val = 0; a_val < a_mod; ++a_val) {
            for (int b_mod = 1; b_mod <= 10; ++b_mod) {
                for (int b_val = 0; b_val < b_mod; ++b_val) {
                    if (a_val % std::gcd(a_mod, b_mod) != b_val % std::gcd(a_mod, b_mod)) continue;

                    mod_constraint<int> a{a_val, a_mod};
                    mod_constraint<int> b{b_val, b_mod};

                    mod_constraint<int> r = a & b;

                    // Check that r.mod is the LCM of a.mod and b.mod
                    int lcm_ab = std::lcm(a.mod, b.mod);
                    REQUIRE(r.mod == lcm_ab);

                    // Check that r.v % a.mod == a.v (and likewise for b)
                    REQUIRE(r.v % a.mod == a.v);
                    REQUIRE(r.v % b.mod == b.v);

                    // Check that r.v is between 0 and r.mod
                    REQUIRE(r.v >= 0);
                    REQUIRE(r.v < r.mod);
                }
            }
        }
    }
}
