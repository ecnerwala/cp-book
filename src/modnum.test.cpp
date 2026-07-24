#include "modnum.hpp"
#include <numeric> // Include for std::lcm and std::gcd
#include <random>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_get_random_seed.hpp>

TEST_CASE("mod_goldilocks test", "[mod_goldilocks]") {
    std::mt19937_64 mt(Catch::getSeed());
    for (int i = -100; i < 100; i++) {
        for (int j = -100; j < 100; j++) {
            mod_goldilocks a(i);
            mod_goldilocks b(j);
            CAPTURE(a, b);
            REQUIRE(uint64_t((__uint128_t(uint64_t(a)) + __uint128_t(uint64_t(b))) % mod_goldilocks::MOD) == uint64_t(a+b));
            REQUIRE(uint64_t((__uint128_t(uint64_t(a)) - __uint128_t(uint64_t(b)) + mod_goldilocks::MOD) % mod_goldilocks::MOD) == uint64_t(a-b));
            REQUIRE(uint64_t((__uint128_t(uint64_t(a)) * __uint128_t(uint64_t(b))) % mod_goldilocks::MOD) == uint64_t(a*b));
        }
    }
    for (int z = 0; z < 1000; z++) {
        mod_goldilocks a(mt());
        mod_goldilocks b(mt());
        CAPTURE(a, b);
        REQUIRE(uint64_t((__uint128_t(uint64_t(a)) + __uint128_t(uint64_t(b))) % mod_goldilocks::MOD) == uint64_t(a+b));
        REQUIRE(uint64_t((__uint128_t(uint64_t(a)) - __uint128_t(uint64_t(b)) + mod_goldilocks::MOD) % mod_goldilocks::MOD) == uint64_t(a-b));
        REQUIRE(uint64_t((__uint128_t(uint64_t(a)) * __uint128_t(uint64_t(b))) % mod_goldilocks::MOD) == uint64_t(a*b));
    }
}

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
