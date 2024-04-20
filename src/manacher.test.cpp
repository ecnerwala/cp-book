#define CATCH_CONFIG_MAIN
#include "manacher.hpp"
#include <catch2/catch.hpp>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <iostream>

// Naive implementation for checking correctness of manacher algorithm
std::vector<int> naive_manacher(const std::string& s) {
    int n = s.size();
    std::vector<int> res(2 * n + 1, 0);
    for (int i = 0; i < 2 * n + 1; ++i) {
        if (i == 0) {
            res[i] = 0; // Index 0 should always be 0 as there's no palindrome before the start of the string
            continue;
        }
        int len = (i % 2 == 0) ? 0 : 1; // Start len at 0 for even indices, 1 for odd
        while (i / 2 - len >= 0 && i / 2 + len < n) {
            if (i % 2 == 0) {
                if (i / 2 - len - 1 < 0 || s[i / 2 - len - 1] != s[i / 2 + len]) break;
            } else {
                if (s[i / 2 - len] != s[i / 2 + len]) break;
            }
            ++len;
        }
        res[i] = 2 * len - (i % 2); // Store the length of the palindrome centered at index i
    }
    return res;
}

// Function to test a single test case of manacher algorithm
void check_stuff(const std::string& test_case) {
    // Call the manacher function
    auto result = manacher(test_case);
    // Call the naive implementation
    auto expected = naive_manacher(test_case);

    // Check the postconditions of the output
    REQUIRE(result.size() == 2 * test_case.size() + 1);

    // Check that the reported palindromes are correct and match the naive implementation
    for (size_t i = 0; i < result.size(); ++i) {
        if (result[i] != expected[i]) {
            std::cout << "Mismatch found on test case: '" << test_case << "' at index " << i << std::endl;
            std::cout << "manacher() result: " << result[i] << ", naive_manacher() expected: " << expected[i] << std::endl;
        }
        REQUIRE(result[i] == expected[i]);
    }
}

// Define the test case for the manacher algorithm
TEST_CASE("Manacher's Algorithm", "[manacher]") {
    // Test with small cases
    std::vector<std::string> small_cases = {"a", "aa", "aba", "abc", "aaa"};
    for (const auto& test_case : small_cases) {
        check_stuff(test_case);
    }

    // Test with larger, randomized cases
    std::mt19937 rng(48); // Standard mersenne_twister_engine seeded with 48
    for (int len : {8, 13, 21, 34, 55, 89}) { // Fibonacci numbers as lengths
        for (int i = 0; i < 5; ++i) { // 5 random test cases for each size
            std::string test_case;
            for (int j = 0; j < len; ++j) {
                char c = std::uniform_int_distribution<char>('a', 'z')(rng);
                test_case += c;
            }
            check_stuff(test_case);
        }
    }
}
