#pragma once

#include <cassert>
#include <concepts>
#include <ranges>

#include "containers/vec.hpp"

namespace wala {

template <typename V> concept palindrome_text = std::ranges::random_access_range<V> && std::ranges::sized_range<V> && std::equality_comparable<std::ranges::range_reference_t<V>>;

/**
 * manacher(S): return the maximum palindromic substring of S centered at each point
 *
 * Input: string (or vector) of length N (no restrictions on character-set)
 * Output: vector res of length 2*N+1
 *   For any 0 <= i <= 2*N:
 *   * i % 2 == res[i] % 2
 *   * the half-open substring S[(i-res[i])/2, (i+res[i])/2) is a palindrome of length res[i]
 *   * For odd palindromes, take odd i, and vice versa
 */
template <palindrome_text V> vec<int> manacher(const V& S) {
	int N = int(std::ranges::size(S));
	auto s = std::ranges::begin(S);
	vec<int> res(2*N+1, 0);
	for (int i = 1, j = -1, r = 0; i < 2*N; i++, j--) {
		if (i > r) {
			r = i+1, res[i] = 1;
		} else {
			res[i] = res[j];
		}
		if (i+res[i] >= r) {
			int b = r>>1, a = i-b;
			while (a > 0 && b < N && s[a-1] == s[b]) {
				a--, b++;
			}
			res[i] = b-a, j = i, r = b<<1;
		}
	}
	return res;
}

/**
 * manacher_odd(S): return the maximum palindromic substring of S centered at each point
 *
 * Input: string (or vector) of length N (no restrictions on character-set)
 * Output: vector res of length N
 *   For any 0 <= i < N:
 *   * the half-open substring S[i-res[i], i+res[i]] is a palindrome of length 2*res[i]+1
 */
template <palindrome_text V> vec<int> manacher_odd(const V& S) {
	int N = int(std::ranges::size(S));
	auto s = std::ranges::begin(S);
	vec<int> res(N);
	for (int i = 1, j = -1, r = 0; i < N; i++, j--) {
		if (i > r) {
			r = i, res[i] = 0;
		} else {
			res[i] = res[j];
		}
		if (i+res[i] >= r) {
			int b = r, a = 2*i-r;
			while (a-1 >= 0 && b+1 < N && s[a-1] == s[b+1]) {
				a--, b++;
			}
			res[i] = b-i, j = i, r = b;
		}
	}
	return res;
}

} // namespace wala
