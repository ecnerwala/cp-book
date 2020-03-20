#pragma once

#include <vector>
#include <cassert>

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
template <typename V> std::vector<int> manacher(const V& S) {
	int N = int(S.size());
	std::vector<int> res(2*N+1, 0);
	for (int i = 1, j = -1, r = 0; i < 2*N; i++, j--) {
		if (i > r) {
			r = i+1, res[i] = 1;
		} else {
			res[i] = res[j];
		}
		if (i+res[i] >= r) {
			int b = r>>1, a = i-b;
			while (a > 0 && b < N && S[a-1] == S[b]) {
				a--, b++;
			}
			res[i] = b-a, j = i, r = b<<1;
		}
	}
	return res;
}
