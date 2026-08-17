#pragma once

#include <cassert>
#include <utility>

// Computes (n on m) == 1 using the binary-gcd method
// m must be positive and odd, and n must be relatively prime
template <typename T> bool is_qr_jacobi(T n, T m) {
	bool r = true;
	assert(m & 1);
	assert(m > 0);
	if (n < 0) {
		if (m & 2) r = !r;
		n = -n;
	}
	while (m > 1) {
		assert(n > 0);
		int t = __builtin_ctzll(n);
		n >>= t;
		if ((t & 1) && (((m & 7) == 3) || ((m & 7) == 5))) {
			r = !r;
		}
		// n and m both odd
		if (n < m) {
			if ((n & 2) && (m & 2)) {
				r = !r;
			}
			using std::swap;
			swap(n, m);
		}
		n -= m;
	}
	return r;
}
