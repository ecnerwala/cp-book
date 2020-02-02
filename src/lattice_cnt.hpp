#pragma once

#include <utility>
#include <cassert>

// number of integer solutions to Ax + By <= C and x,y >= 0
long long lattice_cnt(long long A, long long B, long long C) {
	using ll = long long;

	assert(A >= 0 && B >= 0);
	if (C < 0) return 0;

	assert(A > 0 && B > 0);
	if (A > B) std::swap(A, B);
	assert(A <= B);

	ll ans = 0;
	while (C >= 0) {
		assert(0 < A && A <= B);

		ll k = B/A;
		ll l = B%A;
		assert(B == k * A + l);

		ll f = C/B;
		ll e = C%B / A;
		ll g = C%B % A;
		assert(C == f * B + e * A + g);
		assert(C == (f * k + e) * A + f * l + g);

		// either x + ky <= f*k+e
		// i.e. 0 <= x <= (f-y) * k + e
		// or x >= fk + e + 1 - ky
		// and Ax + (Ak+l) y <= C = (fk + e + 1) A + fl - A + g
		// Let z = x - (fk + e + 1 - ky)
		// Az + A(fk + e + 1 - ky) + Aky + ly <= C = A (fk + e + 1) + fl - A + g
		// Az + ly <= fl - A + g

		ans += (f+1) * (e+1) + (f+1) * f / 2 * k;

		C = f*l - A + g;
		B = A;
		A = l;
	}
	return ans;
}

// count the number of 0 <= (a * x % m) < c for 0 <= x < n
long long mod_count(long long a, long long m, long long c, long long n) {
	assert(a >= 0 && m > 0 && n >= 0 && c >= 0);
	if (c >= m) return n;
	assert(c < m);

	if (n == 0) return 0;

	a %= m;

	// we want solutions to 0 <= (a+m)x - my < c with 0 <= x <= N-1
	// iff 0 <= (a+m)x - my < c with 0 <= x and y <= (a+m) * (n-1) / m
	return lattice_cnt(m, a+m, (a+m) * (n-1) / m * m + c - 1) - lattice_cnt(m, a+m, (a+m) * (n-1) / m * m - 1);
}
