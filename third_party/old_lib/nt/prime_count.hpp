#pragma once

#include <bits/stdc++.h>

namespace nt { namespace prime_count {

using namespace std;

typedef long long ll;

// reference implementation of a counter
struct counter {
	ll v;

	counter() : v(0) {}
	counter(ll) : v(1) {}
	static counter sum_all(ll x) {
		assert(x >= 0);
		counter res;
		res.v = x;
		return res;
	}

	counter& operator += (const counter &other) {
		v += other.v;
		return *this;
	}
	counter& operator -= (const counter &other) {
		v -= other.v;
		return *this;
	}
	counter& operator *= (const counter &other) {
		v *= other.v;
		return *this;
	}

	counter operator + (const counter &other) const {
		counter res = *this;
		res += other;
		return res;
	}
	counter operator - (const counter &other) const {
		counter res = *this;
		res -= other;
		return res;
	}
	counter operator * (const counter &other) const {
		counter res = *this;
		res *= other;
		return res;
	}
};
inline std::ostream& operator << (std::ostream &out, const counter &cnt) {
	return out << cnt.v;
}

// Compute pi(x) = \sum_{prime p<=x} T(p) for all x = \floor(N/a)
// T must implement all ring methods
// As a constructor function, T must be fully multiplicative and support efficiently finding things like sum_{i=1}^n T(i)
// T(x) = 1, T(x) = x^c supports this, as do some more complicated functions like T(x) = pair(x%4==1, x%4==3)
template <typename T, ll N> struct prime_count {

	static constexpr ll cbrt(ll a) {
		assert(a >= 0);
		ll r = 0;
		while (r * r * r <= a) r ++;
		r--;
		return r;
	}

	static constexpr ll A = cbrt(N);
	static constexpr ll B = N / A;

	ll min_prime[B+1];
	T pi_small[B+1]; // precompute the phis

	ll P;
	ll primes[B+1];

	ll PA;

	ll to_inv(ll x) {
		ll n = N / x;
		assert(x == N / n);
		assert(n <= A);
		return n;
	}

	T phi[A+1]; // phi[i] = phi(N/i, PA-1)
	// where phi(x, p) = sum of values <= x with all prime factors > primes[p]

	/* Equivalent, slow calculations of phi
	T phi_slow(ll x, ll i) {
		if (x <= 0) return T();
		if (i < 0) {
			return T::sum_all(x);
		}
		T res = phi_slow(x, i-1) - T(primes[i]) * phi_slow(x / primes[i], i-1);
		return res;
	}

	T phi_slow2(ll x, ll i) {
		ll p = (i >= 0) ? primes[i] : 0;
		T res;
		for (ll v = 1; v <= x; v++) {
			if (min_prime[v] == 0 || min_prime[v] > p) {
				res += T(v);
			}
		}
		return res;
	}
	*/

	void compute_phi() {
		// binary-indexed tree
		T (*dp)[A+1] = new T[A+1][A+1]; // dp[j][p] = sum T(i) for i <= N/j and min_prime[i] > p

		T *bit = new T[A+1];
		for (int i = 0, j = B; j > A; j--) {
			assert(N/j <= B);
			assert(i <= N/j);
			while (i < N/j) {
				i++;
				// insert i into the bit for all values <= min_prime[i]
				T val(i);
				for (ll v = (min_prime[i] == 0) ? A : min(min_prime[i], A); v; v -= v & (-v)) {
					bit[v] += val;
				}
			}
			assert(i == N/j);

			// now we iterate over prime factors of j
			for (ll cur = j; cur > 1; ) {
				ll p = min_prime[cur];
				assert(j % p == 0);
				assert(j > A);
				if (p > A) break;
				if (j / p <= A) {
					// number of values with min_prime >= p
					T phi_cur; // phi(N/j, p)
					for (ll v = p; v <= A; v += v & (-v)) {
						phi_cur += bit[v];
					}
					dp[j/p][p] -= T(p) * phi_cur;
				}
				while (cur % p == 0) cur /= p;
			}
		}
		delete[] bit;

		cerr << "second sieve done\n";

		for (int j = 1; j <= A; j++) {
			dp[j][0] = T::sum_all(N/j);
		}

		for (int i = 0; i < PA; i++) {
			ll q = (i > 0) ? primes[i-1] : 0;
			ll p = primes[i];
			for (int j = 1; j <= A; j++) {
				dp[j][p] += dp[j][q];
				if (j * p <= A) {
					dp[j][p] -= T(p) * dp[j*p][q];
				} else {
					// our first pass covered this
				}
			}
		}

		for (int j = 1; j <= A; j++) {
			phi[j] = dp[j][primes[PA-1]];
		}

		delete[] dp;
	}

	T p2_A(ll x) {
		static unordered_map<ll, T> mp;
		if (mp.count(x)) return mp[x];

		assert(x <= N);
		T res;
		for (ll i = PA; i < P && primes[i] * primes[i] <= x; i++) {
			ll p = primes[i];
			assert(x/p >= p);
			res += T(p) * (pi(x/p) - pi(p-1));
		}

		return mp[x] = res;
	}

	T pi(ll x) {
		if (x <= B) {
			return pi_small[x];
		}
		T res = pi(A);
		//res += phi_slow(x, PA-1);
		res += phi[to_inv(x)];
		res -= T(1);
		res -= p2_A(x);
		return res;
	}

	T count(ll x) {
		return pi(x);
	}

	prime_count() {
		for (ll i = 0; i <= B; i++) {
			if (i <= 1) min_prime[i] = 0;
			else if (i % 2 == 0) min_prime[i] = 2;
			else min_prime[i] = i;
		}

		for (ll p = 3; p * p <= B; p++) {
			if (min_prime[p] != p) continue;
			for (ll j = p * p; j <= B; j += p) {
				if (min_prime[j] == j) min_prime[j] = p;
			}
		}
		cerr << "sieve done\n";

		pi_small[0] = T();
		P = 0;
		for (ll p = 1; p <= B; p++) {
			pi_small[p] = pi_small[p-1];

			if (min_prime[p] != p) continue;
			primes[P++] = p;
			pi_small[p] += T(p);
		}

		PA = 0;
		while (PA < P && primes[PA] <= A) PA++;
		assert(PA <= P);

		cerr << "primes gathered\n";

		compute_phi();

		cerr << "precompute done\n";
	}
};

}} // namespace nt::prime_count
