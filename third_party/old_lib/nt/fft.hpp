#pragma once

#include <cassert>

namespace pe { namespace fft {

using ll = long long;

template <typename T> void bitReversal(T *a, int n) {
	// i and j are reverses
	for (int j = 1, i = 0; j < n-1; j++) {
		// do j++ except with i
		for (int k = n/2; k > (i^=k); k /= 2);
		if (j < i) swap(a[i], a[j]);
	}
}

constexpr int powmod(int a, int b, int m) {
	a %= m;
	int r = 1;
	while (b) {
		if (b & 1) r = int((ll(r)*a) % m);
		a = int((ll(a)*a) % m);
		b /= 2;
	}
	return r;
}

constexpr int inv(int a, int m) {
	return a == 1 ? 1 : m - int(ll(inv(m % a, a)) * m / a);
}

template <int n_param, int p_param, int rt_param> struct NTT {
	static constexpr int n = n_param;
	static constexpr int p = p_param;
	static constexpr int rt = rt_param;
	static constexpr int rtinv = inv(rt, p); // 1/w

	static_assert((n & (n-1)) == 0);
	static_assert(ll(rt) * rtinv % p == 1);
	static_assert(powmod(rt, n, p) == 1);
	static_assert(powmod(rt, n/2, p) == p-1);

	static void fft(int a[n], int sz, bool invert=false) {
		assert((sz & (sz - 1)) == 0);
		assert(sz <= n);
		for (int i = 0; i < sz; i++) {
			a[i] %= p;
			if (a[i] < 0) a[i] += p;
		}
		int w = invert ? rt : rtinv;
		for (int h = n; h > 1; h /= 2, w = int(ll(w) * w % p)) {
			if (h > sz) continue; // multiply w all the way up
			assert(sz % h == 0);
			for (int k = 0; k < sz; k += h) {
				int wi = 1;
				for (int i = 0; i < h/2; i++, wi = int(ll(wi) * w % p)) {
					int l = a[k+i], r = a[k+i+h/2];
					a[k+i] = int((ll(l)+r) % p);
					a[k+i+h/2] = int((ll(l)+ll(p-r)) * wi % p);
				}
				assert(wi == p-1);
			}
		}
		assert(w == 1);

		bitReversal(a, sz);

		if (invert) {
			int sinv = inv(sz, p);
			for (int i = 0; i < sz; i++) a[i] = int(ll(sinv) * a[i] % p);
		}
	}

	static void ifft(int a[n], int sz) {
		fft(a, sz, true);
	}

	static void test() {
		int *a = new int[n];
		// test 1: just 1
		for (int i = 0; i < n; i++) {
			a[i] = (i == 0);
		}
		fft(a, n);
		assert(a[1] == 1);
		assert(a[2] == 1);
		assert(a[3] == 1);
		assert(a[4] == 1);
		assert(a[1234567] == 1);
		ifft(a, n);
		for (int i = 0; i < n; i++) {
			assert(a[i] == (i == 0));
		}
		// test 2: lots of numbers
		for (int i = 0; i < n; i++) {
			a[i] = i;
		}
		fft(a, n);
		ifft(a, n);
		for (int i = 0; i < n; i++) {
			assert(a[i] == i);
		}

		fft(a, 32);
		assert(a[32] == 32); // hasn't been touched
		ifft(a, 32);
		for (int i = 0; i < n; i++) {
			assert(a[i] == i);
		}
		delete[] a;
	}

	// these things can be fully equal, but cannot overlap
	static void polymul(int dst[n], const int a[n], const int b[n], int sz) {
		if (dst == a && dst == b) {
			fft(dst, sz);
			for (int i = 0; i < sz; i++) {
				dst[i] = int(ll(dst[i]) * dst[i] % p);
			}
			ifft(dst, sz);
			return;
		} else {
			assert(dst != a || dst != b);
			if (dst == b) {
				swap(a, b);
			}
			assert(dst != b);
			int *tmp = new int[n];
			memcpy(tmp, b, sizeof(int) * sz);
			fft(tmp, sz);
			if (dst != a) memcpy(dst, a, sizeof(int) * sz);
			fft(dst, sz);
			for (int i = 0; i < sz; i++) {
				dst[i] = int(ll(dst[i]) * tmp[i] % p);
			}
			ifft(dst, sz);
			delete[] tmp;
		}
	}
};

using ntt1 = NTT<1ll << 25, (1 << 25) * 33 + 1, 309>;
using ntt2 = NTT<1ll << 25, (1 << 25) * 51 + 1, 40>;
using ntt3 = NTT<1ll << 25, (1 << 25) * 54 + 1, 103>;

// MODULUS and number of terms to keep
const int MOD = int(1e9) + 7;
const int N = ll(1e7);

int tmp1[1 << 25];
int tmp2[1 << 25];
int tmp3[1 << 25];
void polymul(int dst[1 << 25], const int src1[1 << 25], const int src2[1 << 25], int lim = N) {
	int sz = 1 << 25;
	while (sz > 1 && sz / 2 > 2 * lim) sz /= 2;
	assert(sz > 2 * lim);
	memcpy(dst, src1, sizeof(int) * size_t(sz));
	memcpy(tmp3, src2, sizeof(int) * size_t(sz));
	for (int i = 0; i <= lim; i++) {
		dst[i] %= MOD;
		if (dst[i] < 0) dst[i] += MOD;
		tmp3[i] %= MOD;
		if (tmp3[i] < 0) tmp3[i] += MOD;
	}
	for (int i = lim+1; i < sz; i++) {
		dst[i] = tmp3[i] = 0;
	}
	ntt1::polymul(tmp1, dst, tmp3, sz);
	ntt2::polymul(tmp2, dst, tmp3, sz);
	ntt3::polymul(tmp3, dst, tmp3, sz);
	ll m = ntt1::p, n = ntt2::p, o = ntt3::p;
	ll invmn = inv(int(m), int(n));
	ll invmno = inv(int(m*n%o), int(o));
	for (int i = 0; i <= lim; i++) {
		ll a = tmp1[i], b = tmp2[i], c = tmp3[i];
		assert(0 <= a && a < m);
		assert(0 <= b && b < n);
		assert(0 <= c && c < o);

		ll v = a + ((b-a)%n*invmn) % n * m;
		if (v < 0) v += m*n;
		assert(0 <= v && v < m*n);
		// result is v + ((c-v) * invmno) % o * mn;
		// we should have set it up so that MOD * MOD * n < o * m * n so that the maximum precise input fits in our crt
		ll ofac = (c%o - v%o) % o * invmno % o;
		if (ofac < 0) ofac += o;
		assert(0 <= ofac && ofac < o);
		ll res = v % MOD + (m % MOD) % MOD * (n % MOD) % MOD * (ofac % MOD) % MOD;
		res %= MOD;
		if (res < 0) res += MOD;
		dst[i] = int(res);
	}
	for (int i = lim+1; i < sz; i++) {
		dst[i] = 0;
	}
}

void test_polymul() {
	int *d = new int[1 << 25];
	int *s = new int[1 << 25];
	int *t = new int[1 << 25];
	for (int i = N+1; i < (1 << 25); i++) {
		s[i] = t[i] = 0;
	}

	// test 1: maximum product
	for (int i = 0; i <= N; i++) {
		s[i] = t[i] = MOD-1;
	}
	polymul(d, s, t);
	for (int i = 0; i <= N; i++) {
		assert(s[i] == MOD-1);
		assert(t[i] == MOD-1);
		assert(d[i] == i+1);
	}

	// Test 2: negatives
	for (int i = 0; i <= N; i++) {
		s[i] = MOD-1, t[i] = 1;
	}
	polymul(d, s, t);
	for (int i = 0; i <= N; i++) {
		assert(s[i] == MOD-1);
		assert(t[i] == 1);
		assert(d[i] == MOD-(i+1));
	}
	delete[] d;
	delete[] s;
	delete[] t;
}

void polyinv(int dst[1 << 25], int src[1 << 25]) {
	assert(src[0]);
	dst[0] = inv(src[0], MOD);
	int *tmp = new int[1 << 25];
	for (int l = 1; l < N; l *= 2) {
		polymul(tmp, src, dst, 2*l-1); // first find our current product
		assert(tmp[0] == 1);
		for (int i = 1; i < l; i++) {
			assert(tmp[i] == 0);
		}
		// transform tmp to the 1-remainder
		for (int i = 0; i < l; i++) {
			tmp[i] = -tmp[i+l];
		}
		for (int i = l; i < (1 << 25); i++) {
			tmp[i] = 0;
		}
		polymul(tmp, tmp, dst, l-1); // multiply the remainder by dst to compute the inverse
		for (int i = 0; i < l && i+l <= N; i++) {
			dst[i+l] = tmp[i];
		}
		cerr << "iteration l = " << l << " done" << '\n';
	}
	polymul(tmp, src, dst);
	assert(tmp[0] == 1);
	for (int i = 1; i <= N; i++) {
		assert(tmp[i] == 0);
	}
	delete[] tmp;
}

}} // namespace pe::fft
