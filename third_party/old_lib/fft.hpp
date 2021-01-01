#pragma once

#include <utility>

template <int L_, typename num> struct fft {
	static constexpr int L = L_;
	static constexpr int N = (1 << L);
	using num_t = num;

	using arr = num[N];
	const num root;
	const num invRoot = num(1) / root;
	const num invN = num(1) / num(N);
	num roots[L+1];
	num invRoots[L+1];
	fft(num root_) : root(root_) {
		roots[0] = root;
		invRoots[0] = invRoot;
		for (int i = 1; i <= L; i++) {
			roots[i] = roots[i-1] * roots[i-1];
			invRoots[i] = invRoots[i-1] * invRoots[i-1];
		}
		//assert(roots[L] == 1);
		//assert(invRoots[L] == 1);
	}

	void bitReversal(arr a) {
		for (int i = 1, j = 0; i < N-1; i++) {
			for (int k = N / 2; k > (j ^= k); k /= 2);
			if (j < i) {
				std::swap(a[j], a[i]);
			}
		}
	}

	void operator () (arr a, bool inv = false) {
		bitReversal(a);
		for (int l = 1, p = 0; l < N; l <<= 1, p++) {
			num w = inv ? invRoots[L - p - 1] : roots[L - p - 1];
			for (int k = 0; k < N; k += (2 * l)) {
				num v = 1;
				for (int i = k; i < k + l; i++, v *= w) {
					num x = a[i];
					num y = v * a[i+l];
					a[i] = x + y;
					a[i+l] = x - y;
				}
			}
		}
		if (inv) {
			for (int i = 0; i < N; i++) {
				a[i] *= invN;
			}
		}
	}
};
