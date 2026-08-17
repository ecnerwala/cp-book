#pragma once
#include<bits/stdc++.h>

template <typename num>
std::vector<num> BerlekampMassey(const std::vector<num>& s) {
	int n = int(s.size()), L = 0, m = 0;
	std::vector<num> C(n), B(n), T;
	C[0] = B[0] = 1;

	num b = 1;
	for(int i = 0; i < n; i++) { ++m;
		num d = s[i];
		for (int j = 1; j <= L; j++) d += C[j] * s[i - j];
		if (d == 0) continue;
		T = C; num coef = d / b;
		for (int j = m; j < n; j++) C[j] -= coef * B[j - m];
		if (2 * L > i) continue;
		L = i + 1 - L; B = T; b = d; m = 0;
	}

	C.resize(L + 1); C.erase(C.begin());
	for (auto& x : C) {
		x = -x;
	}
	return C;
}

template <typename num>
num linearRec(const std::vector<num>& S, const std::vector<num>& tr, int64_t k) {
	int n = int(tr.size());
	assert(S.size() >= tr.size());

	auto combine = [&](std::vector<num> a, std::vector<num> b, bool e = false) {
		// multiply a * b * x^e
		std::vector<num> res(int(a.size()) + int(b.size()));
		for (int i = 0; i < int(a.size()); i++) {
			for (int j = 0; j < int(b.size()); j++) {
				res[i + j + e] += a[i] * b[j];
			}
		}
		for (int i = int(res.size())-1; i >= n; --i) {
			for (int j = 0; j < n; j++) {
				res[i - 1 - j] += res[i] * tr[j];
			}
		}
		res.resize(n);
		return res;
	};

	std::vector<num> pol(n);
	if (n > 0) pol[0] = num(1);

	assert(k >= 0);
	for (int i = 64 - 1 - (k == 0 ? 64 : __builtin_clzll(k)); i >= 0; i--) {
		pol = combine(pol, pol, (k >> i) & 1);
	}

	num res = 0;
	for (int i = 0; i < n; i++) res += pol[i] * S[i];
	return res;
}
