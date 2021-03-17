#pragma once

#include <vector>
#include <bitset>
#include <cassert>

// Compute the characteristic polynomial of a square matrix A over some field.
// Not numerically stable at all.
// Takes argument by value, use std::move if you can.
template <typename num> std::vector<num> charPoly(std::vector<std::vector<num>> A) {
	int N = int(A.size());
	std::vector<num> res; res.reserve(N+1);
	res.push_back(num(1));
	for (int i = 0, deg = 0; i < N; i++) {
		auto& Ai = A[i];

		int c = i+1;
		while (c < N && Ai[c] == num(0)) c++;
		if (c == N) {
			res.resize(i+2, num(0));
			for (int x = deg; x >= 0; x--) {
				num v = res[x];
				for (int y = x+1, z = i; z >= deg; z--, y++) {
					res[y] -= v * Ai[z];
				}
			}
			deg = i+1;
			continue;
		}

		num vc = Ai[c];
		num ivc = inv(vc);

		Ai[c] = Ai[i+1];
		Ai[i+1] = 0;

		std::swap(A[i+1], A[c]);
		auto& Ai1 = A[i+1];
		for (int k = deg; k < N; k++) {
			Ai1[k] *= vc;
		}

		for (int k = i+1; k < N; k++) {
			auto& Ak = A[k];
			{
				auto& x = Ak[i+1];
				auto& y = Ak[c];
				num tmp = y;
				y = x;
				x = tmp * ivc;
			}
			{
				num v = Ak[i+1];
				for (int j = deg; j < N; j++) {
					Ak[j] -= v * Ai[j];
				}
			}
			if (k > i+1) {
				num v = Ai[k];
				for (int j = deg; j < N; j++) {
					Ai1[j] += v * Ak[j];
				}
			}
		}

		for (int k = deg; k <= i; k++) {
			Ai1[k+1] += Ai[k];
		}
	}
	reverse(res.begin(), res.end());
	return res;
}

// Compute the characteristic polynomial of a square matrix A over F2.
// Takes argument by value, use std::move if you can.
// Note that MAXS must be at least N+1
template <std::size_t MAXS> std::bitset<MAXS> charPoly(std::vector<std::bitset<MAXS>> A) {
	using bs = std::bitset<MAXS>;
	int N = int(A.size());
	assert(MAXS >= N+1);
	bs ans; ans[0] = 1;
	int deg = 0;
	for (int i = 0; i < N; i++) {
		{
			int j = int(A[i]._Find_next(i));
			if (j >= N) {
				bs nans;
				for (; deg <= i; ans <<= 1, deg++) {
					if (A[i][deg]) nans ^= ans;
				}
				ans ^= nans;
				continue;
			}
			if (j != i+1) {
				swap(A[j], A[i+1]);
				for (auto& a : A) {
					bool tmp = a[j];
					a[j] = a[i+1];
					a[i+1] = tmp;
				}
			}
		}
		assert(A[i][i+1]);
		bs msk = A[i]; msk.flip(i+1);
		for (int k = 0; k < N; k++) {
			if (msk[k]) A[i+1] ^= A[k];
		}
		for (auto& a : A) {
			if (a[i+1]) a ^= msk;
		}
	}
	return ans;
}
