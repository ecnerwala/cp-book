#pragma once

#include <array>
#include <cstddef>

namespace wala {

// Small NxN matrix over num, row-major
template <typename num, int N> struct mat {
	std::array<num, size_t(N) * N> a{};
	num& operator[](std::array<int, 2> rc) { return a[size_t(rc[0]) * N + rc[1]]; }
	const num& operator[](std::array<int, 2> rc) const { return a[size_t(rc[0]) * N + rc[1]]; }
	num* data() { return a.data(); }
	const num* data() const { return a.data(); }
	mat& operator+=(const mat& o) { for (int i = 0; i < N*N; i++) a[i] += o.a[i]; return *this; }
	friend mat operator+(mat x, const mat& y) { x += y; return x; }
	mat& operator-=(const mat& o) { for (int i = 0; i < N*N; i++) a[i] -= o.a[i]; return *this; }
	friend mat operator-(mat x, const mat& y) { x -= y; return x; }
	friend mat operator*(const mat& x, const mat& y) {
		mat r;
		for (int i = 0; i < N; i++) for (int k = 0; k < N; k++) for (int j = 0; j < N; j++)
			r[{i, j}] += x[{i, k}] * y[{k, j}];
		return r;
	}
	mat& operator*=(const mat& o) { return *this = *this * o; }
	friend bool operator==(const mat&, const mat&) = default;
};

} // namespace wala
