#pragma once

#include <array>
#include <cstddef>

namespace wala {

// Truncated polynomial mod x^N over num
template <typename num, int N> struct trunc_series {
	std::array<num, size_t(N)> a{};
	num& operator[](int i) { return a[size_t(i)]; }
	const num& operator[](int i) const { return a[size_t(i)]; }
	num* data() { return a.data(); }
	const num* data() const { return a.data(); }
	trunc_series& operator+=(const trunc_series& o) { for (int i = 0; i < N; i++) a[i] += o.a[i]; return *this; }
	friend trunc_series operator+(trunc_series x, const trunc_series& y) { x += y; return x; }
	trunc_series& operator-=(const trunc_series& o) { for (int i = 0; i < N; i++) a[i] -= o.a[i]; return *this; }
	friend trunc_series operator-(trunc_series x, const trunc_series& y) { x -= y; return x; }
	friend trunc_series operator*(const trunc_series& x, const trunc_series& y) {
		trunc_series r;
		for (int i = 0; i < N; i++) for (int j = 0; j < N - i; j++) r[i + j] += x[i] * y[j];
		return r;
	}
	trunc_series& operator*=(const trunc_series& o) { return *this = *this * o; }
	friend bool operator==(const trunc_series&, const trunc_series&) = default;
};

} // namespace wala
