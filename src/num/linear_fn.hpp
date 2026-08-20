#pragma once

namespace wala {

// Linear function x -> a * x + b over num.
// operator* is composition: (f * g)(x) = f(g(x)).
// Default-constructs to the identity.
template <typename num> struct linear_fn {
	num a = 1;
	num b = 0;
	num operator()(const num& x) const { return a * x + b; }
	friend linear_fn operator*(const linear_fn& f, const linear_fn& g) {
		return {f.a * g.a, f.a * g.b + f.b};
	}
	linear_fn& operator*=(const linear_fn& o) { return *this = *this * o; }
	friend bool operator==(const linear_fn&, const linear_fn&) = default;
};

} // namespace wala
