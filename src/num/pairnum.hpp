#pragma once

#include <iostream>

#include "modnum.hpp"

namespace wala {

template <typename U, typename V> struct pairnum : num_ops<pairnum<U, V>> {
	using Self = pairnum;
	U u;
	V v;

	pairnum() : u(0), v(0) {}
	pairnum(long long val) : u(val), v(val) {}
	pairnum(const U& u_, const V& v_) : u(u_), v(v_) {}

	friend std::ostream& operator << (std::ostream& out, const Self& n) { return out << '(' << n.u << ',' << ' ' << n.v << ')'; }
	friend std::istream& operator >> (std::istream& in, Self& n) { long long val; in >> val; n = Self(val); return in; }

	friend bool operator == (const Self& a, const Self& b) { return a.u == b.u && a.v == b.v; }

	Self inv() const {
		return Self(u.inv(), v.inv());
	}
	Self neg() const {
		return Self(u.neg(), v.neg());
	}

	Self& operator ++ () {
		++u, ++v;
		return *this;
	}
	Self& operator -- () {
		--u, --v;
		return *this;
	}

	Self& operator += (const Self& o) {
		u += o.u;
		v += o.v;
		return *this;
	}
	Self& operator -= (const Self& o) {
		u -= o.u;
		v -= o.v;
		return *this;
	}
	Self& operator *= (const Self& o) {
		u *= o.u;
		v *= o.v;
		return *this;
	}
	Self& operator /= (const Self& o) {
		u /= o.u;
		v /= o.v;
		return *this;
	}
};

} // namespace wala
