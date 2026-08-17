#pragma once

#include <iostream>
#include <numeric>

template <typename T, typename MulT=T> struct fraction_t {
	T numer = 0, denom = 1;

	fraction_t() : numer(0), denom(1) {}
	fraction_t(T v) : numer(v), denom(1) {}
	fraction_t(T n, T d) : numer(n), denom(d) {
		if (denom < 0 || (denom == 0 && numer < 0)) {
			numer = -numer;
			denom = -denom;
		}
	}
	template <typename U, typename V> explicit fraction_t(const fraction_t<U, V> o) : numer(T(o.numer)), denom(T(o.denom)) {}

	friend std::ostream& operator << (std::ostream& o, const fraction_t& f) {
		return o << f.numer << '/' << f.denom;
	}
	friend std::istream& operator >> (std::istream& i, const fraction_t& f) {
		return i >> f.numer >> f.denom;
	}

	friend MulT cross(const fraction_t& a, const fraction_t& b) {
		return MulT(a.numer) * MulT(b.denom) - MulT(b.numer) * MulT(a.denom);
	}

	friend bool operator == (const fraction_t& a, const fraction_t& b) {
		return cross(a, b) == 0;
	}
	friend std::strong_ordering operator <=> (const fraction_t& a, const fraction_t& b) {
		return cross(a, b) <=> 0;
	}

	fraction_t operator + () const { return fraction_t(+numer, denom); }
	fraction_t operator - () const { return fraction_t(-numer, denom); }

	fraction_t& operator *= (const fraction_t& o) {
		numer *= o.numer;
		denom *= o.denom;
		return *this;
	}
	fraction_t& operator /= (const fraction_t& o) {
		numer *= o.denom;
		denom *= o.numer;
		return *this;
	}
	friend fraction_t operator * (const fraction_t& a, const fraction_t& b) {
		return fraction_t(a.numer * b.numer, a.denom * b.denom);
	}
	friend fraction_t operator / (const fraction_t& a, const fraction_t& b) {
		return fraction_t(a.numer * b.denom, a.denom * b.numer);
	}

	friend fraction_t operator + (const fraction_t& a, const fraction_t& b) {
		return {a.numer * b.denom + b.numer * a.denom, a.denom * b.denom};
	}
	friend fraction_t operator - (const fraction_t& a, const fraction_t& b) {
		return {a.numer * b.denom - b.numer * a.denom, a.denom * b.denom};
	}
	fraction_t& operator += (const fraction_t& o) { return *this = *this + o; }
	fraction_t& operator -= (const fraction_t& o) { return *this = *this - o; }

	fraction_t& reduce() {
		using std::gcd;
		T g = gcd(numer, denom);
		numer /= g;
		denom /= g;
		return *this;
	}
};
