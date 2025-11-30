#pragma once

#include<bits/stdc++.h>

const double PI = acos(-1.);
const double TAU = 2 * PI;

template <typename T, typename AreaT=T, typename VolT=T> struct Point3D {
	using P = Point3D;

	T x, y, z;
	Point3D() : x(0), y(0), z(0) {}
	Point3D(T x_, T y_, T z_) : x(x_), y(y_), z(z_) {}

	template <typename U, typename V, typename W> explicit Point3D(const Point3D<U, V, W>& p) : x(T(p.x)), y(T(p.y)), z(T(p.z)) {}

	friend std::istream& operator >> (std::istream& i, P& p) { return i >> p.x >> p.y >> p.z; }
	friend std::ostream& operator << (std::ostream& o, const P& p) { return o << "(" << p.x << "," << p.y << "," << p.z << ")"; }

	friend bool operator == (const P& a, const P& b) { return a.x == b.x && a.y == b.y && a.z == b.z; }
	friend bool operator != (const P& a, const P& b) { return a.x != b.x || a.y != b.y || a.z != b.z; }

	P& operator += (const P& o) { x += o.x, y += o.y, z += o.z; return *this; }
	P& operator -= (const P& o) { x -= o.x, y -= o.y, z -= o.z; return *this; }
	friend P operator + (const P& a, const P& b) { return P(a) += b; }
	friend P operator - (const P& a, const P& b) { return P(a) -= b; }

	P& operator *= (const T& t) { x *= t, y *= t, z *= t; return *this; }
	P& operator /= (const T& t) { x /= t, y /= t, z /= t; return *this; }
	friend P operator * (const P& p, const T& t) { return P(p) *= t; }
	friend P operator * (const T& t, const P& p) { return P(p) *= t; }
	friend P operator / (const P& a, const T& t) { return P(a) /= t; }

	friend P operator + (const P& a) { return P(+a.x, +a.y, +a.z); }
	friend P operator - (const P& a) { return P(-a.x, -a.y, -a.z); }

	friend AreaT dot(const P& a, const P& b) { return AreaT(a.x) * AreaT(b.x) + AreaT(a.y) * AreaT(b.y) + AreaT(a.z) * AreaT(b.z); }
	friend AreaT norm(const P& a) { return dot(a,a); }
	// We're playing a little loose with this type, expliitly cast it if you need
	friend Point3D<AreaT, VolT> cross(const P& a, const P& b) { return Point3D<AreaT, VolT>(AreaT(a.y) * AreaT(b.z) - AreaT(a.z) * AreaT(b.y), AreaT(a.z) * AreaT(b.x) - AreaT(a.x) * AreaT(b.z), AreaT(a.x) * AreaT(b.y) - AreaT(a.y) * AreaT(b.x)); }

	friend T int_norm(const P& p) {
		return std::gcd(std::gcd(abs(p.x), abs(p.y)), abs(p.z));
	}
	friend P int_unit(const P& p) {
		T g = int_norm(p);
		return g ? p / g : p;
	}

	friend T abs(const P& a) { return std::sqrt(std::max(T(0), norm(a))); }
	friend P unit(const P& a) { return a / abs(a); }

	friend VolT vol(const P& a, const P& b, const P& c, const P& d) { return dot(cross(b-a, c-a), Point3D<AreaT, VolT>(d-a)); }

	friend bool lexLess(const P& a, const P& b) { return std::tie(a.x, a.y, a.z) < std::tie(b.x, b.y, b.z); }

	friend bool parallelSame(const P& a, const P& b) {
		assert(a != P());
		assert(b != P());
		return lexLess(a, P()) == lexLess(b, P());
	}
};
