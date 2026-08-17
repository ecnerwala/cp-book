#pragma once

#include <utility>
#include <array>
#include <tuple>
#include <iostream>

template<typename num = int>
struct hurwitz_quaternion {
	// we store the doubled quaternion
	num s,x,y,z;
	hurwitz_quaternion() : s(0), x(0), y(0), z(0) {}
	hurwitz_quaternion(num v) : s(2*v), x(0), y(0), z(0) {}
	hurwitz_quaternion(num s_, num x_, num y_, num z_) : s(2*s_), x(2*x_), y(2*y_), z(2*z_) {}
	struct doubled_coords_tag {};
	hurwitz_quaternion(doubled_coords_tag, num s_, num x_, num y_, num z_) : s(s_), x(x_), y(y_), z(z_) {
		assert((s & 1) == (x & 1) && (s & 1) == (y & 1) && (s & 1) == (z & 1));
	}
	friend std::ostream& operator << (std::ostream& o, const hurwitz_quaternion& q) {
		o << double(q.s)/2;
		{
			std::ios_base::fmtflags f(o.flags());
			o << std::showpos << double(q.x)/2 << "i" << double(q.y)/2 << "j" << double(q.z)/2 << "k";
			o.flags(f);
		}
		return o;
	}

	explicit operator bool() const {
		return s || x || y || z;
	}

	friend bool operator == (const hurwitz_quaternion& a, const hurwitz_quaternion& b) {
		return std::tie(a.s,a.x,a.y,a.z) == std::tie(b.s,b.x,b.y,b.z);
	}
	friend bool operator != (const hurwitz_quaternion& a, const hurwitz_quaternion& b) { return !(a == b); }

	num real_doubled() const {
		return s;
	}
	num real() const {
		assert(!(s & 1));
		return s >> 1;
	}
	std::array<num, 3> imag_doubled() const {
		return {x, y, z};
	}
	std::array<num, 3> imag() const {
		assert(!(s & 1));
		return {x>>1, y>>1, z>>1};
	}
	std::array<num, 4> coords_doubled() const {
		return {s, x, y, z};
	}
	std::array<num, 4> coords() const {
		assert(!(s & 1));
		return {s>>1, x>>1, y>>1, z>>1};
	}

	friend num norm(const hurwitz_quaternion& q) {
		return (q.s * q.s + q.x * q.x + q.y * q.y + q.z * q.z) >> 2;
	}
	friend hurwitz_quaternion conj(const hurwitz_quaternion& q) {
		return hurwitz_quaternion(doubled_coords_tag{}, q.s, -q.x, -q.y, -q.z);
	}

	friend hurwitz_quaternion operator + (const hurwitz_quaternion& q) {
		return hurwitz_quaternion(doubled_coords_tag{}, +q.s, +q.x, +q.y, +q.z);
	}
	friend hurwitz_quaternion operator - (const hurwitz_quaternion& q) {
		return hurwitz_quaternion(doubled_coords_tag{}, -q.s, -q.x, -q.y, -q.z);
	}

	hurwitz_quaternion& operator += (const hurwitz_quaternion& o) {
		s += o.s;
		x += o.x;
		y += o.y;
		z += o.z;
		return *this;
	}
	friend hurwitz_quaternion operator + (const hurwitz_quaternion& a, const hurwitz_quaternion& b) {
		return hurwitz_quaternion(doubled_coords_tag{}, a.s + b.s, a.x + b.x, a.y + b.y, a.z + b.z);
	}
	hurwitz_quaternion& operator -= (const hurwitz_quaternion& o) {
		s -= o.s;
		x -= o.x;
		y -= o.y;
		z -= o.z;
		return *this;
	}
	friend hurwitz_quaternion operator - (const hurwitz_quaternion& a, const hurwitz_quaternion& b) {
		return hurwitz_quaternion(doubled_coords_tag{}, a.s - b.s, a.x - b.x, a.y - b.y, a.z - b.z);
	}

	friend hurwitz_quaternion operator * (const num& a, const hurwitz_quaternion& q) {
		return hurwitz_quaternion(doubled_coords_tag{}, a*q.s, a*q.x, a*q.y, a*q.z);
	}
	friend hurwitz_quaternion operator * (const hurwitz_quaternion& q, const num& a) {
		return hurwitz_quaternion(doubled_coords_tag{}, q.s*a, q.x*a, q.y*a, q.z*a);
	}
	hurwitz_quaternion& operator *= (const num& a) {
		s *= a;
		x *= a;
		y *= a;
		z *= a;
		return *this;
	}

	friend hurwitz_quaternion operator * (const hurwitz_quaternion& a, const hurwitz_quaternion& b) {
		return hurwitz_quaternion(
			doubled_coords_tag{},
			(a.s * b.s - a.x * b.x - a.y * b.y - a.z * b.z) >> 1,
			(a.s * b.x + a.x * b.s + a.y * b.z - a.z * b.y) >> 1,
			(a.s * b.y + a.y * b.s + a.z * b.x - a.x * b.z) >> 1,
			(a.s * b.z + a.z * b.s + a.x * b.y - a.y * b.x) >> 1
		);
	}
	hurwitz_quaternion& operator *= (const hurwitz_quaternion& o) {
		return *this = *this * o;
	}

	struct div_t {
		hurwitz_quaternion quot, rem;
	};
	// a = b * quot + rem
	friend div_t right_div(const hurwitz_quaternion& a, const hurwitz_quaternion& b) {
		hurwitz_quaternion numer = conj(b) * a;
		num denom = norm(b);

		auto floor_div = [](num u, num v) -> num {
			if ((u^v) >= 0) {
				return u/v;
			} else {
				auto res = std::div(u, v);
				return res.quot - bool(res.rem);
			}
		};
		num s = floor_div(numer.s, denom);
		num x = floor_div(numer.x, denom);
		num y = floor_div(numer.y, denom);
		num z = floor_div(numer.z, denom);

		hurwitz_quaternion q_odd(doubled_coords_tag{}, s | 1, x | 1, y | 1, z | 1);
		hurwitz_quaternion r_odd = a - b * q_odd;
		hurwitz_quaternion q_even(doubled_coords_tag{}, (s+1)&~num(1), (x+1)&~num(1), (y+1)&~num(1), (z+1)&~num(1));
		hurwitz_quaternion r_even = a - b * q_even;
		div_t res = norm(r_odd) < norm(r_even) ? div_t{q_odd, r_odd} : div_t{q_even, r_even};
		assert(norm(res.rem) < norm(b));
		return res;
	}

	// a = ga', b = gb'
	friend hurwitz_quaternion right_gcd(hurwitz_quaternion a, hurwitz_quaternion b) {
		while (a) {
			b = right_div(b, a).rem;
			std::swap(a, b);
		}
		return b;
	}
};
