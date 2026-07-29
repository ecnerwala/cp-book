#pragma once

#include <cassert>
#include <concepts>
#include <cstdint>
#include <iostream>
#include <limits>
#include <type_traits>
#include <utility>

template <typename T> T mod_inv_in_range(T a, T m) {
	// assert(0 <= a && a < m);
	T x = a, y = m;
	// abs coeff of a in x and y (they're always opposite sign)
	T vx = 1, vy = 0;
	bool swap = false;
	while (x) {
		T k = y / x;
		y %= x;
		vy += k * vx;
		std::swap(x, y);
		std::swap(vx, vy);
		swap ^= 1;
	}
	assert(y == 1);
	return swap ? vy : m - vy;
}

template <typename T> struct extended_gcd_result {
	T gcd;
	T coeff_a, coeff_b;
};
template <typename T> extended_gcd_result<T> extended_gcd(T a, T b) {
	T x = a, y = b;
	// coeff of a and b in x and y
	T ax = 1, ay = 0;
	T bx = 0, by = 1;
	while (x) {
		T k = y / x;
		y %= x;
		ay -= k * ax;
		by -= k * bx;
		std::swap(x, y);
		std::swap(ax, ay);
		std::swap(bx, by);
	}
	return {y, ay, by};
}

template <typename T> T mod_inv(T a, T m) {
	a %= m;
	a = a < 0 ? a + m : a;
	return mod_inv_in_range(a, m);
}

// Derives the boilerplate operator surface of a number type from its compound
// ops, ==, neg(), and inv().
// Bodies are only instantiated on use, so a type may omit some of the
// underlying pieces if the corresponding derived ops are never called.
template <typename Self>
struct num_ops {
	Self operator+ () const { return static_cast<const Self&>(*this); }
	Self operator- () const { return static_cast<const Self&>(*this).neg(); }

	friend Self operator ++ (Self& a, int) { Self r = a; ++a; return r; }
	friend Self operator -- (Self& a, int) { Self r = a; --a; return r; }
	friend Self operator + (const Self& a, const Self& b) { return Self(a) += b; }
	friend Self operator - (const Self& a, const Self& b) { return Self(a) -= b; }
	friend Self operator * (const Self& a, const Self& b) { return Self(a) *= b; }
	friend Self operator / (const Self& a, const Self& b) { return Self(a) /= b; }

	friend bool operator != (const Self& a, const Self& b) { return !(a == b); }

	friend Self neg(const Self& a) { return a.neg(); }
	friend Self inv(const Self& a) { return a.inv(); }
};

// Storage and arithmetic for numbers mod Self::MOD, as a reduced
// representative v in [0, MOD) of unsigned type V.
// The type provides static MOD (of type V), reduce (value -> representative),
// and *=;
// everything else is derived here, valid for any MOD up to V's full range
// (sums and differences are tracked mod 2^bits, so no headroom is needed).
// Hooks may be overridden in the type's own body (e.g. a faster += / -=).
template <typename Self, typename V>
struct mod_ops : num_ops<Self> {
	static_assert(std::unsigned_integral<V>);
	V v;

	struct is_reduced_tag {};

	mod_ops() : v(0) {}
	mod_ops(V v_, is_reduced_tag) : v(v_) { assert(v < Self::MOD); }
	template <std::integral I> mod_ops(I x) : v(Self::reduce(x)) {}

	static Self from_reduced(V v) { return Self(v, is_reduced_tag{}); }

	// A negative value reduces via its nonnegative complement: x = -1 - ~x.
	static V reduce(std::signed_integral auto x) {
		using U = std::make_unsigned_t<decltype(x)>;
		return x < 0 ? V(Self::MOD - 1 - Self::reduce(U(~x))) : Self::reduce(U(x));
	}

	explicit operator V() const { return v; }
	std::make_signed_t<V> balanced() const {
		return std::make_signed_t<V>(Self::MOD-v > v ? v : v - Self::MOD);
	}

	friend bool operator == (const Self& a, const Self& b) { return a.v == b.v; }
	friend std::ostream& operator << (std::ostream& out, const Self& n) { return out << n.v; }
	friend std::istream& operator >> (std::istream& in, Self& n) { int64_t v_; in >> v_; n = Self(v_); return in; }

	Self& operator ++ () {
		++v;
		if (v == Self::MOD) v = 0;
		return self();
	}
	Self& operator -- () {
		if (v == 0) v = Self::MOD;
		--v;
		return self();
	}
	Self& operator += (const Self& o) { v = Self::sub_mod_raw(v, Self::MOD - o.v); return self(); }
	Self& operator -= (const Self& o) { v = Self::sub_mod_raw(v, o.v); return self(); }
	Self& operator /= (const Self& o) { return self() *= o.inv(); }

	// Returns a - b mod MOD, for b in [0, MOD]; wraparound detects the underflow.
	static V sub_mod_raw(V a, V b) { return a < b ? a - b + Self::MOD : a - b; }

	Self neg() const { return from_reduced(v ? Self::MOD - v : 0); }
	Self inv() const { return from_reduced(mod_inv_in_range(v, Self::MOD)); }

private:
	Self& self() { return static_cast<Self&>(*this); }
};

template <auto MOD_> struct modnum : mod_ops<modnum<MOD_>, std::make_unsigned_t<decltype(MOD_)>> {
	using Self = modnum;
	static_assert(MOD_ > 0, "MOD must be positive");
	using V = std::make_unsigned_t<decltype(MOD_)>;
	static constexpr V MOD = V(MOD_);

	using base = mod_ops<modnum, V>;
	using base::base;
	using base::v;
	using base::reduce;

	static V reduce(std::unsigned_integral auto x) { return V(x % MOD); }

	explicit operator std::make_signed_t<V>() const
		requires (MOD <= V(std::numeric_limits<std::make_signed_t<V>>::max()))
	{
		return std::make_signed_t<V>(v);
	}

	Self& operator *= (const Self& o) {
		if constexpr (sizeof(V) <= 4) v = V(uint64_t(v) * o.v % MOD);
		else v = V(__uint128_t(v) * o.v % MOD);
		return *this;
	}
};

struct mod_goldilocks : mod_ops<mod_goldilocks, uint64_t> {
	using Self = mod_goldilocks;
	static constexpr uint64_t MOD = 0xffffffff00000001ull;
	static constexpr uint64_t EPS = -MOD;
	// We have 2^32 is a primitive 6th root of unity.
	// Note that omega_8 + omega_8^7 == 2^24 - 2^72 == sqrt(2)
	// We'll pick the root so that 2^24 - 2^72 is our primitive 384th root of unity.
	static constexpr uint64_t PRIMITIVE_ROOT = 2717;

	using base = mod_ops<mod_goldilocks, uint64_t>;
	using base::base;
	using base::reduce;
	mod_goldilocks() = default;
	mod_goldilocks(__int128_t a) : base(a < 0 ? uint64_t(MOD - 1 - __uint128_t(~a) % MOD) : uint64_t(__uint128_t(a) % MOD), is_reduced_tag{}) {}
	mod_goldilocks(__uint128_t a) : base(uint64_t(a % MOD), is_reduced_tag{}) {}

	// Avoids the division: any uint64_t is within MOD of reduced.
	static uint64_t reduce(std::unsigned_integral auto x) {
		static_assert(sizeof(x) <= 8);
		uint64_t a = x;
		return a >= MOD ? a - MOD : a;
	}

	// returns a-b, assuming -MOD <= a-b, e.g. b <= MOD
	static uint64_t sub_mod_raw(uint64_t a, uint64_t b) {
#if defined(__x86_64__)
		// TODO: We could try to write this using intrinsics, but GCC sometimes produces the wrong code.
		uint64_t res_wrapped = a;
		uint64_t adjustment = b;
		asm (
			// AT&T syntax: SRC DST
			"sub %[y], %[x]\n\t"
			// Trick from plonky2 implementation:
			// After sub, flag CF is set iff we underflowed. We want to correct by EPS == 2^32 - 1 iff C is set.
			// sbb (subtract with borrow) computes DST <- DST - SRC - CF
			// Thus, we can use the 32-bit form of sbb on a dummy register to load CF ? EPS : 0.
			// Here, we'll just reuse the original register holding b.
			"sbb %k[y], %k[y]\n\t"
			: [x] "+r"(res_wrapped),
			[y] "+r"(adjustment)
			:
			: "cc"
		);
#else
		uint64_t res_wrapped = a - b;
		uint64_t adjustment = (res_wrapped > a) ? EPS : 0;
#endif
		return res_wrapped - adjustment;
	}

	// Reduce lo + 2^64 * mi + 2^96 * hi, where hi <= MOD
	static uint64_t reduce_u160_raw(uint64_t lo, uint32_t mi, uint64_t hi) {
		// result = lo - hi + EPS * mi
		// 0 <= lo <= 2^64 - 1 = MOD + EPS - 1
		// 0 <= EPS * mi <= (2^32 - 1) * EPS = MOD - 1 - EPS
		// 0 <= hi <= MOD
		// -MOD <= lo - hi + EPS * mi <= 2*MOD-2
		// so we do have some leeway
		return sub_mod_raw(sub_mod_raw(lo, hi), MOD-(uint64_t(mi)<<32)+mi);
	}

	static uint64_t reduce_u128_raw(__uint128_t v) {
		uint64_t hi = uint64_t(v >> 64);
		uint64_t lo = uint64_t(v);
		uint32_t hi_hi = uint32_t(hi >> 32);
		uint32_t hi_lo = uint32_t(hi);
		return reduce_u160_raw(lo, hi_lo, hi_hi);
	}

	Self& operator *= (Self o) {
		v = reduce_u128_raw(__uint128_t(v) * __uint128_t(o.v));
		return *this;
	}
};

template <typename T> T power(T a, long long b) {
	assert(b >= 0);
	T r = 1; while (b) { if (b & 1) r *= a; b >>= 1; a *= a; } return r;
}

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

template <typename tag> struct dynamic_modnum : mod_ops<dynamic_modnum<tag>, uint32_t> {
	using Self = dynamic_modnum;
private:
	inline static uint32_t MOD_ = 0;
	inline static uint64_t BARRETT_M = 0;

public:
	// Make only the const-reference public, to force the use of set_mod
	static constexpr uint32_t const& MOD = MOD_;

	using base = mod_ops<dynamic_modnum, uint32_t>;
	using base::base;
	using base::v;
	using base::reduce;

	// Barret reduction taken from KACTL:
	/**
	 * Author: Simon Lindholm
	 * Date: 2020-05-30
	 * License: CC0
	 * Source: https://en.wikipedia.org/wiki/Barrett_reduction
	 * Description: Compute $a \% b$ about 5 times faster than usual, where $b$ is constant but not known at compile time.
	 * Returns a value congruent to $a \pmod b$ in the range $[0, 2b)$.
	 * Status: proven correct, stress-tested
	 * Measured as having 4 times lower latency, and 8 times higher throughput, see stress-test.
	 * Details:
	 * More precisely, it can be proven that the result equals 0 only if $a = 0$,
	 * and otherwise lies in $[1, (1 + a/2^64) * b)$.
	 */
	static void set_mod(int mod) {
		assert(mod > 0);
		MOD_ = uint32_t(mod);
		BARRETT_M = (uint64_t(-1) / MOD);
	}
	static uint32_t barrett_reduce_partial(uint64_t a) {
		return uint32_t(a - uint64_t((__uint128_t(BARRETT_M) * a) >> 64) * MOD);
	}
	static uint32_t barrett_reduce(uint64_t a) {
		int32_t res = int32_t(barrett_reduce_partial(a) - MOD);
		return uint32_t((res < 0) ? res + int32_t(MOD) : res);
	}

	struct mod_reader {
		friend std::istream& operator >> (std::istream& i, mod_reader) {
			int mod; i >> mod;
			Self::set_mod(mod);
			return i;
		}
	};
	static mod_reader MOD_READER() {
		return mod_reader();
	}

	static uint32_t reduce(std::unsigned_integral auto x) {
		static_assert(sizeof(x) <= 8);
		return barrett_reduce(x);
	}

	explicit operator int() const { return int(v); }

	Self& operator *= (const Self& o) {
		v = barrett_reduce(uint64_t(v) * o.v);
		return *this;
	}
};

template <typename T> struct mod_constraint {
	T v, mod;

	friend mod_constraint operator & (mod_constraint a, mod_constraint b) {
		if (a.mod < b.mod) std::swap(a, b);
		if (b.mod == 1) return a;

		extended_gcd_result<T> egcd = extended_gcd<T>(a.mod, b.mod);
		assert(a.v % egcd.gcd == b.v % egcd.gcd);

		T extra = b.v - a.v % b.mod;
		extra /= egcd.gcd;

		extra *= egcd.coeff_a;
		extra %= b.mod / egcd.gcd;
		extra += (extra < 0) ? b.mod / egcd.gcd : 0;

		return mod_constraint{
			a.v + extra * a.mod,
			a.mod * (b.mod / egcd.gcd)
		};
	}
};
