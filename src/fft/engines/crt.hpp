#pragma once

#include <cassert>
#include <cstdint>
#include <span>
#include <utility>

#include "fft/engine.hpp"
#include "fft/engines/ntt.hpp"
#include "num/modnum.hpp"

namespace ecnerwala::fft::engines {

// Multiplies mod `mnum` by running NTTs modulo two FFT-friendly primes and CRT'ing.
// Inputs use balanced representatives (|v| <= MOD/2), so the true integer coefficients
// are bounded by n (MOD/2)^2.
template <typename mnum, typename num1 = mod_goldilocks, typename num2 = modnum<(15 << 27) + 1>>
struct crt {
	static_assert(sizeof(decltype(mnum::MOD)) <= 4, "n (MOD/2)^2 must fit the CRT modulus product");
	using value_type = mnum;
	static constexpr bool commutative = true;
	static constexpr int unit_scale = 1;
	using E1 = ntt<num1>;
	using E2 = ntt<num2>;
	template <int A = 1> struct transformed_t {
		typename E1::transformed t1;
		typename E2::transformed t2;
		int size() const { return t1.size(); }
		transformed_t() = default;
		transformed_t(typename E1::transformed&& t1_, typename E2::transformed&& t2_)
			: t1(std::move(t1_)), t2(std::move(t2_)) {}
		template <int A2> requires (A2 != A) explicit(A2 > A) transformed_t(transformed_t<A2>&& o)
			: t1(std::move(o.t1)), t2(std::move(o.t2)) {}
	};
	using transformed = transformed_t<1>;
	template <int K> struct product_t {
		typename E1::product p1;
		typename E2::product p2;
		int size() const { return sz(p1); }
		product_t() = default;
		product_t(typename E1::product&& p1_, typename E2::product&& p2_)
			: p1(std::move(p1_)), p2(std::move(p2_)) {}
		template <int K2> requires (K2 != K) explicit(K2 > K) product_t(product_t<K2>&& o)
			: p1(std::move(o.p1)), p2(std::move(o.p2)) {}
	};
	using product = product_t<1>;

	static transformed transform(std::span<const mnum> a, int n) {
		assert(sz(a) <= 2 * n);
		auto b1 = buffer_pool<num1>::get(sz(a));
		auto b2 = buffer_pool<num2>::get(sz(a));
		for (int i = 0; i < sz(a); i++) { int64_t v = a[i].balanced(); b1[i] = num1(v); b2[i] = num2(v); }
		return transformed{
			E1::transform(std::span<const num1>(b1.span()), n),
			E2::transform(std::span<const num2>(b2.span()), n),
		};
	}
	static void extend_to(transformed& t, int m, std::span<const mnum> coeffs) {
		if (t.size() >= m) return;
		auto b1 = buffer_pool<num1>::get(sz(coeffs));
		auto b2 = buffer_pool<num2>::get(sz(coeffs));
		for (int i = 0; i < sz(coeffs); i++) { int64_t v = coeffs[i].balanced(); b1[i] = num1(v); b2[i] = num2(v); }
		E1::extend_to(t.t1, m, std::span<const num1>(b1.span()));
		E2::extend_to(t.t2, m, std::span<const num2>(b2.span()));
	}
	template <int A> static transformed_t<A> downsample(const transformed_t<A>& t, int n, bool odd) {
		return transformed_t<A>{E1::downsample(t.t1, n, odd), E2::downsample(t.t2, n, odd)};
	}
	template <int K> static product_t<K> downsample(const product_t<K>& p, int n, bool odd) {
		return product_t<K>{E1::downsample(p.p1, n, odd), E2::downsample(p.p2, n, odd)};
	}
	template <int A> static transformed_t<A> upsample(const transformed_t<A>& t, int n, bool odd) {
		return transformed_t<A>{E1::upsample(t.t1, n, odd), E2::upsample(t.t2, n, odd)};
	}
	template <int K> static product_t<K> upsample(const product_t<K>& p, int n, bool odd) {
		return product_t<K>{E1::upsample(p.p1, n, odd), E2::upsample(p.p2, n, odd)};
	}
	template <int A> static transformed_t<A> negate_arg(const transformed_t<A>& t, int n) {
		return transformed_t<A>{E1::negate_arg(t.t1, n), E2::negate_arg(t.t2, n)};
	}
	// Exact per prime; the scale tracks the true (integer) coefficient growth.
	template <int A, int B> static transformed_t<A + B> add(transformed_t<A>&& a, const transformed_t<B>& b) {
		return transformed_t<A + B>{E1::add(std::move(a.t1), b.t1), E2::add(std::move(a.t2), b.t2)};
	}
	template <int A, int B> static product_t<A * B> mul(const transformed_t<A>& a, const transformed_t<B>& b, int n) {
		return product_t<A * B>{E1::mul(a.t1, b.t1, n), E2::mul(a.t2, b.t2, n)};
	}
	template <int A> static product_t<A * A> sq(const transformed_t<A>& a, int n) { return mul(a, a, n); }
	template <int A1, int B1, int A2, int B2>
	static product_t<A1 * B1 + A2 * B2> mul2(
		const transformed_t<A1>& a1, const transformed_t<B1>& b1,
		const transformed_t<A2>& a2, const transformed_t<B2>& b2,
		int n
	) {
		return product_t<A1 * B1 + A2 * B2>{
			E1::mul2(a1.t1, b1.t1, a2.t1, b2.t1, n),
			E2::mul2(a1.t2, b1.t2, a2.t2, b2.t2, n),
		};
	}
	template <int K1, int K2> static product_t<K1 + K2> add(product_t<K1>&& a, product_t<K2>&& b) {
		return product_t<K1 + K2>{E1::add(std::move(a.p1), b.p1), E2::add(std::move(a.p2), b.p2)};
	}
	template <int K = 1, typename Op = assign_op> static void finish(product_t<K>&& p, std::span<mnum> out, Op op = {}) {
		// The reconstruction needs |c| < whole/2; balanced inputs bound each addend's
		// true coefficients by n (MOD/2)^2, so the safe length is divided by the
		// accumulated scale. K <= 2 is very conservative (~2^35 even for MOD ~ 2^30).
		static_assert(K <= 2, "crt: accumulated scale too large");
		int n = p.size();
		assert(sz(out) <= n);
		auto o1 = buffer_pool<num1>::get(sz(out));
		auto o2 = buffer_pool<num2>::get(sz(out));
		E1::finish(std::move(p.p1), o1.span());
		E2::finish(std::move(p.p2), o2.span());

		// TODO: Could hardcode these
		num1 inv_n2 = inv(num1(num2::MOD));
		num2 inv_n1 = inv(num2(num1::MOD));
		__int128_t whole = __int128_t(num1::MOD) * __int128_t(num2::MOD);

		mnum m1_mod = mnum(num1::MOD);
		mnum m2_mod = mnum(num2::MOD);
		mnum whole_mod = m1_mod * m2_mod;
		for (int i = 0; i < sz(out); i++) {
			num1 v1 = o1[i] * inv_n2;
			num2 v2 = o2[i] * inv_n1;
			mnum o_mod = mnum(uint64_t(v1)) * m2_mod + mnum(int(v2)) * m1_mod;
			__int128_t o_exact = __int128_t(uint64_t(v1)) * __int128_t(num2::MOD) + __int128_t(int(v2)) * __int128_t(num1::MOD);
			if (o_exact >= whole) { o_exact -= whole; o_mod -= whole_mod; }
			// Balanced representatives: |o| <= whole/2
			if (o_exact > whole / 2) o_mod -= whole_mod;
			op(out[i], o_mod);
		}
	}
};

/* namespace ecnerwala::fft::engines */ }
