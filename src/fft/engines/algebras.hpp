#pragma once

#include <array>
#include <concepts>
#include <cstddef>
#include <span>
#include <utility>

#include "fft/engine.hpp"
#include "num/mat.hpp"
#include "num/trunc_series.hpp"

namespace wala::fft::engines {

// ==== wrapper engines ====

// componentwise

// These are "componentwise engines" which model free modules/algebras over the underlying ring.
// Matrices are the canonical example: we can take entry-wise transforms, then multiply/add in transformed space.
//
// We start with a shared `componentwise` base class which handles all linear ops, i.e. not mul.
//
// If the underlying has unit_scale = 1, we may need to avoid accumulation of sums in transformed space;
// then, the product and the transformed data may have different dimensions.
// We'll represent this by an array Ofs of prefix offsets mapping each input/transform-space dimension to a range of product-space dimensions.
// Specifically, out[c] = sum prod[Ofs[c]:Ofs[c+1]]

template <int L> constexpr std::array<int, size_t(L) + 1> componentwise_iota = [] {
	std::array<int, size_t(L) + 1> r{};
	for (int i = 0; i <= L; i++) r[size_t(i)] = i;
	return r;
}();

template <engine E, typename V, int L, std::array<int, size_t(L) + 1> Ofs = componentwise_iota<L>>
struct componentwise {
	using S = typename E::value_type;
	using value_type = V;
	static constexpr int P = Ofs[size_t(L)];  // total product components
	static constexpr int unit_scale = E::unit_scale;
	template <int A = unit_scale> struct transformed_t {
		std::array<typename E::template transformed_t<A>, size_t(L)> t;
		int size() const { return t[0].size(); }
		transformed_t() = default;
		template <int A2> requires (A2 != A) explicit(A2 > A) transformed_t(transformed_t<A2>&& o) {
			for (int c = 0; c < L; c++)
				t[c] = typename E::template transformed_t<A>(std::move(o.t[c]));
		}
	};
	using transformed = transformed_t<>;
	// TODO: if E::product_t == E::transformed_t, mirror that here
	template <int K> struct product_t {
		std::array<typename E::template product_t<K>, size_t(P)> t;
		int size() const { return t[0].size(); }
		product_t() = default;
		template <int K2> requires (K2 != K) explicit(K2 > K) product_t(product_t<K2>&& o) {
			for (int c = 0; c < P; c++)
				t[c] = typename E::template product_t<K>(std::move(o.t[c]));
		}
	};

	static transformed transform(std::span<const V> a, int n) {
		transformed r;
		auto buf = buffer_pool<S>::get(sz(a));
		for (int c = 0; c < L; c++) {
			for (int i = 0; i < sz(a); i++) buf[i] = a[i].data()[c];
			r.t[c] = E::transform(std::span<const S>(buf.span()), n);
		}
		return r;
	}
	static void extend_to(transformed& t, int m, std::span<const V> coeffs) {
		if (t.size() >= m) return;
		auto buf = buffer_pool<S>::get(sz(coeffs));
		for (int c = 0; c < L; c++) {
			for (int i = 0; i < sz(coeffs); i++) buf[i] = coeffs[i].data()[c];
			E::extend_to(t.t[c], m, std::span<const S>(buf.span()));
		}
	}
	template <int A> static transformed_t<A> downsample(const transformed_t<A>& t, int n, bool odd) {
		transformed_t<A> r;
		for (int c = 0; c < L; c++) r.t[c] = E::downsample(t.t[c], n, odd);
		return r;
	}
	template <int K> static product_t<K> downsample(const product_t<K>& p, int n, bool odd) {
		product_t<K> r;
		for (int c = 0; c < P; c++) r.t[c] = E::downsample(p.t[c], n, odd);
		return r;
	}
	template <int A> static transformed_t<A> upsample(const transformed_t<A>& t, int n, bool odd) {
		transformed_t<A> r;
		for (int c = 0; c < L; c++) r.t[c] = E::upsample(t.t[c], n, odd);
		return r;
	}
	template <int K> static product_t<K> upsample(const product_t<K>& p, int n, bool odd) {
		product_t<K> r;
		for (int c = 0; c < P; c++) r.t[c] = E::upsample(p.t[c], n, odd);
		return r;
	}
	template <int A> static transformed_t<A> negate_arg(const transformed_t<A>& t, int n) {
		transformed_t<A> r;
		for (int c = 0; c < L; c++) r.t[c] = E::negate_arg(t.t[c], n);
		return r;
	}
	template <int A, int B> static transformed_t<A + B> add(transformed_t<A>&& a, const transformed_t<B>& b) {
		transformed_t<A + B> r;
		for (int c = 0; c < L; c++) r.t[c] = E::add(std::move(a.t[c]), b.t[c]);
		return r;
	}
	template <int K1, int K2> static product_t<K1 + K2> add(product_t<K1>&& a, product_t<K2>&& b) {
		product_t<K1 + K2> r;
		for (int c = 0; c < P; c++) r.t[c] = E::add(std::move(a.t[c]), std::move(b.t[c]));
		return r;
	}
	template <int K, typename Op = assign_op> static void finish(product_t<K>&& p, std::span<V> out, Op op = {}) {
		auto buf = buffer_pool<S>::get(sz(out));
		auto emit = [&](std::span<V> dst) {
			for (int c = 0; c < L; c++) {
				E::finish(std::move(p.t[Ofs[size_t(c)]]), buf.span());
				for (int j = Ofs[size_t(c)] + 1; j < Ofs[size_t(c) + 1]; j++)
					E::finish(std::move(p.t[j]), buf.span(), add_op{});
				for (int i = 0; i < sz(dst); i++) dst[i].data()[c] = buf[i];
			}
		};
		// Op must see each out element whole, exactly once, so compose
		// non-assign ops through an element buffer.
		if constexpr (std::same_as<Op, assign_op>) {
			emit(out);
		} else {
			auto vbuf = buffer_pool<V>::get(sz(out));
			emit(vbuf.span());
			for (int i = 0; i < sz(out); i++) op(out[i], vbuf.span()[i]);
		}
	}
};

// Convolve mat<N> (NxN matrices), with accumulation in product space
template <engine E, int N>
struct matrix : componentwise<E, mat<typename E::value_type, N>, N * N> {
	using base = componentwise<E, mat<typename E::value_type, N>, N * N>;
	static constexpr bool commutative = false;
	static constexpr int unit_scale = base::unit_scale;
	template <int A = unit_scale> using transformed_t = typename base::template transformed_t<A>;
	template <int K> using product_t = typename base::template product_t<K * N>;
	using transformed = typename base::transformed;
	using product = product_t<unit_scale * unit_scale>;

	// right fold over k so a tracked inner engine's per-addend types line up
	template <int A, int B, int k = 0>
	static auto entry(const transformed_t<A>& a, const transformed_t<B>& b, int r, int c, int n) {
		auto e = E::mul(a.t[size_t(r) * N + k], b.t[size_t(k) * N + c], n);
		if constexpr (k + 1 == N) return e;
		else return E::add(std::move(e), entry<A, B, k + 1>(a, b, r, c, n));
	}
	template <int A, int B>
	static product_t<A * B> mul(const transformed_t<A>& a, const transformed_t<B>& b, int n) {
		product_t<A * B> p;
		for (int r = 0; r < N; r++) for (int c = 0; c < N; c++)
			p.t[size_t(r) * N + c] = entry<A, B>(a, b, r, c, n);
		return p;
	}
	template <int A> static auto sq(const transformed_t<A>& a, int n) { return mul(a, a, n); }
	template <int A1, int B1, int A2, int B2, int k = 0>
	static auto entry2(
		const transformed_t<A1>& a1, const transformed_t<B1>& b1,
		const transformed_t<A2>& a2, const transformed_t<B2>& b2,
		int r, int c, int n
	) {
		auto e = E::mul2(
			a1.t[size_t(r) * N + k], b1.t[size_t(k) * N + c],
			a2.t[size_t(r) * N + k], b2.t[size_t(k) * N + c],
			n
		);
		if constexpr (k + 1 == N) return e;
		else return E::add(std::move(e), entry2<A1, B1, A2, B2, k + 1>(a1, b1, a2, b2, r, c, n));
	}
	template <int A1, int B1, int A2, int B2>
	static product_t<A1 * B1 + A2 * B2> mul2(
		const transformed_t<A1>& a1, const transformed_t<B1>& b1,
		const transformed_t<A2>& a2, const transformed_t<B2>& b2,
		int n
	) {
		product_t<A1 * B1 + A2 * B2> p;
		for (int r = 0; r < N; r++) for (int c = 0; c < N; c++)
			p.t[size_t(r) * N + c] = entry2<A1, B1, A2, B2>(a1, b1, a2, b2, r, c, n);
		return p;
	}
};

// Convolve trunc_series<num, N> (power series truncated at N), with accumulation in product space
template <engine E, int N>
struct trunc : componentwise<E, trunc_series<typename E::value_type, N>, N> {
	using base = componentwise<E, trunc_series<typename E::value_type, N>, N>;
	static constexpr bool commutative = E::commutative;
	static constexpr int unit_scale = base::unit_scale;
	template <int A = unit_scale> using transformed_t = typename base::template transformed_t<A>;
	template <int K> using product_t = typename base::template product_t<K * N>;
	using transformed = typename base::transformed;
	using product = product_t<unit_scale * unit_scale>;

	template <int A, int B, int s, int i = 0>
	static auto entry(const transformed_t<A>& a, const transformed_t<B>& b, int n) {
		auto e = E::mul(a.t[size_t(i)], b.t[size_t(s - i)], n);
		if constexpr (i == s) return e;
		else return E::add(std::move(e), entry<A, B, s, i + 1>(a, b, n));
	}
	template <int A, int B>
	static product_t<A * B> mul(const transformed_t<A>& a, const transformed_t<B>& b, int n) {
		product_t<A * B> p;
		[&]<size_t... s_>(std::index_sequence<s_...>) {
			((p.t[s_] = entry<A, B, int(s_)>(a, b, n)), ...);
		}(std::make_index_sequence<size_t(N)>{});
		return p;
	}
	template <int A> static auto sq(const transformed_t<A>& a, int n) { return mul(a, a, n); }
	template <int A1, int B1, int A2, int B2, int s, int i = 0>
	static auto entry2(
		const transformed_t<A1>& a1, const transformed_t<B1>& b1,
		const transformed_t<A2>& a2, const transformed_t<B2>& b2,
		int n
	) {
		auto e = E::mul2(a1.t[size_t(i)], b1.t[size_t(s - i)], a2.t[size_t(i)], b2.t[size_t(s - i)], n);
		if constexpr (i == s) return e;
		else return E::add(std::move(e), entry2<A1, B1, A2, B2, s, i + 1>(a1, b1, a2, b2, n));
	}
	template <int A1, int B1, int A2, int B2>
	static product_t<A1 * B1 + A2 * B2> mul2(
		const transformed_t<A1>& a1, const transformed_t<B1>& b1,
		const transformed_t<A2>& a2, const transformed_t<B2>& b2,
		int n
	) {
		product_t<A1 * B1 + A2 * B2> p;
		[&]<size_t... s_>(std::index_sequence<s_...>) {
			((p.t[s_] = entry2<A1, B1, A2, B2, int(s_)>(a1, b1, a2, b2, n)), ...);
		}(std::make_index_sequence<size_t(N)>{});
		return p;
	}
};

// Stable variants of the wrapper engines: do not accumulate in product space.
// This costs an extra log factor.

template <int N> constexpr std::array<int, size_t(N) * N + 1> matrix_stable_ofs = [] {
	std::array<int, size_t(N) * N + 1> r{};
	for (int i = 0; i <= N * N; i++) r[size_t(i)] = i * N;
	return r;
}();

template <engine E, int N>
struct matrix_stable
		: componentwise<E, mat<typename E::value_type, N>, N * N, matrix_stable_ofs<N>> {
	using base = componentwise<E, mat<typename E::value_type, N>, N * N, matrix_stable_ofs<N>>;
	static constexpr bool commutative = false;
	static constexpr int unit_scale = base::unit_scale;
	template <int A = unit_scale> using transformed_t = typename base::template transformed_t<A>;
	template <int K> using product_t = typename base::template product_t<K>;
	using transformed = typename base::transformed;
	using product = product_t<unit_scale * unit_scale>;

	template <int A, int B>
	static product_t<A * B> mul(const transformed_t<A>& a, const transformed_t<B>& b, int n) {
		product_t<A * B> p;
		// entry (r, c)'s k-th addend a(r,k)*b(k,c), grouped per the offsets
		for (int r = 0; r < N; r++) for (int c = 0; c < N; c++) for (int k = 0; k < N; k++)
			p.t[(size_t(r) * N + c) * N + k] = E::mul(a.t[size_t(r) * N + k], b.t[size_t(k) * N + c], n);
		return p;
	}
	template <int A> static auto sq(const transformed_t<A>& a, int n) { return mul(a, a, n); }
	template <int A1, int B1, int A2, int B2>
	static product_t<A1 * B1 + A2 * B2> mul2(
		const transformed_t<A1>& a1, const transformed_t<B1>& b1,
		const transformed_t<A2>& a2, const transformed_t<B2>& b2,
		int n
	) {
		product_t<A1 * B1 + A2 * B2> p;
		for (int r = 0; r < N; r++) for (int c = 0; c < N; c++) for (int k = 0; k < N; k++)
			p.t[(size_t(r) * N + c) * N + k] = E::mul2(
				a1.t[size_t(r) * N + k], b1.t[size_t(k) * N + c],
				a2.t[size_t(r) * N + k], b2.t[size_t(k) * N + c],
				n
			);
		return p;
	}
};

template <int N> constexpr std::array<int, size_t(N) + 1> trunc_series_stable_ofs = [] {
	std::array<int, size_t(N) + 1> r{};
	for (int i = 0; i <= N; i++) r[size_t(i)] = i * (i + 1) / 2;
	return r;
}();

template <engine E, int N>
struct trunc_stable
		: componentwise<E, trunc_series<typename E::value_type, N>, N, trunc_series_stable_ofs<N>> {
	using base = componentwise<E, trunc_series<typename E::value_type, N>, N, trunc_series_stable_ofs<N>>;
	static constexpr bool commutative = E::commutative;
	static constexpr int unit_scale = base::unit_scale;
	template <int A = unit_scale> using transformed_t = typename base::template transformed_t<A>;
	template <int K> using product_t = typename base::template product_t<K>;
	using transformed = typename base::transformed;
	using product = product_t<unit_scale * unit_scale>;

	template <int A, int B>
	static product_t<A * B> mul(const transformed_t<A>& a, const transformed_t<B>& b, int n) {
		product_t<A * B> p;
		for (int s = 0; s < N; s++) for (int i = 0; i <= s; i++)
			p.t[size_t(trunc_series_stable_ofs<N>[size_t(s)] + i)] = E::mul(a.t[size_t(i)], b.t[size_t(s - i)], n);
		return p;
	}
	template <int A> static auto sq(const transformed_t<A>& a, int n) { return mul(a, a, n); }
	template <int A1, int B1, int A2, int B2>
	static product_t<A1 * B1 + A2 * B2> mul2(
		const transformed_t<A1>& a1, const transformed_t<B1>& b1,
		const transformed_t<A2>& a2, const transformed_t<B2>& b2,
		int n
	) {
		product_t<A1 * B1 + A2 * B2> p;
		for (int s = 0; s < N; s++) for (int i = 0; i <= s; i++)
			p.t[size_t(trunc_series_stable_ofs<N>[size_t(s)] + i)] = E::mul2(
				a1.t[size_t(i)], b1.t[size_t(s - i)],
				a2.t[size_t(i)], b2.t[size_t(s - i)],
				n
			);
		return p;
	}
};

/* namespace wala::fft::engines */ }
