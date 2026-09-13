#pragma once

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

#include "fft/engine.hpp"

namespace wala::fft {

// ==== multiply layer ====
// These are free functions to convolve spans.
//
// The interfaces will typically take input spans, an output span, and an Op representing how to fold the result into the output.
// Output spans may alias one of the input spans.
// Output spans may be shorter than expected; the output will just be truncated.
//
// Some functions may also take E::transformed& objects associated with the input
// spans. These will be lazily filled (see E::extend_to) and used if available.
//
// A trailing E::transformed& keep receives the pointwise product as the
// transform of the (full, untruncated) result when the engine's product is a
// transform (unit_scale == 0), and is cleared otherwise.
//
// Value types build their products from the staged pieces (conv_size_for,
// E::extend_to, E::mul, finish_linear / finish_middle) so they can size their
// output between transforming the inputs and finishing: the output may then
// alias an input even when it has to grow.

// Trailing tag requesting a result that carries its transform.
struct keep_t {};
inline constexpr keep_t keep{};

// Circular convolution mod n (power of 2)
template <engine E, fold_op<typename E::value_type> Op = assign_op>
void multiply_circular(std::span<const typename E::value_type> a, std::span<const typename E::value_type> b,
		std::span<typename E::value_type> out, int n, Op op = {}) {
	assert(!(n & (n-1)));
	auto ta = E::transform(a, n);
	auto tb = E::transform(b, n);
	E::finish(E::mul(ta, tb, n), out, op);
}

template <engine E, fold_op<typename E::value_type> Op = assign_op>
void square_circular(std::span<const typename E::value_type> a, std::span<typename E::value_type> out, int n, Op op = {}) {
	assert(!(n & (n-1)));
	auto ta = E::transform(a, n);
	E::finish(E::sq(ta, n), out, op);
}

// Arrays of length 2^k + 1 are somewhat common, so we will optimize them by
// multiplying mod 2^k, and fixing up the leading coefficient.
// conv_size_for(s) is the transform size n for a linear product of s
// coefficients, and whether that cut applies (n == s - 1, and the result is
// s <= 2n coefficients, so it is still a valid input to E::extend_to).
struct conv_size { int n; bool cut; };
inline conv_size conv_size_for(int s) {
	int n = nextPow2(s);
	bool cut = (n == 2 * (s - 1));
	return {cut ? n / 2 : n, cut};
}

namespace detail {
// Call op while lazily applying the correction if necessary
template <typename T, typename Op>
void emit_linear(std::span<T> buf, int n, int s, bool cut, T c0, std::span<T> out, Op op) {
	T cn{};
	if (cut) {
		cn = buf[0] - c0;
		buf[0] = c0;
	}
	int lim = min(sz(out), min(s, n));
	for (int i = 0; i < lim; i++) op(out[i], buf[i]);
	if (cut && sz(out) >= s) op(out[s-1], cn);
}

// Applies op, diverting the wrapped leading coefficient of a cut product:
// out[0] receives c0 and the wraparound term is captured into cn for the
// caller to emit at out[s-1].
template <typename T, typename Op>
struct cut_op {
	Op op;
	T* out0;
	T c0;
	T& cn;
	void operator()(T& x, T v) const {
		if (&x == out0) { cn = v - c0; v = c0; }
		op(x, v);
	}
};

/* namespace detail */ }

// Finish the size-n circular product p of inputs whose linear product has s
// coefficients ({n, cut} = conv_size_for(s), c0 = a[0] * b[0]) into out,
// applying the cut correction in place.
template <engine E, typename P, fold_op<typename E::value_type> Op = assign_op>
void finish_linear(
	P&& p, int n, int s, bool cut,
	typename E::value_type c0, std::span<typename E::value_type> out, Op op = {}
) {
	using T = typename E::value_type;
	if (sz(out) == 0) return;
	int lim = min(sz(out), min(s, n));
	if (!cut) {
		E::finish(std::move(p), out.subspan(0, lim), op);
	} else {
		T cn{};
		E::finish(std::move(p), out.subspan(0, lim), detail::cut_op<T, Op>{op, &out[0], c0, cn});
		if (sz(out) >= s) op(out[s-1], cn);
	}
}
template <engine E, typename P>
void finish_linear(
	P&& p, int n, int s, bool cut,
	typename E::value_type c0, std::span<typename E::value_type> out, transformed<E>& keep_
) {
	if constexpr (std::same_as<typename E::product, transformed<E>>) {
		keep_ = p;
	} else {
		keep_ = transformed<E>{};
	}
	finish_linear<E>(std::move(p), n, s, cut, c0, out);
}

template <engine E, fold_op<typename E::value_type> Op = assign_op>
void multiply(std::span<const typename E::value_type> a, std::span<const typename E::value_type> b,
		std::span<typename E::value_type> out, Op op = {}) {
	using T = typename E::value_type;
	if (sz(a) == 0 || sz(b) == 0) return;
	int s = sz(a) + sz(b) - 1;
	auto [n, cut] = conv_size_for(s);
	T c0 = a[0] * b[0];
	auto buf = buffer_pool<T>::get(n);
	multiply_circular<E>(a, b, buf.span(), n);
	detail::emit_linear<T>(buf.span(), n, s, cut, c0, out, op);
}

template <engine E, fold_op<typename E::value_type> Op = assign_op>
void multiply(std::span<const typename E::value_type> a, transformed<E>& ta,
		std::span<const typename E::value_type> b, transformed<E>& tb,
		std::span<typename E::value_type> out, Op op = {}) {
	using T = typename E::value_type;
	if (sz(a) == 0 || sz(b) == 0) return;
	int s = sz(a) + sz(b) - 1;
	auto [n, cut] = conv_size_for(s);
	T c0 = a[0] * b[0];
	E::extend_to(ta, n, a);
	E::extend_to(tb, n, b);
	finish_linear<E>(E::mul(ta, tb, n), n, s, cut, c0, out, op);
}
template <engine E>
void multiply(std::span<const typename E::value_type> a, transformed<E>& ta,
		std::span<const typename E::value_type> b, transformed<E>& tb,
		std::span<typename E::value_type> out, transformed<E>& keep_) {
	using T = typename E::value_type;
	keep_ = transformed<E>{};
	if (sz(a) == 0 || sz(b) == 0) return;
	int s = sz(a) + sz(b) - 1;
	auto [n, cut] = conv_size_for(s);
	T c0 = a[0] * b[0];
	E::extend_to(ta, n, a);
	E::extend_to(tb, n, b);
	finish_linear<E>(E::mul(ta, tb, n), n, s, cut, c0, out, keep_);
}

template <engine E, fold_op<typename E::value_type> Op = assign_op>
void multiply_add2(std::span<const typename E::value_type> a1, transformed<E>& ta1,
		std::span<const typename E::value_type> b1, transformed<E>& tb1,
		std::span<const typename E::value_type> a2, transformed<E>& ta2,
		std::span<const typename E::value_type> b2, transformed<E>& tb2,
		std::span<typename E::value_type> out, Op op = {}) {
	using T = typename E::value_type;
	assert(sz(a1) > 0 && sz(b1) > 0 && sz(a2) > 0 && sz(b2) > 0);
	int s = sz(a1) + sz(b1) - 1;
	assert(sz(a2) + sz(b2) - 1 == s);
	auto [n, cut] = conv_size_for(s);
	T c0 = a1[0] * b1[0] + a2[0] * b2[0];
	E::extend_to(ta1, n, a1); E::extend_to(tb1, n, b1);
	E::extend_to(ta2, n, a2); E::extend_to(tb2, n, b2);
	finish_linear<E>(E::mul2(ta1, tb1, ta2, tb2, n), n, s, cut, c0, out, op);
}
template <engine E>
void multiply_add2(std::span<const typename E::value_type> a1, transformed<E>& ta1,
		std::span<const typename E::value_type> b1, transformed<E>& tb1,
		std::span<const typename E::value_type> a2, transformed<E>& ta2,
		std::span<const typename E::value_type> b2, transformed<E>& tb2,
		std::span<typename E::value_type> out, transformed<E>& keep_) {
	using T = typename E::value_type;
	assert(sz(a1) > 0 && sz(b1) > 0 && sz(a2) > 0 && sz(b2) > 0);
	int s = sz(a1) + sz(b1) - 1;
	assert(sz(a2) + sz(b2) - 1 == s);
	auto [n, cut] = conv_size_for(s);
	T c0 = a1[0] * b1[0] + a2[0] * b2[0];
	E::extend_to(ta1, n, a1); E::extend_to(tb1, n, b1);
	E::extend_to(ta2, n, a2); E::extend_to(tb2, n, b2);
	finish_linear<E>(E::mul2(ta1, tb1, ta2, tb2, n), n, s, cut, c0, out, keep_);
}

template <engine E, fold_op<typename E::value_type> Op = assign_op>
void square(std::span<const typename E::value_type> a, std::span<typename E::value_type> out, Op op = {}) {
	using T = typename E::value_type;
	if (sz(a) == 0) return;
	int s = 2 * sz(a) - 1;
	auto [n, cut] = conv_size_for(s);
	T c0 = a[0] * a[0];
	auto buf = buffer_pool<T>::get(n);
	square_circular<E>(a, buf.span(), n);
	detail::emit_linear<T>(buf.span(), n, s, cut, c0, out, op);
}

template <engine E, fold_op<typename E::value_type> Op = assign_op>
void square(std::span<const typename E::value_type> a, transformed<E>& ta,
		std::span<typename E::value_type> out, Op op = {}) {
	using T = typename E::value_type;
	if (sz(a) == 0) return;
	int s = 2 * sz(a) - 1;
	auto [n, cut] = conv_size_for(s);
	T c0 = a[0] * a[0];
	E::extend_to(ta, n, a);
	finish_linear<E>(E::sq(ta, n), n, s, cut, c0, out, op);
}
template <engine E>
void square(std::span<const typename E::value_type> a, transformed<E>& ta,
		std::span<typename E::value_type> out, transformed<E>& keep_) {
	using T = typename E::value_type;
	keep_ = transformed<E>{};
	if (sz(a) == 0) return;
	int s = 2 * sz(a) - 1;
	auto [n, cut] = conv_size_for(s);
	T c0 = a[0] * a[0];
	E::extend_to(ta, n, a);
	finish_linear<E>(E::sq(ta, n), n, s, cut, c0, out, keep_);
}

template <engine E> vector<typename E::value_type> multiply(
		const vector<typename E::value_type>& a, const vector<typename E::value_type>& b) {
	using T = typename E::value_type;
	if (sz(a) == 0 || sz(b) == 0) return {};
	vector<T> r(sz(a) + sz(b) - 1);
	multiply<E>(std::span<const T>(a), std::span<const T>(b), std::span<T>(r));
	return r;
}

template <engine E> vector<typename E::value_type> square(const vector<typename E::value_type>& a) {
	using T = typename E::value_type;
	if (sz(a) == 0) return {};
	vector<T> r(2 * sz(a) - 1);
	square<E>(std::span<const T>(a), std::span<T>(r));
	return r;
}

namespace detail {
// emit_linear but for middle_product
template <typename T, typename Op>
void emit_middle(std::span<T> buf, bool cut, int la, int lb, T c0, T ctop, std::span<T> out, Op op) {
	int m = la - lb + 1;
	T cn{};
	if (cut) {
		cn = buf[0] - c0; // for lb == 1 these coincide: slot 0 = c_0 + c_n and ctop = c_n
		buf[lb - 1] -= ctop;
	}
	int lim = min(sz(out), cut ? m - 1 : m);
	for (int t = 0; t < lim; t++) op(out[t], buf[lb - 1 + t]);
	if (cut && sz(out) >= m) op(out[m-1], cn);
}
/* namespace detail */ }

// The middle product's dot-product case (sz(a) == sz(b)), which needs no transforms.
template <typename T, fold_op<T> Op = assign_op>
void middle_product_dot(std::span<const T> a, std::span<const T> b, std::span<T> out, Op op = {}) {
	assert(sz(a) == sz(b));
	T r{};
	for (int i = 0; i < sz(a); i++) {
		r += a[i] * b[sz(b) - 1 - i];
	}
	if (sz(out) > 0) op(out[0], r);
}

// Finish the size-n circular product p of a (length la) and b (length lb < la)
// as their middle product ({n, cut} = conv_size_for(la), c0 = a[0] * b[0],
// ctop = a[la-1] * b[lb-1]).
template <engine E, typename P, fold_op<typename E::value_type> Op = assign_op>
void finish_middle(
	P&& p, int n, bool cut, int la, int lb,
	typename E::value_type c0, typename E::value_type ctop,
	std::span<typename E::value_type> out, Op op = {}
) {
	using T = typename E::value_type;
	auto buf = buffer_pool<T>::get(n);
	E::finish(std::move(p), buf.span());
	detail::emit_middle<T>(buf.span(), cut, la, lb, c0, ctop, out, op);
}

// Middle product (the transposed multiplication): takes only coefficients of a * b which include terms from all of b.
// Must have len(a) >= len(b)
template <engine E, fold_op<typename E::value_type> Op = assign_op>
void middle_product(std::span<const typename E::value_type> a, std::span<const typename E::value_type> b,
		std::span<typename E::value_type> out, Op op = {}) {
	using T = typename E::value_type;
	if (sz(a) == 0 || sz(b) == 0) return;
	assert(sz(a) >= sz(b));
	if (sz(a) == sz(b)) return middle_product_dot<T>(a, b, out, op);
	auto [n, cut] = conv_size_for(sz(a));
	auto buf = buffer_pool<T>::get(n);
	multiply_circular<E>(a, b, buf.span(), n);
	detail::emit_middle<T>(buf.span(), cut, sz(a), sz(b),
			a[0] * b[0], a[sz(a) - 1] * b[sz(b) - 1], out, op);
}

// TODO: Let's decide whether to keep vector<> returning forms or not; this
// largely depends on whether we think these functions are a public interface or
// merely convenience for value type implementors.
template <engine E> vector<typename E::value_type> middle_product(
		std::span<const typename E::value_type> a, std::span<const typename E::value_type> b) {
	using T = typename E::value_type;
	if (sz(a) == 0 || sz(b) == 0) return {};
	assert(sz(a) >= sz(b));
	vector<T> r(size_t(sz(a) - sz(b) + 1));
	middle_product<E>(a, b, std::span<T>(r));
	return r;
}

template <engine E, fold_op<typename E::value_type> Op = assign_op>
void middle_product(std::span<const typename E::value_type> a, transformed<E>& ta,
		std::span<const typename E::value_type> b, transformed<E>& tb,
		std::span<typename E::value_type> out, Op op = {}) {
	using T = typename E::value_type;
	if (sz(a) == 0 || sz(b) == 0) return;
	assert(sz(a) >= sz(b));
	if (sz(a) == sz(b)) return middle_product_dot<T>(a, b, out, op);
	auto [n, cut] = conv_size_for(sz(a));
	E::extend_to(ta, n, a);
	E::extend_to(tb, n, b);
	finish_middle<E>(E::mul(ta, tb, n), n, cut, sz(a), sz(b),
			a[0] * b[0], a[sz(a) - 1] * b[sz(b) - 1], out, op);
}

template <engine E>
vector<typename E::value_type> middle_product(std::span<const typename E::value_type> a, transformed<E>& ta,
		std::span<const typename E::value_type> b, transformed<E>& tb) {
	using T = typename E::value_type;
	if (sz(a) == 0 || sz(b) == 0) return {};
	assert(sz(a) >= sz(b));
	vector<T> r(size_t(sz(a) - sz(b) + 1));
	middle_product<E>(a, ta, b, tb, std::span<T>(r));
	return r;
}

/* namespace wala::fft */ }
