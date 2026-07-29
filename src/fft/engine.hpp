#pragma once

#include <concepts>
#include <span>
#include <type_traits>
#include <utility>

#include "fft/common.hpp"

namespace ecnerwala::fft {

// ==== engine concept ====

// Output operations for the finish step to express arbitrary fusion into the output buffer.
struct assign_op { template <typename T> void operator()(T& d, T v) const { d = v; } };
struct add_op { template <typename T> void operator()(T& d, T v) const { d += v; } };
struct sub_op { template <typename T> void operator()(T& d, T v) const { d -= v; } };
struct add_twice_op { template <typename T> void operator()(T& d, T v) const { d += v + v; } };

// `engine` contract
//   engine represents a way of packing/unpacking sequences over an arbitrary ring into FFT-style transforms.
//   We expect transforms/products of transforms to be linear but potentially lossy/imprecise, so we'll track precision
//   at compile-time as a template parameter.
//
//   E::value_type   The ring we operate over
//   E::unit_scale   0 or 1 depending on whether there's error that can accumulate
//   E::commutative  A marker for whether the ring is commutative
//
//   transformed_t<A>  The transform of a sequence. This object owns its data buffer.
//   product_t<A>      The product of 2 transforms. May equal transformed_t, particularly when unit_scale = 0.
//
//   transformed       alias for transformed_t<unit_scale>
//   product           alias for product_t<unit_scale>
//
//   The basic multiplication API is
//      transform(span<const value_type> in, int n) -> transformed_t<unit_scale>
//      mul(transformed_t<A>, transformed_t<B>, int n) -> product_t<A*B>
//      mul2(a1, b1, a2, b2, int n) -> product_t<A1*B1 + A2*B2>, computing a1*b1 + a2*b2 in one pass
//      finish(product_t<A>, span<value_type>& out, Op) -> void
//
//   Input span can be length up to 2n.
//   Output spans can be length up to n; only the prefix that exists is filled.
//   finish applies Op exactly once per out element, in index order, to
//   value_type targets (so ops may be stateful).
//   Transforms can be longer than necessary, and only the relevant prefix is used.
//
//   For non-exact engines, there's some subtlety in whether we wrap before or after packing.
//   We will choose to wrap *after* packing, which hurts error bounds but makes the prefix condition more uniform.
//
//   Additionally, we have APIs to take advantage of linearity in both transformed and product space:
//      add(transformed_t<A>, transformed_t<B>) -> transformed_t<A+B>
//      add(product_t<K1>, product_t<K2>) -> product_t<K1+K2>
//
//   Finally, we expose some additional fast-transform optimization paths.
//   extend_to only operates on transformed_t<unit_scale>; the others are scale-generic.
//   downsample is also defined on product_t (halving before finish saves inverse-transform work).
//      extend_to       build (if empty) or grow a transform to size m by repeated doubling; feed the same coefficients (sz <= 2m) every time
//                      (each doubling step reads only the coefficients that fit, so zero-padded buffers are fine:
//                       coefficients past twice the existing transform's size must be zero, or it couldn't be a prefix)
//      downsample      compute the half-sized transform/product of just the even (odd = false) or odd terms of the input
//      negate_arg      size n transform of A(-x)
template <typename E>
concept engine = requires(
	std::span<const typename E::value_type> in,
	std::span<typename E::value_type> out,
	typename E::transformed& t,
	const typename E::transformed& ct,
	typename E::product& p,
	const typename E::product& cp,
	int n
) {
	typename E::value_type;
	{ E::transform(in, n) } -> std::same_as<typename E::transformed>;
	{ ct.size() } -> std::same_as<int>;
	E::extend_to(t, n, in);
	{ E::downsample(ct, n, false) } -> std::same_as<typename E::transformed>;
	{ E::downsample(cp, n, false) } -> std::same_as<typename E::product>;
	{ E::negate_arg(ct, n) } -> std::same_as<typename E::transformed>;
	{ E::mul(ct, ct, n) } -> std::same_as<typename E::product>;
	{ E::sq(ct, n) } -> std::same_as<typename E::product>;
	{ E::mul2(ct, ct, ct, ct, n) } -> std::same_as<typename E::template product_t<2 * E::unit_scale>>;
	E::finish(std::move(p), out);
	E::finish(std::move(p), out, add_op{});
	E::finish(E::add(std::move(p), std::move(p)), out);
	{ E::add(E::transform(in, n), ct) } -> std::same_as<typename E::template transformed_t<2 * E::unit_scale>>;
	{ E::add(std::move(p), std::move(p)) } -> std::same_as<typename E::template product_t<2 * E::unit_scale>>;
	requires std::same_as<std::remove_cvref_t<decltype(E::commutative)>, bool>;
	requires std::same_as<std::remove_cvref_t<decltype(E::unit_scale)>, int>;
};

// ==== storage/ownership conventions ====
//
// Parameter forms encode clobber permission:
//   const&        borrowed, read-only
//   && / by-value sink: the callee owns it and may clobber (finish works in the product's buffer)
//   & out-param   overwritten, reusing the object's existing capacity
//
// Engines may additionally provide out-parameter forms of the value-returning
// primitives so callers can pre-construct/reserve a transformed/product and
// build into it (contiguous or reused storage across a loop):
//   transform(in, n, transformed& out)
//   mul(a, b, n, product& out)
//   downsample(t, n, odd, out)
//   negate_arg(t, n, transformed& out)
// These are optional; the *_into helpers below dispatch to them when present
// and fall back to move-assignment (fresh allocation) otherwise.
template <typename E>
void transform_into(std::span<const typename E::value_type> in, int n, typename E::transformed& out) {
	if constexpr (requires { E::transform(in, n, out); }) E::transform(in, n, out);
	else out = E::transform(in, n);
}
template <typename E, typename P>
void mul_into(const typename E::transformed& a, const typename E::transformed& b, int n, P& out) {
	if constexpr (requires { E::mul(a, b, n, out); }) E::mul(a, b, n, out);
	else out = E::mul(a, b, n);
}
template <typename E, typename P>
void downsample_into(const P& t, int n, bool odd, P& out) {
	if constexpr (requires { E::downsample(t, n, odd, out); }) E::downsample(t, n, odd, out);
	else out = E::downsample(t, n, odd);
}
template <typename E>
void negate_arg_into(const typename E::transformed& t, int n, typename E::transformed& out) {
	if constexpr (requires { E::negate_arg(t, n, out); }) E::negate_arg(t, n, out);
	else out = E::negate_arg(t, n);
}

// Constrains two engine-parameterized value types to share the same engine.
template <typename A, typename B>
concept same_engine = std::same_as<typename A::engine_t, typename B::engine_t>;

// short spelling for E::transformed at use sites
template <engine E> using transformed = typename E::transformed;

/* namespace ecnerwala::fft */ }
