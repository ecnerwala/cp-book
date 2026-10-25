#pragma once

#include <algorithm>
#include <concepts>
#include <initializer_list>
#include <memory>
#include <ranges>
#include <utility>

#include "containers/vec.hpp"

namespace wala {

// A bounded_vec stored as (base, top, lim) pointers instead of (sz, cap), for hot stacks.
// Shares its interface with bounded_vec through detail::bounded_ops.
template <typename T> struct bounded_stack : detail::bounded_ops {
	static_assert(!std::is_const_v<T>);
	using value_type = T;
	using iterator = detail::contiguous_iterator_impl<bounded_stack, false>;
	using const_iterator = detail::contiguous_iterator_impl<bounded_stack, true>;

	T* base = nullptr;
	T* top = nullptr;
	T* lim = nullptr;

	bounded_stack() = default;

	friend void swap(bounded_stack& a, bounded_stack& b) noexcept { std::swap(a.base, b.base); std::swap(a.top, b.top); std::swap(a.lim, b.lim); }
	bounded_stack(bounded_stack&& o) noexcept : bounded_stack() { swap(*this, o); }
	bounded_stack& operator=(bounded_stack&& o) noexcept { swap(*this, o); return *this; }

	~bounded_stack() { std::destroy(base, top); detail::deallocate(base, int(lim - base)); }

	explicit bounded_stack(with_capacity c) : base(detail::allocate<T>(c.n)), top(base), lim(base + c.n) {}
	template <typename... A> requires detail::uninitialized_constructor<T, A...>
	explicit bounded_stack(with_capacity c, int n, const A&... a) : bounded_stack(c) { grow_to(n, a...); }
	template <typename... A> requires detail::uninitialized_constructor<T, A...>
	explicit bounded_stack(int n, const A&... a) : bounded_stack(with_capacity{n}, n, a...) {}
	explicit bounded_stack(with_capacity c, int n, const T& v) requires std::copy_constructible<T> : bounded_stack(c) { grow_to(n, v); }
	explicit bounded_stack(int n, const T& v) requires std::copy_constructible<T> : bounded_stack(with_capacity{n}, n, v) {}
	template <std::ranges::input_range R>
	requires std::constructible_from<T, std::ranges::range_reference_t<R>>
	bounded_stack(with_capacity c, std::from_range_t, R&& r) : bounded_stack(c) { append_range(std::forward<R>(r)); }
	template <std::ranges::sized_range R>
	requires std::constructible_from<T, std::ranges::range_reference_t<R>>
	bounded_stack(std::from_range_t f, R&& r) : bounded_stack(with_capacity{int(std::ranges::size(r))}, f, r) {}
	bounded_stack(std::initializer_list<T> list) requires std::copy_constructible<T> : bounded_stack(std::from_range, list) {}
	bounded_stack(with_capacity c, std::initializer_list<T> list) requires std::copy_constructible<T> : bounded_stack(c, std::from_range, list) {}
	bounded_stack(const bounded_stack& o) requires std::copy_constructible<T> : bounded_stack(with_capacity{o.capacity()}, std::from_range, o) {}
	bounded_stack& operator=(const bounded_stack& o) requires std::copy_constructible<T> { *this = bounded_stack(o); return *this; }

	explicit bounded_stack(detail::allocation<T> a) : base(a.p), top(a.p + a.sz), lim(a.p + a.cap) {
		debug_assert(0 <= a.sz && a.sz <= a.cap, "adopting allocation with size {} out of range [0, {}]", a.sz, a.cap);
	}
	detail::allocation<T> release() {
		detail::allocation<T> a{base, size(), capacity()};
		base = top = lim = nullptr;
		return a;
	}

	T* end_ptr() const { return top; }
	T* cap_ptr() const { return lim; }
	void set_end(T* p) { top = p; }
	const bounded_stack* owner() const { return this; }
};

} // namespace wala
