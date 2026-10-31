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
// Kept separate from the (sz, base) mixin in vec.hpp so neither representation pays for the other.
template <typename T> struct bounded_stack : detail::insert_ops {
	static_assert(!std::is_const_v<T>);
	using value_type = T;
	using size_type = int;
	using iterator = detail::contiguous_iterator_impl<bounded_stack, false>;
	using const_iterator = detail::contiguous_iterator_impl<bounded_stack, true>;

	T* base = nullptr;
	T* top = nullptr;
	T* lim = nullptr;

	bounded_stack() = default;

	friend void swap(bounded_stack& a, bounded_stack& b) noexcept { std::swap(a.base, b.base); std::swap(a.top, b.top); std::swap(a.lim, b.lim); }
	bounded_stack(bounded_stack&& o) noexcept : bounded_stack() { swap(*this, o); }
	bounded_stack& operator=(bounded_stack&& o) noexcept { swap(*this, o); return *this; }

	// Destruction
	~bounded_stack() { std::destroy(base, top); detail::deallocate(base, int(lim - base)); }

	// Construction
	explicit bounded_stack(with_capacity c) : base(detail::allocate<T>(c.n)), top(base), lim(base + c.n) {}
	template <typename... A> requires detail::uninitialized_constructor<T, A...>
	explicit bounded_stack(with_capacity c, size_type n, const A&... a) : bounded_stack(c) {
		debug_assert(n <= c.n, "size {} exceeds capacity {}", n, c.n);
		detail::uninitialized_construct(base, base + n, a...); top = base + n;
	}
	template <typename... A> requires detail::uninitialized_constructor<T, A...>
	explicit bounded_stack(size_type n, const A&... a) : bounded_stack(with_capacity{n}, n, a...) {}
	explicit bounded_stack(with_capacity c, size_type n, const T& v) requires std::copy_constructible<T> : bounded_stack(c) { grow_by<T>(n, v); }
	explicit bounded_stack(size_type n, const T& v) requires std::copy_constructible<T> : bounded_stack(with_capacity{n}, n, v) {}
	template <std::ranges::input_range R>
	requires std::constructible_from<T, std::ranges::range_reference_t<R>>
	bounded_stack(with_capacity c, std::from_range_t, R&& r) : bounded_stack(c) {
		auto [in, out] = std::ranges::uninitialized_copy(r, std::ranges::subrange(base, lim));
		// Must have drained all of r
		debug_assert(in == std::ranges::end(r), "range size exceeds capacity {}", c.n);
		top = out;
	}
	template <std::ranges::sized_range R>
	requires std::constructible_from<T, std::ranges::range_reference_t<R>>
	bounded_stack(std::from_range_t f, R&& r) : bounded_stack(with_capacity{size_type(std::ranges::size(r))}, f, r) {}
	bounded_stack(std::initializer_list<T> list) requires std::copy_constructible<T> : bounded_stack(std::from_range, list) {}
	bounded_stack(with_capacity c, std::initializer_list<T> list) requires std::copy_constructible<T> : bounded_stack(c, std::from_range, list) {}
	bounded_stack(const bounded_stack& o) requires std::copy_constructible<T> : bounded_stack(with_capacity{o.capacity()}, std::from_range, o) {}
	bounded_stack& operator=(const bounded_stack& o) requires std::copy_constructible<T> { *this = bounded_stack(o); return *this; }

	// Conversions, transferring the allocation
	explicit bounded_stack(detail::allocation<T> a) : base(a.p), top(a.p + a.sz), lim(a.p + a.cap) {
		debug_assert(0 <= a.sz && a.sz <= a.cap, "adopting allocation with size {} out of range [0, {}]", a.sz, a.cap);
	}
	detail::allocation<T> release() {
		detail::allocation<T> a{base, size_type(top - base), size_type(lim - base)};
		base = top = lim = nullptr;
		return a;
	}
	vec<T> into_vec() && {
		debug_assert(top == lim, "converting container with size {} and capacity {} to vec", top - base, lim - base);
		return vec<T>(release());
	}
	bounded_vec<T> into_bounded() && { return bounded_vec<T>(release()); }
	growable_vec<T> into_growable() && { return growable_vec<T>(release()); }

	[[nodiscard]] bool empty() const { return top == base; }
	[[nodiscard]] size_type size() const { return size_type(top - base); }
	[[nodiscard]] size_type capacity() const { return size_type(lim - base); }

	void assert_valid_deref(const T* p) const {
		debug_assert(base <= p && p < top, "iterator position {} out of range [0, {})", p - base, top - base);
	}
	void assert_valid_index(size_type i) const {
		debug_assert(0 <= i && base + i < top, "index {} out of range [0, {})", i, top - base);
	}
	void assert_valid_insert(size_type i) const {
		debug_assert(0 <= i && base + i <= top, "index {} out of range [0, {}]", i, top - base);
	}
	void assert_nonempty() const {
		debug_assert(top > base, "container is empty");
	}

	T* data() { return base; }
	const T* data() const { return base; }
	iterator begin() { return iterator(base, this); }
	const_iterator begin() const { return const_iterator(base, this); }
	iterator end() { return iterator(top, this); }
	const_iterator end() const { return const_iterator(top, this); }
	T& operator [] (size_type i) { assert_valid_index(i); return base[i]; }
	const T& operator [] (size_type i) const { assert_valid_index(i); return base[i]; }
	T& front() { assert_nonempty(); return base[0]; }
	const T& front() const { assert_nonempty(); return base[0]; }
	T& back() { assert_nonempty(); return top[-1]; }
	const T& back() const { assert_nonempty(); return top[-1]; }

	// Resizing ops
	template <typename... Args> requires std::constructible_from<T, Args&&...>
	T& emplace_back(Args&&... args) {
		debug_assert(top < lim, "pushing into full container (capacity {})", lim - base);
		return *std::construct_at(top++, std::forward<Args>(args)...);
	}
	template <std::ranges::input_range R>
	requires std::constructible_from<T, std::ranges::range_reference_t<R>>
	std::ranges::subrange<iterator> append_range(R&& r) {
		auto old_end = end();
		// Use raw pointers since the memory is uninitialized
		auto [in, out] = std::ranges::uninitialized_copy(r, std::ranges::subrange(top, lim));
		// Must have drained all of r
		debug_assert(in == std::ranges::end(r), "range size exceeds capacity {}", lim - base);
		top = out;
		return {old_end, end()};
	}
	template <typename... A> requires detail::uninitialized_constructor<T, A...>
	std::ranges::subrange<iterator> grow_by(size_type n, const A&... a) {
		debug_assert(0 <= n, "growing container by negative {}", n);
		debug_assert(top + n <= lim, "growing container with size {} by {} exceeds capacity {}", top - base, n, lim - base);
		auto old_end = end();
		detail::uninitialized_construct(top, top + n, a...);
		top += n;
		return {old_end, end()};
	}
	std::ranges::subrange<iterator> grow_by(size_type n, const T& v) requires std::copy_constructible<T> { return grow_by<T>(n, v); }
	template <typename... A> requires detail::uninitialized_constructor<T, A...>
	std::ranges::subrange<iterator> grow_to(size_type n, const A&... a) {
		debug_assert(top <= base + n, "growing container with size {} to smaller target {}", top - base, n);
		debug_assert(base + n <= lim, "growing container to {} exceeds capacity {}", n, lim - base);
		auto old_end = end();
		detail::uninitialized_construct(top, base + n, a...);
		top = base + n;
		return {old_end, end()};
	}
	std::ranges::subrange<iterator> grow_to(size_type n, const T& v) requires std::copy_constructible<T> { return grow_to<T>(n, v); }
	template <typename... A> requires detail::uninitialized_constructor<T, A...>
	std::ranges::subrange<iterator> grow_to(iterator it, const A&... a) {
		debug_assert(top <= it.p, "growing container with size {} to smaller target {}", top - base, it.p - base);
		debug_assert(it.p <= lim, "growing container to {} exceeds capacity {}", it.p - base, lim - base);
		auto old_end = end();
		detail::uninitialized_construct(top, it.p, a...);
		top = it.p;
		return {old_end, end()};
	}
	std::ranges::subrange<iterator> grow_to(iterator it, const T& v) requires std::copy_constructible<T> { return grow_to<T>(it, v); }

	// Remove from the back
	void pop_back() { assert_nonempty(); std::destroy_at(--top); }
	void shrink_by(size_type n) {
		debug_assert(0 <= n && base <= top - n, "shrinking container with size {} by {}", top - base, n);
		std::destroy(top - n, top);
		top -= n;
	}
	void shrink_to(size_type n) {
		debug_assert(0 <= n && base + n <= top, "shrinking container with size {} to {}", top - base, n);
		std::destroy(base + n, top);
		top = base + n;
	}
	void shrink_to(iterator it) {
		debug_assert(base <= it.p && it.p <= top, "shrinking container with size {} to {}", top - base, it.p - base);
		std::destroy(it.p, top);
		top = it.p;
	}

	void clear() { shrink_to(0); }

	friend bool operator == (const bounded_stack& a, const bounded_stack& b) { return std::ranges::equal(a, b); }
	friend auto operator <=> (const bounded_stack& a, const bounded_stack& b) {
		return std::lexicographical_compare_three_way(a.begin(), a.end(), b.begin(), b.end());
	}
};

} // namespace wala
