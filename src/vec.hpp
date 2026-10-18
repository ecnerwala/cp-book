#pragma once

#include <algorithm>
#include <concepts>
#include <cstdlib>
#include <initializer_list>
#include <iterator>
#include <memory>
#include <ranges>
#include <span>
#include <variant>

namespace wala {

#ifdef _GLIBCXX_DEBUG
constexpr bool WALA_DEBUG = true;
#else
constexpr bool WALA_DEBUG = false;
#endif

inline void debug_assert(bool b) {
	if (!b) [[unlikely]] {
		if constexpr (WALA_DEBUG) {
			std::abort();
		} else {
			__builtin_unreachable();
		}
	}
}

struct uninit_t{}; inline constexpr uninit_t uninit{};
struct with_capacity{ int n; };

namespace detail {
template <typename C, bool Const>
struct contiguous_iterator_impl {
	using Self = contiguous_iterator_impl;

	using iterator_category = std::random_access_iterator_tag;
	using iterator_concept = std::contiguous_iterator_tag;
	// value_type should be non-cv qualified
	using value_type = C::value_type;
	using T = std::conditional_t<Const, const value_type, value_type>;
	using pointer = T*;
	using reference = T&;
	using difference_type = int;

	pointer p = nullptr;
	[[no_unique_address]] std::conditional_t<WALA_DEBUG, const C*, std::monostate> c{};

	contiguous_iterator_impl() = default;
	contiguous_iterator_impl(T* p_, [[maybe_unused]] const C* c_) : p(p_) { if constexpr (WALA_DEBUG) c = c_; }
	// Cast from non-const to const
	template <bool OConst> requires (Const && !OConst) contiguous_iterator_impl(const contiguous_iterator_impl<C, OConst>& o)  : p(o.p), c(o.c) {}

	void assert_valid_deref() const { if constexpr (WALA_DEBUG) { c->assert_valid_deref(p); } }
	static void assert_same_container([[maybe_unused]] Self a, [[maybe_unused]] Self b) { if constexpr (WALA_DEBUG) { debug_assert(a.c == b.c); } }

	T& operator * () const { assert_valid_deref(); return *p; }
	T* operator -> () const { assert_valid_deref(); return p; }
	T& operator [] (int n) const { return *(*this + n); }

	Self& operator ++ () { ++p; return *this; }
	Self operator ++ (int) { auto tmp = *this; operator ++ (); return tmp; }
	Self& operator -- () { --p; return *this; }
	Self operator -- (int) { auto tmp = *this; operator -- (); return tmp; }
	Self& operator += (int n) { p += n; return *this; }
	friend Self operator + (Self a, int n) { return a += n; }
	friend Self operator + (int n, Self a) { return a + n; }
	Self& operator -= (int n) { p -= n; return *this; }
	friend Self operator - (Self a, int n) { return a -= n; }
	friend int operator - (Self a, Self b) { assert_same_container(a, b); return int(a.p - b.p); }

	friend bool operator == (Self a, Self b) { assert_same_container(a, b); return a.p == b.p; }
	friend auto operator <=> (Self a, Self b) { assert_same_container(a, b); return a.p <=> b.p; }
};
} // namespace detail
} // namespace wala

template <typename C, bool Const> struct std::pointer_traits<wala::detail::contiguous_iterator_impl<C, Const>> {
	using pointer = wala::detail::contiguous_iterator_impl<C, Const>;
	using element_type = pointer::T;
	using difference_type = pointer::difference_type;
	static element_type* to_address(pointer it) noexcept { return it.p; }
};

namespace wala {

namespace detail {
template <typename T> T* allocate(int n) {
	debug_assert(n >= 0);
	return n ? std::allocator<T>().allocate(std::size_t(n)) : nullptr;
}
template <typename T> void deallocate(T* p, int n) {
	if (p) std::allocator<T>().deallocate(p, std::size_t(n));
}
} // namespace detail

template <typename T> struct vec {
	static_assert(!std::is_const_v<T>);
	using value_type = T;
	using size_type = int;
	using difference_type = int;
	using iterator = detail::contiguous_iterator_impl<vec, false>;
	using const_iterator = detail::contiguous_iterator_impl<vec, true>;

	size_type sz = 0;
	T* base = nullptr;

	vec() = default;

	friend void swap(vec& a, vec& b) noexcept { std::swap(a.sz, b.sz); std::swap(a.base, b.base); }
	vec(vec&& o) noexcept : vec() { swap(*this, o); }
	vec& operator=(vec&& o) noexcept { swap(*this, o); return *this; }
	vec(const vec&) = delete;
	vec& operator=(const vec&) = delete;

	~vec() { std::destroy_n(base, sz); detail::deallocate(base, sz); }
private:
	struct alloc_only { size_type n; };
	explicit vec(alloc_only a) : sz(a.n), base(detail::allocate<T>(a.n)) {}
public:
	explicit vec(size_type n, uninit_t) requires std::is_trivially_default_constructible_v<T> : vec(alloc_only{n}) {}
	explicit vec(size_type n) : vec(alloc_only{n}) { std::uninitialized_value_construct_n(base, sz); }
	explicit vec(size_type n, const T& v) : vec(alloc_only{n}) { std::uninitialized_fill_n(base, sz, v); }

	template <std::ranges::sized_range R> requires std::constructible_from<T, std::ranges::range_reference_t<R>>
	vec(std::from_range_t, R&& r) : vec(alloc_only{size_type(std::ranges::size(r))}) {
		std::ranges::uninitialized_copy(r, std::span(base, sz));
	}

	vec(std::initializer_list<T> list) requires std::copy_constructible<T> : vec(std::from_range, list) {}

	[[nodiscard]] vec clone() const requires std::copy_constructible<T> { return vec(std::from_range, *this); }

	[[nodiscard]] bool empty() const { return !sz; }
	[[nodiscard]] size_type size() const { return sz; }

	void assert_valid_deref(const T* p) const { debug_assert(base <= p); debug_assert(p < base + sz); }
	void assert_valid_index(size_type i) const { debug_assert(0 <= i); debug_assert(i < sz); }
	void assert_nonempty() const { debug_assert(sz > 0); }

	T* data() { return base; }
	const T* data() const { return base; }
	iterator begin() { return iterator(base, this); }
	const_iterator begin() const { return const_iterator(base, this); }
	iterator end() { return iterator(base + sz, this); }
	const_iterator end() const { return const_iterator(base + sz, this); }
	T& operator [] (size_type i) { assert_valid_index(i); return base[i]; }
	const T& operator [] (size_type i) const { assert_valid_index(i); return base[i]; }
	T& front() { assert_nonempty(); return base[0]; }
	const T& front() const { assert_nonempty(); return base[0]; }
	T& back() { assert_nonempty(); return base[sz-1]; }
	const T& back() const { assert_nonempty(); return base[sz-1]; }

	friend bool operator == (const vec& a, const vec& b) { return std::ranges::equal(a, b); }
	friend auto operator <=> (const vec& a, const vec& b) {
		return std::lexicographical_compare_three_way(a.begin(), a.end(), b.begin(), b.end());
	}
};

template <typename T> struct bounded_vec {
	static_assert(!std::is_const_v<T>);
	using value_type = T;
	using size_type = int;
	using difference_type = int;
	using iterator = detail::contiguous_iterator_impl<bounded_vec, false>;
	using const_iterator = detail::contiguous_iterator_impl<bounded_vec, true>;

	size_type sz = 0;
	size_type cap = 0;
	T* base = nullptr;

	bounded_vec() = default;

	friend void swap(bounded_vec& a, bounded_vec& b) noexcept { std::swap(a.sz, b.sz); std::swap(a.cap, b.cap); std::swap(a.base, b.base); }
	bounded_vec(bounded_vec&& o) noexcept : bounded_vec() { swap(*this, o); }
	bounded_vec& operator=(bounded_vec&& o) noexcept { swap(*this, o); return *this; }
	bounded_vec(const bounded_vec&) = delete;
	bounded_vec& operator=(const bounded_vec&) = delete;

	// Destruction
	~bounded_vec() { std::destroy_n(base, sz); detail::deallocate(base, cap); }

	// Construction
	explicit bounded_vec(with_capacity c) : cap(c.n), base(detail::allocate<T>(c.n)) {}
	explicit bounded_vec(size_type n, uninit_t, with_capacity c) requires std::is_trivially_default_constructible_v<T> : bounded_vec(c) { sz = n; }
	explicit bounded_vec(size_type n, uninit_t u) requires std::is_trivially_default_constructible_v<T> : bounded_vec(n, u, with_capacity{n}) {}
	explicit bounded_vec(size_type n, with_capacity c) : bounded_vec(c) { std::uninitialized_value_construct_n(base, sz = n); }
	explicit bounded_vec(size_type n) : bounded_vec(n, with_capacity{n}) {}
	explicit bounded_vec(size_type n, const T& v, with_capacity c) : bounded_vec(c) { std::uninitialized_fill_n(base, sz = n, v); }
	explicit bounded_vec(size_type n, const T& v) : bounded_vec(n, v, with_capacity{n}) {}
	template <std::ranges::input_range R>
	requires std::constructible_from<T, std::ranges::range_reference_t<R>>
	explicit bounded_vec(std::from_range_t, R&& r, with_capacity c) : bounded_vec(c) {
		auto [in, out] = std::ranges::uninitialized_copy(r, std::ranges::subrange(base, base + cap));
		// Must have drained all of r
		debug_assert(in == std::ranges::end(r));
		sz = size_type(out - base);
	}
	template <std::ranges::sized_range R>
	requires std::constructible_from<T, std::ranges::range_reference_t<R>>
	explicit bounded_vec(std::from_range_t f, R&& r) : bounded_vec(f, r, with_capacity{size_type(std::ranges::size(r))}) {}
	bounded_vec(std::initializer_list<T> list) requires std::copy_constructible<T> : bounded_vec(std::from_range, list) {}
	bounded_vec(std::initializer_list<T> list, with_capacity c) requires std::copy_constructible<T> : bounded_vec(std::from_range, list, c) {}

	[[nodiscard]] bounded_vec clone() const requires std::copy_constructible<T> { return bounded_vec(std::from_range, *this, with_capacity{cap}); }

	[[nodiscard]] bool empty() const { return !sz; }
	[[nodiscard]] size_type size() const { return sz; }
	[[nodiscard]] size_type capacity() const { return cap; }

	void assert_valid_deref(const T* p) const { debug_assert(base <= p); debug_assert(p < base + sz); }
	void assert_valid_index(size_type i) const { debug_assert(0 <= i); debug_assert(i < sz); }
	void assert_nonempty() const { debug_assert(sz > 0); }

	T* data() { return base; }
	const T* data() const { return base; }
	iterator begin() { return iterator(base, this); }
	const_iterator begin() const { return const_iterator(base, this); }
	iterator end() { return iterator(base + sz, this); }
	const_iterator end() const { return const_iterator(base + sz, this); }
	T& operator [] (size_type i) { assert_valid_index(i); return base[i]; }
	const T& operator [] (size_type i) const { assert_valid_index(i); return base[i]; }
	T& front() { assert_nonempty(); return base[0]; }
	const T& front() const { assert_nonempty(); return base[0]; }
	T& back() { assert_nonempty(); return base[sz-1]; }
	const T& back() const { assert_nonempty(); return base[sz-1]; }

	// Resizing ops
	template <typename... Args> T& emplace_back(Args&&... args) { debug_assert(sz < cap); return *std::construct_at(base + sz++, std::forward<Args>(args)...); }
	void push_back(const T& v) { emplace_back(v); }
	void push_back(T&& v) { emplace_back(std::move(v)); }
	void pop_back() { debug_assert(sz > 0); std::destroy_at(base + --sz); }
	template <class... Args> T& emplace_front(Args&&... args) { emplace_back(std::forward<Args>(args)...); std::ranges::rotate(*this, end() - 1); return front(); }
	void push_front(const T& v) { emplace_front(v); }
	void push_front(T&& v) { emplace_front(std::move(v)); }
	void pop_front() { std::ranges::rotate(*this, begin() + 1); pop_back(); }
	void clear() { shrink_to(0); }
	void shrink_to(size_type n) { debug_assert(0 <= n && n <= sz); std::destroy(base + n, base + sz); sz = n; }
	void grow_to(size_type n, uninit_t) requires std::is_trivially_default_constructible_v<T> { debug_assert(sz <= n && n <= cap); sz = n; }
	void grow_to(size_type n) { debug_assert(sz <= n && n <= cap); std::uninitialized_value_construct(base + sz, base + n); sz = n; }
	void grow_to(size_type n, const T& v) { debug_assert(sz <= n && n <= cap); std::uninitialized_fill(base + sz, base + n, v); sz = n; }

	friend bool operator == (const bounded_vec& a, const bounded_vec& b) { return std::ranges::equal(a, b); }
	friend auto operator <=> (const bounded_vec& a, const bounded_vec& b) {
		return std::lexicographical_compare_three_way(a.begin(), a.end(), b.begin(), b.end());
	}
};

template <typename T> struct bounded_stack {
	static_assert(!std::is_const_v<T>);
	using value_type = T;
	using size_type = int;
	using difference_type = int;
	using iterator = detail::contiguous_iterator_impl<bounded_stack, false>;
	using const_iterator = detail::contiguous_iterator_impl<bounded_stack, true>;

	T* base = nullptr;
	T* top = nullptr;
	T* lim = nullptr;

	bounded_stack() = default;

	friend void swap(bounded_stack& a, bounded_stack& b) noexcept { std::swap(a.base, b.base); std::swap(a.top, b.top); std::swap(a.lim, b.lim); }
	bounded_stack(bounded_stack&& o) noexcept : bounded_stack() { swap(*this, o); }
	bounded_stack& operator=(bounded_stack&& o) noexcept { swap(*this, o); return *this; }
	bounded_stack(const bounded_stack&) = delete;
	bounded_stack& operator=(const bounded_stack&) = delete;

	// Destruction
	~bounded_stack() { std::destroy(base, top); detail::deallocate(base, int(lim - base)); }

	// Construction
	explicit bounded_stack(with_capacity c) : base(detail::allocate<T>(c.n)), top(base), lim(base + c.n) {}
	explicit bounded_stack(size_type n, uninit_t, with_capacity c) requires std::is_trivially_default_constructible_v<T> : bounded_stack(c) { top = base + n; }
	explicit bounded_stack(size_type n, uninit_t u) requires std::is_trivially_default_constructible_v<T> : bounded_stack(n, u, with_capacity{n}) {}
	explicit bounded_stack(size_type n, with_capacity c) : bounded_stack(c) { std::uninitialized_value_construct(base, top = base + n); }
	explicit bounded_stack(size_type n) : bounded_stack(n, with_capacity{n}) {}
	explicit bounded_stack(size_type n, const T& v, with_capacity c) : bounded_stack(c) { std::uninitialized_fill(base, top = base + n, v); }
	explicit bounded_stack(size_type n, const T& v) : bounded_stack(n, v, with_capacity{n}) {}
	template <std::ranges::input_range R>
	requires std::constructible_from<T, std::ranges::range_reference_t<R>>
	explicit bounded_stack(std::from_range_t, R&& r, with_capacity c) : bounded_stack(c) {
		auto [in, out] = std::ranges::uninitialized_copy(r, std::ranges::subrange(base, lim));
		// Must have drained all of r
		debug_assert(in == std::ranges::end(r));
		top = out;
	}
	template <std::ranges::sized_range R>
	requires std::constructible_from<T, std::ranges::range_reference_t<R>>
	explicit bounded_stack(std::from_range_t f, R&& r) : bounded_stack(f, r, with_capacity{size_type(std::ranges::size(r))}) {}
	bounded_stack(std::initializer_list<T> list) requires std::copy_constructible<T> : bounded_stack(std::from_range, list) {}
	bounded_stack(std::initializer_list<T> list, with_capacity c) requires std::copy_constructible<T> : bounded_stack(std::from_range, list, c) {}

	[[nodiscard]] bounded_stack clone() const requires std::copy_constructible<T> { return bounded_stack(std::from_range, *this, with_capacity{capacity()}); }

	[[nodiscard]] bool empty() const { return top == base; }
	[[nodiscard]] size_type size() const { return size_type(top - base); }
	[[nodiscard]] size_type capacity() const { return size_type(lim - base); }

	void assert_valid_deref(const T* p) const { debug_assert(base <= p); debug_assert(p < top); }
	void assert_valid_index(size_type i) const { debug_assert(0 <= i); debug_assert(i < size()); }
	void assert_nonempty() const { debug_assert(top > base); }

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
	template <typename... Args> T& emplace_back(Args&&... args) { debug_assert(top < lim); return *std::construct_at(top++, std::forward<Args>(args)...); }
	void push_back(const T& v) { emplace_back(v); }
	void push_back(T&& v) { emplace_back(std::move(v)); }
	void pop_back() { debug_assert(top > base); std::destroy_at(--top); }
	template <class... Args> T& emplace_front(Args&&... args) { emplace_back(std::forward<Args>(args)...); std::ranges::rotate(*this, end() - 1); return front(); }
	void push_front(const T& v) { emplace_front(v); }
	void push_front(T&& v) { emplace_front(std::move(v)); }
	void pop_front() { std::ranges::rotate(*this, begin() + 1); pop_back(); }
	void clear() { std::destroy(base, top); top = base; }
	void shrink_by(size_type n) { debug_assert(0 <= n && base <= top - n); std::destroy(top - n, top); top -= n; }
	void shrink_to(size_type n) { debug_assert(0 <= n && base + n <= top); std::destroy(base + n, top); top = base + n; }
	void grow_by(size_type n, uninit_t) requires std::is_trivially_default_constructible_v<T> { debug_assert(0 <= n && top + n <= lim); top += n; }
	void grow_by(size_type n) { debug_assert(0 <= n && top + n <= lim); std::uninitialized_value_construct(top, top + n); top += n; }
	void grow_by(size_type n, const T& v) { debug_assert(0 <= n && top + n <= lim); std::uninitialized_fill(top, top + n, v); top += n; }
	void grow_to(size_type n, uninit_t) requires std::is_trivially_default_constructible_v<T> { debug_assert(top <= base + n && base + n <= lim); top = base + n; }
	void grow_to(size_type n) { debug_assert(top <= base + n && base + n <= lim); std::uninitialized_value_construct(top, base + n); top = base + n; }
	void grow_to(size_type n, const T& v) { debug_assert(top <= base + n && base + n <= lim); std::uninitialized_fill(top, base + n, v); top = base + n; }

	friend bool operator == (const bounded_stack& a, const bounded_stack& b) { return std::ranges::equal(a, b); }
	friend auto operator <=> (const bounded_stack& a, const bounded_stack& b) {
		return std::lexicographical_compare_three_way(a.begin(), a.end(), b.begin(), b.end());
	}
};

} // namespace wala
