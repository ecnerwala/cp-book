#pragma once

#include <algorithm>
#include <concepts>
#include <cstdlib>
#include <format>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <memory>
#include <ranges>
#include <source_location>
#include <span>
#include <variant>

namespace wala {

#ifdef _GLIBCXX_DEBUG
constexpr bool WALA_DEBUG = true;
#else
constexpr bool WALA_DEBUG = false;
#endif

template <typename... A> struct fmt_loc {
	std::format_string<A...> fmt;
	std::source_location loc;
	template <typename S> consteval fmt_loc(const S& s, std::source_location l = std::source_location::current()) : fmt(s), loc(l) {}
};
template <typename... A> void debug_assert(bool b, fmt_loc<std::type_identity_t<A>...> f, const A&... a) {
	if (!b) [[unlikely]] {
		if constexpr (WALA_DEBUG) {
			std::cerr << f.loc.file_name() << ':' << f.loc.line() << ": wala: " << std::vformat(f.fmt.get(), std::make_format_args(a...)) << '\n';
			std::abort();
		} else {
			__builtin_unreachable();
		}
	}
}

struct uninit_t{}; inline constexpr uninit_t uninit{};
struct with_capacity{ int n; };

template <typename T> struct vec;
template <typename T> struct bounded_vec;
template <typename T> struct bounded_stack;
template <typename T> struct growable_vec;

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
	static void assert_same_container([[maybe_unused]] Self a, [[maybe_unused]] Self b) {
		if constexpr (WALA_DEBUG) {
			debug_assert(a.c == b.c, "comparing two iterators from different containers");
		}
	}

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
	debug_assert(n >= 0, "allocating negative {} items", n);
	return n ? std::allocator<T>().allocate(std::size_t(n)) : nullptr;
}
template <typename T> void deallocate(T* p, int n) {
	if (p) std::allocator<T>().deallocate(p, std::size_t(n));
}

// Ownership of an existing allocation of cap items, the first sz of which are live.
// Passed between the containers by release() and the allocation constructors.
template <typename T> struct allocation { T* p; int sz; int cap; };

template <typename T> void uninitialized_construct(T* first, T* last) { std::uninitialized_value_construct(first, last); }
template <typename T> void uninitialized_construct(T*, T*, uninit_t) requires std::is_trivially_default_constructible_v<T> {}
template <typename T> void uninitialized_construct(T* first, T* last, const std::type_identity_t<T>& v) { std::uninitialized_fill(first, last, v); }
template <typename T, typename... A> concept uninitialized_constructor = requires (T* p, const A&... a) { detail::uninitialized_construct(p, p, a...); };

// Members shared by vec and bounded_vec, which both store (sz, base).
// Self also provides owner(), the container pointer that debug iterators are checked against.
template <typename S> concept moved = !std::is_lvalue_reference_v<S>;
template <typename S> using pointer_t = std::conditional_t<std::is_const_v<S>, const typename S::value_type*, typename S::value_type*>;
struct contiguous_ops {
	template <typename Self> [[nodiscard]] bool empty(this const Self& self) { return !self.sz; }
	template <typename Self> [[nodiscard]] Self::size_type size(this const Self& self) { return self.sz; }

	template <typename Self> void assert_valid_deref(this const Self& self, const Self::value_type* p) {
		debug_assert(self.base <= p && p < self.base + self.sz, "iterator position {} out of range [0, {})", p - self.base, self.sz);
	}
	template <typename Self> void assert_valid_index(this const Self& self, Self::size_type i) {
		debug_assert(0 <= i && i < self.sz, "index {} out of range [0, {})", i, self.sz);
	}
	template <typename Self> void assert_valid_insert(this const Self& self, Self::size_type i) {
		debug_assert(0 <= i && i <= self.sz, "index {} out of range [0, {}]", i, self.sz);
	}
	template <typename Self> void assert_nonempty(this const Self& self) {
		debug_assert(self.sz > 0, "container is empty");
	}

	template <typename Self> pointer_t<Self> data(this Self& self) { return self.base; }
	template <typename Self> auto begin(this Self& self) {
		return std::conditional_t<std::is_const_v<Self>, typename Self::const_iterator, typename Self::iterator>(self.base, self.owner());
	}
	template <typename Self> auto end(this Self& self) {
		return std::conditional_t<std::is_const_v<Self>, typename Self::const_iterator, typename Self::iterator>(self.base + self.sz, self.owner());
	}
	template <typename Self> auto& operator [] (this Self& self, Self::size_type i) { self.assert_valid_index(i); return self.data()[i]; }
	template <typename Self> auto& front(this Self& self) { self.assert_nonempty(); return self.data()[0]; }
	template <typename Self> auto& back(this Self& self) { self.assert_nonempty(); return self.data()[self.sz-1]; }

	template <typename Self> bool operator == (this const Self& a, const Self& b) { return std::ranges::equal(a, b); }
	template <typename Self> auto operator <=> (this const Self& a, const Self& b) {
		return std::lexicographical_compare_three_way(a.begin(), a.end(), b.begin(), b.end());
	}

	// Conversions, transferring the allocation
	template <moved Self> auto into_bounded(this Self&& self) { return bounded_vec<typename Self::value_type>(self.release()); }
	template <moved Self> auto into_stack(this Self&& self) { return bounded_stack<typename Self::value_type>(self.release()); }
	template <moved Self> auto into_growable(this Self&& self) { return growable_vec<typename Self::value_type>(self.release()); }
};


struct insert_ops {
	template <typename Self>
	void push_back(this Self& self, const Self::value_type& v) requires std::copy_constructible<typename Self::value_type> {
		self.emplace_back(v);
	}
	template <typename Self>
	void push_back(this Self& self, Self::value_type&& v) {
		self.emplace_back(std::move(v));
	}

	// Add/remove to the middle
	template <typename Self, typename... Args>
	requires std::constructible_from<typename Self::value_type, Args&&...>
	Self::value_type& emplace(this Self& self, Self::size_type i, Args&&... args) {
		self.assert_valid_insert(i);
		self.emplace_back(std::forward<Args>(args)...);
		std::rotate(self.begin() + i, self.end() - 1, self.end());
		return self[i];
	}
	template <typename Self>
	void insert(this Self& self, Self::size_type i, const Self::value_type& v) requires std::copy_constructible<typename Self::value_type> {
		self.emplace(i, v);
	}
	template <typename Self>
	void insert(this Self& self, Self::size_type i, Self::value_type&& v) {
		self.emplace(i, std::move(v));
	}
	template <typename Self, typename... A> requires detail::uninitialized_constructor<typename Self::value_type, A...>
	void insert_n(this Self& self, Self::size_type i, Self::size_type n, const A&... a) {
		self.assert_valid_insert(i);
		auto [first, last] = self.grow_by(n, a...);
		std::rotate(self.begin() + i, first, last);
	}
	template <typename Self>
	void insert_n(this Self& self, Self::size_type i, Self::size_type n, const Self::value_type& v) requires std::copy_constructible<typename Self::value_type> {
		self.template insert_n<Self, typename Self::value_type>(i, n, v);
	}
	template <typename Self, std::ranges::input_range R>
	requires std::constructible_from<typename Self::value_type, std::ranges::range_reference_t<R>>
	void insert_range(this Self& self, Self::size_type i, R&& r) {
		self.assert_valid_insert(i);
		auto [first, last] = self.append_range(std::forward<R>(r));
		std::rotate(self.begin() + i, first, last);
	}

	template <typename Self>
	void erase(this Self& self, Self::size_type i) {
		self.assert_valid_index(i);
		std::shift_left(self.begin() + i, self.end(), 1);
		self.pop_back();
	}
	template <typename Self>
	void erase_n(this Self& self, Self::size_type i, Self::size_type n) {
		self.assert_valid_insert(i);
		self.assert_valid_insert(i + n);
		std::shift_left(self.begin() + i, self.end(), n);
		self.shrink_by(n);
	}

	template <typename Self, typename... Args>
	requires std::constructible_from<typename Self::value_type, Args&&...>
	Self::value_type& emplace_front(this Self& self, Args&&... args) {
		return self.emplace(0, std::forward<Args>(args)...);
	}
	template <typename Self> void push_front(this Self& self, const Self::value_type& v) requires std::copy_constructible<typename Self::value_type> { self.emplace_front(v); }
	template <typename Self> void push_front(this Self& self, Self::value_type&& v) { self.emplace_front(std::move(v)); }
	template <typename Self> void pop_front(this Self& self) { self.erase(0); }
};

} // namespace detail

template <typename T> struct vec : detail::contiguous_ops {
	static_assert(!std::is_const_v<T>);
	using value_type = T;
	using size_type = int;
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
	template <typename... A> requires detail::uninitialized_constructor<T, A...>
	explicit vec(size_type n, const A&... a) : vec(alloc_only{n}) { detail::uninitialized_construct(base, base + n, a...); }
	explicit vec(size_type n, const T& v) requires std::copy_constructible<T> : vec(alloc_only{n}) { detail::uninitialized_construct(base, base + n, v); }

	template <std::ranges::sized_range R> requires std::constructible_from<T, std::ranges::range_reference_t<R>>
	vec(std::from_range_t, R&& r) : vec(alloc_only{size_type(std::ranges::size(r))}) {
		std::ranges::uninitialized_copy(r, std::span(base, sz));
	}
	vec(std::initializer_list<T> list) requires std::copy_constructible<T> : vec(std::from_range, list) {}

	vec(const vec& o) requires std::copy_constructible<T> : vec(std::from_range, o) {}
	vec& operator=(const vec& o) requires std::copy_constructible<T> { *this = vec(o); return *this; }

	// Conversions, transferring the allocation
	explicit vec(detail::allocation<T> a) : sz(a.sz), base(a.p) {
		debug_assert(a.sz == a.cap, "adopting allocation with size {} not equal to capacity {}", a.sz, a.cap);
	}
	detail::allocation<T> release() {
		detail::allocation<T> a{base, sz, sz};
		base = nullptr;
		sz = 0;
		return a;
	}
	const vec* owner() const { return this; }
};

template <typename T> struct bounded_vec : detail::contiguous_ops, detail::insert_ops {
	static_assert(!std::is_const_v<T>);
	using value_type = T;
	using size_type = int;
	using iterator = detail::contiguous_iterator_impl<bounded_vec, false>;
	using const_iterator = detail::contiguous_iterator_impl<bounded_vec, true>;

	size_type sz = 0;
	size_type cap = 0;
	T* base = nullptr;

	bounded_vec() = default;

	friend void swap(bounded_vec& a, bounded_vec& b) noexcept { std::swap(a.sz, b.sz); std::swap(a.cap, b.cap); std::swap(a.base, b.base); }
	bounded_vec(bounded_vec&& o) noexcept : bounded_vec() { swap(*this, o); }
	bounded_vec& operator=(bounded_vec&& o) noexcept { swap(*this, o); return *this; }

	// Destruction
	~bounded_vec() { std::destroy_n(base, sz); detail::deallocate(base, cap); }

	// Construction
	explicit bounded_vec(with_capacity c) : cap(c.n), base(detail::allocate<T>(c.n)) {}
	template <typename... A> requires detail::uninitialized_constructor<T, A...>
	explicit bounded_vec(with_capacity c, size_type n, const A&... a) : bounded_vec(c) {
		debug_assert(n <= c.n, "size {} exceeds capacity {}", n, c.n);
		detail::uninitialized_construct(base, base + n, a...); sz = n;
	}
	template <typename... A> requires detail::uninitialized_constructor<T, A...>
	explicit bounded_vec(size_type n, const A&... a) : bounded_vec(with_capacity{n}, n, a...) {}
	explicit bounded_vec(with_capacity c, size_type n, const T& v) requires std::copy_constructible<T> : bounded_vec(c) { grow_by<T>(n, v); }
	explicit bounded_vec(size_type n, const T& v) requires std::copy_constructible<T> : bounded_vec(with_capacity{n}, n, v) {}
	template <std::ranges::input_range R>
	requires std::constructible_from<T, std::ranges::range_reference_t<R>>
	bounded_vec(with_capacity c, std::from_range_t, R&& r) : bounded_vec(c) {
		auto [in, out] = std::ranges::uninitialized_copy(r, std::ranges::subrange(base, base + cap));
		// Must have drained all of r
		debug_assert(in == std::ranges::end(r), "range size exceeds capacity {}", c.n);
		sz = size_type(out - base);
	}
	template <std::ranges::sized_range R>
	requires std::constructible_from<T, std::ranges::range_reference_t<R>>
	bounded_vec(std::from_range_t f, R&& r) : bounded_vec(with_capacity{size_type(std::ranges::size(r))}, f, r) {}
	bounded_vec(std::initializer_list<T> list) requires std::copy_constructible<T> : bounded_vec(std::from_range, list) {}
	bounded_vec(with_capacity c, std::initializer_list<T> list) requires std::copy_constructible<T> : bounded_vec(c, std::from_range, list) {}
	bounded_vec(const bounded_vec& o) requires std::copy_constructible<T> : bounded_vec(with_capacity{o.cap}, std::from_range, o) {}
	bounded_vec& operator=(const bounded_vec& o) requires std::copy_constructible<T> { *this = bounded_vec(o); return *this; }

	// Conversions, transferring the allocation
	explicit bounded_vec(detail::allocation<T> a) : sz(a.sz), cap(a.cap), base(a.p) {
		debug_assert(0 <= a.sz && a.sz <= a.cap, "adopting allocation with size {} out of range [0, {}]", a.sz, a.cap);
	}
	detail::allocation<T> release() {
		detail::allocation<T> a{base, sz, cap};
		base = nullptr;
		sz = cap = 0;
		return a;
	}
	const bounded_vec* owner() const { return this; }
	vec<T> into_vec() && {
		debug_assert(sz == cap, "converting container with size {} and capacity {} to vec", sz, cap);
		return vec<T>(release());
	}

	[[nodiscard]] size_type capacity() const { return cap; }

	// Resizing ops
	template <typename... Args> requires std::constructible_from<T, Args&&...>
	T& emplace_back(Args&&... args) {
		debug_assert(sz < cap, "pushing into full container (capacity {})", cap);
		return *std::construct_at(base + sz++, std::forward<Args>(args)...);
	}
	template <std::ranges::input_range R>
	requires std::constructible_from<T, std::ranges::range_reference_t<R>>
	std::ranges::subrange<iterator> append_range(R&& r) {
		auto old_end = end();
		// Use raw pointers since the memory is uninitialized
		auto [in, out] = std::ranges::uninitialized_copy(r, std::ranges::subrange(base + sz, base + cap));
		// Must have drained all of r
		debug_assert(in == std::ranges::end(r), "range size exceeds capacity {}", cap);
		sz = size_type(out - base);
		return {old_end, end()};
	}
	template <typename... A> requires detail::uninitialized_constructor<T, A...>
	std::ranges::subrange<iterator> grow_by(size_type n, const A&... a) {
		debug_assert(0 <= n, "growing container by negative {}", n);
		debug_assert(sz + n <= cap, "growing container with size {} by {} exceeds capacity {}", sz, n, cap);
		auto old_end = end();
		detail::uninitialized_construct(base + sz, base + sz + n, a...);
		sz += n;
		return {old_end, end()};
	}
	std::ranges::subrange<iterator> grow_by(size_type n, const T& v) requires std::copy_constructible<T> { return grow_by<T>(n, v); }
	template <typename... A> requires detail::uninitialized_constructor<T, A...>
	std::ranges::subrange<iterator> grow_to(size_type n, const A&... a) {
		debug_assert(sz <= n, "growing container with size {} to smaller target {}", sz, n);
		debug_assert(n <= cap, "growing container to {} exceeds capacity {}", n, cap);
		auto old_end = end();
		detail::uninitialized_construct(base + sz, base + n, a...);
		sz = n;
		return {old_end, end()};
	}
	std::ranges::subrange<iterator> grow_to(size_type n, const T& v) requires std::copy_constructible<T> { return grow_to<T>(n, v); }
	template <typename... A> requires detail::uninitialized_constructor<T, A...>
	std::ranges::subrange<iterator> grow_to(iterator it, const A&... a) {
		debug_assert(base + sz <= it.p, "growing container with size {} to smaller target {}", sz, it.p - base);
		debug_assert(it.p <= base + cap, "growing container to {} exceeds capacity {}", it.p - base, cap);
		auto old_end = end();
		detail::uninitialized_construct(base + sz, it.p, a...);
		sz = size_type(it.p - base);
		return {old_end, end()};
	}
	std::ranges::subrange<iterator> grow_to(iterator it, const T& v) requires std::copy_constructible<T> { return grow_to<T>(it, v); }

	// Remove from the back
	void pop_back() { assert_nonempty(); std::destroy_at(base + --sz); }
	void shrink_by(size_type n) {
		debug_assert(0 <= n && n <= sz, "shrinking container with size {} by {}", sz, n);
		std::destroy(base + sz - n, base + sz);
		sz -= n;
	}
	void shrink_to(size_type n) {
		debug_assert(0 <= n && n <= sz, "shrinking container with size {} to {}", sz, n);
		std::destroy(base + n, base + sz);
		sz = n;
	}
	void shrink_to(iterator it) {
		debug_assert(base <= it.p && it.p <= base + sz, "shrinking container with size {} to {}", sz, it.p - base);
		std::destroy(it.p, base + sz);
		sz = size_type(it.p - base);
	}

	void clear() { shrink_to(0); }
};

template <typename T> struct growable_vec : private bounded_vec<T> {
	using bounded = bounded_vec<T>;
	using typename bounded::value_type, typename bounded::size_type, typename bounded::iterator, typename bounded::const_iterator;
	// Inherit constructors
	using bounded::bounded_vec;
	using bounded::empty, bounded::size, bounded::capacity, bounded::data, bounded::begin, bounded::end, bounded::operator[], bounded::front, bounded::back;
	// Needed for insert_ops and other consumers
	using bounded::assert_valid_deref, bounded::assert_valid_index, bounded::assert_valid_insert, bounded::assert_nonempty;
	using bounded::pop_back, bounded::shrink_by, bounded::shrink_to, bounded::clear;
	// From insert_ops
	using bounded::push_back, bounded::emplace, bounded::insert, bounded::insert_n, bounded::insert_range, bounded::erase, bounded::erase_n, bounded::emplace_front, bounded::push_front, bounded::pop_front;

	friend void swap(growable_vec& a, growable_vec& b) noexcept { swap(static_cast<bounded&>(a), static_cast<bounded&>(b)); }

	// Conversions, transferring the allocation
	using bounded::release, bounded::sz, bounded::base, bounded::owner, bounded::operator==, bounded::operator<=>;
	vec<T> into_vec() && {
		shrink_to_fit();
		return std::move(static_cast<bounded&>(*this)).into_vec();
	}
	bounded_vec<T> into_bounded() && { return std::move(static_cast<bounded&>(*this)); }
	bounded_stack<T> into_stack() && { return std::move(static_cast<bounded&>(*this)).into_stack(); }

	// Update capacity when running growing ops
	void set_capacity(size_type n) {
		debug_assert(n >= size(), "reserving less capacity {} than size {}", n, size());
		bounded::operator=(bounded(with_capacity{n}, std::from_range, *this | std::views::as_rvalue));
	}
	void reserve_for(int n) {
		int c = capacity();
		if (n > c) {
			set_capacity(std::max(2 * c, n));
		}
	}
	void shrink_to_fit() { set_capacity(size()); }

	template <typename... Args> requires std::constructible_from<T, Args&&...>
	T& emplace_back(Args&&... args) {
		if (size() < capacity()) return bounded::emplace_back(std::forward<Args>(args)...);
		// args may refer to our own elements, so construct before reallocating
		T x(std::forward<Args>(args)...);
		reserve_for(size() + 1);
		return bounded::emplace_back(std::move(x));
	}
	template <std::ranges::input_range R>
	requires std::constructible_from<T, std::ranges::range_reference_t<R>>
	std::ranges::subrange<iterator> append_range(R&& r) {
		if constexpr (std::ranges::sized_range<R>) {
			reserve_for(size() + size_type(std::ranges::size(r)));
			return bounded::append_range(std::forward<R>(r));
		} else {
			size_type old_sz = size();
			for (auto&& x : r) emplace_back(std::forward<decltype(x)>(x));
			return {begin() + old_sz, end()};
		}
	}
	template <typename... A> requires detail::uninitialized_constructor<T, A...>
	std::ranges::subrange<iterator> grow_by(size_type n, const A&... a) {
		if (size() + n <= capacity()) return bounded::grow_by(n, a...);
		// a may refer to our own elements, so copy before reallocating
		std::tuple<A...> c(a...);
		reserve_for(size() + n);
		return std::apply([&](const A&... x) { return bounded::grow_by(n, x...); }, c);
	}
	std::ranges::subrange<iterator> grow_by(size_type n, const T& v) requires std::copy_constructible<T> { return grow_by<T>(n, v); }
	template <typename... A> requires detail::uninitialized_constructor<T, A...>
	std::ranges::subrange<iterator> grow_to(size_type n, const A&... a) {
		debug_assert(size() <= n, "growing container with size {} to smaller target {}", size(), n);
		return grow_by(n - size(), a...);
	}
	std::ranges::subrange<iterator> grow_to(size_type n, const T& v) requires std::copy_constructible<T> { return grow_to<T>(n, v); }
};

} // namespace wala
