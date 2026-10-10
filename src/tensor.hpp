#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <compare>
#include <concepts>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <initializer_list>
#include <iterator>
#include <memory>
#include <ranges>
#include <span>
#include <type_traits>
#include <utility>

namespace wala {

#ifdef _GLIBCXX_DEBUG
inline constexpr bool WALA_DEBUG = true;
#else
inline constexpr bool WALA_DEBUG = false;
#endif

namespace detail {
[[noreturn, gnu::cold, gnu::noinline]] inline void check_fail(const char* what, std::ptrdiff_t i, std::ptrdiff_t n) {
	std::fprintf(stderr, "wala: %s %td out of range [0, %td)\n", what, i, n);
	std::abort();
}
}

inline void debug_assert(bool b) {
	if (!b) [[unlikely]] {
		if constexpr (WALA_DEBUG) {
			std::abort();
		} else {
			__builtin_unreachable();
		}
	}
}

inline void check_range(const char* what, std::ptrdiff_t i, std::ptrdiff_t n) {
	if (!(0 <= i && i < n)) [[unlikely]] {
		if constexpr (WALA_DEBUG) {
			detail::check_fail(what, i, n);
		} else {
			__builtin_unreachable();
		}
	}
}

struct uninit_t {};
inline constexpr uninit_t uninit{};
struct with_capacity { int n; };

template <typename T> struct vec;
template <typename T> struct bounded_vec;
template <typename T> struct bounded_stack;

namespace detail {

struct empty {};

template <typename T> T* allocate(int n) {
	debug_assert(n >= 0);
	return n ? std::allocator<T>().allocate(std::size_t(n)) : nullptr;
}
template <typename T> void deallocate(T* p, int n) {
	if (p) std::allocator<T>().deallocate(p, std::size_t(n));
}

template <typename C, bool Const> struct contiguous_iterator_impl {
	using Self = contiguous_iterator_impl;

	using T = std::conditional_t<Const, const typename C::value_type, typename C::value_type>;
	using iterator_concept = std::contiguous_iterator_tag;
	using iterator_category = std::random_access_iterator_tag;
	using value_type = std::remove_const_t<T>;
	using difference_type = int;

	T* p = nullptr;
	[[no_unique_address]] std::conditional_t<WALA_DEBUG, const C*, empty> c{};

	contiguous_iterator_impl() = default;
	contiguous_iterator_impl(T* p_, [[maybe_unused]] const C* c_) : p(p_) { if constexpr (WALA_DEBUG) c = c_; }
	contiguous_iterator_impl(const contiguous_iterator_impl&) = default;
	contiguous_iterator_impl& operator=(const contiguous_iterator_impl&) = default;
	contiguous_iterator_impl(const contiguous_iterator_impl<C, false>& o) requires Const : p(o.p), c(o.c) {}

	void check_deref() const { if constexpr (WALA_DEBUG) c->check_iter(p); }
	static void check_same([[maybe_unused]] Self a, [[maybe_unused]] Self b) { if constexpr (WALA_DEBUG) debug_assert(a.c == b.c); }

	T& operator*() const { check_deref(); return *p; }
	T* operator->() const { check_deref(); return p; }
	T& operator[](int n) const { return *(*this + n); }
	Self& operator++() { ++p; return *this; }
	Self operator++(int) { Self o = *this; operator++(); return o; }
	Self& operator--() { --p; return *this; }
	Self operator--(int) { Self o = *this; operator--(); return o; }
	Self& operator+=(int n) { p += n; return *this; }
	Self& operator-=(int n) { p -= n; return *this; }
	friend Self operator+(Self it, int n) { return it += n; }
	friend Self operator+(int n, Self it) { return it += n; }
	friend Self operator-(Self it, int n) { return it -= n; }
	friend int operator-(Self a, Self b) { check_same(a, b); return int(a.p - b.p); }
	friend bool operator==(Self a, Self b) { check_same(a, b); return a.p == b.p; }
	friend std::strong_ordering operator<=>(Self a, Self b) { check_same(a, b); return a.p <=> b.p; }
};

} // namespace detail
} // namespace wala

// std::to_address(it) must not go through operator-> (which is deref-checked); end() is a valid argument.
template <typename C, bool Const> struct std::pointer_traits<wala::detail::contiguous_iterator_impl<C, Const>> {
	using pointer = wala::detail::contiguous_iterator_impl<C, Const>;
	using element_type = typename pointer::T;
	using difference_type = int;
	static element_type* to_address(pointer it) noexcept { return it.p; }
};

namespace wala {
namespace detail {

// Everything derivable from data() and size(); Self owns the storage.
template <typename Self, typename T> struct contiguous_container {
	static_assert(!std::is_const_v<T>);
	using value_type = T;
	using size_type = int;
	using difference_type = int;
	using iterator = contiguous_iterator_impl<Self, false>;
	using const_iterator = contiguous_iterator_impl<Self, true>;

	bool empty() const { return self().size() == 0; }
	iterator begin() { return {self().data(), &self()}; }
	iterator end() { return {self().data() + self().size(), &self()}; }
	const_iterator begin() const { return {self().data(), &self()}; }
	const_iterator end() const { return {self().data() + self().size(), &self()}; }
	T& operator[](int i) { check_index(i); return self().data()[i]; }
	const T& operator[](int i) const { check_index(i); return self().data()[i]; }
	T& front() { check_nonempty(); return self().data()[0]; }
	const T& front() const { check_nonempty(); return self().data()[0]; }
	T& back() { check_nonempty(); return self().data()[self().size() - 1]; }
	const T& back() const { check_nonempty(); return self().data()[self().size() - 1]; }

	friend bool operator==(const Self& a, const Self& b) { return std::ranges::equal(a, b); }

	void check_index(int i) const { check_range("index", i, self().size()); }
	void check_iter(const T* p) const { check_range("iterator", p - self().data(), self().size()); }
	void check_nonempty() const { debug_assert(!empty()); }

private:
	Self& self() { return static_cast<Self&>(*this); }
	const Self& self() const { return static_cast<const Self&>(*this); }
};

}

// Fixed-size vector: always full, no size-changing operations.
template <typename T> struct vec : detail::contiguous_container<vec<T>, T> {
	T* base = nullptr;
	int sz = 0;

	vec() = default;
	~vec() { std::destroy_n(base, sz); detail::deallocate(base, sz); }

	friend void swap(vec& a, vec& b) noexcept { std::swap(a.base, b.base); std::swap(a.sz, b.sz); }
	vec(vec&& o) noexcept : vec() { swap(*this, o); }
	vec& operator=(vec&& o) noexcept { swap(*this, o); return *this; }
	vec(const vec&) = delete;
	vec& operator=(const vec&) = delete;

	explicit vec(int n, uninit_t) requires std::is_trivially_default_constructible_v<T> : base(detail::allocate<T>(n)), sz(n) {}
	explicit vec(int n) : base(detail::allocate<T>(n)), sz(n) { std::uninitialized_value_construct_n(base, n); }
	explicit vec(int n, const T& v) : base(detail::allocate<T>(n)), sz(n) { std::uninitialized_fill_n(base, n, v); }
	vec(std::initializer_list<T> il) : vec(std::from_range, il) {}
	template <std::ranges::sized_range R>
		requires std::constructible_from<T, std::ranges::range_reference_t<R>>
	vec(std::from_range_t, R&& r) : base(detail::allocate<T>(int(std::ranges::size(r)))), sz(int(std::ranges::size(r))) {
		std::ranges::uninitialized_copy(r, std::span(base, sz));
	}

	[[nodiscard]] vec clone() const { return vec(std::from_range, *this); }
	[[nodiscard]] bounded_vec<T> into_bounded() &&;
	[[nodiscard]] bounded_stack<T> into_stack() &&;

	int size() const { return sz; }
	T* data() { return base; }
	const T* data() const { return base; }
};

// Fixed-capacity vector: grow_to/shrink_to within capacity, never reallocates.
template <typename T> struct bounded_vec : detail::contiguous_container<bounded_vec<T>, T> {
	T* base = nullptr;
	int sz = 0;
	int cap = 0;

	bounded_vec() = default;
	~bounded_vec() { std::destroy_n(base, sz); detail::deallocate(base, cap); }

	friend void swap(bounded_vec& a, bounded_vec& b) noexcept { std::swap(a.base, b.base); std::swap(a.sz, b.sz); std::swap(a.cap, b.cap); }
	bounded_vec(bounded_vec&& o) noexcept : bounded_vec() { swap(*this, o); }
	bounded_vec& operator=(bounded_vec&& o) noexcept { swap(*this, o); return *this; }
	bounded_vec(const bounded_vec&) = delete;
	bounded_vec& operator=(const bounded_vec&) = delete;

	explicit bounded_vec(with_capacity c) : base(detail::allocate<T>(c.n)), cap(c.n) {}
	explicit bounded_vec(int n, uninit_t, with_capacity c) requires std::is_trivially_default_constructible_v<T> : bounded_vec(c) { grow_to(n, uninit); }
	explicit bounded_vec(int n, const T& v, with_capacity c) : bounded_vec(c) { grow_to(n, v); }
	bounded_vec(std::initializer_list<T> il) : bounded_vec(std::from_range, il) {}
	bounded_vec(std::initializer_list<T> il, with_capacity c) : bounded_vec(std::from_range, il, c) {}
	template <std::ranges::input_range R>
		requires std::constructible_from<T, std::ranges::range_reference_t<R>>
	bounded_vec(std::from_range_t, R&& r, with_capacity c) : bounded_vec(c) {
		if constexpr (std::ranges::sized_range<R>) {
			int n = int(std::ranges::size(r));
			debug_assert(n <= cap);
			std::ranges::uninitialized_copy_n(std::ranges::begin(r), n, base, base + n);
			sz = n;
		} else {
			for (auto&& x : r) emplace_back(std::forward<decltype(x)>(x));
		}
	}
	template <std::ranges::sized_range R>
		requires std::constructible_from<T, std::ranges::range_reference_t<R>>
	bounded_vec(std::from_range_t, R&& r) : bounded_vec(std::from_range, std::forward<R>(r), with_capacity{int(std::ranges::size(r))}) {}

	[[nodiscard]] bounded_vec clone() const { return bounded_vec(std::from_range, *this, with_capacity{cap}); }
	[[nodiscard]] vec<T> into_vec() &&;
	[[nodiscard]] bounded_stack<T> into_stack() &&;

	int size() const { return sz; }
	int capacity() const { return cap; }
	bool full() const { return sz == cap; }
	T* data() { return base; }
	const T* data() const { return base; }

	template <typename... Args>
	T& emplace_back(Args&&... args) { debug_assert(sz < cap); return *std::construct_at(base + sz++, std::forward<Args>(args)...); }
	void push_back(const T& v) { emplace_back(v); }
	void push_back(T&& v) { emplace_back(std::move(v)); }
	void pop_back() { debug_assert(sz > 0); std::destroy_at(base + --sz); }
	void clear() { shrink_to(0); }
	void shrink_to(int n) { debug_assert(0 <= n && n <= sz); std::destroy(base + n, base + sz); sz = n; }
	void grow_to(int n, uninit_t) requires std::is_trivially_default_constructible_v<T> { debug_assert(sz <= n && n <= cap); sz = n; }
	void grow_to(int n) { debug_assert(sz <= n && n <= cap); std::uninitialized_value_construct(base + sz, base + n); sz = n; }
	void grow_to(int n, const T& v) { debug_assert(sz <= n && n <= cap); std::uninitialized_fill(base + sz, base + n, v); sz = n; }
	void clear_and_set(int n, const T& v) { clear(); grow_to(n, v); }
};

// bounded_vec with end() stored as a pointer; cheaper when only the top is touched.
template <typename T> struct bounded_stack : detail::contiguous_container<bounded_stack<T>, T> {
	using Base = detail::contiguous_container<bounded_stack<T>, T>;
	using typename Base::iterator;
	using typename Base::const_iterator;

	T* base = nullptr;
	T* top = nullptr;
	T* lim = nullptr;

	bounded_stack() = default;
	~bounded_stack() { std::destroy(base, top); detail::deallocate(base, capacity()); }

	friend void swap(bounded_stack& a, bounded_stack& b) noexcept { std::swap(a.base, b.base); std::swap(a.top, b.top); std::swap(a.lim, b.lim); }
	bounded_stack(bounded_stack&& o) noexcept : bounded_stack() { swap(*this, o); }
	bounded_stack& operator=(bounded_stack&& o) noexcept { swap(*this, o); return *this; }
	bounded_stack(const bounded_stack&) = delete;
	bounded_stack& operator=(const bounded_stack&) = delete;

	explicit bounded_stack(with_capacity c) : base(detail::allocate<T>(c.n)), top(base), lim(base + c.n) {}
	explicit bounded_stack(int n, uninit_t, with_capacity c) requires std::is_trivially_default_constructible_v<T> : bounded_stack(bounded_vec<T>(n, uninit, c).into_stack()) {}
	explicit bounded_stack(int n, const T& v, with_capacity c) : bounded_stack(bounded_vec<T>(n, v, c).into_stack()) {}
	bounded_stack(std::initializer_list<T> il) : bounded_stack(bounded_vec<T>(il).into_stack()) {}
	bounded_stack(std::initializer_list<T> il, with_capacity c) : bounded_stack(bounded_vec<T>(il, c).into_stack()) {}
	template <std::ranges::input_range R>
		requires std::constructible_from<T, std::ranges::range_reference_t<R>>
	bounded_stack(std::from_range_t, R&& r, with_capacity c) : bounded_stack(bounded_vec<T>(std::from_range, std::forward<R>(r), c).into_stack()) {}
	template <std::ranges::sized_range R>
		requires std::constructible_from<T, std::ranges::range_reference_t<R>>
	bounded_stack(std::from_range_t, R&& r) : bounded_stack(bounded_vec<T>(std::from_range, std::forward<R>(r)).into_stack()) {}

	[[nodiscard]] bounded_stack clone() const { return bounded_stack(std::from_range, *this, with_capacity{capacity()}); }
	[[nodiscard]] vec<T> into_vec() &&;
	[[nodiscard]] bounded_vec<T> into_bounded() &&;

	int size() const { return int(top - base); }
	int capacity() const { return int(lim - base); }
	bool empty() const { return top == base; }
	bool full() const { return top == lim; }
	T* data() { return base; }
	const T* data() const { return base; }
	iterator end() { return {top, this}; }
	const_iterator end() const { return {top, this}; }

	template <typename... Args>
	T& emplace_back(Args&&... args) { debug_assert(top < lim); return *std::construct_at(top++, std::forward<Args>(args)...); }
	void push_back(const T& v) { emplace_back(v); }
	void push_back(T&& v) { emplace_back(std::move(v)); }
	void pop_back() { debug_assert(top > base); std::destroy_at(--top); }
	void clear() { shrink_to(0); }
	void shrink_to(int n) { debug_assert(0 <= n && n <= size()); std::destroy(base + n, top); top = base + n; }
	void grow_to(int n, uninit_t) requires std::is_trivially_default_constructible_v<T> { debug_assert(size() <= n && n <= capacity()); top = base + n; }
	void grow_to(int n) { debug_assert(size() <= n && n <= capacity()); std::uninitialized_value_construct(top, base + n); top = base + n; }
	void grow_to(int n, const T& v) { debug_assert(size() <= n && n <= capacity()); std::uninitialized_fill(top, base + n, v); top = base + n; }
	void clear_and_set(int n, const T& v) { clear(); grow_to(n, v); }
};

// Conversions steal the allocation; the source is left empty.
template <typename T> bounded_vec<T> vec<T>::into_bounded() && {
	bounded_vec<T> r;
	r.base = std::exchange(base, nullptr);
	r.sz = r.cap = std::exchange(sz, 0);
	return r;
}
template <typename T> bounded_stack<T> vec<T>::into_stack() && {
	bounded_stack<T> r;
	r.base = std::exchange(base, nullptr);
	r.top = r.lim = r.base + std::exchange(sz, 0);
	return r;
}
template <typename T> vec<T> bounded_vec<T>::into_vec() && {
	debug_assert(full());
	vec<T> r;
	r.base = std::exchange(base, nullptr);
	r.sz = std::exchange(sz, 0);
	cap = 0;
	return r;
}
template <typename T> bounded_stack<T> bounded_vec<T>::into_stack() && {
	bounded_stack<T> r;
	r.base = std::exchange(base, nullptr);
	r.top = r.base + std::exchange(sz, 0);
	r.lim = r.base + std::exchange(cap, 0);
	return r;
}
template <typename T> vec<T> bounded_stack<T>::into_vec() && {
	debug_assert(full());
	vec<T> r;
	r.sz = size();
	r.base = std::exchange(base, nullptr);
	top = lim = nullptr;
	return r;
}
template <typename T> bounded_vec<T> bounded_stack<T>::into_bounded() && {
	bounded_vec<T> r;
	r.sz = size();
	r.cap = capacity();
	r.base = std::exchange(base, nullptr);
	top = lim = nullptr;
	return r;
}

template <typename T, int NDIMS> struct tensor_view {
	static_assert(NDIMS >= 0, "NDIMS must be nonnegative");

protected:
	std::array<int, NDIMS> shape;
	std::array<int, NDIMS> strides;
	T* data;

	tensor_view(std::array<int, NDIMS> shape_, std::array<int, NDIMS> strides_, T* data_) : shape(shape_), strides(strides_), data(data_) {}

public:
	tensor_view() : shape{0}, strides{0}, data(nullptr) {}

protected:
	int flatten_index(std::array<int, NDIMS> idx) const {
		int res = 0;
		for (int i = 0; i < NDIMS; i++) { res += idx[i] * strides[i]; }
		return res;
	}
	int flatten_index_checked(std::array<int, NDIMS> idx) const {
		int res = 0;
		for (int i = 0; i < NDIMS; i++) {
			assert(0 <= idx[i] && idx[i] < shape[i]);
			res += idx[i] * strides[i];
		}
		return res;
	}

public:
	T& operator[] (std::array<int, NDIMS> idx) const {
#ifdef _GLIBCXX_DEBUG
		return data[flatten_index_checked(idx)];
#else
		return data[flatten_index(idx)];
#endif
	}
	T& at(std::array<int, NDIMS> idx) const {
		return data[flatten_index_checked(idx)];
	}

	template <int D = NDIMS>
	typename std::enable_if<(0 < D), tensor_view<T, NDIMS-1>>::type operator[] (int idx) const {
		std::array<int, NDIMS-1> nshape; std::copy(shape.begin()+1, shape.end(), nshape.begin());
		std::array<int, NDIMS-1> nstrides; std::copy(strides.begin()+1, strides.end(), nstrides.begin());
		T* ndata = data + (strides[0] * idx);
		return tensor_view<T, NDIMS-1>(nshape, nstrides, ndata);
	}
	template <int D = NDIMS>
	typename std::enable_if<(0 < D), tensor_view<T, NDIMS-1>>::type at(int idx) const {
		assert(0 <= idx && idx < shape[0]);
		return operator[](idx);
	}

	template <int D = NDIMS>
	typename std::enable_if<(0 == D), T&>::type operator * () const {
		return *data;
	}

	template <typename U, int D> friend struct tensor_view;
	template <typename U, int D> friend struct tensor;
};

template <typename T, int NDIMS> struct tensor {
	static_assert(NDIMS >= 0, "NDIMS must be nonnegative");

protected:
	std::array<int, NDIMS> shape;
	std::array<int, NDIMS> strides;
	vec<T> data;

public:
	tensor() : shape{0}, strides{0}, data() {}

	explicit tensor(std::array<int, NDIMS> shape_, const T& t = T()) {
		shape = shape_;
		int len = 1;
		for (int i = NDIMS-1; i >= 0; i--) {
			strides[i] = len;
			len *= shape[i];
		}
		data = vec<T>(len, t);
	}

	tensor(const tensor& o) : shape(o.shape), strides(o.strides), data(o.data.clone()) {}
	tensor& operator=(tensor&& o) noexcept {
		using std::swap;
		swap(shape, o.shape);
		swap(strides, o.strides);
		swap(data, o.data);
		return *this;
	}
	tensor(tensor&& o) noexcept : tensor() {
		*this = std::move(o);
	}
	tensor& operator=(const tensor& o) {
		return *this = tensor(o);
	}

	using view_t = tensor_view<T, NDIMS>;
	view_t view() {
		return tensor_view<T, NDIMS>(shape, strides, data.data());
	}
	operator view_t() {
		return view();
	}

	using const_view_t = tensor_view<const T, NDIMS>;
	const_view_t view() const {
		return tensor_view<const T, NDIMS>(shape, strides, data.data());
	}
	operator const_view_t() const {
		return view();
	}

	T& operator[] (std::array<int, NDIMS> idx) { return view()[idx]; }
	T& at(std::array<int, NDIMS> idx) { return view().at(idx); }
	const T& operator[] (std::array<int, NDIMS> idx) const { return view()[idx]; }
	const T& at(std::array<int, NDIMS> idx) const { return view().at(idx); }

	template <int D = NDIMS>
	typename std::enable_if<(0 < D), tensor_view<T, NDIMS-1>>::type operator[] (int idx) {
		return view()[idx];
	}
	template <int D = NDIMS>
	typename std::enable_if<(0 < D), tensor_view<T, NDIMS-1>>::type at(int idx) {
		return view().at(idx);
	}

	template <int D = NDIMS>
	typename std::enable_if<(0 < D), tensor_view<const T, NDIMS-1>>::type operator[] (int idx) const {
		return view()[idx];
	}
	template <int D = NDIMS>
	typename std::enable_if<(0 < D), tensor_view<const T, NDIMS-1>>::type at(int idx) const {
		return view().at(idx);
	}

	template <int D = NDIMS>
	typename std::enable_if<(0 == D), T&>::type operator * () {
		return *view();
	}
	template <int D = NDIMS>
	typename std::enable_if<(0 == D), const T&>::type operator * () const {
		return *view();
	}
};

} // namespace wala
