#pragma once

#include <array>
#include <ranges>
#include <memory>

namespace wala {

#ifdef _GLIBCXX_DEBUG
inline constexpr bool WALA_DEBUG = true;
#else
inline constexpr bool WALA_DEBUG = false;
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

struct uninit_t {} uninit;
struct with_capacity { int n; };

namespace detail {

template <typename C, bool Const> struct contiguous_iterator_impl {
	using Self = contiguous_iterator_impl;

	using T = std::conditional_t<Const, const typename C::value_type, typename C::value_type>;
	using iterator_category = std::contiguous_iterator_tag;
	using iterator_concept = std::contiguous_iterator_tag;
	using value_type = std::remove_const_t<T>;
	using difference_type = std::ptrdiff_t;

	T* p = nullptr;
	[[no_unique_address]] std::conditional_t<WALA_DEBUG, const C*, std::monostate> c{};
	contiguous_iterator_impl() = default;
	contiguous_iterator_impl(T* p_, [[maybe_unused]] const C* c_) : p(p_) { if constexpr (WALA_DEBUG) c = c_; }
	contiguous_iterator_impl(const contiguous_iterator_impl& o) = default;
	contiguous_iterator_impl(const contiguous_iterator_impl<C, false>& o) requires Const : p(o.p), c(o.c) {}
	T& operator*() const { if constexpr (WALA_DEBUG) { c->check_iter(p); } return *p; }
	T* operator->() const { if constexpr (WALA_DEBUG) { c->check_iter(p); } return p; }
	T& operator[](difference_type n) const { return *(*this + n); }
	Self& operator++() { ++p; return *this; }
	Self operator++(int) { Self o = *this; operator++(); return o; }
	Self& operator--() { --p; return *this; }
	Self operator--(int) { Self o = *this; operator--(); return o; }
	Self& operator+=(difference_type n) { p += n; return *this; }
	friend Self operator+(Self it, difference_type n) { return it += n; }
	friend Self operator+(difference_type n, Self it) { return it += n; }
	Self& operator-=(difference_type n) { p -= n; return *this; }
	friend Self operator-(Self it, difference_type n) { return it -= n; }
	friend difference_type operator-(Self a, Self b) { return difference_type(a.p - b.p); }
	friend auto operator<=>(Self, Self) = default; // TODO: Assert that c is equal?
};
}

template <typename T> struct vec {
	static_assert(!std::is_const_v<T>);
	using value_type = T;
	using size_type = int;
	using difference_type = std::ptrdiff_t;
	using iterator = detail::contiguous_iterator_impl<vec, false>;
	using const_iterator = detail::contiguous_iterator_impl<vec, true>;

	T* base = nullptr;
	int sz = 0;

	vec() = default;
	~vec() { std::destroy_n(base, sz); if (base) std::allocator<T>().deallocate(base, sz); }

	friend void swap(vec& a, vec& b) noexcept { std::swap(a.base, b.base); std::swap(a.sz, b.sz); }
	vec(vec&& o) noexcept : vec() { swap(*this, o); }
	vec& operator= (vec&& o) noexcept { swap(*this, o); return *this; }

	explicit vec(int n) : base(std::allocator<T>().allocate(n)), sz(n) { std::uninitialized_value_construct_n(base, n); }
	explicit vec(int n, const T& v) : base(std::allocator<T>().allocate(n)), sz(n) { std::uninitialized_fill_n(base, n, v); }
	explicit vec(int n, uninit_t) : base(std::allocator<T>().allocate(n)), sz(n) { static_assert(std::is_trivially_default_constructible_v<T>); }
	template <std::ranges::sized_range R>
	vec(std::from_range_t, R&& r) : base(std::allocator<T>().allocate(std::ranges::size(r))), sz(int(std::ranges::size(r))) {
		std::ranges::uninitialized_copy(r, std::span(base, sz));
	}
	[[nodiscard]] vec<T> clone() const { return vec(std::from_range, *this); }
	int size() const { return sz; }
	bool empty() const { return sz == 0; }
	T* data() { return base; }
	const T* data() const { return base; }
	iterator begin() { return {base, this}; }
	iterator end() { return {base + sz, this}; }
	const_iterator begin() const { return {base, this}; }
	const_iterator end() const { return {base + sz, this}; }
	void check_index(int i) const {
		debug_assert(0 <= i);
		debug_assert(i < sz);
	}
	void check_iter(const T* p) const {
		debug_assert(base <= p);
		debug_assert(p < base + sz);
	}
	void check_nonempty() const { debug_assert(!empty()); }
	T& operator[] (int i) { check_index(i); return base[i]; }
	const T& operator[] (int i) const { check_index(i); return base[i]; }
	T& front() { check_nonempty(); return base[0]; }
	const T& front() const { check_nonempty(); return base[0]; }
	T& back() { check_nonempty(); return base[sz-1]; }
	const T& back() const { check_nonempty(); return base[sz-1]; }

	friend bool operator == (const vec& a, const vec& b) { return std::ranges::equal(a, b); }
};

template <typename T> struct bounded_vec {
	T* base = nullptr;
	int sz = 0;
	int cap = 0;

	bounded_vec() = default;
	~bounded_vec() { std::destroy_n(base, sz); if (base) std::allocator<T>().deallocate(base, cap); }

	friend void swap(bounded_vec& a, bounded_vec& b) noexcept { std::swap(a.base, b.base); std::swap(a.sz, b.sz); std::swap(a.cap, b.cap); }
	bounded_vec(bounded_vec&& o) noexcept : bounded_vec() { swap(*this, o); }
	bounded_vec& operator= (bounded_vec&& o) noexcept { swap(*this, o); return *this; }

	bounded_vec(const bounded_vec&) = delete;
	bounded_vec& operator= (const bounded_vec&) = delete;

	explicit bounded_vec(with_capacity c) : base(std::allocator<T>().allocate(c.n)), sz(0), cap(c.n) {}

	[[nodiscard]] int size() const { return sz; }
	[[nodiscard]] bool empty() const { return !sz; }
	[[nodiscard]] bool full() const { return sz == cap; }
	T* data() { return base; }
	const T* data() const { return base; }
	T* begin() { return base; }
	T* end() { return base + sz; }
	const T* begin() const { return base; }
	const T* end() const { return base + sz; }
	T& operator[] (int i) {
		return base[i];
	}
	const T& operator[] (int i) const { return base[i]; }
	T& front() { return base[0]; }
	const T& front() const { return base[0]; }
	T& back() { return base[sz-1]; }
	const T& back() const { return base[sz-1]; }

	template <typename... Args>
	T& emplace_back(Args&&... args) { assert(sz < cap); return *std::construct_at(base + sz++, std::forward<Args>(args)...); }
	void push_back(const T& v) { emplace_back(v); }
	void push_back(T&& v) { emplace_back(std::move(v)); }
	void pop_back() { assert(sz > 0); std::destroy_at(base + --sz); }
	void clear() { std::destroy(base, base + sz); sz = 0; }
	void shrink_to(int n) { assert(0 <= n && n <= sz); std::destroy(base + n, base + sz); sz = n; }
	void grow_to(int n, const T& v) { assert(sz <= n && n <= cap); std::uninitialized_fill(base + sz, base + n, v); sz = n; }
	void grow_to(int n) { assert(sz <= n && n <= cap); std::uninitialized_value_construct(base + sz, base + n); sz = n; }
	void assign(int n, const T& v) { clear(); grow_to(n, v); }
};

template <typename T> struct bounded_stack {
	T* base = nullptr;
	T* top = nullptr;
	T* cap = nullptr;

	bounded_stack() = default;
	~bounded_stack() { std::destroy(base, top); if (base) std::allocator<T>().deallocate(base, cap - base); }

	friend void swap(bounded_stack& a, bounded_stack& b) noexcept { std::swap(a.base, b.base); std::swap(a.top, b.top); std::swap(a.cap, b.cap); }
	bounded_stack(bounded_stack&& o) noexcept : bounded_stack() { swap(*this, o); }
	bounded_stack& operator= (bounded_stack&& o) noexcept { swap(*this, o); return *this; }

	bounded_stack(const bounded_stack&) = delete;
	bounded_stack& operator= (const bounded_stack&) = delete;

	explicit bounded_stack(with_capacity c) : base(std::allocator<T>().allocate(c.n)), top(base), cap(base + c.n) {}

	int size() const { return int(top - base); }
	bool empty() const { return top == base; }
	T* data() { return base; }
	const T* data() const { return base; }
	T* begin() { return base; }
	T* end() { return top; }
	const T* begin() const { return base; }
	const T* end() const { return top; }
	T& operator[] (int i) { return base[i]; }
	const T& operator[] (int i) const { return base[i]; }
	T& front() { return base[0]; }
	const T& front() const { return base[0]; }
	T& back() { return top[-1]; }
	const T& back() const { return top[-1]; }

	template <typename... Args>
	T& emplace_back(Args&&... args) { assert(top < cap); return *std::construct_at(top++, std::forward<Args>(args)...); }
	void push_back(const T& v) { emplace_back(v); }
	void push_back(T&& v) { emplace_back(std::move(v)); }
	void pop_back() { assert(top > base); std::destroy_at(--top); }
	void pop_to(T* ntop) { assert(ntop <= top); std::destroy(ntop, top); top = ntop; }
	// Only support shrinking; bulk delete
	void shrink_to(int n) { pop_to(base + n); }
	void clear() { shrink_to(0); }
};

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
