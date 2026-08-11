#pragma once

#include <cassert>
#include <compare>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iostream>
#include <utility>
#include <vector>

// int_newtype<Tag, T>: a strongly-typed integer wrapper, intended for keeping
// separate index spaces (or other integer-like quantities) as separate types.
//
//   * Distinct Tag => distinct, non-interconvertible type.
//   * Explicit casts between newtypes (and to/from the underlying type).
//   * Implicit conversion from raw integers: from T itself at runtime, and
//     from compile-time constants of any integral type (range-checked) via a
//     consteval constructor. Conversion back to T is explicit.
//   * Ring semantics within a type: +, -, *, /, %, bit ops, and shifts are all
//     closed over the newtype (with 0 and 1 as distinguished points, and raw
//     integers embedding implicitly as operands); operands from other newtype
//     spaces always require an explicit cast.
template <typename Tag, typename T = int> struct int_newtype {
	using underlying_type = T;

	T v;

	int_newtype() : v(0) {}
	template <std::integral U> consteval int_newtype(U v_) : v(T(v_)) {
		assert(std::in_range<T>(v_));
	}
	constexpr int_newtype(T v_) : v(v_) {}
	template <typename Tag2, typename T2> explicit constexpr int_newtype(int_newtype<Tag2, T2> o) : v(T(o.v)) {}

	explicit constexpr operator T() const { return v; }

	friend constexpr bool operator == (int_newtype a, int_newtype b) { return a.v == b.v; }
	friend constexpr std::strong_ordering operator <=> (int_newtype a, int_newtype b) { return a.v <=> b.v; }

	constexpr int_newtype& operator ++ () { ++v; return *this; }
	constexpr int_newtype operator ++ (int) { int_newtype r = *this; ++v; return r; }
	constexpr int_newtype& operator -- () { --v; return *this; }
	constexpr int_newtype operator -- (int) { int_newtype r = *this; --v; return r; }

	constexpr int_newtype operator + () const { return *this; }
	constexpr int_newtype operator - () const { return int_newtype(T(-v)); }
	constexpr int_newtype operator ~ () const { return int_newtype(T(~v)); }

	constexpr int_newtype& operator += (int_newtype o) { v = T(v + o.v); return *this; }
	constexpr int_newtype& operator -= (int_newtype o) { v = T(v - o.v); return *this; }
	constexpr int_newtype& operator *= (int_newtype o) { v = T(v * o.v); return *this; }
	constexpr int_newtype& operator /= (int_newtype o) { v = T(v / o.v); return *this; }
	constexpr int_newtype& operator %= (int_newtype o) { v = T(v % o.v); return *this; }
	constexpr int_newtype& operator &= (int_newtype o) { v = T(v & o.v); return *this; }
	constexpr int_newtype& operator |= (int_newtype o) { v = T(v | o.v); return *this; }
	constexpr int_newtype& operator ^= (int_newtype o) { v = T(v ^ o.v); return *this; }
	constexpr int_newtype& operator <<= (int s) { v = T(v << s); return *this; }
	constexpr int_newtype& operator >>= (int s) { v = T(v >> s); return *this; }

	friend constexpr int_newtype operator + (int_newtype a, int_newtype b) { return int_newtype(T(a.v + b.v)); }
	friend constexpr int_newtype operator - (int_newtype a, int_newtype b) { return int_newtype(T(a.v - b.v)); }
	friend constexpr int_newtype operator * (int_newtype a, int_newtype b) { return int_newtype(T(a.v * b.v)); }
	friend constexpr int_newtype operator / (int_newtype a, int_newtype b) { return int_newtype(T(a.v / b.v)); }
	friend constexpr int_newtype operator % (int_newtype a, int_newtype b) { return int_newtype(T(a.v % b.v)); }
	friend constexpr int_newtype operator & (int_newtype a, int_newtype b) { return int_newtype(T(a.v & b.v)); }
	friend constexpr int_newtype operator | (int_newtype a, int_newtype b) { return int_newtype(T(a.v | b.v)); }
	friend constexpr int_newtype operator ^ (int_newtype a, int_newtype b) { return int_newtype(T(a.v ^ b.v)); }
	friend constexpr int_newtype operator << (int_newtype a, int s) { return int_newtype(T(a.v << s)); }
	friend constexpr int_newtype operator >> (int_newtype a, int s) { return int_newtype(T(a.v >> s)); }

	friend std::ostream& operator << (std::ostream& out, int_newtype n) { return out << n.v; }
	friend std::istream& operator >> (std::istream& in, int_newtype& n) { return in >> n.v; }
};

template <typename Tag, typename T> struct std::hash<int_newtype<Tag, T>> {
	size_t operator () (int_newtype<Tag, T> n) const { return std::hash<T>{}(n.v); }
};

template <typename> inline constexpr bool is_int_newtype_v = false;
template <typename Tag, typename T> inline constexpr bool is_int_newtype_v<int_newtype<Tag, T>> = true;
template <typename K> concept IntNewtype = is_int_newtype_v<K>;

// Keys for typed containers: either a newtype or a plain integer, so that
// structures can be templated on a key type with `int` as the default.
template <typename K> concept IntKey = IntNewtype<K> || std::integral<K>;
template <typename K> struct newtype_underlying { using type = K; };
template <typename Tag, typename T> struct newtype_underlying<int_newtype<Tag, T>> { using type = T; };
template <typename K> using newtype_underlying_t = typename newtype_underlying<K>::type;
template <IntKey K> constexpr newtype_underlying_t<K> newtype_raw(K k) {
	if constexpr (is_int_newtype_v<K>) return k.v;
	else return k;
}

// Declares a newtype with a fresh tag, e.g.:
//   INT_NEWTYPE(node_t);            // int_newtype<..., int>
//   INT_NEWTYPE(edge_t, int64_t);   // int_newtype<..., int64_t>
#define INT_NEWTYPE(name, ...) struct name ## _newtype_tag__; using name = int_newtype<name ## _newtype_tag__ __VA_OPT__(,) __VA_ARGS__>

// vec<K, V>: contiguous storage (std::vector<V>) indexed only by keys of
// type K (a newtype or a plain integer). Sizes and iteration are expressed
// in key-space.
template <IntKey K, typename V> struct vec {
	using key_type = K;
	using value_type = V;

	std::vector<V> data;

	vec() {}
	explicit vec(K n) : data(size_t(newtype_raw(n))) {}
	vec(K n, const V& init) : data(size_t(newtype_raw(n)), init) {}
	explicit vec(std::vector<V> data_) : data(std::move(data_)) {}
	vec(std::initializer_list<V> init) : data(init) {}

	V& operator [] (K k) { return data[size_t(newtype_raw(k))]; }
	const V& operator [] (K k) const { return data[size_t(newtype_raw(k))]; }
	V& at(K k) { return data.at(size_t(newtype_raw(k))); }
	const V& at(K k) const { return data.at(size_t(newtype_raw(k))); }

	// Past-the-end key; loop with: for (K i{}; i < a.size(); ++i)
	K size() const { return K(newtype_underlying_t<K>(data.size())); }
	bool empty() const { return data.empty(); }
	void resize(K n) { data.resize(size_t(newtype_raw(n))); }
	void resize(K n, const V& init) { data.resize(size_t(newtype_raw(n)), init); }
	void assign(K n, const V& init) { data.assign(size_t(newtype_raw(n)), init); }
	void reserve(K n) { data.reserve(size_t(newtype_raw(n))); }
	void clear() { data.clear(); }

	// Returns the key of the inserted element.
	K push_back(const V& val) { data.push_back(val); return K(newtype_underlying_t<K>(data.size() - 1)); }
	K push_back(V&& val) { data.push_back(std::move(val)); return K(newtype_underlying_t<K>(data.size() - 1)); }
	template <typename... Args> K emplace_back(Args&&... args) {
		data.emplace_back(std::forward<Args>(args)...);
		return K(newtype_underlying_t<K>(data.size() - 1));
	}
	void pop_back() { data.pop_back(); }
	V& back() { return data.back(); }
	const V& back() const { return data.back(); }
	V& front() { return data.front(); }
	const V& front() const { return data.front(); }

	// Value iteration (range-for over values, contiguous memory)
	auto begin() { return data.begin(); }
	auto begin() const { return data.begin(); }
	auto end() { return data.end(); }
	auto end() const { return data.end(); }

	// Key iteration: for (K i : a.keys()) { ... }
	struct key_range {
		K lo, hi;
		struct iterator {
			K k;
			K operator * () const { return k; }
			iterator& operator ++ () { ++k; return *this; }
			bool operator == (const iterator&) const = default;
		};
		iterator begin() const { return {lo}; }
		iterator end() const { return {hi}; }
	};
	key_range keys() const { return {K{}, size()}; }

	friend bool operator == (const vec& a, const vec& b) { return a.data == b.data; }
};
