#pragma once

#include <functional>
#include <utility>

template <typename F> struct reverse_comparator : private F {
	template <typename Arg1, typename Arg2> inline bool operator() (Arg1&& arg1, Arg2&& arg2) {
		return F::operator()(std::forward<Arg2>(arg2), std::forward<Arg1>(arg1));
	}
};

template <typename T, typename C = int, typename Comp = std::less<T>> struct cnt_min {
	T v;
	C cnt;

	cnt_min() : v(), cnt(0) {}
	explicit cnt_min(T v_) : v(v_), cnt(1) {}
	cnt_min(T v_, C cnt_) : v(v_), cnt(cnt_) {}

	friend cnt_min operator + (const cnt_min& a, const cnt_min& b) {
		if (!b.cnt) return a;
		else if (!a.cnt) return b;
		else if (Comp().operator()(a.v, b.v)) return a;
		else if (Comp().operator()(b.v, a.v)) return b;
		else return cnt_min(a.v, a.cnt + b.cnt);
	}

	cnt_min& operator += (const cnt_min& o) {
		return *this = (*this + o);
	}
};

template <typename T, typename C = int, typename Comp = std::less<T>> using cnt_max = cnt_min<T, C, reverse_comparator<Comp>>;
