#include <cassert>

namespace seg_tree {

// Floor of log_2(a); index of highest 1-bit
int log_2(int a) {
	return a ? (8 * sizeof(a)) - __builtin_clz(a) : -1;
}

template <typename F> void for_parents_up(int a, F f) {
	a >>= 1;
	for (; a > 0; a >>= 1) {
		f(a);
	}
}

template <typename F> void for_parents_down(int a, F f) {
	a >>= 1;
	for (int L = log_2(a); L >= 0; L--) {
		f(a >> L);
	}
}

template <typename F> void for_range_parents_up(int a, int b, F f) {
	assert(a && b);
	a >>= __builtin_ctz(a);
	b >>= __builtin_ctz(b);

	// TODO: dedup
	for_parents_up(a, f);
	for_parents_up(b, f);
}

template <typename F> void for_range_parents_down(int a, int b, F f) {
	assert(a && b);
	a >>= __builtin_ctz(a);
	b >>= __builtin_ctz(b);

	// TODO: dedup
	for_parents_down(a, f);
	for_parents_down(b, f);
}

template <typename F> void for_point(int a, F f) {
	for (; a > 0; a >>= 1) {
		f(a);
	}
}

template <typename F> void for_range(int a, int b, F f) {
	assert(a <= b && b <= 2*a);
	for (; a < b; a >>= 1, b >>= 1) {
		if (a & 1) f(a++, false);
		if (b & 1) f(--b, true);
	}
}

};
