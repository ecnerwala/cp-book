#include <cassert>
#include <array>

namespace seg_tree {

// Floor of log_2(a); index of highest 1-bit
int log_2(int a) {
	return a ? (8 * sizeof(a)) - 1 - __builtin_clz(a) : -1;
}

int next_pow_2(int a) {
	assert(a > 0);
	return 1 << log_2(2*a-1);
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

template <typename F> void for_range_parents_up(std::array<int, 2> r, F f) {
	auto [a, b] = r;
	assert(a && b);
	a >>= __builtin_ctz(a);
	b >>= __builtin_ctz(b);

	// TODO: dedup
	for_parents_up(a, f);
	for_parents_up(b, f);
}

template <typename F> void for_range_parents_down(std::array<int, 2> r, F f) {
	auto [a, b] = r;
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

template <typename F> void for_range(std::array<int, 2> r, F f) {
	auto [a, b] = r;
	assert(a <= b && b <= 2*a);
	for (; a < b; a >>= 1, b >>= 1) {
		if (a & 1) f(a++, false);
		if (b & 1) f(--b, true);
	}
}

namespace in_order {

int get_point(int N, int a) {
	assert(0 <= a && a < N);
	int S = next_pow_2(N);
	a += S;
	return a >= 2 * N ? a - N : a;
}

std::array<int, 2> get_range(int N, std::array<int, 2> p) {
	auto [a, b] = p;
	assert(0 <= a && a <= b && b <= N);
	int S = next_pow_2(N);
	a += S, b += S;
	return { (a >= 2 * N ? 2*(a-N) : a), (b >= 2 * N ? 2*(b-N) : b) };
}

int get_node_index(int N, int a) {
	assert(N <= a && a < 2 * N);
	int S = next_pow_2(N);
	return (a < S ? a + N : a) - S;
}

std::array<int, 2> get_node_bounds(int N, int a) {
	assert(1 <= a && a < 2 * N);
	int l = __builtin_clz(a) - __builtin_clz(2*N-1);
	int S = next_pow_2(N);
	int x = a << l, y = (a+1) << l;
	assert(S <= x && x < y && y <= 2*S);
	return {(x >= 2 * N ? (x>>1) + N : x) - S, (y >= 2 * N ? (y>>1) + N : y) - S};
}

int get_node_size(int N, int a) {
	assert(1 <= a && a < 2 * N);
	auto [x, y] = get_node_bounds(N, a);
	return y - x;
}

} // namespace in_order

namespace circular {

int get_point(int N, int a) {
	assert(0 <= a && a < N);
	return N + a;
}

std::array<int, 2> get_range(int N, std::array<int, 2> p) {
	auto [a, b] = p;
	assert(0 <= a && a <= b && b <= N);
	return { N + a, N + b };
}

int get_node_index(int N, int a) {
	assert(N <= a && a < 2 * N);
	return a - N;
}

std::array<int, 2> get_node_bounds(int N, int a) {
	assert(1 <= a && a < 2 * N);
	int l = __builtin_clz(a) - __builtin_clz(2*N-1);
	int S = next_pow_2(N);
	int x = a << l, y = (a+1) << l;
	assert(S <= x && x < y && y <= 2*S);
	return {(x >= 2 * N ? x >> 1 : x) - N, (y >= 2 * N ? y >> 1 : y) - N};
}

int get_node_size(int N, int a) {
	assert(1 <= a && a < 2 * N);
	return in_order::get_node_size(N, a);
}

} // namespace circular

} // namespace seg_tree
