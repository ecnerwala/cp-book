#pragma once

#include <vector>
#include <cassert>

/** Binary-indexed tree
 *
 *  A binary indexed tree with N nodes of type T implicitly provides the
 *  following two functions for 0 <= i <= N:
 *
 *      prefix(int i) -> list<T&>
 *      suffix(int i) -> list<T&>
 *
 *  such that size(suffix(i) intersect prefix(j)) = (1 if i < j else 0).
 *  Furthermore, the resulting lists always have size at most log_2(N).
 *
 *  This can be used to implement either point-update/(prefix|suffix)-query or
 *  (prefix|suffix)-update/point-query over a virtual array of size N of a
 *  commutative monoid. This can be generalized to implement
 *  point-update/range-query or range-update/point-query over a virtual array
 *  of size N of a commutative group.
 *
 *  With 0-indexed data, prefixes are more natural:
 *   * For range update/query, use for_prefix for the ranges and for_suffix for the points.
 *   * For prefix update/query, no change.
 *   * For suffix update/query, use for_prefix(point + 1); 1-index the data.
 */
template <typename T> struct bit {
private:
	std::vector<T> dat;
public:
	bit() {}
	explicit bit(size_t N) : dat(N) {}
	bit(size_t N, const T& t) : dat(N, t) {}

	template <typename F> void for_prefix(int i, F f) const {
		assert(0 <= i && i <= int(dat.size()));
		for (int a = i; a; a &= a-1) {
			f(dat[a-1]);
		}
	}

	template <typename F> void for_prefix(int i, F f) {
		assert(0 <= i && i <= int(dat.size()));
		for (int a = i; a; a &= a-1) {
			f(dat[a-1]);
		}
	}

	template <typename F> void for_suffix(int i, F f) const {
		assert(0 <= i && i <= int(dat.size()));
		for (int a = i; a < int(dat.size()); a |= a+1) {
			f(dat[a]);
		}
	}

	template <typename F> void for_suffix(int i, F f) {
		assert(0 <= i && i <= int(dat.size()));
		for (int a = i; a < int(dat.size()); a |= a+1) {
			f(dat[a]);
		}
	}
};
