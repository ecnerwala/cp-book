#pragma once
// WIP!!!!

#include <cassert>
#include <functional>
#include <utility>

template <typename T, int& S, int MAXS, int HIGHS = -1, int NEGS = 0> struct myarray {
	static_assert(NEGS >= 0, "We must have at least zero");
	static_assert(HIGHS >= -1, "We must go up to S-1");
	T v[MAXS];
	T& operator [] (int i) {
		assert(S + (NEGS + HIGHS + 1) <= MAXS);
		assert(i >= -NEGS);
		assert(i <= S + HIGHS);
		return v[i + NEGS];
	}
	const T& operator [] (int i) const {
		assert(S + (NEGS + HIGHS + 1) <= MAXS);
		assert(i >= -NEGS);
		assert(i <= S + HIGHS);
		return v[i + NEGS];
	}

	struct slice {
		myarray* array;

		T* begin() {
			return &(*array)[0];
		}
		T* end() {
			return &(*array)[S];
		}
		const T* cbegin() const {
			return &(*array)[0];
		}
		const T* cend() const {
			return &(*array)[S];
		}
		const T* begin() const { return cbegin(); }
		const T* end() const { return cend(); }
	};
	slice values() {
		return slice{this};
	}
	const value_set values() const {
		return value_set{this};
	}

};
