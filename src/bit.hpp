#pragma once

#include <vector>
#include <cassert>

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
