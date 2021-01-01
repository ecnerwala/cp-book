#pragma once

#include <bits/stdc++.h>

namespace nt { namespace arith_sum {

using ll = long long;

constexpr ll floor_sqrt(ll n) {
	if (n == 0) return 0; assert(n > 0);
	return std::llround(std::trunc(std::sqrt(n)));
}

template <ll N_, typename T> struct div_vector {
	static constexpr ll N = N_;
	static constexpr ll rt = floor_sqrt(N);
	static_assert(N > 0, "N must be positive");
	static_assert(N / rt >= rt, "Sqrt is wrong!");
	static_assert(N / (rt+1) < N / rt, "Sqrt not increasing!");

	static constexpr ll sz = rt + N / (rt+1);
	using container = std::vector<T>;
	container contents;

	div_vector() : contents(sz, T()) { }
	template <typename F> div_vector(const F& func) : contents(sz) {
		for (auto it : *this) {
			it.second = func(it.first);
		}
	}
	template <typename F> void fill(const F& func) {
		for (auto it : *this) {
			it.second = func(it.first);
		}
	}

	size_t to_index(ll v) const {
		assert(0 < v);
		assert(v <= N);

		if (v <= rt) {
			return size_t(v-1);
		} else {
			// v = N // d implies vd <= N, and d might as well be as large as possible, so d = N / v is valid
			ll d = N/v;
			assert(v == N / d);

			assert(d <= N / (rt+1));
			return size_t(contents.size()-d);
		}
	}
	size_t to_bucket(ll v) const {
		assert(0 < v);
		assert(v <= N);

		if (v <= rt) {
			return size_t(v-1);
		} else {
			ll d = N/v;

			assert(d <= N / (rt+1));
			return size_t(contents.size()-d);
		}
	}

	T& operator[] (ll v) {
		return contents[to_index(v)];
	}
	const T& operator[] (ll v) const {
		return contents[to_index(v)];
	}
	T& bucket(ll v) {
		return contents[to_bucket(v)];
	}
	const T& bucket(ll v) const {
		return contents[to_bucket(v)];
	}

	struct iterator {
		div_vector &dv;
		typename container::iterator it;
		ll v;

		std::pair<const ll, T&> operator* () {
			assert(v);
			return {v, *it};
		}

		iterator& operator ++ () {
			assert(v);
			++it;
			if (v == dv.N) {
				v = 0;
			} else if (v < dv.rt) {
				v ++;
			} else {
				v = dv.N / (dv.N / (v + 1));
			}
			return *this;
		}

		iterator operator ++ (int) {
			const_iterator res = *this;
			++(*this);
			return res;
		}

		bool operator == (const iterator &o) const { return &dv == &o.dv && it == o.it; }
		bool operator != (const iterator &o) const { return &dv != &o.dv || it != o.it; }
	};
	struct const_iterator {
		const div_vector &dv;
		typename container::const_iterator it;
		ll v;

		std::pair<const ll, const T&> operator* () {
			assert(v);
			return {v, *it};
		}

		const_iterator& operator ++ () {
			assert(v);
			++it;
			if (v == dv.N) {
				v = 0;
			} else if (v < dv.rt) {
				v ++;
			} else {
				v = dv.N / (dv.N / (v + 1));
			}
			return *this;
		}

		const_iterator operator ++ (int) {
			const_iterator res = *this;
			++(*this);
			return res;
		}

		bool operator == (const const_iterator &o) const { return &dv == &o.dv && it == o.it; }
		bool operator != (const const_iterator &o) const { return &dv != &o.dv || it != o.it; }
	};

	iterator begin() {
		return {*this, contents.begin(), 1};
	}
	iterator end() {
		return {*this, contents.end(), 0};
	}
	const_iterator begin() const {
		return {*this, contents.begin(), 1};
	}
	const_iterator end() const {
		return {*this, contents.end(), 0};
	}
	const_iterator cbegin() const {
		return {*this, contents.begin(), 1};
	}
	const_iterator cend() const {
		return {*this, contents.end(), 0};
	}

private:
	// finalize should finalize r[i] based on a[i] and b[i], or likewise
	// For normal convolution, this should run
	// r[i] += (i > 1) ? a[i] * b[1] + b[i] * a[1] : a[i] * b[i];
	// Should return the new prefix convolution
	template <typename F>
	friend void convolve_helper(div_vector& r, const div_vector& a, const div_vector& b, F finalize) {
		ll hyb = std::max(1ll, std::llround(std::pow(N, 2./3)));
		hyb = N / (N / hyb); // make sure it's a perfect bucket
		assert(hyb >= 1);

		T pa = 0;
		T pb = 0;
		T pr = 0;
		for (auto it : r) {
			if (it.first <= hyb) {
				it.second += pr;

				it.second -= pa * b[1];
				it.second -= pb * a[1];

				pr = finalize(it);

				T va = a[it.first] - pa;
				T vb = b[it.first] - pb;
				pa = a[it.first];
				pb = b[it.first];

				for (ll i = 2; i < it.first && i <= hyb / it.first; i++) {
					r.bucket(it.first * i) += va * (b[i] - (i>1 ? b[i-1] : 0)) + vb * (a[i] - (i>1 ? a[i-1] : 0));
				}
				if (2 <= it.first && it.first <= hyb / it.first) {
					r.bucket(it.first * it.first) += va * vb;
				}

			} else {
				ll curRt = floor_sqrt(it.first);
				assert(curRt <= r.rt);
				pa = a[1];
				for (ll i = 2; i <= curRt; i++) {
					ll j = it.first / i;
					assert(i <= j);

					it.second -= pa * b[j];
					it.second += a[i] * b[j];
					pa = a[i];
				}
				for (ll j = it.first / (curRt+1); j > 1; j--) {
					ll i = it.first / j;
					assert(j < i);

					it.second -= pa * b[j];
					it.second += a[i] * b[j];
					pa = a[i];
				}
				it.second -= pa * b[1];

				finalize(it);
			}
		}
	}

public:
	friend div_vector convolve(const div_vector& a, const div_vector& b) {
		div_vector r;
		convolve_helper(
			r, a, b,
			[&](auto it) -> T {
				if (it.first == 1) {
					it.second += a[1] * b[1];
				} else {
					it.second += a[it.first] * b[1];
					it.second += b[it.first] * a[1];
				}
				return it.second;
			}
		);
		return r;
	}

	friend div_vector operator * (const div_vector& a, const div_vector& b) {
		return convolve(a, b);
	}

	friend div_vector divide(const div_vector& a, const div_vector& b) {
		div_vector r;
		convolve_helper(
			r, r, b,
			[&](auto it) -> T {
				if (it.first > 1) {
					it.second += b[it.first] * r[1];
				}

				T pref = it.second;
				T rem = a[it.first] - pref;
				it.second = rem / b[1];
				assert(a[it.first] == pref + it.second * b[1]);
				return a[it.first];
			}
		);
		return r;
	}

	friend div_vector operator / (const div_vector& a, const div_vector& b) {
		return divide(a, b);
	}

	friend div_vector sqrt(const div_vector& a) {
		div_vector r;
		convolve_helper(
			r, r, r,
			[&](auto it) -> T {
				T pref = it.second;
				T rem = a[it.first] - pref;
				if (it.first == 1) {
					assert(rem == 1);
					it.second = 1;
					assert(a[it.first] == pref + it.second * it.second);
				} else {
					it.second = rem / 2 / r[1];
					assert(a[it.first] == pref + it.second * r[1] + it.second * r[1]);
				}
				return a[it.first];
			}
		);
		return r;
	}

	static div_vector DELTA() {
		return div_vector([](ll) -> T { return T(1); });
	}
	static div_vector ONES() {
		return div_vector([](ll v) -> T { return v; });
	}
	static div_vector ID() {
		return div_vector([](ll v) -> T { return v % 2 == 0 ? (T(v/2) * T(v+1)) : (T(v) * T((v+1)/2)); });
	}
	static div_vector POLY(int k) {
		// we first compute partial differences
		std::vector<T> diff(k+3);
		diff[0] = 0;
		for (int i = 1; i <= k+2; i++) {
			diff[i] = 1;
			for (int e = 0; e < k; e++) diff[i] *= i;
		}
		for (int i = 1; i <= k+1; i++) {
			for (int j = k+2; j > i; j--) {
				diff[j] -= diff[j-1];
			}
		}

		assert(diff[k+2] == 0);
		return div_vector([&](ll v) -> T {
			T res = 0;
			std::vector<ll> factors;
			for (int i = 0; i <= k+1; i++) {
				T binom = 1;
				for (ll factor : factors) {
					binom *= T(factor);
				}
				res += diff[i] * binom;
				factors.push_back(v-i);
				ll den = i+1;
				for (ll& factor : factors) {
					ll g = std::gcd(factor, den);
					den /= g;
					factor /= g;
				}
				assert(den == 1);
			}
			return res;
		});
	}
	static div_vector SIGMA() {
		return ONES() * ONES();
	}
	static div_vector SIGMA1() {
		return ONES() * ID();
	}
	static div_vector SIGMA(ll k) {
		return ONES() * POLY(k);
	}
};

}} // namespace nt::arith_sum
