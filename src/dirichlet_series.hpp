#pragma once

#include <cstdint>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <type_traits>

namespace dirichlet_series {

inline int64_t inv(int64_t v) {
	assert(v == 1);
	return 1;
}

constexpr int64_t floor_sqrt(int64_t N) {
	int64_t a = N;
	while (true) {
		int64_t b = N/a;
		assert(a >= b);
		if (a-b <= 1) return b;
		a = (a+b+1)>>1;
	}
}

class div_vector_layout {
public:
	int64_t N;
	int64_t rt = floor_sqrt(N);
	int64_t rt2 = rt - (rt * (rt + 1) > N);
	int len = int(1 + rt + rt2);

	constexpr div_vector_layout(int64_t N_ = 1) : N(N_) {}

	constexpr int get_value_bucket(int64_t a) const {
		return a <= rt ? int(a) : len - int(N/a);
	}
	constexpr int64_t get_bucket_bound(int i) const {
		return i <= rt ? i : N/(len-i);
	}
};

template <const div_vector_layout& layout, typename T> class div_vector {
public:
	// Let's just make everything public, getters and setters are too much work
	T* st = new T[layout.len+1]{}; // Allocate one extra on each side
	T* en = st + layout.len;

	div_vector() {}

	template <typename F, std::enable_if_t<std::is_invocable_r_v<T, F, int64_t>, bool> = true>
	div_vector(F f) {
		for (int i = 1; i < layout.len; i++) {
			st[i] = f(layout.get_bucket_bound(i));
		}
	}

	/* Rule of 5 declarations */
	div_vector(div_vector const& o) {
		std::copy(o.st, o.en, st);
	}
	div_vector& operator = (div_vector const& o) {
		std::copy(o.st, o.en, st);
		return *this;
	}
	friend void swap(div_vector& a, div_vector& b) {
		std::swap(a.st, b.st);
		std::swap(a.en, b.en);
	}
	div_vector(div_vector && o) : st(nullptr), en(nullptr) {
		swap(*this, o);
	}
	div_vector& operator = (div_vector && o) {
		swap(*this, o);
		return *this;
	}
	~div_vector() { delete[] st; }

	T& operator [] (int64_t v) { return st[layout.get_value_bucket(v)]; }
	T& operator [] (int64_t v) const { return st[layout.get_value_bucket(v)]; }
};

template <div_vector_layout const& layout, typename T, typename Derived> class vectorspace_mixin {
private:
	Derived& underlying() {
		return static_cast<Derived&>(*this);
	}
	Derived const& underlying() const {
		return static_cast<Derived const&>(*this);
	}

public:
	friend Derived operator + (Derived&& a) {
		for (int64_t i = 1; i < layout.len; i++) {
			a.st[i] = +a.st[i];
		}
		return a;
	}
	friend Derived operator + (Derived const& a) { return +Derived(a); }
	friend Derived operator - (Derived && a) {
		for (int64_t i = 1; i < layout.len; i++) {
			a.st[i] = -a.st[i];
		}
		return a;
	}
	friend Derived operator - (Derived const& a) { return -Derived(a); }

	Derived& operator += (Derived const& o) {
		for (int64_t i = 1; i < layout.len; i++) {
			underlying().st[i] += o.st[i];
		}
		return underlying();
	}
	friend Derived operator + (Derived && a, Derived const& b) { return a += b; }
	friend Derived operator + (Derived const& a, Derived && b) {
		for (int64_t i = 1; i < layout.len; i++) {
			b.st[i] = a.st[i] + b.st[i];
		}
		return b;
	}
	friend Derived operator + (Derived && a, Derived && b) { return std::move(a) + b; }
	friend Derived operator + (Derived const& a, Derived const& b) { return Derived(a) + b; }

	template <typename F> Derived& operator += (F f) {
		for (int64_t i = 1; i < layout.len; i++) {
			underlying().st[i] += f(layout.get_bucket_bound(i));
		}
		return underlying();
	}

	Derived& operator -= (Derived const& o) {
		for (int64_t i = 1; i < layout.len; i++) {
			underlying().st[i] -= o.st[i];
		}
		return underlying();
	}
	friend Derived operator - (Derived && a, Derived const& b) { return a -= b; }
	friend Derived operator - (Derived const& a, Derived && b) {
		for (int64_t i = 1; i < layout.len; i++) {
			b.st[i] = a.st[i] - b.st[i];
		}
		return b;
	}
	friend Derived operator - (Derived && a, Derived && b) { return std::move(a) - b; }
	friend Derived operator - (Derived const& a, Derived const& b) { return Derived(a) - b; }

	template <typename F> Derived& operator -= (F f) {
		for (int64_t i = 1; i < layout.len; i++) {
			underlying().st[i] -= f(layout.get_bucket_bound(i));
		}
		return underlying();
	}

	Derived& operator *= (T const& t) {
		for (int64_t i = 1; i < layout.len; i++) {
			underlying().st[i] *= t;
		}
		return underlying();
	}
	friend Derived operator * (Derived && a, T const& t) { return a *= t; }
	friend Derived operator * (Derived const& a, T const& t) { return Derived(a) * t; }
	// Just in case, don't assume multiplication is commutative.
	friend Derived operator * (T const& t, Derived && a) {
		for (int64_t i = 1; i < layout.len; i++) {
			a.st[i] = t * a.st[i];
		}
		return a;
	}
	friend Derived operator * (T const& t, Derived const& a) { return t * Derived(a); }

	Derived& operator /= (T const& t) {
		for (int64_t i = 1; i < layout.len; i++) {
			underlying().st[i] /= t;
		}
		return underlying();
	}
	friend Derived operator / (Derived && a, T const& t) { return a /= t; }
	friend Derived operator / (Derived const& a, T const& t) { return Derived(a) / t; }
};

template <div_vector_layout const& layout, typename T> class dirichlet_series_values : public div_vector<layout, T>, public vectorspace_mixin<layout, T, dirichlet_series_values<layout, T>> {
public:
	using div_vector<layout, T>::div_vector;

	template <typename U> explicit dirichlet_series_values(dirichlet_series_values<layout, U> const& o) {
		for (int i = 1; i < layout.len; i++) {
			this->st[i] = T(o.st[i]);
		}
	}
};

template <div_vector_layout const& layout, typename T> class dirichlet_series_prefix : public div_vector<layout, T>, public vectorspace_mixin<layout, T, dirichlet_series_prefix<layout, T>> {
public:
	using div_vector<layout, T>::div_vector;

	template <typename U> explicit dirichlet_series_prefix(dirichlet_series_prefix<layout, U> const& o) {
		for (int i = 1; i < layout.len; i++) {
			this->st[i] = T(o.st[i]);
		}
	}

	explicit dirichlet_series_prefix(dirichlet_series_values<layout, T> && o) : div_vector<layout, T>(static_cast<div_vector<layout, T>&&>(std::move(o))) {
		for (int i = 2; i < layout.len; i++) {
			this->st[i] += this->st[i-1];
		}
	}
	explicit dirichlet_series_prefix(dirichlet_series_values<layout, T> const& o) {
		T pref = this->st[1] = o.st[1];
		for (int i = 2; i < layout.len; i++) {
			this->st[i] = (pref += o.st[i]);
		}
	}

	explicit operator dirichlet_series_values<layout, T> () && {
		dirichlet_series_values<layout, T> r(static_cast<div_vector<layout, T>&&>(std::move(*this)));
		for (int i = layout.len - 1; i > 1; i--) {
			r.st[i] -= r.st[i-1];
		}
		return r;
	}

	explicit operator dirichlet_series_values<layout, T> () const& {
		dirichlet_series_values<layout, T> r;
		for (int i = layout.len - 1; i > 1; i--) {
			r.st[i] = this->st[i] - this->st[i-1];
		}
		r.st[1] = this->st[1];
		return r;
	}

private:
	// This essentially runs *this += a * b, except it doesn't convolve any
	// terms involving 1*i and leaves those for the user-provided function f.
	// (f is called for each i in [2, layout.len-1].) This allows us to
	// easily implement multiplication or division or sqrt. (Note that a or b
	// are allowed to be equal to this.)
	template <typename F>
	void convolve_helper(dirichlet_series_prefix const& a, dirichlet_series_prefix const& b, F f) {
		T cur_sum = a.st[1] * b.st[1];
		for (int i = 2; i < layout.len; i++) {
			cur_sum += this->st[i];

			// Handle cases with N/(z+1) < x * y <= N/z with max(x,y) > z, where z = layout.len - i
			if (i > layout.rt) {
				int z = int(layout.len - i);
				int64_t lo = layout.N/(z+1), hi = layout.N/z;
				int64_t rt_over_z = layout.rt / z;
				int64_t rt_over_z1 = layout.rt / (z+1);
				int64_t x_max = layout.N / z / (z+1);

				for (int64_t x = 2; x * x <= hi && x <= x_max; x++) {
					int yhi_idx = (x <= rt_over_z) ? layout.len - int(x*z) : int(hi / x);

					bool is_small = x * x <= lo;
					int ylo_idx = is_small
						? ((x <= rt_over_z1) ? layout.len - int(x*(z+1)) : int(lo / x))
						: int(x - 1);

					T ax = a.st[x] - a.st[x-1];
					T bx = b.st[x] - b.st[x-1];

					cur_sum += ax * (b.st[yhi_idx] - b.st[ylo_idx]) + (a.st[yhi_idx] - a.st[ylo_idx]) * bx;
					if (!is_small) {
						cur_sum -= ax * bx;
					}
				}
			}

			this->st[i] = f(i, cur_sum);

			T ai = a.st[i] - a.st[i-1];
			T bi = b.st[i] - b.st[i-1];

			cur_sum += ai * b.st[1] + a.st[1] * bi;

			// Handle cases with x * i with 2 <= x <= i <= N/x/i
			if (i <= layout.rt) {
				int64_t rt_over_i = layout.rt / i;
				int64_t N_over_i = layout.N / i;
				int x_max = int(std::min<int64_t>(N_over_i / i, i));
				for (int x = 2; x <= x_max; x++) {
					T v;
					if (x == i) {
						v = ai * bi;
					} else {
						v = ai * (b.st[x] - b.st[x-1]) + (a.st[x] - a.st[x-1]) * bi;
					}

					if (x <= rt_over_i) {
						this->st[x*i] += v;
					} else {
						this->en[-(N_over_i/x)] += v;
					}
				}
			}
		}
	}

public:
	friend dirichlet_series_prefix operator * (dirichlet_series_prefix const& a, dirichlet_series_prefix const& b) {
		dirichlet_series_prefix r;
		r.st[1] = a.st[1] * b.st[1];
		r.convolve_helper(a, b, [&](int i, T cur_sum) -> T {
			return cur_sum + (a.st[i] - a.st[i-1]) * b.st[1] + a.st[1] * (b.st[i] - b.st[i-1]);
		});
		return r;
	}
	dirichlet_series_prefix& operator *= (const dirichlet_series_prefix& o) { return *this = *this * o; }

	friend dirichlet_series_prefix operator / (dirichlet_series_prefix const& a, dirichlet_series_prefix const& b) {
		dirichlet_series_prefix r;
		T inv_b1 = inv(b.st[1]);
		r.st[1] = a.st[1] * inv_b1;
		r.convolve_helper(r, b, [&](int i, T cur_sum) -> T {
			return (a.st[i] - (cur_sum + r.st[1] * (b.st[i] - b.st[i-1]))) * inv_b1 + r.st[i-1];
		});
		return r;
	}
	dirichlet_series_prefix& operator /= (const dirichlet_series_prefix& o) { return *this = *this / o; }

	friend dirichlet_series_prefix sqrt(const dirichlet_series_prefix& a) {
		dirichlet_series_prefix r;
		// assert(a.st[1] == 1);
		r.st[1] = 1;
		T inv_2 = inv(T(2));
		r.convolve_helper(r, r, [&](int i, T cur_sum) -> T {
			return (a.st[i] - cur_sum) * inv_2 + r.st[i-1];
		});
		return r;
	}

	friend dirichlet_series_values<layout, T> inverse_euler_transform(dirichlet_series_prefix a) {
		// assert(a.st[1] == 1);

		// Phase 1: manually eliminate values up to the 6th root of a
		dirichlet_series_values<layout, T> small_values;

		int64_t x;
		for (x = 2; layout.rt / x / x / x > 0; x++) {
			T v = a.st[x] - T(1);
			if (v == 0) continue; // Small optimization, good for prime counting in particular
			small_values.st[x] = v;
			for (int i = layout.len - 1; i >= x; i--) {
				a.st[i] -= a.st[layout.get_value_bucket(layout.get_bucket_bound(i) / x)] * v;
			}
		}

		for (int i = 1; i < layout.len; i++) {
			a.st[i] -= T(1);
		}

		// Phase 2: now we take log of the remaining thing, using just the first few terms.
		// In particular, we take log_a = a^5 / 5 - a^4 / 4 + a^3 / 3 - a^2 / 2 + a
		dirichlet_series_prefix log_a;
		std::array<T, 6> invs{T{}, T(1), inv(T(2)), inv(T(3)), inv(T(4)), inv(T(5))};
		for (int i = 1; i < layout.len; i++) {
			log_a.st[i] = (a.st[i] - T(1)) * invs[5] - invs[4];
		}
		log_a *= a;
		for (int i = 1; i < layout.len; i++) {
			log_a.st[i] += invs[3];
		}
		log_a *= a;
		for (int i = 1; i < layout.len; i++) {
			log_a.st[i] -= invs[2];
		}
		log_a *= a;
		for (int i = 1; i < layout.len; i++) {
			log_a.st[i] += invs[1];
		}
		log_a *= a;

		// Phase 3: correct log_a; we need to get rid of the extra powers.
		dirichlet_series_values<layout, T> r = dirichlet_series_values<layout, T>(log_a);
		for (; x <= layout.rt; x++) {
			T v = r.st[x];
			int e = 1;
			T pv = v;
			int64_t px = x;
			while (px <= layout.N/x) {
				e++;
				px *= x;
				pv *= v;
				r.st[layout.get_value_bucket(px)] -= pv * invs[e];
			}
		}
		return r += small_values;
	}
};

// TODO: This will be useful for sparse convolution, which is nice for e.g. exp/log/prime counting
template <div_vector_layout const& layout, typename T> class dirichlet_series_seg : public div_vector<layout, T>, public vectorspace_mixin<layout, T, dirichlet_series_seg<layout, T>> {
public:
	using div_vector<layout, T>::div_vector;
};

}
