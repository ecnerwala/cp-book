#pragma once

#include <cstdint>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <type_traits>
#include <array>

namespace dirichlet_series {

inline int inv(int v) {
	assert(v == 1);
	return 1;
}

inline int64_t inv(int64_t v) {
	assert(v == 1);
	return 1;
}

constexpr int64_t floor_sqrt(int64_t N) {
	assert(N >= 0);
	if (N == 0) return 0;
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
	int len = int(2 * rt + (rt * (rt+1) <= N));

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

	div_vector() = default;

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

template <div_vector_layout const& layout, typename T> class dirichlet_series_values;
template <div_vector_layout const& layout, typename T> class dirichlet_series_prefix;
template <div_vector_layout const& layout, typename T> class dirichlet_series_binary_indexed_tree;

template <div_vector_layout const& layout, typename T> class dirichlet_series_values : public div_vector<layout, T>, public vectorspace_mixin<layout, T, dirichlet_series_values<layout, T>> {
public:
	dirichlet_series_values() = default;

	template <typename F, std::enable_if_t<std::is_invocable_r_v<T, F, int64_t, int64_t>, bool> = true>
	dirichlet_series_values(F f) {
		for (int i = 1; i < layout.len; i++) {
			this->st[i] = f(layout.get_bucket_bound(i-1), layout.get_bucket_bound(i));
		}
	}

	template <typename U> explicit dirichlet_series_values(dirichlet_series_values<layout, U> const& o) {
		for (int i = 1; i < layout.len; i++) {
			this->st[i] = T(o.st[i]);
		}
	}

	explicit dirichlet_series_values(dirichlet_series_prefix<layout, T> && o) : div_vector<layout, T>(static_cast<div_vector<layout, T>&&>(std::move(o))) {
		for (int i = layout.len - 1; i > 1; i--) {
			this->st[i] -= this->st[i-1];
		}
	}
	explicit dirichlet_series_values(dirichlet_series_prefix<layout, T> const& o) {
		for (int i = layout.len - 1; i > 1; i--) {
			this->st[i] = o.st[i] - o.st[i-1];
		}
		this->st[1] = o.st[1];
	}
};

template <div_vector_layout const& layout, typename T> class dirichlet_series_prefix : public div_vector<layout, T>, public vectorspace_mixin<layout, T, dirichlet_series_prefix<layout, T>> {
public:
	dirichlet_series_prefix() = default;

	template <typename F, std::enable_if_t<std::is_invocable_r_v<T, F, int64_t>, bool> = true>
	dirichlet_series_prefix(F f) {
		for (int i = 1; i < layout.len; i++) {
			this->st[i] = f(layout.get_bucket_bound(i));
		}
	}

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

private:
	// This essentially runs *this += a * b, except it doesn't convolve any
	// terms involving 1*i and leaves those for the user-provided function f.
	// (f is called for each i in [2, layout.len-1].) This allows us to
	// easily implement multiplication or division or sqrt. (Note that a or b
	// are allowed to be equal to this.)
	template <typename F>
	void convolve_helper(dirichlet_series_prefix const& a, dirichlet_series_prefix const& b, F f) {
		// We roughly want to apply this[N/z] += a_val[x] * b_val[y] for all xyz <= N
		//
		// We'll split into the following cases (WLOG x <= y):
		// 0a. x = 1 or y = 1
		// 0b. x = y > 1
		// 1. x < y <= z <= N/x/y
		// 2. max(x, z) < y <= N/x/z

		T cur_sum = a.st[1] * b.st[1];
		for (int i = 2; i < layout.len; i++) {
			cur_sum += this->st[i];

			// Case 2: max(x, z) < y <= N/x/z
			// x^2 <= N / z
			// x <= N / z / z
			if (i >= layout.len - layout.rt) {
				int z = int(layout.len - i);
				assert(z <= layout.rt);
				int64_t rt_over_z = layout.rt/z;
				int64_t N_over_z = layout.N/z;
				int64_t x_max = N_over_z/(z+1);

				T tot_val = T();
				for (int64_t x = 2; x * (x+1) <= N_over_z && x <= x_max; x++) {
					// ylo = std::max(x, z)
					int ylo_idx = std::max(int(x), z);
					// yhi = N / x / z
					int yhi_idx = int(x <= rt_over_z ? layout.len - x * z : N_over_z / x);
					assert(ylo_idx < yhi_idx);

					T ax = a.st[x] - a.st[x-1];
					T bx = b.st[x] - b.st[x-1];
					T ay = a.st[yhi_idx] - a.st[ylo_idx];
					T by = b.st[yhi_idx] - b.st[ylo_idx];

					T v = ax * by + ay * bx;
					tot_val += v;
				}
				cur_sum += tot_val;
				if (i+1 < layout.len) {
					this->st[i+1] -= tot_val;
				}
			}

			this->st[i] = f(i, cur_sum);

			T ai = a.st[i] - a.st[i-1];
			T bi = b.st[i] - b.st[i-1];

			// Case 0a: x = 1
			cur_sum += ai * b.st[1] + a.st[1] * bi;

			if (i <= layout.rt) {
				// Case 1: x < y <= z <= N/x/y (y = i)
				// xy <= z <= N/y
				int64_t rt_over_i = layout.rt / i;
				int64_t N_over_i = layout.N / i;
				int x_max = int(std::min<int64_t>(N_over_i / i, i-1));
				T tot_sub = T();
				for (int x = 2; x <= x_max; x++) {
					T v;
					v = ai * (b.st[x] - b.st[x-1]) + (a.st[x] - a.st[x-1]) * bi;

					int zlo_idx = int(x <= rt_over_i ? x * i : layout.len - (N_over_i / x));
					this->st[zlo_idx] += v;
					tot_sub += v;
				}
				this->en[-(i-1)] -= tot_sub;

				// Case 0b: x = y > 1
				{
					int zlo_idx = int(i <= rt_over_i ? i * i : layout.len - (N_over_i / i));
					this->st[zlo_idx] += ai * bi;
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

	// This computes a pseudo-Euler transform of the sequence.
	//
	// Formally, given a Dirichlet series
	//   A = sum a_i / i^s,
	// we output the Dirichlet series corresponding to
	//   B = prod 1 / (1 - a_i i^{-s})
	//
	// Note: strictly speaking, the standard Euler transform over a generating function should be
	//   A = sum a_i / i^s -> B = prod 1 / (1 - i^{-s})^a_i
	// but our defintion is better suited for totally multiplicative functions,
	// and always works over general rings. Also, the two definitions match
	// when the a_i are always 0/1.
	//
	// This runs in $O(n^{2/3})$ time, but requires small inverses (up to 1/120).
	friend dirichlet_series_prefix euler_transform_fraction(dirichlet_series_prefix a_pref) {
		dirichlet_series_values<layout, T> a(std::move(a_pref));
		// assert(a.st[1] == 0);

		// Phase 0: stash away values up to the 6th root of N
		int x;
		for (x = 2; layout.rt / x / x / x > 0; x++) { }

		// Phase 1: adjust the values and insert the necessary extra powers
		std::array<T, 6> invs{T{}, T(1), inv(T(2)), inv(T(3)), inv(T(4)), inv(T(5))};
		for (int i = int(layout.rt); i >= x; i--) {
			T v = a.st[i];
			int e = 1;
			T pv = v;
			int64_t pi = i;
			while (pi <= layout.N/i) {
				e++;
				pi *= i;
				pv *= v;
				a.st[layout.get_value_bucket(pi)] += pv * invs[e];
			}
		}

		// Phase 2: now we take exp of the adjusted version
		// In particular, we take e^a = 1 + a + a^2 / 2 + a^3 / 6 + a^4 / 24 + a^5 / 120
		dirichlet_series_prefix v;
		for (int i = x; i < layout.len; i++) {
			v.st[i] = v.st[i-1] + a.st[i];
		}

		dirichlet_series_prefix r = v * v;
		for (int i = x; i < layout.len; i++) {
			r.st[i] = r.st[i] * invs[5] + v.st[i];
		}
		r *= v;
		for (int i = x; i < layout.len; i++) {
			r.st[i] = r.st[i] * invs[4] + v.st[i];
		}
		r *= v;
		for (int i = x; i < layout.len; i++) {
			r.st[i] = r.st[i] * invs[3] + v.st[i];
		}
		r *= v;
		for (int i = x; i < layout.len; i++) {
			r.st[i] = r.st[i] * invs[2] + v.st[i];
		}

		for (int i = 1; i < layout.len; i++) {
			r.st[i] += T(1);
		}

		// Phase 3: apply the extra below x
		for (x--; x >= 2; x--) {
			T ax = a.st[x];
			if (ax == 0) continue;
			for (int i = x; i < layout.len; i++) {
				r.st[i] += r.st[layout.get_value_bucket(layout.get_bucket_bound(i) / x)] * ax;
			}
		}

		return r;
	}

	// This computes the inverse of the pseudo-Euler transformation. See the
	// comment on euler_transform() for more details.
	friend dirichlet_series_prefix inverse_euler_transform_fraction(dirichlet_series_prefix a) {
		dirichlet_series_values<layout, T> r;

		// assert(a.st[1] == 1);

		// Phase 1: manually eliminate values up to the 6th root of a

		int x;
		for (x = 2; layout.rt / x / x / x > 0; x++) {
			T v = a.st[x] - T(1);
			if (v == 0) continue; // Small optimization, good for prime counting in particular
			r.st[x] = v;
			for (int i = layout.len - 1; i > x; i--) {
				a.st[i] -= a.st[layout.get_value_bucket(layout.get_bucket_bound(i) / x)] * v;
			}
			a.st[x] = T(1);
		}

		for (int i = 1; i < x; i++) {
			a.st[i] = T();
		}
		for (int i = x; i < layout.len; i++) {
			a.st[i] -= T(1);
		}

		std::array<T, 6> invs{T{}, T(1), inv(T(2)), inv(T(3)), inv(T(4)), inv(T(5))};

		// Phase 2: now we take log of the remaining thing, using just the first few terms.
		// In particular, we take log_a = a^5 / 5 - a^4 / 4 + a^3 / 3 - a^2 / 2 + a
		dirichlet_series_prefix log_a;
		for (int i = x; i < layout.len; i++) {
			log_a.st[i] = a.st[i] * invs[5];
		}
		log_a *= a;
		for (int i = x; i < layout.len; i++) {
			log_a.st[i] -= a.st[i] * invs[4];
		}
		log_a *= a;
		for (int i = x; i < layout.len; i++) {
			log_a.st[i] += a.st[i] * invs[3];
		}
		log_a *= a;
		for (int i = x; i < layout.len; i++) {
			log_a.st[i] -= a.st[i] * invs[2];
		}
		log_a *= a;
		for (int i = x; i < layout.len; i++) {
			log_a.st[i] += a.st[i] * invs[1];
		}

		// Phase 3: correct log_a; we need to get rid of the extra powers.
		for (int i = x; i < layout.len; i++) {
			r.st[i] = log_a.st[i] - log_a.st[i-1];
		}
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

		return dirichlet_series_prefix(std::move(r));
	}

	friend dirichlet_series_prefix euler_transform_binary_indexed_tree(dirichlet_series_prefix a_pref) {
		int x = 2;
		while (x <= layout.N / x / x) x++;

		dirichlet_series_prefix r_pref;
		for (int i = x; i < layout.len; i++) {
			r_pref.st[i] = a_pref.st[i] - a_pref.st[x-1];
		}

		for (int i = x; i <= layout.rt; i++) {
			T vi = a_pref.st[i] - a_pref.st[i-1];
			if (vi == 0) continue;

			int64_t N_over_i = layout.N / i;
			int64_t rt_over_i = layout.rt / i;
			int64_t max_z = N_over_i / i;
			for (int z = 1; z <= max_z; z++) {
				int jlo_idx = i-1;
				int jhi_idx = int(z <= rt_over_i ? layout.len - i * z : N_over_i / z);
				assert(jlo_idx < jhi_idx);
				T v = vi * (a_pref.st[jhi_idx] - a_pref.st[jlo_idx]);
				r_pref.st[layout.len-z] += v;
			}
		}

		dirichlet_series_binary_indexed_tree<layout, T> r(std::move(r_pref));
		r.increment_bucket_suffix(1, T(1));
		for (int i = x-1; i >= 2; i--) {
			T cur = a_pref.st[i] - a_pref.st[i-1];
			if (cur == 0) continue;
			r.sparse_mul_unlimited(i, cur);
		}
		return dirichlet_series_prefix(std::move(r));
	}

	friend dirichlet_series_prefix inverse_euler_transform_binary_indexed_tree(dirichlet_series_prefix a_pref) {
		// assert(a_pref.st[1] == 1);

		dirichlet_series_binary_indexed_tree<layout, T> a_bit(std::move(a_pref));
		dirichlet_series_values<layout, T> r;

		// First, use the BIT to clear up to N^1/3
		int x;
		for (x = 2; x <= layout.N / x / x; x++) {
			T cur = a_bit.get_bucket_prefix(x) - T(1);
			if (cur == 0) continue;
			r.st[x] = cur;
			a_bit.sparse_div_unlimited(x, cur);
		}
		a_pref = dirichlet_series_prefix<layout, T>(std::move(a_bit));

		// Now, a_pref contains terms of the form r[i] or r[i] * r[j], so let's
		// subtract out the semiprimes.
		// Note that N^1/3 < x <= i <= j, so N/i/j <= N^1/3 < x

		for (int i = x; i < layout.len; i++) {
			T vi = a_pref.st[i] - a_pref.st[i-1];
			r.st[i] = vi;
		}

		// We want i <= j <= N/i/z
		for (int i = x; i <= layout.rt; i++) {
			T vi = r.st[i];
			if (vi == 0) continue;

			int64_t N_over_i = layout.N / i;
			int64_t rt_over_i = layout.rt / i;
			int64_t max_z = N_over_i / i;
			for (int z = 1; z <= max_z; z++) {
				int jlo_idx = i-1;
				int jhi_idx = int(z <= rt_over_i ? layout.len - i * z : N_over_i / z);
				assert(jlo_idx < jhi_idx);
				T v = vi * (a_pref.st[jhi_idx] - a_pref.st[jlo_idx]);
				r.en[-z] -= v;
				if (z > 0) r.en[-(z-1)] += v;
			}
		}

		return dirichlet_series_prefix(std::move(r));
	}
};

// TODO: This will be useful for sparse convolution, which is nice for e.g. exp/log/prime counting
template <div_vector_layout const& layout, typename T> class dirichlet_series_binary_indexed_tree : public div_vector<layout, T>, public vectorspace_mixin<layout, T, dirichlet_series_binary_indexed_tree<layout, T>> {
public:
	dirichlet_series_binary_indexed_tree() = default;

	template <typename U> explicit dirichlet_series_binary_indexed_tree(dirichlet_series_binary_indexed_tree<layout, U> const& o) {
		for (int i = 1; i < layout.len; i++) {
			this->st[i] = T(o.st[i]);
		}
	}

	explicit dirichlet_series_binary_indexed_tree(dirichlet_series_prefix<layout, T> && o) : div_vector<layout, T>(static_cast<div_vector<layout, T>&&>(std::move(o))) {
		for (int i = layout.len - 1; i >= 1; i--) {
			this->st[i] -= this->st[i & (i-1)];
		}
	}
	explicit dirichlet_series_binary_indexed_tree(dirichlet_series_prefix<layout, T> const& o) {
		for (int i = layout.len - 1; i >= 1; i--) {
			this->st[i] = o.st[i] - o.st[i & (i-1)];
		}
	}
	explicit operator dirichlet_series_prefix<layout, T> () && {
		dirichlet_series_prefix<layout, T> r;
		swap(static_cast<div_vector<layout, T>&>(r), static_cast<div_vector<layout, T>&>(*this));
		for (int i = 1; i < layout.len; i++) {
			r.st[i] += r.st[i & (i-1)];
		}
		return r;
	}
	explicit operator dirichlet_series_prefix<layout, T> () const& {
		dirichlet_series_prefix<layout, T> r;
		for (int i = 1; i < layout.len; i++) {
			r.st[i] = this->st[i] + r.st[i & (i-1)];
		}
		return r;
	}

	T get_bucket_prefix(int a) const {
		T r = T();
		for (; a > 0; a -= a & -a) {
			r += this->st[a];
		}
		return r;
	}
	T get_prefix(int64_t v) const {
		return get_bucket_prefix(layout.get_value_bucket(v));
	}
	void increment_bucket_suffix(int a, T d) {
		for (; a < layout.len; a += a & -a) {
			this->st[a] += d;
		}
	}
	void increment_suffix(int64_t v, T d) const {
		return increment_bucket_suffix(layout.get_value_bucket(v), d);
	}

	// These 4 functions facilitate some simple sparse convolution.
	// They each take O(sqrt(N/x) log(N)) time.

	// multiply by (1 + w x^s)
	void sparse_mul_at_most_one(int64_t x, T w) {
		assert(x > 1);
		int64_t j = 1;
		T cur = get_bucket_prefix(int(j <= layout.rt / x ? layout.len - j * x : layout.N / x / j));
		for (; (j+1) <= layout.N / x / (j+1); j++) {
			T nxt = get_bucket_prefix(int(j + 1 <= layout.rt / x ? layout.len - (j+1) * x : layout.N / x / (j+1)));
			if (cur != nxt) {
				increment_bucket_suffix(int(layout.len - j), w * (cur - nxt));
				cur = nxt;
			}
		}
		for (int64_t i = layout.N / x / j; i > 0; i--) {
			T nxt = get_bucket_prefix(int(i-1));
			if (cur != nxt) {
				increment_bucket_suffix(int(i <= layout.rt / x ? i * x : layout.len - layout.N / x / i), w * (cur - nxt));
				cur = nxt;
			}
		}
	}

	// multiply by 1/(1 - w x^s) = 1 + wx^s + w^2 (x^2)^s + ...
	void sparse_mul_unlimited(int64_t x, T w) {
		assert(x > 1);
		T prv = T();
		int64_t i;
		for (i = 1; i <= layout.N / x / i; i++) {
			T cur = get_bucket_prefix(int(i));
			if (cur != prv) {
				increment_bucket_suffix(int(i <= layout.rt / x ? i * x : layout.len - layout.N / x / i), w * (cur - prv));
				prv = cur;
			}
		}
		for (int64_t j = layout.N / x / i; j >= 1; j--) {
			T cur = get_bucket_prefix(int(j <= layout.rt / x ? layout.len - j * x : layout.N / x / j));
			if (cur != prv) {
				increment_bucket_suffix(int(layout.len - j), w * (cur - prv));
				prv = cur;
			}
		}
	}

	// divide by (1 + w x^s)
	void sparse_div_at_most_one(int64_t x, T w) {
		return sparse_mul_unlimited(x, -w);
	}

	// divide by 1/(1 - w x^s) = 1 + wx^s + w^2 (x^2)^s + ...
	void sparse_div_unlimited(int64_t x, T w) {
		return sparse_mul_at_most_one(x, -w);
	}
};

}
